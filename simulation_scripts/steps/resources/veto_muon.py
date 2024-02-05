from __future__ import division
import numpy as np

from I3Tray import I3Tray
from icecube import icetray, dataclasses
from icecube.phys_services import I3Calculator

from ic3_labels.labels.utils import muon as mu_utils
from ic3_labels.labels.utils import detector
from ic3_labels.labels.utils import geometry


class InjectSingleVetoMuon(icetray.I3ConditionalModule):

    """Class to inject an accompanying muon for a provided neutrino event

    This is intended to run after neutrino injection. It is expected that the
    first primary in the provided I3MCTree is the injected neutrino event.
    For this particle, an accompanying muon will be injected at the convex hull
    of the detector.

    Attributes
    ----------
    mctree_name : str
        The name of the I3MCTree key.
    """

    def __init__(self, context):
        """Class to inject accompanying muons

        Parameters
        ----------
        context : TYPE
            Description
        """
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddOutBox('OutBox')
        self.AddParameter(
            'mctree_name', 'The name of the I3MCTree to use.', 'I3MCTree')
        self.AddParameter(
            'n_frames_per_neutrino',
            'The number of coincidence events to inject per imported event.',
            1)
        self.AddParameter(
            'sampling_settings',
            'Settings specifying how the muon energy is sampled.',
            {'method': 'power_law', 'range': [10, 1e7], 'gamma': 2})
        self.AddParameter(
            'uncorrelated_muon_settings',
            'If provided, the injected muon will be an uncorrelated coincident'
            ' muon sampled via the provided settings. The provided keys must '
            'include: anchor_(x/y/z)_range, time_range, '
            'azimuth_range [deg], zenith_range [deg]',
            None)
        self.AddParameter(
            'random_service',
            'The random service or seed to use. If this is an '
            'integer, a numpy random state will be created with '
            'the seed set to `random_service`',
            42)
        self.AddParameter(
            'output_key',
            'The key to which the sampling information is written to.',
            'MCVetoMuonInjectionInfo')

    def Configure(self):
        """Configures MuonLossProfileFilter.
        """
        self.mctree_name = self.GetParameter('mctree_name')
        self.n_frames_per_neutrino = self.GetParameter('n_frames_per_neutrino')
        self.sampling_settings = self.GetParameter('sampling_settings')
        self.uncorrelated_muon_settings = self.GetParameter(
            'uncorrelated_muon_settings')
        self.random_service = self.GetParameter('random_service')
        self.sampling_method = self.sampling_settings['method']
        self.output_key = self.GetParameter('output_key')

        if isinstance(self.random_service, int):
            self.random_service = np.random.RandomState(self.random_service)

    def DAQ(self, frame):
        """Inject accompanying muons
        """

        for i in range(self.n_frames_per_neutrino):

            # get I3MCTree
            mc_tree = dataclasses.I3MCTree(frame[self.mctree_name])

            # primary particle
            primary = mc_tree.get_primaries()[0]

            # inject uncorrelated muon with arbitrary anchor point
            if self.uncorrelated_muon_settings is not None:
                primary = dataclasses.I3Particle(primary)
                x = self.random_service.uniform(
                    *self.uncorrelated_muon_settings['anchor_x_range'])
                y = self.random_service.uniform(
                    *self.uncorrelated_muon_settings['anchor_y_range'])
                z = self.random_service.uniform(
                    *self.uncorrelated_muon_settings['anchor_z_range'])
                t = self.random_service.uniform(
                    *self.uncorrelated_muon_settings['time_range'])
                azimuth = self.random_service.uniform(*np.deg2rad(
                    self.uncorrelated_muon_settings['azimuth_range']))
                zenith = self.random_service.uniform(*np.deg2rad(
                    self.uncorrelated_muon_settings['zenith_range']))

                primary.pos = dataclasses.I3Position(x, y, z)
                primary.dir = dataclasses.I3Direction(zenith, azimuth)
                primary.time = t

            # compute entry point
            intersection_ts = mu_utils.get_muon_convex_hull_intersections(
                primary, convex_hull=detector.icecube_hull)
            if len(intersection_ts) == 0:
                # particle did not hit convex hull, use closest approach
                closest_position = I3Calculator.closest_approach_position(
                    primary, dataclasses.I3Position(0., 0., 0.))
                distance = mu_utils.get_distance_along_track_to_point(
                    primary.pos, primary.dir, closest_position)
            else:
                distance = min(intersection_ts)

            # in order to not land exactly on the convex hull and potentially
            # cause issues for label generation, we will walk back a little
            # further for primary and muon injection
            distance -= 1
            inj_pos = primary.pos + distance * primary.dir
            inj_time = primary.time + distance / dataclasses.I3Constants.c
            inj_dir = dataclasses.I3Direction(primary.dir)

            # inject new muon
            mc_tree, injection_info = self._create_mc_tree(
                inj_pos=inj_pos,
                inj_time=inj_time,
                inj_dir=inj_dir,
            )

            # copy frame
            frame_copy = icetray.I3Frame(frame)

            # replace I3MCTree
            frame_copy[self.mctree_name + 'VetoMuon_preMuonProp'] = mc_tree

            # add info to frame
            injection_info['injection_counter'] = float(i)
            if self.uncorrelated_muon_settings is not None:
                injection_info['is_correlated'] = False
            else:
                injection_info['is_correlated'] = True
            frame_copy[self.output_key] = injection_info

            # push frame on to subsequent modules
            self.PushFrame(frame_copy)

    def _powerlaw_sampler(self, e_min, e_max, gamma):
        """Sample from Powerlaw Distribution

        Sample `num` events from a power law with index gamma between e_min
        and e_max by using the analytic inversion method.
        The power law pdf is given by
        .. math::
           \mathrm{pdf}(\gamma) = x^{-\gamma} / \mathrm{norm}
        where norm ensures an area under curve of one. Positive spectral index
        gamma means a falling spectrum.
        Note: When :math:`\gamma=1` the integral is
        .. math::
           \int 1/x \mathrm{d}x = ln(x) + c
        This case is also handled.

        Parameters
        ----------
        e_min : float
            The lower bound of the PDF, needed for proper normalization.
        e_max : float
            The upper bound of the PDF, needed for proper normalization.
        gamma : float, optional
            Power law index.
        num : int
            Number of random numbers to generate.

        Returns
        -------
        float
            The sampled energy from the specified powerlaw distribution.
        """
        u = self.random_service.uniform(0., 1.)

        if gamma == 1:
            return np.exp(u * np.log(e_max / e_min)) * e_min
        else:
            radicant = (u * (e_max**(1. - gamma) - e_min**(1. - gamma))
                        + e_min**(1. - gamma))
            return radicant**(1. / (1. - gamma))

    def _sample_energy(self):
        """Sample Energy of the injected Muon

        Returns
        -------
        double
            The sampled energy
        """
        if self.sampling_method == 'power_law':
            energy = self._powerlaw_sampler(
                e_min=self.sampling_settings['range'][0],
                e_max=self.sampling_settings['range'][1],
                gamma=self.sampling_settings['gamma'],
            )

        else:
            raise ValueError('Unknown sampling method: {}'.format(
                self.sampling_method))
        return energy

    def _create_mc_tree(self, inj_pos, inj_time, inj_dir):
        """Inject accompanying muon in provided I3MCTree

        Parameters
        ----------
        inj_pos : I3Position
            The position at which to inject the muon.
        inj_time : double
            The time at which to inject the muon.
        inj_dir : I3Direction
            The direction of the injected muon.

        Returns
        -------
        I3MCTree
            The modified I3MCTree with the injected Muon.
        I3MapStringDouble
            Information on the injected muon.
        """
        mc_tree = dataclasses.I3MCTree()
        muon_primary = dataclasses.I3Particle()
        muon_primary.shape = dataclasses.I3Particle.ParticleShape.Primary
        muon_primary.dir = dataclasses.I3Direction(inj_dir)
        muon_primary.pos = dataclasses.I3Position(inj_pos)
        muon_primary.time = inj_time

        muon = dataclasses.I3Particle()
        muon.dir = dataclasses.I3Direction(inj_dir)
        muon.pos = dataclasses.I3Position(inj_pos)
        muon.time = inj_time
        muon.location_type = dataclasses.I3Particle.LocationType.InIce

        # sample type: MuPlus or MuMinus
        u = self.random_service.uniform(0., 1.)
        if u > 0.5:
            pdg_encoding = 13
        else:
            pdg_encoding = -13
        muon.pdg_encoding = pdg_encoding

        # sample energy
        muon.energy = self._sample_energy()

        # add muon primary to I3MCTree
        mc_tree.add_primary(muon_primary)

        # add muon as new child
        mc_tree.append_child(muon_primary, muon)

        # add info
        injection_info = dataclasses.I3MapStringDouble({
            'muon_energy': muon.energy,
            'muon_pdg_encoding': muon.pdg_encoding,
        })

        return mc_tree, injection_info


class CombineMCTrees(icetray.I3ConditionalModule):

    def __init__(self, context):
        """Class to combine two I3MCTrees

        Parameters
        ----------
        context : TYPE
            Description
        """
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddOutBox('OutBox')
        self.AddParameter('tree1', 'The name of the first I3MCTree.')
        self.AddParameter('tree2', 'The name of the first I3MCTree.')
        self.AddParameter(
            'output_key', 'Key to which the combined I3MCTree is written to.')

    def Configure(self):
        """Configures MuonLossProfileFilter.
        """
        self.tree1 = self.GetParameter('tree1')
        self.tree2 = self.GetParameter('tree2')
        self.output_key = self.GetParameter('output_key')

    def DAQ(self, frame):
        """Merge trees
        """
        tree1 = dataclasses.I3MCTree(frame[self.tree1])
        tree2 = dataclasses.I3MCTree(frame[self.tree2])
        tree1.merge(tree2)
        frame[self.output_key] = tree1
        self.PushFrame(frame)
