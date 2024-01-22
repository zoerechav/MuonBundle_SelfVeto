from __future__ import division
import numpy as np
np.set_printoptions(formatter={'all':lambda x: str(x)})
import click
import os
from I3Tray import I3Tray
from icecube import icetray, dataclasses,MuonGun
from icecube.phys_services import I3Calculator

from ic3_labels.labels.utils import muon as mu_utils
from ic3_labels.labels.utils import detector
from ic3_labels.labels.utils import geometry
import matplotlib.path as mpltPath
low_angles= np.linspace(0,0.5,13)
high_angles=np.linspace(0.5,1,13)
# new_low_angles=np.delete(low_angles,0) 
angles_space=np.concatenate( (low_angles,high_angles) )

depth_space = np.linspace(1.5,2.5,26)
flav_dict={}# include antiparticles once you've done

flav_dict=    {12:'NuE',-12:'NuE',14:'NuMu',-14:'NuMu',16:'NuTau',-16:'NuTau'}
deep=2.06
def get_surface_det(gcdFile=None):
        gcdFile=gcdFile
        bound_2D=[]
        surface_det = MuonGun.ExtrudedPolygon.from_file(gcdFile, padding=0)##Build Polygon from I3Geometry
        
        x=[(surface_det.x[i],surface_det.y[i])for i in range(len(surface_det.x))]###getting only x and y
        bound_2D=mpltPath.Path(x)#Projection of detector on x,y plane
        return bound_2D, surface_det

def boundary_check(particle1,gcdFile=None):
        ####checks if particle is inside the detector###
       
        inlimit = False  
        if ((particle1.pos.z <=max(surface_det.z)) and (particle1.pos.z>=min(surface_det.z))):
            if bound_2D.contains_points([(particle1.pos.x, particle1.pos.y)]):
                inlimit=True

        return inlimit
def default_sampler( e_min=1e1, e_max=1e7, gamma=3,multiplicity=1):
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
        u = np.random.uniform(0., 1.,multiplicity)

        if gamma == 1:
            return np.exp(u * np.log(e_max / e_min)) * e_min
        else:
            radicant = (u * (e_max**(1. - gamma) - e_min**(1. - gamma))
                        + e_min**(1. - gamma))
            return radicant**(1. / (1. - gamma))
class InjectVetoMuons(icetray.I3ConditionalModule):

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
        print(self.n_frames_per_neutrino)
        for i in range(self.n_frames_per_neutrino):
            print()
            print(frame['OldRunNumber'].value,';',frame['OldEventID'].value)
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
            global flavour
            flavour=flav_dict[primary.pdg_encoding]
            print('measuring multiplicity')
            try:
                bundle_multiplicity,depth,binnedsampling=self._sample_multiplicity(primary)
                print('multiplicity:{}, depth:{}, binnedsample:{}'.format(bundle_multiplicity,depth,binnedsampling))
                click.echo('multiplicity measured, going to muon injection')
                            # inject new muon
                print('NeutrinoEnergy',primary.energy)
                
                mc_tree, injection_info,energysample = self._create_bundle_mc_tree(
                    inj_energy=primary.energy,
                    inj_pos=inj_pos,
                    inj_time=inj_time,
                    inj_dir=inj_dir,
                    depth=depth,
                    multiplicity=bundle_multiplicity,
                )
            except Exception as e:
                print(e)
                # energysample=False
                # # inject new muon
                # mc_tree, injection_info = self._create_mc_tree(
                #     inj_pos=inj_pos,
                #     inj_time=inj_time,
                #     inj_dir=inj_dir,
                # )

            # copy frame
            frame_copy = icetray.I3Frame(frame)

            # replace I3MCTree
            frame_copy[self.mctree_name + 'VetoMuon_preMuonProp'] = mc_tree

            # add info to frame
            try:
                injection_info['injection_counter'] = float(i)
                if binnedsampling:
                    injection_info['MultiplicitySample'] = True
                else:
                    injection_info['MultiplicitySample'] = False
                if energysample:
                    injection_info['EnergySample'] = True
                else:
                    injection_info['EnergySample'] = False
                frame_copy[self.output_key] = injection_info

            # push frame on to subsequent modules
                self.PushFrame(frame_copy)
            except Exception as e:
                print(e)
            

    
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

    def _sample_multiplicity(self,primary):
        """Sample Multiplicity of the Muon bundle for given neutrino

        Returns
        -------
        int
            The bundle multiplicity
        """
        from random import seed
        import random
        mult_bins=np.linspace(0,100,100+1)
        mult_bin_centers = np.sqrt(mult_bins[:-1] * mult_bins[1:])
        bound_2D,surface_det = get_surface_det(gcdFile='/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz')

        intersection=surface_det.intersection(primary.pos, primary.dir)#points of intersection
        z_inter=primary.pos.z-intersection.first*np.cos(primary.dir.zenith)
        depth=1948.07-z_inter
        
        depth=np.around(depth_space[np.digitize(depth/1000,depth_space)-1],decimals=2)
        coszen=np.around(angles_space[np.digitize(np.cos(primary.dir.zenith),angles_space)-1],decimals=2)
        E_nu_bins=np.logspace(1,7,10+1)
        binnedsample=False
        
        try:
            #sample multiplicity from pdf for specific flavour of neutrinos 
            if os.path.exists('/data/user/vbasu/SelfVetoArrays/'+flavour+'_Mult_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'):
                Hist2D=np.load('/data/user/vbasu/SelfVetoArrays/'+flavour+'_Mult_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy')
                print('MultFile:','/data/user/vbasu/SelfVetoArrays/'+flavour+'_Mult_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy')
            #sample multiplicity from pdf for all neutrinos 
            elif os.path.exists('/data/user/vbasu/SelfVetoArrays/Mult_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'):
                Hist2D=np.load('/data/user/vbasu/SelfVetoArrays/Mult_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy')
                print('MultFile:','/data/user/vbasu/SelfVetoArrays/Mult_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy')
            energy_index=np.digitize(primary.energy,E_nu_bins)
            Hist1D=Hist2D[:,energy_index]
            Hist1D=Hist1D/np.sum(Hist1D)
            
            multiplicity= np.random.choice(np.arange(len(mult_bin_centers)), 1, p=Hist1D)[0]
            if multiplicity<1:
                multiplicity=1
            binnedsample=True
        except Exception as e:
            print(e)
            multiplicity=1
        print('multiplicity:{}, depth:{}, binnedsample:{}'.format(multiplicity,depth,binnedsample))
        return multiplicity,depth,binnedsample
        
    def _sample_bundle_energies(self,multiplicity,neutrino_energy,coszen,depth):
        """Sample Energy of the injected Muon

        Returns
        -------
        double
            The sampled energy
        """
        E_nu_bins=np.logspace(1,7,10+1)
        gridpts=np.logspace(-1,7,20)
        index=(np.linspace(0,19,19+1,dtype=int))
        energy_index=np.digitize(neutrino_energy,E_nu_bins)-1
        inj_muon_energies=[]
        if coszen<=0:
            print('CosZen <0')
            return inj_muon_energies,False
        coszen=angles_space[np.digitize(coszen,angles_space)]
        depth=depth_space[np.digitize(depth,depth_space)-1]
        coszen=np.around(coszen,decimals=2)
        depth=np.around(depth,decimals=2)
        energysample=True
        print('Bundle Energy Sampler')
        try:
            if multiplicity==1:
                filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_SingleSpline_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'
                print(filename)
                Hist2D=np.load('/data/user/vbasu/SelfVetoArrays/'+flavour+'_SingleSpline_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy',allow_pickle=True)
                weight_array=Hist2D[:,energy_index].T/np.sum(Hist2D[:,energy_index])
                gridpts=np.logspace(-1,7,len(Hist2D))
                index=(np.linspace(0,len(Hist2D)-1,len(Hist2D),dtype=int))
                indexarray=np.arange(len(index)**2)
                probarray=(weight_array).ravel()
                sample_indices = np.random.choice(len(probarray), 1, p=probarray)[0]
                weight_index=np.unravel_index(sample_indices, (len(index)))

                inj_muon_energies=[gridpts[weight_index]]
            elif multiplicity==2:
                filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_DoubleSpline_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'
                if os.path.exists(filename):#Spline exists
                    Hist3D=np.load('/data/user/vbasu/SelfVetoArrays/'+flavour+'_DoubleSpline_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy',allow_pickle=True)
                else:#Spline does not exist, using raw histogram instead
                    Hist3D=np.load('/data/user/vbasu/SelfVetoArrays/'+flavour+'_Double_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy',allow_pickle=True)
                weight_array=Hist3D[:,:,energy_index].T/np.sum(Hist3D[:,:,energy_index])
                gridpts=np.logspace(-1,7,len(Hist3D))
                index=(np.linspace(0,len(Hist3D)-1,len(Hist3D),dtype=int))
                indexarray=np.arange(len(index)**2)
                probarray=(weight_array).ravel()
                sample_indices = np.random.choice(len(probarray), 1, p=probarray)[0]
                weight_index=np.unravel_index(sample_indices, (len(index), len(index)))
                
                inj_muon_energies=[gridpts[weight_index[0]],gridpts[weight_index[1]]]
            elif multiplicity==3:
                filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_TripleSpline_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'
                if os.path.exists(filename):#Spline exists
                    Hist4D=np.load('/data/user/vbasu/SelfVetoArrays/'+flavour+'_TripleSpline_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy',allow_pickle=True)
                else:#Spline does not exist, using raw histogram instead
                    Hist4D=np.load('/data/user/vbasu/SelfVetoArrays/'+flavour+'_Triple_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy',allow_pickle=True)
                weight_array=Hist4D[:,:,:,energy_index].T/np.sum(Hist4D[:,:,:,energy_index])
                gridpts=np.logspace(-1,7,len(Hist4D))
                index=(np.linspace(0,len(Hist4D)-1,len(Hist4D),dtype=int))
                indexarray=np.arange(len(index)**3)
                probarray=(weight_array).ravel()
                sample_indices = np.random.choice(len(probarray), 1, p=probarray)[0]
                weight_index=np.unravel_index(sample_indices, (len(index), len(index),len(index)))
                
                inj_muon_energies=[gridpts[weight_index[0]],gridpts[weight_index[1]],gridpts[weight_index[2]]]
            elif multiplicity>=4:
                filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_QuadrupleSpline_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'
                if os.path.exists(filename):#Spline exists
                    Hist5D=np.load('/data/user/vbasu/SelfVetoArrays/'+flavour+'_QuadrupleSpline_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy',allow_pickle=True)
                else:#Spline does not exist, using raw histogram instead
                    Hist5D=np.load('/data/user/vbasu/SelfVetoArrays/Quadruple_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy',allow_pickle=True)
                weight_array=Hist5D[:,:,:,:,energy_index].T/np.sum(Hist5D[:,:,:,:,energy_index])
                gridpts=np.logspace(-1,7,len(Hist5D))
                index=(np.linspace(0,len(Hist5D)-1,len(Hist5D),dtype=int))
                indexarray=np.arange(len(index)**5)
                probarray=(weight_array).ravel()
                sample_indices = np.random.choice(len(probarray), 1, p=probarray)[0]
                weight_index=np.unravel_index(sample_indices, (len(index), len(index),len(index), len(index)))

                inj_muon_energies=[gridpts[weight_index[0]],gridpts[weight_index[1]],gridpts[weight_index[2]],gridpts[weight_index[3]]]
        
        except Exception as e:
            print(e)
            filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_EnergyTotalSplined_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'
            if os.path.exists(filename):#Use splined histogram for total energy
                print('EnergySum used')
                Hist2D=np.load('/data/user/vbasu/SelfVetoArrays/'+flavour+'_EnergyTotalSplined_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy',allow_pickle=True)
                weight_array=Hist2D[:,energy_index].T/np.sum(Hist2D[:,energy_index])
                gridpts=np.logspace(-1,7,len(Hist2D))
                index=(np.linspace(0,len(Hist2D)-1,len(Hist2D),dtype=int))
                indexarray=np.arange(len(index)**2)
                probarray=np.nan_to_num((weight_array).ravel())/np.sum(np.nan_to_num((weight_array).ravel()))
                if not np.isnan(np.sum(probarray)):
                    sample_indices = np.random.choice(len(probarray), 1, p=probarray)[0]
                    weight_index=np.unravel_index(sample_indices, (len(index)))

                    inj_muon_energies=[gridpts[weight_index]]
                    energysample=True
                else:
                    print('nans in probarray??')
                    inj_muon_energies=default_sampler(multiplicity=multiplicity)
                    energysample=False
            else:#use default energy from powerlaw
                inj_muon_energies=default_sampler(multiplicity=multiplicity)
                energysample=False
        print('EnergySample',energysample)
        return inj_muon_energies,energysample

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
        
    def _create_bundle_mc_tree(self, inj_energy,inj_pos, inj_time, inj_dir,depth,multiplicity):
        """Inject accompanying muon in provided I3MCTree

        Parameters
        ----------
        inj_energy : Primary energy 
            Neutrino Energy.
        inj_pos : I3Position
            The position at which to inject the muons.
        inj_time : double
            The time at which to inject the muons.
        inj_dir : I3Direction
            The direction of the injected muons.
        multiplicity : int
            The number of injected muons.

        Returns
        -------
        I3MCTree
            The modified I3MCTree with the injected Muon.
        I3MapStringDouble
            Information on the injected muon.
        """
        mc_tree = dataclasses.I3MCTree()
        print('Empty MC Tree created')
        injectdict={}
        keys=['muon_1_energy','muon_1_pdg_encoding','muon_2_energy','muon_2_pdg_encoding','muon_3_energy','muon_3_pdg_encoding','muon_4_energy','muon_4_pdg_encoding']
        for key in keys:#initializing injection dictionary to zero
            injectdict[key]=0
        print('Starting Sampler')
        inj_muon_energies,energysample=self._sample_bundle_energies(multiplicity,inj_energy,np.cos(inj_dir.zenith),depth)
        #len(inj_muon_energies) is the injected multiplicity
        print('Energies',inj_muon_energies)
        ctr=1
        if not inj_muon_energies:
            injection_info=dataclasses.I3MapStringDouble(injectdict)
            return mc_tree, injection_info,energysample
        for muon_energy in inj_muon_energies:

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
            # muon.energy = self._sample_energy()
            
            muon.energy = muon_energy

            # add muon primary to I3MCTree
            mc_tree.add_primary(muon_primary)

            # add muon as new child
            mc_tree.append_child(muon_primary, muon)

            # add info
            keys=['muon_'+str(ctr)+'_energy','muon_'+str(ctr)+'_pdg_encoding']
            ctr+=1
            injectdict[keys[0]]=muon.energy
            injectdict[keys[1]]=muon.pdg_encoding
            # injectdict={
            #     'muon_energy': muon.energy,
            #     'muon_pdg_encoding': muon.pdg_encoding,
            # }
        injection_info=dataclasses.I3MapStringDouble(injectdict)
        return mc_tree, injection_info,energysample


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
