#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/icetray-start
#METAPROJECT combo/V01-00-00
import click
import yaml

import numpy as np

from icecube.simprod import segments

from I3Tray import I3Tray
from icecube import icetray, dataclasses
from icecube import sim_services, MuonGun

from utils import create_random_services, get_run_folder
from resources.biased_simulation import BaseSimulationBias
from dom_distance_cut import OversizeSplitterNSplits, generate_stream_object


@click.command()
@click.argument('cfg', type=click.Path(exists=True))
@click.argument('run_number', type=int)
@click.option('--scratch/--no-scratch', default=True)
def main(cfg, run_number, scratch):
    with open(cfg, 'r') as stream:
        if int(yaml.__version__[0]) < 5:
            # backwards compatibility for yaml versions before version 5
            cfg = yaml.load(stream)
        else:
            cfg = yaml.full_load(stream)
    cfg['run_number'] = run_number
    cfg['run_folder'] = get_run_folder(run_number)
    if scratch:
        outfile = cfg['scratchfile_pattern'].format(**cfg)
    else:
        outfile = cfg['outfile_pattern'].format(**cfg)
    outfile = outfile.replace(' ', '0')
    click.echo('NEvents: {}'.format(cfg['n_events_per_run']))
    click.echo('EMin: {}'.format(cfg['e_min']))
    click.echo('EMax: {}'.format(cfg['e_max']))
    click.echo('Gamma: {}'.format(cfg['gamma']))
    click.echo('ZenithMin: {}'.format(cfg['zenith_min']))
    click.echo('ZenithMax: {}'.format(cfg['zenith_max']))
    click.echo('AzimuthMin: {}'.format(cfg['azimuth_min']))
    click.echo('AzimuthMax: {}'.format(cfg['azimuth_max']))
    if cfg['neutrino_flavor'] is None:
        click.echo('NeutrinoTypes: {}'.format(cfg['neutrino_types']))
        click.echo('PrimaryTypeRatio: {}'.format(cfg['primary_type_ratio']))
    else:
        click.echo('NeutrinoFlavor: {}'.format(cfg['neutrino_flavor']))
    if 'ApplyBaseSimulationBias' in cfg and cfg['ApplyBaseSimulationBias']:
        click.echo('Apply simulation bias: True')
    else:
        click.echo('Apply simulation bias: True')

    tray = I3Tray()

    if 'ApplyBaseSimulationBias' in cfg and cfg['ApplyBaseSimulationBias']:
        n_services = 3
    else:
        n_services = 2
    random_services, _ = create_random_services(
        dataset_number=cfg['dataset_number'],
        run_number=cfg['run_number'],
        seed=cfg['seed'],
        n_services=n_services,
        use_gslrng=cfg['random_service_use_gslrng'],
    )

    random_service, random_service_prop = random_services[:2]
    tray.context['I3RandomService'] = random_service

    tray.AddModule("I3InfiniteSource",
                   "TheSource",
                   Prefix=cfg['gcd'],
                   Stream=icetray.I3Frame.DAQ)

    tray.AddSegment(
        segments.GenerateNeutrinos, "GenerateNeutrinos",
        RunID=run_number,
        RandomService=random_service,
        NumEvents=cfg['n_events_per_run'],
        FromEnergy=cfg['e_min']*icetray.I3Units.GeV,
        ToEnergy=cfg['e_max']*icetray.I3Units.GeV,
        ZenithRange=[cfg['zenith_min'] * icetray.I3Units.deg,
                     cfg['zenith_max'] * icetray.I3Units.deg],
        AzimuthRange=[cfg['azimuth_min'] * icetray.I3Units.deg,
                      cfg['azimuth_max'] * icetray.I3Units.deg],
        **cfg['additional_GenerateNeutrinos_settings']
    )

    # propagate muons if config exists in config
    # Note: Snowstorm may perform muon propagation internally
    if 'muon_propagation_config' in cfg:
        tray.AddSegment(segments.PropagateMuons,
                        'propagate_muons',
                        RandomService=random_service_prop,
                        **cfg['muon_propagation_config'])
    else:
        # In this case we are not propagating the I3MCTree yet, but
        # are letting this be done by snowstorm propagation
        # We need to add a key named 'I3MCTree', since snowstorm expects this
        # It will propagate the particles for us.
        tray.AddModule('Rename', keys=['I3MCTree_preMuonProp', 'I3MCTree'])

    # Bias simulation if desired
    if 'ApplyBaseSimulationBias' in cfg and cfg['ApplyBaseSimulationBias']:
        tray.AddModule(
            BaseSimulationBias,
            'BaseSimulationBias',
            random_service=random_services[2],
            **cfg['BaseSimulationBiasSettings']
        )

    if cfg['distance_splits'] is not None:
        import dom_distance_cut as dom_cut
        click.echo('Oversizestreams')
        stream_objects = dom_cut.generate_stream_object(
            cut_distances=cfg['distance_splits'],
            dom_limits=cfg['threshold_doms'],
            oversize_factors=cfg['oversize_factors'])
        tray.AddModule(dom_cut.OversizeSplitterNSplits,
                       "OversizeSplitterNSplits",
                       thresholds=cfg['distance_splits'],
                       thresholds_doms=cfg['threshold_doms'],
                       oversize_factors=cfg['oversize_factors'],
                       simulaton_type=cfg['neutrino_flavor'].lower())
        for stream_i in stream_objects:
            outfile_i = stream_i.transform_filepath(outfile)
            click.echo('\t{}'.format(stream_i))
            click.echo('\tOutfile: {}'.format(outfile_i))
            tray.AddModule("I3Writer",
                           "writer_{}".format(stream_i.stream_name),
                           Filename=outfile_i,
                           Streams=[icetray.I3Frame.DAQ,
                                    icetray.I3Frame.Physics,
                                    icetray.I3Frame.Stream('S'),
                                    icetray.I3Frame.Stream('M')],
                           If=stream_i)
    else:
        click.echo('Output: {}'.format(outfile))
        tray.AddModule("I3Writer", "writer",
                       Filename=outfile,
                       Streams=[icetray.I3Frame.DAQ,
                                icetray.I3Frame.Physics,
                                icetray.I3Frame.Stream('S'),
                                icetray.I3Frame.Stream('M')])
    click.echo('Scratch: {}'.format(scratch))
    tray.AddModule("TrashCan", "the can")
    tray.Execute()
    tray.Finish()


if __name__ == '__main__':
    main()
