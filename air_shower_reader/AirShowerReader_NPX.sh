#!/bin/bash

#prep inputs
infile=$1
outfile=$2
#setup env
source /home/vbasu/test/testdev_sh.sh
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`
export HDF5_USE_FILE_LOCKING='FALSE'
#run job                                              
#python /home/vbasu/scripts/MCPropagator.py -i /data/sim/IceCube/2016/filtered/level2/neutrino-generator/21218/0000000-0000999/Level2_IC86.2016_NuE.021218.000999.i3.zst -o /home/vbasu/testoutput/L2_propagated_condor.i3.zst -g /cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz
/cvmfs/icecube.opensciencegrid.org/users/vbasu/meta-projects/combo3/build/env-shell.sh python /home/zrechav/scripts/air_shower_reader/AirShowerReader.py -i ${infile} -o ${outfile}  
