#!/bin/bash

#prep inputs
infile=$1
outfile=$2
#setup env
source /home/vbasu/test/testdev_sh.sh
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`
export HDF5_USE_FILE_LOCKING='FALSE'

/cvmfs/icecube.opensciencegrid.org/users/vbasu/meta-projects/combo3/build/env-shell.sh python /home/zrechav/SelfVeto_Correlation_Tables/scripts/air_shower_reader/AirShowerReader.py -i ${infile} -o ${outfile}  
