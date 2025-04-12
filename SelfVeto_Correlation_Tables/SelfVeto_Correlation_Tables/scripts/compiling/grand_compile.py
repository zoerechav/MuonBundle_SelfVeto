#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/icetray-start
#METAPROJECT /home/zrechav/i3/icetray/build
import pyarrow
from icecube.icetray import I3Tray
from icecube import icetray, phys_services, dataio, dataclasses,MuonGun
import simweights
import pandas as pd
from pathlib import Path
import simweights
import h5py
import numpy as np
import traceback


import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", type=str, default = "/data/user/zrechav/compiled_hdf5s/airshowered_corsika/22803/", dest = "INPUT", help = "input i3 file directory")

parser.add_argument("-o", "--output", action="store", type=str, default ='/data/user/zrechav/compiled_hdf5s/airshowered_corsika/correlation/22803_0000000-0003999_correlation.hdf5', dest = "OUTPUT", help = "Output i3 file")




args = parser.parse_args()

print(args)

input_path = Path(args.INPUT)

input_files = sorted(str(f) for f in input_path.glob("*.hdf5"))
#input_files = sorted(glob.glob(input_path + '*_weighted.hdf5'))
print(input_files)
output_file = args.OUTPUT


# with h5py.File(output_file, 'w') as combined_file:
#     for file in input_files:
#         #file = (file)
#         print(file)
#         with h5py.File(file, 'r') as f:
#             # Iterate over all keys (datasets) in the current file
#             for key in f.keys():
#                 # Copy the dataset to the combined file
#                 f.copy(key, combined_file)


table_keys = [
    "CorsikaWeightMap", 
    "I3EventHeader", 
    "MuonMultiplicity",
    "Muon_Energy_L1",
    "Muon_Energy_L2",
    "Muon_Energy_L3",
    "Muon_Energy_L4",
    "Muon_L1", 
    "Muon_L1_Depth",    
    "Muon_L1_z_inter", 
    "Muon_L2",
    "Muon_L2_Depth",  
    "Muon_L2_z_inter", 
    "Muon_L3",
    "Muon_L3_Depth",    
    "Muon_L3_z_inter",
    "Muon_L4",
    "Muon_L4_Depth",  
    "Muon_L4_z_inter", 
    "NeutrinoParent", 
    "PolyplopiaPrimary",
    "PrimaryMass",
    "ShowerNeutrino", 
    "Total_Muon_Energy",
    "num_neutrinos",
    "shower_neutrino_depth", 
    "shower_neutrino_energy",
    "shower_neutrino_type", 
    "shower_neutrino_zenith",
    "total_neutrino_energy",
    "x_inter",
    "y_inter", 
    "z_inter"
]

with h5py.File(output_file, 'w') as combined_file:
    for file in input_files:
        print(f'Processing: {file}')
        with h5py.File(file, 'r') as f:
            # Iterate over the specified keys in table_keys
            for key in table_keys:
                if key in f:  # Check if the key exists in the current file
                    data = f[key][:]
                    
                    if key in combined_file:
                        # If the dataset already exists, resize and concatenate
                        combined_file[key].resize((combined_file[key].shape[0] + data.shape[0]), axis=0)
                        combined_file[key][-data.shape[0]:] = data
                    else:
                        # If the dataset does not exist, create it
                        combined_file.create_dataset(key, data=data, maxshape=(None,), chunks=True)
                         #combined_file.create_dataset(key, data=data, maxshape=(None,), chunks=True, compression='gzip', compression_opts=9)
print(f'Combined HDF5 file created: {output_file}')
#print(f'Combined {len(hdf5_files)} files into {output_file}')

print("done")
