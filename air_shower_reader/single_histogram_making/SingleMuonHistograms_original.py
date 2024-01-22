#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/vbasu/meta-projects/combo2/build
import numpy as np
import math
import os.path
from os import path

import tables
import h5py
import numpy as np
import math
import os.path
from os import path

import tables
import glob
import h5py
import pandas as pd


import argparse
import os
 
# handling of command line arguments  
#from optparse import OptionParser
parser = argparse.ArgumentParser()
#usage = """%prog [options]"""
#parser.set_usage(usage)


parser.add_argument("-f", "--flavour", action="store", type=str, default='NuE', dest="FLAVOUR", help="Flavour")
parser.add_argument("-a", "--angle", action="store", type=float, default=0.20, dest="ANGLE", help="Input angle")
parser.add_argument("-d", "--depth", action="store", type=float, default=1.75, dest="DEPTH", help="Inputdepth")


# parse cmd line args, bail out if anything is not understood
args = parser.parse_args()

coszen=np.around(args.ANGLE,decimals=2)
depth=np.around(args.DEPTH,decimals=2)
flavour=args.FLAVOUR
##Parameter in file name is lower boundary of bin range (ex: 0.04-0.08 cos(zenith))
SingleMuonMask=np.load('/home/zrechav/scripts/air_shower_reader/picklemaking/masked_pickles/SingleMuonMask.npy',mmap_mode='r',allow_pickle=True)
DoubleMuonMask = np.load('/home/zrechav/scripts/air_shower_reader/picklemaking/masked_pickles/DoubleMuonMask.npy', mmap_mode='r',allow_pickle=True)
angle_mask=np.load('/home/zrechav/scripts/air_shower_reader/picklemaking/masked_pickles/Angle_mask_'+str(coszen)+'.npy',mmap_mode='r',allow_pickle=True)
depth_mask=np.load('/home/zrechav/scripts/air_shower_reader/picklemaking/masked_pickles/Depth_mask_'+str(depth)+'.npy',mmap_mode='r',allow_pickle=True)
Flav_mask=np.load('/home/zrechav/scripts/air_shower_reader/picklemaking/masked_pickles/Flavour_mask_'+str(flavour)+'.npy',mmap_mode='r',allow_pickle=True)
# print('Double Muons:',np.sum(DoubleMuonMask))

E_nu_bins=np.logspace(1,8,10+1)
nu_bin_centers = np.sqrt(E_nu_bins[:-1] * E_nu_bins[1:])
E_mu_bins=np.logspace(-1,8,10+1)
mu_bin_centers = np.sqrt(E_mu_bins[:-1] * E_mu_bins[1:])
##Correlated Muon Energies (complete)
Muon_Corr_Energy_L1, Muon_Corr_Energy_L2, Muon_Corr_Energy_L3, Muon_Corr_Energy_L4, Muon_Corr_Energy_L5 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
##Correlated Muon Depths (complete)
Muon_Corr_Depth_L1, Muon_Corr_Depth_L2, Muon_Corr_Depth_L3, Muon_Corr_Depth_L4, Muon_Corr_Depth_L5 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
##Correlated Muon Zeniths (complete)
Muon_Corr_Zenith_L1, Muon_Corr_Zenith_L2, Muon_Corr_Zenith_L3, Muon_Corr_Zenith_L4, Muon_Corr_Zenith_L5, = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
##Correlated Neutrino Separation (complete)
Muon_Corr_Neutrino_separation_L1, Muon_Corr_Neutrino_separation_L2, Muon_Corr_Neutrino_separation_L3, Muon_Corr_Neutrino_separation_L4, Muon_Corr_Neutrino_separation_L5 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
##Correlated Shower Separation (complete)
Muon_Corr_Shower_separation_L1, Muon_Corr_Shower_separation_L2, Muon_Corr_Shower_separation_L3, Muon_Corr_Shower_separation_L4, Muon_Corr_Shower_separation_L5 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
##Correlated Muon Multiplicities (complete)
MuonMultiplicity_Corr, MuonMultiplicity_Corr_10, MuonMultiplicity_Corr_100, MuonMultiplicity_Corr_200 = np.array([]), np.array([]), np.array([]), np.array([])

##Uncorrelated Muon Information -- THIS IS WHAT I NEED
Muon_UNCorr_Energy_L1, Muon_UNCorr_Energy_L2, Muon_UNCorr_Energy_L3, Muon_UNCorr_Energy_L4, Muon_UNCorr_Energy_L5 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
##Correlated Muon Depths (complete)
Muon_UNCorr_Depth_L1, Muon_UNCorr_Depth_L2, Muon_UNCorr_Depth_L3, Muon_UNCorr_Depth_L4, Muon_UNCorr_Depth_L5 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
##Correlated Muon Zeniths (complete)
Muon_UNCorr_Zenith_L1, Muon_UNCorr_Zenith_L2, Muon_UNCorr_Zenith_L3, Muon_UNCorr_Zenith_L4, Muon_UNCorr_Zenith_L5, = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
##Correlated Neutrino Separation (complete)
Muon_UNCorr_Neutrino_separation_L1, Muon_UNCorr_Neutrino_separation_L2, Muon_UNCorr_Neutrino_separation_L3, Muon_UNCorr_Neutrino_separation_L4, Muon_UNCorr_Neutrino_separation_L5 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
##Correlated Shower Separation (complete)
Muon_UNCorr_Shower_separation_L1, Muon_UNCorr_Shower_separation_L2, Muon_UNCorr_Shower_separation_L3, Muon_UNCorr_Shower_separation_L4, Muon_UNCorr_Shower_separation_L5 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
##Correlated Muon Multiplicities (complete)
MuonMultiplicity_UNCorr, MuonMultiplicity_UNCorr_10, MuonMultiplicity_UNCorr_100, MuonMultiplicity_UNCorr_200 = np.array([]), np.array([]), np.array([]), np.array([])

##Muon Energies (complete)
Muon_Energy_L1, Muon_Energy_L2, Muon_Energy_L3, Muon_Energy_L4, Muon_Energy_L5 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
##Muon Depths (complete)
Muon_Depth_L1, Muon_Depth_L2, Muon_Depth_L3, Muon_Depth_L4, Muon_Depth_L5 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
##Muon Zeniths (complete)
Muon_Zenith_L1, Muon_Zenith_L2, Muon_Zenith_L3, Muon_Zenith_L4, Muon_Zenith_L5, = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
##Neutrino Separation (complete)
Muon_Neutrino_separation_L1, Muon_Neutrino_separation_L2, Muon_Neutrino_separation_L3, Muon_Neutrino_separation_L4, Muon_Neutrino_separation_L5 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
##Shower Separation (complete)
Muon_Shower_separation_L1, Muon_Shower_separation_L2, Muon_Shower_separation_L3, Muon_Shower_separation_L4, Muon_Shower_separation_L5 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
##Muon Multiplicities (complete)
MuonMultiplicity, MuonMultiplicity_10, MuonMultiplicity_100, MuonMultiplicity_200 = np.array([]), np.array([]), np.array([]), np.array([])


Max_Deposited_Energy = np.array([])     # (complete)
Energy_Shower_Neutrino = np.array([])   # (complete)
Zenith_Shower_Neutrino = np.array([])   # (complete)
Depth_Shower_Neutrino = np.array([])    # (complete)


Energy_Brightest_Muon = np.array([])   # (complete)
Energy_BrightestMuonAtDepth = np.array([]) # (complete)
Depth_BrightestMuon = np.array([])     # (complete)
Zenith_BrightestMuon = np.array([])    # (complete)


Total_Muon_energy = np.array([])
Total_Corr_Muon_energy = np.array([])
Total_UNCorr_Muon_energy = np.array([])

Flavour_Shower_Neutrino = np.array([])
GaisserH4a_weight = np.array([])
PrimaryMass, PrimaryEnergyperNucleon= np.array([]),np.array([])
filename = '/data/user/zrechav/compiled_hdf5s/airshowered_corsika/Corsika_20904_sample.hdf5'
#filename='/data/user/zrechav/airshowered_corsika/20904/0000000-0000999/airshowered_corsika_20904_000000.hdf5'
#filename = '/home/zrechav/test/AirShowerReader.hdf5'
print(filename)
# filename="/home/vbasu/testoutput/CorsikaShowerReader2.hdf5"
if (path.exists(filename)):
    try:
        print('I EXIST')
        hdf=h5py.File(filename, 'r')
########ENERGY, ZENITH, AND WEIGHTING INFORMATION#######
        GaisserH4a_weight = np.append(GaisserH4a_weight,np.asarray(hdf.get('GaisserH4a_weight')['item']))
        Zenith_Shower_Neutrino = np.append(Zenith_Shower_Neutrino,np.asarray(hdf.get('shower_neutrino_zenith')['value']))
        Flavour_Shower_Neutrino = np.append(Flavour_Shower_Neutrino,np.asarray(hdf.get('shower_neutrino_type')['value']))
        Energy_Shower_Neutrino = np.append(Energy_Shower_Neutrino,np.asarray(hdf.get('shower_neutrino_energy')['value']))
        Depth_Shower_Neutrino = np.append(Depth_Shower_Neutrino,np.asarray(hdf.get('shower_neutrino_depth')['value']))
        
        Max_Deposited_Energy = np.append(Max_Deposited_Energy, np.asarray(hdf.get('max_depo_energy')['value']))
        Energy_Brightest_Muon = np.append(Energy_Brightest_Muon,np.asarray(hdf.get('BrightestMuon')['energy']))
        Energy_BrightestMuonAtDepth = np.append(Energy_BrightestMuonAtDepth , np.asarray(hdf.get('Energy_BrightestMuonAtDepth')['value']))
        Depth_BrightestMuon = np.append(Depth_BrightestMuon, np.asarray(hdf.get('BrightestMuonDepth')['value']))
        Zenith_BrightestMuon = np.append(Zenith_BrightestMuon, np.asarray(hdf.get('BrightestMuon')['zenith']))
########MUON INFORMATION: MAY OR MAY NOT BE CORRELATED#######      
        MuonMultiplicity = np.append(MuonMultiplicity,np.asarray(hdf.get('MuonMultiplicity')['value']))
        
        MuonMultiplicity_10 = np.append(MuonMultiplicity_10,np.asarray(hdf.get('MuonMultiplicity_10')['value']))
        MuonMultiplicity_100 = np.append(MuonMultiplicity_100,np.asarray(hdf.get('MuonMultiplicity_100')['value']))
        MuonMultiplicity_200 = np.append(MuonMultiplicity_200,np.asarray(hdf.get('MuonMultiplicity_200')['value']))
        
        Muon_Energy_L1=np.append(Muon_Energy_L1,np.asarray(hdf.get('Muon_Energy_L1')['value']))
        Muon_Energy_L2=np.append(Muon_Energy_L2,np.asarray(hdf.get('Muon_Energy_L2')['value']))
        Muon_Energy_L3=np.append(Muon_Energy_L3,np.asarray(hdf.get('Muon_Energy_L3')['value']))
        Muon_Energy_L4=np.append(Muon_Energy_L4,np.asarray(hdf.get('Muon_Energy_L4')['value']))
        Muon_Energy_L5=np.append(Muon_Energy_L5,np.asarray(hdf.get('Muon_Energy_L5')['value']))
        print('DID I MAKE IT HERE?')
        Muon_Depth_L1=np.append(Muon_Depth_L1,np.asarray(hdf.get('Muon_L1_Depth')['value']))
        Muon_Depth_L2=np.append(Muon_Depth_L2,np.asarray(hdf.get('Muon_L2_Depth')['value']))
        Muon_Depth_L3=np.append(Muon_Depth_L3,np.asarray(hdf.get('Muon_L3_Depth')['value']))
        Muon_Depth_L4=np.append(Muon_Depth_L4,np.asarray(hdf.get('Muon_L4_Depth')['value']))
        Muon_Depth_L5=np.append(Muon_Depth_L5,np.asarray(hdf.get('Muon_L5_Depth')['value']))
        
        Muon_Zenith_L1=np.append(Muon_Zenith_L1,np.asarray(hdf.get('Muon_L1')['zenith']))
        Muon_Zenith_L2=np.append(Muon_Zenith_L2,np.asarray(hdf.get('Muon_L2')['zenith']))
        Muon_Zenith_L3=np.append(Muon_Zenith_L3,np.asarray(hdf.get('Muon_L3')['zenith']))
        Muon_Zenith_L4=np.append(Muon_Zenith_L4,np.asarray(hdf.get('Muon_L4')['zenith']))
        Muon_Zenith_L5=np.append(Muon_Zenith_L5,np.asarray(hdf.get('Muon_L5')['zenith']))
        print('DOUBLE CHECKING THIS IS INDEED THE FIX')
        
        Muon_Neutrino_separation_L1 = np.append(Muon_Neutrino_separation_L1,np.asarray(hdf.get('Muon_L1_Neutrino_separation')['value']))
        Muon_Neutrino_separation_L2 = np.append(Muon_Neutrino_separation_L2,np.asarray(hdf.get('Muon_L2_Neutrino_separation')['value']))
        Muon_Neutrino_separation_L3 = np.append(Muon_Neutrino_separation_L3,np.asarray(hdf.get('Muon_L3_Neutrino_separation')['value']))
        Muon_Neutrino_separation_L4 = np.append(Muon_Neutrino_separation_L4,np.asarray(hdf.get('Muon_L4_Neutrino_separation')['value']))
        Muon_Neutrino_separation_L5 = np.append(Muon_Neutrino_separation_L5,np.asarray(hdf.get('Muon_L5_Neutrino_separation')['value']))
        print('NOW WHAT')
        Muon_Shower_separation_L1 = np.append(Muon_Shower_separation_L1,np.asarray(hdf.get('Muon_L1_Shower_separation')['value']))
        Muon_Shower_separation_L2 = np.append(Muon_Shower_separation_L2,np.asarray(hdf.get('Muon_L2_Shower_separation')['value']))
        Muon_Shower_separation_L3 = np.append(Muon_Shower_separation_L3,np.asarray(hdf.get('Muon_L3_Shower_separation')['value']))
        Muon_Shower_separation_L4 = np.append(Muon_Shower_separation_L4,np.asarray(hdf.get('Muon_L4_Shower_separation')['value']))
        Muon_Shower_separation_L5 = np.append(Muon_Shower_separation_L5,np.asarray(hdf.get('Muon_L5_Shower_separation')['value']))
        print('MEH')
        Total_Muon_energy = np.append(Total_Muon_energy,np.asarray(hdf.get('Total_Muon_energy')['value']))
########CORRELATED MUON INFORMATION############
        print('HERRRRRRRRRE')
        MuonMultiplicity_Corr = np.append(MuonMultiplicity_Corr,np.asarray(hdf.get('MuonMultiplicity_Corr')['value']))
        MuonMultiplicity_Corr_10 = np.append(MuonMultiplicity_Corr_10,np.asarray(hdf.get('MuonMultiplicity_Corr_10')['value']))
        MuonMultiplicity_Corr_100 = np.append(MuonMultiplicity_Corr_100,np.asarray(hdf.get('MuonMultiplicity_Corr_100')['value']))
        MuonMultiplicity_Corr_200 = np.append(MuonMultiplicity_Corr_200,np.asarray(hdf.get('MuonMultiplicity_Corr_200')['value']))
        print('EEEEK')
        Muon_Corr_Energy_L1 = np.append(Muon_Corr_Energy_L1, np.asarray(hdf.get('Muon_Corr_Energy_L1')['value']))
        Muon_Corr_Energy_L2 = np.append(Muon_Corr_Energy_L2, np.asarray(hdf.get('Muon_Corr_Energy_L2')['value']))
        Muon_Corr_Energy_L3 = np.append(Muon_Corr_Energy_L3, np.asarray(hdf.get('Muon_Corr_Energy_L3')['value']))
        Muon_Corr_Energy_L4 = np.append(Muon_Corr_Energy_L4, np.asarray(hdf.get('Muon_Corr_Energy_L4')['value']))
        Muon_Corr_Energy_L5 = np.append(Muon_Corr_Energy_L5, np.asarray(hdf.get('Muon_Corr_Energy_L5')['value']))
        
        Muon_Corr_Depth_L1 = np.append(Muon_Corr_Depth_L1, np.asarray(hdf.get('Muon_Corr_L1_Depth')['value']))
        Muon_Corr_Depth_L2 = np.append(Muon_Corr_Depth_L2, np.asarray(hdf.get('Muon_Corr_L2_Depth')['value']))
        Muon_Corr_Depth_L3 = np.append(Muon_Corr_Depth_L3, np.asarray(hdf.get('Muon_Corr_L3_Depth')['value']))
        Muon_Corr_Depth_L4 = np.append(Muon_Corr_Depth_L4, np.asarray(hdf.get('Muon_Corr_L4_Depth')['value']))
        Muon_Corr_Depth_L5 = np.append(Muon_Corr_Depth_L5, np.asarray(hdf.get('Muon_Corr_L5_Depth')['value']))
        print('HOOPLAH')
        Muon_Corr_Zenith_L1=np.append(Muon_Corr_Zenith_L1,np.asarray(hdf.get('Muon_Corr_L1')['zenith']))
        Muon_Corr_Zenith_L2=np.append(Muon_Corr_Zenith_L2,np.asarray(hdf.get('Muon_Corr_L2')['zenith']))
        Muon_Corr_Zenith_L3=np.append(Muon_Corr_Zenith_L3,np.asarray(hdf.get('Muon_Corr_L3')['zenith']))
        Muon_Corr_Zenith_L4=np.append(Muon_Corr_Zenith_L4,np.asarray(hdf.get('Muon_Corr_L4')['zenith']))
        Muon_Corr_Zenith_L5=np.append(Muon_Corr_Zenith_L5,np.asarray(hdf.get('Muon_Corr_L5')['zenith']))
        
        Muon_Corr_Neutrino_separation_L1 = np.append(Muon_Corr_Neutrino_separation_L1, np.asarray(hdf.get('Muon_Corr_L1_Neutrino_separation')['value']))
        Muon_Corr_Neutrino_separation_L2 = np.append(Muon_Corr_Neutrino_separation_L2, np.asarray(hdf.get('Muon_Corr_L2_Neutrino_separation')['value']))
        Muon_Corr_Neutrino_separation_L3 = np.append(Muon_Corr_Neutrino_separation_L3, np.asarray(hdf.get('Muon_Corr_L3_Neutrino_separation')['value']))
        Muon_Corr_Neutrino_separation_L4 = np.append(Muon_Corr_Neutrino_separation_L4, np.asarray(hdf.get('Muon_Corr_L4_Neutrino_separation')['value']))
        Muon_Corr_Neutrino_separation_L5 = np.append(Muon_Corr_Neutrino_separation_L5, np.asarray(hdf.get('Muon_Corr_L5_Neutrino_separation')['value']))
        
        Muon_Corr_Shower_separation_L1 = np.append(Muon_Corr_Shower_separation_L1, np.asarray(hdf.get('Muon_Corr_L1_Shower_separation')['value']))
        Muon_Corr_Shower_separation_L2 = np.append(Muon_Corr_Shower_separation_L2, np.asarray(hdf.get('Muon_Corr_L2_Shower_separation')['value']))
        Muon_Corr_Shower_separation_L3 = np.append(Muon_Corr_Shower_separation_L3, np.asarray(hdf.get('Muon_Corr_L3_Shower_separation')['value']))
        Muon_Corr_Shower_separation_L4 = np.append(Muon_Corr_Shower_separation_L4, np.asarray(hdf.get('Muon_Corr_L4_Shower_separation')['value']))
        Muon_Corr_Shower_separation_L5 = np.append(Muon_Corr_Shower_separation_L5, np.asarray(hdf.get('Muon_Corr_L5_Shower_separation')['value']))
        print('FINAL COUNTDOWN')
        Total_Corr_Muon_energy = np.append(Total_Corr_Muon_energy,np.asarray(hdf.get('Total_Corr_Muon_energy')['value']))
########UNCORRELATED MUON INFORMATION############
        print("AM I ALL THE WAY DOWN HERE?")
        MuonMultiplicity_UNCorr = np.append(MuonMultiplicity_UNCorr,np.asarray(hdf.get('MuonMultiplicity_UnCorr')['value']))
        MuonMultiplicity_UNCorr_10 = np.append(MuonMultiplicity_UNCorr_10,np.asarray(hdf.get('MuonMultiplicity_UnCorr_10')['value']))
        MuonMultiplicity_UNCorr_100 = np.append(MuonMultiplicity_UNCorr_100,np.asarray(hdf.get('MuonMultiplicity_UnCorr_100')['value']))
        MuonMultiplicity_UNCorr_200 = np.append(MuonMultiplicity_UNCorr_200,np.asarray(hdf.get('MuonMultiplicity_UnCorr_200')['value']))

        Muon_UNCorr_Energy_L1 = np.append(Muon_UNCorr_Energy_L1, np.asarray(hdf.get('Muon_UnCorr_Energy_L1')['value']))
        Muon_UNCorr_Energy_L2 = np.append(Muon_UNCorr_Energy_L2, np.asarray(hdf.get('Muon_UnCorr_Energy_L2')['value']))
        Muon_UNCorr_Energy_L3 = np.append(Muon_UNCorr_Energy_L3, np.asarray(hdf.get('Muon_UnCorr_Energy_L3')['value']))
        Muon_UNCorr_Energy_L4 = np.append(Muon_UNCorr_Energy_L4, np.asarray(hdf.get('Muon_UnCorr_Energy_L4')['value']))
        Muon_UNCorr_Energy_L5 = np.append(Muon_UNCorr_Energy_L5, np.asarray(hdf.get('Muon_UnCorr_Energy_L5')['value']))

        Muon_UNCorr_Depth_L1 = np.append(Muon_UNCorr_Depth_L1, np.asarray(hdf.get('Muon_UnCorr_L1_Depth')['value']))
        Muon_UNCorr_Depth_L2 = np.append(Muon_UNCorr_Depth_L2, np.asarray(hdf.get('Muon_UnCorr_L2_Depth')['value']))
        Muon_UNCorr_Depth_L3 = np.append(Muon_UNCorr_Depth_L3, np.asarray(hdf.get('Muon_UnCorr_L3_Depth')['value']))
        Muon_UNCorr_Depth_L4 = np.append(Muon_UNCorr_Depth_L4, np.asarray(hdf.get('Muon_UnCorr_L4_Depth')['value']))
        Muon_UNCorr_Depth_L5 = np.append(Muon_UNCorr_Depth_L5, np.asarray(hdf.get('Muon_UnCorr_L5_Depth')['value']))

        Muon_UNCorr_Zenith_L1=np.append(Muon_UNCorr_Zenith_L1,np.asarray(hdf.get('Muon_UnCorr_L1')['zenith']))
        Muon_UNCorr_Zenith_L2=np.append(Muon_UNCorr_Zenith_L2,np.asarray(hdf.get('Muon_UnCorr_L2')['zenith']))
        Muon_UNCorr_Zenith_L3=np.append(Muon_UNCorr_Zenith_L3,np.asarray(hdf.get('Muon_UnCorr_L3')['zenith']))
        Muon_UNCorr_Zenith_L4=np.append(Muon_UNCorr_Zenith_L4,np.asarray(hdf.get('Muon_UnCorr_L4')['zenith']))
        Muon_UNCorr_Zenith_L5=np.append(Muon_UNCorr_Zenith_L5,np.asarray(hdf.get('Muon_UnCorr_L5')['zenith']))

        Muon_UNCorr_Neutrino_separation_L1 = np.append(Muon_UNCorr_Neutrino_separation_L1, np.asarray(hdf.get('Muon_UnCorr_L1_Neutrino_separation')['value']))
        Muon_UNCorr_Neutrino_separation_L2 = np.append(Muon_UNCorr_Neutrino_separation_L2, np.asarray(hdf.get('Muon_UnCorr_L2_Neutrino_separation')['value']))
        Muon_UNCorr_Neutrino_separation_L3 = np.append(Muon_UNCorr_Neutrino_separation_L3, np.asarray(hdf.get('Muon_UnCorr_L3_Neutrino_separation')['value']))
        Muon_UNCorr_Neutrino_separation_L4 = np.append(Muon_UNCorr_Neutrino_separation_L4, np.asarray(hdf.get('Muon_UnCorr_L4_Neutrino_separation')['value']))
        Muon_UNCorr_Neutrino_separation_L5 = np.append(Muon_UNCorr_Neutrino_separation_L5, np.asarray(hdf.get('Muon_UnCorr_L5_Neutrino_separation')['value']))

        Muon_UNCorr_Shower_separation_L1 = np.append(Muon_UNCorr_Shower_separation_L1, np.asarray(hdf.get('Muon_UnCorr_L1_Shower_separation')['value']))
        Muon_UNCorr_Shower_separation_L2 = np.append(Muon_UNCorr_Shower_separation_L2, np.asarray(hdf.get('Muon_UnCorr_L2_Shower_separation')['value']))
        Muon_UNCorr_Shower_separation_L3 = np.append(Muon_UNCorr_Shower_separation_L3, np.asarray(hdf.get('Muon_UnCorr_L3_Shower_separation')['value']))
        Muon_UNCorr_Shower_separation_L4 = np.append(Muon_UNCorr_Shower_separation_L4, np.asarray(hdf.get('Muon_UnCorr_L4_Shower_separation')['value']))
        Muon_UNCorr_Shower_separation_L5 = np.append(Muon_UNCorr_Shower_separation_L5, np.asarray(hdf.get('Muon_UnCorr_L5_Shower_separation')['value']))
        print('WHAT THE BISCUITS')
        Total_UNCorr_Muon_energy = np.append(Total_UNCorr_Muon_energy,np.asarray(hdf.get('Total_UnCorr_Muon_energy')['value']))
        print('I FINISHED LOADING MY ARRAYS WITH HDF INFO')
        hdf.close()
    except Exception as e:
        print(filename+" Is Faulty")
        print(e)

        hdf.close()


MuonEnergyThreshMask=[energy>10 for energy in Muon_Energy_L1]

Energy_Brightest_Muon=np.float32(Energy_Brightest_Muon[MuonEnergyThreshMask])
Energy_Shower_Neutrino=np.float32(Energy_Shower_Neutrino[MuonEnergyThreshMask])
Flavour_Shower_Neutrino=np.float32(Flavour_Shower_Neutrino[MuonEnergyThreshMask])
Muon_Energy_L1 = Muon_Energy_L1[MuonEnergyThreshMask]
Muon_Energy_L2 = Muon_Energy_L2[MuonEnergyThreshMask]
Muon_UNCorr_Energy_L1 = Muon_UNCorr_Energy_L1[MuonEnergyThreshMask]
#print(len(Muon_Energy_L2))
GaisserH4a_weight = GaisserH4a_weight[MuonEnergyThreshMask]

# MuonEnergyThreshMask=[energy>10 for energy in Energy_L1_muon]
# Energy_Shower_Neutrino=(Energy_Shower_Neutrino[MuonEnergyThreshMask])
# Flavour_Shower_Neutrino=(Flavour_Shower_Neutrino[MuonEnergyThreshMask])
# Energy_L1_muon=(Energy_L1_muon[MuonEnergyThreshMask])
# # Energy_L2_muon=(Energy_L2_muon[MuonEnergyThreshMask])
# Weight_Gaisser_corsika=(Weight_Gaisser_corsika[MuonEnergyThreshMask])

# Energy_L1_muon=(Energy_L1_muon[MuonEnergyThreshMask])

#MultiplicityMasks=[np.logical_and(SingleMuonMask,Flav_mask)]
MultiplicityMasks=[np.logical_and(DoubleMuonMask,Flav_mask)]
print('len mult mask', len(MultiplicityMasks))
print('MultiplicityLen',np.sum(MultiplicityMasks))
print('angle mask len', np.sum(angle_mask))
print('depth mask len', np.sum(depth_mask))
for multmask in MultiplicityMasks:
    ADmask=np.logical_and.reduce((angle_mask,depth_mask,multmask))
    print('AD mask sum',np.sum(ADmask))
    #print('E1 AD',Muon_Energy_L1[ADmask])
    if np.sum(ADmask)>0:
        print('I AM IN HERE')
        Hist2D_trans,xedges,yedges = np.histogram2d(Energy_Shower_Neutrino[ADmask],Muon_Energy_L1[ADmask],bins = [E_nu_bins,E_mu_bins],weights = GaisserH4a_weight[ADmask])
        #Hist2D_trans,xedges,yedges = np.histogram2d(Muon_Energy_L1[ADmask],Muon_Energy_L2[ADmask],bins = [E_mu_bins,E_mu_bins],weights = GaisserH4a_weight[ADmask])
        #print(Hist2D_trans)
        Hist2D_norm_rows = (Hist2D_trans.T/np.sum(Hist2D_trans,axis=1)).T
        #print(Hist2D_norm_rows)
        Hist2D=Hist2D_norm_rows.T
        filename = '/home/zrechav/scripts/air_shower_reader/single_histogram_making/histograms/' + flavour + '_Single_Zen_' + str(coszen) + '_Depth_' + str(depth) + '.npy'
        #filename = '/home/zrechav/scripts/air_shower_reader/single_histogram_making/histograms/' + flavour + '_Double_Zen_' + str(coszen) + '_Depth_' + str(depth) + '.npy'
        print('Zen_'+str(coszen)+'_Depth_'+str(depth)+' '+str(np.sum(Hist2D)))
        np.save(filename, Hist2D)
    else:
        print('Empty array')
    