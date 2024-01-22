#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/vbasu/meta-projects/combo3/build
import numpy as np
import math
import os.path
from os import path
from icecube.tableio import I3TableWriter
from icecube.hdfwriter import I3HDFTableService, I3HDFWriter
import tables
import h5py
import numpy as np
import math
import os.path
from os import path
from icecube.tableio import I3TableWriter
from icecube.hdfwriter import I3HDFTableService, I3HDFWriter
import tables
import glob
import h5py
import pandas as pd
from icecube.dataclasses import I3Double, I3Particle, I3Direction, I3Position, I3VectorI3Particle, I3Constants, I3VectorOMKey

import argparse

# handling of command line arguments  
from optparse import OptionParser
parser = OptionParser()
usage = """%prog [options]"""
parser.set_usage(usage)


parser.add_option("-f", "--flavour", action="store", type="string", default='NuMu', dest="FLAVOUR", help="Flavour")
parser.add_option("-a", "--angle", action="store", type="float", default=0.17, dest="ANGLE", help="Input angle")
parser.add_option("-d", "--depth", action="store", type="float", default=1.5, dest="DEPTH", help="Inputdepth")


# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()

coszen=np.around(options.ANGLE,decimals=2)
depth=np.around(options.DEPTH,decimals=2)
flavour=options.FLAVOUR

QuadMuonMask=np.load('/data/user/vbasu/SelfVetoArrays/PickleFiles/QuadrupleMuonMask.npy',mmap_mode='r',allow_pickle=True)
angle_mask=np.load('/data/user/vbasu/SelfVetoArrays/PickleFiles/Angle_mask_'+str(coszen)+'.npy',mmap_mode='r',allow_pickle=True)
depth_mask=np.load('/data/user/vbasu/SelfVetoArrays/PickleFiles/Depth_mask_'+str(depth)+'.npy',mmap_mode='r',allow_pickle=True)
Flav_mask=np.load('/data/user/vbasu/SelfVetoArrays/PickleFiles/Flav_mask_'+str(flavour)+'.npy',mmap_mode='r',allow_pickle=True)
# print('Double Muons:',np.sum(DoubleMuonMask))

E_nu_bins=np.logspace(1,7,10+1)
nu_bin_centers = np.sqrt(E_nu_bins[:-1] * E_nu_bins[1:])
E_mu_bins=np.logspace(-1,7,10+1)
mu_bin_centers = np.sqrt(E_mu_bins[:-1] * E_mu_bins[1:])
# Energy_L1_muon=np.load('/data/user/vbasu/SelfVetoArrays/PickleFiles/Energy_L1_muon.npz',mmap_mode='r',allow_pickle=True)
# Energy_L2_muon=np.load('/data/user/vbasu/SelfVetoArrays/PickleFiles/Energy_L2_muon.npz',mmap_mode='r',allow_pickle=True)
# Energy_Shower_Neutrino=np.load('/data/user/vbasu/SelfVetoArrays/PickleFiles/Energy_Shower_Neutrino.npz',mmap_mode='r',allow_pickle=True)
# Weight_Gaisser_corsika=np.load('/data/user/vbasu/SelfVetoArrays/PickleFiles/Weight_Gaisser_corsika.npz',mmap_mode='r',allow_pickle=True)
Energy_brightest_muon,Energy_max_muon,max_depo_energy,Energy_Shower_Neutrino=np.array([]),np.array([]),np.array([]),np.array([])
Weight_Gaisser_corsika,Weight_Hoerandel_corsika,Weight_Honda_corsika=np.array([]),np.array([]),np.array([])
Energy_brightest_muon_depth_track,Energy_brightest_muon_at_depth,Energy_HE_muon_at_depth=np.array([]),np.array([]),np.array([])
MuonMultiplicity,brightest_muon_zenith,brightest_muon_depth=np.array([]),np.array([]),np.array([])
MuonMultiplicity_10,MuonMultiplicity_100,MuonMultiplicity_200=np.array([]),np.array([]),np.array([])
PrimaryMass,PrimaryEnergyperNucleon=np.array([]),np.array([])
Flavour_Shower_Neutrino=np.array([])
# TotalMuonEnergy=np.array([])
NEvents_corsika,Weight_corsika,PrimaryDepthMC_corsika,PrimaryZenithMC_corsika,PrimaryEnergyMC_corsika=np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
Energy_L1_muon,Depth_L1_muon,Zenith_L1_muon,NeutrinoSep_L1_muon,ShowerSep_L1_muon=np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
Energy_L2_muon,Depth_L2_muon,Zenith_L2_muon,NeutrinoSep_L2_muon,ShowerSep_L2_muon=np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
Energy_L3_muon,Depth_L3_muon,Zenith_L3_muon,NeutrinoSep_L3_muon,ShowerSep_L3_muon=np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
Energy_L4_muon,Depth_L4_muon,Zenith_L4_muon,NeutrinoSep_L4_muon,ShowerSep_L4_muon=np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
Energy_L5_muon,Depth_L5_muon,Zenith_L5_muon,NeutrinoSep_L5_muon,ShowerSep_L5_muon=np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
filename='/data/user/vbasu/CONDOR_output/Compiled_HDF/CorsikaShowers.hdf5'
print(filename)
# filename="/home/vbasu/testoutput/CorsikaShowerReader2.hdf5"
if (path.exists(filename)):
    try:
        hdf=h5py.File(filename, 'r')
        # Energy_brightest_muon=np.append(Energy_brightest_muon,np.asarray(hdf.get('Energy_brightest_muon')))
        Energy_Shower_Neutrino=np.append(Energy_Shower_Neutrino,np.asarray(hdf.get('Energy_Shower_Neutrino')))
        Flavour_Shower_Neutrino=np.append(Flavour_Shower_Neutrino,np.asarray(hdf.get('Flavour_Shower_Neutrino')))
        # Energy_brightest_muon_depth_track=np.append(Energy_brightest_muon_depth_track,np.asarray(hdf.get('Energy_brightest_muon_depth_track')))
        # Energy_brightest_muon_at_depth=np.append(Energy_brightest_muon_at_depth,np.asarray(hdf.get('Energy_brightest_muon_at_depth')))
        # Energy_HE_muon_at_depth=np.append(Energy_HE_muon_at_depth,np.asarray(hdf.get('Energy_HEmuon_at_depth')))
        # MuonMultiplicity=np.append(MuonMultiplicity,np.asarray(hdf.get('MuonMultiplicity')))
        # MuonMultiplicity_10=np.append(MuonMultiplicity_10,np.asarray(hdf.get('MuonMultiplicity_10')))
        # MuonMultiplicity_100=np.append(MuonMultiplicity_100,np.asarray(hdf.get('MuonMultiplicity_100')))
        # MuonMultiplicity_200=np.append(MuonMultiplicity_200,np.asarray(hdf.get('MuonMultiplicity_200')))
        # brightest_muon_zenith=np.append(brightest_muon_zenith,np.asarray(hdf.get('brightest_muon_zenith')))
        # brightest_muon_depth=np.append(brightest_muon_depth,np.asarray(hdf.get('brightest_muon_depth')))
        # Energy_max_muon=np.append(Energy_max_muon,np.asarray(hdf.get('Energy_max_muon')))
        # max_depo_energy=np.append(max_depo_energy,np.asarray(hdf.get('max_depo_energy')))
        Weight_Gaisser_corsika=np.append(Weight_Gaisser_corsika,np.asarray(hdf.get('Weight_Gaisser_corsika')))
        # Weight_Hoerandel_corsika=np.append(Weight_Hoerandel_corsika,np.asarray(hdf.get('Weight_Hoerandel_corsika')))
        # Weight_Honda_corsika=np.append(Weight_Honda_corsika,np.asarray(hdf.get('Weight_Honda_corsika')))
        print(filename)
        Energy_L1_muon=np.append(Energy_L1_muon,np.asarray(hdf.get('Energy_L1_muon')))
        Energy_L2_muon=np.append(Energy_L2_muon,np.asarray(hdf.get('Energy_L2_muon')))
        Energy_L3_muon=np.append(Energy_L3_muon,np.asarray(hdf.get('Energy_L3_muon')))
        Energy_L4_muon=np.append(Energy_L4_muon,np.asarray(hdf.get('Energy_L4_muon')))
        # Energy_L5_muon=np.append(Energy_L5_muon,np.asarray(hdf.get('Energy_L5_muon')))
        hdf.close()
    except Exception as e:
        print(filename+" Is Faulty")
        print(e)

        hdf.close()
# MuonEnergyThreshMask=[energy>10 for energy in Energy_L1_muon]
# Energy_Shower_Neutrino=(Energy_Shower_Neutrino[MuonEnergyThreshMask])
# Flavour_Shower_Neutrino=(Flavour_Shower_Neutrino[MuonEnergyThreshMask])
# Energy_L1_muon=(Energy_L1_muon[MuonEnergyThreshMask])
# Energy_L2_muon=(Energy_L2_muon[MuonEnergyThreshMask])
# Energy_L3_muon=(Energy_L3_muon[MuonEnergyThreshMask])
# Energy_L4_muon=(Energy_L4_muon[MuonEnergyThreshMask])

# Weight_Gaisser_corsika=(Weight_Gaisser_corsika[MuonEnergyThreshMask])

MultiplicityMasks=[np.logical_and(QuadMuonMask,Flav_mask)]
for multmask in MultiplicityMasks:
    ADmask=np.logical_and.reduce((angle_mask,depth_mask,multmask))
    print('AD mask sum',np.sum(ADmask))
    print('E1 AD',Energy_L1_muon[ADmask])
    if np.sum(ADmask)>0:
        stackarray=np.stack((Energy_L1_muon[ADmask],Energy_L2_muon[ADmask],Energy_L3_muon[ADmask],Energy_L4_muon[ADmask],Energy_Shower_Neutrino[ADmask]),axis=1)
        Hist5D,edges=np.histogramdd(stackarray,bins=(E_mu_bins,E_mu_bins,E_mu_bins,E_mu_bins,E_nu_bins),weights=Weight_Gaisser_corsika[ADmask])
        filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_Quadruple_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'
        print('Zen_'+str(coszen)+'_Depth_'+str(depth)+' '+str(np.sum(Hist5D)))
        np.save(filename, Hist5D)
    else:
        print('Empty Array')