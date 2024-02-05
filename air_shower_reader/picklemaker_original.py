#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/vbasu/meta-projects/combo2/build

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

livetime_year=365.*24*60*60
flav_dict=    {'NuMu': '21217', 'NuE': '21218', 'NuTau': '21219'}
type_dict={'NuMu': I3Particle.NuMu, 'NuE': I3Particle.NuE, 'NuTau': I3Particle.NuTau}

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
        # Energy_L1_muon=np.append(Energy_L1_muon,np.asarray(hdf.get('Energy_L1_muon')))
        # Energy_L2_muon=np.append(Energy_L2_muon,np.asarray(hdf.get('Energy_L2_muon')))
        # Energy_L3_muon=np.append(Energy_L3_muon,np.asarray(hdf.get('Energy_L3_muon')))
        # Energy_L4_muon=np.append(Energy_L4_muon,np.asarray(hdf.get('Energy_L4_muon')))
        # Energy_L5_muon=np.append(Energy_L5_muon,np.asarray(hdf.get('Energy_L5_muon')))
        # Depth_L1_muon=np.append(Depth_L1_muon,np.asarray(hdf.get('Depth_L1_muon')))
        # Depth_L2_muon=np.append(Depth_L2_muon,np.asarray(hdf.get('Depth_L2_muon')))
        # Depth_L3_muon=np.append(Depth_L3_muon,np.asarray(hdf.get('Depth_L3_muon')))
        # Depth_L4_muon=np.append(Depth_L4_muon,np.asarray(hdf.get('Depth_L4_muon')))
        # Depth_L5_muon=np.append(Depth_L5_muon,np.asarray(hdf.get('Depth_L5_muon')))
        # Zenith_L1_muon=np.append(Zenith_L1_muon,np.asarray(hdf.get('Zenith_L1_muon')))
        # Zenith_L2_muon=np.append(Zenith_L2_muon,np.asarray(hdf.get('Zenith_L2_muon')))
        # Zenith_L3_muon=np.append(Zenith_L3_muon,np.asarray(hdf.get('Zenith_L3_muon')))
        # Zenith_L4_muon=np.append(Zenith_L4_muon,np.asarray(hdf.get('Zenith_L4_muon')))
        # Zenith_L5_muon=np.append(Zenith_L5_muon,np.asarray(hdf.get('Zenith_L5_muon')))
#         NeutrinoSep_L1_muon=np.append(NeutrinoSep_L1_muon,np.asarray(hdf.get('NeutrinoSep_L1_muon')))
#         NeutrinoSep_L2_muon=np.append(NeutrinoSep_L2_muon,np.asarray(hdf.get('NeutrinoSep_L2_muon')))
#         NeutrinoSep_L3_muon=np.append(NeutrinoSep_L3_muon,np.asarray(hdf.get('NeutrinoSep_L3_muon')))
#         NeutrinoSep_L4_muon=np.append(NeutrinoSep_L4_muon,np.asarray(hdf.get('NeutrinoSep_L4_muon')))
#         NeutrinoSep_L5_muon=np.append(NeutrinoSep_L5_muon,np.asarray(hdf.get('NeutrinoSep_L5_muon')))
#         ShowerSep_L1_muon=np.append(ShowerSep_L1_muon,np.asarray(hdf.get('ShowerSep_L1_muon')))
#         ShowerSep_L2_muon=np.append(ShowerSep_L2_muon,np.asarray(hdf.get('ShowerSep_L2_muon')))
#         ShowerSep_L3_muon=np.append(ShowerSep_L3_muon,np.asarray(hdf.get('ShowerSep_L3_muon')))
#         ShowerSep_L4_muon=np.append(ShowerSep_L4_muon,np.asarray(hdf.get('ShowerSep_L4_muon')))
#         ShowerSep_L5_muon=np.append(ShowerSep_L5_muon,np.asarray(hdf.get('ShowerSep_L5_muon')))
        
#         TotalMuonEnergy=np.append(TotalMuonEnergy,np.asarray(hdf.get('TotalMuonEnergy')))
        print(filename)
        PrimaryDepthMC_corsika=np.append(PrimaryDepthMC_corsika,np.asarray(hdf.get('PrimaryDepthMC')))
        # PrimaryEnergyMC_corsika=np.append(PrimaryEnergyMC_corsika,np.asarray(hdf.get('PrimaryEnergyMC')))
        PrimaryZenithMC_corsika=np.append(PrimaryZenithMC_corsika,np.asarray(hdf.get('PrimaryZenithMC')))
        # PrimaryMass=np.append(PrimaryMass,np.asarray(hdf.get('PrimaryMass')))
        

        hdf.close()
    except Exception as e:
        print(filename+" Is Faulty")
        print(e)

        hdf.close()

                
PrimaryEnergyperNucleon=PrimaryEnergyMC_corsika/PrimaryMass  
# print(len(Energy_brightest_muon)) 
print(len(PrimaryZenithMC_corsika))
print(np.amin(Energy_Shower_Neutrino),np.amax(Energy_Shower_Neutrino))
# print(np.amin(Energy_brightest_muon_at_depth),np.amax(Energy_brightest_muon_at_depth))
# print(np.amin(Energy_HE_muon_at_depth),np.amax(Energy_HE_muon_at_depth))
# print(np.amin(PrimaryEnergyMC_corsika),np.amax(PrimaryEnergyMC_corsika))
# print(np.amin(max_depo_energy),np.amax(max_depo_energy))
# print(np.amin(MuonMultiplicity),np.amax(MuonMultiplicity))
# MuonEnergyThreshMask=[energy>10 for energy in Energy_L1_muon]
# Energy_brightest_muon=np.float32(Energy_brightest_muon[MuonEnergyThreshMask])
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/Energy_brightest_muon',Energy_brightest_muon)
# Energy_Shower_Neutrino=np.float32(Energy_Shower_Neutrino)
# np.savez_compressed('/data/user/vbasu/SelfVetoArrays/PickleFiles/Energy_Shower_Neutrino',Energy_Shower_Neutrino)
# Flavour_Shower_Neutrino=np.float32(Flavour_Shower_Neutrino[MuonEnergyThreshMask])
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/Flavour_Shower_Neutrino',Flavour_Shower_Neutrino)
# print('Flavour masks')
# Flav_mask_NuE=np.where(np.abs(Flavour_Shower_Neutrino)==12., True, False)
# np.save('/data/user/vbasu/SelfVetoArrays/PickleFiles/Flav_mask_NuE.npy',Flav_mask_NuE)
# Flav_mask_NuMu=np.where(np.abs(Flavour_Shower_Neutrino)==14., True, False)
# np.save('/data/user/vbasu/SelfVetoArrays/PickleFiles/Flav_mask_NuMu.npy',Flav_mask_NuMu)
# Flav_mask_NuTau=np.where(np.abs(Flavour_Shower_Neutrino)==16., True, False)
# np.save('/data/user/vbasu/SelfVetoArrays/PickleFiles/Flav_mask_NuTau.npy',Flav_mask_NuTau)
# Energy_brightest_muon_at_depth=np.float32(Energy_brightest_muon_at_depth[MuonEnergyThreshMask])
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/Energy_brightest_muon_at_depth',Energy_brightest_muon_at_depth)
# Energy_HE_muon_at_depth=np.float32(Energy_HE_muon_at_depth[MuonEnergyThreshMask])
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/Energy_HE_muon_at_depth',Energy_HE_muon_at_depth)

# TotalMuonEnergy=TotalMuonEnergy[MuonEnergyThreshMask]
# MuonMultiplicity=np.float32(MuonMultiplicity[MuonEnergyThreshMask])
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/MuonMultiplicity',MuonMultiplicity)
# MuonMultiplicity_10=np.float32(MuonMultiplicity_10[MuonEnergyThreshMask])
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/MuonMultiplicity_10',MuonMultiplicity_10)
# MuonMultiplicity_100=np.float32(MuonMultiplicity_100[MuonEnergyThreshMask])
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/MuonMultiplicity_100',MuonMultiplicity_100)
# MuonMultiplicity_200=np.float32(MuonMultiplicity_200[MuonEnergyThreshMask])
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/MuonMultiplicity_200',MuonMultiplicity_200)
# print('muon energies')
# # Energy_L1_muon=np.float32(Energy_L1_muon[MuonEnergyThreshMask])
# # np.savez_compressed('/data/user/vbasu/SelfVetoArrays/PickleFiles/Energy_L1_muon',Energy_L1_muon)
# # Energy_L2_muon=np.float32(Energy_L2_muon[MuonEnergyThreshMask])
# np.savez_compressed('/data/user/vbasu/SelfVetoArrays/PickleFiles/Energy_L2_muon',Energy_L2_muon)
# # Energy_L3_muon=np.float32(Energy_L3_muon[MuonEnergyThreshMask])
# np.savez_compressed('/data/user/vbasu/SelfVetoArrays/PickleFiles/Energy_L3_muon',Energy_L3_muon)
# # Energy_L4_muon=np.float32(Energy_L4_muon[MuonEnergyThreshMask])
# np.savez_compressed('/data/user/vbasu/SelfVetoArrays/PickleFiles/Energy_L4_muon',Energy_L4_muon)
# # Energy_L5_muon=np.float32(Energy_L5_muon[MuonEnergyThreshMask])
# np.savez_compressed('/data/user/vbasu/SelfVetoArrays/PickleFiles/Energy_L5_muon',Energy_L5_muon)

# Depth_L1_muon=np.float32(Depth_L1_muon[MuonEnergyThreshMask])
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/Depth_L1_muon',Depth_L1_muon)
# Depth_L2_muon=np.float32(Depth_L2_muon[MuonEnergyThreshMask])
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/Depth_L2_muon',Depth_L2_muon)
# Depth_L3_muon=np.float32(Depth_L3_muon[MuonEnergyThreshMask])
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/Depth_L3_muon',Depth_L3_muon)
# Depth_L4_muon=np.float32(Depth_L4_muon[MuonEnergyThreshMask])
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/Depth_L4_muon',Depth_L4_muon)
# Depth_L5_muon=np.float32(Depth_L5_muon[MuonEnergyThreshMask])
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/Depth_L5_muon',Depth_L5_muon)

# Zenith_L1_muon=np.float32(Zenith_L1_muon[MuonEnergyThreshMask])
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/Zenith_L1_muon',Zenith_L1_muon)
# Zenith_L2_muon=np.float32(Zenith_L2_muon[MuonEnergyThreshMask])
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/Zenith_L2_muon',Zenith_L2_muon)
# Zenith_L3_muon=np.float32(Zenith_L3_muon[MuonEnergyThreshMask])
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/Zenith_L3_muon',Zenith_L3_muon)
# Zenith_L4_muon=np.float32(Zenith_L4_muon[MuonEnergyThreshMask])
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/Zenith_L4_muon',Zenith_L4_muon)
# Zenith_L5_muon=np.float32(Zenith_L5_muon[MuonEnergyThreshMask])
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/Zenith_L5_muon',Zenith_L5_muon)

# brightest_muon_zenith=brightest_muon_zenith[MuonEnergyThreshMask]
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/brightest_muon_zenith',brightest_muon_zenith)
# brightest_muon_depth=brightest_muon_depth[MuonEnergyThreshMask]
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/brightest_muon_depth',brightest_muon_depth)
# Energy_max_muon=Energy_max_muon[MuonEnergyThreshMask]
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/Energy_max_muon',Energy_max_muon)
# max_depo_energy=np.float32(max_depo_energy[MuonEnergyThreshMask])
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/max_depo_energy',max_depo_energy)
# Weight_Gaisser_corsika=np.float32(Weight_Gaisser_corsika[MuonEnergyThreshMask])
#np.savez_compressed('/data/user/vbasu/SelfVetoArrays/PickleFiles/Weight_Gaisser_corsika',Weight_Gaisser_corsika)
# Weight_Hoerandel_corsika=Weight_Hoerandel_corsika[MuonEnergyThreshMask]
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/Weight_Hoerandel_corsika',Weight_Hoerandel_corsika)
# Weight_Honda_corsika=Weight_Honda_corsika[MuonEnergyThreshMask]
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/Weight_Honda_corsika',Weight_Honda_corsika)

# PrimaryEnergyMC_corsika=np.float32(PrimaryEnergyMC_corsika[MuonEnergyThreshMask])
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/PrimaryEnergyMC_corsika',PrimaryEnergyMC_corsika)

# PrimaryZenithMC_corsika=np.float32(PrimaryZenithMC_corsika[MuonEnergyThreshMask])
# print(len(PrimaryZenithMC_corsika))

# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/PrimaryZenithMC_corsika',PrimaryZenithMC_corsika)
# PrimaryDepthMC_corsika=np.float32(PrimaryDepthMC_corsika[MuonEnergyThreshMask])
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/PrimaryDepthMC_corsika',PrimaryDepthMC_corsika)
# PrimaryEnergyperNucleon=np.float32(PrimaryEnergyperNucleon[MuonEnergyThreshMask])
# np.savez_compressed('/home/vbasu/scripts/MESE_Studies/SelfVeto/PickleFiles/PrimaryEnergyperNucleon',PrimaryEnergyperNucleon)
# print('multiplicity masks')
# L2_mask=np.where(Energy_L2_muon==0, True, False).tolist()
# L3_mask=np.where(Energy_L3_muon==0, True, False).tolist()
# L4_mask=np.where(Energy_L4_muon==0, True, False).tolist()
# L5_mask=np.where(Energy_L5_muon==0, True, False).tolist()
# SingleMuonMask=L2_mask
# print('Single Muons:',np.sum(SingleMuonMask))
# np.save('/data/user/vbasu/SelfVetoArrays/PickleFiles/SingleMuonMask.npy',SingleMuonMask)
# DoubleMuonMask=np.logical_and(L3_mask,np.logical_not(L2_mask))
# print('Double Muons:',np.sum(DoubleMuonMask))
# np.save('/data/user/vbasu/SelfVetoArrays/PickleFiles/DoubleMuonMask.npy',DoubleMuonMask)
# TripleMuonMask=np.logical_and.reduce((L4_mask,np.logical_not(L3_mask),np.logical_not(L2_mask)))
# print('Triple Muons:',np.sum(TripleMuonMask))
# np.save('/data/user/vbasu/SelfVetoArrays/PickleFiles/TripleMuonMask.npy',TripleMuonMask)
# QuadMuonMask=np.logical_and.reduce((L5_mask,np.logical_not(L4_mask),np.logical_not(L3_mask),np.logical_not(L2_mask)))
# print('Quadruple Muons:',np.sum(QuadMuonMask))
# np.save('/data/user/vbasu/SelfVetoArrays/PickleFiles/QuadrupleMuonMask.npy',QuadMuonMask)


E_nu_bins=np.logspace(2,7,10+1)
nu_bin_centers = np.sqrt(E_nu_bins[:-1] * E_nu_bins[1:])
E_mu_bins=np.logspace(-1,7,10+1)
mu_bin_centers = np.sqrt(E_mu_bins[:-1] * E_mu_bins[1:])

low_angles= np.linspace(0,0.5,13)
high_angles=np.linspace(0.5,1,13)
# new_low_angles=np.delete(low_angles,0) 
angles_space=np.concatenate( (low_angles,high_angles) )


mult_bins=np.logspace(0,5,20+1)



angle_vals_prim=np.digitize(np.cos(PrimaryZenithMC_corsika),angles_space)-1

depth_space = np.linspace(1.42,2.5,28)
depths_prim=np.digitize(PrimaryDepthMC_corsika/1000,depth_space)

print('angle masks')
angle_masks_prim = [angle_vals_prim == i for i in range(len(angles_space))]
for i,angle in enumerate(angles_space):
    print(np.sum(angle_masks_prim[i]))
    np.save('/data/user/vbasu/SelfVetoArrays/PickleFiles/Angle_mask_'+str(np.around(angle,decimals=2))+'.npy',angle_masks_prim[i])

print('depth masks')
depth_masks_prim = [depths_prim == i for i in range(len(depth_space))]
for i,depth in enumerate(depth_space):
    np.save('/data/user/vbasu/SelfVetoArrays/PickleFiles/Depth_mask_'+str(np.around(depth,decimals=2))+'.npy',depth_masks_prim[i])
