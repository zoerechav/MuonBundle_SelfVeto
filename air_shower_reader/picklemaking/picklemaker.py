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
MuonMultiplicity_UNCorr= np.array([])



Max_Deposited_Energy = np.array([])     # (complete)
Energy_Shower_Neutrino = np.array([])   # (complete)
Zenith_Shower_Neutrino = np.array([])   # (complete)
Depth_Shower_Neutrino = np.array([])    # (complete)


Total_UNCorr_Muon_energy = np.array([])

Flavour_Shower_Neutrino = np.array([])
GaisserH4a_weight = np.array([])
PrimaryMass, PrimaryEnergyperNucleon= np.array([]),np.array([])
#filename = '/home/zrechav/test/AirShowerReader.hdf5'
filename = '/data/user/zrechav/compiled_hdf5s/airshowered_corsika/Corsika_20904_5000.hdf5'
#filename='/data/user/zrechav/airshowered_corsika/20904/0000000-0000999/airshowered_corsika_20904_000000.hdf5'
print(filename)
# filename="/home/vbasu/testoutput/CorsikaShowerReader2.hdf5"
if (path.exists(filename)):
    try:
        print('I EXIST')
        hdf=h5py.File(filename, 'r')
########ENERGY, ZENITH, AND WEIGHTING INFORMATION#######
        GaisserH4a_weight = np.append(GaisserH4a_weight,np.asarray(hdf.get('GaisserH4a_weight')['item']))
        Zenith_Shower_Neutrino = np.append(Zenith_Shower_Neutrino,np.asarray((hdf.get('shower_neutrino_zenith')['value'])))
            
        Flavour_Shower_Neutrino = np.append(Flavour_Shower_Neutrino,np.asarray(hdf.get('shower_neutrino_type')['value']))
        Energy_Shower_Neutrino = np.append(Energy_Shower_Neutrino,np.asarray(hdf.get('shower_neutrino_energy')['value']))
        Depth_Shower_Neutrino = np.append(Depth_Shower_Neutrino,np.asarray(hdf.get('shower_neutrino_depth')['value']))
        
        

########UNCORRELATED MUON INFORMATION############
        MuonMultiplicity_UNCorr = np.append(MuonMultiplicity_UNCorr,np.asarray(hdf.get('MuonMultiplicity_UnCorr')['value']))

        Muon_UNCorr_Energy_L1 = np.append(Muon_UNCorr_Energy_L1, np.asarray(hdf.get('Muon_UnCorr_Energy_L1')['value']))
        Muon_UNCorr_Energy_L2 = np.append(Muon_UNCorr_Energy_L2, np.asarray(hdf.get('Muon_UnCorr_Energy_L2')['value']))
        Muon_UNCorr_Energy_L3 = np.append(Muon_UNCorr_Energy_L3, np.asarray(hdf.get('Muon_UnCorr_Energy_L3')['value']))
        Muon_UNCorr_Energy_L4 = np.append(Muon_UNCorr_Energy_L4, np.asarray(hdf.get('Muon_UnCorr_Energy_L4')['value']))
        Muon_UNCorr_Energy_L5 = np.append(Muon_UNCorr_Energy_L5, np.asarray(hdf.get('Muon_UnCorr_Energy_L5')['value']))
        print('making progress')
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
        print('UAU')
        Muon_UNCorr_Neutrino_separation_L1 = np.append(Muon_UNCorr_Neutrino_separation_L1, np.asarray(hdf.get('Muon_UnCorr_L1_Neutrino_separation')['value']))
        Muon_UNCorr_Neutrino_separation_L2 = np.append(Muon_UNCorr_Neutrino_separation_L2, np.asarray(hdf.get('Muon_UnCorr_L2_Neutrino_separation')['value']))
        Muon_UNCorr_Neutrino_separation_L3 = np.append(Muon_UNCorr_Neutrino_separation_L3, np.asarray(hdf.get('Muon_UnCorr_L3_Neutrino_separation')['value']))
        Muon_UNCorr_Neutrino_separation_L4 = np.append(Muon_UNCorr_Neutrino_separation_L4, np.asarray(hdf.get('Muon_UnCorr_L4_Neutrino_separation')['value']))
        Muon_UNCorr_Neutrino_separation_L5 = np.append(Muon_UNCorr_Neutrino_separation_L5, np.asarray(hdf.get('Muon_UnCorr_L5_Neutrino_separation')['value']))
        print('almost therrrrr')
        Muon_UNCorr_Shower_separation_L1 = np.append(Muon_UNCorr_Shower_separation_L1, np.asarray(hdf.get('Muon_UnCorr_L1_Shower_separation')['value']))
        Muon_UNCorr_Shower_separation_L2 = np.append(Muon_UNCorr_Shower_separation_L2, np.asarray(hdf.get('Muon_UnCorr_L2_Shower_separation')['value']))
        Muon_UNCorr_Shower_separation_L3 = np.append(Muon_UNCorr_Shower_separation_L3, np.asarray(hdf.get('Muon_UnCorr_L3_Shower_separation')['value']))
        Muon_UNCorr_Shower_separation_L4 = np.append(Muon_UNCorr_Shower_separation_L4, np.asarray(hdf.get('Muon_UnCorr_L4_Shower_separation')['value']))
        Muon_UNCorr_Shower_separation_L5 = np.append(Muon_UNCorr_Shower_separation_L5, np.asarray(hdf.get('Muon_UnCorr_L5_Shower_separation')['value']))

        Total_UNCorr_Muon_energy = np.append(Total_UNCorr_Muon_energy,np.asarray(hdf.get('Total_UnCorr_Muon_energy')['value']))
        print('I FINISHED LOADING MY ARRAYS WITH HDF INFO')
        hdf.close()
    except Exception as e:
        print(filename+" Is Faulty")
        print(e)

        hdf.close()

##MASKING AND SAVING INFORMATION TO SAVE TIME LATER ON                
            
MuonEnergyThreshMask=[energy>1e1 for energy in Muon_UNCorr_Energy_L1]


Energy_Shower_Neutrino=np.float32(Energy_Shower_Neutrino)
np.savez_compressed('/home/zrechav/scripts/air_shower_reader/picklemaking/pickles/Energy_Shower_Neutrino',Energy_Shower_Neutrino)
Flavour_Shower_Neutrino=np.float32(Flavour_Shower_Neutrino[MuonEnergyThreshMask])
np.savez_compressed('/home/zrechav/scripts/air_shower_reader/picklemaking/pickles/Flavour_Shower_Neutrino',Flavour_Shower_Neutrino)

# print('Flavour masks')
Flav_mask_NuE=np.where(np.abs(Flavour_Shower_Neutrino)==12., True, False)
np.save('/home/zrechav/scripts/air_shower_reader/picklemaking/masked_pickles/Flavour_mask_NuE.npy',Flav_mask_NuE)
Flav_mask_NuMu=np.where(np.abs(Flavour_Shower_Neutrino)==14., True, False)
np.save('/home/zrechav/scripts/air_shower_reader/picklemaking/masked_pickles/Flavour_mask_NuMu.npy',Flav_mask_NuMu)
Flav_mask_NuTau=np.where(np.abs(Flavour_Shower_Neutrino)==16., True, False)
np.save('/home/zrechav/scripts/air_shower_reader/picklemaking/masked_pickles/Flavour_mask_NuTau.npy',Flav_mask_NuTau)

##my plots do not express correlated muons by energy when I use these values -- try expressly using uncorrelated muons
Muon_Energy_L2 = Muon_UNCorr_Energy_L2[MuonEnergyThreshMask]
Muon_Energy_L3 = Muon_UNCorr_Energy_L3[MuonEnergyThreshMask]
Muon_Energy_L4 = Muon_UNCorr_Energy_L4[MuonEnergyThreshMask]
Muon_Energy_L5 = Muon_UNCorr_Energy_L5[MuonEnergyThreshMask]

L2_mask=np.where(Muon_Energy_L2==0, True, False).tolist()
L3_mask=np.where(Muon_Energy_L3==0, True, False).tolist()
L4_mask=np.where(Muon_Energy_L4==0, True, False).tolist()
L5_mask=np.where(Muon_Energy_L5==0, True, False).tolist()

SingleMuonMask=L2_mask
# print('Single Muons:',np.sum(SingleMuonMask))
np.save('/home/zrechav/scripts/air_shower_reader/picklemaking/masked_pickles/SingleMuonMask.npy',SingleMuonMask)
DoubleMuonMask=np.logical_and(L3_mask,np.logical_not(L2_mask))
# print('Double Muons:',np.sum(DoubleMuonMask))
np.save('/home/zrechav/scripts/air_shower_reader/picklemaking/masked_pickles/DoubleMuonMask.npy',DoubleMuonMask)
TripleMuonMask=np.logical_and.reduce((L4_mask,np.logical_not(L3_mask),np.logical_not(L2_mask)))
# print('Triple Muons:',np.sum(TripleMuonMask))
np.save('/home/zrechav/scripts/air_shower_reader/picklemaking/masked_pickles/TripleMuonMask.npy',TripleMuonMask)
QuadMuonMask=np.logical_and.reduce((L5_mask,np.logical_not(L4_mask),np.logical_not(L3_mask),np.logical_not(L2_mask)))
# print('Quadruple Muons:',np.sum(QuadMuonMask))
np.save('/home/zrechav/scripts/air_shower_reader/picklemaking/masked_pickles/QuadrupleMuonMask.npy',QuadMuonMask)

##DEFINE ENERGY, ZENITH, AND DEPTH BINS HERE 
E_nu_bins=np.logspace(1,7,12+1)
print(E_nu_bins)
nu_bin_centers = np.sqrt(E_nu_bins[:-1] * E_nu_bins[1:])
E_mu_bins=np.logspace(1,7,12+1)
print(E_mu_bins)
mu_bin_centers = np.sqrt(E_mu_bins[:-1] * E_mu_bins[1:])

angles_space= np.linspace(0.,1.0,6)
#high_angles=np.linspace(0.5,1,11)
# new_low_angles=np.delete(low_angles,0) 
#angles_space=np.concatenate( (low_angles,high_angles) )

print(angles_space)
#mult_bins=np.logspace(0,5,20+1)

#########REPLACE WITH PRIMARY DEPTHMC ZENITH AND DEPTH ASAP
Zenith_Shower_Neutrino = Zenith_Shower_Neutrino[MuonEnergyThreshMask]
angle_vals_prim=np.digitize(np.cos(Zenith_Shower_Neutrino),angles_space)-1

#Depth_BrightestMuon = Depth_BrightestMuon[MuonEnergyThreshMask]
Depth_Shower_Neutrino = Depth_Shower_Neutrino[MuonEnergyThreshMask]
dust_layer = np.linspace(2.,2.1,1)
upp = np.linspace(1.40,2.0,2)
low = np.linspace(2.1,2.5,2)
eh = [3.0]
depth_space= np.concatenate([upp,low])
print(depth_space)
#depth_space = np.linspace(1.42,2.5,28)
depths_prim=np.digitize(Depth_Shower_Neutrino/1000,depth_space)

print('angle masks')
angle_masks_prim = [angle_vals_prim == i for i in range(len(angles_space))]
for i,angle in enumerate(angles_space):
    print(np.sum(angle_masks_prim[i]),',',angle) 
    np.save('/home/zrechav/scripts/air_shower_reader/picklemaking/masked_pickles/Angle_mask_' + str(np.around(angle,decimals=2))+'.npy',angle_masks_prim[i])

print('depth masks')
depth_masks_prim = [depths_prim == i for i in range(len(depth_space))]
for i,depth in enumerate(depth_space):
    print(np.sum(depth_masks_prim[i]),',',depth)
    np.save('/home/zrechav/scripts/air_shower_reader/picklemaking/masked_pickles/Depth_mask_' +str(np.around(depth,decimals=2))+'.npy',depth_masks_prim[i])
