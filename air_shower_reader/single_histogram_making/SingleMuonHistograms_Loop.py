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
 



E_nu_bins=np.logspace(1,7,12+1)
nu_bin_centers = np.sqrt(E_nu_bins[:-1] * E_nu_bins[1:])
E_mu_bins=np.logspace(1,7,12+1)
mu_bin_centers = np.sqrt(E_mu_bins[:-1] * E_mu_bins[1:])




##Parameter in file name is lower boundary of bin range (ex: 0.04-0.08 cos(zenith))
SingleMuonMask=np.load('/home/zrechav/scripts/air_shower_reader/picklemaking/masked_pickles/SingleMuonMask.npy',mmap_mode='r',allow_pickle=True)


Muon_UNCorr_Energy_L1 = np.array([])
Energy_Shower_Neutrino = np.array([])   
Zenith_Shower_Neutrino = np.array([])   
Depth_Shower_Neutrino = np.array([])   
Flavour_Shower_Neutrino = np.array([])
GaisserH4a_weight = np.array([])

filename = '/data/user/zrechav/compiled_hdf5s/airshowered_corsika/Corsika_20904_5000.hdf5'
if (path.exists(filename)):
    try:
        print('I EXIST')
        hdf=h5py.File(filename, 'r')

        GaisserH4a_weight = np.append(GaisserH4a_weight,np.asarray(hdf.get('GaisserH4a_weight')['item']))
        Zenith_Shower_Neutrino = np.append(Zenith_Shower_Neutrino,np.asarray(hdf.get('shower_neutrino_zenith')['value']))
        Flavour_Shower_Neutrino = np.append(Flavour_Shower_Neutrino,np.asarray(hdf.get('shower_neutrino_type')['value']))
        Energy_Shower_Neutrino = np.append(Energy_Shower_Neutrino,np.asarray(hdf.get('shower_neutrino_energy')['value']))
        Depth_Shower_Neutrino = np.append(Depth_Shower_Neutrino,np.asarray(hdf.get('shower_neutrino_depth')['value']))
                    
        
        Muon_UNCorr_Energy_L1=np.append(Muon_UNCorr_Energy_L1,np.asarray(hdf.get('Muon_UnCorr_Energy_L1')['value']))
       
        


        hdf.close()
    except Exception as e:
        print(filename+" Is Faulty")
        print(e)

        hdf.close()
        
MuonEnergyThreshMask=[energy>10. for energy in Muon_UNCorr_Energy_L1]
Energy_Shower_Neutrino=np.float32(Energy_Shower_Neutrino[MuonEnergyThreshMask])
Flavour_Shower_Neutrino=np.float32(Flavour_Shower_Neutrino[MuonEnergyThreshMask])
Muon_UNCorr_Energy_L1 = Muon_UNCorr_Energy_L1[MuonEnergyThreshMask]
GaisserH4a_weight = GaisserH4a_weight[MuonEnergyThreshMask]



flavours = ['NuE','NuMu','NuTau']


angles= np.linspace(0,1.0,6)

upp = np.linspace(1.40,2.0,2)
low = np.linspace(2.1,2.5,2)
depths = np.concatenate([upp,low])



for flavour in flavours:
    Flav_mask = np.load('/home/zrechav/scripts/air_shower_reader/picklemaking/masked_pickles/Flavour_mask_' + str(flavour) +'.npy', mmap_mode = 'r', allow_pickle=True)
    for angle in angles:
        angle_mask = np.load('/home/zrechav/scripts/air_shower_reader/picklemaking/masked_pickles/Angle_mask_' + str(np.round(angle,2)) + '.npy', mmap_mode = 'r', allow_pickle=True)
        for depth in depths:
            depth_mask = np.load('/home/zrechav/scripts/air_shower_reader/picklemaking/masked_pickles/Depth_mask_' + str(np.round(depth,2)) + '.npy', mmap_mode='r',allow_pickle=True)



            print(flavour,angle,depth)

            MultiplicityMasks=[np.logical_and(SingleMuonMask,Flav_mask)]
            print('len mult mask', len(MultiplicityMasks))
            print('MultiplicityLen',np.sum(MultiplicityMasks))
            print('angle mask len', np.sum(angle_mask))
            print('depth mask len', np.sum(depth_mask))
            for multmask in MultiplicityMasks:
                ADmask=np.logical_and.reduce((angle_mask,depth_mask,multmask))
                print('AD mask sum',np.sum(ADmask))

                if np.sum(ADmask)>0:
                    print('I AM IN HERE')
                    Hist2D_trans,xedges,yedges = np.histogram2d(Energy_Shower_Neutrino[ADmask], Muon_UNCorr_Energy_L1[ADmask], bins = [E_nu_bins,E_mu_bins], weights = GaisserH4a_weight[ADmask])

                    #Hist2D_norm_rows = (Hist2D_trans.T/np.sum(Hist2D_trans,axis=1)).T ##why this step?

                    #Hist2D=Hist2D_norm_rows.T ##why this step?
                    Hist2D= Hist2D_trans/np.sum(Hist2D_trans)
                    filename = '/home/zrechav/scripts/air_shower_reader/single_histogram_making/histograms/' + flavour + '_Single_Zen_' + str(np.round(angle,2)) + '_Depth_' + str(np.round(depth,2)) + '.npy'
                    np.save(filename, Hist2D)
                else:
                    print('Empty array')
    