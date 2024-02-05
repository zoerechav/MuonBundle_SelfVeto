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
parser.add_argument("-a", "--angle", action="store", type=float, default=0.2, dest="ANGLE", help="Input angle")
parser.add_argument("-d", "--depth", action="store", type=float, default=2.1, dest="DEPTH", help="Inputdepth")


# parse cmd line args, bail out if anything is not understood
args = parser.parse_args()

coszen=np.around(args.ANGLE,decimals=2)
depth=np.around(args.DEPTH,decimals=2)
flavour=args.FLAVOUR
##Parameter in file name is lower boundary of bin range (ex: 0.04-0.08 cos(zenith))
DoubleMuonMask = np.load('/home/zrechav/scripts/air_shower_reader/picklemaking/masked_pickles/DoubleMuonMask.npy', mmap_mode='r',allow_pickle=True)

E_nu_bins=np.logspace(1,7,12+1)
nu_bin_centers = np.sqrt(E_nu_bins[:-1] * E_nu_bins[1:])
E_mu_bins=np.logspace(1,7,12+1)
mu_bin_centers = np.sqrt(E_mu_bins[:-1] * E_mu_bins[1:])
print(E_nu_bins)

Muon_UNCorr_Energy_L1 = np.array([])
Muon_UNCorr_Energy_L2 = np.array([])
Energy_Shower_Neutrino = np.array([])   
Zenith_Shower_Neutrino = np.array([])   
Depth_Shower_Neutrino = np.array([])   
Flavour_Shower_Neutrino = np.array([])
GaisserH4a_weight = np.array([])
Total_Muon_energy = np.array([])

filename = '/data/user/zrechav/compiled_hdf5s/airshowered_corsika/Corsika_20904_5000.hdf5'
print(filename)

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

        Muon_UNCorr_Energy_L1=np.append(Muon_UNCorr_Energy_L1,np.asarray(hdf.get('Muon_UnCorr_Energy_L1')['value']))
        Muon_UNCorr_Energy_L2=np.append(Muon_UNCorr_Energy_L2,np.asarray(hdf.get('Muon_UnCorr_Energy_L2')['value']))
        print(Muon_UNCorr_Energy_L2)
        print('DID I MAKE IT HERE?')
        

        print('DOUBLE CHECKING THIS IS INDEED THE FIX')
        
        
        hdf.close()
    except Exception as e:
        print(filename+" Is Faulty")
        print(e)

        hdf.close()


MuonEnergyThreshMask=[energy>1e1 for energy in Muon_UNCorr_Energy_L1]
Energy_Shower_Neutrino=np.float32(Energy_Shower_Neutrino[MuonEnergyThreshMask])
Flavour_Shower_Neutrino=np.float32(Flavour_Shower_Neutrino[MuonEnergyThreshMask])
Muon_Energy_L1 = Muon_UNCorr_Energy_L1[MuonEnergyThreshMask]
Muon_Energy_L2 = Muon_UNCorr_Energy_L2[MuonEnergyThreshMask]
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
            MultiplicityMasks=[np.logical_and(DoubleMuonMask,Flav_mask)]
            #print('Thresh mask',np.sum(MuonEnergyThreshMask))

            #print('Angle Mask',np.shape(angle_mask))
            #print('Depth Mask',np.shape(depth_mask))
            #print('DoubleMuonMask',np.shape(DoubleMuonMask))
            #print('Flavour mask',np.shape(Flav_mask))
            #print('Energy:',np.shape(Muon_Energy_L1))
            print('angle mask len', np.sum(angle_mask))
            print('depth mask len', np.sum(depth_mask))
            for multmask in MultiplicityMasks:
            
            
                ADmask=np.logical_and.reduce((angle_mask,depth_mask,multmask))
                print('AD mask sum',np.sum(ADmask))
                #print('E1 AD',Muon_Energy_L1[ADmask])
                #print('E2 AD',Muon_Energy_L2[ADmask])
                if np.sum(ADmask)>0:

                    stackarray = np.stack((Muon_Energy_L1[ADmask], Muon_Energy_L2[ADmask], Energy_Shower_Neutrino[ADmask]), axis=1)
                    #print(stackarray.shape)
                    Hist3D,edges=np.histogramdd(stackarray,bins=(E_mu_bins,E_mu_bins,E_nu_bins),weights=GaisserH4a_weight[ADmask])

                    #print(Hist3D.shape)
                    filename = '/home/zrechav/scripts/air_shower_reader/double_histogram_making/histograms/' + flavour + '_Double_Zen_' + str(np.around(angle,2)) + '_Depth_' + str(np.around(depth,2)) + '.npy'   #filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_Double_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'
                    print('Zen_'+str(angle)+'_Depth_'+str(depth)+' '+str(np.sum(Hist3D)))
                    np.save(filename, Hist3D)
                else:
                    print('Empty array')