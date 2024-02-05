#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/vbasu/meta-projects/combo2/build

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib import colors as clrs
import matplotlib as mpl
from scipy.ndimage import zoom
import os
import argparse  
from optparse import OptionParser

parser = OptionParser()
usage = """%prog [options]"""
parser.set_usage(usage)

parser.add_option("-f", "--flavour", action="store", type="string", default='NuMu', dest="FLAVOUR", help="Flavour")
parser.add_option("-a", "--angle", action="store", type="float", default=0.2, dest="ANGLE", help="Input angle")
parser.add_option("-d", "--depth", action="store", type="float", default=2.1, dest="DEPTH", help="Inputdepth")

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()

E_mu_bins = np.logspace(1,7,12+1)
E_nu_bins = np.logspace(1,7,12+1)
mu_bin_centers = np.sqrt(E_mu_bins[:-1] * E_mu_bins[1:])
nu_bin_centers = np.sqrt(E_nu_bins[:-1] * E_nu_bins[1:])

gridpts_mu = np.logspace(1,7,31+1)
gridpts_nu = np.logspace(1,7,31+1)

coszen=np.around(options.ANGLE,decimals=2)
depth=np.around(options.DEPTH,decimals=2)
flavour=options.FLAVOUR




def Make_3D_Spline(  filename,
                     E_mu_bins = np.logspace(1,7,12+1),
                     E_nu_bins = np.logspace(1,7,12+1),
                     gridpts_mu = np.logspace(1,7,31+1),
                     gridpts_nu = np.logspace(1,7,31+1),
                     spline_method = 'linear',
                     
                    ) :    
    
    mu_bin_centers = np.sqrt(E_mu_bins[:-1] * E_mu_bins[1:])
    nu_bin_centers = np.sqrt(E_nu_bins[:-1] * E_nu_bins[1:])
    
    Hist4D=np.load(filename)
    #print(len(E_nu_bins))
    #print(Hist4D.shape)
    
    Hist4D_splined=np.empty((len(gridpts_mu),len(gridpts_mu),len(gridpts_mu),len(E_nu_bins)))
    for test_index in range(len(E_nu_bins)-1):
        test_nu_energy=E_nu_bins[test_index]
        Hist3D=Hist4D[:,:,:,test_index]
        Hist3D=Hist3D/np.sum(Hist3D)
        Hist3D=np.nan_to_num(Hist3D)

        xv, yv,zv=np.meshgrid((mu_bin_centers),(mu_bin_centers),(mu_bin_centers))
        plotx,ploty,plotz=np.meshgrid((gridpts_mu),(gridpts_mu),(gridpts_mu))
        mask=np.where(Hist3D!=0)

        try:
            grid_z1 = griddata((xv[mask],yv[mask],zv[mask]), Hist3D[mask], (plotx,ploty,plotz), method=spline_method)      
            grid_z1=grid_z1/np.sum(grid_z1)
            grid_z1=np.nan_to_num(grid_z1)
            Hist4D_splined[:,:,:,test_index]=grid_z1.T
            
        except Exception as e:
            #print(e)
            finearray=(Hist3D.T).repeat(2, 0).repeat(2, 1).repeat(2,2)
            length = finearray.shape[0]
            # Resize the array to shape (32, 32) using zoom, so Hist3D has same dimensions
            reshaped_finearray = zoom(finearray, (len(gridpts_mu)/length, len(gridpts_mu)/length,len(gridpts_mu)/length), order=1)
            #print('I have been reshaped ',reshaped_finearray.shape)
            reshaped_finearray=np.nan_to_num(reshaped_finearray) 
            Hist4D_splined[:,:,:,test_index]=reshaped_finearray
    #print(Hist4D_splined)
    return Hist4D_splined


flavours = ['NuE','NuMu','NuTau']


angles= np.linspace(0,1.0,6)

upp = np.linspace(1.40,2.0,2)
low = np.linspace(2.1,2.5,2)
depths = np.concatenate([upp,low])

for flavour in flavours:
    for angle in angles:
        for depth in depths:
            filename = '/home/zrechav/scripts/air_shower_reader/triple_histogram_making/histograms/' + flavour + '_Triple_Zen_' + str(np.around(angle,2)) + '_Depth_' + str(depth) + '.npy'
            print(filename)
            Splined_3D_Histogram = np.empty((len(gridpts_mu),len(gridpts_mu),len(gridpts_mu),len(E_nu_bins)))
            
            if os.path.exists(filename):
                try:
                    Splined_3D_Histogram = Make_3D_Spline(filename)

                    outfile = '/home/zrechav/scripts/air_shower_reader/spline_fitting/splines/3D/' + flavour + '_Triple_Splined_Zen_' +str(np.around(angle,2)) + '_Depth_' + str(depth) + '.npy' 
                    np.save(outfile, Splined_3D_Histogram)
                except Exception as e:
                    print(flavour,angle,depth)
                    print(e)
                    continue
        
