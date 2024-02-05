#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/vbasu/meta-projects/combo2/build

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib import colors as clrs
import matplotlib as mpl
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




def Make_1D_Spline(  Hist2D_str,
                     E_mu_bins = np.logspace(1,7,12+1),
                     E_nu_bins = np.logspace(1,7,12+1),
                     gridpts_mu = np.logspace(1,7,31+1),
                     gridpts_nu = np.logspace(1,7,31+1),
                     spline_method = 'linear',
                     
                    ) :
    #Hist2D_str: string of .npy file location
    #E_mu_bins: muon energy bins
    #E_nu_bins: neutrino energy bins
    #gridpts_mu: muon energy splining bins
    #gridpts_nu: neutrino energy splining bins
    #spline_method: griddata spline method input: linear, cubic, nearest
 
    Hist2D = np.load(Hist2D_str,mmap_mode='r')
    E_mu_bins=E_mu_bins
    E_nu_bins=E_nu_bins
    mu_bin_centers = np.sqrt(E_mu_bins[:-1] * E_mu_bins[1:])
    nu_bin_centers = np.sqrt(E_nu_bins[:-1] * E_nu_bins[1:])
    
    mask=np.where(Hist2D!=0)

    ##Spline Linear Interpolation##
    xv, yv=np.meshgrid((mu_bin_centers),(mu_bin_centers))
    gridpts_nu= gridpts_mu
    gridpts_mu = gridpts_mu
    plotx,ploty=np.meshgrid((len(gridpts_mu)),len(gridpts_mu))
    grid_z1 = griddata((xv[mask],yv[mask]), Hist2D[mask], (plotx,ploty), method=spline_method)
    grid_z1 = np.nan_to_num(grid_z1) ##why this step? (eliminating nans in arrays)
    grid_z1 = grid_z1/np.sum(grid_z1) ##why this step? (normalizing the splined histogram)
    Hist2D_splined = grid_z1.T
    return Hist2D_splined

flavours = ['NuE','NuMu','NuTau']


angles= np.linspace(0,1.0,6)

upp = np.linspace(1.40,2.0,2)
low = np.linspace(2.1,2.5,2)
depths = np.concatenate([upp,low])

for flavour in flavours:
    for angle in angles:
        for depth in depths:
            filename = '/home/zrechav/scripts/air_shower_reader/single_histogram_making/histograms/' + flavour + '_Single_Zen_' + str(np.around(angle,2)) + '_Depth_' + str(depth) + '.npy'
            print(filename)
            Splined_Histogram = None
            if os.path.exists(filename):
                try:
                    Splined_Histogram = Make_1D_Spline(filename)

                    outfile = '/home/zrechav/scripts/air_shower_reader/spline_fitting/splines/1D/' + flavour + '_Single_Splined_Zen_' +str(np.around(angle,2)) + '_Depth_' + str(depth) + '.npy' 
                    np.save(outfile, Splined_Histogram)
                except Exception as e:
                    print(flavour,angle,depth)
                    print(e)
                    continue
        
