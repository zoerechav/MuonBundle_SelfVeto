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
import pylab
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, CubicSpline, pchip_interpolate, interpn
import scipy.special as sp
from scipy import integrate
from random import seed
import sys, math
import copy
import glob
import numpy as np

import collections
import os

import random
# seed(1234)
from photospline import glam_fit, ndsparse, bspline
from mpl_toolkits import mplot3d
from scipy.optimize import curve_fit
def polygaussian(x,c, gauss_mean, gauss_sigma,d ):#
    normpdf=c*np.exp(-(x-gauss_mean)**2/(2*gauss_sigma**2))-d
    return normpdf
mult_bins=np.logspace(0,5,20+1)


low_angles= np.linspace(0,0.5,11)
high_angles=np.linspace(0.5,1,11)
new_low_angles=np.delete(low_angles,0) 
angles_space=np.concatenate( (new_low_angles,high_angles) )
E_nu_bins=np.logspace(1,8,10+1)
nu_bin_centers = np.sqrt(E_nu_bins[:-1] * E_nu_bins[1:])
E_mu_bins=np.logspace(-1,8,10+1)
mu_bin_centers = np.sqrt(E_mu_bins[:-1] * E_mu_bins[1:])
#test_index=4
#test_nu_energy=E_nu_bins[test_index]
ctr=0
import argparse

# handling of command line arguments  
from optparse import OptionParser
parser = OptionParser()
usage = """%prog [options]"""
parser.set_usage(usage)

parser.add_option("-f", "--flavour", action="store", type="string", default='NuMu', dest="FLAVOUR", help="Flavour")
parser.add_option("-a", "--angle", action="store", type="float", default=0.35, dest="ANGLE", help="Input angle")
parser.add_option("-d", "--depth", action="store", type="float", default=1.95, dest="DEPTH", help="Inputdepth")


# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()

gridpts=np.logspace(-1,7,20)
import os
plotx=np.log10(gridpts)
from scipy.interpolate import griddata
# for angle in angles_space:
#     for depth in depth_space:
coszen=np.around(options.ANGLE,decimals=2)
depth=np.around(options.DEPTH,decimals=2)
flavour=options.FLAVOUR
#filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_EnergyTotal_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'

# filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_Single_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'
# filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_Quadruple_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'
filename = '/data/user/vbasu/SelfVetoArrays/'+flavour+'_Double_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'
# filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_Triple_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'

# filename='/data/user/vbasu/SelfVetoArrays/Quintuple_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'
print(filename)
#test_index=4

if os.path.exists(filename):
    Hist2D=np.load(filename, mmap_mode="r")
    print(np.shape(Hist2D))
    Hist2D_splined=np.empty((len(gridpts),len(E_nu_bins)))
    for test_index in range(len(E_nu_bins)-1):
        
        test_nu_energy=E_nu_bins[test_index]
        
        print('Test Energy:',np.around(test_nu_energy,decimals=2),'GeV')
        Hist1D=Hist2D[:,test_index]
        print(Hist1D)
        Hist1D=Hist1D/np.sum(Hist1D)
        print(Hist1D)
        xv=(mu_bin_centers)
        mask=np.where(Hist1D!=0)

        try:
            grid_z1 = griddata((xv[mask]), Hist1D[mask], (plotx), method='linear')
            grid_z1=np.nan_to_num(grid_z1)
            grid_z1=grid_z1/np.sum(grid_z1)

            Hist2D_splined[:,test_index]=grid_z1.T
        except Exception as e:
            finearray=(Hist1D.T).repeat(2, 0)
    #                 print(np.shape(finearray))
            Hist2D_splined[:,test_index]=finearray
    filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_EnergyTotalSplined_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'
    print('Zen_'+str(coszen)+'_Depth_'+str(depth))
    np.save(filename, Hist2D_splined)
# 