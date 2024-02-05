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


low_angles= np.linspace(0,0.5,13)
high_angles=np.linspace(0.5,1,13)
new_low_angles=np.delete(low_angles,0) 
angles_space=np.concatenate( (new_low_angles,high_angles) )
E_nu_bins=np.logspace(1,7,10+1)
nu_bin_centers = np.sqrt(E_nu_bins[:-1] * E_nu_bins[1:])
E_mu_bins=np.logspace(-1,7,10+1)
mu_bin_centers = np.sqrt(E_mu_bins[:-1] * E_mu_bins[1:])
test_index=4
test_nu_energy=E_nu_bins[test_index]
ctr=0
import argparse

# handling of command line arguments  
from optparse import OptionParser
parser = OptionParser()
usage = """%prog [options]"""
parser.set_usage(usage)

parser.add_option("-f", "--flavour", action="store", type="string", default='NuMu', dest="FLAVOUR", help="Flavour")
parser.add_option("-a", "--angle", action="store", type="float", default=0.83, dest="ANGLE", help="Input angle")
parser.add_option("-d", "--depth", action="store", type="float", default=1.66, dest="DEPTH", help="Inputdepth")


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
filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_EnergyTotal_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'

# filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_Single_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'
# filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_Quadruple_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'
# filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_Double_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'
# filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_Triple_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'

# filename='/data/user/vbasu/SelfVetoArrays/Quintuple_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'
print(filename)
test_index=4
# if os.path.exists(filename):
# 	Hist5D=np.load(filename, mmap_mode="r")
# 	print(np.shape(Hist5D))
# 	# Hist5D.close()
# 	Hist5D_splined=np.empty((len(gridpts),len(gridpts),len(gridpts),len(gridpts),len(E_nu_bins)))
# 	for test_index in range(len(E_nu_bins)-1):
# 	    test_nu_energy=E_nu_bins[test_index]
# 	    print('Test Energy:',np.around(test_nu_energy,decimals=2),'GeV')
# 	    Hist4D=Hist5D[:,:,:,:,test_index]
# 	    Hist4D=Hist4D/np.sum(Hist4D)
# 	    xv, yv,zv,tv=np.meshgrid(np.log10(mu_bin_centers),np.log10(mu_bin_centers),np.log10(mu_bin_centers),np.log10(mu_bin_centers))
# 	    mask=np.where(Hist4D!=0)

# 	    try:
# 	        grid_z1 = griddata((xv[mask],yv[mask],zv[mask],tv[mask]), Hist4D[mask], (plotx,ploty,plotz,plott), method='linear')
# 	        grid_z1=np.nan_to_num(grid_z1)
# 	        grid_z1=grid_z1/np.sum(grid_z1)

# 	        Hist5D_splined[:,:,:,:,test_index]=grid_z1.T
# 	    except Exception as e:
# 	        finearray=(Hist4D.T).repeat(2, 0).repeat(2, 1).repeat(2, 2).repeat(2,3)
# 	#                 print(np.shape(finearray))
# 	        Hist5D_splined[:,:,:,:,test_index]=finearray
# 	filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_QuadrupleSpline_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'
# 	# filename='SelfVeto/QuinSplines/QuintupleSpline_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'
# 	print('Zen_'+str(coszen)+'_Depth_'+str(depth))
# 	np.save(filename, Hist5D_splined)
# if os.path.exists(filename):
# 	Hist4D=np.load(filename, mmap_mode="r")
# 	print(np.shape(Hist4D))
# 	# Hist5D.close()
# 	Hist4D_splined=np.empty((len(gridpts),len(gridpts),len(gridpts),len(E_nu_bins)))
# 	for test_index in range(len(E_nu_bins)-1):
# 	    test_nu_energy=E_nu_bins[test_index]
# 	    print('Test Energy:',np.around(test_nu_energy,decimals=2),'GeV')
# 	    Hist3D=Hist4D[:,:,:,test_index]
# 	    Hist3D=Hist3D/np.sum(Hist3D)
# 	    xv, yv,zv=np.meshgrid(np.log10(mu_bin_centers),np.log10(mu_bin_centers),np.log10(mu_bin_centers))
# 	    mask=np.where(Hist3D!=0)

# 	    try:
# 	        grid_z1 = griddata((xv[mask],yv[mask],zv[mask]), Hist3D[mask], (plotx,ploty,plotz), method='linear')
# 	        grid_z1=np.nan_to_num(grid_z1)
# 	        grid_z1=grid_z1/np.sum(grid_z1)

# 	        Hist4D_splined[:,:,:,test_index]=grid_z1.T
# 	    except Exception as e:
# 	        finearray=(Hist3D.T).repeat(2, 0).repeat(2, 1).repeat(2, 2)
# 	#                 print(np.shape(finearray))
# 	        Hist4D_splined[:,:,:,test_index]=finearray
# 	filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_TripleSpline_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'
# 	print('Zen_'+str(coszen)+'_Depth_'+str(depth))
# 	np.save(filename, Hist4D_splined)

# if os.path.exists(filename):
# 	Hist3D=np.load(filename, mmap_mode="r")
# 	print(np.shape(Hist3D))
# 	# Hist5D.close()
# 	Hist3D_splined=np.empty((len(gridpts),len(gridpts),len(E_nu_bins)))
# 	for test_index in range(len(E_nu_bins)-1):
# 	    test_nu_energy=E_nu_bins[test_index]
# 	    print('Test Energy:',np.around(test_nu_energy,decimals=2),'GeV')
# 	    Hist2D=Hist3D[:,:,test_index]
# 	    Hist2D=Hist2D/np.sum(Hist2D)
# 	    xv, yv=np.meshgrid(np.log10(mu_bin_centers),np.log10(mu_bin_centers))
# 	    mask=np.where(Hist2D!=0)

# 	    try:
# 	        grid_z1 = griddata((xv[mask],yv[mask]), Hist2D[mask], (plotx,ploty), method='linear')
# 	        grid_z1=np.nan_to_num(grid_z1)
# 	        grid_z1=grid_z1/np.sum(grid_z1)

# 	        Hist3D_splined[:,:,test_index]=grid_z1.T
# 	    except Exception as e:
# 	        finearray=(Hist2D.T).repeat(2, 0).repeat(2, 1)
# 	#                 print(np.shape(finearray))
# 	        Hist3D_splined[:,:,test_index]=finearray
# 	filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_DoubleSpline_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'
# 	print('Zen_'+str(coszen)+'_Depth_'+str(depth))
# 	np.save(filename, Hist3D_splined)
# if os.path.exists(filename):
# 	FaultyArrays=[]
# 	Hist2D=np.load(filename, mmap_mode="r")
# 	print(np.shape(Hist2D))
# 	# Hist5D.close()
# 	for test_index in range(len(mu_bin_centers)):
#         try:
#             test_index=int(test_index)
#             test_nu_energy=E_nu_bins[test_index]
#             if np.isnan(Hist2D.T[:][test_index]).all():continue            
#             xdata,ydata=mu_bin_centers[Hist2D.T[:][test_index]!=0],Hist2D.T[:][test_index][Hist2D.T[:][test_index]!=0]

#             SampleMu=np.logspace(0,7,16+1)
 
#             logx=np.log10(xdata)
#             logy=np.log10(ydata)
#             n=len(xdata)
#             sigma = np.abs(sum(logy*(logx-np.median(np.log10(mu_bin_centers)))**2)/n)      #note this correction
#             popt_mod, pcov = curve_fit(polygaussian, logx, logy,
#                                        p0=[15,np.median(np.log10(mu_bin_centers)),sigma,10], 
#                                        method='dogbox',
#                                        bounds=([0,1,0.01,5],[15, 10,100,30])
#                                       )
#             fitPDF = polygaussian(np.log10(mu_bin_centers), *popt_mod                                  
#                                  )
#             print(popt_mod)
#             print('X values:',logx)
#             print('Target values:',logy)
#             print('Fit values:',fitPDF)
# #             print()
#             filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_SingleSplined_Zen_'+str(coszen)+'_Depth_'+str(depth)+'_Enu_'+str(np.around(test_nu_energy,decimals=2))+'.npy'
#             print('Zen_'+str(coszen)+'_Depth_'+str(depth)+'_Enu_'+str(np.around(test_nu_energy,decimals=2))+' '+str(ctr))
#             np.save(filename, popt_mod)

# #             plt.plot(mu_bin_centers,smooth(Hist2D.T[:][test_index],2),linestyle='dashdot',label=r'Smoothed', color=colors[ctr%4])
#             ctr+=1
#         except Exception as e:
#             filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_SingleSplined_Zen_'+str(coszen)+'_Depth_'+str(depth)+'_Enu_'+str(np.around(test_nu_energy,decimals=2))+'.npy'
#             FaultyArrays.append(filename)
#             if np.isnan(Hist2D.T[:][test_index]).all():continue
#             print(Hist2D.T[:][test_index])
#             zero_index=np.nonzero(Hist2D.T[:][test_index])[0]
#             e_left=np.log10(E_mu_bins[zero_index])
#             e_right=np.log10(E_mu_bins[(zero_index+1)])
#             e_center=np.log10(mu_bin_centers[zero_index])
#             trialsigma=0.5*(e_right-e_left)
#             popt=np.asarray([1,e_center,trialsigma,1])
#             print(filename)
#             np.save(filename, popt)
#     print(len(FaultyArrays))
# 	print(FaultyArrays)

# if os.path.exists(filename):
#     FaultyArrays=[]
#     Hist2D=np.load(filename, mmap_mode="r")
#     print(np.shape(Hist2D))
#     popt_array=[]
#     # Hist5D.close()
#     for test_index in range(len(mu_bin_centers)):
#         try:
#             test_index=int(test_index)
#             test_nu_energy=E_nu_bins[test_index]
#             if np.isnan(Hist2D.T[:][test_index]).all():continue            
#             xdata,ydata=mu_bin_centers[Hist2D.T[:][test_index]!=0],Hist2D.T[:][test_index][Hist2D.T[:][test_index]!=0]

#             SampleMu=np.logspace(0,7,16+1)
 
#             logx=np.log10(xdata)
#             logy=np.log10(ydata)
#             n=len(xdata)
#             sigma = np.abs(sum(logy*(logx-np.median(np.log10(mu_bin_centers)))**2)/n)      #note this correction
#             popt_mod, pcov = curve_fit(polygaussian, logx, logy,
#                                        p0=[15,np.median(np.log10(mu_bin_centers)),sigma,10], 
#                                        method='dogbox',
#                                        bounds=([0,1,0.01,5],[15, 10,100,30])
#                                       )
#             fitPDF = polygaussian(np.log10(mu_bin_centers), *popt_mod                                  
#                                  )
#             print(popt_mod)
#             print('X values:',logx)
#             print('Target values:',logy)
#             print('Fit values:',fitPDF)
# #             print()
#             # filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_EnergyTotalSplined_Zen_'+str(coszen)+'_Depth_'+str(depth)+'_Enu_'+str(np.around(test_nu_energy,decimals=2))+'.npy'
#             # print('Zen_'+str(coszen)+'_Depth_'+str(depth)+'_Enu_'+str(np.around(test_nu_energy,decimals=2))+' '+str(ctr))
#             # np.save(filename, popt_mod)
#             popt_array.append(popt_mod)

# #             plt.plot(mu_bin_centers,smooth(Hist2D.T[:][test_index],2),linestyle='dashdot',label=r'Smoothed', color=colors[ctr%4])
#             ctr+=1
#         except Exception as e:
#             filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_EnergyTotalSplined_Zen_'+str(coszen)+'_Depth_'+str(depth)+'_Enu_'+str(np.around(test_nu_energy,decimals=2))+'.npy'
#             FaultyArrays.append(filename)
#             if np.isnan(Hist2D.T[:][test_index]).all():continue
#             print(Hist2D.T[:][test_index])
#             zero_index=np.nonzero(Hist2D.T[:][test_index])[0]
#             e_left=np.log10(E_mu_bins[zero_index])
#             e_right=np.log10(E_mu_bins[(zero_index+1)])
#             e_center=np.log10(mu_bin_centers[zero_index])
#             trialsigma=0.5*(e_right-e_left)
#             if (len(e_center)>0 and len(trialsigma>0)):
#                 popt=np.asarray([1,e_center[0],trialsigma[0],1])
#                 print('Faulty Files')
#                 print(filename)
#                 print(popt)
#                 print()
#                 popt_array.append(popt)
#             else:
#                 print('DoubleFault')
#                 print('E_left:',e_left)
#                 print('E_center:',e_center)
#                 print('E_right:',e_right)
#                 print(zero_index)
#                 popt_array.append(np.array([]))
#                 print()
#     print(len(FaultyArrays))
#     print(FaultyArrays)
#     filename='/data/user/vbasu/SelfVetoArrays/'+flavour+'_EnergyTotalSplined_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'
#     print('Zen_'+str(coszen)+'_Depth_'+str(depth))
#     np.save(filename, np.array(popt_array))
#     print(popt_array)   
if os.path.exists(filename):
    Hist2D=np.load(filename, mmap_mode="r")
    print(np.shape(Hist2D))
    # Hist5D.close()
    Hist2D_splined=np.empty((len(gridpts),len(E_nu_bins)))
    for test_index in range(len(E_nu_bins)-1):
        test_nu_energy=E_nu_bins[test_index]
        print('Test Energy:',np.around(test_nu_energy,decimals=2),'GeV')
        Hist1D=Hist2D[:,test_index]
        Hist1D=Hist1D/np.sum(Hist1D)
        xv=np.log10(mu_bin_centers)
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