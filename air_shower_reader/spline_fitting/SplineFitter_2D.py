#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/vbasu/meta-projects/combo3/build


import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.ndimage import zoom
from matplotlib import colors as clrs
import matplotlib as mpl
import os
import argparse  
from optparse import OptionParser


parser = OptionParser()
usage = """%prog [options]"""
parser.set_usage(usage)

parser.add_option("-f", "--flavour", action="store", type="string", default='NuMu', dest="FLAVOUR", help="Flavour")
parser.add_option("-a", "--angle", action="store", type="float", default=0.17, dest="ANGLE", help="Input angle")
parser.add_option("-d", "--depth", action="store", type="float", default=1.5, dest="DEPTH", help="Inputdepth")


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

def make_spline_for_2D_bin(
                    Hist3D,
                    test_index,
                    Hist3D_splined,
                    E_mu_bins = np.logspace(1,7,12+1),
                    E_nu_bins = np.logspace(1,7,12+1),
                    gridpts_mu = np.logspace(1,7,31+1),
                    gridpts_nu = np.logspace(1,7,31+1),
                    spline_method = 'linear',
                     
                    ) :    

    #Hist3D: Histogram for 2D case
    #test_index: neutrino energy bin
    #Hist3D_splined: 3D Histogram being filled
    #E_mu_bins: muon energy bins
    #E_nu_bins: neutrino energy bins
    #gridpts_mu: muon energy splining bins
    #gridpts_nu: neutrino energy splining bins
    #spline_method: griddata spline method input: linear, cubic, nearest
    Hist3D_splined = Hist3D_splined
    test_nu_energy = E_nu_bins[test_index]
    Hist2D=Hist3D[:,:,test_index]
    Hist2D=Hist2D/np.sum(Hist2D) ##normalizing the histogram 
    xv, yv=np.meshgrid(mu_bin_centers, mu_bin_centers)
    plotx,ploty=np.meshgrid(gridpts_mu,gridpts_mu)
    mask=np.where(Hist2D!=0)
    

    try:
        grid_z1 = griddata((xv[mask],yv[mask]), Hist2D[mask], (plotx,ploty), method='linear')
        grid_z1=np.nan_to_num(grid_z1)
        grid_z1=grid_z1/np.sum(grid_z1) ##normalizing the splined hist

        Hist3D_splined[:,:,test_index]=grid_z1.T ##adding modified array to HIST3D
    except Exception as e:
        print('splining didnt work')
        print(e)
        ##instead of splining, just make hist more dense
        finearray=(Hist2D.T).repeat(2, 0).repeat(2, 1)
        length = finearray.shape[0]
        # Resize the array to shape (32, 32) using zoom, so Hist3D has same dimensions
        reshaped_finearray = zoom(finearray, (len(gridpts_mu)/length, len(gridpts_mu)/length), order=1)
        print(reshaped_finearray.shape)
        reshaped_finearray=np.nan_to_num(reshaped_finearray) 
        Hist3D_splined[:,:,test_index]=reshaped_finearray
        print('I have been reshaped') 
    return Hist3D_splined

def Make_2D_Spline( Hist3D_str,
                     E_mu_bins = np.logspace(1,7,12+1),
                     E_nu_bins = np.logspace(1,7,12+1),
                     gridpts_mu = np.logspace(1,7,31+1),
                     gridpts_nu = np.logspace(1,7,31+1),
                     spline_method = 'linear',
                     
                    ) :    

    #Hist3D_str: string of .npy file location
    #E_mu_bins: muon energy bins
    #E_nu_bins: neutrino energy bins
    #gridpts_mu: muon energy splining bins
    #gridpts_nu: neutrino energy splining bins
    #spline_method: griddata spline method input: linear, cubic, nearest 
    Hist3D = np.load(Hist3D_str)
    Hist3d_splined=np.empty((len(gridpts_mu),len(gridpts_mu),len(E_nu_bins)))
    E_mu_bins=E_mu_bins
    E_nu_bins=E_nu_bins
    mu_bin_centers = np.sqrt(E_mu_bins[:-1] * E_mu_bins[1:])
    nu_bin_centers = np.sqrt(E_nu_bins[:-1] * E_nu_bins[1:])

    for test_index in range(len(E_nu_bins)-1):
        make_spline_for_2D_bin(Hist3D,test_index,Hist3d_splined)
        
#        # Plot the result for each test_index
#         plt.pcolormesh(gridpts_mu, gridpts_nu, Hist3d_splined[:, :, test_index],
#                        cmap='turbo', norm=clrs.LogNorm(vmin=10**-4, vmax=1))
#         plt.title(f'Index {test_index} in E_nu_bins')
#         plt.xlabel('Muon 1 Energy')
#         plt.ylabel('Muon 2 Energy')
#         plt.xscale('log')
#         plt.yscale('log')
#         plt.colorbar()
#         plt.show()
        
    return (Hist3d_splined)

flavours = ['NuE','NuMu','NuTau']


angles= np.linspace(0,1.0,6)

upp = np.linspace(1.40,2.0,2)
low = np.linspace(2.1,2.5,2)
depths = np.concatenate([upp,low])


for flavour in flavours:
    for angle in angles:
        for depth in depths:
            filename = '/home/zrechav/scripts/air_shower_reader/double_histogram_making/histograms/' + flavour + '_Double_Zen_' + str(np.around(angle,2)) + '_Depth_' + str(depth)+'.npy'
            print(filename)

            Splined_2D_Histogram = np.empty((len(gridpts_mu),len(gridpts_mu),len(E_nu_bins)))
            
            if os.path.exists(filename):
                try:
                    Splined_2D_Histogram = Make_2D_Spline(filename)
                    outfile = '/home/zrechav/scripts/air_shower_reader/spline_fitting/splines/2D/' + flavour + '_DoubleSpline_Zen_' + str(np.around(angle,2)) + '_Depth_' + str(depth) + '.npy'
                    print('Zen_'+str(np.around(angle,2))+'_Depth_'+str(depth))
                    np.save(outfile, Splined_2D_Histogram)
                except Exception as e:
                    print(flavour,angle,depth)
                    print(e)
                    continue
