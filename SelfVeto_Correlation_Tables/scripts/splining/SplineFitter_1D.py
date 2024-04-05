#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/vbasu/meta-projects/combo2/build

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib import colors as clrs
import matplotlib as mpl
import argparse
import os
import yaml
import sys
 
parser = argparse.ArgumentParser()


with open('/home/zrechav/SelfVeto_Correlation_Tables/scripts/config.yaml', 'r') as yaml_file:
    config = yaml.safe_load(yaml_file)


E_nu_bins=eval(config['E_nu_bins'])
nu_bin_centers = eval(config['nu_bin_centers'])
E_mu_bins=eval(config['E_mu_bins'])
mu_bin_centers = eval(config['mu_bin_centers'])

gridpts_mu = eval(config['gridpts_mu'])
gridpts_nu = eval(config['gridpts_nu'])

print(gridpts_mu)

def Make_1D_Spline(HistND_str,
                     E_mu_bins = E_mu_bins,
                     E_nu_bins = E_nu_bins,
                     gridpts_mu = gridpts_mu,
                     gridpts_nu = gridpts_nu,
                     spline_method = 'linear'):
     #HistND_str: string of 2D *.npy file location
     #E_mu_bins: muon energy bins
     #E_nu_bins: neutrino energy bins
     #gridpts_mu: muon energy splining bins
     #gridpts_nu: neutrino energy splining bins
     #spline_method: griddata spline method input: linear, cubic, nearest
    try:
        HistND = np.load(HistND_str, mmap_mode='r')
        HistND = np.nan_to_num(HistND)

        grid_dims = (E_nu_bin_centers, E_mu_bin_centers)
        grid_arrays = np.meshgrid(gridpts_nu_centers, gridpts_mu_centers, indexing='ij')

        interp = RegularGridInterpolator(grid_dims, HistND, method=spline_method, bounds_error=False, fill_value=0)
        interpolated_data = interp(np.stack(grid_arrays, axis=-1))

        return interpolated_data

    except Exception as e:
        print(e)

# def Make_1D_Spline(  Hist2D_str,
#                      E_mu_bins = E_mu_bins,
#                      E_nu_bins = E_nu_bins,
#                      gridpts_mu = gridpts_mu,
#                      gridpts_nu = gridpts_nu,
#                      spline_method = 'linear',
                     
#                     ) :
#     #Hist2D_str: string of .npy file location
#     #E_mu_bins: muon energy bins
#     #E_nu_bins: neutrino energy bins
#     #gridpts_mu: muon energy splining bins
#     #gridpts_nu: neutrino energy splining bins
#     #spline_method: griddata spline method input: linear, cubic, nearest
 
#     Hist2D = np.load(Hist2D_str,mmap_mode='r')
#     print(Hist2D.shape)
#     try:
#         mask=np.where(Hist2D!=0)

#         ##Spline Linear Interpolation##
#         xv, yv=np.meshgrid((E_nu_bins),(E_mu_bins))
#         plotx,ploty=np.meshgrid((gridpts_nu),(gridpts_mu))
#         #print(plotx)
#         print(ploty)

#         grid_z1 = griddata((xv[mask],yv[mask]), Hist2D[mask], (plotx,ploty), method=spline_method)
#         print(grid_z1.shape)
#         grid_z1 = np.nan_to_num(grid_z1)
#         Hist2D_splined = grid_z1
#         return Hist2D_splined
#     except Exception as e:
#         #print(flavour,angle,depth)
#         print(e)
#         print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e).__name__, e)
#         #continue

#         # Normalize Spline
#     #     grid_z1 = grid_z1 / grid_z1.sum(axis=1,keepdims=True)
#     #     grid_z1 = np.nan_to_num(grid_z1)

        
    #print(Hist2D_splined)
    

flavours = config['flavours']
angles= eval(config['angles'])
depths = eval(config['depths'])

for flavour in flavours:
    for angle in angles:
        for depth in depths:
            filename = config['single_hists_base'] + flavour + '_Single_Zen_' + str(np.around(angle,2)) + '_Depth_' + str(depth) + '.npy'
            
            Splined_Histogram = None
            if os.path.exists(filename):
                try:
                    Splined_Histogram = Make_1D_Spline(filename)
                    print(Splined_Histogram.shape)
                    print(' I was splined')
                    outfile = config['single_splines_base'] + flavour + '_Single_Splined_Zen_' +str(np.around(angle,2)) + '_Depth_' + str(depth) + '.npy' 
                    np.save(outfile, Splined_Histogram)
                except Exception as e:
                    print(flavour,angle,depth)
                    print(e)
                    print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e).__name__, e)
                    continue
        
