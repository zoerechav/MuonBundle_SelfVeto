#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/vbasu/meta-projects/combo2/build

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib import colors as clrs
import matplotlib as mpl
from scipy.ndimage import zoom
import argparse
import os
import yaml
parser = argparse.ArgumentParser()


with open('/home/zrechav/SelfVeto_Correlation_Tables/scripts/config.yaml', 'r') as yaml_file:
    config = yaml.safe_load(yaml_file)


E_nu_bins=eval(config['E_nu_bins'])
nu_bin_centers = eval(config['nu_bin_centers'])
E_mu_bins=eval(config['E_mu_bins'])
mu_bin_centers = eval(config['mu_bin_centers'])

gridpts_mu = eval(config['gridpts_mu'])
gridpts_nu = eval(config['gridpts_nu'])


def Make_3D_Spline(  filename,
                     E_mu_bins = E_mu_bins,
                     E_nu_bins = E_nu_bins,
                     gridpts_mu = gridpts_mu,
                     gridpts_nu = gridpts_nu,
                     spline_method = 'linear',
                     
                    ) :    
    
    
    
    Hist4D=np.load(filename)
    
    Hist4D_splined=np.empty((len(gridpts_mu),len(gridpts_mu),len(gridpts_mu),len(E_nu_bins)))
    for test_index in range(len(E_nu_bins)-1):
        test_nu_energy=E_nu_bins[test_index]
        Hist3D=Hist4D[:,:,:,test_index]
        #Hist3D=Hist3D/np.sum(Hist3D)   
        Hist3D=np.nan_to_num(Hist3D)

        xv, yv,zv=np.meshgrid((mu_bin_centers),(mu_bin_centers),(mu_bin_centers))
        plotx,ploty,plotz=np.meshgrid((gridpts_mu),(gridpts_mu),(gridpts_mu))
        mask=np.where(Hist3D!=0)

        try:
            grid_z1 = griddata((xv[mask],yv[mask],zv[mask]), Hist3D[mask], (plotx,ploty,plotz), method=spline_method)      
            #grid_z1=grid_z1/np.sum(grid_z1)
            ##new line of normalization code
            #grid_z1 = grid_z1/np.sum(Hist4D[:,:,:,energy_index]
            grid_z1 = grid_z1 / np.sum(grid_z1, axis=-1, keepdims=True)
            grid_z1=np.nan_to_num(grid_z1)
            Hist4D_splined[:,:,:,test_index]=grid_z1
            
        except Exception as e:
            #print(e)
            finearray=(Hist3D).repeat(2, 0).repeat(2, 1).repeat(2,2)
            length = finearray.shape[0]
            # Resize the array to shape (32, 32) using zoom, so Hist3D has same dimensions
            reshaped_finearray = zoom(finearray, (len(gridpts_mu)/length, len(gridpts_mu)/length,len(gridpts_mu)/length), order=1)
            ##new line for normalization
            reshaped_finearray = reshaped_finearray / np.sum(reshaped_finearray, axis=-1, keepdims=True)
            reshaped_finearray=np.nan_to_num(reshaped_finearray) 
            Hist4D_splined[:,:,:,test_index]=reshaped_finearray
    return Hist4D_splined


flavours = config['flavours']
angles= eval(config['angles'])
depths = eval(config['depths'])

for flavour in flavours:
    for angle in angles:
        for depth in depths:
            filename = config['triple_hists_base'] + flavour + '_Triple_Zen_' + str(np.around(angle,2)) + '_Depth_' + str(depth) + '.npy'
            Splined_3D_Histogram = np.empty((len(gridpts_mu),len(gridpts_mu),len(gridpts_mu),len(E_nu_bins)))
            
            if os.path.exists(filename):
                try:
                    Splined_3D_Histogram = Make_3D_Spline(filename)

                    outfile = config['triple_splines_base'] + flavour + '_Triple_Splined_Zen_' +str(np.around(angle,2)) + '_Depth_' + str(depth) + '.npy' 
                    np.save(outfile, Splined_3D_Histogram)
                except Exception as e:
                    print(flavour,angle,depth)
                    print(e)
                    continue
        
