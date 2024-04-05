#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/vbasu/meta-projects/combo3/build


import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.ndimage import zoom
from matplotlib import colors as clrs
import matplotlib as mpl
import yaml

import argparse
import os
 
parser = argparse.ArgumentParser()


with open('/home/zrechav/SelfVeto_Correlation_Tables/scripts/config.yaml', 'r') as yaml_file:
    config = yaml.safe_load(yaml_file)


E_nu_bins=eval(config['E_nu_bins'])
nu_bin_centers = eval(config['nu_bin_centers'])
E_mu_bins=eval(config['E_mu_bins'])
mu_bin_centers = eval(config['mu_bin_centers'])

gridpts_mu = eval(config['gridpts_mu'])
gridpts_nu = eval(config['gridpts_nu'])

def make_spline_for_2D_bin(
                    Hist3D,
                    test_index,
                    Hist3D_splined,
                    E_mu_bins = E_mu_bins,
                    E_nu_bins = E_nu_bins,
                    gridpts_mu = gridpts_mu,
                    gridpts_nu = gridpts_nu,
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
   
    xv, yv=np.meshgrid(E_mu_bins, E_mu_bins)
    plotx,ploty=np.meshgrid(gridpts_mu,gridpts_mu)
    mask=np.where(Hist2D!=0)
    #print(Hist2D)

    try:
        grid_z1 = griddata((xv[mask],yv[mask]), Hist2D[mask], (plotx,ploty), method='linear')
        grid_z1=np.nan_to_num(grid_z1)
        ##new normalization line
        #grid_z1=grid_z1/np.sum(grid_z1) ##normalizing the splined hist
        
        print(grid_z1)
        print('i DID work')
        Hist3D_splined[:,:,test_index]=grid_z1 ##adding modified array to HIST3D
    except Exception as e:
        print('splining didnt work')
        print(e)
        ##instead of splining, just make hist more dense
        finearray=(Hist2D).repeat(2, 0).repeat(2, 1)
        length = finearray.shape[0]
        # Resize the array to shape (32, 32) using zoom, so Hist3D has same dimensions
        reshaped_finearray = zoom(finearray, (len(gridpts_mu)/length, len(gridpts_mu)/length), order=1)
        #print(reshaped_finearray.shape)
        reshaped_finearray=np.nan_to_num(reshaped_finearray) 
        Hist3D_splined[:,:,test_index]=reshaped_finearray
        print('I have been reshaped') 
    return Hist3D_splined

def Make_2D_Spline( Hist3D_str,
                     E_mu_bins = E_mu_bins,
                     E_nu_bins = E_nu_bins,
                     gridpts_mu = gridpts_mu,
                     gridpts_nu = gridpts_nu,
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

    for test_index in range(len(E_nu_bins)-1):
        make_spline_for_2D_bin(Hist3D,test_index,Hist3d_splined)
        
        
    return (Hist3d_splined)

flavours = config['flavours']
angles= eval(config['angles'])
depths = eval(config['depths'])


for flavour in flavours:
    for angle in angles:
        for depth in depths:
            filename = config['double_hists_base'] + flavour + '_Double_Zen_' + str(np.around(angle,2)) + '_Depth_' + str(depth)+'.npy'
            print(filename)

            Splined_2D_Histogram = np.empty((len(gridpts_mu),len(gridpts_mu),len(E_nu_bins)))
            
            if os.path.exists(filename):
                try:
                    Splined_2D_Histogram = Make_2D_Spline(filename)
                    outfile = config['double_splines_base'] + flavour + '_DoubleSpline_Zen_' + str(np.around(angle,2)) + '_Depth_' + str(depth) + '.npy'                    
                    np.save(outfile, Splined_2D_Histogram)
                except Exception as e:
                    print(flavour,angle,depth)
                    print(e)
                    continue
