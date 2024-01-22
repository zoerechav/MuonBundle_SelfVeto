#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/vbasu/meta-projects/combo2/build
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib import colors as clrs
import matplotlib as mpl
import os


script ='/home/zrechav/scripts/air_shower_reader/single_histogram_making/SingleMuonHistograms.py'

low_angles = np.linspace(0,0.5,11)
high_angles = np.linspace(0.5,1,11)
angles_space = np.concatenate( (low_angles,high_angles) )

upp = np.linspace(1.40,2.0,12 + 1)
low = np.linspace(2.1,2.5,8 +1)
depth_space = np.concatenate([upp,low])

flavours = ['NuE','NuMu','NuTau']
count = 0
if os.path.exists(script):
    print("File exists!")
    try:
        for flavour in flavours:
        #print(flavour)
            for angle in angles_space:
                #print(angle)
                for depth in depth_space:
                    if count > 5:
                        break
                    count +=1
                  
                    #print(depth)
                    command = command = f"{script} --flavour {flavour} --angle {np.around(angle, 2)} --depth {np.around(depth, 2)}"
                    print(command)
                    os.system(command)
                    print('one down')
    except Exception as e:
        print(e)
else:
    print("File does not exist.")

            #os.system("python myOtherScript.py arg1 arg2 arg3")