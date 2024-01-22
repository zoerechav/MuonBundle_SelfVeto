#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/vbasu/meta-projects/combo3/build


from I3Tray import *
from icecube import icetray, phys_services, dataio, dataclasses,MuonGun
from icecube.simprod import segments
import sys, math
import copy
import glob
import numpy as np
import nuflux
import collections
import os
import random
import ROOT

import time
from icecube.hdfwriter import I3HDFTableService, I3HDFWriter

import argparse
sys.path.append('/cvmfs/icecube.opensciencegrid.org/users/vbasu/scripts/NewIce/')
sys.path.append('/home/vbasu/scripts/')
import meseUtils
import MCWeightUtils_Airshowers

from icecube import dataclasses
mass_dict=    {'2212': 1, #PPlus
               '1000020040': 4,#Alpha
               '1000070140': 14,#N14
               '1000130270': 27,#Al27
               '1000260560': 56#Fe56
              }

start_time = time.asctime()
print ('Started:', start_time)

# handling of command line arguments  
from optparse import OptionParser
parser = OptionParser()
usage = """%prog [options]"""
parser.set_usage(usage)
corsika_filelist_length = 500

#parser.add_option("-i", "--input", action="store", type="string", default='/home/zrechav/test/test_mctree_prop.i3.zst', dest="INPUT", help="Output i3 file")
# parser.add_option("-o", "--output", action="store", type="string", default="/home/vbasu/testoutput/L2_selfveto", dest="OUTPUT", help="Output i3 file")

parser.add_option("-i", "--input", action="store", type="string", default="/data/user/zrechav/propagated_corsika/20904/0000000-0000999/propagated_corsika_20904_000000.i3.zst", dest="INPUT", help="Output i3 file")
parser.add_option("-o", "--output", action="store", type="string", default="/home/zrechav/test/AirShowerReader", dest="OUTPUT", help="Output i3 file")
parser.add_option("-g", "--gcd", action="store", type="string", default="/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz", dest="GCD", help="GCD file for input i3 file")

#parser.add_option("-f","--grid", action="store", type="int", default=0, dest="GridFTP", help="to use gridftp or not. 0 is false and 1 is true")

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()

infile = sorted(glob.glob(options.INPUT))
outfile = options.OUTPUT
gcdFile=options.GCD 
tray = I3Tray()

print ("Reading input file...",  options.INPUT)
tray.AddSegment(dataio.I3Reader, 'reader', FilenameList=infile)
#print ("Reading input file...",  options.INPUT)  
#print(options.GCD)


#########################################
## GATHER SOME MC TRUTH
#########################################   


NuTypes=[                    
                          dataclasses.I3Particle.ParticleType.NuE,
                          dataclasses.I3Particle.ParticleType.NuEBar,
                          dataclasses.I3Particle.ParticleType.NuMu,
                          dataclasses.I3Particle.ParticleType.NuMuBar,
                          dataclasses.I3Particle.ParticleType.NuTau,
                          dataclasses.I3Particle.ParticleType.NuTauBar,
                          ]
MuTypes=[         
                          dataclasses.I3Particle.ParticleType.MuPlus,
                          dataclasses.I3Particle.ParticleType.MuMinus,
                        
                          ]

import matplotlib.path as mpltPath
def get_surface_det(gcdFile=None):
    
    gcdFile=gcdFile
    #print(gcdFile,'IM HERE')
    bound_2D=[]
    surface_det = MuonGun.ExtrudedPolygon.from_file(gcdFile, padding=200)##Build Polygon from I3Geometry
    #surface_det = MuonGun.Cylinder(1000,700) ##Cylinder 1000m high with 700m radius
    #print(surface_det)
    x=[(surface_det.x[i],surface_det.y[i])for i in range(len(surface_det.x))]###getting only x and y
    bound_2D=mpltPath.Path(x)#Projection of detector on x,y plane
    #bound_2D = surface_det
    return bound_2D, surface_det

##call here to reduce computation time
bound_2D,surface_det = get_surface_det(gcdFile=gcdFile)

def boundary_check(particle1,gcdFile=None):
    ####checks if particle is inside the detector###
    inlimit = False  
    if ((particle1.pos.z <=max(surface_det.z)) and (particle1.pos.z>=min(surface_det.z))):
        if bound_2D.contains_points([(particle1.pos.x, particle1.pos.y)]):
            inlimit=True

    return inlimit


def NeutrinoSelector(frame):
    #print('IM in NEUTRINO SELECTOR')
    if 'I3MCTree' in frame and 'Homogenized_QTot' in frame:
        frame['I3MCTree_copy']=frame['I3MCTree']
        mctree = frame['I3MCTree_copy']
        #mctree = frame['I3MCTree_preMuonProp']
        e_dep_total=0
        neutrino = None
        neutrinos=[]
        primary=None
        for p in mctree:#search for neutrino
            if mctree.depth(p)==0 and primary==None:
                primary=p
                if 'PolyplopiaPrimary' in frame:
                    frame['PolyplopiaPrimaryCopy'] = frame['PolyplopiaPrimary']
                    del frame['PolyplopiaPrimary']
                frame["PolyplopiaPrimary"]=dataclasses.I3Particle(p)
                frame["PrimaryMass"]=dataclasses.I3Double(mass_dict[str(p.pdg_encoding)])
            #if p.type in NuTypes:
            if (p.type in NuTypes) and (p.energy > 10.):
                neutrinos.append(p)
        if neutrinos==[]:
            return False        
        neutrino=random.choice(neutrinos)
        neutrino_parent=mctree.parent(neutrino)
        frame['ShowerNeutrino'] =dataclasses.I3Particle(neutrino)
        intersection=surface_det.intersection(neutrino.pos, neutrino.dir)#points of intersection
        z_inter=neutrino.pos.z-intersection.first*np.cos(neutrino.dir.zenith)
        #print(z_inter)
        depth = 1948.07 - z_inter
        #print(depth)
        frame['z_inter'] = dataclasses.I3Double(z_inter)
        frame['shower_neutrino_depth'] = dataclasses.I3Double(1948.07 - z_inter)
        frame['shower_neutrino_energy']=dataclasses.I3Double(neutrino.energy)
        frame['shower_neutrino_zenith']=dataclasses.I3Double(neutrino.dir.zenith)
        frame['shower_neutrino_type']=dataclasses.I3Double(neutrino.type)
        frame['NeutrinoParent'] =dataclasses.I3Particle(neutrino_parent)
        mctree.erase_children(neutrino)
                
        

        
tray.Add(NeutrinoSelector)
#required because not interested in NuMu cc muons :/
def issameparticle(particle1, particle2):
        isequal = True
        if (particle1.pos-particle2.pos).magnitude > 1e-3:
            isequal = False
        if np.abs(particle1.dir.zenith-particle2.dir.zenith) > 1e-3:
            isequal = False
        if np.abs(particle1.dir.azimuth-particle2.dir.azimuth) > 1e-3:
            isequal = False
        if np.abs(particle1.energy-particle2.energy) > 1e-3:
            isequal = False
        if np.abs(particle1.time-particle2.time) > 1e-3:
            isequal = False
        if particle1.type != particle2.type:
            isequal = False        
        return isequal
def get_lateral_separation(muon, particle2):

    # c = 2.99792458e8 #* I3Units::m / (I3Units::second)
    gamma_muon=muon.energy/0.10566 #MuonMass=105.66MeV/c^2
    p_muon=np.sqrt((muon.energy)**2-(0.10566)**2)/dataclasses.I3Constants.c
    theta=muon.dir.angle(particle2.dir)
    pt_muon=p_muon*np.sin(theta)
    muon_path_length=muon.pos.z/np.cos(muon.dir.zenith)
    lateral_separation=pt_muon*dataclasses.I3Constants.c*muon_path_length/muon.energy
    
    return lateral_separation


def get_uncorr_deposit_energy(frame):
    if 'I3MCTree_copy' in frame:
        mctree=frame['I3MCTree_copy']
        #mctree=frame['I3MCTree_preMuonProp']
        muon_multiplicity=0
#         muon_multiplicity_10=0
#         muon_multiplicity_100=0
#         muon_multiplicity_200=0
        e_muon_total=0
        muon_list=[]
        muon_energy_list=[]
        for track in MuonGun.Track.harvest(mctree, frame['MMCTrackList']):
            if (boundary_check(track,gcdFile)): continue
            #muon_parent=mctree.parent(track)
            #if issameparticle(muon_parent,frame['NeutrinoParent']): continue
            track_energy_at_det=track.get_energy(surface_det.intersection(track.pos, track.dir).first)
            if track_energy_at_det>10.: ##10 GeV minimum energy
                    e_muon_total+=track_energy_at_det
                    muon_list.append(track)
                    muon_energy_list.append(track_energy_at_det)
                    muon_multiplicity+=1
#             if track_energy_at_det>10:
#                     # print(muon.energy,muon.pos)
#                     muon_multiplicity_10+=1
#             if track_energy_at_det>100:
#                     muon_multiplicity_100+=1
#             if track_energy_at_det>200:
#                     muon_multiplicity_200+=1
        
        frame['Total_UnCorr_Muon_energy'] =dataclasses.I3Double(e_muon_total)
        frame['MuonMultiplicity_UnCorr'] =dataclasses.I3Double(muon_multiplicity)
#         frame['MuonMultiplicity_UnCorr_10'] =dataclasses.I3Double(muon_multiplicity_10)
#         frame['MuonMultiplicity_UnCorr_100'] =dataclasses.I3Double(muon_multiplicity_100)
#         frame['MuonMultiplicity_UnCorr_200'] =dataclasses.I3Double(muon_multiplicity_200)
        
        sorted_energy_particles=sorted(zip(muon_energy_list, muon_list),reverse=True)[:5]
        ctr=0
        for energy,muon in sorted_energy_particles:
            ctr+=1
            energy_name="Muon_UnCorr_Energy_L{}".format(ctr)
            muon_name="Muon_UnCorr_L{}".format(ctr)
            muon_depth_name="Muon_UnCorr_L{}_Depth".format(ctr)
            frame[energy_name]=dataclasses.I3Double(energy)
            frame[muon_name]=dataclasses.I3Particle(muon)
            
            intersection=surface_det.intersection(muon.pos, muon.dir)#points of intersection
            z_inter=muon.pos.z-intersection.first*np.cos(muon.dir.zenith)
            depth=1948.07-z_inter
            frame[muon_depth_name]=dataclasses.I3Double(depth)
            muon_sep_neut_name="Muon_UnCorr_L{}_Neutrino_separation".format(ctr)
            muon_sep_shower_name="Muon_UnCorr_L{}_Shower_separation".format(ctr)
   
            sep1=get_lateral_separation(muon,frame["ShowerNeutrino"])
            frame[muon_sep_neut_name]=dataclasses.I3Double(sep1)
            sep2=get_lateral_separation(muon,frame["PolyplopiaPrimary"])
            frame[muon_sep_shower_name]=dataclasses.I3Double(sep2)
            
tray.Add(get_uncorr_deposit_energy)

import simweights

def corsika_weighter(frame):
    weight_keys = [
    "CylinderLength",
    "CylinderRadius",
    "EnergyPrimaryMax",
    "EnergyPrimaryMin",
    "NEvents",
    "OverSampling",
    "ParticleType",
    "PrimaryEnergy",
    "PrimarySpectralIndex",
    "PrimaryType",
    "ThetaMax",
    "ThetaMin",
    "Weight",
    ]

    particle_keys = ["type", "energy", "zenith"]

    CorsikaWeightMap: dict = {k: [] for k in weight_keys}
    PolyplopiaPrimary: dict = {k: [] for k in ["type", "energy", "zenith"]}
    MCtype_corsika = np.array([])
    MCenergy_corsika = np.array([])

    
    if "FilterMask" in frame:
        # Frame may contain coincident events so select injected primary shower 'PolyplopiaPrimary'
        MCtype_corsika = np.append(MCtype_corsika, frame["PolyplopiaPrimary"].type)
        MCenergy_corsika = np.append(MCenergy_corsika, frame["PolyplopiaPrimary"].energy)

        for k in weight_keys:
            CorsikaWeightMap[k].append(frame["CorsikaWeightMap"][k])
        #print(frame["CorsikaWeightMap"])
        PolyplopiaPrimary["zenith"].append(frame["PolyplopiaPrimary"].dir.zenith)
        PolyplopiaPrimary["type"].append(frame["PolyplopiaPrimary"].type)
        PolyplopiaPrimary["energy"].append(frame["PolyplopiaPrimary"].energy)

    fobj = {"CorsikaWeightMap": CorsikaWeightMap, "PolyplopiaPrimary": PolyplopiaPrimary}
    wobj = simweights.CorsikaWeighter(fobj, nfiles=corsika_filelist_length)
    Weights_GaisserH4a = wobj.get_weights(simweights.GaisserH4a())
    frame["GaisserH4a_weight"] = dataclasses.I3VectorDouble(Weights_GaisserH4a)
    
tray.AddModule(corsika_weighter,'corsika_weighter')

###########book keys for hdf file###########
table_keys=[
            'GaisserH4a_weight',
            "HomogenizedQTot",
            "OneWeight",
            "NEvents",
            "PrimaryMass",
            "PolyplopiaPrimary",
            "I3MCWeightDict",
            "CorsikaWeightMap",
            'PrimaryFlavour',
            'Interaction_Type',
            'OneWeight',
            'MuonWeight',
            'MuonWeight2',
            "MCPrimaryType",
            "MCPrimary",
            'IsUpgoingMuon_L4',
            'Muon_UnCorr_Energy_L1',
            'Muon_UnCorr_Energy_L2',
            'Muon_UnCorr_Energy_L3',
            'Muon_UnCorr_Energy_L4',
            'Muon_UnCorr_Energy_L5',
            'Muon_UnCorr_L1',
            'Muon_UnCorr_L2',
            'Muon_UnCorr_L3',
            'Muon_UnCorr_L4',
            'Muon_UnCorr_L5',
            'Muon_UnCorr_L1_Depth',
            'Muon_UnCorr_L2_Depth',
            'Muon_UnCorr_L3_Depth',
            'Muon_UnCorr_L4_Depth',
            'Muon_UnCorr_L5_Depth',
            'Muon_UnCorr_L1_Neutrino_separation',
            'Muon_UnCorr_L2_Neutrino_separation',
            'Muon_UnCorr_L3_Neutrino_separation',
            'Muon_UnCorr_L4_Neutrino_separation',
            'Muon_UnCorr_L5_Neutrino_separation',
            'Muon_UnCorr_L1_Shower_separation',
            'Muon_UnCorr_L2_Shower_separation',
            'Muon_UnCorr_L3_Shower_separation',
            'Muon_UnCorr_L4_Shower_separation',
            'Muon_UnCorr_L5_Shower_separation',
            'MuonMultiplicity_UnCorr',
            
            'ShowerNeutrino',
            'Total_UnCorr_Muon_energy',
            'BrightestMuonAtDepth',
            'Energy_BrightestMuonAtDepth',
            'SPEFit4_rlogL' ,
            'SPEFitSingle_rlogL',
            'SPEFit2_rlogL',
            "I3EventHeader",
            "NEvents",
            'PolyplopiaPrimaryCopy',
            'z_inter',
            'shower_neutrino_depth',
            'shower_neutrino_energy',
            'shower_neutrino_zenith',
            'shower_neutrino_type'
           ]

tray.Add(I3HDFWriter,'HDFwriter',
              Output=outfile+".hdf5",
              keys         = table_keys,
              SubEventStreams = ['InIceSplit']
)
# tray.AddModule('I3Writer', 'writer',
#     DropOrphanStreams=[icetray.I3Frame.DAQ],
#     Streams=[ icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
#     filename=outfile+".i3.zst")

tray.AddModule('TrashCan', 'thecan')

tray.Execute()
tray.Finish()

del tray

stop_time = time.asctime()

print ('Started:', start_time)
print ('Ended:', stop_time)

