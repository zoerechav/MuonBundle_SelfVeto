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

parser.add_option("-g", "--gcd", action="store", type="string", default="/cvmfs/icecube.opensciencegrid.org/users/vbasu/GCD_files/Level2pass2_IC86.2016_data_Run00129060_0117_89_292_GCD.i3.zst", dest="GCD", help="GCD file for input i3 file")
parser.add_option("-i", "--input", action="store", type="string", default="/data/sim/IceCube/2016/filtered/level2/CORSIKA-in-ice/20904/0000000-0000999/Level2_IC86.2016_corsika.020904.000000.i3.zst", dest="INPUT", help="Output i3 file")
# parser.add_option("-o", "--output", action="store", type="string", default="/home/vbasu/testoutput/L2_selfveto", dest="OUTPUT", help="Output i3 file")

#parser.add_option("-i", "--input", action="store", type="string", default="/data/sim/IceCube/2011/filtered/level2/CORSIKA-in-ice/10661/00000-00999/Level2_IC86.2011_corsika.010661.000069.00.i3.bz2", dest="INPUT", help="Output i3 file")
parser.add_option("-o", "--output", action="store", type="string", default="/home/vbasu/testoutput/CorsikaShowerReader2", dest="OUTPUT", help="Output i3 file")

parser.add_option("-f","--grid", action="store", type="int", default=0, dest="GridFTP", help="to use gridftp or not. 0 is false and 1 is true")

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()

gridftp=bool(options.GridFTP)
if gridftp==True:
    grid = 'gsiftp://gridftp-users.icecube.wisc.edu'
    outfile = grid+options.OUTPUT
    gcdFile=options.GCD
    infiles=[options.GCD, grid+options.INPUT]
    
else:
    infile = options.INPUT
    outfile = options.OUTPUT
    gcdFile=options.GCD
        
tray = I3Tray()
if gridftp==True:
    tray.context['I3FileStager'] = dataio.get_stagers()

    tray.AddModule('I3Reader', 'reader', FilenameList=infiles)
        
    print ("Reading input file...",  grid+options.INPUT)
else:
    tray.AddSegment(dataio.I3Reader, 'reader', FilenameList=[options.GCD, infile])
    print ("Reading input file...",  options.INPUT)    


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
    bound_2D=[]
    surface_det = MuonGun.ExtrudedPolygon.from_file(gcdFile, padding=0)##Build Polygon from I3Geometry
    
    x=[(surface_det.x[i],surface_det.y[i])for i in range(len(surface_det.x))]###getting only x and y
    bound_2D=mpltPath.Path(x)#Projection of detector on x,y plane
    return bound_2D, surface_det

bound_2D,surface_det = get_surface_det(gcdFile=gcdFile)
def boundary_check(particle1,gcdFile=None):
    ####checks if particle is inside the detector###
   
    inlimit = False  
    if ((particle1.pos.z <=max(surface_det.z)) and (particle1.pos.z>=min(surface_det.z))):
        if bound_2D.contains_points([(particle1.pos.x, particle1.pos.y)]):
            inlimit=True

    return inlimit

def NeutrinoSelector(frame):
    if 'I3MCTree' in frame and 'LineFit' in frame:
        frame['I3MCTree_copy']=frame['I3MCTree']
        mctree = frame['I3MCTree_copy']
        e_dep_total=0
        bound_2D,surface_det = get_surface_det(gcdFile=gcdFile)
        neutrino = None
        neutrinos=[]
        primary=None
        for p in mctree:#search for neutrino
            if mctree.depth(p)==0 and primary==None:
                primary=p
                frame["PolyplopiaPrimary"]=dataclasses.I3Particle(p)
                frame["PrimaryMass"]=dataclasses.I3Double(mass_dict[str(p.pdg_encoding)])
            if p.type in NuTypes:
                neutrinos.append(p)
                neutrino=p
                neutrino_parent=mctree.parent(p)
                frame['ShowerNeutrino'] =dataclasses.I3Particle(p)
                intersection=surface_det.intersection(p.pos, p.dir)#points of intersection
                z_inter=p.pos.z-intersection.first*np.cos(p.dir.zenith)
                depth=1948.07-z_inter
                frame["PrimaryDepthMC"]=dataclasses.I3Double(depth)
                frame['NeutrinoParent'] =dataclasses.I3Particle(neutrino_parent)
                mctree.erase_children(p)
                break
        if neutrino==None:
            return False

        
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

def get_deposit_energy(frame):
    energy_max=0
    e_max_muon=None
    losses_max=0
    
    brightest_muon=None
    bound_2D,surface_det = get_surface_det(gcdFile=gcdFile)
    if 'I3MCTree_copy' in frame:
                            
        mctree = frame['I3MCTree_copy'] 
        for muon in mctree:#get muon energies
            if (boundary_check(muon)): continue
            if muon.type in MuTypes:
                losses = 0
                daughters = mctree.get_daughters(muon)
                
                if len(daughters)==0:continue
                
                for p in daughters:

                    if not p.is_cascade: continue
                    if not p.location_type == dataclasses.I3Particle.InIce: continue
                    if p.shape == p.Dark: continue
                    if p.type in [p.Hadrons, p.PiPlus, p.PiMinus, p.NuclInt]:
                        #hadlosses += p.energy
                        if p.energy < 1*I3Units.GeV:
                            losses += 0.8*p.energy
                        else:
                            energyScalingFactor = 1.0 + ((p.energy/I3Units.GeV/0.399)**-0.130)*(0.467 - 1)
                            losses += energyScalingFactor*p.energy
                    else:
                        #emlosses += p.energy
                        losses += p.energy
                
                if losses>losses_max:
                    losses_max=losses
                    brightest_muon=muon
        if brightest_muon == None: return False
        frame['BrightestMuon'] =dataclasses.I3Particle(brightest_muon)
        intersection=surface_det.intersection(brightest_muon.pos, brightest_muon.dir)#points of intersection
        z_inter=brightest_muon.pos.z-intersection.first*np.cos(brightest_muon.dir.zenith)
        depth=1948.07-z_inter
        frame["BrightestMuonDepth"]=dataclasses.I3Double(depth)
        e0=0
        for track in MuonGun.Track.harvest(mctree, frame['MMCTrackList']):
            if (boundary_check(track)): continue
            intersections = surface_det.intersection(track.pos, track.dir)
            # Get the corresponding energies
            e0, e1 = track.get_energy(intersections.first), track.get_energy(intersections.second)            
            if e0>energy_max:
                energy_max=e0
                if 'HEMuon' in frame:
                    del frame['HEMuon']
                    del frame['Energy_HEMuonAtDepth']
                    del frame['HEMuonDepth']
                frame['HEMuon'] =dataclasses.I3Particle(track)
                frame['Energy_HEMuonAtDepth']=dataclasses.I3Double(e0)
                intersection=surface_det.intersection(frame['HEMuon'].pos, frame['HEMuon'].dir)#points of intersection
                z_inter=frame['HEMuon'].pos.z-intersection.first*np.cos(frame['HEMuon'].dir.zenith)
                depth=1948.07-z_inter
                frame["HEMuonDepth"]=dataclasses.I3Double(depth)
            if not issameparticle(track, brightest_muon):continue
            intersections = surface_det.intersection(track.pos, track.dir)
            # Get the corresponding energies
            e0, e1 = track.get_energy(intersections.first), track.get_energy(intersections.second)
            frame['Energy_BrightestMuonAtDepth']=dataclasses.I3Double(e0)
            frame['BrightestMuonAtDepth']=dataclasses.I3Particle(track)
            # break
        muon_multiplicity=0
        muon_multiplicity_10=0
        muon_multiplicity_100=0
        muon_multiplicity_200=0
        e_muon_total=0
        muon_list=[]
        muon_energy_list=[]
        for track in MuonGun.Track.harvest(mctree, frame['MMCTrackList']):
            if (boundary_check(track)): continue
            track_energy_at_det=track.get_energy(surface_det.intersection(track.pos, track.dir).first)
            if track_energy_at_det>0:
                    e_muon_total+=track_energy_at_det
                    muon_list.append(track)
                    muon_energy_list.append(track_energy_at_det)
                    muon_multiplicity+=1
            if track_energy_at_det>10:
                    # print(muon.energy,muon.pos)
                    muon_multiplicity_10+=1
            if track_energy_at_det>100:
                    muon_multiplicity_100+=1
            if track_energy_at_det>200:
                    muon_multiplicity_200+=1
        frame['max_depo_energy'] =dataclasses.I3Double(losses_max)
        frame['Total_Muon_energy'] =dataclasses.I3Double(e_muon_total)
        frame['MuonMultiplicity'] =dataclasses.I3Double(muon_multiplicity)
        frame['MuonMultiplicity_10'] =dataclasses.I3Double(muon_multiplicity_10)
        frame['MuonMultiplicity_100'] =dataclasses.I3Double(muon_multiplicity_100)
        frame['MuonMultiplicity_200'] =dataclasses.I3Double(muon_multiplicity_200)
        

        sorted_energy_particles=sorted(zip(muon_energy_list, muon_list),reverse=True)[:5]
        ctr=0
        for energy,muon in sorted_energy_particles:
            ctr+=1
            energy_name="Muon_Energy_L{}".format(ctr)
            muon_name="Muon_L{}".format(ctr)
            muon_depth_name="Muon_L{}_Depth".format(ctr)
            
            frame[energy_name]=dataclasses.I3Double(energy)
            frame[muon_name]=dataclasses.I3Particle(muon)
            
            intersection=surface_det.intersection(muon.pos, muon.dir)#points of intersection
            z_inter=muon.pos.z-intersection.first*np.cos(muon.dir.zenith)
            depth=1948.07-z_inter
            muon_sep_neut_name="Muon_L{}_Neutrino_separation".format(ctr)
            muon_sep_shower_name="Muon_L{}_Shower_separation".format(ctr)
            frame[muon_depth_name]=dataclasses.I3Double(depth)
            sep1=get_lateral_separation(muon,frame["ShowerNeutrino"])
            frame[muon_sep_neut_name]=dataclasses.I3Double(sep1)
            sep2=get_lateral_separation(muon,frame["PolyplopiaPrimary"])
            frame[muon_sep_shower_name]=dataclasses.I3Double(sep2)

tray.Add(get_deposit_energy)

def get_corr_deposit_energy(frame):
    bound_2D,surface_det = get_surface_det(gcdFile=gcdFile)
    if 'I3MCTree_copy' in frame:
        mctree=frame['I3MCTree_copy']
        muon_multiplicity=0
        muon_multiplicity_10=0
        muon_multiplicity_100=0
        muon_multiplicity_200=0
        e_muon_total=0
        muon_list=[]
        muon_energy_list=[]
        for track in MuonGun.Track.harvest(mctree, frame['MMCTrackList']):
            if (boundary_check(track)): continue
            muon_parent=mctree.parent(track)
            if not issameparticle(muon_parent,frame['NeutrinoParent']): continue
            track_energy_at_det=track.get_energy(surface_det.intersection(track.pos, track.dir).first)
            if track_energy_at_det>0:
                    e_muon_total+=track_energy_at_det
                    muon_list.append(track)
                    muon_energy_list.append(track_energy_at_det)
                    muon_multiplicity+=1
            if track_energy_at_det>10:
                    # print(muon.energy,muon.pos)
                    muon_multiplicity_10+=1
            if track_energy_at_det>100:
                    muon_multiplicity_100+=1
            if track_energy_at_det>200:
                    muon_multiplicity_200+=1
        
        frame['Total_Corr_Muon_energy'] =dataclasses.I3Double(e_muon_total)
        frame['MuonMultiplicity_Corr'] =dataclasses.I3Double(muon_multiplicity)
        frame['MuonMultiplicity_Corr_10'] =dataclasses.I3Double(muon_multiplicity_10)
        frame['MuonMultiplicity_Corr_100'] =dataclasses.I3Double(muon_multiplicity_100)
        frame['MuonMultiplicity_Corr_200'] =dataclasses.I3Double(muon_multiplicity_200)
        
        sorted_energy_particles=sorted(zip(muon_energy_list, muon_list),reverse=True)[:5]
        ctr=0
        for energy,muon in sorted_energy_particles:
            ctr+=1
            energy_name="Muon_Corr_Energy_L{}".format(ctr)
            muon_name="Muon_Corr_L{}".format(ctr)
            muon_depth_name="Muon_Corr_L{}_Depth".format(ctr)
            frame[energy_name]=dataclasses.I3Double(energy)
            frame[muon_name]=dataclasses.I3Particle(muon)
            
            intersection=surface_det.intersection(muon.pos, muon.dir)#points of intersection
            z_inter=muon.pos.z-intersection.first*np.cos(muon.dir.zenith)
            depth=1948.07-z_inter
            frame[muon_depth_name]=dataclasses.I3Double(depth)
            muon_sep_neut_name="Muon_Corr_L{}_Neutrino_separation".format(ctr)
            muon_sep_shower_name="Muon_Corr_L{}_Shower_separation".format(ctr)
            
            sep1=get_lateral_separation(muon,frame["ShowerNeutrino"])
            frame[muon_sep_neut_name]=dataclasses.I3Double(sep1)
            sep2=get_lateral_separation(muon,frame["PolyplopiaPrimary"])
            frame[muon_sep_shower_name]=dataclasses.I3Double(sep2)
            
tray.Add(get_corr_deposit_energy)
def get_uncorr_deposit_energy(frame):
    bound_2D,surface_det = get_surface_det(gcdFile=gcdFile)
    if 'I3MCTree_copy' in frame:
        mctree=frame['I3MCTree_copy']
        muon_multiplicity=0
        muon_multiplicity_10=0
        muon_multiplicity_100=0
        muon_multiplicity_200=0
        e_muon_total=0
        muon_list=[]
        muon_energy_list=[]
        for track in MuonGun.Track.harvest(mctree, frame['MMCTrackList']):
            if (boundary_check(track)): continue
            muon_parent=mctree.parent(track)
            if issameparticle(muon_parent,frame['NeutrinoParent']): continue
            track_energy_at_det=track.get_energy(surface_det.intersection(track.pos, track.dir).first)
            if track_energy_at_det>0:
                    e_muon_total+=track_energy_at_det
                    muon_list.append(track)
                    muon_energy_list.append(track_energy_at_det)
                    muon_multiplicity+=1
            if track_energy_at_det>10:
                    # print(muon.energy,muon.pos)
                    muon_multiplicity_10+=1
            if track_energy_at_det>100:
                    muon_multiplicity_100+=1
            if track_energy_at_det>200:
                    muon_multiplicity_200+=1
        
        frame['Total_UnCorr_Muon_energy'] =dataclasses.I3Double(e_muon_total)
        frame['MuonMultiplicity_UnCorr'] =dataclasses.I3Double(muon_multiplicity)
        frame['MuonMultiplicity_UnCorr_10'] =dataclasses.I3Double(muon_multiplicity_10)
        frame['MuonMultiplicity_UnCorr_100'] =dataclasses.I3Double(muon_multiplicity_100)
        frame['MuonMultiplicity_UnCorr_200'] =dataclasses.I3Double(muon_multiplicity_200)
        
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
from icecube.weighting.fluxes import GaisserH4a, Hoerandel5,Honda2004
from icecube.weighting import weighting, get_weighted_primary   
#weight Corsika sample now!
corsika_dids = [10661,20904] 
simprod_nfiles={"21002":9979,"21217":21851,"21218":11991,"21219":15376,
                "21220":9999,"21221":10000,
                # "20904":740335,
                "10661":10000, 
                "20904":50000, 
                "20891": 497612, "20881": 99904,
                "20852":99094, "20849":9948, "20848":99782, "20789":99998,
                "20788":99998, "20787":99743,
                "scat5":9992,"scat-5":9994,"abs5":9993,"abs-5":9988,
                "domeff90":9831,"domeff95":9829,"domeff105":9813,"domeff110":9813,
                "p0=-2.0_p1=-0.2":9822,"p0=-2.0_p1=0.0":9825,"p0=-2.0_p1=0.2":9828,
                "p0=-1.0_p1=-0.2":9840,"p0=-1.0_p1=0.0":9839,"p0=-1.0_p1=0.2":9830,
                "p0=0.0_p1=-0.2":9830,"p0=0.0_p1=0.0":9834,"p0=0.0_p1=0.2":9847,
                "p0=1.0_p1=-0.2":9831,"p0=1.0_p1=0.0":9832,"p0=1.0_p1=0.2":9817
                }
def weighter_corsika(frame,flux=None):
    if not frame.Has("MCPrimary"):
        get_weighted_primary(frame, MCPrimary="MCPrimary")
    MCPrimary = frame["MCPrimary"]
    energy = MCPrimary.energy
    
    ptype = MCPrimary.type
    placeholder = True
    for DID in corsika_dids:
        if(placeholder): 
            generator = weighting.from_simprod(DID) * simprod_nfiles["%5.0f"%DID];
            placeholder = False
        else: 
            generator += weighting.from_simprod(DID) * simprod_nfiles["%5.0f"%DID];
    if flux is 'GaisserH4a':
        flux=GaisserH4a()
        weights = flux(energy, ptype) / generator(energy, ptype)
        if np.isnan(weights): weights=0.0;
        frame['Weight_GaisserH4a']=dataclasses.I3Double(weights)
    elif flux is 'Hoerandel5':
        flux=Hoerandel5()
        weights = flux(energy, ptype) / generator(energy, ptype)
        if np.isnan(weights): weights=0.0;
        frame['Weight_Hoerandel5']=dataclasses.I3Double(weights)
    elif flux is 'Honda2004':
        flux=Honda2004()
        weights = flux(energy, ptype) / generator(energy, ptype)
        if np.isnan(weights): weights=0.0;
        frame['Weight_Honda2004']=dataclasses.I3Double(weights)
    #do a little cleanup now
    flux = 0; generator =0; weights = 0; MCPrimary = 0;
tray.AddModule(weighter_corsika, "WeighterCorsika_Gaisser",flux='GaisserH4a')
tray.AddModule(weighter_corsika, "WeighterCorsika_Hoerandel",flux='Hoerandel5')
tray.AddModule(weighter_corsika, "WeighterCorsika_Honda",flux='Honda2004')
###########book keys for hdf file###########
table_keys=["IsCascade_true","IsTrack_true",
            "HomogenizedQTot","HomogenizedQTot_toposplit",'HomogenizedQTot_TWTS',
            "OneWeight",'EntryMuon','Entry_Energies','Entry_Zenith','Entry_Azimuth',
            "NEvents","PrimaryMass",
            "VetoLayerQTot","VetoLayer0", "VetoLayer1","IsHESE_ck","IsHESE_jvs","IsMESEL3",'EnteringMuon_1',
            "PolyplopiaPrimary", "I3MCWeightDict","CorsikaWeightMap","I3PrimaryInjectorInfo",
            'VisibleEnergyMC','ReconstructedEnergyMonopod_L5',
                'ReconstructedDirectionMonopod_L5','ReconstructedTypeMonopod_L5',
                'MuEXAngularEnergy','MuEXAngularDirection',
                'TrackFitEnergy','TrackFitDirection',
                'TrackFit_AvgDistQ','L5MonopodFit4_AvgDistQ','PrimaryFlavour','PrimaryDepthMC','Interaction_Type',
                'MillipedeDepositedEnergy_TWTS',
                'OneWeight','MuonWeight','MuonWeight2',
                'Weight_Honda2004','Weight_Hoerandel5','Weight_GaisserH4a',
                'PrimaryEnergyMC','PrimaryZenithMC',"MCPrimaryType","MCPrimary",
                'IsUpgoingMuon_L4',
                'Muon_Energy_L1','Muon_Energy_L2','Muon_Energy_L3','Muon_Energy_L4','Muon_Energy_L5',
                'Muon_L1','Muon_L2','Muon_L3','Muon_L4','Muon_L5',
                'Muon_L1_Neutrino_separation','Muon_L2_Neutrino_separation','Muon_L3_Neutrino_separation','Muon_L4_Neutrino_separation','Muon_L5_Neutrino_separation',
                'Muon_L1_Shower_separation','Muon_L2_Shower_separation','Muon_L3_Shower_separation','Muon_L4_Shower_separation','Muon_L5_Shower_separation',
                'Muon_L1_Depth','Muon_L2_Depth','Muon_L3_Depth','Muon_L4_Depth','Muon_L5_Depth',
                'MuonEnergy','MuonMultiplicity','MuonMultiplicity_10','MuonMultiplicity_100','MuonMultiplicity_200',
                'Muon_Corr_Energy_L1','Muon_Corr_Energy_L2','Muon_Corr_Energy_L3','Muon_Corr_Energy_L4','Muon_Corr_Energy_L5',
                'Muon_Corr_L1','Muon_Corr_L2','Muon_Corr_L3','Muon_Corr_L4','Muon_Corr_L5',
                'Muon_Corr_L1_Depth','Muon_Corr_L2_Depth','Muon_Corr_L3_Depth','Muon_Corr_L4_Depth','Muon_Corr_L5_Depth',
                'Muon_Corr_L1_Neutrino_separation','Muon_Corr_L2_Neutrino_separation','Muon_Corr_L3_Neutrino_separation','Muon_Corr_L4_Neutrino_separation','Muon_Corr_L5_Neutrino_separation',
                'Muon_Corr_L1_Shower_separation','Muon_Corr_L2_Shower_separation','Muon_Corr_L3_Shower_separation','Muon_Corr_L4_Shower_separation','Muon_Corr_L5_Shower_separation',
                'MuonMultiplicity_Corr','MuonMultiplicity_Corr_10','MuonMultiplicity_Corr_100','MuonMultiplicity_Corr_200',

                'Muon_UnCorr_Energy_L1','Muon_UnCorr_Energy_L2','Muon_UnCorr_Energy_L3','Muon_UnCorr_Energy_L4','Muon_UnCorr_Energy_L5',
                'Muon_UnCorr_L1','Muon_UnCorr_L2','Muon_UnCorr_L3','Muon_UnCorr_L4','Muon_UnCorr_L5',
                'Muon_UnCorr_L1_Depth','Muon_UnCorr_L2_Depth','Muon_UnCorr_L3_Depth','Muon_UnCorr_L4_Depth','Muon_UnCorr_L5_Depth',
                'Muon_UnCorr_L1_Neutrino_separation','Muon_UnCorr_L2_Neutrino_separation','Muon_UnCorr_L3_Neutrino_separation','Muon_UnCorr_L4_Neutrino_separation','Muon_UnCorr_L5_Neutrino_separation',
                'Muon_UnCorr_L1_Shower_separation','Muon_UnCorr_L2_Shower_separation','Muon_UnCorr_L3_Shower_separation','Muon_UnCorr_L4_Shower_separation','Muon_UnCorr_L5_Shower_separation',
                'MuonMultiplicity_UnCorr','MuonMultiplicity_UnCorr_10','MuonMultiplicity_UnCorr_100','MuonMultiplicity_UnCorr_200',
                
                'ShowerNeutrino','max_depo_energy','BrightestMuon','HEMuon','BrightestMuonDepth','Energy_HEMuonAtDepth','HEMuonDepth',
                'IsHese','IsHESE_jvs',"IsHESE_ck",'IsMESEL3',
                'IsCascade_dnn','IsTrack_dnn','Total_Muon_energy',
                'BrightestMuonAtDepth','Energy_BrightestMuonAtDepth',
                'IsStartingEvent_L4','UpgoingMuon_CoincCut',
                'IsCascade_recoL4','IsTrack_recoL4','SPEFit4_rlogL' ,'SPEFitSingle_rlogL','SPEFit2_rlogL',
                'IsVetoTrack',"I3EventHeader",'SPEFit4_offlineFitParams' ,'SPEFitSingle_offlineFitParams','SPEFit2_offlineFitParams',
                'IsCascade_true','IsTrack_true',
                "NEvents"]

tray.Add(I3HDFWriter,'HDFwriter',
              Output=outfile+".hdf5",
              keys         = table_keys,
              SubEventStreams = ['I3NullSplitter_0000','topological_split','in_ice','InIceSplit','NullSplit']
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
