#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# usage: basf2 MakeNtuple_multi.py "./20210402/evt-0.mdst"
# last: 2021-07-19-1

import os
import sys

import basf2
import modularAnalysis as ma
import variables as va
import variables.collections as vc
import variables.utils as vu
from vertex import kFit
from b2biiConversion import convertBelleMdstToBelleIIMdst
from b2biiMonitors import addBeamParamsConversionMonitors
from b2biiMonitors import addTrackConversionMonitors
from b2biiMonitors import addNeutralsConversionMonitors

import fei
from glob import glob

os.environ['USE_GRAND_REPROCESS_DATA'] = '1'
os.environ['PGUSER'] = 'g0db'

debug = False

basf2.set_debug_level(10)

#basf2.conditions.prepend_globaltag("analysis_tools_light-2104-poseidon data_reprocessing_proc9 online")
basf2.conditions.prepend_globaltag(ma.getAnalysisGlobaltag())

# set random seed
basf2.set_random_seed(42)

# If you want to create the monitoring histograms (recommended in the beginning of your analysis), just provide any argument
monitoring = False
if len(sys.argv) == 1:
    basf2.B2WARNING("In the beginning of an analysis it is recommended to study the monitoring histograms.\n"
                    "These tell you if the conversion works as expected.\n"
                    "If you want to create them, just provide any argument to this script.")
else:
    monitoring = False

fName = sys.argv[1]
destination = sys.argv[2]

if not fName.endswith(".mdst"): sys.exit(1)
else:
    basename = os.path.basename(fName)
    name = os.path.splitext(basename)[0]
print("=================="+fName+" is conducted=================")


# Convert
my_path = basf2.create_path()
inputfile=fName

convertBelleMdstToBelleIIMdst(inputfile, applySkim=True, enableLocalDB=False, path=my_path)
#ma.inputMdst(environmentType='default',filename=inputfile,path=my_path)

# Reconstruct
if monitoring:
    # Create monitoring histograms if requested
    addBeamParamsConversionMonitors(path=my_path)
    addTrackConversionMonitors(path=my_path)
    addNeutralsConversionMonitors(path=my_path)

# fill particle list
stdPhotons.loadStdGoodBellePhoton(path=my_path)

ma.fillParticleList(decayString="K+:mychargedKaon", cut="atcPIDBelle(3,2)>0.6 and eIDBelle < 0.9 and muIDBelle < 0.9 and dr < 2 and dz < 4",path=my_path)
ma.fillParticleList(decayString="K+:all", cut="", path=my_path)
ma.cutAndCopyList('K_S0:good', 'K_S0:mdst', cut='goodBelleKshort', path=my_path)
ma.cutAndCopyList("pi0:good", "pi0:mdst", cut="M>0.118 and M<0.15", path=my_path)
#ma.cutAndCopyList("gamma:good", "gamma:mdst", cut="E>0.05",path=my_path)
#ma.fillParticleList(decayString="pi+:all", cut="", path=my_path)
#ma.fillParticleList(decayString="gamma:all", cut="", path=my_path)

#ma.reconstructDecay("pi0:all -> gamma:all gamma:all", cut="", path = my_path)

#ma.reconstructDecay("K_S0:ch0 -> pi+:all pi-:all", cut ="", dmID=0, path=my_path)
#ma.reconstructDecay("K_S0:ch1 -> pi0:all pi0:all", cut ="", dmID=1, path=my_path)
#ma.copyLists(outputListName="K_S0:myKaonShort", inputListNames=["K_S0:ch0", "K_S0:ch1"], path=my_path)

# ==============================
# FEI
# ==============================

particles = fei.get_default_channels()
configuration = fei.config.FeiConfiguration(prefix="FEI_B2BII_light-2012-minos",training=False, monitor=False)
feistate = fei.get_path(particles, configuration)
my_path.add_path(feistate.path)

# --- Debugger ---
if debug:
    ma.printList('B+:generic',full=False,path=path)
    printVariableValues('B+:generic',['isSignal','mcPDG','genMotherID','genMotherPDG'],path=path)

# ma.applyCuts("B+:generic","Mbc>5.27 and abs(deltaE)<0.1", path=my_path)

# tag side
track_selection = "dr < 2 and dz < 4"
cluster_selection = '[[E > 0.10 and formula(theta/3.14*180) < 31.4] or \
                      [E > 0.05 and formula(theta/3.14*180) > 31.4 and formula(theta/3.14*180) < 130.7] or \
                      [E > 0.15 and formula(theta/3.14*180) > 130.7]]'
cleanMask = ("cleanMask",track_selection, cluster_selection)
ma.buildRestOfEvent("B+:generic",path=my_path)
ma.appendROEMasks("B+:generic",[cleanMask],path=my_path)
ma.buildContinuumSuppression(list_name="B+:generic",roe_mask = "cleanMask", path=my_path)
# To Do: background suppresion using ROE

# signal side
# B+ -> K+ nu nu
ma.reconstructDecay("B+:sig_ch0 -> K+:all ?nu", cut="", dmID=0, path = my_path)
# B0 -> K0 nu nu & K0 -> K0S
# K0S -> pi+ pi- or pi0 pi0
#ma.reconstructDecay("B0:neutral_sig_ch0 -> K_S0:myKaonShort ?nu", cut="", dmID=0, path=my_path)

# To Do: add other decay modes
ma.copyLists(outputListName="B+:sig", inputListNames=["B+:sig_ch0"], path=my_path)
#ma.copyLists(outputListName="B0:neutral_sig", inputListNames=["B0:neutral_sig_ch0"],path=my_path)

# Upsilon(4S) -> Btag K
ma.reconstructDecay("Upsilon(4S):withoutneutrino_charged -> B+:generic B-:sig",cut ="", dmID = 0, path=my_path)
#ma.reconstructDecay("Upsilon(4S):withoutneutrino_neutral -> B0:generic anti-B0:neutral_sig",cut ="", dmID = 1, path=my_path)
#ma.copyLists(outputListName="Upsilon(4S):withoutneutrino", inputListNames=["Upsilon(4S):withoutneutrino_charged", "Upsilon(4S):withoutneutrino_neutral"], path=my_path)
ma.copyLists(outputListName="Upsilon(4S):withoutneutrino", inputListNames=["Upsilon(4S):withoutneutrino_charged"], path=my_path)
ma.buildRestOfEvent("Upsilon(4S):withoutneutrino", inputParticlelists=["K_S0:good", "pi0:good", "gamma:mdst"], path=my_path)
    
# apply cut: no charged track, E, cluster on ROE
ma.appendROEMasks("Upsilon(4S):withoutneutrino",[cleanMask],path=my_path)
# ma.applyCuts("Upsilon(4S):withoutneutrino","nROE_Tracks(cleanMask)==0",path=my_path)

# best candidate selection
va.variables.addAlias("SignalProbofBtag_rank","daughter(0,extraInfo(SignalProbability))")
ma.rankByHighest(particleList="Upsilon(4S):withoutneutrino", variable="SignalProbofBtag_rank",allowMultiRank=True,outputVariable="Upsilon_rank",path=my_path)
# ma.applyCuts("Upsilon(4S):withoutneutrino","extraInfo(Upsilon_rank) == 1",path=my_path)
    
# --- MC truth ---
ma.looseMCTruth(list_name = "B+:generic", path = my_path)
ma.looseMCTruth(list_name = "B+:sig", path = my_path)
#ma.looseMCTruth(list_name = "B0:generic", path = my_path)
#ma.looseMCTruth(list_name = "B0:neutral_sig", path = my_path)
ma.looseMCTruth(list_name = "Upsilon(4S):withoutneutrino", path = my_path)

# --- decay string ---
my_path.add_module('ParticleMCDecayString', listName="B+:generic", fileName=destination + "/" + name+'hashmap_Btag.root')
my_path.add_module('ParticleMCDecayString', listName="B+:sig", fileName=destination + "/" + name+'hashmap_Bsig.root')
#my_path.add_module('ParticleMCDecayString', listName="B0:generic", fileName=name+'hashmap_Btag_neutral.root')
#my_path.add_module('ParticleMCDecayString', listName="B0:neutral_sig", fileName=name+'hashmap_Bsig_neutral.root')
my_path.add_module('ParticleMCDecayString', listName='Upsilon(4S):withoutneutrino', fileName=destination + "/" + name+'hashmap_Upsilon.root')

# get variables
Kinematics = ["E", "InvM", "M", "p", "phi", "theta"]
Kinematics_CMS = ["useCMSFrame(E)", "useCMSFrame(p)", "useCMSFrame(phi)", "useCMSFrame(theta)"]
Kinematics_RecoilRestFrame = ["useTagSideRecoilRestFrame(daughter(1,E),0)", "useTagSideRecoilRestFrame(daughter(1,p),0)", "useTagSideRecoilRestFrame(daughter(1,phi),0)", "useTagSideRecoilRestFrame(daughter(1,theta),0)"]
Btag_cut = ["Mbc", "deltaE"]
mcvar = ["mcPDG", "mcE", "mcP", "mcPhi", "mcTheta", "genParticleID", "isSignal", "isSignalAcceptMissingNeutrino", "isSignalAcceptMissingGamma"]
loosemcvar = ["extraInfo(looseMCMotherPDG)", "extraInfo(looseMCWrongDaughterN)", "extraInfo(looseMCWrongDaughterPDG)", "extraInfo(looseMCWrongDaughterBiB)"]
decayhash = ['extraInfo(DecayHash)', 'extraInfo(DecayHashExtended)']
othervar = ["PDG", "extraInfo(decayModeID)"]
continuumsup_vars = ["R2", "thrustBm", "thrustOm", "cosTBTO", "cosTBz", "KSFWVariables(et)", "KSFWVariables(mm2)", "KSFWVariables(hso00)", "KSFWVariables(hso02)", "KSFWVariables(hso04)", "KSFWVariables(hso10)", "KSFWVariables(hso12)", "KSFWVariables(hso14)", "KSFWVariables(hso20)", "KSFWVariables(hso22)", "KSFWVariables(hso24)", "CleoConeCS(1)", "CleoConeCS(2)", "CleoConeCS(3)", "CleoConeCS(4)", "CleoConeCS(5)", "CleoConeCS(6)", "CleoConeCS(7)", "CleoConeCS(8)", "CleoConeCS(9)"]

Btag_vars = vu.create_aliases(list_of_variables = Kinematics + Btag_cut + Kinematics_CMS + mcvar + loosemcvar + decayhash + othervar + continuumsup_vars + ["extraInfo(SignalProbability)"], wrapper = "daughter(0,{variable})",prefix="Btag")

Bsig_vars = vu.create_aliases(list_of_variables = Kinematics + Kinematics_CMS + mcvar + loosemcvar + decayhash + othervar + ["daughter(0, atcPIDBelle(3,2))", "daughter(0, eIDBelle)", "daughter(0, muIDBelle)", "daughter(0, dr)", "daughter(0, dz)", "daughter(0, mcPDG)", "daughter(0, genParticleID)"], wrapper = "daughter(1,{variable})", prefix="Bsig")

U_vars = Kinematics + Kinematics_CMS + Kinematics_RecoilRestFrame + mcvar + decayhash + othervar + ["extraInfo(Upsilon_rank)", "nROE_ECLClusters(cleanMask)", "nROE_NeutralECLClusters(cleanMask)", "nROE_KLMClusters", "roeE(cleanMask)", "roeMC_E", "nROE_Tracks(cleanMask)", "weMissE(cleanMask,0)", "weMissE(cleanMask,5)", "roeEextra(cleanMask)", "roeNeextra(cleanMask)", "useCMSFrame(roeE(cleanMask))", "useCMSFrame(roeEextra(cleanMask))", "useCMSFrame(roeNeextra(cleanMask))", "nROE_ParticlesInList(gamma:mdst)", "nROE_ParticlesInList(K_S0:good)", "nROE_ParticlesInList(pi0:good)"]

output_file = destination + "/" + name+"_Ntuple.root"
ma.variablesToNtuple(decayString="Upsilon(4S):withoutneutrino",variables=Btag_vars,filename=output_file,treename="Btag",path=my_path)
ma.variablesToNtuple(decayString="Upsilon(4S):withoutneutrino",variables=Bsig_vars,filename=output_file,treename="Bsig",path=my_path)
ma.variablesToNtuple(decayString="Upsilon(4S):withoutneutrino",variables=U_vars,filename=output_file,treename="Upsilon",path=my_path)    

# progress
basf2.process(my_path)

# Print call statistics
print(basf2.statistics)
