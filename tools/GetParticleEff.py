#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# usage: basf2 MakeNtuple_multi.py "./20210402/evt-0.mdst"
# last: 2022-07-14-00

import os
import sys

import basf2
import modularAnalysis as ma
import variables as va
import variables.collections as vc
import variables.utils as vu
import vertex

import stdPi0s
import stdV0s
import stdPhotons

from glob import glob

import argparse

# get data type

parser = argparse.ArgumentParser(description='Sample type')
parser.add_argument('--inputfile', type=str, help='destination of output files for local run. ex) /home/belle2/junewoo/storage_b1/input_Jpsi/B02XJpsi_SKIM_00012_job227326565_00.udst.root')
parser.add_argument('--destination', type=str, help='destination of output files for local run. ex) /home/belle2/junewoo/storage_b1/output')

args = parser.parse_args()

basf2.conditions.prepend_globaltag(ma.getAnalysisGlobaltag())

# set random seed
basf2.set_random_seed(42)

fName = ""
if args.inputfile is not None:
    fName = args.inputfile
    if not fName.endswith(".root"): sys.exit(2)
    else:
        basename = os.path.basename(fName)
        name = os.path.splitext(basename)[0]
    print("=================="+fName+" is conducted=================")

destination = ""

if args.destination is not None:
    destination = args.destination

output_file = "Ntuple_" + "nominal"  + ".root"
output_file = destination + "/Ntuple/" + name + "_" + output_file

# Read uDST
my_path = basf2.create_path()

ma.inputMdst(environmentType='default',filename=fName,path=my_path)

# fill Upsilon(4S) MC particle
ma.fillParticleListFromMC('Upsilon(4S):PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)

# fill kaon only / pion only
ma.fillParticleList(decayString="K+:Kaon_base", cut="nCDCHits > 20 and dr < 0.5 and abs(dz) < 2 and nPXDHits > 0", path=my_path)
ma.fillParticleList(decayString="pi+:Pion_base", cut="nCDCHits > 20 and dr < 0.5 and abs(dz) < 2 and nPXDHits > 0", path=my_path)
ma.matchMCTruth("K+:Kaon_base", path=my_path)
ma.matchMCTruth("pi+:Pion_base", path=my_path)

# fill my charged particles
ma.cutAndCopyList("K+:mychargedKaon", "K+:Kaon_base", cut="kaonID>0.6 and nCDCHits > 20 and dr < 0.5 and abs(dz) < 2 and nPXDHits > 0", path=my_path)
ma.cutAndCopyList("pi+:mychargedPion", "pi+:Pion_base", cut="pionID>0.6 and nCDCHits > 20 and dr < 0.5 and abs(dz) < 2 and nPXDHits > 0", path=my_path)

# fill kaon/pion ID particle
ma.cutAndCopyList("K+:KaonID", "K+:Kaon_base", cut="kaonID>0.6", path=my_path)
ma.cutAndCopyList("pi+:PionID", "pi+:Pion_base", cut="pionID>0.6", path=my_path)

# select true particles
ma.cutAndCopyList("K+:Kaon_true", "K+:Kaon_base", cut="abs(mcPDG) == 321", path=my_path)
ma.cutAndCopyList("pi+:Pion_true", "pi+:Pion_base", cut="abs(mcPDG) == 211", path=my_path)

ma.cutAndCopyList("K+:mychargedKaon_true", "K+:mychargedKaon", cut="abs(mcPDG) == 321", path=my_path)
ma.cutAndCopyList("pi+:mychargedPion_true", "pi+:mychargedPion", cut="abs(mcPDG) == 211", path=my_path)

ma.cutAndCopyList("K+:KaonID_true", "K+:KaonID", cut="abs(mcPDG) == 321", path=my_path)
ma.cutAndCopyList("pi+:PionID_true", "pi+:PionID", cut="abs(mcPDG) == 211", path=my_path)

# for fragmentation correction
ma.fillParticleListFromMC('K+:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
ma.fillParticleListFromMC('K*+:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
ma.fillParticleListFromMC('pi+:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
ma.fillParticleListFromMC('e+:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
ma.fillParticleListFromMC('nu_e:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
ma.fillParticleListFromMC('nu_mu:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
ma.fillParticleListFromMC('nu_tau:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
ma.fillParticleListFromMC('gamma:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
ma.fillParticleListFromMC('pi0:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
ma.fillParticleListFromMC('Xsu:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
ma.fillParticleListFromMC('Xsd:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
ma.fillParticleListFromMC("K0:PrimaryMC", cut = "mcPrimary", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
ma.fillParticleListFromMC("K*0:PrimaryMC", cut = "mcPrimary", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
ma.fillParticleListFromMC("B+:PrimaryMC", cut = "mcPrimary", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
ma.fillParticleListFromMC("B0:PrimaryMC", cut = "mcPrimary", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
    
ma.reconstructMCDecay('K*+:MC_ch1 =direct=> K0:PrimaryMC pi+:PrimaryMC ?gamma', cut='', dmID = 0, path=my_path)
ma.reconstructMCDecay('K*+:MC_ch2 =direct=> K+:PrimaryMC pi0:PrimaryMC ?gamma', cut='', dmID = 0, path=my_path)
ma.reconstructMCDecay('K*0:MC_ch1 =direct=> K+:PrimaryMC pi-:PrimaryMC ?gamma', cut='', dmID = 0, path=my_path)
ma.reconstructMCDecay('K*0:MC_ch2 =direct=> K0:PrimaryMC pi0:PrimaryMC ?gamma', cut='', dmID = 0, path=my_path)
    
# B plus
ma.reconstructMCDecay('Xsu:MCch1 -> K+:PrimaryMC', cut='', dmID = 1, path=my_path)
ma.reconstructMCDecay('Xsu:MCch2 -> K+:PrimaryMC pi0:PrimaryMC', cut='', dmID = 2, path=my_path)
ma.reconstructMCDecay('Xsu:MCch3 -> K0:PrimaryMC pi+:PrimaryMC', cut='', dmID = 3, path=my_path)
ma.reconstructMCDecay('Xsu:MCch4 -> K+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC', cut='', dmID = 4, path=my_path)
ma.reconstructMCDecay('Xsu:MCch5 -> K0:PrimaryMC pi+:PrimaryMC pi0:PrimaryMC', cut='', dmID = 5, path=my_path)
ma.reconstructMCDecay('Xsu:MCch6 -> K+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC pi0:PrimaryMC', cut='', dmID = 6, path=my_path)
ma.reconstructMCDecay('Xsu:MCch7 -> K0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC', cut='', dmID = 7, path=my_path)
ma.reconstructMCDecay('Xsu:MCch8 -> K+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC', cut='', dmID = 8, path=my_path)
ma.reconstructMCDecay('Xsu:MCch9 -> K0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC pi0:PrimaryMC', cut='', dmID = 9, path=my_path)
# 2pi0
ma.reconstructMCDecay('Xsu:MCch10 -> K+:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 10, path=my_path)
ma.reconstructMCDecay('Xsu:MCch11 -> K0:PrimaryMC pi+:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 11, path=my_path)
ma.reconstructMCDecay('Xsu:MCch12 -> K+:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 12, path=my_path)
# 3K
ma.reconstructMCDecay('Xsu:MCch13 -> K+:PrimaryMC K-:PrimaryMC K+:PrimaryMC', cut='', dmID = 13, path=my_path)
ma.reconstructMCDecay('Xsu:MCch14 -> K+:PrimaryMC K-:PrimaryMC K0:PrimaryMC pi+:PrimaryMC', cut='', dmID = 14, path=my_path)
ma.reconstructMCDecay('Xsu:MCch15 -> K+:PrimaryMC K-:PrimaryMC K+:PrimaryMC pi0:PrimaryMC', cut='', dmID = 15, path=my_path)
    
ma.copyLists(outputListName="Xsu:MCcomb", inputListNames=["Xsu:MCch1", "Xsu:MCch2","Xsu:MCch3","Xsu:MCch4","Xsu:MCch5","Xsu:MCch6","Xsu:MCch7","Xsu:MCch8","Xsu:MCch9","Xsu:MCch10", "Xsu:MCch11","Xsu:MCch12","Xsu:MCch13","Xsu:MCch14","Xsu:MCch15"], path=my_path)
    
ma.reconstructMCDecay('B+:Kcharge_e =direct=> K+:PrimaryMC nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
ma.reconstructMCDecay('B+:Kcharge_mu =direct=> K+:PrimaryMC nu_mu:PrimaryMC anti-nu_mu:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
ma.reconstructMCDecay('B+:Kcharge_tau =direct=> K+:PrimaryMC nu_tau:PrimaryMC anti-nu_tau:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
ma.copyLists(outputListName="B+:Kcharge_total", inputListNames=["B+:Kcharge_e", "B+:Kcharge_mu", "B+:Kcharge_tau"], path=my_path)
    
ma.reconstructMCDecay('B+:Kstarcharge_ch1_e =direct=> K*+:MC_ch1 nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.reconstructMCDecay('B+:Kstarcharge_ch1_mu =direct=> K*+:MC_ch1 nu_mu:PrimaryMC anti-nu_mu:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.reconstructMCDecay('B+:Kstarcharge_ch1_tau =direct=> K*+:MC_ch1 nu_tau:PrimaryMC anti-nu_tau:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.copyLists(outputListName="B+:Kstarcharge_ch1_total", inputListNames=["B+:Kstarcharge_ch1_e", "B+:Kstarcharge_ch1_mu", "B+:Kstarcharge_ch1_tau"], path=my_path)
    
ma.reconstructMCDecay('B+:Kstarcharge_ch2_e =direct=> K*+:MC_ch2 nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.reconstructMCDecay('B+:Kstarcharge_ch2_mu =direct=> K*+:MC_ch2 nu_mu:PrimaryMC anti-nu_mu:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.reconstructMCDecay('B+:Kstarcharge_ch2_tau =direct=> K*+:MC_ch2 nu_tau:PrimaryMC anti-nu_tau:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.copyLists(outputListName="B+:Kstarcharge_ch2_total", inputListNames=["B+:Kstarcharge_ch2_e", "B+:Kstarcharge_ch2_mu", "B+:Kstarcharge_ch2_tau"], path=my_path)
    
# B zero
ma.reconstructMCDecay('Xsd:MCch16 -> K0:PrimaryMC', cut='', dmID = 16, path=my_path)
ma.reconstructMCDecay('Xsd:MCch17 -> K+:PrimaryMC pi-:PrimaryMC', cut='', dmID = 17, path=my_path)
ma.reconstructMCDecay('Xsd:MCch18 -> K0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 18, path=my_path)
ma.reconstructMCDecay('Xsd:MCch19 -> K+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC', cut='', dmID = 19, path=my_path)
ma.reconstructMCDecay('Xsd:MCch20 -> K0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC', cut='', dmID = 20, path=my_path)
ma.reconstructMCDecay('Xsd:MCch21 -> K+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC', cut='', dmID = 21, path=my_path)
ma.reconstructMCDecay('Xsd:MCch22 -> K0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC', cut='', dmID = 22, path=my_path)
ma.reconstructMCDecay('Xsd:MCch23 -> K+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC', cut='', dmID = 23, path=my_path)
ma.reconstructMCDecay('Xsd:MCch24 -> K0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC', cut='', dmID = 24, path=my_path)
# 2pi0
ma.reconstructMCDecay('Xsd:MCch25 -> K0:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 25, path=my_path)
ma.reconstructMCDecay('Xsd:MCch26 -> K+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 26, path=my_path)
ma.reconstructMCDecay('Xsd:MCch27 -> K0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 27, path=my_path)
# 3K
ma.reconstructMCDecay('Xsd:MCch28 -> K+:PrimaryMC K-:PrimaryMC K0:PrimaryMC', cut='', dmID = 28, path=my_path)
ma.reconstructMCDecay('Xsd:MCch29 -> K+:PrimaryMC K-:PrimaryMC K+:PrimaryMC pi-:PrimaryMC', cut='', dmID = 29, path=my_path)
ma.reconstructMCDecay('Xsd:MCch30 -> K+:PrimaryMC K-:PrimaryMC K0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 30, path=my_path)
    
ma.copyLists(outputListName="Xsd:MCcomb", inputListNames=["Xsd:MCch16", "Xsd:MCch17","Xsd:MCch18","Xsd:MCch19","Xsd:MCch20","Xsd:MCch21","Xsd:MCch22","Xsd:MCch23","Xsd:MCch24","Xsd:MCch25", "Xsd:MCch26","Xsd:MCch27","Xsd:MCch28","Xsd:MCch29","Xsd:MCch30"], path=my_path)
    
ma.reconstructMCDecay('B0:Kneutral_e =direct=> K0:PrimaryMC nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
ma.reconstructMCDecay('B0:Kneutral_mu =direct=> K0:PrimaryMC nu_mu:PrimaryMC anti-nu_mu:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
ma.reconstructMCDecay('B0:Kneutral_tau =direct=> K0:PrimaryMC nu_tau:PrimaryMC anti-nu_tau:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
ma.reconstructMCDecay('B0:Kneutral_e_oppose =direct=> anti-K0:PrimaryMC nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
ma.reconstructMCDecay('B0:Kneutral_mu_oppose =direct=> anti-K0:PrimaryMC nu_mu:PrimaryMC anti-nu_mu:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
ma.reconstructMCDecay('B0:Kneutral_tau_oppose =direct=> anti-K0:PrimaryMC nu_tau:PrimaryMC anti-nu_tau:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
ma.copyLists(outputListName="B0:Kneutral_total", inputListNames=["B0:Kneutral_e", "B0:Kneutral_mu", "B0:Kneutral_tau", "B0:Kneutral_e_oppose", "B0:Kneutral_mu_oppose", "B0:Kneutral_tau_oppose"], path=my_path)
    
ma.reconstructMCDecay('B0:Kstarneutral_ch1_e =direct=> K*0:MC_ch1 nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.reconstructMCDecay('B0:Kstarneutral_ch1_mu =direct=> K*0:MC_ch1 nu_mu:PrimaryMC anti-nu_mu:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.reconstructMCDecay('B0:Kstarneutral_ch1_tau =direct=> K*0:MC_ch1 nu_tau:PrimaryMC anti-nu_tau:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.copyLists(outputListName="B0:Kstarneutral_ch1_total", inputListNames=["B0:Kstarneutral_ch1_e", "B0:Kstarneutral_ch1_mu", "B0:Kstarneutral_ch1_tau"], path=my_path)
    
ma.reconstructMCDecay('B0:Kstarneutral_ch2_e =direct=> K*0:MC_ch2 nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.reconstructMCDecay('B0:Kstarneutral_ch2_mu =direct=> K*0:MC_ch2 nu_mu:PrimaryMC anti-nu_mu:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.reconstructMCDecay('B0:Kstarneutral_ch2_tau =direct=> K*0:MC_ch2 nu_tau:PrimaryMC anti-nu_tau:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.copyLists(outputListName="B0:Kstarneutral_ch2_total", inputListNames=["B0:Kstarneutral_ch2_e", "B0:Kstarneutral_ch2_mu", "B0:Kstarneutral_ch2_tau"], path=my_path)
    
    
# calculate q2 and costheta for systematic uncertainty
# B+:Kcharge_e
ma.reconstructMCDecay('B+:MC_Kstarcharge_e =direct=> K*+:PrimaryMC nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', path=my_path)
ma.reconstructMCDecay('B+:MC_Xsu_e =direct=> Xsu:PrimaryMC nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', path=my_path)
ma.variablesToExtraInfo("B+:Kcharge_e", {"constant(0)":"helicityangle"}, path=my_path)
ma.variablesToExtraInfo("B+:MC_Kstarcharge_e", {"cosHelicityAngle(0,0)":"helicityangle"}, path=my_path)
ma.variablesToExtraInfo("B+:MC_Xsu_e", {"constant(0)":"helicityangle"}, path=my_path)
ma.copyLists(outputListName="B+:MC_signal_total_e", inputListNames=["B+:Kcharge_e", "B+:MC_Kstarcharge_e", "B+:MC_Xsu_e"], path=my_path)
    
ma.copyLists(outputListName="B0:MC_Kneutral_total_e", inputListNames=["B0:Kneutral_e", "B0:Kneutral_e_oppose"], path=my_path)
ma.reconstructMCDecay('B0:MC_Kstarneutral_e =direct=> K*0:PrimaryMC nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', path=my_path)
ma.reconstructMCDecay('B0:MC_Xsd_e =direct=> Xsd:PrimaryMC nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', path=my_path)
ma.variablesToExtraInfo("B0:MC_Kneutral_total_e", {"constant(0)":"helicityangle"}, path=my_path)
ma.variablesToExtraInfo("B0:MC_Kstarneutral_e", {"cosHelicityAngle(0,0)":"helicityangle"}, path=my_path)
ma.variablesToExtraInfo("B0:MC_Xsd_e", {"constant(0)":"helicityangle"}, path=my_path)
ma.copyLists(outputListName="B0:MC_signal_total_e", inputListNames=["B0:MC_Kneutral_total_e", "B0:MC_Kstarneutral_e", "B0:MC_Xsd_e"], path=my_path)
    
ma.fillParticleListFromMC('nu_e:MC_signal', cut = 'mcPrimary and [isDescendantOfList(B+:MC_signal_total_e,1) or isDescendantOfList(B0:MC_signal_total_e,1)]', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path) 
    
decay_num = ["nParticlesInList(B+:Kcharge_total)", "nParticlesInList(B+:Kstarcharge_ch1_total)", "nParticlesInList(B+:Kstarcharge_ch2_total)", "nParticlesInList(Xsu:MCcomb)", "nParticlesInList(Xsu:MCch1)", "nParticlesInList(Xsu:MCch2)", "nParticlesInList(Xsu:MCch3)", "nParticlesInList(Xsu:MCch4)", "nParticlesInList(Xsu:MCch5)", "nParticlesInList(Xsu:MCch6)", "nParticlesInList(Xsu:MCch7)", "nParticlesInList(Xsu:MCch8)", "nParticlesInList(Xsu:MCch9)", "nParticlesInList(Xsu:MCch10)", "nParticlesInList(Xsu:MCch11)", "nParticlesInList(Xsu:MCch12)", "nParticlesInList(Xsu:MCch13)", "nParticlesInList(Xsu:MCch14)", "nParticlesInList(Xsu:MCch15)", "nParticlesInList(B0:Kneutral_total)", "nParticlesInList(B0:Kstarneutral_ch1_total)", "nParticlesInList(B0:Kstarneutral_ch2_total)", "nParticlesInList(Xsd:MCcomb)", "nParticlesInList(Xsd:MCch16)", "nParticlesInList(Xsd:MCch17)", "nParticlesInList(Xsd:MCch18)", "nParticlesInList(Xsd:MCch19)", "nParticlesInList(Xsd:MCch20)", "nParticlesInList(Xsd:MCch21)", "nParticlesInList(Xsd:MCch22)", "nParticlesInList(Xsd:MCch23)", "nParticlesInList(Xsd:MCch24)", "nParticlesInList(Xsd:MCch25)", "nParticlesInList(Xsd:MCch26)", "nParticlesInList(Xsd:MCch27)", "nParticlesInList(Xsd:MCch28)", "nParticlesInList(Xsd:MCch29)", "nParticlesInList(Xsd:MCch30)"]
systematics = ["nParticlesInList(nu_e:MC_signal)", "nParticlesInList(B+:MC_signal_total_e)", "nParticlesInList(B0:MC_signal_total_e)", "invMassInLists(nu_e:MC_signal)", "averageValueInList(B+:MC_signal_total_e, extraInfo(helicityangle))", "averageValueInList(B0:MC_signal_total_e, extraInfo(helicityangle))", "averageValueInList(B+:MC_signal_total_e, daughter(0, M))", "averageValueInList(B0:MC_signal_total_e, daughter(0, M))", "averageValueInList(B+:MC_signal_total_e, M)", "averageValueInList(B0:MC_signal_total_e, M)", "nParticlesInList(B+:PrimaryMC)", "nParticlesInList(B0:PrimaryMC)"]  

# make Ntuple
NMyChargedParticle = ["nParticlesInList(K+:mychargedKaon)", "nParticlesInList(pi+:mychargedPion)", "nParticlesInList(K+:mychargedKaon_true)", "nParticlesInList(pi+:mychargedPion_true)"]
NPIDChargedParticle = ["nParticlesInList(K+:KaonID)", "nParticlesInList(pi+:PionID)", "nParticlesInList(K+:KaonID_true)", "nParticlesInList(pi+:PionID_true)"]
NMCChargedParticle = ["nParticlesInList(K+:Kaon_true)", "nParticlesInList(pi+:Pion_true)"]
va.variables.addAlias("qsquared","formula(\
(mcDaughter(0, mcDaughter(1, E))+mcDaughter(0, mcDaughter(2, E)))^2\
- (mcDaughter(0, mcDaughter(1, px))+mcDaughter(0, mcDaughter(2, px)))^2\
- (mcDaughter(0, mcDaughter(1, py))+mcDaughter(0, mcDaughter(2, py)))^2\
- (mcDaughter(0, mcDaughter(1, pz))+mcDaughter(0, mcDaughter(2, pz)))^2\
)")
XsKinematics = ["mcDaughter(0, mcDaughter(0, M))", "qsquared"]
ma.variablesToNtuple(decayString="Upsilon(4S):PrimaryMC",variables=NMyChargedParticle + NPIDChargedParticle + NMCChargedParticle + XsKinematics + decay_num + systematics, filename=output_file, treename="Upsilon", path=my_path)

ma.variablesToNtuple(decayString="K+:Kaon_true", variables = ["p", "cosTheta", "kaonID", "pionID"] + ["averageValueInList(Upsilon(4S):PrimaryMC, qsquared)", "averageValueInList(Upsilon(4S):PrimaryMC, mcDaughter(0, mcDaughter(0, M)))"] + decay_num + systematics, filename=output_file, treename="Kaon", path=my_path)
ma.variablesToNtuple(decayString="pi+:Pion_true", variables = ["p", "cosTheta", "kaonID", "pionID"] + ["averageValueInList(Upsilon(4S):PrimaryMC, qsquared)", "averageValueInList(Upsilon(4S):PrimaryMC, mcDaughter(0, mcDaughter(0, M)))"] + decay_num + systematics, filename=output_file, treename="Pion", path=my_path)

# for debug
#kinematics = ["mcDaughter(0, E)", "mcDaughter(0, px)", "mcDaughter(0, py)", "mcDaughter(0, pz)", "mcDaughter(1, E)", "mcDaughter(1, px)", "mcDaughter(1, py)", "mcDaughter(1, pz)"]
#my_path.add_module("ParticlePrinter",listName="Upsilon(4S):PrimaryMC", variables=kinematics)
#ma.printMCParticles(onlyPrimaries=True, showProperties=True, showMomenta=True, path=my_path)

    
# progress
basf2.process(my_path)

# Print call statistics
print(basf2.statistics)
