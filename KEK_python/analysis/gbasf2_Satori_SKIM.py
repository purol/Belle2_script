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
parser.add_argument('--inputsample', required=True, type=str, help='sample type. list) default, MC_SIGNAL, MC_BKG_BB, MC_Jpsi, MC_Jpsi_SIGNAL, MC_BKG_BB_for_Jpsi, data_Jpsi, run_dependent')
args = parser.parse_args()
if(args.inputsample!="default" and args.inputsample!="MC_BKG_BB" and args.inputsample!="MC_Jpsi" and args.inputsample!="MC_BKG_BB_for_Jpsi" and args.inputsample!="MC_SIGNAL" and args.inputsample!="MC_Jpsi_SIGNAL" and args.inputsample!="data_Jpsi" and args.inputsample!="run_dependent"):
    exit(1)

basf2.conditions.prepend_globaltag(ma.getAnalysisGlobaltag())

# only for proc13 & prompt data (on and off both)
#basf2.conditions.prepend_globaltag("data_beam_conditions_proc13prompt")
#basf2.conditions.prepend_globaltag('data_trigger_conditions_proc13')

# set random seed
basf2.set_random_seed(42)

hashmapName = "hashmap_Upsilon.root"
JpsiName = "Jpsi.root"
output_file = "SKIM.root"

# Read uDST
my_path = basf2.create_path()

inputfile="/group/belle2/dataprod/MC/SkimTraining/mixed_BGx1.mdst_000001_prod00009434_task10020000001.root"
ma.inputMdst(environmentType='default',filename=inputfile,path=my_path)

# only for proc13 & prompt data (on and off both)
#ma.correctEnergyBias(inputListNames=["gamma:all"], tableName="PhotonEnergyBiasCorrection_MC15ri_Nov2022", path=my_path)
#ma.correctEnergyBias(inputListNames=["gamma:all"], tableName="PhotonEnergyBiasCorrection_MC15rd_June2023", path=my_path)

if(args.inputsample=="MC_BKG_BB"):
    # reject signal sample
    # Load particles from MCParticle at first
    ma.fillParticleListFromMC('K+:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
    ma.fillParticleListFromMC('K*+:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
    ma.fillParticleListFromMC('nu_e:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
    ma.fillParticleListFromMC('nu_mu:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
    ma.fillParticleListFromMC('nu_tau:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
    ma.fillParticleListFromMC('Xsu:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
    ma.fillParticleListFromMC('Xsd:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
    ma.fillParticleListFromMC("K0:PrimaryMC", cut = "mcPrimary", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
    ma.fillParticleListFromMC("K*0:PrimaryMC", cut = "mcPrimary", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
    ma.fillParticleListFromMC("Upsilon(4S):PrimaryMC", cut = "mcPrimary", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
    
    ma.reconstructMCDecay('B+:Kcharge_e =direct=> K+:PrimaryMC nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', path=my_path)
    ma.reconstructMCDecay('B+:Kcharge_mu =direct=> K+:PrimaryMC nu_mu:PrimaryMC anti-nu_mu:PrimaryMC ?gamma', cut='', path=my_path)
    ma.reconstructMCDecay('B+:Kcharge_tau =direct=> K+:PrimaryMC nu_tau:PrimaryMC anti-nu_tau:PrimaryMC ?gamma', cut='', path=my_path)
    ma.copyLists(outputListName="B+:Kcharge_total", inputListNames=["B+:Kcharge_e", "B+:Kcharge_mu", "B+:Kcharge_tau"], path=my_path)
    
    ma.reconstructMCDecay('B+:Kstarcharge_e =direct=> K*+:PrimaryMC nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', path=my_path)
    ma.reconstructMCDecay('B+:Kstarcharge_mu =direct=> K*+:PrimaryMC nu_mu:PrimaryMC anti-nu_mu:PrimaryMC ?gamma', cut='', path=my_path)
    ma.reconstructMCDecay('B+:Kstarcharge_tau =direct=> K*+:PrimaryMC nu_tau:PrimaryMC anti-nu_tau:PrimaryMC ?gamma', cut='', path=my_path)
    ma.copyLists(outputListName="B+:Kstarcharge_total", inputListNames=["B+:Kstarcharge_e", "B+:Kstarcharge_mu", "B+:Kstarcharge_tau"], path=my_path)
    
    ma.reconstructMCDecay('B+:Xsu_e =direct=> Xsu:PrimaryMC nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', path=my_path)
    ma.reconstructMCDecay('B+:Xsu_mu =direct=> Xsu:PrimaryMC nu_mu:PrimaryMC anti-nu_mu:PrimaryMC ?gamma', cut='', path=my_path)
    ma.reconstructMCDecay('B+:Xsu_tau =direct=> Xsu:PrimaryMC nu_tau:PrimaryMC anti-nu_tau:PrimaryMC ?gamma', cut='', path=my_path)
    ma.copyLists(outputListName="B+:Xsu_total", inputListNames=["B+:Xsu_e", "B+:Xsu_mu", "B+:Xsu_tau"], path=my_path)
    
    ma.copyLists(outputListName="B+:SIGNAL_total", inputListNames=["B+:Kcharge_total", "B+:Kstarcharge_total", "B+:Xsu_total"], path=my_path)
    
    ma.reconstructMCDecay('B0:Kneutral_e =direct=> K0:PrimaryMC nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', path=my_path)
    ma.reconstructMCDecay('B0:Kneutral_mu =direct=> K0:PrimaryMC nu_mu:PrimaryMC anti-nu_mu:PrimaryMC ?gamma', cut='', path=my_path)
    ma.reconstructMCDecay('B0:Kneutral_tau =direct=> K0:PrimaryMC nu_tau:PrimaryMC anti-nu_tau:PrimaryMC ?gamma', cut='', path=my_path)
    ma.reconstructMCDecay('B0:Kneutral_e_oppose =direct=> anti-K0:PrimaryMC nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', path=my_path)
    ma.reconstructMCDecay('B0:Kneutral_mu_oppose =direct=> anti-K0:PrimaryMC nu_mu:PrimaryMC anti-nu_mu:PrimaryMC ?gamma', cut='', path=my_path)
    ma.reconstructMCDecay('B0:Kneutral_tau_oppose =direct=> anti-K0:PrimaryMC nu_tau:PrimaryMC anti-nu_tau:PrimaryMC ?gamma', cut='', path=my_path)
    ma.copyLists(outputListName="B0:Kneutral_total", inputListNames=["B0:Kneutral_e", "B0:Kneutral_mu", "B0:Kneutral_tau", "B0:Kneutral_e_oppose", "B0:Kneutral_mu_oppose", "B0:Kneutral_tau_oppose"], path=my_path)
    
    ma.reconstructMCDecay('B0:Kstarneutral_e =direct=> K*0:PrimaryMC nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', path=my_path)
    ma.reconstructMCDecay('B0:Kstarneutral_mu =direct=> K*0:PrimaryMC nu_mu:PrimaryMC anti-nu_mu:PrimaryMC ?gamma', cut='', path=my_path)
    ma.reconstructMCDecay('B0:Kstarneutral_tau =direct=> K*0:PrimaryMC nu_tau:PrimaryMC anti-nu_tau:PrimaryMC ?gamma', cut='', path=my_path)
    ma.copyLists(outputListName="B0:Kstarneutral_total", inputListNames=["B0:Kstarneutral_e", "B0:Kstarneutral_mu", "B0:Kstarneutral_tau"], path=my_path)
    
    ma.reconstructMCDecay('B0:Xsd_e =direct=> Xsd:PrimaryMC nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', path=my_path)
    ma.reconstructMCDecay('B0:Xsd_mu =direct=> Xsd:PrimaryMC nu_mu:PrimaryMC anti-nu_mu:PrimaryMC ?gamma', cut='', path=my_path)
    ma.reconstructMCDecay('B0:Xsd_tau =direct=> Xsd:PrimaryMC nu_tau:PrimaryMC anti-nu_tau:PrimaryMC ?gamma', cut='', path=my_path)
    ma.copyLists(outputListName="B0:Xsd_total", inputListNames=["B0:Xsd_e", "B0:Xsd_mu", "B0:Xsd_tau"], path=my_path)
    
    ma.copyLists(outputListName="B0:SIGNAL_total", inputListNames=["B0:Kneutral_total", "B0:Kstarneutral_total", "B0:Xsd_total"], path=my_path)
    
    ma.applyEventCuts("nParticlesInList(B+:SIGNAL_total) < 1 and nParticlesInList(B0:SIGNAL_total) < 1", path=my_path)

# remove B->Xs Jpsi(->mu mu) from BB background
if(args.inputsample == "MC_BKG_BB_for_Jpsi"):
    ma.fillParticleListFromMC('B+:fromUpsilon','genMotherPDG==300553',addDaughters=True,path=my_path)
    ma.fillParticleListFromMC('B0:fromUpsilon','genMotherPDG==300553',addDaughters=True,path=my_path)
    ma.fillParticleListFromMC('mu+:fromJpsi','genMotherPDG==443',addDaughters=True, path=my_path)
    ma.fillParticleListFromMC('J/psi:direct','[isMCDescendantOfList(B+:fromUpsilon, -1) == 1] or [isMCDescendantOfList(B0:fromUpsilon, -1) == 1]',addDaughters=True,path=my_path)
    ma.fillParticleListFromMC('mu+:fromJpsi_confirm','isMCDescendantOfList(J/psi:direct,1)',addDaughters=True, path=my_path)
    ma.reconstructMCDecay('J/psi:mumu_BKG_temp =direct=> mu+:fromJpsi_confirm mu-:fromJpsi_confirm', '', path=my_path)
    ma.applyEventCuts("nParticlesInList(J/psi:mumu_BKG_temp) < 1", path=my_path)

# remove non B->Xs Jpsi(->mu mu) from BB background
if(args.inputsample == "MC_Jpsi_SIGNAL"):
    ma.fillParticleListFromMC('B+:fromUpsilon','genMotherPDG==300553',addDaughters=True,path=my_path)
    ma.fillParticleListFromMC('B0:fromUpsilon','genMotherPDG==300553',addDaughters=True,path=my_path)
    ma.fillParticleListFromMC('mu+:fromJpsi','genMotherPDG==443',addDaughters=True, path=my_path)
    ma.fillParticleListFromMC('J/psi:direct','[isMCDescendantOfList(B+:fromUpsilon, -1) == 1] or [isMCDescendantOfList(B0:fromUpsilon, -1) == 1]',addDaughters=True,path=my_path)
    ma.fillParticleListFromMC('mu+:fromJpsi_confirm','isMCDescendantOfList(J/psi:direct,1)',addDaughters=True, path=my_path)
    ma.reconstructMCDecay('J/psi:mumu_SIGNAL_temp =direct=> mu+:fromJpsi_confirm mu-:fromJpsi_confirm', '', path=my_path)
    ma.applyEventCuts("nParticlesInList(J/psi:mumu_SIGNAL_temp) > 0", path=my_path)

# fix wrong dec file
if(args.inputsample == "MC_BKG_BB_for_Jpsi" or args.inputsample == "MC_Jpsi_SIGNAL"):
    ma.fillParticleListFromMC('eta_c:Jpsifix_temp', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
    ma.fillParticleListFromMC('psi(2S):Jpsifix_temp', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
    ma.fillParticleListFromMC('chi_c0:Jpsifix_temp', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
    ma.fillParticleListFromMC("K0:Jpsifix_temp", cut = "mcPrimary", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
    ma.fillParticleListFromMC("J/psi:Jpsifix_temp", cut = "mcPrimary", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
    ma.fillParticleListFromMC('K+:Jpsifix_temp', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
    ma.fillParticleListFromMC('pi+:Jpsifix_temp', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)

    # following decay modes should be removed
    ma.reconstructMCDecay('B0:Jpsifix_temp_ch1 =direct=> eta_c:Jpsifix_temp K0:Jpsifix_temp', cut='', path=my_path)
    ma.reconstructMCDecay('B0:Jpsifix_temp_ch2 =direct=> psi(2S):Jpsifix_temp K0:Jpsifix_temp', cut='', path=my_path)
    ma.reconstructMCDecay('B0:Jpsifix_temp_ch3 =direct=> chi_c0:Jpsifix_temp K0:Jpsifix_temp', cut='', path=my_path)
    ma.reconstructMCDecay('anti-B0:Jpsifix_temp_ch4 =direct=> eta_c:Jpsifix_temp anti-K0:Jpsifix_temp', cut='', path=my_path)
    ma.reconstructMCDecay('anti-B0:Jpsifix_temp_ch5 =direct=> psi(2S):Jpsifix_temp anti-K0:Jpsifix_temp', cut='', path=my_path)
    ma.reconstructMCDecay('anti-B0:Jpsifix_temp_ch6 =direct=> chi_c0:Jpsifix_temp anti-K0:Jpsifix_temp', cut='', path=my_path)
    ma.applyEventCuts("[nParticlesInList(B0:Jpsifix_temp_ch1) < 1] and [nParticlesInList(B0:Jpsifix_temp_ch2) < 1] and [nParticlesInList(B0:Jpsifix_temp_ch3) < 1] and [nParticlesInList(anti-B0:Jpsifix_temp_ch4) < 1] and [nParticlesInList(anti-B0:Jpsifix_temp_ch5) < 1] and [nParticlesInList(anti-B0:Jpsifix_temp_ch6) < 1]", path=my_path)

    # following decay mode should be scaled (0.000811420 -> 0.000055137)
    ma.reconstructMCDecay('B+:Jpsifix_temp_ch7 =direct=> J/psi:Jpsifix_temp K+:Jpsifix_temp pi+:Jpsifix_temp pi-:Jpsifix_temp', cut='', path=my_path)
    ma.applyEventCuts("eventRandom < formula(0.06795124596^nParticlesInList(B+:Jpsifix_temp_ch7))", path=my_path)

# primary particle list
stdPhotons.stdPhotons(listtype="all", beamBackgroundMVAWeight="MC15ri", fakePhotonMVAWeight="MC15ri", path=my_path)

ma.fillParticleList(decayString="K+:mychargedKaon", cut="kaonID>0.6 and nCDCHits > 20 and dr < 0.5 and abs(dz) < 2 and nPXDHits > 0",path=my_path)
ma.fillParticleList(decayString="pi+:mychargedPion", cut="pionID>0.6 and nCDCHits > 20 and dr < 0.5 and abs(dz) < 2 and nPXDHits > 0",path=my_path)

stdV0s.stdKshorts(path=my_path)
ma.cutAndCopyList("K_S0:myKaonshort", "K_S0:merged", cut="[0.4876 < M < 0.5076] and significanceOfDistance > 50", path=my_path)

stdPi0s.stdPi0s(listtype="eff30_May2020",path=my_path)
ma.cutAndCopyList("pi0:myneutralPion", "pi0:eff30_May2020", cut="p > 0.4", path=my_path)

stdPhotons.stdPhotons(listtype="pi0eff30_May2020", path=my_path)
ma.cutAndCopyList("gamma:mygamma", "gamma:pi0eff30_May2020", cut="", path=my_path)

# Xsu
ma.reconstructDecay(decayString="@Xsu:ch0 -> K+:mychargedKaon", cut="M < 2.8", dmID=0, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch1 -> K+:mychargedKaon pi0:myneutralPion", cut="M < 2.8", dmID=1, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch2 -> K_S0:myKaonshort pi+:mychargedPion", cut="M < 2.8", dmID=2, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch3 -> K+:mychargedKaon pi+:mychargedPion pi-:mychargedPion", cut="M < 2.8", dmID=3, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch4 -> K_S0:myKaonshort pi+:mychargedPion pi0:myneutralPion", cut="M < 2.8", dmID=4, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch5 -> K+:mychargedKaon pi-:mychargedPion pi+:mychargedPion pi0:myneutralPion", cut="M < 2.8", dmID=5, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch6 -> K_S0:myKaonshort pi+:mychargedPion pi-:mychargedPion pi+:mychargedPion", cut="M < 2.8", dmID=6, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch7 -> K+:mychargedKaon pi-:mychargedPion pi+:mychargedPion pi-:mychargedPion pi+:mychargedPion", cut="M < 2.8", dmID=7, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch8 -> K_S0:myKaonshort pi+:mychargedPion pi-:mychargedPion pi+:mychargedPion pi0:myneutralPion", cut="M < 2.8", dmID=8, path=my_path)
# 2pi0
ma.reconstructDecay(decayString="@Xsu:ch9 -> K+:mychargedKaon pi0:myneutralPion  pi0:myneutralPion", cut="M < 2.8", dmID=9, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch10 -> K_S0:myKaonshort pi+:mychargedPion pi0:myneutralPion pi0:myneutralPion", cut="M < 2.8", dmID=10, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch11 -> K+:mychargedKaon pi+:mychargedPion pi-:mychargedPion pi0:myneutralPion pi0:myneutralPion", cut="M < 2.8", dmID=11, path=my_path)
# 3K
ma.reconstructDecay(decayString="@Xsu:ch12 -> K+:mychargedKaon K-:mychargedKaon K+:mychargedKaon", cut="M < 2.8", dmID=12, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch13 -> K+:mychargedKaon K-:mychargedKaon K_S0:myKaonshort pi+:mychargedPion", cut="M < 2.8", dmID=13, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch14 -> K+:mychargedKaon K-:mychargedKaon K+:mychargedKaon pi0:myneutralPion", cut="M < 2.8", dmID=14, path=my_path)


# Xsd
ma.reconstructDecay(decayString="@Xsd:ch0 -> K_S0:myKaonshort", cut="M < 2.8", dmID=0, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch1 -> K+:mychargedKaon pi-:mychargedPion", cut="M < 2.8", dmID=1, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch2 -> K_S0:myKaonshort pi0:myneutralPion", cut="M < 2.8", dmID=2, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch3 -> K+:mychargedKaon pi-:mychargedPion pi0:myneutralPion", cut="M < 2.8", dmID=3, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch4 -> K_S0:myKaonshort pi+:mychargedPion pi-:mychargedPion", cut="M < 2.8", dmID=4, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch5 -> K+:mychargedKaon pi-:mychargedPion pi+:mychargedPion pi-:mychargedPion", cut="M < 2.8", dmID=5, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch6 -> K_S0:myKaonshort pi+:mychargedPion pi-:mychargedPion pi0:myneutralPion", cut="M < 2.8", dmID=6, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch7 -> K+:mychargedKaon pi-:mychargedPion pi+:mychargedPion pi-:mychargedPion pi0:myneutralPion", cut="M < 2.8", dmID=7, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch8 -> K_S0:myKaonshort pi+:mychargedPion pi-:mychargedPion pi+:mychargedPion pi-:mychargedPion", cut="M < 2.8", dmID=8, path=my_path)
# 2pi0
ma.reconstructDecay(decayString="@Xsd:ch9 -> K_S0:myKaonshort pi0:myneutralPion pi0:myneutralPion", cut="M < 2.8", dmID=9, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch10 -> K+:mychargedKaon pi-:mychargedPion pi0:myneutralPion pi0:myneutralPion", cut="M < 2.8", dmID=10, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch11 -> K_S0:myKaonshort pi+:mychargedPion pi-:mychargedPion pi0:myneutralPion pi0:myneutralPion", cut="M < 2.8", dmID=11, path=my_path)
# 3K
ma.reconstructDecay(decayString="@Xsd:ch12 -> K+:mychargedKaon K-:mychargedKaon K_S0:myKaonshort", cut="M < 2.8", dmID=12, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch13 -> K+:mychargedKaon K-:mychargedKaon K+:mychargedKaon pi-:mychargedPion", cut="M < 2.8", dmID=13, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch14 -> K+:mychargedKaon K-:mychargedKaon K_S0:myKaonshort pi0:myneutralPion", cut="M < 2.8", dmID=14, path=my_path)

# uncertainty of KS0
ma.variablesToExtraInfo("Xsu:ch0", {"constant(0)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch1", {"constant(0)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch2", {"daughter(0,distance)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch3", {"constant(0)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch4", {"daughter(0,distance)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch5", {"constant(0)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch6", {"daughter(0,distance)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch7", {"constant(0)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch8", {"daughter(0,distance)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch9", {"constant(0)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch10", {"daughter(0,distance)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch11", {"constant(0)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch12", {"constant(0)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch13", {"daughter(2,distance)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch14", {"constant(0)" : "KS0_3D_distance"}, path=my_path)

ma.variablesToExtraInfo("Xsd:ch0", {"daughter(0,distance)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch1", {"constant(0)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch2", {"daughter(0,distance)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch3", {"constant(0)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch4", {"daughter(0,distance)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch5", {"constant(0)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch6", {"daughter(0,distance)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch7", {"constant(0)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch8", {"daughter(0,distance)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch9", {"daughter(0,distance)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch10", {"constant(0)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch11", {"daughter(0,distance)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch12", {"daughter(2,distance)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch13", {"constant(0)" : "KS0_3D_distance"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch14", {"daughter(2,distance)" : "KS0_3D_distance"}, path=my_path)

# combine
ma.copyLists(outputListName="Xsu:comb", inputListNames=["Xsu:ch0", "Xsu:ch1", "Xsu:ch2", "Xsu:ch3", "Xsu:ch4", "Xsu:ch5", "Xsu:ch6", "Xsu:ch7", "Xsu:ch8", "Xsu:ch9", "Xsu:ch10", "Xsu:ch11", "Xsu:ch12", "Xsu:ch13", "Xsu:ch14"], path=my_path)
ma.copyLists(outputListName="Xsd:comb", inputListNames=["Xsd:ch0", "Xsd:ch1", "Xsd:ch2", "Xsd:ch3", "Xsd:ch4", "Xsd:ch5", "Xsd:ch6", "Xsd:ch7", "Xsd:ch8", "Xsd:ch9", "Xsd:ch10", "Xsd:ch11", "Xsd:ch12", "Xsd:ch13", "Xsd:ch14"], path=my_path)

# Ntuple for preselection of Xs
#if(args.inputsample=="MC_SIGNAL") or (args.inputsample=="MC_Jpsi_SIGNAL"):
#    ma.cutAndCopyList('Xsu:random', 'Xsu:comb','eventRandom > -1.0', path=my_path)
#    ma.cutAndCopyList('Xsd:random', 'Xsd:comb','eventRandom > -1.0', path=my_path)
#else:
#    ma.cutAndCopyList('Xsu:random', 'Xsu:comb','eventRandom > 0.95', path=my_path)
#    ma.cutAndCopyList('Xsd:random', 'Xsd:comb','eventRandom > 0.95', path=my_path)
#ma.variablesToNtuple(decayString="Xsu:random",variables=["M"],filename=output_file,treename="Xsu_preselection",path=my_path)
#ma.variablesToNtuple(decayString="Xsd:random",variables=["M"],filename=output_file,treename="Xsd_preselection",path=my_path)

# Cut on Mxs (preselection)
ma.applyCuts("Xsu:comb", "[M < 2.0]", path=my_path)
ma.applyCuts("Xsd:comb", "[M < 2.0]", path=my_path)

# B->X Jpsi reconstruction
if(args.inputsample == "MC_Jpsi" or args.inputsample == "MC_BKG_BB_for_Jpsi" or args.inputsample == "MC_Jpsi_SIGNAL" or args.inputsample == "data_Jpsi"):
    ma.fillParticleList(decayString="mu+:mychargedMuon", cut="[muonID>0.1] and [dr < 2] and [abs(dz) < 4]", path=my_path)
    ma.reconstructDecay("J/psi:temp_for_plot -> mu-:mychargedMuon mu+:mychargedMuon", cut="abs(dM) < 0.1", path = my_path)
    ma.cutAndCopyList("J/psi:temp", "J/psi:temp_for_plot", "abs(dM) < 0.05",path=my_path)

    # apply event cut to obtan the proper number of event for XsJpsi analysis (without signal embedding)
    ma.applyEventCuts("nParticlesInList(J/psi:temp) > 0", path=my_path)

    ma.reconstructDecay("B+:temp_sig -> Xsu:comb J/psi:temp", cut="[Mbc > 5.25] and [abs(deltaE) < 0.1]", path = my_path)
    ma.reconstructDecay("B0:temp_sig -> Xsd:comb J/psi:temp", cut="[Mbc > 5.25] and [abs(deltaE) < 0.1]", path = my_path)
    ma.reconstructDecay("Upsilon(4S):temp_withoutneutrino_charged -> B+:feiHadronic B-:temp_sig",cut ="", dmID = 0, path=my_path)
    ma.reconstructDecay("Upsilon(4S):temp_withoutneutrino_neutral_opposite_cp -> B0:feiHadronic anti-B0:temp_sig",cut ="", dmID = 1, path=my_path)
    ma.reconstructDecay("Upsilon(4S):temp_withoutneutrino_neutral_same_cp -> B0:feiHadronic B0:temp_sig",cut ="", dmID = 1, path=my_path)
    ma.copyLists(outputListName="Upsilon(4S):temp_withoutneutrino", inputListNames=["Upsilon(4S):temp_withoutneutrino_charged", "Upsilon(4S):temp_withoutneutrino_neutral_opposite_cp", "Upsilon(4S):temp_withoutneutrino_neutral_same_cp"], path=my_path)
    ma.cutAndCopyList("J/psi:temp_daughterUpsilon", "J/psi:temp", "isDescendantOfList(Upsilon(4S):temp_withoutneutrino,-1)",path=my_path)
    ma.rankByLowest(particleList="J/psi:temp_daughterUpsilon", variable="random",allowMultiRank=True,outputVariable="Jpsi_random_rank",path=my_path) 
    ma.cutAndCopyList("J/psi:temp_BCS", "J/psi:temp_daughterUpsilon", "extraInfo(Jpsi_random_rank) == 1",path=my_path)
    ma.cutAndCopyList("mu+:fromJpsimuchargedMuon", "mu+:mychargedMuon", "isDescendantOfList(J/psi:temp_BCS,1) == 1",path=my_path)
    ma.applyEventCuts("nParticlesInList(J/psi:temp_BCS) > 0", path=my_path)

    # remove Xsu/Xsd which include muons
    ma.reconstructDecay("B+:temp_sig_JpsiBCS -> Xsu:comb J/psi:temp_BCS", cut="", path = my_path)
    ma.reconstructDecay("B0:temp_sig_JpsiBCS -> Xsd:comb J/psi:temp_BCS", cut="", path = my_path)
    ma.applyCuts("Xsu:comb","isDescendantOfList(B+:temp_sig_JpsiBCS,1)",path=my_path)
    ma.applyCuts("Xsd:comb","isDescendantOfList(B0:temp_sig_JpsiBCS,1)",path=my_path)

    # remove Btag which include muons
    ma.reconstructDecay("B+:feiHadronic_JpsiBCS -> B+:feiHadronic J/psi:temp_BCS", cut="", path = my_path)
    ma.reconstructDecay("B0:feiHadronic_JpsiBCS -> B0:feiHadronic J/psi:temp_BCS", cut="", path = my_path)
    ma.applyCuts("B+:feiHadronic","isDescendantOfList(B+:feiHadronic_JpsiBCS,1)",path=my_path)
    ma.applyCuts("B0:feiHadronic","isDescendantOfList(B0:feiHadronic_JpsiBCS,1)",path=my_path)

    # save Jpsi info into Ntuple
    #ma.variablesToNtuple(decayString="J/psi:temp_for_plot",variables=["dM", "M"],filename=JpsiName,treename="Jpsiforplot",path=my_path)
    #ma.variablesToNtuple(decayString="J/psi:temp_BCS",variables=["dM", "M", "nParticlesInList(mu+:fromJpsimuchargedMuon)", "nParticlesInList(J/psi:temp_daughterUpsilon)", "nParticlesInList(J/psi:temp)"],filename=JpsiName,treename="JpsiBCS",path=my_path)
    #ma.variablesToNtuple(decayString="Upsilon(4S):temp_withoutneutrino",variables=["extraInfo(decayModeID)", "daughter(1, daughter(1, dM))", "daughter(1, Mbc)", "daughter(1, deltaE)"],filename=JpsiName,treename="Upsilon_Jpsi",path=my_path)

# mask definition
track_selection = "[dr < 0.5] and [abs(dz) < 2]"
cluster_selection = '[[clusterReg==1 and E>0.08] or \
                      [clusterReg==2 and E>0.05] or \
                      [clusterReg==3 and E>0.06]] and [abs(clusterTiming) < formula(2.0 * clusterErrorTiming)] and [abs(clusterTiming) < 200] and [minC2TDist > 20]'
cluster_selection_v133 = '[[clusterReg==1 and E>0.08] or \
                           [clusterReg==2 and E>0.05] or \
                           [clusterReg==3 and E>0.06]] and [beamBackgroundSuppression > 0.5] and [fakePhotonSuppression > 0.3]'
cluster_selection_v200 = '[[clusterReg==1 and E>0.10] or \
                          [clusterReg==2 and E>0.06] or \
                          [clusterReg==3 and E>0.15]] and [minC2TDist > 50] and [inCDCAcceptance]'
cluster_selection_v200_800 = '[[clusterReg==1 and E>formula(0.10/0.800)] or \
                               [clusterReg==2 and E>formula(0.06/0.800)] or \
                               [clusterReg==3 and E>formula(0.15/0.800)]] and [minC2TDist > 50] and [inCDCAcceptance]'
cluster_selection_v200_825 = '[[clusterReg==1 and E>formula(0.10/0.825)] or \
                               [clusterReg==2 and E>formula(0.06/0.825)] or \
                               [clusterReg==3 and E>formula(0.15/0.825)]] and [minC2TDist > 50] and [inCDCAcceptance]'
cluster_selection_v200_850 = '[[clusterReg==1 and E>formula(0.10/0.850)] or \
                               [clusterReg==2 and E>formula(0.06/0.850)] or \
                               [clusterReg==3 and E>formula(0.15/0.850)]] and [minC2TDist > 50] and [inCDCAcceptance]'
cluster_selection_v200_875 = '[[clusterReg==1 and E>formula(0.10/0.875)] or \
                               [clusterReg==2 and E>formula(0.06/0.875)] or \
                               [clusterReg==3 and E>formula(0.15/0.875)]] and [minC2TDist > 50] and [inCDCAcceptance]'
cluster_selection_v200_900 = '[[clusterReg==1 and E>formula(0.10/0.900)] or \
                               [clusterReg==2 and E>formula(0.06/0.900)] or \
                               [clusterReg==3 and E>formula(0.15/0.900)]] and [minC2TDist > 50] and [inCDCAcceptance]'
cluster_selection_v200_925 = '[[clusterReg==1 and E>formula(0.10/0.925)] or \
                               [clusterReg==2 and E>formula(0.06/0.925)] or \
                               [clusterReg==3 and E>formula(0.15/0.925)]] and [minC2TDist > 50] and [inCDCAcceptance]'
cluster_selection_v200_950 = '[[clusterReg==1 and E>formula(0.10/0.950)] or \
                               [clusterReg==2 and E>formula(0.06/0.950)] or \
                               [clusterReg==3 and E>formula(0.15/0.950)]] and [minC2TDist > 50] and [inCDCAcceptance]'
cluster_selection_v200_975 = '[[clusterReg==1 and E>formula(0.10/0.975)] or \
                               [clusterReg==2 and E>formula(0.06/0.975)] or \
                               [clusterReg==3 and E>formula(0.15/0.975)]] and [minC2TDist > 50] and [inCDCAcceptance]'
cluster_selection_v200_025 = '[[clusterReg==1 and E>formula(0.10/1.025)] or \
                               [clusterReg==2 and E>formula(0.06/1.025)] or \
                               [clusterReg==3 and E>formula(0.15/1.025)]] and [minC2TDist > 50] and [inCDCAcceptance]'
cluster_selection_v200_050 = '[[clusterReg==1 and E>formula(0.10/1.050)] or \
                               [clusterReg==2 and E>formula(0.06/1.050)] or \
                               [clusterReg==3 and E>formula(0.15/1.050)]] and [minC2TDist > 50] and [inCDCAcceptance]'
cluster_selection_v200_075 = '[[clusterReg==1 and E>formula(0.10/1.075)] or \
                               [clusterReg==2 and E>formula(0.06/1.075)] or \
                               [clusterReg==3 and E>formula(0.15/1.075)]] and [minC2TDist > 50] and [inCDCAcceptance]'
cluster_selection_v200_100 = '[[clusterReg==1 and E>formula(0.10/1.100)] or \
                               [clusterReg==2 and E>formula(0.06/1.100)] or \
                               [clusterReg==3 and E>formula(0.15/1.100)]] and [minC2TDist > 50] and [inCDCAcceptance]'
cluster_minimum_energy ='[[clusterReg==1 and E>0.08] or \
                           [clusterReg==2 and E>0.05] or \
                           [clusterReg==3 and E>0.06]]'
if(args.inputsample == "MC_Jpsi" or args.inputsample == "MC_BKG_BB_for_Jpsi" or args.inputsample == "MC_Jpsi_SIGNAL" or args.inputsample == "data_Jpsi"):
    track_selection += " and [sourceObjectIsInList(mu+:fromJpsimuchargedMuon)==0]"
    cluster_selection += " and [sourceObjectIsInList(mu+:fromJpsimuchargedMuon)==0]"
    cluster_selection_v133 += " and [sourceObjectIsInList(mu+:fromJpsimuchargedMuon)==0]"
    cluster_selection_v200 += " and [sourceObjectIsInList(mu+:fromJpsimuchargedMuon)==0]"
    cluster_selection_v200_800 += " and [sourceObjectIsInList(mu+:fromJpsimuchargedMuon)==0]"
    cluster_selection_v200_825 += " and [sourceObjectIsInList(mu+:fromJpsimuchargedMuon)==0]"
    cluster_selection_v200_850 += " and [sourceObjectIsInList(mu+:fromJpsimuchargedMuon)==0]"
    cluster_selection_v200_875 += " and [sourceObjectIsInList(mu+:fromJpsimuchargedMuon)==0]"
    cluster_selection_v200_900 += " and [sourceObjectIsInList(mu+:fromJpsimuchargedMuon)==0]"
    cluster_selection_v200_925 += " and [sourceObjectIsInList(mu+:fromJpsimuchargedMuon)==0]"
    cluster_selection_v200_950 += " and [sourceObjectIsInList(mu+:fromJpsimuchargedMuon)==0]"
    cluster_selection_v200_975 += " and [sourceObjectIsInList(mu+:fromJpsimuchargedMuon)==0]"
    cluster_selection_v200_025 += " and [sourceObjectIsInList(mu+:fromJpsimuchargedMuon)==0]"
    cluster_selection_v200_050 += " and [sourceObjectIsInList(mu+:fromJpsimuchargedMuon)==0]"
    cluster_selection_v200_075 += " and [sourceObjectIsInList(mu+:fromJpsimuchargedMuon)==0]"
    cluster_selection_v200_100 += " and [sourceObjectIsInList(mu+:fromJpsimuchargedMuon)==0]"
cleanMask = ("cleanMask",track_selection, cluster_selection_v200)

loose_track_selection = "[dr < 2] and [abs(dz) < 4]"
if(args.inputsample == "MC_Jpsi" or args.inputsample == "MC_BKG_BB_for_Jpsi" or args.inputsample == "MC_Jpsi_SIGNAL" or args.inputsample == "data_Jpsi"):
    loose_track_selection += " and [sourceObjectIsInList(mu+:fromJpsimuchargedMuon)==0]"
looseMask = ("looseMask", loose_track_selection, cluster_selection_v200)

# tag side
ma.buildRestOfEvent("B+:feiHadronic",path=my_path)
ma.appendROEMasks("B+:feiHadronic",[cleanMask],path=my_path)
ma.buildContinuumSuppression(list_name="B+:feiHadronic",roe_mask = "cleanMask", path=my_path)
ma.buildRestOfEvent("B0:feiHadronic",path=my_path)
ma.appendROEMasks("B0:feiHadronic",[cleanMask],path=my_path)
ma.buildContinuumSuppression(list_name="B0:feiHadronic",roe_mask = "cleanMask", path=my_path)

# signal side
ma.reconstructDecay("B+:sig -> Xsu:comb ?nu", cut="", dmID=0, path = my_path)
ma.buildRestOfEvent("B+:sig",path=my_path)
ma.appendROEMasks("B+:sig",[cleanMask],path=my_path)
ma.buildContinuumSuppression(list_name="B+:sig",roe_mask = "cleanMask", path=my_path)
ma.reconstructDecay("B0:sig -> Xsd:comb ?nu", cut="", dmID=0, path = my_path)
ma.buildRestOfEvent("B0:sig",path=my_path)
ma.appendROEMasks("B0:sig",[cleanMask],path=my_path)
ma.buildContinuumSuppression(list_name="B0:sig",roe_mask = "cleanMask", path=my_path)

# Upsilon(4S) -> Btag K
ma.reconstructDecay("Upsilon(4S):withoutneutrino_charged -> B+:feiHadronic B-:sig",cut ="", dmID = 0, path=my_path)
ma.reconstructDecay("Upsilon(4S):withoutneutrino_neutral_opposite_cp -> B0:feiHadronic anti-B0:sig",cut ="", dmID = 1, path=my_path)
ma.reconstructDecay("Upsilon(4S):withoutneutrino_neutral_same_cp -> B0:feiHadronic B0:sig",cut ="", dmID = 1, path=my_path)
ma.copyLists(outputListName="Upsilon(4S):withoutneutrino", inputListNames=["Upsilon(4S):withoutneutrino_charged", "Upsilon(4S):withoutneutrino_neutral_opposite_cp", "Upsilon(4S):withoutneutrino_neutral_same_cp"], path=my_path)
ma.buildRestOfEvent("Upsilon(4S):withoutneutrino", inputParticlelists=["K_S0:myKaonshort", "pi0:myneutralPion", "gamma:mygamma"], path=my_path)
    
# apply cut: no charged track, E, cluster on ROE
ma.appendROEMasks("Upsilon(4S):withoutneutrino",[cleanMask, looseMask],path=my_path)
# ma.applyCuts("Upsilon(4S):withoutneutrino","nROE_Tracks(cleanMask)==0",path=my_path)

# Get (psig + pmiss)^2 variable
va.variables.addAlias("MsquaredBsig_op0","formula( (weMissE(cleanMask,0)+daughter(1, useCMSFrame(E)))^2 - (weMissPx(cleanMask,0)+daughter(1, useCMSFrame(px)))^2 - (weMissPy(cleanMask,0)+daughter(1, useCMSFrame(py)))^2 - (weMissPz(cleanMask,0)+daughter(1, useCMSFrame(pz)))^2 )")
va.variables.addAlias("MsquaredBsig_op1","formula( (weMissE(cleanMask,1)+daughter(1, useCMSFrame(E)))^2 - (weMissPx(cleanMask,1)+daughter(1, useCMSFrame(px)))^2 - (weMissPy(cleanMask,1)+daughter(1, useCMSFrame(py)))^2 - (weMissPz(cleanMask,1)+daughter(1, useCMSFrame(pz)))^2 )")
va.variables.addAlias("MsquaredBsig_op2","formula( (weMissE(cleanMask,2)+daughter(1, useCMSFrame(E)))^2 - (weMissPx(cleanMask,2)+daughter(1, useCMSFrame(px)))^2 - (weMissPy(cleanMask,2)+daughter(1, useCMSFrame(py)))^2 - (weMissPz(cleanMask,2)+daughter(1, useCMSFrame(pz)))^2 )")
va.variables.addAlias("MsquaredBsig_op3","formula( (weMissE(cleanMask,3)+daughter(1, useCMSFrame(E)))^2 - (weMissPx(cleanMask,3)+daughter(1, useCMSFrame(px)))^2 - (weMissPy(cleanMask,3)+daughter(1, useCMSFrame(py)))^2 - (weMissPz(cleanMask,3)+daughter(1, useCMSFrame(pz)))^2 )")
va.variables.addAlias("MsquaredBsig_op4","formula( (weMissE(cleanMask,4)+daughter(1, useCMSFrame(E)))^2 - (weMissPx(cleanMask,4)+daughter(1, useCMSFrame(px)))^2 - (weMissPy(cleanMask,4)+daughter(1, useCMSFrame(py)))^2 - (weMissPz(cleanMask,4)+daughter(1, useCMSFrame(pz)))^2 )")
va.variables.addAlias("MsquaredBsig_op7","formula( (weMissE(cleanMask,7)+daughter(1, useCMSFrame(E)))^2 - (weMissPx(cleanMask,7)+daughter(1, useCMSFrame(px)))^2 - (weMissPy(cleanMask,7)+daughter(1, useCMSFrame(py)))^2 - (weMissPz(cleanMask,7)+daughter(1, useCMSFrame(pz)))^2 )")
va.variables.addAlias("qsquared","formula( daughter(1,M)*daughter(1,M)-constant(2)*useTagSideRecoilRestFrame(daughter(0, eRecoil), 0)*useTagSideRecoilRestFrame(daughter(1, E), 0) )")

# best candidate selection
va.variables.addAlias("SignalProbofBtag_rank","daughter(0,extraInfo(SignalProbability))")
ma.rankByHighest(particleList="Upsilon(4S):withoutneutrino", variable="SignalProbofBtag_rank",allowMultiRank=True,outputVariable="Upsilon_rank",path=my_path)
# ma.applyCuts("Upsilon(4S):withoutneutrino","extraInfo(Upsilon_rank) == 1",path=my_path)

# Ntuple for preselection
#Btag_preselection_vars = vu.create_aliases(list_of_variables = ["extraInfo(SignalProbability)"], wrapper = "daughter(0,{variable})",prefix="Btag")
#Bsig_preselection_vars = vu.create_aliases(list_of_variables = [ "daughter(0, extraInfo(decayModeID))"], wrapper = "daughter(1,{variable})", prefix="Bsig")
#U_preselection_vars = ["extraInfo(decayModeID)", "nROE_Tracks(looseMask)", "nROE_ParticlesInList(pi0:myneutralPion)"]
#if(args.inputsample=="MC_SIGNAL") or (args.inputsample=="MC_Jpsi_SIGNAL"):
#    ma.cutAndCopyList('Upsilon(4S):random', 'Upsilon(4S):withoutneutrino','eventRandom > -1.0', path=my_path)
#else:
#    ma.cutAndCopyList('Upsilon(4S):random', 'Upsilon(4S):withoutneutrino','eventRandom > 0.95', path=my_path)
#ma.variablesToNtuple(decayString="Upsilon(4S):random",variables=Btag_preselection_vars,filename=output_file,treename="Btag_preselection",path=my_path)
#ma.variablesToNtuple(decayString="Upsilon(4S):random",variables=Bsig_preselection_vars,filename=output_file,treename="Bsig_preselection",path=my_path)
#ma.variablesToNtuple(decayString="Upsilon(4S):random",variables=U_preselection_vars,filename=output_file,treename="Upsilon_preselection",path=my_path)

# Least Cuts (preselection)
ma.applyCuts("Upsilon(4S):withoutneutrino", "[nROE_Tracks(cleanMask) < 0.5] and [nROE_ParticlesInList(pi0:myneutralPion) < 0.5] and [nROE_ParticlesInList(K_S0:myKaonshort) < 0.5]", path=my_path)

# Lets Skim
ma.applyEventCuts("nParticlesInList(Upsilon(4S):withoutneutrino) > 0", path=my_path)
ma.removeParticlesNotInLists(lists_to_keep=["B+:feiHadronic", "B0:feiHadronic"], path=my_path)
ma.outputUdst(filename=output_file,particleLists=["B+:feiHadronic", "B0:feiHadronic"],path=my_path)

# progress
basf2.process(my_path)

# Print call statistics
print(basf2.statistics)
