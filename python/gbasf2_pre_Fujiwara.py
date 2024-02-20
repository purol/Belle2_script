#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# usage: basf2 MakeNtuple_multi.py "./20210402/evt-0.mdst"
# last: 2022-06-24-00

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
parser.add_argument('--inputsample', required=True, type=str, help='sample type. list) default, MC_SIGNAL, MC_BKG_BB, MC_Jpsi, MC_Jpsi_SIGNAL, MC_BKG_BB_for_Jpsi, run_dependent')
args = parser.parse_args()
if(args.inputsample!="default" and args.inputsample!="MC_BKG_BB" and args.inputsample!="MC_Jpsi" and args.inputsample!="MC_BKG_BB_for_Jpsi" and args.inputsample!="MC_SIGNAL" and args.inputsample!="MC_Jpsi_SIGNAL" and args.inputsample!="run_dependent"):
    exit(1)

basf2.conditions.prepend_globaltag(ma.getAnalysisGlobaltag())
if args.inputsample != "run_dependent":
    basf2.conditions.prepend_globaltag("hadronid_Fall2021_Official_v0")

# set random seed
basf2.set_random_seed(42)

# Read uDST
my_path = basf2.create_path()

inputfile="/group/belle2/dataprod/MC/SkimTraining/mixed_BGx1.mdst_000001_prod00009434_task10020000001.root"
ma.inputMdst(environmentType='default',filename=inputfile,path=my_path)

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

if(args.inputsample == "MC_BKG_BB_for_Jpsi"):
    ma.fillParticleListFromMC('B+:fromUpsilon','genMotherPDG==300553',addDaughters=True,path=my_path)
    ma.fillParticleListFromMC('B0:fromUpsilon','genMotherPDG==300553',addDaughters=True,path=my_path)
    ma.fillParticleListFromMC('mu+:fromJpsi','genMotherPDG==443',addDaughters=True, path=my_path)
    ma.fillParticleListFromMC('J/psi:direct','[isMCDescendantOfList(B+:fromUpsilon, -1) == 1] or [isMCDescendantOfList(B0:fromUpsilon, -1) == 1]',addDaughters=True,path=my_path)
    ma.fillParticleListFromMC('mu+:fromJpsi_confirm','isMCDescendantOfList(J/psi:direct,1)',addDaughters=True, path=my_path)
    ma.reconstructMCDecay('J/psi:mumu =direct=> mu+:fromJpsi_confirm mu-:fromJpsi_confirm', '', path=my_path)
    ma.applyEventCuts("nParticlesInList(J/psi:mumu) < 1", path=my_path)

# primary particle list
ma.fillParticleList(decayString="K+:mychargedKaon", cut="kaonID>0.6 and nCDCHits > 20 and dr < 0.5 and abs(dz) < 2 and nPXDHits > 0",path=my_path)
ma.fillParticleList(decayString="pi+:mychargedPion", cut="pionID>0.6 and nCDCHits > 20 and dr < 0.5 and abs(dz) < 2 and nPXDHits > 0",path=my_path)

if args.inputsample != "run_dependent":
    # get uncertainty and correction factors
    ma.matchMCTruth("K+:mychargedKaon", my_path)
    reweighter_K = basf2.register_module("ParticleWeighting")
    reweighter_K.set_name("ParticleWeighting_Kaon")
    reweighter_K.param("tableName", "ParticleReweighting:kaonID_eff_dstardkpipi_FixedThresh06")
    reweighter_K.param("particleList", "K+:mychargedKaon")
    my_path.add_module(reweighter_K)
    ma.matchMCTruth("pi+:mychargedPion", my_path)
    reweighter_Pi = basf2.register_module("ParticleWeighting")
    reweighter_Pi.set_name("ParticleWeighting_Pion")
    reweighter_Pi.param("tableName", "ParticleReweighting:pionID_eff_dstardkpipi_FixedThresh06")
    reweighter_Pi.param("particleList", "pi+:mychargedPion")
    my_path.add_module(reweighter_Pi)
    va.variables.addAlias("Kaon_mean","extraInfo(ParticleReweighting:kaonID_eff_dstardkpipi_FixedThresh06_data_MC_ratio)")
    va.variables.addAlias("Kaon_stat_up","extraInfo(ParticleReweighting:kaonID_eff_dstardkpipi_FixedThresh06_data_MC_uncertainty_stat_up)")
    va.variables.addAlias("Kaon_stat_dn","extraInfo(ParticleReweighting:kaonID_eff_dstardkpipi_FixedThresh06_data_MC_uncertainty_stat_dn)")
    va.variables.addAlias("Kaon_sys_up","extraInfo(ParticleReweighting:kaonID_eff_dstardkpipi_FixedThresh06_data_MC_uncertainty_sys_up)")
    va.variables.addAlias("Kaon_sys_dn","extraInfo(ParticleReweighting:kaonID_eff_dstardkpipi_FixedThresh06_data_MC_uncertainty_sys_dn)")
    va.variables.addAlias("Kaon_total_up","formula((Kaon_sys_up*Kaon_sys_up + Kaon_stat_up*Kaon_stat_up)^0.5)")
    va.variables.addAlias("Kaon_total_dn","formula((Kaon_sys_dn*Kaon_sys_dn + Kaon_stat_dn*Kaon_stat_dn)^0.5)")
    va.variables.addAlias("Kaon_rel_unce_up","formula(Kaon_total_up/Kaon_mean)")
    va.variables.addAlias("Kaon_rel_unce_dn","formula(Kaon_total_dn/Kaon_mean)")
    va.variables.addAlias("Pion_mean","extraInfo(ParticleReweighting:pionID_eff_dstardkpipi_FixedThresh06_data_MC_ratio)")
    va.variables.addAlias("Pion_stat_up","extraInfo(ParticleReweighting:pionID_eff_dstardkpipi_FixedThresh06_data_MC_uncertainty_stat_up)")
    va.variables.addAlias("Pion_stat_dn","extraInfo(ParticleReweighting:pionID_eff_dstardkpipi_FixedThresh06_data_MC_uncertainty_stat_dn)")
    va.variables.addAlias("Pion_sys_up","extraInfo(ParticleReweighting:pionID_eff_dstardkpipi_FixedThresh06_data_MC_uncertainty_sys_up)")
    va.variables.addAlias("Pion_sys_dn","extraInfo(ParticleReweighting:pionID_eff_dstardkpipi_FixedThresh06_data_MC_uncertainty_sys_dn)")
    va.variables.addAlias("Pion_total_up","formula((Pion_sys_up*Pion_sys_up + Pion_stat_up*Pion_stat_up)^0.5)")
    va.variables.addAlias("Pion_total_dn","formula((Pion_sys_dn*Pion_sys_dn + Pion_stat_dn*Pion_stat_dn)^0.5)")
    va.variables.addAlias("Pion_rel_unce_up","formula(Pion_total_up/Pion_mean)")
    va.variables.addAlias("Pion_rel_unce_dn","formula(Pion_total_dn/Pion_mean)")
else:
    va.variables.addAlias("Kaon_mean","constant(1)")
    va.variables.addAlias("Kaon_stat_up","constant(1)")
    va.variables.addAlias("Kaon_stat_dn","constant(1)")
    va.variables.addAlias("Kaon_sys_up","constant(1)")
    va.variables.addAlias("Kaon_sys_dn","constant(1)")
    va.variables.addAlias("Kaon_total_up","constant(1)")
    va.variables.addAlias("Kaon_total_dn","constant(1)")
    va.variables.addAlias("Kaon_rel_unce_up","constant(1)")
    va.variables.addAlias("Kaon_rel_unce_dn","constant(1)")
    va.variables.addAlias("Pion_mean","constant(1)")
    va.variables.addAlias("Pion_stat_up","constant(1)")
    va.variables.addAlias("Pion_stat_dn","constant(1)")
    va.variables.addAlias("Pion_sys_up","constant(1)")
    va.variables.addAlias("Pion_sys_dn","constant(1)")
    va.variables.addAlias("Pion_total_up","constant(1)")
    va.variables.addAlias("Pion_total_dn","constant(1)")
    va.variables.addAlias("Pion_rel_unce_up","constant(1)")
    va.variables.addAlias("Pion_rel_unce_dn","constant(1)")

stdV0s.stdKshorts(path=my_path)
ma.cutAndCopyList("K_S0:myKaonshort", "K_S0:merged", cut="[0.4876 < M < 0.5076] and significanceOfDistance > 50", path=my_path)

stdPi0s.stdPi0s(listtype="eff30_May2020",path=my_path)
ma.cutAndCopyList("pi0:myneutralPion", "pi0:eff30_May2020", cut="p > 0.4", path=my_path)

stdPhotons.stdPhotons(listtype="pi0eff30_May2020", path=my_path)
ma.cutAndCopyList("gamma:mygamma", "gamma:pi0eff30_May2020", cut="", path=my_path)

# non-primary particle
#ma.reconstructDecay(decayString="K*0:myneutralKaonstar -> K+:mychargedKaon pi-:mychargedPion", cut="",path=my_path)
#vertex.kFit("K*0:myneutralKaonstar", 0, path=my_path)

#ma.reconstructDecay(decayString="K*+:mychargedKaonstar -> K_S0:myKaonshort pi+:mychargedPion", cut="", dmID=1, path=my_path)


# Xsu
ma.reconstructDecay(decayString="@Xsu:ch0 -> K+:mychargedKaon", cut="M < 2", dmID=0, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch1 -> K+:mychargedKaon pi0:myneutralPion", cut="M < 2", dmID=1, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch2 -> K_S0:myKaonshort pi+:mychargedPion", cut="M < 2", dmID=2, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch3 -> K+:mychargedKaon pi+:mychargedPion pi-:mychargedPion", cut="M < 2", dmID=3, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch4 -> K_S0:myKaonshort pi+:mychargedPion pi0:myneutralPion", cut="M < 2", dmID=4, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch5 -> K+:mychargedKaon pi-:mychargedPion pi+:mychargedPion pi0:myneutralPion", cut="M < 2", dmID=5, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch6 -> K_S0:myKaonshort pi+:mychargedPion pi-:mychargedPion pi+:mychargedPion", cut="M < 2", dmID=6, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch7 -> K+:mychargedKaon pi-:mychargedPion pi+:mychargedPion pi-:mychargedPion pi+:mychargedPion", cut="M < 2", dmID=7, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch8 -> K_S0:myKaonshort pi+:mychargedPion pi-:mychargedPion pi+:mychargedPion pi0:myneutralPion", cut="M < 2", dmID=8, path=my_path)
# 2pi0
ma.reconstructDecay(decayString="@Xsu:ch9 -> K+:mychargedKaon pi0:myneutralPion  pi0:myneutralPion", cut="M < 2", dmID=9, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch10 -> K_S0:myKaonshort pi+:mychargedPion pi0:myneutralPion pi0:myneutralPion", cut="M < 2", dmID=10, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch11 -> K+:mychargedKaon pi+:mychargedPion pi-:mychargedPion pi0:myneutralPion pi0:myneutralPion", cut="M < 2", dmID=11, path=my_path)
# 3K
ma.reconstructDecay(decayString="@Xsu:ch12 -> K+:mychargedKaon K-:mychargedKaon K+:mychargedKaon", cut="M < 2", dmID=12, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch13 -> K+:mychargedKaon K-:mychargedKaon K_S0:myKaonshort pi+:mychargedPion", cut="M < 2", dmID=13, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch14 -> K+:mychargedKaon K-:mychargedKaon K+:mychargedKaon pi0:myneutralPion", cut="M < 2", dmID=14, path=my_path)


# Xsd
ma.reconstructDecay(decayString="@Xsd:ch0 -> K_S0:myKaonshort", cut="M < 2", dmID=0, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch1 -> K+:mychargedKaon pi-:mychargedPion", cut="M < 2", dmID=1, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch2 -> K_S0:myKaonshort pi0:myneutralPion", cut="M < 2", dmID=2, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch3 -> K+:mychargedKaon pi-:mychargedPion pi0:myneutralPion", cut="M < 2", dmID=3, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch4 -> K_S0:myKaonshort pi+:mychargedPion pi-:mychargedPion", cut="M < 2", dmID=4, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch5 -> K+:mychargedKaon pi-:mychargedPion pi+:mychargedPion pi-:mychargedPion", cut="M < 2", dmID=5, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch6 -> K_S0:myKaonshort pi+:mychargedPion pi-:mychargedPion pi0:myneutralPion", cut="M < 2", dmID=6, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch7 -> K+:mychargedKaon pi-:mychargedPion pi+:mychargedPion pi-:mychargedPion pi0:myneutralPion", cut="M < 2", dmID=7, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch8 -> K_S0:myKaonshort pi+:mychargedPion pi-:mychargedPion pi+:mychargedPion pi-:mychargedPion", cut="M < 2", dmID=8, path=my_path)
# 2pi0
ma.reconstructDecay(decayString="@Xsd:ch9 -> K_S0:myKaonshort pi0:myneutralPion pi0:myneutralPion", cut="M < 2", dmID=9, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch10 -> K+:mychargedKaon pi-:mychargedPion pi0:myneutralPion pi0:myneutralPion", cut="M < 2", dmID=10, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch11 -> K_S0:myKaonshort pi+:mychargedPion pi-:mychargedPion pi0:myneutralPion pi0:myneutralPion", cut="M < 2", dmID=11, path=my_path)
# 3K
ma.reconstructDecay(decayString="@Xsd:ch12 -> K+:mychargedKaon K-:mychargedKaon K_S0:myKaonshort", cut="M < 2", dmID=12, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch13 -> K+:mychargedKaon K-:mychargedKaon K+:mychargedKaon pi-:mychargedPion", cut="M < 2", dmID=13, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch14 -> K+:mychargedKaon K-:mychargedKaon K_S0:myKaonshort pi0:myneutralPion", cut="M < 2", dmID=14, path=my_path)

# uncertainty
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

# correction factor from Kaon PID
ma.variablesToExtraInfo("Xsu:ch0", {"daughter(0,Kaon_mean)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch1", {"daughter(0,Kaon_mean)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch2", {"constant(1)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch3", {"daughter(0,Kaon_mean)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch4", {"constant(1)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch5", {"daughter(0,Kaon_mean)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch6", {"constant(1)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch7", {"daughter(0,Kaon_mean)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch8", {"constant(1)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch9", {"daughter(0,Kaon_mean)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch10", {"constant(1)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch11", {"daughter(0,Kaon_mean)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch12", {"formula(daughter(0,Kaon_mean)*daughter(1,Kaon_mean)*daughter(2,Kaon_mean))" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch13", {"formula(daughter(0,Kaon_mean)*daughter(1,Kaon_mean))" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch14", {"formula(daughter(0,Kaon_mean)*daughter(1,Kaon_mean)*daughter(2,Kaon_mean))" : "Kaon_PID_correction"}, path=my_path)

ma.variablesToExtraInfo("Xsd:ch0", {"constant(1)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch1", {"daughter(0,Kaon_mean)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch2", {"constant(1)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch3", {"daughter(0,Kaon_mean)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch4", {"constant(1)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch5", {"daughter(0,Kaon_mean)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch6", {"constant(1)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch7", {"daughter(0,Kaon_mean)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch8", {"constant(1)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch9", {"constant(1)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch10", {"daughter(0,Kaon_mean)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch11", {"constant(1)" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch12", {"formula(daughter(0,Kaon_mean)*daughter(1,Kaon_mean))" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch13", {"formula(daughter(0,Kaon_mean)*daughter(1,Kaon_mean)*daughter(2,Kaon_mean))" : "Kaon_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch14", {"formula(daughter(0,Kaon_mean)*daughter(1,Kaon_mean))" : "Kaon_PID_correction"}, path=my_path)

# uncertainty from Kaon PID
ma.variablesToExtraInfo("Xsu:ch0", {"daughter(0,Kaon_rel_unce_up)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch1", {"daughter(0,Kaon_rel_unce_up)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch2", {"constant(0)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch3", {"daughter(0,Kaon_rel_unce_up)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch4", {"constant(0)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch5", {"daughter(0,Kaon_rel_unce_up)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch6", {"constant(0)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch7", {"daughter(0,Kaon_rel_unce_up)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch8", {"constant(0)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch9", {"daughter(0,Kaon_rel_unce_up)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch10", {"constant(0)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch11", {"daughter(0,Kaon_rel_unce_up)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch12", {"formula(daughter(0,Kaon_rel_unce_up)+daughter(1,Kaon_rel_unce_up)+daughter(2,Kaon_rel_unce_up))" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch13", {"formula(daughter(0,Kaon_rel_unce_up)+daughter(1,Kaon_rel_unce_up))" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch14", {"formula(daughter(0,Kaon_rel_unce_up)+daughter(1,Kaon_rel_unce_up)+daughter(2,Kaon_rel_unce_up))" : "Kaon_PID_rel_uncer_up"}, path=my_path)

ma.variablesToExtraInfo("Xsd:ch0", {"constant(0)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch1", {"daughter(0,Kaon_rel_unce_up)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch2", {"constant(0)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch3", {"daughter(0,Kaon_rel_unce_up)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch4", {"constant(0)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch5", {"daughter(0,Kaon_rel_unce_up)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch6", {"constant(0)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch7", {"daughter(0,Kaon_rel_unce_up)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch8", {"constant(0)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch9", {"constant(0)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch10", {"daughter(0,Kaon_rel_unce_up)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch11", {"constant(0)" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch12", {"formula(daughter(0,Kaon_rel_unce_up)+daughter(1,Kaon_rel_unce_up))" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch13", {"formula(daughter(0,Kaon_rel_unce_up)+daughter(1,Kaon_rel_unce_up)+daughter(2,Kaon_rel_unce_up))" : "Kaon_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch14", {"formula(daughter(0,Kaon_rel_unce_up)+daughter(1,Kaon_rel_unce_up))" : "Kaon_PID_rel_uncer_up"}, path=my_path)

ma.variablesToExtraInfo("Xsu:ch0", {"daughter(0,Kaon_rel_unce_dn)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch1", {"daughter(0,Kaon_rel_unce_dn)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch2", {"constant(0)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch3", {"daughter(0,Kaon_rel_unce_dn)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch4", {"constant(0)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch5", {"daughter(0,Kaon_rel_unce_dn)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch6", {"constant(0)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch7", {"daughter(0,Kaon_rel_unce_dn)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch8", {"constant(0)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch9", {"daughter(0,Kaon_rel_unce_dn)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch10", {"constant(0)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch11", {"daughter(0,Kaon_rel_unce_dn)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch12", {"formula(daughter(0,Kaon_rel_unce_dn)+daughter(1,Kaon_rel_unce_dn)+daughter(2,Kaon_rel_unce_dn))" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch13", {"formula(daughter(0,Kaon_rel_unce_dn)+daughter(1,Kaon_rel_unce_dn))" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch14", {"formula(daughter(0,Kaon_rel_unce_dn)+daughter(1,Kaon_rel_unce_dn)+daughter(2,Kaon_rel_unce_dn))" : "Kaon_PID_rel_uncer_dn"}, path=my_path)

ma.variablesToExtraInfo("Xsd:ch0", {"constant(0)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch1", {"daughter(0,Kaon_rel_unce_dn)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch2", {"constant(0)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch3", {"daughter(0,Kaon_rel_unce_dn)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch4", {"constant(0)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch5", {"daughter(0,Kaon_rel_unce_dn)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch6", {"constant(0)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch7", {"daughter(0,Kaon_rel_unce_dn)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch8", {"constant(0)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch9", {"constant(0)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch10", {"daughter(0,Kaon_rel_unce_dn)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch11", {"constant(0)" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch12", {"formula(daughter(0,Kaon_rel_unce_dn)+daughter(1,Kaon_rel_unce_dn))" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch13", {"formula(daughter(0,Kaon_rel_unce_dn)+daughter(1,Kaon_rel_unce_dn)+daughter(2,Kaon_rel_unce_dn))" : "Kaon_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch14", {"formula(daughter(0,Kaon_rel_unce_dn)+daughter(1,Kaon_rel_unce_dn))" : "Kaon_PID_rel_uncer_dn"}, path=my_path)

# correction factor from Pion PID
ma.variablesToExtraInfo("Xsu:ch0", {"constant(1)" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch1", {"constant(1)" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch2", {"daughter(1,Pion_mean)" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch3", {"formula(daughter(1,Pion_mean)*daughter(2,Pion_mean))" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch4", {"daughter(1,Pion_mean)" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch5", {"formula(daughter(1,Pion_mean)*daughter(2,Pion_mean))" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch6", {"formula(daughter(1,Pion_mean)*daughter(2,Pion_mean)*daughter(3,Pion_mean))" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch7", {"formula(daughter(1,Pion_mean)*daughter(2,Pion_mean)*daughter(3,Pion_mean)*daughter(4,Pion_mean))" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch8", {"formula(daughter(1,Pion_mean)*daughter(2,Pion_mean)*daughter(3,Pion_mean))" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch9", {"constant(1)" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch10", {"daughter(1,Pion_mean)" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch11", {"formula(daughter(1,Pion_mean)*daughter(2,Pion_mean))" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch12", {"constant(1)" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch13", {"daughter(3,Pion_mean)" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch14", {"constant(1)" : "Pion_PID_correction"}, path=my_path)

ma.variablesToExtraInfo("Xsd:ch0", {"constant(1)" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch1", {"daughter(1,Pion_mean)" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch2", {"constant(1)" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch3", {"daughter(1,Pion_mean)" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch4", {"formula(daughter(1,Pion_mean)*daughter(2,Pion_mean))" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch5", {"formula(daughter(1,Pion_mean)*daughter(2,Pion_mean)*daughter(3,Pion_mean))" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch6", {"formula(daughter(1,Pion_mean)*daughter(2,Pion_mean))" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch7", {"formula(daughter(1,Pion_mean)*daughter(2,Pion_mean)*daughter(3,Pion_mean))" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch8", {"formula(daughter(1,Pion_mean)*daughter(2,Pion_mean)*daughter(3,Pion_mean)*daughter(4,Pion_mean))" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch9", {"constant(1)" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch10", {"daughter(1,Pion_mean)" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch11", {"formula(daughter(1,Pion_mean)*daughter(2,Pion_mean))" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch12", {"constant(1)" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch13", {"daughter(3,Pion_mean)" : "Pion_PID_correction"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch14", {"constant(1)" : "Pion_PID_correction"}, path=my_path)

# uncertainty from Pion PID
ma.variablesToExtraInfo("Xsu:ch0", {"constant(0)" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch1", {"constant(0)" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch2", {"daughter(1,Pion_rel_unce_up)" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch3", {"formula(daughter(1,Pion_rel_unce_up)+daughter(2,Pion_rel_unce_up))" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch4", {"daughter(1,Pion_rel_unce_up)" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch5", {"formula(daughter(1,Pion_rel_unce_up)+daughter(2,Pion_rel_unce_up))" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch6", {"formula(daughter(1,Pion_rel_unce_up)+daughter(2,Pion_rel_unce_up)+daughter(3,Pion_rel_unce_up))" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch7", {"formula(daughter(1,Pion_rel_unce_up)+daughter(2,Pion_rel_unce_up)+daughter(3,Pion_rel_unce_up)+daughter(4,Pion_rel_unce_up))" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch8", {"formula(daughter(1,Pion_rel_unce_up)+daughter(2,Pion_rel_unce_up)+daughter(3,Pion_rel_unce_up))" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch9", {"constant(0)" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch10", {"daughter(1,Pion_rel_unce_up)" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch11", {"formula(daughter(1,Pion_rel_unce_up)+daughter(2,Pion_rel_unce_up))" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch12", {"constant(0)" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch13", {"daughter(3,Pion_rel_unce_up)" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch14", {"constant(0)" : "Pion_PID_rel_uncer_up"}, path=my_path)

ma.variablesToExtraInfo("Xsd:ch0", {"constant(0)" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch1", {"daughter(1,Pion_rel_unce_up)" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch2", {"constant(0)" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch3", {"daughter(1,Pion_rel_unce_up)" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch4", {"formula(daughter(1,Pion_rel_unce_up)+daughter(2,Pion_rel_unce_up))" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch5", {"formula(daughter(1,Pion_rel_unce_up)+daughter(2,Pion_rel_unce_up)+daughter(3,Pion_rel_unce_up))" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch6", {"formula(daughter(1,Pion_rel_unce_up)+daughter(2,Pion_rel_unce_up))" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch7", {"formula(daughter(1,Pion_rel_unce_up)+daughter(2,Pion_rel_unce_up)+daughter(3,Pion_rel_unce_up))" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch8", {"formula(daughter(1,Pion_rel_unce_up)+daughter(2,Pion_rel_unce_up)+daughter(3,Pion_rel_unce_up)+daughter(4,Pion_rel_unce_up))" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch9", {"constant(0)" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch10", {"daughter(1,Pion_rel_unce_up)" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch11", {"formula(daughter(1,Pion_rel_unce_up)+daughter(2,Pion_rel_unce_up))" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch12", {"constant(0)" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch13", {"daughter(3,Pion_rel_unce_up)" : "Pion_PID_rel_uncer_up"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch14", {"constant(0)" : "Pion_PID_rel_uncer_up"}, path=my_path)

ma.variablesToExtraInfo("Xsu:ch0", {"constant(0)" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch1", {"constant(0)" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch2", {"daughter(1,Pion_rel_unce_dn)" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch3", {"formula(daughter(1,Pion_rel_unce_dn)+daughter(2,Pion_rel_unce_dn))" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch4", {"daughter(1,Pion_rel_unce_dn)" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch5", {"formula(daughter(1,Pion_rel_unce_dn)+daughter(2,Pion_rel_unce_dn))" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch6", {"formula(daughter(1,Pion_rel_unce_dn)+daughter(2,Pion_rel_unce_dn)+daughter(3,Pion_rel_unce_dn))" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch7", {"formula(daughter(1,Pion_rel_unce_dn)+daughter(2,Pion_rel_unce_dn)+daughter(3,Pion_rel_unce_dn)+daughter(4,Pion_rel_unce_dn))" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch8", {"formula(daughter(1,Pion_rel_unce_dn)+daughter(2,Pion_rel_unce_dn)+daughter(3,Pion_rel_unce_dn))" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch9", {"constant(0)" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch10", {"daughter(1,Pion_rel_unce_dn)" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch11", {"formula(daughter(1,Pion_rel_unce_dn)+daughter(2,Pion_rel_unce_dn))" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch12", {"constant(0)" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch13", {"daughter(3,Pion_rel_unce_dn)" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsu:ch14", {"constant(0)" : "Pion_PID_rel_uncer_dn"}, path=my_path)

ma.variablesToExtraInfo("Xsd:ch0", {"constant(0)" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch1", {"daughter(1,Pion_rel_unce_dn)" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch2", {"constant(0)" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch3", {"daughter(1,Pion_rel_unce_dn)" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch4", {"formula(daughter(1,Pion_rel_unce_dn)+daughter(2,Pion_rel_unce_dn))" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch5", {"formula(daughter(1,Pion_rel_unce_dn)+daughter(2,Pion_rel_unce_dn)+daughter(3,Pion_rel_unce_dn))" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch6", {"formula(daughter(1,Pion_rel_unce_dn)+daughter(2,Pion_rel_unce_dn))" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch7", {"formula(daughter(1,Pion_rel_unce_dn)+daughter(2,Pion_rel_unce_dn)+daughter(3,Pion_rel_unce_dn))" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch8", {"formula(daughter(1,Pion_rel_unce_dn)+daughter(2,Pion_rel_unce_dn)+daughter(3,Pion_rel_unce_dn)+daughter(4,Pion_rel_unce_dn))" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch9", {"constant(0)" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch10", {"daughter(1,Pion_rel_unce_dn)" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch11", {"formula(daughter(1,Pion_rel_unce_dn)+daughter(2,Pion_rel_unce_dn))" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch12", {"constant(0)" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch13", {"daughter(3,Pion_rel_unce_dn)" : "Pion_PID_rel_uncer_dn"}, path=my_path)
ma.variablesToExtraInfo("Xsd:ch14", {"constant(0)" : "Pion_PID_rel_uncer_dn"}, path=my_path)


# combine
ma.copyLists(outputListName="Xsu:comb", inputListNames=["Xsu:ch0", "Xsu:ch1", "Xsu:ch2", "Xsu:ch3", "Xsu:ch4", "Xsu:ch5", "Xsu:ch6", "Xsu:ch7", "Xsu:ch8", "Xsu:ch9", "Xsu:ch10", "Xsu:ch11", "Xsu:ch12", "Xsu:ch13", "Xsu:ch14"], path=my_path)
ma.copyLists(outputListName="Xsd:comb", inputListNames=["Xsd:ch0", "Xsd:ch1", "Xsd:ch2", "Xsd:ch3", "Xsd:ch4", "Xsd:ch5", "Xsd:ch6", "Xsd:ch7", "Xsd:ch8", "Xsd:ch9", "Xsd:ch10", "Xsd:ch11", "Xsd:ch12", "Xsd:ch13", "Xsd:ch14"], path=my_path)

# B->X Jpsi reconstruction
if(args.inputsample == "MC_Jpsi" or args.inputsample == "MC_BKG_BB_for_Jpsi" or args.inputsample == "MC_Jpsi_SIGNAL"):
    ma.fillParticleList(decayString="mu+:mychargedMuon", cut="muonID>0.6 and dr < 0.5 and abs(dz) < 2", path=my_path)
    ma.reconstructDecay("J/psi:temp -> mu-:mychargedMuon mu+:mychargedMuon", cut="abs(dM) < 0.05", path = my_path)
    ma.reconstructDecay("B+:temp_sig -> Xsu:comb J/psi:temp", cut="Mbc > 5.25 and abs(deltaE) < 0.1", dmID=0, path = my_path)
    ma.reconstructDecay("B0:temp_sig -> Xsd:comb J/psi:temp", cut="Mbc > 5.25 and abs(deltaE) < 0.1", dmID=0, path = my_path)
    ma.reconstructDecay("Upsilon(4S):temp_withoutneutrino_charged -> B+:feiHadronic B-:temp_sig",cut ="", dmID = 0, path=my_path)
    ma.reconstructDecay("Upsilon(4S):temp_withoutneutrino_neutral_opposite_cp -> B0:feiHadronic anti-B0:temp_sig",cut ="", dmID = 1, path=my_path)
    ma.reconstructDecay("Upsilon(4S):temp_withoutneutrino_neutral_same_cp -> B0:feiHadronic B0:temp_sig",cut ="", dmID = 1, path=my_path)
    ma.copyLists(outputListName="Upsilon(4S):temp_withoutneutrino", inputListNames=["Upsilon(4S):temp_withoutneutrino_charged", "Upsilon(4S):temp_withoutneutrino_neutral_opposite_cp", "Upsilon(4S):temp_withoutneutrino_neutral_same_cp"], path=my_path)
    va.variables.addAlias("dM_Jpsi_rank","daughter(1,daughter(1, abs(dM)))")
    ma.rankByLowest(particleList="Upsilon(4S):temp_withoutneutrino", variable="dM_Jpsi_rank",allowMultiRank=True,outputVariable="Upsilon_Jpsi_rank",path=my_path)
    ma.cutAndCopyList('Upsilon(4S):temp_withoutneutrino_BCS', 'Upsilon(4S):temp_withoutneutrino','extraInfo(Upsilon_Jpsi_rank) == 1', path=my_path)
    
    ma.cutAndCopyList('mu+:fromUpsilonmychargedMuon', 'mu+:mychargedMuon','isDescendantOfList(Upsilon(4S):temp_withoutneutrino_BCS,3)', path=my_path)
    ma.cutAndCopyList('J/psi:fromUpsilontemp', 'J/psi:temp','isDescendantOfList(Upsilon(4S):temp_withoutneutrino_BCS,2)', path=my_path)
    ma.variablesToNtuple(decayString="Upsilon(4S):temp_withoutneutrino_BCS",variables=["extraInfo(decayModeID)", "daughter(0, extraInfo(decayModeID))", "daughter(1, daughter(0, extraInfo(decayModeID)))", "daughter(1, daughter(1, dM))", "daughter(1, Mbc)", "daughter(1, deltaE)", "nParticlesInList(mu+:fromUpsilonmychargedMuon)", "nParticlesInList(J/psi:fromUpsilontemp)", "nParticlesInList(Upsilon(4S):temp_withoutneutrino_BCS)", "nParticlesInList(Upsilon(4S):temp_withoutneutrino)"],filename="Jpsi.root",treename="Jpsi",path=my_path)
    ma.applyEventCuts("nParticlesInList(Upsilon(4S):temp_withoutneutrino_BCS) > 0", path=my_path)
    ma.applyCuts("B+:feiHadronic","isDescendantOfList(Upsilon(4S):temp_withoutneutrino_BCS,1)",path=my_path)
    ma.applyCuts("B0:feiHadronic","isDescendantOfList(Upsilon(4S):temp_withoutneutrino_BCS,1)",path=my_path)
    ma.applyCuts("Xsu:comb","isDescendantOfList(Upsilon(4S):temp_withoutneutrino_BCS,2)",path=my_path)
    ma.applyCuts("Xsd:comb","isDescendantOfList(Upsilon(4S):temp_withoutneutrino_BCS,2)",path=my_path)

# tag side
track_selection = "dr < 0.5 and abs(dz) < 2"
cluster_selection = '[[clusterReg==1 and E>0.08] or \
                      [clusterReg==2 and E>0.03] or \
                      [clusterReg==3 and E>0.06]]'
if(args.inputsample == "MC_Jpsi" or args.inputsample == "MC_BKG_BB_for_Jpsi" or args.inputsample == "MC_Jpsi_SIGNAL"):
    track_selection += " and sourceObjectIsInList(mu+:fromUpsilonmychargedMuon)==0"
    cluster_selection += " and sourceObjectIsInList(mu+:fromUpsilonmychargedMuon)==0"
cleanMask = ("cleanMask",track_selection, cluster_selection)
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
ma.appendROEMasks("Upsilon(4S):withoutneutrino",[cleanMask],path=my_path)
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

# Least Cuts
# ma.applyCuts("Upsilon(4S):withoutneutrino", "[nROE_Tracks(cleanMask) < 0.5] and [roeEextra(cleanMask) < 5] and [daughter(0, chiProb) > -0.5]", path=my_path)

Btag_vars = vu.create_aliases(list_of_variables = ["chiProb", "extraInfo(SignalProbability)"], wrapper = "daughter(0,{variable})",prefix="Btag")

Bsig_vars = vu.create_aliases(list_of_variables = [ "daughter(0, extraInfo(decayModeID))"], wrapper = "daughter(1,{variable})", prefix="Bsig")

U_vars = ["extraInfo(decayModeID)", "nROE_Tracks(cleanMask)", "roeEextra(cleanMask)"]

# Ntuple output
output_file = "Ntuple.root"
ma.variablesToNtuple(decayString="Upsilon(4S):withoutneutrino",variables=Btag_vars,filename=output_file,treename="Btag",path=my_path)
ma.variablesToNtuple(decayString="Upsilon(4S):withoutneutrino",variables=Bsig_vars,filename=output_file,treename="Bsig",path=my_path)
ma.variablesToNtuple(decayString="Upsilon(4S):withoutneutrino",variables=U_vars,filename=output_file,treename="Upsilon",path=my_path)

# progress
basf2.process(my_path)

# Print call statistics
print(basf2.statistics)
