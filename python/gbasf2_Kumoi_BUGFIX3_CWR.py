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

# side bands
parser.add_argument("--nominal", action="store_true") # run with nominal analysis
parser.add_argument("--MbcSide", action="store_true") # run Mbc < 5.26 GeB
parser.add_argument("--chargeSide", action="store_true") # run wrong charge sideband
parser.add_argument("--MuonSide", action="store_true") # run MuonID side band
parser.add_argument("--ElectronSide", action="store_true") # run EID side band
parser.add_argument("--Jpsi", action="store_true") # run J/psi analysis
parser.add_argument("--ppbar", action="store_true") # run Xs p pbar analysis
parser.add_argument("--KSKS", action="store_true") # run Xs KS0 KS0 analysis

# data or MC
parser.add_argument('--sampletype', required=True, type=str, help='sample type. list) data, MC_qq, MC_BB, MC_SIGNAL')

args = parser.parse_args()

if((args.sampletype!="data") and (args.sampletype!="MC_qq") and (args.sampletype!="MC_BB") and (args.sampletype!="MC_SIGNAL")):
    print("please use proper sampletype")
    exit(1)

# set run option
Options = []
if(args.nominal):
    Options.append("nominal")
if(args.MbcSide):
    Options.append("MbcSide")
if(args.chargeSide):
    Options.append("chargeSide")
if(args.MuonSide):
    Options.append("MuonSide")
if(args.ElectronSide):
    Options.append("ElectronSide")
if(args.Jpsi):
    Options.append("Jpsi")
if(args.ppbar):
    Options.append("ppbar")
if(args.KSKS):
    Options.append("KSKS")

basf2.conditions.prepend_globaltag(ma.getAnalysisGlobaltag())

# only for proc13 & prompt data (on and off both)
#basf2.conditions.prepend_globaltag("data_beam_conditions_proc13prompt")
#basf2.conditions.prepend_globaltag('data_trigger_conditions_proc13')

# set random seed
basf2.set_random_seed(42)

for analysistype in Options:
    
    hashmapName = "hashmap_Upsilon.root"
    JpsiName = "Jpsi.root"
    output_file = "Ntuple_" + analysistype  + ".root"
    
    # Read uDST
    my_path = basf2.create_path()
    
    inputfile="/group/belle2/dataprod/MC/SkimTraining/mixed_BGx1.mdst_000001_prod00009434_task10020000001.root"
    ma.inputMdst(environmentType='default',filename=inputfile,path=my_path)
    
    # only for proc13 & prompt data (on and off both)
    #ma.fillParticleList(decayString="gamma:all", cut="",path=my_path)
    #ma.correctEnergyBias(inputListNames=["gamma:all"], tableName="PhotonEnergyBiasCorrection_MC15ri_Nov2022", path=my_path)
    #ma.correctEnergyBias(inputListNames=["gamma:all"], tableName="PhotonEnergyBiasCorrection_MC15rd_June2023", path=my_path)
    
    if((args.sampletype=="MC_BB") and (analysistype=="nominal")):
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
    if((args.sampletype=="MC_BB") and (analysistype=="Jpsi")):
        ma.fillParticleListFromMC('B+:fromUpsilon','genMotherPDG==300553',addDaughters=True,path=my_path)
        ma.fillParticleListFromMC('B0:fromUpsilon','genMotherPDG==300553',addDaughters=True,path=my_path)
        ma.fillParticleListFromMC('mu+:fromJpsi','genMotherPDG==443',addDaughters=True, path=my_path)
        ma.fillParticleListFromMC('J/psi:direct','[isMCDescendantOfList(B+:fromUpsilon, -1) == 1] or [isMCDescendantOfList(B0:fromUpsilon, -1) == 1]',addDaughters=True,path=my_path)
        ma.fillParticleListFromMC('mu+:fromJpsi_confirm','isMCDescendantOfList(J/psi:direct,1)',addDaughters=True, path=my_path)
        ma.reconstructMCDecay('J/psi:mumu_BKG_temp =direct=> mu+:fromJpsi_confirm mu-:fromJpsi_confirm', '', path=my_path)
        ma.applyEventCuts("nParticlesInList(J/psi:mumu_BKG_temp) < 1", path=my_path)
    
    # remove non B->Xs Jpsi(->mu mu) from BB background
    if((args.sampletype=="MC_SIGNAL") and (analysistype=="Jpsi")):
        ma.fillParticleListFromMC('B+:fromUpsilon','genMotherPDG==300553',addDaughters=True,path=my_path)
        ma.fillParticleListFromMC('B0:fromUpsilon','genMotherPDG==300553',addDaughters=True,path=my_path)
        ma.fillParticleListFromMC('mu+:fromJpsi','genMotherPDG==443',addDaughters=True, path=my_path)
        ma.fillParticleListFromMC('J/psi:direct','[isMCDescendantOfList(B+:fromUpsilon, -1) == 1] or [isMCDescendantOfList(B0:fromUpsilon, -1) == 1]',addDaughters=True,path=my_path)
        ma.fillParticleListFromMC('mu+:fromJpsi_confirm','isMCDescendantOfList(J/psi:direct,1)',addDaughters=True, path=my_path)
        ma.reconstructMCDecay('J/psi:mumu_SIGNAL_temp =direct=> mu+:fromJpsi_confirm mu-:fromJpsi_confirm', '', path=my_path)
        ma.applyEventCuts("nParticlesInList(J/psi:mumu_SIGNAL_temp) > 0", path=my_path)
    
    # fix wrong dec file
    if( (analysistype=="Jpsi") and ((args.sampletype=="MC_BB") or (args.sampletype=="MC_SIGNAL")) ):
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
    
    if(analysistype=="MuonSide"):
        ma.fillParticleList(decayString="K+:mychargedKaon", cut="muonID_noSVD>0.9 and nCDCHits > 20 and dr < 0.5 and abs(dz) < 2 and nPXDHits > 0",path=my_path)
        ma.fillParticleList(decayString="pi+:mychargedPion", cut="pionID>0.6 and nCDCHits > 20 and dr < 0.5 and abs(dz) < 2 and nPXDHits > 0",path=my_path)
    elif(analysistype=="ElectronSide"):
        ma.fillParticleList(decayString="K+:mychargedKaon", cut="electronID_noSVD_noTOP>0.9 and nCDCHits > 20 and dr < 0.5 and abs(dz) < 2 and nPXDHits > 0",path=my_path)
        ma.fillParticleList(decayString="pi+:mychargedPion", cut="pionID>0.6 and nCDCHits > 20 and dr < 0.5 and abs(dz) < 2 and nPXDHits > 0",path=my_path)
    else:
        ma.fillParticleList(decayString="K+:mychargedKaon", cut="kaonID>0.6 and nCDCHits > 20 and dr < 0.5 and abs(dz) < 2 and nPXDHits > 0",path=my_path)
        ma.fillParticleList(decayString="pi+:mychargedPion", cut="pionID>0.6 and nCDCHits > 20 and dr < 0.5 and abs(dz) < 2 and nPXDHits > 0",path=my_path)
    
    stdV0s.stdKshorts(path=my_path)
    ma.cutAndCopyList("K_S0:myKaonshort", "K_S0:merged", cut="[0.4876 < M < 0.5076] and significanceOfDistance > 50", path=my_path)
    
    stdPi0s.stdPi0s(listtype="eff30_May2020",path=my_path)
    ma.cutAndCopyList("pi0:myneutralPion", "pi0:eff30_May2020", cut="p > 0.4", path=my_path)
    
    stdPhotons.stdPhotons(listtype="pi0eff30_May2020", path=my_path)
    ma.cutAndCopyList("gamma:mygamma", "gamma:pi0eff30_May2020", cut="", path=my_path)

    if(analysistype=="ppbar"):
         ma.fillParticleList(decayString="p+:mychargedProton", cut="protonID>0.6 and nCDCHits > 20 and dr < 0.5 and abs(dz) < 2 and nPXDHits > 0",path=my_path)

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
    
    # remove decay modes for PID sideband
    if(analysistype=="MuonSide"):
        ma.applyCuts("Xsu:comb", "[extraInfo(decayModeID) == 0] or [extraInfo(decayModeID) == 1] or [extraInfo(decayModeID) == 3] or [extraInfo(decayModeID) == 5] or [extraInfo(decayModeID) == 7] or [extraInfo(decayModeID) == 9] or [extraInfo(decayModeID) == 11] or [extraInfo(decayModeID) == 12] or [extraInfo(decayModeID) == 13] or [extraInfo(decayModeID) == 14]", path=my_path)
        ma.applyCuts("Xsd:comb", "[extraInfo(decayModeID) == 1] or [extraInfo(decayModeID) == 3] or [extraInfo(decayModeID) == 5] or [extraInfo(decayModeID) == 7] or [extraInfo(decayModeID) == 10] or [extraInfo(decayModeID) == 12] or [extraInfo(decayModeID) == 13] or [extraInfo(decayModeID) == 14]", path=my_path)
    elif(analysistype=="ElectronSide"):
        ma.applyCuts("Xsu:comb", "[extraInfo(decayModeID) == 0] or [extraInfo(decayModeID) == 1] or [extraInfo(decayModeID) == 3] or [extraInfo(decayModeID) == 5] or [extraInfo(decayModeID) == 7] or [extraInfo(decayModeID) == 9] or [extraInfo(decayModeID) == 11] or [extraInfo(decayModeID) == 12] or [extraInfo(decayModeID) == 13] or [extraInfo(decayModeID) == 14]", path=my_path)
        ma.applyCuts("Xsd:comb", "[extraInfo(decayModeID) == 1] or [extraInfo(decayModeID) == 3] or [extraInfo(decayModeID) == 5] or [extraInfo(decayModeID) == 7] or [extraInfo(decayModeID) == 10] or [extraInfo(decayModeID) == 12] or [extraInfo(decayModeID) == 13] or [extraInfo(decayModeID) == 14]", path=my_path)
    
    # Ntuple for preselection of Xs
    if(analysistype=="nominal"):
        if(args.sampletype=="MC_SIGNAL"):
            ma.cutAndCopyList('Xsu:random', 'Xsu:comb','eventRandom > -1.0', path=my_path)
            ma.cutAndCopyList('Xsd:random', 'Xsd:comb','eventRandom > -1.0', path=my_path)
        else:
            ma.cutAndCopyList('Xsu:random', 'Xsu:comb','eventRandom > 0.95', path=my_path)
            ma.cutAndCopyList('Xsd:random', 'Xsd:comb','eventRandom > 0.95', path=my_path)
        ma.variablesToNtuple(decayString="Xsu:random",variables=["M"],filename=output_file,treename="Xsu_preselection",path=my_path)
        ma.variablesToNtuple(decayString="Xsd:random",variables=["M"],filename=output_file,treename="Xsd_preselection",path=my_path)
    
    # Cut on Mxs (preselection)
    ma.applyCuts("Xsu:comb", "[M < 2.0]", path=my_path)
    ma.applyCuts("Xsd:comb", "[M < 2.0]", path=my_path)
    
    # B->X Jpsi reconstruction
    if(analysistype=="Jpsi"):
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
        ma.variablesToNtuple(decayString="J/psi:temp_for_plot",variables=["dM", "M"],filename=JpsiName,treename="Jpsiforplot",path=my_path)
        ma.variablesToNtuple(decayString="J/psi:temp_BCS",variables=["dM", "M", "nParticlesInList(mu+:fromJpsimuchargedMuon)", "nParticlesInList(J/psi:temp_daughterUpsilon)", "nParticlesInList(J/psi:temp)"],filename=JpsiName,treename="JpsiBCS",path=my_path)
        ma.variablesToNtuple(decayString="Upsilon(4S):temp_withoutneutrino",variables=["extraInfo(decayModeID)", "daughter(1, daughter(1, dM))", "daughter(1, Mbc)", "daughter(1, deltaE)"],filename=JpsiName,treename="Upsilon_Jpsi",path=my_path)
    
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
    if(analysistype=="Jpsi"):
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
    if(analysistype=="Jpsi"):
        loose_track_selection += " and [sourceObjectIsInList(mu+:fromJpsimuchargedMuon)==0]"
    looseMask = ("looseMask", loose_track_selection, cluster_selection_v200)
    
    # Mbc side band
    if(analysistype=="MbcSide"):
        ma.applyCuts("B+:feiHadronic","Mbc < 5.26",path=my_path)
        ma.applyCuts("B0:feiHadronic","Mbc < 5.26",path=my_path)
    
    # tag side
    ma.buildRestOfEvent("B+:feiHadronic",path=my_path)
    ma.appendROEMasks("B+:feiHadronic",[cleanMask],path=my_path)
    ma.buildContinuumSuppression(list_name="B+:feiHadronic",roe_mask = "cleanMask", path=my_path)
    ma.buildRestOfEvent("B0:feiHadronic",path=my_path)
    ma.appendROEMasks("B0:feiHadronic",[cleanMask],path=my_path)
    ma.buildContinuumSuppression(list_name="B0:feiHadronic",roe_mask = "cleanMask", path=my_path)
    
    # signal side
    if(analysistype=="ppbar"):
        ma.reconstructDecay("B+:sig -> Xsu:comb p+:mychargedProton anti-p-:mychargedProton", cut="", dmID=0, path = my_path)
    elif(analysistype=="KSKS"):
        ma.reconstructDecay("B+:sig -> Xsu:comb K_S0:myKaonshort K_S0:myKaonshort", cut="", dmID=0, path = my_path)
    else:
        ma.reconstructDecay("B+:sig -> Xsu:comb ?nu", cut="", dmID=0, path = my_path)
    ma.buildRestOfEvent("B+:sig",path=my_path)
    ma.appendROEMasks("B+:sig",[cleanMask],path=my_path)
    ma.buildContinuumSuppression(list_name="B+:sig",roe_mask = "cleanMask", path=my_path)
    if(analysistype=="ppbar"):
        ma.reconstructDecay("B0:sig -> Xsd:comb p+:mychargedProton anti-p-:mychargedProton", cut="", dmID=0, path = my_path)
    elif(analysistype=="KSKS"):
        ma.reconstructDecay("B0:sig -> Xsd:comb K_S0:myKaonshort K_S0:myKaonshort", cut="", dmID=0, path = my_path)
    else:
        ma.reconstructDecay("B0:sig -> Xsd:comb ?nu", cut="", dmID=0, path = my_path)
    ma.buildRestOfEvent("B0:sig",path=my_path)
    ma.appendROEMasks("B0:sig",[cleanMask],path=my_path)
    ma.buildContinuumSuppression(list_name="B0:sig",roe_mask = "cleanMask", path=my_path)
    
    if(analysistype=="chargeSide"):
        # Upsilon(4S) -> Btag Xs wrong charge side band
        ma.reconstructDecay("Upsilon(4S):withoutneutrino_charged_opposite_cp -> B+:feiHadronic anti-B0:sig",cut ="", dmID = 0, path=my_path, allowChargeViolation=True)
        ma.reconstructDecay("Upsilon(4S):withoutneutrino_charged_same_cp -> B+:feiHadronic B0:sig",cut ="", dmID = 0, path=my_path, allowChargeViolation=True)
        ma.reconstructDecay("Upsilon(4S):withoutneutrino_neutral_opposite_cp -> B0:feiHadronic B-:sig",cut ="", dmID = 1, path=my_path, allowChargeViolation=True)
        ma.reconstructDecay("Upsilon(4S):withoutneutrino_neutral_same_cp -> B0:feiHadronic B+:sig",cut ="", dmID = 1, path=my_path, allowChargeViolation=True)
        ma.copyLists(outputListName="Upsilon(4S):withoutneutrino", inputListNames=["Upsilon(4S):withoutneutrino_charged_opposite_cp", "Upsilon(4S):withoutneutrino_charged_same_cp", "Upsilon(4S):withoutneutrino_neutral_opposite_cp", "Upsilon(4S):withoutneutrino_neutral_same_cp"], path=my_path)
        ma.buildRestOfEvent("Upsilon(4S):withoutneutrino", inputParticlelists=["K_S0:myKaonshort", "pi0:myneutralPion", "gamma:mygamma"], path=my_path)
    else:
        # Upsilon(4S) -> Btag Xs
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
    if(analysistype=="nominal"):
        Btag_preselection_vars = vu.create_aliases(list_of_variables = ["extraInfo(SignalProbability)"], wrapper = "daughter(0,{variable})",prefix="Btag")
        Bsig_preselection_vars = vu.create_aliases(list_of_variables = [ "daughter(0, extraInfo(decayModeID))"], wrapper = "daughter(1,{variable})", prefix="Bsig")
        U_preselection_vars = ["extraInfo(decayModeID)", "nROE_Tracks(looseMask)", "nROE_ParticlesInList(pi0:myneutralPion)", "nROE_Tracks(cleanMask)"]
        if(args.sampletype=="MC_SIGNAL"):
            ma.cutAndCopyList('Upsilon(4S):random', 'Upsilon(4S):withoutneutrino','eventRandom > -1.0', path=my_path)
        else:
            ma.cutAndCopyList('Upsilon(4S):random', 'Upsilon(4S):withoutneutrino','eventRandom > 0.95', path=my_path)
        ma.variablesToNtuple(decayString="Upsilon(4S):random",variables=Btag_preselection_vars,filename=output_file,treename="Btag_preselection",path=my_path)
        ma.variablesToNtuple(decayString="Upsilon(4S):random",variables=Bsig_preselection_vars,filename=output_file,treename="Bsig_preselection",path=my_path)
        ma.variablesToNtuple(decayString="Upsilon(4S):random",variables=U_preselection_vars,filename=output_file,treename="Upsilon_preselection",path=my_path)
    
    # Least Cuts (preselection)
    ma.applyCuts("Upsilon(4S):withoutneutrino", "[nROE_Tracks(cleanMask) < 0.5] and [nROE_ParticlesInList(pi0:myneutralPion) < 0.5]", path=my_path)
    
    # systematic works only for MC
    if((args.sampletype=="MC_qq") or (args.sampletype=="MC_BB") or (args.sampletype=="MC_SIGNAL")):
    
        # systematic for B->[Xc -> KL0 X] anything
        ma.findMCDecay("D+:DecayIntoKL0", "D+ -> K_L0 ... ?nu ?gamma", path=my_path)
        ma.findMCDecay("D0:DecayIntoKL0", "D0 -> K_L0 ... ?nu ?gamma", path=my_path)
    
        # systematic for B->K(*) nnbar
        ma.fillParticleListFromMC('n0:B2Knnbar', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
    
        ma.findMCDecay("B+:Knn", "B+ =direct=> K+ n0 anti-n0", appendAllDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.rankByLowest(particleList="B+:Knn", variable="random",allowMultiRank=True,outputVariable="random_rank",path=my_path)
        ma.cutAndCopyList("B+:Knn_BCS", "B+:Knn", "extraInfo(random_rank) == 1",path=my_path)
    
        ma.findMCDecay("B+:Kstarnn", "B+ =direct=> K*+ n0 anti-n0", appendAllDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.rankByLowest(particleList="B+:Kstarnn", variable="random",allowMultiRank=True,outputVariable="random_rank",path=my_path)
        ma.cutAndCopyList("B+:Kstarnn_BCS", "B+:Kstarnn", "extraInfo(random_rank) == 1",path=my_path)
    
        ma.findMCDecay("B0:K0nn", "B0 =direct=> K0 n0 anti-n0", appendAllDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.rankByLowest(particleList="B0:K0nn", variable="random",allowMultiRank=True,outputVariable="random_rank",path=my_path)
        ma.cutAndCopyList("B0:K0nn_BCS", "B0:K0nn", "extraInfo(random_rank) == 1",path=my_path)
    
        ma.findMCDecay("B0:Kstar0nn", "B0 =direct=> K*0 n0 anti-n0", appendAllDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.rankByLowest(particleList="B0:Kstar0nn", variable="random",allowMultiRank=True,outputVariable="random_rank",path=my_path)
        ma.cutAndCopyList("B0:Kstar0nn_BCS", "B0:Kstar0nn", "extraInfo(random_rank) == 1",path=my_path)
    
        ma.cutAndCopyList("n0:Knn", "n0:B2Knnbar", cut="isDescendantOfList(B+:Knn_BCS, 1)", path=my_path)
        ma.cutAndCopyList("n0:Kstarnn", "n0:B2Knnbar", cut="isDescendantOfList(B+:Kstarnn_BCS, 1)", path=my_path)
        ma.cutAndCopyList("n0:K0nn", "n0:B2Knnbar", cut="isDescendantOfList(B0:K0nn_BCS, 1)", path=my_path)
        ma.cutAndCopyList("n0:Kstar0nn", "n0:B2Knnbar", cut="isDescendantOfList(B0:Kstar0nn_BCS, 1)", path=my_path)

        # systematic for B->X nnbar
        ma.findMCDecay("B+:Xnn", "B+ -> n0 anti-n0 ... ?gamma", appendAllDaughters=False, path=my_path)
        ma.findMCDecay("B0:Xnn", "B0 -> n0 anti-n0 ... ?gamma", appendAllDaughters=False, path=my_path)
        ma.fillParticleListFromMC("n0:Xnn", cut = "[isMCDescendantOfList(B+:Xnn, -1) == 1] or [isMCDescendantOfList(B0:Xnn, -1) == 1]", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)

        # systematic for B+ -> K+ KL0 KL0
        ma.findMCDecay("B+:KpKLKL_all", "B+ -> K+ K_L0 K_L0 ?gamma", appendAllDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.findMCDecay("B+:KpKLKL_NR", "B+ =direct=> K+ [K0 =direct=> K_L0] [anti-K0 =direct=> K_L0] ?gamma", appendAllDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.rankByLowest(particleList="B+:KpKLKL_NR", variable="random",allowMultiRank=True,outputVariable="random_rank",path=my_path)
        ma.applyCuts("B+:KpKLKL_NR", "extraInfo(random_rank) == 1", path=my_path)
    
        # systematic for B0 -> KS0 KL0 KL0
        ma.findMCDecay("B0:KSKLKL_all", "B0 -> K_S0 K_L0 K_L0 ?gamma", appendAllDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.findMCDecay("B0:KSKLKL_NR_1", "B0 =direct=> [K0 =direct=> K_S0] [K0 =direct=> K_L0] [anti-K0 =direct=> K_L0] ?gamma", appendAllDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.findMCDecay("B0:KSKLKL_NR_2", "B0 =direct=> [anti-K0 =direct=> K_S0] [K0 =direct=> K_L0] [K0 =direct=> K_L0] ?gamma", appendAllDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.copyLists(outputListName="B0:KSKLKL_NR", inputListNames=["B0:KSKLKL_NR_1", "B0:KSKLKL_NR_2"], path=my_path)
        ma.rankByLowest(particleList="B0:KSKLKL_NR", variable="random",allowMultiRank=True,outputVariable="random_rank",path=my_path)
        ma.applyCuts("B0:KSKLKL_NR", "extraInfo(random_rank) == 1", path=my_path)

        # systematic for B -> X KL0 KL0
        ma.findMCDecay("B+:XKLKL", "B+ -> K_L0 K_L0 ... ?gamma", appendAllDaughters=False, path=my_path)
        ma.findMCDecay("B0:XKLKL", "B0 -> K_L0 K_L0 ... ?gamma", appendAllDaughters=False, path=my_path)
        ma.fillParticleListFromMC("K_L0:XKLKL", cut = "[isMCDescendantOfList(B+:XKLKL, -1) == 1] or [isMCDescendantOfList(B0:XKLKL, -1) == 1]", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.fillParticleListFromMC("K_S0:PrimaryMC_XKLKL", cut = "mcPrimary", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.reconstructMCDecay('B+:XKLKLKS -> K_L0:XKLKL K_L0:XKLKL K_S0:PrimaryMC_XKLKL ... ?gamma', cut='', path=my_path)
        ma.reconstructMCDecay('B0:XKLKLKS -> K_L0:XKLKL K_L0:XKLKL K_S0:PrimaryMC_XKLKL ... ?gamma', cut='', path=my_path)
        ma.fillParticleListFromMC("K_S0:XKLKL", cut = "[isMCDescendantOfList(B+:XKLKLKS, -1) == 1] or [isMCDescendantOfList(B0:XKLKLKS, -1) == 1]", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.rankByHighest(particleList="K_L0:XKLKL", variable="random",allowMultiRank=False,outputVariable="KL_XKLKL_rank",path=my_path)
        ma.cutAndCopyList("K_L0:XKLKL_1st", "K_L0:XKLKL", "extraInfo(KL_XKLKL_rank) == 1",path=my_path)
        ma.cutAndCopyList("K_L0:XKLKL_2nd", "K_L0:XKLKL", "extraInfo(KL_XKLKL_rank) == 2",path=my_path)

    # lepton information
    ma.fillParticleList(decayString="e+:ElectronFBDT", cut=track_selection + " and electronID>0.6",path=my_path)
    ma.fillParticleList(decayString="e+:ElectronFBDT_loose", cut=track_selection + " and electronID>0.1",path=my_path)
    ma.fillParticleList(decayString="e+:ElectronFBDT_tight", cut=track_selection + " and electronID>0.9",path=my_path)
    ma.fillParticleList(decayString="mu+:MuonFBDT", cut=track_selection + " and muonID>0.6", path=my_path)
    ma.fillParticleList(decayString="mu+:MuonFBDT_loose", cut=track_selection + " and muonID>0.1", path=my_path)
    ma.fillParticleList(decayString="mu+:MuonFBDT_tight", cut=track_selection + " and muonID>0.9", path=my_path)
    
    # Apply cut for D veto (save computer resource)
    ma.applyCuts("Xsu:comb", "isDescendantOfList(Upsilon(4S):withoutneutrino, 2) == 1", path=my_path)
    ma.applyCuts("Xsd:comb", "isDescendantOfList(Upsilon(4S):withoutneutrino, 2) == 1", path=my_path)

    # add variables for B->X KS0 KS0 analysis
    if(analysistype=="KSKS"):
        ma.applyCuts("B+:sig", "isDescendantOfList(Upsilon(4S):withoutneutrino, 1) == 1", path=my_path)
        ma.applyCuts("B0:sig", "isDescendantOfList(Upsilon(4S):withoutneutrino, 1) == 1", path=my_path)

        # for B+
        KSKS_path_one = basf2.Path()
        KSKS_deadEndPath_one = basf2.Path()
        ma.signalSideParticleFilter('B+:sig', '', KSKS_path_one, KSKS_deadEndPath_one)
        ma.fillSignalSideParticleList('B+:sig_KSKS', '^B+:sig', KSKS_path_one)
        ma.cutAndCopyList("K_S0:sig_KSKS_p", "K_S0:myKaonshort", cut="isDescendantOfList(B+:sig_KSKS, -1) == 1", path=KSKS_path_one)
        ma.rankByHighest(particleList="K_S0:sig_KSKS_p", variable="random",allowMultiRank=False,outputVariable="KS_KSKS_p_rank",path=KSKS_path_one)
        ma.cutAndCopyList("K_S0:sig_KSKS_p_1st", "K_S0:sig_KSKS_p", "extraInfo(KS_KSKS_p_rank) == 1",path=KSKS_path_one)
        ma.cutAndCopyList("K_S0:sig_KSKS_p_2nd", "K_S0:sig_KSKS_p", "extraInfo(KS_KSKS_p_rank) == 2",path=KSKS_path_one)
        ma.variableToSignalSideExtraInfo('B+:sig_KSKS', {"nParticlesInList(K_S0:sig_KSKS_p)": "N_KS_KSKS"}, path=KSKS_path_one)
        ma.variableToSignalSideExtraInfo('B+:sig_KSKS', {"averageValueInList(K_S0:sig_KSKS_p,E)": "avg_KS_KSKS_E"}, path=KSKS_path_one)
        ma.variableToSignalSideExtraInfo('B+:sig_KSKS', {"averageValueInList(K_S0:sig_KSKS_p,px)": "avg_KS_KSKS_px"}, path=KSKS_path_one)
        ma.variableToSignalSideExtraInfo('B+:sig_KSKS', {"averageValueInList(K_S0:sig_KSKS_p,py)": "avg_KS_KSKS_py"}, path=KSKS_path_one)
        ma.variableToSignalSideExtraInfo('B+:sig_KSKS', {"averageValueInList(K_S0:sig_KSKS_p,pz)": "avg_KS_KSKS_pz"}, path=KSKS_path_one)
        ma.variableToSignalSideExtraInfo('K_S0:sig_KSKS_p_1st', {"E": "KS_KSKS_1st_E"}, path=KSKS_path_one)
        ma.variableToSignalSideExtraInfo('K_S0:sig_KSKS_p_1st', {"px": "KS_KSKS_1st_px"}, path=KSKS_path_one)
        ma.variableToSignalSideExtraInfo('K_S0:sig_KSKS_p_1st', {"py": "KS_KSKS_1st_py"}, path=KSKS_path_one)
        ma.variableToSignalSideExtraInfo('K_S0:sig_KSKS_p_1st', {"pz": "KS_KSKS_1st_pz"}, path=KSKS_path_one)
        ma.variableToSignalSideExtraInfo('K_S0:sig_KSKS_p_2nd', {"E": "KS_KSKS_2nd_E"}, path=KSKS_path_one)
        ma.variableToSignalSideExtraInfo('K_S0:sig_KSKS_p_2nd', {"px": "KS_KSKS_2nd_px"}, path=KSKS_path_one)
        ma.variableToSignalSideExtraInfo('K_S0:sig_KSKS_p_2nd', {"py": "KS_KSKS_2nd_py"}, path=KSKS_path_one)
        ma.variableToSignalSideExtraInfo('K_S0:sig_KSKS_p_2nd', {"pz": "KS_KSKS_2nd_pz"}, path=KSKS_path_one)
        my_path.for_each('RestOfEvent', 'RestOfEvents', KSKS_path_one)

        # for B0
        KSKS_path_two = basf2.Path()
        KSKS_deadEndPath_two = basf2.Path()
        ma.signalSideParticleFilter('B0:sig', '', KSKS_path_two, KSKS_deadEndPath_two)
        ma.fillSignalSideParticleList('B0:sig_KSKS', '^B0:sig', KSKS_path_two)
        ma.cutAndCopyList("K_S0:sig_KSKS_n", "K_S0:myKaonshort", cut="isDescendantOfList(B0:sig_KSKS, -1) == 1", path=KSKS_path_two)
        ma.rankByHighest(particleList="K_S0:sig_KSKS_n", variable="random",allowMultiRank=False,outputVariable="KS_KSKS_n_rank",path=KSKS_path_two)
        ma.cutAndCopyList("K_S0:sig_KSKS_n_1st", "K_S0:sig_KSKS_n", "extraInfo(KS_KSKS_n_rank) == 1",path=KSKS_path_two)
        ma.cutAndCopyList("K_S0:sig_KSKS_n_2nd", "K_S0:sig_KSKS_n", "extraInfo(KS_KSKS_n_rank) == 2",path=KSKS_path_two)
        ma.variableToSignalSideExtraInfo('B0:sig_KSKS', {"nParticlesInList(K_S0:sig_KSKS_n)": "N_KS_KSKS"}, path=KSKS_path_two)
        ma.variableToSignalSideExtraInfo('B0:sig_KSKS', {"averageValueInList(K_S0:sig_KSKS_n,E)": "avg_KS_KSKS_E"}, path=KSKS_path_two)
        ma.variableToSignalSideExtraInfo('B0:sig_KSKS', {"averageValueInList(K_S0:sig_KSKS_n,px)": "avg_KS_KSKS_px"}, path=KSKS_path_two)
        ma.variableToSignalSideExtraInfo('B0:sig_KSKS', {"averageValueInList(K_S0:sig_KSKS_n,py)": "avg_KS_KSKS_py"}, path=KSKS_path_two)
        ma.variableToSignalSideExtraInfo('B0:sig_KSKS', {"averageValueInList(K_S0:sig_KSKS_n,pz)": "avg_KS_KSKS_pz"}, path=KSKS_path_two)
        ma.variableToSignalSideExtraInfo('K_S0:sig_KSKS_n_1st', {"E": "KS_KSKS_1st_E"}, path=KSKS_path_two)
        ma.variableToSignalSideExtraInfo('K_S0:sig_KSKS_n_1st', {"px": "KS_KSKS_1st_px"}, path=KSKS_path_two)
        ma.variableToSignalSideExtraInfo('K_S0:sig_KSKS_n_1st', {"py": "KS_KSKS_1st_py"}, path=KSKS_path_two)
        ma.variableToSignalSideExtraInfo('K_S0:sig_KSKS_n_1st', {"pz": "KS_KSKS_1st_pz"}, path=KSKS_path_two)
        ma.variableToSignalSideExtraInfo('K_S0:sig_KSKS_n_2nd', {"E": "KS_KSKS_2nd_E"}, path=KSKS_path_two)
        ma.variableToSignalSideExtraInfo('K_S0:sig_KSKS_n_2nd', {"px": "KS_KSKS_2nd_px"}, path=KSKS_path_two)
        ma.variableToSignalSideExtraInfo('K_S0:sig_KSKS_n_2nd', {"py": "KS_KSKS_2nd_py"}, path=KSKS_path_two)
        ma.variableToSignalSideExtraInfo('K_S0:sig_KSKS_n_2nd', {"pz": "KS_KSKS_2nd_pz"}, path=KSKS_path_two)
        my_path.for_each('RestOfEvent', 'RestOfEvents', KSKS_path_two)

    # D veto
    ma.buildRestOfEvent("Xsu:comb", path=my_path)
    roe_path_one = basf2.Path()
    deadEndPath_one = basf2.Path()
    ma.signalSideParticleFilter('Xsu:comb', '', roe_path_one, deadEndPath_one)
    ma.fillSignalSideParticleList('Xsu:sig', '^Xsu:comb', roe_path_one)
    ma.cutAndCopyList("K+:sig_su", "K+:mychargedKaon", cut="isDescendantOfList(Xsu:sig, 1) == 1", path=roe_path_one)
    ma.cutAndCopyList("K_S0:sig_su", "K_S0:myKaonshort", cut="isDescendantOfList(Xsu:sig, 1) == 1", path=roe_path_one)
    ma.cutAndCopyList("pi+:sig_su", "pi+:mychargedPion", cut="isDescendantOfList(Xsu:sig, 1) == 1", path=roe_path_one)
    ma.cutAndCopyList("pi0:sig_su", "pi0:myneutralPion", cut="isDescendantOfList(Xsu:sig, 1) == 1", path=roe_path_one)
    ma.matchMCTruth("K+:sig_su", path=roe_path_one)
    ma.matchMCTruth("pi+:sig_su", path=roe_path_one)
    
    # uncertainty of slow track
    ma.cutAndCopyList("K+:slow1_su", "K+:mychargedKaon", cut="isDescendantOfList(Xsu:sig, 1) == 1 and 0.05 < p < 0.12", path=roe_path_one)
    ma.cutAndCopyList("K+:slow2_su", "K+:mychargedKaon", cut="isDescendantOfList(Xsu:sig, 1) == 1 and 0.12 < p < 0.16", path=roe_path_one)
    ma.cutAndCopyList("K+:slow3_su", "K+:mychargedKaon", cut="isDescendantOfList(Xsu:sig, 1) == 1 and 0.16 < p < 0.20", path=roe_path_one)
    ma.cutAndCopyList("pi+:slow1_su", "pi+:mychargedPion", cut="isDescendantOfList(Xsu:sig, 1) == 1 and 0.05 < p < 0.12", path=roe_path_one)
    ma.cutAndCopyList("pi+:slow2_su", "pi+:mychargedPion", cut="isDescendantOfList(Xsu:sig, 1) == 1 and 0.12 < p < 0.16", path=roe_path_one)
    ma.cutAndCopyList("pi+:slow3_su", "pi+:mychargedPion", cut="isDescendantOfList(Xsu:sig, 1) == 1 and 0.16 < p < 0.20", path=roe_path_one)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'nParticlesInList(K+:slow1_su)': 'nKslow1'}, path=roe_path_one)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'nParticlesInList(K+:slow2_su)': 'nKslow2'}, path=roe_path_one)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'nParticlesInList(K+:slow3_su)': 'nKslow3'}, path=roe_path_one)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'nParticlesInList(pi+:slow1_su)': 'nPislow1'}, path=roe_path_one)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'nParticlesInList(pi+:slow2_su)': 'nPislow2'}, path=roe_path_one)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'nParticlesInList(pi+:slow3_su)': 'nPislow3'}, path=roe_path_one)
    
    # uncertainty of PID
    pbins = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
    cosThetabins = [-0.866, -0.682, -0.4226, -0.1045, 0.225, 0.5, 0.766, 0.8829, 0.9563]
    index = 0
    for i_p in range(len(pbins) - 1):
        for i_cosTheta in range(len(cosThetabins) - 1):
            pcut_string = "[" + str(pbins[i_p]) + " < p < " + str(pbins[i_p+1]) + "]"
            cosThetacut_string = "[" + str(cosThetabins[i_cosTheta]) + " < cosTheta < " + str(cosThetabins[i_cosTheta+1]) + "]"
            ma.cutAndCopyList("K+:truePIDbin" + str(index), "K+:sig_su", cut="[abs(mcPDG) == 321] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + pcut_string + " and " + cosThetacut_string, path=roe_path_one)
            ma.cutAndCopyList("K+:misPIDbin" + str(index), "K+:sig_su", cut="[abs(mcPDG) == 211] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + pcut_string + " and " + cosThetacut_string, path=roe_path_one)
            ma.cutAndCopyList("pi+:truePIDbin" + str(index), "pi+:sig_su", cut="[abs(mcPDG) == 211] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + pcut_string + " and " + cosThetacut_string, path=roe_path_one)
            ma.cutAndCopyList("pi+:misPIDbin" + str(index), "pi+:sig_su", cut="[abs(mcPDG) == 321] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + pcut_string + " and " + cosThetacut_string, path=roe_path_one)
            index = index + 1
    pcut_string = "[ p < " + str(pbins[0]) + "]" + " or " + "[ p > " + str(pbins[-1]) + "]"
    cosThetacut_string = "[ cosTheta < " + str(cosThetabins[0]) + "]" + " or " + "[ cosTheta > " + str(cosThetabins[-1]) + "]"
    ma.cutAndCopyList("K+:truePIDbin" + str(index), "K+:sig_su", cut="[abs(mcPDG) == 321] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + cosThetacut_string + "]", path=roe_path_one)
    ma.cutAndCopyList("K+:misPIDbin" + str(index), "K+:sig_su", cut="[abs(mcPDG) == 211] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + cosThetacut_string + "]", path=roe_path_one)
    ma.cutAndCopyList("pi+:truePIDbin" + str(index), "pi+:sig_su", cut="[abs(mcPDG) == 211] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + cosThetacut_string + "]", path=roe_path_one)
    ma.cutAndCopyList("pi+:misPIDbin" + str(index), "pi+:sig_su", cut="[abs(mcPDG) == 321] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + cosThetacut_string + "]", path=roe_path_one)
    index = 0
    for i_p in range(len(pbins) - 1):
        for i_cosTheta in range(len(cosThetabins) - 1):
            ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(K+:truePIDbin" + str(index) + ")": "nKtruebin"+str(index)}, path=roe_path_one)
            ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(K+:misPIDbin" + str(index) + ")": "nKmisbin"+str(index)}, path=roe_path_one)
            ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(pi+:truePIDbin" + str(index) + ")": "npitruebin"+str(index)}, path=roe_path_one)
            ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(pi+:misPIDbin" + str(index) + ")": "npimisbin"+str(index)}, path=roe_path_one)
            index = index + 1
    ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(K+:truePIDbin" + str(index) + ")": "nKtruebin"+str(index)}, path=roe_path_one)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(K+:misPIDbin" + str(index) + ")": "nKmisbin"+str(index)}, path=roe_path_one)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(pi+:truePIDbin" + str(index) + ")": "npitruebin"+str(index)}, path=roe_path_one)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(pi+:misPIDbin" + str(index) + ")": "npimisbin"+str(index)}, path=roe_path_one)
    
    pbins_e = [0.4, 0.5, 1.0, 1.5, 2, 2.5, 3]
    thetabins_e = [0.22, 0.56, 1.13, 1.57, 1.88, 2.23, 2.71]
    index = 0
    for i_p in range(len(pbins_e) - 1):
        for i_theta in range(len(thetabins_e) - 1):
            pcut_string = "[" + str(pbins_e[i_p]) + " < p < " + str(pbins_e[i_p+1]) + "]"
            thetacut_string = "[" + str(thetabins_e[i_theta]) + " < theta < " + str(thetabins_e[i_theta+1]) + "]"
            ma.cutAndCopyList("K+:fakeEbin_n" + str(index), "K+:sig_su", cut="[mcPDG == 11] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + pcut_string + " and " + thetacut_string, path=roe_path_one)
            ma.cutAndCopyList("K+:fakeEbin_p" + str(index), "K+:sig_su", cut="[mcPDG == -11] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + pcut_string + " and " + thetacut_string, path=roe_path_one)
            ma.cutAndCopyList("pi+:fakeEbin_n" + str(index), "pi+:sig_su", cut="[mcPDG == 11] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + pcut_string + " and " + thetacut_string, path=roe_path_one)
            ma.cutAndCopyList("pi+:fakeEbin_p" + str(index), "pi+:sig_su", cut="[mcPDG == -11] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + pcut_string + " and " + thetacut_string, path=roe_path_one)
            index = index + 1
    pcut_string = "[ p < " + str(pbins_e[0]) + "]" + " or " + "[ p > " + str(pbins_e[-1]) + "]"
    thetacut_string = "[ theta < " + str(thetabins_e[0]) + "]" + " or " + "[ theta > " + str(thetabins_e[-1]) + "]"
    ma.cutAndCopyList("K+:fakeEbin_n" + str(index), "K+:sig_su", cut="[mcPDG == 11] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + thetacut_string + "]", path=roe_path_one)
    ma.cutAndCopyList("K+:fakeEbin_p" + str(index), "K+:sig_su", cut="[mcPDG == -11] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + thetacut_string + "]", path=roe_path_one)
    ma.cutAndCopyList("pi+:fakeEbin_n" + str(index), "pi+:sig_su", cut="[mcPDG == 11] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + thetacut_string + "]", path=roe_path_one)
    ma.cutAndCopyList("pi+:fakeEbin_p" + str(index), "pi+:sig_su", cut="[mcPDG == -11] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + thetacut_string + "]", path=roe_path_one)
    index = 0
    for i_p in range(len(pbins_e) - 1):
        for i_theta in range(len(thetabins_e) - 1):
            ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(K+:fakeEbin_n" + str(index) + ")": "nKfakeEbin_n"+str(index)}, path=roe_path_one)
            ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(K+:fakeEbin_p" + str(index) + ")": "nKfakeEbin_p"+str(index)}, path=roe_path_one)
            ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(pi+:fakeEbin_n" + str(index) + ")": "npifakeEbin_n"+str(index)}, path=roe_path_one)
            ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(pi+:fakeEbin_p" + str(index) + ")": "npifakeEbin_p"+str(index)}, path=roe_path_one)
            index = index + 1
    ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(K+:fakeEbin_n" + str(index) + ")": "nKfakeEbin_n"+str(index)}, path=roe_path_one)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(K+:fakeEbin_p" + str(index) + ")": "nKfakeEbin_p"+str(index)}, path=roe_path_one)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(pi+:fakeEbin_n" + str(index) + ")": "npifakeEbin_n"+str(index)}, path=roe_path_one)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(pi+:fakeEbin_p" + str(index) + ")": "npifakeEbin_p"+str(index)}, path=roe_path_one)
    
    pbins_mu = [0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5]
    thetabins_mu = [0.4, 0.64, 0.82, 1.16, 1.46, 1.78, 2.13, 2.22, 2.6]
    index = 0
    for i_p in range(len(pbins_mu) - 1):
        for i_theta in range(len(thetabins_mu) - 1):
            pcut_string = "[" + str(pbins_mu[i_p]) + " < p < " + str(pbins_mu[i_p+1]) + "]"
            thetacut_string = "[" + str(thetabins_mu[i_theta]) + " < theta < " + str(thetabins_mu[i_theta+1]) + "]"
            ma.cutAndCopyList("K+:fakeMUbin_n" + str(index), "K+:sig_su", cut="[mcPDG == 13] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + pcut_string + " and " + thetacut_string, path=roe_path_one)
            ma.cutAndCopyList("K+:fakeMUbin_p" + str(index), "K+:sig_su", cut="[mcPDG == -13] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + pcut_string + " and " + thetacut_string, path=roe_path_one)
            ma.cutAndCopyList("pi+:fakeMUbin_n" + str(index), "pi+:sig_su", cut="[mcPDG == 13] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + pcut_string + " and " + thetacut_string, path=roe_path_one)
            ma.cutAndCopyList("pi+:fakeMUbin_p" + str(index), "pi+:sig_su", cut="[mcPDG == -13] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + pcut_string + " and " + thetacut_string, path=roe_path_one)
            index = index + 1
    pcut_string = "[ p < " + str(pbins_mu[0]) + "]" + " or " + "[ p > " + str(pbins_mu[-1]) + "]"
    thetacut_string = "[ theta < " + str(thetabins_mu[0]) + "]" + " or " + "[ theta > " + str(thetabins_mu[-1]) + "]"
    ma.cutAndCopyList("K+:fakeMUbin_n" + str(index), "K+:sig_su", cut="[mcPDG == 13] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + thetacut_string + "]", path=roe_path_one)
    ma.cutAndCopyList("K+:fakeMUbin_p" + str(index), "K+:sig_su", cut="[mcPDG == -13] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + thetacut_string + "]", path=roe_path_one)
    ma.cutAndCopyList("pi+:fakeMUbin_n" + str(index), "pi+:sig_su", cut="[mcPDG == 13] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + thetacut_string + "]", path=roe_path_one)
    ma.cutAndCopyList("pi+:fakeMUbin_p" + str(index), "pi+:sig_su", cut="[mcPDG == -13] and [isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + thetacut_string + "]", path=roe_path_one)
    index = 0
    for i_p in range(len(pbins_mu) - 1):
        for i_theta in range(len(thetabins_mu) - 1):
            ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(K+:fakeMUbin_n" + str(index) + ")": "nKfakeMUbin_n"+str(index)}, path=roe_path_one)
            ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(K+:fakeMUbin_p" + str(index) + ")": "nKfakeMUbin_p"+str(index)}, path=roe_path_one)
            ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(pi+:fakeMUbin_n" + str(index) + ")": "npifakeMUbin_n"+str(index)}, path=roe_path_one)
            ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(pi+:fakeMUbin_p" + str(index) + ")": "npifakeMUbin_p"+str(index)}, path=roe_path_one)
            index = index + 1
    ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(K+:fakeMUbin_n" + str(index) + ")": "nKfakeMUbin_n"+str(index)}, path=roe_path_one)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(K+:fakeMUbin_p" + str(index) + ")": "nKfakeMUbin_p"+str(index)}, path=roe_path_one)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(pi+:fakeMUbin_n" + str(index) + ")": "npifakeMUbin_n"+str(index)}, path=roe_path_one)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(pi+:fakeMUbin_p" + str(index) + ")": "npifakeMUbin_p"+str(index)}, path=roe_path_one)
    
    # uncertainty of pi0
    pi0_pbins = [0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0]
    index = 0
    for i_p in range(len(pi0_pbins) - 1):
        pcut_string = "[" + str(pi0_pbins[i_p]) + " < p < " + str(pi0_pbins[i_p+1]) + "]"
        ma.cutAndCopyList("pi0:pbin" + str(index), "pi0:sig_su", cut="[isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + pcut_string, path=roe_path_one)
        index = index + 1
    pcut_string = "[ p < " + str(pi0_pbins[0]) + "]" + " or " + "[ p > " + str(pi0_pbins[-1]) + "]"
    ma.cutAndCopyList("pi0:pbin" + str(index), "pi0:sig_su", cut="[isDescendantOfList(Xsu:sig, 1) == 1]" + " and " + "[" + pcut_string + "]", path=roe_path_one)
    index = 0
    for i_p in range(len(pi0_pbins) - 1):
        ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(pi0:pbin" + str(index) + ")": "npi0bin"+str(index)}, path=roe_path_one)
        index = index + 1
    ma.variableToSignalSideExtraInfo('Xsu:sig', {"nParticlesInList(pi0:pbin" + str(index) + ")": "npi0bin"+str(index)}, path=roe_path_one)
    
    # D veto all
    ma.fillParticleList(decayString="pi+:simple", cut="dr < 0.5 and abs(dz) < 2 and nPXDHits > 0",path=roe_path_one)
    ma.cutAndCopyList("pi0:simple", "pi0:eff30_May2020", cut="p > 0.4", path=roe_path_one)
    if(analysistype=="Jpsi"):
        ma.applyCuts("pi+:simple", "sourceObjectIsInList(mu+:fromJpsimuchargedMuon)==0", path=roe_path_one)
    ma.reconstructDecay(decayString="D+:simpleveto1_su -> K+:sig_su pi-:simple pi+:simple", cut="M < 2.1", dmID=1, path=roe_path_one) # DCS
    ma.reconstructDecay(decayString="D+:simpleveto2_su -> K_S0:sig_su pi+:simple", cut="M < 2.1", dmID=2, path=roe_path_one)
    ma.reconstructDecay(decayString="D+:simpleveto3_su -> K-:sig_su pi+:simple pi+:simple", cut="M < 2.1", dmID=3, path=roe_path_one)
    ma.reconstructDecay(decayString="D+:simpleveto4_su -> K_S0:sig_su pi+:simple pi0:simple", cut="M < 2.1", dmID=4, path=roe_path_one) # pi0
    ma.reconstructDecay(decayString="D+:simpleveto5_su -> K-:sig_su pi+:simple pi+:simple pi0:simple", cut="M < 2.1", dmID=5, path=roe_path_one) # pi0
    ma.reconstructDecay(decayString="D+:simpleveto6_su -> K_S0:sig_su pi+:simple pi0:simple pi0:simple", cut="M < 2.1", dmID=6, path=roe_path_one) # pi0
    vertex.treeFit(list_name="D+:simpleveto1_su", conf_level=0, updateAllDaughters=False, path=roe_path_one)
    vertex.treeFit(list_name="D+:simpleveto2_su", conf_level=0, updateAllDaughters=False, path=roe_path_one)
    vertex.treeFit(list_name="D+:simpleveto3_su", conf_level=0, updateAllDaughters=False, path=roe_path_one)
    vertex.treeFit(list_name="D+:simpleveto4_su", conf_level=0, updateAllDaughters=False, path=roe_path_one)
    vertex.treeFit(list_name="D+:simpleveto5_su", conf_level=0, updateAllDaughters=False, path=roe_path_one)
    vertex.treeFit(list_name="D+:simpleveto6_su", conf_level=0, updateAllDaughters=False, path=roe_path_one)
    ma.copyLists(outputListName="D+:simplevetocomb_su", inputListNames=["D+:simpleveto1_su", "D+:simpleveto2_su", "D+:simpleveto3_su"], path=roe_path_one)
    ma.copyLists(outputListName="D+:simplevetocomb_su_noDCS", inputListNames=["D+:simpleveto2_su", "D+:simpleveto3_su"], path=roe_path_one)
    ma.copyLists(outputListName="D+:simplevetocomb_su_yespizero", inputListNames=["D+:simpleveto1_su", "D+:simpleveto2_su", "D+:simpleveto3_su", "D+:simpleveto4_su", "D+:simpleveto5_su", "D+:simpleveto6_su"], path=roe_path_one)
    
    Dcsimple_yes_su = basf2.create_path()
    Dcsimple_no_su = basf2.create_path()
    ma.rankByHighest(particleList="D+:simplevetocomb_su", variable="chiProb", outputVariable="Dsimple_rank",path=Dcsimple_yes_su)
    ma.cutAndCopyList("D+:simplevetocomb_BCS_su", "D+:simplevetocomb_su", cut="extraInfo(Dsimple_rank)==1", path=Dcsimple_yes_su)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_su", {'countInList(D+:simplevetocomb_su, chiProb>0.1)': 'nDc'}, path=Dcsimple_yes_su)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_su", {'medianValueInList(D+:simplevetocomb_su,chiProb)': 'Dc_pValue_med'}, path=Dcsimple_yes_su)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_su", {'formula(averageValueInList(D+:simplevetocomb_su, formula((chiProb-averageValueInList(D+:simplevetocomb_su,chiProb))**2))**0.5)': 'Dc_pValue_std'}, path=Dcsimple_yes_su)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_su", {'chiProb': 'Dcsimpleveto_chiProb'}, path=Dcsimple_yes_su)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_su", {'dr': 'Dcsimpleveto_dr'}, path=Dcsimple_yes_su)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_su", {'dz': 'Dcsimpleveto_dz'}, path=Dcsimple_yes_su)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_su", {'M': 'Dcsimpleveto_M'}, path=Dcsimple_yes_su)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(0)': 'nDc'}, path=Dcsimple_no_su)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'Dc_pValue_med'}, path=Dcsimple_no_su)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'Dc_pValue_std'}, path=Dcsimple_no_su)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'Dcsimpleveto_chiProb'}, path=Dcsimple_no_su)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'Dcsimpleveto_dr'}, path=Dcsimple_no_su)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'Dcsimpleveto_dz'}, path=Dcsimple_no_su)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'Dcsimpleveto_M'}, path=Dcsimple_no_su)
    nDcsimple_su_module = roe_path_one.add_module("VariableToReturnValue", variable="nParticlesInList(D+:simplevetocomb_su)")
    nDcsimple_su_module.if_value(">=1", Dcsimple_yes_su, basf2.AfterConditionPath.CONTINUE)
    nDcsimple_su_module.if_value("<1", Dcsimple_no_su, basf2.AfterConditionPath.CONTINUE)

    Dcsimple_yes_su_noDCS = basf2.create_path()
    Dcsimple_no_su_noDCS = basf2.create_path()
    ma.rankByHighest(particleList="D+:simplevetocomb_su_noDCS", variable="chiProb", outputVariable="Dsimple_rank_noDCS",path=Dcsimple_yes_su_noDCS)
    ma.cutAndCopyList("D+:simplevetocomb_BCS_su_noDCS", "D+:simplevetocomb_su_noDCS", cut="extraInfo(Dsimple_rank_noDCS)==1", path=Dcsimple_yes_su_noDCS)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_su_noDCS", {'countInList(D+:simplevetocomb_su_noDCS, chiProb>0.1)': 'nDc_noDCS'}, path=Dcsimple_yes_su_noDCS)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_su_noDCS", {'medianValueInList(D+:simplevetocomb_su_noDCS,chiProb)': 'Dc_pValue_med_noDCS'}, path=Dcsimple_yes_su_noDCS)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_su_noDCS", {'formula(averageValueInList(D+:simplevetocomb_su_noDCS, formula((chiProb-averageValueInList(D+:simplevetocomb_su_noDCS,chiProb))**2))**0.5)': 'Dc_pValue_std_noDCS'}, path=Dcsimple_yes_su_noDCS)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_su_noDCS", {'chiProb': 'Dcsimpleveto_chiProb_noDCS'}, path=Dcsimple_yes_su_noDCS)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_su_noDCS", {'dr': 'Dcsimpleveto_dr_noDCS'}, path=Dcsimple_yes_su_noDCS)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_su_noDCS", {'dz': 'Dcsimpleveto_dz_noDCS'}, path=Dcsimple_yes_su_noDCS)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_su_noDCS", {'M': 'Dcsimpleveto_M_noDCS'}, path=Dcsimple_yes_su_noDCS)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(0)': 'nDc_noDCS'}, path=Dcsimple_no_su_noDCS)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'Dc_pValue_med_noDCS'}, path=Dcsimple_no_su_noDCS)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'Dc_pValue_std_noDCS'}, path=Dcsimple_no_su_noDCS)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'Dcsimpleveto_chiProb_noDCS'}, path=Dcsimple_no_su_noDCS)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'Dcsimpleveto_dr_noDCS'}, path=Dcsimple_no_su_noDCS)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'Dcsimpleveto_dz_noDCS'}, path=Dcsimple_no_su_noDCS)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'Dcsimpleveto_M_noDCS'}, path=Dcsimple_no_su_noDCS)
    nDcsimple_su_module = roe_path_one.add_module("VariableToReturnValue", variable="nParticlesInList(D+:simplevetocomb_su_noDCS)")
    nDcsimple_su_module.if_value(">=1", Dcsimple_yes_su_noDCS, basf2.AfterConditionPath.CONTINUE)
    nDcsimple_su_module.if_value("<1", Dcsimple_no_su_noDCS, basf2.AfterConditionPath.CONTINUE)

    Dcsimple_yes_su_yespizero = basf2.create_path()
    Dcsimple_no_su_yespizero = basf2.create_path()
    ma.rankByHighest(particleList="D+:simplevetocomb_su_yespizero", variable="chiProb", outputVariable="Dsimple_rank_yespizero",path=Dcsimple_yes_su_yespizero)
    ma.cutAndCopyList("D+:simplevetocomb_BCS_su_yespizero", "D+:simplevetocomb_su_yespizero", cut="extraInfo(Dsimple_rank_yespizero)==1", path=Dcsimple_yes_su_yespizero)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_su_yespizero", {'countInList(D+:simplevetocomb_su_yespizero, chiProb>0.1)': 'nDc_yespizero'}, path=Dcsimple_yes_su_yespizero)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_su_yespizero", {'medianValueInList(D+:simplevetocomb_su_yespizero,chiProb)': 'Dc_pValue_med_yespizero'}, path=Dcsimple_yes_su_yespizero)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_su_yespizero", {'formula(averageValueInList(D+:simplevetocomb_su_yespizero, formula((chiProb-averageValueInList(D+:simplevetocomb_su_yespizero,chiProb))**2))**0.5)': 'Dc_pValue_std_yespizero'}, path=Dcsimple_yes_su_yespizero)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_su_yespizero", {'chiProb': 'Dcsimpleveto_chiProb_yespizero'}, path=Dcsimple_yes_su_yespizero)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_su_yespizero", {'dr': 'Dcsimpleveto_dr_yespizero'}, path=Dcsimple_yes_su_yespizero)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_su_yespizero", {'dz': 'Dcsimpleveto_dz_yespizero'}, path=Dcsimple_yes_su_yespizero)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_su_yespizero", {'M': 'Dcsimpleveto_M_yespizero'}, path=Dcsimple_yes_su_yespizero)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(0)': 'nDc_yespizero'}, path=Dcsimple_no_su_yespizero)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'Dc_pValue_med_yespizero'}, path=Dcsimple_no_su_yespizero)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'Dc_pValue_std_yespizero'}, path=Dcsimple_no_su_yespizero)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'Dcsimpleveto_chiProb_yespizero'}, path=Dcsimple_no_su_yespizero)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'Dcsimpleveto_dr_yespizero'}, path=Dcsimple_no_su_yespizero)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'Dcsimpleveto_dz_yespizero'}, path=Dcsimple_no_su_yespizero)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'Dcsimpleveto_M_yespizero'}, path=Dcsimple_no_su_yespizero)
    nDcsimple_su_module = roe_path_one.add_module("VariableToReturnValue", variable="nParticlesInList(D+:simplevetocomb_su_yespizero)")
    nDcsimple_su_module.if_value(">=1", Dcsimple_yes_su_yespizero, basf2.AfterConditionPath.CONTINUE)
    nDcsimple_su_module.if_value("<1", Dcsimple_no_su_yespizero, basf2.AfterConditionPath.CONTINUE)
    
    ma.reconstructDecay(decayString="D0:simpleveto1_su -> K+:sig_su pi-:simple", cut="M < 2.1", dmID=1, path=roe_path_one)
    ma.reconstructDecay(decayString="D0:simpleveto2_su -> K_S0:sig_su pi+:simple pi-:simple", cut="M < 2.1", dmID=2, path=roe_path_one)
    ma.reconstructDecay(decayString="D0:simpleveto3_su -> K+:sig_su pi-:simple pi+:simple pi-:simple", cut="M < 2.1", dmID=3, path=roe_path_one)
    ma.reconstructDecay(decayString="D0:simpleveto4_su -> K-:sig_su pi+:simple pi0:simple", cut="M < 2.1", dmID=4, path=roe_path_one)
    ma.reconstructDecay(decayString="D0:simpleveto5_su -> K_S0:sig_su pi+:simple pi-:simple pi0:simple", cut="M < 2.1", dmID=5, path=roe_path_one)
    ma.reconstructDecay(decayString="D0:simpleveto6_su -> K-:sig_su pi+:simple pi0:simple pi0:simple", cut="M < 2.1", dmID=6, path=roe_path_one)
    vertex.treeFit(list_name="D0:simpleveto1_su", conf_level=0, updateAllDaughters=False, path=roe_path_one)
    vertex.treeFit(list_name="D0:simpleveto2_su", conf_level=0, updateAllDaughters=False, path=roe_path_one)
    vertex.treeFit(list_name="D0:simpleveto3_su", conf_level=0, updateAllDaughters=False, path=roe_path_one)
    vertex.treeFit(list_name="D0:simpleveto4_su", conf_level=0, updateAllDaughters=False, path=roe_path_one)
    vertex.treeFit(list_name="D0:simpleveto5_su", conf_level=0, updateAllDaughters=False, path=roe_path_one)
    vertex.treeFit(list_name="D0:simpleveto6_su", conf_level=0, updateAllDaughters=False, path=roe_path_one)
    ma.copyLists(outputListName="D0:simplevetocomb_su", inputListNames=["D0:simpleveto1_su", "D0:simpleveto2_su", "D0:simpleveto3_su"], path=roe_path_one)
    ma.copyLists(outputListName="D0:simplevetocomb_su_yespizero", inputListNames=["D0:simpleveto1_su", "D0:simpleveto2_su", "D0:simpleveto3_su", "D0:simpleveto4_su", "D0:simpleveto5_su", "D0:simpleveto6_su"], path=roe_path_one)
    
    Dnsimple_yes_su = basf2.create_path()
    Dnsimple_no_su = basf2.create_path()
    ma.rankByHighest(particleList="D0:simplevetocomb_su", variable="chiProb", outputVariable="Dsimple_rank",path=Dnsimple_yes_su)
    ma.cutAndCopyList("D0:simplevetocomb_BCS_su", "D0:simplevetocomb_su", cut="extraInfo(Dsimple_rank)==1", path=Dnsimple_yes_su)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_su", {'countInList(D0:simplevetocomb_su, chiProb>0.1)': 'nD0'}, path=Dnsimple_yes_su)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_su", {'medianValueInList(D0:simplevetocomb_su,chiProb)': 'D0_pValue_med'}, path=Dnsimple_yes_su)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_su", {'formula(averageValueInList(D0:simplevetocomb_su, formula((chiProb-averageValueInList(D0:simplevetocomb_su,chiProb))**2))**0.5)': 'D0_pValue_std'}, path=Dnsimple_yes_su)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_su", {'chiProb': 'D0simpleveto_chiProb'}, path=Dnsimple_yes_su)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_su", {'dr': 'D0simpleveto_dr'}, path=Dnsimple_yes_su)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_su", {'dz': 'D0simpleveto_dz'}, path=Dnsimple_yes_su)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_su", {'M': 'D0simpleveto_M'}, path=Dnsimple_yes_su)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(0)': 'nD0'}, path=Dnsimple_no_su)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'D0_pValue_med'}, path=Dnsimple_no_su)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'D0_pValue_std'}, path=Dnsimple_no_su)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'D0simpleveto_chiProb'}, path=Dnsimple_no_su)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'D0simpleveto_dr'}, path=Dnsimple_no_su)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'D0simpleveto_dz'}, path=Dnsimple_no_su)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'D0simpleveto_M'}, path=Dnsimple_no_su)
    nDnsimple_su_module = roe_path_one.add_module("VariableToReturnValue", variable="nParticlesInList(D0:simplevetocomb_su)")
    nDnsimple_su_module.if_value(">=1", Dnsimple_yes_su, basf2.AfterConditionPath.CONTINUE)
    nDnsimple_su_module.if_value("<1", Dnsimple_no_su, basf2.AfterConditionPath.CONTINUE)

    Dnsimple_yes_su_yespizero = basf2.create_path()
    Dnsimple_no_su_yespizero = basf2.create_path()
    ma.rankByHighest(particleList="D0:simplevetocomb_su_yespizero", variable="chiProb", outputVariable="Dsimple_rank_yespizero",path=Dnsimple_yes_su_yespizero)
    ma.cutAndCopyList("D0:simplevetocomb_BCS_su_yespizero", "D0:simplevetocomb_su_yespizero", cut="extraInfo(Dsimple_rank_yespizero)==1", path=Dnsimple_yes_su_yespizero)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_su_yespizero", {'countInList(D0:simplevetocomb_su_yespizero, chiProb>0.1)': 'nD0_yespizero'}, path=Dnsimple_yes_su_yespizero)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_su_yespizero", {'medianValueInList(D0:simplevetocomb_su_yespizero,chiProb)': 'D0_pValue_med_yespizero'}, path=Dnsimple_yes_su_yespizero)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_su_yespizero", {'formula(averageValueInList(D0:simplevetocomb_su_yespizero, formula((chiProb-averageValueInList(D0:simplevetocomb_su_yespizero,chiProb))**2))**0.5)': 'D0_pValue_std_yespizero'}, path=Dnsimple_yes_su_yespizero)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_su_yespizero", {'chiProb': 'D0simpleveto_chiProb_yespizero'}, path=Dnsimple_yes_su_yespizero)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_su_yespizero", {'dr': 'D0simpleveto_dr_yespizero'}, path=Dnsimple_yes_su_yespizero)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_su_yespizero", {'dz': 'D0simpleveto_dz_yespizero'}, path=Dnsimple_yes_su_yespizero)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_su_yespizero", {'M': 'D0simpleveto_M_yespizero'}, path=Dnsimple_yes_su_yespizero)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(0)': 'nD0_yespizero'}, path=Dnsimple_no_su_yespizero)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'D0_pValue_med_yespizero'}, path=Dnsimple_no_su_yespizero)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'D0_pValue_std_yespizero'}, path=Dnsimple_no_su_yespizero)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'D0simpleveto_chiProb_yespizero'}, path=Dnsimple_no_su_yespizero)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'D0simpleveto_dr_yespizero'}, path=Dnsimple_no_su_yespizero)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'D0simpleveto_dz_yespizero'}, path=Dnsimple_no_su_yespizero)
    ma.variableToSignalSideExtraInfo('Xsu:sig', {'constant(-1)': 'D0simpleveto_M_yespizero'}, path=Dnsimple_no_su_yespizero)
    nDnsimple_su_module = roe_path_one.add_module("VariableToReturnValue", variable="nParticlesInList(D0:simplevetocomb_su_yespizero)")
    nDnsimple_su_module.if_value(">=1", Dnsimple_yes_su_yespizero, basf2.AfterConditionPath.CONTINUE)
    nDnsimple_su_module.if_value("<1", Dnsimple_no_su_yespizero, basf2.AfterConditionPath.CONTINUE)
    
    my_path.for_each('RestOfEvent', 'RestOfEvents', roe_path_one)
    
    ma.buildRestOfEvent("Xsd:comb", path=my_path)
    roe_path_two = basf2.Path()
    deadEndPath_two = basf2.Path()
    ma.signalSideParticleFilter('Xsd:comb', '', roe_path_two, deadEndPath_two)
    ma.fillSignalSideParticleList('Xsd:sig', '^Xsd:comb', roe_path_two)
    ma.cutAndCopyList("K+:sig_sd", "K+:mychargedKaon", cut="isDescendantOfList(Xsd:sig, 1) == 1", path=roe_path_two)
    ma.cutAndCopyList("K_S0:sig_sd", "K_S0:myKaonshort", cut="isDescendantOfList(Xsd:sig, 1) == 1", path=roe_path_two)
    ma.cutAndCopyList("pi+:sig_sd", "pi+:mychargedPion", cut="isDescendantOfList(Xsd:sig, 1) == 1", path=roe_path_two)
    ma.cutAndCopyList("pi0:sig_sd", "pi0:myneutralPion", cut="isDescendantOfList(Xsd:sig, 1) == 1", path=roe_path_two)
    ma.matchMCTruth("K+:sig_sd", path=roe_path_two)
    ma.matchMCTruth("pi+:sig_sd", path=roe_path_two)
    
    # uncertainty of slow track
    ma.cutAndCopyList("K+:slow1_sd", "K+:mychargedKaon", cut="isDescendantOfList(Xsd:sig, 1) == 1 and 0.05 < p < 0.12", path=roe_path_two)
    ma.cutAndCopyList("K+:slow2_sd", "K+:mychargedKaon", cut="isDescendantOfList(Xsd:sig, 1) == 1 and 0.12 < p < 0.16", path=roe_path_two)
    ma.cutAndCopyList("K+:slow3_sd", "K+:mychargedKaon", cut="isDescendantOfList(Xsd:sig, 1) == 1 and 0.16 < p < 0.20", path=roe_path_two)
    ma.cutAndCopyList("pi+:slow1_sd", "pi+:mychargedPion", cut="isDescendantOfList(Xsd:sig, 1) == 1 and 0.05 < p < 0.12", path=roe_path_two)
    ma.cutAndCopyList("pi+:slow2_sd", "pi+:mychargedPion", cut="isDescendantOfList(Xsd:sig, 1) == 1 and 0.12 < p < 0.16", path=roe_path_two)
    ma.cutAndCopyList("pi+:slow3_sd", "pi+:mychargedPion", cut="isDescendantOfList(Xsd:sig, 1) == 1 and 0.16 < p < 0.20", path=roe_path_two)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'nParticlesInList(K+:slow1_sd)': 'nKslow1'}, path=roe_path_two)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'nParticlesInList(K+:slow2_sd)': 'nKslow2'}, path=roe_path_two)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'nParticlesInList(K+:slow3_sd)': 'nKslow3'}, path=roe_path_two)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'nParticlesInList(pi+:slow1_sd)': 'nPislow1'}, path=roe_path_two)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'nParticlesInList(pi+:slow2_sd)': 'nPislow2'}, path=roe_path_two)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'nParticlesInList(pi+:slow3_sd)': 'nPislow3'}, path=roe_path_two)
    
    # uncertainty of PID
    pbins = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
    cosThetabins = [-0.866, -0.682, -0.4226, -0.1045, 0.225, 0.5, 0.766, 0.8829, 0.9563]
    index = 0
    for i_p in range(len(pbins) - 1):
        for i_cosTheta in range(len(cosThetabins) - 1):
            pcut_string = "[" + str(pbins[i_p]) + " < p < " + str(pbins[i_p+1]) + "]"
            cosThetacut_string = "[" + str(cosThetabins[i_cosTheta]) + " < cosTheta < " + str(cosThetabins[i_cosTheta+1]) + "]"
            ma.cutAndCopyList("K+:truePIDbin" + str(index), "K+:sig_sd", cut="[abs(mcPDG) == 321] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + pcut_string + " and " + cosThetacut_string, path=roe_path_two)
            ma.cutAndCopyList("K+:misPIDbin" + str(index), "K+:sig_sd", cut="[abs(mcPDG) == 211] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + pcut_string + " and " + cosThetacut_string, path=roe_path_two)
            ma.cutAndCopyList("pi+:truePIDbin" + str(index), "pi+:sig_sd", cut="[abs(mcPDG) == 211] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + pcut_string + " and " + cosThetacut_string, path=roe_path_two)
            ma.cutAndCopyList("pi+:misPIDbin" + str(index), "pi+:sig_sd", cut="[abs(mcPDG) == 321] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + pcut_string + " and " + cosThetacut_string, path=roe_path_two)
            index = index + 1
    pcut_string = "[ p < " + str(pbins[0]) + "]" + " or " + "[ p > " + str(pbins[-1]) + "]"
    cosThetacut_string = "[ cosTheta < " + str(cosThetabins[0]) + "]" + " or " + "[ cosTheta > " + str(cosThetabins[-1]) + "]"
    ma.cutAndCopyList("K+:truePIDbin" + str(index), "K+:sig_sd", cut="[abs(mcPDG) == 321] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + cosThetacut_string + "]", path=roe_path_two)
    ma.cutAndCopyList("K+:misPIDbin" + str(index), "K+:sig_sd", cut="[abs(mcPDG) == 211] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + cosThetacut_string + "]", path=roe_path_two)
    ma.cutAndCopyList("pi+:truePIDbin" + str(index), "pi+:sig_sd", cut="[abs(mcPDG) == 211] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + cosThetacut_string + "]", path=roe_path_two)
    ma.cutAndCopyList("pi+:misPIDbin" + str(index), "pi+:sig_sd", cut="[abs(mcPDG) == 321] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + cosThetacut_string + "]", path=roe_path_two)
    index = 0
    for i_p in range(len(pbins) - 1):
        for i_cosTheta in range(len(cosThetabins) - 1):
            ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(K+:truePIDbin" + str(index) + ")": "nKtruebin"+str(index)}, path=roe_path_two)
            ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(K+:misPIDbin" + str(index) + ")": "nKmisbin"+str(index)}, path=roe_path_two)
            ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(pi+:truePIDbin" + str(index) + ")": "npitruebin"+str(index)}, path=roe_path_two)
            ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(pi+:misPIDbin" + str(index) + ")": "npimisbin"+str(index)}, path=roe_path_two)
            index = index + 1
    ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(K+:truePIDbin" + str(index) + ")": "nKtruebin"+str(index)}, path=roe_path_two)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(K+:misPIDbin" + str(index) + ")": "nKmisbin"+str(index)}, path=roe_path_two)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(pi+:truePIDbin" + str(index) + ")": "npitruebin"+str(index)}, path=roe_path_two)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(pi+:misPIDbin" + str(index) + ")": "npimisbin"+str(index)}, path=roe_path_two)
    
    pbins_e = [0.4, 0.5, 1.0, 1.5, 2, 2.5, 3]
    thetabins_e = [0.22, 0.56, 1.13, 1.57, 1.88, 2.23, 2.71]
    index = 0
    for i_p in range(len(pbins_e) - 1):
        for i_theta in range(len(thetabins_e) - 1):
            pcut_string = "[" + str(pbins_e[i_p]) + " < p < " + str(pbins_e[i_p+1]) + "]"
            thetacut_string = "[" + str(thetabins_e[i_theta]) + " < theta < " + str(thetabins_e[i_theta+1]) + "]"
            ma.cutAndCopyList("K+:fakeEbin_n" + str(index), "K+:sig_sd", cut="[mcPDG == 11] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + pcut_string + " and " + thetacut_string, path=roe_path_two)
            ma.cutAndCopyList("K+:fakeEbin_p" + str(index), "K+:sig_sd", cut="[mcPDG == -11] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + pcut_string + " and " + thetacut_string, path=roe_path_two)
            ma.cutAndCopyList("pi+:fakeEbin_n" + str(index), "pi+:sig_sd", cut="[mcPDG == 11] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + pcut_string + " and " + thetacut_string, path=roe_path_two)
            ma.cutAndCopyList("pi+:fakeEbin_p" + str(index), "pi+:sig_sd", cut="[mcPDG == -11] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + pcut_string + " and " + thetacut_string, path=roe_path_two)
            index = index + 1
    pcut_string = "[ p < " + str(pbins_e[0]) + "]" + " or " + "[ p > " + str(pbins_e[-1]) + "]"
    thetacut_string = "[ theta < " + str(thetabins_e[0]) + "]" + " or " + "[ theta > " + str(thetabins_e[-1]) + "]"
    ma.cutAndCopyList("K+:fakeEbin_n" + str(index), "K+:sig_sd", cut="[mcPDG == 11] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + thetacut_string + "]", path=roe_path_two)
    ma.cutAndCopyList("K+:fakeEbin_p" + str(index), "K+:sig_sd", cut="[mcPDG == -11] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + thetacut_string + "]", path=roe_path_two)
    ma.cutAndCopyList("pi+:fakeEbin_n" + str(index), "pi+:sig_sd", cut="[mcPDG == 11] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + thetacut_string + "]", path=roe_path_two)
    ma.cutAndCopyList("pi+:fakeEbin_p" + str(index), "pi+:sig_sd", cut="[mcPDG == -11] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + thetacut_string + "]", path=roe_path_two)
    index = 0
    for i_p in range(len(pbins_e) - 1):
        for i_theta in range(len(thetabins_e) - 1):
            ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(K+:fakeEbin_n" + str(index) + ")": "nKfakeEbin_n"+str(index)}, path=roe_path_two)
            ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(K+:fakeEbin_p" + str(index) + ")": "nKfakeEbin_p"+str(index)}, path=roe_path_two)
            ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(pi+:fakeEbin_n" + str(index) + ")": "npifakeEbin_n"+str(index)}, path=roe_path_two)
            ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(pi+:fakeEbin_p" + str(index) + ")": "npifakeEbin_p"+str(index)}, path=roe_path_two)
            index = index + 1
    ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(K+:fakeEbin_n" + str(index) + ")": "nKfakeEbin_n"+str(index)}, path=roe_path_two)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(K+:fakeEbin_p" + str(index) + ")": "nKfakeEbin_p"+str(index)}, path=roe_path_two)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(pi+:fakeEbin_n" + str(index) + ")": "npifakeEbin_n"+str(index)}, path=roe_path_two)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(pi+:fakeEbin_p" + str(index) + ")": "npifakeEbin_p"+str(index)}, path=roe_path_two)
    
    pbins_mu = [0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5]
    thetabins_mu = [0.4, 0.64, 0.82, 1.16, 1.46, 1.78, 2.13, 2.22, 2.6]
    index = 0
    for i_p in range(len(pbins_mu) - 1):
        for i_theta in range(len(thetabins_mu) - 1):
            pcut_string = "[" + str(pbins_mu[i_p]) + " < p < " + str(pbins_mu[i_p+1]) + "]"
            thetacut_string = "[" + str(thetabins_mu[i_theta]) + " < theta < " + str(thetabins_mu[i_theta+1]) + "]"
            ma.cutAndCopyList("K+:fakeMUbin_n" + str(index), "K+:sig_sd", cut="[mcPDG == 13] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + pcut_string + " and " + thetacut_string, path=roe_path_two)
            ma.cutAndCopyList("K+:fakeMUbin_p" + str(index), "K+:sig_sd", cut="[mcPDG == -13] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + pcut_string + " and " + thetacut_string, path=roe_path_two)
            ma.cutAndCopyList("pi+:fakeMUbin_n" + str(index), "pi+:sig_sd", cut="[mcPDG == 13] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + pcut_string + " and " + thetacut_string, path=roe_path_two)
            ma.cutAndCopyList("pi+:fakeMUbin_p" + str(index), "pi+:sig_sd", cut="[mcPDG == -13] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + pcut_string + " and " + thetacut_string, path=roe_path_two)
            index = index + 1
    pcut_string = "[ p < " + str(pbins_mu[0]) + "]" + " or " + "[ p > " + str(pbins_mu[-1]) + "]"
    thetacut_string = "[ theta < " + str(thetabins_mu[0]) + "]" + " or " + "[ theta > " + str(thetabins_mu[-1]) + "]"
    ma.cutAndCopyList("K+:fakeMUbin_n" + str(index), "K+:sig_sd", cut="[mcPDG == 13] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + thetacut_string + "]", path=roe_path_two)
    ma.cutAndCopyList("K+:fakeMUbin_p" + str(index), "K+:sig_sd", cut="[mcPDG == -13] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + thetacut_string + "]", path=roe_path_two)
    ma.cutAndCopyList("pi+:fakeMUbin_n" + str(index), "pi+:sig_sd", cut="[mcPDG == 13] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + thetacut_string + "]", path=roe_path_two)
    ma.cutAndCopyList("pi+:fakeMUbin_p" + str(index), "pi+:sig_sd", cut="[mcPDG == -13] and [isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + "[" +  pcut_string + " or " + thetacut_string + "]", path=roe_path_two)
    index = 0
    for i_p in range(len(pbins_mu) - 1):
        for i_theta in range(len(thetabins_mu) - 1):
            ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(K+:fakeMUbin_n" + str(index) + ")": "nKfakeMUbin_n"+str(index)}, path=roe_path_two)
            ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(K+:fakeMUbin_p" + str(index) + ")": "nKfakeMUbin_p"+str(index)}, path=roe_path_two)
            ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(pi+:fakeMUbin_n" + str(index) + ")": "npifakeMUbin_n"+str(index)}, path=roe_path_two)
            ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(pi+:fakeMUbin_p" + str(index) + ")": "npifakeMUbin_p"+str(index)}, path=roe_path_two)
            index = index + 1
    ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(K+:fakeMUbin_n" + str(index) + ")": "nKfakeMUbin_n"+str(index)}, path=roe_path_two)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(K+:fakeMUbin_p" + str(index) + ")": "nKfakeMUbin_p"+str(index)}, path=roe_path_two)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(pi+:fakeMUbin_n" + str(index) + ")": "npifakeMUbin_n"+str(index)}, path=roe_path_two)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(pi+:fakeMUbin_p" + str(index) + ")": "npifakeMUbin_p"+str(index)}, path=roe_path_two)
    
    # uncertainty of pi0
    pi0_pbins = [0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0]
    index = 0
    for i_p in range(len(pi0_pbins) - 1):
        pcut_string = "[" + str(pi0_pbins[i_p]) + " < p < " + str(pi0_pbins[i_p+1]) + "]"
        ma.cutAndCopyList("pi0:pbin" + str(index), "pi0:sig_sd", cut="[isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + pcut_string, path=roe_path_two)
        index = index + 1
    pcut_string = "[ p < " + str(pi0_pbins[0]) + "]" + " or " + "[ p > " + str(pi0_pbins[-1]) + "]"
    ma.cutAndCopyList("pi0:pbin" + str(index), "pi0:sig_sd", cut="[isDescendantOfList(Xsd:sig, 1) == 1]" + " and " + "[" + pcut_string + "]", path=roe_path_two)
    index = 0
    for i_p in range(len(pi0_pbins) - 1):
        ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(pi0:pbin" + str(index) + ")": "npi0bin"+str(index)}, path=roe_path_two)
        index = index + 1
    ma.variableToSignalSideExtraInfo('Xsd:sig', {"nParticlesInList(pi0:pbin" + str(index) + ")": "npi0bin"+str(index)}, path=roe_path_two)
    
    # D veto all
    ma.fillParticleList(decayString="pi+:simple", cut="dr < 0.5 and abs(dz) < 2 and nPXDHits > 0",path=roe_path_two)
    ma.cutAndCopyList("pi0:simple", "pi0:eff30_May2020", cut="p > 0.4", path=roe_path_two)
    if(analysistype=="Jpsi"):
        ma.applyCuts("pi+:simple", "sourceObjectIsInList(mu+:fromJpsimuchargedMuon)==0", path=roe_path_two)
    ma.reconstructDecay(decayString="D+:simpleveto1_sd -> K+:sig_sd pi-:simple pi+:simple", cut="M < 2.1", dmID=1, path=roe_path_two) # DCS
    ma.reconstructDecay(decayString="D+:simpleveto2_sd -> K_S0:sig_sd pi+:simple", cut="M < 2.1", dmID=2, path=roe_path_two)
    ma.reconstructDecay(decayString="D+:simpleveto3_sd -> K-:sig_sd pi+:simple pi+:simple", cut="M < 2.1", dmID=3, path=roe_path_two)
    ma.reconstructDecay(decayString="D+:simpleveto4_sd -> K_S0:sig_sd pi+:simple pi0:simple", cut="M < 2.1", dmID=4, path=roe_path_two) # pi0
    ma.reconstructDecay(decayString="D+:simpleveto5_sd -> K-:sig_sd pi+:simple pi+:simple pi0:simple", cut="M < 2.1", dmID=5, path=roe_path_two) # pi0
    ma.reconstructDecay(decayString="D+:simpleveto6_sd -> K_S0:sig_sd pi+:simple pi0:simple pi0:simple", cut="M < 2.1", dmID=6, path=roe_path_two) # pi0
    vertex.treeFit(list_name="D+:simpleveto1_sd", conf_level=0, updateAllDaughters=False, path=roe_path_two)
    vertex.treeFit(list_name="D+:simpleveto2_sd", conf_level=0, updateAllDaughters=False, path=roe_path_two)
    vertex.treeFit(list_name="D+:simpleveto3_sd", conf_level=0, updateAllDaughters=False, path=roe_path_two)
    vertex.treeFit(list_name="D+:simpleveto4_sd", conf_level=0, updateAllDaughters=False, path=roe_path_two)
    vertex.treeFit(list_name="D+:simpleveto5_sd", conf_level=0, updateAllDaughters=False, path=roe_path_two)
    vertex.treeFit(list_name="D+:simpleveto6_sd", conf_level=0, updateAllDaughters=False, path=roe_path_two)
    ma.copyLists(outputListName="D+:simplevetocomb_sd", inputListNames=["D+:simpleveto1_sd", "D+:simpleveto2_sd", "D+:simpleveto3_sd"], path=roe_path_two)
    ma.copyLists(outputListName="D+:simplevetocomb_sd_noDCS", inputListNames=["D+:simpleveto2_sd", "D+:simpleveto3_sd"], path=roe_path_two)
    ma.copyLists(outputListName="D+:simplevetocomb_sd_yespizero", inputListNames=["D+:simpleveto1_sd", "D+:simpleveto2_sd", "D+:simpleveto3_sd", "D+:simpleveto4_sd", "D+:simpleveto5_sd", "D+:simpleveto6_sd"], path=roe_path_two)
    
    Dcsimple_yes_sd = basf2.create_path()
    Dcsimple_no_sd = basf2.create_path()
    ma.rankByHighest(particleList="D+:simplevetocomb_sd", variable="chiProb", outputVariable="Dsimple_rank",path=Dcsimple_yes_sd)
    ma.cutAndCopyList("D+:simplevetocomb_BCS_sd", "D+:simplevetocomb_sd", cut="extraInfo(Dsimple_rank)==1", path=Dcsimple_yes_sd)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_sd", {'countInList(D+:simplevetocomb_sd, chiProb>0.1)': 'nDc'}, path=Dcsimple_yes_sd)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_sd", {'medianValueInList(D+:simplevetocomb_sd,chiProb)': 'Dc_pValue_med'}, path=Dcsimple_yes_sd)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_sd", {'formula(averageValueInList(D+:simplevetocomb_sd, formula((chiProb-averageValueInList(D+:simplevetocomb_sd,chiProb))**2))**0.5)': 'Dc_pValue_std'}, path=Dcsimple_yes_sd)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_sd", {'chiProb': 'Dcsimpleveto_chiProb'}, path=Dcsimple_yes_sd)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_sd", {'dr': 'Dcsimpleveto_dr'}, path=Dcsimple_yes_sd)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_sd", {'dz': 'Dcsimpleveto_dz'}, path=Dcsimple_yes_sd)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_sd", {'M': 'Dcsimpleveto_M'}, path=Dcsimple_yes_sd)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(0)': 'nDc'}, path=Dcsimple_no_sd)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'Dc_pValue_med'}, path=Dcsimple_no_sd)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'Dc_pValue_std'}, path=Dcsimple_no_sd)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'Dcsimpleveto_chiProb'}, path=Dcsimple_no_sd)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'Dcsimpleveto_dr'}, path=Dcsimple_no_sd)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'Dcsimpleveto_dz'}, path=Dcsimple_no_sd)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'Dcsimpleveto_M'}, path=Dcsimple_no_sd)
    nDcsimple_sd_module = roe_path_two.add_module("VariableToReturnValue", variable="nParticlesInList(D+:simplevetocomb_sd)")
    nDcsimple_sd_module.if_value(">=1", Dcsimple_yes_sd, basf2.AfterConditionPath.CONTINUE)
    nDcsimple_sd_module.if_value("<1", Dcsimple_no_sd, basf2.AfterConditionPath.CONTINUE)

    Dcsimple_yes_sd_noDCS = basf2.create_path()
    Dcsimple_no_sd_noDCS = basf2.create_path()
    ma.rankByHighest(particleList="D+:simplevetocomb_sd_noDCS", variable="chiProb", outputVariable="Dsimple_rank_noDCS",path=Dcsimple_yes_sd_noDCS)
    ma.cutAndCopyList("D+:simplevetocomb_BCS_sd_noDCS", "D+:simplevetocomb_sd_noDCS", cut="extraInfo(Dsimple_rank_noDCS)==1", path=Dcsimple_yes_sd_noDCS)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_sd_noDCS", {'countInList(D+:simplevetocomb_sd_noDCS, chiProb>0.1)': 'nDc_noDCS'}, path=Dcsimple_yes_sd_noDCS)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_sd_noDCS", {'medianValueInList(D+:simplevetocomb_sd_noDCS,chiProb)': 'Dc_pValue_med_noDCS'}, path=Dcsimple_yes_sd_noDCS)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_sd_noDCS", {'formula(averageValueInList(D+:simplevetocomb_sd_noDCS, formula((chiProb-averageValueInList(D+:simplevetocomb_sd_noDCS,chiProb))**2))**0.5)': 'Dc_pValue_std_noDCS'}, path=Dcsimple_yes_sd_noDCS)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_sd_noDCS", {'chiProb': 'Dcsimpleveto_chiProb_noDCS'}, path=Dcsimple_yes_sd_noDCS)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_sd_noDCS", {'dr': 'Dcsimpleveto_dr_noDCS'}, path=Dcsimple_yes_sd_noDCS)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_sd_noDCS", {'dz': 'Dcsimpleveto_dz_noDCS'}, path=Dcsimple_yes_sd_noDCS)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_sd_noDCS", {'M': 'Dcsimpleveto_M_noDCS'}, path=Dcsimple_yes_sd_noDCS)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(0)': 'nDc_noDCS'}, path=Dcsimple_no_sd_noDCS)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'Dc_pValue_med_noDCS'}, path=Dcsimple_no_sd_noDCS)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'Dc_pValue_std_noDCS'}, path=Dcsimple_no_sd_noDCS)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'Dcsimpleveto_chiProb_noDCS'}, path=Dcsimple_no_sd_noDCS)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'Dcsimpleveto_dr_noDCS'}, path=Dcsimple_no_sd_noDCS)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'Dcsimpleveto_dz_noDCS'}, path=Dcsimple_no_sd_noDCS)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'Dcsimpleveto_M_noDCS'}, path=Dcsimple_no_sd_noDCS)
    nDcsimple_sd_module = roe_path_two.add_module("VariableToReturnValue", variable="nParticlesInList(D+:simplevetocomb_sd_noDCS)")
    nDcsimple_sd_module.if_value(">=1", Dcsimple_yes_sd_noDCS, basf2.AfterConditionPath.CONTINUE)
    nDcsimple_sd_module.if_value("<1", Dcsimple_no_sd_noDCS, basf2.AfterConditionPath.CONTINUE)

    Dcsimple_yes_sd_yespizero = basf2.create_path()
    Dcsimple_no_sd_yespizero = basf2.create_path()
    ma.rankByHighest(particleList="D+:simplevetocomb_sd_yespizero", variable="chiProb", outputVariable="Dsimple_rank_yespizero",path=Dcsimple_yes_sd_yespizero)
    ma.cutAndCopyList("D+:simplevetocomb_BCS_sd_yespizero", "D+:simplevetocomb_sd_yespizero", cut="extraInfo(Dsimple_rank_yespizero)==1", path=Dcsimple_yes_sd_yespizero)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_sd_yespizero", {'countInList(D+:simplevetocomb_sd_yespizero, chiProb>0.1)': 'nDc_yespizero'}, path=Dcsimple_yes_sd_yespizero)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_sd_yespizero", {'medianValueInList(D+:simplevetocomb_sd_yespizero,chiProb)': 'Dc_pValue_med_yespizero'}, path=Dcsimple_yes_sd_yespizero)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_sd_yespizero", {'formula(averageValueInList(D+:simplevetocomb_sd_yespizero, formula((chiProb-averageValueInList(D+:simplevetocomb_sd_yespizero,chiProb))**2))**0.5)': 'Dc_pValue_std_yespizero'}, path=Dcsimple_yes_sd_yespizero)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_sd_yespizero", {'chiProb': 'Dcsimpleveto_chiProb_yespizero'}, path=Dcsimple_yes_sd_yespizero)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_sd_yespizero", {'dr': 'Dcsimpleveto_dr_yespizero'}, path=Dcsimple_yes_sd_yespizero)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_sd_yespizero", {'dz': 'Dcsimpleveto_dz_yespizero'}, path=Dcsimple_yes_sd_yespizero)
    ma.variableToSignalSideExtraInfo("D+:simplevetocomb_BCS_sd_yespizero", {'M': 'Dcsimpleveto_M_yespizero'}, path=Dcsimple_yes_sd_yespizero)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(0)': 'nDc_yespizero'}, path=Dcsimple_no_sd_yespizero)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'Dc_pValue_med_yespizero'}, path=Dcsimple_no_sd_yespizero)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'Dc_pValue_std_yespizero'}, path=Dcsimple_no_sd_yespizero)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'Dcsimpleveto_chiProb_yespizero'}, path=Dcsimple_no_sd_yespizero)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'Dcsimpleveto_dr_yespizero'}, path=Dcsimple_no_sd_yespizero)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'Dcsimpleveto_dz_yespizero'}, path=Dcsimple_no_sd_yespizero)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'Dcsimpleveto_M_yespizero'}, path=Dcsimple_no_sd_yespizero)
    nDcsimple_sd_module = roe_path_two.add_module("VariableToReturnValue", variable="nParticlesInList(D+:simplevetocomb_sd_yespizero)")
    nDcsimple_sd_module.if_value(">=1", Dcsimple_yes_sd_yespizero, basf2.AfterConditionPath.CONTINUE)
    nDcsimple_sd_module.if_value("<1", Dcsimple_no_sd_yespizero, basf2.AfterConditionPath.CONTINUE)

    ma.reconstructDecay(decayString="D0:simpleveto1_sd -> K+:sig_sd pi-:simple", cut="M < 2.1", dmID=1, path=roe_path_two)
    ma.reconstructDecay(decayString="D0:simpleveto2_sd -> K_S0:sig_sd pi+:simple pi-:simple", cut="M < 2.1", dmID=2, path=roe_path_two)
    ma.reconstructDecay(decayString="D0:simpleveto3_sd -> K+:sig_sd pi-:simple pi+:simple pi-:simple", cut="M < 2.1", dmID=3, path=roe_path_two)
    ma.reconstructDecay(decayString="D0:simpleveto4_sd -> K-:sig_sd pi+:simple pi0:simple", cut="M < 2.1", dmID=4, path=roe_path_two)
    ma.reconstructDecay(decayString="D0:simpleveto5_sd -> K_S0:sig_sd pi+:simple pi-:simple pi0:simple", cut="M < 2.1", dmID=5, path=roe_path_two)
    ma.reconstructDecay(decayString="D0:simpleveto6_sd -> K-:sig_sd pi+:simple pi0:simple pi0:simple", cut="M < 2.1", dmID=6, path=roe_path_two)
    vertex.treeFit(list_name="D0:simpleveto1_sd", conf_level=0, updateAllDaughters=False, path=roe_path_two)
    vertex.treeFit(list_name="D0:simpleveto2_sd", conf_level=0, updateAllDaughters=False, path=roe_path_two)
    vertex.treeFit(list_name="D0:simpleveto3_sd", conf_level=0, updateAllDaughters=False, path=roe_path_two)
    vertex.treeFit(list_name="D0:simpleveto4_sd", conf_level=0, updateAllDaughters=False, path=roe_path_two)
    vertex.treeFit(list_name="D0:simpleveto5_sd", conf_level=0, updateAllDaughters=False, path=roe_path_two)
    vertex.treeFit(list_name="D0:simpleveto6_sd", conf_level=0, updateAllDaughters=False, path=roe_path_two)
    ma.copyLists(outputListName="D0:simplevetocomb_sd", inputListNames=["D0:simpleveto1_sd", "D0:simpleveto2_sd", "D0:simpleveto3_sd"], path=roe_path_two)
    ma.copyLists(outputListName="D0:simplevetocomb_sd_yespizero", inputListNames=["D0:simpleveto1_sd", "D0:simpleveto2_sd", "D0:simpleveto3_sd", "D0:simpleveto4_sd", "D0:simpleveto5_sd", "D0:simpleveto6_sd"], path=roe_path_two)    
    
    Dnsimple_yes_sd = basf2.create_path()
    Dnsimple_no_sd = basf2.create_path()
    ma.rankByHighest(particleList="D0:simplevetocomb_sd", variable="chiProb", outputVariable="Dsimple_rank",path=Dnsimple_yes_sd)
    ma.cutAndCopyList("D0:simplevetocomb_BCS_sd", "D0:simplevetocomb_sd", cut="extraInfo(Dsimple_rank)==1", path=Dnsimple_yes_sd)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_sd", {'countInList(D0:simplevetocomb_sd, chiProb>0.1)': 'nD0'}, path=Dnsimple_yes_sd)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_sd", {'medianValueInList(D0:simplevetocomb_sd,chiProb)': 'D0_pValue_med'}, path=Dnsimple_yes_sd)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_sd", {'formula(averageValueInList(D0:simplevetocomb_sd, formula((chiProb-averageValueInList(D0:simplevetocomb_sd,chiProb))**2))**0.5)': 'D0_pValue_std'}, path=Dnsimple_yes_sd)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_sd", {'chiProb': 'D0simpleveto_chiProb'}, path=Dnsimple_yes_sd)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_sd", {'dr': 'D0simpleveto_dr'}, path=Dnsimple_yes_sd)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_sd", {'dz': 'D0simpleveto_dz'}, path=Dnsimple_yes_sd)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_sd", {'M': 'D0simpleveto_M'}, path=Dnsimple_yes_sd)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(0)': 'nD0'}, path=Dnsimple_no_sd)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'D0_pValue_med'}, path=Dnsimple_no_sd)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'D0_pValue_std'}, path=Dnsimple_no_sd)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'D0simpleveto_chiProb'}, path=Dnsimple_no_sd)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'D0simpleveto_dr'}, path=Dnsimple_no_sd)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'D0simpleveto_dz'}, path=Dnsimple_no_sd)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'D0simpleveto_M'}, path=Dnsimple_no_sd)
    nDnsimple_sd_module = roe_path_two.add_module("VariableToReturnValue", variable="nParticlesInList(D0:simplevetocomb_sd)")
    nDnsimple_sd_module.if_value(">=1", Dnsimple_yes_sd, basf2.AfterConditionPath.CONTINUE)
    nDnsimple_sd_module.if_value("<1", Dnsimple_no_sd, basf2.AfterConditionPath.CONTINUE)

    Dnsimple_yes_sd_yespizero = basf2.create_path()
    Dnsimple_no_sd_yespizero = basf2.create_path()
    ma.rankByHighest(particleList="D0:simplevetocomb_sd_yespizero", variable="chiProb", outputVariable="Dsimple_rank_yespizero",path=Dnsimple_yes_sd_yespizero)
    ma.cutAndCopyList("D0:simplevetocomb_BCS_sd_yespizero", "D0:simplevetocomb_sd_yespizero", cut="extraInfo(Dsimple_rank_yespizero)==1", path=Dnsimple_yes_sd_yespizero)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_sd_yespizero", {'countInList(D0:simplevetocomb_sd_yespizero, chiProb>0.1)': 'nD0_yespizero'}, path=Dnsimple_yes_sd_yespizero)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_sd_yespizero", {'medianValueInList(D0:simplevetocomb_sd_yespizero,chiProb)': 'D0_pValue_med_yespizero'}, path=Dnsimple_yes_sd_yespizero)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_sd_yespizero", {'formula(averageValueInList(D0:simplevetocomb_sd_yespizero, formula((chiProb-averageValueInList(D0:simplevetocomb_sd_yespizero,chiProb))**2))**0.5)': 'D0_pValue_std_yespizero'}, path=Dnsimple_yes_sd_yespizero)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_sd_yespizero", {'chiProb': 'D0simpleveto_chiProb_yespizero'}, path=Dnsimple_yes_sd_yespizero)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_sd_yespizero", {'dr': 'D0simpleveto_dr_yespizero'}, path=Dnsimple_yes_sd_yespizero)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_sd_yespizero", {'dz': 'D0simpleveto_dz_yespizero'}, path=Dnsimple_yes_sd_yespizero)
    ma.variableToSignalSideExtraInfo("D0:simplevetocomb_BCS_sd_yespizero", {'M': 'D0simpleveto_M_yespizero'}, path=Dnsimple_yes_sd_yespizero)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(0)': 'nD0_yespizero'}, path=Dnsimple_no_sd_yespizero)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'D0_pValue_med_yespizero'}, path=Dnsimple_no_sd_yespizero)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'D0_pValue_std_yespizero'}, path=Dnsimple_no_sd_yespizero)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'D0simpleveto_chiProb_yespizero'}, path=Dnsimple_no_sd_yespizero)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'D0simpleveto_dr_yespizero'}, path=Dnsimple_no_sd_yespizero)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'D0simpleveto_dz_yespizero'}, path=Dnsimple_no_sd_yespizero)
    ma.variableToSignalSideExtraInfo('Xsd:sig', {'constant(-1)': 'D0simpleveto_M_yespizero'}, path=Dnsimple_no_sd_yespizero)
    nDnsimple_sd_module = roe_path_two.add_module("VariableToReturnValue", variable="nParticlesInList(D0:simplevetocomb_sd_yespizero)")
    nDnsimple_sd_module.if_value(">=1", Dnsimple_yes_sd_yespizero, basf2.AfterConditionPath.CONTINUE)
    nDnsimple_sd_module.if_value("<1", Dnsimple_no_sd_yespizero, basf2.AfterConditionPath.CONTINUE)
    
    my_path.for_each('RestOfEvent', 'RestOfEvents', roe_path_two)
    
    
    # treefit and get vertex information
    roe_path_three = basf2.Path()
    deadEndPath_three = basf2.Path()
    ma.signalSideParticleFilter('Xsu:comb', '[extraInfo(decayModeID) != 0] and [extraInfo(decayModeID) != 1] and [extraInfo(decayModeID) != 9]', roe_path_three, deadEndPath_three)
    ma.fillSignalSideParticleList('Xsu:sig_tree', '^Xsu:comb', roe_path_three)
    vertex.treeFit(list_name="Xsu:sig_tree", conf_level=-1, updateAllDaughters=False, path=roe_path_three)
    ma.variableToSignalSideExtraInfo("Xsu:sig_tree", {'dr': 'mydr'}, path=roe_path_three)
    ma.variableToSignalSideExtraInfo("Xsu:sig_tree", {'dz': 'mydz'}, path=roe_path_three)
    ma.variableToSignalSideExtraInfo("Xsu:sig_tree", {'chiProb': 'mychiProb'}, path=roe_path_three)
    my_path.for_each('RestOfEvent', 'RestOfEvents', roe_path_three)
    
    roe_path_four = basf2.Path()
    deadEndPath_four = basf2.Path()
    ma.signalSideParticleFilter('Xsu:comb', '[extraInfo(decayModeID) == 0] or [extraInfo(decayModeID) == 1] or [extraInfo(decayModeID) == 9]', roe_path_four, deadEndPath_four)
    ma.fillSignalSideParticleList('Xsu:sig_tree', '^Xsu:comb', roe_path_four)
    ma.variableToSignalSideExtraInfo("Xsu:sig_tree", {'daughter(0, dr)': 'mydr'}, path=roe_path_four)
    ma.variableToSignalSideExtraInfo("Xsu:sig_tree", {'daughter(0, dz)': 'mydz'}, path=roe_path_four)
    ma.variableToSignalSideExtraInfo("Xsu:sig_tree", {'constant(2)': 'mychiProb'}, path=roe_path_four)
    my_path.for_each('RestOfEvent', 'RestOfEvents', roe_path_four)
    
    roe_path_five = basf2.Path()
    deadEndPath_five = basf2.Path()
    ma.signalSideParticleFilter('Xsd:comb', '[extraInfo(decayModeID) != 0] and [extraInfo(decayModeID) != 2] and [extraInfo(decayModeID) != 9]', roe_path_five, deadEndPath_five)
    ma.fillSignalSideParticleList('Xsd:sig_tree', '^Xsd:comb', roe_path_five)
    vertex.treeFit(list_name="Xsd:sig_tree", conf_level=-1, updateAllDaughters=False, path=roe_path_five)
    ma.variableToSignalSideExtraInfo("Xsd:sig_tree", {'dr': 'mydr'}, path=roe_path_five)
    ma.variableToSignalSideExtraInfo("Xsd:sig_tree", {'dz': 'mydz'}, path=roe_path_five)
    ma.variableToSignalSideExtraInfo("Xsd:sig_tree", {'chiProb': 'mychiProb'}, path=roe_path_five)
    my_path.for_each('RestOfEvent', 'RestOfEvents', roe_path_five)
    
    roe_path_six = basf2.Path()
    deadEndPath_six = basf2.Path()
    ma.signalSideParticleFilter('Xsd:comb', '[extraInfo(decayModeID) == 0] or [extraInfo(decayModeID) == 2] or [extraInfo(decayModeID) == 9]', roe_path_six, deadEndPath_six)
    ma.fillSignalSideParticleList('Xsd:sig_tree', '^Xsd:comb', roe_path_six)
    ma.variableToSignalSideExtraInfo("Xsd:sig_tree", {'daughter(0, dr)': 'mydr'}, path=roe_path_six)
    ma.variableToSignalSideExtraInfo("Xsd:sig_tree", {'daughter(0, dz)': 'mydz'}, path=roe_path_six)
    ma.variableToSignalSideExtraInfo("Xsd:sig_tree", {'constant(2)': 'mychiProb'}, path=roe_path_six)
    my_path.for_each('RestOfEvent', 'RestOfEvents', roe_path_six)
    
    # --- build Event Kinematics ---
    track_default = 'pt > 0.1 and thetaInCDCAcceptance and abs(dz) < 3 and dr < 0.5'
    gamma_default = 'E > 0.05 and thetaInCDCAcceptance and abs(clusterTiming) < 200'
    if(analysistype=="Jpsi"):
        track_default += " and sourceObjectIsInList(mu+:fromJpsimuchargedMuon)==0"
        gamma_default += " and sourceObjectIsInList(mu+:fromJpsimuchargedMuon)==0"
    ma.buildEventKinematics(custom_cuts=(track_default, gamma_default), path=my_path)
    
    # --- build Event shape ---
    ma.fillParticleList('pi+:evtshape_default', track_default, path=my_path)
    ma.fillParticleList('gamma:evtshape_default', gamma_default, path=my_path)
    ma.buildEventShape(inputListNames=['pi+:evtshape_default', 'gamma:evtshape_default'], path=my_path)
    
    # --- entire tree fit
    vertex.treeFit(list_name="Upsilon(4S):withoutneutrino", conf_level=-1, updateAllDaughters=False, path=my_path)
    
    # --- calculate Eecl + check unmatched photon ---
    roe_path_upsilon = basf2.Path()
    deadEndPath_upsilon = basf2.Path()
    ma.signalSideParticleFilter('Upsilon(4S):withoutneutrino', '', roe_path_upsilon, deadEndPath_upsilon)
    ma.fillSignalSideParticleList('Upsilon(4S):sig', '^Upsilon(4S):withoutneutrino', roe_path_upsilon)
    ma.fillParticleList("gamma:clean", cut=cluster_selection + " and [isInRestOfEvent == 1]", path=roe_path_upsilon)
    ma.fillParticleList("gamma:v133", cut=cluster_minimum_energy + " and [isInRestOfEvent == 1]", path=roe_path_upsilon)
    ma.fillParticleList("gamma:v200", cut=cluster_minimum_energy + " and [isInRestOfEvent == 1]", path=roe_path_upsilon)
    ma.fillParticleList("gamma:v200_800", cut=cluster_selection_v200_800 + " and [isInRestOfEvent == 1]", path=roe_path_upsilon)
    ma.fillParticleList("gamma:v200_825", cut=cluster_selection_v200_825 + " and [isInRestOfEvent == 1]", path=roe_path_upsilon)
    ma.fillParticleList("gamma:v200_850", cut=cluster_selection_v200_850 + " and [isInRestOfEvent == 1]", path=roe_path_upsilon)
    ma.fillParticleList("gamma:v200_875", cut=cluster_selection_v200_875 + " and [isInRestOfEvent == 1]", path=roe_path_upsilon)
    ma.fillParticleList("gamma:v200_900", cut=cluster_selection_v200_900 + " and [isInRestOfEvent == 1]", path=roe_path_upsilon)
    ma.fillParticleList("gamma:v200_925", cut=cluster_selection_v200_925 + " and [isInRestOfEvent == 1]", path=roe_path_upsilon)
    ma.fillParticleList("gamma:v200_950", cut=cluster_selection_v200_950 + " and [isInRestOfEvent == 1]", path=roe_path_upsilon)
    ma.fillParticleList("gamma:v200_975", cut=cluster_selection_v200_975 + " and [isInRestOfEvent == 1]", path=roe_path_upsilon)
    ma.fillParticleList("gamma:v200_025", cut=cluster_selection_v200_025 + " and [isInRestOfEvent == 1]", path=roe_path_upsilon)
    ma.fillParticleList("gamma:v200_050", cut=cluster_selection_v200_050 + " and [isInRestOfEvent == 1]", path=roe_path_upsilon)
    ma.fillParticleList("gamma:v200_075", cut=cluster_selection_v200_075 + " and [isInRestOfEvent == 1]", path=roe_path_upsilon)
    ma.fillParticleList("gamma:v200_100", cut=cluster_selection_v200_100 + " and [isInRestOfEvent == 1]", path=roe_path_upsilon)
    ma.applyCuts("gamma:v133", cut = cluster_selection_v133, path=roe_path_upsilon)
    ma.applyCuts("gamma:v200", cut = cluster_selection_v200, path=roe_path_upsilon)
    ma.matchMCTruth("gamma:clean", path=roe_path_upsilon)
    ma.matchMCTruth("gamma:v133", path=roe_path_upsilon)
    ma.matchMCTruth("gamma:v200", path=roe_path_upsilon)
    ma.matchMCTruth("gamma:v200_800", path=roe_path_upsilon)
    ma.matchMCTruth("gamma:v200_825", path=roe_path_upsilon)
    ma.matchMCTruth("gamma:v200_850", path=roe_path_upsilon)
    ma.matchMCTruth("gamma:v200_875", path=roe_path_upsilon)
    ma.matchMCTruth("gamma:v200_900", path=roe_path_upsilon)
    ma.matchMCTruth("gamma:v200_925", path=roe_path_upsilon)
    ma.matchMCTruth("gamma:v200_950", path=roe_path_upsilon)
    ma.matchMCTruth("gamma:v200_975", path=roe_path_upsilon)
    ma.matchMCTruth("gamma:v200_025", path=roe_path_upsilon)
    ma.matchMCTruth("gamma:v200_050", path=roe_path_upsilon)
    ma.matchMCTruth("gamma:v200_075", path=roe_path_upsilon)
    ma.matchMCTruth("gamma:v200_100", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:clean_matched", "gamma:clean", cut="mcPDG == 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:clean_unmatched", "gamma:clean", cut="mcPDG != 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v133_matched", "gamma:v133", cut="mcPDG == 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v133_unmatched", "gamma:v133", cut="mcPDG != 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_matched", "gamma:v200", cut="mcPDG == 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_unmatched", "gamma:v200", cut="mcPDG != 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_800_matched", "gamma:v200_800", cut="mcPDG == 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_800_unmatched", "gamma:v200_800", cut="mcPDG != 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_825_matched", "gamma:v200_825", cut="mcPDG == 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_825_unmatched", "gamma:v200_825", cut="mcPDG != 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_850_matched", "gamma:v200_850", cut="mcPDG == 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_850_unmatched", "gamma:v200_850", cut="mcPDG != 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_875_matched", "gamma:v200_875", cut="mcPDG == 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_875_unmatched", "gamma:v200_875", cut="mcPDG != 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_900_matched", "gamma:v200_900", cut="mcPDG == 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_900_unmatched", "gamma:v200_900", cut="mcPDG != 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_925_matched", "gamma:v200_925", cut="mcPDG == 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_925_unmatched", "gamma:v200_925", cut="mcPDG != 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_950_matched", "gamma:v200_950", cut="mcPDG == 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_950_unmatched", "gamma:v200_950", cut="mcPDG != 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_975_matched", "gamma:v200_975", cut="mcPDG == 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_975_unmatched", "gamma:v200_975", cut="mcPDG != 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_025_matched", "gamma:v200_025", cut="mcPDG == 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_025_unmatched", "gamma:v200_025", cut="mcPDG != 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_050_matched", "gamma:v200_050", cut="mcPDG == 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_050_unmatched", "gamma:v200_050", cut="mcPDG != 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_075_matched", "gamma:v200_075", cut="mcPDG == 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_075_unmatched", "gamma:v200_075", cut="mcPDG != 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_100_matched", "gamma:v200_100", cut="mcPDG == 22", path=roe_path_upsilon)
    ma.cutAndCopyList("gamma:v200_100_unmatched", "gamma:v200_100", cut="mcPDG != 22", path=roe_path_upsilon)
    ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'nParticlesInList(gamma:clean)': 'Ngamma'}, path=roe_path_upsilon)
    ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'totalECLEnergyOfParticlesInList(gamma:clean)': 'Eecl'}, path=roe_path_upsilon)
    ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'totalECLEnergyOfParticlesInList(gamma:clean_matched)': 'Eecl_matched'}, path=roe_path_upsilon)
    ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'totalECLEnergyOfParticlesInList(gamma:clean_unmatched)': 'Eecl_unmatched'}, path=roe_path_upsilon)
    ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'nParticlesInList(gamma:v133)': 'Ngammav133'}, path=roe_path_upsilon)
    ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'totalECLEnergyOfParticlesInList(gamma:v133)': 'Eeclv133'}, path=roe_path_upsilon)
    ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'totalECLEnergyOfParticlesInList(gamma:v133_matched)': 'Eeclv133_matched'}, path=roe_path_upsilon)
    ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'totalECLEnergyOfParticlesInList(gamma:v133_unmatched)': 'Eeclv133_unmatched'}, path=roe_path_upsilon)
    ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'nParticlesInList(gamma:v200)': 'Ngammav200'}, path=roe_path_upsilon)
    ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'nParticlesInList(gamma:v200_matched)': 'Ngammav200_matched'}, path=roe_path_upsilon)
    ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'nParticlesInList(gamma:v200_unmatched)': 'Ngammav200_unmatched'}, path=roe_path_upsilon)
    ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'totalECLEnergyOfParticlesInList(gamma:v200)': 'Eeclv200'}, path=roe_path_upsilon)
    ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'totalECLEnergyOfParticlesInList(gamma:v200_matched)': 'Eeclv200_matched'}, path=roe_path_upsilon)
    ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'totalECLEnergyOfParticlesInList(gamma:v200_unmatched)': 'Eeclv200_unmatched'}, path=roe_path_upsilon)
    ROEECL_scale = []
    for numbering in ["800", "825", "850", "875", "900", "925", "950", "975", "025", "050", "075", "100"]:
        ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'nParticlesInList(gamma:v200_' + numbering + ')': 'Ngammav200_' + numbering}, path=roe_path_upsilon)
        ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'nParticlesInList(gamma:v200_' + numbering + '_matched)': 'Ngammav200_' + numbering + '_matched'}, path=roe_path_upsilon)
        ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'nParticlesInList(gamma:v200_' + numbering + '_unmatched)': 'Ngammav200_' + numbering + '_unmatched'}, path=roe_path_upsilon)
        ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'totalECLEnergyOfParticlesInList(gamma:v200_' + numbering + ')': 'Eeclv200_' + numbering}, path=roe_path_upsilon)
        ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'totalECLEnergyOfParticlesInList(gamma:v200_' + numbering + '_matched)': 'Eeclv200_' + numbering + '_matched'}, path=roe_path_upsilon)
        ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'totalECLEnergyOfParticlesInList(gamma:v200_' + numbering + '_unmatched)': 'Eeclv200_' + numbering + '_unmatched'}, path=roe_path_upsilon)
    # get gamma from KL0
    ma.cutAndCopyList("gamma:v200_KL0", "gamma:v200", "mcPDG == 130", path=roe_path_upsilon)
    ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'nParticlesInList(gamma:v200_KL0)': 'Ngammav200_KL0'}, path=roe_path_upsilon)
    # ECL from tracks which do not satisfy loose track
    ma.fillParticleList("pi+:mismatch", cut="[[dr > 2] or [abs(dz) > 4]] and " + "[clusterE > 0] and " + "[isInRestOfEvent == 1]", path=roe_path_upsilon)
    ma.fillParticleList("pi+:mismatch_Ecut", cut="[[dr > 2] or [abs(dz) > 4]] and " + "[[clusterReg==1 and clusterE>0.10] or [clusterReg==2 and clusterE>0.06] or [clusterReg==3 and clusterE>0.15]] and " + "[isInRestOfEvent == 1]", path=roe_path_upsilon)
    ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'totalECLEnergyOfParticlesInList(pi+:mismatch)': 'Eecl_mismatch'}, path=roe_path_upsilon)
    ma.variableToSignalSideExtraInfo('Upsilon(4S):sig', {'totalECLEnergyOfParticlesInList(pi+:mismatch_Ecut)': 'Eecl_mismatch_Ecut'}, path=roe_path_upsilon)
    my_path.for_each('RestOfEvent', 'RestOfEvents', roe_path_upsilon)
        
    # --- MC truth ---
    ma.looseMCTruth(list_name = "B+:feiHadronic", path = my_path)
    ma.looseMCTruth(list_name = "B+:sig", path = my_path)
    ma.looseMCTruth(list_name = "B0:feiHadronic", path = my_path)
    ma.looseMCTruth(list_name = "B0:sig", path = my_path)
    ma.looseMCTruth(list_name = "Upsilon(4S):withoutneutrino", path = my_path)
    
    # --- decay string ---
    if(analysistype=="nominal"):
        my_path.add_module('ParticleMCDecayString', listName='Upsilon(4S):withoutneutrino', fileName=hashmapName)
    
    # get variables
    Kinematics = ["E", "InvM", "M", "p", "phi", "theta", "dM", "pt"]
    Kinematics_CMS = ["useCMSFrame(E)", "useCMSFrame(p)", "useCMSFrame(phi)", "useCMSFrame(theta)", "useCMSFrame(pt)"]
    Kinematics_RecoilRestFrame = ["useTagSideRecoilRestFrame(daughter(1,E),0)", "useTagSideRecoilRestFrame(daughter(1,p),0)", "useTagSideRecoilRestFrame(daughter(1,phi),0)", "useTagSideRecoilRestFrame(daughter(1,theta),0)"]
    Btag_cut = ["Mbc", "deltaE"]
    mcvar = ["mcPDG", "mcE", "mcP", "mcPhi", "mcTheta", "genParticleID", "isSignal", "isSignalAcceptMissingNeutrino", "isSignalAcceptMissingGamma"]
    loosemcvar = ["extraInfo(looseMCMotherPDG)", "extraInfo(looseMCWrongDaughterN)", "extraInfo(looseMCWrongDaughterPDG)", "extraInfo(looseMCWrongDaughterBiB)"]
    decayhash = ['extraInfo(DecayHash)', 'extraInfo(DecayHashExtended)']
    othervar = ["PDG", "extraInfo(decayModeID)"]
    continuumsup_vars = ["R2", "thrustBm", "thrustOm", "cosTBTO", "cosTBz", "KSFWVariables(et)", "KSFWVariables(mm2)", "KSFWVariables(hso00)", "KSFWVariables(hso01)", "KSFWVariables(hso02)", "KSFWVariables(hso03)", "KSFWVariables(hso04)", "KSFWVariables(hso10)", "KSFWVariables(hso12)", "KSFWVariables(hso14)", "KSFWVariables(hso20)", "KSFWVariables(hso22)", "KSFWVariables(hso24)", "KSFWVariables(hoo0)", "KSFWVariables(hoo1)", "KSFWVariables(hoo2)", "KSFWVariables(hoo3)", "KSFWVariables(hoo4)", "CleoConeCS(1)", "CleoConeCS(2)", "CleoConeCS(3)", "CleoConeCS(4)", "CleoConeCS(5)", "CleoConeCS(6)", "CleoConeCS(7)", "CleoConeCS(8)", "CleoConeCS(9)"]
    EvtKinematics = ["missingMomentumOfEvent", "missingMomentumOfEvent_theta", "missingMass2OfEvent", "visibleEnergyOfEventCMS", "missingEnergyOfEventCMS"]
    simpleDvetovar = ["daughter(0, extraInfo(nDc))", "daughter(0, extraInfo(Dc_pValue_med))", "daughter(0, extraInfo(Dc_pValue_std))", "daughter(0, extraInfo(Dcsimpleveto_chiProb))", "daughter(0, extraInfo(Dcsimpleveto_dr))", "daughter(0, extraInfo(Dcsimpleveto_dz))", "daughter(0, extraInfo(Dcsimpleveto_M))", "daughter(0, extraInfo(nD0))", "daughter(0, extraInfo(D0_pValue_med))", "daughter(0, extraInfo(D0_pValue_std))", "daughter(0, extraInfo(D0simpleveto_chiProb))", "daughter(0, extraInfo(D0simpleveto_dr))", "daughter(0, extraInfo(D0simpleveto_dz))", "daughter(0, extraInfo(D0simpleveto_M))"]
    simpleDvetovar = simpleDvetovar + ["daughter(0, extraInfo(nDc_noDCS))", "daughter(0, extraInfo(Dc_pValue_med_noDCS))", "daughter(0, extraInfo(Dc_pValue_std_noDCS))", "daughter(0, extraInfo(Dcsimpleveto_chiProb_noDCS))", "daughter(0, extraInfo(Dcsimpleveto_dr_noDCS))", "daughter(0, extraInfo(Dcsimpleveto_dz_noDCS))", "daughter(0, extraInfo(Dcsimpleveto_M_noDCS))"]
    simpleDvetovar = simpleDvetovar + ["daughter(0, extraInfo(nDc_yespizero))", "daughter(0, extraInfo(Dc_pValue_med_yespizero))", "daughter(0, extraInfo(Dc_pValue_std_yespizero))", "daughter(0, extraInfo(Dcsimpleveto_chiProb_yespizero))", "daughter(0, extraInfo(Dcsimpleveto_dr_yespizero))", "daughter(0, extraInfo(Dcsimpleveto_dz_yespizero))", "daughter(0, extraInfo(Dcsimpleveto_M_yespizero))", "daughter(0, extraInfo(nD0_yespizero))", "daughter(0, extraInfo(D0_pValue_med_yespizero))", "daughter(0, extraInfo(D0_pValue_std_yespizero))", "daughter(0, extraInfo(D0simpleveto_chiProb_yespizero))", "daughter(0, extraInfo(D0simpleveto_dr_yespizero))", "daughter(0, extraInfo(D0simpleveto_dz_yespizero))", "daughter(0, extraInfo(D0simpleveto_M_yespizero))"]
    vertexXs = ["daughter(0, extraInfo(mydr))", "daughter(0, extraInfo(mydz))", "daughter(0, extraInfo(mychiProb))"]
    leptonInfo = ["nParticlesInList(e+:ElectronFBDT)", "nParticlesInList(mu+:MuonFBDT)", "nParticlesInList(e+:ElectronFBDT_loose)", "nParticlesInList(mu+:MuonFBDT_loose)", "nParticlesInList(e+:ElectronFBDT_tight)", "nParticlesInList(mu+:MuonFBDT_tight)"]
    PIDsyst = []
    index = 0
    pbins = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
    cosThetabins = [-0.866, -0.682, -0.4226, -0.1045, 0.225, 0.5, 0.766, 0.8829, 0.9563]
    for i_p in range(len(pbins) - 1):
        for i_cosTheta in range(len(cosThetabins) - 1):
            PIDsyst.append("daughter(0,extraInfo(nKtruebin"+str(index)+"))")
            PIDsyst.append("daughter(0,extraInfo(nKmisbin"+str(index)+"))")
            PIDsyst.append("daughter(0,extraInfo(npitruebin"+str(index)+"))")
            PIDsyst.append("daughter(0,extraInfo(npimisbin"+str(index)+"))")
            index = index + 1
    PIDsyst.append("daughter(0,extraInfo(nKtruebin"+str(index)+"))")
    PIDsyst.append("daughter(0,extraInfo(nKmisbin"+str(index)+"))")
    PIDsyst.append("daughter(0,extraInfo(npitruebin"+str(index)+"))")
    PIDsyst.append("daughter(0,extraInfo(npimisbin"+str(index)+"))")
    
    index = 0
    pbins_e = [0.4, 0.5, 1.0, 1.5, 2, 2.5, 3]
    thetabins_e = [0.22, 0.56, 1.13, 1.57, 1.88, 2.23, 2.71]
    for i_p in range(len(pbins_e) - 1):
        for i_theta in range(len(thetabins_e) - 1):
            PIDsyst.append("daughter(0,extraInfo(nKfakeEbin_n"+str(index)+"))")
            PIDsyst.append("daughter(0,extraInfo(nKfakeEbin_p"+str(index)+"))")
            PIDsyst.append("daughter(0,extraInfo(npifakeEbin_n"+str(index)+"))")
            PIDsyst.append("daughter(0,extraInfo(npifakeEbin_p"+str(index)+"))")
            index = index + 1
    PIDsyst.append("daughter(0,extraInfo(nKfakeEbin_n"+str(index)+"))")
    PIDsyst.append("daughter(0,extraInfo(nKfakeEbin_p"+str(index)+"))")
    PIDsyst.append("daughter(0,extraInfo(npifakeEbin_n"+str(index)+"))")
    PIDsyst.append("daughter(0,extraInfo(npifakeEbin_p"+str(index)+"))")
    
    index = 0
    pbins_mu = [0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5]
    thetabins_mu = [0.4, 0.64, 0.82, 1.16, 1.46, 1.78, 2.13, 2.22, 2.6]
    for i_p in range(len(pbins_mu) - 1):
        for i_theta in range(len(thetabins_mu) - 1):
            PIDsyst.append("daughter(0,extraInfo(nKfakeMUbin_n"+str(index)+"))")
            PIDsyst.append("daughter(0,extraInfo(nKfakeMUbin_p"+str(index)+"))")
            PIDsyst.append("daughter(0,extraInfo(npifakeMUbin_n"+str(index)+"))")
            PIDsyst.append("daughter(0,extraInfo(npifakeMUbin_p"+str(index)+"))")
            index = index + 1
    PIDsyst.append("daughter(0,extraInfo(nKfakeMUbin_n"+str(index)+"))")
    PIDsyst.append("daughter(0,extraInfo(nKfakeMUbin_p"+str(index)+"))")
    PIDsyst.append("daughter(0,extraInfo(npifakeMUbin_n"+str(index)+"))")
    PIDsyst.append("daughter(0,extraInfo(npifakeMUbin_p"+str(index)+"))")
    
    pi0syst = []
    index = 0
    pi0_pbins = [0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0]
    for i_p in range(len(pi0_pbins) - 1):
        pi0syst.append("daughter(0,extraInfo(npi0bin"+str(index)+"))")
        index = index + 1
    pi0syst.append("daughter(0,extraInfo(npi0bin"+str(index)+"))")
    
    ROEECL_scale = []
    for numbering in ["800", "825", "850", "875", "900", "925", "950", "975", "025", "050", "075", "100"]:
        ROEECL_scale.append("extraInfo(Ngammav200_" + numbering + ")")
        ROEECL_scale.append("extraInfo(Ngammav200_" + numbering + "_matched)")
        ROEECL_scale.append("extraInfo(Ngammav200_" + numbering + "_unmatched)")
        ROEECL_scale.append("extraInfo(Eeclv200_" + numbering + ")")
        ROEECL_scale.append("extraInfo(Eeclv200_" + numbering + "_matched)")
        ROEECL_scale.append("extraInfo(Eeclv200_" + numbering + "_unmatched)")
    
    ROEECL = ["extraInfo(Ngamma)", "extraInfo(Eecl)", "extraInfo(Eecl_matched)", "extraInfo(Eecl_unmatched)", "extraInfo(Ngammav133)", "extraInfo(Eeclv133)", "extraInfo(Eeclv133_matched)", "extraInfo(Eeclv133_unmatched)", "extraInfo(Ngammav200)", "extraInfo(Ngammav200_matched)", "extraInfo(Ngammav200_unmatched)", "extraInfo(Eeclv200)", "extraInfo(Eeclv200_matched)", "extraInfo(Eeclv200_unmatched)", "extraInfo(Ngammav200_KL0)", "extraInfo(Eecl_mismatch)", "extraInfo(Eecl_mismatch_Ecut)"] + ROEECL_scale
    
    DTOKL0 = []
    KNN = []
    XNN = []
    KKLKL = []
    XKLKL = []
    if((args.sampletype=="MC_qq") or (args.sampletype=="MC_BB") or (args.sampletype=="MC_SIGNAL")):
        DTOKL0.append("nParticlesInList(D+:DecayIntoKL0)")
        DTOKL0.append("nParticlesInList(D0:DecayIntoKL0)")
    
        KNN.append("nParticlesInList(B+:Knn)")
        KNN.append("invMassInLists(n0:Knn)")
        KNN.append("nParticlesInList(B+:Kstarnn)")
        KNN.append("invMassInLists(n0:Kstarnn)")
        KNN.append("nParticlesInList(B0:K0nn)")
        KNN.append("invMassInLists(n0:K0nn)")
        KNN.append("nParticlesInList(B0:Kstar0nn)")
        KNN.append("invMassInLists(n0:Kstar0nn)")

        XNN.append("nParticlesInList(B+:Xnn)")
        XNN.append("nParticlesInList(B0:Xnn)")
        XNN.append("invMassInLists(n0:Xnn)")

        KKLKL.append("nParticlesInList(B+:KpKLKL_all)")
        KKLKL.append("nParticlesInList(B+:KpKLKL_NR)")
        KKLKL.append("averageValueInList(B+:KpKLKL_NR, daughterInvariantMass(0, 1))")
        KKLKL.append("averageValueInList(B+:KpKLKL_NR, daughterInvariantMass(0, 2))")
        KKLKL.append("averageValueInList(B+:KpKLKL_NR, daughterInvariantMass(1, 2))")
    
        KKLKL.append("nParticlesInList(B0:KSKLKL_all)")
        KKLKL.append("nParticlesInList(B0:KSKLKL_NR)")
        KKLKL.append("averageValueInList(B0:KSKLKL_NR, daughterInvariantMass(0, 1))")
        KKLKL.append("averageValueInList(B0:KSKLKL_NR, daughterInvariantMass(0, 2))")
        KKLKL.append("averageValueInList(B0:KSKLKL_NR, daughterInvariantMass(1, 2))")

        XKLKL.append("nParticlesInList(K_L0:XKLKL)")
        XKLKL.append("averageValueInList(K_L0:XKLKL, E)")
        XKLKL.append("averageValueInList(K_L0:XKLKL, px)")
        XKLKL.append("averageValueInList(K_L0:XKLKL, py)")
        XKLKL.append("averageValueInList(K_L0:XKLKL, pz)")
        XKLKL.append("averageValueInList(K_L0:XKLKL_1st, E)")
        XKLKL.append("averageValueInList(K_L0:XKLKL_1st, px)")
        XKLKL.append("averageValueInList(K_L0:XKLKL_1st, py)")
        XKLKL.append("averageValueInList(K_L0:XKLKL_1st, pz)")
        XKLKL.append("averageValueInList(K_L0:XKLKL_2nd, E)")
        XKLKL.append("averageValueInList(K_L0:XKLKL_2nd, px)")
        XKLKL.append("averageValueInList(K_L0:XKLKL_2nd, py)")
        XKLKL.append("averageValueInList(K_L0:XKLKL_2nd, pz)")
        XKLKL.append("nParticlesInList(K_S0:XKLKL)")
        XKLKL.append("averageValueInList(K_S0:XKLKL, E)")
        XKLKL.append("averageValueInList(K_S0:XKLKL, px)")
        XKLKL.append("averageValueInList(K_S0:XKLKL, py)")
        XKLKL.append("averageValueInList(K_S0:XKLKL, pz)")

    else:
        va.variables.addAlias("nParticlesInList__boD__pl__clDecayIntoKL0__bc", "constant(0)")
        va.variables.addAlias("nParticlesInList__boD0__clDecayIntoKL0__bc", "constant(0)")
        DTOKL0.append("nParticlesInList__boD__pl__clDecayIntoKL0__bc")
        DTOKL0.append("nParticlesInList__boD0__clDecayIntoKL0__bc")
    
        va.variables.addAlias("nParticlesInList__boB__pl__clKnn__bc", "constant(0)")
        va.variables.addAlias("invMassInLists__bon0__clKnn__bc", "constant(0)")
        va.variables.addAlias("nParticlesInList__boB__pl__clKstarnn__bc", "constant(0)")
        va.variables.addAlias("invMassInLists__bon0__clKstarnn__bc", "constant(0)")
        va.variables.addAlias("nParticlesInList__boB0__clK0nn__bc", "constant(0)")
        va.variables.addAlias("invMassInLists__bon0__clK0nn__bc", "constant(0)")
        va.variables.addAlias("nParticlesInList__boB0__clKstar0nn__bc", "constant(0)")
        va.variables.addAlias("invMassInLists__bon0__clKstar0nn__bc", "constant(0)")
        KNN.append("nParticlesInList__boB__pl__clKnn__bc")
        KNN.append("invMassInLists__bon0__clKnn__bc")
        KNN.append("nParticlesInList__boB__pl__clKstarnn__bc")
        KNN.append("invMassInLists__bon0__clKstarnn__bc")
        KNN.append("nParticlesInList__boB0__clK0nn__bc")
        KNN.append("invMassInLists__bon0__clK0nn__bc")
        KNN.append("nParticlesInList__boB0__clKstar0nn__bc")
        KNN.append("invMassInLists__bon0__clKstar0nn__bc")

        va.variables.addAlias("nParticlesInList__boB__pl__clXnn__bc", "constant(0)")
        va.variables.addAlias("nParticlesInList__boB0__clXnn__bc", "constant(0)")
        va.variables.addAlias("invMassInLists__bon0__clXnn__bc", "constant(0)")
        XNN.append("nParticlesInList__boB__pl__clXnn__bc")
        XNN.append("nParticlesInList__boB0__clXnn__bc")
        XNN.append("invMassInLists__bon0__clXnn__bc")

        va.variables.addAlias("nParticlesInList__boB__pl__clKpKLKL_all__bc", "constant(0)")
        va.variables.addAlias("nParticlesInList__boB__pl__clKpKLKL_NR__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", "constant(0)")
        KKLKL.append("nParticlesInList__boB__pl__clKpKLKL_all__bc")
        KKLKL.append("nParticlesInList__boB__pl__clKpKLKL_NR__bc")
        KKLKL.append("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc")
        KKLKL.append("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc")
        KKLKL.append("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc")
    
        va.variables.addAlias("nParticlesInList__boB0__clKSKLKL_all__bc", "constant(0)")
        va.variables.addAlias("nParticlesInList__boB0__clKSKLKL_NR__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", "constant(0)")
        KKLKL.append("nParticlesInList__boB0__clKSKLKL_all__bc")
        KKLKL.append("nParticlesInList__boB0__clKSKLKL_NR__bc")
        KKLKL.append("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc")
        KKLKL.append("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc")
        KKLKL.append("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc")

        va.variables.addAlias("nParticlesInList__boK_L0__clXKLKL__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boK_L0__clXKLKL__cm__spE__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boK_L0__clXKLKL__cm__sppx__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boK_L0__clXKLKL__cm__sppy__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boK_L0__clXKLKL__cm__sppz__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boK_L0__clXKLKL_1st__cm__spE__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boK_L0__clXKLKL_1st__cm__sppx__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boK_L0__clXKLKL_1st__cm__sppy__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boK_L0__clXKLKL_1st__cm__sppz__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boK_L0__clXKLKL_2nd__cm__spE__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppx__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppy__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppz__bc", "constant(0)")
        va.variables.addAlias("nParticlesInList__boK_S0__clXKLKL__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boK_S0__clXKLKL__cm__spE__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boK_S0__clXKLKL__cm__sppx__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boK_S0__clXKLKL__cm__sppy__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boK_S0__clXKLKL__cm__sppz__bc", "constant(0)")
        XKLKL.append("nParticlesInList__boK_L0__clXKLKL__bc")
        XKLKL.append("averageValueInList__boK_L0__clXKLKL__cm__spE__bc")
        XKLKL.append("averageValueInList__boK_L0__clXKLKL__cm__sppx__bc")
        XKLKL.append("averageValueInList__boK_L0__clXKLKL__cm__sppy__bc")
        XKLKL.append("averageValueInList__boK_L0__clXKLKL__cm__sppz__bc")
        XKLKL.append("averageValueInList__boK_L0__clXKLKL_1st__cm__spE__bc")
        XKLKL.append("averageValueInList__boK_L0__clXKLKL_1st__cm__sppx__bc")
        XKLKL.append("averageValueInList__boK_L0__clXKLKL_1st__cm__sppy__bc")
        XKLKL.append("averageValueInList__boK_L0__clXKLKL_1st__cm__sppz__bc")
        XKLKL.append("averageValueInList__boK_L0__clXKLKL_2nd__cm__spE__bc")
        XKLKL.append("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppx__bc")
        XKLKL.append("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppy__bc")
        XKLKL.append("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppz__bc")
        XKLKL.append("nParticlesInList__boK_S0__clXKLKL__bc")
        XKLKL.append("averageValueInList__boK_S0__clXKLKL__cm__spE__bc")
        XKLKL.append("averageValueInList__boK_S0__clXKLKL__cm__sppx__bc")
        XKLKL.append("averageValueInList__boK_S0__clXKLKL__cm__sppy__bc")
        XKLKL.append("averageValueInList__boK_S0__clXKLKL__cm__sppz__bc")

    Btag_vars = vu.create_aliases(list_of_variables = Kinematics + Btag_cut + Kinematics_CMS + othervar + mcvar + loosemcvar + continuumsup_vars + ["chiProb", "extraInfo(SignalProbability)", "dr", "dz"], wrapper = "daughter(0,{variable})",prefix="Btag")
    
    Bsig_vars = vu.create_aliases(list_of_variables = Kinematics + Btag_cut + Kinematics_CMS + simpleDvetovar + othervar + mcvar + loosemcvar + continuumsup_vars + vertexXs + [ "daughter(0, M)", "daughter(0, extraInfo(decayModeID))", "daughter(0, extraInfo(KS0_3D_distance))"] + PIDsyst + pi0syst + ["daughter(0, extraInfo(nKslow1))", "daughter(0, extraInfo(nKslow2))", "daughter(0, extraInfo(nKslow3))", "daughter(0, extraInfo(nPislow1))", "daughter(0, extraInfo(nPislow2))", "daughter(0, extraInfo(nPislow3))", "dr", "dz", "daughter(0, daughter(0, dr))", "daughter(0, daughter(0, dz))"], wrapper = "daughter(1,{variable})", prefix="Bsig")
    
    U_vars = Kinematics + Kinematics_CMS + Kinematics_RecoilRestFrame + EvtKinematics + decayhash + othervar + mcvar + loosemcvar + vc.event_shape + leptonInfo + ["extraInfo(Upsilon_rank)",  "nROE_ECLClusters(cleanMask)", "nROE_NeutralECLClusters(cleanMask)", "nROE_KLMClusters", "nROE_Tracks(cleanMask)", "roeEextra(cleanMask)", "roeNeextra(cleanMask)", "useCMSFrame(roeNeextra(cleanMask))", "nROE_ParticlesInList(K_S0:myKaonshort)", "nROE_ParticlesInList(pi0:myneutralPion)", "nROE_ParticlesInList(gamma:mygamma)", "nRemainingTracksInEvent", "MsquaredBsig_op0", "MsquaredBsig_op1", "MsquaredBsig_op2", "MsquaredBsig_op3", "MsquaredBsig_op4", "MsquaredBsig_op7", "roeP(cleanMask)", "roeM(cleanMask)", "roePTheta(cleanMask)", "qsquared", "extraInfo(chiSquared)", "extraInfo(ndf)", "chiProb", "dr", "dz", "beamE", "nROE_Tracks(looseMask)", "Ecms", "m2RecoilSignalSide"] + ROEECL + DTOKL0 + KNN + XNN + KKLKL + XKLKL
    if(analysistype=="Jpsi"):
        U_vars = U_vars + ["nParticlesInList(mu+:fromJpsimuchargedMuon)", "nParticlesInList(J/psi:temp_BCS)", "nParticlesInList(Upsilon(4S):temp_withoutneutrino)"]
    elif(analysistype=="ppbar"):
        U_vars = U_vars + ["daughter(1, daughterInvariantMass(1, 2))"]
    elif(analysistype=="KSKS"):
        U_vars = U_vars + ["daughter(1,extraInfo(N_KS_KSKS))", "daughter(1,extraInfo(avg_KS_KSKS_E))", "daughter(1,extraInfo(avg_KS_KSKS_px))", "daughter(1,extraInfo(avg_KS_KSKS_py))", "daughter(1,extraInfo(avg_KS_KSKS_pz))", "daughter(1,extraInfo(KS_KSKS_1st_E))", "daughter(1,extraInfo(KS_KSKS_1st_px))", "daughter(1,extraInfo(KS_KSKS_1st_py))", "daughter(1,extraInfo(KS_KSKS_1st_pz))", "daughter(1,extraInfo(KS_KSKS_2nd_E))", "daughter(1,extraInfo(KS_KSKS_2nd_px))", "daughter(1,extraInfo(KS_KSKS_2nd_py))", "daughter(1,extraInfo(KS_KSKS_2nd_pz))"]


    # Ntuple output
    ma.variablesToNtuple(decayString="Upsilon(4S):withoutneutrino",variables=Btag_vars,filename=output_file,treename="Btag",path=my_path)
    ma.variablesToNtuple(decayString="Upsilon(4S):withoutneutrino",variables=Bsig_vars,filename=output_file,treename="Bsig",path=my_path)
    ma.variablesToNtuple(decayString="Upsilon(4S):withoutneutrino",variables=U_vars,filename=output_file,treename="Upsilon",path=my_path)
    
    if((args.sampletype=="MC_SIGNAL") and (analysistype=="nominal")):
        # ==========================================MC data for calculate efficiency====================================
        # Load particles from MCParticle at first
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
        ma.fillParticleListFromMC("Upsilon(4S):PrimaryMC", cut = "mcPrimary", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
    
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

        # detail investigation for K pi pi
        ma.reconstructMCDecay('Xsu:Kpipi_ch1 =direct=> rho0:PrimaryMC K+:PrimaryMC', cut='', dmID = 1, path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipi_ch2 =direct=> pi+:PrimaryMC K*0:PrimaryMC', cut='', dmID = 2, path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipi_ch3 =direct=> K+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC', cut='', dmID = 3, path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipi_ch4 =direct=> K+:PrimaryMC omega:PrimaryMC', cut='', dmID = 4, path=my_path)
        ma.copyLists(outputListName="Xsu:Kppimpip", inputListNames=["Xsu:Kpipi_ch1", "Xsu:Kpipi_ch2", "Xsu:Kpipi_ch3", "Xsu:Kpipi_ch4"], path=my_path)

        ma.reconstructMCDecay('Xsu:Kpipi_ch5 =direct=> rho+:PrimaryMC K0:PrimaryMC', cut='', dmID = 5, path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipi_ch6 =direct=> pi+:PrimaryMC K*0:PrimaryMC', cut='', dmID = 6, path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipi_ch7 =direct=> K0:PrimaryMC pi0:PrimaryMC pi+:PrimaryMC', cut='', dmID = 7, path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipi_ch8 =direct=> pi0:PrimaryMC K*+:PrimaryMC', cut='', dmID = 8, path=my_path)
        ma.copyLists(outputListName="Xsu:Kzpizpip", inputListNames=["Xsu:Kpipi_ch5", "Xsu:Kpipi_ch6", "Xsu:Kpipi_ch7", "Xsu:Kpipi_ch8"], path=my_path)

        ma.reconstructMCDecay('Xsd:Kpipi_ch1 =direct=> rho-:PrimaryMC K+:PrimaryMC', cut='', dmID = 1, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipi_ch2 =direct=> pi-:PrimaryMC K*+:PrimaryMC', cut='', dmID = 2, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipi_ch3 =direct=> K+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC', cut='', dmID = 3, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipi_ch4 =direct=> pi0:PrimaryMC K*0:PrimaryMC', cut='', dmID = 4, path=my_path)
        ma.copyLists(outputListName="Xsd:Kppimpiz", inputListNames=["Xsd:Kpipi_ch1", "Xsd:Kpipi_ch2", "Xsd:Kpipi_ch3", "Xsd:Kpipi_ch4"], path=my_path)

        ma.reconstructMCDecay('Xsd:Kpipi_ch5 =direct=> rho0:PrimaryMC K0:PrimaryMC', cut='', dmID = 5, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipi_ch6 =direct=> pi-:PrimaryMC K*+:PrimaryMC', cut='', dmID = 6, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipi_ch7 =direct=> pi-:PrimaryMC pi+:PrimaryMC K0:PrimaryMC', cut='', dmID = 7, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipi_ch8 =direct=> omega:PrimaryMC K0:PrimaryMC', cut='', dmID = 8, path=my_path)
        ma.copyLists(outputListName="Xsd:Kzpippim", inputListNames=["Xsd:Kpipi_ch5", "Xsd:Kpipi_ch6", "Xsd:Kpipi_ch7", "Xsd:Kpipi_ch8"], path=my_path)

        # detail investigation for K pi pi pi
        ma.reconstructMCDecay('Xsu:Kpipipi_ch1 =direct=> rho+:PrimaryMC pi-:PrimaryMC K+:PrimaryMC', cut='',dmID = 1,  path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipipi_ch2 =direct=> K+:PrimaryMC eta:PrimaryMC', cut='', dmID = 2, path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipipi_ch3 =direct=> omega:PrimaryMC K+:PrimaryMC', cut='', dmID = 3, path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipipi_ch4 =direct=> pi-:PrimaryMC pi+:PrimaryMC pi0:PrimaryMC K+:PrimaryMC', cut='', dmID = 4, path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipipi_ch5 =direct=> K+:PrimaryMC phi:PrimaryMC', cut='', dmID = 5, path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipipi_ch6 =direct=> rho+:PrimaryMC K*0:PrimaryMC', cut='', dmID = 6, path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipipi_ch7 =direct=> pi+:PrimaryMC pi0:PrimaryMC K*0:PrimaryMC', cut='', dmID = 7, path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipipi_ch8 =direct=> rho-:PrimaryMC pi+:PrimaryMC K+:PrimaryMC', cut='',dmID = 8, path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipipi_ch9 =direct=> rho0:PrimaryMC pi0:PrimaryMC K+:PrimaryMC', cut='',dmID = 9,  path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipipi_ch10 =direct=> rho0:PrimaryMC K*+:PrimaryMC', cut='',dmID = 10,  path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipipi_ch11 =direct=> pi-:PrimaryMC pi+:PrimaryMC K*+:PrimaryMC', cut='',dmID = 11,  path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipipi_ch12 =direct=> omega:PrimaryMC K+:PrimaryMC pi0:PrimaryMC', cut='', dmID = 12, path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipipi_ch13 =direct=> omega:PrimaryMC K*+:PrimaryMC', cut='', dmID = 13, path=my_path)
        ma.copyLists(outputListName="Xsu:Kppippippiz", inputListNames=["Xsu:Kpipipi_ch1", "Xsu:Kpipipi_ch2", "Xsu:Kpipipi_ch3", "Xsu:Kpipipi_ch4", "Xsu:Kpipipi_ch5", "Xsu:Kpipipi_ch6", "Xsu:Kpipipi_ch7", "Xsu:Kpipipi_ch8", "Xsu:Kpipipi_ch9", "Xsu:Kpipipi_ch10", "Xsu:Kpipipi_ch11", "Xsu:Kpipipi_ch12", "Xsu:Kpipipi_ch13"], path=my_path)

        ma.reconstructMCDecay('Xsu:Kpipipi_ch14 =direct=> K*+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC', cut='',dmID = 14,  path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipipi_ch15 =direct=> K0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC', cut='', dmID = 15, path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipipi_ch16 =direct=> K*+:PrimaryMC rho0:PrimaryMC', cut='', dmID = 16, path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipipi_ch17 =direct=> K*+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC', cut='', dmID = 17, path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipipi_ch18 =direct=> K0:PrimaryMC pi+:PrimaryMC rho0:PrimaryMC', cut='', dmID = 18, path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipipi_ch19 =direct=> K*+:PrimaryMC omega:PrimaryMC', cut='', dmID = 19, path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipipi_ch20 =direct=> K0:PrimaryMC pi+:PrimaryMC omega:PrimaryMC', cut='', dmID = 20, path=my_path)
        ma.copyLists(outputListName="Xsu:Kzpippimpip", inputListNames=["Xsu:Kpipipi_ch14", "Xsu:Kpipipi_ch15", "Xsu:Kpipipi_ch16", "Xsu:Kpipipi_ch17", "Xsu:Kpipipi_ch18", "Xsu:Kpipipi_ch19", "Xsu:Kpipipi_ch20"], path=my_path)
        
        ma.reconstructMCDecay('Xsu:Kpipipi_ch21 =direct=> pi0:PrimaryMC rho+:PrimaryMC K0:PrimaryMC', cut='',dmID = 21,  path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipipi_ch22 =direct=> K*0:PrimaryMC pi+:PrimaryMC pi0:PrimaryMC', cut='', dmID = 22, path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipipi_ch23 =direct=> K*0:PrimaryMC rho+:PrimaryMC', cut='', dmID = 23, path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipipi_ch24 =direct=> K*+:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 24, path=my_path)
        ma.reconstructMCDecay('Xsu:Kpipipi_ch25 =direct=> K0:PrimaryMC pi+:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 25, path=my_path)
        ma.copyLists(outputListName="Xsu:Kzpippizpiz", inputListNames=["Xsu:Kpipipi_ch21", "Xsu:Kpipipi_ch22", "Xsu:Kpipipi_ch23", "Xsu:Kpipipi_ch24", "Xsu:Kpipipi_ch25"], path=my_path)

        ma.reconstructMCDecay('Xsd:Kpipipi_ch1 =direct=> rho0:PrimaryMC K+:PrimaryMC pi-:PrimaryMC', cut='',dmID = 1,  path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipipi_ch2 =direct=> K*0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC', cut='', dmID = 2, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipipi_ch3 =direct=> K*0:PrimaryMC rho0:PrimaryMC', cut='', dmID = 3, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipipi_ch4 =direct=> pi-:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC K+:PrimaryMC', cut='', dmID = 4, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipipi_ch5 =direct=> omega:PrimaryMC K*0:PrimaryMC', cut='', dmID = 5, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipipi_ch6 =direct=> omega:PrimaryMC K+:PrimaryMC pi-:PrimaryMC', cut='', dmID = 6, path=my_path)
        ma.copyLists(outputListName="Xsd:Kppimpippim", inputListNames=["Xsd:Kpipipi_ch1", "Xsd:Kpipipi_ch2", "Xsd:Kpipipi_ch3", "Xsd:Kpipipi_ch4", "Xsd:Kpipipi_ch5", "Xsd:Kpipipi_ch6"], path=my_path)

        ma.reconstructMCDecay('Xsd:Kpipipi_ch7 =direct=> K0:PrimaryMC omega:PrimaryMC', cut='', dmID = 7, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipipi_ch8 =direct=> K0:PrimaryMC pi0:PrimaryMC rho0:PrimaryMC', cut='', dmID = 8, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipipi_ch9 =direct=> K0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC', cut='', dmID = 9, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipipi_ch10 =direct=> K*0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC', cut='', dmID = 10, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipipi_ch11 =direct=> K0:PrimaryMC pi-:PrimaryMC rho+:PrimaryMC', cut='', dmID = 11, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipipi_ch12 =direct=> K*+:PrimaryMC rho-:PrimaryMC', cut='', dmID = 12, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipipi_ch13 =direct=> K0:PrimaryMC phi:PrimaryMC', cut='', dmID = 13, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipipi_ch14 =direct=> K*+:PrimaryMC pi0:PrimaryMC pi-:PrimaryMC', cut='', dmID = 14, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipipi_ch15 =direct=> K*0:PrimaryMC rho0:PrimaryMC', cut='', dmID = 15, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipipi_ch16 =direct=> K0:PrimaryMC eta:PrimaryMC', cut='', dmID = 16, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipipi_ch17 =direct=> K0:PrimaryMC pi+:PrimaryMC rho-:PrimaryMC', cut='', dmID = 17, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipipi_ch18 =direct=> K0:PrimaryMC pi0:PrimaryMC omega:PrimaryMC', cut='', dmID = 18, path=my_path)
        ma.copyLists(outputListName="Xsd:Kzpippimpiz", inputListNames=["Xsd:Kpipipi_ch7", "Xsd:Kpipipi_ch8", "Xsd:Kpipipi_ch9", "Xsd:Kpipipi_ch10", "Xsd:Kpipipi_ch11", "Xsd:Kpipipi_ch12", "Xsd:Kpipipi_ch13", "Xsd:Kpipipi_ch14", "Xsd:Kpipipi_ch15", "Xsd:Kpipipi_ch16", "Xsd:Kpipipi_ch17", "Xsd:Kpipipi_ch18"], path=my_path)

        ma.reconstructMCDecay('Xsd:Kpipipi_ch19 =direct=> K+:PrimaryMC rho-:PrimaryMC pi0:PrimaryMC', cut='', dmID = 19, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipipi_ch20 =direct=> K*+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC', cut='', dmID = 20, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipipi_ch21 =direct=> K+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 21, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipipi_ch22 =direct=> K*+:PrimaryMC rho-:PrimaryMC', cut='', dmID = 22, path=my_path)
        ma.reconstructMCDecay('Xsd:Kpipipi_ch23 =direct=> K*0:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 23, path=my_path)
        ma.copyLists(outputListName="Xsd:Kppimpizpiz", inputListNames=["Xsd:Kpipipi_ch19", "Xsd:Kpipipi_ch20", "Xsd:Kpipipi_ch21", "Xsd:Kpipipi_ch22", "Xsd:Kpipipi_ch23"], path=my_path)

        ma.reconstructMCDecay('Xsu:1270 -> K_1+:PrimaryMC ...', cut='',dmID = 0,  path=my_path)
        ma.reconstructMCDecay('Xsu:1430 -> K_2*+:PrimaryMC ...', cut='',dmID = 0,  path=my_path)

        variables = ["nParticlesInList(Xsu:MCch4)", "nParticlesInList(Xsu:MCch5)", "nParticlesInList(Xsu:MCch6)", "nParticlesInList(Xsu:MCch7)", "nParticlesInList(Xsu:MCch11)"]
        variables = variables + ["nParticlesInList(Xsd:MCch19)", "nParticlesInList(Xsd:MCch20)", "nParticlesInList(Xsd:MCch21)", "nParticlesInList(Xsd:MCch22)", "nParticlesInList(Xsd:MCch26)"]

        for index in range(1,9):
            variables = variables + ["nParticlesInList(Xsu:Kpipi_ch" + str(index) + ")"]
        for index in range(1,26):
            variables = variables + ["nParticlesInList(Xsu:Kpipipi_ch" + str(index) + ")"]
        for index in range(1,9):
            variables = variables + ["nParticlesInList(Xsd:Kpipi_ch" + str(index) + ")"]
        for index in range(1,24):
            variables = variables + ["nParticlesInList(Xsd:Kpipipi_ch" + str(index) + ")"]

        variables = variables + ["nParticlesInList(Xsu:1270)", "nParticlesInList(Xsu:1430)"]

        variables = variables + ["averageValueInList(Xsu:PrimaryMC, M)", "averageValueInList(Xsd:PrimaryMC, M)", "daughter(0, daughter(0, M))"]
    
        ma.variablesToNtuple(decayString="Upsilon(4S):withoutneutrino",variables=decay_num+systematics+variables,filename=output_file,treename="Xs",path=my_path)
        # ======================================================end=====================================================
    
    if((args.sampletype=="MC_SIGNAL") and (analysistype=="Jpsi")):
        # ==========================================MC data for calculate efficiency====================================
        # Load particles from MCParticle at first
        ma.fillParticleListFromMC('K+:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.fillParticleListFromMC('K*+:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.fillParticleListFromMC('pi+:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.fillParticleListFromMC('e+:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.fillParticleListFromMC('nu_e:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.fillParticleListFromMC('nu_mu:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.fillParticleListFromMC('nu_tau:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.fillParticleListFromMC('gamma:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.fillParticleListFromMC('pi0:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.fillParticleListFromMC("K0:PrimaryMC", cut = "mcPrimary", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.fillParticleListFromMC("K_S0:PrimaryMC", cut = "mcPrimary", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.fillParticleListFromMC("K_L0:PrimaryMC", cut = "mcPrimary", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.fillParticleListFromMC("K*0:PrimaryMC", cut = "mcPrimary", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.fillParticleListFromMC("B+:PrimaryMC", cut = "mcPrimary", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.fillParticleListFromMC("B0:PrimaryMC", cut = "mcPrimary", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.fillParticleListFromMC("Upsilon(4S):PrimaryMC", cut = "mcPrimary", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
    
        ma.reconstructMCDecay('K*+:MC_ch1 =direct=> K0:PrimaryMC pi+:PrimaryMC ?gamma', cut='', dmID = 0, path=my_path)
        ma.reconstructMCDecay('K*+:MC_ch2 =direct=> K+:PrimaryMC pi0:PrimaryMC ?gamma', cut='', dmID = 0, path=my_path)
        ma.reconstructMCDecay('K*0:MC_ch1 =direct=> K+:PrimaryMC pi-:PrimaryMC ?gamma', cut='', dmID = 0, path=my_path)
        ma.reconstructMCDecay('K*0:MC_ch2 =direct=> K0:PrimaryMC pi0:PrimaryMC ?gamma', cut='', dmID = 0, path=my_path)
        ma.reconstructMCDecay('K*0:MC_ch2_short =direct=> K_S0:PrimaryMC pi0:PrimaryMC ?gamma', cut='', dmID = 0, path=my_path)
        ma.reconstructMCDecay('K*0:MC_ch2_long =direct=> K_L0:PrimaryMC pi0:PrimaryMC ?gamma', cut='', dmID = 0, path=my_path)
    
        # Jpsi
        ma.fillParticleListFromMC('mu+:PrimaryMC', cut = 'mcPrimary', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
        ma.reconstructMCDecay('J/psi:mumu =direct=> mu+:PrimaryMC mu-:PrimaryMC', '', path=my_path)
    
        # select Jpsi from B+/0
        ma.cutAndCopyList("J/psi:mumu_Bplus", "J/psi:mumu", "isMCDescendantOfList(B+:PrimaryMC,-1)", path=my_path)
        ma.cutAndCopyList("J/psi:mumu_Bzero", "J/psi:mumu", "isMCDescendantOfList(B0:PrimaryMC,-1)", path=my_path)
    
        # B plus
        ma.reconstructMCDecay('B+:MCch1 -> J/psi:mumu K+:PrimaryMC', cut='', dmID = 1, path=my_path)
        ma.reconstructMCDecay('B+:MCch2 -> J/psi:mumu K+:PrimaryMC pi0:PrimaryMC', cut='', dmID = 2, path=my_path)
        ma.reconstructMCDecay('B+:MCch3 -> J/psi:mumu K0:PrimaryMC pi+:PrimaryMC', cut='', dmID = 3, path=my_path)
        ma.reconstructMCDecay('B+:MCch4 -> J/psi:mumu K+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC', cut='', dmID = 4, path=my_path)
        ma.reconstructMCDecay('B+:MCch5 -> J/psi:mumu K0:PrimaryMC pi+:PrimaryMC pi0:PrimaryMC', cut='', dmID = 5, path=my_path)
        ma.reconstructMCDecay('B+:MCch6 -> J/psi:mumu K+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC pi0:PrimaryMC', cut='', dmID = 6, path=my_path)
        ma.reconstructMCDecay('B+:MCch7 -> J/psi:mumu K0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC', cut='', dmID = 7, path=my_path)
        ma.reconstructMCDecay('B+:MCch8 -> J/psi:mumu K+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC', cut='', dmID = 8, path=my_path)
        ma.reconstructMCDecay('B+:MCch9 -> J/psi:mumu K0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC pi0:PrimaryMC', cut='', dmID = 9, path=my_path)
        # 2pi0
        ma.reconstructMCDecay('B+:MCch10 -> J/psi:mumu K+:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 10, path=my_path)
        ma.reconstructMCDecay('B+:MCch11 -> J/psi:mumu K0:PrimaryMC pi+:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 11, path=my_path)
        ma.reconstructMCDecay('B+:MCch12 -> J/psi:mumu K+:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 12, path=my_path)
        # 3K
        ma.reconstructMCDecay('B+:MCch13 -> J/psi:mumu K+:PrimaryMC K-:PrimaryMC K+:PrimaryMC', cut='', dmID = 13, path=my_path)
        ma.reconstructMCDecay('B+:MCch14 -> J/psi:mumu K+:PrimaryMC K-:PrimaryMC K0:PrimaryMC pi+:PrimaryMC', cut='', dmID = 14, path=my_path)
        ma.reconstructMCDecay('B+:MCch15 -> J/psi:mumu K+:PrimaryMC K-:PrimaryMC K+:PrimaryMC pi0:PrimaryMC', cut='', dmID = 15, path=my_path)
    
        ma.copyLists(outputListName="B+:MCcomb", inputListNames=["B+:MCch1", "B+:MCch2","B+:MCch3","B+:MCch4","B+:MCch5","B+:MCch6","B+:MCch7","B+:MCch8","B+:MCch9","B+:MCch10", "B+:MCch11","B+:MCch12","B+:MCch13","B+:MCch14","B+:MCch15"], path=my_path)
    
        ma.reconstructMCDecay('B+:Kcharge_total =direct=> J/psi:mumu K+:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
    
        ma.reconstructMCDecay('B+:Kstarcharge_ch1_total =direct=> J/psi:mumu K*+:MC_ch1 ?gamma', cut='', dmID = 2, path=my_path)
    
        ma.reconstructMCDecay('B+:Kstarcharge_ch2_total =direct=> J/psi:mumu K*+:MC_ch2 ?gamma', cut='', dmID = 2, path=my_path)
    
        # B zero
        ma.reconstructMCDecay('B0:MCch16 -> J/psi:mumu K0:PrimaryMC', cut='', dmID = 16, path=my_path)
        ma.reconstructMCDecay('B0:MCch17 -> J/psi:mumu K+:PrimaryMC pi-:PrimaryMC', cut='', dmID = 17, path=my_path)
        ma.reconstructMCDecay('B0:MCch18 -> J/psi:mumu K0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 18, path=my_path)
        ma.reconstructMCDecay('B0:MCch19 -> J/psi:mumu K+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC', cut='', dmID = 19, path=my_path)
        ma.reconstructMCDecay('B0:MCch20 -> J/psi:mumu K0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC', cut='', dmID = 20, path=my_path)
        ma.reconstructMCDecay('B0:MCch21 -> J/psi:mumu K+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC', cut='', dmID = 21, path=my_path)
        ma.reconstructMCDecay('B0:MCch22 -> J/psi:mumu K0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC', cut='', dmID = 22, path=my_path)
        ma.reconstructMCDecay('B0:MCch23 -> J/psi:mumu K+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC', cut='', dmID = 23, path=my_path)
        ma.reconstructMCDecay('B0:MCch24 -> J/psi:mumu K0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC', cut='', dmID = 24, path=my_path)
        # 2pi0
        ma.reconstructMCDecay('B0:MCch25 -> J/psi:mumu K0:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 25, path=my_path)
        ma.reconstructMCDecay('B0:MCch26 -> J/psi:mumu K+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 26, path=my_path)
        ma.reconstructMCDecay('B0:MCch27 -> J/psi:mumu K0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 27, path=my_path)
        # 3K
        ma.reconstructMCDecay('B0:MCch28 -> J/psi:mumu K+:PrimaryMC K-:PrimaryMC K0:PrimaryMC', cut='', dmID = 28, path=my_path)
        ma.reconstructMCDecay('B0:MCch29 -> J/psi:mumu K+:PrimaryMC K-:PrimaryMC K+:PrimaryMC pi-:PrimaryMC', cut='', dmID = 29, path=my_path)
        ma.reconstructMCDecay('B0:MCch30 -> J/psi:mumu K+:PrimaryMC K-:PrimaryMC K0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 30, path=my_path)
    
        ma.copyLists(outputListName="B0:MCcomb", inputListNames=["B0:MCch16", "B0:MCch17","B0:MCch18","B0:MCch19","B0:MCch20","B0:MCch21","B0:MCch22","B0:MCch23","B0:MCch24","B0:MCch25", "B0:MCch26","B0:MCch27","B0:MCch28","B0:MCch29","B0:MCch30"], path=my_path)
    
        ma.reconstructMCDecay('B0:Kneutral =direct=> J/psi:mumu K0:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
        ma.reconstructMCDecay('B0:Kneutral_oppose =direct=> J/psi:mumu anti-K0:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
        ma.reconstructMCDecay('B0:Kneutral_short =direct=> J/psi:mumu K_S0:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
        ma.reconstructMCDecay('B0:Kneutral_long =direct=> J/psi:mumu K_L0:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
        ma.copyLists(outputListName="B0:Kneutral_total", inputListNames=["B0:Kneutral", "B0:Kneutral_oppose", "B0:Kneutral_short", "B0:Kneutral_long"], path=my_path)
    
        ma.reconstructMCDecay('B0:Kstarneutral_ch1 =direct=> J/psi:mumu K*0:MC_ch1 ?gamma', cut='', dmID = 2, path=my_path)
        ma.reconstructMCDecay('B0:Kstarneutral_ch1_oppose =direct=> J/psi:mumu anti-K*0:MC_ch1 ?gamma', cut='', dmID = 2, path=my_path)
        ma.copyLists(outputListName="B0:Kstarneutral_ch1_total", inputListNames=["B0:Kstarneutral_ch1", "B0:Kstarneutral_ch1_oppose"], path=my_path)
    
        ma.reconstructMCDecay('B0:Kstarneutral_ch2 =direct=> J/psi:mumu K*0:MC_ch2 ?gamma', cut='', dmID = 2, path=my_path)
        ma.reconstructMCDecay('B0:Kstarneutral_ch2_short =direct=> J/psi:mumu K*0:MC_ch2_short ?gamma', cut='', dmID = 2, path=my_path)
        ma.reconstructMCDecay('B0:Kstarneutral_ch2_long =direct=> J/psi:mumu K*0:MC_ch2_long ?gamma', cut='', dmID = 2, path=my_path)
        ma.copyLists(outputListName="B0:Kstarneutral_ch2_total", inputListNames=["B0:Kstarneutral_ch2", "B0:Kstarneutral_ch2_short", "B0:Kstarneutral_ch2_long"], path=my_path)
    
    
        # Add alias to reuse the code
        va.variables.addAlias("nParticlesInList__boXsu__clMCcomb__bc", "nParticlesInList(B+:MCcomb)")  
        va.variables.addAlias("nParticlesInList__boXsu__clMCch1__bc", "nParticlesInList(B+:MCch1)")  
        va.variables.addAlias("nParticlesInList__boXsu__clMCch2__bc", "nParticlesInList(B+:MCch2)")
        va.variables.addAlias("nParticlesInList__boXsu__clMCch3__bc", "nParticlesInList(B+:MCch3)")
        va.variables.addAlias("nParticlesInList__boXsu__clMCch4__bc", "nParticlesInList(B+:MCch4)")
        va.variables.addAlias("nParticlesInList__boXsu__clMCch5__bc", "nParticlesInList(B+:MCch5)")
        va.variables.addAlias("nParticlesInList__boXsu__clMCch6__bc", "nParticlesInList(B+:MCch6)")
        va.variables.addAlias("nParticlesInList__boXsu__clMCch7__bc", "nParticlesInList(B+:MCch7)")
        va.variables.addAlias("nParticlesInList__boXsu__clMCch8__bc", "nParticlesInList(B+:MCch8)")
        va.variables.addAlias("nParticlesInList__boXsu__clMCch9__bc", "nParticlesInList(B+:MCch9)")
        va.variables.addAlias("nParticlesInList__boXsu__clMCch10__bc", "nParticlesInList(B+:MCch10)")
        va.variables.addAlias("nParticlesInList__boXsu__clMCch11__bc", "nParticlesInList(B+:MCch11)")
        va.variables.addAlias("nParticlesInList__boXsu__clMCch12__bc", "nParticlesInList(B+:MCch12)")
        va.variables.addAlias("nParticlesInList__boXsu__clMCch13__bc", "nParticlesInList(B+:MCch13)")
        va.variables.addAlias("nParticlesInList__boXsu__clMCch14__bc", "nParticlesInList(B+:MCch14)")
        va.variables.addAlias("nParticlesInList__boXsu__clMCch15__bc", "nParticlesInList(B+:MCch15)")
        va.variables.addAlias("nParticlesInList__boXsd__clMCcomb__bc", "nParticlesInList(B0:MCcomb)")
        va.variables.addAlias("nParticlesInList__boXsd__clMCch16__bc", "nParticlesInList(B0:MCch16)")
        va.variables.addAlias("nParticlesInList__boXsd__clMCch17__bc", "nParticlesInList(B0:MCch17)")
        va.variables.addAlias("nParticlesInList__boXsd__clMCch18__bc", "nParticlesInList(B0:MCch18)")
        va.variables.addAlias("nParticlesInList__boXsd__clMCch19__bc", "nParticlesInList(B0:MCch19)")
        va.variables.addAlias("nParticlesInList__boXsd__clMCch20__bc", "nParticlesInList(B0:MCch20)")
        va.variables.addAlias("nParticlesInList__boXsd__clMCch21__bc", "nParticlesInList(B0:MCch21)")
        va.variables.addAlias("nParticlesInList__boXsd__clMCch22__bc", "nParticlesInList(B0:MCch22)")
        va.variables.addAlias("nParticlesInList__boXsd__clMCch23__bc", "nParticlesInList(B0:MCch23)")
        va.variables.addAlias("nParticlesInList__boXsd__clMCch24__bc", "nParticlesInList(B0:MCch24)")
        va.variables.addAlias("nParticlesInList__boXsd__clMCch25__bc", "nParticlesInList(B0:MCch25)")
        va.variables.addAlias("nParticlesInList__boXsd__clMCch26__bc", "nParticlesInList(B0:MCch26)")
        va.variables.addAlias("nParticlesInList__boXsd__clMCch27__bc", "nParticlesInList(B0:MCch27)")
        va.variables.addAlias("nParticlesInList__boXsd__clMCch28__bc", "nParticlesInList(B0:MCch28)")
        va.variables.addAlias("nParticlesInList__boXsd__clMCch29__bc", "nParticlesInList(B0:MCch29)")
        va.variables.addAlias("nParticlesInList__boXsd__clMCch30__bc", "nParticlesInList(B0:MCch30)")
        va.variables.addAlias("nParticlesInList__bonu_e__clMC_signal__bc", "constant(0)")
        va.variables.addAlias("nParticlesInList__boB__pl__clMC_signal_total_e__bc", "nParticlesInList(J/psi:mumu_Bplus)")
        va.variables.addAlias("nParticlesInList__boB0__clMC_signal_total_e__bc", "nParticlesInList(J/psi:mumu_Bzero)")
        va.variables.addAlias("invMassInLists__bonu_e__clMC_signal__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boB__pl__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boB0__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boB__pl__clMC_signal_total_e__cm__spM__bc", "constant(0)")
        va.variables.addAlias("averageValueInList__boB0__clMC_signal_total_e__cm__spM__bc", "constant(0)")
    
        decay_num = []
        decay_num.append("nParticlesInList(B+:Kcharge_total)")
        decay_num.append("nParticlesInList(B+:Kstarcharge_ch1_total)")
        decay_num.append("nParticlesInList(B+:Kstarcharge_ch2_total)")
        decay_num.append("nParticlesInList__boXsu__clMCcomb__bc")
        for i in range(1,16):
            decay_num.append("nParticlesInList__boXsu__clMCch" + str(i) + "__bc")
        decay_num.append("nParticlesInList__boXsd__clMCcomb__bc")
        decay_num.append("nParticlesInList(B0:Kneutral_total)")
        decay_num.append("nParticlesInList(B0:Kstarneutral_ch1_total)")
        decay_num.append("nParticlesInList(B0:Kstarneutral_ch2_total)")
        for i in range(16,31):
            decay_num.append("nParticlesInList__boXsd__clMCch" + str(i) + "__bc")
    
        systematics = ["nParticlesInList__bonu_e__clMC_signal__bc", "nParticlesInList__boB__pl__clMC_signal_total_e__bc", "nParticlesInList__boB0__clMC_signal_total_e__bc", "invMassInLists__bonu_e__clMC_signal__bc", "averageValueInList__boB__pl__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", "averageValueInList__boB0__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", "averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", "averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", "averageValueInList__boB__pl__clMC_signal_total_e__cm__spM__bc", "averageValueInList__boB0__clMC_signal_total_e__cm__spM__bc", "nParticlesInList(B+:PrimaryMC)", "nParticlesInList(B0:PrimaryMC)"]
    
        ma.variablesToNtuple(decayString="Upsilon(4S):withoutneutrino",variables=decay_num+systematics,filename=output_file,treename="Xs",path=my_path)
    
    # progress
    basf2.process(my_path)

# Print call statistics
print(basf2.statistics)
