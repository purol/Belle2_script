#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# usage: basf2 MakeNtuple_multi.py "./20210402/evt-0.mdst"
# last: 2021-09-03-0

import os
import sys

import basf2
import modularAnalysis as ma
import variables as va
import variables.collections as vc
import variables.utils as vu

from glob import glob

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

if not fName.endswith(".root"): sys.exit(1)
else:
    basename = os.path.basename(fName)
    name = os.path.splitext(basename)[0]
print("=================="+fName+" is conducted=================")

my_path = basf2.create_path()
inputfile=fName

ma.inputMdst(environmentType='default',filename=inputfile,path=my_path)

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

ma.fillParticleListFromMC("eta:PrimaryMC", cut = "mcPrimary", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)

# for excited K
ma.fillParticleListFromMC("rho+:PrimaryMC", cut = "mcPrimary", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
ma.fillParticleListFromMC("rho0:PrimaryMC", cut = "mcPrimary", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)
ma.fillParticleListFromMC("omega:PrimaryMC", cut = "mcPrimary", addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)

ma.reconstructMCDecay('K*+:PrimaryMC_ch1 =direct=> K0:PrimaryMC pi+:PrimaryMC ?gamma', cut='', dmID = 0, path=my_path)
ma.reconstructMCDecay('K*+:PrimaryMC_ch2 =direct=> K+:PrimaryMC pi0:PrimaryMC ?gamma', cut='', dmID = 0, path=my_path)
ma.reconstructMCDecay('K*0:PrimaryMC_ch1 =direct=> K+:PrimaryMC pi-:PrimaryMC ?gamma', cut='', dmID = 0, path=my_path)
ma.reconstructMCDecay('K*0:PrimaryMC_ch2 =direct=> K0:PrimaryMC pi0:PrimaryMC ?gamma', cut='', dmID = 0, path=my_path)

# B plus
ma.reconstructMCDecay('Xsu:ch1 -> K+:PrimaryMC', cut='', dmID = 1, path=my_path)
ma.reconstructMCDecay('Xsu:ch2 -> K+:PrimaryMC pi0:PrimaryMC', cut='', dmID = 2, path=my_path)
ma.reconstructMCDecay('Xsu:ch3 -> K0:PrimaryMC pi+:PrimaryMC', cut='', dmID = 3, path=my_path)
ma.reconstructMCDecay('Xsu:ch4 -> K+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC', cut='', dmID = 4, path=my_path)
ma.reconstructMCDecay('Xsu:ch5 -> K0:PrimaryMC pi+:PrimaryMC pi0:PrimaryMC', cut='', dmID = 5, path=my_path)
ma.reconstructMCDecay('Xsu:ch6 -> K+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC pi0:PrimaryMC', cut='', dmID = 6, path=my_path)
ma.reconstructMCDecay('Xsu:ch7 -> K0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC', cut='', dmID = 7, path=my_path)
ma.reconstructMCDecay('Xsu:ch8 -> K+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC', cut='', dmID = 8, path=my_path)
ma.reconstructMCDecay('Xsu:ch9 -> K0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC pi0:PrimaryMC', cut='', dmID = 9, path=my_path)
# 2pi0
ma.reconstructMCDecay('Xsu:ch10 -> K+:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 10, path=my_path)
ma.reconstructMCDecay('Xsu:ch11 -> K0:PrimaryMC pi+:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 11, path=my_path)
ma.reconstructMCDecay('Xsu:ch12 -> K+:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 12, path=my_path)
# 3K
ma.reconstructMCDecay('Xsu:ch13 -> K+:PrimaryMC K-:PrimaryMC K+:PrimaryMC', cut='', dmID = 13, path=my_path)
ma.reconstructMCDecay('Xsu:ch14 -> K+:PrimaryMC K-:PrimaryMC K0:PrimaryMC pi+:PrimaryMC', cut='', dmID = 14, path=my_path)
ma.reconstructMCDecay('Xsu:ch15 -> K+:PrimaryMC K-:PrimaryMC K+:PrimaryMC pi0:PrimaryMC', cut='', dmID = 15, path=my_path)

# eta
ma.reconstructMCDecay('Xsu:eta_ch1 -> K+:PrimaryMC eta:PrimaryMC', cut='', path=my_path)
ma.reconstructMCDecay('Xsu:eta_ch2 -> K0:PrimaryMC eta:PrimaryMC pi+:PrimaryMC', cut='', path=my_path)
ma.reconstructMCDecay('Xsu:eta_ch3 -> K+:PrimaryMC eta:PrimaryMC pi0:PrimaryMC', cut='', path=my_path)
ma.reconstructMCDecay('Xsu:eta_ch4 -> K+:PrimaryMC eta:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC', cut='', path=my_path)
ma.reconstructMCDecay('Xsu:eta_ch5 -> K0:PrimaryMC eta:PrimaryMC pi+:PrimaryMC pi0:PrimaryMC', cut='', path=my_path)
Xsu_etas = ["Xsu:eta_ch1", "Xsu:eta_ch2", "Xsu:eta_ch3", "Xsu:eta_ch4", "Xsu:eta_ch5"]
ma.copyLists(outputListName="Xsu:comb", inputListNames=["Xsu:ch1", "Xsu:ch2","Xsu:ch3","Xsu:ch4","Xsu:ch5","Xsu:ch6","Xsu:ch7","Xsu:ch8","Xsu:ch9","Xsu:ch10", "Xsu:ch11","Xsu:ch12","Xsu:ch13","Xsu:ch14","Xsu:ch15"]+Xsu_etas, path=my_path)

ma.reconstructMCDecay('B+:Kcharge_e =direct=> K+:PrimaryMC nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
ma.reconstructMCDecay('B+:Kcharge_mu =direct=> K+:PrimaryMC nu_mu:PrimaryMC anti-nu_mu:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
ma.reconstructMCDecay('B+:Kcharge_tau =direct=> K+:PrimaryMC nu_tau:PrimaryMC anti-nu_tau:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
ma.copyLists(outputListName="B+:Kcharge_total", inputListNames=["B+:Kcharge_e", "B+:Kcharge_mu", "B+:Kcharge_tau"], path=my_path)

ma.reconstructMCDecay('B+:Kstarcharge_ch1_e =direct=> K*+:PrimaryMC_ch1 nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.reconstructMCDecay('B+:Kstarcharge_ch1_mu =direct=> K*+:PrimaryMC_ch1 nu_mu:PrimaryMC anti-nu_mu:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.reconstructMCDecay('B+:Kstarcharge_ch1_tau =direct=> K*+:PrimaryMC_ch1 nu_tau:PrimaryMC anti-nu_tau:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.copyLists(outputListName="B+:Kstarcharge_ch1_total", inputListNames=["B+:Kstarcharge_ch1_e", "B+:Kstarcharge_ch1_mu", "B+:Kstarcharge_ch1_tau"], path=my_path)

ma.reconstructMCDecay('B+:Kstarcharge_ch2_e =direct=> K*+:PrimaryMC_ch2 nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.reconstructMCDecay('B+:Kstarcharge_ch2_mu =direct=> K*+:PrimaryMC_ch2 nu_mu:PrimaryMC anti-nu_mu:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.reconstructMCDecay('B+:Kstarcharge_ch2_tau =direct=> K*+:PrimaryMC_ch2 nu_tau:PrimaryMC anti-nu_tau:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.copyLists(outputListName="B+:Kstarcharge_ch2_total", inputListNames=["B+:Kstarcharge_ch2_e", "B+:Kstarcharge_ch2_mu", "B+:Kstarcharge_ch2_tau"], path=my_path)

# B zero
ma.reconstructMCDecay('Xsd:ch16 -> K0:PrimaryMC', cut='', dmID = 16, path=my_path)
ma.reconstructMCDecay('Xsd:ch17 -> K+:PrimaryMC pi-:PrimaryMC', cut='', dmID = 17, path=my_path)
ma.reconstructMCDecay('Xsd:ch18 -> K0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 18, path=my_path)
ma.reconstructMCDecay('Xsd:ch19 -> K+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC', cut='', dmID = 19, path=my_path)
ma.reconstructMCDecay('Xsd:ch20 -> K0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC', cut='', dmID = 20, path=my_path)
ma.reconstructMCDecay('Xsd:ch21 -> K+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC', cut='', dmID = 21, path=my_path)
ma.reconstructMCDecay('Xsd:ch22 -> K0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC', cut='', dmID = 22, path=my_path)
ma.reconstructMCDecay('Xsd:ch23 -> K+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC', cut='', dmID = 23, path=my_path)
ma.reconstructMCDecay('Xsd:ch24 -> K0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC', cut='', dmID = 24, path=my_path)
# 2pi0
ma.reconstructMCDecay('Xsd:ch25 -> K0:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 25, path=my_path)
ma.reconstructMCDecay('Xsd:ch26 -> K+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 26, path=my_path)
ma.reconstructMCDecay('Xsd:ch27 -> K0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 27, path=my_path)
# 3K
ma.reconstructMCDecay('Xsd:ch28 -> K+:PrimaryMC K-:PrimaryMC K0:PrimaryMC', cut='', dmID = 28, path=my_path)
ma.reconstructMCDecay('Xsd:ch29 -> K+:PrimaryMC K-:PrimaryMC K+:PrimaryMC pi-:PrimaryMC', cut='', dmID = 29, path=my_path)
ma.reconstructMCDecay('Xsd:ch30 -> K+:PrimaryMC K-:PrimaryMC K0:PrimaryMC pi0:PrimaryMC', cut='', dmID = 30, path=my_path)

# eta
ma.reconstructMCDecay('Xsd:eta_ch1 -> K0:PrimaryMC eta:PrimaryMC', cut='', path=my_path) 
ma.reconstructMCDecay('Xsd:eta_ch2 -> K+:PrimaryMC eta:PrimaryMC pi-:PrimaryMC', cut='', path=my_path)
ma.reconstructMCDecay('Xsd:eta_ch3 -> K0:PrimaryMC eta:PrimaryMC pi0:PrimaryMC', cut='', path=my_path)
ma.reconstructMCDecay('Xsd:eta_ch4 -> K0:PrimaryMC eta:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC', cut='', path=my_path)
ma.reconstructMCDecay('Xsd:eta_ch5 -> K+:PrimaryMC eta:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC', cut='', path=my_path)
Xsd_etas = ["Xsd:eta_ch1", "Xsd:eta_ch2", "Xsd:eta_ch3", "Xsd:eta_ch4", "Xsd:eta_ch5"]
ma.copyLists(outputListName="Xsd:comb", inputListNames=["Xsd:ch16", "Xsd:ch17","Xsd:ch18","Xsd:ch19","Xsd:ch20","Xsd:ch21","Xsd:ch22","Xsd:ch23","Xsd:ch24","Xsd:ch25", "Xsd:ch26","Xsd:ch27","Xsd:ch28","Xsd:ch29","Xsd:ch30"]+Xsd_etas, path=my_path)

ma.reconstructMCDecay('B0:Kneutral_e =direct=> K0:PrimaryMC nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
ma.reconstructMCDecay('B0:Kneutral_mu =direct=> K0:PrimaryMC nu_mu:PrimaryMC anti-nu_mu:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
ma.reconstructMCDecay('B0:Kneutral_tau =direct=> K0:PrimaryMC nu_tau:PrimaryMC anti-nu_tau:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
ma.reconstructMCDecay('B0:Kneutral_e_oppose =direct=> anti-K0:PrimaryMC nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
ma.reconstructMCDecay('B0:Kneutral_mu_oppose =direct=> anti-K0:PrimaryMC nu_mu:PrimaryMC anti-nu_mu:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
ma.reconstructMCDecay('B0:Kneutral_tau_oppose =direct=> anti-K0:PrimaryMC nu_tau:PrimaryMC anti-nu_tau:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
ma.copyLists(outputListName="B0:Kneutral_total", inputListNames=["B0:Kneutral_e", "B0:Kneutral_mu", "B0:Kneutral_tau", "B0:Kneutral_e_oppose", "B0:Kneutral_mu_oppose", "B0:Kneutral_tau_oppose"], path=my_path)

ma.reconstructMCDecay('B0:Kstarneutral_ch1_e =direct=> K*0:PrimaryMC_ch1 nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.reconstructMCDecay('B0:Kstarneutral_ch1_mu =direct=> K*0:PrimaryMC_ch1 nu_mu:PrimaryMC anti-nu_mu:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.reconstructMCDecay('B0:Kstarneutral_ch1_tau =direct=> K*0:PrimaryMC_ch1 nu_tau:PrimaryMC anti-nu_tau:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.copyLists(outputListName="B0:Kstarneutral_ch1_total", inputListNames=["B0:Kstarneutral_ch1_e", "B0:Kstarneutral_ch1_mu", "B0:Kstarneutral_ch1_tau"], path=my_path)

ma.reconstructMCDecay('B0:Kstarneutral_ch2_e =direct=> K*0:PrimaryMC_ch2 nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.reconstructMCDecay('B0:Kstarneutral_ch2_mu =direct=> K*0:PrimaryMC_ch2 nu_mu:PrimaryMC anti-nu_mu:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.reconstructMCDecay('B0:Kstarneutral_ch2_tau =direct=> K*0:PrimaryMC_ch2 nu_tau:PrimaryMC anti-nu_tau:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.copyLists(outputListName="B0:Kstarneutral_ch2_total", inputListNames=["B0:Kstarneutral_ch2_e", "B0:Kstarneutral_ch2_mu", "B0:Kstarneutral_ch2_tau"], path=my_path)

# calculate q2 and costheta for systematic uncertainty
# B+:Kcharge_e
ma.reconstructMCDecay('B+:MC_Kstarcharge_e =direct=> K*+:PrimaryMC nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', path=my_path)
ma.reconstructMCDecay('B+:MC_Xsu_e =direct=> Xsu:PrimaryMC nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', path=my_path)
ma.copyLists(outputListName="B+:MC_signal_total_e", inputListNames=["B+:Kcharge_e", "B+:MC_Kstarcharge_e", "B+:MC_Xsu_e"], path=my_path)

ma.copyLists(outputListName="B0:MC_Kneutral_total_e", inputListNames=["B0:Kneutral_e", "B0:Kneutral_e_oppose"], path=my_path)
ma.reconstructMCDecay('B0:MC_Kstarneutral_e =direct=> K*0:PrimaryMC nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', path=my_path)
ma.reconstructMCDecay('B0:MC_Xsd_e =direct=> Xsd:PrimaryMC nu_e:PrimaryMC anti-nu_e:PrimaryMC ?gamma', cut='', path=my_path)
ma.copyLists(outputListName="B0:MC_signal_total_e", inputListNames=["B0:MC_Kneutral_total_e", "B0:MC_Kstarneutral_e", "B0:MC_Xsd_e"], path=my_path)

ma.fillParticleListFromMC('nu_e:MC_signal', cut = 'mcPrimary and [isDescendantOfList(B+:MC_signal_total_e,1) or isDescendantOfList(B0:MC_signal_total_e,1)]', addDaughters=True, skipNonPrimaryDaughters=True, path=my_path)

# for excited K
ma.reconstructMCDecay('Xsu:excited_ch1 =direct=> rho+:PrimaryMC K0:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
ma.reconstructMCDecay('Xsu:excited_ch2 =direct=> rho0:PrimaryMC K+:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.reconstructMCDecay('Xsu:excited_ch3 =direct=> K*0:PrimaryMC pi+:PrimaryMC ?gamma', cut='', dmID = 3, path=my_path)
ma.reconstructMCDecay('Xsu:excited_ch4 =direct=> K*+:PrimaryMC pi0:PrimaryMC ?gamma', cut='', dmID = 4, path=my_path)
ma.reconstructMCDecay('Xsu:excited_ch5 =direct=> omega:PrimaryMC K+:PrimaryMC ?gamma', cut='', dmID = 5, path=my_path)
ma.reconstructMCDecay('Xsu:excited_ch6 =direct=> K+:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC ?gamma', cut='', dmID = 6, path=my_path)
ma.reconstructMCDecay('Xsu:excited_ch7 =direct=> K0:PrimaryMC pi+:PrimaryMC pi0:PrimaryMC ?gamma', cut='', dmID = 7, path=my_path)
ma.reconstructMCDecay('Xsu:excited_ch8 =direct=> K+:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC ?gamma', cut='', dmID = 8, path=my_path)
ma.reconstructMCDecay('Xsu:excited_ch9 =direct=> K0:PrimaryMC pi+:PrimaryMC ?gamma', cut='', dmID = 9, path=my_path)
ma.reconstructMCDecay('Xsu:excited_ch10 =direct=> K+:PrimaryMC pi0:PrimaryMC ?gamma', cut='', dmID = 10, path=my_path)

ma.reconstructMCDecay('Xsd:excited_ch1 =direct=> rho-:PrimaryMC K+:PrimaryMC ?gamma', cut='', dmID = 1, path=my_path)
ma.reconstructMCDecay('Xsd:excited_ch2 =direct=> rho0:PrimaryMC K0:PrimaryMC ?gamma', cut='', dmID = 2, path=my_path)
ma.reconstructMCDecay('Xsd:excited_ch3 =direct=> K*+:PrimaryMC pi-:PrimaryMC ?gamma', cut='', dmID = 3, path=my_path)
ma.reconstructMCDecay('Xsd:excited_ch4 =direct=> K*0:PrimaryMC pi0:PrimaryMC ?gamma', cut='', dmID = 4, path=my_path)
ma.reconstructMCDecay('Xsd:excited_ch5 =direct=> omega:PrimaryMC K0:PrimaryMC ?gamma', cut='', dmID = 5, path=my_path)
ma.reconstructMCDecay('Xsd:excited_ch6 =direct=> K0:PrimaryMC pi+:PrimaryMC pi-:PrimaryMC ?gamma', cut='', dmID = 6, path=my_path)
ma.reconstructMCDecay('Xsd:excited_ch7 =direct=> K+:PrimaryMC pi-:PrimaryMC pi0:PrimaryMC ?gamma', cut='', dmID = 7, path=my_path)
ma.reconstructMCDecay('Xsd:excited_ch8 =direct=> K0:PrimaryMC pi0:PrimaryMC pi0:PrimaryMC ?gamma', cut='', dmID = 8, path=my_path)
ma.reconstructMCDecay('Xsd:excited_ch9 =direct=> K+:PrimaryMC pi-:PrimaryMC ?gamma', cut='', dmID = 9, path=my_path)
ma.reconstructMCDecay('Xsd:excited_ch10 =direct=> K0:PrimaryMC pi0:PrimaryMC ?gamma', cut='', dmID = 10, path=my_path)

# Upsilon

#ma.reconstructMCDecay('Xsu:ch2 -> K+:PrimaryMC pi0:PrimaryMC', cut='', dmID = 2, path=my_path)
#ma.reconstructMCDecay('Xsu:ch3 -> K0:PrimaryMC pi+:PrimaryMC', cut='', dmID = 3, path=my_path)
#ma.reconstructMCDecay('Xsu:ch4 -> K+:PrimaryMC pi-:PrimaryMC pi+:PrimaryMC', cut='', dmID = 4, path=my_path)
#ma.reconstructMCDecay('Xsu:ch5 -> K0:PrimaryMC pi+:PrimaryMC pi0:PrimaryMC', cut='', dmID = 5, path=my_path)

#ma.copyLists(outputListName="Xsu:total", inputListNames=["Xsu:ch1", "Xsu:ch2", "Xsu:ch3", "Xsu:ch4", "Xsu:ch5"], path=my_path)

#my_path.add_module("ParticlePrinter",listName="Xsu:total", variables=["M"])

#ma.applyEventCuts('[nParticlesInList(B+:ch_total) == 1 and nParticlesInList(B0:ch_total) == 0] or [nParticlesInList(B+:ch_total) == 0 and nParticlesInList(B0:ch_total) == 1]',path=my_path)

Tree_shape = ["daughter(0, nDaughters)", "daughter(0, PDG)", "daughter(0, daughter(0,PDG))", "daughter(0, daughter(1,PDG))", "daughter(0, daughter(2,PDG))", "daughter(0, daughter(0,M))", "daughter(1, nDaughters)", "daughter(1, PDG)", "daughter(1, daughter(0,PDG))", "daughter(1, daughter(1,PDG))", "daughter(1, daughter(2,PDG))", "daughter(1, daughter(0,M))"]
decay_num = ["nParticlesInList(B+:Kcharge_total)", "nParticlesInList(B+:Kstarcharge_ch1_total)", "nParticlesInList(B+:Kstarcharge_ch2_total)", "nParticlesInList(Xsu:comb)", "nParticlesInList(Xsu:ch1)", "nParticlesInList(Xsu:ch2)", "nParticlesInList(Xsu:ch3)", "nParticlesInList(Xsu:ch4)", "nParticlesInList(Xsu:ch5)", "nParticlesInList(Xsu:ch6)", "nParticlesInList(Xsu:ch7)", "nParticlesInList(Xsu:ch8)", "nParticlesInList(Xsu:ch9)", "nParticlesInList(Xsu:ch10)", "nParticlesInList(Xsu:ch11)", "nParticlesInList(Xsu:ch12)", "nParticlesInList(Xsu:ch13)", "nParticlesInList(Xsu:ch14)", "nParticlesInList(Xsu:ch15)", "nParticlesInList(B0:Kneutral_total)", "nParticlesInList(B0:Kstarneutral_ch1_total)", "nParticlesInList(B0:Kstarneutral_ch2_total)", "nParticlesInList(Xsd:comb)", "nParticlesInList(Xsd:ch16)", "nParticlesInList(Xsd:ch17)", "nParticlesInList(Xsd:ch18)", "nParticlesInList(Xsd:ch19)", "nParticlesInList(Xsd:ch20)", "nParticlesInList(Xsd:ch21)", "nParticlesInList(Xsd:ch22)", "nParticlesInList(Xsd:ch23)", "nParticlesInList(Xsd:ch24)", "nParticlesInList(Xsd:ch25)", "nParticlesInList(Xsd:ch26)", "nParticlesInList(Xsd:ch27)", "nParticlesInList(Xsd:ch28)", "nParticlesInList(Xsd:ch29)", "nParticlesInList(Xsd:ch30)"]
q2 = ["invMassInLists(nu_e:MC_signal)"]
systematics = ["averageValueInList(B+:MC_signal_total_e, daughter(0, M))", "averageValueInList(B0:MC_signal_total_e, daughter(0, M))", "averageValueInList(B+:MC_signal_total_e, M)", "averageValueInList(B0:MC_signal_total_e, M)"]

output_file = destination + "/" + name+"_Decay_Info.root"
ma.variablesToNtuple(decayString="Upsilon(4S):PrimaryMC",variables=Tree_shape+decay_num+q2+systematics,filename=output_file,treename="Xs",path=my_path)

#ma.applyEventCuts("nParticlesInList(Xsu:1410) > 0", path=my_path)
#ma.printMCParticles(onlyPrimaries=True, path=my_path)

#my_path.add_module("ParticlePrinter",listName="Upsilon(4S):PrimaryMC", variables=[])
#my_path.add_module("ParticlePrinter",listName="Xsu:comb", variables=["M"])
#my_path.add_module("ParticlePrinter",listName="B+:ch_total", variables=["M"])
#my_path.add_module("ParticlePrinter",listName="B+:PrimaryMC", variables=["M"])
#my_path.add_module("ParticlePrinter",listName="Xsu:ch5", variables=["M"])

#ma.printMCParticles(onlyPrimaries=True, path=my_path)
#my_path.add_module("ParticlePrinter",listName="Xsu:ch3", variables=["M"])
#my_path.add_module("ParticlePrinter",listName="K0:PrimaryMC", variables=["M"])

# progress
basf2.process(my_path)

# Print call statistics
print(basf2.statistics)
