#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# usage: basf2 MakeNtuple_multi.py "./20210402/evt-0.mdst"
# last: 2021-10-29-00

import os
import sys

import basf2
import modularAnalysis as ma
import variables as va
import variables.collections as vc
import variables.utils as vu
import vertex

from glob import glob

basf2.conditions.prepend_globaltag(ma.getAnalysisGlobaltag())

# set random seed
basf2.set_random_seed(42)

# Read uDST
my_path = basf2.create_path()

inputfile="/group/belle2/dataprod/MC/SkimTraining/mixed_BGx1.mdst_000001_prod00009434_task10020000001.root"
ma.inputMdst(environmentType='default',filename=inputfile,path=my_path)

# all particle list
ma.fillParticleList(decayString="K+:all", cut="",path=my_path)
ma.fillParticleList(decayString="pi+:all", cut="",path=my_path)
ma.fillParticleList(decayString="gamma:all", cut="",path=my_path)

# primary particle list
ma.fillParticleList(decayString="K+:mychargedKaon", cut="kaonID>0.6 and nCDCHits > 20 and dr < 0.5 and abs(dz) < 2",path=my_path)
ma.fillParticleList(decayString="pi+:mychargedPion", cut="pionID>0.4 and nCDCHits > 20 and dr < 0.5 and abs(dz) < 2",path=my_path)

ma.reconstructDecay(decayString="K_S0:myKaonshort -> pi+:all pi-:all", cut="M > 0.3 and M < 0.7",path=my_path)
vertex.kFit('K_S0:myKaonshort', 0, path=my_path)
ma.applyCuts("K_S0:myKaonshort", "abs(dM) < 0.01 and significanceOfDistance > 50", path=my_path)

ma.fillParticleList(decayString="gamma:mygamma", cut="inCDCAcceptance and clusterNHits > 1.5 and [[clusterReg==1 and E>0.08] or [clusterReg==2 and E>0.03] or [clusterReg==3 and E>0.06]]",path=my_path)

ma.reconstructDecay(decayString="pi0:myneutralPion -> gamma:mygamma gamma:mygamma", cut="M > 0.12 and M < 0.145 and p > 0.4 and abs(daughterAngle(0,1))<1.4 and abs(daughterDiffOfPhi(0,1))<1.5", path=my_path)

# non-primary particle
#ma.reconstructDecay(decayString="K*0:myneutralKaonstar -> K+:mychargedKaon pi-:mychargedPion", cut="",path=my_path)
#vertex.kFit("K*0:myneutralKaonstar", 0, path=my_path)

#ma.reconstructDecay(decayString="K*+:mychargedKaonstar -> K_S0:myKaonshort pi+:mychargedPion", cut="", dmID=1, path=my_path)


# Xsu
ma.reconstructDecay(decayString="@Xsu:ch0 -> K+:mychargedKaon", cut="M < 2.2", dmID=0, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch1 -> K+:mychargedKaon pi0:myneutralPion", cut="M < 2.2", dmID=1, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch2 -> K_S0:myKaonshort pi+:mychargedPion", cut="M < 2.2", dmID=2, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch3 -> K+:mychargedKaon pi+:mychargedPion pi-:mychargedPion", cut="M < 2.2", dmID=3, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch4 -> K_S0:myKaonshort pi+:mychargedPion pi0:myneutralPion", cut="M < 2.2", dmID=4, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch5 -> K+:mychargedKaon pi-:mychargedPion pi+:mychargedPion pi0:myneutralPion", cut="M < 2.2", dmID=5, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch6 -> K_S0:myKaonshort pi+:mychargedPion pi-:mychargedPion pi+:mychargedPion", cut="M < 2.2", dmID=6, path=my_path)
#ma.reconstructDecay(decayString="@Xsu:ch7 -> K+:mychargedKaon pi-:mychargedPion pi+:mychargedPion pi-:mychargedPion pi+:mychargedPion", cut="M < 2.2", dmID=7, path=my_path)
#ma.reconstructDecay(decayString="@Xsu:ch8 -> K_S0:myKaonshort pi+:mychargedPion pi-:mychargedPion pi+:mychargedPion pi0:myneutralPion", cut="M < 2.2", dmID=8, path=my_path)
# 2pi0
#ma.reconstructDecay(decayString="@Xsu:ch9 -> K+:mychargedKaon pi0:myneutralPion  pi0:myneutralPion", cut="M < 2.2", dmID=9, path=my_path)
#ma.reconstructDecay(decayString="@Xsu:ch10 -> K_S0:myKaonshort pi+:mychargedPion pi0:myneutralPion pi0:myneutralPion", cut="M < 2.2", dmID=10, path=my_path)
#ma.reconstructDecay(decayString="@Xsu:ch11 -> K+:mychargedKaon pi+:mychargedPion pi-:mychargedPion pi0:myneutralPion pi0:myneutralPion", cut="M < 2.2", dmID=11, path=my_path)
# 3K
ma.reconstructDecay(decayString="@Xsu:ch12 -> K+:mychargedKaon K-:mychargedKaon K+:mychargedKaon", cut="M < 2.2", dmID=12, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch13 -> K+:mychargedKaon K-:mychargedKaon K_S0:myKaonshort pi+:mychargedPion", cut="M < 2.2", dmID=13, path=my_path)
ma.reconstructDecay(decayString="@Xsu:ch14 -> K+:mychargedKaon K-:mychargedKaon K+:mychargedKaon pi0:myneutralPion", cut="M < 2.2", dmID=14, path=my_path)
ma.copyLists(outputListName="Xsu:comb", inputListNames=["Xsu:ch0", "Xsu:ch1", "Xsu:ch2", "Xsu:ch3", "Xsu:ch4", "Xsu:ch5", "Xsu:ch6", "Xsu:ch12", "Xsu:ch13", "Xsu:ch14"], path=my_path)


# Xsd
ma.reconstructDecay(decayString="@Xsd:ch0 -> K_S0:myKaonshort", cut="M < 2.2", dmID=0, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch1 -> K+:mychargedKaon pi-:mychargedPion", cut="M < 2.2", dmID=1, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch2 -> K_S0:myKaonshort pi0:myneutralPion", cut="M < 2.2", dmID=2, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch3 -> K+:mychargedKaon pi-:mychargedPion pi0:myneutralPion", cut="M < 2.2", dmID=3, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch4 -> K_S0:myKaonshort pi+:mychargedPion pi-:mychargedPion", cut="M < 2.2", dmID=4, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch5 -> K+:mychargedKaon pi-:mychargedPion pi+:mychargedPion pi-:mychargedPion", cut="M < 2.2", dmID=5, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch6 -> K_S0:myKaonshort pi+:mychargedPion pi-:mychargedPion pi0:myneutralPion", cut="M < 2.2", dmID=6, path=my_path)
#ma.reconstructDecay(decayString="@Xsd:ch7 -> K+:mychargedKaon pi-:mychargedPion pi+:mychargedPion pi-:mychargedPion pi0:myneutralPion", cut="M < 2.2", dmID=7, path=my_path)
#ma.reconstructDecay(decayString="@Xsd:ch8 -> K_S0:myKaonshort pi+:mychargedPion pi-:mychargedPion pi+:mychargedPion pi-:mychargedPion", cut="M < 2.2", dmID=8, path=my_path)
# 2pi0
#ma.reconstructDecay(decayString="@Xsd:ch9 -> K_S0:myKaonshort pi0:myneutralPion pi0:myneutralPion", cut="M < 2.2", dmID=9, path=my_path)
#ma.reconstructDecay(decayString="@Xsd:ch10 -> K+:mychargedKaon pi-:mychargedPion pi0:myneutralPion pi0:myneutralPion", cut="M < 2.2", dmID=10, path=my_path)
#ma.reconstructDecay(decayString="@Xsd:ch11 -> K_S0:myKaonshort pi+:mychargedPion pi-:mychargedPion pi0:myneutralPion pi0:myneutralPion", cut="M < 2.2", dmID=11, path=my_path)
# 3K
ma.reconstructDecay(decayString="@Xsd:ch12 -> K+:mychargedKaon K-:mychargedKaon K_S0:myKaonshort", cut="M < 2.2", dmID=12, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch13 -> K+:mychargedKaon K-:mychargedKaon K+:mychargedKaon pi-:mychargedPion", cut="M < 2.2", dmID=13, path=my_path)
ma.reconstructDecay(decayString="@Xsd:ch14 -> K+:mychargedKaon K-:mychargedKaon K_S0:myKaonshort pi0:myneutralPion", cut="M < 2.2", dmID=14, path=my_path)
ma.copyLists(outputListName="Xsd:comb", inputListNames=["Xsd:ch0", "Xsd:ch1", "Xsd:ch2", "Xsd:ch3", "Xsd:ch4", "Xsd:ch5", "Xsd:ch6", "Xsd:ch12", "Xsd:ch13", "Xsd:ch14"], path=my_path)


# reconstruct J/psi and pick appropriate one
ma.fillParticleList(decayString="mu+:mychargedMuon", cut="muonID>0.9 and p>0.7 and dr < 0.5 and abs(dz) < 2",path=my_path)
ma.reconstructDecay(decayString="J/psi:myJpsi_beforedMcut -> mu+:mychargedMuon mu-:mychargedMuon", cut="",path=my_path)
ma.cutAndCopyList('J/psi:myJpsi', 'J/psi:myJpsi_beforedMcut','abs(dM)<0.05', path=my_path)

ma.reconstructDecay("B+:sig_Jpsi -> Xsu:comb J/psi:myJpsi", cut="", dmID=0, path = my_path)
ma.reconstructDecay("B0:sig_Jpsi -> Xsd:comb J/psi:myJpsi", cut="", dmID=0, path = my_path)
#ma.applyCuts('B+:sig_Jpsi', 'Mbc > 5.25 and abs(deltaE)<0.10', path=my_path)
#ma.applyCuts('B0:sig_Jpsi', 'Mbc > 5.25 and abs(deltaE)<0.10', path=my_path)

ma.reconstructDecay("Upsilon(4S):withoutneutrino_charged_Jpsi -> B+:feiHadronic B-:sig_Jpsi",cut ="", dmID = 0, path=my_path)
ma.reconstructDecay("Upsilon(4S):withoutneutrino_neutral_opposite_cp_Jpsi -> B0:feiHadronic anti-B0:sig_Jpsi",cut ="", dmID = 1, path=my_path)
ma.reconstructDecay("Upsilon(4S):withoutneutrino_neutral_same_cp_Jpsi -> B0:feiHadronic B0:sig_Jpsi",cut ="", dmID = 1, path=my_path)
ma.copyLists(outputListName="Upsilon(4S):withoutneutrino_Jpsi", inputListNames=["Upsilon(4S):withoutneutrino_charged_Jpsi", "Upsilon(4S):withoutneutrino_neutral_opposite_cp_Jpsi", "Upsilon(4S):withoutneutrino_neutral_same_cp_Jpsi"], path=my_path)

ma.cutAndCopyList('J/psi:myJpsi_proper', 'J/psi:myJpsi','isDescendantOfList(Upsilon(4S):withoutneutrino_Jpsi,2)', path=my_path)
ma.cutAndCopyList('mu+:mychargedMuon_proper', 'mu+:mychargedMuon','isDescendantOfList(J/psi:myJpsi_proper,1)', path=my_path)

# tag side
track_selection = "dr < 0.5 and abs(dz) < 2 and sourceObjectIsInList(mu+:mychargedMuon_proper)==0"
cluster_selection = '[[clusterReg==1 and E>0.08] or \
                      [clusterReg==2 and E>0.03] or \
                      [clusterReg==3 and E>0.06]] and sourceObjectIsInList(mu+:mychargedMuon_proper)==0'
cleanMask = ("cleanMask",track_selection, cluster_selection)
ma.buildRestOfEvent("B+:feiHadronic",path=my_path)
ma.appendROEMasks("B+:feiHadronic",[cleanMask],path=my_path)
ma.buildContinuumSuppression(list_name="B+:feiHadronic",roe_mask = "cleanMask", path=my_path)
ma.buildRestOfEvent("B0:feiHadronic",path=my_path)
ma.appendROEMasks("B0:feiHadronic",[cleanMask],path=my_path)
ma.buildContinuumSuppression(list_name="B0:feiHadronic",roe_mask = "cleanMask", path=my_path)

# signal side
ma.reconstructDecay("B+:sig -> Xsu:comb ?nu", cut="", dmID=0, path = my_path)
ma.reconstructDecay("B0:sig -> Xsd:comb ?nu", cut="", dmID=0, path = my_path)

# Upsilon(4S) -> Btag K
ma.reconstructDecay("Upsilon(4S):withoutneutrino_charged -> B+:feiHadronic B-:sig",cut ="", dmID = 0, path=my_path)
ma.reconstructDecay("Upsilon(4S):withoutneutrino_neutral_opposite_cp -> B0:feiHadronic anti-B0:sig",cut ="", dmID = 1, path=my_path)
ma.reconstructDecay("Upsilon(4S):withoutneutrino_neutral_same_cp -> B0:feiHadronic B0:sig",cut ="", dmID = 1, path=my_path)
ma.copyLists(outputListName="Upsilon(4S):withoutneutrino", inputListNames=["Upsilon(4S):withoutneutrino_charged", "Upsilon(4S):withoutneutrino_neutral_opposite_cp", "Upsilon(4S):withoutneutrino_neutral_same_cp"], path=my_path)
ma.buildRestOfEvent("Upsilon(4S):withoutneutrino", inputParticlelists=["K_S0:myKaonshort", "pi0:myneutralPion", "J/psi:myJpsi_proper"], path=my_path)
ma.applyCuts('Upsilon(4S):withoutneutrino', 'nROE_ParticlesInList(J/psi:myJpsi_proper)==1 and nRemainingTracksInEvent==2', path=my_path)

# apply cut: no charged track, E, cluster on ROE
ma.appendROEMasks("Upsilon(4S):withoutneutrino",[cleanMask],path=my_path)
# ma.applyCuts("Upsilon(4S):withoutneutrino","nROE_Tracks(cleanMask)==0",path=my_path)

# best candidate selection
va.variables.addAlias("SignalProbofBtag_rank","daughter(0,extraInfo(SignalProbability))")
ma.rankByHighest(particleList="Upsilon(4S):withoutneutrino", variable="SignalProbofBtag_rank",allowMultiRank=True,outputVariable="Upsilon_rank",path=my_path)
# ma.applyCuts("Upsilon(4S):withoutneutrino","extraInfo(Upsilon_rank) == 1",path=my_path)

# --- build Event Kinematics ---
# for release-05-02-17
trackCuts = 'pt > 0.1'
trackCuts += ' and thetaInCDCAcceptance'
trackCuts += ' and abs(dz) < 3'
trackCuts += ' and dr < 0.5'
trackCuts += ' and sourceObjectIsInList(mu+:mychargedMuon_proper)==0'

gammaCuts = 'E > 0.05'
gammaCuts += ' and thetaInCDCAcceptance'
gammaCuts += ' and sourceObjectIsInList(mu+:mychargedMuon_proper)==0'

ma.fillParticleList('pi+:evtkin', '', path=my_path)
ma.fillParticleList('gamma:evtkin', '', path=my_path)

ma.applyCuts('pi+:evtkin', trackCuts, path=my_path)
ma.applyCuts('gamma:evtkin', gammaCuts, path=my_path)

particlelistsforevtkinematics = ['pi+:evtkin', 'gamma:evtkin']

ma.buildEventKinematics(path=my_path, inputListNames=particlelistsforevtkinematics)

    
# --- MC truth ---
ma.looseMCTruth(list_name = "B+:feiHadronic", path = my_path)
ma.looseMCTruth(list_name = "B+:sig", path = my_path)
ma.looseMCTruth(list_name = "B0:feiHadronic", path = my_path)
ma.looseMCTruth(list_name = "B0:sig", path = my_path)
ma.looseMCTruth(list_name = "Upsilon(4S):withoutneutrino", path = my_path)

# --- decay string ---
my_path.add_module('ParticleMCDecayString', listName='Upsilon(4S):withoutneutrino', fileName='hashmap_Upsilon.root')

va.variables.addAlias("energy_asymmetry", "abs(formula((daughter(0,E)-daughter(1,E))/(daughter(0,E)+daughter(1,E))))")

# get variables
Kinematics = ["E", "InvM", "M", "p", "phi", "theta", "dM"]
Kinematics_CMS = ["useCMSFrame(E)", "useCMSFrame(p)", "useCMSFrame(phi)", "useCMSFrame(theta)"]
Kinematics_RecoilRestFrame = ["useTagSideRecoilRestFrame(daughter(1,E),0)", "useTagSideRecoilRestFrame(daughter(1,p),0)", "useTagSideRecoilRestFrame(daughter(1,phi),0)", "useTagSideRecoilRestFrame(daughter(1,theta),0)"]
Btag_cut = ["Mbc", "deltaE"]
mcvar = ["mcPDG", "mcE", "mcP", "mcPhi", "mcTheta", "genParticleID", "isSignal", "isSignalAcceptMissingNeutrino", "isSignalAcceptMissingGamma"]
loosemcvar = ["extraInfo(looseMCMotherPDG)", "extraInfo(looseMCWrongDaughterN)", "extraInfo(looseMCWrongDaughterPDG)", "extraInfo(looseMCWrongDaughterBiB)"]
decayhash = ['extraInfo(DecayHash)', 'extraInfo(DecayHashExtended)']
othervar = ["PDG", "extraInfo(decayModeID)"]
continuumsup_vars = ["R2", "thrustBm", "thrustOm", "cosTBTO", "cosTBz", "KSFWVariables(et)", "KSFWVariables(mm2)", "KSFWVariables(hso00)", "KSFWVariables(hso01)", "KSFWVariables(hso02)", "KSFWVariables(hso03)", "KSFWVariables(hso04)", "KSFWVariables(hso10)", "KSFWVariables(hso12)", "KSFWVariables(hso14)", "KSFWVariables(hso20)", "KSFWVariables(hso22)", "KSFWVariables(hso24)", "KSFWVariables(hoo0)", "KSFWVariables(hoo1)", "KSFWVariables(hoo2)", "KSFWVariables(hoo3)", "KSFWVariables(hoo4)", "CleoConeCS(1)", "CleoConeCS(2)", "CleoConeCS(3)", "CleoConeCS(4)", "CleoConeCS(5)", "CleoConeCS(6)", "CleoConeCS(7)", "CleoConeCS(8)", "CleoConeCS(9)"]
EvtKinematics = ["missingMomentumOfEvent", "missingMomentumOfEvent_theta", "missingMass2OfEvent", "visibleEnergyOfEventCMS", "missingEnergyOfEventCMS"]

Btag_vars = vu.create_aliases(list_of_variables = Kinematics + Btag_cut + Kinematics_CMS + othervar + continuumsup_vars + ["chiProb", "extraInfo(SignalProbability)"], wrapper = "daughter(0,{variable})",prefix="Btag")

Bsig_vars = vu.create_aliases(list_of_variables = Kinematics + Kinematics_CMS + othervar + [ "daughter(0, M)", "daughter(0, extraInfo(decayModeID))"], wrapper = "daughter(1,{variable})", prefix="Bsig")

U_vars = Kinematics + Kinematics_CMS + Kinematics_RecoilRestFrame + EvtKinematics + decayhash + othervar + ["extraInfo(Upsilon_rank)",  "nROE_ECLClusters(cleanMask)", "nROE_NeutralECLClusters(cleanMask)", "nROE_KLMClusters", "nROE_Tracks(cleanMask)", "roeEextra(cleanMask)", "roeNeextra(cleanMask)", "useCMSFrame(roeNeextra(cleanMask))", "nROE_ParticlesInList(K_S0:myKaonshort)", "nROE_ParticlesInList(pi0:myneutralPion)", "nRemainingTracksInEvent"]

# Ntuple output
output_file = "Ntuple.root"
ma.variablesToNtuple(decayString="Upsilon(4S):withoutneutrino",variables=Btag_vars,filename=output_file,treename="Btag",path=my_path)
ma.variablesToNtuple(decayString="Upsilon(4S):withoutneutrino",variables=Bsig_vars,filename=output_file,treename="Bsig",path=my_path)
ma.variablesToNtuple(decayString="Upsilon(4S):withoutneutrino",variables=U_vars+["nParticlesInList(J/psi:myJpsi_proper)", "nROE_ParticlesInList(J/psi:myJpsi_proper)"],filename=output_file,treename="Upsilon",path=my_path)
ma.variablesToNtuple(decayString="J/psi:myJpsi_beforedMcut",variables=["M"],filename=output_file,treename="Jpsi",path=my_path)

# for plot of primary particles
#output_file = destination + "/Kaon/" + name+"_Kaon_Ntuple.root"
#ma.variablesToNtuple(decayString="K+:all",variables=["kaonID", "nCDCHits", "dr", "dz", "M", "dM"],filename=output_file,treename="Kaon",path=my_path)

#output_file = destination + "/Pion/" + name+"_Pion_Ntuple.root"
#ma.variablesToNtuple(decayString="pi+:all",variables=["pionID", "nCDCHits", "dr", "dz", "M", "dM"],filename=output_file,treename="Pion",path=my_path)

#output_file = destination + "/Gamma/" + name+"_gamma_Ntuple.root"
#ma.variablesToNtuple(decayString="gamma:all",variables=["inCDCAcceptance", "clusterNHits", "E"],filename=output_file,treename="Gamma",path=my_path)

#ma.fillParticleList(decayString="gamma:temp", cut="inCDCAcceptance and clusterNHits > 1.5 and [[clusterReg==1 and E>0.08] or [clusterReg==2 and E>0.03] or [clusterReg==3 and E>0.06]]",path=my_path)
#ma.reconstructDecay(decayString="pi0:temp -> gamma:temp gamma:temp", cut="", path=my_path)
#output_file = destination + "/Pion0/" + name+"_Pi0_Ntuple.root"
#ma.variablesToNtuple(decayString="pi0:temp",variables=["M", "p", "daughterAngle(0,1)", "daughterDiffOfPhi(0,1)", "dM"],filename=output_file,treename="Pion0",path=my_path)

#ma.reconstructDecay(decayString="K_S0:temp -> pi+:all pi-:all", cut="M > 0.3 and M < 0.7",path=my_path)
#vertex.kFit('K_S0:temp', 0, path=my_path)
#output_file = destination + "/KaonS0/" + name+"_KS0_Ntuple.root"
#ma.variablesToNtuple(decayString="K_S0:temp",variables=["M", "dM", "significanceOfDistance"],filename=output_file,treename="KaonS0",path=my_path)

# progress
basf2.process(my_path)

# Print call statistics
print(basf2.statistics)
