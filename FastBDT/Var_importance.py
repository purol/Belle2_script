#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

import operator

unsorted_list = [('aplanarity',0.00124887),
('Bsig_daughter_0_extraInfo_mydz',0.00257936),
('Btag_CleoConeCS_1',0.00261798),
('Btag_CleoConeCS_2',0.00320847),
('Btag_CleoConeCS_3',0.00242604),
('Btag_CleoConeCS_4',0.00175438),
('Btag_CleoConeCS_5',0.00120218),
('Btag_CleoConeCS_6',0.00107955),
('Btag_CleoConeCS_7',0.00089965),
('Btag_CleoConeCS_8',0.000581416),
('Btag_cosTBTO',0.031304),
('Btag_cosTBz',0.00516187),
('Btag_deltaE',0.0114402),
('Btag_KSFWVariables_hso01',0.00167718),
('Btag_KSFWVariables_hso03',0.000880145),
('Btag_KSFWVariables_hso04',0.00159479),
('Btag_KSFWVariables_hso14',0.00247632),
('Btag_KSFWVariables_hso24',0.00440373),
('Btag_thrustBm',0.00285456),
('Btag_thrustOm',0.0189037),
('Btag_useCMSFrame_theta',0.00581198),
('cleoConeThrust0',0.000870397),
('cleoConeThrust1',0.000932222),
('cleoConeThrust2',0.00118071),
('cleoConeThrust3',0.00122411),
('cleoConeThrust4',0.00149616),
('cleoConeThrust5',0.00181101),
('cleoConeThrust6',0.00193436),
('cleoConeThrust7',0.00216326),
('cleoConeThrust8',0.00201599),
('harmonicMomentThrust2',0.1024),
('harmonicMomentThrust3',0.00107157),
('harmonicMomentThrust4',0.0126387),
('Btag_extraInfo_SignalProbability',0.0888454),
('missingMomentumOfEvent',0.0924589),
('missingMomentumOfEvent_theta',0.0176805),
('nParticlesInList__bomu__pl__clMuonFBDT_tight__bc',0.000352504),
('roeEextra__bocleanMask__bc',0.53449),
('roePTheta__bocleanMask__bc',0.00616496),
('Bsig_daughter_0_extraInfo_Dc_pValue_std',0.00517205),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_dr',0.00782621),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_dz',0.00105794),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_M',0.0012041),
('Bsig_daughter_0_extraInfo_D0simpleveto_chiProb',0.00794192),
('Bsig_daughter_0_extraInfo_D0simpleveto_dz',0.00133605),
('Bsig_daughter_0_extraInfo_D0simpleveto_M',0.0016248)]

sorted_list = sorted(unsorted_list, key=operator.itemgetter(1), reverse=True)

features_sorted = []
importance_sorted = []

for i in sorted_list:
    features_sorted += [i[0]]
    importance_sorted += [i[1]]

features_sorted.reverse()
importance_sorted.reverse()

plt.rcParams['figure.figsize'] = [15, 12]

fig, ax = plt.subplots()
bars = ax.barh(features_sorted,importance_sorted)
ax.bar_label(bars, fmt = '%.1e')

#plt.barh(y, Nevt)
#plt.yticks(y, decay)

plt.xlabel("Variable Importance")
plt.ylabel("Variable")
plt.grid(True, axis = "x", linestyle="--")

plt.savefig("Var_importance.png", bbox_inches='tight')
