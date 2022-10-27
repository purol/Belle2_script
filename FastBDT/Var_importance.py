#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

import operator

unsorted_list = [('aplanarity',0.000971157),
('Bsig_daughter_0_extraInfo_mydz',0.00257489),
('Btag_CleoConeCS_1',0.00219841),
('Btag_CleoConeCS_2',0.00279873),
('Btag_CleoConeCS_3',0.00196512),
('Btag_CleoConeCS_4',0.000903873),
('Btag_CleoConeCS_5',0.00062769),
('Btag_cosTBTO',0.0639741),
('Btag_deltaE',0.0169642),
('Btag_KSFWVariables_et',0.00303896),
('Btag_KSFWVariables_hoo2',0.02285),
('Btag_KSFWVariables_hoo4',0.00852612),
('Btag_KSFWVariables_hso03',0.000621794),
('Btag_KSFWVariables_hso14',0.00172002),
('Btag_KSFWVariables_hso20',0.0167142),
('Btag_KSFWVariables_hso22',0.0181347),
('Btag_KSFWVariables_hso24',0.00410001),
('Btag_R2',0.0330778),
('Btag_thrustBm',0.0058272),
('Btag_thrustOm',0.00743986),
('Btag_useCMSFrame_theta',0.0076057),
('cleoConeThrust0',0.000765439),
('cleoConeThrust1',0.00031088),
('cleoConeThrust2',0.00031526),
('cleoConeThrust3',0.000312257),
('cleoConeThrust4',0.000294896),
('cleoConeThrust5',0.000265724),
('cleoConeThrust6',0.000315318),
('cleoConeThrust8',0.000398141),
('foxWolframR2',0.00197906),
('foxWolframR3',0.00687666),
('harmonicMomentThrust3',0.000860169),
('harmonicMomentThrust4',0.0611736),
('Btag_extraInfo_SignalProbability',0.127427),
('missingEnergyOfEventCMS',0.0622492),
('missingMomentumOfEvent',0.0713292),
('missingMomentumOfEvent_theta',0.0134532),
('nParticlesInList__bomu__pl__clMuonFBDT_tight__bc',0.00151885),
('roeEextra__bocleanMask__bc',0.410296),
('roePTheta__bocleanMask__bc',0.00317799),
('sphericity',0.00765201),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_dz',0.00200869),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_M',0.00105562),
('Bsig_daughter_0_extraInfo_D0simpleveto_dz',0.00159548),
('Bsig_daughter_0_extraInfo_D0simpleveto_M',0.00173521)]

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
