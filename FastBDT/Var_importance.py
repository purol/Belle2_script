#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

import operator

unsorted_list = [('aplanarity',0.000830879),
('Bsig_daughter_0_extraInfo_mydz',0.0020833),
('Btag_CleoConeCS_1',0.00161476),
('Btag_CleoConeCS_2',0.00326914),
('Btag_CleoConeCS_3',0.002454),
('Btag_CleoConeCS_4',0.00105397),
('Btag_CleoConeCS_5',0.000701435),
('Btag_cosTBTO',0.0890277),
('Btag_deltaE',0.0202895),
('Btag_KSFWVariables_hso03',0.000916079),
('Btag_KSFWVariables_hso14',0.00204479),
('Btag_KSFWVariables_hso24',0.00455722),
('Btag_thrustBm',0.0149471),
('Btag_thrustOm',0.021214),
('Btag_useCMSFrame_theta',0.00833965),
('cleoConeThrust0',0.000630619),
('cleoConeThrust1',0.000461592),
('cleoConeThrust2',0.000799499),
('cleoConeThrust3',0.000526804),
('cleoConeThrust4',0.000761998),
('cleoConeThrust5',0.000957679),
('cleoConeThrust6',0.00140398),
('cleoConeThrust7',0.000842434),
('cleoConeThrust8',0.00027346),
('harmonicMomentThrust3',0.000853729),
('harmonicMomentThrust4',0.0940929),
('Btag_extraInfo_SignalProbability',0.132696),
('missingEnergyOfEventCMS',0.0697817),
('missingMomentumOfEvent',0.0770795),
('missingMomentumOfEvent_theta',0.0166417),
('nParticlesInList__bomu__pl__clMuonFBDT_tight__bc',0.00273364),
('roeEextra__bocleanMask__bc',0.408655),
('roePTheta__bocleanMask__bc',0.00307813),
('thrustAxisCosTheta',0.00596359),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_dz',0.00235419),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_M',0.00249352),
('Bsig_daughter_0_extraInfo_D0simpleveto_dz',0.00178637),
('Bsig_daughter_0_extraInfo_D0simpleveto_M',0.00178797)]

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
