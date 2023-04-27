#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

import operator

unsorted_list = [('aplanarity',0.00116713),
('Bsig_cosTBTO',0.0267418),
('Bsig_KSFWVariables_hso01',0.0105757),
('Bsig_useCMSFrame_p',0.165944),
('Btag_chiProb',0.00066659),
('Btag_CleoConeCS_1',0.00270384),
('Btag_CleoConeCS_4',0.00129182),
('Btag_CleoConeCS_5',0.00116388),
('Btag_KSFWVariables_hoo1',0.0620366),
('Btag_KSFWVariables_hoo2',0.0414117),
('Btag_KSFWVariables_hoo3',0.0407185),
('Btag_KSFWVariables_hoo4',0.139524),
('Btag_KSFWVariables_hso00',0.0020964),
('Btag_KSFWVariables_hso02',0.0249966),
('Btag_KSFWVariables_hso24',0.00162291),
('Btag_useCMSFrame_theta',0.00119127),
('extraInfo__boEeclv133__bc',0.0591407),
('extraInfo__boNgammav133__bc',0.132281),
('foxWolframR3',0.00381777),
('foxWolframR4',0.00937896),
('harmonicMomentThrust1',0.00185911),
('harmonicMomentThrust2',0.00120147),
('Btag_extraInfo_SignalProbability',0.117097),
('missingMomentumOfEvent',0.031783),
('missingMomentumOfEvent_theta',0.00238295),
('nRemainingTracksInEvent',0.0344527),
('roePTheta__bocleanMask__bc',0.0514964),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_std',0.00374422),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_M',0.00236372),
('Bsig_daughter_0_extraInfo_D0simpleveto_std',0.0232902),
('Bsig_daughter_0_extraInfo_D0simpleveto_M',0.00185735),
]

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
