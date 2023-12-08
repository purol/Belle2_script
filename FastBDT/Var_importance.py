#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

import operator

unsorted_list = [('Bsig_cosTBTO',0.0718164),
('Bsig_KSFWVariables_hso01',0.00286603),
('Bsig_KSFWVariables_hso04',0.0037269),
('Bsig_thrustBm',0.0428134),
('Bsig_useCMSFrame_p',0.0305719),
('Btag_CleoConeCS_1',0.0119502),
('Btag_CleoConeCS_2',0.00179437),
('Btag_CleoConeCS_3',0.000650742),
('Btag_cosTBTO',0.00543584),
('Btag_KSFWVariables_hoo1',0.00343226),
('Btag_KSFWVariables_hoo2',0.0444208),
('Btag_KSFWVariables_hoo3',0.00142186),
('Btag_KSFWVariables_hoo4',0.00252565),
('Btag_KSFWVariables_hso02',0.00356423),
('Btag_KSFWVariables_hso24',0.00298083),
('Btag_useCMSFrame_theta',0.0113786),
('extraInfo__boEeclv200__bc',0.348061),
('extraInfo__boNgammav200__bc',0.00394138),
('foxWolframR3',0.00835863),
('foxWolframR4',0.00313178),
('harmonicMomentThrust1',0.00300595),
('harmonicMomentThrust2',0.225677),
('missingMomentumOfEvent',0.105818),
('missingMomentumOfEvent_theta',0.011467),
('nRemainingTracksInEvent',0.0249758),
('roePTheta__bocleanMask__bc',0.00264763),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_std',0.0104871),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_M',0.00204328),
('Bsig_daughter_0_extraInfo_D0simpleveto_std',0.00609654),
('Bsig_daughter_0_extraInfo_D0simpleveto_M',0.00293917)]

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
