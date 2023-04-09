#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

import operator

unsorted_list = [('Bsig_KSFWVariables_hso04',0.00525997),
('Bsig_cosTBTO',0.0617019),
('Bsig_useCMSFrame_p',0.0263803),
('Btag_CleoConeCS_1',0.00319137),
('Btag_CleoConeCS_2',0.00406663),
('Btag_CleoConeCS_3',0.00236327),
('Btag_KSFWVariables_hoo2',0.0154766),
('Btag_KSFWVariables_hoo4',0.0134585),
('Btag_cosTBTO',0.0265327),
('Btag_extraInfo_SignalProbability',0.130505),
('Btag_useCMSFrame_theta',0.00830994),
('extraInfo__boEeclv133__bc',0.3376),
('extraInfo__boNgammav133__bc',0.0535455),
('foxWolframR3',0.0146698),
('harmonicMomentThrust2',0.147624),
('missingMomentumOfEvent',0.0874751),
('missingMomentumOfEvent_theta',0.0130407),
('nRemainingTracksInEvent',0.0194831),
('roePTheta__bocleanMask__bc',0.00801542),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_std',0.0112886),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_M',0.0021517),
('Bsig_daughter_0_extraInfo_D0simpleveto_std',0.00568092),
('Bsig_daughter_0_extraInfo_D0simpleveto_M',0.00217899)]

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
