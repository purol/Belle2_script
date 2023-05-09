#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

import operator

unsorted_list = [('Bsig_cosTBTO',0.0751316),
('Bsig_KSFWVariables_hso01',0.00421535),
('Bsig_KSFWVariables_hso04',0.00605932),
('Bsig_thrustBm',0.0361318),
('Bsig_useCMSFrame_p',0.0295386),
('Btag_CleoConeCS_1',0.00551212),
('Btag_CleoConeCS_2',0.00550448),
('Btag_CleoConeCS_3',0.00285481),
('Btag_cosTBTO',0.0336962),
('Btag_KSFWVariables_hoo1',0.00382561),
('Btag_KSFWVariables_hoo2',0.0124959),
('Btag_KSFWVariables_hoo3',0.00406581),
('Btag_KSFWVariables_hoo4',0.00473457),
('Btag_KSFWVariables_hso02',0.00346009),
('Btag_KSFWVariables_hso24',0.0058533),
('Btag_useCMSFrame_theta',0.0149139),
('extraInfo__boEeclv133__bc',0.311495),
('extraInfo__boNgammav133__bc',0.111438),
('foxWolframR3',0.0180351),
('foxWolframR4',0.00406593),
('harmonicMomentThrust1',0.00379932),
('harmonicMomentThrust2',0.135232),
('missingMomentumOfEvent',0.100741),
('missingMomentumOfEvent_theta',0.0146081),
('nRemainingTracksInEvent',0.0222412),
('roePTheta__bocleanMask__bc',0.0103653),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_std',0.0072559),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_M',0.00190923),
('Bsig_daughter_0_extraInfo_D0simpleveto_std',0.00695781),
('Bsig_daughter_0_extraInfo_D0simpleveto_M',0.00386383)]

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
