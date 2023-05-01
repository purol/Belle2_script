#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

import operator

unsorted_list = [('Bsig_cosTBTO',0.0572256),
('Bsig_KSFWVariables_hso01',0.0150911),
('Bsig_KSFWVariables_hso04',0.024265),
('Bsig_useCMSFrame_p',0.1499),
('Btag_CleoConeCS_1',0.00449115),
('Btag_CleoConeCS_2',0.0105743),
('Btag_CleoConeCS_3',0.0115667),
('Btag_cosTBTO',0.0143197),
('Btag_KSFWVariables_hoo1',0.0185975),
('Btag_KSFWVariables_hoo2',0.0351522),
('Btag_KSFWVariables_hoo3',0.0316465),
('Btag_KSFWVariables_hoo4',0.0224693),
('Btag_KSFWVariables_hso02',0.0380922),
('Btag_KSFWVariables_hso24',0.0092871),
('Btag_useCMSFrame_theta',0.0158178),
('extraInfo__boEeclv133__bc',0.0726524),
('extraInfo__boNgammav133__bc',0.0292424),
('foxWolframR3',0.00635379),
('foxWolframR4',0.0223789),
('harmonicMomentThrust1',0.0204177),
('harmonicMomentThrust2',0.0143004),
('missingMomentumOfEvent',0.035924),
('missingMomentumOfEvent_theta',0.0230698),
('nRemainingTracksInEvent',0.0286217),
('roePTheta__bocleanMask__bc',0.0463659),
('Bsig_thrustBm',0.14596),
('Btag_useCMSFrame_phi',0.0109773),
('harmonicMomentThrust4',0.0191746),
('nParticlesInList__bomu__pl__clMuonFBDT_loose__bc',0.00615088),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_std',0.0101079),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_M',0.00663294),
('Bsig_daughter_0_extraInfo_D0simpleveto_std',0.0155866),
('Bsig_daughter_0_extraInfo_D0simpleveto_M',0.0275865)]

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
