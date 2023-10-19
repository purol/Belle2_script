#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

import operator

unsorted_list = [('Bsig_cosTBTO',0.0348497),
('Bsig_KSFWVariables_hso01',0.0180582),
('Bsig_KSFWVariables_hso04',0.0189861),
('Bsig_thrustBm',0.0384106),
('Bsig_useCMSFrame_p',0.0756693),
('Btag_CleoConeCS_1',0.00859337),
('Btag_CleoConeCS_2',0.00989422),
('Btag_CleoConeCS_3',0.0206479),
('Btag_cosTBTO',0.0220002),
('Btag_KSFWVariables_hoo1',0.0118673),
('Btag_KSFWVariables_hoo2',0.0341887),
('Btag_KSFWVariables_hoo3',0.0139298),
('Btag_KSFWVariables_hoo4',0.0580678),
('Btag_KSFWVariables_hso02',0.0447381),
('Btag_KSFWVariables_hso24',0.0227099),
('Btag_useCMSFrame_theta',0.0257925),
('extraInfo__boEeclv200__bc',0.0512747),
('extraInfo__boNgammav200__bc',0.00309053),
('foxWolframR3',0.0187597),
('foxWolframR4',0.0189489),
('harmonicMomentThrust1',0.0295685),
('harmonicMomentThrust2',0.0154538),
('missingMomentumOfEvent',0.0247489),
('missingMomentumOfEvent_theta',0.0231131),
('nRemainingTracksInEvent',0.017295),
('roePTheta__bocleanMask__bc',0.0150202),
('Bsig_KSFWVariables_hso12',0.0201688),
('Bsig_KSFWVariables_hso14',0.0192086),
('Bsig_KSFWVariables_hso24',0.0603728),
('missingEnergyOfEventCMS',0.0856621),
('foxWolframR1',0.0589581),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_std',0.0200389),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_M',0.0167502),
('Bsig_daughter_0_extraInfo_D0simpleveto_std',0.0255495),
('Bsig_daughter_0_extraInfo_D0simpleveto_M',0.0176142)]

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
