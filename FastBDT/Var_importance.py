#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

import operator

unsorted_list = [('Bsig_useCMSFrame_p',0.0114112),
('Btag_chiProb',0.0203474),
('Btag_CleoConeCS_1',0.00630871),
('Btag_CleoConeCS_3',0.00279814),
('Btag_CleoConeCS_4',0.00157121),
('Btag_CleoConeCS_5',0.00303248),
('Btag_cosTBz',0.00311278),
('Btag_deltaE',0.00557656),
('Btag_KSFWVariables_hso00',0.011668),
('Btag_KSFWVariables_hso01',0.00402076),
('Btag_KSFWVariables_hso03',0.0026794),
('Btag_KSFWVariables_hso04',0.00374816),
('Btag_KSFWVariables_hso10',0.0301568),
('Btag_KSFWVariables_hso14',0.00521375),
('Btag_KSFWVariables_hso24',0.013396),
('cleoConeThrust2',0.00397513),
('cleoConeThrust3',0.00768356),
('cleoConeThrust4',0.0126673),
('cleoConeThrust5',0.0149928),
('harmonicMomentThrust1',0.027301),
('harmonicMomentThrust3',0.00766604),
('Btag_extraInfo_SignalProbability',0.0300164),
('roeEextra__bocleanMask__bc',0.536936),
('roePextra__bocleanMask__bc',0.00398646),
('Btag_useCMSFrame_theta',0.00301324),
('missingMomentumOfEvent_theta',0.0106413),
('missingEnergyOfEventCMS',0.172411),
('thrustAxisCosTheta',0.00275628),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_chiProb',0.00997659),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_dr',0.0136406),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_M',0.00839424),
('Bsig_daughter_0_extraInfo_D0simpleveto_chiProb',0.00455441),
('Bsig_daughter_0_extraInfo_D0simpleveto_dr',0.00152385),
('Bsig_daughter_0_extraInfo_D0simpleveto_M',0.00282261)]
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
