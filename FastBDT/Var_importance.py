#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

import operator

unsorted_list = [('Bsig_cosTBTO',0.0096505),
('Bsig_KSFWVariables_et',0.0304734),
('Bsig_KSFWVariables_hoo0',0.00802777),
('Bsig_KSFWVariables_hoo1',0.00855805),
('Bsig_KSFWVariables_hoo2',0.016853),
('Bsig_KSFWVariables_hoo3',0.0105011),
('Bsig_KSFWVariables_hoo4',0.0129169),
('Bsig_KSFWVariables_hso00',0.00901067),
('Bsig_KSFWVariables_hso02',0.019749),
('Bsig_KSFWVariables_hso03',0.00994932),
('Bsig_KSFWVariables_hso04',0.00807879),
('Bsig_KSFWVariables_hso01',0.0119061),
('Bsig_KSFWVariables_hso10',0.00559536),
('Bsig_KSFWVariables_hso12',0.0118382),
('Bsig_KSFWVariables_hso14',0.00849833),
('Bsig_KSFWVariables_hso20',0.00870862),
('Bsig_KSFWVariables_hso22',0.0098036),
('Bsig_KSFWVariables_hso24',0.015724),
('Bsig_KSFWVariables_mm2',0.0116566),
('Bsig_thrustBm',0.0135392),
('Bsig_useCMSFrame_p',0.0312252),
('Btag_CleoConeCS_1',0.00480608),
('Btag_CleoConeCS_2',0.00536518),
('Btag_CleoConeCS_3',0.00683842),
('Btag_CleoConeCS_4',0.0090295),
('Btag_CleoConeCS_5',0.00857488),
('Btag_CleoConeCS_6',0.00789467),
('Btag_CleoConeCS_7',0.00697891),
('Btag_CleoConeCS_8',0.00250599),
('Btag_CleoConeCS_9',0.00372308),
('Btag_cosTBTO',0.00935646),
('Btag_KSFWVariables_et',0.0102488),
('Btag_KSFWVariables_hoo0',0.00605639),
('Btag_KSFWVariables_hoo1',0.00240131),
('Btag_KSFWVariables_hoo2',0.0160859),
('Btag_KSFWVariables_hoo3',0.00364906),
('Btag_KSFWVariables_hoo4',0.0471804),
('Btag_KSFWVariables_hso00',0.0163136),
('Btag_KSFWVariables_hso01',0.011839),
('Btag_KSFWVariables_hso02',0.0153719),
('Btag_KSFWVariables_hso03',0.0118221),
('Btag_KSFWVariables_hso04',0.0163075),
('Btag_KSFWVariables_hso10',0.0402382),
('Btag_KSFWVariables_hso12',0.00623092),
('Btag_KSFWVariables_hso14',0.0054562),
('Btag_KSFWVariables_hso20',0.0140618),
('Btag_KSFWVariables_hso22',0.0217588),
('Btag_KSFWVariables_hso24',0.00812044),
('Btag_KSFWVariables_mm2',0.00871863),
('Btag_R2',0.00640313),
('Btag_thrustOm',0.0159889),
('Btag_useCMSFrame_phi',0.0135634),
('Btag_useCMSFrame_theta',0.0143532),
('cleoConeThrust0',0.00220345),
('cleoConeThrust1',0.00952132),
('cleoConeThrust2',0.00831641),
('cleoConeThrust3',0.00705928),
('cleoConeThrust4',0.00833232),
('cleoConeThrust5',0.00524683),
('cleoConeThrust6',0.00703462),
('cleoConeThrust7',0.00702498),
('cleoConeThrust8',0.00388191),
('extraInfo__boEeclv200__bc',0.0138044),
('extraInfo__boNgammav200__bc',0.00171965),
('foxWolframR1',0.0393216),
('foxWolframR2',0.00619402),
('foxWolframR3',0.00925667),
('foxWolframR4',0.00767796),
('harmonicMomentThrust1',0.0144304),
('harmonicMomentThrust2',0.00873244),
('harmonicMomentThrust3',0.0122654),
('harmonicMomentThrust4',0.00740637),
('missingEnergyOfEventCMS',0.0286505),
('missingMass2OfEvent',0.0116804),
('missingMomentumOfEvent',0.00948775),
('missingMomentumOfEvent_theta',0.0114257),
('roePTheta__bocleanMask__bc',0.00774016),
('useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc',0.0136015),
('useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc',0.0278167),
('visibleEnergyOfEventCMS',0.0151187),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_std',0.00695976),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_M',0.00561571),
('Bsig_daughter_0_extraInfo_D0simpleveto_std',0.0147528),
('Bsig_daughter_0_extraInfo_D0simpleveto_M',0.0182138),
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
