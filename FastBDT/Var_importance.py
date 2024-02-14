#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

import operator

unsorted_list = [('Bsig_cosTBTO',0.0101986),
('Bsig_KSFWVariables_et',0.0214735),
('Bsig_KSFWVariables_hoo0',0.0082472),
('Bsig_KSFWVariables_hoo1',0.00950944),
('Bsig_KSFWVariables_hoo2',0.0147349),
('Bsig_KSFWVariables_hoo3',0.0123265),
('Bsig_KSFWVariables_hoo4',0.0156334),
('Bsig_KSFWVariables_hso00',0.0125506),
('Bsig_KSFWVariables_hso02',0.0151078),
('Bsig_KSFWVariables_hso03',0.00627845),
('Bsig_KSFWVariables_hso04',0.0101944),
('Bsig_KSFWVariables_hso01',0.0145051),
('Bsig_KSFWVariables_hso10',0.0121506),
('Bsig_KSFWVariables_hso12',0.0127805),
('Bsig_KSFWVariables_hso14',0.0144483),
('Bsig_KSFWVariables_hso20',0.0142836),
('Bsig_KSFWVariables_hso22',0.0118136),
('Bsig_KSFWVariables_hso24',0.0162346),
('Bsig_KSFWVariables_mm2',0.0166195),
('Bsig_thrustBm',0.0127134),
('Bsig_useCMSFrame_p',0.0235164),
('Btag_CleoConeCS_1',0.00473604),
('Btag_CleoConeCS_2',0.0100519),
('Btag_CleoConeCS_3',0.00586553),
('Btag_CleoConeCS_4',0.00588283),
('Btag_CleoConeCS_5',0.00542218),
('Btag_CleoConeCS_6',0.00812838),
('Btag_CleoConeCS_7',0.0124335),
('Btag_CleoConeCS_8',0.00715929),
('Btag_CleoConeCS_9',0.00948522),
('Btag_cosTBTO',0.00391538),
('Btag_KSFWVariables_et',0.0133768),
('Btag_KSFWVariables_hoo0',0.00678459),
('Btag_KSFWVariables_hoo1',0.0118354),
('Btag_KSFWVariables_hoo2',0.0183833),
('Btag_KSFWVariables_hoo3',0.0132279),
('Btag_KSFWVariables_hoo4',0.0157104),
('Btag_KSFWVariables_hso00',0.0195857),
('Btag_KSFWVariables_hso01',0.0127986),
('Btag_KSFWVariables_hso02',0.0154957),
('Btag_KSFWVariables_hso03',0.0104841),
('Btag_KSFWVariables_hso04',0.0191875),
('Btag_KSFWVariables_hso10',0.0167085),
('Btag_KSFWVariables_hso12',0.00631565),
('Btag_KSFWVariables_hso14',0.00797747),
('Btag_KSFWVariables_hso20',0.0223386),
('Btag_KSFWVariables_hso22',0.0141938),
('Btag_KSFWVariables_hso24',0.0109401),
('Btag_KSFWVariables_mm2',0.00696028),
('Btag_R2',0.00837291),
('Btag_thrustOm',0.012666),
('Btag_useCMSFrame_phi',0.00981679),
('Btag_useCMSFrame_theta',0.0147375),
('cleoConeThrust0',0.00597869),
('cleoConeThrust1',0.00753636),
('cleoConeThrust2',0.00867155),
('cleoConeThrust3',0.00770377),
('cleoConeThrust4',0.00789413),
('cleoConeThrust5',0.00880025),
('cleoConeThrust6',0.00799151),
('cleoConeThrust7',0.00548562),
('cleoConeThrust8',0.00636177),
('extraInfo__boEeclv200__bc',0.012725),
('extraInfo__boNgammav200__bc',0.00166066),
('foxWolframR1',0.0334162),
('foxWolframR2',0.00521502),
('foxWolframR3',0.0128183),
('foxWolframR4',0.0117364),
('harmonicMomentThrust1',0.0115549),
('harmonicMomentThrust2',0.00589356),
('harmonicMomentThrust3',0.00946979),
('harmonicMomentThrust4',0.00492167),
('missingEnergyOfEventCMS',0.0253359),
('missingMass2OfEvent',0.0118035),
('missingMomentumOfEvent',0.0110275),
('missingMomentumOfEvent_theta',0.0162576),
('roePTheta__bocleanMask__bc',0.00664797),
('useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc',0.016991),
('useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc',0.0258085),
('visibleEnergyOfEventCMS',0.0137891),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_std',0.0123737),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_M',0.010653),
('Bsig_daughter_0_extraInfo_D0simpleveto_std',0.00904251),
('Bsig_daughter_0_extraInfo_D0simpleveto_M',0.018136),
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
