#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

import operator

unsorted_list = [('Bsig_cosTBTO',0.00969085),
('Bsig_KSFWVariables_et',0.0363833),
('Bsig_KSFWVariables_hoo0',0.00911118),
('Bsig_KSFWVariables_hoo1',0.00862658),
('Bsig_KSFWVariables_hoo2',0.0159697),
('Bsig_KSFWVariables_hoo3',0.0082178),
('Bsig_KSFWVariables_hoo4',0.00742839),
('Bsig_KSFWVariables_hso00',0.00664268),
('Bsig_KSFWVariables_hso02',0.024038),
('Bsig_KSFWVariables_hso03',0.00629923),
('Bsig_KSFWVariables_hso04',0.00460015),
('Bsig_KSFWVariables_hso01',0.00744825),
('Bsig_KSFWVariables_hso10',0.00381509),
('Bsig_KSFWVariables_hso12',0.00832681),
('Bsig_KSFWVariables_hso14',0.00576359),
('Bsig_KSFWVariables_hso20',0.0106512),
('Bsig_KSFWVariables_hso22',0.00879741),
('Bsig_KSFWVariables_hso24',0.0205911),
('Bsig_KSFWVariables_mm2',0.00784157),
('Bsig_thrustBm',0.0245765),
('Bsig_useCMSFrame_p',0.0468444),
('Btag_CleoConeCS_1',0.00429733),
('Btag_CleoConeCS_2',0.00618423),
('Btag_CleoConeCS_3',0.00557909),
('Btag_CleoConeCS_4',0.00593665),
('Btag_CleoConeCS_5',0.00481497),
('Btag_CleoConeCS_6',0.00570199),
('Btag_CleoConeCS_7',0.00730373),
('Btag_CleoConeCS_8',0.00124661),
('Btag_CleoConeCS_9',0.00542528),
('Btag_cosTBTO',0.00872757),
('Btag_KSFWVariables_et',0.0115984),
('Btag_KSFWVariables_hoo0',0.00480045),
('Btag_KSFWVariables_hoo1',0.00169186),
('Btag_KSFWVariables_hoo2',0.0125954),
('Btag_KSFWVariables_hoo3',0.00359275),
('Btag_KSFWVariables_hoo4',0.0577087),
('Btag_KSFWVariables_hso00',0.016236),
('Btag_KSFWVariables_hso01',0.00900101),
('Btag_KSFWVariables_hso02',0.0144893),
('Btag_KSFWVariables_hso03',0.00794642),
('Btag_KSFWVariables_hso04',0.00976658),
('Btag_KSFWVariables_hso10',0.0596972),
('Btag_KSFWVariables_hso12',0.00356017),
('Btag_KSFWVariables_hso14',0.00494943),
('Btag_KSFWVariables_hso20',0.0146912),
('Btag_KSFWVariables_hso22',0.02428),
('Btag_KSFWVariables_hso24',0.00707394),
('Btag_KSFWVariables_mm2',0.00710154),
('Btag_R2',0.00364589),
('Btag_thrustOm',0.0204551),
('Btag_useCMSFrame_phi',0.00714375),
('Btag_useCMSFrame_theta',0.017039),
('cleoConeThrust0',0.00144683),
('cleoConeThrust1',0.00738887),
('cleoConeThrust2',0.0055581),
('cleoConeThrust3',0.00575368),
('cleoConeThrust4',0.00550703),
('cleoConeThrust5',0.003986),
('cleoConeThrust6',0.00552016),
('cleoConeThrust7',0.00705725),
('cleoConeThrust8',0.00275188),
('extraInfo__boEeclv200__bc',0.0187943),
('extraInfo__boNgammav200__bc',0.000416692),
('foxWolframR1',0.0481476),
('foxWolframR2',0.00447422),
('foxWolframR3',0.00526971),
('foxWolframR4',0.00609454),
('harmonicMomentThrust1',0.00984696),
('harmonicMomentThrust2',0.00508836),
('harmonicMomentThrust3',0.00420304),
('harmonicMomentThrust4',0.00806868),
('missingEnergyOfEventCMS',0.0549729),
('missingMass2OfEvent',0.0136928),
('missingMomentumOfEvent',0.0066552),
('missingMomentumOfEvent_theta',0.0087888),
('roePTheta__bocleanMask__bc',0.00680361),
('useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc',0.0191019),
('useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc',0.0261817),
('visibleEnergyOfEventCMS',0.0130431),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_std',0.00593163),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_M',0.00437127),
('Bsig_daughter_0_extraInfo_D0simpleveto_std',0.0129757),
('Bsig_daughter_0_extraInfo_D0simpleveto_M',0.0161624),
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
