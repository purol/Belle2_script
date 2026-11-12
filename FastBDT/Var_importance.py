#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

import operator

unsorted_list = [('Bsig_cosTBTO',0.0734967),
('Bsig_KSFWVariables_hso01',0.00295744),
('Bsig_KSFWVariables_hso04',0.00439809),
('Bsig_thrustBm',0.0435055),
('Bsig_useCMSFrame_p',0.0309568),
('Btag_CleoConeCS_1',0.0121274),
('Btag_CleoConeCS_2',0.00213454),
('Btag_CleoConeCS_3',0.00077828),
('Btag_cosTBTO',0.00652719),
('Btag_KSFWVariables_hoo1',0.00358606),
('Btag_KSFWVariables_hoo2',0.0465571),
('Btag_KSFWVariables_hoo3',0.00159162),
('Btag_KSFWVariables_hoo4',0.00199616),
('Btag_KSFWVariables_hso02',0.00458749),
('Btag_KSFWVariables_hso24',0.00377825),
('Btag_useCMSFrame_theta',0.012093),
(r'$E_{\mathrm{ecl}}$', 0.333803),
('extraInfo__boNgammav200__bc',0.00413708),
('foxWolframR3',0.0086905),
('foxWolframR4',0.00358842),
('harmonicMomentThrust1',0.00367437),
('harmonicMomentThrust2',0.225409),
('missingMomentumOfEvent',0.107491),
('missingMomentumOfEvent_theta',0.0124422),
('nRemainingTracksInEvent',0.0261395),
('roePTheta__bocleanMask__bc',0.0027493),
('Bsig_daughter_0_extraInfo_Dc_pValue_std_noDCS',0.00399743),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_M_noDCS',0.00162415),
('Bsig_daughter_0_extraInfo_D0_pValue_std',0.00765407),
('Bsig_daughter_0_extraInfo_D0simpleveto_M',0.00243661),
('Bsig_daughter_0_extraInfo_D0_pValue_std_yespizero',0.00241626),
('Bsig_daughter_0_extraInfo_D0simpleveto_M_yespizero',0.00267526)]

sorted_list = sorted(unsorted_list, key=operator.itemgetter(1), reverse=True)

features_sorted = []
importance_sorted = []

for i in sorted_list:
    features_sorted += [i[0]]
    importance_sorted += [i[1]]

features_sorted.reverse()
importance_sorted.reverse()

plt.rcParams.update({
    "font.size": 9,
    "pdf.fonttype": 42,
})

# 논문의 0.95\textwidth에 해당하는 폭으로 생성
fig, ax = plt.subplots(
    figsize=(15.124 / 2.54, 20.5 / 2.54)
)

# 긴 변수명을 위한 왼쪽 공간 확보
fig.subplots_adjust(
    left=0.61,
    right=0.99,
    bottom=0.065,
    top=0.99,
)

bars = ax.barh(
    features_sorted,
    importance_sorted,
    height=0.7,
)

ax.bar_label(
    bars,
    fmt="%.1e",
    padding=3,
    fontsize=9,
)

ax.set_xlim(0, 0.45)
ax.set_xticks([0.0, 0.1, 0.2, 0.3, 0.4])
ax.set_xlabel("Variable importance", fontsize=10)

ax.tick_params(axis="both", labelsize=9)
ax.tick_params(axis="y", length=0)
ax.set_axisbelow(True)
ax.grid(axis="x", linestyle="--", alpha=0.4)

# 지정한 물리적 크기를 유지
fig.savefig("Var_importance.pdf", bbox_inches=None)
plt.close(fig)
