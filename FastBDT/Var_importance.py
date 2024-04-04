#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

import operator

unsorted_list = [('Bsig_cosTBTO',0.0697935),
('Bsig_KSFWVariables_hso01',0.00280361),
('Bsig_KSFWVariables_hso04',0.00425243),
('Bsig_thrustBm',0.0434154),
('Bsig_useCMSFrame_p',0.0305483),
('Btag_CleoConeCS_1',0.0118515),
('Btag_CleoConeCS_2',0.0017198),
('Btag_CleoConeCS_3',0.000640396),
('Btag_cosTBTO',0.00697709),
('Btag_KSFWVariables_hoo1',0.00338822),
('Btag_KSFWVariables_hoo2',0.0455039),
('Btag_KSFWVariables_hoo3',0.00143683),
('Btag_KSFWVariables_hoo4',0.0016711),
('Btag_KSFWVariables_hso02',0.00404667),
('Btag_KSFWVariables_hso24',0.00350526),
('Btag_useCMSFrame_theta',0.0117302),
('extraInfo__boEeclv200__bc',0.345956),
('extraInfo__boNgammav200__bc',0.00391882),
('foxWolframR3',0.00837098),
('foxWolframR4',0.00340001),
('harmonicMomentThrust1',0.00322416),
('harmonicMomentThrust2',0.22617),
('missingMomentumOfEvent',0.106035),
('missingMomentumOfEvent_theta',0.012132),
('nRemainingTracksInEvent',0.0250172),
('roePTheta__bocleanMask__bc',0.00276636),
('Bsig_daughter_0_extraInfo_Dc_pValue_std_noDCS',0.00379848),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_M_noDCS',0.00147984),
('Bsig_daughter_0_extraInfo_D0_pValue_std',0.00777465),
('Bsig_daughter_0_extraInfo_D0simpleveto_M',0.0023973),
('Bsig_daughter_0_extraInfo_D0_pValue_std_yespizero',0.00211022),
('Bsig_daughter_0_extraInfo_D0simpleveto_M_yespizero',0.00216408)]


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
