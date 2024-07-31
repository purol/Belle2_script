#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

import operator

unsorted_list = [('Bsig_cosTBTO',0.0720262),
('Bsig_KSFWVariables_hso01',0.00298101),
('Bsig_KSFWVariables_hso04',0.00446538),
('Bsig_thrustBm',0.0435347),
('Bsig_useCMSFrame_p',0.0303834),
('Btag_CleoConeCS_1',0.0125584),
('Btag_CleoConeCS_2',0.00214436),
('Btag_CleoConeCS_3',0.000799951),
('Btag_cosTBTO',0.00591209),
('Btag_KSFWVariables_hoo1',0.0036583),
('Btag_KSFWVariables_hoo2',0.045386),
('Btag_KSFWVariables_hoo3',0.00174413),
('Btag_KSFWVariables_hoo4',0.00194335),
('Btag_KSFWVariables_hso02',0.00422073),
('Btag_KSFWVariables_hso24',0.00387719),
('Btag_useCMSFrame_theta',0.011977),
('extraInfo__boEeclv200__bc',0.335058),
('extraInfo__boNgammav200__bc',0.00425715),
('foxWolframR3',0.00862743),
('foxWolframR4',0.00383641),
('harmonicMomentThrust1',0.00348617),
('harmonicMomentThrust2',0.226642),
('missingMomentumOfEvent',0.108332),
('missingMomentumOfEvent_theta',0.0130551),
('nRemainingTracksInEvent',0.025718),
('roePTheta__bocleanMask__bc',0.00298187),
('Bsig_daughter_0_extraInfo_Dc_pValue_std_noDCS',0.00408601),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_M_noDCS',0.00160124),
('Bsig_daughter_0_extraInfo_D0_pValue_std',0.00780082),
('Bsig_daughter_0_extraInfo_D0simpleveto_M',0.00240111),
('Bsig_daughter_0_extraInfo_D0_pValue_std_yespizero',0.00193056),
('Bsig_daughter_0_extraInfo_D0simpleveto_M_yespizero',0.0025739)]


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
