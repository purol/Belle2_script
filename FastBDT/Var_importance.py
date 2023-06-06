#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

import operator

unsorted_list = [('Bsig_cosTBTO',0.078748),
('Bsig_KSFWVariables_hso01',0.00353357),
('Bsig_KSFWVariables_hso04',0.00462359),
('Bsig_thrustBm',0.0333034),
('Bsig_useCMSFrame_p',0.036978),
('Btag_CleoConeCS_1',0.0055565),
('Btag_CleoConeCS_2',0.00678546),
('Btag_CleoConeCS_3',0.00298066),
('Btag_cosTBTO',0.0433491),
('Btag_KSFWVariables_hoo1',0.00301216),
('Btag_KSFWVariables_hoo2',0.0159541),
('Btag_KSFWVariables_hoo3',0.00192277),
('Btag_KSFWVariables_hoo4',0.00685493),
('Btag_KSFWVariables_hso02',0.00444346),
('Btag_KSFWVariables_hso24',0.0041799),
('Btag_useCMSFrame_theta',0.0125288),
('extraInfo__boEeclv133__bc',0.339821),
('extraInfo__boNgammav133__bc',0.0774836),
('foxWolframR3',0.0154621),
('foxWolframR4',0.00376471),
('harmonicMomentThrust1',0.00236999),
('harmonicMomentThrust2',0.131556),
('missingMomentumOfEvent',0.0933369),
('missingMomentumOfEvent_theta',0.013647),
('nRemainingTracksInEvent',0.023571),
('roePTheta__bocleanMask__bc',0.0101465),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_std',0.0132783),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_M',0.00175872),
('Bsig_daughter_0_extraInfo_D0simpleveto_std',0.00646523),
('Bsig_daughter_0_extraInfo_D0simpleveto_M',0.00258386)]

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
