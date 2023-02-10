#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

import operator

unsorted_list = [('aplanarity',0.00165872),
('Bsig_cosTBTO',0.0874969),
('Bsig_KSFWVariables_hoo1',0.00183806),
('Bsig_KSFWVariables_hso01',0.00275083),
('Bsig_KSFWVariables_hso04',0.0146183),
('Bsig_useCMSFrame_p',0.0220495),
('Btag_chiProb',0.00337943),
('Btag_CleoConeCS_1',0.0217006),
('Btag_CleoConeCS_2',0.0274642),
('Btag_CleoConeCS_3',0.0117423),
('Btag_CleoConeCS_4',0.00559127),
('Btag_CleoConeCS_5',0.00281039),
('Btag_cosTBTO',0.0644747),
('Btag_KSFWVariables_hoo1',0.00361151),
('Btag_KSFWVariables_hoo2',0.0194036),
('Btag_KSFWVariables_hoo3',0.00333259),
('Btag_KSFWVariables_hoo4',0.00658666),
('Btag_KSFWVariables_hso00',0.00370252),
('Btag_KSFWVariables_hso01',0.00185574),
('Btag_KSFWVariables_hso02',0.012249),
('Btag_KSFWVariables_hso24',0.00405913),
('Btag_useCMSFrame_theta',0.00883227),
('extraInfo__boEeclv133__bc',0.322912),
('extraInfo__boNgammav133__bc',0.051052),
('foxWolframR3',0.0346456),
('foxWolframR4',0.00341475),
('harmonicMomentThrust1',0.00216114),
('harmonicMomentThrust3',0.00257841),
('Btag_extraInfo_SignalProbability',0.125475),
('missingMomentumOfEvent',0.0882783),
('missingMomentumOfEvent_theta',0.0120194),
('nParticlesInList__bomu__pl__clMuonFBDT_tight__bc',0.00196066),
('roePTheta__bocleanMask__bc',0.00745622),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_std',0.00886326),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_M',0.00182483),
('Bsig_daughter_0_extraInfo_D0simpleveto_std',0.00408337),
('Bsig_daughter_0_extraInfo_D0simpleveto_M',0.00206687)]

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
