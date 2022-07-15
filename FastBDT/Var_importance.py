#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

import operator

unsorted_list = [('nRemainingTracksInEvent',0.0788328),
('Btag_chiProb',0.00566332),
('Btag_extraInfo_SignalProbability',0.0753352),
('thrustAxisCosTheta',0.00316563),
('missingMomentumOfEvent_theta',0.0127736),
('Btag_deltaE',0.0124639),
('Btag_useCMSFrame_theta',0.0064941),
('Btag_cosTBTO',0.0698603),
('cleoConeThrust0',0.011121),
('cleoConeThrust1',0.029143),
('cleoConeThrust2',0.00420551),
('cleoConeThrust3',0.00803929),
('cleoConeThrust4',0.00639699),
('cleoConeThrust5',0.00769844),
('cleoConeThrust6',0.0112332),
('cleoConeThrust7',0.0158574),
('cleoConeThrust8',0.0220275),
('harmonicMomentThrust1',0.142754),
('harmonicMomentThrust2',0.252297),
('harmonicMomentThrust3',0.00260351),
('harmonicMomentThrust4',0.0253447),
('sphericity',0.00688973),
('aplanarity',0.0133771),
('Btag_thrustBm',0.0114822),
('roePTheta__bocleanMask__bc',0.00729793),
('Btag_cosTBz',0.00939708),
('Btag_CleoConeCS_1',0.00592369),
('Btag_CleoConeCS_2',0.00899569),
('Btag_CleoConeCS_3',0.0103833),
('Btag_CleoConeCS_4',0.0109701),
('Btag_CleoConeCS_5',0.0141836),
('Btag_CleoConeCS_6',0.0148591),
('Btag_CleoConeCS_7',0.0160743),
('Btag_CleoConeCS_8',0.0164257),
('Btag_CleoConeCS_9',0.0234189),
('Bsig_daughter_0_extraInfo_mychiProb',0.00195392),
('Bsig_daughter_0_extraInfo_mydr',0.00880824),
('Bsig_daughter_0_extraInfo_mydz',0.00172997),
('Bsig_M',0.0145187)]
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
