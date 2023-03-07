#!/usr/bin/env python
import numpy as np
import sys, os
import math

import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt

import operator

unsorted_list = [('aplanarity',0.00263161),
('Bsig_cosTBTO',0.0573381),
('Bsig_KSFWVariables_hso01',0.00246914),
('Bsig_KSFWVariables_hso04',0.00483575),
('Bsig_useCMSFrame_p',0.0225347),
('Btag_chiProb',0.00337531),
('Btag_CleoConeCS_1',0.00524673),
('Btag_CleoConeCS_2',0.00706932),
('Btag_CleoConeCS_3',0.00425641),
('Btag_CleoConeCS_4',0.00255),
('Btag_CleoConeCS_5',0.0017607),
('Btag_cosTBTO',0.0278043),
('Btag_KSFWVariables_hoo1',0.00286039),
('Btag_KSFWVariables_hoo2',0.0118566),
('Btag_KSFWVariables_hoo3',0.00343706),
('Btag_KSFWVariables_hoo4',0.0108072),
('Btag_KSFWVariables_hso00',0.00302907),
('Btag_KSFWVariables_hso02',0.00346785),
('Btag_KSFWVariables_hso24',0.00389539),
('Btag_useCMSFrame_theta',0.00820106),
('extraInfo__boEeclv133__bc',0.301237),
('extraInfo__boNgammav133__bc',0.059394),
('foxWolframR3',0.0138706),
('foxWolframR4',0.00246731),
('harmonicMomentThrust1',0.00203473),
('harmonicMomentThrust2',0.150477),
('harmonicMomentThrust3',0.00195767),
('Btag_extraInfo_SignalProbability',0.129111),
('missingMomentumOfEvent',0.0875447),
('missingMomentumOfEvent_theta',0.0135882),
('nParticlesInList__bomu__pl__clMuonFBDT_tight__bc',0.00234836),
('nRemainingTracksInEvent',0.0198082),
('roePTheta__bocleanMask__bc',0.00678574),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_std',0.0113299),
('Bsig_daughter_0_extraInfo_Dcsimpleveto_M',0.00184014),
('Bsig_daughter_0_extraInfo_D0simpleveto_std',0.00456584),
('Bsig_daughter_0_extraInfo_D0simpleveto_M',0.00221283)]

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
