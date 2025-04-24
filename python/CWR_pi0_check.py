#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# usage: basf2 MakeNtuple_multi.py "./20210402/evt-0.mdst"
# last: 2022-07-14-00

import os
import sys

import basf2
import modularAnalysis as ma
import variables as va
import variables.collections as vc
import variables.utils as vu
import vertex

import stdPi0s
import stdV0s
import stdPhotons

from glob import glob

my_path = basf2.create_path()

inputfile="/home/belle2/junewoo/storage_ghi/20220929_SIGNAL_gsim/output/B2Xsnunu_324.root"
ma.inputMdst(environmentType='default',filename=inputfile,path=my_path)

stdPhotons.stdPhotons(listtype="pi0eff30_May2020", path=my_path)
ma.reconstructDecay('pi0:eff30_May2020 -> gamma:pi0eff30_May2020 gamma:pi0eff30_May2020', '', 1, True, my_path)
ma.matchMCTruth('pi0:eff30_May2020', my_path)

ma.variablesToNtuple(decayString="pi0:eff30_May2020",variables=["M", "isSignal"],filename="./pi0_Ntuple.root",treename="pi0",path=my_path)
    
# progress
basf2.process(my_path)

# Print call statistics
print(basf2.statistics)
