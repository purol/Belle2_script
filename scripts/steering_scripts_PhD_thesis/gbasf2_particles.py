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

import argparse

# get data type
parser = argparse.ArgumentParser(description='Sample type')

parser.add_argument('--inputpath', required=True, type=str, help='Directory containing input ROOT files')
parser.add_argument('--destination', required=False, type=str, help='destination of output files for local run. ex) /home/belle2/junewoo/storage_b1/output')
parser.add_argument( '--filter', required=False, type=str, default="", help='Only read ROOT files containing this string' )

args = parser.parse_args()

basf2.conditions.prepend_globaltag(ma.getAnalysisGlobaltag())

# set random seed
basf2.set_random_seed(42)

# set output name
inputpath = args.inputpath

# Get all ROOT files in the directory
inputfiles = sorted(glob(os.path.join(inputpath, f"*{args.filter}*.root")))

outputpath = args.destination

output_file = outputpath + "/output_Ntuple.root"

# Read uDST
my_path = basf2.create_path()
ma.inputMdstList(environmentType='default', filelist=inputfiles, path=my_path)

# select only 1/100
ma.applyEventCuts("eventRandom < 0.01", path=my_path)

# define basic particles
ma.fillParticleList(decayString="K+:all", cut="", path=my_path)
ma.fillParticleList(decayString="pi+:all", cut="", path=my_path)

ma.reconstructDecay('K_S0:preliminary -> pi+:all pi-:all', '', 1, path=my_path)
vertex.treeFit("K_S0:preliminary", conf_level=0.0, path=my_path, updateAllDaughters=False)

ma.fillParticleList(decayString="gamma:all", cut="", path=my_path)
ma.reconstructDecay('pi0:preliminary -> gamma:all gamma:all', '', 1, path=my_path)

# MC match
ma.matchMCTruth("K+:all", path=my_path)
ma.matchMCTruth("pi+:all", path=my_path)
ma.matchMCTruth("K_S0:preliminary", path=my_path)
ma.matchMCTruth("gamma:all", path=my_path)
ma.matchMCTruth("pi0:preliminary", path=my_path)

# for pion/kaon
ma.variablesToNtuple(decayString="K+:all",variables = ["kaonID", "nCDCHits", "nPXDHits", "dr", "dz", "isSignal", "mcPDG"], filename=output_file, treename="kaon", path=my_path)
ma.variablesToNtuple(decayString="pi+:all",variables = ["pionID", "nCDCHits", "nPXDHits", "dr", "dz", "isSignal", "mcPDG"], filename=output_file, treename="pion", path=my_path)

# for KS0
ma.variablesToNtuple(decayString="K_S0:preliminary",variables = ["M", "significanceOfDistance", "isSignal", "mcPDG"], filename=output_file, treename="KSzero", path=my_path)

# for pi0
var1 = ["daughter(0," + v + ")" for v in ["clusterNHits", "thetaInCDCAcceptance", "clusterReg", "E"]]
var2 = ["daughter(1," + v + ")" for v in ["clusterNHits", "thetaInCDCAcceptance", "clusterReg", "E"]]
ma.variablesToNtuple(decayString="pi0:preliminary",variables = ["InvM", "daughterDiffOf(0,1,phi)", "daughterAngle(0,1)", "p", "isSignal", "mcPDG"] + var1 + var2, filename=output_file, treename="pizero", path=my_path)

# progress
basf2.process(my_path)

# Print call statistics
print(basf2.statistics)
