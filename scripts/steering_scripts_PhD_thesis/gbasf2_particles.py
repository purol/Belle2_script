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

parser.add_argument('--inputfile', required=False, type=str, help='destination of output files for local run. ex) /home/belle2/junewoo/storage_b1/input_Jpsi/B02XJpsi_SKIM_00012_job227326565_00.udst.root')
parser.add_argument('--destination', required=False, type=str, help='destination of output files for local run. ex) /home/belle2/junewoo/storage_b1/output')

args = parser.parse_args()

basf2.conditions.prepend_globaltag(ma.getAnalysisGlobaltag())

# set random seed
basf2.set_random_seed(42)

# set output name
inputfile = args.inputfile

basename = os.path.basename(inputfile)
name = os.path.splitext(basename)[0]

outputpath = args.destination
output_file = "Ntuple"
output_file = outputpath + "/" + name + "_" + output_file + ".root"

# Read uDST
my_path = basf2.create_path()
ma.inputMdst(environmentType='default', filename=inputfile, path=my_path)

# select only 1/100
ma.applyEventCuts("eventRandom < 0.1", path=my_path)

# define basic particles
ma.fillParticleList(decayString="K+:all", cut="", path=my_path)
ma.fillParticleList(decayString="pi+:all", cut="", path=my_path)

ma.reconstructDecay('K_S0:preliminary -> pi+:all pi-:all', cut = '', ignoreIfTooManyCandidates=False, path=my_path)
vertex.treeFit("K_S0:preliminary", conf_level=0.0, path=my_path, updateAllDaughters=False)

ma.fillParticleList(decayString="gamma:all", cut="", path=my_path)

stdPhotons.stdPhotons('pi0eff30_May2020', path=my_path)
ma.reconstructDecay('pi0:preliminary -> gamma:pi0eff30_May2020 gamma:pi0eff30_May2020', cut = '', ignoreIfTooManyCandidates=False, path=my_path)

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

# for gamma
ma.variablesToNtuple(decayString="gamma:all",variables = ["clusterNHits", "thetaInCDCAcceptance", "clusterReg", "E", "theta", "isSignal", "mcPDG"], filename=output_file, treename="gamma", path=my_path)

# for pi0
ma.variablesToNtuple(decayString="pi0:preliminary",variables = ["InvM", "daughterDiffOf(0,1,phi)", "daughterAngle(0,1)", "p", "isSignal", "mcPDG"], filename=output_file, treename="pizero", path=my_path)

# progress
basf2.process(my_path)

# Print call statistics
print(basf2.statistics)