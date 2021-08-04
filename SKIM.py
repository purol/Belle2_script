#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# usage: basf2 MakeNtuple_multi.py "./20210402/evt-0.mdst"
# last: 2021-08-03-0

import os
import sys

import basf2
import modularAnalysis as ma
import variables as va
import variables.collections as vc
import variables.utils as vu
from vertex import kFit
from b2biiConversion import convertBelleMdstToBelleIIMdst
from b2biiMonitors import addBeamParamsConversionMonitors
from b2biiMonitors import addTrackConversionMonitors
from b2biiMonitors import addNeutralsConversionMonitors

import fei
from glob import glob

os.environ['USE_GRAND_REPROCESS_DATA'] = '1'
os.environ['PGUSER'] = 'g0db'

#basf2.conditions.prepend_globaltag("analysis_tools_light-2104-poseidon data_reprocessing_proc9 online")
basf2.conditions.prepend_globaltag(ma.getAnalysisGlobaltag())

# set random seed
basf2.set_random_seed(42)

# If you want to create the monitoring histograms (recommended in the beginning of your analysis), just provide any argument
monitoring = False
if len(sys.argv) == 1:
    basf2.B2WARNING("In the beginning of an analysis it is recommended to study the monitoring histograms.\n"
                    "These tell you if the conversion works as expected.\n"
                    "If you want to create them, just provide any argument to this script.")
else:
    monitoring = False

fName = sys.argv[1]
destination = sys.argv[2]

if not fName.endswith(".mdst"): sys.exit(1)
else:
    basename = os.path.basename(fName)
    name = os.path.splitext(basename)[0]
print("=================="+fName+" is conducted=================")


# Convert
my_path = basf2.create_path()
inputfile=fName

convertBelleMdstToBelleIIMdst(inputfile, applySkim=True, enableLocalDB=False, path=my_path)

# Reconstruct
if monitoring:
    # Create monitoring histograms if requested
    addBeamParamsConversionMonitors(path=my_path)
    addTrackConversionMonitors(path=my_path)
    addNeutralsConversionMonitors(path=my_path)

# ==============================
# FEI
# ==============================

particles = fei.get_default_channels(semileptonic=False)
configuration = fei.config.FeiConfiguration(prefix="FEI_B2BII_light-2012-minos",training=False, monitor=False)
feistate = fei.get_path(particles, configuration)
my_path.add_path(feistate.path)

# tag side
ma.applyCuts("B+:generic","Mbc>5.2 and abs(deltaE)<0.5", path=my_path)
ma.applyCuts("B+:generic","extraInfo(SignalProbability) > 0.001",path=my_path)
ma.applyCuts("B0:generic","Mbc>5.2 and abs(deltaE)<0.5", path=my_path)
ma.applyCuts("B0:generic","extraInfo(SignalProbability) > 0.001",path=my_path)

# clean mask
track_selection = "dr < 2 and abs(dz) < 4"
cluster_selection = '[[E > 0.10 and formula(theta/3.14*180) < 31.4] or \
                      [E > 0.05 and formula(theta/3.14*180) > 31.4 and formula(theta/3.14*180) < 130.7] or \
                      [E > 0.15 and formula(theta/3.14*180) > 130.7]]'
cleanMask = ("cleanMask",track_selection, cluster_selection)

ma.buildRestOfEvent("B+:generic", path=my_path)
ma.appendROEMasks("B+:generic",[cleanMask],path=my_path)
ma.buildRestOfEvent("B0:generic", path=my_path)
ma.appendROEMasks("B0:generic",[cleanMask],path=my_path)

ma.applyCuts("B+:generic","nROE_Tracks(cleanMask) < 3",path=my_path)
ma.applyCuts("B0:generic","nROE_Tracks(cleanMask) < 3",path=my_path)

ma.applyEventCuts("nParticlesInList(B+:generic) > 0 or nParticlesInList(B0:generic) > 0", path=my_path)

ma.removeParticlesNotInLists(lists_to_keep=["B+:generic", "B0:generic", "K_S0:mdst", "pi0:mdst", "gamma:mdst", "K_L0:mdst", "Lambda0:mdst", "gamma:v0mdst"], path=my_path)

output_file = destination + "/" + name+"_SKIM.root"
ma.outputUdst(filename=output_file,particleLists=["B+:generic", "B0:generic", "K_S0:mdst", "pi0:mdst", "gamma:mdst", "K_L0:mdst", "Lambda0:mdst", "gamma:v0mdst"],path=my_path)
# progress
basf2.process(my_path)

# Print call statistics
print(basf2.statistics)
