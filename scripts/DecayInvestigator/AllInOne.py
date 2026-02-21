#!/usr/bin/env python3

import root_pandas
import decayHash
import basf2 as b2
from decayHash import DecayHashMap
import sys
import os
import glob

fFiles_path = sys.argv[1] # path which includes named root file
fhash_path = sys.argv[2] # path which includes hashmap
file_list = os.listdir(fFiles_path)

#data = root_pandas.read_root("/home/belle2/junewoo/storage_ghi/20210903_B2Xsnunu_SKIM/Ntuple_for_hashmap/output/Ntuple/B2Xsnunu_669_gsim_SKIM_Ntuple.root", "Upsilon")
#hashmap = DecayHashMap(data, removeRadiativeGammaFlag=False)
#hashmappath="/home/belle2/junewoo/storage_ghi/20210903_B2Xsnunu_SKIM/Ntuple_for_hashmap/output/hashmap/"

for i in range(0, len(file_list)):
    if not file_list[i].endswith(".root"):
        continue
    data = root_pandas.read_root(os.path.join(fFiles_path, file_list[i]),"Upsilon")
    if len(data) == 0:
        continue

    hashfile = glob.glob(fhash_path+"*_"+file_list[i].split('_')[4]+"_*.root")[0]
    hashmap = DecayHashMap(hashfile, removeRadiativeGammaFlag=True)

    for j in range(0, len(data)):
        candidate = data.iloc[j][["extraInfo__boDecayHash__bc", "extraInfo__boDecayHashExtended__bc"]].values
        print(data.iloc[j]["__experiment__"])
        print(data.iloc[j]["__run__"])
        print(data.iloc[j]["__event__"])
        print(data.iloc[j]["__candidate__"])
        print(data.iloc[j]["__ncandidates__"])

        # print the reconstruced decay
        print("Reconstructed Decay: ")
        rec = hashmap.get_reconstructed_decay(*candidate)
        print(rec.to_string())

        # print the original decay as simulated in MC
        print("Monte Carlo Decay: ")
        org = hashmap.get_original_decay(*candidate)
        print(org.to_string())

        print("===============================")

# search for a specific decay (sub-decay)
#print("Search for decay:")
#search_decay = decayHash.Belle2.DecayTree('511 (-> 130 (-> -11 11 22) 443)')
#print(search_decay.to_string())
#found = hashmap.get_original_decay(data.iloc[42]["extraInfo__boDecayHash__bc"],
#                                   data.iloc[42]["extraInfo__boDecayHashExtended__bc"]).find_decay(search_decay)
#print("Found: ", found)
sys.exit(0)
