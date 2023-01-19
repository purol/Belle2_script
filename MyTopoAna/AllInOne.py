#!/usr/bin/env python3

import root_pandas
import decayHash
import basf2 as b2
from decayHash import DecayHashMap
import sys
import os
import glob

fName = sys.argv[1] # unnamed root file (merged file)
fFiles_path = sys.argv[2] # path which includes named root file
fhash_path = sys.argv[3] # path which includes hashmap
file_list = os.listdir(fFiles_path)
data = root_pandas.read_root(fName,"Upsilon")

#data = root_pandas.read_root("/home/belle2/junewoo/storage_ghi/20210903_B2Xsnunu_SKIM/Ntuple_for_hashmap/output/Ntuple/B2Xsnunu_669_gsim_SKIM_Ntuple.root", "Upsilon")
#hashmap = DecayHashMap(data, removeRadiativeGammaFlag=False)
#hashmappath="/home/belle2/junewoo/storage_ghi/20210903_B2Xsnunu_SKIM/Ntuple_for_hashmap/output/hashmap/"


if len(data) == 0:
    sys.exit(0)

for i in range(0, len(data)):
    # get one reconstructed J/psi
    find = False
    candidate = data.iloc[i][["extraInfo__boDecayHash__bc", "extraInfo__boDecayHashExtended__bc"]].values
    if not (data.iloc[i]["MVA_Continuum"] > 0.99):
        continue
    print(data.iloc[i]["MVA_Continuum"])   
    print(data.iloc[i]["nRemainingTracksInEvent"])
    print(data.iloc[i]["__experiment__"])
    print(data.iloc[i]["__run__"])
    print(data.iloc[i]["__event__"])
    print(data.iloc[i]["__candidate__"])
    print(data.iloc[i]["__ncandidates__"])

    for j in range(0, len(file_list)): # Ntuple_00000_job212085751_00_before_Eecl_cut_final_output.root
        if find == True:
            break
        data_temp = root_pandas.read_root(os.path.join(sys.argv[2], file_list[j]),"Upsilon")
        if len(data_temp) == 0:
            continue
        for k in range(0, len(data_temp)):
            if data_temp.iloc[k]["__event__"] == data.iloc[i]["__event__"] and data_temp.iloc[k]["__candidate__"] == data.iloc[i]["__candidate__"] and data_temp.iloc[k]["__ncandidates__"] == data.iloc[i]["__ncandidates__"]:
                print(file_list[j])
                for hash_file in glob.glob(fhash_path+"hashmap_Upsilon_"+file_list[j].split('_')[1]+"_*.root"):
                    hashmap = DecayHashMap(hash_file, removeRadiativeGammaFlag=False)

                    # print the reconstruced decay
                    print("Reconstructed Decay: ")
                    rec = hashmap.get_reconstructed_decay(*candidate)
                    print(rec.to_string())

                    # print the original decay as simulated in MC
                    print("Monte Carlo Decay: ")
                    org = hashmap.get_original_decay(*candidate)
                    print(org.to_string())

                print("===============================")
                find = True
                break

# search for a specific decay (sub-decay)
#print("Search for decay:")
#search_decay = decayHash.Belle2.DecayTree('511 (-> 130 (-> -11 11 22) 443)')
#print(search_decay.to_string())
#found = hashmap.get_original_decay(data.iloc[42]["extraInfo__boDecayHash__bc"],
#                                   data.iloc[42]["extraInfo__boDecayHashExtended__bc"]).find_decay(search_decay)
#print("Found: ", found)
sys.exit(0)
