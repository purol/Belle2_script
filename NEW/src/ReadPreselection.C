// last update: 2021-10-08
// for Belle2 data

/*
when you add new variables:
revise # define N_Needed_info ...
revise typedef struct data
revise void Loader::GetData(TFile* input_file)
revise void Loader::PrintRootFile(std::string output_name)
revise void Loader::PrintSeparateRootFile(std::string output_name)
revise void Loader::ConvertIntoSeparateDataFile(std::string output_name, double flag = 0)
*/

#include <iostream>
#include <stdio.h>
#include <string>
#include <algorithm>
#include <vector>
#include <numeric>
#include <time.h>
#include <random>
#include <sstream>
#include <queue>

#include "TH1.h"
#include "TH2.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TCollection.h"
#include "TPad.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TMath.h"

#include "constants.h"
#include "base.h"
#include "ObtainWeight.h"
#include "correctors.h"

/*
when you add new variables:
revise # define N_Needed_info ...
revise typedef struct dataLoader::ConvertIntoSeparateDataFile(std::string output_name, double flag = 0)
revise void Loader::GetData(TFile* input_file)
check MakeShiftDoubleToInt if it affects the code
revise void Loader::PrintRootFile(std::string output_name)
revise void Loader::PrintSeparateRootFile(std::string output_name)
revise void Loader::ConvertIntoSeparateDataFile(std::string output_name, double flag = 0)
*/

typedef struct data {
    int __experiment__;
    int __run__;
    unsigned int __event__;
    int __candidate__;
    int __ncandidates__;

    double npi0;

} Data;

std::queue<Data> TotalData;

void Initialize() {
    std::queue<Data> empty;
    std::swap(TotalData, empty);
}

void ReadFiles(TFile* input_file) {
    TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon_preselection");

    Data temp = { 0 };

    tree_upsilon->SetBranchAddress("__experiment__", &temp.__experiment__);
    tree_upsilon->SetBranchAddress("__run__", &temp.__run__);
    tree_upsilon->SetBranchAddress("__event__", &temp.__event__);
    tree_upsilon->SetBranchAddress("__candidate__", &temp.__candidate__);
    tree_upsilon->SetBranchAddress("__ncandidates__", &temp.__ncandidates__);

    tree_upsilon->SetBranchAddress("nROE_ParticlesInList__bopi0__clmyneutralPion__bc", &temp.npi0);

    printf("%lld entries...\n", tree_upsilon->GetEntries());
    for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
        tree_upsilon->GetEntry(j);

        TotalData.push(temp);
    }
    input_file->Close();
}

void CountNum(double* Nevt, std::string filename, const char* type, const char* MC_version, const char* category, bool smartmode) {
    typedef struct labels {
        int __experiment__;
        int __run__;
        unsigned int __event__;
        int __ncandidates__;
    } Labels;
    std::vector<Labels> label_list;

    std::queue<Data> temp_queue;
    temp_queue.swap(TotalData);
    while (!temp_queue.empty()) {
        Data temp = temp_queue.front();
        temp_queue.pop();

        bool overlap = false;
        for (unsigned int i = 0; i < label_list.size(); i++) {
            if (label_list.at(i).__experiment__ == temp.__experiment__ && label_list.at(i).__run__ == temp.__run__ && label_list.at(i).__event__ == temp.__event__ && label_list.at(i).__ncandidates__ == temp.__ncandidates__) {
                overlap = true;
                break;
            }
        }
        if (overlap == false) {
            // Number of event
            if (smartmode == false) (*Nevt) = (*Nevt) + 1;
            else {
                if (strcmp(type, "SIGNAL") == 0) {
                    printf("SIGNAL type is not available in this code!\n");
                    exit(1);
                }
                else (*Nevt) = (*Nevt) + ObtainWeight(type, MC_version, category, filename) * corrector_Knn.GetCorrectionFactorCancelOutObtainWeight(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, filename, MC_version, false);
            }
            Labels temp_Labels;
            temp_Labels.__experiment__ = temp.__experiment__;
            temp_Labels.__run__ = temp.__run__;
            temp_Labels.__event__ = temp.__event__;
            temp_Labels.__ncandidates__ = temp.__ncandidates__;

            label_list.push_back(temp_Labels);

        }

        TotalData.push(temp);
    }

}

void CustomizedCut() {

    std::queue<Data> temp_queue;
    temp_queue.swap(TotalData);
    while (!temp_queue.empty()) {
        Data temp = temp_queue.front();
        temp_queue.pop();

        if (temp.npi0 > 0.5) continue;

        TotalData.push(temp);
    }

}

void End() {
    std::queue<Data> empty;
    std::swap(TotalData, empty);
}

void ReadPreselection() {

    const char* dirname = "";

    std::vector<std::string> names;

    load_files(dirname, &names);

    double Nevt_1 = 0.0;
    double Nevt_2 = 0.0;

    for (unsigned int i = 0; i < names.size(); i++) {

        Initialize();

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        ReadFiles(input_file);

        CountNum(&Nevt_1, names.at(i), "CHG", "MC15rd", "validation", true);

        CustomizedCut();

        CountNum(&Nevt_2, names.at(i), "CHG", "MC15rd", "validation", true);

        End();

    }

    printf("Nevt_1: %lf\n", Nevt_1);
    printf("Nevt_2: %lf\n", Nevt_2);
    
}
