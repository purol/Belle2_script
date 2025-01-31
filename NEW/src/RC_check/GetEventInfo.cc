// last update: 2023-07-14
// for Belle2 data

#include <stdio.h>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"

#include "constants.h"
#include "base.h"

/*
when you add new variables:
revise # define N_Needed_info ...
revise typedef struct data
revise void Loader::GetData(TFile* input_file)
check MakeShiftDoubleToInt if it affects the code
revise void Loader::PrintRootFile(std::string output_name)
revise void Loader::PrintSeparateRootFile(std::string output_name)
revise void Loader::ConvertIntoSeparateDataFile(std::string output_name, double flag = 0)
*/

/* ====================================== */

int main(int argc, char* argv[]) {
    /*
    * argv[1]: dirname
    */

    std::vector<std::string> names;

    load_files(argv[1], &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        int __experiment__;
        int __run__;
        unsigned int __event__;
        int __candidate__;
        int __ncandidates__;

        double MXs;
        double BDT_output;

        TFile* input_file = new TFile((argv[1] + std::string("/") + names.at(i)).c_str(), "read");
        //printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* temp_tree = (TTree*)input_file->Get("Upsilon");

        temp_tree->SetBranchAddress("__experiment__", &__experiment__);
        temp_tree->SetBranchAddress("__run__", &__run__);
        temp_tree->SetBranchAddress("__event__", &__event__);
        temp_tree->SetBranchAddress("__candidate__", &__candidate__);
        temp_tree->SetBranchAddress("__ncandidates__", &__ncandidates__);

        temp_tree->SetBranchAddress("Bsig_M", &MXs);

        temp_tree->SetBranchAddress("MVA_BB", &BDT_output);

        for (unsigned int j = 0; j < temp_tree->GetEntries(); j++) { // Fill
            temp_tree->GetEntry(j);

            if ((MXs < 0.6) && (BDT_output > 0.979) && (BDT_output < 0.986)) {
                printf("%d %d %u %d %d\n", __experiment__, __run__, __event__, __candidate__, __ncandidates__);
            }
        }


        input_file->Close();
        delete input_file;

    }


    return 0;
}
