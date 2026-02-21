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
        float BDT_output;

        TFile* input_file = new TFile((argv[1] + std::string("/") + names.at(i)).c_str(), "read");
        //printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* temp_Upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* temp_Bsig = (TTree*)input_file->Get("Bsig");

        temp_Upsilon->SetBranchAddress("__experiment__", &__experiment__);
        temp_Upsilon->SetBranchAddress("__run__", &__run__);
        temp_Upsilon->SetBranchAddress("__event__", &__event__);
        temp_Upsilon->SetBranchAddress("__candidate__", &__candidate__);
        temp_Upsilon->SetBranchAddress("__ncandidates__", &__ncandidates__);

        temp_Bsig->SetBranchAddress("Bsig_M", &MXs);

        temp_Upsilon->SetBranchAddress("MVA_BB", &BDT_output);

        for (unsigned int j = 0; j < temp_Upsilon->GetEntries(); j++) { // Fill
            temp_Upsilon->GetEntry(j);
            temp_Bsig->GetEntry(j);

            if((__experiment__ == 16) && (__run__ == 1096) && (__event__ == 37017320) && (__candidate__ == 0) && (__ncandidates__ == 2)) printf("%lf %f\n", MXs, BDT_output);
            else if ((__experiment__ == 26) && (__run__ == 1527) && (__event__ == 6440078) && (__candidate__ == 0) && (__ncandidates__ == 1)) printf("%lf %f\n", MXs, BDT_output);
            else if ((__experiment__ == 24) && (__run__ == 830) && (__event__ == 36530887) && (__candidate__ == 0) && (__ncandidates__ == 1)) printf("%lf %f\n", MXs, BDT_output);
            else if ((__experiment__ == 26) && (__run__ == 923) && (__event__ == 19393310) && (__candidate__ == 0) && (__ncandidates__ == 1)) printf("%lf %f\n", MXs, BDT_output);
            else if ((__experiment__ == 18) && (__run__ == 346) && (__event__ == 11598646) && (__candidate__ == 0) && (__ncandidates__ == 1)) printf("%lf %f\n", MXs, BDT_output);
        }


        input_file->Close();
        delete input_file;

    }


    return 0;
}
