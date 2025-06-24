// last update: 2023-07-14
// for Belle2 data

#include <stdio.h>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

#include "constants.h"
#include "base.h"
#include "ObtainWeight.h"
#include "correctors.h"
#include "TCanvas.h"
#include "TPaveText.h"

typedef struct info {
    double Bsig_M = -1;
    double FBDToutput = -1;
    double Mbc_tag = -1;
} Info;

void Fill(const char* dirname, const char* included_string, std::vector<Info>* data) { // get nominal PDF with appropriate correction

    int __experiment__;
    int __run__;
    unsigned int __event__;
    int __candidate__;
    int __ncandidates__;
    double Bsig_M = -1;

    float FBDToutput = -1;
    double Mbc_tag = -1;

    std::vector<std::string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("__experiment__", &__experiment__);
        tree_upsilon->SetBranchAddress("__run__", &__run__);
        tree_upsilon->SetBranchAddress("__event__", &__event__);
        tree_upsilon->SetBranchAddress("__candidate__", &__candidate__);
        tree_upsilon->SetBranchAddress("__ncandidates__", &__ncandidates__);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);
        tree_Btag->SetBranchAddress("Btag_Mbc", &Mbc_tag);
        tree_upsilon->SetBranchAddress("MVA_BB", &FBDToutput);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            data->push_back({ Bsig_M, FBDToutput, Mbc_tag });

        }
        input_file->Close();

        printf("%s has %lf events (with correction)\n", dirname, Nevt);

    }

    return;
}

int main(int argc, char* argv[]) {
    // print CSV file to check the correlation between Mbc and FBDT output

    std::vector<Info> data;

    const char* dirname_CHG = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD_LS_MC_ABCD/CHG_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut";
    const char* dirname_MIX = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD_LS_MC_ABCD/MIX_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut";
    const char* dirname_UUBAR = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD_LS_MC_ABCD/UUBAR_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut";
    const char* dirname_DDBAR = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD_LS_MC_ABCD/DDBAR_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut";
    const char* dirname_SSBAR = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD_LS_MC_ABCD/SSBAR_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut";
    const char* dirname_CHARM = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD_LS_MC_ABCD/CHARM_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut";


    Fill(dirname_CHG, "root", &data);
    Fill(dirname_MIX, "root", &data);
    Fill(dirname_UUBAR, "root", &data);
    Fill(dirname_DDBAR, "root", &data);
    Fill(dirname_SSBAR, "root", &data);
    Fill(dirname_CHARM, "root", &data);

    FILE* fp = fopen("ABCD_method.csv", "w");

    fprintf(fp, "MB,Mbc,FBDT\n");
    for (int i = 0; i < data.size(); i++) fprintf(fp, "%lf,%lf,%lf\n", data.at(i).Bsig_M, data.at(i).Mbc_tag, data.at(i).FBDToutput);

    fclose(fp);

    return 0;
}
