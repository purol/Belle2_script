#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include <cmath>
#include <float.h>
#include <string>
#include <string.h>
#include "TFile.h"
#include "TTree.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "RooGenericPdf.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "TH1.h"
#include "TF1.h"
#include "RooDataHist.h"
#include <time.h>
#include <random>
#include <algorithm>
#include "TH2.h"
#include "RooRandom.h"
#include "template.h"
#include "constants.h"
#include "ObtainWeight.h"
#include "correctors.h"
#include "base.h"

typedef struct meta {
    int __experiment__;
    int __run__;
    unsigned int __event__;
    int __candidate__;
    int __ncandidates__;
    std::string filename;
} Meta;

void PrintMeta(const char* filename, std::vector<Meta> metas) {
    FILE* fp;

    fp = fopen(filename, "w");
    fprintf(fp, "exp,run,event,filename\n");
    for (int i = 0; i < metas.size(); i++) {
        fprintf(fp, "%d,%d,%u,%s\n", metas.at(i).__experiment__, metas.at(i).__run__, metas.at(i).__event__, metas.at(i).filename.c_str());
    }
    fclose(fp);
}

void GetNominalNevt(const char* dirname, const char* included_string, std::vector<Meta>* metas) { // get nominal PDF with appropriate correction

    int __experiment__;
    int __run__;
    unsigned int __event__;
    int __candidate__;
    int __ncandidates__;

    std::vector<std::string> names;
    load_files(dirname, &names, included_string);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");

        tree_upsilon->SetBranchAddress("__experiment__", &__experiment__);
        tree_upsilon->SetBranchAddress("__run__", &__run__);
        tree_upsilon->SetBranchAddress("__event__", &__event__);
        tree_upsilon->SetBranchAddress("__candidate__", &__candidate__);
        tree_upsilon->SetBranchAddress("__ncandidates__", &__ncandidates__);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            Meta meta = { __experiment__, __run__, __event__, __candidate__, __ncandidates__, names.at(i) };
            metas->push_back(meta);

        }
        input_file->Close();

    }

    return;
}

int main(int argc, char* argv[])
{

    std::vector<Meta> metas_SIGNAL;
    std::vector<Meta> metas_CHG;
    std::vector<Meta> metas_MIX;
    std::vector<Meta> metas_UUBAR;
    std::vector<Meta> metas_DDBAR;
    std::vector<Meta> metas_SSBAR;
    std::vector<Meta> metas_CHARM;


    /* ====================================== */
    // define path for Ntuple
    const char* MC_dirname_SIGNAL = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/SIGNAL_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut";

    const char* MC_dirname_CHG = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/CHG_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut";
    const char* MC_dirname_MIX = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/MIX_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut";
    const char* MC_dirname_UUBAR = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/UUBAR_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut";
    const char* MC_dirname_DDBAR = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/DDBAR_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut";
    const char* MC_dirname_SSBAR = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/SSBAR_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut";
    const char* MC_dirname_CHARM = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/CHARM_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut";
    /* ====================================== */



    /* ====================================== */
    // get nominal Nevt
    GetNominalNevt(MC_dirname_SIGNAL, "B2Knunu", &metas_SIGNAL);
    GetNominalNevt(MC_dirname_SIGNAL, "B2Kstarnunu", &metas_SIGNAL);
    GetNominalNevt(MC_dirname_SIGNAL, "B2Xsnunu", &metas_SIGNAL);
    GetNominalNevt(MC_dirname_SIGNAL, "B02K0nunu", &metas_SIGNAL);
    GetNominalNevt(MC_dirname_SIGNAL, "B02Kstar0nunu", &metas_SIGNAL);
    GetNominalNevt(MC_dirname_SIGNAL, "B02Xsnunu", &metas_SIGNAL);

    GetNominalNevt(MC_dirname_CHG, "root", &metas_CHG);
    GetNominalNevt(MC_dirname_MIX, "root", &metas_MIX);
    GetNominalNevt(MC_dirname_UUBAR, "root", &metas_UUBAR);
    GetNominalNevt(MC_dirname_DDBAR, "root", &metas_DDBAR);
    GetNominalNevt(MC_dirname_SSBAR, "root", &metas_SSBAR);
    GetNominalNevt(MC_dirname_CHARM, "root", &metas_CHARM);
    /* ====================================== */



    /* ====================================== */
    // file output
    PrintMeta("SIGNAL_meta.csv", metas_SIGNAL);

    PrintMeta("CHG_meta.csv", metas_CHG);
    PrintMeta("MIX_meta.csv", metas_MIX);
    PrintMeta("UUBAR_meta.csv", metas_UUBAR);
    PrintMeta("DDBAR_meta.csv", metas_DDBAR);
    PrintMeta("SSBAR_meta.csv", metas_SSBAR);
    PrintMeta("CHARM_meta.csv", metas_CHARM);
    /* ====================================== */

}
