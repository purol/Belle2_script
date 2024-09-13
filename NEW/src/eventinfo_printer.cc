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
#include "constants.h"
#include "base.h"

/* ====================================== */

void ReadROOT(const char* dirname, const char* included_string, std::vector<int>* experiments_, std::vector<int>* runs_, std::vector<unsigned int>* events_) {

    int __experiment__;
    int __run__;
    unsigned int __event__;
    int __candidate__;
    int __ncandidates__;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;

    double Bsig_M = -1;

    float MVA_BB;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

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

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);

        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_BB);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            // select 0.0 < MXs^reco < 0.6 GeV
            if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) {}
            else continue;

            // select K+
            if ((std::abs(Upsilon_ID) < MyEPSILON) && (std::abs(Bsig_ID) < MyEPSILON)) {}
            else continue;

            // select the last 3 bins
            if (MVA_BB > 0.979) {}
            else continue;

            experiments_->push_back(__experiment__);
            runs_->push_back(__run__);
            events_->push_back(__event__);

        }
    }

}

int main(int argc, char* argv[])
{

    std::vector<int> experiments;
    std::vector<int> runs;
    std::vector<unsigned int> events;

    /* ====================================== */
    // define path for Ntuple
    const char* ROOT_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD_LS_data/SIGNAL_analysis/validation_v008/final_output_root_after_MVA_Application_after_cut";
    /* ====================================== */



    ReadROOT(ROOT_dirname, "root", &experiments, &runs, &events);



    /* ====================================== */
    // file output
    FILE* fp;

    fp = fopen("event_info.csv", "w");
    fprintf(fp, "experiment,run,event\n");
    for (int i = 0; i < experiments.size(); i++) {
        fprintf(fp, "%lf,%lf,%lf\n", experiments.at(i), runs.at(i), events.at(i));
    }
    fclose(fp);
    /* ====================================== */

}
