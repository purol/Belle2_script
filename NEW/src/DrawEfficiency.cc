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

Corrector corrector;
Corrector_Knn corrector_Knn;
Corrector_Fragmentation corrector_Fragmentation;

# define MCTYPE "MC15rd"

# define Scale_Kplus_here (N_Kplus_nunubar_LS1/1000000.0)
# define Scale_Kplusstar_here (N_Kplusstar_nunubar_LS1/1000000.0)
# define Scale_Xsu_nonresonant_here (N_Xsu_nonresonant_nunubar_LS1/5000000.0)
# define Scale_K0_here (N_K0_nunubar_LS1/1000000.0)
# define Scale_K0star_here (N_K0star_nunubar_LS1/1000000.0)
# define Scale_Xsd_nonresonant_here (N_Xsd_nunubar_LS1/5000000.0)

/* ====================================== */

TH1D* Ncandidate_BKG_after_preselection = new TH1D("Ncandidate_BKG_after_preselection", ";M_{X_{s}}^{reco} [GeV/c^{2}]", 9, 0.45, 3.0);
TH1D* Nevt_BKG_before_FBDT_cut = new TH1D("Nevt_BKG_before_FBDT_cut", ";M_{X_{s}}^{reco} [GeV/c^{2}]", 9, 0.45, 3.0);
TH1D* Nevt_BKG_after_selection = new TH1D("Nevt_BKG_after_selection", ";M_{X_{s}}^{reco} [GeV/c^{2}]", 9, 0.45, 3.0);

TH1D* Ncandidate_SIGNAL_after_preselection = new TH1D("Ncandidate_SIGNAL_after_preselection", ";M_{X_{s}}^{reco} [GeV/c^{2}]", 9, 0.45, 3.0);
TH1D* Nevt_SIGNAL_before_FBDT_cut = new TH1D("Nevt_SIGNAL_before_FBDT_cut", ";M_{X_{s}}^{reco} [GeV/c^{2}]", 9, 0.45, 3.0);
TH1D* Nevt_SIGNAL_after_selection = new TH1D("Nevt_SIGNAL_after_selection", ";M_{X_{s}}^{reco} [GeV/c^{2}]", 9, 0.45, 3.0);

TH1D* Nevt_SIGNAL_initial_true = new TH1D("Nevt_SIGNAL_initial_true", ";M_{X_{s}}^{true} [GeV/c^{2}]", 9, 0.45, 3.0);
TH1D* Nevt_SIGNAL_after_preselection_true = new TH1D("Nevt_SIGNAL_after_preselection_true", ";M_{X_{s}}^{true} [GeV/c^{2}]", 9, 0.45, 3.0);
TH1D* Nevt_SIGNAL_before_FBDT_cut_true = new TH1D("Nevt_SIGNAL_before_FBDT_cut_true", ";M_{X_{s}}^{true} [GeV/c^{2}]", 9, 0.45, 3.0);
TH1D* Nevt_SIGNAL_after_selection_true = new TH1D("Nevt_SIGNAL_after_selection_true", ";M_{X_{s}}^{true} [GeV/c^{2}]", 9, 0.45, 3.0);

TH1D* SIGNAL_pre_and_selection_efficiency_true = new TH1D("SIGNAL_pre_and_selection_efficiency_true", ";M_{X_{s}}^{true} [GeV/c^{2}];Efficiency", 9, 0.45, 3.0);
TH1D* SIGNAL_selection_efficiency = new TH1D("SIGNAL_selection_efficiency", ";M_{X_{s}}^{reco} [GeV/c^{2}];Efficiency", 9, 0.45, 3.0);
TH1D* SIGNAL_FBDT_efficiency = new TH1D("SIGNAL_FBDT_efficiency", ";M_{X_{s}}^{reco} [GeV/c^{2}];FBDT Efficiency", 9, 0.45, 3.0);
TH1D* SIGNAL_FBDT_efficiency_true = new TH1D("SIGNAL_FBDT_efficiency_true", ";M_{X_{s}}^{true} [GeV/c^{2}];FBDT Efficiency", 9, 0.45, 3.0);

TH1D* BKG_selection_efficiency = new TH1D("BKG_selection_efficiency", ";M_{X_{s}}^{reco} [GeV/c^{2}];1 - (rejection rate)", 9, 0.45, 3.0);
TH1D* BKG_FBDT_efficiency = new TH1D("BKG_FBDT_efficiency", ";M_{X_{s}}^{reco} [GeV/c^{2}];FBDT 1 - (rejection rate)", 9, 0.45, 3.0);

TH1D* Nevt_SIGNAL_initial_q2_true = new TH1D("Nevt_SIGNAL_initial_q2_true", ";q^{2} [GeV^{2}/c^{4}]", 10, 0.0, 22.9);
TH1D* Nevt_SIGNAL_after_selection_q2_true = new TH1D("Nevt_SIGNAL_after_selection_q2_true", ";q^{2} [GeV^{2}/c^{4}]", 10, 0.0, 22.9);
TH1D* SIGNAL_efficiency_q2_true = new TH1D("SIGNAL_efficiency_q2_true", ";q^{2} [GeV^{2}/c^{4}];Efficiency", 10, 0.0, 22.9);

void FillBKG(const char* dirname, const char* included_string, TH1D* hist, double weight_var = 1.0) { // get nominal PDF with appropriate correction

    int __experiment__;
    int __run__;
    unsigned int __event__;
    int __candidate__;
    int __ncandidates__;
    double Bsig_M = -1;

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

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            double total_weight = weight_var * corrector_Knn.GetCorrectionFactorCancelOutObtainWeight(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, names.at(i), MCTYPE, false);

            hist->Fill(Bsig_M, total_weight);

            Nevt = Nevt + total_weight;

        }
        input_file->Close();

        printf("%s has %lf events (with correction)\n", dirname, Nevt);

    }

    return;
}

void FillSIGNAL_q2(const char* dirname, const char* included_string, const char* type, TH1D* hist, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get nominal PDF with appropriate correction
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    int __experiment__;
    int __run__;
    unsigned int __event__;
    int __candidate__;
    int __ncandidates__;
    double Bsig_M = -1.0;
    int Decay[N_decay] = { 0 };
    double Mxs_MC = -1.0;
    double invM = -1.0;

    std::vector<std::string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");
        TTree* tree_Xs = (TTree*)input_file->Get("Xs");

        tree_upsilon->SetBranchAddress("__experiment__", &__experiment__);
        tree_upsilon->SetBranchAddress("__run__", &__run__);
        tree_upsilon->SetBranchAddress("__event__", &__event__);
        tree_upsilon->SetBranchAddress("__candidate__", &__candidate__);
        tree_upsilon->SetBranchAddress("__ncandidates__", &__ncandidates__);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKcharge_total__bc", &Decay[0]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &Decay[1]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &Decay[2]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCcomb__bc", &Decay[3]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch1__bc", &Decay[4]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch2__bc", &Decay[5]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch3__bc", &Decay[6]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch4__bc", &Decay[7]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch5__bc", &Decay[8]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch6__bc", &Decay[9]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch7__bc", &Decay[10]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch8__bc", &Decay[11]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch9__bc", &Decay[12]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch10__bc", &Decay[13]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch11__bc", &Decay[14]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch12__bc", &Decay[15]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch13__bc", &Decay[16]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch14__bc", &Decay[17]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch15__bc", &Decay[18]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKneutral_total__bc", &Decay[19]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &Decay[20]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &Decay[21]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCcomb__bc", &Decay[22]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch16__bc", &Decay[23]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch17__bc", &Decay[24]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch18__bc", &Decay[25]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch19__bc", &Decay[26]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch20__bc", &Decay[27]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch21__bc", &Decay[28]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch22__bc", &Decay[29]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch23__bc", &Decay[30]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch24__bc", &Decay[31]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch25__bc", &Decay[32]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch26__bc", &Decay[33]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch27__bc", &Decay[34]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch28__bc", &Decay[35]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch29__bc", &Decay[36]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch30__bc", &Decay[37]);
        if (strcmp(type, "Bplus") == 0) tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_MC);
        else tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_MC);
        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            tree_Xs->GetEntry(j);

            double total_weight = weight_var;

            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            if ((invM * invM >= 0) && (invM * invM <= 30.0)) {}
            else continue;

            hist->Fill(invM * invM, total_weight);

            Nevt = Nevt + total_weight;

        }
        input_file->Close();

        printf("%s has %lf events (with correction)\n", dirname, Nevt);

    }

    return;
}

void FillSIGNAL(const char* dirname, const char* included_string, const char* type, TH1D* hist_reco, TH1D* hist_true, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get nominal PDF with appropriate correction
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    int __experiment__;
    int __run__;
    unsigned int __event__;
    int __candidate__;
    int __ncandidates__;
    double Bsig_M = -1.0;
    int Decay[N_decay] = { 0 };
    double Mxs_MC = -1.0;
    double invM = -1.0;

    std::vector<std::string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");
        TTree* tree_Xs = (TTree*)input_file->Get("Xs");

        tree_upsilon->SetBranchAddress("__experiment__", &__experiment__);
        tree_upsilon->SetBranchAddress("__run__", &__run__);
        tree_upsilon->SetBranchAddress("__event__", &__event__);
        tree_upsilon->SetBranchAddress("__candidate__", &__candidate__);
        tree_upsilon->SetBranchAddress("__ncandidates__", &__ncandidates__);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKcharge_total__bc", &Decay[0]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &Decay[1]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &Decay[2]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCcomb__bc", &Decay[3]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch1__bc", &Decay[4]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch2__bc", &Decay[5]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch3__bc", &Decay[6]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch4__bc", &Decay[7]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch5__bc", &Decay[8]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch6__bc", &Decay[9]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch7__bc", &Decay[10]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch8__bc", &Decay[11]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch9__bc", &Decay[12]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch10__bc", &Decay[13]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch11__bc", &Decay[14]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch12__bc", &Decay[15]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch13__bc", &Decay[16]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch14__bc", &Decay[17]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch15__bc", &Decay[18]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKneutral_total__bc", &Decay[19]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &Decay[20]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &Decay[21]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCcomb__bc", &Decay[22]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch16__bc", &Decay[23]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch17__bc", &Decay[24]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch18__bc", &Decay[25]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch19__bc", &Decay[26]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch20__bc", &Decay[27]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch21__bc", &Decay[28]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch22__bc", &Decay[29]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch23__bc", &Decay[30]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch24__bc", &Decay[31]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch25__bc", &Decay[32]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch26__bc", &Decay[33]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch27__bc", &Decay[34]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch28__bc", &Decay[35]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch29__bc", &Decay[36]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch30__bc", &Decay[37]);
        if(strcmp(type, "Bplus") == 0) tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_MC);
        else tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_MC);
        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            tree_Xs->GetEntry(j);

            double total_weight = weight_var;

            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            hist_reco->Fill(Bsig_M, total_weight);
            hist_true->Fill(Mxs_MC, total_weight);

            Nevt = Nevt + total_weight;

        }
        input_file->Close();

        printf("%s has %lf events (with correction)\n", dirname, Nevt);

    }

    return;
}

void ReadDecayInfo(const char* dirname, const char* included_string, const char* type, TH1D* hist_true, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get nominal PDF with appropriate correction
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    int Decay[N_decay] = { 0 };
    double Mxs_MC = -1.0;
    double invM = -1.0;

    std::vector<std::string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_Xs = (TTree*)input_file->Get("Xs");

        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKcharge_total__bc", &Decay[0]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &Decay[1]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &Decay[2]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clcomb__bc", &Decay[3]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch1__bc", &Decay[4]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch2__bc", &Decay[5]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch3__bc", &Decay[6]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch4__bc", &Decay[7]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch5__bc", &Decay[8]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch6__bc", &Decay[9]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch7__bc", &Decay[10]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch8__bc", &Decay[11]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch9__bc", &Decay[12]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch10__bc", &Decay[13]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch11__bc", &Decay[14]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch12__bc", &Decay[15]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch13__bc", &Decay[16]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch14__bc", &Decay[17]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch15__bc", &Decay[18]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKneutral_total__bc", &Decay[19]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &Decay[20]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &Decay[21]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clcomb__bc", &Decay[22]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch16__bc", &Decay[23]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch17__bc", &Decay[24]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch18__bc", &Decay[25]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch19__bc", &Decay[26]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch20__bc", &Decay[27]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch21__bc", &Decay[28]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch22__bc", &Decay[29]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch23__bc", &Decay[30]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch24__bc", &Decay[31]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch25__bc", &Decay[32]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch26__bc", &Decay[33]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch27__bc", &Decay[34]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch28__bc", &Decay[35]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch29__bc", &Decay[36]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch30__bc", &Decay[37]);
        if (strcmp(type, "Bplus") == 0) tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_MC);
        else tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_MC);
        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

        printf("%lld entries...\n", tree_Xs->GetEntries());
        for (unsigned int j = 0; j < tree_Xs->GetEntries(); j++) { // Fill
            tree_Xs->GetEntry(j);

            double total_weight = weight_var;

            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            if (std::isnan(total_weight)) {
                total_weight = weight_var;
            }

            hist_true->Fill(Mxs_MC, total_weight);

            Nevt = Nevt + total_weight;

        }
        input_file->Close();

        printf("%s has %lf events (with correction)\n", dirname, Nevt);

    }

    return;
}

void ReadDecayInfo_q2(const char* dirname, const char* included_string, const char* type, TH1D* hist, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get nominal PDF with appropriate correction
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    int Decay[N_decay] = { 0 };
    double Mxs_MC = -1.0;
    double invM = -1.0;

    std::vector<std::string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_Xs = (TTree*)input_file->Get("Xs");

        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKcharge_total__bc", &Decay[0]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &Decay[1]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &Decay[2]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clcomb__bc", &Decay[3]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch1__bc", &Decay[4]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch2__bc", &Decay[5]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch3__bc", &Decay[6]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch4__bc", &Decay[7]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch5__bc", &Decay[8]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch6__bc", &Decay[9]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch7__bc", &Decay[10]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch8__bc", &Decay[11]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch9__bc", &Decay[12]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch10__bc", &Decay[13]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch11__bc", &Decay[14]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch12__bc", &Decay[15]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch13__bc", &Decay[16]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch14__bc", &Decay[17]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch15__bc", &Decay[18]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKneutral_total__bc", &Decay[19]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &Decay[20]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &Decay[21]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clcomb__bc", &Decay[22]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch16__bc", &Decay[23]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch17__bc", &Decay[24]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch18__bc", &Decay[25]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch19__bc", &Decay[26]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch20__bc", &Decay[27]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch21__bc", &Decay[28]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch22__bc", &Decay[29]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch23__bc", &Decay[30]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch24__bc", &Decay[31]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch25__bc", &Decay[32]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch26__bc", &Decay[33]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch27__bc", &Decay[34]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch28__bc", &Decay[35]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch29__bc", &Decay[36]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch30__bc", &Decay[37]);
        if (strcmp(type, "Bplus") == 0) tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_MC);
        else tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_MC);
        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

        printf("%lld entries...\n", tree_Xs->GetEntries());
        for (unsigned int j = 0; j < tree_Xs->GetEntries(); j++) { // Fill
            tree_Xs->GetEntry(j);

            double total_weight = weight_var;

            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            if (std::isnan(total_weight)) {
                total_weight = weight_var;
            }

            if ((invM * invM >= 0) && (invM * invM <= 30.0)) {}
            else continue;

            hist->Fill(invM * invM, total_weight);

            Nevt = Nevt + total_weight;

        }
        input_file->Close();

        printf("%s has %lf events (with correction)\n", dirname, Nevt);

    }

    return;
}

int main(int argc, char* argv[]) {
    // function to calculate efficiency as a function of MXs

    const char* dirname_CHG_after_preselection = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD/CHG_analysis/validation_v008/before_Mbc_cut";
    const char* dirname_MIX_after_preselection = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD/MIX_analysis/validation_v008/before_Mbc_cut";
    const char* dirname_UUBAR_after_preselection = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD/UUBAR_analysis/validation_v008/before_Mbc_cut";
    const char* dirname_DDBAR_after_preselection = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD/DDBAR_analysis/validation_v008/before_Mbc_cut";
    const char* dirname_SSBAR_after_preselection = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD/SSBAR_analysis/validation_v008/before_Mbc_cut";
    const char* dirname_CHARM_after_preselection = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD/CHARM_analysis/validation_v008/before_Mbc_cut";

    const char* dirname_CHG_after_before_FBDT_cut = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD/CHG_analysis/validation_v008/final_output";
    const char* dirname_MIX_after_before_FBDT_cut = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD/MIX_analysis/validation_v008/final_output";
    const char* dirname_UUBAR_after_before_FBDT_cut = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD/UUBAR_analysis/validation_v008/final_output";
    const char* dirname_DDBAR_after_before_FBDT_cut = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD/DDBAR_analysis/validation_v008/final_output";
    const char* dirname_SSBAR_after_before_FBDT_cut = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD/SSBAR_analysis/validation_v008/final_output";
    const char* dirname_CHARM_after_before_FBDT_cut = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD/CHARM_analysis/validation_v008/final_output";

    const char* dirname_CHG_after_selection = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD/CHG_analysis/validation_v008/final_output_root_after_MVA_Application_after_cut/Merge";
    const char* dirname_MIX_after_selection = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD/MIX_analysis/validation_v008/final_output_root_after_MVA_Application_after_cut/Merge";
    const char* dirname_UUBAR_after_selection = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD/UUBAR_analysis/validation_v008/final_output_root_after_MVA_Application_after_cut/Merge";
    const char* dirname_DDBAR_after_selection = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD/DDBAR_analysis/validation_v008/final_output_root_after_MVA_Application_after_cut/Merge";
    const char* dirname_SSBAR_after_selection = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD/SSBAR_analysis/validation_v008/final_output_root_after_MVA_Application_after_cut/Merge";
    const char* dirname_CHARM_after_selection = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD/CHARM_analysis/validation_v008/final_output_root_after_MVA_Application_after_cut/Merge";

    const char* dirname_SIGNAL_initial = "/home/belle2/junewoo/storage_ghi/20220929_SIGNAL_decayInfo_again/small";
    const char* dirname_SIGNAL_after_preselection = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD/SIGNAL_analysis/validation_v008/before_Mbc_cut";
    const char* dirname_SIGNAL_after_before_FBDT_cut = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD/SIGNAL_analysis/validation_v008/final_output";
    const char* dirname_SIGNAL_after_selection = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD/SIGNAL_analysis/validation_v008/final_output_root_after_MVA_Application_after_cut/Merge";

    FillBKG(dirname_CHG_after_preselection, "root", Ncandidate_BKG_after_preselection, ObtainWeight("CHG", MCTYPE, "validation", "CHG"));
    FillBKG(dirname_MIX_after_preselection, "root", Ncandidate_BKG_after_preselection, ObtainWeight("MIX", MCTYPE, "validation", "MIX"));
    FillBKG(dirname_UUBAR_after_preselection, "root", Ncandidate_BKG_after_preselection, ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"));
    FillBKG(dirname_DDBAR_after_preselection, "root", Ncandidate_BKG_after_preselection, ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"));
    FillBKG(dirname_SSBAR_after_preselection, "root", Ncandidate_BKG_after_preselection, ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"));
    FillBKG(dirname_CHARM_after_preselection, "root", Ncandidate_BKG_after_preselection, ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"));

    FillBKG(dirname_CHG_after_before_FBDT_cut, "root", Nevt_BKG_before_FBDT_cut, ObtainWeight("CHG", MCTYPE, "validation", "CHG"));
    FillBKG(dirname_MIX_after_before_FBDT_cut, "root", Nevt_BKG_before_FBDT_cut, ObtainWeight("MIX", MCTYPE, "validation", "MIX"));
    FillBKG(dirname_UUBAR_after_before_FBDT_cut, "root", Nevt_BKG_before_FBDT_cut, ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"));
    FillBKG(dirname_DDBAR_after_before_FBDT_cut, "root", Nevt_BKG_before_FBDT_cut, ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"));
    FillBKG(dirname_SSBAR_after_before_FBDT_cut, "root", Nevt_BKG_before_FBDT_cut, ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"));
    FillBKG(dirname_CHARM_after_before_FBDT_cut, "root", Nevt_BKG_before_FBDT_cut, ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"));

    FillBKG(dirname_CHG_after_selection, "root", Nevt_BKG_after_selection, ObtainWeight("CHG", MCTYPE, "validation", "CHG"));
    FillBKG(dirname_MIX_after_selection, "root", Nevt_BKG_after_selection, ObtainWeight("MIX", MCTYPE, "validation", "MIX"));
    FillBKG(dirname_UUBAR_after_selection, "root", Nevt_BKG_after_selection, ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"));
    FillBKG(dirname_DDBAR_after_selection, "root", Nevt_BKG_after_selection, ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"));
    FillBKG(dirname_SSBAR_after_selection, "root", Nevt_BKG_after_selection, ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"));
    FillBKG(dirname_CHARM_after_selection, "root", Nevt_BKG_after_selection, ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"));

    FillSIGNAL(dirname_SIGNAL_after_preselection, "B2Knunu", "Bplus", Ncandidate_SIGNAL_after_preselection, Nevt_SIGNAL_after_preselection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    FillSIGNAL(dirname_SIGNAL_after_preselection, "B2Kstarnunu", "Bplus", Ncandidate_SIGNAL_after_preselection, Nevt_SIGNAL_after_preselection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    FillSIGNAL(dirname_SIGNAL_after_preselection, "B2Xsnunu", "Bplus", Ncandidate_SIGNAL_after_preselection, Nevt_SIGNAL_after_preselection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    FillSIGNAL(dirname_SIGNAL_after_preselection, "B02K0nunu", "Bzero", Ncandidate_SIGNAL_after_preselection, Nevt_SIGNAL_after_preselection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    FillSIGNAL(dirname_SIGNAL_after_preselection, "B02Kstar0nunu", "Bzero", Ncandidate_SIGNAL_after_preselection, Nevt_SIGNAL_after_preselection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    FillSIGNAL(dirname_SIGNAL_after_preselection, "B02Xsnunu", "Bzero", Ncandidate_SIGNAL_after_preselection, Nevt_SIGNAL_after_preselection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    FillSIGNAL(dirname_SIGNAL_after_before_FBDT_cut, "B2Knunu", "Bplus", Nevt_SIGNAL_before_FBDT_cut, Nevt_SIGNAL_before_FBDT_cut_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    FillSIGNAL(dirname_SIGNAL_after_before_FBDT_cut, "B2Kstarnunu", "Bplus", Nevt_SIGNAL_before_FBDT_cut, Nevt_SIGNAL_before_FBDT_cut_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    FillSIGNAL(dirname_SIGNAL_after_before_FBDT_cut, "B2Xsnunu", "Bplus", Nevt_SIGNAL_before_FBDT_cut, Nevt_SIGNAL_before_FBDT_cut_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    FillSIGNAL(dirname_SIGNAL_after_before_FBDT_cut, "B02K0nunu", "Bzero", Nevt_SIGNAL_before_FBDT_cut, Nevt_SIGNAL_before_FBDT_cut_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    FillSIGNAL(dirname_SIGNAL_after_before_FBDT_cut, "B02Kstar0nunu", "Bzero", Nevt_SIGNAL_before_FBDT_cut, Nevt_SIGNAL_before_FBDT_cut_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    FillSIGNAL(dirname_SIGNAL_after_before_FBDT_cut, "B02Xsnunu", "Bzero", Nevt_SIGNAL_before_FBDT_cut, Nevt_SIGNAL_before_FBDT_cut_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    FillSIGNAL(dirname_SIGNAL_after_selection, "B2Knunu", "Bplus", Nevt_SIGNAL_after_selection, Nevt_SIGNAL_after_selection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    FillSIGNAL(dirname_SIGNAL_after_selection, "B2Kstarnunu", "Bplus", Nevt_SIGNAL_after_selection, Nevt_SIGNAL_after_selection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    FillSIGNAL(dirname_SIGNAL_after_selection, "B2Xsnunu", "Bplus", Nevt_SIGNAL_after_selection, Nevt_SIGNAL_after_selection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    FillSIGNAL(dirname_SIGNAL_after_selection, "B02K0nunu", "Bzero", Nevt_SIGNAL_after_selection, Nevt_SIGNAL_after_selection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    FillSIGNAL(dirname_SIGNAL_after_selection, "B02Kstar0nunu", "Bzero", Nevt_SIGNAL_after_selection, Nevt_SIGNAL_after_selection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    FillSIGNAL(dirname_SIGNAL_after_selection, "B02Xsnunu", "Bzero", Nevt_SIGNAL_after_selection, Nevt_SIGNAL_after_selection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    FillSIGNAL_q2(dirname_SIGNAL_after_selection, "B2Knunu", "Bplus", Nevt_SIGNAL_after_selection_q2_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    FillSIGNAL_q2(dirname_SIGNAL_after_selection, "B2Kstarnunu", "Bplus", Nevt_SIGNAL_after_selection_q2_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    FillSIGNAL_q2(dirname_SIGNAL_after_selection, "B2Xsnunu", "Bplus", Nevt_SIGNAL_after_selection_q2_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    FillSIGNAL_q2(dirname_SIGNAL_after_selection, "B02K0nunu", "Bzero", Nevt_SIGNAL_after_selection_q2_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    FillSIGNAL_q2(dirname_SIGNAL_after_selection, "B02Kstar0nunu", "Bzero", Nevt_SIGNAL_after_selection_q2_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    FillSIGNAL_q2(dirname_SIGNAL_after_selection, "B02Xsnunu", "Bzero", Nevt_SIGNAL_after_selection_q2_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    ReadDecayInfo(dirname_SIGNAL_initial, "B2Knunu", "Bplus", Nevt_SIGNAL_initial_true, Scale_Kplus_here, "B2Knunu");
    ReadDecayInfo(dirname_SIGNAL_initial, "B2Kstarnunu", "Bplus", Nevt_SIGNAL_initial_true, Scale_Kplusstar_here, "otherwise");
    ReadDecayInfo(dirname_SIGNAL_initial, "B2Xsnunu", "Bplus", Nevt_SIGNAL_initial_true, Scale_Xsu_nonresonant_here, "B2Xsnunu");
    ReadDecayInfo(dirname_SIGNAL_initial, "B02K0nunu", "Bzero", Nevt_SIGNAL_initial_true, Scale_K0_here, "B02K0nunu");
    ReadDecayInfo(dirname_SIGNAL_initial, "B02Kstar0nunu", "Bzero", Nevt_SIGNAL_initial_true, Scale_K0star_here, "otherwise");
    ReadDecayInfo(dirname_SIGNAL_initial, "B02Xsnunu", "Bzero", Nevt_SIGNAL_initial_true, Scale_Xsd_nonresonant_here, "B02Xsnunu");

    ReadDecayInfo_q2(dirname_SIGNAL_initial, "B2Knunu", "Bplus", Nevt_SIGNAL_initial_q2_true, Scale_Kplus_here, "B2Knunu");
    ReadDecayInfo_q2(dirname_SIGNAL_initial, "B2Kstarnunu", "Bplus", Nevt_SIGNAL_initial_q2_true, Scale_Kplusstar_here, "otherwise");
    ReadDecayInfo_q2(dirname_SIGNAL_initial, "B2Xsnunu", "Bplus", Nevt_SIGNAL_initial_q2_true, Scale_Xsu_nonresonant_here, "B2Xsnunu");
    ReadDecayInfo_q2(dirname_SIGNAL_initial, "B02K0nunu", "Bzero", Nevt_SIGNAL_initial_q2_true, Scale_K0_here, "B02K0nunu");
    ReadDecayInfo_q2(dirname_SIGNAL_initial, "B02Kstar0nunu", "Bzero", Nevt_SIGNAL_initial_q2_true, Scale_K0star_here, "otherwise");
    ReadDecayInfo_q2(dirname_SIGNAL_initial, "B02Xsnunu", "Bzero", Nevt_SIGNAL_initial_q2_true, Scale_Xsd_nonresonant_here, "B02Xsnunu");

    SIGNAL_pre_and_selection_efficiency_true->Divide(Nevt_SIGNAL_after_selection_true, Nevt_SIGNAL_initial_true, 1.0, 1.0, "B");
    SIGNAL_selection_efficiency->Divide(Nevt_SIGNAL_after_selection, Ncandidate_SIGNAL_after_preselection, 1.0, 1.0, "B");
    SIGNAL_FBDT_efficiency->Divide(Nevt_SIGNAL_after_selection, Nevt_SIGNAL_before_FBDT_cut, 1.0, 1.0, "B");

    BKG_selection_efficiency->Divide(Nevt_BKG_after_selection, Ncandidate_BKG_after_preselection, 1.0, 1.0, "B");
    BKG_FBDT_efficiency->Divide(Nevt_BKG_after_selection, Nevt_BKG_before_FBDT_cut, 1.0, 1.0, "B");

    SIGNAL_FBDT_efficiency_true->Divide(Nevt_SIGNAL_after_selection_true, Nevt_SIGNAL_before_FBDT_cut_true, 1.0, 1.0, "B");
    SIGNAL_FBDT_efficiency_true->SetStats(false);

    SIGNAL_efficiency_q2_true->Divide(Nevt_SIGNAL_after_selection_q2_true, Nevt_SIGNAL_initial_q2_true, 1.0, 1.0, "B");
    SIGNAL_efficiency_q2_true->SetStats(false);

    // define TCanvas
    TCanvas* c_temp = new TCanvas("c", "", 1200, 1000);

    // set draw option
    SIGNAL_pre_and_selection_efficiency_true->SetMarkerSize(0.0);
    SIGNAL_pre_and_selection_efficiency_true->SetLineWidth(2);
    SIGNAL_pre_and_selection_efficiency_true->SetMarkerColor(4);
    SIGNAL_pre_and_selection_efficiency_true->SetMarkerStyle(21);
    SIGNAL_pre_and_selection_efficiency_true->SetStats(false);

    SIGNAL_FBDT_efficiency->SetMarkerSize(0.0);
    SIGNAL_FBDT_efficiency->SetLineWidth(2);
    SIGNAL_FBDT_efficiency->SetMarkerColor(4);
    SIGNAL_FBDT_efficiency->SetMarkerStyle(21);
    SIGNAL_FBDT_efficiency->SetStats(false);

    BKG_FBDT_efficiency->SetMarkerSize(0.0);
    BKG_FBDT_efficiency->SetLineWidth(2);
    BKG_FBDT_efficiency->SetMarkerColor(4);
    BKG_FBDT_efficiency->SetMarkerStyle(21);
    BKG_FBDT_efficiency->SetStats(false);

    BKG_selection_efficiency->SetMarkerSize(0.0);
    BKG_selection_efficiency->SetLineWidth(2);
    BKG_selection_efficiency->SetMarkerColor(4);
    BKG_selection_efficiency->SetMarkerStyle(21);
    BKG_selection_efficiency->SetStats(false);

    TPaveText* pt;

    // draw SIGNAL preselection + selection efficiency true
    SIGNAL_pre_and_selection_efficiency_true->Draw("P");
    pt = new TPaveText(0.135, 0.88, 0.5, 1.0, "NDC NB"); pt->SetFillStyle(0); pt->SetLineWidth(0); pt->AddText("efficiency for preselection + selection"); pt->Draw();
    c_temp->SaveAs("Efficiency.png");

    // draw SIGNAL FBDT efficiency reco
    SIGNAL_FBDT_efficiency->Draw("P");
    pt = new TPaveText(0.135, 0.88, 0.5, 1.0, "NDC NB"); pt->SetFillStyle(0); pt->SetLineWidth(0); pt->AddText("efficiency for FBDT"); pt->Draw();
    c_temp->SaveAs("Signal_FBDT_reco.png");

    // draw SIGNAL FBDT efficiency true
    SIGNAL_FBDT_efficiency_true->Draw("P");
    pt = new TPaveText(0.135, 0.88, 0.5, 1.0, "NDC NB"); pt->SetFillStyle(0); pt->SetLineWidth(0); pt->AddText("efficiency for FBDT"); pt->Draw();
    c_temp->SaveAs("FBDT_Efficiency.png");

    // draw SIGNAL efficiency reco
    SIGNAL_selection_efficiency->Draw("P");
    pt = new TPaveText(0.135, 0.88, 0.5, 1.0, "NDC NB"); pt->SetFillStyle(0); pt->SetLineWidth(0); pt->AddText("efficiency for selection"); pt->Draw();
    c_temp->SaveAs("Signal_eff_reco.png");

    // draw BKG FBDT efficiency
    BKG_FBDT_efficiency->Draw("P");
    pt = new TPaveText(0.135, 0.88, 0.5, 1.0, "NDC NB"); pt->SetFillStyle(0); pt->SetLineWidth(0); pt->AddText("(1 - rejection rate) for FBDT"); pt->Draw();
    c_temp->SaveAs("Background_FBDT.png");

    // draw BKG 1 - (rejection)
    BKG_selection_efficiency->Draw("P");
    pt = new TPaveText(0.135, 0.88, 0.5, 1.0, "NDC NB"); pt->SetFillStyle(0); pt->SetLineWidth(0); pt->AddText("(1 - rejection rate) for selection"); pt->Draw();
    c_temp->SaveAs("Background_eff_reco.png");

    // draw efficiency as a function of q2
    SIGNAL_efficiency_q2_true->Draw("P");
    pt = new TPaveText(0.135, 0.88, 0.5, 1.0, "NDC NB"); pt->SetFillStyle(0); pt->SetLineWidth(0); pt->AddText("efficiency"); pt->Draw();
    c_temp->SaveAs("Efficiency_q2.png");

    return 0;
}
