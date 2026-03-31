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
Corrector_FEI corrector_FEI;
Corrector_PID corrector_PID;
Corrector_pi0 corrector_pi0;
Corrector_FakePID corrector_FakePID;
Corrector_Knn corrector_Knn;
Corrector_Xsnn corrector_Xsnn;
Corrector_Multiplicity corrector_Multiplicity;
Corrector_KpKLKL corrector_KpKLKL;
Corrector_KSKLKL corrector_KSKLKL;
Corrector_phiKL corrector_phiKL;
Corrector_KstarKLKL corrector_KstarKLKL;
Corrector_XsKLKL corrector_XsKLKL;
Corrector_BtoDtoXKL corrector_BtoDtoXKL;
Corrector_KS0 corrector_KS0;
Corrector_Fragmentation corrector_Fragmentation;

# define MCTYPE "MC15rd"

# define Scale_Kplus_here (N_Kplus_nunubar_LS1/1000000.0)
# define Scale_Kplusstar_here (N_Kplusstar_nunubar_LS1/1000000.0)
# define Scale_Xsu_nonresonant_here (N_Xsu_nonresonant_nunubar_LS1/5000000.0)
# define Scale_K0_here (N_K0_nunubar_LS1/1000000.0)
# define Scale_K0star_here (N_K0star_nunubar_LS1/1000000.0)
# define Scale_Xsd_nonresonant_here (N_Xsd_nunubar_LS1/5000000.0)

/* ====================================== */

std::vector<double> boundary = { 0.45, 0.75, 1.05, 1.35, 1.65, 1.95, 2.15, 2.45, 2.75, 3.05, 5.28 };

TH1D* Ncandidate_BKG_after_preselection = new TH1D("Ncandidate_BKG_after_preselection", ";M_{X_{s}}^{reco} [GeV/c^{2}]", 9, 0.45, 3.0);
TH1D* Nevt_BKG_before_FBDT_cut = new TH1D("Nevt_BKG_before_FBDT_cut", ";M_{X_{s}}^{reco} [GeV/c^{2}]", 9, 0.45, 3.0);
TH1D* Nevt_BKG_after_selection = new TH1D("Nevt_BKG_after_selection", ";M_{X_{s}}^{reco} [GeV/c^{2}]", 9, 0.45, 3.0);

TH1D* Ncandidate_SIGNAL_after_preselection = new TH1D("Ncandidate_SIGNAL_after_preselection", ";M_{X_{s}}^{reco} [GeV/c^{2}]", 9, 0.45, 3.0);
TH1D* Nevt_SIGNAL_before_FBDT_cut = new TH1D("Nevt_SIGNAL_before_FBDT_cut", ";M_{X_{s}}^{reco} [GeV/c^{2}]", 9, 0.45, 3.0);
TH1D* Nevt_SIGNAL_after_selection = new TH1D("Nevt_SIGNAL_after_selection", ";M_{X_{s}}^{reco} [GeV/c^{2}]", 9, 0.45, 3.0);

TH1D* Nevt_SIGNAL_initial_true = new TH1D("Nevt_SIGNAL_initial_true", ";M_{X_{s}}^{true} [GeV/c^{2}]", 10, 0, 10);
TH1D* Nevt_SIGNAL_after_preselection_true = new TH1D("Nevt_SIGNAL_after_preselection_true", ";M_{X_{s}}^{true} [GeV/c^{2}]", 9, 0.45, 3.0);
TH1D* Nevt_SIGNAL_before_FBDT_cut_true = new TH1D("Nevt_SIGNAL_before_FBDT_cut_true", ";M_{X_{s}}^{true} [GeV/c^{2}]", 9, 0.45, 3.0);
TH1D* Nevt_SIGNAL_after_selection_true = new TH1D("Nevt_SIGNAL_after_selection_true", ";M_{X_{s}}^{true} [GeV/c^{2}]", 10, 0, 10);

TH1D* SIGNAL_pre_and_selection_efficiency_true = new TH1D("SIGNAL_pre_and_selection_efficiency_true", ";M_{X_{s}}^{true} [GeV/c^{2}];Efficiency", 10, 0, 10);
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

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_pi0_syst_MC15ri[N_pi0_syst_MC15ri] = { 0.0 };
    double temp_N_pi0_syst_MC15rd[N_pi0_syst_MC15rd] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    double KS0_flight_distance = 0;
    double KS0_costheta = 0;
    double KS0_p = 0;
    double KS0_D1p = 0;
    double KS0_D2p = 0;

    double Bsig_M = -1.0;
    int Decay[N_decay] = { 0 };
    double Mxs_MC = -1.0;
    double invM = -1.0;

    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double invM_Xnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;
    double N_Xplusnn = 0;
    double N_Xzeronn = 0;

    double Ngamma_v200 = -1;

    double s13_KpKLKL = -1;
    double s23_KpKLKL = -1;
    double nB2KpKLKL_all_KpKLKL = -1;
    double nB2KpKLKL_NR_KpKLKL = -1;

    double s13_KSKLKL = -1;
    double s23_KSKLKL = -1;
    double s12_KSKLKL = -1;
    double nB2KSKLKL_all_KSKLKL = -1;
    double nB2KSKLKL_NR_KSKLKL = -1;
    double nB02KLphi2KSKL = -1;

    double nKL_XKLKL = -1;
    double XKLKL_E_1st = -1;
    double XKLKL_px_1st = -1;
    double XKLKL_py_1st = -1;
    double XKLKL_pz_1st = -1;
    double XKLKL_E_2nd = -1;
    double XKLKL_px_2nd = -1;
    double XKLKL_py_2nd = -1;
    double XKLKL_pz_2nd = -1;
    double nB2KstarKLKL = -1;
    double nB02KstarKLKL = -1;

    double nDptoXKL = -1;
    double nD0toXKL = -1;

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

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);

        tree_upsilon->SetBranchAddress("__experiment__", &__experiment__);
        tree_upsilon->SetBranchAddress("__run__", &__run__);
        tree_upsilon->SetBranchAddress("__event__", &__event__);
        tree_upsilon->SetBranchAddress("__candidate__", &__candidate__);
        tree_upsilon->SetBranchAddress("__ncandidates__", &__ncandidates__);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);

        for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKtruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[0][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKmisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[1][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npitruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[2][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npimisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[3][i_PID]);
        }
        for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npi0bin" + std::to_string(i_pi0)).c_str(), &temp_N_bin_pi0[i_pi0]);
        for (int i_pi0 = 0; i_pi0 < N_pi0_syst_MC15ri; i_pi0++) tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npi0MC15ribin" + std::to_string(i_pi0)).c_str(), &temp_N_pi0_syst_MC15ri[i_pi0]);
        for (int i_pi0 = 0; i_pi0 < N_pi0_syst_MC15rd; i_pi0++) tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npi0MC15rdbin" + std::to_string(i_pi0)).c_str(), &temp_N_pi0_syst_MC15rd[i_pi0]);
        for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeEbin_n" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeE[0][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeEbin_p" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeE[1][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeEbin_n" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeE[2][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeEbin_p" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeE[3][i_fake]);
        }
        for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeMUbin_n" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeMU[0][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeMUbin_p" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeMU[1][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeMUbin_n" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeMU[2][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeMUbin_p" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeMU[3][i_fake]);
        }
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_3D_distance", &KS0_flight_distance);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_costheta", &KS0_costheta);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_p", &KS0_p);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_D1p", &KS0_D1p);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_D2p", &KS0_D2p);

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


        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clXnn__bc", &N_Xplusnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clXnn__bc", &N_Xzeronn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clXnn__bc", &invM_Xnn);

        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s13_KpKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s23_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_all__bc", &nB2KpKLKL_all_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_NR__bc", &nB2KpKLKL_NR_KpKLKL);

        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s13_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &s23_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s12_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_all__bc", &nB2KSKLKL_all_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_NR__bc", &nB2KSKLKL_NR_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_phi__bc", &nB02KLphi2KSKL);

        tree_upsilon->SetBranchAddress("nParticlesInList__boK_L0__clXKLKL__bc", &nKL_XKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__spE__bc", &XKLKL_E_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__sppx__bc", &XKLKL_px_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__sppy__bc", &XKLKL_py_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__sppz__bc", &XKLKL_pz_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__spE__bc", &XKLKL_E_2nd);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppx__bc", &XKLKL_px_2nd);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppy__bc", &XKLKL_py_2nd);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppz__bc", &XKLKL_pz_2nd);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarKLKL__bc", &nB2KstarKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstarKLKL__bc", &nB02KstarKLKL);

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            tree_Xs->GetEntry(j);

            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            double Correction_KS0 = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE), temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE), temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            if (std::string(MCTYPE) == "MC15ri") for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            else if (std::string(MCTYPE) == "MC15rd") for (int i_pi0 = 0; i_pi0 < N_pi0_syst_MC15rd; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_pi0_syst_MC15rd[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(0, i_fake, MCTYPE), temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(1, i_fake, MCTYPE), temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(2, i_fake, MCTYPE), temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(3, i_fake, MCTYPE), temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(0, i_fake, MCTYPE), temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(1, i_fake, MCTYPE), temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(2, i_fake, MCTYPE), temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(3, i_fake, MCTYPE), temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }
            if (std::string(MCTYPE) == "MC15ri") Correction_KS0 = Correction_KS0 * corrector_KS0.GetCorrectionFactor(KS0_flight_distance, MCTYPE);
            else if (std::string(MCTYPE) == "MC15rd") Correction_KS0 = Correction_KS0 * corrector_KS0.GetCorrectionFactor(KS0_p, KS0_costheta, KS0_flight_distance, KS0_D1p, KS0_D2p, MCTYPE);

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorCancelOutObtainWeight(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, names.at(i), MCTYPE, true);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor, it is not applied now. it is for a systematic uncertainty
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = 1.0;
            if (nB02KLphi2KSKL < MyEPSILON) Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);
            else Correction_KSKLKL = corrector_phiKL.GetCorrectionFactorAtGeneric(nB02KLphi2KSKL);

            // B --> K* KL KL correction factor
            double Correction_KstarKLKL = corrector_KstarKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KstarKLKL + nB02KstarKLKL);

            // B --> X KL KL correction factor
            double Correction_XKLKL = corrector_XsKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KpKLKL_all_KpKLKL, nB2KSKLKL_all_KSKLKL, nB2KstarKLKL + nB02KstarKLKL, nKL_XKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_KS0 * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL;

            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            // apply FBDT correction
            if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) total_weight = total_weight * 1.0045;
            else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) total_weight = total_weight * 1.0460;
            else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) total_weight = total_weight * 0.9717;

            int myIndex = static_cast<int>(std::upper_bound(boundary.begin(), boundary.end(), Mxs_MC) - boundary.begin());

            hist_reco->Fill(Bsig_M, total_weight);
            hist_true->Fill(myIndex - 0.5, total_weight);

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

            int myIndex = static_cast<int>(std::upper_bound(boundary.begin(), boundary.end(), Mxs_MC) - boundary.begin());

            hist_true->Fill(myIndex - 0.5, total_weight);

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

    const char* dirname_CHG_after_preselection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/CHG_analysis/validation_v009/before_Mbc_cut";
    const char* dirname_MIX_after_preselection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/MIX_analysis/validation_v009/before_Mbc_cut";
    const char* dirname_UUBAR_after_preselection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/UUBAR_analysis/validation_v009/before_Mbc_cut";
    const char* dirname_DDBAR_after_preselection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/DDBAR_analysis/validation_v009/before_Mbc_cut";
    const char* dirname_SSBAR_after_preselection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/SSBAR_analysis/validation_v009/before_Mbc_cut";
    const char* dirname_CHARM_after_preselection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/CHARM_analysis/validation_v009/before_Mbc_cut";

    const char* dirname_CHG_after_before_FBDT_cut = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/CHG_analysis/validation_v009/final_output";
    const char* dirname_MIX_after_before_FBDT_cut = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/MIX_analysis/validation_v009/final_output";
    const char* dirname_UUBAR_after_before_FBDT_cut = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/UUBAR_analysis/validation_v009/final_output";
    const char* dirname_DDBAR_after_before_FBDT_cut = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/DDBAR_analysis/validation_v009/final_output";
    const char* dirname_SSBAR_after_before_FBDT_cut = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/SSBAR_analysis/validation_v009/final_output";
    const char* dirname_CHARM_after_before_FBDT_cut = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/CHARM_analysis/validation_v009/final_output";

    const char* dirname_CHG_after_selection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/CHG_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut/Merge";
    const char* dirname_MIX_after_selection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/MIX_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut/Merge";
    const char* dirname_UUBAR_after_selection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/UUBAR_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut/Merge";
    const char* dirname_DDBAR_after_selection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/DDBAR_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut/Merge";
    const char* dirname_SSBAR_after_selection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/SSBAR_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut/Merge";
    const char* dirname_CHARM_after_selection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/CHARM_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut/Merge";

    const char* dirname_SIGNAL_initial = "/home/belle2/junewoo/storage_ghi/20220929_SIGNAL_decayInfo_again/small";
    const char* dirname_SIGNAL_after_preselection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/SIGNAL_analysis/validation_v009/before_Mbc_cut";
    const char* dirname_SIGNAL_after_before_FBDT_cut = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/SIGNAL_analysis/validation_v009/final_output";
    const char* dirname_SIGNAL_after_selection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/SIGNAL_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut/Merge";

    //FillBKG(dirname_CHG_after_preselection, "root", Ncandidate_BKG_after_preselection, ObtainWeight("CHG", MCTYPE, "validation", "CHG"));
    //FillBKG(dirname_MIX_after_preselection, "root", Ncandidate_BKG_after_preselection, ObtainWeight("MIX", MCTYPE, "validation", "MIX"));
    //FillBKG(dirname_UUBAR_after_preselection, "root", Ncandidate_BKG_after_preselection, ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"));
    //FillBKG(dirname_DDBAR_after_preselection, "root", Ncandidate_BKG_after_preselection, ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"));
    //FillBKG(dirname_SSBAR_after_preselection, "root", Ncandidate_BKG_after_preselection, ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"));
    //FillBKG(dirname_CHARM_after_preselection, "root", Ncandidate_BKG_after_preselection, ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"));

    //FillBKG(dirname_CHG_after_before_FBDT_cut, "root", Nevt_BKG_before_FBDT_cut, ObtainWeight("CHG", MCTYPE, "validation", "CHG"));
    //FillBKG(dirname_MIX_after_before_FBDT_cut, "root", Nevt_BKG_before_FBDT_cut, ObtainWeight("MIX", MCTYPE, "validation", "MIX"));
    //FillBKG(dirname_UUBAR_after_before_FBDT_cut, "root", Nevt_BKG_before_FBDT_cut, ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"));
    //FillBKG(dirname_DDBAR_after_before_FBDT_cut, "root", Nevt_BKG_before_FBDT_cut, ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"));
    //FillBKG(dirname_SSBAR_after_before_FBDT_cut, "root", Nevt_BKG_before_FBDT_cut, ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"));
    //FillBKG(dirname_CHARM_after_before_FBDT_cut, "root", Nevt_BKG_before_FBDT_cut, ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"));

    //FillBKG(dirname_CHG_after_selection, "root", Nevt_BKG_after_selection, ObtainWeight("CHG", MCTYPE, "validation", "CHG"));
    //FillBKG(dirname_MIX_after_selection, "root", Nevt_BKG_after_selection, ObtainWeight("MIX", MCTYPE, "validation", "MIX"));
    //FillBKG(dirname_UUBAR_after_selection, "root", Nevt_BKG_after_selection, ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"));
    //FillBKG(dirname_DDBAR_after_selection, "root", Nevt_BKG_after_selection, ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"));
    //FillBKG(dirname_SSBAR_after_selection, "root", Nevt_BKG_after_selection, ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"));
    //FillBKG(dirname_CHARM_after_selection, "root", Nevt_BKG_after_selection, ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"));

    //FillSIGNAL(dirname_SIGNAL_after_preselection, "B2Knunu", "Bplus", Ncandidate_SIGNAL_after_preselection, Nevt_SIGNAL_after_preselection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    //FillSIGNAL(dirname_SIGNAL_after_preselection, "B2Kstarnunu", "Bplus", Ncandidate_SIGNAL_after_preselection, Nevt_SIGNAL_after_preselection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    //FillSIGNAL(dirname_SIGNAL_after_preselection, "B2Xsnunu", "Bplus", Ncandidate_SIGNAL_after_preselection, Nevt_SIGNAL_after_preselection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    //FillSIGNAL(dirname_SIGNAL_after_preselection, "B02K0nunu", "Bzero", Ncandidate_SIGNAL_after_preselection, Nevt_SIGNAL_after_preselection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    //FillSIGNAL(dirname_SIGNAL_after_preselection, "B02Kstar0nunu", "Bzero", Ncandidate_SIGNAL_after_preselection, Nevt_SIGNAL_after_preselection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    //FillSIGNAL(dirname_SIGNAL_after_preselection, "B02Xsnunu", "Bzero", Ncandidate_SIGNAL_after_preselection, Nevt_SIGNAL_after_preselection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    //FillSIGNAL(dirname_SIGNAL_after_before_FBDT_cut, "B2Knunu", "Bplus", Nevt_SIGNAL_before_FBDT_cut, Nevt_SIGNAL_before_FBDT_cut_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    //FillSIGNAL(dirname_SIGNAL_after_before_FBDT_cut, "B2Kstarnunu", "Bplus", Nevt_SIGNAL_before_FBDT_cut, Nevt_SIGNAL_before_FBDT_cut_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    //FillSIGNAL(dirname_SIGNAL_after_before_FBDT_cut, "B2Xsnunu", "Bplus", Nevt_SIGNAL_before_FBDT_cut, Nevt_SIGNAL_before_FBDT_cut_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    //FillSIGNAL(dirname_SIGNAL_after_before_FBDT_cut, "B02K0nunu", "Bzero", Nevt_SIGNAL_before_FBDT_cut, Nevt_SIGNAL_before_FBDT_cut_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    //FillSIGNAL(dirname_SIGNAL_after_before_FBDT_cut, "B02Kstar0nunu", "Bzero", Nevt_SIGNAL_before_FBDT_cut, Nevt_SIGNAL_before_FBDT_cut_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    //FillSIGNAL(dirname_SIGNAL_after_before_FBDT_cut, "B02Xsnunu", "Bzero", Nevt_SIGNAL_before_FBDT_cut, Nevt_SIGNAL_before_FBDT_cut_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    FillSIGNAL(dirname_SIGNAL_after_selection, "B2Knunu", "Bplus", Nevt_SIGNAL_after_selection, Nevt_SIGNAL_after_selection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    FillSIGNAL(dirname_SIGNAL_after_selection, "B2Kstarnunu", "Bplus", Nevt_SIGNAL_after_selection, Nevt_SIGNAL_after_selection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    FillSIGNAL(dirname_SIGNAL_after_selection, "B2Xsnunu", "Bplus", Nevt_SIGNAL_after_selection, Nevt_SIGNAL_after_selection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    FillSIGNAL(dirname_SIGNAL_after_selection, "B02K0nunu", "Bzero", Nevt_SIGNAL_after_selection, Nevt_SIGNAL_after_selection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    FillSIGNAL(dirname_SIGNAL_after_selection, "B02Kstar0nunu", "Bzero", Nevt_SIGNAL_after_selection, Nevt_SIGNAL_after_selection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    FillSIGNAL(dirname_SIGNAL_after_selection, "B02Xsnunu", "Bzero", Nevt_SIGNAL_after_selection, Nevt_SIGNAL_after_selection_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    //FillSIGNAL_q2(dirname_SIGNAL_after_selection, "B2Knunu", "Bplus", Nevt_SIGNAL_after_selection_q2_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    //FillSIGNAL_q2(dirname_SIGNAL_after_selection, "B2Kstarnunu", "Bplus", Nevt_SIGNAL_after_selection_q2_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    //FillSIGNAL_q2(dirname_SIGNAL_after_selection, "B2Xsnunu", "Bplus", Nevt_SIGNAL_after_selection_q2_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    //FillSIGNAL_q2(dirname_SIGNAL_after_selection, "B02K0nunu", "Bzero", Nevt_SIGNAL_after_selection_q2_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    //FillSIGNAL_q2(dirname_SIGNAL_after_selection, "B02Kstar0nunu", "Bzero", Nevt_SIGNAL_after_selection_q2_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    //FillSIGNAL_q2(dirname_SIGNAL_after_selection, "B02Xsnunu", "Bzero", Nevt_SIGNAL_after_selection_q2_true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    ReadDecayInfo(dirname_SIGNAL_initial, "B2Knunu", "Bplus", Nevt_SIGNAL_initial_true, Scale_Kplus_here, "B2Knunu");
    ReadDecayInfo(dirname_SIGNAL_initial, "B2Kstarnunu", "Bplus", Nevt_SIGNAL_initial_true, Scale_Kplusstar_here, "otherwise");
    ReadDecayInfo(dirname_SIGNAL_initial, "B2Xsnunu", "Bplus", Nevt_SIGNAL_initial_true, Scale_Xsu_nonresonant_here, "B2Xsnunu");
    ReadDecayInfo(dirname_SIGNAL_initial, "B02K0nunu", "Bzero", Nevt_SIGNAL_initial_true, Scale_K0_here, "B02K0nunu");
    ReadDecayInfo(dirname_SIGNAL_initial, "B02Kstar0nunu", "Bzero", Nevt_SIGNAL_initial_true, Scale_K0star_here, "otherwise");
    ReadDecayInfo(dirname_SIGNAL_initial, "B02Xsnunu", "Bzero", Nevt_SIGNAL_initial_true, Scale_Xsd_nonresonant_here, "B02Xsnunu");

    //ReadDecayInfo_q2(dirname_SIGNAL_initial, "B2Knunu", "Bplus", Nevt_SIGNAL_initial_q2_true, Scale_Kplus_here, "B2Knunu");
    //ReadDecayInfo_q2(dirname_SIGNAL_initial, "B2Kstarnunu", "Bplus", Nevt_SIGNAL_initial_q2_true, Scale_Kplusstar_here, "otherwise");
    //ReadDecayInfo_q2(dirname_SIGNAL_initial, "B2Xsnunu", "Bplus", Nevt_SIGNAL_initial_q2_true, Scale_Xsu_nonresonant_here, "B2Xsnunu");
    //ReadDecayInfo_q2(dirname_SIGNAL_initial, "B02K0nunu", "Bzero", Nevt_SIGNAL_initial_q2_true, Scale_K0_here, "B02K0nunu");
    //ReadDecayInfo_q2(dirname_SIGNAL_initial, "B02Kstar0nunu", "Bzero", Nevt_SIGNAL_initial_q2_true, Scale_K0star_here, "otherwise");
    //ReadDecayInfo_q2(dirname_SIGNAL_initial, "B02Xsnunu", "Bzero", Nevt_SIGNAL_initial_q2_true, Scale_Xsd_nonresonant_here, "B02Xsnunu");

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

    for (int i = 0; i < 10; i++) {
        printf("%.9lf,", SIGNAL_pre_and_selection_efficiency_true->GetBinContent(i + 1));
    }
    printf("\n");
    for (int i = 0; i < 10; i++) {
        printf("%.9lf,", SIGNAL_pre_and_selection_efficiency_true->GetBinError(i + 1));
    }
    printf("\n");

    return 0;
}
