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
Corrector_Multiplicity corrector_Multiplicity;
Corrector_Fragmentation corrector_Fragmentation;

# define MCTYPE "MC15ri"

# define Scale_Kplus_here (N_Kplus_nunubar_LS1/1000000.0)
# define Scale_Kplusstar_here (N_Kplusstar_nunubar_LS1/1000000.0)
# define Scale_Xsu_nonresonant_here (N_Xsu_nonresonant_nunubar_LS1/5000000.0)
# define Scale_K0_here (N_K0_nunubar_LS1/1000000.0)
# define Scale_K0star_here (N_K0star_nunubar_LS1/1000000.0)
# define Scale_Xsd_nonresonant_here (N_Xsd_nunubar_LS1/5000000.0)

/* ====================================== */

bool AcceptIt(int Decay[N_decay], int Xsu_KpipiToTree[N_Xsu_Kpipi], int Xsu_KpipipiToTree[N_Xsu_Kpipipi], int Xsd_KpipiToTree[N_Xsd_Kpipi], int Xsd_KpipipiToTree[N_Xsd_Kpipipi], int index) {

    // accept decay with index. Only works for non-resonant Xs sample

    if (Decay[7] > 0.5) { // B+ -> K+ pi- pi+

        if (index == 63) return true;

        if (Xsu_KpipiToTree[0] > 0.5) { // B+ -> rho0 K+
            if (index == 0) return true;
            else return false;
        }
        else if (Xsu_KpipiToTree[1] > 0.5) { // B+ -> pi+ K*0
            if (index == 1) return true;
            else return false;
        }
        else if (Xsu_KpipiToTree[2] > 0.5) { // B+ -> K+ pi- pi+
            if (index == 2) return true;
            else return false;
        }
        else if (Xsu_KpipiToTree[3] > 0.5) { // B+ -> K+ omega
            if (index == 3) return true;
            else return false;
        }
    }
    else if (Decay[8] > 0.5) { // B+ -> K0 pi- pi0

        if (index == 64) return true;

        if (Xsu_KpipiToTree[4] > 0.5) { // B+ -> rho+ K0
            if (index == 4) return true;
            else return false;
        }
        else if (Xsu_KpipiToTree[5] > 0.5) { // B+ -> pi+ K*0
            if (index == 5) return true;
            else return false;
        }
        else if (Xsu_KpipiToTree[6] > 0.5) { // B+ -> K0 pi0 pi+
            if (index == 6) return true;
            else return false;
        }
        else if (Xsu_KpipiToTree[7] > 0.5) { // B+ -> pi0 K*+
            if (index == 7) return true;
            else return false;
        }
    }
    else if (Decay[27] > 0.5) { // B0 -> K0 pi+ pi-

        if (index == 65) return true;

        if (Xsd_KpipiToTree[4] > 0.5) { // B+ -> rho0 K0
            if (index == 8) return true;
            else return false;
        }
        else if (Xsd_KpipiToTree[5] > 0.5) { // B+ -> pi- K*+
            if (index == 9) return true;
            else return false;
        }
        else if (Xsd_KpipiToTree[6] > 0.5) { // B+ -> pi- pi+ K0
            if (index == 10) return true;
            else return false;
        }
        else if (Xsd_KpipiToTree[7] > 0.5) { // B+ -> omega K0
            if (index == 11) return true;
            else return false;
        }
    }
    else if (Decay[26] > 0.5) { // B0 -> K+ pi- pi0

        if (index == 66) return true;

        if (Xsd_KpipiToTree[0] > 0.5) { // B+ -> rho- K+
            if (index == 12) return true;
            else return false;
        }
        else if (Xsd_KpipiToTree[1] > 0.5) { // B+ -> pi- K*+
            if (index == 13) return true;
            else return false;
        }
        else if (Xsd_KpipiToTree[2] > 0.5) { // B+ -> K+ pi- pi0
            if (index == 14) return true;
            else return false;
        }
        else if (Xsd_KpipiToTree[3] > 0.5) { // B+ -> pi0 K*0:
            if (index == 15) return true;
            else return false;
        }
    }
    else if (Decay[9] > 0.5) { // B+ -> K+ pi- pi+ pi0

        if (index == 67) return true;

        if (Xsu_KpipipiToTree[0] > 0.5) { // B+ -> rho+ pi- K+
            if (index == 16) return true;
            else return false;
        }
        else if (Xsu_KpipipiToTree[1] > 0.5) { // B+ -> K+ eta
            if (index == 17) return true;
            else return false;
        }
        else if (Xsu_KpipipiToTree[2] > 0.5) { // B+ -> omega K+
            if (index == 18) return true;
            else return false;
        }
        else if (Xsu_KpipipiToTree[3] > 0.5) { // B+ -> pi- pi+ pi0 K+
            if (index == 19) return true;
            else return false;
        }
        else if (Xsu_KpipipiToTree[4] > 0.5) { // B+ -> K+ phi
            if (index == 20) return true;
            else return false;
        }
        else if (Xsu_KpipipiToTree[5] > 0.5) { // B+ -> rho+ K*0
            if (index == 21) return true;
            else return false;
        }
        else if (Xsu_KpipipiToTree[6] > 0.5) { // B+ -> pi+ pi0 K*0
            if (index == 22) return true;
            else return false;
        }
        else if (Xsu_KpipipiToTree[7] > 0.5) { // B+ -> rho- pi+ K+
            if (index == 23) return true;
            else return false;
        }
        else if (Xsu_KpipipiToTree[8] > 0.5) { // B+ -> rho0 pi0 K+
            if (index == 24) return true;
            else return false;
        }
        else if (Xsu_KpipipiToTree[9] > 0.5) { // B+ -> rho0 K*+
            if (index == 25) return true;
            else return false;
        }
        else if (Xsu_KpipipiToTree[10] > 0.5) { // B+ -> pi- pi+ K*+
            if (index == 26) return true;
            else return false;
        }
        else if (Xsu_KpipipiToTree[11] > 0.5) { // B+ -> omega K+ pi0
            if (index == 27) return true;
            else return false;
        }
        else if (Xsu_KpipipiToTree[12] > 0.5) { // B+ -> omega K*+
            if (index == 28) return true;
            else return false;
        }
    }
    else if (Decay[10] > 0.5) { // B+ -> K0 pi+ pi- pi+

        if (index == 68) return true;

        if (Xsu_KpipipiToTree[13] > 0.5) { // B+ -> K*+ pi- pi+
            if (index == 29) return true;
            else return false;
        }
        else if (Xsu_KpipipiToTree[14] > 0.5) { // B+ -> K0 pi+ pi- pi+
            if (index == 30) return true;
            else return false;
        }
        else if (Xsu_KpipipiToTree[15] > 0.5) { // B+ -> K*+ rho0
            if (index == 31) return true;
            else return false;
        }
        else if (Xsu_KpipipiToTree[17] > 0.5) { // B+ -> K0 pi+ rho0
            if (index == 32) return true;
            else return false;
        }
        else if (Xsu_KpipipiToTree[18] > 0.5) { // B+ -> K*+ omega
            if (index == 33) return true;
            else return false;
        }
        else if (Xsu_KpipipiToTree[19] > 0.5) { // B+ -> K0 pi+ omega
            if (index == 34) return true;
            else return false;
        }
    }
    else if (Decay[14] > 0.5) { // B+ -> K0 pi+ pi0 pi0

        if (index == 69) return true;

        if (Xsu_KpipipiToTree[20] > 0.5) { // B+ -> pi0 rho+ K0
            if (index == 35) return true;
            else return false;
        }
        else if (Xsu_KpipipiToTree[21] > 0.5) { // B+ -> K*0 pi+ pi0
            if (index == 36) return true;
            else return false;
        }
        else if (Xsu_KpipipiToTree[22] > 0.5) { // B+ -> K*0 rho+
            if (index == 37) return true;
            else return false;
        }
        else if (Xsu_KpipipiToTree[23] > 0.5) { // B+ -> K*+ pi0 pi0
            if (index == 38) return true;
            else return false;
        }
        else if (Xsu_KpipipiToTree[24] > 0.5) { // B+ -> K0 pi+ pi0 pi0
            if (index == 39) return true;
            else return false;
        }
    }
    else if (Decay[29] > 0.5) { // B0 -> K0 pi+ pi- pi0

        if (index == 70) return true;

        if (Xsd_KpipipiToTree[6] > 0.5) { // B0 -> K0 omega
            if (index == 40) return true;
            else return false;
        }
        else if (Xsd_KpipipiToTree[7] > 0.5) { // B0 -> K0 pi0 rho0
            if (index == 41) return true;
            else return false;
        }
        else if (Xsd_KpipipiToTree[8] > 0.5) { // B0 -> K0 pi+ pi- pi0
            if (index == 42) return true;
            else return false;
        }
        else if (Xsd_KpipipiToTree[9] > 0.5) { // B0 -> K*0 pi+ pi-
            if (index == 43) return true;
            else return false;
        }
        else if (Xsd_KpipipiToTree[10] > 0.5) { // B0 -> K0 pi- rho+
            if (index == 44) return true;
            else return false;
        }
        else if (Xsd_KpipipiToTree[11] > 0.5) { // B0 -> K*+ rho-
            if (index == 45) return true;
            else return false;
        }
        else if (Xsd_KpipipiToTree[12] > 0.5) { // B0 -> K0 phi
            if (index == 46) return true;
            else return false;
        }
        else if (Xsd_KpipipiToTree[13] > 0.5) { // B0 -> K*+ pi0 pi-
            if (index == 47) return true;
            else return false;
        }
        else if (Xsd_KpipipiToTree[14] > 0.5) { // B0 -> K*0 rho0
            if (index == 48) return true;
            else return false;
        }
        else if (Xsd_KpipipiToTree[15] > 0.5) { // B0 -> K0 eta
            if (index == 49) return true;
            else return false;
        }
        else if (Xsd_KpipipiToTree[16] > 0.5) { // B0 -> K0 pi+ rho-
            if (index == 50) return true;
            else return false;
        }
        else if (Xsd_KpipipiToTree[17] > 0.5) { // B0 -> K0 pi0 omega
            if (index == 51) return true;
            else return false;
        }
    }
    else if (Decay[28] > 0.5) { // B0 -> K+ pi- pi+ pi-

        if (index == 71) return true;

        if (Xsd_KpipipiToTree[0] > 0.5) { // B0 -> rho0 K+ pi-
            if (index == 52) return true;
            else return false;
        }
        else if (Xsd_KpipipiToTree[1] > 0.5) { // B0 -> K*0 pi+ pi-
            if (index == 53) return true;
            else return false;
        }
        else if (Xsd_KpipipiToTree[2] > 0.5) { // B0 -> K*0 rho0
            if (index == 54) return true;
            else return false;
        }
        else if (Xsd_KpipipiToTree[3] > 0.5) { // B0 -> pi- pi+ pi- K+
            if (index == 55) return true;
            else return false;
        }
        else if (Xsd_KpipipiToTree[4] > 0.5) { // B0 -> omega K*0
            if (index == 56) return true;
            else return false;
        }
        else if (Xsd_KpipipiToTree[5] > 0.5) { // B0 -> omega K+ pi-
            if (index == 57) return true;
            else return false;
        }
    }
    else if (Decay[33] > 0.5) { // B0 -> K+ pi- pi0 pi0

        if (index == 72) return true;

        if (Xsd_KpipipiToTree[18] > 0.5) { // B0 -> K+ rho- pi0
            if (index == 58) return true;
            else return false;
        }
        else if (Xsd_KpipipiToTree[19] > 0.5) { // B0 -> K*+ pi- pi0
            if (index == 59) return true;
            else return false;
        }
        else if (Xsd_KpipipiToTree[20] > 0.5) { // B0 -> K+ pi- pi0 pi0
            if (index == 60) return true;
            else return false;
        }
        else if (Xsd_KpipipiToTree[21] > 0.5) { // B0 -> K*+ rho-
            if (index == 61) return true;
            else return false;
        }
        else if (Xsd_KpipipiToTree[22] > 0.5) { // B0 -> K*0 pi0 pi0
            if (index == 62) return true;
            else return false;
        }
    }

    return false;
}

int INDEX = -1;

TH1D* Ncandidate_BKG_after_preselection = new TH1D("Ncandidate_BKG_after_preselection", ";M_{X_{s}}^{reco} [GeV/c^{2}]", 9, 0.45, 3.0);
TH1D* Nevt_BKG_before_FBDT_cut = new TH1D("Nevt_BKG_before_FBDT_cut", ";M_{X_{s}}^{reco} [GeV/c^{2}]", 9, 0.45, 3.0);
TH1D* Nevt_BKG_after_selection = new TH1D("Nevt_BKG_after_selection", ";M_{X_{s}}^{reco} [GeV/c^{2}]", 9, 0.45, 3.0);

TH1D* Ncandidate_SIGNAL_after_preselection = new TH1D("Ncandidate_SIGNAL_after_preselection", ";M_{X_{s}}^{reco} [GeV/c^{2}]", 9, 0.45, 3.0);
TH1D* Nevt_SIGNAL_before_FBDT_cut = new TH1D("Nevt_SIGNAL_before_FBDT_cut", ";M_{X_{s}}^{reco} [GeV/c^{2}]", 9, 0.45, 3.0);
TH1D* Nevt_SIGNAL_after_selection = new TH1D("Nevt_SIGNAL_after_selection", ";M_{X_{s}}^{reco} [GeV/c^{2}]", 9, 0.45, 3.0);

TH1D* Nevt_SIGNAL_initial_true = new TH1D("Nevt_SIGNAL_initial_true", ";M_{X_{s}}^{true} [GeV/c^{2}]", 9, 1.1, 4.0);
TH1D* Nevt_SIGNAL_after_preselection_true = new TH1D("Nevt_SIGNAL_after_preselection_true", ";M_{X_{s}}^{true} [GeV/c^{2}]", 9, 1.1, 4.0);
TH1D* Nevt_SIGNAL_before_FBDT_cut_true = new TH1D("Nevt_SIGNAL_before_FBDT_cut_true", ";M_{X_{s}}^{true} [GeV/c^{2}]", 9, 1.1, 4.0);
TH1D* Nevt_SIGNAL_after_selection_true = new TH1D("Nevt_SIGNAL_after_selection_true", ";M_{X_{s}}^{true} [GeV/c^{2}]", 9, 1.1, 4.0);

TH1D* SIGNAL_pre_and_selection_efficiency_true = new TH1D("SIGNAL_pre_and_selection_efficiency_true", ";M_{X_{s}}^{true} [GeV/c^{2}];Efficiency", 9, 1.1, 4.0);
TH1D* SIGNAL_selection_efficiency = new TH1D("SIGNAL_selection_efficiency", ";M_{X_{s}}^{reco} [GeV/c^{2}];Efficiency", 9, 1.1, 4.0);
TH1D* SIGNAL_FBDT_efficiency = new TH1D("SIGNAL_FBDT_efficiency", ";M_{X_{s}}^{reco} [GeV/c^{2}];FBDT Efficiency", 9, 1.1, 4.0);
TH1D* SIGNAL_FBDT_efficiency_true = new TH1D("SIGNAL_FBDT_efficiency_true", ";M_{X_{s}}^{true} [GeV/c^{2}];FBDT Efficiency", 9, 1.1, 4.0);

TH1D* BKG_selection_efficiency = new TH1D("BKG_selection_efficiency", ";M_{X_{s}}^{reco} [GeV/c^{2}];1 - (rejection rate)", 9, 0.45, 3.0);
TH1D* BKG_FBDT_efficiency = new TH1D("BKG_FBDT_efficiency", ";M_{X_{s}}^{reco} [GeV/c^{2}];FBDT 1 - (rejection rate)", 9, 0.45, 3.0);

TH1D* Nevt_SIGNAL_initial_q2_true = new TH1D("Nevt_SIGNAL_initial_q2_true", ";q^{2} [GeV^{2}/c^{4}]", 10, 0.0, 22.9);
TH1D* Nevt_SIGNAL_after_selection_q2_true = new TH1D("Nevt_SIGNAL_after_selection_q2_true", ";q^{2} [GeV^{2}/c^{4}]", 10, 0.0, 22.9);
TH1D* SIGNAL_efficiency_q2_true = new TH1D("SIGNAL_efficiency_q2_true", ";q^{2} [GeV^{2}/c^{4}];Efficiency", 10, 0.0, 22.9);

TH1D* Nevt_SIGNAL_initial_true_all = new TH1D("Nevt_SIGNAL_initial_true_all", "", 1, 0.0, 100.0);
TH1D* Nevt_SIGNAL_after_selection_true_all = new TH1D("Nevt_SIGNAL_after_selection_true_all", "", 1, 0.0, 100.0);
TH1D* SIGNAL_pre_and_selection_efficiency_true_all = new TH1D("SIGNAL_pre_and_selection_efficiency_true_all", "", 1, 0.0, 100.0);

void FillBKG(const char* dirname, const char* included_string, TH1D* hist, const char* type, double weight_var = 1.0) { // get nominal PDF with appropriate correction
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    int __experiment__;
    int __run__;
    unsigned int __event__;
    int __candidate__;
    int __ncandidates__;
    double Bsig_M = -1;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

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
        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);
        for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKtruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[0][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKmisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[1][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npitruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[2][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npimisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[3][i_PID]);
        }
        for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npi0bin" + std::to_string(i_pi0)).c_str(), &temp_N_bin_pi0[i_pi0]);
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

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE), temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE), temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
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

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * corrector_Knn.GetCorrectionFactorCancelOutObtainWeight(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, names.at(i), MCTYPE, false);

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
    int Xsu_KpipiToTree[N_Xsu_Kpipi] = { 0 };
    int Xsu_KpipipiToTree[N_Xsu_Kpipipi] = { 0 };
    int Xsd_KpipiToTree[N_Xsd_Kpipi] = { 0 };
    int Xsd_KpipipiToTree[N_Xsd_Kpipipi] = { 0 };
    double Mxs_MC = -1.0;
    double invM = -1.0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

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
        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);
        for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKtruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[0][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKmisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[1][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npitruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[2][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npimisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[3][i_PID]);
        }
        for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npi0bin" + std::to_string(i_pi0)).c_str(), &temp_N_bin_pi0[i_pi0]);
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
        for (int i = 0; i < N_Xsu_Kpipi; i++) tree_Xs->SetBranchAddress(("nParticlesInList__boXsu__clKpipi_ch" + std::to_string(i + 1) + "__bc").c_str(), &Xsu_KpipiToTree[i]);
        for (int i = 0; i < N_Xsu_Kpipipi; i++) tree_Xs->SetBranchAddress(("nParticlesInList__boXsu__clKpipipi_ch" + std::to_string(i + 1) + "__bc").c_str(), &Xsu_KpipipiToTree[i]);
        for (int i = 0; i < N_Xsd_Kpipi; i++) tree_Xs->SetBranchAddress(("nParticlesInList__boXsd__clKpipi_ch" + std::to_string(i + 1) + "__bc").c_str(), &Xsd_KpipiToTree[i]);
        for (int i = 0; i < N_Xsd_Kpipipi; i++) tree_Xs->SetBranchAddress(("nParticlesInList__boXsd__clKpipipi_ch" + std::to_string(i + 1) + "__bc").c_str(), &Xsd_KpipipiToTree[i]);
        if (strcmp(type, "Bplus") == 0) tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_MC);
        else tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_MC);
        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            tree_Xs->GetEntry(j);

            // Here we get only 
            if (AcceptIt(Decay, Xsu_KpipiToTree, Xsu_KpipipiToTree, Xsd_KpipiToTree, Xsd_KpipipiToTree, INDEX)) {}
            else continue;

            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE), temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE), temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
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

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake;

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
    int Xsu_KpipiToTree[N_Xsu_Kpipi] = { 0 };
    int Xsu_KpipipiToTree[N_Xsu_Kpipipi] = { 0 };
    int Xsd_KpipiToTree[N_Xsd_Kpipi] = { 0 };
    int Xsd_KpipipiToTree[N_Xsd_Kpipipi] = { 0 };
    double Mxs_MC = -1.0;
    double invM = -1.0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

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
        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);
        for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKtruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[0][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKmisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[1][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npitruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[2][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npimisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[3][i_PID]);
        }
        for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npi0bin" + std::to_string(i_pi0)).c_str(), &temp_N_bin_pi0[i_pi0]);
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
        for (int i = 0; i < N_Xsu_Kpipi; i++) tree_Xs->SetBranchAddress(("nParticlesInList__boXsu__clKpipi_ch" + std::to_string(i + 1) + "__bc").c_str(), &Xsu_KpipiToTree[i]);
        for (int i = 0; i < N_Xsu_Kpipipi; i++) tree_Xs->SetBranchAddress(("nParticlesInList__boXsu__clKpipipi_ch" + std::to_string(i + 1) + "__bc").c_str(), &Xsu_KpipipiToTree[i]);
        for (int i = 0; i < N_Xsd_Kpipi; i++) tree_Xs->SetBranchAddress(("nParticlesInList__boXsd__clKpipi_ch" + std::to_string(i + 1) + "__bc").c_str(), &Xsd_KpipiToTree[i]);
        for (int i = 0; i < N_Xsd_Kpipipi; i++) tree_Xs->SetBranchAddress(("nParticlesInList__boXsd__clKpipipi_ch" + std::to_string(i + 1) + "__bc").c_str(), &Xsd_KpipipiToTree[i]);
        if(strcmp(type, "Bplus") == 0) tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_MC);
        else tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_MC);
        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            tree_Xs->GetEntry(j);

            // Here we get only 
            if (AcceptIt(Decay, Xsu_KpipiToTree, Xsu_KpipipiToTree, Xsd_KpipiToTree, Xsd_KpipipiToTree, INDEX)) {}
            else continue;

            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE), temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE), temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
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

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake;

            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            hist_reco->Fill(Bsig_M, total_weight);
            hist_true->Fill(Mxs_MC, total_weight);

            Nevt = Nevt + total_weight;

            Nevt_SIGNAL_after_selection_true_all->Fill(Mxs_MC, total_weight);

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
    int Xsu_KpipiToTree[N_Xsu_Kpipi] = { 0 };
    int Xsu_KpipipiToTree[N_Xsu_Kpipipi] = { 0 };
    int Xsd_KpipiToTree[N_Xsd_Kpipi] = { 0 };
    int Xsd_KpipipiToTree[N_Xsd_Kpipipi] = { 0 };
    double Mxs_MC = -1.0;
    double invM = -1.0;

    std::vector<std::string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_Xs = (TTree*)input_file->Get("decay");

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
        for (int i = 0; i < N_Xsu_Kpipi; i++) tree_Xs->SetBranchAddress(("nParticlesInList__boXsu__clKpipi_ch" + std::to_string(i + 1) + "__bc").c_str(), &Xsu_KpipiToTree[i]);
        for (int i = 0; i < N_Xsu_Kpipipi; i++) tree_Xs->SetBranchAddress(("nParticlesInList__boXsu__clKpipipi_ch" + std::to_string(i + 1) + "__bc").c_str(), &Xsu_KpipipiToTree[i]);
        for (int i = 0; i < N_Xsd_Kpipi; i++) tree_Xs->SetBranchAddress(("nParticlesInList__boXsd__clKpipi_ch" + std::to_string(i + 1) + "__bc").c_str(), &Xsd_KpipiToTree[i]);
        for (int i = 0; i < N_Xsd_Kpipipi; i++) tree_Xs->SetBranchAddress(("nParticlesInList__boXsd__clKpipipi_ch" + std::to_string(i + 1) + "__bc").c_str(), &Xsd_KpipipiToTree[i]);
        if (strcmp(type, "Bplus") == 0) tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_MC);
        else tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_MC);
        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

        printf("%lld entries...\n", tree_Xs->GetEntries());
        for (unsigned int j = 0; j < tree_Xs->GetEntries(); j++) { // Fill
            tree_Xs->GetEntry(j);

            // Here we get only 
            if (AcceptIt(Decay, Xsu_KpipiToTree, Xsu_KpipipiToTree, Xsd_KpipiToTree, Xsd_KpipipiToTree, INDEX)) {}
            else continue;

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

            Nevt_SIGNAL_initial_true_all->Fill(Mxs_MC, total_weight);

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
    int Xsu_KpipiToTree[N_Xsu_Kpipi] = { 0 };
    int Xsu_KpipipiToTree[N_Xsu_Kpipipi] = { 0 };
    int Xsd_KpipiToTree[N_Xsd_Kpipi] = { 0 };
    int Xsd_KpipipiToTree[N_Xsd_Kpipipi] = { 0 };
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
        for (int i = 0; i < N_Xsu_Kpipi; i++) tree_Xs->SetBranchAddress(("nParticlesInList__boXsu__clKpipi_ch" + std::to_string(i + 1) + "__bc").c_str(), &Xsu_KpipiToTree[i]);
        for (int i = 0; i < N_Xsu_Kpipipi; i++) tree_Xs->SetBranchAddress(("nParticlesInList__boXsu__clKpipipi_ch" + std::to_string(i + 1) + "__bc").c_str(), &Xsu_KpipipiToTree[i]);
        for (int i = 0; i < N_Xsd_Kpipi; i++) tree_Xs->SetBranchAddress(("nParticlesInList__boXsd__clKpipi_ch" + std::to_string(i + 1) + "__bc").c_str(), &Xsd_KpipiToTree[i]);
        for (int i = 0; i < N_Xsd_Kpipipi; i++) tree_Xs->SetBranchAddress(("nParticlesInList__boXsd__clKpipipi_ch" + std::to_string(i + 1) + "__bc").c_str(), &Xsd_KpipipiToTree[i]);
        if (strcmp(type, "Bplus") == 0) tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_MC);
        else tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_MC);
        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

        printf("%lld entries...\n", tree_Xs->GetEntries());
        for (unsigned int j = 0; j < tree_Xs->GetEntries(); j++) { // Fill
            tree_Xs->GetEntry(j);

            // Here we get only 
            if (AcceptIt(Decay, Xsu_KpipiToTree, Xsu_KpipipiToTree, Xsd_KpipiToTree, Xsd_KpipipiToTree, INDEX)) {}
            else continue;

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
    INDEX = std::atoi(argv[1]);

    const char* dirname_CHG_after_preselection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/CWR/CHG_analysis/test_v000/before_Mbc_cut";
    const char* dirname_MIX_after_preselection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/CWR/MIX_analysis/test_v000/before_Mbc_cut";
    const char* dirname_UUBAR_after_preselection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/CWR/UUBAR_analysis/test_v000/before_Mbc_cut";
    const char* dirname_DDBAR_after_preselection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/CWR/DDBAR_analysis/test_v000/before_Mbc_cut";
    const char* dirname_SSBAR_after_preselection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/CWR/SSBAR_analysis/test_v000/before_Mbc_cut";
    const char* dirname_CHARM_after_preselection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/CWR/CHARM_analysis/test_v000/before_Mbc_cut";

    const char* dirname_CHG_after_before_FBDT_cut = "/home/belle2/junewoo/storage_b1/bsub/Analysis/CWR/CHG_analysis/test_v000/final_output";
    const char* dirname_MIX_after_before_FBDT_cut = "/home/belle2/junewoo/storage_b1/bsub/Analysis/CWR/MIX_analysis/test_v000/final_output";
    const char* dirname_UUBAR_after_before_FBDT_cut = "/home/belle2/junewoo/storage_b1/bsub/Analysis/CWR/UUBAR_analysis/test_v000/final_output";
    const char* dirname_DDBAR_after_before_FBDT_cut = "/home/belle2/junewoo/storage_b1/bsub/Analysis/CWR/DDBAR_analysis/test_v000/final_output";
    const char* dirname_SSBAR_after_before_FBDT_cut = "/home/belle2/junewoo/storage_b1/bsub/Analysis/CWR/SSBAR_analysis/test_v000/final_output";
    const char* dirname_CHARM_after_before_FBDT_cut = "/home/belle2/junewoo/storage_b1/bsub/Analysis/CWR/CHARM_analysis/test_v000/final_output";

    const char* dirname_CHG_after_selection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/CWR/CHG_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/Merge";
    const char* dirname_MIX_after_selection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/CWR/MIX_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/Merge";
    const char* dirname_UUBAR_after_selection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/CWR/UUBAR_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/Merge";
    const char* dirname_DDBAR_after_selection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/CWR/DDBAR_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/Merge";
    const char* dirname_SSBAR_after_selection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/CWR/SSBAR_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/Merge";
    const char* dirname_CHARM_after_selection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/CWR/CHARM_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/Merge";

    const char* dirname_SIGNAL_initial = "/home/belle2/junewoo/storage_ghi/20220929_SIGNAL_gsim_CWR_Karim/output/Ntuple";
    const char* dirname_SIGNAL_after_preselection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/CWR/SIGNAL_analysis/test_v000/before_Mbc_cut";
    const char* dirname_SIGNAL_after_before_FBDT_cut = "/home/belle2/junewoo/storage_b1/bsub/Analysis/CWR/SIGNAL_analysis/test_v000/final_output";
    const char* dirname_SIGNAL_after_selection = "/home/belle2/junewoo/storage_b1/bsub/Analysis/CWR/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/Merge";
    /*
    FillBKG(dirname_CHG_after_preselection, "root", Ncandidate_BKG_after_preselection, ObtainWeight("CHG", MCTYPE, "test", "CHG"));
    FillBKG(dirname_MIX_after_preselection, "root", Ncandidate_BKG_after_preselection, ObtainWeight("MIX", MCTYPE, "test", "MIX"));
    FillBKG(dirname_UUBAR_after_preselection, "root", Ncandidate_BKG_after_preselection, ObtainWeight("UUBAR", MCTYPE, "test", "UUBAR"));
    FillBKG(dirname_DDBAR_after_preselection, "root", Ncandidate_BKG_after_preselection, ObtainWeight("DDBAR", MCTYPE, "test", "DDBAR"));
    FillBKG(dirname_SSBAR_after_preselection, "root", Ncandidate_BKG_after_preselection, ObtainWeight("SSBAR", MCTYPE, "test", "SSBAR"));
    FillBKG(dirname_CHARM_after_preselection, "root", Ncandidate_BKG_after_preselection, ObtainWeight("CHARM", MCTYPE, "test", "CHARM"));

    FillBKG(dirname_CHG_after_before_FBDT_cut, "root", Nevt_BKG_before_FBDT_cut, ObtainWeight("CHG", MCTYPE, "test", "CHG"));
    FillBKG(dirname_MIX_after_before_FBDT_cut, "root", Nevt_BKG_before_FBDT_cut, ObtainWeight("MIX", MCTYPE, "test", "MIX"));
    FillBKG(dirname_UUBAR_after_before_FBDT_cut, "root", Nevt_BKG_before_FBDT_cut, ObtainWeight("UUBAR", MCTYPE, "test", "UUBAR"));
    FillBKG(dirname_DDBAR_after_before_FBDT_cut, "root", Nevt_BKG_before_FBDT_cut, ObtainWeight("DDBAR", MCTYPE, "test", "DDBAR"));
    FillBKG(dirname_SSBAR_after_before_FBDT_cut, "root", Nevt_BKG_before_FBDT_cut, ObtainWeight("SSBAR", MCTYPE, "test", "SSBAR"));
    FillBKG(dirname_CHARM_after_before_FBDT_cut, "root", Nevt_BKG_before_FBDT_cut, ObtainWeight("CHARM", MCTYPE, "test", "CHARM"));

    FillBKG(dirname_CHG_after_selection, "root", Nevt_BKG_after_selection, ObtainWeight("CHG", MCTYPE, "test", "CHG"));
    FillBKG(dirname_MIX_after_selection, "root", Nevt_BKG_after_selection, ObtainWeight("MIX", MCTYPE, "test", "MIX"));
    FillBKG(dirname_UUBAR_after_selection, "root", Nevt_BKG_after_selection, ObtainWeight("UUBAR", MCTYPE, "test", "UUBAR"));
    FillBKG(dirname_DDBAR_after_selection, "root", Nevt_BKG_after_selection, ObtainWeight("DDBAR", MCTYPE, "test", "DDBAR"));
    FillBKG(dirname_SSBAR_after_selection, "root", Nevt_BKG_after_selection, ObtainWeight("SSBAR", MCTYPE, "test", "SSBAR"));
    FillBKG(dirname_CHARM_after_selection, "root", Nevt_BKG_after_selection, ObtainWeight("CHARM", MCTYPE, "test", "CHARM"));

    FillSIGNAL(dirname_SIGNAL_after_preselection, "B2Knunu", "Bplus", Ncandidate_SIGNAL_after_preselection, Nevt_SIGNAL_after_preselection_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B2Knunu"), "B2Knunu");
    FillSIGNAL(dirname_SIGNAL_after_preselection, "B2Kstarnunu", "Bplus", Ncandidate_SIGNAL_after_preselection, Nevt_SIGNAL_after_preselection_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B2Kstarnunu"), "otherwise");
    FillSIGNAL(dirname_SIGNAL_after_preselection, "B2Xsnunu", "Bplus", Ncandidate_SIGNAL_after_preselection, Nevt_SIGNAL_after_preselection_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B2Xsnunu"), "B2Xsnunu");
    FillSIGNAL(dirname_SIGNAL_after_preselection, "B02K0nunu", "Bzero", Ncandidate_SIGNAL_after_preselection, Nevt_SIGNAL_after_preselection_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B02K0nunu"), "B02K0nunu");
    FillSIGNAL(dirname_SIGNAL_after_preselection, "B02Kstar0nunu", "Bzero", Ncandidate_SIGNAL_after_preselection, Nevt_SIGNAL_after_preselection_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B02Kstar0nunu"), "otherwise");
    FillSIGNAL(dirname_SIGNAL_after_preselection, "B02Xsnunu", "Bzero", Ncandidate_SIGNAL_after_preselection, Nevt_SIGNAL_after_preselection_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B02Xsnunu"), "B02Xsnunu");

    FillSIGNAL(dirname_SIGNAL_after_before_FBDT_cut, "B2Knunu", "Bplus", Nevt_SIGNAL_before_FBDT_cut, Nevt_SIGNAL_before_FBDT_cut_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B2Knunu"), "B2Knunu");
    FillSIGNAL(dirname_SIGNAL_after_before_FBDT_cut, "B2Kstarnunu", "Bplus", Nevt_SIGNAL_before_FBDT_cut, Nevt_SIGNAL_before_FBDT_cut_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B2Kstarnunu"), "otherwise");
    FillSIGNAL(dirname_SIGNAL_after_before_FBDT_cut, "B2Xsnunu", "Bplus", Nevt_SIGNAL_before_FBDT_cut, Nevt_SIGNAL_before_FBDT_cut_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B2Xsnunu"), "B2Xsnunu");
    FillSIGNAL(dirname_SIGNAL_after_before_FBDT_cut, "B02K0nunu", "Bzero", Nevt_SIGNAL_before_FBDT_cut, Nevt_SIGNAL_before_FBDT_cut_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B02K0nunu"), "B02K0nunu");
    FillSIGNAL(dirname_SIGNAL_after_before_FBDT_cut, "B02Kstar0nunu", "Bzero", Nevt_SIGNAL_before_FBDT_cut, Nevt_SIGNAL_before_FBDT_cut_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B02Kstar0nunu"), "otherwise");
    FillSIGNAL(dirname_SIGNAL_after_before_FBDT_cut, "B02Xsnunu", "Bzero", Nevt_SIGNAL_before_FBDT_cut, Nevt_SIGNAL_before_FBDT_cut_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B02Xsnunu"), "B02Xsnunu");
    */
    FillSIGNAL(dirname_SIGNAL_after_selection, "B2Knunu", "Bplus", Nevt_SIGNAL_after_selection, Nevt_SIGNAL_after_selection_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B2Knunu"), "B2Knunu");
    FillSIGNAL(dirname_SIGNAL_after_selection, "B2Kstarnunu", "Bplus", Nevt_SIGNAL_after_selection, Nevt_SIGNAL_after_selection_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B2Kstarnunu"), "otherwise");
    FillSIGNAL(dirname_SIGNAL_after_selection, "B2Xsnunu", "Bplus", Nevt_SIGNAL_after_selection, Nevt_SIGNAL_after_selection_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B2Xsnunu"), "B2Xsnunu");
    FillSIGNAL(dirname_SIGNAL_after_selection, "B02K0nunu", "Bzero", Nevt_SIGNAL_after_selection, Nevt_SIGNAL_after_selection_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B02K0nunu"), "B02K0nunu");
    FillSIGNAL(dirname_SIGNAL_after_selection, "B02Kstar0nunu", "Bzero", Nevt_SIGNAL_after_selection, Nevt_SIGNAL_after_selection_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B02Kstar0nunu"), "otherwise");
    FillSIGNAL(dirname_SIGNAL_after_selection, "B02Xsnunu", "Bzero", Nevt_SIGNAL_after_selection, Nevt_SIGNAL_after_selection_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B02Xsnunu"), "B02Xsnunu");
    /*
    FillSIGNAL_q2(dirname_SIGNAL_after_selection, "B2Knunu", "Bplus", Nevt_SIGNAL_after_selection_q2_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B2Knunu"), "B2Knunu");
    FillSIGNAL_q2(dirname_SIGNAL_after_selection, "B2Kstarnunu", "Bplus", Nevt_SIGNAL_after_selection_q2_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B2Kstarnunu"), "otherwise");
    FillSIGNAL_q2(dirname_SIGNAL_after_selection, "B2Xsnunu", "Bplus", Nevt_SIGNAL_after_selection_q2_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B2Xsnunu"), "B2Xsnunu");
    FillSIGNAL_q2(dirname_SIGNAL_after_selection, "B02K0nunu", "Bzero", Nevt_SIGNAL_after_selection_q2_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B02K0nunu"), "B02K0nunu");
    FillSIGNAL_q2(dirname_SIGNAL_after_selection, "B02Kstar0nunu", "Bzero", Nevt_SIGNAL_after_selection_q2_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B02Kstar0nunu"), "otherwise");
    FillSIGNAL_q2(dirname_SIGNAL_after_selection, "B02Xsnunu", "Bzero", Nevt_SIGNAL_after_selection_q2_true, ObtainWeight("SIGNAL", MCTYPE, "test", "B02Xsnunu"), "B02Xsnunu");
    */
    ReadDecayInfo(dirname_SIGNAL_initial, "B2Knunu", "Bplus", Nevt_SIGNAL_initial_true, Scale_Kplus_here, "B2Knunu");
    ReadDecayInfo(dirname_SIGNAL_initial, "B2Kstarnunu", "Bplus", Nevt_SIGNAL_initial_true, Scale_Kplusstar_here, "otherwise");
    ReadDecayInfo(dirname_SIGNAL_initial, "B2Xsnunu", "Bplus", Nevt_SIGNAL_initial_true, Scale_Xsu_nonresonant_here, "B2Xsnunu");
    ReadDecayInfo(dirname_SIGNAL_initial, "B02K0nunu", "Bzero", Nevt_SIGNAL_initial_true, Scale_K0_here, "B02K0nunu");
    ReadDecayInfo(dirname_SIGNAL_initial, "B02Kstar0nunu", "Bzero", Nevt_SIGNAL_initial_true, Scale_K0star_here, "otherwise");
    ReadDecayInfo(dirname_SIGNAL_initial, "B02Xsnunu", "Bzero", Nevt_SIGNAL_initial_true, Scale_Xsd_nonresonant_here, "B02Xsnunu");
    /*
    ReadDecayInfo_q2(dirname_SIGNAL_initial, "B2Knunu", "Bplus", Nevt_SIGNAL_initial_q2_true, Scale_Kplus_here, "B2Knunu");
    ReadDecayInfo_q2(dirname_SIGNAL_initial, "B2Kstarnunu", "Bplus", Nevt_SIGNAL_initial_q2_true, Scale_Kplusstar_here, "otherwise");
    ReadDecayInfo_q2(dirname_SIGNAL_initial, "B2Xsnunu", "Bplus", Nevt_SIGNAL_initial_q2_true, Scale_Xsu_nonresonant_here, "B2Xsnunu");
    ReadDecayInfo_q2(dirname_SIGNAL_initial, "B02K0nunu", "Bzero", Nevt_SIGNAL_initial_q2_true, Scale_K0_here, "B02K0nunu");
    ReadDecayInfo_q2(dirname_SIGNAL_initial, "B02Kstar0nunu", "Bzero", Nevt_SIGNAL_initial_q2_true, Scale_K0star_here, "otherwise");
    ReadDecayInfo_q2(dirname_SIGNAL_initial, "B02Xsnunu", "Bzero", Nevt_SIGNAL_initial_q2_true, Scale_Xsd_nonresonant_here, "B02Xsnunu");
    */
    SIGNAL_pre_and_selection_efficiency_true->Divide(Nevt_SIGNAL_after_selection_true, Nevt_SIGNAL_initial_true, 1.0, 1.0, "B");
    /*
    SIGNAL_selection_efficiency->Divide(Nevt_SIGNAL_after_selection, Ncandidate_SIGNAL_after_preselection, 1.0, 1.0, "B");
    SIGNAL_FBDT_efficiency->Divide(Nevt_SIGNAL_after_selection, Nevt_SIGNAL_before_FBDT_cut, 1.0, 1.0, "B");

    BKG_selection_efficiency->Divide(Nevt_BKG_after_selection, Ncandidate_BKG_after_preselection, 1.0, 1.0, "B");
    BKG_FBDT_efficiency->Divide(Nevt_BKG_after_selection, Nevt_BKG_before_FBDT_cut, 1.0, 1.0, "B");

    SIGNAL_FBDT_efficiency_true->Divide(Nevt_SIGNAL_after_selection_true, Nevt_SIGNAL_before_FBDT_cut_true, 1.0, 1.0, "B");
    SIGNAL_FBDT_efficiency_true->SetStats(false);

    SIGNAL_efficiency_q2_true->Divide(Nevt_SIGNAL_after_selection_q2_true, Nevt_SIGNAL_initial_q2_true, 1.0, 1.0, "B");
    SIGNAL_efficiency_q2_true->SetStats(false);
    */
    // define TCanvas
    TCanvas* c_temp = new TCanvas("c", "", 1200, 1000);

    // set draw option
    SIGNAL_pre_and_selection_efficiency_true->SetMarkerSize(0.0);
    SIGNAL_pre_and_selection_efficiency_true->SetLineWidth(2);
    SIGNAL_pre_and_selection_efficiency_true->SetMarkerColor(4);
    SIGNAL_pre_and_selection_efficiency_true->SetMarkerStyle(21);
    SIGNAL_pre_and_selection_efficiency_true->SetStats(false);
    /*
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
    */
    TPaveText* pt;

    // draw SIGNAL preselection + selection efficiency true
    SIGNAL_pre_and_selection_efficiency_true->Draw("P");
    pt = new TPaveText(0.135, 0.88, 0.5, 1.0, "NDC NB"); pt->SetFillStyle(0); pt->SetLineWidth(0); pt->AddText("efficiency for preselection + selection"); pt->Draw();
    c_temp->SaveAs(("Efficiency_index_" + std::to_string(INDEX) + ".png").c_str());
    /*
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
    */

    // calculate entire eff
    SIGNAL_pre_and_selection_efficiency_true_all->Divide(Nevt_SIGNAL_after_selection_true_all, Nevt_SIGNAL_initial_true_all, 1.0, 1.0, "B");
    printf("eff: %lf +- %lf\n", SIGNAL_pre_and_selection_efficiency_true_all->GetBinContent(1), SIGNAL_pre_and_selection_efficiency_true_all->GetBinError(1));

    return 0;
}
