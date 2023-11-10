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
using namespace RooFit;
using std::string;
using std::to_string;

# define N_decay 38 // five decay mode + others

# define MyEPSILON 0.000001
# define MCTYPE "MC15ri"

// arXiv:1409.4557v2, PhysRevD.107.014511
# define TB0 1.5195 // (Table. 1)
# define TBp 1.6384 // (Table. 1)
# define BR_Kplus_nunubar 0.000005044 // Table VI = (5.044 +- 0.402) * 10^{-6}
# define BR_K0star_nunubar 0.00000919 // (eq. 11)
# define BR_K0_nunubar (BR_Kplus_nunubar*TB0/TBp) // under (eq. 15). In Table VI = (4.6669 +- 0.3707) * 10^{-6}
# define BR_Kplusstar_nunubar (BR_K0star_nunubar*TBp/TB0) // under (eq. 15)
# define BR_Xs_nunubar 0.000029 // (eq. 23)
# define BR_Xsu_nonresonant_nunubar (BR_Xs_nunubar - BR_Kplus_nunubar - BR_Kplusstar_nunubar)
# define BR_Xsd_nonresonant_nunubar (BR_Xs_nunubar - BR_K0_nunubar - BR_K0star_nunubar)

// according to DIRAC
# define N_BpBp_1invab 540000000.0
# define N_B0B0_1invab 510000000.0

# define BR_BpBp 0.514
# define BR_B0B0 0.486

// https://confluence.desy.de/pages/viewpage.action?spaceKey=BI&title=Conference+readiness
# define N_BB_LS1 387100000.0 // NBB = (387.1 +/- 5.6) x 10^6

# define N_Kplus_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_BpBp/(BR_BpBp+BR_B0B0)) * BR_Kplus_nunubar)
# define N_Kplusstar_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_BpBp/(BR_BpBp+BR_B0B0)) * BR_Kplusstar_nunubar)
# define N_Xsu_nonresonant_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_BpBp/(BR_BpBp+BR_B0B0)) * BR_Xsu_nonresonant_nunubar)
# define N_K0_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_B0B0/(BR_BpBp+BR_B0B0)) * BR_K0_nunubar)
# define N_K0star_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_B0B0/(BR_BpBp+BR_B0B0)) * BR_K0star_nunubar)
# define N_Xsd_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_B0B0/(BR_BpBp+BR_B0B0)) * BR_Xsd_nonresonant_nunubar)

// SIGNAL MC sample number before skimming
# define N_Kplus_train 7039000.0
# define N_K0_train 7166624.0
# define N_Kplusstar_train 7039000.0
# define N_K0star_train 7166624.0
# define N_Xsu_nonresonant_train 35195000.0
# define N_Xsd_nonresonant_train 34940430.0
# define N_Kplus_test 2961000.0
# define N_K0_test 2833376.0
# define N_Kplusstar_test 2961000.0
# define N_K0star_test 2833376.0
# define N_Xsu_nonresonant_test 14805000.0
# define N_Xsd_nonresonant_test 15059570.0

// scale factor for SIGNAL MC sample until LS1
# define Scale_Kplus_train (N_Kplus_nunubar_LS1/N_Kplus_train)
# define Scale_Kplusstar_train (N_Kplusstar_nunubar_LS1/N_Kplusstar_train)
# define Scale_Xsu_nonresonant_train (N_Xsu_nonresonant_nunubar_LS1/N_Xsu_nonresonant_train)
# define Scale_K0_train (N_K0_nunubar_LS1/N_K0_train)
# define Scale_K0star_train (N_K0star_nunubar_LS1/N_K0star_train)
# define Scale_Xsd_nonresonant_train (N_Xsd_nunubar_LS1/N_Xsd_nonresonant_train)
# define Scale_Kplus_test (N_Kplus_nunubar_LS1/N_Kplus_test)
# define Scale_Kplusstar_test (N_Kplusstar_nunubar_LS1/N_Kplusstar_test)
# define Scale_Xsu_nonresonant_test (N_Xsu_nonresonant_nunubar_LS1/N_Xsu_nonresonant_test)
# define Scale_K0_test (N_K0_nunubar_LS1/N_K0_test)
# define Scale_K0star_test (N_K0star_nunubar_LS1/N_K0star_test)
# define Scale_Xsd_nonresonant_test (N_Xsd_nunubar_LS1/N_Xsd_nonresonant_test)

// BKG MC sample number (2.8/ab for BB, 1.0/ab for qq)
# define N_CHG_test 139768443.0
# define N_MIX_test 104591690.0
# define N_UUBAR_test 107641351.0
# define N_DDBAR_test 25588403.0
# define N_SSBAR_test 20668794.0
# define N_CHARM_test 116108850.0
# define N_CHG_train 151915459.0
# define N_MIX_train 120136353.0
# define N_UUBAR_train 158509639.0
# define N_DDBAR_train 38644413.0
# define N_SSBAR_train 31001866.0
# define N_CHARM_train 174901296.0

// new scale factor for BKG MC sample with additional 1/ab (364.436 - 2.763 = 361.673/fb), until LS1
# define Scale_CHG_train ((N_BB_LS1* (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab * (N_CHG_train / (N_CHG_train + N_CHG_test)) ))
# define Scale_MIX_train ((N_BB_LS1* (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab * (N_MIX_train / (N_MIX_train + N_MIX_test)) ))
# define Scale_UUBAR_train (0.361673/((N_UUBAR_train/(N_UUBAR_train + N_UUBAR_test))*1.0))
# define Scale_DDBAR_train (0.361673/((N_DDBAR_train/(N_DDBAR_train + N_DDBAR_test))*1.0))
# define Scale_SSBAR_train (0.361673/((N_SSBAR_train/(N_SSBAR_train + N_SSBAR_test))*1.0))
# define Scale_CHARM_train (0.361673/((N_CHARM_train/(N_CHARM_train + N_CHARM_test))*1.0))
# define Scale_CHG_test ((N_BB_LS1* (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab * (N_CHG_test / (N_CHG_train + N_CHG_test)) ))
# define Scale_MIX_test ((N_BB_LS1* (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab * (N_MIX_test / (N_MIX_train + N_MIX_test)) ))
# define Scale_UUBAR_test (0.361673/((N_UUBAR_test/(N_UUBAR_train + N_UUBAR_test))*1.0))
# define Scale_DDBAR_test (0.361673/((N_DDBAR_test/(N_DDBAR_train + N_DDBAR_test))*1.0))
# define Scale_SSBAR_test (0.361673/((N_SSBAR_test/(N_SSBAR_train + N_SSBAR_test))*1.0))
# define Scale_CHARM_test (0.361673/((N_CHARM_test/(N_CHARM_train + N_CHARM_test))*1.0))

// scale facto for systematic MC sample
# define N_K0star_nunubar_syst 10000000.0
# define N_Kplusstar_nunubar_syst 10000000.0
# define N_Xsu_nonresonant_nunubar_syst 10000000.0
# define N_Xsd_nonresonant_nunubar_syst 10000000.0

// scale factor for each systematic MC sample (364.436 - 2.763 = 361.673/fb)
# define Scale_Kplusstar_syst (N_Kplusstar_nunubar_LS1/N_Kplusstar_nunubar_syst)
# define Scale_Xsu_nonresonant_syst (N_Xsu_nonresonant_nunubar_LS1/N_Xsu_nonresonant_nunubar_syst)
# define Scale_K0star_syst (N_K0star_nunubar_LS1/N_K0star_nunubar_syst)
# define Scale_Xsd_nonresonant_syst (N_Xsd_nunubar_LS1/N_Xsd_nonresonant_nunubar_syst)

# define KS0_rel_uncertainty 0.6 // %/cm
# define track_rel_uncertainty 0.69 // %
# define Kaon_PID_max_uncertainty 0.1 // not percentage. relative uncertainty

# define FEI_cal_Bc_num 12
# define FEI_cal_B0_num 11
class Corrector_FEI {
private:
    static constexpr double FEI_cal_Bc_MC15ri[FEI_cal_Bc_num] = { 1.04, 0.79, 0.69, 0.56, 0.97, 0.95, 0.74, 0.57, 0.91, 0.51, 0.34, 0.59 };
    static constexpr double FEI_cal_Bc_uncertainty_MC15ri[FEI_cal_Bc_num] = { 0.03, 0.03, 0.05, 0.11, 0.03, 0.03, 0.02, 0.06, 0.1, 0.13, 0.07, 0.02 }; // not relative uncertainty. absolute uncertainty
    static constexpr double FEI_cal_Bc_modeID_MC15ri[FEI_cal_Bc_num] = { 0.0, 1.0, 3.0, 4.0, 15.0, 16.0, 18.0, 19.0, 23.0, 24.0, 30.0, -1.0 };
    static constexpr double FEI_cal_B0_MC15ri[FEI_cal_B0_num] = { 1.16, 0.94, 0.81, 0.79, 0.99, 1.03, 0.67, 0.66, 0.69, 0.49, 0.79 };
    static constexpr double FEI_cal_B0_uncertainty_MC15ri[FEI_cal_B0_num] = { 0.04, 0.05, 0.06, 0.02, 0.03, 0.06, 0.02, 0.03, 0.02, 0.02, 0.12 }; // not relative uncertainty. absolute uncertainty
    static constexpr double FEI_cal_B0_modeID_MC15ri[FEI_cal_B0_num] = { 0.0, 1.0, 3.0, 4.0, 5.0, 15.0, 16.0, 18.0, 19.0, 26.0, -1.0 };
public:
    Corrector_FEI();
    double GetFEICalFactor(double UpsilonID, double BtagID, std::string type);
    double GetFEICalFactor(int index, bool IsItCharged, std::string type);
    double GetFEICalFactorUncer(double UpsilonID, double BtagID, std::string type);
    double GetmodeID(int index, bool IsItCharged, std::string type);
};

Corrector_FEI corrector_FEI;

Corrector_FEI::Corrector_FEI() {}

double Corrector_FEI::GetFEICalFactor(double UpsilonID, double BtagID, std::string type) {
    // UpsilonID => charged: 0, mixed: 1

    if (type == "MC15ri") {
        if (UpsilonID > -0.5 && UpsilonID < 0.5) { // charged
            for (int i = 0; i < FEI_cal_Bc_num - 1; i++) {
                if (BtagID > FEI_cal_Bc_modeID_MC15ri[i] - 0.5 && BtagID < FEI_cal_Bc_modeID_MC15ri[i] + 0.5) return FEI_cal_Bc_MC15ri[i];
            }
            return FEI_cal_Bc_MC15ri[FEI_cal_Bc_num - 1];
        }
        else if (UpsilonID > 0.5 && UpsilonID < 1.5) { // mixed
            for (int i = 0; i < FEI_cal_B0_num - 1; i++) {
                if (BtagID > FEI_cal_B0_modeID_MC15ri[i] - 0.5 && BtagID < FEI_cal_B0_modeID_MC15ri[i] + 0.5) return FEI_cal_B0_MC15ri[i];
            }
            return FEI_cal_B0_MC15ri[FEI_cal_B0_num - 1];
        }
    }
    else {
        printf("[Corrector_FEI] Invalid type!\n");
        exit(1);
    }

    printf("[Corrector_FEI] error! unexpected decay ID\n");
    exit(1);
    return 0;

}

double Corrector_FEI::GetFEICalFactor(int index, bool IsItCharged, std::string type) {
    if (type == "MC15ri") {
        if (IsItCharged) return FEI_cal_Bc_MC15ri[index];
        else FEI_cal_B0_MC15ri[index];
    }
    else {
        printf("[Corrector_FEI] Invalid type!\n");
        exit(1);
        return 0;
    }
}

double Corrector_FEI::GetFEICalFactorUncer(double UpsilonID, double BtagID, std::string type) {
    // UpsilonID => charged: 0, mixed: 1

    if (type == "MC15ri") {
        if (UpsilonID > -0.5 && UpsilonID < 0.5) { // charged
            for (int i = 0; i < FEI_cal_Bc_num - 1; i++) {
                if (BtagID > FEI_cal_Bc_modeID_MC15ri[i] - 0.5 && BtagID < FEI_cal_Bc_modeID_MC15ri[i] + 0.5) return FEI_cal_Bc_uncertainty_MC15ri[i];
            }
            return FEI_cal_Bc_uncertainty_MC15ri[FEI_cal_Bc_num - 1];
        }
        else if (UpsilonID > 0.5 && UpsilonID < 1.5) { // mixed
            for (int i = 0; i < FEI_cal_B0_num - 1; i++) {
                if (BtagID > FEI_cal_B0_modeID_MC15ri[i] - 0.5 && BtagID < FEI_cal_B0_modeID_MC15ri[i] + 0.5) return FEI_cal_B0_uncertainty_MC15ri[i];
            }
            return FEI_cal_B0_uncertainty_MC15ri[FEI_cal_B0_num - 1];
        }
    }
    else {
        printf("[Corrector_FEI] Invalid type!\n");
        exit(1);
    }

    printf("[Corrector_FEI] error! unexpected decay ID\n");
    exit(1);
    return 0;

}

double Corrector_FEI::GetmodeID(int index, bool IsItCharged, std::string type) {
    if (type == "MC15ri") {
        if (IsItCharged) return FEI_cal_Bc_modeID_MC15ri[index];
        else FEI_cal_B0_modeID_MC15ri[index];
    }
    else {
        printf("[Corrector_FEI] Invalid type!\n");
        exit(1);
        return 0;
    }
}

# define Xsu_frag_decay1 531.4
# define Xsu_frag_decay2 1062.1
# define Xsu_frag_decay3 673.9
# define Xsu_frag_decay4 825.3
# define Xsu_frag_decay5 473.3
# define Xsu_frag_decay6 119.1
# define Xsu_frag_decay7 24.6
# define Xsu_frag_decay8 148.2
# define Xsu_frag_decay9 127.6
# define Xsu_frag_decay10 102.9
# define Xsu_frag_decay11 91.8
# define Xsu_frag_decay12 54.8
# define Xsu_frag_decay13 21.7
# define Xsu_frag_decay14 11.0
# define Xsu_frag_decay30 858.4

# define Xsd_frag_decay16 1072.6
# define Xsd_frag_decay17 536.2
# define Xsd_frag_decay18 835.5
# define Xsd_frag_decay19 677.4
# define Xsd_frag_decay20 120.5
# define Xsd_frag_decay21 475.9
# define Xsd_frag_decay22 149.8
# define Xsd_frag_decay23 25.0
# define Xsd_frag_decay24 129.2
# define Xsd_frag_decay25 104.4
# define Xsd_frag_decay26 92.4
# define Xsd_frag_decay27 54.9
# define Xsd_frag_decay28 14.5
# define Xsd_frag_decay29 10.9
# define Xsd_frag_decay30 871.0

// uncertainty of BR from parametric reason
# define Sigma_BR_Kplus_nunubar 0.000000402
# define Sigma_BR_K0star_nunubar 0.00000050
# define Sigma_BR_K0_nunubar (Sigma_BR_Kplus_nunubar*TB0/TBp)
# define Sigma_BR_Kplusstar_nunubar (Sigma_BR_K0star_nunubar*TBp/TB0)
# define Sigma_BR_Xs_nunubar 0.000003

# define N_PID_syst 73
class Corrector_PID {
private:
    double PID_correction_MC15ri[4][N_PID_syst]; // K-true, K-mis, pi-true, pi-miss
    double PID_correction_stat_uncer_MC15ri[4][N_PID_syst]; // K-true, K-mis, pi-true, pi-miss
    double PID_correction_sys_uncer_MC15ri[4][N_PID_syst]; // K-true, K-mis, pi-true, pi-miss
    double PID_correction_uncer_MC15ri[4][N_PID_syst]; // K-true, K-mis, pi-true, pi-miss
    void ReadPIDFile_MC15ri();
public:
    Corrector_PID();
    double GetCorrectionFactor(int PID_type, int bin_PID, std::string type);
    double GetUncertainty(int PID_type, int bin_PID, std::string type);
};

Corrector_PID corrector_PID;

Corrector_PID::Corrector_PID() {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < N_PID_syst; j++) {
            PID_correction_MC15ri[i][j] = 0.0;
            PID_correction_stat_uncer_MC15ri[i][j] = 0.0;
            PID_correction_sys_uncer_MC15ri[i][j] = 0.0;
            PID_correction_uncer_MC15ri[i][j] = 0.0;
        }
    }
    ReadPIDFile_MC15ri();
}

void Corrector_PID::ReadPIDFile_MC15ri() {
    const char* KID_true_file = "/home/belle2/junewoo/storage_b1/bsub/systematic/MC15ri_PID/KaonEff.csv";
    const char* KID_mis_file = "/home/belle2/junewoo/storage_b1/bsub/systematic/MC15ri_PID/Kaonmis.csv";
    const char* PID_true_file = "/home/belle2/junewoo/storage_b1/bsub/systematic/MC15ri_PID/PionEff.csv";
    const char* PID_mis_file = "/home/belle2/junewoo/storage_b1/bsub/systematic/MC15ri_PID/Pionmis.csv";

    FILE* fp_KID_true = fopen(KID_true_file, "r");
    FILE* fp_KID_mis = fopen(KID_mis_file, "r");
    FILE* fp_PID_true = fopen(PID_true_file, "r");
    FILE* fp_PID_mis = fopen(PID_mis_file, "r");

    fscanf(fp_KID_true, "p_min,p_max,cosTheta_min,cosTheta_max,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,threshold,variable\n");
    fscanf(fp_KID_mis, "p_min,p_max,cosTheta_min,cosTheta_max,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,threshold,variable\n");
    fscanf(fp_PID_true, "p_min,p_max,cosTheta_min,cosTheta_max,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,threshold,variable\n");
    fscanf(fp_PID_mis, "p_min,p_max,cosTheta_min,cosTheta_max,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,threshold,variable\n");

    double temp_p_min;
    double temp_p_max;
    double temp_cosTheta_min;
    double temp_cosTheta_max;
    double temp_data_MC_ratio;
    double temp_data_MC_uncertainty_stat_up;
    double temp_data_MC_uncertainty_stat_dn;
    double temp_data_MC_uncertainty_sys_up;
    double temp_data_MC_uncertainty_sys_dn;
    double temp_data_efficiency;
    double temp_data_uncertainty_stat_up;
    double temp_data_uncertainty_stat_dn;
    double temp_data_uncertainty_sys_up;
    double temp_data_uncertainty_sys_dn;
    double temp_MC_efficiency;
    double temp_MC_uncertainty_stat_up;
    double temp_MC_uncertainty_stat_dn;
    double temp_MC_uncertainty_sys_up;
    double temp_MC_uncertainty_sys_dn;
    double temp_threshold;

    for (int i = 0; i < N_PID_syst - 1; i++) {
        fscanf(fp_KID_true, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,kaonID\n", &temp_p_min, &temp_p_max, &temp_cosTheta_min, &temp_cosTheta_max, &PID_correction_MC15ri[0][i], &PID_correction_stat_uncer_MC15ri[0][i], &temp_data_MC_uncertainty_stat_dn, &PID_correction_sys_uncer_MC15ri[0][i], &temp_data_MC_uncertainty_sys_dn, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_threshold);
        fscanf(fp_KID_mis, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,kaonID\n", &temp_p_min, &temp_p_max, &temp_cosTheta_min, &temp_cosTheta_max, &PID_correction_MC15ri[1][i], &PID_correction_stat_uncer_MC15ri[1][i], &temp_data_MC_uncertainty_stat_dn, &PID_correction_sys_uncer_MC15ri[1][i], &temp_data_MC_uncertainty_sys_dn, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_threshold);
        fscanf(fp_PID_true, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,pionID\n", &temp_p_min, &temp_p_max, &temp_cosTheta_min, &temp_cosTheta_max, &PID_correction_MC15ri[2][i], &PID_correction_stat_uncer_MC15ri[2][i], &temp_data_MC_uncertainty_stat_dn, &PID_correction_sys_uncer_MC15ri[2][i], &temp_data_MC_uncertainty_sys_dn, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_threshold);
        fscanf(fp_PID_mis, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,pionID\n", &temp_p_min, &temp_p_max, &temp_cosTheta_min, &temp_cosTheta_max, &PID_correction_MC15ri[3][i], &PID_correction_stat_uncer_MC15ri[3][i], &temp_data_MC_uncertainty_stat_dn, &PID_correction_sys_uncer_MC15ri[3][i], &temp_data_MC_uncertainty_sys_dn, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_threshold);

        if ((std::abs(PID_correction_MC15ri[0][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer_MC15ri[0][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer_MC15ri[0][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction_MC15ri[0][i]) > 10000.0)) {
            PID_correction_MC15ri[0][i] = 1.0;
            PID_correction_stat_uncer_MC15ri[0][i] = 0.0;
            PID_correction_sys_uncer_MC15ri[0][i] = 0.0;
        }
        if ((std::abs(PID_correction_MC15ri[1][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer_MC15ri[1][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer_MC15ri[1][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction_MC15ri[1][i]) > 10000.0)) {
            PID_correction_MC15ri[1][i] = 1.0;
            PID_correction_stat_uncer_MC15ri[1][i] = 0.0;
            PID_correction_sys_uncer_MC15ri[1][i] = 0.0;
        }
        if ((std::abs(PID_correction_MC15ri[2][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer_MC15ri[2][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer_MC15ri[2][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction_MC15ri[2][i]) > 10000.0)) {
            PID_correction_MC15ri[2][i] = 1.0;
            PID_correction_stat_uncer_MC15ri[2][i] = 0.0;
            PID_correction_sys_uncer_MC15ri[2][i] = 0.0;
        }
        if ((std::abs(PID_correction_MC15ri[3][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer_MC15ri[3][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer_MC15ri[3][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction_MC15ri[3][i]) > 10000.0)) {
            PID_correction_MC15ri[3][i] = 1.0;
            PID_correction_stat_uncer_MC15ri[3][i] = 0.0;
            PID_correction_sys_uncer_MC15ri[3][i] = 0.0;
        }

        PID_correction_uncer_MC15ri[0][i] = std::sqrt(PID_correction_stat_uncer_MC15ri[0][i] * PID_correction_stat_uncer_MC15ri[0][i] + PID_correction_sys_uncer_MC15ri[0][i] * PID_correction_sys_uncer_MC15ri[0][i]);
        PID_correction_uncer_MC15ri[1][i] = std::sqrt(PID_correction_stat_uncer_MC15ri[1][i] * PID_correction_stat_uncer_MC15ri[1][i] + PID_correction_sys_uncer_MC15ri[1][i] * PID_correction_sys_uncer_MC15ri[1][i]);
        PID_correction_uncer_MC15ri[2][i] = std::sqrt(PID_correction_stat_uncer_MC15ri[2][i] * PID_correction_stat_uncer_MC15ri[2][i] + PID_correction_sys_uncer_MC15ri[2][i] * PID_correction_sys_uncer_MC15ri[2][i]);
        PID_correction_uncer_MC15ri[3][i] = std::sqrt(PID_correction_stat_uncer_MC15ri[3][i] * PID_correction_stat_uncer_MC15ri[3][i] + PID_correction_sys_uncer_MC15ri[3][i] * PID_correction_sys_uncer_MC15ri[3][i]);
    }

    PID_correction_MC15ri[0][N_PID_syst - 1] = 1.0;
    PID_correction_stat_uncer_MC15ri[0][N_PID_syst - 1] = 0.0;
    PID_correction_sys_uncer_MC15ri[0][N_PID_syst - 1] = 0.0;
    PID_correction_uncer_MC15ri[0][N_PID_syst - 1] = 0.0;

    PID_correction_MC15ri[1][N_PID_syst - 1] = 1.0;
    PID_correction_stat_uncer_MC15ri[1][N_PID_syst - 1] = 0.0;
    PID_correction_sys_uncer_MC15ri[1][N_PID_syst - 1] = 0.0;
    PID_correction_uncer_MC15ri[1][N_PID_syst - 1] = 0.0;

    PID_correction_MC15ri[2][N_PID_syst - 1] = 1.0;
    PID_correction_stat_uncer_MC15ri[2][N_PID_syst - 1] = 0.0;
    PID_correction_sys_uncer_MC15ri[2][N_PID_syst - 1] = 0.0;
    PID_correction_uncer_MC15ri[2][N_PID_syst - 1] = 0.0;

    PID_correction_MC15ri[3][N_PID_syst - 1] = 1.0;
    PID_correction_stat_uncer_MC15ri[3][N_PID_syst - 1] = 0.0;
    PID_correction_sys_uncer_MC15ri[3][N_PID_syst - 1] = 0.0;
    PID_correction_uncer_MC15ri[3][N_PID_syst - 1] = 0.0;

    fclose(fp_KID_true);
    fclose(fp_KID_mis);
    fclose(fp_PID_true);
    fclose(fp_PID_mis);
}

double Corrector_PID::GetCorrectionFactor(int PID_type, int bin_PID, std::string type) {
    if (type == "MC15ri") return PID_correction_MC15ri[PID_type][bin_PID];
    else {
        printf("[Corrector_pi0] Invalid type!\n");
        exit(1);
    }
}

double Corrector_PID::GetUncertainty(int PID_type, int bin_PID, std::string type) {
    if (type == "MC15ri") return PID_correction_uncer_MC15ri[PID_type][bin_PID];
    else {
        printf("[Corrector_pi0] Invalid type!\n");
        exit(1);
    }
}

# define N_fakeE_syst 37
# define N_fakeMU_syst 49
class Corrector_FakePID {
private:
    double PID_fakeE_correction_MC15ri[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double PID_fakeE_uncer_MC15ri[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double PID_fakeMU_correction_MC15ri[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double PID_fakeMU_uncer_MC15ri[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    void ReadPIDFile_MC15ri();
public:
    Corrector_FakePID();
    double GetCorrectionFactorfakeE(int PID_type, int bin_PID, std::string type);
    double GetCorrectionFactorfakeMU(int PID_type, int bin_PID, std::string type);
    double GetUncertaintyfakeE(int PID_type, int bin_PID, std::string type);
    double GetUncertaintyfakeMU(int PID_type, int bin_PID, std::string type);
};

Corrector_FakePID corrector_FakePID;

Corrector_FakePID::Corrector_FakePID() {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < N_fakeE_syst; j++) {
            PID_fakeE_correction_MC15ri[i][j] = 0.0;
            PID_fakeE_uncer_MC15ri[i][j] = 0.0;
        }
    }
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < N_fakeMU_syst; j++) {
            PID_fakeMU_correction_MC15ri[i][j] = 0.0;
            PID_fakeMU_uncer_MC15ri[i][j] = 0.0;
        }
    }
    ReadPIDFile_MC15ri();
}

void Corrector_FakePID::ReadPIDFile_MC15ri() {
    // initialization
    for (int i = 0; i < N_fakeE_syst; i++) {
        PID_fakeE_correction_MC15ri[0][i] = 1.0; //  K-, K+, pi-, pi+
        PID_fakeE_correction_MC15ri[1][i] = 1.0;
        PID_fakeE_correction_MC15ri[2][i] = 1.0;
        PID_fakeE_correction_MC15ri[3][i] = 1.0;

        PID_fakeE_uncer_MC15ri[0][i] = 0.0;
        PID_fakeE_uncer_MC15ri[1][i] = 0.0;
        PID_fakeE_uncer_MC15ri[2][i] = 0.0;
        PID_fakeE_uncer_MC15ri[3][i] = 0.0;
    }

    for (int i = 0; i < N_fakeMU_syst; i++) {
        PID_fakeMU_correction_MC15ri[0][i] = 1.0;
        PID_fakeMU_correction_MC15ri[1][i] = 1.0;
        PID_fakeMU_correction_MC15ri[2][i] = 1.0;
        PID_fakeMU_correction_MC15ri[3][i] = 1.0;

        PID_fakeMU_uncer_MC15ri[0][i] = 0.0;
        PID_fakeMU_uncer_MC15ri[1][i] = 0.0;
        PID_fakeMU_uncer_MC15ri[2][i] = 0.0;
        PID_fakeMU_uncer_MC15ri[3][i] = 0.0;
    }

    const char* K_fromE_file = "/home/belle2/junewoo/storage_b1/bsub/systematic/MC15ri_PID/kaonID_efficiency_electron.csv";
    const char* K_fromMU_file = "/home/belle2/junewoo/storage_b1/bsub/systematic/MC15ri_PID/kaonID_efficiency_muon.csv";
    const char* pi_fromE_file = "/home/belle2/junewoo/storage_b1/bsub/systematic/MC15ri_PID/pionID_efficiency_electron.csv";
    const char* pi_fromMU_file = "/home/belle2/junewoo/storage_b1/bsub/systematic/MC15ri_PID/pionID_efficiency_muon.csv";

    FILE* fp_K_fromE = fopen(K_fromE_file, "r");
    FILE* fp_K_fromMU = fopen(K_fromMU_file, "r");
    FILE* fp_pi_fromE = fopen(pi_fromE_file, "r");
    FILE* fp_pi_fromMU = fopen(pi_fromMU_file, "r");

    fscanf(fp_K_fromE, "variable,charge,p_min,p_max,theta_min,theta_max,iso_score_min,iso_score_max,working_point,threshold,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn\n");
    fscanf(fp_K_fromMU, "variable,charge,p_min,p_max,theta_min,theta_max,iso_score_min,iso_score_max,working_point,threshold,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn\n");
    fscanf(fp_pi_fromE, "variable,charge,p_min,p_max,theta_min,theta_max,iso_score_min,iso_score_max,working_point,threshold,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn\n");
    fscanf(fp_pi_fromMU, "variable,charge,p_min,p_max,theta_min,theta_max,iso_score_min,iso_score_max,working_point,threshold,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn\n");

    char temp_charge;
    double temp_p_min;
    double temp_p_max;
    double temp_theta_min;
    double temp_theta_max;
    double temp_data_MC_ratio;
    double temp_data_MC_uncertainty_stat_up;
    double temp_data_MC_uncertainty_stat_dn;
    double temp_data_MC_uncertainty_sys_up;
    double temp_data_MC_uncertainty_sys_dn;
    double temp_data_efficiency;
    double temp_data_uncertainty_stat_up;
    double temp_data_uncertainty_stat_dn;
    double temp_data_uncertainty_sys_up;
    double temp_data_uncertainty_sys_dn;
    double temp_MC_efficiency;
    double temp_MC_uncertainty_stat_up;
    double temp_MC_uncertainty_stat_dn;
    double temp_MC_uncertainty_sys_up;
    double temp_MC_uncertainty_sys_dn;
    double temp_threshold;

    // Kaon from fake electron
    while (fscanf(fp_K_fromE, "kaonID,%c,%lf,%lf,%lf,%lf,0.0,1.0,FixedThresh06,0.6,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &temp_charge, &temp_p_min, &temp_p_max, &temp_theta_min, &temp_theta_max, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_data_MC_ratio, &temp_data_MC_uncertainty_stat_up, &temp_data_MC_uncertainty_stat_dn, &temp_data_MC_uncertainty_sys_up, &temp_data_MC_uncertainty_sys_dn) != EOF) {
        int p_bin = -1;
        int theta_bin = -1;

        if (std::abs(temp_p_min - 0.4) < MyEPSILON && std::abs(temp_p_max - 0.5) < MyEPSILON) p_bin = 0;
        else if (std::abs(temp_p_min - 0.5) < MyEPSILON && std::abs(temp_p_max - 1.0) < MyEPSILON) p_bin = 1;
        else if (std::abs(temp_p_min - 1.0) < MyEPSILON && std::abs(temp_p_max - 1.5) < MyEPSILON) p_bin = 2;
        else if (std::abs(temp_p_min - 1.5) < MyEPSILON && std::abs(temp_p_max - 2.0) < MyEPSILON) p_bin = 3;
        else if (std::abs(temp_p_min - 2.0) < MyEPSILON && std::abs(temp_p_max - 2.5) < MyEPSILON) p_bin = 4;
        else if (std::abs(temp_p_min - 2.5) < MyEPSILON && std::abs(temp_p_max - 3.0) < MyEPSILON) p_bin = 5;
        else {
            printf("[ERROR] unknown p bin!\n");
            exit(1);
        }

        if (std::abs(temp_theta_min - 0.22) < MyEPSILON && std::abs(temp_theta_max - 0.56) < MyEPSILON) theta_bin = 0;
        else if (std::abs(temp_theta_min - 0.56) < MyEPSILON && std::abs(temp_theta_max - 1.13) < MyEPSILON) theta_bin = 1;
        else if (std::abs(temp_theta_min - 1.13) < MyEPSILON && std::abs(temp_theta_max - 1.57) < MyEPSILON) theta_bin = 2;
        else if (std::abs(temp_theta_min - 1.57) < MyEPSILON && std::abs(temp_theta_max - 1.88) < MyEPSILON) theta_bin = 3;
        else if (std::abs(temp_theta_min - 1.88) < MyEPSILON && std::abs(temp_theta_max - 2.23) < MyEPSILON) theta_bin = 4;
        else if (std::abs(temp_theta_min - 2.23) < MyEPSILON && std::abs(temp_theta_max - 2.71) < MyEPSILON) theta_bin = 5;
        else {
            printf("[ERROR] unknown theta bin!\n");
            exit(1);
        }

        int bin = theta_bin + 6 * p_bin;

        if (temp_charge == '+') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeE_correction_MC15ri[1][bin] = temp_data_MC_ratio;
            PID_fakeE_uncer_MC15ri[1][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeE_correction_MC15ri[1][bin]) < MyEPSILON) {
                PID_fakeE_correction_MC15ri[1][bin] = 1.0;
                PID_fakeE_uncer_MC15ri[1][bin] = 0.0;
            }

        }
        else if (temp_charge == '-') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeE_correction_MC15ri[0][bin] = temp_data_MC_ratio;
            PID_fakeE_uncer_MC15ri[0][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeE_correction_MC15ri[0][bin]) < MyEPSILON) {
                PID_fakeE_correction_MC15ri[0][bin] = 1.0;
                PID_fakeE_uncer_MC15ri[0][bin] = 0.0;
            }

        }
        else {
            printf("[ERROR] unknown charge!\n");
            exit(1);
        }
    }
    fclose(fp_K_fromE);

    // Pion from fake electron
    while (fscanf(fp_pi_fromE, "pionID,%c,%lf,%lf,%lf,%lf,0.0,1.0,FixedThresh06,0.6,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &temp_charge, &temp_p_min, &temp_p_max, &temp_theta_min, &temp_theta_max, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_data_MC_ratio, &temp_data_MC_uncertainty_stat_up, &temp_data_MC_uncertainty_stat_dn, &temp_data_MC_uncertainty_sys_up, &temp_data_MC_uncertainty_sys_dn) != EOF) {
        int p_bin = -1;
        int theta_bin = -1;

        if (std::abs(temp_p_min - 0.4) < MyEPSILON && std::abs(temp_p_max - 0.5) < MyEPSILON) p_bin = 0;
        else if (std::abs(temp_p_min - 0.5) < MyEPSILON && std::abs(temp_p_max - 1.0) < MyEPSILON) p_bin = 1;
        else if (std::abs(temp_p_min - 1.0) < MyEPSILON && std::abs(temp_p_max - 1.5) < MyEPSILON) p_bin = 2;
        else if (std::abs(temp_p_min - 1.5) < MyEPSILON && std::abs(temp_p_max - 2.0) < MyEPSILON) p_bin = 3;
        else if (std::abs(temp_p_min - 2.0) < MyEPSILON && std::abs(temp_p_max - 2.5) < MyEPSILON) p_bin = 4;
        else if (std::abs(temp_p_min - 2.5) < MyEPSILON && std::abs(temp_p_max - 3.0) < MyEPSILON) p_bin = 5;
        else {
            printf("[ERROR] unknown p bin!\n");
            exit(1);
        }

        if (std::abs(temp_theta_min - 0.22) < MyEPSILON && std::abs(temp_theta_max - 0.56) < MyEPSILON) theta_bin = 0;
        else if (std::abs(temp_theta_min - 0.56) < MyEPSILON && std::abs(temp_theta_max - 1.13) < MyEPSILON) theta_bin = 1;
        else if (std::abs(temp_theta_min - 1.13) < MyEPSILON && std::abs(temp_theta_max - 1.57) < MyEPSILON) theta_bin = 2;
        else if (std::abs(temp_theta_min - 1.57) < MyEPSILON && std::abs(temp_theta_max - 1.88) < MyEPSILON) theta_bin = 3;
        else if (std::abs(temp_theta_min - 1.88) < MyEPSILON && std::abs(temp_theta_max - 2.23) < MyEPSILON) theta_bin = 4;
        else if (std::abs(temp_theta_min - 2.23) < MyEPSILON && std::abs(temp_theta_max - 2.71) < MyEPSILON) theta_bin = 5;
        else {
            printf("[ERROR] unknown theta bin!\n");
            exit(1);
        }

        int bin = theta_bin + 6 * p_bin;

        if (temp_charge == '+') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeE_correction_MC15ri[3][bin] = temp_data_MC_ratio;
            PID_fakeE_uncer_MC15ri[3][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeE_correction_MC15ri[3][bin]) < MyEPSILON) {
                PID_fakeE_correction_MC15ri[3][bin] = 1.0;
                PID_fakeE_uncer_MC15ri[3][bin] = 0.0;
            }

        }
        else if (temp_charge == '-') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeE_correction_MC15ri[2][bin] = temp_data_MC_ratio;
            PID_fakeE_uncer_MC15ri[2][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeE_correction_MC15ri[2][bin]) < MyEPSILON) {
                PID_fakeE_correction_MC15ri[2][bin] = 1.0;
                PID_fakeE_uncer_MC15ri[2][bin] = 0.0;
            }

        }
        else {
            printf("[ERROR] unknown charge!\n");
            exit(1);
        }
    }
    fclose(fp_pi_fromE);

    // Kaon from fake muon
    while (fscanf(fp_K_fromMU, "kaonID,%c,%lf,%lf,%lf,%lf,0.0,1.0,FixedThresh06,0.6,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &temp_charge, &temp_p_min, &temp_p_max, &temp_theta_min, &temp_theta_max, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_data_MC_ratio, &temp_data_MC_uncertainty_stat_up, &temp_data_MC_uncertainty_stat_dn, &temp_data_MC_uncertainty_sys_up, &temp_data_MC_uncertainty_sys_dn) != EOF) {
        int p_bin = -1;
        int theta_bin = -1;

        if (std::abs(temp_p_min - 0.4) < MyEPSILON && std::abs(temp_p_max - 0.5) < MyEPSILON) p_bin = 0;
        else if (std::abs(temp_p_min - 0.5) < MyEPSILON && std::abs(temp_p_max - 0.7) < MyEPSILON) p_bin = 1;
        else if (std::abs(temp_p_min - 0.7) < MyEPSILON && std::abs(temp_p_max - 1.0) < MyEPSILON) p_bin = 2;
        else if (std::abs(temp_p_min - 1.0) < MyEPSILON && std::abs(temp_p_max - 1.5) < MyEPSILON) p_bin = 3;
        else if (std::abs(temp_p_min - 1.5) < MyEPSILON && std::abs(temp_p_max - 2.0) < MyEPSILON) p_bin = 4;
        else if (std::abs(temp_p_min - 2.0) < MyEPSILON && std::abs(temp_p_max - 2.5) < MyEPSILON) p_bin = 5;
        else {
            printf("[ERROR] unknown p bin!\n");
            exit(1);
        }

        if (std::abs(temp_theta_min - 0.4) < MyEPSILON && std::abs(temp_theta_max - 0.64) < MyEPSILON) theta_bin = 0;
        else if (std::abs(temp_theta_min - 0.64) < MyEPSILON && std::abs(temp_theta_max - 0.82) < MyEPSILON) theta_bin = 1;
        else if (std::abs(temp_theta_min - 0.82) < MyEPSILON && std::abs(temp_theta_max - 1.16) < MyEPSILON) theta_bin = 2;
        else if (std::abs(temp_theta_min - 1.16) < MyEPSILON && std::abs(temp_theta_max - 1.46) < MyEPSILON) theta_bin = 3;
        else if (std::abs(temp_theta_min - 1.46) < MyEPSILON && std::abs(temp_theta_max - 1.78) < MyEPSILON) theta_bin = 4;
        else if (std::abs(temp_theta_min - 1.78) < MyEPSILON && std::abs(temp_theta_max - 2.13) < MyEPSILON) theta_bin = 5;
        else if (std::abs(temp_theta_min - 2.13) < MyEPSILON && std::abs(temp_theta_max - 2.22) < MyEPSILON) theta_bin = 6;
        else if (std::abs(temp_theta_min - 2.22) < MyEPSILON && std::abs(temp_theta_max - 2.6) < MyEPSILON) theta_bin = 7;
        else {
            printf("[ERROR] unknown theta bin!\n");
            exit(1);
        }

        int bin = theta_bin + 8 * p_bin;

        if (temp_charge == '+') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeMU_correction_MC15ri[1][bin] = temp_data_MC_ratio;
            PID_fakeMU_uncer_MC15ri[1][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeMU_correction_MC15ri[1][bin]) < MyEPSILON) {
                PID_fakeMU_correction_MC15ri[1][bin] = 1.0;
                PID_fakeMU_uncer_MC15ri[1][bin] = 0.0;
            }

        }
        else if (temp_charge == '-') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeMU_correction_MC15ri[0][bin] = temp_data_MC_ratio;
            PID_fakeMU_uncer_MC15ri[0][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeMU_correction_MC15ri[0][bin]) < MyEPSILON) {
                PID_fakeMU_correction_MC15ri[0][bin] = 1.0;
                PID_fakeMU_uncer_MC15ri[0][bin] = 0.0;
            }

        }
        else {
            printf("[ERROR] unknown charge!\n");
            exit(1);
        }
    }
    fclose(fp_K_fromMU);

    // Pion from fake muon
    while (fscanf(fp_pi_fromMU, "pionID,%c,%lf,%lf,%lf,%lf,0.0,1.0,FixedThresh06,0.6,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &temp_charge, &temp_p_min, &temp_p_max, &temp_theta_min, &temp_theta_max, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_data_MC_ratio, &temp_data_MC_uncertainty_stat_up, &temp_data_MC_uncertainty_stat_dn, &temp_data_MC_uncertainty_sys_up, &temp_data_MC_uncertainty_sys_dn) != EOF) {
        int p_bin = -1;
        int theta_bin = -1;

        if (std::abs(temp_p_min - 0.4) < MyEPSILON && std::abs(temp_p_max - 0.5) < MyEPSILON) p_bin = 0;
        else if (std::abs(temp_p_min - 0.5) < MyEPSILON && std::abs(temp_p_max - 0.7) < MyEPSILON) p_bin = 1;
        else if (std::abs(temp_p_min - 0.7) < MyEPSILON && std::abs(temp_p_max - 1.0) < MyEPSILON) p_bin = 2;
        else if (std::abs(temp_p_min - 1.0) < MyEPSILON && std::abs(temp_p_max - 1.5) < MyEPSILON) p_bin = 3;
        else if (std::abs(temp_p_min - 1.5) < MyEPSILON && std::abs(temp_p_max - 2.0) < MyEPSILON) p_bin = 4;
        else if (std::abs(temp_p_min - 2.0) < MyEPSILON && std::abs(temp_p_max - 2.5) < MyEPSILON) p_bin = 5;
        else {
            printf("[ERROR] unknown p bin!\n");
            exit(1);
        }

        if (std::abs(temp_theta_min - 0.4) < MyEPSILON && std::abs(temp_theta_max - 0.64) < MyEPSILON) theta_bin = 0;
        else if (std::abs(temp_theta_min - 0.64) < MyEPSILON && std::abs(temp_theta_max - 0.82) < MyEPSILON) theta_bin = 1;
        else if (std::abs(temp_theta_min - 0.82) < MyEPSILON && std::abs(temp_theta_max - 1.16) < MyEPSILON) theta_bin = 2;
        else if (std::abs(temp_theta_min - 1.16) < MyEPSILON && std::abs(temp_theta_max - 1.46) < MyEPSILON) theta_bin = 3;
        else if (std::abs(temp_theta_min - 1.46) < MyEPSILON && std::abs(temp_theta_max - 1.78) < MyEPSILON) theta_bin = 4;
        else if (std::abs(temp_theta_min - 1.78) < MyEPSILON && std::abs(temp_theta_max - 2.13) < MyEPSILON) theta_bin = 5;
        else if (std::abs(temp_theta_min - 2.13) < MyEPSILON && std::abs(temp_theta_max - 2.22) < MyEPSILON) theta_bin = 6;
        else if (std::abs(temp_theta_min - 2.22) < MyEPSILON && std::abs(temp_theta_max - 2.6) < MyEPSILON) theta_bin = 7;
        else {
            printf("[ERROR] unknown theta bin!\n");
            exit(1);
        }

        int bin = theta_bin + 8 * p_bin;

        if (temp_charge == '+') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeMU_correction_MC15ri[3][bin] = temp_data_MC_ratio;
            PID_fakeMU_uncer_MC15ri[3][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeMU_correction_MC15ri[3][bin]) < MyEPSILON) {
                PID_fakeMU_correction_MC15ri[3][bin] = 1.0;
                PID_fakeMU_uncer_MC15ri[3][bin] = 0.0;
            }

        }
        else if (temp_charge == '-') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeMU_correction_MC15ri[2][bin] = temp_data_MC_ratio;
            PID_fakeMU_uncer_MC15ri[2][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeMU_correction_MC15ri[2][bin]) < MyEPSILON) {
                PID_fakeMU_correction_MC15ri[2][bin] = 1.0;
                PID_fakeMU_uncer_MC15ri[2][bin] = 0.0;
            }

        }
        else {
            printf("[ERROR] unknown charge!\n");
            exit(1);
        }
    }
    fclose(fp_pi_fromMU);
}

double Corrector_FakePID::GetCorrectionFactorfakeE(int PID_type, int bin_PID, std::string type) {
    if (type == "MC15ri") return PID_fakeE_correction_MC15ri[PID_type][bin_PID];
    else {
        printf("[Corrector_FakePID] Invalid type!\n");
        exit(1);
    }
}

double Corrector_FakePID::GetCorrectionFactorfakeMU(int PID_type, int bin_PID, std::string type) {
    if (type == "MC15ri") return PID_fakeMU_correction_MC15ri[PID_type][bin_PID];
    else {
        printf("[Corrector_FakePID] Invalid type!\n");
        exit(1);
    }
}

double Corrector_FakePID::GetUncertaintyfakeE(int PID_type, int bin_PID, std::string type) {
    if (type == "MC15ri") return PID_fakeE_uncer_MC15ri[PID_type][bin_PID];
    else {
        printf("[Corrector_FakePID] Invalid type!\n");
        exit(1);
    }
}

double Corrector_FakePID::GetUncertaintyfakeMU(int PID_type, int bin_PID, std::string type) {
    if (type == "MC15ri") return PID_fakeMU_uncer_MC15ri[PID_type][bin_PID];
    else {
        printf("[Corrector_FakePID] Invalid type!\n");
        exit(1);
    }
}

# define N_pi0_syst 8
class Corrector_pi0 {
private:
    static constexpr double pi0_correction_MC15ri[N_pi0_syst] = { 0.917, 0.965, 0.988, 1.013, 1.042, 1.044, 1.011, 1.0 };
    static constexpr double pi0_stat_uncer_MC15ri[N_pi0_syst] = { 0.004, 0.004, 0.004, 0.005, 0.004, 0.005, 0.005, 0.0 };
    static constexpr double pi0_sys_uncer1_MC15ri[N_pi0_syst] = { 0.049, 0.036, 0.079, 0.058, 0.045, 0.041, 0.040, 0.0 };
    static constexpr double pi0_sys_uncer2_MC15ri[N_pi0_syst] = { 0.0, 0.0, 0.0, 0.0, 0.039, 0.051, 0.030, 0.0 };
public:
    Corrector_pi0();
    double GetCorrectionFactor(int bin_pi0, std::string type);
    double GetStatUncertainty(int bin_pi0, std::string type);
    double GetSystUncertainty1(int bin_pi0, std::string type);
    double GetSystUncertainty2(int bin_pi0, std::string type);
};

Corrector_pi0 corrector_pi0;

Corrector_pi0::Corrector_pi0() {}

double Corrector_pi0::GetCorrectionFactor(int bin_pi0, std::string type) {
    if (type == "MC15ri") return pi0_correction_MC15ri[bin_pi0];
    else {
        printf("[Corrector_pi0] Invalid type!\n");
        exit(1);
    }
}

double Corrector_pi0::GetStatUncertainty(int bin_pi0, std::string type) {
    if (type == "MC15ri") return pi0_stat_uncer_MC15ri[bin_pi0];
    else {
        printf("[Corrector_pi0] Invalid type!\n");
        exit(1);
    }
}

double Corrector_pi0::GetSystUncertainty1(int bin_pi0, std::string type) {
    if (type == "MC15ri") return pi0_sys_uncer1_MC15ri[bin_pi0];
    else {
        printf("[Corrector_pi0] Invalid type!\n");
        exit(1);
    }
}

double Corrector_pi0::GetSystUncertainty2(int bin_pi0, std::string type) {
    if (type == "MC15ri") return pi0_sys_uncer2_MC15ri[bin_pi0];
    else {
        printf("[Corrector_pi0] Invalid type!\n");
        exit(1);
    }
}

# define N_Knn_type 4 //  B2Knn B2Kstarnn B02K0nn B02K0starnn
const double B2Knn_up_uncer[N_Knn_type] = { // relative uncertainty
    0.32 / 2.66, 0.28 / 1.24, 0.5 / 5.9, 0.8 / 3.6
};
const double B2Knn_dn_uncer[N_Knn_type] = { // relative uncertainty
    0.32 / 2.66, 0.25 / 1.24, 0.5 / 5.9, 0.7 / 3.6
};

int N_correction_multiplicity_weightfile = 0; // maximum value of Ngamma in weightfile
int N_correction_multiplicity_MCsample = -1; // maximum value of Ngamma in MC samples
TH1D* weights_Ngamma_uncer; // multiplication factor. If this number is multiplied to the number of event, result is the number of event with 1 sigma (nominal + 1 sigma)

class Corrector_Fragmentation {
private:
    static constexpr int N_Category_gamma = 9;
    static constexpr int N_Bin_gamma = 5;
    static constexpr double Bins_gamma[N_Bin_gamma - 1] = { 1.15, 1.5, 2.0, 2.4 };
    static constexpr double Fragmentation_Xsgamma[N_Bin_gamma][N_Category_gamma] = { // unit: percentage with respect to "un-missing" mode
        {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00},
        {10.6, 5.88, 23.2, 44.5, 0.46, 9.99, 0.52, 4.78, 0.00},
        {3.12, 1.13, 15.7, 20.6, 9.48, 26.9, 5.29, 15.5, 2.20},
        {1.65, 0.51, 9.02, 7.98, 16.8, 26.7, 10.3, 21.8, 5.38},
        {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}
    };
    static constexpr double Fragmentation_Uncertainty_Xsgamma[N_Bin_gamma][N_Category_gamma] = { // relative uncertainty
        {0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000},
        {0.0651, 0.0990, 0.0554, 0.0400, 0.6304, 0.1175, 0.7115, 0.2395, 0.0000},
        {0.1053, 0.2727, 0.0686, 0.0741, 0.0942, 0.0733, 0.2777, 0.1715, 0.1168},
        {0.4250, 1.4490, 0.2646, 0.4003, 0.2509, 0.2671, 0.8409, 0.6682, 0.2898},
        {0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000}
    };

    static constexpr double Nevt_Nominal_before_Xsgamma_MC15[N_Bin_gamma][N_Category_gamma] = {
        {228.010067, 113.889034, 97.106455, 99.597443, 2.350444, 17.694690, 0.265738, 25.726503, 0.000112}, // [-inf, 1.15]
        {1217.857056, 608.996474, 683.434991, 821.802920, 54.319919, 320.575746, 29.362053, 194.574254, 0.829998}, // [1.15, 1.5]
        {745.156145, 372.620424, 554.609584, 717.744130, 132.646240, 507.837684, 190.020473, 312.058874, 120.779580}, // [1.5, 2.0]
        {125.738141, 62.841816, 119.699133, 149.610098, 58.094970, 157.286832, 123.180903, 139.932048, 49.973636}, // [2.0, 2.4]
        {22.820347, 11.451874, 24.859499, 30.438579, 15.176837, 37.050655, 37.906534, 38.685844, 11.884932} // [2.4, inf]
    };
    static constexpr double Nevt_Hmb_before_Xsgamma_MC15[N_Bin_gamma][N_Category_gamma] = {
        {234.717370, 116.931851, 100.045961, 102.811438, 2.443726, 18.213362, 0.283348, 26.352459, 0.000000}, // [-inf, 1.15]
        {1237.919103, 618.966253, 694.727432, 833.619261, 54.873905, 324.784821, 29.598600, 196.698641, 0.862629}, // [1.15, 1.5]
        {739.765455, 368.897914, 550.204743, 712.129243, 131.851781, 504.040211, 187.911583, 309.873035, 119.610335}, // [1.5, 2.0]
        {122.975491, 61.761023, 116.640122, 146.325654, 56.851972, 153.571830, 120.538083, 137.029154, 49.297575}, // [2.0, 2.4]
        {22.077625, 11.336836, 24.298764, 29.398639, 14.646905, 35.886317, 36.466259, 37.697931, 11.619447} // [2.4, inf]
    };
    static constexpr double Nevt_Hpf_before_Xsgamma_MC15[N_Bin_gamma][N_Category_gamma] = {
        {211.414635, 105.941646, 89.984692, 92.845350, 2.132717, 16.210070, 0.250627, 23.805551, 0.000569}, // [-inf, 1.15]
        {1145.424694, 573.572296, 643.905532, 775.309971, 51.161825, 303.573493, 27.884826, 183.335197, 0.805736}, // [1.15, 1.5]
        {736.485305, 369.574161, 550.673698, 711.747889, 132.674662, 505.825855, 190.843751, 312.993403, 121.016612}, // [1.5, 2.0]
        {137.460490, 68.603354, 131.880622, 163.688524, 63.757702, 172.478521, 135.887581, 153.778110, 54.689402}, // [2.0, 2.4]
        {29.109281, 14.773702, 32.144158, 39.024643, 19.707985, 47.791272, 49.197184, 50.159756, 15.393167} // [2.4, inf]
    };
    static constexpr double Nevt_Htransition_before_Xsgamma_MC15[N_Bin_gamma][N_Category_gamma] = {
        {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, // [-inf, 1.15]
        {1126.011312, 562.651975, 658.733369, 809.726662, 55.294719, 337.131227, 31.789385, 187.359070, 0.920546}, // [1.15, 1.5]
        {838.951829, 419.735045, 624.498178, 809.052632, 149.312766, 571.803869, 214.105973, 351.790571, 136.066782}, // [1.5, 2.0]
        {141.771711, 71.022628, 135.144794, 168.801789, 65.504295, 176.864279, 139.131693, 157.662450, 56.155419}, // [2.0, 2.4]
        {25.757474, 12.884694, 27.851586, 34.061168, 16.946366, 41.803329, 42.925015, 43.435334, 13.297957} // [2.4, inf]
    };
    static constexpr double Nevt_Lmb_before_Xsgamma_MC15[N_Bin_gamma][N_Category_gamma] = {
        {219.785163, 110.184748, 93.899905, 96.479643, 2.274163, 17.081852, 0.246506, 24.692561, 0.000000}, // [-inf, 1.15]
        {1195.125885, 599.131943, 672.175942, 810.139571, 53.494062, 316.373461, 29.121914, 191.567680, 0.804038}, // [1.15, 1.5]
        {750.803012, 374.158691, 558.976621, 723.091791, 134.075592, 512.438695, 192.030343, 314.970543, 122.127927}, // [1.5, 2.0]
        {128.850611, 64.181819, 121.971830, 153.567926, 59.265244, 160.822175, 126.599064, 143.704766, 51.333129}, // [2.0, 2.4]
        {23.591953, 11.740537, 25.576218, 31.342805, 15.603211, 38.154400, 39.437907, 39.746047, 12.199234} // [2.4, inf]
    };
    static constexpr double Nevt_Lpf_before_Xsgamma_MC15[N_Bin_gamma][N_Category_gamma] = {
        {250.065120, 125.068681, 106.236894, 109.408332, 2.586983, 19.616770, 0.297757, 28.249527, 0.000000}, // [-inf, 1.15]
        {1307.816445, 653.786860, 732.097197, 879.672426, 58.091454, 342.326888, 31.063100, 208.314601, 0.883302}, // [1.15, 1.5]
        {747.709095, 373.918779, 554.600603, 718.457434, 131.129261, 504.948261, 186.036870, 308.211200, 119.650476}, // [1.5, 2.0]
        {109.711774, 55.152876, 104.872451, 130.564821, 50.850965, 137.527349, 107.684824, 121.713044, 43.870882}, // [2.0, 2.4]
        {16.301694, 8.096028, 17.764354, 21.644725, 10.719722, 26.442162, 26.669031, 27.403075, 8.479496} // [2.4, inf]
    };
    static constexpr double Nevt_Ltransition_before_Xsgamma_MC15[N_Bin_gamma][N_Category_gamma] = {
        {687.014043, 343.100287, 227.441995, 234.500737, 5.644878, 35.769345, 0.466491, 61.879282, 0.000000}, // [-inf, 1.15]
        {1095.635246, 549.051012, 615.383772, 741.741024, 48.963482, 288.728398, 26.536985, 175.580393, 0.736949}, // [1.15, 1.5]
        {671.592637, 336.051687, 499.387200, 646.831974, 119.008094, 456.176166, 171.088734, 281.417112, 108.985261}, // [1.5, 2.0]
        {113.378578, 56.429701, 107.697531, 134.644214, 52.336623, 141.679228, 110.495689, 126.585786, 45.203748}, // [2.0, 2.4]
        {20.750444, 10.246957, 22.414418, 27.261420, 13.458927, 33.490579, 34.139560, 34.988253, 10.795822} // [2.4, inf]
    };

    static constexpr double Nevt_Nominal_missing_before_Xsgamma_MC15[N_Bin_gamma] = { 48.015884, 360.363091, 676.491505, 518.892840, 292.939844 }; // number of missing mode (MC)
    static constexpr double Nevt_Hmb_missing_before_Xsgamma_MC15[N_Bin_gamma] = { 49.450546, 365.252264, 669.873825, 507.353413, 284.711591 }; // number of missing mode (MC)
    static constexpr double Nevt_Hpf_missing_before_Xsgamma_MC15[N_Bin_gamma] = { 44.686813, 340.083522, 676.515003, 574.044827, 388.939028 }; // number of missing mode (MC)
    static constexpr double Nevt_Htransition_missing_before_Xsgamma_MC15[N_Bin_gamma] = { 0.000000, 348.935210, 762.920563, 585.338309, 330.039376 }; // number of missing mode (MC)
    static constexpr double Nevt_Lmb_missing_before_Xsgamma_MC15[N_Bin_gamma] = { 46.502106, 354.423602, 683.769325, 533.204382, 302.346798 }; // number of missing mode (MC)
    static constexpr double Nevt_Lpf_missing_before_Xsgamma_MC15[N_Bin_gamma] = { 52.743632, 386.345000, 669.316781, 448.947976, 200.139315 }; // number of missing mode (MC)
    static constexpr double Nevt_Ltransition_missing_before_Xsgamma_MC15[N_Bin_gamma] = { 93.606503, 323.683107, 609.133190, 468.259061, 263.808865 }; // number of missing mode (MC)

    double Total_Nevt_Nominal_before_Xsgamma_MC15[N_Bin_gamma]; // without missing mode
    double Total_Nevt_Hmb_before_Xsgamma_MC15[N_Bin_gamma];
    double Total_Nevt_Hpf_before_Xsgamma_MC15[N_Bin_gamma];
    double Total_Nevt_Htransition_before_Xsgamma_MC15[N_Bin_gamma];
    double Total_Nevt_Lmb_before_Xsgamma_MC15[N_Bin_gamma];
    double Total_Nevt_Lpf_before_Xsgamma_MC15[N_Bin_gamma];
    double Total_Nevt_Ltransition_before_Xsgamma_MC15[N_Bin_gamma];
public:
    enum class SystType { // what kind of Xs sample is corrected?
        Nominal,
        Hmb,
        Hpf,
        Htransition,
        Lmb,
        Lpf,
        Ltransition
    };
    enum class Sample {
        gamma // B -> Xs gamma is used
    };

    Corrector_Fragmentation();
    int Classify(int Decay[N_decay], Sample sample);
    int GetMxBin(double MXs, Sample sample);
    double GetCorrectionFactor(int Decay[N_decay], double MXs, SystType systtype, Sample sample, std::string type);
    double FluctuateCorrection(int Decay[N_decay], double MXs, SystType systtype, int TargetMxsBin, int TargetCategory, bool IsTargetCategoryUp, Sample sample, std::string type);
    int GetNMxsBin(Sample sample);
    int GetNCategory(Sample sample);
};

Corrector_Fragmentation corrector_Fragmentation;

Corrector_Fragmentation::Corrector_Fragmentation() {
    for (int i = 0; i < N_Bin_gamma; i++) {
        Total_Nevt_Nominal_before_Xsgamma_MC15[i] = 0.0;
        Total_Nevt_Hmb_before_Xsgamma_MC15[i] = 0.0;
        Total_Nevt_Hpf_before_Xsgamma_MC15[i] = 0.0;
        Total_Nevt_Htransition_before_Xsgamma_MC15[i] = 0.0;
        Total_Nevt_Lmb_before_Xsgamma_MC15[i] = 0.0;
        Total_Nevt_Lpf_before_Xsgamma_MC15[i] = 0.0;
        Total_Nevt_Ltransition_before_Xsgamma_MC15[i] = 0.0;
    }

    for (int i = 0; i < N_Bin_gamma; i++) {
        for (int j = 0; j < N_Category_gamma; j++) {
            Total_Nevt_Nominal_before_Xsgamma_MC15[i] = Total_Nevt_Nominal_before_Xsgamma_MC15[i] + Nevt_Nominal_before_Xsgamma_MC15[i][j];
            Total_Nevt_Hmb_before_Xsgamma_MC15[i] = Total_Nevt_Hmb_before_Xsgamma_MC15[i] + Nevt_Hmb_before_Xsgamma_MC15[i][j];
            Total_Nevt_Hpf_before_Xsgamma_MC15[i] = Total_Nevt_Hpf_before_Xsgamma_MC15[i] + Nevt_Hpf_before_Xsgamma_MC15[i][j];
            Total_Nevt_Htransition_before_Xsgamma_MC15[i] = Total_Nevt_Htransition_before_Xsgamma_MC15[i] + Nevt_Htransition_before_Xsgamma_MC15[i][j];
            Total_Nevt_Lmb_before_Xsgamma_MC15[i] = Total_Nevt_Lmb_before_Xsgamma_MC15[i] + Nevt_Lmb_before_Xsgamma_MC15[i][j];
            Total_Nevt_Lpf_before_Xsgamma_MC15[i] = Total_Nevt_Lpf_before_Xsgamma_MC15[i] + Nevt_Lpf_before_Xsgamma_MC15[i][j];
            Total_Nevt_Ltransition_before_Xsgamma_MC15[i] = Total_Nevt_Ltransition_before_Xsgamma_MC15[i] + Nevt_Ltransition_before_Xsgamma_MC15[i][j];
        }
    }
}

int Corrector_Fragmentation::Classify(int Decay[N_decay], Sample sample) {
    if (sample == Corrector_Fragmentation::Sample::gamma) {
        /*
            0: Xs2Kpi_wopi0_MC // Xs -> K pi w/o pi0
            1: Xs2Kpi_wpi0_MC, // Xs -> K pi w/ 1pi0
            2: Xs2Kpipi_wopi0_MC, // Xs -> K pi pi w/o pi0
            3: Xs2Kpipi_wpi0_MC, // Xs -> K pi pi w/ 1pi0
            4: Xs2Kpipipi_wopi0_MC, // Xs -> K pi pi pi w/o pi0
            5: Xs2Kpipipi_wpi0_MC, // Xs -> K pi pi pi w/ 1pi0
            6: Xs2Kpipipipi_atmost1pi0_MC, // Xs -> K pi pi pi pi at most 1pi0
            7: Xs2pi0pi0_MC, // including 2pi0
            8: Xs2KKK_MC, // including 3K
            9: MC_other, // other non-resonant Xs decay
        */
        if (Decay[6] > 0 || Decay[24] > 0) return 0;
        else if (Decay[5] > 0 || Decay[25] > 0) return 1;
        else if (Decay[7] > 0 || Decay[27] > 0) return 2;
        else if (Decay[8] > 0 || Decay[26] > 0) return 3;
        else if (Decay[10] > 0 || Decay[28] > 0) return 4;
        else if (Decay[9] > 0 || Decay[29] > 0) return 5;
        else if (Decay[11] > 0 || Decay[12] > 0 || Decay[30] > 0 || Decay[31] > 0) return 6;
        else if (Decay[13] > 0 || Decay[14] > 0 || Decay[15] > 0 || Decay[32] > 0 || Decay[33] > 0 || Decay[34] > 0) return 7;
        else if (Decay[16] > 0 || Decay[17] > 0 || Decay[18] > 0 || Decay[35] > 0 || Decay[36] > 0 || Decay[37] > 0) return 8;
        else return N_Category_gamma;
    }
    else {
        printf("[Corrector_Fragmentation::Classify] cannot find sample");
        exit(1);
        return 1;
    }
}

int Corrector_Fragmentation::GetMxBin(double MXs, Sample sample) {
    if (sample == Corrector_Fragmentation::Sample::gamma) { // [-inf, 1.15, 1.5, 2.0, 2.4, inf]
        for (int i = 0; i < N_Bin_gamma - 1; i++) {
            if (MXs < Bins_gamma[i]) return i;
        }
        return N_Bin_gamma - 1;
    }
    else {
        printf("[Corrector_Fragmentation::GetMxBin] cannot find sample");
        exit(1);
        return 1;
    }
}

double Corrector_Fragmentation::GetCorrectionFactor(int Decay[N_decay], double MXs, SystType systtype, Sample sample, std::string type) {
    int Category = Classify(Decay, sample);
    int MxsBin = GetMxBin(MXs, sample);

    if (sample == Corrector_Fragmentation::Sample::gamma) {

        if (Category == N_Category_gamma) return 1.0; // no correction if it is missing mode

        if ((MxsBin == 0) || (MxsBin == N_Bin_gamma - 1)) return 1.0; // no correction if Mxs bin is out of range

        if (type == "MC15ri") {
            if (systtype == Corrector_Fragmentation::SystType::Nominal) {
                if (Nevt_Nominal_before_Xsgamma_MC15[MxsBin][Category] < MyEPSILON) return 1.0; // no correction if Nevt is 0 at MC
                return (0.01 * Fragmentation_Xsgamma[MxsBin][Category]) / (Nevt_Nominal_before_Xsgamma_MC15[MxsBin][Category] / Total_Nevt_Nominal_before_Xsgamma_MC15[MxsBin]);
            }
            else if (systtype == Corrector_Fragmentation::SystType::Hmb) {
                if (Nevt_Hmb_before_Xsgamma_MC15[MxsBin][Category] < MyEPSILON) return 1.0; // no correction if Nevt is 0 at MC
                return (0.01 * Fragmentation_Xsgamma[MxsBin][Category]) / (Nevt_Hmb_before_Xsgamma_MC15[MxsBin][Category] / Total_Nevt_Hmb_before_Xsgamma_MC15[MxsBin]);
            }
            else if (systtype == Corrector_Fragmentation::SystType::Hpf) {
                if (Nevt_Hpf_before_Xsgamma_MC15[MxsBin][Category] < MyEPSILON) return 1.0; // no correction if Nevt is 0 at MC
                return (0.01 * Fragmentation_Xsgamma[MxsBin][Category]) / (Nevt_Hpf_before_Xsgamma_MC15[MxsBin][Category] / Total_Nevt_Hpf_before_Xsgamma_MC15[MxsBin]);
            }
            else if (systtype == Corrector_Fragmentation::SystType::Htransition) {
                if (Nevt_Htransition_before_Xsgamma_MC15[MxsBin][Category] < MyEPSILON) return 1.0; // no correction if Nevt is 0 at MC
                return (0.01 * Fragmentation_Xsgamma[MxsBin][Category]) / (Nevt_Htransition_before_Xsgamma_MC15[MxsBin][Category] / Total_Nevt_Htransition_before_Xsgamma_MC15[MxsBin]);
            }
            else if (systtype == Corrector_Fragmentation::SystType::Lmb) {
                if (Nevt_Lmb_before_Xsgamma_MC15[MxsBin][Category] < MyEPSILON) return 1.0; // no correction if Nevt is 0 at MC
                return (0.01 * Fragmentation_Xsgamma[MxsBin][Category]) / (Nevt_Lmb_before_Xsgamma_MC15[MxsBin][Category] / Total_Nevt_Lmb_before_Xsgamma_MC15[MxsBin]);
            }
            else if (systtype == Corrector_Fragmentation::SystType::Lpf) {
                if (Nevt_Lpf_before_Xsgamma_MC15[MxsBin][Category] < MyEPSILON) return 1.0; // no correction if Nevt is 0 at MC
                return (0.01 * Fragmentation_Xsgamma[MxsBin][Category]) / (Nevt_Lpf_before_Xsgamma_MC15[MxsBin][Category] / Total_Nevt_Lpf_before_Xsgamma_MC15[MxsBin]);
            }
            else if (systtype == Corrector_Fragmentation::SystType::Ltransition) {
                if (Nevt_Ltransition_before_Xsgamma_MC15[MxsBin][Category] < MyEPSILON) return 1.0; // no correction if Nevt is 0 at MC
                return (0.01 * Fragmentation_Xsgamma[MxsBin][Category]) / (Nevt_Ltransition_before_Xsgamma_MC15[MxsBin][Category] / Total_Nevt_Ltransition_before_Xsgamma_MC15[MxsBin]);
            }
        }
        else {
            printf("[Corrector_Fragmentation::GetCorrectionFactor] Invalid type!\n");
            exit(1);
            return 0;
        }
    }
    else {
        printf("[Corrector_Fragmentation::GetCorrectionFactor] cannot find sample");
        exit(1);
        return 1;
    }
}

double Corrector_Fragmentation::FluctuateCorrection(int Decay[N_decay], double MXs, SystType systtype, int TargetMxsBin, int TargetCategory, bool IsTargetCategoryUp, Sample sample, std::string type) {
    /*
    Nevt of TargetCategory in TargetMxsBin is fluctuated.
    we calculate the change of Decay[N_decay] (Category)
    */
    int Category = Classify(Decay, sample);
    int MxsBin = GetMxBin(MXs, sample);

    if (TargetMxsBin != MxsBin) return 1.0; // no fluctuation if it is not in target MXs region

    const double RelativeUncertainty = 1.0;

    if (sample == Corrector_Fragmentation::Sample::gamma) {

        double TotalNevtAtMxsBinWithMissing = Total_Nevt_Nominal_before_Xsgamma_MC15[TargetMxsBin] + Nevt_Nominal_missing_before_Xsgamma_MC15[TargetMxsBin];
        double TargetNevtAtMxsBin;
        if (TargetCategory == N_Category_gamma) TargetNevtAtMxsBin = Nevt_Nominal_missing_before_Xsgamma_MC15[TargetMxsBin]; // if it is missing mode
        else TargetNevtAtMxsBin = Total_Nevt_Nominal_before_Xsgamma_MC15[TargetMxsBin] * (Fragmentation_Xsgamma[TargetMxsBin][TargetCategory] * 0.01);
        double TotalNevtAtMxsBinWithMissingWithoutTargetCategory = TotalNevtAtMxsBinWithMissing - TargetNevtAtMxsBin;

        if (Category == TargetCategory) { // it is target category
            if (IsTargetCategoryUp) {
                if ((TotalNevtAtMxsBinWithMissingWithoutTargetCategory - TargetNevtAtMxsBin * RelativeUncertainty) > 0) return (1.0 + RelativeUncertainty); // the number of event in target category is not large
                else return (1.0 + (TotalNevtAtMxsBinWithMissingWithoutTargetCategory / TargetNevtAtMxsBin)); // the number of event in target category is large
            }
            else return (1.0 - RelativeUncertainty);
        }
        else {
            if (IsTargetCategoryUp) {
                double output = (TotalNevtAtMxsBinWithMissingWithoutTargetCategory - TargetNevtAtMxsBin * RelativeUncertainty) / TotalNevtAtMxsBinWithMissingWithoutTargetCategory;
                if (output < 0) return 0.0; // the number of event in target category is large
                else return output; // the number of event in target category is not large
            }
            else {
                double output = (TotalNevtAtMxsBinWithMissingWithoutTargetCategory + TargetNevtAtMxsBin * RelativeUncertainty) / TotalNevtAtMxsBinWithMissingWithoutTargetCategory;
                if (output < 0) {
                    printf("[Corrector_Fragmentation::FluctuateCorrection] minus Nevt!\n");
                    exit(1);
                }
                else return output;
            }
        }
    }
    else {
        printf("[Corrector_Fragmentation::GetCorrectionFactor] cannot find sample");
        exit(1);
        return 1;
    }
}

int Corrector_Fragmentation::GetNMxsBin(Sample sample) {
    if (sample == Corrector_Fragmentation::Sample::gamma) return N_Bin_gamma;
    else {
        printf("[Corrector_Fragmentation::GetMxBin] cannot find sample");
        exit(1);
        return 1;
    }
}

int Corrector_Fragmentation::GetNCategory(Sample sample) {
    if (sample == Corrector_Fragmentation::Sample::gamma) return (N_Category_gamma + 1); // visible + missing
    else {
        printf("[Corrector_Fragmentation::GetMxBin] cannot find sample");
        exit(1);
        return 1;
    }
}

TH1D* Xsu_Hmb_weight;
TH1D* Xsu_Lmb_weight;
TH1D* Xsu_Hpf_weight;
TH1D* Xsu_Lpf_weight;
TH1D* Xsu_Htransition_weight;
TH1D* Xsu_Ltransition_weight;
TH1D* Xsd_Hmb_weight;
TH1D* Xsd_Lmb_weight;
TH1D* Xsd_Hpf_weight;
TH1D* Xsd_Lpf_weight;
TH1D* Xsd_Htransition_weight;
TH1D* Xsd_Ltransition_weight;
TH1D* Kstar_delta_weight;
TH1D* K0star_delta_weight;
double ReadWeightHist(TH1D* hist, double value);

enum DecayMode { // reco level
    B2Kc = 0,
    B2KcPi0,
    B2Ks0Pic,
    B2KcPicPic,
    B2Ks0PicPi0,
    B2KcPicPicPi0,
    B2Ks0PicPicPic,
    B2KcPicPicPicPic,
    B2Ks0PicPicPicPi0,
    B2KcPi0Pi0,
    B2Ks0PicPi0Pi0,
    B2KcPicPicPi0Pi0,
    B2KcKcKc,
    B2KcKcKs0Pic,
    B2KcKcKcPi0,
    B02Ks0,
    B02KcPic,
    B02Ks0Pi0,
    B02KcPicPi0,
    B02Ks0PicPic,
    B02KcPicPicPic,
    B02Ks0PicPicPi0,
    B02KcPicPicPicPi0,
    B02Ks0PicPicPicPic,
    B02Ks0Pi0Pi0,
    B02KcPicPi0Pi0,
    B02Ks0PicPicPi0Pi0,
    B02KcKcKs0,
    B02KcKcKcPic,
    B02KcKcKs0Pi0,
    MAX_NUM_DECAYMODE
};
enum DecayModeMC { // MC level
    Xsu2Kc_MC = 0,
    Xsu2Kcstar2KcPi0_MC,
    Xsu2Kcstar2K0Pic_MC,
    Xsu2KcPi0_MC,
    Xsu2K0Pic_MC,
    Xsu2KcPicPic_MC,
    Xsu2K0PicPi0_MC,
    Xsu2KcPicPicPi0_MC,
    Xsu2K0PicPicPic_MC,
    Xsu2KcPicPicPicPic_MC,
    Xsu2K0PicPicPicPi0_MC,
    Xsu2KcPi0Pi0_MC,
    Xsu2K0PicPi0Pi0_MC,
    Xsu2KcPicPicPi0Pi0_MC,
    Xsu2KcKcKc_MC,
    Xsu2KcKcK0Pic_MC,
    Xsu2KcKcKcPi0_MC,
    Xsd2K0_MC,
    Xsd2K0star2KcPic_MC,
    Xsd2K0star2K0Pi0_MC,
    Xsd2KcPic_MC,
    Xsd2K0Pi0_MC,
    Xsd2KcPicPi0_MC,
    Xsd2K0PicPic_MC,
    Xsd2KcPicPicPic_MC,
    Xsd2K0PicPicPi0_MC,
    Xsd2KcPicPicPicPi0_MC,
    Xsd2K0PicPicPicPic_MC,
    Xsd2K0Pi0Pi0_MC,
    Xsd2KcPicPi0Pi0_MC,
    Xsd2K0PicPicPi0Pi0_MC,
    Xsd2KcKcK0_MC,
    Xsd2KcKcKcPic_MC,
    Xsd2KcKcK0Pi0_MC,
    other,
    MAX_NUM_DECAYMODE_MC
};

void DecayArrayToXsOutputDecay(int Decay[N_decay], double output_Decay[MAX_NUM_DECAYMODE_MC]) {
    output_Decay[Xsu2Kc_MC] = 0;
    output_Decay[Xsu2Kcstar2KcPi0_MC] = 0;
    output_Decay[Xsu2Kcstar2K0Pic_MC] = 0;
    output_Decay[Xsu2KcPi0_MC] = Decay[5];
    output_Decay[Xsu2K0Pic_MC] = Decay[6];
    output_Decay[Xsu2KcPicPic_MC] = Decay[7];
    output_Decay[Xsu2K0PicPi0_MC] = Decay[8];
    output_Decay[Xsu2KcPicPicPi0_MC] = Decay[9];
    output_Decay[Xsu2K0PicPicPic_MC] = Decay[10];
    output_Decay[Xsu2KcPicPicPicPic_MC] = Decay[11];
    output_Decay[Xsu2K0PicPicPicPi0_MC] = Decay[12];
    output_Decay[Xsu2KcPi0Pi0_MC] = Decay[13];
    output_Decay[Xsu2K0PicPi0Pi0_MC] = Decay[14];
    output_Decay[Xsu2KcPicPicPi0Pi0_MC] = Decay[15];
    output_Decay[Xsu2KcKcKc_MC] = Decay[16];
    output_Decay[Xsu2KcKcK0Pic_MC] = Decay[17];
    output_Decay[Xsu2KcKcKcPi0_MC] = Decay[18];
    output_Decay[Xsd2K0_MC] = 0;
    output_Decay[Xsd2K0star2KcPic_MC] = 0;
    output_Decay[Xsd2K0star2K0Pi0_MC] = 0;
    output_Decay[Xsd2KcPic_MC] = Decay[24];
    output_Decay[Xsd2K0Pi0_MC] = Decay[25];
    output_Decay[Xsd2KcPicPi0_MC] = Decay[26];
    output_Decay[Xsd2K0PicPic_MC] = Decay[27];
    output_Decay[Xsd2KcPicPicPic_MC] = Decay[28];
    output_Decay[Xsd2K0PicPicPi0_MC] = Decay[29];
    output_Decay[Xsd2KcPicPicPicPi0_MC] = Decay[30];
    output_Decay[Xsd2K0PicPicPicPic_MC] = Decay[31];
    output_Decay[Xsd2K0Pi0Pi0_MC] = Decay[32];
    output_Decay[Xsd2KcPicPi0Pi0_MC] = Decay[33];
    output_Decay[Xsd2K0PicPicPi0Pi0_MC] = Decay[34];
    output_Decay[Xsd2KcKcK0_MC] = Decay[35];
    output_Decay[Xsd2KcKcKcPic_MC] = Decay[36];
    output_Decay[Xsd2KcKcK0Pi0_MC] = Decay[37];

    output_Decay[other] = 0.0;
    for (int i = 0; i < MAX_NUM_DECAYMODE_MC; i++) {
        if (output_Decay[i] > 0.5) return;
    }

    output_Decay[other] = 1.0;
}

double GetFragmentationWeight(const char* type, DecayModeMC SelectedDecayMode, bool IsItUp) {
    if (strcmp(type, "Bplus") == 0 && ((SelectedDecayMode >= Xsu2KcPi0_MC && SelectedDecayMode <= Xsu2KcKcKcPi0_MC) || (SelectedDecayMode == other))) {
        double Decay_Num[MAX_NUM_DECAYMODE_MC] = { 0.0, 0.0, 0.0, 
            Xsu_frag_decay1, Xsu_frag_decay2, Xsu_frag_decay3, Xsu_frag_decay4, 
            Xsu_frag_decay5, Xsu_frag_decay6, Xsu_frag_decay7, Xsu_frag_decay8, 
            Xsu_frag_decay9, Xsu_frag_decay10, Xsu_frag_decay11, Xsu_frag_decay12, 
            Xsu_frag_decay13, Xsu_frag_decay14, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Xsu_frag_decay30 };
        double Total_Decay_Num = 0;
        for (int i = 0; i < MAX_NUM_DECAYMODE_MC; i++) Total_Decay_Num = Total_Decay_Num + Decay_Num[i];

        double Selected_Decay_Num = 0;
        switch (SelectedDecayMode) {
        case Xsu2KcPi0_MC:
            Selected_Decay_Num = Xsu_frag_decay1;
            break;
        case Xsu2K0Pic_MC:
            Selected_Decay_Num = Xsu_frag_decay2;
            break;
        case Xsu2KcPicPic_MC:
            Selected_Decay_Num = Xsu_frag_decay3;
            break;
        case Xsu2K0PicPi0_MC:
            Selected_Decay_Num = Xsu_frag_decay4;
            break;
        case Xsu2KcPicPicPi0_MC:
            Selected_Decay_Num = Xsu_frag_decay5;
            break;
        case Xsu2K0PicPicPic_MC:
            Selected_Decay_Num = Xsu_frag_decay6;
            break;
        case Xsu2KcPicPicPicPic_MC:
            Selected_Decay_Num = Xsu_frag_decay7;
            break;
        case Xsu2K0PicPicPicPi0_MC:
            Selected_Decay_Num = Xsu_frag_decay8;
            break;
        case Xsu2KcPi0Pi0_MC:
            Selected_Decay_Num = Xsu_frag_decay9;
            break;
        case Xsu2K0PicPi0Pi0_MC:
            Selected_Decay_Num = Xsu_frag_decay10;
            break;
        case Xsu2KcPicPicPi0Pi0_MC:
            Selected_Decay_Num = Xsu_frag_decay11;
            break;
        case Xsu2KcKcKc_MC:
            Selected_Decay_Num = Xsu_frag_decay12;
            break;
        case Xsu2KcKcK0Pic_MC:
            Selected_Decay_Num = Xsu_frag_decay13;
            break;
        case Xsu2KcKcKcPi0_MC:
            Selected_Decay_Num = Xsu_frag_decay14;
            break;
        case other:
            Selected_Decay_Num = Xsu_frag_decay30;
            break;
        default:
            printf("[ERROR] Input value of SelectedDecayMode is not appropriate!\n");
            exit(1);
        }

        if (IsItUp) return ((Total_Decay_Num - 1.5 * Selected_Decay_Num) / (Total_Decay_Num - Selected_Decay_Num));
        else return ((Total_Decay_Num - 0.5 * Selected_Decay_Num) / (Total_Decay_Num - Selected_Decay_Num));
    
    }
    else if (strcmp(type, "Bzero") == 0 && (SelectedDecayMode >= Xsd2KcPic_MC && SelectedDecayMode <= other)) {

        double Decay_Num[MAX_NUM_DECAYMODE_MC] = { 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, Xsd_frag_decay16,
          Xsd_frag_decay17, Xsd_frag_decay18, Xsd_frag_decay19, Xsd_frag_decay20, Xsd_frag_decay21, Xsd_frag_decay22,
          Xsd_frag_decay23, Xsd_frag_decay24, Xsd_frag_decay25, Xsd_frag_decay26, Xsd_frag_decay27,
          Xsd_frag_decay28, Xsd_frag_decay29, Xsd_frag_decay30 };
        double Total_Decay_Num = 0;
        for (int i = 0; i < MAX_NUM_DECAYMODE_MC; i++) Total_Decay_Num = Total_Decay_Num + Decay_Num[i];

        double Selected_Decay_Num = 0;
        switch (SelectedDecayMode) {
        case Xsd2KcPic_MC:
            Selected_Decay_Num = Xsd_frag_decay16;
            break;
        case Xsd2K0Pi0_MC:
            Selected_Decay_Num = Xsd_frag_decay17;
            break;
        case Xsd2KcPicPi0_MC:
            Selected_Decay_Num = Xsd_frag_decay18;
            break;
        case Xsd2K0PicPic_MC:
            Selected_Decay_Num = Xsd_frag_decay19;
            break;
        case Xsd2KcPicPicPic_MC:
            Selected_Decay_Num = Xsd_frag_decay20;
            break;
        case Xsd2K0PicPicPi0_MC:
            Selected_Decay_Num = Xsd_frag_decay21;
            break;
        case Xsd2KcPicPicPicPi0_MC:
            Selected_Decay_Num = Xsd_frag_decay22;
            break;
        case Xsd2K0PicPicPicPic_MC:
            Selected_Decay_Num = Xsd_frag_decay23;
            break;
        case Xsd2K0Pi0Pi0_MC:
            Selected_Decay_Num = Xsd_frag_decay24;
            break;
        case Xsd2KcPicPi0Pi0_MC:
            Selected_Decay_Num = Xsd_frag_decay25;
            break;
        case Xsd2K0PicPicPi0Pi0_MC:
            Selected_Decay_Num = Xsd_frag_decay26;
            break;
        case Xsd2KcKcK0_MC:
            Selected_Decay_Num = Xsd_frag_decay27;
            break;
        case Xsd2KcKcKcPic_MC:
            Selected_Decay_Num = Xsd_frag_decay28;
            break;
        case Xsd2KcKcK0Pi0_MC:
            Selected_Decay_Num = Xsd_frag_decay29;
            break;
        case other:
            Selected_Decay_Num = Xsd_frag_decay30;
            break;
        default:
            printf("[ERROR] Input value of SelectedDecayMode is not appropriate!\n");
            exit(1);
        }

        if (IsItUp) return ((Total_Decay_Num - 1.5 * Selected_Decay_Num) / (Total_Decay_Num - Selected_Decay_Num));
        else return ((Total_Decay_Num - 0.5 * Selected_Decay_Num) / (Total_Decay_Num - Selected_Decay_Num));

    }
    else if (strcmp(type, "Continuum") == 0) {
        printf("[ERROR] unexpected type name\n");
        printf("[ERROR] Continuum type cannot be selected for fragmentation pdf\n");
        exit(1);
    }
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    printf("[ERROR] unexpected error\n");
    exit(1);
    return 0;
}


//# define Htransition_Xsu_change (-0.142)
//# define Htransition_Xsd_change (-0.099)
//# define Ltransition_Xsu_change (-0.002)
//# define Ltransition_Xsd_change (0.067)
//# define Hmb_Xsu_change (-0.067)
//# define Hmb_Xsd_change (-0.058)
//# define Lmb_Xsu_change (-0.062)
//# define Lmb_Xsd_change (-0.046)
//# define Hpf_Xsu_change (-0.120)
//# define Hpf_Xsd_change (-0.073)
//# define Lpf_Xsu_change (-0.014)
//# define Lpf_Xsd_change (0.019)

//# define FEI_cal_Bc 1.0
//# define FEI_cal_Bc_uncertainty 0.0
//# define FEI_cal_B0 1.0
//# define FEI_cal_B0_uncertainty 0.0

//# define Htransition_Xsu_change 0.0
//# define Htransition_Xsd_change 0.0
//# define Ltransition_Xsu_change 0.0
//# define Ltransition_Xsd_change 0.0
//# define Hmb_Xsu_change 0.0
//# define Hmb_Xsd_change 0.0
//# define Lmb_Xsu_change 0.0
//# define Lmb_Xsd_change 0.0
//# define Hpf_Xsu_change 0.0
//# define Hpf_Xsd_change 0.0
//# define Lpf_Xsu_change 0.0
//# define Lpf_Xsd_change 0.0

# define RarityBins 6
# define BinMIN 0.96
# define BinMAX 1.0

/* ====================================== */
// Rarity module
bool sorting(pair<double, double> a, pair<double, double> b) {

    return a.first < b.first;

}

class CDF {
private:
    std::vector<double> m_values;
    std::vector<double> m_weights;
    double m_weights_sum;

    std::vector<double> m_signal_output;
    std::vector<double> m_signal_weight;
    std::vector<double> m_background_output;
    std::vector<double> m_background_weight;

    const double m_max = 1.0;
    const double m_min = 0.92;

    template <typename A, typename B>
    void zip(const std::vector<A>& a, const std::vector<B>& b, std::vector<std::pair<A, B>>& zipped)
    {
        for (size_t i = 0; i < a.size(); ++i)
        {
            zipped.push_back(std::make_pair(a[i], b[i]));
        }
    }

    template <typename A, typename B>
    void unzip(const std::vector<std::pair<A, B>>& zipped, std::vector<A>& a, std::vector<B>& b)
    {
        for (size_t i = 0; i < a.size(); i++)
        {
            a[i] = zipped[i].first;
            b[i] = zipped[i].second;
        }
    }

    void m_load_files(const char* dirname, std::vector<std::string>* names) {
        TSystemDirectory dir(dirname, dirname);
        TList* files = dir.GetListOfFiles();
        if (files) {
            TSystemFile* file;
            TString fname;
            TIter next(files);
            while ((file = (TSystemFile*)next())) {
                fname = file->GetName();
                if (!file->IsDirectory() && fname.EndsWith(".root")) {
                    names->push_back(fname.Data());
                }
            }
        }
    }


public:
    void init(std::vector<double> values, std::vector<double> weights) {
        // Zip the vectors together
        std::vector<std::pair<double, double>> zipped;
        zip(values, weights, zipped);

        sort(zipped.begin(), zipped.end(), sorting);

        // Write the sorted pairs back to the vectors
        unzip(zipped, m_values, m_weights);
        m_weights_sum = std::accumulate(m_weights.begin(), m_weights.end(), 0.0);
    }
    void initbypath(const char* dirname, double weight) {
        std::vector<string> names;

       m_load_files(dirname, &names);

        for (unsigned int i = 0; i < names.size(); i++) {
            float MVA_BB = 0;

            TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");

            TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
            tree_upsilon->SetBranchAddress("MVA_BB", &MVA_BB);

            for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
                tree_upsilon->GetEntry(j);

                m_values.push_back(MVA_BB);
                m_weights.push_back(weight);

            }
            input_file->Close();
        }

        // Zip the vectors together
        std::vector<std::pair<double, double>> zipped;
        zip(m_values, m_weights, zipped);

        sort(zipped.begin(), zipped.end(), sorting);

        // Write the sorted pairs back to the vectors
        unzip(zipped, m_values, m_weights);
        m_weights_sum = std::accumulate(m_weights.begin(), m_weights.end(), 0.0);
    }
    double GetCDFvalue(float value) {
        const double max = 1.0;
        const double min = 0.85;

        auto larger_or_equal = std::lower_bound(m_values.begin(), m_values.end(), value);
        auto larger = std::upper_bound(m_values.begin(), m_values.end(), value);

        int until_index = -1;
        if (larger_or_equal - m_values.begin() == larger - m_values.begin()) { // input value is not same with one of value which is saved
            until_index = larger_or_equal - m_values.begin() - 1;
        }
        else { // input value is same with one of value which is saved
            until_index = larger_or_equal - m_values.begin();
            double weights_partial_sum = std::accumulate(m_weights.begin(), until_index + m_weights.begin(), 0.0);
            return weights_partial_sum / m_weights_sum;
        }

        if (value < m_values.at(0)) {
            int temp_index = std::upper_bound(m_values.begin(), m_values.end(), m_values.at(0)) - m_values.begin() - 1;
            double weights_partial_sum = (value - m_min) * std::accumulate(m_weights.begin(), temp_index + m_weights.begin(), 0.0) / (m_values.at(0) - m_min);
            return weights_partial_sum / m_weights_sum;
        }
        else if (value >= *(m_values.end() - 1)) {
            return 1.0;
        }

        double weights_partial_sum_back = std::accumulate(m_weights.begin(), until_index + m_weights.begin(), 0.0);
        double weights_partial_sum_front = std::accumulate(m_weights.begin(), until_index + m_weights.begin() + 1, 0.0);

        int index_back = until_index;
        int index_front = until_index + 1;

        double weights_partial_sum = weights_partial_sum_back + (value - m_values.at(index_back)) * (weights_partial_sum_front - weights_partial_sum_back) / (m_values.at(index_front) - m_values.at(index_back));

        return weights_partial_sum / m_weights_sum;

    }
    void print() {
        for (unsigned int i = 0; i < m_values.size(); i++) {
            printf("[%lf, %lf] ", m_values.at(i), m_weights.at(i));
        }
        printf("\n");
    }
    void Evaluate(const char* dirname, double weight, bool IsSignal) {
        std::vector<string> names;

        m_load_files(dirname, &names);

        for (unsigned int i = 0; i < names.size(); i++) {
            float MVA_BB = 0;

            TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");

            TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
            tree_upsilon->SetBranchAddress("MVA_BB", &MVA_BB);

            for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
                tree_upsilon->GetEntry(j);

                if (IsSignal) {
                    m_signal_output.push_back(GetCDFvalue(MVA_BB));
                    m_signal_weight.push_back(weight);
                }
                else {
                    m_background_output.push_back(GetCDFvalue(MVA_BB));
                    m_background_weight.push_back(weight);
                }
            }
            input_file->Close();

        }
    }
    void PrintPNG() {
        TH1D* temp_signal_hist = new TH1D("CDF_signal", ";CDF;arbitrary unit", 100, 0.0, 1.0);
        TH1D* temp_background_hist = new TH1D("CDF_background", ";CDF;arbitrary unit", 100, 0.0, 1.0);

        for (unsigned int i = 0; i < m_signal_output.size(); i++) temp_signal_hist->Fill(m_signal_output.at(i), m_signal_weight.at(i));
        for (unsigned int i = 0; i < m_background_output.size(); i++) temp_background_hist->Fill(m_background_output.at(i), m_background_weight.at(i));

        temp_signal_hist->Scale(1.0 / temp_signal_hist->Integral(), "width");
        temp_background_hist->Scale(1.0 / temp_background_hist->Integral(), "width");

        TCanvas* c_temp = new TCanvas("c", "", 800, 800);

        temp_signal_hist->Draw("hist"); c_temp->SaveAs("signal_cdf.png");
        temp_background_hist->Draw("hist"); c_temp->SaveAs("background_cdf.png");

        delete c_temp;
    }
};

CDF cdf;

class Corrector {
private:
    const int STEP;

    const double mKp;
    const double mK0;

    const double mBp;
    const double mB0;

    const double alpha0_old;
    const double alpha1_old;
    const double alpha2_old;
    const double mp_Bp_old;
    const double mp_B0_old;

    const double alpha0_new; // 0.2545 +- 0.0090
    const double alpha1_new; // -0.71 +- 0.14
    const double alpha2_new; // 0.32 +- 0.59
    const double mp_B_new; //  5.4158 +- 0.0015

    double Total_Bp_old;
    double Total_B0_old;
    double Total_Bp_new;
    double Total_B0_new;
public:
    Corrector();
    double GetCorrectionFactor(double q2, std::string type);
    double GetUncertaintyWeight(double q2, std::string type, double alpha0_other, double alpha1_other, double alpha2_other);
    void Drawq2Plot(bool IsItNew, std::string type);
};

Corrector corrector;

Corrector::Corrector() :
    STEP(2000),
    mKp(0.493677),
    mK0(0.497611),
    mBp(5.27934),
    mB0(5.27965),
    alpha0_old(0.432),
    alpha1_old(-0.664),
    alpha2_old(-1.2),
    mp_Bp_old(5.27934 + 0.046),
    mp_B0_old(5.27965 + 0.046),
    alpha0_new(0.2545),
    alpha1_new(-0.71),
    alpha2_new(0.32),
    mp_B_new(5.4158)
{
    Total_Bp_old = 0;
    Total_B0_old = 0;
    Total_Bp_new = 0;
    Total_B0_new = 0;

    // Bp old
    for (int i = 0; i < STEP; i++) {
        double q2 = ((double)i) * (mBp - mKp) * (mBp - mKp) / ((double)STEP);
        double delq2 = (mBp - mKp) * (mBp - mKp) / ((double)STEP);

        double tp = (mBp + mKp) * (mBp + mKp);
        double tm = (mBp - mKp) * (mBp - mKp);
        double t0 = tp * (1 - std::sqrt(1 - tm / tp));
        double z = (std::sqrt(tp - q2) - std::sqrt(tp - t0)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0));
        double fp = (1 / (1 - q2 / (mp_Bp_old * mp_Bp_old))) * (alpha0_old + alpha1_old * z + alpha2_old * z * z + (-alpha1_old + 2 * alpha2_old) * z * z * z / 3);
        double lambda = std::pow(mBp, 4) + std::pow(mKp, 4) + std::pow(q2, 2) - 2 * (std::pow(mBp, 2) * std::pow(mKp, 2) + std::pow(mKp, 2) * q2 + std::pow(mBp, 2) * q2);

        Total_Bp_old = Total_Bp_old + std::pow(lambda, 1.5) * std::pow(fp, 2) * delq2;
    }

    // B0 old
    for (int i = 0; i < STEP; i++) {
        double q2 = ((double)i) * (mB0 - mK0) * (mB0 - mK0) / ((double)STEP);
        double delq2 = (mB0 - mK0) * (mB0 - mK0) / ((double)STEP);

        double tp = (mB0 + mK0) * (mB0 + mK0);
        double tm = (mB0 - mK0) * (mB0 - mK0);
        double t0 = tp * (1 - std::sqrt(1 - tm / tp));
        double z = (std::sqrt(tp - q2) - std::sqrt(tp - t0)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0));
        double fp = (1 / (1 - q2 / (mp_B0_old * mp_B0_old))) * (alpha0_old + alpha1_old * z + alpha2_old * z * z + (-alpha1_old + 2 * alpha2_old) * z * z * z / 3);
        double lambda = std::pow(mB0, 4) + std::pow(mK0, 4) + std::pow(q2, 2) - 2 * (std::pow(mB0, 2) * std::pow(mK0, 2) + std::pow(mK0, 2) * q2 + std::pow(mB0, 2) * q2);

        Total_B0_old = Total_B0_old + std::pow(lambda, 1.5) * std::pow(fp, 2) * delq2;
    }

    // Bp new
    for (int i = 0; i < STEP; i++) {
        double q2 = ((double)i) * (mBp - mKp) * (mBp - mKp) / ((double)STEP);
        double delq2 = (mBp - mKp) * (mBp - mKp) / ((double)STEP);

        double tp = (mBp + mKp) * (mBp + mKp);
        double tm = (mBp - mKp) * (mBp - mKp);
        double t0 = 0.0;
        double z = (std::sqrt(tp - q2) - std::sqrt(tp - t0)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0));
        double fp = (1 / (1 - q2 / (mp_B_new * mp_B_new))) * (alpha0_new + alpha1_new * z + alpha2_new * z * z + (-alpha1_new + 2 * alpha2_new) * z * z * z / 3);
        double lambda = std::pow(mBp, 4) + std::pow(mKp, 4) + std::pow(q2, 2) - 2 * (std::pow(mBp, 2) * std::pow(mKp, 2) + std::pow(mKp, 2) * q2 + std::pow(mBp, 2) * q2);

        Total_Bp_new = Total_Bp_new + std::pow(lambda, 1.5) * std::pow(fp, 2) * delq2;
    }

    // B0 new
    for (int i = 0; i < STEP; i++) {
        double q2 = ((double)i) * (mB0 - mK0) * (mB0 - mK0) / ((double)STEP);
        double delq2 = (mB0 - mK0) * (mB0 - mK0) / ((double)STEP);

        double tp = (mB0 + mK0) * (mB0 + mK0);
        double tm = (mB0 - mK0) * (mB0 - mK0);
        double t0 = 0.0;
        double z = (std::sqrt(tp - q2) - std::sqrt(tp - t0)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0));
        double fp = (1 / (1 - q2 / (mp_B_new * mp_B_new))) * (alpha0_new + alpha1_new * z + alpha2_new * z * z + (-alpha1_new + 2 * alpha2_new) * z * z * z / 3);
        double lambda = std::pow(mB0, 4) + std::pow(mK0, 4) + std::pow(q2, 2) - 2 * (std::pow(mB0, 2) * std::pow(mK0, 2) + std::pow(mK0, 2) * q2 + std::pow(mB0, 2) * q2);

        Total_B0_new = Total_B0_new + std::pow(lambda, 1.5) * std::pow(fp, 2) * delq2;
    }
}

double Corrector::GetCorrectionFactor(double q2, std::string type) {
    if (type == "Bplus") {
        double tp = (mBp + mKp) * (mBp + mKp);
        double tm = (mBp - mKp) * (mBp - mKp);
        double t0_old = tp * (1 - std::sqrt(1 - tm / tp));
        double t0_new = 0.0;
        double z_old = (std::sqrt(tp - q2) - std::sqrt(tp - t0_old)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0_old));
        double z_new = (std::sqrt(tp - q2) - std::sqrt(tp - t0_new)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0_new));

        double fp_old = (1 / (1 - q2 / (mp_Bp_old * mp_Bp_old))) * (alpha0_old + alpha1_old * z_old + alpha2_old * z_old * z_old + (-alpha1_old + 2 * alpha2_old) * z_old * z_old * z_old / 3);
        double fp_new = (1 / (1 - q2 / (mp_B_new * mp_B_new))) * (alpha0_new + alpha1_new * z_new + alpha2_new * z_new * z_new + (-alpha1_new + 2 * alpha2_new) * z_new * z_new * z_new / 3);
        double lambda = std::pow(mBp, 4) + std::pow(mKp, 4) + std::pow(q2, 2) - 2 * (std::pow(mBp, 2) * std::pow(mKp, 2) + std::pow(mKp, 2) * q2 + std::pow(mBp, 2) * q2);

        if ((std::pow(lambda, 1.5) * std::pow(fp_old, 2) / Total_Bp_old) == 0) {
            printf("[ERROR] try to divide by 0!\n");
            exit(1);
        }
        return (std::pow(lambda, 1.5) * std::pow(fp_new, 2) / Total_Bp_new) / (std::pow(lambda, 1.5) * std::pow(fp_old, 2) / Total_Bp_old);
    }
    else if (type == "Bzero") {
        double tp = (mB0 + mK0) * (mB0 + mK0);
        double tm = (mB0 - mK0) * (mB0 - mK0);
        double t0_old = tp * (1 - std::sqrt(1 - tm / tp));
        double t0_new = 0.0;
        double z_old = (std::sqrt(tp - q2) - std::sqrt(tp - t0_old)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0_old));
        double z_new = (std::sqrt(tp - q2) - std::sqrt(tp - t0_new)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0_new));

        double fp_old = (1 / (1 - q2 / (mp_B0_old * mp_B0_old))) * (alpha0_old + alpha1_old * z_old + alpha2_old * z_old * z_old + (-alpha1_old + 2 * alpha2_old) * z_old * z_old * z_old / 3);
        double fp_new = (1 / (1 - q2 / (mp_B_new * mp_B_new))) * (alpha0_new + alpha1_new * z_new + alpha2_new * z_new * z_new + (-alpha1_new + 2 * alpha2_new) * z_new * z_new * z_new / 3);
        double lambda = std::pow(mB0, 4) + std::pow(mK0, 4) + std::pow(q2, 2) - 2 * (std::pow(mB0, 2) * std::pow(mK0, 2) + std::pow(mK0, 2) * q2 + std::pow(mB0, 2) * q2);

        if ((std::pow(lambda, 1.5) * std::pow(fp_old, 2) / Total_B0_old) == 0) {
            printf("[ERROR] try to divide by 0!\n");
            exit(1);
        }
        return (std::pow(lambda, 1.5) * std::pow(fp_new, 2) / Total_B0_new) / (std::pow(lambda, 1.5) * std::pow(fp_old, 2) / Total_B0_old);
    }

    printf("[ERROR] undefined type name");
    exit(1);

    return 0;
}

double Corrector::GetUncertaintyWeight(double q2, std::string type, double alpha0_other, double alpha1_other, double alpha2_other) {
    if (type == "Bplus") {
        double tp = (mBp + mKp) * (mBp + mKp);
        double tm = (mBp - mKp) * (mBp - mKp);
        double t0_new = 0.0;
        double z_new = (std::sqrt(tp - q2) - std::sqrt(tp - t0_new)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0_new));

        double fp_new = (1 / (1 - q2 / (mp_B_new * mp_B_new))) * (alpha0_new + alpha1_new * z_new + alpha2_new * z_new * z_new + (-alpha1_new + 2 * alpha2_new) * z_new * z_new * z_new / 3);
        double fp_other = (1 / (1 - q2 / (mp_B_new * mp_B_new))) * (alpha0_other + alpha1_other * z_new + alpha2_other * z_new * z_new + (-alpha1_other + 2 * alpha2_other) * z_new * z_new * z_new / 3);
        double lambda = std::pow(mBp, 4) + std::pow(mKp, 4) + std::pow(q2, 2) - 2 * (std::pow(mBp, 2) * std::pow(mKp, 2) + std::pow(mKp, 2) * q2 + std::pow(mBp, 2) * q2);

        if (std::pow(lambda, 1.5) * std::pow(fp_new, 2) == 0) {
            printf("[ERROR] try to divide by 0!\n");
            exit(1);
        }
        return (std::pow(lambda, 1.5) * std::pow(fp_other, 2)) / (std::pow(lambda, 1.5) * std::pow(fp_new, 2));
    }
    else if (type == "Bzero") {
        double tp = (mB0 + mK0) * (mB0 + mK0);
        double tm = (mB0 - mK0) * (mB0 - mK0);
        double t0_new = 0.0;
   
        double z_new = (std::sqrt(tp - q2) - std::sqrt(tp - t0_new)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0_new));

        double fp_new = (1 / (1 - q2 / (mp_B_new * mp_B_new))) * (alpha0_new + alpha1_new * z_new + alpha2_new * z_new * z_new + (-alpha1_new + 2 * alpha2_new) * z_new * z_new * z_new / 3);
        double fp_other = (1 / (1 - q2 / (mp_B_new * mp_B_new))) * (alpha0_other + alpha1_other * z_new + alpha2_other * z_new * z_new + (-alpha1_other + 2 * alpha2_other) * z_new * z_new * z_new / 3);
        double lambda = std::pow(mB0, 4) + std::pow(mK0, 4) + std::pow(q2, 2) - 2 * (std::pow(mB0, 2) * std::pow(mK0, 2) + std::pow(mK0, 2) * q2 + std::pow(mB0, 2) * q2);

        if (std::pow(lambda, 1.5) * std::pow(fp_new, 2) == 0) {
            printf("[ERROR] try to divide by 0!\n");
            exit(1);
        }
        return (std::pow(lambda, 1.5) * std::pow(fp_other, 2)) / (std::pow(lambda, 1.5) * std::pow(fp_new, 2));
    }
}

class Corrector_Knn {
private:

    // K+nn
    int STEP_Knn;
    double mininvM_Knn;
    double maxinvM_Knn;
    TH1D* weights_Knn;
    const double DECAY_DEC_BR_Knn;
    const double new_BR_K0pp;
    double new_BR_Knn;
    const double Nraw_initial_Knn;
    double Nscale_initial_Knn;

    // K*+nn
    int STEP_Kstarnn;
    double mininvM_Kstarnn;
    double maxinvM_Kstarnn;
    TH1D* weights_Kstarnn;
    const double DECAY_DEC_BR_Kstarnn;
    const double new_BR_K0starpp;
    double new_BR_Kstarnn;
    const double Nraw_initial_Kstarnn;
    double Nscale_initial_Kstarnn;

    // K0nn
    int STEP_K0nn;
    double mininvM_K0nn;
    double maxinvM_K0nn;
    TH1D* weights_K0nn;
    const double DECAY_DEC_BR_K0nn;
    const double new_BR_Kpp;
    double new_BR_K0nn;
    const double Nraw_initial_K0nn;
    double Nscale_initial_K0nn;

    // K0*nn
    int STEP_K0starnn;
    double mininvM_K0starnn;
    double maxinvM_K0starnn;
    TH1D* weights_K0starnn;
    const double DECAY_DEC_BR_K0starnn;
    const double new_BR_Kstarpp;
    double new_BR_K0starnn;
    const double Nraw_initial_K0starnn;
    double Nscale_initial_K0starnn;

    const double N_EPSILON;
    const double CUTOFF;

    const double tau_Bp;
    const double tau_B0;

public:
    Corrector_Knn();
    double GetCorrectionFactor(double invM_Knn, double invM_Kstarnn, double invM_K0nn, double invM_K0starnn, double N_Knn, double N_Kstarnn, double N_K0nn, double N_K0starnn);
    double GetCorrectionFactorAtGeneric(double invM_Knn, double invM_Kstarnn, double invM_K0nn, double invM_K0starnn, double N_Knn, double N_Kstarnn, double N_K0nn, double N_K0starnn);
};

Corrector_Knn corrector_Knn;

Corrector_Knn::Corrector_Knn() :
    DECAY_DEC_BR_Knn(0.0000057),
    new_BR_K0pp(0.00000266),
    Nraw_initial_Knn(1000000.0),
    DECAY_DEC_BR_Kstarnn(0.0000057),
    new_BR_K0starpp(0.00000124),
    Nraw_initial_Kstarnn(1000000.0),
    DECAY_DEC_BR_K0nn(0.000002),
    new_BR_Kpp(0.0000059),
    Nraw_initial_K0nn(1000000.0),
    DECAY_DEC_BR_K0starnn(0.0000056),
    new_BR_Kstarpp(0.0000036),
    Nraw_initial_K0starnn(1000000.0),
    N_EPSILON(0.01),
    CUTOFF(50.0),
    tau_Bp(1.6384), // ps
    tau_B0(1.5195) // ps
{
    FILE* fp;

    // read Knn weights
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/Knn_weight/Knn_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &STEP_Knn, &mininvM_Knn, &maxinvM_Knn);
    weights_Knn = new TH1D("Knn_weights", ";;", STEP_Knn, mininvM_Knn, maxinvM_Knn);
    for (int i = 0; i < STEP_Knn; i++) {
        double temp;
        fscanf(fp, "%lf\n", &temp);
        if (temp < CUTOFF) weights_Knn->SetBinContent(i + 1, temp);
        else weights_Knn->SetBinContent(i + 1, CUTOFF);
    }
    fclose(fp);

    // read Kstarnn weights
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/Knn_weight/Kstarnn_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &STEP_Kstarnn, &mininvM_Kstarnn, &maxinvM_Kstarnn);
    weights_Kstarnn = new TH1D("Kstarnn_weights", ";;", STEP_Kstarnn, mininvM_Kstarnn, maxinvM_Kstarnn);
    for (int i = 0; i < STEP_Kstarnn; i++) {
        double temp;
        fscanf(fp, "%lf\n", &temp);
        if (temp < CUTOFF) weights_Kstarnn->SetBinContent(i + 1, temp);
        else weights_Kstarnn->SetBinContent(i + 1, CUTOFF);
    }
    fclose(fp);

    // read K0nn weights
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/Knn_weight/K0nn_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &STEP_K0nn, &mininvM_K0nn, &maxinvM_K0nn);
    weights_K0nn = new TH1D("K0nn_weights", ";;", STEP_K0nn, mininvM_K0nn, maxinvM_K0nn);
    for (int i = 0; i < STEP_K0nn; i++) {
        double temp;
        fscanf(fp, "%lf\n", &temp);
        if (temp < CUTOFF) weights_K0nn->SetBinContent(i + 1, temp);
        else weights_K0nn->SetBinContent(i + 1, CUTOFF);
    }
    fclose(fp);

    // read K0starnn weights
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/Knn_weight/K0starnn_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &STEP_K0starnn, &mininvM_K0starnn, &maxinvM_K0starnn);
    weights_K0starnn = new TH1D("K0starnn_weights", ";;", STEP_K0starnn, mininvM_K0starnn, maxinvM_K0starnn);
    for (int i = 0; i < STEP_K0starnn; i++) {
        double temp;
        fscanf(fp, "%lf\n", &temp);
        if (temp < CUTOFF) weights_K0starnn->SetBinContent(i + 1, temp);
        else weights_K0starnn->SetBinContent(i + 1, CUTOFF);
    }
    fclose(fp);

    // calculate the expected number of event
    // use proper isospin
    new_BR_Knn = new_BR_K0pp * (tau_Bp / tau_B0);
    new_BR_Kstarnn = new_BR_K0starpp * (tau_Bp / tau_B0);
    new_BR_K0nn = new_BR_Kpp * (tau_B0 / tau_Bp);
    new_BR_K0starnn = new_BR_Kstarpp * (tau_B0 / tau_Bp);

    Nscale_initial_Knn = (2.0 * N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0)) * new_BR_Knn);
    Nscale_initial_Kstarnn = (2.0 * N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0)) * new_BR_Kstarnn);
    Nscale_initial_K0nn = (2.0 * N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0)) * new_BR_K0nn);
    Nscale_initial_K0starnn = (2.0 * N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0)) * new_BR_K0starnn);
}

double Corrector_Knn::GetCorrectionFactor(double invM_Knn, double invM_Kstarnn, double invM_K0nn, double invM_K0starnn, double N_Knn, double N_Kstarnn, double N_K0nn, double N_K0starnn) {

    double Correction_Knn = 1;
    double Correction_Kstarnn = 1;
    double Correction_K0nn = 1;
    double Correction_K0starnn = 1;

    if (N_Knn < N_EPSILON) Correction_Knn = 1;
    else {
        int Bin = weights_Knn->FindBin(invM_Knn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_Knn) Bin = STEP_Knn;
        Correction_Knn = std::pow((Nscale_initial_Knn / Nraw_initial_Knn) * weights_Knn->GetBinContent(Bin), N_Knn); // BR correction * invM correction
    }

    if (N_Kstarnn < N_EPSILON) Correction_Kstarnn = 1;
    else {
        int Bin = weights_Kstarnn->FindBin(invM_Kstarnn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_Kstarnn) Bin = STEP_Kstarnn;
        Correction_Kstarnn = std::pow((Nscale_initial_Kstarnn / Nraw_initial_Kstarnn) * weights_Kstarnn->GetBinContent(Bin), N_Kstarnn); // BR correction * invM correction
    }

    if (N_K0nn < N_EPSILON) Correction_K0nn = 1;
    else {
        int Bin = weights_K0nn->FindBin(invM_K0nn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_K0nn) Bin = STEP_K0nn;
        Correction_K0nn = std::pow((Nscale_initial_K0nn / Nraw_initial_K0nn) * weights_K0nn->GetBinContent(Bin), N_K0nn); // BR correction * invM correction
    }

    if (N_K0starnn < N_EPSILON) Correction_K0starnn = 1;
    else {
        int Bin = weights_K0starnn->FindBin(invM_K0starnn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_K0starnn) Bin = STEP_K0starnn;
        Correction_K0starnn = std::pow((Nscale_initial_K0starnn / Nraw_initial_K0starnn) * weights_K0starnn->GetBinContent(Bin), N_K0starnn); // BR correction * invM correction
    }

    return Correction_Knn * Correction_Kstarnn * Correction_K0nn * Correction_K0starnn;
}

double Corrector_Knn::GetCorrectionFactorAtGeneric(double invM_Knn, double invM_Kstarnn, double invM_K0nn, double invM_K0starnn, double N_Knn, double N_Kstarnn, double N_K0nn, double N_K0starnn) {

    double Correction_Knn = 1;
    double Correction_Kstarnn = 1;
    double Correction_K0nn = 1;
    double Correction_K0starnn = 1;

    if (N_Knn < N_EPSILON) Correction_Knn = 1;
    else {
        int Bin = weights_Knn->FindBin(invM_Knn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_Knn) Bin = STEP_Knn;
        Correction_Knn = std::pow((new_BR_Knn / DECAY_DEC_BR_Knn) * weights_Knn->GetBinContent(Bin), N_Knn); // BR correction * invM correction
    }

    if (N_Kstarnn < N_EPSILON) Correction_Kstarnn = 1;
    else {
        int Bin = weights_Kstarnn->FindBin(invM_Kstarnn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_Kstarnn) Bin = STEP_Kstarnn;
        Correction_Kstarnn = std::pow((new_BR_Kstarnn / DECAY_DEC_BR_Kstarnn) * weights_Kstarnn->GetBinContent(Bin), N_Kstarnn); // BR correction * invM correction
    }

    if (N_K0nn < N_EPSILON) Correction_K0nn = 1;
    else {
        int Bin = weights_K0nn->FindBin(invM_K0nn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_K0nn) Bin = STEP_K0nn;
        Correction_K0nn = std::pow((new_BR_K0nn / DECAY_DEC_BR_K0nn) * weights_K0nn->GetBinContent(Bin), N_K0nn); // BR correction * invM correction
    }

    if (N_K0starnn < N_EPSILON) Correction_K0starnn = 1;
    else {
        int Bin = weights_K0starnn->FindBin(invM_K0starnn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_K0starnn) Bin = STEP_K0starnn;
        Correction_K0starnn = std::pow((new_BR_K0starnn / DECAY_DEC_BR_K0starnn) * weights_K0starnn->GetBinContent(Bin), N_K0starnn); // BR correction * invM correction
    }

    return Correction_Knn * Correction_Kstarnn * Correction_K0nn * Correction_K0starnn;
}

class Corrector_Multiplicity {
private:

    int NgammaMAX;
    TH1D* weights_Ngamma;
    const double CUTOFF;

public:
    Corrector_Multiplicity();
    double GetCorrectionFactor(double Ngamma);
};

Corrector_Multiplicity corrector_Multiplicity;

Corrector_Multiplicity::Corrector_Multiplicity() :
    CUTOFF(50.0)
{
    FILE* fp;

    // read Knn weights
    fp = fopen("./multiplicity_weight.txt", "r");
    fscanf(fp, "%d\n", &NgammaMAX);
    weights_Ngamma = new TH1D("weights_Ngamma", ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);
    for (int i = 0; i < NgammaMAX + 1; i++) {
        double temp1;
        double temp2;
        double temp3;
        fscanf(fp, "%lf %lf %lf\n", &temp1, &temp2, &temp3);
        if (temp3 < CUTOFF) weights_Ngamma->SetBinContent(i + 1, temp3);
        else weights_Ngamma->SetBinContent(i + 1, CUTOFF);
    }
    fclose(fp);

}

double Corrector_Multiplicity::GetCorrectionFactor(double Ngamma) {
    int Bin = weights_Ngamma->FindBin(Ngamma);
    if (Bin < 1) {
        printf("[ERROR] Ngamma is smaller than 0!\n");
        exit(1);
    }
    else if (Bin > NgammaMAX + 1) return 1.0;

    return weights_Ngamma->GetBinContent(Bin);
}

class Corrector_KpKLKL {
private:

    int s13_NBin;
    double s13_min;
    double s13_max;

    int s23_NBin;
    double s23_min;
    double s23_max;

    TH2D* weights_KpKLKL;

    const double N_EPSILON;

public:
    Corrector_KpKLKL();
    double GetCorrectionFactorAtGeneric(double s13, double s23, double nB2KpKLKL_all, double nB2KpKLKL_NR);
};

Corrector_KpKLKL corrector_KpKLKL;

Corrector_KpKLKL::Corrector_KpKLKL() :
    N_EPSILON(0.01)
{
    FILE* fp;

    // read KpKLKL weights
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/KpKLKL/KpKLKL_weight.txt", "r");
    fscanf(fp, "s13: %d %lf %lf\n", &s13_NBin, &s13_min, &s13_max);
    fscanf(fp, "s23: %d %lf %lf\n", &s23_NBin, &s23_min, &s23_max);
    weights_KpKLKL = new TH2D("KpKLKL_weights", ";;", s13_NBin, s13_min, s13_max, s23_NBin, s23_min, s23_max);
    for (int i = 0; i < s13_NBin; i++) {
        for (int j = 0; j < s23_NBin; j++) {
            double s13;
            double s23;
            double weight;
            fscanf(fp, "%lf %lf %lf\n", &s13, &s23, &weight);
            weights_KpKLKL->Fill(s13, s23, weight);
        }
    }
    fclose(fp);

}

double Corrector_KpKLKL::GetCorrectionFactorAtGeneric(double s13, double s23, double nB2KpKLKL_all, double nB2KpKLKL_NR) {

    if (nB2KpKLKL_all < N_EPSILON) return 1.0; // no correction needed
    if (nB2KpKLKL_all - nB2KpKLKL_NR > N_EPSILON) return 0.0; // remove B+ --> K+ [X --> KL0 KL0]
    if (nB2KpKLKL_all < 0 || nB2KpKLKL_NR < 0) {
        printf("[Corrector_KpKLKL] number of decay is smaller than 0!\n");
        exit(1);
    }

    // check s13 and s23
    double s13_ = std::min(s13, s23);
    double s23_ = std::max(s13, s23);

    double Correction = 1;

    int GLobalBin_weight = weights_KpKLKL->FindBin(s13_, s23_);
    Correction = weights_KpKLKL->GetBinContent(GLobalBin_weight);

    return Correction;
}

class Corrector_KSKLKL {
private:

    int smax_NBin;
    double smax_min;
    double smax_max;

    int smin_NBin;
    double smin_min;
    double smin_max;

    TH2D* weights_KSKLKL;

    const double N_EPSILON;

public:
    Corrector_KSKLKL();
    double GetCorrectionFactorAtGeneric(double smax, double smin, double nB2KSKLKL_all, double nB2KSKLKL_NR);
};

Corrector_KSKLKL corrector_KSKLKL;

Corrector_KSKLKL::Corrector_KSKLKL() :
    N_EPSILON(0.01)
{
    FILE* fp;

    // read KSKLKL weights
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/KSKLKL/KSKLKL_weight.txt", "r");
    fscanf(fp, "smax: %d %lf %lf\n", &smax_NBin, &smax_min, &smax_max);
    fscanf(fp, "smin: %d %lf %lf\n", &smin_NBin, &smin_min, &smin_max);
    weights_KSKLKL = new TH2D("KSKLKL_weights", ";;", smax_NBin, smax_min, smax_max, smin_NBin, smin_min, smin_max);
    for (int i = 0; i < smax_NBin; i++) {
        for (int j = 0; j < smin_NBin; j++) {
            double smax;
            double smin;
            double weight;
            fscanf(fp, "%lf %lf %lf\n", &smax, &smin, &weight);
            weights_KSKLKL->Fill(smax, smin, weight);
        }
    }
    fclose(fp);

}

double Corrector_KSKLKL::GetCorrectionFactorAtGeneric(double smax, double smin, double nB2KSKLKL_all, double nB2KSKLKL_NR) {

    if (nB2KSKLKL_all < N_EPSILON) return 1.0; // no correction needed
    if (nB2KSKLKL_all - nB2KSKLKL_NR > N_EPSILON) return 0.0; // remove B+ --> K+ [X --> KL0 KL0]
    if (nB2KSKLKL_all < 0 || nB2KSKLKL_NR < 0) {
        printf("[Corrector_KSKLKL] number of decay is smaller than 0!\n");
        exit(1);
    }

    // check smax and smin
    double smax_ = std::max(smax, smin);
    double smin_ = std::min(smax, smin);

    double Correction = 1;

    int GLobalBin_weight = weights_KSKLKL->FindBin(smax_, smin_);
    Correction = weights_KSKLKL->GetBinContent(GLobalBin_weight);

    return Correction;
}

class Corrector_BtoDtoXKL {
private:

    const double Nominal_correction;
    const double relative_uncertainty_correction; // relative uncertainty

    const double N_EPSILON;

public:
    Corrector_BtoDtoXKL();
    double GetCorrectionFactorAtGeneric(double nBtoDtoXKL);
    double GetRelativeUncertainty(double nBtoDtoXKL);
};

Corrector_BtoDtoXKL corrector_BtoDtoXKL;

Corrector_BtoDtoXKL::Corrector_BtoDtoXKL() :
    N_EPSILON(0.01),
    Nominal_correction(1.3),
    relative_uncertainty_correction(0.1 / 1.3)
{

}

double Corrector_BtoDtoXKL::GetCorrectionFactorAtGeneric(double nBtoDtoXKL) {
    if (nBtoDtoXKL < N_EPSILON) return 1.0; // no correction needed
    if (nBtoDtoXKL < 0) {
        printf("[Corrector_BtoDtoXKL] number of decay is smaller than 0!\n");
        exit(1);
    }

    double Correction = std::pow(Nominal_correction, nBtoDtoXKL);

    return Correction;
}

double Corrector_BtoDtoXKL::GetRelativeUncertainty(double nBtoDtoXKL) {
    if (nBtoDtoXKL < N_EPSILON) return 0.0; // no uncertainty needed
    if (nBtoDtoXKL < 0) {
        printf("[Corrector_BtoDtoXKL] number of decay is smaller than 0!\n");
        exit(1);
    }

    double RelativeUncertainty = nBtoDtoXKL * relative_uncertainty_correction;

    return RelativeUncertainty;
}

/* ====================================== */

// global variables to calculate uncertainties
std::vector<double> Ns;
std::vector<int> ntracks;
std::vector<int> npi0s;
std::vector<int> dmID;
std::vector<double> KS0_3D_distance;
std::vector<double> KaonID_correction;
std::vector<double> KaonID_rel_up;
std::vector<double> KaonID_rel_dn;

void load_files(const char* dirname, std::vector<string>* names);

void load_files(const char* dirname, std::vector<string>* names) {
    TSystemDirectory dir(dirname, dirname);
    TList* files = dir.GetListOfFiles();
    if (files) {
        TSystemFile* file;
        TString fname;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(".root")) {
                names->push_back(fname.Data());
            }
        }
    }
}

void load_files(const char* dirname, std::vector<string>* names, const char* included_string) {
    TSystemDirectory dir(dirname, dirname);
    TList* files = dir.GetListOfFiles();
    if (files) {
        TSystemFile* file;
        TString fname;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(".root") && fname.Contains(included_string)) {
                names->push_back(fname.Data());
            }
        }
    }
}

double GetEvtNum(const char* dirname, const char* included_string, double weight_var = 1.0) {

    double Total_Evt_Num = 0;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        Total_Evt_Num = Total_Evt_Num + weight_var * tree_upsilon->GetEntries();

        input_file->Close();

    }

    return Total_Evt_Num;

}

double GetNpi0(double Upsilon_ID, double Bsig_ID) {

    if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) return 0.0; // B2Kc
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) return 1.0; // B2KcPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) return 0.0; // B2Ks0Pic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) return 0.0; // B2KcPicPic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) return 1.0; // B2Ks0PicPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) return 1.0; // B2KcPicPicPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) return 0.0; // B2Ks0PicPicPic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) return 0.0; // B2KcPicPicPicPic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) return 1.0; // B2Ks0PicPicPicPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) return 2.0; // B2KcPi0Pi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) return 2.0; // B2Ks0PicPi0Pi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) return 2.0; // B2KcPicPicPi0Pi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) return 0.0; // B2KcKcKc
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) return 0.0; // B2KcKcKs0Pic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) return 1.0; // B2KcKcKcPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) return 0.0; // B02Ks0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) return 0.0; // B02KcPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) return 1.0; // B02Ks0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) return 1.0; // B02KcPicPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) return 0.0; // B02Ks0PicPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) return 0.0; // B02KcPicPicPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) return 1.0; // B02Ks0PicPicPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) return 1.0; // B02KcPicPicPicPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) return 0.0; // B02Ks0PicPicPicPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) return 2.0; // B02Ks0Pi0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) return 2.0; // B02KcPicPi0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) return 2.0; // B02Ks0PicPicPi0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) return 0.0; // B02KcKcKs0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) return 0.0; // B02KcKcKcPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) return 1.0; // B02KcKcKs0Pi0
    else {
        printf("[ERROR] unexpected decay ID\n");
        exit(1);
    }

    return -1.0;
}

double GetNtracks(double Upsilon_ID, double Bsig_ID) {

    if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) return 1.0; // B2Kc
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) return 1.0; // B2KcPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) return 1.0; // B2Ks0Pic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) return 3.0; // B2KcPicPic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) return 1.0; // B2Ks0PicPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) return 3.0; // B2KcPicPicPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) return 3.0; // B2Ks0PicPicPic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) return 5.0; // B2KcPicPicPicPic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) return 3.0; // B2Ks0PicPicPicPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) return 1.0; // B2KcPi0Pi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) return 1.0; // B2Ks0PicPi0Pi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) return 3.0; // B2KcPicPicPi0Pi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) return 3.0; // B2KcKcKc
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) return 3.0; // B2KcKcKs0Pic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) return 3.0; // B2KcKcKcPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) return 0.0; // B02Ks0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) return 2.0; // B02KcPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) return 0.0; // B02Ks0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) return 2.0; // B02KcPicPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) return 2.0; // B02Ks0PicPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) return 4.0; // B02KcPicPicPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) return 2.0; // B02Ks0PicPicPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) return 4.0; // B02KcPicPicPicPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) return 4.0; // B02Ks0PicPicPicPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) return 0.0; // B02Ks0Pi0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) return 2.0; // B02KcPicPi0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) return 2.0; // B02Ks0PicPicPi0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) return 2.0; // B02KcKcKs0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) return 4.0; // B02KcKcKcPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) return 2.0; // B02KcKcKs0Pi0
    else {
        printf("[ERROR] unexpected decay ID\n");
        exit(1);
    }

    return -1.0;
}

double BDTcToWeight(double BDTc) {

    if (BDTc > (5.0 / 6.0)) return std::sqrt(5.0);
    else return std::sqrt(BDTc / (1.0 - BDTc));

}

double GetNominalPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get nominal PDF with appropriate correction
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;

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

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

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

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

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

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if(CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double GetTrackPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, bool IsItUp, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get track efficiency uncertainty PDF with appropriate correction
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;

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

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

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

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);
            double Ntrack = GetNtracks(Upsilon_ID, Bsig_ID);

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

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double Correction_track = 1.0;
            if (IsItUp == true) Correction_track = (1 + (Ntrack * track_rel_uncertainty / 100.0));
            else Correction_track = (1 - (Ntrack * track_rel_uncertainty / 100.0));

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_BtoDtoXKL * Correction_track;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double GetKS0PDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, bool IsItUp, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get KS0 uncertainty PDF with appropriate correction
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_KS0_3D_distance = -1;

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;

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

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
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
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_3D_distance", &temp_KS0_3D_distance);
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

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

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

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

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double KS0_correction = 1.0;
            if(IsItUp == true) KS0_correction = 1 + (KS0_rel_uncertainty * temp_KS0_3D_distance / 100.0);
            else KS0_correction = 1 - (KS0_rel_uncertainty * temp_KS0_3D_distance / 100.0);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_BtoDtoXKL * KS0_correction;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

int GetFEIcorrelatedPDFs(const char* dirname, TH1D* CHG_nominal_hist, TH1D* MIX_nominal_hist, TH1D* Signal_nominal_hist, TH1D*** CHG_hists, TH1D*** MIX_hists, TH1D*** Signal_hists) { // get shape sys histogram from txt file
    int Nentry = 0; // number of eigen values/vectors
    double eigen_value = 0; // eigen value
    double weight_sys[RarityBins * 3] = { 0.0 }; // eigen vector

    FILE* fp;
    fp = fopen(dirname, "r");
    while (true) {
        if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
        for (int i = 0; i < RarityBins * 3; i++) {
            if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
        }
        Nentry++;
    }
    fclose(fp);

    *CHG_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *MIX_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *Signal_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);

    for (int i = 0; i < Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_FEI_correlated" + std::to_string(i) + "_p").c_str(), ("CHG_FEI_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_FEI_correlated" + std::to_string(i) + "_p").c_str(), ("MIX_FEI_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists)[i] = new TH1D(("Signal_FEI_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_FEI_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
    }

    for (int i = Nentry; i < 2 * Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_FEI_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHG_FEI_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_FEI_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("MIX_FEI_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists)[i] = new TH1D(("Signal_FEI_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_FEI_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
    }

    fp = fopen(dirname, "r");
    for (int i = 0; i < Nentry; i++) {
        fscanf(fp, "%lf\n", &eigen_value);
        for (int j = 0; j < RarityBins * 3; j++) fscanf(fp, "%lf\n", &weight_sys[j]);

        for (int k = 0; k < RarityBins; k++) {
            (*CHG_hists)[i]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 1 * RarityBins]));
            (*Signal_hists)[i]->SetBinContent(k + 1, Signal_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 2 * RarityBins]));

            (*CHG_hists)[i + Nentry]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i + Nentry]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 1 * RarityBins]));
            (*Signal_hists)[i + Nentry]->SetBinContent(k + 1, Signal_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 2 * RarityBins]));
        }

    }

    fclose(fp);

    return Nentry;
}

void GetFEIUncorrelatedPDFs(const char* dirname, TH1D* CHG_hist, TH1D* MIX_hist, TH1D* Signal_hist) { // get shape sys histogram from txt file
    FILE* fp;
    fp = fopen(dirname, "r");

    double weight_sys[RarityBins * 3] = { 0.0 };
    for (int i = 0; i < RarityBins * 3; i++) fscanf(fp, "%lf\n", &weight_sys[i]);
    fclose(fp);

    for (int i = 0; i < RarityBins * 3; i++) weight_sys[i] = std::sqrt(weight_sys[i]);

    for (int i = 0; i < RarityBins; i++) CHG_hist->SetBinContent(i + 1, weight_sys[i]);
    for (int i = 0; i < RarityBins; i++) MIX_hist->SetBinContent(i + 1, weight_sys[RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) Signal_hist->SetBinContent(i + 1, weight_sys[2 * RarityBins + i]);
}

int GetPIDcorrelatedPDFs(const char* dirname, TH1D* CHG_nominal_hist, TH1D* MIX_nominal_hist, TH1D* UUBAR_nominal_hist, TH1D* DDBAR_nominal_hist, TH1D* SSBAR_nominal_hist, TH1D* CHARM_nominal_hist, TH1D* SIGNAL_nominal_hist, TH1D*** CHG_hists, TH1D*** MIX_hists, TH1D*** UUBAR_hists, TH1D*** DDBAR_hists, TH1D*** SSBAR_hists, TH1D*** CHARM_hists, TH1D*** SIGNAL_hists, bool IsItKID) { // get shape sys histogram from txt file
    int Nentry = 0; // number of eigen values/vectors
    double eigen_value = 0; // eigen value
    double weight_sys[RarityBins * 7] = { 0.0 }; // eigen vector

    FILE* fp;
    fp = fopen(dirname, "r");
    while (true) {
        if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
        for (int i = 0; i < RarityBins * 7; i++) {
            if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
        }
        Nentry++;
    }
    fclose(fp);

    *CHG_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *MIX_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *UUBAR_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *DDBAR_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *SSBAR_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *CHARM_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *SIGNAL_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);

    for (int i = 0; i < Nentry; i++) {
        if (IsItKID) {
            (*CHG_hists)[i] = new TH1D(("CHG_KID_correlated" + std::to_string(i) + "_p").c_str(), ("CHG_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*MIX_hists)[i] = new TH1D(("MIX_KID_correlated" + std::to_string(i) + "_p").c_str(), ("MIX_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*UUBAR_hists)[i] = new TH1D(("UUBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), ("UUBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*DDBAR_hists)[i] = new TH1D(("DDBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), ("DDBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*SSBAR_hists)[i] = new TH1D(("SSBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), ("SSBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*CHARM_hists)[i] = new TH1D(("CHARM_KID_correlated" + std::to_string(i) + "_p").c_str(), ("CHARM_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*SIGNAL_hists)[i] = new TH1D(("Signal_KID_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        }
        else {
            (*CHG_hists)[i] = new TH1D(("CHG_PID_correlated" + std::to_string(i) + "_p").c_str(), ("CHG_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*MIX_hists)[i] = new TH1D(("MIX_PID_correlated" + std::to_string(i) + "_p").c_str(), ("MIX_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*UUBAR_hists)[i] = new TH1D(("UUBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), ("UUBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*DDBAR_hists)[i] = new TH1D(("DDBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), ("DDBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*SSBAR_hists)[i] = new TH1D(("SSBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), ("SSBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*CHARM_hists)[i] = new TH1D(("CHARM_PID_correlated" + std::to_string(i) + "_p").c_str(), ("CHARM_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*SIGNAL_hists)[i] = new TH1D(("Signal_PID_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        }
    }

    for (int i = Nentry; i < 2* Nentry; i++) {
        if (IsItKID) {
            (*CHG_hists)[i] = new TH1D(("CHG_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHG_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*MIX_hists)[i] = new TH1D(("MIX_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("MIX_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*UUBAR_hists)[i] = new TH1D(("UUBAR_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("UUBAR_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*DDBAR_hists)[i] = new TH1D(("DDBAR_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("DDBAR_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*SSBAR_hists)[i] = new TH1D(("SSBAR_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("SSBAR_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*CHARM_hists)[i] = new TH1D(("CHARM_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHARM_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*SIGNAL_hists)[i] = new TH1D(("Signal_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        }
        else {
            (*CHG_hists)[i] = new TH1D(("CHG_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHG_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*MIX_hists)[i] = new TH1D(("MIX_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("MIX_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*UUBAR_hists)[i] = new TH1D(("UUBAR_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("UUBAR_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*DDBAR_hists)[i] = new TH1D(("DDBAR_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("DDBAR_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*SSBAR_hists)[i] = new TH1D(("SSBAR_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("SSBAR_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*CHARM_hists)[i] = new TH1D(("CHARM_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHARM_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*SIGNAL_hists)[i] = new TH1D(("Signal_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        }
    }

    fp = fopen(dirname, "r");
    for (int i = 0; i < Nentry; i++) {
        fscanf(fp, "%lf\n", &eigen_value);
        for (int j = 0; j < RarityBins * 7; j++) fscanf(fp, "%lf\n", &weight_sys[j]);

        for (int k = 0; k < RarityBins; k++) {
            (*CHG_hists)[i]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 1 * RarityBins]));
            (*UUBAR_hists)[i]->SetBinContent(k + 1, UUBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 2 * RarityBins]));
            (*DDBAR_hists)[i]->SetBinContent(k + 1, DDBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 3 * RarityBins]));
            (*SSBAR_hists)[i]->SetBinContent(k + 1, SSBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 4 * RarityBins]));
            (*CHARM_hists)[i]->SetBinContent(k + 1, CHARM_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 5 * RarityBins]));
            (*SIGNAL_hists)[i]->SetBinContent(k + 1, SIGNAL_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 6 * RarityBins]));

            (*CHG_hists)[i + Nentry]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i + Nentry]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 1 * RarityBins]));
            (*UUBAR_hists)[i + Nentry]->SetBinContent(k + 1, UUBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 2 * RarityBins]));
            (*DDBAR_hists)[i + Nentry]->SetBinContent(k + 1, DDBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 3 * RarityBins]));
            (*SSBAR_hists)[i + Nentry]->SetBinContent(k + 1, SSBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 4 * RarityBins]));
            (*CHARM_hists)[i + Nentry]->SetBinContent(k + 1, CHARM_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 5 * RarityBins]));
            (*SIGNAL_hists)[i + Nentry]->SetBinContent(k + 1, SIGNAL_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 6 * RarityBins]));
        }

    }

    fclose(fp);

    return Nentry;
}

void GetPIDUncorrelatedPDFs(const char* dirname, TH1D* CHG_hist, TH1D* MIX_hist, TH1D* UUBAR_hist, TH1D* DDBAR_hist, TH1D* SSBAR_hist, TH1D* CHARM_hist, TH1D* SIGNAL_hist) { // get shape sys histogram from txt file
    FILE* fp;
    fp = fopen(dirname, "r");

    double weight_sys[RarityBins * 7] = { 0.0 };
    for (int i = 0; i < RarityBins * 7; i++) fscanf(fp, "%lf\n", &weight_sys[i]);
    fclose(fp);

    for (int i = 0; i < RarityBins * 7; i++) weight_sys[i] = std::sqrt(weight_sys[i]);

    for (int i = 0; i < RarityBins; i++) CHG_hist->SetBinContent(i + 1, weight_sys[i]);
    for (int i = 0; i < RarityBins; i++) MIX_hist->SetBinContent(i + 1, weight_sys[RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) UUBAR_hist->SetBinContent(i + 1, weight_sys[2 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) DDBAR_hist->SetBinContent(i + 1, weight_sys[3 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) SSBAR_hist->SetBinContent(i + 1, weight_sys[4 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) CHARM_hist->SetBinContent(i + 1, weight_sys[5 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) SIGNAL_hist->SetBinContent(i + 1, weight_sys[6 * RarityBins + i]);
}

int GetBRcorrelatedPDFs(const char* dirname, TH1D* CHG_nominal_hist, TH1D* MIX_nominal_hist, TH1D* Signal_nominal_hist, TH1D*** CHG_hists, TH1D*** MIX_hists, TH1D*** Signal_hists) { // get shape sys histogram from txt file
    int Nentry = 0; // number of eigen values/vectors
    double eigen_value = 0; // eigen value
    double weight_sys[RarityBins * 3] = { 0.0 }; // eigen vector

    FILE* fp;
    fp = fopen(dirname, "r");
    while (true) {
        if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
        for (int i = 0; i < RarityBins * 3; i++) {
            if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
        }
        Nentry++;
    }
    fclose(fp);

    *CHG_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *MIX_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *Signal_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);

    for (int i = 0; i < Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_BR_correlated" + std::to_string(i) + "_p").c_str(), ("CHG_BR_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_BR_correlated" + std::to_string(i) + "_p").c_str(), ("MIX_BR_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists)[i] = new TH1D(("Signal_BR_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_BR_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
    }

    for (int i = Nentry; i < 2 * Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_BR_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHG_BR_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_BR_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("MIX_BR_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists)[i] = new TH1D(("Signal_BR_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_BR_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
    }

    fp = fopen(dirname, "r");
    for (int i = 0; i < Nentry; i++) {
        fscanf(fp, "%lf\n", &eigen_value);
        for (int j = 0; j < RarityBins * 3; j++) fscanf(fp, "%lf\n", &weight_sys[j]);

        for (int k = 0; k < RarityBins; k++) {
            (*CHG_hists)[i]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 1 * RarityBins]));
            (*Signal_hists)[i]->SetBinContent(k + 1, Signal_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 2 * RarityBins]));

            (*CHG_hists)[i + Nentry]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i + Nentry]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 1 * RarityBins]));
            (*Signal_hists)[i + Nentry]->SetBinContent(k + 1, Signal_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 2 * RarityBins]));
        }

    }

    fclose(fp);

    return Nentry;
}

void GetBRUncorrelatedPDFs(const char* dirname, TH1D* CHG_hist, TH1D* MIX_hist, TH1D* Signal_hist) { // get shape sys histogram from txt file
    FILE* fp;
    fp = fopen(dirname, "r");

    double weight_sys[RarityBins * 3] = { 0.0 };
    for (int i = 0; i < RarityBins * 3; i++) fscanf(fp, "%lf\n", &weight_sys[i]);
    fclose(fp);

    for (int i = 0; i < RarityBins * 3; i++) weight_sys[i] = std::sqrt(weight_sys[i]);

    for (int i = 0; i < RarityBins; i++) CHG_hist->SetBinContent(i + 1, weight_sys[i]);
    for (int i = 0; i < RarityBins; i++) MIX_hist->SetBinContent(i + 1, weight_sys[RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) Signal_hist->SetBinContent(i + 1, weight_sys[2 * RarityBins + i]);
}

int Getpi0correlatedPDFs(const char* dirname, TH1D* CHG_nominal_hist, TH1D* MIX_nominal_hist, TH1D* UUBAR_nominal_hist, TH1D* DDBAR_nominal_hist, TH1D* SSBAR_nominal_hist, TH1D* CHARM_nominal_hist, TH1D* SIGNAL_nominal_hist, TH1D*** CHG_hists, TH1D*** MIX_hists, TH1D*** UUBAR_hists, TH1D*** DDBAR_hists, TH1D*** SSBAR_hists, TH1D*** CHARM_hists, TH1D*** SIGNAL_hists) { // get shape sys histogram from txt file
    int Nentry = 0; // number of eigen values/vectors
    double eigen_value = 0; // eigen value
    double weight_sys[RarityBins * 7] = { 0.0 }; // eigen vector

    FILE* fp;
    fp = fopen(dirname, "r");
    while (true) {
        if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
        for (int i = 0; i < RarityBins * 7; i++) {
            if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
        }
        Nentry++;
    }
    fclose(fp);

    *CHG_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *MIX_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *UUBAR_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *DDBAR_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *SSBAR_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *CHARM_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *SIGNAL_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);

    for (int i = 0; i < Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("CHG_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("MIX_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*UUBAR_hists)[i] = new TH1D(("UUBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("UUBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*DDBAR_hists)[i] = new TH1D(("DDBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("DDBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*SSBAR_hists)[i] = new TH1D(("SSBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("SSBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*CHARM_hists)[i] = new TH1D(("CHARM_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("CHARM_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*SIGNAL_hists)[i] = new TH1D(("Signal_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
    }

    for (int i = Nentry; i < 2 * Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHG_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("MIX_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*UUBAR_hists)[i] = new TH1D(("UUBAR_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("UUBAR_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*DDBAR_hists)[i] = new TH1D(("DDBAR_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("DDBAR_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*SSBAR_hists)[i] = new TH1D(("SSBAR_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("SSBAR_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*CHARM_hists)[i] = new TH1D(("CHARM_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHARM_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*SIGNAL_hists)[i] = new TH1D(("Signal_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
    }

    fp = fopen(dirname, "r");
    for (int i = 0; i < Nentry; i++) {
        fscanf(fp, "%lf\n", &eigen_value);
        for (int j = 0; j < RarityBins * 7; j++) fscanf(fp, "%lf\n", &weight_sys[j]);

        for (int k = 0; k < RarityBins; k++) {
            (*CHG_hists)[i]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 1 * RarityBins]));
            (*UUBAR_hists)[i]->SetBinContent(k + 1, UUBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 2 * RarityBins]));
            (*DDBAR_hists)[i]->SetBinContent(k + 1, DDBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 3 * RarityBins]));
            (*SSBAR_hists)[i]->SetBinContent(k + 1, SSBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 4 * RarityBins]));
            (*CHARM_hists)[i]->SetBinContent(k + 1, CHARM_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 5 * RarityBins]));
            (*SIGNAL_hists)[i]->SetBinContent(k + 1, SIGNAL_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 6 * RarityBins]));

            (*CHG_hists)[i + Nentry]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i + Nentry]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 1 * RarityBins]));
            (*UUBAR_hists)[i + Nentry]->SetBinContent(k + 1, UUBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 2 * RarityBins]));
            (*DDBAR_hists)[i + Nentry]->SetBinContent(k + 1, DDBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 3 * RarityBins]));
            (*SSBAR_hists)[i + Nentry]->SetBinContent(k + 1, SSBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 4 * RarityBins]));
            (*CHARM_hists)[i + Nentry]->SetBinContent(k + 1, CHARM_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 5 * RarityBins]));
            (*SIGNAL_hists)[i + Nentry]->SetBinContent(k + 1, SIGNAL_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 6 * RarityBins]));
        }

    }

    fclose(fp);

    return Nentry;
}

void Getpi0UncorrelatedPDFs(const char* dirname, TH1D* CHG_hist, TH1D* MIX_hist, TH1D* UUBAR_hist, TH1D* DDBAR_hist, TH1D* SSBAR_hist, TH1D* CHARM_hist, TH1D* SIGNAL_hist) { // get shape sys histogram from txt file
    FILE* fp;
    fp = fopen(dirname, "r");

    double weight_sys[RarityBins * 7] = { 0.0 };
    for (int i = 0; i < RarityBins * 7; i++) fscanf(fp, "%lf\n", &weight_sys[i]);
    fclose(fp);

    for (int i = 0; i < RarityBins * 7; i++) weight_sys[i] = std::sqrt(weight_sys[i]);

    for (int i = 0; i < RarityBins; i++) CHG_hist->SetBinContent(i + 1, weight_sys[i]);
    for (int i = 0; i < RarityBins; i++) MIX_hist->SetBinContent(i + 1, weight_sys[RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) UUBAR_hist->SetBinContent(i + 1, weight_sys[2 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) DDBAR_hist->SetBinContent(i + 1, weight_sys[3 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) SSBAR_hist->SetBinContent(i + 1, weight_sys[4 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) CHARM_hist->SetBinContent(i + 1, weight_sys[5 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) SIGNAL_hist->SetBinContent(i + 1, weight_sys[6 * RarityBins + i]);
}

void GetKffPDFs(const char* dirname, const char* included_string, TH1D* hist[7], double Correction_factor_BR[7], const char* type, int charge, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get Kff uncertainty PDF with appropriate correction
    // Be careful! You should combine this function with `GetNominalPDFs` function!
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    double m_b = -1; // B meson mass
    double m_k = -1; // Kaon mass

    const double alpha0 = 0.2545; // 0.2545 +- 0.0090
    const double alpha1 = -0.71; // -0.71 +- 0.14
    const double alpha2 = 0.32; // 0.32 +- 0.59

    /*
    <PCA>
    e0 = 0.00090 alpha0 + 0.19987 alpha1 + 0.97982 alpha2
    e1 = 0.06182 alpha0 + 0.97794 alpha1 - 0.19954 alpha2
    e2 = 0.99809 alpha0 - 0.06075 alpha2 + 0.01148 alpha2
    lambda0 = 0.60197
    lambda1 = 0.07314
    lambda2 = 0.00778
    */
    const double LinearCoefficients[3][3] = {
        {0.000898287544607,   0.199871550874459,   0.979821696141969},
        {0.061817247836003,   0.977937066768286, - 0.199543782941545},
        {0.998087080844888, - 0.060749128320321,   0.011477040509285}
    };
    const double Lambdas[7][3] = {
        { 0.0,                0.0,                0.0              },
        { 0.601972208713257,  0.0,                0.0              },
        {-0.601972208713257,  0.0,                0.0              },
        { 0.0,                0.073143483537865,  0.0              },
        { 0.0,               -0.073143483537865,  0.0              },
        { 0.0,                0.0,                0.007777700923050},
        { 0.0,                0.0,               -0.007777700923050}
    };
    double fluctuations[7][3] = { 0.0 };
    for (int i = 0; i < 7; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) fluctuations[i][j] = fluctuations[i][j] + Lambdas[i][k] * LinearCoefficients[k][j];
        }
    }
    /*
    fluctuations[0][0] = 0.0; // alpha0 change
    fluctuations[0][1] = 0.0; // alpha1 change
    fluctuations[0][2] = 0.0; // alpha2 change
    fluctuations[1][0] = 0.60197 * 0.00090; // alpha0 change Lambdas[1][0]*LinearCoefficients[0][0] + Lambdas[1][1]*LinearCoefficients[1][0] + Lambdas[1][2]*LinearCoefficients[2][0]
    fluctuations[1][1] = 0.60197 * 0.19987; // alpha1 change Lambdas[1][0]*LinearCoefficients[0][1] + Lambdas[1][1]*LinearCoefficients[1][1] + Lambdas[1][2]*LinearCoefficients[2][1]
    fluctuations[1][2] = 0.60197 * 0.97982; // alpha2 change
    fluctuations[2][0] = (-0.60197) * 0.00090; // alpha0 change
    fluctuations[2][1] = (-0.60197) * 0.19987; // alpha1 change
    fluctuations[2][2] = (-0.60197) * 0.97982; // alpha2 change
    fluctuations[3][0] = 0.07314 * 0.06182; // alpha0 change
    fluctuations[3][1] = 0.07314 * 0.97794; // alpha1 change
    fluctuations[3][2] = 0.07314 * (-0.19954); // alpha2 change
    fluctuations[4][0] // alpha0 change
    fluctuations[4][1] // alpha1 change
    fluctuations[4][2] // alpha2 change
    fluctuations[5][0] // alpha0 change
    fluctuations[5][1] // alpha1 change
    fluctuations[5][2] // alpha2 change
    fluctuations[6][0] // alpha0 change
    fluctuations[6][1] // alpha1 change
    fluctuations[6][2] // alpha2 change
    */
    double value[7] = { 0.0 }; // value of lambda^1.5 * fp*fp
    double Nevts[7] = { 0.0 }; // number of events at each fluctuations
    double tot_value[7] = { 0.0 }; // value of integral lambda^1.5 * fp*fp

    double q2 = -1;

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;

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

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
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
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

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

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        TTree* tree_Xs = (TTree*)input_file->Get("Xs");

        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &q2);
        if (charge == 0) {
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spM__bc", &m_b);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &m_k);
        }
        else {
            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spM__bc", &m_b);
            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &m_k);
        }

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            tree_Xs->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

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

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            q2 = q2 * q2;

            for (int k = 0; k < 7; k++) {
                const double alpha0_fluc = alpha0 + fluctuations[k][0];
                const double alpha1_fluc = alpha1 + fluctuations[k][1];
                const double alpha2_fluc = alpha2 + fluctuations[k][2];

                double mp = 5.4158;
                double tp = (m_b + m_k) * (m_b + m_k);
                double tm = (m_b - m_k) * (m_b - m_k);
                double t0 = 0.0;
                double z = (sqrt(tp - q2) - sqrt(tp - t0)) / (sqrt(tp - q2) + sqrt(tp - t0));
                double fp = (1 / (1 - q2 / (mp * mp))) * (alpha0_fluc + alpha1_fluc * z + alpha2_fluc * z * z + (-alpha1_fluc + 2 * alpha2_fluc) * z * z * z / 3);
                double lambda = (m_b * m_b * m_b * m_b) + (m_k * m_k * m_k * m_k) + (q2 * q2) - 2 * (m_b * m_b * m_k * m_k + m_b * m_b * q2 + m_k * m_k * q2);

                value[k] = std::pow(lambda, 1.5) * fp * fp;
                double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_BtoDtoXKL * (value[k] / value[0]);
                if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(q2, "Bplus");
                else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(q2, "Bzero");
                hist[k]->Fill(MVA_var, total_weight );
                Nevts[k] = Nevts[k] + total_weight;
            }

        }
        input_file->Close();
    }


    // calculate fluctuation of total BR
    for (int k = 0; k < 7; k++) {
        const int step = 2000;
        double m_b_fix = -1;
        double m_k_fix = -1;
        if (charge == 0) {
            m_b_fix = 5.27965;
            m_k_fix = 0.497611;
        }
        else if (charge == 1 || charge == -1) {
            m_b_fix = 5.27934;
            m_k_fix = 0.493677;
        }
        const double q2_min = 0;
        const double q2_max = (m_b_fix - m_k_fix) * (m_b_fix - m_k_fix);
        const double delq2 = (q2_max - q2_min) / step; // max: (mb-mk)^2

        const double alpha0_fluc = alpha0 + fluctuations[k][0];
        const double alpha1_fluc = alpha1 + fluctuations[k][1];
        const double alpha2_fluc = alpha2 + fluctuations[k][2];

        for (int i = 0; i < step; i++) {
            const double q2_integral = q2_min + delq2 * i;

            double mp = 5.4158;
            double tp = (m_b_fix + m_k_fix) * (m_b_fix + m_k_fix);
            double tm = (m_b_fix - m_k_fix) * (m_b_fix - m_k_fix);
            double t0 = 0.0;
            double z = (sqrt(tp - q2_integral) - sqrt(tp - t0)) / (sqrt(tp - q2_integral) + sqrt(tp - t0));
            double fp = (1 / (1 - q2_integral / (mp * mp))) * (alpha0_fluc + alpha1_fluc * z + alpha2_fluc * z * z + (-alpha1_fluc + 2 * alpha2_fluc) * z * z * z / 3);
            double lambda = (m_b_fix * m_b_fix * m_b_fix * m_b_fix) + (m_k_fix * m_k_fix * m_k_fix * m_k_fix) + (q2_integral * q2_integral) - 2 * (m_b_fix * m_b_fix * m_k_fix * m_k_fix + m_b_fix * m_b_fix * q2_integral + m_k_fix * m_k_fix * q2_integral);

            tot_value[k] = tot_value[k] + std::pow(lambda, 1.5) * fp * fp * delq2;
        }
        if (charge == 0) {
            Correction_factor_BR[k] = tot_value[k] / tot_value[0];
        }
        else if (charge == 1 || charge == -1) {
            Correction_factor_BR[k] = tot_value[k] / tot_value[0];
        }
    }


    return;
}

void GetKstarffPDFs(const char* dirname, const char* included_string, TH1D* hist[19], double Correction_factor_BR[19], const char* type, int charge, double weight_var = 1.0) { // get Kstarff uncertainty PDF with appropriate correction
    // Be careful! You should combine this function with `GetNominalPDFs` function!
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    double m_b = -1; // B meson mass
    double m_k = -1; // Kaon star mass
    double costheta = -100; // costheta

    const double alpha0_A1 = 0.3;
    const double alpha1_A1 = 0.39;
    const double alpha2_A1 = 1.19;
    const double alpha0_A12 = 0.27;
    const double alpha1_A12 = 0.53;
    const double alpha2_A12 = 0.48;
    const double alpha0_v0 = 0.38;
    const double alpha1_v0 = -1.17;
    const double alpha2_v0 = 2.42;
    const double mR_A1 = 5.829;
    const double mR_A12 = 5.829;
    const double mR_v0 = 5.415;
    /*
    <PCA>
    e0 = - 0.00477 A1a0 + 0.01377 A1a1 + 0.29659 A1a2 + 0.00314 A12a0 + 0.02574 A12a1 + 0.11862 A12a2 - 0.00662 Va0 + 0.02683 Va1 + 0.94674 Va2
    e1 = 0.00955 A1a0 + 0.14336 A1a1 + 0.83303 A1a2 + 0.00253 A12a0 + 0.04200 A12a1 + 0.40190 A12a2 + 0.01359 Va0 + 0.14322 Va1 - 0.31847 Va2
    e2 = 0.00490 A1a0 - 0.02719 A1a1 - 0.43464 A1a2 - 0.00490 A12a0 + 0.14067 A12a1 + 0.88604 A12a2 + 0.00321 Va0 + 0.07099 Va1 + 0.01977 Va2
    e = 0.A1a0 + A1a1 + A1a2 + A12a0 + A12a1 + A12a2 + Va0 + Va1 + Va2
    e = A1a0 + A1a1 + A1a2 + A12a0 + A12a1 + A12a2 + Va0 + Va1 + Va2
    e = A1a0 + A1a1 + A1a2 + A12a0 + A12a1 + A12a2 + Va0 + Va1 + Va2
    e = A1a0 + A1a1 + A1a2 + A12a0 + A12a1 + A12a2 + Va0 + Va1 + Va2
    e = A1a0 + A1a1 + A1a2 + A12a0 + A12a1 + A12a2 + Va0 + Va1 + Va2
    e = A1a0 + A1a1 + A1a2 + A12a0 + A12a1 + A12a2 + Va0 + Va1 + Va2
    lambda0 = 1.57728
    lambda1 = 1.06016
    lambda2 = 0.52088
    lambda3 = 0.23284
    lambda4 = 0.08996
    lambda5 = 0.05088
    lambda6 = 0.01842
    lambda7 = 0.00740
    lambda8 = 0.00177
    */
    const double LinearCoefficients[9][9] = {
    {-0.004767833929074,   0.013772238042338,   0.296589616507393,   0.003139670856264,   0.025743553369441,   0.118620792973715, -0.006615160985686,   0.026825078669141,   0.946739334135628},
    {0.009545743410425,   0.143364495729883,   0.833028724155009,   0.002528027950777,   0.041995533293814,   0.401904805723855,   0.013592144331223,   0.143215823688685, -0.318473920017488},
    {0.004903939104258, -0.027192276054885, -0.434636497785463, -0.004895872506074,   0.140668610588867,   0.886044604252808,   0.003208697626475,   0.070992606469561,   0.019766944252755},
    {0.081726107544427,   0.418943528790113, -0.148580744618415,   0.008936060075137,   0.016121321335068, -0.134221136343898,   0.101230277583689,   0.875232126303976,   0.033121258010682},
    {-0.047393076509321, -0.020252979416420,   0.020465460529877,   0.194484738990744,   0.966871985468078, -0.140890698347323, -0.062525020565583, -0.015984271025904, -0.015622549984120},
    {0.130012725270740,   0.880660653289362, -0.078326566909849, -0.023912474389986,   0.027560190412870,   0.018608679240451,   0.008992181636485, -0.446264603521217,   0.022087208404295},
    {0.512122358575286, -0.109791833791786,   0.012308677320526,   0.445594862213020, -0.020341462354386,   0.009165569689839,   0.721186808802633, -0.079569219353840,   0.005540910930872},
    {-0.519183885001133,   0.030474942310095,   0.001530220137484, -0.500080935622905,   0.116801392077100, -0.015739303946084,   0.680871786694491, -0.044168374846196,   0.002926015416780},
    {-0.664990460261152,   0.118522184960495, -0.020121877756987,   0.716134125194759, -0.166939357181738,   0.027983481968456,   0.042738598277210, -0.002952732689141,   0.000271198372049}
    };
    const double Lambdas[19][9] = {
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {1.577276362747197, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {-1.577276362747197, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 1.060164618950133, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, -1.060164618950133, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.520878893100167, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, -0.520878893100167, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.232840750162789, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, -0.232840750162789, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.089960782159514, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, -0.089960782159514, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.050880752840377, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, -0.050880752840377, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.018416666017616, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.018416666017616, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.007399540696482, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.007399540696482, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001770848832452},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.001770848832452}
    };
    double fluctuations[19][9] = { 0.0 };
    for (int i = 0; i < 19; i++) {
        for (int j = 0; j < 9; j++) {
            for (int k = 0; k < 9; k++) fluctuations[i][j] = fluctuations[i][j] + Lambdas[i][k] * LinearCoefficients[k][j];
        }
    }

    double value[19] = { 0.0 }; // 
    double Nevts[19] = { 0.0 }; // number of events at each fluctuations
    double tot_value[19] = { 0.0 }; // value of integral dBR/dsdcostheta

    double q2 = -1;

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;

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

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
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
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

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

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        TTree* tree_Xs = (TTree*)input_file->Get("Xs");

        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &q2);
        if (charge == 0) {
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spM__bc", &m_b);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &m_k);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", &costheta);
        }
        else {
            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spM__bc", &m_b);
            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &m_k);
            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", &costheta);
        }

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            tree_Xs->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

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

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            q2 = q2 * q2;

            for (int k = 0; k < 19; k++) {
                const double alpha0_A1_fluc = alpha0_A1 + fluctuations[k][0];
                const double alpha1_A1_fluc = alpha1_A1 + fluctuations[k][1];
                const double alpha2_A1_fluc = alpha2_A1 + fluctuations[k][2];
                const double alpha0_A12_fluc = alpha0_A12 + fluctuations[k][3];
                const double alpha1_A12_fluc = alpha1_A12 + fluctuations[k][4];
                const double alpha2_A12_fluc = alpha2_A12 + fluctuations[k][5];
                const double alpha0_v0_fluc = alpha0_v0 + fluctuations[k][6];
                const double alpha1_v0_fluc = alpha1_v0 + fluctuations[k][7];
                const double alpha2_v0_fluc = alpha2_v0 + fluctuations[k][8];

                double tp = (m_b + m_k) * (m_b + m_k);
                double tm = (m_b - m_k) * (m_b - m_k);
                double t0 = tp * (1 - sqrt(1 - tm / tp));
                double z = (sqrt(tp - q2) - sqrt(tp - t0)) / (sqrt(tp - q2) + sqrt(tp - t0));
                double z0 = (sqrt(tp) - sqrt(tp - t0)) / (sqrt(tp) + sqrt(tp - t0));

                double v0 = (1 / (1 - q2 / (mR_v0 * mR_v0))) * (alpha0_v0_fluc + alpha1_v0_fluc * (z - z0) + alpha2_v0_fluc * (z - z0) * (z - z0));
                double A1 = (1 / (1 - q2 / (mR_A1 * mR_A1))) * (alpha0_A1_fluc + alpha1_A1_fluc * (z - z0) + alpha2_A1_fluc * (z - z0) * (z - z0));
                double A12 = (1 / (1 - q2 / (mR_A12 * mR_A12))) * (alpha0_A12_fluc + alpha1_A12_fluc * (z - z0) + alpha2_A12_fluc * (z - z0) * (z - z0));
                double lambda = (tp - q2) * (tm - q2);
                double A2 = ((m_b + m_k) * (m_b + m_k) * (m_b * m_b - m_k * m_k - q2) * A1 - A12 * 16 * m_b * m_k * m_k * (m_b + m_k)) / lambda;

                double sB = q2 / (m_b * m_b);
                double m_k_tilda = m_k / m_b;
                double Lambda = 1 + std::pow(m_k_tilda, 4) + sB * sB - 2 * (m_k_tilda * m_k_tilda + sB + sB * m_k_tilda * m_k_tilda);

                double Amp_parallel = -2 * (std::sqrt(sB)) * (std::pow(Lambda, 1.0 / 4.0)) * (std::sqrt(2)) * (1 + m_k_tilda) * A1;
                double Amp_vertical = 2 * (std::sqrt(sB)) * (std::pow(Lambda, 1.0 / 4.0)) * (std::sqrt(2)) * (std::sqrt(Lambda)) * v0 / (1 + m_k_tilda);
                double Amp_0 = -1 * (std::sqrt(sB)) * (std::pow(Lambda, 1.0 / 4.0)) * (1.0 / m_k_tilda) * (1.0 / std::pow(sB, 0.5)) * ((1 - m_k_tilda * m_k_tilda - sB) * (1 + m_k_tilda) * A1 - Lambda * A2 / (1 + m_k_tilda));

                value[k] = (3.0 / 4.0) * (Amp_vertical * Amp_vertical + Amp_parallel * Amp_parallel) * (1 - costheta * costheta) + (3.0 / 2.0) * Amp_0 * Amp_0 * costheta * costheta;
                double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_BtoDtoXKL * (value[k] / value[0]);
                if (q2 < MyEPSILON) total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_BtoDtoXKL; // Makeshift
                hist[k]->Fill(MVA_var, total_weight);
                Nevts[k] = Nevts[k] + total_weight;
            }
        }
        input_file->Close();
    }


    // calculate fluctuation of total BR
    for (int k = 0; k < 19; k++) {
        const int step = 100;
        double m_b_fix = -1;
        double m_k_fix = -1;
        if (charge == 0) {
            m_b_fix = 5.27965;
            m_k_fix = 0.89555;
        }
        else if (charge == 1 || charge == -1) {
            m_b_fix = 5.27934;
            m_k_fix = 0.89166;
        }
        const double q2_min = 0;
        const double q2_max = (m_b_fix - m_k_fix) * (m_b_fix - m_k_fix);
        const double delq2 = (q2_max - q2_min) / step; // max: (mb-mk)^2

        const double costheta_min = -1;
        const double costheta_max = 1;
        const double delcostheta = (costheta_max - costheta_min) / step;

        const double alpha0_A1_fluc = alpha0_A1 + fluctuations[k][0];
        const double alpha1_A1_fluc = alpha1_A1 + fluctuations[k][1];
        const double alpha2_A1_fluc = alpha2_A1 + fluctuations[k][2];
        const double alpha0_A12_fluc = alpha0_A12 + fluctuations[k][3];
        const double alpha1_A12_fluc = alpha1_A12 + fluctuations[k][4];
        const double alpha2_A12_fluc = alpha2_A12 + fluctuations[k][5];
        const double alpha0_v0_fluc = alpha0_v0 + fluctuations[k][6];
        const double alpha1_v0_fluc = alpha1_v0 + fluctuations[k][7];
        const double alpha2_v0_fluc = alpha2_v0 + fluctuations[k][8];

        for (int i = 0; i < step; i++) {
            const double q2_integral = q2_min + delq2 * i;

            for (int j = 0; j < step; j++) {
                const double costheta_integral = costheta_min + delcostheta * j;

                double tp = (m_b_fix + m_k_fix) * (m_b_fix + m_k_fix);
                double tm = (m_b_fix - m_k_fix) * (m_b_fix - m_k_fix);
                double t0 = tp * (1 - sqrt(1 - tm / tp));
                double z = (sqrt(tp - q2_integral) - sqrt(tp - t0)) / (sqrt(tp - q2_integral) + sqrt(tp - t0));
                double z0 = (sqrt(tp) - sqrt(tp - t0)) / (sqrt(tp) + sqrt(tp - t0));

                double v0 = (1 / (1 - q2_integral / (mR_v0 * mR_v0))) * (alpha0_v0_fluc + alpha1_v0_fluc * (z - z0) + alpha2_v0_fluc * (z - z0) * (z - z0));
                double A1 = (1 / (1 - q2_integral / (mR_A1 * mR_A1))) * (alpha0_A1_fluc + alpha1_A1_fluc * (z - z0) + alpha2_A1_fluc * (z - z0) * (z - z0));
                double A12 = (1 / (1 - q2_integral / (mR_A12 * mR_A12))) * (alpha0_A12_fluc + alpha1_A12_fluc * (z - z0) + alpha2_A12_fluc * (z - z0) * (z - z0));
                double lambda = (tp - q2_integral) * (tm - q2_integral);
                double A2 = ((m_b_fix + m_k_fix) * (m_b_fix + m_k_fix) * (m_b_fix * m_b_fix - m_k_fix * m_k_fix - q2_integral) * A1 - A12 * 16 * m_b_fix * m_k_fix * m_k_fix * (m_b_fix + m_k_fix)) / lambda;

                double sB = q2_integral / (m_b_fix * m_b_fix);
                double m_k_fix_tilda = m_k_fix / m_b_fix;
                double Lambda = 1 + std::pow(m_k_fix_tilda, 4) + sB * sB - 2 * (m_k_fix_tilda * m_k_fix_tilda + sB + sB * m_k_fix_tilda * m_k_fix_tilda);

                double Amp_parallel = -2 * (std::sqrt(sB)) * (std::pow(Lambda, 1.0 / 4.0)) * (std::sqrt(2)) * (1 + m_k_fix_tilda) * A1;
                double Amp_vertical = 2 * (std::sqrt(sB)) * (std::pow(Lambda, 1.0 / 4.0)) * (std::sqrt(2)) * (std::sqrt(Lambda)) * v0 / (1 + m_k_fix_tilda);
                double Amp_0 = -1 * (std::pow(Lambda, 1.0 / 4.0)) * (1.0 / m_k_fix_tilda) * ((1 - m_k_fix_tilda * m_k_fix_tilda - sB) * (1 + m_k_fix_tilda) * A1 - Lambda * A2 / (1 + m_k_fix_tilda));

                tot_value[k] = tot_value[k] + ((3.0 / 4.0) * (Amp_vertical * Amp_vertical + Amp_parallel * Amp_parallel) * (1 - costheta_integral * costheta_integral) + (3.0 / 2.0) * Amp_0 * Amp_0 * costheta_integral * costheta_integral) * delcostheta * (delq2 / (m_b_fix * m_b_fix));
            }
        }
        if (charge == 0) {
            Correction_factor_BR[k] = tot_value[k] / tot_value[0];
        }
        else if (charge == 1 || charge == -1) {
            Correction_factor_BR[k] = tot_value[k] / tot_value[0];
        }
    }


    return Nevt;
}

double GetFragmentationPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, int TargetMxsBin, int TargetCategory, bool IsItUp, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get fragmentation PDF with with appropriate correction
    // Be careful! You should combine this function with `GetNominalPDFs` function!
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;

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

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

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

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

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

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE) * corrector_Fragmentation.FluctuateCorrection(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, TargetMxsBin, TargetCategory, IsItUp, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE) * corrector_Fragmentation.FluctuateCorrection(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, TargetMxsBin, TargetCategory, IsItUp, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double GetNevtWithBDTc(const char* dirname, const char* included_string, const char* type, const char* sample, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get the number of event after BDTc correction
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    float MVA_var = 0;
    float BDTc_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;

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

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

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

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        tree_upsilon->SetBranchAddress("MVA_Continuum", &BDTc_var);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

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

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double BDTc_weight = BDTcToWeight(BDTc_var);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_BtoDtoXKL * BDTc_weight;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            Nevt = Nevt + total_weight;

        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double GetBDTcPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get BDTc PDF with appropriate correction
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    float MVA_var = 0;
    float BDTc_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;

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

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

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

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        tree_upsilon->SetBranchAddress("MVA_Continuum", &BDTc_var);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // apply BDTc weight with appropritate norm factor
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

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

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double BDTc_weight = BDTcToWeight(BDTc_var);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_BtoDtoXKL * BDTc_weight;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }

        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

void GetNegativeChangePDFs(TH1D* nominal_hist, TH1D* positive_hist, TH1D* negative_hist) {
    for (int i = 0; i < RarityBins; i++) {
        double nominal_value = nominal_hist->GetBinContent(i + 1);
        double positive_value = positive_hist->GetBinContent(i + 1);
        double deviation = positive_value - nominal_value;

        negative_hist->SetBinContent(i + 1, nominal_value - deviation);
    }
}

double GetmbPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, bool IsIthigh, bool IsItXsu, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // apply correction for mb weight files
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {
        printf("[ERROR] unexpected type name\n");
        printf("[ERROR] Continuum type cannot be selected for mb pdf\n");
        exit(1);
    }
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;

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

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

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

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

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

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            if (IsIthigh == true && IsItXsu == true) total_weight = total_weight * ReadWeightHist(Xsu_Hmb_weight, invM * invM);
            else if (IsIthigh == false && IsItXsu == true) total_weight = total_weight * ReadWeightHist(Xsu_Lmb_weight, invM * invM);
            else if (IsIthigh == true && IsItXsu == false) total_weight = total_weight * ReadWeightHist(Xsd_Hmb_weight, invM * invM);
            else if (IsIthigh == false && IsItXsu == false) total_weight = total_weight * ReadWeightHist(Xsd_Lmb_weight, invM * invM);

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double GetpfPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, bool IsIthigh, bool IsItXsu, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // apply correction for pf weight files
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {
        printf("[ERROR] unexpected type name\n");
        printf("[ERROR] Continuum type cannot be selected for pf pdf\n");
        exit(1);
    }
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;
    double MXs = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;

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

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            if (IsItXsu) tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &MXs);
            else tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &MXs);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

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

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

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

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, MXs, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, MXs, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            if (IsIthigh == true && IsItXsu == true) total_weight = total_weight * ReadWeightHist(Xsu_Hpf_weight, MXs);
            else if (IsIthigh == false && IsItXsu == true) total_weight = total_weight * ReadWeightHist(Xsu_Lpf_weight, MXs);
            else if (IsIthigh == true && IsItXsu == false) total_weight = total_weight * ReadWeightHist(Xsd_Hpf_weight, MXs);
            else if (IsIthigh == false && IsItXsu == false) total_weight = total_weight * ReadWeightHist(Xsd_Lpf_weight, MXs);

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double GetTransitionPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, bool IsIthigh, bool IsItXsu, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // apply correction for transition weight files
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;

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

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

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

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

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

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            if (IsIthigh == true && IsItXsu == true) total_weight = total_weight * ReadWeightHist(Xsu_Htransition_weight, MVA_var);
            else if (IsIthigh == false && IsItXsu == true) total_weight = total_weight * ReadWeightHist(Xsu_Ltransition_weight, MVA_var);
            else if (IsIthigh == true && IsItXsu == false) total_weight = total_weight * ReadWeightHist(Xsd_Htransition_weight, MVA_var);
            else if (IsIthigh == false && IsItXsu == false) total_weight = total_weight * ReadWeightHist(Xsd_Ltransition_weight, MVA_var);

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double GetKstardeltaPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, bool IsItXsu, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // apply correction for pf weight files
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    otherwise
    */
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {
        printf("[ERROR] unexpected type name\n");
        printf("[ERROR] Continuum type cannot be selected for pf pdf\n");
        exit(1);
    }
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    double invM = -1.0;
    double MXs = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;

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

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");
        TTree* tree_Xs = (TTree*)input_file->Get("Xs");

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu")) tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

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

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        if (IsItXsu) tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &MXs);
        else tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &MXs);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            tree_Xs->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

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

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");

            if (IsItXsu == true) total_weight = total_weight * ReadWeightHist(Kstar_delta_weight, MXs);
            else total_weight = total_weight * ReadWeightHist(K0star_delta_weight, MXs);

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

void GetMultiplicityPDFs(const char* dirname, const char* included_string, TH1D** hist, const char* type, const char* sample, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get nominal PDF with uncertainties on multiplicity correction
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;

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

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    int N_hist = -1;

    if (N_correction_multiplicity_weightfile > N_correction_multiplicity_MCsample) N_hist = N_correction_multiplicity_MCsample + 1;
    else N_hist = N_correction_multiplicity_weightfile + 1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

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

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

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

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            for (int k = 0; k < N_hist; k++) {
                double Correction_multiplicity_uncertainty = 1.0;
                if (std::abs(k - Ngamma_v200) < MyEPSILON) Correction_multiplicity_uncertainty = weights_Ngamma_uncer->GetBinContent(k + 1);
                hist[k]->Fill(MVA_var, total_weight * Correction_multiplicity_uncertainty);
            }
            
        }
        input_file->Close();

    }

}

double GetKnnBRPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, int KnnType, bool IsItUp, double weight_var = 1.0, std::string CorrectionType = "otherwise") {  // get nominal PDF with uncertainties on B->Knn decays
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
    /*
    * KnnType
    * 1: B+ -> K+ n nbar
    * 2: B+ -> Kstar+ n nbar
    * 3: B0 -> K0 n nbar
    * 4: B0 -> Kstar0 n nbar
    */
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;

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

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

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

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

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

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);
            if (KnnType == 1) { // K+ n nbar
                if (IsItUp) Correction_Knn = Correction_Knn * std::pow(1 + B2Knn_up_uncer[0], N_Knn);
                else Correction_Knn = Correction_Knn * std::pow(1 - B2Knn_dn_uncer[0], N_Knn);
            }
            else if (KnnType == 2) { // K*+ n nbar
                if (IsItUp) Correction_Knn = Correction_Knn * std::pow(1 + B2Knn_up_uncer[1], N_Kstarnn);
                else Correction_Knn = Correction_Knn * std::pow(1 - B2Knn_dn_uncer[1], N_Kstarnn);
            }
            else if (KnnType == 3) { // K0 n nbar
                if (IsItUp) Correction_Knn = Correction_Knn * std::pow(1 + B2Knn_up_uncer[2], N_K0nn);
                else Correction_Knn = Correction_Knn * std::pow(1 - B2Knn_dn_uncer[2], N_K0nn);
            }
            else if (KnnType == 4) { // K*0 n nbar
                if (IsItUp) Correction_Knn = Correction_Knn * std::pow(1 + B2Knn_up_uncer[3], N_K0starnn);
                else Correction_Knn = Correction_Knn * std::pow(1 - B2Knn_dn_uncer[3], N_K0starnn);
            }

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double GetBtoDtoXKLPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, bool IsItUp, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get nominal PDF with uncertainties on B-> [D -> X KL0] anything decays
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;

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

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

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

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

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

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) {
                Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);
                if (IsItUp) Correction_BtoDtoXKL = Correction_BtoDtoXKL * (1 + corrector_BtoDtoXKL.GetRelativeUncertainty(nDptoXKL + nD0toXKL));
                else Correction_BtoDtoXKL = Correction_BtoDtoXKL * (1 - corrector_BtoDtoXKL.GetRelativeUncertainty(nDptoXKL + nD0toXKL));
            }

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double ReadWeightHist(TH1D* hist, double value) {
    int Bin = hist->FindBin(value);
    if (Bin < 1) Bin = 1;
    else if (Bin > hist->GetNbinsX()) Bin = hist->GetNbinsX();
    return hist->GetBinContent(Bin);
}

void ReadSignalModelingFile() {

    FILE* fp;
    int Nentry = 0;
    double RangeMIN = -1;
    double RangeMAX = -1;

    // Xsu_Hmb
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/signal_modeling/Xsu_Hmb_weight.txt","r");
    fscanf(fp, "%d %lf %lf\n", &Nentry, &RangeMIN, &RangeMAX);
    Xsu_Hmb_weight = new TH1D("Xsu_Hmb_weight", ";;", Nentry, RangeMIN, RangeMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsu_Hmb_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Xsu_Lmb_weight
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/signal_modeling/Xsu_Lmb_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &Nentry, &RangeMIN, &RangeMAX);
    Xsu_Lmb_weight = new TH1D("Xsu_Lmb_weight", ";;", Nentry, RangeMIN, RangeMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsu_Lmb_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Xsu_Hpf_weight
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/signal_modeling/Xsu_Hpf_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &Nentry, &RangeMIN, &RangeMAX);
    Xsu_Hpf_weight = new TH1D("Xsu_Hpf_weight", ";;", Nentry, RangeMIN, RangeMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsu_Hpf_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Xsu_Lpf_weight
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/signal_modeling/Xsu_Lpf_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &Nentry, &RangeMIN, &RangeMAX);
    Xsu_Lpf_weight = new TH1D("Xsu_Lpf_weight", ";;", Nentry, RangeMIN, RangeMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsu_Lpf_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Xsu_Htransition_weight
    fp = fopen("Xsu_Htransition_weight.txt", "r");
    fscanf(fp, "%d\n", &Nentry);
    if (Nentry != RarityBins) {
        printf("[ReadSignalModelingFile] Nentry and RarityBins are not the same!\n");
        exit(1);
    }
    Xsu_Htransition_weight = new TH1D("Xsu_Htransition_weight", ";;", Nentry, BinMIN, BinMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsu_Htransition_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Xsu_Ltransition_weight
    fp = fopen("Xsu_Ltransition_weight.txt", "r");
    fscanf(fp, "%d\n", &Nentry);
    if (Nentry != RarityBins) {
        printf("[ReadSignalModelingFile] Nentry and RarityBins are not the same!\n");
        exit(1);
    }
    Xsu_Ltransition_weight = new TH1D("Xsu_Ltransition_weight", ";;", Nentry, BinMIN, BinMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsu_Ltransition_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Xsd_Hmb
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/signal_modeling/Xsd_Hmb_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &Nentry, &RangeMIN, &RangeMAX);
    Xsd_Hmb_weight = new TH1D("Xsd_Hmb_weight", ";;", Nentry, RangeMIN, RangeMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsd_Hmb_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Xsd_Lmb_weight
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/signal_modeling/Xsd_Lmb_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &Nentry, &RangeMIN, &RangeMAX);
    Xsd_Lmb_weight = new TH1D("Xsd_Lmb_weight", ";;", Nentry, RangeMIN, RangeMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsd_Lmb_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Xsd_Hpf_weight
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/signal_modeling/Xsd_Hpf_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &Nentry, &RangeMIN, &RangeMAX);
    Xsd_Hpf_weight = new TH1D("Xsd_Hpf_weight", ";;", Nentry, RangeMIN, RangeMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsd_Hpf_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Xsd_Lpf_weight
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/signal_modeling/Xsd_Lpf_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &Nentry, &RangeMIN, &RangeMAX);
    Xsd_Lpf_weight = new TH1D("Xsd_Lpf_weight", ";;", Nentry, RangeMIN, RangeMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsd_Lpf_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Xsd_Htransition_weight
    fp = fopen("Xsd_Htransition_weight.txt", "r");
    fscanf(fp, "%d\n", &Nentry);
    if (Nentry != RarityBins) {
        printf("[ReadSignalModelingFile] Nentry and RarityBins are not the same!\n");
        exit(1);
    }
    Xsd_Htransition_weight = new TH1D("Xsd_Htransition_weight", ";;", Nentry, BinMIN, BinMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsd_Htransition_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Xsd_Ltransition_weight
    fp = fopen("Xsd_Ltransition_weight.txt", "r");
    fscanf(fp, "%d\n", &Nentry);
    if (Nentry != RarityBins) {
        printf("[ReadSignalModelingFile] Nentry and RarityBins are not the same!\n");
        exit(1);
    }
    Xsd_Ltransition_weight = new TH1D("Xsd_Ltransition_weight", ";;", Nentry, BinMIN, BinMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsd_Ltransition_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Kstar_delta_weight
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/signal_modeling/MKstar_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &Nentry, &RangeMIN, &RangeMAX);
    Kstar_delta_weight = new TH1D("Kstar_delta_weight", ";;", Nentry, RangeMIN, RangeMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Kstar_delta_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // K0star_delta_weight
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/signal_modeling/MK0star_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &Nentry, &RangeMIN, &RangeMAX);
    K0star_delta_weight = new TH1D("K0star_delta_weight", ";;", Nentry, RangeMIN, RangeMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        K0star_delta_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

}

void ReadMultiplicityFile() {

    FILE* fp;
    int NgammaMAX = 0;
    const double CUTOFF = 50;

    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/multiplicity/multiplicity_weight_uncertainty.txt", "r");
    fscanf(fp, "%d\n", &NgammaMAX);
    weights_Ngamma_uncer = new TH1D("weights_Ngamma_uncer", ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);
    for (int i = 0; i < NgammaMAX + 1; i++) {
        double temp1;
        double temp2;
        double temp3;
        fscanf(fp, "%lf %lf %lf\n", &temp1, &temp2, &temp3);
        if (temp3 < CUTOFF) weights_Ngamma_uncer->SetBinContent(i + 1, temp3);
        else weights_Ngamma_uncer->SetBinContent(i + 1, CUTOFF);

        if(temp2 > MyEPSILON) N_correction_multiplicity_weightfile = i;
    }

    fclose(fp);

}

void CheckMaxNgamma(const char* dirname, const char* included_string) {

    double Ngamma_v200 = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");

        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);

            if (N_correction_multiplicity_MCsample < std::lround(Ngamma_v200)) N_correction_multiplicity_MCsample = std::lround(Ngamma_v200);

        }
        input_file->Close();


    }

}

void AllocateMultiplicityPDFs(TH1D*** hist, const char* prefix, bool IsItPositive) {
    // prefix: Signal CHG MIX UUBAR DDBAR SSBAR CHARM

    int N_hist = -1;

    if (N_correction_multiplicity_weightfile > N_correction_multiplicity_MCsample) N_hist = N_correction_multiplicity_MCsample + 1;
    else N_hist = N_correction_multiplicity_weightfile + 1;

    *hist = (TH1D**)malloc(sizeof(TH1D*) * N_hist);

    for (int i = 0; i < N_hist; i++) {
        if(IsItPositive) (*hist)[i] = new TH1D((std::string(prefix) + "_multiplicity" + std::to_string(i) + "_p").c_str(), (std::string(prefix) + "_multiplicity" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        else (*hist)[i] = new TH1D((std::string(prefix) + "_multiplicity" + std::to_string(i) + "_m").c_str(), (std::string(prefix) + "_multiplicity" + std::to_string(i) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
    }

    return N_hist;
}

int PrintInfoMultiplicity() {

    int N_hist = -1;

    if (N_correction_multiplicity_weightfile > N_correction_multiplicity_MCsample) N_hist = N_correction_multiplicity_MCsample + 1;
    else N_hist = N_correction_multiplicity_weightfile + 1;

    FILE* fp;
    fp = fopen("Multiplicity_info.txt", "w");
    fprintf(fp, "%d", N_hist);
    fclose(fp);

    return N_hist;
}

void ClearHist(TH1D* Signal_hist, TH1D* CHG_hist, TH1D* MIX_hist, TH1D* UUBAR_hist, TH1D* DDBAR_hist, TH1D* SSBAR_hist, TH1D* CHARM_hist) {
    Signal_hist->Reset();
    CHG_hist->Reset();
    MIX_hist->Reset();
    UUBAR_hist->Reset();
    DDBAR_hist->Reset();
    SSBAR_hist->Reset();
    CHARM_hist->Reset();
}

void AddSQRTHist(TH1D* output_hist, TH1D* input_hist, int Nbin) {
    for (int i = 0; i < Nbin; i++) {
        double temp_1 = output_hist->GetBinContent(i + 1);
        double temp_2 = input_hist->GetBinContent(i + 1);
        double temp_3 = std::sqrt(temp_1 * temp_1 + temp_2 * temp_2);
        output_hist->SetBinContent(i + 1, temp_3);
    }
}

void Signal_yield_fit_BDT_Rarity_HistFactory()
{

    ReadSignalModelingFile();
    ReadMultiplicityFile();



    /* ====================================== */
    // Define PDFs for HistFactory
    // nominal PDFs
    TH1D* Signal_nominal = new TH1D("Signal_nominal", "Signal_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_nominal = new TH1D("CHG_nominal", "CHG_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_nominal = new TH1D("MIX_nominal", "MIX_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_nominal = new TH1D("UUBAR_nominal", "UUBAR_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_nominal = new TH1D("DDBAR_nominal", "DDBAR_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_nominal = new TH1D("SSBAR_nominal", "SSBAR_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_nominal = new TH1D("CHARM_nominal", "CHARM_nominal", RarityBins, BinMIN, BinMAX);

    // track uncertainty
    TH1D* Signal_track_p = new TH1D("Signal_track_p", "Signal_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_track_p = new TH1D("CHG_track_p", "CHG_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_track_p = new TH1D("MIX_track_p", "MIX_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_track_p = new TH1D("UUBAR_track_p", "UUBAR_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_track_p = new TH1D("DDBAR_track_p", "DDBAR_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_track_p = new TH1D("SSBAR_track_p", "SSBAR_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_track_p = new TH1D("CHARM_track_p", "CHARM_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_track_m = new TH1D("Signal_track_m", "Signal_track_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_track_m = new TH1D("CHG_track_m", "CHG_track_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_track_m = new TH1D("MIX_track_m", "MIX_track_m", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_track_m = new TH1D("UUBAR_track_m", "UUBAR_track_m", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_track_m = new TH1D("DDBAR_track_m", "DDBAR_track_m", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_track_m = new TH1D("SSBAR_track_m", "SSBAR_track_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_track_m = new TH1D("CHARM_track_m", "CHARM_track_m", RarityBins, BinMIN, BinMAX);

    // KS0 uncertainty
    TH1D* Signal_KS0_p = new TH1D("Signal_KS0_p", "Signal_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_KS0_p = new TH1D("CHG_KS0_p", "CHG_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_KS0_p = new TH1D("MIX_KS0_p", "MIX_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_KS0_p = new TH1D("UUBAR_KS0_p", "UUBAR_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_KS0_p = new TH1D("DDBAR_KS0_p", "DDBAR_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_KS0_p = new TH1D("SSBAR_KS0_p", "SSBAR_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_KS0_p = new TH1D("CHARM_KS0_p", "CHARM_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_KS0_m = new TH1D("Signal_KS0_m", "Signal_KS0_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_KS0_m = new TH1D("CHG_KS0_m", "CHG_KS0_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_KS0_m = new TH1D("MIX_KS0_m", "MIX_KS0_m", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_KS0_m = new TH1D("UUBAR_KS0_m", "UUBAR_KS0_m", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_KS0_m = new TH1D("DDBAR_KS0_m", "DDBAR_KS0_m", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_KS0_m = new TH1D("SSBAR_KS0_m", "SSBAR_KS0_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_KS0_m = new TH1D("CHARM_KS0_m", "CHARM_KS0_m", RarityBins, BinMIN, BinMAX);

    // FEI uncertainty (correlated)
    TH1D** Signal_FEI_correlated;
    TH1D** CHG_FEI_correlated;
    TH1D** MIX_FEI_correlated;

    // FEI uncertainty (uncorrelated)
    TH1D* Signal_FEI_uncorrelated = new TH1D("Signal_FEI_uncorrelated", "Signal_FEI_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_FEI_uncorrelated = new TH1D("CHG_FEI_uncorrelated", "CHG_FEI_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_FEI_uncorrelated = new TH1D("MIX_FEI_uncorrelated", "MIX_FEI_uncorrelated", RarityBins, BinMIN, BinMAX);

    // Kaon PID uncertainty (correlated)
    TH1D** Signal_KID_correlated;
    TH1D** CHG_KID_correlated;
    TH1D** MIX_KID_correlated;
    TH1D** UUBAR_KID_correlated;
    TH1D** DDBAR_KID_correlated;
    TH1D** SSBAR_KID_correlated;
    TH1D** CHARM_KID_correlated;

    // Kaon PID uncertainty (uncorrelated)
    TH1D* Signal_KID_uncorrelated = new TH1D("Signal_KID_uncorrelated", "Signal_KID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_KID_uncorrelated = new TH1D("CHG_KID_uncorrelated", "CHG_KID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_KID_uncorrelated = new TH1D("MIX_KID_uncorrelated", "MIX_KID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_KID_uncorrelated = new TH1D("UUBAR_KID_uncorrelated", "UUBAR_KID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_KID_uncorrelated = new TH1D("DDBAR_KID_uncorrelated", "DDBAR_KID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_KID_uncorrelated = new TH1D("SSBAR_KID_uncorrelated", "SSBAR_KID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_KID_uncorrelated = new TH1D("CHARM_KID_uncorrelated", "CHARM_KID_uncorrelated", RarityBins, BinMIN, BinMAX);

    // Pion PID uncertainty (correlated)
    TH1D** Signal_PID_correlated;
    TH1D** CHG_PID_correlated;
    TH1D** MIX_PID_correlated;
    TH1D** UUBAR_PID_correlated;
    TH1D** DDBAR_PID_correlated;
    TH1D** SSBAR_PID_correlated;
    TH1D** CHARM_PID_correlated;

    // Pion PID uncertainty (uncorrelated)
    TH1D* Signal_PID_uncorrelated = new TH1D("Signal_PID_uncorrelated", "Signal_PID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_PID_uncorrelated = new TH1D("CHG_PID_uncorrelated", "CHG_PID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_PID_uncorrelated = new TH1D("MIX_PID_uncorrelated", "MIX_PID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_PID_uncorrelated = new TH1D("UUBAR_PID_uncorrelated", "UUBAR_PID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_PID_uncorrelated = new TH1D("DDBAR_PID_uncorrelated", "DDBAR_PID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_PID_uncorrelated = new TH1D("SSBAR_PID_uncorrelated", "SSBAR_PID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_PID_uncorrelated = new TH1D("CHARM_PID_uncorrelated", "CHARM_PID_uncorrelated", RarityBins, BinMIN, BinMAX);

    // BB BR uncertainty(correlated)
    TH1D** Signal_BR_correlated;
    TH1D** CHG_BR_correlated;
    TH1D** MIX_BR_correlated;

    // BB BR uncertainty (uncorrelated)
    TH1D* Signal_BR_uncorrelated = new TH1D("Signal_BR_uncorrelated", "Signal_BR_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_BR_uncorrelated = new TH1D("CHG_BR_uncorrelated", "CHG_BR_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_BR_uncorrelated = new TH1D("MIX_BR_uncorrelated", "MIX_BR_uncorrelated", RarityBins, BinMIN, BinMAX);

    // pi0 uncertainty (correlated)
    TH1D** Signal_pi0_correlated;
    TH1D** CHG_pi0_correlated;
    TH1D** MIX_pi0_correlated;
    TH1D** UUBAR_pi0_correlated;
    TH1D** DDBAR_pi0_correlated;
    TH1D** SSBAR_pi0_correlated;
    TH1D** CHARM_pi0_correlated;

    // pi0 uncertainty (uncorrelated)
    TH1D* Signal_pi0_uncorrelated = new TH1D("Signal_pi0_uncorrelated", "Signal_pi0_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_pi0_uncorrelated = new TH1D("CHG_pi0_uncorrelated", "CHG_pi0_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_pi0_uncorrelated = new TH1D("MIX_pi0_uncorrelated", "MIX_pi0_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_pi0_uncorrelated = new TH1D("UUBAR_pi0_uncorrelated", "UUBAR_pi0_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_pi0_uncorrelated = new TH1D("DDBAR_pi0_uncorrelated", "DDBAR_pi0_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_pi0_uncorrelated = new TH1D("SSBAR_pi0_uncorrelated", "SSBAR_pi0_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_pi0_uncorrelated = new TH1D("CHARM_pi0_uncorrelated", "CHARM_pi0_uncorrelated", RarityBins, BinMIN, BinMAX);

    // K nu nubar form factor
    TH1D* Signal_Kff1_p = new TH1D("Signal_Kff1_p", "Signal_Kff1_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kff1_m = new TH1D("Signal_Kff1_m", "Signal_Kff1_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kff2_p = new TH1D("Signal_Kff2_p", "Signal_Kff2_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kff2_m = new TH1D("Signal_Kff2_m", "Signal_Kff2_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kff3_p = new TH1D("Signal_Kff3_p", "Signal_Kff3_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kff3_m = new TH1D("Signal_Kff3_m", "Signal_Kff3_m", RarityBins, BinMIN, BinMAX);

    // Kstar nu nubar form factor
    TH1D* Signal_Kstarff1_p = new TH1D("Signal_Kstarff1_p", "Signal_Kstarff1_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff1_m = new TH1D("Signal_Kstarff1_m", "Signal_Kstarff1_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff2_p = new TH1D("Signal_Kstarff2_p", "Signal_Kstarff2_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff2_m = new TH1D("Signal_Kstarff2_m", "Signal_Kstarff2_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff3_p = new TH1D("Signal_Kstarff3_p", "Signal_Kstarff3_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff3_m = new TH1D("Signal_Kstarff3_m", "Signal_Kstarff3_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff4_p = new TH1D("Signal_Kstarff4_p", "Signal_Kstarff4_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff4_m = new TH1D("Signal_Kstarff4_m", "Signal_Kstarff4_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff5_p = new TH1D("Signal_Kstarff5_p", "Signal_Kstarff5_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff5_m = new TH1D("Signal_Kstarff5_m", "Signal_Kstarff5_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff6_p = new TH1D("Signal_Kstarff6_p", "Signal_Kstarff6_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff6_m = new TH1D("Signal_Kstarff6_m", "Signal_Kstarff6_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff7_p = new TH1D("Signal_Kstarff7_p", "Signal_Kstarff7_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff7_m = new TH1D("Signal_Kstarff7_m", "Signal_Kstarff7_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff8_p = new TH1D("Signal_Kstarff8_p", "Signal_Kstarff8_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff8_m = new TH1D("Signal_Kstarff8_m", "Signal_Kstarff8_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff9_p = new TH1D("Signal_Kstarff9_p", "Signal_Kstarff9_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff9_m = new TH1D("Signal_Kstarff9_m", "Signal_Kstarff9_m", RarityBins, BinMIN, BinMAX);

    // faction
    TH1D* Signal_Kfrac_p = new TH1D("Signal_Kfrac_p", "Signal_Kfrac_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kfrac_m = new TH1D("Signal_Kfrac_m", "Signal_Kfrac_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarfrac_p = new TH1D("Signal_Kstarfrac_p", "Signal_Kstarfrac_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarfrac_m = new TH1D("Signal_Kstarfrac_m", "Signal_Kstarfrac_m", RarityBins, BinMIN, BinMAX);

    // fragmentation
    std::vector<std::string> name_Signal_Fragmentation;
    for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
        for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
            name_Signal_Fragmentation.push_back(std::string("Signal_Xs_frag_decay_") + std::to_string(MxsBin) + std::string("_") + std::to_string(Category) + std::string("_p"));
        }
    }
    for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
        for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
            name_Signal_Fragmentation.push_back(std::string("Signal_Xs_frag_decay_") + std::to_string(MxsBin) + std::string("_") + std::to_string(Category) + std::string("_m"));
        }
    }
    std::vector<TH1D*> Signal_Fragmentaions;
    for (unsigned int i = 0; i < name_Signal_Fragmentation.size(); i++) {
        TH1D* temp = new TH1D(name_Signal_Fragmentation.at(i).c_str(), name_Signal_Fragmentation.at(i).c_str(), RarityBins, BinMIN, BinMAX);
        Signal_Fragmentaions.push_back(temp);
    }

    // pf
    TH1D* Signal_pf_p = new TH1D("Signal_pf_p", "Signal_pf_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_pf_m = new TH1D("Signal_pf_m", "Signal_pf_m", RarityBins, BinMIN, BinMAX);

    // mb
    TH1D* Signal_mb_p = new TH1D("Signal_mb_p", "Signal_mb_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_mb_m = new TH1D("Signal_mb_m", "Signal_mb_m", RarityBins, BinMIN, BinMAX);

    // transition
    TH1D* Signal_transition_p = new TH1D("Signal_transition_p", "Signal_transition_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_transition_m = new TH1D("Signal_transition_m", "Signal_transition_m", RarityBins, BinMIN, BinMAX);

    // mK*
    TH1D* Signal_mKstar_p = new TH1D("Signal_mKstar_p", "Signal_mKstar_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_mKstar_m = new TH1D("Signal_mKstar_m", "Signal_mKstar_m", RarityBins, BinMIN, BinMAX);

    // BDTc
    //TH1D* CHG_BDTc_p = new TH1D("CHG_BDTc_p", "CHG_BDTc_p", RarityBins, BinMIN, BinMAX);
    //TH1D* MIX_BDTc_p = new TH1D("MIX_BDTc_p", "MIX_BDTc_p", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_BDTc_p = new TH1D("UUBAR_BDTc_p", "UUBAR_BDTc_p", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_BDTc_p = new TH1D("DDBAR_BDTc_p", "DDBAR_BDTc_p", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_BDTc_p = new TH1D("SSBAR_BDTc_p", "SSBAR_BDTc_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_BDTc_p = new TH1D("CHARM_BDTc_p", "CHARM_BDTc_p", RarityBins, BinMIN, BinMAX);
    //TH1D* CHG_BDTc_m = new TH1D("CHG_BDTc_m", "CHG_BDTc_m", RarityBins, BinMIN, BinMAX);
    //TH1D* MIX_BDTc_m = new TH1D("MIX_BDTc_m", "MIX_BDTc_m", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_BDTc_m = new TH1D("UUBAR_BDTc_m", "UUBAR_BDTc_m", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_BDTc_m = new TH1D("DDBAR_BDTc_m", "DDBAR_BDTc_m", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_BDTc_m = new TH1D("SSBAR_BDTc_m", "SSBAR_BDTc_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_BDTc_m = new TH1D("CHARM_BDTc_m", "CHARM_BDTc_m", RarityBins, BinMIN, BinMAX);

    // B2Knn BR
    TH1D* CHG_Knn_p = new TH1D("CHG_Knn_p", "CHG_Knn_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Knn_p = new TH1D("MIX_Knn_p", "MIX_Knn_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Knn_m = new TH1D("CHG_Knn_m", "CHG_Knn_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Knn_m = new TH1D("MIX_Knn_m", "MIX_Knn_m", RarityBins, BinMIN, BinMAX);

    TH1D* CHG_Kstarnn_p = new TH1D("CHG_Kstarnn_p", "CHG_Kstarnn_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kstarnn_p = new TH1D("MIX_Kstarnn_p", "MIX_Kstarnn_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kstarnn_m = new TH1D("CHG_Kstarnn_m", "CHG_Kstarnn_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kstarnn_m = new TH1D("MIX_Kstarnn_m", "MIX_Kstarnn_m", RarityBins, BinMIN, BinMAX);

    TH1D* CHG_K0nn_p = new TH1D("CHG_K0nn_p", "CHG_K0nn_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_K0nn_p = new TH1D("MIX_K0nn_p", "MIX_K0nn_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_K0nn_m = new TH1D("CHG_K0nn_m", "CHG_K0nn_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_K0nn_m = new TH1D("MIX_K0nn_m", "MIX_K0nn_m", RarityBins, BinMIN, BinMAX);

    TH1D* CHG_K0starnn_p = new TH1D("CHG_K0starnn_p", "CHG_K0starnn_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_K0starnn_p = new TH1D("MIX_K0starnn_p", "MIX_K0starnn_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_K0starnn_m = new TH1D("CHG_K0starnn_m", "CHG_K0starnn_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_K0starnn_m = new TH1D("MIX_K0starnn_m", "MIX_K0starnn_m", RarityBins, BinMIN, BinMAX);

    // multiplicity
    TH1D** Signal_multiplicity_p;
    TH1D** CHG_multiplicity_p;
    TH1D** MIX_multiplicity_p;
    TH1D** UUBAR_multiplicity_p;
    TH1D** DDBAR_multiplicity_p;
    TH1D** SSBAR_multiplicity_p;
    TH1D** CHARM_multiplicity_p;
    TH1D** Signal_multiplicity_m;
    TH1D** CHG_multiplicity_m;
    TH1D** MIX_multiplicity_m;
    TH1D** UUBAR_multiplicity_m;
    TH1D** DDBAR_multiplicity_m;
    TH1D** SSBAR_multiplicity_m;
    TH1D** CHARM_multiplicity_m;

    // B -> [D -> X KL0] anything
    TH1D* Signal_BtoDtoXKL_p = new TH1D("Signal_BtoDtoXKL_p", "Signal_BtoDtoXKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_BtoDtoXKL_p = new TH1D("CHG_BtoDtoXKL_p", "CHG_BtoDtoXKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_BtoDtoXKL_p = new TH1D("MIX_BtoDtoXKL_p", "MIX_BtoDtoXKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_BtoDtoXKL_m = new TH1D("Signal_BtoDtoXKL_m", "Signal_BtoDtoXKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_BtoDtoXKL_m = new TH1D("CHG_BtoDtoXKL_m", "CHG_BtoDtoXKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_BtoDtoXKL_m = new TH1D("MIX_BtoDtoXKL_m", "MIX_BtoDtoXKL_m", RarityBins, BinMIN, BinMAX);

    // all of uncorrelated uncertainties
    TH1D* Signal_all_uncorrelated = new TH1D("Signal_all_uncorrelated", "Signal_all_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_all_uncorrelated = new TH1D("CHG_all_uncorrelated", "CHG_all_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_all_uncorrelated = new TH1D("MIX_all_uncorrelated", "MIX_all_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_all_uncorrelated = new TH1D("UUBAR_all_uncorrelated", "UUBAR_all_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_all_uncorrelated = new TH1D("DDBAR_all_uncorrelated", "DDBAR_all_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_all_uncorrelated = new TH1D("SSBAR_all_uncorrelated", "SSBAR_all_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_all_uncorrelated = new TH1D("CHARM_all_uncorrelated", "CHARM_all_uncorrelated", RarityBins, BinMIN, BinMAX);
    /* ====================================== */



    /* ====================================== */
    // define path for Ntuple
    const char* MC_dirname_SIGNAL = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut";

    const char* MC_dirname_CHG = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori/CHG_analysis/test_v000/final_output_root_after_MVA_Application_after_cut";
    const char* MC_dirname_MIX = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori/MIX_analysis/test_v000/final_output_root_after_MVA_Application_after_cut";
    const char* MC_dirname_UUBAR = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori/UUBAR_analysis/test_v000/final_output_root_after_MVA_Application_after_cut";
    const char* MC_dirname_DDBAR = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori/DDBAR_analysis/test_v000/final_output_root_after_MVA_Application_after_cut";
    const char* MC_dirname_SSBAR = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori/SSBAR_analysis/test_v000/final_output_root_after_MVA_Application_after_cut";
    const char* MC_dirname_CHARM = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori/CHARM_analysis/test_v000/final_output_root_after_MVA_Application_after_cut";

    // for signal modeling
    const char* MC_dirname_syst = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_syst/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut";
 
    // for FEI
    const char* FEI_correlated_info = "./FEI_selected.txt";
    const char* FEI_uncorrelated_info = "./FEI_cov_remain_truncated.txt";

    // for PID
    const char* KID_correlated_info = "./KID_selected.txt";
    const char* PID_correlated_info = "./PID_selected.txt";
    const char* KID_uncorrelated_info = "./KID_cov_remain_truncated.txt";
    const char* PID_uncorrelated_info = "./PID_cov_remain_truncated.txt";

    // for BB BR
    const char* BR_correlated_info = "./BR_selected.txt";
    const char* BR_uncorrelated_info = "./BR_cov_remain_truncated.txt";

    // for pi0
    const char* pi0_correlated_info = "./pi0_selected.txt";
    const char* pi0_uncorrelated_info = "./pi0_cov_remain_truncated.txt";

    /* ====================================== */



    /* ====================================== */
    // Get PDFs
    // get nominal pdfs
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_nominal, "Bplus", "SIGNAL", Scale_Kplus_test, "B2Knunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_nominal, "Bplus", "SIGNAL", Scale_Kplusstar_test, "otherwise");
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_nominal, "Bplus", "SIGNAL", Scale_Xsu_nonresonant_test, "B2Xsnunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_nominal, "Bzero", "SIGNAL", Scale_K0_test, "B02K0nunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_nominal, "Bzero", "SIGNAL", Scale_K0star_test, "otherwise");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_nominal, "Bzero", "SIGNAL", Scale_Xsd_nonresonant_test, "B02Xsnunu");

    double NevtCHGWithoutBDTc = GetNominalPDFs(MC_dirname_CHG, "root", CHG_nominal, "Bplus", "CHG", Scale_CHG_test, "otherwise");
    double NevtMIXWithoutBDTc = GetNominalPDFs(MC_dirname_MIX, "root", MIX_nominal, "Bzero", "MIX", Scale_MIX_test, "otherwise");
    double NevtUUBARWithoutBDTc = GetNominalPDFs(MC_dirname_UUBAR, "root", UUBAR_nominal, "Continuum", "UUBAR", Scale_UUBAR_test, "otherwise");
    double NevtDDBARWithoutBDTc = GetNominalPDFs(MC_dirname_DDBAR, "root", DDBAR_nominal, "Continuum", "DDBAR", Scale_DDBAR_test, "otherwise");
    double NevtSSBARWithoutBDTc = GetNominalPDFs(MC_dirname_SSBAR, "root", SSBAR_nominal, "Continuum", "SSBAR", Scale_SSBAR_test, "otherwise");
    double NevtCHARMWithoutBDTc = GetNominalPDFs(MC_dirname_CHARM, "root", CHARM_nominal, "Continuum", "CHARM", Scale_CHARM_test, "otherwise");

    // get track uncertainty pdfs
    GetTrackPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_track_p, "Bplus", "SIGNAL", true, Scale_Kplus_test, "B2Knunu");
    GetTrackPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_track_p, "Bplus", "SIGNAL", true, Scale_Kplusstar_test, "otherwise");
    GetTrackPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_track_p, "Bplus", "SIGNAL", true, Scale_Xsu_nonresonant_test, "B2Xsnunu");
    GetTrackPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_track_p, "Bzero", "SIGNAL", true, Scale_K0_test, "B02K0nunu");
    GetTrackPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_track_p, "Bzero", "SIGNAL", true, Scale_K0star_test, "otherwise");
    GetTrackPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_track_p, "Bzero", "SIGNAL", true, Scale_Xsd_nonresonant_test, "B02Xsnunu");

    GetTrackPDFs(MC_dirname_CHG, "root", CHG_track_p, "Bplus", "CHG", true, Scale_CHG_test, "otherwise");
    GetTrackPDFs(MC_dirname_MIX, "root", MIX_track_p, "Bzero", "MIX", true, Scale_MIX_test, "otherwise");
    GetTrackPDFs(MC_dirname_UUBAR, "root", UUBAR_track_p, "Continuum", "UUBAR", true, Scale_UUBAR_test, "otherwise");
    GetTrackPDFs(MC_dirname_DDBAR, "root", DDBAR_track_p, "Continuum", "DDBAR", true, Scale_DDBAR_test, "otherwise");
    GetTrackPDFs(MC_dirname_SSBAR, "root", SSBAR_track_p, "Continuum", "SSBAR", true, Scale_SSBAR_test, "otherwise");
    GetTrackPDFs(MC_dirname_CHARM, "root", CHARM_track_p, "Continuum", "CHARM", true, Scale_CHARM_test, "otherwise");

    GetTrackPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_track_m, "Bplus", "SIGNAL", false, Scale_Kplus_test, "B2Knunu");
    GetTrackPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_track_m, "Bplus", "SIGNAL", false, Scale_Kplusstar_test, "otherwise");
    GetTrackPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_track_m, "Bplus", "SIGNAL", false, Scale_Xsu_nonresonant_test, "B2Xsnunu");
    GetTrackPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_track_m, "Bzero", "SIGNAL", false, Scale_K0_test, "B02K0nunu");
    GetTrackPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_track_m, "Bzero", "SIGNAL", false, Scale_K0star_test, "otherwise");
    GetTrackPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_track_m, "Bzero", "SIGNAL", false, Scale_Xsd_nonresonant_test, "B02Xsnunu");

    GetTrackPDFs(MC_dirname_CHG, "root", CHG_track_m, "Bplus", "CHG", false, Scale_CHG_test, "otherwise");
    GetTrackPDFs(MC_dirname_MIX, "root", MIX_track_m, "Bzero", "MIX", false, Scale_MIX_test, "otherwise");
    GetTrackPDFs(MC_dirname_UUBAR, "root", UUBAR_track_m, "Continuum", "UUBAR", false, Scale_UUBAR_test, "otherwise");
    GetTrackPDFs(MC_dirname_DDBAR, "root", DDBAR_track_m, "Continuum", "DDBAR", false, Scale_DDBAR_test, "otherwise");
    GetTrackPDFs(MC_dirname_SSBAR, "root", SSBAR_track_m, "Continuum", "SSBAR", false, Scale_SSBAR_test, "otherwise");
    GetTrackPDFs(MC_dirname_CHARM, "root", CHARM_track_m, "Continuum", "CHARM", false, Scale_CHARM_test, "otherwise");

    // get KS0 uncertainty pdfs
    GetKS0PDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_KS0_p, "Bplus", "SIGNAL", true, Scale_Kplus_test, "B2Knunu");
    GetKS0PDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_KS0_p, "Bplus", "SIGNAL", true, Scale_Kplusstar_test, "otherwise");
    GetKS0PDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_KS0_p, "Bplus", "SIGNAL", true, Scale_Xsu_nonresonant_test, "B2Xsnunu");
    GetKS0PDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_KS0_p, "Bzero", "SIGNAL", true, Scale_K0_test, "B02K0nunu");
    GetKS0PDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_KS0_p, "Bzero", "SIGNAL", true, Scale_K0star_test, "otherwise");
    GetKS0PDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_KS0_p, "Bzero", "SIGNAL", true, Scale_Xsd_nonresonant_test, "B02Xsnunu");

    GetKS0PDFs(MC_dirname_CHG, "root", CHG_KS0_p, "Bplus", "CHG", true, Scale_CHG_test, "otherwise");
    GetKS0PDFs(MC_dirname_MIX, "root", MIX_KS0_p, "Bzero", "MIX", true, Scale_MIX_test, "otherwise");
    GetKS0PDFs(MC_dirname_UUBAR, "root", UUBAR_KS0_p, "Continuum", "UUBAR", true, Scale_UUBAR_test, "otherwise");
    GetKS0PDFs(MC_dirname_DDBAR, "root", DDBAR_KS0_p, "Continuum", "DDBAR", true, Scale_DDBAR_test, "otherwise");
    GetKS0PDFs(MC_dirname_SSBAR, "root", SSBAR_KS0_p, "Continuum", "SSBAR", true, Scale_SSBAR_test, "otherwise");
    GetKS0PDFs(MC_dirname_CHARM, "root", CHARM_KS0_p, "Continuum", "CHARM", true, Scale_CHARM_test, "otherwise");

    GetKS0PDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_KS0_m, "Bplus", "SIGNAL", false, Scale_Kplus_test, "B2Knunu");
    GetKS0PDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_KS0_m, "Bplus", "SIGNAL", false, Scale_Kplusstar_test, "otherwise");
    GetKS0PDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_KS0_m, "Bplus", "SIGNAL", false, Scale_Xsu_nonresonant_test, "B2Xsnunu");
    GetKS0PDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_KS0_m, "Bzero", "SIGNAL", false, Scale_K0_test, "B02K0nunu");
    GetKS0PDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_KS0_m, "Bzero", "SIGNAL", false, Scale_K0star_test, "otherwise");
    GetKS0PDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_KS0_m, "Bzero", "SIGNAL", false, Scale_Xsd_nonresonant_test, "B02Xsnunu");

    GetKS0PDFs(MC_dirname_CHG, "root", CHG_KS0_m, "Bplus", "CHG", false, Scale_CHG_test, "otherwise");
    GetKS0PDFs(MC_dirname_MIX, "root", MIX_KS0_m, "Bzero", "MIX", false, Scale_MIX_test, "otherwise");
    GetKS0PDFs(MC_dirname_UUBAR, "root", UUBAR_KS0_m, "Continuum", "UUBAR", false, Scale_UUBAR_test, "otherwise");
    GetKS0PDFs(MC_dirname_DDBAR, "root", DDBAR_KS0_m, "Continuum", "DDBAR", false, Scale_DDBAR_test, "otherwise");
    GetKS0PDFs(MC_dirname_SSBAR, "root", SSBAR_KS0_m, "Continuum", "SSBAR", false, Scale_SSBAR_test, "otherwise");
    GetKS0PDFs(MC_dirname_CHARM, "root", CHARM_KS0_m, "Continuum", "CHARM", false, Scale_CHARM_test, "otherwise");

    // get FEI uncertainty pdfs (correlated)
    int NPDFs_FEI = GetFEIcorrelatedPDFs(FEI_correlated_info, CHG_nominal, MIX_nominal, Signal_nominal, &CHG_FEI_correlated, &MIX_FEI_correlated, &Signal_FEI_correlated);

    // get FEI uncertainty pdfs (uncorrelated)
    GetFEIUncorrelatedPDFs(FEI_uncorrelated_info, CHG_FEI_uncorrelated, MIX_FEI_uncorrelated, Signal_FEI_uncorrelated);

    // get KID uncertainty pdfs (correlated)
    int NPDFs_KID = GetPIDcorrelatedPDFs(KID_correlated_info, CHG_nominal, MIX_nominal, UUBAR_nominal, DDBAR_nominal, SSBAR_nominal, CHARM_nominal, Signal_nominal, &CHG_KID_correlated, &MIX_KID_correlated, &UUBAR_KID_correlated, &DDBAR_KID_correlated, &SSBAR_KID_correlated, &CHARM_KID_correlated, &Signal_KID_correlated, true);
    int NPDFs_PID = GetPIDcorrelatedPDFs(PID_correlated_info, CHG_nominal, MIX_nominal, UUBAR_nominal, DDBAR_nominal, SSBAR_nominal, CHARM_nominal, Signal_nominal, &CHG_PID_correlated, &MIX_PID_correlated, &UUBAR_PID_correlated, &DDBAR_PID_correlated, &SSBAR_PID_correlated, &CHARM_PID_correlated, &Signal_PID_correlated, false);

    // get KID uncertainty pdfs (uncorrelated)
    GetPIDUncorrelatedPDFs(KID_uncorrelated_info, CHG_KID_uncorrelated, MIX_KID_uncorrelated, UUBAR_KID_uncorrelated, DDBAR_KID_uncorrelated, SSBAR_KID_uncorrelated, CHARM_KID_uncorrelated, Signal_KID_uncorrelated);
    GetPIDUncorrelatedPDFs(PID_uncorrelated_info, CHG_PID_uncorrelated, MIX_PID_uncorrelated, UUBAR_PID_uncorrelated, DDBAR_PID_uncorrelated, SSBAR_PID_uncorrelated, CHARM_PID_uncorrelated, Signal_PID_uncorrelated);

    // get BB BR uncertainty pdfs (correlated)
    int NPDFs_BR = GetBRcorrelatedPDFs(BR_correlated_info, CHG_nominal, MIX_nominal, Signal_nominal, &CHG_BR_correlated, &MIX_BR_correlated, &Signal_BR_correlated);

    // get BB BR uncertainty pdfs (uncorrelated)
    GetBRUncorrelatedPDFs(BR_uncorrelated_info, CHG_BR_uncorrelated, MIX_BR_uncorrelated, Signal_BR_uncorrelated);

    // get pi0 uncertainty pdfs (correlated)
    int NPDFs_pi0 = Getpi0correlatedPDFs(pi0_correlated_info, CHG_nominal, MIX_nominal, UUBAR_nominal, DDBAR_nominal, SSBAR_nominal, CHARM_nominal, Signal_nominal, &CHG_pi0_correlated, &MIX_pi0_correlated, &UUBAR_pi0_correlated, &DDBAR_pi0_correlated, &SSBAR_pi0_correlated, &CHARM_pi0_correlated, &Signal_pi0_correlated);

    // get pi0 uncertainty pdfs (uncorrelated)
    Getpi0UncorrelatedPDFs(pi0_uncorrelated_info, CHG_pi0_uncorrelated, MIX_pi0_uncorrelated, UUBAR_pi0_uncorrelated, DDBAR_pi0_uncorrelated, SSBAR_pi0_uncorrelated, CHARM_pi0_uncorrelated, Signal_pi0_uncorrelated);

    // get Kff uncertainty pdfs
    TH1D* Signal_Kff_temp = new TH1D("Signal_Kff_temp", "Signal_Kff_temp", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kff_array[7] = { Signal_Kff_temp, Signal_Kff1_p, Signal_Kff1_m, Signal_Kff2_p, Signal_Kff2_m, Signal_Kff3_p, Signal_Kff3_m };
    double Correction_factor_BR_Kplus[7] = { 0.0 };
    double Correction_factor_BR_Kzero[7] = { 0.0 };

    GetKffPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_Kff_array, Correction_factor_BR_Kplus, "Bplus", 1, Scale_Kplus_test, "B2Knunu");
    for (int i = 0; i < 7; i++) GetNominalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_Kff_array[i], "Bplus", "SIGNAL", Scale_Kplusstar_test, "otherwise");
    for (int i = 0; i < 7; i++) {
        const double non_res_Xsu_correction_factor = (BR_Xsu_nonresonant_nunubar + BR_Kplus_nunubar - BR_Kplus_nunubar * Correction_factor_BR_Kplus[i]) / BR_Xsu_nonresonant_nunubar;
        GetNominalPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_Kff_array[i], "Bplus", "SIGNAL", Scale_Xsu_nonresonant_test * non_res_Xsu_correction_factor, "B2Xsnunu");
    }
    GetKffPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_Kff_array, Correction_factor_BR_Kzero, "Bzero", 0, Scale_K0_test, "B02K0nunu");
    for (int i = 0; i < 7; i++) GetNominalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_Kff_array[i], "Bzero", "SIGNAL", Scale_K0star_test, "otherwise");
    for (int i = 0; i < 7; i++) {
        const double non_res_Xsd_correction_factor = (BR_Xsd_nonresonant_nunubar + BR_K0_nunubar - BR_K0_nunubar * Correction_factor_BR_Kzero[i]) / BR_Xsd_nonresonant_nunubar;
        GetNominalPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_Kff_array[i], "Bzero", "SIGNAL", Scale_Xsd_nonresonant_test * non_res_Xsd_correction_factor, "B02Xsnunu");
    }

    // get Kstar ff uncertainty pdfs
    TH1D* Signal_Kstarff_temp = new TH1D("Signal_Kstarff_temp", "Signal_Kstarff_temp", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff_array[19] = { Signal_Kstarff_temp,
        Signal_Kstarff1_p, Signal_Kstarff1_m,
        Signal_Kstarff2_p, Signal_Kstarff2_m,
        Signal_Kstarff3_p, Signal_Kstarff3_m,
        Signal_Kstarff4_p, Signal_Kstarff4_m,
        Signal_Kstarff5_p, Signal_Kstarff5_m,
        Signal_Kstarff6_p, Signal_Kstarff6_m,
        Signal_Kstarff7_p, Signal_Kstarff7_m,
        Signal_Kstarff8_p, Signal_Kstarff8_m,
        Signal_Kstarff9_p, Signal_Kstarff9_m
    };
    double Correction_factor_BR_Kstarplus[19] = { 0.0 };
    double Correction_factor_BR_Kstarzero[19] = { 0.0 };

    for (int i = 0; i < 19; i++) GetNominalPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_Kstarff_array[i], "Bplus", "SIGNAL", Scale_Kplus_test, "B2Knunu");
    GetKstarffPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_Kstarff_array, Correction_factor_BR_Kstarplus, "Bplus", 1, Scale_Kplusstar_test);
    for (int i = 0; i < 19; i++) {
        const double non_res_Xsu_correction_factor = (BR_Xsu_nonresonant_nunubar + BR_Kplusstar_nunubar - BR_Kplusstar_nunubar * Correction_factor_BR_Kstarplus[i]) / BR_Xsu_nonresonant_nunubar;
        GetNominalPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_Kstarff_array[i], "Bplus", "SIGNAL", Scale_Xsu_nonresonant_test * non_res_Xsu_correction_factor, "B2Xsnunu");
    }
    for (int i = 0; i < 19; i++) GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_Kstarff_array[i], "Bzero", "SIGNAL", Scale_K0_test, "B02K0nunu");
    GetKstarffPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_Kstarff_array, Correction_factor_BR_Kstarzero, "Bzero", 0, Scale_K0star_test);
    for (int i = 0; i < 19; i++) {
        const double non_res_Xsd_correction_factor = (BR_Xsd_nonresonant_nunubar + BR_K0star_nunubar - BR_K0star_nunubar * Correction_factor_BR_Kstarzero[i]) / BR_Xsd_nonresonant_nunubar;
        GetNominalPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_Kstarff_array[i], "Bzero", "SIGNAL", Scale_Xsd_nonresonant_test * non_res_Xsd_correction_factor, "B02Xsnunu");
    }

    GetNominalPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_Kfrac_p, "Bplus", "SIGNAL", Scale_Kplus_test * ((BR_Kplus_nunubar + Sigma_BR_Kplus_nunubar) / BR_Kplus_nunubar), "B2Knunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_Kfrac_p, "Bplus", "SIGNAL", Scale_Kplusstar_test, "otherwise");
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_Kfrac_p, "Bplus", "SIGNAL", Scale_Xsu_nonresonant_test * ((BR_Xsu_nonresonant_nunubar - Sigma_BR_Kplus_nunubar) / BR_Xsu_nonresonant_nunubar), "B2Xsnunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_Kfrac_p, "Bzero", "SIGNAL", Scale_K0_test * ((BR_K0_nunubar + Sigma_BR_K0_nunubar) / BR_K0_nunubar), "B02K0nunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_Kfrac_p, "Bzero", "SIGNAL", Scale_K0star_test, "otherwise");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_Kfrac_p, "Bzero", "SIGNAL", Scale_Xsd_nonresonant_test * ((BR_Xsd_nonresonant_nunubar - Sigma_BR_K0_nunubar) / BR_Xsd_nonresonant_nunubar), "B02Xsnunu");

    GetNominalPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_Kfrac_m, "Bplus", "SIGNAL", Scale_Kplus_test * ((BR_Kplus_nunubar - Sigma_BR_Kplus_nunubar) / BR_Kplus_nunubar), "B2Knunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_Kfrac_m, "Bplus", "SIGNAL", Scale_Kplusstar_test, "otherwise");
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_Kfrac_m, "Bplus", "SIGNAL", Scale_Xsu_nonresonant_test * ((BR_Xsu_nonresonant_nunubar + Sigma_BR_Kplus_nunubar) / BR_Xsu_nonresonant_nunubar), "B2Xsnunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_Kfrac_m, "Bzero", "SIGNAL", Scale_K0_test * ((BR_K0_nunubar - Sigma_BR_K0_nunubar) / BR_K0_nunubar), "B02K0nunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_Kfrac_m, "Bzero", "SIGNAL", Scale_K0star_test, "otherwise");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_Kfrac_m, "Bzero", "SIGNAL", Scale_Xsd_nonresonant_test * ((BR_Xsd_nonresonant_nunubar + Sigma_BR_K0_nunubar) / BR_Xsd_nonresonant_nunubar), "B02Xsnunu");

    GetNominalPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_Kstarfrac_p, "Bplus", "SIGNAL", Scale_Kplus_test, "B2Knunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_Kstarfrac_p, "Bplus", "SIGNAL", Scale_Kplusstar_test * ((BR_Kplusstar_nunubar + Sigma_BR_Kplusstar_nunubar) / BR_Kplusstar_nunubar), "otherwise");
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_Kstarfrac_p, "Bplus", "SIGNAL", Scale_Xsu_nonresonant_test * ((BR_Xsu_nonresonant_nunubar - Sigma_BR_Kplusstar_nunubar) / BR_Xsu_nonresonant_nunubar), "B2Xsnunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_Kstarfrac_p, "Bzero", "SIGNAL", Scale_K0_test, "B02K0nunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_Kstarfrac_p, "Bzero", "SIGNAL", Scale_K0star_test * ((BR_K0star_nunubar + Sigma_BR_K0star_nunubar) / BR_K0star_nunubar), "otherwise");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_Kstarfrac_p, "Bzero", "SIGNAL", Scale_Xsd_nonresonant_test * ((BR_Xsd_nonresonant_nunubar - Sigma_BR_K0star_nunubar) / BR_Xsd_nonresonant_nunubar), "B02Xsnunu");

    GetNominalPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_Kstarfrac_m, "Bplus", "SIGNAL", Scale_Kplus_test, "B2Knunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_Kstarfrac_m, "Bplus", "SIGNAL", Scale_Kplusstar_test * ((BR_Kplusstar_nunubar - Sigma_BR_Kplusstar_nunubar) / BR_Kplusstar_nunubar), "otherwise");
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_Kstarfrac_m, "Bplus", "SIGNAL", Scale_Xsu_nonresonant_test * ((BR_Xsu_nonresonant_nunubar + Sigma_BR_Kplusstar_nunubar) / BR_Xsu_nonresonant_nunubar), "B2Xsnunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_Kstarfrac_m, "Bzero", "SIGNAL", Scale_K0_test, "B02K0nunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_Kstarfrac_m, "Bzero", "SIGNAL", Scale_K0star_test * ((BR_K0star_nunubar - Sigma_BR_K0star_nunubar) / BR_K0star_nunubar), "otherwise");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_Kstarfrac_m, "Bzero", "SIGNAL", Scale_Xsd_nonresonant_test * ((BR_Xsd_nonresonant_nunubar + Sigma_BR_K0star_nunubar) / BR_Xsd_nonresonant_nunubar), "B02Xsnunu");

    // get fragmentation uncertainty pdfs
    for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
        for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
            int temp_index = MxsBin * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma) + Category;
            bool IsItUp = true;
            if (name_Signal_Fragmentation.at(temp_index).find("_p") != string::npos) IsItUp = true;
            else if (name_Signal_Fragmentation.at(temp_index).find("_m") != string::npos) IsItUp = false;
            else {
                printf("[ERROR] unexpected error!\n");
                exit(1);
            }

            GetNominalPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_Fragmentaions.at(temp_index), "Bplus", "SIGNAL", Scale_Kplus_test, "B2Knunu");
            GetNominalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_Fragmentaions.at(temp_index), "Bplus", "SIGNAL", Scale_Kplusstar_test, "otherwise");
            GetFragmentationPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_Fragmentaions.at(temp_index), "Bplus", "SIGNAL", MxsBin, Category, IsItUp, Scale_Xsu_nonresonant_test, "B2Xsnunu");
            GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_Fragmentaions.at(temp_index), "Bzero", "SIGNAL", Scale_K0_test, "B02K0nunu");
            GetNominalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_Fragmentaions.at(temp_index), "Bzero", "SIGNAL", Scale_K0star_test, "otherwise");
            GetFragmentationPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_Fragmentaions.at(temp_index), "Bzero", "SIGNAL", MxsBin, Category, IsItUp, Scale_Xsd_nonresonant_test, "B02Xsnunu");

        }
    }


    // get pf uncertainty pdfs
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_pf_p, "Bplus", "SIGNAL", Scale_Kplus_test, "B2Knunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_pf_p, "Bplus", "SIGNAL", Scale_Kplusstar_test, "otherwise");
    GetpfPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_pf_p, "Bplus", "SIGNAL", true, true, Scale_Xsu_nonresonant_test, "B2Xsnunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_pf_p, "Bzero", "SIGNAL", Scale_K0_test, "B02K0nunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_pf_p, "Bzero", "SIGNAL", Scale_K0star_test, "otherwise");
    GetpfPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_pf_p, "Bzero", "SIGNAL", true, false, Scale_Xsd_nonresonant_test, "B02Xsnunu");

    GetNominalPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_pf_m, "Bplus", "SIGNAL", Scale_Kplus_test, "B2Knunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_pf_m, "Bplus", "SIGNAL", Scale_Kplusstar_test, "otherwise");
    GetpfPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_pf_m, "Bplus", "SIGNAL", false, true, Scale_Xsu_nonresonant_test, "B2Xsnunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_pf_m, "Bzero", "SIGNAL", Scale_K0_test, "B02K0nunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_pf_m, "Bzero", "SIGNAL", Scale_K0star_test, "otherwise");
    GetpfPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_pf_m, "Bzero", "SIGNAL", false, false, Scale_Xsd_nonresonant_test, "B02Xsnunu");

    // mb
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_mb_p, "Bplus", "SIGNAL", Scale_Kplus_test, "B2Knunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_mb_p, "Bplus", "SIGNAL", Scale_Kplusstar_test, "otherwise");
    GetmbPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_mb_p, "Bplus", "SIGNAL", true, true, Scale_Xsu_nonresonant_test, "B2Xsnunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_mb_p, "Bzero", "SIGNAL", Scale_K0_test, "B02K0nunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_mb_p, "Bzero", "SIGNAL", Scale_K0star_test, "otherwise");
    GetmbPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_mb_p, "Bzero", "SIGNAL", true, false, Scale_Xsd_nonresonant_test, "B02Xsnunu");

    GetNominalPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_mb_m, "Bplus", "SIGNAL", Scale_Kplus_test, "B2Knunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_mb_m, "Bplus", "SIGNAL", Scale_Kplusstar_test, "otherwise");
    GetmbPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_mb_m, "Bplus", "SIGNAL", false, true, Scale_Xsu_nonresonant_test, "B2Xsnunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_mb_m, "Bzero", "SIGNAL", Scale_K0_test, "B02K0nunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_mb_m, "Bzero", "SIGNAL", Scale_K0star_test, "otherwise");
    GetmbPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_mb_m, "Bzero", "SIGNAL", false, false, Scale_Xsd_nonresonant_test, "B02Xsnunu");

    // transition
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_transition_p, "Bplus", "SIGNAL", Scale_Kplus_test, "B2Knunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_transition_p, "Bplus", "SIGNAL", Scale_Kplusstar_test, "otherwise");
    GetTransitionPDFs(MC_dirname_syst, "B2Xsnunu_Htransition", Signal_transition_p, "Bplus", "SIGNAL", true, true, Scale_Xsu_nonresonant_syst, "B2Xsnunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_transition_p, "Bzero", "SIGNAL", Scale_K0_test, "B02K0nunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_transition_p, "Bzero", "SIGNAL", Scale_K0star_test, "otherwise");
    GetTransitionPDFs(MC_dirname_syst, "B02Xsnunu_Htransition", Signal_transition_p, "Bzero", "SIGNAL", true, false, Scale_Xsd_nonresonant_syst, "B02Xsnunu");

    GetNominalPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_transition_m, "Bplus", "SIGNAL", Scale_Kplus_test, "B2Knunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_transition_m, "Bplus", "SIGNAL", Scale_Kplusstar_test, "otherwise");
    GetTransitionPDFs(MC_dirname_syst, "B2Xsnunu_Ltransition", Signal_transition_m, "Bplus", "SIGNAL", false, true, Scale_Xsu_nonresonant_syst, "B2Xsnunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_transition_m, "Bzero", "SIGNAL", Scale_K0_test, "B02K0nunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_transition_m, "Bzero", "SIGNAL", Scale_K0star_test, "otherwise");
    GetTransitionPDFs(MC_dirname_syst, "B02Xsnunu_Ltransition", Signal_transition_m, "Bzero", "SIGNAL", false, false, Scale_Xsd_nonresonant_syst, "B02Xsnunu");

    // mKstar fixed
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_mKstar_p, "Bplus", "SIGNAL", Scale_Kplus_test, "B2Knunu");
    GetKstardeltaPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_mKstar_p, "Bplus", "SIGNAL", true, Scale_Kplusstar_test, "otherwise");
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_mKstar_p, "Bplus", "SIGNAL", Scale_Xsu_nonresonant_test, "B2Xsnunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_mKstar_p, "Bzero", "SIGNAL", Scale_K0_test, "B02K0nunu");
    GetKstardeltaPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_mKstar_p, "Bzero", "SIGNAL", false, Scale_K0star_test, "otherwise");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_mKstar_p, "Bzero", "SIGNAL", Scale_Xsd_nonresonant_test, "B02Xsnunu");

    GetNegativeChangePDFs(Signal_nominal, Signal_mKstar_p, Signal_mKstar_m);

    // get BDTc uncertainty pdfs
    //double NevtCHGWithBDTc = GetNevtWithBDTc(MC_dirname_CHG, "root", "Bplus", "CHG", Scale_CHG_test, "otherwise");
    //double NevtMIXWithBDTc = GetNevtWithBDTc(MC_dirname_MIX, "root", "Bzero", "MIX", Scale_MIX_test, "otherwise");
    double NevtUUBARWithBDTc = GetNevtWithBDTc(MC_dirname_UUBAR, "root", "Continuum", "UUBAR", Scale_UUBAR_test, "otherwise");
    double NevtDDBARWithBDTc = GetNevtWithBDTc(MC_dirname_DDBAR, "root", "Continuum", "DDBAR", Scale_DDBAR_test, "otherwise");
    double NevtSSBARWithBDTc = GetNevtWithBDTc(MC_dirname_SSBAR, "root", "Continuum", "SSBAR", Scale_SSBAR_test, "otherwise");
    double NevtCHARMWithBDTc = GetNevtWithBDTc(MC_dirname_CHARM, "root", "Continuum", "CHARM", Scale_CHARM_test, "otherwise");

    double NormFactor = (NevtUUBARWithoutBDTc + NevtDDBARWithoutBDTc + NevtSSBARWithoutBDTc + NevtCHARMWithoutBDTc) / (NevtUUBARWithBDTc + NevtDDBARWithBDTc + NevtSSBARWithBDTc + NevtCHARMWithBDTc);

    //GetBDTcPDFs(MC_dirname_CHG, "root", CHG_BDTc_p, "Bplus", "CHG", Scale_CHG_test, "otherwise");
    //GetBDTcPDFs(MC_dirname_MIX, "root", MIX_BDTc_p, "Bzero", "MIX", Scale_MIX_test, "otherwise");
    double NevtUUBARWithBDTcWithNorm = GetBDTcPDFs(MC_dirname_UUBAR, "root", UUBAR_BDTc_p, "Continuum", "UUBAR", Scale_UUBAR_test * NormFactor, "otherwise");
    double NevtDDBARWithBDTcWithNorm = GetBDTcPDFs(MC_dirname_DDBAR, "root", DDBAR_BDTc_p, "Continuum", "DDBAR", Scale_DDBAR_test * NormFactor, "otherwise");
    double NevtSSBARWithBDTcWithNorm = GetBDTcPDFs(MC_dirname_SSBAR, "root", SSBAR_BDTc_p, "Continuum", "SSBAR", Scale_SSBAR_test * NormFactor, "otherwise");
    double NevtCHARMWithBDTcWithNorm = GetBDTcPDFs(MC_dirname_CHARM, "root", CHARM_BDTc_p, "Continuum", "CHARM", Scale_CHARM_test * NormFactor, "otherwise");

    //GetNegativeChangePDFs(CHG_nominal, CHG_BDTc_p, CHG_BDTc_m);
    //GetNegativeChangePDFs(MIX_nominal, MIX_BDTc_p, MIX_BDTc_m);
    GetNegativeChangePDFs(UUBAR_nominal, UUBAR_BDTc_p, UUBAR_BDTc_m);
    GetNegativeChangePDFs(DDBAR_nominal, DDBAR_BDTc_p, DDBAR_BDTc_m);
    GetNegativeChangePDFs(SSBAR_nominal, SSBAR_BDTc_p, SSBAR_BDTc_m);
    GetNegativeChangePDFs(CHARM_nominal, CHARM_BDTc_p, CHARM_BDTc_m);

    if (std::abs((NevtUUBARWithoutBDTc + NevtDDBARWithoutBDTc + NevtSSBARWithoutBDTc + NevtCHARMWithoutBDTc) - (NevtUUBARWithBDTcWithNorm + NevtDDBARWithBDTcWithNorm + NevtSSBARWithBDTcWithNorm + NevtCHARMWithBDTcWithNorm)) > 0.1) { // final check
        printf("ERROR! before and after BDTc correction shows different number of event!\n");
        exit(1);
    }

    // B2Knn BR
    GetKnnBRPDFs(MC_dirname_CHG, "root", CHG_Knn_p, "Bplus", "CHG", 1, true, Scale_CHG_test, "otherwise");
    GetKnnBRPDFs(MC_dirname_MIX, "root", MIX_Knn_p, "Bzero", "MIX", 1, true, Scale_MIX_test, "otherwise");

    GetKnnBRPDFs(MC_dirname_CHG, "root", CHG_Knn_m, "Bplus", "CHG", 1, false, Scale_CHG_test, "otherwise");
    GetKnnBRPDFs(MC_dirname_MIX, "root", MIX_Knn_m, "Bzero", "MIX", 1, false, Scale_MIX_test, "otherwise");
    
    GetKnnBRPDFs(MC_dirname_CHG, "root", CHG_Kstarnn_p, "Bplus", "CHG", 2, true, Scale_CHG_test, "otherwise");
    GetKnnBRPDFs(MC_dirname_MIX, "root", MIX_Kstarnn_p, "Bzero", "MIX", 2, true, Scale_MIX_test, "otherwise");

    GetKnnBRPDFs(MC_dirname_CHG, "root", CHG_Kstarnn_m, "Bplus", "CHG", 2, false, Scale_CHG_test, "otherwise");
    GetKnnBRPDFs(MC_dirname_MIX, "root", MIX_Kstarnn_m, "Bzero", "MIX", 2, false, Scale_MIX_test, "otherwise");

    GetKnnBRPDFs(MC_dirname_CHG, "root", CHG_K0nn_p, "Bplus", "CHG", 3, true, Scale_CHG_test, "otherwise");
    GetKnnBRPDFs(MC_dirname_MIX, "root", MIX_K0nn_p, "Bzero", "MIX", 3, true, Scale_MIX_test, "otherwise");

    GetKnnBRPDFs(MC_dirname_CHG, "root", CHG_K0nn_m, "Bplus", "CHG", 3, false, Scale_CHG_test, "otherwise");
    GetKnnBRPDFs(MC_dirname_MIX, "root", MIX_K0nn_m, "Bzero", "MIX", 3, false, Scale_MIX_test, "otherwise");

    GetKnnBRPDFs(MC_dirname_CHG, "root", CHG_K0starnn_p, "Bplus", "CHG", 4, true, Scale_CHG_test, "otherwise");
    GetKnnBRPDFs(MC_dirname_MIX, "root", MIX_K0starnn_p, "Bzero", "MIX", 4, true, Scale_MIX_test, "otherwise");

    GetKnnBRPDFs(MC_dirname_CHG, "root", CHG_K0starnn_m, "Bplus", "CHG", 4, false, Scale_CHG_test, "otherwise");
    GetKnnBRPDFs(MC_dirname_MIX, "root", MIX_K0starnn_m, "Bzero", "MIX", 4, false, Scale_MIX_test, "otherwise");

    // multiplicity
    CheckMaxNgamma(MC_dirname_SIGNAL, "B2Knunu");
    CheckMaxNgamma(MC_dirname_SIGNAL, "B2Kstarnunu");
    CheckMaxNgamma(MC_dirname_SIGNAL, "B2Xsnunu");
    CheckMaxNgamma(MC_dirname_SIGNAL, "B02K0nunu");
    CheckMaxNgamma(MC_dirname_SIGNAL, "B02Kstar0nunu");
    CheckMaxNgamma(MC_dirname_SIGNAL, "B02Xsnunu");

    CheckMaxNgamma(MC_dirname_CHG, "root");
    CheckMaxNgamma(MC_dirname_MIX, "root");
    CheckMaxNgamma(MC_dirname_UUBAR, "root");
    CheckMaxNgamma(MC_dirname_DDBAR, "root");
    CheckMaxNgamma(MC_dirname_SSBAR, "root");
    CheckMaxNgamma(MC_dirname_CHARM, "root");
    
    AllocateMultiplicityPDFs(&Signal_multiplicity_p, "Signal", true);
    AllocateMultiplicityPDFs(&CHG_multiplicity_p, "CHG", true);
    AllocateMultiplicityPDFs(&MIX_multiplicity_p, "MIX", true);
    AllocateMultiplicityPDFs(&UUBAR_multiplicity_p, "UUBAR", true);
    AllocateMultiplicityPDFs(&DDBAR_multiplicity_p, "DDBAR", true);
    AllocateMultiplicityPDFs(&SSBAR_multiplicity_p, "SSBAR", true);
    AllocateMultiplicityPDFs(&CHARM_multiplicity_p, "CHARM", true);

    AllocateMultiplicityPDFs(&Signal_multiplicity_m, "Signal", false);
    AllocateMultiplicityPDFs(&CHG_multiplicity_m, "CHG", false);
    AllocateMultiplicityPDFs(&MIX_multiplicity_m, "MIX", false);
    AllocateMultiplicityPDFs(&UUBAR_multiplicity_m, "UUBAR", false);
    AllocateMultiplicityPDFs(&DDBAR_multiplicity_m, "DDBAR", false);
    AllocateMultiplicityPDFs(&SSBAR_multiplicity_m, "SSBAR", false);
    AllocateMultiplicityPDFs(&CHARM_multiplicity_m, "CHARM", false);

    GetMultiplicityPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_multiplicity_p, "Bplus", "SIGNAL", Scale_Kplus_test, "B2Knunu");
    GetMultiplicityPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_multiplicity_p, "Bplus", "SIGNAL", Scale_Kplusstar_test, "otherwise");
    GetMultiplicityPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_multiplicity_p, "Bplus", "SIGNAL", Scale_Xsu_nonresonant_test, "B2Xsnunu");
    GetMultiplicityPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_multiplicity_p, "Bzero", "SIGNAL", Scale_K0_test, "B02K0nunu");
    GetMultiplicityPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_multiplicity_p, "Bzero", "SIGNAL", Scale_K0star_test, "otherwise");
    GetMultiplicityPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_multiplicity_p, "Bzero", "SIGNAL", Scale_Xsd_nonresonant_test, "B02Xsnunu");

    GetMultiplicityPDFs(MC_dirname_CHG, "root", CHG_multiplicity_p, "Bplus", "CHG", Scale_CHG_test, "otherwise");
    GetMultiplicityPDFs(MC_dirname_MIX, "root", MIX_multiplicity_p, "Bzero", "MIX", Scale_MIX_test, "otherwise");
    GetMultiplicityPDFs(MC_dirname_UUBAR, "root", UUBAR_multiplicity_p, "Continuum", "UUBAR", Scale_UUBAR_test, "otherwise");
    GetMultiplicityPDFs(MC_dirname_DDBAR, "root", DDBAR_multiplicity_p, "Continuum", "DDBAR", Scale_DDBAR_test, "otherwise");
    GetMultiplicityPDFs(MC_dirname_SSBAR, "root", SSBAR_multiplicity_p, "Continuum", "SSBAR", Scale_SSBAR_test, "otherwise");
    GetMultiplicityPDFs(MC_dirname_CHARM, "root", CHARM_multiplicity_p, "Continuum", "CHARM", Scale_CHARM_test, "otherwise");

    int N_hist_multiplicity = PrintInfoMultiplicity();

    for (int i = 0; i < N_hist_multiplicity; i++) {
        GetNegativeChangePDFs(Signal_nominal, Signal_multiplicity_p[i], Signal_multiplicity_m[i]);
        GetNegativeChangePDFs(CHG_nominal, CHG_multiplicity_p[i], CHG_multiplicity_m[i]);
        GetNegativeChangePDFs(MIX_nominal, MIX_multiplicity_p[i], MIX_multiplicity_m[i]);
        GetNegativeChangePDFs(UUBAR_nominal, UUBAR_multiplicity_p[i], UUBAR_multiplicity_m[i]);
        GetNegativeChangePDFs(DDBAR_nominal, DDBAR_multiplicity_p[i], DDBAR_multiplicity_m[i]);
        GetNegativeChangePDFs(SSBAR_nominal, SSBAR_multiplicity_p[i], SSBAR_multiplicity_m[i]);
        GetNegativeChangePDFs(CHARM_nominal, CHARM_multiplicity_p[i], CHARM_multiplicity_m[i]);
    }

    // B -> [D -> X KL0] anything
    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_BtoDtoXKL_p, "Bplus", "SIGNAL", true, Scale_Kplus_test, "B2Knunu");
    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_BtoDtoXKL_p, "Bplus", "SIGNAL", true, Scale_Kplusstar_test, "otherwise");
    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_BtoDtoXKL_p, "Bplus", "SIGNAL", true, Scale_Xsu_nonresonant_test, "B2Xsnunu");
    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_BtoDtoXKL_p, "Bzero", "SIGNAL", true, Scale_K0_test, "B02K0nunu");
    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_BtoDtoXKL_p, "Bzero", "SIGNAL", true, Scale_K0star_test, "otherwise");
    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_BtoDtoXKL_p, "Bzero", "SIGNAL", true, Scale_Xsd_nonresonant_test, "B02Xsnunu");

    GetBtoDtoXKLPDFs(MC_dirname_CHG, "root", CHG_BtoDtoXKL_p, "Bplus", "CHG", true, Scale_CHG_test, "otherwise");
    GetBtoDtoXKLPDFs(MC_dirname_MIX, "root", MIX_BtoDtoXKL_p, "Bzero", "MIX", true, Scale_MIX_test, "otherwise");

    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_BtoDtoXKL_m, "Bplus", "SIGNAL", false, Scale_Kplus_test, "B2Knunu");
    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_BtoDtoXKL_m, "Bplus", "SIGNAL", false, Scale_Kplusstar_test, "otherwise");
    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_BtoDtoXKL_m, "Bplus", "SIGNAL", false, Scale_Xsu_nonresonant_test, "B2Xsnunu");
    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_BtoDtoXKL_m, "Bzero", "SIGNAL", false, Scale_K0_test, "B02K0nunu");
    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_BtoDtoXKL_m, "Bzero", "SIGNAL", false, Scale_K0star_test, "otherwise");
    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_BtoDtoXKL_m, "Bzero", "SIGNAL", false, Scale_Xsd_nonresonant_test, "B02Xsnunu");

    GetBtoDtoXKLPDFs(MC_dirname_CHG, "root", CHG_BtoDtoXKL_m, "Bplus", "CHG", false, Scale_CHG_test, "otherwise");
    GetBtoDtoXKLPDFs(MC_dirname_MIX, "root", MIX_BtoDtoXKL_m, "Bzero", "MIX", false, Scale_MIX_test, "otherwise");

    // calculate all uncorrelated pdfs
    ClearHist(Signal_all_uncorrelated, CHG_all_uncorrelated, MIX_all_uncorrelated, UUBAR_all_uncorrelated, DDBAR_all_uncorrelated, SSBAR_all_uncorrelated, CHARM_all_uncorrelated);

    AddSQRTHist(Signal_all_uncorrelated, Signal_FEI_uncorrelated, RarityBins);
    AddSQRTHist(CHG_all_uncorrelated, CHG_FEI_uncorrelated, RarityBins);
    AddSQRTHist(MIX_all_uncorrelated, MIX_FEI_uncorrelated, RarityBins);

    AddSQRTHist(Signal_all_uncorrelated, Signal_KID_uncorrelated, RarityBins);
    AddSQRTHist(CHG_all_uncorrelated, CHG_KID_uncorrelated, RarityBins);
    AddSQRTHist(MIX_all_uncorrelated, MIX_KID_uncorrelated, RarityBins);
    AddSQRTHist(UUBAR_all_uncorrelated, UUBAR_KID_uncorrelated, RarityBins);
    AddSQRTHist(DDBAR_all_uncorrelated, DDBAR_KID_uncorrelated, RarityBins);
    AddSQRTHist(SSBAR_all_uncorrelated, SSBAR_KID_uncorrelated, RarityBins);
    AddSQRTHist(CHARM_all_uncorrelated, CHARM_KID_uncorrelated, RarityBins);

    AddSQRTHist(Signal_all_uncorrelated, Signal_PID_uncorrelated, RarityBins);
    AddSQRTHist(CHG_all_uncorrelated, CHG_PID_uncorrelated, RarityBins);
    AddSQRTHist(MIX_all_uncorrelated, MIX_PID_uncorrelated, RarityBins);
    AddSQRTHist(UUBAR_all_uncorrelated, UUBAR_PID_uncorrelated, RarityBins);
    AddSQRTHist(DDBAR_all_uncorrelated, DDBAR_PID_uncorrelated, RarityBins);
    AddSQRTHist(SSBAR_all_uncorrelated, SSBAR_PID_uncorrelated, RarityBins);
    AddSQRTHist(CHARM_all_uncorrelated, CHARM_PID_uncorrelated, RarityBins);

    AddSQRTHist(Signal_all_uncorrelated, Signal_BR_uncorrelated, RarityBins);
    AddSQRTHist(CHG_all_uncorrelated, CHG_BR_uncorrelated, RarityBins);
    AddSQRTHist(MIX_all_uncorrelated, MIX_BR_uncorrelated, RarityBins);

    AddSQRTHist(Signal_all_uncorrelated, Signal_pi0_uncorrelated, RarityBins);
    AddSQRTHist(CHG_all_uncorrelated, CHG_pi0_uncorrelated, RarityBins);
    AddSQRTHist(MIX_all_uncorrelated, MIX_pi0_uncorrelated, RarityBins);
    AddSQRTHist(UUBAR_all_uncorrelated, UUBAR_pi0_uncorrelated, RarityBins);
    AddSQRTHist(DDBAR_all_uncorrelated, DDBAR_pi0_uncorrelated, RarityBins);
    AddSQRTHist(SSBAR_all_uncorrelated, SSBAR_pi0_uncorrelated, RarityBins);
    AddSQRTHist(CHARM_all_uncorrelated, CHARM_pi0_uncorrelated, RarityBins);
    /* ====================================== */



    /* ====================================== */
    // define path for Data
    const char* DATA_dirname_Knunu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Knunu";
    const char* DATA_dirname_Kstarnunu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Kstarnunu";
    const char* DATA_dirname_Xsununu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Xsnunu";
    const char* DATA_dirname_K0nunu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0nunu";
    const char* DATA_dirname_K0starnunu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0starnunu";
    const char* DATA_dirname_Xsdnunu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02Xsnunu";

    const char* DATA_dirname_CHG = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/CHG_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_MIX = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/MIX_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_UUBAR = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/UUBAR_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_DDBAR = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/DDBAR_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_SSBAR = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SSBAR_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_CHARM = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/CHARM_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    /* ====================================== */



    /* ====================================== */
    // Define histograms for Data
    TH1D* total_DATA = new TH1D("total_DATA", "total_DATA", RarityBins, BinMIN, BinMAX);
    /* ====================================== */



    /* ====================================== */
    // Save data
    GetNominalPDFs(DATA_dirname_Knunu, "root", total_DATA, "Bplus", "DATA", Scale_Kplus_test);
    GetNominalPDFs(DATA_dirname_Kstarnunu, "root", total_DATA, "Bplus", "DATA", Scale_Kplusstar_test);
    GetNominalPDFs(DATA_dirname_Xsununu, "root", total_DATA, "Bplus", "DATA", Scale_Xsu_nonresonant_test);
    GetNominalPDFs(DATA_dirname_K0nunu, "root", total_DATA, "Bzero", "DATA", Scale_K0_test);
    GetNominalPDFs(DATA_dirname_K0starnunu, "root", total_DATA, "Bzero", "DATA", Scale_K0star_test);
    GetNominalPDFs(DATA_dirname_Xsdnunu, "root", total_DATA, "Bzero", "DATA", Scale_Xsd_nonresonant_test);

    GetNominalPDFs(DATA_dirname_CHG, "root", total_DATA, "Bplus", "DATA", 1.0);
    GetNominalPDFs(DATA_dirname_MIX, "root", total_DATA, "Bzero", "DATA", 1.0);
    GetNominalPDFs(DATA_dirname_UUBAR, "root", total_DATA, "Continuum", "DATA", 1.0);
    GetNominalPDFs(DATA_dirname_DDBAR, "root", total_DATA, "Continuum", "DATA", 1.0);
    GetNominalPDFs(DATA_dirname_SSBAR, "root", total_DATA, "Continuum", "DATA", 1.0);
    GetNominalPDFs(DATA_dirname_CHARM, "root", total_DATA, "Continuum", "DATA", 1.0);
    /* ====================================== */



    /* ====================================== */
    // Save histograms
    TFile* file = new TFile("PDFandDATA.root", "RECREATE");

    Signal_nominal->Write();
    CHG_nominal->Write();
    MIX_nominal->Write();
    UUBAR_nominal->Write();
    DDBAR_nominal->Write();
    SSBAR_nominal->Write();
    CHARM_nominal->Write();

    // track uncertainty
    Signal_track_p->Write();
    CHG_track_p->Write();
    MIX_track_p->Write();
    UUBAR_track_p->Write();
    DDBAR_track_p->Write();
    SSBAR_track_p->Write();
    CHARM_track_p->Write();
    Signal_track_m->Write();
    CHG_track_m->Write();
    MIX_track_m->Write();
    UUBAR_track_m->Write();
    DDBAR_track_m->Write();
    SSBAR_track_m->Write();
    CHARM_track_m->Write();

    // KS0 uncertainty
    Signal_KS0_p->Write();
    CHG_KS0_p->Write();
    MIX_KS0_p->Write();
    UUBAR_KS0_p->Write();
    DDBAR_KS0_p->Write();
    SSBAR_KS0_p->Write();
    CHARM_KS0_p->Write();
    Signal_KS0_m->Write();
    CHG_KS0_m->Write();
    MIX_KS0_m->Write();
    UUBAR_KS0_m->Write();
    DDBAR_KS0_m->Write();
    SSBAR_KS0_m->Write();
    CHARM_KS0_m->Write();

    // FEI uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_FEI; i++) {
        Signal_FEI_correlated[i]->Write();
        CHG_FEI_correlated[i]->Write();
        MIX_FEI_correlated[i]->Write();
    }

    // FEI uncertainty (uncorrelated)
    Signal_FEI_uncorrelated->Write();
    CHG_FEI_uncorrelated->Write();
    MIX_FEI_uncorrelated->Write();

    // Kaon PID uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_KID; i++) {
        Signal_KID_correlated[i]->Write();
        CHG_KID_correlated[i]->Write();
        MIX_KID_correlated[i]->Write();
        UUBAR_KID_correlated[i]->Write();
        DDBAR_KID_correlated[i]->Write();
        SSBAR_KID_correlated[i]->Write();
        CHARM_KID_correlated[i]->Write();
    }

    // Kaon PID uncertainty (uncorrelated)
    Signal_KID_uncorrelated->Write();
    CHG_KID_uncorrelated->Write();
    MIX_KID_uncorrelated->Write();
    UUBAR_KID_uncorrelated->Write();
    DDBAR_KID_uncorrelated->Write();
    SSBAR_KID_uncorrelated->Write();
    CHARM_KID_uncorrelated->Write();

    // Pion PID uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_PID; i++) {
        Signal_PID_correlated[i]->Write();
        CHG_PID_correlated[i]->Write();
        MIX_PID_correlated[i]->Write();
        UUBAR_PID_correlated[i]->Write();
        DDBAR_PID_correlated[i]->Write();
        SSBAR_PID_correlated[i]->Write();
        CHARM_PID_correlated[i]->Write();
    }

    // Pion PID uncertainty (uncorrelated)
    Signal_PID_uncorrelated->Write();
    CHG_PID_uncorrelated->Write();
    MIX_PID_uncorrelated->Write();
    UUBAR_PID_uncorrelated->Write();
    DDBAR_PID_uncorrelated->Write();
    SSBAR_PID_uncorrelated->Write();
    CHARM_PID_uncorrelated->Write();

    // BB BR uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_BR; i++) {
        Signal_BR_correlated[i]->Write();
        CHG_BR_correlated[i]->Write();
        MIX_BR_correlated[i]->Write();
    }

    // BB BR uncertainty (uncorrelated)
    Signal_BR_uncorrelated->Write();
    CHG_BR_uncorrelated->Write();
    MIX_BR_uncorrelated->Write();

    // pi0 uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_pi0; i++) {
        Signal_pi0_correlated[i]->Write();
        CHG_pi0_correlated[i]->Write();
        MIX_pi0_correlated[i]->Write();
        UUBAR_pi0_correlated[i]->Write();
        DDBAR_pi0_correlated[i]->Write();
        SSBAR_pi0_correlated[i]->Write();
        CHARM_pi0_correlated[i]->Write();
    }

    // pi0 uncertainty (uncorrelated)
    Signal_pi0_uncorrelated->Write();
    CHG_pi0_uncorrelated->Write();
    MIX_pi0_uncorrelated->Write();
    UUBAR_pi0_uncorrelated->Write();
    DDBAR_pi0_uncorrelated->Write();
    SSBAR_pi0_uncorrelated->Write();
    CHARM_pi0_uncorrelated->Write();

    // K nu nubar form factor
    Signal_Kff1_p->Write();
    Signal_Kff1_m->Write();
    Signal_Kff2_p->Write();
    Signal_Kff2_m->Write();
    Signal_Kff3_p->Write();
    Signal_Kff3_m->Write();

    // Kstar nu nubar form factor
    Signal_Kstarff1_p->Write();
    Signal_Kstarff1_m->Write();
    Signal_Kstarff2_p->Write();
    Signal_Kstarff2_m->Write();
    Signal_Kstarff3_p->Write();
    Signal_Kstarff3_m->Write();
    Signal_Kstarff4_p->Write();
    Signal_Kstarff4_m->Write();
    Signal_Kstarff5_p->Write();
    Signal_Kstarff5_m->Write();
    Signal_Kstarff6_p->Write();
    Signal_Kstarff6_m->Write();
    Signal_Kstarff7_p->Write();
    Signal_Kstarff7_m->Write();
    Signal_Kstarff8_p->Write();
    Signal_Kstarff8_m->Write();
    Signal_Kstarff9_p->Write();
    Signal_Kstarff9_m->Write();

    // faction
    Signal_Kfrac_p->Write();
    Signal_Kfrac_m->Write();
    Signal_Kstarfrac_p->Write();
    Signal_Kstarfrac_m->Write();

    // fragmentation uncertainty
    for (unsigned int i = 0; i < Signal_Fragmentaions.size(); i++) Signal_Fragmentaions.at(i)->Write();

    // pf uncertainty
    Signal_pf_p->Write();
    Signal_pf_m->Write();

    // mb uncertainty
    Signal_mb_p->Write();
    Signal_mb_m->Write();

    // transition uncertainty
    Signal_transition_p->Write();
    Signal_transition_m->Write();

    // fixed mKstar uncertainty
    Signal_mKstar_p->Write();
    Signal_mKstar_m->Write();

    // BDTc uncertainty
    //CHG_BDTc_p->Write();
    //MIX_BDTc_p->Write();
    UUBAR_BDTc_p->Write();
    DDBAR_BDTc_p->Write();
    SSBAR_BDTc_p->Write();
    CHARM_BDTc_p->Write();
    //CHG_BDTc_m->Write();
    //MIX_BDTc_m->Write();
    UUBAR_BDTc_m->Write();
    DDBAR_BDTc_m->Write();
    SSBAR_BDTc_m->Write();
    CHARM_BDTc_m->Write();

    // B2Knn BR uncertainty
    CHG_Knn_p->Write();
    MIX_Knn_p->Write();
    CHG_Knn_m->Write();
    MIX_Knn_m->Write();
    CHG_Kstarnn_p->Write();
    MIX_Kstarnn_p->Write();
    CHG_Kstarnn_m->Write();
    MIX_Kstarnn_m->Write();
    CHG_K0nn_p->Write();
    MIX_K0nn_p->Write();
    CHG_K0nn_m->Write();
    MIX_K0nn_m->Write();
    CHG_K0starnn_p->Write();
    MIX_K0starnn_p->Write();
    CHG_K0starnn_m->Write();
    MIX_K0starnn_m->Write();

    // multiplicity uncertainty
    for (int i = 0; i < N_hist_multiplicity; i++) {
        Signal_multiplicity_p[i]->Write();
        CHG_multiplicity_p[i]->Write();
        MIX_multiplicity_p[i]->Write();
        UUBAR_multiplicity_p[i]->Write();
        DDBAR_multiplicity_p[i]->Write();
        SSBAR_multiplicity_p[i]->Write();
        CHARM_multiplicity_p[i]->Write();
        Signal_multiplicity_m[i]->Write();
        CHG_multiplicity_m[i]->Write();
        MIX_multiplicity_m[i]->Write();
        UUBAR_multiplicity_m[i]->Write();
        DDBAR_multiplicity_m[i]->Write();
        SSBAR_multiplicity_m[i]->Write();
        CHARM_multiplicity_m[i]->Write();
    }

    // B -> [D -> X KL0] anything uncertainties
    Signal_BtoDtoXKL_p->Write();
    CHG_BtoDtoXKL_p->Write();
    MIX_BtoDtoXKL_p->Write();
    Signal_BtoDtoXKL_m->Write();
    CHG_BtoDtoXKL_m->Write();
    MIX_BtoDtoXKL_m->Write();

    // all uncorrelated uncertainties
    Signal_all_uncorrelated->Write();
    CHG_all_uncorrelated->Write(); 
    MIX_all_uncorrelated->Write(); 
    UUBAR_all_uncorrelated->Write(); 
    DDBAR_all_uncorrelated->Write(); 
    SSBAR_all_uncorrelated->Write(); 
    CHARM_all_uncorrelated->Write();

    total_DATA->Write();

    file->Close();
    /* ====================================== */
}
