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

# include <algorithm>
# include <float.h>

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

// SIGNAL MC sample number before skimming for MC15rd
# define N_Kplus_validation_MC15rd 3962022.0
# define N_K0_validation_MC15rd 3962022.0
# define N_Kplusstar_validation_MC15rd 3962022.0
# define N_K0star_validation_MC15rd 3962022.0
# define N_Xsu_nonresonant_validation_MC15rd 15846594.0
# define N_Xsd_nonresonant_validation_MC15rd 15846594.0

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

// scale factor for SIGNAL MC sample until LS1 for MC15rd
# define Scale_Kplus_validation_MC15rd (N_Kplus_nunubar_LS1/N_Kplus_validation_MC15rd)
# define Scale_Kplusstar_validation_MC15rd (N_Kplusstar_nunubar_LS1/N_Kplusstar_validation_MC15rd)
# define Scale_Xsu_nonresonant_validation_MC15rd (N_Xsu_nonresonant_nunubar_LS1/N_Xsu_nonresonant_validation_MC15rd)
# define Scale_K0_validation_MC15rd (N_K0_nunubar_LS1/N_K0_validation_MC15rd)
# define Scale_K0star_validation_MC15rd (N_K0star_nunubar_LS1/N_K0star_validation_MC15rd)
# define Scale_Xsd_nonresonant_validation_MC15rd (N_Xsd_nunubar_LS1/N_Xsd_nonresonant_validation_MC15rd)

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
# define Scale_CHG_validation ((N_BB_LS1* (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab))
# define Scale_MIX_validation ((N_BB_LS1* (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab))
# define Scale_UUBAR_validation (0.361673)
# define Scale_DDBAR_validation (0.361673)
# define Scale_SSBAR_validation (0.361673)
# define Scale_CHARM_validation (0.361673)

// BKG MC sample number for MC15rd
# define N_CHG_validation_MC15rd 785108449.0 // 1458.959/fb
# define N_MIX_validation_MC15rd 741492304.0 // 1458.959/fb
# define N_UUBAR_validation_MC15rd 2306265848.0 // 1458.959/fb
# define N_DDBAR_validation_MC15rd 576209482.0 // 1458.959/fb
# define N_SSBAR_validation_MC15rd 526874294.0 // 1458.959/fb
# define N_CHARM_validation_MC15rd 1889822323.0 // 1458.959/fb

# define Scale_CHG_validation_MC15rd (0.361673/1.458959)
# define Scale_MIX_validation_MC15rd (0.361673/1.458959)
# define Scale_UUBAR_validation_MC15rd (0.361673/1.458959)
# define Scale_DDBAR_validation_MC15rd (0.361673/1.458959)
# define Scale_SSBAR_validation_MC15rd (0.361673/1.458959)
# define Scale_CHARM_validation_MC15rd (0.361673/1.458959)

double ObtainWeight(const char* type, const char* MC_version, const char* category, std::string filename);

# define KS0_rel_uncertainty 0.6 // %/cm
# define track_rel_uncertainty 0.24 // %
# define Kaon_PID_max_uncertainty 0.1 // not percentage. relative uncertainty

# define Nvar_num 1

# define CAL 1.0
# define CAL_qq 1.0
# define Stream 0.25

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

    // it is MC15ri correction factor. It should be fixed
    static constexpr double FEI_cal_Bc_MC15rd[FEI_cal_Bc_num] = { 1.04, 0.79, 0.69, 0.56, 0.97, 0.95, 0.74, 0.57, 0.91, 0.51, 0.34, 0.59 };
    static constexpr double FEI_cal_Bc_uncertainty_MC15rd[FEI_cal_Bc_num] = { 0.03, 0.03, 0.05, 0.11, 0.03, 0.03, 0.02, 0.06, 0.1, 0.13, 0.07, 0.02 }; // not relative uncertainty. absolute uncertainty
    static constexpr double FEI_cal_Bc_modeID_MC15rd[FEI_cal_Bc_num] = { 0.0, 1.0, 3.0, 4.0, 15.0, 16.0, 18.0, 19.0, 23.0, 24.0, 30.0, -1.0 };
    static constexpr double FEI_cal_B0_MC15rd[FEI_cal_B0_num] = { 1.16, 0.94, 0.81, 0.79, 0.99, 1.03, 0.67, 0.66, 0.69, 0.49, 0.79 };
    static constexpr double FEI_cal_B0_uncertainty_MC15rd[FEI_cal_B0_num] = { 0.04, 0.05, 0.06, 0.02, 0.03, 0.06, 0.02, 0.03, 0.02, 0.02, 0.12 }; // not relative uncertainty. absolute uncertainty
    static constexpr double FEI_cal_B0_modeID_MC15rd[FEI_cal_B0_num] = { 0.0, 1.0, 3.0, 4.0, 5.0, 15.0, 16.0, 18.0, 19.0, 26.0, -1.0 };
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
    else if (type == "MC15rd") {
        if (UpsilonID > -0.5 && UpsilonID < 0.5) { // charged
            for (int i = 0; i < FEI_cal_Bc_num - 1; i++) {
                if (BtagID > FEI_cal_Bc_modeID_MC15rd[i] - 0.5 && BtagID < FEI_cal_Bc_modeID_MC15rd[i] + 0.5) return FEI_cal_Bc_MC15rd[i];
            }
            return FEI_cal_Bc_MC15rd[FEI_cal_Bc_num - 1];
        }
        else if (UpsilonID > 0.5 && UpsilonID < 1.5) { // mixed
            for (int i = 0; i < FEI_cal_B0_num - 1; i++) {
                if (BtagID > FEI_cal_B0_modeID_MC15rd[i] - 0.5 && BtagID < FEI_cal_B0_modeID_MC15rd[i] + 0.5) return FEI_cal_B0_MC15rd[i];
            }
            return FEI_cal_B0_MC15rd[FEI_cal_B0_num - 1];
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
        else return FEI_cal_B0_MC15ri[index];
    }
    else if (type == "MC15rd") {
        if (IsItCharged) return FEI_cal_Bc_MC15rd[index];
        else return FEI_cal_B0_MC15rd[index];
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
    else if (type == "MC15rd") {
        if (UpsilonID > -0.5 && UpsilonID < 0.5) { // charged
            for (int i = 0; i < FEI_cal_Bc_num - 1; i++) {
                if (BtagID > FEI_cal_Bc_modeID_MC15rd[i] - 0.5 && BtagID < FEI_cal_Bc_modeID_MC15rd[i] + 0.5) return FEI_cal_Bc_uncertainty_MC15rd[i];
            }
            return FEI_cal_Bc_uncertainty_MC15rd[FEI_cal_Bc_num - 1];
        }
        else if (UpsilonID > 0.5 && UpsilonID < 1.5) { // mixed
            for (int i = 0; i < FEI_cal_B0_num - 1; i++) {
                if (BtagID > FEI_cal_B0_modeID_MC15rd[i] - 0.5 && BtagID < FEI_cal_B0_modeID_MC15rd[i] + 0.5) return FEI_cal_B0_uncertainty_MC15rd[i];
            }
            return FEI_cal_B0_uncertainty_MC15rd[FEI_cal_B0_num - 1];
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
        else return FEI_cal_B0_modeID_MC15ri[index];
    }
    else if (type == "MC15rd") {
        if (IsItCharged) return FEI_cal_Bc_modeID_MC15rd[index];
        else return FEI_cal_B0_modeID_MC15rd[index];
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

    double PID_correction_MC15rd[4][N_PID_syst]; // K-true, K-mis, pi-true, pi-miss
    double PID_correction_stat_uncer_MC15rd[4][N_PID_syst]; // K-true, K-mis, pi-true, pi-miss
    double PID_correction_sys_uncer_MC15rd[4][N_PID_syst]; // K-true, K-mis, pi-true, pi-miss
    double PID_correction_uncer_MC15rd[4][N_PID_syst]; // K-true, K-mis, pi-true, pi-miss

    void ReadPIDFile_MC15ri();
    void ReadPIDFile_MC15rd();
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
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < N_PID_syst; j++) {
            PID_correction_MC15rd[i][j] = 0.0;
            PID_correction_stat_uncer_MC15rd[i][j] = 0.0;
            PID_correction_sys_uncer_MC15rd[i][j] = 0.0;
            PID_correction_uncer_MC15rd[i][j] = 0.0;
        }
    }

    ReadPIDFile_MC15ri();
    ReadPIDFile_MC15rd();
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

void Corrector_PID::ReadPIDFile_MC15rd() {
    const char* KID_true_file = "/home/belle2/junewoo/storage_b1/bsub/systematic/MC15rd_PID/KaonEff.csv";
    const char* KID_mis_file = "/home/belle2/junewoo/storage_b1/bsub/systematic/MC15rd_PID/Kaonmis.csv";
    const char* PID_true_file = "/home/belle2/junewoo/storage_b1/bsub/systematic/MC15rd_PID/PionEff.csv";
    const char* PID_mis_file = "/home/belle2/junewoo/storage_b1/bsub/systematic/MC15rd_PID/Pionmis.csv";

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
        fscanf(fp_KID_true, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,kaonID\n", &temp_p_min, &temp_p_max, &temp_cosTheta_min, &temp_cosTheta_max, &PID_correction_MC15rd[0][i], &PID_correction_stat_uncer_MC15rd[0][i], &temp_data_MC_uncertainty_stat_dn, &PID_correction_sys_uncer_MC15rd[0][i], &temp_data_MC_uncertainty_sys_dn, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_threshold);
        fscanf(fp_KID_mis, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,kaonID\n", &temp_p_min, &temp_p_max, &temp_cosTheta_min, &temp_cosTheta_max, &PID_correction_MC15rd[1][i], &PID_correction_stat_uncer_MC15rd[1][i], &temp_data_MC_uncertainty_stat_dn, &PID_correction_sys_uncer_MC15rd[1][i], &temp_data_MC_uncertainty_sys_dn, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_threshold);
        fscanf(fp_PID_true, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,pionID\n", &temp_p_min, &temp_p_max, &temp_cosTheta_min, &temp_cosTheta_max, &PID_correction_MC15rd[2][i], &PID_correction_stat_uncer_MC15rd[2][i], &temp_data_MC_uncertainty_stat_dn, &PID_correction_sys_uncer_MC15rd[2][i], &temp_data_MC_uncertainty_sys_dn, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_threshold);
        fscanf(fp_PID_mis, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,pionID\n", &temp_p_min, &temp_p_max, &temp_cosTheta_min, &temp_cosTheta_max, &PID_correction_MC15rd[3][i], &PID_correction_stat_uncer_MC15rd[3][i], &temp_data_MC_uncertainty_stat_dn, &PID_correction_sys_uncer_MC15rd[3][i], &temp_data_MC_uncertainty_sys_dn, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_threshold);

        if ((std::abs(PID_correction_MC15rd[0][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer_MC15rd[0][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer_MC15rd[0][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction_MC15rd[0][i]) > 10000.0)) {
            PID_correction_MC15rd[0][i] = 1.0;
            PID_correction_stat_uncer_MC15rd[0][i] = 0.0;
            PID_correction_sys_uncer_MC15rd[0][i] = 0.0;
        }
        if ((std::abs(PID_correction_MC15rd[1][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer_MC15rd[1][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer_MC15rd[1][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction_MC15rd[1][i]) > 10000.0)) {
            PID_correction_MC15rd[1][i] = 1.0;
            PID_correction_stat_uncer_MC15rd[1][i] = 0.0;
            PID_correction_sys_uncer_MC15rd[1][i] = 0.0;
        }
        if ((std::abs(PID_correction_MC15rd[2][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer_MC15rd[2][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer_MC15rd[2][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction_MC15rd[2][i]) > 10000.0)) {
            PID_correction_MC15rd[2][i] = 1.0;
            PID_correction_stat_uncer_MC15rd[2][i] = 0.0;
            PID_correction_sys_uncer_MC15rd[2][i] = 0.0;
        }
        if ((std::abs(PID_correction_MC15rd[3][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer_MC15rd[3][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer_MC15rd[3][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction_MC15rd[3][i]) > 10000.0)) {
            PID_correction_MC15rd[3][i] = 1.0;
            PID_correction_stat_uncer_MC15rd[3][i] = 0.0;
            PID_correction_sys_uncer_MC15rd[3][i] = 0.0;
        }

        PID_correction_uncer_MC15rd[0][i] = std::sqrt(PID_correction_stat_uncer_MC15rd[0][i] * PID_correction_stat_uncer_MC15rd[0][i] + PID_correction_sys_uncer_MC15rd[0][i] * PID_correction_sys_uncer_MC15rd[0][i]);
        PID_correction_uncer_MC15rd[1][i] = std::sqrt(PID_correction_stat_uncer_MC15rd[1][i] * PID_correction_stat_uncer_MC15rd[1][i] + PID_correction_sys_uncer_MC15rd[1][i] * PID_correction_sys_uncer_MC15rd[1][i]);
        PID_correction_uncer_MC15rd[2][i] = std::sqrt(PID_correction_stat_uncer_MC15rd[2][i] * PID_correction_stat_uncer_MC15rd[2][i] + PID_correction_sys_uncer_MC15rd[2][i] * PID_correction_sys_uncer_MC15rd[2][i]);
        PID_correction_uncer_MC15rd[3][i] = std::sqrt(PID_correction_stat_uncer_MC15rd[3][i] * PID_correction_stat_uncer_MC15rd[3][i] + PID_correction_sys_uncer_MC15rd[3][i] * PID_correction_sys_uncer_MC15rd[3][i]);
    }

    PID_correction_MC15rd[0][N_PID_syst - 1] = 1.0;
    PID_correction_stat_uncer_MC15rd[0][N_PID_syst - 1] = 0.0;
    PID_correction_sys_uncer_MC15rd[0][N_PID_syst - 1] = 0.0;
    PID_correction_uncer_MC15rd[0][N_PID_syst - 1] = 0.0;

    PID_correction_MC15rd[1][N_PID_syst - 1] = 1.0;
    PID_correction_stat_uncer_MC15rd[1][N_PID_syst - 1] = 0.0;
    PID_correction_sys_uncer_MC15rd[1][N_PID_syst - 1] = 0.0;
    PID_correction_uncer_MC15rd[1][N_PID_syst - 1] = 0.0;

    PID_correction_MC15rd[2][N_PID_syst - 1] = 1.0;
    PID_correction_stat_uncer_MC15rd[2][N_PID_syst - 1] = 0.0;
    PID_correction_sys_uncer_MC15rd[2][N_PID_syst - 1] = 0.0;
    PID_correction_uncer_MC15rd[2][N_PID_syst - 1] = 0.0;

    PID_correction_MC15rd[3][N_PID_syst - 1] = 1.0;
    PID_correction_stat_uncer_MC15rd[3][N_PID_syst - 1] = 0.0;
    PID_correction_sys_uncer_MC15rd[3][N_PID_syst - 1] = 0.0;
    PID_correction_uncer_MC15rd[3][N_PID_syst - 1] = 0.0;

    fclose(fp_KID_true);
    fclose(fp_KID_mis);
    fclose(fp_PID_true);
    fclose(fp_PID_mis);
}

double Corrector_PID::GetCorrectionFactor(int PID_type, int bin_PID, std::string type) {
    if (type == "MC15ri") return PID_correction_MC15ri[PID_type][bin_PID];
    else if (type == "MC15rd") return PID_correction_MC15rd[PID_type][bin_PID];
    else {
        printf("[Corrector_pi0] Invalid type!\n");
        exit(1);
    }
}

double Corrector_PID::GetUncertainty(int PID_type, int bin_PID, std::string type) {
    if (type == "MC15ri") return PID_correction_uncer_MC15ri[PID_type][bin_PID];
    else if (type == "MC15rd") return PID_correction_uncer_MC15rd[PID_type][bin_PID];
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

    double PID_fakeE_correction_MC15rd[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double PID_fakeE_uncer_MC15rd[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double PID_fakeMU_correction_MC15rd[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double PID_fakeMU_uncer_MC15rd[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    void ReadPIDFile_MC15ri();
    void ReadPIDFile_MC15rd();
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
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < N_fakeE_syst; j++) {
            PID_fakeE_correction_MC15rd[i][j] = 0.0;
            PID_fakeE_uncer_MC15rd[i][j] = 0.0;
        }
    }
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < N_fakeMU_syst; j++) {
            PID_fakeMU_correction_MC15rd[i][j] = 0.0;
            PID_fakeMU_uncer_MC15rd[i][j] = 0.0;
        }
    }

    ReadPIDFile_MC15ri();
    ReadPIDFile_MC15rd();
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

void Corrector_FakePID::ReadPIDFile_MC15rd() {
    // initialization
    for (int i = 0; i < N_fakeE_syst; i++) {
        PID_fakeE_correction_MC15rd[0][i] = 1.0; //  K-, K+, pi-, pi+
        PID_fakeE_correction_MC15rd[1][i] = 1.0;
        PID_fakeE_correction_MC15rd[2][i] = 1.0;
        PID_fakeE_correction_MC15rd[3][i] = 1.0;

        PID_fakeE_uncer_MC15rd[0][i] = 0.0;
        PID_fakeE_uncer_MC15rd[1][i] = 0.0;
        PID_fakeE_uncer_MC15rd[2][i] = 0.0;
        PID_fakeE_uncer_MC15rd[3][i] = 0.0;
    }

    for (int i = 0; i < N_fakeMU_syst; i++) {
        PID_fakeMU_correction_MC15rd[0][i] = 1.0;
        PID_fakeMU_correction_MC15rd[1][i] = 1.0;
        PID_fakeMU_correction_MC15rd[2][i] = 1.0;
        PID_fakeMU_correction_MC15rd[3][i] = 1.0;

        PID_fakeMU_uncer_MC15rd[0][i] = 0.0;
        PID_fakeMU_uncer_MC15rd[1][i] = 0.0;
        PID_fakeMU_uncer_MC15rd[2][i] = 0.0;
        PID_fakeMU_uncer_MC15rd[3][i] = 0.0;
    }

    const char* K_fromE_file = "/home/belle2/junewoo/storage_b1/bsub/systematic/MC15rd_PID/kaonID_efficiency_electron.csv";
    const char* K_fromMU_file = "/home/belle2/junewoo/storage_b1/bsub/systematic/MC15rd_PID/kaonID_efficiency_muon.csv";
    const char* pi_fromE_file = "/home/belle2/junewoo/storage_b1/bsub/systematic/MC15rd_PID/pionID_efficiency_electron.csv";
    const char* pi_fromMU_file = "/home/belle2/junewoo/storage_b1/bsub/systematic/MC15rd_PID/pionID_efficiency_muon.csv";

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

            PID_fakeE_correction_MC15rd[1][bin] = temp_data_MC_ratio;
            PID_fakeE_uncer_MC15rd[1][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeE_correction_MC15rd[1][bin]) < MyEPSILON) {
                PID_fakeE_correction_MC15rd[1][bin] = 1.0;
                PID_fakeE_uncer_MC15rd[1][bin] = 0.0;
            }

        }
        else if (temp_charge == '-') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeE_correction_MC15rd[0][bin] = temp_data_MC_ratio;
            PID_fakeE_uncer_MC15rd[0][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeE_correction_MC15rd[0][bin]) < MyEPSILON) {
                PID_fakeE_correction_MC15rd[0][bin] = 1.0;
                PID_fakeE_uncer_MC15rd[0][bin] = 0.0;
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

            PID_fakeE_correction_MC15rd[3][bin] = temp_data_MC_ratio;
            PID_fakeE_uncer_MC15rd[3][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeE_correction_MC15rd[3][bin]) < MyEPSILON) {
                PID_fakeE_correction_MC15rd[3][bin] = 1.0;
                PID_fakeE_uncer_MC15rd[3][bin] = 0.0;
            }

        }
        else if (temp_charge == '-') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeE_correction_MC15rd[2][bin] = temp_data_MC_ratio;
            PID_fakeE_uncer_MC15rd[2][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeE_correction_MC15rd[2][bin]) < MyEPSILON) {
                PID_fakeE_correction_MC15rd[2][bin] = 1.0;
                PID_fakeE_uncer_MC15rd[2][bin] = 0.0;
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

            PID_fakeMU_correction_MC15rd[1][bin] = temp_data_MC_ratio;
            PID_fakeMU_uncer_MC15rd[1][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeMU_correction_MC15rd[1][bin]) < MyEPSILON) {
                PID_fakeMU_correction_MC15rd[1][bin] = 1.0;
                PID_fakeMU_uncer_MC15rd[1][bin] = 0.0;
            }

        }
        else if (temp_charge == '-') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeMU_correction_MC15rd[0][bin] = temp_data_MC_ratio;
            PID_fakeMU_uncer_MC15rd[0][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeMU_correction_MC15rd[0][bin]) < MyEPSILON) {
                PID_fakeMU_correction_MC15rd[0][bin] = 1.0;
                PID_fakeMU_uncer_MC15rd[0][bin] = 0.0;
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

            PID_fakeMU_correction_MC15rd[3][bin] = temp_data_MC_ratio;
            PID_fakeMU_uncer_MC15rd[3][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeMU_correction_MC15rd[3][bin]) < MyEPSILON) {
                PID_fakeMU_correction_MC15rd[3][bin] = 1.0;
                PID_fakeMU_uncer_MC15rd[3][bin] = 0.0;
            }

        }
        else if (temp_charge == '-') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeMU_correction_MC15rd[2][bin] = temp_data_MC_ratio;
            PID_fakeMU_uncer_MC15rd[2][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeMU_correction_MC15rd[2][bin]) < MyEPSILON) {
                PID_fakeMU_correction_MC15rd[2][bin] = 1.0;
                PID_fakeMU_uncer_MC15rd[2][bin] = 0.0;
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
    else if (type == "MC15rd") return PID_fakeE_correction_MC15rd[PID_type][bin_PID];
    else {
        printf("[Corrector_FakePID] Invalid type!\n");
        exit(1);
    }
}

double Corrector_FakePID::GetCorrectionFactorfakeMU(int PID_type, int bin_PID, std::string type) {
    if (type == "MC15ri") return PID_fakeMU_correction_MC15ri[PID_type][bin_PID];
    else if (type == "MC15rd") return PID_fakeMU_correction_MC15rd[PID_type][bin_PID];
    else {
        printf("[Corrector_FakePID] Invalid type!\n");
        exit(1);
    }
}

double Corrector_FakePID::GetUncertaintyfakeE(int PID_type, int bin_PID, std::string type) {
    if (type == "MC15ri") return PID_fakeE_uncer_MC15ri[PID_type][bin_PID];
    else if (type == "MC15rd") return PID_fakeE_uncer_MC15rd[PID_type][bin_PID];
    else {
        printf("[Corrector_FakePID] Invalid type!\n");
        exit(1);
    }
}

double Corrector_FakePID::GetUncertaintyfakeMU(int PID_type, int bin_PID, std::string type) {
    if (type == "MC15ri") return PID_fakeMU_uncer_MC15ri[PID_type][bin_PID];
    else if (type == "MC15rd") return PID_fakeMU_uncer_MC15rd[PID_type][bin_PID];
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

    // it is MC15ri correction factor. It should be fixed
    static constexpr double pi0_correction_MC15rd[N_pi0_syst] = { 0.917, 0.965, 0.988, 1.013, 1.042, 1.044, 1.011, 1.0 };
    static constexpr double pi0_stat_uncer_MC15rd[N_pi0_syst] = { 0.004, 0.004, 0.004, 0.005, 0.004, 0.005, 0.005, 0.0 };
    static constexpr double pi0_sys_uncer1_MC15rd[N_pi0_syst] = { 0.049, 0.036, 0.079, 0.058, 0.045, 0.041, 0.040, 0.0 };
    static constexpr double pi0_sys_uncer2_MC15rd[N_pi0_syst] = { 0.0, 0.0, 0.0, 0.0, 0.039, 0.051, 0.030, 0.0 };
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
    else if (type == "MC15rd") return pi0_correction_MC15rd[bin_pi0];
    else {
        printf("[Corrector_pi0] Invalid type!\n");
        exit(1);
    }
}

double Corrector_pi0::GetStatUncertainty(int bin_pi0, std::string type) {
    if (type == "MC15ri") return pi0_stat_uncer_MC15ri[bin_pi0];
    else if (type == "MC15rd") return pi0_stat_uncer_MC15rd[bin_pi0];
    else {
        printf("[Corrector_pi0] Invalid type!\n");
        exit(1);
    }
}

double Corrector_pi0::GetSystUncertainty1(int bin_pi0, std::string type) {
    if (type == "MC15ri") return pi0_sys_uncer1_MC15ri[bin_pi0];
    else if (type == "MC15rd") return pi0_sys_uncer1_MC15rd[bin_pi0];
    else {
        printf("[Corrector_pi0] Invalid type!\n");
        exit(1);
    }
}

double Corrector_pi0::GetSystUncertainty2(int bin_pi0, std::string type) {
    if (type == "MC15ri") return pi0_sys_uncer2_MC15ri[bin_pi0];
    else if (type == "MC15rd") return pi0_sys_uncer2_MC15rd[bin_pi0];
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

        if ((type == "MC15ri") || (type == "MC15rd")) {
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
            else {
                printf("[Corrector_Fragmentation::GetCorrectionFactor] unexpected syst type!\n");
                exit(1);
                return 1;
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
bool sorting(std::pair<double, double> a, std::pair<double, double> b) {

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
    Corrector_Multiplicity(const char* filename);
    double GetCorrectionFactor(double Ngamma);
};

Corrector_Multiplicity corrector_Multiplicity;

Corrector_Multiplicity::Corrector_Multiplicity() :
    CUTOFF(50.0)
{
    FILE* fp;

    // read multiplicity weights
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

Corrector_Multiplicity::Corrector_Multiplicity(const char* filename) :
    CUTOFF(50.0)
{
    FILE* fp;

    // read multiplicity weights
    fp = fopen(filename, "r");
    fscanf(fp, "%d\n", &NgammaMAX);
    weights_Ngamma = new TH1D(("weights_Ngamma_" + std::string(filename)).c_str(), ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);
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

bool hasEnding(std::string const& fullString, std::string const& ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
    }
    else {
        return false;
    }
}


void load_files(const char *dirname, std::vector<string>* names){
   TSystemDirectory dir(dirname, dirname);
   TList *files = dir.GetListOfFiles();
   if (files) {
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next())) {
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

void LetsFillNominal_ri(const char* dirname, std::vector<std::string> variable_names, std::vector<std::string> branch_names, std::vector<double> variable_values[Nvar_num], const char* included_string, std::vector<int>* numberings, std::vector<double>* weights, std::string SampleName, const char* type, std::string CorrectionType = "otherwise") {
    /*
    SampleName for Knn
    CHG
    MIX
    UUBAR
    DDBAR
    SSBAR
    CHARM
    SIGNAL
    */
    /*
    type
    Bplus
    Bzero
    Continuum
    */
    /*
    0: charged
    1: mixed
    2: uubar
    3: ddbar
    4: ssbar
    5: ccbar
    6: tautau
    7: mumu
    8: gg
    9: ee
    10: eeee
    11: eemumu
    12: llXX
    13: hhISR
    */

    double var[Nvar_num] = { 0.0 };
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

    double FEI_calibration_factor = -1;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;

    int Ngamma_v200_index = -1;
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

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        for (int k = 0; k < (int)variable_names.size(); k++) {
            if (branch_names.at(k) == std::string("Upsilon")) tree_upsilon->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else if (branch_names.at(k) == std::string("Bsig")) tree_Bsig->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else if (branch_names.at(k) == std::string("Btag")) tree_Btag->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else {
                printf("ERROR! \n");
                exit(1);
            }
        }

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID); // charged: 0, mixed: 1
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

        Ngamma_v200_index = std::find(variable_names.begin(), variable_names.end(), std::string("extraInfo__boNgammav200__bc")) - variable_names.begin();

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

            for (int k = 0; k < (int)variable_names.size(); k++) variable_values[k].push_back(var[k]);

            // FEI calibration factor
            if (strcmp(type, "Bplus") == 0) FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) FEI_calibration_factor = CAL_qq;

            // Fill numberings
            double weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            if (SampleName == "CHG") {
                numberings->push_back(0);
            }
            else if (SampleName == "MIX") {
                numberings->push_back(1);
            }
            else if (SampleName == "UUBAR") {
                numberings->push_back(2);
            }
            else if (SampleName == "DDBAR") {
                numberings->push_back(3);
            }
            else if (SampleName == "SSBAR") {
                numberings->push_back(4);
            }
            else if (SampleName == "CHARM") {
                numberings->push_back(5);
            }
            else if (SampleName == "SIGNAL") {
                numberings->push_back(14);
            }
            //else if (job_id >= 256846858 && job_id <= 256847295) numberings->push_back(6);
            //else if (job_id >= 256847296 && job_id <= 256847807) numberings->push_back(7);
            //else if (job_id >= 256847808 && job_id <= 256848291) numberings->push_back(8);
            //else if (job_id >= 256848292 && job_id <= 256848743) numberings->push_back(9);
            //else if (job_id >= 256848744 && job_id <= 256849128) numberings->push_back(10);
            //else if (job_id >= 256849129 && job_id <= 256849396) numberings->push_back(11);
            else {
                printf("undefined job id!\n");
                exit(1);
            }

            // Fill calibration factors
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

            // Multiplicity correction factor, it is not applied now. it is for a systematic uncertainty
            Ngamma_v200 = var[Ngamma_v200_index];
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if (SampleName == "CHG" || SampleName == "MIX" || SampleName == "SIGNAL") Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double total_weight = FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            weights->push_back(total_weight);

        }
        input_file->Close();

    }

}

void LetsFill(const char* dirname, std::vector<std::string> variable_names, std::vector<std::string> branch_names, std::vector<double> variable_values[Nvar_num], const char* included_string) {
    double var[Nvar_num] = { 0.0 };

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        for (int k = 0; k < (int)variable_names.size(); k++) {
            if (branch_names.at(k) == std::string("Upsilon")) tree_upsilon->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else if (branch_names.at(k) == std::string("Bsig")) tree_Bsig->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else if (branch_names.at(k) == std::string("Btag")) tree_Btag->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else {
                printf("ERROR! \n");
                exit(1);
            }
        }

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            for (int k = 0; k < (int)variable_names.size(); k++) variable_values[k].push_back(var[k]);
        }
        input_file->Close();

    }

}

void LetsFillNominal_ri_correction(const char* dirname, std::vector<std::string> variable_names, std::vector<std::string> branch_names, std::vector<double> variable_values[Nvar_num], const char* included_string, std::vector<int>* numberings, std::vector<double>* weights, std::string SampleName, const char* type, std::string CorrectionType = "otherwise", double NormFactor = 1.0) {
    /*
    SampleName for Knn
    CHG
    MIX
    UUBAR
    DDBAR
    SSBAR
    CHARM
    SIGNAL
    */
    /*
    type
    Bplus
    Bzero
    Continuum
    */
    /*
    0: charged
    1: mixed
    2: uubar
    3: ddbar
    4: ssbar
    5: ccbar
    6: tautau
    7: mumu
    8: gg
    9: ee
    10: eeee
    11: eemumu
    12: llXX
    13: hhISR
    */

    double var[Nvar_num] = { 0.0 };
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

    double FEI_calibration_factor = -1;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;

    int Ngamma_v200_index = -1;
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

    float BDTc = -1;
    double BDTc_correction = -1;

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

        for (int k = 0; k < (int)variable_names.size(); k++) {
            if (branch_names.at(k) == std::string("Upsilon")) tree_upsilon->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else if (branch_names.at(k) == std::string("Bsig")) tree_Bsig->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else if (branch_names.at(k) == std::string("Btag")) tree_Btag->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else {
                printf("ERROR! \n");
                exit(1);
            }
        }

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID); // charged: 0, mixed: 1
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

        Ngamma_v200_index = std::find(variable_names.begin(), variable_names.end(), std::string("extraInfo__boNgammav200__bc")) - variable_names.begin();

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

        tree_upsilon->SetBranchAddress("MVA_Continuum", &BDTc);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            // BDTc correction factor
            if (BDTc > (5.0 / 6.0)) BDTc_correction = 5.0;
            else BDTc_correction = (BDTc / (1.0 - BDTc));
            BDTc_correction = BDTc_correction * NormFactor;

            for (int k = 0; k < (int)variable_names.size(); k++) variable_values[k].push_back(var[k]);

            // FEI calibration factor
            if (strcmp(type, "Bplus") == 0) FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) FEI_calibration_factor = CAL_qq;

            // Fill numberings
            double weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            if (SampleName == "CHG") {
                numberings->push_back(0);
            }
            else if (SampleName == "MIX") {
                numberings->push_back(1);
            }
            else if (SampleName == "UUBAR") {
                numberings->push_back(2);
            }
            else if (SampleName == "DDBAR") {
                numberings->push_back(3);
            }
            else if (SampleName == "SSBAR") {
                numberings->push_back(4);
            }
            else if (SampleName == "CHARM") {
                numberings->push_back(5);
            }
            else if (SampleName == "SIGNAL") {
                numberings->push_back(14);
            }
            //else if (job_id >= 256846858 && job_id <= 256847295) numberings->push_back(6);
            //else if (job_id >= 256847296 && job_id <= 256847807) numberings->push_back(7);
            //else if (job_id >= 256847808 && job_id <= 256848291) numberings->push_back(8);
            //else if (job_id >= 256848292 && job_id <= 256848743) numberings->push_back(9);
            //else if (job_id >= 256848744 && job_id <= 256849128) numberings->push_back(10);
            //else if (job_id >= 256849129 && job_id <= 256849396) numberings->push_back(11);
            else {
                printf("undefined job id!\n");
                exit(1);
            }

            // Fill calibration factors
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

            // Multiplicity correction factor, it is not applied now. it is for a systematic uncertainty
            Ngamma_v200 = var[Ngamma_v200_index];
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if (SampleName == "CHG" || SampleName == "MIX" || SampleName == "SIGNAL") Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double total_weight = FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_BtoDtoXKL * BDTc_correction;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            weights->push_back(total_weight);

        }
        input_file->Close();

    }

}

typedef struct _Nevt {
    double NevtwithoutCorrection;
    double NevtwithCorrection;
} Nevt;

void NevtCount_ri(const char* dirname, std::string SampleName, Nevt* nevt, const char* included_string, const char* type, std::string CorrectionType = "otherwise") {
    /*
    SampleName for Knn
    CHG
    MIX
    UUBAR
    DDBAR
    SSBAR
    CHARM
    SIGNAL
    */
    /*
    type
    Bplus
    Bzero
    Continuum
    */
    /*
    0: charged
    1: mixed
    2: uubar
    3: ddbar
    4: ssbar
    5: ccbar
    6: tautau
    7: mumu
    8: gg
    9: ee
    10: eeee
    11: eemumu
    12: llXX
    13: hhISR
    */

    double var[Nvar_num] = { 0.0 };
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

    double FEI_calibration_factor = -1;

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

    float BDTc = -1;
    double BDTc_correction = -1;

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

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID); // charged: 0, mixed: 1
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

        tree_upsilon->SetBranchAddress("MVA_Continuum", &BDTc);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            // BDTc correction factor
            if (BDTc > (5.0 / 6.0)) BDTc_correction = 5.0;
            else BDTc_correction = (BDTc / (1.0 - BDTc));

            // FEI calibration factor
            if (strcmp(type, "Bplus") == 0) FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) FEI_calibration_factor = CAL_qq;

            // Fill numberings
            double weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));

            // Fill calibration factors
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

            // Multiplicity correction factor, it is not applied now. it is for a systematic uncertainty
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if (SampleName == "CHG" || SampleName == "MIX" || SampleName == "SIGNAL") Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double total_weight = FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + total_weight;
            nevt->NevtwithCorrection = nevt->NevtwithCorrection + total_weight * BDTc_correction;


        }
        input_file->Close();

    }

}

double ObtainWeight(const char* type, const char* MC_version, const char* category, std::string filename) {
    if (strcmp(MC_version, "data") == 0) return 1.0; // no weight if it is data no matter what other values are
    else if ((strcmp(MC_version, "MC15ri") == 0) && (strcmp(category, "test") == 0)) { // MC15ri test
        if ((strcmp(type, "SIGNAL") == 0)) {
            if (filename.find("B2Knunu") != std::string::npos) return Scale_Kplus_test;
            else if (filename.find("B2Kstarnunu") != std::string::npos) return Scale_Kplusstar_test;
            else if (filename.find("B2Xsnunu") != std::string::npos) return Scale_Xsu_nonresonant_test;
            else if (filename.find("B02K0nunu") != std::string::npos) return Scale_K0_test;
            else if (filename.find("B02Kstar0nunu") != std::string::npos) return Scale_K0star_test;
            else if (filename.find("B02Xsnunu") != std::string::npos) return Scale_Xsd_nonresonant_test;
            else { printf("[ObtainWeight] undefined type for SIGNAL\n"); exit(1); }
        }
        else if ((strcmp(type, "CHG") == 0)) return Scale_CHG_test;
        else if ((strcmp(type, "MIX") == 0)) return Scale_MIX_test;
        else if ((strcmp(type, "UUBAR") == 0)) return Scale_UUBAR_test;
        else if ((strcmp(type, "DDBAR") == 0)) return Scale_DDBAR_test;
        else if ((strcmp(type, "SSBAR") == 0)) return Scale_SSBAR_test;
        else if ((strcmp(type, "CHARM") == 0)) return Scale_CHARM_test;
        else { printf("[ObtainWeight] undefined type for MC15ri test\n"); exit(1); }
    }
    else if ((strcmp(MC_version, "MC15ri") == 0) && (strcmp(category, "train") == 0)) { // MC15ri train
        if ((strcmp(type, "SIGNAL") == 0)) {
            if (filename.find("B2Knunu") != std::string::npos) return Scale_Kplus_train;
            else if (filename.find("B2Kstarnunu") != std::string::npos) return Scale_Kplusstar_train;
            else if (filename.find("B2Xsnunu") != std::string::npos) return Scale_Xsu_nonresonant_train;
            else if (filename.find("B02K0nunu") != std::string::npos) return Scale_K0_train;
            else if (filename.find("B02Kstar0nunu") != std::string::npos) return Scale_K0star_train;
            else if (filename.find("B02Xsnunu") != std::string::npos) return Scale_Xsd_nonresonant_train;
            else { printf("[ObtainWeight] undefined type for SIGNAL\n"); exit(1); }
        }
        else if ((strcmp(type, "CHG") == 0)) return Scale_CHG_train;
        else if ((strcmp(type, "MIX") == 0)) return Scale_MIX_train;
        else if ((strcmp(type, "UUBAR") == 0)) return Scale_UUBAR_train;
        else if ((strcmp(type, "DDBAR") == 0)) return Scale_DDBAR_train;
        else if ((strcmp(type, "SSBAR") == 0)) return Scale_SSBAR_train;
        else if ((strcmp(type, "CHARM") == 0)) return Scale_CHARM_train;
        else { printf("[ObtainWeight] undefined type for MC15ri train\n"); exit(1); }
    }
    else if ((strcmp(MC_version, "MC15ri") == 0) && (strcmp(category, "validation") == 0)) { // MC15ri validation
        if ((strcmp(type, "SIGNAL") == 0)) {
            if (filename.find("CHG") != std::string::npos) return Scale_CHG_validation; // it is Jpsi signal analysis with generic CHG sample
            else if (filename.find("MIX") != std::string::npos) return Scale_MIX_validation; // it is Jpsi signal analysis with generic CHG sample
            return 1.0;
        }
        else if ((strcmp(type, "CHG") == 0)) return Scale_CHG_validation;
        else if ((strcmp(type, "MIX") == 0)) return Scale_MIX_validation;
        else if ((strcmp(type, "UUBAR") == 0)) return Scale_UUBAR_validation;
        else if ((strcmp(type, "DDBAR") == 0)) return Scale_DDBAR_validation;
        else if ((strcmp(type, "SSBAR") == 0)) return Scale_SSBAR_validation;
        else if ((strcmp(type, "CHARM") == 0)) return Scale_CHARM_validation;
        else { printf("[ObtainWeight] undefined type for MC15ri validation\n"); exit(1); }
    }
    else if ((strcmp(MC_version, "MC15rd") == 0) && (strcmp(category, "validation") == 0)) { // MC15rd validation
        if ((strcmp(type, "SIGNAL") == 0)) {
            if (filename.find("B2Knunu") != std::string::npos) return Scale_Kplus_validation_MC15rd;
            else if (filename.find("B2Kstarnunu") != std::string::npos) return Scale_Kplusstar_validation_MC15rd;
            else if (filename.find("B2Xsnunu") != std::string::npos) return Scale_Xsu_nonresonant_validation_MC15rd;
            else if (filename.find("B02K0nunu") != std::string::npos) return Scale_K0_validation_MC15rd;
            else if (filename.find("B02Kstar0nunu") != std::string::npos) return Scale_K0star_validation_MC15rd;
            else if (filename.find("B02Xsnunu") != std::string::npos) return Scale_Xsd_nonresonant_validation_MC15rd;
            else if (filename.find("CHG") != std::string::npos) return Scale_CHG_validation_MC15rd; // it is Jpsi signal analysis with generic CHG sample
            else if (filename.find("MIX") != std::string::npos) return Scale_MIX_validation_MC15rd; // it is Jpsi signal analysis with generic CHG sample
            return 1.0; // just data
        }
        else if ((strcmp(type, "CHG") == 0)) return Scale_CHG_validation_MC15rd;
        else if ((strcmp(type, "MIX") == 0)) return Scale_MIX_validation_MC15rd;
        else if ((strcmp(type, "UUBAR") == 0)) return Scale_UUBAR_validation_MC15rd;
        else if ((strcmp(type, "DDBAR") == 0)) return Scale_DDBAR_validation_MC15rd;
        else if ((strcmp(type, "SSBAR") == 0)) return Scale_SSBAR_validation_MC15rd;
        else if ((strcmp(type, "CHARM") == 0)) return Scale_CHARM_validation_MC15rd;
        else { printf("[ObtainWeight] undefined type for MC15ri validation\n"); exit(1); }
    }

    printf("[ObtainWeight] no matched case!\n");
    exit(1);
    return 1.0;
}

void THStack_plot_nominal() {

    Nevt nevt_CHG = { 0.0, 0.0 };
    Nevt nevt_MIX = { 0.0, 0.0 };
    Nevt nevt_UUBAR = { 0.0, 0.0 };
    Nevt nevt_DDBAR = { 0.0, 0.0 };
    Nevt nevt_SSBAR = { 0.0, 0.0 };
    Nevt nevt_CHARM = { 0.0, 0.0 };
    Nevt nevt_SIGNAL = { 0.0, 0.0 };

    const char* Nominal_MC_CHG_validation_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/CHG_analysis/validation_v000/final_output";
    const char* Nominal_MC_MIX_validation_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/MIX_analysis/validation_v000/final_output";
    const char* Nominal_MC_UUBAR_validation_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/UUBAR_analysis/validation_v000/final_output";
    const char* Nominal_MC_DDBAR_validation_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/DDBAR_analysis/validation_v000/final_output";
    const char* Nominal_MC_SSBAR_validation_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SSBAR_analysis/validation_v000/final_output";
    const char* Nominal_MC_CHARM_validation_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/CHARM_analysis/validation_v000/final_output";
    const char* Nominal_MC_SIGNAL_validation_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/validation_v000/final_output";

    const char* Nominal_data_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD_LS_data/SIGNAL_analysis/validation_v000/final_output";

    NevtCount_ri(Nominal_MC_CHG_validation_dirname, "CHG", &nevt_CHG, "root", "Bplus", "otherwise");
    NevtCount_ri(Nominal_MC_MIX_validation_dirname, "MIX", &nevt_MIX, "root", "Bzero", "otherwise");
    NevtCount_ri(Nominal_MC_UUBAR_validation_dirname, "UUBAR", &nevt_UUBAR, "root", "Continuum", "otherwise");
    NevtCount_ri(Nominal_MC_DDBAR_validation_dirname, "DDBAR", &nevt_DDBAR, "root", "Continuum", "otherwise");
    NevtCount_ri(Nominal_MC_SSBAR_validation_dirname, "SSBAR", &nevt_SSBAR, "root", "Continuum", "otherwise");
    NevtCount_ri(Nominal_MC_CHARM_validation_dirname, "CHARM", &nevt_CHARM, "root", "Continuum", "otherwise");
    NevtCount_ri(Nominal_MC_SIGNAL_validation_dirname, "SIGNAL", &nevt_SIGNAL, "B2Knunu", "Bplus", "B2Knunu");
    NevtCount_ri(Nominal_MC_SIGNAL_validation_dirname, "SIGNAL", &nevt_SIGNAL, "B2Kstarnunu", "Bplus", "otherwise");
    NevtCount_ri(Nominal_MC_SIGNAL_validation_dirname, "SIGNAL", &nevt_SIGNAL, "B2Xsnunu", "Bplus", "B2Xsnunu");
    NevtCount_ri(Nominal_MC_SIGNAL_validation_dirname, "SIGNAL", &nevt_SIGNAL, "B02K0nunu", "Bzero", "B02K0nunu");
    NevtCount_ri(Nominal_MC_SIGNAL_validation_dirname, "SIGNAL", &nevt_SIGNAL, "B02Kstar0nunu", "Bzero", "otherwise");
    NevtCount_ri(Nominal_MC_SIGNAL_validation_dirname, "SIGNAL", &nevt_SIGNAL, "B02Xsnunu", "Bzero", "B02Xsnunu");

    double NormFactor_CHG = nevt_CHG.NevtwithoutCorrection / nevt_CHG.NevtwithCorrection;
    double NormFactor_MIX = nevt_MIX.NevtwithoutCorrection / nevt_MIX.NevtwithCorrection;
    double NormFactor_UUBAR = nevt_UUBAR.NevtwithoutCorrection / nevt_UUBAR.NevtwithCorrection;
    double NormFactor_DDBAR = nevt_DDBAR.NevtwithoutCorrection / nevt_DDBAR.NevtwithCorrection;
    double NormFactor_SSBAR = nevt_SSBAR.NevtwithoutCorrection / nevt_SSBAR.NevtwithCorrection;
    double NormFactor_CHARM = nevt_CHARM.NevtwithoutCorrection / nevt_CHARM.NevtwithCorrection;
    double NormFactor_SIGNAL = nevt_SIGNAL.NevtwithoutCorrection / nevt_SIGNAL.NevtwithCorrection;

    std::vector<std::string> variable_names;
    std::vector<std::string> branch_names;
    variable_names.push_back("nRemainingTracksInEvent"); branch_names.push_back("Upsilon");
    variable_names.push_back("Btag_chiProb"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_extraInfo_SignalProbability"); branch_names.push_back("Btag");
    variable_names.push_back("thrustAxisCosTheta"); branch_names.push_back("Upsilon");
    variable_names.push_back("missingMomentumOfEvent_theta"); branch_names.push_back("Upsilon");
    variable_names.push_back("missingEnergyOfEventCMS"); branch_names.push_back("Upsilon");
    variable_names.push_back("missingMomentumOfEvent"); branch_names.push_back("Upsilon");
    variable_names.push_back("Btag_deltaE"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_Mbc"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_useCMSFrame_theta"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_R2"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_cosTBTO"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_et"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_mm2"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso00"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso01"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso02"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso03"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso04"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso10"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso12"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso14"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso20"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso22"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso24"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hoo0"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hoo1"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hoo2"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hoo3"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hoo4"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_dr"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_dz"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_useCMSFrame_p"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_useCMSFrame_phi"); branch_names.push_back("Btag");
    variable_names.push_back("roeEextra__bocleanMask__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("extraInfo__boEeclv133__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("extraInfo__boNgammav133__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("extraInfo__boEeclv200__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("extraInfo__boNgammav200__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("Btag_thrustOm"); branch_names.push_back("Btag");
    variable_names.push_back("nParticlesInList__boe__pl__clElectronFBDT__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("nParticlesInList__bomu__pl__clMuonFBDT__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("nParticlesInList__boe__pl__clElectronFBDT_loose__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("nParticlesInList__bomu__pl__clMuonFBDT_loose__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("nParticlesInList__boe__pl__clElectronFBDT_tight__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("nParticlesInList__bomu__pl__clMuonFBDT_tight__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("Bsig_M"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_useCMSFrame_p"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dc_pValue_med"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dc_pValue_std"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dcsimpleveto_chiProb"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dcsimpleveto_dr"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dcsimpleveto_dz"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dcsimpleveto_M"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_D0_pValue_med"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_D0_pValue_std"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_D0simpleveto_chiProb"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_D0simpleveto_dr"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_D0simpleveto_dz"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_D0simpleveto_M"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_thrustBm"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_thrustOm"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_cosTBTO"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_cosTBz"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_et"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_mm2"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hso00"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hso01"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hso02"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hso03"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hso04"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hso10"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hso12"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hso14"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hso20"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hso22"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hso24"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hoo0"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hoo1"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hoo2"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hoo3"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hoo4"); branch_names.push_back("Bsig");
    variable_names.push_back("foxWolframR1"); branch_names.push_back("Upsilon");
    variable_names.push_back("foxWolframR2"); branch_names.push_back("Upsilon");
    variable_names.push_back("foxWolframR3"); branch_names.push_back("Upsilon");
    variable_names.push_back("foxWolframR4"); branch_names.push_back("Upsilon");
    variable_names.push_back("cleoConeThrust0"); branch_names.push_back("Upsilon");
    variable_names.push_back("cleoConeThrust1"); branch_names.push_back("Upsilon");
    variable_names.push_back("cleoConeThrust2"); branch_names.push_back("Upsilon");
    variable_names.push_back("cleoConeThrust3"); branch_names.push_back("Upsilon");
    variable_names.push_back("cleoConeThrust4"); branch_names.push_back("Upsilon");
    variable_names.push_back("cleoConeThrust5"); branch_names.push_back("Upsilon");
    variable_names.push_back("cleoConeThrust6"); branch_names.push_back("Upsilon");
    variable_names.push_back("cleoConeThrust7"); branch_names.push_back("Upsilon");
    variable_names.push_back("cleoConeThrust8"); branch_names.push_back("Upsilon");
    variable_names.push_back("harmonicMomentThrust0"); branch_names.push_back("Upsilon");
    variable_names.push_back("harmonicMomentThrust1"); branch_names.push_back("Upsilon");
    variable_names.push_back("harmonicMomentThrust2"); branch_names.push_back("Upsilon");
    variable_names.push_back("harmonicMomentThrust3"); branch_names.push_back("Upsilon");
    variable_names.push_back("harmonicMomentThrust4"); branch_names.push_back("Upsilon");
    variable_names.push_back("sphericity"); branch_names.push_back("Upsilon");
    variable_names.push_back("aplanarity"); branch_names.push_back("Upsilon");
    variable_names.push_back("Btag_thrustBm"); branch_names.push_back("Btag");
    variable_names.push_back("roePTheta__bocleanMask__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("Btag_cosTBz"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_CleoConeCS_1"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_CleoConeCS_2"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_CleoConeCS_3"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_CleoConeCS_4"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_CleoConeCS_5"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_CleoConeCS_6"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_CleoConeCS_7"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_CleoConeCS_8"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_CleoConeCS_9"); branch_names.push_back("Btag");
    variable_names.push_back("Bsig_daughter_0_extraInfo_mychiProb"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_mydr"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_mydz"); branch_names.push_back("Bsig");
    variable_names.push_back("nROE_KLMClusters"); branch_names.push_back("Upsilon");
    variable_names.push_back("nROE_ECLClusters__bocleanMask__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("nROE_ParticlesInList__bopi0__clmyneutralPion__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("nROE_ParticlesInList__bogamma__clmygamma__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("missingMass2OfEvent"); branch_names.push_back("Upsilon");
    variable_names.push_back("visibleEnergyOfEventCMS"); branch_names.push_back("Upsilon");
    variable_names.push_back("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("chiProb"); branch_names.push_back("Upsilon");
    variable_names.push_back("dr"); branch_names.push_back("Upsilon");
    variable_names.push_back("dz"); branch_names.push_back("Upsilon");

    int Nvar = static_cast<int>(variable_names.size());
    if (Nvar != Nvar_num) exit(1);
    std::vector<double> Nominal_MC_values[Nvar_num];
    std::vector<double> charged_values[Nvar_num];
    std::vector<double> mixed_values[Nvar_num];
    std::vector<double> uubar_values[Nvar_num];
    std::vector<double> ddbar_values[Nvar_num];
    std::vector<double> ssbar_values[Nvar_num];
    std::vector<double> ccbar_values[Nvar_num];
    std::vector<double> taupair_values[Nvar_num];
    std::vector<double> mumu_values[Nvar_num];
    std::vector<double> gg_values[Nvar_num];
    std::vector<double> ee_values[Nvar_num];
    std::vector<double> eeee_values[Nvar_num];
    std::vector<double> eemumu_values[Nvar_num];
    std::vector<double> llXX_values[Nvar_num];
    std::vector<double> hhISR_values[Nvar_num];
    std::vector<double> SIGNAL_values[Nvar_num];
    std::vector<int> Nominal_MC_numbering;

    std::vector<double> Nominal_data_values[Nvar_num];

    std::vector<double> weights;
    std::vector<double> charged_weights;
    std::vector<double> mixed_weights;
    std::vector<double> uubar_weights;
    std::vector<double> ddbar_weights;
    std::vector<double> ssbar_weights;
    std::vector<double> ccbar_weights;
    std::vector<double> taupair_weights;
    std::vector<double> mumu_weights;
    std::vector<double> gg_weights;
    std::vector<double> ee_weights;
    std::vector<double> eeee_weights;
    std::vector<double> eemumu_weights;
    std::vector<double> llXX_weights;
    std::vector<double> hhISR_weights;
    std::vector<double> SIGNAL_weights;

    /*
    LetsFillNominal_ri(Nominal_MC_CHG_validation_dirname, variable_names, branch_names, Nominal_MC_values, &Nominal_MC_numbering, &weights, "CHG");
    LetsFillNominal_ri(Nominal_MC_MIX_validation_dirname, variable_names, branch_names, Nominal_MC_values, &Nominal_MC_numbering, &weights, "MIX");
    LetsFillNominal_ri_correction(Nominal_MC_UUBAR_validation_dirname, variable_names, branch_names, Nominal_MC_values, &Nominal_MC_numbering, &weights, "UUBAR", NormFactor_UUBAR);
    LetsFillNominal_ri_correction(Nominal_MC_DDBAR_validation_dirname, variable_names, branch_names, Nominal_MC_values, &Nominal_MC_numbering, &weights, "DDBAR", NormFactor_DDBAR);
    LetsFillNominal_ri_correction(Nominal_MC_SSBAR_validation_dirname, variable_names, branch_names, Nominal_MC_values, &Nominal_MC_numbering, &weights, "SSBAR", NormFactor_SSBAR);
    LetsFillNominal_ri_correction(Nominal_MC_CHARM_validation_dirname, variable_names, branch_names, Nominal_MC_values, &Nominal_MC_numbering, &weights, "CHARM", NormFactor_CHARM);
    */
    LetsFillNominal_ri(Nominal_MC_CHG_validation_dirname, variable_names, branch_names, Nominal_MC_values, "root", &Nominal_MC_numbering, &weights, "CHG", "Bplus", "otherwise");
    LetsFillNominal_ri(Nominal_MC_MIX_validation_dirname, variable_names, branch_names, Nominal_MC_values, "root", &Nominal_MC_numbering, &weights, "MIX", "Bzero", "otherwise");
    LetsFillNominal_ri(Nominal_MC_UUBAR_validation_dirname, variable_names, branch_names, Nominal_MC_values, "root", &Nominal_MC_numbering, &weights, "UUBAR", "Continuum", "otherwise");
    LetsFillNominal_ri(Nominal_MC_DDBAR_validation_dirname, variable_names, branch_names, Nominal_MC_values, "root", &Nominal_MC_numbering, &weights, "DDBAR", "Continuum", "otherwise");
    LetsFillNominal_ri(Nominal_MC_SSBAR_validation_dirname, variable_names, branch_names, Nominal_MC_values, "root", &Nominal_MC_numbering, &weights, "SSBAR", "Continuum", "otherwise");
    LetsFillNominal_ri(Nominal_MC_CHARM_validation_dirname, variable_names, branch_names, Nominal_MC_values, "root", &Nominal_MC_numbering, &weights, "CHARM", "Continuum", "otherwise");
    LetsFillNominal_ri(Nominal_MC_SIGNAL_validation_dirname, variable_names, branch_names, Nominal_MC_values, "B2Knunu", &Nominal_MC_numbering, &weights, "SIGNAL", "Bplus", "B2Knunu");
    LetsFillNominal_ri(Nominal_MC_SIGNAL_validation_dirname, variable_names, branch_names, Nominal_MC_values, "B2Kstarnunu", &Nominal_MC_numbering, &weights, "SIGNAL", "Bplus", "otherwise");
    LetsFillNominal_ri(Nominal_MC_SIGNAL_validation_dirname, variable_names, branch_names, Nominal_MC_values, "B2Xsnunu", &Nominal_MC_numbering, &weights, "SIGNAL", "Bplus", "B2Xsnunu");
    LetsFillNominal_ri(Nominal_MC_SIGNAL_validation_dirname, variable_names, branch_names, Nominal_MC_values, "B02K0nunu", &Nominal_MC_numbering, &weights, "SIGNAL", "Bzero", "B02K0nunu");
    LetsFillNominal_ri(Nominal_MC_SIGNAL_validation_dirname, variable_names, branch_names, Nominal_MC_values, "B02Kstar0nunu", &Nominal_MC_numbering, &weights, "SIGNAL", "Bzero", "otherwise");
    LetsFillNominal_ri(Nominal_MC_SIGNAL_validation_dirname, variable_names, branch_names, Nominal_MC_values, "B02Xsnunu", &Nominal_MC_numbering, &weights, "SIGNAL", "Bzero", "B02Xsnunu");
    
    LetsFill(Nominal_data_dirname, variable_names, branch_names, Nominal_data_values, "root");

    // sort variables
    for (int k = 0; k < (int)Nominal_MC_numbering.size(); k++) {
        if (Nominal_MC_numbering.at(k) == 0) {
            for (int l = 0; l < (int)variable_names.size(); l++) charged_values[l].push_back(Nominal_MC_values[l].at(k));
            charged_weights.push_back(weights.at(k));
        }
        else if (Nominal_MC_numbering.at(k) == 1) {
            for (int l = 0; l < (int)variable_names.size(); l++) mixed_values[l].push_back(Nominal_MC_values[l].at(k));
            mixed_weights.push_back(weights.at(k));
        }
        else if (Nominal_MC_numbering.at(k) == 2) {
            for (int l = 0; l < (int)variable_names.size(); l++) uubar_values[l].push_back(Nominal_MC_values[l].at(k));
            uubar_weights.push_back(weights.at(k));
        }
        else if (Nominal_MC_numbering.at(k) == 3) {
            for (int l = 0; l < (int)variable_names.size(); l++) ddbar_values[l].push_back(Nominal_MC_values[l].at(k));
            ddbar_weights.push_back(weights.at(k));
        }
        else if (Nominal_MC_numbering.at(k) == 4) {
            for (int l = 0; l < (int)variable_names.size(); l++) ssbar_values[l].push_back(Nominal_MC_values[l].at(k));
            ssbar_weights.push_back(weights.at(k));
        }
        else if (Nominal_MC_numbering.at(k) == 5) {
            for (int l = 0; l < (int)variable_names.size(); l++) ccbar_values[l].push_back(Nominal_MC_values[l].at(k));
            ccbar_weights.push_back(weights.at(k));
        }
        else if (Nominal_MC_numbering.at(k) == 6) {
            for (int l = 0; l < (int)variable_names.size(); l++) taupair_values[l].push_back(Nominal_MC_values[l].at(k));
            taupair_weights.push_back(weights.at(k));
        }
        else if (Nominal_MC_numbering.at(k) == 7) {
            for (int l = 0; l < (int)variable_names.size(); l++) mumu_values[l].push_back(Nominal_MC_values[l].at(k));
            mumu_weights.push_back(weights.at(k));
        }
        else if (Nominal_MC_numbering.at(k) == 8) {
            for (int l = 0; l < (int)variable_names.size(); l++) gg_values[l].push_back(Nominal_MC_values[l].at(k));
            gg_weights.push_back(weights.at(k));
        }
        else if (Nominal_MC_numbering.at(k) == 9) {
            for (int l = 0; l < (int)variable_names.size(); l++) ee_values[l].push_back(Nominal_MC_values[l].at(k));
            ee_weights.push_back(weights.at(k));
        }
        else if (Nominal_MC_numbering.at(k) == 10) {
            for (int l = 0; l < (int)variable_names.size(); l++) eeee_values[l].push_back(Nominal_MC_values[l].at(k));
            eeee_weights.push_back(weights.at(k));
        }
        else if (Nominal_MC_numbering.at(k) == 11) {
            for (int l = 0; l < (int)variable_names.size(); l++) eemumu_values[l].push_back(Nominal_MC_values[l].at(k));
            eemumu_weights.push_back(weights.at(k));
        }
        else if (Nominal_MC_numbering.at(k) == 12) {
            for (int l = 0; l < (int)variable_names.size(); l++) llXX_values[l].push_back(Nominal_MC_values[l].at(k));
            llXX_weights.push_back(weights.at(k));
        }
        else if (Nominal_MC_numbering.at(k) == 13) {
            for (int l = 0; l < (int)variable_names.size(); l++) hhISR_values[l].push_back(Nominal_MC_values[l].at(k));
            hhISR_weights.push_back(weights.at(k));
        }
        else if (Nominal_MC_numbering.at(k) == 14) {
            for (int l = 0; l < (int)variable_names.size(); l++) SIGNAL_values[l].push_back(Nominal_MC_values[l].at(k));
            SIGNAL_weights.push_back(weights.at(k));
        }
        else {
            printf("undefined numbering!\n");
            exit(1);
        }
    }

    THStack* Stack[Nvar_num];
    TH1D* charged_hist[Nvar_num];
    TH1D* mixed_hist[Nvar_num];
    TH1D* uubar_hist[Nvar_num];
    TH1D* ddbar_hist[Nvar_num];
    TH1D* ssbar_hist[Nvar_num];
    TH1D* ccbar_hist[Nvar_num];
    TH1D* taupair_hist[Nvar_num];
    TH1D* mumu_hist[Nvar_num];
    TH1D* gg_hist[Nvar_num];
    TH1D* ee_hist[Nvar_num];
    TH1D* eeee_hist[Nvar_num];
    TH1D* eemumu_hist[Nvar_num];
    TH1D* llXX_hist[Nvar_num];
    TH1D* hhISR_hist[Nvar_num];
    TH1D* SIGNAL_hist[Nvar_num];
    TH1D* stat_error_hist[Nvar_num];
    TH1D* data_hist[Nvar_num];
    TH1D* Ratio_hist[Nvar_num];

    for (int k = 0; k < (int)variable_names.size(); k++) { // malloc TH1D
        std::vector<double> temp_v;
        temp_v.insert(temp_v.end(), charged_values[k].begin(), charged_values[k].end());
        temp_v.insert(temp_v.end(), mixed_values[k].begin(), mixed_values[k].end());
        temp_v.insert(temp_v.end(), uubar_values[k].begin(), uubar_values[k].end());
        temp_v.insert(temp_v.end(), ddbar_values[k].begin(), ddbar_values[k].end());
        temp_v.insert(temp_v.end(), ssbar_values[k].begin(), ssbar_values[k].end());
        temp_v.insert(temp_v.end(), ccbar_values[k].begin(), ccbar_values[k].end());
        temp_v.insert(temp_v.end(), taupair_values[k].begin(), taupair_values[k].end());
        temp_v.insert(temp_v.end(), mumu_values[k].begin(), mumu_values[k].end());
        temp_v.insert(temp_v.end(), gg_values[k].begin(), gg_values[k].end());
        temp_v.insert(temp_v.end(), ee_values[k].begin(), ee_values[k].end());
        temp_v.insert(temp_v.end(), eeee_values[k].begin(), eeee_values[k].end());
        temp_v.insert(temp_v.end(), eemumu_values[k].begin(), eemumu_values[k].end());
        temp_v.insert(temp_v.end(), llXX_values[k].begin(), llXX_values[k].end());
        temp_v.insert(temp_v.end(), hhISR_values[k].begin(), hhISR_values[k].end());
        temp_v.insert(temp_v.end(), SIGNAL_values[k].begin(), SIGNAL_values[k].end());
        temp_v.insert(temp_v.end(), Nominal_data_values[k].begin(), Nominal_data_values[k].end());


        double min = *min_element(temp_v.begin(), temp_v.end());
        double max = *max_element(temp_v.begin(), temp_v.end());
        int bins = 100;

        if (hasEnding(variable_names.at(k), std::string("dr"))) { // exceptions
            max = 0.2;
            min = 0.0;
        }
        else if (hasEnding(variable_names.at(k), std::string("dz"))) {
            max = 0.2;
            min = -0.2;
        }
        else if (hasEnding(variable_names.at(k), std::string("M"))) {
            max = 2.0;
            min = 0.0;
        }
        else if (hasEnding(variable_names.at(k), std::string("chiProb"))) {
            max = 1.0;
            min = 0.0;
        }
        else if (hasEnding(variable_names.at(k), std::string("harmonicMomentThrust1"))) {
            max = 0.6;
            min = -0.6;
        }
        else if (hasEnding(variable_names.at(k), std::string("harmonicMomentThrust2"))) {
            max = 1.0;
            min = 0.0;
        }
        else if (hasEnding(variable_names.at(k), std::string("harmonicMomentThrust3"))) {
            max = 1.0;
            min = -1.0;
        }
        else if (hasEnding(variable_names.at(k), std::string("harmonicMomentThrust4"))) {
            max = 1.0;
            min = -0.5;
        }
        else if (hasEnding(variable_names.at(k), std::string("missingMomentumOfEvent"))) {
            max = 5.0;
        }
        else if (hasEnding(variable_names.at(k), std::string("missingEnergyOfEventCMS"))) {
            min = -1.5;
        }
        else if (hasEnding(variable_names.at(k), std::string("Btag_extraInfo_SignalProbability"))) {
            max = 0;
            min = -3;
            variable_names.at(k) = std::string("log_{10}SignalProbability");
        }
        else if (hasEnding(variable_names.at(k), std::string("Btag_thrustOm"))) {
            min = 0.5;
        }
        else if (hasEnding(variable_names.at(k), std::string("foxWolframR1"))) {
            max = 0.25;
        }
        else if (hasEnding(variable_names.at(k), std::string("foxWolframR3"))) {
            max = 0.4;
        }

        Stack[k] = new THStack(variable_names.at(k).c_str(), (";" + variable_names.at(k) + ";number of candidates").c_str());
        charged_hist[k] = new TH1D("charged", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        mixed_hist[k] = new TH1D("mixed", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        uubar_hist[k] = new TH1D("u#bar{u}", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        ddbar_hist[k] = new TH1D("d#bar{d}", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        ssbar_hist[k] = new TH1D("s#bar{s}", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        ccbar_hist[k] = new TH1D("c#bar{c}", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        taupair_hist[k] = new TH1D("#tau#bar{#tau}", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        mumu_hist[k] = new TH1D("#mu#bar{#mu}", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        gg_hist[k] = new TH1D("gg", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        ee_hist[k] = new TH1D("e#bar{e}", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        eeee_hist[k] = new TH1D("ee#bar{e}#bar{e}", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        eemumu_hist[k] = new TH1D("e#bar{e}#mu#bar{#mu}", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        llXX_hist[k] = new TH1D("\ell#bar{\ell}XX", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        hhISR_hist[k] = new TH1D("hhISR", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        SIGNAL_hist[k] = new TH1D("SIGNAL", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        stat_error_hist[k] = new TH1D("MC stat error", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        data_hist[k] = new TH1D("data", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        Ratio_hist[k] = new TH1D((variable_names.at(k) + "_ratio").c_str(), ";;data/MC", bins, min, max);
    }

    int index = std::find(variable_names.begin(), variable_names.end(), std::string("log_{10}SignalProbability")) - variable_names.begin();
    for (int i = 0; i < (int)charged_values[index].size(); i++) charged_values[index].at(i) = log10l(charged_values[index].at(i));
    for (int i = 0; i < (int)mixed_values[index].size(); i++) mixed_values[index].at(i) = log10l(mixed_values[index].at(i));
    for (int i = 0; i < (int)uubar_values[index].size(); i++) uubar_values[index].at(i) = log10l(uubar_values[index].at(i));
    for (int i = 0; i < (int)ddbar_values[index].size(); i++) ddbar_values[index].at(i) = log10l(ddbar_values[index].at(i));
    for (int i = 0; i < (int)ssbar_values[index].size(); i++) ssbar_values[index].at(i) = log10l(ssbar_values[index].at(i));
    for (int i = 0; i < (int)ccbar_values[index].size(); i++) ccbar_values[index].at(i) = log10l(ccbar_values[index].at(i));
    for (int i = 0; i < (int)taupair_values[index].size(); i++) taupair_values[index].at(i) = log10l(taupair_values[index].at(i));
    for (int i = 0; i < (int)mumu_values[index].size(); i++) mumu_values[index].at(i) = log10l(mumu_values[index].at(i));
    for (int i = 0; i < (int)gg_values[index].size(); i++) gg_values[index].at(i) = log10l(gg_values[index].at(i));
    for (int i = 0; i < (int)ee_values[index].size(); i++) ee_values[index].at(i) = log10l(ee_values[index].at(i));
    for (int i = 0; i < (int)eeee_values[index].size(); i++) eeee_values[index].at(i) = log10l(eeee_values[index].at(i));
    for (int i = 0; i < (int)eemumu_values[index].size(); i++) eemumu_values[index].at(i) = log10l(eemumu_values[index].at(i));
    for (int i = 0; i < (int)llXX_values[index].size(); i++) llXX_values[index].at(i) = log10l(llXX_values[index].at(i));
    for (int i = 0; i < (int)hhISR_values[index].size(); i++) hhISR_values[index].at(i) = log10l(hhISR_values[index].at(i));
    for (int i = 0; i < (int)SIGNAL_values[index].size(); i++) SIGNAL_values[index].at(i) = log10l(SIGNAL_values[index].at(i));
    for (int i = 0; i < (int)Nominal_MC_values[index].size(); i++) Nominal_MC_values[index].at(i) = log10l(Nominal_MC_values[index].at(i));
    for (int i = 0; i < (int)Nominal_data_values[index].size(); i++) Nominal_data_values[index].at(i) = log10l(Nominal_data_values[index].at(i));

    for (int k = 0; k < (int)variable_names.size(); k++) { // fill
        for (int i = 0; i < (int)charged_values[k].size(); i++) charged_hist[k]->Fill(charged_values[k].at(i), charged_weights.at(i));
        for (int i = 0; i < (int)mixed_values[k].size(); i++) mixed_hist[k]->Fill(mixed_values[k].at(i), mixed_weights.at(i));
        for (int i = 0; i < (int)uubar_values[k].size(); i++) uubar_hist[k]->Fill(uubar_values[k].at(i), uubar_weights.at(i));
        for (int i = 0; i < (int)ddbar_values[k].size(); i++) ddbar_hist[k]->Fill(ddbar_values[k].at(i), ddbar_weights.at(i));
        for (int i = 0; i < (int)ssbar_values[k].size(); i++) ssbar_hist[k]->Fill(ssbar_values[k].at(i), ssbar_weights.at(i));
        for (int i = 0; i < (int)ccbar_values[k].size(); i++) ccbar_hist[k]->Fill(ccbar_values[k].at(i), ccbar_weights.at(i));
        for (int i = 0; i < (int)taupair_values[k].size(); i++) taupair_hist[k]->Fill(taupair_values[k].at(i), taupair_weights.at(i));
        for (int i = 0; i < (int)mumu_values[k].size(); i++) mumu_hist[k]->Fill(mumu_values[k].at(i), mumu_weights.at(i));
        for (int i = 0; i < (int)gg_values[k].size(); i++) gg_hist[k]->Fill(gg_values[k].at(i), gg_weights.at(i));
        for (int i = 0; i < (int)ee_values[k].size(); i++) ee_hist[k]->Fill(ee_values[k].at(i), ee_weights.at(i));
        for (int i = 0; i < (int)eeee_values[k].size(); i++) eeee_hist[k]->Fill(eeee_values[k].at(i), eeee_weights.at(i));
        for (int i = 0; i < (int)eemumu_values[k].size(); i++) eemumu_hist[k]->Fill(eemumu_values[k].at(i), eemumu_weights.at(i));
        for (int i = 0; i < (int)llXX_values[k].size(); i++) llXX_hist[k]->Fill(llXX_values[k].at(i), llXX_weights.at(i));
        for (int i = 0; i < (int)hhISR_values[k].size(); i++) hhISR_hist[k]->Fill(hhISR_values[k].at(i), hhISR_weights.at(i));
        for (int i = 0; i < (int)SIGNAL_values[k].size(); i++) SIGNAL_hist[k]->Fill(SIGNAL_values[k].at(i), SIGNAL_weights.at(i));
        for (int i = 0; i < (int)Nominal_MC_values[k].size(); i++) stat_error_hist[k]->Fill(Nominal_MC_values[k].at(i), weights.at(i));
        for (int i = 0; i < (int)Nominal_data_values[k].size(); i++) data_hist[k]->Fill(Nominal_data_values[k].at(i));
    }

    printf("charged: %d\n", (int)charged_values[0].size());
    printf("mixed: %d\n", (int)mixed_values[0].size());
    printf("uubar: %d\n", (int)uubar_values[0].size());
    printf("ddbar: %d\n", (int)ddbar_values[0].size());
    printf("ssbar: %d\n", (int)ssbar_values[0].size());
    printf("ccbar: %d\n", (int)ccbar_values[0].size());
    printf("taupair: %d\n", (int)taupair_values[0].size());
    printf("mumu: %d\n", (int)mumu_values[0].size());
    printf("gg: %d\n", (int)gg_values[0].size());
    printf("ee: %d\n", (int)ee_values[0].size());
    printf("eeee: %d\n", (int)eeee_values[0].size());
    printf("eemuu: %d\n", (int)eemumu_values[0].size());
    printf("llXX: %d\n", (int)llXX_values[0].size());
    printf("hhISR: %d\n", (int)hhISR_values[0].size());
    printf("SIGNAL: %d\n", (int)SIGNAL_values[0].size());
    printf("data: %d\n", (int)Nominal_data_values[0].size());

    for (int k = 0; k < (int)variable_names.size(); k++) { // draw
        Stack[k]->Add(charged_hist[k]);
        Stack[k]->Add(mixed_hist[k]);
        Stack[k]->Add(uubar_hist[k]);
        Stack[k]->Add(ddbar_hist[k]);
        Stack[k]->Add(ssbar_hist[k]);
        Stack[k]->Add(ccbar_hist[k]);
        Stack[k]->Add(SIGNAL_hist[k]);
        //Stack[k]->Add(taupair_hist[k]);
        //Stack[k]->Add(mumu_hist[k]);
        //Stack[k]->Add(gg_hist[k]);
        //Stack[k]->Add(ee_hist[k]);
        //Stack[k]->Add(eeee_hist[k]);
        //Stack[k]->Add(eemumu_hist[k]);
        //Stack[k]->Add(llXX_hist[k]);
        //Stack[k]->Add(hhISR_hist[k]);

        Ratio_hist[k]->SetLineColor(kBlack); Ratio_hist[k]->SetMarkerStyle(21); Ratio_hist[k]->Sumw2(); Ratio_hist[k]->SetStats(0);
        Ratio_hist[k]->Divide(data_hist[k], stat_error_hist[k]);

        TCanvas* c_temp = new TCanvas("c", "", 800, 800); c_temp->cd();

        TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.35, 1.0, 1.0);
        pad1->SetBottomMargin(0.08); pad1->SetLeftMargin(0.15);
        pad1->SetGridx(); pad1->Draw(); pad1->cd();

        gStyle->SetPalette(kPastel);

        Float_t ymax_1 = Stack[k]->GetMaximum();
        Float_t ymax_2 = data_hist[k]->GetMaximum();
        double real_max = 0;
        if (ymax_1 > ymax_2) real_max = ymax_1;
        else real_max = ymax_2;

        Stack[k]->SetMaximum(real_max * 1.1);

        Stack[k]->Draw("pfc Hist");
        stat_error_hist[k]->SetFillColor(12); stat_error_hist[k]->SetLineWidth(0); stat_error_hist[k]->SetFillStyle(3004); stat_error_hist[k]->Draw("e2 SAME");
        data_hist[k]->SetLineWidth(2); data_hist[k]->SetLineColor(kBlack); data_hist[k]->SetMarkerStyle(8); data_hist[k]->Draw("SAME eP");
        TLegend* legend = pad1->BuildLegend(0.9, 0.9, 0.7, 0.7);
        legend->SetFillStyle(0); legend->SetLineWidth(0);
        TPaveText* pt = new TPaveText(0.135, 0.88, 0.5, 1.0, "NDC NB"); pt->SetFillStyle(0); pt->SetLineWidth(0); pt->AddText(("MC scaled to data, Data/MC= " + std::to_string(CAL)).c_str()); pt->Draw();

        c_temp->cd();
        TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1, 0.3); pad2->SetBottomMargin(0.15); pad2->SetLeftMargin(0.15); pad2->SetGridx(); pad2->Draw(); pad2->cd();
        Ratio_hist[k]->SetMinimum(0.5); Ratio_hist[k]->SetMaximum(1.5); Ratio_hist[k]->SetLineWidth(2);
        Ratio_hist[k]->GetYaxis()->SetTitleSize(0.08); Ratio_hist[k]->GetYaxis()->SetTitleOffset(0.5);
        Ratio_hist[k]->GetXaxis()->SetLabelSize(0.08); Ratio_hist[k]->GetYaxis()->SetLabelSize(0.08);
        Ratio_hist[k]->Draw("e0p");
        TLine* line = new TLine(Ratio_hist[k]->GetXaxis()->GetXmin(), 1.0, Ratio_hist[k]->GetXaxis()->GetXmax(), 1.0);
        line->SetLineColor(kRed);
        line->SetLineStyle(1); line->SetLineWidth(3);
        line->Draw();

        c_temp->SetBottomMargin(0.0);
        c_temp->SaveAs((variable_names.at(k) + "_nominal.png").c_str());

        delete c_temp;
    }

    // Print data-MC discrepancy
    double MC_sum = 0;
    for (int i = 0; i < (int)Nominal_MC_values[0].size(); i++) MC_sum = MC_sum + weights.at(i);
    printf("data num: %ld\n", Nominal_data_values[0].size());
    printf("MC num with calibration: %lf\n", MC_sum);
}
