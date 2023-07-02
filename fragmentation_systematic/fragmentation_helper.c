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

// SIGNAL MC sample number befor skimming
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

# define KS0_rel_uncertainty 0.6 // %/cm
# define track_rel_uncertainty 0.69 // %
// # define pi0_correction 0.932
# define pi0_rel_uncertainty ((0.0369 / 0.932) * 100.0) // %
# define Kaon_PID_max_uncertainty 0.1 // not percentage. relative uncertainty
// https://indico.belle2.org/event/6872/contributions/37447/attachments/17127/25504/FEIperformance_B2GM.pdf
# define FEI_cal_Bc 0.679
# define FEI_cal_Bc_uncertainty (0.017/FEI_cal_Bc) // not percentage. relative uncertainty
# define FEI_cal_B0 0.713
# define FEI_cal_B0_uncertainty (0.019/FEI_cal_B0) // not percentage. relative uncertainty

# define Nvar_num 126

# define CAL 1.1728
# define CAL_qq 1.0
# define Stream 0.25

# define N_PID_syst 73
# define N_fakeE_syst 37
# define N_fakeMU_syst 49
double PID_correction[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
double PID_correction_stat_uncer[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
double PID_correction_sys_uncer[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
double PID_correction_uncer[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss

double PID_fakeE_correction[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
double PID_fakeE_uncer[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+

double PID_fakeMU_correction[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+
double PID_fakeMU_uncer[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

# define N_pi0_syst 8
const double pi0_correction[N_pi0_syst] = {
    0.917, 0.965, 0.988, 1.013, 1.042, 1.044, 1.011, 1.0
};
const double pi0_stat_uncer[N_pi0_syst] = {
    0.004, 0.004, 0.004, 0.005, 0.004, 0.005, 0.005, 0.0
};
const double pi0_sys_uncer1[N_pi0_syst] = {
    0.049, 0.036, 0.079, 0.058, 0.045, 0.041, 0.040, 0.0
};
const double pi0_sys_uncer2[N_pi0_syst] = {
    0.0, 0.0, 0.0, 0.0, 0.039, 0.051, 0.030, 0.0
};

void ReadPIDFile() {
    const char* KID_true_file = "KaonEff.csv";
    const char* KID_mis_file = "Kaonmis.csv";
    const char* PID_true_file = "PionEff.csv";
    const char* PID_mis_file = "Pionmis.csv";

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
        fscanf(fp_KID_true, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,kaonID\n", &temp_p_min, &temp_p_max, &temp_cosTheta_min, &temp_cosTheta_max, &PID_correction[0][i], &PID_correction_stat_uncer[0][i], &temp_data_MC_uncertainty_stat_dn, &PID_correction_sys_uncer[0][i], &temp_data_MC_uncertainty_sys_dn, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_threshold);
        fscanf(fp_KID_mis, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,kaonID\n", &temp_p_min, &temp_p_max, &temp_cosTheta_min, &temp_cosTheta_max, &PID_correction[1][i], &PID_correction_stat_uncer[1][i], &temp_data_MC_uncertainty_stat_dn, &PID_correction_sys_uncer[1][i], &temp_data_MC_uncertainty_sys_dn, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_threshold);
        fscanf(fp_PID_true, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,pionID\n", &temp_p_min, &temp_p_max, &temp_cosTheta_min, &temp_cosTheta_max, &PID_correction[2][i], &PID_correction_stat_uncer[2][i], &temp_data_MC_uncertainty_stat_dn, &PID_correction_sys_uncer[2][i], &temp_data_MC_uncertainty_sys_dn, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_threshold);
        fscanf(fp_PID_mis, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,pionID\n", &temp_p_min, &temp_p_max, &temp_cosTheta_min, &temp_cosTheta_max, &PID_correction[3][i], &PID_correction_stat_uncer[3][i], &temp_data_MC_uncertainty_stat_dn, &PID_correction_sys_uncer[3][i], &temp_data_MC_uncertainty_sys_dn, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_threshold);

        if ((std::abs(PID_correction[0][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer[0][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer[0][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction[0][i]) > 10000.0)) {
            PID_correction[0][i] = 1.0;
            PID_correction_stat_uncer[0][i] = 0.0;
            PID_correction_sys_uncer[0][i] = 0.0;
        }
        if ((std::abs(PID_correction[1][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer[1][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer[1][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction[1][i]) > 10000.0)) {
            PID_correction[1][i] = 1.0;
            PID_correction_stat_uncer[1][i] = 0.0;
            PID_correction_sys_uncer[1][i] = 0.0;
        }
        if ((std::abs(PID_correction[2][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer[2][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer[2][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction[2][i]) > 10000.0)) {
            PID_correction[2][i] = 1.0;
            PID_correction_stat_uncer[2][i] = 0.0;
            PID_correction_sys_uncer[2][i] = 0.0;
        }
        if ((std::abs(PID_correction[3][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer[3][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer[3][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction[3][i]) > 10000.0)) {
            PID_correction[3][i] = 1.0;
            PID_correction_stat_uncer[3][i] = 0.0;
            PID_correction_sys_uncer[3][i] = 0.0;
        }

        PID_correction_uncer[0][i] = std::sqrt(PID_correction_stat_uncer[0][i] * PID_correction_stat_uncer[0][i] + PID_correction_sys_uncer[0][i] * PID_correction_sys_uncer[0][i]);
        PID_correction_uncer[1][i] = std::sqrt(PID_correction_stat_uncer[1][i] * PID_correction_stat_uncer[1][i] + PID_correction_sys_uncer[1][i] * PID_correction_sys_uncer[1][i]);
        PID_correction_uncer[2][i] = std::sqrt(PID_correction_stat_uncer[2][i] * PID_correction_stat_uncer[2][i] + PID_correction_sys_uncer[2][i] * PID_correction_sys_uncer[2][i]);
        PID_correction_uncer[3][i] = std::sqrt(PID_correction_stat_uncer[3][i] * PID_correction_stat_uncer[3][i] + PID_correction_sys_uncer[3][i] * PID_correction_sys_uncer[3][i]);
    }

    PID_correction[0][N_PID_syst - 1] = 1.0;
    PID_correction_stat_uncer[0][N_PID_syst - 1] = 0.0;
    PID_correction_sys_uncer[0][N_PID_syst - 1] = 0.0;
    PID_correction_uncer[0][N_PID_syst - 1] = 0.0;

    PID_correction[1][N_PID_syst - 1] = 1.0;
    PID_correction_stat_uncer[1][N_PID_syst - 1] = 0.0;
    PID_correction_sys_uncer[1][N_PID_syst - 1] = 0.0;
    PID_correction_uncer[1][N_PID_syst - 1] = 0.0;

    PID_correction[2][N_PID_syst - 1] = 1.0;
    PID_correction_stat_uncer[2][N_PID_syst - 1] = 0.0;
    PID_correction_sys_uncer[2][N_PID_syst - 1] = 0.0;
    PID_correction_uncer[2][N_PID_syst - 1] = 0.0;

    PID_correction[3][N_PID_syst - 1] = 1.0;
    PID_correction_stat_uncer[3][N_PID_syst - 1] = 0.0;
    PID_correction_sys_uncer[3][N_PID_syst - 1] = 0.0;
    PID_correction_uncer[3][N_PID_syst - 1] = 0.0;

    fclose(fp_KID_true);
    fclose(fp_KID_mis);
    fclose(fp_PID_true);
    fclose(fp_PID_mis);
}

void ReadFakePIDFile() {
    // initialization
    for (int i = 0; i < N_fakeE_syst; i++) {
        PID_fakeE_correction[0][i] = 1.0; //  K-, K+, pi-, pi+
        PID_fakeE_correction[1][i] = 1.0;
        PID_fakeE_correction[2][i] = 1.0;
        PID_fakeE_correction[3][i] = 1.0;

        PID_fakeE_uncer[0][i] = 0.0;
        PID_fakeE_uncer[1][i] = 0.0;
        PID_fakeE_uncer[2][i] = 0.0;
        PID_fakeE_uncer[3][i] = 0.0;
    }

    for (int i = 0; i < N_fakeMU_syst; i++) {
        PID_fakeMU_correction[0][i] = 1.0;
        PID_fakeMU_correction[1][i] = 1.0;
        PID_fakeMU_correction[2][i] = 1.0;
        PID_fakeMU_correction[3][i] = 1.0;

        PID_fakeMU_uncer[0][i] = 0.0;
        PID_fakeMU_uncer[1][i] = 0.0;
        PID_fakeMU_uncer[2][i] = 0.0;
        PID_fakeMU_uncer[3][i] = 0.0;
    }

    const char* K_fromE_file = "kaonID_efficiency_electron.csv";
    const char* K_fromMU_file = "kaonID_efficiency_muon.csv";
    const char* pi_fromE_file = "pionID_efficiency_electron.csv";
    const char* pi_fromMU_file = "pionID_efficiency_muon.csv";

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

            PID_fakeE_correction[1][bin] = temp_data_MC_ratio;
            PID_fakeE_uncer[1][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeE_correction[1][bin]) < MyEPSILON) {
                PID_fakeE_correction[1][bin] = 1.0;
                PID_fakeE_uncer[1][bin] = 0.0;
            }

        }
        else if (temp_charge == '-') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeE_correction[0][bin] = temp_data_MC_ratio;
            PID_fakeE_uncer[0][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeE_correction[0][bin]) < MyEPSILON) {
                PID_fakeE_correction[0][bin] = 1.0;
                PID_fakeE_uncer[0][bin] = 0.0;
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

            PID_fakeE_correction[3][bin] = temp_data_MC_ratio;
            PID_fakeE_uncer[3][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeE_correction[3][bin]) < MyEPSILON) {
                PID_fakeE_correction[3][bin] = 1.0;
                PID_fakeE_uncer[3][bin] = 0.0;
            }

        }
        else if (temp_charge == '-') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeE_correction[2][bin] = temp_data_MC_ratio;
            PID_fakeE_uncer[2][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeE_correction[2][bin]) < MyEPSILON) {
                PID_fakeE_correction[2][bin] = 1.0;
                PID_fakeE_uncer[2][bin] = 0.0;
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

            PID_fakeMU_correction[1][bin] = temp_data_MC_ratio;
            PID_fakeMU_uncer[1][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeMU_correction[1][bin]) < MyEPSILON) {
                PID_fakeMU_correction[1][bin] = 1.0;
                PID_fakeMU_uncer[1][bin] = 0.0;
            }

        }
        else if (temp_charge == '-') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeMU_correction[0][bin] = temp_data_MC_ratio;
            PID_fakeMU_uncer[0][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeMU_correction[0][bin]) < MyEPSILON) {
                PID_fakeMU_correction[0][bin] = 1.0;
                PID_fakeMU_uncer[0][bin] = 0.0;
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

            PID_fakeMU_correction[3][bin] = temp_data_MC_ratio;
            PID_fakeMU_uncer[3][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeMU_correction[3][bin]) < MyEPSILON) {
                PID_fakeMU_correction[3][bin] = 1.0;
                PID_fakeMU_uncer[3][bin] = 0.0;
            }

        }
        else if (temp_charge == '-') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeMU_correction[2][bin] = temp_data_MC_ratio;
            PID_fakeMU_uncer[2][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeMU_correction[2][bin]) < MyEPSILON) {
                PID_fakeMU_correction[2][bin] = 1.0;
                PID_fakeMU_uncer[2][bin] = 0.0;
            }

        }
        else {
            printf("[ERROR] unknown charge!\n");
            exit(1);
        }
    }
    fclose(fp_pi_fromMU);
}

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
enum SimpleDecayMode { // reco level
    B2Kpi = 0,
    B2Kpipi,
    reco_other,
    MAX_NUM_SIMPLE_DECAYMODE
};
enum SimpleDecayModeMC { // MC level
    Xs2Kpi_MC = 0, // Xs -> K pi nonresonant
    Xs2Kpipi_MC, // Xs -> K pi pi nonresonant
    MC_other, // other non-resonant Xs decay
    MAX_NUM_SIMPLE_DECAYMODE_MC
};
enum DecayModeReduced { // reco level
    rB2KcPi0 = 0,
    rB2Ks0Pic,
    rB2KcPicPic,
    rB2Ks0PicPi0,
    rB2KcPicPicPi0,
    rB2Ks0PicPicPic,
    rB2KcPicPicPicPic,
    rB2Ks0PicPicPicPi0,
    rB2KcPi0Pi0,
    rB2Ks0PicPi0Pi0,
    rB2KcPicPicPi0Pi0,
    rB2KcKcKc,
    rB2KcKcKs0Pic,
    rB2KcKcKcPi0,
    rB02KcPic,
    rB02Ks0Pi0,
    rB02KcPicPi0,
    rB02Ks0PicPic,
    rB02KcPicPicPic,
    rB02Ks0PicPicPi0,
    rB02KcPicPicPicPi0,
    rB02Ks0PicPicPicPic,
    rB02Ks0Pi0Pi0,
    rB02KcPicPi0Pi0,
    rB02Ks0PicPicPi0Pi0,
    rB02KcKcKs0,
    rB02KcKcKcPic,
    rB02KcKcKs0Pi0,
    rMAX_NUM_DECAYMODE
};
enum DecayModeReducedMC { // MC level
    rXsu2KcPi0_MC = 0, // non-res
    rXsu2K0Pic_MC,
    rXsu2KcPicPic_MC,
    rXsu2K0PicPi0_MC,
    rXsu2KcPicPicPi0_MC,
    rXsu2K0PicPicPic_MC,
    rXsu2KcPicPicPicPic_MC,
    rXsu2K0PicPicPicPi0_MC,
    rXsu2KcPi0Pi0_MC,
    rXsu2K0PicPi0Pi0_MC,
    rXsu2KcPicPicPi0Pi0_MC,
    rXsu2KcKcKc_MC,
    rXsu2KcKcK0Pic_MC,
    rXsu2KcKcKcPi0_MC,
    rXsd2KcPic_MC,
    rXsd2K0Pi0_MC,
    rXsd2KcPicPi0_MC,
    rXsd2K0PicPic_MC,
    rXsd2KcPicPicPic_MC,
    rXsd2K0PicPicPi0_MC,
    rXsd2KcPicPicPicPi0_MC,
    rXsd2K0PicPicPicPic_MC,
    rXsd2K0Pi0Pi0_MC,
    rXsd2KcPicPi0Pi0_MC,
    rXsd2K0PicPicPi0Pi0_MC,
    rXsd2KcKcK0_MC,
    rXsd2KcKcKcPic_MC,
    rXsd2KcKcK0Pi0_MC,
    rother,
    rMAX_NUM_DECAYMODE_MC
};

const double N_initial_MC[DecayModeMC::MAX_NUM_DECAYMODE_MC] = {
    23733.733263,
    11217.329797,
    22435.028056,
    4230.433342,
    8386.683031,
    34488.397795,
    38705.566510,
    21384.911778,
    8143.375388,
    5025.082736,
    6635.015387,
    7975.970886,
    5662.030457,
    5989.224581,
    1379.152711,
    747.977559,
    393.762899,
    19391.715296,
    18894.214376,
    9498.627921,
    10499.285975,
    5240.605341,
    34712.309001,
    46080.930494,
    17038.115418,
    20442.541511,
    6439.845695,
    3004.187464,
    4427.954926,
    5339.097244,
    4745.440684,
    1097.305539,
    503.403056,
    366.670228,
    104680.699924
};
const double N_BKG_initial = 1704693583.076672;

const double rN_initial_MC[DecayModeReducedMC::rMAX_NUM_DECAYMODE_MC] = {
    4230.433342,
    8386.683031,
    34488.397795,
    38705.566510,
    21384.911778,
    8143.375388,
    5025.082736,
    6635.015387,
    7975.970886,
    5662.030457,
    5989.224581,
    1379.152711,
    747.977559,
    393.762899,
    10499.285975,
    5240.605341,
    34712.309001,
    46080.930494,
    17038.115418,
    20442.541511,
    6439.845695,
    3004.187464,
    4427.954926,
    5339.097244,
    4745.440684,
    1097.305539,
    503.403056,
    366.670228,
    104680.699924
};

bool TrueIfDecayModeMatch(double Upsilon_decayID, double Bsig_decayID, DecayMode decaymode) {
    switch (decaymode) {
    case DecayMode::B2Kc:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > -0.5 && Bsig_decayID < 0.5) return true;
        return false;
        break;
    case DecayMode::B2KcPi0:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 0.5 && Bsig_decayID < 1.5) return true;
        return false;
        break;
    case DecayMode::B2Ks0Pic:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 1.5 && Bsig_decayID < 2.5) return true;
        return false;
        break;
    case DecayMode::B2KcPicPic:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 2.5 && Bsig_decayID < 3.5) return true;
        return false;
        break;
    case DecayMode::B2Ks0PicPi0:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 3.5 && Bsig_decayID < 4.5) return true;
        return false;
        break;
    case DecayMode::B2KcPicPicPi0:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 4.5 && Bsig_decayID < 5.5) return true;
        return false;
        break;
    case DecayMode::B2Ks0PicPicPic:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 5.5 && Bsig_decayID < 6.5) return true;
        return false;
        break;
    case DecayMode::B2KcPicPicPicPic:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 6.5 && Bsig_decayID < 7.5) return true;
        return false;
        break;
    case DecayMode::B2Ks0PicPicPicPi0:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 7.5 && Bsig_decayID < 8.5) return true;
        return false;
        break;
    case DecayMode::B2KcPi0Pi0:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 8.5 && Bsig_decayID < 9.5) return true;
        return false;
        break;
    case DecayMode::B2Ks0PicPi0Pi0:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 9.5 && Bsig_decayID < 10.5) return true;
        return false;
        break;
    case DecayMode::B2KcPicPicPi0Pi0:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 10.5 && Bsig_decayID < 11.5) return true;
        return false;
        break;
    case DecayMode::B2KcKcKc:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 11.5 && Bsig_decayID < 12.5) return true;
        return false;
        break;
    case DecayMode::B2KcKcKs0Pic:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 12.5 && Bsig_decayID < 13.5) return true;
        return false;
        break;
    case DecayMode::B2KcKcKcPi0:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 13.5 && Bsig_decayID < 14.5) return true;
        return false;
        break;
    case DecayMode::B02Ks0:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > -0.5 && Bsig_decayID < 0.5) return true;
        return false;
        break;
    case DecayMode::B02KcPic:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 0.5 && Bsig_decayID < 1.5) return true;
        return false;
        break;
    case DecayMode::B02Ks0Pi0:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 1.5 && Bsig_decayID < 2.5) return true;
        return false;
        break;
    case DecayMode::B02KcPicPi0:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 2.5 && Bsig_decayID < 3.5) return true;
        return false;
        break;
    case DecayMode::B02Ks0PicPic:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 3.5 && Bsig_decayID < 4.5) return true;
        return false;
        break;
    case DecayMode::B02KcPicPicPic:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 4.5 && Bsig_decayID < 5.5) return true;
        return false;
        break;
    case DecayMode::B02Ks0PicPicPi0:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 5.5 && Bsig_decayID < 6.5) return true;
        return false;
        break;
    case DecayMode::B02KcPicPicPicPi0:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 6.5 && Bsig_decayID < 7.5) return true;
        return false;
        break;
    case DecayMode::B02Ks0PicPicPicPic:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 7.5 && Bsig_decayID < 8.5) return true;
        return false;
        break;
    case DecayMode::B02Ks0Pi0Pi0:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 8.5 && Bsig_decayID < 9.5) return true;
        return false;
        break;
    case DecayMode::B02KcPicPi0Pi0:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 9.5 && Bsig_decayID < 10.5) return true;
        return false;
        break;
    case DecayMode::B02Ks0PicPicPi0Pi0:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 10.5 && Bsig_decayID < 11.5) return true;
        return false;
        break;
    case DecayMode::B02KcKcKs0:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 11.5 && Bsig_decayID < 12.5) return true;
        return false;
        break;
    case DecayMode::B02KcKcKcPic:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 12.5 && Bsig_decayID < 13.5) return true;
        return false;
        break;
    case DecayMode::B02KcKcKs0Pi0:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 13.5 && Bsig_decayID < 14.5) return true;
        return false;
        break;
    default:
        printf("ERROR! 036\n");
        exit(1);
        break;
    }

    printf("ERROR! 037\n");
    exit(1);
    return false;
}

bool TrueIfDecayModeMatch_MC(int Decay[N_decay], DecayModeMC decaymodeMC) {

    switch (decaymodeMC) {
    case DecayModeMC::Xsu2Kc_MC:
        if (Decay[0] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2Kcstar2KcPi0_MC:
        if (Decay[2] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2Kcstar2K0Pic_MC:
        if (Decay[1] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2KcPi0_MC:
        if (Decay[5] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2K0Pic_MC:
        if (Decay[6] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2KcPicPic_MC:
        if (Decay[7] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2K0PicPi0_MC:
        if (Decay[8] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2KcPicPicPi0_MC:
        if (Decay[9] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2K0PicPicPic_MC:
        if (Decay[10] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2KcPicPicPicPic_MC:
        if (Decay[11] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2K0PicPicPicPi0_MC:
        if (Decay[12] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2KcPi0Pi0_MC:
        if (Decay[13] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2K0PicPi0Pi0_MC:
        if (Decay[14] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2KcPicPicPi0Pi0_MC:
        if (Decay[15] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2KcKcKc_MC:
        if (Decay[16] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2KcKcK0Pic_MC:
        if (Decay[17] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2KcKcKcPi0_MC:
        if (Decay[18] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2K0_MC:
        if (Decay[19] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2K0star2KcPic_MC:
        if (Decay[20] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2K0star2K0Pi0_MC:
        if (Decay[21] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2KcPic_MC:
        if (Decay[24] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2K0Pi0_MC:
        if (Decay[25] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2KcPicPi0_MC:
        if (Decay[26] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2K0PicPic_MC:
        if (Decay[27] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2KcPicPicPic_MC:
        if (Decay[28] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2K0PicPicPi0_MC:
        if (Decay[29] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2KcPicPicPicPi0_MC:
        if (Decay[30] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2K0PicPicPicPic_MC:
        if (Decay[31] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2K0Pi0Pi0_MC:
        if (Decay[32] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2KcPicPi0Pi0_MC:
        if (Decay[33] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2K0PicPicPi0Pi0_MC:
        if (Decay[34] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2KcKcK0_MC:
        if (Decay[35] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2KcKcKcPic_MC:
        if (Decay[36] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2KcKcK0Pi0_MC:
        if (Decay[37] > 0) return true;
        return false;
        break;
    case DecayModeMC::other:
        return true;
        break;
    default:
        printf("ERROR! Input value of TrueIfDecayModeMatch_MC is not appropriate\n");
        exit(1);
        break;
    }

    printf("ERROR! Input value of TrueIfDecayModeMatch_MC is not appropriate\n");
    exit(1);
    return false;
}

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

void LetsFillJpsi_ri(const char* dirname, double Confusion[DecayMode::MAX_NUM_DECAYMODE][DecayModeMC::MAX_NUM_DECAYMODE_MC], std::string SampleName) {
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

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    double FEI_calibration_factor = -1;

    double Mxs = -1;

    int nBp = -1;
    int nB0 = -1;

    int Decay[N_decay] = { 0 }; // MC level info

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");
        TTree* tree_Xs = (TTree*)input_file->Get("Xs");

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID); // charged: 0, mixed: 1
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKtruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[0][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKmisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[1][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npitruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[2][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npimisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[3][i_PID]);
        }
        for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npi0bin" + std::to_string(i_pi0)).c_str(), &temp_N_bin_pi0[i_pi0]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clPrimaryMC__bc", &nBp);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clPrimaryMC__bc", &nB0);
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

        tree_Bsig->SetBranchAddress("Bsig_M", &Mxs);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            tree_Xs->GetEntry(j);

            // fix a double counting
            if (Decay[1] > 0) Decay[6] = 0;
            if (Decay[2] > 0) Decay[5] = 0;

            if (Mxs < 1.1) continue;

            double weight_ri = 0.0;
            if (SampleName == "SIGNAL") {
                FEI_calibration_factor = CAL_qq;
                if (nBp > 0) {
                    FEI_calibration_factor = FEI_cal_Bc;
                    weight_ri = ((N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab)); // total 2.8/ab
                }
                else if (nB0 > 0) {
                    FEI_calibration_factor = FEI_cal_B0;
                    weight_ri = ((N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab)); // total 2.8/ab
                }
                else {
                    printf("ERROR 255");
                    exit(1);
                }
            }
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
                Correction_KID = Correction_KID * std::pow(PID_correction[0][i_PID], temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(PID_correction[1][i_PID], temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(PID_correction[2][i_PID], temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(PID_correction[3][i_PID], temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(pi0_correction[i_pi0], temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(PID_fakeE_correction[0][i_fake], temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_fake = Correction_fake * std::pow(PID_fakeE_correction[1][i_fake], temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(PID_fakeE_correction[2][i_fake], temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(PID_fakeE_correction[3][i_fake], temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(PID_fakeMU_correction[0][i_fake], temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_fake = Correction_fake * std::pow(PID_fakeMU_correction[1][i_fake], temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(PID_fakeMU_correction[2][i_fake], temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(PID_fakeMU_correction[3][i_fake], temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }

            int decaymodeid = -1;
            int decaymodeid_MC = -1;

            for (int i = 0; i < DecayMode::MAX_NUM_DECAYMODE; i++) { // find reco decay mode
                if (TrueIfDecayModeMatch(Upsilon_ID, Bsig_ID, static_cast<DecayMode>(i))) {
                    decaymodeid = i;
                    break;
                }
            }
            if (decaymodeid == DecayMode::MAX_NUM_DECAYMODE) {
                printf("ERROR! Reco decay id cannot be found\n");
                exit(1);
            }

            for (int i = 0; i < DecayModeMC::MAX_NUM_DECAYMODE_MC; i++) { // find MC decay mode
                if (TrueIfDecayModeMatch_MC(Decay, static_cast<DecayModeMC>(i))) {
                    decaymodeid_MC = i;
                    break;
                }
            }
            if (decaymodeid_MC == DecayModeMC::MAX_NUM_DECAYMODE_MC) {
                printf("ERROR! MC decay id cannot be found\n");
                exit(1);
            }

            Confusion[decaymodeid][decaymodeid_MC] = Confusion[decaymodeid][decaymodeid_MC] + FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake;


        }
        input_file->Close();

    }

}

void LetsFillJpsi_ri(const char* dirname, double OneDConfusion[DecayMode::MAX_NUM_DECAYMODE], std::string SampleName) {
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

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    double FEI_calibration_factor = -1;

    double Mxs = -1;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID); // charged: 0, mixed: 1
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
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

        tree_Bsig->SetBranchAddress("Bsig_M", &Mxs);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            if (Mxs < 1.1) continue;

            double weight_ri = 0.0;
            if (SampleName == "CHG") {
                FEI_calibration_factor = FEI_cal_Bc;
                weight_ri = ((N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab)); // total 2.8/ab for BB
            }
            else if (SampleName == "MIX") {
                FEI_calibration_factor = FEI_cal_B0;
                weight_ri = ((N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab)); // total 2.8/ab for BB
            }
            else if (SampleName == "UUBAR") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = ((0.364436 - 0.002763) / 1.0); // total 1.0/ab for qq
            }
            else if (SampleName == "DDBAR") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = ((0.364436 - 0.002763) / 1.0); // total 1.0/ab for qq
            }
            else if (SampleName == "SSBAR") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = ((0.364436 - 0.002763) / 1.0); // total 1.0/ab for qq
            }
            else if (SampleName == "CHARM") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = ((0.364436 - 0.002763) / 1.0); // total 1.0/ab for qq
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
                Correction_KID = Correction_KID * std::pow(PID_correction[0][i_PID], temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(PID_correction[1][i_PID], temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(PID_correction[2][i_PID], temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(PID_correction[3][i_PID], temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(pi0_correction[i_pi0], temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(PID_fakeE_correction[0][i_fake], temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_fake = Correction_fake * std::pow(PID_fakeE_correction[1][i_fake], temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(PID_fakeE_correction[2][i_fake], temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(PID_fakeE_correction[3][i_fake], temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(PID_fakeMU_correction[0][i_fake], temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_fake = Correction_fake * std::pow(PID_fakeMU_correction[1][i_fake], temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(PID_fakeMU_correction[2][i_fake], temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(PID_fakeMU_correction[3][i_fake], temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }

            int decaymodeid = -1;

            for (int i = 0; i < DecayMode::MAX_NUM_DECAYMODE; i++) { // find reco decay mode
                if (TrueIfDecayModeMatch(Upsilon_ID, Bsig_ID, static_cast<DecayMode>(i))) {
                    decaymodeid = i;
                    break;
                }
            }
            if (decaymodeid == DecayMode::MAX_NUM_DECAYMODE) {
                printf("ERROR! Reco decay id cannot be found\n");
                exit(1);
            }

            OneDConfusion[decaymodeid] = OneDConfusion[decaymodeid] + FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake;


        }
        input_file->Close();

    }

}

void LetsFill(const char* dirname, double OneDConfusion[DecayMode::MAX_NUM_DECAYMODE]) {

    std::vector<string> names;
    load_files(dirname, &names);

    double Upsilon_decayID;
    double Bsig_decayID;

    double Mxs = -1;

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_decayID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_decayID);

        tree_Bsig->SetBranchAddress("Bsig_M", &Mxs);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            if (Mxs < 1.1) continue;

            int decaymodeid = -1;
            for (int i = 0; i < DecayMode::MAX_NUM_DECAYMODE; i++) { // find reco decay mode
                if (TrueIfDecayModeMatch(Upsilon_decayID, Bsig_decayID, static_cast<DecayMode>(i))) {
                    decaymodeid = i;
                    break;
                }
            }
            if (decaymodeid == DecayMode::MAX_NUM_DECAYMODE) {
                printf("ERROR! Reco decay id cannot be found\n");
                exit(1);
            }

            OneDConfusion[decaymodeid] = OneDConfusion[decaymodeid] + 1.0;

        }
        input_file->Close();

    }

}

void LetsFill(const char* dirname, double OneDConfusion[DecayMode::MAX_NUM_DECAYMODE], const char* included_string) {

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Upsilon_decayID;
    double Bsig_decayID;

    double Mxs = -1;

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_decayID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_decayID);

        tree_Bsig->SetBranchAddress("Bsig_M", &Mxs);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            if (Mxs < 1.1) continue;

            int decaymodeid = -1;
            for (int i = 0; i < DecayMode::MAX_NUM_DECAYMODE; i++) { // find reco decay mode
                if (TrueIfDecayModeMatch(Upsilon_decayID, Bsig_decayID, static_cast<DecayMode>(i))) {
                    decaymodeid = i;
                    break;
                }
            }
            if (decaymodeid == DecayMode::MAX_NUM_DECAYMODE) {
                printf("ERROR! Reco decay id cannot be found\n");
                exit(1);
            }

            OneDConfusion[decaymodeid] = OneDConfusion[decaymodeid] + 1.0;

        }
        input_file->Close();

    }

}

void DecayMatrixToSimpleDecayMatrix(double Confusion[DecayMode::MAX_NUM_DECAYMODE][DecayModeMC::MAX_NUM_DECAYMODE_MC], double SimpleConfusion[DecayMode::MAX_NUM_DECAYMODE][SimpleDecayModeMC::MAX_NUM_SIMPLE_DECAYMODE_MC]) {
    for (int i = 0; i < DecayMode::MAX_NUM_DECAYMODE; i++) { 

        SimpleConfusion[i][SimpleDecayModeMC::Xs2Kpi_MC]
            = Confusion[i][DecayModeMC::Xsu2KcPi0_MC]
            + Confusion[i][DecayModeMC::Xsu2K0Pic_MC]
            + Confusion[i][DecayModeMC::Xsd2KcPic_MC]
            + Confusion[i][DecayModeMC::Xsd2K0Pi0_MC];

        SimpleConfusion[i][SimpleDecayModeMC::Xs2Kpipi_MC]
            = Confusion[i][DecayModeMC::Xsu2KcPicPic_MC]
            + Confusion[i][DecayModeMC::Xsu2K0PicPi0_MC]
            + Confusion[i][DecayModeMC::Xsu2KcPi0Pi0_MC]
            + Confusion[i][DecayModeMC::Xsd2KcPicPi0_MC]
            + Confusion[i][DecayModeMC::Xsd2K0PicPic_MC]
            + Confusion[i][DecayModeMC::Xsd2K0Pi0Pi0_MC];

        SimpleConfusion[i][SimpleDecayModeMC::MC_other]
            = Confusion[i][DecayModeMC::Xsu2KcPicPicPi0_MC]
            + Confusion[i][DecayModeMC::Xsu2K0PicPicPic_MC]
            + Confusion[i][DecayModeMC::Xsu2KcPicPicPicPic_MC]
            + Confusion[i][DecayModeMC::Xsu2K0PicPicPicPi0_MC]
            + Confusion[i][DecayModeMC::Xsu2K0PicPi0Pi0_MC]
            + Confusion[i][DecayModeMC::Xsu2KcPicPicPi0Pi0_MC]
            + Confusion[i][DecayModeMC::Xsu2KcKcKc_MC]
            + Confusion[i][DecayModeMC::Xsu2KcKcK0Pic_MC]
            + Confusion[i][DecayModeMC::Xsu2KcKcKcPi0_MC]
            + Confusion[i][DecayModeMC::Xsd2KcPicPicPic_MC]
            + Confusion[i][DecayModeMC::Xsd2K0PicPicPi0_MC]
            + Confusion[i][DecayModeMC::Xsd2KcPicPicPicPi0_MC]
            + Confusion[i][DecayModeMC::Xsd2K0PicPicPicPic_MC]
            + Confusion[i][DecayModeMC::Xsd2KcPicPi0Pi0_MC]
            + Confusion[i][DecayModeMC::Xsd2K0PicPicPi0Pi0_MC]
            + Confusion[i][DecayModeMC::Xsd2KcKcK0_MC]
            + Confusion[i][DecayModeMC::Xsd2KcKcKcPic_MC]
            + Confusion[i][DecayModeMC::Xsd2KcKcK0Pi0_MC]
            + Confusion[i][DecayModeMC::other];

    }
}

void DecayMatrixToReducedDecayMatrix(double Confusion[DecayMode::MAX_NUM_DECAYMODE][DecayModeMC::MAX_NUM_DECAYMODE_MC], double rConfusion[DecayModeReduced::rMAX_NUM_DECAYMODE][DecayModeReducedMC::rMAX_NUM_DECAYMODE_MC]) {
    for (int i = 0; i < DecayMode::MAX_NUM_DECAYMODE; i++) {
        for (int j = 0; j < DecayModeMC::MAX_NUM_DECAYMODE_MC; j++) {

            DecayModeReduced temp_DecayModeReduced;
            DecayModeReducedMC temp_DecayModeReducedMC;

            switch (i) {
            case DecayMode::B2Kc:
                break;
            case DecayMode::B2KcPi0:
                temp_DecayModeReduced = DecayModeReduced::rB2KcPi0;
                break;
            case DecayMode::B2Ks0Pic:
                temp_DecayModeReduced = DecayModeReduced::rB2Ks0Pic;
                break;
            case DecayMode::B2KcPicPic:
                temp_DecayModeReduced = DecayModeReduced::rB2KcPicPic;
                break;
            case DecayMode::B2Ks0PicPi0:
                temp_DecayModeReduced = DecayModeReduced::rB2Ks0PicPi0;
                break;
            case DecayMode::B2KcPicPicPi0:
                temp_DecayModeReduced = DecayModeReduced::rB2KcPicPicPi0;
                break;
            case DecayMode::B2Ks0PicPicPic:
                temp_DecayModeReduced = DecayModeReduced::rB2Ks0PicPicPic;
                break;
            case DecayMode::B2KcPicPicPicPic:
                temp_DecayModeReduced = DecayModeReduced::rB2KcPicPicPicPic;
                break;
            case DecayMode::B2Ks0PicPicPicPi0:
                temp_DecayModeReduced = DecayModeReduced::rB2Ks0PicPicPicPi0;
                break;
            case DecayMode::B2KcPi0Pi0:
                temp_DecayModeReduced = DecayModeReduced::rB2KcPi0Pi0;
                break;
            case DecayMode::B2Ks0PicPi0Pi0:
                temp_DecayModeReduced = DecayModeReduced::rB2Ks0PicPi0Pi0;
                break;
            case DecayMode::B2KcPicPicPi0Pi0:
                temp_DecayModeReduced = DecayModeReduced::rB2KcPicPicPi0Pi0;
                break;
            case DecayMode::B2KcKcKc:
                temp_DecayModeReduced = DecayModeReduced::rB2KcKcKc;
                break;
            case DecayMode::B2KcKcKs0Pic:
                temp_DecayModeReduced = DecayModeReduced::rB2KcKcKs0Pic;
                break;
            case DecayMode::B2KcKcKcPi0:
                temp_DecayModeReduced = DecayModeReduced::rB2KcKcKcPi0;
                break;
            case DecayMode::B02Ks0:
                break;
            case DecayMode::B02KcPic:
                temp_DecayModeReduced = DecayModeReduced::rB02KcPic;
                break;
            case DecayMode::B02Ks0Pi0:
                temp_DecayModeReduced = DecayModeReduced::rB02Ks0Pi0;
                break;
            case DecayMode::B02KcPicPi0:
                temp_DecayModeReduced = DecayModeReduced::rB02KcPicPi0;
                break;
            case DecayMode::B02Ks0PicPic:
                temp_DecayModeReduced = DecayModeReduced::rB02Ks0PicPic;
                break;
            case DecayMode::B02KcPicPicPic:
                temp_DecayModeReduced = DecayModeReduced::rB02KcPicPicPic;
                break;
            case DecayMode::B02Ks0PicPicPi0:
                temp_DecayModeReduced = DecayModeReduced::rB02Ks0PicPicPi0;
                break;
            case DecayMode::B02KcPicPicPicPi0:
                temp_DecayModeReduced = DecayModeReduced::rB02KcPicPicPicPi0;
                break;
            case DecayMode::B02Ks0PicPicPicPic:
                temp_DecayModeReduced = DecayModeReduced::rB02Ks0PicPicPicPic;
                break;
            case DecayMode::B02Ks0Pi0Pi0:
                temp_DecayModeReduced = DecayModeReduced::rB02Ks0Pi0Pi0;
                break;
            case DecayMode::B02KcPicPi0Pi0:
                temp_DecayModeReduced = DecayModeReduced::rB02KcPicPi0Pi0;
                break;
            case DecayMode::B02Ks0PicPicPi0Pi0:
                temp_DecayModeReduced = DecayModeReduced::rB02Ks0PicPicPi0Pi0;
                break;
            case DecayMode::B02KcKcKs0:
                temp_DecayModeReduced = DecayModeReduced::rB02KcKcKs0;
                break;
            case DecayMode::B02KcKcKcPic:
                temp_DecayModeReduced = DecayModeReduced::rB02KcKcKcPic;
                break;
            case DecayMode::B02KcKcKs0Pi0:
                temp_DecayModeReduced = DecayModeReduced::rB02KcKcKs0Pi0;
                break;
            }


            switch (j) {
            case DecayModeMC::Xsu2Kc_MC:
                break;
            case DecayModeMC::Xsu2Kcstar2KcPi0_MC:
                    break;
            case DecayModeMC::Xsu2Kcstar2K0Pic_MC:
                    break;
            case DecayModeMC::Xsu2KcPi0_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsu2KcPi0_MC;
                    break;
            case DecayModeMC::Xsu2K0Pic_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsu2K0Pic_MC;
                    break;
            case DecayModeMC::Xsu2KcPicPic_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsu2KcPicPic_MC;
                    break;
            case DecayModeMC::Xsu2K0PicPi0_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsu2K0PicPi0_MC;
                    break;
            case DecayModeMC::Xsu2KcPicPicPi0_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsu2KcPicPicPi0_MC;
                    break;
            case DecayModeMC::Xsu2K0PicPicPic_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsu2K0PicPicPic_MC;
                    break;
            case DecayModeMC::Xsu2KcPicPicPicPic_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsu2KcPicPicPicPic_MC;
                    break;
            case DecayModeMC::Xsu2K0PicPicPicPi0_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsu2K0PicPicPicPi0_MC;
                    break;
            case DecayModeMC::Xsu2KcPi0Pi0_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsu2KcPi0Pi0_MC;
                    break;
            case DecayModeMC::Xsu2K0PicPi0Pi0_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsu2K0PicPi0Pi0_MC;
                    break;
            case DecayModeMC::Xsu2KcPicPicPi0Pi0_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsu2KcPicPicPi0Pi0_MC;
                    break;
            case DecayModeMC::Xsu2KcKcKc_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsu2KcKcKc_MC
                    break;
            case DecayModeMC::Xsu2KcKcK0Pic_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsu2KcKcK0Pic_MC;
                    break;
            case DecayModeMC::Xsu2KcKcKcPi0_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsu2KcKcKcPi0_MC;
                    break;
            case DecayModeMC::Xsd2K0_MC:
                    break;
            case DecayModeMC::Xsd2K0star2KcPic_MC:
                    break;
            case DecayModeMC::Xsd2K0star2K0Pi0_MC:
                    break;
            case DecayModeMC::Xsd2KcPic_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsd2KcPic_MC;
                    break;
            case DecayModeMC::Xsd2K0Pi0_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsd2K0Pi0_MC;
                    break;
            case DecayModeMC::Xsd2KcPicPi0_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsd2KcPicPi0_MC;
                    break;
            case DecayModeMC::Xsd2K0PicPic_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsd2K0PicPic_MC;
                    break;
            case DecayModeMC::Xsd2KcPicPicPic_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsd2KcPicPicPic_MC;
                    break;
            case DecayModeMC::Xsd2K0PicPicPi0_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsd2K0PicPicPi0_MC;
                    break;
            case DecayModeMC::Xsd2KcPicPicPicPi0_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsd2KcPicPicPicPi0_MC;
                    break;
            case DecayModeMC::Xsd2K0PicPicPicPic_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsd2K0PicPicPicPic_MC
                    break;
            case DecayModeMC::Xsd2K0Pi0Pi0_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsd2K0Pi0Pi0_MC;
                    break;
            case DecayModeMC::Xsd2KcPicPi0Pi0_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsd2KcPicPi0Pi0_MC;
                    break;
            case DecayModeMC::Xsd2K0PicPicPi0Pi0_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsd2K0PicPicPi0Pi0_MC;
                    break;
            case DecayModeMC::Xsd2KcKcK0_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsd2KcKcK0_MC;
                    break;
            case DecayModeMC::Xsd2KcKcKcPic_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsd2KcKcKcPic_MC;
                    break;
            case DecayModeMC::Xsd2KcKcK0Pi0_MC:
                DecayModeReducedMC = DecayModeReducedMC::rXsd2KcKcK0Pi0_MC;
                    break;
            case DecayModeMC::other:
                DecayModeReducedMC = DecayModeReducedMC::rother;
                    break;
            }

            rConfusion[temp_DecayModeReduced][DecayModeReducedMC] = Confusion[i][j];

        }


    }
}

void DecayMatrixToReducedDecayMatrix(double OneDConfusion[DecayMode::MAX_NUM_DECAYMODE], double rOneDConfusion[DecayModeReduced::rMAX_NUM_DECAYMODE]) {
    for (int i = 0; i < DecayMode::MAX_NUM_DECAYMODE; i++) {

            DecayModeReduced temp_DecayModeReduced;

            switch (i) {
            case DecayMode::B2Kc:
                break;
            case DecayMode::B2KcPi0:
                temp_DecayModeReduced = DecayModeReduced::rB2KcPi0;
                break;
            case DecayMode::B2Ks0Pic:
                temp_DecayModeReduced = DecayModeReduced::rB2Ks0Pic;
                break;
            case DecayMode::B2KcPicPic:
                temp_DecayModeReduced = DecayModeReduced::rB2KcPicPic;
                break;
            case DecayMode::B2Ks0PicPi0:
                temp_DecayModeReduced = DecayModeReduced::rB2Ks0PicPi0;
                break;
            case DecayMode::B2KcPicPicPi0:
                temp_DecayModeReduced = DecayModeReduced::rB2KcPicPicPi0;
                break;
            case DecayMode::B2Ks0PicPicPic:
                temp_DecayModeReduced = DecayModeReduced::rB2Ks0PicPicPic;
                break;
            case DecayMode::B2KcPicPicPicPic:
                temp_DecayModeReduced = DecayModeReduced::rB2KcPicPicPicPic;
                break;
            case DecayMode::B2Ks0PicPicPicPi0:
                temp_DecayModeReduced = DecayModeReduced::rB2Ks0PicPicPicPi0;
                break;
            case DecayMode::B2KcPi0Pi0:
                temp_DecayModeReduced = DecayModeReduced::rB2KcPi0Pi0;
                break;
            case DecayMode::B2Ks0PicPi0Pi0:
                temp_DecayModeReduced = DecayModeReduced::rB2Ks0PicPi0Pi0;
                break;
            case DecayMode::B2KcPicPicPi0Pi0:
                temp_DecayModeReduced = DecayModeReduced::rB2KcPicPicPi0Pi0;
                break;
            case DecayMode::B2KcKcKc:
                temp_DecayModeReduced = DecayModeReduced::rB2KcKcKc;
                break;
            case DecayMode::B2KcKcKs0Pic:
                temp_DecayModeReduced = DecayModeReduced::rB2KcKcKs0Pic;
                break;
            case DecayMode::B2KcKcKcPi0:
                temp_DecayModeReduced = DecayModeReduced::rB2KcKcKcPi0;
                break;
            case DecayMode::B02Ks0:
                break;
            case DecayMode::B02KcPic:
                temp_DecayModeReduced = DecayModeReduced::rB02KcPic;
                break;
            case DecayMode::B02Ks0Pi0:
                temp_DecayModeReduced = DecayModeReduced::rB02Ks0Pi0;
                break;
            case DecayMode::B02KcPicPi0:
                temp_DecayModeReduced = DecayModeReduced::rB02KcPicPi0;
                break;
            case DecayMode::B02Ks0PicPic:
                temp_DecayModeReduced = DecayModeReduced::rB02Ks0PicPic;
                break;
            case DecayMode::B02KcPicPicPic:
                temp_DecayModeReduced = DecayModeReduced::rB02KcPicPicPic;
                break;
            case DecayMode::B02Ks0PicPicPi0:
                temp_DecayModeReduced = DecayModeReduced::rB02Ks0PicPicPi0;
                break;
            case DecayMode::B02KcPicPicPicPi0:
                temp_DecayModeReduced = DecayModeReduced::rB02KcPicPicPicPi0;
                break;
            case DecayMode::B02Ks0PicPicPicPic:
                temp_DecayModeReduced = DecayModeReduced::rB02Ks0PicPicPicPic;
                break;
            case DecayMode::B02Ks0Pi0Pi0:
                temp_DecayModeReduced = DecayModeReduced::rB02Ks0Pi0Pi0;
                break;
            case DecayMode::B02KcPicPi0Pi0:
                temp_DecayModeReduced = DecayModeReduced::rB02KcPicPi0Pi0;
                break;
            case DecayMode::B02Ks0PicPicPi0Pi0:
                temp_DecayModeReduced = DecayModeReduced::rB02Ks0PicPicPi0Pi0;
                break;
            case DecayMode::B02KcKcKs0:
                temp_DecayModeReduced = DecayModeReduced::rB02KcKcKs0;
                break;
            case DecayMode::B02KcKcKcPic:
                temp_DecayModeReduced = DecayModeReduced::rB02KcKcKcPic;
                break;
            case DecayMode::B02KcKcKs0Pi0:
                temp_DecayModeReduced = DecayModeReduced::rB02KcKcKs0Pi0;
                break;
            }

            rOneDConfusion[temp_DecayModeReduced] = OneDConfusion[i];

    }
}

void ExtractReducedElement(double Confusion[DecayMode::MAX_NUM_DECAYMODE][DecayModeMC::MAX_NUM_DECAYMODE_MC], double rOneDConfusion_KKstar[DecayModeReduced::rMAX_NUM_DECAYMODE]) {
    for (int i = 0; i < DecayMode::MAX_NUM_DECAYMODE; i++) {

        DecayModeReduced temp_DecayModeReduced;

        switch (i) {
        case DecayMode::B2Kc:
            break;
        case DecayMode::B2KcPi0:
            temp_DecayModeReduced = DecayModeReduced::rB2KcPi0;
            break;
        case DecayMode::B2Ks0Pic:
            temp_DecayModeReduced = DecayModeReduced::rB2Ks0Pic;
            break;
        case DecayMode::B2KcPicPic:
            temp_DecayModeReduced = DecayModeReduced::rB2KcPicPic;
            break;
        case DecayMode::B2Ks0PicPi0:
            temp_DecayModeReduced = DecayModeReduced::rB2Ks0PicPi0;
            break;
        case DecayMode::B2KcPicPicPi0:
            temp_DecayModeReduced = DecayModeReduced::rB2KcPicPicPi0;
            break;
        case DecayMode::B2Ks0PicPicPic:
            temp_DecayModeReduced = DecayModeReduced::rB2Ks0PicPicPic;
            break;
        case DecayMode::B2KcPicPicPicPic:
            temp_DecayModeReduced = DecayModeReduced::rB2KcPicPicPicPic;
            break;
        case DecayMode::B2Ks0PicPicPicPi0:
            temp_DecayModeReduced = DecayModeReduced::rB2Ks0PicPicPicPi0;
            break;
        case DecayMode::B2KcPi0Pi0:
            temp_DecayModeReduced = DecayModeReduced::rB2KcPi0Pi0;
            break;
        case DecayMode::B2Ks0PicPi0Pi0:
            temp_DecayModeReduced = DecayModeReduced::rB2Ks0PicPi0Pi0;
            break;
        case DecayMode::B2KcPicPicPi0Pi0:
            temp_DecayModeReduced = DecayModeReduced::rB2KcPicPicPi0Pi0;
            break;
        case DecayMode::B2KcKcKc:
            temp_DecayModeReduced = DecayModeReduced::rB2KcKcKc;
            break;
        case DecayMode::B2KcKcKs0Pic:
            temp_DecayModeReduced = DecayModeReduced::rB2KcKcKs0Pic;
            break;
        case DecayMode::B2KcKcKcPi0:
            temp_DecayModeReduced = DecayModeReduced::rB2KcKcKcPi0;
            break;
        case DecayMode::B02Ks0:
            break;
        case DecayMode::B02KcPic:
            temp_DecayModeReduced = DecayModeReduced::rB02KcPic;
            break;
        case DecayMode::B02Ks0Pi0:
            temp_DecayModeReduced = DecayModeReduced::rB02Ks0Pi0;
            break;
        case DecayMode::B02KcPicPi0:
            temp_DecayModeReduced = DecayModeReduced::rB02KcPicPi0;
            break;
        case DecayMode::B02Ks0PicPic:
            temp_DecayModeReduced = DecayModeReduced::rB02Ks0PicPic;
            break;
        case DecayMode::B02KcPicPicPic:
            temp_DecayModeReduced = DecayModeReduced::rB02KcPicPicPic;
            break;
        case DecayMode::B02Ks0PicPicPi0:
            temp_DecayModeReduced = DecayModeReduced::rB02Ks0PicPicPi0;
            break;
        case DecayMode::B02KcPicPicPicPi0:
            temp_DecayModeReduced = DecayModeReduced::rB02KcPicPicPicPi0;
            break;
        case DecayMode::B02Ks0PicPicPicPic:
            temp_DecayModeReduced = DecayModeReduced::rB02Ks0PicPicPicPic;
            break;
        case DecayMode::B02Ks0Pi0Pi0:
            temp_DecayModeReduced = DecayModeReduced::rB02Ks0Pi0Pi0;
            break;
        case DecayMode::B02KcPicPi0Pi0:
            temp_DecayModeReduced = DecayModeReduced::rB02KcPicPi0Pi0;
            break;
        case DecayMode::B02Ks0PicPicPi0Pi0:
            temp_DecayModeReduced = DecayModeReduced::rB02Ks0PicPicPi0Pi0;
            break;
        case DecayMode::B02KcKcKs0:
            temp_DecayModeReduced = DecayModeReduced::rB02KcKcKs0;
            break;
        case DecayMode::B02KcKcKcPic:
            temp_DecayModeReduced = DecayModeReduced::rB02KcKcKcPic;
            break;
        case DecayMode::B02KcKcKs0Pi0:
            temp_DecayModeReduced = DecayModeReduced::rB02KcKcKs0Pi0;
            break;
        }

        rOneDConfusion_KKstar[temp_DecayModeReduced]
            = Confusion[i][DecayModeMC::Xsu2Kc_MC]
            + Confusion[i][DecayModeMC::Xsu2Kcstar2KcPi0_MC]
            + Confusion[i][DecayModeMC::Xsu2Kcstar2K0Pic_MC]
            + Confusion[i][DecayModeMC::Xsd2K0_MC]
            + Confusion[i][DecayModeMC::Xsd2K0star2KcPic_MC]
            + Confusion[i][DecayModeMC::Xsd2K0star2K0Pi0_MC];

    }
}

void SetError(double OneDConfusion_data[DecayMode::MAX_NUM_DECAYMODE], double Error_OneDConfusion_data[DecayMode::MAX_NUM_DECAYMODE]) {
    for (int i = 0; i < DecayMode::MAX_NUM_DECAYMODE; i++) {
        if (OneDConfusion_data[i] < MyEPSILON) Error_OneDConfusion_data[i] = 1.0;
        else Error_OneDConfusion_data[i] = std::sqrt(OneDConfusion_data[i]);
    }
}

void PrintConfusionMatrix(double Confusion[DecayMode::MAX_NUM_DECAYMODE][DecayModeMC::MAX_NUM_DECAYMODE_MC], double OneDConfusion[DecayMode::MAX_NUM_DECAYMODE], double OneDConfusion_data[DecayMode::MAX_NUM_DECAYMODE]) {
    printf("--------------- confusion matrix for XsJ/psi signal ---------------\n");
    printf("[");
    for (int i = 0; i < DecayMode::MAX_NUM_DECAYMODE; i++) {
        printf("[");
        for (int j = 0; j < DecayModeMC::MAX_NUM_DECAYMODE_MC; j++) {
            if (j != DecayModeMC::MAX_NUM_DECAYMODE_MC - 1) printf("%f,", Confusion[i][j]);
            else printf("%f", Confusion[i][j]);
        }
        if (i != DecayMode::MAX_NUM_DECAYMODE - 1) printf("], ");
        else printf("]");
    }
    printf("]\n");
    printf("--------------- confusion matrix for XsJ/psi signal ---------------\n");

    printf("--------------- 1D confusion matrix for background ---------------\n");
    printf("[");
    for (int i = 0; i < DecayMode::MAX_NUM_DECAYMODE; i++) {
        if (i != DecayMode::MAX_NUM_DECAYMODE) printf("[%f],", OneDConfusion[i]);
        else printf("[%f]", OneDConfusion[i]);
    }
    printf("]\n");
    printf("--------------- 1D confusion matrix for background ---------------\n");

    printf("--------------- 1D confusion matrix for data ---------------\n");
    printf("[");
    for (int i = 0; i < DecayMode::MAX_NUM_DECAYMODE; i++) {
        if (i != DecayMode::MAX_NUM_DECAYMODE) printf("[%f],", OneDConfusion_data[i]);
        else printf("[%f]", OneDConfusion_data[i]);
    }
    printf("]\n");
    printf("--------------- 1D confusion matrix for data ---------------\n");
}

void PrintReducedConfusionMatrix(double rConfusion[DecayModeReduced::rMAX_NUM_DECAYMODE][DecayModeReducedMC::rMAX_NUM_DECAYMODE_MC], double rOneDConfusion[DecayModeReduced::rMAX_NUM_DECAYMODE], double rOneDConfusion_data[DecayModeReduced::rMAX_NUM_DECAYMODE]) {
    printf("--------------- reduced confusion matrix for XsJ/psi signal ---------------\n");
    printf("[");
    for (int i = 0; i < DecayModeReduced::rMAX_NUM_DECAYMODE; i++) {
        printf("[");
        for (int j = 0; j < DecayModeReducedMC::rMAX_NUM_DECAYMODE_MC; j++) {
            if (j != DecayModeReducedMC::rMAX_NUM_DECAYMODE_MC - 1) printf("%f,", rConfusion[i][j]);
            else printf("%f", rConfusion[i][j]);
        }
        if (i != DecayModeReduced::rMAX_NUM_DECAYMODE - 1) printf("], ");
        else printf("]");
    }
    printf("]\n");
    printf("--------------- reduced confusion matrix for XsJ/psi signal ---------------\n");

    printf("--------------- reduced 1D confusion matrix for background ---------------\n");
    printf("[");
    for (int i = 0; i < DecayModeReduced::rMAX_NUM_DECAYMODE; i++) {
        if (i != DecayModeReduced::rMAX_NUM_DECAYMODE) printf("[%f],", rOneDConfusion[i]);
        else printf("[%f]", rOneDConfusion[i]);
    }
    printf("]\n");
    printf("--------------- reduced 1D confusion matrix for background ---------------\n");

    printf("--------------- reduced 1D confusion matrix for data ---------------\n");
    printf("[");
    for (int i = 0; i < DecayModeReduced::rMAX_NUM_DECAYMODE; i++) {
        if (i != DecayModeReduced::rMAX_NUM_DECAYMODE) printf("[%f],", rOneDConfusion_data[i]);
        else printf("[%f]", rOneDConfusion_data[i]);
    }
    printf("]\n");
    printf("--------------- reduced 1D confusion matrix for data ---------------\n");
}

void PrintFinalResults(double Eff[DecayMode::MAX_NUM_DECAYMODE][DecayModeMC::MAX_NUM_DECAYMODE_MC], double OneDConfusion[DecayMode::MAX_NUM_DECAYMODE], double OneDConfusion_data[DecayMode::MAX_NUM_DECAYMODE], double Error_OneDConfusion_data[DecayMode::MAX_NUM_DECAYMODE]) {

    printf("--------------- efficiency matrix for XsJ/psi signal ---------------\n");
    printf("[");
    for (int i = 0; i < DecayMode::MAX_NUM_DECAYMODE; i++) {
        printf("[");
        for (int j = 0; j < DecayModeMC::MAX_NUM_DECAYMODE_MC; j++) {
            if (j != DecayModeMC::MAX_NUM_DECAYMODE_MC - 1) printf("%.9lf,", Eff[i][j]);
            else printf("%.9lf", Eff[i][j]);
        }
        if (i != DecayMode::MAX_NUM_DECAYMODE - 1) printf("], ");
        else printf("]");
    }
    printf("]\n");
    printf("--------------- efficiency matrix for XsJ/psi signal ---------------\n");

    printf("--------------- 1D confusion matrix for data - MC BKG ---------------\n");
    printf("[");
    for (int i = 0; i < DecayMode::MAX_NUM_DECAYMODE; i++) {
        if (i != DecayMode::MAX_NUM_DECAYMODE - 1) printf("[%f],", OneDConfusion_data[i] - OneDConfusion[i]);
        else printf("[%f]", OneDConfusion_data[i] - OneDConfusion[i]);
    }
    printf("]\n");
    printf("--------------- 1D confusion matrix for data - MC BKG ---------------\n");

    printf("--------------- 1D confusion matrix for data error ---------------\n");
    printf("[");
    for (int i = 0; i < DecayMode::MAX_NUM_DECAYMODE; i++) {
        if (i != DecayMode::MAX_NUM_DECAYMODE - 1) printf("[%f],", Error_OneDConfusion_data[i]);
        else printf("[%f]", Error_OneDConfusion_data[i]);
    }
    printf("]\n");
    printf("--------------- 1D confusion matrix for data error ---------------\n");

}

void PrintReducedFinalResults(double rEff[DecayModeReduced::rMAX_NUM_DECAYMODE][DecayModeReducedMC::rMAX_NUM_DECAYMODE_MC], double rOneDConfusion_KKstar[DecayModeReduced::rMAX_NUM_DECAYMODE], double rOneDConfusion[DecayModeReduced::rMAX_NUM_DECAYMODE], double rOneDConfusion_data[DecayModeReduced::rMAX_NUM_DECAYMODE], double rError_OneDConfusion_data[DecayModeReduced::rMAX_NUM_DECAYMODE]) {

    printf("--------------- reduced efficiency matrix for XsJ/psi signal ---------------\n");
    printf("[");
    for (int i = 0; i < DecayModeReduced::rMAX_NUM_DECAYMODE; i++) {
        printf("[");
        for (int j = 0; j < DecayModeReducedMC::rMAX_NUM_DECAYMODE_MC; j++) {
            if (j != DecayModeReducedMC::rMAX_NUM_DECAYMODE_MC - 1) printf("%.9lf,", rEff[i][j]);
            else printf("%.9lf", rEff[i][j]);
        }
        if (i != DecayModeReduced::rMAX_NUM_DECAYMODE - 1) printf("], ");
        else printf("]");
    }
    printf("]\n");
    printf("--------------- reduced efficiency matrix for XsJ/psi signal ---------------\n");

    printf("--------------- reduced 1D confusion matrix for data - MC BKG - KKstar ---------------\n");
    printf("[");
    for (int i = 0; i < DecayModeReduced::rMAX_NUM_DECAYMODE; i++) {
        if (i != DecayModeReduced::rMAX_NUM_DECAYMODE - 1) printf("[%f],", rOneDConfusion_data[i] - rOneDConfusion[i] - rOneDConfusion_KKstar[i]);
        else printf("[%f]", rOneDConfusion_data[i] - rOneDConfusion[i] - rOneDConfusion_KKstar[i]);
    }
    printf("]\n");
    printf("--------------- reduced 1D confusion matrix for data - MC BKG - KKstar ---------------\n");

    printf("--------------- reduced 1D confusion matrix for data error ---------------\n");
    printf("[");
    for (int i = 0; i < DecayModeReduced::rMAX_NUM_DECAYMODE; i++) {
        if (i != DecayModeReduced::rMAX_NUM_DECAYMODE - 1) printf("[%f],", rError_OneDConfusion_data[i]);
        else printf("[%f]", rError_OneDConfusion_data[i]);
    }
    printf("]\n");
    printf("--------------- reduced 1D confusion matrix for data error ---------------\n");

}

void fragmentation_helper() {

    ReadPIDFile();
    ReadFakePIDFile();

    // define confusion matrix
    double Confusion[DecayMode::MAX_NUM_DECAYMODE][DecayModeMC::MAX_NUM_DECAYMODE_MC] = { 0.0 }; // [reco][MC truth] SIGNAL
    double OneDConfusion[DecayMode::MAX_NUM_DECAYMODE] = { 0.0 }; // [reco] CHG MIX UUBAR DDBAR SSBAR CHARM
    double OneDConfusion_data[DecayMode::MAX_NUM_DECAYMODE] = { 0.0 }; // [reco] data
    double Error_OneDConfusion_data[DecayMode::MAX_NUM_DECAYMODE] = { 0.0 }; // [reco] data error

    double rConfusion[DecayModeReduced::rMAX_NUM_DECAYMODE][DecayModeReducedMC::rMAX_NUM_DECAYMODE_MC] = { 0.0 }; // [reco][MC truth] SIGNAL, remove resonant component
    double rOneDConfusion[DecayModeReduced::rMAX_NUM_DECAYMODE] = { 0.0 }; // [reco] CHG MIX UUBAR DDBAR SSBAR CHARM
    double rOneDConfusion_data[DecayModeReduced::rMAX_NUM_DECAYMODE] = { 0.0 }; // [reco] data
    double rError_OneDConfusion_data[DecayModeReduced::rMAX_NUM_DECAYMODE] = { 0.0 }; // [reco] data error
    double rOneDConfusion_KKstar[DecayModeReduced::rMAX_NUM_DECAYMODE] = { 0.0 }; // [reco][MC truth] SIGNAL, only resonant component

    // dirnames
    const char* Jpsi_MC_SIGNAL_dirname = "/home/jwpark/storage/BKG_gbasf2/Aunn_LS_MC_Jpsi/SIGNAL_analysis/validation_v002/final_output";
    const char* Jpsi_MC_CHG_dirname = "/home/jwpark/storage/BKG_gbasf2/Aunn_LS_MC_Jpsi/CHG_analysis/validation_v002/final_output";
    const char* Jpsi_MC_MIX_dirname = "/home/jwpark/storage/BKG_gbasf2/Aunn_LS_MC_Jpsi/MIX_analysis/validation_v002/final_output";
    const char* Jpsi_MC_UUBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Aunn_LS_MC_Jpsi/UUBAR_analysis/validation_v002/final_output";
    const char* Jpsi_MC_DDBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Aunn_LS_MC_Jpsi/DDBAR_analysis/validation_v002/final_output";
    const char* Jpsi_MC_SSBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Aunn_LS_MC_Jpsi/SSBAR_analysis/validation_v002/final_output";
    const char* Jpsi_MC_CHARM_dirname = "/home/jwpark/storage/BKG_gbasf2/Aunn_LS_MC_Jpsi/CHARM_analysis/validation_v002/final_output";
    const char* Jpsi_data_dirname = "/home/jwpark/storage/BKG_gbasf2/Aunn_LS_data_Jpsi/SIGNAL_analysis/validation_v002/final_output";

    LetsFillJpsi_ri(Jpsi_MC_SIGNAL_dirname, Confusion, "SIGNAL");
    LetsFillJpsi_ri(Jpsi_MC_CHG_dirname, OneDConfusion, "CHG");
    LetsFillJpsi_ri(Jpsi_MC_MIX_dirname, OneDConfusion, "MIX");
    LetsFillJpsi_ri(Jpsi_MC_UUBAR_dirname, OneDConfusion, "UUBAR");
    LetsFillJpsi_ri(Jpsi_MC_DDBAR_dirname, OneDConfusion, "DDBAR");
    LetsFillJpsi_ri(Jpsi_MC_SSBAR_dirname, OneDConfusion, "SSBAR");
    LetsFillJpsi_ri(Jpsi_MC_CHARM_dirname, OneDConfusion, "CHARM");
    LetsFill(Jpsi_data_dirname, OneDConfusion_data);

    SetError(OneDConfusion_data, Error_OneDConfusion_data);

    DecayMatrixToReducedDecayMatrix(Confusion, rConfusion);
    DecayMatrixToReducedDecayMatrix(OneDConfusion, rOneDConfusion);
    DecayMatrixToReducedDecayMatrix(OneDConfusion_data, rOneDConfusion_data);
    DecayMatrixToReducedDecayMatrix(Error_OneDConfusion_data, rError_OneDConfusion_data);

    ExtractReducedElement(Confusion, rOneDConfusion_KKstar);

    //PrintConfusionMatrix(Confusion, OneDConfusion, OneDConfusion_data);
    PrintReducedConfusionMatrix(rConfusion, rOneDConfusion, rOneDConfusion_data);

    // get efficiency matrix
    double Eff[DecayMode::MAX_NUM_DECAYMODE][DecayModeMC::MAX_NUM_DECAYMODE_MC] = { 0.0 }; // [reco][MC truth]
    double rEff[DecayModeReduced::rMAX_NUM_DECAYMODE][DecayModeReducedMC::rMAX_NUM_DECAYMODE_MC] = { 0.0 }; // [reco][MC truth]

    for (int i = 0; i < DecayMode::MAX_NUM_DECAYMODE; i++) {
        for (int j = 0; j < DecayModeMC::MAX_NUM_DECAYMODE_MC; j++) {
            Eff[i][j] = Confusion[i][j] / N_initial_MC[j];
        }
    }

    for (int i = 0; i < DecayModeReduced::rMAX_NUM_DECAYMODE; i++) {
        for (int j = 0; j < DecayModeReducedMC::rMAX_NUM_DECAYMODE_MC; j++) {
            rEff[i][j] = rConfusion[i][j] / rN_initial_MC[j];
        }
    }

    //PrintFinalResults(Eff, OneDConfusion, OneDConfusion_data, Error_OneDConfusion_data);
    PrintReducedFinalResults(rEff, rOneDConfusion_KKstar, rOneDConfusion, rOneDConfusion_data, rError_OneDConfusion_data);

}
