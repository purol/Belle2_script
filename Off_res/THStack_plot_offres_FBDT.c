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

# define FEI_cal_Bc_num 12
# define FEI_cal_B0_num 11
double FEI_cal_Bc[FEI_cal_Bc_num] = { 1.04, 0.79, 0.69, 0.56, 0.97, 0.95, 0.74, 0.57, 0.91, 0.51, 0.34, 0.59 };
double FEI_cal_Bc_uncertainty[FEI_cal_Bc_num] = { 0.03, 0.03, 0.05, 0.11, 0.03, 0.03, 0.02, 0.06, 0.1, 0.13, 0.07, 0.02 }; // not relative uncertainty. absolute uncertainty
double FEI_cal_Bc_modeID[FEI_cal_Bc_num] = { 0.0, 1.0, 3.0, 4.0, 15.0, 16.0, 18.0, 19.0, 23.0, 24.0, 30.0, -1.0 };
double FEI_cal_B0[FEI_cal_B0_num] = { 1.16, 0.94, 0.81, 0.79, 0.99, 1.03, 0.67, 0.66, 0.69, 0.49, 0.79 };
double FEI_cal_B0_uncertainty[FEI_cal_B0_num] = { 0.04, 0.05, 0.06, 0.02, 0.03, 0.06, 0.02, 0.03, 0.02, 0.02, 0.12 }; // not relative uncertainty. absolute uncertainty
double FEI_cal_B0_modeID[FEI_cal_B0_num] = { 0.0, 1.0, 3.0, 4.0, 5.0, 15.0, 16.0, 18.0, 19.0, 26.0, -1.0 };

# define Nvar_num 1

# define CAL 1.0
# define CAL_qq 1.0
# define Stream 0.25

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

double GetFEICalFactor(double UpsilonID, double BtagID) {
    // UpsilonID => charged: 0, mixed: 1

    if (UpsilonID > -0.5 && UpsilonID < 0.5) { // charged
        for (int i = 0; i < FEI_cal_Bc_num - 1; i++) {
            if (BtagID > FEI_cal_Bc_modeID[i] - 0.5 && BtagID < FEI_cal_Bc_modeID[i] + 0.5) return FEI_cal_Bc[i];
        }
        return FEI_cal_Bc[FEI_cal_Bc_num - 1];
    }
    else if (UpsilonID > 0.5 && UpsilonID < 1.5) { // mixed
        for (int i = 0; i < FEI_cal_B0_num - 1; i++) {
            if (BtagID > FEI_cal_B0_modeID[i] - 0.5 && BtagID < FEI_cal_B0_modeID[i] + 0.5) return FEI_cal_B0[i];
        }
        return FEI_cal_B0[FEI_cal_B0_num - 1];
    }

    printf("[GetFEICalFactor] error! unexpected decay ID\n");
    exit(1);
    return 0;

}

double GetFEICalFactorUncer(double UpsilonID, double BtagID) {
    // UpsilonID => charged: 0, mixed: 1

    if (UpsilonID > -0.5 && UpsilonID < 0.5) { // charged
        for (int i = 0; i < FEI_cal_Bc_num - 1; i++) {
            if (BtagID > FEI_cal_Bc_modeID[i] - 0.5 && BtagID < FEI_cal_Bc_modeID[i] + 0.5) return FEI_cal_Bc_uncertainty[i];
        }
        return FEI_cal_Bc_uncertainty[FEI_cal_Bc_num - 1];
    }
    else if (UpsilonID > 0.5 && UpsilonID < 1.5) { // mixed
        for (int i = 0; i < FEI_cal_B0_num - 1; i++) {
            if (BtagID > FEI_cal_B0_modeID[i] - 0.5 && BtagID < FEI_cal_B0_modeID[i] + 0.5) return FEI_cal_B0[i];
        }
        return FEI_cal_B0_uncertainty[FEI_cal_B0_num - 1];
    }

    printf("[GetFEICalFactor] error! unexpected decay ID\n");
    exit(1);
    return 0;

}

double BDTcToWeight(double BDTc) {

    if (BDTc > (5.0 / 6.0)) return 5.0;
    else return (BDTc / (1.0 - BDTc));

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

void LetsFillOffres_ri(const char* dirname, std::vector<std::string> variable_names, std::vector<std::string> branch_names, std::vector<double> variable_values[Nvar_num], std::vector<int>* numberings, std::vector<double>* weights, std::string SampleName) {
    /*
    0: uubar
    1: ddbar
    2: ssbar
    3: ccbar
    4: tautau
    5: mumu
    6: gg
    7: ee
    8: eeee
    9: eemumu
    10: llXX
    11: hhISR
    */

    float var[Nvar_num] = { 0.0 };
    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    double FEI_calibration_factor = -1;

    double Ngamma_v200 = -1;

    std::vector<string> names;
    load_files(dirname, &names);

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
        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            for (int k = 0; k < (int)variable_names.size(); k++) variable_values[k].push_back(var[k]);

            // Fill numberings
            double weight_ri = 0.0;
            if (SampleName == "UUBAR") {
                numberings->push_back(0);
                FEI_calibration_factor = CAL_qq;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i)) * (0.042329 / 0.361673); // scale to 364.436 and scale to 42.329
            }
            else if (SampleName == "DDBAR") {
                numberings->push_back(1);
                FEI_calibration_factor = CAL_qq;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i)) * (0.042329 / 0.361673); // scale to 364.436 and scale to 42.329
            }
            else if (SampleName == "SSBAR") {
                numberings->push_back(2);
                FEI_calibration_factor = CAL_qq;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i)) * (0.042329 / 0.361673); // scale to 364.436 and scale to 42.329
            }
            else if (SampleName == "CHARM") {
                numberings->push_back(3);
                FEI_calibration_factor = CAL_qq;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i)) * (0.042329 / 0.361673); // scale to 364.436 and scale to 42.329
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

            // Multiplicity correction factor, it is not applied now. it is for a systematic uncertainty
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_multiplicity);


        }
        input_file->Close();

    }

}

void LetsFill(const char* dirname, std::vector<std::string> variable_names, std::vector<std::string> branch_names, std::vector<double> variable_values[Nvar_num]) {
    float var[Nvar_num] = { 0.0 };

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        for (int k = 0; k < (int)variable_names.size(); k++) {
            if(branch_names.at(k) == std::string("Upsilon")) tree_upsilon->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
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

void LetsFill(const char* dirname, std::vector<std::string> variable_names, std::vector<std::string> branch_names, std::vector<double> variable_values[Nvar_num], const char* included_string) {
    float var[Nvar_num] = { 0.0 };

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

void LetsFillOffres_ri_correction(const char* dirname, std::vector<std::string> variable_names, std::vector<std::string> branch_names, std::vector<double> variable_values[Nvar_num], std::vector<int>* numberings, std::vector<double>* weights, std::string SampleName, double NormFactor = 1.0) {
    /*
    0: uubar
    1: ddbar
    2: ssbar
    3: ccbar
    4: tautau
    5: mumu
    6: gg
    7: ee
    8: eeee
    9: eemumu
    10: llXX
    11: hhISR
    */

    float var[Nvar_num] = { 0.0 };
    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    double FEI_calibration_factor = -1;

    double Ngamma_v200 = -1;

    float BDTc = -1;
    double BDTc_correction = -1;

    std::vector<string> names;
    load_files(dirname, &names);

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
        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);
        tree_upsilon->SetBranchAddress("MVA_Continuum", &BDTc);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            // BDTc correction factor
            BDTc_correction = BDTcToWeight(BDTc) * NormFactor;

            for (int k = 0; k < (int)variable_names.size(); k++) variable_values[k].push_back(var[k]);

            // Fill numberings
            double weight_ri = 0.0;
            if (SampleName == "UUBAR") {
                numberings->push_back(0);
                FEI_calibration_factor = CAL_qq;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i)) * (0.042329 / 0.361673); // scale to 364.436 and scale to 42.329
            }
            else if (SampleName == "DDBAR") {
                numberings->push_back(1);
                FEI_calibration_factor = CAL_qq;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i)) * (0.042329 / 0.361673); // scale to 364.436 and scale to 42.329
            }
            else if (SampleName == "SSBAR") {
                numberings->push_back(2);
                FEI_calibration_factor = CAL_qq;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i)) * (0.042329 / 0.361673); // scale to 364.436 and scale to 42.329
            }
            else if (SampleName == "CHARM") {
                numberings->push_back(3);
                FEI_calibration_factor = CAL_qq;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i)) * (0.042329 / 0.361673); // scale to 364.436 and scale to 42.329
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

            // Multiplicity correction factor, it is not applied now. it is for a systematic uncertainty
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_multiplicity * BDTc_correction);


        }
        input_file->Close();

    }

}

typedef struct _Nevt {
    double NevtwithoutCorrection;
    double NevtwithCorrection;
} Nevt;

void NevtCount_ri(const char* dirname, std::string SampleName, Nevt* nevt) {
    /*
    0: uubar
    1: ddbar
    2: ssbar
    3: ccbar
    4: tautau
    5: mumu
    6: gg
    7: ee
    8: eeee
    9: eemumu
    10: llXX
    11: hhISR
    */

    float var[Nvar_num] = { 0.0 };
    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    double FEI_calibration_factor = -1;

    double Ngamma_v200 = -1;

    float BDTc = -1;
    double BDTc_correction = -1;

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
        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);
        tree_upsilon->SetBranchAddress("MVA_Continuum", &BDTc);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            // BDTc correction factor
            BDTc_correction = BDTcToWeight(BDTc);

            // Fill numberings
            double weight_ri = 0.0;
            if (SampleName == "UUBAR") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i)) * (0.042329 / 0.361673); // scale to 364.436 and scale to 42.329
            }
            else if (SampleName == "DDBAR") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i)) * (0.042329 / 0.361673); // scale to 364.436 and scale to 42.329
            }
            else if (SampleName == "SSBAR") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i)) * (0.042329 / 0.361673); // scale to 364.436 and scale to 42.329
            }
            else if (SampleName == "CHARM") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i)) * (0.042329 / 0.361673); // scale to 364.436 and scale to 42.329
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

            // Multiplicity correction factor, it is not applied now. it is for a systematic uncertainty
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_multiplicity;
            nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_multiplicity * BDTc_correction;


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

void THStack_plot_offres_FBDT() {


    Nevt nevt_UUBAR = { 0.0, 0.0 };
    Nevt nevt_DDBAR = { 0.0, 0.0 };
    Nevt nevt_SSBAR = { 0.0, 0.0 };
    Nevt nevt_CHARM = { 0.0, 0.0 };

    const char* Offres_MC_UUBAR_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_LS_MC_off/UUBAR_analysis/test_v000/final_output";
    const char* Offres_MC_DDBAR_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_LS_MC_off/DDBAR_analysis/test_v000/final_output";
    const char* Offres_MC_SSBAR_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_LS_MC_off/SSBAR_analysis/test_v000/final_output";
    const char* Offres_MC_CHARM_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_LS_MC_off/CHARM_analysis/test_v000/final_output";
    const char* Offres_data_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_LS_data_off/SIGNAL_analysis/validation_v000/final_output";

    NevtCount_ri(Offres_MC_UUBAR_dirname, "UUBAR", &nevt_UUBAR);
    NevtCount_ri(Offres_MC_DDBAR_dirname, "DDBAR", &nevt_DDBAR);
    NevtCount_ri(Offres_MC_SSBAR_dirname, "SSBAR", &nevt_SSBAR);
    NevtCount_ri(Offres_MC_CHARM_dirname, "CHARM", &nevt_CHARM);

    double NormFactor_UUBAR = nevt_UUBAR.NevtwithoutCorrection / nevt_UUBAR.NevtwithCorrection;
    double NormFactor_DDBAR = nevt_DDBAR.NevtwithoutCorrection / nevt_DDBAR.NevtwithCorrection;
    double NormFactor_SSBAR = nevt_SSBAR.NevtwithoutCorrection / nevt_SSBAR.NevtwithCorrection;
    double NormFactor_CHARM = nevt_CHARM.NevtwithoutCorrection / nevt_CHARM.NevtwithCorrection;

    std::vector<std::string> variable_names;
    std::vector<std::string> branch_names;
    variable_names.push_back("MVA_BB"); branch_names.push_back("Upsilon");

    int Nvar = static_cast<int>(variable_names.size());
    if (Nvar != Nvar_num) exit(1);
    std::vector<double> Offres_MC_values[Nvar_num];
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
    std::vector<int> Offres_MC_numbering;

    std::vector<double> Offres_data_values[Nvar_num];

    std::vector<double> weights;
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

    /*
    LetsFillOffres_ri_correction(Offres_MC_UUBAR_dirname, variable_names, branch_names, Offres_MC_values, &Offres_MC_numbering, &weights, "UUBAR", NormFactor_UUBAR);
    LetsFillOffres_ri_correction(Offres_MC_DDBAR_dirname, variable_names, branch_names, Offres_MC_values, &Offres_MC_numbering, &weights, "DDBAR", NormFactor_DDBAR);
    LetsFillOffres_ri_correction(Offres_MC_SSBAR_dirname, variable_names, branch_names, Offres_MC_values, &Offres_MC_numbering, &weights, "SSBAR", NormFactor_SSBAR);
    LetsFillOffres_ri_correction(Offres_MC_CHARM_dirname, variable_names, branch_names, Offres_MC_values, &Offres_MC_numbering, &weights, "CHARM", NormFactor_CHARM);
    */
    LetsFillOffres_ri(Offres_MC_UUBAR_dirname, variable_names, branch_names, Offres_MC_values, &Offres_MC_numbering, &weights, "UUBAR");
    LetsFillOffres_ri(Offres_MC_DDBAR_dirname, variable_names, branch_names, Offres_MC_values, &Offres_MC_numbering, &weights, "DDBAR");
    LetsFillOffres_ri(Offres_MC_SSBAR_dirname, variable_names, branch_names, Offres_MC_values, &Offres_MC_numbering, &weights, "SSBAR");
    LetsFillOffres_ri(Offres_MC_CHARM_dirname, variable_names, branch_names, Offres_MC_values, &Offres_MC_numbering, &weights, "CHARM");
    LetsFill(Offres_data_dirname, variable_names, branch_names, Offres_data_values);

    // sort variables
    for (int k = 0; k < (int)Offres_MC_numbering.size(); k++) {
        if (Offres_MC_numbering.at(k) == 0) {
            for (int l = 0; l < (int)variable_names.size(); l++) uubar_values[l].push_back(Offres_MC_values[l].at(k));
            uubar_weights.push_back(weights.at(k));
        }
        else if (Offres_MC_numbering.at(k) == 1) {
            for (int l = 0; l < (int)variable_names.size(); l++) ddbar_values[l].push_back(Offres_MC_values[l].at(k));
            ddbar_weights.push_back(weights.at(k));
        }
        else if (Offres_MC_numbering.at(k) == 2) {
            for (int l = 0; l < (int)variable_names.size(); l++) ssbar_values[l].push_back(Offres_MC_values[l].at(k));
            ssbar_weights.push_back(weights.at(k));
        }
        else if (Offres_MC_numbering.at(k) == 3) {
            for (int l = 0; l < (int)variable_names.size(); l++) ccbar_values[l].push_back(Offres_MC_values[l].at(k));
            ccbar_weights.push_back(weights.at(k));
        }
        else if (Offres_MC_numbering.at(k) == 4) {
            for (int l = 0; l < (int)variable_names.size(); l++) taupair_values[l].push_back(Offres_MC_values[l].at(k));
            taupair_weights.push_back(weights.at(k));
        }
        else if (Offres_MC_numbering.at(k) == 5) {
            for (int l = 0; l < (int)variable_names.size(); l++) mumu_values[l].push_back(Offres_MC_values[l].at(k));
            mumu_weights.push_back(weights.at(k));
        }
        else if (Offres_MC_numbering.at(k) == 6) {
            for (int l = 0; l < (int)variable_names.size(); l++) gg_values[l].push_back(Offres_MC_values[l].at(k));
            gg_weights.push_back(weights.at(k));
        }
        else if (Offres_MC_numbering.at(k) == 7) {
            for (int l = 0; l < (int)variable_names.size(); l++) ee_values[l].push_back(Offres_MC_values[l].at(k));
            ee_weights.push_back(weights.at(k));
        }
        else if (Offres_MC_numbering.at(k) == 8) {
            for (int l = 0; l < (int)variable_names.size(); l++) eeee_values[l].push_back(Offres_MC_values[l].at(k));
            eeee_weights.push_back(weights.at(k));
        }
        else if (Offres_MC_numbering.at(k) == 9) {
            for (int l = 0; l < (int)variable_names.size(); l++) eemumu_values[l].push_back(Offres_MC_values[l].at(k));
            eemumu_weights.push_back(weights.at(k));
        }
        else if (Offres_MC_numbering.at(k) == 10) {
            for (int l = 0; l < (int)variable_names.size(); l++) llXX_values[l].push_back(Offres_MC_values[l].at(k));
            llXX_weights.push_back(weights.at(k));
        }
        else if (Offres_MC_numbering.at(k) == 11) {
            for (int l = 0; l < (int)variable_names.size(); l++) hhISR_values[l].push_back(Offres_MC_values[l].at(k));
            hhISR_weights.push_back(weights.at(k));
        }
        else {
            printf("undefined numbering!\n");
            exit(1);
        }
    }

    THStack* Stack[Nvar_num];
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
    TH1D* stat_error_hist[Nvar_num];
    TH1D* data_hist[Nvar_num];
    TH1D* Ratio_hist[Nvar_num];

    TH1D* MC_one_bin[Nvar_num];
    TH1D* data_one_bin[Nvar_num];
    TH1D* Ratio_one_bin[Nvar_num];

    for (int k = 0; k < (int)variable_names.size(); k++) { // malloc TH1D
        std::vector<double> temp_v;
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
        temp_v.insert(temp_v.end(), Offres_data_values[k].begin(), Offres_data_values[k].end());


        double min = *min_element(temp_v.begin(), temp_v.end());
        double max = *max_element(temp_v.begin(), temp_v.end());
        int bins = 10;
        min = 0.9;
        max = 1.0;

        if (hasEnding(variable_names.at(k),std::string("dr"))) { // exceptions
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
        else if(hasEnding(variable_names.at(k), std::string("harmonicMomentThrust1"))){
            max = 1.0;
            min = -1.0;
        }
        else if(hasEnding(variable_names.at(k), std::string("harmonicMomentThrust2"))){
            max = 1.0;
            min = 0.0;
        }
        else if(hasEnding(variable_names.at(k), std::string("harmonicMomentThrust3"))){
            max = 1.0;
            min = -1.0;
        }
        else if(hasEnding(variable_names.at(k), std::string("harmonicMomentThrust4"))){
            max = 1.0;
            min = -0.5;
        }
        else if(hasEnding(variable_names.at(k), std::string("missingMomentumOfEvent"))){
            max = 5.0;
        }
        else if(hasEnding(variable_names.at(k), std::string("missingEnergyOfEventCMS"))){
            min = -1.5;
        }
        else if(hasEnding(variable_names.at(k), std::string("Btag_extraInfo_SignalProbability"))){
            max = 0;
            min = -3;
            variable_names.at(k) = std::string("log_{10}SignalProbability");
        }
        else if(hasEnding(variable_names.at(k), std::string("Btag_thrustOm"))){
            min = 0.5;
        }

        Stack[k] = new THStack(variable_names.at(k).c_str(), (";"+ variable_names.at(k) + ";number of candidates").c_str());
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
        stat_error_hist[k] = new TH1D("MC stat error", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        data_hist[k] = new TH1D("data", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        Ratio_hist[k] = new TH1D((variable_names.at(k) + "_ratio").c_str(), ";;data/MC", bins, min, max);

        MC_one_bin[k] = new TH1D((variable_names.at(k) + "_MC_one_bin").c_str(), ";number of candidates", 1, min, max);
        data_one_bin[k] = new TH1D((variable_names.at(k) + "_data_one_bin").c_str(), ";number of candidates", 1, min, max);
        Ratio_one_bin[k] = new TH1D((variable_names.at(k) + "_ratio_one_bin").c_str(), ";number of candidates", 1, min, max);
    }

    for (int k = 0; k < (int)variable_names.size(); k++) { // fill
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
        for (int i = 0; i < (int)Offres_MC_values[k].size(); i++) stat_error_hist[k]->Fill(Offres_MC_values[k].at(i), weights.at(i));
        for (int i = 0; i < (int)Offres_data_values[k].size(); i++) data_hist[k]->Fill(Offres_data_values[k].at(i));

        for (int i = 0; i < (int)Offres_MC_values[k].size(); i++) MC_one_bin[k]->Fill(Offres_MC_values[k].at(i), weights.at(i));
        for (int i = 0; i < (int)Offres_data_values[k].size(); i++) data_one_bin[k]->Fill(Offres_data_values[k].at(i));
    }

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
    printf("data: %d\n", (int)Offres_data_values[0].size());

    for (int k = 0; k < (int)variable_names.size(); k++) { // draw
        Stack[k]->Add(uubar_hist[k]);
        Stack[k]->Add(ddbar_hist[k]);
        Stack[k]->Add(ssbar_hist[k]);
        Stack[k]->Add(ccbar_hist[k]);
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

        Ratio_one_bin[k]->Divide(data_one_bin[k], MC_one_bin[k]);

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
        TPaveText *pt = new TPaveText(0.135,0.88,0.5, 1.0,"NDC NB"); pt->SetFillStyle(0); pt->SetLineWidth(0); pt->AddText( ("MC scaled to data, Data/MC= " + std::to_string(CAL)).c_str() ); pt->Draw();

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
        c_temp->SaveAs((variable_names.at(k)+"_offres.png").c_str());

        delete c_temp;
    }

    // Print data-MC discrepancy
    double MC_sum = 0;
    for (int i = 0; i < (int)Offres_MC_values[0].size(); i++) MC_sum = MC_sum + weights.at(i);
    printf("data num: %ld\n", Offres_data_values[0].size());
    printf("MC num with calibration: %lf\n\n", MC_sum);
    printf("MC with calibration: %lf +- %lf\n", MC_one_bin[0]->GetBinContent(1), MC_one_bin[0]->GetBinError(1));
    printf("data with calibration: %lf +- %lf\n", data_one_bin[0]->GetBinContent(1), data_one_bin[0]->GetBinError(1));
    printf("data/MC with calibration: %lf +- %lf\n", Ratio_one_bin[0]->GetBinContent(1), Ratio_one_bin[0]->GetBinError(1));
}
