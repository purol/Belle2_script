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
# include <cmath>

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
# define N_CHG_test 135437652.0
# define N_MIX_test 104176964.0
# define N_UUBAR_test 106465323.0
# define N_DDBAR_test 25669347.0
# define N_SSBAR_test 20666484.0
# define N_CHARM_test 116381155.0
# define N_CHG_train 156250364.0
# define N_MIX_train 120553880.0
# define N_UUBAR_train 159685667.0
# define N_DDBAR_train 38563469.0
# define N_SSBAR_train 31004176.0
# define N_CHARM_train 174628991.0

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

# define Nvar_num 1

# define CAL 1.1728
# define CAL_qq 1.0
# define Stream 0.25

# define N_PID_syst 73
double PID_correction[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
double PID_correction_stat_uncer[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
double PID_correction_sys_uncer[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
double PID_correction_uncer[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss

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

void LetsFillJpsi(const char* dirname, std::vector<std::string> variable_names, std::vector<std::string> branch_names, std::vector<double> variable_values[Nvar_num], std::vector<int>* numberings, std::vector<double>* weights) {
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

    float var[Nvar_num] = { 0.0 };
    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };

    double FEI_calibration_factor = -1;

    double mu0_weight = -1;
    double mu0_mispi_weight = -1;
    double mu0_misK_weight = -1;
    double mu0_pid = -1;
    double mu1_weight = -1;
    double mu1_mispi_weight = -1;
    double mu1_misK_weight = -1;
    double mu1_pid = -1;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");
        TTree* tree_lepton_syst = (TTree*)input_file->Get("lepton_syst");

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
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_0_weight_muonID_noSVD_eff_FixedThresh05__bc", &mu0_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_0_weight_muonID_noSVD_misid_pi_FixedThresh05__bc", &mu0_mispi_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_0_weight_muonID_noSVD_misid_K_FixedThresh05__bc", &mu0_misK_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_0_mcPDG__bc", &mu0_pid);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_1_weight_muonID_noSVD_eff_FixedThresh05__bc", &mu1_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_1_weight_muonID_noSVD_misid_pi_FixedThresh05__bc", &mu1_mispi_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_1_weight_muonID_noSVD_misid_K_FixedThresh05__bc", &mu1_misK_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_1_mcPDG__bc", &mu1_pid);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            tree_lepton_syst->GetEntry(j);

            for (int k = 0; k < (int)variable_names.size(); k++) variable_values[k].push_back(var[k]);

            // Fill numberings
            int job_id = stoi(names.at(i).substr(16, 9));
            if (job_id >= 265736574 && job_id <= 265736629) {
                numberings->push_back(0);
                FEI_calibration_factor = FEI_cal_Bc;
            }
            else if (job_id >= 265736630 && job_id <= 265736675) {
                numberings->push_back(1);
                FEI_calibration_factor = FEI_cal_B0;
            }
            else if (job_id >= 265736722 && job_id <= 265736767) {
                numberings->push_back(2);
                FEI_calibration_factor = CAL_qq;
            }
            else if (job_id >= 265736676 && job_id <= 265736721) {
                numberings->push_back(3);
                FEI_calibration_factor = CAL_qq;
            }
            else if (job_id >= 256846090 && job_id <= 256846478) {
                numberings->push_back(4);
                FEI_calibration_factor = CAL_qq;
            }
            else if (job_id >= 256846479 && job_id <= 256846857) {
                numberings->push_back(5);
                FEI_calibration_factor = CAL_qq;
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
            double Correction_mu0 = 1;
            double Correction_mu1 = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(PID_correction[0][i_PID], temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(PID_correction[1][i_PID], temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(PID_correction[2][i_PID], temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(PID_correction[3][i_PID], temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(pi0_correction[i_pi0], temp_N_bin_pi0[i_pi0]);
            if (std::abs(mu0_pid) > 12.5 && std::abs(mu0_pid) < 13.5) Correction_mu0 = mu0_weight;
            else if (std::abs(mu0_pid) > 210.5 && std::abs(mu0_pid) < 211.5) Correction_mu0 = mu0_mispi_weight;
            else if (std::abs(mu0_pid) > 320.5 && std::abs(mu0_pid) < 321.5) Correction_mu0 = mu0_misK_weight;
            if (std::abs(mu1_pid) > 12.5 && std::abs(mu1_pid) < 13.5) Correction_mu1 = mu1_weight;
            else if (std::abs(mu1_pid) > 210.5 && std::abs(mu1_pid) < 211.5) Correction_mu1 = mu1_mispi_weight;
            else if (std::abs(mu1_pid) > 320.5 && std::abs(mu1_pid) < 321.5) Correction_mu1 = mu1_misK_weight;

            weights->push_back(FEI_calibration_factor * CAL * Stream * Correction_pi0 * Correction_KID * Correction_PID * Correction_mu0 * Correction_mu1);


        }
        input_file->Close();

    }

}

void LetsFillJpsi_ri(const char* dirname, std::vector<std::string> variable_names, std::vector<std::string> branch_names, std::vector<double> variable_values[Nvar_num], std::vector<int>* numberings, std::vector<double>* weights, std::string SampleName) {
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

    float var[Nvar_num] = { 0.0 };
    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };

    double FEI_calibration_factor = -1;

    double mu0_weight = -1;
    double mu0_mispi_weight = -1;
    double mu0_misK_weight = -1;
    double mu0_pid = -1;
    double mu1_weight = -1;
    double mu1_mispi_weight = -1;
    double mu1_misK_weight = -1;
    double mu1_pid = -1;

    int nBp = -1;
    int nB0 = -1;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");
        TTree* tree_Xs;
        if(SampleName == "SIGNAL") tree_Xs = (TTree*)input_file->Get("Xs");
        TTree* tree_lepton_syst = (TTree*)input_file->Get("lepton_syst");

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
        if (SampleName == "SIGNAL") {
            tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clPrimaryMC__bc", &nBp);
            tree_Xs->SetBranchAddress("nParticlesInList__boB0__clPrimaryMC__bc", &nB0);
        }
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_0_weight_muonID_noSVD_eff_FixedThresh05__bc", &mu0_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_0_weight_muonID_noSVD_misid_pi_FixedThresh05__bc", &mu0_mispi_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_0_weight_muonID_noSVD_misid_K_FixedThresh05__bc", &mu0_misK_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_0_mcPDG__bc", &mu0_pid);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_1_weight_muonID_noSVD_eff_FixedThresh05__bc", &mu1_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_1_weight_muonID_noSVD_misid_pi_FixedThresh05__bc", &mu1_mispi_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_1_weight_muonID_noSVD_misid_K_FixedThresh05__bc", &mu1_misK_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_1_mcPDG__bc", &mu1_pid);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if (SampleName == "SIGNAL") tree_Xs->GetEntry(j);
            tree_lepton_syst->GetEntry(j);

            for (int k = 0; k < (int)variable_names.size(); k++) variable_values[k].push_back(var[k]);

            // Fill numberings
            double weight_ri = 0.0;
            if (SampleName == "CHG") {
                numberings->push_back(0);
                FEI_calibration_factor = FEI_cal_Bc;
                weight_ri = ((N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab)); // total 2.8/ab for BB
            }
            else if (SampleName == "MIX") {
                numberings->push_back(1);
                FEI_calibration_factor = FEI_cal_B0;
                weight_ri = ((N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab)); // total 2.8/ab for BB
            }
            else if (SampleName == "UUBAR") {
                numberings->push_back(2);
                FEI_calibration_factor = CAL_qq;
                weight_ri = ((0.364436 - 0.002763) / 1.0); // total 1.0/ab for qq
            }
            else if (SampleName == "DDBAR") {
                numberings->push_back(3);
                FEI_calibration_factor = CAL_qq;
                weight_ri = ((0.364436 - 0.002763) / 1.0); // total 1.0/ab for qq
            }
            else if (SampleName == "SSBAR") {
                numberings->push_back(4);
                FEI_calibration_factor = CAL_qq;
                weight_ri = ((0.364436 - 0.002763) / 1.0); // total 1.0/ab for qq
            }
            else if (SampleName == "CHARM") {
                numberings->push_back(5);
                FEI_calibration_factor = CAL_qq;
                weight_ri = ((0.364436 - 0.002763) / 1.0); // total 1.0/ab for qq
            }
            //else if (job_id >= 256846858 && job_id <= 256847295) numberings->push_back(6);
            //else if (job_id >= 256847296 && job_id <= 256847807) numberings->push_back(7);
            //else if (job_id >= 256847808 && job_id <= 256848291) numberings->push_back(8);
            //else if (job_id >= 256848292 && job_id <= 256848743) numberings->push_back(9);
            //else if (job_id >= 256848744 && job_id <= 256849128) numberings->push_back(10);
            //else if (job_id >= 256849129 && job_id <= 256849396) numberings->push_back(11);
            else if (SampleName == "SIGNAL") {
                numberings->push_back(14);
                FEI_calibration_factor = CAL_qq;
                if (nBp > 0) {
                    FEI_calibration_factor = FEI_cal_Bc;
                    weight_ri = ((N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab)); // total 1.8/ab
                }
                else if (nB0 > 0) {
                    FEI_calibration_factor = FEI_cal_B0;
                    weight_ri = ((N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab)); // total 1.8/ab
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
            double Correction_mu0 = 1;
            double Correction_mu1 = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(PID_correction[0][i_PID], temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(PID_correction[1][i_PID], temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(PID_correction[2][i_PID], temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(PID_correction[3][i_PID], temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(pi0_correction[i_pi0], temp_N_bin_pi0[i_pi0]);
            if (std::abs(mu0_pid) > 12.5 && std::abs(mu0_pid) < 13.5) Correction_mu0 = mu0_weight;
            else if (std::abs(mu0_pid) > 210.5 && std::abs(mu0_pid) < 211.5) Correction_mu0 = mu0_mispi_weight;
            else if (std::abs(mu0_pid) > 320.5 && std::abs(mu0_pid) < 321.5) Correction_mu0 = mu0_misK_weight;
            if (std::abs(mu1_pid) > 12.5 && std::abs(mu1_pid) < 13.5) Correction_mu1 = mu1_weight;
            else if (std::abs(mu1_pid) > 210.5 && std::abs(mu1_pid) < 211.5) Correction_mu1 = mu1_mispi_weight;
            else if (std::abs(mu1_pid) > 320.5 && std::abs(mu1_pid) < 321.5) Correction_mu1 = mu1_misK_weight;

            weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID* Correction_mu0* Correction_mu1);


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

void LetsFillJpsi_ri_correction(const char* dirname, std::vector<std::string> variable_names, std::vector<std::string> branch_names, std::vector<double> variable_values[Nvar_num], std::vector<int>* numberings, std::vector<double>* weights, std::string SampleName, double NormFactor = 1.0) {
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

    float var[Nvar_num] = { 0.0 };
    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };

    double FEI_calibration_factor = -1;

    double mu0_weight = -1;
    double mu0_mispi_weight = -1;
    double mu0_misK_weight = -1;
    double mu0_pid = -1;
    double mu1_weight = -1;
    double mu1_mispi_weight = -1;
    double mu1_misK_weight = -1;
    double mu1_pid = -1;

    int nBp = -1;
    int nB0 = -1;

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
        TTree* tree_Xs;
        if (SampleName == "SIGNAL") tree_Xs = (TTree*)input_file->Get("Xs");
        TTree* tree_lepton_syst = (TTree*)input_file->Get("lepton_syst");

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
        if (SampleName == "SIGNAL") {
            tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clPrimaryMC__bc", &nBp);
            tree_Xs->SetBranchAddress("nParticlesInList__boB0__clPrimaryMC__bc", &nB0);
        }
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_0_weight_muonID_noSVD_eff_FixedThresh05__bc", &mu0_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_0_weight_muonID_noSVD_misid_pi_FixedThresh05__bc", &mu0_mispi_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_0_weight_muonID_noSVD_misid_K_FixedThresh05__bc", &mu0_misK_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_0_mcPDG__bc", &mu0_pid);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_1_weight_muonID_noSVD_eff_FixedThresh05__bc", &mu1_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_1_weight_muonID_noSVD_misid_pi_FixedThresh05__bc", &mu1_mispi_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_1_weight_muonID_noSVD_misid_K_FixedThresh05__bc", &mu1_misK_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_1_mcPDG__bc", &mu1_pid);
        tree_upsilon->SetBranchAddress("MVA_Continuum", &BDTc);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            tree_lepton_syst->GetEntry(j);

            // BDTc correction factor
            if (BDTc > (5.0 / 6.0)) BDTc_correction = 5.0;
            else BDTc_correction = (BDTc / (1.0 - BDTc));
            BDTc_correction = BDTc_correction / NormFactor;

            for (int k = 0; k < (int)variable_names.size(); k++) variable_values[k].push_back(var[k]);

            // Fill numberings
            double weight_ri = 0.0;
            if (SampleName == "CHG") {
                numberings->push_back(0);
                FEI_calibration_factor = FEI_cal_Bc;
                weight_ri = ((N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab)); // total 2.8/ab for BB
            }
            else if (SampleName == "MIX") {
                numberings->push_back(1);
                FEI_calibration_factor = FEI_cal_B0;
                weight_ri = ((N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab)); // total 2.8/ab for BB
            }
            else if (SampleName == "UUBAR") {
                numberings->push_back(2);
                FEI_calibration_factor = CAL_qq;
                weight_ri = ((0.364436 - 0.002763) / 1.0); // total 1.0/ab for qq
            }
            else if (SampleName == "DDBAR") {
                numberings->push_back(3);
                FEI_calibration_factor = CAL_qq;
                weight_ri = ((0.364436 - 0.002763) / 1.0); // total 1.0/ab for qq
            }
            else if (SampleName == "SSBAR") {
                numberings->push_back(4);
                FEI_calibration_factor = CAL_qq;
                weight_ri = ((0.364436 - 0.002763) / 1.0); // total 1.0/ab for qq
            }
            else if (SampleName == "CHARM") {
                numberings->push_back(5);
                FEI_calibration_factor = CAL_qq;
                weight_ri = ((0.364436 - 0.002763) / 1.0); // total 1.0/ab for qq
            }
            //else if (job_id >= 256846858 && job_id <= 256847295) numberings->push_back(6);
            //else if (job_id >= 256847296 && job_id <= 256847807) numberings->push_back(7);
            //else if (job_id >= 256847808 && job_id <= 256848291) numberings->push_back(8);
            //else if (job_id >= 256848292 && job_id <= 256848743) numberings->push_back(9);
            //else if (job_id >= 256848744 && job_id <= 256849128) numberings->push_back(10);
            //else if (job_id >= 256849129 && job_id <= 256849396) numberings->push_back(11);
            else if (SampleName == "SIGNAL") {
                numberings->push_back(14);
                FEI_calibration_factor = CAL_qq;
                if (nBp > 0) {
                    FEI_calibration_factor = FEI_cal_Bc;
                    weight_ri = ((N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab)); // total 1.8/ab
                }
                else if (nB0 > 0) {
                    FEI_calibration_factor = FEI_cal_B0;
                    weight_ri = ((N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab)); // total 1.8/ab
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
            double Correction_mu0 = 1;
            double Correction_mu1 = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(PID_correction[0][i_PID], temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(PID_correction[1][i_PID], temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(PID_correction[2][i_PID], temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(PID_correction[3][i_PID], temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(pi0_correction[i_pi0], temp_N_bin_pi0[i_pi0]);
            if (std::abs(mu0_pid) > 12.5 && std::abs(mu0_pid) < 13.5) Correction_mu0 = mu0_weight;
            else if (std::abs(mu0_pid) > 210.5 && std::abs(mu0_pid) < 211.5) Correction_mu0 = mu0_mispi_weight;
            else if (std::abs(mu0_pid) > 320.5 && std::abs(mu0_pid) < 321.5) Correction_mu0 = mu0_misK_weight;
            if (std::abs(mu1_pid) > 12.5 && std::abs(mu1_pid) < 13.5) Correction_mu1 = mu1_weight;
            else if (std::abs(mu1_pid) > 210.5 && std::abs(mu1_pid) < 211.5) Correction_mu1 = mu1_mispi_weight;
            else if (std::abs(mu1_pid) > 320.5 && std::abs(mu1_pid) < 321.5) Correction_mu1 = mu1_misK_weight;

            weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * BDTc_correction);


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
    14: signal (B-->Xs J/psi)
    */

    float var[Nvar_num] = { 0.0 };
    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };

    double FEI_calibration_factor = -1;

    double mu0_weight = -1;
    double mu0_mispi_weight = -1;
    double mu0_misK_weight = -1;
    double mu0_pid = -1;
    double mu1_weight = -1;
    double mu1_mispi_weight = -1;
    double mu1_misK_weight = -1;
    double mu1_pid = -1;

    int nBp = -1;
    int nB0 = -1;

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
        TTree* tree_Xs;
        if (SampleName == "SIGNAL") tree_Xs = (TTree*)input_file->Get("Xs");
        TTree* tree_lepton_syst = (TTree*)input_file->Get("lepton_syst");

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID); // charged: 0, mixed: 1
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKtruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[0][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKmisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[1][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npitruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[2][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npimisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[3][i_PID]);
        }
        for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npi0bin" + std::to_string(i_pi0)).c_str(), &temp_N_bin_pi0[i_pi0]);
        if (SampleName == "SIGNAL") {
            tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clPrimaryMC__bc", &nBp);
            tree_Xs->SetBranchAddress("nParticlesInList__boB0__clPrimaryMC__bc", &nB0);
        }
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_0_weight_muonID_noSVD_eff_FixedThresh05__bc", &mu0_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_0_weight_muonID_noSVD_misid_pi_FixedThresh05__bc", &mu0_mispi_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_0_weight_muonID_noSVD_misid_K_FixedThresh05__bc", &mu0_misK_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_0_mcPDG__bc", &mu0_pid);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_1_weight_muonID_noSVD_eff_FixedThresh05__bc", &mu1_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_1_weight_muonID_noSVD_misid_pi_FixedThresh05__bc", &mu1_mispi_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_1_weight_muonID_noSVD_misid_K_FixedThresh05__bc", &mu1_misK_weight);
        tree_lepton_syst->SetBranchAddress("eventExtraInfo__bomu_1_mcPDG__bc", &mu1_pid);
        tree_upsilon->SetBranchAddress("MVA_Continuum", &BDTc);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            tree_lepton_syst->GetEntry(j);

            // BDTc correction factor
            if (BDTc > (5.0 / 6.0)) BDTc_correction = 5.0;
            else BDTc_correction = (BDTc / (1.0 - BDTc));

            // Fill numberings
            double weight_ri = 0.0;
            if (SampleName == "CHG") {
                FEI_calibration_factor = FEI_cal_Bc;
                weight_ri = ((N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab)); // total 0.8/ab for BB
            }
            else if (SampleName == "MIX") {
                FEI_calibration_factor = FEI_cal_B0;
                weight_ri = ((N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab)); // total 0.8/ab for BB
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
            else if (SampleName == "SIGNAL") {
                FEI_calibration_factor = CAL_qq;
                if (nBp > 0) {
                    FEI_calibration_factor = FEI_cal_Bc;
                    weight_ri = ((N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab)); // total 1.8/ab
                }
                else if (nB0 > 0) {
                    FEI_calibration_factor = FEI_cal_B0;
                    weight_ri = ((N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab)); // total 1.8/ab
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
            double Correction_mu0 = 1;
            double Correction_mu1 = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(PID_correction[0][i_PID], temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(PID_correction[1][i_PID], temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(PID_correction[2][i_PID], temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(PID_correction[3][i_PID], temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(pi0_correction[i_pi0], temp_N_bin_pi0[i_pi0]);
            if (std::abs(mu0_pid) > 12.5 && std::abs(mu0_pid) < 13.5) Correction_mu0 = mu0_weight;
            else if (std::abs(mu0_pid) > 210.5 && std::abs(mu0_pid) < 211.5) Correction_mu0 = mu0_mispi_weight;
            else if (std::abs(mu0_pid) > 320.5 && std::abs(mu0_pid) < 321.5) Correction_mu0 = mu0_misK_weight;
            if (std::abs(mu1_pid) > 12.5 && std::abs(mu1_pid) < 13.5) Correction_mu1 = mu1_weight;
            else if (std::abs(mu1_pid) > 210.5 && std::abs(mu1_pid) < 211.5) Correction_mu1 = mu1_mispi_weight;
            else if (std::abs(mu1_pid) > 320.5 && std::abs(mu1_pid) < 321.5) Correction_mu1 = mu1_misK_weight;

            nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_mu0 * Correction_mu1;
            nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_mu0 * Correction_mu1 * BDTc_correction;


        }
        input_file->Close();

    }

}

void THStack_plot_Jpsi_FBDT() {

    ReadPIDFile();

    // dirnames
    const char* Jpsi_MC_SIGNAL_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_Jpsi/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut";
    const char* Jpsi_MC_CHG_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_Jpsi/CHG_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut";
    const char* Jpsi_MC_MIX_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_Jpsi/MIX_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut";
    const char* Jpsi_MC_UUBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_Jpsi/UUBAR_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut";
    const char* Jpsi_MC_DDBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_Jpsi/DDBAR_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut";
    const char* Jpsi_MC_SSBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_Jpsi/SSBAR_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut";
    const char* Jpsi_MC_CHARM_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_Jpsi/CHARM_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut";
    const char* Jpsi_data_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_data_Jpsi/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut";

    std::vector<std::string> variable_names;
    std::vector<std::string> branch_names;
    variable_names.push_back("MVA_BB"); branch_names.push_back("Upsilon");

    int Nvar = static_cast<int>(variable_names.size());
    if (Nvar != Nvar_num) exit(1);
    std::vector<double> Jpsi_MC_values[Nvar_num];
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
    std::vector<double> signal_values[Nvar_num];
    std::vector<int> Jpsi_MC_numbering;

    std::vector<double> Jpsi_data_values[Nvar_num];

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
    std::vector<double> signal_weights;

    LetsFillJpsi_ri(Jpsi_MC_SIGNAL_dirname, variable_names, branch_names, Jpsi_MC_values, &Jpsi_MC_numbering, &weights, "SIGNAL");
    LetsFillJpsi_ri(Jpsi_MC_CHG_dirname, variable_names, branch_names, Jpsi_MC_values, &Jpsi_MC_numbering, &weights, "CHG");
    LetsFillJpsi_ri(Jpsi_MC_MIX_dirname, variable_names, branch_names, Jpsi_MC_values, &Jpsi_MC_numbering, &weights, "MIX");
    LetsFillJpsi_ri(Jpsi_MC_UUBAR_dirname, variable_names, branch_names, Jpsi_MC_values, &Jpsi_MC_numbering, &weights, "UUBAR");
    LetsFillJpsi_ri(Jpsi_MC_DDBAR_dirname, variable_names, branch_names, Jpsi_MC_values, &Jpsi_MC_numbering, &weights, "DDBAR");
    LetsFillJpsi_ri(Jpsi_MC_SSBAR_dirname, variable_names, branch_names, Jpsi_MC_values, &Jpsi_MC_numbering, &weights, "SSBAR");
    LetsFillJpsi_ri(Jpsi_MC_CHARM_dirname, variable_names, branch_names, Jpsi_MC_values, &Jpsi_MC_numbering, &weights, "CHARM");
    LetsFill(Jpsi_data_dirname, variable_names, branch_names, Jpsi_data_values);

    // sort variables
    for (int k = 0; k < (int)Jpsi_MC_numbering.size(); k++) {
        if (Jpsi_MC_numbering.at(k) == 0) {
            for (int l = 0; l < (int)variable_names.size(); l++) charged_values[l].push_back(Jpsi_MC_values[l].at(k));
            charged_weights.push_back(weights.at(k));
        }
        else if (Jpsi_MC_numbering.at(k) == 1) {
            for (int l = 0; l < (int)variable_names.size(); l++) mixed_values[l].push_back(Jpsi_MC_values[l].at(k));
            mixed_weights.push_back(weights.at(k));
        }
        else if (Jpsi_MC_numbering.at(k) == 2) {
            for (int l = 0; l < (int)variable_names.size(); l++) uubar_values[l].push_back(Jpsi_MC_values[l].at(k));
            uubar_weights.push_back(weights.at(k));
        }
        else if (Jpsi_MC_numbering.at(k) == 3) {
            for (int l = 0; l < (int)variable_names.size(); l++) ddbar_values[l].push_back(Jpsi_MC_values[l].at(k));
            ddbar_weights.push_back(weights.at(k));
        }
        else if (Jpsi_MC_numbering.at(k) == 4) {
            for (int l = 0; l < (int)variable_names.size(); l++) ssbar_values[l].push_back(Jpsi_MC_values[l].at(k));
            ssbar_weights.push_back(weights.at(k));
        }
        else if (Jpsi_MC_numbering.at(k) == 5) {
            for (int l = 0; l < (int)variable_names.size(); l++) ccbar_values[l].push_back(Jpsi_MC_values[l].at(k));
            ccbar_weights.push_back(weights.at(k));
        }
        else if (Jpsi_MC_numbering.at(k) == 6) {
            for (int l = 0; l < (int)variable_names.size(); l++) taupair_values[l].push_back(Jpsi_MC_values[l].at(k));
            taupair_weights.push_back(weights.at(k));
        }
        else if (Jpsi_MC_numbering.at(k) == 7) {
            for (int l = 0; l < (int)variable_names.size(); l++) mumu_values[l].push_back(Jpsi_MC_values[l].at(k));
            mumu_weights.push_back(weights.at(k));
        }
        else if (Jpsi_MC_numbering.at(k) == 8) {
            for (int l = 0; l < (int)variable_names.size(); l++) gg_values[l].push_back(Jpsi_MC_values[l].at(k));
            gg_weights.push_back(weights.at(k));
        }
        else if (Jpsi_MC_numbering.at(k) == 9) {
            for (int l = 0; l < (int)variable_names.size(); l++) ee_values[l].push_back(Jpsi_MC_values[l].at(k));
            ee_weights.push_back(weights.at(k));
        }
        else if (Jpsi_MC_numbering.at(k) == 10) {
            for (int l = 0; l < (int)variable_names.size(); l++) eeee_values[l].push_back(Jpsi_MC_values[l].at(k));
            eeee_weights.push_back(weights.at(k));
        }
        else if (Jpsi_MC_numbering.at(k) == 11) {
            for (int l = 0; l < (int)variable_names.size(); l++) eemumu_values[l].push_back(Jpsi_MC_values[l].at(k));
            eemumu_weights.push_back(weights.at(k));
        }
        else if (Jpsi_MC_numbering.at(k) == 12) {
            for (int l = 0; l < (int)variable_names.size(); l++) llXX_values[l].push_back(Jpsi_MC_values[l].at(k));
            llXX_weights.push_back(weights.at(k));
        }
        else if (Jpsi_MC_numbering.at(k) == 13) {
            for (int l = 0; l < (int)variable_names.size(); l++) hhISR_values[l].push_back(Jpsi_MC_values[l].at(k));
            hhISR_weights.push_back(weights.at(k));
        }
        else if (Jpsi_MC_numbering.at(k) == 14) {
            for (int l = 0; l < (int)variable_names.size(); l++) signal_values[l].push_back(Jpsi_MC_values[l].at(k));
            signal_weights.push_back(weights.at(k));
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
    TH1D* signal_hist[Nvar_num];
    TH1D* stat_error_hist[Nvar_num];
    TH1D* data_hist[Nvar_num];
    TH1D* Ratio_hist[Nvar_num];

    TH1D* MC_one_bin[Nvar_num];
    TH1D* data_one_bin[Nvar_num];
    TH1D* Ratio_one_bin[Nvar_num];

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
        temp_v.insert(temp_v.end(), signal_values[k].begin(), signal_values[k].end());
        temp_v.insert(temp_v.end(), Jpsi_data_values[k].begin(), Jpsi_data_values[k].end());


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
        signal_hist[k] = new TH1D("signal", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        stat_error_hist[k] = new TH1D("MC stat error", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        data_hist[k] = new TH1D("data", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        Ratio_hist[k] = new TH1D((variable_names.at(k) + "_ratio").c_str(), ";;data/MC", bins, min, max);

        MC_one_bin[k] = new TH1D((variable_names.at(k) + "_MC_one_bin").c_str(), ";number of candidates", 1, min, max);
        data_one_bin[k] = new TH1D((variable_names.at(k) + "_data_one_bin").c_str(), ";number of candidates", 1, min, max);
        Ratio_one_bin[k] = new TH1D((variable_names.at(k) + "_ratio_one_bin").c_str(), ";number of candidates", 1, min, max);
    }
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
        for (int i = 0; i < (int)signal_values[k].size(); i++) signal_hist[k]->Fill(signal_values[k].at(i), signal_weights.at(i));
        for (int i = 0; i < (int)Jpsi_MC_values[k].size(); i++) stat_error_hist[k]->Fill(Jpsi_MC_values[k].at(i), weights.at(i));
        for (int i = 0; i < (int)Jpsi_data_values[k].size(); i++) data_hist[k]->Fill(Jpsi_data_values[k].at(i));

        for (int i = 0; i < (int)Jpsi_MC_values[k].size(); i++) MC_one_bin[k]->Fill(Jpsi_MC_values[k].at(i), weights.at(i));
        for (int i = 0; i < (int)Jpsi_data_values[k].size(); i++) data_one_bin[k]->Fill(Jpsi_data_values[k].at(i));
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
    printf("signal: %d\n", (int)signal_values[0].size());
    printf("data: %d\n", (int)Jpsi_data_values[0].size());

    for (int k = 0; k < (int)variable_names.size(); k++) { // draw
        Stack[k]->Add(charged_hist[k]);
        Stack[k]->Add(mixed_hist[k]);
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
        Stack[k]->Add(signal_hist[k]);

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
        c_temp->SaveAs((variable_names.at(k) + "_Jpsi.png").c_str());

        delete c_temp;
    }

    // Print data-MC discrepancy
    double MC_sum = 0;
    for (int i = 0; i < (int)Jpsi_MC_values[0].size(); i++) MC_sum = MC_sum + weights.at(i);
    printf("data num: %ld\n", Jpsi_data_values[0].size());
    printf("MC num with calibration: %lf\n", MC_sum);
    printf("MC with calibration: %lf +- %lf\n", MC_one_bin[0]->GetBinContent(1), MC_one_bin[0]->GetBinError(1));
    printf("data with calibration: %lf +- %lf\n", data_one_bin[0]->GetBinContent(1), data_one_bin[0]->GetBinError(1));
    printf("data/MC with calibration: %lf +- %lf\n", Ratio_one_bin[0]->GetBinContent(1), Ratio_one_bin[0]->GetBinError(1));
}
