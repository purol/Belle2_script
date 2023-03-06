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

// BKG MC sample number (0.8/ab for BB, 1.0/ab for qq)
# define N_CHG_test 32042497.0
# define N_MIX_test 24693710.0
# define N_UUBAR_test 94447089.0
# define N_DDBAR_test 22664556.0
# define N_SSBAR_test 19244661.0
# define N_CHARM_test 107541168.0
# define N_CHG_train 48052238.0
# define N_MIX_train 37030486.0
# define N_UUBAR_train 141671998.0
# define N_DDBAR_train 34114182.0
# define N_SSBAR_train 28859338.0
# define N_CHARM_train 161280679.0

// new scale factor for BKG MC sample with additional 1/ab (364.436 - 2.763 = 361.673/fb), until LS1
# define Scale_CHG_train ((N_BB_LS1* (BR_BpBp / (BR_BpBp + BR_B0B0))) / (0.8 * N_BpBp_1invab * (N_CHG_train / (N_CHG_train + N_CHG_test)) + N_BpBp_1invab))
# define Scale_MIX_train ((N_BB_LS1* (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (0.8 * N_B0B0_1invab * (N_MIX_train / (N_MIX_train + N_MIX_test)) + N_B0B0_1invab))
# define Scale_UUBAR_train (0.361673/((N_UUBAR_train/(N_UUBAR_train + N_UUBAR_test))*1.0))
# define Scale_DDBAR_train (0.361673/((N_DDBAR_train/(N_DDBAR_train + N_DDBAR_test))*1.0))
# define Scale_SSBAR_train (0.361673/((N_SSBAR_train/(N_SSBAR_train + N_SSBAR_test))*1.0))
# define Scale_CHARM_train (0.361673/((N_CHARM_train/(N_CHARM_train + N_CHARM_test))*1.0))
# define Scale_CHG_test ((N_BB_LS1* (BR_BpBp / (BR_BpBp + BR_B0B0))) / (0.8 * N_BpBp_1invab * (N_CHG_test / (N_CHG_train + N_CHG_test)) + N_BpBp_1invab))
# define Scale_MIX_test ((N_BB_LS1* (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (0.8 * N_B0B0_1invab * (N_MIX_test / (N_MIX_train + N_MIX_test)) + N_B0B0_1invab))
# define Scale_UUBAR_test (0.361673/((N_UUBAR_test/(N_UUBAR_train + N_UUBAR_test))*1.0))
# define Scale_DDBAR_test (0.361673/((N_DDBAR_test/(N_DDBAR_train + N_DDBAR_test))*1.0))
# define Scale_SSBAR_test (0.361673/((N_SSBAR_test/(N_SSBAR_train + N_SSBAR_test))*1.0))
# define Scale_CHARM_test (0.361673/((N_CHARM_test/(N_CHARM_train + N_CHARM_test))*1.0))

# define KS0_rel_uncertainty 0.6 // %/cm
# define track_rel_uncertainty 0.69 // %
# define pi0_correction 0.932
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

    double FEI_calibration_factor = -1;

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

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

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
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(PID_correction[0][i_PID], temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(PID_correction[1][i_PID], temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(PID_correction[2][i_PID], temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(PID_correction[3][i_PID], temp_N_bin_PID[3][i_PID]); // mis PID
            }
            if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) { // B2Kc
                weights->push_back(FEI_calibration_factor *CAL * Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) { // B2KcPi0
                weights->push_back(FEI_calibration_factor * pi0_correction*CAL * Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) { // B2Ks0Pic
                weights->push_back(FEI_calibration_factor *CAL * Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) { // B2KcPicPic
                weights->push_back(FEI_calibration_factor *CAL * Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) { // B2Ks0PicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction*CAL * Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) { // B2KcPicPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction*CAL * Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) { // B2Ks0PicPicPic
                weights->push_back(FEI_calibration_factor *CAL * Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) { // B2KcPicPicPicPic
                weights->push_back(FEI_calibration_factor *CAL * Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) { // B2Ks0PicPicPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction*CAL * Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) { // B2KcPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction*CAL * Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) { // B2Ks0PicPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction*CAL * Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) { // B2KcPicPicPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction*CAL * Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) { // B2KcKcKc
                weights->push_back(FEI_calibration_factor *CAL * Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) { // B2KcKcKs0Pic
                weights->push_back(FEI_calibration_factor *CAL * Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) { // B2KcKcKcPi0
                weights->push_back(FEI_calibration_factor* pi0_correction*CAL* Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) { // B02Ks0
                weights->push_back(FEI_calibration_factor*CAL* Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) { // B02KcPic
                weights->push_back(FEI_calibration_factor*CAL* Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) { // B02Ks0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction*CAL * Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) { // B02KcPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction*CAL * Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) { // B02Ks0PicPic
                weights->push_back(FEI_calibration_factor*CAL * Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) { // B02KcPicPicPic
                weights->push_back(FEI_calibration_factor*CAL * Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) { // B02Ks0PicPicPi0
                weights->push_back(FEI_calibration_factor* pi0_correction*CAL * Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) { // B02KcPicPicPicPi0
                weights->push_back(FEI_calibration_factor* pi0_correction*CAL* Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) { // B02Ks0PicPicPicPic
                weights->push_back(FEI_calibration_factor* CAL* Stream* Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) { // B02Ks0Pi0Pi0
                weights->push_back(FEI_calibration_factor* pi0_correction* pi0_correction* CAL* Stream* Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) { // B02KcPicPi0Pi0
                weights->push_back(FEI_calibration_factor* pi0_correction* pi0_correction* CAL* Stream* Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) { // B02Ks0PicPicPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) { // B02KcKcKs0
                weights->push_back(FEI_calibration_factor*CAL* Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) { // B02KcKcKcPic
                weights->push_back(FEI_calibration_factor*CAL* Stream * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) { // B02KcKcKs0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction*CAL * Stream * Correction_KID * Correction_PID);
            }
            else {
                printf("[ERROR] unexpected decay ID\n");
                exit(1);
            }


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

    double FEI_calibration_factor = -1;

    int nXsu = -1;
    int nXsd = -1;

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
        if (SampleName == "SIGNAL") {
            tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCcomb__bc", &nXsu);
            tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCcomb__bc", &nXsd);
        }

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if (SampleName == "SIGNAL") tree_Xs->GetEntry(j);

            for (int k = 0; k < (int)variable_names.size(); k++) variable_values[k].push_back(var[k]);

            // Fill numberings
            double weight_ri = 0.0;
            if (SampleName == "CHG") {
                numberings->push_back(0);
                FEI_calibration_factor = FEI_cal_Bc;
                weight_ri = ((N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0))) / (0.8 * N_BpBp_1invab)); // total 0.8/ab for BB
            }
            else if (SampleName == "MIX") {
                numberings->push_back(1);
                FEI_calibration_factor = FEI_cal_B0;
                weight_ri = ((N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (0.8 * N_B0B0_1invab)); // total 0.8/ab for BB
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
                if (nXsu > 0) {
                    FEI_calibration_factor = FEI_cal_Bc;
                    weight_ri = ((N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0))) / (1.8 * N_BpBp_1invab)); // total 1.8/ab
                }
                else if (nXsd > 0) {
                    FEI_calibration_factor = FEI_cal_B0;
                    weight_ri = ((N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (1.8 * N_B0B0_1invab)); // total 1.8/ab
                }
                else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5) {
                    FEI_calibration_factor = FEI_cal_Bc;
                    weight_ri = ((N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0))) / (1.8 * N_BpBp_1invab)); // total 1.8/ab
                }
                else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5) {
                    FEI_calibration_factor = FEI_cal_B0;
                    weight_ri = ((N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (1.8 * N_B0B0_1invab)); // total 1.8/ab
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
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(PID_correction[0][i_PID], temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(PID_correction[1][i_PID], temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(PID_correction[2][i_PID], temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(PID_correction[3][i_PID], temp_N_bin_PID[3][i_PID]); // mis PID
            }
            if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) { // B2Kc
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) { // B2KcPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) { // B2Ks0Pic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) { // B2KcPicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) { // B2Ks0PicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) { // B2KcPicPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) { // B2Ks0PicPicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) { // B2KcPicPicPicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) { // B2Ks0PicPicPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) { // B2KcPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) { // B2Ks0PicPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) { // B2KcPicPicPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) { // B2KcKcKc
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) { // B2KcKcKs0Pic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) { // B2KcKcKcPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) { // B02Ks0
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) { // B02KcPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) { // B02Ks0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) { // B02KcPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) { // B02Ks0PicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) { // B02KcPicPicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) { // B02Ks0PicPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) { // B02KcPicPicPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) { // B02Ks0PicPicPicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) { // B02Ks0Pi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) { // B02KcPicPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) { // B02Ks0PicPicPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) { // B02KcKcKs0
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) { // B02KcKcKcPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) { // B02KcKcKs0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID);
            }
            else {
                printf("[ERROR] unexpected decay ID\n");
                exit(1);
            }


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

    double FEI_calibration_factor = -1;

    int nXsu = -1;
    int nXsd = -1;

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
        if (SampleName == "SIGNAL") {
            tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCcomb__bc", &nXsu);
            tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCcomb__bc", &nXsd);
        }
        tree_upsilon->SetBranchAddress("MVA_Continuum", &BDTc);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

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
                weight_ri = ((N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0))) / (0.8 * N_BpBp_1invab)); // total 0.8/ab for BB
            }
            else if (SampleName == "MIX") {
                numberings->push_back(1);
                FEI_calibration_factor = FEI_cal_B0;
                weight_ri = ((N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (0.8 * N_B0B0_1invab)); // total 0.8/ab for BB
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
                if (nXsu > 0) {
                    FEI_calibration_factor = FEI_cal_Bc;
                    weight_ri = ((N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0))) / (1.8 * N_BpBp_1invab)); // total 1.8/ab
                }
                else if (nXsd > 0) {
                    FEI_calibration_factor = FEI_cal_B0;
                    weight_ri = ((N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (1.8 * N_B0B0_1invab)); // total 1.8/ab
                }
                else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5) {
                    FEI_calibration_factor = FEI_cal_Bc;
                    weight_ri = ((N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0))) / (1.8 * N_BpBp_1invab)); // total 1.8/ab
                }
                else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5) {
                    FEI_calibration_factor = FEI_cal_B0;
                    weight_ri = ((N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (1.8 * N_B0B0_1invab)); // total 1.8/ab
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
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(PID_correction[0][i_PID], temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(PID_correction[1][i_PID], temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(PID_correction[2][i_PID], temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(PID_correction[3][i_PID], temp_N_bin_PID[3][i_PID]); // mis PID
            }
            if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) { // B2Kc
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) { // B2KcPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) { // B2Ks0Pic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) { // B2KcPicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) { // B2Ks0PicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) { // B2KcPicPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) { // B2Ks0PicPicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) { // B2KcPicPicPicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) { // B2Ks0PicPicPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) { // B2KcPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) { // B2Ks0PicPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) { // B2KcPicPicPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) { // B2KcKcKc
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) { // B2KcKcKs0Pic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) { // B2KcKcKcPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) { // B02Ks0
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) { // B02KcPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) { // B02Ks0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) { // B02KcPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) { // B02Ks0PicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) { // B02KcPicPicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) { // B02Ks0PicPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) { // B02KcPicPicPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) { // B02Ks0PicPicPicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) { // B02Ks0Pi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) { // B02KcPicPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) { // B02Ks0PicPicPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) { // B02KcKcKs0
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) { // B02KcKcKcPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) { // B02KcKcKs0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction);
            }
            else {
                printf("[ERROR] unexpected decay ID\n");
                exit(1);
            }


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

    double FEI_calibration_factor = -1;

    int nXsu = -1;
    int nXsd = -1;

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

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID); // charged: 0, mixed: 1
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKtruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[0][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKmisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[1][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npitruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[2][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npimisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[3][i_PID]);
        }
        if (SampleName == "SIGNAL") {
            tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCcomb__bc", &nXsu);
            tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCcomb__bc", &nXsd);
        }
        tree_upsilon->SetBranchAddress("MVA_Continuum", &BDTc);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            // BDTc correction factor
            if (BDTc > (5.0 / 6.0)) BDTc_correction = 5.0;
            else BDTc_correction = (BDTc / (1.0 - BDTc));

            // Fill numberings
            double weight_ri = 0.0;
            if (SampleName == "CHG") {
                FEI_calibration_factor = FEI_cal_Bc;
                weight_ri = ((N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0))) / (0.8 * N_BpBp_1invab)); // total 0.8/ab for BB
            }
            else if (SampleName == "MIX") {
                FEI_calibration_factor = FEI_cal_B0;
                weight_ri = ((N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (0.8 * N_B0B0_1invab)); // total 0.8/ab for BB
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
                if (nXsu > 0) {
                    FEI_calibration_factor = FEI_cal_Bc;
                    weight_ri = ((N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0))) / (1.8 * N_BpBp_1invab)); // total 1.8/ab
                }
                else if (nXsd > 0) {
                    FEI_calibration_factor = FEI_cal_B0;
                    weight_ri = ((N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (1.8 * N_B0B0_1invab)); // total 1.8/ab
                }
                else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5) {
                    FEI_calibration_factor = FEI_cal_Bc;
                    weight_ri = ((N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0))) / (1.8 * N_BpBp_1invab)); // total 1.8/ab
                }
                else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5) {
                    FEI_calibration_factor = FEI_cal_B0;
                    weight_ri = ((N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (1.8 * N_B0B0_1invab)); // total 1.8/ab
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
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(PID_correction[0][i_PID], temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(PID_correction[1][i_PID], temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(PID_correction[2][i_PID], temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(PID_correction[3][i_PID], temp_N_bin_PID[3][i_PID]); // mis PID
            }
            if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) { // B2Kc
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) { // B2KcPi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) { // B2Ks0Pic
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) { // B2KcPicPic
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) { // B2Ks0PicPi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) { // B2KcPicPicPi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) { // B2Ks0PicPicPic
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) { // B2KcPicPicPicPic
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) { // B2Ks0PicPicPicPi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) { // B2KcPi0Pi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) { // B2Ks0PicPi0Pi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) { // B2KcPicPicPi0Pi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) { // B2KcKcKc
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) { // B2KcKcKs0Pic
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) { // B2KcKcKcPi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) { // B02Ks0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) { // B02KcPic
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) { // B02Ks0Pi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) { // B02KcPicPi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) { // B02Ks0PicPic
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) { // B02KcPicPicPic
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) { // B02Ks0PicPicPi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) { // B02KcPicPicPicPi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) { // B02Ks0PicPicPicPic
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) { // B02Ks0Pi0Pi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) { // B02KcPicPi0Pi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) { // B02Ks0PicPicPi0Pi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) { // B02KcKcKs0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) { // B02KcKcKcPic
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) { // B02KcKcKs0Pi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * Correction_PID * BDTc_correction;
            }
            else {
                printf("[ERROR] unexpected decay ID\n");
                exit(1);
            }


        }
        input_file->Close();

    }

}

void THStack_Jpsi_FBDT_efficiency() {

    ReadPIDFile();

    // dirnames
    const char* Jpsi_MC_SIGNAL_before_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_Jpsi/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application";
    const char* Jpsi_MC_CHG_before_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_Jpsi/CHG_analysis/validation_v000/final_output_root_after_MVA_Application";
    const char* Jpsi_MC_MIX_before_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_Jpsi/MIX_analysis/validation_v000/final_output_root_after_MVA_Application";
    const char* Jpsi_MC_UUBAR_before_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_Jpsi/UUBAR_analysis/validation_v000/final_output_root_after_MVA_Application";
    const char* Jpsi_MC_DDBAR_before_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_Jpsi/DDBAR_analysis/validation_v000/final_output_root_after_MVA_Application";
    const char* Jpsi_MC_SSBAR_before_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_Jpsi/SSBAR_analysis/validation_v000/final_output_root_after_MVA_Application";
    const char* Jpsi_MC_CHARM_before_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_Jpsi/CHARM_analysis/validation_v000/final_output_root_after_MVA_Application";
    const char* Jpsi_data_before_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_data_Jpsi/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application";

    const char* Jpsi_MC_SIGNAL_after_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_Jpsi/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut";
    const char* Jpsi_MC_CHG_after_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_Jpsi/CHG_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut";
    const char* Jpsi_MC_MIX_after_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_Jpsi/MIX_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut";
    const char* Jpsi_MC_UUBAR_after_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_Jpsi/UUBAR_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut";
    const char* Jpsi_MC_DDBAR_after_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_Jpsi/DDBAR_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut";
    const char* Jpsi_MC_SSBAR_after_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_Jpsi/SSBAR_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut";
    const char* Jpsi_MC_CHARM_after_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_Jpsi/CHARM_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut";
    const char* Jpsi_data_after_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_data_Jpsi/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut";

    std::vector<std::string> variable_names;
    std::vector<std::string> branch_names;
    variable_names.push_back("MVA_BB"); branch_names.push_back("Upsilon");

    int Nvar = static_cast<int>(variable_names.size());
    if (Nvar != Nvar_num) exit(1);
    std::vector<double> Jpsi_MC_values_before[Nvar_num];
    std::vector<double> Jpsi_MC_values_after[Nvar_num];
    std::vector<int> Jpsi_MC_numbering_before;
    std::vector<int> Jpsi_MC_numbering_after;

    std::vector<double> Jpsi_data_values_before[Nvar_num];
    std::vector<double> Jpsi_data_values_after[Nvar_num];

    std::vector<double> weights_before;
    std::vector<double> weights_after;

    LetsFillJpsi_ri(Jpsi_MC_SIGNAL_before_dirname, variable_names, branch_names, Jpsi_MC_values_before, &Jpsi_MC_numbering_before, &weights_before, "SIGNAL");
    LetsFillJpsi_ri(Jpsi_MC_CHG_before_dirname, variable_names, branch_names, Jpsi_MC_values_before, &Jpsi_MC_numbering_before, &weights_before, "CHG");
    LetsFillJpsi_ri(Jpsi_MC_MIX_before_dirname, variable_names, branch_names, Jpsi_MC_values_before, &Jpsi_MC_numbering_before, &weights_before, "MIX");
    LetsFillJpsi_ri(Jpsi_MC_UUBAR_before_dirname, variable_names, branch_names, Jpsi_MC_values_before, &Jpsi_MC_numbering_before, &weights_before, "UUBAR");
    LetsFillJpsi_ri(Jpsi_MC_DDBAR_before_dirname, variable_names, branch_names, Jpsi_MC_values_before, &Jpsi_MC_numbering_before, &weights_before, "DDBAR");
    LetsFillJpsi_ri(Jpsi_MC_SSBAR_before_dirname, variable_names, branch_names, Jpsi_MC_values_before, &Jpsi_MC_numbering_before, &weights_before, "SSBAR");
    LetsFillJpsi_ri(Jpsi_MC_CHARM_before_dirname, variable_names, branch_names, Jpsi_MC_values_before, &Jpsi_MC_numbering_before, &weights_before, "CHARM");
    LetsFillJpsi_ri(Jpsi_MC_SIGNAL_after_dirname, variable_names, branch_names, Jpsi_MC_values_after, &Jpsi_MC_numbering_after, &weights_after, "SIGNAL");
    LetsFillJpsi_ri(Jpsi_MC_CHG_after_dirname, variable_names, branch_names, Jpsi_MC_values_after, &Jpsi_MC_numbering_after, &weights_after, "CHG");
    LetsFillJpsi_ri(Jpsi_MC_MIX_after_dirname, variable_names, branch_names, Jpsi_MC_values_after, &Jpsi_MC_numbering_after, &weights_after, "MIX");
    LetsFillJpsi_ri(Jpsi_MC_UUBAR_after_dirname, variable_names, branch_names, Jpsi_MC_values_after, &Jpsi_MC_numbering_after, &weights_after, "UUBAR");
    LetsFillJpsi_ri(Jpsi_MC_DDBAR_after_dirname, variable_names, branch_names, Jpsi_MC_values_after, &Jpsi_MC_numbering_after, &weights_after, "DDBAR");
    LetsFillJpsi_ri(Jpsi_MC_SSBAR_after_dirname, variable_names, branch_names, Jpsi_MC_values_after, &Jpsi_MC_numbering_after, &weights_after, "SSBAR");
    LetsFillJpsi_ri(Jpsi_MC_CHARM_after_dirname, variable_names, branch_names, Jpsi_MC_values_after, &Jpsi_MC_numbering_after, &weights_after, "CHARM");
    LetsFill(Jpsi_data_before_dirname, variable_names, branch_names, Jpsi_data_values_before);
    LetsFill(Jpsi_data_after_dirname, variable_names, branch_names, Jpsi_data_values_after);

    TH1D* MC_SIGNAL_before_one_bin = new TH1D("MC_SIGNAL_before_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* MC_BKG_before_one_bin = new TH1D("MC_BKG_before_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* data_before_one_bin = new TH1D("data_before_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);

    TH1D* MC_SIGNAL_after_one_bin = new TH1D("MC_SIGNAL_after_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* MC_BKG_after_one_bin = new TH1D("MC_BKG_after_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* data_after_one_bin = new TH1D("data_after_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);

    for (int k = 0; k < (int)Jpsi_MC_numbering_before.size(); k++) {
        if (Jpsi_MC_numbering_before.at(k) >= 0 && Jpsi_MC_numbering_before.at(k) <= 13) { // BKG
            for (int l = 0; l < (int)variable_names.size(); l++) MC_BKG_before_one_bin->Fill(Jpsi_MC_values_before[l].at(k), weights_before.at(k));
        }
        else if (Jpsi_MC_numbering_before.at(k) == 14) { // signal
            for (int l = 0; l < (int)variable_names.size(); l++) MC_SIGNAL_before_one_bin->Fill(Jpsi_MC_values_before[l].at(k), weights_before.at(k));
        }
        else {
            printf("undefined numbering!\n");
            exit(1);
        }
    }

    for (int k = 0; k < (int)Jpsi_data_values_before[0].size(); k++) {
        for (int l = 0; l < (int)variable_names.size(); l++) data_before_one_bin->Fill(Jpsi_data_values_before[l].at(k));
    }

    for (int k = 0; k < (int)Jpsi_MC_numbering_after.size(); k++) {
        if (Jpsi_MC_numbering_after.at(k) >= 0 && Jpsi_MC_numbering_after.at(k) <= 13) { // BKG
            for (int l = 0; l < (int)variable_names.size(); l++) MC_BKG_after_one_bin->Fill(Jpsi_MC_values_after[l].at(k), weights_after.at(k));
        }
        else if (Jpsi_MC_numbering_after.at(k) == 14) { // signal
            for (int l = 0; l < (int)variable_names.size(); l++) MC_SIGNAL_after_one_bin->Fill(Jpsi_MC_values_after[l].at(k), weights_after.at(k));
        }
        else {
            printf("undefined numbering!\n");
            exit(1);
        }
    }

    for (int k = 0; k < (int)Jpsi_data_values_after[0].size(); k++) {
        for (int l = 0; l < (int)variable_names.size(); l++) data_after_one_bin->Fill(Jpsi_data_values_after[l].at(k));
    }

    double Nsig_before = MC_SIGNAL_before_one_bin->GetBinContent(1);
    double NBKG_before = MC_BKG_before_one_bin->GetBinContent(1);
    double Ndata_before = data_before_one_bin->GetBinContent(1);
    double Nsig_after = MC_SIGNAL_after_one_bin->GetBinContent(1);
    double NBKG_after = MC_BKG_after_one_bin->GetBinContent(1);
    double Ndata_after = data_after_one_bin->GetBinContent(1);

    double Nsig_before_uncer = MC_SIGNAL_before_one_bin->GetBinError(1);
    double NBKG_before_uncer = MC_BKG_before_one_bin->GetBinError(1);
    double Ndata_before_uncer = data_before_one_bin->GetBinError(1);
    double Nsig_after_uncer = MC_SIGNAL_after_one_bin->GetBinError(1);
    double NBKG_after_uncer = MC_BKG_after_one_bin->GetBinError(1);
    double Ndata_after_uncer = data_after_one_bin->GetBinError(1);

    // Print data-MC discrepancy
    printf("Nsig before: %lf +- %lf\n", Nsig_before, Nsig_before_uncer);
    printf("NBKG before: %lf +- %lf\n", NBKG_before, NBKG_before_uncer);
    printf("Ndata before: %lf +- %lf\n", Ndata_before, Ndata_before_uncer);
    printf("Nsig after: %lf +- %lf\n", Nsig_after, Nsig_after_uncer);
    printf("NBKG after: %lf +- %lf\n", NBKG_after, NBKG_after_uncer);
    printf("Ndata after: %lf +- %lf\n", Ndata_after, Ndata_after_uncer);

    double efficiency_data = (Ndata_after - NBKG_after) / (Ndata_before - NBKG_before);
    double efficiency_MC = Nsig_after / Nsig_before;

    /* ============================== assume no correlation ============================== */
    /* calculate uncertainty of efficiency data*/
    double efficiency_data_uncer_Ndata_after = Ndata_after_uncer / (Ndata_after - NBKG_after);
    double efficiency_data_uncer_NBKG_after = NBKG_after_uncer / (Ndata_after - NBKG_after);
    double efficiency_data_uncer_Ndata_before = Ndata_before_uncer / (Ndata_before - NBKG_before);
    double efficiency_data_uncer_NBKG_before = NBKG_before_uncer / (Ndata_before - NBKG_before);
    double efficiency_data_uncer = std::sqrt(
        efficiency_data_uncer_Ndata_after * efficiency_data_uncer_Ndata_after +
        efficiency_data_uncer_NBKG_after * efficiency_data_uncer_NBKG_after +
        efficiency_data_uncer_Ndata_before * efficiency_data_uncer_Ndata_before +
        efficiency_data_uncer_NBKG_before * efficiency_data_uncer_NBKG_before
    ); // it is relative uncertainty

    /* calculate uncertainty of efficiency MC*/
    double efficiency_MC_uncer_Nsig_after = Nsig_after_uncer / Nsig_after;
    double efficiency_MC_uncer_Nsig_before = Nsig_before_uncer / Nsig_before;
    double efficiency_MC_uncer = std::sqrt(
        efficiency_MC_uncer_Nsig_after * efficiency_MC_uncer_Nsig_after +
        efficiency_MC_uncer_Nsig_before * efficiency_MC_uncer_Nsig_before
    ); // it is relative uncertainty

    /* calculate total uncertainty */
    double efficiency_ratio = efficiency_data / efficiency_MC;
    double efficiency_ratio_uncer = efficiency_ratio * std::sqrt(
        efficiency_data_uncer * efficiency_data_uncer +
        efficiency_MC_uncer * efficiency_MC_uncer
    ); // it is relative uncertainty
    /* ============================== assume no correlation ============================== */

    /* ============================== assume 100% correlation ============================== */
    /* calculate uncertainty of efficiency data*/
    //double efficiency_data_uncer_Ndata_after = Ndata_after_uncer / (Ndata_after - NBKG_after);
    //double efficiency_data_uncer_NBKG_after = NBKG_after_uncer / (Ndata_after - NBKG_after);
    //double efficiency_data_uncer_Ndata_before = Ndata_before_uncer / (Ndata_before - NBKG_before);
    //double efficiency_data_uncer_NBKG_before = NBKG_before_uncer / (Ndata_before - NBKG_before);
    //double efficiency_data_uncer = std::sqrt(
    //    (efficiency_data_uncer_Ndata_after - efficiency_data_uncer_Ndata_before) *
    //    (efficiency_data_uncer_Ndata_after - efficiency_data_uncer_Ndata_before));

    /* calculate uncertainty of efficiency MC*/
    //double efficiency_MC_uncer = std::sqrt(
    //    (-efficiency_data_uncer_NBKG_after + efficiency_data_uncer_NBKG_before) *
    //    (-efficiency_data_uncer_NBKG_after + efficiency_data_uncer_NBKG_before));

    /* calculate total uncertainty */
    //double efficiency_ratio = efficiency_data / efficiency_MC;
    //double efficiency_ratio_uncer = efficiency_ratio * std::sqrt(
    //    efficiency_data_uncer * efficiency_data_uncer +
    //    efficiency_MC_uncer * efficiency_MC_uncer
    //); // it is relative uncertainty
    /* ============================== assume 100% correlation ============================== */

    /* ============================== Uncertainty for all cut ============================== */
    /* calculate ratio between data and MC */
    double Nevt_ratio = Ndata_after / (NBKG_after + Nsig_after);

    /* calculate uncertainty of the ratio */
    double N_data_rel_uncer = Ndata_after_uncer / Ndata_after;
    double N_MC_rel_uncer = std::sqrt(Nsig_after_uncer * Nsig_after_uncer + NBKG_after_uncer * NBKG_after_uncer) / (Nsig_after + NBKG_after);

    double Nevt_ratio_uncer = std::sqrt(
        N_data_rel_uncer * N_data_rel_uncer +
        N_MC_rel_uncer * N_MC_rel_uncer
    ); // it is relative uncertainty
    /* ============================== Uncertainty for all cut ============================== */

    printf("eps_data/eps_MC for FBDT = %lf +- %lf\n", efficiency_ratio, efficiency_ratio * efficiency_ratio_uncer);
    printf("N_{data}/N_{MC}          = %lf +- %lf\n", Nevt_ratio, Nevt_ratio * Nevt_ratio_uncer);
}
