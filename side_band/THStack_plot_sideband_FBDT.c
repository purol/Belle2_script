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

void LetsFillSideBand(const char* dirname, std::vector<std::string> variable_names, std::vector<std::string> branch_names, std::vector<double> variable_values[Nvar_num], std::vector<int>* numberings, std::vector<double>* weights) {
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

void LetsFillSideBand_ri(const char* dirname, std::vector<std::string> variable_names, std::vector<std::string> branch_names, std::vector<double> variable_values[Nvar_num], std::vector<int>* numberings, std::vector<double>* weights, std::string SampleName) {
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
    double var[Nvar_num] = { 0.0 };

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

void LetsFillSideBand_ri_correction(const char* dirname, std::vector<std::string> variable_names, std::vector<std::string> branch_names, std::vector<double> variable_values[Nvar_num], std::vector<int>* numberings, std::vector<double>* weights, std::string SampleName, double NormFactor = 1.0) {
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
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss

    double FEI_calibration_factor = -1;

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
    */

    double var[Nvar_num] = { 0.0 };
    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss

    double FEI_calibration_factor = -1;

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

void THStack_plot_sideband_FBDT() {

    ReadPIDFile();

    Nevt nevt_CHG = { 0.0, 0.0 };
    Nevt nevt_MIX = { 0.0, 0.0 };
    Nevt nevt_UUBAR = { 0.0, 0.0 };
    Nevt nevt_DDBAR = { 0.0, 0.0 };
    Nevt nevt_SSBAR = { 0.0, 0.0 };
    Nevt nevt_CHARM = { 0.0, 0.0 };

    const char* Sideband_MC_CHG_train_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_side/CHG_analysis/train_v000/final_output";
    const char* Sideband_MC_MIX_train_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_side/MIX_analysis/train_v000/final_output";
    const char* Sideband_MC_UUBAR_train_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_side/UUBAR_analysis/train_v000/final_output";
    const char* Sideband_MC_DDBAR_train_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_side/DDBAR_analysis/train_v000/final_output";
    const char* Sideband_MC_SSBAR_train_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_side/SSBAR_analysis/train_v000/final_output";
    const char* Sideband_MC_CHARM_train_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_side/CHARM_analysis/train_v000/final_output";

    const char* Sideband_MC_CHG_test_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_side/CHG_analysis/test_v000/final_output";
    const char* Sideband_MC_MIX_test_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_side/MIX_analysis/test_v000/final_output";
    const char* Sideband_MC_UUBAR_test_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_side/UUBAR_analysis/test_v000/final_output";
    const char* Sideband_MC_DDBAR_test_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_side/DDBAR_analysis/test_v000/final_output";
    const char* Sideband_MC_SSBAR_test_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_side/SSBAR_analysis/test_v000/final_output";
    const char* Sideband_MC_CHARM_test_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_MC_side/CHARM_analysis/test_v000/final_output";

    const char* Sideband_data_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_LS_data_side/SIGNAL_analysis/validation_v000/final_output";

    NevtCount_ri(Sideband_MC_CHG_train_dirname, "CHG", &nevt_CHG);
    NevtCount_ri(Sideband_MC_MIX_train_dirname, "MIX", &nevt_MIX);
    NevtCount_ri(Sideband_MC_UUBAR_train_dirname, "UUBAR", &nevt_UUBAR);
    NevtCount_ri(Sideband_MC_DDBAR_train_dirname, "DDBAR", &nevt_DDBAR);
    NevtCount_ri(Sideband_MC_SSBAR_train_dirname, "SSBAR", &nevt_SSBAR);
    NevtCount_ri(Sideband_MC_CHARM_train_dirname, "CHARM", &nevt_CHARM);

    NevtCount_ri(Sideband_MC_CHG_test_dirname, "CHG", &nevt_CHG);
    NevtCount_ri(Sideband_MC_MIX_test_dirname, "MIX", &nevt_MIX);
    NevtCount_ri(Sideband_MC_UUBAR_test_dirname, "UUBAR", &nevt_UUBAR);
    NevtCount_ri(Sideband_MC_DDBAR_test_dirname, "DDBAR", &nevt_DDBAR);
    NevtCount_ri(Sideband_MC_SSBAR_test_dirname, "SSBAR", &nevt_SSBAR);
    NevtCount_ri(Sideband_MC_CHARM_test_dirname, "CHARM", &nevt_CHARM);

    double NormFactor_CHG = nevt_CHG.NevtwithoutCorrection / nevt_CHG.NevtwithCorrection;
    double NormFactor_MIX = nevt_MIX.NevtwithoutCorrection / nevt_MIX.NevtwithCorrection;
    double NormFactor_UUBAR = nevt_UUBAR.NevtwithoutCorrection / nevt_UUBAR.NevtwithCorrection;
    double NormFactor_DDBAR = nevt_DDBAR.NevtwithoutCorrection / nevt_DDBAR.NevtwithCorrection;
    double NormFactor_SSBAR = nevt_SSBAR.NevtwithoutCorrection / nevt_SSBAR.NevtwithCorrection;
    double NormFactor_CHARM = nevt_CHARM.NevtwithoutCorrection / nevt_CHARM.NevtwithCorrection;

    std::vector<std::string> variable_names;
    std::vector<std::string> branch_names;
    variable_names.push_back("MVA_BB"); branch_names.push_back("Upsilon");

    int Nvar = static_cast<int>(variable_names.size());
    if (Nvar != Nvar_num) exit(1);
    std::vector<double> Sideband_MC_values[Nvar_num];
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
    std::vector<int> Sideband_MC_numbering;

    std::vector<double> Sideband_data_values[Nvar_num];

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

    LetsFillSideBand_ri(Sideband_MC_CHG_train_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "CHG");
    LetsFillSideBand_ri(Sideband_MC_MIX_train_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "MIX");
    LetsFillSideBand_ri_correction(Sideband_MC_UUBAR_train_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "UUBAR", NormFactor_UUBAR);
    LetsFillSideBand_ri_correction(Sideband_MC_DDBAR_train_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "DDBAR", NormFactor_DDBAR);
    LetsFillSideBand_ri_correction(Sideband_MC_SSBAR_train_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "SSBAR", NormFactor_SSBAR);
    LetsFillSideBand_ri_correction(Sideband_MC_CHARM_train_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "CHARM", NormFactor_CHARM);
    LetsFillSideBand_ri(Sideband_MC_CHG_test_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "CHG");
    LetsFillSideBand_ri(Sideband_MC_MIX_test_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "MIX");
    LetsFillSideBand_ri_correction(Sideband_MC_UUBAR_test_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "UUBAR", NormFactor_UUBAR);
    LetsFillSideBand_ri_correction(Sideband_MC_DDBAR_test_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "DDBAR", NormFactor_DDBAR);
    LetsFillSideBand_ri_correction(Sideband_MC_SSBAR_test_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "SSBAR", NormFactor_SSBAR);
    LetsFillSideBand_ri_correction(Sideband_MC_CHARM_test_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "CHARM", NormFactor_CHARM);
    /*
    LetsFillSideBand_ri(Sideband_MC_CHG_train_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "CHG");
    LetsFillSideBand_ri(Sideband_MC_MIX_train_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "MIX");
    LetsFillSideBand_ri(Sideband_MC_UUBAR_train_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "UUBAR");
    LetsFillSideBand_ri(Sideband_MC_DDBAR_train_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "DDBAR");
    LetsFillSideBand_ri(Sideband_MC_SSBAR_train_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "SSBAR");
    LetsFillSideBand_ri(Sideband_MC_CHARM_train_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "CHARM");
    LetsFillSideBand_ri(Sideband_MC_CHG_test_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "CHG");
    LetsFillSideBand_ri(Sideband_MC_MIX_test_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "MIX");
    LetsFillSideBand_ri(Sideband_MC_UUBAR_test_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "UUBAR");
    LetsFillSideBand_ri(Sideband_MC_DDBAR_test_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "DDBAR");
    LetsFillSideBand_ri(Sideband_MC_SSBAR_test_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "SSBAR");
    LetsFillSideBand_ri(Sideband_MC_CHARM_test_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "CHARM");
    */
    LetsFill(Sideband_data_dirname, variable_names, branch_names, Sideband_data_values);

    // sort variables
    for (int k = 0; k < (int)Sideband_MC_numbering.size(); k++) {
        if (Sideband_MC_numbering.at(k) == 0) {
            for (int l = 0; l < (int)variable_names.size(); l++) charged_values[l].push_back(Sideband_MC_values[l].at(k));
            charged_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == 1) {
            for (int l = 0; l < (int)variable_names.size(); l++) mixed_values[l].push_back(Sideband_MC_values[l].at(k));
            mixed_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == 2) {
            for (int l = 0; l < (int)variable_names.size(); l++) uubar_values[l].push_back(Sideband_MC_values[l].at(k));
            uubar_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == 3) {
            for (int l = 0; l < (int)variable_names.size(); l++) ddbar_values[l].push_back(Sideband_MC_values[l].at(k));
            ddbar_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == 4) {
            for (int l = 0; l < (int)variable_names.size(); l++) ssbar_values[l].push_back(Sideband_MC_values[l].at(k));
            ssbar_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == 5) {
            for (int l = 0; l < (int)variable_names.size(); l++) ccbar_values[l].push_back(Sideband_MC_values[l].at(k));
            ccbar_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == 6) {
            for (int l = 0; l < (int)variable_names.size(); l++) taupair_values[l].push_back(Sideband_MC_values[l].at(k));
            taupair_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == 7) {
            for (int l = 0; l < (int)variable_names.size(); l++) mumu_values[l].push_back(Sideband_MC_values[l].at(k));
            mumu_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == 8) {
            for (int l = 0; l < (int)variable_names.size(); l++) gg_values[l].push_back(Sideband_MC_values[l].at(k));
            gg_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == 9) {
            for (int l = 0; l < (int)variable_names.size(); l++) ee_values[l].push_back(Sideband_MC_values[l].at(k));
            ee_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == 10) {
            for (int l = 0; l < (int)variable_names.size(); l++) eeee_values[l].push_back(Sideband_MC_values[l].at(k));
            eeee_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == 11) {
            for (int l = 0; l < (int)variable_names.size(); l++) eemumu_values[l].push_back(Sideband_MC_values[l].at(k));
            eemumu_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == 12) {
            for (int l = 0; l < (int)variable_names.size(); l++) llXX_values[l].push_back(Sideband_MC_values[l].at(k));
            llXX_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == 13) {
            for (int l = 0; l < (int)variable_names.size(); l++) hhISR_values[l].push_back(Sideband_MC_values[l].at(k));
            hhISR_weights.push_back(weights.at(k));
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
        temp_v.insert(temp_v.end(), Sideband_data_values[k].begin(), Sideband_data_values[k].end());


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
        for (int i = 0; i < (int)Sideband_MC_values[k].size(); i++) stat_error_hist[k]->Fill(Sideband_MC_values[k].at(i), weights.at(i));
        for (int i = 0; i < (int)Sideband_data_values[k].size(); i++) data_hist[k]->Fill(Sideband_data_values[k].at(i));

        for (int i = 0; i < (int)Sideband_MC_values[k].size(); i++) MC_one_bin[k]->Fill(Sideband_MC_values[k].at(i), weights.at(i));
        for (int i = 0; i < (int)Sideband_data_values[k].size(); i++) data_one_bin[k]->Fill(Sideband_data_values[k].at(i));
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
    printf("data: %d\n", (int)Sideband_data_values[0].size());

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
        c_temp->SaveAs((variable_names.at(k) + "_sideband.png").c_str());

        delete c_temp;
    }

    // Print data-MC discrepancy
    double MC_sum = 0;
    for (int i = 0; i < (int)Sideband_MC_values[0].size(); i++) MC_sum = MC_sum + weights.at(i);
    printf("data num: %ld\n", Sideband_data_values[0].size());
    printf("MC num with calibration: %lf\n", MC_sum);
    printf("MC with calibration: %lf +- %lf\n", MC_one_bin[0]->GetBinContent(1), MC_one_bin[0]->GetBinError(1));
    printf("data with calibration: %lf +- %lf\n", data_one_bin[0]->GetBinContent(1), data_one_bin[0]->GetBinError(1));
    printf("data/MC with calibration: %lf +- %lf\n", Ratio_one_bin[0]->GetBinContent(1), Ratio_one_bin[0]->GetBinError(1));
}
