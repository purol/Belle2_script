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

# define Nvar_num 116

# define CAL 1.1728
# define CAL_qq 1.0
# define Stream 0.25

# define Nvar_num 116

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

double LetsCount_ri(const char* dirname, std::string SampleName, double* Ncounts) {

    double Nevt = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };

    double temp_nKaonfrome = 0;
    double temp_nKaonfrommu = 0;
    double temp_nPionfrome = 0;
    double temp_nPionfrommu = 0;

    double FEI_calibration_factor = -1;

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
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nKfrome", &temp_nKaonfrome);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nKfrommu", &temp_nKaonfrommu);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_npifrome", &temp_nPionfrome);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_npifrommu", &temp_nPionfrommu);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);


            // get weight from luminosity and FEI
            double weight_ri = 0.0;
            if (SampleName == "CHG") {
                FEI_calibration_factor = FEI_cal_Bc;
                weight_ri = Scale_CHG_test;
            }
            else if (SampleName == "MIX") {
                FEI_calibration_factor = FEI_cal_B0;
                weight_ri = Scale_MIX_test;
            }
            else if (SampleName == "UUBAR") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = Scale_UUBAR_test; // total 1.0/ab for qq
            }
            else if (SampleName == "DDBAR") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = Scale_DDBAR_test; // total 1.0/ab for qq
            }
            else if (SampleName == "SSBAR") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = Scale_SSBAR_test; // total 1.0/ab for qq
            }
            else if (SampleName == "CHARM") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = Scale_CHARM_test; // total 1.0/ab for qq
            }
            else if (SampleName == "B2Knunu") {
                FEI_calibration_factor = FEI_cal_Bc;
                weight_ri = Scale_Kplus_test;
            }
            else if (SampleName == "B2Kstarnunu") {
                FEI_calibration_factor = FEI_cal_Bc;
                weight_ri = Scale_Kplusstar_test;
            }
            else if (SampleName == "B2Xsnunu") {
                FEI_calibration_factor = FEI_cal_Bc;
                weight_ri = Scale_Xsu_nonresonant_test;
            }
            else if (SampleName == "B02K0nunu") {
                FEI_calibration_factor = FEI_cal_B0;
                weight_ri = Scale_K0_test;
            }
            else if (SampleName == "B02Kstar0nunu") {
                FEI_calibration_factor = FEI_cal_B0;
                weight_ri = Scale_K0star_test;
            }
            else if (SampleName == "B02Xsnunu") {
                FEI_calibration_factor = FEI_cal_B0;
                weight_ri = Scale_Xsd_nonresonant_test;
            }
            else {
                printf("undefined job id!\n");
                exit(1);
            }

            // Fill calibration factors
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(PID_correction[0][i_PID], temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(PID_correction[1][i_PID], temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(PID_correction[2][i_PID], temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(PID_correction[3][i_PID], temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(pi0_correction[i_pi0], temp_N_bin_pi0[i_pi0]);

            int nKaon = 0;
            int nPion = 0;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                nKaon = nKaon + temp_N_bin_PID[0][i_PID];
                nKaon = nKaon + temp_N_bin_PID[1][i_PID];
                nPion = nPion + temp_N_bin_PID[2][i_PID];
                nPion = nPion + temp_N_bin_PID[3][i_PID];
            }

            Nevt = Nevt + FEI_calibration_factor * weight_ri * Correction_pi0 * Correction_KID * Correction_PID;
            Ncounts[0] = Ncounts[0] + nKaon * (FEI_calibration_factor * weight_ri * Correction_pi0 * Correction_KID * Correction_PID);
            Ncounts[1] = Ncounts[1] + nPion * (FEI_calibration_factor * weight_ri * Correction_pi0 * Correction_KID * Correction_PID);
            Ncounts[2] = Ncounts[2] + temp_nKaonfrome * (FEI_calibration_factor * weight_ri * Correction_pi0 * Correction_KID * Correction_PID);
            Ncounts[3] = Ncounts[3] + temp_nKaonfrommu * (FEI_calibration_factor * weight_ri * Correction_pi0 * Correction_KID * Correction_PID);
            Ncounts[4] = Ncounts[4] + temp_nPionfrome * (FEI_calibration_factor * weight_ri * Correction_pi0 * Correction_KID * Correction_PID);
            Ncounts[5] = Ncounts[5] + temp_nPionfrommu * (FEI_calibration_factor * weight_ri * Correction_pi0 * Correction_KID * Correction_PID);
        }
        input_file->Close();

    }

    return Nevt;

}

double LetsCount_ri(const char* dirname, std::string SampleName, const char* included_string, double* Ncounts) {

    double Nevt = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };

    double temp_nKaonfrome = 0;
    double temp_nKaonfrommu = 0;
    double temp_nPionfrome = 0;
    double temp_nPionfrommu = 0;

    double FEI_calibration_factor = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

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
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nKfrome", &temp_nKaonfrome);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nKfrommu", &temp_nKaonfrommu);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_npifrome", &temp_nPionfrome);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_npifrommu", &temp_nPionfrommu);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            // get weight from luminosity and FEI
            double weight_ri = 0.0;
            if (SampleName == "CHG") {
                FEI_calibration_factor = FEI_cal_Bc;
                weight_ri = Scale_CHG_test;
            }
            else if (SampleName == "MIX") {
                FEI_calibration_factor = FEI_cal_B0;
                weight_ri = Scale_MIX_test;
            }
            else if (SampleName == "UUBAR") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = Scale_UUBAR_test; // total 1.0/ab for qq
            }
            else if (SampleName == "DDBAR") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = Scale_DDBAR_test; // total 1.0/ab for qq
            }
            else if (SampleName == "SSBAR") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = Scale_SSBAR_test; // total 1.0/ab for qq
            }
            else if (SampleName == "CHARM") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = Scale_CHARM_test; // total 1.0/ab for qq
            }
            else if (SampleName == "B2Knunu") {
                FEI_calibration_factor = FEI_cal_Bc;
                weight_ri = Scale_Kplus_test;
            }
            else if (SampleName == "B2Kstarnunu") {
                FEI_calibration_factor = FEI_cal_Bc;
                weight_ri = Scale_Kplusstar_test;
            }
            else if (SampleName == "B2Xsnunu") {
                FEI_calibration_factor = FEI_cal_Bc;
                weight_ri = Scale_Xsu_nonresonant_test;
            }
            else if (SampleName == "B02K0nunu") {
                FEI_calibration_factor = FEI_cal_B0;
                weight_ri = Scale_K0_test;
            }
            else if (SampleName == "B02Kstar0nunu") {
                FEI_calibration_factor = FEI_cal_B0;
                weight_ri = Scale_K0star_test;
            }
            else if (SampleName == "B02Xsnunu") {
                FEI_calibration_factor = FEI_cal_B0;
                weight_ri = Scale_Xsd_nonresonant_test;
            }
            else {
                printf("undefined job id!\n");
                exit(1);
            }

            // Fill calibration factors
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(PID_correction[0][i_PID], temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(PID_correction[1][i_PID], temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(PID_correction[2][i_PID], temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(PID_correction[3][i_PID], temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(pi0_correction[i_pi0], temp_N_bin_pi0[i_pi0]);

            int nKaon = 0;
            int nPion = 0;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                nKaon = nKaon + temp_N_bin_PID[0][i_PID];
                nKaon = nKaon + temp_N_bin_PID[1][i_PID];
                nPion = nPion + temp_N_bin_PID[2][i_PID];
                nPion = nPion + temp_N_bin_PID[3][i_PID];
            }

            Nevt = Nevt + FEI_calibration_factor * weight_ri * Correction_pi0 * Correction_KID * Correction_PID;
            Ncounts[0] = Ncounts[0] + nKaon * (FEI_calibration_factor * weight_ri * Correction_pi0 * Correction_KID * Correction_PID);
            Ncounts[1] = Ncounts[1] + nPion * (FEI_calibration_factor * weight_ri * Correction_pi0 * Correction_KID * Correction_PID);
            Ncounts[2] = Ncounts[2] + temp_nKaonfrome * (FEI_calibration_factor * weight_ri * Correction_pi0 * Correction_KID * Correction_PID);
            Ncounts[3] = Ncounts[3] + temp_nKaonfrommu * (FEI_calibration_factor * weight_ri * Correction_pi0 * Correction_KID * Correction_PID);
            Ncounts[4] = Ncounts[4] + temp_nPionfrome * (FEI_calibration_factor * weight_ri * Correction_pi0 * Correction_KID * Correction_PID);
            Ncounts[5] = Ncounts[5] + temp_nPionfrommu * (FEI_calibration_factor * weight_ri * Correction_pi0 * Correction_KID * Correction_PID);
        }
        input_file->Close();

    }

    return Nevt;

}

void THStack_plot() {

    ReadPIDFile();

    const char* Knunu_dirname = "/home/jwpark/storage/BKG_gbasf2/Suwako/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut";
    const char* Kstarnunu_dirname = "/home/jwpark/storage/BKG_gbasf2/Suwako/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut";
    const char* Xsununu_dirname = "/home/jwpark/storage/BKG_gbasf2/Suwako/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut";
    const char* K0nunu_dirname = "/home/jwpark/storage/BKG_gbasf2/Suwako/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut";
    const char* K0starnunu_dirname = "/home/jwpark/storage/BKG_gbasf2/Suwako/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut";
    const char* Xsdnunu_dirname = "/home/jwpark/storage/BKG_gbasf2/Suwako/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut";
    const char* CHG_dirname = "/home/jwpark/storage/BKG_gbasf2/Suwako/CHG_analysis/test_v000/final_output_root_after_MVA_Application_after_cut";
    const char* MIX_dirname = "/home/jwpark/storage/BKG_gbasf2/Suwako/MIX_analysis/test_v000/final_output_root_after_MVA_Application_after_cut";
    const char* UUBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Suwako/UUBAR_analysis/test_v000/final_output_root_after_MVA_Application_after_cut";
    const char* DDBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Suwako/DDBAR_analysis/test_v000/final_output_root_after_MVA_Application_after_cut";
    const char* SSBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Suwako/SSBAR_analysis/test_v000/final_output_root_after_MVA_Application_after_cut";
    const char* CHARM_dirname = "/home/jwpark/storage/BKG_gbasf2/Suwako/CHARM_analysis/test_v000/final_output_root_after_MVA_Application_after_cut";
    
    /*
    * Ncounts:
    * [0]: the number of total kaon candidates
    * [1]: the number of total pion candidates
    * [2]: the number of kaon candidates from e
    * [3]: the number of kaon candidates from mu
    * [4]: the number of pion candidates from e
    * [5]: the number of pion candidates from mu
    */

    double* Ncounts_B2Knunu = (double*)malloc(sizeof(double) * 6);
    double* Ncounts_B2Kstarnunu = (double*)malloc(sizeof(double) * 6);
    double* Ncounts_B2Xsnunu = (double*)malloc(sizeof(double) * 6);
    double* Ncounts_B02K0nunu = (double*)malloc(sizeof(double) * 6);
    double* Ncounts_B02Kstar0nunu = (double*)malloc(sizeof(double) * 6);
    double* Ncounts_B02Xsnunu = (double*)malloc(sizeof(double) * 6);

    double* Ncounts_CHG = (double*)malloc(sizeof(double) * 6);
    double* Ncounts_MIX = (double*)malloc(sizeof(double) * 6);
    double* Ncounts_UUBAR = (double*)malloc(sizeof(double) * 6);
    double* Ncounts_DDBAR = (double*)malloc(sizeof(double) * 6);
    double* Ncounts_SSBAR = (double*)malloc(sizeof(double) * 6);
    double* Ncounts_CHARM = (double*)malloc(sizeof(double) * 6);

    for (int i = 0; i < 6; i++) {
        Ncounts_B2Knunu[i] = 0;
        Ncounts_B2Kstarnunu[i] = 0;
        Ncounts_B2Xsnunu[i] = 0;
        Ncounts_B02K0nunu[i] = 0;
        Ncounts_B02Kstar0nunu[i] = 0;
        Ncounts_B02Xsnunu[i] = 0;

        Ncounts_CHG[i] = 0;
        Ncounts_MIX[i] = 0;
        Ncounts_UUBAR[i] = 0;
        Ncounts_DDBAR[i] = 0;
        Ncounts_SSBAR[i] = 0;
        Ncounts_CHARM[i] = 0;
    }

    double Nevt_B2Knunu = LetsCount_ri(Knunu_dirname, "B2Knunu", "B2Knunu", Ncounts_B2Knunu);
    double Nevt_B2Kstarnunu = LetsCount_ri(Kstarnunu_dirname, "B2Kstarnunu", "B2Kstarnunu", Ncounts_B2Kstarnunu);
    double Nevt_B2Xsnunu = LetsCount_ri(Xsununu_dirname, "B2Xsnunu", "B2Xsnunu", Ncounts_B2Xsnunu);
    double Nevt_B02K0nunu = LetsCount_ri(K0nunu_dirname, "B02K0nunu", "B02K0nunu", Ncounts_B02K0nunu);
    double Nevt_B02Kstar0nunu = LetsCount_ri(K0starnunu_dirname, "B02Kstar0nunu", "B02Kstar0nunu", Ncounts_B02Kstar0nunu);
    double Nevt_B02Xsnunu = LetsCount_ri(Xsdnunu_dirname, "B02Xsnunu", "B02Xsnunu", Ncounts_B02Xsnunu);

    double Nevt_CHG = LetsCount_ri(CHG_dirname, "CHG", Ncounts_CHG);
    double Nevt_MIX = LetsCount_ri(MIX_dirname, "MIX", Ncounts_MIX);
    double Nevt_UUBAR = LetsCount_ri(UUBAR_dirname, "UUBAR", Ncounts_UUBAR);
    double Nevt_DDBAR = LetsCount_ri(DDBAR_dirname, "DDBAR", Ncounts_DDBAR);
    double Nevt_SSBAR = LetsCount_ri(SSBAR_dirname, "SSBAR", Ncounts_SSBAR);
    double Nevt_CHARM = LetsCount_ri(CHARM_dirname, "CHARM", Ncounts_CHARM);

    double nKaon_signal_all =
        Ncounts_B2Knunu[0] +
        Ncounts_B2Kstarnunu[0] +
        Ncounts_B2Xsnunu[0] +
        Ncounts_B02K0nunu[0] +
        Ncounts_B02Kstar0nunu[0] +
        Ncounts_B02Xsnunu[0];

    double nKaon_bkg_all =
        Ncounts_CHG[0] +
        Ncounts_MIX[0] +
        Ncounts_UUBAR[0] +
        Ncounts_DDBAR[0] +
        Ncounts_SSBAR[0] +
        Ncounts_CHARM[0];

    double nPion_signal_all =
        Ncounts_B2Knunu[1] +
        Ncounts_B2Kstarnunu[1] +
        Ncounts_B2Xsnunu[1] +
        Ncounts_B02K0nunu[1] +
        Ncounts_B02Kstar0nunu[1] +
        Ncounts_B02Xsnunu[1];

    double nPion_bkg_all =
        Ncounts_CHG[1] +
        Ncounts_MIX[1] +
        Ncounts_UUBAR[1] +
        Ncounts_DDBAR[1] +
        Ncounts_SSBAR[1] +
        Ncounts_CHARM[1];

    double nKaon_signal_frome =
        Ncounts_B2Knunu[2] +
        Ncounts_B2Kstarnunu[2] +
        Ncounts_B2Xsnunu[2] +
        Ncounts_B02K0nunu[2] +
        Ncounts_B02Kstar0nunu[2] +
        Ncounts_B02Xsnunu[2];

    double nKaon_bkg_frome =
        Ncounts_CHG[2] +
        Ncounts_MIX[2] +
        Ncounts_UUBAR[2] +
        Ncounts_DDBAR[2] +
        Ncounts_SSBAR[2] +
        Ncounts_CHARM[2];

    double nKaon_signal_frommu =
        Ncounts_B2Knunu[3] +
        Ncounts_B2Kstarnunu[3] +
        Ncounts_B2Xsnunu[3] +
        Ncounts_B02K0nunu[3] +
        Ncounts_B02Kstar0nunu[3] +
        Ncounts_B02Xsnunu[3];

    double nKaon_bkg_frommu =
        Ncounts_CHG[3] +
        Ncounts_MIX[3] +
        Ncounts_UUBAR[3] +
        Ncounts_DDBAR[3] +
        Ncounts_SSBAR[3] +
        Ncounts_CHARM[3];

    double nPion_signal_frome =
        Ncounts_B2Knunu[4] +
        Ncounts_B2Kstarnunu[4] +
        Ncounts_B2Xsnunu[4] +
        Ncounts_B02K0nunu[4] +
        Ncounts_B02Kstar0nunu[4] +
        Ncounts_B02Xsnunu[4];

    double nPion_bkg_frome =
        Ncounts_CHG[4] +
        Ncounts_MIX[4] +
        Ncounts_UUBAR[4] +
        Ncounts_DDBAR[4] +
        Ncounts_SSBAR[4] +
        Ncounts_CHARM[4];

    double nPion_signal_frommu =
        Ncounts_B2Knunu[5] +
        Ncounts_B2Kstarnunu[5] +
        Ncounts_B2Xsnunu[5] +
        Ncounts_B02K0nunu[5] +
        Ncounts_B02Kstar0nunu[5] +
        Ncounts_B02Xsnunu[5];

    double nPion_bkg_frommu =
        Ncounts_CHG[5] +
        Ncounts_MIX[5] +
        Ncounts_UUBAR[5] +
        Ncounts_DDBAR[5] +
        Ncounts_SSBAR[5] +
        Ncounts_CHARM[5];

    printf("all Kaon at signal: %lf\n", nKaon_signal_all);
    printf("all Kaon at background: %lf\n", nKaon_bkg_all);
    printf("all Pion at signal: %lf\n", nPion_signal_all);
    printf("all Pion at background: %lf\n", nPion_bkg_all);
    printf("Kaon from e at signal: %lf\n", nKaon_signal_frome);
    printf("Kaon from e at background: %lf\n", nKaon_bkg_frome);
    printf("Kaon from mu at signal: %lf\n", nKaon_signal_frommu);
    printf("Kaon from mu at background: %lf\n", nKaon_bkg_frommu);
    printf("Pion from e at signal: %lf\n", nPion_signal_frome);
    printf("Pion from e at background: %lf\n", nPion_bkg_frome);
    printf("Pion from mu at signal: %lf\n", nPion_signal_frommu);
    printf("Pion from mu at background: %lf\n", nPion_bkg_frommu);

    printf("Kaon from lepton/all Kaon at signal: %lf\n", (nKaon_signal_frome+ nKaon_signal_frommu)/ nKaon_signal_all);
    printf("Pion from lepton/all Pion at signal: %lf\n", (nPion_signal_frome + nPion_signal_frommu) / nPion_signal_all);
    printf("Kaon from lepton/all Kaon at bkg: %lf\n", (nKaon_bkg_frome + nKaon_bkg_frommu) / nKaon_bkg_all);
    printf("Pion from lepton/all Pion at bkg: %lf\n", (nPion_bkg_frome + nPion_bkg_frommu) / nPion_bkg_all);
}
