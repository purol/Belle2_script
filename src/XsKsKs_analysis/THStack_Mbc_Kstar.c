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

#include "constants.h"

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

#include <algorithm>
#include <random>
#include <cmath>

std::vector<double> calculateMKSKS(double N_KSKS, double avg_E, double avg_px, double avg_py, double avg_pz, double first_E, double first_px, double first_py, double first_pz, double second_E, double second_px, double second_py, double second_pz, double Bsig_E, double Bsig_p, double Bsig_pt, double Bsig_phi){

if((N_KSKS > 1.5) && (N_KSKS < 2.5)){ // only one KSKS pair
    double sum_E = first_E + second_E;
    double sum_px = first_px + second_px;
    double sum_py = first_py + second_py;
    double sum_pz = first_pz + second_pz;
    std::vector<double> temp;
    temp.push_back(std::sqrt(sum_E*sum_E - sum_px*sum_px - sum_py*sum_py - sum_pz*sum_pz));
    return temp;
}
else if ((N_KSKS > 2.5) && (N_KSKS < 3.5)){

    double third_E = 3.0 * avg_E - (first_E + second_E);
    double third_px = 3.0 * avg_px - (first_px + second_px);
    double third_py = 3.0 * avg_py - (first_py + second_py);
    double third_pz = 3.0 * avg_pz - (first_pz + second_pz);

    double Bsig_pz = std::sqrt(Bsig_p * Bsig_p - Bsig_pt * Bsig_pt);
    double Bsig_py = Bsig_pt * std::sin(Bsig_phi);
    double Bsig_px = Bsig_pt * std::cos(Bsig_phi);

    double unknown_E = Bsig_E - 3.0 * avg_E;
    double unknown_px = Bsig_px - 3.0 * avg_px;
    double unknown_py = Bsig_py - 3.0 * avg_py;
    double unknown_pz = Bsig_pz - 3.0 * avg_pz;
    printf("Mass of unkown particle: %lf\n", std::sqrt(unknown_E*unknown_E-unknown_px*unknown_px-unknown_py*unknown_py-unknown_pz*unknown_pz));
    printf("Mass of B meson: %lf\n", std::sqrt(Bsig_E * Bsig_E - Bsig_px * Bsig_px - Bsig_py * Bsig_py - Bsig_pz * Bsig_pz));

    // variable preparation
    double spec_E = -1; // MXs_E
    double spec_px = -1;
    double spec_py = -1;
    double spec_pz = -1;
    double spec_M = -1;

    double KSKS_E = -1; // KSKS_E
    double KSKS_px = -1;
    double KSKS_py = -1;
    double KSKS_pz = -1;
    double KSKS_M = -1;

    std::vector<double> MKSKSs;

    // combine with first
    spec_E = Bsig_E - second_E - third_E;
    spec_px = Bsig_px - second_px - third_px;
    spec_py = Bsig_py - second_py - third_py;
    spec_pz = Bsig_pz - second_pz - third_pz;
    spec_M = std::sqrt(spec_E * spec_E - spec_px * spec_px - spec_py * spec_py - spec_pz * spec_pz);
    KSKS_E = second_E + third_E;
    KSKS_px = second_px + third_px;
    KSKS_py = second_py + third_py;
    KSKS_pz = second_pz + third_pz;
    KSKS_M = std::sqrt(KSKS_E * KSKS_E - KSKS_px * KSKS_px - KSKS_py * KSKS_py - KSKS_pz * KSKS_pz);
    if((spec_M > 0.6) && (spec_M < 1.0)) {
        MKSKSs.push_back(KSKS_M);
    }
printf("MXs: %lf\n", spec_M);
printf("MKSKS: %lf\n", KSKS_M);

    // combine with second
    spec_E = Bsig_E - first_E - third_E;
    spec_px = Bsig_px - first_px - third_px;
    spec_py = Bsig_py - first_py - third_py;
    spec_pz = Bsig_pz - first_pz - third_pz;
    spec_M = std::sqrt(spec_E * spec_E - spec_px * spec_px - spec_py * spec_py - spec_pz * spec_pz);
    KSKS_E = first_E + third_E;
    KSKS_px = first_px + third_px;
    KSKS_py = first_py + third_py;
    KSKS_pz = first_pz + third_pz;
    KSKS_M = std::sqrt(KSKS_E * KSKS_E - KSKS_px * KSKS_px - KSKS_py * KSKS_py - KSKS_pz * KSKS_pz);
    if((spec_M > 0.6) && (spec_M < 1.0)) {
        MKSKSs.push_back(KSKS_M);
    }
printf("MXs: %lf\n", spec_M);
printf("MKSKS: %lf\n", KSKS_M);

    // combine with third
    spec_E = Bsig_E - first_E - second_E;
    spec_px = Bsig_px - first_px - second_px;
    spec_py = Bsig_py - first_py - second_py;
    spec_pz = Bsig_pz - first_pz - second_pz;
    spec_M = std::sqrt(spec_E * spec_E - spec_px * spec_px - spec_py * spec_py - spec_pz * spec_pz);
    KSKS_E = first_E + second_E;
    KSKS_px = first_px + second_px;
    KSKS_py = first_py + second_py;
    KSKS_pz = first_pz + second_pz;
    KSKS_M = std::sqrt(KSKS_E * KSKS_E - KSKS_px * KSKS_px - KSKS_py * KSKS_py - KSKS_pz * KSKS_pz);
    if((spec_M > 0.6) && (spec_M < 1.0)) {
        MKSKSs.push_back(KSKS_M);
    }
printf("MXs: %lf\n", spec_M);
printf("MKSKS: %lf\n", KSKS_M);

    return MKSKSs;

}
else{
    printf("improper number of KS\n");
    exit(1);
}

}

THStack* Stack = new THStack("thstack", ";M_{bc} [GeV/c^{2}];number of candidates");
TH1D* SIGNAL_hist = new TH1D("signal", ";M_{bc} [GeV/c^{2}];number of candidates", 30, 5.24, 5.3);
TH1D* CHG_hist = new TH1D("charged", ";M_{bc} [GeV/c^{2}];number of candidates", 30, 5.24, 5.3);
TH1D* MIX_hist = new TH1D("mixed", ";M_{bc} [GeV/c^{2}];number of candidates", 30, 5.24, 5.3);
TH1D* UUBAR_hist = new TH1D("u#bar{u}", ";M_{bc} [GeV/c^{2}];number of candidates", 30, 5.24, 5.3);
TH1D* DDBAR_hist = new TH1D("d#bar{d}", ";M_{bc} [GeV/c^{2}];number of candidates", 30, 5.24, 5.3);
TH1D* SSBAR_hist = new TH1D("s#bar{s}", ";M_{bc} [GeV/c^{2}];number of candidates", 30, 5.24, 5.3);
TH1D* CHARM_hist = new TH1D("c#bar{c}", ";M_{bc} [GeV/c^{2}];number of candidates", 30, 5.24, 5.3);

void LetsFill(const char* dirname, TH1D* hist, double weight = 1, bool isSignal = false) {
    double var = 0;
    double isSignal_var = -1;
    double nBplusToXKSKS_var = -1;
    double nBzeroToXKSKS_var = -1;

    double NKS0 = -1;

    double KS0_avg_E = -1;
    double KS0_avg_px = -1;
    double KS0_avg_py = -1;
    double KS0_avg_pz = -1;

    double KS0_1st_E = -1;
    double KS0_1st_px = -1;
    double KS0_1st_py = -1;
    double KS0_1st_pz = -1;

    double KS0_2nd_E = -1;
    double KS0_2nd_px = -1;
    double KS0_2nd_py = -1;
    double KS0_2nd_pz = -1;

    double Bsig_E = -1;
    double Bsig_p = -1;
    double Bsig_pt = -1;
    double Bsig_phi = -1;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_Bsig->SetBranchAddress("Bsig_Mbc", &var);
        tree_Bsig->SetBranchAddress("Bsig_isSignal", &isSignal_var);

        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boN_KS_KSKS__bc__bc", &NKS0);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boavg_KS_KSKS_E__bc__bc", &KS0_avg_E);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boavg_KS_KSKS_px__bc__bc", &KS0_avg_px);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boavg_KS_KSKS_py__bc__bc", &KS0_avg_py);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boavg_KS_KSKS_pz__bc__bc", &KS0_avg_pz);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boKS_KSKS_1st_E__bc__bc", &KS0_1st_E);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boKS_KSKS_1st_px__bc__bc", &KS0_1st_px);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boKS_KSKS_1st_py__bc__bc", &KS0_1st_py);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boKS_KSKS_1st_pz__bc__bc", &KS0_1st_pz);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boKS_KSKS_2nd_E__bc__bc", &KS0_2nd_E);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boKS_KSKS_2nd_px__bc__bc", &KS0_2nd_px);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boKS_KSKS_2nd_py__bc__bc", &KS0_2nd_py);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boKS_KSKS_2nd_pz__bc__bc", &KS0_2nd_pz);

        tree_Bsig->SetBranchAddress("Bsig_E", &Bsig_E);
        tree_Bsig->SetBranchAddress("Bsig_p", &Bsig_p);
        tree_Bsig->SetBranchAddress("Bsig_pt", &Bsig_pt);
        tree_Bsig->SetBranchAddress("Bsig_phi", &Bsig_phi);

        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clXKSKS__bc", &nBplusToXKSKS_var);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clXKSKS__bc", &nBzeroToXKSKS_var);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            if(isSignal == true){
                if((names.at(i).find("CHG") == std::string::npos) && (names.at(i).find("MIX") == std::string::npos)) continue;
                else{
                    if((nBplusToXKSKS_var > 0.5) || (nBzeroToXKSKS_var > 0.5)) {}
                    else continue;
                }
            }
            else{
                if((names.at(i).find("CHG") == std::string::npos) && (names.at(i).find("MIX") == std::string::npos)) {}
                else{
                    if((nBplusToXKSKS_var > 0.5) || (nBzeroToXKSKS_var > 0.5)) continue;
                    else {}
                }
            }

            std::vector<double> MKSKSs = calculateMKSKS(NKS0, KS0_avg_E, KS0_avg_px, KS0_avg_py, KS0_avg_pz, KS0_1st_E, KS0_1st_px, KS0_1st_py, KS0_1st_pz, KS0_2nd_E, KS0_2nd_px, KS0_2nd_py, KS0_2nd_pz, Bsig_E, Bsig_p, Bsig_pt, Bsig_phi);

            //double equal_weight = 1.0 / MKSKSs.size();

            for(int i = 0; i < MKSKSs.size(); i++){
                double invM = MKSKSs.at(i);
                hist->Fill(var, weight);
            }
        }
        input_file->Close();

    }

}

void LetsFill(const char* dirname, TH1D* hist, const char* included_string, double weight = 1, bool isSignal = false) {
    double var = 0;
    double isSignal_var = -1;
    double nBplusToXKSKS_var = -1;
    double nBzeroToXKSKS_var = -1;

    double NKS0 = -1;

    double KS0_avg_E = -1;
    double KS0_avg_px = -1;
    double KS0_avg_py = -1;
    double KS0_avg_pz = -1;

    double KS0_1st_E = -1;
    double KS0_1st_px = -1;
    double KS0_1st_py = -1;
    double KS0_1st_pz = -1;

    double KS0_2nd_E = -1;
    double KS0_2nd_px = -1;
    double KS0_2nd_py = -1;
    double KS0_2nd_pz = -1;

    double Bsig_E = -1;
    double Bsig_p = -1;
    double Bsig_pt = -1;
    double Bsig_phi = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_Bsig->SetBranchAddress("Bsig_Mbc", &var);
        tree_Bsig->SetBranchAddress("Bsig_isSignal", &isSignal_var);

        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boN_KS_KSKS__bc__bc", &NKS0);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boavg_KS_KSKS_E__bc__bc", &KS0_avg_E);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boavg_KS_KSKS_px__bc__bc", &KS0_avg_px);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boavg_KS_KSKS_py__bc__bc", &KS0_avg_py);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boavg_KS_KSKS_pz__bc__bc", &KS0_avg_pz);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boKS_KSKS_1st_E__bc__bc", &KS0_1st_E);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boKS_KSKS_1st_px__bc__bc", &KS0_1st_px);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boKS_KSKS_1st_py__bc__bc", &KS0_1st_py);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boKS_KSKS_1st_pz__bc__bc", &KS0_1st_pz);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boKS_KSKS_2nd_E__bc__bc", &KS0_2nd_E);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boKS_KSKS_2nd_px__bc__bc", &KS0_2nd_px);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boKS_KSKS_2nd_py__bc__bc", &KS0_2nd_py);
        tree_upsilon->SetBranchAddress("daughter__bo1__cmextraInfo__boKS_KSKS_2nd_pz__bc__bc", &KS0_2nd_pz);

        tree_Bsig->SetBranchAddress("Bsig_E", &Bsig_E);
        tree_Bsig->SetBranchAddress("Bsig_p", &Bsig_p);
        tree_Bsig->SetBranchAddress("Bsig_pt", &Bsig_pt);
        tree_Bsig->SetBranchAddress("Bsig_phi", &Bsig_phi);

        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clXKSKS__bc", &nBplusToXKSKS_var);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clXKSKS__bc", &nBzeroToXKSKS_var);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            if(isSignal == true){
                if((names.at(i).find("CHG") == std::string::npos) && (names.at(i).find("MIX") == std::string::npos)) continue;
                else{
                    if((nBplusToXKSKS_var > 0.5) || (nBzeroToXKSKS_var > 0.5)) {}
                    else continue;
                }
            }
            else{
                if((names.at(i).find("CHG") == std::string::npos) && (names.at(i).find("MIX") == std::string::npos)) {}
                else{
                    if((nBplusToXKSKS_var > 0.5) || (nBzeroToXKSKS_var > 0.5)) continue;
                    else {}
                }
            }

            std::vector<double> MKSKSs = calculateMKSKS(NKS0, KS0_avg_E, KS0_avg_px, KS0_avg_py, KS0_avg_pz, KS0_1st_E, KS0_1st_px, KS0_1st_py, KS0_1st_pz, KS0_2nd_E, KS0_2nd_px, KS0_2nd_py, KS0_2nd_pz, Bsig_E, Bsig_p, Bsig_pt, Bsig_phi);

            //double equal_weight = 1.0 / MKSKSs.size();

            for(int i = 0; i < MKSKSs.size(); i++){
                double invM = MKSKSs.at(i);
                hist->Fill(var, weight);
            }

        }
        input_file->Close();

    }

}

void THStack_Mbc() {

    const char* CHG_dirname = "/home/belle2/junewoo/storage_ghi/pp_SS_systematic_BCS_bug_fix/SS/SS_sample_output_Kstar/CHG";
    const char* MIX_dirname = "/home/belle2/junewoo/storage_ghi/pp_SS_systematic_BCS_bug_fix/SS/SS_sample_output_Kstar/MIX";
    const char* UUBAR_dirname = "/home/belle2/junewoo/storage_ghi/pp_SS_systematic_BCS_bug_fix/SS/SS_sample_output_Kstar/UUBAR";
    const char* DDBAR_dirname = "/home/belle2/junewoo/storage_ghi/pp_SS_systematic_BCS_bug_fix/SS/SS_sample_output_Kstar/DDBAR";
    const char* SSBAR_dirname = "/home/belle2/junewoo/storage_ghi/pp_SS_systematic_BCS_bug_fix/SS/SS_sample_output_Kstar/SSBAR";
    const char* CHARM_dirname = "/home/belle2/junewoo/storage_ghi/pp_SS_systematic_BCS_bug_fix/SS/SS_sample_output_Kstar/CHARM";

    LetsFill(CHG_dirname, SIGNAL_hist, "final_final", Scale_CHG_validation_MC15rd, true);
    LetsFill(MIX_dirname, SIGNAL_hist, "final_final", Scale_MIX_validation_MC15rd, true);

    LetsFill(CHG_dirname, CHG_hist, "final_final", Scale_CHG_validation_MC15rd, false);
    LetsFill(MIX_dirname, MIX_hist, "final_final", Scale_MIX_validation_MC15rd, false);
    LetsFill(UUBAR_dirname, UUBAR_hist, "final_final", Scale_UUBAR_validation_MC15rd, false);
    LetsFill(DDBAR_dirname, DDBAR_hist, "final_final", Scale_DDBAR_validation_MC15rd, false);
    LetsFill(SSBAR_dirname, SSBAR_hist, "final_final", Scale_SSBAR_validation_MC15rd, false);
    LetsFill(CHARM_dirname, CHARM_hist, "final_final", Scale_CHARM_validation_MC15rd, false);


    double CHG_int = CHG_hist->Integral();
    double MIX_int = MIX_hist->Integral();
    double UUBAR_int = UUBAR_hist->Integral();
    double DDBAR_int = DDBAR_hist->Integral();
    double SSBAR_int = SSBAR_hist->Integral();
    double CHARM_int = CHARM_hist->Integral();

    double BKG_int = CHG_int + MIX_int + UUBAR_int + DDBAR_int + SSBAR_int + CHARM_int;

    /*
    CHG_hist->Scale(1.0 / BKG_int, "width");
    MIX_hist->Scale(1.0 / BKG_int, "width");
    UUBAR_hist->Scale(1.0 / BKG_int, "width");
    DDBAR_hist->Scale(1.0 / BKG_int, "width");
    SSBAR_hist->Scale(1.0 / BKG_int, "width");
    CHARM_hist->Scale(1.0 / BKG_int, "width");
    */

    TGaxis::SetMaxDigits(3);

    Stack->Add(CHG_hist);
    Stack->Add(MIX_hist);
    Stack->Add(UUBAR_hist);
    Stack->Add(DDBAR_hist);
    Stack->Add(SSBAR_hist);
    Stack->Add(CHARM_hist);
    Stack->SetMaximum(4.0);

    //SIGNAL_hist->Scale(10000.0/BKG_int, "width");
    //SIGNAL_hist->Scale(1.0 / SIGNAL_int, "width");
    //SIGNAL_hist->Scale(20.0);
    SIGNAL_hist->SetLineWidth(3);
    SIGNAL_hist->SetLineColor(2);
    SIGNAL_hist->SetFillStyle(0);

    TCanvas* c_temp = new TCanvas("c", "", 1500, 1200);
    c_temp->SetCanvasSize(1500, 1200);
    c_temp->SetLeftMargin(0.17);
    c_temp->SetRightMargin(0.04);
    c_temp->SetBottomMargin(0.16);
    c_temp->SetTopMargin(0.05);
    c_temp->cd();

    gStyle->SetPalette(kGistEarth);
    gStyle->SetOptStat(0);

    Stack->SetMinimum(0.0);
    Stack->Draw("pfc Hist");
    Stack->GetXaxis()->SetTitle("M_{bc} [GeV/c^{2}]");
    Stack->GetYaxis()->SetTitle("Number of candidates");
    Stack->GetXaxis()->SetLabelFont(42);
    Stack->GetYaxis()->SetLabelFont(42);
    Stack->GetXaxis()->SetTitleFont(42);
    Stack->GetYaxis()->SetTitleFont(42);
    Stack->GetXaxis()->SetLabelSize(0.040);
    Stack->GetYaxis()->SetLabelSize(0.040);
    Stack->GetXaxis()->SetTitleSize(0.045);
    Stack->GetYaxis()->SetTitleSize(0.045);
    Stack->GetXaxis()->SetTitleOffset(1.25);
    Stack->GetYaxis()->SetTitleOffset(1.40);
    Stack->GetXaxis()->CenterTitle();
    Stack->GetYaxis()->CenterTitle();
    Stack->GetXaxis()->SetNdivisions(505);
    Stack->GetYaxis()->SetNdivisions(505);
    SIGNAL_hist->Draw("HistSAME");
    c_temp->Modified();
    c_temp->Update();

    TLegend* legend = new TLegend(0.23, 0.76, 0.94, 0.93);
    legend->SetNColumns(3);
    legend->SetTextFont(42);
    legend->SetTextSize(0.040);
    legend->SetMargin(0.18);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->AddEntry(CHG_hist, "charged", "f");
    legend->AddEntry(MIX_hist, "mixed", "f");
    legend->AddEntry(UUBAR_hist, "u#bar{u}", "f");
    legend->AddEntry(DDBAR_hist, "d#bar{d}", "f");
    legend->AddEntry(SSBAR_hist, "s#bar{s}", "f");
    legend->AddEntry(CHARM_hist, "c#bar{c}", "f");
    legend->AddEntry(SIGNAL_hist, "B #rightarrow K#scale[0.93]{*}K_{S}^{0}K_{S}^{0}", "l");
    legend->Draw();

    c_temp->RedrawAxis();
    c_temp->Modified();
    c_temp->Update();
    c_temp->SaveAs("Plot_Mbc.png");
    c_temp->Print("Plot_Mbc.pdf", "pdf Portrait");

//    free(line);
//    gPad->BuildLegend();
//    c_temp->SaveAs("Plot_with_legend.png");
}
