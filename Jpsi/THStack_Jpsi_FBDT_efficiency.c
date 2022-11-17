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

// arXiv:1409.4557v2
# define TB0 1.5195 // (Table. 1)
# define TBp 1.6384 // (Table. 1)
# define BR_Kplus_nunubar 0.00000398 // (eq. 10)
# define BR_K0star_nunubar 0.00000919 // (eq. 11)
# define BR_K0_nunubar (BR_Kplus_nunubar*TB0/TBp) // under (eq. 15)
# define BR_Kplusstar_nunubar (BR_K0star_nunubar*TBp/TB0) // under (eq. 15)
# define BR_Xs_nunubar 0.000029 // (eq. 23)
# define BR_Xsu_nonresonant_nunubar (BR_Xs_nunubar - BR_Kplus_nunubar - BR_Kplusstar_nunubar)
# define BR_Xsd_nonresonant_nunubar (BR_Xs_nunubar - BR_K0_nunubar - BR_K0star_nunubar)

// https://confluence.desy.de/pages/viewpage.action?pageId=107054222
# define N_BpBp_1invab 565400000.0
# define N_B0B0_1invab 534600000.0

# define N_Kplus_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Kplus_nunubar)
# define N_Kplusstar_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Kplusstar_nunubar)
# define N_Xsu_nonresonant_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Xsu_nonresonant_nunubar)
# define N_K0_nunubar_1invab (2.0 * N_B0B0_1invab * BR_K0_nunubar)
# define N_K0star_nunubar_1invab (2.0 * N_B0B0_1invab * BR_K0star_nunubar)
# define N_Xsd_nunubar_1invab (2.0 * N_B0B0_1invab * BR_Xsd_nonresonant_nunubar)

// my MC sample number
# define N_Kplus_nunubar 10000000.0
# define N_K0_nunubar 10000000.0
# define N_Kplusstar_nunubar 10000000.0
# define N_K0star_nunubar 10000000.0
# define N_Xsu_nonresonant_nunubar 50000000.0
# define N_Xsd_nonresonant_nunubar 50000000.0

// scale factor for each MC sample
# define Scale_Kplus (N_Kplus_nunubar_1invab/N_Kplus_nunubar)
# define Scale_Kplusstar (N_Kplusstar_nunubar_1invab/N_Kplusstar_nunubar)
# define Scale_Xsu_nonresonant (N_Xsu_nonresonant_nunubar_1invab/N_Xsu_nonresonant_nunubar)
# define Scale_K0 (N_K0_nunubar_1invab/N_K0_nunubar)
# define Scale_K0star (N_K0star_nunubar_1invab/N_K0star_nunubar)
# define Scale_Xsd_nonresonant (N_Xsd_nunubar_1invab/N_Xsd_nonresonant_nunubar)

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

    double var[Nvar_num] = { 0.0 };
    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_KaonID_correction = -1;
    double temp_nKaon_excep = -1;

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
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_correction", &temp_KaonID_correction);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nKexcep", &temp_nKaon_excep);

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
            double Correction_KID = temp_KaonID_correction * std::pow(-1, temp_nKaon_excep);
            if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) { // B2Kc
                weights->push_back(FEI_calibration_factor *CAL * Stream * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) { // B2KcPi0
                weights->push_back(FEI_calibration_factor * pi0_correction*CAL * Stream * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) { // B2Ks0Pic
                weights->push_back(FEI_calibration_factor *CAL * Stream * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) { // B2KcPicPic
                weights->push_back(FEI_calibration_factor *CAL * Stream * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) { // B2Ks0PicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction*CAL * Stream * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) { // B2KcPicPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction*CAL * Stream * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) { // B2Ks0PicPicPic
                weights->push_back(FEI_calibration_factor *CAL * Stream * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) { // B2KcPicPicPicPic
                weights->push_back(FEI_calibration_factor *CAL * Stream * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) { // B2Ks0PicPicPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction*CAL * Stream * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) { // B2KcPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction*CAL * Stream * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) { // B2Ks0PicPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction*CAL * Stream * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) { // B2KcPicPicPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction*CAL * Stream * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) { // B2KcKcKc
                weights->push_back(FEI_calibration_factor *CAL * Stream * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) { // B2KcKcKs0Pic
                weights->push_back(FEI_calibration_factor *CAL * Stream * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) { // B2KcKcKcPi0
                weights->push_back(FEI_calibration_factor* pi0_correction*CAL* Stream * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) { // B02Ks0
                weights->push_back(FEI_calibration_factor*CAL* Stream * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) { // B02KcPic
                weights->push_back(FEI_calibration_factor*CAL* Stream * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) { // B02Ks0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction*CAL * Stream * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) { // B02KcPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction*CAL * Stream * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) { // B02Ks0PicPic
                weights->push_back(FEI_calibration_factor*CAL * Stream * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) { // B02KcPicPicPic
                weights->push_back(FEI_calibration_factor*CAL * Stream * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) { // B02Ks0PicPicPi0
                weights->push_back(FEI_calibration_factor* pi0_correction*CAL * Stream * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) { // B02KcPicPicPicPi0
                weights->push_back(FEI_calibration_factor* pi0_correction*CAL* Stream * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) { // B02Ks0PicPicPicPic
                weights->push_back(FEI_calibration_factor* CAL* Stream* Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) { // B02Ks0Pi0Pi0
                weights->push_back(FEI_calibration_factor* pi0_correction* pi0_correction* CAL* Stream* Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) { // B02KcPicPi0Pi0
                weights->push_back(FEI_calibration_factor* pi0_correction* pi0_correction* CAL* Stream* Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) { // B02Ks0PicPicPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * Stream * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) { // B02KcKcKs0
                weights->push_back(FEI_calibration_factor*CAL* Stream * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) { // B02KcKcKcPic
                weights->push_back(FEI_calibration_factor*CAL* Stream * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) { // B02KcKcKs0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction*CAL * Stream * Correction_KID);
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

    double var[Nvar_num] = { 0.0 };
    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_KaonID_correction = -1;
    double temp_nKaon_excep = -1;

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
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_correction", &temp_KaonID_correction);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nKexcep", &temp_nKaon_excep);
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
                weight_ri = (0.364436 / 0.8); // total 0.8/ab for BB
            }
            else if (SampleName == "MIX") {
                numberings->push_back(1);
                FEI_calibration_factor = FEI_cal_B0;
                weight_ri = (0.364436 / 0.8); // total 0.8/ab for BB
            }
            else if (SampleName == "UUBAR") {
                numberings->push_back(2);
                FEI_calibration_factor = CAL_qq;
                weight_ri = (0.364436 / 1.0); // total 1.0/ab for qq
            }
            else if (SampleName == "DDBAR") {
                numberings->push_back(3);
                FEI_calibration_factor = CAL_qq;
                weight_ri = (0.364436 / 1.0); // total 1.0/ab for qq
            }
            else if (SampleName == "SSBAR") {
                numberings->push_back(4);
                FEI_calibration_factor = CAL_qq;
                weight_ri = (0.364436 / 1.0); // total 1.0/ab for qq
            }
            else if (SampleName == "CHARM") {
                numberings->push_back(5);
                FEI_calibration_factor = CAL_qq;
                weight_ri = (0.364436 / 1.0); // total 1.0/ab for qq
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
                if(nXsu > 0) FEI_calibration_factor = FEI_cal_Bc;
                else if(nXsd > 0) FEI_calibration_factor = FEI_cal_B0;
                else if(Upsilon_ID > -0.5 && Upsilon_ID < 0.5) FEI_calibration_factor = FEI_cal_Bc;
                else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5) FEI_calibration_factor = FEI_cal_B0;
                else {
                    printf("ERROR 255");
                    exit(1);
                }
                weight_ri = (0.364436 / 25.0); // total 25.0/ab
            }
            else {
                printf("undefined job id!\n");
                exit(1);
            }

            // Fill calibration factors
            double Correction_KID = temp_KaonID_correction * std::pow(-1, temp_nKaon_excep);
            if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) { // B2Kc
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) { // B2KcPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) { // B2Ks0Pic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) { // B2KcPicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) { // B2Ks0PicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) { // B2KcPicPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) { // B2Ks0PicPicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) { // B2KcPicPicPicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) { // B2Ks0PicPicPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) { // B2KcPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) { // B2Ks0PicPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) { // B2KcPicPicPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) { // B2KcKcKc
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) { // B2KcKcKs0Pic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) { // B2KcKcKcPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) { // B02Ks0
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) { // B02KcPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) { // B02Ks0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) { // B02KcPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) { // B02Ks0PicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) { // B02KcPicPicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) { // B02Ks0PicPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) { // B02KcPicPicPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) { // B02Ks0PicPicPicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) { // B02Ks0Pi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) { // B02KcPicPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) { // B02Ks0PicPicPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) { // B02KcKcKs0
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) { // B02KcKcKcPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) { // B02KcKcKs0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID);
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

    double var[Nvar_num] = { 0.0 };
    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_KaonID_correction = -1;
    double temp_nKaon_excep = -1;

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
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_correction", &temp_KaonID_correction);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nKexcep", &temp_nKaon_excep);
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
                weight_ri = (0.364436 / 0.8); // total 0.8/ab for BB
            }
            else if (SampleName == "MIX") {
                numberings->push_back(1);
                FEI_calibration_factor = FEI_cal_B0;
                weight_ri = (0.364436 / 0.8); // total 0.8/ab for BB
            }
            else if (SampleName == "UUBAR") {
                numberings->push_back(2);
                FEI_calibration_factor = CAL_qq;
                weight_ri = (0.364436 / 1.0); // total 1.0/ab for qq
            }
            else if (SampleName == "DDBAR") {
                numberings->push_back(3);
                FEI_calibration_factor = CAL_qq;
                weight_ri = (0.364436 / 1.0); // total 1.0/ab for qq
            }
            else if (SampleName == "SSBAR") {
                numberings->push_back(4);
                FEI_calibration_factor = CAL_qq;
                weight_ri = (0.364436 / 1.0); // total 1.0/ab for qq
            }
            else if (SampleName == "CHARM") {
                numberings->push_back(5);
                FEI_calibration_factor = CAL_qq;
                weight_ri = (0.364436 / 1.0); // total 1.0/ab for qq
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
            double Correction_KID = temp_KaonID_correction * std::pow(-1, temp_nKaon_excep);
            if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) { // B2Kc
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) { // B2KcPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) { // B2Ks0Pic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) { // B2KcPicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) { // B2Ks0PicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) { // B2KcPicPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) { // B2Ks0PicPicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) { // B2KcPicPicPicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) { // B2Ks0PicPicPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) { // B2KcPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) { // B2Ks0PicPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) { // B2KcPicPicPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) { // B2KcKcKc
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) { // B2KcKcKs0Pic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) { // B2KcKcKcPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) { // B02Ks0
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) { // B02KcPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) { // B02Ks0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) { // B02KcPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) { // B02Ks0PicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) { // B02KcPicPicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) { // B02Ks0PicPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) { // B02KcPicPicPicPi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) { // B02Ks0PicPicPicPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) { // B02Ks0Pi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) { // B02KcPicPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) { // B02Ks0PicPicPi0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) { // B02KcKcKs0
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) { // B02KcKcKcPic
                weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) { // B02KcKcKs0Pi0
                weights->push_back(FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction);
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

    double var[Nvar_num] = { 0.0 };
    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_KaonID_correction = -1;
    double temp_nKaon_excep = -1;

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
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_correction", &temp_KaonID_correction);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nKexcep", &temp_nKaon_excep);
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
                weight_ri = (0.364436 / 0.8); // total 0.8/ab for BB
            }
            else if (SampleName == "MIX") {
                FEI_calibration_factor = FEI_cal_B0;
                weight_ri = (0.364436 / 0.8); // total 0.8/ab for BB
            }
            else if (SampleName == "UUBAR") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = (0.364436 / 1.0); // total 1.0/ab for qq
            }
            else if (SampleName == "DDBAR") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = (0.364436 / 1.0); // total 1.0/ab for qq
            }
            else if (SampleName == "SSBAR") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = (0.364436 / 1.0); // total 1.0/ab for qq
            }
            else if (SampleName == "CHARM") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = (0.364436 / 1.0); // total 1.0/ab for qq
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
            double Correction_KID = temp_KaonID_correction * std::pow(-1, temp_nKaon_excep);
            if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) { // B2Kc
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) { // B2KcPi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) { // B2Ks0Pic
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) { // B2KcPicPic
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) { // B2Ks0PicPi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) { // B2KcPicPicPi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) { // B2Ks0PicPicPic
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) { // B2KcPicPicPicPic
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) { // B2Ks0PicPicPicPi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) { // B2KcPi0Pi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) { // B2Ks0PicPi0Pi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) { // B2KcPicPicPi0Pi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) { // B2KcKcKc
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) { // B2KcKcKs0Pic
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) { // B2KcKcKcPi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) { // B02Ks0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) { // B02KcPic
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) { // B02Ks0Pi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) { // B02KcPicPi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) { // B02Ks0PicPic
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) { // B02KcPicPicPic
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) { // B02Ks0PicPicPi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) { // B02KcPicPicPicPi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) { // B02Ks0PicPicPicPic
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) { // B02Ks0Pi0Pi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) { // B02KcPicPi0Pi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) { // B02Ks0PicPicPi0Pi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) { // B02KcKcKs0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) { // B02KcKcKcPic
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_KID * BDTc_correction;
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) { // B02KcKcKs0Pi0
                nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID;
                nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * pi0_correction * CAL * weight_ri * Correction_KID * BDTc_correction;
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

    TH1D* MC_SIGNAL_before_one_bin = new TH1D("MC_SIGNAL_before_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* MC_CHG_before_one_bin = new TH1D("MC_CHG_before_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* MC_MIX_before_one_bin = new TH1D("MC_MIX_before_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* MC_UUBAR_before_one_bin = new TH1D("MC_UUBAR_before_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* MC_DDBAR_before_one_bin = new TH1D("MC_DDBAR_before_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* MC_SSBAR_before_one_bin = new TH1D("MC_SSBAR_before_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* MC_CHARM_before_one_bin = new TH1D("MC_CHARM_before_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);

    TH1D* MC_SIGNAL_after_one_bin = new TH1D("MC_SIGNAL_after_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* MC_CHG_after_one_bin = new TH1D("MC_CHG_after_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* MC_MIX_after_one_bin = new TH1D("MC_MIX_after_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* MC_UUBAR_after_one_bin = new TH1D("MC_UUBAR_after_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* MC_DDBAR_after_one_bin = new TH1D("MC_DDBAR_after_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* MC_SSBAR_after_one_bin = new TH1D("MC_SSBAR_after_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* MC_CHARM_after_one_bin = new TH1D("MC_CHARM_after_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);

    // Print data-MC discrepancy
    double MC_sum = 0;
    for (int i = 0; i < (int)Jpsi_MC_values[0].size(); i++) MC_sum = MC_sum + weights.at(i);
    printf("data num: %ld\n", Jpsi_data_values[0].size());
    printf("MC num with calibration: %lf\n", MC_sum);
    printf("MC with calibration: %lf +- %lf\n", MC_one_bin[0]->GetBinContent(1), MC_one_bin[0]->GetBinError(1));
    printf("data with calibration: %lf +- %lf\n", data_one_bin[0]->GetBinContent(1), data_one_bin[0]->GetBinError(1));
    printf("MC/data with calibration: %lf +- %lf\n", Ratio_one_bin[0]->GetBinContent(1), Ratio_one_bin[0]->GetBinError(1));
}
