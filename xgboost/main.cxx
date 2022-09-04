#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
# include <vector>
#include <fstream>

#include <TMath.h>
#include <TColor.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TString.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THStack.h>
#include <TPaveText.h>
#include <TKey.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <numeric>
#include <cfloat>

#include <xgboost/c_api.h>

# define N_Needed_info 37
# define N_event_info 15
# define N_Upsilon_info 47
# define N_Bsig_info 78
# define N_Btag_info 7
# define N_decay 38 // five decay mode + others

# define Nvar 45
# define DvetoNvar 7

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

using std::string;

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

int GetEntryNum(const char* dirname, const char* included_string) {
    unsigned int N_Entry = 0;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    for (unsigned int i = 0; i < names.size(); ++i) {
        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        TTree* tree_data = (TTree*)input_file->Get("data");
        N_Entry = N_Entry + tree_data->GetEntries();
        input_file->Close();
    }

    return (int)N_Entry;

}

int GetEntryNum(const char* dirname, bool tempissignal) {
    unsigned int N_Entry = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    double Decay_Kplus = -1;
    double Decay_Kplusstar_ch1 = -1;
    double Decay_Kplusstar_ch2 = -1;
    double Decay_K0 = -1;
    double Decay_K0star_ch1 = -1;
    double Decay_K0star_ch2 = -1;

    double Mxs = -1;

    for (unsigned int i = 0; i < names.size(); ++i) {
        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        TTree* tree_data = (TTree*)input_file->Get("data");

        if (tempissignal) {
            tree_data->SetBranchAddress("nParticlesInList__boB__pl__clKcharge_total__bc", &Decay_Kplus);
            tree_data->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &Decay_Kplusstar_ch1);
            tree_data->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &Decay_Kplusstar_ch2);
            tree_data->SetBranchAddress("nParticlesInList__boB0__clKneutral_total__bc", &Decay_K0);
            tree_data->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &Decay_K0star_ch1);
            tree_data->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &Decay_K0star_ch2);
            tree_data->SetBranchAddress("Bsig_M", &Mxs);

            for (unsigned int j = 0; j < tree_data->GetEntries(); j++) {
                tree_data->GetEntry(j);

                if ((Decay_Kplus > 0.5 || Decay_Kplusstar_ch1 > 0.5 || Decay_Kplusstar_ch2 > 0.5 || Decay_K0 > 0.5 || Decay_K0star_ch1 > 0.5 || Decay_K0star_ch2 > 0.5) && Mxs > 1.1) continue;
                else if ((Decay_Kplus < 0.5 && Decay_Kplusstar_ch1 < 0.5 && Decay_Kplusstar_ch2 < 0.5 && Decay_K0 < 0.5 && Decay_K0star_ch1 < 0.5 && Decay_K0star_ch2 < 0.5) && Mxs < 1.1) continue;

                N_Entry++;
            }
        }
        else {
            N_Entry = N_Entry + tree_data->GetEntries();
        }

        input_file->Close();
    }

    return (int)N_Entry;

}

void FillVariables(const char* filename, float * input_data, float* IsSignal, float* weight, bool tempissignal, int* indicator, float weight_N = 1.0) {
    TFile* input_file = new TFile(filename, "read");

    TTree* tree_data = (TTree*)input_file->Get("data");

    double Vars[Nvar];
    int flag;

    double Dc_chiProb;
    double Dc_pvalue_med;
    double Dc_pvalue_std;
    double Dc_dr;
    double Dc_dz;
    double Dc_M;
    double D0_chiProb;
    double D0_pvalue_med;
    double D0_pvalue_std;
    double D0_dr;
    double D0_dz;
    double D0_M;

    double Decay_Kplus = -1;
    double Decay_Kplusstar_ch1 = -1;
    double Decay_Kplusstar_ch2 = -1;
    double Decay_K0 = -1;
    double Decay_K0star_ch1 = -1;
    double Decay_K0star_ch2 = -1;

    double Mxs = -1;
    double Pcms = -1;

    tree_data->SetBranchAddress("aplanarity", &Vars[0]);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_mydz", &Vars[1]);
    tree_data->SetBranchAddress("Btag_CleoConeCS_1", &Vars[2]);
    tree_data->SetBranchAddress("Btag_CleoConeCS_2", &Vars[3]);
    tree_data->SetBranchAddress("Btag_CleoConeCS_3", &Vars[4]);
    tree_data->SetBranchAddress("Btag_CleoConeCS_4", &Vars[5]);
    tree_data->SetBranchAddress("Btag_CleoConeCS_5", &Vars[6]);
    tree_data->SetBranchAddress("Btag_CleoConeCS_6", &Vars[7]);
    tree_data->SetBranchAddress("Btag_CleoConeCS_7", &Vars[8]);
    tree_data->SetBranchAddress("Btag_CleoConeCS_8", &Vars[9]);
    tree_data->SetBranchAddress("Btag_cosTBTO", &Vars[10]);
    tree_data->SetBranchAddress("Btag_cosTBz", &Vars[11]);
    tree_data->SetBranchAddress("Btag_deltaE", &Vars[12]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hso01", &Vars[13]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hso03", &Vars[14]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hso04", &Vars[15]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hso14", &Vars[16]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hso24", &Vars[17]);
    tree_data->SetBranchAddress("Btag_thrustBm", &Vars[18]);
    tree_data->SetBranchAddress("Btag_thrustOm", &Vars[19]);
    tree_data->SetBranchAddress("Btag_useCMSFrame_theta", &Vars[20]);
    tree_data->SetBranchAddress("cleoConeThrust0", &Vars[21]);
    tree_data->SetBranchAddress("cleoConeThrust1", &Vars[22]);
    tree_data->SetBranchAddress("cleoConeThrust2", &Vars[23]);
    tree_data->SetBranchAddress("cleoConeThrust3", &Vars[24]);
    tree_data->SetBranchAddress("cleoConeThrust4", &Vars[25]);
    tree_data->SetBranchAddress("cleoConeThrust5", &Vars[26]);
    tree_data->SetBranchAddress("cleoConeThrust6", &Vars[27]);
    tree_data->SetBranchAddress("cleoConeThrust7", &Vars[28]);
    tree_data->SetBranchAddress("cleoConeThrust8", &Vars[29]);
    tree_data->SetBranchAddress("harmonicMomentThrust2", &Vars[30]);
    tree_data->SetBranchAddress("harmonicMomentThrust3", &Vars[31]);
    tree_data->SetBranchAddress("harmonicMomentThrust4", &Vars[32]);
    tree_data->SetBranchAddress("Btag_extraInfo_SignalProbability", &Vars[33]);
    tree_data->SetBranchAddress("missingMomentumOfEvent_theta", &Vars[34]);
    tree_data->SetBranchAddress("nParticlesInList__bomu__pl__clMuonFBDT_tight__bc", &Vars[35]);
    tree_data->SetBranchAddress("roeEextra__bocleanMask__bc", &Vars[36]);
    tree_data->SetBranchAddress("roePTheta__bocleanMask__bc", &Vars[37]);

    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dc_pValue_std", &Dc_pvalue_std);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_dr", &Dc_dr);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_dz", &Dc_dz);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_M", &Dc_M);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_D0simpleveto_chiProb", &D0_chiProb);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_D0simpleveto_dz", &D0_dz);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_D0simpleveto_M", &D0_M);

    tree_data->SetBranchAddress("flag", &flag);

    if (tempissignal) {
        tree_data->SetBranchAddress("nParticlesInList__boB__pl__clKcharge_total__bc", &Decay_Kplus);
        tree_data->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &Decay_Kplusstar_ch1);
        tree_data->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &Decay_Kplusstar_ch2);
        tree_data->SetBranchAddress("nParticlesInList__boB0__clKneutral_total__bc", &Decay_K0);
        tree_data->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &Decay_K0star_ch1);
        tree_data->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &Decay_K0star_ch2);
    }
    tree_data->SetBranchAddress("Bsig_M", &Mxs);
    tree_data->SetBranchAddress("Bsig_useCMSFrame_p", &Pcms);

    //printf("%lld entries...\n", tree_data->GetEntries());
    for (unsigned int j = 0; j < tree_data->GetEntries(); j++) { // Fill
        tree_data->GetEntry(j);
        if (tempissignal == true && (Decay_Kplus > 0.5 || Decay_Kplusstar_ch1 > 0.5 || Decay_Kplusstar_ch2 > 0.5 || Decay_K0 > 0.5 || Decay_K0star_ch1 > 0.5 || Decay_K0star_ch2 > 0.5) && Mxs > 1.1) continue;
        else if (tempissignal == true && (Decay_Kplus < 0.5 && Decay_Kplusstar_ch1 < 0.5 && Decay_Kplusstar_ch2 < 0.5 && Decay_K0 < 0.5 && Decay_K0star_ch1 < 0.5 && Decay_K0star_ch2 < 0.5) && Mxs < 1.1) continue;

        for (unsigned int k = 0; k < Nvar - DvetoNvar; k++) {
            input_data[*indicator] = (float)Vars[k];
            *indicator = *indicator + 1;
        }

        if (Dc_chiProb > -0.5) {
            input_data[*indicator] = (float)Dc_pvalue_std;
            *indicator = *indicator + 1;
            input_data[*indicator] = (float)Dc_dr;
            *indicator = *indicator + 1;
            input_data[*indicator] = (float)Dc_dz;
            *indicator = *indicator + 1;
            input_data[*indicator] = (float)Dc_M;
            *indicator = *indicator + 1;
        }
        else {
            input_data[*indicator] = 0.0f;
            *indicator = *indicator + 1;
            input_data[*indicator] = -1.0f;
            *indicator = *indicator + 1;
            input_data[*indicator] = -100.0f;
            *indicator = *indicator + 1;
            input_data[*indicator] = 0.0f;
            *indicator = *indicator + 1;
        }
        if (D0_chiProb > -0.5) {
            input_data[*indicator] = (float)D0_chiProb;
            *indicator = *indicator + 1;
            input_data[*indicator] = (float)D0_dz;
            *indicator = *indicator + 1;
            input_data[*indicator] = (float)D0_M;
            *indicator = *indicator + 1;
        }
        else {
            input_data[*indicator] = 0.0f;
            *indicator = *indicator + 1;
            input_data[*indicator] = -100.0f;
            *indicator = *indicator + 1;
            input_data[*indicator] = 0.0f;
            *indicator = *indicator + 1;
        }

        if (tempissignal) IsSignal[(*indicator - 1) / Nvar] = 1.0f;
        else IsSignal[(*indicator - 1) / (int) Nvar] = 0.0f;
        weight[(*indicator - 1) / (int) Nvar] = weight_N;

    }

    input_file->Close();
}

double CalculateAUC(bst_ulong out_len, const float* f, float* IsSignal, float* weight) {
    const int step = 100;
    double AUC = 0;
    double NBKG_total = 0;
    double NSIG_total = 0;
    std::vector<double> TPRs;
    std::vector<double> FPRs;

    for (unsigned int i = 0; i < out_len; i++) {
        if (IsSignal[i] > 0.5f) NSIG_total = NSIG_total + weight[i];
        else NBKG_total = NBKG_total + weight[i];
    }

    for (int i = 0; i < step; i++) {
        float value = ((float)i) / ((float)step);
        double NBKG = 0;
        double NSIG = 0;

        for (unsigned int k = 0; k < out_len; k++) {
            float p = f[k];
            if (p >= value) {
                if (IsSignal[k] > 0.5f) NSIG = NSIG + weight[k];
                else NBKG = NBKG + weight[k];
            }
        }

        double TPR = NSIG / NSIG_total;
        double FPR = NBKG / NBKG_total;

        TPRs.push_back(TPR);
        FPRs.push_back(FPR);
    }

    for (unsigned int i = 0; i < TPRs.size(); ++i) {
        if (i != TPRs.size() - 1) {
            double del_FPR = FPRs.at(i) - FPRs.at(i + 1);
            double avg_TPR = (TPRs.at(i) + TPRs.at(i + 1)) / 2.0;
            AUC = AUC + del_FPR * avg_TPR;
        }
        else {
            double del_FPR = FPRs.at(i) - 0.0;
            double avg_TPR = (TPRs.at(i) + 0.0) / 2.0;
            AUC = AUC + del_FPR * avg_TPR;
        }
    }

    return AUC;

}

#define safe_xgboost(call) {  \
  int err = (call); \
  if (err != 0) { \
    fprintf(stderr, "%s:%d: error in %s: %s\n", __FILE__, __LINE__, #call, XGBGetLastError());  \
    exit(1); \
  } \
}

int main(int argc, char** argv) {
    // set name
    std::string model_name = std::string(argv[1]) + "_" + std::string(argv[2]) + "_" + std::string(argv[3]) + "_" + std::string(argv[4]) + "_" + std::string(argv[5]) + "_" + std::string(argv[6]) + "_" + std::string(argv[7]);



    // AUCs
    double AUC_train = 0.0;
    double AUC_test = 0.0;



    // input file
    const char* SIGNAL_input_train = "/home/belle2/junewoo/storage_b1/GridSearch/SIGNAL_analysis/train_v001/final_output/DataFile";
    const char* CHG_input_train = "/home/belle2/junewoo/storage_b1/GridSearch/CHG_analysis/train_v001/final_output/DataFile";
    const char* MIX_input_train = "/home/belle2/junewoo/storage_b1/GridSearch/MIX_analysis/train_v001/final_output/DataFile";
    const char* UUBAR_input_train = "/home/belle2/junewoo/storage_b1/GridSearch/UUBAR_analysis/train_v001/final_output/DataFile";
    const char* DDBAR_input_train = "/home/belle2/junewoo/storage_b1/GridSearch/DDBAR_analysis/train_v001/final_output/DataFile";
    const char* SSBAR_input_train = "/home/belle2/junewoo/storage_b1/GridSearch/SSBAR_analysis/train_v001/final_output/DataFile";
    const char* CHARM_input_train = "/home/belle2/junewoo/storage_b1/GridSearch/CHARM_analysis/train_v001/final_output/DataFile";

    const char* SIGNAL_input_test = "/home/belle2/junewoo/storage_b1/GridSearch/SIGNAL_analysis/test_v001/final_output/DataFile";
    const char* CHG_input_test = "/home/belle2/junewoo/storage_b1/GridSearch/CHG_analysis/test_v001/final_output/DataFile";
    const char* MIX_input_test = "/home/belle2/junewoo/storage_b1/GridSearch/MIX_analysis/test_v001/final_output/DataFile";
    const char* UUBAR_input_test = "/home/belle2/junewoo/storage_b1/GridSearch/UUBAR_analysis/test_v001/final_output/DataFile";
    const char* DDBAR_input_test = "/home/belle2/junewoo/storage_b1/GridSearch/DDBAR_analysis/test_v001/final_output/DataFile";
    const char* SSBAR_input_test = "/home/belle2/junewoo/storage_b1/GridSearch/SSBAR_analysis/test_v001/final_output/DataFile";
    const char* CHARM_input_test = "/home/belle2/junewoo/storage_b1/GridSearch/CHARM_analysis/test_v001/final_output/DataFile";

    int N_entry_train = GetEntryNum(SIGNAL_input_train, true) +
        GetEntryNum(CHG_input_train, false) +
        GetEntryNum(MIX_input_train, false) +
        GetEntryNum(UUBAR_input_train, false) +
        GetEntryNum(DDBAR_input_train, false) +
        GetEntryNum(SSBAR_input_train, false) +
        GetEntryNum(CHARM_input_train, false);

    int N_entry_test = GetEntryNum(SIGNAL_input_test, true) +
        GetEntryNum(CHG_input_test, false) +
        GetEntryNum(MIX_input_test, false) +
        GetEntryNum(UUBAR_input_test, false) +
        GetEntryNum(DDBAR_input_test, false) +
        GetEntryNum(SSBAR_input_test, false) +
        GetEntryNum(CHARM_input_test, false);




    // set training sample
    float* train_set; int train_indicator = 0; float* train_weight; float* train_IsSignal;
    train_set = (float*)malloc(sizeof(float) * N_entry_train * Nvar);
    train_weight = (float*)malloc(sizeof(float) * N_entry_train);
    train_IsSignal = (float*)malloc(sizeof(float) * N_entry_train);

    // fill training sample
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B2Knunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), train_set, train_IsSignal, train_weight, true, &train_indicator, 3.0 / 7.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B2Kstarnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), train_set, train_IsSignal, train_weight, true, &train_indicator, 3.0 / 7.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B2Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), train_set, train_IsSignal, train_weight, true, &train_indicator, 3.0 / 7.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B02K0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), train_set, train_IsSignal, train_weight, true, &train_indicator, 3.0 / 7.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B02Kstar0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), train_set, train_IsSignal, train_weight, true, &train_indicator, 3.0 / 7.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B02Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), train_set, train_IsSignal, train_weight, true, &train_indicator, 3.0 / 7.0);
        }
    }
    {
        std::vector<string> names;
        load_files(CHG_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHG_input_train + std::string("/") + names.at(i)).c_str(), train_set, train_IsSignal, train_weight, false, &train_indicator, (0.3 / 1.0));
        }
    }
    {
        std::vector<string> names;
        load_files(MIX_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((MIX_input_train + std::string("/") + names.at(i)).c_str(), train_set, train_IsSignal, train_weight, false, &train_indicator, (0.3 / 1.0));
        }
    }
    {
        std::vector<string> names;
        load_files(UUBAR_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((UUBAR_input_train + std::string("/") + names.at(i)).c_str(), train_set, train_IsSignal, train_weight, false, &train_indicator, (0.3 / 1.1));
        }
    }
    {
        std::vector<string> names;
        load_files(DDBAR_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((DDBAR_input_train + std::string("/") + names.at(i)).c_str(), train_set, train_IsSignal, train_weight, false, &train_indicator, (0.3 / 1.1));
        }
    }
    {
        std::vector<string> names;
        load_files(SSBAR_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SSBAR_input_train + std::string("/") + names.at(i)).c_str(), train_set, train_IsSignal, train_weight, false, &train_indicator, (0.3 / 1.1));
        }
    }
    {
        std::vector<string> names;
        load_files(CHARM_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHARM_input_train + std::string("/") + names.at(i)).c_str(), train_set, train_IsSignal, train_weight, false, &train_indicator, (0.3 / 1.1));
        }
    }

    // convert to DMatrix
    DMatrixHandle h_train;
    safe_xgboost(XGDMatrixCreateFromMat(train_set, N_entry_train, Nvar, FLT_MAX, &h_train));

    // load the labels
    safe_xgboost(XGDMatrixSetFloatInfo(h_train, "label", train_IsSignal, N_entry_train));


    // create the booster and load some parameters
    BoosterHandle h_booster;
    safe_xgboost(XGBoosterCreate(&h_train, 1, &h_booster));
    safe_xgboost(XGBoosterSetParam(h_booster, "booster", "gbtree"));
    safe_xgboost(XGBoosterSetParam(h_booster, "objective", "binary:logistic"));
    safe_xgboost(XGBoosterSetParam(h_booster, "eta", argv[1])); // 0.1 0.2 0.3 0.4 0.5 | 0.3
    safe_xgboost(XGBoosterSetParam(h_booster, "gamma", argv[2])); // 0 1 10 | 0
    safe_xgboost(XGBoosterSetParam(h_booster, "max_depth", argv[3])); // 3 4 5 6 7 | 6
    safe_xgboost(XGBoosterSetParam(h_booster, "min_child_weight", argv[4])); // 0.5 1.0 1.5 | 1
    safe_xgboost(XGBoosterSetParam(h_booster, "max_delta_step", "0"));
    safe_xgboost(XGBoosterSetParam(h_booster, "subsample", argv[5])); // 0.3 0.6 1.0 | 1
    safe_xgboost(XGBoosterSetParam(h_booster, "sampling_method", "uniform"));
    safe_xgboost(XGBoosterSetParam(h_booster, "colsample_bytree", argv[6])); // 0.3 0.5 1.0 | 1
    safe_xgboost(XGBoosterSetParam(h_booster, "colsample_bylevel", "1"));
    safe_xgboost(XGBoosterSetParam(h_booster, "colsample_bynode", "1"));
    safe_xgboost(XGBoosterSetParam(h_booster, "lambda", "1"));
    safe_xgboost(XGBoosterSetParam(h_booster, "alpha", "0"));
    safe_xgboost(XGBoosterSetParam(h_booster, "num_parallel_tree", "1"));
    const int Niter = atoi(argv[7]); // 100 500 1000 1500 2000



    // perform 200 learning iterations
    for (int iter = 0; iter < Niter; iter++) safe_xgboost(XGBoosterUpdateOneIter(h_booster, iter, h_train));




    // free training sample
    safe_xgboost(XGDMatrixFree(h_train));
    free(train_IsSignal);
    free(train_weight);
    free(train_set);




    // set training sample2
    float* train_set2; int train_indicator2 = 0; float* train_weight2; float* train_IsSignal2;
    train_set2 = (float*)malloc(sizeof(float) * N_entry_train * Nvar);
    train_weight2 = (float*)malloc(sizeof(float) * N_entry_train);
    train_IsSignal2 = (float*)malloc(sizeof(float) * N_entry_train);

    // fill training sample2
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B2Knunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), train_set2, train_IsSignal2, train_weight2, true, &train_indicator2, Scale_Kplus * (3.0 / 7.0));
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B2Kstarnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), train_set2, train_IsSignal2, train_weight2, true, &train_indicator2, Scale_Kplusstar * (3.0 / 7.0));
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B2Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), train_set2, train_IsSignal2, train_weight2, true, &train_indicator2, Scale_Xsu_nonresonant * (3.0 / 7.0));
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B02K0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), train_set2, train_IsSignal2, train_weight2, true, &train_indicator2, Scale_K0 * (3.0 / 7.0));
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B02Kstar0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), train_set2, train_IsSignal2, train_weight2, true, &train_indicator2, Scale_K0star * (3.0 / 7.0));
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B02Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), train_set2, train_IsSignal2, train_weight2, true, &train_indicator2, Scale_Xsd_nonresonant * (3.0 / 7.0));
        }
    }
    {
        std::vector<string> names;
        load_files(CHG_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHG_input_train + std::string("/") + names.at(i)).c_str(), train_set2, train_IsSignal2, train_weight2, false, &train_indicator2, (0.3 / 1.0));
        }
    }
    {
        std::vector<string> names;
        load_files(MIX_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((MIX_input_train + std::string("/") + names.at(i)).c_str(), train_set2, train_IsSignal2, train_weight2, false, &train_indicator2, (0.3 / 1.0));
        }
    }
    {
        std::vector<string> names;
        load_files(UUBAR_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((UUBAR_input_train + std::string("/") + names.at(i)).c_str(), train_set2, train_IsSignal2, train_weight2, false, &train_indicator2, (0.3 / 1.1));
        }
    }
    {
        std::vector<string> names;
        load_files(DDBAR_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((DDBAR_input_train + std::string("/") + names.at(i)).c_str(), train_set2, train_IsSignal2, train_weight2, false, &train_indicator2, (0.3 / 1.1));
        }
    }
    {
        std::vector<string> names;
        load_files(SSBAR_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SSBAR_input_train + std::string("/") + names.at(i)).c_str(), train_set2, train_IsSignal2, train_weight2, false, &train_indicator2, (0.3 / 1.1));
        }
    }
    {
        std::vector<string> names;
        load_files(CHARM_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHARM_input_train + std::string("/") + names.at(i)).c_str(), train_set2, train_IsSignal2, train_weight2, false, &train_indicator2, (0.3 / 1.1));
        }
    }

    // convert to DMatrix
    DMatrixHandle h_train2;
    safe_xgboost(XGDMatrixCreateFromMat((float*)train_set2, N_entry_train, Nvar, FLT_MAX, &h_train2));



    // perform AUC calculation
    bst_ulong out_len;
    const float* f;
    safe_xgboost(XGBoosterPredict(h_booster, h_train2, 0, Niter, 0, &out_len, &f));
    AUC_train = CalculateAUC(out_len, f, train_IsSignal2, train_weight2);




    // free training sample2
    safe_xgboost(XGDMatrixFree(h_train2));
    free(train_IsSignal2);
    free(train_weight2);
    free(train_set2);




    // set testing sample
    float* test_set; int test_indicator = 0; float* test_weight; float* test_IsSignal;
    test_set = (float*)malloc(sizeof(float) * N_entry_test * Nvar);
    test_weight = (float*)malloc(sizeof(float) * N_entry_test);
    test_IsSignal = (float*)malloc(sizeof(float) * N_entry_test);

    // fill testing sample
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B2Knunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), test_set, test_IsSignal, test_weight, true, &test_indicator, Scale_Kplus);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B2Kstarnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), test_set, test_IsSignal, test_weight, true, &test_indicator, Scale_Kplusstar);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B2Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), test_set, test_IsSignal, test_weight, true, &test_indicator, Scale_Xsu_nonresonant);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B02K0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), test_set, test_IsSignal, test_weight, true, &test_indicator, Scale_K0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B02Kstar0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), test_set, test_IsSignal, test_weight, true, &test_indicator, Scale_K0star);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B02Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), test_set, test_IsSignal, test_weight, true, &test_indicator, Scale_Xsd_nonresonant);
        }
    }
    {
        std::vector<string> names;
        load_files(CHG_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHG_input_test + std::string("/") + names.at(i)).c_str(), test_set, test_IsSignal, test_weight, false, &test_indicator, 0.5);
        }
    }
    {
        std::vector<string> names;
        load_files(MIX_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((MIX_input_test + std::string("/") + names.at(i)).c_str(), test_set, test_IsSignal, test_weight, false, &test_indicator, 0.5);
        }
    }
    {
        std::vector<string> names;
        load_files(UUBAR_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((UUBAR_input_test + std::string("/") + names.at(i)).c_str(), test_set, test_IsSignal, test_weight, false, &test_indicator, 0.5);
        }
    }
    {
        std::vector<string> names;
        load_files(DDBAR_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((DDBAR_input_test + std::string("/") + names.at(i)).c_str(), test_set, test_IsSignal, test_weight, false, &test_indicator, 0.5);
        }
    }
    {
        std::vector<string> names;
        load_files(SSBAR_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SSBAR_input_test + std::string("/") + names.at(i)).c_str(), test_set, test_IsSignal, test_weight, false, &test_indicator, 0.5);
        }
    }
    {
        std::vector<string> names;
        load_files(CHARM_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHARM_input_test + std::string("/") + names.at(i)).c_str(), test_set, test_IsSignal, test_weight, false, &test_indicator, 0.5);
        }
    }

    // convert to DMatrix
    DMatrixHandle h_test;
    safe_xgboost(XGDMatrixCreateFromMat((float*)test_set, N_entry_test, Nvar, FLT_MAX, &h_test));


    // perform AUC calculation
    bst_ulong out_len;
    const float* f;
    safe_xgboost(XGBoosterPredict(h_booster, h_test, 0, Niter, 0, &out_len, &f));
    AUC_test = CalculateAUC(out_len, f, test_IsSignal, test_weight);




    // free testing sample
    safe_xgboost(XGDMatrixFree(h_test));
    free(test_IsSignal);
    free(test_weight);
    free(test_set);




    // save model
    safe_xgboost(XGBoosterSaveModel(h_booster, ("/home/belle2/junewoo/storage_b1/GridSearch/out/xgboost_" + model_name + ".model").c_str()));
    FILE* fp;
    fp = fopen(("/home/belle2/junewoo/storage_b1/GridSearch/out/Result_" + std::string(argv[1]) + "_" + std::string(argv[2]) + "_" + std::string(argv[3]) + "_" + std::string(argv[4]) + "_" + std::string(argv[5]) + "_" + std::string(argv[6])).c_str(), "w");
    fprintf(fp, "%s_%s_%s_%s_%s_%s_%s %lf %lf %d\n", argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], AUC_train, AUC_test);
    fclose(fp);




    // free xgboost internal structures
    safe_xgboost(XGBoosterFree(h_booster));

    return 0;
}
