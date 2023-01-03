
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

#include "Classifier.h"

# define N_Needed_info 37
//# define N_event_info 15
# define N_Upsilon_info 57
# define N_Bsig_info 81
# define N_Btag_info 9
# define N_decay 38 // five decay mode + others
# define N_decay_nparticles 3 // # of nu_e, B, B0
# define N_decay_syst_ff 7 // helicity angle + q2

# define Nvar 44
# define DvetoNvar 4

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

/*
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
*/

// SIGNAL MC sample number
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

// scale factor for SIGNAL MC sample (364.436 - 2.763 = 361.673/fb)
# define Scale_Kplus_train (0.361673 * N_Kplus_nunubar_1invab/N_Kplus_train)
# define Scale_Kplusstar_train (0.361673 * N_Kplusstar_nunubar_1invab/N_Kplusstar_train)
# define Scale_Xsu_nonresonant_train (0.361673 * N_Xsu_nonresonant_nunubar_1invab/N_Xsu_nonresonant_train)
# define Scale_K0_train (0.361673 * N_K0_nunubar_1invab/N_K0_train)
# define Scale_K0star_train (0.361673 * N_K0star_nunubar_1invab/N_K0star_train)
# define Scale_Xsd_nonresonant_train (0.361673 * N_Xsd_nunubar_1invab/N_Xsd_nonresonant_train)
# define Scale_Kplus_test (0.361673 * N_Kplus_nunubar_1invab/N_Kplus_test)
# define Scale_Kplusstar_test (0.361673 * N_Kplusstar_nunubar_1invab/N_Kplusstar_test)
# define Scale_Xsu_nonresonant_test (0.361673 * N_Xsu_nonresonant_nunubar_1invab/N_Xsu_nonresonant_test)
# define Scale_K0_test (0.361673 * N_K0_nunubar_1invab/N_K0_test)
# define Scale_K0star_test (0.361673 * N_K0star_nunubar_1invab/N_K0star_test)
# define Scale_Xsd_nonresonant_test (0.361673 * N_Xsd_nunubar_1invab/N_Xsd_nonresonant_test)

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

// new scale factor for BKG MC sample with additional 1/ab (364.436 - 2.763 = 361.673/fb)
# define Scale_CHG_train (0.361673/((N_CHG_train/(N_CHG_train + N_CHG_test))*0.8+1.0))
# define Scale_MIX_train (0.361673/((N_MIX_train/(N_MIX_train + N_MIX_test))*0.8+1.0))
# define Scale_UUBAR_train (0.361673/((N_UUBAR_train/(N_UUBAR_train + N_UUBAR_test))*1.0))
# define Scale_DDBAR_train (0.361673/((N_DDBAR_train/(N_DDBAR_train + N_DDBAR_test))*1.0))
# define Scale_SSBAR_train (0.361673/((N_SSBAR_train/(N_SSBAR_train + N_SSBAR_test))*1.0))
# define Scale_CHARM_train (0.361673/((N_CHARM_train/(N_CHARM_train + N_CHARM_test))*1.0))
# define Scale_CHG_test (0.361673/((N_CHG_test/(N_CHG_train + N_CHG_test))*0.8+1.0))
# define Scale_MIX_test (0.361673/((N_MIX_test/(N_MIX_train + N_MIX_test))*0.8+1.0))
# define Scale_UUBAR_test (0.361673/((N_UUBAR_test/(N_UUBAR_train + N_UUBAR_test))*1.0))
# define Scale_DDBAR_test (0.361673/((N_DDBAR_test/(N_DDBAR_train + N_DDBAR_test))*1.0))
# define Scale_SSBAR_test (0.361673/((N_SSBAR_test/(N_SSBAR_train + N_SSBAR_test))*1.0))
# define Scale_CHARM_test (0.361673/((N_CHARM_test/(N_CHARM_train + N_CHARM_test))*1.0))

using std::string;

std::vector<std::string> var_names;

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

void FillVariables(const char * filename, std::vector<float> input_vars[Nvar], std::vector<bool>* IsSignal, std::vector<float>* weight, bool tempissignal, float weight_N = 1.0) {
    TFile* input_file = new TFile(filename, "read");

    TTree* tree_data = (TTree*)input_file->Get("data");

    double Vars[Nvar];
    int flag;

    double Dc_chiProb; // 0.0
    double Dc_pvalue_med;
    double Dc_pvalue_std; // 0.0
    double Dc_dr; // -1.0
    double Dc_dz; // -100.0
    double Dc_M; // 0.0
    double D0_chiProb;
    double D0_pvalue_med;
    double D0_pvalue_std;
    double D0_dr;
    double D0_dz;
    double D0_M;

    int Decay_Kplus = -1;
    int Decay_Kplusstar_ch1 = -1;
    int Decay_Kplusstar_ch2 = -1;
    int Decay_K0 = -1;
    int Decay_K0star_ch1 = -1;
    int Decay_K0star_ch2 = -1;

    double Mxs = -1;
    double Pcms = -1;

    tree_data->SetBranchAddress("aplanarity", &Vars[0]);
    tree_data->SetBranchAddress("Bsig_cosTBTO", &Vars[1]);
    tree_data->SetBranchAddress("Bsig_KSFWVariables_hoo1", &Vars[2]);
    tree_data->SetBranchAddress("Bsig_KSFWVariables_hoo3", &Vars[3]);
    tree_data->SetBranchAddress("Bsig_KSFWVariables_hso01", &Vars[4]);
    tree_data->SetBranchAddress("Bsig_KSFWVariables_hso03", &Vars[5]);
    tree_data->SetBranchAddress("Bsig_KSFWVariables_hso04", &Vars[6]);
    tree_data->SetBranchAddress("Bsig_useCMSFrame_p", &Vars[7]);
    tree_data->SetBranchAddress("Btag_chiProb", &Vars[8]);
    tree_data->SetBranchAddress("Btag_CleoConeCS_1", &Vars[9]);
    tree_data->SetBranchAddress("Btag_CleoConeCS_2", &Vars[10]);
    tree_data->SetBranchAddress("Btag_CleoConeCS_3", &Vars[11]);
    tree_data->SetBranchAddress("Btag_CleoConeCS_4", &Vars[12]);
    tree_data->SetBranchAddress("Btag_CleoConeCS_5", &Vars[13]);
    tree_data->SetBranchAddress("Btag_cosTBTO", &Vars[14]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hoo1", &Vars[15]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hoo2", &Vars[16]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hoo3", &Vars[17]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hoo4", &Vars[18]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hso00", &Vars[19]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hso01", &Vars[20]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hso02", &Vars[21]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hso22", &Vars[22]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hso24", &Vars[23]);
    tree_data->SetBranchAddress("Btag_useCMSFrame_theta", &Vars[24]);
    tree_data->SetBranchAddress("cleoConeThrust1", &Vars[25]);
    tree_data->SetBranchAddress("cleoConeThrust2", &Vars[26]);
    tree_data->SetBranchAddress("cleoConeThrust3", &Vars[27]);
    tree_data->SetBranchAddress("extraInfo__boEeclv133__bc", &Vars[28]);
    tree_data->SetBranchAddress("foxWolframR3", &Vars[29]);
    tree_data->SetBranchAddress("foxWolframR4", &Vars[30]);
    tree_data->SetBranchAddress("harmonicMomentThrust1", &Vars[31]);
    tree_data->SetBranchAddress("harmonicMomentThrust2", &Vars[32]);
    tree_data->SetBranchAddress("harmonicMomentThrust3", &Vars[33]);
    tree_data->SetBranchAddress("harmonicMomentThrust4", &Vars[34]);
    tree_data->SetBranchAddress("Btag_extraInfo_SignalProbability", &Vars[35]);
    tree_data->SetBranchAddress("missingMomentumOfEvent", &Vars[36]);
    tree_data->SetBranchAddress("missingMomentumOfEvent_theta", &Vars[37]);
    tree_data->SetBranchAddress("nParticlesInList__bomu__pl__clMuonFBDT_tight__bc", &Vars[38]);
    tree_data->SetBranchAddress("roePTheta__bocleanMask__bc", &Vars[39]);

    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dc_pValue_med", &Dc_pvalue_med);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dc_pValue_std", &Dc_pvalue_std);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_chiProb", &Dc_chiProb);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_dz", &Dc_dz);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_M", &Dc_M);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_D0_pValue_med", &D0_pvalue_med);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_D0_pValue_std", &D0_pvalue_std);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_D0simpleveto_chiProb", &D0_chiProb);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_D0simpleveto_dz", &D0_dz);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_D0simpleveto_M", &D0_M);

    tree_data->SetBranchAddress("flag", &flag);

    if(tempissignal){
        tree_data->SetBranchAddress("nParticlesInList__boB__pl__clKcharge_total__bc", &Decay_Kplus);
        tree_data->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &Decay_Kplusstar_ch1);
        tree_data->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &Decay_Kplusstar_ch2);
        tree_data->SetBranchAddress("nParticlesInList__boB0__clKneutral_total__bc", &Decay_K0);
        tree_data->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &Decay_K0star_ch1);
        tree_data->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &Decay_K0star_ch2);
    }
    tree_data->SetBranchAddress("Bsig_M", &Mxs);

    int Nevt = 0;
    //printf("%lld entries...\n", tree_data->GetEntries());
    for (unsigned int j = 0; j < tree_data->GetEntries(); j++) { // Fill
        tree_data->GetEntry(j);
        if(tempissignal == true && (Decay_Kplus > 0.5 || Decay_Kplusstar_ch1 > 0.5 || Decay_Kplusstar_ch2 > 0.5 || Decay_K0 > 0.5 || Decay_K0star_ch1 > 0.5 || Decay_K0star_ch2 > 0.5) && Mxs > 1.1) continue;
        else if(tempissignal == true && (Decay_Kplus < 0.5 && Decay_Kplusstar_ch1 < 0.5 && Decay_Kplusstar_ch2 < 0.5 && Decay_K0 < 0.5 &&  Decay_K0star_ch1 < 0.5 && Decay_K0star_ch2 < 0.5) && Mxs < 1.1) continue;
        Nevt++;

        for (unsigned int k = 0; k < Nvar - DvetoNvar; k++) input_vars[k].push_back((float) Vars[k]); 

        if(Dc_chiProb > -0.5){
            input_vars[Nvar - DvetoNvar + 0].push_back((float)Dc_pvalue_std);
            input_vars[Nvar - DvetoNvar + 1].push_back((float)Dc_M);
        }
        else {
            input_vars[Nvar - DvetoNvar + 0].push_back((float)0.0);
            input_vars[Nvar - DvetoNvar + 1].push_back((float) 0.0);
        }
        if(D0_chiProb > -0.5){
            input_vars[Nvar - DvetoNvar + 2].push_back((float)D0_pvalue_std);
            input_vars[Nvar - DvetoNvar + 3].push_back((float) D0_M);
        }
        else {
            input_vars[Nvar - DvetoNvar + 2].push_back((float) 0.0);
            input_vars[Nvar - DvetoNvar + 3].push_back((float) 0.0);
        }

        IsSignal->push_back(tempissignal);

        weight->push_back(weight_N);

    }

    input_file->Close();
    //printf("==> Total %d events survive...\n", Nevt);
}

double PrintMaximumFOM(const FastBDT::Classifier& classifier, std::vector<std::vector<float>> InputVariables, std::vector<bool> IsSignal, std::vector<float> weight) {
    const int step = 100;
    double FOM_max = -1;

    for (int i = 0; i < step; i++) {
        float value = ((float)i) / ((float)step);
        double NBKG = 0;
        double NSIG = 0;

        for (unsigned int i = 0; i < IsSignal.size(); ++i) {
            std::vector<float> temp;
            for (int j = 0; j < Nvar; j++) temp.push_back(InputVariables.at(j).at(i));
            float p = classifier.predict(temp);
            if (p > value) {
                if(IsSignal[i]) NSIG = NSIG + weight[i];
                else NBKG = NBKG + weight[i];
            }
        }

        double FOM = NSIG / std::sqrt(NBKG + NSIG);
        if (FOM > FOM_max) FOM_max = FOM;

    }

    return FOM_max;
}

double PrintAUC(const FastBDT::Classifier& classifier, std::vector<std::vector<float>> InputVariables, std::vector<bool> IsSignal, std::vector<float> weight) {
    const int step = 100;
    double AUC = 0;
    double NBKG_total = 0;
    double NSIG_total = 0;
    std::vector<double> TPRs;
    std::vector<double> FPRs;

    for (unsigned int i = 0; i < IsSignal.size(); ++i) {
        if (IsSignal[i]) NSIG_total = NSIG_total + weight[i];
        else NBKG_total = NBKG_total + weight[i];
    }

    for (int i = 0; i < step; i++) {
        float value = ((float)i) / ((float)step);
        double NBKG = 0;
        double NSIG = 0;

        for (unsigned int i = 0; i < IsSignal.size(); ++i) {
            std::vector<float> temp;
            for (int j = 0; j < Nvar; j++) temp.push_back(InputVariables.at(j).at(i));
            float p = classifier.predict(temp);
            if (p >= value) {
                if (IsSignal[i]) NSIG = NSIG + weight[i];
                else NBKG = NBKG + weight[i];
            }
        }

        double TPR = NSIG / NSIG_total;
        double FPR = NBKG / NBKG_total;

        TPRs.push_back(TPR);
        FPRs.push_back(FPR);
    }

    for (unsigned int i = 0; i < TPRs.size(); ++i) {
        if ( i != TPRs.size() - 1) {
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

double PrintAVG(const FastBDT::Classifier& classifier, std::vector<std::vector<float>> InputVariables, std::vector<bool> IsSignal, std::vector<float> weight, bool SelectSignal) {
    double NBKG_total = 0;
    double NSIG_total = 0;

    double NBKG_AVG = 0;
    double NSIG_AVG = 0;

    for (unsigned int i = 0; i < IsSignal.size(); ++i) {
        if (IsSignal[i]) NSIG_total = NSIG_total + weight[i];
        else NBKG_total = NBKG_total + weight[i];
    }

    for (unsigned int i = 0; i < IsSignal.size(); ++i) {
        std::vector<float> temp;
        for (int j = 0; j < Nvar; j++) temp.push_back(InputVariables.at(j).at(i));
        float p = classifier.predict(temp);

        if (IsSignal[i]) NSIG_AVG = NSIG_AVG + p * weight[i];
        else NBKG_AVG = NBKG_AVG + p * weight[i];

    }

    NSIG_AVG = NSIG_AVG / NSIG_total;
    NBKG_AVG = NBKG_AVG / NBKG_total;

    if (SelectSignal) return NSIG_AVG;
    else return NBKG_AVG;
}

int main(int argc, char* argv[])
{
    // grid search
    // unsigned int nTrees[5] = { 100, 500, 1000, 1500, 2000 };  default is 100
    // unsigned int depth[3] = { 2, 3, 4 };  default is 3 
    // double shrinkage[4] = { 0.05, 0.1, 0.15, 0.2 };  default is 0.1
    // double subsample[5] = { 0.3, 0.4, 0.5, 0.6, 0.7 };  default is 0.5
    // unsigned int binning[4] = { 6, 7, 8, 9 };  default is 2^8 bins per feature

    unsigned int nTrees = (unsigned int)atoi(argv[1]);
    unsigned int depth = (unsigned int)atoi(argv[2]);
    double shrinkage = atof(argv[3]);
    double subsample = atof(argv[4]);
    unsigned int binning_num = (unsigned int)atoi(argv[5]);

    // set classifier option
    FastBDT::Classifier classifier;
    classifier.SetNTrees(nTrees);
    classifier.SetDepth(depth);
    classifier.SetShrinkage(shrinkage);
    classifier.SetSubsample(subsample);
    std::vector<unsigned int> binning(Nvar, binning_num); classifier.SetBinning(binning);



    // input file
    const char* SIGNAL_input_train = "/home/belle2/junewoo/storage_b1/GridSearch/SIGNAL_analysis/train_v000/final_output/DataFile";
    const char* CHG_input_train = "/home/belle2/junewoo/storage_b1/GridSearch/CHG_analysis/train_v000/final_output/DataFile";
    const char* MIX_input_train = "/home/belle2/junewoo/storage_b1/GridSearch/MIX_analysis/train_v000/final_output/DataFile";
    const char* UUBAR_input_train = "/home/belle2/junewoo/storage_b1/GridSearch/UUBAR_analysis/train_v000/final_output/DataFile";
    const char* DDBAR_input_train = "/home/belle2/junewoo/storage_b1/GridSearch/DDBAR_analysis/train_v000/final_output/DataFile";
    const char* SSBAR_input_train = "/home/belle2/junewoo/storage_b1/GridSearch/SSBAR_analysis/train_v000/final_output/DataFile";
    const char* CHARM_input_train = "/home/belle2/junewoo/storage_b1/GridSearch/CHARM_analysis/train_v000/final_output/DataFile";

    const char* SIGNAL_input_test = "/home/belle2/junewoo/storage_b1/GridSearch/SIGNAL_analysis/test_v000/final_output/DataFile";
    const char* CHG_input_test = "/home/belle2/junewoo/storage_b1/GridSearch/CHG_analysis/test_v000/final_output/DataFile";
    const char* MIX_input_test = "/home/belle2/junewoo/storage_b1/GridSearch/MIX_analysis/test_v000/final_output/DataFile";
    const char* UUBAR_input_test = "/home/belle2/junewoo/storage_b1/GridSearch/UUBAR_analysis/test_v000/final_output/DataFile";
    const char* DDBAR_input_test = "/home/belle2/junewoo/storage_b1/GridSearch/DDBAR_analysis/test_v000/final_output/DataFile";
    const char* SSBAR_input_test = "/home/belle2/junewoo/storage_b1/GridSearch/SSBAR_analysis/test_v000/final_output/DataFile";
    const char* CHARM_input_test = "/home/belle2/junewoo/storage_b1/GridSearch/CHARM_analysis/test_v000/final_output/DataFile";



    // define input of the classifier
    std::vector<std::vector<float>> InputVariables;
    std::vector<bool> IsSignal;
    std::vector<float> weight;

    // define input variables
    std::vector<float> input_vars[Nvar];

    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B2Knunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B2Kstarnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B2Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B02K0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B02Kstar0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B02Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(CHG_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHG_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, Scale_CHG_train);
        }
    }
    {
        std::vector<string> names;
        load_files(MIX_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((MIX_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, Scale_MIX_train);
        }
    }
    {
        std::vector<string> names;
        load_files(UUBAR_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((UUBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, Scale_UUBAR_train);
        }
    }
    {
        std::vector<string> names;
        load_files(DDBAR_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((DDBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, Scale_DDBAR_train);
        }
    }
    {
        std::vector<string> names;
        load_files(SSBAR_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SSBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, Scale_SSBAR_train);
        }
    }
    {
        std::vector<string> names;
        load_files(CHARM_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHARM_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, Scale_CHARM_train);
        }
    }

    // convert to double vector
    for (int i = 0; i < Nvar; i++) {
        InputVariables.push_back(input_vars[i]);
    }



    // fit
    classifier.fit(InputVariables, IsSignal, weight);

    // clear vector to save memory
    for (unsigned int i = 0; i < InputVariables.size(); ++i) std::vector<float>().swap(InputVariables.at(i));
    std::vector<std::vector<float>>().swap(InputVariables);
    std::vector<bool>().swap(IsSignal);
    std::vector<float>().swap(weight);
    for (int i = 0; i < Nvar; i++) std::vector<float>().swap(input_vars[i]);



    // test sample
    std::vector<std::vector<float>> InputVariables2;
    std::vector<bool> IsSignal2;
    std::vector<float> weight2;

    std::vector<float> input_vars2[Nvar];

    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B2Knunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, Scale_Kplus_test);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B2Kstarnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, Scale_Kplusstar_test);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B2Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, Scale_Xsu_nonresonant_test);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B02K0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, Scale_K0_test);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B02Kstar0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, Scale_K0star_test);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B02Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, Scale_Xsd_nonresonant_test);
        }
    }
    {
        std::vector<string> names;
        load_files(CHG_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((CHG_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, Scale_CHG_test);
        }
    }
    {
        std::vector<string> names;
        load_files(MIX_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((MIX_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, Scale_MIX_test);
        }
    }
    {
        std::vector<string> names;
        load_files(UUBAR_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((UUBAR_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, Scale_UUBAR_test);
        }
    }
    {
        std::vector<string> names;
        load_files(DDBAR_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((DDBAR_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, Scale_DDBAR_test);
        }
    }
    {
        std::vector<string> names;
        load_files(SSBAR_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SSBAR_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, Scale_SSBAR_test);
        }
    }
    {
        std::vector<string> names;
        load_files(CHARM_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((CHARM_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, Scale_CHARM_test);
        }
    }

    // convert to double vector
    for (int i = 0; i < Nvar; i++) {
        InputVariables2.push_back(input_vars2[i]);
    }



    // get FOM for testing sample
    // double test_FOM = PrintMaximumFOM(classifier, InputVariables2, IsSignal2, weight2);

    // get AUC for testing sample
    double test_AUC = PrintAUC(classifier, InputVariables2, IsSignal2, weight2);

    // clear vector to save memory
    for (unsigned int i = 0; i < InputVariables2.size(); ++i) std::vector<float>().swap(InputVariables2.at(i));
    std::vector<std::vector<float>>().swap(InputVariables2);
    std::vector<bool>().swap(IsSignal2);
    std::vector<float>().swap(weight2);
    for (int i = 0; i < Nvar; i++) std::vector<float>().swap(input_vars2[i]);



    // train sample with proper weight
    std::vector<std::vector<float>> InputVariables3;
    std::vector<bool> IsSignal3;
    std::vector<float> weight3;

    std::vector<float> input_vars3[Nvar];

    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B2Knunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, Scale_Kplus_train);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B2Kstarnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, Scale_Kplusstar_train);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B2Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, Scale_Xsu_nonresonant_train);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B02K0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, Scale_K0_train);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B02Kstar0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, Scale_K0star_train);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B02Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, Scale_Xsd_nonresonant_train);
        }
    }
    {
        std::vector<string> names;
        load_files(CHG_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHG_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, Scale_CHG_train);
        }
    }
    {
        std::vector<string> names;
        load_files(MIX_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((MIX_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, Scale_MIX_train);
        }
    }
    {
        std::vector<string> names;
        load_files(UUBAR_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((UUBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, Scale_UUBAR_train);
        }
    }
    {
        std::vector<string> names;
        load_files(DDBAR_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((DDBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, Scale_DDBAR_train);
        }
    }
    {
        std::vector<string> names;
        load_files(SSBAR_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SSBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, Scale_SSBAR_train);
        }
    }
    {
        std::vector<string> names;
        load_files(CHARM_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHARM_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, Scale_CHARM_train);
        }
    }

    // convert to double vector
    for (int i = 0; i < Nvar; i++) {
        InputVariables3.push_back(input_vars3[i]);
    }



    // get FOM for training sample
    // double train_FOM = PrintMaximumFOM(classifier, InputVariables3, IsSignal3, weight3);

    // get AUC for training sample
    double train_AUC = PrintAUC(classifier, InputVariables3, IsSignal3, weight3);

    // clear vector to save memory
    for (unsigned int i = 0; i < InputVariables3.size(); ++i) std::vector<float>().swap(InputVariables3.at(i));
    std::vector<std::vector<float>>().swap(InputVariables3);
    std::vector<bool>().swap(IsSignal3);
    std::vector<float>().swap(weight3);
    for (int i = 0; i < Nvar; i++) std::vector<float>().swap(input_vars3[i]);



    printf("%u_%u_%lf_%lf_%u %lf %lf\n", nTrees, depth, shrinkage, subsample, binning_num, train_AUC, test_AUC);

    FILE* fp;
    fp = fopen(("/home/belle2/junewoo/storage_b1/GridSearch/out/Result_" + std::string(argv[1]) + "_" + std::string(argv[2]) + "_" + std::string(argv[3]) + "_" + std::string(argv[4]) + "_" + std::string(argv[5])).c_str(), "w");
    fprintf(fp, "%u_%u_%lf_%lf_%u %lf %lf\n", nTrees, depth, shrinkage, subsample, binning_num, train_AUC, test_AUC);
    fclose(fp);



    // save model
    std::fstream out_stream(("/home/belle2/junewoo/storage_b1/GridSearch/out/classifier_" + std::string(argv[1]) + "_" + std::string(argv[2]) + "_" + std::string(argv[3]) + "_" + std::string(argv[4]) + "_" + std::string(argv[5])+".weightfile").c_str(), std::ios_base::out | std::ios_base::trunc);
    out_stream << classifier << std::endl;
    out_stream.close();

    return 0;
}
