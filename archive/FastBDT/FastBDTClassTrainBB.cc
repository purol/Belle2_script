
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
# define N_event_info 15
# define N_Upsilon_info 47
# define N_Bsig_info 78
# define N_Btag_info 7
# define N_decay 38 // five decay mode + others

# define Nvar 46
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
    tree_data->SetBranchAddress("missingMomentumOfEvent", &Vars[34]);
    tree_data->SetBranchAddress("missingMomentumOfEvent_theta", &Vars[35]);
    tree_data->SetBranchAddress("nParticlesInList__bomu__pl__clMuonFBDT_tight__bc", &Vars[36]);
    tree_data->SetBranchAddress("roeEextra__bocleanMask__bc", &Vars[37]);
    tree_data->SetBranchAddress("roePTheta__bocleanMask__bc", &Vars[38]);

    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dc_pValue_std", &Dc_pvalue_std);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_dr", &Dc_dr);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_dz", &Dc_dz);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_M", &Dc_M);
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

        if (Dc_chiProb > -0.5) {
            input_vars[Nvar - DvetoNvar + 0].push_back((float)Dc_pvalue_std);
            input_vars[Nvar - DvetoNvar + 1].push_back((float)Dc_dr);
            input_vars[Nvar - DvetoNvar + 2].push_back((float)Dc_dz);
            input_vars[Nvar - DvetoNvar + 3].push_back((float)Dc_M);
        }
        else {
            input_vars[Nvar - DvetoNvar + 0].push_back((float)0.0);
            input_vars[Nvar - DvetoNvar + 1].push_back((float)-1.0);
            input_vars[Nvar - DvetoNvar + 2].push_back((float)-100.0);
            input_vars[Nvar - DvetoNvar + 3].push_back((float)0.0);
        }
        if (D0_chiProb > -0.5) {
            input_vars[Nvar - DvetoNvar + 4].push_back((float)D0_chiProb);
            input_vars[Nvar - DvetoNvar + 5].push_back((float)D0_dz);
            input_vars[Nvar - DvetoNvar + 6].push_back((float)D0_M);
        }
        else {
            input_vars[Nvar - DvetoNvar + 4].push_back((float)0.0);
            input_vars[Nvar - DvetoNvar + 5].push_back((float)-100.0);
            input_vars[Nvar - DvetoNvar + 6].push_back((float)0.0);
        }

        IsSignal->push_back(tempissignal);

        weight->push_back(weight_N);

    }

    input_file->Close();
    //printf("==> Total %d events survive...\n", Nevt);
}

float GetScore(const FastBDT::Classifier& classifier, std::vector<std::vector<float>> InputVariables, std::vector<bool> IsSignal) {
    float sum = 0;
    for (unsigned int i = 0; i < IsSignal.size(); ++i) {
        std::vector<float> temp;
        for (int j = 0; j < Nvar; j++) temp.push_back(InputVariables.at(j).at(i));
        float p = classifier.predict(temp);
        sum += (static_cast<int>(IsSignal[i]) - p) * (static_cast<int>(IsSignal[i]) - p);
    }
    return sum / IsSignal.size();
}

float GetWeightedScore(const FastBDT::Classifier& classifier, std::vector<std::vector<float>> InputVariables, std::vector<bool> IsSignal, std::vector<float> weight) {
    float sum = 0;
    float N_sum = 0;
    for (unsigned int i = 0; i < IsSignal.size(); ++i) {
        std::vector<float> temp;
        for (int j = 0; j < Nvar; j++) temp.push_back(InputVariables.at(j).at(i));
        float p = classifier.predict(temp);
        sum += (static_cast<int>(IsSignal[i]) - p) * (static_cast<int>(IsSignal[i]) - p) * weight[i];
        N_sum = N_sum + weight[i];
    }
    return sum / N_sum;
}

void KSTest(const FastBDT::Classifier& classifier, std::vector<std::vector<float>> InputVariables_train, std::vector<bool> IsSignal_train, std::vector<float> weight_train, std::vector<std::vector<float>> InputVariables_test, std::vector<bool> IsSignal_test, std::vector<float> weight_test){

    TH1D* OBB_BKG_train = new TH1D("FastBDT_{2} BKG train", ";FastBDT_{2};", 40, 0, 1.0);
    TH1D* OBB_BKG_test = new TH1D("FastBDT_{2} BKG test", ";FastBDT_{2};", 40, 0, 1.0);
    TH1D* OBB_SIGNAL_train = new TH1D("FastBDT_{2} SIGNAL train", ";FastBDT_{2};", 40, 0, 1.0);
    TH1D* OBB_SIGNAL_test = new TH1D("FastBDT_{2} SIGNAL test", ";FastBDT_{2};", 40, 0, 1.0);

    for (unsigned int i = 0; i < IsSignal_train.size(); ++i) {
        if(IsSignal_train.at(i) == true) {
            std::vector<float> temp;
            for (int j = 0; j < Nvar; j++) temp.push_back(InputVariables_train.at(j).at(i));
            float p = classifier.predict(temp);
            OBB_SIGNAL_train->Fill(p, weight_train.at(i));
        }
        else if(IsSignal_train.at(i) == false) {
            std::vector<float> temp;
            for (int j = 0; j < Nvar; j++) temp.push_back(InputVariables_train.at(j).at(i));
            float p = classifier.predict(temp);
            OBB_BKG_train->Fill(p, weight_train.at(i));
        }
    }

    for (unsigned int i = 0; i < IsSignal_test.size(); ++i) {
        if(IsSignal_test.at(i) == true) {
            std::vector<float> temp;
            for (int j = 0; j < Nvar; j++) temp.push_back(InputVariables_test.at(j).at(i));
            float p = classifier.predict(temp);
            OBB_SIGNAL_test->Fill(p, weight_test.at(i));
        }
        else if(IsSignal_test.at(i) == false) {
            std::vector<float> temp;
            for (int j = 0; j < Nvar; j++) temp.push_back(InputVariables_test.at(j).at(i));
            float p = classifier.predict(temp);
            OBB_BKG_test->Fill(p, weight_test.at(i));
        }
    }

    OBB_SIGNAL_test->SetMarkerStyle(kFullCircle);
    OBB_SIGNAL_test->SetLineColor(kBlue);
    OBB_SIGNAL_test->SetMarkerColor(kBlue);
    OBB_SIGNAL_test->SetLineWidth(1);

    OBB_BKG_test->SetMarkerStyle(kFullCircle);
    OBB_BKG_test->SetLineColor(kRed);
    OBB_BKG_test->SetMarkerColor(kRed);
    OBB_BKG_test->SetLineWidth(1);

    OBB_SIGNAL_train->SetFillStyle(3004);
    OBB_SIGNAL_train->SetLineColor(kBlue);
    OBB_SIGNAL_train->SetFillColor(kBlue);

    OBB_BKG_train->SetFillStyle(3005);
    OBB_BKG_train->SetLineColor(kRed);
    OBB_BKG_train->SetFillColor(kRed);

    double factor = 1.0;

    // normalization
    OBB_SIGNAL_train->Scale(factor / OBB_SIGNAL_train->Integral(), "width");
    OBB_SIGNAL_test->Scale(factor / OBB_SIGNAL_test->Integral(), "width");
    OBB_BKG_train->Scale(factor / OBB_BKG_train->Integral(), "width");
    OBB_BKG_test->Scale(factor / OBB_BKG_test->Integral(), "width");

    double p_value_SIGNAL = OBB_SIGNAL_test->KolmogorovTest(OBB_SIGNAL_train);
    double p_value_BKG = OBB_BKG_test->KolmogorovTest(OBB_BKG_train);
    printf("p value SIGNAL: %lf\n", p_value_SIGNAL);
    printf("p value BKG: %lf\n", p_value_BKG);

    gStyle->SetOptStat(0);

    TCanvas* c_temp = new TCanvas("c", "", 600, 600); c_temp->cd();
    double OBB_BKG_train_max = OBB_BKG_train->GetMaximum(); double OBB_SIGNAL_train_max = OBB_SIGNAL_train->GetMaximum();
    if(OBB_BKG_train_max > OBB_SIGNAL_train_max) OBB_BKG_train->SetMaximum(1.05 * OBB_BKG_train_max);
    else OBB_BKG_train->SetMaximum(1.05 * OBB_SIGNAL_train_max);
    OBB_BKG_train->Draw("Hist"); OBB_SIGNAL_train->Draw("HistSAME");
    OBB_BKG_test->Draw("AP SAME"); OBB_SIGNAL_test->Draw("AP SAME");
    TLegend* legend = gPad->BuildLegend(0.9, 0.9, 0.6, 0.6); legend->SetFillStyle(0);
    c_temp->SaveAs("OBB_Plot.png"); c_temp->SetLogy(); c_temp->SaveAs("OBB_Plot_log.png");
    delete c_temp;
}

int main()
{
    var_names.push_back("aplanarity");
    var_names.push_back("Bsig_daughter_0_extraInfo_mydz");
    var_names.push_back("Btag_CleoConeCS_1");
    var_names.push_back("Btag_CleoConeCS_2");
    var_names.push_back("Btag_CleoConeCS_3");
    var_names.push_back("Btag_CleoConeCS_4");
    var_names.push_back("Btag_CleoConeCS_5");
    var_names.push_back("Btag_CleoConeCS_6");
    var_names.push_back("Btag_CleoConeCS_7");
    var_names.push_back("Btag_CleoConeCS_8");
    var_names.push_back("Btag_cosTBTO");
    var_names.push_back("Btag_cosTBz");
    var_names.push_back("Btag_deltaE");
    var_names.push_back("Btag_KSFWVariables_hso01");
    var_names.push_back("Btag_KSFWVariables_hso03");
    var_names.push_back("Btag_KSFWVariables_hso04");
    var_names.push_back("Btag_KSFWVariables_hso14");
    var_names.push_back("Btag_KSFWVariables_hso24");
    var_names.push_back("Btag_thrustBm");
    var_names.push_back("Btag_thrustOm");
    var_names.push_back("Btag_useCMSFrame_theta");
    var_names.push_back("cleoConeThrust0");
    var_names.push_back("cleoConeThrust1");
    var_names.push_back("cleoConeThrust2");
    var_names.push_back("cleoConeThrust3");
    var_names.push_back("cleoConeThrust4");
    var_names.push_back("cleoConeThrust5");
    var_names.push_back("cleoConeThrust6");
    var_names.push_back("cleoConeThrust7");
    var_names.push_back("cleoConeThrust8");
    var_names.push_back("harmonicMomentThrust2");
    var_names.push_back("harmonicMomentThrust3");
    var_names.push_back("harmonicMomentThrust4");
    var_names.push_back("Btag_extraInfo_SignalProbability");
    var_names.push_back("missingMomentumOfEvent");
    var_names.push_back("missingMomentumOfEvent_theta");
    var_names.push_back("nParticlesInList__bomu__pl__clMuonFBDT_tight__bc");
    var_names.push_back("roeEextra__bocleanMask__bc");
    var_names.push_back("roePTheta__bocleanMask__bc");

    var_names.push_back("Bsig_daughter_0_extraInfo_Dc_pValue_std");
    var_names.push_back("Bsig_daughter_0_extraInfo_Dcsimpleveto_dr");
    var_names.push_back("Bsig_daughter_0_extraInfo_Dcsimpleveto_dz");
    var_names.push_back("Bsig_daughter_0_extraInfo_Dcsimpleveto_M");
    var_names.push_back("Bsig_daughter_0_extraInfo_D0simpleveto_chiProb");
    var_names.push_back("Bsig_daughter_0_extraInfo_D0simpleveto_dz");
    var_names.push_back("Bsig_daughter_0_extraInfo_D0simpleveto_M");

    if (Nvar != var_names.size()) {
        printf("ERROR!\n");
        exit(1);
    }

    // define classifier and set options
    FastBDT::Classifier classifier;
    classifier.SetNTrees(1000);

    // define input of the classifier
    std::vector<std::vector<float>> InputVariables;
    std::vector<bool> IsSignal;
    std::vector<float> weight;

    // define input variables
    std::vector<float> input_vars[Nvar];

    // input file
    const char* SIGNAL_input_train = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SIGNAL_analysis/train_v001/final_output/DataFile";
    const char* CHG_input_train = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/CHG_analysis/train_v001/final_output/DataFile";
    const char* MIX_input_train = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/MIX_analysis/train_v001/final_output/DataFile";
    const char* UUBAR_input_train = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/UUBAR_analysis/train_v001/final_output/DataFile";
    const char* DDBAR_input_train = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/DDBAR_analysis/train_v001/final_output/DataFile";
    const char* SSBAR_input_train = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SSBAR_analysis/train_v001/final_output/DataFile";
    const char* CHARM_input_train = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/CHARM_analysis/train_v001/final_output/DataFile";

    const char* SIGNAL_input_test = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SIGNAL_analysis/test_v001/final_output/DataFile";
    const char* CHG_input_test = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/CHG_analysis/test_v001/final_output/DataFile";
    const char* MIX_input_test = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/MIX_analysis/test_v001/final_output/DataFile";
    const char* UUBAR_input_test = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/UUBAR_analysis/test_v001/final_output/DataFile";
    const char* DDBAR_input_test = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/DDBAR_analysis/test_v001/final_output/DataFile";
    const char* SSBAR_input_test = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SSBAR_analysis/test_v001/final_output/DataFile";
    const char* CHARM_input_test = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/CHARM_analysis/test_v001/final_output/DataFile";

    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B2Knunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 3.0 / 4.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B2Kstarnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 3.0 / 4.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B2Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 3.0 / 4.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B02K0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 3.0 / 4.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B02Kstar0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 3.0 / 4.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B02Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 3.0 / 4.0);
        }
    }
    {
        std::vector<string> names;
        load_files(CHG_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHG_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, (0.3 / 1.6));
        }
    }
    {
        std::vector<string> names;
        load_files(MIX_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((MIX_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, (0.3 / 1.6));
        }
    }
    {
        std::vector<string> names;
        load_files(UUBAR_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((UUBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, (0.3 / 1.7));
        }
    }
    {
        std::vector<string> names;
        load_files(DDBAR_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((DDBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, (0.3 / 1.7));
        }
    }
    {
        std::vector<string> names;
        load_files(SSBAR_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SSBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, (0.3 / 1.7));
        }
    }
    {
        std::vector<string> names;
        load_files(CHARM_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHARM_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, (0.3 / 1.7));
        }
    }

    // convert to double vector
    for (int i = 0; i < Nvar; i++) {
        InputVariables.push_back(input_vars[i]);
    }

    // fit
    classifier.fit(InputVariables, IsSignal, weight);

    // print result
    std::cout << "Train Score " << GetScore(classifier, InputVariables, IsSignal) << std::endl;
    std::cout << "Weighted Train Score " << GetWeightedScore(classifier, InputVariables, IsSignal, weight) << std::endl;

    std::fstream out_stream("BB.weightfile", std::ios_base::out | std::ios_base::trunc);
    out_stream << classifier << std::endl;
    out_stream.close();

    classifier.Print();

    std::map<unsigned int, double> rank;
    rank = classifier.GetVariableRanking();
    printf("Variable importance:\n");
    for (auto iter = rank.begin(); iter != rank.end(); iter++)
    {
        std::cout << "[" << iter->first << ", " << iter->second << "]" << " ";
    }
    printf("\n\n");
    printf("Variable importance for plot:\n");
    printf("[");
    for (auto iter = rank.begin(); iter != rank.end(); iter++)
    {
        int index = std::distance(rank.begin(), iter);
        std::cout << "(\'" << var_names.at(index) << "\'," << iter->second << ")";
        if (index == Nvar - 1) {}
        else {
            std::cout << "," << std::endl;
        }
    }
    printf("]");
    printf("\n\n");


    // print result of testing sample
    std::fstream in_stream("BB.weightfile", std::ios_base::in);
    FastBDT::Classifier classifier2(in_stream);

    std::vector<std::vector<float>> InputVariables2;
    std::vector<bool> IsSignal2;
    std::vector<float> weight2;

    std::vector<float> input_vars2[Nvar];

    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B2Knunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B2Kstarnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B2Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B02K0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B02Kstar0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B02Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(CHG_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHG_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false);
        }
    }
    {
        std::vector<string> names;
        load_files(MIX_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((MIX_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false);
        }
    }
    {
        std::vector<string> names;
        load_files(UUBAR_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((UUBAR_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false);
        }
    }
    {
        std::vector<string> names;
        load_files(DDBAR_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((DDBAR_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false);
        }
    }
    {
        std::vector<string> names;
        load_files(SSBAR_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SSBAR_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false);
        }
    }
    {
        std::vector<string> names;
        load_files(CHARM_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHARM_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false);
        }
    }

    // convert to double vector
    for (int i = 0; i < Nvar; i++) {
        InputVariables2.push_back(input_vars2[i]);
    }

    std::cout << "Test Score " << GetScore(classifier2, InputVariables2, IsSignal2) << std::endl;
    std::cout << "Weighted Test Score " << GetWeightedScore(classifier2, InputVariables2, IsSignal2, weight2) << std::endl;

    classifier2.Print();

    KSTest(classifier2, InputVariables, IsSignal, weight, InputVariables2, IsSignal2, weight2);

    return 0;
}
