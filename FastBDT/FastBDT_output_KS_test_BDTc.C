#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"

#include <string>
# include <vector>

using namespace RooFit ;

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

TH1F* FBDTc_MC_train = new TH1F("FastBDT_{c} MC train", ";FastBDT_{c};", 40, 0, 1.0);
TH1F* FBDTc_MC_test = new TH1F("FastBDT_{c} MC test", ";FastBDT_{c};", 40, 0, 1.0);

TH1F* FBDTc_MC_train_correction = new TH1F("FastBDT_{c} MC train with c", ";FastBDT_{c};", 40, 0, 1.0);
TH1F* FBDTc_MC_test_correction = new TH1F("FastBDT_{c} MC test with c", ";FastBDT_{c};", 40, 0, 1.0);

TH1F* FBDTc_data_train = new TH1F("FastBDT_{c} data train", ";FastBDT_{c};", 40, 0, 1.0);
TH1F* FBDTc_data_test = new TH1F("FastBDT_{c} data test", ";FastBDT_{c};", 40, 0, 1.0);

typedef struct _Nevt {
    double NevtwithoutCorrection;
    double NevtwithCorrection;
} Nevt;

double LetsFill(const char* filename, TH1F* FBDTc_hist, Nevt* nevt, double weight_var = 1.0) {
    float FBDT_var = 0;
    double NormFactor = 0;

    TFile* input_file = new TFile(filename, "read");

    TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
    tree_upsilon->SetBranchAddress("MVA_Continuum", &FBDT_var);

    for (unsigned int i = 0; i < tree_upsilon->GetEntries(); i++) {
        tree_upsilon->GetEntry(i);
        FBDTc_hist->Fill(FBDT_var, weight_var);

        if (FBDT_var > (5.0 / 6.0)) NormFactor = NormFactor + 5.0 * weight_var;
        else NormFactor = NormFactor + (FBDT_var / (1.0 - FBDT_var)) * weight_var;

        nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + weight_var;
        if (FBDT_var > (5.0 / 6.0)) nevt->NevtwithCorrection = nevt->NevtwithCorrection + 5.0 * weight_var;
        else nevt->NevtwithCorrection = nevt->NevtwithCorrection + (FBDT_var / (1.0 - FBDT_var)) * weight_var;

    }
    input_file->Close();

    return NormFactor;
}

double LetsFill(const char* filename, TH1F* FBDTc_hist, double weight_var = 1.0) {
    float FBDT_var = 0;
    double NormFactor = 0;

    TFile* input_file = new TFile(filename, "read");

    TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
    tree_upsilon->SetBranchAddress("MVA_Continuum", &FBDT_var);

    for (unsigned int i = 0; i < tree_upsilon->GetEntries(); i++) {
        tree_upsilon->GetEntry(i);
        FBDTc_hist->Fill(FBDT_var, weight_var);

        if (FBDT_var > (5.0 / 6.0)) NormFactor = NormFactor + 5.0 * weight_var;
        else NormFactor = NormFactor + (FBDT_var / (1.0 - FBDT_var)) * weight_var;

    }
    input_file->Close();

    return NormFactor;
}

double LetsFillwithCorrection(const char* filename, TH1F* FBDTc_hist, double weight_var = 1.0, double NormFactor = 1.0) {
    float FBDT_var = 0;

    TFile* input_file = new TFile(filename, "read");

    TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
    tree_upsilon->SetBranchAddress("MVA_Continuum", &FBDT_var);

    for (unsigned int i = 0; i < tree_upsilon->GetEntries(); i++) {
        tree_upsilon->GetEntry(i);

        if (FBDT_var > (5.0 / 6.0)) FBDTc_hist->Fill(FBDT_var, 5.0 * weight_var / NormFactor);
        else FBDTc_hist->Fill(FBDT_var, (FBDT_var / (1.0 - FBDT_var)) * weight_var / NormFactor);

    }
    input_file->Close();

    return NormFactor;
}

void FastBDT_output_KS_test_BDTc()
{
    Nevt nevt_train = { 0.0, 0.0 };
    Nevt nevt_test = { 0.0, 0.0 };

    // dirnames
    const char* side_data = "./temp_v000_validation/SIGNAL";

    const char* side_MC_CHG = "./temp_v000_validation/CHG";
    const char* side_MC_MIX = "./temp_v000_validation/MIX";
    const char* side_MC_UUBAR = "./temp_v000_validation/UUBAR";
    const char* side_MC_DDBAR = "./temp_v000_validation/DDBAR";
    const char* side_MC_SSBAR = "./temp_v000_validation/SSBAR";
    const char* side_MC_CHARM = "./temp_v000_validation/CHARM";


    // train sample
    {
        std::vector<string> names;
        load_files(side_data, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2 == 0) continue;
            LetsFill((side_data + std::string("/") + names.at(i)).c_str(), FBDTc_data_train, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_CHG, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2 == 0) continue;
            LetsFill((side_MC_CHG + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train, &nevt_train, 0.364436 / 0.9);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_MIX, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2 == 0) continue;
            LetsFill((side_MC_MIX + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train, &nevt_train, 0.364436 / 0.9);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_UUBAR, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2 == 0) continue;
            LetsFill((side_MC_UUBAR + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train, &nevt_train, 0.364436 / 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_DDBAR, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2 == 0) continue;
            LetsFill((side_MC_DDBAR + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train, &nevt_train, 0.364436 / 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_SSBAR, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2 == 0) continue;
            LetsFill((side_MC_SSBAR + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train, &nevt_train, 0.364436 / 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_CHARM, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2 == 0) continue;
            LetsFill((side_MC_CHARM + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train, &nevt_train, 0.364436 / 1.0);
        }
    }

    // test sample
    {
        std::vector<string> names;
        load_files(side_data, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2 == 1) continue;
            LetsFill((side_data + std::string("/") + names.at(i)).c_str(), FBDTc_data_test, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_CHG, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2 == 1) continue;
            LetsFill((side_MC_CHG + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test, &nevt_test, 0.364436 / 0.9);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_MIX, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2 == 1) continue;
            LetsFill((side_MC_MIX + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test, &nevt_test, 0.364436 / 0.9);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_UUBAR, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2 == 1) continue;
            LetsFill((side_MC_UUBAR + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test, &nevt_test, 0.364436 / 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_DDBAR, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2 == 1) continue;
            LetsFill((side_MC_DDBAR + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test, &nevt_test, 0.364436 / 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_SSBAR, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2 == 1) continue;
            LetsFill((side_MC_SSBAR + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test, &nevt_test, 0.364436 / 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_CHARM, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2 == 1) continue;
            LetsFill((side_MC_CHARM + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test, &nevt_test, 0.364436 / 1.0);
        }
    }


    double NormFactorForTrain = nevt_train.NevtwithoutCorrection / nevt_train.NevtwithCorrection;
    double NormFactorForTest = nevt_test.NevtwithoutCorrection / nevt_test.NevtwithCorrection;


    // Corrected MC train
    {
        std::vector<string> names;
        load_files(side_MC_CHG, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2 == 0) continue;
            LetsFillwithCorrection((side_MC_CHG + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train_correction, 0.364436 / 0.9, NormFactorForTrain);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_MIX, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2 == 0) continue;
            LetsFillwithCorrection((side_MC_MIX + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train_correction, 0.364436 / 0.9, NormFactorForTrain);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_UUBAR, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2 == 0) continue;
            LetsFillwithCorrection((side_MC_UUBAR + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train_correction, 0.364436 / 1.0, NormFactorForTrain);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_DDBAR, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2 == 0) continue;
            LetsFillwithCorrection((side_MC_DDBAR + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train_correction, 0.364436 / 1.0, NormFactorForTrain);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_SSBAR, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2 == 0) continue;
            LetsFillwithCorrection((side_MC_SSBAR + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train_correction, 0.364436 / 1.0, NormFactorForTrain);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_CHARM, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2 == 0) continue;
            LetsFillwithCorrection((side_MC_CHARM + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train_correction, 0.364436 / 1.0, NormFactorForTrain);
        }
    }


    // Corrected MC test
    {
        std::vector<string> names;
        load_files(side_MC_CHG, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2  == 1) continue;
            LetsFillwithCorrection((side_MC_CHG + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test_correction, 0.364436 / 0.9, NormFactorForTest);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_MIX, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2  == 1) continue;
            LetsFillwithCorrection((side_MC_MIX + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test_correction, 0.364436 / 0.9, NormFactorForTest);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_UUBAR, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2  == 1) continue;
            LetsFillwithCorrection((side_MC_UUBAR + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test_correction, 0.364436 / 1.0, NormFactorForTest);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_DDBAR, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2  == 1) continue;
            LetsFillwithCorrection((side_MC_DDBAR + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test_correction, 0.364436 / 1.0, NormFactorForTest);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_SSBAR, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2  == 1) continue;
            LetsFillwithCorrection((side_MC_SSBAR + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test_correction, 0.364436 / 1.0, NormFactorForTest);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_CHARM, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 2 == 1) continue;
            LetsFillwithCorrection((side_MC_CHARM + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test_correction, 0.364436 / 1.0, NormFactorForTest);
        }
    }



    // normalization
    double factor = 1.0;
    FBDTc_MC_train->Scale(factor / FBDTc_MC_train->Integral(), "width");
    FBDTc_MC_test->Scale(factor / FBDTc_MC_test->Integral(), "width");
    FBDTc_MC_train_correction->Scale(factor / FBDTc_MC_train_correction->Integral(), "width");
    FBDTc_MC_test_correction->Scale(factor / FBDTc_MC_test_correction->Integral(), "width");
    FBDTc_data_train->Scale(factor / FBDTc_data_train->Integral(), "width");
    FBDTc_data_test->Scale(factor / FBDTc_data_test->Integral(), "width");


    // set color (MC: kRed, data: kBlue)
    // solid line: test, histogram: traing
    FBDTc_data_test->SetMarkerStyle(kFullCircle);
    FBDTc_data_test->SetLineColor(kBlue);
    FBDTc_data_test->SetMarkerColor(kBlue);
    FBDTc_data_test->SetLineWidth(1);

    FBDTc_MC_test->SetMarkerStyle(kFullCircle);
    FBDTc_MC_test->SetLineColor(kRed);
    FBDTc_MC_test->SetMarkerColor(kRed);
    FBDTc_MC_test->SetLineWidth(1);

    FBDTc_data_train->SetFillStyle(3004);
    FBDTc_data_train->SetLineColor(kBlue);
    FBDTc_data_train->SetFillColor(kBlue);

    FBDTc_MC_train->SetFillStyle(3005);
    FBDTc_MC_train->SetLineColor(kRed);
    FBDTc_MC_train->SetFillColor(kRed);

    FBDTc_MC_test_correction->SetMarkerStyle(kFullCircle);
    FBDTc_MC_test_correction->SetLineColor(kRed);
    FBDTc_MC_test_correction->SetMarkerColor(kRed);
    FBDTc_MC_test_correction->SetLineWidth(1);

    FBDTc_MC_train_correction->SetFillStyle(3005);
    FBDTc_MC_train_correction->SetLineColor(kRed);
    FBDTc_MC_train_correction->SetFillColor(kRed);

    gStyle->SetOptStat(0);

    TCanvas* c_temp = new TCanvas("c", "", 600, 600); c_temp->cd(); gPad->SetLogy();
    //double OBB_BKG_train_max = OBB_BKG_train->GetMaximum(); double OBB_SIGNAL_train_max = OBB_SIGNAL_train->GetMaximum();
    //if(OBB_BKG_train_max > OBB_SIGNAL_train_max) OBB_BKG_train->SetMaximum(1.05 * OBB_BKG_train_max);
    //else OBB_BKG_train->SetMaximum(1.05 * OBB_SIGNAL_train_max);
    FBDTc_MC_train->Draw("Hist"); FBDTc_data_train->Draw("HistSAME");
    FBDTc_MC_test->Draw("AP SAME"); FBDTc_data_test->Draw("AP SAME");
    gPad->BuildLegend(0.9, 0.9, 0.6, 0.6);
    c_temp->SaveAs("BDTc_BeforeCorrection.png");

    TCanvas* c_temp_2 = new TCanvas("c2", "", 600, 600); c_temp_2->cd(); gPad->SetLogy();
    //double Oqq_BKG_train_max = Oqq_BKG_train->GetMaximum(); double Oqq_SIGNAL_train_max = Oqq_SIGNAL_train->GetMaximum();
    //if(Oqq_BKG_train_max > Oqq_SIGNAL_train_max) Oqq_BKG_train->SetMaximum(1.05 * Oqq_BKG_train_max);
    //else Oqq_BKG_train->SetMaximum(1.05 * Oqq_SIGNAL_train_max);
    FBDTc_MC_train_correction->Draw("Hist"); FBDTc_data_train->Draw("HistSAME");
    FBDTc_MC_test_correction->Draw("AP SAME"); FBDTc_data_test->Draw("AP SAME");
    gPad->BuildLegend(0.9, 0.9, 0.6, 0.6);
    c_temp_2->SaveAs("BDTc_AfterCorrection.png");

}
