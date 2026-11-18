#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"
#include "RooPlot.h"

#include <string>
#include <algorithm>
# include <vector>

using namespace RooFit ;

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

TH1F* FBDTc_MC_train = new TH1F("MC train", ";BDT_{c};", 40, 0.0, 1.0);
TH1F* FBDTc_MC_test = new TH1F("MC test", ";BDT_{c};", 40, 0.0, 1.0);

TH1F* FBDTc_MC_train_correction = new TH1F("MC train with correction", ";BDT_{c};", 40, 0.0, 1.0);
TH1F* FBDTc_MC_test_correction = new TH1F("MC test with correction", ";BDT_{c};", 40, 0.0, 1.0);

TH1F* FBDTc_data_train = new TH1F("data train", ";BDT_{c};", 40, 0.0, 1.0);
TH1F* FBDTc_data_test = new TH1F("data test", ";BDT_{c};", 40, 0.0, 1.0);

typedef struct _Nevt {
    double NevtwithoutCorrection;
    double NevtwithCorrection;
    double BDTcSum;
} Nevt;

double BDTcToWeight(double BDTc) {

    if ( BDTc > (5.0 / 6.0) ) return std::sqrt(5.0);
    else return std::sqrt( BDTc / (1.0 - BDTc) );

    //return 2*BDTc;

}

double LetsFill(const char* filename, TH1F* FBDTc_hist, Nevt* nevt, double weight_var = 1.0) {
    float FBDT_var = 0;
    double NormFactor = 0;

    TFile* input_file = new TFile(filename, "read");

    TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
    tree_upsilon->SetBranchAddress("MVA_Continuum", &FBDT_var);

    for (unsigned int i = 0; i < tree_upsilon->GetEntries(); i++) {
        tree_upsilon->GetEntry(i);
        FBDTc_hist->Fill(FBDT_var, weight_var);

        NormFactor = NormFactor + BDTcToWeight(FBDT_var) * weight_var;

        nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + weight_var;
        nevt->NevtwithCorrection = nevt->NevtwithCorrection + BDTcToWeight(FBDT_var) * weight_var;
        nevt->BDTcSum = nevt->BDTcSum + FBDT_var * weight_var;

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

        NormFactor = NormFactor + BDTcToWeight(FBDT_var) * weight_var;

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

        FBDTc_hist->Fill(FBDT_var, BDTcToWeight(FBDT_var) * weight_var * NormFactor);

    }
    input_file->Close();

    return NormFactor;
}

void FastBDT_output_KS_test_BDTc_off()
{
    Nevt nevt_train = { 0.0, 0.0, 0.0 };
    Nevt nevt_test = { 0.0, 0.0, 0.0 };

    Nevt data_all = { 0.0, 0.0, 0.0 };

    // dirnames
    const char* side_data = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD_LS_data_off/temp_v009/SIGNAL";

    const char* side_MC_UUBAR_validation = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD_LS_MC_off/temp_v009/UUBAR";
    const char* side_MC_DDBAR_validation = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD_LS_MC_off/temp_v009/DDBAR";
    const char* side_MC_SSBAR_validation = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD_LS_MC_off/temp_v009/SSBAR";
    const char* side_MC_CHARM_validation = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD_LS_MC_off/temp_v009/CHARM";



    // train sample
    {
        std::vector<string> names;
        load_files(side_data, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 == 0) continue;
            LetsFill((side_data + std::string("/") + names.at(i)).c_str(), FBDTc_data_train, &data_all, 1.0 / 0.83152260744);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_UUBAR_validation, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 == 0) continue;
            LetsFill((side_MC_UUBAR_validation + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train, &nevt_train, 0.25 / 0.87399429016);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_DDBAR_validation, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 == 0) continue;
            LetsFill((side_MC_DDBAR_validation + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train, &nevt_train, 0.25 / 0.86799468791);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_SSBAR_validation, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 == 0) continue;
            LetsFill((side_MC_SSBAR_validation + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train, &nevt_train, 0.25 / 0.83277800916);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_CHARM_validation, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 == 0) continue;
            LetsFill((side_MC_CHARM_validation + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train, &nevt_train, 0.25 / 0.85138558986);
        }
    }

    // test sample
    {
        std::vector<string> names;
        load_files(side_data, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 != 0) continue;
            LetsFill((side_data + std::string("/") + names.at(i)).c_str(), FBDTc_data_test, &data_all, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_UUBAR_validation, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 != 0) continue;
            LetsFill((side_MC_UUBAR_validation + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test, &nevt_test, 0.25 / 0.12600570983);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_DDBAR_validation, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 != 0) continue;
            LetsFill((side_MC_DDBAR_validation + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test, &nevt_test, 0.25 / 0.13200531208);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_SSBAR_validation, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 != 0) continue;
            LetsFill((side_MC_SSBAR_validation + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test, &nevt_test, 0.25 / 0.16722199083);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_CHARM_validation, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 != 0) continue;
            LetsFill((side_MC_CHARM_validation + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test, &nevt_test, 0.25 / 0.14861441013);
        }
    }


    double NormFactorForTrain = nevt_train.NevtwithoutCorrection / nevt_train.NevtwithCorrection;
    double NormFactorForTest = nevt_test.NevtwithoutCorrection / nevt_test.NevtwithCorrection;


    // Corrected MC train
    {
        std::vector<string> names;
        load_files(side_MC_UUBAR_validation, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 == 0) continue;
            LetsFillwithCorrection((side_MC_UUBAR_validation + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train_correction, 0.25 / 0.87399429016, NormFactorForTrain);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_DDBAR_validation, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 == 0) continue;
            LetsFillwithCorrection((side_MC_DDBAR_validation + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train_correction, 0.25 / 0.86799468791, NormFactorForTrain);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_SSBAR_validation, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 == 0) continue;
            LetsFillwithCorrection((side_MC_SSBAR_validation + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train_correction, 0.25 / 0.83277800916, NormFactorForTrain);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_CHARM_validation, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 == 0) continue;
            LetsFillwithCorrection((side_MC_CHARM_validation + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train_correction, 0.25 / 0.85138558986, NormFactorForTrain);
        }
    }


    // Corrected MC test
    {
        std::vector<string> names;
        load_files(side_MC_UUBAR_validation, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 != 0) continue;
            LetsFillwithCorrection((side_MC_UUBAR_validation + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test_correction, 0.25 / 0.12600570983, NormFactorForTest);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_DDBAR_validation, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 != 0) continue;
            LetsFillwithCorrection((side_MC_DDBAR_validation + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test_correction, 0.25 / 0.13200531208, NormFactorForTest);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_SSBAR_validation, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 != 0) continue;
            LetsFillwithCorrection((side_MC_SSBAR_validation + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test_correction, 0.25 / 0.16722199083, NormFactorForTest);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_CHARM_validation, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 != 0) continue;
            LetsFillwithCorrection((side_MC_CHARM_validation + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test_correction, 0.25 / 0.14861441013, NormFactorForTest);
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

    // Use the same vertical scale for both panels.
    const double plotMax = 1.15 * std::max({
        FBDTc_MC_train->GetMaximum(),
        FBDTc_MC_test->GetMaximum(),
        FBDTc_MC_train_correction->GetMaximum(),
        FBDTc_MC_test_correction->GetMaximum(),
        FBDTc_data_train->GetMaximum(),
        FBDTc_data_test->GetMaximum()
    });

    auto drawBDTc = [&](const char* canvasName,
                        TH1F* mcTrain, TH1F* mcTest,
                        const char* outputName) {
        auto* canvas = new TCanvas(canvasName, "", 600, 600);
        canvas->SetCanvasSize(600, 600);
        canvas->SetLeftMargin(0.16);
        canvas->SetRightMargin(0.04);
        canvas->SetBottomMargin(0.14);
        canvas->SetTopMargin(0.05);
        canvas->cd();

        mcTrain->SetMinimum(0.0);
        mcTrain->SetMaximum(plotMax);

        for (TAxis* axis :
             {mcTrain->GetXaxis(), mcTrain->GetYaxis()}) {
            axis->SetLabelFont(43);
            axis->SetLabelSize(26);
            axis->SetTitleFont(43);
            axis->SetTitleSize(28);
            axis->SetLabelOffset(0.008);
        }

        mcTrain->GetXaxis()->SetTitle("BDT_{c}");
        mcTrain->GetXaxis()->SetTitleOffset(1.1);
        mcTrain->GetXaxis()->SetNdivisions(505);
        mcTrain->GetYaxis()->SetTitle("Arbitrary units");
        mcTrain->GetYaxis()->SetTitleOffset(1.15);

        mcTrain->Draw("HIST");
        FBDTc_data_train->Draw("HIST SAME");
        mcTest->Draw("E1 SAME");
        FBDTc_data_test->Draw("E1 SAME");

        auto* legend = new TLegend(0.62, 0.67, 0.95, 0.94);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextFont(43);
        legend->SetTextSize(26);
        legend->SetMargin(0.20);
        legend->AddEntry(mcTrain,          "MC (train)",   "f");
        legend->AddEntry(FBDTc_data_train, "Data (train)", "f");
        legend->AddEntry(mcTest,           "MC (test)",    "lep");
        legend->AddEntry(FBDTc_data_test,  "Data (test)",  "lep");
        legend->Draw();

        canvas->RedrawAxis();
        canvas->Modified();
        canvas->Update();
        canvas->Print((std::string(outputName) + ".pdf").c_str(),
                      "pdf Portrait");
        canvas->SaveAs((std::string(outputName) + ".png").c_str());
    };

    drawBDTc("c_BDTc_before",
             FBDTc_MC_train, FBDTc_MC_test,
             "BDTc_BeforeCorrection");
    drawBDTc("c_BDTc_after",
             FBDTc_MC_train_correction, FBDTc_MC_test_correction,
             "BDTc_AfterCorrection");

    printf("average BDTc of data: %lf\n", data_all.BDTcSum / data_all.NevtwithoutCorrection);
}
