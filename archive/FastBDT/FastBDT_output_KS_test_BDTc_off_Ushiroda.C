#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
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

// Export four comparisons separately and together as a 2x2 figure.
void DrawBDTcComparisons()
{
    gStyle->SetOptStat(0);

    // Share the axis ranges and leave room above the error bars for the legend.
    double plotMax = 0.0;
    for (TH1F* histogram : {FBDTc_MC_train, FBDTc_MC_test,
                            FBDTc_data_train, FBDTc_data_test,
                            FBDTc_MC_test_correction}) {
        for (int bin = 1; bin <= histogram->GetNbinsX(); ++bin) {
            plotMax = std::max(plotMax, histogram->GetBinContent(bin)
                                       + histogram->GetBinError(bin));
        }
    }
    plotMax *= 1.30;

    auto drawBDTc = [&](const char* canvasName,
                        TH1F* referenceInput, TH1F* pointsInput,
                        const char* referenceLabel, const char* pointsLabel,
                        bool showReferenceError, const char* outputName) -> TCanvas* {
        auto* canvas = new TCanvas(canvasName, "", 600, 600);
        canvas->SetCanvasSize(600, 600);
        canvas->SetLeftMargin(0.16);
        canvas->SetRightMargin(0.04);
        canvas->SetBottomMargin(0.14);
        canvas->SetTopMargin(0.05);
        canvas->cd();

        // MC test is a point series in (a), but a histogram in (c).
        // Clone each series so styling one panel cannot alter another panel.
        auto* reference = static_cast<TH1F*>(referenceInput->Clone(
            (std::string(canvasName) + "_reference").c_str()));
        auto* points = static_cast<TH1F*>(pointsInput->Clone(
            (std::string(canvasName) + "_points").c_str()));
        reference->SetDirectory(nullptr);
        points->SetDirectory(nullptr);

        // Use the same visual roles in all panels: blue histogram, black points.
        reference->SetLineColor(kBlue + 1);
        reference->SetLineWidth(2);
        reference->SetFillStyle(0);
        reference->SetMinimum(0.0);
        reference->SetMaximum(plotMax);
        points->SetMarkerStyle(kFullCircle);
        points->SetMarkerColor(kBlack);
        points->SetLineColor(kBlack);
        points->SetLineWidth(1);
        points->SetFillStyle(0);

        for (TAxis* axis : {reference->GetXaxis(), reference->GetYaxis()}) {
            axis->SetLabelFont(43);
            axis->SetLabelSize(26);
            axis->SetTitleFont(43);
            axis->SetTitleSize(28);
            axis->SetLabelOffset(0.008);
        }
        reference->GetXaxis()->SetTitle("BDT_{c}");
        reference->GetXaxis()->SetTitleOffset(1.1);
        reference->GetXaxis()->SetNdivisions(505);
        reference->GetYaxis()->SetTitle("Arbitrary units");
        reference->GetYaxis()->SetTitleOffset(1.15);

        reference->Draw("HIST");
        TH1F* referenceBand = nullptr;
        if (showReferenceError) {
            // Show MC errors as a band, without a second set of points.
            referenceBand = static_cast<TH1F*>(reference->Clone(
                (std::string(canvasName) + "_uncertainty").c_str()));
            referenceBand->SetDirectory(nullptr);
            referenceBand->SetFillColorAlpha(kBlue + 1, 0.20);
            referenceBand->SetFillStyle(1001);
            referenceBand->Draw("E2 SAME");
            reference->Draw("HIST SAME");
        }
        points->Draw("E1 SAME");

        auto* legend = new TLegend(0.44, 0.81, 0.95, 0.94);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextFont(43);
        legend->SetTextSize(24);
        legend->SetMargin(0.16);
        legend->AddEntry(referenceBand ? referenceBand : reference,
                         referenceLabel, showReferenceError ? "lf" : "l");
        legend->AddEntry(points, pointsLabel, "lep");
        legend->Draw();

        canvas->RedrawAxis();
        canvas->Modified();
        canvas->Update();
        canvas->Print((std::string(outputName) + ".pdf").c_str(),
                      "pdf Portrait");
        canvas->SaveAs((std::string(outputName) + ".png").c_str());
        return canvas;
    };

    TCanvas* panels[] = {
        drawBDTc("c_BDTc_mc_train_test",
                 FBDTc_MC_train, FBDTc_MC_test,
                 "MC (train)", "MC (test)", false,
                 "BDTc_MCTrainTest"),
        drawBDTc("c_BDTc_data_train_test",
                 FBDTc_data_train, FBDTc_data_test,
                 "Data (train)", "Data (test)", false,
                 "BDTc_DataTrainTest"),
        // The same test data are used before and after the MC correction.
        drawBDTc("c_BDTc_test_data_uncorrected_mc",
                 FBDTc_MC_test, FBDTc_data_test,
                 "Uncorrected MC (test)", "Data (test)", true,
                 "BDTc_TestDataUncorrectedMC"),
        drawBDTc("c_BDTc_test_data_corrected_mc",
                 FBDTc_MC_test_correction, FBDTc_data_test,
                 "Corrected MC (test)", "Data (test)", true,
                 "BDTc_TestDataCorrectedMC")
    };

    auto* comparison = new TCanvas("c_BDTc_comparison_2x2", "", 1200, 1200);
    comparison->SetCanvasSize(1200, 1200);
    comparison->Divide(2, 2, 0.005, 0.005);
    const char* panelLabels[] = {"(a)", "(b)", "(c)", "(d)"};
    for (int panel = 0; panel < 4; ++panel) {
        comparison->cd(panel + 1);
        panels[panel]->DrawClonePad();
        comparison->cd(panel + 1);
        TLatex label;
        label.SetNDC();
        label.SetTextFont(43);
        label.SetTextSize(26);
        label.DrawLatex(0.19, 0.91, panelLabels[panel]);
    }
    comparison->Modified();
    comparison->Update();
    comparison->Print("BDTc_Comparison_2x2.pdf", "pdf Portrait");
    comparison->SaveAs("BDTc_Comparison_2x2.png");
}

void FastBDT_output_KS_test_BDTc_off_Ushiroda()
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


    DrawBDTcComparisons();

    printf("average BDTc of data: %lf\n", data_all.BDTcSum / data_all.NevtwithoutCorrection);
}
