#include <algorithm>
#include <string>
#include "TPaveStats.h"

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

void ReadToyRootFile_MXs2(){

    RooRealVar mu_pull_roorealvar("mu_pull", ";pull of #mu;Toys", -4.0, 4.0);
    RooDataSet mu_pull_RooDataSet("mu_pull_RooDataSet", "mu_pull_RooDataSet", RooArgSet(mu_pull_roorealvar));

    TH1F* ToyMCmu = new TH1F("ToyMCmu", ";#mu;Toys", 40, -25, 25);
    TH1F* ToyMCmuerror = new TH1F("ToyMCmuerror", ";error of #mu;Toys", 50, 2, 8);
    TH1F* ToyMCmuHIerror = new TH1F("ToyMCmuHIerror", ";error of #mu;Toys", 50, 2, 8);
    TH1F* ToyMCmuLOerror = new TH1F("ToyMCmuLOerror", ";error of #mu;Toys", 50, -8, -2);
    TH1F* ToyMCmupull = new TH1F("ToyMCmupull", ";pull of #mu;Toys", 41, -4, 4);

    ToyMCmu->SetMarkerStyle(kFullCircle);
    ToyMCmu->SetLineColor(kBlack);
    ToyMCmu->SetMarkerColor(kBlack);

    ToyMCmuerror->SetMarkerStyle(kFullCircle);
    ToyMCmuerror->SetLineColor(kBlack);
    ToyMCmuerror->SetMarkerColor(kBlack);

    ToyMCmuHIerror->SetMarkerStyle(kFullCircle);
    ToyMCmuHIerror->SetLineColor(kBlack);
    ToyMCmuHIerror->SetMarkerColor(kBlack);

    ToyMCmuLOerror->SetMarkerStyle(kFullCircle);
    ToyMCmuLOerror->SetLineColor(kBlack);
    ToyMCmuLOerror->SetMarkerColor(kBlack);

    ToyMCmupull->SetMarkerStyle(kFullCircle);
    ToyMCmupull->SetLineColor(kBlack);
    ToyMCmupull->SetMarkerColor(kBlack);

    std::vector<string> names;
    const char* dirname = "./";
    load_files(dirname, &names, "TOY_result");
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        TTree* temp_tree = (TTree*)input_file->Get("TOY_result");

        double temp_mu_true = -1;
        double temp_mu_fitting = -1;
        double temp_mu_error = -1;
        double temp_mu_HIerror = -1;
        double temp_mu_LOerror = -1;

        int temp_covQual = -1;
        int temp_status = -1;
        double temp_edm = -1;

        temp_tree->SetBranchAddress("mu_MXs2_true", &temp_mu_true);
        temp_tree->SetBranchAddress("mu_MXs2_value", &temp_mu_fitting);
        temp_tree->SetBranchAddress("mu_MXs2_error", &temp_mu_error);
        temp_tree->SetBranchAddress("mu_MXs2_HIerror", &temp_mu_HIerror);
        temp_tree->SetBranchAddress("mu_MXs2_LOerror", &temp_mu_LOerror);

        temp_tree->SetBranchAddress("covQual", &temp_covQual);
        temp_tree->SetBranchAddress("status", &temp_status);
        temp_tree->SetBranchAddress("edm", &temp_edm);

        for (unsigned int j = 0; j < temp_tree->GetEntries(); j++) { // Fill
            temp_tree->GetEntry(j);
            ToyMCmu->Fill(temp_mu_fitting);
            ToyMCmuerror->Fill(temp_mu_error);
            ToyMCmuHIerror->Fill(temp_mu_HIerror);
            ToyMCmuLOerror->Fill(temp_mu_LOerror);
            if (temp_mu_true >= temp_mu_fitting) {
                ToyMCmupull->Fill((temp_mu_true - temp_mu_fitting) / temp_mu_HIerror);
                mu_pull_roorealvar = (temp_mu_true - temp_mu_fitting) / temp_mu_HIerror;
                mu_pull_RooDataSet.add(RooArgSet(mu_pull_roorealvar));
            }
            else {
                ToyMCmupull->Fill((temp_mu_fitting - temp_mu_true) / temp_mu_LOerror);
                mu_pull_roorealvar = (temp_mu_fitting - temp_mu_true) / temp_mu_LOerror;
                mu_pull_RooDataSet.add(RooArgSet(mu_pull_roorealvar));
            }
        }
        input_file->Close();
    }

    RooRealVar gausmean("mean", "", 0.0, -1.0, 1.0);
    RooRealVar gauswidth("sigma", "", 1.0, 0.9, 1.1);
    RooGaussian gauss("gauss", "gauss", mu_pull_roorealvar, gausmean, gauswidth);
    RooRealVar nentry("nentry", "number of entries", 10000, 9500, 10500);
    RooExtendPdf egauss("egauss", "extended gauss", gauss, nentry);
    egauss.fitTo(mu_pull_RooDataSet, RooFit::Save());

    // Shared display settings for figures placed at 0.49 text width in the thesis.
    auto styleToyAxes = [](TAxis* xaxis, TAxis* yaxis) {
        TAxis* axes[] = {xaxis, yaxis};
        for (TAxis* axis : axes) {
            axis->SetTitleFont(42);
            axis->SetLabelFont(42);
            axis->SetTitleSize(0.050);
            axis->SetLabelSize(0.045);
            axis->SetLabelOffset(0.008);
            axis->SetNdivisions(505);
        }
        xaxis->SetTitleOffset(1.10);
        yaxis->SetTitleOffset(1.65);
    };

    auto makeToyCanvas = [](const char* name) {
        TCanvas* canvas = new TCanvas(name, "", 800, 800);
        canvas->SetLeftMargin(0.20);
        canvas->SetRightMargin(0.05);
        canvas->SetBottomMargin(0.16);
        canvas->SetTopMargin(0.06);
        canvas->SetTicks(1, 1);
        return canvas;
    };

    auto saveToyCanvas = [](TCanvas* canvas, const char* basename) {
        canvas->Modified();
        canvas->Update();
        canvas->Print((std::string(basename) + ".pdf").c_str(), "pdf Portrait");
        canvas->SaveAs((std::string(basename) + ".png").c_str());
    };

    auto drawToyHistogram = [&](TH1* histogram, const char* basename, bool hasFit) {
        TCanvas* canvas = makeToyCanvas(basename);
        if (hasFit) histogram->Fit("gaus");
        styleToyAxes(histogram->GetXaxis(), histogram->GetYaxis());
        histogram->SetLineWidth(2);
        histogram->SetMarkerSize(0.8);
        histogram->SetStats(kTRUE);
        histogram->SetMinimum(0.0);
        double ymax = 0.0;
        for (int bin = 1; bin <= histogram->GetNbinsX(); ++bin) {
            ymax = std::max(ymax, histogram->GetBinContent(bin) + histogram->GetBinError(bin));
        }
        // Keep the enlarged statistics/fit box above the data and error bars.
        histogram->SetMaximum(ymax > 0.0 ? ymax * (hasFit ? 2.0 : 1.5) : 1.0);
        histogram->Draw("PE1");
        canvas->Update();
        TPaveStats* stats = dynamic_cast<TPaveStats*>(histogram->FindObject("stats"));
        if (stats) {
            stats->SetX1NDC(hasFit ? 0.40 : 0.51);
            stats->SetX2NDC(0.94);
            stats->SetY1NDC(hasFit ? 0.60 : 0.72);
            stats->SetY2NDC(0.94);
            stats->SetTextFont(42);
            stats->SetTextSize(0.040);
            stats->SetBorderSize(1);
            stats->SetFillColor(kWhite);
            stats->SetFillStyle(1001);
        }
        saveToyCanvas(canvas, basename);
        delete canvas;
    };

    RooPlot* muframe = mu_pull_roorealvar.frame(RooFit::Bins(41), RooFit::Title(" "));
    mu_pull_RooDataSet.plotOn(muframe);
    egauss.plotOn(muframe, RooFit::LineColor(kRed));
    egauss.paramOn(muframe, RooFit::Layout(0.43, 0.94, 0.94));
    muframe->getAttText()->SetTextFont(42);
    muframe->getAttText()->SetTextSize(0.040);
    styleToyAxes(muframe->GetXaxis(), muframe->GetYaxis());
    muframe->GetXaxis()->SetTitle("pull of #mu");
    muframe->GetYaxis()->SetTitle("Toys");
    const double pullFrameMaximum = muframe->GetMaximum();
    muframe->SetMinimum(0.0);
    muframe->SetMaximum(pullFrameMaximum > 0.0 ? pullFrameMaximum * 1.3 : 1.0);

    TCanvas* c = makeToyCanvas("mu_gauss_fit");
    muframe->Draw();
    saveToyCanvas(c, "TOYMC_mu_MXs2_pull_unbinned");
    delete c;

    const int previousOptStat = gStyle->GetOptStat();
    const int previousOptFit = gStyle->GetOptFit();
    const std::string previousFitFormat = gStyle->GetFitFormat();
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(11);
    gStyle->SetFitFormat(".3g");

    drawToyHistogram(ToyMCmu, "TOYMC_mu_MXs2", false);
    drawToyHistogram(ToyMCmuerror, "TOYMC_mu_MXs2_error", false);
    drawToyHistogram(ToyMCmuHIerror, "TOYMC_mu_MXs2_HIerror", false);
    drawToyHistogram(ToyMCmuLOerror, "TOYMC_mu_MXs2_LOerror", false);

    drawToyHistogram(ToyMCmupull, "TOYMC_mu_MXs2_pull", true);
    gStyle->SetOptStat(previousOptStat);
    gStyle->SetOptFit(previousOptFit);
    gStyle->SetFitFormat(previousFitFormat.c_str());

    TF1* fit_pull_gauss = ToyMCmupull->GetFunction("gaus");
    const double constant_gauss = fit_pull_gauss->GetParameter(0);
    const double mean_gauss = fit_pull_gauss->GetParameter(1);
    const double sigma_gauss = fit_pull_gauss->GetParameter(2);

    delete ToyMCmu;
    delete ToyMCmuerror;
    delete ToyMCmuHIerror;
    delete ToyMCmuLOerror;
    delete ToyMCmupull;

}
