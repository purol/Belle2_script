#include <algorithm>
#include <cmath>
#include "TError.h"
#include "TRandom3.h"
#include <string>
#include <vector>
#include "TPaveStats.h"

#define LT_number 11

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

void ReadLinearityRootFile_median(){

    // Nonparametric bootstrap of the fitted values, separately at each injected mu.
    const unsigned int bootstrapReplicates = 5000;
    const unsigned int bootstrapSeed = 20261120;
    TRandom3 bootstrapRandom(bootstrapSeed);
    if (bootstrapReplicates < 2) {
        Error("ReadLinearityRootFile_median", "At least two bootstrap replicates are required.");
        return;
    }

    std::vector<double> mus;
    std::vector<double> mus_err;

    double Inputmu[LT_number] = { 0 };
    double Inputmuerror[LT_number] = { 0 };
    double outputmu[LT_number] = { 0 };
    double outputmuerror[LT_number] = { 0 };

    for (int i = 0; i < LT_number; i++) {
        double injected_mu = i * 1.0;
        Inputmu[i] = injected_mu;

        std::vector<double> out_mus;
        unsigned int Nentry = 0;

        std::vector<string> names;
        const char* dirname = "./";
        load_files(dirname, &names, ("LT_result_" + std::to_string(injected_mu)).c_str() );
        for (unsigned int k = 0; k < names.size(); k++) {

            TFile* input_file = new TFile((dirname + std::string("/") + names.at(k)).c_str(), "read");
            if (input_file->IsZombie()) {
                Error("ReadLinearityRootFile_median", "Cannot open toy file: %s", names.at(k).c_str());
                delete input_file;
                return;
            }
            TTree* temp_tree = (TTree*)input_file->Get("LT_result");
            if (!temp_tree) {
                Error("ReadLinearityRootFile_median", "Missing LT_result tree in %s", names.at(k).c_str());
                delete input_file;
                return;
            }

            double out_mu;
            Nentry = temp_tree->GetEntries();

            if (temp_tree->SetBranchAddress("mu_MXs3_value", &out_mu) < 0) {
                Error("ReadLinearityRootFile_median", "Cannot read mu_MXs3_value in %s", names.at(k).c_str());
                delete input_file;
                return;
            }

            for (unsigned int j = 0; j < Nentry; j++) { // Fill
                if (temp_tree->GetEntry(j) <= 0 || !std::isfinite(out_mu)) {
                    Error("ReadLinearityRootFile_median", "Invalid fitted value in %s at entry %u", names.at(k).c_str(), j);
                    delete input_file;
                    return;
                }
                out_mus.push_back(out_mu);
            }
            input_file->Close();
            delete input_file;
        }
        if (out_mus.size() < 2) {
            Error("ReadLinearityRootFile_median", "At least two fitted values are required at injected mu = %g", injected_mu);
            return;
        }

        RooRealVar  mu_roorealvar("mu_roorealvar", "mu", injected_mu - 30.0, injected_mu + 30.0);
        RooDataSet mu_RooDataSet("mu_RooDataSet", "mu_RooDataSet", RooArgSet(mu_roorealvar));
        for (int j = 0; j < out_mus.size(); j++) {
            mu_roorealvar = out_mus.at(j);
            mu_RooDataSet.add(RooArgSet(mu_roorealvar));
        }

        // get median
        std::sort(out_mus.begin(), out_mus.end());
        size_t size = out_mus.size();

        outputmu[i] = out_mus[size / 2];

        // Resample N fitted values with replacement. Keep the existing median
        // convention: the upper middle value when N is even.
        std::vector<double> bootstrapSample(size);
        double bootstrapMean = 0.0;
        double bootstrapM2 = 0.0;
        for (unsigned int b = 0; b < bootstrapReplicates; ++b) {
            for (std::size_t j = 0; j < size; ++j) {
                bootstrapSample[j] = out_mus[bootstrapRandom.Integer(static_cast<UInt_t>(size))];
            }
            auto middle = bootstrapSample.begin() + size / 2;
            std::sort(bootstrapSample.begin(), bootstrapSample.end());
            const double bootstrapMedian = *middle;
            // Welford's algorithm for the variance across bootstrap medians.
            const double delta = bootstrapMedian - bootstrapMean;
            bootstrapMean += delta / (b + 1);
            bootstrapM2 += delta * (bootstrapMedian - bootstrapMean);
        }
        // The spread of the bootstrap medians estimates the median SE.
        outputmuerror[i] = std::sqrt(bootstrapM2 / (bootstrapReplicates - 1));
        if (!std::isfinite(outputmuerror[i]) || outputmuerror[i] <= 0.0) {
            Error("ReadLinearityRootFile_median", "Invalid bootstrap standard error at injected mu = %g", injected_mu);
            return;
        }
        Printf("Injected mu = %g, N = %zu, median = %.8g, bootstrap SE = %.8g",
               injected_mu, size, outputmu[i], outputmuerror[i]);

        RooPlot* muframe = mu_roorealvar.frame(RooFit::Bins(400), RooFit::Title(" "));
        mu_RooDataSet.plotOn(muframe);
        TCanvas* c = new TCanvas("mu_gauss_fit", "mu_gauss_fit", 600, 600);
        gPad->SetLeftMargin(0.15); muframe->GetYaxis()->SetTitleOffset(1.4); muframe->Draw(); c->SaveAs(("mu_" + to_string(injected_mu) + "_MXs3_distribution.png").c_str());

    }

    // Match the toy plots when placed at 0.49 text width in the thesis.
    const int previousOptFit = gStyle->GetOptFit();
    const std::string previousFitFormat = gStyle->GetFitFormat();
    gStyle->SetOptFit(11);
    gStyle->SetFitFormat(".3g");

    TCanvas* c = new TCanvas("Linearity test canvas", "", 800, 800);
    c->SetLeftMargin(0.20);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.16);
    c->SetTopMargin(0.06);
    c->SetTicks(1, 1);

    TGraphErrors* gr = new TGraphErrors(LT_number, Inputmu, outputmu, Inputmuerror, outputmuerror);
    gr->SetMarkerStyle(21);
    gr->SetMarkerSize(1.1);
    gr->SetLineWidth(2);
    gr->SetTitle(";Injected #mu;Median fitted #mu");
    gr->Fit("pol1");
    TF1* linearFit = gr->GetFunction("pol1");
    if (linearFit) {
        linearFit->SetParName(0, "Intercept");
        linearFit->SetParName(1, "Slope");
        linearFit->SetLineColor(kRed);
        linearFit->SetLineWidth(2);
    }

    TAxis* axes[] = {gr->GetXaxis(), gr->GetYaxis()};
    for (TAxis* axis : axes) {
        axis->SetTitleFont(42);
        axis->SetLabelFont(42);
        axis->SetTitleSize(0.050);
        axis->SetLabelSize(0.045);
        axis->SetLabelOffset(0.008);
        axis->SetNdivisions(505);
    }
    gr->GetXaxis()->SetTitleOffset(1.10);
    gr->GetYaxis()->SetTitleOffset(1.65);
    gr->Draw("AP");
    c->Update();

    // The upper-left corner is clear of the fitted line and data points.
    TPaveStats* stats = dynamic_cast<TPaveStats*>(gr->GetListOfFunctions()->FindObject("stats"));
    if (stats) {
        stats->SetX1NDC(0.23);
        stats->SetX2NDC(0.74);
        stats->SetY1NDC(0.78);
        stats->SetY2NDC(0.92);
        stats->SetTextFont(42);
        stats->SetTextSize(0.040);
        stats->SetBorderSize(1);
        stats->SetFillColor(kWhite);
        stats->SetFillStyle(1001);
    }
    c->Modified();
    c->Update();
    c->Print("Linearity_MXs3_test_median.pdf", "pdf Portrait");
    c->SaveAs("Linearity_MXs3_test_median.png");

    /*
    gr->SetMinimum(-1);
    gr->SetMaximum(12.0);
    gr->Draw("AP");
    c->Modified(); c->Update(); c->SaveAs("Linearity_test_magnify.png");
    */

    delete gr;
    delete c;

    gStyle->SetOptFit(previousOptFit);
    gStyle->SetFitFormat(previousFitFormat.c_str());

}
