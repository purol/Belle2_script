#include <algorithm>
#include <cmath>
#include <limits>
#include "TError.h"
#include "TRandom3.h"
#include <string>
#include <vector>
#include "TPaveStats.h"

namespace {

struct ToyBootstrapStats {
    double mean = 0.0;
    double stddev = 0.0;
    double median = 0.0;
    double meanError = 0.0;
    double stddevError = 0.0;
    double medianError = 0.0;
};

struct ToyRunningMoments {
    std::size_t count = 0;
    double mean = 0.0;
    double m2 = 0.0;

    void add(double value) {
        ++count;
        const double delta = value - mean;
        mean += delta / static_cast<double>(count);
        m2 += delta * (value - mean);
    }

    double stddev() const {
        // Match the ROOT histogram Std Dev convention (denominator N).
        return std::sqrt(std::max(0.0, m2 / static_cast<double>(count)));
    }

    double standardError() const {
        // Standard deviation across bootstrap replicates, not SD / sqrt(B).
        return std::sqrt(std::max(0.0, m2 / static_cast<double>(count - 1)));
    }
};

double toyMedian(std::vector<double>& values) {
    std::sort(values.begin(), values.end());
    const std::size_t middle = values.size() / 2;
    return values.size() % 2 == 0
        ? 0.5 * values[middle - 1] + 0.5 * values[middle]
        : values[middle];
}

bool bootstrapToyStats(std::vector<double> values, bool includeMedian,
                       unsigned int replicates, UInt_t seed,
                       ToyBootstrapStats& result) {
    if (values.size() < 2 || replicates < 2
        || values.size() > std::numeric_limits<UInt_t>::max()) {
        Error("bootstrapToyStats", "Need at least two values and two replicates, with N fitting in UInt_t.");
        return false;
    }
    for (double value : values) {
        if (!std::isfinite(value)) {
            Error("bootstrapToyStats", "A toy value is not finite.");
            return false;
        }
    }

    // Use individual toy values, including histogram under/overflow entries.
    // Sorting also makes the fixed-seed result independent of input-file order.
    std::sort(values.begin(), values.end());
    ToyRunningMoments original;
    for (double value : values) original.add(value);
    result = ToyBootstrapStats();
    result.mean = original.mean;
    result.stddev = original.stddev();
    if (includeMedian) result.median = toyMedian(values);

    TRandom3 random(seed);
    std::vector<double> sample(includeMedian ? values.size() : 0);
    ToyRunningMoments means, stddevs, medians;
    const UInt_t count = static_cast<UInt_t>(values.size());
    for (unsigned int b = 0; b < replicates; ++b) {
        ToyRunningMoments resampled;
        for (UInt_t j = 0; j < count; ++j) {
            const double value = values[random.Integer(count)];
            resampled.add(value);
            if (includeMedian) sample[j] = value;
        }
        means.add(resampled.mean);
        stddevs.add(resampled.stddev());
        if (includeMedian) medians.add(toyMedian(sample));
    }
    result.meanError = means.standardError();
    result.stddevError = stddevs.standardError();
    if (includeMedian) result.medianError = medians.standardError();
    if (!std::isfinite(result.mean) || !std::isfinite(result.stddev)
        || !std::isfinite(result.meanError) || !std::isfinite(result.stddevError)
        || (includeMedian && (!std::isfinite(result.median) || !std::isfinite(result.medianError)))) {
        Error("bootstrapToyStats", "The summary or bootstrap standard error is not finite.");
        return false;
    }
    return true;
}

} // namespace


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

void ReadToyRootFile_MXs3(){

    const unsigned int bootstrapReplicates = 5000;
    const UInt_t bootstrapSeed = 20261120;

    RooRealVar mu_pull_roorealvar("mu_pull", ";pull of #mu;Toys", -4.0, 4.0);
    RooDataSet mu_pull_RooDataSet("mu_pull_RooDataSet", "mu_pull_RooDataSet", RooArgSet(mu_pull_roorealvar));

    TH1F* ToyMCmu = new TH1F("ToyMCmu", ";#mu;Toys", 40, -40, 40);
    TH1F* ToyMCmuerror = new TH1F("ToyMCmuerror", ";error of #mu;Toys", 50, 6, 16);
    TH1F* ToyMCmuHIerror = new TH1F("ToyMCmuHIerror", ";error of #mu;Toys", 50, 6, 16);
    TH1F* ToyMCmuLOerror = new TH1F("ToyMCmuLOerror", ";error of #mu;Toys", 50, -16, -6);
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

    std::vector<double> fittedMuValues;
    std::vector<double> positiveErrorValues;
    std::vector<double> negativeErrorValues;
    std::vector<string> names;
    const char* dirname = "./";
    load_files(dirname, &names, "TOY_result");
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        if (input_file->IsZombie()) {
            Error("ReadToyRootFile", "Cannot open %s", names.at(i).c_str());
            delete input_file;
            return;
        }
        TTree* temp_tree = (TTree*)input_file->Get("TOY_result");
        if (!temp_tree) {
            Error("ReadToyRootFile", "Missing TOY_result tree in %s", names.at(i).c_str());
            delete input_file;
            return;
        }

        double temp_mu_true = -1;
        double temp_mu_fitting = -1;
        double temp_mu_error = -1;
        double temp_mu_HIerror = -1;
        double temp_mu_LOerror = -1;

        int temp_covQual = -1;
        int temp_status = -1;
        double temp_edm = -1;

        temp_tree->SetBranchAddress("mu_MXs3_true", &temp_mu_true);
        if (temp_tree->SetBranchAddress("mu_MXs3_value", &temp_mu_fitting) < 0) {
            Error("ReadToyRootFile", "Missing mu_MXs3_value branch in %s", names.at(i).c_str());
            delete input_file;
            return;
        }
        temp_tree->SetBranchAddress("mu_MXs3_error", &temp_mu_error);
        if (temp_tree->SetBranchAddress("mu_MXs3_HIerror", &temp_mu_HIerror) < 0) {
            Error("ReadToyRootFile", "Missing mu_MXs3_HIerror branch in %s", names.at(i).c_str());
            delete input_file;
            return;
        }
        if (temp_tree->SetBranchAddress("mu_MXs3_LOerror", &temp_mu_LOerror) < 0) {
            Error("ReadToyRootFile", "Missing mu_MXs3_LOerror branch in %s", names.at(i).c_str());
            delete input_file;
            return;
        }

        temp_tree->SetBranchAddress("covQual", &temp_covQual);
        temp_tree->SetBranchAddress("status", &temp_status);
        temp_tree->SetBranchAddress("edm", &temp_edm);

        for (unsigned int j = 0; j < temp_tree->GetEntries(); j++) { // Fill
            if (temp_tree->GetEntry(j) <= 0 || !std::isfinite(temp_mu_fitting)
                || !std::isfinite(temp_mu_HIerror) || !std::isfinite(temp_mu_LOerror)) {
                Error("ReadToyRootFile", "Invalid toy value in %s at entry %u", names.at(i).c_str(), j);
                delete input_file;
                return;
            }
            ToyMCmu->Fill(temp_mu_fitting);
            fittedMuValues.push_back(temp_mu_fitting);
            ToyMCmuerror->Fill(temp_mu_error);
            ToyMCmuHIerror->Fill(temp_mu_HIerror);
            positiveErrorValues.push_back(temp_mu_HIerror);
            ToyMCmuLOerror->Fill(temp_mu_LOerror);
            negativeErrorValues.push_back(temp_mu_LOerror);
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
        delete input_file;
    }

    ToyBootstrapStats muStats, positiveErrorStats, negativeErrorStats;
    if (!bootstrapToyStats(fittedMuValues, true, bootstrapReplicates, bootstrapSeed, muStats)
        || !bootstrapToyStats(positiveErrorValues, false, bootstrapReplicates, bootstrapSeed + 1, positiveErrorStats)
        || !bootstrapToyStats(negativeErrorValues, false, bootstrapReplicates, bootstrapSeed + 2, negativeErrorStats)) {
        return;
    }
    auto reportBootstrap = [&](const char* quantity, const ToyBootstrapStats& summary, bool includeMedian) {
        Printf("Bootstrap MXs3 %s: N=%zu mean=%.10g meanSE=%.10g stddev=%.10g stddevSE=%.10g",
               quantity, fittedMuValues.size(), summary.mean, summary.meanError, summary.stddev, summary.stddevError);
        if (includeMedian) Printf("Bootstrap MXs3 median=%.10g medianSE=%.10g", summary.median, summary.medianError);
    };
    reportBootstrap("value", muStats, true);
    reportBootstrap("HIerror", positiveErrorStats, false);
    reportBootstrap("LOerror", negativeErrorStats, false);

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
        yaxis->SetTitleOffset(1.25);
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
        const ToyBootstrapStats* summary = histogram == ToyMCmu ? &muStats
            : histogram == ToyMCmuHIerror ? &positiveErrorStats
            : histogram == ToyMCmuLOerror ? &negativeErrorStats : nullptr;
        const bool showMedian = histogram == ToyMCmu;
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
        histogram->SetMaximum(ymax > 0.0 ? ymax * (hasFit ? 2.0 : (showMedian ? 1.7 : 1.5)) : 1.0);
        histogram->Draw("PE1");
        canvas->Update();
        TPaveStats* stats = dynamic_cast<TPaveStats*>(histogram->FindObject("stats"));
        if (stats) {
            stats->SetX1NDC(summary ? 0.37 : 0.51);
            stats->SetX2NDC(0.94);
            stats->SetY1NDC(hasFit ? 0.60 : (showMedian ? 0.665 : 0.72));
            stats->SetY2NDC(0.94);
            stats->SetTextFont(42);
            stats->SetTextSize(0.040);
            stats->SetBorderSize(1);
            stats->SetFillColor(kWhite);
            stats->SetFillStyle(1001);
            if (summary) {
                // All central values and errors use the same unbinned toy sample.
                // Preserve the custom entries when the canvas is updated or saved.
                stats->SetName(TString::Format("%s_bootstrap_stats", histogram->GetName()).Data());
                stats->Clear();
                stats->AddText(histogram->GetName());
                stats->AddText(TString::Format("Entries = %.0f", histogram->GetEntries()).Data());
                stats->AddText(TString::Format("Mean = %.4g #pm %.2g", summary->mean, summary->meanError).Data());
                stats->AddText(TString::Format("Std Dev = %.4g #pm %.2g", summary->stddev, summary->stddevError).Data());
                if (showMedian) {
                    stats->AddText(TString::Format("Median = %.4g #pm %.2g", summary->median, summary->medianError).Data());
                }
                histogram->SetStats(kFALSE);
            }
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
    muframe->SetMaximum(pullFrameMaximum > 0.0 ? pullFrameMaximum * 1.6 : 1.0);

    TCanvas* c = makeToyCanvas("mu_gauss_fit");
    muframe->Draw();
    saveToyCanvas(c, "TOYMC_mu_MXs3_pull_unbinned");
    delete c;

    const int previousOptStat = gStyle->GetOptStat();
    const int previousOptFit = gStyle->GetOptFit();
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(11);

    drawToyHistogram(ToyMCmu, "TOYMC_mu_MXs3", false);
    drawToyHistogram(ToyMCmuerror, "TOYMC_mu_MXs3_error", false);
    drawToyHistogram(ToyMCmuHIerror, "TOYMC_mu_MXs3_HIerror", false);
    drawToyHistogram(ToyMCmuLOerror, "TOYMC_mu_MXs3_LOerror", false);

    drawToyHistogram(ToyMCmupull, "TOYMC_mu_MXs3_pull", true);
    gStyle->SetOptStat(previousOptStat);
    gStyle->SetOptFit(previousOptFit);

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
