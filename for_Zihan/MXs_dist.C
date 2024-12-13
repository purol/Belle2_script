#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"

bool hasEnding(std::string const& fullString, std::string const& ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
    }
    else {
        return false;
    }
}


void load_files(const char* dirname, std::vector<std::string>* names) {
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

void load_files(const char* dirname, std::vector<std::string>* names, const char* included_string) {
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

TH1D* MXs_11 = new TH1D("MXs cut = 1.1 GeV", ";M_{Xs} [GeV];arbitrary unit", 30, 1.1, 4.0);
TH1D* MXs_10 = new TH1D("MXs cut = 1.0 GeV", ";M_{Xs} [GeV];arbitrary unit", 30, 1.1, 4.0);

void LetsFill(const char* dirname, const char* included_string, TH1D* MXs_hist, double weight_var = 1.0) {
    double Mx = 0;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");

        tree_upsilon->SetBranchAddress("mcDaughter__bo0__cm__spmcDaughter__bo0__cm__spM__bc__bc", &Mx);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);

            // cut
            if (Mx < 1.1) continue;

            MXs_hist->Fill(Mx, weight_var);
        }
        input_file->Close();

    }

}

void MXs_dist()
{

    // get data from root files
    const char* dirname = "./";

    LetsFill(dirname, "B02Xsll_11", MXs_11, 1.0);
    LetsFill(dirname, "B02Xsll_10", MXs_10, 1.0);

    double factor = 1.0;

    // normalization
    MXs_11->Scale(factor / MXs_11->Integral(), "width");
    MXs_10->Scale(factor / MXs_10->Integral(), "width");

    MXs_11->SetFillStyle(3004);
    MXs_11->SetLineColor(kBlue);
    MXs_11->SetFillColor(kBlue);

    MXs_10->SetFillStyle(3005);
    MXs_10->SetLineColor(kRed);
    MXs_10->SetFillColor(kRed);

    // draw plot and get p-value
    gStyle->SetOptStat(0);

    TCanvas* c_temp = new TCanvas("c", "", 600, 600); c_temp->cd();
    MXs_11->Draw("Hist"); MXs_10->Draw("HistSAME");
    TLegend* legend = gPad->BuildLegend(0.9, 0.9, 0.6, 0.6); legend->SetFillStyle(0); legend->SetLineWidth(0);
    c_temp->SaveAs("MXsll.png");

}
