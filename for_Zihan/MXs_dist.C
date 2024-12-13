#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"

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

TH1D* MXs_11 = new TH1D("MXs_11", ";M_{Xs} [GeV];arbitrary unit", 30, 1.1, 4.0);
TH1D* MXs_10 = new TH1D("MXs_10", ";M_{Xs} [GeV];arbitrary unit", 30, 1.1, 4.0);

void LetsFill(const char* dirname, TH1D* MXs_hist, double weight_var = 1.0) {
    double Mx = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");

        tree_upsilon->SetBranchAddress("test", &Mx);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);

            MXs_hist->Fill(Mx, weight_var);
        }
        input_file->Close();

    }

}

void MXs_dist()
{

    // get data from root files
    const char* dirname = "";

    LetsFill(dirname, Oqq_SIGNAL_train, 1.0);
    LetsFill(dirname, OBB_BKG_train, 1.0);

    double factor = 1.0;

    // normalization
    Oqq_SIGNAL_train->Scale(factor / Oqq_SIGNAL_train->Integral(), "width");
    OBB_BKG_train->Scale(factor / OBB_BKG_train->Integral(), "width");

    Oqq_SIGNAL_train->SetFillStyle(3004);
    Oqq_SIGNAL_train->SetLineColor(kBlue);
    Oqq_SIGNAL_train->SetFillColor(kBlue);

    OBB_BKG_train->SetFillStyle(3005);
    OBB_BKG_train->SetLineColor(kRed);
    OBB_BKG_train->SetFillColor(kRed);

    // draw plot and get p-value
    gStyle->SetOptStat(0);

    TCanvas* c_temp = new TCanvas("c", "", 600, 600); c_temp->cd(); gPad->SetLogy();
    Oqq_SIGNAL_train->Draw("Hist"); OBB_BKG_train->Draw("HistSAME");
    TLegend* legend = gPad->BuildLegend(0.9, 0.9, 0.6, 0.6); legend->SetFillStyle(0); legend->SetLineWidth(0);
    c_temp->SaveAs("MXsll.png");

}
