#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
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

RooRealVar  Mbc("Mbc", "Mbc", 5.27, 5.29);
RooRealVar  Eecl("Eecl", "Eecl", 0, 2);
RooDataSet info("2Dinfo", "2Dinfo", RooArgSet(Mbc, Eecl));

RooRealVar  Mbc_test("Mbc_test", "Mbc", 5.27, 5.29);
RooRealVar  Eecl_test("Eecl_test", "Eecl", 0, 2);
RooDataSet info_test("2Dinfo_test", "2Dinfo", RooArgSet(Mbc_test, Eecl_test));

void LetsAdd(const char* dirname, RooRealVar* Mbc_, RooRealVar*  Eecl_, RooDataSet* info_) {
    double Eecl_var = 0;
    double Mbc_var = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("roeEextra__bocleanMask__bc", &Eecl_var); // Eecl
        tree_Btag->SetBranchAddress("Btag_Mbc", &Mbc_var); // Mbc

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            
            *Eecl_ = Eecl_var;
            *Mbc_ = Mbc_var;
            info_->add(RooArgSet(*Mbc_, *Eecl_));
        }
        input_file->Close();

    }

}

void sideband_peaking_fit()
{
    // for sideband fitting
    // to determine mean and width of peaking background

    const char* dirname = "./rootfile/";

    LetsAdd(dirname, &Mbc, &Eecl, &info);

    
    RooPlot* Mbcframe = Mbc.frame(Bins(10), Title("Operations on binned datasets"));

    RooDataSet* d_Mbc = (RooDataSet*)info.reduce(RooArgSet(Mbc));


    // ARGUS (non-peaking component)
    RooRealVar End("End", "End of Argus", 5.29);
    RooRealVar c("curvature", "curvature of ARGUS", 1, 0.5, 1.2);
    RooRealVar power("power", "power of ARGUS", 0.5, 0.4, 0.53);

    RooArgusBG pdfARGUS("ARGUS", "ARGUS PDF", Mbc, End, c, power);


    // Gaussian (peaking component)
    RooRealVar mean("mean", "mean of gaussian", 5.28, 5.27, 5.29);
    RooRealVar sigma("sigma", "width of gaussian", 0.3, 0.2, 0.4);

    RooGaussian pdfGaussian("Gaussian", "Gaussian PDF", Mbc, mean, sigma);


    // total distribution
    RooRealVar Npeaking("Npeaking", "number of peaking events", 500, 0., 10000);
    RooRealVar Nnonpeaking("Nnonpeaking", "number of non-peaking events", 500, 0, 10000);
    RooAddPdf  model("model", "nonpeak + peak", RooArgList(pdfARGUS, pdfGaussian), RooArgList(Nnonpeaking, Npeaking));


    // fit
    model.fitTo(Mbc);


    // Draw result
    d_Mbc->plotOn(Mbcframe);
    model.plotOn(xframe, LineColor(kRed), Normalization(1.0, RooAbsReal::RelativeExpected));
    model.plotOn(xframe, Components(pdfARGUS), LineColor(kBlue), LineStyle(kDashed), Normalization(1.0, RooAbsReal::RelativeExpected));
    model.plotOn(xframe, Components(pdfGaussian), LineColor(kViolet), LineStyle(kDashed), Normalization(1.0, RooAbsReal::RelativeExpected));

    new TCanvas("Mbc", "Mbc", 600, 600);
    gPad->SetLeftMargin(0.15); Mbcframe->GetYaxis()->SetTitleOffset(1.4); Mbcframe->Draw();
}
