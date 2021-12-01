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

RooRealVar  Mbc_DATA("Mbc_DATA", "Mbc", 5.27, 5.29);
RooRealVar  Eecl_DATA("Eecl_DATA", "Eecl", 0, 1.2);
RooDataSet info_DATA("2Dinfo_DATA", "2Dinfo", RooArgSet(Mbc_DATA, Eecl_DATA));

RooRealVar  Mbc_MC_signal("Mbc_MC_signal", "Mbc", 5.27, 5.29);
RooRealVar  Eecl_MC_signal("Eecl_MC_signal", "Eecl", 0, 1.2);
RooDataSet info_MC_signal("2Dinfo_MC_signal", "2Dinfo", RooArgSet(Mbc_MC_signal, Eecl_MC_signal));

RooRealVar  Mbc_MC_background("Mbc_MC_background", "Mbc", 5.27, 5.29);
RooRealVar  Eecl_MC_background("Eecl_MC_background", "Eecl", 0, 1.2);
RooDataSet info_MC_background("2Dinfo_MC_background", "2Dinfo", RooArgSet(Mbc_MC_background, Eecl_MC_background));

void LetsAdd(const char* dirname, RooRealVar* Mbc_, RooRealVar*  Eecl_, RooDataSet* info_, double weight_ = 1.0) {
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
            info_->add(RooArgSet(*Mbc_, *Eecl_), weight_);
        }
        input_file->Close();

    }

}

void Signal_yield_fit()
{
    // to extract signal yield

    // get data from root files
    const char* MC_dirname_signal = "./output_after_TMVA_test/SIGNAL";
    LetsAdd(MC_dirname_signal, &Mbc_MC_signal, &Eecl_MC_signal, &info_MC_signal);

    const char* MC_dirname_CHG = "./output_after_TMVA_test/CHG";
    const char* MC_dirname_MIX = "./output_after_TMVA_test/MIX";
    const char* MC_dirname_UUBAR = "./output_after_TMVA_test/UUBAR";
    const char* MC_dirname_DDBAR = "./output_after_TMVA_test/DDBAR";
    const char* MC_dirname_SSBAR = "./output_after_TMVA_test/SSBAR";
    const char* MC_dirname_CHARM = "./output_after_TMVA_test/CHARM";
    LetsAdd(MC_dirname_CHG, &Mbc_MC_background, &Eecl_MC_background, &info_MC_background);
    LetsAdd(MC_dirname_MIX, &Mbc_MC_background, &Eecl_MC_background, &info_MC_background);
    LetsAdd(MC_dirname_UUBAR, &Mbc_MC_background, &Eecl_MC_background, &info_MC_background);
    LetsAdd(MC_dirname_DDBAR, &Mbc_MC_background, &Eecl_MC_background, &info_MC_background);
    LetsAdd(MC_dirname_SSBAR, &Mbc_MC_background, &Eecl_MC_background, &info_MC_background);
    LetsAdd(MC_dirname_CHARM, &Mbc_MC_background, &Eecl_MC_background, &info_MC_background);

    const char* DATA_dirname_BKG = "./output_after_TMVA_train/BKG";
    const char* DATA_dirname_SIGNAL = "./output_after_TMVA_train/SIGNAL";
    LetsAdd(DATA_dirname_BKG, &Mbc_DATA, &Eecl_DATA, &info_DATA);
    LetsAdd(DATA_dirname_SIGNAL, &Mbc_DATA, &Eecl_DATA, &info_DATA, 0.0);


    // define frame and get ready to make pdfs
    Eecl_DATA.setBins(12);
    RooPlot* Eeclframe = Eecl_DATA.frame(Bins(12), Title("Operations on binned datasets"));
    RooDataSet* d_Eecl = (RooDataSet*)info_DATA.reduce(RooArgSet(Eecl_DATA));

    Eecl_MC_signal.setBins(12);
    RooDataSet* dataset_Eecl_MC_signal = (RooDataSet*)info_MC_signal.reduce(RooArgSet(Eecl_MC_signal));
    RooDataHist* hist_Eecl_MC_signal = dataset_Eecl_MC_signal->binnedClone();
    Eecl_MC_background.setBins(12);
    RooDataSet* dataset_Eecl_MC_background = (RooDataSet*)info_MC_background.reduce(RooArgSet(Eecl_MC_background));
    RooDataHist* hist_Eecl_MC_background = dataset_Eecl_MC_background->binnedClone();


    // define pdf and extended pdf
    RooHistPdf histpdf_Eecl_signal("histpdf_Eecl_signal", "histpdf_Eecl_signal", Eecl_DATA, *hist_Eecl_MC_signal, 0);
    RooRealVar nsig("nsig", "number of signal events", 10, 0, 20);
    RooExtendPdf esig("esignal", "extended signal p.d.f", histpdf_Eecl_signal, nsig);

    RooHistPdf histpdf_Eecl_background("histpdf_Eecl_background", "histpdf_Eecl_background", Eecl_DATA, *hist_Eecl_MC_background, 0);
    RooRealVar nbkg("nbkg", "number of background events", 130, 0, 250);
    RooExtendPdf ebkg("ebkg", "extended background p.d.f", histpdf_Eecl_background, nbkg);

    RooAddPdf  totalpdf("model", "b+n", RooArgList(ebkg, esig));

    // fit
    totalpdf.fitTo(*d_Eecl);


    // Draw result
    d_Eecl->plotOn(Eeclframe);
    totalpdf.plotOn(Eeclframe, LineColor(kRed), Normalization(1.0, RooAbsReal::RelativeExpected));
    totalpdf.plotOn(Eeclframe, Components(esig), LineColor(kBlue), LineStyle(kDashed), Normalization(1.0, RooAbsReal::RelativeExpected));
    totalpdf.plotOn(Eeclframe, Components(ebkg), LineColor(kViolet), LineStyle(kDashed), Normalization(1.0, RooAbsReal::RelativeExpected));

    new TCanvas("Eecl", "Eecl", 600, 600);
    gPad->SetLeftMargin(0.15); Eeclframe->GetYaxis()->SetTitleOffset(1.4); Eeclframe->Draw();
}
