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

RooRealVar  Mbc_DATA("Mbc", "Mbc_DATA", 5.27, 5.29);
RooRealVar  Eecl_DATA("Eecl", "Eecl_DATA", 0, 1.2);
RooRealVar weight_DATA("weight", "weight_DATA", 0.0, 1.0);
RooDataSet info_DATA("2Dinfo", "2Dinfo_DATA", RooArgSet(Mbc_DATA, Eecl_DATA, weight_DATA), WeightVar("weight"));

RooRealVar  Mbc_MC_signal("Mbc", "Mbc_MC_signal", 5.27, 5.29);
RooRealVar  Eecl_MC_signal("Eecl", "Eecl_MC_signal", 0, 1.2);
RooRealVar weight_MC_signal("weight", "weight_MC_signal", 0.0, 1.0);
RooDataSet info_MC_signal("2Dinfo", "2Dinfo_MC_signal", RooArgSet(Mbc_MC_signal, Eecl_MC_signal, weight_MC_signal), WeightVar("weight"));

RooRealVar  Mbc_MC_background("Mbc", "Mbc_MC_background", 5.27, 5.29);
RooRealVar  Eecl_MC_background("Eecl", "Eecl_MC_background", 0, 1.2);
RooRealVar weight_MC_background("weight", "weight_MC_background", 0.0, 1.0);
RooDataSet info_MC_background("2Dinfo", "2Dinfo_MC_background", RooArgSet(Mbc_MC_background, Eecl_MC_background, weight_MC_background), WeightVar("weight"));

void LetsAdd(const char* dirname, RooRealVar* Mbc_, RooRealVar*  Eecl_, RooRealVar* weight_, RooDataSet* info_, double weight_var = 1.0) {
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
            *weight_ = weight_var;
            info_->add(RooArgSet(*Mbc_, *Eecl_), weight_->getVal());
        }
        input_file->Close();

    }

}

void Signal_yield_fit()
{
    // to extract signal yield
    RooRealVar EeclFit("Eecl", "Eecl", 0, 1.2, "GeV");
    EeclFit.setBins(12);
    RooPlot* Eeclframe = EeclFit.frame(Bins(12), Title(" "));

    // get data from root files
    const char* MC_dirname_signal = "./output_after_TMVA_test/SIGNAL";
    LetsAdd(MC_dirname_signal, &Mbc_MC_signal, &Eecl_MC_signal, &weight_MC_signal, &info_MC_signal);

    const char* MC_dirname_CHG = "./output_after_TMVA_test/CHG";
    const char* MC_dirname_MIX = "./output_after_TMVA_test/MIX";
    const char* MC_dirname_UUBAR = "./output_after_TMVA_test/UUBAR";
    const char* MC_dirname_DDBAR = "./output_after_TMVA_test/DDBAR";
    const char* MC_dirname_SSBAR = "./output_after_TMVA_test/SSBAR";
    const char* MC_dirname_CHARM = "./output_after_TMVA_test/CHARM";
    LetsAdd(MC_dirname_CHG, &Mbc_MC_background, &Eecl_MC_background, &weight_MC_background, &info_MC_background);
    LetsAdd(MC_dirname_MIX, &Mbc_MC_background, &Eecl_MC_background, &weight_MC_background, &info_MC_background);
    LetsAdd(MC_dirname_UUBAR, &Mbc_MC_background, &Eecl_MC_background, &weight_MC_background, &info_MC_background);
    LetsAdd(MC_dirname_DDBAR, &Mbc_MC_background, &Eecl_MC_background, &weight_MC_background, &info_MC_background);
    LetsAdd(MC_dirname_SSBAR, &Mbc_MC_background, &Eecl_MC_background, &weight_MC_background, &info_MC_background);
    LetsAdd(MC_dirname_CHARM, &Mbc_MC_background, &Eecl_MC_background, &weight_MC_background, &info_MC_background);

    const char* DATA_dirname_BKG = "./output_after_TMVA_test/BKG";
    const char* DATA_dirname_SIGNAL = "./output_after_TMVA_test/SIGNAL";
    LetsAdd(DATA_dirname_BKG, &Mbc_DATA, &Eecl_DATA, &weight_DATA, &info_DATA);
    LetsAdd(DATA_dirname_SIGNAL, &Mbc_DATA, &Eecl_DATA, &weight_DATA, &info_DATA, 0.003385);


    // define frame and get ready to make pdfs
    Eecl_DATA.setBins(12);
    //RooPlot* Eeclframe = Eecl_DATA.frame(Bins(12), Title("Operations on binned datasets"));
    RooDataSet* d_Eecl = (RooDataSet*)info_DATA.reduce(RooArgSet(Eecl_DATA));

    Eecl_MC_signal.setBins(12);
    RooDataSet* dataset_Eecl_MC_signal = (RooDataSet*)info_MC_signal.reduce(RooArgSet(Eecl_MC_signal));
    RooDataHist hist_Eecl_MC_signal("hist_Eecl_MC_signal", "histogram for Eecl of MC signal samples", EeclFit, *dataset_Eecl_MC_signal);
    Eecl_MC_background.setBins(12);
    RooDataSet* dataset_Eecl_MC_background = (RooDataSet*)info_MC_background.reduce(RooArgSet(Eecl_MC_background));
    RooDataHist hist_Eecl_MC_background("hist_Eecl_MC_background", "histogram for Eecl of MC background samples", EeclFit, *dataset_Eecl_MC_background);


    // define pdf and extended pdf
    RooHistPdf histpdf_Eecl_signal("histpdf_Eecl_signal", "histpdf_Eecl_signal", EeclFit, hist_Eecl_MC_signal, 0);
    RooRealVar nsig("nsig", "number of signal events", 10, -100, 100);
    RooExtendPdf esig("esignal", "extended signal p.d.f", histpdf_Eecl_signal, nsig);

    RooHistPdf histpdf_Eecl_background("histpdf_Eecl_background", "histpdf_Eecl_background", EeclFit, hist_Eecl_MC_background, 0);
    RooRealVar nbkg("nbkg", "number of background events", 130, -300, 300);
    RooExtendPdf ebkg("ebkg", "extended background p.d.f", histpdf_Eecl_background, nbkg);

    RooAddPdf  totalpdf("model", "b+n", RooArgList(ebkg, esig));

    // fit
    totalpdf.fitTo(*d_Eecl);

    //info_DATA.get(0)->Print("V");
    //info_DATA.get(1)->Print("V");

    // Draw result
    d_Eecl->plotOn(Eeclframe);
    totalpdf.plotOn(Eeclframe, LineColor(kRed));
    totalpdf.plotOn(Eeclframe, Components(esig), LineColor(kBlue), LineStyle(kDashed));
    totalpdf.plotOn(Eeclframe, Components(ebkg), LineColor(kViolet), LineStyle(kDashed));
    //dataset_Eecl_MC_signal->plotOn(Eeclframe);
    //histpdf_Eecl_signal.plotOn(Eeclframe);
    //totalpdf.plotOn(Eeclframe, LineColor(kRed), Normalization(1.0, RooAbsReal::RelativeExpected));

    TCanvas* c = new TCanvas("Eecl", "Eecl", 600, 600);
    gPad->SetLeftMargin(0.15); Eeclframe->GetYaxis()->SetTitleOffset(1.4); Eeclframe->Draw(); c->SaveAs("Eecl_distribution.png");



    /* ============== toy MC study ============== */
    RooRealVar  Eecl_TOY("Eecl", "Eecl_TOY", 0, 1.2);
    Eecl_TOY.setBins(12);

    RooMCStudy* mcstudy = new RooMCStudy(totalpdf, Eecl_TOY, Binned(kTRUE), Silence(), Extended(),FitOptions(Save(kTRUE), PrintEvalErrors(0)));
    mcstudy->generateAndFit(1000);

    // Make plots of the distributions of mean, the error on mean and the pull of mean
    RooPlot* frame1 = mcstudy->plotParam(nsig, Bins(40));
    RooPlot* frame2 = mcstudy->plotError(nsig, Bins(40));
    RooPlot* frame3 = mcstudy->plotPull(nsig, Bins(40), FitGauss(kTRUE));


    // Draw all plots on a canvas
    gStyle->SetOptStat(0);
    TCanvas* cf = new TCanvas("rf801_mcstudy", "rf801_mcstudy", 1200, 400);
    cf->Divide(3, 1);
    cf->cd(1); gPad->SetLeftMargin(0.15); frame1->GetYaxis()->SetTitleOffset(1.4); frame1->Draw();
    cf->cd(2); gPad->SetLeftMargin(0.15); frame2->GetYaxis()->SetTitleOffset(1.4); frame2->Draw();
    cf->cd(3); gPad->SetLeftMargin(0.15); frame3->GetYaxis()->SetTitleOffset(1.4); frame3->Draw();
}
