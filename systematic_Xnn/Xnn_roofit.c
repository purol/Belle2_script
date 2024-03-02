// last update: 2021-10-13
// for Belle2 data

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit;

void load_files(const char *dirname, std::vector<string>* names){
   TSystemDirectory dir(dirname, dirname);
   TList *files = dir.GetListOfFiles();
   if (files) {
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next())) {
         fname = file->GetName();
         if (!file->IsDirectory() && fname.EndsWith(".root")) {
            names->push_back(fname.Data());
         }
      }
   }
}

void LetsAdd(const char* dirname, RooRealVar* Mbc_, RooDataSet* info_) {
    double Mbc_var = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_Bsig->SetBranchAddress("Bsig_Mbc", &Mbc_var); // Mbc

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            double total_weight = 1.0;

            *Mbc_ = Mbc_var;
            info_->add(RooArgSet(*Mbc_), total_weight);
        }
        input_file->Close();

    }

}

void Xnn_roofit(){

    // Observable:
    RooRealVar  Mbc("Mbc", "Mbc", 5.2, 5.3);
    RooDataSet info("info", "info", RooArgSet(Mbc));

    const char* dirname = "./rootfile/";
    LetsAdd(dirname, &Mbc, &info);
    RooPlot* Mbcframe = Mbc.frame(Bins(50), Title("Operations on binned datasets"));
    RooDataSet* d_Mbc = (RooDataSet*)info.reduce(RooArgSet(Mbc));

    // ARGUS background
    RooRealVar End("End", "End of Argus", 5.29);
    RooRealVar c("curvature", "curvature of ARGUS", 1, 0.9, 1.1);
    RooRealVar power("power", "power of ARGUS", 0.5);

    RooArgusBG pdfARGUS("ARGUS", "ARGUS PDF", Mbc, End, c, power);

    // crystal ball signal
    RooRealVar mean("mean", "mean of Crystal ball", 5.28);
    RooRealVar width("width", "width of Crystal ball", 0.0027, 0.0024, 0.0030);
    RooRealVar alpha("alpha", "alpha of Crystal ball", 2.42, 2.30, 2.50);
    RooRealVar n("n", "n of Crystal ball", 1.117, 1.0, 1.2);

    RooCrystalBall pdfCRYSTAL("CRYSTAL", "CRYSTAL PDF", Mbc, mean, width, alpha, n);

    // Construct a signal and background PDF:
    RooRealVar nsig("nsig", "nsig", 10, 0, 20);
    RooRealVar nbkg("nbkg", "nkbkg", 800, 0, 10000);
    RooExtendPdf esig("esignal", "extended signal p.d.f", pdfCRYSTAL, nsig);
    RooExtendPdf ebkg("ebkg", "extended background p.d.f", pdfARGUS, nbkg);

    RooAddPdf  model("model", "b+n", RooArgList(ebkg, esig));

    // fit
    model.fitTo(d_Mbc, Extended());

    // Draw result
    d_Mbc->plotOn(Mbcframe);
    model.plotOn(Mbcframe, LineColor(kRed), Normalization(1.0, RooAbsReal::RelativeExpected));
    model.plotOn(Mbcframe, Components(pdfARGUS), LineColor(kBlue), LineStyle(kDashed), Normalization(1.0, RooAbsReal::RelativeExpected));
    model.plotOn(Mbcframe, Components(pdfCRYSTAL), LineColor(kViolet), LineStyle(kDashed), Normalization(1.0, RooAbsReal::RelativeExpected));

    new TCanvas("Mbc", "Mbc", 600, 600);
    gPad->SetLeftMargin(0.15); Mbcframe->GetYaxis()->SetTitleOffset(1.4); Mbcframe->Draw();

}
