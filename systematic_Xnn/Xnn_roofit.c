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

double LetsAdd(const char* dirname, RooRealVar* Mbc_, RooRealVar* weight_, RooDataSet* info_) {
    double Nevt = 0;
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
            *weight_ = total_weight;
            info_->add(RooArgSet(*Mbc_), weight_->getVal());
            Nevt = Nevt + total_weight;
        }
        input_file->Close();

    }

    return Nevt;
}

void Xnn_roofit(){

    // Observable:
    RooRealVar Mbc("Mbc", "Mbc", 5.2, 5.3);
    RooRealVar weight("weight", "weight", 0.0, 10.0);
    RooDataSet info("info", "info", RooArgSet(Mbc));

    const char* dirname = "./";
    LetsAdd(dirname, &Mbc, &weight, &info);
    RooPlot* Mbcframe = Mbc.frame(Bins(50), Title("M_{bc}^{sig} fit"));
    RooDataSet* d_Mbc = (RooDataSet*)info.reduce(RooArgSet(Mbc));

    // ARGUS background
    RooRealVar End("End", "End of Argus", 5.29);
    RooRealVar c("curvature", "curvature of ARGUS", -40, -100, 20);
    RooRealVar power("power", "power of ARGUS", 0.5);

    RooArgusBG pdfARGUS("ARGUS", "ARGUS PDF", Mbc, End, c, power);

    // crystal ball signal
    RooRealVar mean("mean", "mean of Crystal ball", 5.28);
    RooRealVar width("width", "width of Crystal ball", 0.003, 0.002, 0.004);
    RooRealVar alpha("alpha", "alpha of Crystal ball", 2.3, 1.9, 2.5);
    RooRealVar n("n", "n of Crystal ball", 1.0, 0.9, 1.1);

    RooCrystalBall pdfCRYSTAL("CRYSTAL", "CRYSTAL PDF", Mbc, mean, width, alpha, n);

    // Construct a signal and background PDF:
    RooRealVar nsig("nsig", "nsig", 20, 10, 130);
    RooRealVar nbkg("nbkg", "nkbkg", 1300, 800, 1800);
    RooExtendPdf esig("esignal", "extended signal p.d.f", pdfCRYSTAL, nsig);
    RooExtendPdf ebkg("ebkg", "extended background p.d.f", pdfARGUS, nbkg);

    RooAddPdf  model("model", "b+n", RooArgList(ebkg, esig));

    // fit
    model.fitTo(*d_Mbc, Extended(), SumW2Error(true));

    // Draw result
    d_Mbc->plotOn(Mbcframe, DataError(RooAbsData::SumW2));
    model.plotOn(Mbcframe, LineColor(kRed), Normalization(1.0, RooAbsReal::RelativeExpected));
    model.plotOn(Mbcframe, Components(pdfARGUS), LineColor(kBlue), LineStyle(kDashed), Normalization(1.0, RooAbsReal::RelativeExpected));
    model.plotOn(Mbcframe, Components(pdfCRYSTAL), LineColor(kViolet), LineStyle(kDashed), Normalization(1.0, RooAbsReal::RelativeExpected));

    TCanvas* c_temp = new TCanvas("Mbc", "Mbc", 600, 600);
    gPad->SetLeftMargin(0.15); Mbcframe->GetYaxis()->SetTitleOffset(1.4); Mbcframe->Draw();
    c_temp->SaveAs("fit_result.png");
}
