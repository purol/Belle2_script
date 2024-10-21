#include <iostream>
#include <stdio.h>

#include "TRandom3.h"
#include "TCanvas.h"
#include "TDatime.h"
#include "TStopwatch.h"
#include "TLegend.h"
#include "TIterator.h"
#include "TH3.h"
#include "TLatex.h"
#include "TLine.h"
#include "TFile.h"

#include "RooChi2Var.h"
#include "RooAbsData.h"
#include "RooRealSumPdf.h"
#include "RooPoisson.h"
#include "RooGaussian.h"
#include "RooRealVar.h"
#include "RooMCStudy.h"
//#include "RooMinuit.h"
#include "RooCategory.h"
#include "RooHistPdf.h"
#include "RooSimultaneous.h"
#include "RooExtendPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "RooParamHistFunc.h"
#include "RooHist.h"
#include "RooRandom.h"
#include "RooErrorVar.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/MinNLLTestStat.h"

#include "RooStats/HistFactory/FlexibleInterpVar.h"
#include "RooStats/HistFactory/PiecewiseInterpolation.h"
//#include "RooStats/HistFactory/HistFactorySimultaneous.h"
#include "RooStats/HistFactory/Channel.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/ParamHistFunc.h"
#include "RooStats/HistFactory/HistFactoryModelUtils.h"
#include "RooStats/HistFactory/RooBarlowBeestonLL.h"

#include "template.h"
#include "Fitter.h"

using namespace RooFit;
using namespace RooStats;
using namespace HistFactory;

double weight_KIDsys[RarityBins * 7] = { 0.0 };
double weight_PIDsys[RarityBins * 7] = { 0.0 };
double weight_BRsys[RarityBins * 3] = { 0.0 };
double weight_pi0sys[RarityBins * 7] = { 0.0 };
double weight_FEIsys[RarityBins * 3] = { 0.0 };

std::vector<std::string> split(std::string str, char Delimiter) {
    std::istringstream iss(str);
    std::string buffer;

    std::vector<std::string> result;

    while (getline(iss, buffer, Delimiter)) {
        result.push_back(buffer);
    }

    return result;
}

int main() {

    ::ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Minimize"); // default: Minuit Migrad
    ::ROOT::Math::MinimizerOptions::SetDefaultStrategy(1); // default 1

    RooStats::UseNLLOffset(true); // default off

    const char* fname = "./PDFandDATA_workspace.root";

    TFile* f = TFile::Open(fname);

    RooWorkspace* w = (RooWorkspace*)f->Get("combined");

    w->Print();
    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    OPTIONS* options = (OPTIONS*)malloc(sizeof(OPTIONS));
    Initialize_options(options, "none");

    w->loadSnapshot("NominalParamValues");
    FixParameters(w, options);

    // test
    //RooRealVar* alpha = w->var("nom_gamma_stat_channel_bin_0");
    //printf("%lf", alpha->getValV());

    // get Category and data
    //RooAbsData* data = (RooAbsData*)w->data("obsData");
    RooDataSet* data = (RooDataSet*)w->data("asimovData");

    // fit
    //RooFitResult* fitres = model->fitTo(*data, RooFit::Minimizer("Minuit2"), RooFit::Extended(false), RooFit::Minos(RooArgSet(*w->var("mu_MXs1"), *w->var("mu_MXs2"), *w->var("mu_MXs3"))), RooFit::SumW2Error(false), Save());
    //double eps = ::ROOT::Math::MinimizerOptions::DefaultTolerance();
    double eps = 0.001;
    RooAbsReal* nll;
    RooFitResult* fitres = MyMinimizeNLL(w, data, &nll, eps);

    RooArgSet fitargs = fitres->floatParsFinal();
    TIterator* iter(fitargs.createIterator());

    for (TObject* a = iter->Next(); a != 0; a = iter->Next()) {
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(a);
        std::string name = rrv->GetName();
        double val = rrv->getVal();
        double err = rrv->getError();
        double HIerr = rrv->getAsymErrorHi();
        double LOerr = rrv->getAsymErrorLo();

        printf("fit result for %s = %lf +- %lf\n", name.c_str(), val, err);
        printf("MINOS error: %lf %lf\n", HIerr, LOerr);

    }

    // draw pull
    iter->Reset();
    Drawpull(w, iter, "param_pull_type0.png", 0);
    Drawpull(w, iter, "param_pull_type1.png", 1);

    // define frame
    // RooPlot* x_frame = x->frame(Title("fit result"));

    // draw
    // GetPlotTemplate(w, data);

    // draw profile likelihood
    /*
    * it is obsolete. Please use `Draw_PLL.cc`
    * 
    RooRealVar* mu_MXs1 = w->var("mu_MXs1");
    mu_MXs1->setRange(mu_MXs1->getValV() - 4, mu_MXs1->getValV() + 4);

    RooRealVar* mu_MXs2 = w->var("mu_MXs2");
    mu_MXs2->setRange(mu_MXs2->getValV() - 4, mu_MXs2->getValV() + 4);

    RooRealVar* mu_MXs3 = w->var("mu_MXs3");
    mu_MXs3->setRange(mu_MXs3->getValV() - 4, mu_MXs3->getValV() + 4);

    RooPlot* mu_frame = mu_MXs3->frame();
    RooAbsReal* pll = nll->createProfile(*mu_MXs3);
    pll->plotOn(mu_frame, RooFit::Precision(-1));

    TCanvas* cmu = new TCanvas("pllPlot", "pllPlot", 700, 700);
    mu_frame->Draw();
    cmu->SaveAs("profile_likelihood_mu.png");
    delete cmu;
    */

    return 0;
}
