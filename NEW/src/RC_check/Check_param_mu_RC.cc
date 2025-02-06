#include <iostream>
#include <stdio.h>
#include <string>
#include <cstring>

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

int main(int argc, char* argv[]) {
    /*
    * argv[1]: version (ex. v020)
    * argv[2]: output path
    */

    ::ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Minimize"); // default: Minuit Migrad
    ::ROOT::Math::MinimizerOptions::SetDefaultStrategy(1); // default 1

    RooStats::UseNLLOffset(true); // default off

    const char* fname = strdup(("./KumoiRD_" + std::string(argv[1]) + "/set1/PDFandDATA_workspace_mu.root").c_str());

    TFile* f = TFile::Open(fname);

    RooWorkspace* w = (RooWorkspace*)f->Get("combined");

    w->Print();
    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    OPTIONS* options = (OPTIONS*)malloc(sizeof(OPTIONS));
    Initialize_options(options, "none", ("./KumoiRD_" + std::string(argv[1]) + "/set1/").c_str());
    w->loadSnapshot("NominalParamValues");
    FixParameters(w, options);

    PrintNevtFile(w, ("./KumoiRD_" + std::string(argv[1]) + "/set1/Nominal_Nevt_mu.txt").c_str());

    // test
    //RooRealVar* alpha = w->var("nom_gamma_stat_channel_bin_0");
    //printf("%lf", alpha->getValV());

    // get Category and data
    RooDataSet* data = (RooDataSet*)w->data("obsData");
    //RooDataSet* data = (RooDataSet*)w->data("asimovData");

    // draw
    GetPlotTemplate(w, data, ("./KumoiRD_" + std::string(argv[1]) + "/set1/prefit_mu.png").c_str());

    // fit
    double eps = 0.1;
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

    FILE* fp;
    fp = fopen(("./KumoiRD_" + std::string(argv[1]) + "/set1/mu_result_mu.csv").c_str(), "w");
    for (TObject* a = iter->Next(); a != 0; a = iter->Next()) {
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(a);
        std::string name = rrv->GetName();
        double val = rrv->getVal();
        double err = rrv->getError();
        double HIerr = rrv->getAsymErrorHi();
        double LOerr = rrv->getAsymErrorLo();

        if (name == "mu") fprintf(fp, "%lf,%lf,%lf,", val, HIerr, LOerr);
    }
    for (TObject* a = iter->Next(); a != 0; a = iter->Next()) {
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(a);
        std::string name = rrv->GetName();
        double val = rrv->getVal();
        double err = rrv->getError();
        double HIerr = rrv->getAsymErrorHi();
        double LOerr = rrv->getAsymErrorLo();

        if (name == "mu_MXs2") fprintf(fp, "%lf,%lf,%lf,", val, HIerr, LOerr);
    }
    for (TObject* a = iter->Next(); a != 0; a = iter->Next()) {
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(a);
        std::string name = rrv->GetName();
        double val = rrv->getVal();
        double err = rrv->getError();
        double HIerr = rrv->getAsymErrorHi();
        double LOerr = rrv->getAsymErrorLo();

        if (name == "mu_MXs3") fprintf(fp, "%lf,%lf,%lf\n", val, HIerr, LOerr);
    }
    fclose(fp);

    // draw pull
    iter->Reset();
    Drawpull(w, iter, ("./KumoiRD_" + std::string(argv[1]) + "/set1/param_pull_type0_mu.png").c_str(), 0);
    Drawpull(w, iter, ("./KumoiRD_" + std::string(argv[1]) + "/set1/param_pull_type1_mu.png").c_str(), 1);

    // define frame
    // RooPlot* x_frame = x->frame(Title("fit result"));

    // draw
    GetPlotTemplate(w, data, ("./KumoiRD_" + std::string(argv[1]) + "/set1/postfit_mu.png").c_str());
    double NSIGNAL = GetNumEvts(w, "Signal_MX1") + GetNumEvts(w, "Signal_MX2") + GetNumEvts(w, "Signal_MX3");
    double NBKG = GetNumEvts(w, "CHG_MX1") + GetNumEvts(w, "CHG_MX2") + GetNumEvts(w, "CHG_MX3") +
GetNumEvts(w, "MIX_MX1") + GetNumEvts(w, "MIX_MX2") + GetNumEvts(w, "MIX_MX3") +
GetNumEvts(w, "UUBAR_MX1") + GetNumEvts(w, "UUBAR_MX2") + GetNumEvts(w, "UUBAR_MX3") +
GetNumEvts(w, "DDBAR_MX1") + GetNumEvts(w, "DDBAR_MX2") + GetNumEvts(w, "DDBAR_MX3") +
GetNumEvts(w, "SSBAR_MX1") + GetNumEvts(w, "SSBAR_MX2") + GetNumEvts(w, "SSBAR_MX3") +
GetNumEvts(w, "CHARM_MX1") + GetNumEvts(w, "CHARM_MX2") + GetNumEvts(w, "CHARM_MX3");

    printf("NSIG: %lf\nNBKG: %lf\n", NSIGNAL, NBKG);

    // get stat uncertainty only
    Initialize_options(options, "all", ("./KumoiRD_" + std::string(argv[1]) + "/set1/").c_str());
    FixParameters(w, options);

    fitres = MyMinimizeNLL(w, data, &nll, eps);

    fitargs = fitres->floatParsFinal();
    TIterator* iter_again(fitargs.createIterator());

    for (TObject* a = iter_again->Next(); a != 0; a = iter_again->Next()) {
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(a);
        std::string name = rrv->GetName();
        double val = rrv->getVal();
        double err = rrv->getError();
        double HIerr = rrv->getAsymErrorHi();
        double LOerr = rrv->getAsymErrorLo();

        printf("fit result for %s = %lf +- %lf\n", name.c_str(), val, err);
        printf("MINOS error: %lf %lf\n", HIerr, LOerr);

    }

    PrintNevtFile(w, ("./KumoiRD_" + std::string(argv[1]) + "/set1/Fit_Nevt_mu.txt").c_str());
    PrintNuisanceParameters(&fitargs);

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
