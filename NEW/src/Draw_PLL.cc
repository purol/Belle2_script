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

#include "RooChi2Var.h"
#include "RooAbsData.h"
#include "RooRealSumPdf.h"
#include "RooPoisson.h"
#include "RooGaussian.h"
#include "RooRealVar.h"
#include "RooMCStudy.h"
#include "RooMinuit.h"
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

#include "RooStats/ModelConfig.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/MinNLLTestStat.h"

#include "RooStats/HistFactory/FlexibleInterpVar.h"
#include "RooStats/HistFactory/PiecewiseInterpolation.h"
#include "RooStats/HistFactory/HistFactorySimultaneous.h"
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
    // argv[1]: {mu_MXs1|mu_MXs2|mu_MXs3}
    // argv[2]: MXs value

    ::ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2"); // default: Minuit Migrad
    ::ROOT::Math::MinimizerOptions::SetDefaultStrategy(1); // default 1

    const char* fname = "./PDFandDATA_workspace.root";

    TFile* f = TFile::Open(fname);

    RooWorkspace* w = (RooWorkspace*)f->Get("combined");

    w->Print();
    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    OPTIONS* options = (OPTIONS*)malloc(sizeof(OPTIONS));
    Initialize_options(options, "none");
    FixParameters(w, options);

    // test
    //RooRealVar* alpha = w->var("nom_gamma_stat_channel_bin_0");
    //printf("%lf", alpha->getValV());

    // Lets tell roofit the right names for our histogram variables //
    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x = (RooRealVar*)obs->find("obs_x_channel");
    x->SetTitle("FBDT output");
    x->setUnit("");

    // get Category and data
    RooCategory* idx = (RooCategory*)obs->find("channelCat");
    //RooAbsData* data = (RooAbsData*)w->data("obsData");
    RooDataSet* data = (RooDataSet*)w->data("asimovData");

    // fit
    //RooFitResult* fitres = model->fitTo(*data, RooFit::Minimizer("Minuit2"), RooFit::Extended(false), RooFit::Minos(RooArgSet(*w->var("mu_MXs1"), *w->var("mu_MXs2"), *w->var("mu_MXs3"))), RooFit::SumW2Error(false), Save());
    //RooAbsReal* nll;
    //ObtainNLL(w, data, &nll);
    //double eps = ::ROOT::Math::MinimizerOptions::DefaultTolerance();
    double eps = 0.001;
    RooAbsReal* nll;
    RooFitResult* fitres = MyMinimizeNLL(w, data, &nll, eps);

    // get PLL value
    RooPlot* mu_frame;
    RooAbsReal* pll;
    double set_value = std::atof(argv[2]);
    double PLL_value = -1;
    if (std::string(argv[1]) == std::string("mu_MXs1")) {
        RooRealVar* mu_MXs1 = w->var("mu_MXs1");
        mu_frame = mu_MXs1->frame();
        RooAbsReal* pll = nll->createProfile(*mu_MXs1);
        mu_MXs1->setVal(set_value);
        PLL_value = pll->getVal();
    }
    else if (std::string(argv[1]) == std::string("mu_MXs2")) {
        RooRealVar* mu_MXs2 = w->var("mu_MXs2");
        mu_frame = mu_MXs2->frame();
        RooAbsReal* pll = nll->createProfile(*mu_MXs2);
        mu_MXs2->setVal(set_value);
        PLL_value = pll->getVal();
    }
    else if (std::string(argv[1]) == std::string("mu_MXs3")) {
        RooRealVar* mu_MXs3 = w->var("mu_MXs3");
        mu_frame = mu_MXs3->frame();
        RooAbsReal* pll = nll->createProfile(*mu_MXs3);
        mu_MXs3->setVal(set_value);
        PLL_value = pll->getVal();
    }
    else {
        printf("choose among {mu_MXs1|mu_MXs2|mu_MXs3}\n");
        exit(1);
    }
    
    // save the result
    FILE* fp = fopen(("./" + std::string(argv[1]) + "_" + std::string(argv[2])).c_str(), "w");
    fprintf(fp, "%s %s %lf", argv[1], argv[2], PLL_value);
    fclose(fp);

    return 0;
}
