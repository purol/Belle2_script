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

    ::ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Minimize"); // default: Minuit Migrad
    ::ROOT::Math::MinimizerOptions::SetDefaultStrategy(1); // default 1

    RooStats::UseNLLOffset(true); // default off

    const char* fname = "./PDFandDATA_workspace.root";

    TFile* f = TFile::Open(fname);

    RooWorkspace* w = (RooWorkspace*)f->Get("combined");

    // w->Print();
    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    OPTIONS* options = (OPTIONS*)malloc(sizeof(OPTIONS));
    Initialize_options(options, "none");

    w->loadSnapshot("NominalParamValues");
    FixParameters(w, options);

    RooDataSet* data = (RooDataSet*)w->data("obsData");

    // fit
    double eps = 0.1;
    RooAbsReal* nll;
    RooFitResult* fitres = MyMinimizeNLL(w, data, &nll, eps, false);

    RooRealVar* x_val_MXs1 = w->var("obs_x_channel_MXs1");
    RooRealVar* x_val_MXs2 = w->var("obs_x_channel_MXs2");
    RooRealVar* x_val_MXs3 = w->var("obs_x_channel_MXs3");
    std::unique_ptr<RooArgSet> params{model->getParameters(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3))};

    w->saveSnapshot("GlobalMinimumParamValues", *params, true);

    // get PLL value
    RooPlot* mu_frame;
    double PLL_value = -1;
    RooAbsReal* pll = nll->createProfile(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3));

    const int NStep = 100;
    for (int i = 0; i < NStep; i++) {
        for (int j = 0; j < NStep; j++) {
            for (int k = 0; k < NStep; k++) {

                w->loadSnapshot("GlobalMinimumParamValues");

                double mu1_local = -100.0 + i * (200.0 / NStep);
                double mu2_local = -100.0 + j * (200.0 / NStep);
                double mu3_local = -100.0 + k * (200.0 / NStep);

                x_val_MXs1->setVal(mu1_local);
                x_val_MXs2->setVal(mu2_local);
                x_val_MXs3->setVal(mu3_local);

                //x_val_MXs1->setConstant(true);
                //x_val_MXs2->setConstant(true);
                //x_val_MXs3->setConstant(true);

                //MyMinimizeNLLReuse(w, data, &nll, eps, false);

                PLL_value = pll->getVal();

                std::unique_ptr<RooArgSet> test_vars{model->getVariables()};
                test_vars->Print("v");

                printf("values: %lf, %lf %lf %lf\n", PLL_value, mu1_local, mu2_local, mu3_local);

            }
        }
    }

    return 0;
}
