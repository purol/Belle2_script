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
    // argv[1]: step index for mu_MXs3.

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
    std::unique_ptr<RooArgSet> params{model->getVariables()};

    w->saveSnapshot("GlobalMinimumParamValues", *params, true);

    // get PLL value
    RooRealVar* mu_MXs1 = w->var("mu_MXs1");
    RooRealVar* mu_MXs2 = w->var("mu_MXs2");
    RooRealVar* mu_MXs3 = w->var("mu_MXs3");
    const double mu_MXs1_global = w->var("mu_MXs1")->getVal();
    const double mu_MXs2_global = w->var("mu_MXs2")->getVal();
    const double mu_MXs3_global = w->var("mu_MXs3")->getVal();
    const double mu_MXs1_err = w->var("mu_MXs1")->getError();
    const double mu_MXs2_err = w->var("mu_MXs2")->getError();
    const double mu_MXs3_err = w->var("mu_MXs3")->getError();
    double PLL_value = -1;
    RooAbsReal* pll = nll->createProfile(RooArgSet(*mu_MXs1, *mu_MXs2, *mu_MXs3));

    FILE* fp = fopen(("scan_result_" + std::string(argv[1]) + ".csv").c_str(), "w");

    const int NStep = 30;
    for (int i = 0; i < NStep; i++) {
        for (int j = 0; j < NStep; j++) {
//            for (int k = 0; k < NStep; k++) {

                int k = std::atoi(argv[1]);

                w->loadSnapshot("GlobalMinimumParamValues");

                double mu1_local = mu_MXs1_global - 3 * mu_MXs1_err + i * (6 * mu_MXs1_err / NStep);
                double mu2_local = mu_MXs2_global - 3 * mu_MXs2_err + j * (6 * mu_MXs2_err / NStep);
                double mu3_local = mu_MXs3_global - 3 * mu_MXs3_err + k * (6 * mu_MXs3_err / NStep);

                mu_MXs1->setVal(mu1_local);
                mu_MXs2->setVal(mu2_local);
                mu_MXs3->setVal(mu3_local);

                //x_val_MXs1->setConstant(true);
                //x_val_MXs2->setConstant(true);
                //x_val_MXs3->setConstant(true);

                //MyMinimizeNLLReuse(w, data, &nll, eps, false);

                PLL_value = pll->getVal();

                //std::unique_ptr<RooArgSet> test_vars{model->getVariables()};
                //test_vars->Print("v");

                fprintf(fp, "%lf,%lf,%lf,%lf", PLL_value, mu_MXs1, mu_MXs2, mu_MXs3);

//            }
        }
    }

    fclose(fp);

    return 0;
}
