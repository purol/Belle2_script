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
    // argv[1]: total mu

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

    // constant variables
    const double BR_1 = 0.0000048514;
    const double BR_2 = 0.0000085024;
    const double BR_3 = 0.0000156653;
    const double BR_total = 0.000029;

    // target BR
    double target_mu = std::atof(argv[1]);
    double target_BR = target_mu * BR_total;

    // fit
    double eps = 0.1;
    RooAbsReal* nll;
    RooFitResult* fitres = MyMinimizeNLL(w, data, &nll, eps, false);

    RooRealVar* x_val_MXs1 = w->var("obs_x_channel_MXs1");
    RooRealVar* x_val_MXs2 = w->var("obs_x_channel_MXs2");
    RooRealVar* x_val_MXs3 = w->var("obs_x_channel_MXs3");
    std::unique_ptr<RooArgSet> params{model->getVariables()};

    w->saveSnapshot("GlobalMinimumParamValues", *params, true);

    // get global value
    RooRealVar* mu_MXs1 = w->var("mu_MXs1");
    RooRealVar* mu_MXs2 = w->var("mu_MXs2");
    RooRealVar* mu_MXs3 = w->var("mu_MXs3");
    const double mu_MXs1_global = w->var("mu_MXs1")->getVal();
    const double mu_MXs2_global = w->var("mu_MXs2")->getVal();
    const double mu_MXs3_global = w->var("mu_MXs3")->getVal();
    const double mu_MXs1_err = w->var("mu_MXs1")->getError();
    const double mu_MXs2_err = w->var("mu_MXs2")->getError();
    const double mu_MXs3_err = w->var("mu_MXs3")->getError();
    const double NLL_global = nll->getVal(); // global minimum of -log(L)
    //double PLL_value = -1;
    RooAbsReal* pll = nll->createProfile(RooArgSet(*mu_MXs1, *mu_MXs2, *mu_MXs3));

    FILE* fp = fopen(("scan_result_" + std::string(argv[1]) + ".csv").c_str(), "w");

    const double step = 0.01;
    double mu_MXs1_conditional = 0;
    double mu_MXs2_conditional = 0;
    double mu_MXs3_conditional = 0;
    double NLL_conditional = DBL_MAX;
    for (double mu1_local = mu_MXs1_global - 3 * mu_MXs1_err; mu1_local < mu_MXs1_global + 3 * mu_MXs1_err; mu1_local = mu1_local + mu_MXs1_err * step) {
        for (double mu2_local = mu_MXs2_global - 3 * mu_MXs2_err; mu2_local < mu_MXs2_global + 3 * mu_MXs2_err; mu2_local = mu2_local + mu_MXs2_err * step) {

            double mu3_local = (target_BR - mu1_local * BR_1 - mu2_local * BR_2) / BR_3;
            if ((mu3_local >= mu_MXs3_global - 3 * mu_MXs3_err) && (mu3_local < mu_MXs3_global + 3 * mu_MXs3_err)) {

                int k = std::atoi(argv[1]);

                w->loadSnapshot("GlobalMinimumParamValues");

                mu_MXs1->setVal(mu1_local);
                mu_MXs2->setVal(mu2_local);
                mu_MXs3->setVal(mu3_local);

                mu_MXs1->setConstant(true);
                mu_MXs2->setConstant(true);
                mu_MXs3->setConstant(true);

                RooFitResult* fitres = MyMinimizeNLLReuse(w, data, &nll, eps, false);

                /* test */
                //PLL_value = pll->getVal();

                std::unique_ptr<RooArgSet> test_vars{model->getVariables()};
                test_vars->Print("v");

                printf("status: %d\n", fitres->status());
                printf("PLL: %lf\n", nll->getVal() - NLL_global);

                w->loadSnapshot("GlobalMinimumParamValues");

                mu_MXs1->setVal(mu1_local);
                mu_MXs2->setVal(mu2_local);
                mu_MXs3->setVal(mu3_local);

                //mu_MXs1->setConstant(true);
                //mu_MXs2->setConstant(true);
                //mu_MXs3->setConstant(true);

                printf("ROOT PLL: %lf\n", pll->getVal());

                /* test */

                if (NLL_conditional > nll->getVal()) {
                    NLL_conditional = nll->getVal();
                    mu_MXs1_conditional = mu1_local;
                    mu_MXs2_conditional = mu2_local;
                    mu_MXs3_conditional = mu3_local;
                }

            }

        }
    }

    fprintf(fp, "%lf,%lf,%lf,%lf,%lf,%lf", target_mu, NLL_global, NLL_conditional, mu_MXs1_conditional, mu_MXs2_conditional, mu_MXs3_conditional);

    fclose(fp);

    return 0;
}
