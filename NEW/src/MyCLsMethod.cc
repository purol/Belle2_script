#include <iostream>
#include <stdio.h>
#include <string>
#include <algorithm>
#include <vector>
#include <numeric>
#include <time.h>
#include <random>
#include <sstream>

#include "TRandom3.h"
#include "TCanvas.h"
#include "TDatime.h"
#include "TStopwatch.h"
#include "TLegend.h"
#include "TIterator.h"
#include "TH3.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"

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
#include "RooMsgService.h"
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

#include "RooStats/FrequentistCalculator.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterPlot.h"

#include "Fitter.h"
#include "constants.h"

std::random_device rd;
std::default_random_engine generator(rd());

/* =========================== configuration =========================== */
const double eps = 0.1;

const double step = 0.2;
const double boundary = 3.0;

const double BR_1 = 0.0000048514;
const double BR_2 = 0.0000085024;
const double BR_3 = 0.0000156653;
const double BR_total = 0.000029;
/* =========================== configuration =========================== */

typedef struct MyFitResult_ {
    double mu1_value;
    double mu2_value;
    double mu3_value;

    double mu1_error;
    double mu2_error;
    double mu3_error;

    double status;

    bool OneSideFlag;

    double mu_value;
} MyFitResult;

double ConditionalFit(RooWorkspace* w, RooAbsReal** nll, double target_mu, double eps, MyFitResult MyGlobalFitResult, const char* snapshot_name, MyFitResult* MyconditionalFitResult) {
    // just manually scan. Too slow. Is there any better way??

    double target_BR = target_mu * BR_total;

    RooRealVar* mu_MXs1 = w->var("mu_MXs1");
    RooRealVar* mu_MXs2 = w->var("mu_MXs2");
    RooRealVar* mu_MXs3 = w->var("mu_MXs3");
    const double mu_MXs1_global = MyGlobalFitResult.mu1_value;
    const double mu_MXs2_global = MyGlobalFitResult.mu2_value;
    const double mu_MXs3_global = MyGlobalFitResult.mu3_value;
    const double mu_MXs1_err = MyGlobalFitResult.mu1_error;
    const double mu_MXs2_err = MyGlobalFitResult.mu2_error;
    const double mu_MXs3_err = MyGlobalFitResult.mu3_error;

    // my ansatz: start from intersection point between \mu=const plane and perpendicular line from global minimum
    const double A_ = BR_1 / BR_total;
    const double B_ = BR_2 / BR_total;
    const double C_ = BR_3 / BR_total;
    const double D_ = -target_mu;
    const double mu_MXs1_initial = mu_MXs1_global - A_ * (A_ * mu_MXs1_global + B_ * mu_MXs2_global + C_ * mu_MXs3_global + D_) / (A_ * A_ + B_ * B_ + C_ * C_);
    const double mu_MXs2_initial = mu_MXs2_global - B_ * (A_ * mu_MXs1_global + B_ * mu_MXs2_global + C_ * mu_MXs3_global + D_) / (A_ * A_ + B_ * B_ + C_ * C_);
    const double mu_MXs3_initial = mu_MXs3_global - C_ * (A_ * mu_MXs1_global + B_ * mu_MXs2_global + C_ * mu_MXs3_global + D_) / (A_ * A_ + B_ * B_ + C_ * C_);

    double PLL_value = -1;
    RooAbsReal* pll = (*nll)->createProfile(RooArgSet(*mu_MXs1, *mu_MXs2, *mu_MXs3));

    double mu_MXs1_conditional = 0;
    double mu_MXs2_conditional = 0;
    double mu_MXs3_conditional = 0;
    double PLL_conditional = DBL_MAX;
    for (double mu1_local = mu_MXs1_initial - boundary * mu_MXs1_err; mu1_local < mu_MXs1_initial + boundary * mu_MXs1_err; mu1_local = mu1_local + mu_MXs1_err * step) {
        for (double mu2_local = mu_MXs2_initial - boundary * mu_MXs2_err; mu2_local < mu_MXs2_initial + boundary * mu_MXs2_err; mu2_local = mu2_local + mu_MXs2_err * step) {

            double mu3_local = (target_BR - mu1_local * BR_1 - mu2_local * BR_2) / BR_3;
            if ((mu3_local >= mu_MXs3_initial - boundary * mu_MXs3_err) && (mu3_local < mu_MXs3_initial + boundary * mu_MXs3_err)) {

                w->loadSnapshot(snapshot_name);

                mu_MXs1->setVal(mu1_local);
                mu_MXs2->setVal(mu2_local);
                mu_MXs3->setVal(mu3_local);

                PLL_value = pll->getVal();
                if ((PLL_conditional > PLL_value) && (PLL_value > 0.0)) {
                    PLL_conditional = pll->getVal();
                    mu_MXs1_conditional = mu1_local;
                    mu_MXs2_conditional = mu2_local;
                    mu_MXs3_conditional = mu3_local;

                    MyconditionalFitResult->mu1_value = mu_MXs1_conditional;
                    MyconditionalFitResult->mu2_value = mu_MXs2_conditional;
                    MyconditionalFitResult->mu3_value = mu_MXs3_conditional;
                    MyconditionalFitResult->mu1_error = mu_MXs1_err;
                    MyconditionalFitResult->mu2_error = mu_MXs2_err;
                    MyconditionalFitResult->mu3_error = mu_MXs3_err;
                    MyconditionalFitResult->status = PLL_conditional;
                    MyconditionalFitResult->OneSideFlag = true;
                    MyconditionalFitResult->mu_value = target_mu;
                }

            }

        }
    }

    return PLL_conditional;
}

double UnConditionalFit(RooWorkspace* w, RooAbsReal** nll, double target_mu, double eps, MyFitResult* MyUnconditionalFitResult) {
    double target_BR = target_mu * BR_total;

    RooFitResult* fitres = MyMinimizeNLLReuse(w, nll, eps, false);

    double mu_MXs1_conditional = w->var("mu_MXs1")->getVal();
    double mu_MXs2_conditional = w->var("mu_MXs2")->getVal();
    double mu_MXs3_conditional = w->var("mu_MXs3")->getVal();

    double mu_MXs1_error_conditional = w->var("mu_MXs1")->getError();
    double mu_MXs2_error_conditional = w->var("mu_MXs2")->getError();
    double mu_MXs3_error_conditional = w->var("mu_MXs3")->getError();

    double mu_conditional = (mu_MXs1_conditional * BR_1 + mu_MXs2_conditional * BR_2 + mu_MXs3_conditional * BR_3) / target_BR;

    MyUnconditionalFitResult->mu1_value = mu_MXs1_conditional;
    MyUnconditionalFitResult->mu2_value = mu_MXs2_conditional;
    MyUnconditionalFitResult->mu3_value = mu_MXs3_conditional;

    MyUnconditionalFitResult->mu1_error = mu_MXs1_error_conditional;
    MyUnconditionalFitResult->mu2_error = mu_MXs2_error_conditional;
    MyUnconditionalFitResult->mu3_error = mu_MXs3_error_conditional;

    MyUnconditionalFitResult->mu_value = mu_conditional;

    MyUnconditionalFitResult->status = fitres->status();

    // SetOneSided(true) condition
    if (mu_conditional >= target_mu) MyUnconditionalFitResult->OneSideFlag = true; // do not conduct conditional fit!
    else MyUnconditionalFitResult->OneSideFlag = false; // please conduct conditional fit

    RooAbsReal::setHideOffset(false);
    return (*nll)->getVal();
}

double GenerateToy() {

}

double FluctuateGlobalObservables() {

}

int main(int argc, char* argv[]) {
    /*
    * This code is a homemade CLs method. This code is based on Frequentist calculator and SetOneSided(true)
    * argv[1]: scanned \mu
    * argv[2]: Ntoy for s+b hypothesis
    * argv[3]: Ntoy for b only hypothesis
    */

    // arguments
    const double scanned_mu = std::atof(argv[1]);
    const int Nsig = std::atoi(argv[2]);
    const int Nbkg = std::atoi(argv[3]);

    // initial setting
    ::ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(2);
    ::ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Minimize"); // default: Minuit Migrad
    ::ROOT::Math::MinimizerOptions::SetDefaultStrategy(1); // default 1

    RooStats::UseNLLOffset(true); // default off

    RooMsgService::instance().setStreamStatus(1, false);
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

    RooRandom::randomGenerator()->SetSeed(rd());

    // get model and pdf
    const char* fname = "./PDFandDATA_workspace.root";
    TFile* f = TFile::Open(fname);
    RooWorkspace* w = (RooWorkspace*)f->Get("combined");

    RooStats::ModelConfig* model = (RooStats::ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* pdf = (RooSimultaneous*)model->GetPdf();

    RooArgSet fGlobalObs = *model->GetGlobalObservables();
    RooArgSet fConditionalObs;

    // get data test statistics
    RooDataSet* data = (RooDataSet*)w->data("obsData");

    RooArgSet* allParams = pdf->getParameters(*data);
    RooStats::RemoveConstantParameters(allParams);

    RooAbsReal* nll = = pdf->createNLL(*data, RooFit::CloneData(false), RooFit::Constrain(*allParams), RooFit::GlobalObservables(fGlobalObs), RooFit::ConditionalObservables(fConditionalObs), RooFit::Offset("bin"));
    
    double data_test_statistic = 0.0;
    MyFitResult MyUnconditionalFitResult;

    double data_uncond_NLL = UnConditionalFit(w, &nll, scanned_mu, eps, &MyUnconditionalFitResult);
    std::unique_ptr<RooArgSet> params{pdf->getVariables()};
    w->saveSnapshot("GlobalMinimumParamValues", *params, true);

    if (MyUnconditionalFitResult.OneSideFlag) data_test_statistic = 0.0;
    else {
        w->loadSnapshot("GlobalMinimumParamValues");
        MyFitResult MyconditionalFitResult;
        double data_cond_PLL = ConditionalFit(w, &nll, scanned_mu, eps, MyUnconditionalFitResult, "GlobalMinimumParamValues", &MyconditionalFitResult);
        data_test_statistic = data_cond_PLL;

        printf("%lf %lf %lf %lf\n", MyconditionalFitResult.mu1_value, MyconditionalFitResult.mu2_value, MyconditionalFitResult.mu3_value, MyconditionalFitResult.status);
    }

    // get s+b test statistics
    //w->loadSnapshot("GlobalMinimumParamValues");
    //ConditionalFit(w, &nll, scanned_mu, eps);

    return 0;
}
