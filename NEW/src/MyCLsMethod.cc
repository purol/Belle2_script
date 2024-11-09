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

const double step = 0.1;
const double boundary = 5.0;

const double BR_1 = 0.0000048514;
const double BR_2 = 0.0000085024;
const double BR_3 = 0.0000156653;
const double BR_total = 0.000029;
/* =========================== configuration =========================== */

double ConditionalFit(RooWorkspace* w, RooAbsReal** nll, double target_mu, double eps) {
    // just manually scan. Too slow. Is there any better way??

    double target_BR = target_mu * BR_total;

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

    double mu_MXs1_conditional = 0;
    double mu_MXs2_conditional = 0;
    double mu_MXs3_conditional = 0;
    double PLL_conditional = DBL_MAX;
    for (double mu1_local = mu_MXs1_global - boundary * mu_MXs1_err; mu1_local < mu_MXs1_global + boundary * mu_MXs1_err; mu1_local = mu1_local + mu_MXs1_err * step) {
        for (double mu2_local = mu_MXs2_global - boundary * mu_MXs2_err; mu2_local < mu_MXs2_global + boundary * mu_MXs2_err; mu2_local = mu2_local + mu_MXs2_err * step) {

            double mu3_local = (target_BR - mu1_local * BR_1 - mu2_local * BR_2) / BR_3;
            if ((mu3_local >= mu_MXs3_global - boundary * mu_MXs3_err) && (mu3_local < mu_MXs3_global + boundary * mu_MXs3_err)) {

                w->loadSnapshot("GlobalMinimumParamValues");

                mu_MXs1->setVal(mu1_local);
                mu_MXs2->setVal(mu2_local);
                mu_MXs3->setVal(mu3_local);

                PLL_value = pll->getVal();
                if ((PLL_conditional > PLL_value) && (PLL_value > 0.0)) {
                    PLL_conditional = pll->getVal();
                    mu_MXs1_conditional = mu1_local;
                    mu_MXs2_conditional = mu2_local;
                    mu_MXs3_conditional = mu3_local;
                }

            }

        }
    }

    return PLL_conditional;
}

double UnConditionalFit(RooWorkspace* w, RooAbsReal** nll, double target_mu, double eps, bool* OneSideFlag) {
    MyMinimizeNLLReuse(w, nll, eps, false);

    double mu_MXs1_conditional = w->var("mu_MXs1")->getVal();
    double mu_MXs2_conditional = w->var("mu_MXs2")->getVal();
    double mu_MXs3_conditional = w->var("mu_MXs3")->getVal();

    double mu_conditional = (mu_MXs1_conditional * BR_1 + mu_MXs2_conditional * BR_2 + mu_MXs3_conditional * BR_3) / target_BR;

    // SetOneSided(true) condition
    if (mu_conditional >= target_mu) (*OneSideFlag) = true; // do not conduct conditional fit!
    else (*OneSideFlag) = false; // please conduct conditional fit

    RooAbsReal::setHideOffset(false);
    return nll->getVal();
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
    RooStats::ModelConfig* model = (RooStats::ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooStats::RooSimultaneous* pdf = (RooStats::RooSimultaneous*)model->GetPdf();

    RooArgSet fGlobalObs = *model->GetGlobalObservables();
    RooArgSet fConditionalObs;

    // get data test statistics
    RooDataSet* data = (RooDataSet*)w->data("obsData");

    RooArgSet* allParams = pdf->getParameters(*data);
    RooStats::RemoveConstantParameters(allParams);

    RooAbsReal* nll = = pdf->createNLL(*data, RooFit::CloneData(false), RooFit::Constrain(*allParams), RooFit::GlobalObservables(fGlobalObs), RooFit::ConditionalObservables(fConditionalObs), RooFit::Offset("bin"));
    
    double data_test_statistic = 0.0;
    bool OneSideFlag = false;

    double data_uncond_NLL = UnConditionalFit(w, &nll, scanned_mu, eps, &OneSideFlag);
    w->saveSnapshot("GlobalMinimumParamValues", pdf->getVariables(), true);

    if (OneSideFlag) data_test_statistic = 0.0;
    else {
        w->loadSnapshot("GlobalMinimumParamValues");
        double data_cond_PLL = ConditionalFit(w, &nll, scanned_mu, eps);
        data_test_statistic = data_cond_PLL;
    }

    // get s+b test statistics
    w->loadSnapshot("GlobalMinimumParamValues");
    ConditionalFit(w, &nll, scanned_mu, eps);

    return 0;
}
