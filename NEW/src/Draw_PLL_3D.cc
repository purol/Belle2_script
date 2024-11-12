#include <iostream>
#include <stdio.h>
#include <vector>
#include <cmath>

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

int findClosestIndex(std::vector<double> vec, double target) {
    int closestIndex = 0;
    double minDifference = std::abs(vec[0] - target);

    for (int i = 1; i < vec.size(); ++i) {
        double difference = std::abs(vec[i] - target);
        if (difference < minDifference) {
            minDifference = difference;
            closestIndex = i;
        }
    }

    return closestIndex;
}

typedef struct csv_setting_ {
    const char* filename;
    double step;
    double range;
    std::vector<double> mu_values;
    std::vector<double> mu1_values;
    std::vector<double> mu2_values;
    std::vector<double> mu3_values;

    std::vector<double> mu_meaningful_values;
    std::vector<double> mu1_meaningful_values;
    std::vector<double> mu2_meaningful_values;
    std::vector<double> mu3_meaningful_values;

} csv_setting;

const bool CSVfileON = false;
csv_setting csv_file = { "total_scan_heavy.csv", 0.035, 4 };

void ReadCSVFile() {
    FILE* fp = fopen(csv_file.filename, "r");

    double target_mu;
    double NLL_global;
    double PLL_conditional;
    double mu_MXs1_conditional;
    double mu_MXs2_conditional;
    double mu_MXs3_conditional;

    while (true) {
        if (fscanf(fp, "%lf,%lf,%lf,%lf,%lf,%lf\n", &target_mu, &NLL_global, &PLL_conditional, &mu_MXs1_conditional, &mu_MXs2_conditional, &mu_MXs3_conditional) == EOF) break;
        csv_file.mu_values.push_back(target_mu);
        csv_file.mu1_values.push_back(mu_MXs1_conditional);
        csv_file.mu2_values.push_back(mu_MXs2_conditional);
        csv_file.mu3_values.push_back(mu_MXs3_conditional);

        if (std::abs(PLL_conditional) < 100000000.0) {
            csv_file.mu_meaningful_values.push_back(target_mu);
            csv_file.mu1_meaningful_values.push_back(mu_MXs1_conditional);
            csv_file.mu2_meaningful_values.push_back(mu_MXs2_conditional);
            csv_file.mu3_meaningful_values.push_back(mu_MXs3_conditional);
        }

    }

    fclose(fp);
}

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

    RooMsgService::instance().setStreamStatus(1, false);
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

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

    // read csv file if ON
    if (CSVfileON) ReadCSVFile();

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
    double mu_MXs1_global = w->var("mu_MXs1")->getVal();
    double mu_MXs2_global = w->var("mu_MXs2")->getVal();
    double mu_MXs3_global = w->var("mu_MXs3")->getVal();
    double mu_MXs1_err = w->var("mu_MXs1")->getError();
    double mu_MXs2_err = w->var("mu_MXs2")->getError();
    double mu_MXs3_err = w->var("mu_MXs3")->getError();
    const double NLL_global = nll->getVal(); // global minimum of -log(L)
    double PLL_value = -1;
    RooAbsReal* pll = nll->createProfile(RooArgSet(*mu_MXs1, *mu_MXs2, *mu_MXs3));

    FILE* fp = fopen(("scan_result_" + std::string(argv[1]) + ".csv").c_str(), "w");

    const double step = 0.05;
    const double range = 4.0;
    double mu_MXs1_conditional = 0;
    double mu_MXs2_conditional = 0;
    double mu_MXs3_conditional = 0;
    double PLL_conditional = DBL_MAX;

    // do fancy thing if CSV ON
    if (CSVfileON) {
        int closest_index = findClosestIndex(csv_file.mu_meaningful_values, target_mu);

        // put values
        mu_MXs1_global = csv_file.mu1_meaningful_values.at(closest_index);
        mu_MXs2_global = csv_file.mu2_meaningful_values.at(closest_index);
        mu_MXs3_global = csv_file.mu3_meaningful_values.at(closest_index);

        /*
        * more good way is to use mu_MXs1_err * csv_file.step;
        * however, CSV file made with limited box space.
        * to find the real global minimum, we use mu_MXs1_err
        */
        mu_MXs1_err = mu_MXs1_err / range;
        mu_MXs2_err = mu_MXs2_err / range;
        mu_MXs3_err = mu_MXs3_err / range;
    }

    for (double mu1_local = mu_MXs1_global - range * mu_MXs1_err; mu1_local < mu_MXs1_global + range * mu_MXs1_err; mu1_local = mu1_local + mu_MXs1_err * step) {
        for (double mu2_local = mu_MXs2_global - range * mu_MXs2_err; mu2_local < mu_MXs2_global + range * mu_MXs2_err; mu2_local = mu2_local + mu_MXs2_err * step) {

            double mu3_local = (target_BR - mu1_local * BR_1 - mu2_local * BR_2) / BR_3;
            if ((mu3_local >= mu_MXs3_global - range * mu_MXs3_err) && (mu3_local < mu_MXs3_global + range * mu_MXs3_err)) {

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

    fprintf(fp, "%lf,%lf,%lf,%lf,%lf,%lf", target_mu, NLL_global, PLL_conditional, mu_MXs1_conditional, mu_MXs2_conditional, mu_MXs3_conditional);

    fclose(fp);

    return 0;
}
