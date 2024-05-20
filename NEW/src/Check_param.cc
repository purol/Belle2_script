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
#include "RooErrorVar.h"

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

void Drawpull(RooWorkspace* w, TIterator* iter) {
    std::vector<double> pulls;
    std::vector<double> pull_errors;
    std::vector<std::string> names;

    for (TObject* Tobj = iter->Next(); Tobj != 0; Tobj = iter->Next()) {
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(Tobj);
        std::string name = rrv->GetName();
        double val = rrv->getVal();
        double err = rrv->getError();

        std::cout.width(25);
        std::cout << name;
        std::cout.width(15);
        std::cout << val;
        std::cout.width(5);
        std::cout << "+-" << err << std::endl;

        if (name.find("alpha") != std::string::npos) {
            pulls.push_back((val - 0.0) / 1.0);
            pull_errors.push_back(err / 1.0);
            names.push_back(name);
        }
        else if (name.find("gamma_stat") != std::string::npos) {

            /* poisson */
            //RooRealVar* norm = w->var(("nom_" + names->at(i)).c_str());
            //std::poisson_distribution<> distribution(norm->getValV());
            //w->var(names->at(i).c_str())->setVal(distribution(generator) / norm->getValV());

            /* gaussian */
            RooRealVar* variable = w->var(name.c_str());
            RooErrorVar* err_variable = variable->errorVar();
            double width = err_variable->getValV();

            pulls.push_back((val - 1.0) / width);
            pull_errors.push_back(err / width);
            names.push_back(name);
        }
        else if ((name.find("gamma") != std::string::npos) && (name.find("uncorr") != std::string::npos)) {

            RooRealVar* variable = w->var(name.c_str());
            RooErrorVar* err_variable = variable->errorVar();
            double width = err_variable->getValV();

            pulls.push_back((val - 1.0) / width);
            pull_errors.push_back(err / width);
            names.push_back(name);

        }
    }

    for (unsigned int i = 0; i < names.size(); i++) {
        std::cout << names.at(i) << " " << pulls.at(i) << "+-" << pull_errors.at(i) << std::endl;
    }

    // draw pull
    int size_pull = pulls.size();

    TH1D* pull_ht = new TH1D("pull data hist", "pull of parameters;;", size_pull, 0, size_pull);
    for (int i = 0; i < size_pull; i++) {
        pull_ht->SetBinContent(i + 1, pulls.at(i));
        pull_ht->SetBinError(i + 1, pull_errors.at(i));
    }
    pull_ht->SetLineWidth(2.0);
    pull_ht->SetMarkerColor(1);
    pull_ht->SetMarkerStyle(21);
    pull_ht->SetLineColor(1);
    char** label_name = (char**)malloc(sizeof(char*) * size_pull);
    for (int i = 0; i < size_pull; i++) {
        label_name[i] = (char*)malloc(sizeof(char) * names.at(i).size() + 1);
        memcpy(label_name[i], names.at(i).c_str(), names.at(i).size() + 1);
    }

    TH1D* pull_one_sigma = new TH1D("1sig hist", "1sig;;", size_pull, 0.0, size_pull);
    for (int i = 0; i < size_pull; i++) {
        pull_one_sigma->SetBinContent(i + 1, 0.0);
        pull_one_sigma->SetBinError(i + 1, 1.0);
    }
    pull_one_sigma->SetFillColor(kGreen);
    pull_one_sigma->SetFillStyle(1001);

    TH1D* pull_two_sigma = new TH1D("2sig hist", "2sig;;", size_pull, 0.0, size_pull);
    for (int i = 0; i < size_pull; i++) {
        pull_two_sigma->SetBinContent(i + 1, 0.0);
        pull_two_sigma->SetBinError(i + 1, 2.0);
    }
    pull_two_sigma->SetFillColor(kYellow);
    pull_two_sigma->SetFillStyle(1001);
    for (int i = 0; i < size_pull; i++) {
        pull_two_sigma->GetXaxis()->SetBinLabel(i + 1, names.at(i).c_str());
    }
    pull_two_sigma->SetStats(false);
    pull_two_sigma->GetYaxis()->SetTitle("(#hat{#theta}-#theta)/#Delta#theta");
    pull_two_sigma->GetYaxis()->SetTitleOffset(1.4);
    pull_two_sigma->GetXaxis()->LabelsOption("v");

    TLine* line = new TLine(0.0, 0.0, size_pull, 0.0);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2); line->SetLineWidth(1);

    TCanvas* cpull = new TCanvas("pull_Plot", "pull Plot", 2200, 800); cpull->SetBottomMargin(0.3);
    pull_two_sigma->Draw("E2");
    pull_one_sigma->Draw("E2 same");
    pull_ht->Draw("e1 same");
    line->Draw();

    cpull->SaveAs("param_pull.png");

    delete cpull;
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
    Drawpull(w, iter);

    // define frame
    RooPlot* x_frame = x->frame(Title("fit result"));

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
