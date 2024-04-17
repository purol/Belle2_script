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
#include "RooMsgService.h"
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

#include "RooStats/FrequentistCalculator.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterPlot.h"

#include "Fitter_indep_norm.h"

using namespace RooFit;
using namespace RooStats;
using namespace HistFactory;

using std::cout;
using std::endl;

# define MyDEBUG false
void Debug(RooWorkspace* w, RooFitResult* fitres, RooDataSet* data);

#include "template.h"
#include "correctors.h"

int Toy_iter_num = 0.0;
int LT_iter_num = 0.0;

std::random_device rd;
std::default_random_engine generator(rd());

void GetNameOfParams(RooWorkspace* w, std::vector<std::string>* names) {
    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    w->loadSnapshot("ParamValues");

    RooAbsData* data = (RooAbsData*)w->data("asimovData");
    RooFitResult* fitres = model->fitTo(*data, RooFit::Extended(false), RooFit::SumW2Error(false), PrintLevel(-1), Save());

    RooArgSet fitargs = fitres->floatParsFinal();
    TIterator* iter(fitargs.createIterator());

    for (TObject* a = iter->Next(); a != 0; a = iter->Next()) {
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(a);
        std::string name = rrv->GetName();
        names->push_back(name);
    }

    w->loadSnapshot("ParamValues");
}

class FileSaver {
private:
    double* m_true_param;
    double* m_fitting_param;
    double* m_fitting_param_error;
    double* m_fitting_param_HIerror;
    double* m_fitting_param_LOerror;
    int m_covQual;
    int m_status;
    double m_edm;

    TFile* m_file;
    TTree* m_tree;

public:
    FileSaver();
    void OpenFile(bool IsItToy, std::vector<std::string>* names, double mu_injected, int indicator);
    void CloseFile();
    void GetTrueValues(RooWorkspace* w, std::vector<std::string>* names);
    void GetFittingValues(RooFitResult* fitres, std::vector<std::string>* names);
    void GetFittingStatus(RooFitResult* fitres);
    void WriteIntoBranch();
};

FileSaver::FileSaver() {
    m_true_param = nullptr;
    m_fitting_param = nullptr;
    m_fitting_param_error = nullptr;
    m_fitting_param_HIerror = nullptr;
    m_fitting_param_LOerror = nullptr;
    m_covQual = -1;
    m_status = -1;
    m_edm = -1;
    TFile* m_file = nullptr;
    TTree* m_tree = nullptr;
}

void FileSaver::OpenFile(bool IsItToy, std::vector<std::string>* names, double mu_injected, int indicator) {

    // malloc params
    m_true_param = (double*)malloc(sizeof(double) * names->size());
    m_fitting_param = (double*)malloc(sizeof(double) * names->size());
    m_fitting_param_error = (double*)malloc(sizeof(double) * names->size());
    m_fitting_param_HIerror = (double*)malloc(sizeof(double) * names->size());
    m_fitting_param_LOerror = (double*)malloc(sizeof(double) * names->size());

    // open file
    if (IsItToy) {
        m_file = new TFile(("./out/TOY_result_" + std::to_string(indicator) + ".root").c_str(), "recreate");
        m_file->cd();
        m_tree = new TTree("TOY_result", "");
    }
    else {
        m_file = new TFile(("./out/LT_result_" + std::to_string(mu_injected) + "_" + std::to_string(indicator) + ".root").c_str(), "recreate");
        m_file->cd();
        m_tree = new TTree("LT_result", "");
    }

    // set branches
    for (unsigned int i = 0; i < names->size(); i++) {
        m_tree->Branch((names->at(i) + "_true").c_str(), &m_true_param[i]);
        m_tree->Branch((names->at(i) + "_value").c_str(), &m_fitting_param[i]);
        m_tree->Branch((names->at(i) + "_error").c_str(), &m_fitting_param_error[i]);
        m_tree->Branch((names->at(i) + "_HIerror").c_str(), &m_fitting_param_HIerror[i]);
        m_tree->Branch((names->at(i) + "_LOerror").c_str(), &m_fitting_param_LOerror[i]);
    }
    m_tree->Branch("covQual", &m_covQual);
    m_tree->Branch("status", &m_status);
    m_tree->Branch("edm", &m_edm);
}

void FileSaver::CloseFile() {
    m_file->cd();
    m_tree->Write();
    m_file->Close();
}

void FileSaver::GetTrueValues(RooWorkspace* w, std::vector<std::string>* names) {

    // write into member variables
    for (unsigned int i = 0; i < names->size(); i++) {
        m_true_param[i] = w->var(names->at(i).c_str())->getValV();
    }

}

void FileSaver::GetFittingValues(RooFitResult* fitres, std::vector<std::string>* names) {

    RooArgSet fitargs = fitres->floatParsFinal();
    TIterator* iter(fitargs.createIterator());

    for (TObject* a = iter->Next(); a != 0; a = iter->Next()) {
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(a);
        std::string name = rrv->GetName();
        double val = rrv->getVal();
        double err = rrv->getError();
        double HIerr = rrv->getAsymErrorHi();
        double LOerr = rrv->getAsymErrorLo();

        for (unsigned int i = 0; i < names->size(); i++) {
            if (name == names->at(i)) {
                m_fitting_param[i] = val;
                m_fitting_param_error[i] = err;
                m_fitting_param_HIerror[i] = HIerr;
                m_fitting_param_LOerror[i] = LOerr;
            }
        }

    }

}

void FileSaver::GetFittingStatus(RooFitResult* fitres) {

    m_covQual = fitres->covQual();
    m_status = fitres->status();
    m_edm = fitres->edm();

}

void FileSaver::WriteIntoBranch() {

    m_tree->Fill();

}

FileSaver filesaver;

std::vector<std::string> split(std::string str, char Delimiter) {
    std::istringstream iss(str);
    std::string buffer;

    std::vector<std::string> result;

    while (getline(iss, buffer, Delimiter)) {
        result.push_back(buffer);
    }

    return result;
}

double SetParamsForToy(RooWorkspace* w, std::vector<std::string>* names, double injected_mu) {

    double Nevt = 0.0;

    w->loadSnapshot("ParamValues");
    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();
    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x_val = w->var("obs_x_channel");

    for (unsigned int i = 0; i < names->size(); i++) {
        if (names->at(i).find("alpha") != std::string::npos) {
            std::normal_distribution<double> distribution(0.0, 1.0);
            w->var(names->at(i).c_str())->setVal(distribution(generator));
        }
        else if (names->at(i).find("gamma_stat") != std::string::npos) {

            /* poisson */
            //RooRealVar* norm = w->var(("nom_" + names->at(i)).c_str());
            //std::poisson_distribution<> distribution(norm->getValV());
            //w->var(names->at(i).c_str())->setVal(distribution(generator) / norm->getValV());

            /* gaussian */
            RooRealVar* variable = w->var(names->at(i).c_str());
            RooErrorVar* err_variable = variable->errorVar();

            std::normal_distribution<double> distribution(1.0, err_variable->getValV());
            w->var(names->at(i).c_str())->setVal(distribution(generator));
        }
        else if ( (names->at(i).find("gamma") != std::string::npos) && (names->at(i).find("uncorr") != std::string::npos) ) {

            RooRealVar* variable = w->var(names->at(i).c_str());
            RooErrorVar* err_variable = variable->errorVar();

            std::normal_distribution<double> distribution(1.0, err_variable->getValV());
            w->var(names->at(i).c_str())->setVal(distribution(generator));

        }
    }

    w->var("mu_MXs1")->setVal(injected_mu);
    w->var("mu_MXs2")->setVal(injected_mu);
    w->var("mu_MXs3")->setVal(injected_mu);

    /* ================================ cal Nexpected ================================*/
    RooAbsBinning const& binning = x_val->getBinning();
    const double oldVal = x_val->getVal();

    for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
        double binCenter = binning.binCenter(iBin);
        double binWidth = binning.binWidth(iBin);

        *x_val = binCenter; // set x value

        for (unsigned int j = 0; j < Sample_names.size(); j++) {
            RooAbsReal* temp_func = w->function(Sample_names.at(j).c_str());
            Nevt = Nevt + temp_func->getValV();
            if (temp_func->getValV() < 0) {
                printf("[ERROR] negative count!\n");
                exit(1);
            }
        }

    }

    *x_val = oldVal;

    return Nevt;

}

void MyToyMCStudy(RooWorkspace *w, std::vector<std::string>* names, double eps, int indicator = 0){

        ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
        RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

        RooArgSet* obs = (RooArgSet*)mc->GetObservables();
        RooRealVar* x = (RooRealVar*)obs->find("obs_x_channel");

        for(int i=0; i< Toy_iter_num; i++) { // Do Toy MC study
            
            double Nevt_total = SetParamsForToy(w, names, 1.0);

            filesaver.GetTrueValues(w, names);

            RooDataSet* genData = model->generate(RooArgSet(*x,model->indexCat()), Nevt_total, false, true, "", false, true);

            w->loadSnapshot("ParamValues");
            //RooFitResult* fitres = model->fitTo(*genData, RooFit::Extended(false), RooFit::Minos(RooArgSet(*w->var("mu_MXs1"), *w->var("mu_MXs2"), *w->var("mu_MXs3"))), RooFit::SumW2Error(false), Save());
            RooAbsReal* nll;
            RooFitResult* fitres = MyMinimizeNLL(w, genData, &nll, eps);

            if(MyDEBUG) Debug(w, fitres, genData);

            filesaver.GetFittingValues(fitres, names);
            filesaver.GetFittingStatus(fitres);
            filesaver.WriteIntoBranch();

            delete fitres;

        }
}

void MyLinearityTest(RooWorkspace* w, std::vector<std::string>* names, double mu_injected, double eps, int indicator = 0) {

    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x = (RooRealVar*)obs->find("obs_x_channel");

    for (int i = 0; i < LT_iter_num; i++) { // Do LT MC study

        double Nevt_total = SetParamsForToy(w, names, mu_injected);

        filesaver.GetTrueValues(w, names);

        RooDataSet* genData = model->generate(RooArgSet(*x, model->indexCat()), Nevt_total, false, true, "", false, true);
        w->loadSnapshot("ParamValues");

        //RooFitResult* fitres = model->fitTo(*genData, RooFit::Extended(false), RooFit::Minos(RooArgSet(*w->var("mu_MXs1"), *w->var("mu_MXs2"), *w->var("mu_MXs3"))), RooFit::SumW2Error(false), Save());
        RooAbsReal* nll;
        RooFitResult* fitres = MyMinimizeNLL(w, genData, &nll, eps);

        filesaver.GetFittingValues(fitres, names);
        filesaver.GetFittingStatus(fitres);
        filesaver.WriteIntoBranch();

        delete fitres;

    }
}

void FitToData(RooWorkspace* w, double eps) {

    w->loadSnapshot("ParamValues");

    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();
    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x_val = w->var("obs_x_channel");

    // get Category and data
    RooCategory* idx = (RooCategory*)obs->find("channelCat");
    RooDataSet* data = (RooDataSet*)w->data("asimovData");

    // fit
    //RooFitResult* fitres = model->fitTo(*data, RooFit::Extended(false), RooFit::Minos(RooArgSet(*w->var("mu_MXs1"), *w->var("mu_MXs2"), *w->var("mu_MXs3"))), RooFit::SumW2Error(false), Save());
    RooAbsReal* nll;
    RooFitResult* fitres = MyMinimizeNLL(w, data, &nll, eps);

    // get expected num of evts for PDFs
    double Signal_Nevts = GetNumEvts(w, "Signal_MX1") + GetNumEvts(w, "Signal_MX2") + GetNumEvts(w, "Signal_MX3");
    double CHG_Nevts = GetNumEvts(w, "CHG_MX1") + GetNumEvts(w, "CHG_MX2") + GetNumEvts(w, "CHG_MX3");
    double MIX_Nevts = GetNumEvts(w, "MIX_MX1") + GetNumEvts(w, "MIX_MX2") + GetNumEvts(w, "MIX_MX3");
    double UUBAR_Nevts = GetNumEvts(w, "UUBAR_MX1") + GetNumEvts(w, "UUBAR_MX2") + GetNumEvts(w, "UUBAR_MX3");
    double DDBAR_Nevts = GetNumEvts(w, "DDBAR_MX1") + GetNumEvts(w, "DDBAR_MX2") + GetNumEvts(w, "DDBAR_MX3");
    double SSBAR_Nevts = GetNumEvts(w, "SSBAR_MX1") + GetNumEvts(w, "SSBAR_MX2") + GetNumEvts(w, "SSBAR_MX3");
    double CHARM_Nevts = GetNumEvts(w, "CHARM_MX1") + GetNumEvts(w, "CHARM_MX2") + GetNumEvts(w, "CHARM_MX3");

    // draw
    RooPlot* x_frame = x_val->frame(Title("FBDT"));
    data->plotOn(x_frame, DataError(RooAbsData::Poisson), Cut("channelCat==0"), DrawOption("ZP"), Name("data"));
    //data->plotOn(x_frame, DataError(RooAbsData::Poisson), Cut("channelCat==0"), MarkerSize(0.4), DrawOption("ZP"), Normalization(1, RooAbsReal::ScaleType::NumEvent));
    model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kRed - 6), LineWidth(0), Components("L_x_*Signal*_StatUncert,L_x_*CHG*_StatUncert,L_x_*MIX*_StatUncert,L_x_*UUBAR*_StatUncert,L_x_*DDBAR*_StatUncert,L_x_*SSBAR*_StatUncert,L_x_*CHARM*_StatUncert"), Normalization(Signal_Nevts + CHG_Nevts + MIX_Nevts + UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("signal"));
    model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kBlue - 6), LineWidth(0), Components("L_x_*CHG*_StatUncert,L_x_*MIX*_StatUncert,L_x_*UUBAR*_StatUncert,L_x_*DDBAR*_StatUncert,L_x_*SSBAR*_StatUncert,L_x_*CHARM*_StatUncert"), Normalization(CHG_Nevts + MIX_Nevts + UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("Charged B"));
    model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kCyan - 6), LineWidth(0), Components("L_x_*MIX*_StatUncert,L_x_*UUBAR*_StatUncert,L_x_*DDBAR*_StatUncert,L_x_*SSBAR*_StatUncert,L_x_*CHARM*_StatUncert"), Normalization(MIX_Nevts + UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("Neutral B"));
    model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kOrange - 6), LineWidth(0), Components("L_x_*UUBAR*_StatUncert,L_x_*DDBAR*_StatUncert,L_x_*SSBAR*_StatUncert,L_x_*CHARM*_StatUncert"), Normalization(UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("Continuum"));
    data->plotOn(x_frame, DataError(RooAbsData::Poisson), Cut("channelCat==0"), DrawOption("ZP"), Name("data"));

    TCanvas* canvas = new TCanvas("sPlot", "sPlot demo", 700, 700);

    // draw PDFs
    x_frame->Draw();

    // draw legend
    TLegend* leg = new TLegend(0.9, 0.9, 0.7, 0.7);
    leg->SetFillStyle(0);
    leg->SetLineWidth(0);
    leg->AddEntry("data", "Data");
    leg->AddEntry("signal", "B #rightarrow X_{s} #nu #bar{#nu}");
    leg->AddEntry("Charged B", "Charged B");
    leg->AddEntry("Neutral B", "Neutral B");
    leg->AddEntry("Continuum", "q#bar{q}");
    leg->Draw();

    canvas->SaveAs("fit_plot.png");

    // Get fitting variables
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

    w->loadSnapshot("ParamValues");

    delete canvas;
    delete fitres;
}

void MyToyMCStudyDataPoisson(RooWorkspace* w, std::vector<std::string>* names, double eps) {

    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x = (RooRealVar*)obs->find("obs_x_channel");

    // construct PDfs from data
    RooAbsData* data = (RooAbsData*)w->data("asimovData");
    TH1D* HistData = (TH1D*)data->createHistogram("HistData", *x, Binning(RarityBins, BinMIN, BinMAX));
    RooDataHist DataHist("DataHist", "DataHist", *x, HistData);
    RooHistPdf dataPDF("dataPDF", "dataPDF", RooArgSet(*x, model->indexCat()), DataHist);

    // get total number of event in data
    double Nevt_data = 0;
    int Nevt_data_int = 0;
    for (int i = 0; i < RarityBins; i++) Nevt_data = Nevt_data + HistData->GetBinContent(i + 1);
    Nevt_data_int = (int)std::round(Nevt_data);

    
    for (int i = 0; i < Toy_iter_num; i++) { // Do Toy MC study

        w->loadSnapshot("ParamValues");

        filesaver.GetTrueValues(w, names);

        RooDataSet* genData = dataPDF.generate(RooArgSet(*x, model->indexCat()), Nevt_data_int, false, true, "", false, true);

        w->loadSnapshot("ParamValues");
        //RooFitResult* fitres = model->fitTo(*genData, RooFit::Extended(false), RooFit::Minos(RooArgSet(*w->var("mu_MXs1"), *w->var("mu_MXs2"), *w->var("mu_MXs3"))), RooFit::SumW2Error(false), Save());
        RooAbsReal* nll;
        RooFitResult* fitres = MyMinimizeNLL(w, genData, &nll, eps);

        filesaver.GetFittingValues(fitres, names);
        filesaver.GetFittingStatus(fitres);
        filesaver.WriteIntoBranch();

        delete fitres;

    }
    
}

void Debug(RooWorkspace* w, RooFitResult* fitres, RooDataSet* data) {
    // deprecated. Do not use it before fix the function
    RooArgSet fitargs = fitres->floatParsFinal();
    double fit_status = fitres->status();

    TIterator* iter(fitargs.createIterator());

    for (TObject* a = iter->Next(); a != 0; a = iter->Next()) {
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(a);
        std::string name = rrv->GetName();
        double val = rrv->getVal();
        double err = rrv->getError();

        if (name == "mu_MXs1") {
            if (std::abs(fit_status) > 0.5) {

                ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
                RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();
                RooArgSet* obs = (RooArgSet*)mc->GetObservables();
                RooRealVar* x_val = w->var("obs_x_channel");

                // get Category and data
                RooCategory* idx = (RooCategory*)obs->find("channelCat");

                // get expected num of evts for PDFs
                double Signal_Nevts = GetNumEvts(w, "Signal_MX1") + GetNumEvts(w, "Signal_MX2") + GetNumEvts(w, "Signal_MX3");
                double CHG_Nevts = GetNumEvts(w, "CHG_MX1") + GetNumEvts(w, "CHG_MX2") + GetNumEvts(w, "CHG_MX3");
                double MIX_Nevts = GetNumEvts(w, "MIX_MX1") + GetNumEvts(w, "MIX_MX2") + GetNumEvts(w, "MIX_MX3");
                double UUBAR_Nevts = GetNumEvts(w, "UUBAR_MX1") + GetNumEvts(w, "UUBAR_MX2") + GetNumEvts(w, "UUBAR_MX3");
                double DDBAR_Nevts = GetNumEvts(w, "DDBAR_MX1") + GetNumEvts(w, "DDBAR_MX2") + GetNumEvts(w, "DDBAR_MX3");
                double SSBAR_Nevts = GetNumEvts(w, "SSBAR_MX1") + GetNumEvts(w, "SSBAR_MX2") + GetNumEvts(w, "SSBAR_MX3");
                double CHARM_Nevts = GetNumEvts(w, "CHARM_MX1") + GetNumEvts(w, "CHARM_MX2") + GetNumEvts(w, "CHARM_MX3");

                // draw
                RooPlot* x_frame = x_val->frame(Title("FBDT"));
                data->plotOn(x_frame, DataError(RooAbsData::Poisson), Cut("channelCat==0"), DrawOption("ZP"), Name("data"));
                //data->plotOn(x_frame, DataError(RooAbsData::Poisson), Cut("channelCat==0"), MarkerSize(0.4), DrawOption("ZP"), Normalization(1, RooAbsReal::ScaleType::NumEvent));
                model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kRed - 6), LineWidth(0), Components("L_x_*Signal*_StatUncert,L_x_*CHG*_StatUncert,L_x_*MIX*_StatUncert,L_x_*UUBAR*_StatUncert,L_x_*DDBAR*_StatUncert,L_x_*SSBAR*_StatUncert,L_x_*CHARM*_StatUncert"), Normalization(Signal_Nevts + CHG_Nevts + MIX_Nevts + UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("signal"));
                model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kBlue - 6), LineWidth(0), Components("L_x_*CHG*_StatUncert,L_x_*MIX*_StatUncert,L_x_*UUBAR*_StatUncert,L_x_*DDBAR*_StatUncert,L_x_*SSBAR*_StatUncert,L_x_*CHARM*_StatUncert"), Normalization(CHG_Nevts + MIX_Nevts + UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("Charged B"));
                model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kCyan - 6), LineWidth(0), Components("L_x_*MIX*_StatUncert,L_x_*UUBAR*_StatUncert,L_x_*DDBAR*_StatUncert,L_x_*SSBAR*_StatUncert,L_x_*CHARM*_StatUncert"), Normalization(MIX_Nevts + UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("Neutral B"));
                model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kOrange - 6), LineWidth(0), Components("L_x_*UUBAR*_StatUncert,L_x_*DDBAR*_StatUncert,L_x_*SSBAR*_StatUncert,L_x_*CHARM*_StatUncert"), Normalization(UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("Continuum"));
                data->plotOn(x_frame, DataError(RooAbsData::Poisson), Cut("channelCat==0"), DrawOption("ZP"), Name("data"));

                TCanvas* canvas = new TCanvas("sPlot", "sPlot demo", 700, 700);

                // draw PDFs
                x_frame->Draw();

                // draw legend
                TLegend* leg = new TLegend(0.9, 0.9, 0.7, 0.7);
                leg->SetFillStyle(0);
                leg->SetLineWidth(0);
                leg->AddEntry("data", "Data");
                leg->AddEntry("signal", "B #rightarrow X_{s} #nu #bar{#nu}");
                leg->AddEntry("Charged B", "Charged B");
                leg->AddEntry("Neutral B", "Neutral B");
                leg->AddEntry("Continuum", "q#bar{q}");
                leg->Draw();

                canvas->SaveAs( (std::string("fit_plot_") + "mu-" + std::to_string(val) + "_" + "err-" + std::to_string(err) + "_" + std::to_string(fit_status) + ".png").c_str() );

                delete canvas;

            }
        }

    }
}

int main(int argc, char* argv[]) {

    ::ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit"); // default: Minuit Migrad
    ::ROOT::Math::MinimizerOptions::SetDefaultStrategy(1); // default 1

    // RooStats::UseNLLOffset(true);

    RooMsgService::instance().setStreamStatus(1, false);
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

    RooRandom::randomGenerator()->SetSeed(rd());

    // argv[1]: {ToyMC|LinearityTest}
    // argv[2]: injected mu when Linearity test
    // argv[3]: eps for minimizer
    // argv[4]: indicator
    // argv[5]: Num of sample
    // argv[6]: fixed nuisance parameter
    double injected_mu = -1;
    double eps = -1.0;
    int indicator = 0;
    std::string fixed_param;

    std::vector<std::string> param_names;
    
    if (std::string(argv[1]) == std::string("ToyMC")) {  // main ToyMC
        if (argc == 5) {
            injected_mu = -1;
            eps = std::atof(argv[2]);
            indicator = std::atoi(argv[3]);
            Toy_iter_num = std::atoi(argv[4]);
        }
        else {
            printf("Toy MC requires 4 arguments {eps} {indicator} {Num of sample}\n");
            exit(1);
        }
    }
    else if (std::string(argv[1]) == std::string("LinearityTest")) { // main LinearityTest 12 1234 0.8
        if (argc == 6) {
            injected_mu = std::atof(argv[2]);
            eps = std::atof(argv[3]);
            indicator = std::atoi(argv[4]);
            LT_iter_num = std::atoi(argv[5]);
        }
        else {
            printf("Linearity test requires 5 arguments: {injected mu} {eps} {indicator} {Num of sample}\n");
            exit(1);
        }
    }
    else if (std::string(argv[1]) == std::string("nuisance")) { // check the power of nuisance parameters
        if (argc == 6) {
            injected_mu = -1;
            eps = std::atof(argv[2]);
            indicator = std::atoi(argv[3]);
            Toy_iter_num = std::atoi(argv[4]);
            fixed_param = std::string(argv[5]);
        }
        else {
            printf("nuisance Toy MC requires 5 arguments {eps} {indicator} {Num of sample} {fixed param}\n");
            exit(1);
        }
    }
    else {
        printf("first argument should be {ToyMC|LinearityTest|nuisance}\n");
        exit(1);
    }

    const char* fname = "./PDFandDATA_workspace.root";

    TFile* f = TFile::Open(fname);

    RooWorkspace* w = (RooWorkspace*)f->Get("combined");

    if (std::string(argv[1]) == std::string("ToyMC")) {
        ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
        RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

        w->loadSnapshot("NominalParamValues");

        RooRealVar* x_val = w->var("obs_x_channel");
        std::unique_ptr<RooArgSet> params{model->getParameters(*x_val)};

        // save snapshot
        w->saveSnapshot("ParamValues", *params, true);
    }
    else if (std::string(argv[1]) == std::string("LinearityTest")) {
        ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
        RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

        w->loadSnapshot("NominalParamValues");

        RooRealVar* x_val = w->var("obs_x_channel");
        std::unique_ptr<RooArgSet> params{model->getParameters(*x_val)};

        // save snapshot
        w->saveSnapshot("ParamValues", *params, true);
    }
    else if (std::string(argv[1]) == std::string("nuisance")) {
        OPTIONS* options = (OPTIONS*)malloc(sizeof(OPTIONS));
        Initialize_options(options, fixed_param.c_str());
        FixParameters(w, options);
    }

    GetNameOfParams(w, &param_names);

    if (std::string(argv[1]) == std::string("ToyMC")) {
        filesaver.OpenFile(true, &param_names, injected_mu, indicator);
        MyToyMCStudy(w, &param_names, eps, indicator);
    }
    else if (std::string(argv[1]) == std::string("LinearityTest")) {
        filesaver.OpenFile(false, &param_names, injected_mu, indicator);
        MyLinearityTest(w, &param_names, injected_mu, eps, indicator);
    }
    else if (std::string(argv[1]) == std::string("nuisance")) {
        filesaver.OpenFile(true, &param_names, injected_mu, indicator);
        MyToyMCStudy(w, &param_names, eps, indicator);
    }

    f->Close();

    filesaver.CloseFile();

    return 0;
}
