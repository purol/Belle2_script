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

    RooDataSet* data = (RooDataSet*)w->data("asimovData");
    //RooFitResult* fitres = model->fitTo(*data, RooFit::Minimizer("Minuit2"), RooFit::Extended(false), RooFit::SumW2Error(false), PrintLevel(-1), Save());
    RooAbsReal* nll;
    RooFitResult* fitres = MyMinimizeNLL(w, data, &nll);

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

std::vector<double> MySetParamsForToy(RooWorkspace* w, std::vector<std::string>* names, double injected_mu) {

    std::vector<double> Nevt;
    for (int i = 0; i < RarityBins; i++) Nevt.push_back(0.0);

    w->loadSnapshot("ParamValues");
    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();
    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x_val_MXs1 = w->var("obs_x_channel_MXs1");
    RooRealVar* x_val_MXs2 = w->var("obs_x_channel_MXs2");
    RooRealVar* x_val_MXs3 = w->var("obs_x_channel_MXs3");

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
        else if ((names->at(i).find("gamma") != std::string::npos) && (names->at(i).find("uncorr") != std::string::npos)) {

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
    {
        RooAbsBinning const& binning = x_val_MXs1->getBinning();
        const double oldVal = x_val_MXs1->getVal();

        for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
            double binCenter = binning.binCenter(iBin);
            double binWidth = binning.binWidth(iBin);

            *x_val_MXs1 = binCenter; // set x value

            for (unsigned int j = 0; j < scaleFactors_pdf_names.size(); j++) {
                if (std::strstr(scaleFactors_pdf_names.at(j).c_str(), "channel_MXs1") == nullptr) continue; // skip non-MXs1

                RooAbsReal* temp_func_scaleFactors = w->function(scaleFactors_pdf_names.at(j).c_str());
                RooAbsReal* temp_func_shapes = w->function(shapes_pdf_names.at(j).c_str());
                Nevt.at(iBin) = Nevt.at(iBin) + (temp_func_scaleFactors->getValV() * temp_func_shapes->getValV());
                if ((temp_func_scaleFactors->getValV() * temp_func_shapes->getValV()) < 0) {
                    printf("[ERROR] negative count!\n");
                    exit(1);
                }
            }

        }

        *x_val_MXs1 = oldVal;
    }

    {
        RooAbsBinning const& binning = x_val_MXs2->getBinning();
        const double oldVal = x_val_MXs2->getVal();

        for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
            double binCenter = binning.binCenter(iBin);
            double binWidth = binning.binWidth(iBin);

            *x_val_MXs2 = binCenter; // set x value

            for (unsigned int j = 0; j < scaleFactors_pdf_names.size(); j++) {
                if (std::strstr(scaleFactors_pdf_names.at(j).c_str(), "channel_MXs2") == nullptr) continue; // skip non-MXs2

                RooAbsReal* temp_func_scaleFactors = w->function(scaleFactors_pdf_names.at(j).c_str());
                RooAbsReal* temp_func_shapes = w->function(shapes_pdf_names.at(j).c_str());
                Nevt.at(iBin + RarityBins_MX1) = Nevt.at(iBin + RarityBins_MX1) + (temp_func_scaleFactors->getValV() * temp_func_shapes->getValV());
                if ((temp_func_scaleFactors->getValV() * temp_func_shapes->getValV()) < 0) {
                    printf("[ERROR] negative count!\n");
                    exit(1);
                }
            }

        }

        *x_val_MXs2 = oldVal;
    }

    {
        RooAbsBinning const& binning = x_val_MXs3->getBinning();
        const double oldVal = x_val_MXs3->getVal();

        for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
            double binCenter = binning.binCenter(iBin);
            double binWidth = binning.binWidth(iBin);

            *x_val_MXs3 = binCenter; // set x value

            for (unsigned int j = 0; j < scaleFactors_pdf_names.size(); j++) {
                if (std::strstr(scaleFactors_pdf_names.at(j).c_str(), "channel_MXs3") == nullptr) continue; // skip non-MXs3

                RooAbsReal* temp_func_scaleFactors = w->function(scaleFactors_pdf_names.at(j).c_str());
                RooAbsReal* temp_func_shapes = w->function(shapes_pdf_names.at(j).c_str());
                Nevt.at(iBin + RarityBins_MX1 + RarityBins_MX2) = Nevt.at(iBin + RarityBins_MX1 + RarityBins_MX2) + (temp_func_scaleFactors->getValV() * temp_func_shapes->getValV());
                if ((temp_func_scaleFactors->getValV() * temp_func_shapes->getValV()) < 0) {
                    printf("[ERROR] negative count!\n");
                    exit(1);
                }
            }

        }

        *x_val_MXs3 = oldVal;
    }

    return Nevt;

}

RooDataSet* MyGenerate(RooWorkspace* w, std::vector<double> Nevts, bool extended) {

    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    // get variables and weight
    RooRealVar* x_val_MXs1 = w->var("obs_x_channel_MXs1");
    RooRealVar* x_val_MXs2 = w->var("obs_x_channel_MXs2");
    RooRealVar* x_val_MXs3 = w->var("obs_x_channel_MXs3");
    RooCategory* channelCat = (RooCategory*)(&model->indexCat());
    RooRealVar* weight_ = new RooRealVar("weight_", "", 0.0, 1000.0);

    // define data
    RooDataSet* genData = new RooDataSet("hmaster", "hmaster", RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, *channelCat, *weight_), weight_->GetName());

    for (int j = 0; j < RarityBins_MX1; j++) {
        // channel 1
        x_val_MXs1->setVal(0.5 + j);
        x_val_MXs2->setVal(0.5);
        x_val_MXs3->setVal(0.5);
        channelCat->setLabel("channel_MXs1");

        // generate
        if (Nevts.at(j) > MyEPSILON) {
            if (extended) {
                std::poisson_distribution<int> distribution((int)floor(Nevts.at(j) + 0.5));
                int Nentry_with_fluctuation = distribution(generator);
                genData->add(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, *channelCat), Nentry_with_fluctuation);
            }
            else {
                genData->add(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, *channelCat), (int)floor(Nevts.at(j) + 0.5));
            }
        }
        else { // no event. Maybe because of partial unblind. Just set 0
            genData->add(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, *channelCat), 0);
        }

    }
    for (int j = 0; j < RarityBins_MX2; j++) {
        // channel 2
        x_val_MXs1->setVal(RarityBins_MX1 - 0.5);
        x_val_MXs2->setVal(0.5 + j);
        x_val_MXs3->setVal(0.5);
        channelCat->setLabel("channel_MXs2");

        // generate
        if (Nevts.at(RarityBins_MX1 + j) > MyEPSILON) {
            if (extended) {
                std::poisson_distribution<int> distribution((int)floor(Nevts.at(RarityBins_MX1 + j) + 0.5));
                int Nentry_with_fluctuation = distribution(generator);
                genData->add(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, *channelCat), Nentry_with_fluctuation);
            }
            else {
                genData->add(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, *channelCat), (int)floor(Nevts.at(RarityBins_MX1 + j) + 0.5));
            }
        }
        else { // no event. Maybe because of partial unblind. Just set 0
            genData->add(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, *channelCat), 0);
        }

    }
    for (int j = 0; j < RarityBins_MX3; j++) {
        // channel 3
        x_val_MXs1->setVal(RarityBins_MX1 - 0.5);
        x_val_MXs2->setVal(RarityBins_MX2 - 0.5);
        x_val_MXs3->setVal(0.5 + j);
        channelCat->setLabel("channel_MXs3");

        // generate
        if (Nevts.at(RarityBins_MX1 + RarityBins_MX2 + j) > MyEPSILON) {
            if (extended) {
                std::poisson_distribution<int> distribution((int)floor(Nevts.at(RarityBins_MX1 + RarityBins_MX2 + j) + 0.5));
                int Nentry_with_fluctuation = distribution(generator);
                genData->add(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, *channelCat), Nentry_with_fluctuation);
            }
            else {
                genData->add(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, *channelCat), (int)floor(Nevts.at(RarityBins_MX1 + RarityBins_MX2 + j) + 0.5));
            }
        }
        else { // no event. Maybe because of partial unblind. Just set 0
            genData->add(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, *channelCat), 0);
        }

    }

    delete weight_;

    return genData;
}

double SetParamsForToy(RooWorkspace* w, std::vector<std::string>* names, double injected_mu) {

    double Nevt = 0.0;

    w->loadSnapshot("ParamValues");
    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();
    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x_val_MXs1 = w->var("obs_x_channel_MXs1");
    RooRealVar* x_val_MXs2 = w->var("obs_x_channel_MXs2");
    RooRealVar* x_val_MXs3 = w->var("obs_x_channel_MXs3");

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
    {
        RooAbsBinning const& binning = x_val_MXs1->getBinning();
        const double oldVal = x_val_MXs1->getVal();

        for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
            double binCenter = binning.binCenter(iBin);
            double binWidth = binning.binWidth(iBin);

            *x_val_MXs1 = binCenter; // set x value

            for (unsigned int j = 0; j < scaleFactors_pdf_names.size(); j++) {
                if (std::strstr(scaleFactors_pdf_names.at(j).c_str(), "channel_MXs1") == nullptr) continue; // skip non-MXs1

                RooAbsReal* temp_func_scaleFactors = w->function(scaleFactors_pdf_names.at(j).c_str());
                RooAbsReal* temp_func_shapes = w->function(shapes_pdf_names.at(j).c_str());
                Nevt = Nevt + (temp_func_scaleFactors->getValV() * temp_func_shapes->getValV());
                if ((temp_func_scaleFactors->getValV() * temp_func_shapes->getValV()) < 0) {
                    printf("[ERROR] negative count!\n");
                    exit(1);
                }
            }

        }

        *x_val_MXs1 = oldVal;
    }

    {
        RooAbsBinning const& binning = x_val_MXs2->getBinning();
        const double oldVal = x_val_MXs2->getVal();

        for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
            double binCenter = binning.binCenter(iBin);
            double binWidth = binning.binWidth(iBin);

            *x_val_MXs2 = binCenter; // set x value

            for (unsigned int j = 0; j < scaleFactors_pdf_names.size(); j++) {
                if (std::strstr(scaleFactors_pdf_names.at(j).c_str(), "channel_MXs2") == nullptr) continue; // skip non-MXs2

                RooAbsReal* temp_func_scaleFactors = w->function(scaleFactors_pdf_names.at(j).c_str());
                RooAbsReal* temp_func_shapes = w->function(shapes_pdf_names.at(j).c_str());
                Nevt = Nevt + (temp_func_scaleFactors->getValV() * temp_func_shapes->getValV());
                if ((temp_func_scaleFactors->getValV() * temp_func_shapes->getValV()) < 0) {
                    printf("[ERROR] negative count!\n");
                    exit(1);
                }
            }

        }

        *x_val_MXs2 = oldVal;
    }

    {
        RooAbsBinning const& binning = x_val_MXs3->getBinning();
        const double oldVal = x_val_MXs3->getVal();

        for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
            double binCenter = binning.binCenter(iBin);
            double binWidth = binning.binWidth(iBin);

            *x_val_MXs3 = binCenter; // set x value

            for (unsigned int j = 0; j < scaleFactors_pdf_names.size(); j++) {
                if (std::strstr(scaleFactors_pdf_names.at(j).c_str(), "channel_MXs3") == nullptr) continue; // skip non-MXs3

                RooAbsReal* temp_func_scaleFactors = w->function(scaleFactors_pdf_names.at(j).c_str());
                RooAbsReal* temp_func_shapes = w->function(shapes_pdf_names.at(j).c_str());
                Nevt = Nevt + (temp_func_scaleFactors->getValV() * temp_func_shapes->getValV());
                if ((temp_func_scaleFactors->getValV() * temp_func_shapes->getValV()) < 0) {
                    printf("[ERROR] negative count!\n");
                    exit(1);
                }
            }

        }

        *x_val_MXs3 = oldVal;
    }

    return Nevt;

}

double DoNotSetParamsForToy(RooWorkspace* w, std::vector<std::string>* names, double injected_mu) {

    double Nevt = 0.0;

    w->loadSnapshot("ParamValues");
    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();
    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x_val_MXs1 = w->var("obs_x_channel_MXs1");
    RooRealVar* x_val_MXs2 = w->var("obs_x_channel_MXs2");
    RooRealVar* x_val_MXs3 = w->var("obs_x_channel_MXs3");

    w->var("mu_MXs1")->setVal(injected_mu);
    w->var("mu_MXs2")->setVal(injected_mu);
    w->var("mu_MXs3")->setVal(injected_mu);

    /* ================================ cal Nexpected ================================*/
    {
        RooAbsBinning const& binning = x_val_MXs1->getBinning();
        const double oldVal = x_val_MXs1->getVal();

        for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
            double binCenter = binning.binCenter(iBin);
            double binWidth = binning.binWidth(iBin);

            *x_val_MXs1 = binCenter; // set x value

            for (unsigned int j = 0; j < scaleFactors_pdf_names.size(); j++) {
                if (std::strstr(scaleFactors_pdf_names.at(j).c_str(), "channel_MXs1") == nullptr) continue; // skip non-MXs1

                RooAbsReal* temp_func_scaleFactors = w->function(scaleFactors_pdf_names.at(j).c_str());
                RooAbsReal* temp_func_shapes = w->function(shapes_pdf_names.at(j).c_str());
                Nevt = Nevt + (temp_func_scaleFactors->getValV() * temp_func_shapes->getValV());
                if ((temp_func_scaleFactors->getValV() * temp_func_shapes->getValV()) < 0) {
                    printf("[ERROR] negative count!\n");
                    exit(1);
                }
            }

        }

        *x_val_MXs1 = oldVal;
    }

    {
        RooAbsBinning const& binning = x_val_MXs2->getBinning();
        const double oldVal = x_val_MXs2->getVal();

        for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
            double binCenter = binning.binCenter(iBin);
            double binWidth = binning.binWidth(iBin);

            *x_val_MXs2 = binCenter; // set x value

            for (unsigned int j = 0; j < scaleFactors_pdf_names.size(); j++) {
                if (std::strstr(scaleFactors_pdf_names.at(j).c_str(), "channel_MXs2") == nullptr) continue; // skip non-MXs2

                RooAbsReal* temp_func_scaleFactors = w->function(scaleFactors_pdf_names.at(j).c_str());
                RooAbsReal* temp_func_shapes = w->function(shapes_pdf_names.at(j).c_str());
                Nevt = Nevt + (temp_func_scaleFactors->getValV() * temp_func_shapes->getValV());
                if ((temp_func_scaleFactors->getValV() * temp_func_shapes->getValV()) < 0) {
                    printf("[ERROR] negative count!\n");
                    exit(1);
                }
            }

        }

        *x_val_MXs2 = oldVal;
    }

    {
        RooAbsBinning const& binning = x_val_MXs3->getBinning();
        const double oldVal = x_val_MXs3->getVal();

        for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
            double binCenter = binning.binCenter(iBin);
            double binWidth = binning.binWidth(iBin);

            *x_val_MXs3 = binCenter; // set x value

            for (unsigned int j = 0; j < scaleFactors_pdf_names.size(); j++) {
                if (std::strstr(scaleFactors_pdf_names.at(j).c_str(), "channel_MXs3") == nullptr) continue; // skip non-MXs3

                RooAbsReal* temp_func_scaleFactors = w->function(scaleFactors_pdf_names.at(j).c_str());
                RooAbsReal* temp_func_shapes = w->function(shapes_pdf_names.at(j).c_str());
                Nevt = Nevt + (temp_func_scaleFactors->getValV() * temp_func_shapes->getValV());
                if ((temp_func_scaleFactors->getValV() * temp_func_shapes->getValV()) < 0) {
                    printf("[ERROR] negative count!\n");
                    exit(1);
                }
            }

        }

        *x_val_MXs3 = oldVal;
    }

    return Nevt;

}

std::vector<double> GetPoissonError(RooWorkspace* w) {
    std::vector<double> Nevts;
    std::vector<double> errors;

    RooRealVar* x_val_MXs1 = w->var("obs_x_channel_MXs1");
    RooRealVar* x_val_MXs2 = w->var("obs_x_channel_MXs2");
    RooRealVar* x_val_MXs3 = w->var("obs_x_channel_MXs3");

    {
        RooAbsBinning const& binning = x_val_MXs1->getBinning();
        const double oldVal = x_val_MXs1->getVal();

        for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
            double Nevt = 0.0;
            double binCenter = binning.binCenter(iBin);
            double binWidth = binning.binWidth(iBin);

            *x_val_MXs1 = binCenter; // set x value

            for (unsigned int j = 0; j < scaleFactors_pdf_names.size(); j++) {
                if (std::strstr(scaleFactors_pdf_names.at(j).c_str(), "channel_MXs1") == nullptr) continue; // skip non-MXs1

                RooAbsReal* temp_func_scaleFactors = w->function(scaleFactors_pdf_names.at(j).c_str());
                RooAbsReal* temp_func_shapes = w->function(shapes_pdf_names.at(j).c_str());
                Nevt = Nevt + (temp_func_scaleFactors->getValV() * temp_func_shapes->getValV());
            }
            Nevts.push_back(Nevt);
        }

        *x_val_MXs1 = oldVal;
    }

    {
        RooAbsBinning const& binning = x_val_MXs2->getBinning();
        const double oldVal = x_val_MXs2->getVal();

        for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
            double Nevt = 0.0;
            double binCenter = binning.binCenter(iBin);
            double binWidth = binning.binWidth(iBin);

            *x_val_MXs2 = binCenter; // set x value

            for (unsigned int j = 0; j < scaleFactors_pdf_names.size(); j++) {
                if (std::strstr(scaleFactors_pdf_names.at(j).c_str(), "channel_MXs2") == nullptr) continue; // skip non-MXs2

                RooAbsReal* temp_func_scaleFactors = w->function(scaleFactors_pdf_names.at(j).c_str());
                RooAbsReal* temp_func_shapes = w->function(shapes_pdf_names.at(j).c_str());
                Nevt = Nevt + (temp_func_scaleFactors->getValV() * temp_func_shapes->getValV());
            }
            Nevts.push_back(Nevt);
        }

        *x_val_MXs2 = oldVal;
    }

    {
        RooAbsBinning const& binning = x_val_MXs3->getBinning();
        const double oldVal = x_val_MXs3->getVal();

        for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
            double Nevt = 0.0;
            double binCenter = binning.binCenter(iBin);
            double binWidth = binning.binWidth(iBin);

            *x_val_MXs3 = binCenter; // set x value

            for (unsigned int j = 0; j < scaleFactors_pdf_names.size(); j++) {
                if (std::strstr(scaleFactors_pdf_names.at(j).c_str(), "channel_MXs3") == nullptr) continue; // skip non-MXs3

                RooAbsReal* temp_func_scaleFactors = w->function(scaleFactors_pdf_names.at(j).c_str());
                RooAbsReal* temp_func_shapes = w->function(shapes_pdf_names.at(j).c_str());
                Nevt = Nevt + (temp_func_scaleFactors->getValV() * temp_func_shapes->getValV());
            }
            Nevts.push_back(Nevt);
        }

        *x_val_MXs3 = oldVal;
    }

    // get poisson fluctuation
    for (int i = 0; i < Nevts.size(); i++) {
        std::poisson_distribution<int> distribution((int) floor(Nevts.at(i) + 0.5));
        errors.push_back((double)(distribution(generator) - (int)floor(Nevts.at(i) + 0.5)));
    }

    return errors;
}

std::vector<double> GetNevts(RooWorkspace* w) {
    std::vector<double> Nevts;

    RooRealVar* x_val_MXs1 = w->var("obs_x_channel_MXs1");
    RooRealVar* x_val_MXs2 = w->var("obs_x_channel_MXs2");
    RooRealVar* x_val_MXs3 = w->var("obs_x_channel_MXs3");

    {
        RooAbsBinning const& binning = x_val_MXs1->getBinning();
        const double oldVal = x_val_MXs1->getVal();

        for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
            double Nevt = 0.0;
            double binCenter = binning.binCenter(iBin);
            double binWidth = binning.binWidth(iBin);

            *x_val_MXs1 = binCenter; // set x value

            for (unsigned int j = 0; j < scaleFactors_pdf_names.size(); j++) {
                if (std::strstr(scaleFactors_pdf_names.at(j).c_str(), "channel_MXs1") == nullptr) continue; // skip non-MXs1

                RooAbsReal* temp_func_scaleFactors = w->function(scaleFactors_pdf_names.at(j).c_str());
                RooAbsReal* temp_func_shapes = w->function(shapes_pdf_names.at(j).c_str());
                Nevt = Nevt + (temp_func_scaleFactors->getValV() * temp_func_shapes->getValV());
            }
            Nevts.push_back(Nevt);
        }

        *x_val_MXs1 = oldVal;
    }

    {
        RooAbsBinning const& binning = x_val_MXs2->getBinning();
        const double oldVal = x_val_MXs2->getVal();

        for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
            double Nevt = 0.0;
            double binCenter = binning.binCenter(iBin);
            double binWidth = binning.binWidth(iBin);

            *x_val_MXs2 = binCenter; // set x value

            for (unsigned int j = 0; j < scaleFactors_pdf_names.size(); j++) {
                if (std::strstr(scaleFactors_pdf_names.at(j).c_str(), "channel_MXs2") == nullptr) continue; // skip non-MXs2

                RooAbsReal* temp_func_scaleFactors = w->function(scaleFactors_pdf_names.at(j).c_str());
                RooAbsReal* temp_func_shapes = w->function(shapes_pdf_names.at(j).c_str());
                Nevt = Nevt + (temp_func_scaleFactors->getValV() * temp_func_shapes->getValV());
            }
            Nevts.push_back(Nevt);
        }

        *x_val_MXs2 = oldVal;
    }

    {
        RooAbsBinning const& binning = x_val_MXs3->getBinning();
        const double oldVal = x_val_MXs3->getVal();

        for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
            double Nevt = 0.0;
            double binCenter = binning.binCenter(iBin);
            double binWidth = binning.binWidth(iBin);

            *x_val_MXs3 = binCenter; // set x value

            for (unsigned int j = 0; j < scaleFactors_pdf_names.size(); j++) {
                if (std::strstr(scaleFactors_pdf_names.at(j).c_str(), "channel_MXs3") == nullptr) continue; // skip non-MXs3

                RooAbsReal* temp_func_scaleFactors = w->function(scaleFactors_pdf_names.at(j).c_str());
                RooAbsReal* temp_func_shapes = w->function(shapes_pdf_names.at(j).c_str());
                Nevt = Nevt + (temp_func_scaleFactors->getValV() * temp_func_shapes->getValV());
            }
            Nevts.push_back(Nevt);
        }

        *x_val_MXs3 = oldVal;
    }

    return Nevts;
}

void MyToyMCStudy(RooWorkspace *w, std::vector<std::string>* names, double eps, int indicator = 0){

        ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
        RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

        RooArgSet* obs = (RooArgSet*)mc->GetObservables();
        RooRealVar* x_val_MXs1 = w->var("obs_x_channel_MXs1");
        RooRealVar* x_val_MXs2 = w->var("obs_x_channel_MXs2");
        RooRealVar* x_val_MXs3 = w->var("obs_x_channel_MXs3");

        for(int i=0; i< Toy_iter_num; i++) { // Do Toy MC study
            
            double Nevt_total = SetParamsForToy(w, names, 1.0);
            // std::vector<double> Nevt_total = MySetParamsForToy(w, names, 1.0);

            filesaver.GetTrueValues(w, names);

            RooDataSet* genData = model->generate(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, model->indexCat()), Nevt_total, false, true, "", false, true);
            // RooDataSet* genData = MyGenerate(w, Nevt_total, true);

            w->loadSnapshot("ParamValues");
            //RooFitResult* fitres = model->fitTo(*genData, RooFit::Minimizer("Minuit2"), RooFit::Extended(false), RooFit::Minos(RooArgSet(*w->var("mu_MXs1"), *w->var("mu_MXs2"), *w->var("mu_MXs3"))), RooFit::SumW2Error(false), Save());
            RooAbsReal* nll;
            RooFitResult* fitres = MyMinimizeNLL(w, genData, &nll, eps);

            if(MyDEBUG) Debug(w, fitres, genData);

            filesaver.GetFittingValues(fitres, names);
            filesaver.GetFittingStatus(fitres);
            filesaver.WriteIntoBranch();

            delete fitres;
            delete genData;

        }
}

void MyToyMCRCStudy(RooWorkspace* w, std::vector<std::string>* names, double eps, int indicator = 0) {
    // based on Makus-san's comment, get Poisson fluctuation first, then fluctuate nuisance parameters

    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x_val_MXs1 = w->var("obs_x_channel_MXs1");
    RooRealVar* x_val_MXs2 = w->var("obs_x_channel_MXs2");
    RooRealVar* x_val_MXs3 = w->var("obs_x_channel_MXs3");
    RooCategory* channelCat = (RooCategory*)(&model->indexCat());
    RooRealVar* weight_ = new RooRealVar("weight_", "", 0.0, 1000.0);

    for (int i = 0; i < Toy_iter_num; i++) { // Do Toy MC study

        // get poisson error with nominal configuration
        w->loadSnapshot("ParamValues");
        std::vector<double> PoissonError = GetPoissonError(w);

        double Nevt_total = SetParamsForToy(w, names, 1.0);
        std::vector<double> Nevts = GetNevts(w);

        filesaver.GetTrueValues(w, names);

        // generate RooDataSet
        w->loadSnapshot("ParamValues");
        RooDataSet* genData = new RooDataSet("hmaster", "hmaster", RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, *channelCat, *weight_), weight_->GetName());
        for (int j = 0; j < RarityBins_MX1; j++) {
            // channel 1
            x_val_MXs1->setVal(0.5 + j);
            x_val_MXs2->setVal(0.5);
            x_val_MXs3->setVal(0.5);
            channelCat->setLabel("channel_MXs1");
            genData->add(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, *channelCat), Nevts.at(j) + PoissonError.at(j));
        }
        for (int j = 0; j < RarityBins_MX2; j++) {
            // channel 2
            x_val_MXs1->setVal(RarityBins_MX1 - 0.5);
            x_val_MXs2->setVal(0.5 + j);
            x_val_MXs3->setVal(0.5);
            channelCat->setLabel("channel_MXs2");
            genData->add(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, *channelCat), Nevts.at(j + RarityBins_MX1) + PoissonError.at(j + RarityBins_MX1));
        }
        for (int j = 0; j < RarityBins_MX3; j++) {
            // channel 2
            x_val_MXs1->setVal(RarityBins_MX1 - 0.5);
            x_val_MXs2->setVal(RarityBins_MX2 - 0.5);
            x_val_MXs3->setVal(0.5 + j);
            channelCat->setLabel("channel_MXs3");
            genData->add(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, *channelCat), Nevts.at(j + RarityBins_MX1 + RarityBins_MX2) + PoissonError.at(j + RarityBins_MX1 + RarityBins_MX2));
        }

        w->loadSnapshot("ParamValues");
        //RooFitResult* fitres = model->fitTo(*genData, RooFit::Minimizer("Minuit2"), RooFit::Extended(false), RooFit::Minos(RooArgSet(*w->var("mu_MXs1"), *w->var("mu_MXs2"), *w->var("mu_MXs3"))), RooFit::SumW2Error(false), Save());
        RooAbsReal* nll;
        RooFitResult* fitres = MyMinimizeNLL(w, genData, &nll, eps);

        delete genData;

        if (MyDEBUG) Debug(w, fitres, genData);

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
    RooRealVar* x_val_MXs1 = w->var("obs_x_channel_MXs1");
    RooRealVar* x_val_MXs2 = w->var("obs_x_channel_MXs2");
    RooRealVar* x_val_MXs3 = w->var("obs_x_channel_MXs3");

    for (int i = 0; i < LT_iter_num; i++) { // Do LT MC study

        double Nevt_total = SetParamsForToy(w, names, mu_injected);
        // std::vector<double> Nevt_total = MySetParamsForToy(w, names, mu_injected);

        filesaver.GetTrueValues(w, names);

        RooDataSet* genData = model->generate(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, model->indexCat()), Nevt_total, false, true, "", false, true);
        // RooDataSet* genData = MyGenerate(w, Nevt_total, true);

        w->loadSnapshot("ParamValues");
        //RooFitResult* fitres = model->fitTo(*genData, RooFit::Minimizer("Minuit2"), RooFit::Extended(false), RooFit::Minos(RooArgSet(*w->var("mu_MXs1"), *w->var("mu_MXs2"), *w->var("mu_MXs3"))), RooFit::SumW2Error(false), Save());
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
    RooRealVar* x_val_MXs1 = w->var("obs_x_channel_MXs1");
    RooRealVar* x_val_MXs2 = w->var("obs_x_channel_MXs2");
    RooRealVar* x_val_MXs3 = w->var("obs_x_channel_MXs3");

    // get Category and data
    RooCategory* idx = (RooCategory*)obs->find("channelCat");
    RooDataSet* data = (RooDataSet*)w->data("asimovData");

    // fit
    //RooFitResult* fitres = model->fitTo(*data, RooFit::Minimizer("Minuit2"), RooFit::Extended(false), RooFit::Minos(RooArgSet(*w->var("mu_MXs1"), *w->var("mu_MXs2"), *w->var("mu_MXs3"))), RooFit::SumW2Error(false), Save());
    RooAbsReal* nll;
    RooFitResult* fitres = MyMinimizeNLL(w, data, &nll, eps);

    // draw
    // GetPlotTemplate(w, data);

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

    delete fitres;
}

void MyToyMCStudyDataPoisson(RooWorkspace* w, std::vector<std::string>* names, double eps, RooDataSet* data, int indicator = 0) {

    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x_val_MXs1 = w->var("obs_x_channel_MXs1");
    RooRealVar* x_val_MXs2 = w->var("obs_x_channel_MXs2");
    RooRealVar* x_val_MXs3 = w->var("obs_x_channel_MXs3");
    RooCategory* channelCat = (RooCategory*)(&model->indexCat());

    for (int i = 0; i < Toy_iter_num; i++) { // Do Toy MC study

        RooRealVar* weight_ = new RooRealVar("weight_", "", 0.0, 1000.0);

        // get nominal configuration
        w->loadSnapshot("ParamValues");
        filesaver.GetTrueValues(w, names);

        // generate RooDataSet
        w->loadSnapshot("ParamValues");
        RooDataSet* genData = new RooDataSet("hmaster", "hmaster", RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, *channelCat, *weight_), weight_->GetName());
        for (int j = 0; j < RarityBins_MX1; j++) {
            // channel 1
            x_val_MXs1->setVal(0.5 + j);
            x_val_MXs2->setVal(0.5);
            x_val_MXs3->setVal(0.5);
            channelCat->setLabel("channel_MXs1");

            // load data
            data->get(j);

            // generate data with Poisson fluctuation
            if (data->weight() > MyEPSILON) {
                std::poisson_distribution<int> distribution((int)floor(data->weight() + 0.5));
                int Nentry_with_fluctuation = distribution(generator);
                genData->add(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, *channelCat), Nentry_with_fluctuation);
            }
            else { // no event. Maybe because of partial unblind. Just set 0
                genData->add(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, *channelCat), 0);
            }

        }
        for (int j = 0; j < RarityBins_MX2; j++) {
            // channel 2
            x_val_MXs1->setVal(RarityBins_MX1 - 0.5);
            x_val_MXs2->setVal(0.5 + j);
            x_val_MXs3->setVal(0.5);
            channelCat->setLabel("channel_MXs2");

            // load data
            data->get(RarityBins_MX1 + j);

            // generate data with Poisson fluctuation
            if (data->weight() > MyEPSILON) {
                std::poisson_distribution<int> distribution((int)floor(data->weight() + 0.5));
                int Nentry_with_fluctuation = distribution(generator);
                genData->add(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, *channelCat), Nentry_with_fluctuation);
            }
            else { // no event. Maybe because of partial unblind. Just set 0
                genData->add(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, *channelCat), 0);
            }

        }
        for (int j = 0; j < RarityBins_MX3; j++) {
            // channel 3
            x_val_MXs1->setVal(RarityBins_MX1 - 0.5);
            x_val_MXs2->setVal(RarityBins_MX2 - 0.5);
            x_val_MXs3->setVal(0.5 + j);
            channelCat->setLabel("channel_MXs3");

            // load data
            data->get(RarityBins_MX1 + RarityBins_MX2 + j);

            // generate data with Poisson fluctuation
            if (data->weight() > MyEPSILON) {
                std::poisson_distribution<int> distribution((int)floor(data->weight() + 0.5));
                int Nentry_with_fluctuation = distribution(generator);
                genData->add(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, *channelCat), Nentry_with_fluctuation);
            }
            else {
                genData->add(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3, *channelCat), 0);
            }

        }

        w->loadSnapshot("ParamValues");
        //RooFitResult* fitres = model->fitTo(*genData, RooFit::Minimizer("Minuit2"), RooFit::Extended(false), RooFit::Minos(RooArgSet(*w->var("mu_MXs1"), *w->var("mu_MXs2"), *w->var("mu_MXs3"))), RooFit::SumW2Error(false), Save());
        RooAbsReal* nll;
        RooFitResult* fitres = MyMinimizeNLL(w, genData, &nll, eps);

        delete genData;
        delete weight_;

        if (MyDEBUG) Debug(w, fitres, genData);

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
                RooRealVar* x_val_MXs1 = w->var("obs_x_channel_MXs1");
                RooRealVar* x_val_MXs2 = w->var("obs_x_channel_MXs2");
                RooRealVar* x_val_MXs3 = w->var("obs_x_channel_MXs3");

                // get Category and data
                RooCategory* idx = (RooCategory*)obs->find("channelCat");

                // draw
                // GetPlotTemplate(w, data);


            }
        }

    }
}

std::vector<double> ReadFittedNevt(const char* filename) {
    std::vector<double> Nevts;
    for (int i = 0; i < RarityBins; i++) Nevts.push_back(0.0);

    FILE* fp = fopen(filename, "r");

    fscanf(fp, "Signal_MXs1:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        Nevts.at(i) = Nevts.at(i) + temp;
    }
    fscanf(fp, "\n");

    fscanf(fp, "Signal_MXs2:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        Nevts.at(i) = Nevts.at(i) + temp;
    }
    fscanf(fp, "\n");

    fscanf(fp, "Signal_MXs3:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        Nevts.at(i) = Nevts.at(i) + temp;
    }
    fscanf(fp, "\n");

    fscanf(fp, "CHG:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        Nevts.at(i) = Nevts.at(i) + temp;
    }
    fscanf(fp, "\n");

    fscanf(fp, "MIX:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        Nevts.at(i) = Nevts.at(i) + temp;
    }
    fscanf(fp, "\n");

    fscanf(fp, "UUBAR:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        Nevts.at(i) = Nevts.at(i) + temp;
    }
    fscanf(fp, "\n");

    fscanf(fp, "DDBAR:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        Nevts.at(i) = Nevts.at(i) + temp;
    }
    fscanf(fp, "\n");

    fscanf(fp, "SSBAR:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        Nevts.at(i) = Nevts.at(i) + temp;
    }
    fscanf(fp, "\n");

    fscanf(fp, "CHARM:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        Nevts.at(i) = Nevts.at(i) + temp;
    }
    fscanf(fp, "\n");

    fclose(fp);

    return Nevts;
}

void MyToyMCStudyWithNevts(RooWorkspace* w, std::vector<std::string>* names, double eps, std::vector<double> Nevts, int indicator = 0) {

    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x_val_MXs1 = w->var("obs_x_channel_MXs1");
    RooRealVar* x_val_MXs2 = w->var("obs_x_channel_MXs2");
    RooRealVar* x_val_MXs3 = w->var("obs_x_channel_MXs3");

    for (int i = 0; i < Toy_iter_num; i++) { // Do Toy MC study

        filesaver.GetTrueValues(w, names);

        RooDataSet* genData = MyGenerate(w, Nevts, true);

        w->loadSnapshot("ParamValues");
        //RooFitResult* fitres = model->fitTo(*genData, RooFit::Minimizer("Minuit2"), RooFit::Extended(false), RooFit::Minos(RooArgSet(*w->var("mu_MXs1"), *w->var("mu_MXs2"), *w->var("mu_MXs3"))), RooFit::SumW2Error(false), Save());
        RooAbsReal* nll;
        RooFitResult* fitres = MyMinimizeNLL(w, genData, &nll, eps);

        if (MyDEBUG) Debug(w, fitres, genData);

        filesaver.GetFittingValues(fitres, names);
        filesaver.GetFittingStatus(fitres);
        filesaver.WriteIntoBranch();

        delete fitres;
        delete genData;

    }
}

int main(int argc, char* argv[]) {

    ::ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(2);
    ::ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Minimize"); // default: Minuit Migrad
    ::ROOT::Math::MinimizerOptions::SetDefaultStrategy(1); // default 1

    RooStats::UseNLLOffset(true); // default off

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
    
    if ((std::string(argv[1]) == std::string("ToyMC")) || (std::string(argv[1]) == std::string("ToyMCRC")) || (std::string(argv[1]) == std::string("ToyMCbox"))) {  // main ToyMC|ToyMCRC|ToyMCbox
        if (argc == 5) {
            injected_mu = -1;
            eps = std::atof(argv[2]);
            indicator = std::atoi(argv[3]);
            Toy_iter_num = std::atoi(argv[4]);
        }
        else {
            printf("Toy MC(RC) requires 4 arguments {eps} {indicator} {Num of sample}\n");
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
    else if ((std::string(argv[1]) == std::string("nuisance")) || (std::string(argv[1]) == std::string("nuisanceLOADNevt")) || (std::string(argv[1]) == std::string("nuisanceLOADNevtLOADNuisance"))) { // check the power of nuisance parameters
        if (argc == 6) {
            injected_mu = -1;
            eps = std::atof(argv[2]);
            indicator = std::atoi(argv[3]);
            Toy_iter_num = std::atoi(argv[4]);
            fixed_param = std::string(argv[5]);
        }
        else {
            printf("nuisance(LOAD) Toy MC requires 5 arguments {eps} {indicator} {Num of sample} {fixed param}\n");
            exit(1);
        }
    }
    else {
        printf("first argument should be {ToyMC|LinearityTest|nuisance|nuisanceLOADNevt|nuisanceLOADNevtLOADNuisance|ToyMCRC|ToyMCbox}\n");
        exit(1);
    }

    ::ROOT::Math::MinimizerOptions::SetDefaultTolerance(eps); // default 0.01. but it is better to use 0.001

    const char* fname = "./PDFandDATA_workspace.root";

    TFile* f = TFile::Open(fname);

    RooWorkspace* w = (RooWorkspace*)f->Get("combined");

    if ((std::string(argv[1]) == std::string("ToyMC")) || (std::string(argv[1]) == std::string("ToyMCRC")) || (std::string(argv[1]) == std::string("ToyMCbox"))) {
        ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
        RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

        w->loadSnapshot("NominalParamValues");

        RooRealVar* x_val_MXs1 = w->var("obs_x_channel_MXs1");
        RooRealVar* x_val_MXs2 = w->var("obs_x_channel_MXs2");
        RooRealVar* x_val_MXs3 = w->var("obs_x_channel_MXs3");
        std::unique_ptr<RooArgSet> params{model->getParameters(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3))};

        // save snapshot
        w->saveSnapshot("ParamValues", *params, true);
    }
    else if (std::string(argv[1]) == std::string("LinearityTest")) {
        ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
        RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

        w->loadSnapshot("NominalParamValues");

        RooRealVar* x_val_MXs1 = w->var("obs_x_channel_MXs1");
        RooRealVar* x_val_MXs2 = w->var("obs_x_channel_MXs2");
        RooRealVar* x_val_MXs3 = w->var("obs_x_channel_MXs3");
        std::unique_ptr<RooArgSet> params{model->getParameters(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3))};

        // save snapshot
        w->saveSnapshot("ParamValues", *params, true);
    }
    else if ((std::string(argv[1]) == std::string("nuisance")) || (std::string(argv[1]) == std::string("nuisanceLOADNevt"))) {
        OPTIONS* options = (OPTIONS*)malloc(sizeof(OPTIONS));
        Initialize_options(options, fixed_param.c_str());

        w->loadSnapshot("NominalParamValues");
        FixParameters(w, options);
    }
    else if (std::string(argv[1]) == std::string("nuisanceLOADNevtLOADNuisance")) {
        OPTIONS* options = (OPTIONS*)malloc(sizeof(OPTIONS));
        Initialize_options(options, fixed_param.c_str());

        w->loadSnapshot("NominalParamValues");
        LoadNuisanceParameter(w, "nuisance.txt");
        FixParameters(w, options);
    }

    GetNameOfParams(w, &param_names);

    if (std::string(argv[1]) == std::string("ToyMC")) {
        filesaver.OpenFile(true, &param_names, injected_mu, indicator);
        MyToyMCStudy(w, &param_names, eps, indicator);
    }
    else if (std::string(argv[1]) == std::string("ToyMCRC")) {
        filesaver.OpenFile(true, &param_names, injected_mu, indicator);
        MyToyMCRCStudy(w, &param_names, eps, indicator);
    }
    else if (std::string(argv[1]) == std::string("ToyMCbox")) {
        RooDataSet* data = (RooDataSet*)w->data("obsData");
        filesaver.OpenFile(true, &param_names, injected_mu, indicator);
        MyToyMCStudyDataPoisson(w, &param_names, eps, data, indicator);
    }
    else if (std::string(argv[1]) == std::string("LinearityTest")) {
        filesaver.OpenFile(false, &param_names, injected_mu, indicator);
        MyLinearityTest(w, &param_names, injected_mu, eps, indicator);
    }
    else if (std::string(argv[1]) == std::string("nuisance")) {
        filesaver.OpenFile(true, &param_names, injected_mu, indicator);
        MyToyMCStudy(w, &param_names, eps, indicator);
    }
    else if (std::string(argv[1]) == std::string("nuisanceLOADNevt")) {
        filesaver.OpenFile(true, &param_names, injected_mu, indicator);
        std::vector<double> Nevts = ReadFittedNevt("Nevts.txt");
        MyToyMCStudyWithNevts(w, &param_names, eps, Nevts, indicator);
    }
    else if (std::string(argv[1]) == std::string("nuisanceLOADNevtLOADNuisance")) {
        filesaver.OpenFile(true, &param_names, injected_mu, indicator);
        std::vector<double> Nevts = ReadFittedNevt("Nevts.txt");
        MyToyMCStudyWithNevts(w, &param_names, eps, Nevts, indicator);
    }

    f->Close();

    filesaver.CloseFile();

    return 0;
}
