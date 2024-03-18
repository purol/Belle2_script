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

using namespace RooFit;
using namespace RooStats;
using namespace HistFactory;

using std::string;
using std::to_string;
using std::cout;
using std::endl;

# define MyDEBUG false
void Debug(RooWorkspace* w, RooFitResult* fitres, RooDataSet* data);

#include "template.h"
#include "correctors.h"

int Toy_iter_num = 0.0;
int LT_iter_num = 0.0;

double weight_KIDsys[RarityBins * 7] = { 0.0 };
double weight_PIDsys[RarityBins * 7] = { 0.0 };
double weight_BRsys[RarityBins * 3] = { 0.0 };
double weight_pi0sys[RarityBins * 7] = { 0.0 };
double weight_FEIsys[RarityBins * 3] = { 0.0 };

Corrector_Fragmentation corrector_Fragmentation;

std::random_device rd;
std::default_random_engine generator(rd());

std::vector<std::string> Sample_names = {
    "L_x_Signal_nominal_MXs1_channel_overallSyst_x_StatUncert_x_channel_Signal_all_uncorr_uncer_ShapeSys",
    "L_x_Signal_nominal_MXs2_channel_overallSyst_x_StatUncert_x_channel_Signal_all_uncorr_uncer_ShapeSys",
    "L_x_Signal_nominal_MXs3_channel_overallSyst_x_StatUncert_x_channel_Signal_all_uncorr_uncer_ShapeSys",
    "L_x_CHG_nominal_channel_overallSyst_x_StatUncert_x_channel_CHG_all_uncorr_uncer_ShapeSys",
    "L_x_MIX_nominal_channel_overallSyst_x_StatUncert_x_channel_MIX_all_uncorr_uncer_ShapeSys",
    "L_x_UUBAR_nominal_channel_overallSyst_x_StatUncert_x_channel_UUBAR_all_uncorr_uncer_ShapeSys",
    "L_x_DDBAR_nominal_channel_overallSyst_x_StatUncert_x_channel_DDBAR_all_uncorr_uncer_ShapeSys",
    "L_x_SSBAR_nominal_channel_overallSyst_x_StatUncert_x_channel_SSBAR_all_uncorr_uncer_ShapeSys",
    "L_x_CHARM_nominal_channel_overallSyst_x_StatUncert_x_channel_CHARM_all_uncorr_uncer_ShapeSys"
};

typedef struct Options
{
    // turn on it?
    bool track;
    bool PID;
    bool KID;
    bool KS0;
    bool pi0;
    bool FBDT;
    bool FEI;
    bool qqbar;
    bool multiplicity;
    bool Kff;
    bool Kstarff;
    bool OLD_Kff;
    bool pf;
    bool Transition;
    bool mb;
    bool fraction;
    bool MCstat;
    bool Fragmentation;
    bool mKstar;
    bool BKGNorm;
    bool BDTc;
    bool BBcounting;
    bool BBBR;
    bool BRBtoKKLKL;
    bool BRKnn;
    bool BRDKL0;
    bool uncorrelated;

    int NEntryFEI;
    int NEntryKID;
    int NEntryPID;
    int NEntryBR;
    int NEntrypi0;
    int NEntryMultiplicity;

    int Nsyst = 25;
} OPTIONS;

void GetNameOfParams(RooWorkspace* w, std::vector<std::string>* names) {
    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    w->loadSnapshot("ParamValues");

    RooAbsData* data = (RooAbsData*)w->data("asimovData");
    RooFitResult* fitres = model->fitTo(*data, RooFit::Extended(true), RooFit::SumW2Error(false), PrintLevel(-1), Save());

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

RooFitResult* MyMinimizeNLL(RooWorkspace* w, RooDataSet* data, RooAbsReal* nll, double tolerance = -1.0) {
    // what we have done
    w->loadSnapshot("ParamValues");
    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    // get nll
    RooArgSet* allParams = model->getParameters(*data);
    RooStats::RemoveConstantParameters(allParams);
    RooArgSet fGlobalObs = *mc->GetGlobalObservables();
    RooArgSet fConditionalObs;
    Bool_t fLOffset = RooStats::IsNLLOffset();
    nll = model->createNLL(*data, RooFit::CloneData(kFALSE), RooFit::Constrain(*allParams), RooFit::GlobalObservables(fGlobalObs), RooFit::ConditionalObservables(fConditionalObs), RooFit::Offset(fLOffset));

    // minimizer option
    TString fMinimizer = ::ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str();
    TString minimizer = fMinimizer;

    TString algorithm = ::ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo();

    Int_t fStrategy = ::ROOT::Math::MinimizerOptions::DefaultStrategy();

    Double_t fTolerance;
    if (tolerance < 0) fTolerance = TMath::Max(1., ::ROOT::Math::MinimizerOptions::DefaultTolerance());
    else fTolerance = tolerance;

    Int_t fPrintLevel = ::ROOT::Math::MinimizerOptions::DefaultPrintLevel();
    //LM: RooMinimizer.setPrintLevel has +1 offset - so subtract  here -1 + an extra -1
    int level = (fPrintLevel == 0) ? -1 : fPrintLevel - 2;


    // follow what ProfileLikelihoodTestStat.cxx does
    const auto& config = RooStats::GetGlobalRooStatsConfig();
    RooMinimizer minim(*nll);
    minim.setStrategy(fStrategy);
    minim.setEvalErrorWall(config.useEvalErrorWall);
    minim.setEps(fTolerance);
    minim.setPrintLevel(level);
    // this causes a memory leak
    minim.optimizeConst(2);
    minim.migrad();
    minim.minos(RooArgSet(*w->var("mu_MXs1")));
    minim.minos(RooArgSet(*w->var("mu_MXs2")));
    minim.minos(RooArgSet(*w->var("mu_MXs3")));

    // fit!
    int status;
    status = minim.minimize(minimizer, algorithm);

    return minim.save();
}

RooFitResult* MinimizeNLL(RooWorkspace* w, RooDataSet* data, RooAbsReal* nll, double tolerance = -1.0) { // this function follows the procedure in ProfileLikelihoodTestStat.cxx
    // what we have done
    // deprecated
    w->loadSnapshot("ParamValues");
    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    // get nll
    RooArgSet* allParams = model->getParameters(*data);
    RooStats::RemoveConstantParameters(allParams);
    RooArgSet fGlobalObs;
    RooArgSet fConditionalObs;
    Bool_t fLOffset = RooStats::IsNLLOffset();
    nll = model->createNLL(*data, RooFit::CloneData(kFALSE), RooFit::Constrain(*allParams), RooFit::GlobalObservables(fGlobalObs), RooFit::ConditionalObservables(fConditionalObs), RooFit::Offset(fLOffset));

    // get default values. follow ProfileLikelihoodTestStat.h
    TString fMinimizer = ::ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str();
    Int_t fStrategy = ::ROOT::Math::MinimizerOptions::DefaultStrategy();
    // avoid default tolerance to be too small (1. is default in RooMinimizer)
    Double_t fTolerance;
    if (tolerance < 0) fTolerance = TMath::Max(1., ::ROOT::Math::MinimizerOptions::DefaultTolerance());
    else fTolerance = tolerance;
    Int_t fPrintLevel = ::ROOT::Math::MinimizerOptions::DefaultPrintLevel();

    // follow what ProfileLikelihoodTestStat.cxx does
    const auto& config = RooStats::GetGlobalRooStatsConfig();
    RooMinimizer minim(*nll);
    minim.setStrategy(fStrategy);
    minim.setEvalErrorWall(config.useEvalErrorWall);
    //LM: RooMinimizer.setPrintLevel has +1 offset - so subtract  here -1 + an extra -1
    int level = (fPrintLevel == 0) ? -1 : fPrintLevel - 2;
    minim.setPrintLevel(level);
    minim.setEps(fTolerance);
    // this causes a memory leak
    minim.optimizeConst(2);
    TString minimizer = fMinimizer;
    TString algorithm = ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo();
    if (algorithm == "Migrad") algorithm = "Minimize"; // prefer to use Minimize instead of Migrad
    int status;
    for (int tries = 1, maxtries = 4; tries <= maxtries; ++tries) {
        status = minim.minimize(minimizer, algorithm);
        if (status % 1000 == 0) {  // ignore erros from Improve
            break;
        }
        else if (tries < maxtries) {
            cout << "    ----> Doing a re-scan first" << endl;
            minim.minimize(minimizer, "Scan");
            if (tries == 2) {
                if (fStrategy == 0) {
                    cout << "    ----> trying with strategy = 1" << endl;;
                    minim.setStrategy(1);
                }
                else
                    tries++; // skip this trial if strategy is already 1
            }
            if (tries == 3) {
                cout << "    ----> trying with improve" << endl;;
                minimizer = "Minuit";
                algorithm = "migradimproved";
            }
        }
    }
    /*
    if (MyDEBUG) {

        RooArgSet fitargs = minim.save()->floatParsFinal();
        TIterator* iter(fitargs.createIterator());

        for (TObject* a = iter->Next(); a != 0; a = iter->Next()) {
            RooRealVar* rrv = dynamic_cast<RooRealVar*>(a);
            std::string name = rrv->GetName();
            double val = rrv->getVal();
            double err = rrv->getError();

            if (name == "mu") {
                if (err < 1.5) {
                    RooRealVar* poi = (RooRealVar*)mc->GetParametersOfInterest()->first();
                    RooPlot* x_frame = poi->frame(Title("profile LL"));
                    //nll->plotOn(x_frame);
                    RooAbsReal* pll = nll->createProfile(*poi);
                    pll->plotOn(x_frame);

                    TCanvas* canvas = new TCanvas("sPlot", "sPlot demo", 700, 700);

                    // draw PDFs
                    x_frame->Draw();

                    canvas->SaveAs((std::string("pll_plot_") + "mu-" + std::to_string(val) + "_" + "err-" + std::to_string(err) + ".png").c_str());

                    delete canvas;
                }

            }
        }
    }
    */
    return minim.save();
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
            RooRealVar* norm = w->var(("nom_" + names->at(i)).c_str());
            std::poisson_distribution<> distribution(norm->getValV());
            w->var(names->at(i).c_str())->setVal(distribution(generator) / norm->getValV());
        }
        else if ( (names->at(i).find("gamma") != std::string::npos) && (names->at(i).find("uncorr") != std::string::npos) ) {
            int sample_index = -1;
            int bin_index = -1;

            if (names->at(i).find("CHG") != std::string::npos) sample_index = 0;
            else if(names->at(i).find("MIX") != std::string::npos) sample_index = 1;
            else if (names->at(i).find("UUBAR") != std::string::npos) sample_index = 2;
            else if (names->at(i).find("DDBAR") != std::string::npos) sample_index = 3;
            else if (names->at(i).find("SSBAR") != std::string::npos) sample_index = 4;
            else if (names->at(i).find("CHARM") != std::string::npos) sample_index = 5;
            else if (names->at(i).find("Signal") != std::string::npos) sample_index = 6;

            std::vector<std::string> temp_strings = split(names->at(i), '_');
            bin_index = stoi(temp_strings.back()); // from 0

            if (names->at(i).find("all") != std::string::npos) {
                double KID_uncertainty = weight_KIDsys[RarityBins * sample_index + bin_index];
                double PID_uncertainty = weight_PIDsys[RarityBins * sample_index + bin_index];
                double BR_uncertainty = 0.0;
                double pi0_uncertainty = weight_pi0sys[RarityBins * sample_index + bin_index];
                double FEI_uncertainty = 0.0;
                if ((names->at(i).find("CHG") != std::string::npos) || (names->at(i).find("MIX") != std::string::npos)) {
                    BR_uncertainty = weight_BRsys[RarityBins * sample_index + bin_index];
                    FEI_uncertainty = weight_FEIsys[RarityBins * sample_index + bin_index];
                }
                else if (names->at(i).find("Signal") != std::string::npos) {
                    BR_uncertainty = weight_BRsys[RarityBins * 2 + bin_index]; // exception for signal BB BR uncorrelated uncertainty!
                    FEI_uncertainty = weight_FEIsys[RarityBins * 2 + bin_index]; // exception for signal FEI uncorrelated uncertainty!
                }
                
                double total_uncertainty = std::sqrt(KID_uncertainty * KID_uncertainty + PID_uncertainty * PID_uncertainty + BR_uncertainty * BR_uncertainty + pi0_uncertainty * pi0_uncertainty + FEI_uncertainty * FEI_uncertainty);
                std::normal_distribution<double> distribution(1.0, total_uncertainty);
                w->var(names->at(i).c_str())->setVal(distribution(generator));
            }

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
            //RooFitResult* fitres = model->fitTo(*genData, RooFit::Extended(false), RooFit::Minos(RooArgSet(*w->var("mu"))), RooFit::SumW2Error(false), PrintLevel(-1), Save());
            RooAbsReal* nll;
            RooFitResult* fitres = MyMinimizeNLL(w, genData, nll, eps);

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

        //RooFitResult* fitres = model->fitTo(*genData, RooFit::Extended(false), RooFit::Minos(RooArgSet(*w->var("mu"))), RooFit::SumW2Error(false), PrintLevel(-1), Save());
        RooAbsReal* nll;
        RooFitResult* fitres = MyMinimizeNLL(w, genData, nll, eps);

        filesaver.GetFittingValues(fitres, names);
        filesaver.GetFittingStatus(fitres);
        filesaver.WriteIntoBranch();

        delete fitres;

    }
}

void ReadPIDuncorrsysFile(const char* dirname_KID, const char* dirname_PID) {
    FILE* fp;

    fp = fopen(dirname_KID, "r");
    for (int i = 0; i < RarityBins * 7; i++) fscanf(fp, "%lf\n", &weight_KIDsys[i]);
    fclose(fp);
    for (int i = 0; i < RarityBins * 7; i++) weight_KIDsys[i] = std::sqrt(weight_KIDsys[i]);

    fp = fopen(dirname_PID, "r");
    for (int i = 0; i < RarityBins * 7; i++) fscanf(fp, "%lf\n", &weight_PIDsys[i]);
    fclose(fp);
    for (int i = 0; i < RarityBins * 7; i++) weight_PIDsys[i] = std::sqrt(weight_PIDsys[i]);
}

void ReadBRuncorrsysFile(const char* dirname_BR) {
    FILE* fp;

    fp = fopen(dirname_BR, "r");
    for (int i = 0; i < RarityBins * 3; i++) fscanf(fp, "%lf\n", &weight_BRsys[i]);
    fclose(fp);
    for (int i = 0; i < RarityBins * 3; i++) weight_BRsys[i] = std::sqrt(weight_BRsys[i]);

}

void Readpi0uncorrsysFile(const char* dirname_pi0) {
    FILE* fp;

    fp = fopen(dirname_pi0, "r");
    for (int i = 0; i < RarityBins * 7; i++) fscanf(fp, "%lf\n", &weight_pi0sys[i]);
    fclose(fp);
    for (int i = 0; i < RarityBins * 7; i++) weight_pi0sys[i] = std::sqrt(weight_pi0sys[i]);

}

void ReadFEIuncorrsysFile(const char* dirname_FEI) {
    FILE* fp;

    fp = fopen(dirname_FEI, "r");
    for (int i = 0; i < RarityBins * 3; i++) fscanf(fp, "%lf\n", &weight_FEIsys[i]);
    fclose(fp);
    for (int i = 0; i < RarityBins * 3; i++) weight_FEIsys[i] = std::sqrt(weight_FEIsys[i]);

}

double GetNumEvts(RooWorkspace* w, const char* sample_type) {

    double Nevt = 0.0;

    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();
    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x_val = w->var("obs_x_channel");

    int index = -1;
    if (strcmp(sample_type, "Signal_MX1") == 0) index = 0;
    else if (strcmp(sample_type, "Signal_MX2") == 0) index = 1;
    else if (strcmp(sample_type, "Signal_MX3") == 0) index = 2;
    else if (strcmp(sample_type, "CHG") == 0) index = 3;
    else if (strcmp(sample_type, "MIX") == 0) index = 4;
    else if (strcmp(sample_type, "UUBAR") == 0) index = 5;
    else if (strcmp(sample_type, "DDBAR") == 0) index = 6;
    else if (strcmp(sample_type, "SSBAR") == 0) index = 7;
    else if (strcmp(sample_type, "CHARM") == 0) index = 8;
    else {
        printf("[ERROR] unexpected sample type!\n");
        exit(1);
    }

    /* ================================ cal Nexpected ================================*/
    RooAbsBinning const& binning = x_val->getBinning();
    const double oldVal = x_val->getVal();

    for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
        double binCenter = binning.binCenter(iBin);
        double binWidth = binning.binWidth(iBin);

        *x_val = binCenter; // set x value

        RooAbsReal* temp_func = w->function(Sample_names.at(index).c_str());
        Nevt = Nevt + temp_func->getValV();
        if (temp_func->getValV() < 0) {
            printf("[ERROR] negative count!\n");
            exit(1);
        }

    }

    *x_val = oldVal;

    return Nevt;

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
    //RooFitResult* fitres = model->fitTo(*data, RooFit::Extended(false), RooFit::Minos(RooArgSet(*w->var("mu"))), RooFit::SumW2Error(false), PrintLevel(-1), Save());
    RooAbsReal* nll;
    RooFitResult* fitres = MyMinimizeNLL(w, data, nll, eps);

    // get expected num of evts for PDFs
    double Signal_Nevts = GetNumEvts(w, "Signal");
    double CHG_Nevts = GetNumEvts(w, "CHG");
    double MIX_Nevts = GetNumEvts(w, "MIX");
    double UUBAR_Nevts = GetNumEvts(w, "UUBAR");
    double DDBAR_Nevts = GetNumEvts(w, "DDBAR");
    double SSBAR_Nevts = GetNumEvts(w, "SSBAR");
    double CHARM_Nevts = GetNumEvts(w, "CHARM");

    // draw
    RooPlot* x_frame = x_val->frame(Title("FBDT"));
    data->plotOn(x_frame, DataError(RooAbsData::Poisson), Cut("channelCat==0"), DrawOption("ZP"), Name("data"));
    //data->plotOn(x_frame, DataError(RooAbsData::Poisson), Cut("channelCat==0"), MarkerSize(0.4), DrawOption("ZP"), Normalization(1, RooAbsReal::ScaleType::NumEvent));
    model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kRed - 6), LineWidth(0), Components("L_x_*Signal*_ShapeSys,L_x_*CHG*_ShapeSys,L_x_*MIX*_ShapeSys,L_x_*UUBAR*_ShapeSys,L_x_*DDBAR*_ShapeSys,L_x_*SSBAR*_ShapeSys,L_x_*CHARM*_ShapeSys"), Normalization(Signal_Nevts + CHG_Nevts + MIX_Nevts + UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("signal"));
    model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kBlue - 6), LineWidth(0), Components("L_x_*CHG*_ShapeSys,L_x_*MIX*_ShapeSys,L_x_*UUBAR*_ShapeSys,L_x_*DDBAR*_ShapeSys,L_x_*SSBAR*_ShapeSys,L_x_*CHARM*_ShapeSys"), Normalization(CHG_Nevts + MIX_Nevts + UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("Charged B"));
    model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kCyan - 6), LineWidth(0), Components("L_x_*MIX*_ShapeSys,L_x_*UUBAR*_ShapeSys,L_x_*DDBAR*_ShapeSys,L_x_*SSBAR*_ShapeSys,L_x_*CHARM*_ShapeSys"), Normalization(MIX_Nevts + UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("Neutral B"));
    model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kOrange - 6), LineWidth(0), Components("L_x_*UUBAR*_ShapeSys,L_x_*DDBAR*_ShapeSys,L_x_*SSBAR*_ShapeSys,L_x_*CHARM*_ShapeSys"), Normalization(UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("Continuum"));
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
        //RooFitResult* fitres = model->fitTo(*genData, RooFit::Extended(false), RooFit::Minos(RooArgSet(*w->var("mu"))), RooFit::SumW2Error(false), PrintLevel(-1), Save());
        RooAbsReal* nll;
        RooFitResult* fitres = MyMinimizeNLL(w, genData, nll, eps);

        filesaver.GetFittingValues(fitres, names);
        filesaver.GetFittingStatus(fitres);
        filesaver.WriteIntoBranch();

        delete fitres;

    }
    
}

int ReadNFEIEigenVector(const char* dirname) {
    int Nentry = 0; // number of eigen values/vectors
    double eigen_value = 0; // eigen value
    double weight_sys[RarityBins * 3] = { 0.0 }; // eigen vector

    FILE* fp;
    fp = fopen(dirname, "r");
    while (true) {
        if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
        for (int i = 0; i < RarityBins * 3; i++) {
            if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
        }
        Nentry++;
    }
    fclose(fp);

    return Nentry;
}

int ReadNPIDEigenVector(const char* dirname) {
    int Nentry = 0; // number of eigen values/vectors
    double eigen_value = 0; // eigen value
    double weight_sys[RarityBins * 7] = { 0.0 }; // eigen vector

    FILE* fp;
    fp = fopen(dirname, "r");
    while (true) {
        if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
        for (int i = 0; i < RarityBins * 7; i++) {
            if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
        }
        Nentry++;
    }
    fclose(fp);

    return Nentry;
}

int ReadNBREigenVector(const char* dirname) {
    int Nentry = 0; // number of eigen values/vectors
    double eigen_value = 0; // eigen value
    double weight_sys[RarityBins * 3] = { 0.0 }; // eigen vector

    FILE* fp;
    fp = fopen(dirname, "r");
    while (true) {
        if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
        for (int i = 0; i < RarityBins * 3; i++) {
            if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
        }
        Nentry++;
    }
    fclose(fp);

    return Nentry;
}

int ReadNpi0EigenVector(const char* dirname) {
    int Nentry = 0; // number of eigen values/vectors
    double eigen_value = 0; // eigen value
    double weight_sys[RarityBins * 7] = { 0.0 }; // eigen vector

    FILE* fp;
    fp = fopen(dirname, "r");
    while (true) {
        if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
        for (int i = 0; i < RarityBins * 7; i++) {
            if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
        }
        Nentry++;
    }
    fclose(fp);

    return Nentry;
}

int ReadMultiplicityInfo(const char* dirname) {

    int N_hist = -1;

    FILE* fp;
    fp = fopen(dirname, "r");
    fscanf(fp, "%d", &N_hist);
    fclose(fp);

    return N_hist;
}

void Initialize_options(OPTIONS* options_, const char* tested_param) {
    // initialize
    options_->track = false;
    options_->PID = false;
    options_->KID = false;
    options_->KS0 = false;
    options_->pi0 = false;
    options_->FBDT = false;
    options_->FEI = false;
    options_->qqbar = false;
    options_->multiplicity = false;
    options_->Kff = false;
    options_->Kstarff = false;
    options_->OLD_Kff = false;
    options_->pf = false;
    options_->Transition = false;
    options_->mb = false;
    options_->fraction = false;
    options_->MCstat = false;
    options_->Fragmentation = false;
    options_->mKstar = false;
    options_->BKGNorm = false;
    options_->BDTc = false;
    options_->BBcounting = false;
    options_->BBBR = false;
    options_->BRBtoKKLKL = false;
    options_->BRKnn = false;
    options_->BRDKL0 = false;
    options_->uncorrelated = false;

    if (std::string(tested_param) == std::string("all")) {
        options_->track = true;
        options_->PID = true;
        options_->KID = true;
        options_->KS0 = true;
        options_->pi0 = true;
        options_->FBDT = true;
        options_->FEI = true;
        options_->qqbar = true;
        options_->multiplicity = true;
        options_->Kff = true;
        options_->Kstarff = true;
        options_->OLD_Kff = true;
        options_->pf = true;
        options_->Transition = true;
        options_->mb = true;
        options_->fraction = true;
        options_->MCstat = true;
        options_->Fragmentation = true;
        options_->mKstar = true;
        options_->BKGNorm = true;
        options_->BDTc = true;
        options_->BBcounting = true;
        options_->BBBR = true;
        options_->BRBtoKKLKL = true;
        options_->BRKnn = true;
        options_->BRDKL0 = true;
        options_->uncorrelated = true;
    }
    else if (std::string(tested_param) == std::string("none")) {}
    else if (std::string(tested_param) == std::string("track")) options_->track = true;
    else if (std::string(tested_param) == std::string("PID")) options_->PID = true;
    else if (std::string(tested_param) == std::string("KID")) options_->KID = true;
    else if (std::string(tested_param) == std::string("KS0")) options_->KS0 = true;
    else if (std::string(tested_param) == std::string("pi0")) options_->pi0 = true;
    else if (std::string(tested_param) == std::string("FBDT_efficiency")) options_->FBDT = true;
    else if (std::string(tested_param) == std::string("FEI")) options_->FEI = true;
    else if (std::string(tested_param) == std::string("qqbar_CAL")) options_->qqbar = true;
    else if (std::string(tested_param) == std::string("photon_multiplicity")) options_->multiplicity = true;
    else if (std::string(tested_param) == std::string("Kff")) options_->Kff = true;
    else if (std::string(tested_param) == std::string("Kstarff")) options_->Kstarff = true;
    else if (std::string(tested_param) == std::string("OLDKff")) options_->OLD_Kff = true;
    else if (std::string(tested_param) == std::string("pf")) options_->pf = true;
    else if (std::string(tested_param) == std::string("Transition")) options_->Transition = true;
    else if (std::string(tested_param) == std::string("mb")) options_->mb = true;
    else if (std::string(tested_param) == std::string("fraction")) options_->fraction = true;
    else if (std::string(tested_param) == std::string("MCstat")) options_->MCstat = true;
    else if (std::string(tested_param) == std::string("Fragmentation")) options_->Fragmentation = true;
    else if (std::string(tested_param) == std::string("mKstar")) options_->mKstar = true;
    else if (std::string(tested_param) == std::string("BKGNorm")) options_->BKGNorm = true;
    else if (std::string(tested_param) == std::string("BDTc")) options_->BDTc = true;
    else if (std::string(tested_param) == std::string("BBcounting")) options_->BBcounting = true;
    else if (std::string(tested_param) == std::string("BBBR")) options_->BBBR = true;
    else if (std::string(tested_param) == std::string("BRBtoKKLKL")) options_->BRBtoKKLKL = true;
    else if (std::string(tested_param) == std::string("BRKnn")) options_->BRKnn = true;
    else if (std::string(tested_param) == std::string("BtoDtoXKL")) options_->BRDKL0 = true;
    else if (std::string(tested_param) == std::string("uncorrelated")) options_->uncorrelated = true;
    else {
        printf("inappropriate parameter name\n");
        exit(1);
    }

    // read entry for nuisance parameters
    options_->NEntryFEI = ReadNFEIEigenVector("./FEI_selected.txt");
    options_->NEntryKID = ReadNPIDEigenVector("./KID_selected.txt");
    options_->NEntryPID = ReadNPIDEigenVector("./PID_selected.txt");
    options_->NEntryBR = ReadNBREigenVector("./BR_selected.txt");
    options_->NEntrypi0 = ReadNpi0EigenVector("./pi0_selected.txt");
    options_->NEntryMultiplicity = ReadMultiplicityInfo("./Multiplicity_info.txt");
}

void FixParameters(RooWorkspace* w, OPTIONS* options_) {

    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    w->loadSnapshot("NominalParamValues");

    RooRealVar* x_val = w->var("obs_x_channel");
    std::unique_ptr<RooArgSet> params{model->getParameters(*x_val)};

    // track
    w->var("alpha_track_eff_uncer")->setConstant(options_->track);

    // PID
    for (int i = 0; i < options_->NEntryPID; i++) w->var(("alpha_PID" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->PID);

    // KID
    for (int i = 0; i < options_->NEntryKID; i++) w->var(("alpha_KID" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->KID);

    // KS0
    w->var("alpha_KS0_reco_uncer")->setConstant(options_->KS0);

    // pi0
    for (int i = 0; i < options_->NEntrypi0; i++) w->var(("alpha_pi0" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->pi0);

    // FBDT efficiency
    w->var("alpha_FBDT_efficiency_uncer")->setConstant(options_->FBDT);

    // FEI
    for (int i = 0; i < options_->NEntryFEI; i++) w->var(("alpha_FEI" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->FEI);

    // qqbar normalization
    w->var("alpha_qq_CAL_UUBAR_uncer")->setConstant(options_->qqbar);
    w->var("alpha_qq_CAL_DDBAR_uncer")->setConstant(options_->qqbar);
    w->var("alpha_qq_CAL_SSBAR_uncer")->setConstant(options_->qqbar);
    w->var("alpha_qq_CAL_CHARM_uncer")->setConstant(options_->qqbar);

    // photon multiplicity correction
    for (int i = 0; i < options_->NEntryMultiplicity; i++) w->var(("alpha_multiplicity" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->multiplicity);

    // B->K form factor
    w->var("alpha_Kff1_uncer")->setConstant(options_->Kff);
    w->var("alpha_Kff2_uncer")->setConstant(options_->Kff);
    w->var("alpha_Kff3_uncer")->setConstant(options_->Kff);

    // B->Kstar form factor
    w->var("alpha_Kstarff1_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff2_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff3_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff4_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff5_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff6_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff7_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff8_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff9_uncer")->setConstant(options_->Kstarff);

    // new B->K form factor
    w->var("alpha_Kff_OLD_uncer")->setConstant(options_->OLD_Kff);

    // fermi motion moment
    w->var("alpha_pf_uncer")->setConstant(options_->pf);

    // K*-Xs transition point
    w->var("alpha_transition_uncer")->setConstant(options_->Transition);

    // b-quark mass
    w->var("alpha_mb_uncer")->setConstant(options_->mb);

    // relative fraction
    w->var("alpha_Kfrac_uncer")->setConstant(options_->fraction);
    w->var("alpha_Kstarfrac_uncer")->setConstant(options_->fraction);

    // MC statistics
    for (int i = 0; i < RarityBins; i++) w->var(("gamma_stat_channel_bin_" + std::to_string(i)).c_str())->setConstant(options_->MCstat);

    // Fragmentation
    for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
        for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
            int temp_index = MxsBin * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma) + Category;
            w->var(("alpha_Xs_fragmentation" + std::to_string(temp_index) + "_uncer").c_str())->setConstant(options_->Fragmentation);
        }
    }

    // mKstar
    w->var("alpha_mKstar_uncer")->setConstant(options_->mKstar);

    // background normalization
    w->var("alpha_mu_CHG")->setConstant(options_->BKGNorm);
    w->var("alpha_mu_MIX")->setConstant(options_->BKGNorm);
    w->var("alpha_mu_UUBAR")->setConstant(options_->BKGNorm);
    w->var("alpha_mu_DDBAR")->setConstant(options_->BKGNorm);
    w->var("alpha_mu_SSBAR")->setConstant(options_->BKGNorm);
    w->var("alpha_mu_CHARM")->setConstant(options_->BKGNorm);

    // BDTc
    w->var("alpha_BDTc_shape_qq")->setConstant(options_->BDTc);

    // BB counting
    w->var("alpha_BB_counting_uncer")->setConstant(options_->BBcounting);

    // BB BR
    for (int i = 0; i < options_->NEntryBR; i++) w->var(("alpha_BBBR" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->BBBR);

    // B->K KL KL BR
    w->var("alpha_BRBtoKKLKL_uncer")->setConstant(options_->BRBtoKKLKL);

    // B->Knn BR
    w->var("alpha_Knn_BR_uncer")->setConstant(options_->BRKnn);
    w->var("alpha_Kstarnn_BR_uncer")->setConstant(options_->BRKnn);
    w->var("alpha_K0nn_BR_uncer")->setConstant(options_->BRKnn);
    w->var("alpha_K0starnn_BR_uncer")->setConstant(options_->BRKnn);

    // B->(D->X KL0) + anything
    w->var("alpha_BtoDtoXKL_uncer")->setConstant(options_->BRDKL0);

    // uncorrelated
    for (int i = 0; i < RarityBins; i++) {
        w->var(("gamma_Signal_MXs1_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
        w->var(("gamma_Signal_MXs2_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
        w->var(("gamma_Signal_MXs3_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
        w->var(("gamma_CHG_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
        w->var(("gamma_MIX_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
        w->var(("gamma_UUBAR_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
        w->var(("gamma_DDBAR_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
        w->var(("gamma_SSBAR_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
        w->var(("gamma_CHARM_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
    }

    // save snapshot
    w->saveSnapshot("ParamValues", *params, true);
}

void Debug(RooWorkspace* w, RooFitResult* fitres, RooDataSet* data) {
    // deprecated. Do not use it before fix the function
    RooArgSet fitargs = fitres->floatParsFinal();
    TIterator* iter(fitargs.createIterator());

    for (TObject* a = iter->Next(); a != 0; a = iter->Next()) {
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(a);
        std::string name = rrv->GetName();
        double val = rrv->getVal();
        double err = rrv->getError();

        if (name == "mu") {
            if (err < 1.5) {

                ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
                RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();
                RooArgSet* obs = (RooArgSet*)mc->GetObservables();
                RooRealVar* x_val = w->var("obs_x_channel");

                // get Category and data
                RooCategory* idx = (RooCategory*)obs->find("channelCat");

                // get expected num of evts for PDFs
                double Signal_Nevts = GetNumEvts(w, "Signal");
                double CHG_Nevts = GetNumEvts(w, "CHG");
                double MIX_Nevts = GetNumEvts(w, "MIX");
                double UUBAR_Nevts = GetNumEvts(w, "UUBAR");
                double DDBAR_Nevts = GetNumEvts(w, "DDBAR");
                double SSBAR_Nevts = GetNumEvts(w, "SSBAR");
                double CHARM_Nevts = GetNumEvts(w, "CHARM");

                // draw
                RooPlot* x_frame = x_val->frame(Title("FBDT"));
                data->plotOn(x_frame, DataError(RooAbsData::Poisson), Cut("channelCat==0"), DrawOption("ZP"), Name("data"));
                //data->plotOn(x_frame, DataError(RooAbsData::Poisson), Cut("channelCat==0"), MarkerSize(0.4), DrawOption("ZP"), Normalization(1, RooAbsReal::ScaleType::NumEvent));
                model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kRed - 6), LineWidth(0), Components("L_x_*Signal*_ShapeSys,L_x_*CHG*_ShapeSys,L_x_*MIX*_ShapeSys,L_x_*UUBAR*_ShapeSys,L_x_*DDBAR*_ShapeSys,L_x_*SSBAR*_ShapeSys,L_x_*CHARM*_ShapeSys"), Normalization(Signal_Nevts + CHG_Nevts + MIX_Nevts + UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("signal"));
                model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kBlue - 6), LineWidth(0), Components("L_x_*CHG*_ShapeSys,L_x_*MIX*_ShapeSys,L_x_*UUBAR*_ShapeSys,L_x_*DDBAR*_ShapeSys,L_x_*SSBAR*_ShapeSys,L_x_*CHARM*_ShapeSys"), Normalization(CHG_Nevts + MIX_Nevts + UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("Charged B"));
                model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kCyan - 6), LineWidth(0), Components("L_x_*MIX*_ShapeSys,L_x_*UUBAR*_ShapeSys,L_x_*DDBAR*_ShapeSys,L_x_*SSBAR*_ShapeSys,L_x_*CHARM*_ShapeSys"), Normalization(MIX_Nevts + UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("Neutral B"));
                model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kOrange - 6), LineWidth(0), Components("L_x_*UUBAR*_ShapeSys,L_x_*DDBAR*_ShapeSys,L_x_*SSBAR*_ShapeSys,L_x_*CHARM*_ShapeSys"), Normalization(UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("Continuum"));
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

                canvas->SaveAs( (std::string("fit_plot_") + "mu-" + std::to_string(val) + "_" + "err-" + std::to_string(err) + ".png").c_str() );

                delete canvas;

            }
        }

    }
}

int main(int argc, char* argv[]) {

    ::ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit"); // default: Minuit Migrad
    ::ROOT::Math::MinimizerOptions::SetDefaultStrategy(1); // default 1

    RooMsgService::instance().setStreamStatus(1, false);
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

    RooRandom::randomGenerator()->SetSeed(rd());

    ReadPIDuncorrsysFile("./KID_cov_remain_truncated.txt", "./PID_cov_remain_truncated.txt");
    ReadBRuncorrsysFile("./BR_cov_remain_truncated.txt");
    Readpi0uncorrsysFile("./pi0_cov_remain_truncated.txt");
    ReadFEIuncorrsysFile("./FEI_cov_remain_truncated.txt");

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
