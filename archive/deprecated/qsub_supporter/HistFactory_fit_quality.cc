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

// arXiv:1409.4557v2
# define TB0 1.5195 // (Table. 1)
# define TBp 1.6384 // (Table. 1)
# define BR_Kplus_nunubar 0.00000398 // (eq. 10)
# define BR_K0star_nunubar 0.00000919 // (eq. 11)
# define BR_K0_nunubar (BR_Kplus_nunubar*TB0/TBp) // under (eq. 15)
# define BR_Kplusstar_nunubar (BR_K0star_nunubar*TBp/TB0) // under (eq. 15)
# define BR_Xs_nunubar 0.000029 // (eq. 23)
# define BR_Xsu_nonresonant_nunubar (BR_Xs_nunubar - BR_Kplus_nunubar - BR_Kplusstar_nunubar)
# define BR_Xsd_nonresonant_nunubar (BR_Xs_nunubar - BR_K0_nunubar - BR_K0star_nunubar)

// uncertainty of BR from parametric reason
# define Sigma_BR_Kplus_nunubar 0.00000019
# define Sigma_BR_K0star_nunubar 0.00000050
# define Sigma_BR_K0_nunubar (Sigma_BR_Kplus_nunubar*TB0/TBp)
# define Sigma_BR_Kplusstar_nunubar (Sigma_BR_K0star_nunubar*TBp/TB0)
# define Sigma_BR_Xs_nunubar 0.000003

// https://confluence.desy.de/pages/viewpage.action?pageId=107054222
# define N_BpBp_1invab 565400000.0
# define N_B0B0_1invab 534600000.0

# define N_Kplus_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Kplus_nunubar)
# define N_Kplusstar_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Kplusstar_nunubar)
# define N_Xsu_nonresonant_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Xsu_nonresonant_nunubar)
# define N_K0_nunubar_1invab (2.0 * N_B0B0_1invab * BR_K0_nunubar)
# define N_K0star_nunubar_1invab (2.0 * N_B0B0_1invab * BR_K0star_nunubar)
# define N_Xsd_nunubar_1invab (2.0 * N_B0B0_1invab * BR_Xsd_nonresonant_nunubar)

# define Lumi_validation_MC 0.3 // ab-1

// my MC sample number
# define N_Kplus_nunubar 10000000.0
# define N_K0_nunubar 10000000.0
# define N_Kplusstar_nunubar 10000000.0
# define N_K0star_nunubar 10000000.0
# define N_Xsu_nonresonant_nunubar 50000000.0
# define N_Xsd_nonresonant_nunubar 50000000.0

// scale factor for each MC sample
# define Scale_Kplus (N_Kplus_nunubar_1invab/N_Kplus_nunubar)
# define Scale_Kplusstar (N_Kplusstar_nunubar_1invab/N_Kplusstar_nunubar)
# define Scale_Xsu_nonresonant (N_Xsu_nonresonant_nunubar_1invab/N_Xsu_nonresonant_nunubar)
# define Scale_K0 (N_K0_nunubar_1invab/N_K0_nunubar)
# define Scale_K0star (N_K0star_nunubar_1invab/N_K0star_nunubar)
# define Scale_Xsd_nonresonant (N_Xsd_nunubar_1invab/N_Xsd_nonresonant_nunubar)

# define KS0_rel_uncertainty 0.6 // %/cm
# define track_rel_uncertainty 0.24 // %
# define pi0_correction 0.932
# define pi0_rel_uncertainty ((0.0369 / 0.932) * 100.0) // %
# define Kaon_PID_max_uncertainty 0.1 // not percentage. relative uncertainty

# define FEI_cal_Bc_num 12
# define FEI_cal_B0_num 11
double FEI_cal_Bc[FEI_cal_Bc_num] = { 1.04, 0.79, 0.69, 0.56, 0.97, 0.95, 0.74, 0.57, 0.91, 0.51, 0.34, 0.59 };
double FEI_cal_Bc_uncertainty[FEI_cal_Bc_num] = { 0.03, 0.03, 0.05, 0.11, 0.03, 0.03, 0.02, 0.06, 0.1, 0.13, 0.07, 0.02 }; // not relative uncertainty. absolute uncertainty
double FEI_cal_Bc_modeID[FEI_cal_Bc_num] = { 0.0, 1.0, 3.0, 4.0, 15.0, 16.0, 18.0, 19.0, 23.0, 24.0, 30.0, -1.0 };
double FEI_cal_B0[FEI_cal_B0_num] = { 1.16, 0.94, 0.81, 0.79, 0.99, 1.03, 0.67, 0.66, 0.69, 0.49, 0.79 };
double FEI_cal_B0_uncertainty[FEI_cal_B0_num] = { 0.04, 0.05, 0.06, 0.02, 0.03, 0.06, 0.02, 0.03, 0.02, 0.02, 0.12 }; // not relative uncertainty. absolute uncertainty
double FEI_cal_B0_modeID[FEI_cal_B0_num] = { 0.0, 1.0, 3.0, 4.0, 5.0, 15.0, 16.0, 18.0, 19.0, 26.0, -1.0 };

//# define Htransition_Xsu_change (-0.142)
//# define Htransition_Xsd_change (-0.099)
//# define Ltransition_Xsu_change (-0.002)
//# define Ltransition_Xsd_change (0.067)
//# define Hmb_Xsu_change (-0.067)
//# define Hmb_Xsd_change (-0.058)
//# define Lmb_Xsu_change (-0.062)
//# define Lmb_Xsd_change (-0.046)
//# define Hpf_Xsu_change (-0.120)
//# define Hpf_Xsd_change (-0.073)
//# define Lpf_Xsu_change (-0.014)
//# define Lpf_Xsd_change (0.019)

//# define FEI_cal_Bc 1.0
//# define FEI_cal_Bc_uncertainty 0.0
//# define FEI_cal_B0 1.0
//# define FEI_cal_B0_uncertainty 0.0

//# define Htransition_Xsu_change 0.0
//# define Htransition_Xsd_change 0.0
//# define Ltransition_Xsu_change 0.0
//# define Ltransition_Xsd_change 0.0
//# define Hmb_Xsu_change 0.0
//# define Hmb_Xsd_change 0.0
//# define Lmb_Xsu_change 0.0
//# define Lmb_Xsd_change 0.0
//# define Hpf_Xsu_change 0.0
//# define Hpf_Xsd_change 0.0
//# define Lpf_Xsu_change 0.0
//# define Lpf_Xsd_change 0.0

//# define Toy_iter_num 500
//# define LT_iter_num 500

int Toy_iter_num = 0.0;
int LT_iter_num = 0.0;

# define RarityBins 6
# define BinMIN 0.96
# define BinMAX 1.0
double weight_KIDsys[RarityBins * 7] = { 0.0 };
double weight_PIDsys[RarityBins * 7] = { 0.0 };
double weight_BRsys[RarityBins * 3] = { 0.0 };
double weight_pi0sys[RarityBins * 7] = { 0.0 };
double weight_FEIsys[RarityBins * 3] = { 0.0 };

std::random_device rd;
std::default_random_engine generator(rd());

std::vector<std::string> Sample_names = {
    "L_x_Signal_nominal_channel_overallSyst_x_StatUncert_x_channel_Signal_all_uncorr_uncer_ShapeSys",
    "L_x_CHG_nominal_channel_overallSyst_x_StatUncert_x_channel_CHG_all_uncorr_uncer_ShapeSys",
    "L_x_MIX_nominal_channel_overallSyst_x_StatUncert_x_channel_MIX_all_uncorr_uncer_ShapeSys",
    "L_x_UUBAR_nominal_channel_overallSyst_x_StatUncert_x_channel_UUBAR_all_uncorr_uncer_ShapeSys",
    "L_x_DDBAR_nominal_channel_overallSyst_x_StatUncert_x_channel_DDBAR_all_uncorr_uncer_ShapeSys",
    "L_x_SSBAR_nominal_channel_overallSyst_x_StatUncert_x_channel_SSBAR_all_uncorr_uncer_ShapeSys",
    "L_x_CHARM_nominal_channel_overallSyst_x_StatUncert_x_channel_CHARM_all_uncorr_uncer_ShapeSys"
};

void GetNameOfParams(RooWorkspace* w, std::vector<std::string>* names) {
    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    w->loadSnapshot("NominalParamValues");

    RooAbsData* data = (RooAbsData*)w->data("asimovData");
    RooFitResult* fitres = model->fitTo(*data, RooFit::Extended(true), RooFit::SumW2Error(false), PrintLevel(-1), Save());

    RooArgSet fitargs = fitres->floatParsFinal();
    TIterator* iter(fitargs.createIterator());

    for (TObject* a = iter->Next(); a != 0; a = iter->Next()) {
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(a);
        std::string name = rrv->GetName();
        names->push_back(name);
    }

    w->loadSnapshot("NominalParamValues");
}

class FileSaver {
private:
    double* m_true_param;
    double* m_fitting_param;
    double* m_fitting_param_error;
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
    double* m_true_param = nullptr;
    double* m_fitting_param = nullptr;
    double* m_fitting_param_error = nullptr;
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

    // open file
    if (IsItToy) {
        m_file = new TFile(("TOY_result_" + std::to_string(indicator) + ".root").c_str(), "recreate");
        m_file->cd();
        m_tree = new TTree("TOY_result", "");
    }
    else {
        m_file = new TFile(("LT_result_" + std::to_string(mu_injected) + "_" + std::to_string(indicator) + ".root").c_str(), "recreate");
        m_file->cd();
        m_tree = new TTree("LT_result", "");
    }

    // set branches
    for (unsigned int i = 0; i < names->size(); i++) {
        m_tree->Branch((names->at(i) + "_true").c_str(), &m_true_param[i]);
        m_tree->Branch((names->at(i) + "_value").c_str(), &m_fitting_param[i]);
        m_tree->Branch((names->at(i) + "_error").c_str(), &m_fitting_param_error[i]);
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

        for (unsigned int i = 0; i < names->size(); i++) {
            if (name == names->at(i)) {
                m_fitting_param[i] = val;
                m_fitting_param_error[i] = err;
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

RooFitResult* MinimizeNLL(RooWorkspace* w, RooDataSet* data, RooAbsReal* nll, double tolerance = -1.0) { // this function follows the procedure in ProfileLikelihoodTestStat.cxx
    // what we have done
    w->loadSnapshot("NominalParamValues");
    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    // get nll
    RooArgSet* allParams = model->getParameters(*data);
    RooStats::RemoveConstantParameters(allParams);
    RooArgSet fGlobalObs;
    RooArgSet fConditionalObs;
    Bool_t fLOffset = RooStats::IsNLLOffset();
    nll = model->createNLL(*data, RooFit::CloneData(kFALSE), RooFit::Constrain(*allParams), RooFit::GlobalObservables(fGlobalObs), RooFit::ConditionalObservables(fConditionalObs), RooFit::Offset(fLOffset));

    // get default values
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

    w->loadSnapshot("NominalParamValues");
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

    w->var("mu")->setVal(injected_mu);

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

            w->loadSnapshot("NominalParamValues");
            //RooFitResult* fitres = model->fitTo(*genData, RooFit::Extended(true), RooFit::SumW2Error(false), PrintLevel(-1), Save());
            RooAbsReal* nll;
            RooFitResult* fitres = MinimizeNLL(w, genData, nll, eps);

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
        w->loadSnapshot("NominalParamValues");

        //RooFitResult* fitres = model->fitTo(*genData, RooFit::Extended(true), RooFit::SumW2Error(false), PrintLevel(-1), Save());
        RooAbsReal* nll;
        RooFitResult* fitres = MinimizeNLL(w, genData, nll, eps);

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
    if (strcmp(sample_type, "Signal") == 0) index = 0;
    else if (strcmp(sample_type, "CHG") == 0) index = 1;
    else if (strcmp(sample_type, "MIX") == 0) index = 2;
    else if (strcmp(sample_type, "UUBAR") == 0) index = 3;
    else if (strcmp(sample_type, "DDBAR") == 0) index = 4;
    else if (strcmp(sample_type, "SSBAR") == 0) index = 5;
    else if (strcmp(sample_type, "CHARM") == 0) index = 6;
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

    w->loadSnapshot("NominalParamValues");

    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();
    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x_val = w->var("obs_x_channel");

    // get Category and data
    RooCategory* idx = (RooCategory*)obs->find("channelCat");
    RooDataSet* data = (RooDataSet*)w->data("asimovData");

    // fit
    RooAbsReal* nll;
    RooFitResult* fitres = MinimizeNLL(w, data, nll, eps);

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

        if (name == std::string("mu")) {
            printf("mu: %lf +- %lf\n", val, err);
        }

    }

    w->loadSnapshot("NominalParamValues");

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

        w->loadSnapshot("NominalParamValues");

        filesaver.GetTrueValues(w, names);

        RooDataSet* genData = dataPDF.generate(RooArgSet(*x, model->indexCat()), Nevt_data_int, false, true, "", false, true);

        w->loadSnapshot("NominalParamValues");
        RooAbsReal* nll;
        RooFitResult* fitres = MinimizeNLL(w, genData, nll, eps);

        filesaver.GetFittingValues(fitres, names);
        filesaver.GetFittingStatus(fitres);
        filesaver.WriteIntoBranch();

        delete fitres;

    }
    
}

int main(int argc, char* argv[]) {
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
    double injected_mu = -1;
    double eps = -1.0;
    int indicator = 0;

    std::vector<std::string> param_names;
    
    if (std::string(argv[1]) == std::string("ToyMC")) {  // main ToyMC
        if (argc == 5) {
            injected_mu = -1;
            eps = std::atof(argv[2]);
            indicator = std::atoi(argv[3]);
            Toy_iter_num = std::atoi(argv[4]);
        }
        else {
            printf("Toy MC requires 3 arguments {eps} {indicator} {Num of sample}\n");
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
            printf("Linearity test requires only 4 arguments: {injected mu} {eps} {indicator} {Num of sample}\n");
            exit(1);
        }
    }
    else {
        printf("first argument should be {ToyMC|LinearityTest}\n");
        exit(1);
    }

    const char* fname = "./PDFandDATA_workspace.root";

    TFile* f = TFile::Open(fname);

    RooWorkspace* w = (RooWorkspace*)f->Get("combined");

    GetNameOfParams(w, &param_names);

    if (std::string(argv[1]) == std::string("ToyMC")) {
        filesaver.OpenFile(true, &param_names, injected_mu, indicator);
        MyToyMCStudy(w, &param_names, eps, indicator);
    }
    else if (std::string(argv[1]) == std::string("LinearityTest")) {
        filesaver.OpenFile(false, &param_names, injected_mu, indicator);
        MyLinearityTest(w, &param_names, injected_mu, eps, indicator);
    }

    f->Close();

    filesaver.CloseFile();

    return 0;
}
