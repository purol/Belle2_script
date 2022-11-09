#include <iostream>
#include <stdio.h>
#include <string>
#include <algorithm>
#include <vector>
#include <numeric>
#include <time.h>
#include <random>

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
# define track_rel_uncertainty 0.69 // %
# define pi0_correction 0.932
# define pi0_rel_uncertainty ((0.0369 / 0.932) * 100.0) // %
# define Kaon_PID_max_uncertainty 0.1 // not percentage. relative uncertainty
// https://indico.belle2.org/event/6872/contributions/37447/attachments/17127/25504/FEIperformance_B2GM.pdf
# define FEI_cal_Bc 0.679
# define FEI_cal_Bc_uncertainty (0.017/FEI_cal_Bc) // not percentage. relative uncertainty
# define FEI_cal_B0 0.713
# define FEI_cal_B0_uncertainty (0.019/FEI_cal_B0) // not percentage. relative uncertainty

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

# define Toy_iter_num 5000
# define LT_iter_num 5000

std::default_random_engine generator;

std::vector<std::string> Sample_names = {
    "L_x_Signal_nominal_channel_overallSyst_x_StatUncert",
    "L_x_CHG_nominal_channel_overallSyst_x_StatUncert",
    "L_x_MIX_nominal_channel_overallSyst_x_StatUncert",
    "L_x_UUBAR_nominal_channel_overallSyst_x_StatUncert",
    "L_x_DDBAR_nominal_channel_overallSyst_x_StatUncert",
    "L_x_SSBAR_nominal_channel_overallSyst_x_StatUncert",
    "L_x_CHARM_nominal_channel_overallSyst_x_StatUncert"
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

void MyToyMCStudy(RooWorkspace *w, std::vector<std::string>* names, double eps){

        std::vector<double> mus;
        std::vector<double> mu_errors;
        std::vector<double> mu_pulls;

        std::vector<double> alpha_mu_CHGs;
        std::vector<double> alpha_mu_CHG_errors;
        std::vector<double> alpha_mu_CHG_pulls;
        std::vector<double> mu_CHGs;

        std::vector<double> alpha_mu_MIXs;
        std::vector<double> alpha_mu_MIX_errors;
        std::vector<double> alpha_mu_MIX_pulls;
        std::vector<double> mu_MIXs;

        std::vector<double> alpha_mu_UUBARs;
        std::vector<double> alpha_mu_UUBAR_errors;
        std::vector<double> alpha_mu_UUBAR_pulls;
        std::vector<double> mu_UUBARs;

        std::vector<double> alpha_mu_DDBARs;
        std::vector<double> alpha_mu_DDBAR_errors;
        std::vector<double> alpha_mu_DDBAR_pulls;
        std::vector<double> mu_DDBARs;

        std::vector<double> alpha_mu_SSBARs;
        std::vector<double> alpha_mu_SSBAR_errors;
        std::vector<double> alpha_mu_SSBAR_pulls;
        std::vector<double> mu_SSBARs;

        std::vector<double> alpha_mu_CHARMs;
        std::vector<double> alpha_mu_CHARM_errors;
        std::vector<double> alpha_mu_CHARM_pulls;
        std::vector<double> mu_CHARMs;

        ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
        RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

        RooArgSet* obs = (RooArgSet*)mc->GetObservables();
        RooRealVar* x = (RooRealVar*)obs->find("obs_x_channel");

        for(int i=0; i< Toy_iter_num; i++) { // Do Toy MC study
            
            double Nevt_total = SetParamsForToy(w, names, 1.0);

            RooDataSet* genData = model->generate(RooArgSet(*x,model->indexCat()), Nevt_total, false, true, "", false, true);

            w->loadSnapshot("NominalParamValues");
            //RooFitResult* fitres = model->fitTo(*genData, RooFit::Extended(true), RooFit::SumW2Error(false), PrintLevel(-1), Save());
            RooAbsReal* nll;
            RooFitResult* fitres = MinimizeNLL(w, genData, nll, eps);

            RooArgSet fitargs_TOY = fitres->floatParsFinal();
            TIterator* iter_TOY(fitargs_TOY.createIterator());

            for(TObject* a_TOY = iter_TOY->Next(); a_TOY != 0; a_TOY = iter_TOY->Next()){
                RooRealVar* rrv_TOY = dynamic_cast<RooRealVar*>(a_TOY);
                std::string name_TOY = rrv_TOY->GetName();
                double val_TOY = rrv_TOY->getVal();
                double err_TOY = rrv_TOY->getError();
                if(name_TOY == std::string("mu")){
                    mus.push_back(val_TOY);
                    mu_errors.push_back(err_TOY);
                    mu_pulls.push_back((val_TOY - 1.0)/err_TOY);
                }
                else if (name_TOY == std::string("alpha_mu_CHG")) {
                    alpha_mu_CHGs.push_back(val_TOY);
                    alpha_mu_CHG_errors.push_back(err_TOY);
                    alpha_mu_CHG_pulls.push_back((val_TOY - w->var("nom_alpha_mu_CHG")->getValV()) / err_TOY);
                }
                else if (name_TOY == std::string("alpha_mu_MIX")) {
                    alpha_mu_MIXs.push_back(val_TOY);
                    alpha_mu_MIX_errors.push_back(err_TOY);
                    alpha_mu_MIX_pulls.push_back((val_TOY - w->var("nom_alpha_mu_MIX")->getValV()) / err_TOY);
                }
                else if (name_TOY == std::string("alpha_mu_UUBAR")) {
                    alpha_mu_UUBARs.push_back(val_TOY);
                    alpha_mu_UUBAR_errors.push_back(err_TOY);
                    alpha_mu_UUBAR_pulls.push_back((val_TOY - w->var("nom_alpha_mu_UUBAR")->getValV()) / err_TOY);
                }
                else if (name_TOY == std::string("alpha_mu_DDBAR")) {
                    alpha_mu_DDBARs.push_back(val_TOY);
                    alpha_mu_DDBAR_errors.push_back(err_TOY);
                    alpha_mu_DDBAR_pulls.push_back((val_TOY - w->var("nom_alpha_mu_DDBAR")->getValV()) / err_TOY);
                }
                else if (name_TOY == std::string("alpha_mu_SSBAR")) {
                    alpha_mu_SSBARs.push_back(val_TOY);
                    alpha_mu_SSBAR_errors.push_back(err_TOY);
                    alpha_mu_SSBAR_pulls.push_back((val_TOY - w->var("nom_alpha_mu_SSBAR")->getValV()) / err_TOY);
                }
                else if (name_TOY == std::string("alpha_mu_CHARM")) {
                    alpha_mu_CHARMs.push_back(val_TOY);
                    alpha_mu_CHARM_errors.push_back(err_TOY);
                    alpha_mu_CHARM_pulls.push_back((val_TOY - w->var("nom_alpha_mu_CHARM")->getValV()) / err_TOY);
                }
            }

            mu_CHGs.push_back(w->function("CHG_nominal_channel_epsilon")->getVal());
            mu_MIXs.push_back(w->function("MIX_nominal_channel_epsilon")->getVal());
            mu_UUBARs.push_back(w->function("UUBAR_nominal_channel_epsilon")->getVal());
            mu_DDBARs.push_back(w->function("DDBAR_nominal_channel_epsilon")->getVal());
            mu_SSBARs.push_back(w->function("SSBAR_nominal_channel_epsilon")->getVal());
            mu_CHARMs.push_back(w->function("CHARM_nominal_channel_epsilon")->getVal());

            delete nll;
            delete fitres;

        }

        // define file to save results
        TFile* temp_file = new TFile("TOY_result.root", "recreate");
        temp_file->cd();
        TTree* temp_tree = new TTree("TOY_result", "");

        double mu = -1;
        double mu_error = -1;
        double mu_pull = -1;

        double alpha_mu_CHG = -1;
        double alpha_mu_CHG_error = -1;
        double alpha_mu_CHG_pull = -1;
        double mu_CHG = -1;

        double alpha_mu_MIX = -1;
        double alpha_mu_MIX_error = -1;
        double alpha_mu_MIX_pull = -1;
        double mu_MIX = -1;

        double alpha_mu_UUBAR = -1;
        double alpha_mu_UUBAR_error = -1;
        double alpha_mu_UUBAR_pull = -1;
        double mu_UUBAR = -1;

        double alpha_mu_DDBAR = -1;
        double alpha_mu_DDBAR_error = -1;
        double alpha_mu_DDBAR_pull = -1;
        double mu_DDBAR = -1;

        double alpha_mu_SSBAR = -1;
        double alpha_mu_SSBAR_error = -1;
        double alpha_mu_SSBAR_pull = -1;
        double mu_SSBAR = -1;

        double alpha_mu_CHARM = -1;
        double alpha_mu_CHARM_error = -1;
        double alpha_mu_CHARM_pull = -1;
        double mu_CHARM = -1;

        temp_tree->Branch("mu", &mu);
        temp_tree->Branch("mu_error", &mu_error);
        temp_tree->Branch("mu_pull", &mu_pull);

        temp_tree->Branch("alpha_mu_CHG", &alpha_mu_CHG);
        temp_tree->Branch("alpha_mu_CHG_error", &alpha_mu_CHG_error);
        temp_tree->Branch("alpha_mu_CHG_pull", &alpha_mu_CHG_pull);
        temp_tree->Branch("mu_CHG", &mu_CHG);

        temp_tree->Branch("alpha_mu_MIX", &alpha_mu_MIX);
        temp_tree->Branch("alpha_mu_MIX_error", &alpha_mu_MIX_error);
        temp_tree->Branch("alpha_mu_MIX_pull", &alpha_mu_MIX_pull);
        temp_tree->Branch("mu_MIX", &mu_MIX);

        temp_tree->Branch("alpha_mu_UUBAR", &alpha_mu_UUBAR);
        temp_tree->Branch("alpha_mu_UUBAR_error", &alpha_mu_UUBAR_error);
        temp_tree->Branch("alpha_mu_UUBAR_pull", &alpha_mu_UUBAR_pull);
        temp_tree->Branch("mu_UUBAR", &mu_UUBAR);

        temp_tree->Branch("alpha_mu_DDBAR", &alpha_mu_DDBAR);
        temp_tree->Branch("alpha_mu_DDBAR_error", &alpha_mu_DDBAR_error);
        temp_tree->Branch("alpha_mu_DDBAR_pull", &alpha_mu_DDBAR_pull);
        temp_tree->Branch("mu_DDBAR", &mu_DDBAR);

        temp_tree->Branch("alpha_mu_SSBAR", &alpha_mu_SSBAR);
        temp_tree->Branch("alpha_mu_SSBAR_error", &alpha_mu_SSBAR_error);
        temp_tree->Branch("alpha_mu_SSBAR_pull", &alpha_mu_SSBAR_pull);
        temp_tree->Branch("mu_SSBAR", &mu_SSBAR);

        temp_tree->Branch("alpha_mu_CHARM", &alpha_mu_CHARM);
        temp_tree->Branch("alpha_mu_CHARM_error", &alpha_mu_CHARM_error);
        temp_tree->Branch("alpha_mu_CHARM_pull", &alpha_mu_CHARM_pull);
        temp_tree->Branch("mu_CHARM", &mu_CHARM);

        for (int i = 0; i < Toy_iter_num; i++) {
            mu = mus.at(i);
            mu_error = mu_errors.at(i);
            mu_pull = mu_pulls.at(i);

            alpha_mu_CHG = alpha_mu_CHGs.at(i);
            alpha_mu_CHG_error = alpha_mu_CHG_errors.at(i);
            alpha_mu_CHG_pull = alpha_mu_CHG_pulls.at(i);
            mu_CHG = mu_CHGs.at(i);

            alpha_mu_MIX = alpha_mu_MIXs.at(i);
            alpha_mu_MIX_error = alpha_mu_MIX_errors.at(i);
            alpha_mu_MIX_pull = alpha_mu_MIX_pulls.at(i);
            mu_MIX = mu_MIXs.at(i);

            alpha_mu_UUBAR = alpha_mu_UUBARs.at(i);
            alpha_mu_UUBAR_error = alpha_mu_UUBAR_errors.at(i);
            alpha_mu_UUBAR_pull = alpha_mu_UUBAR_pulls.at(i);
            mu_UUBAR = mu_UUBARs.at(i);

            alpha_mu_DDBAR = alpha_mu_DDBARs.at(i);
            alpha_mu_DDBAR_error = alpha_mu_DDBAR_errors.at(i);
            alpha_mu_DDBAR_pull = alpha_mu_DDBAR_pulls.at(i);
            mu_DDBAR = mu_DDBARs.at(i);

            alpha_mu_SSBAR = alpha_mu_SSBARs.at(i);
            alpha_mu_SSBAR_error = alpha_mu_SSBAR_errors.at(i);
            alpha_mu_SSBAR_pull = alpha_mu_SSBAR_pulls.at(i);
            mu_SSBAR = mu_SSBARs.at(i);

            alpha_mu_CHARM = alpha_mu_CHARMs.at(i);
            alpha_mu_CHARM_error = alpha_mu_CHARM_errors.at(i);
            alpha_mu_CHARM_pull = alpha_mu_CHARM_pulls.at(i);
            mu_CHARM = mu_CHARMs.at(i);

            temp_tree->Fill();
        }

        temp_file->cd();
        temp_tree->Write();
        temp_file->Close();
}

void MyLinearityTest(RooWorkspace* w, std::vector<std::string>* names, double mu_injected, double eps) {

    std::vector<double> mus;
    std::vector<double> mu_errors;
    std::vector<double> mu_pulls;

    std::vector<double> alpha_mu_CHGs;
    std::vector<double> alpha_mu_CHG_errors;
    std::vector<double> alpha_mu_CHG_pulls;
    std::vector<double> mu_CHGs;

    std::vector<double> alpha_mu_MIXs;
    std::vector<double> alpha_mu_MIX_errors;
    std::vector<double> alpha_mu_MIX_pulls;
    std::vector<double> mu_MIXs;

    std::vector<double> alpha_mu_UUBARs;
    std::vector<double> alpha_mu_UUBAR_errors;
    std::vector<double> alpha_mu_UUBAR_pulls;
    std::vector<double> mu_UUBARs;

    std::vector<double> alpha_mu_DDBARs;
    std::vector<double> alpha_mu_DDBAR_errors;
    std::vector<double> alpha_mu_DDBAR_pulls;
    std::vector<double> mu_DDBARs;

    std::vector<double> alpha_mu_SSBARs;
    std::vector<double> alpha_mu_SSBAR_errors;
    std::vector<double> alpha_mu_SSBAR_pulls;
    std::vector<double> mu_SSBARs;

    std::vector<double> alpha_mu_CHARMs;
    std::vector<double> alpha_mu_CHARM_errors;
    std::vector<double> alpha_mu_CHARM_pulls;
    std::vector<double> mu_CHARMs;

    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x = (RooRealVar*)obs->find("obs_x_channel");

    for (int i = 0; i < LT_iter_num; i++) { // Do LT MC study

        double Nevt_total = SetParamsForToy(w, names, mu_injected);

        RooDataSet* genData = model->generate(RooArgSet(*x, model->indexCat()), Nevt_total, false, true, "", false, true);
        w->loadSnapshot("NominalParamValues");

        //RooFitResult* fitres = model->fitTo(*genData, RooFit::Extended(true), RooFit::SumW2Error(false), PrintLevel(-1), Save());
        RooAbsReal* nll
        RooFitResult* fitres = MinimizeNLL(w, genData, nll, eps);
        delete nll;

        RooArgSet fitargs_LT = fitres->floatParsFinal();
        TIterator* iter_LT(fitargs_LT.createIterator());

        for (TObject* a_LT = iter_LT->Next(); a_LT != 0; a_LT = iter_LT->Next()) {
            RooRealVar* rrv_LT = dynamic_cast<RooRealVar*>(a_LT);
            std::string name_LT = rrv_LT->GetName();
            double val_LT = rrv_LT->getVal();
            double err_LT = rrv_LT->getError();
            if (name_LT == std::string("mu")) {
                mus.push_back(val_LT);
                mu_errors.push_back(err_LT);
                mu_pulls.push_back((val_LT - mu_injected) / err_LT);
            }
            else if (name_LT == std::string("alpha_mu_CHG")) {
                alpha_mu_CHGs.push_back(val_LT);
                alpha_mu_CHG_errors.push_back(err_LT);
                alpha_mu_CHG_pulls.push_back((val_LT - w->var("nom_alpha_mu_CHG")->getValV()) / err_LT);
            }
            else if (name_LT == std::string("alpha_mu_MIX")) {
                alpha_mu_MIXs.push_back(val_LT);
                alpha_mu_MIX_errors.push_back(err_LT);
                alpha_mu_MIX_pulls.push_back((val_LT - w->var("nom_alpha_mu_MIX")->getValV()) / err_LT);
            }
            else if (name_LT == std::string("alpha_mu_UUBAR")) {
                alpha_mu_UUBARs.push_back(val_LT);
                alpha_mu_UUBAR_errors.push_back(err_LT);
                alpha_mu_UUBAR_pulls.push_back((val_LT - w->var("nom_alpha_mu_UUBAR")->getValV()) / err_LT);
            }
            else if (name_LT == std::string("alpha_mu_DDBAR")) {
                alpha_mu_DDBARs.push_back(val_LT);
                alpha_mu_DDBAR_errors.push_back(err_LT);
                alpha_mu_DDBAR_pulls.push_back((val_LT - w->var("nom_alpha_mu_DDBAR")->getValV()) / err_LT);
            }
            else if (name_LT == std::string("alpha_mu_SSBAR")) {
                alpha_mu_SSBARs.push_back(val_LT);
                alpha_mu_SSBAR_errors.push_back(err_LT);
                alpha_mu_SSBAR_pulls.push_back((val_LT - w->var("nom_alpha_mu_SSBAR")->getValV()) / err_LT);
            }
            else if (name_LT == std::string("alpha_mu_CHARM")) {
                alpha_mu_CHARMs.push_back(val_LT);
                alpha_mu_CHARM_errors.push_back(err_LT);
                alpha_mu_CHARM_pulls.push_back((val_LT - w->var("nom_alpha_mu_CHARM")->getValV()) / err_LT);
            }
        }

        mu_CHGs.push_back(w->function("CHG_nominal_channel_epsilon")->getVal());
        mu_MIXs.push_back(w->function("MIX_nominal_channel_epsilon")->getVal());
        mu_UUBARs.push_back(w->function("UUBAR_nominal_channel_epsilon")->getVal());
        mu_DDBARs.push_back(w->function("DDBAR_nominal_channel_epsilon")->getVal());
        mu_SSBARs.push_back(w->function("SSBAR_nominal_channel_epsilon")->getVal());
        mu_CHARMs.push_back(w->function("CHARM_nominal_channel_epsilon")->getVal());

        delete nll;
        delete fitres;

    }

    // define file to save results
    TFile* temp_file = new TFile( ("LT_result_" + std::to_string(mu_injected) + ".root").c_str() , "recreate");
    temp_file->cd();
    TTree* temp_tree = new TTree("LT_result", "");

    double mu = -1;
    double mu_error = -1;
    double mu_pull = -1;

    double alpha_mu_CHG = -1;
    double alpha_mu_CHG_error = -1;
    double alpha_mu_CHG_pull = -1;
    double mu_CHG = -1;

    double alpha_mu_MIX = -1;
    double alpha_mu_MIX_error = -1;
    double alpha_mu_MIX_pull = -1;
    double mu_MIX = -1;

    double alpha_mu_UUBAR = -1;
    double alpha_mu_UUBAR_error = -1;
    double alpha_mu_UUBAR_pull = -1;
    double mu_UUBAR = -1;

    double alpha_mu_DDBAR = -1;
    double alpha_mu_DDBAR_error = -1;
    double alpha_mu_DDBAR_pull = -1;
    double mu_DDBAR = -1;

    double alpha_mu_SSBAR = -1;
    double alpha_mu_SSBAR_error = -1;
    double alpha_mu_SSBAR_pull = -1;
    double mu_SSBAR = -1;

    double alpha_mu_CHARM = -1;
    double alpha_mu_CHARM_error = -1;
    double alpha_mu_CHARM_pull = -1;
    double mu_CHARM = -1;

    temp_tree->Branch("mu", &mu);
    temp_tree->Branch("mu_error", &mu_error);
    temp_tree->Branch("mu_pull", &mu_pull);

    temp_tree->Branch("alpha_mu_CHG", &alpha_mu_CHG);
    temp_tree->Branch("alpha_mu_CHG_error", &alpha_mu_CHG_error);
    temp_tree->Branch("alpha_mu_CHG_pull", &alpha_mu_CHG_pull);
    temp_tree->Branch("mu_CHG", &mu_CHG);

    temp_tree->Branch("alpha_mu_MIX", &alpha_mu_MIX);
    temp_tree->Branch("alpha_mu_MIX_error", &alpha_mu_MIX_error);
    temp_tree->Branch("alpha_mu_MIX_pull", &alpha_mu_MIX_pull);
    temp_tree->Branch("mu_MIX", &mu_MIX);

    temp_tree->Branch("alpha_mu_UUBAR", &alpha_mu_UUBAR);
    temp_tree->Branch("alpha_mu_UUBAR_error", &alpha_mu_UUBAR_error);
    temp_tree->Branch("alpha_mu_UUBAR_pull", &alpha_mu_UUBAR_pull);
    temp_tree->Branch("mu_UUBAR", &mu_UUBAR);

    temp_tree->Branch("alpha_mu_DDBAR", &alpha_mu_DDBAR);
    temp_tree->Branch("alpha_mu_DDBAR_error", &alpha_mu_DDBAR_error);
    temp_tree->Branch("alpha_mu_DDBAR_pull", &alpha_mu_DDBAR_pull);
    temp_tree->Branch("mu_DDBAR", &mu_DDBAR);

    temp_tree->Branch("alpha_mu_SSBAR", &alpha_mu_SSBAR);
    temp_tree->Branch("alpha_mu_SSBAR_error", &alpha_mu_SSBAR_error);
    temp_tree->Branch("alpha_mu_SSBAR_pull", &alpha_mu_SSBAR_pull);
    temp_tree->Branch("mu_SSBAR", &mu_SSBAR);

    temp_tree->Branch("alpha_mu_CHARM", &alpha_mu_CHARM);
    temp_tree->Branch("alpha_mu_CHARM_error", &alpha_mu_CHARM_error);
    temp_tree->Branch("alpha_mu_CHARM_pull", &alpha_mu_CHARM_pull);
    temp_tree->Branch("mu_CHARM", &mu_CHARM);

    for (int i = 0; i < LT_iter_num; i++) {
        mu = mus.at(i);
        mu_error = mu_errors.at(i);
        mu_pull = mu_pulls.at(i);

        alpha_mu_CHG = alpha_mu_CHGs.at(i);
        alpha_mu_CHG_error = alpha_mu_CHG_errors.at(i);
        alpha_mu_CHG_pull = alpha_mu_CHG_pulls.at(i);
        mu_CHG = mu_CHGs.at(i);

        alpha_mu_MIX = alpha_mu_MIXs.at(i);
        alpha_mu_MIX_error = alpha_mu_MIX_errors.at(i);
        alpha_mu_MIX_pull = alpha_mu_MIX_pulls.at(i);
        mu_MIX = mu_MIXs.at(i);

        alpha_mu_UUBAR = alpha_mu_UUBARs.at(i);
        alpha_mu_UUBAR_error = alpha_mu_UUBAR_errors.at(i);
        alpha_mu_UUBAR_pull = alpha_mu_UUBAR_pulls.at(i);
        mu_UUBAR = mu_UUBARs.at(i);

        alpha_mu_DDBAR = alpha_mu_DDBARs.at(i);
        alpha_mu_DDBAR_error = alpha_mu_DDBAR_errors.at(i);
        alpha_mu_DDBAR_pull = alpha_mu_DDBAR_pulls.at(i);
        mu_DDBAR = mu_DDBARs.at(i);

        alpha_mu_SSBAR = alpha_mu_SSBARs.at(i);
        alpha_mu_SSBAR_error = alpha_mu_SSBAR_errors.at(i);
        alpha_mu_SSBAR_pull = alpha_mu_SSBAR_pulls.at(i);
        mu_SSBAR = mu_SSBARs.at(i);

        alpha_mu_CHARM = alpha_mu_CHARMs.at(i);
        alpha_mu_CHARM_error = alpha_mu_CHARM_errors.at(i);
        alpha_mu_CHARM_pull = alpha_mu_CHARM_pulls.at(i);
        mu_CHARM = mu_CHARMs.at(i);

        temp_tree->Fill();
    }

    temp_file->cd();
    temp_tree->Write();
    temp_file->Close();
}

int main(int argc, char* argv[]) {
    RooMsgService::instance().setStreamStatus(1, false);
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

    RooRandom::randomGenerator()->SetSeed(time(NULL));

    // argv[1]: {ToyMC|LinearityTest}
    // argv[2]: injected mu when Linearity test
    // argv[3]: eps for minimizer
    double injected_mu = -1;
    double eps = -1.0;

    std::vector<std::string> param_names;
    
    if (std::string(argv[1]) == std::string("ToyMC")) {  // main ToyMC
        if (argc == 3) {
            injected_mu = -1;
            eps = std::atof(argv[2]);
        }
        else {
            printf("Toy MC requires 1 arguments\n");
            exit(1);
        }
    }
    else if (std::string(argv[1]) == std::string("LinearityTest")) { // main LinearityTest 12 1234 0.8
        if (argc == 4) {
            injected_mu = std::atof(argv[2]);
            eps = std::atof(argv[3]);
        }
        else {
            printf("Linearity test requires only 2 arguments: {injected mu} {eps}\n");
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
        MyToyMCStudy(w, &param_names, eps);
    }
    else if (std::string(argv[1]) == std::string("LinearityTest")) {
        MyLinearityTest(w, &param_names, injected_mu, eps);
    }

    f->Close();

    return 0;
}
