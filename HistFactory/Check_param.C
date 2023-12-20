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

using namespace RooFit;
using namespace RooStats;
using namespace HistFactory;

using std::string;
using std::to_string;

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

# define RarityBins 6

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
    minim.minos(RooArgSet(*w->var("mu")));

    // fit!
    int status;
    status = minim.minimize(minimizer, algorithm);

    return minim.save();
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
            RooRealVar* norm = w->var(("nom_" + name).c_str());
            double width = std::pow(norm->getValV(), -0.5);
            pulls.push_back((val - 1.0) / width);
            pull_errors.push_back(err / width);
            names.push_back(name);
        }
        else if ((name.find("gamma") != std::string::npos) && (name.find("uncorr") != std::string::npos)) {
            int sample_index = -1;
            int bin_index = -1;

            if (name.find("CHG") != std::string::npos) sample_index = 0;
            else if (name.find("MIX") != std::string::npos) sample_index = 1;
            else if (name.find("UUBAR") != std::string::npos) sample_index = 2;
            else if (name.find("DDBAR") != std::string::npos) sample_index = 3;
            else if (name.find("SSBAR") != std::string::npos) sample_index = 4;
            else if (name.find("CHARM") != std::string::npos) sample_index = 5;
            else if (name.find("Signal") != std::string::npos) sample_index = 6;

            std::vector<std::string> temp_strings = split(name, '_');
            bin_index = stoi(temp_strings.back()); // from 0

            if (name.find("all") != std::string::npos) {
                double KID_uncertainty = weight_KIDsys[RarityBins * sample_index + bin_index];
                double PID_uncertainty = weight_PIDsys[RarityBins * sample_index + bin_index];
                double BR_uncertainty = 0.0;
                double pi0_uncertainty = weight_pi0sys[RarityBins * sample_index + bin_index];
                double FEI_uncertainty = 0.0;
                if ((name.find("CHG") != std::string::npos) || (name.find("MIX") != std::string::npos)) {
                    BR_uncertainty = weight_BRsys[RarityBins * sample_index + bin_index];
                    FEI_uncertainty = weight_FEIsys[RarityBins * sample_index + bin_index];
                }
                else if (name.find("Signal") != std::string::npos) {
                    BR_uncertainty = weight_BRsys[RarityBins * 2 + bin_index]; // exception for signal BB BR uncorrelated uncertainty!
                    FEI_uncertainty = weight_FEIsys[RarityBins * 2 + bin_index]; // exception for signal FEI uncorrelated uncertainty!
                }

                double total_uncertainty = std::sqrt(KID_uncertainty * KID_uncertainty + PID_uncertainty * PID_uncertainty + BR_uncertainty * BR_uncertainty + pi0_uncertainty * pi0_uncertainty + FEI_uncertainty * FEI_uncertainty);
                pulls.push_back((val - 1.0) / total_uncertainty);
                pull_errors.push_back(err / total_uncertainty);
                names.push_back(name);
            }
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

    TCanvas* cpull = new TCanvas("pull_Plot", "pull Plot", 700, 800); cpull->SetBottomMargin(0.3);
    pull_two_sigma->Draw("E2");
    pull_one_sigma->Draw("E2 same");
    pull_ht->Draw("e1 same");
    line->Draw();
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

int Check_param() {

    ::ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit"); // default: Minuit Migrad
    ::ROOT::Math::MinimizerOptions::SetDefaultStrategy(1); // default 1

    ReadPIDuncorrsysFile("./KID_cov_remain_truncated.txt", "./PID_cov_remain_truncated.txt");
    ReadBRuncorrsysFile("./BR_cov_remain_truncated.txt");
    Readpi0uncorrsysFile("./pi0_cov_remain_truncated.txt");
    ReadFEIuncorrsysFile("./FEI_cov_remain_truncated.txt");

    const char* fname = "./PDFandDATA_workspace.root";

    TFile* f = TFile::Open(fname);

    RooWorkspace* w = (RooWorkspace*)f->Get("combined");

    w->Print();
    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();


    // test
    //RooRealVar* alpha = w->var("nom_gamma_stat_channel_bin_0");
    //printf("%lf", alpha->getValV());

    // Lets tell roofit the right names for our histogram variables //
    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x = (RooRealVar*)obs->find("obs_x_channel");
    x->SetTitle("FBDT output");
    x->setUnit("");

    // get Category and data
    RooCategory* idx = (RooCategory*)obs->find("channelCat");
    //RooAbsData* data = (RooAbsData*)w->data("obsData");
    RooDataSet* data = (RooDataSet*)w->data("asimovData");

    // fit
    double eps = 0.001;
    RooAbsReal* nll;
    RooFitResult* fitres = MyMinimizeNLL(w, data, nll, eps);

    RooArgSet fitargs = fitres->floatParsFinal();
    TIterator* iter(fitargs.createIterator());

    // draw pull
    Drawpull(w, iter);

    // define frame
    RooPlot* x_frame = x->frame(Title("Tramsformed FBDT_{1}"));

    // draw
    data->plotOn(x_frame, Name("data_name"), DataError(RooAbsData::SumW2), Cut("channelCat==0"), MarkerSize(0.4), DrawOption("ZP"));
    model->plotOn(x_frame, Name("CHG_name"), Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kRed + 1), LineWidth(0));
    model->plotOn(x_frame, Name("MIX_name"), Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kViolet + 1), LineWidth(0), Components("*MIX*, *UUBAR*, *DDBAR*, *SSBAR*, *CHARM*, *Signal*"));
    model->plotOn(x_frame, Name("UUBAR_name"), Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kBlue + 1), LineWidth(0), Components("*UUBAR*, *DDBAR*, *SSBAR*, *CHARM*, *Signal*"));
    model->plotOn(x_frame, Name("DDBAR_name"), Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kOrange + 1), LineWidth(0), Components("*DDBAR*, *SSBAR*, *CHARM*, *Signal*"));
    model->plotOn(x_frame, Name("SSBAR_name"), Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kGreen + 1), LineWidth(0), Components("*SSBAR*, *CHARM*, *Signal*"));
    model->plotOn(x_frame, Name("CHARM_name"), Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kYellow + 1), LineWidth(0), Components("*CHARM*, *Signal*"));
    model->plotOn(x_frame, Name("signal_name"), Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kPink + 1), LineWidth(0), Components("*Signal*"));
    data->plotOn(x_frame, DataError(RooAbsData::SumW2), Cut("channelCat==0"), MarkerSize(0.4), DrawOption("ZP"));
    //model->paramOn(x_frame);
    TCanvas* cdata = new TCanvas("sPlot", "sPlot demo", 700, 700);
    x_frame->Draw();

    TLegend* leg = new TLegend(0.7, 0.6, 0.9, 0.9);
    leg->SetFillStyle(0);
    leg->SetLineWidth(0);
    leg->AddEntry("data_name", "Data", "LP");
    leg->AddEntry("CHG_name", "charged", "F");
    leg->AddEntry("MIX_name", "mixed", "F");
    leg->AddEntry("UUBAR_name", "u#bar{u}", "F");
    leg->AddEntry("DDBAR_name", "d#bar{d}", "F");
    leg->AddEntry("SSBAR_name", "s#bar{s}", "F");
    leg->AddEntry("CHARM_name", "c#bar{c}", "F");
    leg->AddEntry("signal_name", "signal", "F");
    leg->Draw();

    cdata->SaveAs("FitResult.png");
    /* ======================== CLS ======================== */

    return 0;
}
