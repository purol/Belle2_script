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

#include "RooStats/HybridCalculator.h"
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

# define RarityBins 10

void GetExpectedCL(RooStats::HypoTestInverterResult* fResults, const char* mu) {
	// get CLs CLb CLs+b
	const int nEntries = fResults->ArraySize();

	double p[7];
	double q[7];
	p[0] = ROOT::Math::normal_cdf(-2);
	p[1] = ROOT::Math::normal_cdf(-1);
	p[2] = 0.5;
	p[3] = ROOT::Math::normal_cdf(1);
	p[4] = ROOT::Math::normal_cdf(2);

	int np = 0;
	for (int j = 0; j < nEntries; ++j) {
		int i = j; // i is the order index
		SamplingDistribution* s = fResults->GetExpectedPValueDist(i);
		if (!s)  continue;
		const std::vector<double>& values = s->GetSamplingDistribution();

		double* x = const_cast<double*>(&values[0]); // need to change TMath::Quantiles
		TMath::Quantiles(values.size(), 5, x, q, p, false);

		FILE* fp;
		fp = fopen(("CLs_hyb_" + std::string(mu) + ".txt").c_str(), "a");
		fprintf(fp, "expected CLs median: %lf\n", q[2]);
		fprintf(fp, "expected CLs +1sigma: %lf\n", q[3] - q[2]);
		fprintf(fp, "expected CLs -1sigma: %lf\n", q[2] - q[1]);
		fprintf(fp, "expected CLs +2sigma: %lf\n", q[4] - q[2]);
		fprintf(fp, "expected CLs -2sigma: %lf\n", q[2] - q[0]);
		fclose(fp);
		if (s) delete s;
		np++;
	}
}

void GetObservedCLs(RooStats::HypoTestInverterResult* fResults, const char* mu, int type = 0) {
	// type 0: CLs
	// type 1: CLb
	// type 2: CLs+b
	if (!(type == 0) && !(type == 1) && !(type == 2)) {
		printf("[ERROR] unvalid type!\n");
		exit(1);
	}

	const int nEntries = fResults->ArraySize();

	std::vector<Double_t> xArray;
	std::vector<Double_t> yArray;
	std::vector<Double_t> yErrArray;

	for (int i = 0; i < nEntries; i++) {
		int index = i;

		double CLVal = 0.0;
		double CLErr = 0.0;
		if (type == 0) {
			CLVal = fResults->GetYValue(index);
			CLErr = fResults->GetYError(index);
		}
		else if (type == 1) {
			CLVal = fResults->CLb(index);
			CLErr = fResults->CLbError(index);
		}
		else if (type == 2) {
			CLVal = fResults->CLsplusb(index);
			CLErr = fResults->CLsplusbError(index);
		}

		if (CLVal < 0.0 || !std::isfinite(CLVal)) {
			printf("Got a confidence level of %f at x=%f (failed fit?). Skipping this point.", CLVal, fResults->GetXValue(index));
			continue;
		}

		FILE* fp;
		fp = fopen(("CLs_hyb_" + std::string(mu) + ".txt").c_str(), "a");
		if (type == 0) {
			fprintf(fp, "observed CLs central value: %lf\n", CLVal);
			fprintf(fp, "observed CLs error: %lf\n", CLErr);
		}
		else if (type == 1) {
			fprintf(fp, "observed CLb central value: %lf\n", CLVal);
			fprintf(fp, "observed CLb error: %lf\n", CLErr);
		}
		else if (type == 2) {
			fprintf(fp, "observed CLs+b central value: %lf\n", CLVal);
			fprintf(fp, "observed CLs+b error: %lf\n", CLErr);
		}
		fclose(fp);

		yArray.push_back(CLVal);
		yErrArray.push_back(CLErr);
		xArray.push_back(fResults->GetXValue(index));
	}
}

int main(int argc, char* argv[]) { // argv[1]: mu value to test

	const char* fname = "/home/belle2/junewoo/storage_b1/HistFactory/PDFandDATA_workspace_10_mu.root";

	TFile* f = TFile::Open(fname);

	RooWorkspace* w = (RooWorkspace*)f->Get("combined");

	w->Print();
	ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
	RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

	// Lets tell roofit the right names for our histogram variables //
	RooArgSet* obs = (RooArgSet*)mc->GetObservables();
	RooRealVar* x = (RooRealVar*)obs->find("obs_x_channel");
	x->SetTitle("Tramsformed FBDT_{1}");
	x->setUnit("");

	// get Category and data
	RooCategory* idx = (RooCategory*)obs->find("channelCat");
	RooAbsData* data = (RooAbsData*)w->data("asimovData");

	// fit
	RooFitResult* fitres = model->fitTo(*data, RooFit::Extended(true), RooFit::SumW2Error(false));

	// define frame
	RooPlot* x_frame = x->frame(Title("Tramsformed FBDT_{1}"));

	// draw
	data->plotOn(x_frame, Name("data_name"), DataError(RooAbsData::Poisson), Cut("channelCat==0"), MarkerSize(0.4), DrawOption("ZP"));
	model->plotOn(x_frame, Name("CHG_name"), Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kRed + 1), LineWidth(0));
	model->plotOn(x_frame, Name("MIX_name"), Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kViolet + 1), LineWidth(0), Components("*MIX*, *UUBAR*, *DDBAR*, *SSBAR*, *CHARM*, *Signal*"));
	model->plotOn(x_frame, Name("UUBAR_name"), Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kBlue + 1), LineWidth(0), Components("*UUBAR*, *DDBAR*, *SSBAR*, *CHARM*, *Signal*"));
	model->plotOn(x_frame, Name("DDBAR_name"), Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kOrange + 1), LineWidth(0), Components("*DDBAR*, *SSBAR*, *CHARM*, *Signal*"));
	model->plotOn(x_frame, Name("SSBAR_name"), Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kGreen + 1), LineWidth(0), Components("*SSBAR*, *CHARM*, *Signal*"));
	model->plotOn(x_frame, Name("CHARM_name"), Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kYellow + 1), LineWidth(0), Components("*CHARM*, *Signal*"));
	model->plotOn(x_frame, Name("signal_name"), Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kPink + 1), LineWidth(0), Components("*Signal*"));
	data->plotOn(x_frame, DataError(RooAbsData::Poisson), Cut("channelCat==0"), MarkerSize(0.4), DrawOption("ZP"));

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

	//cdata->SaveAs("/home/belle2/junewoo/storage_b1/HistFactory/plot/FitResult.png");
	/* ======================== CLS ======================== */

	if (data->isWeighted()) {
		Info("Hypocal", "Data set is weighted, nentries = %d and sum of weights = %8.1f.\n", data->numEntries(), data->sumEntries());
	}

	RooStats::ModelConfig* sbModel = (RooStats::ModelConfig*)w->obj("ModelConfig");
	RooStats::ModelConfig* bModel = (RooStats::ModelConfig*)sbModel->Clone("BonlyModel");
	RooRealVar* poi = (RooRealVar*)bModel->GetParametersOfInterest()->first();
	poi->setVal(0);
	bModel->SetSnapshot(*poi);

	RooStats::HybridCalculator HybCalc(*data, *bModel, *sbModel);
	RooStats::ProfileLikelihoodTestStat* plr = new RooStats::ProfileLikelihoodTestStat(*sbModel->GetPdf());
	plr->SetOneSided(true);

	RooStats::ToyMCSampler* toymcs = (RooStats::ToyMCSampler*)HybCalc.GetTestStatSampler();
	toymcs->SetTestStatistic(plr);
	HybCalc.SetToys(5000, 5000);

	RooStats::HypoTestInverter inverter(HybCalc);
	//inverter.SetConfidenceLevel(0.90);
	inverter.UseCLs(true);
	inverter.SetVerbose(false);
	inverter.SetFixedScan(1, std::stof(argv[1]), std::stof(argv[1])); // set number of points , xmin and xmax

	TStopwatch sw;
	sw.Start();

	RooStats::HypoTestInverterResult* result = inverter.GetInterval();

	sw.Stop();
	printf("consumed time: %lf (s)\n", sw.RealTime());

	GetExpectedCL(result, argv[1]);
	GetObservedCLs(result, argv[1], 0);
	GetObservedCLs(result, argv[1], 1);
	GetObservedCLs(result, argv[1], 2);

	return 0;
}
