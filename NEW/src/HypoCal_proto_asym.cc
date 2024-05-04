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

#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterPlot.h"

using namespace RooFit;
using namespace RooStats;
using namespace HistFactory;

using std::string;
using std::to_string;
using std::cout;
using std::endl;

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
		fp = fopen(("CLs_asym_" + std::string(mu) + ".txt").c_str(), "a");
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
		fp = fopen(("CLs_asym_" + std::string(mu) + ".txt").c_str(), "a");
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

int main(int argc, char* argv[]) { // argv[1]: mu value to test, argv[2]: eps

	::ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2"); // default: Minuit Migrad
	::ROOT::Math::MinimizerOptions::SetDefaultStrategy(1); // default 1
	double eps = std::atof(argv[2]);
	::ROOT::Math::MinimizerOptions::SetDefaultTolerance(eps); // default 0.01. but it is better to use 0.001

	const char* fname = "./PDFandDATA_workspace.root";

	TFile* f = TFile::Open(fname);

	RooWorkspace* w = (RooWorkspace*)f->Get("combined");

	w->Print();
	ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
	RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

	// Lets tell roofit the right names for our histogram variables //
	RooArgSet* obs = (RooArgSet*)mc->GetObservables();
	RooRealVar* x = (RooRealVar*)obs->find("obs_x_channel");
	x->SetTitle("FBDT");
	x->setUnit("");

	// get Category and data
	RooCategory* idx = (RooCategory*)obs->find("channelCat");
	RooAbsData* data = (RooAbsData*)w->data("asimovData");

	// fit
	RooFitResult* fitres = model->fitTo(*data, RooFit::Extended(true), RooFit::SumW2Error(false));

	// define frame
	RooPlot* x_frame = x->frame(Title("FBDT"));

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
	// try to make 0 for mu_MXs1, mu_MXs2, and mu_MXs3
	RooRealVar* mu_MXs1 = (RooRealVar*)bModel->GetNuisanceParameters()->find("mu_MXs1");
	if (mu_MXs1 != nullptr) mu_MXs1->setVal(0);
	RooRealVar* mu_MXs2 = (RooRealVar*)bModel->GetNuisanceParameters()->find("mu_MXs2");
	if (mu_MXs2 != nullptr) mu_MXs2->setVal(0);
	RooRealVar* mu_MXs3 = (RooRealVar*)bModel->GetNuisanceParameters()->find("mu_MXs3");
	if (mu_MXs3 != nullptr) mu_MXs3->setVal(0);
	bModel->SetSnapshot(*poi);

	RooStats::AsymptoticCalculator AsymCalc(*data, *bModel, *sbModel);
	RooStats::ProfileLikelihoodTestStat* plr = new RooStats::ProfileLikelihoodTestStat(*sbModel->GetPdf());
	plr->SetOneSided(true);
	plr->SetMinimizer("Minuit2");
	plr->SetStrategy(1);
	//plr->SetLOffset(true);
	plr->SetTolerance(eps);

	RooStats::ToyMCSampler* toymcs = (RooStats::ToyMCSampler*)AsymCalc.GetTestStatSampler();

	RooStats::HypoTestInverter inverter(AsymCalc);
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
