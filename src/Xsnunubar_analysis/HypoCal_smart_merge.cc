#include <iostream>
#include <stdio.h>
#include <iomanip>

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

#include "RooStats/SamplingDistPlot.h"

#include "base.h"

using namespace RooFit;
using namespace RooStats;
using namespace HistFactory;

using std::string;
using std::to_string;
using std::cout;
using std::endl;

void PrintTestStat(RooStats::HypoTestInverterResult* result, std::string mu_string) {
	TCanvas* c = new TCanvas("c", "c", 696, 472);

	HypoTestInverterPlot* plot = new HypoTestInverterPlot("plot", "plot", result);
	SamplingDistPlot* pl = plot->MakeTestStatPlot(0);
	pl->SetLogYaxis(true);
	pl->Draw();

	c->SaveAs(("TestStat_" + mu_string + ".png").c_str());

	delete plot;
	delete c;

}

int main(int argc, char* argv[]) {

	const double mu_max = 10.0;

	double mu = 0.0;
	RooStats::HypoTestInverterResult* result = nullptr;

	mu = 0.0;
	while(mu < mu_max) {
		std::stringstream stream;
		stream << std::fixed << std::setprecision(1) << mu;
		std::string mu_string = stream.str();

		std::vector<std::string> names;
		load_files("./", &names, ("Hypotestinverter_freq_" + mu_string + "_").c_str());

		for (unsigned int i = 0; i < names.size(); i++) {
			TFile* file = TFile::Open(names.at(i).c_str());
			if (result == nullptr) result = (RooStats::HypoTestInverterResult*)file->Get("result_mu_MXs1");
			else {
				RooStats::HypoTestInverterResult* result_temp = (RooStats::HypoTestInverterResult*)file->Get("result_mu_MXs1");
				result->Add(*result_temp);
			}
			file->Close();
		}

		mu = mu + 0.1;
	}

	result->SetConfidenceLevel(0.9);

	double upperLimit = result->UpperLimit();
	double upperLimitError = result->UpperLimitEstimatedError();

	TCanvas* cHT = new TCanvas("HypoTestInverter Scan", "HypoTestInverter Scan");
	HypoTestInverterPlot* plotHT = new HypoTestInverterPlot("Result_Plot", "HypoTest Scan Result", result);
	cHT->SetLogy(false);
	plotHT->Draw("CLb 2CL");
	cHT->SaveAs("CLs_smart.png");

	std::cout << "Expected upper limits, using the B (alternate) model : " << std::endl;
	std::cout << " expected limit (median) " << result->GetExpectedUpperLimit(0) << std::endl;
	std::cout << " expected limit (-1 sig) " << result->GetExpectedUpperLimit(-1) << std::endl;
	std::cout << " expected limit (+1 sig) " << result->GetExpectedUpperLimit(1) << std::endl;
	std::cout << " expected limit (-2 sig) " << result->GetExpectedUpperLimit(-2) << std::endl;
	std::cout << " expected limit (+2 sig) " << result->GetExpectedUpperLimit(2) << std::endl;
	std::cout << " observed " << upperLimit << std::endl;

	return 0;
}
