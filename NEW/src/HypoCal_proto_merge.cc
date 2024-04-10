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

#include "RooStats/FrequentistCalculator.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterPlot.h"

#include "base.h"

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
		fp = fopen(("CLs_freq_" + std::string(mu) + ".txt").c_str(), "a");
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
		fp = fopen(("CLs_freq_" + std::string(mu) + ".txt").c_str(), "a");
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

int main(int argc, char* argv[]) {

	const double mu_MXs1_max = 10.0;
	const double mu_MXs2_max = 15.0;
	const double mu_MXs3_max = 30.0;

	const double mu_MXs1_step = 0.2;
	const double mu_MXs2_step = 0.2;
	const double mu_MXs3_step = 0.3;

	
	for () {

	}

	std::vector<std::string> names;

	load_files("./", &names, "");

	const char* fname = "./PDFandDATA_workspace.root";

	TFile* f = TFile::Open(fname);

	RooStats::HypoTestInverterResult* result = inverter.GetInterval();


	//GetExpectedCL(result, argv[1]);
	//GetObservedCLs(result, argv[1], 0);
	//GetObservedCLs(result, argv[1], 1);
	//GetObservedCLs(result, argv[1], 2);

	return 0;
}
