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
#include "RooStats/HistFactory/Sample.h"

#include "template.h"
#include "correctors.h"
#include "constants.h"

using namespace RooFit;
using namespace RooStats;
using namespace HistFactory;

using std::string;
using std::to_string;

Corrector_Fragmentation corrector_Fragmentation;

int WorkSpace_modify() {
	// this code make several model with different \mu values

	const char* fname = "./PDFandDATA_workspace.root";
	TFile* f = TFile::Open(fname);
	RooWorkspace* w = (RooWorkspace*)f->Get("combined");

	RooStats::ModelConfig* sbModel = (RooStats::ModelConfig*)w->obj("ModelConfig");

	for (int i = 0; i < 50; i++) {
		RooStats::ModelConfig* temp_Model = (RooStats::ModelConfig*)sbModel->Clone(("ModelConfig_" + std::to_string(i)).c_str());
		printf("print!\n");

		const double BR_1 = 0.0000048514;
		const double BR_2 = 0.0000085024;
		const double BR_3 = 0.0000156653;
		const double BR_total = 0.000029;

		RooRealVar* mu_MXs1 = w->var("mu_MXs1");
		RooRealVar* mu_MXs2 = w->var("mu_MXs2");
		RooRealVar* mu_MXs3 = w->var("mu_MXs3");

		mu_MXs1->setConstant(true);
		RooRealVar* mu = new RooRealVar("mu", "mu", -100.0, 100.0);
		RooRealVar* mu1 = new RooFormulaVar("mu_MXs1_alt", "(@2*0.000029-@0*0.0000085024-@1*0.0000156653)/0.0000048514", RooArgList(*mu_MXs2, *mu_MXs3, *mu));

		RooProduct* MXs1_MXs1_scale = (RooProduct*)w->function("Signal_MXs1_nominal_MXs1_channel_MXs1_scaleFactors");
		MXs1_MXs1_scale->addTerm(mu1);

		RooProduct* MXs1_MXs2_scale = (RooProduct*)w->function("Signal_MXs1_nominal_MXs2_channel_MXs2_scaleFactors");
		MXs1_MXs2_scale->addTerm(mu1);

		RooProduct* MXs1_MXs3_scale = (RooProduct*)w->function("Signal_MXs1_nominal_MXs3_channel_MXs3_scaleFactors");
		MXs1_MXs3_scale->addTerm(mu1);

		w->Print();
	}


	return 0;
}