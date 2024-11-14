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
		RooStats::ModelConfig* temp_Model = (RooStats::ModelConfig*)sbModel->Clone(("BonlyModel_" + std::to_string(i)).c_str());
		w.factory("PROD::Signal_MXs1_nominal_MXs1_channel_MXs1_scaleFactors(Signal_MXs1_nominal_MXs1_channel_MXs1_epsilon, FBDT_efficiency_CAL_MXs1)");
		printf("print!\n");
		w->Print();
	}


	return 0;
}