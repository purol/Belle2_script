#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include <cmath>
#include <float.h>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "RooGenericPdf.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "TH1.h"
#include "TF1.h"
#include "RooDataHist.h"
using namespace RooFit;
using std::string;
using std::to_string;

void Hypocal()
{
	TFile* f = TFile::Open("CLs_test.root");
	RooWorkspace* w = (RooWorkspace*)f->Get("w");
	w->Print();
	RooAbsData* data = w->data("observed_data");

	if (data->isWeighted()) {
		Info("Hypocal", "Data set is weighted, nentries = %d and sum of weights = %8.1f.\n", data->numEntries(), data->sumEntries());
	}

	RooStats::ModelConfig* sbModel = (RooStats::ModelConfig*)w->obj("ModelConfig");
	RooStats::ModelConfig* bModel = (RooStats::ModelConfig*)sbModel->Clone("BonlyModel");
	RooRealVar* poi = (RooRealVar*)bModel->GetParametersOfInterest()->first();
	poi->setVal(0);
	bModel->SetSnapshot(*poi);

	RooStats::HybridCalculator hybCalc(*data, *bModel, *sbModel);
	RooStats::ProfileLikelihoodTestStat* plr = new RooStats::ProfileLikelihoodTestStat(*sbModel->GetPdf());
	plr->SetOneSided(true);

	RooStats::ToyMCSampler* toymcs = (RooStats::ToyMCSampler*)hybCalc.GetTestStatSampler();
	toymcs->SetTestStatistic(plr);
	hybCalc.SetToys(1000, 1000);

	RooStats::HypoTestInverter inverter(hybCalc);
	inverter.SetConfidenceLevel(0.90);
	inverter.UseCLs(true);
	inverter.SetVerbose(false);
	inverter.SetFixedScan(300, 0.0, 0.0002); // set number of points , xmin and xmax

	RooStats::HypoTestInverterResult* result = inverter.GetInterval();

	cout << 100 * inverter.ConfidenceLevel() << "%  upper limit : " << result->UpperLimit() << endl;

	std::cout << "Expected upper limits, using the B (alternate) model : " << std::endl;
	std::cout << " expected limit (median) " << result->GetExpectedUpperLimit(0) << std::endl;
	std::cout << " expected limit (-1 sig) " << result->GetExpectedUpperLimit(-1) << std::endl;
	std::cout << " expected limit (+1 sig) " << result->GetExpectedUpperLimit(1) << std::endl;
	std::cout << " expected limit (-2 sig) " << result->GetExpectedUpperLimit(-2) << std::endl;
	std::cout << " expected limit (+2 sig) " << result->GetExpectedUpperLimit(2) << std::endl;

	TCanvas* c1 = new TCanvas();
	RooStats::HypoTestInverterPlot* plot = new RooStats::HypoTestInverterPlot("HTI_Result_Plot", "HypoTest Scan Result", result);
	plot->Draw("CLb 2CL");  // plot also CLb and CLs+b 
	c1->Draw(); c1->SaveAs("CLs_hyb.png");
}
