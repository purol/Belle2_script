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

int template_draw() {

	const char* fname = "PDFandDATA.root";

	TFile* f = new TFile(fname, "read");


	TH1D* signal_hist = (TH1D*) f->Get("Signal_nominal");
	TH1D* chg_hist = (TH1D*)f->Get("CHG_nominal");
	TH1D* mix_hist = (TH1D*)f->Get("MIX_nominal");
	TH1D* uubar_hist = (TH1D*)f->Get("UUBAR_nominal");
	TH1D* ddbar_hist = (TH1D*)f->Get("DDBAR_nominal");
	TH1D* ssbar_hist = (TH1D*)f->Get("SSBAR_nominal");
	TH1D* charm_hist = (TH1D*)f->Get("CHARM_nominal");

	TCanvas* c_temp = new TCanvas("c", "", 800, 800);
	signal_hist->SetStats(false);
	signal_hist->GetYaxis()->SetRangeUser(0., 2.3);
	signal_hist->SetFillColor(kBlue + 1);
	signal_hist->SetLineWidth(0);
	signal_hist->GetXaxis()->SetTitle("bin index");
	signal_hist->GetYaxis()->SetTitle("number of events");
	signal_hist->Draw("hist");
	c_temp->SaveAs("signal_hist.png");

	THStack* Stack = new THStack("thstack", ";bin index;number of events");

	chg_hist->SetFillColor(kRed + 1);
	chg_hist->SetLineWidth(0);

	mix_hist->SetFillColor(kOrange + 1);
	mix_hist->SetLineWidth(0);

	uubar_hist->SetFillColor(kGreen + 1);
	uubar_hist->SetLineWidth(0);

	ddbar_hist->SetFillColor(kAzure + 1);
	ddbar_hist->SetLineWidth(0);

	ssbar_hist->SetFillColor(kViolet + 1);
	ssbar_hist->SetLineWidth(0);

	charm_hist->SetFillColor(kPink + 1);
	charm_hist->SetLineWidth(0);

	Stack->Add(chg_hist);
	Stack->Add(mix_hist);
	Stack->Add(uubar_hist);
	Stack->Add(ddbar_hist);
	Stack->Add(ssbar_hist);
	Stack->Add(charm_hist);

	Stack->Draw("hist");

	TLegend* legend = gPad->BuildLegend(0.9, 0.9, 0.7, 0.7);
	legend->SetFillStyle(0);
	legend->SetLineWidth(0);

	c_temp->SaveAs("bkgs_hist.png");

	return 0;
}