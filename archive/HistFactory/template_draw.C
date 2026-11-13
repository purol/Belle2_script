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
#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TAxis.h"
#include "TStyle.h"

#include "RooChi2Var.h"
#include "RooAbsData.h"
#include "RooRealSumPdf.h"
#include "RooPoisson.h"
#include "RooGaussian.h"
#include "RooRealVar.h"
#include "RooMCStudy.h"
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

	const char* fname = "PDFandDATA_nominal.root";

	TFile* f = new TFile(fname, "read");


	TH1D* signal_hist = (TH1D*) f->Get("Signal_nominal");
	TH1D* chg_hist = (TH1D*)f->Get("CHG_nominal");
	TH1D* mix_hist = (TH1D*)f->Get("MIX_nominal");
	TH1D* uubar_hist = (TH1D*)f->Get("UUBAR_nominal");
	TH1D* ddbar_hist = (TH1D*)f->Get("DDBAR_nominal");
	TH1D* ssbar_hist = (TH1D*)f->Get("SSBAR_nominal");
	TH1D* charm_hist = (TH1D*)f->Get("CHARM_nominal");

	TCanvas* c_temp = new TCanvas("c", "", 800, 800);
	c_temp->SetCanvasSize(800, 800);
	c_temp->SetLeftMargin(0.17);
	c_temp->SetRightMargin(0.04);
	c_temp->SetBottomMargin(0.16);
	c_temp->SetTopMargin(0.06);
	c_temp->cd();
	gStyle->SetOptStat(0);

	// Both panels are included at half the text width in the thesis.
	const auto style_axes = [](TAxis* x_axis, TAxis* y_axis) {
		x_axis->SetTitle("Bin index");
		y_axis->SetTitle("Events");
		x_axis->SetLabelFont(42);
		y_axis->SetLabelFont(42);
		x_axis->SetTitleFont(42);
		y_axis->SetTitleFont(42);
		x_axis->SetLabelSize(0.040);
		y_axis->SetLabelSize(0.045);
		x_axis->SetTitleSize(0.055);
		y_axis->SetTitleSize(0.055);
		x_axis->SetLabelOffset(0.010);
		y_axis->SetLabelOffset(0.010);
		x_axis->SetTitleOffset(1.20);
		y_axis->SetTitleOffset(1.35);
		x_axis->CenterTitle();
		y_axis->CenterTitle();
		y_axis->SetNdivisions(505);
	};

	signal_hist->SetStats(false);
	signal_hist->SetTitle("");
	signal_hist->GetYaxis()->SetRangeUser(0., 5.0);
	signal_hist->SetFillColor(kBlue + 1);
	signal_hist->SetLineWidth(0);
	style_axes(signal_hist->GetXaxis(), signal_hist->GetYaxis());
	signal_hist->Draw("hist");
	c_temp->RedrawAxis();
	c_temp->Modified();
	c_temp->Update();
	c_temp->SaveAs("signal_hist.png");
	c_temp->Print("signal_hist.pdf", "pdf Portrait");

	c_temp->Clear();
	c_temp->cd();
	THStack* Stack = new THStack("thstack", ";Bin index;Events");

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

	// Reserve space above the distribution for a two-row legend.
	Stack->SetMinimum(0.0);
	Stack->SetMaximum(1.35 * Stack->GetMaximum());
	Stack->Draw("hist");
	style_axes(Stack->GetXaxis(), Stack->GetYaxis());
	c_temp->Modified();
	c_temp->Update();

	TLegend* legend = new TLegend(0.22, 0.79, 0.94, 0.93);
	legend->SetNColumns(3);
	legend->SetTextFont(42);
	legend->SetTextSize(0.045);
	legend->SetMargin(0.18);
	legend->SetFillStyle(0);
	legend->SetBorderSize(0);
	legend->AddEntry(chg_hist, "charged", "f");
	legend->AddEntry(mix_hist, "mixed", "f");
	legend->AddEntry(uubar_hist, "u#bar{u}", "f");
	legend->AddEntry(ddbar_hist, "d#bar{d}", "f");
	legend->AddEntry(ssbar_hist, "s#bar{s}", "f");
	legend->AddEntry(charm_hist, "c#bar{c}", "f");
	legend->Draw();

	c_temp->RedrawAxis();
	c_temp->Modified();
	c_temp->Update();
	c_temp->SaveAs("bkgs_hist.png");
	c_temp->Print("bkgs_hist.pdf", "pdf Portrait");

	return 0;
}
