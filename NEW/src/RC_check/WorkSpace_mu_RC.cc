#include <iostream>
#include <stdio.h>
#include <string>
#include <cstring>

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
#include "TFile.h"

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

const bool KnunuOnly = false;

void ModifyPiecewiseInterpolationForAll(RooWorkspace* ws, int code) {
	RooArgSet funcs = ws->allFunctions();
	TIter it = funcs.createIterator();
	TObject* tempObj = nullptr;
	while ((tempObj = it.Next())) {
		PiecewiseInterpolation* piec = dynamic_cast<PiecewiseInterpolation*>(tempObj);
		if (piec) {
			piec->setAllInterpCodes(code);
		}
	}
}

void ModifyPiecewiseInterpolation(RooWorkspace* ws, const char* name, int code) {
	RooArgSet funcs = ws->allFunctions();
	TIter it = funcs.createIterator();
	TObject* tempObj = nullptr;
	while ((tempObj = it.Next())) {
		PiecewiseInterpolation* piec = dynamic_cast<PiecewiseInterpolation*>(tempObj);
		if (piec) {
			RooRealVar* alpha = ws->var(name);
			piec->setInterpCode(*alpha, code);
		}
	}
}

void ModifyFlexibleInterpolationForAll(RooWorkspace* ws, int code) {
	RooArgSet funcs = ws->allFunctions();
	TIter it = funcs.createIterator();
	TObject* tempObj = nullptr;
	while ((tempObj = it.Next())) {
		FlexibleInterpVar* flex = dynamic_cast<FlexibleInterpVar*>(tempObj);
		if (flex) {
			flex->setAllInterpCodes(code);
		}
	}
}

void ModifyFlexibleInterpolation(RooWorkspace* ws, const char* name, int code) {
	RooArgSet funcs = ws->allFunctions();
	TIter it = funcs.createIterator();
	TObject* tempObj = nullptr;
	while ((tempObj = it.Next())) {
		FlexibleInterpVar* flex = dynamic_cast<FlexibleInterpVar*>(tempObj);
		if (flex) {
			RooRealVar* alpha = ws->var(name);
			flex->setInterpCode(*alpha, code);
		}
	}
}

void CheckInterpolation(RooWorkspace* ws) {
	RooArgSet funcs = ws->allFunctions();
	TIter it = funcs.createIterator();
	TObject* tempObj = 0;
	while ((tempObj = it.Next())) {
		FlexibleInterpVar* flex = dynamic_cast<FlexibleInterpVar*>(tempObj);
		if (flex) {
			flex->printAllInterpCodes();
		}
		PiecewiseInterpolation* piec = dynamic_cast<PiecewiseInterpolation*>(tempObj);
		if (piec) {
			piec->printAllInterpCodes();
		}
	}
}

int ReadNFEIEigenVector(const char* dirname) {
	int Nentry = 0; // number of eigen values/vectors
	double eigen_value = 0; // eigen value
	double weight_sys[RarityBins * 5] = { 0.0 }; // eigen vector

	FILE* fp;
	fp = fopen(dirname, "r");
	while (true) {
		if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
		for (int i = 0; i < RarityBins * 5; i++) {
			if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
		}
		Nentry++;
	}
	fclose(fp);

	return Nentry;
}

int ReadNPIDEigenVector(const char* dirname) {
	int Nentry = 0; // number of eigen values/vectors
	double eigen_value = 0; // eigen value
	double weight_sys[RarityBins * 9] = { 0.0 }; // eigen vector

	FILE* fp;
	fp = fopen(dirname, "r");
	while (true) {
		if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
		for (int i = 0; i < RarityBins * 9; i++) {
			if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
		}
		Nentry++;
	}
	fclose(fp);

	return Nentry;
}

int ReadNBREigenVector(const char* dirname) {
	int Nentry = 0; // number of eigen values/vectors
	double eigen_value = 0; // eigen value
	double weight_sys[RarityBins * 5] = { 0.0 }; // eigen vector

	FILE* fp;
	fp = fopen(dirname, "r");
	while (true) {
		if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
		for (int i = 0; i < RarityBins * 5; i++) {
			if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
		}
		Nentry++;
	}
	fclose(fp);

	return Nentry;
}

int ReadNpi0EigenVector(const char* dirname) {
	int Nentry = 0; // number of eigen values/vectors
	double eigen_value = 0; // eigen value
	double weight_sys[RarityBins * 9] = { 0.0 }; // eigen vector

	FILE* fp;
	fp = fopen(dirname, "r");
	while (true) {
		if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
		for (int i = 0; i < RarityBins * 9; i++) {
			if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
		}
		Nentry++;
	}
	fclose(fp);

	return Nentry;
}

int ReadMultiplicityInfo(const char* dirname) {
	int Nentry = 0; // number of eigen values/vectors
	double eigen_value = 0; // eigen value
	double weight_sys[RarityBins * 9] = { 0.0 }; // eigen vector

	FILE* fp;
	fp = fopen(dirname, "r");
	while (true) {
		if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
		for (int i = 0; i < RarityBins * 9; i++) {
			if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
		}
		Nentry++;
	}
	fclose(fp);

	return Nentry;
}

int ReadNFragmentationEigenVector(const char* dirname) {
	int Nentry = 0; // number of eigen values/vectors
	double eigen_value = 0; // eigen value
	double weight_sys[RarityBins * 5] = { 0.0 }; // eigen vector

	FILE* fp;
	fp = fopen(dirname, "r");
	while (true) {
		if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
		for (int i = 0; i < RarityBins * 5; i++) {
			if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
		}
		Nentry++;
	}
	fclose(fp);

	return Nentry;
}

bool IsThereAnyChange(const char* fname, const char* nominal_hist_name, const char* hist_m_name, const char* hist_p_name) {
	TFile* file = new TFile(fname);
	file->cd();

	TH1D* nominal_hist = (TH1D*)file->Get(nominal_hist_name);
	TH1D* hist_m = (TH1D*)file->Get(hist_m_name);
	TH1D* hist_p = (TH1D*)file->Get(hist_p_name);

	bool ThereIsChange = false;
	for (int i = 0; i < RarityBins; i++) {
		double nominal_value = nominal_hist->GetBinContent(i + 1);
		double m_value = hist_m->GetBinContent(i + 1);
		double p_value = hist_p->GetBinContent(i + 1);

		if (std::abs(m_value - nominal_value) < MyEPSILON) ThereIsChange = false;
		else {
			ThereIsChange = true;
			break;
		}

		if (std::abs(p_value - nominal_value) < MyEPSILON) ThereIsChange = false;
		else {
			ThereIsChange = true;
			break;
		}
	}

	file->Close();

	if (ThereIsChange == false) {
		printf("[IsThereAnyChange] There is not difference among %s, %s, and %s in %s\n", nominal_hist_name, hist_m_name, hist_p_name, fname);
	}

	return ThereIsChange;
}

void AddEmptySample(HistFactory::Channel* channel, int MXs_bin) {
	// MXs_bin: reco bin

	std::string bin_name = "";
	if (MXs_bin == 1) bin_name = std::string("MXs1");
	else if (MXs_bin == 2) bin_name = std::string("MXs2");
	else if (MXs_bin == 3) bin_name = std::string("MXs3");
	else {
		printf("inappropritate bin\n");
		exit(1);
	}

	// make empty root file
	const char* empty_file_name = "empty_file.root";
	TFile* input_file = new TFile(empty_file_name, "RECREATE");
	TH1D* temp_hist;
	if (MXs_bin == 1) {
		temp_hist = new TH1D("empty_hist", "empty_hist", RarityBins_MX1, 0, RarityBins_MX1);
		for (int i = 0; i < RarityBins_MX1; i++) {
			temp_hist->SetBinError(i + 1, 0);
			temp_hist->SetBinContent(i + 1, 0);
		}
	}
	else if (MXs_bin == 2) {
		temp_hist = new TH1D("empty_hist", "empty_hist", RarityBins_MX2, 0, RarityBins_MX2);
		for (int i = 0; i < RarityBins_MX2; i++) {
			temp_hist->SetBinError(i + 1, 0);
			temp_hist->SetBinContent(i + 1, 0);
		}
	}
	else if (MXs_bin == 3) {
		temp_hist = new TH1D("empty_hist", "empty_hist", RarityBins_MX3, 0, RarityBins_MX3);
		for (int i = 0; i < RarityBins_MX3; i++) {
			temp_hist->SetBinError(i + 1, 0);
			temp_hist->SetBinContent(i + 1, 0);
		}
	}
	else {
		printf("inappropritate bin\n");
		exit(1);
	}
	temp_hist->Write();
	input_file->Close();

	// point data
	channel->SetData("empty_hist", empty_file_name);

	// read signal template
	/* ================================ SIGNAL with true MXs1 ================================ */
	RooStats::HistFactory::Sample sig_temp_MXs1(("Signal_MXs1_nominal_" + bin_name).c_str(), "empty_hist", empty_file_name);
	channel->AddSample(sig_temp_MXs1);

	/* ================================ SIGNAL with true MXs2 ================================ */
	RooStats::HistFactory::Sample sig_temp_MXs2(("Signal_MXs2_nominal_" + bin_name).c_str(), "empty_hist", empty_file_name);
	channel->AddSample(sig_temp_MXs2);

	/* ================================ SIGNAL with true MXs3 ================================ */
	RooStats::HistFactory::Sample sig_temp_MXs3(("Signal_MXs3_nominal_" + bin_name).c_str(), "empty_hist", empty_file_name);
	channel->AddSample(sig_temp_MXs3);

	// read background template
	/* ================================ CHG ================================ */
	RooStats::HistFactory::Sample CHG_temp(("CHG_nominal_" + bin_name).c_str(), "empty_hist", empty_file_name);
	channel->AddSample(CHG_temp);
	/* ================================ CHG ================================ */

	/* ================================ MIX ================================ */
	RooStats::HistFactory::Sample MIX_temp(("MIX_nominal_" + bin_name).c_str(), "empty_hist", empty_file_name);
	channel->AddSample(MIX_temp);
	/* ================================ MIX ================================ */

	/* ================================ UUBAR ================================ */
	RooStats::HistFactory::Sample UUBAR_temp(("UUBAR_nominal_" + bin_name).c_str(), "empty_hist", empty_file_name);
	channel->AddSample(UUBAR_temp);
	/* ================================ UUBAR ================================ */

	/* ================================ DDBAR ================================ */
	RooStats::HistFactory::Sample DDBAR_temp(("DDBAR_nominal_" + bin_name).c_str(), "empty_hist", empty_file_name);
	channel->AddSample(DDBAR_temp);
	/* ================================ DDBAR ================================ */

	/* ================================ SSBAR ================================ */
	RooStats::HistFactory::Sample SSBAR_temp(("SSBAR_nominal_" + bin_name).c_str(), "empty_hist", empty_file_name);
	channel->AddSample(SSBAR_temp);
	/* ================================ SSBAR ================================ */

	/* ================================ CHARM ================================ */
	RooStats::HistFactory::Sample CHARM_temp(("CHARM_nominal_" + bin_name).c_str(), "empty_hist", empty_file_name);
	channel->AddSample(CHARM_temp);
	/* ================================ CHARM ================================ */
}

void AddSample(HistFactory::Channel* channel, const char* fname, int MXs_bin, const double expmu, double qq_CAL, double qq_CAL_relativeuncer, double FBDT_CAL, double FBDT_CAL_relativeuncer, double bkg_norm_relativeuncer, const char* version) {
	// MXs_bin: reco bin

	int NEntryFEI = ReadNFEIEigenVector(("./KumoiRD_" + std::string(version) + "/set1/FEI_selected.txt").c_str());
	int NEntryKID = ReadNPIDEigenVector(("./KumoiRD_" + std::string(version) + "/set1/KID_selected.txt").c_str());
	int NEntryPID = ReadNPIDEigenVector(("./KumoiRD_" + std::string(version) + "/set1/PID_selected.txt").c_str());
	int NEntryBR = ReadNBREigenVector(("./KumoiRD_" + std::string(version) + "/set1/BR_selected.txt").c_str());
	int NEntrypi0 = ReadNpi0EigenVector(("./KumoiRD_" + std::string(version) + "/set1/pi0_selected.txt").c_str());
	int NEntryMultiplicity = ReadMultiplicityInfo(("./KumoiRD_" + std::string(version) + "/set1/multiplicity_selected.txt").c_str());
	int NEntryFragmentation = ReadNFragmentationEigenVector(("./KumoiRD_" + std::string(version) + "/set1/Fragmentation_selected.txt").c_str());

	std::string bin_name = "";
	if (MXs_bin == 1) bin_name = std::string("MXs1");
	else if (MXs_bin == 2) bin_name = std::string("MXs2");
	else if (MXs_bin == 3) bin_name = std::string("MXs3");
	else {
		printf("inappropritate bin\n");
		exit(1);
	}

	// point data
	channel->SetData("total_DATA", fname);

	// read signal template
	/* ================================ SIGNAL with true MXs1 ================================ */
	RooStats::HistFactory::Sample sig_temp_MXs1(("Signal_MXs1_nominal_" + bin_name).c_str(), "Signal_MXs1_nominal", fname);
	if (IsThereAnyChange(fname, "Signal_MXs1_nominal", "Signal_MXs1_track_m", "Signal_MXs1_track_p")) sig_temp_MXs1.AddHistoSys("track_eff_uncer", "Signal_MXs1_track_m", fname, "", "Signal_MXs1_track_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs1_nominal", "Signal_MXs1_KS0_m", "Signal_MXs1_KS0_p")) sig_temp_MXs1.AddHistoSys("KS0_reco_uncer", "Signal_MXs1_KS0_m", fname, "", "Signal_MXs1_KS0_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs1_nominal", "Signal_MXs1_BtoDtoXKL_m", "Signal_MXs1_BtoDtoXKL_p")) sig_temp_MXs1.AddHistoSys("BtoDtoXKL_uncer", "Signal_MXs1_BtoDtoXKL_m", fname, "", "Signal_MXs1_BtoDtoXKL_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs1_nominal", "Signal_MXs1_BRBtoXKLKL_m", "Signal_MXs1_BRBtoXKLKL_p")) sig_temp_MXs1.AddHistoSys("BRBtoXKLKL_uncer", "Signal_MXs1_BRBtoXKLKL_m", fname, "", "Signal_MXs1_BRBtoXKLKL_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs1_nominal", "Signal_MXs1_EffECLKL_m", "Signal_MXs1_EffECLKL_p")) sig_temp_MXs1.AddHistoSys("EffECLKL_uncer", "Signal_MXs1_EffECLKL_m", fname, "", "Signal_MXs1_EffECLKL_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs1_nominal", "Signal_MXs1_NEWFEICAL_m", "Signal_MXs1_NEWFEICAL_p")) sig_temp_MXs1.AddHistoSys("NEWFEICAL_uncer", "Signal_MXs1_NEWFEICAL_m", fname, "", "Signal_MXs1_NEWFEICAL_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs1_nominal", "Signal_MXs1_Xnn_m", "Signal_MXs1_Xnn_p")) sig_temp_MXs1.AddHistoSys("Xnn_BR_uncer", "Signal_MXs1_Xnn_m", fname, "", "Signal_MXs1_Xnn_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) if (IsThereAnyChange(fname, "Signal_MXs1_nominal", ("Signal_MXs1_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), ("Signal_MXs1_multiplicity_correlated" + std::to_string(i) + "_p").c_str())) sig_temp_MXs1.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs1_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs1_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryFEI; i++) if (IsThereAnyChange(fname, "Signal_MXs1_nominal", ("Signal_MXs1_FEI_correlated" + std::to_string(i) + "_m").c_str(), ("Signal_MXs1_FEI_correlated" + std::to_string(i) + "_p").c_str())) sig_temp_MXs1.AddHistoSys(("FEI" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs1_FEI_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs1_FEI_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) if (IsThereAnyChange(fname, "Signal_MXs1_nominal", ("Signal_MXs1_KID_correlated" + std::to_string(i) + "_m").c_str(), ("Signal_MXs1_KID_correlated" + std::to_string(i) + "_p").c_str())) sig_temp_MXs1.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs1_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs1_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntryPID; i++) if (IsThereAnyChange(fname, "Signal_MXs1_nominal", ("Signal_MXs1_PID_correlated" + std::to_string(i) + "_m").c_str(), ("Signal_MXs1_PID_correlated" + std::to_string(i) + "_p").c_str())) sig_temp_MXs1.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs1_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs1_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryBR; i++) if (IsThereAnyChange(fname, "Signal_MXs1_nominal", ("Signal_MXs1_BR_correlated" + std::to_string(i) + "_m").c_str(), ("Signal_MXs1_BR_correlated" + std::to_string(i) + "_p").c_str())) sig_temp_MXs1.AddHistoSys(("BBBR" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs1_BR_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs1_BR_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntrypi0; i++) if (IsThereAnyChange(fname, "Signal_MXs1_nominal", ("Signal_MXs1_pi0_correlated" + std::to_string(i) + "_m").c_str(), ("Signal_MXs1_pi0_correlated" + std::to_string(i) + "_p").c_str())) sig_temp_MXs1.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs1_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs1_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs1_nominal", "Signal_MXs1_Kff1_m", "Signal_MXs1_Kff1_p")) sig_temp_MXs1.AddHistoSys("Kff1_uncer", "Signal_MXs1_Kff1_m", fname, "", "Signal_MXs1_Kff1_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs1_nominal", "Signal_MXs1_Kff2_m", "Signal_MXs1_Kff2_p")) sig_temp_MXs1.AddHistoSys("Kff2_uncer", "Signal_MXs1_Kff2_m", fname, "", "Signal_MXs1_Kff2_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs1_nominal", "Signal_MXs1_Kff3_m", "Signal_MXs1_Kff3_p")) sig_temp_MXs1.AddHistoSys("Kff3_uncer", "Signal_MXs1_Kff3_m", fname, "", "Signal_MXs1_Kff3_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs1_nominal", "Signal_MXs1_Kfrac_m", "Signal_MXs1_Kfrac_p")) sig_temp_MXs1.AddHistoSys("Kfrac_uncer", "Signal_MXs1_Kfrac_m", fname, "", "Signal_MXs1_Kfrac_p", fname, "");
	if (MXs_bin == 1) sig_temp_MXs1.AddNormFactor(("FBDT_efficiency_CAL_" + bin_name).c_str(), FBDT_CAL, FBDT_CAL, FBDT_CAL);
	if (MXs_bin == 1) sig_temp_MXs1.AddOverallSys(("FBDT_efficiency_uncer_" + bin_name).c_str(), 1.0 - FBDT_CAL_relativeuncer, 1.0 + FBDT_CAL_relativeuncer);
	sig_temp_MXs1.AddOverallSys("BB_counting_uncer", 0.9855, 1.0145);
	//sig_temp_MXs1.AddOverallSys("mu1_fitter_bias_mean", 0.8615, 1.1385);
	//sig_temp_MXs1.AddOverallSys("mu1_fitter_bias_sigma", 0.9937, 1.0063);
	sig_temp_MXs1.ActivateStatError("Signal_MXs1_all_uncorrelated_MC_stat", fname, "");
	sig_temp_MXs1.SetNormalizeByTheory(kFALSE);
	sig_temp_MXs1.AddNormFactor("mu_MXs1", expmu, -100.0, 100.0);
	channel->AddSample(sig_temp_MXs1);

	/* ================================ SIGNAL with true MXs2 ================================ */
	RooStats::HistFactory::Sample sig_temp_MXs2(("Signal_MXs2_nominal_" + bin_name).c_str(), "Signal_MXs2_nominal", fname);
	if (IsThereAnyChange(fname, "Signal_MXs2_nominal", "Signal_MXs2_track_m", "Signal_MXs2_track_p")) sig_temp_MXs2.AddHistoSys("track_eff_uncer", "Signal_MXs2_track_m", fname, "", "Signal_MXs2_track_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs2_nominal", "Signal_MXs2_KS0_m", "Signal_MXs2_KS0_p")) sig_temp_MXs2.AddHistoSys("KS0_reco_uncer", "Signal_MXs2_KS0_m", fname, "", "Signal_MXs2_KS0_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs2_nominal", "Signal_MXs2_BtoDtoXKL_m", "Signal_MXs2_BtoDtoXKL_p")) sig_temp_MXs2.AddHistoSys("BtoDtoXKL_uncer", "Signal_MXs2_BtoDtoXKL_m", fname, "", "Signal_MXs2_BtoDtoXKL_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs2_nominal", "Signal_MXs2_BRBtoXKLKL_m", "Signal_MXs2_BRBtoXKLKL_p")) sig_temp_MXs2.AddHistoSys("BRBtoXKLKL_uncer", "Signal_MXs2_BRBtoXKLKL_m", fname, "", "Signal_MXs2_BRBtoXKLKL_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs2_nominal", "Signal_MXs2_EffECLKL_m", "Signal_MXs2_EffECLKL_p")) sig_temp_MXs2.AddHistoSys("EffECLKL_uncer", "Signal_MXs2_EffECLKL_m", fname, "", "Signal_MXs2_EffECLKL_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs2_nominal", "Signal_MXs2_NEWFEICAL_m", "Signal_MXs2_NEWFEICAL_p")) sig_temp_MXs2.AddHistoSys("NEWFEICAL_uncer", "Signal_MXs2_NEWFEICAL_m", fname, "", "Signal_MXs2_NEWFEICAL_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs2_nominal", "Signal_MXs2_Xnn_m", "Signal_MXs2_Xnn_p")) sig_temp_MXs2.AddHistoSys("Xnn_BR_uncer", "Signal_MXs2_Xnn_m", fname, "", "Signal_MXs2_Xnn_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) if (IsThereAnyChange(fname, "Signal_MXs2_nominal", ("Signal_MXs2_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), ("Signal_MXs2_multiplicity_correlated" + std::to_string(i) + "_p").c_str())) sig_temp_MXs2.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs2_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs2_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryFEI; i++) if (IsThereAnyChange(fname, "Signal_MXs2_nominal", ("Signal_MXs2_FEI_correlated" + std::to_string(i) + "_m").c_str(), ("Signal_MXs2_FEI_correlated" + std::to_string(i) + "_p").c_str())) sig_temp_MXs2.AddHistoSys(("FEI" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs2_FEI_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs2_FEI_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) if (IsThereAnyChange(fname, "Signal_MXs2_nominal", ("Signal_MXs2_KID_correlated" + std::to_string(i) + "_m").c_str(), ("Signal_MXs2_KID_correlated" + std::to_string(i) + "_p").c_str())) sig_temp_MXs2.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs2_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs2_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntryPID; i++) if (IsThereAnyChange(fname, "Signal_MXs2_nominal", ("Signal_MXs2_PID_correlated" + std::to_string(i) + "_m").c_str(), ("Signal_MXs2_PID_correlated" + std::to_string(i) + "_p").c_str())) sig_temp_MXs2.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs2_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs2_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryBR; i++) if (IsThereAnyChange(fname, "Signal_MXs2_nominal", ("Signal_MXs2_BR_correlated" + std::to_string(i) + "_m").c_str(), ("Signal_MXs2_BR_correlated" + std::to_string(i) + "_p").c_str())) sig_temp_MXs2.AddHistoSys(("BBBR" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs2_BR_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs2_BR_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntrypi0; i++) if (IsThereAnyChange(fname, "Signal_MXs2_nominal", ("Signal_MXs2_pi0_correlated" + std::to_string(i) + "_m").c_str(), ("Signal_MXs2_pi0_correlated" + std::to_string(i) + "_p").c_str())) sig_temp_MXs2.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs2_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs2_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs2_nominal", "Signal_MXs2_Kstarff1_m", "Signal_MXs2_Kstarff1_p")) sig_temp_MXs2.AddHistoSys("Kstarff1_uncer", "Signal_MXs2_Kstarff1_m", fname, "", "Signal_MXs2_Kstarff1_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs2_nominal", "Signal_MXs2_Kstarff2_m", "Signal_MXs2_Kstarff2_p")) sig_temp_MXs2.AddHistoSys("Kstarff2_uncer", "Signal_MXs2_Kstarff2_m", fname, "", "Signal_MXs2_Kstarff2_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs2_nominal", "Signal_MXs2_Kstarff3_m", "Signal_MXs2_Kstarff3_p")) sig_temp_MXs2.AddHistoSys("Kstarff3_uncer", "Signal_MXs2_Kstarff3_m", fname, "", "Signal_MXs2_Kstarff3_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs2_nominal", "Signal_MXs2_Kstarff4_m", "Signal_MXs2_Kstarff4_p")) sig_temp_MXs2.AddHistoSys("Kstarff4_uncer", "Signal_MXs2_Kstarff4_m", fname, "", "Signal_MXs2_Kstarff4_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs2_nominal", "Signal_MXs2_Kstarff5_m", "Signal_MXs2_Kstarff5_p")) sig_temp_MXs2.AddHistoSys("Kstarff5_uncer", "Signal_MXs2_Kstarff5_m", fname, "", "Signal_MXs2_Kstarff5_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs2_nominal", "Signal_MXs2_Kstarff6_m", "Signal_MXs2_Kstarff6_p")) sig_temp_MXs2.AddHistoSys("Kstarff6_uncer", "Signal_MXs2_Kstarff6_m", fname, "", "Signal_MXs2_Kstarff6_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs2_nominal", "Signal_MXs2_Kstarff7_m", "Signal_MXs2_Kstarff7_p")) sig_temp_MXs2.AddHistoSys("Kstarff7_uncer", "Signal_MXs2_Kstarff7_m", fname, "", "Signal_MXs2_Kstarff7_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs2_nominal", "Signal_MXs2_Kstarff8_m", "Signal_MXs2_Kstarff8_p")) sig_temp_MXs2.AddHistoSys("Kstarff8_uncer", "Signal_MXs2_Kstarff8_m", fname, "", "Signal_MXs2_Kstarff8_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs2_nominal", "Signal_MXs2_Kstarff9_m", "Signal_MXs2_Kstarff9_p")) sig_temp_MXs2.AddHistoSys("Kstarff9_uncer", "Signal_MXs2_Kstarff9_m", fname, "", "Signal_MXs2_Kstarff9_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs2_nominal", "Signal_MXs2_Kstarfrac_m", "Signal_MXs2_Kstarfrac_p")) sig_temp_MXs2.AddHistoSys("Kstarfrac_uncer", "Signal_MXs2_Kstarfrac_m", fname, "", "Signal_MXs2_Kstarfrac_p", fname, "");
	// sig_temp_MXs2.AddHistoSys("mKstar_uncer", "Signal_MXs2_mKstar_m", fname, "", "Signal_MXs2_mKstar_p", fname, "");
	if (MXs_bin == 2) sig_temp_MXs2.AddNormFactor(("FBDT_efficiency_CAL_" + bin_name).c_str(), FBDT_CAL, FBDT_CAL, FBDT_CAL);
	if (MXs_bin == 2) sig_temp_MXs2.AddOverallSys(("FBDT_efficiency_uncer_" + bin_name).c_str(), 1.0 - FBDT_CAL_relativeuncer, 1.0 + FBDT_CAL_relativeuncer);
	sig_temp_MXs2.AddOverallSys("BB_counting_uncer", 0.9855, 1.0145);
	//sig_temp_MXs2.AddOverallSys("mu2_fitter_bias_mean", 0.9715, 1.0284);
	//sig_temp_MXs2.AddOverallSys("mu2_fitter_bias_sigma", 0.9830, 1.0170);
	sig_temp_MXs2.ActivateStatError("Signal_MXs2_all_uncorrelated_MC_stat", fname, "");
	sig_temp_MXs2.SetNormalizeByTheory(kFALSE);
	if (KnunuOnly) sig_temp_MXs2.AddNormFactor("mu_MXs2", expmu, 1.0, 1.0);
	else sig_temp_MXs2.AddNormFactor("mu_MXs2", expmu, -100.0, 100.0);
	channel->AddSample(sig_temp_MXs2);

	/* ================================ SIGNAL with true MXs3 ================================ */
	RooStats::HistFactory::Sample sig_temp_MXs3(("Signal_MXs3_nominal_" + bin_name).c_str(), "Signal_MXs3_nominal", fname);
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_track_m", "Signal_MXs3_track_p")) sig_temp_MXs3.AddHistoSys("track_eff_uncer", "Signal_MXs3_track_m", fname, "", "Signal_MXs3_track_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_KS0_m", "Signal_MXs3_KS0_p")) sig_temp_MXs3.AddHistoSys("KS0_reco_uncer", "Signal_MXs3_KS0_m", fname, "", "Signal_MXs3_KS0_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_BtoDtoXKL_m", "Signal_MXs3_BtoDtoXKL_p")) sig_temp_MXs3.AddHistoSys("BtoDtoXKL_uncer", "Signal_MXs3_BtoDtoXKL_m", fname, "", "Signal_MXs3_BtoDtoXKL_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_BRBtoXKLKL_m", "Signal_MXs3_BRBtoXKLKL_p")) sig_temp_MXs3.AddHistoSys("BRBtoXKLKL_uncer", "Signal_MXs3_BRBtoXKLKL_m", fname, "", "Signal_MXs3_BRBtoXKLKL_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_EffECLKL_m", "Signal_MXs3_EffECLKL_p")) sig_temp_MXs3.AddHistoSys("EffECLKL_uncer", "Signal_MXs3_EffECLKL_m", fname, "", "Signal_MXs3_EffECLKL_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_NEWFEICAL_m", "Signal_MXs3_NEWFEICAL_p")) sig_temp_MXs3.AddHistoSys("NEWFEICAL_uncer", "Signal_MXs3_NEWFEICAL_m", fname, "", "Signal_MXs3_NEWFEICAL_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_Xnn_m", "Signal_MXs3_Xnn_p")) sig_temp_MXs3.AddHistoSys("Xnn_BR_uncer", "Signal_MXs3_Xnn_m", fname, "", "Signal_MXs3_Xnn_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) if (IsThereAnyChange(fname, "Signal_MXs3_nominal", ("Signal_MXs3_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), ("Signal_MXs3_multiplicity_correlated" + std::to_string(i) + "_p").c_str())) sig_temp_MXs3.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs3_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs3_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryFEI; i++) if (IsThereAnyChange(fname, "Signal_MXs3_nominal", ("Signal_MXs3_FEI_correlated" + std::to_string(i) + "_m").c_str(), ("Signal_MXs3_FEI_correlated" + std::to_string(i) + "_p").c_str())) sig_temp_MXs3.AddHistoSys(("FEI" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs3_FEI_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs3_FEI_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) if (IsThereAnyChange(fname, "Signal_MXs3_nominal", ("Signal_MXs3_KID_correlated" + std::to_string(i) + "_m").c_str(), ("Signal_MXs3_KID_correlated" + std::to_string(i) + "_p").c_str())) sig_temp_MXs3.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs3_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs3_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntryPID; i++) if (IsThereAnyChange(fname, "Signal_MXs3_nominal", ("Signal_MXs3_PID_correlated" + std::to_string(i) + "_m").c_str(), ("Signal_MXs3_PID_correlated" + std::to_string(i) + "_p").c_str())) sig_temp_MXs3.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs3_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs3_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryBR; i++) if (IsThereAnyChange(fname, "Signal_MXs3_nominal", ("Signal_MXs3_BR_correlated" + std::to_string(i) + "_m").c_str(), ("Signal_MXs3_BR_correlated" + std::to_string(i) + "_p").c_str())) sig_temp_MXs3.AddHistoSys(("BBBR" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs3_BR_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs3_BR_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntrypi0; i++) if (IsThereAnyChange(fname, "Signal_MXs3_nominal", ("Signal_MXs3_pi0_correlated" + std::to_string(i) + "_m").c_str(), ("Signal_MXs3_pi0_correlated" + std::to_string(i) + "_p").c_str())) sig_temp_MXs3.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs3_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs3_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_Kff1_m", "Signal_MXs3_Kff1_p")) sig_temp_MXs3.AddHistoSys("Kff1_uncer", "Signal_MXs3_Kff1_m", fname, "", "Signal_MXs3_Kff1_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_Kff2_m", "Signal_MXs3_Kff2_p")) sig_temp_MXs3.AddHistoSys("Kff2_uncer", "Signal_MXs3_Kff2_m", fname, "", "Signal_MXs3_Kff2_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_Kff3_m", "Signal_MXs3_Kff3_p")) sig_temp_MXs3.AddHistoSys("Kff3_uncer", "Signal_MXs3_Kff3_m", fname, "", "Signal_MXs3_Kff3_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_Kstarff1_m", "Signal_MXs3_Kstarff1_p")) sig_temp_MXs3.AddHistoSys("Kstarff1_uncer", "Signal_MXs3_Kstarff1_m", fname, "", "Signal_MXs3_Kstarff1_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_Kstarff2_m", "Signal_MXs3_Kstarff2_p")) sig_temp_MXs3.AddHistoSys("Kstarff2_uncer", "Signal_MXs3_Kstarff2_m", fname, "", "Signal_MXs3_Kstarff2_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_Kstarff3_m", "Signal_MXs3_Kstarff3_p")) sig_temp_MXs3.AddHistoSys("Kstarff3_uncer", "Signal_MXs3_Kstarff3_m", fname, "", "Signal_MXs3_Kstarff3_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_Kstarff4_m", "Signal_MXs3_Kstarff4_p")) sig_temp_MXs3.AddHistoSys("Kstarff4_uncer", "Signal_MXs3_Kstarff4_m", fname, "", "Signal_MXs3_Kstarff4_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_Kstarff5_m", "Signal_MXs3_Kstarff5_p")) sig_temp_MXs3.AddHistoSys("Kstarff5_uncer", "Signal_MXs3_Kstarff5_m", fname, "", "Signal_MXs3_Kstarff5_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_Kstarff6_m", "Signal_MXs3_Kstarff6_p")) sig_temp_MXs3.AddHistoSys("Kstarff6_uncer", "Signal_MXs3_Kstarff6_m", fname, "", "Signal_MXs3_Kstarff6_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_Kstarff7_m", "Signal_MXs3_Kstarff7_p")) sig_temp_MXs3.AddHistoSys("Kstarff7_uncer", "Signal_MXs3_Kstarff7_m", fname, "", "Signal_MXs3_Kstarff7_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_Kstarff8_m", "Signal_MXs3_Kstarff8_p")) sig_temp_MXs3.AddHistoSys("Kstarff8_uncer", "Signal_MXs3_Kstarff8_m", fname, "", "Signal_MXs3_Kstarff8_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_Kstarff9_m", "Signal_MXs3_Kstarff9_p")) sig_temp_MXs3.AddHistoSys("Kstarff9_uncer", "Signal_MXs3_Kstarff9_m", fname, "", "Signal_MXs3_Kstarff9_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_Kfrac_m", "Signal_MXs3_Kfrac_p")) sig_temp_MXs3.AddHistoSys("Kfrac_uncer", "Signal_MXs3_Kfrac_m", fname, "", "Signal_MXs3_Kfrac_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_Kstarfrac_m", "Signal_MXs3_Kstarfrac_p")) sig_temp_MXs3.AddHistoSys("Kstarfrac_uncer", "Signal_MXs3_Kstarfrac_m", fname, "", "Signal_MXs3_Kstarfrac_p", fname, "");
	for (int i = 0; i < NEntryFragmentation; i++) if (IsThereAnyChange(fname, "Signal_MXs3_nominal", ("Signal_MXs3_Fragmentation_correlated" + std::to_string(i) + "_m").c_str(), ("Signal_MXs3_Fragmentation_correlated" + std::to_string(i) + "_p").c_str())) sig_temp_MXs3.AddHistoSys(("Xs_fragmentation" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs3_Fragmentation_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs3_Fragmentation_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_pf_m", "Signal_MXs3_pf_p")) sig_temp_MXs3.AddHistoSys("pf_uncer", "Signal_MXs3_pf_m", fname, "", "Signal_MXs3_pf_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_mb_m", "Signal_MXs3_mb_p")) sig_temp_MXs3.AddHistoSys("mb_uncer", "Signal_MXs3_mb_m", fname, "", "Signal_MXs3_mb_p", fname, "");
	if (IsThereAnyChange(fname, "Signal_MXs3_nominal", "Signal_MXs3_transition_m", "Signal_MXs3_transition_p")) {
		// read shape change
		sig_temp_MXs3.AddHistoSys("transition_uncer", "Signal_MXs3_transition_m", fname, "", "Signal_MXs3_transition_p", fname, "");

		// read overall norm change
		double scale_p = 1.0;
		double scale_m = 1.0;
		FILE* fp;

		if (MXs_bin == 1) fp = fopen(("./KumoiRD_" + std::string(version) + "/set1/Scale_Signal_MXs3_transition_p_MXs1.txt").c_str(), "r");
		else if (MXs_bin == 2) fp = fopen(("./KumoiRD_" + std::string(version) + "/set1/Scale_Signal_MXs3_transition_p_MXs2.txt").c_str(), "r");
		else if (MXs_bin == 3) fp = fopen(("./KumoiRD_" + std::string(version) + "/set1/Scale_Signal_MXs3_transition_p_MXs3.txt").c_str(), "r");
		fscanf(fp, "%lf", &scale_p);
		fclose(fp);

		if (MXs_bin == 1) fp = fopen(("./KumoiRD_" + std::string(version) + "/set1/Scale_Signal_MXs3_transition_m_MXs1.txt").c_str(), "r");
		else if (MXs_bin == 2) fp = fopen(("./KumoiRD_" + std::string(version) + "/set1/Scale_Signal_MXs3_transition_m_MXs2.txt").c_str(), "r");
		else if (MXs_bin == 3) fp = fopen(("./KumoiRD_" + std::string(version) + "/set1/Scale_Signal_MXs3_transition_m_MXs3.txt").c_str(), "r");
		fscanf(fp, "%lf", &scale_m);
		fclose(fp);

		if (scale_p < MyEPSILON) scale_p = 1.0;
		if (scale_m < MyEPSILON) scale_m = 1.0;

		sig_temp_MXs3.AddOverallSys("transition_uncer", 1.0 / scale_m, 1.0 / scale_p);

		printf("OverallSys %lf and %lf are assigned\n", scale_m, scale_p);
	}
	// sig_temp_MXs3.AddHistoSys("mKstar_uncer", "Signal_MXs3_mKstar_m", fname, "", "Signal_MXs3_mKstar_p", fname, "");
	if (MXs_bin == 3) sig_temp_MXs3.AddNormFactor(("FBDT_efficiency_CAL_" + bin_name).c_str(), FBDT_CAL, FBDT_CAL, FBDT_CAL);
	if (MXs_bin == 3) sig_temp_MXs3.AddOverallSys(("FBDT_efficiency_uncer_" + bin_name).c_str(), 1.0 - FBDT_CAL_relativeuncer, 1.0 + FBDT_CAL_relativeuncer);
	sig_temp_MXs3.AddOverallSys("BB_counting_uncer", 0.9855, 1.0145);
	//sig_temp_MXs3.AddOverallSys("mu3_fitter_bias_mean", 0.9197, 1.0803);
	//sig_temp_MXs3.AddOverallSys("mu3_fitter_bias_sigma", 0.9757, 1.0243);
	sig_temp_MXs3.ActivateStatError("Signal_MXs3_all_uncorrelated_MC_stat", fname, "");
	sig_temp_MXs3.SetNormalizeByTheory(kFALSE);
	if (KnunuOnly) sig_temp_MXs3.AddNormFactor("mu_MXs3", expmu, 1.0, 1.0);
	else sig_temp_MXs3.AddNormFactor("mu_MXs3", expmu, -100.0, 100.0);
	channel->AddSample(sig_temp_MXs3);

	// read background template
	/* ================================ CHG ================================ */
	RooStats::HistFactory::Sample CHG_temp(("CHG_nominal_" + bin_name).c_str(), "CHG_nominal", fname);
	if (IsThereAnyChange(fname, "CHG_nominal", "CHG_track_m", "CHG_track_p")) CHG_temp.AddHistoSys("track_eff_uncer", "CHG_track_m", fname, "", "CHG_track_p", fname, "");
	if (IsThereAnyChange(fname, "CHG_nominal", "CHG_KS0_m", "CHG_KS0_p")) CHG_temp.AddHistoSys("KS0_reco_uncer", "CHG_KS0_m", fname, "", "CHG_KS0_p", fname, "");
	if (IsThereAnyChange(fname, "CHG_nominal", "CHG_BtoDtoXKL_m", "CHG_BtoDtoXKL_p")) CHG_temp.AddHistoSys("BtoDtoXKL_uncer", "CHG_BtoDtoXKL_m", fname, "", "CHG_BtoDtoXKL_p", fname, "");
	if (IsThereAnyChange(fname, "CHG_nominal", "CHG_BRBtoXKLKL_m", "CHG_BRBtoXKLKL_p")) CHG_temp.AddHistoSys("BRBtoXKLKL_uncer", "CHG_BRBtoXKLKL_m", fname, "", "CHG_BRBtoXKLKL_p", fname, "");
	if (IsThereAnyChange(fname, "CHG_nominal", "CHG_EffECLKL_m", "CHG_EffECLKL_p")) CHG_temp.AddHistoSys("EffECLKL_uncer", "CHG_EffECLKL_m", fname, "", "CHG_EffECLKL_p", fname, "");
	if (IsThereAnyChange(fname, "CHG_nominal", "CHG_Xnn_m", "CHG_Xnn_p")) CHG_temp.AddHistoSys("Xnn_BR_uncer", "CHG_Xnn_m", fname, "", "CHG_Xnn_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) if (IsThereAnyChange(fname, "CHG_nominal", ("CHG_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), ("CHG_multiplicity_correlated" + std::to_string(i) + "_p").c_str())) CHG_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("CHG_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryFEI; i++) if (IsThereAnyChange(fname, "CHG_nominal", ("CHG_FEI_correlated" + std::to_string(i) + "_m").c_str(), ("CHG_FEI_correlated" + std::to_string(i) + "_p").c_str())) CHG_temp.AddHistoSys(("FEI" + std::to_string(i) + "_uncer").c_str(), ("CHG_FEI_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_FEI_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) if (IsThereAnyChange(fname, "CHG_nominal", ("CHG_KID_correlated" + std::to_string(i) + "_m").c_str(), ("CHG_KID_correlated" + std::to_string(i) + "_p").c_str())) CHG_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("CHG_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntryPID; i++) if (IsThereAnyChange(fname, "CHG_nominal", ("CHG_PID_correlated" + std::to_string(i) + "_m").c_str(), ("CHG_PID_correlated" + std::to_string(i) + "_p").c_str())) CHG_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("CHG_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryBR; i++) if (IsThereAnyChange(fname, "CHG_nominal", ("CHG_BR_correlated" + std::to_string(i) + "_m").c_str(), ("CHG_BR_correlated" + std::to_string(i) + "_p").c_str())) CHG_temp.AddHistoSys(("BBBR" + std::to_string(i) + "_uncer").c_str(), ("CHG_BR_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_BR_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntrypi0; i++) if (IsThereAnyChange(fname, "CHG_nominal", ("CHG_pi0_correlated" + std::to_string(i) + "_m").c_str(), ("CHG_pi0_correlated" + std::to_string(i) + "_p").c_str())) CHG_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("CHG_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	//CHG_temp.AddHistoSys("BDTc_shape_BB", "CHG_BDTc_m", fname, "", "CHG_BDTc_p", fname, "");
	CHG_temp.AddOverallSys(("mu_CHG_" + bin_name).c_str(), 1.0 - bkg_norm_relativeuncer, 1.0 + bkg_norm_relativeuncer);
	CHG_temp.AddOverallSys("BB_counting_uncer", 0.9855, 1.0145);
	CHG_temp.ActivateStatError("CHG_all_uncorrelated_MC_stat", fname, "");
	//if (MXs_bin == 1) CHG_temp.AddShapeSys(("CHG_additional_uncorr_uncer" + bin_name).c_str(), RooStats::HistFactory::Constraint::Type::Gaussian, "CHG_rel_uncer", fname, "");
	CHG_temp.SetNormalizeByTheory(kFALSE);
	channel->AddSample(CHG_temp);
	/* ================================ CHG ================================ */

	/* ================================ MIX ================================ */
	RooStats::HistFactory::Sample MIX_temp(("MIX_nominal_" + bin_name).c_str(), "MIX_nominal", fname);
	if (IsThereAnyChange(fname, "MIX_nominal", "MIX_track_m", "MIX_track_p")) MIX_temp.AddHistoSys("track_eff_uncer", "MIX_track_m", fname, "", "MIX_track_p", fname, "");
	if (IsThereAnyChange(fname, "MIX_nominal", "MIX_KS0_m", "MIX_KS0_p")) MIX_temp.AddHistoSys("KS0_reco_uncer", "MIX_KS0_m", fname, "", "MIX_KS0_p", fname, "");
	if (IsThereAnyChange(fname, "MIX_nominal", "MIX_BtoDtoXKL_m", "MIX_BtoDtoXKL_p")) MIX_temp.AddHistoSys("BtoDtoXKL_uncer", "MIX_BtoDtoXKL_m", fname, "", "MIX_BtoDtoXKL_p", fname, "");
	if (IsThereAnyChange(fname, "MIX_nominal", "MIX_BRBtoXKLKL_m", "MIX_BRBtoXKLKL_p")) MIX_temp.AddHistoSys("BRBtoXKLKL_uncer", "MIX_BRBtoXKLKL_m", fname, "", "MIX_BRBtoXKLKL_p", fname, "");
	if (IsThereAnyChange(fname, "MIX_nominal", "MIX_EffECLKL_m", "MIX_EffECLKL_p")) MIX_temp.AddHistoSys("EffECLKL_uncer", "MIX_EffECLKL_m", fname, "", "MIX_EffECLKL_p", fname, "");
	if (IsThereAnyChange(fname, "MIX_nominal", "MIX_Xnn_m", "MIX_Xnn_p")) MIX_temp.AddHistoSys("Xnn_BR_uncer", "MIX_Xnn_m", fname, "", "MIX_Xnn_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) if (IsThereAnyChange(fname, "MIX_nominal", ("MIX_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), ("MIX_multiplicity_correlated" + std::to_string(i) + "_p").c_str())) MIX_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("MIX_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryFEI; i++) if (IsThereAnyChange(fname, "MIX_nominal", ("MIX_FEI_correlated" + std::to_string(i) + "_m").c_str(), ("MIX_FEI_correlated" + std::to_string(i) + "_p").c_str())) MIX_temp.AddHistoSys(("FEI" + std::to_string(i) + "_uncer").c_str(), ("MIX_FEI_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_FEI_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) if (IsThereAnyChange(fname, "MIX_nominal", ("MIX_KID_correlated" + std::to_string(i) + "_m").c_str(), ("MIX_KID_correlated" + std::to_string(i) + "_p").c_str())) MIX_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("MIX_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntryPID; i++) if (IsThereAnyChange(fname, "MIX_nominal", ("MIX_PID_correlated" + std::to_string(i) + "_m").c_str(), ("MIX_PID_correlated" + std::to_string(i) + "_p").c_str())) MIX_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("MIX_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryBR; i++) if (IsThereAnyChange(fname, "MIX_nominal", ("MIX_BR_correlated" + std::to_string(i) + "_m").c_str(), ("MIX_BR_correlated" + std::to_string(i) + "_p").c_str())) MIX_temp.AddHistoSys(("BBBR" + std::to_string(i) + "_uncer").c_str(), ("MIX_BR_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_BR_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntrypi0; i++) if (IsThereAnyChange(fname, "MIX_nominal", ("MIX_pi0_correlated" + std::to_string(i) + "_m").c_str(), ("MIX_pi0_correlated" + std::to_string(i) + "_p").c_str())) MIX_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("MIX_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	//MIX_temp.AddHistoSys("BDTc_shape_BB", "MIX_BDTc_m", fname, "", "MIX_BDTc_p", fname, "");
	MIX_temp.AddOverallSys(("mu_MIX_" + bin_name).c_str(), 1.0 - bkg_norm_relativeuncer, 1.0 + bkg_norm_relativeuncer);
	MIX_temp.AddOverallSys("BB_counting_uncer", 0.9855, 1.0145);
	MIX_temp.ActivateStatError("MIX_all_uncorrelated_MC_stat", fname, "");
	//if (MXs_bin == 1) MIX_temp.AddShapeSys(("MIX_additional_uncorr_uncer" + bin_name).c_str(), RooStats::HistFactory::Constraint::Type::Gaussian, "MIX_rel_uncer", fname, "");
	MIX_temp.SetNormalizeByTheory(kFALSE);
	channel->AddSample(MIX_temp);
	/* ================================ MIX ================================ */

	/* ================================ UUBAR ================================ */
	RooStats::HistFactory::Sample UUBAR_temp(("UUBAR_nominal_" + bin_name).c_str(), "UUBAR_nominal", fname);
	if (IsThereAnyChange(fname, "UUBAR_nominal", "UUBAR_track_m", "UUBAR_track_p")) UUBAR_temp.AddHistoSys("track_eff_uncer", "UUBAR_track_m", fname, "", "UUBAR_track_p", fname, "");
	if (IsThereAnyChange(fname, "UUBAR_nominal", "UUBAR_KS0_m", "UUBAR_KS0_p")) UUBAR_temp.AddHistoSys("KS0_reco_uncer", "UUBAR_KS0_m", fname, "", "UUBAR_KS0_p", fname, "");
	if (IsThereAnyChange(fname, "UUBAR_nominal", "UUBAR_EffECLKL_m", "UUBAR_EffECLKL_p")) UUBAR_temp.AddHistoSys("EffECLKL_uncer", "UUBAR_EffECLKL_m", fname, "", "UUBAR_EffECLKL_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) if (IsThereAnyChange(fname, "UUBAR_nominal", ("UUBAR_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), ("UUBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str())) UUBAR_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("UUBAR_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("UUBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) if (IsThereAnyChange(fname, "UUBAR_nominal", ("UUBAR_KID_correlated" + std::to_string(i) + "_m").c_str(), ("UUBAR_KID_correlated" + std::to_string(i) + "_p").c_str())) UUBAR_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("UUBAR_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("UUBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntryPID; i++) if (IsThereAnyChange(fname, "UUBAR_nominal", ("UUBAR_PID_correlated" + std::to_string(i) + "_m").c_str(), ("UUBAR_PID_correlated" + std::to_string(i) + "_p").c_str())) UUBAR_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("UUBAR_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("UUBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntrypi0; i++) if (IsThereAnyChange(fname, "UUBAR_nominal", ("UUBAR_pi0_correlated" + std::to_string(i) + "_m").c_str(), ("UUBAR_pi0_correlated" + std::to_string(i) + "_p").c_str())) UUBAR_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("UUBAR_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("UUBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (IsThereAnyChange(fname, "UUBAR_nominal", "UUBAR_BDTc_m", "UUBAR_BDTc_p")) UUBAR_temp.AddHistoSys("BDTc_shape_qq", "UUBAR_BDTc_m", fname, "", "UUBAR_BDTc_p", fname, "");
	UUBAR_temp.AddNormFactor(("qq_CAL_UUBAR_" + bin_name).c_str(), qq_CAL, qq_CAL, qq_CAL);
	UUBAR_temp.AddOverallSys(("qq_CAL_UUBAR_uncer_" + bin_name).c_str(), 1.0 - qq_CAL_relativeuncer, 1.0 + qq_CAL_relativeuncer);
	UUBAR_temp.AddOverallSys(("mu_UUBAR_" + bin_name).c_str(), 1.0 - bkg_norm_relativeuncer, 1.0 + bkg_norm_relativeuncer);
	UUBAR_temp.ActivateStatError("UUBAR_all_uncorrelated_MC_stat", fname, "");
	//if (MXs_bin == 1) UUBAR_temp.AddShapeSys(("UUBAR_additional_uncorr_uncer" + bin_name).c_str(), RooStats::HistFactory::Constraint::Type::Gaussian, "UUBAR_rel_uncer", fname, "");
	UUBAR_temp.SetNormalizeByTheory(kFALSE);
	channel->AddSample(UUBAR_temp);
	/* ================================ UUBAR ================================ */

	/* ================================ DDBAR ================================ */
	RooStats::HistFactory::Sample DDBAR_temp(("DDBAR_nominal_" + bin_name).c_str(), "DDBAR_nominal", fname);
	if (IsThereAnyChange(fname, "DDBAR_nominal", "DDBAR_track_m", "DDBAR_track_p")) DDBAR_temp.AddHistoSys("track_eff_uncer", "DDBAR_track_m", fname, "", "DDBAR_track_p", fname, "");
	if (IsThereAnyChange(fname, "DDBAR_nominal", "DDBAR_KS0_m", "DDBAR_KS0_p")) DDBAR_temp.AddHistoSys("KS0_reco_uncer", "DDBAR_KS0_m", fname, "", "DDBAR_KS0_p", fname, "");
	if (IsThereAnyChange(fname, "DDBAR_nominal", "DDBAR_EffECLKL_m", "DDBAR_EffECLKL_p")) DDBAR_temp.AddHistoSys("EffECLKL_uncer", "DDBAR_EffECLKL_m", fname, "", "DDBAR_EffECLKL_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) if (IsThereAnyChange(fname, "DDBAR_nominal", ("DDBAR_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), ("DDBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str())) DDBAR_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("DDBAR_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("DDBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) if (IsThereAnyChange(fname, "DDBAR_nominal", ("DDBAR_KID_correlated" + std::to_string(i) + "_m").c_str(), ("DDBAR_KID_correlated" + std::to_string(i) + "_p").c_str())) DDBAR_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("DDBAR_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("DDBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntryPID; i++) if (IsThereAnyChange(fname, "DDBAR_nominal", ("DDBAR_PID_correlated" + std::to_string(i) + "_m").c_str(), ("DDBAR_PID_correlated" + std::to_string(i) + "_p").c_str())) DDBAR_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("DDBAR_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("DDBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntrypi0; i++) if (IsThereAnyChange(fname, "DDBAR_nominal", ("DDBAR_pi0_correlated" + std::to_string(i) + "_m").c_str(), ("DDBAR_pi0_correlated" + std::to_string(i) + "_p").c_str())) DDBAR_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("DDBAR_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("DDBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (IsThereAnyChange(fname, "DDBAR_nominal", "DDBAR_BDTc_m", "DDBAR_BDTc_p")) DDBAR_temp.AddHistoSys("BDTc_shape_qq", "DDBAR_BDTc_m", fname, "", "DDBAR_BDTc_p", fname, "");
	DDBAR_temp.AddNormFactor(("qq_CAL_DDBAR_" + bin_name).c_str(), qq_CAL, qq_CAL, qq_CAL);
	DDBAR_temp.AddOverallSys(("qq_CAL_DDBAR_uncer_" + bin_name).c_str(), 1.0 - qq_CAL_relativeuncer, 1.0 + qq_CAL_relativeuncer);
	DDBAR_temp.AddOverallSys(("mu_DDBAR_" + bin_name).c_str(), 1.0 - bkg_norm_relativeuncer, 1.0 + bkg_norm_relativeuncer);
	DDBAR_temp.ActivateStatError("DDBAR_all_uncorrelated_MC_stat", fname, "");
	//if (MXs_bin == 1) DDBAR_temp.AddShapeSys(("DDBAR_additional_uncorr_uncer" + bin_name).c_str(), RooStats::HistFactory::Constraint::Type::Gaussian, "DDBAR_rel_uncer", fname, "");
	DDBAR_temp.SetNormalizeByTheory(kFALSE);
	channel->AddSample(DDBAR_temp);
	/* ================================ DDBAR ================================ */

	/* ================================ SSBAR ================================ */
	RooStats::HistFactory::Sample SSBAR_temp(("SSBAR_nominal_" + bin_name).c_str(), "SSBAR_nominal", fname);
	if (IsThereAnyChange(fname, "SSBAR_nominal", "SSBAR_track_m", "SSBAR_track_p")) SSBAR_temp.AddHistoSys("track_eff_uncer", "SSBAR_track_m", fname, "", "SSBAR_track_p", fname, "");
	if (IsThereAnyChange(fname, "SSBAR_nominal", "SSBAR_KS0_m", "SSBAR_KS0_p")) SSBAR_temp.AddHistoSys("KS0_reco_uncer", "SSBAR_KS0_m", fname, "", "SSBAR_KS0_p", fname, "");
	if (IsThereAnyChange(fname, "SSBAR_nominal", "SSBAR_EffECLKL_m", "SSBAR_EffECLKL_p")) SSBAR_temp.AddHistoSys("EffECLKL_uncer", "SSBAR_EffECLKL_m", fname, "", "SSBAR_EffECLKL_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) if (IsThereAnyChange(fname, "SSBAR_nominal", ("SSBAR_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), ("SSBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str())) SSBAR_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("SSBAR_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("SSBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) if (IsThereAnyChange(fname, "SSBAR_nominal", ("SSBAR_KID_correlated" + std::to_string(i) + "_m").c_str(), ("SSBAR_KID_correlated" + std::to_string(i) + "_p").c_str())) SSBAR_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("SSBAR_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("SSBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntryPID; i++) if (IsThereAnyChange(fname, "SSBAR_nominal", ("SSBAR_PID_correlated" + std::to_string(i) + "_m").c_str(), ("SSBAR_PID_correlated" + std::to_string(i) + "_p").c_str())) SSBAR_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("SSBAR_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("SSBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntrypi0; i++) if (IsThereAnyChange(fname, "SSBAR_nominal", ("SSBAR_pi0_correlated" + std::to_string(i) + "_m").c_str(), ("SSBAR_pi0_correlated" + std::to_string(i) + "_p").c_str())) SSBAR_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("SSBAR_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("SSBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (IsThereAnyChange(fname, "SSBAR_nominal", "SSBAR_BDTc_m", "SSBAR_BDTc_p")) SSBAR_temp.AddHistoSys("BDTc_shape_qq", "SSBAR_BDTc_m", fname, "", "SSBAR_BDTc_p", fname, "");
	SSBAR_temp.AddNormFactor(("qq_CAL_SSBAR_" + bin_name).c_str(), qq_CAL, qq_CAL, qq_CAL);
	SSBAR_temp.AddOverallSys(("qq_CAL_SSBAR_uncer_" + bin_name).c_str(), 1.0 - qq_CAL_relativeuncer, 1.0 + qq_CAL_relativeuncer);
	SSBAR_temp.AddOverallSys(("mu_SSBAR_" + bin_name).c_str(), 1.0 - bkg_norm_relativeuncer, 1.0 + bkg_norm_relativeuncer);
	SSBAR_temp.ActivateStatError("SSBAR_all_uncorrelated_MC_stat", fname, "");
	//if (MXs_bin == 1) SSBAR_temp.AddShapeSys(("SSBAR_additional_uncorr_uncer" + bin_name).c_str(), RooStats::HistFactory::Constraint::Type::Gaussian, "SSBAR_rel_uncer", fname, "");
	SSBAR_temp.SetNormalizeByTheory(kFALSE);
	channel->AddSample(SSBAR_temp);
	/* ================================ SSBAR ================================ */

	/* ================================ CHARM ================================ */
	RooStats::HistFactory::Sample CHARM_temp(("CHARM_nominal_" + bin_name).c_str(), "CHARM_nominal", fname);
	if (IsThereAnyChange(fname, "CHARM_nominal", "CHARM_track_m", "CHARM_track_p")) CHARM_temp.AddHistoSys("track_eff_uncer", "CHARM_track_m", fname, "", "CHARM_track_p", fname, "");
	if (IsThereAnyChange(fname, "CHARM_nominal", "CHARM_KS0_m", "CHARM_KS0_p")) CHARM_temp.AddHistoSys("KS0_reco_uncer", "CHARM_KS0_m", fname, "", "CHARM_KS0_p", fname, "");
	if (IsThereAnyChange(fname, "CHARM_nominal", "CHARM_EffECLKL_m", "CHARM_EffECLKL_p")) CHARM_temp.AddHistoSys("EffECLKL_uncer", "CHARM_EffECLKL_m", fname, "", "CHARM_EffECLKL_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) if (IsThereAnyChange(fname, "CHARM_nominal", ("CHARM_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), ("CHARM_multiplicity_correlated" + std::to_string(i) + "_p").c_str())) CHARM_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("CHARM_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHARM_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) if (IsThereAnyChange(fname, "CHARM_nominal", ("CHARM_KID_correlated" + std::to_string(i) + "_m").c_str(), ("CHARM_KID_correlated" + std::to_string(i) + "_p").c_str())) CHARM_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("CHARM_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHARM_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntryPID; i++) if (IsThereAnyChange(fname, "CHARM_nominal", ("CHARM_PID_correlated" + std::to_string(i) + "_m").c_str(), ("CHARM_PID_correlated" + std::to_string(i) + "_p").c_str())) CHARM_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("CHARM_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHARM_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntrypi0; i++) if (IsThereAnyChange(fname, "CHARM_nominal", ("CHARM_pi0_correlated" + std::to_string(i) + "_m").c_str(), ("CHARM_pi0_correlated" + std::to_string(i) + "_p").c_str())) CHARM_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("CHARM_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHARM_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (IsThereAnyChange(fname, "CHARM_nominal", "CHARM_BDTc_m", "CHARM_BDTc_p")) CHARM_temp.AddHistoSys("BDTc_shape_qq", "CHARM_BDTc_m", fname, "", "CHARM_BDTc_p", fname, "");
	CHARM_temp.AddNormFactor(("qq_CAL_CHARM_" + bin_name).c_str(), qq_CAL, qq_CAL, qq_CAL);
	CHARM_temp.AddOverallSys(("qq_CAL_CHARM_uncer_" + bin_name).c_str(), 1.0 - qq_CAL_relativeuncer, 1.0 + qq_CAL_relativeuncer);
	CHARM_temp.AddOverallSys(("mu_CHARM_" + bin_name).c_str(), 1.0 - bkg_norm_relativeuncer, 1.0 + bkg_norm_relativeuncer);
	CHARM_temp.ActivateStatError("CHARM_all_uncorrelated_MC_stat", fname, "");
	//if (MXs_bin == 1) CHARM_temp.AddShapeSys(("CHARM_additional_uncorr_uncer" + bin_name).c_str(), RooStats::HistFactory::Constraint::Type::Gaussian, "CHARM_rel_uncer", fname, "");
	CHARM_temp.SetNormalizeByTheory(kFALSE);
	channel->AddSample(CHARM_temp);
	/* ================================ CHARM ================================ */
}

int main(int argc, char* argv[]) {
	/*
	* argv[1]: version (ex. v020)
	* argv[2]: output path
	*/

	// read CAL
	double CAL_off_K = 1.0;
	double CAL_off_Kstar = 1.0;
	double CAL_off_Xs = 1.0;

	double CAL_off_K_uncer = 1.0;
	double CAL_off_Kstar_uncer = 1.0;
	double CAL_off_Xs_uncer = 1.0;

	FILE* fp;
	fp = fopen(("/home/belle2/junewoo/storage_b1/bsub/FastBDTinputdrawtest/KumoiRD_off_RC_check/" + std::string(argv[1]) + "/MC_data_ratio.txt").c_str(), "r");
	fscanf(fp, "%lf %lf\n", &CAL_off_K, &CAL_off_K_uncer);
	fscanf(fp, "%lf %lf\n", &CAL_off_Kstar, &CAL_off_Kstar_uncer);
	fscanf(fp, "%lf %lf\n", &CAL_off_Xs, &CAL_off_Xs_uncer);
	fclose(fp);

	double CAL_FBDT_K = 1.0;
	double CAL_FBDT_Kstar = 1.0;
	double CAL_FBDT_Xs = 1.0;

	double CAL_FBDT_K_uncer = 1.0;
	double CAL_FBDT_Kstar_uncer = 1.0;
	double CAL_FBDT_Xs_uncer = 1.0;

	fp = fopen(("/home/belle2/junewoo/storage_b1/bsub/FastBDTinputdrawtest/KumoiRD_Jpsi_RC_check/" + std::string(argv[1]) + "/MC_data_ratio.txt").c_str(), "r");
	fscanf(fp, "%lf %lf\n", &CAL_FBDT_K, &CAL_FBDT_K_uncer);
	fscanf(fp, "%lf %lf\n", &CAL_FBDT_Kstar, &CAL_FBDT_Kstar_uncer);
	fscanf(fp, "%lf %lf\n", &CAL_FBDT_Xs, &CAL_FBDT_Xs_uncer);
	fclose(fp);


	const double expmu = 1.0;
	const char* fname = strdup(("./KumoiRD_" + std::string(argv[1]) + "/set1/PDFandDATA_nominal.root").c_str());
	const char* fname_MXs1 = strdup(("./KumoiRD_" + std::string(argv[1]) + "/set1/PDFandDATA_MXs1.root").c_str());
	const char* fname_MXs2 = strdup(("./KumoiRD_" + std::string(argv[1]) + "/set1/PDFandDATA_MXs2.root").c_str());
	const char* fname_MXs3 = strdup(("./KumoiRD_" + std::string(argv[1]) + "/set1/PDFandDATA_MXs3.root").c_str());

	int NEntryFragmentation = ReadNFragmentationEigenVector(("./KumoiRD_" + std::string(argv[1]) + "/set1/Fragmentation_selected.txt").c_str());

	RooStats::HistFactory::Measurement meas("my_measurement", "my measurement");
	meas.SetOutputFilePrefix((std::string(argv[2]) + "/PDFandDATA/my_measurement").c_str());
	meas.SetExportOnly(kTRUE);

	// setting measurement
	meas.SetPOI("mu");
	meas.SetLumi(1.0);
	meas.AddConstantParam("Lumi");

	// set constant parameters
	meas.AddConstantParam("FBDT_efficiency_CAL_MXs1");
	meas.AddConstantParam("FBDT_efficiency_CAL_MXs2");
	meas.AddConstantParam("FBDT_efficiency_CAL_MXs3");

	meas.AddConstantParam("qq_CAL_UUBAR_MXs1");
	meas.AddConstantParam("qq_CAL_UUBAR_MXs2");
	meas.AddConstantParam("qq_CAL_UUBAR_MXs3");

	meas.AddConstantParam("qq_CAL_DDBAR_MXs1");
	meas.AddConstantParam("qq_CAL_DDBAR_MXs2");
	meas.AddConstantParam("qq_CAL_DDBAR_MXs3");

	meas.AddConstantParam("qq_CAL_SSBAR_MXs1");
	meas.AddConstantParam("qq_CAL_SSBAR_MXs2");
	meas.AddConstantParam("qq_CAL_SSBAR_MXs3");

	meas.AddConstantParam("qq_CAL_CHARM_MXs1");
	meas.AddConstantParam("qq_CAL_CHARM_MXs2");
	meas.AddConstantParam("qq_CAL_CHARM_MXs3");

	// define channels
	HistFactory::Channel channel_MXs1("channel_MXs1");
	channel_MXs1.SetStatErrorConfig(1e-5, "Gaussian");

	HistFactory::Channel channel_MXs2("channel_MXs2");
	channel_MXs2.SetStatErrorConfig(1e-5, "Gaussian");

	HistFactory::Channel channel_MXs3("channel_MXs3");
	channel_MXs3.SetStatErrorConfig(1e-5, "Gaussian");

	// get MC
	AddSample(&channel_MXs1, fname_MXs1, 1, expmu, CAL_off_K, CAL_off_K_uncer / CAL_off_K, CAL_FBDT_K, CAL_FBDT_K_uncer / CAL_FBDT_K, 0.2, argv[1]);
	if (KnunuOnly) {
		AddEmptySample(&channel_MXs2, 2);
		AddEmptySample(&channel_MXs3, 3);
	}
	else {
		AddSample(&channel_MXs2, fname_MXs2, 2, expmu, CAL_off_Kstar, CAL_off_Kstar_uncer / CAL_off_Kstar, CAL_FBDT_Kstar, CAL_FBDT_Kstar_uncer / CAL_FBDT_Kstar, 0.2, argv[1]);
		AddSample(&channel_MXs3, fname_MXs3, 3, expmu, CAL_off_Xs, CAL_off_Xs_uncer / CAL_off_Xs, CAL_FBDT_Xs, CAL_FBDT_Xs_uncer / CAL_FBDT_Xs, 0.2, argv[1]);
	}

	// add channel to measurement
	meas.AddChannel(channel_MXs1);
	meas.AddChannel(channel_MXs2);
	meas.AddChannel(channel_MXs3);
	meas.CollectHistograms();

	RooWorkspace* w;
	w = RooStats::HistFactory::MakeModelAndMeasurementFast(meas);

	w->Print();
	CheckInterpolation(w);
	w->writeToFile((std::string(argv[2]) + "/PDFandDATA_workspace_mu.root").c_str());

	// meas.PrintXML("PDFandDATA_MXs1");

	return 0;
}