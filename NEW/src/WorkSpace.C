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

#include "template.h"
#include "correctors.h"

using namespace RooFit;
using namespace RooStats;
using namespace HistFactory;

using std::string;
using std::to_string;

Corrector_Fragmentation corrector_Fragmentation;

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

void AddSample(HistFactory::Channel* channel, const char* fname, int MXs_bin, const double expmu, double qq_CAL, double qq_CAL_relativeuncer, double FBDT_CAL, double FBDT_CAL_relativeuncer) {
	// MXs_bin: reco bin

	int NEntryFEI = ReadNFEIEigenVector("./FEI_selected.txt");
	int NEntryKID = ReadNPIDEigenVector("./KID_selected.txt");
	int NEntryPID = ReadNPIDEigenVector("./PID_selected.txt");
	int NEntryBR = ReadNBREigenVector("./BR_selected.txt");
	int NEntrypi0 = ReadNpi0EigenVector("./pi0_selected.txt");
	int NEntryMultiplicity = ReadMultiplicityInfo("./multiplicity_selected.txt");
	int NEntryFragmentation = ReadNFragmentationEigenVector("./Fragmentation_selected.txt");

	std::string bin_name = "";
	if (MXs_bin == 1) bin_name = std::string("MXs1");
	else if (MXs_bin == 2) bin_name = std::string("MXs2");
	else if (MXs_bin == 3) bin_name = std::string("MXs3");
	else {
		printf("inappropritate bin\n");
		exit(1);
	}

	// read signal template
	/* ================================ SIGNAL with true MXs1 ================================ */
	RooStats::HistFactory::Sample sig_temp_MXs1(("Signal_MXs1_nominal_" + bin_name).c_str(), "Signal_MXs1_nominal", fname);
	sig_temp_MXs1.AddHistoSys("track_eff_uncer", "Signal_MXs1_track_m", fname, "", "Signal_MXs1_track_p", fname, "");
	sig_temp_MXs1.AddHistoSys("KS0_reco_uncer", "Signal_MXs1_KS0_m", fname, "", "Signal_MXs1_KS0_p", fname, "");
	sig_temp_MXs1.AddHistoSys("BtoDtoXKL_uncer", "Signal_MXs1_BtoDtoXKL_m", fname, "", "Signal_MXs1_BtoDtoXKL_p", fname, "");
	// sig_temp_MXs1.AddHistoSys("BRBtoXKLKL_uncer", "Signal_MXs1_BRBtoXKLKL_m", fname, "", "Signal_MXs1_BRBtoXKLKL_p", fname, ""); In KumoiRD_v003 set17, it is 0
	// sig_temp_MXs1.AddHistoSys("Xnn_BR_uncer", "Signal_MXs1_Xnn_m", fname, "", "Signal_MXs1_Xnn_p", fname, ""); In KumoiRD_v003 set17, it is 0
	for (int i = 0; i < NEntryMultiplicity; i++) sig_temp_MXs1.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs1_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs1_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryFEI; i++) sig_temp_MXs1.AddHistoSys(("FEI" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs1_FEI_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs1_FEI_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) sig_temp_MXs1.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs1_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs1_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntryPID; i++) sig_temp_MXs1.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs1_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs1_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryBR; i++) sig_temp_MXs1.AddHistoSys(("BBBR" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs1_BR_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs1_BR_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntrypi0; i++) sig_temp_MXs1.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs1_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs1_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	sig_temp_MXs1.AddHistoSys("Kff1_uncer", "Signal_MXs1_Kff1_m", fname, "", "Signal_MXs1_Kff1_p", fname, "");
	sig_temp_MXs1.AddHistoSys("Kff2_uncer", "Signal_MXs1_Kff2_m", fname, "", "Signal_MXs1_Kff2_p", fname, "");
	sig_temp_MXs1.AddHistoSys("Kff3_uncer", "Signal_MXs1_Kff3_m", fname, "", "Signal_MXs1_Kff3_p", fname, "");
	sig_temp_MXs1.AddHistoSys("Kfrac_uncer", "Signal_MXs1_Kfrac_m", fname, "", "Signal_MXs1_Kfrac_p", fname, "");
	if (MXs_bin == 1) sig_temp_MXs1.AddNormFactor(("FBDT_efficiency_CAL_" + bin_name).c_str(), FBDT_CAL, FBDT_CAL, FBDT_CAL, true);
	if (MXs_bin == 1) sig_temp_MXs1.AddOverallSys(("FBDT_efficiency_uncer_" + bin_name).c_str(), 1.0 - FBDT_CAL_relativeuncer, 1.0 + FBDT_CAL_relativeuncer);
	sig_temp_MXs1.AddOverallSys("BB_counting_uncer", 0.9855, 1.0145);
	sig_temp_MXs1.ActivateStatError("Signal_MXs1_all_uncorrelated_MC_stat", fname, "");
	sig_temp_MXs1.SetNormalizeByTheory(kFALSE);
	sig_temp_MXs1.AddNormFactor("mu_MXs1", expmu, -100.0, 100.0);
	channel->AddSample(sig_temp_MXs1);

	/* ================================ SIGNAL with true MXs2 ================================ */
	RooStats::HistFactory::Sample sig_temp_MXs2(("Signal_MXs2_nominal_" + bin_name).c_str(), "Signal_MXs2_nominal", fname);
	sig_temp_MXs2.AddHistoSys("track_eff_uncer", "Signal_MXs2_track_m", fname, "", "Signal_MXs2_track_p", fname, "");
	sig_temp_MXs2.AddHistoSys("KS0_reco_uncer", "Signal_MXs2_KS0_m", fname, "", "Signal_MXs2_KS0_p", fname, "");
	sig_temp_MXs2.AddHistoSys("BtoDtoXKL_uncer", "Signal_MXs2_BtoDtoXKL_m", fname, "", "Signal_MXs2_BtoDtoXKL_p", fname, "");
	// sig_temp_MXs2.AddHistoSys("BRBtoXKLKL_uncer", "Signal_MXs2_BRBtoXKLKL_m", fname, "", "Signal_MXs2_BRBtoXKLKL_p", fname, ""); In KumoiRD_v003 set17, it is 0
	// sig_temp_MXs2.AddHistoSys("Xnn_BR_uncer", "Signal_MXs2_Xnn_m", fname, "", "Signal_MXs2_Xnn_p", fname, ""); In KumoiRD_v003 set17, it is 0
	for (int i = 0; i < NEntryMultiplicity; i++) sig_temp_MXs2.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs2_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs2_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryFEI; i++) sig_temp_MXs2.AddHistoSys(("FEI" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs2_FEI_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs2_FEI_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) sig_temp_MXs2.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs2_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs2_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntryPID; i++) sig_temp_MXs2.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs2_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs2_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryBR; i++) sig_temp_MXs2.AddHistoSys(("BBBR" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs2_BR_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs2_BR_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntrypi0; i++) sig_temp_MXs2.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs2_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs2_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	sig_temp_MXs2.AddHistoSys("Kstarff1_uncer", "Signal_MXs2_Kstarff1_m", fname, "", "Signal_MXs2_Kstarff1_p", fname, "");
	sig_temp_MXs2.AddHistoSys("Kstarff2_uncer", "Signal_MXs2_Kstarff2_m", fname, "", "Signal_MXs2_Kstarff2_p", fname, "");
	sig_temp_MXs2.AddHistoSys("Kstarff3_uncer", "Signal_MXs2_Kstarff3_m", fname, "", "Signal_MXs2_Kstarff3_p", fname, "");
	sig_temp_MXs2.AddHistoSys("Kstarff4_uncer", "Signal_MXs2_Kstarff4_m", fname, "", "Signal_MXs2_Kstarff4_p", fname, "");
	sig_temp_MXs2.AddHistoSys("Kstarff5_uncer", "Signal_MXs2_Kstarff5_m", fname, "", "Signal_MXs2_Kstarff5_p", fname, "");
	sig_temp_MXs2.AddHistoSys("Kstarff6_uncer", "Signal_MXs2_Kstarff6_m", fname, "", "Signal_MXs2_Kstarff6_p", fname, "");
	sig_temp_MXs2.AddHistoSys("Kstarff7_uncer", "Signal_MXs2_Kstarff7_m", fname, "", "Signal_MXs2_Kstarff7_p", fname, "");
	sig_temp_MXs2.AddHistoSys("Kstarff8_uncer", "Signal_MXs2_Kstarff8_m", fname, "", "Signal_MXs2_Kstarff8_p", fname, "");
	sig_temp_MXs2.AddHistoSys("Kstarff9_uncer", "Signal_MXs2_Kstarff9_m", fname, "", "Signal_MXs2_Kstarff9_p", fname, "");
	sig_temp_MXs2.AddHistoSys("Kstarff9_uncer", "Signal_MXs2_Kstarff9_m", fname, "", "Signal_MXs2_Kstarff9_p", fname, "");
	sig_temp_MXs2.AddHistoSys("Kstarfrac_uncer", "Signal_MXs2_Kstarfrac_m", fname, "", "Signal_MXs2_Kstarfrac_p", fname, "");
	// sig_temp_MXs2.AddHistoSys("mKstar_uncer", "Signal_MXs2_mKstar_m", fname, "", "Signal_MXs2_mKstar_p", fname, "");
	if (MXs_bin == 2) sig_temp_MXs2.AddNormFactor(("FBDT_efficiency_CAL_" + bin_name).c_str(), FBDT_CAL, FBDT_CAL, FBDT_CAL, true);
	if (MXs_bin == 2) sig_temp_MXs2.AddOverallSys(("FBDT_efficiency_uncer_" + bin_name).c_str(), 1.0 - FBDT_CAL_relativeuncer, 1.0 + FBDT_CAL_relativeuncer);
	sig_temp_MXs2.AddOverallSys("BB_counting_uncer", 0.9855, 1.0145);
	sig_temp_MXs2.ActivateStatError("Signal_MXs2_all_uncorrelated_MC_stat", fname, "");
	sig_temp_MXs2.SetNormalizeByTheory(kFALSE);
	sig_temp_MXs2.AddNormFactor("mu_MXs2", expmu, -100.0, 100.0);
	channel->AddSample(sig_temp_MXs2);

	/* ================================ SIGNAL with true MXs3 ================================ */
	RooStats::HistFactory::Sample sig_temp_MXs3(("Signal_MXs3_nominal_" + bin_name).c_str(), "Signal_MXs3_nominal", fname);
	sig_temp_MXs3.AddHistoSys("track_eff_uncer", "Signal_MXs3_track_m", fname, "", "Signal_MXs3_track_p", fname, "");
	sig_temp_MXs3.AddHistoSys("KS0_reco_uncer", "Signal_MXs3_KS0_m", fname, "", "Signal_MXs3_KS0_p", fname, "");
	sig_temp_MXs3.AddHistoSys("BtoDtoXKL_uncer", "Signal_MXs3_BtoDtoXKL_m", fname, "", "Signal_MXs3_BtoDtoXKL_p", fname, "");
	// sig_temp_MXs3.AddHistoSys("BRBtoXKLKL_uncer", "Signal_MXs3_BRBtoXKLKL_m", fname, "", "Signal_MXs3_BRBtoXKLKL_p", fname, ""); In KumoiRD_v003 set17, it is 0
	// sig_temp_MXs3.AddHistoSys("Xnn_BR_uncer", "Signal_MXs3_Xnn_m", fname, "", "Signal_MXs3_Xnn_p", fname, ""); In KumoiRD_v003 set17, it is 0
	for (int i = 0; i < NEntryMultiplicity; i++) sig_temp_MXs3.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs3_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs3_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryFEI; i++) sig_temp_MXs3.AddHistoSys(("FEI" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs3_FEI_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs3_FEI_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) sig_temp_MXs3.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs3_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs3_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntryPID; i++) sig_temp_MXs3.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs3_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs3_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryBR; i++) sig_temp_MXs3.AddHistoSys(("BBBR" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs3_BR_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs3_BR_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntrypi0; i++) sig_temp_MXs3.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs3_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs3_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	sig_temp_MXs3.AddHistoSys("Kff1_uncer", "Signal_MXs3_Kff1_m", fname, "", "Signal_MXs3_Kff1_p", fname, "");
	sig_temp_MXs3.AddHistoSys("Kff2_uncer", "Signal_MXs3_Kff2_m", fname, "", "Signal_MXs3_Kff2_p", fname, "");
	sig_temp_MXs3.AddHistoSys("Kff3_uncer", "Signal_MXs3_Kff3_m", fname, "", "Signal_MXs3_Kff3_p", fname, "");
	sig_temp_MXs3.AddHistoSys("Kstarff1_uncer", "Signal_MXs3_Kstarff1_m", fname, "", "Signal_MXs3_Kstarff1_p", fname, "");
	sig_temp_MXs3.AddHistoSys("Kstarff2_uncer", "Signal_MXs3_Kstarff2_m", fname, "", "Signal_MXs3_Kstarff2_p", fname, "");
	sig_temp_MXs3.AddHistoSys("Kstarff3_uncer", "Signal_MXs3_Kstarff3_m", fname, "", "Signal_MXs3_Kstarff3_p", fname, "");
	sig_temp_MXs3.AddHistoSys("Kstarff4_uncer", "Signal_MXs3_Kstarff4_m", fname, "", "Signal_MXs3_Kstarff4_p", fname, "");
	sig_temp_MXs3.AddHistoSys("Kstarff5_uncer", "Signal_MXs3_Kstarff5_m", fname, "", "Signal_MXs3_Kstarff5_p", fname, "");
	sig_temp_MXs3.AddHistoSys("Kstarff6_uncer", "Signal_MXs3_Kstarff6_m", fname, "", "Signal_MXs3_Kstarff6_p", fname, "");
	sig_temp_MXs3.AddHistoSys("Kstarff7_uncer", "Signal_MXs3_Kstarff7_m", fname, "", "Signal_MXs3_Kstarff7_p", fname, "");
	sig_temp_MXs3.AddHistoSys("Kstarff8_uncer", "Signal_MXs3_Kstarff8_m", fname, "", "Signal_MXs3_Kstarff8_p", fname, "");
	sig_temp_MXs3.AddHistoSys("Kstarff9_uncer", "Signal_MXs3_Kstarff9_m", fname, "", "Signal_MXs3_Kstarff9_p", fname, "");
	sig_temp_MXs3.AddHistoSys("Kstarff9_uncer", "Signal_MXs3_Kstarff9_m", fname, "", "Signal_MXs3_Kstarff9_p", fname, "");
	sig_temp_MXs3.AddHistoSys("Kfrac_uncer", "Signal_MXs3_Kfrac_m", fname, "", "Signal_MXs3_Kfrac_p", fname, "");
	sig_temp_MXs3.AddHistoSys("Kstarfrac_uncer", "Signal_MXs3_Kstarfrac_m", fname, "", "Signal_MXs3_Kstarfrac_p", fname, "");
	for (int i = 0; i < NEntryFragmentation; i++) sig_temp_MXs3.AddHistoSys(("Xs_fragmentation" + std::to_string(i) + "_uncer").c_str(), ("Signal_MXs3_Fragmentation_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_MXs3_Fragmentation_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	sig_temp_MXs3.AddHistoSys("pf_uncer", "Signal_MXs3_pf_m", fname, "", "Signal_MXs3_pf_p", fname, "");
	sig_temp_MXs3.AddHistoSys("mb_uncer", "Signal_MXs3_mb_m", fname, "", "Signal_MXs3_mb_p", fname, "");
	sig_temp_MXs3.AddHistoSys("transition_uncer", "Signal_MXs3_transition_m", fname, "", "Signal_MXs3_transition_p", fname, "");
	// sig_temp_MXs3.AddHistoSys("mKstar_uncer", "Signal_MXs3_mKstar_m", fname, "", "Signal_MXs3_mKstar_p", fname, "");
	if (MXs_bin == 3) sig_temp_MXs3.AddNormFactor(("FBDT_efficiency_CAL_" + bin_name).c_str(), FBDT_CAL, FBDT_CAL, FBDT_CAL, true);
	if (MXs_bin == 3) sig_temp_MXs3.AddOverallSys(("FBDT_efficiency_uncer_" + bin_name).c_str(), 1.0 - FBDT_CAL_relativeuncer, 1.0 + FBDT_CAL_relativeuncer);
	sig_temp_MXs3.AddOverallSys("BB_counting_uncer", 0.9855, 1.0145);
	sig_temp_MXs3.ActivateStatError("Signal_MXs3_all_uncorrelated_MC_stat", fname, "");
	sig_temp_MXs3.SetNormalizeByTheory(kFALSE);
	sig_temp_MXs3.AddNormFactor("mu_MXs3", expmu, -100.0, 100.0);
	channel->AddSample(sig_temp_MXs3);

	// read background template
	/* ================================ CHG ================================ */
	RooStats::HistFactory::Sample CHG_temp(("CHG_nominal_" + bin_name).c_str(), "CHG_nominal", fname);
	CHG_temp.AddHistoSys("track_eff_uncer", "CHG_track_m", fname, "", "CHG_track_p", fname, "");
	CHG_temp.AddHistoSys("KS0_reco_uncer", "CHG_KS0_m", fname, "", "CHG_KS0_p", fname, "");
	CHG_temp.AddHistoSys("BtoDtoXKL_uncer", "CHG_BtoDtoXKL_m", fname, "", "CHG_BtoDtoXKL_p", fname, "");
	CHG_temp.AddHistoSys("BRBtoXKLKL_uncer", "CHG_BRBtoXKLKL_m", fname, "", "CHG_BRBtoXKLKL_p", fname, "");
	CHG_temp.AddHistoSys("Xnn_BR_uncer", "CHG_Xnn_m", fname, "", "CHG_Xnn_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) CHG_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("CHG_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryFEI; i++) CHG_temp.AddHistoSys(("FEI" + std::to_string(i) + "_uncer").c_str(), ("CHG_FEI_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_FEI_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) CHG_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("CHG_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntryPID; i++) CHG_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("CHG_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryBR; i++) CHG_temp.AddHistoSys(("BBBR" + std::to_string(i) + "_uncer").c_str(), ("CHG_BR_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_BR_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntrypi0; i++) CHG_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("CHG_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	//CHG_temp.AddHistoSys("BDTc_shape_BB", "CHG_BDTc_m", fname, "", "CHG_BDTc_p", fname, "");
	CHG_temp.AddOverallSys(("mu_CHG_" + bin_name).c_str(), 0.7, 1.3);
	CHG_temp.AddOverallSys("BB_counting_uncer", 0.9855, 1.0145);
	CHG_temp.ActivateStatError("CHG_all_uncorrelated_MC_stat", fname, "");
	CHG_temp.SetNormalizeByTheory(kFALSE);
	channel->AddSample(CHG_temp);
	/* ================================ CHG ================================ */

	/* ================================ MIX ================================ */
	RooStats::HistFactory::Sample MIX_temp(("MIX_nominal_" + bin_name).c_str(), "MIX_nominal", fname);
	MIX_temp.AddHistoSys("track_eff_uncer", "MIX_track_m", fname, "", "MIX_track_p", fname, "");
	MIX_temp.AddHistoSys("KS0_reco_uncer", "MIX_KS0_m", fname, "", "MIX_KS0_p", fname, "");
	MIX_temp.AddHistoSys("BtoDtoXKL_uncer", "MIX_BtoDtoXKL_m", fname, "", "MIX_BtoDtoXKL_p", fname, "");
	MIX_temp.AddHistoSys("BRBtoXKLKL_uncer", "MIX_BRBtoXKLKL_m", fname, "", "MIX_BRBtoXKLKL_p", fname, "");
	MIX_temp.AddHistoSys("Xnn_BR_uncer", "MIX_Xnn_m", fname, "", "MIX_Xnn_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) MIX_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("MIX_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryFEI; i++) MIX_temp.AddHistoSys(("FEI" + std::to_string(i) + "_uncer").c_str(), ("MIX_FEI_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_FEI_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) MIX_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("MIX_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntryPID; i++) MIX_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("MIX_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryBR; i++) MIX_temp.AddHistoSys(("BBBR" + std::to_string(i) + "_uncer").c_str(), ("MIX_BR_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_BR_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntrypi0; i++) MIX_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("MIX_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	//MIX_temp.AddHistoSys("BDTc_shape_BB", "MIX_BDTc_m", fname, "", "MIX_BDTc_p", fname, "");
	MIX_temp.AddOverallSys(("mu_MIX_" + bin_name).c_str(), 0.7, 1.3);
	MIX_temp.AddOverallSys("BB_counting_uncer", 0.9855, 1.0145);
	MIX_temp.ActivateStatError("MIX_all_uncorrelated_MC_stat", fname, "");
	MIX_temp.SetNormalizeByTheory(kFALSE);
	channel->AddSample(MIX_temp);
	/* ================================ MIX ================================ */

	/* ================================ UUBAR ================================ */
	RooStats::HistFactory::Sample UUBAR_temp(("UUBAR_nominal_" + bin_name).c_str(), "UUBAR_nominal", fname);
	UUBAR_temp.AddHistoSys("track_eff_uncer", "UUBAR_track_m", fname, "", "UUBAR_track_p", fname, "");
	UUBAR_temp.AddHistoSys("KS0_reco_uncer", "UUBAR_KS0_m", fname, "", "UUBAR_KS0_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) UUBAR_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("UUBAR_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("UUBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) UUBAR_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("UUBAR_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("UUBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntryPID; i++) UUBAR_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("UUBAR_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("UUBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntrypi0; i++) UUBAR_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("UUBAR_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("UUBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	UUBAR_temp.AddHistoSys("BDTc_shape_qq", "UUBAR_BDTc_m", fname, "", "UUBAR_BDTc_p", fname, "");
	UUBAR_temp.AddNormFactor(("qq_CAL_UUBAR_" + bin_name).c_str(), qq_CAL, qq_CAL, qq_CAL, true);
	//UUBAR_temp.AddOverallSys(("qq_CAL_UUBAR_uncer_" + bin_name).c_str(), 1.0 - qq_CAL_relativeuncer, 1.0 + qq_CAL_relativeuncer);
	UUBAR_temp.AddOverallSys(("mu_UUBAR_" + bin_name).c_str(), 1.0 - std::sqrt(0.3 * 0.3 + qq_CAL_relativeuncer * qq_CAL_relativeuncer), 1.0 + std::sqrt(0.3 * 0.3 + qq_CAL_relativeuncer * qq_CAL_relativeuncer));
	UUBAR_temp.ActivateStatError("UUBAR_all_uncorrelated_MC_stat", fname, "");
	UUBAR_temp.SetNormalizeByTheory(kFALSE);
	channel->AddSample(UUBAR_temp);
	/* ================================ UUBAR ================================ */

	/* ================================ DDBAR ================================ */
	RooStats::HistFactory::Sample DDBAR_temp(("DDBAR_nominal_" + bin_name).c_str(), "DDBAR_nominal", fname);
	DDBAR_temp.AddHistoSys("track_eff_uncer", "DDBAR_track_m", fname, "", "DDBAR_track_p", fname, "");
	DDBAR_temp.AddHistoSys("KS0_reco_uncer", "DDBAR_KS0_m", fname, "", "DDBAR_KS0_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) DDBAR_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("DDBAR_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("DDBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) DDBAR_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("DDBAR_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("DDBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntryPID; i++) DDBAR_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("DDBAR_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("DDBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntrypi0; i++) DDBAR_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("DDBAR_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("DDBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	DDBAR_temp.AddHistoSys("BDTc_shape_qq", "DDBAR_BDTc_m", fname, "", "DDBAR_BDTc_p", fname, "");
	DDBAR_temp.AddNormFactor(("qq_CAL_DDBAR_" + bin_name).c_str(), qq_CAL, qq_CAL, qq_CAL, true);
	//DDBAR_temp.AddOverallSys(("qq_CAL_DDBAR_uncer_" + bin_name).c_str(), 1.0 - qq_CAL_relativeuncer, 1.0 + qq_CAL_relativeuncer);
	DDBAR_temp.AddOverallSys(("mu_DDBAR_" + bin_name).c_str(), 1.0 - std::sqrt(0.3 * 0.3 + qq_CAL_relativeuncer * qq_CAL_relativeuncer), 1.0 + std::sqrt(0.3 * 0.3 + qq_CAL_relativeuncer * qq_CAL_relativeuncer));
	DDBAR_temp.ActivateStatError("DDBAR_all_uncorrelated_MC_stat", fname, "");
	DDBAR_temp.SetNormalizeByTheory(kFALSE);
	channel->AddSample(DDBAR_temp);
	/* ================================ DDBAR ================================ */

	/* ================================ SSBAR ================================ */
	RooStats::HistFactory::Sample SSBAR_temp(("SSBAR_nominal_" + bin_name).c_str(), "SSBAR_nominal", fname);
	SSBAR_temp.AddHistoSys("track_eff_uncer", "SSBAR_track_m", fname, "", "SSBAR_track_p", fname, "");
	SSBAR_temp.AddHistoSys("KS0_reco_uncer", "SSBAR_KS0_m", fname, "", "SSBAR_KS0_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) SSBAR_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("SSBAR_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("SSBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) SSBAR_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("SSBAR_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("SSBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntryPID; i++) SSBAR_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("SSBAR_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("SSBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntrypi0; i++) SSBAR_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("SSBAR_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("SSBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	SSBAR_temp.AddHistoSys("BDTc_shape_qq", "SSBAR_BDTc_m", fname, "", "SSBAR_BDTc_p", fname, "");
	SSBAR_temp.AddNormFactor(("qq_CAL_SSBAR_" + bin_name).c_str(), qq_CAL, qq_CAL, qq_CAL, true);
	//SSBAR_temp.AddOverallSys(("qq_CAL_SSBAR_uncer_" + bin_name).c_str(), 1.0 - qq_CAL_relativeuncer, 1.0 + qq_CAL_relativeuncer);
	SSBAR_temp.AddOverallSys(("mu_SSBAR_" + bin_name).c_str(), 1.0 - std::sqrt(0.3 * 0.3 + qq_CAL_relativeuncer * qq_CAL_relativeuncer), 1.0 + std::sqrt(0.3 * 0.3 + qq_CAL_relativeuncer * qq_CAL_relativeuncer));
	SSBAR_temp.ActivateStatError("SSBAR_all_uncorrelated_MC_stat", fname, "");
	SSBAR_temp.SetNormalizeByTheory(kFALSE);
	channel->AddSample(SSBAR_temp);
	/* ================================ SSBAR ================================ */

	/* ================================ CHARM ================================ */
	RooStats::HistFactory::Sample CHARM_temp(("CHARM_nominal_" + bin_name).c_str(), "CHARM_nominal", fname);
	CHARM_temp.AddHistoSys("track_eff_uncer", "CHARM_track_m", fname, "", "CHARM_track_p", fname, "");
	CHARM_temp.AddHistoSys("KS0_reco_uncer", "CHARM_KS0_m", fname, "", "CHARM_KS0_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) CHARM_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("CHARM_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHARM_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) CHARM_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("CHARM_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHARM_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntryPID; i++) CHARM_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("CHARM_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHARM_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	if (MXs_bin != 1) for (int i = 0; i < NEntrypi0; i++) CHARM_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("CHARM_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHARM_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	CHARM_temp.AddHistoSys("BDTc_shape_qq", "CHARM_BDTc_m", fname, "", "CHARM_BDTc_p", fname, "");
	CHARM_temp.AddNormFactor(("qq_CAL_CHARM_" + bin_name).c_str(), qq_CAL, qq_CAL, qq_CAL, true);
	//CHARM_temp.AddOverallSys(("qq_CAL_CHARM_uncer_" + bin_name).c_str(), 1.0 - qq_CAL_relativeuncer, 1.0 + qq_CAL_relativeuncer);
	CHARM_temp.AddOverallSys(("mu_CHARM_" + bin_name).c_str(), 1.0 - std::sqrt(0.3 * 0.3 + qq_CAL_relativeuncer * qq_CAL_relativeuncer), 1.0 + std::sqrt(0.3 * 0.3 + qq_CAL_relativeuncer * qq_CAL_relativeuncer));
	CHARM_temp.ActivateStatError("CHARM_all_uncorrelated_MC_stat", fname, "");
	CHARM_temp.SetNormalizeByTheory(kFALSE);
	channel->AddSample(CHARM_temp);
	/* ================================ CHARM ================================ */
}

int WorkSpace() {

	const double expmu = 1.0;
	const char* fname = "PDFandDATA_nominal.root";
	const char* fname_MXs1 = "PDFandDATA_MXs1.root";
	const char* fname_MXs2 = "PDFandDATA_MXs2.root";
	const char* fname_MXs3 = "PDFandDATA_MXs3.root";

	RooStats::HistFactory::Measurement meas("my_measurement", "my measurement");
	meas.SetOutputFilePrefix("results/my_measurement");
	meas.SetExportOnly(kTRUE);

	// setting measurement
	meas.SetPOI("mu_MXs1");
	meas.SetLumi(1.0);
	meas.AddConstantParam("Lumi");

	// define channels
	HistFactory::Channel channel_MXs1("channel_MXs1");
	channel_MXs1.SetStatErrorConfig(1e-5, "Gaussian");

	HistFactory::Channel channel_MXs2("channel_MXs2");
	channel_MXs2.SetStatErrorConfig(1e-5, "Gaussian");

	HistFactory::Channel channel_MXs3("channel_MXs3");
	channel_MXs3.SetStatErrorConfig(1e-5, "Gaussian");

	// point data
	channel_MXs1.SetData("total_DATA", fname_MXs1);
	channel_MXs2.SetData("total_DATA", fname_MXs2);
	channel_MXs3.SetData("total_DATA", fname_MXs3);

	// get MC
	AddSample(&channel_MXs1, fname_MXs1, 1, expmu, 1.1652, 0.1105, 1.0062, 0.0381);
	AddSample(&channel_MXs2, fname_MXs2, 2, expmu, 0.8544, 0.2287, 1.1006, 0.0742);
	AddSample(&channel_MXs3, fname_MXs3, 3, expmu, 0.7892, 0.1311, 0.9391, 0.1513);

	// add channel to measurement
	meas.AddChannel(channel_MXs1);
	meas.AddChannel(channel_MXs2);
	meas.AddChannel(channel_MXs3);
	meas.CollectHistograms();

	RooWorkspace* w;
	w = RooStats::HistFactory::MakeModelAndMeasurementFast(meas);

	w->Print();
	w->writeToFile("PDFandDATA_workspace.root");

	return 0;
}