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
	double weight_sys[RarityBins * 3] = { 0.0 }; // eigen vector

	FILE* fp;
	fp = fopen(dirname, "r");
	while (true) {
		if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
		for (int i = 0; i < RarityBins * 3; i++) {
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
	double weight_sys[RarityBins * 7] = { 0.0 }; // eigen vector

	FILE* fp;
	fp = fopen(dirname, "r");
	while (true) {
		if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
		for (int i = 0; i < RarityBins * 7; i++) {
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
	double weight_sys[RarityBins * 3] = { 0.0 }; // eigen vector

	FILE* fp;
	fp = fopen(dirname, "r");
	while (true) {
		if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
		for (int i = 0; i < RarityBins * 3; i++) {
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
	double weight_sys[RarityBins * 7] = { 0.0 }; // eigen vector

	FILE* fp;
	fp = fopen(dirname, "r");
	while (true) {
		if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
		for (int i = 0; i < RarityBins * 7; i++) {
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
	double weight_sys[RarityBins * 7] = { 0.0 }; // eigen vector

	FILE* fp;
	fp = fopen(dirname, "r");
	while (true) {
		if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
		for (int i = 0; i < RarityBins * 7; i++) {
			if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
		}
		Nentry++;
	}
	fclose(fp);

	return Nentry;
}

int ReadNfractionEigenVector(const char* dirname) {
	int Nentry = 0; // number of eigen values/vectors
	double eigen_value = 0; // eigen value
	double weight_sys[RarityBins * 3] = { 0.0 }; // eigen vector

	FILE* fp;
	fp = fopen(dirname, "r");
	while (true) {
		if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
		for (int i = 0; i < RarityBins * 3; i++) {
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
	double weight_sys[RarityBins * 3] = { 0.0 }; // eigen vector

	FILE* fp;
	fp = fopen(dirname, "r");
	while (true) {
		if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
		for (int i = 0; i < RarityBins * 3; i++) {
			if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
		}
		Nentry++;
	}
	fclose(fp);

	return Nentry;
}

int WorkSpace() {

	int NEntryFEI = ReadNFEIEigenVector("./FEI_selected.txt");
	int NEntryKID = ReadNPIDEigenVector("./KID_selected.txt");
	int NEntryPID = ReadNPIDEigenVector("./PID_selected.txt");
	int NEntryBR = ReadNBREigenVector("./BR_selected.txt");
	int NEntrypi0 = ReadNpi0EigenVector("./pi0_selected.txt");
	int NEntryMultiplicity = ReadMultiplicityInfo("./multiplicity_selected.txt");
	int NEntryfraction = ReadNfractionEigenVector("./fraction_selected.txt");
	int NEntryFragmentation = ReadNFragmentationEigenVector("./Fragmentation_selected.txt");

	const double expmu = 1.0;
	const char* fname = "PDFandDATA_nominal.root";
	const char* fname_MXs1 = "PDFandDATA_MXs1.root";
	const char* fname_MXs2 = "PDFandDATA_MXs2.root";
	const char* fname_MXs3 = "PDFandDATA_MXs3.root";

	RooStats::HistFactory::Measurement meas("my_measurement", "my measurement");
	meas.SetOutputFilePrefix("results/my_measurement");
	meas.SetExportOnly(kTRUE);

	// setting measurement
	meas.SetPOI("mu_MXs3");
	meas.SetLumi(1.0);
	meas.AddConstantParam("Lumi");

	// define channel
	HistFactory::Channel channel("channel");
	channel.SetStatErrorConfig(1e-5, "Gaussian");

	// point data
	channel.SetData("total_DATA", fname);

	// read signal template
	/* ================================ SIGNAL MX1 ================================ */
	RooStats::HistFactory::Sample sig_temp_MXs1("Signal_nominal_MXs1", "Signal_nominal", fname_MXs1);
	sig_temp_MXs1.AddHistoSys("track_eff_uncer", "Signal_track_m", fname_MXs1, "", "Signal_track_p", fname_MXs1, "");
	sig_temp_MXs1.AddHistoSys("KS0_reco_uncer", "Signal_KS0_m", fname_MXs1, "", "Signal_KS0_p", fname_MXs1, "");
	sig_temp_MXs1.AddHistoSys("BtoDtoXKL_uncer", "Signal_BtoDtoXKL_m", fname_MXs1, "", "Signal_BtoDtoXKL_p", fname_MXs1, "");
	sig_temp_MXs1.AddHistoSys("BRBtoXKLKL_uncer", "Signal_BRBtoXKLKL_m", fname_MXs1, "", "Signal_BRBtoXKLKL_p", fname_MXs1, "");
	sig_temp_MXs1.AddHistoSys("Xnn_BR_uncer", "Signal_Xnn_m", fname_MXs1, "", "Signal_Xnn_p", fname_MXs1, "");
	for (int i = 0; i < NEntryMultiplicity; i++) sig_temp_MXs1.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("Signal_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs1, "", ("Signal_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs1, "");
	for (int i = 0; i < NEntryFEI; i++) sig_temp_MXs1.AddHistoSys(("FEI" + std::to_string(i) + "_uncer").c_str(), ("Signal_FEI_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs1, "", ("Signal_FEI_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs1, "");
	for (int i = 0; i < NEntryKID; i++) sig_temp_MXs1.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("Signal_KID_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs1, "", ("Signal_KID_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs1, "");
	for (int i = 0; i < NEntryPID; i++) sig_temp_MXs1.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("Signal_PID_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs1, "", ("Signal_PID_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs1, "");
	for (int i = 0; i < NEntryBR; i++) sig_temp_MXs1.AddHistoSys(("BBBR" + std::to_string(i) + "_uncer").c_str(), ("Signal_BR_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs1, "", ("Signal_BR_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs1, "");
	for (int i = 0; i < NEntrypi0; i++) sig_temp_MXs1.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("Signal_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs1, "", ("Signal_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs1, "");
	sig_temp_MXs1.AddHistoSys("Kff1_uncer", "Signal_Kff1_m", fname_MXs1, "", "Signal_Kff1_p", fname_MXs1, "");
	sig_temp_MXs1.AddHistoSys("Kff2_uncer", "Signal_Kff2_m", fname_MXs1, "", "Signal_Kff2_p", fname_MXs1, "");
	sig_temp_MXs1.AddHistoSys("Kff3_uncer", "Signal_Kff3_m", fname_MXs1, "", "Signal_Kff3_p", fname_MXs1, "");
	sig_temp_MXs1.AddHistoSys("Kstarff1_uncer", "Signal_Kstarff1_m", fname_MXs1, "", "Signal_Kstarff1_p", fname_MXs1, "");
	sig_temp_MXs1.AddHistoSys("Kstarff2_uncer", "Signal_Kstarff2_m", fname_MXs1, "", "Signal_Kstarff2_p", fname_MXs1, "");
	sig_temp_MXs1.AddHistoSys("Kstarff3_uncer", "Signal_Kstarff3_m", fname_MXs1, "", "Signal_Kstarff3_p", fname_MXs1, "");
	sig_temp_MXs1.AddHistoSys("Kstarff4_uncer", "Signal_Kstarff4_m", fname_MXs1, "", "Signal_Kstarff4_p", fname_MXs1, "");
	sig_temp_MXs1.AddHistoSys("Kstarff5_uncer", "Signal_Kstarff5_m", fname_MXs1, "", "Signal_Kstarff5_p", fname_MXs1, "");
	sig_temp_MXs1.AddHistoSys("Kstarff6_uncer", "Signal_Kstarff6_m", fname_MXs1, "", "Signal_Kstarff6_p", fname_MXs1, "");
	sig_temp_MXs1.AddHistoSys("Kstarff7_uncer", "Signal_Kstarff7_m", fname_MXs1, "", "Signal_Kstarff7_p", fname_MXs1, "");
	sig_temp_MXs1.AddHistoSys("Kstarff8_uncer", "Signal_Kstarff8_m", fname_MXs1, "", "Signal_Kstarff8_p", fname_MXs1, "");
	sig_temp_MXs1.AddHistoSys("Kstarff9_uncer", "Signal_Kstarff9_m", fname_MXs1, "", "Signal_Kstarff9_p", fname_MXs1, "");
	sig_temp_MXs1.AddHistoSys("Kstarff9_uncer", "Signal_Kstarff9_m", fname_MXs1, "", "Signal_Kstarff9_p", fname_MXs1, "");
	sig_temp_MXs1.AddHistoSys("Kff_OLD_uncer", "Signal_Kff_OLD_m", fname_MXs1, "", "Signal_Kff_OLD_p", fname_MXs1, "");
	for (int i = 0; i < NEntryfraction; i++) sig_temp_MXs1.AddHistoSys(("fraction" + std::to_string(i) + "_uncer").c_str(), ("Signal_fraction_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs1, "", ("Signal_fraction_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs1, "");
	for (int i = 0; i < NEntryFragmentation; i++) sig_temp_MXs1.AddHistoSys(("Xs_fragmentation" + std::to_string(i) + "_uncer").c_str(), ("Signal_Fragmentation_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs1, "", ("Signal_Fragmentation_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs1, "");
	sig_temp_MXs1.AddHistoSys("pf_uncer", "Signal_pf_m", fname_MXs1, "", "Signal_pf_p", fname_MXs1, "");
	sig_temp_MXs1.AddHistoSys("mb_uncer", "Signal_mb_m", fname_MXs1, "", "Signal_mb_p", fname_MXs1, "");
	sig_temp_MXs1.AddHistoSys("transition_uncer", "Signal_transition_m", fname_MXs1, "", "Signal_transition_p", fname_MXs1, "");
	sig_temp_MXs1.AddHistoSys("mKstar_uncer", "Signal_mKstar_m", fname_MXs1, "", "Signal_mKstar_p", fname_MXs1, "");
	sig_temp_MXs1.AddNormFactor("FBDT_efficiency_CAL", 0.8605, 0.8605, 0.8605, true);
	sig_temp_MXs1.AddOverallSys("FBDT_efficiency_uncer", 0.9262, 1.0737);
	sig_temp_MXs1.AddOverallSys("BB_counting_uncer", 0.9855, 1.0145);
	sig_temp_MXs1.ActivateStatError("Signal_all_uncorrelated_MC_stat", fname_MXs1, "");
	sig_temp_MXs1.SetNormalizeByTheory(kFALSE);
	sig_temp_MXs1.AddNormFactor("mu_MXs1", expmu, -100.0, 100.0);
	channel.AddSample(sig_temp_MXs1);

	/* ================================ SIGNAL MX2 ================================ */
	RooStats::HistFactory::Sample sig_temp_MXs2("Signal_nominal_MXs2", "Signal_nominal", fname_MXs2);
	sig_temp_MXs2.AddHistoSys("track_eff_uncer", "Signal_track_m", fname_MXs2, "", "Signal_track_p", fname_MXs2, "");
	sig_temp_MXs2.AddHistoSys("KS0_reco_uncer", "Signal_KS0_m", fname_MXs2, "", "Signal_KS0_p", fname_MXs2, "");
	sig_temp_MXs2.AddHistoSys("BtoDtoXKL_uncer", "Signal_BtoDtoXKL_m", fname_MXs2, "", "Signal_BtoDtoXKL_p", fname_MXs2, "");
	sig_temp_MXs2.AddHistoSys("BRBtoXKLKL_uncer", "Signal_BRBtoXKLKL_m", fname_MXs2, "", "Signal_BRBtoXKLKL_p", fname_MXs2, "");
	sig_temp_MXs2.AddHistoSys("Xnn_BR_uncer", "Signal_Xnn_m", fname_MXs2, "", "Signal_Xnn_p", fname_MXs2, "");
	for (int i = 0; i < NEntryMultiplicity; i++) sig_temp_MXs2.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("Signal_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs2, "", ("Signal_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs2, "");
	for (int i = 0; i < NEntryFEI; i++) sig_temp_MXs2.AddHistoSys(("FEI" + std::to_string(i) + "_uncer").c_str(), ("Signal_FEI_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs2, "", ("Signal_FEI_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs2, "");
	for (int i = 0; i < NEntryKID; i++) sig_temp_MXs2.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("Signal_KID_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs2, "", ("Signal_KID_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs2, "");
	for (int i = 0; i < NEntryPID; i++) sig_temp_MXs2.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("Signal_PID_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs2, "", ("Signal_PID_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs2, "");
	for (int i = 0; i < NEntryBR; i++) sig_temp_MXs2.AddHistoSys(("BBBR" + std::to_string(i) + "_uncer").c_str(), ("Signal_BR_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs2, "", ("Signal_BR_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs2, "");
	for (int i = 0; i < NEntrypi0; i++) sig_temp_MXs2.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("Signal_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs2, "", ("Signal_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs2, "");
	sig_temp_MXs2.AddHistoSys("Kff1_uncer", "Signal_Kff1_m", fname_MXs2, "", "Signal_Kff1_p", fname_MXs2, "");
	sig_temp_MXs2.AddHistoSys("Kff2_uncer", "Signal_Kff2_m", fname_MXs2, "", "Signal_Kff2_p", fname_MXs2, "");
	sig_temp_MXs2.AddHistoSys("Kff3_uncer", "Signal_Kff3_m", fname_MXs2, "", "Signal_Kff3_p", fname_MXs2, "");
	sig_temp_MXs2.AddHistoSys("Kstarff1_uncer", "Signal_Kstarff1_m", fname_MXs2, "", "Signal_Kstarff1_p", fname_MXs2, "");
	sig_temp_MXs2.AddHistoSys("Kstarff2_uncer", "Signal_Kstarff2_m", fname_MXs2, "", "Signal_Kstarff2_p", fname_MXs2, "");
	sig_temp_MXs2.AddHistoSys("Kstarff3_uncer", "Signal_Kstarff3_m", fname_MXs2, "", "Signal_Kstarff3_p", fname_MXs2, "");
	sig_temp_MXs2.AddHistoSys("Kstarff4_uncer", "Signal_Kstarff4_m", fname_MXs2, "", "Signal_Kstarff4_p", fname_MXs2, "");
	sig_temp_MXs2.AddHistoSys("Kstarff5_uncer", "Signal_Kstarff5_m", fname_MXs2, "", "Signal_Kstarff5_p", fname_MXs2, "");
	sig_temp_MXs2.AddHistoSys("Kstarff6_uncer", "Signal_Kstarff6_m", fname_MXs2, "", "Signal_Kstarff6_p", fname_MXs2, "");
	sig_temp_MXs2.AddHistoSys("Kstarff7_uncer", "Signal_Kstarff7_m", fname_MXs2, "", "Signal_Kstarff7_p", fname_MXs2, "");
	sig_temp_MXs2.AddHistoSys("Kstarff8_uncer", "Signal_Kstarff8_m", fname_MXs2, "", "Signal_Kstarff8_p", fname_MXs2, "");
	sig_temp_MXs2.AddHistoSys("Kstarff9_uncer", "Signal_Kstarff9_m", fname_MXs2, "", "Signal_Kstarff9_p", fname_MXs2, "");
	sig_temp_MXs2.AddHistoSys("Kstarff9_uncer", "Signal_Kstarff9_m", fname_MXs2, "", "Signal_Kstarff9_p", fname_MXs2, "");
	sig_temp_MXs2.AddHistoSys("Kff_OLD_uncer", "Signal_Kff_OLD_m", fname_MXs2, "", "Signal_Kff_OLD_p", fname_MXs2, "");
	for (int i = 0; i < NEntryfraction; i++) sig_temp_MXs2.AddHistoSys(("fraction" + std::to_string(i) + "_uncer").c_str(), ("Signal_fraction_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs2, "", ("Signal_fraction_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs2, "");
	for (int i = 0; i < NEntryFragmentation; i++) sig_temp_MXs2.AddHistoSys(("Xs_fragmentation" + std::to_string(i) + "_uncer").c_str(), ("Signal_Fragmentation_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs2, "", ("Signal_Fragmentation_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs2, "");
	sig_temp_MXs2.AddHistoSys("pf_uncer", "Signal_pf_m", fname_MXs2, "", "Signal_pf_p", fname_MXs2, "");
	sig_temp_MXs2.AddHistoSys("mb_uncer", "Signal_mb_m", fname_MXs2, "", "Signal_mb_p", fname_MXs2, "");
	sig_temp_MXs2.AddHistoSys("transition_uncer", "Signal_transition_m", fname_MXs2, "", "Signal_transition_p", fname_MXs2, "");
	sig_temp_MXs2.AddHistoSys("mKstar_uncer", "Signal_mKstar_m", fname_MXs2, "", "Signal_mKstar_p", fname_MXs2, "");
	sig_temp_MXs2.AddNormFactor("FBDT_efficiency_CAL", 0.8605, 0.8605, 0.8605, true);
	sig_temp_MXs2.AddOverallSys("FBDT_efficiency_uncer", 0.9262, 1.0737);
	sig_temp_MXs2.AddOverallSys("BB_counting_uncer", 0.9855, 1.0145);
	sig_temp_MXs2.ActivateStatError("Signal_all_uncorrelated_MC_stat", fname_MXs2, "");
	sig_temp_MXs2.SetNormalizeByTheory(kFALSE);
	sig_temp_MXs2.AddNormFactor("mu_MXs2", expmu, -100.0, 100.0);
	channel.AddSample(sig_temp_MXs2);

	/* ================================ SIGNAL MX3 ================================ */
	RooStats::HistFactory::Sample sig_temp_MXs3("Signal_nominal_MXs3", "Signal_nominal", fname_MXs3);
	sig_temp_MXs3.AddHistoSys("track_eff_uncer", "Signal_track_m", fname_MXs3, "", "Signal_track_p", fname_MXs3, "");
	sig_temp_MXs3.AddHistoSys("KS0_reco_uncer", "Signal_KS0_m", fname_MXs3, "", "Signal_KS0_p", fname_MXs3, "");
	sig_temp_MXs3.AddHistoSys("BtoDtoXKL_uncer", "Signal_BtoDtoXKL_m", fname_MXs3, "", "Signal_BtoDtoXKL_p", fname_MXs3, "");
	sig_temp_MXs3.AddHistoSys("BRBtoXKLKL_uncer", "Signal_BRBtoXKLKL_m", fname_MXs3, "", "Signal_BRBtoXKLKL_p", fname_MXs3, "");
	sig_temp_MXs3.AddHistoSys("Xnn_BR_uncer", "Signal_Xnn_m", fname_MXs3, "", "Signal_Xnn_p", fname_MXs3, "");
	for (int i = 0; i < NEntryMultiplicity; i++) sig_temp_MXs3.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("Signal_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs3, "", ("Signal_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs3, "");
	for (int i = 0; i < NEntryFEI; i++) sig_temp_MXs3.AddHistoSys(("FEI" + std::to_string(i) + "_uncer").c_str(), ("Signal_FEI_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs3, "", ("Signal_FEI_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs3, "");
	for (int i = 0; i < NEntryKID; i++) sig_temp_MXs3.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("Signal_KID_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs3, "", ("Signal_KID_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs3, "");
	for (int i = 0; i < NEntryPID; i++) sig_temp_MXs3.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("Signal_PID_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs3, "", ("Signal_PID_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs3, "");
	for (int i = 0; i < NEntryBR; i++) sig_temp_MXs3.AddHistoSys(("BBBR" + std::to_string(i) + "_uncer").c_str(), ("Signal_BR_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs3, "", ("Signal_BR_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs3, "");
	for (int i = 0; i < NEntrypi0; i++) sig_temp_MXs3.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("Signal_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs3, "", ("Signal_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs3, "");
	sig_temp_MXs3.AddHistoSys("Kff1_uncer", "Signal_Kff1_m", fname_MXs3, "", "Signal_Kff1_p", fname_MXs3, "");
	sig_temp_MXs3.AddHistoSys("Kff2_uncer", "Signal_Kff2_m", fname_MXs3, "", "Signal_Kff2_p", fname_MXs3, "");
	sig_temp_MXs3.AddHistoSys("Kff3_uncer", "Signal_Kff3_m", fname_MXs3, "", "Signal_Kff3_p", fname_MXs3, "");
	sig_temp_MXs3.AddHistoSys("Kstarff1_uncer", "Signal_Kstarff1_m", fname_MXs3, "", "Signal_Kstarff1_p", fname_MXs3, "");
	sig_temp_MXs3.AddHistoSys("Kstarff2_uncer", "Signal_Kstarff2_m", fname_MXs3, "", "Signal_Kstarff2_p", fname_MXs3, "");
	sig_temp_MXs3.AddHistoSys("Kstarff3_uncer", "Signal_Kstarff3_m", fname_MXs3, "", "Signal_Kstarff3_p", fname_MXs3, "");
	sig_temp_MXs3.AddHistoSys("Kstarff4_uncer", "Signal_Kstarff4_m", fname_MXs3, "", "Signal_Kstarff4_p", fname_MXs3, "");
	sig_temp_MXs3.AddHistoSys("Kstarff5_uncer", "Signal_Kstarff5_m", fname_MXs3, "", "Signal_Kstarff5_p", fname_MXs3, "");
	sig_temp_MXs3.AddHistoSys("Kstarff6_uncer", "Signal_Kstarff6_m", fname_MXs3, "", "Signal_Kstarff6_p", fname_MXs3, "");
	sig_temp_MXs3.AddHistoSys("Kstarff7_uncer", "Signal_Kstarff7_m", fname_MXs3, "", "Signal_Kstarff7_p", fname_MXs3, "");
	sig_temp_MXs3.AddHistoSys("Kstarff8_uncer", "Signal_Kstarff8_m", fname_MXs3, "", "Signal_Kstarff8_p", fname_MXs3, "");
	sig_temp_MXs3.AddHistoSys("Kstarff9_uncer", "Signal_Kstarff9_m", fname_MXs3, "", "Signal_Kstarff9_p", fname_MXs3, "");
	sig_temp_MXs3.AddHistoSys("Kstarff9_uncer", "Signal_Kstarff9_m", fname_MXs3, "", "Signal_Kstarff9_p", fname_MXs3, "");
	sig_temp_MXs3.AddHistoSys("Kff_OLD_uncer", "Signal_Kff_OLD_m", fname_MXs3, "", "Signal_Kff_OLD_p", fname_MXs3, "");
	for (int i = 0; i < NEntryfraction; i++) sig_temp_MXs3.AddHistoSys(("fraction" + std::to_string(i) + "_uncer").c_str(), ("Signal_fraction_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs3, "", ("Signal_fraction_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs3, "");
	for (int i = 0; i < NEntryFragmentation; i++) sig_temp_MXs3.AddHistoSys(("Xs_fragmentation" + std::to_string(i) + "_uncer").c_str(), ("Signal_Fragmentation_correlated" + std::to_string(i) + "_m").c_str(), fname_MXs3, "", ("Signal_Fragmentation_correlated" + std::to_string(i) + "_p").c_str(), fname_MXs3, "");
	sig_temp_MXs3.AddHistoSys("pf_uncer", "Signal_pf_m", fname_MXs3, "", "Signal_pf_p", fname_MXs3, "");
	sig_temp_MXs3.AddHistoSys("mb_uncer", "Signal_mb_m", fname_MXs3, "", "Signal_mb_p", fname_MXs3, "");
	sig_temp_MXs3.AddHistoSys("transition_uncer", "Signal_transition_m", fname_MXs3, "", "Signal_transition_p", fname_MXs3, "");
	sig_temp_MXs3.AddHistoSys("mKstar_uncer", "Signal_mKstar_m", fname_MXs3, "", "Signal_mKstar_p", fname_MXs3, "");
	sig_temp_MXs3.AddNormFactor("FBDT_efficiency_CAL", 0.8605, 0.8605, 0.8605, true);
	sig_temp_MXs3.AddOverallSys("FBDT_efficiency_uncer", 0.9262, 1.0737);
	sig_temp_MXs3.AddOverallSys("BB_counting_uncer", 0.9855, 1.0145);
	sig_temp_MXs3.ActivateStatError("Signal_all_uncorrelated_MC_stat", fname_MXs3, "");
	sig_temp_MXs3.SetNormalizeByTheory(kFALSE);
	sig_temp_MXs3.AddNormFactor("mu_MXs3", expmu, -100.0, 100.0);
	channel.AddSample(sig_temp_MXs3);

	// read background template
	/* ================================ CHG ================================ */
	RooStats::HistFactory::Sample CHG_temp("CHG_nominal", "CHG_nominal", fname);
	CHG_temp.AddHistoSys("track_eff_uncer", "CHG_track_m", fname, "", "CHG_track_p", fname, "");
	CHG_temp.AddHistoSys("KS0_reco_uncer", "CHG_KS0_m", fname, "", "CHG_KS0_p", fname, "");
	CHG_temp.AddHistoSys("BtoDtoXKL_uncer", "CHG_BtoDtoXKL_m", fname, "", "CHG_BtoDtoXKL_p", fname, "");
	CHG_temp.AddHistoSys("BRBtoXKLKL_uncer", "CHG_BRBtoXKLKL_m", fname, "", "CHG_BRBtoXKLKL_p", fname, "");
	CHG_temp.AddHistoSys("Xnn_BR_uncer", "CHG_Xnn_m", fname, "", "CHG_Xnn_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) CHG_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("CHG_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryFEI; i++) CHG_temp.AddHistoSys(("FEI" + std::to_string(i) + "_uncer").c_str(), ("CHG_FEI_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_FEI_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) CHG_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("CHG_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryPID; i++) CHG_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("CHG_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryBR; i++) CHG_temp.AddHistoSys(("BBBR" + std::to_string(i) + "_uncer").c_str(), ("CHG_BR_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_BR_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntrypi0; i++) CHG_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("CHG_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	//CHG_temp.AddHistoSys("BDTc_shape_BB", "CHG_BDTc_m", fname, "", "CHG_BDTc_p", fname, "");
	CHG_temp.AddOverallSys("mu_CHG", 0.7, 1.3);
	CHG_temp.AddOverallSys("BB_counting_uncer", 0.9855, 1.0145);
	CHG_temp.ActivateStatError("CHG_all_uncorrelated_MC_stat", fname, "");
	CHG_temp.SetNormalizeByTheory(kFALSE);
	channel.AddSample(CHG_temp);
	/* ================================ CHG ================================ */

	/* ================================ MIX ================================ */
	RooStats::HistFactory::Sample MIX_temp("MIX_nominal", "MIX_nominal", fname);
	MIX_temp.AddHistoSys("track_eff_uncer", "MIX_track_m", fname, "", "MIX_track_p", fname, "");
	MIX_temp.AddHistoSys("KS0_reco_uncer", "MIX_KS0_m", fname, "", "MIX_KS0_p", fname, "");
	MIX_temp.AddHistoSys("BtoDtoXKL_uncer", "MIX_BtoDtoXKL_m", fname, "", "MIX_BtoDtoXKL_p", fname, "");
	MIX_temp.AddHistoSys("BRBtoXKLKL_uncer", "MIX_BRBtoXKLKL_m", fname, "", "MIX_BRBtoXKLKL_p", fname, "");
	MIX_temp.AddHistoSys("Xnn_BR_uncer", "MIX_Xnn_m", fname, "", "MIX_Xnn_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) MIX_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("MIX_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryFEI; i++) MIX_temp.AddHistoSys(("FEI" + std::to_string(i) + "_uncer").c_str(), ("MIX_FEI_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_FEI_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) MIX_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("MIX_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryPID; i++) MIX_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("MIX_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryBR; i++) MIX_temp.AddHistoSys(("BBBR" + std::to_string(i) + "_uncer").c_str(), ("MIX_BR_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_BR_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntrypi0; i++) MIX_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("MIX_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	//MIX_temp.AddHistoSys("BDTc_shape_BB", "MIX_BDTc_m", fname, "", "MIX_BDTc_p", fname, "");
	MIX_temp.AddOverallSys("mu_MIX", 0.7, 1.3);
	MIX_temp.AddOverallSys("BB_counting_uncer", 0.9855, 1.0145);
	MIX_temp.ActivateStatError("MIX_all_uncorrelated_MC_stat", fname, "");
	MIX_temp.SetNormalizeByTheory(kFALSE);
	channel.AddSample(MIX_temp);
	/* ================================ MIX ================================ */

	/* ================================ UUBAR ================================ */
	RooStats::HistFactory::Sample UUBAR_temp("UUBAR_nominal", "UUBAR_nominal", fname);
	UUBAR_temp.AddHistoSys("track_eff_uncer", "UUBAR_track_m", fname, "", "UUBAR_track_p", fname, "");
	UUBAR_temp.AddHistoSys("KS0_reco_uncer", "UUBAR_KS0_m", fname, "", "UUBAR_KS0_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) UUBAR_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("UUBAR_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("UUBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) UUBAR_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("UUBAR_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("UUBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryPID; i++) UUBAR_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("UUBAR_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("UUBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntrypi0; i++) UUBAR_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("UUBAR_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("UUBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	UUBAR_temp.AddHistoSys("BDTc_shape_qq", "UUBAR_BDTc_m", fname, "", "UUBAR_BDTc_p", fname, "");
	UUBAR_temp.AddNormFactor("qq_CAL_UUBAR", 1.1052, 1.1052, 1.1052, true);
	UUBAR_temp.AddOverallSys("qq_CAL_UUBAR_uncer", 0.7932, 1.2068);
	UUBAR_temp.AddOverallSys("mu_UUBAR", 0.7, 1.3);
	UUBAR_temp.ActivateStatError("UUBAR_all_uncorrelated_MC_stat", fname, "");
	UUBAR_temp.SetNormalizeByTheory(kFALSE);
	channel.AddSample(UUBAR_temp);
	/* ================================ UUBAR ================================ */

	/* ================================ DDBAR ================================ */
	RooStats::HistFactory::Sample DDBAR_temp("DDBAR_nominal", "DDBAR_nominal", fname);
	DDBAR_temp.AddHistoSys("track_eff_uncer", "DDBAR_track_m", fname, "", "DDBAR_track_p", fname, "");
	DDBAR_temp.AddHistoSys("KS0_reco_uncer", "DDBAR_KS0_m", fname, "", "DDBAR_KS0_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) DDBAR_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("DDBAR_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("DDBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) DDBAR_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("DDBAR_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("DDBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryPID; i++) DDBAR_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("DDBAR_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("DDBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntrypi0; i++) DDBAR_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("DDBAR_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("DDBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	DDBAR_temp.AddHistoSys("BDTc_shape_qq", "DDBAR_BDTc_m", fname, "", "DDBAR_BDTc_p", fname, "");
	DDBAR_temp.AddNormFactor("qq_CAL_DDBAR", 1.1052, 1.1052, 1.1052, true);
	DDBAR_temp.AddOverallSys("qq_CAL_DDBAR_uncer", 0.7932, 1.2068);
	DDBAR_temp.AddOverallSys("mu_DDBAR", 0.7, 1.3);
	DDBAR_temp.ActivateStatError("DDBAR_all_uncorrelated_MC_stat", fname, "");
	DDBAR_temp.SetNormalizeByTheory(kFALSE);
	channel.AddSample(DDBAR_temp);
	/* ================================ DDBAR ================================ */

	/* ================================ SSBAR ================================ */
	RooStats::HistFactory::Sample SSBAR_temp("SSBAR_nominal", "SSBAR_nominal", fname);
	SSBAR_temp.AddHistoSys("track_eff_uncer", "SSBAR_track_m", fname, "", "SSBAR_track_p", fname, "");
	SSBAR_temp.AddHistoSys("KS0_reco_uncer", "SSBAR_KS0_m", fname, "", "SSBAR_KS0_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) SSBAR_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("SSBAR_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("SSBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) SSBAR_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("SSBAR_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("SSBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryPID; i++) SSBAR_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("SSBAR_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("SSBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntrypi0; i++) SSBAR_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("SSBAR_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("SSBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	SSBAR_temp.AddHistoSys("BDTc_shape_qq", "SSBAR_BDTc_m", fname, "", "SSBAR_BDTc_p", fname, "");
	SSBAR_temp.AddNormFactor("qq_CAL_SSBAR", 1.1052, 1.1052, 1.1052, true);
	SSBAR_temp.AddOverallSys("qq_CAL_SSBAR_uncer", 0.7932, 1.2068);
	SSBAR_temp.AddOverallSys("mu_SSBAR", 0.7, 1.3);
	SSBAR_temp.ActivateStatError("SSBAR_all_uncorrelated_MC_stat", fname, "");
	SSBAR_temp.SetNormalizeByTheory(kFALSE);
	channel.AddSample(SSBAR_temp);
	/* ================================ SSBAR ================================ */

	/* ================================ CHARM ================================ */
	RooStats::HistFactory::Sample CHARM_temp("CHARM_nominal", "CHARM_nominal", fname);
	CHARM_temp.AddHistoSys("track_eff_uncer", "CHARM_track_m", fname, "", "CHARM_track_p", fname, "");
	CHARM_temp.AddHistoSys("KS0_reco_uncer", "CHARM_KS0_m", fname, "", "CHARM_KS0_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) CHARM_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("CHARM_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHARM_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) CHARM_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("CHARM_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHARM_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryPID; i++) CHARM_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("CHARM_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHARM_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntrypi0; i++) CHARM_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("CHARM_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHARM_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	CHARM_temp.AddHistoSys("BDTc_shape_qq", "CHARM_BDTc_m", fname, "", "CHARM_BDTc_p", fname, "");
	CHARM_temp.AddNormFactor("qq_CAL_CHARM", 1.1052, 1.1052, 1.1052, true);
	CHARM_temp.AddOverallSys("qq_CAL_CHARM_uncer", 0.7932, 1.2068);
	CHARM_temp.AddOverallSys("mu_CHARM", 0.7, 1.3);
	CHARM_temp.ActivateStatError("CHARM_all_uncorrelated_MC_stat", fname, "");
	CHARM_temp.SetNormalizeByTheory(kFALSE);
	channel.AddSample(CHARM_temp);
	/* ================================ CHARM ================================ */

	// add channel to measurement
	meas.AddChannel(channel);
	meas.CollectHistograms();

	RooWorkspace* w;
	w = RooStats::HistFactory::MakeModelAndMeasurementFast(meas);

	w->Print();
	w->writeToFile("PDFandDATA_workspace.root");

	return 0;
}