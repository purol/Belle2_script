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

int WorkSpace_Knunu() {

	int NEntryFEI = ReadNFEIEigenVector("./FEI_selected.txt");
	int NEntryKID = ReadNPIDEigenVector("./KID_selected.txt");
	int NEntryBR = ReadNBREigenVector("./BR_selected.txt");
	int NEntryMultiplicity = ReadMultiplicityInfo("./multiplicity_selected.txt");
	int NEntryfraction = ReadNfractionEigenVector("./fraction_selected.txt");
	int NEntryFragmentation = ReadNFragmentationEigenVector("./Fragmentation_selected.txt");

	const double expmu = 1.0;
    const char* fname = "PDFandDATA_nominal.root";

	RooStats::HistFactory::Measurement meas("my_measurement", "my measurement");
	meas.SetOutputFilePrefix("results/my_measurement");
	meas.SetExportOnly(kTRUE);

	// setting measurement
	meas.SetPOI("mu");
	meas.SetLumi(1.0);
	meas.AddConstantParam("Lumi");

	// define channel
	HistFactory::Channel channel("channel");
	channel.SetStatErrorConfig(1e-5, "Gaussian");

	// point data
	channel.SetData("total_DATA", fname);

	// read signal template
    /* ================================ SIGNAL ================================ */
	RooStats::HistFactory::Sample sig_temp("Signal_nominal_MXs1", "Signal_nominal", fname);
	sig_temp.AddHistoSys("track_eff_uncer", "Signal_track_m", fname, "", "Signal_track_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) sig_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("Signal_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryFEI; i++) sig_temp.AddHistoSys(("FEI" + std::to_string(i) + "_uncer").c_str(), ("Signal_FEI_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_FEI_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) sig_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("Signal_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryBR; i++) sig_temp.AddHistoSys(("BBBR" + std::to_string(i) + "_uncer").c_str(), ("Signal_BR_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_BR_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	sig_temp.AddHistoSys("Kff1_uncer", "Signal_Kff1_m", fname, "", "Signal_Kff1_p", fname, "");
	sig_temp.AddHistoSys("Kff2_uncer", "Signal_Kff2_m", fname, "", "Signal_Kff2_p", fname, "");
	sig_temp.AddHistoSys("Kff3_uncer", "Signal_Kff3_m", fname, "", "Signal_Kff3_p", fname, "");
	sig_temp.AddHistoSys("Kff_OLD_uncer", "Signal_Kff_OLD_m", fname, "", "Signal_Kff_OLD_p", fname, "");
	for (int i = 0; i < NEntryfraction; i++) sig_temp.AddHistoSys(("fraction" + std::to_string(i) + "_uncer").c_str(), ("Signal_fraction_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_fraction_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	sig_temp.AddNormFactor("FBDT_efficiency_CAL", 0.8605, 0.8605, 0.8605, true);
	sig_temp.AddOverallSys("FBDT_efficiency_uncer", 0.9262, 1.0737);
	sig_temp.AddOverallSys("BB_counting_uncer", 0.9855, 1.0145);
	sig_temp.ActivateStatError("Signal_all_uncorrelated_MC_stat", fname, "");
	sig_temp.SetNormalizeByTheory(kFALSE);
	sig_temp.AddNormFactor("mu", expmu, -100.0, 100.0);
	channel.AddSample(sig_temp);

	// read background template
	/* ================================ CHG ================================ */
	RooStats::HistFactory::Sample CHG_temp("CHG_nominal", "CHG_nominal", fname);
	CHG_temp.AddHistoSys("track_eff_uncer", "CHG_track_m", fname, "", "CHG_track_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) CHG_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("CHG_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryFEI; i++) CHG_temp.AddHistoSys(("FEI" + std::to_string(i) + "_uncer").c_str(), ("CHG_FEI_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_FEI_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) CHG_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("CHG_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryBR; i++) CHG_temp.AddHistoSys(("BBBR" + std::to_string(i) + "_uncer").c_str(), ("CHG_BR_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_BR_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	CHG_temp.AddHistoSys("Kff1_uncer", "CHG_Kff1_m", fname, "", "CHG_Kff1_p", fname, "");
	CHG_temp.AddHistoSys("Kff2_uncer", "CHG_Kff2_m", fname, "", "CHG_Kff2_p", fname, "");
	CHG_temp.AddHistoSys("Kff3_uncer", "CHG_Kff3_m", fname, "", "CHG_Kff3_p", fname, "");
	CHG_temp.AddHistoSys("Kstarff1_uncer", "CHG_Kstarff1_m", fname, "", "CHG_Kstarff1_p", fname, "");
	CHG_temp.AddHistoSys("Kstarff2_uncer", "CHG_Kstarff2_m", fname, "", "CHG_Kstarff2_p", fname, "");
	CHG_temp.AddHistoSys("Kstarff3_uncer", "CHG_Kstarff3_m", fname, "", "CHG_Kstarff3_p", fname, "");
	CHG_temp.AddHistoSys("Kstarff4_uncer", "CHG_Kstarff4_m", fname, "", "CHG_Kstarff4_p", fname, "");
	CHG_temp.AddHistoSys("Kstarff5_uncer", "CHG_Kstarff5_m", fname, "", "CHG_Kstarff5_p", fname, "");
	CHG_temp.AddHistoSys("Kstarff6_uncer", "CHG_Kstarff6_m", fname, "", "CHG_Kstarff6_p", fname, "");
	CHG_temp.AddHistoSys("Kstarff7_uncer", "CHG_Kstarff7_m", fname, "", "CHG_Kstarff7_p", fname, "");
	CHG_temp.AddHistoSys("Kstarff8_uncer", "CHG_Kstarff8_m", fname, "", "CHG_Kstarff8_p", fname, "");
	CHG_temp.AddHistoSys("Kstarff9_uncer", "CHG_Kstarff9_m", fname, "", "CHG_Kstarff9_p", fname, "");
	CHG_temp.AddHistoSys("Kstarff9_uncer", "CHG_Kstarff9_m", fname, "", "CHG_Kstarff9_p", fname, "");
	CHG_temp.AddHistoSys("Kff_OLD_uncer", "CHG_Kff_OLD_m", fname, "", "CHG_Kff_OLD_p", fname, "");
	for (int i = 0; i < NEntryfraction; i++) CHG_temp.AddHistoSys(("fraction" + std::to_string(i) + "_uncer").c_str(), ("CHG_fraction_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_fraction_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryFragmentation; i++) CHG_temp.AddHistoSys(("Xs_fragmentation" + std::to_string(i) + "_uncer").c_str(), ("CHG_Fragmentation_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_Fragmentation_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	CHG_temp.AddHistoSys("pf_uncer", "CHG_pf_m", fname, "", "CHG_pf_p", fname, "");
	CHG_temp.AddHistoSys("mb_uncer", "CHG_mb_m", fname, "", "CHG_mb_p", fname, "");
	CHG_temp.AddHistoSys("transition_uncer", "CHG_transition_m", fname, "", "CHG_transition_p", fname, "");
	CHG_temp.AddHistoSys("mKstar_uncer", "CHG_mKstar_m", fname, "", "CHG_mKstar_p", fname, "");
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
	for (int i = 0; i < NEntryMultiplicity; i++) MIX_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("MIX_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryFEI; i++) MIX_temp.AddHistoSys(("FEI" + std::to_string(i) + "_uncer").c_str(), ("MIX_FEI_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_FEI_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) MIX_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("MIX_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryBR; i++) MIX_temp.AddHistoSys(("BBBR" + std::to_string(i) + "_uncer").c_str(), ("MIX_BR_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_BR_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	MIX_temp.AddHistoSys("Kff1_uncer", "MIX_Kff1_m", fname, "", "MIX_Kff1_p", fname, "");
	MIX_temp.AddHistoSys("Kff2_uncer", "MIX_Kff2_m", fname, "", "MIX_Kff2_p", fname, "");
	MIX_temp.AddHistoSys("Kff3_uncer", "MIX_Kff3_m", fname, "", "MIX_Kff3_p", fname, "");
	MIX_temp.AddHistoSys("Kstarff1_uncer", "MIX_Kstarff1_m", fname, "", "MIX_Kstarff1_p", fname, "");
	MIX_temp.AddHistoSys("Kstarff2_uncer", "MIX_Kstarff2_m", fname, "", "MIX_Kstarff2_p", fname, "");
	MIX_temp.AddHistoSys("Kstarff3_uncer", "MIX_Kstarff3_m", fname, "", "MIX_Kstarff3_p", fname, "");
	MIX_temp.AddHistoSys("Kstarff4_uncer", "MIX_Kstarff4_m", fname, "", "MIX_Kstarff4_p", fname, "");
	MIX_temp.AddHistoSys("Kstarff5_uncer", "MIX_Kstarff5_m", fname, "", "MIX_Kstarff5_p", fname, "");
	MIX_temp.AddHistoSys("Kstarff6_uncer", "MIX_Kstarff6_m", fname, "", "MIX_Kstarff6_p", fname, "");
	MIX_temp.AddHistoSys("Kstarff7_uncer", "MIX_Kstarff7_m", fname, "", "MIX_Kstarff7_p", fname, "");
	MIX_temp.AddHistoSys("Kstarff8_uncer", "MIX_Kstarff8_m", fname, "", "MIX_Kstarff8_p", fname, "");
	MIX_temp.AddHistoSys("Kstarff9_uncer", "MIX_Kstarff9_m", fname, "", "MIX_Kstarff9_p", fname, "");
	MIX_temp.AddHistoSys("Kstarff9_uncer", "MIX_Kstarff9_m", fname, "", "MIX_Kstarff9_p", fname, "");
	MIX_temp.AddHistoSys("Kff_OLD_uncer", "MIX_Kff_OLD_m", fname, "", "MIX_Kff_OLD_p", fname, "");
	for (int i = 0; i < NEntryfraction; i++) MIX_temp.AddHistoSys(("fraction" + std::to_string(i) + "_uncer").c_str(), ("MIX_fraction_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_fraction_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryFragmentation; i++) MIX_temp.AddHistoSys(("Xs_fragmentation" + std::to_string(i) + "_uncer").c_str(), ("MIX_Fragmentation_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_Fragmentation_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	MIX_temp.AddHistoSys("pf_uncer", "MIX_pf_m", fname, "", "MIX_pf_p", fname, "");
	MIX_temp.AddHistoSys("mb_uncer", "MIX_mb_m", fname, "", "MIX_mb_p", fname, "");
	MIX_temp.AddHistoSys("transition_uncer", "MIX_transition_m", fname, "", "MIX_transition_p", fname, "");
	MIX_temp.AddHistoSys("mKstar_uncer", "MIX_mKstar_m", fname, "", "MIX_mKstar_p", fname, "");
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
	for (int i = 0; i < NEntryMultiplicity; i++) UUBAR_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("UUBAR_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("UUBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) UUBAR_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("UUBAR_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("UUBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
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
	for (int i = 0; i < NEntryMultiplicity; i++) DDBAR_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("DDBAR_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("DDBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) DDBAR_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("DDBAR_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("DDBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
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
	for (int i = 0; i < NEntryMultiplicity; i++) SSBAR_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("SSBAR_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("SSBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) SSBAR_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("SSBAR_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("SSBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
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
	for (int i = 0; i < NEntryMultiplicity; i++) CHARM_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("CHARM_multiplicity_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHARM_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryKID; i++) CHARM_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("CHARM_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHARM_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
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
