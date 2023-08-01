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

// arXiv:1409.4557v2, PhysRevD.107.014511
# define TB0 1.5195 // (Table. 1)
# define TBp 1.6384 // (Table. 1)
# define BR_Kplus_nunubar 0.000005044 // Table VI = (5.044 +- 0.402) * 10^{-6}
# define BR_K0star_nunubar 0.00000919 // (eq. 11)
# define BR_K0_nunubar (BR_Kplus_nunubar*TB0/TBp) // under (eq. 15). In Table VI = (4.6669 +- 0.3707) * 10^{-6}
# define BR_Kplusstar_nunubar (BR_K0star_nunubar*TBp/TB0) // under (eq. 15)
# define BR_Xs_nunubar 0.000029 // (eq. 23)
# define BR_Xsu_nonresonant_nunubar (BR_Xs_nunubar - BR_Kplus_nunubar - BR_Kplusstar_nunubar)
# define BR_Xsd_nonresonant_nunubar (BR_Xs_nunubar - BR_K0_nunubar - BR_K0star_nunubar)

// according to DIRAC
# define N_BpBp_1invab 540000000.0
# define N_B0B0_1invab 510000000.0

# define BR_BpBp 0.514
# define BR_B0B0 0.486

// https://confluence.desy.de/pages/viewpage.action?spaceKey=BI&title=Conference+readiness
# define N_BB_LS1 387100000.0 // NBB = (387.1 +/- 5.6) x 10^6

# define N_Kplus_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_BpBp/(BR_BpBp+BR_B0B0)) * BR_Kplus_nunubar)
# define N_Kplusstar_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_BpBp/(BR_BpBp+BR_B0B0)) * BR_Kplusstar_nunubar)
# define N_Xsu_nonresonant_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_BpBp/(BR_BpBp+BR_B0B0)) * BR_Xsu_nonresonant_nunubar)
# define N_K0_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_B0B0/(BR_BpBp+BR_B0B0)) * BR_K0_nunubar)
# define N_K0star_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_B0B0/(BR_BpBp+BR_B0B0)) * BR_K0star_nunubar)
# define N_Xsd_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_B0B0/(BR_BpBp+BR_B0B0)) * BR_Xsd_nonresonant_nunubar)

// SIGNAL MC sample number before skimming
# define N_Kplus_train 7039000.0
# define N_K0_train 7166624.0
# define N_Kplusstar_train 7039000.0
# define N_K0star_train 7166624.0
# define N_Xsu_nonresonant_train 35195000.0
# define N_Xsd_nonresonant_train 34940430.0
# define N_Kplus_test 2961000.0
# define N_K0_test 2833376.0
# define N_Kplusstar_test 2961000.0
# define N_K0star_test 2833376.0
# define N_Xsu_nonresonant_test 14805000.0
# define N_Xsd_nonresonant_test 15059570.0

// scale factor for SIGNAL MC sample until LS1
# define Scale_Kplus_train (N_Kplus_nunubar_LS1/N_Kplus_train)
# define Scale_Kplusstar_train (N_Kplusstar_nunubar_LS1/N_Kplusstar_train)
# define Scale_Xsu_nonresonant_train (N_Xsu_nonresonant_nunubar_LS1/N_Xsu_nonresonant_train)
# define Scale_K0_train (N_K0_nunubar_LS1/N_K0_train)
# define Scale_K0star_train (N_K0star_nunubar_LS1/N_K0star_train)
# define Scale_Xsd_nonresonant_train (N_Xsd_nunubar_LS1/N_Xsd_nonresonant_train)
# define Scale_Kplus_test (N_Kplus_nunubar_LS1/N_Kplus_test)
# define Scale_Kplusstar_test (N_Kplusstar_nunubar_LS1/N_Kplusstar_test)
# define Scale_Xsu_nonresonant_test (N_Xsu_nonresonant_nunubar_LS1/N_Xsu_nonresonant_test)
# define Scale_K0_test (N_K0_nunubar_LS1/N_K0_test)
# define Scale_K0star_test (N_K0star_nunubar_LS1/N_K0star_test)
# define Scale_Xsd_nonresonant_test (N_Xsd_nunubar_LS1/N_Xsd_nonresonant_test)

// BKG MC sample number (2.8/ab for BB, 1.0/ab for qq)
# define N_CHG_test 139768443.0
# define N_MIX_test 104591690.0
# define N_UUBAR_test 107641351.0
# define N_DDBAR_test 25588403.0
# define N_SSBAR_test 20668794.0
# define N_CHARM_test 116108850.0
# define N_CHG_train 151915459.0
# define N_MIX_train 120136353.0
# define N_UUBAR_train 158509639.0
# define N_DDBAR_train 38644413.0
# define N_SSBAR_train 31001866.0
# define N_CHARM_train 174901296.0

// new scale factor for BKG MC sample with additional 1/ab (364.436 - 2.763 = 361.673/fb), until LS1
# define Scale_CHG_train ((N_BB_LS1* (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab * (N_CHG_train / (N_CHG_train + N_CHG_test)) ))
# define Scale_MIX_train ((N_BB_LS1* (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab * (N_MIX_train / (N_MIX_train + N_MIX_test)) ))
# define Scale_UUBAR_train (0.361673/((N_UUBAR_train/(N_UUBAR_train + N_UUBAR_test))*1.0))
# define Scale_DDBAR_train (0.361673/((N_DDBAR_train/(N_DDBAR_train + N_DDBAR_test))*1.0))
# define Scale_SSBAR_train (0.361673/((N_SSBAR_train/(N_SSBAR_train + N_SSBAR_test))*1.0))
# define Scale_CHARM_train (0.361673/((N_CHARM_train/(N_CHARM_train + N_CHARM_test))*1.0))
# define Scale_CHG_test ((N_BB_LS1* (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab * (N_CHG_test / (N_CHG_train + N_CHG_test)) ))
# define Scale_MIX_test ((N_BB_LS1* (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab * (N_MIX_test / (N_MIX_train + N_MIX_test)) ))
# define Scale_UUBAR_test (0.361673/((N_UUBAR_test/(N_UUBAR_train + N_UUBAR_test))*1.0))
# define Scale_DDBAR_test (0.361673/((N_DDBAR_test/(N_DDBAR_train + N_DDBAR_test))*1.0))
# define Scale_SSBAR_test (0.361673/((N_SSBAR_test/(N_SSBAR_train + N_SSBAR_test))*1.0))
# define Scale_CHARM_test (0.361673/((N_CHARM_test/(N_CHARM_train + N_CHARM_test))*1.0))

// several constant for systematic uncertainty
# define KS0_rel_uncertainty 0.6 // %/cm
# define track_rel_uncertainty 0.69 // %
// # define pi0_correction 0.932
// # define pi0_rel_uncertainty ((0.0369 / 0.932) * 100.0) // %
# define FEI_cal_Bc 0.679
# define FEI_cal_Bc_uncertainty (0.017/FEI_cal_Bc) // not percentage. relative uncertainty
# define FEI_cal_B0 0.713
# define FEI_cal_B0_uncertainty (0.019/FEI_cal_B0) // not percentage. relative uncertainty

# define RarityBins 8

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

int WorkSpace() {

	int NEntryKID = ReadNPIDEigenVector("./KID_selected.txt");
	int NEntryPID = ReadNPIDEigenVector("./PID_selected.txt");
	int NEntryBR = ReadNBREigenVector("./BR_selected.txt");
	int NEntrypi0 = ReadNpi0EigenVector("./pi0_selected.txt");

	const double expmu = 1.0;
	const char* fname = "PDFandDATA.root";

	RooStats::HistFactory::Measurement meas("my_measurement", "my measurement");
	meas.SetOutputFilePrefix("results/my_measurement");
	meas.SetExportOnly(kTRUE);

	// setting measurement
	meas.SetPOI("mu");
	meas.SetLumi(1.0);
	meas.AddConstantParam("Lumi");

	// define channel
	HistFactory::Channel channel("channel");
	channel.SetStatErrorConfig(1e-5, "Poisson");

	// point data
	channel.SetData("total_DATA", fname);

	// read signal template
	RooStats::HistFactory::Sample sig_temp("Signal_nominal", "Signal_nominal", fname);
	sig_temp.AddHistoSys("FEI_charged_uncer", "Signal_FEI_charged_m", fname, "", "Signal_FEI_charged_p", fname, "");
	sig_temp.AddHistoSys("FEI_neutral_uncer", "Signal_FEI_neutral_m", fname, "", "Signal_FEI_neutral_p", fname, "");
	sig_temp.AddHistoSys("track_eff_uncer", "Signal_track_m", fname, "", "Signal_track_p", fname, "");
	sig_temp.AddHistoSys("KS0_reco_uncer", "Signal_KS0_m", fname, "", "Signal_KS0_p", fname, "");
	sig_temp.AddHistoSys("Eecl_uncer", "Signal_Eecl_m", fname, "", "Signal_Eecl_p", fname, "");
	for (int i = 0; i < NEntryKID; i++) sig_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("Signal_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryPID; i++) sig_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("Signal_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryBR; i++) sig_temp.AddHistoSys(("BBBR" + std::to_string(i) + "_uncer").c_str(), ("Signal_BR_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_BR_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntrypi0; i++) sig_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("Signal_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	sig_temp.AddShapeSys("Signal_all_uncorr_uncer", RooStats::HistFactory::Constraint::Type::Gaussian, "Signal_all_uncorrelated", fname, "");
	sig_temp.AddHistoSys("Kff1_uncer", "Signal_Kff1_m", fname, "", "Signal_Kff1_p", fname, "");
	sig_temp.AddHistoSys("Kff2_uncer", "Signal_Kff2_m", fname, "", "Signal_Kff2_p", fname, "");
	sig_temp.AddHistoSys("Kff3_uncer", "Signal_Kff3_m", fname, "", "Signal_Kff3_p", fname, "");
	sig_temp.AddHistoSys("Kstarff1_uncer", "Signal_Kstarff1_m", fname, "", "Signal_Kstarff1_p", fname, "");
	sig_temp.AddHistoSys("Kstarff2_uncer", "Signal_Kstarff2_m", fname, "", "Signal_Kstarff2_p", fname, "");
	sig_temp.AddHistoSys("Kstarff3_uncer", "Signal_Kstarff3_m", fname, "", "Signal_Kstarff3_p", fname, "");
	sig_temp.AddHistoSys("Kstarff4_uncer", "Signal_Kstarff4_m", fname, "", "Signal_Kstarff4_p", fname, "");
	sig_temp.AddHistoSys("Kstarff5_uncer", "Signal_Kstarff5_m", fname, "", "Signal_Kstarff5_p", fname, "");
	sig_temp.AddHistoSys("Kstarff6_uncer", "Signal_Kstarff6_m", fname, "", "Signal_Kstarff6_p", fname, "");
	sig_temp.AddHistoSys("Kstarff7_uncer", "Signal_Kstarff7_m", fname, "", "Signal_Kstarff7_p", fname, "");
	sig_temp.AddHistoSys("Kstarff8_uncer", "Signal_Kstarff8_m", fname, "", "Signal_Kstarff8_p", fname, "");
	sig_temp.AddHistoSys("Kstarff9_uncer", "Signal_Kstarff9_m", fname, "", "Signal_Kstarff9_p", fname, "");
	sig_temp.AddHistoSys("Kfrac_uncer", "Signal_Kfrac_m", fname, "", "Signal_Kfrac_p", fname, "");
	sig_temp.AddHistoSys("Kstarfrac_uncer", "Signal_Kstarfrac_m", fname, "", "Signal_Kstarfrac_p", fname, "");
	sig_temp.AddHistoSys("Xsu_fragmentation1_uncer", "Signal_Xsu_frag_decay1_m", fname, "", "Signal_Xsu_frag_decay1_p", fname, "");
	sig_temp.AddHistoSys("Xsu_fragmentation2_uncer", "Signal_Xsu_frag_decay2_m", fname, "", "Signal_Xsu_frag_decay2_p", fname, "");
	sig_temp.AddHistoSys("Xsu_fragmentation3_uncer", "Signal_Xsu_frag_decay3_m", fname, "", "Signal_Xsu_frag_decay3_p", fname, "");
	sig_temp.AddHistoSys("Xsu_fragmentation4_uncer", "Signal_Xsu_frag_decay4_m", fname, "", "Signal_Xsu_frag_decay4_p", fname, "");
	sig_temp.AddHistoSys("Xsu_fragmentation5_uncer", "Signal_Xsu_frag_decay5_m", fname, "", "Signal_Xsu_frag_decay5_p", fname, "");
	sig_temp.AddHistoSys("Xsu_fragmentation6_uncer", "Signal_Xsu_frag_decay6_m", fname, "", "Signal_Xsu_frag_decay6_p", fname, "");
	sig_temp.AddHistoSys("Xsu_fragmentation7_uncer", "Signal_Xsu_frag_decay7_m", fname, "", "Signal_Xsu_frag_decay7_p", fname, "");
	sig_temp.AddHistoSys("Xsu_fragmentation8_uncer", "Signal_Xsu_frag_decay8_m", fname, "", "Signal_Xsu_frag_decay8_p", fname, "");
	sig_temp.AddHistoSys("Xsu_fragmentation9_uncer", "Signal_Xsu_frag_decay9_m", fname, "", "Signal_Xsu_frag_decay9_p", fname, "");
	sig_temp.AddHistoSys("Xsu_fragmentation10_uncer", "Signal_Xsu_frag_decay10_m", fname, "", "Signal_Xsu_frag_decay10_p", fname, "");
	sig_temp.AddHistoSys("Xsu_fragmentation11_uncer", "Signal_Xsu_frag_decay11_m", fname, "", "Signal_Xsu_frag_decay11_p", fname, "");
	sig_temp.AddHistoSys("Xsu_fragmentation12_uncer", "Signal_Xsu_frag_decay12_m", fname, "", "Signal_Xsu_frag_decay12_p", fname, "");
	sig_temp.AddHistoSys("Xsu_fragmentation13_uncer", "Signal_Xsu_frag_decay13_m", fname, "", "Signal_Xsu_frag_decay13_p", fname, "");
	sig_temp.AddHistoSys("Xsu_fragmentation14_uncer", "Signal_Xsu_frag_decay14_m", fname, "", "Signal_Xsu_frag_decay14_p", fname, "");
	sig_temp.AddHistoSys("Xsu_fragmentation30_uncer", "Signal_Xsu_frag_decay30_m", fname, "", "Signal_Xsu_frag_decay30_p", fname, "");
	sig_temp.AddHistoSys("Xsd_fragmentation16_uncer", "Signal_Xsd_frag_decay16_m", fname, "", "Signal_Xsd_frag_decay16_p", fname, "");
	sig_temp.AddHistoSys("Xsd_fragmentation17_uncer", "Signal_Xsd_frag_decay17_m", fname, "", "Signal_Xsd_frag_decay17_p", fname, "");
	sig_temp.AddHistoSys("Xsd_fragmentation18_uncer", "Signal_Xsd_frag_decay18_m", fname, "", "Signal_Xsd_frag_decay18_p", fname, "");
	sig_temp.AddHistoSys("Xsd_fragmentation19_uncer", "Signal_Xsd_frag_decay19_m", fname, "", "Signal_Xsd_frag_decay19_p", fname, "");
	sig_temp.AddHistoSys("Xsd_fragmentation20_uncer", "Signal_Xsd_frag_decay20_m", fname, "", "Signal_Xsd_frag_decay20_p", fname, "");
	sig_temp.AddHistoSys("Xsd_fragmentation21_uncer", "Signal_Xsd_frag_decay21_m", fname, "", "Signal_Xsd_frag_decay21_p", fname, "");
	sig_temp.AddHistoSys("Xsd_fragmentation22_uncer", "Signal_Xsd_frag_decay22_m", fname, "", "Signal_Xsd_frag_decay22_p", fname, "");
	sig_temp.AddHistoSys("Xsd_fragmentation23_uncer", "Signal_Xsd_frag_decay23_m", fname, "", "Signal_Xsd_frag_decay23_p", fname, "");
	sig_temp.AddHistoSys("Xsd_fragmentation24_uncer", "Signal_Xsd_frag_decay24_m", fname, "", "Signal_Xsd_frag_decay24_p", fname, "");
	sig_temp.AddHistoSys("Xsd_fragmentation25_uncer", "Signal_Xsd_frag_decay25_m", fname, "", "Signal_Xsd_frag_decay25_p", fname, "");
	sig_temp.AddHistoSys("Xsd_fragmentation26_uncer", "Signal_Xsd_frag_decay26_m", fname, "", "Signal_Xsd_frag_decay26_p", fname, "");
	sig_temp.AddHistoSys("Xsd_fragmentation27_uncer", "Signal_Xsd_frag_decay27_m", fname, "", "Signal_Xsd_frag_decay27_p", fname, "");
	sig_temp.AddHistoSys("Xsd_fragmentation28_uncer", "Signal_Xsd_frag_decay28_m", fname, "", "Signal_Xsd_frag_decay28_p", fname, "");
	sig_temp.AddHistoSys("Xsd_fragmentation29_uncer", "Signal_Xsd_frag_decay29_m", fname, "", "Signal_Xsd_frag_decay29_p", fname, "");
	sig_temp.AddHistoSys("Xsd_fragmentation30_uncer", "Signal_Xsd_frag_decay30_m", fname, "", "Signal_Xsd_frag_decay30_p", fname, "");
	sig_temp.AddHistoSys("pf_uncer", "Signal_pf_m", fname, "", "Signal_pf_p", fname, "");
	sig_temp.AddHistoSys("mb_uncer", "Signal_mb_m", fname, "", "Signal_mb_p", fname, "");
	sig_temp.AddHistoSys("transition_uncer", "Signal_transition_m", fname, "", "Signal_transition_p", fname, "");
	sig_temp.AddHistoSys("mKstar_uncer", "Signal_mKstar_m", fname, "", "Signal_mKstar_p", fname, "");
	sig_temp.AddNormFactor("FBDT_efficiency_CAL", 0.8605, 0.8605, 0.8605, true);
	sig_temp.AddOverallSys("FBDT_efficiency_uncer", 0.9262, 1.0737);
	sig_temp.AddOverallSys("BB_counting_uncer", 0.9855, 1.0145);
	sig_temp.ActivateStatError();
	sig_temp.SetNormalizeByTheory(kFALSE);
	sig_temp.AddNormFactor("mu", expmu, -100.0, 100.0);
	channel.AddSample(sig_temp);

	// read background template
	/* ================================ CHG ================================ */
	RooStats::HistFactory::Sample CHG_temp("CHG_nominal", "CHG_nominal", fname);
	CHG_temp.AddOverallSys("FEI_charged_uncer", 1 - FEI_cal_Bc_uncertainty, 1 + FEI_cal_Bc_uncertainty);
	CHG_temp.AddHistoSys("track_eff_uncer", "CHG_track_m", fname, "", "CHG_track_p", fname, "");
	CHG_temp.AddHistoSys("KS0_reco_uncer", "CHG_KS0_m", fname, "", "CHG_KS0_p", fname, "");
	CHG_temp.AddHistoSys("Eecl_uncer", "CHG_Eecl_m", fname, "", "CHG_Eecl_p", fname, "");
	CHG_temp.AddHistoSys("Knn_BR_uncer", "CHG_Knn_m", fname, "", "CHG_Knn_p", fname, "");
	CHG_temp.AddHistoSys("Kstarnn_BR_uncer", "CHG_Kstarnn_m", fname, "", "CHG_Kstarnn_p", fname, "");
	for (int i = 0; i < NEntryKID; i++) CHG_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("CHG_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryPID; i++) CHG_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("CHG_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryBR; i++) CHG_temp.AddHistoSys(("BBBR" + std::to_string(i) + "_uncer").c_str(), ("CHG_BR_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_BR_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntrypi0; i++) CHG_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("CHG_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	CHG_temp.AddShapeSys("CHG_all_uncorr_uncer", RooStats::HistFactory::Constraint::Type::Gaussian, "CHG_all_uncorrelated", fname, "");
	//CHG_temp.AddHistoSys("BDTc_shape_BB", "CHG_BDTc_m", fname, "", "CHG_BDTc_p", fname, "");
    //CHG_temp.AddNormFactor("Norm_CHG_CAL", 1.3, 1.3, 1.3, true);
	CHG_temp.AddOverallSys("mu_CHG", 0.7, 1.3);
	CHG_temp.AddOverallSys("BB_counting_uncer", 0.9855, 1.0145);
	CHG_temp.ActivateStatError();
	CHG_temp.SetNormalizeByTheory(kFALSE);
	channel.AddSample(CHG_temp);
	/* ================================ CHG ================================ */

	/* ================================ MIX ================================ */
	RooStats::HistFactory::Sample MIX_temp("MIX_nominal", "MIX_nominal", fname);
	MIX_temp.AddOverallSys("FEI_neutral_uncer", 1 - FEI_cal_B0_uncertainty, 1 + FEI_cal_B0_uncertainty);
	MIX_temp.AddHistoSys("track_eff_uncer", "MIX_track_m", fname, "", "MIX_track_p", fname, "");
	MIX_temp.AddHistoSys("KS0_reco_uncer", "MIX_KS0_m", fname, "", "MIX_KS0_p", fname, "");
	MIX_temp.AddHistoSys("Eecl_uncer", "MIX_Eecl_m", fname, "", "MIX_Eecl_p", fname, "");
	MIX_temp.AddHistoSys("K0nn_BR_uncer", "MIX_K0nn_m", fname, "", "MIX_K0nn_p", fname, "");
	MIX_temp.AddHistoSys("K0starnn_BR_uncer", "MIX_K0starnn_m", fname, "", "MIX_K0starnn_p", fname, "");
	for (int i = 0; i < NEntryKID; i++) MIX_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("MIX_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryPID; i++) MIX_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("MIX_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryBR; i++) MIX_temp.AddHistoSys(("BBBR" + std::to_string(i) + "_uncer").c_str(), ("MIX_BR_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_BR_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntrypi0; i++) MIX_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("MIX_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	MIX_temp.AddShapeSys("MIX_all_uncorr_uncer", RooStats::HistFactory::Constraint::Type::Gaussian, "MIX_all_uncorrelated", fname, "");
	//MIX_temp.AddHistoSys("BDTc_shape_BB", "MIX_BDTc_m", fname, "", "MIX_BDTc_p", fname, "");
    //MIX_temp.AddNormFactor("Norm_MIX_CAL", 1.3, 1.3, 1.3, true);
	MIX_temp.AddOverallSys("mu_MIX", 0.7, 1.3);
	MIX_temp.AddOverallSys("BB_counting_uncer", 0.9855, 1.0145);
	MIX_temp.ActivateStatError();
	MIX_temp.SetNormalizeByTheory(kFALSE);
	channel.AddSample(MIX_temp);
	/* ================================ MIX ================================ */

	/* ================================ UUBAR ================================ */
	RooStats::HistFactory::Sample UUBAR_temp("UUBAR_nominal", "UUBAR_nominal", fname);
	UUBAR_temp.AddHistoSys("track_eff_uncer", "UUBAR_track_m", fname, "", "UUBAR_track_p", fname, "");
	UUBAR_temp.AddHistoSys("KS0_reco_uncer", "UUBAR_KS0_m", fname, "", "UUBAR_KS0_p", fname, "");
	UUBAR_temp.AddHistoSys("Eecl_uncer", "UUBAR_Eecl_m", fname, "", "UUBAR_Eecl_p", fname, "");
	for (int i = 0; i < NEntryKID; i++) UUBAR_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("UUBAR_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("UUBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryPID; i++) UUBAR_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("UUBAR_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("UUBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntrypi0; i++) UUBAR_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("UUBAR_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("UUBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	UUBAR_temp.AddShapeSys("UUBAR_all_uncorr_uncer", RooStats::HistFactory::Constraint::Type::Gaussian, "UUBAR_all_uncorrelated", fname, "");
	UUBAR_temp.AddHistoSys("BDTc_shape_qq", "UUBAR_BDTc_m", fname, "", "UUBAR_BDTc_p", fname, "");
	UUBAR_temp.AddNormFactor("qq_CAL_UUBAR", 1.1052, 1.1052, 1.1052, true);
	UUBAR_temp.AddOverallSys("qq_CAL_UUBAR_uncer", 0.7932, 1.2068);
    //UUBAR_temp.AddNormFactor("Norm_UUBAR_CAL", 1.3, 1.3, 1.3, true);
	UUBAR_temp.AddOverallSys("mu_UUBAR", 0.7, 1.3);
	UUBAR_temp.ActivateStatError();
	UUBAR_temp.SetNormalizeByTheory(kFALSE);
	channel.AddSample(UUBAR_temp);
	/* ================================ UUBAR ================================ */

	/* ================================ DDBAR ================================ */
	RooStats::HistFactory::Sample DDBAR_temp("DDBAR_nominal", "DDBAR_nominal", fname);
	DDBAR_temp.AddHistoSys("track_eff_uncer", "DDBAR_track_m", fname, "", "DDBAR_track_p", fname, "");
	DDBAR_temp.AddHistoSys("KS0_reco_uncer", "DDBAR_KS0_m", fname, "", "DDBAR_KS0_p", fname, "");
	DDBAR_temp.AddHistoSys("Eecl_uncer", "DDBAR_Eecl_m", fname, "", "DDBAR_Eecl_p", fname, "");
	for (int i = 0; i < NEntryKID; i++) DDBAR_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("DDBAR_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("DDBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryPID; i++) DDBAR_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("DDBAR_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("DDBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntrypi0; i++) DDBAR_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("DDBAR_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("DDBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	DDBAR_temp.AddShapeSys("DDBAR_all_uncorr_uncer", RooStats::HistFactory::Constraint::Type::Gaussian, "DDBAR_all_uncorrelated", fname, "");
	DDBAR_temp.AddHistoSys("BDTc_shape_qq", "DDBAR_BDTc_m", fname, "", "DDBAR_BDTc_p", fname, "");
	DDBAR_temp.AddNormFactor("qq_CAL_DDBAR", 1.1052, 1.1052, 1.1052, true);
	DDBAR_temp.AddOverallSys("qq_CAL_DDBAR_uncer", 0.7932, 1.2068);
    //DDBAR_temp.AddNormFactor("Norm_DDBAR_CAL", 1.3, 1.3, 1.3, true);
	DDBAR_temp.AddOverallSys("mu_DDBAR", 0.7, 1.3);
	DDBAR_temp.ActivateStatError();
	DDBAR_temp.SetNormalizeByTheory(kFALSE);
	channel.AddSample(DDBAR_temp);
	/* ================================ DDBAR ================================ */

	/* ================================ SSBAR ================================ */
	RooStats::HistFactory::Sample SSBAR_temp("SSBAR_nominal", "SSBAR_nominal", fname);
	SSBAR_temp.AddHistoSys("track_eff_uncer", "SSBAR_track_m", fname, "", "SSBAR_track_p", fname, "");
	SSBAR_temp.AddHistoSys("KS0_reco_uncer", "SSBAR_KS0_m", fname, "", "SSBAR_KS0_p", fname, "");
	SSBAR_temp.AddHistoSys("Eecl_uncer", "SSBAR_Eecl_m", fname, "", "SSBAR_Eecl_p", fname, "");
	for (int i = 0; i < NEntryKID; i++) SSBAR_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("SSBAR_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("SSBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryPID; i++) SSBAR_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("SSBAR_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("SSBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntrypi0; i++) SSBAR_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("SSBAR_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("SSBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	SSBAR_temp.AddShapeSys("SSBAR_all_uncorr_uncer", RooStats::HistFactory::Constraint::Type::Gaussian, "SSBAR_all_uncorrelated", fname, "");
	SSBAR_temp.AddHistoSys("BDTc_shape_qq", "SSBAR_BDTc_m", fname, "", "SSBAR_BDTc_p", fname, "");
	SSBAR_temp.AddNormFactor("qq_CAL_SSBAR", 1.1052, 1.1052, 1.1052, true);
	SSBAR_temp.AddOverallSys("qq_CAL_SSBAR_uncer", 0.7932, 1.2068);
    //SSBAR_temp.AddNormFactor("Norm_SSBAR_CAL", 1.3, 1.3, 1.3, true);
	SSBAR_temp.AddOverallSys("mu_SSBAR", 0.7, 1.3);
	SSBAR_temp.ActivateStatError();
	SSBAR_temp.SetNormalizeByTheory(kFALSE);
	channel.AddSample(SSBAR_temp);
	/* ================================ SSBAR ================================ */

	/* ================================ CHARM ================================ */
	RooStats::HistFactory::Sample CHARM_temp("CHARM_nominal", "CHARM_nominal", fname);
	CHARM_temp.AddHistoSys("track_eff_uncer", "CHARM_track_m", fname, "", "CHARM_track_p", fname, "");
	CHARM_temp.AddHistoSys("KS0_reco_uncer", "CHARM_KS0_m", fname, "", "CHARM_KS0_p", fname, "");
	CHARM_temp.AddHistoSys("Eecl_uncer", "CHARM_Eecl_m", fname, "", "CHARM_Eecl_p", fname, "");
	for (int i = 0; i < NEntryKID; i++) CHARM_temp.AddHistoSys(("KID" + std::to_string(i) + "_uncer").c_str(), ("CHARM_KID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHARM_KID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryPID; i++) CHARM_temp.AddHistoSys(("PID" + std::to_string(i) + "_uncer").c_str(), ("CHARM_PID_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHARM_PID_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntrypi0; i++) CHARM_temp.AddHistoSys(("pi0" + std::to_string(i) + "_uncer").c_str(), ("CHARM_pi0_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHARM_pi0_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
	CHARM_temp.AddShapeSys("CHARM_all_uncorr_uncer", RooStats::HistFactory::Constraint::Type::Gaussian, "CHARM_all_uncorrelated", fname, "");
	CHARM_temp.AddHistoSys("BDTc_shape_qq", "CHARM_BDTc_m", fname, "", "CHARM_BDTc_p", fname, "");
	CHARM_temp.AddNormFactor("qq_CAL_CHARM", 1.1052, 1.1052, 1.1052, true);
	CHARM_temp.AddOverallSys("qq_CAL_CHARM_uncer", 0.7932, 1.2068);
    //CHARM_temp.AddNormFactor("Norm_CHARM_CAL", 1.3, 1.3, 1.3, true);
	CHARM_temp.AddOverallSys("mu_CHARM", 0.7, 1.3);
	CHARM_temp.ActivateStatError();
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
