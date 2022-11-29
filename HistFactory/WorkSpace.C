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

// arXiv:1409.4557v2
# define TB0 1.5195 // (Table. 1)
# define TBp 1.6384 // (Table. 1)
# define BR_Kplus_nunubar 0.00000398 // (eq. 10)
# define BR_K0star_nunubar 0.00000919 // (eq. 11)
# define BR_K0_nunubar (BR_Kplus_nunubar*TB0/TBp) // under (eq. 15)
# define BR_Kplusstar_nunubar (BR_K0star_nunubar*TBp/TB0) // under (eq. 15)
# define BR_Xs_nunubar 0.000029 // (eq. 23)
# define BR_Xsu_nonresonant_nunubar (BR_Xs_nunubar - BR_Kplus_nunubar - BR_Kplusstar_nunubar)
# define BR_Xsd_nonresonant_nunubar (BR_Xs_nunubar - BR_K0_nunubar - BR_K0star_nunubar)

// uncertainty of BR from parametric reason
# define Sigma_BR_Kplus_nunubar 0.00000019
# define Sigma_BR_K0star_nunubar 0.00000050
# define Sigma_BR_K0_nunubar (Sigma_BR_Kplus_nunubar*TB0/TBp)
# define Sigma_BR_Kplusstar_nunubar (Sigma_BR_K0star_nunubar*TBp/TB0)
# define Sigma_BR_Xs_nunubar 0.000003

// https://confluence.desy.de/pages/viewpage.action?pageId=107054222
# define N_BpBp_1invab 565400000.0
# define N_B0B0_1invab 534600000.0

# define N_Kplus_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Kplus_nunubar)
# define N_Kplusstar_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Kplusstar_nunubar)
# define N_Xsu_nonresonant_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Xsu_nonresonant_nunubar)
# define N_K0_nunubar_1invab (2.0 * N_B0B0_1invab * BR_K0_nunubar)
# define N_K0star_nunubar_1invab (2.0 * N_B0B0_1invab * BR_K0star_nunubar)
# define N_Xsd_nunubar_1invab (2.0 * N_B0B0_1invab * BR_Xsd_nonresonant_nunubar)

# define Lumi_validation_MC 0.3 // ab-1

// my MC sample number
# define N_Kplus_nunubar 10000000.0
# define N_K0_nunubar 10000000.0
# define N_Kplusstar_nunubar 10000000.0
# define N_K0star_nunubar 10000000.0
# define N_Xsu_nonresonant_nunubar 50000000.0
# define N_Xsd_nonresonant_nunubar 50000000.0

// scale factor for each MC sample
# define Scale_Kplus (N_Kplus_nunubar_1invab/N_Kplus_nunubar)
# define Scale_Kplusstar (N_Kplusstar_nunubar_1invab/N_Kplusstar_nunubar)
# define Scale_Xsu_nonresonant (N_Xsu_nonresonant_nunubar_1invab/N_Xsu_nonresonant_nunubar)
# define Scale_K0 (N_K0_nunubar_1invab/N_K0_nunubar)
# define Scale_K0star (N_K0star_nunubar_1invab/N_K0star_nunubar)
# define Scale_Xsd_nonresonant (N_Xsd_nunubar_1invab/N_Xsd_nonresonant_nunubar)

# define KS0_rel_uncertainty 0.6 // %/cm
# define track_rel_uncertainty 0.69 // %
# define pi0_correction 0.932
# define pi0_rel_uncertainty ((0.0369 / 0.932) * 100.0) // %
# define Kaon_PID_max_uncertainty 0.1 // not percentage. relative uncertainty
// https://indico.belle2.org/event/6872/contributions/37447/attachments/17127/25504/FEIperformance_B2GM.pdf
# define FEI_cal_Bc 0.679
# define FEI_cal_Bc_uncertainty (0.017/FEI_cal_Bc) // not percentage. relative uncertainty
# define FEI_cal_B0 0.713
# define FEI_cal_B0_uncertainty (0.019/FEI_cal_B0) // not percentage. relative uncertainty

//# define Htransition_Xsu_change (-0.142)
//# define Htransition_Xsd_change (-0.099)
//# define Ltransition_Xsu_change (-0.002)
//# define Ltransition_Xsd_change (0.067)
//# define Hmb_Xsu_change (-0.067)
//# define Hmb_Xsd_change (-0.058)
//# define Lmb_Xsu_change (-0.062)
//# define Lmb_Xsd_change (-0.046)
//# define Hpf_Xsu_change (-0.120)
//# define Hpf_Xsd_change (-0.073)
//# define Lpf_Xsu_change (-0.014)
//# define Lpf_Xsd_change (0.019)

//# define FEI_cal_Bc 1.0
//# define FEI_cal_Bc_uncertainty 0.0
//# define FEI_cal_B0 1.0
//# define FEI_cal_B0_uncertainty 0.0

//# define Htransition_Xsu_change 0.0
//# define Htransition_Xsd_change 0.0
//# define Ltransition_Xsu_change 0.0
//# define Ltransition_Xsd_change 0.0
//# define Hmb_Xsu_change 0.0
//# define Hmb_Xsd_change 0.0
//# define Lmb_Xsu_change 0.0
//# define Lmb_Xsd_change 0.0
//# define Hpf_Xsu_change 0.0
//# define Hpf_Xsd_change 0.0
//# define Lpf_Xsu_change 0.0
//# define Lpf_Xsd_change 0.0

# define RarityBins 14

int WorkSpace() {

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
	sig_temp.AddHistoSys("pi0_reco_uncer", "Signal_pi0_m", fname, "", "Signal_pi0_p", fname, "");
	sig_temp.AddHistoSys("track_eff_uncer", "Signal_track_m", fname, "", "Signal_track_p", fname, "");
	sig_temp.AddHistoSys("KS0_reco_uncer", "Signal_KS0_m", fname, "", "Signal_KS0_p", fname, "");
	sig_temp.AddHistoSys("KID_eff_uncer", "Signal_KID_m", fname, "", "Signal_KID_p", fname, "");
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
	sig_temp.AddOverallSys("FBDT_efficiency_uncer", 0.85, 1.15);
	sig_temp.ActivateStatError();
	sig_temp.SetNormalizeByTheory(kFALSE);
	sig_temp.AddNormFactor("mu", expmu, -100.0, 100.0);
	channel.AddSample(sig_temp);

	// read background template
	/* ================================ CHG ================================ */
	RooStats::HistFactory::Sample CHG_temp("CHG_nominal", "CHG_nominal", fname);
	CHG_temp.AddOverallSys("FEI_charged_uncer", 1 - FEI_cal_Bc_uncertainty, 1 + FEI_cal_Bc_uncertainty);
	CHG_temp.AddHistoSys("pi0_reco_uncer", "CHG_pi0_m", fname, "", "CHG_pi0_p", fname, "");
	CHG_temp.AddHistoSys("track_eff_uncer", "CHG_track_m", fname, "", "CHG_track_p", fname, "");
	CHG_temp.AddHistoSys("KS0_reco_uncer", "CHG_KS0_m", fname, "", "CHG_KS0_p", fname, "");
	CHG_temp.AddHistoSys("KID_eff_uncer", "CHG_KID_m", fname, "", "CHG_KID_p", fname, "");
	CHG_temp.AddHistoSys("BDTc_shape_BB", "CHG_BDTc_m", fname, "", "CHG_BDTc_p", fname, "");
    //CHG_temp.AddNormFactor("Norm_CHG_CAL", 1.3, 1.3, 1.3, true);
	CHG_temp.AddOverallSys("mu_CHG", 0.8, 1.2);
	CHG_temp.ActivateStatError();
	CHG_temp.SetNormalizeByTheory(kFALSE);
	channel.AddSample(CHG_temp);
	/* ================================ CHG ================================ */

	/* ================================ MIX ================================ */
	RooStats::HistFactory::Sample MIX_temp("MIX_nominal", "MIX_nominal", fname);
	MIX_temp.AddOverallSys("FEI_neutral_uncer", 1 - FEI_cal_B0_uncertainty, 1 + FEI_cal_B0_uncertainty);
	MIX_temp.AddHistoSys("pi0_reco_uncer", "MIX_pi0_m", fname, "", "MIX_pi0_p", fname, "");
	MIX_temp.AddHistoSys("track_eff_uncer", "MIX_track_m", fname, "", "MIX_track_p", fname, "");
	MIX_temp.AddHistoSys("KS0_reco_uncer", "MIX_KS0_m", fname, "", "MIX_KS0_p", fname, "");
	MIX_temp.AddHistoSys("KID_eff_uncer", "MIX_KID_m", fname, "", "MIX_KID_p", fname, "");
	MIX_temp.AddHistoSys("BDTc_shape_BB", "MIX_BDTc_m", fname, "", "MIX_BDTc_p", fname, "");
    //MIX_temp.AddNormFactor("Norm_MIX_CAL", 1.3, 1.3, 1.3, true);
	MIX_temp.AddOverallSys("mu_MIX", 0.8, 1.2);
	MIX_temp.ActivateStatError();
	MIX_temp.SetNormalizeByTheory(kFALSE);
	channel.AddSample(MIX_temp);
	/* ================================ MIX ================================ */

	/* ================================ UUBAR ================================ */
	RooStats::HistFactory::Sample UUBAR_temp("UUBAR_nominal", "UUBAR_nominal", fname);
	UUBAR_temp.AddHistoSys("pi0_reco_uncer", "UUBAR_pi0_m", fname, "", "UUBAR_pi0_p", fname, "");
	UUBAR_temp.AddHistoSys("track_eff_uncer", "UUBAR_track_m", fname, "", "UUBAR_track_p", fname, "");
	UUBAR_temp.AddHistoSys("KS0_reco_uncer", "UUBAR_KS0_m", fname, "", "UUBAR_KS0_p", fname, "");
	UUBAR_temp.AddHistoSys("KID_eff_uncer", "UUBAR_KID_m", fname, "", "UUBAR_KID_p", fname, "");
	UUBAR_temp.AddHistoSys("BDTc_shape_qq", "UUBAR_BDTc_m", fname, "", "UUBAR_BDTc_p", fname, "");
	UUBAR_temp.AddNormFactor("qq_CAL_UUBAR", 1.2323, 1.2323, 1.2323, true);
	UUBAR_temp.AddOverallSys("qq_CAL_UUBAR_uncer", 0.9355, 1.0645);
    //UUBAR_temp.AddNormFactor("Norm_UUBAR_CAL", 1.3, 1.3, 1.3, true);
	UUBAR_temp.AddOverallSys("mu_UUBAR", 0.8, 1.2);
	UUBAR_temp.ActivateStatError();
	UUBAR_temp.SetNormalizeByTheory(kFALSE);
	channel.AddSample(UUBAR_temp);
	/* ================================ UUBAR ================================ */

	/* ================================ DDBAR ================================ */
	RooStats::HistFactory::Sample DDBAR_temp("DDBAR_nominal", "DDBAR_nominal", fname);
	DDBAR_temp.AddHistoSys("pi0_reco_uncer", "DDBAR_pi0_m", fname, "", "DDBAR_pi0_p", fname, "");
	DDBAR_temp.AddHistoSys("track_eff_uncer", "DDBAR_track_m", fname, "", "DDBAR_track_p", fname, "");
	DDBAR_temp.AddHistoSys("KS0_reco_uncer", "DDBAR_KS0_m", fname, "", "DDBAR_KS0_p", fname, "");
	DDBAR_temp.AddHistoSys("KID_eff_uncer", "DDBAR_KID_m", fname, "", "DDBAR_KID_p", fname, "");
	DDBAR_temp.AddHistoSys("BDTc_shape_qq", "DDBAR_BDTc_m", fname, "", "DDBAR_BDTc_p", fname, "");
	DDBAR_temp.AddNormFactor("qq_CAL_DDBAR", 1.2323, 1.2323, 1.2323, true);
	DDBAR_temp.AddOverallSys("qq_CAL_DDBAR_uncer", 0.9355, 1.0645);
    //DDBAR_temp.AddNormFactor("Norm_DDBAR_CAL", 1.3, 1.3, 1.3, true);
	DDBAR_temp.AddOverallSys("mu_DDBAR", 0.8, 1.2);
	DDBAR_temp.ActivateStatError();
	DDBAR_temp.SetNormalizeByTheory(kFALSE);
	channel.AddSample(DDBAR_temp);
	/* ================================ DDBAR ================================ */

	/* ================================ SSBAR ================================ */
	RooStats::HistFactory::Sample SSBAR_temp("SSBAR_nominal", "SSBAR_nominal", fname);
	SSBAR_temp.AddHistoSys("pi0_reco_uncer", "SSBAR_pi0_m", fname, "", "SSBAR_pi0_p", fname, "");
	SSBAR_temp.AddHistoSys("track_eff_uncer", "SSBAR_track_m", fname, "", "SSBAR_track_p", fname, "");
	SSBAR_temp.AddHistoSys("KS0_reco_uncer", "SSBAR_KS0_m", fname, "", "SSBAR_KS0_p", fname, "");
	SSBAR_temp.AddHistoSys("KID_eff_uncer", "SSBAR_KID_m", fname, "", "SSBAR_KID_p", fname, "");
	SSBAR_temp.AddHistoSys("BDTc_shape_qq", "SSBAR_BDTc_m", fname, "", "SSBAR_BDTc_p", fname, "");
	SSBAR_temp.AddNormFactor("qq_CAL_SSBAR", 1.2323, 1.2323, 1.2323, true);
	SSBAR_temp.AddOverallSys("qq_CAL_SSBAR_uncer", 0.9355, 1.0645);
    //SSBAR_temp.AddNormFactor("Norm_SSBAR_CAL", 1.3, 1.3, 1.3, true);
	SSBAR_temp.AddOverallSys("mu_SSBAR", 0.8, 1.2);
	SSBAR_temp.ActivateStatError();
	SSBAR_temp.SetNormalizeByTheory(kFALSE);
	channel.AddSample(SSBAR_temp);
	/* ================================ SSBAR ================================ */

	/* ================================ CHARM ================================ */
	RooStats::HistFactory::Sample CHARM_temp("CHARM_nominal", "CHARM_nominal", fname);
	CHARM_temp.AddHistoSys("pi0_reco_uncer", "CHARM_pi0_m", fname, "", "CHARM_pi0_p", fname, "");
	CHARM_temp.AddHistoSys("track_eff_uncer", "CHARM_track_m", fname, "", "CHARM_track_p", fname, "");
	CHARM_temp.AddHistoSys("KS0_reco_uncer", "CHARM_KS0_m", fname, "", "CHARM_KS0_p", fname, "");
	CHARM_temp.AddHistoSys("KID_eff_uncer", "CHARM_KID_m", fname, "", "CHARM_KID_p", fname, "");
	CHARM_temp.AddHistoSys("BDTc_shape_qq", "CHARM_BDTc_m", fname, "", "CHARM_BDTc_p", fname, "");
	CHARM_temp.AddNormFactor("qq_CAL_CHARM", 1.2323, 1.2323, 1.2323, true);
	CHARM_temp.AddOverallSys("qq_CAL_CHARM_uncer", 0.9355, 1.0645);
    //CHARM_temp.AddNormFactor("Norm_CHARM_CAL", 1.3, 1.3, 1.3, true);
	CHARM_temp.AddOverallSys("mu_CHARM", 0.8, 1.2);
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
