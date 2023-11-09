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

# define N_decay 38 // five decay mode + others

# define MyEPSILON 0.000001
# define MCTYPE "MC15ri"

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

# define FEI_cal_Bc_num 12
# define FEI_cal_B0_num 11
double FEI_cal_Bc[FEI_cal_Bc_num] = { 1.04, 0.79, 0.69, 0.56, 0.97, 0.95, 0.74, 0.57, 0.91, 0.51, 0.34, 0.59 };
double FEI_cal_Bc_uncertainty[FEI_cal_Bc_num] = { 0.03, 0.03, 0.05, 0.11, 0.03, 0.03, 0.02, 0.06, 0.1, 0.13, 0.07, 0.02 }; // not relative uncertainty. absolute uncertainty
double FEI_cal_Bc_modeID[FEI_cal_Bc_num] = { 0.0, 1.0, 3.0, 4.0, 15.0, 16.0, 18.0, 19.0, 23.0, 24.0, 30.0, -1.0 };
double FEI_cal_B0[FEI_cal_B0_num] = { 1.16, 0.94, 0.81, 0.79, 0.99, 1.03, 0.67, 0.66, 0.69, 0.49, 0.79 };
double FEI_cal_B0_uncertainty[FEI_cal_B0_num] = { 0.04, 0.05, 0.06, 0.02, 0.03, 0.06, 0.02, 0.03, 0.02, 0.02, 0.12 }; // not relative uncertainty. absolute uncertainty
double FEI_cal_B0_modeID[FEI_cal_B0_num] = { 0.0, 1.0, 3.0, 4.0, 5.0, 15.0, 16.0, 18.0, 19.0, 26.0, -1.0 };

# define RarityBins 6

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

	int N_hist = -1;

	FILE* fp;
	fp = fopen(dirname, "r");
	fprintf(fp, "%d", N_hist);
	fclose(fp);

	return N_hist;
}

class Corrector_Fragmentation {
private:
	static constexpr int N_Category_gamma = 9;
	static constexpr int N_Bin_gamma = 5;
	static constexpr double Bins_gamma[N_Bin_gamma - 1] = { 1.15, 1.5, 2.0, 2.4 };
	static constexpr double Fragmentation_Xsgamma[N_Bin_gamma][N_Category_gamma] = { // unit: percentage with respect to "un-missing" mode
		{0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00},
		{10.6, 5.88, 23.2, 44.5, 0.46, 9.99, 0.52, 4.78, 0.00},
		{3.12, 1.13, 15.7, 20.6, 9.48, 26.9, 5.29, 15.5, 2.20},
		{1.65, 0.51, 9.02, 7.98, 16.8, 26.7, 10.3, 21.8, 5.38},
		{0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}
	};
	static constexpr double Fragmentation_Uncertainty_Xsgamma[N_Bin_gamma][N_Category_gamma] = { // relative uncertainty
		{0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000},
		{0.0651, 0.0990, 0.0554, 0.0400, 0.6304, 0.1175, 0.7115, 0.2395, 0.0000},
		{0.1053, 0.2727, 0.0686, 0.0741, 0.0942, 0.0733, 0.2777, 0.1715, 0.1168},
		{0.4250, 1.4490, 0.2646, 0.4003, 0.2509, 0.2671, 0.8409, 0.6682, 0.2898},
		{0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000}
	};

	static constexpr double Nevt_Nominal_before_Xsgamma_MC15[N_Bin_gamma][N_Category_gamma] = {
		{228.010067, 113.889034, 97.106455, 99.597443, 2.350444, 17.694690, 0.265738, 25.726503, 0.000112}, // [-inf, 1.15]
		{1217.857056, 608.996474, 683.434991, 821.802920, 54.319919, 320.575746, 29.362053, 194.574254, 0.829998}, // [1.15, 1.5]
		{745.156145, 372.620424, 554.609584, 717.744130, 132.646240, 507.837684, 190.020473, 312.058874, 120.779580}, // [1.5, 2.0]
		{125.738141, 62.841816, 119.699133, 149.610098, 58.094970, 157.286832, 123.180903, 139.932048, 49.973636}, // [2.0, 2.4]
		{22.820347, 11.451874, 24.859499, 30.438579, 15.176837, 37.050655, 37.906534, 38.685844, 11.884932} // [2.4, inf]
	};
	static constexpr double Nevt_Hmb_before_Xsgamma_MC15[N_Bin_gamma][N_Category_gamma] = {
		{234.717370, 116.931851, 100.045961, 102.811438, 2.443726, 18.213362, 0.283348, 26.352459, 0.000000}, // [-inf, 1.15]
		{1237.919103, 618.966253, 694.727432, 833.619261, 54.873905, 324.784821, 29.598600, 196.698641, 0.862629}, // [1.15, 1.5]
		{739.765455, 368.897914, 550.204743, 712.129243, 131.851781, 504.040211, 187.911583, 309.873035, 119.610335}, // [1.5, 2.0]
		{122.975491, 61.761023, 116.640122, 146.325654, 56.851972, 153.571830, 120.538083, 137.029154, 49.297575}, // [2.0, 2.4]
		{22.077625, 11.336836, 24.298764, 29.398639, 14.646905, 35.886317, 36.466259, 37.697931, 11.619447} // [2.4, inf]
	};
	static constexpr double Nevt_Hpf_before_Xsgamma_MC15[N_Bin_gamma][N_Category_gamma] = {
		{211.414635, 105.941646, 89.984692, 92.845350, 2.132717, 16.210070, 0.250627, 23.805551, 0.000569}, // [-inf, 1.15]
		{1145.424694, 573.572296, 643.905532, 775.309971, 51.161825, 303.573493, 27.884826, 183.335197, 0.805736}, // [1.15, 1.5]
		{736.485305, 369.574161, 550.673698, 711.747889, 132.674662, 505.825855, 190.843751, 312.993403, 121.016612}, // [1.5, 2.0]
		{137.460490, 68.603354, 131.880622, 163.688524, 63.757702, 172.478521, 135.887581, 153.778110, 54.689402}, // [2.0, 2.4]
		{29.109281, 14.773702, 32.144158, 39.024643, 19.707985, 47.791272, 49.197184, 50.159756, 15.393167} // [2.4, inf]
	};
	static constexpr double Nevt_Htransition_before_Xsgamma_MC15[N_Bin_gamma][N_Category_gamma] = {
		{0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, // [-inf, 1.15]
		{1126.011312, 562.651975, 658.733369, 809.726662, 55.294719, 337.131227, 31.789385, 187.359070, 0.920546}, // [1.15, 1.5]
		{838.951829, 419.735045, 624.498178, 809.052632, 149.312766, 571.803869, 214.105973, 351.790571, 136.066782}, // [1.5, 2.0]
		{141.771711, 71.022628, 135.144794, 168.801789, 65.504295, 176.864279, 139.131693, 157.662450, 56.155419}, // [2.0, 2.4]
		{25.757474, 12.884694, 27.851586, 34.061168, 16.946366, 41.803329, 42.925015, 43.435334, 13.297957} // [2.4, inf]
	};
	static constexpr double Nevt_Lmb_before_Xsgamma_MC15[N_Bin_gamma][N_Category_gamma] = {
		{219.785163, 110.184748, 93.899905, 96.479643, 2.274163, 17.081852, 0.246506, 24.692561, 0.000000}, // [-inf, 1.15]
		{1195.125885, 599.131943, 672.175942, 810.139571, 53.494062, 316.373461, 29.121914, 191.567680, 0.804038}, // [1.15, 1.5]
		{750.803012, 374.158691, 558.976621, 723.091791, 134.075592, 512.438695, 192.030343, 314.970543, 122.127927}, // [1.5, 2.0]
		{128.850611, 64.181819, 121.971830, 153.567926, 59.265244, 160.822175, 126.599064, 143.704766, 51.333129}, // [2.0, 2.4]
		{23.591953, 11.740537, 25.576218, 31.342805, 15.603211, 38.154400, 39.437907, 39.746047, 12.199234} // [2.4, inf]
	};
	static constexpr double Nevt_Lpf_before_Xsgamma_MC15[N_Bin_gamma][N_Category_gamma] = {
		{250.065120, 125.068681, 106.236894, 109.408332, 2.586983, 19.616770, 0.297757, 28.249527, 0.000000}, // [-inf, 1.15]
		{1307.816445, 653.786860, 732.097197, 879.672426, 58.091454, 342.326888, 31.063100, 208.314601, 0.883302}, // [1.15, 1.5]
		{747.709095, 373.918779, 554.600603, 718.457434, 131.129261, 504.948261, 186.036870, 308.211200, 119.650476}, // [1.5, 2.0]
		{109.711774, 55.152876, 104.872451, 130.564821, 50.850965, 137.527349, 107.684824, 121.713044, 43.870882}, // [2.0, 2.4]
		{16.301694, 8.096028, 17.764354, 21.644725, 10.719722, 26.442162, 26.669031, 27.403075, 8.479496} // [2.4, inf]
	};
	static constexpr double Nevt_Ltransition_before_Xsgamma_MC15[N_Bin_gamma][N_Category_gamma] = {
		{687.014043, 343.100287, 227.441995, 234.500737, 5.644878, 35.769345, 0.466491, 61.879282, 0.000000}, // [-inf, 1.15]
		{1095.635246, 549.051012, 615.383772, 741.741024, 48.963482, 288.728398, 26.536985, 175.580393, 0.736949}, // [1.15, 1.5]
		{671.592637, 336.051687, 499.387200, 646.831974, 119.008094, 456.176166, 171.088734, 281.417112, 108.985261}, // [1.5, 2.0]
		{113.378578, 56.429701, 107.697531, 134.644214, 52.336623, 141.679228, 110.495689, 126.585786, 45.203748}, // [2.0, 2.4]
		{20.750444, 10.246957, 22.414418, 27.261420, 13.458927, 33.490579, 34.139560, 34.988253, 10.795822} // [2.4, inf]
	};

	static constexpr double Nevt_Nominal_missing_before_Xsgamma_MC15[N_Bin_gamma] = { 48.015884, 360.363091, 676.491505, 518.892840, 292.939844 }; // number of missing mode (MC)
	static constexpr double Nevt_Hmb_missing_before_Xsgamma_MC15[N_Bin_gamma] = { 49.450546, 365.252264, 669.873825, 507.353413, 284.711591 }; // number of missing mode (MC)
	static constexpr double Nevt_Hpf_missing_before_Xsgamma_MC15[N_Bin_gamma] = { 44.686813, 340.083522, 676.515003, 574.044827, 388.939028 }; // number of missing mode (MC)
	static constexpr double Nevt_Htransition_missing_before_Xsgamma_MC15[N_Bin_gamma] = { 0.000000, 348.935210, 762.920563, 585.338309, 330.039376 }; // number of missing mode (MC)
	static constexpr double Nevt_Lmb_missing_before_Xsgamma_MC15[N_Bin_gamma] = { 46.502106, 354.423602, 683.769325, 533.204382, 302.346798 }; // number of missing mode (MC)
	static constexpr double Nevt_Lpf_missing_before_Xsgamma_MC15[N_Bin_gamma] = { 52.743632, 386.345000, 669.316781, 448.947976, 200.139315 }; // number of missing mode (MC)
	static constexpr double Nevt_Ltransition_missing_before_Xsgamma_MC15[N_Bin_gamma] = { 93.606503, 323.683107, 609.133190, 468.259061, 263.808865 }; // number of missing mode (MC)

	double Total_Nevt_Nominal_before_Xsgamma_MC15[N_Bin_gamma]; // without missing mode
	double Total_Nevt_Hmb_before_Xsgamma_MC15[N_Bin_gamma];
	double Total_Nevt_Hpf_before_Xsgamma_MC15[N_Bin_gamma];
	double Total_Nevt_Htransition_before_Xsgamma_MC15[N_Bin_gamma];
	double Total_Nevt_Lmb_before_Xsgamma_MC15[N_Bin_gamma];
	double Total_Nevt_Lpf_before_Xsgamma_MC15[N_Bin_gamma];
	double Total_Nevt_Ltransition_before_Xsgamma_MC15[N_Bin_gamma];
public:
	enum class SystType { // what kind of Xs sample is corrected?
		Nominal,
		Hmb,
		Hpf,
		Htransition,
		Lmb,
		Lpf,
		Ltransition
	};
	enum class Sample {
		gamma // B -> Xs gamma is used
	};

	Corrector_Fragmentation();
	int Classify(int Decay[N_decay], Sample sample);
	int GetMxBin(double MXs, Sample sample);
	double GetCorrectionFactor(int Decay[N_decay], double MXs, SystType systtype, Sample sample, std::string type);
	double FluctuateCorrection(int Decay[N_decay], double MXs, SystType systtype, int TargetMxsBin, int TargetCategory, bool IsTargetCategoryUp, Sample sample, std::string type);
	int GetNMxsBin(Sample sample);
	int GetNCategory(Sample sample);
};

Corrector_Fragmentation corrector_Fragmentation;

Corrector_Fragmentation::Corrector_Fragmentation() {
	for (int i = 0; i < N_Bin_gamma; i++) {
		Total_Nevt_Nominal_before_Xsgamma_MC15[i] = 0.0;
		Total_Nevt_Hmb_before_Xsgamma_MC15[i] = 0.0;
		Total_Nevt_Hpf_before_Xsgamma_MC15[i] = 0.0;
		Total_Nevt_Htransition_before_Xsgamma_MC15[i] = 0.0;
		Total_Nevt_Lmb_before_Xsgamma_MC15[i] = 0.0;
		Total_Nevt_Lpf_before_Xsgamma_MC15[i] = 0.0;
		Total_Nevt_Ltransition_before_Xsgamma_MC15[i] = 0.0;
	}

	for (int i = 0; i < N_Bin_gamma; i++) {
		for (int j = 0; j < N_Category_gamma; j++) {
			Total_Nevt_Nominal_before_Xsgamma_MC15[i] = Total_Nevt_Nominal_before_Xsgamma_MC15[i] + Nevt_Nominal_before_Xsgamma_MC15[i][j];
			Total_Nevt_Hmb_before_Xsgamma_MC15[i] = Total_Nevt_Hmb_before_Xsgamma_MC15[i] + Nevt_Hmb_before_Xsgamma_MC15[i][j];
			Total_Nevt_Hpf_before_Xsgamma_MC15[i] = Total_Nevt_Hpf_before_Xsgamma_MC15[i] + Nevt_Hpf_before_Xsgamma_MC15[i][j];
			Total_Nevt_Htransition_before_Xsgamma_MC15[i] = Total_Nevt_Htransition_before_Xsgamma_MC15[i] + Nevt_Htransition_before_Xsgamma_MC15[i][j];
			Total_Nevt_Lmb_before_Xsgamma_MC15[i] = Total_Nevt_Lmb_before_Xsgamma_MC15[i] + Nevt_Lmb_before_Xsgamma_MC15[i][j];
			Total_Nevt_Lpf_before_Xsgamma_MC15[i] = Total_Nevt_Lpf_before_Xsgamma_MC15[i] + Nevt_Lpf_before_Xsgamma_MC15[i][j];
			Total_Nevt_Ltransition_before_Xsgamma_MC15[i] = Total_Nevt_Ltransition_before_Xsgamma_MC15[i] + Nevt_Ltransition_before_Xsgamma_MC15[i][j];
		}
	}
}

int Corrector_Fragmentation::Classify(int Decay[N_decay], Sample sample) {
	if (sample == Corrector_Fragmentation::Sample::gamma) {
		/*
			0: Xs2Kpi_wopi0_MC // Xs -> K pi w/o pi0
			1: Xs2Kpi_wpi0_MC, // Xs -> K pi w/ 1pi0
			2: Xs2Kpipi_wopi0_MC, // Xs -> K pi pi w/o pi0
			3: Xs2Kpipi_wpi0_MC, // Xs -> K pi pi w/ 1pi0
			4: Xs2Kpipipi_wopi0_MC, // Xs -> K pi pi pi w/o pi0
			5: Xs2Kpipipi_wpi0_MC, // Xs -> K pi pi pi w/ 1pi0
			6: Xs2Kpipipipi_atmost1pi0_MC, // Xs -> K pi pi pi pi at most 1pi0
			7: Xs2pi0pi0_MC, // including 2pi0
			8: Xs2KKK_MC, // including 3K
			9: MC_other, // other non-resonant Xs decay
		*/
		if (Decay[6] > 0 || Decay[24] > 0) return 0;
		else if (Decay[5] > 0 || Decay[25] > 0) return 1;
		else if (Decay[7] > 0 || Decay[27] > 0) return 2;
		else if (Decay[8] > 0 || Decay[26] > 0) return 3;
		else if (Decay[10] > 0 || Decay[28] > 0) return 4;
		else if (Decay[9] > 0 || Decay[29] > 0) return 5;
		else if (Decay[11] > 0 || Decay[12] > 0 || Decay[30] > 0 || Decay[31] > 0) return 6;
		else if (Decay[13] > 0 || Decay[14] > 0 || Decay[15] > 0 || Decay[32] > 0 || Decay[33] > 0 || Decay[34] > 0) return 7;
		else if (Decay[16] > 0 || Decay[17] > 0 || Decay[18] > 0 || Decay[35] > 0 || Decay[36] > 0 || Decay[37] > 0) return 8;
		else return N_Category_gamma;
	}
	else {
		printf("[Corrector_Fragmentation::Classify] cannot find sample");
		exit(1);
		return 1;
	}
}

int Corrector_Fragmentation::GetMxBin(double MXs, Sample sample) {
	if (sample == Corrector_Fragmentation::Sample::gamma) { // [-inf, 1.15, 1.5, 2.0, 2.4, inf]
		for (int i = 0; i < N_Bin_gamma - 1; i++) {
			if (MXs < Bins_gamma[i]) return i;
		}
		return N_Bin_gamma - 1;
	}
	else {
		printf("[Corrector_Fragmentation::GetMxBin] cannot find sample");
		exit(1);
		return 1;
	}
}

double Corrector_Fragmentation::GetCorrectionFactor(int Decay[N_decay], double MXs, SystType systtype, Sample sample, std::string type) {
	int Category = Classify(Decay, sample);
	int MxsBin = GetMxBin(MXs, sample);

	if (sample == Corrector_Fragmentation::Sample::gamma) {

		if (Category == N_Category_gamma) return 1.0; // no correction if it is missing mode

		if ((MxsBin == 0) || (MxsBin == N_Bin_gamma - 1)) return 1.0; // no correction if Mxs bin is out of range

		if (type == "MC15ri") {
			if (systtype == Corrector_Fragmentation::SystType::Nominal) {
				if (Nevt_Nominal_before_Xsgamma_MC15[MxsBin][Category] < MyEPSILON) return 1.0; // no correction if Nevt is 0 at MC
				return (0.01 * Fragmentation_Xsgamma[MxsBin][Category]) / (Nevt_Nominal_before_Xsgamma_MC15[MxsBin][Category] / Total_Nevt_Nominal_before_Xsgamma_MC15[MxsBin]);
			}
			else if (systtype == Corrector_Fragmentation::SystType::Hmb) {
				if (Nevt_Hmb_before_Xsgamma_MC15[MxsBin][Category] < MyEPSILON) return 1.0; // no correction if Nevt is 0 at MC
				return (0.01 * Fragmentation_Xsgamma[MxsBin][Category]) / (Nevt_Hmb_before_Xsgamma_MC15[MxsBin][Category] / Total_Nevt_Hmb_before_Xsgamma_MC15[MxsBin]);
			}
			else if (systtype == Corrector_Fragmentation::SystType::Hpf) {
				if (Nevt_Hpf_before_Xsgamma_MC15[MxsBin][Category] < MyEPSILON) return 1.0; // no correction if Nevt is 0 at MC
				return (0.01 * Fragmentation_Xsgamma[MxsBin][Category]) / (Nevt_Hpf_before_Xsgamma_MC15[MxsBin][Category] / Total_Nevt_Hpf_before_Xsgamma_MC15[MxsBin]);
			}
			else if (systtype == Corrector_Fragmentation::SystType::Htransition) {
				if (Nevt_Htransition_before_Xsgamma_MC15[MxsBin][Category] < MyEPSILON) return 1.0; // no correction if Nevt is 0 at MC
				return (0.01 * Fragmentation_Xsgamma[MxsBin][Category]) / (Nevt_Htransition_before_Xsgamma_MC15[MxsBin][Category] / Total_Nevt_Htransition_before_Xsgamma_MC15[MxsBin]);
			}
			else if (systtype == Corrector_Fragmentation::SystType::Lmb) {
				if (Nevt_Lmb_before_Xsgamma_MC15[MxsBin][Category] < MyEPSILON) return 1.0; // no correction if Nevt is 0 at MC
				return (0.01 * Fragmentation_Xsgamma[MxsBin][Category]) / (Nevt_Lmb_before_Xsgamma_MC15[MxsBin][Category] / Total_Nevt_Lmb_before_Xsgamma_MC15[MxsBin]);
			}
			else if (systtype == Corrector_Fragmentation::SystType::Lpf) {
				if (Nevt_Lpf_before_Xsgamma_MC15[MxsBin][Category] < MyEPSILON) return 1.0; // no correction if Nevt is 0 at MC
				return (0.01 * Fragmentation_Xsgamma[MxsBin][Category]) / (Nevt_Lpf_before_Xsgamma_MC15[MxsBin][Category] / Total_Nevt_Lpf_before_Xsgamma_MC15[MxsBin]);
			}
			else if (systtype == Corrector_Fragmentation::SystType::Ltransition) {
				if (Nevt_Ltransition_before_Xsgamma_MC15[MxsBin][Category] < MyEPSILON) return 1.0; // no correction if Nevt is 0 at MC
				return (0.01 * Fragmentation_Xsgamma[MxsBin][Category]) / (Nevt_Ltransition_before_Xsgamma_MC15[MxsBin][Category] / Total_Nevt_Ltransition_before_Xsgamma_MC15[MxsBin]);
			}
		}
		else {
			printf("[Corrector_Fragmentation::GetCorrectionFactor] Invalid type!\n");
			exit(1);
			return 0;
		}
	}
	else {
		printf("[Corrector_Fragmentation::GetCorrectionFactor] cannot find sample");
		exit(1);
		return 1;
	}
}

double Corrector_Fragmentation::FluctuateCorrection(int Decay[N_decay], double MXs, SystType systtype, int TargetMxsBin, int TargetCategory, bool IsTargetCategoryUp, Sample sample, std::string type) {
	/*
	Nevt of TargetCategory in TargetMxsBin is fluctuated.
	we calculate the change of Decay[N_decay] (Category)
	*/
	int Category = Classify(Decay, sample);
	int MxsBin = GetMxBin(MXs, sample);

	if (TargetMxsBin != MxsBin) return 1.0; // no fluctuation if it is not in target MXs region

	const double RelativeUncertainty = 1.0;

	if (sample == Corrector_Fragmentation::Sample::gamma) {

		double TotalNevtAtMxsBinWithMissing = Total_Nevt_Nominal_before_Xsgamma_MC15[TargetMxsBin] + Nevt_Nominal_missing_before_Xsgamma_MC15[TargetMxsBin];
		double TargetNevtAtMxsBin;
		if (TargetCategory == N_Category_gamma) TargetNevtAtMxsBin = Nevt_Nominal_missing_before_Xsgamma_MC15[TargetMxsBin]; // if it is missing mode
		else TargetNevtAtMxsBin = Total_Nevt_Nominal_before_Xsgamma_MC15[TargetMxsBin] * Fragmentation_Xsgamma[TargetMxsBin][TargetCategory];
		double TotalNevtAtMxsBinWithMissingWithoutTargetCategory = TotalNevtAtMxsBinWithMissing - TargetNevtAtMxsBin;

		if (Category == TargetCategory) {
			if (IsTargetCategoryUp) return (1.0 + RelativeUncertainty);
			else return (1.0 - RelativeUncertainty);
		}
		else {
			if (IsTargetCategoryUp) {
				double output = (TotalNevtAtMxsBinWithMissingWithoutTargetCategory - TargetNevtAtMxsBin * RelativeUncertainty) / TotalNevtAtMxsBinWithMissingWithoutTargetCategory;
				if (output < 0) {
					printf("[Corrector_Fragmentation::FluctuateCorrection] minus Nevt!\n");
					exit(1);
				}
				return output;
			}
			else {
				double output = (TotalNevtAtMxsBinWithMissingWithoutTargetCategory + TargetNevtAtMxsBin * RelativeUncertainty) / TotalNevtAtMxsBinWithMissingWithoutTargetCategory;
				if (output < 0) {
					printf("[Corrector_Fragmentation::FluctuateCorrection] minus Nevt!\n");
					exit(1);
				}
				return output;
			}
		}
	}
	else {
		printf("[Corrector_Fragmentation::GetCorrectionFactor] cannot find sample");
		exit(1);
		return 1;
	}
}

int Corrector_Fragmentation::GetNMxsBin(Sample sample) {
	if (sample == Corrector_Fragmentation::Sample::gamma) return N_Bin_gamma;
	else {
		printf("[Corrector_Fragmentation::GetMxBin] cannot find sample");
		exit(1);
		return 1;
	}
}

int Corrector_Fragmentation::GetNCategory(Sample sample) {
	if (sample == Corrector_Fragmentation::Sample::gamma) return (N_Category_gamma + 1); // visible + missing
	else {
		printf("[Corrector_Fragmentation::GetMxBin] cannot find sample");
		exit(1);
		return 1;
	}
}

int WorkSpace() {

	int NEntryFEI = ReadNFEIEigenVector("./FEI_selected.txt");
	int NEntryKID = ReadNPIDEigenVector("./KID_selected.txt");
	int NEntryPID = ReadNPIDEigenVector("./PID_selected.txt");
	int NEntryBR = ReadNBREigenVector("./BR_selected.txt");
	int NEntrypi0 = ReadNpi0EigenVector("./pi0_selected.txt");
	int NEntryMultiplicity = ReadMultiplicityInfo("./Multiplicity_info.txt");

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
	sig_temp.AddHistoSys("track_eff_uncer", "Signal_track_m", fname, "", "Signal_track_p", fname, "");
	sig_temp.AddHistoSys("KS0_reco_uncer", "Signal_KS0_m", fname, "", "Signal_KS0_p", fname, "");
	sig_temp.AddHistoSys("BtoDtoXKL_uncer", "Signal_BtoDtoXKL_m", fname, "", "Signal_BtoDtoXKL_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) sig_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("Signal_multiplicity" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_multiplicity" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryFEI; i++) sig_temp.AddHistoSys(("FEI" + std::to_string(i) + "_uncer").c_str(), ("Signal_FEI_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("Signal_FEI_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
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
	for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
		for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
			int temp_index = MxsBin * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma) + Category;

			std::string uncertainty_name = "Xs_fragmentation" + std::to_string(temp_index) + "_uncer";
			std::string plus_name = std::string("Signal_Xs_frag_decay_") + std::to_string(MxsBin) + std::string("_") + std::to_string(Category) + std::string("_p");
			std::string minus_name = std::string("Signal_Xs_frag_decay_") + std::to_string(MxsBin) + std::string("_") + std::to_string(Category) + std::string("_m");

			sig_temp.AddHistoSys(uncertainty_name.c_str(), plus_name.c_str(), fname, "", minus_name.c_str(), fname, "");
		}
	}
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
	CHG_temp.AddHistoSys("track_eff_uncer", "CHG_track_m", fname, "", "CHG_track_p", fname, "");
	CHG_temp.AddHistoSys("KS0_reco_uncer", "CHG_KS0_m", fname, "", "CHG_KS0_p", fname, "");
	CHG_temp.AddHistoSys("BtoDtoXKL_uncer", "CHG_BtoDtoXKL_m", fname, "", "CHG_BtoDtoXKL_p", fname, "");
	CHG_temp.AddHistoSys("Knn_BR_uncer", "CHG_Knn_m", fname, "", "CHG_Knn_p", fname, "");
	CHG_temp.AddHistoSys("Kstarnn_BR_uncer", "CHG_Kstarnn_m", fname, "", "CHG_Kstarnn_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) CHG_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("CHG_multiplicity" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_multiplicity" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryFEI; i++) CHG_temp.AddHistoSys(("FEI" + std::to_string(i) + "_uncer").c_str(), ("CHG_FEI_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("CHG_FEI_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
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
	MIX_temp.AddHistoSys("track_eff_uncer", "MIX_track_m", fname, "", "MIX_track_p", fname, "");
	MIX_temp.AddHistoSys("KS0_reco_uncer", "MIX_KS0_m", fname, "", "MIX_KS0_p", fname, "");
	MIX_temp.AddHistoSys("BtoDtoXKL_uncer", "MIX_BtoDtoXKL_m", fname, "", "MIX_BtoDtoXKL_p", fname, "");
	MIX_temp.AddHistoSys("K0nn_BR_uncer", "MIX_K0nn_m", fname, "", "MIX_K0nn_p", fname, "");
	MIX_temp.AddHistoSys("K0starnn_BR_uncer", "MIX_K0starnn_m", fname, "", "MIX_K0starnn_p", fname, "");
	for (int i = 0; i < NEntryMultiplicity; i++) MIX_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("MIX_multiplicity" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_multiplicity" + std::to_string(i) + "_p").c_str(), fname, "");
	for (int i = 0; i < NEntryFEI; i++) MIX_temp.AddHistoSys(("FEI" + std::to_string(i) + "_uncer").c_str(), ("MIX_FEI_correlated" + std::to_string(i) + "_m").c_str(), fname, "", ("MIX_FEI_correlated" + std::to_string(i) + "_p").c_str(), fname, "");
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
	for (int i = 0; i < NEntryMultiplicity; i++) UUBAR_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("UUBAR_multiplicity" + std::to_string(i) + "_m").c_str(), fname, "", ("UUBAR_multiplicity" + std::to_string(i) + "_p").c_str(), fname, "");
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
	for (int i = 0; i < NEntryMultiplicity; i++) DDBAR_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("DDBAR_multiplicity" + std::to_string(i) + "_m").c_str(), fname, "", ("DDBAR_multiplicity" + std::to_string(i) + "_p").c_str(), fname, "");
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
	for (int i = 0; i < NEntryMultiplicity; i++) SSBAR_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("SSBAR_multiplicity" + std::to_string(i) + "_m").c_str(), fname, "", ("SSBAR_multiplicity" + std::to_string(i) + "_p").c_str(), fname, "");
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
	for (int i = 0; i < NEntryMultiplicity; i++) CHARM_temp.AddHistoSys(("multiplicity" + std::to_string(i) + "_uncer").c_str(), ("CHARM_multiplicity" + std::to_string(i) + "_m").c_str(), fname, "", ("CHARM_multiplicity" + std::to_string(i) + "_p").c_str(), fname, "");
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
