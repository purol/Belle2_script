#include <iostream>
#include <stdio.h>
#include <string>
#include <algorithm>
#include <vector>
#include <numeric>
#include <time.h>
#include <random>
#include <sstream>

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
#include "RooMsgService.h"

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

#include "RooStats/FrequentistCalculator.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterPlot.h"

using namespace RooFit;
using namespace RooStats;
using namespace HistFactory;

using std::string;
using std::to_string;
using std::cout;
using std::endl;

# define MyDEBUG true
void Debug(RooWorkspace* w, RooFitResult* fitres, RooDataSet* data);

# define N_Needed_info 37
//# define N_event_info 15
# define N_Upsilon_info 164
# define N_Bsig_info 738
# define N_Btag_info 11
# define N_decay 38 // five decay mode + others
# define N_decay_nparticles 5 // # of nu_e, B->Xs nu_e nu_e_bar, B0->Xs nu_e nu_e_bar, B+-, B0
# define N_decay_syst_ff 7 // helicity angle + q2
# define N_PID_syst 73
# define N_fakeE_syst 37
# define N_fakeMU_syst 49
# define N_pi0_syst 8
# define index_q2 0
# define index_MXs_Bc 3
# define index_MXs_B0 4

# define MyEPSILON 0.000001

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

# define FEI_cal_Bc_num 12
# define FEI_cal_B0_num 11
double FEI_cal_Bc[FEI_cal_Bc_num] = { 1.04, 0.79, 0.69, 0.56, 0.97, 0.95, 0.74, 0.57, 0.91, 0.51, 0.34, 0.59 };
double FEI_cal_Bc_uncertainty[FEI_cal_Bc_num] = { 0.03, 0.03, 0.05, 0.11, 0.03, 0.03, 0.02, 0.06, 0.1, 0.13, 0.07, 0.02 }; // not relative uncertainty. absolute uncertainty
double FEI_cal_Bc_modeID[FEI_cal_Bc_num] = { 0.0, 1.0, 3.0, 4.0, 15.0, 16.0, 18.0, 19.0, 23.0, 24.0, 30.0, -1.0 };
double FEI_cal_B0[FEI_cal_B0_num] = { 1.16, 0.94, 0.81, 0.79, 0.99, 1.03, 0.67, 0.66, 0.69, 0.49, 0.79 };
double FEI_cal_B0_uncertainty[FEI_cal_B0_num] = { 0.04, 0.05, 0.06, 0.02, 0.03, 0.06, 0.02, 0.03, 0.02, 0.02, 0.12 }; // not relative uncertainty. absolute uncertainty
double FEI_cal_B0_modeID[FEI_cal_B0_num] = { 0.0, 1.0, 3.0, 4.0, 5.0, 15.0, 16.0, 18.0, 19.0, 26.0, -1.0 };

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

//# define Toy_iter_num 500
//# define LT_iter_num 500

int Toy_iter_num = 0.0;
int LT_iter_num = 0.0;

# define RarityBins 6
# define BinMIN 0.96
# define BinMAX 1.0
double weight_KIDsys[RarityBins * 7] = { 0.0 };
double weight_PIDsys[RarityBins * 7] = { 0.0 };
double weight_BRsys[RarityBins * 3] = { 0.0 };
double weight_pi0sys[RarityBins * 7] = { 0.0 };
double weight_FEIsys[RarityBins * 3] = { 0.0 };

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

        if ((type == "MC15ri") || (type == "MC15rd")) {
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
            else {
                printf("[Corrector_Fragmentation::GetCorrectionFactor] unexpected syst type!\n");
                exit(1);
                return 1;
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
        else TargetNevtAtMxsBin = Total_Nevt_Nominal_before_Xsgamma_MC15[TargetMxsBin] * (Fragmentation_Xsgamma[TargetMxsBin][TargetCategory] * 0.01);
        double TotalNevtAtMxsBinWithMissingWithoutTargetCategory = TotalNevtAtMxsBinWithMissing - TargetNevtAtMxsBin;

        if (Category == TargetCategory) { // it is target category
            if (IsTargetCategoryUp) {
                if ((TotalNevtAtMxsBinWithMissingWithoutTargetCategory - TargetNevtAtMxsBin * RelativeUncertainty) > 0) return (1.0 + RelativeUncertainty); // the number of event in target category is not large
                else return (1.0 + (TotalNevtAtMxsBinWithMissingWithoutTargetCategory / TargetNevtAtMxsBin)); // the number of event in target category is large
            }
            else return (1.0 - RelativeUncertainty);
        }
        else {
            if (IsTargetCategoryUp) {
                double output = (TotalNevtAtMxsBinWithMissingWithoutTargetCategory - TargetNevtAtMxsBin * RelativeUncertainty) / TotalNevtAtMxsBinWithMissingWithoutTargetCategory;
                if (output < 0) return 0.0; // the number of event in target category is large
                else return output; // the number of event in target category is not large
            }
            else {
                double output = (TotalNevtAtMxsBinWithMissingWithoutTargetCategory + TargetNevtAtMxsBin * RelativeUncertainty) / TotalNevtAtMxsBinWithMissingWithoutTargetCategory;
                if (output < 0) {
                    printf("[Corrector_Fragmentation::FluctuateCorrection] minus Nevt!\n");
                    exit(1);
                }
                else return output;
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

std::random_device rd;
std::default_random_engine generator(rd());

std::vector<std::string> Sample_names = {
    "L_x_Signal_nominal_channel_overallSyst_x_StatUncert_x_channel_Signal_all_uncorr_uncer_ShapeSys",
    "L_x_CHG_nominal_channel_overallSyst_x_StatUncert_x_channel_CHG_all_uncorr_uncer_ShapeSys",
    "L_x_MIX_nominal_channel_overallSyst_x_StatUncert_x_channel_MIX_all_uncorr_uncer_ShapeSys",
    "L_x_UUBAR_nominal_channel_overallSyst_x_StatUncert_x_channel_UUBAR_all_uncorr_uncer_ShapeSys",
    "L_x_DDBAR_nominal_channel_overallSyst_x_StatUncert_x_channel_DDBAR_all_uncorr_uncer_ShapeSys",
    "L_x_SSBAR_nominal_channel_overallSyst_x_StatUncert_x_channel_SSBAR_all_uncorr_uncer_ShapeSys",
    "L_x_CHARM_nominal_channel_overallSyst_x_StatUncert_x_channel_CHARM_all_uncorr_uncer_ShapeSys"
};

typedef struct Options
{
    // turn on it?
    bool track;
    bool PID;
    bool KID;
    bool KS0;
    bool pi0;
    bool FBDT;
    bool FEI;
    bool qqbar;
    bool multiplicity;
    bool Kff;
    bool Kstarff;
    bool OLD_Kff;
    bool pf;
    bool Transition;
    bool mb;
    bool fraction;
    bool MCstat;
    bool Fragmentation;
    bool mKstar;
    bool BKGNorm;
    bool BDTc;
    bool BBcounting;
    bool BBBR;
    bool BRBtoKKLKL;
    bool BRKnn;
    bool BRDKL0;
    bool uncorrelated;

    int NEntryFEI;
    int NEntryKID;
    int NEntryPID;
    int NEntryBR;
    int NEntrypi0;
    int NEntryMultiplicity;

    int Nsyst = 25;
} OPTIONS;

void GetNameOfParams(RooWorkspace* w, std::vector<std::string>* names) {
    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    w->loadSnapshot("ParamValues");

    RooAbsData* data = (RooAbsData*)w->data("asimovData");
    RooFitResult* fitres = model->fitTo(*data, RooFit::Extended(true), RooFit::SumW2Error(false), PrintLevel(-1), Save());

    RooArgSet fitargs = fitres->floatParsFinal();
    TIterator* iter(fitargs.createIterator());

    for (TObject* a = iter->Next(); a != 0; a = iter->Next()) {
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(a);
        std::string name = rrv->GetName();
        names->push_back(name);
    }

    w->loadSnapshot("ParamValues");
}

class FileSaver {
private:
    double* m_true_param;
    double* m_fitting_param;
    double* m_fitting_param_error;
    double* m_fitting_param_HIerror;
    double* m_fitting_param_LOerror;
    int m_covQual;
    int m_status;
    double m_edm;

    TFile* m_file;
    TTree* m_tree;

public:
    FileSaver();
    void OpenFile(bool IsItToy, std::vector<std::string>* names, double mu_injected, int indicator);
    void CloseFile();
    void GetTrueValues(RooWorkspace* w, std::vector<std::string>* names);
    void GetFittingValues(RooFitResult* fitres, std::vector<std::string>* names);
    void GetFittingStatus(RooFitResult* fitres);
    void WriteIntoBranch();
};

FileSaver::FileSaver() {
    m_true_param = nullptr;
    m_fitting_param = nullptr;
    m_fitting_param_error = nullptr;
    m_fitting_param_HIerror = nullptr;
    m_fitting_param_LOerror = nullptr;
    m_covQual = -1;
    m_status = -1;
    m_edm = -1;
    TFile* m_file = nullptr;
    TTree* m_tree = nullptr;
}

void FileSaver::OpenFile(bool IsItToy, std::vector<std::string>* names, double mu_injected, int indicator) {

    // malloc params
    m_true_param = (double*)malloc(sizeof(double) * names->size());
    m_fitting_param = (double*)malloc(sizeof(double) * names->size());
    m_fitting_param_error = (double*)malloc(sizeof(double) * names->size());
    m_fitting_param_HIerror = (double*)malloc(sizeof(double) * names->size());
    m_fitting_param_LOerror = (double*)malloc(sizeof(double) * names->size());

    // open file
    if (IsItToy) {
        m_file = new TFile(("./out/TOY_result_" + std::to_string(indicator) + ".root").c_str(), "recreate");
        m_file->cd();
        m_tree = new TTree("TOY_result", "");
    }
    else {
        m_file = new TFile(("./out/LT_result_" + std::to_string(mu_injected) + "_" + std::to_string(indicator) + ".root").c_str(), "recreate");
        m_file->cd();
        m_tree = new TTree("LT_result", "");
    }

    // set branches
    for (unsigned int i = 0; i < names->size(); i++) {
        m_tree->Branch((names->at(i) + "_true").c_str(), &m_true_param[i]);
        m_tree->Branch((names->at(i) + "_value").c_str(), &m_fitting_param[i]);
        m_tree->Branch((names->at(i) + "_error").c_str(), &m_fitting_param_error[i]);
        m_tree->Branch((names->at(i) + "_HIerror").c_str(), &m_fitting_param_HIerror[i]);
        m_tree->Branch((names->at(i) + "_LOerror").c_str(), &m_fitting_param_LOerror[i]);
    }
    m_tree->Branch("covQual", &m_covQual);
    m_tree->Branch("status", &m_status);
    m_tree->Branch("edm", &m_edm);
}

void FileSaver::CloseFile() {
    m_file->cd();
    m_tree->Write();
    m_file->Close();
}

void FileSaver::GetTrueValues(RooWorkspace* w, std::vector<std::string>* names) {

    // write into member variables
    for (unsigned int i = 0; i < names->size(); i++) {
        m_true_param[i] = w->var(names->at(i).c_str())->getValV();
    }

}

void FileSaver::GetFittingValues(RooFitResult* fitres, std::vector<std::string>* names) {

    RooArgSet fitargs = fitres->floatParsFinal();
    TIterator* iter(fitargs.createIterator());

    for (TObject* a = iter->Next(); a != 0; a = iter->Next()) {
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(a);
        std::string name = rrv->GetName();
        double val = rrv->getVal();
        double err = rrv->getError();
        double HIerr = rrv->getAsymErrorHi();
        double LOerr = rrv->getAsymErrorLo();

        for (unsigned int i = 0; i < names->size(); i++) {
            if (name == names->at(i)) {
                m_fitting_param[i] = val;
                m_fitting_param_error[i] = err;
                m_fitting_param_HIerror[i] = HIerr;
                m_fitting_param_LOerror[i] = LOerr;
            }
        }

    }

}

void FileSaver::GetFittingStatus(RooFitResult* fitres) {

    m_covQual = fitres->covQual();
    m_status = fitres->status();
    m_edm = fitres->edm();

}

void FileSaver::WriteIntoBranch() {

    m_tree->Fill();

}

FileSaver filesaver;

RooFitResult* MyMinimizeNLL(RooWorkspace* w, RooDataSet* data, RooAbsReal* nll, double tolerance = -1.0) {
    // what we have done
    w->loadSnapshot("ParamValues");
    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    // get nll
    RooArgSet* allParams = model->getParameters(*data);
    RooStats::RemoveConstantParameters(allParams);
    RooArgSet fGlobalObs = *mc->GetGlobalObservables();
    RooArgSet fConditionalObs;
    Bool_t fLOffset = RooStats::IsNLLOffset();
    nll = model->createNLL(*data, RooFit::CloneData(kFALSE), RooFit::Constrain(*allParams), RooFit::GlobalObservables(fGlobalObs), RooFit::ConditionalObservables(fConditionalObs), RooFit::Offset(fLOffset));

    // minimizer option
    TString fMinimizer = ::ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str();
    TString minimizer = fMinimizer;

    TString algorithm = ::ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo();

    Int_t fStrategy = ::ROOT::Math::MinimizerOptions::DefaultStrategy();

    Double_t fTolerance;
    if (tolerance < 0) fTolerance = TMath::Max(1., ::ROOT::Math::MinimizerOptions::DefaultTolerance());
    else fTolerance = tolerance;

    Int_t fPrintLevel = ::ROOT::Math::MinimizerOptions::DefaultPrintLevel();
    //LM: RooMinimizer.setPrintLevel has +1 offset - so subtract  here -1 + an extra -1
    int level = (fPrintLevel == 0) ? -1 : fPrintLevel - 2;


    // follow what ProfileLikelihoodTestStat.cxx does
    const auto& config = RooStats::GetGlobalRooStatsConfig();
    RooMinimizer minim(*nll);
    minim.setStrategy(fStrategy);
    minim.setEvalErrorWall(config.useEvalErrorWall);
    minim.setEps(fTolerance);
    minim.setPrintLevel(level);
    // this causes a memory leak
    minim.optimizeConst(2);
    minim.migrad();
    minim.minos(RooArgSet(*w->var("mu")));

    // fit!
    int status;
    status = minim.minimize(minimizer, algorithm);

    return minim.save();
}

RooFitResult* MinimizeNLL(RooWorkspace* w, RooDataSet* data, RooAbsReal* nll, double tolerance = -1.0) { // this function follows the procedure in ProfileLikelihoodTestStat.cxx
    // what we have done
    // deprecated
    w->loadSnapshot("ParamValues");
    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    // get nll
    RooArgSet* allParams = model->getParameters(*data);
    RooStats::RemoveConstantParameters(allParams);
    RooArgSet fGlobalObs;
    RooArgSet fConditionalObs;
    Bool_t fLOffset = RooStats::IsNLLOffset();
    nll = model->createNLL(*data, RooFit::CloneData(kFALSE), RooFit::Constrain(*allParams), RooFit::GlobalObservables(fGlobalObs), RooFit::ConditionalObservables(fConditionalObs), RooFit::Offset(fLOffset));

    // get default values. follow ProfileLikelihoodTestStat.h
    TString fMinimizer = ::ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str();
    Int_t fStrategy = ::ROOT::Math::MinimizerOptions::DefaultStrategy();
    // avoid default tolerance to be too small (1. is default in RooMinimizer)
    Double_t fTolerance;
    if (tolerance < 0) fTolerance = TMath::Max(1., ::ROOT::Math::MinimizerOptions::DefaultTolerance());
    else fTolerance = tolerance;
    Int_t fPrintLevel = ::ROOT::Math::MinimizerOptions::DefaultPrintLevel();

    // follow what ProfileLikelihoodTestStat.cxx does
    const auto& config = RooStats::GetGlobalRooStatsConfig();
    RooMinimizer minim(*nll);
    minim.setStrategy(fStrategy);
    minim.setEvalErrorWall(config.useEvalErrorWall);
    //LM: RooMinimizer.setPrintLevel has +1 offset - so subtract  here -1 + an extra -1
    int level = (fPrintLevel == 0) ? -1 : fPrintLevel - 2;
    minim.setPrintLevel(level);
    minim.setEps(fTolerance);
    // this causes a memory leak
    minim.optimizeConst(2);
    TString minimizer = fMinimizer;
    TString algorithm = ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo();
    if (algorithm == "Migrad") algorithm = "Minimize"; // prefer to use Minimize instead of Migrad
    int status;
    for (int tries = 1, maxtries = 4; tries <= maxtries; ++tries) {
        status = minim.minimize(minimizer, algorithm);
        if (status % 1000 == 0) {  // ignore erros from Improve
            break;
        }
        else if (tries < maxtries) {
            cout << "    ----> Doing a re-scan first" << endl;
            minim.minimize(minimizer, "Scan");
            if (tries == 2) {
                if (fStrategy == 0) {
                    cout << "    ----> trying with strategy = 1" << endl;;
                    minim.setStrategy(1);
                }
                else
                    tries++; // skip this trial if strategy is already 1
            }
            if (tries == 3) {
                cout << "    ----> trying with improve" << endl;;
                minimizer = "Minuit";
                algorithm = "migradimproved";
            }
        }
    }
    /*
    if (MyDEBUG) {

        RooArgSet fitargs = minim.save()->floatParsFinal();
        TIterator* iter(fitargs.createIterator());

        for (TObject* a = iter->Next(); a != 0; a = iter->Next()) {
            RooRealVar* rrv = dynamic_cast<RooRealVar*>(a);
            std::string name = rrv->GetName();
            double val = rrv->getVal();
            double err = rrv->getError();

            if (name == "mu") {
                if (err < 1.5) {
                    RooRealVar* poi = (RooRealVar*)mc->GetParametersOfInterest()->first();
                    RooPlot* x_frame = poi->frame(Title("profile LL"));
                    //nll->plotOn(x_frame);
                    RooAbsReal* pll = nll->createProfile(*poi);
                    pll->plotOn(x_frame);

                    TCanvas* canvas = new TCanvas("sPlot", "sPlot demo", 700, 700);

                    // draw PDFs
                    x_frame->Draw();

                    canvas->SaveAs((std::string("pll_plot_") + "mu-" + std::to_string(val) + "_" + "err-" + std::to_string(err) + ".png").c_str());

                    delete canvas;
                }

            }
        }
    }
    */
    return minim.save();
}

std::vector<std::string> split(std::string str, char Delimiter) {
    std::istringstream iss(str);
    std::string buffer;

    std::vector<std::string> result;

    while (getline(iss, buffer, Delimiter)) {
        result.push_back(buffer);
    }

    return result;
}

double SetParamsForToy(RooWorkspace* w, std::vector<std::string>* names, double injected_mu) {

    double Nevt = 0.0;

    w->loadSnapshot("ParamValues");
    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();
    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x_val = w->var("obs_x_channel");

    for (unsigned int i = 0; i < names->size(); i++) {
        if (names->at(i).find("alpha") != std::string::npos) {
            std::normal_distribution<double> distribution(0.0, 1.0);
            w->var(names->at(i).c_str())->setVal(distribution(generator));
        }
        else if (names->at(i).find("gamma_stat") != std::string::npos) {
            RooRealVar* norm = w->var(("nom_" + names->at(i)).c_str());
            std::poisson_distribution<> distribution(norm->getValV());
            w->var(names->at(i).c_str())->setVal(distribution(generator) / norm->getValV());
        }
        else if ( (names->at(i).find("gamma") != std::string::npos) && (names->at(i).find("uncorr") != std::string::npos) ) {
            int sample_index = -1;
            int bin_index = -1;

            if (names->at(i).find("CHG") != std::string::npos) sample_index = 0;
            else if(names->at(i).find("MIX") != std::string::npos) sample_index = 1;
            else if (names->at(i).find("UUBAR") != std::string::npos) sample_index = 2;
            else if (names->at(i).find("DDBAR") != std::string::npos) sample_index = 3;
            else if (names->at(i).find("SSBAR") != std::string::npos) sample_index = 4;
            else if (names->at(i).find("CHARM") != std::string::npos) sample_index = 5;
            else if (names->at(i).find("Signal") != std::string::npos) sample_index = 6;

            std::vector<std::string> temp_strings = split(names->at(i), '_');
            bin_index = stoi(temp_strings.back()); // from 0

            if (names->at(i).find("all") != std::string::npos) {
                double KID_uncertainty = weight_KIDsys[RarityBins * sample_index + bin_index];
                double PID_uncertainty = weight_PIDsys[RarityBins * sample_index + bin_index];
                double BR_uncertainty = 0.0;
                double pi0_uncertainty = weight_pi0sys[RarityBins * sample_index + bin_index];
                double FEI_uncertainty = 0.0;
                if ((names->at(i).find("CHG") != std::string::npos) || (names->at(i).find("MIX") != std::string::npos)) {
                    BR_uncertainty = weight_BRsys[RarityBins * sample_index + bin_index];
                    FEI_uncertainty = weight_FEIsys[RarityBins * sample_index + bin_index];
                }
                else if (names->at(i).find("Signal") != std::string::npos) {
                    BR_uncertainty = weight_BRsys[RarityBins * 2 + bin_index]; // exception for signal BB BR uncorrelated uncertainty!
                    FEI_uncertainty = weight_FEIsys[RarityBins * 2 + bin_index]; // exception for signal FEI uncorrelated uncertainty!
                }
                
                double total_uncertainty = std::sqrt(KID_uncertainty * KID_uncertainty + PID_uncertainty * PID_uncertainty + BR_uncertainty * BR_uncertainty + pi0_uncertainty * pi0_uncertainty + FEI_uncertainty * FEI_uncertainty);
                std::normal_distribution<double> distribution(1.0, total_uncertainty);
                w->var(names->at(i).c_str())->setVal(distribution(generator));
            }

        }
    }

    w->var("mu")->setVal(injected_mu);

    /* ================================ cal Nexpected ================================*/
    RooAbsBinning const& binning = x_val->getBinning();
    const double oldVal = x_val->getVal();

    for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
        double binCenter = binning.binCenter(iBin);
        double binWidth = binning.binWidth(iBin);

        *x_val = binCenter; // set x value

        for (unsigned int j = 0; j < Sample_names.size(); j++) {
            RooAbsReal* temp_func = w->function(Sample_names.at(j).c_str());
            Nevt = Nevt + temp_func->getValV();
            if (temp_func->getValV() < 0) {
                printf("[ERROR] negative count!\n");
                exit(1);
            }
        }

    }

    *x_val = oldVal;

    return Nevt;

}

void MyToyMCStudy(RooWorkspace *w, std::vector<std::string>* names, double eps, int indicator = 0){

        ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
        RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

        RooArgSet* obs = (RooArgSet*)mc->GetObservables();
        RooRealVar* x = (RooRealVar*)obs->find("obs_x_channel");

        for(int i=0; i< Toy_iter_num; i++) { // Do Toy MC study
            
            double Nevt_total = SetParamsForToy(w, names, 1.0);

            filesaver.GetTrueValues(w, names);

            RooDataSet* genData = model->generate(RooArgSet(*x,model->indexCat()), Nevt_total, false, true, "", false, true);

            w->loadSnapshot("ParamValues");
            //RooFitResult* fitres = model->fitTo(*genData, RooFit::Extended(false), RooFit::Minos(RooArgSet(*w->var("mu"))), RooFit::SumW2Error(false), PrintLevel(-1), Save());
            RooAbsReal* nll;
            RooFitResult* fitres = MyMinimizeNLL(w, genData, nll, eps);

            if(MyDEBUG) Debug(w, fitres, genData);

            filesaver.GetFittingValues(fitres, names);
            filesaver.GetFittingStatus(fitres);
            filesaver.WriteIntoBranch();

            delete fitres;

        }
}

void MyLinearityTest(RooWorkspace* w, std::vector<std::string>* names, double mu_injected, double eps, int indicator = 0) {

    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x = (RooRealVar*)obs->find("obs_x_channel");

    for (int i = 0; i < LT_iter_num; i++) { // Do LT MC study

        double Nevt_total = SetParamsForToy(w, names, mu_injected);

        filesaver.GetTrueValues(w, names);

        RooDataSet* genData = model->generate(RooArgSet(*x, model->indexCat()), Nevt_total, false, true, "", false, true);
        w->loadSnapshot("ParamValues");

        //RooFitResult* fitres = model->fitTo(*genData, RooFit::Extended(false), RooFit::Minos(RooArgSet(*w->var("mu"))), RooFit::SumW2Error(false), PrintLevel(-1), Save());
        RooAbsReal* nll;
        RooFitResult* fitres = MyMinimizeNLL(w, genData, nll, eps);

        filesaver.GetFittingValues(fitres, names);
        filesaver.GetFittingStatus(fitres);
        filesaver.WriteIntoBranch();

        delete fitres;

    }
}

void ReadPIDuncorrsysFile(const char* dirname_KID, const char* dirname_PID) {
    FILE* fp;

    fp = fopen(dirname_KID, "r");
    for (int i = 0; i < RarityBins * 7; i++) fscanf(fp, "%lf\n", &weight_KIDsys[i]);
    fclose(fp);
    for (int i = 0; i < RarityBins * 7; i++) weight_KIDsys[i] = std::sqrt(weight_KIDsys[i]);

    fp = fopen(dirname_PID, "r");
    for (int i = 0; i < RarityBins * 7; i++) fscanf(fp, "%lf\n", &weight_PIDsys[i]);
    fclose(fp);
    for (int i = 0; i < RarityBins * 7; i++) weight_PIDsys[i] = std::sqrt(weight_PIDsys[i]);
}

void ReadBRuncorrsysFile(const char* dirname_BR) {
    FILE* fp;

    fp = fopen(dirname_BR, "r");
    for (int i = 0; i < RarityBins * 3; i++) fscanf(fp, "%lf\n", &weight_BRsys[i]);
    fclose(fp);
    for (int i = 0; i < RarityBins * 3; i++) weight_BRsys[i] = std::sqrt(weight_BRsys[i]);

}

void Readpi0uncorrsysFile(const char* dirname_pi0) {
    FILE* fp;

    fp = fopen(dirname_pi0, "r");
    for (int i = 0; i < RarityBins * 7; i++) fscanf(fp, "%lf\n", &weight_pi0sys[i]);
    fclose(fp);
    for (int i = 0; i < RarityBins * 7; i++) weight_pi0sys[i] = std::sqrt(weight_pi0sys[i]);

}

void ReadFEIuncorrsysFile(const char* dirname_FEI) {
    FILE* fp;

    fp = fopen(dirname_FEI, "r");
    for (int i = 0; i < RarityBins * 3; i++) fscanf(fp, "%lf\n", &weight_FEIsys[i]);
    fclose(fp);
    for (int i = 0; i < RarityBins * 3; i++) weight_FEIsys[i] = std::sqrt(weight_FEIsys[i]);

}

double GetNumEvts(RooWorkspace* w, const char* sample_type) {

    double Nevt = 0.0;

    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();
    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x_val = w->var("obs_x_channel");

    int index = -1;
    if (strcmp(sample_type, "Signal") == 0) index = 0;
    else if (strcmp(sample_type, "CHG") == 0) index = 1;
    else if (strcmp(sample_type, "MIX") == 0) index = 2;
    else if (strcmp(sample_type, "UUBAR") == 0) index = 3;
    else if (strcmp(sample_type, "DDBAR") == 0) index = 4;
    else if (strcmp(sample_type, "SSBAR") == 0) index = 5;
    else if (strcmp(sample_type, "CHARM") == 0) index = 6;
    else {
        printf("[ERROR] unexpected sample type!\n");
        exit(1);
    }

    /* ================================ cal Nexpected ================================*/
    RooAbsBinning const& binning = x_val->getBinning();
    const double oldVal = x_val->getVal();

    for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
        double binCenter = binning.binCenter(iBin);
        double binWidth = binning.binWidth(iBin);

        *x_val = binCenter; // set x value

        RooAbsReal* temp_func = w->function(Sample_names.at(index).c_str());
        Nevt = Nevt + temp_func->getValV();
        if (temp_func->getValV() < 0) {
            printf("[ERROR] negative count!\n");
            exit(1);
        }

    }

    *x_val = oldVal;

    return Nevt;

}

void FitToData(RooWorkspace* w, double eps) {

    w->loadSnapshot("ParamValues");

    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();
    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x_val = w->var("obs_x_channel");

    // get Category and data
    RooCategory* idx = (RooCategory*)obs->find("channelCat");
    RooDataSet* data = (RooDataSet*)w->data("asimovData");

    // fit
    //RooFitResult* fitres = model->fitTo(*data, RooFit::Extended(false), RooFit::Minos(RooArgSet(*w->var("mu"))), RooFit::SumW2Error(false), PrintLevel(-1), Save());
    RooAbsReal* nll;
    RooFitResult* fitres = MyMinimizeNLL(w, data, nll, eps);

    // get expected num of evts for PDFs
    double Signal_Nevts = GetNumEvts(w, "Signal");
    double CHG_Nevts = GetNumEvts(w, "CHG");
    double MIX_Nevts = GetNumEvts(w, "MIX");
    double UUBAR_Nevts = GetNumEvts(w, "UUBAR");
    double DDBAR_Nevts = GetNumEvts(w, "DDBAR");
    double SSBAR_Nevts = GetNumEvts(w, "SSBAR");
    double CHARM_Nevts = GetNumEvts(w, "CHARM");

    // draw
    RooPlot* x_frame = x_val->frame(Title("FBDT"));
    data->plotOn(x_frame, DataError(RooAbsData::Poisson), Cut("channelCat==0"), DrawOption("ZP"), Name("data"));
    //data->plotOn(x_frame, DataError(RooAbsData::Poisson), Cut("channelCat==0"), MarkerSize(0.4), DrawOption("ZP"), Normalization(1, RooAbsReal::ScaleType::NumEvent));
    model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kRed - 6), LineWidth(0), Components("L_x_*Signal*_ShapeSys,L_x_*CHG*_ShapeSys,L_x_*MIX*_ShapeSys,L_x_*UUBAR*_ShapeSys,L_x_*DDBAR*_ShapeSys,L_x_*SSBAR*_ShapeSys,L_x_*CHARM*_ShapeSys"), Normalization(Signal_Nevts + CHG_Nevts + MIX_Nevts + UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("signal"));
    model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kBlue - 6), LineWidth(0), Components("L_x_*CHG*_ShapeSys,L_x_*MIX*_ShapeSys,L_x_*UUBAR*_ShapeSys,L_x_*DDBAR*_ShapeSys,L_x_*SSBAR*_ShapeSys,L_x_*CHARM*_ShapeSys"), Normalization(CHG_Nevts + MIX_Nevts + UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("Charged B"));
    model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kCyan - 6), LineWidth(0), Components("L_x_*MIX*_ShapeSys,L_x_*UUBAR*_ShapeSys,L_x_*DDBAR*_ShapeSys,L_x_*SSBAR*_ShapeSys,L_x_*CHARM*_ShapeSys"), Normalization(MIX_Nevts + UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("Neutral B"));
    model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kOrange - 6), LineWidth(0), Components("L_x_*UUBAR*_ShapeSys,L_x_*DDBAR*_ShapeSys,L_x_*SSBAR*_ShapeSys,L_x_*CHARM*_ShapeSys"), Normalization(UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("Continuum"));
    data->plotOn(x_frame, DataError(RooAbsData::Poisson), Cut("channelCat==0"), DrawOption("ZP"), Name("data"));

    TCanvas* canvas = new TCanvas("sPlot", "sPlot demo", 700, 700);

    // draw PDFs
    x_frame->Draw();

    // draw legend
    TLegend* leg = new TLegend(0.9, 0.9, 0.7, 0.7);
    leg->SetFillStyle(0);
    leg->SetLineWidth(0);
    leg->AddEntry("data", "Data");
    leg->AddEntry("signal", "B #rightarrow X_{s} #nu #bar{#nu}");
    leg->AddEntry("Charged B", "Charged B");
    leg->AddEntry("Neutral B", "Neutral B");
    leg->AddEntry("Continuum", "q#bar{q}");
    leg->Draw();

    canvas->SaveAs("fit_plot.png");

    // Get fitting variables
    RooArgSet fitargs = fitres->floatParsFinal();
    TIterator* iter(fitargs.createIterator());

    for (TObject* a = iter->Next(); a != 0; a = iter->Next()) {
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(a);
        std::string name = rrv->GetName();
        double val = rrv->getVal();
        double err = rrv->getError();

        if (name == std::string("mu")) {
            printf("mu: %lf +- %lf\n", val, err);
        }

    }

    w->loadSnapshot("ParamValues");

    delete canvas;
    delete fitres;
}

void MyToyMCStudyDataPoisson(RooWorkspace* w, std::vector<std::string>* names, double eps) {

    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x = (RooRealVar*)obs->find("obs_x_channel");

    // construct PDfs from data
    RooAbsData* data = (RooAbsData*)w->data("asimovData");
    TH1D* HistData = (TH1D*)data->createHistogram("HistData", *x, Binning(RarityBins, BinMIN, BinMAX));
    RooDataHist DataHist("DataHist", "DataHist", *x, HistData);
    RooHistPdf dataPDF("dataPDF", "dataPDF", RooArgSet(*x, model->indexCat()), DataHist);

    // get total number of event in data
    double Nevt_data = 0;
    int Nevt_data_int = 0;
    for (int i = 0; i < RarityBins; i++) Nevt_data = Nevt_data + HistData->GetBinContent(i + 1);
    Nevt_data_int = (int)std::round(Nevt_data);

    
    for (int i = 0; i < Toy_iter_num; i++) { // Do Toy MC study

        w->loadSnapshot("ParamValues");

        filesaver.GetTrueValues(w, names);

        RooDataSet* genData = dataPDF.generate(RooArgSet(*x, model->indexCat()), Nevt_data_int, false, true, "", false, true);

        w->loadSnapshot("ParamValues");
        //RooFitResult* fitres = model->fitTo(*genData, RooFit::Extended(false), RooFit::Minos(RooArgSet(*w->var("mu"))), RooFit::SumW2Error(false), PrintLevel(-1), Save());
        RooAbsReal* nll;
        RooFitResult* fitres = MyMinimizeNLL(w, genData, nll, eps);

        filesaver.GetFittingValues(fitres, names);
        filesaver.GetFittingStatus(fitres);
        filesaver.WriteIntoBranch();

        delete fitres;

    }
    
}

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
    fscanf(fp, "%d", &N_hist);
    fclose(fp);

    return N_hist;
}

void Initialize_options(OPTIONS* options_, const char* tested_param) {
    // initialize
    options_->track = false;
    options_->PID = false;
    options_->KID = false;
    options_->KS0 = false;
    options_->pi0 = false;
    options_->FBDT = false;
    options_->FEI = false;
    options_->qqbar = false;
    options_->multiplicity = false;
    options_->Kff = false;
    options_->Kstarff = false;
    options_->OLD_Kff = false;
    options_->pf = false;
    options_->Transition = false;
    options_->mb = false;
    options_->fraction = false;
    options_->MCstat = false;
    options_->Fragmentation = false;
    options_->mKstar = false;
    options_->BKGNorm = false;
    options_->BDTc = false;
    options_->BBcounting = false;
    options_->BBBR = false;
    options_->BRBtoKKLKL = false;
    options_->BRKnn = false;
    options_->BRDKL0 = false;
    options_->uncorrelated = false;

    if (std::string(tested_param) == std::string("all")) {
        options_->track = true;
        options_->PID = true;
        options_->KID = true;
        options_->KS0 = true;
        options_->pi0 = true;
        options_->FBDT = true;
        options_->FEI = true;
        options_->qqbar = true;
        options_->multiplicity = true;
        options_->Kff = true;
        options_->Kstarff = true;
        options_->OLD_Kff = true;
        options_->pf = true;
        options_->Transition = true;
        options_->mb = true;
        options_->fraction = true;
        options_->MCstat = true;
        options_->Fragmentation = true;
        options_->mKstar = true;
        options_->BKGNorm = true;
        options_->BDTc = true;
        options_->BBcounting = true;
        options_->BBBR = true;
        options_->BRBtoKKLKL = true;
        options_->BRKnn = true;
        options_->BRDKL0 = true;
        options_->uncorrelated = true;
    }
    else if (std::string(tested_param) == std::string("none")) {}
    else if (std::string(tested_param) == std::string("track")) options_->track = true;
    else if (std::string(tested_param) == std::string("PID")) options_->PID = true;
    else if (std::string(tested_param) == std::string("KID")) options_->KID = true;
    else if (std::string(tested_param) == std::string("KS0")) options_->KS0 = true;
    else if (std::string(tested_param) == std::string("pi0")) options_->pi0 = true;
    else if (std::string(tested_param) == std::string("FBDT_efficiency")) options_->FBDT = true;
    else if (std::string(tested_param) == std::string("FEI")) options_->FEI = true;
    else if (std::string(tested_param) == std::string("qqbar_CAL")) options_->qqbar = true;
    else if (std::string(tested_param) == std::string("photon_multiplicity")) options_->multiplicity = true;
    else if (std::string(tested_param) == std::string("Kff")) options_->Kff = true;
    else if (std::string(tested_param) == std::string("Kstarff")) options_->Kstarff = true;
    else if (std::string(tested_param) == std::string("OLDKff")) options_->OLD_Kff = true;
    else if (std::string(tested_param) == std::string("pf")) options_->pf = true;
    else if (std::string(tested_param) == std::string("Transition")) options_->Transition = true;
    else if (std::string(tested_param) == std::string("mb")) options_->mb = true;
    else if (std::string(tested_param) == std::string("fraction")) options_->fraction = true;
    else if (std::string(tested_param) == std::string("MCstat")) options_->MCstat = true;
    else if (std::string(tested_param) == std::string("Fragmentation")) options_->Fragmentation = true;
    else if (std::string(tested_param) == std::string("mKstar")) options_->mKstar = true;
    else if (std::string(tested_param) == std::string("BKGNorm")) options_->BKGNorm = true;
    else if (std::string(tested_param) == std::string("BDTc")) options_->BDTc = true;
    else if (std::string(tested_param) == std::string("BBcounting")) options_->BBcounting = true;
    else if (std::string(tested_param) == std::string("BBBR")) options_->BBBR = true;
    else if (std::string(tested_param) == std::string("BRBtoKKLKL")) options_->BRBtoKKLKL = true;
    else if (std::string(tested_param) == std::string("BRKnn")) options_->BRKnn = true;
    else if (std::string(tested_param) == std::string("BtoDtoXKL")) options_->BRDKL0 = true;
    else if (std::string(tested_param) == std::string("uncorrelated")) options_->uncorrelated = true;
    else {
        printf("inappropriate parameter name\n");
        exit(1);
    }

    // read entry for nuisance parameters
    options_->NEntryFEI = ReadNFEIEigenVector("./FEI_selected.txt");
    options_->NEntryKID = ReadNPIDEigenVector("./KID_selected.txt");
    options_->NEntryPID = ReadNPIDEigenVector("./PID_selected.txt");
    options_->NEntryBR = ReadNBREigenVector("./BR_selected.txt");
    options_->NEntrypi0 = ReadNpi0EigenVector("./pi0_selected.txt");
    options_->NEntryMultiplicity = ReadMultiplicityInfo("./Multiplicity_info.txt");
}

void FixParameters(RooWorkspace* w, OPTIONS* options_) {

    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    w->loadSnapshot("NominalParamValues");

    RooRealVar* x_val = w->var("obs_x_channel");
    std::unique_ptr<RooArgSet> params{model->getParameters(*x_val)};

    // track
    w->var("alpha_track_eff_uncer")->setConstant(options_->track);

    // PID
    for (int i = 0; i < options_->NEntryPID; i++) w->var(("alpha_PID" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->PID);

    // KID
    for (int i = 0; i < options_->NEntryKID; i++) w->var(("alpha_KID" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->KID);

    // KS0
    w->var("alpha_KS0_reco_uncer")->setConstant(options_->KS0);

    // pi0
    for (int i = 0; i < options_->NEntrypi0; i++) w->var(("alpha_pi0" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->pi0);

    // FBDT efficiency
    w->var("alpha_FBDT_efficiency_uncer")->setConstant(options_->FBDT);

    // FEI
    for (int i = 0; i < options_->NEntryFEI; i++) w->var(("alpha_FEI" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->FEI);

    // qqbar normalization
    w->var("alpha_qq_CAL_UUBAR_uncer")->setConstant(options_->qqbar);
    w->var("alpha_qq_CAL_DDBAR_uncer")->setConstant(options_->qqbar);
    w->var("alpha_qq_CAL_SSBAR_uncer")->setConstant(options_->qqbar);
    w->var("alpha_qq_CAL_CHARM_uncer")->setConstant(options_->qqbar);

    // photon multiplicity correction
    for (int i = 0; i < options_->NEntryMultiplicity; i++) w->var(("alpha_multiplicity" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->multiplicity);

    // B->K form factor
    w->var("alpha_Kff1_uncer")->setConstant(options_->Kff);
    w->var("alpha_Kff2_uncer")->setConstant(options_->Kff);
    w->var("alpha_Kff3_uncer")->setConstant(options_->Kff);

    // B->Kstar form factor
    w->var("alpha_Kstarff1_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff2_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff3_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff4_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff5_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff6_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff7_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff8_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff9_uncer")->setConstant(options_->Kstarff);

    // new B->K form factor
    w->var("alpha_Kff_OLD_uncer")->setConstant(options_->OLD_Kff);

    // fermi motion moment
    w->var("alpha_pf_uncer")->setConstant(options_->pf);

    // K*-Xs transition point
    w->var("alpha_transition_uncer")->setConstant(options_->Transition);

    // b-quark mass
    w->var("alpha_mb_uncer")->setConstant(options_->mb);

    // relative fraction
    w->var("alpha_Kfrac_uncer")->setConstant(options_->fraction);
    w->var("alpha_Kstarfrac_uncer")->setConstant(options_->fraction);

    // MC statistics
    for (int i = 0; i < RarityBins; i++) w->var(("gamma_stat_channel_bin_" + std::to_string(i)).c_str())->setConstant(options_->MCstat);

    // Fragmentation
    for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
        for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
            int temp_index = MxsBin * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma) + Category;
            w->var(("alpha_Xs_fragmentation" + std::to_string(temp_index) + "_uncer").c_str())->setConstant(options_->Fragmentation);
        }
    }

    // mKstar
    w->var("alpha_mKstar_uncer")->setConstant(options_->mKstar);

    // background normalization
    w->var("alpha_mu_CHG")->setConstant(options_->BKGNorm);
    w->var("alpha_mu_MIX")->setConstant(options_->BKGNorm);
    w->var("alpha_mu_UUBAR")->setConstant(options_->BKGNorm);
    w->var("alpha_mu_DDBAR")->setConstant(options_->BKGNorm);
    w->var("alpha_mu_SSBAR")->setConstant(options_->BKGNorm);
    w->var("alpha_mu_CHARM")->setConstant(options_->BKGNorm);

    // BDTc
    w->var("alpha_BDTc_shape_qq")->setConstant(options_->BDTc);

    // BB counting
    w->var("alpha_BB_counting_uncer")->setConstant(options_->BBcounting);

    // BB BR
    for (int i = 0; i < options_->NEntryBR; i++) w->var(("alpha_BBBR" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->BBBR);

    // B->K KL KL BR
    w->var("alpha_BRBtoKKLKL_uncer")->setConstant(options_->BRBtoKKLKL);

    // B->Knn BR
    w->var("alpha_Knn_BR_uncer")->setConstant(options_->BRKnn);
    w->var("alpha_Kstarnn_BR_uncer")->setConstant(options_->BRKnn);
    w->var("alpha_K0nn_BR_uncer")->setConstant(options_->BRKnn);
    w->var("alpha_K0starnn_BR_uncer")->setConstant(options_->BRKnn);

    // B->(D->X KL0) + anything
    w->var("alpha_BtoDtoXKL_uncer")->setConstant(options_->BRDKL0);

    // uncorrelated
    for (int i = 0; i < RarityBins; i++) {
        w->var(("gamma_Signal_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
        w->var(("gamma_CHG_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
        w->var(("gamma_MIX_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
        w->var(("gamma_UUBAR_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
        w->var(("gamma_DDBAR_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
        w->var(("gamma_SSBAR_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
        w->var(("gamma_CHARM_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
    }

    // save snapshot
    w->saveSnapshot("ParamValues", *params, true);
}

void Debug(RooWorkspace* w, RooFitResult* fitres, RooDataSet* data) {

    RooArgSet fitargs = fitres->floatParsFinal();
    TIterator* iter(fitargs.createIterator());

    for (TObject* a = iter->Next(); a != 0; a = iter->Next()) {
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(a);
        std::string name = rrv->GetName();
        double val = rrv->getVal();
        double err = rrv->getError();

        if (name == "mu") {
            if (err < 1.5) {

                ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
                RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();
                RooArgSet* obs = (RooArgSet*)mc->GetObservables();
                RooRealVar* x_val = w->var("obs_x_channel");

                // get Category and data
                RooCategory* idx = (RooCategory*)obs->find("channelCat");

                // get expected num of evts for PDFs
                double Signal_Nevts = GetNumEvts(w, "Signal");
                double CHG_Nevts = GetNumEvts(w, "CHG");
                double MIX_Nevts = GetNumEvts(w, "MIX");
                double UUBAR_Nevts = GetNumEvts(w, "UUBAR");
                double DDBAR_Nevts = GetNumEvts(w, "DDBAR");
                double SSBAR_Nevts = GetNumEvts(w, "SSBAR");
                double CHARM_Nevts = GetNumEvts(w, "CHARM");

                // draw
                RooPlot* x_frame = x_val->frame(Title("FBDT"));
                data->plotOn(x_frame, DataError(RooAbsData::Poisson), Cut("channelCat==0"), DrawOption("ZP"), Name("data"));
                //data->plotOn(x_frame, DataError(RooAbsData::Poisson), Cut("channelCat==0"), MarkerSize(0.4), DrawOption("ZP"), Normalization(1, RooAbsReal::ScaleType::NumEvent));
                model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kRed - 6), LineWidth(0), Components("L_x_*Signal*_ShapeSys,L_x_*CHG*_ShapeSys,L_x_*MIX*_ShapeSys,L_x_*UUBAR*_ShapeSys,L_x_*DDBAR*_ShapeSys,L_x_*SSBAR*_ShapeSys,L_x_*CHARM*_ShapeSys"), Normalization(Signal_Nevts + CHG_Nevts + MIX_Nevts + UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("signal"));
                model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kBlue - 6), LineWidth(0), Components("L_x_*CHG*_ShapeSys,L_x_*MIX*_ShapeSys,L_x_*UUBAR*_ShapeSys,L_x_*DDBAR*_ShapeSys,L_x_*SSBAR*_ShapeSys,L_x_*CHARM*_ShapeSys"), Normalization(CHG_Nevts + MIX_Nevts + UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("Charged B"));
                model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kCyan - 6), LineWidth(0), Components("L_x_*MIX*_ShapeSys,L_x_*UUBAR*_ShapeSys,L_x_*DDBAR*_ShapeSys,L_x_*SSBAR*_ShapeSys,L_x_*CHARM*_ShapeSys"), Normalization(MIX_Nevts + UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("Neutral B"));
                model->plotOn(x_frame, Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kOrange - 6), LineWidth(0), Components("L_x_*UUBAR*_ShapeSys,L_x_*DDBAR*_ShapeSys,L_x_*SSBAR*_ShapeSys,L_x_*CHARM*_ShapeSys"), Normalization(UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent), Name("Continuum"));
                data->plotOn(x_frame, DataError(RooAbsData::Poisson), Cut("channelCat==0"), DrawOption("ZP"), Name("data"));

                TCanvas* canvas = new TCanvas("sPlot", "sPlot demo", 700, 700);

                // draw PDFs
                x_frame->Draw();

                // draw legend
                TLegend* leg = new TLegend(0.9, 0.9, 0.7, 0.7);
                leg->SetFillStyle(0);
                leg->SetLineWidth(0);
                leg->AddEntry("data", "Data");
                leg->AddEntry("signal", "B #rightarrow X_{s} #nu #bar{#nu}");
                leg->AddEntry("Charged B", "Charged B");
                leg->AddEntry("Neutral B", "Neutral B");
                leg->AddEntry("Continuum", "q#bar{q}");
                leg->Draw();

                canvas->SaveAs( (std::string("fit_plot_") + "mu-" + std::to_string(val) + "_" + "err-" + std::to_string(err) + ".png").c_str() );

                delete canvas;

            }
        }

    }
}

int main(int argc, char* argv[]) {

    ::ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit"); // default: Minuit Migrad
    ::ROOT::Math::MinimizerOptions::SetDefaultStrategy(1); // default 1

    RooMsgService::instance().setStreamStatus(1, false);
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

    RooRandom::randomGenerator()->SetSeed(rd());

    ReadPIDuncorrsysFile("./KID_cov_remain_truncated.txt", "./PID_cov_remain_truncated.txt");
    ReadBRuncorrsysFile("./BR_cov_remain_truncated.txt");
    Readpi0uncorrsysFile("./pi0_cov_remain_truncated.txt");
    ReadFEIuncorrsysFile("./FEI_cov_remain_truncated.txt");

    // argv[1]: {ToyMC|LinearityTest}
    // argv[2]: injected mu when Linearity test
    // argv[3]: eps for minimizer
    // argv[4]: indicator
    // argv[5]: Num of sample
    // argv[6]: fixed nuisance parameter
    double injected_mu = -1;
    double eps = -1.0;
    int indicator = 0;
    std::string fixed_param;

    std::vector<std::string> param_names;
    
    if (std::string(argv[1]) == std::string("ToyMC")) {  // main ToyMC
        if (argc == 5) {
            injected_mu = -1;
            eps = std::atof(argv[2]);
            indicator = std::atoi(argv[3]);
            Toy_iter_num = std::atoi(argv[4]);
        }
        else {
            printf("Toy MC requires 4 arguments {eps} {indicator} {Num of sample}\n");
            exit(1);
        }
    }
    else if (std::string(argv[1]) == std::string("LinearityTest")) { // main LinearityTest 12 1234 0.8
        if (argc == 6) {
            injected_mu = std::atof(argv[2]);
            eps = std::atof(argv[3]);
            indicator = std::atoi(argv[4]);
            LT_iter_num = std::atoi(argv[5]);
        }
        else {
            printf("Linearity test requires 5 arguments: {injected mu} {eps} {indicator} {Num of sample}\n");
            exit(1);
        }
    }
    else if (std::string(argv[1]) == std::string("nuisance")) { // check the power of nuisance parameters
        if (argc == 6) {
            injected_mu = -1;
            eps = std::atof(argv[2]);
            indicator = std::atoi(argv[3]);
            Toy_iter_num = std::atoi(argv[4]);
            fixed_param = std::string(argv[5]);
        }
        else {
            printf("nuisance Toy MC requires 5 arguments {eps} {indicator} {Num of sample} {fixed param}\n");
            exit(1);
        }
    }
    else {
        printf("first argument should be {ToyMC|LinearityTest|nuisance}\n");
        exit(1);
    }

    const char* fname = "./PDFandDATA_workspace.root";

    TFile* f = TFile::Open(fname);

    RooWorkspace* w = (RooWorkspace*)f->Get("combined");

    if (std::string(argv[1]) == std::string("ToyMC")) {
        ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
        RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

        w->loadSnapshot("NominalParamValues");

        RooRealVar* x_val = w->var("obs_x_channel");
        std::unique_ptr<RooArgSet> params{model->getParameters(*x_val)};

        // save snapshot
        w->saveSnapshot("ParamValues", *params, true);
    }
    else if (std::string(argv[1]) == std::string("LinearityTest")) {
        ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
        RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

        w->loadSnapshot("NominalParamValues");

        RooRealVar* x_val = w->var("obs_x_channel");
        std::unique_ptr<RooArgSet> params{model->getParameters(*x_val)};

        // save snapshot
        w->saveSnapshot("ParamValues", *params, true);
    }
    else if (std::string(argv[1]) == std::string("nuisance")) {
        OPTIONS* options = (OPTIONS*)malloc(sizeof(OPTIONS));
        Initialize_options(options, fixed_param.c_str());
        FixParameters(w, options);
    }

    GetNameOfParams(w, &param_names);

    if (std::string(argv[1]) == std::string("ToyMC")) {
        filesaver.OpenFile(true, &param_names, injected_mu, indicator);
        MyToyMCStudy(w, &param_names, eps, indicator);
    }
    else if (std::string(argv[1]) == std::string("LinearityTest")) {
        filesaver.OpenFile(false, &param_names, injected_mu, indicator);
        MyLinearityTest(w, &param_names, injected_mu, eps, indicator);
    }
    else if (std::string(argv[1]) == std::string("nuisance")) {
        filesaver.OpenFile(true, &param_names, injected_mu, indicator);
        MyToyMCStudy(w, &param_names, eps, indicator);
    }

    f->Close();

    filesaver.CloseFile();

    return 0;
}
