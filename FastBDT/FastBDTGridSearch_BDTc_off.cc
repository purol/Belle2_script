
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
# include <vector>
#include <fstream>

#include <TMath.h>
#include <TColor.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TString.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THStack.h>
#include <TPaveText.h>
#include <TKey.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <numeric>

#include "Classifier.h"

# define N_Needed_info 37
//# define N_event_info 15
# define N_Upsilon_info 57
# define N_Bsig_info 81
# define N_Btag_info 9
# define N_decay 38 // five decay mode + others
# define N_decay_nparticles 3 // # of nu_e, B, B0
# define N_decay_syst_ff 7 // helicity angle + q2

# define Nvar 23
# define DvetoNvar 4

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

// SIGNAL MC sample number befor skimming
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

class Corrector {
private:
    const int STEP;

    const double mKp;
    const double mK0;

    const double mBp;
    const double mB0;

    const double alpha0_old;
    const double alpha1_old;
    const double alpha2_old;
    const double mp_Bp_old;
    const double mp_B0_old;

    const double alpha0_new; // 0.2545 +- 0.0090
    const double alpha1_new; // -0.71 +- 0.14
    const double alpha2_new; // 0.32 +- 0.59
    const double mp_B_new; //  5.4158 +- 0.0015

    double Total_Bp_old;
    double Total_B0_old;
    double Total_Bp_new;
    double Total_B0_new;
public:
    Corrector();
    double GetCorrectionFactor(double q2, std::string type);
    void Drawq2Plot(bool IsItNew, std::string type);
};

Corrector corrector;

Corrector::Corrector() :
    STEP(100),
    mKp(0.493677),
    mK0(0.497611),
    mBp(5.27934),
    mB0(5.27965),
    alpha0_old(0.432),
    alpha1_old(-0.664),
    alpha2_old(-1.2),
    mp_Bp_old(5.27934 + 0.046),
    mp_B0_old(5.27965 + 0.046),
    alpha0_new(0.2545),
    alpha1_new(-0.71),
    alpha2_new(0.32),
    mp_B_new(5.4158)
{
    Total_Bp_old = 0;
    Total_B0_old = 0;
    Total_Bp_new = 0;
    Total_B0_new = 0;

    // Bp old
    for (int i = 0; i < STEP; i++) {
        double q2 = ((double)i) * (mBp - mKp) * (mBp - mKp) / ((double)STEP);
        double delq2 = (mBp - mKp) * (mBp - mKp) / ((double)STEP);

        double tp = (mBp + mKp) * (mBp + mKp);
        double tm = (mBp - mKp) * (mBp - mKp);
        double t0 = tp * (1 - std::sqrt(1 - tm / tp));
        double z = (std::sqrt(tp - q2) - std::sqrt(tp - t0)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0));
        double fp = (1 / (1 - q2 / (mp_Bp_old * mp_Bp_old))) * (alpha0_old + alpha1_old * z + alpha2_old * z * z + (-alpha1_old + 2 * alpha2_old) * z * z * z / 3);
        double lambda = std::pow(mBp, 4) + std::pow(mKp, 4) + std::pow(q2, 2) - 2 * (std::pow(mBp, 2) * std::pow(mKp, 2) + std::pow(mKp, 2) * q2 + std::pow(mBp, 2) * q2);

        Total_Bp_old = Total_Bp_old + std::pow(lambda, 1.5) * std::pow(fp, 2) * delq2;
    }

    // B0 old
    for (int i = 0; i < STEP; i++) {
        double q2 = ((double)i) * (mB0 - mK0) * (mB0 - mK0) / ((double)STEP);
        double delq2 = (mB0 - mK0) * (mB0 - mK0) / ((double)STEP);

        double tp = (mB0 + mK0) * (mB0 + mK0);
        double tm = (mB0 - mK0) * (mB0 - mK0);
        double t0 = tp * (1 - std::sqrt(1 - tm / tp));
        double z = (std::sqrt(tp - q2) - std::sqrt(tp - t0)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0));
        double fp = (1 / (1 - q2 / (mp_B0_old * mp_B0_old))) * (alpha0_old + alpha1_old * z + alpha2_old * z * z + (-alpha1_old + 2 * alpha2_old) * z * z * z / 3);
        double lambda = std::pow(mB0, 4) + std::pow(mK0, 4) + std::pow(q2, 2) - 2 * (std::pow(mB0, 2) * std::pow(mK0, 2) + std::pow(mK0, 2) * q2 + std::pow(mB0, 2) * q2);

        Total_B0_old = Total_B0_old + std::pow(lambda, 1.5) * std::pow(fp, 2) * delq2;
    }

    // Bp new
    for (int i = 0; i < STEP; i++) {
        double q2 = ((double)i) * (mBp - mKp) * (mBp - mKp) / ((double)STEP);
        double delq2 = (mBp - mKp) * (mBp - mKp) / ((double)STEP);

        double tp = (mBp + mKp) * (mBp + mKp);
        double tm = (mBp - mKp) * (mBp - mKp);
        double t0 = 0.0;
        double z = (std::sqrt(tp - q2) - std::sqrt(tp - t0)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0));
        double fp = (1 / (1 - q2 / (mp_B_new * mp_B_new))) * (alpha0_new + alpha1_new * z + alpha2_new * z * z + (-alpha1_new + 2 * alpha2_new) * z * z * z / 3);
        double lambda = std::pow(mBp, 4) + std::pow(mKp, 4) + std::pow(q2, 2) - 2 * (std::pow(mBp, 2) * std::pow(mKp, 2) + std::pow(mKp, 2) * q2 + std::pow(mBp, 2) * q2);

        Total_Bp_new = Total_Bp_new + std::pow(lambda, 1.5) * std::pow(fp, 2) * delq2;
    }

    // B0 new
    for (int i = 0; i < STEP; i++) {
        double q2 = ((double)i) * (mB0 - mK0) * (mB0 - mK0) / ((double)STEP);
        double delq2 = (mB0 - mK0) * (mB0 - mK0) / ((double)STEP);

        double tp = (mB0 + mK0) * (mB0 + mK0);
        double tm = (mB0 - mK0) * (mB0 - mK0);
        double t0 = 0.0;
        double z = (std::sqrt(tp - q2) - std::sqrt(tp - t0)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0));
        double fp = (1 / (1 - q2 / (mp_B_new * mp_B_new))) * (alpha0_new + alpha1_new * z + alpha2_new * z * z + (-alpha1_new + 2 * alpha2_new) * z * z * z / 3);
        double lambda = std::pow(mB0, 4) + std::pow(mK0, 4) + std::pow(q2, 2) - 2 * (std::pow(mB0, 2) * std::pow(mK0, 2) + std::pow(mK0, 2) * q2 + std::pow(mB0, 2) * q2);

        Total_B0_new = Total_B0_new + std::pow(lambda, 1.5) * std::pow(fp, 2) * delq2;
    }
}

double Corrector::GetCorrectionFactor(double q2, std::string type) {
    if (type == "Bplus") {
        double tp = (mBp + mKp) * (mBp + mKp);
        double tm = (mBp - mKp) * (mBp - mKp);
        double t0_old = tp * (1 - std::sqrt(1 - tm / tp));
        double t0_new = 0.0;
        double z_old = (std::sqrt(tp - q2) - std::sqrt(tp - t0_old)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0_old));
        double z_new = (std::sqrt(tp - q2) - std::sqrt(tp - t0_new)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0_new));

        double fp_old = (1 / (1 - q2 / (mp_Bp_old * mp_Bp_old))) * (alpha0_old + alpha1_old * z_old + alpha2_old * z_old * z_old + (-alpha1_old + 2 * alpha2_old) * z_old * z_old * z_old / 3);
        double fp_new = (1 / (1 - q2 / (mp_B_new * mp_B_new))) * (alpha0_new + alpha1_new * z_new + alpha2_new * z_new * z_new + (-alpha1_new + 2 * alpha2_new) * z_new * z_new * z_new / 3);
        double lambda = std::pow(mBp, 4) + std::pow(mKp, 4) + std::pow(q2, 2) - 2 * (std::pow(mBp, 2) * std::pow(mKp, 2) + std::pow(mKp, 2) * q2 + std::pow(mBp, 2) * q2);

        if ((std::pow(lambda, 1.5) * std::pow(fp_old, 2) / Total_Bp_old) == 0) {
            printf("[ERROR] try to divide by 0!\n");
            exit(1);
        }
        return (std::pow(lambda, 1.5) * std::pow(fp_new, 2) / Total_Bp_new) / (std::pow(lambda, 1.5) * std::pow(fp_old, 2) / Total_Bp_old);
    }
    else if (type == "Bzero") {
        double tp = (mB0 + mK0) * (mB0 + mK0);
        double tm = (mB0 - mK0) * (mB0 - mK0);
        double t0_old = tp * (1 - std::sqrt(1 - tm / tp));
        double t0_new = 0.0;
        double z_old = (std::sqrt(tp - q2) - std::sqrt(tp - t0_old)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0_old));
        double z_new = (std::sqrt(tp - q2) - std::sqrt(tp - t0_new)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0_new));

        double fp_old = (1 / (1 - q2 / (mp_B0_old * mp_B0_old))) * (alpha0_old + alpha1_old * z_old + alpha2_old * z_old * z_old + (-alpha1_old + 2 * alpha2_old) * z_old * z_old * z_old / 3);
        double fp_new = (1 / (1 - q2 / (mp_B_new * mp_B_new))) * (alpha0_new + alpha1_new * z_new + alpha2_new * z_new * z_new + (-alpha1_new + 2 * alpha2_new) * z_new * z_new * z_new / 3);
        double lambda = std::pow(mB0, 4) + std::pow(mK0, 4) + std::pow(q2, 2) - 2 * (std::pow(mB0, 2) * std::pow(mK0, 2) + std::pow(mK0, 2) * q2 + std::pow(mB0, 2) * q2);

        if ((std::pow(lambda, 1.5) * std::pow(fp_old, 2) / Total_B0_old) == 0) {
            printf("[ERROR] try to divide by 0!\n");
            exit(1);
        }
        return (std::pow(lambda, 1.5) * std::pow(fp_new, 2) / Total_B0_new) / (std::pow(lambda, 1.5) * std::pow(fp_old, 2) / Total_B0_old);
    }

    printf("[ERROR] undefined type name");
    exit(1);

    return 0;
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
    double FluctuateCorrection(int Decay[N_decay], double MXs, SystType systtype, int TargetCategory, bool IsTargetCategoryUp, Sample sample, std::string type);
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

double Corrector_Fragmentation::FluctuateCorrection(int Decay[N_decay], double MXs, SystType systtype, int TargetCategory, bool IsTargetCategoryUp, Sample sample, std::string type) {
    /*
    Nevt of TargetCategory is fluctuated.
    we calculate the change of Decay[N_decay] (Category)
    */
    int Category = Classify(Decay, sample);
    int MxsBin = GetMxBin(MXs, sample);

    const double RelativeUncertainty = 1.0;
    double TotalNevtAtMxsBinWithMissing = Total_Nevt_Nominal_before_Xsgamma_MC15[MxsBin] + Nevt_Nominal_missing_before_Xsgamma_MC15[MxsBin];
    double TargetNevtAtMxsBin;
    if (TargetCategory == N_Category_gamma) TargetNevtAtMxsBin = Nevt_Nominal_missing_before_Xsgamma_MC15[MxsBin]; // if it is missing mode
    else TargetNevtAtMxsBin = Total_Nevt_Nominal_before_Xsgamma_MC15[MxsBin] * Fragmentation_Xsgamma[MxsBin][TargetCategory];
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

using std::string;

std::vector<std::string> var_names;

void load_files(const char* dirname, std::vector<string>* names) {
    TSystemDirectory dir(dirname, dirname);
    TList* files = dir.GetListOfFiles();
    if (files) {
        TSystemFile* file;
        TString fname;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(".root")) {
                names->push_back(fname.Data());
            }
        }
    }
}

void load_files(const char* dirname, std::vector<string>* names, const char* included_string) {
    TSystemDirectory dir(dirname, dirname);
    TList* files = dir.GetListOfFiles();
    if (files) {
        TSystemFile* file;
        TString fname;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(".root") && fname.Contains(included_string)) {
                names->push_back(fname.Data());
            }
        }
    }
}

void FillVariables(const char * filename, std::vector<float> input_vars[Nvar], std::vector<bool>* IsSignal, std::vector<float>* weight, bool tempissignal, float weight_N = 1.0) {
    TFile* input_file = new TFile(filename, "read");

    TTree* tree_data = (TTree*)input_file->Get("data");

    double Vars[Nvar];
    int flag;

    double Dc_chiProb; // 0.0
    double Dc_pvalue_med;
    double Dc_pvalue_std; // 0.0
    double Dc_dr; // -1.0
    double Dc_dz; // -100.0
    double Dc_M; // 0.0
    double D0_chiProb;
    double D0_pvalue_med;
    double D0_pvalue_std;
    double D0_dr;
    double D0_dz;
    double D0_M;

    int Decay_Kplus = -1;
    int Decay_Kplusstar_ch1 = -1;
    int Decay_Kplusstar_ch2 = -1;
    int Decay_K0 = -1;
    int Decay_K0star_ch1 = -1;
    int Decay_K0star_ch2 = -1;

    double Mxs = -1;
    double Pcms = -1;

    tree_data->SetBranchAddress("Bsig_KSFWVariables_hso04", &Vars[0]);
    tree_data->SetBranchAddress("Bsig_cosTBTO", &Vars[1]);
    tree_data->SetBranchAddress("Bsig_useCMSFrame_p", &Vars[2]);
    tree_data->SetBranchAddress("Btag_CleoConeCS_1", &Vars[3]);
    tree_data->SetBranchAddress("Btag_CleoConeCS_2", &Vars[4]);
    tree_data->SetBranchAddress("Btag_CleoConeCS_3", &Vars[5]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hoo2", &Vars[6]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hoo4", &Vars[7]);
    tree_data->SetBranchAddress("Btag_cosTBTO", &Vars[8]);
    tree_data->SetBranchAddress("Btag_extraInfo_SignalProbability", &Vars[9]);
    tree_data->SetBranchAddress("Btag_useCMSFrame_theta", &Vars[10]);
    tree_data->SetBranchAddress("extraInfo__boEeclv133__bc", &Vars[11]);
    tree_data->SetBranchAddress("extraInfo__boNgammav133__bc", &Vars[12]);
    tree_data->SetBranchAddress("foxWolframR3", &Vars[13]);
    tree_data->SetBranchAddress("harmonicMomentThrust2", &Vars[14]);
    tree_data->SetBranchAddress("missingMomentumOfEvent", &Vars[15]);
    tree_data->SetBranchAddress("missingMomentumOfEvent_theta", &Vars[16]);
    tree_data->SetBranchAddress("nRemainingTracksInEvent", &Vars[17]);
    tree_data->SetBranchAddress("roePTheta__bocleanMask__bc", &Vars[18]);

    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dc_pValue_med", &Dc_pvalue_med);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dc_pValue_std", &Dc_pvalue_std);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_chiProb", &Dc_chiProb);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_dz", &Dc_dz);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_M", &Dc_M);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_D0_pValue_med", &D0_pvalue_med);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_D0_pValue_std", &D0_pvalue_std);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_D0simpleveto_chiProb", &D0_chiProb);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_D0simpleveto_dz", &D0_dz);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_D0simpleveto_M", &D0_M);

    tree_data->SetBranchAddress("flag", &flag);

    tree_data->SetBranchAddress("Bsig_M", &Mxs);

    int Nevt = 0;
    //printf("%lld entries...\n", tree_data->GetEntries());
    for (unsigned int j = 0; j < tree_data->GetEntries(); j++) { // Fill
        tree_data->GetEntry(j);

        Nevt++;

        for (unsigned int k = 0; k < Nvar - DvetoNvar; k++) input_vars[k].push_back((float) Vars[k]); 

        if (Dc_chiProb > -0.5) {
            input_vars[Nvar - DvetoNvar + 0].push_back((float)Dc_pvalue_std);
            input_vars[Nvar - DvetoNvar + 1].push_back((float)Dc_M);
        }
        else {
            input_vars[Nvar - DvetoNvar + 0].push_back((float)0.0);
            input_vars[Nvar - DvetoNvar + 1].push_back((float)0.0);
        }
        if (D0_chiProb > -0.5) {
            input_vars[Nvar - DvetoNvar + 2].push_back((float)D0_pvalue_std);
            input_vars[Nvar - DvetoNvar + 3].push_back((float)D0_M);
        }
        else {
            input_vars[Nvar - DvetoNvar + 2].push_back((float)0.0);
            input_vars[Nvar - DvetoNvar + 3].push_back((float)0.0);
        }

        IsSignal->push_back(tempissignal);

        weight->push_back(weight_N);

    }

    input_file->Close();
    //printf("==> Total %d events survive...\n", Nevt);
}

double PrintMaximumFOM(const FastBDT::Classifier& classifier, std::vector<std::vector<float>> InputVariables, std::vector<bool> IsSignal, std::vector<float> weight) {
    const int step = 100;
    double FOM_max = -1;

    for (int i = 0; i < step; i++) {
        float value = ((float)i) / ((float)step);
        double NBKG = 0;
        double NSIG = 0;

        for (unsigned int i = 0; i < IsSignal.size(); ++i) {
            std::vector<float> temp;
            for (int j = 0; j < Nvar; j++) temp.push_back(InputVariables.at(j).at(i));
            float p = classifier.predict(temp);
            if (p > value) {
                if(IsSignal[i]) NSIG = NSIG + weight[i];
                else NBKG = NBKG + weight[i];
            }
        }

        double FOM = NSIG / std::sqrt(NBKG + NSIG);
        if (FOM > FOM_max) FOM_max = FOM;

    }

    return FOM_max;
}

double PrintAUC(const FastBDT::Classifier& classifier, std::vector<std::vector<float>> InputVariables, std::vector<bool> IsSignal, std::vector<float> weight) {
    const int step = 100;
    double AUC = 0;
    double NBKG_total = 0;
    double NSIG_total = 0;
    std::vector<double> TPRs;
    std::vector<double> FPRs;

    for (unsigned int i = 0; i < IsSignal.size(); ++i) {
        if (IsSignal[i]) NSIG_total = NSIG_total + weight[i];
        else NBKG_total = NBKG_total + weight[i];
    }

    for (int i = 0; i < step; i++) {
        float value = ((float)i) / ((float)step);
        double NBKG = 0;
        double NSIG = 0;

        for (unsigned int i = 0; i < IsSignal.size(); ++i) {
            std::vector<float> temp;
            for (int j = 0; j < Nvar; j++) temp.push_back(InputVariables.at(j).at(i));
            float p = classifier.predict(temp);
            if (p >= value) {
                if (IsSignal[i]) NSIG = NSIG + weight[i];
                else NBKG = NBKG + weight[i];
            }
        }

        double TPR = NSIG / NSIG_total;
        double FPR = NBKG / NBKG_total;

        TPRs.push_back(TPR);
        FPRs.push_back(FPR);
    }

    for (unsigned int i = 0; i < TPRs.size(); ++i) {
        if ( i != TPRs.size() - 1) {
            double del_FPR = FPRs.at(i) - FPRs.at(i + 1);
            double avg_TPR = (TPRs.at(i) + TPRs.at(i + 1)) / 2.0;
            AUC = AUC + del_FPR * avg_TPR;
        }
        else {
            double del_FPR = FPRs.at(i) - 0.0;
            double avg_TPR = (TPRs.at(i) + 0.0) / 2.0;
            AUC = AUC + del_FPR * avg_TPR;
        }
    }

    return AUC;
}

double PrintAVG(const FastBDT::Classifier& classifier, std::vector<std::vector<float>> InputVariables, std::vector<bool> IsSignal, std::vector<float> weight, bool SelectSignal) {
    double NBKG_total = 0;
    double NSIG_total = 0;

    double NBKG_AVG = 0;
    double NSIG_AVG = 0;

    for (unsigned int i = 0; i < IsSignal.size(); ++i) {
        if (IsSignal[i]) NSIG_total = NSIG_total + weight[i];
        else NBKG_total = NBKG_total + weight[i];
    }

    for (unsigned int i = 0; i < IsSignal.size(); ++i) {
        std::vector<float> temp;
        for (int j = 0; j < Nvar; j++) temp.push_back(InputVariables.at(j).at(i));
        float p = classifier.predict(temp);

        if (IsSignal[i]) NSIG_AVG = NSIG_AVG + p * weight[i];
        else NBKG_AVG = NBKG_AVG + p * weight[i];

    }

    NSIG_AVG = NSIG_AVG / NSIG_total;
    NBKG_AVG = NBKG_AVG / NBKG_total;

    if (SelectSignal) return NSIG_AVG;
    else return NBKG_AVG;
}

int main(int argc, char* argv[]) // offres total: 42.329/fb
{
    // grid search
    // unsigned int nTrees[5] = { 100, 500, 1000, 1500, 2000 };  default is 100
    // unsigned int depth[3] = { 2, 3, 4 };  default is 3 
    // double shrinkage[4] = { 0.05, 0.1, 0.15, 0.2 };  default is 0.1
    // double subsample[5] = { 0.3, 0.4, 0.5, 0.6, 0.7 };  default is 0.5
    // unsigned int binning[4] = { 6, 7, 8, 9 };  default is 2^8 bins per feature

    unsigned int nTrees = (unsigned int)atoi(argv[1]);
    unsigned int depth = (unsigned int)atoi(argv[2]);
    double shrinkage = atof(argv[3]);
    double subsample = atof(argv[4]);
    unsigned int binning_num = (unsigned int)atoi(argv[5]);

    // set classifier option
    FastBDT::Classifier classifier;
    classifier.SetNTrees(nTrees);
    classifier.SetDepth(depth);
    classifier.SetShrinkage(shrinkage);
    classifier.SetSubsample(subsample);
    std::vector<unsigned int> binning(Nvar, binning_num); classifier.SetBinning(binning);



    // input file
    const char* off_data = "/home/belle2/junewoo/storage_b1/GridSearch_BDTc/data/SIGNAL_analysis/validation_v000/final_output/DataFile";
   
    const char* off_MC_UUBAR_train = "/home/belle2/junewoo/storage_b1/GridSearch_BDTc/MC/UUBAR_analysis/train_v000/final_output/DataFile";
    const char* off_MC_DDBAR_train = "/home/belle2/junewoo/storage_b1/GridSearch_BDTc/MC/DDBAR_analysis/train_v000/final_output/DataFile";
    const char* off_MC_SSBAR_train = "/home/belle2/junewoo/storage_b1/GridSearch_BDTc/MC/SSBAR_analysis/train_v000/final_output/DataFile";
    const char* off_MC_CHARM_train = "/home/belle2/junewoo/storage_b1/GridSearch_BDTc/MC/CHARM_analysis/train_v000/final_output/DataFile";

    const char* off_MC_UUBAR_test = "/home/belle2/junewoo/storage_b1/GridSearch_BDTc/MC/UUBAR_analysis/test_v000/final_output/DataFile";
    const char* off_MC_DDBAR_test = "/home/belle2/junewoo/storage_b1/GridSearch_BDTc/MC/DDBAR_analysis/test_v000/final_output/DataFile";
    const char* off_MC_SSBAR_test = "/home/belle2/junewoo/storage_b1/GridSearch_BDTc/MC/SSBAR_analysis/test_v000/final_output/DataFile";
    const char* off_MC_CHARM_test = "/home/belle2/junewoo/storage_b1/GridSearch_BDTc/MC/CHARM_analysis/test_v000/final_output/DataFile";

    // define input of the classifier
    std::vector<std::vector<float>> InputVariables;
    std::vector<bool> IsSignal;
    std::vector<float> weight;

    // define input variables
    std::vector<float> input_vars[Nvar];

    {
        std::vector<string> names;
        load_files(off_data, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 == 0) continue; // offresdata: 7758 of 9012 (36.1697/fb)
            FillVariables((off_data + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(off_MC_UUBAR_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((off_MC_UUBAR_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, Scale_UUBAR_train * (0.0364390 / 0.361673));
        }
    }
    {
        std::vector<string> names;
        load_files(off_MC_DDBAR_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((off_MC_DDBAR_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, Scale_DDBAR_train * (0.0364390 / 0.361673));
        }
    }
    {
        std::vector<string> names;
        load_files(off_MC_SSBAR_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((off_MC_SSBAR_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, Scale_SSBAR_train * (0.0364390 / 0.361673));
        }
    }
    {
        std::vector<string> names;
        load_files(off_MC_CHARM_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((off_MC_CHARM_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, Scale_CHARM_train * (0.0364390 / 0.361673));
        }
    }


    // convert to double vector
    for (int i = 0; i < Nvar; i++) {
        InputVariables.push_back(input_vars[i]);
    }



    // fit
    classifier.fit(InputVariables, IsSignal, weight);


    // get FOM for training sample
    // double train_FOM = PrintMaximumFOM(classifier, InputVariables3, IsSignal3, weight3);

    // get AUC for training sample
    double train_AUC = PrintAUC(classifier, InputVariables, IsSignal, weight);


    // get average value for BKG and SIG
    double train_AVG_SIG = PrintAVG(classifier, InputVariables, IsSignal, weight, true);
    double train_AVG_BKG = PrintAVG(classifier, InputVariables, IsSignal, weight, false);


    // clear vector to save memory
    for (unsigned int i = 0; i < InputVariables.size(); ++i) std::vector<float>().swap(InputVariables.at(i));
    std::vector<std::vector<float>>().swap(InputVariables);
    std::vector<bool>().swap(IsSignal);
    std::vector<float>().swap(weight);
    for (int i = 0; i < Nvar; i++) std::vector<float>().swap(input_vars[i]);



    // test sample
    std::vector<std::vector<float>> InputVariables2;
    std::vector<bool> IsSignal2;
    std::vector<float> weight2;

    std::vector<float> input_vars2[Nvar];

    {
        std::vector<string> names;
        load_files(off_data, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 != 0) continue; // offresdata: 1254 of 9012  (6.1593/fb)
            FillVariables((off_data + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(off_MC_UUBAR_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((off_MC_UUBAR_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, Scale_UUBAR_test * (0.0058900 / 0.361673));
        }
    }
    {
        std::vector<string> names;
        load_files(off_MC_DDBAR_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((off_MC_DDBAR_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, Scale_DDBAR_test * (0.0058900 / 0.361673));
        }
    }
    {
        std::vector<string> names;
        load_files(off_MC_SSBAR_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((off_MC_SSBAR_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, Scale_SSBAR_test * (0.0058900 / 0.361673));
        }
    }
    {
        std::vector<string> names;
        load_files(off_MC_CHARM_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((off_MC_CHARM_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, Scale_CHARM_test* (0.0058900 / 0.361673));
        }
    }


    // convert to double vector
    for (int i = 0; i < Nvar; i++) {
        InputVariables2.push_back(input_vars2[i]);
    }



    // get FOM for testing sample
    // double test_FOM = PrintMaximumFOM(classifier, InputVariables2, IsSignal2, weight2);

    // get AUC for testing sample
    double test_AUC = PrintAUC(classifier, InputVariables2, IsSignal2, weight2);

    // get average value for BKG and SIG
    double AVG_SIG = PrintAVG(classifier, InputVariables2, IsSignal2, weight2, true);
    double AVG_BKG = PrintAVG(classifier, InputVariables2, IsSignal2, weight2, false);

    // clear vector to save memory
    for (unsigned int i = 0; i < InputVariables2.size(); ++i) std::vector<float>().swap(InputVariables2.at(i));
    std::vector<std::vector<float>>().swap(InputVariables2);
    std::vector<bool>().swap(IsSignal2);
    std::vector<float>().swap(weight2);
    for (int i = 0; i < Nvar; i++) std::vector<float>().swap(input_vars2[i]);



    printf("%u_%u_%lf_%lf_%u %lf %lf %lf %lf %lf %lf\n", nTrees, depth, shrinkage, subsample, binning_num, train_AUC, test_AUC, train_AVG_SIG, train_AVG_BKG, AVG_SIG, AVG_BKG);

    FILE* fp;
    fp = fopen(("/home/belle2/junewoo/storage_b1/GridSearch_BDTc/out/Result_" + std::string(argv[1]) + "_" + std::string(argv[2]) + "_" + std::string(argv[3]) + "_" + std::string(argv[4]) + "_" + std::string(argv[5])).c_str(), "w");
    fprintf(fp, "%u_%u_%lf_%lf_%u %lf %lf\n", nTrees, depth, shrinkage, subsample, binning_num, train_AUC, test_AUC);
    fclose(fp);



    // save model
    std::fstream out_stream(("/home/belle2/junewoo/storage_b1/GridSearch_BDTc/out/classifier_" + std::string(argv[1]) + "_" + std::string(argv[2]) + "_" + std::string(argv[3]) + "_" + std::string(argv[4]) + "_" + std::string(argv[5])+".weightfile").c_str(), std::ios_base::out | std::ios_base::trunc);
    out_stream << classifier << std::endl;
    out_stream.close();

    return 0;
}
