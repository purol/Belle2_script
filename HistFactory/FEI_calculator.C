#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include <cmath>
#include <float.h>
#include <string>
#include <string.h>
#include "TFile.h"
#include "TTree.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "RooGenericPdf.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "TH1.h"
#include "TF1.h"
#include "RooDataHist.h"
#include <time.h>
#include <random>
#include <algorithm>
using namespace RooFit;
using std::string;
using std::to_string;

# define MyEPSILON 0.000001

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

// scale facto for systematic MC sample
# define N_K0star_nunubar_syst 10000000.0
# define N_Kplusstar_nunubar_syst 10000000.0
# define N_Xsu_nonresonant_nunubar_syst 10000000.0
# define N_Xsd_nonresonant_nunubar_syst 10000000.0

// scale factor for each systematic MC sample (364.436 - 2.763 = 361.673/fb)
# define Scale_Kplusstar_syst (N_Kplusstar_nunubar_LS1/N_Kplusstar_nunubar_syst)
# define Scale_Xsu_nonresonant_syst (N_Xsu_nonresonant_nunubar_LS1/N_Xsu_nonresonant_nunubar_syst)
# define Scale_K0star_syst (N_K0star_nunubar_LS1/N_K0star_nunubar_syst)
# define Scale_Xsd_nonresonant_syst (N_Xsd_nunubar_LS1/N_Xsd_nonresonant_nunubar_syst)

# define KS0_rel_uncertainty 0.6 // %/cm
# define track_rel_uncertainty 0.69 // %
// # define pi0_correction 0.932
// # define pi0_rel_uncertainty ((0.0369 / 0.932) * 100.0) // %
# define Kaon_PID_max_uncertainty 0.1 // not percentage. relative uncertainty

# define FEI_cal_Bc_num 12
# define FEI_cal_B0_num 11
double FEI_cal_Bc[FEI_cal_Bc_num] = { 1.04, 0.79, 0.69, 0.56, 0.97, 0.95, 0.74, 0.57, 0.91, 0.51, 0.34, 0.59 };
double FEI_cal_Bc_uncertainty[FEI_cal_Bc_num] = { 0.03, 0.03, 0.05, 0.11, 0.03, 0.03, 0.02, 0.06, 0.1, 0.13, 0.07, 0.02 }; // not relative uncertainty. absolute uncertainty
double FEI_cal_Bc_modeID[FEI_cal_Bc_num] = { 0.0, 1.0, 3.0, 4.0, 15.0, 16.0, 18.0, 19.0, 23.0, 24.0, 30.0, -1.0 };
double FEI_cal_B0[FEI_cal_B0_num] = { 1.16, 0.95, 0.84, 0.78, 0.99, 1.01, 0.67, 0.65, 0.69, 0.58, 0.81 };
double FEI_cal_B0_uncertainty[FEI_cal_B0_num] = { 0.04, 0.03, 0.02, 0.02, 0.03, 0.03, 0.02, 0.02, 0.02, 0.16, 0.13 }; // not relative uncertainty. absolute uncertainty
double FEI_cal_B0_modeID[FEI_cal_B0_num] = { 0.0, 1.0, 3.0, 4.0, 5.0, 15.0, 16.0, 18.0, 19.0, 26.0, -1.0 };

// FEI_cal_Bc_eigenvalue[0] <---------> FEI_cal_Bc_eigenvector[0][i]
double FEI_cal_Bc_eigenvalue[FEI_cal_Bc_num] = { 0.0 };
double FEI_cal_Bc_eigenvector[FEI_cal_Bc_num][FEI_cal_Bc_num] = { 0.0 };
double FEI_cal_Bc_fluctuated[FEI_cal_Bc_num] = { 0.0 };
double FEI_cal_B0_eigenvalue[FEI_cal_B0_num] = { 0.0 };
double FEI_cal_B0_eigenvector[FEI_cal_B0_num][FEI_cal_B0_num] = { 0.0 };
double FEI_cal_B0_fluctuated[FEI_cal_B0_num] = { 0.0 };

# define Xsu_frag_decay1 531.4
# define Xsu_frag_decay2 1062.1
# define Xsu_frag_decay3 673.9
# define Xsu_frag_decay4 825.3
# define Xsu_frag_decay5 473.3
# define Xsu_frag_decay6 119.1
# define Xsu_frag_decay7 24.6
# define Xsu_frag_decay8 148.2
# define Xsu_frag_decay9 127.6
# define Xsu_frag_decay10 102.9
# define Xsu_frag_decay11 91.8
# define Xsu_frag_decay12 54.8
# define Xsu_frag_decay13 21.7
# define Xsu_frag_decay14 11.0
# define Xsu_frag_decay30 858.4

# define Xsd_frag_decay16 1072.6
# define Xsd_frag_decay17 536.2
# define Xsd_frag_decay18 835.5
# define Xsd_frag_decay19 677.4
# define Xsd_frag_decay20 120.5
# define Xsd_frag_decay21 475.9
# define Xsd_frag_decay22 149.8
# define Xsd_frag_decay23 25.0
# define Xsd_frag_decay24 129.2
# define Xsd_frag_decay25 104.4
# define Xsd_frag_decay26 92.4
# define Xsd_frag_decay27 54.9
# define Xsd_frag_decay28 14.5
# define Xsd_frag_decay29 10.9
# define Xsd_frag_decay30 871.0

std::random_device rd;
std::default_random_engine generator(rd());

# define N_PID_syst 73
# define N_fakeE_syst 37
# define N_fakeMU_syst 49
double PID_correction[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
double PID_correction_stat_uncer[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
double PID_correction_sys_uncer[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
double PID_correction_uncer[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss

double PID_fakeE_correction[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
double PID_fakeE_uncer[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+

double PID_fakeMU_correction[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+
double PID_fakeMU_uncer[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

# define N_pi0_syst 8
const double pi0_correction[N_pi0_syst] = {
    0.917, 0.965, 0.988, 1.013, 1.042, 1.044, 1.011, 1.0
};
const double pi0_stat_uncer[N_pi0_syst] = {
    0.004, 0.004, 0.004, 0.005, 0.004, 0.005, 0.005, 0.0
};
const double pi0_sys_uncer1[N_pi0_syst] = {
    0.049, 0.036, 0.079, 0.058, 0.045, 0.041, 0.040, 0.0
};
const double pi0_sys_uncer2[N_pi0_syst] = {
    0.0, 0.0, 0.0, 0.0, 0.039, 0.051, 0.030, 0.0
};

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

class Corrector_Knn {
private:

    // K+nn
    int STEP_Knn;
    double mininvM_Knn;
    double maxinvM_Knn;
    TH1D* weights_Knn;
    const double DECAY_DEC_BR_Knn;
    const double new_BR_Knn;
    const double Nraw_initial_Knn;
    double Nscale_initial_Knn;

    // K*+nn
    int STEP_Kstarnn;
    double mininvM_Kstarnn;
    double maxinvM_Kstarnn;
    TH1D* weights_Kstarnn;
    const double DECAY_DEC_BR_Kstarnn;
    const double new_BR_Kstarnn;
    const double Nraw_initial_Kstarnn;
    double Nscale_initial_Kstarnn;

    // K0nn
    int STEP_K0nn;
    double mininvM_K0nn;
    double maxinvM_K0nn;
    TH1D* weights_K0nn;
    const double DECAY_DEC_BR_K0nn;
    const double new_BR_K0nn;
    const double Nraw_initial_K0nn;
    double Nscale_initial_K0nn;

    // K0*nn
    int STEP_K0starnn;
    double mininvM_K0starnn;
    double maxinvM_K0starnn;
    TH1D* weights_K0starnn;
    const double DECAY_DEC_BR_K0starnn;
    const double new_BR_K0starnn;
    const double Nraw_initial_K0starnn;
    double Nscale_initial_K0starnn;

    const double N_EPSILON;
    const double CUTOFF;

public:
    Corrector_Knn();
    double GetCorrectionFactor(double invM_Knn, double invM_Kstarnn, double invM_K0nn, double invM_K0starnn, double N_Knn, double N_Kstarnn, double N_K0nn, double N_K0starnn);
    double GetCorrectionFactorAtGeneric(double invM_Knn, double invM_Kstarnn, double invM_K0nn, double invM_K0starnn, double N_Knn, double N_Kstarnn, double N_K0nn, double N_K0starnn);
};

Corrector_Knn corrector_Knn;

Corrector_Knn::Corrector_Knn() :
    DECAY_DEC_BR_Knn(0.0000057),
    new_BR_Knn(0.0000059),
    Nraw_initial_Knn(1000000.0),
    DECAY_DEC_BR_Kstarnn(0.0000057),
    new_BR_Kstarnn(0.0000036),
    Nraw_initial_Kstarnn(1000000.0),
    DECAY_DEC_BR_K0nn(0.000002),
    new_BR_K0nn(0.00000266),
    Nraw_initial_K0nn(1000000.0),
    DECAY_DEC_BR_K0starnn(0.0000056),
    new_BR_K0starnn(0.00000124),
    Nraw_initial_K0starnn(1000000.0),
    N_EPSILON(0.01),
    CUTOFF(50.0)
{
    FILE* fp;

    // read Knn weights
    fp = fopen("/home/jwpark/storage/BKG_gbasf2/Knn_weight/Knn_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &STEP_Knn, &mininvM_Knn, &maxinvM_Knn);
    weights_Knn = new TH1D("Knn_weights", ";;", STEP_Knn, mininvM_Knn, maxinvM_Knn);
    for (int i = 0; i < STEP_Knn; i++) {
        double temp;
        fscanf(fp, "%lf\n", &temp);
        if (temp < CUTOFF) weights_Knn->SetBinContent(i + 1, temp);
        else weights_Knn->SetBinContent(i + 1, CUTOFF);
    }
    fclose(fp);

    // read Kstarnn weights
    fp = fopen("/home/jwpark/storage/BKG_gbasf2/Knn_weight/Kstarnn_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &STEP_Kstarnn, &mininvM_Kstarnn, &maxinvM_Kstarnn);
    weights_Kstarnn = new TH1D("Kstarnn_weights", ";;", STEP_Kstarnn, mininvM_Kstarnn, maxinvM_Kstarnn);
    for (int i = 0; i < STEP_Kstarnn; i++) {
        double temp;
        fscanf(fp, "%lf\n", &temp);
        if (temp < CUTOFF) weights_Kstarnn->SetBinContent(i + 1, temp);
        else weights_Kstarnn->SetBinContent(i + 1, CUTOFF);
    }
    fclose(fp);

    // read K0nn weights
    fp = fopen("/home/jwpark/storage/BKG_gbasf2/Knn_weight/K0nn_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &STEP_K0nn, &mininvM_K0nn, &maxinvM_K0nn);
    weights_K0nn = new TH1D("K0nn_weights", ";;", STEP_K0nn, mininvM_K0nn, maxinvM_K0nn);
    for (int i = 0; i < STEP_K0nn; i++) {
        double temp;
        fscanf(fp, "%lf\n", &temp);
        if (temp < CUTOFF) weights_K0nn->SetBinContent(i + 1, temp);
        else weights_K0nn->SetBinContent(i + 1, CUTOFF);
    }
    fclose(fp);

    // read K0starnn weights
    fp = fopen("/home/jwpark/storage/BKG_gbasf2/Knn_weight/K0starnn_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &STEP_K0starnn, &mininvM_K0starnn, &maxinvM_K0starnn);
    weights_K0starnn = new TH1D("K0starnn_weights", ";;", STEP_K0starnn, mininvM_K0starnn, maxinvM_K0starnn);
    for (int i = 0; i < STEP_K0starnn; i++) {
        double temp;
        fscanf(fp, "%lf\n", &temp);
        if (temp < CUTOFF) weights_K0starnn->SetBinContent(i + 1, temp);
        else weights_K0starnn->SetBinContent(i + 1, CUTOFF);
    }
    fclose(fp);

    // calculate the expected number of event
    Nscale_initial_Knn = (2.0 * N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0)) * new_BR_Knn);
    Nscale_initial_Kstarnn = (2.0 * N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0)) * new_BR_Kstarnn);
    Nscale_initial_K0nn = (2.0 * N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0)) * new_BR_K0nn);
    Nscale_initial_K0starnn = (2.0 * N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0)) * new_BR_K0starnn);
}

double Corrector_Knn::GetCorrectionFactor(double invM_Knn, double invM_Kstarnn, double invM_K0nn, double invM_K0starnn, double N_Knn, double N_Kstarnn, double N_K0nn, double N_K0starnn) {

    double Correction_Knn = 1;
    double Correction_Kstarnn = 1;
    double Correction_K0nn = 1;
    double Correction_K0starnn = 1;

    if (N_Knn < N_EPSILON) Correction_Knn = 1;
    else {
        int Bin = weights_Knn->FindBin(invM_Knn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_Knn) Bin = STEP_Knn;
        Correction_Knn = std::pow((Nscale_initial_Knn / Nraw_initial_Knn) * weights_Knn->GetBinContent(Bin), N_Knn); // BR correction * invM correction
    }

    if (N_Kstarnn < N_EPSILON) Correction_Kstarnn = 1;
    else {
        int Bin = weights_Kstarnn->FindBin(invM_Kstarnn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_Kstarnn) Bin = STEP_Kstarnn;
        Correction_Kstarnn = std::pow((Nscale_initial_Kstarnn / Nraw_initial_Kstarnn) * weights_Kstarnn->GetBinContent(Bin), N_Kstarnn); // BR correction * invM correction
    }

    if (N_K0nn < N_EPSILON) Correction_K0nn = 1;
    else {
        int Bin = weights_K0nn->FindBin(invM_K0nn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_K0nn) Bin = STEP_K0nn;
        Correction_K0nn = std::pow((Nscale_initial_K0nn / Nraw_initial_K0nn) * weights_K0nn->GetBinContent(Bin), N_K0nn); // BR correction * invM correction
    }

    if (N_K0starnn < N_EPSILON) Correction_K0starnn = 1;
    else {
        int Bin = weights_K0starnn->FindBin(invM_K0starnn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_K0starnn) Bin = STEP_K0starnn;
        Correction_K0starnn = std::pow((Nscale_initial_K0starnn / Nraw_initial_K0starnn) * weights_K0starnn->GetBinContent(Bin), N_K0starnn); // BR correction * invM correction
    }

    return Correction_Knn * Correction_Kstarnn * Correction_K0nn * Correction_K0starnn;
}

double Corrector_Knn::GetCorrectionFactorAtGeneric(double invM_Knn, double invM_Kstarnn, double invM_K0nn, double invM_K0starnn, double N_Knn, double N_Kstarnn, double N_K0nn, double N_K0starnn) {

    double Correction_Knn = 1;
    double Correction_Kstarnn = 1;
    double Correction_K0nn = 1;
    double Correction_K0starnn = 1;

    if (N_Knn < N_EPSILON) Correction_Knn = 1;
    else {
        int Bin = weights_Knn->FindBin(invM_Knn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_Knn) Bin = STEP_Knn;
        Correction_Knn = std::pow((new_BR_Knn / DECAY_DEC_BR_Knn) * weights_Knn->GetBinContent(Bin), N_Knn); // BR correction * invM correction
    }

    if (N_Kstarnn < N_EPSILON) Correction_Kstarnn = 1;
    else {
        int Bin = weights_Kstarnn->FindBin(invM_Kstarnn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_Kstarnn) Bin = STEP_Kstarnn;
        Correction_Kstarnn = std::pow((new_BR_Kstarnn / DECAY_DEC_BR_Kstarnn) * weights_Kstarnn->GetBinContent(Bin), N_Kstarnn); // BR correction * invM correction
    }

    if (N_K0nn < N_EPSILON) Correction_K0nn = 1;
    else {
        int Bin = weights_K0nn->FindBin(invM_K0nn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_K0nn) Bin = STEP_K0nn;
        Correction_K0nn = std::pow((new_BR_K0nn / DECAY_DEC_BR_K0nn) * weights_K0nn->GetBinContent(Bin), N_K0nn); // BR correction * invM correction
    }

    if (N_K0starnn < N_EPSILON) Correction_K0starnn = 1;
    else {
        int Bin = weights_K0starnn->FindBin(invM_K0starnn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_K0starnn) Bin = STEP_K0starnn;
        Correction_K0starnn = std::pow((new_BR_K0starnn / DECAY_DEC_BR_K0starnn) * weights_K0starnn->GetBinContent(Bin), N_K0starnn); // BR correction * invM correction
    }

    return Correction_Knn * Correction_Kstarnn * Correction_K0nn * Correction_K0starnn;
}

std::vector<EvtInfo> EvtInfos;

struct _BRuncertainty {
    std::vector<int> DMID;
    std::vector<double> RelativeUncertainty;
    std::vector<double> BR_correction_fluctuated;
} BRuncertainty = {
    {
        -1,

        0, 1, 2, 3, 4,
        5, 6, 7, 8, 9,
        10, 11, 12, 13, 14,
        15, 16, 17, 18, 19,
        20, 21, 22, 23, 24,
        25, 26,

        100, 101, 102, 103, 104,
        105, 106, 107, 108, 109,
        110, 111, 112, 113, 114,
        115, 116, 117, 118, 119,
        120, 121, 122
    },
    {
        0,

        0.03943, 0.03913, 0.03913, 0.13433, 0.03943,
        0.02778, 0.26316, 0.46667, 0.02896, 0.80000,
        0.71429, 0.17347, 0.32468, 0.10638, 0.22222,
        1.00000, 0.10000, 0.03809, 0.22018, 0.04336,
        0.11650, 0.29091, 0.21053, 0.07381, 0.16296,
        0.01863, 0.20732,

        0.02414, 0.02414, 0.04018, 0.33333, 0.04018,
        0.20769, 0.15789, 0.13235, 0.15341, 0.56140,
        0.04745, 0.21905, 0.03187, 0.05696, 0.48718,
        0.55000, 0.13750, 0.08333, 0.07910, 0.21622,
        0.04545, 0.11111, 0.07910
    },
    {}
};

enum DecayMode { // reco level
    B2Kc = 0,
    B2KcPi0,
    B2Ks0Pic,
    B2KcPicPic,
    B2Ks0PicPi0,
    B2KcPicPicPi0,
    B2Ks0PicPicPic,
    B2KcPicPicPicPic,
    B2Ks0PicPicPicPi0,
    B2KcPi0Pi0,
    B2Ks0PicPi0Pi0,
    B2KcPicPicPi0Pi0,
    B2KcKcKc,
    B2KcKcKs0Pic,
    B2KcKcKcPi0,
    B02Ks0,
    B02KcPic,
    B02Ks0Pi0,
    B02KcPicPi0,
    B02Ks0PicPic,
    B02KcPicPicPic,
    B02Ks0PicPicPi0,
    B02KcPicPicPicPi0,
    B02Ks0PicPicPicPic,
    B02Ks0Pi0Pi0,
    B02KcPicPi0Pi0,
    B02Ks0PicPicPi0Pi0,
    B02KcKcKs0,
    B02KcKcKcPic,
    B02KcKcKs0Pi0,
    MAX_NUM_DECAYMODE
};
enum DecayModeMC { // MC level
    Xsu2Kc_MC = 0,
    Xsu2Kcstar2KcPi0_MC,
    Xsu2Kcstar2K0Pic_MC,
    Xsu2KcPi0_MC,
    Xsu2K0Pic_MC,
    Xsu2KcPicPic_MC,
    Xsu2K0PicPi0_MC,
    Xsu2KcPicPicPi0_MC,
    Xsu2K0PicPicPic_MC,
    Xsu2KcPicPicPicPic_MC,
    Xsu2K0PicPicPicPi0_MC,
    Xsu2KcPi0Pi0_MC,
    Xsu2K0PicPi0Pi0_MC,
    Xsu2KcPicPicPi0Pi0_MC,
    Xsu2KcKcKc_MC,
    Xsu2KcKcK0Pic_MC,
    Xsu2KcKcKcPi0_MC,
    Xsd2K0_MC,
    Xsd2K0star2KcPic_MC,
    Xsd2K0star2K0Pi0_MC,
    Xsd2KcPic_MC,
    Xsd2K0Pi0_MC,
    Xsd2KcPicPi0_MC,
    Xsd2K0PicPic_MC,
    Xsd2KcPicPicPic_MC,
    Xsd2K0PicPicPi0_MC,
    Xsd2KcPicPicPicPi0_MC,
    Xsd2K0PicPicPicPic_MC,
    Xsd2K0Pi0Pi0_MC,
    Xsd2KcPicPi0Pi0_MC,
    Xsd2K0PicPicPi0Pi0_MC,
    Xsd2KcKcK0_MC,
    Xsd2KcKcKcPic_MC,
    Xsd2KcKcK0Pi0_MC,
    other,
    MAX_NUM_DECAYMODE_MC
};

# define N_decay 38 // five decay mode + others

# define RarityBins 8
# define NToys 1000

/* ====================================== */
// Rarity module
bool sorting(pair<double, double> a, pair<double, double> b) {

    return a.first < b.first;

}

class CDF {
private:
    std::vector<double> m_values;
    std::vector<double> m_weights;
    double m_weights_sum;

    std::vector<double> m_signal_output;
    std::vector<double> m_signal_weight;
    std::vector<double> m_background_output;
    std::vector<double> m_background_weight;

    const double m_max = 1.0;
    const double m_min = 0.92;

    template <typename A, typename B>
    void zip(const std::vector<A>& a, const std::vector<B>& b, std::vector<std::pair<A, B>>& zipped)
    {
        for (size_t i = 0; i < a.size(); ++i)
        {
            zipped.push_back(std::make_pair(a[i], b[i]));
        }
    }

    template <typename A, typename B>
    void unzip(const std::vector<std::pair<A, B>>& zipped, std::vector<A>& a, std::vector<B>& b)
    {
        for (size_t i = 0; i < a.size(); i++)
        {
            a[i] = zipped[i].first;
            b[i] = zipped[i].second;
        }
    }

    void m_load_files(const char* dirname, std::vector<std::string>* names) {
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


public:
    void init(std::vector<double> values, std::vector<double> weights) {
        // Zip the vectors together
        std::vector<std::pair<double, double>> zipped;
        zip(values, weights, zipped);

        sort(zipped.begin(), zipped.end(), sorting);

        // Write the sorted pairs back to the vectors
        unzip(zipped, m_values, m_weights);
        m_weights_sum = std::accumulate(m_weights.begin(), m_weights.end(), 0.0);
    }
    void initbypath(const char* dirname, double weight) {
        std::vector<string> names;

       m_load_files(dirname, &names);

        for (unsigned int i = 0; i < names.size(); i++) {
            float MVA_BB = 0;

            TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");

            TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
            tree_upsilon->SetBranchAddress("MVA_BB", &MVA_BB);

            for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
                tree_upsilon->GetEntry(j);

                m_values.push_back(MVA_BB);
                m_weights.push_back(weight);

            }
            input_file->Close();
        }

        // Zip the vectors together
        std::vector<std::pair<double, double>> zipped;
        zip(m_values, m_weights, zipped);

        sort(zipped.begin(), zipped.end(), sorting);

        // Write the sorted pairs back to the vectors
        unzip(zipped, m_values, m_weights);
        m_weights_sum = std::accumulate(m_weights.begin(), m_weights.end(), 0.0);
    }
    double GetCDFvalue(float value) {
        const double max = 1.0;
        const double min = 0.85;

        auto larger_or_equal = std::lower_bound(m_values.begin(), m_values.end(), value);
        auto larger = std::upper_bound(m_values.begin(), m_values.end(), value);

        int until_index = -1;
        if (larger_or_equal - m_values.begin() == larger - m_values.begin()) { // input value is not same with one of value which is saved
            until_index = larger_or_equal - m_values.begin() - 1;
        }
        else { // input value is same with one of value which is saved
            until_index = larger_or_equal - m_values.begin();
            double weights_partial_sum = std::accumulate(m_weights.begin(), until_index + m_weights.begin(), 0.0);
            return weights_partial_sum / m_weights_sum;
        }

        if (value < m_values.at(0)) {
            int temp_index = std::upper_bound(m_values.begin(), m_values.end(), m_values.at(0)) - m_values.begin() - 1;
            double weights_partial_sum = (value - m_min) * std::accumulate(m_weights.begin(), temp_index + m_weights.begin(), 0.0) / (m_values.at(0) - m_min);
            return weights_partial_sum / m_weights_sum;
        }
        else if (value >= *(m_values.end() - 1)) {
            return 1.0;
        }

        double weights_partial_sum_back = std::accumulate(m_weights.begin(), until_index + m_weights.begin(), 0.0);
        double weights_partial_sum_front = std::accumulate(m_weights.begin(), until_index + m_weights.begin() + 1, 0.0);

        int index_back = until_index;
        int index_front = until_index + 1;

        double weights_partial_sum = weights_partial_sum_back + (value - m_values.at(index_back)) * (weights_partial_sum_front - weights_partial_sum_back) / (m_values.at(index_front) - m_values.at(index_back));

        return weights_partial_sum / m_weights_sum;

    }
    void print() {
        for (unsigned int i = 0; i < m_values.size(); i++) {
            printf("[%lf, %lf] ", m_values.at(i), m_weights.at(i));
        }
        printf("\n");
    }
    void Evaluate(const char* dirname, double weight, bool IsSignal) {
        std::vector<string> names;

        m_load_files(dirname, &names);

        for (unsigned int i = 0; i < names.size(); i++) {
            float MVA_BB = 0;

            TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");

            TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
            tree_upsilon->SetBranchAddress("MVA_BB", &MVA_BB);

            for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
                tree_upsilon->GetEntry(j);

                if (IsSignal) {
                    m_signal_output.push_back(GetCDFvalue(MVA_BB));
                    m_signal_weight.push_back(weight);
                }
                else {
                    m_background_output.push_back(GetCDFvalue(MVA_BB));
                    m_background_weight.push_back(weight);
                }
            }
            input_file->Close();

        }
    }
    void PrintPNG() {
        TH1D* temp_signal_hist = new TH1D("CDF_signal", ";CDF;arbitrary unit", 100, 0.0, 1.0);
        TH1D* temp_background_hist = new TH1D("CDF_background", ";CDF;arbitrary unit", 100, 0.0, 1.0);

        for (unsigned int i = 0; i < m_signal_output.size(); i++) temp_signal_hist->Fill(m_signal_output.at(i), m_signal_weight.at(i));
        for (unsigned int i = 0; i < m_background_output.size(); i++) temp_background_hist->Fill(m_background_output.at(i), m_background_weight.at(i));

        temp_signal_hist->Scale(1.0 / temp_signal_hist->Integral(), "width");
        temp_background_hist->Scale(1.0 / temp_background_hist->Integral(), "width");

        TCanvas* c_temp = new TCanvas("c", "", 800, 800);

        temp_signal_hist->Draw("hist"); c_temp->SaveAs("signal_cdf.png");
        temp_background_hist->Draw("hist"); c_temp->SaveAs("background_cdf.png");

        delete c_temp;
    }
};

CDF cdf;

/* ====================================== */

double GetFEICalFactor(double UpsilonID, double BtagID) {
    // UpsilonID => charged: 0, mixed: 1

    if (UpsilonID > -0.5 && UpsilonID < 0.5) { // charged
        for (int i = 0; i < FEI_cal_Bc_num - 1; i++) {
            if (BtagID > FEI_cal_Bc_modeID[i] - 0.5 && BtagID < FEI_cal_Bc_modeID[i] + 0.5) return FEI_cal_Bc[i];
        }
        return FEI_cal_Bc[FEI_cal_Bc_num - 1];
    }
    else if (UpsilonID > 0.5 && UpsilonID < 1.5) { // mixed
        for (int i = 0; i < FEI_cal_B0_num - 1; i++) {
            if (BtagID > FEI_cal_B0_modeID[i] - 0.5 && BtagID < FEI_cal_B0_modeID[i] + 0.5) return FEI_cal_B0[i];
        }
        return FEI_cal_B0[FEI_cal_B0_num - 1];
    }

    printf("[GetFEICalFactor] error! unexpected decay ID\n");
    exit(1);
    return 0;

}

double GetFEICalFactorUncer(double UpsilonID, double BtagID) {
    // UpsilonID => charged: 0, mixed: 1

    if (UpsilonID > -0.5 && UpsilonID < 0.5) { // charged
        for (int i = 0; i < FEI_cal_Bc_num - 1; i++) {
            if (BtagID > FEI_cal_Bc_modeID[i] - 0.5 && BtagID < FEI_cal_Bc_modeID[i] + 0.5) return FEI_cal_Bc_uncertainty[i];
        }
        return FEI_cal_Bc_uncertainty[FEI_cal_Bc_num - 1];
    }
    else if (UpsilonID > 0.5 && UpsilonID < 1.5) { // mixed
        for (int i = 0; i < FEI_cal_B0_num - 1; i++) {
            if (BtagID > FEI_cal_B0_modeID[i] - 0.5 && BtagID < FEI_cal_B0_modeID[i] + 0.5) return FEI_cal_B0[i];
        }
        return FEI_cal_B0_uncertainty[FEI_cal_B0_num - 1];
    }

    printf("[GetFEICalFactor] error! unexpected decay ID\n");
    exit(1);
    return 0;

}

double GetFEICalFactor_fluctuated(double UpsilonID, double BtagID) {
    // UpsilonID => charged: 0, mixed: 1

    if (UpsilonID > -0.5 && UpsilonID < 0.5) { // charged
        for (int i = 0; i < FEI_cal_Bc_num - 1; i++) {
            if (BtagID > FEI_cal_Bc_modeID[i] - 0.5 && BtagID < FEI_cal_Bc_modeID[i] + 0.5) return FEI_cal_Bc_fluctuated[i];
        }
        return FEI_cal_Bc_fluctuated[FEI_cal_Bc_num - 1];
    }
    else if (UpsilonID > 0.5 && UpsilonID < 1.5) { // mixed
        for (int i = 0; i < FEI_cal_B0_num - 1; i++) {
            if (BtagID > FEI_cal_B0_modeID[i] - 0.5 && BtagID < FEI_cal_B0_modeID[i] + 0.5) return FEI_cal_B0_fluctuated[i];
        }
        return FEI_cal_B0_fluctuated[FEI_cal_B0_num - 1];
    }

    printf("[GetFEICalFactor_fluctuated] error! unexpected decay ID\n");
    exit(1);
    return 0;

}

// global variables to calculate uncertainties
std::vector<double> Ns;
std::vector<int> ntracks;
std::vector<int> npi0s;
std::vector<int> dmID;
std::vector<double> KS0_3D_distance;
std::vector<double> KaonID_correction;
std::vector<double> KaonID_rel_up;
std::vector<double> KaonID_rel_dn;

void load_files(const char* dirname, std::vector<string>* names);

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

double GetNpi0(double Upsilon_ID, double Bsig_ID) {

    if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) return 0.0; // B2Kc
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) return 1.0; // B2KcPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) return 0.0; // B2Ks0Pic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) return 0.0; // B2KcPicPic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) return 1.0; // B2Ks0PicPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) return 1.0; // B2KcPicPicPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) return 0.0; // B2Ks0PicPicPic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) return 0.0; // B2KcPicPicPicPic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) return 1.0; // B2Ks0PicPicPicPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) return 2.0; // B2KcPi0Pi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) return 2.0; // B2Ks0PicPi0Pi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) return 2.0; // B2KcPicPicPi0Pi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) return 0.0; // B2KcKcKc
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) return 0.0; // B2KcKcKs0Pic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) return 1.0; // B2KcKcKcPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) return 0.0; // B02Ks0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) return 0.0; // B02KcPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) return 1.0; // B02Ks0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) return 1.0; // B02KcPicPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) return 0.0; // B02Ks0PicPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) return 0.0; // B02KcPicPicPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) return 1.0; // B02Ks0PicPicPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) return 1.0; // B02KcPicPicPicPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) return 0.0; // B02Ks0PicPicPicPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) return 2.0; // B02Ks0Pi0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) return 2.0; // B02KcPicPi0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) return 2.0; // B02Ks0PicPicPi0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) return 0.0; // B02KcKcKs0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) return 0.0; // B02KcKcKcPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) return 1.0; // B02KcKcKs0Pi0
    else {
        printf("[ERROR] unexpected decay ID\n");
        exit(1);
    }

    return -1.0;
}

double GetNtracks(double Upsilon_ID, double Bsig_ID) {

    if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) return 1.0; // B2Kc
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) return 1.0; // B2KcPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) return 1.0; // B2Ks0Pic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) return 3.0; // B2KcPicPic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) return 1.0; // B2Ks0PicPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) return 3.0; // B2KcPicPicPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) return 3.0; // B2Ks0PicPicPic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) return 5.0; // B2KcPicPicPicPic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) return 3.0; // B2Ks0PicPicPicPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) return 1.0; // B2KcPi0Pi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) return 1.0; // B2Ks0PicPi0Pi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) return 3.0; // B2KcPicPicPi0Pi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) return 3.0; // B2KcKcKc
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) return 3.0; // B2KcKcKs0Pic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) return 3.0; // B2KcKcKcPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) return 0.0; // B02Ks0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) return 2.0; // B02KcPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) return 0.0; // B02Ks0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) return 2.0; // B02KcPicPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) return 2.0; // B02Ks0PicPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) return 4.0; // B02KcPicPicPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) return 2.0; // B02Ks0PicPicPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) return 4.0; // B02KcPicPicPicPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) return 4.0; // B02Ks0PicPicPicPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) return 0.0; // B02Ks0Pi0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) return 2.0; // B02KcPicPi0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) return 2.0; // B02Ks0PicPicPi0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) return 2.0; // B02KcKcKs0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) return 4.0; // B02KcKcKcPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) return 2.0; // B02KcKcKs0Pi0
    else {
        printf("[ERROR] unexpected decay ID\n");
        exit(1);
    }

    return -1.0;
}

void GetNominalNevt(const char* dirname, TH1D* hist, const char* type, const char* sample, double Nevt_nominal[RarityBins * 3], double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get nominal PDF with appropriate correction
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Knn
    B2Kstarnn
    B02K0nn
    B02K0starnn
    otherwise
    */
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    if (strcmp(sample, "CHG") == 0) {}
    else if (strcmp(sample, "MIX") == 0) {}
    else if (strcmp(sample, "SIGNAL") == 0) {}
    else if (strcmp(sample, "Knn") == 0) {}
    else if (strcmp(sample, "Kstarnn") == 0) {}
    else if (strcmp(sample, "K0nn") == 0) {}
    else if (strcmp(sample, "K0starnn") == 0) {}
    else {
        printf("[ERROR] unexpected sample name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKtruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[0][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKmisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[1][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npitruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[2][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npimisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[3][i_PID]);
        }
        for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npi0bin" + std::to_string(i_pi0)).c_str(), &temp_N_bin_pi0[i_pi0]);
        for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeEbin_n" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeE[0][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeEbin_p" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeE[1][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeEbin_n" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeE[2][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeEbin_p" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeE[3][i_fake]);
        }
        for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeMUbin_n" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeMU[0][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeMUbin_p" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeMU[1][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeMUbin_n" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeMU[2][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeMUbin_p" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeMU[3][i_fake]);
        }
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu")) tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);
        if ((CorrectionType == "B2Knn") || (CorrectionType == "B2Kstarnn") || (CorrectionType == "B02K0nn") || (CorrectionType == "B02K0starnn")) {
            tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
            tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
            tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
            tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
            tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
            tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
            tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
            tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);
        }

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu")) tree_Xs->GetEntry(j);

            // do not remove Knn bkg from CHG/MIX generic background. Because of technical hardness, special Knn background is not included for BR syst calculator
            // if ((strcmp(sample, "CHG") == 0) && (N_Knn > MyEPSILON)) continue;
            // else if ((strcmp(sample, "CHG") == 0) && (N_Kstarnn > MyEPSILON)) continue;
            // else if ((strcmp(sample, "MIX") == 0) && (N_K0nn > MyEPSILON)) continue;
            // else if ((strcmp(sample, "MIX") == 0) && (N_K0starnn > MyEPSILON)) continue;

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = GetFEICalFactor(Upsilon_ID, Btag_ID);
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = GetFEICalFactor(Upsilon_ID, Btag_ID);
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(PID_correction[0][i_PID], temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(PID_correction[1][i_PID], temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(PID_correction[2][i_PID], temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(PID_correction[3][i_PID], temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(pi0_correction[i_pi0], temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(PID_fakeE_correction[0][i_fake], temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_fake = Correction_fake * std::pow(PID_fakeE_correction[1][i_fake], temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(PID_fakeE_correction[2][i_fake], temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(PID_fakeE_correction[3][i_fake], temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(PID_fakeMU_correction[0][i_fake], temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_fake = Correction_fake * std::pow(PID_fakeMU_correction[1][i_fake], temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(PID_fakeMU_correction[2][i_fake], temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(PID_fakeMU_correction[3][i_fake], temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }

            double total_weight = weight_var * Correction_pi0 * Correction_FEI * Correction_KID * Correction_PID * Correction_fake;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            if ((CorrectionType == "B2Knn") || (CorrectionType == "B2Kstarnn") || (CorrectionType == "B02K0nn") || (CorrectionType == "B02K0starnn")) total_weight = total_weight * corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

        printf("%s has %lf events (with correction)\n", dirname, Nevt);

    }

    int ArrayBinID = -1;
    
    if (strcmp(sample, "CHG") == 0) ArrayBinID = 0;
    else if (strcmp(sample, "MIX") == 0) ArrayBinID = 1;
    else if (strcmp(sample, "SIGNAL") == 0) ArrayBinID = 2;
    else if (strcmp(sample, "Knn") == 0) ArrayBinID = 0;
    else if (strcmp(sample, "Kstarnn") == 0) ArrayBinID = 0;
    else if (strcmp(sample, "K0nn") == 0) ArrayBinID = 1;
    else if (strcmp(sample, "K0starnn") == 0) ArrayBinID = 1;

    for (int i = 0; i < RarityBins; i++) {
        Nevt_nominal[ArrayBinID * RarityBins + i] = Nevt_nominal[ArrayBinID * RarityBins + i] + hist->GetBinContent(i + 1);
    }

    return;
}

void GetFlucNevt(const char* dirname, TH1D* hist, const char* type, const char* sample, double Nevt_fluc[NToys][RarityBins * 3], int ToyNum, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get nominal PDF with appropriate correction
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Knn
    B2Kstarnn
    B02K0nn
    B02K0starnn
    otherwise
    */
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    if (strcmp(sample, "CHG") == 0) {}
    else if (strcmp(sample, "MIX") == 0) {}
    else if (strcmp(sample, "SIGNAL") == 0) {}
    else if (strcmp(sample, "Knn") == 0) {}
    else if (strcmp(sample, "Kstarnn") == 0) {}
    else if (strcmp(sample, "K0nn") == 0) {}
    else if (strcmp(sample, "K0starnn") == 0) {}
    else {
        printf("[ERROR] unexpected sample name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;

    int __experiment__ = 0;
    int __run__ = 0;
    unsigned int __event__ = 0;
    int __candidate__ = 0;
    int __ncandidates__ = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKtruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[0][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKmisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[1][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npitruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[2][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npimisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[3][i_PID]);
        }
        for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npi0bin" + std::to_string(i_pi0)).c_str(), &temp_N_bin_pi0[i_pi0]);
        for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeEbin_n" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeE[0][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeEbin_p" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeE[1][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeEbin_n" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeE[2][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeEbin_p" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeE[3][i_fake]);
        }
        for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeMUbin_n" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeMU[0][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeMUbin_p" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeMU[1][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeMUbin_n" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeMU[2][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeMUbin_p" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeMU[3][i_fake]);
        }
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu")) tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);
        if ((CorrectionType == "B2Knn") || (CorrectionType == "B2Kstarnn") || (CorrectionType == "B02K0nn") || (CorrectionType == "B02K0starnn")) {
            tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
            tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
            tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
            tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
            tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
            tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
            tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
            tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);
        }

        tree_upsilon->SetBranchAddress("__experiment__", &__experiment__);
        tree_upsilon->SetBranchAddress("__run__", &__run__);
        tree_upsilon->SetBranchAddress("__event__", &__event__);
        tree_upsilon->SetBranchAddress("__candidate__", &__candidate__);
        tree_upsilon->SetBranchAddress("__ncandidates__", &__ncandidates__);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu")) tree_Xs->GetEntry(j);

            // do not remove Knn bkg from CHG/MIX generic background. Because of technical hardness, special Knn background is not included for BR syst calculator
            // if ((strcmp(sample, "CHG") == 0) && (N_Knn > MyEPSILON)) continue;
            // else if ((strcmp(sample, "CHG") == 0) && (N_Kstarnn > MyEPSILON)) continue;
            // else if ((strcmp(sample, "MIX") == 0) && (N_K0nn > MyEPSILON)) continue;
            // else if ((strcmp(sample, "MIX") == 0) && (N_K0starnn > MyEPSILON)) continue;

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = GetFEICalFactor_fluctuated(Upsilon_ID, Btag_ID);
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = GetFEICalFactor_fluctuated(Upsilon_ID, Btag_ID);
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(PID_correction[0][i_PID], temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(PID_correction[1][i_PID], temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(PID_correction[2][i_PID], temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(PID_correction[3][i_PID], temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(pi0_correction[i_pi0], temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(PID_fakeE_correction[0][i_fake], temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_fake = Correction_fake * std::pow(PID_fakeE_correction[1][i_fake], temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(PID_fakeE_correction[2][i_fake], temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(PID_fakeE_correction[3][i_fake], temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(PID_fakeMU_correction[0][i_fake], temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_fake = Correction_fake * std::pow(PID_fakeMU_correction[1][i_fake], temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(PID_fakeMU_correction[2][i_fake], temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(PID_fakeMU_correction[3][i_fake], temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }

            double total_weight = weight_var * Correction_pi0 * Correction_FEI * Correction_KID * Correction_PID * Correction_fake;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            if ((CorrectionType == "B2Knn") || (CorrectionType == "B2Kstarnn") || (CorrectionType == "B02K0nn") || (CorrectionType == "B02K0starnn")) total_weight = total_weight * corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

        printf("%s has %lf events (with correction)\n", dirname, Nevt);

    }

    int ArrayBinID = -1;

    if (strcmp(sample, "CHG") == 0) ArrayBinID = 0;
    else if (strcmp(sample, "MIX") == 0) ArrayBinID = 1;
    else if (strcmp(sample, "SIGNAL") == 0) ArrayBinID = 2;
    else if (strcmp(sample, "Knn") == 0) ArrayBinID = 0;
    else if (strcmp(sample, "Kstarnn") == 0) ArrayBinID = 0;
    else if (strcmp(sample, "K0nn") == 0) ArrayBinID = 1;
    else if (strcmp(sample, "K0starnn") == 0) ArrayBinID = 1;

    for (int i = 0; i < RarityBins; i++) {
        Nevt_fluc[ToyNum][ArrayBinID * RarityBins + i] = Nevt_fluc[ToyNum][ArrayBinID * RarityBins + i] + hist->GetBinContent(i + 1);
    }

    return;
}

void ReadFEIcalFile() {
    const char* FEI_cal_Bp_file = "FEI_cal_Bp_eigen.txt";
    const char* FEI_cal_B0_file = "FEI_cal_B0_eigen.txt";

    FILE* fp_FEI_cal_Bp = fopen(FEI_cal_Bp_file, "r");
    FILE* fp_FEI_cal_B0 = fopen(FEI_cal_B0_file, "r");

    for (int i = 0; i < FEI_cal_Bc_num; i++) {
        fscanf(fp_FEI_cal_Bp, "%lf\n", &FEI_cal_Bc_eigenvalue[i]);
        for (int j = 0; j < FEI_cal_Bc_num; j++) fscanf(fp_FEI_cal_Bp, "%lf\n", &FEI_cal_Bc_eigenvector[i][j]);
    }

    for (int i = 0; i < FEI_cal_B0_num; i++) {
        fscanf(fp_FEI_cal_B0, "%lf\n", &FEI_cal_B0_eigenvalue[i]);
        for (int j = 0; j < FEI_cal_B0_num; j++) fscanf(fp_FEI_cal_B0, "%lf\n", &FEI_cal_B0_eigenvector[i][j]);
    }

    fclose(fp_FEI_cal_Bp);
    fclose(fp_FEI_cal_B0);
}

void ReadPIDFile() {
    const char* KID_true_file = "KaonEff.csv";
    const char* KID_mis_file = "Kaonmis.csv";
    const char* PID_true_file = "PionEff.csv";
    const char* PID_mis_file = "Pionmis.csv";

    FILE* fp_KID_true = fopen(KID_true_file, "r");
    FILE* fp_KID_mis = fopen(KID_mis_file, "r");
    FILE* fp_PID_true = fopen(PID_true_file, "r");
    FILE* fp_PID_mis = fopen(PID_mis_file, "r");

    fscanf(fp_KID_true, "p_min,p_max,cosTheta_min,cosTheta_max,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,threshold,variable\n");
    fscanf(fp_KID_mis, "p_min,p_max,cosTheta_min,cosTheta_max,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,threshold,variable\n");
    fscanf(fp_PID_true, "p_min,p_max,cosTheta_min,cosTheta_max,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,threshold,variable\n");
    fscanf(fp_PID_mis, "p_min,p_max,cosTheta_min,cosTheta_max,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,threshold,variable\n");

    double temp_p_min;
    double temp_p_max;
    double temp_cosTheta_min;
    double temp_cosTheta_max;
    double temp_data_MC_ratio;
    double temp_data_MC_uncertainty_stat_up;
    double temp_data_MC_uncertainty_stat_dn;
    double temp_data_MC_uncertainty_sys_up;
    double temp_data_MC_uncertainty_sys_dn;
    double temp_data_efficiency;
    double temp_data_uncertainty_stat_up;
    double temp_data_uncertainty_stat_dn;
    double temp_data_uncertainty_sys_up;
    double temp_data_uncertainty_sys_dn;
    double temp_MC_efficiency;
    double temp_MC_uncertainty_stat_up;
    double temp_MC_uncertainty_stat_dn;
    double temp_MC_uncertainty_sys_up;
    double temp_MC_uncertainty_sys_dn;
    double temp_threshold;

    for (int i = 0; i < N_PID_syst - 1; i++) {
        fscanf(fp_KID_true, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,kaonID\n", &temp_p_min, &temp_p_max, &temp_cosTheta_min, &temp_cosTheta_max, &PID_correction[0][i], &PID_correction_stat_uncer[0][i], &temp_data_MC_uncertainty_stat_dn, &PID_correction_sys_uncer[0][i], &temp_data_MC_uncertainty_sys_dn, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_threshold);
        fscanf(fp_KID_mis, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,kaonID\n", &temp_p_min, &temp_p_max, &temp_cosTheta_min, &temp_cosTheta_max, &PID_correction[1][i], &PID_correction_stat_uncer[1][i], &temp_data_MC_uncertainty_stat_dn, &PID_correction_sys_uncer[1][i], &temp_data_MC_uncertainty_sys_dn, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_threshold);
        fscanf(fp_PID_true, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,pionID\n", &temp_p_min, &temp_p_max, &temp_cosTheta_min, &temp_cosTheta_max, &PID_correction[2][i], &PID_correction_stat_uncer[2][i], &temp_data_MC_uncertainty_stat_dn, &PID_correction_sys_uncer[2][i], &temp_data_MC_uncertainty_sys_dn, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_threshold);
        fscanf(fp_PID_mis, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,pionID\n", &temp_p_min, &temp_p_max, &temp_cosTheta_min, &temp_cosTheta_max, &PID_correction[3][i], &PID_correction_stat_uncer[3][i], &temp_data_MC_uncertainty_stat_dn, &PID_correction_sys_uncer[3][i], &temp_data_MC_uncertainty_sys_dn, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_threshold);

        if ((std::abs(PID_correction[0][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer[0][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer[0][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction[0][i]) > 10000.0)) {
            PID_correction[0][i] = 1.0;
            PID_correction_stat_uncer[0][i] = 0.0;
            PID_correction_sys_uncer[0][i] = 0.0;
        }
        if ((std::abs(PID_correction[1][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer[1][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer[1][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction[1][i]) > 10000.0)) {
            PID_correction[1][i] = 1.0;
            PID_correction_stat_uncer[1][i] = 0.0;
            PID_correction_sys_uncer[1][i] = 0.0;
        }
        if ((std::abs(PID_correction[2][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer[2][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer[2][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction[2][i]) > 10000.0)) {
            PID_correction[2][i] = 1.0;
            PID_correction_stat_uncer[2][i] = 0.0;
            PID_correction_sys_uncer[2][i] = 0.0;
        }
        if ((std::abs(PID_correction[3][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer[3][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer[3][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction[3][i]) > 10000.0)) {
            PID_correction[3][i] = 1.0;
            PID_correction_stat_uncer[3][i] = 0.0;
            PID_correction_sys_uncer[3][i] = 0.0;
        }

        PID_correction_uncer[0][i] = std::sqrt(PID_correction_stat_uncer[0][i] * PID_correction_stat_uncer[0][i] + PID_correction_sys_uncer[0][i] * PID_correction_sys_uncer[0][i]);
        PID_correction_uncer[1][i] = std::sqrt(PID_correction_stat_uncer[1][i] * PID_correction_stat_uncer[1][i] + PID_correction_sys_uncer[1][i] * PID_correction_sys_uncer[1][i]);
        PID_correction_uncer[2][i] = std::sqrt(PID_correction_stat_uncer[2][i] * PID_correction_stat_uncer[2][i] + PID_correction_sys_uncer[2][i] * PID_correction_sys_uncer[2][i]);
        PID_correction_uncer[3][i] = std::sqrt(PID_correction_stat_uncer[3][i] * PID_correction_stat_uncer[3][i] + PID_correction_sys_uncer[3][i] * PID_correction_sys_uncer[3][i]);
    }

    PID_correction[0][N_PID_syst - 1] = 1.0;
    PID_correction_stat_uncer[0][N_PID_syst - 1] = 0.0;
    PID_correction_sys_uncer[0][N_PID_syst - 1] = 0.0;
    PID_correction_uncer[0][N_PID_syst - 1] = 0.0;

    PID_correction[1][N_PID_syst - 1] = 1.0;
    PID_correction_stat_uncer[1][N_PID_syst - 1] = 0.0;
    PID_correction_sys_uncer[1][N_PID_syst - 1] = 0.0;
    PID_correction_uncer[1][N_PID_syst - 1] = 0.0;

    PID_correction[2][N_PID_syst - 1] = 1.0;
    PID_correction_stat_uncer[2][N_PID_syst - 1] = 0.0;
    PID_correction_sys_uncer[2][N_PID_syst - 1] = 0.0;
    PID_correction_uncer[2][N_PID_syst - 1] = 0.0;

    PID_correction[3][N_PID_syst - 1] = 1.0;
    PID_correction_stat_uncer[3][N_PID_syst - 1] = 0.0;
    PID_correction_sys_uncer[3][N_PID_syst - 1] = 0.0;
    PID_correction_uncer[3][N_PID_syst - 1] = 0.0;

    fclose(fp_KID_true);
    fclose(fp_KID_mis);
    fclose(fp_PID_true);
    fclose(fp_PID_mis);
}

void ReadFakePIDFile() {
    // initialization
    for (int i = 0; i < N_fakeE_syst; i++) {
        PID_fakeE_correction[0][i] = 1.0; //  K-, K+, pi-, pi+
        PID_fakeE_correction[1][i] = 1.0;
        PID_fakeE_correction[2][i] = 1.0;
        PID_fakeE_correction[3][i] = 1.0;

        PID_fakeE_uncer[0][i] = 0.0;
        PID_fakeE_uncer[1][i] = 0.0;
        PID_fakeE_uncer[2][i] = 0.0;
        PID_fakeE_uncer[3][i] = 0.0;
    }

    for (int i = 0; i < N_fakeMU_syst; i++) {
        PID_fakeMU_correction[0][i] = 1.0;
        PID_fakeMU_correction[1][i] = 1.0;
        PID_fakeMU_correction[2][i] = 1.0;
        PID_fakeMU_correction[3][i] = 1.0;

        PID_fakeMU_uncer[0][i] = 0.0;
        PID_fakeMU_uncer[1][i] = 0.0;
        PID_fakeMU_uncer[2][i] = 0.0;
        PID_fakeMU_uncer[3][i] = 0.0;
    }

    const char* K_fromE_file = "kaonID_efficiency_electron.csv";
    const char* K_fromMU_file = "kaonID_efficiency_muon.csv";
    const char* pi_fromE_file = "pionID_efficiency_electron.csv";
    const char* pi_fromMU_file = "pionID_efficiency_muon.csv";

    FILE* fp_K_fromE = fopen(K_fromE_file, "r");
    FILE* fp_K_fromMU = fopen(K_fromMU_file, "r");
    FILE* fp_pi_fromE = fopen(pi_fromE_file, "r");
    FILE* fp_pi_fromMU = fopen(pi_fromMU_file, "r");

    fscanf(fp_K_fromE, "variable,charge,p_min,p_max,theta_min,theta_max,iso_score_min,iso_score_max,working_point,threshold,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn\n");
    fscanf(fp_K_fromMU, "variable,charge,p_min,p_max,theta_min,theta_max,iso_score_min,iso_score_max,working_point,threshold,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn\n");
    fscanf(fp_pi_fromE, "variable,charge,p_min,p_max,theta_min,theta_max,iso_score_min,iso_score_max,working_point,threshold,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn\n");
    fscanf(fp_pi_fromMU, "variable,charge,p_min,p_max,theta_min,theta_max,iso_score_min,iso_score_max,working_point,threshold,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn\n");

    char temp_charge;
    double temp_p_min;
    double temp_p_max;
    double temp_theta_min;
    double temp_theta_max;
    double temp_data_MC_ratio;
    double temp_data_MC_uncertainty_stat_up;
    double temp_data_MC_uncertainty_stat_dn;
    double temp_data_MC_uncertainty_sys_up;
    double temp_data_MC_uncertainty_sys_dn;
    double temp_data_efficiency;
    double temp_data_uncertainty_stat_up;
    double temp_data_uncertainty_stat_dn;
    double temp_data_uncertainty_sys_up;
    double temp_data_uncertainty_sys_dn;
    double temp_MC_efficiency;
    double temp_MC_uncertainty_stat_up;
    double temp_MC_uncertainty_stat_dn;
    double temp_MC_uncertainty_sys_up;
    double temp_MC_uncertainty_sys_dn;
    double temp_threshold;

    // Kaon from fake electron
    while (fscanf(fp_K_fromE, "kaonID,%c,%lf,%lf,%lf,%lf,0.0,1.0,FixedThresh06,0.6,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &temp_charge, &temp_p_min, &temp_p_max, &temp_theta_min, &temp_theta_max, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_data_MC_ratio, &temp_data_MC_uncertainty_stat_up, &temp_data_MC_uncertainty_stat_dn, &temp_data_MC_uncertainty_sys_up, &temp_data_MC_uncertainty_sys_dn) != EOF) {
        int p_bin = -1;
        int theta_bin = -1;

        if (std::abs(temp_p_min - 0.4) < MyEPSILON && std::abs(temp_p_max - 0.5) < MyEPSILON) p_bin = 0;
        else if (std::abs(temp_p_min - 0.5) < MyEPSILON && std::abs(temp_p_max - 1.0) < MyEPSILON) p_bin = 1;
        else if (std::abs(temp_p_min - 1.0) < MyEPSILON && std::abs(temp_p_max - 1.5) < MyEPSILON) p_bin = 2;
        else if (std::abs(temp_p_min - 1.5) < MyEPSILON && std::abs(temp_p_max - 2.0) < MyEPSILON) p_bin = 3;
        else if (std::abs(temp_p_min - 2.0) < MyEPSILON && std::abs(temp_p_max - 2.5) < MyEPSILON) p_bin = 4;
        else if (std::abs(temp_p_min - 2.5) < MyEPSILON && std::abs(temp_p_max - 3.0) < MyEPSILON) p_bin = 5;
        else {
            printf("[ERROR] unknown p bin!\n");
            exit(1);
        }

        if (std::abs(temp_theta_min - 0.22) < MyEPSILON && std::abs(temp_theta_max - 0.56) < MyEPSILON) theta_bin = 0;
        else if (std::abs(temp_theta_min - 0.56) < MyEPSILON && std::abs(temp_theta_max - 1.13) < MyEPSILON) theta_bin = 1;
        else if (std::abs(temp_theta_min - 1.13) < MyEPSILON && std::abs(temp_theta_max - 1.57) < MyEPSILON) theta_bin = 2;
        else if (std::abs(temp_theta_min - 1.57) < MyEPSILON && std::abs(temp_theta_max - 1.88) < MyEPSILON) theta_bin = 3;
        else if (std::abs(temp_theta_min - 1.88) < MyEPSILON && std::abs(temp_theta_max - 2.23) < MyEPSILON) theta_bin = 4;
        else if (std::abs(temp_theta_min - 2.23) < MyEPSILON && std::abs(temp_theta_max - 2.71) < MyEPSILON) theta_bin = 5;
        else {
            printf("[ERROR] unknown theta bin!\n");
            exit(1);
        }

        int bin = theta_bin + 6 * p_bin;

        if (temp_charge == '+') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeE_correction[1][bin] = temp_data_MC_ratio;
            PID_fakeE_uncer[1][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeE_correction[1][bin]) < MyEPSILON) {
                PID_fakeE_correction[1][bin] = 1.0;
                PID_fakeE_uncer[1][bin] = 0.0;
            }

        }
        else if (temp_charge == '-') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeE_correction[0][bin] = temp_data_MC_ratio;
            PID_fakeE_uncer[0][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeE_correction[0][bin]) < MyEPSILON) {
                PID_fakeE_correction[0][bin] = 1.0;
                PID_fakeE_uncer[0][bin] = 0.0;
            }

        }
        else {
            printf("[ERROR] unknown charge!\n");
            exit(1);
        }
    }
    fclose(fp_K_fromE);

    // Pion from fake electron
    while (fscanf(fp_pi_fromE, "pionID,%c,%lf,%lf,%lf,%lf,0.0,1.0,FixedThresh06,0.6,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &temp_charge, &temp_p_min, &temp_p_max, &temp_theta_min, &temp_theta_max, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_data_MC_ratio, &temp_data_MC_uncertainty_stat_up, &temp_data_MC_uncertainty_stat_dn, &temp_data_MC_uncertainty_sys_up, &temp_data_MC_uncertainty_sys_dn) != EOF) {
        int p_bin = -1;
        int theta_bin = -1;

        if (std::abs(temp_p_min - 0.4) < MyEPSILON && std::abs(temp_p_max - 0.5) < MyEPSILON) p_bin = 0;
        else if (std::abs(temp_p_min - 0.5) < MyEPSILON && std::abs(temp_p_max - 1.0) < MyEPSILON) p_bin = 1;
        else if (std::abs(temp_p_min - 1.0) < MyEPSILON && std::abs(temp_p_max - 1.5) < MyEPSILON) p_bin = 2;
        else if (std::abs(temp_p_min - 1.5) < MyEPSILON && std::abs(temp_p_max - 2.0) < MyEPSILON) p_bin = 3;
        else if (std::abs(temp_p_min - 2.0) < MyEPSILON && std::abs(temp_p_max - 2.5) < MyEPSILON) p_bin = 4;
        else if (std::abs(temp_p_min - 2.5) < MyEPSILON && std::abs(temp_p_max - 3.0) < MyEPSILON) p_bin = 5;
        else {
            printf("[ERROR] unknown p bin!\n");
            exit(1);
        }

        if (std::abs(temp_theta_min - 0.22) < MyEPSILON && std::abs(temp_theta_max - 0.56) < MyEPSILON) theta_bin = 0;
        else if (std::abs(temp_theta_min - 0.56) < MyEPSILON && std::abs(temp_theta_max - 1.13) < MyEPSILON) theta_bin = 1;
        else if (std::abs(temp_theta_min - 1.13) < MyEPSILON && std::abs(temp_theta_max - 1.57) < MyEPSILON) theta_bin = 2;
        else if (std::abs(temp_theta_min - 1.57) < MyEPSILON && std::abs(temp_theta_max - 1.88) < MyEPSILON) theta_bin = 3;
        else if (std::abs(temp_theta_min - 1.88) < MyEPSILON && std::abs(temp_theta_max - 2.23) < MyEPSILON) theta_bin = 4;
        else if (std::abs(temp_theta_min - 2.23) < MyEPSILON && std::abs(temp_theta_max - 2.71) < MyEPSILON) theta_bin = 5;
        else {
            printf("[ERROR] unknown theta bin!\n");
            exit(1);
        }

        int bin = theta_bin + 6 * p_bin;

        if (temp_charge == '+') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeE_correction[3][bin] = temp_data_MC_ratio;
            PID_fakeE_uncer[3][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeE_correction[3][bin]) < MyEPSILON) {
                PID_fakeE_correction[3][bin] = 1.0;
                PID_fakeE_uncer[3][bin] = 0.0;
            }

        }
        else if (temp_charge == '-') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeE_correction[2][bin] = temp_data_MC_ratio;
            PID_fakeE_uncer[2][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeE_correction[2][bin]) < MyEPSILON) {
                PID_fakeE_correction[2][bin] = 1.0;
                PID_fakeE_uncer[2][bin] = 0.0;
            }

        }
        else {
            printf("[ERROR] unknown charge!\n");
            exit(1);
        }
    }
    fclose(fp_pi_fromE);

    // Kaon from fake muon
    while (fscanf(fp_K_fromMU, "kaonID,%c,%lf,%lf,%lf,%lf,0.0,1.0,FixedThresh06,0.6,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &temp_charge, &temp_p_min, &temp_p_max, &temp_theta_min, &temp_theta_max, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_data_MC_ratio, &temp_data_MC_uncertainty_stat_up, &temp_data_MC_uncertainty_stat_dn, &temp_data_MC_uncertainty_sys_up, &temp_data_MC_uncertainty_sys_dn) != EOF) {
        int p_bin = -1;
        int theta_bin = -1;

        if (std::abs(temp_p_min - 0.4) < MyEPSILON && std::abs(temp_p_max - 0.5) < MyEPSILON) p_bin = 0;
        else if (std::abs(temp_p_min - 0.5) < MyEPSILON && std::abs(temp_p_max - 0.7) < MyEPSILON) p_bin = 1;
        else if (std::abs(temp_p_min - 0.7) < MyEPSILON && std::abs(temp_p_max - 1.0) < MyEPSILON) p_bin = 2;
        else if (std::abs(temp_p_min - 1.0) < MyEPSILON && std::abs(temp_p_max - 1.5) < MyEPSILON) p_bin = 3;
        else if (std::abs(temp_p_min - 1.5) < MyEPSILON && std::abs(temp_p_max - 2.0) < MyEPSILON) p_bin = 4;
        else if (std::abs(temp_p_min - 2.0) < MyEPSILON && std::abs(temp_p_max - 2.5) < MyEPSILON) p_bin = 5;
        else {
            printf("[ERROR] unknown p bin!\n");
            exit(1);
        }

        if (std::abs(temp_theta_min - 0.4) < MyEPSILON && std::abs(temp_theta_max - 0.64) < MyEPSILON) theta_bin = 0;
        else if (std::abs(temp_theta_min - 0.64) < MyEPSILON && std::abs(temp_theta_max - 0.82) < MyEPSILON) theta_bin = 1;
        else if (std::abs(temp_theta_min - 0.82) < MyEPSILON && std::abs(temp_theta_max - 1.16) < MyEPSILON) theta_bin = 2;
        else if (std::abs(temp_theta_min - 1.16) < MyEPSILON && std::abs(temp_theta_max - 1.46) < MyEPSILON) theta_bin = 3;
        else if (std::abs(temp_theta_min - 1.46) < MyEPSILON && std::abs(temp_theta_max - 1.78) < MyEPSILON) theta_bin = 4;
        else if (std::abs(temp_theta_min - 1.78) < MyEPSILON && std::abs(temp_theta_max - 2.13) < MyEPSILON) theta_bin = 5;
        else if (std::abs(temp_theta_min - 2.13) < MyEPSILON && std::abs(temp_theta_max - 2.22) < MyEPSILON) theta_bin = 6;
        else if (std::abs(temp_theta_min - 2.22) < MyEPSILON && std::abs(temp_theta_max - 2.6) < MyEPSILON) theta_bin = 7;
        else {
            printf("[ERROR] unknown theta bin!\n");
            exit(1);
        }

        int bin = theta_bin + 8 * p_bin;

        if (temp_charge == '+') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeMU_correction[1][bin] = temp_data_MC_ratio;
            PID_fakeMU_uncer[1][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeMU_correction[1][bin]) < MyEPSILON) {
                PID_fakeMU_correction[1][bin] = 1.0;
                PID_fakeMU_uncer[1][bin] = 0.0;
            }

        }
        else if (temp_charge == '-') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeMU_correction[0][bin] = temp_data_MC_ratio;
            PID_fakeMU_uncer[0][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeMU_correction[0][bin]) < MyEPSILON) {
                PID_fakeMU_correction[0][bin] = 1.0;
                PID_fakeMU_uncer[0][bin] = 0.0;
            }

        }
        else {
            printf("[ERROR] unknown charge!\n");
            exit(1);
        }
    }
    fclose(fp_K_fromMU);

    // Pion from fake muon
    while (fscanf(fp_pi_fromMU, "pionID,%c,%lf,%lf,%lf,%lf,0.0,1.0,FixedThresh06,0.6,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &temp_charge, &temp_p_min, &temp_p_max, &temp_theta_min, &temp_theta_max, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_data_MC_ratio, &temp_data_MC_uncertainty_stat_up, &temp_data_MC_uncertainty_stat_dn, &temp_data_MC_uncertainty_sys_up, &temp_data_MC_uncertainty_sys_dn) != EOF) {
        int p_bin = -1;
        int theta_bin = -1;

        if (std::abs(temp_p_min - 0.4) < MyEPSILON && std::abs(temp_p_max - 0.5) < MyEPSILON) p_bin = 0;
        else if (std::abs(temp_p_min - 0.5) < MyEPSILON && std::abs(temp_p_max - 0.7) < MyEPSILON) p_bin = 1;
        else if (std::abs(temp_p_min - 0.7) < MyEPSILON && std::abs(temp_p_max - 1.0) < MyEPSILON) p_bin = 2;
        else if (std::abs(temp_p_min - 1.0) < MyEPSILON && std::abs(temp_p_max - 1.5) < MyEPSILON) p_bin = 3;
        else if (std::abs(temp_p_min - 1.5) < MyEPSILON && std::abs(temp_p_max - 2.0) < MyEPSILON) p_bin = 4;
        else if (std::abs(temp_p_min - 2.0) < MyEPSILON && std::abs(temp_p_max - 2.5) < MyEPSILON) p_bin = 5;
        else {
            printf("[ERROR] unknown p bin!\n");
            exit(1);
        }

        if (std::abs(temp_theta_min - 0.4) < MyEPSILON && std::abs(temp_theta_max - 0.64) < MyEPSILON) theta_bin = 0;
        else if (std::abs(temp_theta_min - 0.64) < MyEPSILON && std::abs(temp_theta_max - 0.82) < MyEPSILON) theta_bin = 1;
        else if (std::abs(temp_theta_min - 0.82) < MyEPSILON && std::abs(temp_theta_max - 1.16) < MyEPSILON) theta_bin = 2;
        else if (std::abs(temp_theta_min - 1.16) < MyEPSILON && std::abs(temp_theta_max - 1.46) < MyEPSILON) theta_bin = 3;
        else if (std::abs(temp_theta_min - 1.46) < MyEPSILON && std::abs(temp_theta_max - 1.78) < MyEPSILON) theta_bin = 4;
        else if (std::abs(temp_theta_min - 1.78) < MyEPSILON && std::abs(temp_theta_max - 2.13) < MyEPSILON) theta_bin = 5;
        else if (std::abs(temp_theta_min - 2.13) < MyEPSILON && std::abs(temp_theta_max - 2.22) < MyEPSILON) theta_bin = 6;
        else if (std::abs(temp_theta_min - 2.22) < MyEPSILON && std::abs(temp_theta_max - 2.6) < MyEPSILON) theta_bin = 7;
        else {
            printf("[ERROR] unknown theta bin!\n");
            exit(1);
        }

        int bin = theta_bin + 8 * p_bin;

        if (temp_charge == '+') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeMU_correction[3][bin] = temp_data_MC_ratio;
            PID_fakeMU_uncer[3][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeMU_correction[3][bin]) < MyEPSILON) {
                PID_fakeMU_correction[3][bin] = 1.0;
                PID_fakeMU_uncer[3][bin] = 0.0;
            }

        }
        else if (temp_charge == '-') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeMU_correction[2][bin] = temp_data_MC_ratio;
            PID_fakeMU_uncer[2][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeMU_correction[2][bin]) < MyEPSILON) {
                PID_fakeMU_correction[2][bin] = 1.0;
                PID_fakeMU_uncer[2][bin] = 0.0;
            }

        }
        else {
            printf("[ERROR] unknown charge!\n");
            exit(1);
        }
    }
    fclose(fp_pi_fromMU);
}

void FluctuateFEIcal() {

    std::normal_distribution<double> FEI_cal_distribution(0.0, 1.0);

    for (int i = 0; i < FEI_cal_Bc_num; i++) {
        double temp_normal = FEI_cal_distribution(generator);

        FEI_cal_Bc_fluctuated[i] = FEI_cal_Bc[i];
        for (int j = 0; j < FEI_cal_Bc_num; j++) FEI_cal_Bc_fluctuated[i] = FEI_cal_Bc_fluctuated[i] + FEI_cal_Bc_eigenvalue[j] * FEI_cal_Bc_eigenvector[j][i] * temp_normal;
    }

    for (int i = 0; i < FEI_cal_B0_num; i++) {
        double temp_normal = FEI_cal_distribution(generator);

        FEI_cal_B0_fluctuated[i] = FEI_cal_B0[i];
        for (int j = 0; j < FEI_cal_B0_num; j++) FEI_cal_B0_fluctuated[i] = FEI_cal_B0_fluctuated[i] + FEI_cal_B0_eigenvalue[j] * FEI_cal_B0_eigenvector[j][i] * temp_normal;
    }

}

void FEI_calculator()
{
    ReadPIDFile();
    ReadFakePIDFile();
    ReadFEIcalFile();

    RooRandom::randomGenerator()->SetSeed(rd());

    double Nevt_nominal[RarityBins * 3] = { 0.0 }; // CHG MIX SIGNAL
    double Nevt_fluc[NToys][RarityBins * 3] = { 0.0 }; // CHG MIX SIGNAL
    double Relative_Uncertainty[NToys][RarityBins * 3] = { 0.0 };
    double Covariance[RarityBins * 3][RarityBins * 3] = { 0.0 };

    /* ====================================== */
    // Seting CDF module
    cdf.initbypath("/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Knunu", Scale_Kplus_test);
    cdf.initbypath("/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Kstarnunu", Scale_Kplusstar_test);
    cdf.initbypath("/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Xsnunu", Scale_Xsu_nonresonant_test);
    cdf.initbypath("/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0nunu", Scale_K0_test);
    cdf.initbypath("/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0starnunu", Scale_K0star_test);
    cdf.initbypath("/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02Xsnunu", Scale_Xsd_nonresonant_test);
    /* ====================================== */



    /* ====================================== */
    // define TH1D for temporary usage
    TH1D* temp_hist = new TH1D("temp_hist", "temp_hist", RarityBins, 0.96, 1.0);
    temp_hist->Reset();
    /* ====================================== */



    /* ====================================== */
    // define path for Ntuple
    const char* MC_dirname_Knunu = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Knunu";
    const char* MC_dirname_Kstarnunu = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Kstarnunu";
    const char* MC_dirname_Xsununu = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Xsnunu";
    const char* MC_dirname_K0nunu = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0nunu";
    const char* MC_dirname_K0starnunu = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0starnunu";
    const char* MC_dirname_Xsdnunu = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02Xsnunu";

    const char* MC_dirname_CHG = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/CHG_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_MIX = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/MIX_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_UUBAR = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/UUBAR_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_DDBAR = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/DDBAR_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_SSBAR = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/SSBAR_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_CHARM = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/CHARM_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    /* ====================================== */



    /* ====================================== */
    // get nominal Nevt
    GetNominalNevt(MC_dirname_Knunu, temp_hist, "Bplus", "SIGNAL", Nevt_nominal, Scale_Kplus_test, "B2Knunu");
    temp_hist->Reset();
    GetNominalNevt(MC_dirname_Kstarnunu, temp_hist, "Bplus", "SIGNAL", Nevt_nominal, Scale_Kplusstar_test, "otherwise");
    temp_hist->Reset();
    GetNominalNevt(MC_dirname_Xsununu, temp_hist, "Bplus", "SIGNAL", Nevt_nominal, Scale_Xsu_nonresonant_test, "otherwise");
    temp_hist->Reset();
    GetNominalNevt(MC_dirname_K0nunu, temp_hist, "Bzero", "SIGNAL", Nevt_nominal, Scale_K0_test, "B02K0nunu");
    temp_hist->Reset();
    GetNominalNevt(MC_dirname_K0starnunu, temp_hist, "Bzero", "SIGNAL", Nevt_nominal, Scale_K0star_test, "otherwise");
    temp_hist->Reset();
    GetNominalNevt(MC_dirname_Xsdnunu, temp_hist, "Bzero", "SIGNAL", Nevt_nominal, Scale_Xsd_nonresonant_test, "otherwise");
    temp_hist->Reset();

    GetNominalNevt(MC_dirname_CHG, temp_hist, "Bplus", "CHG", Nevt_nominal, Scale_CHG_test, "B2Knn");
    temp_hist->Reset();
    GetNominalNevt(MC_dirname_MIX, temp_hist, "Bzero", "MIX", Nevt_nominal, Scale_MIX_test, "B2Knn");
    temp_hist->Reset();
    /* ====================================== */



    /* ====================================== */
    // get fluctuated Nevt for BR
    for (int i = 0; i < NToys; i++) {
        FluctuateFEIcal();

        GetFlucNevt(MC_dirname_Knunu, temp_hist, "Bplus", "SIGNAL", Nevt_fluc, i, Scale_Kplus_test, "B2Knunu");
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_Kstarnunu, temp_hist, "Bplus", "SIGNAL", Nevt_fluc, i, Scale_Kplusstar_test, "otherwise");
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_Xsununu, temp_hist, "Bplus", "SIGNAL", Nevt_fluc, i, Scale_Xsu_nonresonant_test, "otherwise");
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_K0nunu, temp_hist, "Bzero", "SIGNAL", Nevt_fluc, i, Scale_K0_test, "B02K0nunu");
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_K0starnunu, temp_hist, "Bzero", "SIGNAL", Nevt_fluc, i, Scale_K0star_test, "otherwise");
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_Xsdnunu, temp_hist, "Bzero", "SIGNAL", Nevt_fluc, i, Scale_Xsd_nonresonant_test, "otherwise");
        temp_hist->Reset();

        GetFlucNevt(MC_dirname_CHG, temp_hist, "Bplus", "CHG", Nevt_fluc, i, Scale_CHG_test, "B2Knn");
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_MIX, temp_hist, "Bzero", "MIX", Nevt_fluc, i, Scale_MIX_test, "B2Knn");
        temp_hist->Reset();
    }
    /* ====================================== */



    /* ====================================== */
    // get relative uncertainty
    for (int i = 0; i < NToys; i++) {
        for (int j = 0; j < RarityBins * 3; j++) {
            if (std::abs(Nevt_nominal[j]) < MyEPSILON) Relative_Uncertainty[i][j] = 1.0;
            else Relative_Uncertainty[i][j] = Nevt_fluc[i][j] / Nevt_nominal[j];
        }
    }
    /* ====================================== */



    /* ====================================== */
    // file output
    FILE* fp;
    
    fp = fopen("FEI_toys.txt","w");
    for (int i = 0; i < NToys; i++) {
        for (int j = 0; j < RarityBins * 3; j++) {
            fprintf(fp, "%lf ", Relative_Uncertainty[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    /* ====================================== */

    delete temp_hist;
}
