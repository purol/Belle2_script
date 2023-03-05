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

// BKG MC sample number (0.8/ab for BB, 1.0/ab for qq)
# define N_CHG_test 32042497.0
# define N_MIX_test 24693710.0
# define N_UUBAR_test 94447089.0
# define N_DDBAR_test 22664556.0
# define N_SSBAR_test 19244661.0
# define N_CHARM_test 107541168.0
# define N_CHG_train 48052238.0
# define N_MIX_train 37030486.0
# define N_UUBAR_train 141671998.0
# define N_DDBAR_train 34114182.0
# define N_SSBAR_train 28859338.0
# define N_CHARM_train 161280679.0

// new scale factor for BKG MC sample with additional 1/ab (364.436 - 2.763 = 361.673/fb), until LS1
# define Scale_CHG_train ((N_BB_LS1* (BR_BpBp / (BR_BpBp + BR_B0B0))) / (0.8 * N_BpBp_1invab * (N_CHG_train / (N_CHG_train + N_CHG_test)) + N_BpBp_1invab))
# define Scale_MIX_train ((N_BB_LS1* (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (0.8 * N_B0B0_1invab * (N_MIX_train / (N_MIX_train + N_MIX_test)) + N_B0B0_1invab))
# define Scale_UUBAR_train (0.361673/((N_UUBAR_train/(N_UUBAR_train + N_UUBAR_test))*1.0))
# define Scale_DDBAR_train (0.361673/((N_DDBAR_train/(N_DDBAR_train + N_DDBAR_test))*1.0))
# define Scale_SSBAR_train (0.361673/((N_SSBAR_train/(N_SSBAR_train + N_SSBAR_test))*1.0))
# define Scale_CHARM_train (0.361673/((N_CHARM_train/(N_CHARM_train + N_CHARM_test))*1.0))
# define Scale_CHG_test ((N_BB_LS1* (BR_BpBp / (BR_BpBp + BR_B0B0))) / (0.8 * N_BpBp_1invab * (N_CHG_test / (N_CHG_train + N_CHG_test)) + N_BpBp_1invab))
# define Scale_MIX_test ((N_BB_LS1* (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (0.8 * N_B0B0_1invab * (N_MIX_test / (N_MIX_train + N_MIX_test)) + N_B0B0_1invab))
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
# define pi0_correction 0.932
# define pi0_rel_uncertainty ((0.0369 / 0.932) * 100.0) // %
# define Kaon_PID_max_uncertainty 0.1 // not percentage. relative uncertainty
// https://indico.belle2.org/event/6872/contributions/37447/attachments/17127/25504/FEIperformance_B2GM.pdf
# define FEI_cal_Bc 0.679
# define FEI_cal_Bc_uncertainty (0.017/FEI_cal_Bc) // not percentage. relative uncertainty
# define FEI_cal_B0 0.713
# define FEI_cal_B0_uncertainty (0.019/FEI_cal_B0) // not percentage. relative uncertainty

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
double PID_correction[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
double PID_correction_stat_uncer[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
double PID_correction_sys_uncer[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
double PID_correction_uncer[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss

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

std::vector<double> Evt_experiments;
std::vector<double> Evt_runs;
std::vector<double> Evt_events;
std::vector<double> Evt_candidate;
std::vector<double> Evt_ncandidatess;
std::vector<int> Evt_DMID1s;
std::vector<int> Evt_DMID2s;

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

double GetBRRelativeUncertainty(int experiment, int run, int event, int candidate, int ncandidates);

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

void GetNominalNevt(const char* dirname, TH1D* hist, const char* type, const char* sample, double Nevt_nominal[RarityBins * 2], double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get nominal PDF with appropriate correction
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
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
    else {
        printf("[ERROR] unexpected sample name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss

    double invM = -1.0;

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
        for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKtruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[0][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKmisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[1][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npitruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[2][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npimisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[3][i_PID]);
        }
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu")) tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu")) tree_Xs->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

            double Correction_pi0 = std::pow(pi0_correction, Npi0);
            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = FEI_cal_Bc;
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = FEI_cal_B0;
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(PID_correction[0][i_PID], temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(PID_correction[1][i_PID], temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(PID_correction[2][i_PID], temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(PID_correction[3][i_PID], temp_N_bin_PID[3][i_PID]); // mis PID
            }

            double total_weight = weight_var * Correction_pi0 * Correction_FEI * Correction_KID * Correction_PID;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

        printf("%s has %lf events (with correction)\n", dirname, Nevt);

    }

    int ArrayBinID = -1;
    
    if (strcmp(sample, "CHG") == 0) ArrayBinID = 0;
    else if (strcmp(sample, "MIX") == 0) ArrayBinID = 1;

    for (int i = 0; i < RarityBins; i++) {
        Nevt_nominal[ArrayBinID * RarityBins + i] = Nevt_nominal[ArrayBinID * RarityBins + i] + hist->GetBinContent(i + 1);
    }

    return;
}

void GetFlucNevt(const char* dirname, TH1D* hist, const char* type, const char* sample, double Nevt_fluc[NToys][RarityBins * 2], int ToyNum, bool IsItKID, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get nominal PDF with appropriate correction
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
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
    else {
        printf("[ERROR] unexpected sample name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss

    double invM = -1.0;

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
        for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKtruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[0][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKmisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[1][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npitruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[2][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npimisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[3][i_PID]);
        }
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu")) tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

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

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

            double Correction_pi0 = std::pow(pi0_correction, Npi0);
            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = FEI_cal_Bc;
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = FEI_cal_B0;
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(PID_correction[0][i_PID], temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(PID_correction[1][i_PID], temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(PID_correction[2][i_PID], temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(PID_correction[3][i_PID], temp_N_bin_PID[3][i_PID]); // mis PID
            }

            double Correction_BR = GetBRRelativeUncertainty(__experiment__, __run__, __event__, __candidate__, __ncandidates__);

            double total_weight = weight_var * Correction_pi0 * Correction_FEI * Correction_KID * Correction_PID * Correction_BR;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

        printf("%s has %lf events (with correction)\n", dirname, Nevt);

    }

    int ArrayBinID = -1;

    if (strcmp(sample, "CHG") == 0) ArrayBinID = 0;
    else if (strcmp(sample, "MIX") == 0) ArrayBinID = 1;

    for (int i = 0; i < RarityBins; i++) {
        Nevt_fluc[ToyNum][ArrayBinID * RarityBins + i] = Nevt_fluc[ToyNum][ArrayBinID * RarityBins + i] + hist->GetBinContent(i + 1);
    }

    return;
}

void ReadEvtFile() {
    const char* CHG_Evt_file = "CHG_Evt";
    const char* MIX_Evt_file = "MIX_Evt";

    FILE* fp_CHG_Evt = fopen(CHG_Evt_file, "r");
    FILE* fp_MIX_Evt = fopen(MIX_Evt_file, "r");

    double temp_experiment = -1;
    double temp_run = -1;
    double temp_event = -1;
    double temp_candidate = -1;
    double temp_ncandidates = -1;
    int temp_DMID1 = -1;
    int temp_DMID2 = -1;

    while (true) {
        if (fscanf(fp_CHG_Evt, "%lf\n", &temp_experiment) == EOF) break;
        if (fscanf(fp_CHG_Evt, "%lf\n", &temp_run) == EOF) break;
        if (fscanf(fp_CHG_Evt, "%lf\n", &temp_event) == EOF) break;
        if (fscanf(fp_CHG_Evt, "%lf\n", &temp_candidate) == EOF) break;
        if (fscanf(fp_CHG_Evt, "%lf\n", &temp_ncandidates) == EOF) break;
        if (fscanf(fp_CHG_Evt, "%d\n", &temp_DMID1) == EOF) break;
        if (fscanf(fp_CHG_Evt, "%d\n", &temp_DMID2) == EOF) break;

        if ((temp_DMID1 > 99) || (temp_DMID2 > 99)) {
            printf("[ReadEvtFile] DMID for CHG is larger than 99!\n");
            exit(1);
        }

        Evt_experiments.push_back(temp_experiment);
        Evt_runs.push_back(temp_run);
        Evt_events.push_back(temp_event);
        Evt_candidate.push_back(temp_candidate);
        Evt_ncandidatess.push_back(temp_ncandidates);
        Evt_DMID1s.push_back(temp_DMID1);
        Evt_DMID2s.push_back(temp_DMID2);
    }

    while (true) {
        if (fscanf(fp_MIX_Evt, "%lf\n", &temp_experiment) == EOF) break;
        if (fscanf(fp_MIX_Evt, "%lf\n", &temp_run) == EOF) break;
        if (fscanf(fp_MIX_Evt, "%lf\n", &temp_event) == EOF) break;
        if (fscanf(fp_MIX_Evt, "%lf\n", &temp_candidate) == EOF) break;
        if (fscanf(fp_MIX_Evt, "%lf\n", &temp_ncandidates) == EOF) break;
        if (fscanf(fp_MIX_Evt, "%d\n", &temp_DMID1) == EOF) break;
        if (fscanf(fp_MIX_Evt, "%d\n", &temp_DMID2) == EOF) break;

        if (((temp_DMID1 < 100) && (temp_DMID1 >= 0)) || ((temp_DMID2 < 100) && (temp_DMID2 >= 0)) ) {
            printf("[ReadEvtFile] DMID for MIX is smaller than 100!\n");
            exit(1);
        }

        Evt_experiments.push_back(temp_experiment);
        Evt_runs.push_back(temp_run);
        Evt_events.push_back(temp_event);
        Evt_candidate.push_back(temp_candidate);
        Evt_ncandidatess.push_back(temp_ncandidates);
        Evt_DMID1s.push_back(temp_DMID1);
        Evt_DMID2s.push_back(temp_DMID2);
    }

    fclose(fp_CHG_Evt);
    fclose(fp_CHG_Evt);
}

double GetBRRelativeUncertainty(int experiment, int run, int event, int candidate, int ncandidates) {
    int temp_Evt_DMID1 = -100;
    int temp_Evt_DMID2 = -100;

    for (unsigned int i = 0; i < Evt_experiments.size(); i++) {
        if ((experiment == Evt_experiments.at(i)) && (run == Evt_runs.at(i)) && (event == Evt_events.at(i)) && (candidate == Evt_candidate.at(i)) && (ncandidates == Evt_ncandidatess.at(i))) {
            temp_Evt_DMID1 = Evt_DMID1s.at(i);
            temp_Evt_DMID2 = Evt_DMID2s.at(i);
        }
    }

    if ((temp_Evt_DMID1 == -100) || (temp_Evt_DMID2 == -100)) {
        printf("[GetBRRelativeUncertainty] improper DM ID\n");
        exit(1);
    }

    double temp_Relative_uncertainty1 = BRuncertainty.BR_correction_fluctuated.at(std::find(BRuncertainty.DMID.begin(), BRuncertainty.DMID.end(), temp_Evt_DMID1) - BRuncertainty.DMID.begin());
    double temp_Relative_uncertainty2 = BRuncertainty.BR_correction_fluctuated.at(std::find(BRuncertainty.DMID.begin(), BRuncertainty.DMID.end(), temp_Evt_DMID2) - BRuncertainty.DMID.begin());

    return temp_Relative_uncertainty1 * temp_Relative_uncertainty2;

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

void FluctuateBBBR() {

    BRuncertainty.BR_correction_fluctuated.clear();

    for (unsigned int i = 0; i < BRuncertainty.RelativeUncertainty.size(); i++) {
        double Correction_BR = 1.0;

        if (std::abs(BRuncertainty.RelativeUncertainty.at(i)) < MyEPSILON) Correction_BR = 1.0;
        else {
            std::lognormal_distribution<double> BR_distribution(0.0, BRuncertainty.RelativeUncertainty.at(i));
            Correction_BR = BR_distribution(generator);
        }

        BRuncertainty.BR_correction_fluctuated.push_back(Correction_BR);
    }

}

void BR_calculator()
{
    ReadPIDFile();
    ReadEvtFile();

    RooRandom::randomGenerator()->SetSeed(rd());

    double Nevt_nominal[RarityBins * 2] = { 0.0 }; // CHG MIX
    double Nevt_fluc[NToys][RarityBins * 2] = { 0.0 }; // CHG MIX
    double Relative_Uncertainty[NToys][RarityBins * 2] = { 0.0 };
    double Covariance[RarityBins * 2][RarityBins * 2] = { 0.0 };

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
    GetNominalNevt(MC_dirname_CHG, temp_hist, "Bplus", "CHG", Nevt_nominal, Scale_CHG_test, "otherwise");
    temp_hist->Reset();
    GetNominalNevt(MC_dirname_MIX, temp_hist, "Bzero", "MIX", Nevt_nominal, Scale_MIX_test, "otherwise");
    temp_hist->Reset();
    /* ====================================== */



    /* ====================================== */
    // get fluctuated Nevt for BR
    for (int i = 0; i < NToys; i++) {
        FluctuateBBBR();

        GetFlucNevt(MC_dirname_CHG, temp_hist, "Bplus", "CHG", Nevt_fluc, i, true, Scale_CHG_test, "otherwise");
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_MIX, temp_hist, "Bzero", "MIX", Nevt_fluc, i, true, Scale_MIX_test, "otherwise");
        temp_hist->Reset();
    }
    /* ====================================== */



    /* ====================================== */
    // get relative uncertainty
    for (int i = 0; i < NToys; i++) {
        for (int j = 0; j < RarityBins * 2; j++) {
            if (std::abs(Nevt_nominal[j]) < MyEPSILON) Relative_Uncertainty[i][j] = 1.0;
            else Relative_Uncertainty[i][j] = Nevt_fluc[i][j] / Nevt_nominal[j];
        }
    }
    /* ====================================== */



    /* ====================================== */
    // file output
    FILE* fp;
    
    fp = fopen("BR_toys.txt","w");
    for (int i = 0; i < NToys; i++) {
        for (int j = 0; j < RarityBins * 2; j++) {
            fprintf(fp, "%lf ", Relative_Uncertainty[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    /* ====================================== */

    delete temp_hist;
}
