
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

# define Nvar 23
# define DvetoNvar 4

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
# define Scale_CHG_validation ((N_BB_LS1* (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab))
# define Scale_MIX_validation ((N_BB_LS1* (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab))
# define Scale_UUBAR_validation (0.361673)
# define Scale_DDBAR_validation (0.361673)
# define Scale_SSBAR_validation (0.361673)
# define Scale_CHARM_validation (0.361673)

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

    double invM = -1;

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

    if(tempissignal){
        tree_data->SetBranchAddress("nParticlesInList__boB__pl__clKcharge_total__bc", &Decay_Kplus);
        tree_data->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &Decay_Kplusstar_ch1);
        tree_data->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &Decay_Kplusstar_ch2);
        tree_data->SetBranchAddress("nParticlesInList__boB0__clKneutral_total__bc", &Decay_K0);
        tree_data->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &Decay_K0star_ch1);
        tree_data->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &Decay_K0star_ch2);
    }
    tree_data->SetBranchAddress("Bsig_M", &Mxs);

    if(tempissignal) tree_data->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

    int Nevt = 0;
    //printf("%lld entries...\n", tree_data->GetEntries());
    for (unsigned int j = 0; j < tree_data->GetEntries(); j++) { // Fill
        tree_data->GetEntry(j);
        if(tempissignal == true && (Decay_Kplus > 0.5 || Decay_Kplusstar_ch1 > 0.5 || Decay_Kplusstar_ch2 > 0.5 || Decay_K0 > 0.5 || Decay_K0star_ch1 > 0.5 || Decay_K0star_ch2 > 0.5) && Mxs > 1.1) continue;
        else if(tempissignal == true && (Decay_Kplus < 0.5 && Decay_Kplusstar_ch1 < 0.5 && Decay_Kplusstar_ch2 < 0.5 && Decay_K0 < 0.5 &&  Decay_K0star_ch1 < 0.5 && Decay_K0star_ch2 < 0.5) && Mxs < 1.1) continue;
        Nevt++;

        for (unsigned int k = 0; k < Nvar - DvetoNvar; k++) input_vars[k].push_back((float) Vars[k]); 

        if(Dc_chiProb > -0.5){
            input_vars[Nvar - DvetoNvar + 0].push_back((float)Dc_pvalue_std);
            input_vars[Nvar - DvetoNvar + 1].push_back((float)Dc_M);
        }
        else {
            input_vars[Nvar - DvetoNvar + 0].push_back((float)0.0);
            input_vars[Nvar - DvetoNvar + 1].push_back((float) 0.0);
        }
        if(D0_chiProb > -0.5){
            input_vars[Nvar - DvetoNvar + 2].push_back((float)D0_pvalue_std);
            input_vars[Nvar - DvetoNvar + 3].push_back((float) D0_M);
        }
        else {
            input_vars[Nvar - DvetoNvar + 2].push_back((float) 0.0);
            input_vars[Nvar - DvetoNvar + 3].push_back((float) 0.0);
        }

        IsSignal->push_back(tempissignal);

        double FF_reweight = 1.0;
        if (Decay_Kplus > 0.5) FF_reweight = FF_reweight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
        else if (Decay_K0 > 0.5) FF_reweight = FF_reweight * corrector.GetCorrectionFactor(invM * invM, "Bzero");

        weight->push_back(weight_N * FF_reweight);

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

double ObtainWeight(const char* type, const char* MC_version, const char* category, std::string filename) {
    if (strcmp(MC_version, "data") == 0) return 1.0; // no weight if it is data no matter what other values are
    else if ((strcmp(MC_version, "MC15ri") == 0) && (strcmp(category, "test") == 0)) { // MC15ri test
        if ((strcmp(type, "SIGNAL") == 0)) {
            if (filename.find("B2Knunu") != std::string::npos) return Scale_Kplus_test;
            else if (filename.find("B2Kstarnunu") != std::string::npos) return Scale_Kplusstar_test;
            else if (filename.find("B2Xsnunu") != std::string::npos) return Scale_Xsu_nonresonant_test;
            else if (filename.find("B02K0nunu") != std::string::npos) return Scale_Xsu_nonresonant_test;
            else if (filename.find("B02Kstar0nunu") != std::string::npos) return Scale_K0star_test;
            else if (filename.find("B02Xsnunu") != std::string::npos) return Scale_Xsd_nonresonant_test;
            else { printf("[ObtainWeight] undefined type for SIGNAL\n"); exit(1); }
        }
        else if ((strcmp(type, "CHG") == 0)) return Scale_CHG_test;
        else if ((strcmp(type, "MIX") == 0)) return Scale_MIX_test;
        else if ((strcmp(type, "UUBAR") == 0)) return Scale_UUBAR_test;
        else if ((strcmp(type, "DDBAR") == 0)) return Scale_DDBAR_test;
        else if ((strcmp(type, "SSBAR") == 0)) return Scale_SSBAR_test;
        else if ((strcmp(type, "CHARM") == 0)) return Scale_CHARM_test;
        else { printf("[ObtainWeight] undefined type for MC15ri test\n"); exit(1); }
    }
    else if ((strcmp(MC_version, "MC15ri") == 0) && (strcmp(category, "train") == 0)) { // MC15ri train
        if ((strcmp(type, "SIGNAL") == 0)) {
            if (filename.find("B2Knunu") != std::string::npos) return Scale_Kplus_train;
            else if (filename.find("B2Kstarnunu") != std::string::npos) return Scale_Kplusstar_train;
            else if (filename.find("B2Xsnunu") != std::string::npos) return Scale_Xsu_nonresonant_train;
            else if (filename.find("B02K0nunu") != std::string::npos) return Scale_Xsu_nonresonant_train;
            else if (filename.find("B02Kstar0nunu") != std::string::npos) return Scale_K0star_train;
            else if (filename.find("B02Xsnunu") != std::string::npos) return Scale_Xsd_nonresonant_train;
            else { printf("[ObtainWeight] undefined type for SIGNAL\n"); exit(1); }
        }
        else if ((strcmp(type, "CHG") == 0)) return Scale_CHG_train;
        else if ((strcmp(type, "MIX") == 0)) return Scale_MIX_train;
        else if ((strcmp(type, "UUBAR") == 0)) return Scale_UUBAR_train;
        else if ((strcmp(type, "DDBAR") == 0)) return Scale_DDBAR_train;
        else if ((strcmp(type, "SSBAR") == 0)) return Scale_SSBAR_train;
        else if ((strcmp(type, "CHARM") == 0)) return Scale_CHARM_train;
        else { printf("[ObtainWeight] undefined type for MC15ri train\n"); exit(1); }
    }
    else if ((strcmp(MC_version, "MC15ri") == 0) && (strcmp(category, "validation") == 0)) { // MC15ri validation
        if ((strcmp(type, "SIGNAL") == 0)) return 1.0;
        else if ((strcmp(type, "CHG") == 0)) return Scale_CHG_validation;
        else if ((strcmp(type, "MIX") == 0)) return Scale_MIX_validation;
        else if ((strcmp(type, "UUBAR") == 0)) return Scale_UUBAR_validation;
        else if ((strcmp(type, "DDBAR") == 0)) return Scale_DDBAR_validation;
        else if ((strcmp(type, "SSBAR") == 0)) return Scale_SSBAR_validation;
        else if ((strcmp(type, "CHARM") == 0)) return Scale_CHARM_validation;
        else { printf("[ObtainWeight] undefined type for MC15ri validation\n"); exit(1); }
    }

    printf("[ObtainWeight] no matched case!\n");
    exit(1);
    return 1.0;
}

int main(int argc, char* argv[])
{
    // grid search
    // unsigned int nTrees[5] = { 100, 500, 1000, 1500, 2000 };  default is 100
    // unsigned int depth[3] = { 2, 3, 4 };  default is 3 
    // double shrinkage[4] = { 0.05, 0.1, 0.15, 0.2 };  default is 0.1
    // double subsample[5] = { 0.3, 0.4, 0.5, 0.6, 0.7 };  default is 0.5
    // unsigned int binning[4] = { 6, 7, 8, 9 };  default is 2^8 bins per feature

    /*
    * argv[1]: nTrees
    * argv[2]: depth
    * argv[3]: shrinkage path
    * argv[4]: subsample type
    * argv[5]: binning
    * argv[6]: version name (ex. Aqua, Kokoro, Satori, ...)
    * argb[7]: dirname (ex. v000, v001, ...)
    * argv[8]: MC version: {MC15ri|MC15rd}
    */

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
    const char* SIGNAL_input_train = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/SIGNAL_analysis/train_" + std::string(argv[7]) + "/final_output_data";
    const char* CHG_input_train = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/CHG_analysis/train_" + std::string(argv[7]) + "/final_output_data";
    const char* MIX_input_train = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/MIX_analysis/train_" + std::string(argv[7]) + "/final_output_data";
    const char* UUBAR_input_train = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/UUBAR_analysis/train_" + std::string(argv[7]) + "/final_output_data";
    const char* DDBAR_input_train = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/DDBAR_analysis/train_" + std::string(argv[7]) + "/final_output_data";
    const char* SSBAR_input_train = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/SSBAR_analysis/train_" + std::string(argv[7]) + "/final_output_data";
    const char* CHARM_input_train = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/CHARM_analysis/train_" + std::string(argv[7]) + "/final_output_data";

    const char* SIGNAL_input_test = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/SIGNAL_analysis/test_" + std::string(argv[7]) + "/final_output_data";
    const char* CHG_input_test = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/CHG_analysis/test_" + std::string(argv[7]) + "/final_output_data";
    const char* MIX_input_test = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/MIX_analysis/test_" + std::string(argv[7]) + "/final_output_data";
    const char* UUBAR_input_test = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/UUBAR_analysis/test_" + std::string(argv[7]) + "/final_output_data";
    const char* DDBAR_input_test = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/DDBAR_analysis/test_" + std::string(argv[7]) + "/final_output_data";
    const char* SSBAR_input_test = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/SSBAR_analysis/test_" + std::string(argv[7]) + "/final_output_data";
    const char* CHARM_input_test = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/CHARM_analysis/test_" + std::string(argv[7]) + "/final_output_data";



    // define input of the classifier
    std::vector<std::vector<float>> InputVariables;
    std::vector<bool> IsSignal;
    std::vector<float> weight;

    // define input variables
    std::vector<float> input_vars[Nvar];

    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B2Knunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B2Kstarnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B2Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B02K0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B02Kstar0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B02Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(CHG_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHG_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, ObtainWeight("CHG", argv[8], "train", std::string("")));
        }
    }
    {
        std::vector<string> names;
        load_files(MIX_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((MIX_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, ObtainWeight("MIX", argv[8], "train", std::string("")));
        }
    }
    {
        std::vector<string> names;
        load_files(UUBAR_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((UUBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, ObtainWeight("UUBAR", argv[8], "train", std::string("")));
        }
    }
    {
        std::vector<string> names;
        load_files(DDBAR_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((DDBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, ObtainWeight("DDBAR", argv[8], "train", std::string("")));
        }
    }
    {
        std::vector<string> names;
        load_files(SSBAR_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SSBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, ObtainWeight("SSBAR", argv[8], "train", std::string("")));
        }
    }
    {
        std::vector<string> names;
        load_files(CHARM_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHARM_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, ObtainWeight("CHARM", argv[8], "train", std::string("")));
        }
    }

    // convert to double vector
    for (int i = 0; i < Nvar; i++) {
        InputVariables.push_back(input_vars[i]);
    }



    // fit
    classifier.fit(InputVariables, IsSignal, weight);

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
        load_files(SIGNAL_input_test, &names, "B2Knunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, ObtainWeight("SIGNAL", argv[8], "test", std::string("B2Knunu")));
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B2Kstarnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, ObtainWeight("SIGNAL", argv[8], "test", std::string("B2Kstarnunu")));
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B2Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, ObtainWeight("SIGNAL", argv[8], "test", std::string("B2Xsnunu")));
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B02K0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, ObtainWeight("SIGNAL", argv[8], "test", std::string("B02K0nunu")));
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B02Kstar0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, ObtainWeight("SIGNAL", argv[8], "test", std::string("B02Kstar0nunu")));
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test, &names, "B02Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, ObtainWeight("SIGNAL", argv[8], "test", std::string("B02Xsnunu")));
        }
    }
    {
        std::vector<string> names;
        load_files(CHG_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((CHG_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, ObtainWeight("CHG", argv[8], "test", std::string("")));
        }
    }
    {
        std::vector<string> names;
        load_files(MIX_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((MIX_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, ObtainWeight("MIX", argv[8], "test", std::string("")));
        }
    }
    {
        std::vector<string> names;
        load_files(UUBAR_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((UUBAR_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, ObtainWeight("UUBAR", argv[8], "test", std::string("")));
        }
    }
    {
        std::vector<string> names;
        load_files(DDBAR_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((DDBAR_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, ObtainWeight("DDBAR", argv[8], "test", std::string("")));
        }
    }
    {
        std::vector<string> names;
        load_files(SSBAR_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SSBAR_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, ObtainWeight("SSBAR", argv[8], "test", std::string("")));
        }
    }
    {
        std::vector<string> names;
        load_files(CHARM_input_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((CHARM_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, ObtainWeight("CHARM", argv[8], "test", std::string("")));
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

    // clear vector to save memory
    for (unsigned int i = 0; i < InputVariables2.size(); ++i) std::vector<float>().swap(InputVariables2.at(i));
    std::vector<std::vector<float>>().swap(InputVariables2);
    std::vector<bool>().swap(IsSignal2);
    std::vector<float>().swap(weight2);
    for (int i = 0; i < Nvar; i++) std::vector<float>().swap(input_vars2[i]);



    // train sample with proper weight
    std::vector<std::vector<float>> InputVariables3;
    std::vector<bool> IsSignal3;
    std::vector<float> weight3;

    std::vector<float> input_vars3[Nvar];

    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B2Knunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, ObtainWeight("SIGNAL", argv[8], "train", std::string("B2Knunu")));
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B2Kstarnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, ObtainWeight("SIGNAL", argv[8], "train", std::string("B2Kstarnunu")));
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B2Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, ObtainWeight("SIGNAL", argv[8], "train", std::string("B2Xsnunu")));
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B02K0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, ObtainWeight("SIGNAL", argv[8], "train", std::string("B02K0nunu")));
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B02Kstar0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, ObtainWeight("SIGNAL", argv[8], "train", std::string("B02Kstar0nunu")));
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train, &names, "B02Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, ObtainWeight("SIGNAL", argv[8], "train", std::string("B02Xsnunu")));
        }
    }
    {
        std::vector<string> names;
        load_files(CHG_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHG_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, ObtainWeight("CHG", argv[8], "train", std::string("")));
        }
    }
    {
        std::vector<string> names;
        load_files(MIX_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((MIX_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, ObtainWeight("MIX", argv[8], "train", std::string("")));
        }
    }
    {
        std::vector<string> names;
        load_files(UUBAR_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((UUBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, ObtainWeight("UUBAR", argv[8], "train", std::string("")));
        }
    }
    {
        std::vector<string> names;
        load_files(DDBAR_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((DDBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, ObtainWeight("DDBAR", argv[8], "train", std::string("")));
        }
    }
    {
        std::vector<string> names;
        load_files(SSBAR_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SSBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, ObtainWeight("SSBAR", argv[8], "train", std::string("")));
        }
    }
    {
        std::vector<string> names;
        load_files(CHARM_input_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHARM_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, ObtainWeight("CHARM", argv[8], "train", std::string("")));
        }
    }

    // convert to double vector
    for (int i = 0; i < Nvar; i++) {
        InputVariables3.push_back(input_vars3[i]);
    }



    // get FOM for training sample
    // double train_FOM = PrintMaximumFOM(classifier, InputVariables3, IsSignal3, weight3);

    // get AUC for training sample
    double train_AUC = PrintAUC(classifier, InputVariables3, IsSignal3, weight3);

    // clear vector to save memory
    for (unsigned int i = 0; i < InputVariables3.size(); ++i) std::vector<float>().swap(InputVariables3.at(i));
    std::vector<std::vector<float>>().swap(InputVariables3);
    std::vector<bool>().swap(IsSignal3);
    std::vector<float>().swap(weight3);
    for (int i = 0; i < Nvar; i++) std::vector<float>().swap(input_vars3[i]);



    printf("%u_%u_%lf_%lf_%u %lf %lf\n", nTrees, depth, shrinkage, subsample, binning_num, train_AUC, test_AUC);

    FILE* fp;
    fp = fopen(("/home/belle2/junewoo/storage_b1/bsub/Analysis/GridSearch/out/Result_" + std::string(argv[1]) + "_" + std::string(argv[2]) + "_" + std::string(argv[3]) + "_" + std::string(argv[4]) + "_" + std::string(argv[5])).c_str(), "w");
    fprintf(fp, "%u_%u_%lf_%lf_%u %lf %lf\n", nTrees, depth, shrinkage, subsample, binning_num, train_AUC, test_AUC);
    fclose(fp);



    // save model
    std::fstream out_stream(("/home/belle2/junewoo/storage_b1/bsub/Analysis/GridSearch/out/classifier_" + std::string(argv[1]) + "_" + std::string(argv[2]) + "_" + std::string(argv[3]) + "_" + std::string(argv[4]) + "_" + std::string(argv[5])+".weightfile").c_str(), std::ios_base::out | std::ios_base::trunc);
    out_stream << classifier << std::endl;
    out_stream.close();

    return 0;
}
