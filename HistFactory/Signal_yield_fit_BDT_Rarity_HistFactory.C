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
using namespace RooFit;
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

# define Lumi_validation_MC 0.364436 // ab-1

/*
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
*/

// SIGNAL MC sample number
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

// scale factor for SIGNAL MC sample (364.436/fb)
# define Scale_Kplus_train (0.364436 * N_Kplus_nunubar_1invab/N_Kplus_train)
# define Scale_Kplusstar_train (0.364436 * N_Kplusstar_nunubar_1invab/N_Kplusstar_train)
# define Scale_Xsu_nonresonant_train (0.364436 * N_Xsu_nonresonant_nunubar_1invab/N_Xsu_nonresonant_train)
# define Scale_K0_train (0.364436 * N_K0_nunubar_1invab/N_K0_train)
# define Scale_K0star_train (0.364436 * N_K0star_nunubar_1invab/N_K0star_train)
# define Scale_Xsd_nonresonant_train (0.364436 * N_Xsd_nunubar_1invab/N_Xsd_nonresonant_train)
# define Scale_Kplus_test (0.364436 * N_Kplus_nunubar_1invab/N_Kplus_test)
# define Scale_Kplusstar_test (0.364436 * N_Kplusstar_nunubar_1invab/N_Kplusstar_test)
# define Scale_Xsu_nonresonant_test (0.364436 * N_Xsu_nonresonant_nunubar_1invab/N_Xsu_nonresonant_test)
# define Scale_K0_test (0.364436 * N_K0_nunubar_1invab/N_K0_test)
# define Scale_K0star_test (0.364436 * N_K0star_nunubar_1invab/N_K0star_test)
# define Scale_Xsd_nonresonant_test (0.364436 * N_Xsd_nunubar_1invab/N_Xsd_nonresonant_test)

// BKG MC sample number (0.8/ab for BB, 1.0/ab for qq)
# define N_CHG_train 32042497.0
# define N_MIX_train 24693710.0
# define N_UUBAR_train 94447089.0
# define N_DDBAR_train 22664556.0
# define N_SSBAR_train 19244661.0
# define N_CHARM_train 107541168.0
# define N_CHG_test 48052238.0
# define N_MIX_test 37030486.0
# define N_UUBAR_test 141671998.0
# define N_DDBAR_test 34114182.0
# define N_SSBAR_test 28859338.0
# define N_CHARM_test 161280679.0

// scale factor for BKG MC sample (364.436/fb)
# define Scale_CHG_train ((0.364436/0.8)*(N_CHG_train + N_CHG_test)/ N_CHG_train)
# define Scale_MIX_train ((0.364436/0.8)*(N_MIX_train + N_MIX_test)/ N_MIX_train)
# define Scale_UUBAR_train ((0.364436/1.0)*(N_UUBAR_train + N_UUBAR_test)/ N_UUBAR_train)
# define Scale_DDBAR_train ((0.364436/1.0)*(N_DDBAR_train + N_DDBAR_test)/ N_DDBAR_train)
# define Scale_SSBAR_train ((0.364436/1.0)*(N_SSBAR_train + N_SSBAR_test)/ N_SSBAR_train)
# define Scale_CHARM_train ((0.364436/1.0)*(N_CHARM_train + N_CHARM_test)/ N_CHARM_train)
# define Scale_CHG_test ((0.364436/0.8)*(N_CHG_train + N_CHG_test)/ N_CHG_test)
# define Scale_MIX_test ((0.364436/0.8)*(N_MIX_train + N_MIX_test)/ N_MIX_test)
# define Scale_UUBAR_test ((0.364436/1.0)*(N_UUBAR_train + N_UUBAR_test)/ N_UUBAR_test)
# define Scale_DDBAR_test ((0.364436/1.0)*(N_DDBAR_train + N_DDBAR_test)/ N_DDBAR_test)
# define Scale_SSBAR_test ((0.364436/1.0)*(N_SSBAR_train + N_SSBAR_test)/ N_SSBAR_test)
# define Scale_CHARM_test ((0.364436/1.0)*(N_CHARM_train + N_CHARM_test)/ N_CHARM_test)

// scale facto for systematic MC sample
# define N_K0star_nunubar_syst 10000000.0
# define N_Kplusstar_nunubar_syst 10000000.0
# define N_Xsu_nonresonant_nunubar_syst 10000000.0
# define N_Xsd_nonresonant_nunubar_syst 10000000.0

// scale factor for each systematic MC sample (364.436/fb)
# define Scale_Kplusstar_syst (0.364436 * N_Kplusstar_nunubar_1invab/N_Kplusstar_nunubar_syst)
# define Scale_Xsu_nonresonant_syst (0.364436 * N_Xsu_nonresonant_nunubar_1invab/N_Xsu_nonresonant_nunubar_syst)
# define Scale_K0star_syst (0.364436 * N_K0star_nunubar_1invab/N_K0star_nunubar_syst)
# define Scale_Xsd_nonresonant_syst (0.364436 * N_Xsd_nunubar_1invab/N_Xsd_nonresonant_nunubar_syst)

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

void DecayArrayToXsOutputDecay(int Decay[N_decay], double output_Decay[MAX_NUM_DECAYMODE_MC]) {
    output_Decay[Xsu2Kc_MC] = 0;
    output_Decay[Xsu2Kcstar2KcPi0_MC] = 0;
    output_Decay[Xsu2Kcstar2K0Pic_MC] = 0;
    output_Decay[Xsu2KcPi0_MC] = Decay[5];
    output_Decay[Xsu2K0Pic_MC] = Decay[6];
    output_Decay[Xsu2KcPicPic_MC] = Decay[7];
    output_Decay[Xsu2K0PicPi0_MC] = Decay[8];
    output_Decay[Xsu2KcPicPicPi0_MC] = Decay[9];
    output_Decay[Xsu2K0PicPicPic_MC] = Decay[10];
    output_Decay[Xsu2KcPicPicPicPic_MC] = Decay[11];
    output_Decay[Xsu2K0PicPicPicPi0_MC] = Decay[12];
    output_Decay[Xsu2KcPi0Pi0_MC] = Decay[13];
    output_Decay[Xsu2K0PicPi0Pi0_MC] = Decay[14];
    output_Decay[Xsu2KcPicPicPi0Pi0_MC] = Decay[15];
    output_Decay[Xsu2KcKcKc_MC] = Decay[16];
    output_Decay[Xsu2KcKcK0Pic_MC] = Decay[17];
    output_Decay[Xsu2KcKcKcPi0_MC] = Decay[18];
    output_Decay[Xsd2K0_MC] = 0;
    output_Decay[Xsd2K0star2KcPic_MC] = 0;
    output_Decay[Xsd2K0star2K0Pi0_MC] = 0;
    output_Decay[Xsd2KcPic_MC] = Decay[24];
    output_Decay[Xsd2K0Pi0_MC] = Decay[25];
    output_Decay[Xsd2KcPicPi0_MC] = Decay[26];
    output_Decay[Xsd2K0PicPic_MC] = Decay[27];
    output_Decay[Xsd2KcPicPicPic_MC] = Decay[28];
    output_Decay[Xsd2K0PicPicPi0_MC] = Decay[29];
    output_Decay[Xsd2KcPicPicPicPi0_MC] = Decay[30];
    output_Decay[Xsd2K0PicPicPicPic_MC] = Decay[31];
    output_Decay[Xsd2K0Pi0Pi0_MC] = Decay[32];
    output_Decay[Xsd2KcPicPi0Pi0_MC] = Decay[33];
    output_Decay[Xsd2K0PicPicPi0Pi0_MC] = Decay[34];
    output_Decay[Xsd2KcKcK0_MC] = Decay[35];
    output_Decay[Xsd2KcKcKcPic_MC] = Decay[36];
    output_Decay[Xsd2KcKcK0Pi0_MC] = Decay[37];

    output_Decay[other] = 0.0;
    for (int i = 0; i < MAX_NUM_DECAYMODE_MC; i++) {
        if (output_Decay[i] > 0.5) return;
    }

    output_Decay[other] = 1.0;
}

double GetFragmentationWeight(const char* type, DecayModeMC SelectedDecayMode, bool IsItUp) {
    if (strcmp(type, "Bplus") == 0 && ((SelectedDecayMode >= Xsu2KcPi0_MC && SelectedDecayMode <= Xsu2KcKcKcPi0_MC) || (SelectedDecayMode == other))) {
        double Decay_Num[MAX_NUM_DECAYMODE_MC] = { 0.0, 0.0, 0.0, 
            Xsu_frag_decay1, Xsu_frag_decay2, Xsu_frag_decay3, Xsu_frag_decay4, 
            Xsu_frag_decay5, Xsu_frag_decay6, Xsu_frag_decay7, Xsu_frag_decay8, 
            Xsu_frag_decay9, Xsu_frag_decay10, Xsu_frag_decay11, Xsu_frag_decay12, 
            Xsu_frag_decay13, Xsu_frag_decay14, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Xsu_frag_decay30 };
        double Total_Decay_Num = 0;
        for (int i = 0; i < MAX_NUM_DECAYMODE_MC; i++) Total_Decay_Num = Total_Decay_Num + Decay_Num[i];

        double Selected_Decay_Num = 0;
        switch (SelectedDecayMode) {
        case Xsu2KcPi0_MC:
            Selected_Decay_Num = Xsu_frag_decay1;
            break;
        case Xsu2K0Pic_MC:
            Selected_Decay_Num = Xsu_frag_decay2;
            break;
        case Xsu2KcPicPic_MC:
            Selected_Decay_Num = Xsu_frag_decay3;
            break;
        case Xsu2K0PicPi0_MC:
            Selected_Decay_Num = Xsu_frag_decay4;
            break;
        case Xsu2KcPicPicPi0_MC:
            Selected_Decay_Num = Xsu_frag_decay5;
            break;
        case Xsu2K0PicPicPic_MC:
            Selected_Decay_Num = Xsu_frag_decay6;
            break;
        case Xsu2KcPicPicPicPic_MC:
            Selected_Decay_Num = Xsu_frag_decay7;
            break;
        case Xsu2K0PicPicPicPi0_MC:
            Selected_Decay_Num = Xsu_frag_decay8;
            break;
        case Xsu2KcPi0Pi0_MC:
            Selected_Decay_Num = Xsu_frag_decay9;
            break;
        case Xsu2K0PicPi0Pi0_MC:
            Selected_Decay_Num = Xsu_frag_decay10;
            break;
        case Xsu2KcPicPicPi0Pi0_MC:
            Selected_Decay_Num = Xsu_frag_decay11;
            break;
        case Xsu2KcKcKc_MC:
            Selected_Decay_Num = Xsu_frag_decay12;
            break;
        case Xsu2KcKcK0Pic_MC:
            Selected_Decay_Num = Xsu_frag_decay13;
            break;
        case Xsu2KcKcKcPi0_MC:
            Selected_Decay_Num = Xsu_frag_decay14;
            break;
        case other:
            Selected_Decay_Num = Xsu_frag_decay30;
            break;
        default:
            printf("[ERROR] Input value of SelectedDecayMode is not appropriate!\n");
            exit(1);
        }

        if (IsItUp) return ((Total_Decay_Num - 1.5 * Selected_Decay_Num) / (Total_Decay_Num - Selected_Decay_Num));
        else return ((Total_Decay_Num - 0.5 * Selected_Decay_Num) / (Total_Decay_Num - Selected_Decay_Num));
    
    }
    else if (strcmp(type, "Bzero") == 0 && (SelectedDecayMode >= Xsd2KcPic_MC && SelectedDecayMode <= other)) {

        double Decay_Num[MAX_NUM_DECAYMODE_MC] = { 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, Xsd_frag_decay16,
          Xsd_frag_decay17, Xsd_frag_decay18, Xsd_frag_decay19, Xsd_frag_decay20, Xsd_frag_decay21, Xsd_frag_decay22,
          Xsd_frag_decay23, Xsd_frag_decay24, Xsd_frag_decay25, Xsd_frag_decay26, Xsd_frag_decay27,
          Xsd_frag_decay28, Xsd_frag_decay29, Xsd_frag_decay30 };
        double Total_Decay_Num = 0;
        for (int i = 0; i < MAX_NUM_DECAYMODE_MC; i++) Total_Decay_Num = Total_Decay_Num + Decay_Num[i];

        double Selected_Decay_Num = 0;
        switch (SelectedDecayMode) {
        case Xsd2KcPic_MC:
            Selected_Decay_Num = Xsd_frag_decay16;
            break;
        case Xsd2K0Pi0_MC:
            Selected_Decay_Num = Xsd_frag_decay17;
            break;
        case Xsd2KcPicPi0_MC:
            Selected_Decay_Num = Xsd_frag_decay18;
            break;
        case Xsd2K0PicPic_MC:
            Selected_Decay_Num = Xsd_frag_decay19;
            break;
        case Xsd2KcPicPicPic_MC:
            Selected_Decay_Num = Xsd_frag_decay20;
            break;
        case Xsd2K0PicPicPi0_MC:
            Selected_Decay_Num = Xsd_frag_decay21;
            break;
        case Xsd2KcPicPicPicPi0_MC:
            Selected_Decay_Num = Xsd_frag_decay22;
            break;
        case Xsd2K0PicPicPicPic_MC:
            Selected_Decay_Num = Xsd_frag_decay23;
            break;
        case Xsd2K0Pi0Pi0_MC:
            Selected_Decay_Num = Xsd_frag_decay24;
            break;
        case Xsd2KcPicPi0Pi0_MC:
            Selected_Decay_Num = Xsd_frag_decay25;
            break;
        case Xsd2K0PicPicPi0Pi0_MC:
            Selected_Decay_Num = Xsd_frag_decay26;
            break;
        case Xsd2KcKcK0_MC:
            Selected_Decay_Num = Xsd_frag_decay27;
            break;
        case Xsd2KcKcKcPic_MC:
            Selected_Decay_Num = Xsd_frag_decay28;
            break;
        case Xsd2KcKcK0Pi0_MC:
            Selected_Decay_Num = Xsd_frag_decay29;
            break;
        case other:
            Selected_Decay_Num = Xsd_frag_decay30;
            break;
        default:
            printf("[ERROR] Input value of SelectedDecayMode is not appropriate!\n");
            exit(1);
        }

        if (IsItUp) return ((Total_Decay_Num - 1.5 * Selected_Decay_Num) / (Total_Decay_Num - Selected_Decay_Num));
        else return ((Total_Decay_Num - 0.5 * Selected_Decay_Num) / (Total_Decay_Num - Selected_Decay_Num));

    }
    else if (strcmp(type, "Continuum") == 0) {
        printf("[ERROR] unexpected type name\n");
        printf("[ERROR] Continuum type cannot be selected for fragmentation pdf\n");
        exit(1);
    }
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    printf("[ERROR] unexpected error\n");
    exit(1);
    return 0;
}


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

# define RarityBins 10

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
            float MVA_Continuum = 0;

            TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");

            TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
            tree_upsilon->SetBranchAddress("MVA_Continuum", &MVA_Continuum);

            for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
                tree_upsilon->GetEntry(j);

                m_values.push_back(MVA_Continuum);
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
            float MVA_Continuum = 0;

            TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");

            TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
            tree_upsilon->SetBranchAddress("MVA_Continuum", &MVA_Continuum);

            for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
                tree_upsilon->GetEntry(j);

                if (IsSignal) {
                    m_signal_output.push_back(GetCDFvalue(MVA_Continuum));
                    m_signal_weight.push_back(weight);
                }
                else {
                    m_background_output.push_back(GetCDFvalue(MVA_Continuum));
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

void LetsCalculateUncertainties(const char* dirname, double weight) {

    double Upsilon_ID = -1;
    double Bsig_ID = -1;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        double temp_KS0_3D_distance = -1;
        double temp_KaonID_correction = -1;
        double temp_KaonID_rel_up = -1;
        double temp_KaonID_rel_dn = -1;
        double temp_nKaon_excep = -1;

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID); // charged: 0, mixed: 1
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_3D_distance", &temp_KS0_3D_distance);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_correction", &temp_KaonID_correction);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_rel_uncer_up", &temp_KaonID_rel_up);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_rel_uncer_dn", &temp_KaonID_rel_dn);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nKexcep", &temp_nKaon_excep);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) { // B2Kc
                ntracks.push_back(1);
                npi0s.push_back(0);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) { // B2KcPi0
                ntracks.push_back(1);
                npi0s.push_back(1);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) { // B2Ks0Pic
                ntracks.push_back(1);
                npi0s.push_back(0);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) { // B2KcPicPic
                ntracks.push_back(3);
                npi0s.push_back(0);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) { // B2Ks0PicPi0
                ntracks.push_back(1);
                npi0s.push_back(1);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) { // B2KcPicPicPi0
                ntracks.push_back(3);
                npi0s.push_back(1);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) { // B2Ks0PicPicPic
                ntracks.push_back(3);
                npi0s.push_back(0);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) { // B2KcPicPicPicPic
                ntracks.push_back(5);
                npi0s.push_back(0);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) { // B2Ks0PicPicPicPi0
                ntracks.push_back(3);
                npi0s.push_back(1);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) { // B2KcPi0Pi0
                ntracks.push_back(1);
                npi0s.push_back(2);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) { // B2Ks0PicPi0Pi0
                ntracks.push_back(1);
                npi0s.push_back(2);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) { // B2KcPicPicPi0Pi0
                ntracks.push_back(3);
                npi0s.push_back(2);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) { // B2KcKcKc
                ntracks.push_back(3);
                npi0s.push_back(0);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) { // B2KcKcKs0Pic
                ntracks.push_back(3);
                npi0s.push_back(0);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) { // B2KcKcKcPi0
                ntracks.push_back(3);
                npi0s.push_back(1);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) { // B02Ks0
                ntracks.push_back(0);
                npi0s.push_back(0);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) { // B02KcPic
                ntracks.push_back(2);
                npi0s.push_back(0);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) { // B02Ks0Pi0
                ntracks.push_back(0);
                npi0s.push_back(1);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) { // B02KcPicPi0
                ntracks.push_back(2);
                npi0s.push_back(1);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) { // B02Ks0PicPic
                ntracks.push_back(2);
                npi0s.push_back(0);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) { // B02KcPicPicPic
                ntracks.push_back(4);
                npi0s.push_back(0);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) { // B02Ks0PicPicPi0
                ntracks.push_back(2);
                npi0s.push_back(1);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) { // B02KcPicPicPicPi0
                ntracks.push_back(4);
                npi0s.push_back(1);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) { // B02Ks0PicPicPicPic
                ntracks.push_back(4);
                npi0s.push_back(0);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) { // B02Ks0Pi0Pi0
                ntracks.push_back(0);
                npi0s.push_back(2);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) { // B02KcPicPi0Pi0
                ntracks.push_back(2);
                npi0s.push_back(2);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) { // B02Ks0PicPicPi0Pi0
                ntracks.push_back(2);
                npi0s.push_back(2);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) { // B02KcKcKs0
                ntracks.push_back(2);
                npi0s.push_back(0);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) { // B02KcKcKcPic
                ntracks.push_back(4);
                npi0s.push_back(0);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) { // B02KcKcKs0Pi0
                ntracks.push_back(2);
                npi0s.push_back(1);
                dmID.push_back(1);
            }
            else {
                printf("[ERROR] unexpected decay ID\n");
                exit(1);
            }
            KS0_3D_distance.push_back(temp_KS0_3D_distance);

            // correct low momentum region for KaonPID
            KaonID_correction.push_back(temp_KaonID_correction * std::pow(-1, temp_nKaon_excep));
            KaonID_rel_up.push_back(temp_KaonID_rel_up + std::sqrt(2) * temp_nKaon_excep + temp_nKaon_excep * Kaon_PID_max_uncertainty);
            KaonID_rel_dn.push_back(temp_KaonID_rel_dn + std::sqrt(2) * temp_nKaon_excep + temp_nKaon_excep * Kaon_PID_max_uncertainty);

            Ns.push_back(weight);
        }
        input_file->Close();

    }
}

double GetEvtNum(const char* dirname, double weight_var = 1.0) {

    double Total_Evt_Num = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        Total_Evt_Num = Total_Evt_Num + weight_var * tree_upsilon->GetEntries();

        input_file->Close();

    }

    return Total_Evt_Num;

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

double GetNominalPDFs(const char* dirname, TH1D* hist, const char* type, double weight_var = 1.0) { // get nominal PDF with appropriate correction
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_KaonID_correction = -1;
    double temp_nKaon_excep = -1;

    std::vector<string> names;
    load_files(dirname, &names);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("MVA_Continuum", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_correction", &temp_KaonID_correction);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nKexcep", &temp_nKaon_excep);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

            double Correction_pi0 = std::pow(pi0_correction, Npi0);
            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = FEI_cal_Bc;
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = FEI_cal_B0;
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = temp_KaonID_correction * std::pow(-1, temp_nKaon_excep);

            double total_weight = weight_var * Correction_pi0 * Correction_FEI * Correction_KID;

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

        printf("%s has %lf events (with correction)\n", dirname, Nevt);

    }
    return Nevt;
}

double GetFEIPDFs(const char* dirname, TH1D* hist, const char* type, bool IsItUp, double weight_var = 1.0) { // get FEI uncertainty PDF with appropriate correction
    // Be careful! You should combine this function with `GetNominalPDFs` function!
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_KaonID_correction = -1;
    double temp_nKaon_excep = -1;

    std::vector<string> names;
    load_files(dirname, &names);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("MVA_Continuum", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_correction", &temp_KaonID_correction);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nKexcep", &temp_nKaon_excep);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

            double Correction_pi0 = std::pow(pi0_correction, Npi0);
            double Correction_FEI = 1.0;
            if (IsItUp == true) { // + sigma value
                if (strcmp(type, "Bplus") == 0) Correction_FEI = FEI_cal_Bc * (1 + FEI_cal_Bc_uncertainty);
                else if (strcmp(type, "Bzero") == 0) Correction_FEI = FEI_cal_B0 * (1 + FEI_cal_B0_uncertainty);
                else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            }
            else { // - sigma value
                if (strcmp(type, "Bplus") == 0) Correction_FEI = FEI_cal_Bc * (1 - FEI_cal_Bc_uncertainty);
                else if (strcmp(type, "Bzero") == 0) Correction_FEI = FEI_cal_B0 * (1 - FEI_cal_B0_uncertainty);
                else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            }
            double Correction_KID = temp_KaonID_correction * std::pow(-1, temp_nKaon_excep);

            double total_weight = weight_var * Correction_pi0 * Correction_FEI * Correction_KID;

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

        printf("%s has %lf events (with correction)\n", dirname, Nevt);

    }
    return Nevt;
}

double GetPi0PDFs(const char* dirname, TH1D* hist, const char* type, bool IsItUp, double weight_var = 1.0) { // get pi0 uncertainty PDF with appropriate correction
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_KaonID_correction = -1;
    double temp_nKaon_excep = -1;

    std::vector<string> names;
    load_files(dirname, &names);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("MVA_Continuum", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_correction", &temp_KaonID_correction);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nKexcep", &temp_nKaon_excep);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

            double Correction_pi0 = 1.0;
            if (IsItUp == true) Correction_pi0 = std::pow(pi0_correction * (1 + (pi0_rel_uncertainty/100.0) ), Npi0); // + sigma value
            else Correction_pi0 = std::pow(pi0_correction * (1 - (pi0_rel_uncertainty / 100.0)), Npi0); // - sigma value
            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = FEI_cal_Bc;
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = FEI_cal_B0;
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = temp_KaonID_correction * std::pow(-1, temp_nKaon_excep);

            double total_weight = weight_var * Correction_pi0 * Correction_FEI * Correction_KID;

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

        printf("%s has %lf events (with correction)\n", dirname, Nevt);

    }
    return Nevt;
}

double GetTrackPDFs(const char* dirname, TH1D* hist, const char* type, bool IsItUp, double weight_var = 1.0) { // get track efficiency uncertainty PDF with appropriate correction
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_KaonID_correction = -1;
    double temp_nKaon_excep = -1;

    std::vector<string> names;
    load_files(dirname, &names);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("MVA_Continuum", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_correction", &temp_KaonID_correction);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nKexcep", &temp_nKaon_excep);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);
            double Ntrack = GetNtracks(Upsilon_ID, Bsig_ID);

            double Correction_pi0 = std::pow(pi0_correction, Npi0);
            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = FEI_cal_Bc;
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = FEI_cal_B0;
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = temp_KaonID_correction * std::pow(-1, temp_nKaon_excep);
            double Correction_track = 1.0;
            if (IsItUp == true) Correction_track = (1 + (Ntrack * track_rel_uncertainty / 100.0));
            else Correction_track = (1 - (Ntrack * track_rel_uncertainty / 100.0));

            double total_weight = weight_var * Correction_pi0 * Correction_FEI * Correction_KID * Correction_track;

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

        printf("%s has %lf events (with correction)\n", dirname, Nevt);

    }
    return Nevt;
}

double GetKS0PDFs(const char* dirname, TH1D* hist, const char* type, bool IsItUp, double weight_var = 1.0) { // get KS0 uncertainty PDF with appropriate correction
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_KaonID_correction = -1;
    double temp_nKaon_excep = -1;
    double temp_KS0_3D_distance = -1;

    std::vector<string> names;
    load_files(dirname, &names);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("MVA_Continuum", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_correction", &temp_KaonID_correction);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nKexcep", &temp_nKaon_excep);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_3D_distance", &temp_KS0_3D_distance);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

            double Correction_pi0 = std::pow(pi0_correction, Npi0);
            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = FEI_cal_Bc;
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = FEI_cal_B0;
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = temp_KaonID_correction * std::pow(-1, temp_nKaon_excep);
            double KS0_correction = 1.0;
            if(IsItUp == true) KS0_correction = 1 + (KS0_rel_uncertainty * temp_KS0_3D_distance / 100.0);
            else KS0_correction = 1 - (KS0_rel_uncertainty * temp_KS0_3D_distance / 100.0);

            double total_weight = weight_var * Correction_pi0 * Correction_FEI * Correction_KID * KS0_correction;

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

        printf("%s has %lf events (with correction)\n", dirname, Nevt);

    }
    return Nevt;
}

double GetKIDPDFs(const char* dirname, TH1D* hist, const char* type, bool IsItUp, double weight_var = 1.0) { // get nominal PDF with appropriate correction
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_KaonID_correction = -1;
    double temp_nKaon_excep = -1;
    double temp_KaonID_rel_up = -1;
    double temp_KaonID_rel_dn = -1;

    std::vector<string> names;
    load_files(dirname, &names);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("MVA_Continuum", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_correction", &temp_KaonID_correction);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nKexcep", &temp_nKaon_excep);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_rel_uncer_up", &temp_KaonID_rel_up);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_rel_uncer_dn", &temp_KaonID_rel_dn);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

            double Correction_pi0 = std::pow(pi0_correction, Npi0);
            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = FEI_cal_Bc;
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = FEI_cal_B0;
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = temp_KaonID_correction * std::pow(-1, temp_nKaon_excep);
            if (IsItUp == true) {
                double KaonID_rel_up = ( temp_KaonID_rel_up + std::sqrt(2) * temp_nKaon_excep ) + temp_nKaon_excep * Kaon_PID_max_uncertainty;
                Correction_KID = Correction_KID + Correction_KID * KaonID_rel_up;
            }
            else {
                double KaonID_rel_dn = (temp_KaonID_rel_dn + std::sqrt(2) * temp_nKaon_excep) + temp_nKaon_excep * Kaon_PID_max_uncertainty;
                Correction_KID = Correction_KID - Correction_KID * KaonID_rel_dn;
            }

            double total_weight = weight_var * Correction_pi0 * Correction_FEI * Correction_KID;

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

        printf("%s has %lf events (with correction)\n", dirname, Nevt);

    }
    return Nevt;
}

void GetKffPDFs(const char* dirname, TH1D* hist[7], double Correction_factor_BR[7], const char* type, int charge, double weight_var = 1.0) { // get Kff uncertainty PDF with appropriate correction
    // Be careful! You should combine this function with `GetNominalPDFs` function!
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    double m_b = -1; // B meson mass
    double m_k = -1; // Kaon mass

    const double alpha0 = 0.432; // +-0.011
    const double alpha1 = -0.664; // +-0.096
    const double alpha2 = -1.2; // +-0.69
    const double alpha0_sigma = 0.011;
    const double alpha1_sigma = -0.664;
    const double alpha2_sigma = -1.2;
    /*
    <PCA>
    e0 = -0.00586 alpha0 + 0.03681 alpha1 + 0.99931 alpha2
    e1 = 0.05181 alpha0 + 0.99799 alpha1 - 0.03646 alpha2
    e2 = 0.99864 alpha0 -0.05156 alpha2 + 0.00776 alpha2
    lambda0 = 0.69047
    lambda1 = 0.09274
    lambda2 = 0.01007
    */
    const double LinearCoefficients[3][3] = {
        {-0.005860659921244,   0.051812539908648,   0.998639631385367},
        {0.036811794399807,   0.997991063981209, -0.051562854912436},
        {0.999305030738039, -0.036459524432847,   0.007756198796883}
    };
    const double Lambdas[7][3] = {
        {0.0, 0.0, 0.0},
        {0.690465567619637, 0.0, 0.0},
        {-0.690465567619637, 0.0, 0.0},
        {0.0, 0.092737087984295, 0.0},
        {0.0, -0.092737087984295, 0.0},
        {0.0, 0.0, 0.010065727287981},
        {0.0, 0.0, -0.010065727287981}
    };
    double fluctuations[7][3] = { 0.0 };
    for (int i = 0; i < 7; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) fluctuations[i][j] = fluctuations[i][j] + Lambdas[i][k] * LinearCoefficients[k][j];
        }
    }
    /*
    fluctuations[0][0] = 0.0; // alpha0 change
    fluctuations[0][1] = 0.0; // alpha1 change
    fluctuations[0][2] = 0.0; // alpha2 change
    fluctuations[1][0] = 0.69047 * (-0.00586); // alpha0 change Lambdas[1][0]*LinearCoefficients[0][0] + Lambdas[1][1]*LinearCoefficients[1][0] + Lambdas[1][2]*LinearCoefficients[2][0]
    fluctuations[1][1] = 0.69047 * 0.03681; // alpha1 change Lambdas[1][0]*LinearCoefficients[0][1] + Lambdas[1][1]*LinearCoefficients[1][1] + Lambdas[1][2]*LinearCoefficients[2][1]
    fluctuations[1][2] = 0.69047 * 0.99931; // alpha2 change
    fluctuations[2][0] = (-0.69047) * (-0.00586); // alpha0 change
    fluctuations[2][1] = (-0.69047) * 0.03681; // alpha1 change
    fluctuations[2][2] = (-0.69047) * 0.99931; // alpha2 change
    fluctuations[3][0] = 0.09274 * 0.05181; // alpha0 change
    fluctuations[3][1] = 0.09274 * 0.99799; // alpha1 change
    fluctuations[3][2] = 0.09274 * (-0.03646); // alpha2 change
    fluctuations[4][0] = (-0.09274) * 0.05181; // alpha0 change
    fluctuations[4][1] = (-0.09274) * 0.99799; // alpha1 change
    fluctuations[4][2] = (-0.09274) * (-0.03646); // alpha2 change
    fluctuations[5][0] = 0.01007 * 0.99864; // alpha0 change
    fluctuations[5][1] = 0.01007 * (-0.05156); // alpha1 change
    fluctuations[5][2] = 0.01007 * 0.00775; // alpha2 change
    fluctuations[6][0] = 0.01007 * 0.99864; // alpha0 change
    fluctuations[6][1] = 0.01007 * (-0.05156); // alpha1 change
    fluctuations[6][2] = 0.01007 * 0.00775; // alpha2 change
    */
    double value[7] = { 0.0 }; // value of lambda^1.5 * fp*fp
    double Nevts[7] = { 0.0 }; // number of events at each fluctuations
    double tot_value[7] = { 0.0 }; // value of integral lambda^1.5 * fp*fp

    double q2 = -1;

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_KaonID_correction = -1;
    double temp_nKaon_excep = -1;

    std::vector<string> names;
    load_files(dirname, &names);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("MVA_Continuum", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_correction", &temp_KaonID_correction);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nKexcep", &temp_nKaon_excep);

        TTree* tree_Xs = (TTree*)input_file->Get("Xs");

        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &q2);
        if (charge == 0) {
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spM__bc", &m_b);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &m_k);
        }
        else {
            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spM__bc", &m_b);
            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &m_k);
        }

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            tree_Xs->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

            double Correction_pi0 = std::pow(pi0_correction, Npi0);
            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = FEI_cal_Bc;
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = FEI_cal_B0;
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = temp_KaonID_correction * std::pow(-1, temp_nKaon_excep);

            q2 = q2 * q2;

            for (int k = 0; k < 7; k++) {
                const double alpha0_fluc = alpha0 + fluctuations[k][0];
                const double alpha1_fluc = alpha1 + fluctuations[k][1];
                const double alpha2_fluc = alpha2 + fluctuations[k][2];

                double mp = m_b + 0.046;
                double tp = (m_b + m_k) * (m_b + m_k);
                double tm = (m_b - m_k) * (m_b - m_k);
                double t0 = tp * (1 - sqrt(1 - tm / tp));
                double z = (sqrt(tp - q2) - sqrt(tp - t0)) / (sqrt(tp - q2) + sqrt(tp - t0));
                double fp = (1 / (1 - q2 / (mp * mp))) * (alpha0_fluc + alpha1_fluc * z + alpha2_fluc * z * z + (-alpha1_fluc + 2 * alpha2_fluc) * z * z * z / 3);
                double lambda = (m_b * m_b * m_b * m_b) + (m_k * m_k * m_k * m_k) + (q2 * q2) - 2 * (m_b * m_b * m_k * m_k + m_b * m_b * q2 + m_k * m_k * q2);

                value[k] = std::pow(lambda, 1.5) * fp * fp;
                double total_weight = weight_var * Correction_pi0 * Correction_FEI * Correction_KID * (value[k] / value[0]);
                hist[k]->Fill(MVA_var, total_weight );
                Nevts[k] = Nevts[k] + total_weight;
            }

        }
        input_file->Close();
    }


    // calculate fluctuation of total BR
    for (int k = 0; k < 7; k++) {
        const int step = 2000;
        double m_b_fix = -1;
        double m_k_fix = -1;
        if (charge == 0) {
            m_b_fix = 5.27965;
            m_k_fix = 0.497611;
        }
        else if (charge == 1 || charge == -1) {
            m_b_fix = 5.27934;
            m_k_fix = 0.493677;
        }
        const double q2_min = 0;
        const double q2_max = (m_b_fix - m_k_fix) * (m_b_fix - m_k_fix);
        const double delq2 = (q2_max - q2_min) / step; // max: (mb-mk)^2

        const double alpha0_fluc = alpha0 + fluctuations[k][0];
        const double alpha1_fluc = alpha1 + fluctuations[k][1];
        const double alpha2_fluc = alpha2 + fluctuations[k][2];

        for (int i = 0; i < step; i++) {
            const double q2_integral = q2_min + delq2 * i;

            double mp = m_b_fix + 0.046;
            double tp = (m_b_fix + m_k_fix) * (m_b_fix + m_k_fix);
            double tm = (m_b_fix - m_k_fix) * (m_b_fix - m_k_fix);
            double t0 = tp * (1 - sqrt(1 - tm / tp));
            double z = (sqrt(tp - q2_integral) - sqrt(tp - t0)) / (sqrt(tp - q2_integral) + sqrt(tp - t0));
            double fp = (1 / (1 - q2_integral / (mp * mp))) * (alpha0_fluc + alpha1_fluc * z + alpha2_fluc * z * z + (-alpha1_fluc + 2 * alpha2_fluc) * z * z * z / 3);
            double lambda = (m_b_fix * m_b_fix * m_b_fix * m_b_fix) + (m_k_fix * m_k_fix * m_k_fix * m_k_fix) + (q2_integral * q2_integral) - 2 * (m_b_fix * m_b_fix * m_k_fix * m_k_fix + m_b_fix * m_b_fix * q2_integral + m_k_fix * m_k_fix * q2_integral);

            tot_value[k] = tot_value[k] + std::pow(lambda, 1.5) * fp * fp * delq2;
        }
        if (charge == 0) {
            Correction_factor_BR[k] = tot_value[k] / tot_value[0];
        }
        else if (charge == 1 || charge == -1) {
            Correction_factor_BR[k] = tot_value[k] / tot_value[0];
        }
    }


    return;
}

void GetKstarffPDFs(const char* dirname, TH1D* hist[19], double Correction_factor_BR[19], const char* type, int charge, double weight_var = 1.0) { // get Kstarff uncertainty PDF with appropriate correction
    // Be careful! You should combine this function with `GetNominalPDFs` function!
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    double m_b = -1; // B meson mass
    double m_k = -1; // Kaon star mass
    double costheta = -100; // costheta

    const double alpha0_A1 = 0.3;
    const double alpha1_A1 = 0.39;
    const double alpha2_A1 = 1.19;
    const double alpha0_A12 = 0.27;
    const double alpha1_A12 = 0.53;
    const double alpha2_A12 = 0.48;
    const double alpha0_v0 = 0.38;
    const double alpha1_v0 = -1.17;
    const double alpha2_v0 = 2.42;
    const double mR_A1 = 5.829;
    const double mR_A12 = 5.829;
    const double mR_v0 = 5.415;
    /*
    <PCA>
    e0 = - 0.00477 A1a0 + 0.01377 A1a1 + 0.29659 A1a2 + 0.00314 A12a0 + 0.02574 A12a1 + 0.11862 A12a2 - 0.00662 Va0 + 0.02683 Va1 + 0.94674 Va2
    e1 = 0.00955 A1a0 + 0.14336 A1a1 + 0.83303 A1a2 + 0.00253 A12a0 + 0.04200 A12a1 + 0.40190 A12a2 + 0.01359 Va0 + 0.14322 Va1 - 0.31847 Va2
    e2 = 0.00490 A1a0 - 0.02719 A1a1 - 0.43464 A1a2 - 0.00490 A12a0 + 0.14067 A12a1 + 0.88604 A12a2 + 0.00321 Va0 + 0.07099 Va1 + 0.01977 Va2
    e = 0.A1a0 + A1a1 + A1a2 + A12a0 + A12a1 + A12a2 + Va0 + Va1 + Va2
    e = A1a0 + A1a1 + A1a2 + A12a0 + A12a1 + A12a2 + Va0 + Va1 + Va2
    e = A1a0 + A1a1 + A1a2 + A12a0 + A12a1 + A12a2 + Va0 + Va1 + Va2
    e = A1a0 + A1a1 + A1a2 + A12a0 + A12a1 + A12a2 + Va0 + Va1 + Va2
    e = A1a0 + A1a1 + A1a2 + A12a0 + A12a1 + A12a2 + Va0 + Va1 + Va2
    e = A1a0 + A1a1 + A1a2 + A12a0 + A12a1 + A12a2 + Va0 + Va1 + Va2
    lambda0 = 1.57728
    lambda1 = 1.06016
    lambda2 = 0.52088
    lambda3 = 0.23284
    lambda4 = 0.08996
    lambda5 = 0.05088
    lambda6 = 0.01842
    lambda7 = 0.00740
    lambda8 = 0.00177
    */
    const double LinearCoefficients[9][9] = {
    {-0.004767833929074,   0.013772238042338,   0.296589616507393,   0.003139670856264,   0.025743553369441,   0.118620792973715, -0.006615160985686,   0.026825078669141,   0.946739334135628},
    {0.009545743410425,   0.143364495729883,   0.833028724155009,   0.002528027950777,   0.041995533293814,   0.401904805723855,   0.013592144331223,   0.143215823688685, -0.318473920017488},
    {0.004903939104258, -0.027192276054885, -0.434636497785463, -0.004895872506074,   0.140668610588867,   0.886044604252808,   0.003208697626475,   0.070992606469561,   0.019766944252755},
    {0.081726107544427,   0.418943528790113, -0.148580744618415,   0.008936060075137,   0.016121321335068, -0.134221136343898,   0.101230277583689,   0.875232126303976,   0.033121258010682},
    {-0.047393076509321, -0.020252979416420,   0.020465460529877,   0.194484738990744,   0.966871985468078, -0.140890698347323, -0.062525020565583, -0.015984271025904, -0.015622549984120},
    {0.130012725270740,   0.880660653289362, -0.078326566909849, -0.023912474389986,   0.027560190412870,   0.018608679240451,   0.008992181636485, -0.446264603521217,   0.022087208404295},
    {0.512122358575286, -0.109791833791786,   0.012308677320526,   0.445594862213020, -0.020341462354386,   0.009165569689839,   0.721186808802633, -0.079569219353840,   0.005540910930872},
    {-0.519183885001133,   0.030474942310095,   0.001530220137484, -0.500080935622905,   0.116801392077100, -0.015739303946084,   0.680871786694491, -0.044168374846196,   0.002926015416780},
    {-0.664990460261152,   0.118522184960495, -0.020121877756987,   0.716134125194759, -0.166939357181738,   0.027983481968456,   0.042738598277210, -0.002952732689141,   0.000271198372049}
    };
    const double Lambdas[19][9] = {
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {1.577276362747197, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {-1.577276362747197, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 1.060164618950133, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, -1.060164618950133, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.520878893100167, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, -0.520878893100167, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.232840750162789, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, -0.232840750162789, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.089960782159514, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, -0.089960782159514, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.050880752840377, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, -0.050880752840377, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.018416666017616, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.018416666017616, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.007399540696482, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.007399540696482, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001770848832452},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.001770848832452}
    };
    double fluctuations[19][9] = { 0.0 };
    for (int i = 0; i < 19; i++) {
        for (int j = 0; j < 9; j++) {
            for (int k = 0; k < 9; k++) fluctuations[i][j] = fluctuations[i][j] + Lambdas[i][k] * LinearCoefficients[k][j];
        }
    }

    double value[19] = { 0.0 }; // 
    double Nevts[19] = { 0.0 }; // number of events at each fluctuations
    double tot_value[19] = { 0.0 }; // value of integral dBR/dsdcostheta

    double q2 = -1;

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_KaonID_correction = -1;
    double temp_nKaon_excep = -1;

    std::vector<string> names;
    load_files(dirname, &names);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("MVA_Continuum", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_correction", &temp_KaonID_correction);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nKexcep", &temp_nKaon_excep);

        TTree* tree_Xs = (TTree*)input_file->Get("Xs");

        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &q2);
        if (charge == 0) {
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spM__bc", &m_b);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &m_k);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", &costheta);
        }
        else {
            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spM__bc", &m_b);
            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &m_k);
            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", &costheta);
        }

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            tree_Xs->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

            double Correction_pi0 = std::pow(pi0_correction, Npi0);
            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = FEI_cal_Bc;
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = FEI_cal_B0;
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = temp_KaonID_correction * std::pow(-1, temp_nKaon_excep);

            q2 = q2 * q2;

            for (int k = 0; k < 19; k++) {
                const double alpha0_A1_fluc = alpha0_A1 + fluctuations[k][0];
                const double alpha1_A1_fluc = alpha1_A1 + fluctuations[k][1];
                const double alpha2_A1_fluc = alpha2_A1 + fluctuations[k][2];
                const double alpha0_A12_fluc = alpha0_A12 + fluctuations[k][3];
                const double alpha1_A12_fluc = alpha1_A12 + fluctuations[k][4];
                const double alpha2_A12_fluc = alpha2_A12 + fluctuations[k][5];
                const double alpha0_v0_fluc = alpha0_v0 + fluctuations[k][6];
                const double alpha1_v0_fluc = alpha1_v0 + fluctuations[k][7];
                const double alpha2_v0_fluc = alpha2_v0 + fluctuations[k][8];

                double tp = (m_b + m_k) * (m_b + m_k);
                double tm = (m_b - m_k) * (m_b - m_k);
                double t0 = tp * (1 - sqrt(1 - tm / tp));
                double z = (sqrt(tp - q2) - sqrt(tp - t0)) / (sqrt(tp - q2) + sqrt(tp - t0));
                double z0 = (sqrt(tp) - sqrt(tp - t0)) / (sqrt(tp) + sqrt(tp - t0));

                double v0 = (1 / (1 - q2 / (mR_v0 * mR_v0))) * (alpha0_v0_fluc + alpha1_v0_fluc * (z - z0) + alpha2_v0_fluc * (z - z0) * (z - z0));
                double A1 = (1 / (1 - q2 / (mR_A1 * mR_A1))) * (alpha0_A1_fluc + alpha1_A1_fluc * (z - z0) + alpha2_A1_fluc * (z - z0) * (z - z0));
                double A12 = (1 / (1 - q2 / (mR_A12 * mR_A12))) * (alpha0_A12_fluc + alpha1_A12_fluc * (z - z0) + alpha2_A12_fluc * (z - z0) * (z - z0));
                double lambda = (tp - q2) * (tm - q2);
                double A2 = ((m_b + m_k) * (m_b + m_k) * (m_b * m_b - m_k * m_k - q2) * A1 - A12 * 16 * m_b * m_k * m_k * (m_b + m_k)) / lambda;

                double sB = q2 / (m_b * m_b);
                double m_k_tilda = m_k / m_b;
                double Lambda = 1 + std::pow(m_k_tilda, 4) + sB * sB - 2 * (m_k_tilda * m_k_tilda + sB + sB * m_k_tilda * m_k_tilda);

                double Amp_parallel = -2 * (std::sqrt(sB)) * (std::pow(Lambda, 1.0 / 4.0)) * (std::sqrt(2)) * (1 + m_k_tilda) * A1;
                double Amp_vertical = 2 * (std::sqrt(sB)) * (std::pow(Lambda, 1.0 / 4.0)) * (std::sqrt(2)) * (std::sqrt(Lambda)) * v0 / (1 + m_k_tilda);
                double Amp_0 = -1 * (std::sqrt(sB)) * (std::pow(Lambda, 1.0 / 4.0)) * (1.0 / m_k_tilda) * (1.0 / std::pow(sB, 0.5)) * ((1 - m_k_tilda * m_k_tilda - sB) * (1 + m_k_tilda) * A1 - Lambda * A2 / (1 + m_k_tilda));

                value[k] = (3.0 / 4.0) * (Amp_vertical * Amp_vertical + Amp_parallel * Amp_parallel) * (1 - costheta * costheta) + (3.0 / 2.0) * Amp_0 * Amp_0 * costheta * costheta;
                double total_weight = weight_var * Correction_pi0 * Correction_FEI * Correction_KID * (value[k] / value[0]);
                hist[k]->Fill(MVA_var, total_weight);
                Nevts[k] = Nevts[k] + total_weight;
            }
        }
        input_file->Close();
    }


    // calculate fluctuation of total BR
    for (int k = 0; k < 19; k++) {
        const int step = 100;
        double m_b_fix = -1;
        double m_k_fix = -1;
        if (charge == 0) {
            m_b_fix = 5.27965;
            m_k_fix = 0.89555;
        }
        else if (charge == 1 || charge == -1) {
            m_b_fix = 5.27934;
            m_k_fix = 0.89166;
        }
        const double q2_min = 0;
        const double q2_max = (m_b_fix - m_k_fix) * (m_b_fix - m_k_fix);
        const double delq2 = (q2_max - q2_min) / step; // max: (mb-mk)^2

        const double costheta_min = -1;
        const double costheta_max = 1;
        const double delcostheta = (costheta_max - costheta_min) / step;

        const double alpha0_A1_fluc = alpha0_A1 + fluctuations[k][0];
        const double alpha1_A1_fluc = alpha1_A1 + fluctuations[k][1];
        const double alpha2_A1_fluc = alpha2_A1 + fluctuations[k][2];
        const double alpha0_A12_fluc = alpha0_A12 + fluctuations[k][3];
        const double alpha1_A12_fluc = alpha1_A12 + fluctuations[k][4];
        const double alpha2_A12_fluc = alpha2_A12 + fluctuations[k][5];
        const double alpha0_v0_fluc = alpha0_v0 + fluctuations[k][6];
        const double alpha1_v0_fluc = alpha1_v0 + fluctuations[k][7];
        const double alpha2_v0_fluc = alpha2_v0 + fluctuations[k][8];

        for (int i = 0; i < step; i++) {
            const double q2_integral = q2_min + delq2 * i;

            for (int j = 0; j < step; j++) {
                const double costheta_integral = costheta_min + delcostheta * j;

                double tp = (m_b_fix + m_k_fix) * (m_b_fix + m_k_fix);
                double tm = (m_b_fix - m_k_fix) * (m_b_fix - m_k_fix);
                double t0 = tp * (1 - sqrt(1 - tm / tp));
                double z = (sqrt(tp - q2_integral) - sqrt(tp - t0)) / (sqrt(tp - q2_integral) + sqrt(tp - t0));
                double z0 = (sqrt(tp) - sqrt(tp - t0)) / (sqrt(tp) + sqrt(tp - t0));

                double v0 = (1 / (1 - q2_integral / (mR_v0 * mR_v0))) * (alpha0_v0_fluc + alpha1_v0_fluc * (z - z0) + alpha2_v0_fluc * (z - z0) * (z - z0));
                double A1 = (1 / (1 - q2_integral / (mR_A1 * mR_A1))) * (alpha0_A1_fluc + alpha1_A1_fluc * (z - z0) + alpha2_A1_fluc * (z - z0) * (z - z0));
                double A12 = (1 / (1 - q2_integral / (mR_A12 * mR_A12))) * (alpha0_A12_fluc + alpha1_A12_fluc * (z - z0) + alpha2_A12_fluc * (z - z0) * (z - z0));
                double lambda = (tp - q2_integral) * (tm - q2_integral);
                double A2 = ((m_b_fix + m_k_fix) * (m_b_fix + m_k_fix) * (m_b_fix * m_b_fix - m_k_fix * m_k_fix - q2_integral) * A1 - A12 * 16 * m_b_fix * m_k_fix * m_k_fix * (m_b_fix + m_k_fix)) / lambda;

                double sB = q2_integral / (m_b_fix * m_b_fix);
                double m_k_fix_tilda = m_k_fix / m_b_fix;
                double Lambda = 1 + std::pow(m_k_fix_tilda, 4) + sB * sB - 2 * (m_k_fix_tilda * m_k_fix_tilda + sB + sB * m_k_fix_tilda * m_k_fix_tilda);

                double Amp_parallel = -2 * (std::sqrt(sB)) * (std::pow(Lambda, 1.0 / 4.0)) * (std::sqrt(2)) * (1 + m_k_fix_tilda) * A1;
                double Amp_vertical = 2 * (std::sqrt(sB)) * (std::pow(Lambda, 1.0 / 4.0)) * (std::sqrt(2)) * (std::sqrt(Lambda)) * v0 / (1 + m_k_fix_tilda);
                double Amp_0 = -1 * (std::pow(Lambda, 1.0 / 4.0)) * (1.0 / m_k_fix_tilda) * ((1 - m_k_fix_tilda * m_k_fix_tilda - sB) * (1 + m_k_fix_tilda) * A1 - Lambda * A2 / (1 + m_k_fix_tilda));

                tot_value[k] = tot_value[k] + ((3.0 / 4.0) * (Amp_vertical * Amp_vertical + Amp_parallel * Amp_parallel) * (1 - costheta_integral * costheta_integral) + (3.0 / 2.0) * Amp_0 * Amp_0 * costheta_integral * costheta_integral) * delcostheta * (delq2 / (m_b_fix * m_b_fix));
            }
        }
        if (charge == 0) {
            Correction_factor_BR[k] = tot_value[k] / tot_value[0];
        }
        else if (charge == 1 || charge == -1) {
            Correction_factor_BR[k] = tot_value[k] / tot_value[0];
        }
    }


    return Nevt;
}

double GetFragmentationPDFs(const char* dirname, TH1D* hist, const char* type, DecayModeMC SelectedDecayMode, bool IsItUp, double weight_var = 1.0) { // get fragmentation uncertainty PDF with appropriate correction
    // Be careful! You should combine this function with `GetNominalPDFs` function!
    if (strcmp(type, "Bplus") == 0 && ( (SelectedDecayMode >= Xsu2KcPi0_MC && SelectedDecayMode <= Xsu2KcKcKcPi0_MC) || (SelectedDecayMode == other) ) ) {}
    else if (strcmp(type, "Bzero") == 0 && (SelectedDecayMode >= Xsd2KcPic_MC && SelectedDecayMode <= other) ) {}
    else if (strcmp(type, "Continuum") == 0) {
        printf("[ERROR] unexpected type name\n");
        printf("[ERROR] Continuum type cannot be selected for fragmentation pdf\n");
        exit(1);
    }
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_KaonID_correction = -1;
    double temp_nKaon_excep = -1;
    double temp_KaonID_rel_up = -1;
    double temp_KaonID_rel_dn = -1;
    int Decay[N_decay] = { 0 };

    double weight_not_Selected = GetFragmentationWeight(type, SelectedDecayMode, IsItUp);

    std::vector<string> names;
    load_files(dirname, &names);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");
        TTree* tree_Xs = (TTree*)input_file->Get("Xs");

        tree_upsilon->SetBranchAddress("MVA_Continuum", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_correction", &temp_KaonID_correction);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nKexcep", &temp_nKaon_excep);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_rel_uncer_up", &temp_KaonID_rel_up);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_rel_uncer_dn", &temp_KaonID_rel_dn);

        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKcharge_total__bc", &Decay[0]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &Decay[1]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &Decay[2]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCcomb__bc", &Decay[3]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch1__bc", &Decay[4]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch2__bc", &Decay[5]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch3__bc", &Decay[6]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch4__bc", &Decay[7]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch5__bc", &Decay[8]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch6__bc", &Decay[9]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch7__bc", &Decay[10]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch8__bc", &Decay[11]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch9__bc", &Decay[12]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch10__bc", &Decay[13]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch11__bc", &Decay[14]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch12__bc", &Decay[15]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch13__bc", &Decay[16]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch14__bc", &Decay[17]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch15__bc", &Decay[18]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKneutral_total__bc", &Decay[19]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &Decay[20]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &Decay[21]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCcomb__bc", &Decay[22]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch16__bc", &Decay[23]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch17__bc", &Decay[24]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch18__bc", &Decay[25]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch19__bc", &Decay[26]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch20__bc", &Decay[27]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch21__bc", &Decay[28]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch22__bc", &Decay[29]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch23__bc", &Decay[30]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch24__bc", &Decay[31]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch25__bc", &Decay[32]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch26__bc", &Decay[33]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch27__bc", &Decay[34]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch28__bc", &Decay[35]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch29__bc", &Decay[36]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch30__bc", &Decay[37]);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            tree_Xs->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

            double Correction_pi0 = std::pow(pi0_correction, Npi0);
            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = FEI_cal_Bc;
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = FEI_cal_B0;
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = temp_KaonID_correction * std::pow(-1, temp_nKaon_excep);

            double output_Decay[MAX_NUM_DECAYMODE_MC] = { 0.0 };
            DecayArrayToXsOutputDecay(Decay, output_Decay);
            double fragmentation_weight = 1.0;
            for (int k = 0; k < MAX_NUM_DECAYMODE_MC; k++) {
                if (output_Decay[k] > 0.5) {
                    if (k == static_cast<int>(SelectedDecayMode)) {
                        if (IsItUp) fragmentation_weight = fragmentation_weight * std::pow(1.5, output_Decay[k]);
                        else fragmentation_weight = fragmentation_weight * std::pow(0.5, output_Decay[k]);
                    }
                    else {
                        fragmentation_weight = fragmentation_weight * std::pow(weight_not_Selected, output_Decay[k]);
                    }
                }
            }

            double total_weight = weight_var * Correction_pi0 * Correction_FEI * Correction_KID * fragmentation_weight;

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

        printf("%s has %lf events (with correction)\n", dirname, Nevt);

    }
    return Nevt;
}

void Signal_yield_fit_BDT_Rarity_HistFactory()
{
    /* ====================================== */
    // Seting CDF module
    cdf.initbypath("/home/jwpark/storage/BKG_gbasf2/Izayoi_release5/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Knunu", Scale_Kplus_test);
    cdf.initbypath("/home/jwpark/storage/BKG_gbasf2/Izayoi_release5/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Kstarnunu", Scale_Kplusstar_test);
    cdf.initbypath("/home/jwpark/storage/BKG_gbasf2/Izayoi_release5/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Xsnunu", Scale_Xsu_nonresonant_test);
    cdf.initbypath("/home/jwpark/storage/BKG_gbasf2/Izayoi_release5/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0nunu", Scale_K0_test);
    cdf.initbypath("/home/jwpark/storage/BKG_gbasf2/Izayoi_release5/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0starnunu", Scale_K0star_test);
    cdf.initbypath("/home/jwpark/storage/BKG_gbasf2/Izayoi_release5/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02Xsnunu", Scale_Xsd_nonresonant_test);
    /* ====================================== */



    /* ====================================== */
    // Define PDFs for HistFactory
    // nominal PDFs
    TH1D* Signal_nominal = new TH1D("Signal_nominal", "Signal_nominal", RarityBins, 0.7, 1.0);
    TH1D* CHG_nominal = new TH1D("CHG_nominal", "CHG_nominal", RarityBins, 0.7, 1.0);
    TH1D* MIX_nominal = new TH1D("MIX_nominal", "MIX_nominal", RarityBins, 0.7, 1.0);
    TH1D* UUBAR_nominal = new TH1D("UUBAR_nominal", "UUBAR_nominal", RarityBins, 0.7, 1.0);
    TH1D* DDBAR_nominal = new TH1D("DDBAR_nominal", "DDBAR_nominal", RarityBins, 0.7, 1.0);
    TH1D* SSBAR_nominal = new TH1D("SSBAR_nominal", "SSBAR_nominal", RarityBins, 0.7, 1.0);
    TH1D* CHARM_nominal = new TH1D("CHARM_nominal", "CHARM_nominal", RarityBins, 0.7, 1.0);

    // FEI uncertainty, BKGs will be cared by overall syst
    TH1D* Signal_FEI_charged_p = new TH1D("Signal_FEI_charged_p", "Signal_FEI_charged_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_FEI_charged_m = new TH1D("Signal_FEI_charged_m", "Signal_FEI_charged_m", RarityBins, 0.7, 1.0);
    TH1D* Signal_FEI_neutral_p = new TH1D("Signal_FEI_neutral_p", "Signal_FEI_neutral_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_FEI_neutral_m = new TH1D("Signal_FEI_neutral_m", "Signal_FEI_neutral_m", RarityBins, 0.7, 1.0);

    // pi0 uncertainty
    TH1D* Signal_pi0_p = new TH1D("Signal_pi0_p", "Signal_pi0_p", RarityBins, 0.7, 1.0);
    TH1D* CHG_pi0_p = new TH1D("CHG_pi0_p", "CHG_pi0_p", RarityBins, 0.7, 1.0);
    TH1D* MIX_pi0_p = new TH1D("MIX_pi0_p", "MIX_pi0_p", RarityBins, 0.7, 1.0);
    TH1D* UUBAR_pi0_p = new TH1D("UUBAR_pi0_p", "UUBAR_pi0_p", RarityBins, 0.7, 1.0);
    TH1D* DDBAR_pi0_p = new TH1D("DDBAR_pi0_p", "DDBAR_pi0_p", RarityBins, 0.7, 1.0);
    TH1D* SSBAR_pi0_p = new TH1D("SSBAR_pi0_p", "SSBAR_pi0_p", RarityBins, 0.7, 1.0);
    TH1D* CHARM_pi0_p = new TH1D("CHARM_pi0_p", "CHARM_pi0_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_pi0_m = new TH1D("Signal_pi0_m", "Signal_pi0_m", RarityBins, 0.7, 1.0);
    TH1D* CHG_pi0_m = new TH1D("CHG_pi0_m", "CHG_pi0_m", RarityBins, 0.7, 1.0);
    TH1D* MIX_pi0_m = new TH1D("MIX_pi0_m", "MIX_pi0_m", RarityBins, 0.7, 1.0);
    TH1D* UUBAR_pi0_m = new TH1D("UUBAR_pi0_m", "UUBAR_pi0_m", RarityBins, 0.7, 1.0);
    TH1D* DDBAR_pi0_m = new TH1D("DDBAR_pi0_m", "DDBAR_pi0_m", RarityBins, 0.7, 1.0);
    TH1D* SSBAR_pi0_m = new TH1D("SSBAR_pi0_m", "SSBAR_pi0_m", RarityBins, 0.7, 1.0);
    TH1D* CHARM_pi0_m = new TH1D("CHARM_pi0_m", "CHARM_pi0_m", RarityBins, 0.7, 1.0);

    // track uncertainty
    TH1D* Signal_track_p = new TH1D("Signal_track_p", "Signal_track_p", RarityBins, 0.7, 1.0);
    TH1D* CHG_track_p = new TH1D("CHG_track_p", "CHG_track_p", RarityBins, 0.7, 1.0);
    TH1D* MIX_track_p = new TH1D("MIX_track_p", "MIX_track_p", RarityBins, 0.7, 1.0);
    TH1D* UUBAR_track_p = new TH1D("UUBAR_track_p", "UUBAR_track_p", RarityBins, 0.7, 1.0);
    TH1D* DDBAR_track_p = new TH1D("DDBAR_track_p", "DDBAR_track_p", RarityBins, 0.7, 1.0);
    TH1D* SSBAR_track_p = new TH1D("SSBAR_track_p", "SSBAR_track_p", RarityBins, 0.7, 1.0);
    TH1D* CHARM_track_p = new TH1D("CHARM_track_p", "CHARM_track_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_track_m = new TH1D("Signal_track_m", "Signal_track_m", RarityBins, 0.7, 1.0);
    TH1D* CHG_track_m = new TH1D("CHG_track_m", "CHG_track_m", RarityBins, 0.7, 1.0);
    TH1D* MIX_track_m = new TH1D("MIX_track_m", "MIX_track_m", RarityBins, 0.7, 1.0);
    TH1D* UUBAR_track_m = new TH1D("UUBAR_track_m", "UUBAR_track_m", RarityBins, 0.7, 1.0);
    TH1D* DDBAR_track_m = new TH1D("DDBAR_track_m", "DDBAR_track_m", RarityBins, 0.7, 1.0);
    TH1D* SSBAR_track_m = new TH1D("SSBAR_track_m", "SSBAR_track_m", RarityBins, 0.7, 1.0);
    TH1D* CHARM_track_m = new TH1D("CHARM_track_m", "CHARM_track_m", RarityBins, 0.7, 1.0);

    // KS0 uncertainty
    TH1D* Signal_KS0_p = new TH1D("Signal_KS0_p", "Signal_KS0_p", RarityBins, 0.7, 1.0);
    TH1D* CHG_KS0_p = new TH1D("CHG_KS0_p", "CHG_KS0_p", RarityBins, 0.7, 1.0);
    TH1D* MIX_KS0_p = new TH1D("MIX_KS0_p", "MIX_KS0_p", RarityBins, 0.7, 1.0);
    TH1D* UUBAR_KS0_p = new TH1D("UUBAR_KS0_p", "UUBAR_KS0_p", RarityBins, 0.7, 1.0);
    TH1D* DDBAR_KS0_p = new TH1D("DDBAR_KS0_p", "DDBAR_KS0_p", RarityBins, 0.7, 1.0);
    TH1D* SSBAR_KS0_p = new TH1D("SSBAR_KS0_p", "SSBAR_KS0_p", RarityBins, 0.7, 1.0);
    TH1D* CHARM_KS0_p = new TH1D("CHARM_KS0_p", "CHARM_KS0_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_KS0_m = new TH1D("Signal_KS0_m", "Signal_KS0_m", RarityBins, 0.7, 1.0);
    TH1D* CHG_KS0_m = new TH1D("CHG_KS0_m", "CHG_KS0_m", RarityBins, 0.7, 1.0);
    TH1D* MIX_KS0_m = new TH1D("MIX_KS0_m", "MIX_KS0_m", RarityBins, 0.7, 1.0);
    TH1D* UUBAR_KS0_m = new TH1D("UUBAR_KS0_m", "UUBAR_KS0_m", RarityBins, 0.7, 1.0);
    TH1D* DDBAR_KS0_m = new TH1D("DDBAR_KS0_m", "DDBAR_KS0_m", RarityBins, 0.7, 1.0);
    TH1D* SSBAR_KS0_m = new TH1D("SSBAR_KS0_m", "SSBAR_KS0_m", RarityBins, 0.7, 1.0);
    TH1D* CHARM_KS0_m = new TH1D("CHARM_KS0_m", "CHARM_KS0_m", RarityBins, 0.7, 1.0);

    // Kaon PID uncertainty
    TH1D* Signal_KID_p = new TH1D("Signal_KID_p", "Signal_KID_p", RarityBins, 0.7, 1.0);
    TH1D* CHG_KID_p = new TH1D("CHG_KID_p", "CHG_KID_p", RarityBins, 0.7, 1.0);
    TH1D* MIX_KID_p = new TH1D("MIX_KID_p", "MIX_KID_p", RarityBins, 0.7, 1.0);
    TH1D* UUBAR_KID_p = new TH1D("UUBAR_KID_p", "UUBAR_KID_p", RarityBins, 0.7, 1.0);
    TH1D* DDBAR_KID_p = new TH1D("DDBAR_KID_p", "DDBAR_KID_p", RarityBins, 0.7, 1.0);
    TH1D* SSBAR_KID_p = new TH1D("SSBAR_KID_p", "SSBAR_KID_p", RarityBins, 0.7, 1.0);
    TH1D* CHARM_KID_p = new TH1D("CHARM_KID_p", "CHARM_KID_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_KID_m = new TH1D("Signal_KID_m", "Signal_KID_m", RarityBins, 0.7, 1.0);
    TH1D* CHG_KID_m = new TH1D("CHG_KID_m", "CHG_KID_m", RarityBins, 0.7, 1.0);
    TH1D* MIX_KID_m = new TH1D("MIX_KID_m", "MIX_KID_m", RarityBins, 0.7, 1.0);
    TH1D* UUBAR_KID_m = new TH1D("UUBAR_KID_m", "UUBAR_KID_m", RarityBins, 0.7, 1.0);
    TH1D* DDBAR_KID_m = new TH1D("DDBAR_KID_m", "DDBAR_KID_m", RarityBins, 0.7, 1.0);
    TH1D* SSBAR_KID_m = new TH1D("SSBAR_KID_m", "SSBAR_KID_m", RarityBins, 0.7, 1.0);
    TH1D* CHARM_KID_m = new TH1D("CHARM_KID_m", "CHARM_KID_m", RarityBins, 0.7, 1.0);

    // K nu nubar form factor
    TH1D* Signal_Kff1_p = new TH1D("Signal_Kff1_p", "Signal_Kff1_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kff1_m = new TH1D("Signal_Kff1_m", "Signal_Kff1_m", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kff2_p = new TH1D("Signal_Kff2_p", "Signal_Kff2_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kff2_m = new TH1D("Signal_Kff2_m", "Signal_Kff2_m", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kff3_p = new TH1D("Signal_Kff3_p", "Signal_Kff3_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kff3_m = new TH1D("Signal_Kff3_m", "Signal_Kff3_m", RarityBins, 0.7, 1.0);

    // Kstar nu nubar form factor
    TH1D* Signal_Kstarff1_p = new TH1D("Signal_Kstarff1_p", "Signal_Kstarff1_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kstarff1_m = new TH1D("Signal_Kstarff1_m", "Signal_Kstarff1_m", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kstarff2_p = new TH1D("Signal_Kstarff2_p", "Signal_Kstarff2_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kstarff2_m = new TH1D("Signal_Kstarff2_m", "Signal_Kstarff2_m", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kstarff3_p = new TH1D("Signal_Kstarff3_p", "Signal_Kstarff3_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kstarff3_m = new TH1D("Signal_Kstarff3_m", "Signal_Kstarff3_m", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kstarff4_p = new TH1D("Signal_Kstarff4_p", "Signal_Kstarff4_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kstarff4_m = new TH1D("Signal_Kstarff4_m", "Signal_Kstarff4_m", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kstarff5_p = new TH1D("Signal_Kstarff5_p", "Signal_Kstarff5_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kstarff5_m = new TH1D("Signal_Kstarff5_m", "Signal_Kstarff5_m", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kstarff6_p = new TH1D("Signal_Kstarff6_p", "Signal_Kstarff6_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kstarff6_m = new TH1D("Signal_Kstarff6_m", "Signal_Kstarff6_m", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kstarff7_p = new TH1D("Signal_Kstarff7_p", "Signal_Kstarff7_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kstarff7_m = new TH1D("Signal_Kstarff7_m", "Signal_Kstarff7_m", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kstarff8_p = new TH1D("Signal_Kstarff8_p", "Signal_Kstarff8_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kstarff8_m = new TH1D("Signal_Kstarff8_m", "Signal_Kstarff8_m", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kstarff9_p = new TH1D("Signal_Kstarff9_p", "Signal_Kstarff9_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kstarff9_m = new TH1D("Signal_Kstarff9_m", "Signal_Kstarff9_m", RarityBins, 0.7, 1.0);

    // faction
    TH1D* Signal_Kfrac_p = new TH1D("Signal_Kfrac_p", "Signal_Kfrac_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kfrac_m = new TH1D("Signal_Kfrac_m", "Signal_Kfrac_m", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kstarfrac_p = new TH1D("Signal_Kstarfrac_p", "Signal_Kstarfrac_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kstarfrac_m = new TH1D("Signal_Kstarfrac_m", "Signal_Kstarfrac_m", RarityBins, 0.7, 1.0);

    // fragmentation
    std::vector<std::string> name_Signal_Fragmentation = {
        "Signal_Xsu_frag_decay1_p", "Signal_Xsu_frag_decay2_p", "Signal_Xsu_frag_decay3_p",
        "Signal_Xsu_frag_decay4_p", "Signal_Xsu_frag_decay5_p", "Signal_Xsu_frag_decay6_p",
        "Signal_Xsu_frag_decay7_p", "Signal_Xsu_frag_decay8_p", "Signal_Xsu_frag_decay9_p",
        "Signal_Xsu_frag_decay10_p", "Signal_Xsu_frag_decay11_p", "Signal_Xsu_frag_decay12_p",
        "Signal_Xsu_frag_decay13_p", "Signal_Xsu_frag_decay14_p", "Signal_Xsu_frag_decay30_p",
        "Signal_Xsd_frag_decay16_p", "Signal_Xsd_frag_decay17_p", "Signal_Xsd_frag_decay18_p",
        "Signal_Xsd_frag_decay19_p", "Signal_Xsd_frag_decay20_p", "Signal_Xsd_frag_decay21_p",
        "Signal_Xsd_frag_decay22_p", "Signal_Xsd_frag_decay23_p", "Signal_Xsd_frag_decay24_p",
        "Signal_Xsd_frag_decay25_p", "Signal_Xsd_frag_decay26_p", "Signal_Xsd_frag_decay27_p",
        "Signal_Xsd_frag_decay28_p", "Signal_Xsd_frag_decay29_p", "Signal_Xsd_frag_decay30_p",
        "Signal_Xsu_frag_decay1_m", "Signal_Xsu_frag_decay2_m", "Signal_Xsu_frag_decay3_m",
        "Signal_Xsu_frag_decay4_m", "Signal_Xsu_frag_decay5_m", "Signal_Xsu_frag_decay6_m",
        "Signal_Xsu_frag_decay7_m", "Signal_Xsu_frag_decay8_m", "Signal_Xsu_frag_decay9_m",
        "Signal_Xsu_frag_decay10_m", "Signal_Xsu_frag_decay11_m", "Signal_Xsu_frag_decay12_m",
        "Signal_Xsu_frag_decay13_m", "Signal_Xsu_frag_decay14_m", "Signal_Xsu_frag_decay30_m",
        "Signal_Xsd_frag_decay16_m", "Signal_Xsd_frag_decay17_m", "Signal_Xsd_frag_decay18_m",
        "Signal_Xsd_frag_decay19_m", "Signal_Xsd_frag_decay20_m", "Signal_Xsd_frag_decay21_m",
        "Signal_Xsd_frag_decay22_m", "Signal_Xsd_frag_decay23_m", "Signal_Xsd_frag_decay24_m",
        "Signal_Xsd_frag_decay25_m", "Signal_Xsd_frag_decay26_m", "Signal_Xsd_frag_decay27_m",
        "Signal_Xsd_frag_decay28_m", "Signal_Xsd_frag_decay29_m", "Signal_Xsd_frag_decay30_m",
    };
    std::vector<TH1D*> Signal_Fragmentaions;
    for (unsigned int i = 0; i < name_Signal_Fragmentation.size(); i++) {
        TH1D* temp = new TH1D(name_Signal_Fragmentation.at(i).c_str(), name_Signal_Fragmentation.at(i).c_str(), RarityBins, 0.7, 1.0);
        Signal_Fragmentaions.push_back(temp);
    }

    // pf
    TH1D* Signal_pf_p = new TH1D("Signal_pf_p", "Signal_pf_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_pf_m = new TH1D("Signal_pf_m", "Signal_pf_m", RarityBins, 0.7, 1.0);

    // mb
    TH1D* Signal_mb_p = new TH1D("Signal_mb_p", "Signal_mb_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_mb_m = new TH1D("Signal_mb_m", "Signal_mb_m", RarityBins, 0.7, 1.0);

    // transition
    TH1D* Signal_transition_p = new TH1D("Signal_transition_p", "Signal_transition_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_transition_m = new TH1D("Signal_transition_m", "Signal_transition_m", RarityBins, 0.7, 1.0);

    // mK*
    TH1D* Signal_mKstar_p = new TH1D("Signal_mKstar_p", "Signal_mKstar_p", RarityBins, 0.7, 1.0);
    TH1D* Signal_mKstar_m = new TH1D("Signal_mKstar_m", "Signal_mKstar_m", RarityBins, 0.7, 1.0);
    /* ====================================== */



    /* ====================================== */
    // define path for Ntuple
    const char* MC_dirname_Knunu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Knunu";
    const char* MC_dirname_Kstarnunu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Kstarnunu";
    const char* MC_dirname_Xsununu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Xsnunu";
    const char* MC_dirname_K0nunu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0nunu";
    const char* MC_dirname_K0starnunu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0starnunu";
    const char* MC_dirname_Xsdnunu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02Xsnunu";

    const char* MC_dirname_CHG = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5/CHG_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_MIX = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5/MIX_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_UUBAR = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5/UUBAR_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_DDBAR = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5/DDBAR_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_SSBAR = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5/SSBAR_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_CHARM = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5/CHARM_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";

    // for signal modeling
    const char* MC_dirname_Xsununu_Hpf = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5_syst/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/Xsu_Hpf";
    const char* MC_dirname_Xsdnunu_Hpf = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5_syst/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/Xsd_Hpf";
    const char* MC_dirname_Xsununu_Lpf = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5_syst/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/Xsu_Lpf";
    const char* MC_dirname_Xsdnunu_Lpf = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5_syst/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/Xsd_Lpf";
    const char* MC_dirname_Xsununu_Hmb = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5_syst/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/Xsu_Hmb";
    const char* MC_dirname_Xsdnunu_Hmb = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5_syst/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/Xsd_Hmb";
    const char* MC_dirname_Xsununu_Lmb = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5_syst/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/Xsu_Lmb";
    const char* MC_dirname_Xsdnunu_Lmb = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5_syst/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/Xsd_Lmb";
    const char* MC_dirname_Xsununu_Htransition = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5_syst/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/Xsu_Htransition";
    const char* MC_dirname_Xsdnunu_Htransition = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5_syst/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/Xsd_Htransition";
    const char* MC_dirname_Xsununu_Ltransition = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5_syst/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/Xsu_Ltransition";
    const char* MC_dirname_Xsdnunu_Ltransition = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5_syst/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/Xsd_Ltransition";

    // for fixed mKstar
    const char* MC_dirname_Kstarnunu_mKstarfixed = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5_syst/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/B2Kstarnunu";
    const char* MC_dirname_K0starnunu_mKstarfixed = "/home/jwpark/storage/BKG_gbasf2/Izayoi_release5_syst/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/B02Kstar0nunu";
    /* ====================================== */



    /* ====================================== */
    // Get PDFs
    // get nominal pdfs
    GetNominalPDFs(MC_dirname_Knunu, Signal_nominal, "Bplus", Scale_Kplus_test);
    GetNominalPDFs(MC_dirname_Kstarnunu, Signal_nominal, "Bplus", Scale_Kplusstar_test);
    GetNominalPDFs(MC_dirname_Xsununu, Signal_nominal, "Bplus", Scale_Xsu_nonresonant_test);
    GetNominalPDFs(MC_dirname_K0nunu, Signal_nominal, "Bzero", Scale_K0_test);
    GetNominalPDFs(MC_dirname_K0starnunu, Signal_nominal, "Bzero", Scale_K0star_test);
    GetNominalPDFs(MC_dirname_Xsdnunu, Signal_nominal, "Bzero", Scale_Xsd_nonresonant_test);

    GetNominalPDFs(MC_dirname_CHG, CHG_nominal, "Bplus", Scale_CHG_test);
    GetNominalPDFs(MC_dirname_MIX, MIX_nominal, "Bzero", Scale_MIX_test);
    GetNominalPDFs(MC_dirname_UUBAR, UUBAR_nominal, "Continuum", Scale_UUBAR_test);
    GetNominalPDFs(MC_dirname_DDBAR, DDBAR_nominal, "Continuum", Scale_DDBAR_test);
    GetNominalPDFs(MC_dirname_SSBAR, SSBAR_nominal, "Continuum", Scale_SSBAR_test);
    GetNominalPDFs(MC_dirname_CHARM, CHARM_nominal, "Continuum", Scale_CHARM_test);

    // get FEI uncertainty pdfs
    GetFEIPDFs(MC_dirname_Knunu, Signal_FEI_charged_p, "Bplus", true, Scale_Kplus_test);
    GetFEIPDFs(MC_dirname_Kstarnunu, Signal_FEI_charged_p, "Bplus", true, Scale_Kplusstar_test);
    GetFEIPDFs(MC_dirname_Xsununu, Signal_FEI_charged_p, "Bplus", true, Scale_Xsu_nonresonant_test);
    GetNominalPDFs(MC_dirname_K0nunu, Signal_FEI_charged_p, "Bzero", Scale_K0_test);
    GetNominalPDFs(MC_dirname_K0starnunu, Signal_FEI_charged_p, "Bzero", Scale_K0star_test);
    GetNominalPDFs(MC_dirname_Xsdnunu, Signal_FEI_charged_p, "Bzero", Scale_Xsd_nonresonant_test);

    GetFEIPDFs(MC_dirname_Knunu, Signal_FEI_charged_m, "Bplus", false, Scale_Kplus_test);
    GetFEIPDFs(MC_dirname_Kstarnunu, Signal_FEI_charged_m, "Bplus", false, Scale_Kplusstar_test);
    GetFEIPDFs(MC_dirname_Xsununu, Signal_FEI_charged_m, "Bplus", false, Scale_Xsu_nonresonant_test);
    GetNominalPDFs(MC_dirname_K0nunu, Signal_FEI_charged_m, "Bzero", Scale_K0_test);
    GetNominalPDFs(MC_dirname_K0starnunu, Signal_FEI_charged_m, "Bzero", Scale_K0star_test);
    GetNominalPDFs(MC_dirname_Xsdnunu, Signal_FEI_charged_m, "Bzero", Scale_Xsd_nonresonant_test);

    GetNominalPDFs(MC_dirname_Knunu, Signal_FEI_neutral_p, "Bplus", Scale_Kplus_test);
    GetNominalPDFs(MC_dirname_Kstarnunu, Signal_FEI_neutral_p, "Bplus", Scale_Kplusstar_test);
    GetNominalPDFs(MC_dirname_Xsununu, Signal_FEI_neutral_p, "Bplus", Scale_Xsu_nonresonant_test);
    GetFEIPDFs(MC_dirname_K0nunu, Signal_FEI_neutral_p, "Bzero", true, Scale_K0_test);
    GetFEIPDFs(MC_dirname_K0starnunu, Signal_FEI_neutral_p, "Bzero", true, Scale_K0star_test);
    GetFEIPDFs(MC_dirname_Xsdnunu, Signal_FEI_neutral_p, "Bzero", true, Scale_Xsd_nonresonant_test);

    GetNominalPDFs(MC_dirname_Knunu, Signal_FEI_neutral_m, "Bplus", Scale_Kplus_test);
    GetNominalPDFs(MC_dirname_Kstarnunu, Signal_FEI_neutral_m, "Bplus", Scale_Kplusstar_test);
    GetNominalPDFs(MC_dirname_Xsununu, Signal_FEI_neutral_m, "Bplus", Scale_Xsu_nonresonant_test);
    GetFEIPDFs(MC_dirname_K0nunu, Signal_FEI_neutral_m, "Bzero", false, Scale_K0_test);
    GetFEIPDFs(MC_dirname_K0starnunu, Signal_FEI_neutral_m, "Bzero", false, Scale_K0star_test);
    GetFEIPDFs(MC_dirname_Xsdnunu, Signal_FEI_neutral_m, "Bzero", false, Scale_Xsd_nonresonant_test);

    // get pi0 uncertainty pdfs
    GetPi0PDFs(MC_dirname_Knunu, Signal_pi0_p, "Bplus", true, Scale_Kplus_test);
    GetPi0PDFs(MC_dirname_Kstarnunu, Signal_pi0_p, "Bplus", true, Scale_Kplusstar_test);
    GetPi0PDFs(MC_dirname_Xsununu, Signal_pi0_p, "Bplus", true, Scale_Xsu_nonresonant_test);
    GetPi0PDFs(MC_dirname_K0nunu, Signal_pi0_p, "Bzero", true, Scale_K0_test);
    GetPi0PDFs(MC_dirname_K0starnunu, Signal_pi0_p, "Bzero", true, Scale_K0star_test);
    GetPi0PDFs(MC_dirname_Xsdnunu, Signal_pi0_p, "Bzero", true, Scale_Xsd_nonresonant_test);

    GetPi0PDFs(MC_dirname_CHG, CHG_pi0_p, "Bplus", true, Scale_CHG_test);
    GetPi0PDFs(MC_dirname_MIX, MIX_pi0_p, "Bzero", true, Scale_MIX_test);
    GetPi0PDFs(MC_dirname_UUBAR, UUBAR_pi0_p, "Continuum", true, Scale_UUBAR_test);
    GetPi0PDFs(MC_dirname_DDBAR, DDBAR_pi0_p, "Continuum", true, Scale_DDBAR_test);
    GetPi0PDFs(MC_dirname_SSBAR, SSBAR_pi0_p, "Continuum", true, Scale_SSBAR_test);
    GetPi0PDFs(MC_dirname_CHARM, CHARM_pi0_p, "Continuum", true, Scale_CHARM_test);

    GetPi0PDFs(MC_dirname_Knunu, Signal_pi0_m, "Bplus", false, Scale_Kplus_test);
    GetPi0PDFs(MC_dirname_Kstarnunu, Signal_pi0_m, "Bplus", false, Scale_Kplusstar_test);
    GetPi0PDFs(MC_dirname_Xsununu, Signal_pi0_m, "Bplus", false, Scale_Xsu_nonresonant_test);
    GetPi0PDFs(MC_dirname_K0nunu, Signal_pi0_m, "Bzero", false, Scale_K0_test);
    GetPi0PDFs(MC_dirname_K0starnunu, Signal_pi0_m, "Bzero", false, Scale_K0star_test);
    GetPi0PDFs(MC_dirname_Xsdnunu, Signal_pi0_m, "Bzero", false, Scale_Xsd_nonresonant_test);

    GetPi0PDFs(MC_dirname_CHG, CHG_pi0_m, "Bplus", false, Scale_CHG_test);
    GetPi0PDFs(MC_dirname_MIX, MIX_pi0_m, "Bzero", false, Scale_MIX_test);
    GetPi0PDFs(MC_dirname_UUBAR, UUBAR_pi0_m, "Continuum", false, Scale_UUBAR_test);
    GetPi0PDFs(MC_dirname_DDBAR, DDBAR_pi0_m, "Continuum", false, Scale_DDBAR_test);
    GetPi0PDFs(MC_dirname_SSBAR, SSBAR_pi0_m, "Continuum", false, Scale_SSBAR_test);
    GetPi0PDFs(MC_dirname_CHARM, CHARM_pi0_m, "Continuum", false, Scale_CHARM_test);

    // get track uncertainty pdfs
    GetTrackPDFs(MC_dirname_Knunu, Signal_track_p, "Bplus", true, Scale_Kplus_test);
    GetTrackPDFs(MC_dirname_Kstarnunu, Signal_track_p, "Bplus", true, Scale_Kplusstar_test);
    GetTrackPDFs(MC_dirname_Xsununu, Signal_track_p, "Bplus", true, Scale_Xsu_nonresonant_test);
    GetTrackPDFs(MC_dirname_K0nunu, Signal_track_p, "Bzero", true, Scale_K0_test);
    GetTrackPDFs(MC_dirname_K0starnunu, Signal_track_p, "Bzero", true, Scale_K0star_test);
    GetTrackPDFs(MC_dirname_Xsdnunu, Signal_track_p, "Bzero", true, Scale_Xsd_nonresonant_test);

    GetTrackPDFs(MC_dirname_CHG, CHG_track_p, "Bplus", true, Scale_CHG_test);
    GetTrackPDFs(MC_dirname_MIX, MIX_track_p, "Bzero", true, Scale_MIX_test);
    GetTrackPDFs(MC_dirname_UUBAR, UUBAR_track_p, "Continuum", true, Scale_UUBAR_test);
    GetTrackPDFs(MC_dirname_DDBAR, DDBAR_track_p, "Continuum", true, Scale_DDBAR_test);
    GetTrackPDFs(MC_dirname_SSBAR, SSBAR_track_p, "Continuum", true, Scale_SSBAR_test);
    GetTrackPDFs(MC_dirname_CHARM, CHARM_track_p, "Continuum", true, Scale_CHARM_test);

    GetTrackPDFs(MC_dirname_Knunu, Signal_track_m, "Bplus", false, Scale_Kplus_test);
    GetTrackPDFs(MC_dirname_Kstarnunu, Signal_track_m, "Bplus", false, Scale_Kplusstar_test);
    GetTrackPDFs(MC_dirname_Xsununu, Signal_track_m, "Bplus", false, Scale_Xsu_nonresonant_test);
    GetTrackPDFs(MC_dirname_K0nunu, Signal_track_m, "Bzero", false, Scale_K0_test);
    GetTrackPDFs(MC_dirname_K0starnunu, Signal_track_m, "Bzero", false, Scale_K0star_test);
    GetTrackPDFs(MC_dirname_Xsdnunu, Signal_track_m, "Bzero", false, Scale_Xsd_nonresonant_test);

    GetTrackPDFs(MC_dirname_CHG, CHG_track_m, "Bplus", false, Scale_CHG_test);
    GetTrackPDFs(MC_dirname_MIX, MIX_track_m, "Bzero", false, Scale_MIX_test);
    GetTrackPDFs(MC_dirname_UUBAR, UUBAR_track_m, "Continuum", false, Scale_UUBAR_test);
    GetTrackPDFs(MC_dirname_DDBAR, DDBAR_track_m, "Continuum", false, Scale_DDBAR_test);
    GetTrackPDFs(MC_dirname_SSBAR, SSBAR_track_m, "Continuum", false, Scale_SSBAR_test);
    GetTrackPDFs(MC_dirname_CHARM, CHARM_track_m, "Continuum", false, Scale_CHARM_test);

    // get KS0 uncertainty pdfs
    GetKS0PDFs(MC_dirname_Knunu, Signal_KS0_p, "Bplus", true, Scale_Kplus_test);
    GetKS0PDFs(MC_dirname_Kstarnunu, Signal_KS0_p, "Bplus", true, Scale_Kplusstar_test);
    GetKS0PDFs(MC_dirname_Xsununu, Signal_KS0_p, "Bplus", true, Scale_Xsu_nonresonant_test);
    GetKS0PDFs(MC_dirname_K0nunu, Signal_KS0_p, "Bzero", true, Scale_K0_test);
    GetKS0PDFs(MC_dirname_K0starnunu, Signal_KS0_p, "Bzero", true, Scale_K0star_test);
    GetKS0PDFs(MC_dirname_Xsdnunu, Signal_KS0_p, "Bzero", true, Scale_Xsd_nonresonant_test);

    GetKS0PDFs(MC_dirname_CHG, CHG_KS0_p, "Bplus", true, Scale_CHG_test);
    GetKS0PDFs(MC_dirname_MIX, MIX_KS0_p, "Bzero", true, Scale_MIX_test);
    GetKS0PDFs(MC_dirname_UUBAR, UUBAR_KS0_p, "Continuum", true, Scale_UUBAR_test);
    GetKS0PDFs(MC_dirname_DDBAR, DDBAR_KS0_p, "Continuum", true, Scale_DDBAR_test);
    GetKS0PDFs(MC_dirname_SSBAR, SSBAR_KS0_p, "Continuum", true, Scale_SSBAR_test);
    GetKS0PDFs(MC_dirname_CHARM, CHARM_KS0_p, "Continuum", true, Scale_CHARM_test);

    GetKS0PDFs(MC_dirname_Knunu, Signal_KS0_m, "Bplus", false, Scale_Kplus_test);
    GetKS0PDFs(MC_dirname_Kstarnunu, Signal_KS0_m, "Bplus", false, Scale_Kplusstar_test);
    GetKS0PDFs(MC_dirname_Xsununu, Signal_KS0_m, "Bplus", false, Scale_Xsu_nonresonant_test);
    GetKS0PDFs(MC_dirname_K0nunu, Signal_KS0_m, "Bzero", false, Scale_K0_test);
    GetKS0PDFs(MC_dirname_K0starnunu, Signal_KS0_m, "Bzero", false, Scale_K0star_test);
    GetKS0PDFs(MC_dirname_Xsdnunu, Signal_KS0_m, "Bzero", false, Scale_Xsd_nonresonant_test);

    GetKS0PDFs(MC_dirname_CHG, CHG_KS0_m, "Bplus", false, Scale_CHG_test);
    GetKS0PDFs(MC_dirname_MIX, MIX_KS0_m, "Bzero", false, Scale_MIX_test);
    GetKS0PDFs(MC_dirname_UUBAR, UUBAR_KS0_m, "Continuum", false, Scale_UUBAR_test);
    GetKS0PDFs(MC_dirname_DDBAR, DDBAR_KS0_m, "Continuum", false, Scale_DDBAR_test);
    GetKS0PDFs(MC_dirname_SSBAR, SSBAR_KS0_m, "Continuum", false, Scale_SSBAR_test);
    GetKS0PDFs(MC_dirname_CHARM, CHARM_KS0_m, "Continuum", false, Scale_CHARM_test);

    // get KID uncertainty pdfs
    GetKIDPDFs(MC_dirname_Knunu, Signal_KID_p, "Bplus", true, Scale_Kplus_test);
    GetKIDPDFs(MC_dirname_Kstarnunu, Signal_KID_p, "Bplus", true, Scale_Kplusstar_test);
    GetKIDPDFs(MC_dirname_Xsununu, Signal_KID_p, "Bplus", true, Scale_Xsu_nonresonant_test);
    GetKIDPDFs(MC_dirname_K0nunu, Signal_KID_p, "Bzero", true, Scale_K0_test);
    GetKIDPDFs(MC_dirname_K0starnunu, Signal_KID_p, "Bzero", true, Scale_K0star_test);
    GetKIDPDFs(MC_dirname_Xsdnunu, Signal_KID_p, "Bzero", true, Scale_Xsd_nonresonant_test);

    GetKIDPDFs(MC_dirname_CHG, CHG_KID_p, "Bplus", true, Scale_CHG_test);
    GetKIDPDFs(MC_dirname_MIX, MIX_KID_p, "Bzero", true, Scale_MIX_test);
    GetKIDPDFs(MC_dirname_UUBAR, UUBAR_KID_p, "Continuum", true, Scale_UUBAR_test);
    GetKIDPDFs(MC_dirname_DDBAR, DDBAR_KID_p, "Continuum", true, Scale_DDBAR_test);
    GetKIDPDFs(MC_dirname_SSBAR, SSBAR_KID_p, "Continuum", true, Scale_SSBAR_test);
    GetKIDPDFs(MC_dirname_CHARM, CHARM_KID_p, "Continuum", true, Scale_CHARM_test);

    GetKIDPDFs(MC_dirname_Knunu, Signal_KID_m, "Bplus", false, Scale_Kplus_test);
    GetKIDPDFs(MC_dirname_Kstarnunu, Signal_KID_m, "Bplus", false, Scale_Kplusstar_test);
    GetKIDPDFs(MC_dirname_Xsununu, Signal_KID_m, "Bplus", false, Scale_Xsu_nonresonant_test);
    GetKIDPDFs(MC_dirname_K0nunu, Signal_KID_m, "Bzero", false, Scale_K0_test);
    GetKIDPDFs(MC_dirname_K0starnunu, Signal_KID_m, "Bzero", false, Scale_K0star_test);
    GetKIDPDFs(MC_dirname_Xsdnunu, Signal_KID_m, "Bzero", false, Scale_Xsd_nonresonant_test);

    GetKIDPDFs(MC_dirname_CHG, CHG_KID_m, "Bplus", false, Scale_CHG_test);
    GetKIDPDFs(MC_dirname_MIX, MIX_KID_m, "Bzero", false, Scale_MIX_test);
    GetKIDPDFs(MC_dirname_UUBAR, UUBAR_KID_m, "Continuum", false, Scale_UUBAR_test);
    GetKIDPDFs(MC_dirname_DDBAR, DDBAR_KID_m, "Continuum", false, Scale_DDBAR_test);
    GetKIDPDFs(MC_dirname_SSBAR, SSBAR_KID_m, "Continuum", false, Scale_SSBAR_test);
    GetKIDPDFs(MC_dirname_CHARM, CHARM_KID_m, "Continuum", false, Scale_CHARM_test);

    // get Kff uncertainty pdfs
    TH1D* Signal_Kff_temp = new TH1D("Signal_Kff_temp", "Signal_Kff_temp", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kff_array[7] = { Signal_Kff_temp, Signal_Kff1_p, Signal_Kff1_m, Signal_Kff2_p, Signal_Kff2_m, Signal_Kff3_p, Signal_Kff3_m };
    double Correction_factor_BR_Kplus[7] = { 0.0 };
    double Correction_factor_BR_Kzero[7] = { 0.0 };

    GetKffPDFs(MC_dirname_Knunu, Signal_Kff_array, Correction_factor_BR_Kplus, "Bplus", 1, Scale_Kplus_test);
    for (int i = 0; i < 7; i++) GetNominalPDFs(MC_dirname_Kstarnunu, Signal_Kff_array[i], "Bplus", Scale_Kplusstar_test);
    for (int i = 0; i < 7; i++) {
        const double non_res_Xsu_correction_factor = (BR_Xsu_nonresonant_nunubar + BR_Kplus_nunubar - BR_Kplus_nunubar * Correction_factor_BR_Kplus[i]) / BR_Xsu_nonresonant_nunubar;
        GetNominalPDFs(MC_dirname_Xsununu, Signal_Kff_array[i], "Bplus", Scale_Xsu_nonresonant_test * non_res_Xsu_correction_factor);
    }
    GetKffPDFs(MC_dirname_K0nunu, Signal_Kff_array, Correction_factor_BR_Kzero, "Bzero", 0, Scale_K0_test);
    for (int i = 0; i < 7; i++) GetNominalPDFs(MC_dirname_K0starnunu, Signal_Kff_array[i], "Bzero", Scale_K0star_test);
    for (int i = 0; i < 7; i++) {
        const double non_res_Xsd_correction_factor = (BR_Xsd_nonresonant_nunubar + BR_K0_nunubar - BR_K0_nunubar * Correction_factor_BR_Kzero[i]) / BR_Xsd_nonresonant_nunubar;
        GetNominalPDFs(MC_dirname_Xsdnunu, Signal_Kff_array[i], "Bzero", Scale_Xsd_nonresonant_test * non_res_Xsd_correction_factor);
    }

    // get Kstar ff uncertainty pdfs
    TH1D* Signal_Kstarff_temp = new TH1D("Signal_Kstarff_temp", "Signal_Kstarff_temp", RarityBins, 0.7, 1.0);
    TH1D* Signal_Kstarff_array[19] = { Signal_Kstarff_temp,
        Signal_Kstarff1_p, Signal_Kstarff1_m,
        Signal_Kstarff2_p, Signal_Kstarff2_m,
        Signal_Kstarff3_p, Signal_Kstarff3_m,
        Signal_Kstarff4_p, Signal_Kstarff4_m,
        Signal_Kstarff5_p, Signal_Kstarff5_m,
        Signal_Kstarff6_p, Signal_Kstarff6_m,
        Signal_Kstarff7_p, Signal_Kstarff7_m,
        Signal_Kstarff8_p, Signal_Kstarff8_m,
        Signal_Kstarff9_p, Signal_Kstarff9_m
    };
    double Correction_factor_BR_Kstarplus[19] = { 0.0 };
    double Correction_factor_BR_Kstarzero[19] = { 0.0 };

    for (int i = 0; i < 19; i++) GetNominalPDFs(MC_dirname_Knunu, Signal_Kstarff_array[i], "Bplus", Scale_Kplus_test);
    GetKstarffPDFs(MC_dirname_Kstarnunu, Signal_Kstarff_array, Correction_factor_BR_Kstarplus, "Bplus", 1, Scale_Kplusstar_test);
    for (int i = 0; i < 19; i++) {
        const double non_res_Xsu_correction_factor = (BR_Xsu_nonresonant_nunubar + BR_Kplusstar_nunubar - BR_Kplusstar_nunubar * Correction_factor_BR_Kstarplus[i]) / BR_Xsu_nonresonant_nunubar;
        GetNominalPDFs(MC_dirname_Xsununu, Signal_Kstarff_array[i], "Bplus", Scale_Xsu_nonresonant_test * non_res_Xsu_correction_factor);
    }
    for (int i = 0; i < 19; i++) GetNominalPDFs(MC_dirname_K0nunu, Signal_Kstarff_array[i], "Bzero", Scale_K0_test);
    GetKstarffPDFs(MC_dirname_K0starnunu, Signal_Kstarff_array, Correction_factor_BR_Kstarzero, "Bzero", 0, Scale_K0star_test);
    for (int i = 0; i < 19; i++) {
        const double non_res_Xsd_correction_factor = (BR_Xsd_nonresonant_nunubar + BR_K0star_nunubar - BR_K0star_nunubar * Correction_factor_BR_Kstarzero[i]) / BR_Xsd_nonresonant_nunubar;
        GetNominalPDFs(MC_dirname_Xsdnunu, Signal_Kstarff_array[i], "Bzero", Scale_Xsd_nonresonant_test * non_res_Xsd_correction_factor);
    }

    // get fraction uncertainty pdfs
    const double total_Xs_Num_plus = Lumi_validation_MC * (2.0 * N_BpBp_1invab * (BR_Xs_nunubar + Sigma_BR_Xs_nunubar) + 2.0 * N_B0B0_1invab * (BR_Xs_nunubar + Sigma_BR_Xs_nunubar));
    const double total_Xs_Num_minus = Lumi_validation_MC * (2.0 * N_BpBp_1invab * (BR_Xs_nunubar - Sigma_BR_Xs_nunubar) + 2.0 * N_B0B0_1invab * (BR_Xs_nunubar - Sigma_BR_Xs_nunubar));

    GetNominalPDFs(MC_dirname_Knunu, Signal_Kfrac_p, "Bplus", Scale_Kplus_test * ((BR_Kplus_nunubar + Sigma_BR_Kplus_nunubar) / BR_Kplus_nunubar));
    GetNominalPDFs(MC_dirname_Kstarnunu, Signal_Kfrac_p, "Bplus", Scale_Kplusstar_test);
    GetNominalPDFs(MC_dirname_Xsununu, Signal_Kfrac_p, "Bplus", Scale_Xsu_nonresonant_test * ((BR_Xsu_nonresonant_nunubar - Sigma_BR_Kplus_nunubar) / BR_Xsu_nonresonant_nunubar));
    GetNominalPDFs(MC_dirname_K0nunu, Signal_Kfrac_p, "Bzero", Scale_K0_test * ((BR_K0_nunubar + Sigma_BR_K0_nunubar) / BR_K0_nunubar));
    GetNominalPDFs(MC_dirname_K0starnunu, Signal_Kfrac_p, "Bzero", Scale_K0star_test);
    GetNominalPDFs(MC_dirname_Xsdnunu, Signal_Kfrac_p, "Bzero", Scale_Xsd_nonresonant_test * ((BR_Xsd_nonresonant_nunubar - Sigma_BR_K0_nunubar) / BR_Xsd_nonresonant_nunubar));

    GetNominalPDFs(MC_dirname_Knunu, Signal_Kfrac_m, "Bplus", Scale_Kplus_test * ((BR_Kplus_nunubar - Sigma_BR_Kplus_nunubar) / BR_Kplus_nunubar));
    GetNominalPDFs(MC_dirname_Kstarnunu, Signal_Kfrac_m, "Bplus", Scale_Kplusstar_test);
    GetNominalPDFs(MC_dirname_Xsununu, Signal_Kfrac_m, "Bplus", Scale_Xsu_nonresonant_test * ((BR_Xsu_nonresonant_nunubar + Sigma_BR_Kplus_nunubar) / BR_Xsu_nonresonant_nunubar));
    GetNominalPDFs(MC_dirname_K0nunu, Signal_Kfrac_m, "Bzero", Scale_K0_test * ((BR_K0_nunubar - Sigma_BR_K0_nunubar) / BR_K0_nunubar));
    GetNominalPDFs(MC_dirname_K0starnunu, Signal_Kfrac_m, "Bzero", Scale_K0star_test);
    GetNominalPDFs(MC_dirname_Xsdnunu, Signal_Kfrac_m, "Bzero", Scale_Xsd_nonresonant_test * ((BR_Xsd_nonresonant_nunubar + Sigma_BR_K0_nunubar) / BR_Xsd_nonresonant_nunubar));

    GetNominalPDFs(MC_dirname_Knunu, Signal_Kstarfrac_p, "Bplus", Scale_Kplus_test);
    GetNominalPDFs(MC_dirname_Kstarnunu, Signal_Kstarfrac_p, "Bplus", Scale_Kplusstar_test * ((BR_Kplusstar_nunubar + Sigma_BR_Kplusstar_nunubar) / BR_Kplusstar_nunubar));
    GetNominalPDFs(MC_dirname_Xsununu, Signal_Kstarfrac_p, "Bplus", Scale_Xsu_nonresonant_test * ((BR_Xsu_nonresonant_nunubar - Sigma_BR_Kplusstar_nunubar) / BR_Xsu_nonresonant_nunubar));
    GetNominalPDFs(MC_dirname_K0nunu, Signal_Kstarfrac_p, "Bzero", Scale_K0_test);
    GetNominalPDFs(MC_dirname_K0starnunu, Signal_Kstarfrac_p, "Bzero", Scale_K0star_test * ((BR_K0star_nunubar + Sigma_BR_K0star_nunubar) / BR_K0star_nunubar));
    GetNominalPDFs(MC_dirname_Xsdnunu, Signal_Kstarfrac_p, "Bzero", Scale_Xsd_nonresonant_test * ((BR_Xsd_nonresonant_nunubar - Sigma_BR_K0star_nunubar) / BR_Xsd_nonresonant_nunubar));

    GetNominalPDFs(MC_dirname_Knunu, Signal_Kstarfrac_m, "Bplus", Scale_Kplus_test);
    GetNominalPDFs(MC_dirname_Kstarnunu, Signal_Kstarfrac_m, "Bplus", Scale_Kplusstar_test * ((BR_Kplusstar_nunubar - Sigma_BR_Kplusstar_nunubar) / BR_Kplusstar_nunubar));
    GetNominalPDFs(MC_dirname_Xsununu, Signal_Kstarfrac_m, "Bplus", Scale_Xsu_nonresonant_test * ((BR_Xsu_nonresonant_nunubar + Sigma_BR_Kplusstar_nunubar) / BR_Xsu_nonresonant_nunubar));
    GetNominalPDFs(MC_dirname_K0nunu, Signal_Kstarfrac_m, "Bzero", Scale_K0_test);
    GetNominalPDFs(MC_dirname_K0starnunu, Signal_Kstarfrac_m, "Bzero", Scale_K0star_test * ((BR_K0star_nunubar - Sigma_BR_K0star_nunubar) / BR_K0star_nunubar));
    GetNominalPDFs(MC_dirname_Xsdnunu, Signal_Kstarfrac_m, "Bzero", Scale_Xsd_nonresonant_test * ((BR_Xsd_nonresonant_nunubar + Sigma_BR_K0star_nunubar) / BR_Xsd_nonresonant_nunubar));

    // get fragmentation uncertainty pdfs
    DecayModeMC Fragmentation_types[2 * (MAX_NUM_DECAYMODE_MC - 6 + 1)] = {
        Xsu2KcPi0_MC, Xsu2K0Pic_MC, Xsu2KcPicPic_MC, Xsu2K0PicPi0_MC,
        Xsu2KcPicPicPi0_MC, Xsu2K0PicPicPic_MC, Xsu2KcPicPicPicPic_MC, Xsu2K0PicPicPicPi0_MC,
        Xsu2KcPi0Pi0_MC, Xsu2K0PicPi0Pi0_MC, Xsu2KcPicPicPi0Pi0_MC, Xsu2KcKcKc_MC,
        Xsu2KcKcK0Pic_MC, Xsu2KcKcKcPi0_MC, other,
        Xsd2KcPic_MC, Xsd2K0Pi0_MC, Xsd2KcPicPi0_MC, Xsd2K0PicPic_MC,
        Xsd2KcPicPicPic_MC, Xsd2K0PicPicPi0_MC, Xsd2KcPicPicPicPi0_MC, Xsd2K0PicPicPicPic_MC,
        Xsd2K0Pi0Pi0_MC, Xsd2KcPicPi0Pi0_MC, Xsd2K0PicPicPi0Pi0_MC, Xsd2KcKcK0_MC,
        Xsd2KcKcKcPic_MC, Xsd2KcKcK0Pi0_MC, other,
        Xsu2KcPi0_MC, Xsu2K0Pic_MC, Xsu2KcPicPic_MC, Xsu2K0PicPi0_MC,
        Xsu2KcPicPicPi0_MC, Xsu2K0PicPicPic_MC, Xsu2KcPicPicPicPic_MC, Xsu2K0PicPicPicPi0_MC,
        Xsu2KcPi0Pi0_MC, Xsu2K0PicPi0Pi0_MC, Xsu2KcPicPicPi0Pi0_MC, Xsu2KcKcKc_MC,
        Xsu2KcKcK0Pic_MC, Xsu2KcKcKcPi0_MC, other,
        Xsd2KcPic_MC, Xsd2K0Pi0_MC, Xsd2KcPicPi0_MC, Xsd2K0PicPic_MC,
        Xsd2KcPicPicPic_MC, Xsd2K0PicPicPi0_MC, Xsd2KcPicPicPicPi0_MC, Xsd2K0PicPicPicPic_MC,
        Xsd2K0Pi0Pi0_MC, Xsd2KcPicPi0Pi0_MC, Xsd2K0PicPicPi0Pi0_MC, Xsd2KcKcK0_MC,
        Xsd2KcKcKcPic_MC, Xsd2KcKcK0Pi0_MC, other
    };
    for (unsigned int i = 0; i < name_Signal_Fragmentation.size(); i++) {
        bool IsItUp = true;
        if(name_Signal_Fragmentation.at(i).find("_p") != string::npos) IsItUp = true;
        else if(name_Signal_Fragmentation.at(i).find("_m") != string::npos) IsItUp = false;
        else {
            printf("[ERROR] unexpected error!\n");
            exit(1);
        }

        GetNominalPDFs(MC_dirname_Knunu, Signal_Fragmentaions.at(i), "Bplus", Scale_Kplus_test);
        GetNominalPDFs(MC_dirname_Kstarnunu, Signal_Fragmentaions.at(i), "Bplus", Scale_Kplusstar_test);
        if(name_Signal_Fragmentation.at(i).find("Xsu") != string::npos) {
            GetFragmentationPDFs(MC_dirname_Xsununu, Signal_Fragmentaions.at(i), "Bplus", Fragmentation_types[i], IsItUp, Scale_Xsu_nonresonant_test);
        }
        else {
            GetNominalPDFs(MC_dirname_Xsununu, Signal_Fragmentaions.at(i), "Bplus", Scale_Xsu_nonresonant_test);
        }
        GetNominalPDFs(MC_dirname_K0nunu, Signal_Fragmentaions.at(i), "Bzero", Scale_K0_test);
        GetNominalPDFs(MC_dirname_K0starnunu, Signal_Fragmentaions.at(i), "Bzero", Scale_K0star_test);
        if (name_Signal_Fragmentation.at(i).find("Xsd") != string::npos) {
            GetFragmentationPDFs(MC_dirname_Xsdnunu, Signal_Fragmentaions.at(i), "Bzero", Fragmentation_types[i], IsItUp, Scale_Xsd_nonresonant_test);
        }
        else {
            GetNominalPDFs(MC_dirname_Xsdnunu, Signal_Fragmentaions.at(i), "Bzero", Scale_Xsd_nonresonant_test);
        }

    }

    // get pf uncertainty pdfs
    GetNominalPDFs(MC_dirname_Knunu, Signal_pf_p, "Bplus", Scale_Kplus_test);
    GetNominalPDFs(MC_dirname_Kstarnunu, Signal_pf_p, "Bplus", Scale_Kplusstar_test);
    GetNominalPDFs(MC_dirname_Xsununu_Hpf, Signal_pf_p, "Bplus", Scale_Xsu_nonresonant_syst);
    GetNominalPDFs(MC_dirname_K0nunu, Signal_pf_p, "Bzero", Scale_K0_test);
    GetNominalPDFs(MC_dirname_K0starnunu, Signal_pf_p, "Bzero", Scale_K0star_test);
    GetNominalPDFs(MC_dirname_Xsdnunu_Hpf, Signal_pf_p, "Bzero", Scale_Xsd_nonresonant_syst);

    GetNominalPDFs(MC_dirname_Knunu, Signal_pf_m, "Bplus", Scale_Kplus_test);
    GetNominalPDFs(MC_dirname_Kstarnunu, Signal_pf_m, "Bplus", Scale_Kplusstar_test);
    GetNominalPDFs(MC_dirname_Xsununu_Lpf, Signal_pf_m, "Bplus", Scale_Xsu_nonresonant_syst);
    GetNominalPDFs(MC_dirname_K0nunu, Signal_pf_m, "Bzero", Scale_K0_test);
    GetNominalPDFs(MC_dirname_K0starnunu, Signal_pf_m, "Bzero", Scale_K0star_test);
    GetNominalPDFs(MC_dirname_Xsdnunu_Lpf, Signal_pf_m, "Bzero", Scale_Xsd_nonresonant_syst);

    // mb
    GetNominalPDFs(MC_dirname_Knunu, Signal_mb_p, "Bplus", Scale_Kplus_test);
    GetNominalPDFs(MC_dirname_Kstarnunu, Signal_mb_p, "Bplus", Scale_Kplusstar_test);
    GetNominalPDFs(MC_dirname_Xsununu_Hmb, Signal_mb_p, "Bplus", Scale_Xsu_nonresonant_syst);
    GetNominalPDFs(MC_dirname_K0nunu, Signal_mb_p, "Bzero", Scale_K0_test);
    GetNominalPDFs(MC_dirname_K0starnunu, Signal_mb_p, "Bzero", Scale_K0star_test);
    GetNominalPDFs(MC_dirname_Xsdnunu_Hmb, Signal_mb_p, "Bzero", Scale_Xsd_nonresonant_syst);

    GetNominalPDFs(MC_dirname_Knunu, Signal_mb_m, "Bplus", Scale_Kplus_test);
    GetNominalPDFs(MC_dirname_Kstarnunu, Signal_mb_m, "Bplus", Scale_Kplusstar_test);
    GetNominalPDFs(MC_dirname_Xsununu_Lmb, Signal_mb_m, "Bplus", Scale_Xsu_nonresonant_syst);
    GetNominalPDFs(MC_dirname_K0nunu, Signal_mb_m, "Bzero", Scale_K0_test);
    GetNominalPDFs(MC_dirname_K0starnunu, Signal_mb_m, "Bzero", Scale_K0star_test);
    GetNominalPDFs(MC_dirname_Xsdnunu_Lmb, Signal_mb_m, "Bzero", Scale_Xsd_nonresonant_syst);

    // transition
    GetNominalPDFs(MC_dirname_Knunu, Signal_transition_p, "Bplus", Scale_Kplus_test);
    GetNominalPDFs(MC_dirname_Kstarnunu, Signal_transition_p, "Bplus", Scale_Kplusstar_test);
    GetNominalPDFs(MC_dirname_Xsununu_Htransition, Signal_transition_p, "Bplus", Scale_Xsu_nonresonant_syst);
    GetNominalPDFs(MC_dirname_K0nunu, Signal_transition_p, "Bzero", Scale_K0_test);
    GetNominalPDFs(MC_dirname_K0starnunu, Signal_transition_p, "Bzero", Scale_K0star_test);
    GetNominalPDFs(MC_dirname_Xsdnunu_Htransition, Signal_transition_p, "Bzero", Scale_Xsd_nonresonant_syst);

    GetNominalPDFs(MC_dirname_Knunu, Signal_transition_m, "Bplus", Scale_Kplus_test);
    GetNominalPDFs(MC_dirname_Kstarnunu, Signal_transition_m, "Bplus", Scale_Kplusstar_test);
    GetNominalPDFs(MC_dirname_Xsununu_Ltransition, Signal_transition_m, "Bplus", Scale_Xsu_nonresonant_syst);
    GetNominalPDFs(MC_dirname_K0nunu, Signal_transition_m, "Bzero", Scale_K0_test);
    GetNominalPDFs(MC_dirname_K0starnunu, Signal_transition_m, "Bzero", Scale_K0star_test);
    GetNominalPDFs(MC_dirname_Xsdnunu_Ltransition, Signal_transition_m, "Bzero", Scale_Xsd_nonresonant_syst);

    // mKstar fixed
    GetNominalPDFs(MC_dirname_Knunu, Signal_mKstar_p, "Bplus", Scale_Kplus_test);
    GetNominalPDFs(MC_dirname_Kstarnunu_mKstarfixed, Signal_mKstar_p, "Bplus", Scale_Kplusstar_syst);
    GetNominalPDFs(MC_dirname_Xsununu, Signal_mKstar_p, "Bplus", Scale_Xsu_nonresonant_test);
    GetNominalPDFs(MC_dirname_K0nunu, Signal_mKstar_p, "Bzero", Scale_K0_test);
    GetNominalPDFs(MC_dirname_K0starnunu_mKstarfixed, Signal_mKstar_p, "Bzero", Scale_K0star_syst);
    GetNominalPDFs(MC_dirname_Xsdnunu, Signal_mKstar_p, "Bzero", Scale_Xsd_nonresonant_test);

    GetNominalPDFs(MC_dirname_Knunu, Signal_mKstar_m, "Bplus", Scale_Kplus_test);
    GetNominalPDFs(MC_dirname_Kstarnunu_mKstarfixed, Signal_mKstar_m, "Bplus", Scale_Kplusstar_syst);
    GetNominalPDFs(MC_dirname_Xsununu, Signal_mKstar_m, "Bplus", Scale_Xsu_nonresonant_test);
    GetNominalPDFs(MC_dirname_K0nunu, Signal_mKstar_m, "Bzero", Scale_K0_test);
    GetNominalPDFs(MC_dirname_K0starnunu_mKstarfixed, Signal_mKstar_m, "Bzero", Scale_K0star_syst);
    GetNominalPDFs(MC_dirname_Xsdnunu, Signal_mKstar_m, "Bzero", Scale_Xsd_nonresonant_test);
    /* ====================================== */



    /* ====================================== */
    // define path for Data
    const char* DATA_dirname_Knunu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SIGNAL_analysis/validation_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Knunu";
    const char* DATA_dirname_Kstarnunu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SIGNAL_analysis/validation_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Kstarnunu";
    const char* DATA_dirname_Xsununu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SIGNAL_analysis/validation_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Xsnunu";
    const char* DATA_dirname_K0nunu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SIGNAL_analysis/validation_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0nunu";
    const char* DATA_dirname_K0starnunu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SIGNAL_analysis/validation_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0starnunu";
    const char* DATA_dirname_Xsdnunu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SIGNAL_analysis/validation_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02Xsnunu";

    const char* DATA_dirname_CHG = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/CHG_analysis/validation_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_MIX = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/MIX_analysis/validation_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_UUBAR = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/UUBAR_analysis/validation_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_DDBAR = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/DDBAR_analysis/validation_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_SSBAR = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SSBAR_analysis/validation_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_CHARM = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/CHARM_analysis/validation_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    /* ====================================== */



    /* ====================================== */
    // Define histograms for Data
    TH1D* total_DATA = new TH1D("total_DATA", "total_DATA", RarityBins, 0.7, 1.0);
    /* ====================================== */



    /* ====================================== */
    // Save data
    GetNominalPDFs(DATA_dirname_Knunu, total_DATA, "Bplus", Scale_Kplus_test);
    GetNominalPDFs(DATA_dirname_Kstarnunu, total_DATA, "Bplus", Scale_Kplusstar_test);
    GetNominalPDFs(DATA_dirname_Xsununu, total_DATA, "Bplus", Scale_Xsu_nonresonant_test);
    GetNominalPDFs(DATA_dirname_K0nunu, total_DATA, "Bzero", Scale_K0_test);
    GetNominalPDFs(DATA_dirname_K0starnunu, total_DATA, "Bzero", Scale_K0star_test);
    GetNominalPDFs(DATA_dirname_Xsdnunu, total_DATA, "Bzero", Scale_Xsd_nonresonant_test);

    GetNominalPDFs(DATA_dirname_CHG, total_DATA, "Bplus", 1.0);
    GetNominalPDFs(DATA_dirname_MIX, total_DATA, "Bzero", 1.0);
    GetNominalPDFs(DATA_dirname_UUBAR, total_DATA, "Continuum", 1.0);
    GetNominalPDFs(DATA_dirname_DDBAR, total_DATA, "Continuum", 1.0);
    GetNominalPDFs(DATA_dirname_SSBAR, total_DATA, "Continuum", 1.0);
    GetNominalPDFs(DATA_dirname_CHARM, total_DATA, "Continuum", 1.0);
    /* ====================================== */



    /* ====================================== */
    // Save histograms
    TFile* file = new TFile("PDFandDATA.root", "RECREATE");

    Signal_nominal->Write();
    CHG_nominal->Write();
    MIX_nominal->Write();
    UUBAR_nominal->Write();
    DDBAR_nominal->Write();
    SSBAR_nominal->Write();
    CHARM_nominal->Write();

    // FEI uncertainty, BKGs will be cared by overall syst
    Signal_FEI_charged_p->Write();
    Signal_FEI_charged_m->Write();
    Signal_FEI_neutral_p->Write();
    Signal_FEI_neutral_m->Write();

    // pi0 uncertainty
    Signal_pi0_p->Write();
    CHG_pi0_p->Write();
    MIX_pi0_p->Write();
    UUBAR_pi0_p->Write();
    DDBAR_pi0_p->Write();
    SSBAR_pi0_p->Write();
    CHARM_pi0_p->Write();
    Signal_pi0_m->Write();
    CHG_pi0_m->Write();
    MIX_pi0_m->Write();
    UUBAR_pi0_m->Write();
    DDBAR_pi0_m->Write();
    SSBAR_pi0_m->Write();
    CHARM_pi0_m->Write();

    // track uncertainty
    Signal_track_p->Write();
    CHG_track_p->Write();
    MIX_track_p->Write();
    UUBAR_track_p->Write();
    DDBAR_track_p->Write();
    SSBAR_track_p->Write();
    CHARM_track_p->Write();
    Signal_track_m->Write();
    CHG_track_m->Write();
    MIX_track_m->Write();
    UUBAR_track_m->Write();
    DDBAR_track_m->Write();
    SSBAR_track_m->Write();
    CHARM_track_m->Write();

    // KS0 uncertainty
    Signal_KS0_p->Write();
    CHG_KS0_p->Write();
    MIX_KS0_p->Write();
    UUBAR_KS0_p->Write();
    DDBAR_KS0_p->Write();
    SSBAR_KS0_p->Write();
    CHARM_KS0_p->Write();
    Signal_KS0_m->Write();
    CHG_KS0_m->Write();
    MIX_KS0_m->Write();
    UUBAR_KS0_m->Write();
    DDBAR_KS0_m->Write();
    SSBAR_KS0_m->Write();
    CHARM_KS0_m->Write();

    // Kaon PID uncertainty
    Signal_KID_p->Write();
    CHG_KID_p->Write();
    MIX_KID_p->Write();
    UUBAR_KID_p->Write();
    DDBAR_KID_p->Write();
    SSBAR_KID_p->Write();
    CHARM_KID_p->Write();
    Signal_KID_m->Write();
    CHG_KID_m->Write();
    MIX_KID_m->Write();
    UUBAR_KID_m->Write();
    DDBAR_KID_m->Write();
    SSBAR_KID_m->Write();
    CHARM_KID_m->Write();

    // K nu nubar form factor
    Signal_Kff1_p->Write();
    Signal_Kff1_m->Write();
    Signal_Kff2_p->Write();
    Signal_Kff2_m->Write();
    Signal_Kff3_p->Write();
    Signal_Kff3_m->Write();

    // Kstar nu nubar form factor
    Signal_Kstarff1_p->Write();
    Signal_Kstarff1_m->Write();
    Signal_Kstarff2_p->Write();
    Signal_Kstarff2_m->Write();
    Signal_Kstarff3_p->Write();
    Signal_Kstarff3_m->Write();
    Signal_Kstarff4_p->Write();
    Signal_Kstarff4_m->Write();
    Signal_Kstarff5_p->Write();
    Signal_Kstarff5_m->Write();
    Signal_Kstarff6_p->Write();
    Signal_Kstarff6_m->Write();
    Signal_Kstarff7_p->Write();
    Signal_Kstarff7_m->Write();
    Signal_Kstarff8_p->Write();
    Signal_Kstarff8_m->Write();
    Signal_Kstarff9_p->Write();
    Signal_Kstarff9_m->Write();

    // faction
    Signal_Kfrac_p->Write();
    Signal_Kfrac_m->Write();
    Signal_Kstarfrac_p->Write();
    Signal_Kstarfrac_m->Write();

    // fragmentation uncertainty
    for (unsigned int i = 0; i < Signal_Fragmentaions.size(); i++) Signal_Fragmentaions.at(i)->Write();

    // pf uncertainty
    Signal_pf_p->Write();
    Signal_pf_m->Write();

    // mb uncertainty
    Signal_mb_p->Write();
    Signal_mb_m->Write();

    // transition uncertainty
    Signal_transition_p->Write();
    Signal_transition_m->Write();

    // fixed mKstar uncertainty
    Signal_mKstar_p->Write();
    Signal_mKstar_m->Write();

    total_DATA->Write();

    file->Close();
    /* ====================================== */
}
