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
#include <numeric>
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
#include "TH2.h"
#include "TF1.h"
#include "RooDataHist.h"
#include "template.h"
#include "constants.h"
#include "ObtainWeight.h"
#include "correctors.h"
#include "base.h"
using namespace RooFit;
using std::string;
using std::to_string;

# define Knunu_only false

Corrector corrector;
Corrector_FEI corrector_FEI;
Corrector_PID corrector_PID;
Corrector_pi0 corrector_pi0;
Corrector_FakePID corrector_FakePID;
Corrector_Knn corrector_Knn;
Corrector_Xsnn corrector_Xsnn;
Corrector_Multiplicity corrector_Multiplicity;
Corrector_KpKLKL corrector_KpKLKL;
Corrector_KSKLKL corrector_KSKLKL;
Corrector_KstarKLKL corrector_KstarKLKL;
Corrector_XsKLKL corrector_XsKLKL;
Corrector_BtoDtoXKL corrector_BtoDtoXKL;
Corrector_Fragmentation corrector_Fragmentation;

# define MCTYPE "MC15rd"

TH1D* Xsu_Hmb_weight;
TH1D* Xsu_Lmb_weight;
TH1D* Xsu_Hpf_weight;
TH1D* Xsu_Lpf_weight;
TH1D* Xsu_Htransition_weight;
TH1D* Xsu_Ltransition_weight;
TH1D* Xsd_Hmb_weight;
TH1D* Xsd_Lmb_weight;
TH1D* Xsd_Hpf_weight;
TH1D* Xsd_Lpf_weight;
TH1D* Xsd_Htransition_weight;
TH1D* Xsd_Ltransition_weight;
TH1D* Kstar_delta_weight;
TH1D* K0star_delta_weight;
double ReadWeightHist(TH1D* hist, double value);

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

/* ====================================== */

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

double BDTcToWeight(double BDTc) {

    if (BDTc > (5.0 / 6.0)) return std::sqrt(5.0);
    else return std::sqrt(BDTc / (1.0 - BDTc));

}

double GetNominalPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get nominal PDF with appropriate correction
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
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
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double invM_Xnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;
    double N_Xplusnn = 0;
    double N_Xzeronn = 0;

    double Ngamma_v200 = -1;

    double s13_KpKLKL = -1;
    double s23_KpKLKL = -1;
    double nB2KpKLKL_all_KpKLKL = -1;
    double nB2KpKLKL_NR_KpKLKL = -1;

    double s13_KSKLKL = -1;
    double s23_KSKLKL = -1;
    double s12_KSKLKL = -1;
    double nB2KSKLKL_all_KSKLKL = -1;
    double nB2KSKLKL_NR_KSKLKL = -1;

    double nKL_XKLKL = -1;
    double XKLKL_E_1st = -1;
    double XKLKL_px_1st = -1;
    double XKLKL_py_1st = -1;
    double XKLKL_pz_1st = -1;
    double XKLKL_E_2nd = -1;
    double XKLKL_px_2nd = -1;
    double XKLKL_py_2nd = -1;
    double XKLKL_pz_2nd = -1;
    double nB2KstarKLKL = -1;
    double nB02KstarKLKL = -1;

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    double Bsig_M = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clXnn__bc", &N_Xplusnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clXnn__bc", &N_Xzeronn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clXnn__bc", &invM_Xnn);

        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s13_KpKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s23_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_all__bc", &nB2KpKLKL_all_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_NR__bc", &nB2KpKLKL_NR_KpKLKL);

        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s13_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &s23_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s12_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_all__bc", &nB2KSKLKL_all_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_NR__bc", &nB2KSKLKL_NR_KSKLKL);

        tree_upsilon->SetBranchAddress("nParticlesInList__boK_L0__clXKLKL__bc", &nKL_XKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__spE__bc", &XKLKL_E_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__sppx__bc", &XKLKL_px_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__sppy__bc", &XKLKL_py_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__sppz__bc", &XKLKL_pz_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__spE__bc", &XKLKL_E_2nd);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppx__bc", &XKLKL_px_2nd);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppy__bc", &XKLKL_py_2nd);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppz__bc", &XKLKL_pz_2nd);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarKLKL__bc", &nB2KstarKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstarKLKL__bc", &nB02KstarKLKL);

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            // select B+ --> K+ nu nubar reconstruction only
            if ((std::abs(Upsilon_ID) < MyEPSILON) && (std::abs(Bsig_ID) < MyEPSILON)) {}
            else continue;

            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE), temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE), temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(0, i_fake, MCTYPE), temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(1, i_fake, MCTYPE), temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(2, i_fake, MCTYPE), temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(3, i_fake, MCTYPE), temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(0, i_fake, MCTYPE), temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(1, i_fake, MCTYPE), temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(2, i_fake, MCTYPE), temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(3, i_fake, MCTYPE), temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B --> K* KL KL correction factor
            double Correction_KstarKLKL = corrector_KstarKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KstarKLKL + nB02KstarKLKL);

            // B --> X KL KL correction factor
            double Correction_XKLKL = corrector_XsKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KpKLKL_all_KpKLKL, nB2KSKLKL_all_KSKLKL, nB2KstarKLKL + nB02KstarKLKL, nKL_XKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if(CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            Nevt = Nevt + FillTemplate(hist, MVA_var, total_weight, Bsig_M);

        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double GetTrackPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, bool IsItUp, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get track efficiency uncertainty PDF with appropriate correction
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
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
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double invM_Xnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;
    double N_Xplusnn = 0;
    double N_Xzeronn = 0;

    double Ngamma_v200 = -1;

    double s13_KpKLKL = -1;
    double s23_KpKLKL = -1;
    double nB2KpKLKL_all_KpKLKL = -1;
    double nB2KpKLKL_NR_KpKLKL = -1;

    double s13_KSKLKL = -1;
    double s23_KSKLKL = -1;
    double s12_KSKLKL = -1;
    double nB2KSKLKL_all_KSKLKL = -1;
    double nB2KSKLKL_NR_KSKLKL = -1;

    double nKL_XKLKL = -1;
    double XKLKL_E_1st = -1;
    double XKLKL_px_1st = -1;
    double XKLKL_py_1st = -1;
    double XKLKL_pz_1st = -1;
    double XKLKL_E_2nd = -1;
    double XKLKL_px_2nd = -1;
    double XKLKL_py_2nd = -1;
    double XKLKL_pz_2nd = -1;
    double nB2KstarKLKL = -1;
    double nB02KstarKLKL = -1;

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    double Bsig_M = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clXnn__bc", &N_Xplusnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clXnn__bc", &N_Xzeronn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clXnn__bc", &invM_Xnn);

        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s13_KpKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s23_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_all__bc", &nB2KpKLKL_all_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_NR__bc", &nB2KpKLKL_NR_KpKLKL);

        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s13_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &s23_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s12_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_all__bc", &nB2KSKLKL_all_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_NR__bc", &nB2KSKLKL_NR_KSKLKL);

        tree_upsilon->SetBranchAddress("nParticlesInList__boK_L0__clXKLKL__bc", &nKL_XKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__spE__bc", &XKLKL_E_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__sppx__bc", &XKLKL_px_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__sppy__bc", &XKLKL_py_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__sppz__bc", &XKLKL_pz_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__spE__bc", &XKLKL_E_2nd);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppx__bc", &XKLKL_px_2nd);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppy__bc", &XKLKL_py_2nd);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppz__bc", &XKLKL_pz_2nd);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarKLKL__bc", &nB2KstarKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstarKLKL__bc", &nB02KstarKLKL);

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            // select B+ --> K+ nu nubar reconstruction only
            if ((std::abs(Upsilon_ID) < MyEPSILON) && (std::abs(Bsig_ID) < MyEPSILON)) {}
            else continue;

            double Ntrack = GetNtracks(Upsilon_ID, Bsig_ID);

            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE), temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE), temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(0, i_fake, MCTYPE), temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(1, i_fake, MCTYPE), temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(2, i_fake, MCTYPE), temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(3, i_fake, MCTYPE), temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(0, i_fake, MCTYPE), temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(1, i_fake, MCTYPE), temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(2, i_fake, MCTYPE), temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(3, i_fake, MCTYPE), temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B --> K* KL KL correction factor
            double Correction_KstarKLKL = corrector_KstarKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KstarKLKL + nB02KstarKLKL);

            // B --> X KL KL correction factor
            double Correction_XKLKL = corrector_XsKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KpKLKL_all_KpKLKL, nB2KSKLKL_all_KSKLKL, nB2KstarKLKL + nB02KstarKLKL, nKL_XKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double Correction_track = 1.0;
            if (IsItUp == true) Correction_track = (1 + (Ntrack * track_rel_uncertainty / 100.0));
            else Correction_track = (1 - (Ntrack * track_rel_uncertainty / 100.0));

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * Correction_track;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            Nevt = Nevt + FillTemplate(hist, MVA_var, total_weight, Bsig_M);

        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double GetKS0PDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, bool IsItUp, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get KS0 uncertainty PDF with appropriate correction
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
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
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_KS0_3D_distance = -1;

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double invM_Xnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;
    double N_Xplusnn = 0;
    double N_Xzeronn = 0;

    double Ngamma_v200 = -1;

    double s13_KpKLKL = -1;
    double s23_KpKLKL = -1;
    double nB2KpKLKL_all_KpKLKL = -1;
    double nB2KpKLKL_NR_KpKLKL = -1;

    double s13_KSKLKL = -1;
    double s23_KSKLKL = -1;
    double s12_KSKLKL = -1;
    double nB2KSKLKL_all_KSKLKL = -1;
    double nB2KSKLKL_NR_KSKLKL = -1;

    double nKL_XKLKL = -1;
    double XKLKL_E_1st = -1;
    double XKLKL_px_1st = -1;
    double XKLKL_py_1st = -1;
    double XKLKL_pz_1st = -1;
    double XKLKL_E_2nd = -1;
    double XKLKL_px_2nd = -1;
    double XKLKL_py_2nd = -1;
    double XKLKL_pz_2nd = -1;
    double nB2KstarKLKL = -1;
    double nB02KstarKLKL = -1;

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    double Bsig_M = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);
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
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_3D_distance", &temp_KS0_3D_distance);
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clXnn__bc", &N_Xplusnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clXnn__bc", &N_Xzeronn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clXnn__bc", &invM_Xnn);

        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s13_KpKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s23_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_all__bc", &nB2KpKLKL_all_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_NR__bc", &nB2KpKLKL_NR_KpKLKL);

        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s13_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &s23_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s12_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_all__bc", &nB2KSKLKL_all_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_NR__bc", &nB2KSKLKL_NR_KSKLKL);

        tree_upsilon->SetBranchAddress("nParticlesInList__boK_L0__clXKLKL__bc", &nKL_XKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__spE__bc", &XKLKL_E_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__sppx__bc", &XKLKL_px_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__sppy__bc", &XKLKL_py_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__sppz__bc", &XKLKL_pz_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__spE__bc", &XKLKL_E_2nd);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppx__bc", &XKLKL_px_2nd);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppy__bc", &XKLKL_py_2nd);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppz__bc", &XKLKL_pz_2nd);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarKLKL__bc", &nB2KstarKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstarKLKL__bc", &nB02KstarKLKL);

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            // select B+ --> K+ nu nubar reconstruction only
            if ((std::abs(Upsilon_ID) < MyEPSILON) && (std::abs(Bsig_ID) < MyEPSILON)) {}
            else continue;

            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE), temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE), temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(0, i_fake, MCTYPE), temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(1, i_fake, MCTYPE), temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(2, i_fake, MCTYPE), temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(3, i_fake, MCTYPE), temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(0, i_fake, MCTYPE), temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(1, i_fake, MCTYPE), temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(2, i_fake, MCTYPE), temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(3, i_fake, MCTYPE), temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B --> K* KL KL correction factor
            double Correction_KstarKLKL = corrector_KstarKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KstarKLKL + nB02KstarKLKL);

            // B --> X KL KL correction factor
            double Correction_XKLKL = corrector_XsKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KpKLKL_all_KpKLKL, nB2KSKLKL_all_KSKLKL, nB2KstarKLKL + nB02KstarKLKL, nKL_XKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double KS0_correction = 1.0;
            if(IsItUp == true) KS0_correction = 1 + (KS0_rel_uncertainty * temp_KS0_3D_distance / 100.0);
            else KS0_correction = 1 - (KS0_rel_uncertainty * temp_KS0_3D_distance / 100.0);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * KS0_correction;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            Nevt = Nevt + FillTemplate(hist, MVA_var, total_weight, Bsig_M);

        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

int GetFEIcorrelatedPDFs(const char* dirname, TH1D* CHG_nominal_hist, TH1D* MIX_nominal_hist, TH1D* Signal_nominal_hist, TH1D*** CHG_hists, TH1D*** MIX_hists, TH1D*** Signal_hists) { // get shape sys histogram from txt file
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

    *CHG_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *MIX_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *Signal_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);

    for (int i = 0; i < Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_FEI_correlated" + std::to_string(i) + "_p").c_str(), ("CHG_FEI_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_FEI_correlated" + std::to_string(i) + "_p").c_str(), ("MIX_FEI_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists)[i] = new TH1D(("Signal_FEI_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_FEI_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
    }

    for (int i = Nentry; i < 2 * Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_FEI_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHG_FEI_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_FEI_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("MIX_FEI_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists)[i] = new TH1D(("Signal_FEI_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_FEI_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
    }

    fp = fopen(dirname, "r");
    for (int i = 0; i < Nentry; i++) {
        fscanf(fp, "%lf\n", &eigen_value);
        for (int j = 0; j < RarityBins * 3; j++) fscanf(fp, "%lf\n", &weight_sys[j]);

        for (int k = 0; k < RarityBins; k++) {
            (*CHG_hists)[i]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 1 * RarityBins]));
            (*Signal_hists)[i]->SetBinContent(k + 1, Signal_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 2 * RarityBins]));

            (*CHG_hists)[i + Nentry]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i + Nentry]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 1 * RarityBins]));
            (*Signal_hists)[i + Nentry]->SetBinContent(k + 1, Signal_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 2 * RarityBins]));
        }

    }

    fclose(fp);

    return Nentry;
}

void GetFEIUncorrelatedPDFs(const char* dirname, TH1D* CHG_hist, TH1D* MIX_hist, TH1D* Signal_hist) { // get shape sys histogram from txt file
    FILE* fp;
    fp = fopen(dirname, "r");

    double weight_sys[RarityBins * 3] = { 0.0 };
    for (int i = 0; i < RarityBins * 3; i++) fscanf(fp, "%lf\n", &weight_sys[i]);
    fclose(fp);

    for (int i = 0; i < RarityBins * 3; i++) weight_sys[i] = std::sqrt(weight_sys[i]);

    for (int i = 0; i < RarityBins; i++) CHG_hist->SetBinContent(i + 1, weight_sys[i]);
    for (int i = 0; i < RarityBins; i++) MIX_hist->SetBinContent(i + 1, weight_sys[RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) Signal_hist->SetBinContent(i + 1, weight_sys[2 * RarityBins + i]);
}

int GetPIDcorrelatedPDFs(const char* dirname, TH1D* CHG_nominal_hist, TH1D* MIX_nominal_hist, TH1D* UUBAR_nominal_hist, TH1D* DDBAR_nominal_hist, TH1D* SSBAR_nominal_hist, TH1D* CHARM_nominal_hist, TH1D* SIGNAL_nominal_hist, TH1D*** CHG_hists, TH1D*** MIX_hists, TH1D*** UUBAR_hists, TH1D*** DDBAR_hists, TH1D*** SSBAR_hists, TH1D*** CHARM_hists, TH1D*** SIGNAL_hists, bool IsItKID) { // get shape sys histogram from txt file
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

    *CHG_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *MIX_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *UUBAR_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *DDBAR_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *SSBAR_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *CHARM_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *SIGNAL_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);

    for (int i = 0; i < Nentry; i++) {
        if (IsItKID) {
            (*CHG_hists)[i] = new TH1D(("CHG_KID_correlated" + std::to_string(i) + "_p").c_str(), ("CHG_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*MIX_hists)[i] = new TH1D(("MIX_KID_correlated" + std::to_string(i) + "_p").c_str(), ("MIX_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*UUBAR_hists)[i] = new TH1D(("UUBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), ("UUBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*DDBAR_hists)[i] = new TH1D(("DDBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), ("DDBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*SSBAR_hists)[i] = new TH1D(("SSBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), ("SSBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*CHARM_hists)[i] = new TH1D(("CHARM_KID_correlated" + std::to_string(i) + "_p").c_str(), ("CHARM_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*SIGNAL_hists)[i] = new TH1D(("Signal_KID_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        }
        else {
            (*CHG_hists)[i] = new TH1D(("CHG_PID_correlated" + std::to_string(i) + "_p").c_str(), ("CHG_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*MIX_hists)[i] = new TH1D(("MIX_PID_correlated" + std::to_string(i) + "_p").c_str(), ("MIX_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*UUBAR_hists)[i] = new TH1D(("UUBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), ("UUBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*DDBAR_hists)[i] = new TH1D(("DDBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), ("DDBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*SSBAR_hists)[i] = new TH1D(("SSBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), ("SSBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*CHARM_hists)[i] = new TH1D(("CHARM_PID_correlated" + std::to_string(i) + "_p").c_str(), ("CHARM_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*SIGNAL_hists)[i] = new TH1D(("Signal_PID_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        }
    }

    for (int i = Nentry; i < 2* Nentry; i++) {
        if (IsItKID) {
            (*CHG_hists)[i] = new TH1D(("CHG_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHG_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*MIX_hists)[i] = new TH1D(("MIX_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("MIX_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*UUBAR_hists)[i] = new TH1D(("UUBAR_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("UUBAR_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*DDBAR_hists)[i] = new TH1D(("DDBAR_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("DDBAR_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*SSBAR_hists)[i] = new TH1D(("SSBAR_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("SSBAR_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*CHARM_hists)[i] = new TH1D(("CHARM_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHARM_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*SIGNAL_hists)[i] = new TH1D(("Signal_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        }
        else {
            (*CHG_hists)[i] = new TH1D(("CHG_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHG_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*MIX_hists)[i] = new TH1D(("MIX_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("MIX_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*UUBAR_hists)[i] = new TH1D(("UUBAR_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("UUBAR_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*DDBAR_hists)[i] = new TH1D(("DDBAR_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("DDBAR_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*SSBAR_hists)[i] = new TH1D(("SSBAR_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("SSBAR_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*CHARM_hists)[i] = new TH1D(("CHARM_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHARM_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*SIGNAL_hists)[i] = new TH1D(("Signal_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        }
    }

    fp = fopen(dirname, "r");
    for (int i = 0; i < Nentry; i++) {
        fscanf(fp, "%lf\n", &eigen_value);
        for (int j = 0; j < RarityBins * 7; j++) fscanf(fp, "%lf\n", &weight_sys[j]);

        for (int k = 0; k < RarityBins; k++) {
            (*CHG_hists)[i]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 1 * RarityBins]));
            (*UUBAR_hists)[i]->SetBinContent(k + 1, UUBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 2 * RarityBins]));
            (*DDBAR_hists)[i]->SetBinContent(k + 1, DDBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 3 * RarityBins]));
            (*SSBAR_hists)[i]->SetBinContent(k + 1, SSBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 4 * RarityBins]));
            (*CHARM_hists)[i]->SetBinContent(k + 1, CHARM_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 5 * RarityBins]));
            (*SIGNAL_hists)[i]->SetBinContent(k + 1, SIGNAL_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 6 * RarityBins]));

            (*CHG_hists)[i + Nentry]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i + Nentry]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 1 * RarityBins]));
            (*UUBAR_hists)[i + Nentry]->SetBinContent(k + 1, UUBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 2 * RarityBins]));
            (*DDBAR_hists)[i + Nentry]->SetBinContent(k + 1, DDBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 3 * RarityBins]));
            (*SSBAR_hists)[i + Nentry]->SetBinContent(k + 1, SSBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 4 * RarityBins]));
            (*CHARM_hists)[i + Nentry]->SetBinContent(k + 1, CHARM_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 5 * RarityBins]));
            (*SIGNAL_hists)[i + Nentry]->SetBinContent(k + 1, SIGNAL_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 6 * RarityBins]));
        }

    }

    fclose(fp);

    return Nentry;
}

void GetPIDUncorrelatedPDFs(const char* dirname, TH1D* CHG_hist, TH1D* MIX_hist, TH1D* UUBAR_hist, TH1D* DDBAR_hist, TH1D* SSBAR_hist, TH1D* CHARM_hist, TH1D* SIGNAL_hist) { // get shape sys histogram from txt file
    FILE* fp;
    fp = fopen(dirname, "r");

    double weight_sys[RarityBins * 7] = { 0.0 };
    for (int i = 0; i < RarityBins * 7; i++) fscanf(fp, "%lf\n", &weight_sys[i]);
    fclose(fp);

    for (int i = 0; i < RarityBins * 7; i++) weight_sys[i] = std::sqrt(weight_sys[i]);

    for (int i = 0; i < RarityBins; i++) CHG_hist->SetBinContent(i + 1, weight_sys[i]);
    for (int i = 0; i < RarityBins; i++) MIX_hist->SetBinContent(i + 1, weight_sys[RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) UUBAR_hist->SetBinContent(i + 1, weight_sys[2 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) DDBAR_hist->SetBinContent(i + 1, weight_sys[3 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) SSBAR_hist->SetBinContent(i + 1, weight_sys[4 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) CHARM_hist->SetBinContent(i + 1, weight_sys[5 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) SIGNAL_hist->SetBinContent(i + 1, weight_sys[6 * RarityBins + i]);
}

int GetBRcorrelatedPDFs(const char* dirname, TH1D* CHG_nominal_hist, TH1D* MIX_nominal_hist, TH1D* Signal_nominal_hist, TH1D*** CHG_hists, TH1D*** MIX_hists, TH1D*** Signal_hists) { // get shape sys histogram from txt file
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

    *CHG_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *MIX_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *Signal_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);

    for (int i = 0; i < Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_BR_correlated" + std::to_string(i) + "_p").c_str(), ("CHG_BR_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_BR_correlated" + std::to_string(i) + "_p").c_str(), ("MIX_BR_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists)[i] = new TH1D(("Signal_BR_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_BR_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
    }

    for (int i = Nentry; i < 2 * Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_BR_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHG_BR_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_BR_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("MIX_BR_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists)[i] = new TH1D(("Signal_BR_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_BR_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
    }

    fp = fopen(dirname, "r");
    for (int i = 0; i < Nentry; i++) {
        fscanf(fp, "%lf\n", &eigen_value);
        for (int j = 0; j < RarityBins * 3; j++) fscanf(fp, "%lf\n", &weight_sys[j]);

        for (int k = 0; k < RarityBins; k++) {
            (*CHG_hists)[i]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 1 * RarityBins]));
            (*Signal_hists)[i]->SetBinContent(k + 1, Signal_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 2 * RarityBins]));

            (*CHG_hists)[i + Nentry]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i + Nentry]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 1 * RarityBins]));
            (*Signal_hists)[i + Nentry]->SetBinContent(k + 1, Signal_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 2 * RarityBins]));
        }

    }

    fclose(fp);

    return Nentry;
}

void GetBRUncorrelatedPDFs(const char* dirname, TH1D* CHG_hist, TH1D* MIX_hist, TH1D* Signal_hist) { // get shape sys histogram from txt file
    FILE* fp;
    fp = fopen(dirname, "r");

    double weight_sys[RarityBins * 3] = { 0.0 };
    for (int i = 0; i < RarityBins * 3; i++) fscanf(fp, "%lf\n", &weight_sys[i]);
    fclose(fp);

    for (int i = 0; i < RarityBins * 3; i++) weight_sys[i] = std::sqrt(weight_sys[i]);

    for (int i = 0; i < RarityBins; i++) CHG_hist->SetBinContent(i + 1, weight_sys[i]);
    for (int i = 0; i < RarityBins; i++) MIX_hist->SetBinContent(i + 1, weight_sys[RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) Signal_hist->SetBinContent(i + 1, weight_sys[2 * RarityBins + i]);
}

int GetfractioncorrelatedPDFs(const char* dirname, TH1D* CHG_nominal_hist, TH1D* MIX_nominal_hist, TH1D* Signal_nominal_hist, TH1D*** CHG_hists, TH1D*** MIX_hists, TH1D*** Signal_hists) { // get shape sys histogram from txt file
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

    *CHG_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *MIX_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *Signal_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);

    for (int i = 0; i < Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_fraction_correlated" + std::to_string(i) + "_p").c_str(), ("CHG_fraction_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_fraction_correlated" + std::to_string(i) + "_p").c_str(), ("MIX_fraction_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists)[i] = new TH1D(("Signal_fraction_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_fraction_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
    }

    for (int i = Nentry; i < 2 * Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_fraction_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHG_fraction_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_fraction_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("MIX_fraction_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists)[i] = new TH1D(("Signal_fraction_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_fraction_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
    }

    fp = fopen(dirname, "r");
    for (int i = 0; i < Nentry; i++) {
        fscanf(fp, "%lf\n", &eigen_value);
        for (int j = 0; j < RarityBins * 3; j++) fscanf(fp, "%lf\n", &weight_sys[j]);

        for (int k = 0; k < RarityBins; k++) {
            (*CHG_hists)[i]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 1 * RarityBins]));
            (*Signal_hists)[i]->SetBinContent(k + 1, Signal_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 2 * RarityBins]));

            (*CHG_hists)[i + Nentry]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i + Nentry]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 1 * RarityBins]));
            (*Signal_hists)[i + Nentry]->SetBinContent(k + 1, Signal_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 2 * RarityBins]));
        }

    }

    fclose(fp);

    return Nentry;
}

void GetfractionUncorrelatedPDFs(const char* dirname, TH1D* CHG_hist, TH1D* MIX_hist, TH1D* Signal_hist) { // get shape sys histogram from txt file
    FILE* fp;
    fp = fopen(dirname, "r");

    double weight_sys[RarityBins * 3] = { 0.0 };
    for (int i = 0; i < RarityBins * 3; i++) fscanf(fp, "%lf\n", &weight_sys[i]);
    fclose(fp);

    for (int i = 0; i < RarityBins * 3; i++) weight_sys[i] = std::sqrt(weight_sys[i]);

    for (int i = 0; i < RarityBins; i++) CHG_hist->SetBinContent(i + 1, weight_sys[i]);
    for (int i = 0; i < RarityBins; i++) MIX_hist->SetBinContent(i + 1, weight_sys[RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) Signal_hist->SetBinContent(i + 1, weight_sys[2 * RarityBins + i]);
}

int GetFragmentationcorrelatedPDFs(const char* dirname, TH1D* CHG_nominal_hist, TH1D* MIX_nominal_hist, TH1D* Signal_nominal_hist, TH1D*** CHG_hists, TH1D*** MIX_hists, TH1D*** Signal_hists) { // get shape sys histogram from txt file
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

    *CHG_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *MIX_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *Signal_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);

    for (int i = 0; i < Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_Fragmentation_correlated" + std::to_string(i) + "_p").c_str(), ("CHG_Fragmentation_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_Fragmentation_correlated" + std::to_string(i) + "_p").c_str(), ("MIX_Fragmentation_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists)[i] = new TH1D(("Signal_Fragmentation_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_Fragmentation_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
    }

    for (int i = Nentry; i < 2 * Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_Fragmentation_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHG_Fragmentation_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_Fragmentation_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("MIX_Fragmentation_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists)[i] = new TH1D(("Signal_Fragmentation_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_Fragmentation_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
    }

    fp = fopen(dirname, "r");
    for (int i = 0; i < Nentry; i++) {
        fscanf(fp, "%lf\n", &eigen_value);
        for (int j = 0; j < RarityBins * 3; j++) fscanf(fp, "%lf\n", &weight_sys[j]);

        for (int k = 0; k < RarityBins; k++) {
            (*CHG_hists)[i]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 1 * RarityBins]));
            (*Signal_hists)[i]->SetBinContent(k + 1, Signal_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 2 * RarityBins]));

            (*CHG_hists)[i + Nentry]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i + Nentry]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 1 * RarityBins]));
            (*Signal_hists)[i + Nentry]->SetBinContent(k + 1, Signal_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 2 * RarityBins]));
        }

    }

    fclose(fp);

    return Nentry;
}

void GetFragmentationUncorrelatedPDFs(const char* dirname, TH1D* CHG_hist, TH1D* MIX_hist, TH1D* Signal_hist) { // get shape sys histogram from txt file
    FILE* fp;
    fp = fopen(dirname, "r");

    double weight_sys[RarityBins * 3] = { 0.0 };
    for (int i = 0; i < RarityBins * 3; i++) fscanf(fp, "%lf\n", &weight_sys[i]);
    fclose(fp);

    for (int i = 0; i < RarityBins * 3; i++) weight_sys[i] = std::sqrt(weight_sys[i]);

    for (int i = 0; i < RarityBins; i++) CHG_hist->SetBinContent(i + 1, weight_sys[i]);
    for (int i = 0; i < RarityBins; i++) MIX_hist->SetBinContent(i + 1, weight_sys[RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) Signal_hist->SetBinContent(i + 1, weight_sys[2 * RarityBins + i]);
}

int Getpi0correlatedPDFs(const char* dirname, TH1D* CHG_nominal_hist, TH1D* MIX_nominal_hist, TH1D* UUBAR_nominal_hist, TH1D* DDBAR_nominal_hist, TH1D* SSBAR_nominal_hist, TH1D* CHARM_nominal_hist, TH1D* SIGNAL_nominal_hist, TH1D*** CHG_hists, TH1D*** MIX_hists, TH1D*** UUBAR_hists, TH1D*** DDBAR_hists, TH1D*** SSBAR_hists, TH1D*** CHARM_hists, TH1D*** SIGNAL_hists) { // get shape sys histogram from txt file
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

    *CHG_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *MIX_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *UUBAR_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *DDBAR_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *SSBAR_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *CHARM_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *SIGNAL_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);

    for (int i = 0; i < Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("CHG_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("MIX_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*UUBAR_hists)[i] = new TH1D(("UUBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("UUBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*DDBAR_hists)[i] = new TH1D(("DDBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("DDBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*SSBAR_hists)[i] = new TH1D(("SSBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("SSBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*CHARM_hists)[i] = new TH1D(("CHARM_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("CHARM_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*SIGNAL_hists)[i] = new TH1D(("Signal_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
    }

    for (int i = Nentry; i < 2 * Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHG_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("MIX_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*UUBAR_hists)[i] = new TH1D(("UUBAR_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("UUBAR_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*DDBAR_hists)[i] = new TH1D(("DDBAR_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("DDBAR_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*SSBAR_hists)[i] = new TH1D(("SSBAR_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("SSBAR_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*CHARM_hists)[i] = new TH1D(("CHARM_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHARM_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*SIGNAL_hists)[i] = new TH1D(("Signal_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
    }

    fp = fopen(dirname, "r");
    for (int i = 0; i < Nentry; i++) {
        fscanf(fp, "%lf\n", &eigen_value);
        for (int j = 0; j < RarityBins * 7; j++) fscanf(fp, "%lf\n", &weight_sys[j]);

        for (int k = 0; k < RarityBins; k++) {
            (*CHG_hists)[i]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 1 * RarityBins]));
            (*UUBAR_hists)[i]->SetBinContent(k + 1, UUBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 2 * RarityBins]));
            (*DDBAR_hists)[i]->SetBinContent(k + 1, DDBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 3 * RarityBins]));
            (*SSBAR_hists)[i]->SetBinContent(k + 1, SSBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 4 * RarityBins]));
            (*CHARM_hists)[i]->SetBinContent(k + 1, CHARM_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 5 * RarityBins]));
            (*SIGNAL_hists)[i]->SetBinContent(k + 1, SIGNAL_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 6 * RarityBins]));

            (*CHG_hists)[i + Nentry]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i + Nentry]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 1 * RarityBins]));
            (*UUBAR_hists)[i + Nentry]->SetBinContent(k + 1, UUBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 2 * RarityBins]));
            (*DDBAR_hists)[i + Nentry]->SetBinContent(k + 1, DDBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 3 * RarityBins]));
            (*SSBAR_hists)[i + Nentry]->SetBinContent(k + 1, SSBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 4 * RarityBins]));
            (*CHARM_hists)[i + Nentry]->SetBinContent(k + 1, CHARM_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 5 * RarityBins]));
            (*SIGNAL_hists)[i + Nentry]->SetBinContent(k + 1, SIGNAL_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 6 * RarityBins]));
        }

    }

    fclose(fp);

    return Nentry;
}

void Getpi0UncorrelatedPDFs(const char* dirname, TH1D* CHG_hist, TH1D* MIX_hist, TH1D* UUBAR_hist, TH1D* DDBAR_hist, TH1D* SSBAR_hist, TH1D* CHARM_hist, TH1D* SIGNAL_hist) { // get shape sys histogram from txt file
    FILE* fp;
    fp = fopen(dirname, "r");

    double weight_sys[RarityBins * 7] = { 0.0 };
    for (int i = 0; i < RarityBins * 7; i++) fscanf(fp, "%lf\n", &weight_sys[i]);
    fclose(fp);

    for (int i = 0; i < RarityBins * 7; i++) weight_sys[i] = std::sqrt(weight_sys[i]);

    for (int i = 0; i < RarityBins; i++) CHG_hist->SetBinContent(i + 1, weight_sys[i]);
    for (int i = 0; i < RarityBins; i++) MIX_hist->SetBinContent(i + 1, weight_sys[RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) UUBAR_hist->SetBinContent(i + 1, weight_sys[2 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) DDBAR_hist->SetBinContent(i + 1, weight_sys[3 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) SSBAR_hist->SetBinContent(i + 1, weight_sys[4 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) CHARM_hist->SetBinContent(i + 1, weight_sys[5 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) SIGNAL_hist->SetBinContent(i + 1, weight_sys[6 * RarityBins + i]);
}

void GetKffPDFs(const char* dirname, const char* included_string, TH1D* hist[7], double Correction_factor_BR[7], const char* type, int charge, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get Kff uncertainty PDF with appropriate correction
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

    const double alpha0 = 0.2545; // 0.2545 +- 0.0090
    const double alpha1 = -0.71; // -0.71 +- 0.14
    const double alpha2 = 0.32; // 0.32 +- 0.59

    /*
    <PCA>
    e0 = 0.00090 alpha0 + 0.19987 alpha1 + 0.97982 alpha2
    e1 = 0.06182 alpha0 + 0.97794 alpha1 - 0.19954 alpha2
    e2 = 0.99809 alpha0 - 0.06075 alpha2 + 0.01148 alpha2
    lambda0 = 0.60197
    lambda1 = 0.07314
    lambda2 = 0.00778
    */
    const double LinearCoefficients[3][3] = {
        {0.000898287544607,   0.199871550874459,   0.979821696141969},
        {0.061817247836003,   0.977937066768286, - 0.199543782941545},
        {0.998087080844888, - 0.060749128320321,   0.011477040509285}
    };
    const double Lambdas[7][3] = {
        { 0.0,                0.0,                0.0              },
        { 0.601972208713257,  0.0,                0.0              },
        {-0.601972208713257,  0.0,                0.0              },
        { 0.0,                0.073143483537865,  0.0              },
        { 0.0,               -0.073143483537865,  0.0              },
        { 0.0,                0.0,                0.007777700923050},
        { 0.0,                0.0,               -0.007777700923050}
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
    fluctuations[1][0] = 0.60197 * 0.00090; // alpha0 change Lambdas[1][0]*LinearCoefficients[0][0] + Lambdas[1][1]*LinearCoefficients[1][0] + Lambdas[1][2]*LinearCoefficients[2][0]
    fluctuations[1][1] = 0.60197 * 0.19987; // alpha1 change Lambdas[1][0]*LinearCoefficients[0][1] + Lambdas[1][1]*LinearCoefficients[1][1] + Lambdas[1][2]*LinearCoefficients[2][1]
    fluctuations[1][2] = 0.60197 * 0.97982; // alpha2 change
    fluctuations[2][0] = (-0.60197) * 0.00090; // alpha0 change
    fluctuations[2][1] = (-0.60197) * 0.19987; // alpha1 change
    fluctuations[2][2] = (-0.60197) * 0.97982; // alpha2 change
    fluctuations[3][0] = 0.07314 * 0.06182; // alpha0 change
    fluctuations[3][1] = 0.07314 * 0.97794; // alpha1 change
    fluctuations[3][2] = 0.07314 * (-0.19954); // alpha2 change
    fluctuations[4][0] // alpha0 change
    fluctuations[4][1] // alpha1 change
    fluctuations[4][2] // alpha2 change
    fluctuations[5][0] // alpha0 change
    fluctuations[5][1] // alpha1 change
    fluctuations[5][2] // alpha2 change
    fluctuations[6][0] // alpha0 change
    fluctuations[6][1] // alpha1 change
    fluctuations[6][2] // alpha2 change
    */
    double value[7] = { 0.0 }; // value of lambda^1.5 * fp*fp
    double Nevts[7] = { 0.0 }; // number of events at each fluctuations
    double tot_value[7] = { 0.0 }; // value of integral lambda^1.5 * fp*fp

    double q2 = -1;

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double invM_Xnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;
    double N_Xplusnn = 0;
    double N_Xzeronn = 0;

    double Ngamma_v200 = -1;

    double s13_KpKLKL = -1;
    double s23_KpKLKL = -1;
    double nB2KpKLKL_all_KpKLKL = -1;
    double nB2KpKLKL_NR_KpKLKL = -1;

    double s13_KSKLKL = -1;
    double s23_KSKLKL = -1;
    double s12_KSKLKL = -1;
    double nB2KSKLKL_all_KSKLKL = -1;
    double nB2KSKLKL_NR_KSKLKL = -1;

    double nKL_XKLKL = -1;
    double XKLKL_E_1st = -1;
    double XKLKL_px_1st = -1;
    double XKLKL_py_1st = -1;
    double XKLKL_pz_1st = -1;
    double XKLKL_E_2nd = -1;
    double XKLKL_px_2nd = -1;
    double XKLKL_py_2nd = -1;
    double XKLKL_pz_2nd = -1;
    double nB2KstarKLKL = -1;
    double nB02KstarKLKL = -1;

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    double Bsig_M = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);
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
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s13_KpKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s23_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_all__bc", &nB2KpKLKL_all_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_NR__bc", &nB2KpKLKL_NR_KpKLKL);

        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s13_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &s23_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s12_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_all__bc", &nB2KSKLKL_all_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_NR__bc", &nB2KSKLKL_NR_KSKLKL);

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

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

            // select B+ --> K+ nu nubar reconstruction only
            if ((std::abs(Upsilon_ID) < MyEPSILON) && (std::abs(Bsig_ID) < MyEPSILON)) {}
            else continue;

            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE), temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE), temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(0, i_fake, MCTYPE), temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(1, i_fake, MCTYPE), temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(2, i_fake, MCTYPE), temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(3, i_fake, MCTYPE), temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(0, i_fake, MCTYPE), temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(1, i_fake, MCTYPE), temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(2, i_fake, MCTYPE), temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(3, i_fake, MCTYPE), temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B --> K* KL KL correction factor
            double Correction_KstarKLKL = corrector_KstarKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KstarKLKL + nB02KstarKLKL);

            // B --> X KL KL correction factor
            double Correction_XKLKL = corrector_XsKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KpKLKL_all_KpKLKL, nB2KSKLKL_all_KSKLKL, nB2KstarKLKL + nB02KstarKLKL, nKL_XKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            q2 = q2 * q2;

            for (int k = 0; k < 7; k++) {
                const double alpha0_fluc = alpha0 + fluctuations[k][0];
                const double alpha1_fluc = alpha1 + fluctuations[k][1];
                const double alpha2_fluc = alpha2 + fluctuations[k][2];

                double mp = 5.4158;
                double tp = (m_b + m_k) * (m_b + m_k);
                double tm = (m_b - m_k) * (m_b - m_k);
                double t0 = 0.0;
                double z = (sqrt(tp - q2) - sqrt(tp - t0)) / (sqrt(tp - q2) + sqrt(tp - t0));
                double fp = (1 / (1 - q2 / (mp * mp))) * (alpha0_fluc + alpha1_fluc * z + alpha2_fluc * z * z + (-alpha1_fluc + 2 * alpha2_fluc) * z * z * z / 3);
                double lambda = (m_b * m_b * m_b * m_b) + (m_k * m_k * m_k * m_k) + (q2 * q2) - 2 * (m_b * m_b * m_k * m_k + m_b * m_b * q2 + m_k * m_k * q2);

                value[k] = std::pow(lambda, 1.5) * fp * fp;
                double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * (value[k] / value[0]);
                if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(q2, "Bplus");
                else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(q2, "Bzero");
                FillTemplate(hist[k], MVA_var, total_weight, Bsig_M);
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

            double mp = 5.4158;
            double tp = (m_b_fix + m_k_fix) * (m_b_fix + m_k_fix);
            double tm = (m_b_fix - m_k_fix) * (m_b_fix - m_k_fix);
            double t0 = 0.0;
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

double GetKstarffPDFs(const char* dirname, const char* included_string, TH1D* hist[19], double Correction_factor_BR[19], const char* type, int charge, double weight_var = 1.0) { // get Kstarff uncertainty PDF with appropriate correction
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
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double invM_Xnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;
    double N_Xplusnn = 0;
    double N_Xzeronn = 0;

    double Ngamma_v200 = -1;

    double s13_KpKLKL = -1;
    double s23_KpKLKL = -1;
    double nB2KpKLKL_all_KpKLKL = -1;
    double nB2KpKLKL_NR_KpKLKL = -1;

    double s13_KSKLKL = -1;
    double s23_KSKLKL = -1;
    double s12_KSKLKL = -1;
    double nB2KSKLKL_all_KSKLKL = -1;
    double nB2KSKLKL_NR_KSKLKL = -1;

    double nKL_XKLKL = -1;
    double XKLKL_E_1st = -1;
    double XKLKL_px_1st = -1;
    double XKLKL_py_1st = -1;
    double XKLKL_pz_1st = -1;
    double XKLKL_E_2nd = -1;
    double XKLKL_px_2nd = -1;
    double XKLKL_py_2nd = -1;
    double XKLKL_pz_2nd = -1;
    double nB2KstarKLKL = -1;
    double nB02KstarKLKL = -1;

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    double Bsig_M = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);
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
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s13_KpKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s23_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_all__bc", &nB2KpKLKL_all_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_NR__bc", &nB2KpKLKL_NR_KpKLKL);

        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s13_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &s23_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s12_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_all__bc", &nB2KSKLKL_all_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_NR__bc", &nB2KSKLKL_NR_KSKLKL);

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

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

            // select B+ --> K+ nu nubar reconstruction only
            if ((std::abs(Upsilon_ID) < MyEPSILON) && (std::abs(Bsig_ID) < MyEPSILON)) {}
            else continue;

            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE), temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE), temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(0, i_fake, MCTYPE), temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(1, i_fake, MCTYPE), temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(2, i_fake, MCTYPE), temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(3, i_fake, MCTYPE), temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(0, i_fake, MCTYPE), temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(1, i_fake, MCTYPE), temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(2, i_fake, MCTYPE), temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(3, i_fake, MCTYPE), temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B --> K* KL KL correction factor
            double Correction_KstarKLKL = corrector_KstarKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KstarKLKL + nB02KstarKLKL);

            // B --> X KL KL correction factor
            double Correction_XKLKL = corrector_XsKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KpKLKL_all_KpKLKL, nB2KSKLKL_all_KSKLKL, nB2KstarKLKL + nB02KstarKLKL, nKL_XKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

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
                double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn* Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * (value[k] / value[0]);
                if (q2 < MyEPSILON) total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL; // Makeshift
                FillTemplate(hist[k], MVA_var, total_weight, Bsig_M);
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

double GetOLDKffPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get PDF with old Kff from NominalPDF
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
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
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double invM_Xnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;
    double N_Xplusnn = 0;
    double N_Xzeronn = 0;

    double Ngamma_v200 = -1;

    double s13_KpKLKL = -1;
    double s23_KpKLKL = -1;
    double nB2KpKLKL_all_KpKLKL = -1;
    double nB2KpKLKL_NR_KpKLKL = -1;

    double s13_KSKLKL = -1;
    double s23_KSKLKL = -1;
    double s12_KSKLKL = -1;
    double nB2KSKLKL_all_KSKLKL = -1;
    double nB2KSKLKL_NR_KSKLKL = -1;

    double nKL_XKLKL = -1;
    double XKLKL_E_1st = -1;
    double XKLKL_px_1st = -1;
    double XKLKL_py_1st = -1;
    double XKLKL_pz_1st = -1;
    double XKLKL_E_2nd = -1;
    double XKLKL_px_2nd = -1;
    double XKLKL_py_2nd = -1;
    double XKLKL_pz_2nd = -1;
    double nB2KstarKLKL = -1;
    double nB02KstarKLKL = -1;

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    double Bsig_M = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s13_KpKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s23_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_all__bc", &nB2KpKLKL_all_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_NR__bc", &nB2KpKLKL_NR_KpKLKL);

        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s13_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &s23_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s12_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_all__bc", &nB2KSKLKL_all_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_NR__bc", &nB2KSKLKL_NR_KSKLKL);

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            // select B+ --> K+ nu nubar reconstruction only
            if ((std::abs(Upsilon_ID) < MyEPSILON) && (std::abs(Bsig_ID) < MyEPSILON)) {}
            else continue;

            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE), temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE), temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(0, i_fake, MCTYPE), temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(1, i_fake, MCTYPE), temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(2, i_fake, MCTYPE), temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(3, i_fake, MCTYPE), temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(0, i_fake, MCTYPE), temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(1, i_fake, MCTYPE), temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(2, i_fake, MCTYPE), temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(3, i_fake, MCTYPE), temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B --> K* KL KL correction factor
            double Correction_KstarKLKL = corrector_KstarKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KstarKLKL + nB02KstarKLKL);

            // B --> X KL KL correction factor
            double Correction_XKLKL = corrector_XsKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KpKLKL_all_KpKLKL, nB2KSKLKL_all_KSKLKL, nB2KstarKLKL + nB02KstarKLKL, nKL_XKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * (BR_Kplus_nunubar_OLD / BR_Kplus_nunubar);
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * (BR_K0_nunubar_OLD / BR_K0_nunubar);
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE) * (BR_Xsu_nonresonant_nunubar_OLD / BR_Xsu_nonresonant_nunubar);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE) * (BR_Xsd_nonresonant_nunubar_OLD / BR_Xsd_nonresonant_nunubar);

            Nevt = Nevt + FillTemplate(hist, MVA_var, total_weight, Bsig_M);

        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double GetNevtWithBDTc(const char* dirname, const char* included_string, const char* type, const char* sample, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get the number of event after BDTc correction
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    float MVA_var = 0;
    float BDTc_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double invM_Xnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;
    double N_Xplusnn = 0;
    double N_Xzeronn = 0;

    double Ngamma_v200 = -1;

    double s13_KpKLKL = -1;
    double s23_KpKLKL = -1;
    double nB2KpKLKL_all_KpKLKL = -1;
    double nB2KpKLKL_NR_KpKLKL = -1;

    double s13_KSKLKL = -1;
    double s23_KSKLKL = -1;
    double s12_KSKLKL = -1;
    double nB2KSKLKL_all_KSKLKL = -1;
    double nB2KSKLKL_NR_KSKLKL = -1;

    double nKL_XKLKL = -1;
    double XKLKL_E_1st = -1;
    double XKLKL_px_1st = -1;
    double XKLKL_py_1st = -1;
    double XKLKL_pz_1st = -1;
    double XKLKL_E_2nd = -1;
    double XKLKL_px_2nd = -1;
    double XKLKL_py_2nd = -1;
    double XKLKL_pz_2nd = -1;
    double nB2KstarKLKL = -1;
    double nB02KstarKLKL = -1;

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    double Bsig_M = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s13_KpKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s23_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_all__bc", &nB2KpKLKL_all_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_NR__bc", &nB2KpKLKL_NR_KpKLKL);

        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s13_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &s23_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s12_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_all__bc", &nB2KSKLKL_all_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_NR__bc", &nB2KSKLKL_NR_KSKLKL);

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        tree_upsilon->SetBranchAddress("MVA_Continuum", &BDTc_var);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            // select B+ --> K+ nu nubar reconstruction only
            if ((std::abs(Upsilon_ID) < MyEPSILON) && (std::abs(Bsig_ID) < MyEPSILON)) {}
            else continue;

            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE), temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE), temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(0, i_fake, MCTYPE), temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(1, i_fake, MCTYPE), temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(2, i_fake, MCTYPE), temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(3, i_fake, MCTYPE), temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(0, i_fake, MCTYPE), temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(1, i_fake, MCTYPE), temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(2, i_fake, MCTYPE), temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(3, i_fake, MCTYPE), temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B --> K* KL KL correction factor
            double Correction_KstarKLKL = corrector_KstarKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KstarKLKL + nB02KstarKLKL);

            // B --> X KL KL correction factor
            double Correction_XKLKL = corrector_XsKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KpKLKL_all_KpKLKL, nB2KSKLKL_all_KSKLKL, nB2KstarKLKL + nB02KstarKLKL, nKL_XKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double BDTc_weight = BDTcToWeight(BDTc_var);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * BDTc_weight;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            Nevt = Nevt + DoesNotFillTemplate(hist, MVA_var, total_weight, Bsig_M);

        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double GetBDTcPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get BDTc PDF with appropriate correction
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    float MVA_var = 0;
    float BDTc_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double invM_Xnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;
    double N_Xplusnn = 0;
    double N_Xzeronn = 0;

    double Ngamma_v200 = -1;

    double s13_KpKLKL = -1;
    double s23_KpKLKL = -1;
    double nB2KpKLKL_all_KpKLKL = -1;
    double nB2KpKLKL_NR_KpKLKL = -1;

    double s13_KSKLKL = -1;
    double s23_KSKLKL = -1;
    double s12_KSKLKL = -1;
    double nB2KSKLKL_all_KSKLKL = -1;
    double nB2KSKLKL_NR_KSKLKL = -1;

    double nKL_XKLKL = -1;
    double XKLKL_E_1st = -1;
    double XKLKL_px_1st = -1;
    double XKLKL_py_1st = -1;
    double XKLKL_pz_1st = -1;
    double XKLKL_E_2nd = -1;
    double XKLKL_px_2nd = -1;
    double XKLKL_py_2nd = -1;
    double XKLKL_pz_2nd = -1;
    double nB2KstarKLKL = -1;
    double nB02KstarKLKL = -1;

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    double Bsig_M = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s13_KpKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s23_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_all__bc", &nB2KpKLKL_all_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_NR__bc", &nB2KpKLKL_NR_KpKLKL);

        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s13_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &s23_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s12_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_all__bc", &nB2KSKLKL_all_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_NR__bc", &nB2KSKLKL_NR_KSKLKL);

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        tree_upsilon->SetBranchAddress("MVA_Continuum", &BDTc_var);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // apply BDTc weight with appropritate norm factor
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            // select B+ --> K+ nu nubar reconstruction only
            if ((std::abs(Upsilon_ID) < MyEPSILON) && (std::abs(Bsig_ID) < MyEPSILON)) {}
            else continue;

            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE), temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE), temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(0, i_fake, MCTYPE), temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(1, i_fake, MCTYPE), temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(2, i_fake, MCTYPE), temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(3, i_fake, MCTYPE), temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(0, i_fake, MCTYPE), temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(1, i_fake, MCTYPE), temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(2, i_fake, MCTYPE), temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(3, i_fake, MCTYPE), temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B --> K* KL KL correction factor
            double Correction_KstarKLKL = corrector_KstarKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KstarKLKL + nB02KstarKLKL);

            // B --> X KL KL correction factor
            double Correction_XKLKL = corrector_XsKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KpKLKL_all_KpKLKL, nB2KSKLKL_all_KSKLKL, nB2KstarKLKL + nB02KstarKLKL, nKL_XKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double BDTc_weight = BDTcToWeight(BDTc_var);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * BDTc_weight;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            Nevt = Nevt + FillTemplate(hist, MVA_var, total_weight, Bsig_M);

        }

        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

void GetNegativeChangePDFs(TH1D* nominal_hist, TH1D* positive_hist, TH1D* negative_hist) {
    for (int i = 0; i < RarityBins; i++) {
        double nominal_value = nominal_hist->GetBinContent(i + 1);
        double positive_value = positive_hist->GetBinContent(i + 1);
        double deviation = positive_value - nominal_value;

        negative_hist->SetBinContent(i + 1, nominal_value - deviation);
    }
}

double GetmbPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, bool IsIthigh, bool IsItXsu, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // apply correction for mb weight files
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {
        printf("[ERROR] unexpected type name\n");
        printf("[ERROR] Continuum type cannot be selected for mb pdf\n");
        exit(1);
    }
    else {
        printf("[ERROR] unexpected type name\n");
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

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double invM_Xnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;
    double N_Xplusnn = 0;
    double N_Xzeronn = 0;

    double Ngamma_v200 = -1;

    double s13_KpKLKL = -1;
    double s23_KpKLKL = -1;
    double nB2KpKLKL_all_KpKLKL = -1;
    double nB2KpKLKL_NR_KpKLKL = -1;

    double s13_KSKLKL = -1;
    double s23_KSKLKL = -1;
    double s12_KSKLKL = -1;
    double nB2KSKLKL_all_KSKLKL = -1;
    double nB2KSKLKL_NR_KSKLKL = -1;

    double nKL_XKLKL = -1;
    double XKLKL_E_1st = -1;
    double XKLKL_px_1st = -1;
    double XKLKL_py_1st = -1;
    double XKLKL_pz_1st = -1;
    double XKLKL_E_2nd = -1;
    double XKLKL_px_2nd = -1;
    double XKLKL_py_2nd = -1;
    double XKLKL_pz_2nd = -1;
    double nB2KstarKLKL = -1;
    double nB02KstarKLKL = -1;

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    double Bsig_M = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s13_KpKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s23_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_all__bc", &nB2KpKLKL_all_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_NR__bc", &nB2KpKLKL_NR_KpKLKL);

        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s13_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &s23_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s12_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_all__bc", &nB2KSKLKL_all_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_NR__bc", &nB2KSKLKL_NR_KSKLKL);

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            // select B+ --> K+ nu nubar reconstruction only
            if ((std::abs(Upsilon_ID) < MyEPSILON) && (std::abs(Bsig_ID) < MyEPSILON)) {}
            else continue;

            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE), temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE), temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(0, i_fake, MCTYPE), temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(1, i_fake, MCTYPE), temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(2, i_fake, MCTYPE), temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(3, i_fake, MCTYPE), temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(0, i_fake, MCTYPE), temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(1, i_fake, MCTYPE), temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(2, i_fake, MCTYPE), temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(3, i_fake, MCTYPE), temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B --> K* KL KL correction factor
            double Correction_KstarKLKL = corrector_KstarKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KstarKLKL + nB02KstarKLKL);

            // B --> X KL KL correction factor
            double Correction_XKLKL = corrector_XsKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KpKLKL_all_KpKLKL, nB2KSKLKL_all_KSKLKL, nB2KstarKLKL + nB02KstarKLKL, nKL_XKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            if (IsIthigh == true && IsItXsu == true) total_weight = total_weight * ReadWeightHist(Xsu_Hmb_weight, invM * invM);
            else if (IsIthigh == false && IsItXsu == true) total_weight = total_weight * ReadWeightHist(Xsu_Lmb_weight, invM * invM);
            else if (IsIthigh == true && IsItXsu == false) total_weight = total_weight * ReadWeightHist(Xsd_Hmb_weight, invM * invM);
            else if (IsIthigh == false && IsItXsu == false) total_weight = total_weight * ReadWeightHist(Xsd_Lmb_weight, invM * invM);

            Nevt = Nevt + FillTemplate(hist, MVA_var, total_weight, Bsig_M);

        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double GetpfPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, bool IsIthigh, bool IsItXsu, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // apply correction for pf weight files
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {
        printf("[ERROR] unexpected type name\n");
        printf("[ERROR] Continuum type cannot be selected for pf pdf\n");
        exit(1);
    }
    else {
        printf("[ERROR] unexpected type name\n");
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

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;
    double MXs = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double invM_Xnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;
    double N_Xplusnn = 0;
    double N_Xzeronn = 0;

    double Ngamma_v200 = -1;

    double s13_KpKLKL = -1;
    double s23_KpKLKL = -1;
    double nB2KpKLKL_all_KpKLKL = -1;
    double nB2KpKLKL_NR_KpKLKL = -1;

    double s13_KSKLKL = -1;
    double s23_KSKLKL = -1;
    double s12_KSKLKL = -1;
    double nB2KSKLKL_all_KSKLKL = -1;
    double nB2KSKLKL_NR_KSKLKL = -1;

    double nKL_XKLKL = -1;
    double XKLKL_E_1st = -1;
    double XKLKL_px_1st = -1;
    double XKLKL_py_1st = -1;
    double XKLKL_pz_1st = -1;
    double XKLKL_E_2nd = -1;
    double XKLKL_px_2nd = -1;
    double XKLKL_py_2nd = -1;
    double XKLKL_pz_2nd = -1;
    double nB2KstarKLKL = -1;
    double nB02KstarKLKL = -1;

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    double Bsig_M = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            if (IsItXsu) tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &MXs);
            else tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &MXs);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s13_KpKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s23_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_all__bc", &nB2KpKLKL_all_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_NR__bc", &nB2KpKLKL_NR_KpKLKL);

        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s13_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &s23_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s12_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_all__bc", &nB2KSKLKL_all_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_NR__bc", &nB2KSKLKL_NR_KSKLKL);

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            // select B+ --> K+ nu nubar reconstruction only
            if ((std::abs(Upsilon_ID) < MyEPSILON) && (std::abs(Bsig_ID) < MyEPSILON)) {}
            else continue;

            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE), temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE), temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(0, i_fake, MCTYPE), temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(1, i_fake, MCTYPE), temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(2, i_fake, MCTYPE), temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(3, i_fake, MCTYPE), temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(0, i_fake, MCTYPE), temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(1, i_fake, MCTYPE), temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(2, i_fake, MCTYPE), temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(3, i_fake, MCTYPE), temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B --> K* KL KL correction factor
            double Correction_KstarKLKL = corrector_KstarKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KstarKLKL + nB02KstarKLKL);

            // B --> X KL KL correction factor
            double Correction_XKLKL = corrector_XsKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KpKLKL_all_KpKLKL, nB2KSKLKL_all_KSKLKL, nB2KstarKLKL + nB02KstarKLKL, nKL_XKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, MXs, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, MXs, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            if (IsIthigh == true && IsItXsu == true) total_weight = total_weight * ReadWeightHist(Xsu_Hpf_weight, MXs);
            else if (IsIthigh == false && IsItXsu == true) total_weight = total_weight * ReadWeightHist(Xsu_Lpf_weight, MXs);
            else if (IsIthigh == true && IsItXsu == false) total_weight = total_weight * ReadWeightHist(Xsd_Hpf_weight, MXs);
            else if (IsIthigh == false && IsItXsu == false) total_weight = total_weight * ReadWeightHist(Xsd_Lpf_weight, MXs);

            Nevt = Nevt + FillTemplate(hist, MVA_var, total_weight, Bsig_M);

        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double GetTransitionPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, bool IsIthigh, bool IsItXsu, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // apply correction for transition weight files
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
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
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double invM_Xnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;
    double N_Xplusnn = 0;
    double N_Xzeronn = 0;

    double Ngamma_v200 = -1;

    double s13_KpKLKL = -1;
    double s23_KpKLKL = -1;
    double nB2KpKLKL_all_KpKLKL = -1;
    double nB2KpKLKL_NR_KpKLKL = -1;

    double s13_KSKLKL = -1;
    double s23_KSKLKL = -1;
    double s12_KSKLKL = -1;
    double nB2KSKLKL_all_KSKLKL = -1;
    double nB2KSKLKL_NR_KSKLKL = -1;

    double nKL_XKLKL = -1;
    double XKLKL_E_1st = -1;
    double XKLKL_px_1st = -1;
    double XKLKL_py_1st = -1;
    double XKLKL_pz_1st = -1;
    double XKLKL_E_2nd = -1;
    double XKLKL_px_2nd = -1;
    double XKLKL_py_2nd = -1;
    double XKLKL_pz_2nd = -1;
    double nB2KstarKLKL = -1;
    double nB02KstarKLKL = -1;

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    double Bsig_M = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s13_KpKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s23_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_all__bc", &nB2KpKLKL_all_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_NR__bc", &nB2KpKLKL_NR_KpKLKL);

        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s13_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &s23_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s12_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_all__bc", &nB2KSKLKL_all_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_NR__bc", &nB2KSKLKL_NR_KSKLKL);

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            // select B+ --> K+ nu nubar reconstruction only
            if ((std::abs(Upsilon_ID) < MyEPSILON) && (std::abs(Bsig_ID) < MyEPSILON)) {}
            else continue;

            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE), temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE), temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(0, i_fake, MCTYPE), temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(1, i_fake, MCTYPE), temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(2, i_fake, MCTYPE), temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(3, i_fake, MCTYPE), temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(0, i_fake, MCTYPE), temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(1, i_fake, MCTYPE), temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(2, i_fake, MCTYPE), temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(3, i_fake, MCTYPE), temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B --> K* KL KL correction factor
            double Correction_KstarKLKL = corrector_KstarKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KstarKLKL + nB02KstarKLKL);

            // B --> X KL KL correction factor
            double Correction_XKLKL = corrector_XsKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KpKLKL_all_KpKLKL, nB2KSKLKL_all_KSKLKL, nB2KstarKLKL + nB02KstarKLKL, nKL_XKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            if (IsIthigh == true && IsItXsu == true) total_weight = total_weight * ReadWeightHist(Xsu_Htransition_weight, MVA_var);
            else if (IsIthigh == false && IsItXsu == true) total_weight = total_weight * ReadWeightHist(Xsu_Ltransition_weight, MVA_var);
            else if (IsIthigh == true && IsItXsu == false) total_weight = total_weight * ReadWeightHist(Xsd_Htransition_weight, MVA_var);
            else if (IsIthigh == false && IsItXsu == false) total_weight = total_weight * ReadWeightHist(Xsd_Ltransition_weight, MVA_var);

            Nevt = Nevt + FillTemplate(hist, MVA_var, total_weight, Bsig_M);

        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double GetKstardeltaPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, bool IsItXsu, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // apply correction for pf weight files
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    otherwise
    */
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {
        printf("[ERROR] unexpected type name\n");
        printf("[ERROR] Continuum type cannot be selected for pf pdf\n");
        exit(1);
    }
    else {
        printf("[ERROR] unexpected type name\n");
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
    double MXs = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double invM_Xnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;
    double N_Xplusnn = 0;
    double N_Xzeronn = 0;

    double Ngamma_v200 = -1;

    double s13_KpKLKL = -1;
    double s23_KpKLKL = -1;
    double nB2KpKLKL_all_KpKLKL = -1;
    double nB2KpKLKL_NR_KpKLKL = -1;

    double s13_KSKLKL = -1;
    double s23_KSKLKL = -1;
    double s12_KSKLKL = -1;
    double nB2KSKLKL_all_KSKLKL = -1;
    double nB2KSKLKL_NR_KSKLKL = -1;

    double nKL_XKLKL = -1;
    double XKLKL_E_1st = -1;
    double XKLKL_px_1st = -1;
    double XKLKL_py_1st = -1;
    double XKLKL_pz_1st = -1;
    double XKLKL_E_2nd = -1;
    double XKLKL_px_2nd = -1;
    double XKLKL_py_2nd = -1;
    double XKLKL_pz_2nd = -1;
    double nB2KstarKLKL = -1;
    double nB02KstarKLKL = -1;

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    double Bsig_M = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");
        TTree* tree_Xs = (TTree*)input_file->Get("Xs");

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);
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
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s13_KpKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s23_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_all__bc", &nB2KpKLKL_all_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_NR__bc", &nB2KpKLKL_NR_KpKLKL);

        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s13_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &s23_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s12_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_all__bc", &nB2KSKLKL_all_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_NR__bc", &nB2KSKLKL_NR_KSKLKL);

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        if (IsItXsu) tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &MXs);
        else tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &MXs);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            tree_Xs->GetEntry(j);

            // select B+ --> K+ nu nubar reconstruction only
            if ((std::abs(Upsilon_ID) < MyEPSILON) && (std::abs(Bsig_ID) < MyEPSILON)) {}
            else continue;

            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE), temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE), temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(0, i_fake, MCTYPE), temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(1, i_fake, MCTYPE), temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(2, i_fake, MCTYPE), temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(3, i_fake, MCTYPE), temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(0, i_fake, MCTYPE), temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(1, i_fake, MCTYPE), temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(2, i_fake, MCTYPE), temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(3, i_fake, MCTYPE), temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B --> K* KL KL correction factor
            double Correction_KstarKLKL = corrector_KstarKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KstarKLKL + nB02KstarKLKL);

            // B --> X KL KL correction factor
            double Correction_XKLKL = corrector_XsKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KpKLKL_all_KpKLKL, nB2KSKLKL_all_KSKLKL, nB2KstarKLKL + nB02KstarKLKL, nKL_XKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");

            if (IsItXsu == true) total_weight = total_weight * ReadWeightHist(Kstar_delta_weight, invM * invM);
            else total_weight = total_weight * ReadWeightHist(K0star_delta_weight, invM * invM);

            Nevt = Nevt + FillTemplate(hist, MVA_var, total_weight, Bsig_M);

        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double GetXnnBRPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, bool IsItUp, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get nominal PDF with uncertainties on all BR(B->X n nbar) decay
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
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
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double invM_Xnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;
    double N_Xplusnn = 0;
    double N_Xzeronn = 0;

    double Ngamma_v200 = -1;

    double s13_KpKLKL = -1;
    double s23_KpKLKL = -1;
    double nB2KpKLKL_all_KpKLKL = -1;
    double nB2KpKLKL_NR_KpKLKL = -1;

    double s13_KSKLKL = -1;
    double s23_KSKLKL = -1;
    double s12_KSKLKL = -1;
    double nB2KSKLKL_all_KSKLKL = -1;
    double nB2KSKLKL_NR_KSKLKL = -1;

    double nKL_XKLKL = -1;
    double XKLKL_E_1st = -1;
    double XKLKL_px_1st = -1;
    double XKLKL_py_1st = -1;
    double XKLKL_pz_1st = -1;
    double XKLKL_E_2nd = -1;
    double XKLKL_px_2nd = -1;
    double XKLKL_py_2nd = -1;
    double XKLKL_pz_2nd = -1;
    double nB2KstarKLKL = -1;
    double nB02KstarKLKL = -1;

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    double Bsig_M = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clXnn__bc", &N_Xplusnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clXnn__bc", &N_Xzeronn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clXnn__bc", &invM_Xnn);

        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s13_KpKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s23_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_all__bc", &nB2KpKLKL_all_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_NR__bc", &nB2KpKLKL_NR_KpKLKL);

        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s13_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &s23_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s12_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_all__bc", &nB2KSKLKL_all_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_NR__bc", &nB2KSKLKL_NR_KSKLKL);

        tree_upsilon->SetBranchAddress("nParticlesInList__boK_L0__clXKLKL__bc", &nKL_XKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__spE__bc", &XKLKL_E_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__sppx__bc", &XKLKL_px_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__sppy__bc", &XKLKL_py_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__sppz__bc", &XKLKL_pz_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__spE__bc", &XKLKL_E_2nd);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppx__bc", &XKLKL_px_2nd);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppy__bc", &XKLKL_py_2nd);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppz__bc", &XKLKL_pz_2nd);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarKLKL__bc", &nB2KstarKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstarKLKL__bc", &nB02KstarKLKL);

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            // select B+ --> K+ nu nubar reconstruction only
            if ((std::abs(Upsilon_ID) < MyEPSILON) && (std::abs(Bsig_ID) < MyEPSILON)) {}
            else continue;

            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE), temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE), temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(0, i_fake, MCTYPE), temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(1, i_fake, MCTYPE), temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(2, i_fake, MCTYPE), temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(3, i_fake, MCTYPE), temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(0, i_fake, MCTYPE), temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(1, i_fake, MCTYPE), temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(2, i_fake, MCTYPE), temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(3, i_fake, MCTYPE), temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);
            if ((N_Xplusnn + N_Xzeronn) > 0.5) { // uncertainty for all B->X n nbar, including K n nbar and Kstar n nbar
                if (IsItUp) Correction_Xnn = Correction_Xnn * 2.0;
                else Correction_Xnn = 0.0;
            }

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B --> K* KL KL correction factor
            double Correction_KstarKLKL = corrector_KstarKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KstarKLKL + nB02KstarKLKL);

            // B --> X KL KL correction factor
            double Correction_XKLKL = corrector_XsKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KpKLKL_all_KpKLKL, nB2KSKLKL_all_KSKLKL, nB2KstarKLKL + nB02KstarKLKL, nKL_XKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            Nevt = Nevt + FillTemplate(hist, MVA_var, total_weight, Bsig_M);

        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

int GetmultiplicitycorrelatedPDFs(const char* dirname, TH1D* CHG_nominal_hist, TH1D* MIX_nominal_hist, TH1D* UUBAR_nominal_hist, TH1D* DDBAR_nominal_hist, TH1D* SSBAR_nominal_hist, TH1D* CHARM_nominal_hist, TH1D* SIGNAL_nominal_hist, TH1D*** CHG_hists, TH1D*** MIX_hists, TH1D*** UUBAR_hists, TH1D*** DDBAR_hists, TH1D*** SSBAR_hists, TH1D*** CHARM_hists, TH1D*** SIGNAL_hists) { // get shape sys histogram from txt file
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

    *CHG_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *MIX_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *UUBAR_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *DDBAR_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *SSBAR_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *CHARM_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *SIGNAL_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);

    for (int i = 0; i < Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), ("CHG_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), ("MIX_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*UUBAR_hists)[i] = new TH1D(("UUBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), ("UUBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*DDBAR_hists)[i] = new TH1D(("DDBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), ("DDBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*SSBAR_hists)[i] = new TH1D(("SSBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), ("SSBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*CHARM_hists)[i] = new TH1D(("CHARM_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), ("CHARM_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*SIGNAL_hists)[i] = new TH1D(("Signal_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
    }

    for (int i = Nentry; i < 2 * Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHG_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("MIX_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*UUBAR_hists)[i] = new TH1D(("UUBAR_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("UUBAR_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*DDBAR_hists)[i] = new TH1D(("DDBAR_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("DDBAR_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*SSBAR_hists)[i] = new TH1D(("SSBAR_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("SSBAR_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*CHARM_hists)[i] = new TH1D(("CHARM_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHARM_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*SIGNAL_hists)[i] = new TH1D(("Signal_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
    }

    fp = fopen(dirname, "r");
    for (int i = 0; i < Nentry; i++) {
        fscanf(fp, "%lf\n", &eigen_value);
        for (int j = 0; j < RarityBins * 7; j++) fscanf(fp, "%lf\n", &weight_sys[j]);

        for (int k = 0; k < RarityBins; k++) {
            (*CHG_hists)[i]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 1 * RarityBins]));
            (*UUBAR_hists)[i]->SetBinContent(k + 1, UUBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 2 * RarityBins]));
            (*DDBAR_hists)[i]->SetBinContent(k + 1, DDBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 3 * RarityBins]));
            (*SSBAR_hists)[i]->SetBinContent(k + 1, SSBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 4 * RarityBins]));
            (*CHARM_hists)[i]->SetBinContent(k + 1, CHARM_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 5 * RarityBins]));
            (*SIGNAL_hists)[i]->SetBinContent(k + 1, SIGNAL_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 6 * RarityBins]));

            (*CHG_hists)[i + Nentry]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i + Nentry]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 1 * RarityBins]));
            (*UUBAR_hists)[i + Nentry]->SetBinContent(k + 1, UUBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 2 * RarityBins]));
            (*DDBAR_hists)[i + Nentry]->SetBinContent(k + 1, DDBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 3 * RarityBins]));
            (*SSBAR_hists)[i + Nentry]->SetBinContent(k + 1, SSBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 4 * RarityBins]));
            (*CHARM_hists)[i + Nentry]->SetBinContent(k + 1, CHARM_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 5 * RarityBins]));
            (*SIGNAL_hists)[i + Nentry]->SetBinContent(k + 1, SIGNAL_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 6 * RarityBins]));
        }

    }

    fclose(fp);

    return Nentry;
}

void GetmultiplicityUncorrelatedPDFs(const char* dirname, TH1D* CHG_hist, TH1D* MIX_hist, TH1D* UUBAR_hist, TH1D* DDBAR_hist, TH1D* SSBAR_hist, TH1D* CHARM_hist, TH1D* SIGNAL_hist) { // get shape sys histogram from txt file
    FILE* fp;
    fp = fopen(dirname, "r");

    double weight_sys[RarityBins * 7] = { 0.0 };
    for (int i = 0; i < RarityBins * 7; i++) fscanf(fp, "%lf\n", &weight_sys[i]);
    fclose(fp);

    for (int i = 0; i < RarityBins * 7; i++) weight_sys[i] = std::sqrt(weight_sys[i]);

    for (int i = 0; i < RarityBins; i++) CHG_hist->SetBinContent(i + 1, weight_sys[i]);
    for (int i = 0; i < RarityBins; i++) MIX_hist->SetBinContent(i + 1, weight_sys[RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) UUBAR_hist->SetBinContent(i + 1, weight_sys[2 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) DDBAR_hist->SetBinContent(i + 1, weight_sys[3 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) SSBAR_hist->SetBinContent(i + 1, weight_sys[4 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) CHARM_hist->SetBinContent(i + 1, weight_sys[5 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) SIGNAL_hist->SetBinContent(i + 1, weight_sys[6 * RarityBins + i]);
}

double GetBtoDtoXKLPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, bool IsItUp, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get nominal PDF with uncertainties on B-> [D -> X KL0] anything decays
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
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
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double invM_Xnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;
    double N_Xplusnn = 0;
    double N_Xzeronn = 0;

    double Ngamma_v200 = -1;

    double s13_KpKLKL = -1;
    double s23_KpKLKL = -1;
    double nB2KpKLKL_all_KpKLKL = -1;
    double nB2KpKLKL_NR_KpKLKL = -1;

    double s13_KSKLKL = -1;
    double s23_KSKLKL = -1;
    double s12_KSKLKL = -1;
    double nB2KSKLKL_all_KSKLKL = -1;
    double nB2KSKLKL_NR_KSKLKL = -1;

    double nKL_XKLKL = -1;
    double XKLKL_E_1st = -1;
    double XKLKL_px_1st = -1;
    double XKLKL_py_1st = -1;
    double XKLKL_pz_1st = -1;
    double XKLKL_E_2nd = -1;
    double XKLKL_px_2nd = -1;
    double XKLKL_py_2nd = -1;
    double XKLKL_pz_2nd = -1;
    double nB2KstarKLKL = -1;
    double nB02KstarKLKL = -1;

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    double Bsig_M = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);

        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s13_KpKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s23_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_all__bc", &nB2KpKLKL_all_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_NR__bc", &nB2KpKLKL_NR_KpKLKL);

        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s13_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &s23_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s12_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_all__bc", &nB2KSKLKL_all_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_NR__bc", &nB2KSKLKL_NR_KSKLKL);

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            // select B+ --> K+ nu nubar reconstruction only
            if ((std::abs(Upsilon_ID) < MyEPSILON) && (std::abs(Bsig_ID) < MyEPSILON)) {}
            else continue;

            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE), temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE), temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(0, i_fake, MCTYPE), temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(1, i_fake, MCTYPE), temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(2, i_fake, MCTYPE), temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(3, i_fake, MCTYPE), temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(0, i_fake, MCTYPE), temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(1, i_fake, MCTYPE), temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(2, i_fake, MCTYPE), temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(3, i_fake, MCTYPE), temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B --> K* KL KL correction factor
            double Correction_KstarKLKL = corrector_KstarKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KstarKLKL + nB02KstarKLKL);

            // B --> X KL KL correction factor
            double Correction_XKLKL = corrector_XsKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KpKLKL_all_KpKLKL, nB2KSKLKL_all_KSKLKL, nB2KstarKLKL + nB02KstarKLKL, nKL_XKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) {
                Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);
                if (IsItUp) Correction_BtoDtoXKL = Correction_BtoDtoXKL * (1 + corrector_BtoDtoXKL.GetRelativeUncertainty(nDptoXKL + nD0toXKL));
                else Correction_BtoDtoXKL = Correction_BtoDtoXKL * (1 - corrector_BtoDtoXKL.GetRelativeUncertainty(nDptoXKL + nD0toXKL));
            }

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            Nevt = Nevt + FillTemplate(hist, MVA_var, total_weight, Bsig_M);

        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double GetBRXKLKLPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, bool IsItUp, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get nominal PDF with uncertainties on all BR(B->X KL KL) decay
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
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
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = -1.0;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double invM_Xnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;
    double N_Xplusnn = 0;
    double N_Xzeronn = 0;

    double Ngamma_v200 = -1;

    double s13_KpKLKL = -1;
    double s23_KpKLKL = -1;
    double nB2KpKLKL_all_KpKLKL = -1;
    double nB2KpKLKL_NR_KpKLKL = -1;

    double s13_KSKLKL = -1;
    double s23_KSKLKL = -1;
    double s12_KSKLKL = -1;
    double nB2KSKLKL_all_KSKLKL = -1;
    double nB2KSKLKL_NR_KSKLKL = -1;

    double nKL_XKLKL = -1;
    double XKLKL_E_1st = -1;
    double XKLKL_px_1st = -1;
    double XKLKL_py_1st = -1;
    double XKLKL_pz_1st = -1;
    double XKLKL_E_2nd = -1;
    double XKLKL_px_2nd = -1;
    double XKLKL_py_2nd = -1;
    double XKLKL_pz_2nd = -1;
    double nB2KstarKLKL = -1;
    double nB02KstarKLKL = -1;

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    double Bsig_M = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        TTree* tree_Xs;
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);
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
        if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) {
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

            tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
        }
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clXnn__bc", &N_Xplusnn);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clXnn__bc", &N_Xzeronn);
        tree_upsilon->SetBranchAddress("invMassInLists__bon0__clXnn__bc", &invM_Xnn);

        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s13_KpKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s23_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_all__bc", &nB2KpKLKL_all_KpKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_NR__bc", &nB2KpKLKL_NR_KpKLKL);

        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &s13_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &s23_KSKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &s12_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_all__bc", &nB2KSKLKL_all_KSKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_NR__bc", &nB2KSKLKL_NR_KSKLKL);

        tree_upsilon->SetBranchAddress("nParticlesInList__boK_L0__clXKLKL__bc", &nKL_XKLKL);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__spE__bc", &XKLKL_E_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__sppx__bc", &XKLKL_px_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__sppy__bc", &XKLKL_py_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_1st__cm__sppz__bc", &XKLKL_pz_1st);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__spE__bc", &XKLKL_E_2nd);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppx__bc", &XKLKL_px_2nd);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppy__bc", &XKLKL_py_2nd);
        tree_upsilon->SetBranchAddress("averageValueInList__boK_L0__clXKLKL_2nd__cm__sppz__bc", &XKLKL_pz_2nd);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarKLKL__bc", &nB2KstarKLKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstarKLKL__bc", &nB02KstarKLKL);

        tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &nDptoXKL);
        tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &nD0toXKL);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if ((CorrectionType == "B2Knunu") || (CorrectionType == "B02K0nunu") || (CorrectionType == "B2Xsnunu") || (CorrectionType == "B02Xsnunu")) tree_Xs->GetEntry(j);

            // select B+ --> K+ nu nubar reconstruction only
            if ((std::abs(Upsilon_ID) < MyEPSILON) && (std::abs(Bsig_ID) < MyEPSILON)) {}
            else continue;

            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE), temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE), temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(0, i_fake, MCTYPE), temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(1, i_fake, MCTYPE), temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(2, i_fake, MCTYPE), temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(3, i_fake, MCTYPE), temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(0, i_fake, MCTYPE), temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(1, i_fake, MCTYPE), temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(2, i_fake, MCTYPE), temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(3, i_fake, MCTYPE), temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B --> K* KL KL correction factor
            double Correction_KstarKLKL = corrector_KstarKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KstarKLKL + nB02KstarKLKL);

            // B --> X KL KL correction factor
            double Correction_XKLKL = corrector_XsKLKL.GetCorrectionFactorAtGeneric(XKLKL_E_1st, XKLKL_px_1st, XKLKL_py_1st, XKLKL_pz_1st, XKLKL_E_2nd, XKLKL_px_2nd, XKLKL_py_2nd, XKLKL_pz_2nd, nB2KpKLKL_all_KpKLKL, nB2KSKLKL_all_KSKLKL, nB2KstarKLKL + nB02KstarKLKL, nKL_XKLKL);
            if (nKL_XKLKL > 1.5) { // uncertainty for all B->X KL KL, including K KL KL and Kstar KL KL
                if (IsItUp) Correction_XKLKL = Correction_XKLKL * 2.0;
                else Correction_XKLKL = 0.0;
            }

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            Nevt = Nevt + FillTemplate(hist, MVA_var, total_weight, Bsig_M);

        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double ReadWeightHist(TH1D* hist, double value) {
    int Bin = hist->FindBin(value);
    if (Bin < 1) Bin = 1;
    else if (Bin > hist->GetNbinsX()) Bin = hist->GetNbinsX();
    return hist->GetBinContent(Bin);
}

void ReadSignalModelingFile() {

    FILE* fp;
    int Nentry = 0;
    double RangeMIN = -1;
    double RangeMAX = -1;

    // Xsu_Hmb
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/signal_modeling/Xsu_Hmb_weight.txt","r");
    fscanf(fp, "%d %lf %lf\n", &Nentry, &RangeMIN, &RangeMAX);
    Xsu_Hmb_weight = new TH1D("Xsu_Hmb_weight", ";;", Nentry, RangeMIN, RangeMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsu_Hmb_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Xsu_Lmb_weight
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/signal_modeling/Xsu_Lmb_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &Nentry, &RangeMIN, &RangeMAX);
    Xsu_Lmb_weight = new TH1D("Xsu_Lmb_weight", ";;", Nentry, RangeMIN, RangeMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsu_Lmb_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Xsu_Hpf_weight
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/signal_modeling/Xsu_Hpf_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &Nentry, &RangeMIN, &RangeMAX);
    Xsu_Hpf_weight = new TH1D("Xsu_Hpf_weight", ";;", Nentry, RangeMIN, RangeMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsu_Hpf_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Xsu_Lpf_weight
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/signal_modeling/Xsu_Lpf_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &Nentry, &RangeMIN, &RangeMAX);
    Xsu_Lpf_weight = new TH1D("Xsu_Lpf_weight", ";;", Nentry, RangeMIN, RangeMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsu_Lpf_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Xsu_Htransition_weight
    fp = fopen("Xsu_Htransition_weight.txt", "r");
    fscanf(fp, "%d\n", &Nentry);
    if (Nentry != RarityBins) {
        printf("[ReadSignalModelingFile] Nentry and RarityBins are not the same!\n");
        exit(1);
    }
    Xsu_Htransition_weight = new TH1D("Xsu_Htransition_weight", ";;", Nentry, BinMIN, BinMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsu_Htransition_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Xsu_Ltransition_weight
    fp = fopen("Xsu_Ltransition_weight.txt", "r");
    fscanf(fp, "%d\n", &Nentry);
    if (Nentry != RarityBins) {
        printf("[ReadSignalModelingFile] Nentry and RarityBins are not the same!\n");
        exit(1);
    }
    Xsu_Ltransition_weight = new TH1D("Xsu_Ltransition_weight", ";;", Nentry, BinMIN, BinMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsu_Ltransition_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Xsd_Hmb
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/signal_modeling/Xsd_Hmb_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &Nentry, &RangeMIN, &RangeMAX);
    Xsd_Hmb_weight = new TH1D("Xsd_Hmb_weight", ";;", Nentry, RangeMIN, RangeMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsd_Hmb_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Xsd_Lmb_weight
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/signal_modeling/Xsd_Lmb_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &Nentry, &RangeMIN, &RangeMAX);
    Xsd_Lmb_weight = new TH1D("Xsd_Lmb_weight", ";;", Nentry, RangeMIN, RangeMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsd_Lmb_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Xsd_Hpf_weight
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/signal_modeling/Xsd_Hpf_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &Nentry, &RangeMIN, &RangeMAX);
    Xsd_Hpf_weight = new TH1D("Xsd_Hpf_weight", ";;", Nentry, RangeMIN, RangeMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsd_Hpf_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Xsd_Lpf_weight
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/signal_modeling/Xsd_Lpf_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &Nentry, &RangeMIN, &RangeMAX);
    Xsd_Lpf_weight = new TH1D("Xsd_Lpf_weight", ";;", Nentry, RangeMIN, RangeMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsd_Lpf_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Xsd_Htransition_weight
    fp = fopen("Xsd_Htransition_weight.txt", "r");
    fscanf(fp, "%d\n", &Nentry);
    if (Nentry != RarityBins) {
        printf("[ReadSignalModelingFile] Nentry and RarityBins are not the same!\n");
        exit(1);
    }
    Xsd_Htransition_weight = new TH1D("Xsd_Htransition_weight", ";;", Nentry, BinMIN, BinMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsd_Htransition_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Xsd_Ltransition_weight
    fp = fopen("Xsd_Ltransition_weight.txt", "r");
    fscanf(fp, "%d\n", &Nentry);
    if (Nentry != RarityBins) {
        printf("[ReadSignalModelingFile] Nentry and RarityBins are not the same!\n");
        exit(1);
    }
    Xsd_Ltransition_weight = new TH1D("Xsd_Ltransition_weight", ";;", Nentry, BinMIN, BinMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Xsd_Ltransition_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // Kstar_delta_weight
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/NEW_mKstar/Kstar_mKstar_uncertainty.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &Nentry, &RangeMIN, &RangeMAX);
    Kstar_delta_weight = new TH1D("Kstar_delta_weight", ";;", Nentry, RangeMIN, RangeMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Kstar_delta_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // K0star_delta_weight
    fp = fopen("/home/belle2/junewoo/storage_b1/bsub/systematic/NEW_mKstar/K0star_mKstar_uncertainty.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &Nentry, &RangeMIN, &RangeMAX);
    K0star_delta_weight = new TH1D("K0star_delta_weight", ";;", Nentry, RangeMIN, RangeMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        K0star_delta_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

}

void ClearHist(TH1D* Signal_hist, TH1D* CHG_hist, TH1D* MIX_hist, TH1D* UUBAR_hist, TH1D* DDBAR_hist, TH1D* SSBAR_hist, TH1D* CHARM_hist) {
    Signal_hist->Reset();
    CHG_hist->Reset();
    MIX_hist->Reset();
    UUBAR_hist->Reset();
    DDBAR_hist->Reset();
    SSBAR_hist->Reset();
    CHARM_hist->Reset();
}

void AddSQRTHist(TH1D* output_hist, TH1D* input_hist, int Nbin) {
    for (int i = 0; i < Nbin; i++) {
        double temp_1 = output_hist->GetBinContent(i + 1);
        double temp_2 = input_hist->GetBinContent(i + 1);
        double temp_3 = std::sqrt(temp_1 * temp_1 + temp_2 * temp_2);
        output_hist->SetBinContent(i + 1, temp_3);
    }
}

void AddSQRTHist(TH1D* output_hist, TH1D* input_hist1, TH1D* input_hist2, int Nbin) {
    for (int i = 0; i < Nbin; i++) {
        double input_1 = input_hist1->GetBinContent(i + 1);
        double input_2 = input_hist2->GetBinContent(i + 1);
        double output = std::sqrt(input_1 * input_1 + input_2 * input_2);
        output_hist->SetBinContent(i + 1, output);
    }
}

void GetMCstatisticalRelativeError(TH1D* nominal_hist, TH1D * MCstat_hist, int Nbin) {
    for (int i = 0; i < Nbin; i++) {
        double value = nominal_hist->GetBinContent(i + 1);
        double error = nominal_hist->GetBinError(i + 1);
        if(value < MyEPSILON) MCstat_hist->SetBinContent(i + 1, 0.0);
        else {
            double relative_error = error / value;
            MCstat_hist->SetBinContent(i + 1, relative_error);
        }
    }
}

void SaveSpecificMXsBin(TH1D* hist, int MXsBin) {
    if (MXsBin == 1) {
        for (int i = RarityBins_MX1; i < RarityBins; i++) {
            hist->SetBinError(i + 1, 0);
            hist->SetBinContent(i + 1, 0);
        }
    }
    else if (MXsBin == 2) {
        for (int i = 0; i < RarityBins_MX1; i++) {
            hist->SetBinError(i + 1, 0);
            hist->SetBinContent(i + 1, 0);
        }
        for (int i = RarityBins_MX1 + RarityBins_MX2; i < RarityBins; i++) {
            hist->SetBinError(i + 1, 0);
            hist->SetBinContent(i + 1, 0);
        }
    }
    else if (MXsBin == 3) {
        for (int i = 0; i < RarityBins_MX1 + RarityBins_MX2; i++) {
            hist->SetBinError(i + 1, 0);
            hist->SetBinContent(i + 1, 0);
        }
    }
    else if (MXsBin == 0) {}
    else {
        printf("Improper Bin index.\n");
        exit(1);
    }
}

int main()
{

    ReadSignalModelingFile();



    /* ====================================== */
    // Define PDFs for HistFactory
    // nominal PDFs
    TH1D* Signal_nominal = new TH1D("Signal_nominal", "Signal_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_nominal = new TH1D("CHG_nominal", "CHG_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_nominal = new TH1D("MIX_nominal", "MIX_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_nominal = new TH1D("UUBAR_nominal", "UUBAR_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_nominal = new TH1D("DDBAR_nominal", "DDBAR_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_nominal = new TH1D("SSBAR_nominal", "SSBAR_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_nominal = new TH1D("CHARM_nominal", "CHARM_nominal", RarityBins, BinMIN, BinMAX);

    // track uncertainty
    TH1D* Signal_track_p = new TH1D("Signal_track_p", "Signal_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_track_p = new TH1D("CHG_track_p", "CHG_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_track_p = new TH1D("MIX_track_p", "MIX_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_track_p = new TH1D("UUBAR_track_p", "UUBAR_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_track_p = new TH1D("DDBAR_track_p", "DDBAR_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_track_p = new TH1D("SSBAR_track_p", "SSBAR_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_track_p = new TH1D("CHARM_track_p", "CHARM_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_track_m = new TH1D("Signal_track_m", "Signal_track_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_track_m = new TH1D("CHG_track_m", "CHG_track_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_track_m = new TH1D("MIX_track_m", "MIX_track_m", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_track_m = new TH1D("UUBAR_track_m", "UUBAR_track_m", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_track_m = new TH1D("DDBAR_track_m", "DDBAR_track_m", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_track_m = new TH1D("SSBAR_track_m", "SSBAR_track_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_track_m = new TH1D("CHARM_track_m", "CHARM_track_m", RarityBins, BinMIN, BinMAX);

    // KS0 uncertainty
    TH1D* Signal_KS0_p = new TH1D("Signal_KS0_p", "Signal_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_KS0_p = new TH1D("CHG_KS0_p", "CHG_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_KS0_p = new TH1D("MIX_KS0_p", "MIX_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_KS0_p = new TH1D("UUBAR_KS0_p", "UUBAR_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_KS0_p = new TH1D("DDBAR_KS0_p", "DDBAR_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_KS0_p = new TH1D("SSBAR_KS0_p", "SSBAR_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_KS0_p = new TH1D("CHARM_KS0_p", "CHARM_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_KS0_m = new TH1D("Signal_KS0_m", "Signal_KS0_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_KS0_m = new TH1D("CHG_KS0_m", "CHG_KS0_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_KS0_m = new TH1D("MIX_KS0_m", "MIX_KS0_m", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_KS0_m = new TH1D("UUBAR_KS0_m", "UUBAR_KS0_m", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_KS0_m = new TH1D("DDBAR_KS0_m", "DDBAR_KS0_m", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_KS0_m = new TH1D("SSBAR_KS0_m", "SSBAR_KS0_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_KS0_m = new TH1D("CHARM_KS0_m", "CHARM_KS0_m", RarityBins, BinMIN, BinMAX);

    // FEI uncertainty (correlated)
    TH1D** Signal_FEI_correlated;
    TH1D** CHG_FEI_correlated;
    TH1D** MIX_FEI_correlated;

    // FEI uncertainty (uncorrelated)
    TH1D* Signal_FEI_uncorrelated = new TH1D("Signal_FEI_uncorrelated", "Signal_FEI_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_FEI_uncorrelated = new TH1D("CHG_FEI_uncorrelated", "CHG_FEI_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_FEI_uncorrelated = new TH1D("MIX_FEI_uncorrelated", "MIX_FEI_uncorrelated", RarityBins, BinMIN, BinMAX);

    // Kaon PID uncertainty (correlated)
    TH1D** Signal_KID_correlated;
    TH1D** CHG_KID_correlated;
    TH1D** MIX_KID_correlated;
    TH1D** UUBAR_KID_correlated;
    TH1D** DDBAR_KID_correlated;
    TH1D** SSBAR_KID_correlated;
    TH1D** CHARM_KID_correlated;

    // Kaon PID uncertainty (uncorrelated)
    TH1D* Signal_KID_uncorrelated = new TH1D("Signal_KID_uncorrelated", "Signal_KID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_KID_uncorrelated = new TH1D("CHG_KID_uncorrelated", "CHG_KID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_KID_uncorrelated = new TH1D("MIX_KID_uncorrelated", "MIX_KID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_KID_uncorrelated = new TH1D("UUBAR_KID_uncorrelated", "UUBAR_KID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_KID_uncorrelated = new TH1D("DDBAR_KID_uncorrelated", "DDBAR_KID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_KID_uncorrelated = new TH1D("SSBAR_KID_uncorrelated", "SSBAR_KID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_KID_uncorrelated = new TH1D("CHARM_KID_uncorrelated", "CHARM_KID_uncorrelated", RarityBins, BinMIN, BinMAX);

    // Pion PID uncertainty (correlated)
    TH1D** Signal_PID_correlated;
    TH1D** CHG_PID_correlated;
    TH1D** MIX_PID_correlated;
    TH1D** UUBAR_PID_correlated;
    TH1D** DDBAR_PID_correlated;
    TH1D** SSBAR_PID_correlated;
    TH1D** CHARM_PID_correlated;

    // Pion PID uncertainty (uncorrelated)
    TH1D* Signal_PID_uncorrelated = new TH1D("Signal_PID_uncorrelated", "Signal_PID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_PID_uncorrelated = new TH1D("CHG_PID_uncorrelated", "CHG_PID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_PID_uncorrelated = new TH1D("MIX_PID_uncorrelated", "MIX_PID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_PID_uncorrelated = new TH1D("UUBAR_PID_uncorrelated", "UUBAR_PID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_PID_uncorrelated = new TH1D("DDBAR_PID_uncorrelated", "DDBAR_PID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_PID_uncorrelated = new TH1D("SSBAR_PID_uncorrelated", "SSBAR_PID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_PID_uncorrelated = new TH1D("CHARM_PID_uncorrelated", "CHARM_PID_uncorrelated", RarityBins, BinMIN, BinMAX);

    // BB BR uncertainty(correlated)
    TH1D** Signal_BR_correlated;
    TH1D** CHG_BR_correlated;
    TH1D** MIX_BR_correlated;

    // BB BR uncertainty (uncorrelated)
    TH1D* Signal_BR_uncorrelated = new TH1D("Signal_BR_uncorrelated", "Signal_BR_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_BR_uncorrelated = new TH1D("CHG_BR_uncorrelated", "CHG_BR_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_BR_uncorrelated = new TH1D("MIX_BR_uncorrelated", "MIX_BR_uncorrelated", RarityBins, BinMIN, BinMAX);

    // pi0 uncertainty (correlated)
    TH1D** Signal_pi0_correlated;
    TH1D** CHG_pi0_correlated;
    TH1D** MIX_pi0_correlated;
    TH1D** UUBAR_pi0_correlated;
    TH1D** DDBAR_pi0_correlated;
    TH1D** SSBAR_pi0_correlated;
    TH1D** CHARM_pi0_correlated;

    // pi0 uncertainty (uncorrelated)
    TH1D* Signal_pi0_uncorrelated = new TH1D("Signal_pi0_uncorrelated", "Signal_pi0_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_pi0_uncorrelated = new TH1D("CHG_pi0_uncorrelated", "CHG_pi0_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_pi0_uncorrelated = new TH1D("MIX_pi0_uncorrelated", "MIX_pi0_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_pi0_uncorrelated = new TH1D("UUBAR_pi0_uncorrelated", "UUBAR_pi0_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_pi0_uncorrelated = new TH1D("DDBAR_pi0_uncorrelated", "DDBAR_pi0_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_pi0_uncorrelated = new TH1D("SSBAR_pi0_uncorrelated", "SSBAR_pi0_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_pi0_uncorrelated = new TH1D("CHARM_pi0_uncorrelated", "CHARM_pi0_uncorrelated", RarityBins, BinMIN, BinMAX);

    // K nu nubar form factor
    TH1D* Signal_Kff1_p = new TH1D("Signal_Kff1_p", "Signal_Kff1_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kff1_m = new TH1D("Signal_Kff1_m", "Signal_Kff1_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kff2_p = new TH1D("Signal_Kff2_p", "Signal_Kff2_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kff2_m = new TH1D("Signal_Kff2_m", "Signal_Kff2_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kff3_p = new TH1D("Signal_Kff3_p", "Signal_Kff3_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kff3_m = new TH1D("Signal_Kff3_m", "Signal_Kff3_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kff1_p = new TH1D("CHG_Kff1_p", "CHG_Kff1_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kff1_m = new TH1D("CHG_Kff1_m", "CHG_Kff1_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kff2_p = new TH1D("CHG_Kff2_p", "CHG_Kff2_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kff2_m = new TH1D("CHG_Kff2_m", "CHG_Kff2_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kff3_p = new TH1D("CHG_Kff3_p", "CHG_Kff3_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kff3_m = new TH1D("CHG_Kff3_m", "CHG_Kff3_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kff1_p = new TH1D("MIX_Kff1_p", "MIX_Kff1_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kff1_m = new TH1D("MIX_Kff1_m", "MIX_Kff1_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kff2_p = new TH1D("MIX_Kff2_p", "MIX_Kff2_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kff2_m = new TH1D("MIX_Kff2_m", "MIX_Kff2_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kff3_p = new TH1D("MIX_Kff3_p", "MIX_Kff3_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kff3_m = new TH1D("MIX_Kff3_m", "MIX_Kff3_m", RarityBins, BinMIN, BinMAX);

    // Kstar nu nubar form factor
    TH1D* CHG_Kstarff1_p = new TH1D("CHG_Kstarff1_p", "CHG_Kstarff1_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kstarff1_m = new TH1D("CHG_Kstarff1_m", "CHG_Kstarff1_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kstarff2_p = new TH1D("CHG_Kstarff2_p", "CHG_Kstarff2_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kstarff2_m = new TH1D("CHG_Kstarff2_m", "CHG_Kstarff2_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kstarff3_p = new TH1D("CHG_Kstarff3_p", "CHG_Kstarff3_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kstarff3_m = new TH1D("CHG_Kstarff3_m", "CHG_Kstarff3_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kstarff4_p = new TH1D("CHG_Kstarff4_p", "CHG_Kstarff4_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kstarff4_m = new TH1D("CHG_Kstarff4_m", "CHG_Kstarff4_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kstarff5_p = new TH1D("CHG_Kstarff5_p", "CHG_Kstarff5_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kstarff5_m = new TH1D("CHG_Kstarff5_m", "CHG_Kstarff5_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kstarff6_p = new TH1D("CHG_Kstarff6_p", "CHG_Kstarff6_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kstarff6_m = new TH1D("CHG_Kstarff6_m", "CHG_Kstarff6_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kstarff7_p = new TH1D("CHG_Kstarff7_p", "CHG_Kstarff7_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kstarff7_m = new TH1D("CHG_Kstarff7_m", "CHG_Kstarff7_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kstarff8_p = new TH1D("CHG_Kstarff8_p", "CHG_Kstarff8_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kstarff8_m = new TH1D("CHG_Kstarff8_m", "CHG_Kstarff8_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kstarff9_p = new TH1D("CHG_Kstarff9_p", "CHG_Kstarff9_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kstarff9_m = new TH1D("CHG_Kstarff9_m", "CHG_Kstarff9_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kstarff1_p = new TH1D("MIX_Kstarff1_p", "MIX_Kstarff1_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kstarff1_m = new TH1D("MIX_Kstarff1_m", "MIX_Kstarff1_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kstarff2_p = new TH1D("MIX_Kstarff2_p", "MIX_Kstarff2_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kstarff2_m = new TH1D("MIX_Kstarff2_m", "MIX_Kstarff2_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kstarff3_p = new TH1D("MIX_Kstarff3_p", "MIX_Kstarff3_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kstarff3_m = new TH1D("MIX_Kstarff3_m", "MIX_Kstarff3_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kstarff4_p = new TH1D("MIX_Kstarff4_p", "MIX_Kstarff4_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kstarff4_m = new TH1D("MIX_Kstarff4_m", "MIX_Kstarff4_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kstarff5_p = new TH1D("MIX_Kstarff5_p", "MIX_Kstarff5_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kstarff5_m = new TH1D("MIX_Kstarff5_m", "MIX_Kstarff5_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kstarff6_p = new TH1D("MIX_Kstarff6_p", "MIX_Kstarff6_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kstarff6_m = new TH1D("MIX_Kstarff6_m", "MIX_Kstarff6_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kstarff7_p = new TH1D("MIX_Kstarff7_p", "MIX_Kstarff7_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kstarff7_m = new TH1D("MIX_Kstarff7_m", "MIX_Kstarff7_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kstarff8_p = new TH1D("MIX_Kstarff8_p", "MIX_Kstarff8_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kstarff8_m = new TH1D("MIX_Kstarff8_m", "MIX_Kstarff8_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kstarff9_p = new TH1D("MIX_Kstarff9_p", "MIX_Kstarff9_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kstarff9_m = new TH1D("MIX_Kstarff9_m", "MIX_Kstarff9_m", RarityBins, BinMIN, BinMAX);

    // K nu nubar OLD form factor
    TH1D* Signal_Kff_OLD_p = new TH1D("Signal_Kff_OLD_p", "Signal_Kff_OLD_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kff_OLD_p = new TH1D("CHG_Kff_OLD_p", "CHG_Kff_OLD_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kff_OLD_p = new TH1D("MIX_Kff_OLD_p", "MIX_Kff_OLD_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kff_OLD_m = new TH1D("Signal_Kff_OLD_m", "Signal_Kff_OLD_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kff_OLD_m = new TH1D("CHG_Kff_OLD_m", "CHG_Kff_OLD_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kff_OLD_m = new TH1D("MIX_Kff_OLD_m", "MIX_Kff_OLD_m", RarityBins, BinMIN, BinMAX);

    // fraction uncertainty(correlated)
    TH1D** Signal_fraction_correlated;
    TH1D** CHG_fraction_correlated;
    TH1D** MIX_fraction_correlated;

    // fraction uncertainty (uncorrelated)
    TH1D* Signal_fraction_uncorrelated = new TH1D("Signal_fraction_uncorrelated", "Signal_fraction_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_fraction_uncorrelated = new TH1D("CHG_fraction_uncorrelated", "CHG_fraction_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_fraction_uncorrelated = new TH1D("MIX_fraction_uncorrelated", "MIX_fraction_uncorrelated", RarityBins, BinMIN, BinMAX);

    // fragmentation uncertainty(correlated)
    TH1D** Signal_Fragmentation_correlated;
    TH1D** CHG_Fragmentation_correlated;
    TH1D** MIX_Fragmentation_correlated;

    // fragmentation uncertainty (uncorrelated)
    TH1D* Signal_Fragmentation_uncorrelated = new TH1D("Signal_Fragmentation_uncorrelated", "Signal_Fragmentation_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Fragmentation_uncorrelated = new TH1D("CHG_Fragmentation_uncorrelated", "CHG_Fragmentation_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Fragmentation_uncorrelated = new TH1D("MIX_Fragmentation_uncorrelated", "MIX_Fragmentation_uncorrelated", RarityBins, BinMIN, BinMAX);

    // pf
    TH1D* CHG_pf_p = new TH1D("CHG_pf_p", "CHG_pf_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_pf_p = new TH1D("MIX_pf_p", "MIX_pf_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_pf_m = new TH1D("CHG_pf_m", "CHG_pf_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_pf_m = new TH1D("MIX_pf_m", "MIX_pf_m", RarityBins, BinMIN, BinMAX);

    // mb
    TH1D* CHG_mb_p = new TH1D("CHG_mb_p", "CHG_mb_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_mb_p = new TH1D("MIX_mb_p", "MIX_mb_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_mb_m = new TH1D("CHG_mb_m", "CHG_mb_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_mb_m = new TH1D("MIX_mb_m", "MIX_mb_m", RarityBins, BinMIN, BinMAX);

    // transition
    TH1D* CHG_transition_p = new TH1D("CHG_transition_p", "CHG_transition_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_transition_p = new TH1D("MIX_transition_p", "MIX_transition_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_transition_m = new TH1D("CHG_transition_m", "CHG_transition_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_transition_m = new TH1D("MIX_transition_m", "MIX_transition_m", RarityBins, BinMIN, BinMAX);

    // mK*
    TH1D* CHG_mKstar_p = new TH1D("CHG_mKstar_p", "CHG_mKstar_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_mKstar_p = new TH1D("MIX_mKstar_p", "MIX_mKstar_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_mKstar_m = new TH1D("CHG_mKstar_m", "CHG_mKstar_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_mKstar_m = new TH1D("MIX_mKstar_m", "MIX_mKstar_m", RarityBins, BinMIN, BinMAX);

    // BDTc
    //TH1D* CHG_BDTc_p = new TH1D("CHG_BDTc_p", "CHG_BDTc_p", RarityBins, BinMIN, BinMAX);
    //TH1D* MIX_BDTc_p = new TH1D("MIX_BDTc_p", "MIX_BDTc_p", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_BDTc_p = new TH1D("UUBAR_BDTc_p", "UUBAR_BDTc_p", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_BDTc_p = new TH1D("DDBAR_BDTc_p", "DDBAR_BDTc_p", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_BDTc_p = new TH1D("SSBAR_BDTc_p", "SSBAR_BDTc_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_BDTc_p = new TH1D("CHARM_BDTc_p", "CHARM_BDTc_p", RarityBins, BinMIN, BinMAX);
    //TH1D* CHG_BDTc_m = new TH1D("CHG_BDTc_m", "CHG_BDTc_m", RarityBins, BinMIN, BinMAX);
    //TH1D* MIX_BDTc_m = new TH1D("MIX_BDTc_m", "MIX_BDTc_m", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_BDTc_m = new TH1D("UUBAR_BDTc_m", "UUBAR_BDTc_m", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_BDTc_m = new TH1D("DDBAR_BDTc_m", "DDBAR_BDTc_m", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_BDTc_m = new TH1D("SSBAR_BDTc_m", "SSBAR_BDTc_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_BDTc_m = new TH1D("CHARM_BDTc_m", "CHARM_BDTc_m", RarityBins, BinMIN, BinMAX);

    // B2Knn BR
    TH1D* Signal_Xnn_p = new TH1D("Signal_Xnn_p", "Signal_Xnn_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Xnn_p = new TH1D("CHG_Xnn_p", "CHG_Xnn_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Xnn_p = new TH1D("MIX_Xnn_p", "MIX_Xnn_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Xnn_m = new TH1D("Signal_Xnn_m", "Signal_Xnn_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Xnn_m = new TH1D("CHG_Xnn_m", "CHG_Xnn_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Xnn_m = new TH1D("MIX_Xnn_m", "MIX_Xnn_m", RarityBins, BinMIN, BinMAX);

    // multiplicity uncertainty (correlated)
    TH1D** Signal_multiplicity_correlated;
    TH1D** CHG_multiplicity_correlated;
    TH1D** MIX_multiplicity_correlated;
    TH1D** UUBAR_multiplicity_correlated;
    TH1D** DDBAR_multiplicity_correlated;
    TH1D** SSBAR_multiplicity_correlated;
    TH1D** CHARM_multiplicity_correlated;

    // multiplicity uncertainty (uncorrelated)
    TH1D* Signal_multiplicity_uncorrelated = new TH1D("Signal_multiplicity_uncorrelated", "Signal_multiplicity_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_multiplicity_uncorrelated = new TH1D("CHG_multiplicity_uncorrelated", "CHG_multiplicity_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_multiplicity_uncorrelated = new TH1D("MIX_multiplicity_uncorrelated", "MIX_multiplicity_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_multiplicity_uncorrelated = new TH1D("UUBAR_multiplicity_uncorrelated", "UUBAR_multiplicity_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_multiplicity_uncorrelated = new TH1D("DDBAR_multiplicity_uncorrelated", "DDBAR_multiplicity_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_multiplicity_uncorrelated = new TH1D("SSBAR_multiplicity_uncorrelated", "SSBAR_multiplicity_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_multiplicity_uncorrelated = new TH1D("CHARM_multiplicity_uncorrelated", "CHARM_multiplicity_uncorrelated", RarityBins, BinMIN, BinMAX);

    // B -> [D -> X KL0] anything
    TH1D* Signal_BtoDtoXKL_p = new TH1D("Signal_BtoDtoXKL_p", "Signal_BtoDtoXKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_BtoDtoXKL_p = new TH1D("CHG_BtoDtoXKL_p", "CHG_BtoDtoXKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_BtoDtoXKL_p = new TH1D("MIX_BtoDtoXKL_p", "MIX_BtoDtoXKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_BtoDtoXKL_m = new TH1D("Signal_BtoDtoXKL_m", "Signal_BtoDtoXKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_BtoDtoXKL_m = new TH1D("CHG_BtoDtoXKL_m", "CHG_BtoDtoXKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_BtoDtoXKL_m = new TH1D("MIX_BtoDtoXKL_m", "MIX_BtoDtoXKL_m", RarityBins, BinMIN, BinMAX);

    // BR(B -> X KL KL)
    TH1D* Signal_BRBtoXKLKL_p = new TH1D("Signal_BRBtoXKLKL_p", "Signal_BRBtoXKLKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_BRBtoXKLKL_p = new TH1D("CHG_BRBtoXKLKL_p", "CHG_BRBtoXKLKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_BRBtoXKLKL_p = new TH1D("MIX_BRBtoXKLKL_p", "MIX_BRBtoXKLKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_BRBtoXKLKL_m = new TH1D("Signal_BRBtoXKLKL_m", "Signal_BRBtoXKLKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_BRBtoXKLKL_m = new TH1D("CHG_BRBtoXKLKL_m", "CHG_BRBtoXKLKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_BRBtoXKLKL_m = new TH1D("MIX_BRBtoXKLKL_m", "MIX_BRBtoXKLKL_m", RarityBins, BinMIN, BinMAX);

    // all of uncorrelated uncertainties
    TH1D* Signal_all_uncorrelated = new TH1D("Signal_all_uncorrelated", "Signal_all_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_all_uncorrelated = new TH1D("CHG_all_uncorrelated", "CHG_all_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_all_uncorrelated = new TH1D("MIX_all_uncorrelated", "MIX_all_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_all_uncorrelated = new TH1D("UUBAR_all_uncorrelated", "UUBAR_all_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_all_uncorrelated = new TH1D("DDBAR_all_uncorrelated", "DDBAR_all_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_all_uncorrelated = new TH1D("SSBAR_all_uncorrelated", "SSBAR_all_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_all_uncorrelated = new TH1D("CHARM_all_uncorrelated", "CHARM_all_uncorrelated", RarityBins, BinMIN, BinMAX);

    // MC statistical uncertainties
    TH1D* Signal_MC_stat = new TH1D("Signal_MC_stat", "Signal_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_MC_stat = new TH1D("CHG_MC_stat", "CHG_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_MC_stat = new TH1D("MIX_MC_stat", "MIX_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_MC_stat = new TH1D("UUBAR_MC_stat", "UUBAR_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_MC_stat = new TH1D("DDBAR_MC_stat", "DDBAR_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_MC_stat = new TH1D("SSBAR_MC_stat", "SSBAR_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_MC_stat = new TH1D("CHARM_MC_stat", "CHARM_MC_stat", RarityBins, BinMIN, BinMAX);

    // all of uncorrelated uncertainties + MC statistical uncertainties
    TH1D* Signal_all_uncorrelated_MC_stat = new TH1D("Signal_all_uncorrelated_MC_stat", "Signal_all_uncorrelated_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_all_uncorrelated_MC_stat = new TH1D("CHG_all_uncorrelated_MC_stat", "CHG_all_uncorrelated_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_all_uncorrelated_MC_stat = new TH1D("MIX_all_uncorrelated_MC_stat", "MIX_all_uncorrelated_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_all_uncorrelated_MC_stat = new TH1D("UUBAR_all_uncorrelated_MC_stat", "UUBAR_all_uncorrelated_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_all_uncorrelated_MC_stat = new TH1D("DDBAR_all_uncorrelated_MC_stat", "DDBAR_all_uncorrelated_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_all_uncorrelated_MC_stat = new TH1D("SSBAR_all_uncorrelated_MC_stat", "SSBAR_all_uncorrelated_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_all_uncorrelated_MC_stat = new TH1D("CHARM_all_uncorrelated_MC_stat", "CHARM_all_uncorrelated_MC_stat", RarityBins, BinMIN, BinMAX);
    /* ====================================== */



    /* ====================================== */
    // define path for Ntuple
    const char* MC_dirname_SIGNAL = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/Merge";

    const char* MC_dirname_CHG = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori/CHG_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/Merge";
    const char* MC_dirname_MIX = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori/MIX_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/Merge";
    const char* MC_dirname_UUBAR = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori/UUBAR_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/Merge";
    const char* MC_dirname_DDBAR = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori/DDBAR_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/Merge";
    const char* MC_dirname_SSBAR = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori/SSBAR_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/Merge";
    const char* MC_dirname_CHARM = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori/CHARM_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/Merge";
 
    // for FEI
    const char* FEI_correlated_info = "./FEI_selected.txt";
    const char* FEI_uncorrelated_info = "./FEI_cov_remain_truncated.txt";

    // for PID
    const char* KID_correlated_info = "./KID_selected.txt";
    const char* PID_correlated_info = "./PID_selected.txt";
    const char* KID_uncorrelated_info = "./KID_cov_remain_truncated.txt";
    const char* PID_uncorrelated_info = "./PID_cov_remain_truncated.txt";

    // for BB BR
    const char* BR_correlated_info = "./BR_selected.txt";
    const char* BR_uncorrelated_info = "./BR_cov_remain_truncated.txt";

    // for pi0
    const char* pi0_correlated_info = "./pi0_selected.txt";
    const char* pi0_uncorrelated_info = "./pi0_cov_remain_truncated.txt";

    // for multiplicity
    const char* multiplicity_correlated_info = "./multiplicity_selected.txt";
    const char* multiplicity_uncorrelated_info = "./multiplicity_cov_remain_truncated.txt";

    // for fraction
    const char* fraction_correlated_info = "./fraction_selected.txt";
    const char* fraction_uncorrelated_info = "./fraction_cov_remain_truncated.txt";

    // for fragmentation
    const char* Fragmentation_correlated_info = "./Fragmentation_selected.txt";
    const char* Fragmentation_uncorrelated_info = "./Fragmentation_cov_remain_truncated.txt";

    /* ====================================== */



    /* ====================================== */
    // Get PDFs
    // get nominal pdfs
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_nominal, "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", CHG_nominal, "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Xsnunu", CHG_nominal, "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", MIX_nominal, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", MIX_nominal, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Xsnunu", MIX_nominal, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    double NevtCHGWithoutBDTc = GetNominalPDFs(MC_dirname_CHG, "root", CHG_nominal, "Bplus", "CHG", ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    double NevtMIXWithoutBDTc = GetNominalPDFs(MC_dirname_MIX, "root", MIX_nominal, "Bzero", "MIX", ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");
    double NevtUUBARWithoutBDTc = GetNominalPDFs(MC_dirname_UUBAR, "root", UUBAR_nominal, "Continuum", "UUBAR", ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"), "otherwise");
    double NevtDDBARWithoutBDTc = GetNominalPDFs(MC_dirname_DDBAR, "root", DDBAR_nominal, "Continuum", "DDBAR", ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"), "otherwise");
    double NevtSSBARWithoutBDTc = GetNominalPDFs(MC_dirname_SSBAR, "root", SSBAR_nominal, "Continuum", "SSBAR", ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"), "otherwise");
    double NevtCHARMWithoutBDTc = GetNominalPDFs(MC_dirname_CHARM, "root", CHARM_nominal, "Continuum", "CHARM", ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"), "otherwise");

    // get track uncertainty pdfs
    GetTrackPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_track_p, "Bplus", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    GetTrackPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", CHG_track_p, "Bplus", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    GetTrackPDFs(MC_dirname_SIGNAL, "B2Xsnunu", CHG_track_p, "Bplus", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    GetTrackPDFs(MC_dirname_SIGNAL, "B02K0nunu", MIX_track_p, "Bzero", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetTrackPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", MIX_track_p, "Bzero", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    GetTrackPDFs(MC_dirname_SIGNAL, "B02Xsnunu", MIX_track_p, "Bzero", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    GetTrackPDFs(MC_dirname_CHG, "root", CHG_track_p, "Bplus", "CHG", true, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    GetTrackPDFs(MC_dirname_MIX, "root", MIX_track_p, "Bzero", "MIX", true, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");
    GetTrackPDFs(MC_dirname_UUBAR, "root", UUBAR_track_p, "Continuum", "UUBAR", true, ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"), "otherwise");
    GetTrackPDFs(MC_dirname_DDBAR, "root", DDBAR_track_p, "Continuum", "DDBAR", true, ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"), "otherwise");
    GetTrackPDFs(MC_dirname_SSBAR, "root", SSBAR_track_p, "Continuum", "SSBAR", true, ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"), "otherwise");
    GetTrackPDFs(MC_dirname_CHARM, "root", CHARM_track_p, "Continuum", "CHARM", true, ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"), "otherwise");

    GetTrackPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_track_m, "Bplus", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    GetTrackPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", CHG_track_m, "Bplus", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    GetTrackPDFs(MC_dirname_SIGNAL, "B2Xsnunu", CHG_track_m, "Bplus", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    GetTrackPDFs(MC_dirname_SIGNAL, "B02K0nunu", MIX_track_m, "Bzero", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetTrackPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", MIX_track_m, "Bzero", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    GetTrackPDFs(MC_dirname_SIGNAL, "B02Xsnunu", MIX_track_m, "Bzero", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    GetTrackPDFs(MC_dirname_CHG, "root", CHG_track_m, "Bplus", "CHG", false, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    GetTrackPDFs(MC_dirname_MIX, "root", MIX_track_m, "Bzero", "MIX", false, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");
    GetTrackPDFs(MC_dirname_UUBAR, "root", UUBAR_track_m, "Continuum", "UUBAR", false, ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"), "otherwise");
    GetTrackPDFs(MC_dirname_DDBAR, "root", DDBAR_track_m, "Continuum", "DDBAR", false, ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"), "otherwise");
    GetTrackPDFs(MC_dirname_SSBAR, "root", SSBAR_track_m, "Continuum", "SSBAR", false, ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"), "otherwise");
    GetTrackPDFs(MC_dirname_CHARM, "root", CHARM_track_m, "Continuum", "CHARM", false, ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"), "otherwise");

    // get KS0 uncertainty pdfs
    GetKS0PDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_KS0_p, "Bplus", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    GetKS0PDFs(MC_dirname_SIGNAL, "B2Kstarnunu", CHG_KS0_p, "Bplus", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    GetKS0PDFs(MC_dirname_SIGNAL, "B2Xsnunu", CHG_KS0_p, "Bplus", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    GetKS0PDFs(MC_dirname_SIGNAL, "B02K0nunu", MIX_KS0_p, "Bzero", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetKS0PDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", MIX_KS0_p, "Bzero", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    GetKS0PDFs(MC_dirname_SIGNAL, "B02Xsnunu", MIX_KS0_p, "Bzero", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    GetKS0PDFs(MC_dirname_CHG, "root", CHG_KS0_p, "Bplus", "CHG", true, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    GetKS0PDFs(MC_dirname_MIX, "root", MIX_KS0_p, "Bzero", "MIX", true, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");
    GetKS0PDFs(MC_dirname_UUBAR, "root", UUBAR_KS0_p, "Continuum", "UUBAR", true, ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"), "otherwise");
    GetKS0PDFs(MC_dirname_DDBAR, "root", DDBAR_KS0_p, "Continuum", "DDBAR", true, ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"), "otherwise");
    GetKS0PDFs(MC_dirname_SSBAR, "root", SSBAR_KS0_p, "Continuum", "SSBAR", true, ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"), "otherwise");
    GetKS0PDFs(MC_dirname_CHARM, "root", CHARM_KS0_p, "Continuum", "CHARM", true, ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"), "otherwise");

    GetKS0PDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_KS0_m, "Bplus", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    GetKS0PDFs(MC_dirname_SIGNAL, "B2Kstarnunu", CHG_KS0_m, "Bplus", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    GetKS0PDFs(MC_dirname_SIGNAL, "B2Xsnunu", CHG_KS0_m, "Bplus", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    GetKS0PDFs(MC_dirname_SIGNAL, "B02K0nunu", MIX_KS0_m, "Bzero", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetKS0PDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", MIX_KS0_m, "Bzero", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    GetKS0PDFs(MC_dirname_SIGNAL, "B02Xsnunu", MIX_KS0_m, "Bzero", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    GetKS0PDFs(MC_dirname_CHG, "root", CHG_KS0_m, "Bplus", "CHG", false, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    GetKS0PDFs(MC_dirname_MIX, "root", MIX_KS0_m, "Bzero", "MIX", false, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");
    GetKS0PDFs(MC_dirname_UUBAR, "root", UUBAR_KS0_m, "Continuum", "UUBAR", false, ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"), "otherwise");
    GetKS0PDFs(MC_dirname_DDBAR, "root", DDBAR_KS0_m, "Continuum", "DDBAR", false, ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"), "otherwise");
    GetKS0PDFs(MC_dirname_SSBAR, "root", SSBAR_KS0_m, "Continuum", "SSBAR", false, ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"), "otherwise");
    GetKS0PDFs(MC_dirname_CHARM, "root", CHARM_KS0_m, "Continuum", "CHARM", false, ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"), "otherwise");

    // get FEI uncertainty pdfs (correlated)
    int NPDFs_FEI = GetFEIcorrelatedPDFs(FEI_correlated_info, CHG_nominal, MIX_nominal, Signal_nominal, &CHG_FEI_correlated, &MIX_FEI_correlated, &Signal_FEI_correlated);

    // get FEI uncertainty pdfs (uncorrelated)
    GetFEIUncorrelatedPDFs(FEI_uncorrelated_info, CHG_FEI_uncorrelated, MIX_FEI_uncorrelated, Signal_FEI_uncorrelated);

    // get KID uncertainty pdfs (correlated)
    int NPDFs_KID = GetPIDcorrelatedPDFs(KID_correlated_info, CHG_nominal, MIX_nominal, UUBAR_nominal, DDBAR_nominal, SSBAR_nominal, CHARM_nominal, Signal_nominal, &CHG_KID_correlated, &MIX_KID_correlated, &UUBAR_KID_correlated, &DDBAR_KID_correlated, &SSBAR_KID_correlated, &CHARM_KID_correlated, &Signal_KID_correlated, true);
    int NPDFs_PID = GetPIDcorrelatedPDFs(PID_correlated_info, CHG_nominal, MIX_nominal, UUBAR_nominal, DDBAR_nominal, SSBAR_nominal, CHARM_nominal, Signal_nominal, &CHG_PID_correlated, &MIX_PID_correlated, &UUBAR_PID_correlated, &DDBAR_PID_correlated, &SSBAR_PID_correlated, &CHARM_PID_correlated, &Signal_PID_correlated, false);

    // get KID uncertainty pdfs (uncorrelated)
    GetPIDUncorrelatedPDFs(KID_uncorrelated_info, CHG_KID_uncorrelated, MIX_KID_uncorrelated, UUBAR_KID_uncorrelated, DDBAR_KID_uncorrelated, SSBAR_KID_uncorrelated, CHARM_KID_uncorrelated, Signal_KID_uncorrelated);
    GetPIDUncorrelatedPDFs(PID_uncorrelated_info, CHG_PID_uncorrelated, MIX_PID_uncorrelated, UUBAR_PID_uncorrelated, DDBAR_PID_uncorrelated, SSBAR_PID_uncorrelated, CHARM_PID_uncorrelated, Signal_PID_uncorrelated);

    // get BB BR uncertainty pdfs (correlated)
    int NPDFs_BR = GetBRcorrelatedPDFs(BR_correlated_info, CHG_nominal, MIX_nominal, Signal_nominal, &CHG_BR_correlated, &MIX_BR_correlated, &Signal_BR_correlated);

    // get BB BR uncertainty pdfs (uncorrelated)
    GetBRUncorrelatedPDFs(BR_uncorrelated_info, CHG_BR_uncorrelated, MIX_BR_uncorrelated, Signal_BR_uncorrelated);

    // get pi0 uncertainty pdfs (correlated)
    int NPDFs_pi0 = Getpi0correlatedPDFs(pi0_correlated_info, CHG_nominal, MIX_nominal, UUBAR_nominal, DDBAR_nominal, SSBAR_nominal, CHARM_nominal, Signal_nominal, &CHG_pi0_correlated, &MIX_pi0_correlated, &UUBAR_pi0_correlated, &DDBAR_pi0_correlated, &SSBAR_pi0_correlated, &CHARM_pi0_correlated, &Signal_pi0_correlated);

    // get pi0 uncertainty pdfs (uncorrelated)
    Getpi0UncorrelatedPDFs(pi0_uncorrelated_info, CHG_pi0_uncorrelated, MIX_pi0_uncorrelated, UUBAR_pi0_uncorrelated, DDBAR_pi0_uncorrelated, SSBAR_pi0_uncorrelated, CHARM_pi0_uncorrelated, Signal_pi0_uncorrelated);

    // get Kff uncertainty pdfs
    TH1D* Signal_Kff_temp = new TH1D("Signal_Kff_temp", "Signal_Kff_temp", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kff_array[7] = { Signal_Kff_temp, Signal_Kff1_p, Signal_Kff1_m, Signal_Kff2_p, Signal_Kff2_m, Signal_Kff3_p, Signal_Kff3_m };
    TH1D* CHG_Kff_temp = new TH1D("CHG_Kff_temp", "CHG_Kff_temp", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kff_array[7] = { CHG_Kff_temp, CHG_Kff1_p, CHG_Kff1_m, CHG_Kff2_p, CHG_Kff2_m, CHG_Kff3_p, CHG_Kff3_m };
    TH1D* MIX_Kff_temp = new TH1D("MIX_Kff_temp", "MIX_Kff_temp", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kff_array[7] = { MIX_Kff_temp, MIX_Kff1_p, MIX_Kff1_m, MIX_Kff2_p, MIX_Kff2_m, MIX_Kff3_p, MIX_Kff3_m };
    double Correction_factor_BR_Kplus[7] = { 0.0 };
    double Correction_factor_BR_Kzero[7] = { 0.0 };

    GetKffPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_Kff_array, Correction_factor_BR_Kplus, "Bplus", 1, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    for (int i = 0; i < 7; i++) GetNominalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", CHG_Kff_array[i], "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    for (int i = 0; i < 7; i++) {
        const double non_res_Xsu_correction_factor = (BR_Xsu_nonresonant_nunubar + BR_Kplus_nunubar - BR_Kplus_nunubar * Correction_factor_BR_Kplus[i]) / BR_Xsu_nonresonant_nunubar;
        GetNominalPDFs(MC_dirname_SIGNAL, "B2Xsnunu", CHG_Kff_array[i], "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu") * non_res_Xsu_correction_factor, "B2Xsnunu");
    }
    GetKffPDFs(MC_dirname_SIGNAL, "B02K0nunu", MIX_Kff_array, Correction_factor_BR_Kzero, "Bzero", 0, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    for (int i = 0; i < 7; i++) GetNominalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", MIX_Kff_array[i], "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    for (int i = 0; i < 7; i++) {
        const double non_res_Xsd_correction_factor = (BR_Xsd_nonresonant_nunubar + BR_K0_nunubar - BR_K0_nunubar * Correction_factor_BR_Kzero[i]) / BR_Xsd_nonresonant_nunubar;
        GetNominalPDFs(MC_dirname_SIGNAL, "B02Xsnunu", MIX_Kff_array[i], "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu") * non_res_Xsd_correction_factor, "B02Xsnunu");
    }

    for (int i = 0; i < 7; i++) {
        GetNominalPDFs(MC_dirname_CHG, "root", CHG_Kff_array[i], "Bplus", "CHG", ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
        GetNominalPDFs(MC_dirname_MIX, "root", MIX_Kff_array[i], "Bzero", "MIX", ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");
    }

    // get Kstar ff uncertainty pdfs
    TH1D* CHG_Kstarff_temp = new TH1D("CHG_Kstarff_temp", "CHG_Kstarff_temp", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Kstarff_array[19] = { CHG_Kstarff_temp,
        CHG_Kstarff1_p, CHG_Kstarff1_m,
        CHG_Kstarff2_p, CHG_Kstarff2_m,
        CHG_Kstarff3_p, CHG_Kstarff3_m,
        CHG_Kstarff4_p, CHG_Kstarff4_m,
        CHG_Kstarff5_p, CHG_Kstarff5_m,
        CHG_Kstarff6_p, CHG_Kstarff6_m,
        CHG_Kstarff7_p, CHG_Kstarff7_m,
        CHG_Kstarff8_p, CHG_Kstarff8_m,
        CHG_Kstarff9_p, CHG_Kstarff9_m
    };
    TH1D* MIX_Kstarff_temp = new TH1D("MIX_Kstarff_temp", "MIX_Kstarff_temp", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Kstarff_array[19] = { MIX_Kstarff_temp,
        MIX_Kstarff1_p, MIX_Kstarff1_m,
        MIX_Kstarff2_p, MIX_Kstarff2_m,
        MIX_Kstarff3_p, MIX_Kstarff3_m,
        MIX_Kstarff4_p, MIX_Kstarff4_m,
        MIX_Kstarff5_p, MIX_Kstarff5_m,
        MIX_Kstarff6_p, MIX_Kstarff6_m,
        MIX_Kstarff7_p, MIX_Kstarff7_m,
        MIX_Kstarff8_p, MIX_Kstarff8_m,
        MIX_Kstarff9_p, MIX_Kstarff9_m
    };
    double Correction_factor_BR_Kstarplus[19] = { 0.0 };
    double Correction_factor_BR_Kstarzero[19] = { 0.0 };

    GetKstarffPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", CHG_Kstarff_array, Correction_factor_BR_Kstarplus, "Bplus", 1, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"));
    for (int i = 0; i < 19; i++) {
        const double non_res_Xsu_correction_factor = (BR_Xsu_nonresonant_nunubar + BR_Kplusstar_nunubar - BR_Kplusstar_nunubar * Correction_factor_BR_Kstarplus[i]) / BR_Xsu_nonresonant_nunubar;
        GetNominalPDFs(MC_dirname_SIGNAL, "B2Xsnunu", CHG_Kstarff_array[i], "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu") * non_res_Xsu_correction_factor, "B2Xsnunu");
    }
    for (int i = 0; i < 19; i++) GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", MIX_Kstarff_array[i], "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetKstarffPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", MIX_Kstarff_array, Correction_factor_BR_Kstarzero, "Bzero", 0, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"));
    for (int i = 0; i < 19; i++) {
        const double non_res_Xsd_correction_factor = (BR_Xsd_nonresonant_nunubar + BR_K0star_nunubar - BR_K0star_nunubar * Correction_factor_BR_Kstarzero[i]) / BR_Xsd_nonresonant_nunubar;
        GetNominalPDFs(MC_dirname_SIGNAL, "B02Xsnunu", MIX_Kstarff_array[i], "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu") * non_res_Xsd_correction_factor, "B02Xsnunu");
    }

    for (int i = 0; i < 7; i++) {
        GetNominalPDFs(MC_dirname_CHG, "root", CHG_Kstarff_array[i], "Bplus", "CHG", ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
        GetNominalPDFs(MC_dirname_MIX, "root", MIX_Kstarff_array[i], "Bzero", "MIX", ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");
    }

    // get OLD Kff uncertainty PDFs
    GetOLDKffPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_Kff_OLD_p, "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    GetOLDKffPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", CHG_Kff_OLD_p, "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    GetOLDKffPDFs(MC_dirname_SIGNAL, "B2Xsnunu", CHG_Kff_OLD_p, "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    GetOLDKffPDFs(MC_dirname_SIGNAL, "B02K0nunu", MIX_Kff_OLD_p, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetOLDKffPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", MIX_Kff_OLD_p, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    GetOLDKffPDFs(MC_dirname_SIGNAL, "B02Xsnunu", MIX_Kff_OLD_p, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    GetNominalPDFs(MC_dirname_CHG, "root", CHG_Kff_OLD_p, "Bplus", "CHG", ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    GetNominalPDFs(MC_dirname_MIX, "root", MIX_Kff_OLD_p, "Bzero", "MIX", ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");

    GetNegativeChangePDFs(Signal_nominal, Signal_Kff_OLD_p, Signal_Kff_OLD_m);
    GetNegativeChangePDFs(CHG_nominal, CHG_Kff_OLD_p, CHG_Kff_OLD_m);
    GetNegativeChangePDFs(MIX_nominal, MIX_Kff_OLD_p, MIX_Kff_OLD_m);


    // get fraction uncertainty pdfs (correlated)
    int NPDFs_fraction = GetfractioncorrelatedPDFs(fraction_correlated_info, CHG_nominal, MIX_nominal, Signal_nominal, &CHG_fraction_correlated, &MIX_fraction_correlated, &Signal_fraction_correlated);

    // get fraction uncertainty pdfs (uncorrelated)
    GetfractionUncorrelatedPDFs(fraction_uncorrelated_info, CHG_fraction_uncorrelated, MIX_fraction_uncorrelated, Signal_fraction_uncorrelated);

    // get fragmentation uncertainty pdfs (correlated)
    int NPDFs_Fragmentation = GetFragmentationcorrelatedPDFs(Fragmentation_correlated_info, CHG_nominal, MIX_nominal, Signal_nominal, &CHG_Fragmentation_correlated, &MIX_Fragmentation_correlated, &Signal_Fragmentation_correlated);

    // get fragmentation uncertainty pdfs (uncorrelated)
    GetFragmentationUncorrelatedPDFs(Fragmentation_uncorrelated_info, CHG_Fragmentation_uncorrelated, MIX_Fragmentation_uncorrelated, Signal_Fragmentation_uncorrelated);

    // get pf uncertainty pdfs
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", CHG_pf_p, "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    GetpfPDFs(MC_dirname_SIGNAL, "B2Xsnunu", CHG_pf_p, "Bplus", "SIGNAL", true, true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", MIX_pf_p, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", MIX_pf_p, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    GetpfPDFs(MC_dirname_SIGNAL, "B02Xsnunu", MIX_pf_p, "Bzero", "SIGNAL", true, false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    GetNominalPDFs(MC_dirname_CHG, "root", CHG_pf_p, "Bplus", "CHG", ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    GetNominalPDFs(MC_dirname_MIX, "root", MIX_pf_p, "Bzero", "MIX", ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");

    GetNominalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", CHG_pf_m, "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    GetpfPDFs(MC_dirname_SIGNAL, "B2Xsnunu", CHG_pf_m, "Bplus", "SIGNAL", false, true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", MIX_pf_m, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", MIX_pf_m, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    GetpfPDFs(MC_dirname_SIGNAL, "B02Xsnunu", MIX_pf_m, "Bzero", "SIGNAL", false, false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    GetNominalPDFs(MC_dirname_CHG, "root", CHG_pf_m, "Bplus", "CHG", ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    GetNominalPDFs(MC_dirname_MIX, "root", MIX_pf_m, "Bzero", "MIX", ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");

    // mb
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", CHG_mb_p, "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    GetmbPDFs(MC_dirname_SIGNAL, "B2Xsnunu", CHG_mb_p, "Bplus", "SIGNAL", true, true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", MIX_mb_p, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", MIX_mb_p, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    GetmbPDFs(MC_dirname_SIGNAL, "B02Xsnunu", MIX_mb_p, "Bzero", "SIGNAL", true, false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    GetNominalPDFs(MC_dirname_CHG, "root", CHG_mb_p, "Bplus", "CHG", ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    GetNominalPDFs(MC_dirname_MIX, "root", MIX_mb_p, "Bzero", "MIX", ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");

    GetNominalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", CHG_mb_m, "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    GetmbPDFs(MC_dirname_SIGNAL, "B2Xsnunu", CHG_mb_m, "Bplus", "SIGNAL", false, true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", MIX_mb_m, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", MIX_mb_m, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    GetmbPDFs(MC_dirname_SIGNAL, "B02Xsnunu", MIX_mb_m, "Bzero", "SIGNAL", false, false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    GetNominalPDFs(MC_dirname_CHG, "root", CHG_mb_m, "Bplus", "CHG", ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    GetNominalPDFs(MC_dirname_MIX, "root", MIX_mb_m, "Bzero", "MIX", ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");

    // transition
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", CHG_transition_p, "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    GetTransitionPDFs(MC_dirname_SIGNAL, "B2Xsnunu", CHG_transition_p, "Bplus", "SIGNAL", true, true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", MIX_transition_p, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", MIX_transition_p, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    GetTransitionPDFs(MC_dirname_SIGNAL, "B02Xsnunu", MIX_transition_p, "Bzero", "SIGNAL", true, false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    GetNominalPDFs(MC_dirname_CHG, "root", CHG_transition_p, "Bplus", "CHG", ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    GetNominalPDFs(MC_dirname_MIX, "root", MIX_transition_p, "Bzero", "MIX", ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");

    GetNominalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", CHG_transition_m, "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    GetTransitionPDFs(MC_dirname_SIGNAL, "B2Xsnunu", CHG_transition_m, "Bplus", "SIGNAL", false, true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", MIX_transition_m, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", MIX_transition_m, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    GetTransitionPDFs(MC_dirname_SIGNAL, "B02Xsnunu", MIX_transition_m, "Bzero", "SIGNAL", false, false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    GetNominalPDFs(MC_dirname_CHG, "root", CHG_transition_m, "Bplus", "CHG", ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    GetNominalPDFs(MC_dirname_MIX, "root", MIX_transition_m, "Bzero", "MIX", ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");


    // mKstar fixed
    GetKstardeltaPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", CHG_mKstar_p, "Bplus", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Xsnunu", CHG_mKstar_p, "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", MIX_mKstar_p, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetKstardeltaPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", MIX_mKstar_p, "Bzero", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Xsnunu", MIX_mKstar_p, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    GetNegativeChangePDFs(CHG_nominal, CHG_mKstar_p, CHG_mKstar_m);
    GetNegativeChangePDFs(MIX_nominal, MIX_mKstar_p, MIX_mKstar_m);


    // get BDTc uncertainty pdfs
    //double NevtCHGWithBDTc = GetNevtWithBDTc(MC_dirname_CHG, "root", "Bplus", "CHG", ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    //double NevtMIXWithBDTc = GetNevtWithBDTc(MC_dirname_MIX, "root", "Bzero", "MIX", ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");
    double NevtUUBARWithBDTc = GetNevtWithBDTc(MC_dirname_UUBAR, "root", "Continuum", "UUBAR", ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"), "otherwise");
    double NevtDDBARWithBDTc = GetNevtWithBDTc(MC_dirname_DDBAR, "root", "Continuum", "DDBAR", ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"), "otherwise");
    double NevtSSBARWithBDTc = GetNevtWithBDTc(MC_dirname_SSBAR, "root", "Continuum", "SSBAR", ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"), "otherwise");
    double NevtCHARMWithBDTc = GetNevtWithBDTc(MC_dirname_CHARM, "root", "Continuum", "CHARM", ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"), "otherwise");

    double NormFactor = (NevtUUBARWithoutBDTc + NevtDDBARWithoutBDTc + NevtSSBARWithoutBDTc + NevtCHARMWithoutBDTc) / (NevtUUBARWithBDTc + NevtDDBARWithBDTc + NevtSSBARWithBDTc + NevtCHARMWithBDTc);

    //GetBDTcPDFs(MC_dirname_CHG, "root", CHG_BDTc_p, "Bplus", "CHG", ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    //GetBDTcPDFs(MC_dirname_MIX, "root", MIX_BDTc_p, "Bzero", "MIX", ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");
    double NevtUUBARWithBDTcWithNorm = GetBDTcPDFs(MC_dirname_UUBAR, "root", UUBAR_BDTc_p, "Continuum", "UUBAR", ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR") * NormFactor, "otherwise");
    double NevtDDBARWithBDTcWithNorm = GetBDTcPDFs(MC_dirname_DDBAR, "root", DDBAR_BDTc_p, "Continuum", "DDBAR", ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR") * NormFactor, "otherwise");
    double NevtSSBARWithBDTcWithNorm = GetBDTcPDFs(MC_dirname_SSBAR, "root", SSBAR_BDTc_p, "Continuum", "SSBAR", ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR") * NormFactor, "otherwise");
    double NevtCHARMWithBDTcWithNorm = GetBDTcPDFs(MC_dirname_CHARM, "root", CHARM_BDTc_p, "Continuum", "CHARM", ObtainWeight("CHARM", MCTYPE, "validation", "CHARM") * NormFactor, "otherwise");

    //GetNegativeChangePDFs(CHG_nominal, CHG_BDTc_p, CHG_BDTc_m);
    //GetNegativeChangePDFs(MIX_nominal, MIX_BDTc_p, MIX_BDTc_m);
    GetNegativeChangePDFs(UUBAR_nominal, UUBAR_BDTc_p, UUBAR_BDTc_m);
    GetNegativeChangePDFs(DDBAR_nominal, DDBAR_BDTc_p, DDBAR_BDTc_m);
    GetNegativeChangePDFs(SSBAR_nominal, SSBAR_BDTc_p, SSBAR_BDTc_m);
    GetNegativeChangePDFs(CHARM_nominal, CHARM_BDTc_p, CHARM_BDTc_m);

    if (std::abs((NevtUUBARWithoutBDTc + NevtDDBARWithoutBDTc + NevtSSBARWithoutBDTc + NevtCHARMWithoutBDTc) - (NevtUUBARWithBDTcWithNorm + NevtDDBARWithBDTcWithNorm + NevtSSBARWithBDTcWithNorm + NevtCHARMWithBDTcWithNorm)) > 0.1) { // final check
        printf("ERROR! before and after BDTc correction shows different number of event!\n");
        exit(1);
    }

    // B2Xnn BR
    GetXnnBRPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_Xnn_p, "Bplus", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    GetXnnBRPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", CHG_Xnn_p, "Bplus", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    GetXnnBRPDFs(MC_dirname_SIGNAL, "B2Xsnunu", CHG_Xnn_p, "Bplus", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    GetXnnBRPDFs(MC_dirname_SIGNAL, "B02K0nunu", MIX_Xnn_p, "Bzero", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetXnnBRPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", MIX_Xnn_p, "Bzero", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    GetXnnBRPDFs(MC_dirname_SIGNAL, "B02Xsnunu", MIX_Xnn_p, "Bzero", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    GetXnnBRPDFs(MC_dirname_CHG, "root", CHG_Xnn_p, "Bplus", "CHG", true, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    GetXnnBRPDFs(MC_dirname_MIX, "root", MIX_Xnn_p, "Bzero", "MIX", true, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");

    GetXnnBRPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_Xnn_m, "Bplus", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    GetXnnBRPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", CHG_Xnn_m, "Bplus", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    GetXnnBRPDFs(MC_dirname_SIGNAL, "B2Xsnunu", CHG_Xnn_m, "Bplus", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    GetXnnBRPDFs(MC_dirname_SIGNAL, "B02K0nunu", MIX_Xnn_m, "Bzero", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetXnnBRPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", MIX_Xnn_m, "Bzero", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    GetXnnBRPDFs(MC_dirname_SIGNAL, "B02Xsnunu", MIX_Xnn_m, "Bzero", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    GetXnnBRPDFs(MC_dirname_CHG, "root", CHG_Xnn_m, "Bplus", "CHG", false, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    GetXnnBRPDFs(MC_dirname_MIX, "root", MIX_Xnn_m, "Bzero", "MIX", false, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");

    // get multiplicity uncertainty pdfs (correlated)
    int NPDFs_multiplicity = GetmultiplicitycorrelatedPDFs(multiplicity_correlated_info, CHG_nominal, MIX_nominal, UUBAR_nominal, DDBAR_nominal, SSBAR_nominal, CHARM_nominal, Signal_nominal, &CHG_multiplicity_correlated, &MIX_multiplicity_correlated, &UUBAR_multiplicity_correlated, &DDBAR_multiplicity_correlated, &SSBAR_multiplicity_correlated, &CHARM_multiplicity_correlated, &Signal_multiplicity_correlated);

    // get multiplicity uncertainty pdfs (uncorrelated)
    GetmultiplicityUncorrelatedPDFs(multiplicity_uncorrelated_info, CHG_multiplicity_uncorrelated, MIX_multiplicity_uncorrelated, UUBAR_multiplicity_uncorrelated, DDBAR_multiplicity_uncorrelated, SSBAR_multiplicity_uncorrelated, CHARM_multiplicity_uncorrelated, Signal_multiplicity_uncorrelated);

    // B -> [D -> X KL0] anything
    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_BtoDtoXKL_p, "Bplus", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", CHG_BtoDtoXKL_p, "Bplus", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B2Xsnunu", CHG_BtoDtoXKL_p, "Bplus", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B02K0nunu", MIX_BtoDtoXKL_p, "Bzero", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", MIX_BtoDtoXKL_p, "Bzero", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B02Xsnunu", MIX_BtoDtoXKL_p, "Bzero", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    GetBtoDtoXKLPDFs(MC_dirname_CHG, "root", CHG_BtoDtoXKL_p, "Bplus", "CHG", true, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    GetBtoDtoXKLPDFs(MC_dirname_MIX, "root", MIX_BtoDtoXKL_p, "Bzero", "MIX", true, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");

    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_BtoDtoXKL_m, "Bplus", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", CHG_BtoDtoXKL_m, "Bplus", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B2Xsnunu", CHG_BtoDtoXKL_m, "Bplus", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B02K0nunu", MIX_BtoDtoXKL_m, "Bzero", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", MIX_BtoDtoXKL_m, "Bzero", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    GetBtoDtoXKLPDFs(MC_dirname_SIGNAL, "B02Xsnunu", MIX_BtoDtoXKL_m, "Bzero", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    GetBtoDtoXKLPDFs(MC_dirname_CHG, "root", CHG_BtoDtoXKL_m, "Bplus", "CHG", false, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    GetBtoDtoXKLPDFs(MC_dirname_MIX, "root", MIX_BtoDtoXKL_m, "Bzero", "MIX", false, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");

    // BR(B -> X KL KL)
    GetBRXKLKLPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_BRBtoXKLKL_p, "Bplus", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    GetBRXKLKLPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", CHG_BRBtoXKLKL_p, "Bplus", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    GetBRXKLKLPDFs(MC_dirname_SIGNAL, "B2Xsnunu", CHG_BRBtoXKLKL_p, "Bplus", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    GetBRXKLKLPDFs(MC_dirname_SIGNAL, "B02K0nunu", MIX_BRBtoXKLKL_p, "Bzero", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetBRXKLKLPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", MIX_BRBtoXKLKL_p, "Bzero", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    GetBRXKLKLPDFs(MC_dirname_SIGNAL, "B02Xsnunu", MIX_BRBtoXKLKL_p, "Bzero", "SIGNAL", true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    GetBRXKLKLPDFs(MC_dirname_CHG, "root", CHG_BRBtoXKLKL_p, "Bplus", "CHG", true, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    GetBRXKLKLPDFs(MC_dirname_MIX, "root", MIX_BRBtoXKLKL_p, "Bzero", "MIX", true, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");

    GetBRXKLKLPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_BRBtoXKLKL_m, "Bplus", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    GetBRXKLKLPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", CHG_BRBtoXKLKL_m, "Bplus", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    GetBRXKLKLPDFs(MC_dirname_SIGNAL, "B2Xsnunu", CHG_BRBtoXKLKL_m, "Bplus", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    GetBRXKLKLPDFs(MC_dirname_SIGNAL, "B02K0nunu", MIX_BRBtoXKLKL_m, "Bzero", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetBRXKLKLPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", MIX_BRBtoXKLKL_m, "Bzero", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    GetBRXKLKLPDFs(MC_dirname_SIGNAL, "B02Xsnunu", MIX_BRBtoXKLKL_m, "Bzero", "SIGNAL", false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    GetBRXKLKLPDFs(MC_dirname_CHG, "root", CHG_BRBtoXKLKL_m, "Bplus", "CHG", false, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    GetBRXKLKLPDFs(MC_dirname_MIX, "root", MIX_BRBtoXKLKL_m, "Bzero", "MIX", false, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");

    // calculate all uncorrelated pdfs
    ClearHist(Signal_all_uncorrelated, CHG_all_uncorrelated, MIX_all_uncorrelated, UUBAR_all_uncorrelated, DDBAR_all_uncorrelated, SSBAR_all_uncorrelated, CHARM_all_uncorrelated);

    AddSQRTHist(Signal_all_uncorrelated, Signal_FEI_uncorrelated, RarityBins);
    AddSQRTHist(CHG_all_uncorrelated, CHG_FEI_uncorrelated, RarityBins);
    AddSQRTHist(MIX_all_uncorrelated, MIX_FEI_uncorrelated, RarityBins);

    AddSQRTHist(Signal_all_uncorrelated, Signal_KID_uncorrelated, RarityBins);
    AddSQRTHist(CHG_all_uncorrelated, CHG_KID_uncorrelated, RarityBins);
    AddSQRTHist(MIX_all_uncorrelated, MIX_KID_uncorrelated, RarityBins);
    AddSQRTHist(UUBAR_all_uncorrelated, UUBAR_KID_uncorrelated, RarityBins);
    AddSQRTHist(DDBAR_all_uncorrelated, DDBAR_KID_uncorrelated, RarityBins);
    AddSQRTHist(SSBAR_all_uncorrelated, SSBAR_KID_uncorrelated, RarityBins);
    AddSQRTHist(CHARM_all_uncorrelated, CHARM_KID_uncorrelated, RarityBins);

    AddSQRTHist(Signal_all_uncorrelated, Signal_PID_uncorrelated, RarityBins);
    AddSQRTHist(CHG_all_uncorrelated, CHG_PID_uncorrelated, RarityBins);
    AddSQRTHist(MIX_all_uncorrelated, MIX_PID_uncorrelated, RarityBins);
    AddSQRTHist(UUBAR_all_uncorrelated, UUBAR_PID_uncorrelated, RarityBins);
    AddSQRTHist(DDBAR_all_uncorrelated, DDBAR_PID_uncorrelated, RarityBins);
    AddSQRTHist(SSBAR_all_uncorrelated, SSBAR_PID_uncorrelated, RarityBins);
    AddSQRTHist(CHARM_all_uncorrelated, CHARM_PID_uncorrelated, RarityBins);

    AddSQRTHist(Signal_all_uncorrelated, Signal_BR_uncorrelated, RarityBins);
    AddSQRTHist(CHG_all_uncorrelated, CHG_BR_uncorrelated, RarityBins);
    AddSQRTHist(MIX_all_uncorrelated, MIX_BR_uncorrelated, RarityBins);

    AddSQRTHist(Signal_all_uncorrelated, Signal_pi0_uncorrelated, RarityBins);
    AddSQRTHist(CHG_all_uncorrelated, CHG_pi0_uncorrelated, RarityBins);
    AddSQRTHist(MIX_all_uncorrelated, MIX_pi0_uncorrelated, RarityBins);
    AddSQRTHist(UUBAR_all_uncorrelated, UUBAR_pi0_uncorrelated, RarityBins);
    AddSQRTHist(DDBAR_all_uncorrelated, DDBAR_pi0_uncorrelated, RarityBins);
    AddSQRTHist(SSBAR_all_uncorrelated, SSBAR_pi0_uncorrelated, RarityBins);
    AddSQRTHist(CHARM_all_uncorrelated, CHARM_pi0_uncorrelated, RarityBins);

    AddSQRTHist(Signal_all_uncorrelated, Signal_multiplicity_uncorrelated, RarityBins);
    AddSQRTHist(CHG_all_uncorrelated, CHG_multiplicity_uncorrelated, RarityBins);
    AddSQRTHist(MIX_all_uncorrelated, MIX_multiplicity_uncorrelated, RarityBins);
    AddSQRTHist(UUBAR_all_uncorrelated, UUBAR_multiplicity_uncorrelated, RarityBins);
    AddSQRTHist(DDBAR_all_uncorrelated, DDBAR_multiplicity_uncorrelated, RarityBins);
    AddSQRTHist(SSBAR_all_uncorrelated, SSBAR_multiplicity_uncorrelated, RarityBins);
    AddSQRTHist(CHARM_all_uncorrelated, CHARM_multiplicity_uncorrelated, RarityBins);

    AddSQRTHist(Signal_all_uncorrelated, Signal_fraction_uncorrelated, RarityBins);
    AddSQRTHist(CHG_all_uncorrelated, CHG_fraction_uncorrelated, RarityBins);
    AddSQRTHist(MIX_all_uncorrelated, MIX_fraction_uncorrelated, RarityBins);

    AddSQRTHist(Signal_all_uncorrelated, Signal_Fragmentation_uncorrelated, RarityBins);
    AddSQRTHist(CHG_all_uncorrelated, CHG_Fragmentation_uncorrelated, RarityBins);
    AddSQRTHist(MIX_all_uncorrelated, MIX_Fragmentation_uncorrelated, RarityBins);

    // calculate MC statistical uncertainties (relative errors)
    GetMCstatisticalRelativeError(Signal_nominal, Signal_MC_stat, RarityBins);
    GetMCstatisticalRelativeError(CHG_nominal, CHG_MC_stat, RarityBins);
    GetMCstatisticalRelativeError(MIX_nominal, MIX_MC_stat, RarityBins);
    GetMCstatisticalRelativeError(UUBAR_nominal, UUBAR_MC_stat, RarityBins);
    GetMCstatisticalRelativeError(DDBAR_nominal, DDBAR_MC_stat, RarityBins);
    GetMCstatisticalRelativeError(SSBAR_nominal, SSBAR_MC_stat, RarityBins);
    GetMCstatisticalRelativeError(CHARM_nominal, CHARM_MC_stat, RarityBins);

    // all of uncorrelated uncertainties + MC statistical uncertainties
    AddSQRTHist(Signal_all_uncorrelated_MC_stat, Signal_all_uncorrelated, Signal_MC_stat, RarityBins);
    AddSQRTHist(CHG_all_uncorrelated_MC_stat, CHG_all_uncorrelated, CHG_MC_stat, RarityBins);
    AddSQRTHist(MIX_all_uncorrelated_MC_stat, MIX_all_uncorrelated, MIX_MC_stat, RarityBins);
    AddSQRTHist(UUBAR_all_uncorrelated_MC_stat, UUBAR_all_uncorrelated, UUBAR_MC_stat, RarityBins);
    AddSQRTHist(DDBAR_all_uncorrelated_MC_stat, DDBAR_all_uncorrelated, DDBAR_MC_stat, RarityBins);
    AddSQRTHist(SSBAR_all_uncorrelated_MC_stat, SSBAR_all_uncorrelated, SSBAR_MC_stat, RarityBins);
    AddSQRTHist(CHARM_all_uncorrelated_MC_stat, CHARM_all_uncorrelated, CHARM_MC_stat, RarityBins);
    /* ====================================== */



    /* ====================================== */
    // define path for Data
    const char* DATA_dirname_Knunu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Knunu";
    const char* DATA_dirname_Kstarnunu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Kstarnunu";
    const char* DATA_dirname_Xsununu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Xsnunu";
    const char* DATA_dirname_K0nunu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0nunu";
    const char* DATA_dirname_K0starnunu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0starnunu";
    const char* DATA_dirname_Xsdnunu = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02Xsnunu";

    const char* DATA_dirname_CHG = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/CHG_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_MIX = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/MIX_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_UUBAR = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/UUBAR_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_DDBAR = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/DDBAR_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_SSBAR = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SSBAR_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_CHARM = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/CHARM_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    /* ====================================== */



    /* ====================================== */
    // Define histograms for Data
    TH1D* total_DATA = new TH1D("total_DATA", "total_DATA", RarityBins, BinMIN, BinMAX);
    /* ====================================== */



    /* ====================================== */
    // Save data
    GetNominalPDFs(DATA_dirname_Knunu, "root", total_DATA, "Bplus", "DATA", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"));
    GetNominalPDFs(DATA_dirname_Kstarnunu, "root", total_DATA, "Bplus", "DATA", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"));
    GetNominalPDFs(DATA_dirname_Xsununu, "root", total_DATA, "Bplus", "DATA", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"));
    GetNominalPDFs(DATA_dirname_K0nunu, "root", total_DATA, "Bzero", "DATA", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"));
    GetNominalPDFs(DATA_dirname_K0starnunu, "root", total_DATA, "Bzero", "DATA", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"));
    GetNominalPDFs(DATA_dirname_Xsdnunu, "root", total_DATA, "Bzero", "DATA", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"));

    GetNominalPDFs(DATA_dirname_CHG, "root", total_DATA, "Bplus", "DATA", 1.0);
    GetNominalPDFs(DATA_dirname_MIX, "root", total_DATA, "Bzero", "DATA", 1.0);
    GetNominalPDFs(DATA_dirname_UUBAR, "root", total_DATA, "Continuum", "DATA", 1.0);
    GetNominalPDFs(DATA_dirname_DDBAR, "root", total_DATA, "Continuum", "DATA", 1.0);
    GetNominalPDFs(DATA_dirname_SSBAR, "root", total_DATA, "Continuum", "DATA", 1.0);
    GetNominalPDFs(DATA_dirname_CHARM, "root", total_DATA, "Continuum", "DATA", 1.0);
    /* ====================================== */



    /* ====================================== */
    // Select specific MXs bin
    const int MXsBin = 1;

    SaveSpecificMXsBin(Signal_nominal, MXsBin);
    SaveSpecificMXsBin(CHG_nominal, MXsBin);
    SaveSpecificMXsBin(MIX_nominal, MXsBin);
    SaveSpecificMXsBin(UUBAR_nominal, MXsBin);
    SaveSpecificMXsBin(DDBAR_nominal, MXsBin);
    SaveSpecificMXsBin(SSBAR_nominal, MXsBin);
    SaveSpecificMXsBin(CHARM_nominal, MXsBin);

    // track uncertainty
    SaveSpecificMXsBin(Signal_track_p, MXsBin);
    SaveSpecificMXsBin(CHG_track_p, MXsBin);
    SaveSpecificMXsBin(MIX_track_p, MXsBin);
    SaveSpecificMXsBin(UUBAR_track_p, MXsBin);
    SaveSpecificMXsBin(DDBAR_track_p, MXsBin);
    SaveSpecificMXsBin(SSBAR_track_p, MXsBin);
    SaveSpecificMXsBin(CHARM_track_p, MXsBin);
    SaveSpecificMXsBin(Signal_track_m, MXsBin);
    SaveSpecificMXsBin(CHG_track_m, MXsBin);
    SaveSpecificMXsBin(MIX_track_m, MXsBin);
    SaveSpecificMXsBin(UUBAR_track_m, MXsBin);
    SaveSpecificMXsBin(DDBAR_track_m, MXsBin);
    SaveSpecificMXsBin(SSBAR_track_m, MXsBin);
    SaveSpecificMXsBin(CHARM_track_m, MXsBin);

    // KS0 uncertainty
    SaveSpecificMXsBin(Signal_KS0_p, MXsBin);
    SaveSpecificMXsBin(CHG_KS0_p, MXsBin);
    SaveSpecificMXsBin(MIX_KS0_p, MXsBin);
    SaveSpecificMXsBin(UUBAR_KS0_p, MXsBin);
    SaveSpecificMXsBin(DDBAR_KS0_p, MXsBin);
    SaveSpecificMXsBin(SSBAR_KS0_p, MXsBin);
    SaveSpecificMXsBin(CHARM_KS0_p, MXsBin);
    SaveSpecificMXsBin(Signal_KS0_m, MXsBin);
    SaveSpecificMXsBin(CHG_KS0_m, MXsBin);
    SaveSpecificMXsBin(MIX_KS0_m, MXsBin);
    SaveSpecificMXsBin(UUBAR_KS0_m, MXsBin);
    SaveSpecificMXsBin(DDBAR_KS0_m, MXsBin);
    SaveSpecificMXsBin(SSBAR_KS0_m, MXsBin);
    SaveSpecificMXsBin(CHARM_KS0_m, MXsBin);

    // FEI uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_FEI; i++) {
        SaveSpecificMXsBin(Signal_FEI_correlated[i], MXsBin);
        SaveSpecificMXsBin(CHG_FEI_correlated[i], MXsBin);
        SaveSpecificMXsBin(MIX_FEI_correlated[i], MXsBin);
    }

    // FEI uncertainty (uncorrelated)
    SaveSpecificMXsBin(Signal_FEI_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHG_FEI_uncorrelated, MXsBin);
    SaveSpecificMXsBin(MIX_FEI_uncorrelated, MXsBin);

    // Kaon PID uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_KID; i++) {
        SaveSpecificMXsBin(Signal_KID_correlated[i], MXsBin);
        SaveSpecificMXsBin(CHG_KID_correlated[i], MXsBin);
        SaveSpecificMXsBin(MIX_KID_correlated[i], MXsBin);
        SaveSpecificMXsBin(UUBAR_KID_correlated[i], MXsBin);
        SaveSpecificMXsBin(DDBAR_KID_correlated[i], MXsBin);
        SaveSpecificMXsBin(SSBAR_KID_correlated[i], MXsBin);
        SaveSpecificMXsBin(CHARM_KID_correlated[i], MXsBin);
    }

    // Kaon PID uncertainty (uncorrelated)
    SaveSpecificMXsBin(Signal_KID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHG_KID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(MIX_KID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(UUBAR_KID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(DDBAR_KID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(SSBAR_KID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHARM_KID_uncorrelated, MXsBin);

    // Pion PID uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_PID; i++) {
        SaveSpecificMXsBin(Signal_PID_correlated[i], MXsBin);
        SaveSpecificMXsBin(CHG_PID_correlated[i], MXsBin);
        SaveSpecificMXsBin(MIX_PID_correlated[i], MXsBin);
        SaveSpecificMXsBin(UUBAR_PID_correlated[i], MXsBin);
        SaveSpecificMXsBin(DDBAR_PID_correlated[i], MXsBin);
        SaveSpecificMXsBin(SSBAR_PID_correlated[i], MXsBin);
        SaveSpecificMXsBin(CHARM_PID_correlated[i], MXsBin);
    }

    // Pion PID uncertainty (uncorrelated)
    SaveSpecificMXsBin(Signal_PID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHG_PID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(MIX_PID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(UUBAR_PID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(DDBAR_PID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(SSBAR_PID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHARM_PID_uncorrelated, MXsBin);

    // BB BR uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_BR; i++) {
        SaveSpecificMXsBin(Signal_BR_correlated[i], MXsBin);
        SaveSpecificMXsBin(CHG_BR_correlated[i], MXsBin);
        SaveSpecificMXsBin(MIX_BR_correlated[i], MXsBin);
    }

    // BB BR uncertainty (uncorrelated)
    SaveSpecificMXsBin(Signal_BR_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHG_BR_uncorrelated, MXsBin);
    SaveSpecificMXsBin(MIX_BR_uncorrelated, MXsBin);

    // pi0 uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_pi0; i++) {
        SaveSpecificMXsBin(Signal_pi0_correlated[i], MXsBin);
        SaveSpecificMXsBin(CHG_pi0_correlated[i], MXsBin);
        SaveSpecificMXsBin(MIX_pi0_correlated[i], MXsBin);
        SaveSpecificMXsBin(UUBAR_pi0_correlated[i], MXsBin);
        SaveSpecificMXsBin(DDBAR_pi0_correlated[i], MXsBin);
        SaveSpecificMXsBin(SSBAR_pi0_correlated[i], MXsBin);
        SaveSpecificMXsBin(CHARM_pi0_correlated[i], MXsBin);
    }

    // pi0 uncertainty (uncorrelated)
    SaveSpecificMXsBin(Signal_pi0_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHG_pi0_uncorrelated, MXsBin);
    SaveSpecificMXsBin(MIX_pi0_uncorrelated, MXsBin);
    SaveSpecificMXsBin(UUBAR_pi0_uncorrelated, MXsBin);
    SaveSpecificMXsBin(DDBAR_pi0_uncorrelated, MXsBin);
    SaveSpecificMXsBin(SSBAR_pi0_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHARM_pi0_uncorrelated, MXsBin);

    // K nu nubar form factor
    SaveSpecificMXsBin(Signal_Kff1_p, MXsBin);
    SaveSpecificMXsBin(Signal_Kff1_m, MXsBin);
    SaveSpecificMXsBin(Signal_Kff2_p, MXsBin);
    SaveSpecificMXsBin(Signal_Kff2_m, MXsBin);
    SaveSpecificMXsBin(Signal_Kff3_p, MXsBin);
    SaveSpecificMXsBin(Signal_Kff3_m, MXsBin);
    SaveSpecificMXsBin(CHG_Kff1_p, MXsBin);
    SaveSpecificMXsBin(CHG_Kff1_m, MXsBin);
    SaveSpecificMXsBin(CHG_Kff2_p, MXsBin);
    SaveSpecificMXsBin(CHG_Kff2_m, MXsBin);
    SaveSpecificMXsBin(CHG_Kff3_p, MXsBin);
    SaveSpecificMXsBin(CHG_Kff3_m, MXsBin);
    SaveSpecificMXsBin(MIX_Kff1_p, MXsBin);
    SaveSpecificMXsBin(MIX_Kff1_m, MXsBin);
    SaveSpecificMXsBin(MIX_Kff2_p, MXsBin);
    SaveSpecificMXsBin(MIX_Kff2_m, MXsBin);
    SaveSpecificMXsBin(MIX_Kff3_p, MXsBin);
    SaveSpecificMXsBin(MIX_Kff3_m, MXsBin);

    // Kstar nu nubar form factor
    SaveSpecificMXsBin(CHG_Kstarff1_p, MXsBin);
    SaveSpecificMXsBin(CHG_Kstarff1_m, MXsBin);
    SaveSpecificMXsBin(CHG_Kstarff2_p, MXsBin);
    SaveSpecificMXsBin(CHG_Kstarff2_m, MXsBin);
    SaveSpecificMXsBin(CHG_Kstarff3_p, MXsBin);
    SaveSpecificMXsBin(CHG_Kstarff3_m, MXsBin);
    SaveSpecificMXsBin(CHG_Kstarff4_p, MXsBin);
    SaveSpecificMXsBin(CHG_Kstarff4_m, MXsBin);
    SaveSpecificMXsBin(CHG_Kstarff5_p, MXsBin);
    SaveSpecificMXsBin(CHG_Kstarff5_m, MXsBin);
    SaveSpecificMXsBin(CHG_Kstarff6_p, MXsBin);
    SaveSpecificMXsBin(CHG_Kstarff6_m, MXsBin);
    SaveSpecificMXsBin(CHG_Kstarff7_p, MXsBin);
    SaveSpecificMXsBin(CHG_Kstarff7_m, MXsBin);
    SaveSpecificMXsBin(CHG_Kstarff8_p, MXsBin);
    SaveSpecificMXsBin(CHG_Kstarff8_m, MXsBin);
    SaveSpecificMXsBin(CHG_Kstarff9_p, MXsBin);
    SaveSpecificMXsBin(CHG_Kstarff9_m, MXsBin);
    SaveSpecificMXsBin(MIX_Kstarff1_p, MXsBin);
    SaveSpecificMXsBin(MIX_Kstarff1_m, MXsBin);
    SaveSpecificMXsBin(MIX_Kstarff2_p, MXsBin);
    SaveSpecificMXsBin(MIX_Kstarff2_m, MXsBin);
    SaveSpecificMXsBin(MIX_Kstarff3_p, MXsBin);
    SaveSpecificMXsBin(MIX_Kstarff3_m, MXsBin);
    SaveSpecificMXsBin(MIX_Kstarff4_p, MXsBin);
    SaveSpecificMXsBin(MIX_Kstarff4_m, MXsBin);
    SaveSpecificMXsBin(MIX_Kstarff5_p, MXsBin);
    SaveSpecificMXsBin(MIX_Kstarff5_m, MXsBin);
    SaveSpecificMXsBin(MIX_Kstarff6_p, MXsBin);
    SaveSpecificMXsBin(MIX_Kstarff6_m, MXsBin);
    SaveSpecificMXsBin(MIX_Kstarff7_p, MXsBin);
    SaveSpecificMXsBin(MIX_Kstarff7_m, MXsBin);
    SaveSpecificMXsBin(MIX_Kstarff8_p, MXsBin);
    SaveSpecificMXsBin(MIX_Kstarff8_m, MXsBin);
    SaveSpecificMXsBin(MIX_Kstarff9_p, MXsBin);
    SaveSpecificMXsBin(MIX_Kstarff9_m, MXsBin);

    // K nu nubar OLD form factor
    SaveSpecificMXsBin(Signal_Kff_OLD_p, MXsBin);
    SaveSpecificMXsBin(CHG_Kff_OLD_p, MXsBin);
    SaveSpecificMXsBin(MIX_Kff_OLD_p, MXsBin);
    SaveSpecificMXsBin(Signal_Kff_OLD_m, MXsBin);
    SaveSpecificMXsBin(CHG_Kff_OLD_m, MXsBin);
    SaveSpecificMXsBin(MIX_Kff_OLD_m, MXsBin);

    // fraction uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_fraction; i++) {
        SaveSpecificMXsBin(Signal_fraction_correlated[i], MXsBin);
        SaveSpecificMXsBin(CHG_fraction_correlated[i], MXsBin);
        SaveSpecificMXsBin(MIX_fraction_correlated[i], MXsBin);
    }

    // fraction uncertainty (uncorrelated)
    SaveSpecificMXsBin(Signal_fraction_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHG_fraction_uncorrelated, MXsBin);
    SaveSpecificMXsBin(MIX_fraction_uncorrelated, MXsBin);

    // fragmentation uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_Fragmentation; i++) {
        SaveSpecificMXsBin(Signal_Fragmentation_correlated[i], MXsBin);
        SaveSpecificMXsBin(CHG_Fragmentation_correlated[i], MXsBin);
        SaveSpecificMXsBin(MIX_Fragmentation_correlated[i], MXsBin);
    }

    // fragmentation uncertainty (uncorrelated)
    SaveSpecificMXsBin(Signal_Fragmentation_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHG_Fragmentation_uncorrelated, MXsBin);
    SaveSpecificMXsBin(MIX_Fragmentation_uncorrelated, MXsBin);

    // pf uncertainty
    SaveSpecificMXsBin(CHG_pf_p, MXsBin);
    SaveSpecificMXsBin(MIX_pf_p, MXsBin);
    SaveSpecificMXsBin(CHG_pf_m, MXsBin);
    SaveSpecificMXsBin(MIX_pf_m, MXsBin);

    // mb uncertainty
    SaveSpecificMXsBin(CHG_mb_p, MXsBin);
    SaveSpecificMXsBin(MIX_mb_p, MXsBin);
    SaveSpecificMXsBin(CHG_mb_m, MXsBin);
    SaveSpecificMXsBin(MIX_mb_m, MXsBin);

    // transition uncertainty
    SaveSpecificMXsBin(CHG_transition_p, MXsBin);
    SaveSpecificMXsBin(MIX_transition_p, MXsBin);
    SaveSpecificMXsBin(CHG_transition_m, MXsBin);
    SaveSpecificMXsBin(MIX_transition_m, MXsBin);

    // fixed mKstar uncertainty
    SaveSpecificMXsBin(CHG_mKstar_p, MXsBin);
    SaveSpecificMXsBin(MIX_mKstar_p, MXsBin);
    SaveSpecificMXsBin(CHG_mKstar_m, MXsBin);
    SaveSpecificMXsBin(MIX_mKstar_m, MXsBin);

    // BDTc
    SaveSpecificMXsBin(UUBAR_BDTc_p, MXsBin);
    SaveSpecificMXsBin(DDBAR_BDTc_p, MXsBin);
    SaveSpecificMXsBin(SSBAR_BDTc_p, MXsBin);
    SaveSpecificMXsBin(CHARM_BDTc_p, MXsBin);
    SaveSpecificMXsBin(UUBAR_BDTc_m, MXsBin);
    SaveSpecificMXsBin(DDBAR_BDTc_m, MXsBin);
    SaveSpecificMXsBin(SSBAR_BDTc_m, MXsBin);
    SaveSpecificMXsBin(CHARM_BDTc_m, MXsBin);

    // B2Knn BR
    SaveSpecificMXsBin(Signal_Xnn_p, MXsBin);
    SaveSpecificMXsBin(CHG_Xnn_p, MXsBin);
    SaveSpecificMXsBin(MIX_Xnn_p, MXsBin);
    SaveSpecificMXsBin(Signal_Xnn_m, MXsBin);
    SaveSpecificMXsBin(CHG_Xnn_m, MXsBin);
    SaveSpecificMXsBin(MIX_Xnn_m, MXsBin);

    // multiplicity uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_multiplicity; i++) {
        SaveSpecificMXsBin(Signal_multiplicity_correlated[i], MXsBin);
        SaveSpecificMXsBin(CHG_multiplicity_correlated[i], MXsBin);
        SaveSpecificMXsBin(MIX_multiplicity_correlated[i], MXsBin);
        SaveSpecificMXsBin(UUBAR_multiplicity_correlated[i], MXsBin);
        SaveSpecificMXsBin(DDBAR_multiplicity_correlated[i], MXsBin);
        SaveSpecificMXsBin(SSBAR_multiplicity_correlated[i], MXsBin);
        SaveSpecificMXsBin(CHARM_multiplicity_correlated[i], MXsBin);
    }

    // multiplicity uncertainty (uncorrelated)
    SaveSpecificMXsBin(Signal_multiplicity_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHG_multiplicity_uncorrelated, MXsBin);
    SaveSpecificMXsBin(MIX_multiplicity_uncorrelated, MXsBin);
    SaveSpecificMXsBin(UUBAR_multiplicity_uncorrelated, MXsBin);
    SaveSpecificMXsBin(DDBAR_multiplicity_uncorrelated, MXsBin);
    SaveSpecificMXsBin(SSBAR_multiplicity_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHARM_multiplicity_uncorrelated, MXsBin);

    // B -> [D -> X KL0] anything uncertainties
    SaveSpecificMXsBin(Signal_BtoDtoXKL_p, MXsBin);
    SaveSpecificMXsBin(CHG_BtoDtoXKL_p, MXsBin);
    SaveSpecificMXsBin(MIX_BtoDtoXKL_p, MXsBin);
    SaveSpecificMXsBin(Signal_BtoDtoXKL_m, MXsBin);
    SaveSpecificMXsBin(CHG_BtoDtoXKL_m, MXsBin);
    SaveSpecificMXsBin(MIX_BtoDtoXKL_m, MXsBin);

    // BR(B -> X KL KL)
    SaveSpecificMXsBin(Signal_BRBtoXKLKL_p, MXsBin);
    SaveSpecificMXsBin(CHG_BRBtoXKLKL_p, MXsBin);
    SaveSpecificMXsBin(MIX_BRBtoXKLKL_p, MXsBin);
    SaveSpecificMXsBin(Signal_BRBtoXKLKL_m, MXsBin);
    SaveSpecificMXsBin(CHG_BRBtoXKLKL_m, MXsBin);
    SaveSpecificMXsBin(MIX_BRBtoXKLKL_m, MXsBin);

    // all uncorrelated uncertainties
    SaveSpecificMXsBin(Signal_all_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHG_all_uncorrelated, MXsBin);
    SaveSpecificMXsBin(MIX_all_uncorrelated, MXsBin);
    SaveSpecificMXsBin(UUBAR_all_uncorrelated, MXsBin);
    SaveSpecificMXsBin(DDBAR_all_uncorrelated, MXsBin);
    SaveSpecificMXsBin(SSBAR_all_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHARM_all_uncorrelated, MXsBin);

    // MC statistical uncertainties
    SaveSpecificMXsBin(Signal_MC_stat, MXsBin);
    SaveSpecificMXsBin(CHG_MC_stat, MXsBin);
    SaveSpecificMXsBin(MIX_MC_stat, MXsBin);
    SaveSpecificMXsBin(UUBAR_MC_stat, MXsBin);
    SaveSpecificMXsBin(DDBAR_MC_stat, MXsBin);
    SaveSpecificMXsBin(SSBAR_MC_stat, MXsBin);
    SaveSpecificMXsBin(CHARM_MC_stat, MXsBin);

    // all of uncorrelated uncertainties + MC statistical uncertainties
    SaveSpecificMXsBin(Signal_all_uncorrelated_MC_stat, MXsBin);
    SaveSpecificMXsBin(CHG_all_uncorrelated_MC_stat, MXsBin);
    SaveSpecificMXsBin(MIX_all_uncorrelated_MC_stat, MXsBin);
    SaveSpecificMXsBin(UUBAR_all_uncorrelated_MC_stat, MXsBin);
    SaveSpecificMXsBin(DDBAR_all_uncorrelated_MC_stat, MXsBin);
    SaveSpecificMXsBin(SSBAR_all_uncorrelated_MC_stat, MXsBin);
    SaveSpecificMXsBin(CHARM_all_uncorrelated_MC_stat, MXsBin);
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

    // FEI uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_FEI; i++) {
        Signal_FEI_correlated[i]->Write();
        CHG_FEI_correlated[i]->Write();
        MIX_FEI_correlated[i]->Write();
    }

    // FEI uncertainty (uncorrelated)
    Signal_FEI_uncorrelated->Write();
    CHG_FEI_uncorrelated->Write();
    MIX_FEI_uncorrelated->Write();

    // Kaon PID uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_KID; i++) {
        Signal_KID_correlated[i]->Write();
        CHG_KID_correlated[i]->Write();
        MIX_KID_correlated[i]->Write();
        UUBAR_KID_correlated[i]->Write();
        DDBAR_KID_correlated[i]->Write();
        SSBAR_KID_correlated[i]->Write();
        CHARM_KID_correlated[i]->Write();
    }

    // Kaon PID uncertainty (uncorrelated)
    Signal_KID_uncorrelated->Write();
    CHG_KID_uncorrelated->Write();
    MIX_KID_uncorrelated->Write();
    UUBAR_KID_uncorrelated->Write();
    DDBAR_KID_uncorrelated->Write();
    SSBAR_KID_uncorrelated->Write();
    CHARM_KID_uncorrelated->Write();

    // Pion PID uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_PID; i++) {
        Signal_PID_correlated[i]->Write();
        CHG_PID_correlated[i]->Write();
        MIX_PID_correlated[i]->Write();
        UUBAR_PID_correlated[i]->Write();
        DDBAR_PID_correlated[i]->Write();
        SSBAR_PID_correlated[i]->Write();
        CHARM_PID_correlated[i]->Write();
    }

    // Pion PID uncertainty (uncorrelated)
    Signal_PID_uncorrelated->Write();
    CHG_PID_uncorrelated->Write();
    MIX_PID_uncorrelated->Write();
    UUBAR_PID_uncorrelated->Write();
    DDBAR_PID_uncorrelated->Write();
    SSBAR_PID_uncorrelated->Write();
    CHARM_PID_uncorrelated->Write();

    // BB BR uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_BR; i++) {
        Signal_BR_correlated[i]->Write();
        CHG_BR_correlated[i]->Write();
        MIX_BR_correlated[i]->Write();
    }

    // BB BR uncertainty (uncorrelated)
    Signal_BR_uncorrelated->Write();
    CHG_BR_uncorrelated->Write();
    MIX_BR_uncorrelated->Write();

    // pi0 uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_pi0; i++) {
        Signal_pi0_correlated[i]->Write();
        CHG_pi0_correlated[i]->Write();
        MIX_pi0_correlated[i]->Write();
        UUBAR_pi0_correlated[i]->Write();
        DDBAR_pi0_correlated[i]->Write();
        SSBAR_pi0_correlated[i]->Write();
        CHARM_pi0_correlated[i]->Write();
    }

    // pi0 uncertainty (uncorrelated)
    Signal_pi0_uncorrelated->Write();
    CHG_pi0_uncorrelated->Write();
    MIX_pi0_uncorrelated->Write();
    UUBAR_pi0_uncorrelated->Write();
    DDBAR_pi0_uncorrelated->Write();
    SSBAR_pi0_uncorrelated->Write();
    CHARM_pi0_uncorrelated->Write();

    // K nu nubar form factor
    Signal_Kff1_p->Write();
    Signal_Kff1_m->Write();
    Signal_Kff2_p->Write();
    Signal_Kff2_m->Write();
    Signal_Kff3_p->Write();
    Signal_Kff3_m->Write();
    CHG_Kff1_p->Write();
    CHG_Kff1_m->Write();
    CHG_Kff2_p->Write();
    CHG_Kff2_m->Write();
    CHG_Kff3_p->Write();
    CHG_Kff3_m->Write();
    MIX_Kff1_p->Write();
    MIX_Kff1_m->Write();
    MIX_Kff2_p->Write();
    MIX_Kff2_m->Write();
    MIX_Kff3_p->Write();
    MIX_Kff3_m->Write();

    // Kstar nu nubar form factor
    CHG_Kstarff1_p->Write();
    CHG_Kstarff1_m->Write();
    CHG_Kstarff2_p->Write();
    CHG_Kstarff2_m->Write();
    CHG_Kstarff3_p->Write();
    CHG_Kstarff3_m->Write();
    CHG_Kstarff4_p->Write();
    CHG_Kstarff4_m->Write();
    CHG_Kstarff5_p->Write();
    CHG_Kstarff5_m->Write();
    CHG_Kstarff6_p->Write();
    CHG_Kstarff6_m->Write();
    CHG_Kstarff7_p->Write();
    CHG_Kstarff7_m->Write();
    CHG_Kstarff8_p->Write();
    CHG_Kstarff8_m->Write();
    CHG_Kstarff9_p->Write();
    CHG_Kstarff9_m->Write();
    MIX_Kstarff1_p->Write();
    MIX_Kstarff1_m->Write();
    MIX_Kstarff2_p->Write();
    MIX_Kstarff2_m->Write();
    MIX_Kstarff3_p->Write();
    MIX_Kstarff3_m->Write();
    MIX_Kstarff4_p->Write();
    MIX_Kstarff4_m->Write();
    MIX_Kstarff5_p->Write();
    MIX_Kstarff5_m->Write();
    MIX_Kstarff6_p->Write();
    MIX_Kstarff6_m->Write();
    MIX_Kstarff7_p->Write();
    MIX_Kstarff7_m->Write();
    MIX_Kstarff8_p->Write();
    MIX_Kstarff8_m->Write();
    MIX_Kstarff9_p->Write();
    MIX_Kstarff9_m->Write();

    // K nu nubar OLD form factor
    Signal_Kff_OLD_p->Write();
    CHG_Kff_OLD_p->Write();
    MIX_Kff_OLD_p->Write();
    Signal_Kff_OLD_m->Write();
    CHG_Kff_OLD_m->Write();
    MIX_Kff_OLD_m->Write();

    // fraction uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_fraction; i++) {
        Signal_fraction_correlated[i]->Write();
        CHG_fraction_correlated[i]->Write();
        MIX_fraction_correlated[i]->Write();
    }

    // fraction uncertainty (uncorrelated)
    Signal_fraction_uncorrelated->Write();
    CHG_fraction_uncorrelated->Write();
    MIX_fraction_uncorrelated->Write();

    // fragmentation uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_Fragmentation; i++) {
        Signal_Fragmentation_correlated[i]->Write();
        CHG_Fragmentation_correlated[i]->Write();
        MIX_Fragmentation_correlated[i]->Write();
    }

    // fragmentation uncertainty (uncorrelated)
    Signal_Fragmentation_uncorrelated->Write();
    CHG_Fragmentation_uncorrelated->Write();
    MIX_Fragmentation_uncorrelated->Write();

    // pf uncertainty
    CHG_pf_p->Write();
    MIX_pf_p->Write();
    CHG_pf_m->Write();
    MIX_pf_m->Write();

    // mb uncertainty
    CHG_mb_p->Write();
    MIX_mb_p->Write();
    CHG_mb_m->Write();
    MIX_mb_m->Write();

    // transition uncertainty
    CHG_transition_p->Write();
    MIX_transition_p->Write();
    CHG_transition_m->Write();
    MIX_transition_m->Write();

    // fixed mKstar uncertainty
    CHG_mKstar_p->Write();
    MIX_mKstar_p->Write();
    CHG_mKstar_m->Write();
    MIX_mKstar_m->Write();

    // BDTc uncertainty
    //CHG_BDTc_p->Write();
    //MIX_BDTc_p->Write();
    UUBAR_BDTc_p->Write();
    DDBAR_BDTc_p->Write();
    SSBAR_BDTc_p->Write();
    CHARM_BDTc_p->Write();
    //CHG_BDTc_m->Write();
    //MIX_BDTc_m->Write();
    UUBAR_BDTc_m->Write();
    DDBAR_BDTc_m->Write();
    SSBAR_BDTc_m->Write();
    CHARM_BDTc_m->Write();

    // B2Xnn BR uncertainty
    Signal_Xnn_p->Write();
    CHG_Xnn_p->Write();
    MIX_Xnn_p->Write();
    Signal_Xnn_m->Write();
    CHG_Xnn_m->Write();
    MIX_Xnn_m->Write();

    // multiplicity uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_multiplicity; i++) {
        Signal_multiplicity_correlated[i]->Write();
        CHG_multiplicity_correlated[i]->Write();
        MIX_multiplicity_correlated[i]->Write();
        UUBAR_multiplicity_correlated[i]->Write();
        DDBAR_multiplicity_correlated[i]->Write();
        SSBAR_multiplicity_correlated[i]->Write();
        CHARM_multiplicity_correlated[i]->Write();
    }

    // multiplicity uncertainty (uncorrelated)
    Signal_multiplicity_uncorrelated->Write();
    CHG_multiplicity_uncorrelated->Write();
    MIX_multiplicity_uncorrelated->Write();
    UUBAR_multiplicity_uncorrelated->Write();
    DDBAR_multiplicity_uncorrelated->Write();
    SSBAR_multiplicity_uncorrelated->Write();
    CHARM_multiplicity_uncorrelated->Write();

    // B -> [D -> X KL0] anything uncertainties
    Signal_BtoDtoXKL_p->Write();
    CHG_BtoDtoXKL_p->Write();
    MIX_BtoDtoXKL_p->Write();
    Signal_BtoDtoXKL_m->Write();
    CHG_BtoDtoXKL_m->Write();
    MIX_BtoDtoXKL_m->Write();

    // BR(B -> X KL KL)
    Signal_BRBtoXKLKL_p->Write();
    CHG_BRBtoXKLKL_p->Write();
    MIX_BRBtoXKLKL_p->Write();
    Signal_BRBtoXKLKL_m->Write();
    CHG_BRBtoXKLKL_m->Write();
    MIX_BRBtoXKLKL_m->Write();

    // all uncorrelated uncertainties
    Signal_all_uncorrelated->Write();
    CHG_all_uncorrelated->Write(); 
    MIX_all_uncorrelated->Write(); 
    UUBAR_all_uncorrelated->Write(); 
    DDBAR_all_uncorrelated->Write(); 
    SSBAR_all_uncorrelated->Write(); 
    CHARM_all_uncorrelated->Write();

    // MC statistical uncertainties
    Signal_MC_stat->Write();
    CHG_MC_stat->Write();
    MIX_MC_stat->Write();
    UUBAR_MC_stat->Write();
    DDBAR_MC_stat->Write();
    SSBAR_MC_stat->Write();
    CHARM_MC_stat->Write();

    // all of uncorrelated uncertainties + MC statistical uncertainties
    Signal_all_uncorrelated_MC_stat->Write();
    CHG_all_uncorrelated_MC_stat->Write();
    MIX_all_uncorrelated_MC_stat->Write();
    UUBAR_all_uncorrelated_MC_stat->Write();
    DDBAR_all_uncorrelated_MC_stat->Write();
    SSBAR_all_uncorrelated_MC_stat->Write();
    CHARM_all_uncorrelated_MC_stat->Write();

    total_DATA->Write();

    file->Close();
    /* ====================================== */

    return 0;
}
