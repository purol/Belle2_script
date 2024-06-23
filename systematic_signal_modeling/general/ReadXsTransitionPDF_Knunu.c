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

#include "constants.h"
#include "correctors.h"
#include "base.h"
#include "template.h"

# define MCTYPE "MC15ri"

// scale facto for systematic MC sample
# define N_Xsu_nonresonant_nunubar_syst 10000000.0
# define N_Xsd_nonresonant_nunubar_syst 10000000.0

// scale factor for each systematic MC sample (364.436 - 2.763 = 361.673/fb)
# define Scale_Xsu_nonresonant_syst (N_Xsu_nonresonant_nunubar_LS1/N_Xsu_nonresonant_nunubar_syst)
# define Scale_Xsd_nonresonant_syst (N_Xsd_nunubar_LS1/N_Xsd_nonresonant_nunubar_syst)

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
Corrector_phiKL corrector_phiKL;
Corrector_KstarKLKL corrector_KstarKLKL;
Corrector_XsKLKL corrector_XsKLKL;
Corrector_BtoDtoXKL corrector_BtoDtoXKL;
Corrector_Fragmentation corrector_Fragmentation;

TH1D* Xsu_Hmb_weight;
TH1D* Xsu_Lmb_weight;
TH1D* Xsu_Hpf_weight;
TH1D* Xsu_Lpf_weight;
TH1D* Xsd_Hmb_weight;
TH1D* Xsd_Lmb_weight;
TH1D* Xsd_Hpf_weight;
TH1D* Xsd_Lpf_weight;
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
// Rarity module
bool sorting(pair<double, double> a, pair<double, double> b) {

    return a.first < b.first;

}

/* ====================================== */

double GetEvtNum(const char* dirname, const char* included_string, double weight_var = 1.0) {

    double Total_Evt_Num = 0;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

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

double BDTcToWeight(double BDTc) {

    if (BDTc > (5.0 / 6.0)) return std::sqrt(5.0);
    else return std::sqrt(BDTc / (1.0 - BDTc));

}

double GetNominalPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, Corrector_Fragmentation::SystType systtype, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get nominal PDF with appropriate correction
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
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;

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
    double nB02KLphi2KSKL = -1;

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
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_phi__bc", &nB02KLphi2KSKL);

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

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

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
            double Correction_Knn = corrector_Knn.GetCorrectionFactorCancelOutObtainWeight(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, names.at(i), MCTYPE, false);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = 1.0;
            if (nB02KLphi2KSKL < MyEPSILON) Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);
            else Correction_KSKLKL = corrector_phiKL.GetCorrectionFactorAtGeneric(nB02KLphi2KSKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if(CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, systtype, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, systtype, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            Nevt = Nevt + FillTemplate(hist, MVA_var, total_weight, Bsig_M);
        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

void ReadXsTransitionPDF_Knunu()
{

    /* ====================================== */
    // Define PDFs for HistFactory
    // nominal PDFs
    TH1D* Xsu_nominal = new TH1D("Xsu_nominal", "Xsu_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* Xsd_nominal = new TH1D("Xsd_nominal", "Xsd_nominal", RarityBins, BinMIN, BinMAX);

    // transition
    TH1D* Xsu_transition_p = new TH1D("Xsu_transition_p", "Xsu_transition_p", RarityBins, BinMIN, BinMAX);
    TH1D* Xsu_transition_m = new TH1D("Xsu_transition_m", "Xsu_transition_m", RarityBins, BinMIN, BinMAX);
    TH1D* Xsd_transition_p = new TH1D("Xsd_transition_p", "Xsd_transition_p", RarityBins, BinMIN, BinMAX);
    TH1D* Xsd_transition_m = new TH1D("Xsd_transition_m", "Xsd_transition_m", RarityBins, BinMIN, BinMAX);
    /* ====================================== */



    /* ====================================== */
    // define path for Ntuple
    const char* MC_dirname_SIGNAL = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut";

    // for signal modeling
    const char* MC_dirname_syst = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_syst/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut";
    /* ====================================== */



    /* ====================================== */
    // Get PDFs
    // get nominal pdfs
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Xsu_nominal, "Bplus", "SIGNAL", Corrector_Fragmentation::SystType::Nominal, Scale_Xsu_nonresonant_test, "B2Xsnunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Xsd_nominal, "Bzero", "SIGNAL", Corrector_Fragmentation::SystType::Nominal, Scale_Xsd_nonresonant_test, "B02Xsnunu");

    // transition
    GetNominalPDFs(MC_dirname_syst, "B2Xsnunu_Htransition", Xsu_transition_p, "Bplus", "SIGNAL", Corrector_Fragmentation::SystType::Htransition, Scale_Xsu_nonresonant_syst, "B2Xsnunu");
    GetNominalPDFs(MC_dirname_syst, "B02Xsnunu_Htransition", Xsd_transition_p, "Bzero", "SIGNAL", Corrector_Fragmentation::SystType::Htransition, Scale_Xsd_nonresonant_syst, "B02Xsnunu");

    GetNominalPDFs(MC_dirname_syst, "B2Xsnunu_Ltransition", Xsu_transition_m, "Bplus", "SIGNAL", Corrector_Fragmentation::SystType::Ltransition, Scale_Xsu_nonresonant_syst, "B2Xsnunu");
    GetNominalPDFs(MC_dirname_syst, "B02Xsnunu_Ltransition", Xsd_transition_m, "Bzero", "SIGNAL", Corrector_Fragmentation::SystType::Ltransition, Scale_Xsd_nonresonant_syst, "B02Xsnunu");
    /* ====================================== */



    /* ====================================== */
    FILE* fp;

    fp = fopen("Xsu_Htransition_weight.txt", "w");
    fprintf(fp, "%d\n", RarityBins);
    for (int i = 0; i < RarityBins; i++) { // Xsu up
        double correction_factor = 1.0;
        if (Xsu_nominal->GetBinContent(i + 1) > MyEPSILON) correction_factor = Xsu_transition_p->GetBinContent(i + 1) / Xsu_nominal->GetBinContent(i + 1);
        else correction_factor = 1.0;
        fprintf(fp, "%lf\n", correction_factor);
    }
    fclose(fp);

    fp = fopen("Xsu_Ltransition_weight.txt", "w");
    fprintf(fp, "%d\n", RarityBins);
    for (int i = 0; i < RarityBins; i++) { // Xsu down
        double correction_factor = 1.0;
        if (Xsu_nominal->GetBinContent(i + 1) > MyEPSILON) correction_factor = Xsu_transition_m->GetBinContent(i + 1) / Xsu_nominal->GetBinContent(i + 1);
        else correction_factor = 1.0;
        fprintf(fp, "%lf\n", correction_factor);
    }
    fclose(fp);

    fp = fopen("Xsd_Htransition_weight.txt", "w");
    fprintf(fp, "%d\n", RarityBins);
    for (int i = 0; i < RarityBins; i++) { // Xsd up
        double correction_factor = 1.0;
        if (Xsd_nominal->GetBinContent(i + 1) > MyEPSILON)correction_factor = Xsd_transition_p->GetBinContent(i + 1) / Xsd_nominal->GetBinContent(i + 1);
        else correction_factor = 1.0;
        fprintf(fp, "%lf\n", correction_factor);
    }
    fclose(fp);

    fp = fopen("Xsd_Ltransition_weight.txt", "w");
    fprintf(fp, "%d\n", RarityBins);
    for (int i = 0; i < RarityBins; i++) { // Xsd down
        double correction_factor = 1.0;
        if (Xsd_nominal->GetBinContent(i + 1) > MyEPSILON) correction_factor = Xsd_transition_m->GetBinContent(i + 1) / Xsd_nominal->GetBinContent(i + 1);
        else correction_factor = 1.0;
        fprintf(fp, "%lf\n", correction_factor);
    }
    fclose(fp);
    /* ====================================== */

    delete Xsu_nominal;
    delete Xsd_nominal;

    delete Xsu_transition_p;
    delete Xsu_transition_m;
    delete Xsd_transition_p;
    delete Xsd_transition_m;
}
