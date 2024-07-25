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
#include "systematic.h"
using namespace RooFit;
using std::string;
using std::to_string;

# define Knunu_only false

enum PDFtype { // reco level
    nominal = 0,
    trackUP,
    trackDOWN,
    KS0UP,
    KS0DOWN,
    BDTc,
    mbXsuUP,
    mbXsuDOWN,
    mbXsdUP,
    mbXsdDOWN,
    pfXsuUP,
    pfXsuDOWN,
    pfXsdUP,
    pfXsdDOWN,
    TransitionXsuUP,
    TransitionXsuDOWN,
    TransitionXsdUP,
    TransitionXsdDOWN,
    KstardeltaKstarplus,
    KstardeltaKstarzero,
    XnnBRUP,
    XnnBRDOWN,
    BRBtoDtoXKLUP,
    BRBtoDtoXKLDOWN,
    BRXKLKLUP,
    BRXKLKLDOWN,
    EffKLECLUP,
    EffKLECLDOWN,
    ControlFEICAL
};

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

double GetDataPDF(const char* dirname, const char* included_string, TH1D* hist, double weight_var = 1.0) { // get data PDF

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;

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

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            double total_weight = weight_var;

            Nevt = Nevt + FillTemplate(hist, MVA_var, total_weight, Bsig_M);

        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double GetPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, PDFtype pdftype, double weight_var = 1.0, std::string CorrectionType = "otherwise", int true_MXs_region = 0) { // get nominal PDF with appropriate correction
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

    if ((true_MXs_region != 0) && (strcmp(sample, "SIGNAL") != 0)) {
        printf("selecting true MXs region only can be done for signal sample\n");
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
    double Ngamma_v200_KL = -1;

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

    double Btag_isSignal = -1;

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
        if (strcmp(sample, "SIGNAL") == 0) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Btag->SetBranchAddress("Btag_isSignal", &Btag_isSignal);
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
        if (strcmp(sample, "SIGNAL") == 0) {
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
        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_KL0__bc", &Ngamma_v200_KL);

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

        tree_upsilon->SetBranchAddress("MVA_Continuum", &BDTc_var);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if (strcmp(sample, "SIGNAL") == 0) tree_Xs->GetEntry(j);

            // select the specific true MXs region
            if (strcmp(sample, "SIGNAL") == 0) {
                double MC_MXs = -1;
                if (strcmp(type, "Bplus") == 0) MC_MXs = Mxs_Bc_MC;
                else if (strcmp(type, "Bzero") == 0) MC_MXs = Mxs_B0_MC;

                // sanity check
                if ((MC_MXs > 0.0) && (MC_MXs < 6.0)) {}
                else { // mass is NaN. try to find true mass region by file name
                    if ((strcmp(included_string, "B2Knunu") == 0) || (strcmp(included_string, "B02K0nunu") == 0)) MC_MXs = 0.4868;
                    else if ((strcmp(included_string, "B2Kstarnunu") == 0) || (strcmp(included_string, "B02Kstar0nunu") == 0)) MC_MXs = 0.8916;
                    else if ((strcmp(included_string, "B2Xsnunu") == 0) || (strcmp(included_string, "B02Xsnunu") == 0)) MC_MXs = 1.5;
                    else {
                        printf("MC Mass of Xs cannot be found and the file name is not expected\n");
                        exit(1);
                    }
                }

                if (true_MXs_region == 0) {}
                else if (true_MXs_region == 1) {
                    if ((MC_MXs > 0.0) && (MC_MXs < 0.6)) {}
                    else continue;
                }
                else if (true_MXs_region == 2) {
                    if ((MC_MXs > 0.6) && (MC_MXs < 1.0)) {}
                    else continue;
                }
                else if (true_MXs_region == 3) {
                    if (MC_MXs > 1.0) {}
                    else continue;
                }
            }

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
            double Correction_Knn = corrector_Knn.GetCorrectionFactorCancelOutObtainWeight(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, names.at(i), MCTYPE, true);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = 1.0;
            if (nB02KLphi2KSKL < MyEPSILON) Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);
            else Correction_KSKLKL = corrector_phiKL.GetCorrectionFactorAtGeneric(nB02KLphi2KSKL);

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

            // lets care about systematics
            if (pdftype == PDFtype::nominal) {}
            else if (pdftype == PDFtype::trackUP) {
                double Ntrack = GetNtracks(Upsilon_ID, Bsig_ID);
                double track_correction = (1 + (Ntrack * track_rel_uncertainty / 100.0));
                total_weight = total_weight * track_correction;
            }
            else if (pdftype == PDFtype::trackDOWN) {
                double Ntrack = GetNtracks(Upsilon_ID, Bsig_ID);
                double track_correction = (1 - (Ntrack * track_rel_uncertainty / 100.0));
                total_weight = total_weight * track_correction;
            }
            else if (pdftype == PDFtype::KS0UP) {
                double KS0_correction = 1 + (KS0_rel_uncertainty * temp_KS0_3D_distance / 100.0);
                total_weight = total_weight * KS0_correction;
            }
            else if (pdftype == PDFtype::KS0DOWN) {
                double KS0_correction = 1 - (KS0_rel_uncertainty * temp_KS0_3D_distance / 100.0);
                total_weight = total_weight * KS0_correction;
            }
            else if (pdftype == PDFtype::BDTc) {
                double BDTc_weight = BDTcToWeight(BDTc_var);
                total_weight = total_weight * BDTc_weight;
            }
            else if (pdftype == PDFtype::mbXsuUP) {
                double mb_correction = ReadWeightHist(Xsu_Hmb_weight, invM * invM);
                total_weight = total_weight * mb_correction;
            }
            else if (pdftype == PDFtype::mbXsuDOWN) {
                double mb_correction = ReadWeightHist(Xsu_Lmb_weight, invM * invM);
                total_weight = total_weight * mb_correction;
            }
            else if (pdftype == PDFtype::mbXsdUP) {
                double mb_correction = ReadWeightHist(Xsd_Hmb_weight, invM * invM);
                total_weight = total_weight * mb_correction;
            }
            else if (pdftype == PDFtype::mbXsdDOWN) {
                double mb_correction = ReadWeightHist(Xsd_Lmb_weight, invM * invM);
                total_weight = total_weight * mb_correction;
            }
            else if (pdftype == PDFtype::pfXsuUP) {
                double pf_correction = ReadWeightHist(Xsu_Hpf_weight, Mxs_Bc_MC);
                total_weight = total_weight * pf_correction;
            }
            else if (pdftype == PDFtype::pfXsuDOWN) {
                double pf_correction = ReadWeightHist(Xsu_Lpf_weight, Mxs_Bc_MC);
                total_weight = total_weight * pf_correction;
            }
            else if (pdftype == PDFtype::pfXsdUP) {
                double pf_correction = ReadWeightHist(Xsd_Hpf_weight, Mxs_B0_MC);
                total_weight = total_weight * pf_correction;
            }
            else if (pdftype == PDFtype::pfXsdDOWN) {
                double pf_correction = ReadWeightHist(Xsd_Lpf_weight, Mxs_B0_MC);
                total_weight = total_weight * pf_correction;
            }
            else if (pdftype == PDFtype::TransitionXsuUP) {
                double bin_index = ReturnBinIndex(MVA_var, Bsig_M);
                double Transition_correction = ReadWeightHist(Xsu_Htransition_weight, bin_index);
                total_weight = total_weight * Transition_correction;
            }
            else if (pdftype == PDFtype::TransitionXsuDOWN) {
                double bin_index = ReturnBinIndex(MVA_var, Bsig_M);
                double Transition_correction = ReadWeightHist(Xsu_Ltransition_weight, bin_index);
                total_weight = total_weight * Transition_correction;
            }
            else if (pdftype == PDFtype::TransitionXsdUP) {
                double bin_index = ReturnBinIndex(MVA_var, Bsig_M);
                double Transition_correction = ReadWeightHist(Xsd_Htransition_weight, bin_index);
                total_weight = total_weight * Transition_correction;
            }
            else if (pdftype == PDFtype::TransitionXsdDOWN) {
                double bin_index = ReturnBinIndex(MVA_var, Bsig_M);
                double Transition_correction = ReadWeightHist(Xsd_Ltransition_weight, bin_index);
                total_weight = total_weight * Transition_correction;
            }
            else if (pdftype == PDFtype::KstardeltaKstarplus) {
                double mKstar_correction = ReadWeightHist(Kstar_delta_weight, invM * invM);
                total_weight = total_weight * mKstar_correction;
            }
            else if (pdftype == PDFtype::KstardeltaKstarzero) {
                double mKstar_correction = ReadWeightHist(K0star_delta_weight, invM * invM);
                total_weight = total_weight * mKstar_correction;
            }
            else if (pdftype == PDFtype::XnnBRUP) {
                double XnnBR_correction = 1.0; // uncertainty for all B->X n nbar, including K n nbar and Kstar n nbar
                if ((N_Xplusnn + N_Xzeronn) > 0.5) XnnBR_correction = 2.0;
                total_weight = total_weight * XnnBR_correction;
            }
            else if (pdftype == PDFtype::XnnBRDOWN) {
                double XnnBR_correction = 1.0; // uncertainty for all B->X n nbar, including K n nbar and Kstar n nbar
                if ((N_Xplusnn + N_Xzeronn) > 0.5) XnnBR_correction = 0.0;
                total_weight = total_weight * XnnBR_correction;
            }
            else if (pdftype == PDFtype::BRBtoDtoXKLUP) {
                if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) {
                    double BRBtoKtoXKL_correction = (1 + corrector_BtoDtoXKL.GetRelativeUncertainty(nDptoXKL + nD0toXKL));
                    total_weight = total_weight * BRBtoKtoXKL_correction;
                }
            }
            else if (pdftype == PDFtype::BRBtoDtoXKLDOWN) {
                if ((strcmp(sample, "CHG") == 0) || (strcmp(sample, "MIX") == 0) || (strcmp(sample, "SIGNAL") == 0)) {
                    double BRBtoKtoXKL_correction = (1 - corrector_BtoDtoXKL.GetRelativeUncertainty(nDptoXKL + nD0toXKL));
                    total_weight = total_weight * BRBtoKtoXKL_correction;
                }
            }
            else if (pdftype == PDFtype::BRXKLKLUP) {
                double XKLKL_correction = 1.0; // uncertainty for all B->X KL KL, including K KL KL and Kstar KL KL
                if (nKL_XKLKL > 1.5) {
                    XKLKL_correction = 2.0;
                }
                total_weight = total_weight * XKLKL_correction;
            }
            else if (pdftype == PDFtype::BRXKLKLDOWN) {
                double XKLKL_correction = 1.0; // uncertainty for all B->X KL KL, including K KL KL and Kstar KL KL
                if (nKL_XKLKL > 1.5) { 
                    XKLKL_correction = 0.0;
                }
                total_weight = total_weight * XKLKL_correction;
            }
            else if (pdftype == PDFtype::EffKLECLUP) {
                double ECLKL_correction = 1.0;
                if (Ngamma_v200_KL > 0.5) ECLKL_correction = (1 + 0.17);
                total_weight = total_weight * ECLKL_correction;
            }
            else if (pdftype == PDFtype::EffKLECLDOWN) {
                double ECLKL_correction = 1.0;
                if (Ngamma_v200_KL > 0.5) ECLKL_correction = (1 - 0.17);
                total_weight = total_weight * ECLKL_correction;
            }
            else if (pdftype == PDFtype::ControlFEICAL) {
                if (strcmp(sample, "SIGNAL") == 0) { // it works only when the sample is signal
                    if (std::abs(Btag_isSignal - 1.0) < MyEPSILON) { // it is applied only when Btag is correct
                        double New_Correction_FEI = 1.0;

                        if (strcmp(type, "Bplus") == 0) {
                            if (std::abs(Btag_ID - 0.0) < MyEPSILON) New_Correction_FEI = 0.880748;
                            else if (std::abs(Btag_ID - 1.0) < MyEPSILON) New_Correction_FEI = 0.880748;
                            else if (std::abs(Btag_ID - 3.0) < MyEPSILON) New_Correction_FEI = 0.880748;
                            else if (std::abs(Btag_ID - 4.0) < MyEPSILON) New_Correction_FEI = 0.880748;
                            else if (std::abs(Btag_ID - 15.0) < MyEPSILON) New_Correction_FEI = 0.586322;
                            else if (std::abs(Btag_ID - 16.0) < MyEPSILON) New_Correction_FEI = 0.586322;
                            else if (std::abs(Btag_ID - 18.0) < MyEPSILON) New_Correction_FEI = 0.586322;
                            else if (std::abs(Btag_ID - 19.0) < MyEPSILON) New_Correction_FEI = 0.586322;
                            else if (std::abs(Btag_ID - 23.0) < MyEPSILON) New_Correction_FEI = 0.230299;
                            else if (std::abs(Btag_ID - 24.0) < MyEPSILON) New_Correction_FEI = 0.230299;
                            else if (std::abs(Btag_ID - 30.0) < MyEPSILON) New_Correction_FEI = 0.230299;
                            else New_Correction_FEI = 0.230299;
                        }
                        else if (strcmp(type, "Bzero") == 0) New_Correction_FEI = 0.968066;

                        // we also need to cancel out the previous FEI correction factor
                        total_weight = total_weight * (New_Correction_FEI / Correction_FEI);
                    }
                }
            }
            else {
                printf("Undefined pdf type!\n");
                exit(1);
            }

            Nevt = Nevt + FillTemplate(hist, MVA_var, total_weight, Bsig_M);

        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double GetFragmentationPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, int TargetMxsBin, int TargetCategory, bool IsItUp, double weight_var = 1.0, std::string CorrectionType = "otherwise", int true_MXs_region = 0) { // get fragmentation PDF with with appropriate correction
    // Be careful! You should combine this function with `GetNominalPDFs` function!
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

    if ((true_MXs_region != 0) && (strcmp(sample, "SIGNAL") != 0)) {
        printf("selecting true MXs region only can be done for signal sample\n");
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
    double Ngamma_v200_KL = -1;

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

    double Btag_isSignal = -1;

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
        if (strcmp(sample, "SIGNAL") == 0) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Btag->SetBranchAddress("Btag_isSignal", &Btag_isSignal);
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
        if (strcmp(sample, "SIGNAL") == 0) {
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
        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_KL0__bc", &Ngamma_v200_KL);

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

        tree_upsilon->SetBranchAddress("MVA_Continuum", &BDTc_var);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if (strcmp(sample, "SIGNAL") == 0) tree_Xs->GetEntry(j);

            // select the specific true MXs region
            if (strcmp(sample, "SIGNAL") == 0) {
                double MC_MXs = -1;
                if (strcmp(type, "Bplus") == 0) MC_MXs = Mxs_Bc_MC;
                else if (strcmp(type, "Bzero") == 0) MC_MXs = Mxs_B0_MC;

                // sanity check
                if ((MC_MXs > 0.0) && (MC_MXs < 6.0)) {}
                else { // mass is NaN. try to find true mass region by file name
                    if ((strcmp(included_string, "B2Knunu") == 0) || (strcmp(included_string, "B02K0nunu") == 0)) MC_MXs = 0.4868;
                    else if ((strcmp(included_string, "B2Kstarnunu") == 0) || (strcmp(included_string, "B02Kstar0nunu") == 0)) MC_MXs = 0.8916;
                    else if ((strcmp(included_string, "B2Xsnunu") == 0) || (strcmp(included_string, "B02Xsnunu") == 0)) MC_MXs = 1.5;
                    else {
                        printf("MC Mass of Xs cannot be found and the file name is not expected\n");
                        exit(1);
                    }
                }

                if (true_MXs_region == 0) {}
                else if (true_MXs_region == 1) {
                    if ((MC_MXs > 0.0) && (MC_MXs < 0.6)) {}
                    else continue;
                }
                else if (true_MXs_region == 2) {
                    if ((MC_MXs > 0.6) && (MC_MXs < 1.0)) {}
                    else continue;
                }
                else if (true_MXs_region == 3) {
                    if (MC_MXs > 1.0) {}
                    else continue;
                }
            }

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
            double Correction_Knn = corrector_Knn.GetCorrectionFactorCancelOutObtainWeight(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, names.at(i), MCTYPE, true);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = 1.0;
            if (nB02KLphi2KSKL < MyEPSILON) Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);
            else Correction_KSKLKL = corrector_phiKL.GetCorrectionFactorAtGeneric(nB02KLphi2KSKL);

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

            // fragmentation fluctuation
            if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.FluctuateCorrection(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, TargetMxsBin, TargetCategory, IsItUp, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.FluctuateCorrection(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, TargetMxsBin, TargetCategory, IsItUp, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else {
                printf("[GetFragmentationPDFs] I think it is not B2Xsnunubar or B02Xsnunubar. Just abort.\n");
                exit(1);
            }

            Nevt = Nevt + FillTemplate(hist, MVA_var, total_weight, Bsig_M);

        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double GetBBBRPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, int TargetdmID, bool IsItUp, double weight_var = 1.0, std::string CorrectionType = "otherwise", int true_MXs_region = 0) { // get BBBR PDF with appropriate correction
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

    if ((true_MXs_region != 0) && (strcmp(sample, "SIGNAL") != 0)) {
        printf("selecting true MXs region only can be done for signal sample\n");
        exit(1);
    }

    float MVA_var = 0;
    float BDTc_var = 0;

    int upsilon_experiment = 0;
    int upsilon_run = 0;
    unsigned int upsilon_event = 0;
    int upsilon_candidate = 0;
    int upsilon_ncandidates = 0;

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
    double Ngamma_v200_KL = -1;

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

    double Btag_isSignal = -1;

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
        if (strcmp(sample, "SIGNAL") == 0) tree_Xs = (TTree*)input_file->Get("Xs");
        else tree_Xs = nullptr;

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("__experiment__", &upsilon_experiment);
        tree_upsilon->SetBranchAddress("__run__", &upsilon_run);
        tree_upsilon->SetBranchAddress("__event__", &upsilon_event);
        tree_upsilon->SetBranchAddress("__candidate__", &upsilon_candidate);
        tree_upsilon->SetBranchAddress("__ncandidates__", &upsilon_ncandidates);

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Btag->SetBranchAddress("Btag_isSignal", &Btag_isSignal);
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
        if (strcmp(sample, "SIGNAL") == 0) {
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
        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_KL0__bc", &Ngamma_v200_KL);

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

        tree_upsilon->SetBranchAddress("MVA_Continuum", &BDTc_var);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if (strcmp(sample, "SIGNAL") == 0) tree_Xs->GetEntry(j);

            // select the specific true MXs region
            if (strcmp(sample, "SIGNAL") == 0) {
                double MC_MXs = -1;
                if (strcmp(type, "Bplus") == 0) MC_MXs = Mxs_Bc_MC;
                else if (strcmp(type, "Bzero") == 0) MC_MXs = Mxs_B0_MC;

                // sanity check
                if ((MC_MXs > 0.0) && (MC_MXs < 6.0)) {}
                else { // mass is NaN. try to find true mass region by file name
                    if ((strcmp(included_string, "B2Knunu") == 0) || (strcmp(included_string, "B02K0nunu") == 0)) MC_MXs = 0.4868;
                    else if ((strcmp(included_string, "B2Kstarnunu") == 0) || (strcmp(included_string, "B02Kstar0nunu") == 0)) MC_MXs = 0.8916;
                    else if ((strcmp(included_string, "B2Xsnunu") == 0) || (strcmp(included_string, "B02Xsnunu") == 0)) MC_MXs = 1.5;
                    else {
                        printf("MC Mass of Xs cannot be found and the file name is not expected\n");
                        exit(1);
                    }
                }

                if (true_MXs_region == 0) {}
                else if (true_MXs_region == 1) {
                    if ((MC_MXs > 0.0) && (MC_MXs < 0.6)) {}
                    else continue;
                }
                else if (true_MXs_region == 2) {
                    if ((MC_MXs > 0.6) && (MC_MXs < 1.0)) {}
                    else continue;
                }
                else if (true_MXs_region == 3) {
                    if (MC_MXs > 1.0) {}
                    else continue;
                }
            }

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
            double Correction_Knn = corrector_Knn.GetCorrectionFactorCancelOutObtainWeight(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, names.at(i), MCTYPE, true);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = 1.0;
            if (nB02KLphi2KSKL < MyEPSILON) Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);
            else Correction_KSKLKL = corrector_phiKL.GetCorrectionFactorAtGeneric(nB02KLphi2KSKL);

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

            // BBBR fluctuation
            total_weight = total_weight * GetBRFluctuation(upsilon_experiment, upsilon_run, upsilon_event, upsilon_candidate, upsilon_ncandidates, names.at(i), TargetdmID, IsItUp);

            Nevt = Nevt + FillTemplate(hist, MVA_var, total_weight, Bsig_M);

        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

int GetFEIcorrelatedPDFs(const char* dirname, TH1D* CHG_nominal_hist, TH1D* MIX_nominal_hist, TH1D* Signal_nominal_hist_MXs1, TH1D* Signal_nominal_hist_MXs2, TH1D* Signal_nominal_hist_MXs3, TH1D*** CHG_hists, TH1D*** MIX_hists, TH1D*** Signal_hists_MXs1, TH1D*** Signal_hists_MXs2, TH1D*** Signal_hists_MXs3) { // get shape sys histogram from txt file
    int Nentry = 0; // number of eigen values/vectors
    double eigen_value = 0; // eigen value
    double weight_sys[RarityBins * 5] = { 0.0 }; // eigen vector

    FILE* fp;
    fp = fopen(dirname, "r");
    while (true) {
        if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
        for (int i = 0; i < RarityBins * 5; i++) {
            if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
        }
        Nentry++;
    }
    fclose(fp);

    *CHG_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *MIX_hists = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *Signal_hists_MXs1 = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *Signal_hists_MXs2 = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *Signal_hists_MXs3 = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);

    for (int i = 0; i < Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_FEI_correlated" + std::to_string(i) + "_p").c_str(), ("CHG_FEI_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_FEI_correlated" + std::to_string(i) + "_p").c_str(), ("MIX_FEI_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists_MXs1)[i] = new TH1D(("Signal_MXs1_FEI_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_MXs1_FEI_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists_MXs2)[i] = new TH1D(("Signal_MXs2_FEI_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_MXs2_FEI_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists_MXs3)[i] = new TH1D(("Signal_MXs3_FEI_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_MXs3_FEI_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
    }

    for (int i = Nentry; i < 2 * Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_FEI_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHG_FEI_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_FEI_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("MIX_FEI_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists_MXs1)[i] = new TH1D(("Signal_MXs1_FEI_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_MXs1_FEI_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists_MXs2)[i] = new TH1D(("Signal_MXs2_FEI_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_MXs2_FEI_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists_MXs3)[i] = new TH1D(("Signal_MXs3_FEI_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_MXs3_FEI_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
    }

    fp = fopen(dirname, "r");
    for (int i = 0; i < Nentry; i++) {
        fscanf(fp, "%lf\n", &eigen_value);
        for (int j = 0; j < RarityBins * 5; j++) fscanf(fp, "%lf\n", &weight_sys[j]);

        for (int k = 0; k < RarityBins; k++) {
            (*CHG_hists)[i]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 1 * RarityBins]));
            (*Signal_hists_MXs1)[i]->SetBinContent(k + 1, Signal_nominal_hist_MXs1->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 2 * RarityBins]));
            (*Signal_hists_MXs2)[i]->SetBinContent(k + 1, Signal_nominal_hist_MXs2->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 3 * RarityBins]));
            (*Signal_hists_MXs3)[i]->SetBinContent(k + 1, Signal_nominal_hist_MXs3->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 4 * RarityBins]));

            (*CHG_hists)[i + Nentry]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i + Nentry]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 1 * RarityBins]));
            (*Signal_hists_MXs1)[i + Nentry]->SetBinContent(k + 1, Signal_nominal_hist_MXs1->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 2 * RarityBins]));
            (*Signal_hists_MXs2)[i + Nentry]->SetBinContent(k + 1, Signal_nominal_hist_MXs2->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 3 * RarityBins]));
            (*Signal_hists_MXs3)[i + Nentry]->SetBinContent(k + 1, Signal_nominal_hist_MXs3->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 4 * RarityBins]));
        }

    }

    fclose(fp);

    return Nentry;
}

void GetFEIUncorrelatedPDFs(const char* dirname, TH1D* CHG_hist, TH1D* MIX_hist, TH1D* Signal_hist_MXs1, TH1D* Signal_hist_MXs2, TH1D* Signal_hist_MXs3) { // get shape sys histogram from txt file
    FILE* fp;
    fp = fopen(dirname, "r");

    double weight_sys[RarityBins * 5] = { 0.0 };
    for (int i = 0; i < RarityBins * 5; i++) fscanf(fp, "%lf\n", &weight_sys[i]);
    fclose(fp);

    for (int i = 0; i < RarityBins * 5; i++) weight_sys[i] = std::sqrt(weight_sys[i]);

    for (int i = 0; i < RarityBins; i++) CHG_hist->SetBinContent(i + 1, weight_sys[i]);
    for (int i = 0; i < RarityBins; i++) MIX_hist->SetBinContent(i + 1, weight_sys[RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) Signal_hist_MXs1->SetBinContent(i + 1, weight_sys[2 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) Signal_hist_MXs2->SetBinContent(i + 1, weight_sys[3 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) Signal_hist_MXs3->SetBinContent(i + 1, weight_sys[4 * RarityBins + i]);
}

int GetPIDcorrelatedPDFs(const char* dirname, TH1D* CHG_nominal_hist, TH1D* MIX_nominal_hist, TH1D* UUBAR_nominal_hist, TH1D* DDBAR_nominal_hist, TH1D* SSBAR_nominal_hist, TH1D* CHARM_nominal_hist, TH1D* Signal_nominal_hist_MXs1, TH1D* Signal_nominal_hist_MXs2, TH1D* Signal_nominal_hist_MXs3, TH1D*** CHG_hists, TH1D*** MIX_hists, TH1D*** UUBAR_hists, TH1D*** DDBAR_hists, TH1D*** SSBAR_hists, TH1D*** CHARM_hists, TH1D*** Signal_hists_MXs1, TH1D*** Signal_hists_MXs2, TH1D*** Signal_hists_MXs3, bool IsItKID) { // get shape sys histogram from txt file
    int Nentry = 0; // number of eigen values/vectors
    double eigen_value = 0; // eigen value
    double weight_sys[RarityBins * 9] = { 0.0 }; // eigen vector

    FILE* fp;
    fp = fopen(dirname, "r");
    while (true) {
        if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
        for (int i = 0; i < RarityBins * 9; i++) {
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
    *Signal_hists_MXs1 = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *Signal_hists_MXs2 = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *Signal_hists_MXs3 = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);

    for (int i = 0; i < Nentry; i++) {
        if (IsItKID) {
            (*CHG_hists)[i] = new TH1D(("CHG_KID_correlated" + std::to_string(i) + "_p").c_str(), ("CHG_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*MIX_hists)[i] = new TH1D(("MIX_KID_correlated" + std::to_string(i) + "_p").c_str(), ("MIX_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*UUBAR_hists)[i] = new TH1D(("UUBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), ("UUBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*DDBAR_hists)[i] = new TH1D(("DDBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), ("DDBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*SSBAR_hists)[i] = new TH1D(("SSBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), ("SSBAR_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*CHARM_hists)[i] = new TH1D(("CHARM_KID_correlated" + std::to_string(i) + "_p").c_str(), ("CHARM_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*Signal_hists_MXs1)[i] = new TH1D(("Signal_MXs1_KID_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_MXs1_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*Signal_hists_MXs2)[i] = new TH1D(("Signal_MXs2_KID_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_MXs2_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*Signal_hists_MXs3)[i] = new TH1D(("Signal_MXs3_KID_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_MXs3_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        }
        else {
            (*CHG_hists)[i] = new TH1D(("CHG_PID_correlated" + std::to_string(i) + "_p").c_str(), ("CHG_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*MIX_hists)[i] = new TH1D(("MIX_PID_correlated" + std::to_string(i) + "_p").c_str(), ("MIX_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*UUBAR_hists)[i] = new TH1D(("UUBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), ("UUBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*DDBAR_hists)[i] = new TH1D(("DDBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), ("DDBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*SSBAR_hists)[i] = new TH1D(("SSBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), ("SSBAR_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*CHARM_hists)[i] = new TH1D(("CHARM_PID_correlated" + std::to_string(i) + "_p").c_str(), ("CHARM_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*Signal_hists_MXs1)[i] = new TH1D(("Signal_MXs1_PID_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_MXs1_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*Signal_hists_MXs2)[i] = new TH1D(("Signal_MXs2_PID_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_MXs2_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
            (*Signal_hists_MXs3)[i] = new TH1D(("Signal_MXs3_PID_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_MXs3_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
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
            (*Signal_hists_MXs1)[i] = new TH1D(("Signal_MXs1_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_MXs1_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*Signal_hists_MXs2)[i] = new TH1D(("Signal_MXs2_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_MXs2_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*Signal_hists_MXs3)[i] = new TH1D(("Signal_MXs3_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_MXs3_KID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        }
        else {
            (*CHG_hists)[i] = new TH1D(("CHG_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHG_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*MIX_hists)[i] = new TH1D(("MIX_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("MIX_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*UUBAR_hists)[i] = new TH1D(("UUBAR_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("UUBAR_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*DDBAR_hists)[i] = new TH1D(("DDBAR_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("DDBAR_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*SSBAR_hists)[i] = new TH1D(("SSBAR_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("SSBAR_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*CHARM_hists)[i] = new TH1D(("CHARM_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHARM_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*Signal_hists_MXs1)[i] = new TH1D(("Signal_MXs1_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_MXs1_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*Signal_hists_MXs2)[i] = new TH1D(("Signal_MXs2_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_MXs2_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
            (*Signal_hists_MXs3)[i] = new TH1D(("Signal_MXs3_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_MXs3_PID_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        }
    }

    fp = fopen(dirname, "r");
    for (int i = 0; i < Nentry; i++) {
        fscanf(fp, "%lf\n", &eigen_value);
        for (int j = 0; j < RarityBins * 9; j++) fscanf(fp, "%lf\n", &weight_sys[j]);

        for (int k = 0; k < RarityBins; k++) {
            (*CHG_hists)[i]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 1 * RarityBins]));
            (*UUBAR_hists)[i]->SetBinContent(k + 1, UUBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 2 * RarityBins]));
            (*DDBAR_hists)[i]->SetBinContent(k + 1, DDBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 3 * RarityBins]));
            (*SSBAR_hists)[i]->SetBinContent(k + 1, SSBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 4 * RarityBins]));
            (*CHARM_hists)[i]->SetBinContent(k + 1, CHARM_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 5 * RarityBins]));
            (*Signal_hists_MXs1)[i]->SetBinContent(k + 1, Signal_nominal_hist_MXs1->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 6 * RarityBins]));
            (*Signal_hists_MXs2)[i]->SetBinContent(k + 1, Signal_nominal_hist_MXs2->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 7 * RarityBins]));
            (*Signal_hists_MXs3)[i]->SetBinContent(k + 1, Signal_nominal_hist_MXs3->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 8 * RarityBins]));

            (*CHG_hists)[i + Nentry]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i + Nentry]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 1 * RarityBins]));
            (*UUBAR_hists)[i + Nentry]->SetBinContent(k + 1, UUBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 2 * RarityBins]));
            (*DDBAR_hists)[i + Nentry]->SetBinContent(k + 1, DDBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 3 * RarityBins]));
            (*SSBAR_hists)[i + Nentry]->SetBinContent(k + 1, SSBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 4 * RarityBins]));
            (*CHARM_hists)[i + Nentry]->SetBinContent(k + 1, CHARM_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 5 * RarityBins]));
            (*Signal_hists_MXs1)[i + Nentry]->SetBinContent(k + 1, Signal_nominal_hist_MXs1->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 6 * RarityBins]));
            (*Signal_hists_MXs2)[i + Nentry]->SetBinContent(k + 1, Signal_nominal_hist_MXs2->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 7 * RarityBins]));
            (*Signal_hists_MXs3)[i + Nentry]->SetBinContent(k + 1, Signal_nominal_hist_MXs3->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 8 * RarityBins]));
        }

    }

    fclose(fp);

    return Nentry;
}

void GetPIDUncorrelatedPDFs(const char* dirname, TH1D* CHG_hist, TH1D* MIX_hist, TH1D* UUBAR_hist, TH1D* DDBAR_hist, TH1D* SSBAR_hist, TH1D* CHARM_hist, TH1D* Signal_hist_MXs1, TH1D* Signal_hist_MXs2, TH1D* Signal_hist_MXs3) { // get shape sys histogram from txt file
    FILE* fp;
    fp = fopen(dirname, "r");

    double weight_sys[RarityBins * 9] = { 0.0 };
    for (int i = 0; i < RarityBins * 9; i++) fscanf(fp, "%lf\n", &weight_sys[i]);
    fclose(fp);

    for (int i = 0; i < RarityBins * 9; i++) weight_sys[i] = std::sqrt(weight_sys[i]);

    for (int i = 0; i < RarityBins; i++) CHG_hist->SetBinContent(i + 1, weight_sys[i]);
    for (int i = 0; i < RarityBins; i++) MIX_hist->SetBinContent(i + 1, weight_sys[RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) UUBAR_hist->SetBinContent(i + 1, weight_sys[2 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) DDBAR_hist->SetBinContent(i + 1, weight_sys[3 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) SSBAR_hist->SetBinContent(i + 1, weight_sys[4 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) CHARM_hist->SetBinContent(i + 1, weight_sys[5 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) Signal_hist_MXs1->SetBinContent(i + 1, weight_sys[6 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) Signal_hist_MXs2->SetBinContent(i + 1, weight_sys[7 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) Signal_hist_MXs3->SetBinContent(i + 1, weight_sys[8 * RarityBins + i]);
}

int Getpi0correlatedPDFs(const char* dirname, TH1D* CHG_nominal_hist, TH1D* MIX_nominal_hist, TH1D* UUBAR_nominal_hist, TH1D* DDBAR_nominal_hist, TH1D* SSBAR_nominal_hist, TH1D* CHARM_nominal_hist, TH1D* Signal_nominal_hist_MXs1, TH1D* Signal_nominal_hist_MXs2, TH1D* Signal_nominal_hist_MXs3, TH1D*** CHG_hists, TH1D*** MIX_hists, TH1D*** UUBAR_hists, TH1D*** DDBAR_hists, TH1D*** SSBAR_hists, TH1D*** CHARM_hists, TH1D*** Signal_hists_MXs1, TH1D*** Signal_hists_MXs2, TH1D*** Signal_hists_MXs3) { // get shape sys histogram from txt file
    int Nentry = 0; // number of eigen values/vectors
    double eigen_value = 0; // eigen value
    double weight_sys[RarityBins * 9] = { 0.0 }; // eigen vector

    FILE* fp;
    fp = fopen(dirname, "r");
    while (true) {
        if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
        for (int i = 0; i < RarityBins * 9; i++) {
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
    *Signal_hists_MXs1 = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *Signal_hists_MXs2 = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *Signal_hists_MXs3 = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);

    for (int i = 0; i < Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("CHG_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("MIX_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*UUBAR_hists)[i] = new TH1D(("UUBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("UUBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*DDBAR_hists)[i] = new TH1D(("DDBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("DDBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*SSBAR_hists)[i] = new TH1D(("SSBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("SSBAR_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*CHARM_hists)[i] = new TH1D(("CHARM_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("CHARM_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists_MXs1)[i] = new TH1D(("Signal_MXs1_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_MXs1_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists_MXs2)[i] = new TH1D(("Signal_MXs2_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_MXs2_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists_MXs3)[i] = new TH1D(("Signal_MXs3_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_MXs3_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
    }

    for (int i = Nentry; i < 2 * Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHG_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("MIX_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*UUBAR_hists)[i] = new TH1D(("UUBAR_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("UUBAR_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*DDBAR_hists)[i] = new TH1D(("DDBAR_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("DDBAR_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*SSBAR_hists)[i] = new TH1D(("SSBAR_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("SSBAR_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*CHARM_hists)[i] = new TH1D(("CHARM_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHARM_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists_MXs1)[i] = new TH1D(("Signal_MXs1_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_MXs1_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists_MXs2)[i] = new TH1D(("Signal_MXs2_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_MXs2_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists_MXs3)[i] = new TH1D(("Signal_MXs3_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_MXs3_pi0_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
    }

    fp = fopen(dirname, "r");
    for (int i = 0; i < Nentry; i++) {
        fscanf(fp, "%lf\n", &eigen_value);
        for (int j = 0; j < RarityBins * 9; j++) fscanf(fp, "%lf\n", &weight_sys[j]);

        for (int k = 0; k < RarityBins; k++) {
            (*CHG_hists)[i]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 1 * RarityBins]));
            (*UUBAR_hists)[i]->SetBinContent(k + 1, UUBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 2 * RarityBins]));
            (*DDBAR_hists)[i]->SetBinContent(k + 1, DDBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 3 * RarityBins]));
            (*SSBAR_hists)[i]->SetBinContent(k + 1, SSBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 4 * RarityBins]));
            (*CHARM_hists)[i]->SetBinContent(k + 1, CHARM_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 5 * RarityBins]));
            (*Signal_hists_MXs1)[i]->SetBinContent(k + 1, Signal_nominal_hist_MXs1->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 6 * RarityBins]));
            (*Signal_hists_MXs2)[i]->SetBinContent(k + 1, Signal_nominal_hist_MXs2->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 7 * RarityBins]));
            (*Signal_hists_MXs3)[i]->SetBinContent(k + 1, Signal_nominal_hist_MXs3->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 8 * RarityBins]));

            (*CHG_hists)[i + Nentry]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i + Nentry]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 1 * RarityBins]));
            (*UUBAR_hists)[i + Nentry]->SetBinContent(k + 1, UUBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 2 * RarityBins]));
            (*DDBAR_hists)[i + Nentry]->SetBinContent(k + 1, DDBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 3 * RarityBins]));
            (*SSBAR_hists)[i + Nentry]->SetBinContent(k + 1, SSBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 4 * RarityBins]));
            (*CHARM_hists)[i + Nentry]->SetBinContent(k + 1, CHARM_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 5 * RarityBins]));
            (*Signal_hists_MXs1)[i + Nentry]->SetBinContent(k + 1, Signal_nominal_hist_MXs1->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 6 * RarityBins]));
            (*Signal_hists_MXs2)[i + Nentry]->SetBinContent(k + 1, Signal_nominal_hist_MXs2->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 7 * RarityBins]));
            (*Signal_hists_MXs3)[i + Nentry]->SetBinContent(k + 1, Signal_nominal_hist_MXs3->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 8 * RarityBins]));
        }

    }

    fclose(fp);

    return Nentry;
}

void Getpi0UncorrelatedPDFs(const char* dirname, TH1D* CHG_hist, TH1D* MIX_hist, TH1D* UUBAR_hist, TH1D* DDBAR_hist, TH1D* SSBAR_hist, TH1D* CHARM_hist, TH1D* Signal_hist_MXs1, TH1D* Signal_hist_MXs2, TH1D* Signal_hist_MXs3) { // get shape sys histogram from txt file
    FILE* fp;
    fp = fopen(dirname, "r");

    double weight_sys[RarityBins * 9] = { 0.0 };
    for (int i = 0; i < RarityBins * 9; i++) fscanf(fp, "%lf\n", &weight_sys[i]);
    fclose(fp);

    for (int i = 0; i < RarityBins * 9; i++) weight_sys[i] = std::sqrt(weight_sys[i]);

    for (int i = 0; i < RarityBins; i++) CHG_hist->SetBinContent(i + 1, weight_sys[i]);
    for (int i = 0; i < RarityBins; i++) MIX_hist->SetBinContent(i + 1, weight_sys[RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) UUBAR_hist->SetBinContent(i + 1, weight_sys[2 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) DDBAR_hist->SetBinContent(i + 1, weight_sys[3 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) SSBAR_hist->SetBinContent(i + 1, weight_sys[4 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) CHARM_hist->SetBinContent(i + 1, weight_sys[5 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) Signal_hist_MXs1->SetBinContent(i + 1, weight_sys[6 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) Signal_hist_MXs2->SetBinContent(i + 1, weight_sys[7 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) Signal_hist_MXs3->SetBinContent(i + 1, weight_sys[8 * RarityBins + i]);
}

void GetKffPDFs(const char* dirname, const char* included_string, TH1D* hist[7], double Correction_factor_BR[7], const char* type, int charge, double weight_var = 1.0, std::string CorrectionType = "otherwise", int true_MXs_region = 0) { // get Kff uncertainty PDF with appropriate correction
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
        {0.061817247836003,   0.977937066768286,  -0.199543782941545},
        {0.998087080844888,  -0.060749128320321,   0.011477040509285}
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
    double nB02KLphi2KSKL = -1;

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
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_phi__bc", &nB02KLphi2KSKL);

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

            // select the specific true MXs region
            // sanity check
            if ((m_k > 0.0) && (m_k < 6.0)) {}
            else { // mass is NaN. try to find true mass region by file name
                if ((strcmp(included_string, "B2Knunu") == 0) || (strcmp(included_string, "B02K0nunu") == 0)) m_k = 0.4868;
                else if ((strcmp(included_string, "B2Kstarnunu") == 0) || (strcmp(included_string, "B02Kstar0nunu") == 0)) m_k = 0.8916;
                else if ((strcmp(included_string, "B2Xsnunu") == 0) || (strcmp(included_string, "B02Xsnunu") == 0)) m_k = 1.5;
                else {
                    printf("MC Mass of Xs cannot be found and the file name is not expected\n");
                    exit(1);
                }
            }

            if (true_MXs_region == 0) {}
            else if (true_MXs_region == 1) {
                if ((m_k > 0.0) && (m_k < 0.6)) {}
                else continue;
            }
            else if (true_MXs_region == 2) {
                if ((m_k > 0.6) && (m_k < 1.0)) {}
                else continue;
            }
            else if (true_MXs_region == 3) {
                if (m_k > 1.0) {}
                else continue;
            }

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
            double Correction_Knn = corrector_Knn.GetCorrectionFactorCancelOutObtainWeight(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, names.at(i), MCTYPE, true);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = 1.0;
            if (nB02KLphi2KSKL < MyEPSILON) Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);
            else Correction_KSKLKL = corrector_phiKL.GetCorrectionFactorAtGeneric(nB02KLphi2KSKL);

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

double GetKstarffPDFs(const char* dirname, const char* included_string, TH1D* hist[19], double Correction_factor_BR[19], const char* type, int charge, double weight_var = 1.0, int true_MXs_region = 0) { // get Kstarff uncertainty PDF with appropriate correction
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
    double nB02KLphi2KSKL = -1;

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
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_phi__bc", &nB02KLphi2KSKL);

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

            // select the specific true MXs region
            // sanity check
            if ((m_k > 0.0) && (m_k < 6.0)) {}
            else { // mass is NaN. try to find true mass region by file name
                if ((strcmp(included_string, "B2Knunu") == 0) || (strcmp(included_string, "B02K0nunu") == 0)) m_k = 0.4868;
                else if ((strcmp(included_string, "B2Kstarnunu") == 0) || (strcmp(included_string, "B02Kstar0nunu") == 0)) m_k = 0.8916;
                else if ((strcmp(included_string, "B2Xsnunu") == 0) || (strcmp(included_string, "B02Xsnunu") == 0)) m_k = 1.5;
                else {
                    printf("MC Mass of Xs cannot be found and the file name is not expected\n");
                    exit(1);
                }
            }

            if (true_MXs_region == 0) {}
            else if (true_MXs_region == 1) {
                if ((m_k > 0.0) && (m_k < 0.6)) {}
                else continue;
            }
            else if (true_MXs_region == 2) {
                if ((m_k > 0.6) && (m_k < 1.0)) {}
                else continue;
            }
            else if (true_MXs_region == 3) {
                if (m_k > 1.0) {}
                else continue;
            }

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
            double Correction_Knn = corrector_Knn.GetCorrectionFactorCancelOutObtainWeight(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, names.at(i), MCTYPE, true);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = 1.0;
            if (nB02KLphi2KSKL < MyEPSILON) Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);
            else Correction_KSKLKL = corrector_phiKL.GetCorrectionFactorAtGeneric(nB02KLphi2KSKL);

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

double GetNevtWithBDTc(const char* dirname, const char* included_string, const char* type, const char* sample, double weight_var = 1.0, std::string CorrectionType = "otherwise", int true_MXs_region = 0) { // get the number of event after BDTc correction
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

    if ((true_MXs_region != 0) && (strcmp(sample, "SIGNAL") != 0)) {
        printf("selecting true MXs region only can be done for signal sample\n");
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
    double nB02KLphi2KSKL = -1;

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
        if (strcmp(sample, "SIGNAL") == 0) tree_Xs = (TTree*)input_file->Get("Xs");
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
        if (strcmp(sample, "SIGNAL") == 0) {
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

        tree_upsilon->SetBranchAddress("MVA_Continuum", &BDTc_var);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if (strcmp(sample, "SIGNAL") == 0) tree_Xs->GetEntry(j);

            // select the specific true MXs region
            if (strcmp(sample, "SIGNAL") == 0) {
                double MC_MXs = -1;
                if (strcmp(type, "Bplus") == 0) MC_MXs = Mxs_Bc_MC;
                else if (strcmp(type, "Bzero") == 0) MC_MXs = Mxs_B0_MC;

                // sanity check
                if ((MC_MXs > 0.0) && (MC_MXs < 6.0)) {}
                else { // mass is NaN. try to find true mass region by file name
                    if ((strcmp(included_string, "B2Knunu") == 0) || (strcmp(included_string, "B02K0nunu") == 0)) MC_MXs = 0.4868;
                    else if ((strcmp(included_string, "B2Kstarnunu") == 0) || (strcmp(included_string, "B02Kstar0nunu") == 0)) MC_MXs = 0.8916;
                    else if ((strcmp(included_string, "B2Xsnunu") == 0) || (strcmp(included_string, "B02Xsnunu") == 0)) MC_MXs = 1.5;
                    else {
                        printf("MC Mass of Xs cannot be found and the file name is not expected\n");
                        exit(1);
                    }
                }

                if (true_MXs_region == 0) {}
                else if (true_MXs_region == 1) {
                    if ((MC_MXs > 0.0) && (MC_MXs < 0.6)) {}
                    else continue;
                }
                else if (true_MXs_region == 2) {
                    if ((MC_MXs > 0.6) && (MC_MXs < 1.0)) {}
                    else continue;
                }
                else if (true_MXs_region == 3) {
                    if (MC_MXs > 1.0) {}
                    else continue;
                }
            }

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
            double Correction_Knn = corrector_Knn.GetCorrectionFactorCancelOutObtainWeight(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, names.at(i), MCTYPE, true);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = 1.0;
            if (nB02KLphi2KSKL < MyEPSILON) Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);
            else Correction_KSKLKL = corrector_phiKL.GetCorrectionFactorAtGeneric(nB02KLphi2KSKL);

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

            Nevt = Nevt + DoesNotFillTemplate(MVA_var, total_weight, Bsig_M);

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

int GetmultiplicitycorrelatedPDFs(const char* dirname, TH1D* CHG_nominal_hist, TH1D* MIX_nominal_hist, TH1D* UUBAR_nominal_hist, TH1D* DDBAR_nominal_hist, TH1D* SSBAR_nominal_hist, TH1D* CHARM_nominal_hist, TH1D* Signal_nominal_hist_MXs1, TH1D* Signal_nominal_hist_MXs2, TH1D* Signal_nominal_hist_MXs3, TH1D*** CHG_hists, TH1D*** MIX_hists, TH1D*** UUBAR_hists, TH1D*** DDBAR_hists, TH1D*** SSBAR_hists, TH1D*** CHARM_hists, TH1D*** Signal_hists_MXs1, TH1D*** Signal_hists_MXs2, TH1D*** Signal_hists_MXs3) { // get shape sys histogram from txt file
    int Nentry = 0; // number of eigen values/vectors
    double eigen_value = 0; // eigen value
    double weight_sys[RarityBins * 9] = { 0.0 }; // eigen vector

    FILE* fp;
    fp = fopen(dirname, "r");
    while (true) {
        if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
        for (int i = 0; i < RarityBins * 9; i++) {
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
    *Signal_hists_MXs1 = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *Signal_hists_MXs2 = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);
    *Signal_hists_MXs3 = (TH1D**)malloc(sizeof(TH1D*) * Nentry * 2);

    for (int i = 0; i < Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), ("CHG_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), ("MIX_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*UUBAR_hists)[i] = new TH1D(("UUBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), ("UUBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*DDBAR_hists)[i] = new TH1D(("DDBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), ("DDBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*SSBAR_hists)[i] = new TH1D(("SSBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), ("SSBAR_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*CHARM_hists)[i] = new TH1D(("CHARM_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), ("CHARM_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists_MXs1)[i] = new TH1D(("Signal_MXs1_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_MXs1_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists_MXs2)[i] = new TH1D(("Signal_MXs2_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_MXs2_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists_MXs3)[i] = new TH1D(("Signal_MXs3_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_MXs3_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
    }

    for (int i = Nentry; i < 2 * Nentry; i++) {
        (*CHG_hists)[i] = new TH1D(("CHG_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHG_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*MIX_hists)[i] = new TH1D(("MIX_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("MIX_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*UUBAR_hists)[i] = new TH1D(("UUBAR_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("UUBAR_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*DDBAR_hists)[i] = new TH1D(("DDBAR_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("DDBAR_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*SSBAR_hists)[i] = new TH1D(("SSBAR_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("SSBAR_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*CHARM_hists)[i] = new TH1D(("CHARM_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("CHARM_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists_MXs1)[i] = new TH1D(("Signal_MXs1_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_MXs1_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists_MXs2)[i] = new TH1D(("Signal_MXs2_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_MXs2_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
        (*Signal_hists_MXs3)[i] = new TH1D(("Signal_MXs3_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), ("Signal_MXs3_multiplicity_correlated" + std::to_string(i - Nentry) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
    }

    fp = fopen(dirname, "r");
    for (int i = 0; i < Nentry; i++) {
        fscanf(fp, "%lf\n", &eigen_value);
        for (int j = 0; j < RarityBins * 9; j++) fscanf(fp, "%lf\n", &weight_sys[j]);

        for (int k = 0; k < RarityBins; k++) {
            (*CHG_hists)[i]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 1 * RarityBins]));
            (*UUBAR_hists)[i]->SetBinContent(k + 1, UUBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 2 * RarityBins]));
            (*DDBAR_hists)[i]->SetBinContent(k + 1, DDBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 3 * RarityBins]));
            (*SSBAR_hists)[i]->SetBinContent(k + 1, SSBAR_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 4 * RarityBins]));
            (*CHARM_hists)[i]->SetBinContent(k + 1, CHARM_nominal_hist->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 5 * RarityBins]));
            (*Signal_hists_MXs1)[i]->SetBinContent(k + 1, Signal_nominal_hist_MXs1->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 6 * RarityBins]));
            (*Signal_hists_MXs2)[i]->SetBinContent(k + 1, Signal_nominal_hist_MXs2->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 7 * RarityBins]));
            (*Signal_hists_MXs3)[i]->SetBinContent(k + 1, Signal_nominal_hist_MXs3->GetBinContent(k + 1) * (1 + eigen_value * weight_sys[k + 8 * RarityBins]));

            (*CHG_hists)[i + Nentry]->SetBinContent(k + 1, CHG_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 0 * RarityBins]));
            (*MIX_hists)[i + Nentry]->SetBinContent(k + 1, MIX_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 1 * RarityBins]));
            (*UUBAR_hists)[i + Nentry]->SetBinContent(k + 1, UUBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 2 * RarityBins]));
            (*DDBAR_hists)[i + Nentry]->SetBinContent(k + 1, DDBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 3 * RarityBins]));
            (*SSBAR_hists)[i + Nentry]->SetBinContent(k + 1, SSBAR_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 4 * RarityBins]));
            (*CHARM_hists)[i + Nentry]->SetBinContent(k + 1, CHARM_nominal_hist->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 5 * RarityBins]));
            (*Signal_hists_MXs1)[i + Nentry]->SetBinContent(k + 1, Signal_nominal_hist_MXs1->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 6 * RarityBins]));
            (*Signal_hists_MXs2)[i + Nentry]->SetBinContent(k + 1, Signal_nominal_hist_MXs2->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 7 * RarityBins]));
            (*Signal_hists_MXs3)[i + Nentry]->SetBinContent(k + 1, Signal_nominal_hist_MXs3->GetBinContent(k + 1) * (1 - eigen_value * weight_sys[k + 8 * RarityBins]));
        }

    }

    fclose(fp);

    return Nentry;
}

void GetmultiplicityUncorrelatedPDFs(const char* dirname, TH1D* CHG_hist, TH1D* MIX_hist, TH1D* UUBAR_hist, TH1D* DDBAR_hist, TH1D* SSBAR_hist, TH1D* CHARM_hist, TH1D* Signal_hist_MXs1, TH1D* Signal_hist_MXs2, TH1D* Signal_hist_MXs3) { // get shape sys histogram from txt file
    FILE* fp;
    fp = fopen(dirname, "r");

    double weight_sys[RarityBins * 9] = { 0.0 };
    for (int i = 0; i < RarityBins * 9; i++) fscanf(fp, "%lf\n", &weight_sys[i]);
    fclose(fp);

    for (int i = 0; i < RarityBins * 9; i++) weight_sys[i] = std::sqrt(weight_sys[i]);

    for (int i = 0; i < RarityBins; i++) CHG_hist->SetBinContent(i + 1, weight_sys[i]);
    for (int i = 0; i < RarityBins; i++) MIX_hist->SetBinContent(i + 1, weight_sys[RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) UUBAR_hist->SetBinContent(i + 1, weight_sys[2 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) DDBAR_hist->SetBinContent(i + 1, weight_sys[3 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) SSBAR_hist->SetBinContent(i + 1, weight_sys[4 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) CHARM_hist->SetBinContent(i + 1, weight_sys[5 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) Signal_hist_MXs1->SetBinContent(i + 1, weight_sys[6 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) Signal_hist_MXs2->SetBinContent(i + 1, weight_sys[7 * RarityBins + i]);
    for (int i = 0; i < RarityBins; i++) Signal_hist_MXs3->SetBinContent(i + 1, weight_sys[8 * RarityBins + i]);
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

void ClearHist(TH1D* hist) {
    hist->Reset();
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

void GetRelativeError(const char* fname, TH1D* relative_error_hist, int Nbin) {
    int Nentry = -1;

    FILE* fp;
    fp = fopen(fname, "r");
    fscanf(fp, "%d\n", &Nentry);
    if (Nentry != Nbin) {
        printf("[GetRelativeError] Nbin and Nentry are different\n");
        exit(1);
    }
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        weight = std::abs(1.0 - weight);
        relative_error_hist->SetBinContent(i + 1, weight);
    }
    fclose(fp);
}

void CalculateEntireRelativeUncertainty(TH1D* entire_rel_error,
    TH1D* signal_MXs1_nominal, TH1D* signal_MXs1_rel_error,
    TH1D* signal_MXs2_nominal, TH1D* signal_MXs2_rel_error,
    TH1D* signal_MXs3_nominal, TH1D* signal_MXs3_rel_error,
    TH1D* CHG_nominal, TH1D* CHG_rel_error,
    TH1D* MIX_nominal, TH1D* MIX_rel_error,
    TH1D* UUBAR_nominal, TH1D* UUBAR_rel_error,
    TH1D* DDBAR_nominal, TH1D* DDBAR_rel_error,
    TH1D* SSBAR_nominal, TH1D* SSBAR_rel_error,
    TH1D* CHARM_nominal, TH1D* CHARM_rel_error, int Nbin) {
    for (int i = 0; i < Nbin; i++) {
        double signal_MXs1_error = signal_MXs1_nominal->GetBinContent(i + 1) * signal_MXs1_rel_error->GetBinContent(i + 1);
        double signal_MXs2_error = signal_MXs2_nominal->GetBinContent(i + 1) * signal_MXs2_rel_error->GetBinContent(i + 1);
        double signal_MXs3_error = signal_MXs3_nominal->GetBinContent(i + 1) * signal_MXs3_rel_error->GetBinContent(i + 1);
        double CHG_error = CHG_nominal->GetBinContent(i + 1) * CHG_rel_error->GetBinContent(i + 1);
        double MIX_error = MIX_nominal->GetBinContent(i + 1) * MIX_rel_error->GetBinContent(i + 1);
        double UUBAR_error = UUBAR_nominal->GetBinContent(i + 1) * UUBAR_rel_error->GetBinContent(i + 1);
        double DDBAR_error = DDBAR_nominal->GetBinContent(i + 1) * DDBAR_rel_error->GetBinContent(i + 1);
        double SSBAR_error = SSBAR_nominal->GetBinContent(i + 1) * SSBAR_rel_error->GetBinContent(i + 1);
        double CHARM_error = CHARM_nominal->GetBinContent(i + 1) * CHARM_rel_error->GetBinContent(i + 1);

        double signal_MXs1_val = signal_MXs1_nominal->GetBinContent(i + 1);
        double signal_MXs2_val = signal_MXs2_nominal->GetBinContent(i + 1);
        double signal_MXs3_val = signal_MXs3_nominal->GetBinContent(i + 1);
        double CHG_val = CHG_nominal->GetBinContent(i + 1);
        double MIX_val = MIX_nominal->GetBinContent(i + 1);
        double UUBAR_val = UUBAR_nominal->GetBinContent(i + 1);
        double DDBAR_val = DDBAR_nominal->GetBinContent(i + 1);
        double SSBAR_val = SSBAR_nominal->GetBinContent(i + 1);
        double CHARM_val = CHARM_nominal->GetBinContent(i + 1);

        entire_rel_error->SetBinContent(i + 1,
            std::sqrt(signal_MXs1_error * signal_MXs1_error +
                signal_MXs2_error * signal_MXs2_error +
                signal_MXs3_error * signal_MXs3_error +
                CHG_error * CHG_error +
                MIX_error * MIX_error +
                UUBAR_error * UUBAR_error +
                DDBAR_error * DDBAR_error +
                SSBAR_error * SSBAR_error +
                CHARM_error * CHARM_error) /
            (signal_MXs1_val + signal_MXs2_val + signal_MXs3_val + CHG_val + MIX_val + UUBAR_val + DDBAR_val + SSBAR_val + CHARM_val)
        );
    }
}

void SaveSpecificMXsBin(TH1D*& hist, int MXsBin) {

    TH1D* replace_hist;
    if (MXsBin == 1) replace_hist = new TH1D("", "", RarityBins_MX1, 0, RarityBins_MX1);
    else if (MXsBin == 2) replace_hist = new TH1D("", "", RarityBins_MX2, 0, RarityBins_MX2);
    else if (MXsBin == 3) replace_hist = new TH1D("", "", RarityBins_MX3, 0, RarityBins_MX3);
    else if (MXsBin == 0) {
        return;
    }
    else {
        printf("Improper Bin index.\n");
        exit(1);
    }

    std::string name_original = std::string(hist->GetName());
    std::string title_original = std::string(hist->GetTitle());

    if (MXsBin == 1) {
        for (int i = 0; i < RarityBins_MX1; i++) {
            replace_hist->SetBinError(i + 1, hist->GetBinError(i + 1));
            replace_hist->SetBinContent(i + 1, hist->GetBinContent(i + 1));
        }
    }
    else if (MXsBin == 2) {
        for (int i = 0; i < RarityBins_MX2; i++) {
            replace_hist->SetBinError(i + 1, hist->GetBinError(i + 1 + RarityBins_MX1));
            replace_hist->SetBinContent(i + 1, hist->GetBinContent(i + 1 + RarityBins_MX1));
        }
    }
    else if (MXsBin == 3) {
        for (int i = 0; i < RarityBins_MX3; i++) {
            replace_hist->SetBinError(i + 1, hist->GetBinError(i + 1 + RarityBins_MX1 + RarityBins_MX2));
            replace_hist->SetBinContent(i + 1, hist->GetBinContent(i + 1 + RarityBins_MX1 + RarityBins_MX2));
        }
    }
    else if (MXsBin == 0) {}
    else {
        printf("Improper Bin index.\n");
        exit(1);
    }

    delete hist;

    replace_hist->SetNameTitle(name_original.c_str(), title_original.c_str());
    hist = replace_hist;
}

void AddPDFs(TH1D* output_hist, TH1D* input_hist) {
    output_hist->Add(output_hist, input_hist);
}

void AddPDFsWithRelativeUncertainty(TH1D* output_hist, TH1D* input_hist, TH1D* output_reluncer_hist, TH1D* input_reluncer_hist, int Nbin) {
    for (int i = 0; i < Nbin; i++) {
        double temp_1 = output_hist->GetBinContent(i + 1);
        double temp_2 = input_hist->GetBinContent(i + 1);
        double temp_3 = output_reluncer_hist->GetBinContent(i + 1);
        double temp_4 = input_reluncer_hist->GetBinContent(i + 1);
        double absolute_uncertainty = std::sqrt((temp_1 * temp_3) * (temp_1 * temp_3) + (temp_2 * temp_4) * (temp_2 * temp_4));
        double relative_uncertainty = absolute_uncertainty / (temp_1 + temp_2);
        output_reluncer_hist->SetBinContent(i + 1, relative_uncertainty);
    }

    output_hist->Add(output_hist, input_hist);
}

int main()
{

    ReadSignalModelingFile();
    ReadEvtFile();


    /* ====================================== */
    // Define PDFs for HistFactory
    // nominal PDFs
    TH1D* Signal_nominal = new TH1D("Signal_nominal", "Signal_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_nominal = new TH1D("Signal_MXs1_nominal", "Signal_MXs1_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_nominal = new TH1D("Signal_MXs2_nominal", "Signal_MXs2_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_nominal = new TH1D("Signal_MXs3_nominal", "Signal_MXs3_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_nominal = new TH1D("CHG_nominal", "CHG_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_nominal = new TH1D("MIX_nominal", "MIX_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_nominal = new TH1D("UUBAR_nominal", "UUBAR_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_nominal = new TH1D("DDBAR_nominal", "DDBAR_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_nominal = new TH1D("SSBAR_nominal", "SSBAR_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_nominal = new TH1D("CHARM_nominal", "CHARM_nominal", RarityBins, BinMIN, BinMAX);

    // track uncertainty
    TH1D* Signal_track_p = new TH1D("Signal_track_p", "Signal_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_track_p = new TH1D("Signal_MXs1_track_p", "Signal_MXs1_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_track_p = new TH1D("Signal_MXs2_track_p", "Signal_MXs2_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_track_p = new TH1D("Signal_MXs3_track_p", "Signal_MXs3_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_track_p = new TH1D("CHG_track_p", "CHG_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_track_p = new TH1D("MIX_track_p", "MIX_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_track_p = new TH1D("UUBAR_track_p", "UUBAR_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_track_p = new TH1D("DDBAR_track_p", "DDBAR_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_track_p = new TH1D("SSBAR_track_p", "SSBAR_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_track_p = new TH1D("CHARM_track_p", "CHARM_track_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_track_m = new TH1D("Signal_track_m", "Signal_track_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_track_m = new TH1D("Signal_MXs1_track_m", "Signal_MXs1_track_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_track_m = new TH1D("Signal_MXs2_track_m", "Signal_MXs2_track_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_track_m = new TH1D("Signal_MXs3_track_m", "Signal_MXs3_track_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_track_m = new TH1D("CHG_track_m", "CHG_track_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_track_m = new TH1D("MIX_track_m", "MIX_track_m", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_track_m = new TH1D("UUBAR_track_m", "UUBAR_track_m", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_track_m = new TH1D("DDBAR_track_m", "DDBAR_track_m", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_track_m = new TH1D("SSBAR_track_m", "SSBAR_track_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_track_m = new TH1D("CHARM_track_m", "CHARM_track_m", RarityBins, BinMIN, BinMAX);

    // KS0 uncertainty
    TH1D* Signal_KS0_p = new TH1D("Signal_KS0_p", "Signal_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_KS0_p = new TH1D("Signal_MXs1_KS0_p", "Signal_MXs1_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_KS0_p = new TH1D("Signal_MXs2_KS0_p", "Signal_MXs2_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_KS0_p = new TH1D("Signal_MXs3_KS0_p", "Signal_MXs3_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_KS0_p = new TH1D("CHG_KS0_p", "CHG_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_KS0_p = new TH1D("MIX_KS0_p", "MIX_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_KS0_p = new TH1D("UUBAR_KS0_p", "UUBAR_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_KS0_p = new TH1D("DDBAR_KS0_p", "DDBAR_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_KS0_p = new TH1D("SSBAR_KS0_p", "SSBAR_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_KS0_p = new TH1D("CHARM_KS0_p", "CHARM_KS0_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_KS0_m = new TH1D("Signal_KS0_m", "Signal_KS0_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_KS0_m = new TH1D("Signal_MXs1_KS0_m", "Signal_MXs1_KS0_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_KS0_m = new TH1D("Signal_MXs2_KS0_m", "Signal_MXs2_KS0_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_KS0_m = new TH1D("Signal_MXs3_KS0_m", "Signal_MXs3_KS0_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_KS0_m = new TH1D("CHG_KS0_m", "CHG_KS0_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_KS0_m = new TH1D("MIX_KS0_m", "MIX_KS0_m", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_KS0_m = new TH1D("UUBAR_KS0_m", "UUBAR_KS0_m", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_KS0_m = new TH1D("DDBAR_KS0_m", "DDBAR_KS0_m", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_KS0_m = new TH1D("SSBAR_KS0_m", "SSBAR_KS0_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_KS0_m = new TH1D("CHARM_KS0_m", "CHARM_KS0_m", RarityBins, BinMIN, BinMAX);

    // FEI uncertainty (correlated)
    TH1D** Signal_FEI_correlated;
    TH1D** Signal_MXs1_FEI_correlated;
    TH1D** Signal_MXs2_FEI_correlated;
    TH1D** Signal_MXs3_FEI_correlated;
    TH1D** CHG_FEI_correlated;
    TH1D** MIX_FEI_correlated;

    // FEI uncertainty (uncorrelated)
    TH1D* Signal_MXs1_FEI_uncorrelated = new TH1D("Signal_MXs1_FEI_uncorrelated", "Signal_MXs1_FEI_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_FEI_uncorrelated = new TH1D("Signal_MXs2_FEI_uncorrelated", "Signal_MXs2_FEI_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_FEI_uncorrelated = new TH1D("Signal_MXs3_FEI_uncorrelated", "Signal_MXs3_FEI_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_FEI_uncorrelated = new TH1D("CHG_FEI_uncorrelated", "CHG_FEI_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_FEI_uncorrelated = new TH1D("MIX_FEI_uncorrelated", "MIX_FEI_uncorrelated", RarityBins, BinMIN, BinMAX);

    // Kaon PID uncertainty (correlated)
    TH1D** Signal_KID_correlated;
    TH1D** Signal_MXs1_KID_correlated;
    TH1D** Signal_MXs2_KID_correlated;
    TH1D** Signal_MXs3_KID_correlated;
    TH1D** CHG_KID_correlated;
    TH1D** MIX_KID_correlated;
    TH1D** UUBAR_KID_correlated;
    TH1D** DDBAR_KID_correlated;
    TH1D** SSBAR_KID_correlated;
    TH1D** CHARM_KID_correlated;

    // Kaon PID uncertainty (uncorrelated)
    TH1D* Signal_MXs1_KID_uncorrelated = new TH1D("Signal_MXs1_KID_uncorrelated", "Signal_MXs1_KID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_KID_uncorrelated = new TH1D("Signal_MXs2_KID_uncorrelated", "Signal_MXs2_KID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_KID_uncorrelated = new TH1D("Signal_MXs3_KID_uncorrelated", "Signal_MXs3_KID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_KID_uncorrelated = new TH1D("CHG_KID_uncorrelated", "CHG_KID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_KID_uncorrelated = new TH1D("MIX_KID_uncorrelated", "MIX_KID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_KID_uncorrelated = new TH1D("UUBAR_KID_uncorrelated", "UUBAR_KID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_KID_uncorrelated = new TH1D("DDBAR_KID_uncorrelated", "DDBAR_KID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_KID_uncorrelated = new TH1D("SSBAR_KID_uncorrelated", "SSBAR_KID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_KID_uncorrelated = new TH1D("CHARM_KID_uncorrelated", "CHARM_KID_uncorrelated", RarityBins, BinMIN, BinMAX);

    // Pion PID uncertainty (correlated)
    TH1D** Signal_PID_correlated;
    TH1D** Signal_MXs1_PID_correlated;
    TH1D** Signal_MXs2_PID_correlated;
    TH1D** Signal_MXs3_PID_correlated;
    TH1D** CHG_PID_correlated;
    TH1D** MIX_PID_correlated;
    TH1D** UUBAR_PID_correlated;
    TH1D** DDBAR_PID_correlated;
    TH1D** SSBAR_PID_correlated;
    TH1D** CHARM_PID_correlated;

    // Pion PID uncertainty (uncorrelated)
    TH1D* Signal_MXs1_PID_uncorrelated = new TH1D("Signal_MXs1_PID_uncorrelated", "Signal_MXs1_PID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_PID_uncorrelated = new TH1D("Signal_MXs2_PID_uncorrelated", "Signal_MXs2_PID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_PID_uncorrelated = new TH1D("Signal_MXs3_PID_uncorrelated", "Signal_MXs3_PID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_PID_uncorrelated = new TH1D("CHG_PID_uncorrelated", "CHG_PID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_PID_uncorrelated = new TH1D("MIX_PID_uncorrelated", "MIX_PID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_PID_uncorrelated = new TH1D("UUBAR_PID_uncorrelated", "UUBAR_PID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_PID_uncorrelated = new TH1D("DDBAR_PID_uncorrelated", "DDBAR_PID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_PID_uncorrelated = new TH1D("SSBAR_PID_uncorrelated", "SSBAR_PID_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_PID_uncorrelated = new TH1D("CHARM_PID_uncorrelated", "CHARM_PID_uncorrelated", RarityBins, BinMIN, BinMAX);

    // BB BR uncertainty
    std::vector<std::string> name_Signal_BR;
    for (int dmIndex = 0; dmIndex < NBRdmID(); dmIndex++) {
        name_Signal_BR.push_back(std::string("Signal_BR_") + std::to_string(dmIndex) + std::string("_p"));
    }
    for (int dmIndex = 0; dmIndex < NBRdmID(); dmIndex++) {
        name_Signal_BR.push_back(std::string("Signal_BR_") + std::to_string(dmIndex) + std::string("_m"));
    }
    std::vector<TH1D*> Signal_BRs;
    for (unsigned int i = 0; i < name_Signal_BR.size(); i++) {
        TH1D* temp = new TH1D(name_Signal_BR.at(i).c_str(), name_Signal_BR.at(i).c_str(), RarityBins, BinMIN, BinMAX);
        Signal_BRs.push_back(temp);
    }

    std::vector<std::string> name_Signal_MXs1_BR;
    for (int dmIndex = 0; dmIndex < NBRdmID(); dmIndex++) {
        name_Signal_MXs1_BR.push_back(std::string("Signal_MXs1_BR_") + std::to_string(dmIndex) + std::string("_p"));
    }
    for (int dmIndex = 0; dmIndex < NBRdmID(); dmIndex++) {
        name_Signal_MXs1_BR.push_back(std::string("Signal_MXs1_BR_") + std::to_string(dmIndex) + std::string("_m"));
    }
    std::vector<TH1D*> Signal_MXs1_BRs;
    for (unsigned int i = 0; i < name_Signal_MXs1_BR.size(); i++) {
        TH1D* temp = new TH1D(name_Signal_MXs1_BR.at(i).c_str(), name_Signal_MXs1_BR.at(i).c_str(), RarityBins, BinMIN, BinMAX);
        Signal_MXs1_BRs.push_back(temp);
    }

    std::vector<std::string> name_Signal_MXs2_BR;
    for (int dmIndex = 0; dmIndex < NBRdmID(); dmIndex++) {
        name_Signal_MXs2_BR.push_back(std::string("Signal_MXs2_BR_") + std::to_string(dmIndex) + std::string("_p"));
    }
    for (int dmIndex = 0; dmIndex < NBRdmID(); dmIndex++) {
        name_Signal_MXs2_BR.push_back(std::string("Signal_MXs2_BR_") + std::to_string(dmIndex) + std::string("_m"));
    }
    std::vector<TH1D*> Signal_MXs2_BRs;
    for (unsigned int i = 0; i < name_Signal_MXs2_BR.size(); i++) {
        TH1D* temp = new TH1D(name_Signal_MXs2_BR.at(i).c_str(), name_Signal_MXs2_BR.at(i).c_str(), RarityBins, BinMIN, BinMAX);
        Signal_MXs2_BRs.push_back(temp);
    }

    std::vector<std::string> name_Signal_MXs3_BR;
    for (int dmIndex = 0; dmIndex < NBRdmID(); dmIndex++) {
        name_Signal_MXs3_BR.push_back(std::string("Signal_MXs3_BR_") + std::to_string(dmIndex) + std::string("_p"));
    }
    for (int dmIndex = 0; dmIndex < NBRdmID(); dmIndex++) {
        name_Signal_MXs3_BR.push_back(std::string("Signal_MXs3_BR_") + std::to_string(dmIndex) + std::string("_m"));
    }
    std::vector<TH1D*> Signal_MXs3_BRs;
    for (unsigned int i = 0; i < name_Signal_MXs3_BR.size(); i++) {
        TH1D* temp = new TH1D(name_Signal_MXs3_BR.at(i).c_str(), name_Signal_MXs3_BR.at(i).c_str(), RarityBins, BinMIN, BinMAX);
        Signal_MXs3_BRs.push_back(temp);
    }

    std::vector<std::string> name_CHG_BR;
    for (int dmIndex = 0; dmIndex < NBRdmID(); dmIndex++) {
        name_CHG_BR.push_back(std::string("CHG_BR_") + std::to_string(dmIndex) + std::string("_p"));
    }
    for (int dmIndex = 0; dmIndex < NBRdmID(); dmIndex++) {
        name_CHG_BR.push_back(std::string("CHG_BR_") + std::to_string(dmIndex) + std::string("_m"));
    }
    std::vector<TH1D*> CHG_BRs;
    for (unsigned int i = 0; i < name_CHG_BR.size(); i++) {
        TH1D* temp = new TH1D(name_CHG_BR.at(i).c_str(), name_CHG_BR.at(i).c_str(), RarityBins, BinMIN, BinMAX);
        CHG_BRs.push_back(temp);
    }

    std::vector<std::string> name_MIX_BR;
    for (int dmIndex = 0; dmIndex < NBRdmID(); dmIndex++) {
        name_MIX_BR.push_back(std::string("MIX_BR_") + std::to_string(dmIndex) + std::string("_p"));
    }
    for (int dmIndex = 0; dmIndex < NBRdmID(); dmIndex++) {
        name_MIX_BR.push_back(std::string("MIX_BR_") + std::to_string(dmIndex) + std::string("_m"));
    }
    std::vector<TH1D*> MIX_BRs;
    for (unsigned int i = 0; i < name_MIX_BR.size(); i++) {
        TH1D* temp = new TH1D(name_MIX_BR.at(i).c_str(), name_MIX_BR.at(i).c_str(), RarityBins, BinMIN, BinMAX);
        MIX_BRs.push_back(temp);
    }

    // pi0 uncertainty (correlated)
    TH1D** Signal_pi0_correlated;
    TH1D** Signal_MXs1_pi0_correlated;
    TH1D** Signal_MXs2_pi0_correlated;
    TH1D** Signal_MXs3_pi0_correlated;
    TH1D** CHG_pi0_correlated;
    TH1D** MIX_pi0_correlated;
    TH1D** UUBAR_pi0_correlated;
    TH1D** DDBAR_pi0_correlated;
    TH1D** SSBAR_pi0_correlated;
    TH1D** CHARM_pi0_correlated;

    // pi0 uncertainty (uncorrelated)
    TH1D* Signal_MXs1_pi0_uncorrelated = new TH1D("Signal_MXs1_pi0_uncorrelated", "Signal_MXs1_pi0_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_pi0_uncorrelated = new TH1D("Signal_MXs2_pi0_uncorrelated", "Signal_MXs2_pi0_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_pi0_uncorrelated = new TH1D("Signal_MXs3_pi0_uncorrelated", "Signal_MXs3_pi0_uncorrelated", RarityBins, BinMIN, BinMAX);
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
    TH1D* Signal_MXs1_Kff1_p = new TH1D("Signal_MXs1_Kff1_p", "Signal_MXs1_Kff1_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_Kff1_m = new TH1D("Signal_MXs1_Kff1_m", "Signal_MXs1_Kff1_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_Kff2_p = new TH1D("Signal_MXs1_Kff2_p", "Signal_MXs1_Kff2_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_Kff2_m = new TH1D("Signal_MXs1_Kff2_m", "Signal_MXs1_Kff2_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_Kff3_p = new TH1D("Signal_MXs1_Kff3_p", "Signal_MXs1_Kff3_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_Kff3_m = new TH1D("Signal_MXs1_Kff3_m", "Signal_MXs1_Kff3_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kff1_p = new TH1D("Signal_MXs3_Kff1_p", "Signal_MXs3_Kff1_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kff1_m = new TH1D("Signal_MXs3_Kff1_m", "Signal_MXs3_Kff1_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kff2_p = new TH1D("Signal_MXs3_Kff2_p", "Signal_MXs3_Kff2_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kff2_m = new TH1D("Signal_MXs3_Kff2_m", "Signal_MXs3_Kff2_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kff3_p = new TH1D("Signal_MXs3_Kff3_p", "Signal_MXs3_Kff3_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kff3_m = new TH1D("Signal_MXs3_Kff3_m", "Signal_MXs3_Kff3_m", RarityBins, BinMIN, BinMAX);

    // Kstar nu nubar form factor
    TH1D* Signal_Kstarff1_p = new TH1D("Signal_Kstarff1_p", "Signal_Kstarff1_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff1_m = new TH1D("Signal_Kstarff1_m", "Signal_Kstarff1_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff2_p = new TH1D("Signal_Kstarff2_p", "Signal_Kstarff2_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff2_m = new TH1D("Signal_Kstarff2_m", "Signal_Kstarff2_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff3_p = new TH1D("Signal_Kstarff3_p", "Signal_Kstarff3_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff3_m = new TH1D("Signal_Kstarff3_m", "Signal_Kstarff3_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff4_p = new TH1D("Signal_Kstarff4_p", "Signal_Kstarff4_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff4_m = new TH1D("Signal_Kstarff4_m", "Signal_Kstarff4_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff5_p = new TH1D("Signal_Kstarff5_p", "Signal_Kstarff5_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff5_m = new TH1D("Signal_Kstarff5_m", "Signal_Kstarff5_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff6_p = new TH1D("Signal_Kstarff6_p", "Signal_Kstarff6_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff6_m = new TH1D("Signal_Kstarff6_m", "Signal_Kstarff6_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff7_p = new TH1D("Signal_Kstarff7_p", "Signal_Kstarff7_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff7_m = new TH1D("Signal_Kstarff7_m", "Signal_Kstarff7_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff8_p = new TH1D("Signal_Kstarff8_p", "Signal_Kstarff8_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff8_m = new TH1D("Signal_Kstarff8_m", "Signal_Kstarff8_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff9_p = new TH1D("Signal_Kstarff9_p", "Signal_Kstarff9_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarff9_m = new TH1D("Signal_Kstarff9_m", "Signal_Kstarff9_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Kstarff1_p = new TH1D("Signal_MXs2_Kstarff1_p", "Signal_MXs2_Kstarff1_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Kstarff1_m = new TH1D("Signal_MXs2_Kstarff1_m", "Signal_MXs2_Kstarff1_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Kstarff2_p = new TH1D("Signal_MXs2_Kstarff2_p", "Signal_MXs2_Kstarff2_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Kstarff2_m = new TH1D("Signal_MXs2_Kstarff2_m", "Signal_MXs2_Kstarff2_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Kstarff3_p = new TH1D("Signal_MXs2_Kstarff3_p", "Signal_MXs2_Kstarff3_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Kstarff3_m = new TH1D("Signal_MXs2_Kstarff3_m", "Signal_MXs2_Kstarff3_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Kstarff4_p = new TH1D("Signal_MXs2_Kstarff4_p", "Signal_MXs2_Kstarff4_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Kstarff4_m = new TH1D("Signal_MXs2_Kstarff4_m", "Signal_MXs2_Kstarff4_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Kstarff5_p = new TH1D("Signal_MXs2_Kstarff5_p", "Signal_MXs2_Kstarff5_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Kstarff5_m = new TH1D("Signal_MXs2_Kstarff5_m", "Signal_MXs2_Kstarff5_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Kstarff6_p = new TH1D("Signal_MXs2_Kstarff6_p", "Signal_MXs2_Kstarff6_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Kstarff6_m = new TH1D("Signal_MXs2_Kstarff6_m", "Signal_MXs2_Kstarff6_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Kstarff7_p = new TH1D("Signal_MXs2_Kstarff7_p", "Signal_MXs2_Kstarff7_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Kstarff7_m = new TH1D("Signal_MXs2_Kstarff7_m", "Signal_MXs2_Kstarff7_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Kstarff8_p = new TH1D("Signal_MXs2_Kstarff8_p", "Signal_MXs2_Kstarff8_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Kstarff8_m = new TH1D("Signal_MXs2_Kstarff8_m", "Signal_MXs2_Kstarff8_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Kstarff9_p = new TH1D("Signal_MXs2_Kstarff9_p", "Signal_MXs2_Kstarff9_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Kstarff9_m = new TH1D("Signal_MXs2_Kstarff9_m", "Signal_MXs2_Kstarff9_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kstarff1_p = new TH1D("Signal_MXs3_Kstarff1_p", "Signal_MXs3_Kstarff1_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kstarff1_m = new TH1D("Signal_MXs3_Kstarff1_m", "Signal_MXs3_Kstarff1_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kstarff2_p = new TH1D("Signal_MXs3_Kstarff2_p", "Signal_MXs3_Kstarff2_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kstarff2_m = new TH1D("Signal_MXs3_Kstarff2_m", "Signal_MXs3_Kstarff2_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kstarff3_p = new TH1D("Signal_MXs3_Kstarff3_p", "Signal_MXs3_Kstarff3_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kstarff3_m = new TH1D("Signal_MXs3_Kstarff3_m", "Signal_MXs3_Kstarff3_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kstarff4_p = new TH1D("Signal_MXs3_Kstarff4_p", "Signal_MXs3_Kstarff4_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kstarff4_m = new TH1D("Signal_MXs3_Kstarff4_m", "Signal_MXs3_Kstarff4_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kstarff5_p = new TH1D("Signal_MXs3_Kstarff5_p", "Signal_MXs3_Kstarff5_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kstarff5_m = new TH1D("Signal_MXs3_Kstarff5_m", "Signal_MXs3_Kstarff5_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kstarff6_p = new TH1D("Signal_MXs3_Kstarff6_p", "Signal_MXs3_Kstarff6_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kstarff6_m = new TH1D("Signal_MXs3_Kstarff6_m", "Signal_MXs3_Kstarff6_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kstarff7_p = new TH1D("Signal_MXs3_Kstarff7_p", "Signal_MXs3_Kstarff7_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kstarff7_m = new TH1D("Signal_MXs3_Kstarff7_m", "Signal_MXs3_Kstarff7_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kstarff8_p = new TH1D("Signal_MXs3_Kstarff8_p", "Signal_MXs3_Kstarff8_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kstarff8_m = new TH1D("Signal_MXs3_Kstarff8_m", "Signal_MXs3_Kstarff8_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kstarff9_p = new TH1D("Signal_MXs3_Kstarff9_p", "Signal_MXs3_Kstarff9_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kstarff9_m = new TH1D("Signal_MXs3_Kstarff9_m", "Signal_MXs3_Kstarff9_m", RarityBins, BinMIN, BinMAX);

    // fraction
    TH1D* Signal_Kfrac_p = new TH1D("Signal_Kfrac_p", "Signal_Kfrac_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kfrac_m = new TH1D("Signal_Kfrac_m", "Signal_Kfrac_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_Kfrac_p = new TH1D("Signal_MXs1_Kfrac_p", "Signal_MXs1_Kfrac_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_Kfrac_m = new TH1D("Signal_MXs1_Kfrac_m", "Signal_MXs1_Kfrac_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kfrac_p = new TH1D("Signal_MXs3_Kfrac_p", "Signal_MXs3_Kfrac_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kfrac_m = new TH1D("Signal_MXs3_Kfrac_m", "Signal_MXs3_Kfrac_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarfrac_p = new TH1D("Signal_Kstarfrac_p", "Signal_Kstarfrac_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Kstarfrac_m = new TH1D("Signal_Kstarfrac_m", "Signal_Kstarfrac_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Kstarfrac_p = new TH1D("Signal_MXs2_Kstarfrac_p", "Signal_MXs2_Kstarfrac_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Kstarfrac_m = new TH1D("Signal_MXs2_Kstarfrac_m", "Signal_MXs2_Kstarfrac_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kstarfrac_p = new TH1D("Signal_MXs3_Kstarfrac_p", "Signal_MXs3_Kstarfrac_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kstarfrac_m = new TH1D("Signal_MXs3_Kstarfrac_m", "Signal_MXs3_Kstarfrac_m", RarityBins, BinMIN, BinMAX);

    // fragmentation uncertainty
    std::vector<std::string> name_Signal_Fragmentation;
    for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
        for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
            name_Signal_Fragmentation.push_back(std::string("Signal_Xs_frag_decay_") + std::to_string(MxsBin) + std::string("_") + std::to_string(Category) + std::string("_p"));
        }
    }
    for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
        for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
            name_Signal_Fragmentation.push_back(std::string("Signal_Xs_frag_decay_") + std::to_string(MxsBin) + std::string("_") + std::to_string(Category) + std::string("_m"));
        }
    }
    std::vector<TH1D*> Signal_Fragmentaions;
    for (unsigned int i = 0; i < name_Signal_Fragmentation.size(); i++) {
        TH1D* temp = new TH1D(name_Signal_Fragmentation.at(i).c_str(), name_Signal_Fragmentation.at(i).c_str(), RarityBins, BinMIN, BinMAX);
        Signal_Fragmentaions.push_back(temp);
    }

    std::vector<std::string> name_Signal_MXs1_Fragmentation;
    for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
        for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
            name_Signal_MXs1_Fragmentation.push_back(std::string("Signal_MXs1_Xs_frag_decay_") + std::to_string(MxsBin) + std::string("_") + std::to_string(Category) + std::string("_p"));
        }
    }
    for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
        for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
            name_Signal_MXs1_Fragmentation.push_back(std::string("Signal_MXs1_Xs_frag_decay_") + std::to_string(MxsBin) + std::string("_") + std::to_string(Category) + std::string("_m"));
        }
    }
    std::vector<TH1D*> Signal_MXs1_Fragmentaions;
    for (unsigned int i = 0; i < name_Signal_MXs1_Fragmentation.size(); i++) {
        TH1D* temp = new TH1D(name_Signal_MXs1_Fragmentation.at(i).c_str(), name_Signal_MXs1_Fragmentation.at(i).c_str(), RarityBins, BinMIN, BinMAX);
        Signal_MXs1_Fragmentaions.push_back(temp);
    }

    std::vector<std::string> name_Signal_MXs2_Fragmentation;
    for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
        for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
            name_Signal_MXs2_Fragmentation.push_back(std::string("Signal_MXs2_Xs_frag_decay_") + std::to_string(MxsBin) + std::string("_") + std::to_string(Category) + std::string("_p"));
        }
    }
    for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
        for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
            name_Signal_MXs2_Fragmentation.push_back(std::string("Signal_MXs2_Xs_frag_decay_") + std::to_string(MxsBin) + std::string("_") + std::to_string(Category) + std::string("_m"));
        }
    }
    std::vector<TH1D*> Signal_MXs2_Fragmentaions;
    for (unsigned int i = 0; i < name_Signal_MXs2_Fragmentation.size(); i++) {
        TH1D* temp = new TH1D(name_Signal_MXs2_Fragmentation.at(i).c_str(), name_Signal_MXs2_Fragmentation.at(i).c_str(), RarityBins, BinMIN, BinMAX);
        Signal_MXs2_Fragmentaions.push_back(temp);
    }

    std::vector<std::string> name_Signal_MXs3_Fragmentation;
    for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
        for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
            name_Signal_MXs3_Fragmentation.push_back(std::string("Signal_MXs3_Xs_frag_decay_") + std::to_string(MxsBin) + std::string("_") + std::to_string(Category) + std::string("_p"));
        }
    }
    for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
        for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
            name_Signal_MXs3_Fragmentation.push_back(std::string("Signal_MXs3_Xs_frag_decay_") + std::to_string(MxsBin) + std::string("_") + std::to_string(Category) + std::string("_m"));
        }
    }
    std::vector<TH1D*> Signal_MXs3_Fragmentaions;
    for (unsigned int i = 0; i < name_Signal_MXs3_Fragmentation.size(); i++) {
        TH1D* temp = new TH1D(name_Signal_MXs3_Fragmentation.at(i).c_str(), name_Signal_MXs3_Fragmentation.at(i).c_str(), RarityBins, BinMIN, BinMAX);
        Signal_MXs3_Fragmentaions.push_back(temp);
    }

    std::vector<std::string> name_CHG_Fragmentation;
    for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
        for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
            name_CHG_Fragmentation.push_back(std::string("CHG_Xs_frag_decay_") + std::to_string(MxsBin) + std::string("_") + std::to_string(Category) + std::string("_p"));
        }
    }
    for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
        for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
            name_CHG_Fragmentation.push_back(std::string("CHG_Xs_frag_decay_") + std::to_string(MxsBin) + std::string("_") + std::to_string(Category) + std::string("_m"));
        }
    }
    std::vector<TH1D*> CHG_Fragmentaions;
    for (unsigned int i = 0; i < name_CHG_Fragmentation.size(); i++) {
        TH1D* temp = new TH1D(name_CHG_Fragmentation.at(i).c_str(), name_CHG_Fragmentation.at(i).c_str(), RarityBins, BinMIN, BinMAX);
        CHG_Fragmentaions.push_back(temp);
    }

    std::vector<std::string> name_MIX_Fragmentation;
    for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
        for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
            name_MIX_Fragmentation.push_back(std::string("MIX_Xs_frag_decay_") + std::to_string(MxsBin) + std::string("_") + std::to_string(Category) + std::string("_p"));
        }
    }
    for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
        for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
            name_MIX_Fragmentation.push_back(std::string("MIX_Xs_frag_decay_") + std::to_string(MxsBin) + std::string("_") + std::to_string(Category) + std::string("_m"));
        }
    }
    std::vector<TH1D*> MIX_Fragmentaions;
    for (unsigned int i = 0; i < name_MIX_Fragmentation.size(); i++) {
        TH1D* temp = new TH1D(name_MIX_Fragmentation.at(i).c_str(), name_MIX_Fragmentation.at(i).c_str(), RarityBins, BinMIN, BinMAX);
        MIX_Fragmentaions.push_back(temp);
    }

    // pf
    TH1D* Signal_pf_p = new TH1D("Signal_pf_p", "Signal_pf_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_pf_p = new TH1D("Signal_MXs3_pf_p", "Signal_MXs3_pf_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_pf_m = new TH1D("Signal_pf_m", "Signal_pf_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_pf_m = new TH1D("Signal_MXs3_pf_m", "Signal_MXs3_pf_m", RarityBins, BinMIN, BinMAX);

    // mb
    TH1D* Signal_mb_p = new TH1D("Signal_mb_p", "Signal_mb_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_mb_p = new TH1D("Signal_MXs3_mb_p", "Signal_MXs3_mb_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_mb_m = new TH1D("Signal_mb_m", "Signal_mb_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_mb_m = new TH1D("Signal_MXs3_mb_m", "Signal_MXs3_mb_m", RarityBins, BinMIN, BinMAX);

    // transition
    TH1D* Signal_transition_p = new TH1D("Signal_transition_p", "Signal_transition_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_transition_p = new TH1D("Signal_MXs3_transition_p", "Signal_MXs3_transition_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_transition_m = new TH1D("Signal_transition_m", "Signal_transition_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_transition_m = new TH1D("Signal_MXs3_transition_m", "Signal_MXs3_transition_m", RarityBins, BinMIN, BinMAX);

    // mK*
    TH1D* Signal_mKstar_p = new TH1D("Signal_mKstar_p", "Signal_mKstar_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_mKstar_p = new TH1D("Signal_MXs2_mKstar_p", "Signal_MXs2_mKstar_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_mKstar_p = new TH1D("Signal_MXs3_mKstar_p", "Signal_MXs3_mKstar_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_mKstar_m = new TH1D("Signal_mKstar_m", "Signal_mKstar_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_mKstar_m = new TH1D("Signal_MXs2_mKstar_m", "Signal_MXs2_mKstar_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_mKstar_m = new TH1D("Signal_MXs3_mKstar_m", "Signal_MXs3_mKstar_m", RarityBins, BinMIN, BinMAX);

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
    TH1D* Signal_MXs1_Xnn_p = new TH1D("Signal_MXs1_Xnn_p", "Signal_MXs1_Xnn_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Xnn_p = new TH1D("Signal_MXs2_Xnn_p", "Signal_MXs2_Xnn_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Xnn_p = new TH1D("Signal_MXs3_Xnn_p", "Signal_MXs3_Xnn_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Xnn_p = new TH1D("CHG_Xnn_p", "CHG_Xnn_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Xnn_p = new TH1D("MIX_Xnn_p", "MIX_Xnn_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_Xnn_m = new TH1D("Signal_Xnn_m", "Signal_Xnn_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_Xnn_m = new TH1D("Signal_MXs1_Xnn_m", "Signal_MXs1_Xnn_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Xnn_m = new TH1D("Signal_MXs2_Xnn_m", "Signal_MXs2_Xnn_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Xnn_m = new TH1D("Signal_MXs3_Xnn_m", "Signal_MXs3_Xnn_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_Xnn_m = new TH1D("CHG_Xnn_m", "CHG_Xnn_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_Xnn_m = new TH1D("MIX_Xnn_m", "MIX_Xnn_m", RarityBins, BinMIN, BinMAX);

    // multiplicity uncertainty (correlated)
    TH1D** Signal_multiplicity_correlated;
    TH1D** Signal_MXs1_multiplicity_correlated;
    TH1D** Signal_MXs2_multiplicity_correlated;
    TH1D** Signal_MXs3_multiplicity_correlated;
    TH1D** CHG_multiplicity_correlated;
    TH1D** MIX_multiplicity_correlated;
    TH1D** UUBAR_multiplicity_correlated;
    TH1D** DDBAR_multiplicity_correlated;
    TH1D** SSBAR_multiplicity_correlated;
    TH1D** CHARM_multiplicity_correlated;

    // multiplicity uncertainty (uncorrelated)
    TH1D* Signal_MXs1_multiplicity_uncorrelated = new TH1D("Signal_MXs1_multiplicity_uncorrelated", "Signal_MXs1_multiplicity_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_multiplicity_uncorrelated = new TH1D("Signal_MXs2_multiplicity_uncorrelated", "Signal_MXs2_multiplicity_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_multiplicity_uncorrelated = new TH1D("Signal_MXs3_multiplicity_uncorrelated", "Signal_MXs3_multiplicity_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_multiplicity_uncorrelated = new TH1D("CHG_multiplicity_uncorrelated", "CHG_multiplicity_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_multiplicity_uncorrelated = new TH1D("MIX_multiplicity_uncorrelated", "MIX_multiplicity_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_multiplicity_uncorrelated = new TH1D("UUBAR_multiplicity_uncorrelated", "UUBAR_multiplicity_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_multiplicity_uncorrelated = new TH1D("DDBAR_multiplicity_uncorrelated", "DDBAR_multiplicity_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_multiplicity_uncorrelated = new TH1D("SSBAR_multiplicity_uncorrelated", "SSBAR_multiplicity_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_multiplicity_uncorrelated = new TH1D("CHARM_multiplicity_uncorrelated", "CHARM_multiplicity_uncorrelated", RarityBins, BinMIN, BinMAX);

    // B -> [D -> X KL0] anything
    TH1D* Signal_BtoDtoXKL_p = new TH1D("Signal_BtoDtoXKL_p", "Signal_BtoDtoXKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_BtoDtoXKL_p = new TH1D("Signal_MXs1_BtoDtoXKL_p", "Signal_MXs1_BtoDtoXKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_BtoDtoXKL_p = new TH1D("Signal_MXs2_BtoDtoXKL_p", "Signal_MXs2_BtoDtoXKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_BtoDtoXKL_p = new TH1D("Signal_MXs3_BtoDtoXKL_p", "Signal_MXs3_BtoDtoXKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_BtoDtoXKL_p = new TH1D("CHG_BtoDtoXKL_p", "CHG_BtoDtoXKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_BtoDtoXKL_p = new TH1D("MIX_BtoDtoXKL_p", "MIX_BtoDtoXKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_BtoDtoXKL_m = new TH1D("Signal_BtoDtoXKL_m", "Signal_BtoDtoXKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_BtoDtoXKL_m = new TH1D("Signal_MXs1_BtoDtoXKL_m", "Signal_MXs1_BtoDtoXKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_BtoDtoXKL_m = new TH1D("Signal_MXs2_BtoDtoXKL_m", "Signal_MXs2_BtoDtoXKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_BtoDtoXKL_m = new TH1D("Signal_MXs3_BtoDtoXKL_m", "Signal_MXs3_BtoDtoXKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_BtoDtoXKL_m = new TH1D("CHG_BtoDtoXKL_m", "CHG_BtoDtoXKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_BtoDtoXKL_m = new TH1D("MIX_BtoDtoXKL_m", "MIX_BtoDtoXKL_m", RarityBins, BinMIN, BinMAX);

    // BR(B -> X KL KL)
    TH1D* Signal_BRBtoXKLKL_p = new TH1D("Signal_BRBtoXKLKL_p", "Signal_BRBtoXKLKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_BRBtoXKLKL_p = new TH1D("Signal_MXs1_BRBtoXKLKL_p", "Signal_MXs1_BRBtoXKLKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_BRBtoXKLKL_p = new TH1D("Signal_MXs2_BRBtoXKLKL_p", "Signal_MXs2_BRBtoXKLKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_BRBtoXKLKL_p = new TH1D("Signal_MXs3_BRBtoXKLKL_p", "Signal_MXs3_BRBtoXKLKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_BRBtoXKLKL_p = new TH1D("CHG_BRBtoXKLKL_p", "CHG_BRBtoXKLKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_BRBtoXKLKL_p = new TH1D("MIX_BRBtoXKLKL_p", "MIX_BRBtoXKLKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_BRBtoXKLKL_m = new TH1D("Signal_BRBtoXKLKL_m", "Signal_BRBtoXKLKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_BRBtoXKLKL_m = new TH1D("Signal_MXs1_BRBtoXKLKL_m", "Signal_MXs1_BRBtoXKLKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_BRBtoXKLKL_m = new TH1D("Signal_MXs2_BRBtoXKLKL_m", "Signal_MXs2_BRBtoXKLKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_BRBtoXKLKL_m = new TH1D("Signal_MXs3_BRBtoXKLKL_m", "Signal_MXs3_BRBtoXKLKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_BRBtoXKLKL_m = new TH1D("CHG_BRBtoXKLKL_m", "CHG_BRBtoXKLKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_BRBtoXKLKL_m = new TH1D("MIX_BRBtoXKLKL_m", "MIX_BRBtoXKLKL_m", RarityBins, BinMIN, BinMAX);

    // efficiency for ECL cluster from KL0
    TH1D* Signal_EffECLKL_p = new TH1D("Signal_EffECLKL_p", "Signal_EffECLKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_EffECLKL_p = new TH1D("Signal_MXs1_EffECLKL_p", "Signal_MXs1_EffECLKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_EffECLKL_p = new TH1D("Signal_MXs2_EffECLKL_p", "Signal_MXs2_EffECLKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_EffECLKL_p = new TH1D("Signal_MXs3_EffECLKL_p", "Signal_MXs3_EffECLKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_EffECLKL_p = new TH1D("CHG_EffECLKL_p", "CHG_EffECLKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_EffECLKL_p = new TH1D("MIX_EffECLKL_p", "MIX_EffECLKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_EffECLKL_p = new TH1D("UUBAR_EffECLKL_p", "UUBAR_EffECLKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_EffECLKL_p = new TH1D("DDBAR_EffECLKL_p", "DDBAR_EffECLKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_EffECLKL_p = new TH1D("SSBAR_EffECLKL_p", "SSBAR_EffECLKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_EffECLKL_p = new TH1D("CHARM_EffECLKL_p", "CHARM_EffECLKL_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_EffECLKL_m = new TH1D("Signal_EffECLKL_m", "Signal_EffECLKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_EffECLKL_m = new TH1D("Signal_MXs1_EffECLKL_m", "Signal_MXs1_EffECLKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_EffECLKL_m = new TH1D("Signal_MXs2_EffECLKL_m", "Signal_MXs2_EffECLKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_EffECLKL_m = new TH1D("Signal_MXs3_EffECLKL_m", "Signal_MXs3_EffECLKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_EffECLKL_m = new TH1D("CHG_EffECLKL_m", "CHG_EffECLKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_EffECLKL_m = new TH1D("MIX_EffECLKL_m", "MIX_EffECLKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_EffECLKL_m = new TH1D("UUBAR_EffECLKL_m", "UUBAR_EffECLKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_EffECLKL_m = new TH1D("DDBAR_EffECLKL_m", "DDBAR_EffECLKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_EffECLKL_m = new TH1D("SSBAR_EffECLKL_m", "SSBAR_EffECLKL_m", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_EffECLKL_m = new TH1D("CHARM_EffECLKL_m", "CHARM_EffECLKL_m", RarityBins, BinMIN, BinMAX);

    // New FEI CAL
    TH1D* Signal_NEWFEICAL_p = new TH1D("Signal_NEWFEICAL_p", "Signal_NEWFEICAL_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_NEWFEICAL_p = new TH1D("Signal_MXs1_NEWFEICAL_p", "Signal_MXs1_NEWFEICAL_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_NEWFEICAL_p = new TH1D("Signal_MXs2_NEWFEICAL_p", "Signal_MXs2_NEWFEICAL_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_NEWFEICAL_p = new TH1D("Signal_MXs3_NEWFEICAL_p", "Signal_MXs3_NEWFEICAL_p", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_NEWFEICAL_m = new TH1D("Signal_NEWFEICAL_m", "Signal_NEWFEICAL_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_NEWFEICAL_m = new TH1D("Signal_MXs1_NEWFEICAL_m", "Signal_MXs1_NEWFEICAL_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_NEWFEICAL_m = new TH1D("Signal_MXs2_NEWFEICAL_m", "Signal_MXs2_NEWFEICAL_m", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_NEWFEICAL_m = new TH1D("Signal_MXs3_NEWFEICAL_m", "Signal_MXs3_NEWFEICAL_m", RarityBins, BinMIN, BinMAX);

    // all of uncorrelated uncertainties
    TH1D* Signal_MXs1_all_uncorrelated = new TH1D("Signal_MXs1_all_uncorrelated", "Signal_MXs1_all_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_all_uncorrelated = new TH1D("Signal_MXs2_all_uncorrelated", "Signal_MXs2_all_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_all_uncorrelated = new TH1D("Signal_MXs3_all_uncorrelated", "Signal_MXs3_all_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_all_uncorrelated = new TH1D("CHG_all_uncorrelated", "CHG_all_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_all_uncorrelated = new TH1D("MIX_all_uncorrelated", "MIX_all_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_all_uncorrelated = new TH1D("UUBAR_all_uncorrelated", "UUBAR_all_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_all_uncorrelated = new TH1D("DDBAR_all_uncorrelated", "DDBAR_all_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_all_uncorrelated = new TH1D("SSBAR_all_uncorrelated", "SSBAR_all_uncorrelated", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_all_uncorrelated = new TH1D("CHARM_all_uncorrelated", "CHARM_all_uncorrelated", RarityBins, BinMIN, BinMAX);

    // MC statistical uncertainties
    TH1D* Signal_MXs1_MC_stat = new TH1D("Signal_MXs1_MC_stat", "Signal_MXs1_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_MC_stat = new TH1D("Signal_MXs2_MC_stat", "Signal_MXs2_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_MC_stat = new TH1D("Signal_MXs3_MC_stat", "Signal_MXs3_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_MC_stat = new TH1D("CHG_MC_stat", "CHG_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_MC_stat = new TH1D("MIX_MC_stat", "MIX_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_MC_stat = new TH1D("UUBAR_MC_stat", "UUBAR_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_MC_stat = new TH1D("DDBAR_MC_stat", "DDBAR_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_MC_stat = new TH1D("SSBAR_MC_stat", "SSBAR_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_MC_stat = new TH1D("CHARM_MC_stat", "CHARM_MC_stat", RarityBins, BinMIN, BinMAX);

    // additional relative uncertainty
    TH1D* CHG_rel_uncer = new TH1D("CHG_rel_uncer", "CHG_rel_uncer", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_rel_uncer = new TH1D("MIX_rel_uncer", "MIX_rel_uncer", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_rel_uncer = new TH1D("UUBAR_rel_uncer", "UUBAR_rel_uncer", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_rel_uncer = new TH1D("DDBAR_rel_uncer", "DDBAR_rel_uncer", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_rel_uncer = new TH1D("SSBAR_rel_uncer", "SSBAR_rel_uncer", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_rel_uncer = new TH1D("CHARM_rel_uncer", "CHARM_rel_uncer", RarityBins, BinMIN, BinMAX);

    // all of uncorrelated uncertainties + MC statistical uncertainties
    TH1D* Signal_all_uncorrelated_MC_stat = new TH1D("Signal_all_uncorrelated_MC_stat", "Signal_all_uncorrelated_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_all_uncorrelated_MC_stat = new TH1D("Signal_MXs1_all_uncorrelated_MC_stat", "Signal_MXs1_all_uncorrelated_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_all_uncorrelated_MC_stat = new TH1D("Signal_MXs2_all_uncorrelated_MC_stat", "Signal_MXs2_all_uncorrelated_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_all_uncorrelated_MC_stat = new TH1D("Signal_MXs3_all_uncorrelated_MC_stat", "Signal_MXs3_all_uncorrelated_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_all_uncorrelated_MC_stat = new TH1D("CHG_all_uncorrelated_MC_stat", "CHG_all_uncorrelated_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_all_uncorrelated_MC_stat = new TH1D("MIX_all_uncorrelated_MC_stat", "MIX_all_uncorrelated_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_all_uncorrelated_MC_stat = new TH1D("UUBAR_all_uncorrelated_MC_stat", "UUBAR_all_uncorrelated_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_all_uncorrelated_MC_stat = new TH1D("DDBAR_all_uncorrelated_MC_stat", "DDBAR_all_uncorrelated_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_all_uncorrelated_MC_stat = new TH1D("SSBAR_all_uncorrelated_MC_stat", "SSBAR_all_uncorrelated_MC_stat", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_all_uncorrelated_MC_stat = new TH1D("CHARM_all_uncorrelated_MC_stat", "CHARM_all_uncorrelated_MC_stat", RarityBins, BinMIN, BinMAX);

    // all of uncorrelated uncertainties + MC statistical uncertainties for entire sample
    TH1D* entire_all_uncorrelated_MC_stat = new TH1D("entire_all_uncorrelated_MC_stat", "entire_all_uncorrelated_MC_stat", RarityBins, BinMIN, BinMAX);
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

    // for additional relative uncertainty
    const char* relative_uncertainty_file = "./dataMCratio_sideband.txt";

    /* ====================================== */



    /* ====================================== */
    // Get PDFs
    // get nominal pdfs
    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs1_nominal, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs1_nominal, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs1_nominal, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs1_nominal, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs1_nominal, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs1_nominal, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 1);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs2_nominal, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs2_nominal, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs2_nominal, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs2_nominal, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs2_nominal, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs2_nominal, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 2);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_nominal, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_nominal, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_nominal, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_nominal, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_nominal, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_nominal, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

    double NevtCHGWithoutBDTc = GetPDFs(MC_dirname_CHG, "root", CHG_nominal, "Bplus", "CHG", PDFtype::nominal, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise", 0);
    double NevtMIXWithoutBDTc = GetPDFs(MC_dirname_MIX, "root", MIX_nominal, "Bzero", "MIX", PDFtype::nominal, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise", 0);
    double NevtUUBARWithoutBDTc = GetPDFs(MC_dirname_UUBAR, "root", UUBAR_nominal, "Continuum", "UUBAR", PDFtype::nominal, ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"), "otherwise", 0);
    double NevtDDBARWithoutBDTc = GetPDFs(MC_dirname_DDBAR, "root", DDBAR_nominal, "Continuum", "DDBAR", PDFtype::nominal, ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"), "otherwise", 0);
    double NevtSSBARWithoutBDTc = GetPDFs(MC_dirname_SSBAR, "root", SSBAR_nominal, "Continuum", "SSBAR", PDFtype::nominal, ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"), "otherwise", 0);
    double NevtCHARMWithoutBDTc = GetPDFs(MC_dirname_CHARM, "root", CHARM_nominal, "Continuum", "CHARM", PDFtype::nominal, ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"), "otherwise", 0);

    // get track uncertainty pdfs
    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs1_track_p, "Bplus", "SIGNAL", PDFtype::trackUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs1_track_p, "Bplus", "SIGNAL", PDFtype::trackUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs1_track_p, "Bplus", "SIGNAL", PDFtype::trackUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs1_track_p, "Bzero", "SIGNAL", PDFtype::trackUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs1_track_p, "Bzero", "SIGNAL", PDFtype::trackUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs1_track_p, "Bzero", "SIGNAL", PDFtype::trackUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 1);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs2_track_p, "Bplus", "SIGNAL", PDFtype::trackUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs2_track_p, "Bplus", "SIGNAL", PDFtype::trackUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs2_track_p, "Bplus", "SIGNAL", PDFtype::trackUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs2_track_p, "Bzero", "SIGNAL", PDFtype::trackUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs2_track_p, "Bzero", "SIGNAL", PDFtype::trackUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs2_track_p, "Bzero", "SIGNAL", PDFtype::trackUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 2);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_track_p, "Bplus", "SIGNAL", PDFtype::trackUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_track_p, "Bplus", "SIGNAL", PDFtype::trackUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_track_p, "Bplus", "SIGNAL", PDFtype::trackUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_track_p, "Bzero", "SIGNAL", PDFtype::trackUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_track_p, "Bzero", "SIGNAL", PDFtype::trackUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_track_p, "Bzero", "SIGNAL", PDFtype::trackUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

    GetPDFs(MC_dirname_CHG, "root", CHG_track_p, "Bplus", "CHG", PDFtype::trackUP, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise", 0);
    GetPDFs(MC_dirname_MIX, "root", MIX_track_p, "Bzero", "MIX", PDFtype::trackUP, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise", 0);
    GetPDFs(MC_dirname_UUBAR, "root", UUBAR_track_p, "Continuum", "UUBAR", PDFtype::trackUP, ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"), "otherwise", 0);
    GetPDFs(MC_dirname_DDBAR, "root", DDBAR_track_p, "Continuum", "DDBAR", PDFtype::trackUP, ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"), "otherwise", 0);
    GetPDFs(MC_dirname_SSBAR, "root", SSBAR_track_p, "Continuum", "SSBAR", PDFtype::trackUP, ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"), "otherwise", 0);
    GetPDFs(MC_dirname_CHARM, "root", CHARM_track_p, "Continuum", "CHARM", PDFtype::trackUP, ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"), "otherwise", 0);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs1_track_m, "Bplus", "SIGNAL", PDFtype::trackDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs1_track_m, "Bplus", "SIGNAL", PDFtype::trackDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs1_track_m, "Bplus", "SIGNAL", PDFtype::trackDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs1_track_m, "Bzero", "SIGNAL", PDFtype::trackDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs1_track_m, "Bzero", "SIGNAL", PDFtype::trackDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs1_track_m, "Bzero", "SIGNAL", PDFtype::trackDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 1);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs2_track_m, "Bplus", "SIGNAL", PDFtype::trackDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs2_track_m, "Bplus", "SIGNAL", PDFtype::trackDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs2_track_m, "Bplus", "SIGNAL", PDFtype::trackDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs2_track_m, "Bzero", "SIGNAL", PDFtype::trackDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs2_track_m, "Bzero", "SIGNAL", PDFtype::trackDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs2_track_m, "Bzero", "SIGNAL", PDFtype::trackDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 2);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_track_m, "Bplus", "SIGNAL", PDFtype::trackDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_track_m, "Bplus", "SIGNAL", PDFtype::trackDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_track_m, "Bplus", "SIGNAL", PDFtype::trackDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_track_m, "Bzero", "SIGNAL", PDFtype::trackDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_track_m, "Bzero", "SIGNAL", PDFtype::trackDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_track_m, "Bzero", "SIGNAL", PDFtype::trackDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

    GetPDFs(MC_dirname_CHG, "root", CHG_track_m, "Bplus", "CHG", PDFtype::trackDOWN, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise", 0);
    GetPDFs(MC_dirname_MIX, "root", MIX_track_m, "Bzero", "MIX", PDFtype::trackDOWN, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise", 0);
    GetPDFs(MC_dirname_UUBAR, "root", UUBAR_track_m, "Continuum", "UUBAR", PDFtype::trackDOWN, ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"), "otherwise", 0);
    GetPDFs(MC_dirname_DDBAR, "root", DDBAR_track_m, "Continuum", "DDBAR", PDFtype::trackDOWN, ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"), "otherwise", 0);
    GetPDFs(MC_dirname_SSBAR, "root", SSBAR_track_m, "Continuum", "SSBAR", PDFtype::trackDOWN, ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"), "otherwise", 0);
    GetPDFs(MC_dirname_CHARM, "root", CHARM_track_m, "Continuum", "CHARM", PDFtype::trackDOWN, ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"), "otherwise", 0);

    // get KS0 uncertainty pdfs
    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs1_KS0_p, "Bplus", "SIGNAL", PDFtype::KS0UP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs1_KS0_p, "Bplus", "SIGNAL", PDFtype::KS0UP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs1_KS0_p, "Bplus", "SIGNAL", PDFtype::KS0UP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs1_KS0_p, "Bzero", "SIGNAL", PDFtype::KS0UP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs1_KS0_p, "Bzero", "SIGNAL", PDFtype::KS0UP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs1_KS0_p, "Bzero", "SIGNAL", PDFtype::KS0UP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 1);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs2_KS0_p, "Bplus", "SIGNAL", PDFtype::KS0UP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs2_KS0_p, "Bplus", "SIGNAL", PDFtype::KS0UP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs2_KS0_p, "Bplus", "SIGNAL", PDFtype::KS0UP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs2_KS0_p, "Bzero", "SIGNAL", PDFtype::KS0UP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs2_KS0_p, "Bzero", "SIGNAL", PDFtype::KS0UP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs2_KS0_p, "Bzero", "SIGNAL", PDFtype::KS0UP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 2);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_KS0_p, "Bplus", "SIGNAL", PDFtype::KS0UP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_KS0_p, "Bplus", "SIGNAL", PDFtype::KS0UP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_KS0_p, "Bplus", "SIGNAL", PDFtype::KS0UP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_KS0_p, "Bzero", "SIGNAL", PDFtype::KS0UP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_KS0_p, "Bzero", "SIGNAL", PDFtype::KS0UP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_KS0_p, "Bzero", "SIGNAL", PDFtype::KS0UP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

    GetPDFs(MC_dirname_CHG, "root", CHG_KS0_p, "Bplus", "CHG", PDFtype::KS0UP, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise", 0);
    GetPDFs(MC_dirname_MIX, "root", MIX_KS0_p, "Bzero", "MIX", PDFtype::KS0UP, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise", 0);
    GetPDFs(MC_dirname_UUBAR, "root", UUBAR_KS0_p, "Continuum", "UUBAR", PDFtype::KS0UP, ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"), "otherwise", 0);
    GetPDFs(MC_dirname_DDBAR, "root", DDBAR_KS0_p, "Continuum", "DDBAR", PDFtype::KS0UP, ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"), "otherwise", 0);
    GetPDFs(MC_dirname_SSBAR, "root", SSBAR_KS0_p, "Continuum", "SSBAR", PDFtype::KS0UP, ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"), "otherwise", 0);
    GetPDFs(MC_dirname_CHARM, "root", CHARM_KS0_p, "Continuum", "CHARM", PDFtype::KS0UP, ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"), "otherwise", 0);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs1_KS0_m, "Bplus", "SIGNAL", PDFtype::KS0DOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs1_KS0_m, "Bplus", "SIGNAL", PDFtype::KS0DOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs1_KS0_m, "Bplus", "SIGNAL", PDFtype::KS0DOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs1_KS0_m, "Bzero", "SIGNAL", PDFtype::KS0DOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs1_KS0_m, "Bzero", "SIGNAL", PDFtype::KS0DOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs1_KS0_m, "Bzero", "SIGNAL", PDFtype::KS0DOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 1);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs2_KS0_m, "Bplus", "SIGNAL", PDFtype::KS0DOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs2_KS0_m, "Bplus", "SIGNAL", PDFtype::KS0DOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs2_KS0_m, "Bplus", "SIGNAL", PDFtype::KS0DOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs2_KS0_m, "Bzero", "SIGNAL", PDFtype::KS0DOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs2_KS0_m, "Bzero", "SIGNAL", PDFtype::KS0DOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs2_KS0_m, "Bzero", "SIGNAL", PDFtype::KS0DOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 2);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_KS0_m, "Bplus", "SIGNAL", PDFtype::KS0DOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_KS0_m, "Bplus", "SIGNAL", PDFtype::KS0DOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_KS0_m, "Bplus", "SIGNAL", PDFtype::KS0DOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_KS0_m, "Bzero", "SIGNAL", PDFtype::KS0DOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_KS0_m, "Bzero", "SIGNAL", PDFtype::KS0DOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_KS0_m, "Bzero", "SIGNAL", PDFtype::KS0DOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

    GetPDFs(MC_dirname_CHG, "root", CHG_KS0_m, "Bplus", "CHG", PDFtype::KS0DOWN, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise", 0);
    GetPDFs(MC_dirname_MIX, "root", MIX_KS0_m, "Bzero", "MIX", PDFtype::KS0DOWN, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise", 0);
    GetPDFs(MC_dirname_UUBAR, "root", UUBAR_KS0_m, "Continuum", "UUBAR", PDFtype::KS0DOWN, ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"), "otherwise", 0);
    GetPDFs(MC_dirname_DDBAR, "root", DDBAR_KS0_m, "Continuum", "DDBAR", PDFtype::KS0DOWN, ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"), "otherwise", 0);
    GetPDFs(MC_dirname_SSBAR, "root", SSBAR_KS0_m, "Continuum", "SSBAR", PDFtype::KS0DOWN, ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"), "otherwise", 0);
    GetPDFs(MC_dirname_CHARM, "root", CHARM_KS0_m, "Continuum", "CHARM", PDFtype::KS0DOWN, ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"), "otherwise", 0);

    // get FEI uncertainty pdfs (correlated)
    int NPDFs_FEI = GetFEIcorrelatedPDFs(FEI_correlated_info, CHG_nominal, MIX_nominal, Signal_MXs1_nominal, Signal_MXs2_nominal, Signal_MXs3_nominal, &CHG_FEI_correlated, &MIX_FEI_correlated, &Signal_MXs1_FEI_correlated, &Signal_MXs2_FEI_correlated, &Signal_MXs3_FEI_correlated);

    // get FEI uncertainty pdfs (uncorrelated)
    GetFEIUncorrelatedPDFs(FEI_uncorrelated_info, CHG_FEI_uncorrelated, MIX_FEI_uncorrelated, Signal_MXs1_FEI_uncorrelated, Signal_MXs2_FEI_uncorrelated, Signal_MXs3_FEI_uncorrelated);

    // get KID uncertainty pdfs (correlated)
    int NPDFs_KID = GetPIDcorrelatedPDFs(KID_correlated_info, CHG_nominal, MIX_nominal, UUBAR_nominal, DDBAR_nominal, SSBAR_nominal, CHARM_nominal, Signal_MXs1_nominal, Signal_MXs2_nominal, Signal_MXs3_nominal, &CHG_KID_correlated, &MIX_KID_correlated, &UUBAR_KID_correlated, &DDBAR_KID_correlated, &SSBAR_KID_correlated, &CHARM_KID_correlated, &Signal_MXs1_KID_correlated, &Signal_MXs2_KID_correlated, &Signal_MXs3_KID_correlated, true);
    int NPDFs_PID = GetPIDcorrelatedPDFs(PID_correlated_info, CHG_nominal, MIX_nominal, UUBAR_nominal, DDBAR_nominal, SSBAR_nominal, CHARM_nominal, Signal_MXs1_nominal, Signal_MXs2_nominal, Signal_MXs3_nominal, &CHG_PID_correlated, &MIX_PID_correlated, &UUBAR_PID_correlated, &DDBAR_PID_correlated, &SSBAR_PID_correlated, &CHARM_PID_correlated, &Signal_MXs1_PID_correlated, &Signal_MXs2_PID_correlated, &Signal_MXs3_PID_correlated, false);

    // get KID uncertainty pdfs (uncorrelated)
    GetPIDUncorrelatedPDFs(KID_uncorrelated_info, CHG_KID_uncorrelated, MIX_KID_uncorrelated, UUBAR_KID_uncorrelated, DDBAR_KID_uncorrelated, SSBAR_KID_uncorrelated, CHARM_KID_uncorrelated, Signal_MXs1_KID_uncorrelated, Signal_MXs2_KID_uncorrelated, Signal_MXs3_KID_uncorrelated);
    GetPIDUncorrelatedPDFs(PID_uncorrelated_info, CHG_PID_uncorrelated, MIX_PID_uncorrelated, UUBAR_PID_uncorrelated, DDBAR_PID_uncorrelated, SSBAR_PID_uncorrelated, CHARM_PID_uncorrelated, Signal_MXs1_PID_uncorrelated, Signal_MXs2_PID_uncorrelated, Signal_MXs3_PID_uncorrelated);

    // get BB BR uncertainty pdfs
    for (int dmIndex = 0; dmIndex < NBRdmID(); dmIndex++) {
        GetBBBRPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs1_BRs.at(dmIndex), "Bplus", "SIGNAL", GetBRdmID(dmIndex), true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 1);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs1_BRs.at(dmIndex), "Bplus", "SIGNAL", GetBRdmID(dmIndex), true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 1);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs1_BRs.at(dmIndex), "Bplus", "SIGNAL", GetBRdmID(dmIndex), true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 1);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs1_BRs.at(dmIndex), "Bzero", "SIGNAL", GetBRdmID(dmIndex), true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 1);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs1_BRs.at(dmIndex), "Bzero", "SIGNAL", GetBRdmID(dmIndex), true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 1);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs1_BRs.at(dmIndex), "Bzero", "SIGNAL", GetBRdmID(dmIndex), true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 1);

        GetBBBRPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs1_BRs.at(dmIndex + NBRdmID()), "Bplus", "SIGNAL", GetBRdmID(dmIndex), false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 1);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs1_BRs.at(dmIndex + NBRdmID()), "Bplus", "SIGNAL", GetBRdmID(dmIndex), false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 1);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs1_BRs.at(dmIndex + NBRdmID()), "Bplus", "SIGNAL", GetBRdmID(dmIndex), false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 1);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs1_BRs.at(dmIndex + NBRdmID()), "Bzero", "SIGNAL", GetBRdmID(dmIndex), false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 1);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs1_BRs.at(dmIndex + NBRdmID()), "Bzero", "SIGNAL", GetBRdmID(dmIndex), false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 1);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs1_BRs.at(dmIndex + NBRdmID()), "Bzero", "SIGNAL", GetBRdmID(dmIndex), false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 1);
    }

    for (int dmIndex = 0; dmIndex < NBRdmID(); dmIndex++) {
        GetBBBRPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs2_BRs.at(dmIndex), "Bplus", "SIGNAL", GetBRdmID(dmIndex), true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 2);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs2_BRs.at(dmIndex), "Bplus", "SIGNAL", GetBRdmID(dmIndex), true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 2);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs2_BRs.at(dmIndex), "Bplus", "SIGNAL", GetBRdmID(dmIndex), true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 2);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs2_BRs.at(dmIndex), "Bzero", "SIGNAL", GetBRdmID(dmIndex), true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 2);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs2_BRs.at(dmIndex), "Bzero", "SIGNAL", GetBRdmID(dmIndex), true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 2);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs2_BRs.at(dmIndex), "Bzero", "SIGNAL", GetBRdmID(dmIndex), true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 2);

        GetBBBRPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs2_BRs.at(dmIndex + NBRdmID()), "Bplus", "SIGNAL", GetBRdmID(dmIndex), false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 2);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs2_BRs.at(dmIndex + NBRdmID()), "Bplus", "SIGNAL", GetBRdmID(dmIndex), false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 2);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs2_BRs.at(dmIndex + NBRdmID()), "Bplus", "SIGNAL", GetBRdmID(dmIndex), false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 2);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs2_BRs.at(dmIndex + NBRdmID()), "Bzero", "SIGNAL", GetBRdmID(dmIndex), false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 2);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs2_BRs.at(dmIndex + NBRdmID()), "Bzero", "SIGNAL", GetBRdmID(dmIndex), false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 2);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs2_BRs.at(dmIndex + NBRdmID()), "Bzero", "SIGNAL", GetBRdmID(dmIndex), false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 2);
    }

    for (int dmIndex = 0; dmIndex < NBRdmID(); dmIndex++) {
        GetBBBRPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_BRs.at(dmIndex), "Bplus", "SIGNAL", GetBRdmID(dmIndex), true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_BRs.at(dmIndex), "Bplus", "SIGNAL", GetBRdmID(dmIndex), true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_BRs.at(dmIndex), "Bplus", "SIGNAL", GetBRdmID(dmIndex), true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_BRs.at(dmIndex), "Bzero", "SIGNAL", GetBRdmID(dmIndex), true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_BRs.at(dmIndex), "Bzero", "SIGNAL", GetBRdmID(dmIndex), true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_BRs.at(dmIndex), "Bzero", "SIGNAL", GetBRdmID(dmIndex), true, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

        GetBBBRPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_BRs.at(dmIndex + NBRdmID()), "Bplus", "SIGNAL", GetBRdmID(dmIndex), false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_BRs.at(dmIndex + NBRdmID()), "Bplus", "SIGNAL", GetBRdmID(dmIndex), false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_BRs.at(dmIndex + NBRdmID()), "Bplus", "SIGNAL", GetBRdmID(dmIndex), false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_BRs.at(dmIndex + NBRdmID()), "Bzero", "SIGNAL", GetBRdmID(dmIndex), false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_BRs.at(dmIndex + NBRdmID()), "Bzero", "SIGNAL", GetBRdmID(dmIndex), false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
        GetBBBRPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_BRs.at(dmIndex + NBRdmID()), "Bzero", "SIGNAL", GetBRdmID(dmIndex), false, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);
    }

    for (int dmIndex = 0; dmIndex < NBRdmID(); dmIndex++) {
        GetBBBRPDFs(MC_dirname_CHG, "root", CHG_BRs.at(dmIndex), "Bplus", "CHG", GetBRdmID(dmIndex), true, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise", 0);
        GetBBBRPDFs(MC_dirname_CHG, "root", CHG_BRs.at(dmIndex + NBRdmID()), "Bplus", "CHG", GetBRdmID(dmIndex), false, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise", 0);
    }

    for (int dmIndex = 0; dmIndex < NBRdmID(); dmIndex++) {
        GetBBBRPDFs(MC_dirname_MIX, "root", MIX_BRs.at(dmIndex), "Bplus", "MIX", GetBRdmID(dmIndex), true, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise", 0);
        GetBBBRPDFs(MC_dirname_MIX, "root", MIX_BRs.at(dmIndex + NBRdmID()), "Bplus", "MIX", GetBRdmID(dmIndex), false, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise", 0);
    }

    // get pi0 uncertainty pdfs (correlated)
    int NPDFs_pi0 = Getpi0correlatedPDFs(pi0_correlated_info, CHG_nominal, MIX_nominal, UUBAR_nominal, DDBAR_nominal, SSBAR_nominal, CHARM_nominal, Signal_MXs1_nominal, Signal_MXs2_nominal, Signal_MXs3_nominal, &CHG_pi0_correlated, &MIX_pi0_correlated, &UUBAR_pi0_correlated, &DDBAR_pi0_correlated, &SSBAR_pi0_correlated, &CHARM_pi0_correlated, &Signal_MXs1_pi0_correlated, &Signal_MXs2_pi0_correlated, &Signal_MXs3_pi0_correlated);

    // get pi0 uncertainty pdfs (uncorrelated)
    Getpi0UncorrelatedPDFs(pi0_uncorrelated_info, CHG_pi0_uncorrelated, MIX_pi0_uncorrelated, UUBAR_pi0_uncorrelated, DDBAR_pi0_uncorrelated, SSBAR_pi0_uncorrelated, CHARM_pi0_uncorrelated, Signal_MXs1_pi0_uncorrelated, Signal_MXs2_pi0_uncorrelated, Signal_MXs3_pi0_uncorrelated);

    // get Kff uncertainty pdfs
    TH1D* Signal_MXs1_Kff_temp = new TH1D("Signal_MXs1_Kff_temp", "Signal_MXs1_Kff_temp", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_Kff_array[7] = { Signal_MXs1_Kff_temp, Signal_MXs1_Kff1_p, Signal_MXs1_Kff1_m, Signal_MXs1_Kff2_p, Signal_MXs1_Kff2_m, Signal_MXs1_Kff3_p, Signal_MXs1_Kff3_m };
    double Correction_factor_BR_Kplus[7] = { 0.0 };
    double Correction_factor_BR_Kzero[7] = { 0.0 };

    GetKffPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs1_Kff_array, Correction_factor_BR_Kplus, "Bplus", 1, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 1);
    for (int i = 0; i < 7; i++) GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs1_Kff_array[i], "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 1);
    for (int i = 0; i < 7; i++) {
        const double non_res_Xsu_correction_factor = (BR_Xsu_nonresonant_nunubar + BR_Kplus_nunubar - BR_Kplus_nunubar * Correction_factor_BR_Kplus[i]) / BR_Xsu_nonresonant_nunubar;
        GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs1_Kff_array[i], "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu") * non_res_Xsu_correction_factor, "B2Xsnunu", 1);
    }
    GetKffPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs1_Kff_array, Correction_factor_BR_Kzero, "Bzero", 0, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 1);
    for (int i = 0; i < 7; i++) GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs1_Kff_array[i], "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 1);
    for (int i = 0; i < 7; i++) {
        const double non_res_Xsd_correction_factor = (BR_Xsd_nonresonant_nunubar + BR_K0_nunubar - BR_K0_nunubar * Correction_factor_BR_Kzero[i]) / BR_Xsd_nonresonant_nunubar;
        GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs1_Kff_array[i], "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu") * non_res_Xsd_correction_factor, "B02Xsnunu", 1);
    }

    TH1D* Signal_MXs3_Kff_temp = new TH1D("Signal_MXs3_Kff_temp", "Signal_MXs3_Kff_temp", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kff_array[7] = { Signal_MXs3_Kff_temp, Signal_MXs3_Kff1_p, Signal_MXs3_Kff1_m, Signal_MXs3_Kff2_p, Signal_MXs3_Kff2_m, Signal_MXs3_Kff3_p, Signal_MXs3_Kff3_m };
    for (int i = 0; i < 7; i++) {
        Correction_factor_BR_Kplus[i] = 0.0;
        Correction_factor_BR_Kzero[i] = 0.0;
    }

    GetKffPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_Kff_array, Correction_factor_BR_Kplus, "Bplus", 1, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    for (int i = 0; i < 7; i++) GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_Kff_array[i], "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    for (int i = 0; i < 7; i++) {
        const double non_res_Xsu_correction_factor = (BR_Xsu_nonresonant_nunubar + BR_Kplus_nunubar - BR_Kplus_nunubar * Correction_factor_BR_Kplus[i]) / BR_Xsu_nonresonant_nunubar;
        GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_Kff_array[i], "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu") * non_res_Xsu_correction_factor, "B2Xsnunu", 3);
    }
    GetKffPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_Kff_array, Correction_factor_BR_Kzero, "Bzero", 0, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    for (int i = 0; i < 7; i++) GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_Kff_array[i], "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    for (int i = 0; i < 7; i++) {
        const double non_res_Xsd_correction_factor = (BR_Xsd_nonresonant_nunubar + BR_K0_nunubar - BR_K0_nunubar * Correction_factor_BR_Kzero[i]) / BR_Xsd_nonresonant_nunubar;
        GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_Kff_array[i], "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu") * non_res_Xsd_correction_factor, "B02Xsnunu", 3);
    }

    // get Kstar ff uncertainty pdfs
    TH1D* Signal_MXs2_Kstarff_temp = new TH1D("Signal_MXs2_Kstarff_temp", "Signal_MXs2_Kstarff_temp", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_Kstarff_array[19] = { Signal_MXs2_Kstarff_temp,
        Signal_MXs2_Kstarff1_p, Signal_MXs2_Kstarff1_m,
        Signal_MXs2_Kstarff2_p, Signal_MXs2_Kstarff2_m,
        Signal_MXs2_Kstarff3_p, Signal_MXs2_Kstarff3_m,
        Signal_MXs2_Kstarff4_p, Signal_MXs2_Kstarff4_m,
        Signal_MXs2_Kstarff5_p, Signal_MXs2_Kstarff5_m,
        Signal_MXs2_Kstarff6_p, Signal_MXs2_Kstarff6_m,
        Signal_MXs2_Kstarff7_p, Signal_MXs2_Kstarff7_m,
        Signal_MXs2_Kstarff8_p, Signal_MXs2_Kstarff8_m,
        Signal_MXs2_Kstarff9_p, Signal_MXs2_Kstarff9_m
    };
    double Correction_factor_BR_Kstarplus[19] = { 0.0 };
    double Correction_factor_BR_Kstarzero[19] = { 0.0 };

    for (int i = 0; i < 19; i++) GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs2_Kstarff_array[i], "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 2);
    GetKstarffPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs2_Kstarff_array, Correction_factor_BR_Kstarplus, "Bplus", 1, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), 2);
    for (int i = 0; i < 19; i++) {
        const double non_res_Xsu_correction_factor = (BR_Xsu_nonresonant_nunubar + BR_Kplusstar_nunubar - BR_Kplusstar_nunubar * Correction_factor_BR_Kstarplus[i]) / BR_Xsu_nonresonant_nunubar;
        GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs2_Kstarff_array[i], "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu") * non_res_Xsu_correction_factor, "B2Xsnunu", 2);
    }
    for (int i = 0; i < 19; i++) GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs2_Kstarff_array[i], "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 2);
    GetKstarffPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs2_Kstarff_array, Correction_factor_BR_Kstarzero, "Bzero", 0, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), 2);
    for (int i = 0; i < 19; i++) {
        const double non_res_Xsd_correction_factor = (BR_Xsd_nonresonant_nunubar + BR_K0star_nunubar - BR_K0star_nunubar * Correction_factor_BR_Kstarzero[i]) / BR_Xsd_nonresonant_nunubar;
        GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs2_Kstarff_array[i], "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu") * non_res_Xsd_correction_factor, "B02Xsnunu", 2);
    }

    TH1D* Signal_MXs3_Kstarff_temp = new TH1D("Signal_MXs3_Kstarff_temp", "Signal_MXs3_Kstarff_temp", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_Kstarff_array[19] = { Signal_MXs3_Kstarff_temp,
        Signal_MXs3_Kstarff1_p, Signal_MXs3_Kstarff1_m,
        Signal_MXs3_Kstarff2_p, Signal_MXs3_Kstarff2_m,
        Signal_MXs3_Kstarff3_p, Signal_MXs3_Kstarff3_m,
        Signal_MXs3_Kstarff4_p, Signal_MXs3_Kstarff4_m,
        Signal_MXs3_Kstarff5_p, Signal_MXs3_Kstarff5_m,
        Signal_MXs3_Kstarff6_p, Signal_MXs3_Kstarff6_m,
        Signal_MXs3_Kstarff7_p, Signal_MXs3_Kstarff7_m,
        Signal_MXs3_Kstarff8_p, Signal_MXs3_Kstarff8_m,
        Signal_MXs3_Kstarff9_p, Signal_MXs3_Kstarff9_m
    };
    for (int i = 0; i < 19; i++) {
        Correction_factor_BR_Kstarplus[i] = 0.0;
        Correction_factor_BR_Kstarzero[i] = 0.0;
    }

    for (int i = 0; i < 19; i++) GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_Kstarff_array[i], "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetKstarffPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_Kstarff_array, Correction_factor_BR_Kstarplus, "Bplus", 1, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), 3);
    for (int i = 0; i < 19; i++) {
        const double non_res_Xsu_correction_factor = (BR_Xsu_nonresonant_nunubar + BR_Kplusstar_nunubar - BR_Kplusstar_nunubar * Correction_factor_BR_Kstarplus[i]) / BR_Xsu_nonresonant_nunubar;
        GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_Kstarff_array[i], "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu") * non_res_Xsu_correction_factor, "B2Xsnunu", 3);
    }
    for (int i = 0; i < 19; i++) GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_Kstarff_array[i], "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetKstarffPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_Kstarff_array, Correction_factor_BR_Kstarzero, "Bzero", 0, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), 3);
    for (int i = 0; i < 19; i++) {
        const double non_res_Xsd_correction_factor = (BR_Xsd_nonresonant_nunubar + BR_K0star_nunubar - BR_K0star_nunubar * Correction_factor_BR_Kstarzero[i]) / BR_Xsd_nonresonant_nunubar;
        GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_Kstarff_array[i], "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu") * non_res_Xsd_correction_factor, "B02Xsnunu", 3);
    }

    // get fraction uncertainty PDFs
    for (int i = 0; i < 7; i++) {
        if (Correction_factor_BR_Kplus[i] < MyEPSILON) {
            printf("Correction_factor_BR_Kplus is equal to zero. Maybe you forgot to initialize it.\n");
            exit(1);
        }
    }
    double Positive_relative_uncertainty_BR_Kplus[3] = { 0.0 };
    double Negative_relative_uncertainty_BR_Kplus[3] = { 0.0 };
    if ((Correction_factor_BR_Kplus[1] > 1.0) && (Correction_factor_BR_Kplus[2] < 1.0)) {
        Positive_relative_uncertainty_BR_Kplus[0] = Correction_factor_BR_Kplus[1] - 1.0;
        Negative_relative_uncertainty_BR_Kplus[0] = 1.0 - Correction_factor_BR_Kplus[2];
    }
    else if ((Correction_factor_BR_Kplus[1] < 1.0) && (Correction_factor_BR_Kplus[2] > 1.0)) {
        Positive_relative_uncertainty_BR_Kplus[0] = 1.0 - Correction_factor_BR_Kplus[1];
        Negative_relative_uncertainty_BR_Kplus[0] = Correction_factor_BR_Kplus[2] - 1.0;
    }
    else {
        printf("Correction_factor_BR_Kplus[1]: %lf\n", Correction_factor_BR_Kplus[1]);
        printf("Correction_factor_BR_Kplus[2]: %lf\n", Correction_factor_BR_Kplus[2]);
        printf("one should be larger than 1.0 and the other should be smaller than 1.0");
        exit(1);
    }
    if ((Correction_factor_BR_Kplus[3] > 1.0) && (Correction_factor_BR_Kplus[4] < 1.0)) {
        Positive_relative_uncertainty_BR_Kplus[1] = Correction_factor_BR_Kplus[3] - 1.0;
        Negative_relative_uncertainty_BR_Kplus[1] = 1.0 - Correction_factor_BR_Kplus[4];
    }
    else if ((Correction_factor_BR_Kplus[3] < 1.0) && (Correction_factor_BR_Kplus[4] > 1.0)) {
        Positive_relative_uncertainty_BR_Kplus[1] = 1.0 - Correction_factor_BR_Kplus[3];
        Negative_relative_uncertainty_BR_Kplus[1] = Correction_factor_BR_Kplus[4] - 1.0;
    }
    else {
        printf("Correction_factor_BR_Kplus[3]: %lf\n", Correction_factor_BR_Kplus[3]);
        printf("Correction_factor_BR_Kplus[4]: %lf\n", Correction_factor_BR_Kplus[4]);
        printf("one should be larger than 1.0 and the other should be smaller than 1.0");
        exit(1);
    }
    if ((Correction_factor_BR_Kplus[5] > 1.0) && (Correction_factor_BR_Kplus[6] < 1.0)) {
        Positive_relative_uncertainty_BR_Kplus[2] = Correction_factor_BR_Kplus[5] - 1.0;
        Negative_relative_uncertainty_BR_Kplus[2] = 1.0 - Correction_factor_BR_Kplus[6];
    }
    else if ((Correction_factor_BR_Kplus[5] < 1.0) && (Correction_factor_BR_Kplus[6] > 1.0)) {
        Positive_relative_uncertainty_BR_Kplus[2] = 1.0 - Correction_factor_BR_Kplus[5];
        Negative_relative_uncertainty_BR_Kplus[2] = Correction_factor_BR_Kplus[6] - 1.0;
    }
    else {
        printf("Correction_factor_BR_Kplus[5]: %lf\n", Correction_factor_BR_Kplus[5]);
        printf("Correction_factor_BR_Kplus[6]: %lf\n", Correction_factor_BR_Kplus[6]);
        printf("one should be larger than 1.0 and the other should be smaller than 1.0");
        exit(1);
    }
    double total_relative_uncertainty_BR_Kplus = Sigma_BR_Kplus_nunubar / BR_Kplus_nunubar;
    double total_parameter_positive_relative_uncertainty_BR_Kplus = std::sqrt(
        total_relative_uncertainty_BR_Kplus * total_relative_uncertainty_BR_Kplus
        - Positive_relative_uncertainty_BR_Kplus[0] * Positive_relative_uncertainty_BR_Kplus[0]
        - Positive_relative_uncertainty_BR_Kplus[1] * Positive_relative_uncertainty_BR_Kplus[1]
        - Positive_relative_uncertainty_BR_Kplus[2] * Positive_relative_uncertainty_BR_Kplus[2]
    );
    double total_parameter_negative_relative_uncertainty_BR_Kplus = std::sqrt(
        total_relative_uncertainty_BR_Kplus * total_relative_uncertainty_BR_Kplus
        - Negative_relative_uncertainty_BR_Kplus[0] * Negative_relative_uncertainty_BR_Kplus[0]
        - Negative_relative_uncertainty_BR_Kplus[1] * Negative_relative_uncertainty_BR_Kplus[1]
        - Negative_relative_uncertainty_BR_Kplus[2] * Negative_relative_uncertainty_BR_Kplus[2]
    );
    printf("total relative uncertainty: %lf\n", total_relative_uncertainty_BR_Kplus);
    printf("total positive parameter relative uncertainty: %lf\n", total_parameter_positive_relative_uncertainty_BR_Kplus);
    printf("total negative parameter relative uncertainty: %lf\n", total_parameter_negative_relative_uncertainty_BR_Kplus);
    if ((total_parameter_positive_relative_uncertainty_BR_Kplus < 0) || (total_parameter_negative_relative_uncertainty_BR_Kplus < 0)) {
        printf("total parameter relative uncertainy is negative\n");
        exit(1);
    }
    double Sigma_BR_Kplus_nunubar_parameter_positive = BR_Kplus_nunubar * total_parameter_positive_relative_uncertainty_BR_Kplus;
    double Sigma_BR_Kplus_nunubar_parameter_negative = BR_Kplus_nunubar * total_parameter_negative_relative_uncertainty_BR_Kplus;
    double Sigma_BR_K0_nunubar_parameter_positive = BR_K0_nunubar * total_parameter_positive_relative_uncertainty_BR_Kplus;
    double Sigma_BR_K0_nunubar_parameter_negative = BR_K0_nunubar * total_parameter_negative_relative_uncertainty_BR_Kplus;

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs1_Kfrac_p, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu") * ((BR_Kplus_nunubar + Sigma_BR_Kplus_nunubar_parameter_positive) / BR_Kplus_nunubar), "B2Knunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs1_Kfrac_p, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs1_Kfrac_p, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu") * ((BR_Xsu_nonresonant_nunubar - Sigma_BR_Kplus_nunubar_parameter_positive) / BR_Xsu_nonresonant_nunubar), "B2Xsnunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs1_Kfrac_p, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu") * ((BR_K0_nunubar + Sigma_BR_K0_nunubar_parameter_positive) / BR_K0_nunubar), "B02K0nunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs1_Kfrac_p, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs1_Kfrac_p, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu") * ((BR_Xsd_nonresonant_nunubar - Sigma_BR_K0_nunubar_parameter_positive) / BR_Xsd_nonresonant_nunubar), "B02Xsnunu", 1);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs1_Kfrac_m, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu") * ((BR_Kplus_nunubar - Sigma_BR_Kplus_nunubar_parameter_negative) / BR_Kplus_nunubar), "B2Knunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs1_Kfrac_m, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs1_Kfrac_m, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu") * ((BR_Xsu_nonresonant_nunubar + Sigma_BR_Kplus_nunubar_parameter_negative) / BR_Xsu_nonresonant_nunubar), "B2Xsnunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs1_Kfrac_m, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu") * ((BR_K0_nunubar - Sigma_BR_K0_nunubar_parameter_negative) / BR_K0_nunubar), "B02K0nunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs1_Kfrac_m, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs1_Kfrac_m, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu") * ((BR_Xsd_nonresonant_nunubar + Sigma_BR_K0_nunubar_parameter_negative) / BR_Xsd_nonresonant_nunubar), "B02Xsnunu", 1);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_Kfrac_p, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu") * ((BR_Kplus_nunubar + Sigma_BR_Kplus_nunubar_parameter_positive) / BR_Kplus_nunubar), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_Kfrac_p, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_Kfrac_p, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu") * ((BR_Xsu_nonresonant_nunubar - Sigma_BR_Kplus_nunubar_parameter_positive) / BR_Xsu_nonresonant_nunubar), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_Kfrac_p, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu") * ((BR_K0_nunubar + Sigma_BR_K0_nunubar_parameter_positive) / BR_K0_nunubar), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_Kfrac_p, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_Kfrac_p, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu") * ((BR_Xsd_nonresonant_nunubar - Sigma_BR_K0_nunubar_parameter_positive) / BR_Xsd_nonresonant_nunubar), "B02Xsnunu", 3);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_Kfrac_m, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu") * ((BR_Kplus_nunubar - Sigma_BR_Kplus_nunubar_parameter_negative) / BR_Kplus_nunubar), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_Kfrac_m, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_Kfrac_m, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu") * ((BR_Xsu_nonresonant_nunubar + Sigma_BR_Kplus_nunubar_parameter_negative) / BR_Xsu_nonresonant_nunubar), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_Kfrac_m, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu") * ((BR_K0_nunubar - Sigma_BR_K0_nunubar_parameter_negative) / BR_K0_nunubar), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_Kfrac_m, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_Kfrac_m, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu") * ((BR_Xsd_nonresonant_nunubar + Sigma_BR_K0_nunubar_parameter_negative) / BR_Xsd_nonresonant_nunubar), "B02Xsnunu", 3);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs2_Kstarfrac_p, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs2_Kstarfrac_p, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu") * ((BR_Kplusstar_nunubar + Sigma_BR_Kplusstar_nunubar) / BR_Kplusstar_nunubar), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs2_Kstarfrac_p, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu") * ((BR_Xsu_nonresonant_nunubar - Sigma_BR_Kplusstar_nunubar) / BR_Xsu_nonresonant_nunubar), "B2Xsnunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs2_Kstarfrac_p, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs2_Kstarfrac_p, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu") * ((BR_K0star_nunubar + Sigma_BR_K0star_nunubar) / BR_K0star_nunubar), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs2_Kstarfrac_p, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu") * ((BR_Xsd_nonresonant_nunubar - Sigma_BR_K0star_nunubar) / BR_Xsd_nonresonant_nunubar), "B02Xsnunu", 2);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs2_Kstarfrac_m, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs2_Kstarfrac_m, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu") * ((BR_Kplusstar_nunubar - Sigma_BR_Kplusstar_nunubar) / BR_Kplusstar_nunubar), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs2_Kstarfrac_m, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu") * ((BR_Xsu_nonresonant_nunubar + Sigma_BR_Kplusstar_nunubar) / BR_Xsu_nonresonant_nunubar), "B2Xsnunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs2_Kstarfrac_m, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs2_Kstarfrac_m, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu") * ((BR_K0star_nunubar - Sigma_BR_K0star_nunubar) / BR_K0star_nunubar), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs2_Kstarfrac_m, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu") * ((BR_Xsd_nonresonant_nunubar + Sigma_BR_K0star_nunubar) / BR_Xsd_nonresonant_nunubar), "B02Xsnunu", 2);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_Kstarfrac_p, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_Kstarfrac_p, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu") * ((BR_Kplusstar_nunubar + Sigma_BR_Kplusstar_nunubar) / BR_Kplusstar_nunubar), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_Kstarfrac_p, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu") * ((BR_Xsu_nonresonant_nunubar - Sigma_BR_Kplusstar_nunubar) / BR_Xsu_nonresonant_nunubar), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_Kstarfrac_p, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_Kstarfrac_p, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu") * ((BR_K0star_nunubar + Sigma_BR_K0star_nunubar) / BR_K0star_nunubar), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_Kstarfrac_p, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu") * ((BR_Xsd_nonresonant_nunubar - Sigma_BR_K0star_nunubar) / BR_Xsd_nonresonant_nunubar), "B02Xsnunu", 3);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_Kstarfrac_m, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_Kstarfrac_m, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu") * ((BR_Kplusstar_nunubar - Sigma_BR_Kplusstar_nunubar) / BR_Kplusstar_nunubar), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_Kstarfrac_m, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu") * ((BR_Xsu_nonresonant_nunubar + Sigma_BR_Kplusstar_nunubar) / BR_Xsu_nonresonant_nunubar), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_Kstarfrac_m, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_Kstarfrac_m, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu") * ((BR_K0star_nunubar - Sigma_BR_K0star_nunubar) / BR_K0star_nunubar), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_Kstarfrac_m, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu") * ((BR_Xsd_nonresonant_nunubar + Sigma_BR_K0star_nunubar) / BR_Xsd_nonresonant_nunubar), "B02Xsnunu", 3);

    // get fragmentation uncertainty pdfs
    for (int sign = 0; sign < 2; sign++) {
        for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
            for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
                int temp_index = MxsBin * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma) + Category + sign * corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma) * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma);
                bool IsItUp = true;
                if (name_Signal_MXs1_Fragmentation.at(temp_index).find("_p") != string::npos) IsItUp = true;
                else if (name_Signal_MXs1_Fragmentation.at(temp_index).find("_m") != string::npos) IsItUp = false;
                else {
                    printf("[ERROR] unexpected error!\n");
                    exit(1);
                }

                GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs1_Fragmentaions.at(temp_index), "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 1);
                GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs1_Fragmentaions.at(temp_index), "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 1);
                GetFragmentationPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs1_Fragmentaions.at(temp_index), "Bplus", "SIGNAL", MxsBin, Category, IsItUp, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 1);
                GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs1_Fragmentaions.at(temp_index), "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 1);
                GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs1_Fragmentaions.at(temp_index), "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 1);
                GetFragmentationPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs1_Fragmentaions.at(temp_index), "Bzero", "SIGNAL", MxsBin, Category, IsItUp, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 1);

            }
        }
    }

    for (int sign = 0; sign < 2; sign++) {
        for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
            for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
                int temp_index = MxsBin * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma) + Category + sign * corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma) * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma);
                bool IsItUp = true;
                if (name_Signal_MXs2_Fragmentation.at(temp_index).find("_p") != string::npos) IsItUp = true;
                else if (name_Signal_MXs2_Fragmentation.at(temp_index).find("_m") != string::npos) IsItUp = false;
                else {
                    printf("[ERROR] unexpected error!\n");
                    exit(1);
                }

                GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs2_Fragmentaions.at(temp_index), "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 2);
                GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs2_Fragmentaions.at(temp_index), "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 2);
                GetFragmentationPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs2_Fragmentaions.at(temp_index), "Bplus", "SIGNAL", MxsBin, Category, IsItUp, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 2);
                GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs2_Fragmentaions.at(temp_index), "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 2);
                GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs2_Fragmentaions.at(temp_index), "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 2);
                GetFragmentationPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs2_Fragmentaions.at(temp_index), "Bzero", "SIGNAL", MxsBin, Category, IsItUp, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 2);

            }
        }
    }

    for (int sign = 0; sign < 2; sign++) {
        for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
            for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
                int temp_index = MxsBin * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma) + Category + sign * corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma) * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma);
                bool IsItUp = true;
                if (name_Signal_MXs3_Fragmentation.at(temp_index).find("_p") != string::npos) IsItUp = true;
                else if (name_Signal_MXs3_Fragmentation.at(temp_index).find("_m") != string::npos) IsItUp = false;
                else {
                    printf("[ERROR] unexpected error!\n");
                    exit(1);
                }

                GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_Fragmentaions.at(temp_index), "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
                GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_Fragmentaions.at(temp_index), "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
                GetFragmentationPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_Fragmentaions.at(temp_index), "Bplus", "SIGNAL", MxsBin, Category, IsItUp, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
                GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_Fragmentaions.at(temp_index), "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
                GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_Fragmentaions.at(temp_index), "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
                GetFragmentationPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_Fragmentaions.at(temp_index), "Bzero", "SIGNAL", MxsBin, Category, IsItUp, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

            }
        }
    }

    for (int sign = 0; sign < 2; sign++) {
        for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
            for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
                int temp_index = MxsBin * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma) + Category + sign * corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma) * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma);
                bool IsItUp = true;
                if (name_CHG_Fragmentation.at(temp_index).find("_p") != string::npos) IsItUp = true;
                else if (name_CHG_Fragmentation.at(temp_index).find("_m") != string::npos) IsItUp = false;
                else {
                    printf("[ERROR] unexpected error!\n");
                    exit(1);
                }

                GetPDFs(MC_dirname_CHG, "root", CHG_Fragmentaions.at(temp_index), "Bplus", "CHG", PDFtype::nominal, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise", 0);

            }
        }
    }

    for (int sign = 0; sign < 2; sign++) {
        for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
            for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
                int temp_index = MxsBin * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma) + Category + sign * corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma) * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma);
                bool IsItUp = true;
                if (name_MIX_Fragmentation.at(temp_index).find("_p") != string::npos) IsItUp = true;
                else if (name_MIX_Fragmentation.at(temp_index).find("_m") != string::npos) IsItUp = false;
                else {
                    printf("[ERROR] unexpected error!\n");
                    exit(1);
                }

                GetPDFs(MC_dirname_MIX, "root", MIX_Fragmentaions.at(temp_index), "Bzero", "MIX", PDFtype::nominal, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise", 0);

            }
        }
    }

    // get pf uncertainty pdfs
    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_pf_p, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_pf_p, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_pf_p, "Bplus", "SIGNAL", PDFtype::pfXsuUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_pf_p, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_pf_p, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_pf_p, "Bzero", "SIGNAL", PDFtype::pfXsdUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_pf_m, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_pf_m, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_pf_m, "Bplus", "SIGNAL", PDFtype::pfXsuDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_pf_m, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_pf_m, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_pf_m, "Bzero", "SIGNAL", PDFtype::pfXsdDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

    // mb
    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_mb_p, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_mb_p, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_mb_p, "Bplus", "SIGNAL", PDFtype::mbXsuUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_mb_p, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_mb_p, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_mb_p, "Bzero", "SIGNAL", PDFtype::mbXsdUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_mb_m, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_mb_m, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_mb_m, "Bplus", "SIGNAL", PDFtype::mbXsuDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_mb_m, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_mb_m, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_mb_m, "Bzero", "SIGNAL", PDFtype::mbXsdDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

    // transition
    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_transition_p, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_transition_p, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_transition_p, "Bplus", "SIGNAL", PDFtype::TransitionXsuUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_transition_p, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_transition_p, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_transition_p, "Bzero", "SIGNAL", PDFtype::TransitionXsdUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_transition_m, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_transition_m, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_transition_m, "Bplus", "SIGNAL", PDFtype::TransitionXsuDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_transition_m, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_transition_m, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_transition_m, "Bzero", "SIGNAL", PDFtype::TransitionXsdDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

    // mKstar fixed
    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs2_mKstar_p, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs2_mKstar_p, "Bplus", "SIGNAL", PDFtype::KstardeltaKstarplus, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs2_mKstar_p, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs2_mKstar_p, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs2_mKstar_p, "Bzero", "SIGNAL", PDFtype::KstardeltaKstarzero, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs2_mKstar_p, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 2);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_mKstar_p, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_mKstar_p, "Bplus", "SIGNAL", PDFtype::KstardeltaKstarplus, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_mKstar_p, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_mKstar_p, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_mKstar_p, "Bzero", "SIGNAL", PDFtype::KstardeltaKstarzero, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_mKstar_p, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

    GetNegativeChangePDFs(Signal_MXs2_nominal, Signal_MXs2_mKstar_p, Signal_MXs2_mKstar_m);
    GetNegativeChangePDFs(Signal_MXs3_nominal, Signal_MXs3_mKstar_p, Signal_MXs3_mKstar_m);

    // get BDTc uncertainty pdfs
    //double NevtCHGWithBDTc = GetNevtWithBDTc(MC_dirname_CHG, "root", "Bplus", "CHG", ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise", 0);
    //double NevtMIXWithBDTc = GetNevtWithBDTc(MC_dirname_MIX, "root", "Bzero", "MIX", ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise", 0);
    double NevtUUBARWithBDTc = GetNevtWithBDTc(MC_dirname_UUBAR, "root", "Continuum", "UUBAR", ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"), "otherwise", 0);
    double NevtDDBARWithBDTc = GetNevtWithBDTc(MC_dirname_DDBAR, "root", "Continuum", "DDBAR", ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"), "otherwise", 0);
    double NevtSSBARWithBDTc = GetNevtWithBDTc(MC_dirname_SSBAR, "root", "Continuum", "SSBAR", ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"), "otherwise", 0);
    double NevtCHARMWithBDTc = GetNevtWithBDTc(MC_dirname_CHARM, "root", "Continuum", "CHARM", ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"), "otherwise", 0);

    double NormFactor = (NevtUUBARWithoutBDTc + NevtDDBARWithoutBDTc + NevtSSBARWithoutBDTc + NevtCHARMWithoutBDTc) / (NevtUUBARWithBDTc + NevtDDBARWithBDTc + NevtSSBARWithBDTc + NevtCHARMWithBDTc);

    //GetBDTcPDFs(MC_dirname_CHG, "root", CHG_BDTc_p, "Bplus", "CHG", ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise", 0);
    //GetBDTcPDFs(MC_dirname_MIX, "root", MIX_BDTc_p, "Bzero", "MIX", ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise", 0);
    double NevtUUBARWithBDTcWithNorm = GetPDFs(MC_dirname_UUBAR, "root", UUBAR_BDTc_p, "Continuum", "UUBAR", PDFtype::BDTc, ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR") * NormFactor, "otherwise", 0);
    double NevtDDBARWithBDTcWithNorm = GetPDFs(MC_dirname_DDBAR, "root", DDBAR_BDTc_p, "Continuum", "DDBAR", PDFtype::BDTc, ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR") * NormFactor, "otherwise", 0);
    double NevtSSBARWithBDTcWithNorm = GetPDFs(MC_dirname_SSBAR, "root", SSBAR_BDTc_p, "Continuum", "SSBAR", PDFtype::BDTc, ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR") * NormFactor, "otherwise", 0);
    double NevtCHARMWithBDTcWithNorm = GetPDFs(MC_dirname_CHARM, "root", CHARM_BDTc_p, "Continuum", "CHARM", PDFtype::BDTc, ObtainWeight("CHARM", MCTYPE, "validation", "CHARM") * NormFactor, "otherwise", 0);

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
    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs1_Xnn_p, "Bplus", "SIGNAL", PDFtype::XnnBRUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs1_Xnn_p, "Bplus", "SIGNAL", PDFtype::XnnBRUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs1_Xnn_p, "Bplus", "SIGNAL", PDFtype::XnnBRUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs1_Xnn_p, "Bzero", "SIGNAL", PDFtype::XnnBRUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs1_Xnn_p, "Bzero", "SIGNAL", PDFtype::XnnBRUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs1_Xnn_p, "Bzero", "SIGNAL", PDFtype::XnnBRUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 1);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs2_Xnn_p, "Bplus", "SIGNAL", PDFtype::XnnBRUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs2_Xnn_p, "Bplus", "SIGNAL", PDFtype::XnnBRUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs2_Xnn_p, "Bplus", "SIGNAL", PDFtype::XnnBRUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs2_Xnn_p, "Bzero", "SIGNAL", PDFtype::XnnBRUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs2_Xnn_p, "Bzero", "SIGNAL", PDFtype::XnnBRUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs2_Xnn_p, "Bzero", "SIGNAL", PDFtype::XnnBRUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 2);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_Xnn_p, "Bplus", "SIGNAL", PDFtype::XnnBRUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_Xnn_p, "Bplus", "SIGNAL", PDFtype::XnnBRUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_Xnn_p, "Bplus", "SIGNAL", PDFtype::XnnBRUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_Xnn_p, "Bzero", "SIGNAL", PDFtype::XnnBRUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_Xnn_p, "Bzero", "SIGNAL", PDFtype::XnnBRUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_Xnn_p, "Bzero", "SIGNAL", PDFtype::XnnBRUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

    GetPDFs(MC_dirname_CHG, "root", CHG_Xnn_p, "Bplus", "CHG", PDFtype::XnnBRUP, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise", 0);
    GetPDFs(MC_dirname_MIX, "root", MIX_Xnn_p, "Bzero", "MIX", PDFtype::XnnBRUP, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise", 0);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs1_Xnn_m, "Bplus", "SIGNAL", PDFtype::XnnBRDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs1_Xnn_m, "Bplus", "SIGNAL", PDFtype::XnnBRDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs1_Xnn_m, "Bplus", "SIGNAL", PDFtype::XnnBRDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs1_Xnn_m, "Bzero", "SIGNAL", PDFtype::XnnBRDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs1_Xnn_m, "Bzero", "SIGNAL", PDFtype::XnnBRDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs1_Xnn_m, "Bzero", "SIGNAL", PDFtype::XnnBRDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 1);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs2_Xnn_m, "Bplus", "SIGNAL", PDFtype::XnnBRDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs2_Xnn_m, "Bplus", "SIGNAL", PDFtype::XnnBRDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs2_Xnn_m, "Bplus", "SIGNAL", PDFtype::XnnBRDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs2_Xnn_m, "Bzero", "SIGNAL", PDFtype::XnnBRDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs2_Xnn_m, "Bzero", "SIGNAL", PDFtype::XnnBRDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs2_Xnn_m, "Bzero", "SIGNAL", PDFtype::XnnBRDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 2);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_Xnn_m, "Bplus", "SIGNAL", PDFtype::XnnBRDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_Xnn_m, "Bplus", "SIGNAL", PDFtype::XnnBRDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_Xnn_m, "Bplus", "SIGNAL", PDFtype::XnnBRDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_Xnn_m, "Bzero", "SIGNAL", PDFtype::XnnBRDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_Xnn_m, "Bzero", "SIGNAL", PDFtype::XnnBRDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_Xnn_m, "Bzero", "SIGNAL", PDFtype::XnnBRDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

    GetPDFs(MC_dirname_CHG, "root", CHG_Xnn_m, "Bplus", "CHG", PDFtype::XnnBRDOWN, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise", 0);
    GetPDFs(MC_dirname_MIX, "root", MIX_Xnn_m, "Bzero", "MIX", PDFtype::XnnBRDOWN, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise", 0);

    // get multiplicity uncertainty pdfs (correlated)
    int NPDFs_multiplicity = GetmultiplicitycorrelatedPDFs(multiplicity_correlated_info, CHG_nominal, MIX_nominal, UUBAR_nominal, DDBAR_nominal, SSBAR_nominal, CHARM_nominal, Signal_MXs1_nominal, Signal_MXs2_nominal, Signal_MXs3_nominal, &CHG_multiplicity_correlated, &MIX_multiplicity_correlated, &UUBAR_multiplicity_correlated, &DDBAR_multiplicity_correlated, &SSBAR_multiplicity_correlated, &CHARM_multiplicity_correlated, &Signal_MXs1_multiplicity_correlated, &Signal_MXs2_multiplicity_correlated, &Signal_MXs3_multiplicity_correlated);

    // get multiplicity uncertainty pdfs (uncorrelated)
    GetmultiplicityUncorrelatedPDFs(multiplicity_uncorrelated_info, CHG_multiplicity_uncorrelated, MIX_multiplicity_uncorrelated, UUBAR_multiplicity_uncorrelated, DDBAR_multiplicity_uncorrelated, SSBAR_multiplicity_uncorrelated, CHARM_multiplicity_uncorrelated, Signal_MXs1_multiplicity_uncorrelated, Signal_MXs2_multiplicity_uncorrelated, Signal_MXs3_multiplicity_uncorrelated);

    // B -> [D -> X KL0] anything
    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs1_BtoDtoXKL_p, "Bplus", "SIGNAL", PDFtype::BRBtoDtoXKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs1_BtoDtoXKL_p, "Bplus", "SIGNAL", PDFtype::BRBtoDtoXKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs1_BtoDtoXKL_p, "Bplus", "SIGNAL", PDFtype::BRBtoDtoXKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs1_BtoDtoXKL_p, "Bzero", "SIGNAL", PDFtype::BRBtoDtoXKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs1_BtoDtoXKL_p, "Bzero", "SIGNAL", PDFtype::BRBtoDtoXKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs1_BtoDtoXKL_p, "Bzero", "SIGNAL", PDFtype::BRBtoDtoXKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 1);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs2_BtoDtoXKL_p, "Bplus", "SIGNAL", PDFtype::BRBtoDtoXKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs2_BtoDtoXKL_p, "Bplus", "SIGNAL", PDFtype::BRBtoDtoXKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs2_BtoDtoXKL_p, "Bplus", "SIGNAL", PDFtype::BRBtoDtoXKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs2_BtoDtoXKL_p, "Bzero", "SIGNAL", PDFtype::BRBtoDtoXKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs2_BtoDtoXKL_p, "Bzero", "SIGNAL", PDFtype::BRBtoDtoXKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs2_BtoDtoXKL_p, "Bzero", "SIGNAL", PDFtype::BRBtoDtoXKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 2);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_BtoDtoXKL_p, "Bplus", "SIGNAL", PDFtype::BRBtoDtoXKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_BtoDtoXKL_p, "Bplus", "SIGNAL", PDFtype::BRBtoDtoXKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_BtoDtoXKL_p, "Bplus", "SIGNAL", PDFtype::BRBtoDtoXKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_BtoDtoXKL_p, "Bzero", "SIGNAL", PDFtype::BRBtoDtoXKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_BtoDtoXKL_p, "Bzero", "SIGNAL", PDFtype::BRBtoDtoXKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_BtoDtoXKL_p, "Bzero", "SIGNAL", PDFtype::BRBtoDtoXKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

    GetPDFs(MC_dirname_CHG, "root", CHG_BtoDtoXKL_p, "Bplus", "CHG", PDFtype::BRBtoDtoXKLUP, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise", 0);
    GetPDFs(MC_dirname_MIX, "root", MIX_BtoDtoXKL_p, "Bzero", "MIX", PDFtype::BRBtoDtoXKLUP, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise", 0);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs1_BtoDtoXKL_m, "Bplus", "SIGNAL", PDFtype::BRBtoDtoXKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs1_BtoDtoXKL_m, "Bplus", "SIGNAL", PDFtype::BRBtoDtoXKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs1_BtoDtoXKL_m, "Bplus", "SIGNAL", PDFtype::BRBtoDtoXKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs1_BtoDtoXKL_m, "Bzero", "SIGNAL", PDFtype::BRBtoDtoXKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs1_BtoDtoXKL_m, "Bzero", "SIGNAL", PDFtype::BRBtoDtoXKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs1_BtoDtoXKL_m, "Bzero", "SIGNAL", PDFtype::BRBtoDtoXKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 1);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs2_BtoDtoXKL_m, "Bplus", "SIGNAL", PDFtype::BRBtoDtoXKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs2_BtoDtoXKL_m, "Bplus", "SIGNAL", PDFtype::BRBtoDtoXKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs2_BtoDtoXKL_m, "Bplus", "SIGNAL", PDFtype::BRBtoDtoXKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs2_BtoDtoXKL_m, "Bzero", "SIGNAL", PDFtype::BRBtoDtoXKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs2_BtoDtoXKL_m, "Bzero", "SIGNAL", PDFtype::BRBtoDtoXKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs2_BtoDtoXKL_m, "Bzero", "SIGNAL", PDFtype::BRBtoDtoXKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 2);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_BtoDtoXKL_m, "Bplus", "SIGNAL", PDFtype::BRBtoDtoXKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_BtoDtoXKL_m, "Bplus", "SIGNAL", PDFtype::BRBtoDtoXKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_BtoDtoXKL_m, "Bplus", "SIGNAL", PDFtype::BRBtoDtoXKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_BtoDtoXKL_m, "Bzero", "SIGNAL", PDFtype::BRBtoDtoXKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_BtoDtoXKL_m, "Bzero", "SIGNAL", PDFtype::BRBtoDtoXKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_BtoDtoXKL_m, "Bzero", "SIGNAL", PDFtype::BRBtoDtoXKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

    GetPDFs(MC_dirname_CHG, "root", CHG_BtoDtoXKL_m, "Bplus", "CHG", PDFtype::BRBtoDtoXKLDOWN, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise", 0);
    GetPDFs(MC_dirname_MIX, "root", MIX_BtoDtoXKL_m, "Bzero", "MIX", PDFtype::BRBtoDtoXKLDOWN, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise", 0);

    // BR(B -> X KL KL)
    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs1_BRBtoXKLKL_p, "Bplus", "SIGNAL", PDFtype::BRXKLKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs1_BRBtoXKLKL_p, "Bplus", "SIGNAL", PDFtype::BRXKLKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs1_BRBtoXKLKL_p, "Bplus", "SIGNAL", PDFtype::BRXKLKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs1_BRBtoXKLKL_p, "Bzero", "SIGNAL", PDFtype::BRXKLKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs1_BRBtoXKLKL_p, "Bzero", "SIGNAL", PDFtype::BRXKLKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs1_BRBtoXKLKL_p, "Bzero", "SIGNAL", PDFtype::BRXKLKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 1);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs2_BRBtoXKLKL_p, "Bplus", "SIGNAL", PDFtype::BRXKLKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs2_BRBtoXKLKL_p, "Bplus", "SIGNAL", PDFtype::BRXKLKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs2_BRBtoXKLKL_p, "Bplus", "SIGNAL", PDFtype::BRXKLKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs2_BRBtoXKLKL_p, "Bzero", "SIGNAL", PDFtype::BRXKLKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs2_BRBtoXKLKL_p, "Bzero", "SIGNAL", PDFtype::BRXKLKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs2_BRBtoXKLKL_p, "Bzero", "SIGNAL", PDFtype::BRXKLKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 2);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_BRBtoXKLKL_p, "Bplus", "SIGNAL", PDFtype::BRXKLKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_BRBtoXKLKL_p, "Bplus", "SIGNAL", PDFtype::BRXKLKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_BRBtoXKLKL_p, "Bplus", "SIGNAL", PDFtype::BRXKLKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_BRBtoXKLKL_p, "Bzero", "SIGNAL", PDFtype::BRXKLKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_BRBtoXKLKL_p, "Bzero", "SIGNAL", PDFtype::BRXKLKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_BRBtoXKLKL_p, "Bzero", "SIGNAL", PDFtype::BRXKLKLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

    GetPDFs(MC_dirname_CHG, "root", CHG_BRBtoXKLKL_p, "Bplus", "CHG", PDFtype::BRXKLKLUP, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise", 0);
    GetPDFs(MC_dirname_MIX, "root", MIX_BRBtoXKLKL_p, "Bzero", "MIX", PDFtype::BRXKLKLUP, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise", 0);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs1_BRBtoXKLKL_m, "Bplus", "SIGNAL", PDFtype::BRXKLKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs1_BRBtoXKLKL_m, "Bplus", "SIGNAL", PDFtype::BRXKLKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs1_BRBtoXKLKL_m, "Bplus", "SIGNAL", PDFtype::BRXKLKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs1_BRBtoXKLKL_m, "Bzero", "SIGNAL", PDFtype::BRXKLKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs1_BRBtoXKLKL_m, "Bzero", "SIGNAL", PDFtype::BRXKLKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs1_BRBtoXKLKL_m, "Bzero", "SIGNAL", PDFtype::BRXKLKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 1);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs2_BRBtoXKLKL_m, "Bplus", "SIGNAL", PDFtype::BRXKLKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs2_BRBtoXKLKL_m, "Bplus", "SIGNAL", PDFtype::BRXKLKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs2_BRBtoXKLKL_m, "Bplus", "SIGNAL", PDFtype::BRXKLKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs2_BRBtoXKLKL_m, "Bzero", "SIGNAL", PDFtype::BRXKLKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs2_BRBtoXKLKL_m, "Bzero", "SIGNAL", PDFtype::BRXKLKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs2_BRBtoXKLKL_m, "Bzero", "SIGNAL", PDFtype::BRXKLKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 2);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_BRBtoXKLKL_m, "Bplus", "SIGNAL", PDFtype::BRXKLKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_BRBtoXKLKL_m, "Bplus", "SIGNAL", PDFtype::BRXKLKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_BRBtoXKLKL_m, "Bplus", "SIGNAL", PDFtype::BRXKLKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_BRBtoXKLKL_m, "Bzero", "SIGNAL", PDFtype::BRXKLKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_BRBtoXKLKL_m, "Bzero", "SIGNAL", PDFtype::BRXKLKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_BRBtoXKLKL_m, "Bzero", "SIGNAL", PDFtype::BRXKLKLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

    GetPDFs(MC_dirname_CHG, "root", CHG_BRBtoXKLKL_m, "Bplus", "CHG", PDFtype::BRXKLKLDOWN, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise", 0);
    GetPDFs(MC_dirname_MIX, "root", MIX_BRBtoXKLKL_m, "Bzero", "MIX", PDFtype::BRXKLKLDOWN, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise", 0);

    // efficiency for ECL cluster from KL0
    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs1_EffECLKL_p, "Bplus", "SIGNAL", PDFtype::EffKLECLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs1_EffECLKL_p, "Bplus", "SIGNAL", PDFtype::EffKLECLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs1_EffECLKL_p, "Bplus", "SIGNAL", PDFtype::EffKLECLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs1_EffECLKL_p, "Bzero", "SIGNAL", PDFtype::EffKLECLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs1_EffECLKL_p, "Bzero", "SIGNAL", PDFtype::EffKLECLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs1_EffECLKL_p, "Bzero", "SIGNAL", PDFtype::EffKLECLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 1);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs2_EffECLKL_p, "Bplus", "SIGNAL", PDFtype::EffKLECLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs2_EffECLKL_p, "Bplus", "SIGNAL", PDFtype::EffKLECLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs2_EffECLKL_p, "Bplus", "SIGNAL", PDFtype::EffKLECLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs2_EffECLKL_p, "Bzero", "SIGNAL", PDFtype::EffKLECLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs2_EffECLKL_p, "Bzero", "SIGNAL", PDFtype::EffKLECLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs2_EffECLKL_p, "Bzero", "SIGNAL", PDFtype::EffKLECLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 2);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_EffECLKL_p, "Bplus", "SIGNAL", PDFtype::EffKLECLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_EffECLKL_p, "Bplus", "SIGNAL", PDFtype::EffKLECLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_EffECLKL_p, "Bplus", "SIGNAL", PDFtype::EffKLECLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_EffECLKL_p, "Bzero", "SIGNAL", PDFtype::EffKLECLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_EffECLKL_p, "Bzero", "SIGNAL", PDFtype::EffKLECLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_EffECLKL_p, "Bzero", "SIGNAL", PDFtype::EffKLECLUP, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

    GetPDFs(MC_dirname_CHG, "root", CHG_EffECLKL_p, "Bplus", "CHG", PDFtype::EffKLECLUP, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise", 0);
    GetPDFs(MC_dirname_MIX, "root", MIX_EffECLKL_p, "Bzero", "MIX", PDFtype::EffKLECLUP, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise", 0);
    GetPDFs(MC_dirname_UUBAR, "root", UUBAR_EffECLKL_p, "Continuum", "UUBAR", PDFtype::EffKLECLUP, ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"), "otherwise", 0);
    GetPDFs(MC_dirname_DDBAR, "root", DDBAR_EffECLKL_p, "Continuum", "DDBAR", PDFtype::EffKLECLUP, ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"), "otherwise", 0);
    GetPDFs(MC_dirname_SSBAR, "root", SSBAR_EffECLKL_p, "Continuum", "SSBAR", PDFtype::EffKLECLUP, ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"), "otherwise", 0);
    GetPDFs(MC_dirname_CHARM, "root", CHARM_EffECLKL_p, "Continuum", "CHARM", PDFtype::EffKLECLUP, ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"), "otherwise", 0);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs1_EffECLKL_m, "Bplus", "SIGNAL", PDFtype::EffKLECLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs1_EffECLKL_m, "Bplus", "SIGNAL", PDFtype::EffKLECLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs1_EffECLKL_m, "Bplus", "SIGNAL", PDFtype::EffKLECLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs1_EffECLKL_m, "Bzero", "SIGNAL", PDFtype::EffKLECLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs1_EffECLKL_m, "Bzero", "SIGNAL", PDFtype::EffKLECLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs1_EffECLKL_m, "Bzero", "SIGNAL", PDFtype::EffKLECLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 1);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs2_EffECLKL_m, "Bplus", "SIGNAL", PDFtype::EffKLECLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs2_EffECLKL_m, "Bplus", "SIGNAL", PDFtype::EffKLECLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs2_EffECLKL_m, "Bplus", "SIGNAL", PDFtype::EffKLECLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs2_EffECLKL_m, "Bzero", "SIGNAL", PDFtype::EffKLECLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs2_EffECLKL_m, "Bzero", "SIGNAL", PDFtype::EffKLECLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs2_EffECLKL_m, "Bzero", "SIGNAL", PDFtype::EffKLECLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 2);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_EffECLKL_m, "Bplus", "SIGNAL", PDFtype::EffKLECLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_EffECLKL_m, "Bplus", "SIGNAL", PDFtype::EffKLECLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_EffECLKL_m, "Bplus", "SIGNAL", PDFtype::EffKLECLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_EffECLKL_m, "Bzero", "SIGNAL", PDFtype::EffKLECLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_EffECLKL_m, "Bzero", "SIGNAL", PDFtype::EffKLECLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_EffECLKL_m, "Bzero", "SIGNAL", PDFtype::EffKLECLDOWN, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

    GetPDFs(MC_dirname_CHG, "root", CHG_EffECLKL_m, "Bplus", "CHG", PDFtype::EffKLECLDOWN, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise", 0);
    GetPDFs(MC_dirname_MIX, "root", MIX_EffECLKL_m, "Bzero", "MIX", PDFtype::EffKLECLDOWN, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise", 0);
    GetPDFs(MC_dirname_UUBAR, "root", UUBAR_EffECLKL_m, "Continuum", "UUBAR", PDFtype::EffKLECLDOWN, ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"), "otherwise", 0);
    GetPDFs(MC_dirname_DDBAR, "root", DDBAR_EffECLKL_m, "Continuum", "DDBAR", PDFtype::EffKLECLDOWN, ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"), "otherwise", 0);
    GetPDFs(MC_dirname_SSBAR, "root", SSBAR_EffECLKL_m, "Continuum", "SSBAR", PDFtype::EffKLECLDOWN, ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"), "otherwise", 0);
    GetPDFs(MC_dirname_CHARM, "root", CHARM_EffECLKL_m, "Continuum", "CHARM", PDFtype::EffKLECLDOWN, ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"), "otherwise", 0);

    // New FEI CAL
    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs1_NEWFEICAL_p, "Bplus", "SIGNAL", PDFtype::ControlFEICAL, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs1_NEWFEICAL_p, "Bplus", "SIGNAL", PDFtype::ControlFEICAL, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs1_NEWFEICAL_p, "Bplus", "SIGNAL", PDFtype::ControlFEICAL, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs1_NEWFEICAL_p, "Bzero", "SIGNAL", PDFtype::ControlFEICAL, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs1_NEWFEICAL_p, "Bzero", "SIGNAL", PDFtype::ControlFEICAL, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 1);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs1_NEWFEICAL_p, "Bzero", "SIGNAL", PDFtype::ControlFEICAL, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 1);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs2_NEWFEICAL_p, "Bplus", "SIGNAL", PDFtype::ControlFEICAL, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs2_NEWFEICAL_p, "Bplus", "SIGNAL", PDFtype::ControlFEICAL, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs2_NEWFEICAL_p, "Bplus", "SIGNAL", PDFtype::ControlFEICAL, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs2_NEWFEICAL_p, "Bzero", "SIGNAL", PDFtype::ControlFEICAL, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs2_NEWFEICAL_p, "Bzero", "SIGNAL", PDFtype::ControlFEICAL, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 2);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs2_NEWFEICAL_p, "Bzero", "SIGNAL", PDFtype::ControlFEICAL, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 2);

    GetPDFs(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_NEWFEICAL_p, "Bplus", "SIGNAL", PDFtype::ControlFEICAL, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_NEWFEICAL_p, "Bplus", "SIGNAL", PDFtype::ControlFEICAL, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_NEWFEICAL_p, "Bplus", "SIGNAL", PDFtype::ControlFEICAL, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_NEWFEICAL_p, "Bzero", "SIGNAL", PDFtype::ControlFEICAL, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_NEWFEICAL_p, "Bzero", "SIGNAL", PDFtype::ControlFEICAL, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_NEWFEICAL_p, "Bzero", "SIGNAL", PDFtype::ControlFEICAL, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

    GetNegativeChangePDFs(Signal_MXs1_nominal, Signal_MXs1_NEWFEICAL_p, Signal_MXs1_NEWFEICAL_m);
    GetNegativeChangePDFs(Signal_MXs2_nominal, Signal_MXs2_NEWFEICAL_p, Signal_MXs2_NEWFEICAL_m);
    GetNegativeChangePDFs(Signal_MXs3_nominal, Signal_MXs3_NEWFEICAL_p, Signal_MXs3_NEWFEICAL_m);

    // calculate all uncorrelated pdfs
    ClearHist(Signal_MXs1_all_uncorrelated);
    ClearHist(Signal_MXs2_all_uncorrelated);
    ClearHist(Signal_MXs3_all_uncorrelated);
    ClearHist(CHG_all_uncorrelated);
    ClearHist(MIX_all_uncorrelated);
    ClearHist(UUBAR_all_uncorrelated);
    ClearHist(DDBAR_all_uncorrelated);
    ClearHist(SSBAR_all_uncorrelated);
    ClearHist(CHARM_all_uncorrelated);

    AddSQRTHist(Signal_MXs1_all_uncorrelated, Signal_MXs1_FEI_uncorrelated, RarityBins);
    AddSQRTHist(Signal_MXs2_all_uncorrelated, Signal_MXs2_FEI_uncorrelated, RarityBins);
    AddSQRTHist(Signal_MXs3_all_uncorrelated, Signal_MXs3_FEI_uncorrelated, RarityBins);
    AddSQRTHist(CHG_all_uncorrelated, CHG_FEI_uncorrelated, RarityBins);
    AddSQRTHist(MIX_all_uncorrelated, MIX_FEI_uncorrelated, RarityBins);

    AddSQRTHist(Signal_MXs1_all_uncorrelated, Signal_MXs1_KID_uncorrelated, RarityBins);
    AddSQRTHist(Signal_MXs2_all_uncorrelated, Signal_MXs2_KID_uncorrelated, RarityBins);
    AddSQRTHist(Signal_MXs3_all_uncorrelated, Signal_MXs3_KID_uncorrelated, RarityBins);
    AddSQRTHist(CHG_all_uncorrelated, CHG_KID_uncorrelated, RarityBins);
    AddSQRTHist(MIX_all_uncorrelated, MIX_KID_uncorrelated, RarityBins);
    AddSQRTHist(UUBAR_all_uncorrelated, UUBAR_KID_uncorrelated, RarityBins);
    AddSQRTHist(DDBAR_all_uncorrelated, DDBAR_KID_uncorrelated, RarityBins);
    AddSQRTHist(SSBAR_all_uncorrelated, SSBAR_KID_uncorrelated, RarityBins);
    AddSQRTHist(CHARM_all_uncorrelated, CHARM_KID_uncorrelated, RarityBins);

    AddSQRTHist(Signal_MXs1_all_uncorrelated, Signal_MXs1_PID_uncorrelated, RarityBins);
    AddSQRTHist(Signal_MXs2_all_uncorrelated, Signal_MXs2_PID_uncorrelated, RarityBins);
    AddSQRTHist(Signal_MXs3_all_uncorrelated, Signal_MXs3_PID_uncorrelated, RarityBins);
    AddSQRTHist(CHG_all_uncorrelated, CHG_PID_uncorrelated, RarityBins);
    AddSQRTHist(MIX_all_uncorrelated, MIX_PID_uncorrelated, RarityBins);
    AddSQRTHist(UUBAR_all_uncorrelated, UUBAR_PID_uncorrelated, RarityBins);
    AddSQRTHist(DDBAR_all_uncorrelated, DDBAR_PID_uncorrelated, RarityBins);
    AddSQRTHist(SSBAR_all_uncorrelated, SSBAR_PID_uncorrelated, RarityBins);
    AddSQRTHist(CHARM_all_uncorrelated, CHARM_PID_uncorrelated, RarityBins);

    AddSQRTHist(Signal_MXs1_all_uncorrelated, Signal_MXs1_pi0_uncorrelated, RarityBins);
    AddSQRTHist(Signal_MXs2_all_uncorrelated, Signal_MXs2_pi0_uncorrelated, RarityBins);
    AddSQRTHist(Signal_MXs3_all_uncorrelated, Signal_MXs3_pi0_uncorrelated, RarityBins);
    AddSQRTHist(CHG_all_uncorrelated, CHG_pi0_uncorrelated, RarityBins);
    AddSQRTHist(MIX_all_uncorrelated, MIX_pi0_uncorrelated, RarityBins);
    AddSQRTHist(UUBAR_all_uncorrelated, UUBAR_pi0_uncorrelated, RarityBins);
    AddSQRTHist(DDBAR_all_uncorrelated, DDBAR_pi0_uncorrelated, RarityBins);
    AddSQRTHist(SSBAR_all_uncorrelated, SSBAR_pi0_uncorrelated, RarityBins);
    AddSQRTHist(CHARM_all_uncorrelated, CHARM_pi0_uncorrelated, RarityBins);

    AddSQRTHist(Signal_MXs1_all_uncorrelated, Signal_MXs1_multiplicity_uncorrelated, RarityBins);
    AddSQRTHist(Signal_MXs2_all_uncorrelated, Signal_MXs2_multiplicity_uncorrelated, RarityBins);
    AddSQRTHist(Signal_MXs3_all_uncorrelated, Signal_MXs3_multiplicity_uncorrelated, RarityBins);
    AddSQRTHist(CHG_all_uncorrelated, CHG_multiplicity_uncorrelated, RarityBins);
    AddSQRTHist(MIX_all_uncorrelated, MIX_multiplicity_uncorrelated, RarityBins);
    AddSQRTHist(UUBAR_all_uncorrelated, UUBAR_multiplicity_uncorrelated, RarityBins);
    AddSQRTHist(DDBAR_all_uncorrelated, DDBAR_multiplicity_uncorrelated, RarityBins);
    AddSQRTHist(SSBAR_all_uncorrelated, SSBAR_multiplicity_uncorrelated, RarityBins);
    AddSQRTHist(CHARM_all_uncorrelated, CHARM_multiplicity_uncorrelated, RarityBins);

    // read additional relative uncertainty 
    GetRelativeError(relative_uncertainty_file, CHG_rel_uncer, RarityBins);
    GetRelativeError(relative_uncertainty_file, MIX_rel_uncer, RarityBins);
    GetRelativeError(relative_uncertainty_file, UUBAR_rel_uncer, RarityBins);
    GetRelativeError(relative_uncertainty_file, DDBAR_rel_uncer, RarityBins);
    GetRelativeError(relative_uncertainty_file, SSBAR_rel_uncer, RarityBins);
    GetRelativeError(relative_uncertainty_file, CHARM_rel_uncer, RarityBins);

    // calculate MC statistical uncertainties (relative errors)
    GetMCstatisticalRelativeError(Signal_MXs1_nominal, Signal_MXs1_MC_stat, RarityBins);
    GetMCstatisticalRelativeError(Signal_MXs2_nominal, Signal_MXs2_MC_stat, RarityBins);
    GetMCstatisticalRelativeError(Signal_MXs3_nominal, Signal_MXs3_MC_stat, RarityBins);
    GetMCstatisticalRelativeError(CHG_nominal, CHG_MC_stat, RarityBins);
    GetMCstatisticalRelativeError(MIX_nominal, MIX_MC_stat, RarityBins);
    GetMCstatisticalRelativeError(UUBAR_nominal, UUBAR_MC_stat, RarityBins);
    GetMCstatisticalRelativeError(DDBAR_nominal, DDBAR_MC_stat, RarityBins);
    GetMCstatisticalRelativeError(SSBAR_nominal, SSBAR_MC_stat, RarityBins);
    GetMCstatisticalRelativeError(CHARM_nominal, CHARM_MC_stat, RarityBins);

    // all of uncorrelated uncertainties + MC statistical uncertainties
    AddSQRTHist(Signal_MXs1_all_uncorrelated_MC_stat, Signal_MXs1_all_uncorrelated, Signal_MXs1_MC_stat, RarityBins);
    AddSQRTHist(Signal_MXs2_all_uncorrelated_MC_stat, Signal_MXs2_all_uncorrelated, Signal_MXs2_MC_stat, RarityBins);
    AddSQRTHist(Signal_MXs3_all_uncorrelated_MC_stat, Signal_MXs3_all_uncorrelated, Signal_MXs3_MC_stat, RarityBins);
    AddSQRTHist(CHG_all_uncorrelated_MC_stat, CHG_all_uncorrelated, CHG_MC_stat, RarityBins);
    AddSQRTHist(MIX_all_uncorrelated_MC_stat, MIX_all_uncorrelated, MIX_MC_stat, RarityBins);
    AddSQRTHist(UUBAR_all_uncorrelated_MC_stat, UUBAR_all_uncorrelated, UUBAR_MC_stat, RarityBins);
    AddSQRTHist(DDBAR_all_uncorrelated_MC_stat, DDBAR_all_uncorrelated, DDBAR_MC_stat, RarityBins);
    AddSQRTHist(SSBAR_all_uncorrelated_MC_stat, SSBAR_all_uncorrelated, SSBAR_MC_stat, RarityBins);
    AddSQRTHist(CHARM_all_uncorrelated_MC_stat, CHARM_all_uncorrelated, CHARM_MC_stat, RarityBins);

    // calculate the entire relative uncertainty
    CalculateEntireRelativeUncertainty(entire_all_uncorrelated_MC_stat,
        Signal_MXs1_nominal, Signal_MXs1_all_uncorrelated_MC_stat,
        Signal_MXs2_nominal, Signal_MXs2_all_uncorrelated_MC_stat,
        Signal_MXs3_nominal, Signal_MXs3_all_uncorrelated_MC_stat,
        CHG_nominal, CHG_all_uncorrelated_MC_stat,
        MIX_nominal, MIX_all_uncorrelated_MC_stat,
        UUBAR_nominal, UUBAR_all_uncorrelated_MC_stat,
        DDBAR_nominal, DDBAR_all_uncorrelated_MC_stat,
        SSBAR_nominal, SSBAR_all_uncorrelated_MC_stat,
        CHARM_nominal, CHARM_all_uncorrelated_MC_stat, RarityBins);
    /* ====================================== */



    /* ====================================== */
    // define path for Data
    const char* DATA_dirname = "/home/jwpark/storage/BKG_gbasf2/Izayoi_again/SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Knunu";
    /* ====================================== */



    /* ====================================== */
    // Define histograms for Data
    TH1D* total_DATA = new TH1D("total_DATA", "total_DATA", RarityBins, BinMIN, BinMAX);
    /* ====================================== */



    /* ====================================== */
    // Save data
    GetDataPDF(DATA_dirname, "root", total_DATA, 1.0);
    /* ====================================== */



    /* ====================================== */
    // calculate total Signal PDFs
    AddPDFsWithRelativeUncertainty(Signal_nominal, Signal_MXs1_nominal, Signal_all_uncorrelated_MC_stat, Signal_MXs1_all_uncorrelated_MC_stat, RarityBins);
    AddPDFsWithRelativeUncertainty(Signal_nominal, Signal_MXs2_nominal, Signal_all_uncorrelated_MC_stat, Signal_MXs2_all_uncorrelated_MC_stat, RarityBins);
    AddPDFsWithRelativeUncertainty(Signal_nominal, Signal_MXs3_nominal, Signal_all_uncorrelated_MC_stat, Signal_MXs3_all_uncorrelated_MC_stat, RarityBins);

    AddPDFs(Signal_track_p, Signal_MXs1_track_p);
    AddPDFs(Signal_track_p, Signal_MXs2_track_p);
    AddPDFs(Signal_track_p, Signal_MXs3_track_p);
    AddPDFs(Signal_track_m, Signal_MXs1_track_m);
    AddPDFs(Signal_track_m, Signal_MXs2_track_m);
    AddPDFs(Signal_track_m, Signal_MXs3_track_m);

    AddPDFs(Signal_KS0_p, Signal_MXs1_KS0_p);
    AddPDFs(Signal_KS0_p, Signal_MXs2_KS0_p);
    AddPDFs(Signal_KS0_p, Signal_MXs3_KS0_p);
    AddPDFs(Signal_KS0_m, Signal_MXs1_KS0_m);
    AddPDFs(Signal_KS0_m, Signal_MXs2_KS0_m);
    AddPDFs(Signal_KS0_m, Signal_MXs3_KS0_m);

    Signal_FEI_correlated = (TH1D**)malloc(sizeof(TH1D*) * NPDFs_FEI * 2);
    for (int i = 0; i < NPDFs_FEI; i++) Signal_FEI_correlated[i] = new TH1D(("Signal_FEI_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_FEI_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
    for (int i = NPDFs_FEI; i < 2 * NPDFs_FEI; i++) Signal_FEI_correlated[i] = new TH1D(("Signal_FEI_correlated" + std::to_string(i - NPDFs_FEI) + "_m").c_str(), ("Signal_FEI_correlated" + std::to_string(i - NPDFs_FEI) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
    for (int i = 0; i < 2 * NPDFs_FEI; i++) {
        AddPDFs(Signal_FEI_correlated[i], Signal_MXs1_FEI_correlated[i]);
        AddPDFs(Signal_FEI_correlated[i], Signal_MXs2_FEI_correlated[i]);
        AddPDFs(Signal_FEI_correlated[i], Signal_MXs3_FEI_correlated[i]);
    }

    Signal_KID_correlated = (TH1D**)malloc(sizeof(TH1D*) * NPDFs_KID * 2);
    for (int i = 0; i < NPDFs_KID; i++) Signal_KID_correlated[i] = new TH1D(("Signal_KID_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_KID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
    for (int i = NPDFs_KID; i < 2 * NPDFs_KID; i++) Signal_KID_correlated[i] = new TH1D(("Signal_KID_correlated" + std::to_string(i - NPDFs_KID) + "_m").c_str(), ("Signal_KID_correlated" + std::to_string(i - NPDFs_KID) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
    for (int i = 0; i < 2 * NPDFs_KID; i++) {
        AddPDFs(Signal_KID_correlated[i], Signal_MXs1_KID_correlated[i]);
        AddPDFs(Signal_KID_correlated[i], Signal_MXs2_KID_correlated[i]);
        AddPDFs(Signal_KID_correlated[i], Signal_MXs3_KID_correlated[i]);
    }

    Signal_PID_correlated = (TH1D**)malloc(sizeof(TH1D*) * NPDFs_PID * 2);
    for (int i = 0; i < NPDFs_PID; i++) Signal_PID_correlated[i] = new TH1D(("Signal_PID_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_PID_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
    for (int i = NPDFs_PID; i < 2 * NPDFs_PID; i++) Signal_PID_correlated[i] = new TH1D(("Signal_PID_correlated" + std::to_string(i - NPDFs_PID) + "_m").c_str(), ("Signal_PID_correlated" + std::to_string(i - NPDFs_PID) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
    for (int i = 0; i < 2 * NPDFs_PID; i++) {
        AddPDFs(Signal_PID_correlated[i], Signal_MXs1_PID_correlated[i]);
        AddPDFs(Signal_PID_correlated[i], Signal_MXs2_PID_correlated[i]);
        AddPDFs(Signal_PID_correlated[i], Signal_MXs3_PID_correlated[i]);
    }

    for (int dmIndex = 0; dmIndex < 2 * NBRdmID(); dmIndex++) {
        AddPDFs(Signal_BRs.at(dmIndex), Signal_MXs1_BRs.at(dmIndex));
        AddPDFs(Signal_BRs.at(dmIndex), Signal_MXs2_BRs.at(dmIndex));
        AddPDFs(Signal_BRs.at(dmIndex), Signal_MXs3_BRs.at(dmIndex));
    }

    Signal_pi0_correlated = (TH1D**)malloc(sizeof(TH1D*) * NPDFs_pi0 * 2);
    for (int i = 0; i < NPDFs_pi0; i++) Signal_pi0_correlated[i] = new TH1D(("Signal_pi0_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_pi0_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
    for (int i = NPDFs_pi0; i < 2 * NPDFs_pi0; i++) Signal_pi0_correlated[i] = new TH1D(("Signal_pi0_correlated" + std::to_string(i - NPDFs_pi0) + "_m").c_str(), ("Signal_pi0_correlated" + std::to_string(i - NPDFs_pi0) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
    for (int i = 0; i < 2 * NPDFs_pi0; i++) {
        AddPDFs(Signal_pi0_correlated[i], Signal_MXs1_pi0_correlated[i]);
        AddPDFs(Signal_pi0_correlated[i], Signal_MXs2_pi0_correlated[i]);
        AddPDFs(Signal_pi0_correlated[i], Signal_MXs3_pi0_correlated[i]);
    }

    AddPDFs(Signal_Kff1_p, Signal_MXs1_Kff1_p);
    AddPDFs(Signal_Kff1_p, Signal_MXs2_nominal);
    AddPDFs(Signal_Kff1_p, Signal_MXs3_Kff1_p);
    AddPDFs(Signal_Kff1_m, Signal_MXs1_Kff1_m);
    AddPDFs(Signal_Kff1_m, Signal_MXs2_nominal);
    AddPDFs(Signal_Kff1_m, Signal_MXs3_Kff1_m);
    AddPDFs(Signal_Kff2_p, Signal_MXs1_Kff2_p);
    AddPDFs(Signal_Kff2_p, Signal_MXs2_nominal);
    AddPDFs(Signal_Kff2_p, Signal_MXs3_Kff2_p);
    AddPDFs(Signal_Kff2_m, Signal_MXs1_Kff2_m);
    AddPDFs(Signal_Kff2_m, Signal_MXs2_nominal);
    AddPDFs(Signal_Kff2_m, Signal_MXs3_Kff2_m);
    AddPDFs(Signal_Kff3_p, Signal_MXs1_Kff3_p);
    AddPDFs(Signal_Kff3_p, Signal_MXs2_nominal);
    AddPDFs(Signal_Kff3_p, Signal_MXs3_Kff3_p);
    AddPDFs(Signal_Kff3_m, Signal_MXs1_Kff3_m);
    AddPDFs(Signal_Kff3_m, Signal_MXs2_nominal);
    AddPDFs(Signal_Kff3_m, Signal_MXs3_Kff3_m);

    AddPDFs(Signal_Kstarff1_p, Signal_MXs1_nominal);
    AddPDFs(Signal_Kstarff1_p, Signal_MXs2_Kstarff1_p);
    AddPDFs(Signal_Kstarff1_p, Signal_MXs3_Kstarff1_p);
    AddPDFs(Signal_Kstarff1_m, Signal_MXs1_nominal);
    AddPDFs(Signal_Kstarff1_m, Signal_MXs2_Kstarff1_m);
    AddPDFs(Signal_Kstarff1_m, Signal_MXs3_Kstarff1_m);
    AddPDFs(Signal_Kstarff2_p, Signal_MXs1_nominal);
    AddPDFs(Signal_Kstarff2_p, Signal_MXs2_Kstarff2_p);
    AddPDFs(Signal_Kstarff2_p, Signal_MXs3_Kstarff2_p);
    AddPDFs(Signal_Kstarff2_m, Signal_MXs1_nominal);
    AddPDFs(Signal_Kstarff2_m, Signal_MXs2_Kstarff2_m);
    AddPDFs(Signal_Kstarff2_m, Signal_MXs3_Kstarff2_m);
    AddPDFs(Signal_Kstarff3_p, Signal_MXs1_nominal);
    AddPDFs(Signal_Kstarff3_p, Signal_MXs2_Kstarff3_p);
    AddPDFs(Signal_Kstarff3_p, Signal_MXs3_Kstarff3_p);
    AddPDFs(Signal_Kstarff3_m, Signal_MXs1_nominal);
    AddPDFs(Signal_Kstarff3_m, Signal_MXs2_Kstarff3_m);
    AddPDFs(Signal_Kstarff3_m, Signal_MXs3_Kstarff3_m);
    AddPDFs(Signal_Kstarff4_p, Signal_MXs1_nominal);
    AddPDFs(Signal_Kstarff4_p, Signal_MXs2_Kstarff4_p);
    AddPDFs(Signal_Kstarff4_p, Signal_MXs3_Kstarff4_p);
    AddPDFs(Signal_Kstarff4_m, Signal_MXs1_nominal);
    AddPDFs(Signal_Kstarff4_m, Signal_MXs2_Kstarff4_m);
    AddPDFs(Signal_Kstarff4_m, Signal_MXs3_Kstarff4_m);
    AddPDFs(Signal_Kstarff5_p, Signal_MXs1_nominal);
    AddPDFs(Signal_Kstarff5_p, Signal_MXs2_Kstarff5_p);
    AddPDFs(Signal_Kstarff5_p, Signal_MXs3_Kstarff5_p);
    AddPDFs(Signal_Kstarff5_m, Signal_MXs1_nominal);
    AddPDFs(Signal_Kstarff5_m, Signal_MXs2_Kstarff5_m);
    AddPDFs(Signal_Kstarff5_m, Signal_MXs3_Kstarff5_m);
    AddPDFs(Signal_Kstarff6_p, Signal_MXs1_nominal);
    AddPDFs(Signal_Kstarff6_p, Signal_MXs2_Kstarff6_p);
    AddPDFs(Signal_Kstarff6_p, Signal_MXs3_Kstarff6_p);
    AddPDFs(Signal_Kstarff6_m, Signal_MXs1_nominal);
    AddPDFs(Signal_Kstarff6_m, Signal_MXs2_Kstarff6_m);
    AddPDFs(Signal_Kstarff6_m, Signal_MXs3_Kstarff6_m);
    AddPDFs(Signal_Kstarff7_p, Signal_MXs1_nominal);
    AddPDFs(Signal_Kstarff7_p, Signal_MXs2_Kstarff7_p);
    AddPDFs(Signal_Kstarff7_p, Signal_MXs3_Kstarff7_p);
    AddPDFs(Signal_Kstarff7_m, Signal_MXs1_nominal);
    AddPDFs(Signal_Kstarff7_m, Signal_MXs2_Kstarff7_m);
    AddPDFs(Signal_Kstarff7_m, Signal_MXs3_Kstarff7_m);
    AddPDFs(Signal_Kstarff8_p, Signal_MXs1_nominal);
    AddPDFs(Signal_Kstarff8_p, Signal_MXs2_Kstarff8_p);
    AddPDFs(Signal_Kstarff8_p, Signal_MXs3_Kstarff8_p);
    AddPDFs(Signal_Kstarff8_m, Signal_MXs1_nominal);
    AddPDFs(Signal_Kstarff8_m, Signal_MXs2_Kstarff8_m);
    AddPDFs(Signal_Kstarff8_m, Signal_MXs3_Kstarff8_m);
    AddPDFs(Signal_Kstarff9_p, Signal_MXs1_nominal);
    AddPDFs(Signal_Kstarff9_p, Signal_MXs2_Kstarff9_p);
    AddPDFs(Signal_Kstarff9_p, Signal_MXs3_Kstarff9_p);
    AddPDFs(Signal_Kstarff9_m, Signal_MXs1_nominal);
    AddPDFs(Signal_Kstarff9_m, Signal_MXs2_Kstarff9_m);
    AddPDFs(Signal_Kstarff9_m, Signal_MXs3_Kstarff9_m);

    AddPDFs(Signal_Kfrac_p, Signal_MXs1_Kfrac_p);
    AddPDFs(Signal_Kfrac_p, Signal_MXs2_nominal);
    AddPDFs(Signal_Kfrac_p, Signal_MXs3_Kfrac_p);
    AddPDFs(Signal_Kfrac_m, Signal_MXs1_Kfrac_m);
    AddPDFs(Signal_Kfrac_m, Signal_MXs2_nominal);
    AddPDFs(Signal_Kfrac_m, Signal_MXs3_Kfrac_m);

    AddPDFs(Signal_Kstarfrac_p, Signal_MXs1_nominal);
    AddPDFs(Signal_Kstarfrac_p, Signal_MXs2_Kstarfrac_p);
    AddPDFs(Signal_Kstarfrac_p, Signal_MXs3_Kstarfrac_p);
    AddPDFs(Signal_Kstarfrac_m, Signal_MXs1_nominal);
    AddPDFs(Signal_Kstarfrac_m, Signal_MXs2_Kstarfrac_m);
    AddPDFs(Signal_Kstarfrac_m, Signal_MXs3_Kstarfrac_m);

    for (int sign = 0; sign < 2; sign++) {
        for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
            for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
                int temp_index = MxsBin * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma) + Category + sign * corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma) * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma);

                AddPDFs(Signal_Fragmentaions.at(temp_index), Signal_MXs1_Fragmentaions.at(temp_index));
                AddPDFs(Signal_Fragmentaions.at(temp_index), Signal_MXs2_Fragmentaions.at(temp_index));
                AddPDFs(Signal_Fragmentaions.at(temp_index), Signal_MXs3_Fragmentaions.at(temp_index));

            }
        }
    }

    AddPDFs(Signal_pf_p, Signal_MXs1_nominal);
    AddPDFs(Signal_pf_p, Signal_MXs2_nominal);
    AddPDFs(Signal_pf_p, Signal_MXs3_pf_p);
    AddPDFs(Signal_pf_m, Signal_MXs1_nominal);
    AddPDFs(Signal_pf_m, Signal_MXs2_nominal);
    AddPDFs(Signal_pf_m, Signal_MXs3_pf_m);

    AddPDFs(Signal_mb_p, Signal_MXs1_nominal);
    AddPDFs(Signal_mb_p, Signal_MXs2_nominal);
    AddPDFs(Signal_mb_p, Signal_MXs3_mb_p);
    AddPDFs(Signal_mb_m, Signal_MXs1_nominal);
    AddPDFs(Signal_mb_m, Signal_MXs2_nominal);
    AddPDFs(Signal_mb_m, Signal_MXs3_mb_m);

    AddPDFs(Signal_transition_p, Signal_MXs1_nominal);
    AddPDFs(Signal_transition_p, Signal_MXs2_nominal);
    AddPDFs(Signal_transition_p, Signal_MXs3_transition_p);
    AddPDFs(Signal_transition_m, Signal_MXs1_nominal);
    AddPDFs(Signal_transition_m, Signal_MXs2_nominal);
    AddPDFs(Signal_transition_m, Signal_MXs3_transition_m);

    AddPDFs(Signal_mKstar_p, Signal_MXs1_nominal);
    AddPDFs(Signal_mKstar_p, Signal_MXs2_mKstar_p);
    AddPDFs(Signal_mKstar_p, Signal_MXs3_mKstar_p);
    AddPDFs(Signal_mKstar_m, Signal_MXs1_nominal);
    AddPDFs(Signal_mKstar_m, Signal_MXs2_mKstar_m);
    AddPDFs(Signal_mKstar_m, Signal_MXs3_mKstar_m);

    AddPDFs(Signal_Xnn_p, Signal_MXs1_Xnn_p);
    AddPDFs(Signal_Xnn_p, Signal_MXs2_Xnn_p);
    AddPDFs(Signal_Xnn_p, Signal_MXs3_Xnn_p);
    AddPDFs(Signal_Xnn_m, Signal_MXs1_Xnn_m);
    AddPDFs(Signal_Xnn_m, Signal_MXs2_Xnn_m);
    AddPDFs(Signal_Xnn_m, Signal_MXs3_Xnn_m);

    Signal_multiplicity_correlated = (TH1D**)malloc(sizeof(TH1D*) * NPDFs_multiplicity * 2);
    for (int i = 0; i < NPDFs_multiplicity; i++) Signal_multiplicity_correlated[i] = new TH1D(("Signal_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), ("Signal_multiplicity_correlated" + std::to_string(i) + "_p").c_str(), RarityBins, BinMIN, BinMAX);
    for (int i = NPDFs_multiplicity; i < 2 * NPDFs_multiplicity; i++) Signal_multiplicity_correlated[i] = new TH1D(("Signal_multiplicity_correlated" + std::to_string(i - NPDFs_multiplicity) + "_m").c_str(), ("Signal_multiplicity_correlated" + std::to_string(i - NPDFs_multiplicity) + "_m").c_str(), RarityBins, BinMIN, BinMAX);
    for (int i = 0; i < 2 * NPDFs_multiplicity; i++) {
        AddPDFs(Signal_multiplicity_correlated[i], Signal_MXs1_multiplicity_correlated[i]);
        AddPDFs(Signal_multiplicity_correlated[i], Signal_MXs2_multiplicity_correlated[i]);
        AddPDFs(Signal_multiplicity_correlated[i], Signal_MXs3_multiplicity_correlated[i]);
    }

    AddPDFs(Signal_BtoDtoXKL_p, Signal_MXs1_BtoDtoXKL_p);
    AddPDFs(Signal_BtoDtoXKL_p, Signal_MXs2_BtoDtoXKL_p);
    AddPDFs(Signal_BtoDtoXKL_p, Signal_MXs3_BtoDtoXKL_p);
    AddPDFs(Signal_BtoDtoXKL_m, Signal_MXs1_BtoDtoXKL_m);
    AddPDFs(Signal_BtoDtoXKL_m, Signal_MXs2_BtoDtoXKL_m);
    AddPDFs(Signal_BtoDtoXKL_m, Signal_MXs3_BtoDtoXKL_m);

    AddPDFs(Signal_BRBtoXKLKL_p, Signal_MXs1_BRBtoXKLKL_p);
    AddPDFs(Signal_BRBtoXKLKL_p, Signal_MXs2_BRBtoXKLKL_p);
    AddPDFs(Signal_BRBtoXKLKL_p, Signal_MXs3_BRBtoXKLKL_p);
    AddPDFs(Signal_BRBtoXKLKL_m, Signal_MXs1_BRBtoXKLKL_m);
    AddPDFs(Signal_BRBtoXKLKL_m, Signal_MXs2_BRBtoXKLKL_m);
    AddPDFs(Signal_BRBtoXKLKL_m, Signal_MXs3_BRBtoXKLKL_m);

    AddPDFs(Signal_EffECLKL_p, Signal_MXs1_EffECLKL_p);
    AddPDFs(Signal_EffECLKL_p, Signal_MXs2_EffECLKL_p);
    AddPDFs(Signal_EffECLKL_p, Signal_MXs3_EffECLKL_p);
    AddPDFs(Signal_EffECLKL_m, Signal_MXs1_EffECLKL_m);
    AddPDFs(Signal_EffECLKL_m, Signal_MXs2_EffECLKL_m);
    AddPDFs(Signal_EffECLKL_m, Signal_MXs3_EffECLKL_m);

    AddPDFs(Signal_NEWFEICAL_p, Signal_MXs1_NEWFEICAL_p);
    AddPDFs(Signal_NEWFEICAL_p, Signal_MXs2_NEWFEICAL_p);
    AddPDFs(Signal_NEWFEICAL_p, Signal_MXs3_NEWFEICAL_p);
    AddPDFs(Signal_NEWFEICAL_m, Signal_MXs1_NEWFEICAL_m);
    AddPDFs(Signal_NEWFEICAL_m, Signal_MXs2_NEWFEICAL_m);
    AddPDFs(Signal_NEWFEICAL_m, Signal_MXs3_NEWFEICAL_m);
    /* ====================================== */



    /* ====================================== */
    // Select specific MXs bin
    const int MXsBin = 1;

    SaveSpecificMXsBin(Signal_nominal, MXsBin);
    SaveSpecificMXsBin(Signal_MXs1_nominal, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_nominal, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_nominal, MXsBin);
    SaveSpecificMXsBin(CHG_nominal, MXsBin);
    SaveSpecificMXsBin(MIX_nominal, MXsBin);
    SaveSpecificMXsBin(UUBAR_nominal, MXsBin);
    SaveSpecificMXsBin(DDBAR_nominal, MXsBin);
    SaveSpecificMXsBin(SSBAR_nominal, MXsBin);
    SaveSpecificMXsBin(CHARM_nominal, MXsBin);

    // track uncertainty
    SaveSpecificMXsBin(Signal_track_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs1_track_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_track_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_track_p, MXsBin);
    SaveSpecificMXsBin(CHG_track_p, MXsBin);
    SaveSpecificMXsBin(MIX_track_p, MXsBin);
    SaveSpecificMXsBin(UUBAR_track_p, MXsBin);
    SaveSpecificMXsBin(DDBAR_track_p, MXsBin);
    SaveSpecificMXsBin(SSBAR_track_p, MXsBin);
    SaveSpecificMXsBin(CHARM_track_p, MXsBin);
    SaveSpecificMXsBin(Signal_track_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs1_track_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_track_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_track_m, MXsBin);
    SaveSpecificMXsBin(CHG_track_m, MXsBin);
    SaveSpecificMXsBin(MIX_track_m, MXsBin);
    SaveSpecificMXsBin(UUBAR_track_m, MXsBin);
    SaveSpecificMXsBin(DDBAR_track_m, MXsBin);
    SaveSpecificMXsBin(SSBAR_track_m, MXsBin);
    SaveSpecificMXsBin(CHARM_track_m, MXsBin);

    // KS0 uncertainty
    SaveSpecificMXsBin(Signal_KS0_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs1_KS0_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_KS0_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_KS0_p, MXsBin);
    SaveSpecificMXsBin(CHG_KS0_p, MXsBin);
    SaveSpecificMXsBin(MIX_KS0_p, MXsBin);
    SaveSpecificMXsBin(UUBAR_KS0_p, MXsBin);
    SaveSpecificMXsBin(DDBAR_KS0_p, MXsBin);
    SaveSpecificMXsBin(SSBAR_KS0_p, MXsBin);
    SaveSpecificMXsBin(CHARM_KS0_p, MXsBin);
    SaveSpecificMXsBin(Signal_KS0_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs1_KS0_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_KS0_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_KS0_m, MXsBin);
    SaveSpecificMXsBin(CHG_KS0_m, MXsBin);
    SaveSpecificMXsBin(MIX_KS0_m, MXsBin);
    SaveSpecificMXsBin(UUBAR_KS0_m, MXsBin);
    SaveSpecificMXsBin(DDBAR_KS0_m, MXsBin);
    SaveSpecificMXsBin(SSBAR_KS0_m, MXsBin);
    SaveSpecificMXsBin(CHARM_KS0_m, MXsBin);

    // FEI uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_FEI; i++) {
        SaveSpecificMXsBin(Signal_FEI_correlated[i], MXsBin);
        SaveSpecificMXsBin(Signal_MXs1_FEI_correlated[i], MXsBin);
        SaveSpecificMXsBin(Signal_MXs2_FEI_correlated[i], MXsBin);
        SaveSpecificMXsBin(Signal_MXs3_FEI_correlated[i], MXsBin);
        SaveSpecificMXsBin(CHG_FEI_correlated[i], MXsBin);
        SaveSpecificMXsBin(MIX_FEI_correlated[i], MXsBin);
    }

    // FEI uncertainty (uncorrelated)
    SaveSpecificMXsBin(Signal_MXs1_FEI_uncorrelated, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_FEI_uncorrelated, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_FEI_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHG_FEI_uncorrelated, MXsBin);
    SaveSpecificMXsBin(MIX_FEI_uncorrelated, MXsBin);

    // Kaon PID uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_KID; i++) {
        SaveSpecificMXsBin(Signal_KID_correlated[i], MXsBin);
        SaveSpecificMXsBin(Signal_MXs1_KID_correlated[i], MXsBin);
        SaveSpecificMXsBin(Signal_MXs2_KID_correlated[i], MXsBin);
        SaveSpecificMXsBin(Signal_MXs3_KID_correlated[i], MXsBin);
        SaveSpecificMXsBin(CHG_KID_correlated[i], MXsBin);
        SaveSpecificMXsBin(MIX_KID_correlated[i], MXsBin);
        SaveSpecificMXsBin(UUBAR_KID_correlated[i], MXsBin);
        SaveSpecificMXsBin(DDBAR_KID_correlated[i], MXsBin);
        SaveSpecificMXsBin(SSBAR_KID_correlated[i], MXsBin);
        SaveSpecificMXsBin(CHARM_KID_correlated[i], MXsBin);
    }

    // Kaon PID uncertainty (uncorrelated)
    SaveSpecificMXsBin(Signal_MXs1_KID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_KID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_KID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHG_KID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(MIX_KID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(UUBAR_KID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(DDBAR_KID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(SSBAR_KID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHARM_KID_uncorrelated, MXsBin);

    // Pion PID uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_PID; i++) {
        SaveSpecificMXsBin(Signal_PID_correlated[i], MXsBin);
        SaveSpecificMXsBin(Signal_MXs1_PID_correlated[i], MXsBin);
        SaveSpecificMXsBin(Signal_MXs2_PID_correlated[i], MXsBin);
        SaveSpecificMXsBin(Signal_MXs3_PID_correlated[i], MXsBin);
        SaveSpecificMXsBin(CHG_PID_correlated[i], MXsBin);
        SaveSpecificMXsBin(MIX_PID_correlated[i], MXsBin);
        SaveSpecificMXsBin(UUBAR_PID_correlated[i], MXsBin);
        SaveSpecificMXsBin(DDBAR_PID_correlated[i], MXsBin);
        SaveSpecificMXsBin(SSBAR_PID_correlated[i], MXsBin);
        SaveSpecificMXsBin(CHARM_PID_correlated[i], MXsBin);
    }

    // Pion PID uncertainty (uncorrelated)
    SaveSpecificMXsBin(Signal_MXs1_PID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_PID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_PID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHG_PID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(MIX_PID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(UUBAR_PID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(DDBAR_PID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(SSBAR_PID_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHARM_PID_uncorrelated, MXsBin);

    // BB BR uncertainty
    for (int dmIndex = 0; dmIndex < 2 * NBRdmID(); dmIndex++) {
        SaveSpecificMXsBin(Signal_BRs.at(dmIndex), MXsBin);
        SaveSpecificMXsBin(Signal_MXs1_BRs.at(dmIndex), MXsBin);
        SaveSpecificMXsBin(Signal_MXs2_BRs.at(dmIndex), MXsBin);
        SaveSpecificMXsBin(Signal_MXs3_BRs.at(dmIndex), MXsBin);
        SaveSpecificMXsBin(CHG_BRs.at(dmIndex), MXsBin);
        SaveSpecificMXsBin(MIX_BRs.at(dmIndex), MXsBin);
    }

    // pi0 uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_pi0; i++) {
        SaveSpecificMXsBin(Signal_pi0_correlated[i], MXsBin);
        SaveSpecificMXsBin(Signal_MXs1_pi0_correlated[i], MXsBin);
        SaveSpecificMXsBin(Signal_MXs2_pi0_correlated[i], MXsBin);
        SaveSpecificMXsBin(Signal_MXs3_pi0_correlated[i], MXsBin);
        SaveSpecificMXsBin(CHG_pi0_correlated[i], MXsBin);
        SaveSpecificMXsBin(MIX_pi0_correlated[i], MXsBin);
        SaveSpecificMXsBin(UUBAR_pi0_correlated[i], MXsBin);
        SaveSpecificMXsBin(DDBAR_pi0_correlated[i], MXsBin);
        SaveSpecificMXsBin(SSBAR_pi0_correlated[i], MXsBin);
        SaveSpecificMXsBin(CHARM_pi0_correlated[i], MXsBin);
    }

    // pi0 uncertainty (uncorrelated)
    SaveSpecificMXsBin(Signal_MXs1_pi0_uncorrelated, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_pi0_uncorrelated, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_pi0_uncorrelated, MXsBin);
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
    SaveSpecificMXsBin(Signal_MXs1_Kff1_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs1_Kff1_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs1_Kff2_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs1_Kff2_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs1_Kff3_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs1_Kff3_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kff1_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kff1_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kff2_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kff2_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kff3_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kff3_m, MXsBin);

    // Kstar nu nubar form factor
    SaveSpecificMXsBin(Signal_Kstarff1_p, MXsBin);
    SaveSpecificMXsBin(Signal_Kstarff1_m, MXsBin);
    SaveSpecificMXsBin(Signal_Kstarff2_p, MXsBin);
    SaveSpecificMXsBin(Signal_Kstarff2_m, MXsBin);
    SaveSpecificMXsBin(Signal_Kstarff3_p, MXsBin);
    SaveSpecificMXsBin(Signal_Kstarff3_m, MXsBin);
    SaveSpecificMXsBin(Signal_Kstarff4_p, MXsBin);
    SaveSpecificMXsBin(Signal_Kstarff4_m, MXsBin);
    SaveSpecificMXsBin(Signal_Kstarff5_p, MXsBin);
    SaveSpecificMXsBin(Signal_Kstarff5_m, MXsBin);
    SaveSpecificMXsBin(Signal_Kstarff6_p, MXsBin);
    SaveSpecificMXsBin(Signal_Kstarff6_m, MXsBin);
    SaveSpecificMXsBin(Signal_Kstarff7_p, MXsBin);
    SaveSpecificMXsBin(Signal_Kstarff7_m, MXsBin);
    SaveSpecificMXsBin(Signal_Kstarff8_p, MXsBin);
    SaveSpecificMXsBin(Signal_Kstarff8_m, MXsBin);
    SaveSpecificMXsBin(Signal_Kstarff9_p, MXsBin);
    SaveSpecificMXsBin(Signal_Kstarff9_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_Kstarff1_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_Kstarff1_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_Kstarff2_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_Kstarff2_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_Kstarff3_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_Kstarff3_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_Kstarff4_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_Kstarff4_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_Kstarff5_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_Kstarff5_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_Kstarff6_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_Kstarff6_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_Kstarff7_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_Kstarff7_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_Kstarff8_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_Kstarff8_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_Kstarff9_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_Kstarff9_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kstarff1_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kstarff1_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kstarff2_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kstarff2_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kstarff3_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kstarff3_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kstarff4_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kstarff4_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kstarff5_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kstarff5_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kstarff6_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kstarff6_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kstarff7_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kstarff7_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kstarff8_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kstarff8_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kstarff9_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kstarff9_m, MXsBin);

    // faction
    SaveSpecificMXsBin(Signal_Kfrac_p, MXsBin);
    SaveSpecificMXsBin(Signal_Kfrac_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs1_Kfrac_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs1_Kfrac_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kfrac_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kfrac_m, MXsBin);
    SaveSpecificMXsBin(Signal_Kstarfrac_p, MXsBin);
    SaveSpecificMXsBin(Signal_Kstarfrac_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_Kstarfrac_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_Kstarfrac_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kstarfrac_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Kstarfrac_m, MXsBin);

    // fragmentation uncertainty
    for (int sign = 0; sign < 2; sign++) {
        for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
            for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
                int temp_index = MxsBin * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma) + Category + sign * corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma) * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma);

                SaveSpecificMXsBin(Signal_Fragmentaions.at(temp_index), MXsBin);
                SaveSpecificMXsBin(Signal_MXs1_Fragmentaions.at(temp_index), MXsBin);
                SaveSpecificMXsBin(Signal_MXs2_Fragmentaions.at(temp_index), MXsBin);
                SaveSpecificMXsBin(Signal_MXs3_Fragmentaions.at(temp_index), MXsBin);
                SaveSpecificMXsBin(CHG_Fragmentaions.at(temp_index), MXsBin);
                SaveSpecificMXsBin(MIX_Fragmentaions.at(temp_index), MXsBin);

            }
        }
    }

    // pf uncertainty
    SaveSpecificMXsBin(Signal_pf_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_pf_p, MXsBin);
    SaveSpecificMXsBin(Signal_pf_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_pf_m, MXsBin);

    // mb uncertainty
    SaveSpecificMXsBin(Signal_mb_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_mb_p, MXsBin);
    SaveSpecificMXsBin(Signal_mb_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_mb_m, MXsBin);

    // transition uncertainty
    SaveSpecificMXsBin(Signal_transition_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_transition_p, MXsBin);
    SaveSpecificMXsBin(Signal_transition_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_transition_m, MXsBin);

    // fixed mKstar uncertainty
    SaveSpecificMXsBin(Signal_mKstar_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_mKstar_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_mKstar_p, MXsBin);
    SaveSpecificMXsBin(Signal_mKstar_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_mKstar_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_mKstar_m, MXsBin);

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
    SaveSpecificMXsBin(Signal_MXs1_Xnn_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_Xnn_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Xnn_p, MXsBin);
    SaveSpecificMXsBin(CHG_Xnn_p, MXsBin);
    SaveSpecificMXsBin(MIX_Xnn_p, MXsBin);
    SaveSpecificMXsBin(Signal_Xnn_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs1_Xnn_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_Xnn_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_Xnn_m, MXsBin);
    SaveSpecificMXsBin(CHG_Xnn_m, MXsBin);
    SaveSpecificMXsBin(MIX_Xnn_m, MXsBin);

    // multiplicity uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_multiplicity; i++) {
        SaveSpecificMXsBin(Signal_multiplicity_correlated[i], MXsBin);
        SaveSpecificMXsBin(Signal_MXs1_multiplicity_correlated[i], MXsBin);
        SaveSpecificMXsBin(Signal_MXs2_multiplicity_correlated[i], MXsBin);
        SaveSpecificMXsBin(Signal_MXs3_multiplicity_correlated[i], MXsBin);
        SaveSpecificMXsBin(CHG_multiplicity_correlated[i], MXsBin);
        SaveSpecificMXsBin(MIX_multiplicity_correlated[i], MXsBin);
        SaveSpecificMXsBin(UUBAR_multiplicity_correlated[i], MXsBin);
        SaveSpecificMXsBin(DDBAR_multiplicity_correlated[i], MXsBin);
        SaveSpecificMXsBin(SSBAR_multiplicity_correlated[i], MXsBin);
        SaveSpecificMXsBin(CHARM_multiplicity_correlated[i], MXsBin);
    }

    // multiplicity uncertainty (uncorrelated)
    SaveSpecificMXsBin(Signal_MXs1_multiplicity_uncorrelated, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_multiplicity_uncorrelated, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_multiplicity_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHG_multiplicity_uncorrelated, MXsBin);
    SaveSpecificMXsBin(MIX_multiplicity_uncorrelated, MXsBin);
    SaveSpecificMXsBin(UUBAR_multiplicity_uncorrelated, MXsBin);
    SaveSpecificMXsBin(DDBAR_multiplicity_uncorrelated, MXsBin);
    SaveSpecificMXsBin(SSBAR_multiplicity_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHARM_multiplicity_uncorrelated, MXsBin);

    // B -> [D -> X KL0] anything uncertainties
    SaveSpecificMXsBin(Signal_BtoDtoXKL_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs1_BtoDtoXKL_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_BtoDtoXKL_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_BtoDtoXKL_p, MXsBin);
    SaveSpecificMXsBin(CHG_BtoDtoXKL_p, MXsBin);
    SaveSpecificMXsBin(MIX_BtoDtoXKL_p, MXsBin);
    SaveSpecificMXsBin(Signal_BtoDtoXKL_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs1_BtoDtoXKL_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_BtoDtoXKL_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_BtoDtoXKL_m, MXsBin);
    SaveSpecificMXsBin(CHG_BtoDtoXKL_m, MXsBin);
    SaveSpecificMXsBin(MIX_BtoDtoXKL_m, MXsBin);

    // BR(B -> X KL KL)
    SaveSpecificMXsBin(Signal_BRBtoXKLKL_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs1_BRBtoXKLKL_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_BRBtoXKLKL_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_BRBtoXKLKL_p, MXsBin);
    SaveSpecificMXsBin(CHG_BRBtoXKLKL_p, MXsBin);
    SaveSpecificMXsBin(MIX_BRBtoXKLKL_p, MXsBin);
    SaveSpecificMXsBin(Signal_BRBtoXKLKL_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs1_BRBtoXKLKL_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_BRBtoXKLKL_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_BRBtoXKLKL_m, MXsBin);
    SaveSpecificMXsBin(CHG_BRBtoXKLKL_m, MXsBin);
    SaveSpecificMXsBin(MIX_BRBtoXKLKL_m, MXsBin);

    // efficiency for ECL cluster from KL0
    SaveSpecificMXsBin(Signal_EffECLKL_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs1_EffECLKL_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_EffECLKL_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_EffECLKL_p, MXsBin);
    SaveSpecificMXsBin(CHG_EffECLKL_p, MXsBin);
    SaveSpecificMXsBin(MIX_EffECLKL_p, MXsBin);
    SaveSpecificMXsBin(UUBAR_EffECLKL_p, MXsBin);
    SaveSpecificMXsBin(DDBAR_EffECLKL_p, MXsBin);
    SaveSpecificMXsBin(SSBAR_EffECLKL_p, MXsBin);
    SaveSpecificMXsBin(CHARM_EffECLKL_p, MXsBin);
    SaveSpecificMXsBin(Signal_EffECLKL_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs1_EffECLKL_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_EffECLKL_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_EffECLKL_m, MXsBin);
    SaveSpecificMXsBin(CHG_EffECLKL_m, MXsBin);
    SaveSpecificMXsBin(MIX_EffECLKL_m, MXsBin);
    SaveSpecificMXsBin(UUBAR_EffECLKL_m, MXsBin);
    SaveSpecificMXsBin(DDBAR_EffECLKL_m, MXsBin);
    SaveSpecificMXsBin(SSBAR_EffECLKL_m, MXsBin);
    SaveSpecificMXsBin(CHARM_EffECLKL_m, MXsBin);

    // New FEI CAL
    SaveSpecificMXsBin(Signal_NEWFEICAL_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs1_NEWFEICAL_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_NEWFEICAL_p, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_NEWFEICAL_p, MXsBin);
    SaveSpecificMXsBin(Signal_NEWFEICAL_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs1_NEWFEICAL_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_NEWFEICAL_m, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_NEWFEICAL_m, MXsBin);

    // all uncorrelated uncertainties
    SaveSpecificMXsBin(Signal_MXs1_all_uncorrelated, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_all_uncorrelated, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_all_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHG_all_uncorrelated, MXsBin);
    SaveSpecificMXsBin(MIX_all_uncorrelated, MXsBin);
    SaveSpecificMXsBin(UUBAR_all_uncorrelated, MXsBin);
    SaveSpecificMXsBin(DDBAR_all_uncorrelated, MXsBin);
    SaveSpecificMXsBin(SSBAR_all_uncorrelated, MXsBin);
    SaveSpecificMXsBin(CHARM_all_uncorrelated, MXsBin);

    // MC statistical uncertainties
    SaveSpecificMXsBin(Signal_MXs1_MC_stat, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_MC_stat, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_MC_stat, MXsBin);
    SaveSpecificMXsBin(CHG_MC_stat, MXsBin);
    SaveSpecificMXsBin(MIX_MC_stat, MXsBin);
    SaveSpecificMXsBin(UUBAR_MC_stat, MXsBin);
    SaveSpecificMXsBin(DDBAR_MC_stat, MXsBin);
    SaveSpecificMXsBin(SSBAR_MC_stat, MXsBin);
    SaveSpecificMXsBin(CHARM_MC_stat, MXsBin);

    // additional relative uncertainty
    SaveSpecificMXsBin(CHG_rel_uncer, MXsBin);
    SaveSpecificMXsBin(MIX_rel_uncer, MXsBin);
    SaveSpecificMXsBin(UUBAR_rel_uncer, MXsBin);
    SaveSpecificMXsBin(DDBAR_rel_uncer, MXsBin);
    SaveSpecificMXsBin(SSBAR_rel_uncer, MXsBin);
    SaveSpecificMXsBin(CHARM_rel_uncer, MXsBin);

    // all of uncorrelated uncertainties + MC statistical uncertainties
    SaveSpecificMXsBin(Signal_all_uncorrelated_MC_stat, MXsBin);
    SaveSpecificMXsBin(Signal_MXs1_all_uncorrelated_MC_stat, MXsBin);
    SaveSpecificMXsBin(Signal_MXs2_all_uncorrelated_MC_stat, MXsBin);
    SaveSpecificMXsBin(Signal_MXs3_all_uncorrelated_MC_stat, MXsBin);
    SaveSpecificMXsBin(CHG_all_uncorrelated_MC_stat, MXsBin);
    SaveSpecificMXsBin(MIX_all_uncorrelated_MC_stat, MXsBin);
    SaveSpecificMXsBin(UUBAR_all_uncorrelated_MC_stat, MXsBin);
    SaveSpecificMXsBin(DDBAR_all_uncorrelated_MC_stat, MXsBin);
    SaveSpecificMXsBin(SSBAR_all_uncorrelated_MC_stat, MXsBin);
    SaveSpecificMXsBin(CHARM_all_uncorrelated_MC_stat, MXsBin);

    // the entire relative uncertainty
    SaveSpecificMXsBin(entire_all_uncorrelated_MC_stat, MXsBin);

    // data
    SaveSpecificMXsBin(total_DATA, MXsBin);
    /* ====================================== */



    /* ====================================== */
    // Save histograms
    TFile* file = new TFile("PDFandDATA.root", "RECREATE");

    Signal_nominal->Write();
    Signal_MXs1_nominal->Write();
    Signal_MXs2_nominal->Write();
    Signal_MXs3_nominal->Write();
    CHG_nominal->Write();
    MIX_nominal->Write();
    UUBAR_nominal->Write();
    DDBAR_nominal->Write();
    SSBAR_nominal->Write();
    CHARM_nominal->Write();

    // track uncertainty
    Signal_track_p->Write();
    Signal_MXs1_track_p->Write();
    Signal_MXs2_track_p->Write();
    Signal_MXs3_track_p->Write();
    CHG_track_p->Write();
    MIX_track_p->Write();
    UUBAR_track_p->Write();
    DDBAR_track_p->Write();
    SSBAR_track_p->Write();
    CHARM_track_p->Write();
    Signal_track_m->Write();
    Signal_MXs1_track_m->Write();
    Signal_MXs2_track_m->Write();
    Signal_MXs3_track_m->Write();
    CHG_track_m->Write();
    MIX_track_m->Write();
    UUBAR_track_m->Write();
    DDBAR_track_m->Write();
    SSBAR_track_m->Write();
    CHARM_track_m->Write();

    // KS0 uncertainty
    Signal_KS0_p->Write();
    Signal_MXs1_KS0_p->Write();
    Signal_MXs2_KS0_p->Write();
    Signal_MXs3_KS0_p->Write();
    CHG_KS0_p->Write();
    MIX_KS0_p->Write();
    UUBAR_KS0_p->Write();
    DDBAR_KS0_p->Write();
    SSBAR_KS0_p->Write();
    CHARM_KS0_p->Write();
    Signal_KS0_m->Write();
    Signal_MXs1_KS0_m->Write();
    Signal_MXs2_KS0_m->Write();
    Signal_MXs3_KS0_m->Write();
    CHG_KS0_m->Write();
    MIX_KS0_m->Write();
    UUBAR_KS0_m->Write();
    DDBAR_KS0_m->Write();
    SSBAR_KS0_m->Write();
    CHARM_KS0_m->Write();

    // FEI uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_FEI; i++) {
        Signal_FEI_correlated[i]->Write();
        Signal_MXs1_FEI_correlated[i]->Write();
        Signal_MXs2_FEI_correlated[i]->Write();
        Signal_MXs3_FEI_correlated[i]->Write();
        CHG_FEI_correlated[i]->Write();
        MIX_FEI_correlated[i]->Write();
    }

    // FEI uncertainty (uncorrelated)
    Signal_MXs1_FEI_uncorrelated->Write();
    Signal_MXs2_FEI_uncorrelated->Write();
    Signal_MXs3_FEI_uncorrelated->Write();
    CHG_FEI_uncorrelated->Write();
    MIX_FEI_uncorrelated->Write();

    // Kaon PID uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_KID; i++) {
        Signal_KID_correlated[i]->Write();
        Signal_MXs1_KID_correlated[i]->Write();
        Signal_MXs2_KID_correlated[i]->Write();
        Signal_MXs3_KID_correlated[i]->Write();
        CHG_KID_correlated[i]->Write();
        MIX_KID_correlated[i]->Write();
        UUBAR_KID_correlated[i]->Write();
        DDBAR_KID_correlated[i]->Write();
        SSBAR_KID_correlated[i]->Write();
        CHARM_KID_correlated[i]->Write();
    }

    // Kaon PID uncertainty (uncorrelated)
    Signal_MXs1_KID_uncorrelated->Write();
    Signal_MXs2_KID_uncorrelated->Write();
    Signal_MXs3_KID_uncorrelated->Write();
    CHG_KID_uncorrelated->Write();
    MIX_KID_uncorrelated->Write();
    UUBAR_KID_uncorrelated->Write();
    DDBAR_KID_uncorrelated->Write();
    SSBAR_KID_uncorrelated->Write();
    CHARM_KID_uncorrelated->Write();

    // Pion PID uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_PID; i++) {
        Signal_PID_correlated[i]->Write();
        Signal_MXs1_PID_correlated[i]->Write();
        Signal_MXs2_PID_correlated[i]->Write();
        Signal_MXs3_PID_correlated[i]->Write();
        CHG_PID_correlated[i]->Write();
        MIX_PID_correlated[i]->Write();
        UUBAR_PID_correlated[i]->Write();
        DDBAR_PID_correlated[i]->Write();
        SSBAR_PID_correlated[i]->Write();
        CHARM_PID_correlated[i]->Write();
    }

    // Pion PID uncertainty (uncorrelated)
    Signal_MXs1_PID_uncorrelated->Write();
    Signal_MXs2_PID_uncorrelated->Write();
    Signal_MXs3_PID_uncorrelated->Write();
    CHG_PID_uncorrelated->Write();
    MIX_PID_uncorrelated->Write();
    UUBAR_PID_uncorrelated->Write();
    DDBAR_PID_uncorrelated->Write();
    SSBAR_PID_uncorrelated->Write();
    CHARM_PID_uncorrelated->Write();

    // BB BR uncertainty
    for (int dmIndex = 0; dmIndex < 2 * NBRdmID(); dmIndex++) {
        Signal_BRs.at(dmIndex)->Write();
        Signal_MXs1_BRs.at(dmIndex)->Write();
        Signal_MXs2_BRs.at(dmIndex)->Write();
        Signal_MXs3_BRs.at(dmIndex)->Write();
        CHG_BRs.at(dmIndex)->Write();
        MIX_BRs.at(dmIndex)->Write();
    }

    // pi0 uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_pi0; i++) {
        Signal_pi0_correlated[i]->Write();
        Signal_MXs1_pi0_correlated[i]->Write();
        Signal_MXs2_pi0_correlated[i]->Write();
        Signal_MXs3_pi0_correlated[i]->Write();
        CHG_pi0_correlated[i]->Write();
        MIX_pi0_correlated[i]->Write();
        UUBAR_pi0_correlated[i]->Write();
        DDBAR_pi0_correlated[i]->Write();
        SSBAR_pi0_correlated[i]->Write();
        CHARM_pi0_correlated[i]->Write();
    }

    // pi0 uncertainty (uncorrelated)
    Signal_MXs1_pi0_uncorrelated->Write();
    Signal_MXs2_pi0_uncorrelated->Write();
    Signal_MXs3_pi0_uncorrelated->Write();
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
    Signal_MXs1_Kff1_p->Write();
    Signal_MXs1_Kff1_m->Write();
    Signal_MXs1_Kff2_p->Write();
    Signal_MXs1_Kff2_m->Write();
    Signal_MXs1_Kff3_p->Write();
    Signal_MXs1_Kff3_m->Write();
    Signal_MXs3_Kff1_p->Write();
    Signal_MXs3_Kff1_m->Write();
    Signal_MXs3_Kff2_p->Write();
    Signal_MXs3_Kff2_m->Write();
    Signal_MXs3_Kff3_p->Write();
    Signal_MXs3_Kff3_m->Write();

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
    Signal_MXs2_Kstarff1_p->Write();
    Signal_MXs2_Kstarff1_m->Write();
    Signal_MXs2_Kstarff2_p->Write();
    Signal_MXs2_Kstarff2_m->Write();
    Signal_MXs2_Kstarff3_p->Write();
    Signal_MXs2_Kstarff3_m->Write();
    Signal_MXs2_Kstarff4_p->Write();
    Signal_MXs2_Kstarff4_m->Write();
    Signal_MXs2_Kstarff5_p->Write();
    Signal_MXs2_Kstarff5_m->Write();
    Signal_MXs2_Kstarff6_p->Write();
    Signal_MXs2_Kstarff6_m->Write();
    Signal_MXs2_Kstarff7_p->Write();
    Signal_MXs2_Kstarff7_m->Write();
    Signal_MXs2_Kstarff8_p->Write();
    Signal_MXs2_Kstarff8_m->Write();
    Signal_MXs2_Kstarff9_p->Write();
    Signal_MXs2_Kstarff9_m->Write();
    Signal_MXs3_Kstarff1_p->Write();
    Signal_MXs3_Kstarff1_m->Write();
    Signal_MXs3_Kstarff2_p->Write();
    Signal_MXs3_Kstarff2_m->Write();
    Signal_MXs3_Kstarff3_p->Write();
    Signal_MXs3_Kstarff3_m->Write();
    Signal_MXs3_Kstarff4_p->Write();
    Signal_MXs3_Kstarff4_m->Write();
    Signal_MXs3_Kstarff5_p->Write();
    Signal_MXs3_Kstarff5_m->Write();
    Signal_MXs3_Kstarff6_p->Write();
    Signal_MXs3_Kstarff6_m->Write();
    Signal_MXs3_Kstarff7_p->Write();
    Signal_MXs3_Kstarff7_m->Write();
    Signal_MXs3_Kstarff8_p->Write();
    Signal_MXs3_Kstarff8_m->Write();
    Signal_MXs3_Kstarff9_p->Write();
    Signal_MXs3_Kstarff9_m->Write();

    // faction
    Signal_Kfrac_p->Write();
    Signal_Kfrac_m->Write();
    Signal_MXs1_Kfrac_p->Write();
    Signal_MXs1_Kfrac_m->Write();
    Signal_MXs3_Kfrac_p->Write();
    Signal_MXs3_Kfrac_m->Write();
    Signal_Kstarfrac_p->Write();
    Signal_Kstarfrac_m->Write();
    Signal_MXs2_Kstarfrac_p->Write();
    Signal_MXs2_Kstarfrac_m->Write();
    Signal_MXs3_Kstarfrac_p->Write();
    Signal_MXs3_Kstarfrac_m->Write();

    // fragmentation uncertainty
    for (int sign = 0; sign < 2; sign++) {
        for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
            for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
                int temp_index = MxsBin * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma) + Category + sign * corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma) * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma);

                Signal_Fragmentaions.at(temp_index)->Write();
                Signal_MXs1_Fragmentaions.at(temp_index)->Write();
                Signal_MXs2_Fragmentaions.at(temp_index)->Write();
                Signal_MXs3_Fragmentaions.at(temp_index)->Write();
                CHG_Fragmentaions.at(temp_index)->Write();
                MIX_Fragmentaions.at(temp_index)->Write();

            }
        }
    }

    // pf uncertainty
    Signal_pf_p->Write();
    Signal_MXs3_pf_p->Write();
    Signal_pf_m->Write();
    Signal_MXs3_pf_m->Write();

    // mb uncertainty
    Signal_mb_p->Write();
    Signal_MXs3_mb_p->Write();
    Signal_mb_m->Write();
    Signal_MXs3_mb_m->Write();

    // transition uncertainty
    Signal_transition_p->Write();
    Signal_MXs3_transition_p->Write();
    Signal_transition_m->Write();
    Signal_MXs3_transition_m->Write();

    // fixed mKstar uncertainty
    Signal_mKstar_p->Write();
    Signal_MXs2_mKstar_p->Write();
    Signal_MXs3_mKstar_p->Write();
    Signal_mKstar_m->Write();
    Signal_MXs2_mKstar_m->Write();
    Signal_MXs3_mKstar_m->Write();

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
    Signal_MXs1_Xnn_p->Write();
    Signal_MXs2_Xnn_p->Write();
    Signal_MXs3_Xnn_p->Write();
    CHG_Xnn_p->Write();
    MIX_Xnn_p->Write();
    Signal_Xnn_m->Write();
    Signal_MXs1_Xnn_m->Write();
    Signal_MXs2_Xnn_m->Write();
    Signal_MXs3_Xnn_m->Write();
    CHG_Xnn_m->Write();
    MIX_Xnn_m->Write();

    // multiplicity uncertainty (correlated)
    for (int i = 0; i < 2 * NPDFs_multiplicity; i++) {
        Signal_multiplicity_correlated[i]->Write();
        Signal_MXs1_multiplicity_correlated[i]->Write();
        Signal_MXs2_multiplicity_correlated[i]->Write();
        Signal_MXs3_multiplicity_correlated[i]->Write();
        CHG_multiplicity_correlated[i]->Write();
        MIX_multiplicity_correlated[i]->Write();
        UUBAR_multiplicity_correlated[i]->Write();
        DDBAR_multiplicity_correlated[i]->Write();
        SSBAR_multiplicity_correlated[i]->Write();
        CHARM_multiplicity_correlated[i]->Write();
    }

    // multiplicity uncertainty (uncorrelated)
    Signal_MXs1_multiplicity_uncorrelated->Write();
    Signal_MXs2_multiplicity_uncorrelated->Write();
    Signal_MXs3_multiplicity_uncorrelated->Write();
    CHG_multiplicity_uncorrelated->Write();
    MIX_multiplicity_uncorrelated->Write();
    UUBAR_multiplicity_uncorrelated->Write();
    DDBAR_multiplicity_uncorrelated->Write();
    SSBAR_multiplicity_uncorrelated->Write();
    CHARM_multiplicity_uncorrelated->Write();

    // B -> [D -> X KL0] anything uncertainties
    Signal_BtoDtoXKL_p->Write();
    Signal_MXs1_BtoDtoXKL_p->Write();
    Signal_MXs2_BtoDtoXKL_p->Write();
    Signal_MXs3_BtoDtoXKL_p->Write();
    CHG_BtoDtoXKL_p->Write();
    MIX_BtoDtoXKL_p->Write();
    Signal_BtoDtoXKL_m->Write();
    Signal_MXs1_BtoDtoXKL_m->Write();
    Signal_MXs2_BtoDtoXKL_m->Write();
    Signal_MXs3_BtoDtoXKL_m->Write();
    CHG_BtoDtoXKL_m->Write();
    MIX_BtoDtoXKL_m->Write();

    // BR(B -> X KL KL)
    Signal_BRBtoXKLKL_p->Write();
    Signal_MXs1_BRBtoXKLKL_p->Write();
    Signal_MXs2_BRBtoXKLKL_p->Write();
    Signal_MXs3_BRBtoXKLKL_p->Write();
    CHG_BRBtoXKLKL_p->Write();
    MIX_BRBtoXKLKL_p->Write();
    Signal_BRBtoXKLKL_m->Write();
    Signal_MXs1_BRBtoXKLKL_m->Write();
    Signal_MXs2_BRBtoXKLKL_m->Write();
    Signal_MXs3_BRBtoXKLKL_m->Write();
    CHG_BRBtoXKLKL_m->Write();
    MIX_BRBtoXKLKL_m->Write();

    // efficiency for ECL cluster from KL0
    Signal_EffECLKL_p->Write();
    Signal_MXs1_EffECLKL_p->Write();
    Signal_MXs2_EffECLKL_p->Write();
    Signal_MXs3_EffECLKL_p->Write();
    CHG_EffECLKL_p->Write();
    MIX_EffECLKL_p->Write();
    UUBAR_EffECLKL_p->Write();
    DDBAR_EffECLKL_p->Write();
    SSBAR_EffECLKL_p->Write();
    CHARM_EffECLKL_p->Write();
    Signal_EffECLKL_m->Write();
    Signal_MXs1_EffECLKL_m->Write();
    Signal_MXs2_EffECLKL_m->Write();
    Signal_MXs3_EffECLKL_m->Write();
    CHG_EffECLKL_m->Write();
    MIX_EffECLKL_m->Write();
    UUBAR_EffECLKL_m->Write();
    DDBAR_EffECLKL_m->Write();
    SSBAR_EffECLKL_m->Write();
    CHARM_EffECLKL_m->Write();

    // New FEI CAL
    Signal_NEWFEICAL_p->Write();
    Signal_MXs1_NEWFEICAL_p->Write();
    Signal_MXs2_NEWFEICAL_p->Write();
    Signal_MXs3_NEWFEICAL_p->Write();
    Signal_NEWFEICAL_m->Write();
    Signal_MXs1_NEWFEICAL_m->Write();
    Signal_MXs2_NEWFEICAL_m->Write();
    Signal_MXs3_NEWFEICAL_m->Write();

    // all uncorrelated uncertainties
    Signal_MXs1_all_uncorrelated->Write();
    Signal_MXs2_all_uncorrelated->Write();
    Signal_MXs3_all_uncorrelated->Write();
    CHG_all_uncorrelated->Write(); 
    MIX_all_uncorrelated->Write(); 
    UUBAR_all_uncorrelated->Write(); 
    DDBAR_all_uncorrelated->Write(); 
    SSBAR_all_uncorrelated->Write(); 
    CHARM_all_uncorrelated->Write();

    // MC statistical uncertainties
    Signal_MXs1_MC_stat->Write();
    Signal_MXs2_MC_stat->Write();
    Signal_MXs3_MC_stat->Write();
    CHG_MC_stat->Write();
    MIX_MC_stat->Write();
    UUBAR_MC_stat->Write();
    DDBAR_MC_stat->Write();
    SSBAR_MC_stat->Write();
    CHARM_MC_stat->Write();

    // additional relative uncertainty
    CHG_rel_uncer->Write();
    MIX_rel_uncer->Write();
    UUBAR_rel_uncer->Write();
    DDBAR_rel_uncer->Write();
    SSBAR_rel_uncer->Write();
    CHARM_rel_uncer->Write();

    // all of uncorrelated uncertainties + MC statistical uncertainties
    Signal_all_uncorrelated_MC_stat->Write();
    Signal_MXs1_all_uncorrelated_MC_stat->Write();
    Signal_MXs2_all_uncorrelated_MC_stat->Write();
    Signal_MXs3_all_uncorrelated_MC_stat->Write();
    CHG_all_uncorrelated_MC_stat->Write();
    MIX_all_uncorrelated_MC_stat->Write();
    UUBAR_all_uncorrelated_MC_stat->Write();
    DDBAR_all_uncorrelated_MC_stat->Write();
    SSBAR_all_uncorrelated_MC_stat->Write();
    CHARM_all_uncorrelated_MC_stat->Write();

    // the entire relative uncertainty
    entire_all_uncorrelated_MC_stat->Write();

    total_DATA->Write();

    file->Close();
    /* ====================================== */

    return 0;
}
