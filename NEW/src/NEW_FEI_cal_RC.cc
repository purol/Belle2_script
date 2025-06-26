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
#include "TLegend.h"
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

const bool KnunuOnly = false;

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
    BRXKLKLDOWN
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
Corrector_KS0 corrector_KS0;
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
    double temp_N_pi0_syst_MC15ri[N_pi0_syst_MC15ri] = { 0.0 };
    double temp_N_pi0_syst_MC15rd[N_pi0_syst_MC15rd] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    double KS0_flight_distance = 0;
    double KS0_costheta = 0;
    double KS0_p = 0;

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
        for (int i_pi0 = 0; i_pi0 < N_pi0_syst_MC15ri; i_pi0++) tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npi0MC15ribin" + std::to_string(i_pi0)).c_str(), &temp_N_pi0_syst_MC15ri[i_pi0]);
        for (int i_pi0 = 0; i_pi0 < N_pi0_syst_MC15rd; i_pi0++) tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npi0MC15rdbin" + std::to_string(i_pi0)).c_str(), &temp_N_pi0_syst_MC15rd[i_pi0]);
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
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_3D_distance", &KS0_flight_distance);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_costheta", &KS0_costheta);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_p", &KS0_p);
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

            if (KnunuOnly) {
                // select B+ --> K+ nu nubar reconstruction only
                if ((std::abs(Upsilon_ID) < MyEPSILON) && (std::abs(Bsig_ID) < MyEPSILON)) {}
                else continue;
            }

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
            double Correction_KS0 = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE), temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE), temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            if (std::string(MCTYPE) == "MC15ri") for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            else if (std::string(MCTYPE) == "MC15rd") for (int i_pi0 = 0; i_pi0 < N_pi0_syst_MC15rd; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_pi0_syst_MC15rd[i_pi0]);
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
            if (std::string(MCTYPE) == "MC15ri") Correction_KS0 = Correction_KS0 * corrector_KS0.GetCorrectionFactor(KS0_flight_distance, MCTYPE);
            else if (std::string(MCTYPE) == "MC15rd") Correction_KS0 = Correction_KS0 * corrector_KS0.GetCorrectionFactor(KS0_p, KS0_costheta, KS0_flight_distance, MCTYPE);

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

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_KS0 * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            // lets care about systematics
            if (pdftype == PDFtype::nominal) {}
            else if (pdftype == PDFtype::trackUP) {
                double Ntrack = GetNtracks(Upsilon_ID, Bsig_ID);
                double track_correction = 1.0;
                if (std::string(MCTYPE) == "MC15ri") track_correction = (1 + (Ntrack * track_rel_uncertainty_MC15ri / 100.0));
                else if (std::string(MCTYPE) == "MC15rd") track_correction = (1 + (Ntrack * track_rel_uncertainty_MC15rd / 100.0));
                total_weight = total_weight * track_correction;
            }
            else if (pdftype == PDFtype::trackDOWN) {
                double Ntrack = GetNtracks(Upsilon_ID, Bsig_ID);
                double track_correction = 1.0;
                if (std::string(MCTYPE) == "MC15ri") track_correction = (1 - (Ntrack * track_rel_uncertainty_MC15ri / 100.0));
                else if (std::string(MCTYPE) == "MC15rd") track_correction = (1 - (Ntrack * track_rel_uncertainty_MC15rd / 100.0));
                total_weight = total_weight * track_correction;
            }
            else if (pdftype == PDFtype::KS0UP) {
                printf("[GetPDFs] do not use `GetPDFs` to calculate the systematic uncertainty from KS0 efficiency calbiration factor anymore.\n");
                exit(1);
            }
            else if (pdftype == PDFtype::KS0DOWN) {
                printf("[GetPDFs] do not use `GetPDFs` to calculate the systematic uncertainty from KS0 efficiency calbiration factor anymore.\n");
                exit(1);
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
                double Transition_correction = ReadWeightHist(Xsu_Htransition_weight, MVA_var);
                total_weight = total_weight * Transition_correction;
            }
            else if (pdftype == PDFtype::TransitionXsuDOWN) {
                double Transition_correction = ReadWeightHist(Xsu_Ltransition_weight, MVA_var);
                total_weight = total_weight * Transition_correction;
            }
            else if (pdftype == PDFtype::TransitionXsdUP) {
                double Transition_correction = ReadWeightHist(Xsd_Htransition_weight, MVA_var);
                total_weight = total_weight * Transition_correction;
            }
            else if (pdftype == PDFtype::TransitionXsdDOWN) {
                double Transition_correction = ReadWeightHist(Xsd_Ltransition_weight, MVA_var);
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

double GetPDFs_NEW_FEI(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, PDFtype pdftype, double weight_var = 1.0, std::string CorrectionType = "otherwise", int true_MXs_region = 0) { // get nominal PDF with appropriate correction
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
    double temp_N_pi0_syst_MC15ri[N_pi0_syst_MC15ri] = { 0.0 };
    double temp_N_pi0_syst_MC15rd[N_pi0_syst_MC15rd] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    double KS0_flight_distance = 0;
    double KS0_costheta = 0;
    double KS0_p = 0;

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
        for (int i_pi0 = 0; i_pi0 < N_pi0_syst_MC15ri; i_pi0++) tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npi0MC15ribin" + std::to_string(i_pi0)).c_str(), &temp_N_pi0_syst_MC15ri[i_pi0]);
        for (int i_pi0 = 0; i_pi0 < N_pi0_syst_MC15rd; i_pi0++) tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npi0MC15rdbin" + std::to_string(i_pi0)).c_str(), &temp_N_pi0_syst_MC15rd[i_pi0]);
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
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_3D_distance", &KS0_flight_distance);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_costheta", &KS0_costheta);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_p", &KS0_p);
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

            if (KnunuOnly) {
                // select B+ --> K+ nu nubar reconstruction only
                if ((std::abs(Upsilon_ID) < MyEPSILON) && (std::abs(Bsig_ID) < MyEPSILON)) {}
                else continue;
            }

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
            if (strcmp(type, "Bplus") == 0) {
                if (std::abs(Btag_isSignal - 1.0) < MyEPSILON) {
                    if (std::abs(Btag_ID - 0.0) < MyEPSILON) Correction_FEI = 0.871866;
                    else if (std::abs(Btag_ID - 1.0) < MyEPSILON) Correction_FEI = 0.871866;
                    else if (std::abs(Btag_ID - 3.0) < MyEPSILON) Correction_FEI = 0.871866;
                    else if (std::abs(Btag_ID - 4.0) < MyEPSILON) Correction_FEI = 0.871866;
                    else if (std::abs(Btag_ID - 15.0) < MyEPSILON) Correction_FEI = 0.580409;
                    else if (std::abs(Btag_ID - 16.0) < MyEPSILON) Correction_FEI = 0.580409;
                    else if (std::abs(Btag_ID - 18.0) < MyEPSILON) Correction_FEI = 0.580409;
                    else if (std::abs(Btag_ID - 19.0) < MyEPSILON) Correction_FEI = 0.580409;
                    else if (std::abs(Btag_ID - 23.0) < MyEPSILON) Correction_FEI = 0.227976;
                    else if (std::abs(Btag_ID - 24.0) < MyEPSILON) Correction_FEI = 0.227976;
                    else if (std::abs(Btag_ID - 30.0) < MyEPSILON) Correction_FEI = 0.227976;
                    else Correction_FEI = 0.227976;
                }
                else Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            }
            else if (strcmp(type, "Bzero") == 0) {
                if(std::abs(Btag_isSignal - 1.0) < MyEPSILON) Correction_FEI = 0.958304; // true Btag
                else Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            }
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            double Correction_KS0 = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE), temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE), temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            if (std::string(MCTYPE) == "MC15ri") for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            else if (std::string(MCTYPE) == "MC15rd") for (int i_pi0 = 0; i_pi0 < N_pi0_syst_MC15rd; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_pi0_syst_MC15rd[i_pi0]);
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
            if (std::string(MCTYPE) == "MC15ri") Correction_KS0 = Correction_KS0 * corrector_KS0.GetCorrectionFactor(KS0_flight_distance, MCTYPE);
            else if (std::string(MCTYPE) == "MC15rd") Correction_KS0 = Correction_KS0 * corrector_KS0.GetCorrectionFactor(KS0_p, KS0_costheta, KS0_flight_distance, MCTYPE);

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

            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_KS0 * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            // lets care about systematics
            if (pdftype == PDFtype::nominal) {}
            else if (pdftype == PDFtype::trackUP) {
                double Ntrack = GetNtracks(Upsilon_ID, Bsig_ID);
                double track_correction = 1.0;
                if (std::string(MCTYPE) == "MC15ri") track_correction = (1 + (Ntrack * track_rel_uncertainty_MC15ri / 100.0));
                else if (std::string(MCTYPE) == "MC15rd") track_correction = (1 + (Ntrack * track_rel_uncertainty_MC15rd / 100.0));
                total_weight = total_weight * track_correction;
            }
            else if (pdftype == PDFtype::trackDOWN) {
                double Ntrack = GetNtracks(Upsilon_ID, Bsig_ID);
                double track_correction = 1.0;
                if (std::string(MCTYPE) == "MC15ri") track_correction = (1 - (Ntrack * track_rel_uncertainty_MC15ri / 100.0));
                else if (std::string(MCTYPE) == "MC15rd") track_correction = (1 - (Ntrack * track_rel_uncertainty_MC15rd / 100.0));
                total_weight = total_weight * track_correction;
            }
            else if (pdftype == PDFtype::KS0UP) {
                printf("[GetPDFs] do not use `GetPDFs` to calculate the systematic uncertainty from KS0 efficiency calbiration factor anymore.\n");
                exit(1);
            }
            else if (pdftype == PDFtype::KS0DOWN) {
                printf("[GetPDFs] do not use `GetPDFs` to calculate the systematic uncertainty from KS0 efficiency calbiration factor anymore.\n");
                exit(1);
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
                double Transition_correction = ReadWeightHist(Xsu_Htransition_weight, MVA_var);
                total_weight = total_weight * Transition_correction;
            }
            else if (pdftype == PDFtype::TransitionXsuDOWN) {
                double Transition_correction = ReadWeightHist(Xsu_Ltransition_weight, MVA_var);
                total_weight = total_weight * Transition_correction;
            }
            else if (pdftype == PDFtype::TransitionXsdUP) {
                double Transition_correction = ReadWeightHist(Xsd_Htransition_weight, MVA_var);
                total_weight = total_weight * Transition_correction;
            }
            else if (pdftype == PDFtype::TransitionXsdDOWN) {
                double Transition_correction = ReadWeightHist(Xsd_Ltransition_weight, MVA_var);
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
    TH1D* BKG_nominal = new TH1D("BKG_nominal", "BKG_nominal", RarityBins, BinMIN, BinMAX);

    TH1D* Signal_NEW = new TH1D("Signal_NEW", "Signal_NEW", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs1_NEW = new TH1D("Signal_MXs1_NEW", "Signal_MXs1_NEW", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs2_NEW = new TH1D("Signal_MXs2_NEW", "Signal_MXs2_NEW", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_MXs3_NEW = new TH1D("Signal_MXs3_NEW", "Signal_MXs3_NEW", RarityBins, BinMIN, BinMAX);
    TH1D* CHG_NEW = new TH1D("CHG_NEW", "CHG_NEW", RarityBins, BinMIN, BinMAX);
    TH1D* MIX_NEW = new TH1D("MIX_NEW", "MIX_NEW", RarityBins, BinMIN, BinMAX);
    TH1D* UUBAR_NEW = new TH1D("UUBAR_NEW", "UUBAR_NEW", RarityBins, BinMIN, BinMAX);
    TH1D* DDBAR_NEW = new TH1D("DDBAR_NEW", "DDBAR_NEW", RarityBins, BinMIN, BinMAX);
    TH1D* SSBAR_NEW = new TH1D("SSBAR_NEW", "SSBAR_NEW", RarityBins, BinMIN, BinMAX);
    TH1D* CHARM_NEW = new TH1D("CHARM_NEW", "CHARM_NEW", RarityBins, BinMIN, BinMAX);
    TH1D* BKG_NEW = new TH1D("BKG_NEW", "BKG_NEW", RarityBins, BinMIN, BinMAX);

    TH1D* Signal_ratio = new TH1D("Signal_ratio", ";;non-official/official", RarityBins, BinMIN, BinMAX);
    TH1D* BKG_ratio = new TH1D("BKG_ratio", ";;non-official/official", RarityBins, BinMIN, BinMAX);
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

    GetPDFs(MC_dirname_CHG, "root", CHG_nominal, "Bplus", "CHG", PDFtype::nominal, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise", 0);
    GetPDFs(MC_dirname_MIX, "root", MIX_nominal, "Bzero", "MIX", PDFtype::nominal, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise", 0);
    GetPDFs(MC_dirname_UUBAR, "root", UUBAR_nominal, "Continuum", "UUBAR", PDFtype::nominal, ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"), "otherwise", 0);
    GetPDFs(MC_dirname_DDBAR, "root", DDBAR_nominal, "Continuum", "DDBAR", PDFtype::nominal, ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"), "otherwise", 0);
    GetPDFs(MC_dirname_SSBAR, "root", SSBAR_nominal, "Continuum", "SSBAR", PDFtype::nominal, ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"), "otherwise", 0);
    GetPDFs(MC_dirname_CHARM, "root", CHARM_nominal, "Continuum", "CHARM", PDFtype::nominal, ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"), "otherwise", 0);


    GetPDFs_NEW_FEI(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs1_NEW, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 1);
    GetPDFs_NEW_FEI(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs1_NEW, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 1);
    GetPDFs_NEW_FEI(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs1_NEW, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 1);
    GetPDFs_NEW_FEI(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs1_NEW, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 1);
    GetPDFs_NEW_FEI(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs1_NEW, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 1);
    GetPDFs_NEW_FEI(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs1_NEW, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 1);

    GetPDFs_NEW_FEI(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs2_NEW, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 2);
    GetPDFs_NEW_FEI(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs2_NEW, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 2);
    GetPDFs_NEW_FEI(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs2_NEW, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 2);
    GetPDFs_NEW_FEI(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs2_NEW, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 2);
    GetPDFs_NEW_FEI(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs2_NEW, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 2);
    GetPDFs_NEW_FEI(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs2_NEW, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 2);

    GetPDFs_NEW_FEI(MC_dirname_SIGNAL, "B2Knunu", Signal_MXs3_NEW, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu", 3);
    GetPDFs_NEW_FEI(MC_dirname_SIGNAL, "B2Kstarnunu", Signal_MXs3_NEW, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise", 3);
    GetPDFs_NEW_FEI(MC_dirname_SIGNAL, "B2Xsnunu", Signal_MXs3_NEW, "Bplus", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu", 3);
    GetPDFs_NEW_FEI(MC_dirname_SIGNAL, "B02K0nunu", Signal_MXs3_NEW, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu", 3);
    GetPDFs_NEW_FEI(MC_dirname_SIGNAL, "B02Kstar0nunu", Signal_MXs3_NEW, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise", 3);
    GetPDFs_NEW_FEI(MC_dirname_SIGNAL, "B02Xsnunu", Signal_MXs3_NEW, "Bzero", "SIGNAL", PDFtype::nominal, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu", 3);

    GetPDFs_NEW_FEI(MC_dirname_CHG, "root", CHG_NEW, "Bplus", "CHG", PDFtype::nominal, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise", 0);
    GetPDFs_NEW_FEI(MC_dirname_MIX, "root", MIX_NEW, "Bzero", "MIX", PDFtype::nominal, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise", 0);
    GetPDFs_NEW_FEI(MC_dirname_UUBAR, "root", UUBAR_NEW, "Continuum", "UUBAR", PDFtype::nominal, ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"), "otherwise", 0);
    GetPDFs_NEW_FEI(MC_dirname_DDBAR, "root", DDBAR_NEW, "Continuum", "DDBAR", PDFtype::nominal, ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"), "otherwise", 0);
    GetPDFs_NEW_FEI(MC_dirname_SSBAR, "root", SSBAR_NEW, "Continuum", "SSBAR", PDFtype::nominal, ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"), "otherwise", 0);
    GetPDFs_NEW_FEI(MC_dirname_CHARM, "root", CHARM_NEW, "Continuum", "CHARM", PDFtype::nominal, ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"), "otherwise", 0);
    /* ====================================== */



    /* ====================================== */
    // calculate total PDFs
    AddPDFs(Signal_nominal, Signal_MXs1_nominal);
    AddPDFs(Signal_nominal, Signal_MXs2_nominal);
    AddPDFs(Signal_nominal, Signal_MXs3_nominal);

    AddPDFs(BKG_nominal, CHG_nominal);
    AddPDFs(BKG_nominal, MIX_nominal);
    AddPDFs(BKG_nominal, UUBAR_nominal);
    AddPDFs(BKG_nominal, DDBAR_nominal);
    AddPDFs(BKG_nominal, SSBAR_nominal);
    AddPDFs(BKG_nominal, CHARM_nominal);

    AddPDFs(Signal_NEW, Signal_MXs1_NEW);
    AddPDFs(Signal_NEW, Signal_MXs2_NEW);
    AddPDFs(Signal_NEW, Signal_MXs3_NEW);

    AddPDFs(BKG_NEW, CHG_NEW);
    AddPDFs(BKG_NEW, MIX_NEW);
    AddPDFs(BKG_NEW, UUBAR_NEW);
    AddPDFs(BKG_NEW, DDBAR_NEW);
    AddPDFs(BKG_NEW, SSBAR_NEW);
    AddPDFs(BKG_NEW, CHARM_NEW);

    Signal_ratio->Divide(Signal_NEW, Signal_nominal);
    BKG_ratio->Divide(BKG_NEW, BKG_nominal);
    /* ====================================== */



    /* ====================================== */
    // Select specific MXs bin
    const int MXsBin = 0;

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
    /* ====================================== */



    /* ====================================== */
    // Draw Hist
    Signal_ratio->SetLineColor(kBlack); Signal_ratio->SetMarkerStyle(21); Signal_ratio->Sumw2(); Signal_ratio->SetStats(0);

    TCanvas* c_temp = new TCanvas("c", "", 800, 800); c_temp->cd();

    TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.35, 1.0, 1.0);
    pad1->SetBottomMargin(0.08); pad1->SetLeftMargin(0.15);
    pad1->SetGridx(); pad1->Draw(); pad1->cd();

    gStyle->SetPalette(kPastel);

    Float_t ymax_1 = Signal_NEW->GetMaximum();
    Float_t ymax_2 = Signal_nominal->GetMaximum();
    double real_max = 0;
    if (ymax_1 > ymax_2) real_max = ymax_1;
    else real_max = ymax_2;

    Signal_NEW->SetMaximum(real_max * 1.1);

    Signal_nominal->SetStats(0);
    Signal_NEW->SetStats(0);

    Signal_nominal->SetFillColor(4); Signal_nominal->SetLineColor(4); Signal_nominal->SetLineWidth(1); Signal_nominal->SetFillStyle(3002); Signal_nominal->Draw("e1 Hist");
    Signal_NEW->SetFillColor(2); Signal_NEW->SetLineColor(2); Signal_NEW->SetLineWidth(1); Signal_NEW->SetFillStyle(3003); Signal_NEW->Draw("e1 Hist SAME");
    TLegend* legend = pad1->BuildLegend(0.95, 0.9, 0.75, 0.6);
    legend->SetFillStyle(0); legend->SetLineWidth(0);

    c_temp->cd();
    TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1, 0.3); pad2->SetBottomMargin(0.15); pad2->SetLeftMargin(0.15); pad2->SetGridx(); pad2->Draw(); pad2->cd();
    Signal_ratio->SetMinimum(0.5); Signal_ratio->SetMaximum(1.5); Signal_ratio->SetLineWidth(2);
    Signal_ratio->GetYaxis()->SetTitleSize(0.08); Signal_ratio->GetYaxis()->SetTitleOffset(0.5);
    Signal_ratio->GetXaxis()->SetLabelSize(0.08); Signal_ratio->GetYaxis()->SetLabelSize(0.08);
    Signal_ratio->Draw("e0p");
    TLine* line = new TLine(Signal_ratio->GetXaxis()->GetXmin(), 1.0, Signal_ratio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineColor(kRed);
    line->SetLineStyle(1); line->SetLineWidth(3);
    line->Draw();

    c_temp->SetBottomMargin(0.0);
    c_temp->SaveAs("FEI_cal_comparison.png");
    /* ====================================== */

    return 0;
}
