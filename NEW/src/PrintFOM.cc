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
#include "TH2.h"
#include "RooRandom.h"
#include "template.h"
#include "constants.h"
#include "ObtainWeight.h"
#include "correctors.h"
#include "base.h"
#include <limits>

# define Nstep 400
# define scan_start 0.5
# define scan_end 1.0

using std::string;
using std::to_string;

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

std::random_device rd;
std::default_random_engine generator(rd());

typedef struct FOM_structure_ {
    double total_weight;
    double MXs;
    float BDT_output;
} FOM_structure;


/* ====================================== */

void GetNominalNevt(const char* dirname, const char* included_string, const char* type, const char* sample, std::vector<FOM_structure>* Nevt_FOM, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get nominal PDF with appropriate correction
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

    if (strcmp(sample, "CHG") == 0) {}
    else if (strcmp(sample, "MIX") == 0) {}
    else if (strcmp(sample, "UUBAR") == 0) {}
    else if (strcmp(sample, "DDBAR") == 0) {}
    else if (strcmp(sample, "SSBAR") == 0) {}
    else if (strcmp(sample, "CHARM") == 0) {}
    else if (strcmp(sample, "SIGNAL") == 0) {}
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

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            if (strcmp(sample, "SIGNAL") == 0) tree_Xs->GetEntry(j);

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

            // Multiplicity correction factor, it is not applied now. it is for a systematic uncertainty
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

            /* precise but not compatible with the previous version
            double total_weight = Correction_FEI * weight_var * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
            */

            // compatible with the previous version
            double total_weight = weight_var * Correction_Knn;
            if (CorrectionType == "B2Knunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
            else if (CorrectionType == "B02K0nunu") total_weight = total_weight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
            else if (CorrectionType == "B2Xsnunu") total_weight = total_weight;
            else if (CorrectionType == "B02Xsnunu") total_weight = total_weight;


            FOM_structure temp_FOM_structure = { total_weight, Bsig_M, MVA_var };
            Nevt_FOM->push_back(temp_FOM_structure);
            Nevt = Nevt + total_weight;
        }
        input_file->Close();

        printf("%s has %lf events (with correction)\n", dirname, Nevt);

    }

    return;
}

void CalculateFOM(
    std::vector<FOM_structure> SIGNAL_FOM_, 
    std::vector<FOM_structure> CHG_FOM_, 
    std::vector<FOM_structure> MIX_FOM_, 
    std::vector<FOM_structure> UUBAR_FOM_,
    std::vector<FOM_structure> DDBAR_FOM_,
    std::vector<FOM_structure> SSBAR_FOM_,
    std::vector<FOM_structure> CHARM_FOM_,
    int mass_region,
    const char* filename) {

    double min_MXs = -1;
    double max_MXs = -1;

    if (mass_region == 0) {
        min_MXs = -1;
        max_MXs = 5.0;
    }
    else if (mass_region == 1) {
        min_MXs = 0.0;
        max_MXs = 0.6;
    }
    else if (mass_region == 2) {
        min_MXs = 0.6;
        max_MXs = 1.0;
    }
    else if (mass_region == 3) {
        min_MXs = 1.0;
        max_MXs = 2.0;
    }

    // mass cut
    for (std::vector<FOM_structure>::iterator it = SIGNAL_FOM_.begin(); it != SIGNAL_FOM_.end();)
    {
        if((min_MXs <= it->MXs) && (it->MXs < max_MXs)) ++it;
        else it = SIGNAL_FOM_.erase(it);
    }

    for (std::vector<FOM_structure>::iterator it = CHG_FOM_.begin(); it != CHG_FOM_.end();)
    {
        if ((min_MXs <= it->MXs) && (it->MXs < max_MXs)) ++it;
        else it = CHG_FOM_.erase(it);
    }

    for (std::vector<FOM_structure>::iterator it = MIX_FOM_.begin(); it != MIX_FOM_.end();)
    {
        if ((min_MXs <= it->MXs) && (it->MXs < max_MXs)) ++it;
        else it = MIX_FOM_.erase(it);
    }

    for (std::vector<FOM_structure>::iterator it = UUBAR_FOM_.begin(); it != UUBAR_FOM_.end();)
    {
        if ((min_MXs <= it->MXs) && (it->MXs < max_MXs)) ++it;
        else it = UUBAR_FOM_.erase(it);
    }

    for (std::vector<FOM_structure>::iterator it = DDBAR_FOM_.begin(); it != DDBAR_FOM_.end();)
    {
        if ((min_MXs <= it->MXs) && (it->MXs < max_MXs)) ++it;
        else it = DDBAR_FOM_.erase(it);
    }

    for (std::vector<FOM_structure>::iterator it = SSBAR_FOM_.begin(); it != SSBAR_FOM_.end();)
    {
        if ((min_MXs <= it->MXs) && (it->MXs < max_MXs)) ++it;
        else it = SSBAR_FOM_.erase(it);
    }

    for (std::vector<FOM_structure>::iterator it = CHARM_FOM_.begin(); it != CHARM_FOM_.end();)
    {
        if ((min_MXs <= it->MXs) && (it->MXs < max_MXs)) ++it;
        else it = CHARM_FOM_.erase(it);
    }

    double FOM_Matrix[Nstep];
    for (int i = 0; i < Nstep; i++) FOM_Matrix[i] = 0;

    double FBDT_cut[Nstep];
    for (int i = 0; i < Nstep; i++) FBDT_cut[i] = 0;

    for (int j = 0; j < Nstep; j++) {
        double BB_output = scan_start + (scan_end - scan_start) * j / Nstep;
        FBDT_cut[j] = BB_output;

        double BKG_num = 0;
        double SIGNAL_num = 0;

        for (int k = 0; k < SIGNAL_FOM_.size(); k++) {
            if(SIGNAL_FOM_.at(k).BDT_output > BB_output) SIGNAL_num = SIGNAL_num + SIGNAL_FOM_.at(k).total_weight;
        }

        for (int k = 0; k < CHG_FOM_.size(); k++) {
            if (CHG_FOM_.at(k).BDT_output > BB_output) BKG_num = BKG_num + CHG_FOM_.at(k).total_weight;
        }

        for (int k = 0; k < MIX_FOM_.size(); k++) {
            if (MIX_FOM_.at(k).BDT_output > BB_output) BKG_num = BKG_num + MIX_FOM_.at(k).total_weight;
        }

        for (int k = 0; k < UUBAR_FOM_.size(); k++) {
            if (UUBAR_FOM_.at(k).BDT_output > BB_output) BKG_num = BKG_num + UUBAR_FOM_.at(k).total_weight;
        }

        for (int k = 0; k < DDBAR_FOM_.size(); k++) {
            if (DDBAR_FOM_.at(k).BDT_output > BB_output) BKG_num = BKG_num + DDBAR_FOM_.at(k).total_weight;
        }

        for (int k = 0; k < SSBAR_FOM_.size(); k++) {
            if (SSBAR_FOM_.at(k).BDT_output > BB_output) BKG_num = BKG_num + SSBAR_FOM_.at(k).total_weight;
        }

        for (int k = 0; k < CHARM_FOM_.size(); k++) {
            if (CHARM_FOM_.at(k).BDT_output > BB_output) BKG_num = BKG_num + CHARM_FOM_.at(k).total_weight;
        }

        FOM_Matrix[j] = SIGNAL_num / sqrt(SIGNAL_num + BKG_num);

    }

    for (int i = 0; i < Nstep; i++) {
        printf("%f ", FOM_Matrix[i]);
        printf("\n");
    }
    
    // get MIN MAX
    double MIN_FOM = std::numeric_limits<double>::max();
    double MAX_FOM = std::numeric_limits<double>::lowest();
    for (int i = 0; i < Nstep; i++) {
        if (MIN_FOM > FOM_Matrix[i]) MIN_FOM = FOM_Matrix[i];
        if (MAX_FOM < FOM_Matrix[i]) MAX_FOM = FOM_Matrix[i];
    }

    // draw
    TCanvas* c3 = new TCanvas("c3", "Graph Draw Options", 200, 10, 600, 600);
    c3->SetLeftMargin(0.15);

    TGraph* gr3 = new TGraph(Nstep, FBDT_cut, FOM_Matrix);
    gr3->SetTitle(";FBDT cut;#frac{S}{#sqrt{S+B}}");
    gr3->SetMarkerStyle(0);
    gr3->SetMinimum(MIN_FOM * 0.7);
    gr3->Draw("");

    TLine* line = new TLine(0.86, MIN_FOM * 0.7, 0.86, MAX_FOM);
    line->SetLineColor(kRed);
    line->SetLineStyle(5);
    line->Draw();

    c3->SaveAs(filename);

    delete c3;
}

void CalculatePunziFOM(
    std::vector<FOM_structure> SIGNAL_FOM_,
    std::vector<FOM_structure> CHG_FOM_,
    std::vector<FOM_structure> MIX_FOM_,
    std::vector<FOM_structure> UUBAR_FOM_,
    std::vector<FOM_structure> DDBAR_FOM_,
    std::vector<FOM_structure> SSBAR_FOM_,
    std::vector<FOM_structure> CHARM_FOM_,
    int mass_region,
    const char* filename,
    double alpha = 1.28) {

    double min_MXs = -1;
    double max_MXs = -1;
    const double Nsig_initial = N_Kplus_nunubar_LS1 + N_Kplusstar_nunubar_LS1 + N_Xsu_nonresonant_nunubar_LS1 + N_K0_nunubar_LS1 + N_K0star_nunubar_LS1 + N_Xsd_nunubar_LS1;

    if (mass_region == 0) {
        min_MXs = -1;
        max_MXs = 5.0;
    }
    else if (mass_region == 1) {
        min_MXs = 0.0;
        max_MXs = 0.6;
    }
    else if (mass_region == 2) {
        min_MXs = 0.6;
        max_MXs = 1.0;
    }
    else if (mass_region == 3) {
        min_MXs = 1.0;
        max_MXs = 2.0;
    }

    // mass cut
    for (std::vector<FOM_structure>::iterator it = SIGNAL_FOM_.begin(); it != SIGNAL_FOM_.end();)
    {
        if ((min_MXs <= it->MXs) && (it->MXs < max_MXs)) ++it;
        else it = SIGNAL_FOM_.erase(it);
    }

    for (std::vector<FOM_structure>::iterator it = CHG_FOM_.begin(); it != CHG_FOM_.end();)
    {
        if ((min_MXs <= it->MXs) && (it->MXs < max_MXs)) ++it;
        else it = CHG_FOM_.erase(it);
    }

    for (std::vector<FOM_structure>::iterator it = MIX_FOM_.begin(); it != MIX_FOM_.end();)
    {
        if ((min_MXs <= it->MXs) && (it->MXs < max_MXs)) ++it;
        else it = MIX_FOM_.erase(it);
    }

    for (std::vector<FOM_structure>::iterator it = UUBAR_FOM_.begin(); it != UUBAR_FOM_.end();)
    {
        if ((min_MXs <= it->MXs) && (it->MXs < max_MXs)) ++it;
        else it = UUBAR_FOM_.erase(it);
    }

    for (std::vector<FOM_structure>::iterator it = DDBAR_FOM_.begin(); it != DDBAR_FOM_.end();)
    {
        if ((min_MXs <= it->MXs) && (it->MXs < max_MXs)) ++it;
        else it = DDBAR_FOM_.erase(it);
    }

    for (std::vector<FOM_structure>::iterator it = SSBAR_FOM_.begin(); it != SSBAR_FOM_.end();)
    {
        if ((min_MXs <= it->MXs) && (it->MXs < max_MXs)) ++it;
        else it = SSBAR_FOM_.erase(it);
    }

    for (std::vector<FOM_structure>::iterator it = CHARM_FOM_.begin(); it != CHARM_FOM_.end();)
    {
        if ((min_MXs <= it->MXs) && (it->MXs < max_MXs)) ++it;
        else it = CHARM_FOM_.erase(it);
    }

    double FOM_Matrix[Nstep];
    for (int i = 0; i < Nstep; i++) FOM_Matrix[i] = 0;

    double FBDT_cut[Nstep];
    for (int i = 0; i < Nstep; i++) FBDT_cut[i] = 0;

    for (int j = 0; j < Nstep; j++) {
        double BB_output = scan_start + (scan_end - scan_start) * j / Nstep;
        FBDT_cut[j] = BB_output;

        double BKG_num = 0;
        double SIGNAL_num = 0;

        for (int k = 0; k < SIGNAL_FOM_.size(); k++) {
            if (SIGNAL_FOM_.at(k).BDT_output > BB_output) SIGNAL_num = SIGNAL_num + SIGNAL_FOM_.at(k).total_weight;
        }

        for (int k = 0; k < CHG_FOM_.size(); k++) {
            if (CHG_FOM_.at(k).BDT_output > BB_output) BKG_num = BKG_num + CHG_FOM_.at(k).total_weight;
        }

        for (int k = 0; k < MIX_FOM_.size(); k++) {
            if (MIX_FOM_.at(k).BDT_output > BB_output) BKG_num = BKG_num + MIX_FOM_.at(k).total_weight;
        }

        for (int k = 0; k < UUBAR_FOM_.size(); k++) {
            if (UUBAR_FOM_.at(k).BDT_output > BB_output) BKG_num = BKG_num + UUBAR_FOM_.at(k).total_weight;
        }

        for (int k = 0; k < DDBAR_FOM_.size(); k++) {
            if (DDBAR_FOM_.at(k).BDT_output > BB_output) BKG_num = BKG_num + DDBAR_FOM_.at(k).total_weight;
        }

        for (int k = 0; k < SSBAR_FOM_.size(); k++) {
            if (SSBAR_FOM_.at(k).BDT_output > BB_output) BKG_num = BKG_num + SSBAR_FOM_.at(k).total_weight;
        }

        for (int k = 0; k < CHARM_FOM_.size(); k++) {
            if (CHARM_FOM_.at(k).BDT_output > BB_output) BKG_num = BKG_num + CHARM_FOM_.at(k).total_weight;
        }

        FOM_Matrix[j] = (SIGNAL_num / Nsig_initial) / (alpha / 2.0 + sqrt(BKG_num));

    }

    for (int i = 0; i < Nstep; i++) {
        printf("%f ", FOM_Matrix[i]);
        printf("\n");
    }

    // get MIN MAX
    double MIN_FOM = std::numeric_limits<double>::max();
    double MAX_FOM = std::numeric_limits<double>::lowest();
    for (int i = 0; i < Nstep; i++) {
        if (MIN_FOM > FOM_Matrix[i]) MIN_FOM = FOM_Matrix[i];
        if (MAX_FOM < FOM_Matrix[i]) MAX_FOM = FOM_Matrix[i];
    }

    // draw
    TCanvas* c3 = new TCanvas("c3", "Graph Draw Options", 200, 10, 600, 600);
    c3->SetLeftMargin(0.15);

    TGraph* gr3 = new TGraph(Nstep, FBDT_cut, FOM_Matrix);
    gr3->SetTitle(";FBDT cut;#frac{S}{#sqrt{S+B}}");
    gr3->SetMarkerStyle(0);
    gr3->SetMinimum(MIN_FOM * 0.7);
    gr3->Draw("");

    TLine* line = new TLine(0.86, MIN_FOM * 0.7, 0.86, MAX_FOM);
    line->SetLineColor(kRed);
    line->SetLineStyle(5);
    line->Draw();

    c3->SaveAs(filename);

    delete c3;
}

int main(int argc, char* argv[])
{

    /* ====================================== */
    // define path for Ntuple
    const char* MC_dirname_SIGNAL = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD/SIGNAL_analysis/validation_v008/final_output_root_after_MVA_Application/Merge/";

    const char* MC_dirname_CHG = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD/CHG_analysis/validation_v008/final_output_root_after_MVA_Application/Merge/";
    const char* MC_dirname_MIX = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD/MIX_analysis/validation_v008/final_output_root_after_MVA_Application/Merge/";
    const char* MC_dirname_UUBAR = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD/UUBAR_analysis/validation_v008/final_output_root_after_MVA_Application/Merge/";
    const char* MC_dirname_DDBAR = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD/DDBAR_analysis/validation_v008/final_output_root_after_MVA_Application/Merge/";
    const char* MC_dirname_SSBAR = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD/SSBAR_analysis/validation_v008/final_output_root_after_MVA_Application/Merge/";
    const char* MC_dirname_CHARM = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD/CHARM_analysis/validation_v008/final_output_root_after_MVA_Application/Merge/";
    /* ====================================== */



    /* ====================================== */
    // define std::vector
    std::vector<FOM_structure> SIGNAL_FOM;

    std::vector<FOM_structure> CHG_FOM;
    std::vector<FOM_structure> MIX_FOM;
    std::vector<FOM_structure> UUBAR_FOM;
    std::vector<FOM_structure> DDBAR_FOM;
    std::vector<FOM_structure> SSBAR_FOM;
    std::vector<FOM_structure> CHARM_FOM;
    /* ====================================== */



    /* ====================================== */
    // get nominal Nevt
    GetNominalNevt(MC_dirname_SIGNAL, "B2Knunu", "Bplus", "SIGNAL", &SIGNAL_FOM, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    GetNominalNevt(MC_dirname_SIGNAL, "B2Kstarnunu", "Bplus", "SIGNAL", &SIGNAL_FOM, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    GetNominalNevt(MC_dirname_SIGNAL, "B2Xsnunu", "Bplus", "SIGNAL", &SIGNAL_FOM, ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    GetNominalNevt(MC_dirname_SIGNAL, "B02K0nunu", "Bzero", "SIGNAL", &SIGNAL_FOM, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetNominalNevt(MC_dirname_SIGNAL, "B02Kstar0nunu", "Bzero", "SIGNAL", &SIGNAL_FOM, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    GetNominalNevt(MC_dirname_SIGNAL, "B02Xsnunu", "Bzero", "SIGNAL", &SIGNAL_FOM, ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    GetNominalNevt(MC_dirname_CHG, "root", "Bplus", "CHG", &CHG_FOM, ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    GetNominalNevt(MC_dirname_MIX, "root", "Bzero", "MIX", &MIX_FOM, ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");
    GetNominalNevt(MC_dirname_UUBAR, "root", "Continuum", "UUBAR", &UUBAR_FOM, ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"), "otherwise");
    GetNominalNevt(MC_dirname_DDBAR, "root", "Continuum", "DDBAR", &DDBAR_FOM, ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"), "otherwise");
    GetNominalNevt(MC_dirname_SSBAR, "root", "Continuum", "SSBAR", &SSBAR_FOM, ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"), "otherwise");
    GetNominalNevt(MC_dirname_CHARM, "root", "Continuum", "CHARM", &CHARM_FOM, ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"), "otherwise");
    /* ====================================== */

    CalculateFOM(SIGNAL_FOM, CHG_FOM, MIX_FOM, UUBAR_FOM, DDBAR_FOM, SSBAR_FOM, CHARM_FOM, 1, "FOM_K.png");
    CalculateFOM(SIGNAL_FOM, CHG_FOM, MIX_FOM, UUBAR_FOM, DDBAR_FOM, SSBAR_FOM, CHARM_FOM, 2, "FOM_Kstar.png");
    CalculateFOM(SIGNAL_FOM, CHG_FOM, MIX_FOM, UUBAR_FOM, DDBAR_FOM, SSBAR_FOM, CHARM_FOM, 3, "FOM_Xs.png");

}
