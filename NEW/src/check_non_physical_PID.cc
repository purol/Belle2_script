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
#include "constants.h"
#include "ObtainWeight.h"
#include "correctors.h"
#include "base.h"
#include "RooDataHist.h"
#include "template.h"
using namespace RooFit;
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
Corrector_KstarKLKL corrector_KstarKLKL;
Corrector_XsKLKL corrector_XsKLKL;
Corrector_BtoDtoXKL corrector_BtoDtoXKL;
Corrector_Fragmentation corrector_Fragmentation;

# define MCTYPE "MC15rd"

/* ====================================== */

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

            Nevt = Nevt + total_weight;

            FillTemplate(hist, MVA_var, total_weight, Bsig_M);
        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

double GetNonPhysicalPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get nominal PDF with non-physical PID correction factor
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

            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;

            bool IsItNonePhysical = false;

            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE), temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE), temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID

                if ((std::abs(corrector_PID.GetCorrectionFactor(0, i_PID, MCTYPE) - 1) < MyEPSILON) && (temp_N_bin_PID[0][i_PID] > 0.5) && (i_PID != N_PID_syst - 1)) IsItNonePhysical = true;
                if ((std::abs(corrector_PID.GetCorrectionFactor(1, i_PID, MCTYPE) - 1) < MyEPSILON) && (temp_N_bin_PID[1][i_PID] > 0.5) && (i_PID != N_PID_syst - 1)) IsItNonePhysical = true;
                if ((std::abs(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE) - 1) < MyEPSILON) && (temp_N_bin_PID[2][i_PID] > 0.5) && (i_PID != N_PID_syst - 1)) IsItNonePhysical = true;
                if ((std::abs(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE) - 1) < MyEPSILON) && (temp_N_bin_PID[3][i_PID] > 0.5) && (i_PID != N_PID_syst - 1)) IsItNonePhysical = true;
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(0, i_fake, MCTYPE), temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(1, i_fake, MCTYPE), temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(2, i_fake, MCTYPE), temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(3, i_fake, MCTYPE), temp_N_bin_fakeE[3][i_fake]); // pi+ from e

                if ((std::abs(corrector_FakePID.GetCorrectionFactorfakeE(0, i_fake, MCTYPE) - 1) < MyEPSILON) && (temp_N_bin_fakeE[0][i_fake] > 0.5) && (i_fake != N_fakeE_syst - 1)) IsItNonePhysical = true;
                if ((std::abs(corrector_FakePID.GetCorrectionFactorfakeE(1, i_fake, MCTYPE) - 1) < MyEPSILON) && (temp_N_bin_fakeE[1][i_fake] > 0.5) && (i_fake != N_fakeE_syst - 1)) IsItNonePhysical = true;
                if ((std::abs(corrector_FakePID.GetCorrectionFactorfakeE(2, i_fake, MCTYPE) - 1) < MyEPSILON) && (temp_N_bin_fakeE[2][i_fake] > 0.5) && (i_fake != N_fakeE_syst - 1)) IsItNonePhysical = true;
                if ((std::abs(corrector_FakePID.GetCorrectionFactorfakeE(3, i_fake, MCTYPE) - 1) < MyEPSILON) && (temp_N_bin_fakeE[3][i_fake] > 0.5) && (i_fake != N_fakeE_syst - 1)) IsItNonePhysical = true;
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(0, i_fake, MCTYPE), temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(1, i_fake, MCTYPE), temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(2, i_fake, MCTYPE), temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(3, i_fake, MCTYPE), temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu

                if ((std::abs(corrector_FakePID.GetCorrectionFactorfakeMU(0, i_fake, MCTYPE) - 1) < MyEPSILON) && (temp_N_bin_fakeMU[0][i_fake] > 0.5) && (i_fake != N_fakeMU_syst - 1)) IsItNonePhysical = true;
                if ((std::abs(corrector_FakePID.GetCorrectionFactorfakeMU(1, i_fake, MCTYPE) - 1) < MyEPSILON) && (temp_N_bin_fakeMU[1][i_fake] > 0.5) && (i_fake != N_fakeMU_syst - 1)) IsItNonePhysical = true;
                if ((std::abs(corrector_FakePID.GetCorrectionFactorfakeMU(2, i_fake, MCTYPE) - 1) < MyEPSILON) && (temp_N_bin_fakeMU[2][i_fake] > 0.5) && (i_fake != N_fakeMU_syst - 1)) IsItNonePhysical = true;
                if ((std::abs(corrector_FakePID.GetCorrectionFactorfakeMU(3, i_fake, MCTYPE) - 1) < MyEPSILON) && (temp_N_bin_fakeMU[3][i_fake] > 0.5) && (i_fake != N_fakeMU_syst - 1)) IsItNonePhysical = true;
            }

            if (IsItNonePhysical == false) continue;

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

            Nevt = Nevt + total_weight;

            FillTemplate(hist, MVA_var, total_weight, Bsig_M);
        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

int main()
{

    /* ====================================== */
    // Define PDFs for HistFactory
    // nominal PDFs
    TH1D* all_nominal = new TH1D("all_nominal", "all_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* non_physicsl_nominal = new TH1D("non_physicsl_nominal", "non_physicsl_nominal", RarityBins, BinMIN, BinMAX);
    TH1D* ratio_nominal = new TH1D("ratio_nominal", "ratio_nominal", RarityBins, BinMIN, BinMAX);
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
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Knunu", all_nominal, "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", all_nominal, "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    GetNominalPDFs(MC_dirname_SIGNAL, "B2Xsnunu", all_nominal, "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02K0nunu", all_nominal, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", all_nominal, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    GetNominalPDFs(MC_dirname_SIGNAL, "B02Xsnunu", all_nominal, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    GetNominalPDFs(MC_dirname_CHG, "root", all_nominal, "Bplus", "CHG", ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    GetNominalPDFs(MC_dirname_MIX, "root", all_nominal, "Bzero", "MIX", ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");
    GetNominalPDFs(MC_dirname_UUBAR, "root", all_nominal, "Continuum", "UUBAR", ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"), "otherwise");
    GetNominalPDFs(MC_dirname_DDBAR, "root", all_nominal, "Continuum", "DDBAR", ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"), "otherwise");
    GetNominalPDFs(MC_dirname_SSBAR, "root", all_nominal, "Continuum", "SSBAR", ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"), "otherwise");
    GetNominalPDFs(MC_dirname_CHARM, "root", all_nominal, "Continuum", "CHARM", ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"), "otherwise");

    // get non-physical pdfs
    GetNonPhysicalPDFs(MC_dirname_SIGNAL, "B2Knunu", non_physicsl_nominal, "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    GetNonPhysicalPDFs(MC_dirname_SIGNAL, "B2Kstarnunu", non_physicsl_nominal, "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    GetNonPhysicalPDFs(MC_dirname_SIGNAL, "B2Xsnunu", non_physicsl_nominal, "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    GetNonPhysicalPDFs(MC_dirname_SIGNAL, "B02K0nunu", non_physicsl_nominal, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    GetNonPhysicalPDFs(MC_dirname_SIGNAL, "B02Kstar0nunu", non_physicsl_nominal, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    GetNonPhysicalPDFs(MC_dirname_SIGNAL, "B02Xsnunu", non_physicsl_nominal, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    GetNonPhysicalPDFs(MC_dirname_CHG, "root", non_physicsl_nominal, "Bplus", "CHG", ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    GetNonPhysicalPDFs(MC_dirname_MIX, "root", non_physicsl_nominal, "Bzero", "MIX", ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");
    GetNonPhysicalPDFs(MC_dirname_UUBAR, "root", non_physicsl_nominal, "Continuum", "UUBAR", ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"), "otherwise");
    GetNonPhysicalPDFs(MC_dirname_DDBAR, "root", non_physicsl_nominal, "Continuum", "DDBAR", ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"), "otherwise");
    GetNonPhysicalPDFs(MC_dirname_SSBAR, "root", non_physicsl_nominal, "Continuum", "SSBAR", ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"), "otherwise");
    GetNonPhysicalPDFs(MC_dirname_CHARM, "root", non_physicsl_nominal, "Continuum", "CHARM", ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"), "otherwise");
    /* ====================================== */

    ratio_nominal->Divide(non_physicsl_nominal, all_nominal);

    printf("print x-asix list for python input:\n");
    printf("[");
    for (int i = 0; i < RarityBins - 1; i++) printf("%d, ", i);
    printf("%d]\n", RarityBins - 1);

    printf("print x-asix error list for python input:\n");
    printf("[");
    for (int i = 0; i < RarityBins - 1; i++) printf("%lf, ", 0.5);
    printf("%lf]\n", 0.5);

    printf("print y-asix list for python input:\n");
    printf("[");
    for (int i = 0; i < RarityBins - 1; i++) printf("%lf, ", ratio_nominal->GetBinContent(i + 1));
    printf("%lf]\n", ratio_nominal->GetBinContent(RarityBins));

    printf("print y-asix error list for python input:\n");
    printf("[");
    for (int i = 0; i < RarityBins - 1; i++) printf("%lf, ", ratio_nominal->GetBinError(i + 1));
    printf("%lf]\n", ratio_nominal->GetBinError(RarityBins));


    /* ====================================== */

    return 0;
}
