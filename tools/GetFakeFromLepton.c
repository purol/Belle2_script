// last update: 2021-10-08
// for Belle2 data

/*
when you add new variables:
revise # define N_Needed_info ...
revise typedef struct data
revise void Loader::GetData(TFile* input_file)
revise void Loader::PrintRootFile(std::string output_name)
revise void Loader::PrintSeparateRootFile(std::string output_name)
revise void Loader::ConvertIntoSeparateDataFile(std::string output_name, double flag = 0)
*/

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

# define MCTYPE "MC15rd"

/* ====================================== */

double LetsCount_ri(const char* dirname, const char* included_string, double* Ncounts, const char* type, const char* sample, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // modified GetPDFs
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

            Nevt = Nevt + total_weight;

            // calculate nKaon & nPion
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                nKaonfromK = nKaonfromK + temp_N_bin_PID[0][i_PID] * total_weight;
                nKaonformpi = nKaonformpi + temp_N_bin_PID[1][i_PID] * total_weight;
                nPionfrompi = nPionfrompi + temp_N_bin_PID[2][i_PID] * total_weight;
                nPionfromK = nPionfromK + temp_N_bin_PID[3][i_PID] * total_weight;
            }
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                nKaonfrome = nKaonfrome + temp_N_bin_fakeE[0][i_fake] * total_weight; // K- from e
                nKaonfrome = nKaonfrome + temp_N_bin_fakeE[1][i_fake] * total_weight; // K+ from e
                nPionfrome = nPionfrome + temp_N_bin_fakeE[2][i_fake] * total_weight; // pi- from e
                nPionfrome = nPionfrome + temp_N_bin_fakeE[3][i_fake] * total_weight; // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                nKaonfrommu = nKaonfrommu + temp_N_bin_fakeMU[0][i_fake] * total_weight; // K- from mu
                nKaonfrommu = nKaonfrommu + temp_N_bin_fakeMU[1][i_fake] * total_weight; // K+ from mu
                nPionfrommu = nPionfrommu + temp_N_bin_fakeMU[2][i_fake] * total_weight; // pi- from mu
                nPionfrommu = nPionfrommu + temp_N_bin_fakeMU[3][i_fake] * total_weight; // pi+ from mu
            }

        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    Ncounts[0] = Ncounts[0] + nKaonfromK;
    Ncounts[1] = Ncounts[1] + nKaonformpi;
    Ncounts[2] = Ncounts[2] + nKaonfrome;
    Ncounts[3] = Ncounts[3] + nKaonfrommu;
    Ncounts[4] = Ncounts[4] + nPionfrompi;
    Ncounts[5] = Ncounts[5] + nPionfromK;
    Ncounts[6] = Ncounts[6] + nPionfrome;
    Ncounts[7] = Ncounts[7] + nPionfrommu;

    return Nevt;
}

void GetFakeFromLepton() {

    const char* Knunu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/validation_v004/final_output_root_after_MVA_Application_after_cut";
    const char* Kstarnunu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/validation_v004/final_output_root_after_MVA_Application_after_cut";
    const char* Xsununu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/validation_v004/final_output_root_after_MVA_Application_after_cut";
    const char* K0nunu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/validation_v004/final_output_root_after_MVA_Application_after_cut";
    const char* K0starnunu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/validation_v004/final_output_root_after_MVA_Application_after_cut";
    const char* Xsdnunu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/validation_v004/final_output_root_after_MVA_Application_after_cut";
    const char* CHG_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/CHG_analysis/validation_v004/final_output_root_after_MVA_Application_after_cut";
    const char* MIX_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/MIX_analysis/validation_v004/final_output_root_after_MVA_Application_after_cut";
    const char* UUBAR_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/UUBAR_analysis/validation_v004/final_output_root_after_MVA_Application_after_cut";
    const char* DDBAR_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/DDBAR_analysis/validation_v004/final_output_root_after_MVA_Application_after_cut";
    const char* SSBAR_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SSBAR_analysis/validation_v004/final_output_root_after_MVA_Application_after_cut";
    const char* CHARM_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/CHARM_analysis/validation_v004/final_output_root_after_MVA_Application_after_cut";
    
    /*
    * Ncounts:
    * [0]: the number of kaon candidates from K
    * [1]: the number of kaon candidates from pi
    * [2]: the number of kaon candidates from e
    * [3]: the number of kaon candidates from mu
    * [4]: the number of pion candidates from pi
    * [5]: the number of pion candidates from K
    * [6]: the number of pion candidates from e
    * [7]: the number of pion candidates from mu
    */

    double* Ncounts_B2Knunu = (double*)malloc(sizeof(double) * 8);
    double* Ncounts_B2Kstarnunu = (double*)malloc(sizeof(double) * 8);
    double* Ncounts_B2Xsnunu = (double*)malloc(sizeof(double) * 8);
    double* Ncounts_B02K0nunu = (double*)malloc(sizeof(double) * 8);
    double* Ncounts_B02Kstar0nunu = (double*)malloc(sizeof(double) * 8);
    double* Ncounts_B02Xsnunu = (double*)malloc(sizeof(double) * 8);

    double* Ncounts_CHG = (double*)malloc(sizeof(double) * 8);
    double* Ncounts_MIX = (double*)malloc(sizeof(double) * 8);
    double* Ncounts_UUBAR = (double*)malloc(sizeof(double) * 8);
    double* Ncounts_DDBAR = (double*)malloc(sizeof(double) * 8);
    double* Ncounts_SSBAR = (double*)malloc(sizeof(double) * 8);
    double* Ncounts_CHARM = (double*)malloc(sizeof(double) * 8);

    for (int i = 0; i < 8; i++) {
        Ncounts_B2Knunu[i] = 0;
        Ncounts_B2Kstarnunu[i] = 0;
        Ncounts_B2Xsnunu[i] = 0;
        Ncounts_B02K0nunu[i] = 0;
        Ncounts_B02Kstar0nunu[i] = 0;
        Ncounts_B02Xsnunu[i] = 0;

        Ncounts_CHG[i] = 0;
        Ncounts_MIX[i] = 0;
        Ncounts_UUBAR[i] = 0;
        Ncounts_DDBAR[i] = 0;
        Ncounts_SSBAR[i] = 0;
        Ncounts_CHARM[i] = 0;
    }

    double Nevt_B2Knunu = LetsCount_ri(Knunu_dirname, "B2Knunu", Ncounts_B2Knunu, "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Knunu"), "B2Knunu");
    double Nevt_B2Kstarnunu = LetsCount_ri(Kstarnunu_dirname, "B2Kstarnunu", Ncounts_B2Kstarnunu, "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Kstarnunu"), "otherwise");
    double Nevt_B2Xsnunu = LetsCount_ri(Xsununu_dirname, "B2Xsnunu", Ncounts_B2Xsnunu, "Bplus", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B2Xsnunu"), "B2Xsnunu");
    double Nevt_B02K0nunu = LetsCount_ri(K0nunu_dirname, "B02K0nunu", Ncounts_B02K0nunu, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02K0nunu"), "B02K0nunu");
    double Nevt_B02Kstar0nunu = LetsCount_ri(K0starnunu_dirname, "B02Kstar0nunu", Ncounts_B02Kstar0nunu, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Kstar0nunu"), "otherwise");
    double Nevt_B02Xsnunu = LetsCount_ri(Xsdnunu_dirname, "B02Xsnunu", Ncounts_B02Xsnunu, "Bzero", "SIGNAL", ObtainWeight("SIGNAL", MCTYPE, "validation", "B02Xsnunu"), "B02Xsnunu");

    double Nevt_CHG = LetsCount_ri(CHG_dirname, "root", Ncounts_CHG, "Bplus", "CHG", ObtainWeight("CHG", MCTYPE, "validation", "CHG"), "otherwise");
    double Nevt_MIX = LetsCount_ri(MIX_dirname, "root", Ncounts_MIX, "Bzero", "MIX", ObtainWeight("MIX", MCTYPE, "validation", "MIX"), "otherwise");
    double Nevt_UUBAR = LetsCount_ri(UUBAR_dirname, "root", Ncounts_UUBAR, "Continuum", "UUBAR", ObtainWeight("UUBAR", MCTYPE, "validation", "UUBAR"), "otherwise");
    double Nevt_DDBAR = LetsCount_ri(DDBAR_dirname, "root", Ncounts_DDBAR, "Continuum", "DDBAR", ObtainWeight("DDBAR", MCTYPE, "validation", "DDBAR"), "otherwise");
    double Nevt_SSBAR = LetsCount_ri(SSBAR_dirname, "root", Ncounts_SSBAR, "Continuum", "SSBAR", ObtainWeight("SSBAR", MCTYPE, "validation", "SSBAR"), "otherwise");
    double Nevt_CHARM = LetsCount_ri(CHARM_dirname, "root", Ncounts_CHARM, "Continuum", "CHARM", ObtainWeight("CHARM", MCTYPE, "validation", "CHARM"), "otherwise");

    double nKaon_signal_fromk =
        Ncounts_B2Knunu[0] +
        Ncounts_B2Kstarnunu[0] +
        Ncounts_B2Xsnunu[0] +
        Ncounts_B02K0nunu[0] +
        Ncounts_B02Kstar0nunu[0] +
        Ncounts_B02Xsnunu[0];

    double nKaon_bkg_fromk =
        Ncounts_CHG[0] +
        Ncounts_MIX[0] +
        Ncounts_UUBAR[0] +
        Ncounts_DDBAR[0] +
        Ncounts_SSBAR[0] +
        Ncounts_CHARM[0];

    double nKaon_signal_frompi =
        Ncounts_B2Knunu[1] +
        Ncounts_B2Kstarnunu[1] +
        Ncounts_B2Xsnunu[1] +
        Ncounts_B02K0nunu[1] +
        Ncounts_B02Kstar0nunu[1] +
        Ncounts_B02Xsnunu[1];

    double nKaon_bkg_frompi =
        Ncounts_CHG[1] +
        Ncounts_MIX[1] +
        Ncounts_UUBAR[1] +
        Ncounts_DDBAR[1] +
        Ncounts_SSBAR[1] +
        Ncounts_CHARM[1];

    double nKaon_signal_frome =
        Ncounts_B2Knunu[2] +
        Ncounts_B2Kstarnunu[2] +
        Ncounts_B2Xsnunu[2] +
        Ncounts_B02K0nunu[2] +
        Ncounts_B02Kstar0nunu[2] +
        Ncounts_B02Xsnunu[2];

    double nKaon_bkg_frome =
        Ncounts_CHG[2] +
        Ncounts_MIX[2] +
        Ncounts_UUBAR[2] +
        Ncounts_DDBAR[2] +
        Ncounts_SSBAR[2] +
        Ncounts_CHARM[2];

    double nKaon_signal_frommu =
        Ncounts_B2Knunu[3] +
        Ncounts_B2Kstarnunu[3] +
        Ncounts_B2Xsnunu[3] +
        Ncounts_B02K0nunu[3] +
        Ncounts_B02Kstar0nunu[3] +
        Ncounts_B02Xsnunu[3];

    double nKaon_bkg_frommu =
        Ncounts_CHG[3] +
        Ncounts_MIX[3] +
        Ncounts_UUBAR[3] +
        Ncounts_DDBAR[3] +
        Ncounts_SSBAR[3] +
        Ncounts_CHARM[3];

    double nPion_signal_frompi =
        Ncounts_B2Knunu[4] +
        Ncounts_B2Kstarnunu[4] +
        Ncounts_B2Xsnunu[4] +
        Ncounts_B02K0nunu[4] +
        Ncounts_B02Kstar0nunu[4] +
        Ncounts_B02Xsnunu[4];

    double nPion_bkg_frompi =
        Ncounts_CHG[4] +
        Ncounts_MIX[4] +
        Ncounts_UUBAR[4] +
        Ncounts_DDBAR[4] +
        Ncounts_SSBAR[4] +
        Ncounts_CHARM[4];

    double nPion_signal_fromk =
        Ncounts_B2Knunu[5] +
        Ncounts_B2Kstarnunu[5] +
        Ncounts_B2Xsnunu[5] +
        Ncounts_B02K0nunu[5] +
        Ncounts_B02Kstar0nunu[5] +
        Ncounts_B02Xsnunu[5];

    double nPion_bkg_fromk =
        Ncounts_CHG[5] +
        Ncounts_MIX[5] +
        Ncounts_UUBAR[5] +
        Ncounts_DDBAR[5] +
        Ncounts_SSBAR[5] +
        Ncounts_CHARM[5];

    double nPion_signal_frome =
        Ncounts_B2Knunu[6] +
        Ncounts_B2Kstarnunu[6] +
        Ncounts_B2Xsnunu[6] +
        Ncounts_B02K0nunu[6] +
        Ncounts_B02Kstar0nunu[6] +
        Ncounts_B02Xsnunu[6];

    double nPion_bkg_frome =
        Ncounts_CHG[6] +
        Ncounts_MIX[6] +
        Ncounts_UUBAR[6] +
        Ncounts_DDBAR[6] +
        Ncounts_SSBAR[6] +
        Ncounts_CHARM[6];

    double nPion_signal_frommu =
        Ncounts_B2Knunu[7] +
        Ncounts_B2Kstarnunu[7] +
        Ncounts_B2Xsnunu[7] +
        Ncounts_B02K0nunu[7] +
        Ncounts_B02Kstar0nunu[7] +
        Ncounts_B02Xsnunu[7];

    double nPion_bkg_frommu =
        Ncounts_CHG[7] +
        Ncounts_MIX[7] +
        Ncounts_UUBAR[7] +
        Ncounts_DDBAR[7] +
        Ncounts_SSBAR[7] +
        Ncounts_CHARM[7];

    double nKaon_signal_all = nKaon_signal_fromk + nKaon_signal_frompi + nKaon_signal_frome + nKaon_signal_frommu;
    double nKaon_bkg_all = nKaon_bkg_fromk + nKaon_bkg_frompi + nKaon_bkg_frome + nKaon_bkg_frommu;
    double nPion_signal_all = nPion_signal_frompi + nPion_signal_fromk + nPion_signal_frome + nPion_signal_frommu;
    double nPion_bkg_all = nPion_bkg_frompi + nPion_bkg_fromk + nPion_bkg_frome + nPion_bkg_frommu;

    printf("all Kaon at signal: %lf\n", nKaon_signal_all);
    printf("all Kaon at background: %lf\n", nKaon_bkg_all);
    printf("all Pion at signal: %lf\n", nPion_signal_all);
    printf("all Pion at background: %lf\n", nPion_bkg_all);

    printf("Kaon from K at signal: %lf\n", nKaon_signal_fromk);
    printf("Kaon from K at background: %lf\n", nKaon_bkg_fromk);
    printf("Kaon from pi at signal: %lf\n", nKaon_signal_frompi);
    printf("Kaon from pi at background: %lf\n", nKaon_bkg_frompi);
    printf("Kaon from e at signal: %lf\n", nKaon_signal_frome);
    printf("Kaon from e at background: %lf\n", nKaon_bkg_frome);
    printf("Kaon from mu at signal: %lf\n", nKaon_signal_frommu);
    printf("Kaon from mu at background: %lf\n", nKaon_bkg_frommu);

    printf("Pion from pi at signal: %lf\n", nPion_signal_frompi);
    printf("Pion from pi at background: %lf\n", nPion_bkg_frompi);
    printf("Pion from K at signal: %lf\n", nPion_signal_fromk);
    printf("Pion from K at background: %lf\n", nPion_bkg_fromk);
    printf("Pion from e at signal: %lf\n", nPion_signal_frome);
    printf("Pion from e at background: %lf\n", nPion_bkg_frome);
    printf("Pion from mu at signal: %lf\n", nPion_signal_frommu);
    printf("Pion from mu at background: %lf\n", nPion_bkg_frommu);

    printf("Kaon from lepton/all Kaon at signal: %lf\n", (nKaon_signal_frome+ nKaon_signal_frommu)/ nKaon_signal_all);
    printf("Pion from lepton/all Pion at signal: %lf\n", (nPion_signal_frome + nPion_signal_frommu) / nPion_signal_all);
    printf("Kaon from lepton/all Kaon at bkg: %lf\n", (nKaon_bkg_frome + nKaon_bkg_frommu) / nKaon_bkg_all);
    printf("Pion from lepton/all Pion at bkg: %lf\n", (nPion_bkg_frome + nPion_bkg_frommu) / nPion_bkg_all);
}
