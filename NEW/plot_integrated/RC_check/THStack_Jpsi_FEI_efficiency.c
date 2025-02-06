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

# define MCTYPE "MC15rd"

double CAL = 1.0; // must be 1.0
# define CAL_qq_K 1.0
# define CAL_qq_Kstar 1.0
# define CAL_qq_Xs 1.0

# include <algorithm>
# include <float.h>
# include <string>
# include <cstring>
# include <vector>

#include "constants.h"
#include "ObtainWeight.h"
#include "correctors.h"
#include "base.h"

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"

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
Corrector_LID corrector_LID;

TH1D* Nevt_MC_Bplus = new TH1D("Nevt_MC_Bplus", "Nevt_MC_Bplus", 36, 0, 36);
TH1D* Nevt_MC_Bzero = new TH1D("Nevt_MC_Bzero", "Nevt_MC_Bzero", 32, 0, 32);

TH1D* Nevt_data_Bplus = new TH1D("Nevt_data_Bplus", "Nevt_data_Bplus", 36, 0, 36);
TH1D* Nevt_data_Bzero = new TH1D("Nevt_data_Bzero", "Nevt_data_Bzero", 32, 0, 32);

TH1D* ratio_Bplus = new TH1D("ratio_Bplus", "ratio_Bplus", 36, 0, 36);
TH1D* ratio_Bzero = new TH1D("ratio_Bzero", "ratio_Bzero", 32, 0, 32);


TH1D* onebin_MC_Bplus = new TH1D("onebin_MC_Bplus", "onebin_MC_Bplus", 1, 0, 1);
TH1D* onebin_MC_Bzero = new TH1D("onebin_MC_Bzero", "onebin_MC_Bzero", 1, 0, 1);

TH1D* onebin_data_Bplus = new TH1D("onebin_data_Bplus", "onebin_data_Bplus", 1, 0, 1);
TH1D* onebin_data_Bzero = new TH1D("onebin_data_Bzero", "onebin_data_Bzero", 1, 0, 1);

TH1D* ratio_onebin_Bplus = new TH1D("ratio_onebin_Bplus", "ratio_Bplus", 1, 0, 1);
TH1D* ratio_onebin_Bzero = new TH1D("ratio_onebin_Bzero", "ratio_Bzero", 1, 0, 1);


TH1D* Nevt_MC_Bplus_compare = new TH1D("Nevt_MC_Bplus_compare", "Nevt_MC_Bplus_compare", 12, 0, 12);
TH1D* Nevt_MC_Bzero_compare = new TH1D("Nevt_MC_Bzero_compare", "Nevt_MC_Bzero_compare", 11, 0, 11);

TH1D* Nevt_data_Bplus_compare = new TH1D("Nevt_data_Bplus_compare", "Nevt_data_Bplus_compare", 12, 0, 12);
TH1D* Nevt_data_Bzero_compare = new TH1D("Nevt_data_Bzero_compare", "Nevt_data_Bzero_compare", 11, 0, 11);

TH1D* ratio_Bplus_compare = new TH1D("ratio_Bplus_compare", "ratio_Bplus_compare", 12, 0, 12);
TH1D* ratio_Bzero_compare = new TH1D("ratio_Bzero_compare", "ratio_Bzero_compare", 11, 0, 11);


TH1D* Nevt_MC_Bplus_My = new TH1D("Nevt_MC_Bplus_My", "Nevt_MC_Bplus_My", 3, 0, 3);
TH1D* Nevt_MC_Bzero_My = new TH1D("Nevt_MC_Bzero_My", "Nevt_MC_Bzero_My", 1, 0, 1);

TH1D* Nevt_data_Bplus_My = new TH1D("Nevt_data_Bplus_My", "Nevt_data_Bplus_My", 3, 0, 3);
TH1D* Nevt_data_Bzero_My = new TH1D("Nevt_data_Bzero_My", "Nevt_data_Bzero_My", 1, 0, 1);

TH1D* ratio_Bplus_My = new TH1D("ratio_Bplus_My", "ratio_Bplus_My", 3, 0, 3);
TH1D* ratio_Bzero_My = new TH1D("ratio_Bzero_My", "ratio_Bzero_My", 1, 0, 1);

void LetsFillKonlyMC(const char* dirname, std::string SampleName, int option = 0, double additional_weight = 1.0, bool IncludeKstar = false) {
    /*
    SampleName for Knn
    CHG
    MIX
    UUBAR
    DDBAR
    SSBAR
    CHARM
    SIGNAL
    */
    /*
    option 0: select all Btag
    option 1: select Btag+
    option 2: select Btag0
    */
    /*
    0: charged
    1: mixed
    2: uubar
    3: ddbar
    4: ssbar
    5: ccbar
    6: tautau
    7: mumu
    8: gg
    9: ee
    10: eeee
    11: eemumu
    12: llXX
    13: hhISR
    */

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    double FEI_calibration_factor = -1;

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

    int Ngamma_v200_index = -1;
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

    int nBp = -1;
    int nB0 = -1;

    double Bsig_M = -1;

    std::vector<std::string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");
        TTree* tree_Xs;
        if (SampleName == "SIGNAL") tree_Xs = (TTree*)input_file->Get("Xs");

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID); // charged: 0, mixed: 1
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
        if (SampleName == "SIGNAL") {
            tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clPrimaryMC__bc", &nBp);
            tree_Xs->SetBranchAddress("nParticlesInList__boB0__clPrimaryMC__bc", &nB0);
        }
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
            if (SampleName == "SIGNAL") tree_Xs->GetEntry(j);

            if (IncludeKstar == false) {
                // select B->K Jpsi only
                if ((Bsig_ID > -0.5) && (Bsig_ID < 0.5)) {}
                else continue;
            }
            else {
                // select B->k(*) Jpsi only
                if ((Bsig_ID > -0.5) && (Bsig_ID < 2.5) && (Bsig_M < 1.0)) {}
                else continue;
            }

            if (option == 1 && Upsilon_ID != 0) continue;
            else if (option == 2 && Upsilon_ID != 1) continue;
            else if (option == 0) {}
            else if ((option != 0) && (option != 1) && (option != 2)) {
                printf("improper option value!\n");
                exit(1);
            }

            double weight_ri = 0.0;
            if (SampleName == "CHG") {
                FEI_calibration_factor = 1.0;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "MIX") {
                FEI_calibration_factor = 1.0;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "UUBAR") {
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) FEI_calibration_factor = CAL_qq_K;
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) FEI_calibration_factor = CAL_qq_Kstar;
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) FEI_calibration_factor = CAL_qq_Xs;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "DDBAR") {
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) FEI_calibration_factor = CAL_qq_K;
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) FEI_calibration_factor = CAL_qq_Kstar;
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) FEI_calibration_factor = CAL_qq_Xs;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "SSBAR") {
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) FEI_calibration_factor = CAL_qq_K;
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) FEI_calibration_factor = CAL_qq_Kstar;
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) FEI_calibration_factor = CAL_qq_Xs;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "CHARM") {
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) FEI_calibration_factor = CAL_qq_K;
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) FEI_calibration_factor = CAL_qq_Kstar;
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) FEI_calibration_factor = CAL_qq_Xs;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "SIGNAL") {
                if (nBp > 0) {
                    FEI_calibration_factor = 1.0;
                    weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
                }
                else if (nB0 > 0) {
                    FEI_calibration_factor = 1.0;
                    weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
                }
                else {
                    printf("ERROR 255");
                    exit(1);
                }
            }
            else {
                printf("undefined job id!\n");
                exit(1);
            }

            // Fill calibration factors
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
            if (SampleName == "CHG" || SampleName == "MIX" || SampleName == "SIGNAL") Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double weights = FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * additional_weight;

            if ((Upsilon_ID > -0.5) && (Upsilon_ID < 0.5)) {
                Nevt_MC_Bplus->Fill(Btag_ID + 0.5, weights);
                onebin_MC_Bplus->Fill(0.5, weights);
                if (std::abs(Btag_ID - 0.0) < MyEPSILON) Nevt_MC_Bplus_compare->Fill(0.5, weights);
                else if (std::abs(Btag_ID - 1.0) < MyEPSILON) Nevt_MC_Bplus_compare->Fill(1.5, weights);
                else if (std::abs(Btag_ID - 3.0) < MyEPSILON) Nevt_MC_Bplus_compare->Fill(2.5, weights);
                else if (std::abs(Btag_ID - 4.0) < MyEPSILON) Nevt_MC_Bplus_compare->Fill(3.5, weights);
                else if (std::abs(Btag_ID - 15.0) < MyEPSILON) Nevt_MC_Bplus_compare->Fill(4.5, weights);
                else if (std::abs(Btag_ID - 16.0) < MyEPSILON) Nevt_MC_Bplus_compare->Fill(5.5, weights);
                else if (std::abs(Btag_ID - 18.0) < MyEPSILON) Nevt_MC_Bplus_compare->Fill(6.5, weights);
                else if (std::abs(Btag_ID - 19.0) < MyEPSILON) Nevt_MC_Bplus_compare->Fill(7.5, weights);
                else if (std::abs(Btag_ID - 23.0) < MyEPSILON) Nevt_MC_Bplus_compare->Fill(8.5, weights);
                else if (std::abs(Btag_ID - 24.0) < MyEPSILON) Nevt_MC_Bplus_compare->Fill(9.5, weights);
                else if (std::abs(Btag_ID - 30.0) < MyEPSILON) Nevt_MC_Bplus_compare->Fill(10.5, weights);
                else Nevt_MC_Bplus_compare->Fill(11.5, weights);
                if (std::abs(Btag_ID - 0.0) < MyEPSILON) Nevt_MC_Bplus_My->Fill(0.5, weights);
                else if (std::abs(Btag_ID - 1.0) < MyEPSILON) Nevt_MC_Bplus_My->Fill(0.5, weights);
                else if (std::abs(Btag_ID - 3.0) < MyEPSILON) Nevt_MC_Bplus_My->Fill(0.5, weights);
                else if (std::abs(Btag_ID - 4.0) < MyEPSILON) Nevt_MC_Bplus_My->Fill(0.5, weights);
                else if (std::abs(Btag_ID - 15.0) < MyEPSILON) Nevt_MC_Bplus_My->Fill(1.5, weights);
                else if (std::abs(Btag_ID - 16.0) < MyEPSILON) Nevt_MC_Bplus_My->Fill(1.5, weights);
                else if (std::abs(Btag_ID - 18.0) < MyEPSILON) Nevt_MC_Bplus_My->Fill(1.5, weights);
                else if (std::abs(Btag_ID - 19.0) < MyEPSILON) Nevt_MC_Bplus_My->Fill(1.5, weights);
                else if (std::abs(Btag_ID - 23.0) < MyEPSILON) Nevt_MC_Bplus_My->Fill(2.5, weights);
                else if (std::abs(Btag_ID - 24.0) < MyEPSILON) Nevt_MC_Bplus_My->Fill(2.5, weights);
                else if (std::abs(Btag_ID - 30.0) < MyEPSILON) Nevt_MC_Bplus_My->Fill(2.5, weights);
                else Nevt_MC_Bplus_My->Fill(2.5, weights);
            }
            else {
                Nevt_MC_Bzero->Fill(Btag_ID + 0.5, weights);
                onebin_MC_Bzero->Fill(0.5, weights);
                if (std::abs(Btag_ID - 0.0) < MyEPSILON) Nevt_MC_Bzero_compare->Fill(0.5, weights);
                else if (std::abs(Btag_ID - 1.0) < MyEPSILON) Nevt_MC_Bzero_compare->Fill(1.5, weights);
                else if (std::abs(Btag_ID - 3.0) < MyEPSILON) Nevt_MC_Bzero_compare->Fill(2.5, weights);
                else if (std::abs(Btag_ID - 4.0) < MyEPSILON) Nevt_MC_Bzero_compare->Fill(3.5, weights);
                else if (std::abs(Btag_ID - 5.0) < MyEPSILON) Nevt_MC_Bzero_compare->Fill(4.5, weights);
                else if (std::abs(Btag_ID - 15.0) < MyEPSILON) Nevt_MC_Bzero_compare->Fill(5.5, weights);
                else if (std::abs(Btag_ID - 16.0) < MyEPSILON) Nevt_MC_Bzero_compare->Fill(6.5, weights);
                else if (std::abs(Btag_ID - 18.0) < MyEPSILON) Nevt_MC_Bzero_compare->Fill(7.5, weights);
                else if (std::abs(Btag_ID - 19.0) < MyEPSILON) Nevt_MC_Bzero_compare->Fill(8.5, weights);
                else if (std::abs(Btag_ID - 26.0) < MyEPSILON) Nevt_MC_Bzero_compare->Fill(9.5, weights);
                else Nevt_MC_Bzero_compare->Fill(10.5, weights);
                Nevt_MC_Bzero_My->Fill(0.5, weights);
            }

        }
        input_file->Close();

    }

}

void LetsFillKonlydata(const char* dirname, const char* included_string = "root", int option = 0, bool IncludeKstar = false) {
    /*
    option 0: select all Btag
    option 1: select Btag+
    option 2: select Btag0
    */
    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;

    double Bsig_M = -1;

    std::vector<std::string> names;
    load_files(dirname, &names, included_string);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID); // charged: 0, mixed: 1
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            if (IncludeKstar == false) {
                // select B->K Jpsi only
                if ((Bsig_ID > -0.5) && (Bsig_ID < 0.5)) {}
                else continue;
            }
            else {
                // select B->k(*) Jpsi only
                if ((Bsig_ID > -0.5) && (Bsig_ID < 2.5) && (Bsig_M < 1.0)) {}
                else continue;
            }

            if (option == 1 && Upsilon_ID != 0) continue;
            else if (option == 2 && Upsilon_ID != 1) continue;
            else if (option == 0) {}
            else if ((option != 0) && (option != 1) && (option != 2)) {
                printf("improper option value!\n");
                exit(1);
            }

            if ((Upsilon_ID > -0.5) && (Upsilon_ID < 0.5)) {
                Nevt_data_Bplus->Fill(Btag_ID + 0.5, 1.0);
                onebin_data_Bplus->Fill(0.5, 1.0);
                if (std::abs(Btag_ID - 0.0) < MyEPSILON) Nevt_data_Bplus_compare->Fill(0.5, 1.0);
                else if (std::abs(Btag_ID - 1.0) < MyEPSILON) Nevt_data_Bplus_compare->Fill(1.5, 1.0);
                else if (std::abs(Btag_ID - 3.0) < MyEPSILON) Nevt_data_Bplus_compare->Fill(2.5, 1.0);
                else if (std::abs(Btag_ID - 4.0) < MyEPSILON) Nevt_data_Bplus_compare->Fill(3.5, 1.0);
                else if (std::abs(Btag_ID - 15.0) < MyEPSILON) Nevt_data_Bplus_compare->Fill(4.5, 1.0);
                else if (std::abs(Btag_ID - 16.0) < MyEPSILON) Nevt_data_Bplus_compare->Fill(5.5, 1.0);
                else if (std::abs(Btag_ID - 18.0) < MyEPSILON) Nevt_data_Bplus_compare->Fill(6.5, 1.0);
                else if (std::abs(Btag_ID - 19.0) < MyEPSILON) Nevt_data_Bplus_compare->Fill(7.5, 1.0);
                else if (std::abs(Btag_ID - 23.0) < MyEPSILON) Nevt_data_Bplus_compare->Fill(8.5, 1.0);
                else if (std::abs(Btag_ID - 24.0) < MyEPSILON) Nevt_data_Bplus_compare->Fill(9.5, 1.0);
                else if (std::abs(Btag_ID - 30.0) < MyEPSILON) Nevt_data_Bplus_compare->Fill(10.5, 1.0);
                else Nevt_data_Bplus_compare->Fill(11.5, 1.0);
                if (std::abs(Btag_ID - 0.0) < MyEPSILON) Nevt_data_Bplus_My->Fill(0.5, 1.0);
                else if (std::abs(Btag_ID - 1.0) < MyEPSILON) Nevt_data_Bplus_My->Fill(0.5, 1.0);
                else if (std::abs(Btag_ID - 3.0) < MyEPSILON) Nevt_data_Bplus_My->Fill(0.5, 1.0);
                else if (std::abs(Btag_ID - 4.0) < MyEPSILON) Nevt_data_Bplus_My->Fill(0.5, 1.0);
                else if (std::abs(Btag_ID - 15.0) < MyEPSILON) Nevt_data_Bplus_My->Fill(1.5, 1.0);
                else if (std::abs(Btag_ID - 16.0) < MyEPSILON) Nevt_data_Bplus_My->Fill(1.5, 1.0);
                else if (std::abs(Btag_ID - 18.0) < MyEPSILON) Nevt_data_Bplus_My->Fill(1.5, 1.0);
                else if (std::abs(Btag_ID - 19.0) < MyEPSILON) Nevt_data_Bplus_My->Fill(1.5, 1.0);
                else if (std::abs(Btag_ID - 23.0) < MyEPSILON) Nevt_data_Bplus_My->Fill(2.5, 1.0);
                else if (std::abs(Btag_ID - 24.0) < MyEPSILON) Nevt_data_Bplus_My->Fill(2.5, 1.0);
                else if (std::abs(Btag_ID - 30.0) < MyEPSILON) Nevt_data_Bplus_My->Fill(2.5, 1.0);
                else Nevt_data_Bplus_My->Fill(2.5, 1.0);
            }
            else {
                Nevt_data_Bzero->Fill(Btag_ID + 0.5, 1.0);
                onebin_data_Bzero->Fill(0.5, 1.0);
                if (std::abs(Btag_ID - 0.0) < MyEPSILON) Nevt_data_Bzero_compare->Fill(0.5, 1.0);
                else if (std::abs(Btag_ID - 1.0) < MyEPSILON) Nevt_data_Bzero_compare->Fill(1.5, 1.0);
                else if (std::abs(Btag_ID - 3.0) < MyEPSILON) Nevt_data_Bzero_compare->Fill(2.5, 1.0);
                else if (std::abs(Btag_ID - 4.0) < MyEPSILON) Nevt_data_Bzero_compare->Fill(3.5, 1.0);
                else if (std::abs(Btag_ID - 5.0) < MyEPSILON) Nevt_data_Bzero_compare->Fill(4.5, 1.0);
                else if (std::abs(Btag_ID - 15.0) < MyEPSILON) Nevt_data_Bzero_compare->Fill(5.5, 1.0);
                else if (std::abs(Btag_ID - 16.0) < MyEPSILON) Nevt_data_Bzero_compare->Fill(6.5, 1.0);
                else if (std::abs(Btag_ID - 18.0) < MyEPSILON) Nevt_data_Bzero_compare->Fill(7.5, 1.0);
                else if (std::abs(Btag_ID - 19.0) < MyEPSILON) Nevt_data_Bzero_compare->Fill(8.5, 1.0);
                else if (std::abs(Btag_ID - 26.0) < MyEPSILON) Nevt_data_Bzero_compare->Fill(9.5, 1.0);
                else Nevt_data_Bzero_compare->Fill(10.5, 1.0);
                Nevt_data_Bzero_My->Fill(0.5, 1.0);
            }

        }
        input_file->Close();

    }

}

void THStack_Jpsi_FEI_efficiency() {
    /*
    * argv[1]: input version (ex. v031)
    * argv[2]: output path
    */


    // dirnames
    const char* Jpsi_MC_SIGNAL_before_dirname = strdup(("/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_MC_Jpsi/SIGNAL_analysis/validation_" + std::string(argv[1]) + "/final_output_root_after_MVA_Application").c_str());
    const char* Jpsi_MC_CHG_before_dirname = strdup(("/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_MC_Jpsi/CHG_analysis/validation_" + std::string(argv[1]) + "/final_output_root_after_MVA_Application").c_str());
    const char* Jpsi_MC_MIX_before_dirname = strdup(("/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_MC_Jpsi/MIX_analysis/validation_" + std::string(argv[1]) + "/final_output_root_after_MVA_Application").c_str());
    const char* Jpsi_MC_UUBAR_before_dirname = strdup(("/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_MC_Jpsi/UUBAR_analysis/validation_" + std::string(argv[1]) + "/final_output_root_after_MVA_Application").c_str());
    const char* Jpsi_MC_DDBAR_before_dirname = strdup(("/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_MC_Jpsi/DDBAR_analysis/validation_" + std::string(argv[1]) + "/final_output_root_after_MVA_Application").c_str());
    const char* Jpsi_MC_SSBAR_before_dirname = strdup(("/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_MC_Jpsi/SSBAR_analysis/validation_" + std::string(argv[1]) + "/final_output_root_after_MVA_Application").c_str());
    const char* Jpsi_MC_CHARM_before_dirname = strdup(("/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_MC_Jpsi/CHARM_analysis/validation_" + std::string(argv[1]) + "/final_output_root_after_MVA_Application").c_str());
    const char* Jpsi_data_before_dirname = strdup(("/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_data_Jpsi/SIGNAL_analysis/validation_" + std::string(argv[1]) + "/final_output_root_after_MVA_Application").c_str());

    LetsFillKonlyMC(Jpsi_MC_SIGNAL_before_dirname, "SIGNAL", 0, 1.0, true);
    LetsFillKonlyMC(Jpsi_MC_CHG_before_dirname, "CHG", 0, 1.0, true);
    LetsFillKonlyMC(Jpsi_MC_MIX_before_dirname, "MIX", 0, 1.0, true);
    LetsFillKonlyMC(Jpsi_MC_UUBAR_before_dirname, "UUBAR", 0, 1.0, true);
    LetsFillKonlyMC(Jpsi_MC_DDBAR_before_dirname, "DDBAR", 0, 1.0, true);
    LetsFillKonlyMC(Jpsi_MC_SSBAR_before_dirname, "SSBAR", 0, 1.0, true);
    LetsFillKonlyMC(Jpsi_MC_CHARM_before_dirname, "CHARM", 0, 1.0, true);
    LetsFillKonlydata(Jpsi_data_before_dirname, "root", 0, true);

    for (int i = 0; i < 36; i++) ratio_Bplus->Divide(Nevt_data_Bplus, Nevt_MC_Bplus);
    for (int i = 0; i < 32; i++) ratio_Bzero->Divide(Nevt_data_Bzero, Nevt_MC_Bzero);

    ratio_onebin_Bplus->Divide(onebin_data_Bplus, onebin_MC_Bplus);
    ratio_onebin_Bzero->Divide(onebin_data_Bzero, onebin_MC_Bzero);

    ratio_Bplus_compare->Divide(Nevt_data_Bplus_compare, Nevt_MC_Bplus_compare);
    ratio_Bzero_compare->Divide(Nevt_data_Bzero_compare, Nevt_MC_Bzero_compare);

    ratio_Bplus_My->Divide(Nevt_data_Bplus_My, Nevt_MC_Bplus_My);
    ratio_Bzero_My->Divide(Nevt_data_Bzero_My, Nevt_MC_Bzero_My);

    printf("=== Nevt MC for Bplus ===\n");
    for (int i = 0; i < 36; i++) printf("channel %d: %lf +- %lf\n", i, Nevt_MC_Bplus->GetBinContent(i + 1), Nevt_MC_Bplus->GetBinError(i + 1));
    printf("=========================\n");

    printf("=== Nevt MC for Bzero ===\n");
    for (int i = 0; i < 32; i++) printf("channel %d: %lf +- %lf\n", i, Nevt_MC_Bzero->GetBinContent(i + 1), Nevt_MC_Bzero->GetBinError(i + 1));
    printf("=========================\n");

    printf("=== Nevt data for Bplus ===\n");
    for (int i = 0; i < 36; i++) printf("channel %d: %lf +- %lf\n", i, Nevt_data_Bplus->GetBinContent(i + 1), Nevt_data_Bplus->GetBinError(i + 1));
    printf("===========================\n");

    printf("=== Nevt data for Bzero ===\n");
    for (int i = 0; i < 32; i++) printf("channel %d: %lf +- %lf\n", i, Nevt_data_Bzero->GetBinContent(i + 1), Nevt_data_Bzero->GetBinError(i + 1));
    printf("===========================\n");

    printf("=== data/MC for Bplus ===\n");
    for (int i = 0; i < 36; i++) printf("channel %d: %lf +- %lf\n", i, ratio_Bplus->GetBinContent(i + 1), ratio_Bplus->GetBinError(i + 1));
    printf("==========================\n");

    printf("=== data/MC for Bzero ===\n");
    for (int i = 0; i < 32; i++) printf("channel %d: %lf +- %lf\n", i, ratio_Bzero->GetBinContent(i + 1), ratio_Bzero->GetBinError(i + 1));
    printf("==========================\n");



    printf("=== Nevt MC for Bplus ===\n");
    printf("all channel: %lf +- %lf\n", onebin_MC_Bplus->GetBinContent(1), onebin_MC_Bplus->GetBinError(1));
    printf("=========================\n");

    printf("=== Nevt MC for Bzero ===\n");
    printf("all channel: %lf +- %lf\n", onebin_MC_Bzero->GetBinContent(1), onebin_MC_Bzero->GetBinError(1));
    printf("=========================\n");

    printf("=== Nevt data for Bplus ===\n");
    printf("all channel: %lf +- %lf\n", onebin_data_Bplus->GetBinContent(1), onebin_data_Bplus->GetBinError(1));
    printf("===========================\n");

    printf("=== Nevt data for Bzero ===\n");
    printf("all channel: %lf +- %lf\n", onebin_data_Bzero->GetBinContent(1), onebin_data_Bzero->GetBinError(1));
    printf("===========================\n");

    printf("=== data/MC for Bplus ===\n");
    printf("all channel: %lf +- %lf\n", ratio_onebin_Bplus->GetBinContent(1), ratio_onebin_Bplus->GetBinError(1));
    printf("==========================\n");

    printf("=== data/MC for Bzero ===\n");
    printf("all channel: %lf +- %lf\n", ratio_onebin_Bzero->GetBinContent(1), ratio_onebin_Bzero->GetBinError(1));
    printf("==========================\n");


    printf("=== Nevt MC for Bplus ===\n");
    for (int i = 0; i < 12; i++) printf("category %d: %lf +- %lf\n", i, Nevt_MC_Bplus_compare->GetBinContent(i + 1), Nevt_MC_Bplus_compare->GetBinError(i + 1));
    printf("=========================\n");

    printf("=== Nevt MC for Bzero ===\n");
    for (int i = 0; i < 11; i++) printf("category %d: %lf +- %lf\n", i, Nevt_MC_Bzero_compare->GetBinContent(i + 1), Nevt_MC_Bzero_compare->GetBinError(i + 1));
    printf("=========================\n");

    printf("=== Nevt data for Bplus ===\n");
    for (int i = 0; i < 12; i++) printf("category %d: %lf +- %lf\n", i, Nevt_data_Bplus_compare->GetBinContent(i + 1), Nevt_data_Bplus_compare->GetBinError(i + 1));
    printf("===========================\n");

    printf("=== Nevt data for Bzero ===\n");
    for (int i = 0; i < 11; i++) printf("category %d: %lf +- %lf\n", i, Nevt_data_Bzero_compare->GetBinContent(i + 1), Nevt_data_Bzero_compare->GetBinError(i + 1));
    printf("===========================\n");

    printf("=== data/MC for Bplus ===\n");
    for (int i = 0; i < 12; i++) printf("category %d: %lf +- %lf\n", i, ratio_Bplus_compare->GetBinContent(i + 1), ratio_Bplus_compare->GetBinError(i + 1));
    printf("==========================\n");

    printf("=== data/MC for Bzero ===\n");
    for (int i = 0; i < 11; i++) printf("category %d: %lf +- %lf\n", i, ratio_Bzero_compare->GetBinContent(i + 1), ratio_Bzero_compare->GetBinError(i + 1));
    printf("==========================\n");


    printf("=== Nevt MC for Bplus ===\n");
    for (int i = 0; i < 3; i++) printf("My category %d: %lf +- %lf\n", i, Nevt_MC_Bplus_My->GetBinContent(i + 1), Nevt_MC_Bplus_My->GetBinError(i + 1));
    printf("=========================\n");

    printf("=== Nevt MC for Bzero ===\n");
    for (int i = 0; i < 1; i++) printf("My category %d: %lf +- %lf\n", i, Nevt_MC_Bzero_My->GetBinContent(i + 1), Nevt_MC_Bzero_My->GetBinError(i + 1));
    printf("=========================\n");

    printf("=== Nevt data for Bplus ===\n");
    for (int i = 0; i < 3; i++) printf("My category %d: %lf +- %lf\n", i, Nevt_data_Bplus_My->GetBinContent(i + 1), Nevt_data_Bplus_My->GetBinError(i + 1));
    printf("===========================\n");

    printf("=== Nevt data for Bzero ===\n");
    for (int i = 0; i < 1; i++) printf("My category %d: %lf +- %lf\n", i, Nevt_data_Bzero_My->GetBinContent(i + 1), Nevt_data_Bzero_My->GetBinError(i + 1));
    printf("===========================\n");

    printf("=== data/MC for Bplus ===\n");
    for (int i = 0; i < 3; i++) printf("My category %d: %lf +- %lf\n", i, ratio_Bplus_My->GetBinContent(i + 1), ratio_Bplus_My->GetBinError(i + 1));
    printf("==========================\n");

    printf("=== data/MC for Bzero ===\n");
    for (int i = 0; i < 1; i++) printf("My category %d: %lf +- %lf\n", i, ratio_Bzero_My->GetBinContent(i + 1), ratio_Bzero_My->GetBinError(i + 1));
    printf("==========================\n");
}
