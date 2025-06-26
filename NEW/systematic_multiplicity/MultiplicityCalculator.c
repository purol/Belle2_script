// last update: 2021-10-13
// for Belle2 data

#include "constants.h"
#include "ObtainWeight.h"
#include "correctors.h"
#include "base.h"

# define NgammaMAX 20

# define MCTYPE "MC15rd"

# define CAL 1.0 // this should be 1.0
# define CAL_qq 1.0

Corrector corrector;
Corrector_FEI corrector_FEI;
Corrector_PID corrector_PID;
Corrector_pi0 corrector_pi0;
Corrector_FakePID corrector_FakePID;
Corrector_Knn corrector_Knn;
Corrector_Xsnn corrector_Xsnn;
Corrector_KpKLKL corrector_KpKLKL;
Corrector_KSKLKL corrector_KSKLKL;
Corrector_phiKL corrector_phiKL;
Corrector_KstarKLKL corrector_KstarKLKL;
Corrector_XsKLKL corrector_XsKLKL;
Corrector_BtoDtoXKL corrector_BtoDtoXKL;
Corrector_KS0 corrector_KS0;

Corrector_Multiplicity* corrector_Multiplicity;

/* ====================================== */

void LetsFillNgamma(const char* dirname, TH1D* hist_Ngamma, std::string SampleName, int option, bool IsMultiplicityCorrectionApplied) {
    /*
    SampleName for Knn
    CHG
    MIX
    UUBAR
    DDBAR
    SSBAR
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

    double var = -1;

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

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &var);

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID); // charged: 0, mixed: 1
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
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

            if (option == 1 && Upsilon_ID != 0) continue;
            else if (option == 2 && Upsilon_ID != 1) continue;

            // Fill numberings
            double weight_ri = 0.0;
            if (SampleName == "CHG") {
                FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "MIX") {
                FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "UUBAR") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "DDBAR") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "SSBAR") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "CHARM") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            //else if (job_id >= 256846858 && job_id <= 256847295) numberings->push_back(6);
            //else if (job_id >= 256847296 && job_id <= 256847807) numberings->push_back(7);
            //else if (job_id >= 256847808 && job_id <= 256848291) numberings->push_back(8);
            //else if (job_id >= 256848292 && job_id <= 256848743) numberings->push_back(9);
            //else if (job_id >= 256848744 && job_id <= 256849128) numberings->push_back(10);
            //else if (job_id >= 256849129 && job_id <= 256849396) numberings->push_back(11);
            else {
                printf("undefined job id!\n");
                exit(1);
            }

            // Fill calibration factors
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

            // Multiplicity correction factor, it is not applied now. it is for a systematic uncertainty
            double Ngamma_v200 = var;
            double Correction_multiplicity = 1.0;
            if (IsMultiplicityCorrectionApplied) Correction_multiplicity = corrector_Multiplicity->GetCorrectionFactor(Ngamma_v200);

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

            double weight = 1.0;
            if (IsMultiplicityCorrectionApplied) weight = FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_KS0 * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL;
            else weight = FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_KS0 * Correction_Knn * Correction_Xnn * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL;

            hist_Ngamma->Fill(var, weight);


        }
        input_file->Close();

    }

}

void LetsFillNgamma(const char* dirname, TH1D* hist_Ngamma, int option) {
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

    double var = -1;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &var);

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID); // charged: 0, mixed: 1
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            if (option == 1 && Upsilon_ID != 0) continue;
            else if (option == 2 && Upsilon_ID != 1) continue;
            hist_Ngamma->Fill(var);

        }
        input_file->Close();

    }

}

void MultiplicityCalculator(){

    
    const char* ChargeSideband_MC_CHG_validation_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD_LS_MC_cside/CHG_analysis/validation_v009/final_output";
    const char* ChargeSideband_MC_MIX_validation_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD_LS_MC_cside/MIX_analysis/validation_v009/final_output";
    const char* ChargeSideband_MC_UUBAR_validation_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD_LS_MC_cside/UUBAR_analysis/validation_v009/final_output";
    const char* ChargeSideband_MC_DDBAR_validation_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD_LS_MC_cside/DDBAR_analysis/validation_v009/final_output";
    const char* ChargeSideband_MC_SSBAR_validation_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD_LS_MC_cside/SSBAR_analysis/validation_v009/final_output";
    const char* ChargeSideband_MC_CHARM_validation_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD_LS_MC_cside/CHARM_analysis/validation_v009/final_output";

    const char* ChargeSideband_data_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD_LS_data_cside/SIGNAL_analysis/validation_v009/final_output";

    // Lets fill!
    TH1D* Ngamma_v200_MC = new TH1D("Ngamma_v200_MC", ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);
    TH1D* Ngamma_v200_data = new TH1D("Ngamma_v200_data", ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);

    LetsFillNgamma(ChargeSideband_MC_CHG_validation_dirname, Ngamma_v200_MC, "CHG", 1, false);
    LetsFillNgamma(ChargeSideband_MC_MIX_validation_dirname, Ngamma_v200_MC, "MIX", 1, false);
    LetsFillNgamma(ChargeSideband_MC_UUBAR_validation_dirname, Ngamma_v200_MC, "UUBAR", 1, false);
    LetsFillNgamma(ChargeSideband_MC_DDBAR_validation_dirname, Ngamma_v200_MC, "DDBAR", 1, false);
    LetsFillNgamma(ChargeSideband_MC_SSBAR_validation_dirname, Ngamma_v200_MC, "SSBAR", 1, false);
    LetsFillNgamma(ChargeSideband_MC_CHARM_validation_dirname, Ngamma_v200_MC, "CHARM", 1, false);

    LetsFillNgamma(ChargeSideband_data_dirname, Ngamma_v200_data, 1);

    // print weight file
    FILE* fp = fopen("multiplicity_weight.txt","w");
    fprintf(fp, "%d\n", NgammaMAX);
    for (int i = 0; i < NgammaMAX + 1; i++) {
        double MC_num_bin = Ngamma_v200_MC->GetBinContent(i + 1);
        double data_num_bin = Ngamma_v200_data->GetBinContent(i + 1);
        if (MC_num_bin > MyEPSILON) fprintf(fp, "%lf %lf %lf\n", data_num_bin, MC_num_bin, data_num_bin / MC_num_bin);
        else fprintf(fp, "%lf %lf %lf\n", data_num_bin, MC_num_bin, 1.0);
    }
    fclose(fp);

    double MC_num = 0;
    double data_num = 0;

    for (int i = 0; i < NgammaMAX + 1; i++) {
        MC_num = MC_num + Ngamma_v200_MC->GetBinContent(i + 1);
        data_num = data_num + Ngamma_v200_data->GetBinContent(i + 1);
    }

    printf("option1 data num: %lf\n", data_num);
    printf("option1 MC num with calibration: %lf\n", MC_num);




    corrector_Multiplicity = new Corrector_Multiplicity();

    // Lets fill!
    TH1D* Ngamma_v200_MC_test = new TH1D("Ngamma_v200_MC_test", ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);
    TH1D* Ngamma_v200_data_test = new TH1D("Ngamma_v200_data_test", ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);

    LetsFillNgamma(ChargeSideband_MC_CHG_validation_dirname, Ngamma_v200_MC_test, "CHG", 2, true);
    LetsFillNgamma(ChargeSideband_MC_MIX_validation_dirname, Ngamma_v200_MC_test, "MIX", 2, true);
    LetsFillNgamma(ChargeSideband_MC_UUBAR_validation_dirname, Ngamma_v200_MC_test, "UUBAR", 2, true);
    LetsFillNgamma(ChargeSideband_MC_DDBAR_validation_dirname, Ngamma_v200_MC_test, "DDBAR", 2, true);
    LetsFillNgamma(ChargeSideband_MC_SSBAR_validation_dirname, Ngamma_v200_MC_test, "SSBAR", 2, true);
    LetsFillNgamma(ChargeSideband_MC_CHARM_validation_dirname, Ngamma_v200_MC_test, "CHARM", 2, true);

    LetsFillNgamma(ChargeSideband_data_dirname, Ngamma_v200_data_test, 2);

    // print weight file
    fp = fopen("multiplicity_weight_uncertainty.txt", "w");
    fprintf(fp, "%d\n", NgammaMAX);
    for (int i = 0; i < NgammaMAX + 1; i++) {
        double MC_num_bin = Ngamma_v200_MC_test->GetBinContent(i + 1);
        double data_num_bin = Ngamma_v200_data_test->GetBinContent(i + 1);
        if (MC_num_bin > MyEPSILON) fprintf(fp, "%lf %lf %lf\n", data_num_bin, MC_num_bin, data_num_bin / MC_num_bin);
        else fprintf(fp, "%lf %lf %lf\n", data_num_bin, MC_num_bin, 1.0);
    }
    fclose(fp);

    MC_num = 0;
    data_num = 0;

    for (int i = 0; i < NgammaMAX + 1; i++) {
        MC_num = MC_num + Ngamma_v200_MC_test->GetBinContent(i + 1);
        data_num = data_num + Ngamma_v200_data_test->GetBinContent(i + 1);
    }

    printf("option2 data num: %lf\n", data_num);
    printf("option2 MC num with calibration: %lf\n", MC_num);
}
