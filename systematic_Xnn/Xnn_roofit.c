// last update: 2021-10-13
// for Belle2 data

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"

#include "constants.h"
#include "ObtainWeight.h"
#include "correctors.h"
#include "base.h"

using namespace RooFit;

# define MCTYPE "MC15rd"

# define CAL 1.0
# define CAL_qq 1.0

Corrector_FEI corrector_FEI;
Corrector_PID corrector_PID;
Corrector_pi0 corrector_pi0;
Corrector_FakePID corrector_FakePID;
Corrector_Knn corrector_Knn;
//Corrector_Multiplicity corrector_Multiplicity;
Corrector_KpKLKL corrector_KpKLKL;
Corrector_KSKLKL corrector_KSKLKL;
Corrector_BtoDtoXKL corrector_BtoDtoXKL;
Corrector_KS0 corrector_KS0;

double LetsAdd(const char* dirname, const char* including_string, RooRealVar* Mbc_, RooRealVar* weight_, RooDataSet* info_, std::string SampleName, double additional_weight = 1.0) {
    double Nevt = 0;
    double Mbc_var = 0;

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

    double nDptoXKL = -1;
    double nD0toXKL = -1;

    int nBp = -1;
    int nB0 = -1;

    double Bsig_M = -1;

    std::vector<string> names;
    load_files(dirname, &names, including_string);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");
        TTree* tree_Xs;
        if (SampleName == "SIGNAL") tree_Xs = (TTree*)input_file->Get("Xs");

        tree_Bsig->SetBranchAddress("Bsig_Mbc", &Mbc_var); // Mbc

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
        for (int i_pi0 = 0; i_pi0 < N_pi0_syst_MC15ri; i_pi0++) tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npi0MC15ribin" + std::to_string(i_pi0)).c_str(), &temp_N_pi0_syst_MC15ri[i_pi0]);
        for (int i_pi0 = 0; i_pi0 < N_pi0_syst_MC15rd; i_pi0++) tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npi0MC15rdbin" + std::to_string(i_pi0)).c_str(), &temp_N_pi0_syst_MC15rd[i_pi0]);
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
            if (SampleName == "SIGNAL") tree_Xs->GetEntry(j);

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
            else if (SampleName == "SIGNAL") {
                if (nBp > 0) {
                    FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
                    weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
                }
                else if (nB0 > 0) {
                    FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
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
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Multiplicity correction factor, it is not applied now. it is for a systematic uncertainty
            //double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // B+ --> K+ KL0 KL0 correction factor
            double Correction_KpKLKL = corrector_KpKLKL.GetCorrectionFactorAtGeneric(s13_KpKLKL, s23_KpKLKL, nB2KpKLKL_all_KpKLKL, nB2KpKLKL_NR_KpKLKL);

            // B0 --> KS0 KL0 KL0 correction factor
            double Correction_KSKLKL = corrector_KSKLKL.GetCorrectionFactorAtGeneric(std::max(std::max(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), std::min(std::min(s13_KSKLKL, s23_KSKLKL), s12_KSKLKL), nB2KSKLKL_all_KSKLKL, nB2KSKLKL_NR_KSKLKL);

            // B-> [D -> KL0 X] anything correction factor
            double Correction_BtoDtoXKL = 1.0;
            if (SampleName == "CHG" || SampleName == "MIX" || SampleName == "SIGNAL") Correction_BtoDtoXKL = corrector_BtoDtoXKL.GetCorrectionFactorAtGeneric(nDptoXKL + nD0toXKL);

            double total_weight = FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_KS0 * Correction_Knn * Correction_KpKLKL * Correction_KSKLKL * Correction_BtoDtoXKL * additional_weight;

            *Mbc_ = Mbc_var;
            *weight_ = total_weight;
            info_->add(RooArgSet(*Mbc_, *weight_), weight_->getVal());
            Nevt = Nevt + total_weight;
        }
        input_file->Close();

    }

    return Nevt;
}

double LetsAdd(const char* dirname, RooRealVar* Mbc_, RooRealVar* weight_, RooDataSet* info_) {
    double Nevt = 0;
    double Mbc_var = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_Bsig->SetBranchAddress("Bsig_Mbc", &Mbc_var); // Mbc

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            double total_weight = 1.0;

            *Mbc_ = Mbc_var;
            *weight_ = total_weight;
            info_->add(RooArgSet(*Mbc_, *weight_), weight_->getVal());
            Nevt = Nevt + total_weight;
        }
        input_file->Close();

    }

    return Nevt;
}

void Xnn_roofit() {

    // Observable:
    RooRealVar Mbc("Mbc", "Mbc", 5.2, 5.3);
    RooRealVar weight("weight", "weight", 0.0, 10.0);
    RooDataSet info("info", "info", RooArgSet(Mbc, weight), WeightVar("weight"));

    const char* dirname = "./";
    LetsAdd(dirname, "CHG",  & Mbc, &weight, &info, "CHG", 1.0);
    LetsAdd(dirname, "MIX", &Mbc, &weight, &info, "MIX", 1.0);
    LetsAdd(dirname, "UUBAR", &Mbc, &weight, &info, "UUBAR", 1.0);
    LetsAdd(dirname, "DDBAR", &Mbc, &weight, &info, "DDBAR", 1.0);
    LetsAdd(dirname, "SSBAR", &Mbc, &weight, &info, "SSBAR", 1.0);
    LetsAdd(dirname, "CHARM", &Mbc, &weight, &info, "CHARM", 1.0);
    RooPlot* Mbcframe = Mbc.frame(Bins(50), Title("M_{bc}^{sig} fit"));
    RooDataSet* d_Mbc = (RooDataSet*)info.reduce(RooArgSet(Mbc));

    // ARGUS background
    RooRealVar End("End", "End of Argus", 5.29); // fixed by experimental setting
    RooRealVar c("curvature", "curvature of ARGUS", -68, -100, 20);
    RooRealVar power("power", "power of ARGUS", 0.5); // just fixed

    RooArgusBG pdfARGUS("ARGUS", "ARGUS PDF", Mbc, End, c, power);

    // crystal ball signal
    RooRealVar mean("mean", "mean of Crystal ball", 5.27981); // fixed by MC
    RooRealVar width("width", "width of Crystal ball", 0.003151); // fixed by MC
    RooRealVar alpha("alpha", "alpha of Crystal ball", 1.327851); // fixed by MC
    RooRealVar n("n", "n of Crystal ball", 0.9, 0.5, 1.5);

    RooCrystalBall pdfCRYSTAL("CRYSTAL", "CRYSTAL PDF", Mbc, mean, width, alpha, n);

    // Construct a signal and background PDF:
    RooRealVar nsig("nsig", "nsig", 20, 0, 50);
    RooRealVar nbkg("nbkg", "nkbkg", 300, 200, 500);
    RooExtendPdf esig("esignal", "extended signal p.d.f", pdfCRYSTAL, nsig);
    RooExtendPdf ebkg("ebkg", "extended background p.d.f", pdfARGUS, nbkg);

    RooAddPdf  model("model", "b+n", RooArgList(ebkg, esig));

    // fit
    RooFitResult* result = model.fitTo(*d_Mbc, Extended(), Save(), SumW2Error(true));

    // print fit result
    RooArgSet fitargs_LT = result->floatParsFinal();
    TIterator* iter_LT(fitargs_LT.createIterator());

    for (TObject* a_LT = iter_LT->Next(); a_LT != 0; a_LT = iter_LT->Next()) {
        RooRealVar* rrv_LT = dynamic_cast<RooRealVar*>(a_LT);
        std::string name_LT = rrv_LT->GetName();
        double val_LT = rrv_LT->getVal();
        double err_LT = rrv_LT->getError();
        printf("%s: %lf +- %lf\n", name_LT.c_str(), val_LT, err_LT);
    }

    // Draw result
    d_Mbc->plotOn(Mbcframe, DataError(RooAbsData::SumW2));
    model.plotOn(Mbcframe, LineColor(kRed), Normalization(1.0, RooAbsReal::RelativeExpected));
    model.plotOn(Mbcframe, Components(pdfARGUS), LineColor(kBlue), LineStyle(kDashed), Normalization(1.0, RooAbsReal::RelativeExpected));
    model.plotOn(Mbcframe, Components(pdfCRYSTAL), LineColor(kViolet), LineStyle(kDashed), Normalization(1.0, RooAbsReal::RelativeExpected));

    TCanvas* c_temp = new TCanvas("Mbc", "Mbc", 600, 600);
    gPad->SetLeftMargin(0.15); Mbcframe->GetYaxis()->SetTitleOffset(1.4); Mbcframe->Draw();
    c_temp->SaveAs("fit_result.png");
}
