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

# include <algorithm>
# include <float.h>

#include "constants.h"
#include "ObtainWeight.h"
#include "correctors.h"
#include "base.h"
#include "template.h"
#include "Fitter.h"

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

# define MCTYPE "MC15rd"

# define KS0_rel_uncertainty 0.6 // %/cm
# define track_rel_uncertainty 0.24 // %
# define Kaon_PID_max_uncertainty 0.1 // not percentage. relative uncertainty

int Nvar_num = -1;

double CAL = 1.0; // must be 1.0
# define CAL_qq_K 1.0
# define CAL_qq_Kstar 1.0
# define CAL_qq_Xs 1.0

bool NormalizeAtEachMXs = false;

bool DrawFitPlots = false;

std::vector<double> Signal_MXs1_Nevts_nominal;
std::vector<double> Signal_MXs2_Nevts_nominal;
std::vector<double> Signal_MXs3_Nevts_nominal;
std::vector<double> CHG_Nevts_nominal;
std::vector<double> MIX_Nevts_nominal;
std::vector<double> UUBAR_Nevts_nominal;
std::vector<double> DDBAR_Nevts_nominal;
std::vector<double> SSBAR_Nevts_nominal;
std::vector<double> CHARM_Nevts_nominal;

std::vector<double> Signal_MXs1_Nevts_fit;
std::vector<double> Signal_MXs2_Nevts_fit;
std::vector<double> Signal_MXs3_Nevts_fit;
std::vector<double> CHG_Nevts_fit;
std::vector<double> MIX_Nevts_fit;
std::vector<double> UUBAR_Nevts_fit;
std::vector<double> DDBAR_Nevts_fit;
std::vector<double> SSBAR_Nevts_fit;
std::vector<double> CHARM_Nevts_fit;

/* ====================================== */

// Histogram for counting
TH1D* MC_one_bin = new TH1D("MC_one_bin", ";number of candidates", 1, -100, 100);
TH1D* data_one_bin = new TH1D("data_one_bin", ";number of candidates", 1, -100, 100);
TH1D* Ratio_one_bin = new TH1D("ratio_one_bin", ";number of candidates", 1, -100, 100);

TH1D* MC_Nevt_MXs1 = new TH1D("MC_Nevt_MXs1", ";number of candidates", 1, -100, 100);
TH1D* MC_Nevt_MXs2 = new TH1D("MC_Nevt_MXs2", ";number of candidates", 1, -100, 100);
TH1D* MC_Nevt_MXs3 = new TH1D("MC_Nevt_MXs3", ";number of candidates", 1, -100, 100);
TH1D* data_Nevt_MXs1 = new TH1D("data_Nevt_MXs1", ";number of candidates", 1, -100, 100);
TH1D* data_Nevt_MXs2 = new TH1D("data_Nevt_MXs2", ";number of candidates", 1, -100, 100);
TH1D* data_Nevt_MXs3 = new TH1D("data_Nevt_MXs3", ";number of candidates", 1, -100, 100);
TH1D* Ratio_Nevt_MXs1 = new TH1D("Ratio_Nevt_MXs1", ";number of candidates", 1, -100, 100);
TH1D* Ratio_Nevt_MXs2 = new TH1D("Ratio_Nevt_MXs2", ";number of candidates", 1, -100, 100);
TH1D* Ratio_Nevt_MXs3 = new TH1D("Ratio_Nevt_MXs3", ";number of candidates", 1, -100, 100);

enum class MCsample {
    CHG,
    MIX,
    UUBAR,
    DDBAR,
    SSBAR,
    CHARM,
    TAU,
    MUMU,
    GG,
    EE,
    EEEE,
    EEMUMU,
    LLXX,
    HHISR,
    SIGNAL
};

double BDTcToWeight(double BDTc) {

    if (BDTc > (5.0 / 6.0)) return std::sqrt(5.0);
    else return std::sqrt(BDTc / (1.0 - BDTc));

}

void ReadNominalNevt(const char* filename) {
    FILE* fp = fopen(filename, "r");

    fscanf(fp, "Signal_MXs1:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        Signal_MXs1_Nevts_nominal.push_back(temp);
    }
    fscanf(fp, "\n");

    fscanf(fp, "Signal_MXs2:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        Signal_MXs2_Nevts_nominal.push_back(temp);
    }
    fscanf(fp, "\n");

    fscanf("Signal_MXs3:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        Signal_MXs3_Nevts_nominal.push_back(temp);
    }
    fscanf(fp, "\n");

    fscanf(fp, "CHG:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        CHG_Nevts_nominal.push_back(temp);
    }
    fscanf(fp, "\n");

    fscanf(fp, "MIX:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        MIX_Nevts_nominal.push_back(temp);
    }
    fscanf(fp, "\n");

    fscanf(fp, "UUBAR:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        UUBAR_Nevts_nominal.push_back(temp);
    }
    fscanf("\n");

    fscanf(fp, "DDBAR:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        DDBAR_Nevts_nominal.push_back(temp);
    }
    fscanf("\n");

    fscanf(fp, "SSBAR:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        SSBAR_Nevts_nominal.push_back(temp);
    }
    fscanf(fp, "\n");

    fscanf(fp, "CHARM:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        CHARM_Nevts_nominal.push_back(temp);
    }
    fscanf(fp, "\n");

    fclose(fp);
}

void ReadFittedNevt(const char* filename) {
    FILE* fp = fopen(filename, "r");

    fscanf(fp, "Signal_MXs1:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        Signal_MXs1_Nevts_fit.push_back(temp);
    }
    fscanf(fp, "\n");

    fscanf(fp, "Signal_MXs2:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        Signal_MXs2_Nevts_fit.push_back(temp);
    }
    fscanf(fp, "\n");

    fscanf(fp, "Signal_MXs3:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        Signal_MXs3_Nevts_fit.push_back(temp);
    }
    fscanf("\n");

    fscanf(fp, "CHG:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        CHG_Nevts_fit.push_back(temp);
    }
    fscanf(fp, "\n");

    fscanf(fp, "MIX:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        MIX_Nevts_fit.push_back(temp);
    }
    fscanf(fp, "\n");

    fscanf(fp, "UUBAR:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        UUBAR_Nevts_fit.push_back(temp);
    }
    fscanf(fp, "\n");

    fscanf(fp, "DDBAR:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        DDBAR_Nevts_fit.push_back(temp);
    }
    fscanf(fp, "\n");

    fscanf(fp, "SSBAR:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        SSBAR_Nevts_fit.push_back(temp);
    }
    fscanf(fp, "\n");

    fscanf(fp, "CHARM:\n");
    for (int i = 0; i < RarityBins; i++) {
        double temp;
        fscanf(fp, "%lf ", &temp);
        CHARM_Nevts_fit.push_back(temp);
    }
    fscanf(fp, "\n");

    fclose(fp);
}

void LetsDrawFitPlot(bool DrawFitPlots_) {

    if (DrawFitPlots_ == true) {
        if (Signal_MXs1_Nevts_nominal.size() != RarityBins) {
            printf("[ERROR] You should call `ReadNominalNevt` first\n");
            exit(1);
        }

        if (Signal_MXs1_Nevts_fit.size() != RarityBins) {
            printf("[ERROR] You should call `ReadFittedNevt` first\n");
            exit(1);
        }
    }

    DrawFitPlots = DrawFitPlots_;
}

double GetFitWeight(double bin_index, double MC_MXs, std::string SampleName) {
    // Currently, I only include this function into `NevtCount_ri`, `LetsFillMC`, and `LetsFillMC_correction`

    // double bin_index starts from 0.5
    int bin_index_int = static_cast <int> (std::floor(bin_index));

    if (SampleName == "CHG") {
        if (CHG_Nevts_nominal.at(bin_index_int) < MyEPSILON) return 1.0;
        else return CHG_Nevts_fit.at(bin_index_int) / CHG_Nevts_nominal.at(bin_index_int);
    }
    else if (SampleName == "MIX") {
        if (MIX_Nevts_nominal.at(bin_index_int) < MyEPSILON) return 1.0;
        else return MIX_Nevts_fit.at(bin_index_int) / MIX_Nevts_nominal.at(bin_index_int);
    }
    else if (SampleName == "UUBAR") {
        if (UUBAR_Nevts_nominal.at(bin_index_int) < MyEPSILON) return 1.0;
        else return UUBAR_Nevts_fit.at(bin_index_int) / UUBAR_Nevts_nominal.at(bin_index_int);
    }
    else if (SampleName == "DDBAR") {
        if (DDBAR_Nevts_nominal.at(bin_index_int) < MyEPSILON) return 1.0;
        else return DDBAR_Nevts_fit.at(bin_index_int) / DDBAR_Nevts_nominal.at(bin_index_int);
    }
    else if (SampleName == "SSBAR") {
        if (SSBAR_Nevts_nominal.at(bin_index_int) < MyEPSILON) return 1.0;
        else return SSBAR_Nevts_fit.at(bin_index_int) / SSBAR_Nevts_nominal.at(bin_index_int);
    }
    else if (SampleName == "CHARM") {
        if (CHARM_Nevts_nominal.at(bin_index_int) < MyEPSILON) return 1.0;
        else return CHARM_Nevts_fit.at(bin_index_int) / CHARM_Nevts_nominal.at(bin_index_int);
    }
    else if (SampleName == "SIGNAL") {

        if ((MC_MXs > 0.0) && (MC_MXs < 0.6)) {
            if (Signal_MXs1_Nevts_nominal.at(bin_index_int) < MyEPSILON) return 1.0;
            else return Signal_MXs1_Nevts_fit.at(bin_index_int) / Signal_MXs1_Nevts_nominal.at(bin_index_int);
        }
        else if ((MC_MXs >= 0.6) && (MC_MXs < 1.0)) {
            if (Signal_MXs2_Nevts_nominal.at(bin_index_int) < MyEPSILON) return 1.0;
            else return Signal_MXs2_Nevts_fit.at(bin_index_int) / Signal_MXs2_Nevts_nominal.at(bin_index_int);
        }
        else if (MC_MXs >= 1.0) {
            if (Signal_MXs3_Nevts_nominal.at(bin_index_int) < MyEPSILON) return 1.0;
            else return Signal_MXs3_Nevts_fit.at(bin_index_int) / Signal_MXs3_Nevts_nominal.at(bin_index_int);
        }
        else {
            printf("[ERROR] unexpected MC_MXs for SIGNAL\n");
            exit(1);
        }

    }

    // never reached
    return -1.0;
}

void LetsFillMC(const char* dirname, std::vector<std::string> variable_names, std::vector<std::string> branch_names, std::vector<double> variable_values[], std::vector<int>* numberings, std::vector<double>* weights, std::string SampleName, int option = 0, double additional_weight = 1.0) {
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

    double* var = (double*)malloc(sizeof(double) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var[i] = 0.0;
    float* var_float = (float*)malloc(sizeof(float) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var_float[i] = 0.0;
    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    double FEI_calibration_factor = -1;

    double MC_MXs = -1;
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

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

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");
        TTree* tree_Xs;
        if (SampleName == "SIGNAL") tree_Xs = (TTree*)input_file->Get("Xs");

        for (int k = 0; k < (int)variable_names.size(); k++) {
            if (variable_names.at(k).find("bin index") != std::string::npos) continue;

            if (branch_names.at(k) == std::string("Upsilon")) {
                if(variable_names.at(k).find("MVA") == std::string::npos) tree_upsilon->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
                else tree_upsilon->SetBranchAddress(variable_names.at(k).c_str(), &var_float[k]);
            }
            else if (branch_names.at(k) == std::string("Bsig")) tree_Bsig->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else if (branch_names.at(k) == std::string("Btag")) tree_Btag->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else {
                printf("ERROR! \n");
                exit(1);
            }
        }

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

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
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

        Ngamma_v200_index = std::find(variable_names.begin(), variable_names.end(), std::string("extraInfo__boNgammav200__bc")) - variable_names.begin();
        if(Ngamma_v200_index == variable_names.size()) tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

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

            if (option == 1 && Upsilon_ID != 0) continue;
            else if (option == 2 && Upsilon_ID != 1) continue;
            else if (option == 0) {}
            else if ((option != 0) && (option != 1) && (option != 2)){
                printf("improper option value!\n");
                exit(1);
            }

            // obtain true MXs
            if (SampleName == "SIGNAL") {

                // set MXs from filename
                if (std::strstr(names.at(i).c_str(), "B2Knunu") != nullptr) MC_MXs = Mxs_Bc_MC;
                else if (std::strstr(names.at(i).c_str(), "B2Kstarnunu") != nullptr) MC_MXs = Mxs_Bc_MC;
                else if (std::strstr(names.at(i).c_str(), "B2Xsnunu") != nullptr) MC_MXs = Mxs_Bc_MC;
                else if (std::strstr(names.at(i).c_str(), "B02K0nunu") != nullptr) MC_MXs = Mxs_B0_MC;
                else if (std::strstr(names.at(i).c_str(), "B02Kstar0nunu") != nullptr) MC_MXs = Mxs_B0_MC;
                else if (std::strstr(names.at(i).c_str(), "B02Xsnunu") != nullptr) MC_MXs = Mxs_B0_MC;
                else {
                    printf("[ERROR] unexpected filename for SIGNAL\n");
                    exit(1);
                }

                // sanity check
                if ((MC_MXs > 0.0) && (MC_MXs < 6.0)) {}
                else {
                    if (std::strstr(names.at(i).c_str(), "B2Knunu") != nullptr) MC_MXs = 0.4868;
                    else if (std::strstr(names.at(i).c_str(), "B2Kstarnunu") != nullptr) MC_MXs = 0.8916;
                    else if (std::strstr(names.at(i).c_str(), "B2Xsnunu") != nullptr) MC_MXs = 1.5;
                    else if (std::strstr(names.at(i).c_str(), "B02K0nunu") != nullptr) MC_MXs = 0.4868;
                    else if (std::strstr(names.at(i).c_str(), "B02Kstar0nunu") != nullptr) MC_MXs = 0.8916;
                    else if (std::strstr(names.at(i).c_str(), "B02Xsnunu") != nullptr) MC_MXs = 1.5;
                }

            }
            else {
                MC_MXs = -1;
                Mxs_Bc_MC = -1;
                Mxs_B0_MC = -1;
            }

            for (int k = 0; k < (int)variable_names.size(); k++) {
                if (variable_names.at(k).find("bin index") != std::string::npos) {
                    auto it = std::find(variable_names.begin(), variable_names.end(), "MVA_BB");
                    int index_FBDT_raw = std::distance(variable_names.begin(), it);
                    variable_values[k].push_back(ReturnBinIndex(var_float[index_FBDT_raw], Bsig_M));
                }
                else if (variable_names.at(k).find("MVA") == std::string::npos) variable_values[k].push_back(var[k]);
                else variable_values[k].push_back((double)var_float[k]);
                //else variable_values[k].push_back(GetBinIndex((double)var_float[k], Bsig_M));
            }

            // Fill numberings
            double weight_ri = 0.0;
            if (SampleName == "CHG") {
                numberings->push_back(static_cast<int>(MCsample::CHG));
                FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "MIX") {
                numberings->push_back(static_cast<int>(MCsample::MIX));
                FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "UUBAR") {
                numberings->push_back(static_cast<int>(MCsample::UUBAR));
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) FEI_calibration_factor = CAL_qq_K;
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) FEI_calibration_factor = CAL_qq_Kstar;
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) FEI_calibration_factor = CAL_qq_Xs;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "DDBAR") {
                numberings->push_back(static_cast<int>(MCsample::DDBAR));
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) FEI_calibration_factor = CAL_qq_K;
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) FEI_calibration_factor = CAL_qq_Kstar;
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) FEI_calibration_factor = CAL_qq_Xs;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "SSBAR") {
                numberings->push_back(static_cast<int>(MCsample::SSBAR));
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) FEI_calibration_factor = CAL_qq_K;
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) FEI_calibration_factor = CAL_qq_Kstar;
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) FEI_calibration_factor = CAL_qq_Xs;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "CHARM") {
                numberings->push_back(static_cast<int>(MCsample::CHARM));
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) FEI_calibration_factor = CAL_qq_K;
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) FEI_calibration_factor = CAL_qq_Kstar;
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) FEI_calibration_factor = CAL_qq_Xs;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "SIGNAL") {
                numberings->push_back(static_cast<int>(MCsample::SIGNAL));
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
            if (Ngamma_v200_index != variable_names.size()) Ngamma_v200 = var[Ngamma_v200_index];
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

            // norm factor for each MXs
            double normfactor_MXs =  1.0;
            if (NormalizeAtEachMXs == true) {
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) normfactor_MXs = data_Nevt_MXs1->GetBinContent(1) / MC_Nevt_MXs1->GetBinContent(1);
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) normfactor_MXs = data_Nevt_MXs2->GetBinContent(1) / MC_Nevt_MXs2->GetBinContent(1);
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) normfactor_MXs = data_Nevt_MXs3->GetBinContent(1) / MC_Nevt_MXs3->GetBinContent(1);
            }

            // scale factor to adjust into fit result
            double fitfactor = 1.0;
            if (DrawFitPlots == true) {
                auto it = std::find(variable_names.begin(), variable_names.end(), "MVA_BB");
                int index_FBDT_raw = std::distance(variable_names.begin(), it);
                double bin_index = ReturnBinIndex(var_float[index_FBDT_raw], Bsig_M);
                fitfactor = GetFitWeight(bin_index, MC_MXs, SampleName);
            }

            weights->push_back(FEI_calibration_factor* CAL* weight_ri* Correction_pi0* Correction_KID* Correction_PID* Correction_fake* Correction_Knn* Correction_Xnn* Correction_multiplicity* Correction_KpKLKL* Correction_KSKLKL* Correction_KstarKLKL* Correction_XKLKL* Correction_BtoDtoXKL* additional_weight * normfactor_MXs * fitfactor);

            MC_one_bin->Fill(1.0, FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * additional_weight * normfactor_MXs * fitfactor);
        }
        input_file->Close();

    }

    free(var);
    free(var_float);

}

void LetsFilldata(const char* dirname, std::vector<std::string> variable_names, std::vector<std::string> branch_names, std::vector<double> variable_values[], const char* included_string = "root", int option = 0) {
    /*
    option 0: select all Btag
    option 1: select Btag+
    option 2: select Btag0
    */
    double* var = (double*)malloc(sizeof(double) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var[i] = 0.0;
    float* var_float = (float*)malloc(sizeof(float) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var_float[i] = 0.0;
    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;

    double Bsig_M = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        for (int k = 0; k < (int)variable_names.size(); k++) {
            if (variable_names.at(k).find("bin index") != std::string::npos) continue;

            if (branch_names.at(k) == std::string("Upsilon")) {
                if (variable_names.at(k).find("MVA") == std::string::npos) tree_upsilon->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
                else tree_upsilon->SetBranchAddress(variable_names.at(k).c_str(), &var_float[k]);
            }
            else if (branch_names.at(k) == std::string("Bsig")) tree_Bsig->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else if (branch_names.at(k) == std::string("Btag")) tree_Btag->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else {
                printf("ERROR! \n");
                exit(1);
            }
        }

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID); // charged: 0, mixed: 1
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            if (option == 1 && Upsilon_ID != 0) continue;
            else if (option == 2 && Upsilon_ID != 1) continue;
            else if (option == 0) {}
            else if ((option != 0) && (option != 1) && (option != 2)) {
                printf("improper option value!\n");
                exit(1);
            }

            for (int k = 0; k < (int)variable_names.size(); k++) {
                if (variable_names.at(k).find("bin index") != std::string::npos) {
                    auto it = std::find(variable_names.begin(), variable_names.end(), "MVA_BB");
                    int index_FBDT_raw = std::distance(variable_names.begin(), it);
                    variable_values[k].push_back(ReturnBinIndex(var_float[index_FBDT_raw], Bsig_M));
                }
                else if (variable_names.at(k).find("MVA") == std::string::npos) variable_values[k].push_back(var[k]);
                else variable_values[k].push_back((double)var_float[k]);
                //else variable_values[k].push_back(GetBinIndex((double)var_float[k], Bsig_M));
            }

            data_one_bin->Fill(1.0);
        }
        input_file->Close();

    }

    free(var);
    free(var_float);

}

void LetsFillembeddedMC(const char* dirname, std::vector<std::string> variable_names, std::vector<std::string> branch_names, std::vector<double> variable_values[], std::vector<double>* weights, std::string SampleName, int option = 0, double additional_weight = 1.0) {
    /*
    option 0: select all Btag
    option 1: select Btag+
    option 2: select Btag0
    */
    double* var = (double*)malloc(sizeof(double) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var[i] = 0.0;
    float* var_float = (float*)malloc(sizeof(float) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var_float[i] = 0.0;
    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;

    double FEI_calibration_factor = -1;

    int Ngamma_v200_index = -1;
    double Ngamma_v200 = -1;

    double Bsig_M = -1;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        for (int k = 0; k < (int)variable_names.size(); k++) {
            if (variable_names.at(k).find("bin index") != std::string::npos) continue;

            if (branch_names.at(k) == std::string("Upsilon")) {
                if (variable_names.at(k).find("MVA") == std::string::npos) tree_upsilon->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
                else tree_upsilon->SetBranchAddress(variable_names.at(k).c_str(), &var_float[k]);
            }
            else if (branch_names.at(k) == std::string("Bsig")) tree_Bsig->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else if (branch_names.at(k) == std::string("Btag")) tree_Btag->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else {
                printf("ERROR! \n");
                exit(1);
            }
        }

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID); // charged: 0, mixed: 1
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        tree_Bsig->SetBranchAddress("Bsig_M", &Bsig_M);

        Ngamma_v200_index = std::find(variable_names.begin(), variable_names.end(), std::string("extraInfo__boNgammav200__bc")) - variable_names.begin();
        if (Ngamma_v200_index == variable_names.size()) tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            if (option == 1 && Upsilon_ID != 0) continue;
            else if (option == 2 && Upsilon_ID != 1) continue;
            else if (option == 0) {}
            else if ((option != 0) && (option != 1) && (option != 2)) {
                printf("improper option value!\n");
                exit(1);
            }

            for (int k = 0; k < (int)variable_names.size(); k++) {
                if (variable_names.at(k).find("bin index") != std::string::npos) {
                    auto it = std::find(variable_names.begin(), variable_names.end(), "MVA_BB");
                    int index_FBDT_raw = std::distance(variable_names.begin(), it);
                    variable_values[k].push_back(ReturnBinIndex(var_float[index_FBDT_raw], Bsig_M));
                }
                else if (variable_names.at(k).find("MVA") == std::string::npos) variable_values[k].push_back(var[k]);
                else variable_values[k].push_back((double)var_float[k]);
                //else variable_values[k].push_back(GetBinIndex((double)var_float[k], Bsig_M));
            }

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
                FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "DDBAR") {
                FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "SSBAR") {
                FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "CHARM") {
                FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else {
                printf("undefined job id!\n");
                exit(1);
            }

            // Multiplicity correction factor, it is not applied now. it is for a systematic uncertainty
            if (Ngamma_v200_index != variable_names.size()) Ngamma_v200 = var[Ngamma_v200_index];
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            // norm factor for each MXs
            double normfactor_MXs = 1.0;
            if (NormalizeAtEachMXs == true) {
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) normfactor_MXs = data_Nevt_MXs1->GetBinContent(1) / MC_Nevt_MXs1->GetBinContent(1);
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) normfactor_MXs = data_Nevt_MXs2->GetBinContent(1) / MC_Nevt_MXs2->GetBinContent(1);
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) normfactor_MXs = data_Nevt_MXs3->GetBinContent(1) / MC_Nevt_MXs3->GetBinContent(1);
            }

            weights->push_back(FEI_calibration_factor * CAL * weight_ri * Correction_multiplicity * additional_weight * normfactor_MXs);

            MC_one_bin->Fill(1.0, FEI_calibration_factor * CAL * weight_ri * Correction_multiplicity * additional_weight * normfactor_MXs);

        }
        input_file->Close();

    }

    free(var);
    free(var_float);

}

void LetsFillMC_correction(const char* dirname, std::vector<std::string> variable_names, std::vector<std::string> branch_names, std::vector<double> variable_values[], std::vector<int>* numberings, std::vector<double>* weights, std::string SampleName, double NormFactor = 1.0, int option = 0, double additional_weight = 1.0) {
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

    double* var = (double*)malloc(sizeof(double) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var[i] = 0.0;
    float* var_float = (float*)malloc(sizeof(float) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var_float[i] = 0.0;
    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    double FEI_calibration_factor = -1;

    double MC_MXs = -1;
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

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

    float BDTc = -1;
    double BDTc_correction = -1;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");
        TTree* tree_Xs;
        if (SampleName == "SIGNAL") tree_Xs = (TTree*)input_file->Get("Xs");

        for (int k = 0; k < (int)variable_names.size(); k++) {
            if (variable_names.at(k).find("bin index") != std::string::npos) continue;

            if (branch_names.at(k) == std::string("Upsilon")) {
                if (variable_names.at(k).find("MVA") == std::string::npos) tree_upsilon->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
                else tree_upsilon->SetBranchAddress(variable_names.at(k).c_str(), &var_float[k]);
            }
            else if (branch_names.at(k) == std::string("Bsig")) tree_Bsig->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else if (branch_names.at(k) == std::string("Btag")) tree_Btag->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else {
                printf("ERROR! \n");
                exit(1);
            }
        }

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

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
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

        Ngamma_v200_index = std::find(variable_names.begin(), variable_names.end(), std::string("extraInfo__boNgammav200__bc")) - variable_names.begin();
        if (Ngamma_v200_index == variable_names.size()) tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

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

        tree_upsilon->SetBranchAddress("MVA_Continuum", &BDTc);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            if (option == 1 && Upsilon_ID != 0) continue;
            else if (option == 2 && Upsilon_ID != 1) continue;
            else if (option == 0) {}
            else if ((option != 0) && (option != 1) && (option != 2)) {
                printf("improper option value!\n");
                exit(1);
            }

            // BDTc correction factor
            BDTc_correction = BDTcToWeight(BDTc) * NormFactor;

            // obtain true MXs
            if (SampleName == "SIGNAL") {

                // set MXs from filename
                if (std::strstr(names.at(i).c_str(), "B2Knunu") != nullptr) MC_MXs = Mxs_Bc_MC;
                else if (std::strstr(names.at(i).c_str(), "B2Kstarnunu") != nullptr) MC_MXs = Mxs_Bc_MC;
                else if (std::strstr(names.at(i).c_str(), "B2Xsnunu") != nullptr) MC_MXs = Mxs_Bc_MC;
                else if (std::strstr(names.at(i).c_str(), "B02K0nunu") != nullptr) MC_MXs = Mxs_B0_MC;
                else if (std::strstr(names.at(i).c_str(), "B02Kstar0nunu") != nullptr) MC_MXs = Mxs_B0_MC;
                else if (std::strstr(names.at(i).c_str(), "B02Xsnunu") != nullptr) MC_MXs = Mxs_B0_MC;
                else {
                    printf("[ERROR] unexpected filename for SIGNAL\n");
                    exit(1);
                }

                // sanity check
                if ((MC_MXs > 0.0) && (MC_MXs < 6.0)) {}
                else {
                    if (std::strstr(names.at(i).c_str(), "B2Knunu") != nullptr) MC_MXs = 0.4868;
                    else if (std::strstr(names.at(i).c_str(), "B2Kstarnunu") != nullptr) MC_MXs = 0.8916;
                    else if (std::strstr(names.at(i).c_str(), "B2Xsnunu") != nullptr) MC_MXs = 1.5;
                    else if (std::strstr(names.at(i).c_str(), "B02K0nunu") != nullptr) MC_MXs = 0.4868;
                    else if (std::strstr(names.at(i).c_str(), "B02Kstar0nunu") != nullptr) MC_MXs = 0.8916;
                    else if (std::strstr(names.at(i).c_str(), "B02Xsnunu") != nullptr) MC_MXs = 1.5;
                }

            }
            else {
                MC_MXs = -1;
                Mxs_Bc_MC = -1;
                Mxs_B0_MC = -1;
            }

            for (int k = 0; k < (int)variable_names.size(); k++) {
                if (variable_names.at(k).find("bin index") != std::string::npos) {
                    auto it = std::find(variable_names.begin(), variable_names.end(), "MVA_BB");
                    int index_FBDT_raw = std::distance(variable_names.begin(), it);
                    variable_values[k].push_back(ReturnBinIndex(var_float[index_FBDT_raw], Bsig_M));
                }
                else if (variable_names.at(k).find("MVA") == std::string::npos) variable_values[k].push_back(var[k]);
                else variable_values[k].push_back((double)var_float[k]);
                //else variable_values[k].push_back(GetBinIndex((double)var_float[k], Bsig_M));
            }

            // Fill numberings
            double weight_ri = 0.0;
            if (SampleName == "CHG") {
                numberings->push_back(static_cast<int>(MCsample::CHG));
                FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "MIX") {
                numberings->push_back(static_cast<int>(MCsample::MIX));
                FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "UUBAR") {
                numberings->push_back(static_cast<int>(MCsample::UUBAR));
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) FEI_calibration_factor = CAL_qq_K;
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) FEI_calibration_factor = CAL_qq_Kstar;
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) FEI_calibration_factor = CAL_qq_Xs;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "DDBAR") {
                numberings->push_back(static_cast<int>(MCsample::DDBAR));
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) FEI_calibration_factor = CAL_qq_K;
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) FEI_calibration_factor = CAL_qq_Kstar;
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) FEI_calibration_factor = CAL_qq_Xs;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "SSBAR") {
                numberings->push_back(static_cast<int>(MCsample::SSBAR));
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) FEI_calibration_factor = CAL_qq_K;
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) FEI_calibration_factor = CAL_qq_Kstar;
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) FEI_calibration_factor = CAL_qq_Xs;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "CHARM") {
                numberings->push_back(static_cast<int>(MCsample::CHARM));
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) FEI_calibration_factor = CAL_qq_K;
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) FEI_calibration_factor = CAL_qq_Kstar;
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) FEI_calibration_factor = CAL_qq_Xs;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "SIGNAL") {
                numberings->push_back(static_cast<int>(MCsample::SIGNAL));
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
            if (Ngamma_v200_index != variable_names.size()) Ngamma_v200 = var[Ngamma_v200_index];
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

            // norm factor for each MXs
            double normfactor_MXs = 1.0;
            if (NormalizeAtEachMXs == true) {
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) normfactor_MXs = data_Nevt_MXs1->GetBinContent(1) / MC_Nevt_MXs1->GetBinContent(1);
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) normfactor_MXs = data_Nevt_MXs2->GetBinContent(1) / MC_Nevt_MXs2->GetBinContent(1);
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) normfactor_MXs = data_Nevt_MXs3->GetBinContent(1) / MC_Nevt_MXs3->GetBinContent(1);
            }

            // scale factor to adjust into fit result
            double fitfactor = 1.0;
            if (DrawFitPlots == true) {
                auto it = std::find(variable_names.begin(), variable_names.end(), "MVA_BB");
                int index_FBDT_raw = std::distance(variable_names.begin(), it);
                double bin_index = ReturnBinIndex(var_float[index_FBDT_raw], Bsig_M);
                fitfactor = GetFitWeight(bin_index, MC_MXs, SampleName);
            }

            weights->push_back(FEI_calibration_factor* CAL* weight_ri* Correction_pi0* Correction_KID* Correction_PID* Correction_fake* Correction_Knn* Correction_Xnn* Correction_multiplicity* Correction_KpKLKL* Correction_KSKLKL* Correction_KstarKLKL* Correction_XKLKL* Correction_BtoDtoXKL* BDTc_correction* additional_weight * normfactor_MXs * fitfactor);

            MC_one_bin->Fill(1.0, FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * BDTc_correction * additional_weight * normfactor_MXs * fitfactor);

        }
        input_file->Close();

    }

    free(var);
    free(var_float);

}

void LetsFillMC_ESide(const char* dirname, std::vector<std::string> variable_names, std::vector<std::string> branch_names, std::vector<double> variable_values[], std::vector<int>* numberings, std::vector<double>* weights, std::string SampleName, int option = 0, double additional_weight = 1.0) {
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

    double* var = (double*)malloc(sizeof(double) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var[i] = 0.0;
    float* var_float = (float*)malloc(sizeof(float) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var_float[i] = 0.0;
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

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");
        TTree* tree_Xs;
        if (SampleName == "SIGNAL") tree_Xs = (TTree*)input_file->Get("Xs");

        for (int k = 0; k < (int)variable_names.size(); k++) {
            if (variable_names.at(k).find("bin index") != std::string::npos) continue;

            if (branch_names.at(k) == std::string("Upsilon")) {
                if (variable_names.at(k).find("MVA") == std::string::npos) tree_upsilon->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
                else tree_upsilon->SetBranchAddress(variable_names.at(k).c_str(), &var_float[k]);
            }
            else if (branch_names.at(k) == std::string("Bsig")) tree_Bsig->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else if (branch_names.at(k) == std::string("Btag")) tree_Btag->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else {
                printf("ERROR! \n");
                exit(1);
            }
        }

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

        Ngamma_v200_index = std::find(variable_names.begin(), variable_names.end(), std::string("extraInfo__boNgammav200__bc")) - variable_names.begin();
        if (Ngamma_v200_index == variable_names.size()) tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

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

            if (option == 1 && Upsilon_ID != 0) continue;
            else if (option == 2 && Upsilon_ID != 1) continue;
            else if (option == 0) {}
            else if ((option != 0) && (option != 1) && (option != 2)) {
                printf("improper option value!\n");
                exit(1);
            }

            for (int k = 0; k < (int)variable_names.size(); k++) {
                if (variable_names.at(k).find("bin index") != std::string::npos) {
                    auto it = std::find(variable_names.begin(), variable_names.end(), "MVA_BB");
                    int index_FBDT_raw = std::distance(variable_names.begin(), it);
                    variable_values[k].push_back(ReturnBinIndex(var_float[index_FBDT_raw], Bsig_M));
                }
                else if (variable_names.at(k).find("MVA") == std::string::npos) variable_values[k].push_back(var[k]);
                else variable_values[k].push_back((double)var_float[k]);
                //else variable_values[k].push_back(GetBinIndex((double)var_float[k], Bsig_M));
            }

            // Fill numberings
            double weight_ri = 0.0;
            if (SampleName == "CHG") {
                numberings->push_back(static_cast<int>(MCsample::CHG));
                FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "MIX") {
                numberings->push_back(static_cast<int>(MCsample::MIX));
                FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "UUBAR") {
                numberings->push_back(static_cast<int>(MCsample::UUBAR));
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) FEI_calibration_factor = CAL_qq_K;
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) FEI_calibration_factor = CAL_qq_Kstar;
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) FEI_calibration_factor = CAL_qq_Xs;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "DDBAR") {
                numberings->push_back(static_cast<int>(MCsample::DDBAR));
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) FEI_calibration_factor = CAL_qq_K;
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) FEI_calibration_factor = CAL_qq_Kstar;
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) FEI_calibration_factor = CAL_qq_Xs;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "SSBAR") {
                numberings->push_back(static_cast<int>(MCsample::SSBAR));
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) FEI_calibration_factor = CAL_qq_K;
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) FEI_calibration_factor = CAL_qq_Kstar;
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) FEI_calibration_factor = CAL_qq_Xs;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "CHARM") {
                numberings->push_back(static_cast<int>(MCsample::CHARM));
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) FEI_calibration_factor = CAL_qq_K;
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) FEI_calibration_factor = CAL_qq_Kstar;
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) FEI_calibration_factor = CAL_qq_Xs;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "SIGNAL") {
                numberings->push_back(static_cast<int>(MCsample::SIGNAL));
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
            double Correction_LID = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_LID = Correction_LID * std::pow(corrector_LID.GetCorrectionFactortrueE(0, i_fake, MCTYPE), temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_LID = Correction_LID * std::pow(corrector_LID.GetCorrectionFactortrueE(1, i_fake, MCTYPE), temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(2, i_fake, MCTYPE), temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(3, i_fake, MCTYPE), temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(2, i_fake, MCTYPE), temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(3, i_fake, MCTYPE), temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorCancelOutObtainWeight(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, names.at(i), MCTYPE, true);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor, it is not applied now. it is for a systematic uncertainty
            if (Ngamma_v200_index != variable_names.size()) Ngamma_v200 = var[Ngamma_v200_index];
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

            // norm factor for each MXs
            double normfactor_MXs = 1.0;
            if (NormalizeAtEachMXs == true) {
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) normfactor_MXs = data_Nevt_MXs1->GetBinContent(1) / MC_Nevt_MXs1->GetBinContent(1);
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) normfactor_MXs = data_Nevt_MXs2->GetBinContent(1) / MC_Nevt_MXs2->GetBinContent(1);
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) normfactor_MXs = data_Nevt_MXs3->GetBinContent(1) / MC_Nevt_MXs3->GetBinContent(1);
            }

            weights->push_back(FEI_calibration_factor* CAL* weight_ri* Correction_pi0* Correction_KID* Correction_PID* Correction_fake* Correction_LID* Correction_Knn* Correction_Xnn* Correction_multiplicity* Correction_KpKLKL* Correction_KSKLKL* Correction_KstarKLKL* Correction_XKLKL* Correction_BtoDtoXKL* additional_weight * normfactor_MXs);

            MC_one_bin->Fill(1.0, FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_LID * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * additional_weight * normfactor_MXs);

        }
        input_file->Close();

    }

    free(var);
    free(var_float);

}

void LetsFillMC_MUSide(const char* dirname, std::vector<std::string> variable_names, std::vector<std::string> branch_names, std::vector<double> variable_values[], std::vector<int>* numberings, std::vector<double>* weights, std::string SampleName, int option = 0, double additional_weight = 1.0) {
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

    double* var = (double*)malloc(sizeof(double) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var[i] = 0.0;
    float* var_float = (float*)malloc(sizeof(float) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var_float[i] = 0.0;
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

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");
        TTree* tree_Xs;
        if (SampleName == "SIGNAL") tree_Xs = (TTree*)input_file->Get("Xs");

        for (int k = 0; k < (int)variable_names.size(); k++) {
            if (variable_names.at(k).find("bin index") != std::string::npos) continue;

            if (branch_names.at(k) == std::string("Upsilon")) {
                if (variable_names.at(k).find("MVA") == std::string::npos) tree_upsilon->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
                else tree_upsilon->SetBranchAddress(variable_names.at(k).c_str(), &var_float[k]);
            }
            else if (branch_names.at(k) == std::string("Bsig")) tree_Bsig->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else if (branch_names.at(k) == std::string("Btag")) tree_Btag->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else {
                printf("ERROR! \n");
                exit(1);
            }
        }

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

        Ngamma_v200_index = std::find(variable_names.begin(), variable_names.end(), std::string("extraInfo__boNgammav200__bc")) - variable_names.begin();
        if (Ngamma_v200_index == variable_names.size()) tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

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

            if (option == 1 && Upsilon_ID != 0) continue;
            else if (option == 2 && Upsilon_ID != 1) continue;
            else if (option == 0) {}
            else if ((option != 0) && (option != 1) && (option != 2)) {
                printf("improper option value!\n");
                exit(1);
            }

            for (int k = 0; k < (int)variable_names.size(); k++) {
                if (variable_names.at(k).find("bin index") != std::string::npos) {
                    auto it = std::find(variable_names.begin(), variable_names.end(), "MVA_BB");
                    int index_FBDT_raw = std::distance(variable_names.begin(), it);
                    variable_values[k].push_back(ReturnBinIndex(var_float[index_FBDT_raw], Bsig_M));
                }
                else if (variable_names.at(k).find("MVA") == std::string::npos) variable_values[k].push_back(var[k]);
                else variable_values[k].push_back((double)var_float[k]);
                //else variable_values[k].push_back(GetBinIndex((double)var_float[k], Bsig_M));
            }

            // Fill numberings
            double weight_ri = 0.0;
            if (SampleName == "CHG") {
                numberings->push_back(static_cast<int>(MCsample::CHG));
                FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "MIX") {
                numberings->push_back(static_cast<int>(MCsample::MIX));
                FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "UUBAR") {
                numberings->push_back(static_cast<int>(MCsample::UUBAR));
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) FEI_calibration_factor = CAL_qq_K;
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) FEI_calibration_factor = CAL_qq_Kstar;
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) FEI_calibration_factor = CAL_qq_Xs;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "DDBAR") {
                numberings->push_back(static_cast<int>(MCsample::DDBAR));
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) FEI_calibration_factor = CAL_qq_K;
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) FEI_calibration_factor = CAL_qq_Kstar;
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) FEI_calibration_factor = CAL_qq_Xs;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "SSBAR") {
                numberings->push_back(static_cast<int>(MCsample::SSBAR));
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) FEI_calibration_factor = CAL_qq_K;
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) FEI_calibration_factor = CAL_qq_Kstar;
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) FEI_calibration_factor = CAL_qq_Xs;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "CHARM") {
                numberings->push_back(static_cast<int>(MCsample::CHARM));
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) FEI_calibration_factor = CAL_qq_K;
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) FEI_calibration_factor = CAL_qq_Kstar;
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) FEI_calibration_factor = CAL_qq_Xs;
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "SIGNAL") {
                numberings->push_back(static_cast<int>(MCsample::SIGNAL));
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
            double Correction_LID = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(2, i_fake, MCTYPE), temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(3, i_fake, MCTYPE), temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_LID = Correction_LID * std::pow(corrector_LID.GetCorrectionFactortrueMU(0, i_fake, MCTYPE), temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_LID = Correction_LID * std::pow(corrector_LID.GetCorrectionFactortrueMU(1, i_fake, MCTYPE), temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(2, i_fake, MCTYPE), temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeMU(3, i_fake, MCTYPE), temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorCancelOutObtainWeight(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, names.at(i), MCTYPE, true);

            // Xsnn correction factor
            double Correction_Xnn = corrector_Xsnn.GetCorrectionFactorAtGeneric(invM_Xnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn, N_Xplusnn + N_Xzeronn);

            // Multiplicity correction factor, it is not applied now. it is for a systematic uncertainty
            if (Ngamma_v200_index != variable_names.size()) Ngamma_v200 = var[Ngamma_v200_index];
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

            // norm factor for each MXs
            double normfactor_MXs = 1.0;
            if (NormalizeAtEachMXs == true) {
                if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) normfactor_MXs = data_Nevt_MXs1->GetBinContent(1) / MC_Nevt_MXs1->GetBinContent(1);
                else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) normfactor_MXs = data_Nevt_MXs2->GetBinContent(1) / MC_Nevt_MXs2->GetBinContent(1);
                else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) normfactor_MXs = data_Nevt_MXs3->GetBinContent(1) / MC_Nevt_MXs3->GetBinContent(1);
            }

            weights->push_back(FEI_calibration_factor* CAL* weight_ri* Correction_pi0* Correction_KID* Correction_PID* Correction_fake* Correction_LID* Correction_Knn* Correction_Xnn* Correction_multiplicity* Correction_KpKLKL* Correction_KSKLKL* Correction_KstarKLKL* Correction_XKLKL* Correction_BtoDtoXKL* additional_weight* normfactor_MXs);

            MC_one_bin->Fill(1.0, FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_LID * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * additional_weight * normfactor_MXs);

        }
        input_file->Close();

    }

    free(var);
    free(var_float);

}

void LetsCountMC(const char* dirname, std::string SampleName, int option = 0, double additional_weight = 1.0) {
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

    if (NormalizeAtEachMXs == false) {
        printf("`NormalizeAtEachMXs` is false. Just skip to count Nevt\n");
    }

    double* var = (double*)malloc(sizeof(double) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var[i] = 0.0;
    float* var_float = (float*)malloc(sizeof(float) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var_float[i] = 0.0;
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

    std::vector<string> names;
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

            if (option == 1 && Upsilon_ID != 0) continue;
            else if (option == 2 && Upsilon_ID != 1) continue;
            else if (option == 0) {}
            else if ((option != 0) && (option != 1) && (option != 2)) {
                printf("improper option value!\n");
                exit(1);
            }

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

            if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) MC_Nevt_MXs1->Fill(1.0, FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * additional_weight);
            else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) MC_Nevt_MXs2->Fill(1.0, FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * additional_weight);
            else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) MC_Nevt_MXs3->Fill(1.0, FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * additional_weight);

        }
        input_file->Close();

    }

    free(var);
    free(var_float);

}

void LetsCountdata(const char* dirname, const char* included_string = "root", int option = 0) {
    /*
    option 0: select all Btag
    option 1: select Btag+
    option 2: select Btag0
    */

    if (NormalizeAtEachMXs == false) {
        printf("`NormalizeAtEachMXs` is false. Just skip to count Nevt\n");
    }

    double* var = (double*)malloc(sizeof(double) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var[i] = 0.0;
    float* var_float = (float*)malloc(sizeof(float) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var_float[i] = 0.0;
    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;

    double Bsig_M = -1;

    std::vector<string> names;
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

            if (option == 1 && Upsilon_ID != 0) continue;
            else if (option == 2 && Upsilon_ID != 1) continue;
            else if (option == 0) {}
            else if ((option != 0) && (option != 1) && (option != 2)) {
                printf("improper option value!\n");
                exit(1);
            }

            if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) data_Nevt_MXs1->Fill(1.0);
            else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) data_Nevt_MXs2->Fill(1.0);
            else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) data_Nevt_MXs3->Fill(1.0);

        }
        input_file->Close();

    }

    free(var);
    free(var_float);

}

void LetsCountembeddedMC(const char* dirname, std::string SampleName, int option = 0, double additional_weight = 1.0) {
    /*
    option 0: select all Btag
    option 1: select Btag+
    option 2: select Btag0
    */

    if (NormalizeAtEachMXs == false) {
        printf("`NormalizeAtEachMXs` is false. Just skip to count Nevt\n");
    }

    double* var = (double*)malloc(sizeof(double) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var[i] = 0.0;
    float* var_float = (float*)malloc(sizeof(float) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var_float[i] = 0.0;
    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;

    double FEI_calibration_factor = -1;

    int Ngamma_v200_index = -1;
    double Ngamma_v200 = -1;

    double Bsig_M = -1;

    std::vector<string> names;
    load_files(dirname, &names);

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

        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &Ngamma_v200);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            if (option == 1 && Upsilon_ID != 0) continue;
            else if (option == 2 && Upsilon_ID != 1) continue;
            else if (option == 0) {}
            else if ((option != 0) && (option != 1) && (option != 2)) {
                printf("improper option value!\n");
                exit(1);
            }

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
                FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "DDBAR") {
                FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "SSBAR") {
                FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else if (SampleName == "CHARM") {
                FEI_calibration_factor = corrector_FEI.GetFEICalFactor(Upsilon_ID, Btag_ID, MCTYPE);
                weight_ri = ObtainWeight(SampleName.c_str(), MCTYPE, "validation", names.at(i));
            }
            else {
                printf("undefined job id!\n");
                exit(1);
            }

            // Multiplicity correction factor, it is not applied now. it is for a systematic uncertainty
            double Correction_multiplicity = corrector_Multiplicity.GetCorrectionFactor(Ngamma_v200);

            if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) MC_Nevt_MXs1->Fill(1.0, FEI_calibration_factor * CAL * weight_ri * Correction_multiplicity * additional_weight);
            else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) MC_Nevt_MXs2->Fill(1.0, FEI_calibration_factor * CAL * weight_ri * Correction_multiplicity * additional_weight);
            else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) MC_Nevt_MXs3->Fill(1.0, FEI_calibration_factor * CAL * weight_ri * Correction_multiplicity * additional_weight);

        }
        input_file->Close();

    }

    free(var);
    free(var_float);

}

void LetsCountMC_correction(const char* dirname, std::string SampleName, double NormFactor = 1.0, int option = 0, double additional_weight = 1.0) {
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

    if (NormalizeAtEachMXs == false) {
        printf("`NormalizeAtEachMXs` is false. Just skip to count Nevt\n");
    }

    double* var = (double*)malloc(sizeof(double) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var[i] = 0.0;
    float* var_float = (float*)malloc(sizeof(float) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var_float[i] = 0.0;
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

    float BDTc = -1;
    double BDTc_correction = -1;

    std::vector<string> names;
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

        tree_upsilon->SetBranchAddress("MVA_Continuum", &BDTc);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            if (option == 1 && Upsilon_ID != 0) continue;
            else if (option == 2 && Upsilon_ID != 1) continue;
            else if (option == 0) {}
            else if ((option != 0) && (option != 1) && (option != 2)) {
                printf("improper option value!\n");
                exit(1);
            }

            // BDTc correction factor
            BDTc_correction = BDTcToWeight(BDTc) * NormFactor;

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

            if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) MC_Nevt_MXs1->Fill(1.0, FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * BDTc_correction * additional_weight);
            else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) MC_Nevt_MXs2->Fill(1.0, FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * BDTc_correction * additional_weight);
            else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) MC_Nevt_MXs3->Fill(1.0, FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * BDTc_correction * additional_weight);

        }
        input_file->Close();

    }

    free(var);
    free(var_float);

}

void LetsCountMC_ESide(const char* dirname, std::string SampleName, int option = 0, double additional_weight = 1.0) {
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

    if (NormalizeAtEachMXs == false) {
        printf("`NormalizeAtEachMXs` is false. Just skip to count Nevt\n");
    }

    double* var = (double*)malloc(sizeof(double) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var[i] = 0.0;
    float* var_float = (float*)malloc(sizeof(float) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var_float[i] = 0.0;
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

    std::vector<string> names;
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

            if (option == 1 && Upsilon_ID != 0) continue;
            else if (option == 2 && Upsilon_ID != 1) continue;
            else if (option == 0) {}
            else if ((option != 0) && (option != 1) && (option != 2)) {
                printf("improper option value!\n");
                exit(1);
            }

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
            double Correction_LID = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_LID = Correction_LID * std::pow(corrector_LID.GetCorrectionFactortrueE(0, i_fake, MCTYPE), temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_LID = Correction_LID * std::pow(corrector_LID.GetCorrectionFactortrueE(1, i_fake, MCTYPE), temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(2, i_fake, MCTYPE), temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(3, i_fake, MCTYPE), temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
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

            if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) MC_Nevt_MXs1->Fill(1.0, FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_LID * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * additional_weight);
            else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) MC_Nevt_MXs2->Fill(1.0, FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_LID * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * additional_weight);
            else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) MC_Nevt_MXs3->Fill(1.0, FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_LID * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * additional_weight);

        }
        input_file->Close();

    }

    free(var);
    free(var_float);

}

void LetsCountMC_MUSide(const char* dirname, std::string SampleName, int option = 0, double additional_weight = 1.0) {
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

    if (NormalizeAtEachMXs == false) {
        printf("`NormalizeAtEachMXs` is false. Just skip to count Nevt\n");
    }

    double* var = (double*)malloc(sizeof(double) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var[i] = 0.0;
    float* var_float = (float*)malloc(sizeof(float) * Nvar_num); for (int i = 0; i < Nvar_num; i++) var_float[i] = 0.0;
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

    std::vector<string> names;
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

            if (option == 1 && Upsilon_ID != 0) continue;
            else if (option == 2 && Upsilon_ID != 1) continue;
            else if (option == 0) {}
            else if ((option != 0) && (option != 1) && (option != 2)) {
                printf("improper option value!\n");
                exit(1);
            }

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
            double Correction_LID = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(2, i_PID, MCTYPE), temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(corrector_PID.GetCorrectionFactor(3, i_PID, MCTYPE), temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(corrector_pi0.GetCorrectionFactor(i_pi0, MCTYPE), temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(2, i_fake, MCTYPE), temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(corrector_FakePID.GetCorrectionFactorfakeE(3, i_fake, MCTYPE), temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_LID = Correction_LID * std::pow(corrector_LID.GetCorrectionFactortrueMU(0, i_fake, MCTYPE), temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_LID = Correction_LID * std::pow(corrector_LID.GetCorrectionFactortrueMU(1, i_fake, MCTYPE), temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
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

            if ((Bsig_M > 0.0) && (Bsig_M < 0.6)) MC_Nevt_MXs1->Fill(1.0, FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_LID * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * additional_weight);
            else if ((Bsig_M > 0.6) && (Bsig_M < 1.0)) MC_Nevt_MXs2->Fill(1.0, FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_LID * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * additional_weight);
            else if ((Bsig_M > 1.0) && (Bsig_M < 2.0)) MC_Nevt_MXs3->Fill(1.0, FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_LID * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * additional_weight);

        }
        input_file->Close();

    }

    free(var);
    free(var_float);

}


typedef struct _Nevt {
    double NevtwithoutCorrection;
    double NevtwithCorrection;
} Nevt;

void NevtCount_ri(const char* dirname, std::string SampleName, Nevt* nevt, double additional_weight = 1.0) {
    /*
    SampleName for Knn
    CHG
    MIX
    UUBAR
    DDBAR
    SSBAR
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

    double MC_MXs = -1;
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

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

    int nBp = -1;
    int nB0 = -1;

    double Bsig_M = -1;

    float MVA_BB = -1;
    float BDTc = -1;
    double BDTc_correction = -1;

    std::vector<string> names;
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

            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
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

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_BB);
        tree_upsilon->SetBranchAddress("MVA_Continuum", &BDTc);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            // BDTc correction factor
            BDTc_correction = BDTcToWeight(BDTc);

            // obtain true MXs
            if (SampleName == "SIGNAL") {

                // set MXs from filename
                if (std::strstr(names.at(i).c_str(), "B2Knunu") != nullptr) MC_MXs = Mxs_Bc_MC;
                else if (std::strstr(names.at(i).c_str(), "B2Kstarnunu") != nullptr) MC_MXs = Mxs_Bc_MC;
                else if (std::strstr(names.at(i).c_str(), "B2Xsnunu") != nullptr) MC_MXs = Mxs_Bc_MC;
                else if (std::strstr(names.at(i).c_str(), "B02K0nunu") != nullptr) MC_MXs = Mxs_B0_MC;
                else if (std::strstr(names.at(i).c_str(), "B02Kstar0nunu") != nullptr) MC_MXs = Mxs_B0_MC;
                else if (std::strstr(names.at(i).c_str(), "B02Xsnunu") != nullptr) MC_MXs = Mxs_B0_MC;
                else {
                    printf("[ERROR] unexpected filename for SIGNAL\n");
                    exit(1);
                }

                // sanity check
                if ((MC_MXs > 0.0) && (MC_MXs < 6.0)) {}
                else {
                    if (std::strstr(names.at(i).c_str(), "B2Knunu") != nullptr) MC_MXs = 0.4868;
                    else if (std::strstr(names.at(i).c_str(), "B2Kstarnunu") != nullptr) MC_MXs = 0.8916;
                    else if (std::strstr(names.at(i).c_str(), "B2Xsnunu") != nullptr) MC_MXs = 1.5;
                    else if (std::strstr(names.at(i).c_str(), "B02K0nunu") != nullptr) MC_MXs = 0.4868;
                    else if (std::strstr(names.at(i).c_str(), "B02Kstar0nunu") != nullptr) MC_MXs = 0.8916;
                    else if (std::strstr(names.at(i).c_str(), "B02Xsnunu") != nullptr) MC_MXs = 1.5;
                }

            }
            else {
                MC_MXs = -1;
                Mxs_Bc_MC = -1;
                Mxs_B0_MC = -1;
            }

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

            // scale factor to adjust into fit result
            double fitfactor = 1.0;
            if (DrawFitPlots == true) {
                double bin_index = ReturnBinIndex(MVA_BB, Bsig_M);
                fitfactor = GetFitWeight(bin_index, MC_MXs, SampleName);
            }

            nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * additional_weight * fitfactor;
            nevt->NevtwithCorrection = nevt->NevtwithCorrection + FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_Xnn * Correction_multiplicity * Correction_KpKLKL * Correction_KSKLKL * Correction_KstarKLKL * Correction_XKLKL * Correction_BtoDtoXKL * BDTc_correction * additional_weight * fitfactor;


        }
        input_file->Close();

    }

}

void PrintDataMCRatio(THStack* MC_stack, TH1D* data_hist, TH1D* ratio_hist, const char* fname = "dataMCratio.txt") {

    FILE* fp = fopen(fname, "w");
    fprintf(fp, "%d\n", RarityBins);
    for (int i = 0; i < RarityBins; i++) { 
        double relative_uncertainty = 1.0;
        if (ratio_hist->GetBinContent(i + 1) > MyEPSILON) relative_uncertainty = ratio_hist->GetBinError(i + 1) / ratio_hist->GetBinContent(i + 1);
        else relative_uncertainty = 1.0;
        fprintf(fp, "%lf\n", relative_uncertainty);
    }
    fclose(fp);
}

void Printchi2(THStack* MC_stack, TH1D* data_hist) {

    TH1* temp = (TH1*) MC_stack->GetStack()->Last();

    if (temp->GetNbinsX() != data_hist->GetNbinsX()) {
        printf("[Printchi2] Nbin is different\n");
        exit(1);
    }

    double chi2 = 0.0;
    for (int i = 0; i < temp->GetNbinsX(); i++) {
        double MC_value = temp->GetBinContent(i + 1);
        double MC_error = temp->GetBinError(i + 1);
        double data_value = data_hist->GetBinContent(i + 1);
        double data_error = data_hist->GetBinError(i + 1);

        chi2 = chi2 + ((MC_value - data_value) * (MC_value - data_value) / (data_error * data_error));
    }

    printf("chi2: %lf\n", chi2);

    delete temp;

}

void PrintDataMC(THStack* MC_stack, TH1D* data_hist) {

    TH1* temp = (TH1*)MC_stack->GetStack()->Last();

    if (temp->GetNbinsX() != data_hist->GetNbinsX()) {
        printf("[PrintDataMC] Nbin is different\n");
        exit(1);
    }

    printf("===== MC =====\n");
    for (int i = 0; i < temp->GetNbinsX(); i++) {
        double MC_value = temp->GetBinContent(i + 1);
        double MC_error = temp->GetBinError(i + 1);
        printf("%d th bin: %lf +- %lf\n", i, MC_value, MC_error);
    }
    printf("===== data =====\n");
    for (int i = 0; i < data_hist->GetNbinsX(); i++) {
        double data_value = data_hist->GetBinContent(i + 1);
        double data_error = data_hist->GetBinError(i + 1);
        printf("%d th bin: %lf +- %lf\n", i, data_value, data_error);
    }

    delete temp;

}