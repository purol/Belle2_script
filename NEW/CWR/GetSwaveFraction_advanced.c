// last update: 2021-10-13
// for Belle2 data

#include <vector>
#include "template.h"
#include "constants.h"
#include "ObtainWeight.h"
#include "correctors.h"
#include "base.h"

# define Nentry 2000
# define MIN_Q2 0.0
# define MAX_Q2 20.0
# define Nfine 1000000

# define Scale_Kplusstar_here (N_Kplusstar_nunubar_LS1/2000000.0)
# define Scale_K0star_here (N_K0star_nunubar_LS1/2000000.0)

# define MCTYPE "MC15ri"

# define KnunuOnly false

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

/*
* [q2 bin][mass bin]
* q2 bin: [-inf,0.1], [0.1,0.98], [0.98,1.1], [1.1,2.5], [2.5,4.0], [4.0,6.0], [6.0,8.0], [8.0,11.0], [11.0,12.5], [12.5,15.0], [15.0,17.0], [17.0,19.0], [19.0, inf]
* mass bin: [-inf,0.644], [0.644,0.796], [0.796,0.996], [0.996,1.2], [1.2,inf]
*/
double N_Kplusstar[13][5] = { 0.0 }; 
double N_K0star[13][5] = { 0.0 };
double N_kappaplus[13][5] = { 0.0 };
double N_kappa0[13][5] = { 0.0 };
double scale[13][5] = { 0.0 };

typedef struct {
    double q_square;
    double MKstar;
} Data;

void load_files_here(const char *dirname, std::vector<string>* names){
   TSystemDirectory dir(dirname, dirname);
   TList *files = dir.GetListOfFiles();
   if (files) {
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next())) {
         fname = file->GetName();
         if (!file->IsDirectory() && fname.EndsWith(".root")) {
            names->push_back(fname.Data());
         }
      }
   }
}

void load_files_here(const char* dirname, std::vector<std::string>* names, const char* included_string) {
    TSystemDirectory dir(dirname, dirname);
    TList* files = dir.GetListOfFiles();
    if (files) {
        TSystemFile* file;
        TString fname;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(".root") && fname.Contains(included_string)) {
                names->push_back(fname.Data());
            }
        }
    }
}

double** GetPDFs(const char* dirname, const char* included_string, const char* type, const char* sample, bool Get06GeVhigher1GeVunder_, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get nominal PDF with appropriate correction
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
    double temp_N_pi0_syst_MC15ri[N_pi0_syst_MC15ri] = { 0.0 };
    double temp_N_pi0_syst_MC15rd[N_pi0_syst_MC15rd] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    double KS0_flight_distance = 0;
    double KS0_costheta = 0;
    double KS0_p = 0;
    double KS0_D1p = 0;
    double KS0_D2p = 0;

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

    double total_Nevt = 0;
    double** Nevt = nullptr;
    Nevt = (double**)malloc(sizeof(double*) * 13);
    for (int i = 0; i < 13; i++) Nevt[i] = (double*)malloc(sizeof(double) * 5);
    for (int i = 0; i < 13; i++) for (int j = 0; j < 5; j++) Nevt[i][j] = 0.0;
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
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_D1p", &KS0_D1p);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_D2p", &KS0_D2p);
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
            double MC_MXs = -1;
            if (strcmp(sample, "SIGNAL") == 0) {
                if (strcmp(type, "Bplus") == 0) MC_MXs = Mxs_Bc_MC;
                else if (strcmp(type, "Bzero") == 0) MC_MXs = Mxs_B0_MC;

                // select 0.6 < MXs <1.0 GeV only
                if ((0.6 < MC_MXs) && (MC_MXs < 1.0)) {}
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
            else if (std::string(MCTYPE) == "MC15rd") Correction_KS0 = Correction_KS0 * corrector_KS0.GetCorrectionFactor(KS0_p, KS0_costheta, KS0_flight_distance, KS0_D1p, KS0_D2p, MCTYPE);

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
            else if (CorrectionType == "kappa") {
                int index_qsq = -1;
                int index_mass = -1;

                if ((invM * invM < 0.1)) index_qsq = 0;
                else if ((0.1 < invM * invM) && (invM * invM < 0.98)) index_qsq = 1;
                else if ((0.98 < invM * invM) && (invM * invM < 1.1)) index_qsq = 2;
                else if ((1.1 < invM * invM) && (invM * invM < 2.5)) index_qsq = 3;
                else if ((2.5 < invM * invM) && (invM * invM < 4.0)) index_qsq = 4;
                else if ((4.0 < invM * invM) && (invM * invM < 6.0)) index_qsq = 5;
                else if ((6.0 < invM * invM) && (invM * invM < 8.0)) index_qsq = 6;
                else if ((8.0 < invM * invM) && (invM * invM < 11.0)) index_qsq = 7;
                else if ((11.0 < invM * invM) && (invM * invM < 12.5)) index_qsq = 8;
                else if ((12.5 < invM * invM) && (invM * invM < 15.0)) index_qsq = 9;
                else if ((15.0 < invM * invM) && (invM * invM < 17.0)) index_qsq = 10;
                else if ((17.0 < invM * invM) && (invM * invM < 19.0)) index_qsq = 11;
                else if (19.0 < invM * invM) index_qsq = 12;

                if ((MC_MXs < 0.644)) index_mass = 0;
                else if ((0.644 < MC_MXs) && (MC_MXs < 0.796)) index_mass = 1;
                else if ((0.796 < MC_MXs) && (MC_MXs < 0.996)) index_mass = 2;
                else if ((0.996 < MC_MXs) && (MC_MXs < 1.2)) index_mass = 3;
                else if ((1.2 < MC_MXs)) index_mass = 4;

                if(strcmp(type, "Bzero") == 0) total_weight = total_weight * scale[index_qsq][index_mass];
                else if (strcmp(type, "Bplus") == 0) total_weight = total_weight * scale[index_qsq][index_mass] * (TBp / TB0);
            }

            int index_qsq = -1;
            int index_mass = -1;

            if ((invM * invM < 0.1)) index_qsq = 0;
            else if ((0.1 < invM * invM) && (invM * invM < 0.98)) index_qsq = 1;
            else if ((0.98 < invM * invM) && (invM * invM < 1.1)) index_qsq = 2;
            else if ((1.1 < invM * invM) && (invM * invM < 2.5)) index_qsq = 3;
            else if ((2.5 < invM * invM) && (invM * invM < 4.0)) index_qsq = 4;
            else if ((4.0 < invM * invM) && (invM * invM < 6.0)) index_qsq = 5;
            else if ((6.0 < invM * invM) && (invM * invM < 8.0)) index_qsq = 6;
            else if ((8.0 < invM * invM) && (invM * invM < 11.0)) index_qsq = 7;
            else if ((11.0 < invM * invM) && (invM * invM < 12.5)) index_qsq = 8;
            else if ((12.5 < invM * invM) && (invM * invM < 15.0)) index_qsq = 9;
            else if ((15.0 < invM * invM) && (invM * invM < 17.0)) index_qsq = 10;
            else if ((17.0 < invM * invM) && (invM * invM < 19.0)) index_qsq = 11;
            else if (19.0 < invM * invM) index_qsq = 12;

            if ((MC_MXs < 0.644)) index_mass = 0;
            else if ((0.644 < MC_MXs) && (MC_MXs < 0.796)) index_mass = 1;
            else if ((0.796 < MC_MXs) && (MC_MXs < 0.996)) index_mass = 2;
            else if ((0.996 < MC_MXs) && (MC_MXs < 1.2)) index_mass = 3;
            else if ((1.2 < MC_MXs)) index_mass = 4;

            if (Get06GeVhigher1GeVunder_) {
                if ((0.6 < MC_MXs) && (MC_MXs < 1.0)) {}
                else continue;
            }

            total_Nevt = total_Nevt + total_weight * projection_multiplication;
            if((index_qsq != -1) && (index_mass != -1)) Nevt[index_qsq][index_mass] = Nevt[index_qsq][index_mass] + total_weight * projection_multiplication;
        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, total_Nevt);

    return Nevt;
}

double GetPDFs(const char* dirname, const char* included_string, TH1D* hist, const char* type, const char* sample, bool Get06GeVhigher1GeVunder_, double weight_var = 1.0, std::string CorrectionType = "otherwise") { // get nominal PDF with appropriate correction
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
    double temp_N_pi0_syst_MC15ri[N_pi0_syst_MC15ri] = { 0.0 };
    double temp_N_pi0_syst_MC15rd[N_pi0_syst_MC15rd] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    double KS0_flight_distance = 0;
    double KS0_costheta = 0;
    double KS0_p = 0;
    double KS0_D1p = 0;
    double KS0_D2p = 0;

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
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_D1p", &KS0_D1p);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_D2p", &KS0_D2p);
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
            double MC_MXs = -1;
            if (strcmp(sample, "SIGNAL") == 0) {
                if (strcmp(type, "Bplus") == 0) MC_MXs = Mxs_Bc_MC;
                else if (strcmp(type, "Bzero") == 0) MC_MXs = Mxs_B0_MC;

                // select 0.6 < MXs <1.0 GeV only
                if ((0.6 < MC_MXs) && (MC_MXs < 1.0)) {}
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
            else if (std::string(MCTYPE) == "MC15rd") Correction_KS0 = Correction_KS0 * corrector_KS0.GetCorrectionFactor(KS0_p, KS0_costheta, KS0_flight_distance, KS0_D1p, KS0_D2p, MCTYPE);

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
            else if (CorrectionType == "kappa") {
                int index_qsq = -1;
                int index_mass = -1;

                if ((invM * invM < 0.1)) index_qsq = 0;
                else if ((0.1 < invM * invM) && (invM * invM < 0.98)) index_qsq = 1;
                else if ((0.98 < invM * invM) && (invM * invM < 1.1)) index_qsq = 2;
                else if ((1.1 < invM * invM) && (invM * invM < 2.5)) index_qsq = 3;
                else if ((2.5 < invM * invM) && (invM * invM < 4.0)) index_qsq = 4;
                else if ((4.0 < invM * invM) && (invM * invM < 6.0)) index_qsq = 5;
                else if ((6.0 < invM * invM) && (invM * invM < 8.0)) index_qsq = 6;
                else if ((8.0 < invM * invM) && (invM * invM < 11.0)) index_qsq = 7;
                else if ((11.0 < invM * invM) && (invM * invM < 12.5)) index_qsq = 8;
                else if ((12.5 < invM * invM) && (invM * invM < 15.0)) index_qsq = 9;
                else if ((15.0 < invM * invM) && (invM * invM < 17.0)) index_qsq = 10;
                else if ((17.0 < invM * invM) && (invM * invM < 19.0)) index_qsq = 11;
                else if (19.0 < invM * invM) index_qsq = 12;

                if ((MC_MXs < 0.644)) index_mass = 0;
                else if ((0.644 < MC_MXs) && (MC_MXs < 0.796)) index_mass = 1;
                else if ((0.796 < MC_MXs) && (MC_MXs < 0.996)) index_mass = 2;
                else if ((0.996 < MC_MXs) && (MC_MXs < 1.2)) index_mass = 3;
                else if ((1.2 < MC_MXs)) index_mass = 4;

                if (strcmp(type, "Bzero") == 0) total_weight = total_weight * scale[index_qsq][index_mass];
                else if (strcmp(type, "Bplus") == 0) total_weight = total_weight * scale[index_qsq][index_mass] * (TBp / TB0);
            }

            int index_qsq = -1;
            int index_mass = -1;

            if ((invM * invM < 0.1)) index_qsq = 0;
            else if ((0.1 < invM * invM) && (invM * invM < 0.98)) index_qsq = 1;
            else if ((0.98 < invM * invM) && (invM * invM < 1.1)) index_qsq = 2;
            else if ((1.1 < invM * invM) && (invM * invM < 2.5)) index_qsq = 3;
            else if ((2.5 < invM * invM) && (invM * invM < 4.0)) index_qsq = 4;
            else if ((4.0 < invM * invM) && (invM * invM < 6.0)) index_qsq = 5;
            else if ((6.0 < invM * invM) && (invM * invM < 8.0)) index_qsq = 6;
            else if ((8.0 < invM * invM) && (invM * invM < 11.0)) index_qsq = 7;
            else if ((11.0 < invM * invM) && (invM * invM < 12.5)) index_qsq = 8;
            else if ((12.5 < invM * invM) && (invM * invM < 15.0)) index_qsq = 9;
            else if ((15.0 < invM * invM) && (invM * invM < 17.0)) index_qsq = 10;
            else if ((17.0 < invM * invM) && (invM * invM < 19.0)) index_qsq = 11;
            else if (19.0 < invM * invM) index_qsq = 12;

            if ((MC_MXs < 0.644)) index_mass = 0;
            else if ((0.644 < MC_MXs) && (MC_MXs < 0.796)) index_mass = 1;
            else if ((0.796 < MC_MXs) && (MC_MXs < 0.996)) index_mass = 2;
            else if ((0.996 < MC_MXs) && (MC_MXs < 1.2)) index_mass = 3;
            else if ((1.2 < MC_MXs)) index_mass = 4;

            if (Get06GeVhigher1GeVunder_) {
                if ((0.6 < MC_MXs) && (MC_MXs < 1.0)) {}
                else continue;
            }

            Nevt = Nevt + FillTemplate(hist, MVA_var, total_weight * projection_multiplication, Bsig_M);

        }
        input_file->Close();

    }
    printf("%s has %lf events (with correction)\n", dirname, Nevt);

    return Nevt;
}

void LetsFillMK0star(const char* dirname, std::vector<Data>* list, const char* included_string, bool Get06GeVhigher1GeVunder_) {

    double invM = 0.0;
    double MXs = 0.0;

    std::vector<string> names;
    load_files_here(dirname, &names, included_string);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_Xs = (TTree*)input_file->Get("Xs");

        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);
        tree_Xs->SetBranchAddress("daughter__bo0__cm__spdaughter__bo0__cmM__bc__bc", &MXs);

        printf("%lld entries...\n", tree_Xs->GetEntries());
        for (unsigned int j = 0; j < tree_Xs->GetEntries(); j++) { // Fill
            tree_Xs->GetEntry(j);

            if (Get06GeVhigher1GeVunder_) {
                if ((0.6 < MXs) && (MXs < 1.0)) {}
                else continue;
            }

            Data temp_data = { invM * invM , MXs };
            list->push_back(temp_data);

        }
        input_file->Close();

    }

}

void FillNKstar(bool Get06GeVhigher1GeVunder_) {
    const char* K0star_dirname = "/home/belle2/junewoo/storage_ghi/20251024_Kappa_decayInfo/output/";
    const char* Kplusstar_dirname = "/home/belle2/junewoo/storage_ghi/20251024_Kappa_decayInfo/output/";
    const char* kappa0_dirname = "/home/belle2/junewoo/storage_ghi/20251024_Kappa_decayInfo/output/";
    const char* kappaplus_dirname = "/home/belle2/junewoo/storage_ghi/20251024_Kappa_decayInfo/output/";

    std::vector<Data> list;

    LetsFillMK0star(K0star_dirname, &list, "B02Kstar", Get06GeVhigher1GeVunder_);
    for (int i = 0; i < list.size(); i++) {
        int index_qsq = -1;
        int index_mass = -1;

        if ((list.at(i).q_square < 0.1)) index_qsq = 0;
        else if ((0.1 < list.at(i).q_square) && (list.at(i).q_square < 0.98)) index_qsq = 1;
        else if ((0.98 < list.at(i).q_square) && (list.at(i).q_square < 1.1)) index_qsq = 2;
        else if ((1.1 < list.at(i).q_square) && (list.at(i).q_square < 2.5)) index_qsq = 3;
        else if ((2.5 < list.at(i).q_square) && (list.at(i).q_square < 4.0)) index_qsq = 4;
        else if ((4.0 < list.at(i).q_square) && (list.at(i).q_square < 6.0)) index_qsq = 5;
        else if ((6.0 < list.at(i).q_square) && (list.at(i).q_square < 8.0)) index_qsq = 6;
        else if ((8.0 < list.at(i).q_square) && (list.at(i).q_square < 11.0)) index_qsq = 7;
        else if ((11.0 < list.at(i).q_square) && (list.at(i).q_square < 12.5)) index_qsq = 8;
        else if ((12.5 < list.at(i).q_square) && (list.at(i).q_square < 15.0)) index_qsq = 9;
        else if ((15.0 < list.at(i).q_square) && (list.at(i).q_square < 17.0)) index_qsq = 10;
        else if ((17.0 < list.at(i).q_square) && (list.at(i).q_square < 19.0)) index_qsq = 11;
        else if (19.0 < list.at(i).q_square) index_qsq = 12;

        if ((list.at(i).MKstar < 0.644)) index_mass = 0;
        else if ((0.644 < list.at(i).MKstar) && (list.at(i).MKstar < 0.796)) index_mass = 1;
        else if ((0.796 < list.at(i).MKstar) && (list.at(i).MKstar < 0.996)) index_mass = 2;
        else if ((0.996 < list.at(i).MKstar) && (list.at(i).MKstar < 1.2)) index_mass = 3;
        else if ((1.2 < list.at(i).MKstar)) index_mass = 4;

        N_K0star[index_qsq][index_mass] = N_K0star[index_qsq][index_mass] + Scale_K0star_here;
    }
    list.clear();

    LetsFillMK0star(Kplusstar_dirname, &list, "B2Kstar", Get06GeVhigher1GeVunder_);
    for (int i = 0; i < list.size(); i++) {
        int index_qsq = -1;
        int index_mass = -1;

        if ((list.at(i).q_square < 0.1)) index_qsq = 0;
        else if ((0.1 < list.at(i).q_square) && (list.at(i).q_square < 0.98)) index_qsq = 1;
        else if ((0.98 < list.at(i).q_square) && (list.at(i).q_square < 1.1)) index_qsq = 2;
        else if ((1.1 < list.at(i).q_square) && (list.at(i).q_square < 2.5)) index_qsq = 3;
        else if ((2.5 < list.at(i).q_square) && (list.at(i).q_square < 4.0)) index_qsq = 4;
        else if ((4.0 < list.at(i).q_square) && (list.at(i).q_square < 6.0)) index_qsq = 5;
        else if ((6.0 < list.at(i).q_square) && (list.at(i).q_square < 8.0)) index_qsq = 6;
        else if ((8.0 < list.at(i).q_square) && (list.at(i).q_square < 11.0)) index_qsq = 7;
        else if ((11.0 < list.at(i).q_square) && (list.at(i).q_square < 12.5)) index_qsq = 8;
        else if ((12.5 < list.at(i).q_square) && (list.at(i).q_square < 15.0)) index_qsq = 9;
        else if ((15.0 < list.at(i).q_square) && (list.at(i).q_square < 17.0)) index_qsq = 10;
        else if ((17.0 < list.at(i).q_square) && (list.at(i).q_square < 19.0)) index_qsq = 11;
        else if (19.0 < list.at(i).q_square) index_qsq = 12;

        if ((list.at(i).MKstar < 0.644)) index_mass = 0;
        else if ((0.644 < list.at(i).MKstar) && (list.at(i).MKstar < 0.796)) index_mass = 1;
        else if ((0.796 < list.at(i).MKstar) && (list.at(i).MKstar < 0.996)) index_mass = 2;
        else if ((0.996 < list.at(i).MKstar) && (list.at(i).MKstar < 1.2)) index_mass = 3;
        else if ((1.2 < list.at(i).MKstar)) index_mass = 4;

        N_Kplusstar[index_qsq][index_mass] = N_Kplusstar[index_qsq][index_mass] + Scale_Kplusstar_here;
    }
    list.clear();

    LetsFillMK0star(kappa0_dirname, &list, "B02Kappa", Get06GeVhigher1GeVunder_);
    for (int i = 0; i < list.size(); i++) {
        int index_qsq = -1;
        int index_mass = -1;

        if ((list.at(i).q_square < 0.1)) index_qsq = 0;
        else if ((0.1 < list.at(i).q_square) && (list.at(i).q_square < 0.98)) index_qsq = 1;
        else if ((0.98 < list.at(i).q_square) && (list.at(i).q_square < 1.1)) index_qsq = 2;
        else if ((1.1 < list.at(i).q_square) && (list.at(i).q_square < 2.5)) index_qsq = 3;
        else if ((2.5 < list.at(i).q_square) && (list.at(i).q_square < 4.0)) index_qsq = 4;
        else if ((4.0 < list.at(i).q_square) && (list.at(i).q_square < 6.0)) index_qsq = 5;
        else if ((6.0 < list.at(i).q_square) && (list.at(i).q_square < 8.0)) index_qsq = 6;
        else if ((8.0 < list.at(i).q_square) && (list.at(i).q_square < 11.0)) index_qsq = 7;
        else if ((11.0 < list.at(i).q_square) && (list.at(i).q_square < 12.5)) index_qsq = 8;
        else if ((12.5 < list.at(i).q_square) && (list.at(i).q_square < 15.0)) index_qsq = 9;
        else if ((15.0 < list.at(i).q_square) && (list.at(i).q_square < 17.0)) index_qsq = 10;
        else if ((17.0 < list.at(i).q_square) && (list.at(i).q_square < 19.0)) index_qsq = 11;
        else if (19.0 < list.at(i).q_square) index_qsq = 12;

        if ((list.at(i).MKstar < 0.644)) index_mass = 0;
        else if ((0.644 < list.at(i).MKstar) && (list.at(i).MKstar < 0.796)) index_mass = 1;
        else if ((0.796 < list.at(i).MKstar) && (list.at(i).MKstar < 0.996)) index_mass = 2;
        else if ((0.996 < list.at(i).MKstar) && (list.at(i).MKstar < 1.2)) index_mass = 3;
        else if ((1.2 < list.at(i).MKstar)) index_mass = 4;

        N_kappa0[index_qsq][index_mass] = N_kappa0[index_qsq][index_mass] + 1.0;
    }
    list.clear();

    LetsFillMK0star(kappaplus_dirname, &list, "B2Kappa", Get06GeVhigher1GeVunder_);
    for (int i = 0; i < list.size(); i++) {
        int index_qsq = -1;
        int index_mass = -1;

        if ((list.at(i).q_square < 0.1)) index_qsq = 0;
        else if ((0.1 < list.at(i).q_square) && (list.at(i).q_square < 0.98)) index_qsq = 1;
        else if ((0.98 < list.at(i).q_square) && (list.at(i).q_square < 1.1)) index_qsq = 2;
        else if ((1.1 < list.at(i).q_square) && (list.at(i).q_square < 2.5)) index_qsq = 3;
        else if ((2.5 < list.at(i).q_square) && (list.at(i).q_square < 4.0)) index_qsq = 4;
        else if ((4.0 < list.at(i).q_square) && (list.at(i).q_square < 6.0)) index_qsq = 5;
        else if ((6.0 < list.at(i).q_square) && (list.at(i).q_square < 8.0)) index_qsq = 6;
        else if ((8.0 < list.at(i).q_square) && (list.at(i).q_square < 11.0)) index_qsq = 7;
        else if ((11.0 < list.at(i).q_square) && (list.at(i).q_square < 12.5)) index_qsq = 8;
        else if ((12.5 < list.at(i).q_square) && (list.at(i).q_square < 15.0)) index_qsq = 9;
        else if ((15.0 < list.at(i).q_square) && (list.at(i).q_square < 17.0)) index_qsq = 10;
        else if ((17.0 < list.at(i).q_square) && (list.at(i).q_square < 19.0)) index_qsq = 11;
        else if (19.0 < list.at(i).q_square) index_qsq = 12;

        if ((list.at(i).MKstar < 0.644)) index_mass = 0;
        else if ((0.644 < list.at(i).MKstar) && (list.at(i).MKstar < 0.796)) index_mass = 1;
        else if ((0.796 < list.at(i).MKstar) && (list.at(i).MKstar < 0.996)) index_mass = 2;
        else if ((0.996 < list.at(i).MKstar) && (list.at(i).MKstar < 1.2)) index_mass = 3;
        else if ((1.2 < list.at(i).MKstar)) index_mass = 4;

        N_kappaplus[index_qsq][index_mass] = N_kappaplus[index_qsq][index_mass] + 1.0;
    }
    list.clear();
}

void CalculateScale() {
    // calculate scale factor
    scale[0][2] = (N_K0star[0][2] * (0.021 / (1 - 0.021))) / N_kappa0[0][2];
    scale[1][2] = (N_K0star[1][2] * (0.021 / (1 - 0.021))) / N_kappa0[1][2];
    scale[2][2] = (N_K0star[2][2] * (0.0825 / (1 - 0.0825))) / N_kappa0[2][2];
    scale[3][2] = (N_K0star[3][2] * (0.144 / (1 - 0.144))) / N_kappa0[3][2];
    scale[4][2] = (N_K0star[4][2] * (0.029 / (1 - 0.029))) / N_kappa0[4][2];
    scale[5][2] = (N_K0star[5][2] * (0.117 / (1 - 0.117))) / N_kappa0[5][2];
    scale[6][2] = (N_K0star[6][2] * (0.033 / (1 - 0.033))) / N_kappa0[6][2];
    scale[7][2] = (N_K0star[7][2] * (0.027 / (1 - 0.027))) / N_kappa0[7][2];
    scale[8][2] = (N_K0star[8][2] * (0.021 / (1 - 0.021))) / N_kappa0[8][2];
    scale[9][2] = (N_K0star[9][2] * (0.0065 / (1 - 0.0065))) / N_kappa0[9][2];
    scale[10][2] = (N_K0star[10][2] * ((-0.008) / (1 + 0.008))) / N_kappa0[10][2];
    scale[11][2] = (N_K0star[11][2] * (0.018 / (1 - 0.018))) / N_kappa0[11][2];
    scale[12][2] = (N_K0star[12][2] * (0.018 / (1 - 0.018))) / N_kappa0[12][2];

    scale[0][1] = scale[0][3] = (N_K0star[0][1] + N_K0star[0][2] + N_K0star[0][3]) / ((N_kappa0[0][1] + N_kappa0[0][3]) * (1.0 / 0.052 - 1)) - (N_kappa0[0][2] * scale[0][2]) / (N_kappa0[0][1] + N_kappa0[0][3]);
    scale[1][1] = scale[1][3] = (N_K0star[1][1] + N_K0star[1][2] + N_K0star[1][3]) / ((N_kappa0[1][1] + N_kappa0[1][3]) * (1.0 / 0.052 - 1)) - (N_kappa0[1][2] * scale[1][2]) / (N_kappa0[1][1] + N_kappa0[1][3]);
    scale[2][1] = scale[2][3] = (N_K0star[2][1] + N_K0star[2][2] + N_K0star[2][3]) / ((N_kappa0[2][1] + N_kappa0[2][3]) * (1.0 / 0.178 - 1)) - (N_kappa0[2][2] * scale[2][2]) / (N_kappa0[2][1] + N_kappa0[2][3]);
    scale[3][1] = scale[3][3] = (N_K0star[3][1] + N_K0star[3][2] + N_K0star[3][3]) / ((N_kappa0[3][1] + N_kappa0[3][3]) * (1.0 / 0.304 - 1)) - (N_kappa0[3][2] * scale[3][2]) / (N_kappa0[3][1] + N_kappa0[3][3]);
    scale[4][1] = scale[4][3] = (N_K0star[4][1] + N_K0star[4][2] + N_K0star[4][3]) / ((N_kappa0[4][1] + N_kappa0[4][3]) * (1.0 / 0.071 - 1)) - (N_kappa0[4][2] * scale[4][2]) / (N_kappa0[4][1] + N_kappa0[4][3]);
    scale[5][1] = scale[5][3] = (N_K0star[5][1] + N_K0star[5][2] + N_K0star[5][3]) / ((N_kappa0[5][1] + N_kappa0[5][3]) * (1.0 / 0.254 - 1)) - (N_kappa0[5][2] * scale[5][2]) / (N_kappa0[5][1] + N_kappa0[5][3]);
    scale[6][1] = scale[6][3] = (N_K0star[6][1] + N_K0star[6][2] + N_K0star[6][3]) / ((N_kappa0[6][1] + N_kappa0[6][3]) * (1.0 / 0.082 - 1)) - (N_kappa0[6][2] * scale[6][2]) / (N_kappa0[6][1] + N_kappa0[6][3]);
    scale[7][1] = scale[7][3] = (N_K0star[7][1] + N_K0star[7][2] + N_K0star[7][3]) / ((N_kappa0[7][1] + N_kappa0[7][3]) * (1.0 / 0.0655 - 1)) - (N_kappa0[7][2] * scale[7][2]) / (N_kappa0[7][1] + N_kappa0[7][3]);
    scale[8][1] = scale[8][3] = (N_K0star[8][1] + N_K0star[8][2] + N_K0star[8][3]) / ((N_kappa0[8][1] + N_kappa0[8][3]) * (1.0 / 0.049 - 1)) - (N_kappa0[8][2] * scale[8][2]) / (N_kappa0[8][1] + N_kappa0[8][3]);
    scale[9][1] = scale[9][3] = (N_K0star[9][1] + N_K0star[9][2] + N_K0star[9][3]) / ((N_kappa0[9][1] + N_kappa0[9][3]) * (1.0 / 0.0165 - 1)) - (N_kappa0[9][2] * scale[9][2]) / (N_kappa0[9][1] + N_kappa0[9][3]);
    scale[10][1] = scale[10][3] = (N_K0star[10][1] + N_K0star[10][2] + N_K0star[10][3]) / ((N_kappa0[10][1] + N_kappa0[10][3]) * (1.0 / -0.016 - 1)) - (N_kappa0[10][2] * scale[10][2]) / (N_kappa0[10][1] + N_kappa0[10][3]);
    scale[11][1] = scale[11][3] = (N_K0star[11][1] + N_K0star[11][2] + N_K0star[11][3]) / ((N_kappa0[11][1] + N_kappa0[11][3]) * (1.0 / 0.034 - 1)) - (N_kappa0[11][2] * scale[11][2]) / (N_kappa0[11][1] + N_kappa0[11][3]);
    scale[12][1] = scale[12][3] = (N_K0star[12][1] + N_K0star[12][2] + N_K0star[12][3]) / ((N_kappa0[12][1] + N_kappa0[12][3]) * (1.0 / 0.034 - 1)) - (N_kappa0[12][2] * scale[12][2]) / (N_kappa0[12][1] + N_kappa0[12][3]);

    for (int i = 0; i < 13; i++) {
        scale[i][0] = scale[i][1];
        scale[i][4] = scale[i][3];
    }
}

void InitNKstar() {
    for (int i = 0; i < 13; i++) {
        for (int j = 0; j < 5; j++) {
            N_Kplusstar[i][j] = 0.0;
            N_K0star[i][j] = 0.0;
            N_kappaplus[i][j] = 0.0;
            N_kappa0[i][j] = 0.0;
        }
    }
}

void GetRatioPlot(TH1D* first_, TH1D* second_, const char* name_) {
    TH1D* Signal_ratio = new TH1D("Signal_ratio", ";;non-official/official", RarityBins, BinMIN, BinMAX);

    Signal_ratio->Divide(first_, second_);

    Signal_ratio->SetLineColor(kBlack); Signal_ratio->SetMarkerStyle(21); Signal_ratio->Sumw2(); Signal_ratio->SetStats(0);

    TCanvas* c_temp = new TCanvas("c", "", 800, 800); c_temp->cd();

    TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.35, 1.0, 1.0);
    pad1->SetBottomMargin(0.08); pad1->SetLeftMargin(0.15);
    pad1->SetGridx(); pad1->Draw(); pad1->cd();

    gStyle->SetPalette(kPastel);

    Float_t ymax_1 = first_->GetMaximum();
    Float_t ymax_2 = second_->GetMaximum();
    double real_max = 0;
    if (ymax_1 > ymax_2) real_max = ymax_1;
    else real_max = ymax_2;

    first_->SetMaximum(real_max * 1.1);

    second_->SetStats(0);
    first_->SetStats(0);

    second_->SetFillColor(4); second_->SetLineColor(4); second_->SetLineWidth(1); second_->SetFillStyle(3002); second_->Draw("e1 Hist");
    first_->SetFillColor(2); first_->SetLineColor(2); first_->SetLineWidth(1); first_->SetFillStyle(3003); first_->Draw("e1 Hist SAME");
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
    c_temp->SaveAs(name_);

    printf("ratio: ");
    for (int i = 0; i < RarityBins; i++) {
        printf("%f,", 1.0 / Signal_ratio->GetBinContent(i + 1));
    }
    printf("\n");
}

void GetSwaveFraction_advanced(){

    FillNKstar(false);

    CalculateScale();

    InitNKstar();

    FillNKstar(true);

    double sum_K0 = 0.0;
    double sum_kappa0 = 0.0;

    for (int i = 0; i < 13; i++) {
        sum_K0 = sum_K0 + N_K0star[i][3];
        sum_kappa0 = sum_kappa0 + scale[i][3] * N_kappa0[i][3];
    }

    printf("%lf\n", (sum_kappa0)/(sum_kappa0 + sum_K0));

    const char* dir_Ntuple = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Swave/SIGNAL_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut";
    double** N_Kstar_after = GetPDFs(dir_Ntuple, "B2Kstar", "Bplus", "SIGNAL", true, Scale_Kplusstar_here, "otherwise");
    double** N_K0star_after = GetPDFs(dir_Ntuple, "B02Kstar", "Bzero", "SIGNAL", true, Scale_K0star_here, "otherwise");
    double** N_kappa_after = GetPDFs(dir_Ntuple, "B2Kappa", "Bplus", "SIGNAL", true, 1.0, "kappa");
    double** N_kappa0_after = GetPDFs(dir_Ntuple, "B02Kappa", "Bzero", "SIGNAL", true, 1.0, "kappa");

    /* =============================== 2D distribution =============================== */
    for (int i = 0; i < 13; i++) {
        for (int j = 0; j < 5; j++) {
            printf("%lf", (N_Kstar_after[i][j] + N_K0star_after[i][j])/(N_Kplusstar[i][j] + N_K0star[i][j]));
            if (j != 4) printf(",");
        }
        printf("\n");
    }

    printf("\n");

    for (int i = 0; i < 13; i++) {
        for (int j = 0; j < 5; j++) {
            double p = (N_Kstar_after[i][j] + N_K0star_after[i][j]) / (N_Kplusstar[i][j] + N_K0star[i][j]);
            printf("%lf", std::sqrt(p* (1 - p) / ((N_Kplusstar[i][j] / Scale_Kplusstar_here) + (N_K0star[i][j] / Scale_K0star_here))));
            if (j != 4) printf(",");
        }
        printf("\n");
    }

    printf("\n");

    for (int i = 0; i < 13; i++) {
        for (int j = 0; j < 5; j++) {
            printf("%lf", (N_kappa_after[i][j] * (TBp / TB0) + N_kappa0_after[i][j]) / (scale[i][j] * N_kappaplus[i][j] * (TBp / TB0) + scale[i][j] * N_kappa0[i][j]));
            if (j != 4) printf(",");
        }
        printf("\n");
    }

    printf("\n");

    for (int i = 0; i < 13; i++) {
        for (int j = 0; j < 5; j++) {
            double p = (N_kappa_after[i][j] * (TBp / TB0) + N_kappa0_after[i][j]) / (scale[i][j] * N_kappaplus[i][j] * (TBp / TB0) + scale[i][j] * N_kappa0[i][j]);
            printf("%lf", std::sqrt(p* (1 - p) / (scale[i][j] * N_kappaplus[i][j] * (TBp / TB0) + scale[i][j] * N_kappa0[i][j])));
            if (j != 4) printf(",");
        }
        printf("\n");
    }

    printf("\n");

    /* =============================== 1D distribution =============================== */
    for (int i = 0; i < 13; i++) {
        double N_sum_Kstar_after = 0.0;
        double N_sum_K0_after = 0.0;
        double N_sum_Kstar = 0.0;
        double N_sum_K0 = 0.0;
        for (int j = 0; j < 5; j++) {
            N_sum_Kstar_after = N_sum_Kstar_after + N_Kstar_after[i][j];
            N_sum_K0_after = N_sum_K0_after + N_K0star_after[i][j];
            N_sum_Kstar = N_sum_Kstar + N_Kplusstar[i][j];
            N_sum_K0 = N_sum_K0 + N_K0star[i][j];
        }
        printf("%lf", (N_sum_Kstar_after + N_sum_K0_after) / (N_sum_Kstar + N_sum_K0));
        if (i != 12) printf(",");
        else printf("\n");
    }

    printf("\n");

    for (int i = 0; i < 13; i++) {
        double N_sum_Kstar_after = 0.0;
        double N_sum_K0_after = 0.0;
        double N_sum_Kstar = 0.0;
        double N_sum_K0 = 0.0;
        for (int j = 0; j < 5; j++) {
            N_sum_Kstar_after = N_sum_Kstar_after + N_Kstar_after[i][j];
            N_sum_K0_after = N_sum_K0_after + N_K0star_after[i][j];
            N_sum_Kstar = N_sum_Kstar + N_Kplusstar[i][j];
            N_sum_K0 = N_sum_K0 + N_K0star[i][j];
        }
        double p = (N_sum_Kstar_after + N_sum_K0_after) / (N_sum_Kstar + N_sum_K0);
        printf("%lf", std::sqrt(p * (1 - p) / ((N_sum_Kstar / Scale_Kplusstar_here) + (N_sum_K0 / Scale_K0star_here))));
        if (i != 12) printf(",");
        else printf("\n");
    }

    printf("\n");

    for (int i = 0; i < 13; i++) {
        double N_sum_kappa_after = 0.0;
        double N_sum_kappa0_after = 0.0;
        double N_sum_kappa = 0.0;
        double N_sum_kappa0 = 0.0;
        for (int j = 0; j < 5; j++) {
            N_sum_kappa_after = N_sum_kappa_after + N_kappa_after[i][j];
            N_sum_kappa0_after = N_sum_kappa0_after + N_kappa0_after[i][j];
            N_sum_kappa = N_sum_kappa + scale[i][j] * N_kappaplus[i][j] * (TBp / TB0);
            N_sum_kappa0 = N_sum_kappa0 + scale[i][j] * N_kappa0[i][j];
        }
        printf("%lf", (N_sum_kappa_after + N_sum_kappa0_after) / (N_sum_kappa + N_sum_kappa0));
        if (i != 12) printf(",");
        else printf("\n");
    }

    printf("\n");

    for (int i = 0; i < 13; i++) {
        double N_sum_kappa_after = 0.0;
        double N_sum_kappa0_after = 0.0;
        double N_sum_kappa = 0.0;
        double N_sum_kappa0 = 0.0;
        double N_sum_kappa_no_scale = 0.0;
        double N_sum_kappa0_no_scale = 0.0;
        for (int j = 0; j < 5; j++) {
            N_sum_kappa_after = N_sum_kappa_after + N_kappa_after[i][j];
            N_sum_kappa0_after = N_sum_kappa0_after + N_kappa0_after[i][j];
            N_sum_kappa = N_sum_kappa + scale[i][j] * N_kappaplus[i][j] * (TBp / TB0);
            N_sum_kappa0 = N_sum_kappa0 + scale[i][j] * N_kappa0[i][j];
            N_sum_kappa_no_scale = N_sum_kappa_no_scale + N_kappaplus[i][j] * (TBp / TB0);
            N_sum_kappa0_no_scale = N_sum_kappa0_no_scale + N_kappa0[i][j];
        }
        double p = (N_sum_kappa_after + N_sum_kappa0_after) / (N_sum_kappa + N_sum_kappa0);
        printf("%lf", std::sqrt(p * (1 - p) / (N_sum_kappa_no_scale + N_sum_kappa0_no_scale)));
        if (i != 12) printf(",");
        else printf("\n");
    }

    printf("\n");

    /* =============================== fraction =============================== */
    {
        double N_total_Kstar_after = 0.0;
        double N_total_kappa_after = 0.0;
        for (int i = 0; i < 13; i++) {
            for (int j = 0; j < 5; j++) {
                N_total_Kstar_after = N_total_Kstar_after + N_Kstar_after[i][j] + N_K0star_after[i][j];
                N_total_kappa_after = N_total_kappa_after + N_kappa_after[i][j] * (TBp / TB0) + N_kappa0_after[i][j];
            }
        }
        printf("Nkappa/NKstar = %lf\n", N_total_kappa_after / N_total_Kstar_after);
    }

    /* =============================== 2D distribution for Nevt =============================== */
    for (int i = 0; i < 13; i++) {
        for (int j = 0; j < 5; j++) {
            printf("%lf", (N_Kstar_after[i][j] + N_K0star_after[i][j]));
            if (j != 4) printf(",");
        }
        printf("\n");
    }

    printf("\n");

    for (int i = 0; i < 13; i++) {
        for (int j = 0; j < 5; j++) {
            printf("%lf", (N_Kplusstar[i][j] + N_K0star[i][j]));
            if (j != 4) printf(",");
        }
        printf("\n");
    }

    printf("\n");

    for (int i = 0; i < 13; i++) {
        for (int j = 0; j < 5; j++) {
            printf("%lf", (N_kappa_after[i][j] * (TBp / TB0) + N_kappa0_after[i][j]));
            if (j != 4) printf(",");
        }
        printf("\n");
    }

    printf("\n");

    for (int i = 0; i < 13; i++) {
        for (int j = 0; j < 5; j++) {
            printf("%lf", (scale[i][j] * N_kappaplus[i][j] * (TBp / TB0) + scale[i][j] * N_kappa0[i][j]));
            if (j != 4) printf(",");
        }
        printf("\n");
    }

    printf("\n");

    /* =============================== eff calculation =============================== */
    double total_N_sum_Kstar_after = 0.0;
    double total_N_sum_K0_after = 0.0;
    double total_N_sum_Kstar = 0.0;
    double total_N_sum_K0 = 0.0;
    for (int i = 0; i < 13; i++) {
        for (int j = 0; j < 5; j++) {
            total_N_sum_Kstar_after = total_N_sum_Kstar_after + N_Kstar_after[i][j];
            total_N_sum_K0_after = total_N_sum_K0_after + N_K0star_after[i][j];
            total_N_sum_Kstar = total_N_sum_Kstar + N_Kplusstar[i][j];
            total_N_sum_K0 = total_N_sum_K0 + N_K0star[i][j];
        }
    }
    printf("eff_Kstaronly = %.9lf, %lf, %lf\n", (total_N_sum_Kstar_after + total_N_sum_K0_after) / (total_N_sum_Kstar + total_N_sum_K0), (total_N_sum_Kstar_after + total_N_sum_K0_after), (total_N_sum_Kstar + total_N_sum_K0));

    double total_N_sum_kappa_after = 0.0;
    double total_N_sum_kappa0_after = 0.0;
    double total_N_sum_kappa = 0.0;
    double total_N_sum_kappa0 = 0.0;
    for (int i = 0; i < 13; i++) {
        for (int j = 0; j < 5; j++) {
            total_N_sum_kappa_after = total_N_sum_kappa_after + N_kappa_after[i][j];
            total_N_sum_kappa0_after = total_N_sum_kappa0_after + N_kappa0_after[i][j];
            total_N_sum_kappa = total_N_sum_kappa + scale[i][j] * N_kappaplus[i][j] * (TBp / TB0);
            total_N_sum_kappa0 = total_N_sum_kappa0 + scale[i][j] * N_kappa0[i][j];
        }
    }
    printf("eff_KstarANDkappa = %.9lf\n", (total_N_sum_Kstar_after + total_N_sum_K0_after + total_N_sum_kappa_after + total_N_sum_kappa0_after) / (total_N_sum_Kstar + total_N_sum_K0 + total_N_sum_kappa + total_N_sum_kappa0));

    /* =============================== get THD =============================== */
    TH1D* Kstar892 = new TH1D("Kstar892", "Kstar892", RarityBins, BinMIN, BinMAX);
    TH1D* Kstar700 = new TH1D("Kstar700+892", "Kstar700", RarityBins, BinMIN, BinMAX);
    TH1D* Signal_ratio = new TH1D("Signal_ratio", ";;896/700", RarityBins, BinMIN, BinMAX);

    GetPDFs(dir_Ntuple, "B2Kstar", Kstar892, "Bplus", "SIGNAL", true, Scale_Kplusstar_here, "otherwise");
    GetPDFs(dir_Ntuple, "B02Kstar", Kstar892, "Bzero", "SIGNAL", true, Scale_K0star_here, "otherwise");
    GetPDFs(dir_Ntuple, "B2Kappa", Kstar700, "Bplus", "SIGNAL", true, 1.0, "kappa");
    GetPDFs(dir_Ntuple, "B02Kappa", Kstar700, "Bzero", "SIGNAL", true, 1.0, "kappa");

    Kstar700->Add(Kstar700, Kstar892);

    double sum_892 = 0.0;
    double sum_700 = 0.0;
    for (int i = 0; i < RarityBins; i++) {
        sum_892 = sum_892 + Kstar892->GetBinContent(i + 1);
        sum_700 = sum_700 + Kstar700->GetBinContent(i + 1);
    }

    Kstar700->Scale(sum_892 / sum_700);

    GetRatioPlot(Kstar892, Kstar700, "ratio.png");

    printf("892 should be scaled to %lf\n", sum_700 / sum_892);
}
