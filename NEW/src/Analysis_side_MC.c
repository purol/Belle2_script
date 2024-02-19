// last update: 2023-07-14
// for Belle2 data

#include <iostream>
#include <stdio.h>
#include <string>
#include <algorithm>
#include <vector>
#include <numeric>
#include <time.h>
#include <random>
#include <sstream>
#include <queue>

#include "TH1.h"
#include "TH2.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TCollection.h"
#include "TPad.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TMath.h"

#include "correctors.h"
#include "constants.h"

/*
when you add new variables:
revise # define N_Needed_info ...
revise typedef struct dataLoader::ConvertIntoSeparateDataFile(std::string output_name, double flag = 0)
revise void Loader::GetData(TFile* input_file)
check MakeShiftDoubleToInt if it affects the code
revise void Loader::PrintRootFile(std::string output_name)
revise void Loader::PrintSeparateRootFile(std::string output_name)
revise void Loader::ConvertIntoSeparateDataFile(std::string output_name, double flag = 0)
*/

# define MCTYPE "MC15ri"

double ObtainWeight(const char* type, const char* MC_version, const char* category, std::string filename);

/* ====================================== */

void load_files(const char *dirname, std::vector<std::string>* names){
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

typedef struct data{
    int __experiment__;
    int __run__;
    unsigned int __event__;
    int __candidate__;
    int __ncandidates__;

    double Upsilon_decayID;
    double Bsig_decayID;

    int upsilon_experiment;
    int upsilon_run;
    unsigned int upsilon_event;
    int upsilon_candidate;
    int upsilon_ncandidates;
    int Bsig_experiment;
    int Bsig_run;
    unsigned int Bsig_event;
    int Bsig_candidate;
    int Bsig_ncandidates;
    int Btag_experiment;
    int Btag_run;
    unsigned int Btag_event;
    int Btag_candidate;
    int Btag_ncandidates;

    double Upsilon_info[N_Upsilon_info];
    // 0: number of ECL clusters in ROE(cleanMask), 1: number of KLM clusters in ROE
    // 2: number of tracks in ROE(cleanMask), 3: roeEextra(cleanMask)
    // 4: nROE_NeutralECLClusters(cleanMask), 5: nROE_K_S0, 6: nROE_pi0
    // 7: missing momentum of event theta, 8: missing momentum
    // 9: missing energy at CMS, 10: number of remaining tracks
    // 11: roeNeextra(cleanMask), 12: useCMSFrame__boroeNeextra__bocleanMask__bc__bc
    // 13: nROE_ParticlesInList__bogamma__clmygamma__bc
    // 14: foxWolframR1, 15: foxWolframR2, 16: foxWolframR3, 17: foxWolframR4
    // 18: harmonicMomentThrust0, 19: harmonicMomentThrust1, 20: harmonicMomentThrust2
    // 21: harmonicMomentThrust3, 22: harmonicMomentThrust4, 23: cleoConeThrust0
    // 24: cleoConeThrust1, 25: cleoConeThrust2
    // 26: cleoConeThrust3, 27: cleoConeThrust4, 28: cleoConeThrust5
    // 29: cleoConeThrust6, 30: cleoConeThrust7, 31: cleoConeThrust8
    // 32: sphericity, 33: aplanarity, 34: thrust, 35: thrustAxisCosTheta
    // 36: MsquaredBsig_op0, 37: MsquaredBsig_op1, 38: MsquaredBsig_op2
    // 39: MsquaredBsig_op3, 40: MsquaredBsig_op4, 41: MsquaredBsig_op7
    // 42: roeP__bocleanMask__bc, 43: roeM__bocleanMask__bc, 44: roePTheta__bocleanMask__bc
    // 45: qsquared - MB^2, 46: chiProb, 47: dr, 48: dz
    // 49: nElectron, 50: nMuon, 51: nElectronloose, 52: nMuonloose
    // 53: nElectrontight, 54: nMuontight, 55: beamE, 56: number of tracks in ROE(looseMask)
    // 57: Ecms, 58: Ngamma_v111, 59: Eecl_v111, 60: Eecl_v111_matched, 61: Eecl_v111_unmatched
    // 62: Ngamma_v133, 63: Eecl_v133, 64: Eecl_v133_matched, 65:Eecl_v133_unmatched
    // 66: Ngamma_v200, 67: Ngamma_v200_matched, 68: Ngamma_v200_unmatched
    // 69: Eecl_v200, 70: Eecl_v200_matched, 71: Eecl_v200_unmatched
    // 72-143: Eecl_systematics
    // 144: nDplustoKL0, 145: nDzerotoKL0
    // 146: nBplustoKnn, 147: inv_Knn, 148: nBplustoK*nn, 149: inv_K*nn, 150: nBzerotoK0nn, 151: inv_K0nn, 152: nBzerotoK0*nn, 153: inv_K0*nn
    // 154: nBplustoKpKLKL_all, 155: nBplustoKpKLKL_NR, 156: BplustoKpKLKL_M01, 157: BplustoKpKLKL_M02, 158: BplustoKpKLKL_M12
    // 159: nBzerotoKSKLKL_all, 160: nBzerotoKSKLKL_NR, 161: BzerotoKSKLKL_M01, 162: BzerotoKSKLKL_M02, 163: BzerotoKSKLKL_M12

    double Bsig_info[N_Bsig_info];
    // 0: Bsig_E, 1: Bsig_E_CMS, 2: Bsig_E_Recoil
    // 3: Bsig_p, 4: Bsig_p_CMS, 5: Bsig_p_Recoil
    // 6: M
    // 7: KS0_distance
    // 8: nKslow1, 9: nKslow2
    // 10: nKslow3, 11: nPislow1, 12: nPislow2
    // 13: nPislow3
    // 14: Bsig_R2, 15: Bsig_thrustBm, 16: Bsig_thrustOm
    // 17: Bsig_cosTBTO, 18: Bsig_cosTBz, 19: Bsig_KSFWVariables_et
    // 20: Bsig_KSFWVariables_mm2, 21: Bsig_KSFWVariables_hso00, 22: Bsig_KSFWVariables_hso01
    // 23: Bsig_KSFWVariables_hso02, 24: Bsig_KSFWVariables_hso03: 25: Bsig_KSFWVariables_hso04
    // 26: Bsig_KSFWVariables_hso10, 27: Bsig_KSFWVariables_hso12, 28: Bsig_KSFWVariables_hso14
    // 29: Bsig_KSFWVariables_hso20, 30: Bsig_KSFWVariables_hso22, 31: Bsig_KSFWVariables_hso24
    // 32: Bsig_KSFWVariables_hoo0, 33: Bsig_KSFWVariables_hoo1, 34: Bsig_KSFWVariables_hoo2
    // 35: Bsig_KSFWVariables_hoo3, 36: Bsig_KSFWVariables_hoo4, 37: Bsig_CleoConeCS_1
    // 38: Bsig_CleoConeCS_2, 39: Bsig_CleoConeCS_3, 40: Bsig_CleoConeCS_4
    // 41: Bsig_CleoConeCS_5, 42: Bsig_CleoConeCS_6, 43: Bsig_CleoConeCS_7
    // 44: Bsig_CleoConeCS_8, 45: Bsig_CleoConeCS_9, 46: Bsig_pt
    // 47: Bsig_useCMSFrame_pt, 48: Bsig_theta, 49: Bsig_useCMSFrame_theta
    // 50: nDc, 51: Dc_pValue_med, 52: Dc_pValue:std, 53: Dcsimpleveto_chiProb
    // 54: Dcsimpleveto_dr, 55: Dcsimpleveto_dz, 56: Dcsimpleveto_M
    // 57: nD0, 58: D0_pValue_med, 59: D0_pValue:std, 60: D0simpleveto_chiProb
    // 61: D0simpleveto_dr, 62: D0simpleveto_dz, 63: D0simpleveto_M
    // 64: mychiProb, 65: mydr, 66: mydz
    // 67-358: PIDs
    // 359-366: pi0s
    // 367-514: fake from E, 542-737: fake from MU

    double Btag_info[N_Btag_info];
    // 0: Btag_dmID, 1: Btag_Mbc, 2: Btag_deltaE
    // 3: Btag_E, 4: Btag_E_CMS, 5: Btag_signalprobability, 6: chiProb_tag
    // 7: dr, 8: dz, 9: Btag_p_CMS, 10: Btag_CMS_phi

    double Needed_info[N_Needed_info];
    // 0: Btag_R2, 1: Btag_thrustBm, 2: Btag_thrustOm. 3: Btag_cosTBTO
    // 4: Btag_cosTBz, 5: Btag_KSFWVariables(et), 6: Btag_KSFWVariables(mm2)
    // 7: Btag_KSFWVariables(hso00), 8: Btag_KSFWVariables(hso01), 9: Btag_KSFWVariables(hso02)
    // 10: Btag_KSFWVariables(hso03), 11: Btag_KSFWVariables(hso04)
    // 12: Btag_KSFWVariables(hso10), 13: Btag_KSFWVariables(hso12), 14: Btag_KSFWVariables(hso14)
    // 15: Btag_KSFWVariables(hso20), 16: Btag_KSFWVariables(hso22), 17: Btag_KSFWVariables(hso24)
    // 18: Btag_KSFWVariables(hoo0), 19: Btag_KSFWVariables(hoo1), 20: Btag_KSFWVariables(hoo2)
    // 21: Btag_KSFWVariables(hoo3), 22: Btag_KSFWVariables(hoo4)
    // 23: Btag_CleoConeCS(1), 24: Btag_CleoConeCS(2), 25: Btag_CleoConeCS(3)
    // 26: Btag_CleoConeCS(4), 27: Btag_CleoConeCS(5), 28: Btag_CleoConeCS(6)
    // 29: Btag_CleoConeCS(7), 30: Btag_CleoConeCS(8), 31:Btag_CleoConeCS(9)
    // 32: missing mass^2, 33: visible energy, 34: theta_CMS_of_tag
    // 35: decayhash, 36: decayhash_extended

    int Decay[N_decay]; // MC level info
    int Decay_nparticles[N_decay_nparticles]; // MC level info for # particle
    double Decay_syst_ff[N_decay_syst_ff]; // MC level info for helicity angle + q2

    float MVA_BB;
    float MVA_Continuum;

    int nROE_mu;
    int nROE_Jpsi;
    int nROE_Upsilon_BCS;
    int nROE_Upsilon;

} Data; 

class Loader {
public:
    enum Variable
    {
        Upsilon = 0,
        Bsig,
        Btag
    };
    enum Inequality
    {
        larger_than = 0,
        smaller_than
    };
    enum Arithmetic
    {
        plus = 0,
        minus,
        times,
        dividedby
    };
    enum BCS_criterion {
        Highest = 0,
        Lowest
    };
    enum Debug
    {
        DebugOn = 0,
        DebugOff
    };
    enum ValueOption {
        Linear = 0,
        Log
    };
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
    enum Qualifier {
        when = 0,
        except
    };
    enum Dvetotype {
        Dchargedwithoutpizero = 0,
        Dchargedwithpizero,
        Dneutralwithoutpizero,
        Dneutralwithpizero
    };
    enum MassRegion {
        SmallMass = 0,
        LargeMass
    };

private:
    std::queue<Data> TotalData;
    std::vector<TH1F*> TH1Fs;
    std::vector<TH2F*> TH2Fs;
    int current_TH1F;
    int current_TH2F;

    std::vector<double> N_events;
    std::vector<double> N_candidates;
    std::vector<double> N_candidates_modes[Loader::MAX_NUM_DECAYMODE];
    std::vector<double> N_MC_modes[Loader::MAX_NUM_DECAYMODE_MC];
    std::vector<std::string> titles;
    int current_N_event;
    int current_N_candidate;

    std::vector<int> experiment_indices;
    std::vector<int> run_indices;
    std::vector<unsigned int> event_indices;
    std::vector<int> candidate_indices;
    std::vector<int> ncandidates_indices;
    int current_N_experiment_index;
    bool DebugIsOn;

    std::vector<TFile*> files;
    std::vector<TTree*> trees_upsilon;
    std::vector<TTree*> trees_Bsig;
    std::vector<TTree*> trees_Btag;
    std::vector<TTree*> trees_Xs;
    int current_file;
    bool AllOfThemHaveXsBranch;
    bool AllOfThemHaveMVAOutput;
    bool AllOfThemHaveJpsiOutput;

    std::vector<THStack*> THStacks;
    std::vector<TH1F*> TH1Fs_THStack[Loader::MAX_NUM_DECAYMODE];
    int current_THStack;

    int current_Confusion_matrix;
    double Confusion[Loader::MAX_NUM_DECAYMODE][Loader::MAX_NUM_DECAYMODE_MC]; // [reco][MC truth]
    double Confusion_square[Loader::MAX_NUM_DECAYMODE][Loader::MAX_NUM_DECAYMODE + 1]; // [reco][MC truth]
    bool Confusion_matrixIsOn;

    int current_FOM;
    double FOM_Matrix[Nstep][Nstep];
    bool FOMIsOn;

    int current_MCcount;
    double MCcount[Loader::MAX_NUM_DECAYMODE_MC];
    bool MCcountOn;

    int UpsilonExperimentToTree;
    int UpsilonRunToTree;
    unsigned int UpsilonEventToTree;
    int UpsilonCandidateToTree;
    int UpsilonNcandidatesToTree;
    int BsigExperimentToTree;
    int BsigRunToTree;
    unsigned int BsigEventToTree;
    int BsigCandidateToTree;
    int BsigNcandidatesToTree;
    int BtagExperimentToTree;
    int BtagRunToTree;
    unsigned int BtagEventToTree;
    int BtagCandidateToTree;
    int BtagNcandidatesToTree;
    double UpsilonDataToTree[N_Upsilon_info];
    double BsigDataToTree[N_Bsig_info];
    double BtagDataToTree[N_Btag_info];
    int DecayDataToTree[N_decay];
    int DecayNparticlesDataToTree[N_decay_nparticles];
    double DecaySystFFDataToTree[N_decay_syst_ff];
    double DataToTree[N_Needed_info];
    double Upsilon_decayIDToTree;
    double Bsig_decayIDToTree;
    float MVA_BB_DataToTree;
    float MVA_Continuum_DataToTree;
    int nROE_mu_DataToTree;
    int nROE_Jpsi_DataToTree;
    int nROE_Upsilon_BCS_DataToTree;
    int nROE_Upsilon_DataToTree;

    bool DoesItHaveXsBranch;
    bool DoesItHaveMVAOutput;
    bool DoesItHaveJpsiOutput;

    bool TrueIfDecayModeMatch(Data temp_data, Loader::DecayMode decaymode);
    bool TrueIfDecayModeMatch_MC(Data temp_data, Loader::DecayModeMC decaymodeMC);

public:
    Loader();
    void initialize();
    void GetData(TFile* input_file);
    bool event_info_is_valid();
    void DrawTH1F(const char* name, const char* title, int nbins, double x_low, double x_high, Loader::Variable variable, int i, Loader::ValueOption dr = Loader::Linear);
    void DrawTH1F(const char* name, const char* title, int nbins, double x_low, double x_high, Loader::Variable variable, int i, Loader::Qualifier qualifier, Loader::DecayMode decaymode, Loader::ValueOption dr = Loader::Linear);
    void DrawTH2F(const char* name, const char* title, int nbinsx, double xlow, double xup, int nbinsy, double ylow, double yup, Loader::Variable variable_1, int i, Loader::Variable variable_2, int j);
    void DrawTH2F(const char* name, const char* title, int nbinsx, double xlow, double xup, int nbinsy, double ylow, double yup, Loader::Variable variable_1, int i, Loader::Variable variable_2, int j, Loader::Qualifier qualifier, Loader::DecayMode decaymode);
    void DrawTHStack(const char* name, const char* title, int nbins, double x_low, double x_high, Loader::Variable variable, int i, Loader::ValueOption dr = Loader::Linear);
    void PrintInformation(std::string title, std::string filename, const char* type, const char* MC_version, const char* category, bool smartmode = true);
    void Cut(Loader::Variable variable, int i, Loader::Inequality inq, double value);
    void Cut(Loader::Variable variable_1, int i_1, Loader::Arithmetic ari, Loader::Variable variable_2, int i_2, Loader::Inequality inq, double value);
    void Cut(const char* bracket_1, double const_1, Loader::Arithmetic ari_1, Loader::Variable variable_1, int i_1, const char* bracket_2, Loader::Arithmetic ari, const char* bracket_3, double const_2, Loader::Arithmetic ari_2, Loader::Variable variable_2, int i_2, const char* bracket_4, Loader::Inequality inq, double value);
    void Cut(Loader::Variable variable, int i, Loader::Inequality inq, double value, Loader::Qualifier qualifier, Loader::DecayMode decaymode);
    void Cut(Loader::Variable variable, int i, Loader::Inequality inq, double value, Loader::Qualifier qualifier, Loader::Variable variable_qual, int i_qual, Loader::Inequality inq_qual, double value_qual);
    void BCS(Loader::Variable variable, int index, Loader::BCS_criterion crit);
    bool IsBCSValid();
    void End();
    void PrintRootFile(std::string output_name);
    void PrintSeparateRootFile(std::string output_name);
    void ConvertIntoSeparateDataFile(std::string output_name, int flag = 0);
    void PrintDebugLogIf(Loader::Variable variable, int i, Loader::Inequality inq, double value);
    void PrintConfusionMatrix(std::string filename, const char* type, const char* MC_version, const char* category, bool smartmode = true);
    void DvetoFor(Loader::Variable variable, int i, double min, double max);
    void BsigFitConvergeFor(Loader::Variable variable, int i);
    void OnlySelectDvetoTypeFor(Loader::Variable variable, int Dchargedvetomassindex, int DchargedvetodmIDindex, int Dneutralvetomassindex, int DneutralvetodmIDindex, Loader::Dvetotype type);
    void DvetoAboutSpecificTypeFor(Loader::Variable variable, int Dchargedvetomassindex, int DchargedvetodmIDindex, int Dneutralvetomassindex, int DneutralvetodmIDindex, Loader::Dvetotype type, double minM, double maxM);
    void PrintFOM(std::string filename, const char* type, const char* MC_version, const char* category, bool smartmode = true);
    void PrintFOM1D(std::string filename, const char* type, const char* MC_version, const char* category, bool smartmode = true);
    void MVACut(double OBB, double Oqq, Loader::MassRegion massRegion);
    void CountMCEvent(std::string filename, const char* type, const char* MC_version, const char* category, bool smartmode = true);
    void SelectDecayModeOf(Loader::DecayMode decaymode);
    void RejectDecayModeOf(Loader::DecayMode decaymode);
    void BeamEnergyCorrectionFromDeltaE(int index_pBcms, int index_EBcms, int index_Mbc, int index_deltaE, double targetEbeamstar, bool IsItBtag);
};

Loader::Loader() {
    std::queue<Data> empty;
    std::swap(TotalData, empty);
    current_TH1F = 0;
    current_TH2F = 0;
    current_N_event = 0;
    current_N_candidate = 0;
    current_file = 0;
    current_N_experiment_index = 0;
    DebugIsOn = false;
    current_THStack = 0;
    DoesItHaveXsBranch = false;
    DoesItHaveMVAOutput = false;
    DoesItHaveJpsiOutput = false;
    AllOfThemHaveXsBranch = true;
    AllOfThemHaveMVAOutput = true;
    AllOfThemHaveJpsiOutput = true;
    current_Confusion_matrix = 0;
    for (int i = 0; i < Loader::MAX_NUM_DECAYMODE; i++) { // initialization
        for (int j = 0; j < Loader::MAX_NUM_DECAYMODE_MC; j++) {
            Confusion[i][j] = 0;
        }
    }
    for (int i = 0; i < Loader::MAX_NUM_DECAYMODE; i++) { // initialization
        for (int j = 0; j < Loader::MAX_NUM_DECAYMODE + 1; j++) {
            Confusion_square[i][j] = 0;
        }
    }
    Confusion_matrixIsOn = false;
    current_FOM = 0;
    for (int i = 0; i < Nstep; i++) for (int j = 0; j < Nstep; j++) FOM_Matrix[i][j] = 0.0; // initialization
    FOMIsOn = false;
    current_MCcount = 0;
    for (int i = 0; i < Loader::MAX_NUM_DECAYMODE_MC; i++) MCcount[i] = 0;
    MCcountOn = false;
}

void Loader::initialize() {
    std::queue<Data> empty;
    std::swap(TotalData, empty);
    current_TH1F = 0;
    current_TH2F = 0;
    current_N_event = 0;
    current_N_candidate = 0;
    current_file = 0;
    current_N_experiment_index = 0;
    DebugIsOn = false;
    current_THStack = 0;
    DoesItHaveXsBranch = false;
    DoesItHaveMVAOutput = false;
    DoesItHaveJpsiOutput = false;
    current_Confusion_matrix = 0;
    Confusion_matrixIsOn = false;
    current_FOM = 0;
    FOMIsOn = false;
    current_MCcount = 0;
    MCcountOn = false;
}

void Loader::GetData(TFile* input_file) {
    DoesItHaveXsBranch = false;
    TList* branches = input_file->GetListOfKeys();
    for (int i = 0; i < input_file->GetNkeys(); i++) {
        TKey* branch = (TKey*)branches->At(i);
        if (std::string("Xs") == branch->GetName()) { DoesItHaveXsBranch = true; }
    }
    if (DoesItHaveXsBranch == false) AllOfThemHaveXsBranch = false;

    TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
    TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
    TTree* tree_Btag = (TTree*)input_file->Get("Btag");

    TTree* tree_Xs;
    if (DoesItHaveXsBranch) tree_Xs = (TTree*)input_file->Get("Xs");
    else tree_Xs = nullptr;

    DoesItHaveMVAOutput = false;
    if (tree_upsilon->FindLeaf("MVA_BB") == 0 || tree_upsilon->FindLeaf("MVA_Continuum") == 0) DoesItHaveMVAOutput = false;
    else DoesItHaveMVAOutput = true;
    if (DoesItHaveMVAOutput == false) AllOfThemHaveMVAOutput = false;

    DoesItHaveJpsiOutput = false;
    if (tree_upsilon->FindLeaf("nParticlesInList__bomu__pl__clfromUpsilonmychargedMuon__bc") == 0 || tree_upsilon->FindLeaf("nParticlesInList__boJ__slpsi__clfromUpsilontemp__bc") == 0 || tree_upsilon->FindLeaf("nParticlesInList__boUpsilon__bo4S__bc__cltemp_withoutneutrino_BCS__bc") == 0 || tree_upsilon->FindLeaf("nParticlesInList__boUpsilon__bo4S__bc__cltemp_withoutneutrino__bc") == 0) DoesItHaveJpsiOutput = false;
    else DoesItHaveJpsiOutput = true;
    if (DoesItHaveJpsiOutput == false) AllOfThemHaveJpsiOutput = false;

    Data temp = { 0 };

    int MakeShiftDoubleToInt[17] = { 0 }; // intermediate variable to convert from int to double

    // get event_info
    tree_upsilon->SetBranchAddress("__experiment__", &temp.upsilon_experiment);
    tree_upsilon->SetBranchAddress("__run__", &temp.upsilon_run);
    tree_upsilon->SetBranchAddress("__event__", &temp.upsilon_event);
    tree_upsilon->SetBranchAddress("__candidate__", &temp.upsilon_candidate);
    tree_upsilon->SetBranchAddress("__ncandidates__", &temp.upsilon_ncandidates);
    tree_Bsig->SetBranchAddress("__experiment__", &temp.Bsig_experiment);
    tree_Bsig->SetBranchAddress("__run__", &temp.Bsig_run);
    tree_Bsig->SetBranchAddress("__event__", &temp.Bsig_event);
    tree_Bsig->SetBranchAddress("__candidate__", &temp.Bsig_candidate);
    tree_Bsig->SetBranchAddress("__ncandidates__", &temp.Bsig_ncandidates);
    tree_Btag->SetBranchAddress("__experiment__", &temp.Btag_experiment);
    tree_Btag->SetBranchAddress("__run__", &temp.Btag_run);
    tree_Btag->SetBranchAddress("__event__", &temp.Btag_event);
    tree_Btag->SetBranchAddress("__candidate__", &temp.Btag_candidate);
    tree_Btag->SetBranchAddress("__ncandidates__", &temp.Btag_ncandidates);

    // get decaymodeID
    tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &temp.Upsilon_decayID);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &temp.Bsig_decayID);

    // get Upsilon_info
    tree_upsilon->SetBranchAddress("nROE_ECLClusters__bocleanMask__bc", &temp.Upsilon_info[0]);
    tree_upsilon->SetBranchAddress("nROE_KLMClusters", &temp.Upsilon_info[1]);
    tree_upsilon->SetBranchAddress("nROE_Tracks__bocleanMask__bc", &temp.Upsilon_info[2]);
    tree_upsilon->SetBranchAddress("roeEextra__bocleanMask__bc", &temp.Upsilon_info[3]);
    tree_upsilon->SetBranchAddress("nROE_NeutralECLClusters__bocleanMask__bc", &temp.Upsilon_info[4]);
    tree_upsilon->SetBranchAddress("nROE_ParticlesInList__boK_S0__clmyKaonshort__bc", &temp.Upsilon_info[5]);
    tree_upsilon->SetBranchAddress("nROE_ParticlesInList__bopi0__clmyneutralPion__bc", &temp.Upsilon_info[6]);
    tree_upsilon->SetBranchAddress("missingMomentumOfEvent_theta", &temp.Upsilon_info[7]);
    tree_upsilon->SetBranchAddress("missingMomentumOfEvent", &temp.Upsilon_info[8]);
    tree_upsilon->SetBranchAddress("missingEnergyOfEventCMS", &temp.Upsilon_info[9]);
    tree_upsilon->SetBranchAddress("nRemainingTracksInEvent", &MakeShiftDoubleToInt[0]); // temp.Upsilon_info[10]
    tree_upsilon->SetBranchAddress("roeNeextra__bocleanMask__bc", &temp.Upsilon_info[11]);
    tree_upsilon->SetBranchAddress("useCMSFrame__boroeNeextra__bocleanMask__bc__bc", &temp.Upsilon_info[12]);
    tree_upsilon->SetBranchAddress("nROE_ParticlesInList__bogamma__clmygamma__bc", &temp.Upsilon_info[13]);
    tree_upsilon->SetBranchAddress("foxWolframR1", &temp.Upsilon_info[14]);
    tree_upsilon->SetBranchAddress("foxWolframR2", &temp.Upsilon_info[15]);
    tree_upsilon->SetBranchAddress("foxWolframR3", &temp.Upsilon_info[16]);
    tree_upsilon->SetBranchAddress("foxWolframR4", &temp.Upsilon_info[17]);
    tree_upsilon->SetBranchAddress("harmonicMomentThrust0", &temp.Upsilon_info[18]);
    tree_upsilon->SetBranchAddress("harmonicMomentThrust1", &temp.Upsilon_info[19]);
    tree_upsilon->SetBranchAddress("harmonicMomentThrust2", &temp.Upsilon_info[20]);
    tree_upsilon->SetBranchAddress("harmonicMomentThrust3", &temp.Upsilon_info[21]);
    tree_upsilon->SetBranchAddress("harmonicMomentThrust4", &temp.Upsilon_info[22]);
    tree_upsilon->SetBranchAddress("cleoConeThrust0", &temp.Upsilon_info[23]);
    tree_upsilon->SetBranchAddress("cleoConeThrust1", &temp.Upsilon_info[24]);
    tree_upsilon->SetBranchAddress("cleoConeThrust2", &temp.Upsilon_info[25]);
    tree_upsilon->SetBranchAddress("cleoConeThrust3", &temp.Upsilon_info[26]);
    tree_upsilon->SetBranchAddress("cleoConeThrust4", &temp.Upsilon_info[27]);
    tree_upsilon->SetBranchAddress("cleoConeThrust5", &temp.Upsilon_info[28]);
    tree_upsilon->SetBranchAddress("cleoConeThrust6", &temp.Upsilon_info[29]);
    tree_upsilon->SetBranchAddress("cleoConeThrust7", &temp.Upsilon_info[30]);
    tree_upsilon->SetBranchAddress("cleoConeThrust8", &temp.Upsilon_info[31]);
    tree_upsilon->SetBranchAddress("sphericity", &temp.Upsilon_info[32]);
    tree_upsilon->SetBranchAddress("aplanarity", &temp.Upsilon_info[33]);
    tree_upsilon->SetBranchAddress("thrust", &temp.Upsilon_info[34]);
    tree_upsilon->SetBranchAddress("thrustAxisCosTheta", &temp.Upsilon_info[35]);
    tree_upsilon->SetBranchAddress("MsquaredBsig_op0", &temp.Upsilon_info[36]);
    tree_upsilon->SetBranchAddress("MsquaredBsig_op1", &temp.Upsilon_info[37]);
    tree_upsilon->SetBranchAddress("MsquaredBsig_op2", &temp.Upsilon_info[38]);
    tree_upsilon->SetBranchAddress("MsquaredBsig_op3", &temp.Upsilon_info[39]);
    tree_upsilon->SetBranchAddress("MsquaredBsig_op4", &temp.Upsilon_info[40]);
    tree_upsilon->SetBranchAddress("MsquaredBsig_op7", &temp.Upsilon_info[41]);
    tree_upsilon->SetBranchAddress("roeP__bocleanMask__bc", &temp.Upsilon_info[42]);
    tree_upsilon->SetBranchAddress("roeM__bocleanMask__bc", &temp.Upsilon_info[43]);
    tree_upsilon->SetBranchAddress("roePTheta__bocleanMask__bc", &temp.Upsilon_info[44]);
    tree_upsilon->SetBranchAddress("qsquared", &temp.Upsilon_info[45]);
    tree_upsilon->SetBranchAddress("chiProb", &temp.Upsilon_info[46]);
    tree_upsilon->SetBranchAddress("dr", &temp.Upsilon_info[47]);
    tree_upsilon->SetBranchAddress("dz", &temp.Upsilon_info[48]);
    tree_upsilon->SetBranchAddress("nParticlesInList__boe__pl__clElectronFBDT__bc", &MakeShiftDoubleToInt[1]); // temp.Upsilon_info[49]
    tree_upsilon->SetBranchAddress("nParticlesInList__bomu__pl__clMuonFBDT__bc", &MakeShiftDoubleToInt[2]); // temp.Upsilon_info[50]
    tree_upsilon->SetBranchAddress("nParticlesInList__boe__pl__clElectronFBDT_loose__bc", &MakeShiftDoubleToInt[3]); // temp.Upsilon_info[51]
    tree_upsilon->SetBranchAddress("nParticlesInList__bomu__pl__clMuonFBDT_loose__bc", &MakeShiftDoubleToInt[4]); // temp.Upsilon_info[52]
    tree_upsilon->SetBranchAddress("nParticlesInList__boe__pl__clElectronFBDT_tight__bc", &MakeShiftDoubleToInt[5]); // temp.Upsilon_info[53]
    tree_upsilon->SetBranchAddress("nParticlesInList__bomu__pl__clMuonFBDT_tight__bc", &MakeShiftDoubleToInt[6]); // temp.Upsilon_info[54]
    tree_upsilon->SetBranchAddress("beamE", &temp.Upsilon_info[55]);
    tree_upsilon->SetBranchAddress("nROE_Tracks__bolooseMask__bc", &temp.Upsilon_info[56]);
    tree_upsilon->SetBranchAddress("Ecms", &temp.Upsilon_info[57]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgamma__bc", &temp.Upsilon_info[58]);
    tree_upsilon->SetBranchAddress("extraInfo__boEecl__bc", &temp.Upsilon_info[59]);
    tree_upsilon->SetBranchAddress("extraInfo__boEecl_matched__bc", &temp.Upsilon_info[60]);
    tree_upsilon->SetBranchAddress("extraInfo__boEecl_unmatched__bc", &temp.Upsilon_info[61]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav133__bc", &temp.Upsilon_info[62]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv133__bc", &temp.Upsilon_info[63]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv133_matched__bc", &temp.Upsilon_info[64]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv133_unmatched__bc", &temp.Upsilon_info[65]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &temp.Upsilon_info[66]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_matched__bc", &temp.Upsilon_info[67]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_unmatched__bc", &temp.Upsilon_info[68]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200__bc", &temp.Upsilon_info[69]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_matched__bc", &temp.Upsilon_info[70]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_unmatched__bc", &temp.Upsilon_info[71]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_800__bc", &temp.Upsilon_info[72]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_800_matched__bc", &temp.Upsilon_info[73]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_800_unmatched__bc", &temp.Upsilon_info[74]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_800__bc", &temp.Upsilon_info[75]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_800_matched__bc", &temp.Upsilon_info[76]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_800_unmatched__bc", &temp.Upsilon_info[77]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_825__bc", &temp.Upsilon_info[78]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_825_matched__bc", &temp.Upsilon_info[79]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_825_unmatched__bc", &temp.Upsilon_info[80]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_825__bc", &temp.Upsilon_info[81]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_825_matched__bc", &temp.Upsilon_info[82]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_825_unmatched__bc", &temp.Upsilon_info[83]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_850__bc", &temp.Upsilon_info[84]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_850_matched__bc", &temp.Upsilon_info[85]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_850_unmatched__bc", &temp.Upsilon_info[86]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_850__bc", &temp.Upsilon_info[87]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_850_matched__bc", &temp.Upsilon_info[88]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_850_unmatched__bc", &temp.Upsilon_info[89]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_875__bc", &temp.Upsilon_info[90]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_875_matched__bc", &temp.Upsilon_info[91]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_875_unmatched__bc", &temp.Upsilon_info[92]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_875__bc", &temp.Upsilon_info[93]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_875_matched__bc", &temp.Upsilon_info[94]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_875_unmatched__bc", &temp.Upsilon_info[95]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_900__bc", &temp.Upsilon_info[96]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_900_matched__bc", &temp.Upsilon_info[97]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_900_unmatched__bc", &temp.Upsilon_info[98]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_900__bc", &temp.Upsilon_info[99]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_900_matched__bc", &temp.Upsilon_info[100]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_900_unmatched__bc", &temp.Upsilon_info[101]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_925__bc", &temp.Upsilon_info[102]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_925_matched__bc", &temp.Upsilon_info[103]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_925_unmatched__bc", &temp.Upsilon_info[104]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_925__bc", &temp.Upsilon_info[105]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_925_matched__bc", &temp.Upsilon_info[106]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_925_unmatched__bc", &temp.Upsilon_info[107]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_950__bc", &temp.Upsilon_info[108]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_950_matched__bc", &temp.Upsilon_info[109]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_950_unmatched__bc", &temp.Upsilon_info[110]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_950__bc", &temp.Upsilon_info[111]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_950_matched__bc", &temp.Upsilon_info[112]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_950_unmatched__bc", &temp.Upsilon_info[113]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_975__bc", &temp.Upsilon_info[114]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_975_matched__bc", &temp.Upsilon_info[115]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_975_unmatched__bc", &temp.Upsilon_info[116]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_975__bc", &temp.Upsilon_info[117]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_975_matched__bc", &temp.Upsilon_info[118]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_975_unmatched__bc", &temp.Upsilon_info[119]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_025__bc", &temp.Upsilon_info[120]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_025_matched__bc", &temp.Upsilon_info[121]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_025_unmatched__bc", &temp.Upsilon_info[122]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_025__bc", &temp.Upsilon_info[123]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_025_matched__bc", &temp.Upsilon_info[124]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_025_unmatched__bc", &temp.Upsilon_info[125]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_050__bc", &temp.Upsilon_info[126]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_050_matched__bc", &temp.Upsilon_info[127]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_050_unmatched__bc", &temp.Upsilon_info[128]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_050__bc", &temp.Upsilon_info[129]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_050_matched__bc", &temp.Upsilon_info[130]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_050_unmatched__bc", &temp.Upsilon_info[131]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_075__bc", &temp.Upsilon_info[132]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_075_matched__bc", &temp.Upsilon_info[133]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_075_unmatched__bc", &temp.Upsilon_info[134]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_075__bc", &temp.Upsilon_info[135]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_075_matched__bc", &temp.Upsilon_info[136]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_075_unmatched__bc", &temp.Upsilon_info[137]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_100__bc", &temp.Upsilon_info[138]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_100_matched__bc", &temp.Upsilon_info[139]);
    tree_upsilon->SetBranchAddress("extraInfo__boNgammav200_100_unmatched__bc", &temp.Upsilon_info[140]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_100__bc", &temp.Upsilon_info[141]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_100_matched__bc", &temp.Upsilon_info[142]);
    tree_upsilon->SetBranchAddress("extraInfo__boEeclv200_100_unmatched__bc", &temp.Upsilon_info[143]);
    tree_upsilon->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &MakeShiftDoubleToInt[7]); // temp.Upsilon_info[144]
    tree_upsilon->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &MakeShiftDoubleToInt[8]); // temp.Upsilon_info[145]
    tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &MakeShiftDoubleToInt[9]); // temp.Upsilon_info[146]
    tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &temp.Upsilon_info[147]);
    tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &MakeShiftDoubleToInt[10]); // temp.Upsilon_info[148]
    tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &temp.Upsilon_info[149]);
    tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &MakeShiftDoubleToInt[11]); // temp.Upsilon_info[150]
    tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &temp.Upsilon_info[151]);
    tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &MakeShiftDoubleToInt[12]); // temp.Upsilon_info[152]
    tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &temp.Upsilon_info[153]);
    tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_all__bc", &MakeShiftDoubleToInt[13]); // temp.Upsilon_info[154]
    tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_NR__bc", &MakeShiftDoubleToInt[14]); // temp.Upsilon_info[155]
    tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &temp.Upsilon_info[156]);
    tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &temp.Upsilon_info[157]);
    tree_upsilon->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &temp.Upsilon_info[158]);
    tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_all__bc", &MakeShiftDoubleToInt[15]); // temp.Upsilon_info[159]
    tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_NR__bc", &MakeShiftDoubleToInt[16]); // temp.Upsilon_info[160]
    tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &temp.Upsilon_info[161]);
    tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &temp.Upsilon_info[162]);
    tree_upsilon->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &temp.Upsilon_info[163]);

    // get Bsig_info
    tree_Bsig->SetBranchAddress("Bsig_E", &temp.Bsig_info[0]);
    tree_Bsig->SetBranchAddress("Bsig_useCMSFrame_E", &temp.Bsig_info[1]);
    tree_upsilon->SetBranchAddress("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc", &temp.Bsig_info[2]);
    tree_Bsig->SetBranchAddress("Bsig_p", &temp.Bsig_info[3]);
    tree_Bsig->SetBranchAddress("Bsig_useCMSFrame_p", &temp.Bsig_info[4]);
    tree_upsilon->SetBranchAddress("useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc", &temp.Bsig_info[5]);
    tree_Bsig->SetBranchAddress("Bsig_M", &temp.Bsig_info[6]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_3D_distance", &temp.Bsig_info[7]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nKslow1", &temp.Bsig_info[8]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nKslow2", &temp.Bsig_info[9]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nKslow3", &temp.Bsig_info[10]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nPislow1", &temp.Bsig_info[11]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nPislow2", &temp.Bsig_info[12]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nPislow3", &temp.Bsig_info[13]);
    tree_Bsig->SetBranchAddress("Bsig_R2", &temp.Bsig_info[14]);
    tree_Bsig->SetBranchAddress("Bsig_thrustBm", &temp.Bsig_info[15]);
    tree_Bsig->SetBranchAddress("Bsig_thrustOm", &temp.Bsig_info[16]);
    tree_Bsig->SetBranchAddress("Bsig_cosTBTO", &temp.Bsig_info[17]);
    tree_Bsig->SetBranchAddress("Bsig_cosTBz", &temp.Bsig_info[18]);
    tree_Bsig->SetBranchAddress("Bsig_KSFWVariables_et", &temp.Bsig_info[19]);
    tree_Bsig->SetBranchAddress("Bsig_KSFWVariables_mm2", &temp.Bsig_info[20]);
    tree_Bsig->SetBranchAddress("Bsig_KSFWVariables_hso00", &temp.Bsig_info[21]);
    tree_Bsig->SetBranchAddress("Bsig_KSFWVariables_hso01", &temp.Bsig_info[22]);
    tree_Bsig->SetBranchAddress("Bsig_KSFWVariables_hso02", &temp.Bsig_info[23]);
    tree_Bsig->SetBranchAddress("Bsig_KSFWVariables_hso03", &temp.Bsig_info[24]);
    tree_Bsig->SetBranchAddress("Bsig_KSFWVariables_hso04", &temp.Bsig_info[25]);
    tree_Bsig->SetBranchAddress("Bsig_KSFWVariables_hso10", &temp.Bsig_info[26]);
    tree_Bsig->SetBranchAddress("Bsig_KSFWVariables_hso12", &temp.Bsig_info[27]);
    tree_Bsig->SetBranchAddress("Bsig_KSFWVariables_hso14", &temp.Bsig_info[28]);
    tree_Bsig->SetBranchAddress("Bsig_KSFWVariables_hso20", &temp.Bsig_info[29]);
    tree_Bsig->SetBranchAddress("Bsig_KSFWVariables_hso22", &temp.Bsig_info[30]);
    tree_Bsig->SetBranchAddress("Bsig_KSFWVariables_hso24", &temp.Bsig_info[31]);
    tree_Bsig->SetBranchAddress("Bsig_KSFWVariables_hoo0", &temp.Bsig_info[32]);
    tree_Bsig->SetBranchAddress("Bsig_KSFWVariables_hoo1", &temp.Bsig_info[33]);
    tree_Bsig->SetBranchAddress("Bsig_KSFWVariables_hoo2", &temp.Bsig_info[34]);
    tree_Bsig->SetBranchAddress("Bsig_KSFWVariables_hoo3", &temp.Bsig_info[35]);
    tree_Bsig->SetBranchAddress("Bsig_KSFWVariables_hoo4", &temp.Bsig_info[36]);
    tree_Bsig->SetBranchAddress("Bsig_CleoConeCS_1", &temp.Bsig_info[37]);
    tree_Bsig->SetBranchAddress("Bsig_CleoConeCS_2", &temp.Bsig_info[38]);
    tree_Bsig->SetBranchAddress("Bsig_CleoConeCS_3", &temp.Bsig_info[39]);
    tree_Bsig->SetBranchAddress("Bsig_CleoConeCS_4", &temp.Bsig_info[40]);
    tree_Bsig->SetBranchAddress("Bsig_CleoConeCS_5", &temp.Bsig_info[41]);
    tree_Bsig->SetBranchAddress("Bsig_CleoConeCS_6", &temp.Bsig_info[42]);
    tree_Bsig->SetBranchAddress("Bsig_CleoConeCS_7", &temp.Bsig_info[43]);
    tree_Bsig->SetBranchAddress("Bsig_CleoConeCS_8", &temp.Bsig_info[44]);
    tree_Bsig->SetBranchAddress("Bsig_CleoConeCS_9", &temp.Bsig_info[45]);
    tree_Bsig->SetBranchAddress("Bsig_pt", &temp.Bsig_info[46]);
    tree_Bsig->SetBranchAddress("Bsig_useCMSFrame_pt", &temp.Bsig_info[47]);
    tree_Bsig->SetBranchAddress("Bsig_theta", &temp.Bsig_info[48]);
    tree_Bsig->SetBranchAddress("Bsig_useCMSFrame_theta", &temp.Bsig_info[49]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nDc", &temp.Bsig_info[50]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Dc_pValue_med", &temp.Bsig_info[51]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Dc_pValue_std", &temp.Bsig_info[52]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_chiProb", &temp.Bsig_info[53]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_dr", &temp.Bsig_info[54]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_dz", &temp.Bsig_info[55]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_M", &temp.Bsig_info[56]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nD0", &temp.Bsig_info[57]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_D0_pValue_med", &temp.Bsig_info[58]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_D0_pValue_std", &temp.Bsig_info[59]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_D0simpleveto_chiProb", &temp.Bsig_info[60]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_D0simpleveto_dr", &temp.Bsig_info[61]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_D0simpleveto_dz", &temp.Bsig_info[62]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_D0simpleveto_M", &temp.Bsig_info[63]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_mychiProb", &temp.Bsig_info[64]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_mydr", &temp.Bsig_info[65]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_mydz", &temp.Bsig_info[66]);
    for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
        tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKtruebin" + std::to_string(i_PID)).c_str(), &temp.Bsig_info[67 + 4 * i_PID]);
        tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKmisbin" + std::to_string(i_PID)).c_str(), &temp.Bsig_info[67 + 4 * i_PID + 1]);
        tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npitruebin" + std::to_string(i_PID)).c_str(), &temp.Bsig_info[67 + 4 * i_PID + 2]);
        tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npimisbin" + std::to_string(i_PID)).c_str(), &temp.Bsig_info[67 + 4 * i_PID + 3]);
    }
    for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npi0bin" + std::to_string(i_pi0)).c_str(), &temp.Bsig_info[359 + i_pi0]);
    for (int i_PID = 0; i_PID < N_fakeE_syst; i_PID++) {
        tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeEbin_n" + std::to_string(i_PID)).c_str(), &temp.Bsig_info[367 + 4 * i_PID]);
        tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeEbin_p" + std::to_string(i_PID)).c_str(), &temp.Bsig_info[367 + 4 * i_PID + 1]);
        tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeEbin_n" + std::to_string(i_PID)).c_str(), &temp.Bsig_info[367 + 4 * i_PID + 2]);
        tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeEbin_p" + std::to_string(i_PID)).c_str(), &temp.Bsig_info[367 + 4 * i_PID + 3]);
    }
    for (int i_PID = 0; i_PID < N_fakeMU_syst; i_PID++) {
        tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeMUbin_n" + std::to_string(i_PID)).c_str(), &temp.Bsig_info[542 + 4 * i_PID]);
        tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeMUbin_p" + std::to_string(i_PID)).c_str(), &temp.Bsig_info[542 + 4 * i_PID + 1]);
        tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeMUbin_n" + std::to_string(i_PID)).c_str(), &temp.Bsig_info[542 + 4 * i_PID + 2]);
        tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeMUbin_p" + std::to_string(i_PID)).c_str(), &temp.Bsig_info[542 + 4 * i_PID + 3]);
    }

    // get Btag_info
    tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &temp.Btag_info[0]);
    tree_Btag->SetBranchAddress("Btag_Mbc", &temp.Btag_info[1]);
    tree_Btag->SetBranchAddress("Btag_deltaE", &temp.Btag_info[2]);
    tree_Btag->SetBranchAddress("Btag_E", &temp.Btag_info[3]);
    tree_Btag->SetBranchAddress("Btag_useCMSFrame_E", &temp.Btag_info[4]);
    tree_Btag->SetBranchAddress("Btag_extraInfo_SignalProbability", &temp.Btag_info[5]);
    tree_Btag->SetBranchAddress("Btag_chiProb", &temp.Btag_info[6]);
    tree_Btag->SetBranchAddress("Btag_dr", &temp.Btag_info[7]);
    tree_Btag->SetBranchAddress("Btag_dz", &temp.Btag_info[8]);
    tree_Btag->SetBranchAddress("Btag_useCMSFrame_p", &temp.Btag_info[9]);
    tree_Btag->SetBranchAddress("Btag_useCMSFrame_phi", &temp.Btag_info[10]);

    // other information I need
    tree_Btag->SetBranchAddress("Btag_R2", &temp.Needed_info[0]);
    tree_Btag->SetBranchAddress("Btag_thrustBm", &temp.Needed_info[1]);
    tree_Btag->SetBranchAddress("Btag_thrustOm", &temp.Needed_info[2]);
    tree_Btag->SetBranchAddress("Btag_cosTBTO", &temp.Needed_info[3]);
    tree_Btag->SetBranchAddress("Btag_cosTBz", &temp.Needed_info[4]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_et", &temp.Needed_info[5]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_mm2", &temp.Needed_info[6]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hso00", &temp.Needed_info[7]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hso01", &temp.Needed_info[8]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hso02", &temp.Needed_info[9]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hso03", &temp.Needed_info[10]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hso04", &temp.Needed_info[11]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hso10", &temp.Needed_info[12]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hso12", &temp.Needed_info[13]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hso14", &temp.Needed_info[14]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hso20", &temp.Needed_info[15]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hso22", &temp.Needed_info[16]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hso24", &temp.Needed_info[17]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hoo0", &temp.Needed_info[18]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hoo1", &temp.Needed_info[19]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hoo2", &temp.Needed_info[20]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hoo3", &temp.Needed_info[21]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hoo4", &temp.Needed_info[22]);
    tree_Btag->SetBranchAddress("Btag_CleoConeCS_1", &temp.Needed_info[23]);
    tree_Btag->SetBranchAddress("Btag_CleoConeCS_2", &temp.Needed_info[24]);
    tree_Btag->SetBranchAddress("Btag_CleoConeCS_3", &temp.Needed_info[25]);
    tree_Btag->SetBranchAddress("Btag_CleoConeCS_4", &temp.Needed_info[26]);
    tree_Btag->SetBranchAddress("Btag_CleoConeCS_5", &temp.Needed_info[27]);
    tree_Btag->SetBranchAddress("Btag_CleoConeCS_6", &temp.Needed_info[28]);
    tree_Btag->SetBranchAddress("Btag_CleoConeCS_7", &temp.Needed_info[29]);
    tree_Btag->SetBranchAddress("Btag_CleoConeCS_8", &temp.Needed_info[30]);
    tree_Btag->SetBranchAddress("Btag_CleoConeCS_9", &temp.Needed_info[31]);
    tree_upsilon->SetBranchAddress("missingMass2OfEvent", &temp.Needed_info[32]);
    tree_upsilon->SetBranchAddress("visibleEnergyOfEventCMS", &temp.Needed_info[33]);
    tree_Btag->SetBranchAddress("Btag_useCMSFrame_theta", &temp.Needed_info[34]);
    tree_upsilon->SetBranchAddress("extraInfo__boDecayHash__bc", &temp.Needed_info[35]);
    tree_upsilon->SetBranchAddress("extraInfo__boDecayHashExtended__bc", &temp.Needed_info[36]);

    if (DoesItHaveXsBranch) {
        // decay mode (MC level)
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKcharge_total__bc", &temp.Decay[0]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &temp.Decay[1]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &temp.Decay[2]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCcomb__bc", &temp.Decay[3]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch1__bc", &temp.Decay[4]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch2__bc", &temp.Decay[5]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch3__bc", &temp.Decay[6]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch4__bc", &temp.Decay[7]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch5__bc", &temp.Decay[8]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch6__bc", &temp.Decay[9]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch7__bc", &temp.Decay[10]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch8__bc", &temp.Decay[11]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch9__bc", &temp.Decay[12]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch10__bc", &temp.Decay[13]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch11__bc", &temp.Decay[14]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch12__bc", &temp.Decay[15]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch13__bc", &temp.Decay[16]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch14__bc", &temp.Decay[17]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch15__bc", &temp.Decay[18]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKneutral_total__bc", &temp.Decay[19]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &temp.Decay[20]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &temp.Decay[21]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCcomb__bc", &temp.Decay[22]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch16__bc", &temp.Decay[23]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch17__bc", &temp.Decay[24]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch18__bc", &temp.Decay[25]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch19__bc", &temp.Decay[26]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch20__bc", &temp.Decay[27]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch21__bc", &temp.Decay[28]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch22__bc", &temp.Decay[29]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch23__bc", &temp.Decay[30]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch24__bc", &temp.Decay[31]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch25__bc", &temp.Decay[32]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch26__bc", &temp.Decay[33]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch27__bc", &temp.Decay[34]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch28__bc", &temp.Decay[35]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch29__bc", &temp.Decay[36]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch30__bc", &temp.Decay[37]);
        tree_Xs->SetBranchAddress("nParticlesInList__bonu_e__clMC_signal__bc", &temp.Decay_nparticles[0]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clMC_signal_total_e__bc", &temp.Decay_nparticles[1]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clMC_signal_total_e__bc", &temp.Decay_nparticles[2]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clPrimaryMC__bc", &temp.Decay_nparticles[3]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clPrimaryMC__bc", &temp.Decay_nparticles[4]);
        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &temp.Decay_syst_ff[index_q2]);
        tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", &temp.Decay_syst_ff[1]);
        tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", &temp.Decay_syst_ff[2]);
        tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &temp.Decay_syst_ff[index_MXs_Bc]);
        tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &temp.Decay_syst_ff[index_MXs_B0]);
        tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spM__bc", &temp.Decay_syst_ff[5]);
        tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spM__bc", &temp.Decay_syst_ff[6]);
    }
    else {
        for (int i = 0; i < N_decay; i++) temp.Decay[i] = -1;
        for (int i = 0; i < N_decay_nparticles; i++) temp.Decay_nparticles[i] = -1;
        for (int i = 0; i < N_decay_syst_ff; i++) temp.Decay_syst_ff[i] = -1;
    }

    if (DoesItHaveMVAOutput) {
        tree_upsilon->SetBranchAddress("MVA_BB", &temp.MVA_BB);
        tree_upsilon->SetBranchAddress("MVA_Continuum", &temp.MVA_Continuum);
    }
    else {
        temp.MVA_BB = -1.0f;
        temp.MVA_Continuum = -1.0f;
    }

    if (DoesItHaveJpsiOutput) {
        tree_upsilon->SetBranchAddress("nParticlesInList__bomu__pl__clfromUpsilonmychargedMuon__bc", &temp.nROE_mu);
        tree_upsilon->SetBranchAddress("nParticlesInList__boJ__slpsi__clfromUpsilontemp__bc", &temp.nROE_Jpsi);
        tree_upsilon->SetBranchAddress("nParticlesInList__boUpsilon__bo4S__bc__cltemp_withoutneutrino_BCS__bc", &temp.nROE_Upsilon_BCS);
        tree_upsilon->SetBranchAddress("nParticlesInList__boUpsilon__bo4S__bc__cltemp_withoutneutrino__bc", &temp.nROE_Upsilon);
    }
    else {
        temp.nROE_mu = -1;
        temp.nROE_Jpsi = -1;
        temp.nROE_Upsilon_BCS = -1;
        temp.nROE_Upsilon = -1;
    }

    printf("%lld entries...\n", tree_upsilon->GetEntries());
    for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
        tree_upsilon->GetEntry(j);
        tree_Bsig->GetEntry(j);
        tree_Btag->GetEntry(j);
        if (DoesItHaveXsBranch) tree_Xs->GetEntry(j);

        // Just Makeshift!
        temp.Upsilon_info[10] = static_cast<double>(MakeShiftDoubleToInt[0]);
        temp.Upsilon_info[49] = static_cast<double>(MakeShiftDoubleToInt[1]);
        temp.Upsilon_info[50] = static_cast<double>(MakeShiftDoubleToInt[2]);
        temp.Upsilon_info[51] = static_cast<double>(MakeShiftDoubleToInt[3]);
        temp.Upsilon_info[52] = static_cast<double>(MakeShiftDoubleToInt[4]);
        temp.Upsilon_info[53] = static_cast<double>(MakeShiftDoubleToInt[5]);
        temp.Upsilon_info[54] = static_cast<double>(MakeShiftDoubleToInt[6]);
        temp.Upsilon_info[144] = static_cast<double>(MakeShiftDoubleToInt[7]);
        temp.Upsilon_info[145] = static_cast<double>(MakeShiftDoubleToInt[8]);
        temp.Upsilon_info[146] = static_cast<double>(MakeShiftDoubleToInt[9]);
        temp.Upsilon_info[148] = static_cast<double>(MakeShiftDoubleToInt[10]);
        temp.Upsilon_info[150] = static_cast<double>(MakeShiftDoubleToInt[11]);
        temp.Upsilon_info[152] = static_cast<double>(MakeShiftDoubleToInt[12]);
        temp.Upsilon_info[154] = static_cast<double>(MakeShiftDoubleToInt[13]);
        temp.Upsilon_info[155] = static_cast<double>(MakeShiftDoubleToInt[14]);
        temp.Upsilon_info[159] = static_cast<double>(MakeShiftDoubleToInt[15]);
        temp.Upsilon_info[160] = static_cast<double>(MakeShiftDoubleToInt[16]);

        TotalData.push(temp);
    }
    input_file->Close();
}

bool Loader::event_info_is_valid() {
    std::queue<Data> temp_queue;
    while (!TotalData.empty()) {
        Data temp_data = TotalData.front();
        TotalData.pop();
        if (temp_data.upsilon_experiment != temp_data.Bsig_experiment || temp_data.upsilon_experiment != temp_data.Btag_experiment || temp_data.Bsig_experiment != temp_data.Btag_experiment) return false;
        if (temp_data.upsilon_run != temp_data.Bsig_run || temp_data.upsilon_run != temp_data.Btag_run || temp_data.Bsig_run != temp_data.Btag_run) return false;
        if (temp_data.upsilon_event != temp_data.Bsig_event || temp_data.upsilon_event != temp_data.Btag_event || temp_data.Bsig_event != temp_data.Btag_event) return false;
        if (temp_data.upsilon_candidate != temp_data.Bsig_candidate || temp_data.upsilon_candidate != temp_data.Btag_candidate || temp_data.Bsig_candidate != temp_data.Btag_candidate) return false;
        if (temp_data.upsilon_ncandidates != temp_data.Bsig_ncandidates || temp_data.upsilon_ncandidates != temp_data.Btag_ncandidates || temp_data.Bsig_ncandidates != temp_data.Btag_ncandidates) return false;
        temp_data.__experiment__ = temp_data.upsilon_experiment;
        temp_data.__run__ = temp_data.upsilon_run;
        temp_data.__event__ = temp_data.upsilon_event;
        temp_data.__candidate__ = temp_data.upsilon_candidate;
        temp_data.__ncandidates__ = temp_data.upsilon_ncandidates;
        temp_queue.push(temp_data);
    }
    TotalData.swap(temp_queue);
    return true;
}

void Loader::DrawTH1F(const char* name, const char* title, int nbins, double x_low, double x_high, Loader::Variable variable, int i, Loader::ValueOption dr) {
    if (TH1Fs.size() == current_TH1F) { // allocate new hist
        TH1F* hist = new TH1F(name, title, nbins, x_low, x_high);
        TH1Fs.push_back(hist);
    }
    else if (TH1Fs.size() > current_TH1F) { // use what I have
    }
    else { // error
        printf("ERROR! 000\n");
        exit(1);
    }

    TH1F* temp_hist = TH1Fs.at(current_TH1F);
    std::queue<Data> temp_queue;
    temp_queue.swap(TotalData);
    while (!temp_queue.empty()) {
        Data temp_data = temp_queue.front();
        temp_queue.pop();
        if (variable == Loader::Upsilon) {
            if(dr == Loader::Linear) temp_hist->Fill(temp_data.Upsilon_info[i]);
            else if (dr == Loader::Log) temp_hist->Fill(TMath::Log10(temp_data.Upsilon_info[i]));
            else { printf("ERROR! 001\n"); exit(1); }
        }
        else if(variable == Loader::Bsig) {
            if (dr == Loader::Linear) temp_hist->Fill(temp_data.Bsig_info[i]);
            else if (dr == Loader::Log) temp_hist->Fill(TMath::Log10(temp_data.Bsig_info[i]));
            else { printf("ERROR! 002\n"); exit(1); }
        }
        else if(variable == Loader::Btag) {
            if (dr == Loader::Linear) temp_hist->Fill(temp_data.Btag_info[i]);
            else if (dr == Loader::Log) temp_hist->Fill(TMath::Log10(temp_data.Btag_info[i]));
            else { printf("ERROR! 003\n"); exit(1); }
        }
        else {
            printf("ERROR! 004\n");
            exit(1);
        }
        TotalData.push(temp_data);
    }

    current_TH1F++;
}

void Loader::DrawTH1F(const char* name, const char* title, int nbins, double x_low, double x_high, Loader::Variable variable, int i, Loader::Qualifier qualifier, Loader::DecayMode decaymode, Loader::ValueOption dr) {
    if (TH1Fs.size() == current_TH1F) { // allocate new hist
        TH1F* hist = new TH1F(name, title, nbins, x_low, x_high);
        TH1Fs.push_back(hist);
    }
    else if (TH1Fs.size() > current_TH1F) { // use what I have
    }
    else { // error
        printf("ERROR! 005\n");
        exit(1);
    }

    TH1F* temp_hist = TH1Fs.at(current_TH1F);
    std::queue<Data> temp_queue;
    temp_queue.swap(TotalData);
    while (!temp_queue.empty()) {
        Data temp_data = temp_queue.front();
        temp_queue.pop();

        if (qualifier == Loader::when) {
            if (!TrueIfDecayModeMatch(temp_data, decaymode)) continue;
        }
        else if (qualifier == Loader::except) {
            if (TrueIfDecayModeMatch(temp_data, decaymode)) continue;
        }
        else {
            printf("ERROR! 006\n");
            exit(1);
        }

        if (variable == Loader::Upsilon) {
            if (dr == Loader::Linear) temp_hist->Fill(temp_data.Upsilon_info[i]);
            else if (dr == Loader::Log) temp_hist->Fill(TMath::Log10(temp_data.Upsilon_info[i]));
            else { printf("ERROR! 007\n"); exit(1); }
        }
        else if (variable == Loader::Bsig) {
            if (dr == Loader::Linear) temp_hist->Fill(temp_data.Bsig_info[i]);
            else if (dr == Loader::Log) temp_hist->Fill(TMath::Log10(temp_data.Bsig_info[i]));
            else { printf("ERROR! 008\n"); exit(1); }
        }
        else if (variable == Loader::Btag) {
            if (dr == Loader::Linear) temp_hist->Fill(temp_data.Btag_info[i]);
            else if (dr == Loader::Log) temp_hist->Fill(TMath::Log10(temp_data.Btag_info[i]));
            else { printf("ERROR! 009\n"); exit(1); }
        }
        else {
            printf("ERROR! 010\n");
            exit(1);
        }
        TotalData.push(temp_data);
    }

    current_TH1F++;
}

void Loader::DrawTH2F(const char* name, const char* title, int nbinsx, double xlow, double xup, int nbinsy, double ylow, double yup, Loader::Variable variable_1, int i, Loader::Variable variable_2, int j) {
    if (TH2Fs.size() == current_TH2F) { // allocate new hist
        TH2F* hist = new TH2F(name, title, nbinsx, xlow, xup, nbinsy, ylow, yup);
        TH2Fs.push_back(hist);
    }
    else if (TH2Fs.size() > current_TH2F) { // use what I have
    }
    else { // error
        printf("ERROR! 011\n");
        exit(1);
    }

    TH2F* temp_hist = TH2Fs.at(current_TH2F);
    std::queue<Data> temp_queue;
    temp_queue.swap(TotalData);
    while (!temp_queue.empty()) {
        Data temp_data = temp_queue.front();
        temp_queue.pop();
        if (variable_1 == Loader::Upsilon) {
            if (variable_2 == Loader::Upsilon) {
                temp_hist->Fill(temp_data.Upsilon_info[i], temp_data.Upsilon_info[j]);
            }
            else if (variable_2 == Loader::Bsig) {
                temp_hist->Fill(temp_data.Upsilon_info[i], temp_data.Bsig_info[j]);
            }
            else if (variable_2 == Loader::Btag) {
                temp_hist->Fill(temp_data.Upsilon_info[i], temp_data.Btag_info[j]);
            }
            else {
                printf("ERROR! 012\n");
                exit(1);
            }
        }
        else if (variable_1 == Loader::Bsig) {
            if (variable_2 == Loader::Upsilon) {
                temp_hist->Fill(temp_data.Bsig_info[i], temp_data.Upsilon_info[j]);
            }
            else if (variable_2 == Loader::Bsig) {
                temp_hist->Fill(temp_data.Bsig_info[i], temp_data.Bsig_info[j]);
            }
            else if (variable_2 == Loader::Btag) {
                temp_hist->Fill(temp_data.Bsig_info[i], temp_data.Btag_info[j]);
            }
            else {
                printf("ERROR! 013\n");
                exit(1);
            }
        }
        else if (variable_1 == Loader::Btag) {
            if (variable_2 == Loader::Upsilon) {
                temp_hist->Fill(temp_data.Btag_info[i], temp_data.Upsilon_info[j]);
            }
            else if (variable_2 == Loader::Bsig) {
                temp_hist->Fill(temp_data.Btag_info[i], temp_data.Bsig_info[j]);
            }
            else if (variable_2 == Loader::Btag) {
                temp_hist->Fill(temp_data.Btag_info[i], temp_data.Btag_info[j]);
            }
            else {
                printf("ERROR! 014\n");
                exit(1);
            }
        }
        else {
            printf("ERROR! 015\n");
            exit(1);
        }
        TotalData.push(temp_data);
    }

    current_TH2F++;
}

void Loader::DrawTH2F(const char* name, const char* title, int nbinsx, double xlow, double xup, int nbinsy, double ylow, double yup, Loader::Variable variable_1, int i, Loader::Variable variable_2, int j, Loader::Qualifier qualifier, Loader::DecayMode decaymode) {
    if (TH2Fs.size() == current_TH2F) { // allocate new hist
        TH2F* hist = new TH2F(name, title, nbinsx, xlow, xup, nbinsy, ylow, yup);
        TH2Fs.push_back(hist);
    }
    else if (TH2Fs.size() > current_TH2F) { // use what I have
    }
    else { // error
        printf("ERROR! 016\n");
        exit(1);
    }

    TH2F* temp_hist = TH2Fs.at(current_TH2F);
    std::queue<Data> temp_queue;
    temp_queue.swap(TotalData);
    while (!temp_queue.empty()) {
        Data temp_data = temp_queue.front();
        temp_queue.pop();

        if (qualifier == Loader::when) {
            if (!TrueIfDecayModeMatch(temp_data, decaymode)) continue;
        }
        else if (qualifier == Loader::except) {
            if (TrueIfDecayModeMatch(temp_data, decaymode)) continue;
        }
        else {
            printf("ERROR! 017\n");
            exit(1);
        }

        if (variable_1 == Loader::Upsilon) {
            if (variable_2 == Loader::Upsilon) {
                temp_hist->Fill(temp_data.Upsilon_info[i], temp_data.Upsilon_info[j]);
            }
            else if (variable_2 == Loader::Bsig) {
                temp_hist->Fill(temp_data.Upsilon_info[i], temp_data.Bsig_info[j]);
            }
            else if (variable_2 == Loader::Btag) {
                temp_hist->Fill(temp_data.Upsilon_info[i], temp_data.Btag_info[j]);
            }
            else {
                printf("ERROR! 018\n");
                exit(1);
            }
        }
        else if (variable_1 == Loader::Bsig) {
            if (variable_2 == Loader::Upsilon) {
                temp_hist->Fill(temp_data.Bsig_info[i], temp_data.Upsilon_info[j]);
            }
            else if (variable_2 == Loader::Bsig) {
                temp_hist->Fill(temp_data.Bsig_info[i], temp_data.Bsig_info[j]);
            }
            else if (variable_2 == Loader::Btag) {
                temp_hist->Fill(temp_data.Bsig_info[i], temp_data.Btag_info[j]);
            }
            else {
                printf("ERROR! 019\n");
                exit(1);
            }
        }
        else if (variable_1 == Loader::Btag) {
            if (variable_2 == Loader::Upsilon) {
                temp_hist->Fill(temp_data.Btag_info[i], temp_data.Upsilon_info[j]);
            }
            else if (variable_2 == Loader::Bsig) {
                temp_hist->Fill(temp_data.Btag_info[i], temp_data.Bsig_info[j]);
            }
            else if (variable_2 == Loader::Btag) {
                temp_hist->Fill(temp_data.Btag_info[i], temp_data.Btag_info[j]);
            }
            else {
                printf("ERROR! 020\n");
                exit(1);
            }
        }
        else {
            printf("ERROR! 021\n");
            exit(1);
        }
        TotalData.push(temp_data);
    }

    current_TH2F++;
}

void Loader::DrawTHStack(const char* name, const char* title, int nbins, double x_low, double x_high, Loader::Variable variable, int i, Loader::ValueOption dr) {
    if (THStacks.size() == current_THStack) { // allocate new thstacks
        THStack* stack = new THStack(name, title);
        THStacks.push_back(stack);
        for (int i = 0; i < Loader::MAX_NUM_DECAYMODE; i++) {
            TH1F* hist = new TH1F((std::string(name) + std::string("_") + std::to_string(i)).c_str(), title, nbins, x_low, x_high);
            TH1Fs_THStack[i].push_back(hist);
        }
    }
    else if (THStacks.size() > current_THStack) { // use what I have
    }
    else { // error
        printf("ERROR! 022\n");
        exit(1);
    }

    TH1F* temp_hist[Loader::MAX_NUM_DECAYMODE];
    for (int i = 0; i < MAX_NUM_DECAYMODE;i++) {
        temp_hist[i] = TH1Fs_THStack[i].at(current_THStack);
    }

    std::queue<Data> temp_queue;
    temp_queue.swap(TotalData);
    while (!temp_queue.empty()) {
        Data temp_data = temp_queue.front();
        temp_queue.pop();

        Loader::DecayMode decaymodeid = Loader::MAX_NUM_DECAYMODE;
        for (int i = 0; i < MAX_NUM_DECAYMODE; i++) {
            if (TrueIfDecayModeMatch(temp_data, static_cast<Loader::DecayMode>(i) )) {
                decaymodeid = static_cast<Loader::DecayMode>(i);
                break;
            }
        }
        if (decaymodeid == Loader::MAX_NUM_DECAYMODE) {
            printf("ERROR! 023\n");
            exit(1);
        }

        if (variable == Loader::Upsilon) {
            if (dr == Loader::Linear) temp_hist[decaymodeid]->Fill(temp_data.Upsilon_info[i]);
            else if (dr == Loader::Log) temp_hist[decaymodeid]->Fill(TMath::Log10(temp_data.Upsilon_info[i]));
            else { printf("ERROR! 024\n"); exit(1); }
        }
        else if (variable == Loader::Bsig) {
            if (dr == Loader::Linear) temp_hist[decaymodeid]->Fill(temp_data.Bsig_info[i]);
            else if (dr == Loader::Log) temp_hist[decaymodeid]->Fill(TMath::Log10(temp_data.Bsig_info[i]));
            else { printf("ERROR! 025\n"); exit(1); }
        }
        else if (variable == Loader::Btag) {
            if (dr == Loader::Linear) temp_hist[decaymodeid]->Fill(temp_data.Btag_info[i]);
            else if (dr == Loader::Log) temp_hist[decaymodeid]->Fill(TMath::Log10(temp_data.Btag_info[i]));
            else { printf("ERROR! 026\n"); exit(1); }
        }
        else {
            printf("ERROR! 027\n");
            exit(1);
        }
        TotalData.push(temp_data);
    }

    current_THStack++;
}

void Loader::PrintInformation(std::string title, std::string filename, const char* type, const char* MC_version, const char* category, bool smartmode) {
    typedef struct labels {
        int __experiment__;
        int __run__;
        unsigned int __event__;
        int __ncandidates__;
        bool IsThisModeExist[Loader::MAX_NUM_DECAYMODE];
    } Labels;
    std::vector<Labels> label_list;

    if (N_events.size() == current_N_event && N_candidates.size() == current_N_candidate && N_events.size() == N_candidates.size()) { // allocate new int
        N_events.push_back(0);
        N_candidates.push_back(0);
        titles.push_back(title);
        for (int i = 0; i < Loader::MAX_NUM_DECAYMODE; i++) N_candidates_modes[i].push_back(0);
        for (int i = 0; i < Loader::MAX_NUM_DECAYMODE_MC; i++) N_MC_modes[i].push_back(0);
    }
    else if (N_events.size() > current_N_event && N_candidates.size() > current_N_candidate && N_events.size() == N_candidates.size() && current_N_event == current_N_candidate) { // use what I have
    }
    else { // error
        printf("ERROR! 028\n");
        exit(1);
    }

    std::queue<Data> temp_queue;
    temp_queue.swap(TotalData);
    while (!temp_queue.empty()) {
        Data temp = temp_queue.front();
        temp_queue.pop();

        bool overlap = false;
        for (unsigned int i = 0; i < label_list.size(); i++) {
            if (label_list.at(i).__experiment__ == temp.__experiment__ && label_list.at(i).__run__ == temp.__run__ && label_list.at(i).__event__ == temp.__event__ && label_list.at(i).__ncandidates__ == temp.__ncandidates__) {
                Loader::DecayMode decaymodeid = Loader::MAX_NUM_DECAYMODE;
                for (int j = 0; j < Loader::MAX_NUM_DECAYMODE; j++) {
                    if (TrueIfDecayModeMatch(temp, static_cast<Loader::DecayMode>(j))) {
                        decaymodeid = static_cast<Loader::DecayMode>(j);
                        break;
                    }
                }
                if (decaymodeid == Loader::MAX_NUM_DECAYMODE) {
                    printf("ERROR! 029\n");
                    exit(1);
                }
                label_list.at(i).IsThisModeExist[decaymodeid] = true;

                overlap = true;
            }
        }
        if (overlap == false) {
            // Number of event
            if (smartmode == false) N_events.at(current_N_event) = N_events.at(current_N_event) + 1;
            else {
                if (strcmp(type, "SIGNAL") == 0) {
                    if (filename.find("B2Knunu") != std::string::npos) {
                        double correction_weight = corrector.GetCorrectionFactor(temp.Decay_syst_ff[index_q2] * temp.Decay_syst_ff[index_q2], "Bplus");
                        N_events.at(current_N_event) = N_events.at(current_N_event) + ObtainWeight(type, MC_version, category, filename) * correction_weight;
                    }
                    else if (filename.find("B2Kstarnunu") != std::string::npos) N_events.at(current_N_event) = N_events.at(current_N_event) + ObtainWeight(type, MC_version, category, filename);
                    else if (filename.find("B2Xsnunu") != std::string::npos) {
                        double correction_fragmentation = corrector_Fragmentation.GetCorrectionFactor(temp.Decay, temp.Decay_syst_ff[index_MXs_Bc], Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MC_version);
                        N_events.at(current_N_event) = N_events.at(current_N_event) + ObtainWeight(type, MC_version, category, filename) * correction_fragmentation;
                    }
                    else if (filename.find("B02K0nunu") != std::string::npos) {
                        double correction_weight = corrector.GetCorrectionFactor(temp.Decay_syst_ff[index_q2] * temp.Decay_syst_ff[index_q2], "Bzero");
                        N_events.at(current_N_event) = N_events.at(current_N_event) + ObtainWeight(type, MC_version, category, filename) * correction_weight;
                    }
                    else if (filename.find("B02Kstar0nunu") != std::string::npos) N_events.at(current_N_event) = N_events.at(current_N_event) + ObtainWeight(type, MC_version, category, filename);
                    else if (filename.find("B02Xsnunu") != std::string::npos) {
                        double correction_fragmentation = corrector_Fragmentation.GetCorrectionFactor(temp.Decay, temp.Decay_syst_ff[index_MXs_B0], Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MC_version);
                        N_events.at(current_N_event) = N_events.at(current_N_event) + ObtainWeight(type, MC_version, category, filename) * correction_fragmentation;
                    }
                    else { N_events.at(current_N_event) = N_events.at(current_N_event) + ObtainWeight(type, MC_version, category, filename); }
                }
                else N_events.at(current_N_event) = N_events.at(current_N_event) + ObtainWeight(type, MC_version, category, filename);
            }
            Labels temp_Labels;
            temp_Labels.__experiment__ = temp.__experiment__;
            temp_Labels.__run__ = temp.__run__;
            temp_Labels.__event__ = temp.__event__;
            temp_Labels.__ncandidates__ = temp.__ncandidates__;

            for (int i = 0; i < Loader::MAX_NUM_DECAYMODE; i++) temp_Labels.IsThisModeExist[i] = false;
            Loader::DecayMode decaymodeid = Loader::MAX_NUM_DECAYMODE;
            for (int i = 0; i < Loader::MAX_NUM_DECAYMODE; i++) {
                if (TrueIfDecayModeMatch(temp, static_cast<Loader::DecayMode>(i))) {
                    decaymodeid = static_cast<Loader::DecayMode>(i);
                    break;
                }
            }
            if (decaymodeid == Loader::MAX_NUM_DECAYMODE) {
                printf("ERROR! 095\n");
                exit(1);
            }
            temp_Labels.IsThisModeExist[decaymodeid] = true;

            label_list.push_back(temp_Labels);


            Loader::DecayModeMC decaymodeid_MC = Loader::MAX_NUM_DECAYMODE_MC;
            if (DoesItHaveXsBranch) {
                for (int i = 0; i < Loader::MAX_NUM_DECAYMODE_MC; i++) { // find MC decay mode
                    if (TrueIfDecayModeMatch_MC(temp, static_cast<Loader::DecayModeMC>(i))) {
                        decaymodeid_MC = static_cast<Loader::DecayModeMC>(i);
                        break;
                    }
                }
                if (decaymodeid_MC == Loader::MAX_NUM_DECAYMODE_MC) {
                    printf("ERROR! MC decay id cannot be found\n");
                    exit(1);
                }

                // Number of event with MC decayID (scaled)
                double temp_N = -1;
                if (decaymodeid_MC == Loader::Xsu2Kc_MC) {
                    double correction_weight = corrector.GetCorrectionFactor(temp.Decay_syst_ff[index_q2] * temp.Decay_syst_ff[index_q2], "Bplus");
                    temp_N = ObtainWeight(type, MC_version, category, "B2Knunu") * correction_weight;
                }
                else if (decaymodeid_MC == Loader::Xsu2Kcstar2KcPi0_MC || decaymodeid_MC == Loader::Xsu2Kcstar2K0Pic_MC) temp_N = ObtainWeight(type, MC_version, category, "B2Kstarnunu");
                else if (static_cast<int>(Xsu2KcPi0_MC) <= static_cast<int>(decaymodeid_MC) && static_cast<int>(decaymodeid_MC) <= static_cast<int>(Xsu2KcKcKcPi0_MC)) {
                    double correction_fragmentation = corrector_Fragmentation.GetCorrectionFactor(temp.Decay, temp.Decay_syst_ff[index_MXs_Bc], Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MC_version);
                    temp_N = ObtainWeight(type, MC_version, category, "B2Xsnunu") * correction_fragmentation;
                }
                else if (decaymodeid_MC == Loader::Xsd2K0_MC) {
                    double correction_weight = corrector.GetCorrectionFactor(temp.Decay_syst_ff[index_q2] * temp.Decay_syst_ff[index_q2], "Bzero");
                    temp_N = ObtainWeight(type, MC_version, category, "B02K0nunu") * correction_weight;
                }
                else if (decaymodeid_MC == Loader::Xsd2K0star2KcPic_MC || decaymodeid_MC == Loader::Xsd2K0star2K0Pi0_MC) temp_N = ObtainWeight(type, MC_version, category, "B02Kstar0nunu");
                else if (static_cast<int>(Xsd2KcPic_MC) <= static_cast<int>(decaymodeid_MC) && static_cast<int>(decaymodeid_MC) <= static_cast<int>(other)) {
                    double correction_fragmentation = corrector_Fragmentation.GetCorrectionFactor(temp.Decay, temp.Decay_syst_ff[index_MXs_B0], Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MC_version);
                    temp_N = ObtainWeight(type, MC_version, category, "B02Xsnunu") * correction_fragmentation;
                }
                else {
                    printf("ERROR 265\n");
                    exit(1);
                }
                N_MC_modes[decaymodeid_MC].at(current_N_event) = N_MC_modes[decaymodeid_MC].at(current_N_event) + temp_N;
            }
        }

        Loader::DecayMode decaymodeid = Loader::MAX_NUM_DECAYMODE;
        for (int i = 0; i < Loader::MAX_NUM_DECAYMODE; i++) {
            if (TrueIfDecayModeMatch(temp, static_cast<Loader::DecayMode>(i))) {
                decaymodeid = static_cast<Loader::DecayMode>(i);
                break;
            }
        }
        if (decaymodeid == Loader::MAX_NUM_DECAYMODE) {
            printf("ERROR! 029\n");
            exit(1);
        }

        // Number of candidate of decayID & Number of candidate
        if (smartmode == false) {
            N_candidates_modes[decaymodeid].at(current_N_candidate) = N_candidates_modes[decaymodeid].at(current_N_candidate) + 1.0;
            N_candidates.at(current_N_candidate) = N_candidates.at(current_N_candidate) + 1.0;
        }
        else {
            if (strcmp(type, "SIGNAL") == 0) {
                if (filename.find("B2Knunu") != std::string::npos) {
                    double correction_weight = corrector.GetCorrectionFactor(temp.Decay_syst_ff[index_q2] * temp.Decay_syst_ff[index_q2], "Bplus");
                    N_candidates_modes[decaymodeid].at(current_N_candidate) = N_candidates_modes[decaymodeid].at(current_N_candidate) + ObtainWeight(type, MC_version, category, filename) * correction_weight;
                    N_candidates.at(current_N_candidate) = N_candidates.at(current_N_candidate) + ObtainWeight(type, MC_version, category, filename) * correction_weight;
                }
                else if (filename.find("B2Kstarnunu") != std::string::npos) {
                    N_candidates_modes[decaymodeid].at(current_N_candidate) = N_candidates_modes[decaymodeid].at(current_N_candidate) + ObtainWeight(type, MC_version, category, filename);
                    N_candidates.at(current_N_candidate) = N_candidates.at(current_N_candidate) + ObtainWeight(type, MC_version, category, filename);
                }
                else if (filename.find("B2Xsnunu") != std::string::npos) {
                    double correction_fragmentation = corrector_Fragmentation.GetCorrectionFactor(temp.Decay, temp.Decay_syst_ff[index_MXs_Bc], Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MC_version);
                    N_candidates_modes[decaymodeid].at(current_N_candidate) = N_candidates_modes[decaymodeid].at(current_N_candidate) + ObtainWeight(type, MC_version, category, filename) * correction_fragmentation;
                    N_candidates.at(current_N_candidate) = N_candidates.at(current_N_candidate) + ObtainWeight(type, MC_version, category, filename) * correction_fragmentation;
                }
                else if (filename.find("B02K0nunu") != std::string::npos) {
                    double correction_weight = corrector.GetCorrectionFactor(temp.Decay_syst_ff[index_q2] * temp.Decay_syst_ff[index_q2], "Bzero");
                    N_candidates_modes[decaymodeid].at(current_N_candidate) = N_candidates_modes[decaymodeid].at(current_N_candidate) + ObtainWeight(type, MC_version, category, filename) * correction_weight;
                    N_candidates.at(current_N_candidate) = N_candidates.at(current_N_candidate) + ObtainWeight(type, MC_version, category, filename) * correction_weight;
                }
                else if (filename.find("B02Kstar0nunu") != std::string::npos) {
                    N_candidates_modes[decaymodeid].at(current_N_candidate) = N_candidates_modes[decaymodeid].at(current_N_candidate) + ObtainWeight(type, MC_version, category, filename);
                    N_candidates.at(current_N_candidate) = N_candidates.at(current_N_candidate) + ObtainWeight(type, MC_version, category, filename);
                }
                else if (filename.find("B02Xsnunu") != std::string::npos) {
                    double correction_fragmentation = corrector_Fragmentation.GetCorrectionFactor(temp.Decay, temp.Decay_syst_ff[index_MXs_B0], Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MC_version);
                    N_candidates_modes[decaymodeid].at(current_N_candidate) = N_candidates_modes[decaymodeid].at(current_N_candidate) + ObtainWeight(type, MC_version, category, filename) * correction_fragmentation;
                    N_candidates.at(current_N_candidate) = N_candidates.at(current_N_candidate) + ObtainWeight(type, MC_version, category, filename) * correction_fragmentation;
                }
                else {
                    N_candidates_modes[decaymodeid].at(current_N_candidate) = N_candidates_modes[decaymodeid].at(current_N_candidate) + ObtainWeight(type, MC_version, category, filename);
                    N_candidates.at(current_N_candidate) = N_candidates.at(current_N_candidate) + ObtainWeight(type, MC_version, category, filename);
                }
            }
            else {
                N_candidates_modes[decaymodeid].at(current_N_candidate) = N_candidates_modes[decaymodeid].at(current_N_candidate) + ObtainWeight(type, MC_version, category, filename);
                N_candidates.at(current_N_candidate) = N_candidates.at(current_N_candidate) + ObtainWeight(type, MC_version, category, filename);
            }

        }

        TotalData.push(temp);
    }

    current_N_event++;
    current_N_candidate++;
}

void Loader::Cut(Loader::Variable variable, int i, Loader::Inequality inq, double value) {
    std::queue<Data> temp_queue;
    while (!TotalData.empty()) {
        Data temp_data = TotalData.front();
        TotalData.pop();
        if (variable == Loader::Upsilon) {
            if(inq == Loader::larger_than && temp_data.Upsilon_info[i] > value) temp_queue.push(temp_data);
            else if (inq == Loader::smaller_than && temp_data.Upsilon_info[i] < value) temp_queue.push(temp_data);
        }
        else if (variable == Loader::Bsig) {
            if (inq == Loader::larger_than && temp_data.Bsig_info[i] > value) temp_queue.push(temp_data);
            else if (inq == Loader::smaller_than && temp_data.Bsig_info[i] < value) temp_queue.push(temp_data);
        }
        else if (variable == Loader::Btag) {
            if (inq == Loader::larger_than && temp_data.Btag_info[i] > value) temp_queue.push(temp_data);
            else if (inq == Loader::smaller_than && temp_data.Btag_info[i] < value) temp_queue.push(temp_data);
        }
        else {
            printf("ERROR! 030\n");
            exit(1);
        }
    }
    TotalData.swap(temp_queue);
}

void Loader::Cut(Loader::Variable variable_1, int i_1, Loader::Arithmetic ari, Loader::Variable variable_2, int i_2, Loader::Inequality inq, double value) {
    std::queue<Data> temp_queue;
    while (!TotalData.empty()) {
        Data temp_data = TotalData.front();
        TotalData.pop();

        double value_1 = 0;
        double value_2 = 0;
        double value_3 = 0;

        if (variable_1 == Loader::Upsilon) value_1 = (double) temp_data.Upsilon_info[i_1];
        else if (variable_1 == Loader::Bsig) value_1 = (double) temp_data.Bsig_info[i_1];
        else if (variable_1 == Loader::Btag) value_1 = (double) temp_data.Btag_info[i_1];
        else {
            printf("ERROR! 1740\n");
            exit(1);
        }

        if (variable_2 == Loader::Upsilon) value_2 = (double) temp_data.Upsilon_info[i_2];
        else if (variable_2 == Loader::Bsig) value_2 = (double) temp_data.Bsig_info[i_2];
        else if (variable_2 == Loader::Btag) value_2 = (double) temp_data.Btag_info[i_2];
        else {
            printf("ERROR! 1746\n");
            exit(1);
        }

        if (ari == Loader::plus) value_3 = value_1 + value_2;
        else if (ari == Loader::minus) value_3 = value_1 - value_2;
        else if (ari == Loader::times) value_3 = value_1 * value_2;
        else if (ari == Loader::dividedby) {
            if (value_2 == 0) {
                printf("you try to divide by zero!\n");
                exit(1);
            }
            value_3 = value_1 / value_2;
        }
        else {
            printf("ERROR! 1762\n");
            exit(1);
        }

        if (inq == Loader::larger_than && value_3 > value) temp_queue.push(temp_data);
        else if (inq == Loader::smaller_than && value_3 < value) temp_queue.push(temp_data);

    }
    TotalData.swap(temp_queue);
}

void Loader::Cut(const char* bracket_1, double const_1, Loader::Arithmetic ari_1, Loader::Variable variable_1, int i_1, const char* bracket_2, Loader::Arithmetic ari, const char* bracket_3, double const_2, Loader::Arithmetic ari_2, Loader::Variable variable_2, int i_2, const char* bracket_4, Loader::Inequality inq, double value) {
    if ((std::strcmp(bracket_1, "(") == 0) && (std::strcmp(bracket_2, ")") == 0) && (std::strcmp(bracket_3, "(") == 0) && (std::strcmp(bracket_4, ")") == 0)) {}
    else {
        printf("ERROR! 1777 Use proper brackets for Cut module!\n");
        exit(1);
    }

    std::queue<Data> temp_queue;
    while (!TotalData.empty()) {
        Data temp_data = TotalData.front();
        TotalData.pop();

        double value_1 = 0;
        double value_2 = 0;
        double value_3 = 0;

        if (variable_1 == Loader::Upsilon) value_1 = (double)temp_data.Upsilon_info[i_1];
        else if (variable_1 == Loader::Bsig) value_1 = (double)temp_data.Bsig_info[i_1];
        else if (variable_1 == Loader::Btag) value_1 = (double)temp_data.Btag_info[i_1];
        else {
            printf("ERROR! 1740\n");
            exit(1);
        }

        if (ari_1 == Loader::plus) value_1 = const_1 + value_1;
        else if (ari_1 == Loader::minus) value_1 = const_1 - value_1;
        else if (ari_1 == Loader::times) value_1 = const_1 * value_1;
        else if (ari_1 == Loader::dividedby) {
            if (value_1 == 0) {
                printf("you try to divide by zero!\n");
                exit(1);
            }
            value_1 = const_1 / value_1;
        }
        else {
            printf("ERROR! 1809\n");
            exit(1);
        }

        if (variable_2 == Loader::Upsilon) value_2 = (double)temp_data.Upsilon_info[i_2];
        else if (variable_2 == Loader::Bsig) value_2 = (double)temp_data.Bsig_info[i_2];
        else if (variable_2 == Loader::Btag) value_2 = (double)temp_data.Btag_info[i_2];
        else {
            printf("ERROR! 1746\n");
            exit(1);
        }

        if (ari_2 == Loader::plus) value_2 = const_2 + value_2;
        else if (ari_2 == Loader::minus) value_2 = const_2 - value_2;
        else if (ari_2 == Loader::times) value_2 = const_2 * value_2;
        else if (ari_2 == Loader::dividedby) {
            if (value_2 == 0) {
                printf("you try to divide by zero!\n");
                exit(1);
            }
            value_2 = const_2 / value_2;
        }
        else {
            printf("ERROR! 1809\n");
            exit(1);
        }

        if (ari == Loader::plus) value_3 = value_1 + value_2;
        else if (ari == Loader::minus) value_3 = value_1 - value_2;
        else if (ari == Loader::times) value_3 = value_1 * value_2;
        else if (ari == Loader::dividedby) {
            if (value_2 == 0) {
                printf("you try to divide by zero!\n");
                exit(1);
            }
            value_3 = value_1 / value_2;
        }
        else {
            printf("ERROR! 1762\n");
            exit(1);
        }

        if (inq == Loader::larger_than && value_3 > value) temp_queue.push(temp_data);
        else if (inq == Loader::smaller_than && value_3 < value) temp_queue.push(temp_data);

    }
    TotalData.swap(temp_queue);
}

void Loader::Cut(Loader::Variable variable, int i, Loader::Inequality inq, double value, Loader::Qualifier qualifier, Loader::DecayMode decaymode) {
    std::queue<Data> temp_queue;
    while (!TotalData.empty()) {
        Data temp_data = TotalData.front();
        TotalData.pop();

        if (qualifier == Loader::when) {
            if (!TrueIfDecayModeMatch(temp_data, decaymode)) {
                temp_queue.push(temp_data);
                continue;
            }
        }
        else if (qualifier == Loader::except) {
            if (TrueIfDecayModeMatch(temp_data, decaymode)) {
                temp_queue.push(temp_data);
                continue;
            }
        }
        else {
            printf("ERROR! 031\n");
            exit(1);
        }

        if (variable == Loader::Upsilon) {
            if (inq == Loader::larger_than && temp_data.Upsilon_info[i] > value) temp_queue.push(temp_data);
            else if (inq == Loader::smaller_than && temp_data.Upsilon_info[i] < value) temp_queue.push(temp_data);
        }
        else if (variable == Loader::Bsig) {
            if (inq == Loader::larger_than && temp_data.Bsig_info[i] > value) temp_queue.push(temp_data);
            else if (inq == Loader::smaller_than && temp_data.Bsig_info[i] < value) temp_queue.push(temp_data);
        }
        else if (variable == Loader::Btag) {
            if (inq == Loader::larger_than && temp_data.Btag_info[i] > value) temp_queue.push(temp_data);
            else if (inq == Loader::smaller_than && temp_data.Btag_info[i] < value) temp_queue.push(temp_data);
        }
        else {
            printf("ERROR! 032\n");
            exit(1);
        }
    }
    TotalData.swap(temp_queue);
}

void Loader::Cut(Loader::Variable variable, int i, Loader::Inequality inq, double value, Loader::Qualifier qualifier, Loader::Variable variable_qual, int i_qual, Loader::Inequality inq_qual, double value_qual) {
    std::queue<Data> temp_queue;
    while (!TotalData.empty()) {
        Data temp_data = TotalData.front();
        TotalData.pop();

        double qualifier_var = -1;
        bool satisfy = false;
        if (variable_qual == Loader::Upsilon) qualifier_var = temp_data.Upsilon_info[i_qual];
        else if(variable_qual == Loader::Bsig) qualifier_var = temp_data.Bsig_info[i_qual];
        else if(variable_qual == Loader::Btag) qualifier_var = temp_data.Btag_info[i_qual];
        else {
            printf("ERROR! 130\n");
            exit(1);
        }
        if (inq_qual == Loader::larger_than && qualifier_var > value_qual) satisfy = true;
        else if(inq_qual == Loader::smaller_than && qualifier_var < value_qual) satisfy = true;
        else satisfy = false;
        if (qualifier == Loader::when) {
            if (!satisfy) {
                temp_queue.push(temp_data);
                continue;
            }
        }
        else if (qualifier == Loader::except) {
            if (satisfy) {
                temp_queue.push(temp_data);
                continue;
            }
        }
        else {
            printf("ERROR! 131\n");
            exit(1);
        }

        if (variable == Loader::Upsilon) {
            if (inq == Loader::larger_than && temp_data.Upsilon_info[i] > value) temp_queue.push(temp_data);
            else if (inq == Loader::smaller_than && temp_data.Upsilon_info[i] < value) temp_queue.push(temp_data);
        }
        else if (variable == Loader::Bsig) {
            if (inq == Loader::larger_than && temp_data.Bsig_info[i] > value) temp_queue.push(temp_data);
            else if (inq == Loader::smaller_than && temp_data.Bsig_info[i] < value) temp_queue.push(temp_data);
        }
        else if (variable == Loader::Btag) {
            if (inq == Loader::larger_than && temp_data.Btag_info[i] > value) temp_queue.push(temp_data);
            else if (inq == Loader::smaller_than && temp_data.Btag_info[i] < value) temp_queue.push(temp_data);
        }
        else {
            printf("ERROR! 330\n");
            exit(1);
        }
    }
    TotalData.swap(temp_queue);
}

void Loader::PrintDebugLogIf(Loader::Variable variable, int i, Loader::Inequality inq, double value) {
    if (current_N_experiment_index > 0) { // allocate new int
        printf("The number of PrintDebugLogIf should not be larger than 1\n");
        printf("Only first PrintDebugLogIf is accepted\n");
        return;
    }

    std::queue<Data> temp_queue;
    temp_queue.swap(TotalData);

    while (!temp_queue.empty()) {
        Data temp_data = temp_queue.front();
        temp_queue.pop();
        if (variable == Loader::Upsilon) {
            if (inq == Loader::larger_than && temp_data.Upsilon_info[i] > value) {
                experiment_indices.push_back(temp_data.__experiment__);
                run_indices.push_back(temp_data.__run__);
                event_indices.push_back(temp_data.__event__);
                candidate_indices.push_back(temp_data.__candidate__);
                ncandidates_indices.push_back(temp_data.__ncandidates__);
            }
            else if (inq == Loader::smaller_than && temp_data.Upsilon_info[i] < value) {
                experiment_indices.push_back(temp_data.__experiment__);
                run_indices.push_back(temp_data.__run__);
                event_indices.push_back(temp_data.__event__);
                candidate_indices.push_back(temp_data.__candidate__);
                ncandidates_indices.push_back(temp_data.__ncandidates__);
            }
        }
        else if (variable == Loader::Bsig) {
            if (inq == Loader::larger_than && temp_data.Bsig_info[i] > value) {
                experiment_indices.push_back(temp_data.__experiment__);
                run_indices.push_back(temp_data.__run__);
                event_indices.push_back(temp_data.__event__);
                candidate_indices.push_back(temp_data.__candidate__);
                ncandidates_indices.push_back(temp_data.__ncandidates__);
            }
            else if (inq == Loader::smaller_than && temp_data.Bsig_info[i] < value) {
                experiment_indices.push_back(temp_data.__experiment__);
                run_indices.push_back(temp_data.__run__);
                event_indices.push_back(temp_data.__event__);
                candidate_indices.push_back(temp_data.__candidate__);
                ncandidates_indices.push_back(temp_data.__ncandidates__);
            }
        }
        else if (variable == Loader::Btag) {
            if (inq == Loader::larger_than && temp_data.Btag_info[i] > value) {
                experiment_indices.push_back(temp_data.__experiment__);
                run_indices.push_back(temp_data.__run__);
                event_indices.push_back(temp_data.__event__);
                candidate_indices.push_back(temp_data.__candidate__);
                ncandidates_indices.push_back(temp_data.__ncandidates__);
            }
            else if (inq == Loader::smaller_than && temp_data.Btag_info[i] < value) {
                experiment_indices.push_back(temp_data.__experiment__);
                run_indices.push_back(temp_data.__run__);
                event_indices.push_back(temp_data.__event__);
                candidate_indices.push_back(temp_data.__candidate__);
                ncandidates_indices.push_back(temp_data.__ncandidates__);
            }
        }
        else {
            printf("ERROR! 033\n");
            exit(1);
        }
        TotalData.push(temp_data);
    }

    DebugIsOn = true;
    current_N_experiment_index++;
}

void Loader::BCS(Loader::Variable variable, int index, Loader::BCS_criterion crit) {
    std::queue<Data> new_container;

    while (!TotalData.empty()) {
        std::vector<Data> temp;
        Data initial_data = TotalData.front();
        int experiment_ = initial_data.__experiment__;
        int run_ = initial_data.__run__;
        unsigned int event_ = initial_data.__event__;
        int ncandidates_ = initial_data.__ncandidates__;
        while (true) { // I suppose that the order of data exists
            if (TotalData.empty()) break;
            Data temp_data = TotalData.front();
            if (temp_data.__experiment__ == experiment_ && temp_data.__run__ == run_ && temp_data.__event__ == event_ && temp_data.__ncandidates__ == ncandidates_) {
                TotalData.pop();
                temp.push_back(temp_data);
            }
            else break;
        }

        switch (variable) {
        case Loader::Upsilon:
            if (crit == Loader::Highest) {
                double max = -std::numeric_limits<double>::max();
                int best_candidate_index = -1;
                for (unsigned int i = 0; i < temp.size(); i++) {
                    if (temp.at(i).Upsilon_info[index] > max) { max = temp.at(i).Upsilon_info[index]; best_candidate_index = i; }
                }
                if (best_candidate_index == -1) { printf("error!\n"); exit(1); }
                new_container.push(temp[best_candidate_index]);
            }
            else if (crit == Loader::Lowest) {
                double min = std::numeric_limits<double>::max();
                int best_candidate_index = -1;
                for (unsigned int i = 0; i < temp.size(); i++) {
                    if (temp.at(i).Upsilon_info[index] < min) { min = temp.at(i).Upsilon_info[index]; best_candidate_index = i; }
                }
                if (best_candidate_index == -1) { printf("error!\n"); exit(1); }
                new_container.push(temp[best_candidate_index]);
            }
            break;
        case Loader::Bsig:
            if (crit == Loader::Highest) {
                double max = -std::numeric_limits<double>::max();
                int best_candidate_index = -1;
                for (unsigned int i = 0; i < temp.size(); i++) {
                    if (temp.at(i).Bsig_info[index] > max) { max = temp.at(i).Bsig_info[index]; best_candidate_index = i; }
                }
                if (best_candidate_index == -1) { printf("error!\n"); exit(1); }
                new_container.push(temp[best_candidate_index]);
            }
            else if (crit == Loader::Lowest) {
                double min = std::numeric_limits<double>::max();
                int best_candidate_index = -1;
                for (unsigned int i = 0; i < temp.size(); i++) {
                    if (temp.at(i).Bsig_info[index] < min) { min = temp.at(i).Bsig_info[index]; best_candidate_index = i; }
                }
                if (best_candidate_index == -1) { printf("error!\n"); exit(1); }
                new_container.push(temp[best_candidate_index]);
            }
            break;
        case Loader::Btag:
            if (crit == Loader::Highest) {
                double max = -std::numeric_limits<double>::max();
                int best_candidate_index = -1;
                for (unsigned int i = 0; i < temp.size(); i++) {
                    if (temp.at(i).Btag_info[index] > max) { max = temp.at(i).Btag_info[index]; best_candidate_index = i; }
                }
                if (best_candidate_index == -1) { printf("error!\n"); exit(1); }
                new_container.push(temp[best_candidate_index]);
            }
            else if (crit == Loader::Lowest) {
                double min = std::numeric_limits<double>::max();
                int best_candidate_index = -1;
                for (unsigned int i = 0; i < temp.size(); i++) {
                    if (temp.at(i).Btag_info[index] < min) { min = temp.at(i).Btag_info[index]; best_candidate_index = i; }
                }
                if (best_candidate_index == -1) { printf("error!\n"); exit(1); }
                new_container.push(temp[best_candidate_index]);
            }
            break;
        default:
            printf("ERROR! 034\n");
            exit(1);
            break;
        }
    }
    while (!TotalData.empty()) TotalData.pop();
    while (!new_container.empty()) {
        Data temp_data = new_container.front();
        new_container.pop();
        TotalData.push(temp_data);
    }
}

bool Loader::IsBCSValid() { // modified for makeshift!!
    bool IsItValid = true;
    bool Eeclsame = false;

    typedef struct labels {
        int __experiment__;
        int __run__;
        unsigned int __event__;
        int __ncandidates__;
        double Eecl; // Makeshift
    } Labels;

    std::vector<Labels> label_list;
    std::queue<Data> TotalData_;
    TotalData_.swap(TotalData);

    while (!TotalData_.empty()) {
        Data temp = TotalData_.front();
        TotalData_.pop();
        for (unsigned int i = 0; i < label_list.size(); i++) {
            if (label_list.at(i).__experiment__ == temp.__experiment__ && label_list.at(i).__run__ == temp.__run__ && label_list.at(i).__event__ == temp.__event__ && label_list.at(i).__ncandidates__ == temp.__ncandidates__) {
                IsItValid = false;
                if (label_list.at(i).Eecl == temp.Upsilon_info[3]) Eeclsame = true; // makeshift
            }
        }
        Labels temp_Labels;
        temp_Labels.__experiment__ = temp.__experiment__;
        temp_Labels.__run__ = temp.__run__;
        temp_Labels.__event__ = temp.__event__;
        temp_Labels.__ncandidates__ = temp.__ncandidates__;
        temp_Labels.Eecl = temp.Upsilon_info[3];
        label_list.push_back(temp_Labels);

        if (IsItValid == true) TotalData.push(temp);
        else if (IsItValid == false && Eeclsame == true) {
            IsItValid = true; // makeshift
        }
        else if (IsItValid == false && Eeclsame == false) {
            IsItValid = true; // makeshift
            TotalData.push(temp);
        }
    }
    return IsItValid;
}

void Loader::End() {
    for (int i = 0; i < N_events.size();i++) {
        printf("%s\n", titles.at(i).c_str());
        printf("Number of event: %lf\n", N_events.at(i));
        printf("Number of candidate: %lf\n", N_candidates.at(i));
        for (int j = 0; j < Loader::MAX_NUM_DECAYMODE; j++) printf("Number of candidate of decayID %d: %lf\n", j, N_candidates_modes[j].at(i));
        if (AllOfThemHaveXsBranch) for (int j = 0; j < Loader::MAX_NUM_DECAYMODE_MC; j++) printf("Number of event with MC decayID %d(scaled): %lf\n", j, N_MC_modes[j].at(i));
    }


    if (Confusion_matrixIsOn == true) {
        printf("--------------- confusion matrix ---------------\n");
        for (int i = 0; i < Loader::MAX_NUM_DECAYMODE; i++) {
            for (int j = 0; j < Loader::MAX_NUM_DECAYMODE_MC; j++) {
                printf("%f ", Confusion[i][j]);
            }
            printf("\n");
        }
        printf("--------------- confusion matrix ---------------\n");

        double Sum_Each_Reco[Loader::MAX_NUM_DECAYMODE];
        for (int i = 0; i < Loader::MAX_NUM_DECAYMODE; i++) Sum_Each_Reco[i] = 0;

        for (int i = 0; i < Loader::MAX_NUM_DECAYMODE; i++) {
            for (int j = 0; j < Loader::MAX_NUM_DECAYMODE_MC; j++) {
                Sum_Each_Reco[i] = Sum_Each_Reco[i] + Confusion[i][j];
            }
        }

        for (int i = 0; i < Loader::MAX_NUM_DECAYMODE; i++) { // normalization
            for (int j = 0; j < Loader::MAX_NUM_DECAYMODE_MC; j++) {
                if (Sum_Each_Reco[i] != 0) Confusion[i][j] = Confusion[i][j] / Sum_Each_Reco[i];
            }
        }

        printf("--------------- normalized confusion matrix ---------------\n");
        for (int i = 0; i < Loader::MAX_NUM_DECAYMODE; i++) {
            for (int j = 0; j < Loader::MAX_NUM_DECAYMODE_MC; j++) {
                printf("%f ", Confusion[i][j]);
            }
            printf("\n");
        }
        printf("--------------- normalized confusion matrix ---------------\n");

        printf("--------------- square confusion matrix ---------------\n");
        for (int i = 0; i < Loader::MAX_NUM_DECAYMODE; i++) {
            for (int j = 0; j < Loader::MAX_NUM_DECAYMODE + 1; j++) {
                printf("%f ", Confusion_square[i][j]);
            }
            printf("\n");
        }
        printf("--------------- square confusion matrix ---------------\n");

        for (int i = 0; i < Loader::MAX_NUM_DECAYMODE; i++) Sum_Each_Reco[i] = 0;

        for (int i = 0; i < Loader::MAX_NUM_DECAYMODE; i++) {
            for (int j = 0; j < Loader::MAX_NUM_DECAYMODE + 1; j++) {
                Sum_Each_Reco[i] = Sum_Each_Reco[i] + Confusion_square[i][j];
            }
        }

        for (int i = 0; i < Loader::MAX_NUM_DECAYMODE; i++) { // normalization
            for (int j = 0; j < Loader::MAX_NUM_DECAYMODE + 1; j++) {
                if (Sum_Each_Reco[i] != 0) Confusion_square[i][j] = Confusion_square[i][j] / Sum_Each_Reco[i];
            }
        }

        printf("--------------- normalized square confusion matrix ---------------\n");
        for (int i = 0; i < Loader::MAX_NUM_DECAYMODE; i++) {
            for (int j = 0; j < Loader::MAX_NUM_DECAYMODE + 1; j++) {
                printf("%f ", Confusion_square[i][j]);
            }
            printf("\n");
        }
        printf("--------------- normalized square confusion matrix ---------------\n");

    }

    if (FOMIsOn == true) {
        printf("--------------- number of event to get FOM ---------------\n");
        printf("--------------- Oqq -> ---------------\n");
        printf("\n");
        for (int i = 0; i < Nstep; i++) {
            for (int j = 0; j < Nstep; j++) {
                printf("%lf ", FOM_Matrix[i][j]);
            }
            printf("\n");
        }
        printf("--------------- number of event to get FOM ---------------\n");
    }

    if (MCcountOn == true) {
        printf("--------------- number of Decay mode ---------------\n");
        for (int j = 0; j < Loader::MAX_NUM_DECAYMODE_MC; j++) printf("Number of event of MCdecayID %d: %lf\n", j, MCcount[j]);
        printf("--------------- number of Decay mode ---------------\n");
    }

    for (int i = 0; i < TH1Fs.size();i++) {
        TCanvas* c_temp = new TCanvas("c", "", 1500, 1200); c_temp->cd();
        TH1Fs.at(i)->Draw("Hist"); c_temp->SaveAs((std::string(TH1Fs.at(i)->GetName()) + ".png").c_str());
        delete c_temp;
    }

    for (int i = 0; i < TH2Fs.size();i++) {
        TH1D* projX;
        TH1D* projY;

        projX = TH2Fs.at(i)->ProjectionX();
        projY = TH2Fs.at(i)->ProjectionY();

        auto c_temp = new TCanvas("c1", "c1", 1000, 1000); c_temp->cd();
        gStyle->SetOptStat(0);

        TPad* center_pad = new TPad("center_pad", "center_pad", 0.0, 0.0, 0.6, 0.6);
        center_pad->Draw();

        TPad* right_pad = new TPad("right_pad", "right_pad", 0.55, 0.0, 1.0, 0.6);
        right_pad->Draw();

        TPad* top_pad = new TPad("top_pad", "top_pad", 0.0, 0.55, 0.6, 1.0);
        top_pad->Draw();

        center_pad->cd();
        gStyle->SetPalette(1);
        TH2Fs.at(i)->Draw("Hist");

        top_pad->cd();
        projX->SetFillColor(33);
        projX->Draw("bar");

        right_pad->cd();
        projY->SetFillColor(33);
        projY->Draw("hbar");
        gPad->RedrawAxis();

        c_temp->SaveAs((std::string(TH2Fs.at(i)->GetName()) + ".png").c_str());

        delete projX; delete projY; delete c_temp;
    }

    for (int i = 0; i < files.size(); i++) {
        TFile* temp_file = files.at(i);
        TTree* temp_tree_upsilon = trees_upsilon.at(i);
        TTree* temp_tree_Bsig = trees_Bsig.at(i);
        TTree* temp_tree_Btag = trees_Btag.at(i);
        TTree* temp_tree_Xs = trees_Xs.at(i);
        temp_file->cd();
        temp_tree_upsilon->Write();
        temp_tree_Bsig->Write();
        temp_tree_Btag->Write();
        if (AllOfThemHaveXsBranch) temp_tree_Xs->Write();
        else delete temp_tree_Xs;
        temp_file->Close();
    }

    if (DebugIsOn == true) {
        printf("========== Debug Log start ==========");
        printf("---------------------\n");
        for (int i = 0; i < experiment_indices.size(); i++) {
            printf("experiment num: %d\n", experiment_indices.at(i));
            printf("run num: %d\n", run_indices.at(i));
            printf("event num: %u\n", event_indices.at(i));
            printf("candidate num: %d\n", candidate_indices.at(i));
            printf("number of candidate: %d\n", ncandidates_indices.at(i));
            printf("---------------------\n");
        }
        printf("========== Debug Log end ==========");
    }

    for (int i = 0; i < THStacks.size(); i++) {
        TCanvas* c_temp = new TCanvas("c", "", 1500, 1200); c_temp->cd();
        gStyle->SetPalette(kOcean);

        for (int j = 0; j < Loader::MAX_NUM_DECAYMODE; j++) THStacks.at(i)->Add(TH1Fs_THStack[j].at(i));
        THStacks.at(i)->Draw("pfc"); c_temp->SaveAs((std::string(THStacks.at(i)->GetName()) + ".png").c_str());
        delete c_temp;
    }
}

void Loader::PrintRootFile(std::string output_name) {
    if (files.size() == current_file && trees_upsilon.size() == current_file && trees_Bsig.size() == current_file && trees_Btag.size() == current_file && trees_Xs.size() == current_file) { // allocate new TFile and TTree
        TFile* file = new TFile(output_name.c_str(), "recreate");
        file->cd();
        TTree* tree_upsilon = new TTree("Upsilon", "");
        TTree* tree_Bsig = new TTree("Bsig", "");
        TTree* tree_Btag = new TTree("Btag", "");
        TTree* tree_Xs = new TTree("Xs", "");

        /*================================================================*/
        // get event_info
        tree_upsilon->Branch("__experiment__", &UpsilonExperimentToTree);
        tree_upsilon->Branch("__run__", &UpsilonRunToTree);
        tree_upsilon->Branch("__event__", &UpsilonEventToTree);
        tree_upsilon->Branch("__candidate__", &UpsilonCandidateToTree);
        tree_upsilon->Branch("__ncandidates__", &UpsilonNcandidatesToTree);
        tree_Bsig->Branch("__experiment__", &BsigExperimentToTree);
        tree_Bsig->Branch("__run__", &BsigRunToTree);
        tree_Bsig->Branch("__event__", &BsigEventToTree);
        tree_Bsig->Branch("__candidate__", &BsigCandidateToTree);
        tree_Bsig->Branch("__ncandidates__", &BsigNcandidatesToTree);
        tree_Btag->Branch("__experiment__", &BtagExperimentToTree);
        tree_Btag->Branch("__run__", &BtagRunToTree);
        tree_Btag->Branch("__event__", &BtagEventToTree);
        tree_Btag->Branch("__candidate__", &BtagCandidateToTree);
        tree_Btag->Branch("__ncandidates__", &BtagNcandidatesToTree);

        // get decaymodeID
        tree_upsilon->Branch("extraInfo__bodecayModeID__bc", &Upsilon_decayIDToTree);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_decayIDToTree);

        // get Upsilon_info
        tree_upsilon->Branch("nROE_ECLClusters__bocleanMask__bc", &UpsilonDataToTree[0]);
        tree_upsilon->Branch("nROE_KLMClusters", &UpsilonDataToTree[1]);
        tree_upsilon->Branch("nROE_Tracks__bocleanMask__bc", &UpsilonDataToTree[2]);
        tree_upsilon->Branch("roeEextra__bocleanMask__bc", &UpsilonDataToTree[3]);
        tree_upsilon->Branch("nROE_NeutralECLClusters__bocleanMask__bc", &UpsilonDataToTree[4]);
        tree_upsilon->Branch("nROE_ParticlesInList__boK_S0__clmyKaonshort__bc", &UpsilonDataToTree[5]);
        tree_upsilon->Branch("nROE_ParticlesInList__bopi0__clmyneutralPion__bc", &UpsilonDataToTree[6]);
        tree_upsilon->Branch("missingMomentumOfEvent_theta", &UpsilonDataToTree[7]);
        tree_upsilon->Branch("missingMomentumOfEvent", &UpsilonDataToTree[8]);
        tree_upsilon->Branch("missingEnergyOfEventCMS", &UpsilonDataToTree[9]);
        tree_upsilon->Branch("nRemainingTracksInEvent", &UpsilonDataToTree[10]);
        tree_upsilon->Branch("roeNeextra__bocleanMask__bc", &UpsilonDataToTree[11]);
        tree_upsilon->Branch("useCMSFrame__boroeNeextra__bocleanMask__bc__bc", &UpsilonDataToTree[12]);
        tree_upsilon->Branch("nROE_ParticlesInList__bogamma__clmygamma__bc", &UpsilonDataToTree[13]);
        tree_upsilon->Branch("foxWolframR1", &UpsilonDataToTree[14]);
        tree_upsilon->Branch("foxWolframR2", &UpsilonDataToTree[15]);
        tree_upsilon->Branch("foxWolframR3", &UpsilonDataToTree[16]);
        tree_upsilon->Branch("foxWolframR4", &UpsilonDataToTree[17]);
        tree_upsilon->Branch("harmonicMomentThrust0", &UpsilonDataToTree[18]);
        tree_upsilon->Branch("harmonicMomentThrust1", &UpsilonDataToTree[19]);
        tree_upsilon->Branch("harmonicMomentThrust2", &UpsilonDataToTree[20]);
        tree_upsilon->Branch("harmonicMomentThrust3", &UpsilonDataToTree[21]);
        tree_upsilon->Branch("harmonicMomentThrust4", &UpsilonDataToTree[22]);
        tree_upsilon->Branch("cleoConeThrust0", &UpsilonDataToTree[23]);
        tree_upsilon->Branch("cleoConeThrust1", &UpsilonDataToTree[24]);
        tree_upsilon->Branch("cleoConeThrust2", &UpsilonDataToTree[25]);
        tree_upsilon->Branch("cleoConeThrust3", &UpsilonDataToTree[26]);
        tree_upsilon->Branch("cleoConeThrust4", &UpsilonDataToTree[27]);
        tree_upsilon->Branch("cleoConeThrust5", &UpsilonDataToTree[28]);
        tree_upsilon->Branch("cleoConeThrust6", &UpsilonDataToTree[29]);
        tree_upsilon->Branch("cleoConeThrust7", &UpsilonDataToTree[30]);
        tree_upsilon->Branch("cleoConeThrust8", &UpsilonDataToTree[31]);
        tree_upsilon->Branch("sphericity", &UpsilonDataToTree[32]);
        tree_upsilon->Branch("aplanarity", &UpsilonDataToTree[33]);
        tree_upsilon->Branch("thrust", &UpsilonDataToTree[34]);
        tree_upsilon->Branch("thrustAxisCosTheta", &UpsilonDataToTree[35]);
        tree_upsilon->Branch("MsquaredBsig_op0", &UpsilonDataToTree[36]);
        tree_upsilon->Branch("MsquaredBsig_op1", &UpsilonDataToTree[37]);
        tree_upsilon->Branch("MsquaredBsig_op2", &UpsilonDataToTree[38]);
        tree_upsilon->Branch("MsquaredBsig_op3", &UpsilonDataToTree[39]);
        tree_upsilon->Branch("MsquaredBsig_op4", &UpsilonDataToTree[40]);
        tree_upsilon->Branch("MsquaredBsig_op7", &UpsilonDataToTree[41]);
        tree_upsilon->Branch("roeP__bocleanMask__bc", &UpsilonDataToTree[42]);
        tree_upsilon->Branch("roeM__bocleanMask__bc", &UpsilonDataToTree[43]);
        tree_upsilon->Branch("roePTheta__bocleanMask__bc", &UpsilonDataToTree[44]);
        tree_upsilon->Branch("qsquared", &UpsilonDataToTree[45]);
        tree_upsilon->Branch("chiProb", &UpsilonDataToTree[46]);
        tree_upsilon->Branch("dr", &UpsilonDataToTree[47]);
        tree_upsilon->Branch("dz", &UpsilonDataToTree[48]);
        tree_upsilon->Branch("nParticlesInList__boe__pl__clElectronFBDT__bc", &UpsilonDataToTree[49]);
        tree_upsilon->Branch("nParticlesInList__bomu__pl__clMuonFBDT__bc", &UpsilonDataToTree[50]);
        tree_upsilon->Branch("nParticlesInList__boe__pl__clElectronFBDT_loose__bc", &UpsilonDataToTree[51]);
        tree_upsilon->Branch("nParticlesInList__bomu__pl__clMuonFBDT_loose__bc", &UpsilonDataToTree[52]);
        tree_upsilon->Branch("nParticlesInList__boe__pl__clElectronFBDT_tight__bc", &UpsilonDataToTree[53]);
        tree_upsilon->Branch("nParticlesInList__bomu__pl__clMuonFBDT_tight__bc", &UpsilonDataToTree[54]);
        tree_upsilon->Branch("beamE", &UpsilonDataToTree[55]);
        tree_upsilon->Branch("nROE_Tracks__bolooseMask__bc", &UpsilonDataToTree[56]);
        tree_upsilon->Branch("Ecms", &UpsilonDataToTree[57]);
        tree_upsilon->Branch("extraInfo__boNgamma__bc", &UpsilonDataToTree[58]);
        tree_upsilon->Branch("extraInfo__boEecl__bc", &UpsilonDataToTree[59]);
        tree_upsilon->Branch("extraInfo__boEecl_matched__bc", &UpsilonDataToTree[60]);
        tree_upsilon->Branch("extraInfo__boEecl_unmatched__bc", &UpsilonDataToTree[61]);
        tree_upsilon->Branch("extraInfo__boNgammav133__bc", &UpsilonDataToTree[62]);
        tree_upsilon->Branch("extraInfo__boEeclv133__bc", &UpsilonDataToTree[63]);
        tree_upsilon->Branch("extraInfo__boEeclv133_matched__bc", &UpsilonDataToTree[64]);
        tree_upsilon->Branch("extraInfo__boEeclv133_unmatched__bc", &UpsilonDataToTree[65]);
        tree_upsilon->Branch("extraInfo__boNgammav200__bc", &UpsilonDataToTree[66]);
        tree_upsilon->Branch("extraInfo__boNgammav200_matched__bc", &UpsilonDataToTree[67]);
        tree_upsilon->Branch("extraInfo__boNgammav200_unmatched__bc", &UpsilonDataToTree[68]);
        tree_upsilon->Branch("extraInfo__boEeclv200__bc", &UpsilonDataToTree[69]);
        tree_upsilon->Branch("extraInfo__boEeclv200_matched__bc", &UpsilonDataToTree[70]);
        tree_upsilon->Branch("extraInfo__boEeclv200_unmatched__bc", &UpsilonDataToTree[71]);
        tree_upsilon->Branch("extraInfo__boNgammav200_800__bc", &UpsilonDataToTree[72]);
        tree_upsilon->Branch("extraInfo__boNgammav200_800_matched__bc", &UpsilonDataToTree[73]);
        tree_upsilon->Branch("extraInfo__boNgammav200_800_unmatched__bc", &UpsilonDataToTree[74]);
        tree_upsilon->Branch("extraInfo__boEeclv200_800__bc", &UpsilonDataToTree[75]);
        tree_upsilon->Branch("extraInfo__boEeclv200_800_matched__bc", &UpsilonDataToTree[76]);
        tree_upsilon->Branch("extraInfo__boEeclv200_800_unmatched__bc", &UpsilonDataToTree[77]);
        tree_upsilon->Branch("extraInfo__boNgammav200_825__bc", &UpsilonDataToTree[78]);
        tree_upsilon->Branch("extraInfo__boNgammav200_825_matched__bc", &UpsilonDataToTree[79]);
        tree_upsilon->Branch("extraInfo__boNgammav200_825_unmatched__bc", &UpsilonDataToTree[80]);
        tree_upsilon->Branch("extraInfo__boEeclv200_825__bc", &UpsilonDataToTree[81]);
        tree_upsilon->Branch("extraInfo__boEeclv200_825_matched__bc", &UpsilonDataToTree[82]);
        tree_upsilon->Branch("extraInfo__boEeclv200_825_unmatched__bc", &UpsilonDataToTree[83]);
        tree_upsilon->Branch("extraInfo__boNgammav200_850__bc", &UpsilonDataToTree[84]);
        tree_upsilon->Branch("extraInfo__boNgammav200_850_matched__bc", &UpsilonDataToTree[85]);
        tree_upsilon->Branch("extraInfo__boNgammav200_850_unmatched__bc", &UpsilonDataToTree[86]);
        tree_upsilon->Branch("extraInfo__boEeclv200_850__bc", &UpsilonDataToTree[87]);
        tree_upsilon->Branch("extraInfo__boEeclv200_850_matched__bc", &UpsilonDataToTree[88]);
        tree_upsilon->Branch("extraInfo__boEeclv200_850_unmatched__bc", &UpsilonDataToTree[89]);
        tree_upsilon->Branch("extraInfo__boNgammav200_875__bc", &UpsilonDataToTree[90]);
        tree_upsilon->Branch("extraInfo__boNgammav200_875_matched__bc", &UpsilonDataToTree[91]);
        tree_upsilon->Branch("extraInfo__boNgammav200_875_unmatched__bc", &UpsilonDataToTree[92]);
        tree_upsilon->Branch("extraInfo__boEeclv200_875__bc", &UpsilonDataToTree[93]);
        tree_upsilon->Branch("extraInfo__boEeclv200_875_matched__bc", &UpsilonDataToTree[94]);
        tree_upsilon->Branch("extraInfo__boEeclv200_875_unmatched__bc", &UpsilonDataToTree[95]);
        tree_upsilon->Branch("extraInfo__boNgammav200_900__bc", &UpsilonDataToTree[96]);
        tree_upsilon->Branch("extraInfo__boNgammav200_900_matched__bc", &UpsilonDataToTree[97]);
        tree_upsilon->Branch("extraInfo__boNgammav200_900_unmatched__bc", &UpsilonDataToTree[98]);
        tree_upsilon->Branch("extraInfo__boEeclv200_900__bc", &UpsilonDataToTree[99]);
        tree_upsilon->Branch("extraInfo__boEeclv200_900_matched__bc", &UpsilonDataToTree[100]);
        tree_upsilon->Branch("extraInfo__boEeclv200_900_unmatched__bc", &UpsilonDataToTree[101]);
        tree_upsilon->Branch("extraInfo__boNgammav200_925__bc", &UpsilonDataToTree[102]);
        tree_upsilon->Branch("extraInfo__boNgammav200_925_matched__bc", &UpsilonDataToTree[103]);
        tree_upsilon->Branch("extraInfo__boNgammav200_925_unmatched__bc", &UpsilonDataToTree[104]);
        tree_upsilon->Branch("extraInfo__boEeclv200_925__bc", &UpsilonDataToTree[105]);
        tree_upsilon->Branch("extraInfo__boEeclv200_925_matched__bc", &UpsilonDataToTree[106]);
        tree_upsilon->Branch("extraInfo__boEeclv200_925_unmatched__bc", &UpsilonDataToTree[107]);
        tree_upsilon->Branch("extraInfo__boNgammav200_950__bc", &UpsilonDataToTree[108]);
        tree_upsilon->Branch("extraInfo__boNgammav200_950_matched__bc", &UpsilonDataToTree[109]);
        tree_upsilon->Branch("extraInfo__boNgammav200_950_unmatched__bc", &UpsilonDataToTree[110]);
        tree_upsilon->Branch("extraInfo__boEeclv200_950__bc", &UpsilonDataToTree[111]);
        tree_upsilon->Branch("extraInfo__boEeclv200_950_matched__bc", &UpsilonDataToTree[112]);
        tree_upsilon->Branch("extraInfo__boEeclv200_950_unmatched__bc", &UpsilonDataToTree[113]);
        tree_upsilon->Branch("extraInfo__boNgammav200_975__bc", &UpsilonDataToTree[114]);
        tree_upsilon->Branch("extraInfo__boNgammav200_975_matched__bc", &UpsilonDataToTree[115]);
        tree_upsilon->Branch("extraInfo__boNgammav200_975_unmatched__bc", &UpsilonDataToTree[116]);
        tree_upsilon->Branch("extraInfo__boEeclv200_975__bc", &UpsilonDataToTree[117]);
        tree_upsilon->Branch("extraInfo__boEeclv200_975_matched__bc", &UpsilonDataToTree[118]);
        tree_upsilon->Branch("extraInfo__boEeclv200_975_unmatched__bc", &UpsilonDataToTree[119]);
        tree_upsilon->Branch("extraInfo__boNgammav200_025__bc", &UpsilonDataToTree[120]);
        tree_upsilon->Branch("extraInfo__boNgammav200_025_matched__bc", &UpsilonDataToTree[121]);
        tree_upsilon->Branch("extraInfo__boNgammav200_025_unmatched__bc", &UpsilonDataToTree[122]);
        tree_upsilon->Branch("extraInfo__boEeclv200_025__bc", &UpsilonDataToTree[123]);
        tree_upsilon->Branch("extraInfo__boEeclv200_025_matched__bc", &UpsilonDataToTree[124]);
        tree_upsilon->Branch("extraInfo__boEeclv200_025_unmatched__bc", &UpsilonDataToTree[125]);
        tree_upsilon->Branch("extraInfo__boNgammav200_050__bc", &UpsilonDataToTree[126]);
        tree_upsilon->Branch("extraInfo__boNgammav200_050_matched__bc", &UpsilonDataToTree[127]);
        tree_upsilon->Branch("extraInfo__boNgammav200_050_unmatched__bc", &UpsilonDataToTree[128]);
        tree_upsilon->Branch("extraInfo__boEeclv200_050__bc", &UpsilonDataToTree[129]);
        tree_upsilon->Branch("extraInfo__boEeclv200_050_matched__bc", &UpsilonDataToTree[130]);
        tree_upsilon->Branch("extraInfo__boEeclv200_050_unmatched__bc", &UpsilonDataToTree[131]);
        tree_upsilon->Branch("extraInfo__boNgammav200_075__bc", &UpsilonDataToTree[132]);
        tree_upsilon->Branch("extraInfo__boNgammav200_075_matched__bc", &UpsilonDataToTree[133]);
        tree_upsilon->Branch("extraInfo__boNgammav200_075_unmatched__bc", &UpsilonDataToTree[134]);
        tree_upsilon->Branch("extraInfo__boEeclv200_075__bc", &UpsilonDataToTree[135]);
        tree_upsilon->Branch("extraInfo__boEeclv200_075_matched__bc", &UpsilonDataToTree[136]);
        tree_upsilon->Branch("extraInfo__boEeclv200_075_unmatched__bc", &UpsilonDataToTree[137]);
        tree_upsilon->Branch("extraInfo__boNgammav200_100__bc", &UpsilonDataToTree[138]);
        tree_upsilon->Branch("extraInfo__boNgammav200_100_matched__bc", &UpsilonDataToTree[139]);
        tree_upsilon->Branch("extraInfo__boNgammav200_100_unmatched__bc", &UpsilonDataToTree[140]);
        tree_upsilon->Branch("extraInfo__boEeclv200_100__bc", &UpsilonDataToTree[141]);
        tree_upsilon->Branch("extraInfo__boEeclv200_100_matched__bc", &UpsilonDataToTree[142]);
        tree_upsilon->Branch("extraInfo__boEeclv200_100_unmatched__bc", &UpsilonDataToTree[143]);
        tree_upsilon->Branch("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &UpsilonDataToTree[144]);
        tree_upsilon->Branch("nParticlesInList__boD0__clDecayIntoKL0__bc", &UpsilonDataToTree[145]); 
        tree_upsilon->Branch("nParticlesInList__boB__pl__clKnn__bc", &UpsilonDataToTree[146]); 
        tree_upsilon->Branch("invMassInLists__bon0__clKnn__bc", &UpsilonDataToTree[147]);
        tree_upsilon->Branch("nParticlesInList__boB__pl__clKstarnn__bc", &UpsilonDataToTree[148]); 
        tree_upsilon->Branch("invMassInLists__bon0__clKstarnn__bc", &UpsilonDataToTree[149]);
        tree_upsilon->Branch("nParticlesInList__boB0__clK0nn__bc", &UpsilonDataToTree[150]); 
        tree_upsilon->Branch("invMassInLists__bon0__clK0nn__bc", &UpsilonDataToTree[151]);
        tree_upsilon->Branch("nParticlesInList__boB0__clKstar0nn__bc", &UpsilonDataToTree[152]);
        tree_upsilon->Branch("invMassInLists__bon0__clKstar0nn__bc", &UpsilonDataToTree[153]);
        tree_upsilon->Branch("nParticlesInList__boB__pl__clKpKLKL_all__bc", &UpsilonDataToTree[154]);
        tree_upsilon->Branch("nParticlesInList__boB__pl__clKpKLKL_NR__bc", &UpsilonDataToTree[155]);
        tree_upsilon->Branch("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &UpsilonDataToTree[156]);
        tree_upsilon->Branch("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &UpsilonDataToTree[157]);
        tree_upsilon->Branch("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &UpsilonDataToTree[158]);
        tree_upsilon->Branch("nParticlesInList__boB0__clKSKLKL_all__bc", &UpsilonDataToTree[159]);
        tree_upsilon->Branch("nParticlesInList__boB0__clKSKLKL_NR__bc", &UpsilonDataToTree[160]);
        tree_upsilon->Branch("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &UpsilonDataToTree[161]);
        tree_upsilon->Branch("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &UpsilonDataToTree[162]);
        tree_upsilon->Branch("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &UpsilonDataToTree[163]);

        // get Bsig_info
        tree_Bsig->Branch("Bsig_E", &BsigDataToTree[0]);
        tree_Bsig->Branch("Bsig_useCMSFrame_E", &BsigDataToTree[1]);
        tree_upsilon->Branch("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc", &BsigDataToTree[2]);
        tree_Bsig->Branch("Bsig_p", &BsigDataToTree[3]);
        tree_Bsig->Branch("Bsig_useCMSFrame_p", &BsigDataToTree[4]);
        tree_upsilon->Branch("useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc", &BsigDataToTree[5]);
        tree_Bsig->Branch("Bsig_M", &BsigDataToTree[6]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_KS0_3D_distance", &BsigDataToTree[7]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nKslow1", &BsigDataToTree[8]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nKslow2", &BsigDataToTree[9]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nKslow3", &BsigDataToTree[10]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nPislow1", &BsigDataToTree[11]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nPislow2", &BsigDataToTree[12]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nPislow3", &BsigDataToTree[13]);
        tree_Bsig->Branch("Bsig_R2", &BsigDataToTree[14]);
        tree_Bsig->Branch("Bsig_thrustBm", &BsigDataToTree[15]);
        tree_Bsig->Branch("Bsig_thrustOm", &BsigDataToTree[16]);
        tree_Bsig->Branch("Bsig_cosTBTO", &BsigDataToTree[17]);
        tree_Bsig->Branch("Bsig_cosTBz", &BsigDataToTree[18]);
        tree_Bsig->Branch("Bsig_KSFWVariables_et", &BsigDataToTree[19]);
        tree_Bsig->Branch("Bsig_KSFWVariables_mm2", &BsigDataToTree[20]);
        tree_Bsig->Branch("Bsig_KSFWVariables_hso00", &BsigDataToTree[21]);
        tree_Bsig->Branch("Bsig_KSFWVariables_hso01", &BsigDataToTree[22]);
        tree_Bsig->Branch("Bsig_KSFWVariables_hso02", &BsigDataToTree[23]);
        tree_Bsig->Branch("Bsig_KSFWVariables_hso03", &BsigDataToTree[24]);
        tree_Bsig->Branch("Bsig_KSFWVariables_hso04", &BsigDataToTree[25]);
        tree_Bsig->Branch("Bsig_KSFWVariables_hso10", &BsigDataToTree[26]);
        tree_Bsig->Branch("Bsig_KSFWVariables_hso12", &BsigDataToTree[27]);
        tree_Bsig->Branch("Bsig_KSFWVariables_hso14", &BsigDataToTree[28]);
        tree_Bsig->Branch("Bsig_KSFWVariables_hso20", &BsigDataToTree[29]);
        tree_Bsig->Branch("Bsig_KSFWVariables_hso22", &BsigDataToTree[30]);
        tree_Bsig->Branch("Bsig_KSFWVariables_hso24", &BsigDataToTree[31]);
        tree_Bsig->Branch("Bsig_KSFWVariables_hoo0", &BsigDataToTree[32]);
        tree_Bsig->Branch("Bsig_KSFWVariables_hoo1", &BsigDataToTree[33]);
        tree_Bsig->Branch("Bsig_KSFWVariables_hoo2", &BsigDataToTree[34]);
        tree_Bsig->Branch("Bsig_KSFWVariables_hoo3", &BsigDataToTree[35]);
        tree_Bsig->Branch("Bsig_KSFWVariables_hoo4", &BsigDataToTree[36]);
        tree_Bsig->Branch("Bsig_CleoConeCS_1", &BsigDataToTree[37]);
        tree_Bsig->Branch("Bsig_CleoConeCS_2", &BsigDataToTree[38]);
        tree_Bsig->Branch("Bsig_CleoConeCS_3", &BsigDataToTree[39]);
        tree_Bsig->Branch("Bsig_CleoConeCS_4", &BsigDataToTree[40]);
        tree_Bsig->Branch("Bsig_CleoConeCS_5", &BsigDataToTree[41]);
        tree_Bsig->Branch("Bsig_CleoConeCS_6", &BsigDataToTree[42]);
        tree_Bsig->Branch("Bsig_CleoConeCS_7", &BsigDataToTree[43]);
        tree_Bsig->Branch("Bsig_CleoConeCS_8", &BsigDataToTree[44]);
        tree_Bsig->Branch("Bsig_CleoConeCS_9", &BsigDataToTree[45]);
        tree_Bsig->Branch("Bsig_pt", &BsigDataToTree[46]);
        tree_Bsig->Branch("Bsig_useCMSFrame_pt", &BsigDataToTree[47]);
        tree_Bsig->Branch("Bsig_theta", &BsigDataToTree[48]);
        tree_Bsig->Branch("Bsig_useCMSFrame_theta", &BsigDataToTree[49]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nDc", &BsigDataToTree[50]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Dc_pValue_med", &BsigDataToTree[51]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Dc_pValue_std", &BsigDataToTree[52]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Dcsimpleveto_chiProb", &BsigDataToTree[53]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Dcsimpleveto_dr", &BsigDataToTree[54]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Dcsimpleveto_dz", &BsigDataToTree[55]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Dcsimpleveto_M", &BsigDataToTree[56]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nD0", &BsigDataToTree[57]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_D0_pValue_med", &BsigDataToTree[58]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_D0_pValue_std", &BsigDataToTree[59]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_D0simpleveto_chiProb", &BsigDataToTree[60]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_D0simpleveto_dr", &BsigDataToTree[61]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_D0simpleveto_dz", &BsigDataToTree[62]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_D0simpleveto_M", &BsigDataToTree[63]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_mychiProb", &BsigDataToTree[64]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_mydr", &BsigDataToTree[65]);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_mydz", &BsigDataToTree[66]);
        for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
            tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_nKtruebin" + std::to_string(i_PID)).c_str(), &BsigDataToTree[67 + 4 * i_PID]);
            tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_nKmisbin" + std::to_string(i_PID)).c_str(), &BsigDataToTree[67 + 4 * i_PID + 1]);
            tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_npitruebin" + std::to_string(i_PID)).c_str(), &BsigDataToTree[67 + 4 * i_PID + 2]);
            tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_npimisbin" + std::to_string(i_PID)).c_str(), &BsigDataToTree[67 + 4 * i_PID + 3]);
        }
        for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_npi0bin" + std::to_string(i_pi0)).c_str(), &BsigDataToTree[359 + i_pi0]);
        for (int i_PID = 0; i_PID < N_fakeE_syst; i_PID++) {
            tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_nKfakeEbin_n" + std::to_string(i_PID)).c_str(), &BsigDataToTree[367 + 4 * i_PID]);
            tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_nKfakeEbin_p" + std::to_string(i_PID)).c_str(), &BsigDataToTree[367 + 4 * i_PID + 1]);
            tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_npifakeEbin_n" + std::to_string(i_PID)).c_str(), &BsigDataToTree[367 + 4 * i_PID + 2]);
            tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_npifakeEbin_p" + std::to_string(i_PID)).c_str(), &BsigDataToTree[367 + 4 * i_PID + 3]);
        }
        for (int i_PID = 0; i_PID < N_fakeMU_syst; i_PID++) {
            tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_nKfakeMUbin_n" + std::to_string(i_PID)).c_str(), &BsigDataToTree[542 + 4 * i_PID]);
            tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_nKfakeMUbin_p" + std::to_string(i_PID)).c_str(), &BsigDataToTree[542 + 4 * i_PID + 1]);
            tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_npifakeMUbin_n" + std::to_string(i_PID)).c_str(), &BsigDataToTree[542 + 4 * i_PID + 2]);
            tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_npifakeMUbin_p" + std::to_string(i_PID)).c_str(), &BsigDataToTree[542 + 4 * i_PID + 3]);
        }

        // get Btag_info
        tree_Btag->Branch("Btag_extraInfo_decayModeID", &BtagDataToTree[0]);
        tree_Btag->Branch("Btag_Mbc", &BtagDataToTree[1]);
        tree_Btag->Branch("Btag_deltaE", &BtagDataToTree[2]);
        tree_Btag->Branch("Btag_E", &BtagDataToTree[3]);
        tree_Btag->Branch("Btag_useCMSFrame_E", &BtagDataToTree[4]);
        tree_Btag->Branch("Btag_extraInfo_SignalProbability", &BtagDataToTree[5]);
        tree_Btag->Branch("Btag_chiProb", &BtagDataToTree[6]);
        tree_Btag->Branch("Btag_dr", &BtagDataToTree[7]);
        tree_Btag->Branch("Btag_dz", &BtagDataToTree[8]);
        tree_Btag->Branch("Btag_useCMSFrame_p", &BtagDataToTree[9]);
        tree_Btag->Branch("Btag_useCMSFrame_phi", &BtagDataToTree[10]);

        // other information I need
        tree_Btag->Branch("Btag_R2", &DataToTree[0]);
        tree_Btag->Branch("Btag_thrustBm", &DataToTree[1]);
        tree_Btag->Branch("Btag_thrustOm", &DataToTree[2]);
        tree_Btag->Branch("Btag_cosTBTO", &DataToTree[3]);
        tree_Btag->Branch("Btag_cosTBz", &DataToTree[4]);
        tree_Btag->Branch("Btag_KSFWVariables_et", &DataToTree[5]);
        tree_Btag->Branch("Btag_KSFWVariables_mm2", &DataToTree[6]);
        tree_Btag->Branch("Btag_KSFWVariables_hso00", &DataToTree[7]);
        tree_Btag->Branch("Btag_KSFWVariables_hso01", &DataToTree[8]);
        tree_Btag->Branch("Btag_KSFWVariables_hso02", &DataToTree[9]);
        tree_Btag->Branch("Btag_KSFWVariables_hso03", &DataToTree[10]);
        tree_Btag->Branch("Btag_KSFWVariables_hso04", &DataToTree[11]);
        tree_Btag->Branch("Btag_KSFWVariables_hso10", &DataToTree[12]);
        tree_Btag->Branch("Btag_KSFWVariables_hso12", &DataToTree[13]);
        tree_Btag->Branch("Btag_KSFWVariables_hso14", &DataToTree[14]);
        tree_Btag->Branch("Btag_KSFWVariables_hso20", &DataToTree[15]);
        tree_Btag->Branch("Btag_KSFWVariables_hso22", &DataToTree[16]);
        tree_Btag->Branch("Btag_KSFWVariables_hso24", &DataToTree[17]);
        tree_Btag->Branch("Btag_KSFWVariables_hoo0", &DataToTree[18]);
        tree_Btag->Branch("Btag_KSFWVariables_hoo1", &DataToTree[19]);
        tree_Btag->Branch("Btag_KSFWVariables_hoo2", &DataToTree[20]);
        tree_Btag->Branch("Btag_KSFWVariables_hoo3", &DataToTree[21]);
        tree_Btag->Branch("Btag_KSFWVariables_hoo4", &DataToTree[22]);
        tree_Btag->Branch("Btag_CleoConeCS_1", &DataToTree[23]);
        tree_Btag->Branch("Btag_CleoConeCS_2", &DataToTree[24]);
        tree_Btag->Branch("Btag_CleoConeCS_3", &DataToTree[25]);
        tree_Btag->Branch("Btag_CleoConeCS_4", &DataToTree[26]);
        tree_Btag->Branch("Btag_CleoConeCS_5", &DataToTree[27]);
        tree_Btag->Branch("Btag_CleoConeCS_6", &DataToTree[28]);
        tree_Btag->Branch("Btag_CleoConeCS_7", &DataToTree[29]);
        tree_Btag->Branch("Btag_CleoConeCS_8", &DataToTree[30]);
        tree_Btag->Branch("Btag_CleoConeCS_9", &DataToTree[31]);
        tree_upsilon->Branch("missingMass2OfEvent", &DataToTree[32]);
        tree_upsilon->Branch("visibleEnergyOfEventCMS", &DataToTree[33]);
        tree_Btag->Branch("Btag_useCMSFrame_theta", &DataToTree[34]);
        tree_upsilon->Branch("extraInfo__boDecayHash__bc", &DataToTree[35]);
        tree_upsilon->Branch("extraInfo__boDecayHashExtended__bc", &DataToTree[36]);

        // decay mode (MC level)
        tree_Xs->Branch("nParticlesInList__boB__pl__clKcharge_total__bc", &DecayDataToTree[0]);
        tree_Xs->Branch("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &DecayDataToTree[1]);
        tree_Xs->Branch("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &DecayDataToTree[2]);
        tree_Xs->Branch("nParticlesInList__boXsu__clMCcomb__bc", &DecayDataToTree[3]);
        tree_Xs->Branch("nParticlesInList__boXsu__clMCch1__bc", &DecayDataToTree[4]);
        tree_Xs->Branch("nParticlesInList__boXsu__clMCch2__bc", &DecayDataToTree[5]);
        tree_Xs->Branch("nParticlesInList__boXsu__clMCch3__bc", &DecayDataToTree[6]);
        tree_Xs->Branch("nParticlesInList__boXsu__clMCch4__bc", &DecayDataToTree[7]);
        tree_Xs->Branch("nParticlesInList__boXsu__clMCch5__bc", &DecayDataToTree[8]);
        tree_Xs->Branch("nParticlesInList__boXsu__clMCch6__bc", &DecayDataToTree[9]);
        tree_Xs->Branch("nParticlesInList__boXsu__clMCch7__bc", &DecayDataToTree[10]);
        tree_Xs->Branch("nParticlesInList__boXsu__clMCch8__bc", &DecayDataToTree[11]);
        tree_Xs->Branch("nParticlesInList__boXsu__clMCch9__bc", &DecayDataToTree[12]);
        tree_Xs->Branch("nParticlesInList__boXsu__clMCch10__bc", &DecayDataToTree[13]);
        tree_Xs->Branch("nParticlesInList__boXsu__clMCch11__bc", &DecayDataToTree[14]);
        tree_Xs->Branch("nParticlesInList__boXsu__clMCch12__bc", &DecayDataToTree[15]);
        tree_Xs->Branch("nParticlesInList__boXsu__clMCch13__bc", &DecayDataToTree[16]);
        tree_Xs->Branch("nParticlesInList__boXsu__clMCch14__bc", &DecayDataToTree[17]);
        tree_Xs->Branch("nParticlesInList__boXsu__clMCch15__bc", &DecayDataToTree[18]);
        tree_Xs->Branch("nParticlesInList__boB0__clKneutral_total__bc", &DecayDataToTree[19]);
        tree_Xs->Branch("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &DecayDataToTree[20]);
        tree_Xs->Branch("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &DecayDataToTree[21]);
        tree_Xs->Branch("nParticlesInList__boXsd__clMCcomb__bc", &DecayDataToTree[22]);
        tree_Xs->Branch("nParticlesInList__boXsd__clMCch16__bc", &DecayDataToTree[23]);
        tree_Xs->Branch("nParticlesInList__boXsd__clMCch17__bc", &DecayDataToTree[24]);
        tree_Xs->Branch("nParticlesInList__boXsd__clMCch18__bc", &DecayDataToTree[25]);
        tree_Xs->Branch("nParticlesInList__boXsd__clMCch19__bc", &DecayDataToTree[26]);
        tree_Xs->Branch("nParticlesInList__boXsd__clMCch20__bc", &DecayDataToTree[27]);
        tree_Xs->Branch("nParticlesInList__boXsd__clMCch21__bc", &DecayDataToTree[28]);
        tree_Xs->Branch("nParticlesInList__boXsd__clMCch22__bc", &DecayDataToTree[29]);
        tree_Xs->Branch("nParticlesInList__boXsd__clMCch23__bc", &DecayDataToTree[30]);
        tree_Xs->Branch("nParticlesInList__boXsd__clMCch24__bc", &DecayDataToTree[31]);
        tree_Xs->Branch("nParticlesInList__boXsd__clMCch25__bc", &DecayDataToTree[32]);
        tree_Xs->Branch("nParticlesInList__boXsd__clMCch26__bc", &DecayDataToTree[33]);
        tree_Xs->Branch("nParticlesInList__boXsd__clMCch27__bc", &DecayDataToTree[34]);
        tree_Xs->Branch("nParticlesInList__boXsd__clMCch28__bc", &DecayDataToTree[35]);
        tree_Xs->Branch("nParticlesInList__boXsd__clMCch29__bc", &DecayDataToTree[36]);
        tree_Xs->Branch("nParticlesInList__boXsd__clMCch30__bc", &DecayDataToTree[37]);
        tree_Xs->Branch("nParticlesInList__bonu_e__clMC_signal__bc", &DecayNparticlesDataToTree[0]);
        tree_Xs->Branch("nParticlesInList__boB__pl__clMC_signal_total_e__bc", &DecayNparticlesDataToTree[1]);
        tree_Xs->Branch("nParticlesInList__boB0__clMC_signal_total_e__bc", &DecayNparticlesDataToTree[2]);
        tree_Xs->Branch("nParticlesInList__boB__pl__clPrimaryMC__bc", &DecayNparticlesDataToTree[3]);
        tree_Xs->Branch("nParticlesInList__boB0__clPrimaryMC__bc", &DecayNparticlesDataToTree[4]);
        tree_Xs->Branch("invMassInLists__bonu_e__clMC_signal__bc", &DecaySystFFDataToTree[index_q2]);
        tree_Xs->Branch("averageValueInList__boB__pl__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", &DecaySystFFDataToTree[1]);
        tree_Xs->Branch("averageValueInList__boB0__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", &DecaySystFFDataToTree[2]);
        tree_Xs->Branch("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &DecaySystFFDataToTree[3]);
        tree_Xs->Branch("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &DecaySystFFDataToTree[4]);
        tree_Xs->Branch("averageValueInList__boB__pl__clMC_signal_total_e__cm__spM__bc", &DecaySystFFDataToTree[5]);
        tree_Xs->Branch("averageValueInList__boB0__clMC_signal_total_e__cm__spM__bc", &DecaySystFFDataToTree[6]);

        tree_upsilon->Branch("MVA_BB", &MVA_BB_DataToTree);
        tree_upsilon->Branch("MVA_Continuum", &MVA_Continuum_DataToTree);

        tree_upsilon->Branch("nParticlesInList__bomu__pl__clfromUpsilonmychargedMuon__bc", &nROE_mu_DataToTree);
        tree_upsilon->Branch("nParticlesInList__boJ__slpsi__clfromUpsilontemp__bc", &nROE_Jpsi_DataToTree);
        tree_upsilon->Branch("nParticlesInList__boUpsilon__bo4S__bc__cltemp_withoutneutrino_BCS__bc", &nROE_Upsilon_BCS_DataToTree);
        tree_upsilon->Branch("nParticlesInList__boUpsilon__bo4S__bc__cltemp_withoutneutrino__bc", &nROE_Upsilon_DataToTree);
        /*================================================================*/
        files.push_back(file);
        trees_upsilon.push_back(tree_upsilon);
        trees_Bsig.push_back(tree_Bsig);
        trees_Btag.push_back(tree_Btag);
        trees_Xs.push_back(tree_Xs);
    }
    else if (files.size() > current_file && trees_upsilon.size() > current_file && files.size() == trees_upsilon.size() && files.size() == trees_Bsig.size() && files.size() == trees_Btag.size() && files.size() == trees_Xs.size()) { // use what I have
    }
    else { // error
        printf("ERROR! 035\n");
        exit(1);
    }

    TFile* temp_file = files.at(current_file);
    temp_file->cd();
    TTree* temp_tree_upsilon = trees_upsilon.at(current_file);
    TTree* temp_tree_Bsig = trees_Bsig.at(current_file);
    TTree* temp_tree_Btag = trees_Btag.at(current_file);
    TTree* temp_tree_Xs = trees_Xs.at(current_file);
    std::queue<Data> temp_queue;
    temp_queue.swap(TotalData);
    while (!temp_queue.empty()) {
        Data temp = temp_queue.front();
        temp_queue.pop();

        {
            UpsilonExperimentToTree = temp.upsilon_experiment;
            UpsilonRunToTree = temp.upsilon_run;
            UpsilonEventToTree = temp.upsilon_event;
            UpsilonCandidateToTree = temp.upsilon_candidate;
            UpsilonNcandidatesToTree = temp.upsilon_ncandidates;
            BsigExperimentToTree = temp.Bsig_experiment;
            BsigRunToTree = temp.Bsig_run;
            BsigEventToTree = temp.Bsig_event;
            BsigCandidateToTree = temp.Bsig_candidate;
            BsigNcandidatesToTree = temp.Bsig_ncandidates;
            BtagExperimentToTree = temp.Btag_experiment;
            BtagRunToTree = temp.Btag_run;
            BtagEventToTree = temp.Btag_event;
            BtagCandidateToTree = temp.Btag_candidate;
            BtagNcandidatesToTree = temp.Btag_ncandidates;
        }
        for (int i = 0; i < N_Upsilon_info; i++) {
            UpsilonDataToTree[i] = temp.Upsilon_info[i];
        }
        for (int i = 0; i < N_Bsig_info; i++) {
            BsigDataToTree[i] = temp.Bsig_info[i];
        }
        for (int i = 0; i < N_Btag_info; i++) {
            BtagDataToTree[i] = temp.Btag_info[i];
        }
        for (int i = 0; i < N_decay; i++) {
            DecayDataToTree[i] = temp.Decay[i];
        }
        for (int i = 0; i < N_decay_nparticles; i++) {
            DecayNparticlesDataToTree[i] = temp.Decay_nparticles[i];
        }
        for (int i = 0; i < N_decay_syst_ff; i++) {
            DecaySystFFDataToTree[i] = temp.Decay_syst_ff[i];
        }
        for (int i = 0; i < N_Needed_info; i++) {
            DataToTree[i] = temp.Needed_info[i];
        }
        Upsilon_decayIDToTree = temp.Upsilon_decayID;
        Bsig_decayIDToTree = temp.Bsig_decayID;

        MVA_BB_DataToTree = temp.MVA_BB;
        MVA_Continuum_DataToTree = temp.MVA_Continuum;

        nROE_mu_DataToTree = temp.nROE_mu;
        nROE_Jpsi_DataToTree = temp.nROE_Jpsi;
        nROE_Upsilon_BCS_DataToTree = temp.nROE_Upsilon_BCS;
        nROE_Upsilon_DataToTree = temp.nROE_Upsilon;

        temp_tree_upsilon->Fill();
        temp_tree_Bsig->Fill();
        temp_tree_Btag->Fill();
        temp_tree_Xs->Fill();

        TotalData.push(temp);
    }

    current_file++;
}

void Loader::PrintSeparateRootFile(std::string output_name) {

    TFile* temp_file = new TFile(output_name.c_str(), "recreate");
    temp_file->cd();
    TTree* temp_tree_upsilon = new TTree("Upsilon", "");
    TTree* temp_tree_Bsig = new TTree("Bsig", "");
    TTree* temp_tree_Btag = new TTree("Btag", "");

    TTree* temp_tree_Xs;
    if (DoesItHaveXsBranch) temp_tree_Xs = new TTree("Xs", "");
    else temp_tree_Xs = nullptr;

    int temp_UpsilonExperimentToTree;
    int temp_UpsilonRunToTree;
    unsigned int temp_UpsilonEventToTree;
    int temp_UpsilonCandidateToTree;
    int temp_UpsilonNcandidatesToTree;
    int temp_BsigExperimentToTree;
    int temp_BsigRunToTree;
    unsigned int temp_BsigEventToTree;
    int temp_BsigCandidateToTree;
    int temp_BsigNcandidatesToTree;
    int temp_BtagExperimentToTree;
    int temp_BtagRunToTree;
    unsigned int temp_BtagEventToTree;
    int temp_BtagCandidateToTree;
    int temp_BtagNcandidatesToTree;
    double temp_UpsilonDataToTree[N_Upsilon_info];
    double temp_BsigDataToTree[N_Bsig_info];
    double temp_BtagDataToTree[N_Btag_info];
    int temp_DecayDataToTree[N_decay];
    int temp_DecayNparticlesDataToTree[N_decay_nparticles];
    double temp_DecaySystFFDataToTree[N_decay_syst_ff];
    double temp_DataToTree[N_Needed_info];
    double temp_Upsilon_decayIDToTree;
    double temp_Bsig_decayIDToTree;

    float temp_MVA_BB_DataToTree;
    float temp_MVA_Continuum_DataToTree;

    int temp_nROE_mu;
    int temp_nROE_Jpsi;
    int temp_nROE_Upsilon_BCS;
    int temp_nROE_Upsilon;

    /*================================================================*/
    // get event_info
    temp_tree_upsilon->Branch("__experiment__", &temp_UpsilonExperimentToTree);
    temp_tree_upsilon->Branch("__run__", &temp_UpsilonRunToTree);
    temp_tree_upsilon->Branch("__event__", &temp_UpsilonEventToTree);
    temp_tree_upsilon->Branch("__candidate__", &temp_UpsilonCandidateToTree);
    temp_tree_upsilon->Branch("__ncandidates__", &temp_UpsilonNcandidatesToTree);
    temp_tree_Bsig->Branch("__experiment__", &temp_BsigExperimentToTree);
    temp_tree_Bsig->Branch("__run__", &temp_BsigRunToTree);
    temp_tree_Bsig->Branch("__event__", &temp_BsigEventToTree);
    temp_tree_Bsig->Branch("__candidate__", &temp_BsigCandidateToTree);
    temp_tree_Bsig->Branch("__ncandidates__", &temp_BsigNcandidatesToTree);
    temp_tree_Btag->Branch("__experiment__", &temp_BtagExperimentToTree);
    temp_tree_Btag->Branch("__run__", &temp_BtagRunToTree);
    temp_tree_Btag->Branch("__event__", &temp_BtagEventToTree);
    temp_tree_Btag->Branch("__candidate__", &temp_BtagCandidateToTree);
    temp_tree_Btag->Branch("__ncandidates__", &temp_BtagNcandidatesToTree);

    // get decaymodeID
    temp_tree_upsilon->Branch("extraInfo__bodecayModeID__bc", &temp_Upsilon_decayIDToTree);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_decayModeID", &temp_Bsig_decayIDToTree);

    // get Upsilon_info
    temp_tree_upsilon->Branch("nROE_ECLClusters__bocleanMask__bc", &temp_UpsilonDataToTree[0]);
    temp_tree_upsilon->Branch("nROE_KLMClusters", &temp_UpsilonDataToTree[1]);
    temp_tree_upsilon->Branch("nROE_Tracks__bocleanMask__bc", &temp_UpsilonDataToTree[2]);
    temp_tree_upsilon->Branch("roeEextra__bocleanMask__bc", &temp_UpsilonDataToTree[3]);
    temp_tree_upsilon->Branch("nROE_NeutralECLClusters__bocleanMask__bc", &temp_UpsilonDataToTree[4]);
    temp_tree_upsilon->Branch("nROE_ParticlesInList__boK_S0__clmyKaonshort__bc", &temp_UpsilonDataToTree[5]);
    temp_tree_upsilon->Branch("nROE_ParticlesInList__bopi0__clmyneutralPion__bc", &temp_UpsilonDataToTree[6]);
    temp_tree_upsilon->Branch("missingMomentumOfEvent_theta", &temp_UpsilonDataToTree[7]);
    temp_tree_upsilon->Branch("missingMomentumOfEvent", &temp_UpsilonDataToTree[8]);
    temp_tree_upsilon->Branch("missingEnergyOfEventCMS", &temp_UpsilonDataToTree[9]);
    temp_tree_upsilon->Branch("nRemainingTracksInEvent", &temp_UpsilonDataToTree[10]);
    temp_tree_upsilon->Branch("roeNeextra__bocleanMask__bc", &temp_UpsilonDataToTree[11]);
    temp_tree_upsilon->Branch("useCMSFrame__boroeNeextra__bocleanMask__bc__bc", &temp_UpsilonDataToTree[12]);
    temp_tree_upsilon->Branch("nROE_ParticlesInList__bogamma__clmygamma__bc", &temp_UpsilonDataToTree[13]);
    temp_tree_upsilon->Branch("foxWolframR1", &temp_UpsilonDataToTree[14]);
    temp_tree_upsilon->Branch("foxWolframR2", &temp_UpsilonDataToTree[15]);
    temp_tree_upsilon->Branch("foxWolframR3", &temp_UpsilonDataToTree[16]);
    temp_tree_upsilon->Branch("foxWolframR4", &temp_UpsilonDataToTree[17]);
    temp_tree_upsilon->Branch("harmonicMomentThrust0", &temp_UpsilonDataToTree[18]);
    temp_tree_upsilon->Branch("harmonicMomentThrust1", &temp_UpsilonDataToTree[19]);
    temp_tree_upsilon->Branch("harmonicMomentThrust2", &temp_UpsilonDataToTree[20]);
    temp_tree_upsilon->Branch("harmonicMomentThrust3", &temp_UpsilonDataToTree[21]);
    temp_tree_upsilon->Branch("harmonicMomentThrust4", &temp_UpsilonDataToTree[22]);
    temp_tree_upsilon->Branch("cleoConeThrust0", &temp_UpsilonDataToTree[23]);
    temp_tree_upsilon->Branch("cleoConeThrust1", &temp_UpsilonDataToTree[24]);
    temp_tree_upsilon->Branch("cleoConeThrust2", &temp_UpsilonDataToTree[25]);
    temp_tree_upsilon->Branch("cleoConeThrust3", &temp_UpsilonDataToTree[26]);
    temp_tree_upsilon->Branch("cleoConeThrust4", &temp_UpsilonDataToTree[27]);
    temp_tree_upsilon->Branch("cleoConeThrust5", &temp_UpsilonDataToTree[28]);
    temp_tree_upsilon->Branch("cleoConeThrust6", &temp_UpsilonDataToTree[29]);
    temp_tree_upsilon->Branch("cleoConeThrust7", &temp_UpsilonDataToTree[30]);
    temp_tree_upsilon->Branch("cleoConeThrust8", &temp_UpsilonDataToTree[31]);
    temp_tree_upsilon->Branch("sphericity", &temp_UpsilonDataToTree[32]);
    temp_tree_upsilon->Branch("aplanarity", &temp_UpsilonDataToTree[33]);
    temp_tree_upsilon->Branch("thrust", &temp_UpsilonDataToTree[34]);
    temp_tree_upsilon->Branch("thrustAxisCosTheta", &temp_UpsilonDataToTree[35]);
    temp_tree_upsilon->Branch("MsquaredBsig_op0", &temp_UpsilonDataToTree[36]);
    temp_tree_upsilon->Branch("MsquaredBsig_op1", &temp_UpsilonDataToTree[37]);
    temp_tree_upsilon->Branch("MsquaredBsig_op2", &temp_UpsilonDataToTree[38]);
    temp_tree_upsilon->Branch("MsquaredBsig_op3", &temp_UpsilonDataToTree[39]);
    temp_tree_upsilon->Branch("MsquaredBsig_op4", &temp_UpsilonDataToTree[40]);
    temp_tree_upsilon->Branch("MsquaredBsig_op7", &temp_UpsilonDataToTree[41]);
    temp_tree_upsilon->Branch("roeP__bocleanMask__bc", &temp_UpsilonDataToTree[42]);
    temp_tree_upsilon->Branch("roeM__bocleanMask__bc", &temp_UpsilonDataToTree[43]);
    temp_tree_upsilon->Branch("roePTheta__bocleanMask__bc", &temp_UpsilonDataToTree[44]);
    temp_tree_upsilon->Branch("qsquared", &temp_UpsilonDataToTree[45]);
    temp_tree_upsilon->Branch("chiProb", &temp_UpsilonDataToTree[46]);
    temp_tree_upsilon->Branch("dr", &temp_UpsilonDataToTree[47]);
    temp_tree_upsilon->Branch("dz", &temp_UpsilonDataToTree[48]);
    temp_tree_upsilon->Branch("nParticlesInList__boe__pl__clElectronFBDT__bc", &temp_UpsilonDataToTree[49]);
    temp_tree_upsilon->Branch("nParticlesInList__bomu__pl__clMuonFBDT__bc", &temp_UpsilonDataToTree[50]);
    temp_tree_upsilon->Branch("nParticlesInList__boe__pl__clElectronFBDT_loose__bc", &temp_UpsilonDataToTree[51]);
    temp_tree_upsilon->Branch("nParticlesInList__bomu__pl__clMuonFBDT_loose__bc", &temp_UpsilonDataToTree[52]);
    temp_tree_upsilon->Branch("nParticlesInList__boe__pl__clElectronFBDT_tight__bc", &temp_UpsilonDataToTree[53]);
    temp_tree_upsilon->Branch("nParticlesInList__bomu__pl__clMuonFBDT_tight__bc", &temp_UpsilonDataToTree[54]);
    temp_tree_upsilon->Branch("beamE", &temp_UpsilonDataToTree[55]);
    temp_tree_upsilon->Branch("nROE_Tracks__bolooseMask__bc", &temp_UpsilonDataToTree[56]);
    temp_tree_upsilon->Branch("Ecms", &temp_UpsilonDataToTree[57]);
    temp_tree_upsilon->Branch("extraInfo__boNgamma__bc", &temp_UpsilonDataToTree[58]);
    temp_tree_upsilon->Branch("extraInfo__boEecl__bc", &temp_UpsilonDataToTree[59]);
    temp_tree_upsilon->Branch("extraInfo__boEecl_matched__bc", &temp_UpsilonDataToTree[60]);
    temp_tree_upsilon->Branch("extraInfo__boEecl_unmatched__bc", &temp_UpsilonDataToTree[61]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav133__bc", &temp_UpsilonDataToTree[62]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv133__bc", &temp_UpsilonDataToTree[63]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv133_matched__bc", &temp_UpsilonDataToTree[64]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv133_unmatched__bc", &temp_UpsilonDataToTree[65]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200__bc", &temp_UpsilonDataToTree[66]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_matched__bc", &temp_UpsilonDataToTree[67]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_unmatched__bc", &temp_UpsilonDataToTree[68]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200__bc", &temp_UpsilonDataToTree[69]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_matched__bc", &temp_UpsilonDataToTree[70]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_unmatched__bc", &temp_UpsilonDataToTree[71]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_800__bc", &temp_UpsilonDataToTree[72]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_800_matched__bc", &temp_UpsilonDataToTree[73]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_800_unmatched__bc", &temp_UpsilonDataToTree[74]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_800__bc", &temp_UpsilonDataToTree[75]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_800_matched__bc", &temp_UpsilonDataToTree[76]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_800_unmatched__bc", &temp_UpsilonDataToTree[77]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_825__bc", &temp_UpsilonDataToTree[78]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_825_matched__bc", &temp_UpsilonDataToTree[79]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_825_unmatched__bc", &temp_UpsilonDataToTree[80]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_825__bc", &temp_UpsilonDataToTree[81]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_825_matched__bc", &temp_UpsilonDataToTree[82]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_825_unmatched__bc", &temp_UpsilonDataToTree[83]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_850__bc", &temp_UpsilonDataToTree[84]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_850_matched__bc", &temp_UpsilonDataToTree[85]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_850_unmatched__bc", &temp_UpsilonDataToTree[86]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_850__bc", &temp_UpsilonDataToTree[87]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_850_matched__bc", &temp_UpsilonDataToTree[88]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_850_unmatched__bc", &temp_UpsilonDataToTree[89]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_875__bc", &temp_UpsilonDataToTree[90]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_875_matched__bc", &temp_UpsilonDataToTree[91]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_875_unmatched__bc", &temp_UpsilonDataToTree[92]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_875__bc", &temp_UpsilonDataToTree[93]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_875_matched__bc", &temp_UpsilonDataToTree[94]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_875_unmatched__bc", &temp_UpsilonDataToTree[95]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_900__bc", &temp_UpsilonDataToTree[96]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_900_matched__bc", &temp_UpsilonDataToTree[97]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_900_unmatched__bc", &temp_UpsilonDataToTree[98]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_900__bc", &temp_UpsilonDataToTree[99]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_900_matched__bc", &temp_UpsilonDataToTree[100]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_900_unmatched__bc", &temp_UpsilonDataToTree[101]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_925__bc", &temp_UpsilonDataToTree[102]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_925_matched__bc", &temp_UpsilonDataToTree[103]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_925_unmatched__bc", &temp_UpsilonDataToTree[104]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_925__bc", &temp_UpsilonDataToTree[105]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_925_matched__bc", &temp_UpsilonDataToTree[106]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_925_unmatched__bc", &temp_UpsilonDataToTree[107]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_950__bc", &temp_UpsilonDataToTree[108]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_950_matched__bc", &temp_UpsilonDataToTree[109]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_950_unmatched__bc", &temp_UpsilonDataToTree[110]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_950__bc", &temp_UpsilonDataToTree[111]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_950_matched__bc", &temp_UpsilonDataToTree[112]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_950_unmatched__bc", &temp_UpsilonDataToTree[113]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_975__bc", &temp_UpsilonDataToTree[114]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_975_matched__bc", &temp_UpsilonDataToTree[115]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_975_unmatched__bc", &temp_UpsilonDataToTree[116]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_975__bc", &temp_UpsilonDataToTree[117]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_975_matched__bc", &temp_UpsilonDataToTree[118]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_975_unmatched__bc", &temp_UpsilonDataToTree[119]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_025__bc", &temp_UpsilonDataToTree[120]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_025_matched__bc", &temp_UpsilonDataToTree[121]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_025_unmatched__bc", &temp_UpsilonDataToTree[122]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_025__bc", &temp_UpsilonDataToTree[123]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_025_matched__bc", &temp_UpsilonDataToTree[124]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_025_unmatched__bc", &temp_UpsilonDataToTree[125]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_050__bc", &temp_UpsilonDataToTree[126]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_050_matched__bc", &temp_UpsilonDataToTree[127]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_050_unmatched__bc", &temp_UpsilonDataToTree[128]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_050__bc", &temp_UpsilonDataToTree[129]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_050_matched__bc", &temp_UpsilonDataToTree[130]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_050_unmatched__bc", &temp_UpsilonDataToTree[131]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_075__bc", &temp_UpsilonDataToTree[132]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_075_matched__bc", &temp_UpsilonDataToTree[133]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_075_unmatched__bc", &temp_UpsilonDataToTree[134]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_075__bc", &temp_UpsilonDataToTree[135]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_075_matched__bc", &temp_UpsilonDataToTree[136]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_075_unmatched__bc", &temp_UpsilonDataToTree[137]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_100__bc", &temp_UpsilonDataToTree[138]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_100_matched__bc", &temp_UpsilonDataToTree[139]);
    temp_tree_upsilon->Branch("extraInfo__boNgammav200_100_unmatched__bc", &temp_UpsilonDataToTree[140]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_100__bc", &temp_UpsilonDataToTree[141]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_100_matched__bc", &temp_UpsilonDataToTree[142]);
    temp_tree_upsilon->Branch("extraInfo__boEeclv200_100_unmatched__bc", &temp_UpsilonDataToTree[143]);
    temp_tree_upsilon->Branch("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &temp_UpsilonDataToTree[144]);
    temp_tree_upsilon->Branch("nParticlesInList__boD0__clDecayIntoKL0__bc", &temp_UpsilonDataToTree[145]);
    temp_tree_upsilon->Branch("nParticlesInList__boB__pl__clKnn__bc", &temp_UpsilonDataToTree[146]);
    temp_tree_upsilon->Branch("invMassInLists__bon0__clKnn__bc", &temp_UpsilonDataToTree[147]);
    temp_tree_upsilon->Branch("nParticlesInList__boB__pl__clKstarnn__bc", &temp_UpsilonDataToTree[148]);
    temp_tree_upsilon->Branch("invMassInLists__bon0__clKstarnn__bc", &temp_UpsilonDataToTree[149]);
    temp_tree_upsilon->Branch("nParticlesInList__boB0__clK0nn__bc", &temp_UpsilonDataToTree[150]);
    temp_tree_upsilon->Branch("invMassInLists__bon0__clK0nn__bc", &temp_UpsilonDataToTree[151]);
    temp_tree_upsilon->Branch("nParticlesInList__boB0__clKstar0nn__bc", &temp_UpsilonDataToTree[152]);
    temp_tree_upsilon->Branch("invMassInLists__bon0__clKstar0nn__bc", &temp_UpsilonDataToTree[153]);
    temp_tree_upsilon->Branch("nParticlesInList__boB__pl__clKpKLKL_all__bc", &temp_UpsilonDataToTree[154]);
    temp_tree_upsilon->Branch("nParticlesInList__boB__pl__clKpKLKL_NR__bc", &temp_UpsilonDataToTree[155]);
    temp_tree_upsilon->Branch("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &temp_UpsilonDataToTree[156]);
    temp_tree_upsilon->Branch("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &temp_UpsilonDataToTree[157]);
    temp_tree_upsilon->Branch("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &temp_UpsilonDataToTree[158]);
    temp_tree_upsilon->Branch("nParticlesInList__boB0__clKSKLKL_all__bc", &temp_UpsilonDataToTree[159]);
    temp_tree_upsilon->Branch("nParticlesInList__boB0__clKSKLKL_NR__bc", &temp_UpsilonDataToTree[160]);
    temp_tree_upsilon->Branch("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &temp_UpsilonDataToTree[161]);
    temp_tree_upsilon->Branch("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &temp_UpsilonDataToTree[162]);
    temp_tree_upsilon->Branch("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &temp_UpsilonDataToTree[163]);

    // get Bsig_info
    temp_tree_Bsig->Branch("Bsig_E", &temp_BsigDataToTree[0]);
    temp_tree_Bsig->Branch("Bsig_useCMSFrame_E", &temp_BsigDataToTree[1]);
    temp_tree_upsilon->Branch("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc", &temp_BsigDataToTree[2]);
    temp_tree_Bsig->Branch("Bsig_p", &temp_BsigDataToTree[3]);
    temp_tree_Bsig->Branch("Bsig_useCMSFrame_p", &temp_BsigDataToTree[4]);
    temp_tree_upsilon->Branch("useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc", &temp_BsigDataToTree[5]);
    temp_tree_Bsig->Branch("Bsig_M", &temp_BsigDataToTree[6]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_KS0_3D_distance", &temp_BsigDataToTree[7]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nKslow1", &temp_BsigDataToTree[8]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nKslow2", &temp_BsigDataToTree[9]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nKslow3", &temp_BsigDataToTree[10]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nPislow1", &temp_BsigDataToTree[11]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nPislow2", &temp_BsigDataToTree[12]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nPislow3", &temp_BsigDataToTree[13]);
    temp_tree_Bsig->Branch("Bsig_R2", &temp_BsigDataToTree[14]);
    temp_tree_Bsig->Branch("Bsig_thrustBm", &temp_BsigDataToTree[15]);
    temp_tree_Bsig->Branch("Bsig_thrustOm", &temp_BsigDataToTree[16]);
    temp_tree_Bsig->Branch("Bsig_cosTBTO", &temp_BsigDataToTree[17]);
    temp_tree_Bsig->Branch("Bsig_cosTBz", &temp_BsigDataToTree[18]);
    temp_tree_Bsig->Branch("Bsig_KSFWVariables_et", &temp_BsigDataToTree[19]);
    temp_tree_Bsig->Branch("Bsig_KSFWVariables_mm2", &temp_BsigDataToTree[20]);
    temp_tree_Bsig->Branch("Bsig_KSFWVariables_hso00", &temp_BsigDataToTree[21]);
    temp_tree_Bsig->Branch("Bsig_KSFWVariables_hso01", &temp_BsigDataToTree[22]);
    temp_tree_Bsig->Branch("Bsig_KSFWVariables_hso02", &temp_BsigDataToTree[23]);
    temp_tree_Bsig->Branch("Bsig_KSFWVariables_hso03", &temp_BsigDataToTree[24]);
    temp_tree_Bsig->Branch("Bsig_KSFWVariables_hso04", &temp_BsigDataToTree[25]);
    temp_tree_Bsig->Branch("Bsig_KSFWVariables_hso10", &temp_BsigDataToTree[26]);
    temp_tree_Bsig->Branch("Bsig_KSFWVariables_hso12", &temp_BsigDataToTree[27]);
    temp_tree_Bsig->Branch("Bsig_KSFWVariables_hso14", &temp_BsigDataToTree[28]);
    temp_tree_Bsig->Branch("Bsig_KSFWVariables_hso20", &temp_BsigDataToTree[29]);
    temp_tree_Bsig->Branch("Bsig_KSFWVariables_hso22", &temp_BsigDataToTree[30]);
    temp_tree_Bsig->Branch("Bsig_KSFWVariables_hso24", &temp_BsigDataToTree[31]);
    temp_tree_Bsig->Branch("Bsig_KSFWVariables_hoo0", &temp_BsigDataToTree[32]);
    temp_tree_Bsig->Branch("Bsig_KSFWVariables_hoo1", &temp_BsigDataToTree[33]);
    temp_tree_Bsig->Branch("Bsig_KSFWVariables_hoo2", &temp_BsigDataToTree[34]);
    temp_tree_Bsig->Branch("Bsig_KSFWVariables_hoo3", &temp_BsigDataToTree[35]);
    temp_tree_Bsig->Branch("Bsig_KSFWVariables_hoo4", &temp_BsigDataToTree[36]);
    temp_tree_Bsig->Branch("Bsig_CleoConeCS_1", &temp_BsigDataToTree[37]);
    temp_tree_Bsig->Branch("Bsig_CleoConeCS_2", &temp_BsigDataToTree[38]);
    temp_tree_Bsig->Branch("Bsig_CleoConeCS_3", &temp_BsigDataToTree[39]);
    temp_tree_Bsig->Branch("Bsig_CleoConeCS_4", &temp_BsigDataToTree[40]);
    temp_tree_Bsig->Branch("Bsig_CleoConeCS_5", &temp_BsigDataToTree[41]);
    temp_tree_Bsig->Branch("Bsig_CleoConeCS_6", &temp_BsigDataToTree[42]);
    temp_tree_Bsig->Branch("Bsig_CleoConeCS_7", &temp_BsigDataToTree[43]);
    temp_tree_Bsig->Branch("Bsig_CleoConeCS_8", &temp_BsigDataToTree[44]);
    temp_tree_Bsig->Branch("Bsig_CleoConeCS_9", &temp_BsigDataToTree[45]);
    temp_tree_Bsig->Branch("Bsig_pt", &temp_BsigDataToTree[46]);
    temp_tree_Bsig->Branch("Bsig_useCMSFrame_pt", &temp_BsigDataToTree[47]);
    temp_tree_Bsig->Branch("Bsig_theta", &temp_BsigDataToTree[48]);
    temp_tree_Bsig->Branch("Bsig_useCMSFrame_theta", &temp_BsigDataToTree[49]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nDc", &temp_BsigDataToTree[50]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Dc_pValue_med", &temp_BsigDataToTree[51]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Dc_pValue_std", &temp_BsigDataToTree[52]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Dcsimpleveto_chiProb", &temp_BsigDataToTree[53]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Dcsimpleveto_dr", &temp_BsigDataToTree[54]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Dcsimpleveto_dz", &temp_BsigDataToTree[55]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Dcsimpleveto_M", &temp_BsigDataToTree[56]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nD0", &temp_BsigDataToTree[57]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_D0_pValue_med", &temp_BsigDataToTree[58]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_D0_pValue_std", &temp_BsigDataToTree[59]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_D0simpleveto_chiProb", &temp_BsigDataToTree[60]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_D0simpleveto_dr", &temp_BsigDataToTree[61]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_D0simpleveto_dz", &temp_BsigDataToTree[62]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_D0simpleveto_M", &temp_BsigDataToTree[63]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_mychiProb", &temp_BsigDataToTree[64]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_mydr", &temp_BsigDataToTree[65]);
    temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_mydz", &temp_BsigDataToTree[66]);
    for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
        temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_nKtruebin" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[67 + 4 * i_PID]);
        temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_nKmisbin" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[67 + 4 * i_PID + 1]);
        temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_npitruebin" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[67 + 4 * i_PID + 2]);
        temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_npimisbin" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[67 + 4 * i_PID + 3]);
    }
    for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_npi0bin" + std::to_string(i_pi0)).c_str(), &temp_BsigDataToTree[359 + i_pi0]);
    for (int i_PID = 0; i_PID < N_fakeE_syst; i_PID++) {
        temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_nKfakeEbin_n" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[367 + 4 * i_PID]);
        temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_nKfakeEbin_p" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[367 + 4 * i_PID + 1]);
        temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_npifakeEbin_n" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[367 + 4 * i_PID + 2]);
        temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_npifakeEbin_p" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[367 + 4 * i_PID + 3]);
    }
    for (int i_PID = 0; i_PID < N_fakeMU_syst; i_PID++) {
        temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_nKfakeMUbin_n" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[542 + 4 * i_PID]);
        temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_nKfakeMUbin_p" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[542 + 4 * i_PID + 1]);
        temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_npifakeMUbin_n" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[542 + 4 * i_PID + 2]);
        temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_npifakeMUbin_p" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[542 + 4 * i_PID + 3]);
    }

    // get Btag_info
    temp_tree_Btag->Branch("Btag_extraInfo_decayModeID", &temp_BtagDataToTree[0]);
    temp_tree_Btag->Branch("Btag_Mbc", &temp_BtagDataToTree[1]);
    temp_tree_Btag->Branch("Btag_deltaE", &temp_BtagDataToTree[2]);
    temp_tree_Btag->Branch("Btag_E", &temp_BtagDataToTree[3]);
    temp_tree_Btag->Branch("Btag_useCMSFrame_E", &temp_BtagDataToTree[4]);
    temp_tree_Btag->Branch("Btag_extraInfo_SignalProbability", &temp_BtagDataToTree[5]);
    temp_tree_Btag->Branch("Btag_chiProb", &temp_BtagDataToTree[6]);
    temp_tree_Btag->Branch("Btag_dr", &temp_BtagDataToTree[7]);
    temp_tree_Btag->Branch("Btag_dz", &temp_BtagDataToTree[8]);
    temp_tree_Btag->Branch("Btag_useCMSFrame_p", &temp_BtagDataToTree[9]);
    temp_tree_Btag->Branch("Btag_useCMSFrame_phi", &temp_BtagDataToTree[10]);

    // other information I need
    temp_tree_Btag->Branch("Btag_R2", &temp_DataToTree[0]);
    temp_tree_Btag->Branch("Btag_thrustBm", &temp_DataToTree[1]);
    temp_tree_Btag->Branch("Btag_thrustOm", &temp_DataToTree[2]);
    temp_tree_Btag->Branch("Btag_cosTBTO", &temp_DataToTree[3]);
    temp_tree_Btag->Branch("Btag_cosTBz", &temp_DataToTree[4]);
    temp_tree_Btag->Branch("Btag_KSFWVariables_et", &temp_DataToTree[5]);
    temp_tree_Btag->Branch("Btag_KSFWVariables_mm2", &temp_DataToTree[6]);
    temp_tree_Btag->Branch("Btag_KSFWVariables_hso00", &temp_DataToTree[7]);
    temp_tree_Btag->Branch("Btag_KSFWVariables_hso01", &temp_DataToTree[8]);
    temp_tree_Btag->Branch("Btag_KSFWVariables_hso02", &temp_DataToTree[9]);
    temp_tree_Btag->Branch("Btag_KSFWVariables_hso03", &temp_DataToTree[10]);
    temp_tree_Btag->Branch("Btag_KSFWVariables_hso04", &temp_DataToTree[11]);
    temp_tree_Btag->Branch("Btag_KSFWVariables_hso10", &temp_DataToTree[12]);
    temp_tree_Btag->Branch("Btag_KSFWVariables_hso12", &temp_DataToTree[13]);
    temp_tree_Btag->Branch("Btag_KSFWVariables_hso14", &temp_DataToTree[14]);
    temp_tree_Btag->Branch("Btag_KSFWVariables_hso20", &temp_DataToTree[15]);
    temp_tree_Btag->Branch("Btag_KSFWVariables_hso22", &temp_DataToTree[16]);
    temp_tree_Btag->Branch("Btag_KSFWVariables_hso24", &temp_DataToTree[17]);
    temp_tree_Btag->Branch("Btag_KSFWVariables_hoo0", &temp_DataToTree[18]);
    temp_tree_Btag->Branch("Btag_KSFWVariables_hoo1", &temp_DataToTree[19]);
    temp_tree_Btag->Branch("Btag_KSFWVariables_hoo2", &temp_DataToTree[20]);
    temp_tree_Btag->Branch("Btag_KSFWVariables_hoo3", &temp_DataToTree[21]);
    temp_tree_Btag->Branch("Btag_KSFWVariables_hoo4", &temp_DataToTree[22]);
    temp_tree_Btag->Branch("Btag_CleoConeCS_1", &temp_DataToTree[23]);
    temp_tree_Btag->Branch("Btag_CleoConeCS_2", &temp_DataToTree[24]);
    temp_tree_Btag->Branch("Btag_CleoConeCS_3", &temp_DataToTree[25]);
    temp_tree_Btag->Branch("Btag_CleoConeCS_4", &temp_DataToTree[26]);
    temp_tree_Btag->Branch("Btag_CleoConeCS_5", &temp_DataToTree[27]);
    temp_tree_Btag->Branch("Btag_CleoConeCS_6", &temp_DataToTree[28]);
    temp_tree_Btag->Branch("Btag_CleoConeCS_7", &temp_DataToTree[29]);
    temp_tree_Btag->Branch("Btag_CleoConeCS_8", &temp_DataToTree[30]);
    temp_tree_Btag->Branch("Btag_CleoConeCS_9", &temp_DataToTree[31]);
    temp_tree_upsilon->Branch("missingMass2OfEvent", &temp_DataToTree[32]);
    temp_tree_upsilon->Branch("visibleEnergyOfEventCMS", &temp_DataToTree[33]);
    temp_tree_Btag->Branch("Btag_useCMSFrame_theta", &temp_DataToTree[34]);
    temp_tree_upsilon->Branch("extraInfo__boDecayHash__bc", &temp_DataToTree[35]);
    temp_tree_upsilon->Branch("extraInfo__boDecayHashExtended__bc", &temp_DataToTree[36]);

    if (DoesItHaveXsBranch) {
        // decay mode (MC level)
        temp_tree_Xs->Branch("nParticlesInList__boB__pl__clKcharge_total__bc", &temp_DecayDataToTree[0]);
        temp_tree_Xs->Branch("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &temp_DecayDataToTree[1]);
        temp_tree_Xs->Branch("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &temp_DecayDataToTree[2]);
        temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCcomb__bc", &temp_DecayDataToTree[3]);
        temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch1__bc", &temp_DecayDataToTree[4]);
        temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch2__bc", &temp_DecayDataToTree[5]);
        temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch3__bc", &temp_DecayDataToTree[6]);
        temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch4__bc", &temp_DecayDataToTree[7]);
        temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch5__bc", &temp_DecayDataToTree[8]);
        temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch6__bc", &temp_DecayDataToTree[9]);
        temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch7__bc", &temp_DecayDataToTree[10]);
        temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch8__bc", &temp_DecayDataToTree[11]);
        temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch9__bc", &temp_DecayDataToTree[12]);
        temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch10__bc", &temp_DecayDataToTree[13]);
        temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch11__bc", &temp_DecayDataToTree[14]);
        temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch12__bc", &temp_DecayDataToTree[15]);
        temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch13__bc", &temp_DecayDataToTree[16]);
        temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch14__bc", &temp_DecayDataToTree[17]);
        temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch15__bc", &temp_DecayDataToTree[18]);
        temp_tree_Xs->Branch("nParticlesInList__boB0__clKneutral_total__bc", &temp_DecayDataToTree[19]);
        temp_tree_Xs->Branch("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &temp_DecayDataToTree[20]);
        temp_tree_Xs->Branch("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &temp_DecayDataToTree[21]);
        temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCcomb__bc", &temp_DecayDataToTree[22]);
        temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch16__bc", &temp_DecayDataToTree[23]);
        temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch17__bc", &temp_DecayDataToTree[24]);
        temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch18__bc", &temp_DecayDataToTree[25]);
        temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch19__bc", &temp_DecayDataToTree[26]);
        temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch20__bc", &temp_DecayDataToTree[27]);
        temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch21__bc", &temp_DecayDataToTree[28]);
        temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch22__bc", &temp_DecayDataToTree[29]);
        temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch23__bc", &temp_DecayDataToTree[30]);
        temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch24__bc", &temp_DecayDataToTree[31]);
        temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch25__bc", &temp_DecayDataToTree[32]);
        temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch26__bc", &temp_DecayDataToTree[33]);
        temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch27__bc", &temp_DecayDataToTree[34]);
        temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch28__bc", &temp_DecayDataToTree[35]);
        temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch29__bc", &temp_DecayDataToTree[36]);
        temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch30__bc", &temp_DecayDataToTree[37]);
        temp_tree_Xs->Branch("nParticlesInList__bonu_e__clMC_signal__bc", &temp_DecayNparticlesDataToTree[0]);
        temp_tree_Xs->Branch("nParticlesInList__boB__pl__clMC_signal_total_e__bc", &temp_DecayNparticlesDataToTree[1]);
        temp_tree_Xs->Branch("nParticlesInList__boB0__clMC_signal_total_e__bc", &temp_DecayNparticlesDataToTree[2]);
        temp_tree_Xs->Branch("nParticlesInList__boB__pl__clPrimaryMC__bc", &temp_DecayNparticlesDataToTree[3]);
        temp_tree_Xs->Branch("nParticlesInList__boB0__clPrimaryMC__bc", &temp_DecayNparticlesDataToTree[4]);
        temp_tree_Xs->Branch("invMassInLists__bonu_e__clMC_signal__bc", &temp_DecaySystFFDataToTree[index_q2]);
        temp_tree_Xs->Branch("averageValueInList__boB__pl__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", &temp_DecaySystFFDataToTree[1]);
        temp_tree_Xs->Branch("averageValueInList__boB0__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", &temp_DecaySystFFDataToTree[2]);
        temp_tree_Xs->Branch("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &temp_DecaySystFFDataToTree[3]);
        temp_tree_Xs->Branch("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &temp_DecaySystFFDataToTree[4]);
        temp_tree_Xs->Branch("averageValueInList__boB__pl__clMC_signal_total_e__cm__spM__bc", &temp_DecaySystFFDataToTree[5]);
        temp_tree_Xs->Branch("averageValueInList__boB0__clMC_signal_total_e__cm__spM__bc", &temp_DecaySystFFDataToTree[6]);
    }
    else {
        for (int i = 0; i < N_decay; i++)  temp_DecayDataToTree[i] = -1;
        for (int i = 0; i < N_decay_nparticles; i++)  temp_DecayNparticlesDataToTree[i] = -1;
        for (int i = 0; i < N_decay_syst_ff; i++)  temp_DecaySystFFDataToTree[i] = -1;
    }

    if (DoesItHaveMVAOutput) {
        temp_tree_upsilon->Branch("MVA_BB", &temp_MVA_BB_DataToTree);
        temp_tree_upsilon->Branch("MVA_Continuum", &temp_MVA_Continuum_DataToTree);
    }
    else {
        temp_MVA_BB_DataToTree = -1.0f;
        temp_MVA_Continuum_DataToTree = -1.0f;
    }

    if (DoesItHaveJpsiOutput) {
        temp_tree_upsilon->Branch("nParticlesInList__bomu__pl__clfromUpsilonmychargedMuon__bc", &temp_nROE_mu);
        temp_tree_upsilon->Branch("nParticlesInList__boJ__slpsi__clfromUpsilontemp__bc", &temp_nROE_Jpsi);
        temp_tree_upsilon->Branch("nParticlesInList__boUpsilon__bo4S__bc__cltemp_withoutneutrino_BCS__bc", &temp_nROE_Upsilon_BCS);
        temp_tree_upsilon->Branch("nParticlesInList__boUpsilon__bo4S__bc__cltemp_withoutneutrino__bc", &temp_nROE_Upsilon);
    }
    else {
        temp_nROE_mu = -1.0;
        temp_nROE_Jpsi = -1.0;
        temp_nROE_Upsilon_BCS = -1.0;
        temp_nROE_Upsilon = -1.0;
    }
    /*================================================================*/

    std::queue<Data> temp_queue;
    temp_queue.swap(TotalData);
    while (!temp_queue.empty()) {
        Data temp = temp_queue.front();
        temp_queue.pop();

        {
            temp_UpsilonExperimentToTree = temp.upsilon_experiment;
            temp_UpsilonRunToTree = temp.upsilon_run;
            temp_UpsilonEventToTree = temp.upsilon_event;
            temp_UpsilonCandidateToTree = temp.upsilon_candidate;
            temp_UpsilonNcandidatesToTree = temp.upsilon_ncandidates;
            temp_BsigExperimentToTree = temp.Bsig_experiment;
            temp_BsigRunToTree = temp.Bsig_run;
            temp_BsigEventToTree = temp.Bsig_event;
            temp_BsigCandidateToTree = temp.Bsig_candidate;
            temp_BsigNcandidatesToTree = temp.Bsig_ncandidates;
            temp_BtagExperimentToTree = temp.Btag_experiment;
            temp_BtagRunToTree = temp.Btag_run;
            temp_BtagEventToTree = temp.Btag_event;
            temp_BtagCandidateToTree = temp.Btag_candidate;
            temp_BtagNcandidatesToTree = temp.Btag_ncandidates;
        }
        for (int i = 0; i < N_Upsilon_info; i++) {
            temp_UpsilonDataToTree[i] = temp.Upsilon_info[i];
        }
        for (int i = 0; i < N_Bsig_info; i++) {
            temp_BsigDataToTree[i] = temp.Bsig_info[i];
        }
        for (int i = 0; i < N_Btag_info; i++) {
            temp_BtagDataToTree[i] = temp.Btag_info[i];
        }
        for (int i = 0; i < N_decay; i++) {
            if (DoesItHaveXsBranch) temp_DecayDataToTree[i] = temp.Decay[i];
        }
        for (int i = 0; i < N_decay_nparticles; i++) {
            if (DoesItHaveXsBranch) temp_DecayNparticlesDataToTree[i] = temp.Decay_nparticles[i];
        }
        for (int i = 0; i < N_decay_syst_ff; i++) {
            if (DoesItHaveXsBranch) temp_DecaySystFFDataToTree[i] = temp.Decay_syst_ff[i];
        }
        for (int i = 0; i < N_Needed_info; i++) {
            temp_DataToTree[i] = temp.Needed_info[i];
        }
        temp_Upsilon_decayIDToTree = temp.Upsilon_decayID;
        temp_Bsig_decayIDToTree = temp.Bsig_decayID;
        if (DoesItHaveMVAOutput) {
            temp_MVA_BB_DataToTree = temp.MVA_BB;
            temp_MVA_Continuum_DataToTree = temp.MVA_Continuum;
        }
        if (DoesItHaveJpsiOutput) {
            temp_nROE_mu = temp.nROE_mu;
            temp_nROE_Jpsi = temp.nROE_Jpsi;
            temp_nROE_Upsilon_BCS = temp.nROE_Upsilon_BCS;
            temp_nROE_Upsilon = temp.nROE_Upsilon;
        }

        temp_tree_upsilon->Fill();
        temp_tree_Bsig->Fill();
        temp_tree_Btag->Fill();
        if (DoesItHaveXsBranch) temp_tree_Xs->Fill();

        TotalData.push(temp);
    }

    temp_file->cd();
    temp_tree_upsilon->Write();
    temp_tree_Bsig->Write();
    temp_tree_Btag->Write();
    if (DoesItHaveXsBranch) temp_tree_Xs->Write();
    temp_file->Close();

}

void Loader::ConvertIntoSeparateDataFile(std::string output_name, int flag) {

    TFile* temp_file = new TFile(output_name.c_str(), "recreate");
    temp_file->cd();
    TTree* temp_tree = new TTree("data", "");

    int temp_ExperimentToTree;
    int temp_RunToTree;
    unsigned int temp_EventToTree;
    int temp_CandidateToTree;
    int temp_NcandidatesToTree;
    double temp_UpsilonDataToTree[N_Upsilon_info];
    double temp_BsigDataToTree[N_Bsig_info];
    double temp_BtagDataToTree[N_Btag_info];
    double temp_DataToTree[N_Needed_info];
    int temp_DecayDataToTree[N_decay];
    int temp_DecayNparticlesDataToTree[N_decay_nparticles];
    double temp_DecaySystFFDataToTree[N_decay_syst_ff];
    double temp_Upsilon_decayIDToTree;
    double temp_Bsig_decayIDToTree;
    int temp_flag;

    float temp_MVA_BB_DataToTree;
    float temp_MVA_Continuum_DataToTree;

    int temp_nROE_mu_DataToTree;
    int temp_nROE_Jpsi_DataToTree;
    int temp_nROE_Upsilon_BCS_DataToTree;
    int temp_nROE_Upsilon_DataToTree;

    /*================================================================*/
    // get event_info
    temp_tree->Branch("__experiment__", &temp_ExperimentToTree);
    temp_tree->Branch("__run__", &temp_RunToTree);
    temp_tree->Branch("__event__", &temp_EventToTree);
    temp_tree->Branch("__candidate__", &temp_CandidateToTree);
    temp_tree->Branch("__ncandidates__", &temp_NcandidatesToTree);

    // get decaymodeID
    temp_tree->Branch("extraInfo__bodecayModeID__bc", &temp_Upsilon_decayIDToTree);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_decayModeID", &temp_Bsig_decayIDToTree);

    // get Upsilon_info
    temp_tree->Branch("nROE_ECLClusters__bocleanMask__bc", &temp_UpsilonDataToTree[0]);
    temp_tree->Branch("nROE_KLMClusters", &temp_UpsilonDataToTree[1]);
    temp_tree->Branch("nROE_Tracks__bocleanMask__bc", &temp_UpsilonDataToTree[2]);
    temp_tree->Branch("roeEextra__bocleanMask__bc", &temp_UpsilonDataToTree[3]);
    temp_tree->Branch("nROE_NeutralECLClusters__bocleanMask__bc", &temp_UpsilonDataToTree[4]);
    temp_tree->Branch("nROE_ParticlesInList__boK_S0__clmyKaonshort__bc", &temp_UpsilonDataToTree[5]);
    temp_tree->Branch("nROE_ParticlesInList__bopi0__clmyneutralPion__bc", &temp_UpsilonDataToTree[6]);
    temp_tree->Branch("missingMomentumOfEvent_theta", &temp_UpsilonDataToTree[7]);
    temp_tree->Branch("missingMomentumOfEvent", &temp_UpsilonDataToTree[8]);
    temp_tree->Branch("missingEnergyOfEventCMS", &temp_UpsilonDataToTree[9]);
    temp_tree->Branch("nRemainingTracksInEvent", &temp_UpsilonDataToTree[10]);
    temp_tree->Branch("roeNeextra__bocleanMask__bc", &temp_UpsilonDataToTree[11]);
    temp_tree->Branch("useCMSFrame__boroeNeextra__bocleanMask__bc__bc", &temp_UpsilonDataToTree[12]);
    temp_tree->Branch("nROE_ParticlesInList__bogamma__clmygamma__bc", &temp_UpsilonDataToTree[13]);
    temp_tree->Branch("foxWolframR1", &temp_UpsilonDataToTree[14]);
    temp_tree->Branch("foxWolframR2", &temp_UpsilonDataToTree[15]);
    temp_tree->Branch("foxWolframR3", &temp_UpsilonDataToTree[16]);
    temp_tree->Branch("foxWolframR4", &temp_UpsilonDataToTree[17]);
    temp_tree->Branch("harmonicMomentThrust0", &temp_UpsilonDataToTree[18]);
    temp_tree->Branch("harmonicMomentThrust1", &temp_UpsilonDataToTree[19]);
    temp_tree->Branch("harmonicMomentThrust2", &temp_UpsilonDataToTree[20]);
    temp_tree->Branch("harmonicMomentThrust3", &temp_UpsilonDataToTree[21]);
    temp_tree->Branch("harmonicMomentThrust4", &temp_UpsilonDataToTree[22]);
    temp_tree->Branch("cleoConeThrust0", &temp_UpsilonDataToTree[23]);
    temp_tree->Branch("cleoConeThrust1", &temp_UpsilonDataToTree[24]);
    temp_tree->Branch("cleoConeThrust2", &temp_UpsilonDataToTree[25]);
    temp_tree->Branch("cleoConeThrust3", &temp_UpsilonDataToTree[26]);
    temp_tree->Branch("cleoConeThrust4", &temp_UpsilonDataToTree[27]);
    temp_tree->Branch("cleoConeThrust5", &temp_UpsilonDataToTree[28]);
    temp_tree->Branch("cleoConeThrust6", &temp_UpsilonDataToTree[29]);
    temp_tree->Branch("cleoConeThrust7", &temp_UpsilonDataToTree[30]);
    temp_tree->Branch("cleoConeThrust8", &temp_UpsilonDataToTree[31]);
    temp_tree->Branch("sphericity", &temp_UpsilonDataToTree[32]);
    temp_tree->Branch("aplanarity", &temp_UpsilonDataToTree[33]);
    temp_tree->Branch("thrust", &temp_UpsilonDataToTree[34]);
    temp_tree->Branch("thrustAxisCosTheta", &temp_UpsilonDataToTree[35]);
    temp_tree->Branch("MsquaredBsig_op0", &temp_UpsilonDataToTree[36]);
    temp_tree->Branch("MsquaredBsig_op1", &temp_UpsilonDataToTree[37]);
    temp_tree->Branch("MsquaredBsig_op2", &temp_UpsilonDataToTree[38]);
    temp_tree->Branch("MsquaredBsig_op3", &temp_UpsilonDataToTree[39]);
    temp_tree->Branch("MsquaredBsig_op4", &temp_UpsilonDataToTree[40]);
    temp_tree->Branch("MsquaredBsig_op7", &temp_UpsilonDataToTree[41]);
    temp_tree->Branch("roeP__bocleanMask__bc", &temp_UpsilonDataToTree[42]);
    temp_tree->Branch("roeM__bocleanMask__bc", &temp_UpsilonDataToTree[43]);
    temp_tree->Branch("roePTheta__bocleanMask__bc", &temp_UpsilonDataToTree[44]);
    temp_tree->Branch("qsquared", &temp_UpsilonDataToTree[45]);
    temp_tree->Branch("chiProb", &temp_UpsilonDataToTree[46]);
    temp_tree->Branch("dr", &temp_UpsilonDataToTree[47]);
    temp_tree->Branch("dz", &temp_UpsilonDataToTree[48]);
    temp_tree->Branch("nParticlesInList__boe__pl__clElectronFBDT__bc", &temp_UpsilonDataToTree[49]);
    temp_tree->Branch("nParticlesInList__bomu__pl__clMuonFBDT__bc", &temp_UpsilonDataToTree[50]);
    temp_tree->Branch("nParticlesInList__boe__pl__clElectronFBDT_loose__bc", &temp_UpsilonDataToTree[51]);
    temp_tree->Branch("nParticlesInList__bomu__pl__clMuonFBDT_loose__bc", &temp_UpsilonDataToTree[52]);
    temp_tree->Branch("nParticlesInList__boe__pl__clElectronFBDT_tight__bc", &temp_UpsilonDataToTree[53]);
    temp_tree->Branch("nParticlesInList__bomu__pl__clMuonFBDT_tight__bc", &temp_UpsilonDataToTree[54]);
    temp_tree->Branch("beamE", &temp_UpsilonDataToTree[55]);
    temp_tree->Branch("nROE_Tracks__bolooseMask__bc", &temp_UpsilonDataToTree[56]);
    temp_tree->Branch("Ecms", &temp_UpsilonDataToTree[57]);
    temp_tree->Branch("extraInfo__boNgamma__bc", &temp_UpsilonDataToTree[58]);
    temp_tree->Branch("extraInfo__boEecl__bc", &temp_UpsilonDataToTree[59]);
    temp_tree->Branch("extraInfo__boEecl_matched__bc", &temp_UpsilonDataToTree[60]);
    temp_tree->Branch("extraInfo__boEecl_unmatched__bc", &temp_UpsilonDataToTree[61]);
    temp_tree->Branch("extraInfo__boNgammav133__bc", &temp_UpsilonDataToTree[62]);
    temp_tree->Branch("extraInfo__boEeclv133__bc", &temp_UpsilonDataToTree[63]);
    temp_tree->Branch("extraInfo__boEeclv133_matched__bc", &temp_UpsilonDataToTree[64]);
    temp_tree->Branch("extraInfo__boEeclv133_unmatched__bc", &temp_UpsilonDataToTree[65]);
    temp_tree->Branch("extraInfo__boNgammav200__bc", &temp_UpsilonDataToTree[66]);
    temp_tree->Branch("extraInfo__boNgammav200_matched__bc", &temp_UpsilonDataToTree[67]);
    temp_tree->Branch("extraInfo__boNgammav200_unmatched__bc", &temp_UpsilonDataToTree[68]);
    temp_tree->Branch("extraInfo__boEeclv200__bc", &temp_UpsilonDataToTree[69]);
    temp_tree->Branch("extraInfo__boEeclv200_matched__bc", &temp_UpsilonDataToTree[70]);
    temp_tree->Branch("extraInfo__boEeclv200_unmatched__bc", &temp_UpsilonDataToTree[71]);
    temp_tree->Branch("extraInfo__boNgammav200_800__bc", &temp_UpsilonDataToTree[72]);
    temp_tree->Branch("extraInfo__boNgammav200_800_matched__bc", &temp_UpsilonDataToTree[73]);
    temp_tree->Branch("extraInfo__boNgammav200_800_unmatched__bc", &temp_UpsilonDataToTree[74]);
    temp_tree->Branch("extraInfo__boEeclv200_800__bc", &temp_UpsilonDataToTree[75]);
    temp_tree->Branch("extraInfo__boEeclv200_800_matched__bc", &temp_UpsilonDataToTree[76]);
    temp_tree->Branch("extraInfo__boEeclv200_800_unmatched__bc", &temp_UpsilonDataToTree[77]);
    temp_tree->Branch("extraInfo__boNgammav200_825__bc", &temp_UpsilonDataToTree[78]);
    temp_tree->Branch("extraInfo__boNgammav200_825_matched__bc", &temp_UpsilonDataToTree[79]);
    temp_tree->Branch("extraInfo__boNgammav200_825_unmatched__bc", &temp_UpsilonDataToTree[80]);
    temp_tree->Branch("extraInfo__boEeclv200_825__bc", &temp_UpsilonDataToTree[81]);
    temp_tree->Branch("extraInfo__boEeclv200_825_matched__bc", &temp_UpsilonDataToTree[82]);
    temp_tree->Branch("extraInfo__boEeclv200_825_unmatched__bc", &temp_UpsilonDataToTree[83]);
    temp_tree->Branch("extraInfo__boNgammav200_850__bc", &temp_UpsilonDataToTree[84]);
    temp_tree->Branch("extraInfo__boNgammav200_850_matched__bc", &temp_UpsilonDataToTree[85]);
    temp_tree->Branch("extraInfo__boNgammav200_850_unmatched__bc", &temp_UpsilonDataToTree[86]);
    temp_tree->Branch("extraInfo__boEeclv200_850__bc", &temp_UpsilonDataToTree[87]);
    temp_tree->Branch("extraInfo__boEeclv200_850_matched__bc", &temp_UpsilonDataToTree[88]);
    temp_tree->Branch("extraInfo__boEeclv200_850_unmatched__bc", &temp_UpsilonDataToTree[89]);
    temp_tree->Branch("extraInfo__boNgammav200_875__bc", &temp_UpsilonDataToTree[90]);
    temp_tree->Branch("extraInfo__boNgammav200_875_matched__bc", &temp_UpsilonDataToTree[91]);
    temp_tree->Branch("extraInfo__boNgammav200_875_unmatched__bc", &temp_UpsilonDataToTree[92]);
    temp_tree->Branch("extraInfo__boEeclv200_875__bc", &temp_UpsilonDataToTree[93]);
    temp_tree->Branch("extraInfo__boEeclv200_875_matched__bc", &temp_UpsilonDataToTree[94]);
    temp_tree->Branch("extraInfo__boEeclv200_875_unmatched__bc", &temp_UpsilonDataToTree[95]);
    temp_tree->Branch("extraInfo__boNgammav200_900__bc", &temp_UpsilonDataToTree[96]);
    temp_tree->Branch("extraInfo__boNgammav200_900_matched__bc", &temp_UpsilonDataToTree[97]);
    temp_tree->Branch("extraInfo__boNgammav200_900_unmatched__bc", &temp_UpsilonDataToTree[98]);
    temp_tree->Branch("extraInfo__boEeclv200_900__bc", &temp_UpsilonDataToTree[99]);
    temp_tree->Branch("extraInfo__boEeclv200_900_matched__bc", &temp_UpsilonDataToTree[100]);
    temp_tree->Branch("extraInfo__boEeclv200_900_unmatched__bc", &temp_UpsilonDataToTree[101]);
    temp_tree->Branch("extraInfo__boNgammav200_925__bc", &temp_UpsilonDataToTree[102]);
    temp_tree->Branch("extraInfo__boNgammav200_925_matched__bc", &temp_UpsilonDataToTree[103]);
    temp_tree->Branch("extraInfo__boNgammav200_925_unmatched__bc", &temp_UpsilonDataToTree[104]);
    temp_tree->Branch("extraInfo__boEeclv200_925__bc", &temp_UpsilonDataToTree[105]);
    temp_tree->Branch("extraInfo__boEeclv200_925_matched__bc", &temp_UpsilonDataToTree[106]);
    temp_tree->Branch("extraInfo__boEeclv200_925_unmatched__bc", &temp_UpsilonDataToTree[107]);
    temp_tree->Branch("extraInfo__boNgammav200_950__bc", &temp_UpsilonDataToTree[108]);
    temp_tree->Branch("extraInfo__boNgammav200_950_matched__bc", &temp_UpsilonDataToTree[109]);
    temp_tree->Branch("extraInfo__boNgammav200_950_unmatched__bc", &temp_UpsilonDataToTree[110]);
    temp_tree->Branch("extraInfo__boEeclv200_950__bc", &temp_UpsilonDataToTree[111]);
    temp_tree->Branch("extraInfo__boEeclv200_950_matched__bc", &temp_UpsilonDataToTree[112]);
    temp_tree->Branch("extraInfo__boEeclv200_950_unmatched__bc", &temp_UpsilonDataToTree[113]);
    temp_tree->Branch("extraInfo__boNgammav200_975__bc", &temp_UpsilonDataToTree[114]);
    temp_tree->Branch("extraInfo__boNgammav200_975_matched__bc", &temp_UpsilonDataToTree[115]);
    temp_tree->Branch("extraInfo__boNgammav200_975_unmatched__bc", &temp_UpsilonDataToTree[116]);
    temp_tree->Branch("extraInfo__boEeclv200_975__bc", &temp_UpsilonDataToTree[117]);
    temp_tree->Branch("extraInfo__boEeclv200_975_matched__bc", &temp_UpsilonDataToTree[118]);
    temp_tree->Branch("extraInfo__boEeclv200_975_unmatched__bc", &temp_UpsilonDataToTree[119]);
    temp_tree->Branch("extraInfo__boNgammav200_025__bc", &temp_UpsilonDataToTree[120]);
    temp_tree->Branch("extraInfo__boNgammav200_025_matched__bc", &temp_UpsilonDataToTree[121]);
    temp_tree->Branch("extraInfo__boNgammav200_025_unmatched__bc", &temp_UpsilonDataToTree[122]);
    temp_tree->Branch("extraInfo__boEeclv200_025__bc", &temp_UpsilonDataToTree[123]);
    temp_tree->Branch("extraInfo__boEeclv200_025_matched__bc", &temp_UpsilonDataToTree[124]);
    temp_tree->Branch("extraInfo__boEeclv200_025_unmatched__bc", &temp_UpsilonDataToTree[125]);
    temp_tree->Branch("extraInfo__boNgammav200_050__bc", &temp_UpsilonDataToTree[126]);
    temp_tree->Branch("extraInfo__boNgammav200_050_matched__bc", &temp_UpsilonDataToTree[127]);
    temp_tree->Branch("extraInfo__boNgammav200_050_unmatched__bc", &temp_UpsilonDataToTree[128]);
    temp_tree->Branch("extraInfo__boEeclv200_050__bc", &temp_UpsilonDataToTree[129]);
    temp_tree->Branch("extraInfo__boEeclv200_050_matched__bc", &temp_UpsilonDataToTree[130]);
    temp_tree->Branch("extraInfo__boEeclv200_050_unmatched__bc", &temp_UpsilonDataToTree[131]);
    temp_tree->Branch("extraInfo__boNgammav200_075__bc", &temp_UpsilonDataToTree[132]);
    temp_tree->Branch("extraInfo__boNgammav200_075_matched__bc", &temp_UpsilonDataToTree[133]);
    temp_tree->Branch("extraInfo__boNgammav200_075_unmatched__bc", &temp_UpsilonDataToTree[134]);
    temp_tree->Branch("extraInfo__boEeclv200_075__bc", &temp_UpsilonDataToTree[135]);
    temp_tree->Branch("extraInfo__boEeclv200_075_matched__bc", &temp_UpsilonDataToTree[136]);
    temp_tree->Branch("extraInfo__boEeclv200_075_unmatched__bc", &temp_UpsilonDataToTree[137]);
    temp_tree->Branch("extraInfo__boNgammav200_100__bc", &temp_UpsilonDataToTree[138]);
    temp_tree->Branch("extraInfo__boNgammav200_100_matched__bc", &temp_UpsilonDataToTree[139]);
    temp_tree->Branch("extraInfo__boNgammav200_100_unmatched__bc", &temp_UpsilonDataToTree[140]);
    temp_tree->Branch("extraInfo__boEeclv200_100__bc", &temp_UpsilonDataToTree[141]);
    temp_tree->Branch("extraInfo__boEeclv200_100_matched__bc", &temp_UpsilonDataToTree[142]);
    temp_tree->Branch("extraInfo__boEeclv200_100_unmatched__bc", &temp_UpsilonDataToTree[143]);
    temp_tree->Branch("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &temp_UpsilonDataToTree[144]);
    temp_tree->Branch("nParticlesInList__boD0__clDecayIntoKL0__bc", &temp_UpsilonDataToTree[145]);
    temp_tree->Branch("nParticlesInList__boB__pl__clKnn__bc", &temp_UpsilonDataToTree[146]);
    temp_tree->Branch("invMassInLists__bon0__clKnn__bc", &temp_UpsilonDataToTree[147]);
    temp_tree->Branch("nParticlesInList__boB__pl__clKstarnn__bc", &temp_UpsilonDataToTree[148]);
    temp_tree->Branch("invMassInLists__bon0__clKstarnn__bc", &temp_UpsilonDataToTree[149]);
    temp_tree->Branch("nParticlesInList__boB0__clK0nn__bc", &temp_UpsilonDataToTree[150]);
    temp_tree->Branch("invMassInLists__bon0__clK0nn__bc", &temp_UpsilonDataToTree[151]);
    temp_tree->Branch("nParticlesInList__boB0__clKstar0nn__bc", &temp_UpsilonDataToTree[152]);
    temp_tree->Branch("invMassInLists__bon0__clKstar0nn__bc", &temp_UpsilonDataToTree[153]);
    temp_tree->Branch("nParticlesInList__boB__pl__clKpKLKL_all__bc", &temp_UpsilonDataToTree[154]);
    temp_tree->Branch("nParticlesInList__boB__pl__clKpKLKL_NR__bc", &temp_UpsilonDataToTree[155]);
    temp_tree->Branch("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &temp_UpsilonDataToTree[156]);
    temp_tree->Branch("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &temp_UpsilonDataToTree[157]);
    temp_tree->Branch("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &temp_UpsilonDataToTree[158]);
    temp_tree->Branch("nParticlesInList__boB0__clKSKLKL_all__bc", &temp_UpsilonDataToTree[159]);
    temp_tree->Branch("nParticlesInList__boB0__clKSKLKL_NR__bc", &temp_UpsilonDataToTree[160]);
    temp_tree->Branch("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &temp_UpsilonDataToTree[161]);
    temp_tree->Branch("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &temp_UpsilonDataToTree[162]);
    temp_tree->Branch("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &temp_UpsilonDataToTree[163]);

    // get Bsig_info
    temp_tree->Branch("Bsig_E", &temp_BsigDataToTree[0]);
    temp_tree->Branch("Bsig_useCMSFrame_E", &temp_BsigDataToTree[1]);
    temp_tree->Branch("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc", &temp_BsigDataToTree[2]);
    temp_tree->Branch("Bsig_p", &temp_BsigDataToTree[3]);
    temp_tree->Branch("Bsig_useCMSFrame_p", &temp_BsigDataToTree[4]);
    temp_tree->Branch("useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc", &temp_BsigDataToTree[5]);
    temp_tree->Branch("Bsig_M", &temp_BsigDataToTree[6]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_KS0_3D_distance", &temp_BsigDataToTree[7]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_nKslow1", &temp_BsigDataToTree[8]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_nKslow2", &temp_BsigDataToTree[9]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_nKslow3", &temp_BsigDataToTree[10]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_nPislow1", &temp_BsigDataToTree[11]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_nPislow2", &temp_BsigDataToTree[12]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_nPislow3", &temp_BsigDataToTree[13]);
    temp_tree->Branch("Bsig_R2", &temp_BsigDataToTree[14]);
    temp_tree->Branch("Bsig_thrustBm", &temp_BsigDataToTree[15]);
    temp_tree->Branch("Bsig_thrustOm", &temp_BsigDataToTree[16]);
    temp_tree->Branch("Bsig_cosTBTO", &temp_BsigDataToTree[17]);
    temp_tree->Branch("Bsig_cosTBz", &temp_BsigDataToTree[18]);
    temp_tree->Branch("Bsig_KSFWVariables_et", &temp_BsigDataToTree[19]);
    temp_tree->Branch("Bsig_KSFWVariables_mm2", &temp_BsigDataToTree[20]);
    temp_tree->Branch("Bsig_KSFWVariables_hso00", &temp_BsigDataToTree[21]);
    temp_tree->Branch("Bsig_KSFWVariables_hso01", &temp_BsigDataToTree[22]);
    temp_tree->Branch("Bsig_KSFWVariables_hso02", &temp_BsigDataToTree[23]);
    temp_tree->Branch("Bsig_KSFWVariables_hso03", &temp_BsigDataToTree[24]);
    temp_tree->Branch("Bsig_KSFWVariables_hso04", &temp_BsigDataToTree[25]);
    temp_tree->Branch("Bsig_KSFWVariables_hso10", &temp_BsigDataToTree[26]);
    temp_tree->Branch("Bsig_KSFWVariables_hso12", &temp_BsigDataToTree[27]);
    temp_tree->Branch("Bsig_KSFWVariables_hso14", &temp_BsigDataToTree[28]);
    temp_tree->Branch("Bsig_KSFWVariables_hso20", &temp_BsigDataToTree[29]);
    temp_tree->Branch("Bsig_KSFWVariables_hso22", &temp_BsigDataToTree[30]);
    temp_tree->Branch("Bsig_KSFWVariables_hso24", &temp_BsigDataToTree[31]);
    temp_tree->Branch("Bsig_KSFWVariables_hoo0", &temp_BsigDataToTree[32]);
    temp_tree->Branch("Bsig_KSFWVariables_hoo1", &temp_BsigDataToTree[33]);
    temp_tree->Branch("Bsig_KSFWVariables_hoo2", &temp_BsigDataToTree[34]);
    temp_tree->Branch("Bsig_KSFWVariables_hoo3", &temp_BsigDataToTree[35]);
    temp_tree->Branch("Bsig_KSFWVariables_hoo4", &temp_BsigDataToTree[36]);
    temp_tree->Branch("Bsig_CleoConeCS_1", &temp_BsigDataToTree[37]);
    temp_tree->Branch("Bsig_CleoConeCS_2", &temp_BsigDataToTree[38]);
    temp_tree->Branch("Bsig_CleoConeCS_3", &temp_BsigDataToTree[39]);
    temp_tree->Branch("Bsig_CleoConeCS_4", &temp_BsigDataToTree[40]);
    temp_tree->Branch("Bsig_CleoConeCS_5", &temp_BsigDataToTree[41]);
    temp_tree->Branch("Bsig_CleoConeCS_6", &temp_BsigDataToTree[42]);
    temp_tree->Branch("Bsig_CleoConeCS_7", &temp_BsigDataToTree[43]);
    temp_tree->Branch("Bsig_CleoConeCS_8", &temp_BsigDataToTree[44]);
    temp_tree->Branch("Bsig_CleoConeCS_9", &temp_BsigDataToTree[45]);
    temp_tree->Branch("Bsig_pt", &temp_BsigDataToTree[46]);
    temp_tree->Branch("Bsig_useCMSFrame_pt", &temp_BsigDataToTree[47]);
    temp_tree->Branch("Bsig_theta", &temp_BsigDataToTree[48]);
    temp_tree->Branch("Bsig_useCMSFrame_theta", &temp_BsigDataToTree[49]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_nDc", &temp_BsigDataToTree[50]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_Dc_pValue_med", &temp_BsigDataToTree[51]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_Dc_pValue_std", &temp_BsigDataToTree[52]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_Dcsimpleveto_chiProb", &temp_BsigDataToTree[53]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_Dcsimpleveto_dr", &temp_BsigDataToTree[54]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_Dcsimpleveto_dz", &temp_BsigDataToTree[55]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_Dcsimpleveto_M", &temp_BsigDataToTree[56]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_nD0", &temp_BsigDataToTree[57]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_D0_pValue_med", &temp_BsigDataToTree[58]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_D0_pValue_std", &temp_BsigDataToTree[59]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_D0simpleveto_chiProb", &temp_BsigDataToTree[60]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_D0simpleveto_dr", &temp_BsigDataToTree[61]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_D0simpleveto_dz", &temp_BsigDataToTree[62]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_D0simpleveto_M", &temp_BsigDataToTree[63]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_mychiProb", &temp_BsigDataToTree[64]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_mydr", &temp_BsigDataToTree[65]);
    temp_tree->Branch("Bsig_daughter_0_extraInfo_mydz", &temp_BsigDataToTree[66]);
    for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
        temp_tree->Branch(("Bsig_daughter_0_extraInfo_nKtruebin" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[67 + 4 * i_PID]);
        temp_tree->Branch(("Bsig_daughter_0_extraInfo_nKmisbin" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[67 + 4 * i_PID + 1]);
        temp_tree->Branch(("Bsig_daughter_0_extraInfo_npitruebin" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[67 + 4 * i_PID + 2]);
        temp_tree->Branch(("Bsig_daughter_0_extraInfo_npimisbin" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[67 + 4 * i_PID + 3]);
    }
    for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) temp_tree->Branch(("Bsig_daughter_0_extraInfo_npi0bin" + std::to_string(i_pi0)).c_str(), &temp_BsigDataToTree[359 + i_pi0]);
    for (int i_PID = 0; i_PID < N_fakeE_syst; i_PID++) {
        temp_tree->Branch(("Bsig_daughter_0_extraInfo_nKfakeEbin_n" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[367 + 4 * i_PID]);
        temp_tree->Branch(("Bsig_daughter_0_extraInfo_nKfakeEbin_p" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[367 + 4 * i_PID + 1]);
        temp_tree->Branch(("Bsig_daughter_0_extraInfo_npifakeEbin_n" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[367 + 4 * i_PID + 2]);
        temp_tree->Branch(("Bsig_daughter_0_extraInfo_npifakeEbin_p" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[367 + 4 * i_PID + 3]);
    }
    for (int i_PID = 0; i_PID < N_fakeMU_syst; i_PID++) {
        temp_tree->Branch(("Bsig_daughter_0_extraInfo_nKfakeMUbin_n" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[542 + 4 * i_PID]);
        temp_tree->Branch(("Bsig_daughter_0_extraInfo_nKfakeMUbin_p" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[542 + 4 * i_PID + 1]);
        temp_tree->Branch(("Bsig_daughter_0_extraInfo_npifakeMUbin_n" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[542 + 4 * i_PID + 2]);
        temp_tree->Branch(("Bsig_daughter_0_extraInfo_npifakeMUbin_p" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[542 + 4 * i_PID + 3]);
    }


    // get Btag_info
    temp_tree->Branch("Btag_extraInfo_decayModeID", &temp_BtagDataToTree[0]);
    temp_tree->Branch("Btag_Mbc", &temp_BtagDataToTree[1]);
    temp_tree->Branch("Btag_deltaE", &temp_BtagDataToTree[2]);
    temp_tree->Branch("Btag_E", &temp_BtagDataToTree[3]);
    temp_tree->Branch("Btag_useCMSFrame_E", &temp_BtagDataToTree[4]);
    temp_tree->Branch("Btag_extraInfo_SignalProbability", &temp_BtagDataToTree[5]);
    temp_tree->Branch("Btag_chiProb", &temp_BtagDataToTree[6]);
    temp_tree->Branch("Btag_dr", &temp_BtagDataToTree[7]);
    temp_tree->Branch("Btag_dz", &temp_BtagDataToTree[8]);
    temp_tree->Branch("Btag_useCMSFrame_p", &temp_BtagDataToTree[9]);
    temp_tree->Branch("Btag_useCMSFrame_phi", &temp_BtagDataToTree[10]);

    // other information I need
    temp_tree->Branch("Btag_R2", &temp_DataToTree[0]);
    temp_tree->Branch("Btag_thrustBm", &temp_DataToTree[1]);
    temp_tree->Branch("Btag_thrustOm", &temp_DataToTree[2]);
    temp_tree->Branch("Btag_cosTBTO", &temp_DataToTree[3]);
    temp_tree->Branch("Btag_cosTBz", &temp_DataToTree[4]);
    temp_tree->Branch("Btag_KSFWVariables_et", &temp_DataToTree[5]);
    temp_tree->Branch("Btag_KSFWVariables_mm2", &temp_DataToTree[6]);
    temp_tree->Branch("Btag_KSFWVariables_hso00", &temp_DataToTree[7]);
    temp_tree->Branch("Btag_KSFWVariables_hso01", &temp_DataToTree[8]);
    temp_tree->Branch("Btag_KSFWVariables_hso02", &temp_DataToTree[9]);
    temp_tree->Branch("Btag_KSFWVariables_hso03", &temp_DataToTree[10]);
    temp_tree->Branch("Btag_KSFWVariables_hso04", &temp_DataToTree[11]);
    temp_tree->Branch("Btag_KSFWVariables_hso10", &temp_DataToTree[12]);
    temp_tree->Branch("Btag_KSFWVariables_hso12", &temp_DataToTree[13]);
    temp_tree->Branch("Btag_KSFWVariables_hso14", &temp_DataToTree[14]);
    temp_tree->Branch("Btag_KSFWVariables_hso20", &temp_DataToTree[15]);
    temp_tree->Branch("Btag_KSFWVariables_hso22", &temp_DataToTree[16]);
    temp_tree->Branch("Btag_KSFWVariables_hso24", &temp_DataToTree[17]);
    temp_tree->Branch("Btag_KSFWVariables_hoo0", &temp_DataToTree[18]);
    temp_tree->Branch("Btag_KSFWVariables_hoo1", &temp_DataToTree[19]);
    temp_tree->Branch("Btag_KSFWVariables_hoo2", &temp_DataToTree[20]);
    temp_tree->Branch("Btag_KSFWVariables_hoo3", &temp_DataToTree[21]);
    temp_tree->Branch("Btag_KSFWVariables_hoo4", &temp_DataToTree[22]);
    temp_tree->Branch("Btag_CleoConeCS_1", &temp_DataToTree[23]);
    temp_tree->Branch("Btag_CleoConeCS_2", &temp_DataToTree[24]);
    temp_tree->Branch("Btag_CleoConeCS_3", &temp_DataToTree[25]);
    temp_tree->Branch("Btag_CleoConeCS_4", &temp_DataToTree[26]);
    temp_tree->Branch("Btag_CleoConeCS_5", &temp_DataToTree[27]);
    temp_tree->Branch("Btag_CleoConeCS_6", &temp_DataToTree[28]);
    temp_tree->Branch("Btag_CleoConeCS_7", &temp_DataToTree[29]);
    temp_tree->Branch("Btag_CleoConeCS_8", &temp_DataToTree[30]);
    temp_tree->Branch("Btag_CleoConeCS_9", &temp_DataToTree[31]);
    temp_tree->Branch("missingMass2OfEvent", &temp_DataToTree[32]);
    temp_tree->Branch("visibleEnergyOfEventCMS", &temp_DataToTree[33]);
    temp_tree->Branch("Btag_useCMSFrame_theta", &temp_DataToTree[34]);
    temp_tree->Branch("extraInfo__boDecayHash__bc", &temp_DataToTree[35]);
    temp_tree->Branch("extraInfo__boDecayHashExtended__bc", &temp_DataToTree[36]);

    if (DoesItHaveXsBranch) {
        // decay mode (MC level)
        temp_tree->Branch("nParticlesInList__boB__pl__clKcharge_total__bc", &temp_DecayDataToTree[0]);
        temp_tree->Branch("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &temp_DecayDataToTree[1]);
        temp_tree->Branch("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &temp_DecayDataToTree[2]);
        temp_tree->Branch("nParticlesInList__boXsu__clMCcomb__bc", &temp_DecayDataToTree[3]);
        temp_tree->Branch("nParticlesInList__boXsu__clMCch1__bc", &temp_DecayDataToTree[4]);
        temp_tree->Branch("nParticlesInList__boXsu__clMCch2__bc", &temp_DecayDataToTree[5]);
        temp_tree->Branch("nParticlesInList__boXsu__clMCch3__bc", &temp_DecayDataToTree[6]);
        temp_tree->Branch("nParticlesInList__boXsu__clMCch4__bc", &temp_DecayDataToTree[7]);
        temp_tree->Branch("nParticlesInList__boXsu__clMCch5__bc", &temp_DecayDataToTree[8]);
        temp_tree->Branch("nParticlesInList__boXsu__clMCch6__bc", &temp_DecayDataToTree[9]);
        temp_tree->Branch("nParticlesInList__boXsu__clMCch7__bc", &temp_DecayDataToTree[10]);
        temp_tree->Branch("nParticlesInList__boXsu__clMCch8__bc", &temp_DecayDataToTree[11]);
        temp_tree->Branch("nParticlesInList__boXsu__clMCch9__bc", &temp_DecayDataToTree[12]);
        temp_tree->Branch("nParticlesInList__boXsu__clMCch10__bc", &temp_DecayDataToTree[13]);
        temp_tree->Branch("nParticlesInList__boXsu__clMCch11__bc", &temp_DecayDataToTree[14]);
        temp_tree->Branch("nParticlesInList__boXsu__clMCch12__bc", &temp_DecayDataToTree[15]);
        temp_tree->Branch("nParticlesInList__boXsu__clMCch13__bc", &temp_DecayDataToTree[16]);
        temp_tree->Branch("nParticlesInList__boXsu__clMCch14__bc", &temp_DecayDataToTree[17]);
        temp_tree->Branch("nParticlesInList__boXsu__clMCch15__bc", &temp_DecayDataToTree[18]);
        temp_tree->Branch("nParticlesInList__boB0__clKneutral_total__bc", &temp_DecayDataToTree[19]);
        temp_tree->Branch("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &temp_DecayDataToTree[20]);
        temp_tree->Branch("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &temp_DecayDataToTree[21]);
        temp_tree->Branch("nParticlesInList__boXsd__clMCcomb__bc", &temp_DecayDataToTree[22]);
        temp_tree->Branch("nParticlesInList__boXsd__clMCch16__bc", &temp_DecayDataToTree[23]);
        temp_tree->Branch("nParticlesInList__boXsd__clMCch17__bc", &temp_DecayDataToTree[24]);
        temp_tree->Branch("nParticlesInList__boXsd__clMCch18__bc", &temp_DecayDataToTree[25]);
        temp_tree->Branch("nParticlesInList__boXsd__clMCch19__bc", &temp_DecayDataToTree[26]);
        temp_tree->Branch("nParticlesInList__boXsd__clMCch20__bc", &temp_DecayDataToTree[27]);
        temp_tree->Branch("nParticlesInList__boXsd__clMCch21__bc", &temp_DecayDataToTree[28]);
        temp_tree->Branch("nParticlesInList__boXsd__clMCch22__bc", &temp_DecayDataToTree[29]);
        temp_tree->Branch("nParticlesInList__boXsd__clMCch23__bc", &temp_DecayDataToTree[30]);
        temp_tree->Branch("nParticlesInList__boXsd__clMCch24__bc", &temp_DecayDataToTree[31]);
        temp_tree->Branch("nParticlesInList__boXsd__clMCch25__bc", &temp_DecayDataToTree[32]);
        temp_tree->Branch("nParticlesInList__boXsd__clMCch26__bc", &temp_DecayDataToTree[33]);
        temp_tree->Branch("nParticlesInList__boXsd__clMCch27__bc", &temp_DecayDataToTree[34]);
        temp_tree->Branch("nParticlesInList__boXsd__clMCch28__bc", &temp_DecayDataToTree[35]);
        temp_tree->Branch("nParticlesInList__boXsd__clMCch29__bc", &temp_DecayDataToTree[36]);
        temp_tree->Branch("nParticlesInList__boXsd__clMCch30__bc", &temp_DecayDataToTree[37]);
        temp_tree->Branch("nParticlesInList__bonu_e__clMC_signal__bc", &temp_DecayNparticlesDataToTree[0]);
        temp_tree->Branch("nParticlesInList__boB__pl__clMC_signal_total_e__bc", &temp_DecayNparticlesDataToTree[1]);
        temp_tree->Branch("nParticlesInList__boB0__clMC_signal_total_e__bc", &temp_DecayNparticlesDataToTree[2]);
        temp_tree->Branch("nParticlesInList__boB__pl__clPrimaryMC__bc", &temp_DecayNparticlesDataToTree[3]);
        temp_tree->Branch("nParticlesInList__boB0__clPrimaryMC__bc", &temp_DecayNparticlesDataToTree[4]);
        temp_tree->Branch("invMassInLists__bonu_e__clMC_signal__bc", &temp_DecaySystFFDataToTree[index_q2]);
        temp_tree->Branch("averageValueInList__boB__pl__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", &temp_DecaySystFFDataToTree[1]);
        temp_tree->Branch("averageValueInList__boB0__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", &temp_DecaySystFFDataToTree[2]);
        temp_tree->Branch("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &temp_DecaySystFFDataToTree[3]);
        temp_tree->Branch("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &temp_DecaySystFFDataToTree[4]);
        temp_tree->Branch("averageValueInList__boB__pl__clMC_signal_total_e__cm__spM__bc", &temp_DecaySystFFDataToTree[5]);
        temp_tree->Branch("averageValueInList__boB0__clMC_signal_total_e__cm__spM__bc", &temp_DecaySystFFDataToTree[6]);
    }
    else {
        for (int i = 0; i < N_decay; i++)  temp_DecayDataToTree[i] = -1;
    }

    if (DoesItHaveMVAOutput) {
        temp_tree->Branch("MVA_BB", &temp_MVA_BB_DataToTree);
        temp_tree->Branch("MVA_Continuum", &temp_MVA_Continuum_DataToTree);
    }
    else {
        temp_MVA_BB_DataToTree = -1.0f;
        temp_MVA_Continuum_DataToTree = -1.0f;
    }

    if (DoesItHaveJpsiOutput) {
        temp_tree->Branch("nParticlesInList__bomu__pl__clfromUpsilonmychargedMuon__bc", &temp_nROE_mu_DataToTree);
        temp_tree->Branch("nParticlesInList__boJ__slpsi__clfromUpsilontemp__bc", &temp_nROE_Jpsi_DataToTree);
        temp_tree->Branch("nParticlesInList__boUpsilon__bo4S__bc__cltemp_withoutneutrino_BCS__bc", &temp_nROE_Upsilon_BCS_DataToTree);
        temp_tree->Branch("nParticlesInList__boUpsilon__bo4S__bc__cltemp_withoutneutrino__bc", &temp_nROE_Upsilon_DataToTree);
    }
    else {
        temp_nROE_mu_DataToTree = -1.0;
        temp_nROE_Jpsi_DataToTree = -1.0;
        temp_nROE_Upsilon_BCS_DataToTree = -1.0;
        temp_nROE_Upsilon_DataToTree = -1.0;
    }

    // flag
    temp_tree->Branch("flag", &temp_flag);
    /*================================================================*/

    std::queue<Data> temp_queue;
    temp_queue.swap(TotalData);
    while (!temp_queue.empty()) {
        Data temp = temp_queue.front();
        temp_queue.pop();

        {
            temp_ExperimentToTree = temp.__experiment__;
            temp_RunToTree = temp.__run__;
            temp_EventToTree = temp.__event__;
            temp_CandidateToTree = temp.__candidate__;
            temp_NcandidatesToTree = temp.__ncandidates__;
        }
        for (int i = 0; i < N_Upsilon_info; i++) {
            temp_UpsilonDataToTree[i] = temp.Upsilon_info[i];
        }
        for (int i = 0; i < N_Bsig_info; i++) {
            temp_BsigDataToTree[i] = temp.Bsig_info[i];
        }
        for (int i = 0; i < N_Btag_info; i++) {
            temp_BtagDataToTree[i] = temp.Btag_info[i];
        }
        for (int i = 0; i < N_Needed_info; i++) {
            temp_DataToTree[i] = temp.Needed_info[i];
        }
        for (int i = 0; i < N_decay; i++) {
            if (DoesItHaveXsBranch) temp_DecayDataToTree[i] = temp.Decay[i];
        }
        for (int i = 0; i < N_decay_nparticles; i++) {
            if (DoesItHaveXsBranch) temp_DecayNparticlesDataToTree[i] = temp.Decay_nparticles[i];
        }
        for (int i = 0; i < N_decay_syst_ff; i++) {
            if (DoesItHaveXsBranch) temp_DecaySystFFDataToTree[i] = temp.Decay_syst_ff[i];
        }
        temp_Upsilon_decayIDToTree = temp.Upsilon_decayID;
        temp_Bsig_decayIDToTree = temp.Bsig_decayID;

        if (DoesItHaveMVAOutput) {
            temp_MVA_BB_DataToTree = temp.MVA_BB;
            temp_MVA_Continuum_DataToTree = temp.MVA_Continuum;
        }

        if (DoesItHaveJpsiOutput) {
            temp_nROE_mu_DataToTree = temp.nROE_mu;
            temp_nROE_Jpsi_DataToTree = temp.nROE_Jpsi;
            temp_nROE_Upsilon_BCS_DataToTree = temp.nROE_Upsilon_BCS;
            temp_nROE_Upsilon_DataToTree = temp.nROE_Upsilon;
        }

        temp_flag = flag;

        temp_tree->Fill();

        TotalData.push(temp);
    }

    temp_file->cd();
    temp_tree->Write();
    temp_file->Close();

}

bool Loader::TrueIfDecayModeMatch(Data temp_data, Loader::DecayMode decaymode) {
    switch (decaymode) {
    case Loader::B2Kc:
        if (temp_data.Upsilon_decayID > -0.5 && temp_data.Upsilon_decayID < 0.5 && temp_data.Bsig_decayID > -0.5 && temp_data.Bsig_decayID < 0.5) return true;
        return false;
        break;
    case Loader::B2KcPi0:
        if (temp_data.Upsilon_decayID > -0.5 && temp_data.Upsilon_decayID < 0.5 && temp_data.Bsig_decayID > 0.5 && temp_data.Bsig_decayID < 1.5) return true;
        return false;
        break;
    case Loader::B2Ks0Pic:
        if (temp_data.Upsilon_decayID > -0.5 && temp_data.Upsilon_decayID < 0.5 && temp_data.Bsig_decayID > 1.5 && temp_data.Bsig_decayID < 2.5) return true;
        return false;
        break;
    case Loader::B2KcPicPic:
        if (temp_data.Upsilon_decayID > -0.5 && temp_data.Upsilon_decayID < 0.5 && temp_data.Bsig_decayID > 2.5 && temp_data.Bsig_decayID < 3.5) return true;
        return false;
        break;
    case Loader::B2Ks0PicPi0:
        if (temp_data.Upsilon_decayID > -0.5 && temp_data.Upsilon_decayID < 0.5 && temp_data.Bsig_decayID > 3.5 && temp_data.Bsig_decayID < 4.5) return true;
        return false;
        break;
    case Loader::B2KcPicPicPi0:
        if (temp_data.Upsilon_decayID > -0.5 && temp_data.Upsilon_decayID < 0.5 && temp_data.Bsig_decayID > 4.5 && temp_data.Bsig_decayID < 5.5) return true;
        return false;
        break;
    case Loader::B2Ks0PicPicPic:
        if (temp_data.Upsilon_decayID > -0.5 && temp_data.Upsilon_decayID < 0.5 && temp_data.Bsig_decayID > 5.5 && temp_data.Bsig_decayID < 6.5) return true;
        return false;
        break;
    case Loader::B2KcPicPicPicPic:
        if (temp_data.Upsilon_decayID > -0.5 && temp_data.Upsilon_decayID < 0.5 && temp_data.Bsig_decayID > 6.5 && temp_data.Bsig_decayID < 7.5) return true;
        return false;
        break;
    case Loader::B2Ks0PicPicPicPi0:
        if (temp_data.Upsilon_decayID > -0.5 && temp_data.Upsilon_decayID < 0.5 && temp_data.Bsig_decayID > 7.5 && temp_data.Bsig_decayID < 8.5) return true;
        return false;
        break;
    case Loader::B2KcPi0Pi0:
        if (temp_data.Upsilon_decayID > -0.5 && temp_data.Upsilon_decayID < 0.5 && temp_data.Bsig_decayID > 8.5 && temp_data.Bsig_decayID < 9.5) return true;
        return false;
        break;
    case Loader::B2Ks0PicPi0Pi0:
        if (temp_data.Upsilon_decayID > -0.5 && temp_data.Upsilon_decayID < 0.5 && temp_data.Bsig_decayID > 9.5 && temp_data.Bsig_decayID < 10.5) return true;
        return false;
        break;
    case Loader::B2KcPicPicPi0Pi0:
        if (temp_data.Upsilon_decayID > -0.5 && temp_data.Upsilon_decayID < 0.5 && temp_data.Bsig_decayID > 10.5 && temp_data.Bsig_decayID < 11.5) return true;
        return false;
        break;
    case Loader::B2KcKcKc:
        if (temp_data.Upsilon_decayID > -0.5 && temp_data.Upsilon_decayID < 0.5 && temp_data.Bsig_decayID > 11.5 && temp_data.Bsig_decayID < 12.5) return true;
        return false;
        break;
    case Loader::B2KcKcKs0Pic:
        if (temp_data.Upsilon_decayID > -0.5 && temp_data.Upsilon_decayID < 0.5 && temp_data.Bsig_decayID > 12.5 && temp_data.Bsig_decayID < 13.5) return true;
        return false;
        break;
    case Loader::B2KcKcKcPi0:
        if (temp_data.Upsilon_decayID > -0.5 && temp_data.Upsilon_decayID < 0.5 && temp_data.Bsig_decayID > 13.5 && temp_data.Bsig_decayID < 14.5) return true;
        return false;
        break;
    case Loader::B02Ks0:
        if (temp_data.Upsilon_decayID > 0.5 && temp_data.Upsilon_decayID < 1.5 && temp_data.Bsig_decayID > -0.5 && temp_data.Bsig_decayID < 0.5) return true;
        return false;
        break;
    case Loader::B02KcPic:
        if (temp_data.Upsilon_decayID > 0.5 && temp_data.Upsilon_decayID < 1.5 && temp_data.Bsig_decayID > 0.5 && temp_data.Bsig_decayID < 1.5) return true;
        return false;
        break;
    case Loader::B02Ks0Pi0:
        if (temp_data.Upsilon_decayID > 0.5 && temp_data.Upsilon_decayID < 1.5 && temp_data.Bsig_decayID > 1.5 && temp_data.Bsig_decayID < 2.5) return true;
        return false;
        break;
    case Loader::B02KcPicPi0:
        if (temp_data.Upsilon_decayID > 0.5 && temp_data.Upsilon_decayID < 1.5 && temp_data.Bsig_decayID > 2.5 && temp_data.Bsig_decayID < 3.5) return true;
        return false;
        break;
    case Loader::B02Ks0PicPic:
        if (temp_data.Upsilon_decayID > 0.5 && temp_data.Upsilon_decayID < 1.5 && temp_data.Bsig_decayID > 3.5 && temp_data.Bsig_decayID < 4.5) return true;
        return false;
        break;
    case Loader::B02KcPicPicPic:
        if (temp_data.Upsilon_decayID > 0.5 && temp_data.Upsilon_decayID < 1.5 && temp_data.Bsig_decayID > 4.5 && temp_data.Bsig_decayID < 5.5) return true;
        return false;
        break;
    case Loader::B02Ks0PicPicPi0:
        if (temp_data.Upsilon_decayID > 0.5 && temp_data.Upsilon_decayID < 1.5 && temp_data.Bsig_decayID > 5.5 && temp_data.Bsig_decayID < 6.5) return true;
        return false;
        break;
    case Loader::B02KcPicPicPicPi0:
        if (temp_data.Upsilon_decayID > 0.5 && temp_data.Upsilon_decayID < 1.5 && temp_data.Bsig_decayID > 6.5 && temp_data.Bsig_decayID < 7.5) return true;
        return false;
        break;
    case Loader::B02Ks0PicPicPicPic:
        if (temp_data.Upsilon_decayID > 0.5 && temp_data.Upsilon_decayID < 1.5 && temp_data.Bsig_decayID > 7.5 && temp_data.Bsig_decayID < 8.5) return true;
        return false;
        break;
    case Loader::B02Ks0Pi0Pi0:
        if (temp_data.Upsilon_decayID > 0.5 && temp_data.Upsilon_decayID < 1.5 && temp_data.Bsig_decayID > 8.5 && temp_data.Bsig_decayID < 9.5) return true;
        return false;
        break;
    case Loader::B02KcPicPi0Pi0:
        if (temp_data.Upsilon_decayID > 0.5 && temp_data.Upsilon_decayID < 1.5 && temp_data.Bsig_decayID > 9.5 && temp_data.Bsig_decayID < 10.5) return true;
        return false;
        break;
    case Loader::B02Ks0PicPicPi0Pi0:
        if (temp_data.Upsilon_decayID > 0.5 && temp_data.Upsilon_decayID < 1.5 && temp_data.Bsig_decayID > 10.5 && temp_data.Bsig_decayID < 11.5) return true;
        return false;
        break;
    case Loader::B02KcKcKs0:
        if (temp_data.Upsilon_decayID > 0.5 && temp_data.Upsilon_decayID < 1.5 && temp_data.Bsig_decayID > 11.5 && temp_data.Bsig_decayID < 12.5) return true;
        return false;
        break;
    case Loader::B02KcKcKcPic:
        if (temp_data.Upsilon_decayID > 0.5 && temp_data.Upsilon_decayID < 1.5 && temp_data.Bsig_decayID > 12.5 && temp_data.Bsig_decayID < 13.5) return true;
        return false;
        break;
    case Loader::B02KcKcKs0Pi0:
        if (temp_data.Upsilon_decayID > 0.5 && temp_data.Upsilon_decayID < 1.5 && temp_data.Bsig_decayID > 13.5 && temp_data.Bsig_decayID < 14.5) return true;
        return false;
        break;
    default:
        printf("ERROR! 036\n");
        exit(1);
        break;
    }

    printf("ERROR! 037\n");
    exit(1);
    return false;
}

bool Loader::TrueIfDecayModeMatch_MC(Data temp_data, Loader::DecayModeMC decaymodeMC) {
    if (DoesItHaveXsBranch == false) {
        printf("ERROR! TrueIfDecayModeMatch_MC is called when the data does not have Xs branch\n");
        exit(1);
    }

    switch (decaymodeMC) {
    case Loader::Xsu2Kc_MC:
        if (temp_data.Decay[0] > 0) return true;
        return false;
        break;
    case Loader::Xsu2Kcstar2KcPi0_MC:
        if (temp_data.Decay[2] > 0) return true;
        return false;
        break;
    case Loader::Xsu2Kcstar2K0Pic_MC:
        if (temp_data.Decay[1] > 0) return true;
        return false;
        break;
    case Loader::Xsu2KcPi0_MC:
        if (temp_data.Decay[5] > 0) return true;
        return false;
        break;
    case Loader::Xsu2K0Pic_MC:
        if (temp_data.Decay[6] > 0) return true;
        return false;
        break;
    case Loader::Xsu2KcPicPic_MC:
        if (temp_data.Decay[7] > 0) return true;
        return false;
        break;
    case Loader::Xsu2K0PicPi0_MC:
        if (temp_data.Decay[8] > 0) return true;
        return false;
        break;
    case Loader::Xsu2KcPicPicPi0_MC:
        if (temp_data.Decay[9] > 0) return true;
        return false;
        break;
    case Loader::Xsu2K0PicPicPic_MC:
        if (temp_data.Decay[10] > 0) return true;
        return false;
        break;
    case Loader::Xsu2KcPicPicPicPic_MC:
        if (temp_data.Decay[11] > 0) return true;
        return false;
        break;
    case Loader::Xsu2K0PicPicPicPi0_MC:
        if (temp_data.Decay[12] > 0) return true;
        return false;
        break;
    case Loader::Xsu2KcPi0Pi0_MC:
        if (temp_data.Decay[13] > 0) return true;
        return false;
        break;
    case Loader::Xsu2K0PicPi0Pi0_MC:
        if (temp_data.Decay[14] > 0) return true;
        return false;
        break;
    case Loader::Xsu2KcPicPicPi0Pi0_MC:
        if (temp_data.Decay[15] > 0) return true;
        return false;
        break;
    case Loader::Xsu2KcKcKc_MC:
        if (temp_data.Decay[16] > 0) return true;
        return false;
        break;
    case Loader::Xsu2KcKcK0Pic_MC:
        if (temp_data.Decay[17] > 0) return true;
        return false;
        break;
    case Loader::Xsu2KcKcKcPi0_MC:
        if (temp_data.Decay[18] > 0) return true;
        return false;
        break;
    case Loader::Xsd2K0_MC:
        if (temp_data.Decay[19] > 0) return true;
        return false;
        break;
    case Loader::Xsd2K0star2KcPic_MC:
        if (temp_data.Decay[20] > 0) return true;
        return false;
        break;
    case Loader::Xsd2K0star2K0Pi0_MC:
        if (temp_data.Decay[21] > 0) return true;
        return false;
        break;
    case Loader::Xsd2KcPic_MC:
        if (temp_data.Decay[24] > 0) return true;
        return false;
        break;
    case Loader::Xsd2K0Pi0_MC:
        if (temp_data.Decay[25] > 0) return true;
        return false;
        break;
    case Loader::Xsd2KcPicPi0_MC:
        if (temp_data.Decay[26] > 0) return true;
        return false;
        break;
    case Loader::Xsd2K0PicPic_MC:
        if (temp_data.Decay[27] > 0) return true;
        return false;
        break;
    case Loader::Xsd2KcPicPicPic_MC:
        if (temp_data.Decay[28] > 0) return true;
        return false;
        break;
    case Loader::Xsd2K0PicPicPi0_MC:
        if (temp_data.Decay[29] > 0) return true;
        return false;
        break;
    case Loader::Xsd2KcPicPicPicPi0_MC:
        if (temp_data.Decay[30] > 0) return true;
        return false;
        break;
    case Loader::Xsd2K0PicPicPicPic_MC:
        if (temp_data.Decay[31] > 0) return true;
        return false;
        break;
    case Loader::Xsd2K0Pi0Pi0_MC:
        if (temp_data.Decay[32] > 0) return true;
        return false;
        break;
    case Loader::Xsd2KcPicPi0Pi0_MC:
        if (temp_data.Decay[33] > 0) return true;
        return false;
        break;
    case Loader::Xsd2K0PicPicPi0Pi0_MC:
        if (temp_data.Decay[34] > 0) return true;
        return false;
        break;
    case Loader::Xsd2KcKcK0_MC:
        if (temp_data.Decay[35] > 0) return true;
        return false;
        break;
    case Loader::Xsd2KcKcKcPic_MC:
        if (temp_data.Decay[36] > 0) return true;
        return false;
        break;
    case Loader::Xsd2KcKcK0Pi0_MC:
        if (temp_data.Decay[37] > 0) return true;
        return false;
        break;
    case Loader::other:
        return true;
        break;
    default:
        printf("ERROR! Input value of TrueIfDecayModeMatch_MC is not appropriate\n");
        exit(1);
        break;
    }

    printf("ERROR! Input value of TrueIfDecayModeMatch_MC is not appropriate\n");
    exit(1);
    return false;
}


void Loader::PrintConfusionMatrix(std::string filename, const char* type, const char* MC_version, const char* category, bool smartmode) {
    if (current_Confusion_matrix > 0) { // allocate new int
        printf("The number of PrintConfusionMatrix should not be larger than 1\n");
        printf("Only first PrintConfusionMatrix is accepted\n");
        return;
    }
    if (DoesItHaveXsBranch == false) {
        printf("ERROR! PrintConfusionMatrix is called when the data does not have Xs branch\n");
        exit(1);
    }

    std::queue<Data> temp_queue;
    temp_queue.swap(TotalData);
    while (!temp_queue.empty()) {
        Data temp = temp_queue.front();
        temp_queue.pop();

        int decaymodeid = -1;
        int decaymodeid_MC = -1;

        for (int i = 0; i < Loader::MAX_NUM_DECAYMODE; i++) { // find reco decay mode
            if (TrueIfDecayModeMatch(temp, static_cast<Loader::DecayMode>(i))) {
                decaymodeid = i;
                break;
            }
        }
        if (decaymodeid == Loader::MAX_NUM_DECAYMODE) {
            printf("ERROR! Reco decay id cannot be found\n");
            exit(1);
        }

        for (int i = 0; i < Loader::MAX_NUM_DECAYMODE_MC; i++) { // find MC decay mode
            if (TrueIfDecayModeMatch_MC(temp, static_cast<Loader::DecayModeMC>(i))) {
                decaymodeid_MC = i;
                break;
            }
        }
        if (decaymodeid_MC == Loader::MAX_NUM_DECAYMODE_MC) {
            printf("ERROR! MC decay id cannot be found\n");
            exit(1);
        }

        if (smartmode == false) Confusion[decaymodeid][decaymodeid_MC]++;
        else {
            if (filename.find("B2Knunu") != std::string::npos) {
                double correction_weight = corrector.GetCorrectionFactor(temp.Decay_syst_ff[index_q2] * temp.Decay_syst_ff[index_q2], "Bplus");
                Confusion[decaymodeid][decaymodeid_MC] = Confusion[decaymodeid][decaymodeid_MC] + ObtainWeight(type, MC_version, category, filename) * correction_weight;
            }
            else if (filename.find("B2Kstarnunu") != std::string::npos) Confusion[decaymodeid][decaymodeid_MC] = Confusion[decaymodeid][decaymodeid_MC] + ObtainWeight(type, MC_version, category, filename);
            else if (filename.find("B2Xsnunu") != std::string::npos) {
                double correction_fragmentation = corrector_Fragmentation.GetCorrectionFactor(temp.Decay, temp.Decay_syst_ff[index_MXs_Bc], Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MC_version);
                Confusion[decaymodeid][decaymodeid_MC] = Confusion[decaymodeid][decaymodeid_MC] + ObtainWeight(type, MC_version, category, filename) * correction_fragmentation;
            }
            else if (filename.find("B02K0nunu") != std::string::npos) {
                double correction_weight = corrector.GetCorrectionFactor(temp.Decay_syst_ff[index_q2] * temp.Decay_syst_ff[index_q2], "Bzero");
                Confusion[decaymodeid][decaymodeid_MC] = Confusion[decaymodeid][decaymodeid_MC] + ObtainWeight(type, MC_version, category, filename) * correction_weight;
            }
            else if (filename.find("B02Kstar0nunu") != std::string::npos) Confusion[decaymodeid][decaymodeid_MC] = Confusion[decaymodeid][decaymodeid_MC] + ObtainWeight(type, MC_version, category, filename);
            else if (filename.find("B02Xsnunu") != std::string::npos) {
                double correction_fragmentation = corrector_Fragmentation.GetCorrectionFactor(temp.Decay, temp.Decay_syst_ff[index_MXs_B0], Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MC_version);
                Confusion[decaymodeid][decaymodeid_MC] = Confusion[decaymodeid][decaymodeid_MC] + ObtainWeight(type, MC_version, category, filename) * correction_fragmentation;
            }
            else { printf("ERROR 142\n"); exit(1); }
        }

        int decaymodeid_MC_for_square = -1;
        if (decaymodeid_MC == 0) decaymodeid_MC_for_square = 0;
        else if (decaymodeid_MC == 1)decaymodeid_MC_for_square = 1;
        else if (decaymodeid_MC == 2)decaymodeid_MC_for_square = 2;
        else if (decaymodeid_MC == 3)decaymodeid_MC_for_square = 1;
        else if (decaymodeid_MC == 4)decaymodeid_MC_for_square = 2;
        else if (decaymodeid_MC == 5)decaymodeid_MC_for_square = 3;
        else if (decaymodeid_MC == 6)decaymodeid_MC_for_square = 4;
        else if (decaymodeid_MC == 7)decaymodeid_MC_for_square = 5;
        else if (decaymodeid_MC == 8)decaymodeid_MC_for_square = 6;
        else if (decaymodeid_MC == 9)decaymodeid_MC_for_square = 7;
        else if (decaymodeid_MC == 10)decaymodeid_MC_for_square = 8;
        else if (decaymodeid_MC == 11)decaymodeid_MC_for_square = 9;
        else if (decaymodeid_MC == 12)decaymodeid_MC_for_square = 10;
        else if (decaymodeid_MC == 13)decaymodeid_MC_for_square = 11;
        else if (decaymodeid_MC == 14)decaymodeid_MC_for_square = 12;
        else if (decaymodeid_MC == 15)decaymodeid_MC_for_square = 13;
        else if (decaymodeid_MC == 16)decaymodeid_MC_for_square = 14;
        else if (decaymodeid_MC == 17)decaymodeid_MC_for_square = 15;
        else if (decaymodeid_MC == 18)decaymodeid_MC_for_square = 16;
        else if (decaymodeid_MC == 19)decaymodeid_MC_for_square = 17;
        else if (decaymodeid_MC == 20)decaymodeid_MC_for_square = 16;
        else if (decaymodeid_MC == 21)decaymodeid_MC_for_square = 17;
        else if (decaymodeid_MC == 22)decaymodeid_MC_for_square = 18;
        else if (decaymodeid_MC == 23)decaymodeid_MC_for_square = 19;
        else if (decaymodeid_MC == 24)decaymodeid_MC_for_square = 20;
        else if (decaymodeid_MC == 25)decaymodeid_MC_for_square = 21;
        else if (decaymodeid_MC == 26)decaymodeid_MC_for_square = 22;
        else if (decaymodeid_MC == 27)decaymodeid_MC_for_square = 23;
        else if (decaymodeid_MC == 28)decaymodeid_MC_for_square = 24;
        else if (decaymodeid_MC == 29)decaymodeid_MC_for_square = 25;
        else if (decaymodeid_MC == 30)decaymodeid_MC_for_square = 26;
        else if (decaymodeid_MC == 31)decaymodeid_MC_for_square = 27;
        else if (decaymodeid_MC == 32)decaymodeid_MC_for_square = 28;
        else if (decaymodeid_MC == 33)decaymodeid_MC_for_square = 29;
        else if (decaymodeid_MC == 34)decaymodeid_MC_for_square = 30;
        if(smartmode == false) Confusion_square[decaymodeid][decaymodeid_MC_for_square]++;
        else {
            if (filename.find("B2Knunu") != std::string::npos) {
                double correction_weight = corrector.GetCorrectionFactor(temp.Decay_syst_ff[index_q2] * temp.Decay_syst_ff[index_q2], "Bplus");
                Confusion_square[decaymodeid][decaymodeid_MC_for_square] = Confusion_square[decaymodeid][decaymodeid_MC_for_square] + ObtainWeight(type, MC_version, category, filename) * correction_weight;
            }
            else if (filename.find("B2Kstarnunu") != std::string::npos) Confusion_square[decaymodeid][decaymodeid_MC_for_square] = Confusion_square[decaymodeid][decaymodeid_MC_for_square] + ObtainWeight(type, MC_version, category, filename);
            else if (filename.find("B2Xsnunu") != std::string::npos) {
                double correction_fragmentation = corrector_Fragmentation.GetCorrectionFactor(temp.Decay, temp.Decay_syst_ff[index_MXs_Bc], Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MC_version);
                Confusion_square[decaymodeid][decaymodeid_MC_for_square] = Confusion_square[decaymodeid][decaymodeid_MC_for_square] + ObtainWeight(type, MC_version, category, filename) * correction_fragmentation;
            }
            else if (filename.find("B02K0nunu") != std::string::npos) {
                double correction_weight = corrector.GetCorrectionFactor(temp.Decay_syst_ff[index_q2] * temp.Decay_syst_ff[index_q2], "Bzero");
                Confusion_square[decaymodeid][decaymodeid_MC_for_square] = Confusion_square[decaymodeid][decaymodeid_MC_for_square] + ObtainWeight(type, MC_version, category, filename) * correction_weight;
            }
            else if (filename.find("B02Kstar0nunu") != std::string::npos) Confusion_square[decaymodeid][decaymodeid_MC_for_square] = Confusion_square[decaymodeid][decaymodeid_MC_for_square] + ObtainWeight(type, MC_version, category, filename);
            else if (filename.find("B02Xsnunu") != std::string::npos) {
                double correction_fragmentation = corrector_Fragmentation.GetCorrectionFactor(temp.Decay, temp.Decay_syst_ff[index_MXs_B0], Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MC_version);
                Confusion_square[decaymodeid][decaymodeid_MC_for_square] = Confusion_square[decaymodeid][decaymodeid_MC_for_square] + ObtainWeight(type, MC_version, category, filename) * correction_fragmentation;
            }
            else { printf("ERROR 142\n"); exit(1); }
        }

        TotalData.push(temp);
    }

    Confusion_matrixIsOn = true;
    current_Confusion_matrix++;
}

void Loader::DvetoFor(Loader::Variable variable, int i, double min, double max) { // remove min < X < max
    std::queue<Data> temp_queue;
    while (!TotalData.empty()) {
        Data temp_data = TotalData.front();
        TotalData.pop();
        if (variable == Loader::Upsilon) {
            if (min < temp_data.Upsilon_info[i] && temp_data.Upsilon_info[i] < max) {}
            else temp_queue.push(temp_data);
        }
        else if (variable == Loader::Bsig) {
            if (min < temp_data.Bsig_info[i] && temp_data.Bsig_info[i] < max) {}
            else temp_queue.push(temp_data);
        }
        else if (variable == Loader::Btag) {
            if (min < temp_data.Btag_info[i] && temp_data.Btag_info[i] < max) {}
            else temp_queue.push(temp_data);
        }
        else {
            printf("ERROR! 499\n");
            exit(1);
        }
    }
    TotalData.swap(temp_queue);
}

void Loader::BsigFitConvergeFor(Loader::Variable variable, int i) { // remove X < -0.5 depending on decay mode
    std::queue<Data> temp_queue;
    while (!TotalData.empty()) {
        Data temp_data = TotalData.front();
        TotalData.pop();

        if (TrueIfDecayModeMatch(temp_data, Loader::B2Kc) || TrueIfDecayModeMatch(temp_data, Loader::B2KcPi0) || TrueIfDecayModeMatch(temp_data, Loader::B02Ks0) || TrueIfDecayModeMatch(temp_data, Loader::B02Ks0Pi0)) {
            temp_queue.push(temp_data);
            continue;
        }

        if (variable == Loader::Upsilon) {
            if (-0.5 < temp_data.Upsilon_info[i]) temp_queue.push(temp_data);
        }
        else if (variable == Loader::Bsig) {
            if (-0.5 < temp_data.Bsig_info[i]) temp_queue.push(temp_data);
        }
        else if (variable == Loader::Btag) {
            if (-0.5 < temp_data.Btag_info[i]) temp_queue.push(temp_data);
        }
        else {
            printf("ERROR! 500\n");
            exit(1);
        }
    }
    TotalData.swap(temp_queue);
}

void Loader::OnlySelectDvetoTypeFor(Loader::Variable variable, int Dchargedvetomassindex, int DchargedvetodmIDindex, int Dneutralvetomassindex, int DneutralvetodmIDindex, Loader::Dvetotype type) { // select specific D veto type
    std::queue<Data> temp_queue;
    while (!TotalData.empty()) {
        Data temp_data = TotalData.front();
        TotalData.pop();

        double temp_dmID_Dcharged;
        double temp_dmID_Dneutral;
        if (variable == Loader::Upsilon) {
            temp_dmID_Dcharged = temp_data.Upsilon_info[DchargedvetodmIDindex];
            temp_dmID_Dneutral = temp_data.Upsilon_info[DneutralvetodmIDindex];
        }
        else if (variable == Loader::Bsig) {
            temp_dmID_Dcharged = temp_data.Bsig_info[DchargedvetodmIDindex];
            temp_dmID_Dneutral = temp_data.Bsig_info[DneutralvetodmIDindex];
        }
        else if (variable == Loader::Btag) {
            temp_dmID_Dcharged = temp_data.Btag_info[DchargedvetodmIDindex];
            temp_dmID_Dneutral = temp_data.Btag_info[DneutralvetodmIDindex];
        }

        if (type == Loader::Dchargedwithoutpizero) {
            if ((0.5 < temp_dmID_Dcharged && temp_dmID_Dcharged < 3.5) || (5.5 < temp_dmID_Dcharged && temp_dmID_Dcharged < 8.5)) temp_queue.push(temp_data);
        }
        else if (type == Loader::Dchargedwithpizero) {
            if ((3.5 < temp_dmID_Dcharged && temp_dmID_Dcharged  < 5.5) || (8.5 < temp_dmID_Dcharged && temp_dmID_Dcharged < 11.5)) temp_queue.push(temp_data);
        }
        else if (type == Loader::Dneutralwithoutpizero) {
            if ((0.5 < temp_dmID_Dneutral && temp_dmID_Dneutral < 3.5) || (6.5 < temp_dmID_Dneutral && temp_dmID_Dneutral < 7.5)) temp_queue.push(temp_data);
        }
        else if (type == Loader::Dneutralwithpizero) {
            if ((3.5 < temp_dmID_Dneutral && temp_dmID_Dneutral < 6.5) || (7.5 < temp_dmID_Dneutral && temp_dmID_Dneutral < 8.5)) temp_queue.push(temp_data);
        }
        else {
            printf("ERROR! 963\n");
            exit(1);
        }
    }
    TotalData.swap(temp_queue);
}

void Loader::DvetoAboutSpecificTypeFor(Loader::Variable variable, int Dchargedvetomassindex, int DchargedvetodmIDindex, int Dneutralvetomassindex, int DneutralvetodmIDindex, Loader::Dvetotype type, double minM, double maxM) { // reject minM < D < max for specific type
    std::queue<Data> temp_queue;
    while (!TotalData.empty()) {
        Data temp_data = TotalData.front();
        TotalData.pop();

        double temp_dmID_Dcharged;
        double temp_M_Dcharged;
        double temp_dmID_Dneutral;
        double temp_M_Dneutral;
        if (variable == Loader::Upsilon) {
            temp_dmID_Dcharged = temp_data.Upsilon_info[DchargedvetodmIDindex];
            temp_M_Dcharged = temp_data.Upsilon_info[Dchargedvetomassindex];
            temp_dmID_Dneutral = temp_data.Upsilon_info[DneutralvetodmIDindex];
            temp_M_Dneutral = temp_data.Upsilon_info[Dneutralvetomassindex];
        }
        else if (variable == Loader::Bsig) {
            temp_dmID_Dcharged = temp_data.Bsig_info[DchargedvetodmIDindex];
            temp_M_Dcharged = temp_data.Bsig_info[Dchargedvetomassindex];
            temp_dmID_Dneutral = temp_data.Bsig_info[DneutralvetodmIDindex];
            temp_M_Dneutral = temp_data.Bsig_info[Dneutralvetomassindex];
        }
        else if (variable == Loader::Btag) {
            temp_dmID_Dcharged = temp_data.Btag_info[DchargedvetodmIDindex];
            temp_M_Dcharged = temp_data.Btag_info[Dchargedvetomassindex];
            temp_dmID_Dneutral = temp_data.Btag_info[DneutralvetodmIDindex];
            temp_M_Dneutral = temp_data.Btag_info[Dneutralvetomassindex];
        }

        if (type == Loader::Dchargedwithoutpizero) {
            if ((0.5 < temp_dmID_Dcharged && temp_dmID_Dcharged < 3.5) || (5.5 < temp_dmID_Dcharged && temp_dmID_Dcharged < 8.5)) {
                if(!(minM < temp_M_Dcharged && temp_M_Dcharged < maxM)) temp_queue.push(temp_data);
            }
            else temp_queue.push(temp_data);
        }
        else if (type == Loader::Dchargedwithpizero) {
            if ((3.5 < temp_dmID_Dcharged && temp_dmID_Dcharged < 5.5) || (8.5 < temp_dmID_Dcharged && temp_dmID_Dcharged < 11.5)) {
                if (!(minM < temp_M_Dcharged && temp_M_Dcharged < maxM)) temp_queue.push(temp_data);
            }
            else temp_queue.push(temp_data);
        }
        else if (type == Loader::Dneutralwithoutpizero) {
            if ((0.5 < temp_dmID_Dneutral && temp_dmID_Dneutral < 3.5) || (6.5 < temp_dmID_Dneutral && temp_dmID_Dneutral < 7.5)) {
                if (!(minM < temp_M_Dneutral && temp_M_Dneutral < maxM)) temp_queue.push(temp_data);
            }
            else temp_queue.push(temp_data);
        }
        else if (type == Loader::Dneutralwithpizero) {
            if ((3.5 < temp_dmID_Dneutral && temp_dmID_Dneutral < 6.5) || (7.5 < temp_dmID_Dneutral && temp_dmID_Dneutral < 8.5)) {
                if (!(minM < temp_M_Dneutral && temp_M_Dneutral < maxM)) temp_queue.push(temp_data);
            }
            else temp_queue.push(temp_data);
        }
        else {
            printf("ERROR! 964\n");
            exit(1);
        }
    }
    TotalData.swap(temp_queue);
}

void Loader::PrintFOM(std::string filename, const char* type, const char* MC_version, const char* category, bool smartmode) {
    if (current_FOM > 0) { // allocate new int
        printf("The number of PrintFOM should not be larger than 1\n");
        printf("Only first PrintFOM is accepted\n");
        return;
    }
    if (DoesItHaveMVAOutput == false) {
        printf("ERROR! PrintFOM is called when the data does not have MVA output\n");
        exit(1);
    }
    typedef struct labels {
        int __experiment__;
        int __run__;
        unsigned int __event__;
        int __ncandidates__;
    } Labels;

    for (int i = 0; i < Nstep; i++) {
        for (int j = 0; j < Nstep; j++) {
            std::queue<Data> temp_queue;
            temp_queue.swap(TotalData);

            double BB_output = start + (end - start) * i / Nstep;
            double Continuum_output = start + (end - start) * j / Nstep;

            std::vector<Labels> label_list;
            double EVT_num = 0.0;

            while (!temp_queue.empty()) {
                Data temp = temp_queue.front();
                temp_queue.pop();

                if (temp.MVA_BB > BB_output && temp.MVA_Continuum > Continuum_output) {
                    bool overlap = false;
                    for (unsigned int k = 0; k < label_list.size(); k++) {
                        if (label_list.at(k).__experiment__ == temp.__experiment__ && label_list.at(k).__run__ == temp.__run__ && label_list.at(k).__event__ == temp.__event__ && label_list.at(k).__ncandidates__ == temp.__ncandidates__) {
                            overlap = true;
                        }
                    }
                    if (overlap == false) {
                        if (smartmode == false) EVT_num = EVT_num + 1.0;
                        else {
                            if (filename.find("B2Knunu") != std::string::npos) {
                                double correction_weight = corrector.GetCorrectionFactor(temp.Decay_syst_ff[index_q2] * temp.Decay_syst_ff[index_q2], "Bplus");
                                EVT_num = EVT_num + ObtainWeight(type, MC_version, category, filename) * correction_weight;
                            }
                            else if (filename.find("B2Kstarnunu") != std::string::npos) EVT_num = EVT_num + ObtainWeight(type, MC_version, category, filename);
                            else if (filename.find("B2Xsnunu") != std::string::npos) {
                                double correction_fragmentation = corrector_Fragmentation.GetCorrectionFactor(temp.Decay, temp.Decay_syst_ff[index_MXs_Bc], Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MC_version);
                                EVT_num = EVT_num + ObtainWeight(type, MC_version, category, filename) * correction_fragmentation;
                            }
                            else if (filename.find("B02K0nunu") != std::string::npos) {
                                double correction_weight = corrector.GetCorrectionFactor(temp.Decay_syst_ff[index_q2] * temp.Decay_syst_ff[index_q2], "Bzero");
                                EVT_num = EVT_num + ObtainWeight(type, MC_version, category, filename) * correction_weight;
                            }
                            else if (filename.find("B02Kstar0nunu") != std::string::npos) EVT_num = EVT_num + ObtainWeight(type, MC_version, category, filename);
                            else if (filename.find("B02Xsnunu") != std::string::npos) {
                                double correction_fragmentation = corrector_Fragmentation.GetCorrectionFactor(temp.Decay, temp.Decay_syst_ff[index_MXs_B0], Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MC_version);
                                EVT_num = EVT_num + ObtainWeight(type, MC_version, category, filename) * correction_fragmentation;
                            }
                            else { EVT_num = EVT_num + ObtainWeight(type, MC_version, category, filename); }
                        }
                        Labels temp_Labels;
                        temp_Labels.__experiment__ = temp.__experiment__;
                        temp_Labels.__run__ = temp.__run__;
                        temp_Labels.__event__ = temp.__event__;
                        temp_Labels.__ncandidates__ = temp.__ncandidates__;
                        label_list.push_back(temp_Labels);
                    }
                }

                TotalData.push(temp);
            }
            FOM_Matrix[i][j] = FOM_Matrix[i][j] + EVT_num;
        }
    }

    FOMIsOn = true;
    current_FOM++;
}

void Loader::PrintFOM1D(std::string filename, const char* type, const char* MC_version, const char* category, bool smartmode) {
    if (current_FOM > 0) { // allocate new int
        printf("The number of PrintFOM should not be larger than 1\n");
        printf("Only first PrintFOM is accepted\n");
        return;
    }
    if (DoesItHaveMVAOutput == false) {
        printf("ERROR! PrintFOM is called when the data does not have MVA output\n");
        exit(1);
    }
    typedef struct labels {
        int __experiment__;
        int __run__;
        unsigned int __event__;
        int __ncandidates__;
    } Labels;

    for (int j = 0; j < Nstep; j++) {
        std::queue<Data> temp_queue;
        temp_queue.swap(TotalData);

        double BB_output = start + (end - start) * j / Nstep;

        std::vector<Labels> label_list;
        double EVT_num = 0.0;

        while (!temp_queue.empty()) {
            Data temp = temp_queue.front();
            temp_queue.pop();

            if (temp.MVA_BB > BB_output) {
                bool overlap = false;
                for (unsigned int k = 0; k < label_list.size(); k++) {
                    if (label_list.at(k).__experiment__ == temp.__experiment__ && label_list.at(k).__run__ == temp.__run__ && label_list.at(k).__event__ == temp.__event__ && label_list.at(k).__ncandidates__ == temp.__ncandidates__) {
                        overlap = true;
                    }
                }
                if (overlap == false) {
                    if (smartmode == false) EVT_num = EVT_num + 1.0;
                    else {
                        if (filename.find("B2Knunu") != std::string::npos) {
                            double correction_weight = corrector.GetCorrectionFactor(temp.Decay_syst_ff[index_q2] * temp.Decay_syst_ff[index_q2], "Bplus");
                            EVT_num = EVT_num + ObtainWeight(type, MC_version, category, filename) * correction_weight;
                        }
                        else if (filename.find("B2Kstarnunu") != std::string::npos) EVT_num = EVT_num + ObtainWeight(type, MC_version, category, filename);
                        else if (filename.find("B2Xsnunu") != std::string::npos) {
                            double correction_fragmentation = corrector_Fragmentation.GetCorrectionFactor(temp.Decay, temp.Decay_syst_ff[index_MXs_Bc], Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MC_version);
                            EVT_num = EVT_num + ObtainWeight(type, MC_version, category, filename) * correction_fragmentation;
                        }
                        else if (filename.find("B02K0nunu") != std::string::npos) {
                            double correction_weight = corrector.GetCorrectionFactor(temp.Decay_syst_ff[index_q2] * temp.Decay_syst_ff[index_q2], "Bzero");
                            EVT_num = EVT_num + ObtainWeight(type, MC_version, category, filename) * correction_weight;
                        }
                        else if (filename.find("B02Kstar0nunu") != std::string::npos) EVT_num = EVT_num + ObtainWeight(type, MC_version, category, filename);
                        else if (filename.find("B02Xsnunu") != std::string::npos) {
                            double correction_fragmentation = corrector_Fragmentation.GetCorrectionFactor(temp.Decay, temp.Decay_syst_ff[index_MXs_B0], Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MC_version);
                            EVT_num = EVT_num + ObtainWeight(type, MC_version, category, filename) * correction_fragmentation;
                        }
                        else EVT_num = EVT_num + ObtainWeight(type, MC_version, category, filename);
                    }
                    Labels temp_Labels;
                    temp_Labels.__experiment__ = temp.__experiment__;
                    temp_Labels.__run__ = temp.__run__;
                    temp_Labels.__event__ = temp.__event__;
                    temp_Labels.__ncandidates__ = temp.__ncandidates__;
                    label_list.push_back(temp_Labels);
                }
            }

            TotalData.push(temp);
        }
        FOM_Matrix[0][j] = FOM_Matrix[0][j] + EVT_num;
    }

    FOMIsOn = true;
    current_FOM++;
}

void Loader::MVACut(double OBB, double Oqq, Loader::MassRegion massRegion) {
    if (DoesItHaveMVAOutput == false) {
        printf("ERROR! MVACut is called when the data does not have MVA output\n");
        exit(1);
    }

    std::queue<Data> temp_queue;
    while (!TotalData.empty()) {
        Data temp_data = TotalData.front();
        TotalData.pop();

        if (massRegion == Loader::SmallMass) {
            if (temp_data.Bsig_info[6] > 1.1) temp_queue.push(temp_data);
            else {
                if(temp_data.MVA_BB > OBB && temp_data.MVA_Continuum > Oqq) temp_queue.push(temp_data);
            }
        }
        else if (massRegion == Loader::LargeMass) {
            if (temp_data.Bsig_info[6] < 1.1) temp_queue.push(temp_data);
            else {
                if (temp_data.MVA_BB > OBB && temp_data.MVA_Continuum > Oqq) temp_queue.push(temp_data);
            }
        }

    }
    TotalData.swap(temp_queue);
}

void Loader::CountMCEvent(std::string filename, const char* type, const char* MC_version, const char* category, bool smartmode) {
    if (current_MCcount > 0) { // allocate new int
        printf("The number of CountMCEvent should not be larger than 1\n");
        printf("Only first CountMCEvent is accepted\n");
        return;
    }
    if (DoesItHaveXsBranch == false) {
        printf("ERROR! CountMCEvent is called when the data does not have Xs branch\n");
        exit(1);
    }
    typedef struct labels {
        int __experiment__;
        int __run__;
        unsigned int __event__;
        int __ncandidates__;
    } Labels;

    std::queue<Data> temp_queue;
    temp_queue.swap(TotalData);

    std::vector<Labels> label_list;

    while (!temp_queue.empty()) {
        Data temp = temp_queue.front();
        temp_queue.pop();

        bool overlap = false;
        for (unsigned int k = 0; k < label_list.size(); k++) {
            if (label_list.at(k).__experiment__ == temp.__experiment__ && label_list.at(k).__run__ == temp.__run__ && label_list.at(k).__event__ == temp.__event__ && label_list.at(k).__ncandidates__ == temp.__ncandidates__) {
                overlap = true;
            }
        }
        if (overlap == false) {
            Labels temp_Labels;
            temp_Labels.__experiment__ = temp.__experiment__;
            temp_Labels.__run__ = temp.__run__;
            temp_Labels.__event__ = temp.__event__;
            temp_Labels.__ncandidates__ = temp.__ncandidates__;
            label_list.push_back(temp_Labels);

            int decaymodeid_MC = -1;
            for (int i = 0; i < Loader::MAX_NUM_DECAYMODE_MC; i++) { // find MC decay mode
                if (TrueIfDecayModeMatch_MC(temp, static_cast<Loader::DecayModeMC>(i))) {
                    decaymodeid_MC = i;
                    break;
                }
            }
            if (decaymodeid_MC == Loader::MAX_NUM_DECAYMODE_MC) {
                printf("ERROR! MC decay id cannot be found\n");
                exit(1);
            }

            if (smartmode == false) MCcount[decaymodeid_MC] = MCcount[decaymodeid_MC] + 1;
            else {
                if (filename.find("B2Knunu") != std::string::npos) {
                    double correction_weight = corrector.GetCorrectionFactor(temp.Decay_syst_ff[index_q2] * temp.Decay_syst_ff[index_q2], "Bplus");
                    MCcount[decaymodeid_MC] = MCcount[decaymodeid_MC] + ObtainWeight(type, MC_version, category, filename) * correction_weight;
                }
                else if (filename.find("B2Kstarnunu") != std::string::npos) MCcount[decaymodeid_MC] = MCcount[decaymodeid_MC] + ObtainWeight(type, MC_version, category, filename);
                else if (filename.find("B2Xsnunu") != std::string::npos) {
                    double correction_fragmentation = corrector_Fragmentation.GetCorrectionFactor(temp.Decay, temp.Decay_syst_ff[index_MXs_Bc], Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MC_version);
                    MCcount[decaymodeid_MC] = MCcount[decaymodeid_MC] + ObtainWeight(type, MC_version, category, filename) * correction_fragmentation;
                }
                else if (filename.find("B02K0nunu") != std::string::npos) {
                    double correction_weight = corrector.GetCorrectionFactor(temp.Decay_syst_ff[index_q2] * temp.Decay_syst_ff[index_q2], "Bzero");
                    MCcount[decaymodeid_MC] = MCcount[decaymodeid_MC] + ObtainWeight(type, MC_version, category, filename) * correction_weight;
                }
                else if (filename.find("B02Kstar0nunu") != std::string::npos) MCcount[decaymodeid_MC] = MCcount[decaymodeid_MC] + ObtainWeight(type, MC_version, category, filename);
                else if (filename.find("B02Xsnunu") != std::string::npos) {
                    double correction_fragmentation = corrector_Fragmentation.GetCorrectionFactor(temp.Decay, temp.Decay_syst_ff[index_MXs_B0], Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MC_version);
                    MCcount[decaymodeid_MC] = MCcount[decaymodeid_MC] + ObtainWeight(type, MC_version, category, filename) * correction_fragmentation;
                }
                else { printf("ERROR 142\n"); exit(1); }
            }
        }

        TotalData.push(temp);
    }

    MCcountOn = true;
    current_MCcount++;
}

void Loader::SelectDecayModeOf(Loader::DecayMode decaymode) {
    std::queue<Data> temp_queue;
    while (!TotalData.empty()) {
        Data temp_data = TotalData.front();
        TotalData.pop();

        if (TrueIfDecayModeMatch(temp_data, decaymode)) temp_queue.push(temp_data);

    }
    TotalData.swap(temp_queue);
}

void Loader::RejectDecayModeOf(Loader::DecayMode decaymode) {
    std::queue<Data> temp_queue;
    while (!TotalData.empty()) {
        Data temp_data = TotalData.front();
        TotalData.pop();

        if (TrueIfDecayModeMatch(temp_data, decaymode)) {}
        else temp_queue.push(temp_data);

    }
    TotalData.swap(temp_queue);
}

void Loader::BeamEnergyCorrectionFromDeltaE(int index_pBcms, int index_EBcms, int index_Mbc, int index_deltaE, double targetEbeamstar, bool IsItBtag) {
    std::queue<Data> temp_queue;
    while (!TotalData.empty()) {
        Data temp_data = TotalData.front();
        TotalData.pop();

        if (IsItBtag) {
            double Ebeamstar = temp_data.Btag_info[index_EBcms] - temp_data.Btag_info[index_deltaE]; // EBstar - deltaE
            temp_data.Btag_info[index_Mbc] = std::sqrt(targetEbeamstar * targetEbeamstar - targetEbeamstar * targetEbeamstar * temp_data.Btag_info[index_pBcms] * temp_data.Btag_info[index_pBcms] / (Ebeamstar * Ebeamstar)); // correct Mbc
            temp_data.Btag_info[index_deltaE] = (targetEbeamstar / Ebeamstar) * temp_data.Btag_info[index_EBcms] - targetEbeamstar; // correct deltaE
        }
        else {
            double Ebeamstar = temp_data.Bsig_info[index_EBcms] - temp_data.Bsig_info[index_deltaE]; // EBstar - deltaE
            temp_data.Bsig_info[index_Mbc] = std::sqrt(targetEbeamstar * targetEbeamstar - targetEbeamstar * targetEbeamstar * temp_data.Bsig_info[index_pBcms] * temp_data.Bsig_info[index_pBcms] / (Ebeamstar * Ebeamstar)); // correct Mbc
            temp_data.Bsig_info[index_deltaE] = (targetEbeamstar / Ebeamstar) * temp_data.Bsig_info[index_EBcms] - targetEbeamstar; // correct deltaE
        }

        temp_queue.push(temp_data);
    }
    TotalData.swap(temp_queue);
}

double ObtainWeight(const char* type, const char* MC_version, const char* category, std::string filename) {
    if (strcmp(MC_version, "data") == 0) return 1.0; // no weight if it is data no matter what other values are
    else if ((strcmp(MC_version, "MC15ri") == 0) && (strcmp(category, "test") == 0)) { // MC15ri test
        if ((strcmp(type, "SIGNAL") == 0)) {
            if (filename.find("B2Knunu") != std::string::npos) return Scale_Kplus_test;
            else if (filename.find("B2Kstarnunu") != std::string::npos) return Scale_Kplusstar_test;
            else if (filename.find("B2Xsnunu") != std::string::npos) return Scale_Xsu_nonresonant_test;
            else if (filename.find("B02K0nunu") != std::string::npos) return Scale_K0_test;
            else if (filename.find("B02Kstar0nunu") != std::string::npos) return Scale_K0star_test;
            else if (filename.find("B02Xsnunu") != std::string::npos) return Scale_Xsd_nonresonant_test;
            else { printf("[ObtainWeight] undefined type for SIGNAL\n"); exit(1); }
        }
        else if ((strcmp(type, "CHG") == 0)) return Scale_CHG_test;
        else if ((strcmp(type, "MIX") == 0)) return Scale_MIX_test;
        else if ((strcmp(type, "UUBAR") == 0)) return Scale_UUBAR_test;
        else if ((strcmp(type, "DDBAR") == 0)) return Scale_DDBAR_test;
        else if ((strcmp(type, "SSBAR") == 0)) return Scale_SSBAR_test;
        else if ((strcmp(type, "CHARM") == 0)) return Scale_CHARM_test;
        else { printf("[ObtainWeight] undefined type for MC15ri test\n"); exit(1); }
    }
    else if ((strcmp(MC_version, "MC15ri") == 0) && (strcmp(category, "train") == 0)) { // MC15ri train
        if ((strcmp(type, "SIGNAL") == 0)) {
            if (filename.find("B2Knunu") != std::string::npos) return Scale_Kplus_train;
            else if (filename.find("B2Kstarnunu") != std::string::npos) return Scale_Kplusstar_train;
            else if (filename.find("B2Xsnunu") != std::string::npos) return Scale_Xsu_nonresonant_train;
            else if (filename.find("B02K0nunu") != std::string::npos) return Scale_K0_train;
            else if (filename.find("B02Kstar0nunu") != std::string::npos) return Scale_K0star_train;
            else if (filename.find("B02Xsnunu") != std::string::npos) return Scale_Xsd_nonresonant_train;
            else { printf("[ObtainWeight] undefined type for SIGNAL\n"); exit(1); }
        }
        else if ((strcmp(type, "CHG") == 0)) return Scale_CHG_train;
        else if ((strcmp(type, "MIX") == 0)) return Scale_MIX_train;
        else if ((strcmp(type, "UUBAR") == 0)) return Scale_UUBAR_train;
        else if ((strcmp(type, "DDBAR") == 0)) return Scale_DDBAR_train;
        else if ((strcmp(type, "SSBAR") == 0)) return Scale_SSBAR_train;
        else if ((strcmp(type, "CHARM") == 0)) return Scale_CHARM_train;
        else { printf("[ObtainWeight] undefined type for MC15ri train\n"); exit(1); }
    }
    else if ((strcmp(MC_version, "MC15ri") == 0) && (strcmp(category, "validation") == 0)) { // MC15ri validation
        if ((strcmp(type, "SIGNAL") == 0)) {
            if (filename.find("CHG") != std::string::npos) return Scale_CHG_validation; // it is Jpsi signal analysis with generic CHG sample
            else if (filename.find("MIX") != std::string::npos) return Scale_MIX_validation; // it is Jpsi signal analysis with generic CHG sample
            return 1.0;
        }
        else if ((strcmp(type, "CHG") == 0)) return Scale_CHG_validation;
        else if ((strcmp(type, "MIX") == 0)) return Scale_MIX_validation;
        else if ((strcmp(type, "UUBAR") == 0)) return Scale_UUBAR_validation;
        else if ((strcmp(type, "DDBAR") == 0)) return Scale_DDBAR_validation;
        else if ((strcmp(type, "SSBAR") == 0)) return Scale_SSBAR_validation;
        else if ((strcmp(type, "CHARM") == 0)) return Scale_CHARM_validation;
        else { printf("[ObtainWeight] undefined type for MC15ri validation\n"); exit(1); }
    }
    else if ((strcmp(MC_version, "MC15rd") == 0) && (strcmp(category, "validation") == 0)) { // MC15rd validation
        if ((strcmp(type, "SIGNAL") == 0)) {
            if (filename.find("B2Knunu") != std::string::npos) return Scale_Kplus_validation_MC15rd;
            else if (filename.find("B2Kstarnunu") != std::string::npos) return Scale_Kplusstar_validation_MC15rd;
            else if (filename.find("B2Xsnunu") != std::string::npos) return Scale_Xsu_nonresonant_validation_MC15rd;
            else if (filename.find("B02K0nunu") != std::string::npos) return Scale_K0_validation_MC15rd;
            else if (filename.find("B02Kstar0nunu") != std::string::npos) return Scale_K0star_validation_MC15rd;
            else if (filename.find("B02Xsnunu") != std::string::npos) return Scale_Xsd_nonresonant_validation_MC15rd;
            else if (filename.find("CHG") != std::string::npos) return Scale_CHG_validation_MC15rd; // it is Jpsi signal analysis with generic CHG sample
            else if (filename.find("MIX") != std::string::npos) return Scale_MIX_validation_MC15rd; // it is Jpsi signal analysis with generic CHG sample
            return 1.0; // just data
        }
        else if ((strcmp(type, "CHG") == 0)) return Scale_CHG_validation_MC15rd;
        else if ((strcmp(type, "MIX") == 0)) return Scale_MIX_validation_MC15rd;
        else if ((strcmp(type, "UUBAR") == 0)) return Scale_UUBAR_validation_MC15rd;
        else if ((strcmp(type, "DDBAR") == 0)) return Scale_DDBAR_validation_MC15rd;
        else if ((strcmp(type, "SSBAR") == 0)) return Scale_SSBAR_validation_MC15rd;
        else if ((strcmp(type, "CHARM") == 0)) return Scale_CHARM_validation_MC15rd;
        else { printf("[ObtainWeight] undefined type for MC15ri validation\n"); exit(1); }
    }

    printf("[ObtainWeight] no matched case!\n");
    exit(1);
    return 1.0;
}

int main(int argc, char* argv[]) {
    /*
    * argv[1]: dirname
    * argv[2]: outputname
    * argv[3]: output path
    * argv[4]: sample type: {SIGNAL|CHG|MIX|UUBAR|DDBAR|SSBAR|CHARM}
    * argv[5]: MC version: {data|MC15ri|MC15rd}
    * argv[6]: sample category: {test|train|validation}
    */

    std::vector<std::string> names;

    load_files(argv[1], &names);

    Loader loader;

    for (unsigned int i = 0; i < names.size(); i++) {
        loader.initialize();

        TFile* input_file = new TFile((argv[1] + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());
        loader.GetData(input_file);
        if (loader.event_info_is_valid() == false) { printf("error!\n"); return 0; }

        std::string::size_type const p(names.at(i).find_last_of('.'));
        std::string file_without_extension = names.at(i).substr(0, p);

        loader.Cut(Loader::Upsilon, 56, Loader::smaller_than, 0.5); // nloose track
        loader.Cut(Loader::Upsilon, 6, Loader::smaller_than, 0.5); // N pi0
        loader.Cut(Loader::Upsilon, 5, Loader::smaller_than, 0.5); // N KS0
        loader.Cut(Loader::Bsig, 64, Loader::larger_than, -0.5); // Bsig vertex fit
        loader.Cut(Loader::Btag, 6, Loader::larger_than, -0.5); // Btag vertex fit
        loader.PrintInformation(std::string("========== inital =========="), names.at(i), argv[4], argv[5], argv[6], true);

        loader.PrintSeparateRootFile(std::string(argv[3]) + "/before_Mbc_cut/" + file_without_extension + std::string("_before_Mbc_cut.root"));
        loader.Cut(Loader::Btag, 1, Loader::larger_than, 5.2);
        loader.PrintInformation(std::string("========== Mbc > 5.2 =========="), names.at(i), argv[4], argv[5], argv[6], true);
        //loader.DrawTH2F("MbcVSdeltaE_after_Mbc_strict_cut", ";Mbc of B_{tag} [GeV];#DeltaE of B_{tag} [GeV]", 100, 5.24, 5.3, 100, -0.2, 0.2, Loader::Btag, 1, Loader::Btag, 2);

        loader.PrintSeparateRootFile(std::string(argv[3]) + "/before_delE_cut/" + file_without_extension + std::string("_before_delE_cut.root"));
        loader.Cut(Loader::Btag, 2, Loader::larger_than, -0.2);
        loader.Cut(Loader::Btag, 2, Loader::smaller_than, 0.2);
        loader.PrintInformation(std::string("========== abs(deltaE) < 0.2 =========="), names.at(i), argv[4], argv[5], argv[6], true);
        //loader.DrawTH2F("MbcVSdeltaE_after_deltaE_strict_cut", ";Mbc of B_{tag} [GeV];#DeltaE of B_{tag} [GeV]", 100, 5.24, 5.3, 100, -0.2, 0.2, Loader::Btag, 1, Loader::Btag, 2);

        loader.PrintSeparateRootFile(std::string(argv[3]) + "/before_Eecl_cut/" + file_without_extension + std::string("_before_Eecl_cut.root"));
        //loader.DrawTH1F("missing_momentum_theta_after_npi0_cut", "#theta_{missing};#theta_{missing} [rad];evt", 100, 0, 3.2, Loader::Upsilon, 7);
        loader.Cut(Loader::Upsilon, 69, Loader::smaller_than, 1.3);
        loader.PrintInformation(std::string("========== Eecl_v200 < 1.3 =========="), names.at(i), argv[4], argv[5], argv[6], true);

        loader.PrintSeparateRootFile(std::string(argv[3]) + "/before_missing_momentum_theta_cut/" + file_without_extension + std::string("_before_missing_momentum_theta_cut.root"));
        //loader.DrawTH1F("missing_momentum_theta_after_npi0_cut", "#theta_{missing};#theta_{missing} [rad];evt", 100, 0, 3.2, Loader::Upsilon, 7);
        loader.Cut(Loader::Upsilon, 7, Loader::smaller_than, 2.618);
        loader.Cut(Loader::Upsilon, 7, Loader::larger_than, 0.297);
        loader.PrintInformation(std::string("========== 0.297 < missing momentum theta < 2.618 =========="), names.at(i), argv[4], argv[5], argv[6], true);

        loader.PrintSeparateRootFile(std::string(argv[3]) + "/before_psig_cut/" + file_without_extension + std::string("_before_psig_cut.root"));
        //loader.DrawTH1F("momentum_Bsig_after_missing_theta_cut", "momentum of B_{sig} at CMS;momentum [GeV];evt", 100, 0, 3.2, Loader::Bsig, 4);
        loader.Cut(Loader::Bsig, 4, Loader::smaller_than, 2.96);
        loader.Cut(Loader::Bsig, 4, Loader::larger_than, 0.5);
        loader.PrintInformation(std::string("========== 0.5 < momentum of signal side < 2.96 =========="), names.at(i), argv[4], argv[5], argv[6], true);

        loader.PrintSeparateRootFile(std::string(argv[3]) + "/before_Dveto_cut/" + file_without_extension + std::string("_before_Dveto_cut.root"));
        //loader.DrawTH1F("Bsig_M_Xs", "mass of X_{s};M_{Xs} [GeV];evt", 100, 0, 3.5, Loader::Bsig, 6);
        loader.DvetoFor(Loader::Bsig, 6, 1.84, 1.89);
        loader.PrintInformation(std::string("========== D veto =========="), names.at(i), argv[4], argv[5], argv[6], true);

        loader.BCS(Loader::Btag, 5, Loader::Highest);
        if (loader.IsBCSValid() == false) {
            printf("ERROR! it is not valid\n");
            exit(1);
        }
        loader.PrintInformation(std::string("========== BCS =========="), names.at(i), argv[4], argv[5], argv[6], true);

        loader.PrintSeparateRootFile(std::string(argv[3]) + "/final_output/" + file_without_extension + std::string("_") + std::string(argv[2]) + ".root");
        loader.ConvertIntoSeparateDataFile(std::string(argv[3]) + "/final_output_data/" + file_without_extension + std::string("_") + std::string(argv[2]) + "_data.root");
    }
    loader.End();

    return 0;
}
