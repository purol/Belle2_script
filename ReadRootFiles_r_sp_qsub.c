// last update: 2021-11-09
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

# include <iostream>
# include <stdio.h>
# include <vector>
# include <string>
# include <queue>
# include <algorithm>
# include <limits>

#include <TMath.h>
#include <TColor.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TString.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THStack.h>
#include <TPaveText.h>

#include <RooChi2Var.h>
#include <RooConstVar.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooAbsArg.h>
#include <RooMinimizer.h>
#include <RooWorkspace.h>
#include <RooStats/ModelConfig.h>

#include <RooGaussian.h>
#include <RooBreitWigner.h>
#include <RooBifurGauss.h>
#include <RooPolynomial.h>
#include <RooCBShape.h>
#include <RooExponential.h>
#include <RooArgusBG.h>

#include <TCanvas.h>
#include <TCut.h>
#include <RooPlot.h>
#include <TAxis.h>

# define N_Needed_info 37
# define N_event_info 15
# define N_Upsilon_info 13
# define N_Bsig_info 7
# define N_Btag_info 7
# define N_decay 38 // five decay mode + others

void load_files(const char *dirname, std::vector<string>* names){
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
    int __event__;
    int __candidate__;
    int __ncandidates__;

    double Upsilon_decayID;
    double Bsig_decayID;

    int event_info[N_event_info];
    // 0: upsilon_experiment, 1: upsilon_run, 2: upsilon_event, 3: upsilon_candidate, 4: upsilon_ncandidates
    // 5: Bsig_experiment, 6: Bsig_run, 7: Bsig_event, 8: Bsig_candidate, 9: Bsig_ncandidates
    // 10: Btag_experiment, 11: Btag_run, 12: Btag_event, 13: Btag_candidate, 14: Btag_ncandidates

    double Upsilon_info[N_Upsilon_info];
    // 0: number of ECL clusters in ROE(cleanMask), 1: number of KLM clusters in ROE
    // 2: number of tracks in ROE(cleanMask), 3: roeEextra(cleanMask)
    // 4: nROE_NeutralECLClusters(cleanMask), 5: nROE_K_S0, 6: nROE_pi0
    // 7: missing momentum of event theta, 8: missing momentum
    // 9: missing energy at CMS, 10: number of remaining tracks
    // 11: roeNeextra(cleanMask), 12: useCMSFrame__boroeNeextra__bocleanMask__bc__bc

    double Bsig_info[N_Bsig_info];
    // 0: Bsig_E, 1: Bsig_E_CMS, 2: Bsig_E_Recoil
    // 3: Bsig_p, 4: Bsig_p_CMS, 5: Bsig_p_Recoil
    // 6: M

    double Btag_info[N_Btag_info];
    // 0: Btag_dmID, 1: Btag_Mbc, 2: Btag_deltaE
    // 3: Btag_E, 4: Btag_E_CMS, 5: Btag_signalprobability, 6: chiProb_tag

    double Needed_info[N_Needed_info];
    // 0: R2, 1: thrustBm, 2: thrustOm. 3: cosTBTO
    // 4: cosTBz, 5: KSFWVariables(et), 6: KSFWVariables(mm2)
    // 7: KSFWVariables(hso00), 8: KSFWVariables(hso01), 9: KSFWVariables(hso02)
    // 10: KSFWVariables(hso03), 11: KSFWVariables(hso04)
    // 12: KSFWVariables(hso10), 13: KSFWVariables(hso12), 14: KSFWVariables(hso14)
    // 15: KSFWVariables(hso20), 16: KSFWVariables(hso22), 17: KSFWVariables(hso24)
    // 18: KSFWVariables(hoo0), 19: KSFWVariables(hoo1), 20: KSFWVariables(hoo2)
    // 21: KSFWVariables(hoo3), 22: KSFWVariables(hoo4)
    // 23: CleoConeCS(1), 24: CleoConeCS(2), 25: CleoConeCS(3)
    // 26: CleoConeCS(4), 27: CleoConeCS(5), 28: CleoConeCS(6)
    // 29: CleoConeCS(7), 30: CleoConeCS(8), 31:CleoConeCS(9)
    // 32: missing mass^2, 33: visible energy, 34: theta_CMS_of_tag
    // 35: decayhash, 36: decayhash_extended

    double Decay[N_decay]; // MC level info

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

private:
    std::queue<Data> TotalData;
    std::vector<TH1F*> TH1Fs;
    std::vector<TH2F*> TH2Fs;
    int current_TH1F;
    int current_TH2F;

    std::vector<int> N_events;
    std::vector<int> N_candidates;
    std::vector<int> N_candidates_modes[Loader::MAX_NUM_DECAYMODE];
    std::vector<std::string> titles;
    int current_N_event;
    int current_N_candidate;

    std::vector<int> experiment_indices;
    std::vector<int> run_indices;
    std::vector<int> event_indices;
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

    std::vector<THStack*> THStacks;
    std::vector<TH1F*> TH1Fs_THStack[Loader::MAX_NUM_DECAYMODE];
    int current_THStack;

    int current_Confusion_matrix;
    double Confusion[Loader::MAX_NUM_DECAYMODE][Loader::MAX_NUM_DECAYMODE_MC]; // [reco][MC truth]
    double Confusion_square[Loader::MAX_NUM_DECAYMODE][Loader::MAX_NUM_DECAYMODE + 1]; // [reco][MC truth]
    bool Confusion_matrixIsOn;

    int EventDataToTree[N_event_info];
    double UpsilonDataToTree[N_Upsilon_info];
    double BsigDataToTree[N_Bsig_info];
    double BtagDataToTree[N_Btag_info];
    double DecayDataToTree[N_decay];
    double DataToTree[N_Needed_info];
    double Upsilon_decayIDToTree;
    double Bsig_decayIDToTree;

    bool DoesItHaveXsBranch;

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
    void PrintInformation(std::string title);
    void Cut(Loader::Variable variable, int i, Loader::Inequality inq, double value);
    void Cut(Loader::Variable variable, int i, Loader::Inequality inq, double value, Loader::Qualifier qualifier, Loader::DecayMode decaymode);
    void BCS(Loader::Variable variable, int index, Loader::BCS_criterion crit);
    bool IsBCSValid();
    void End();
    void PrintRootFile(std::string output_name);
    void PrintSeparateRootFile(std::string output_name);
    void ConvertIntoSeparateDataFile(std::string output_name, int flag = 0);
    void PrintDebugLogIf(Loader::Variable variable, int i, Loader::Inequality inq, double value);
    void PrintConfusionMatrix();
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
    AllOfThemHaveXsBranch = true;
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
    current_Confusion_matrix = 0;
    Confusion_matrixIsOn = false;
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

    Data temp = { 0 };

    // get event_info
    tree_upsilon->SetBranchAddress("__experiment__", &temp.event_info[0]);
    tree_upsilon->SetBranchAddress("__run__", &temp.event_info[1]);
    tree_upsilon->SetBranchAddress("__event__", &temp.event_info[2]);
    tree_upsilon->SetBranchAddress("__candidate__", &temp.event_info[3]);
    tree_upsilon->SetBranchAddress("__ncandidates__", &temp.event_info[4]);
    tree_Bsig->SetBranchAddress("__experiment__", &temp.event_info[5]);
    tree_Bsig->SetBranchAddress("__run__", &temp.event_info[6]);
    tree_Bsig->SetBranchAddress("__event__", &temp.event_info[7]);
    tree_Bsig->SetBranchAddress("__candidate__", &temp.event_info[8]);
    tree_Bsig->SetBranchAddress("__ncandidates__", &temp.event_info[9]);
    tree_Btag->SetBranchAddress("__experiment__", &temp.event_info[10]);
    tree_Btag->SetBranchAddress("__run__", &temp.event_info[11]);
    tree_Btag->SetBranchAddress("__event__", &temp.event_info[12]);
    tree_Btag->SetBranchAddress("__candidate__", &temp.event_info[13]);
    tree_Btag->SetBranchAddress("__ncandidates__", &temp.event_info[14]);

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
    tree_upsilon->SetBranchAddress("nRemainingTracksInEvent", &temp.Upsilon_info[10]);
    tree_upsilon->SetBranchAddress("roeNeextra__bocleanMask__bc", &temp.Upsilon_info[11]);
    tree_upsilon->SetBranchAddress("useCMSFrame__boroeNeextra__bocleanMask__bc__bc", &temp.Upsilon_info[12]);

    // get Bsig_info
    tree_Bsig->SetBranchAddress("Bsig_E", &temp.Bsig_info[0]);
    tree_Bsig->SetBranchAddress("Bsig_useCMSFrame_E", &temp.Bsig_info[1]);
    tree_upsilon->SetBranchAddress("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc", &temp.Bsig_info[2]);
    tree_Bsig->SetBranchAddress("Bsig_p", &temp.Bsig_info[3]);
    tree_Bsig->SetBranchAddress("Bsig_useCMSFrame_p", &temp.Bsig_info[4]);
    tree_upsilon->SetBranchAddress("useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc", &temp.Bsig_info[5]);
    tree_Bsig->SetBranchAddress("Bsig_M", &temp.Bsig_info[6]);

    // get Btag_info
    tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &temp.Btag_info[0]);
    tree_Btag->SetBranchAddress("Btag_Mbc", &temp.Btag_info[1]);
    tree_Btag->SetBranchAddress("Btag_deltaE", &temp.Btag_info[2]);
    tree_Btag->SetBranchAddress("Btag_E", &temp.Btag_info[3]);
    tree_Btag->SetBranchAddress("Btag_useCMSFrame_E", &temp.Btag_info[4]);
    tree_Btag->SetBranchAddress("Btag_extraInfo_SignalProbability", &temp.Btag_info[5]);
    tree_Btag->SetBranchAddress("Btag_chiProb", &temp.Btag_info[6]);

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
    }
    else {
        for (int i = 0; i < N_decay; i++) temp.Decay[i] = -1;
    }

    printf("%lld entries...\n", tree_upsilon->GetEntries());
    for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
        tree_upsilon->GetEntry(j);
        tree_Bsig->GetEntry(j);
        tree_Btag->GetEntry(j);
        if (DoesItHaveXsBranch) tree_Xs->GetEntry(j);
        TotalData.push(temp);
    }
    input_file->Close();
}

bool Loader::event_info_is_valid() {
    std::queue<Data> temp_queue;
    while (!TotalData.empty()) {
        Data temp_data = TotalData.front();
        TotalData.pop();
        if (temp_data.event_info[0] != temp_data.event_info[5] || temp_data.event_info[0] != temp_data.event_info[10] || temp_data.event_info[5] != temp_data.event_info[10]) return false;
        if (temp_data.event_info[1] != temp_data.event_info[6] || temp_data.event_info[1] != temp_data.event_info[11] || temp_data.event_info[6] != temp_data.event_info[11]) return false;
        if (temp_data.event_info[2] != temp_data.event_info[7] || temp_data.event_info[2] != temp_data.event_info[12] || temp_data.event_info[7] != temp_data.event_info[12]) return false;
        if (temp_data.event_info[3] != temp_data.event_info[8] || temp_data.event_info[3] != temp_data.event_info[13] || temp_data.event_info[8] != temp_data.event_info[13]) return false;
        if (temp_data.event_info[4] != temp_data.event_info[9] || temp_data.event_info[4] != temp_data.event_info[14] || temp_data.event_info[9] != temp_data.event_info[14]) return false;
        temp_data.__experiment__ = temp_data.event_info[0];
        temp_data.__run__ = temp_data.event_info[1];
        temp_data.__event__ = temp_data.event_info[2];
        temp_data.__candidate__ = temp_data.event_info[3];
        temp_data.__ncandidates__ = temp_data.event_info[4];
        temp_queue.push(temp_data);
    }
    TotalData.swap(temp_queue);
    return true;
}

void Loader::DrawTH1F(const char* name, const char* title, int nbins, double x_low, double x_high, Loader::Variable variable, int i, Loader::ValueOption dr = Loader::Linear) {
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

void Loader::DrawTH1F(const char* name, const char* title, int nbins, double x_low, double x_high, Loader::Variable variable, int i, Loader::Qualifier qualifier, Loader::DecayMode decaymode, Loader::ValueOption dr = Loader::Linear) {
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

void Loader::DrawTHStack(const char* name, const char* title, int nbins, double x_low, double x_high, Loader::Variable variable, int i, Loader::ValueOption dr = Loader::Linear) {
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

void Loader::PrintInformation(std::string title) {
    typedef struct labels {
        int __experiment__;
        int __run__;
        int __event__;
        int __ncandidates__;
    } Labels;
    std::vector<Labels> label_list;

    if (N_events.size() == current_N_event && N_candidates.size() == current_N_candidate && N_events.size() == N_candidates.size()) { // allocate new int
        N_events.push_back(0);
        N_candidates.push_back(0);
        titles.push_back(title);
        for (int i = 0; i < Loader::MAX_NUM_DECAYMODE; i++) N_candidates_modes[i].push_back(0);
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
            if (label_list.at(i).__experiment__ == temp.__experiment__ && label_list.at(i).__run__ == temp.__run__ && label_list.at(i).__event__ == temp.__event__ && label_list.at(i).__ncandidates__ == temp.__ncandidates__) { overlap = true; }
        }
        if (overlap == false) {
            N_events.at(current_N_event) = N_events.at(current_N_event) + 1;
            Labels temp_Labels;
            temp_Labels.__experiment__ = temp.__experiment__;
            temp_Labels.__run__ = temp.__run__;
            temp_Labels.__event__ = temp.__event__;
            temp_Labels.__ncandidates__ = temp.__ncandidates__;
            label_list.push_back(temp_Labels);
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
        N_candidates_modes[decaymodeid].at(current_N_candidate)++;

        TotalData.push(temp);
    }
    N_candidates.at(current_N_candidate) = N_candidates.at(current_N_candidate) + TotalData.size();

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
        int event_ = initial_data.__event__;
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

bool Loader::IsBCSValid() {
    bool IsItValid = true;

    typedef struct labels {
        int __experiment__;
        int __run__;
        int __event__;
        int __ncandidates__;
    } Labels;

    std::vector<Labels> label_list;
    std::queue<Data> TotalData_;
    TotalData_.swap(TotalData);

    while (!TotalData_.empty()) {
        Data temp = TotalData_.front();
        TotalData_.pop();
        for (unsigned int i = 0; i < label_list.size(); i++) {
            if (label_list.at(i).__experiment__ == temp.__experiment__ && label_list.at(i).__run__ == temp.__run__ && label_list.at(i).__event__ == temp.__event__ && label_list.at(i).__ncandidates__ == temp.__ncandidates__) IsItValid = false;
        }
        Labels temp_Labels;
        temp_Labels.__experiment__ = temp.__experiment__;
        temp_Labels.__run__ = temp.__run__;
        temp_Labels.__event__ = temp.__event__;
        temp_Labels.__ncandidates__ = temp.__ncandidates__;
        label_list.push_back(temp_Labels);

        TotalData.push(temp);
    }
    return IsItValid;
}

void Loader::End() {
    for (int i = 0; i < N_events.size();i++) {
        printf("%s\n", titles.at(i).c_str());
        printf("Number of event: %d\n", N_events.at(i));
        printf("Number of candidate: %d\n", N_candidates.at(i));
        for (int j = 0; j < Loader::MAX_NUM_DECAYMODE; j++) printf("Number of candidate of decayID %d: %d\n", j, N_candidates_modes[j].at(i));
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
            printf("event num: %d\n", event_indices.at(i));
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
        tree_upsilon->Branch("__experiment__", &EventDataToTree[0]);
        tree_upsilon->Branch("__run__", &EventDataToTree[1]);
        tree_upsilon->Branch("__event__", &EventDataToTree[2]);
        tree_upsilon->Branch("__candidate__", &EventDataToTree[3]);
        tree_upsilon->Branch("__ncandidates__", &EventDataToTree[4]);
        tree_Bsig->Branch("__experiment__", &EventDataToTree[5]);
        tree_Bsig->Branch("__run__", &EventDataToTree[6]);
        tree_Bsig->Branch("__event__", &EventDataToTree[7]);
        tree_Bsig->Branch("__candidate__", &EventDataToTree[8]);
        tree_Bsig->Branch("__ncandidates__", &EventDataToTree[9]);
        tree_Btag->Branch("__experiment__", &EventDataToTree[10]);
        tree_Btag->Branch("__run__", &EventDataToTree[11]);
        tree_Btag->Branch("__event__", &EventDataToTree[12]);
        tree_Btag->Branch("__candidate__", &EventDataToTree[13]);
        tree_Btag->Branch("__ncandidates__", &EventDataToTree[14]);

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

        // get Bsig_info
        tree_Bsig->Branch("Bsig_E", &BsigDataToTree[0]);
        tree_Bsig->Branch("Bsig_useCMSFrame_E", &BsigDataToTree[1]);
        tree_upsilon->Branch("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc", &BsigDataToTree[2]);
        tree_Bsig->Branch("Bsig_p", &BsigDataToTree[3]);
        tree_Bsig->Branch("Bsig_useCMSFrame_p", &BsigDataToTree[4]);
        tree_upsilon->Branch("useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc", &BsigDataToTree[5]);
        tree_Bsig->Branch("Bsig_M", &BsigDataToTree[6]);

        // get Btag_info
        tree_Btag->Branch("Btag_extraInfo_decayModeID", &BtagDataToTree[0]);
        tree_Btag->Branch("Btag_Mbc", &BtagDataToTree[1]);
        tree_Btag->Branch("Btag_deltaE", &BtagDataToTree[2]);
        tree_Btag->Branch("Btag_E", &BtagDataToTree[3]);
        tree_Btag->Branch("Btag_useCMSFrame_E", &BtagDataToTree[4]);
        tree_Btag->Branch("Btag_extraInfo_SignalProbability", &BtagDataToTree[5]);
        tree_Btag->Branch("Btag_chiProb", &BtagDataToTree[6]);

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

        for (int i = 0; i < N_event_info; i++) {
            EventDataToTree[i] = temp.event_info[i];
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
        for (int i = 0; i < N_Needed_info; i++) {
            DataToTree[i] = temp.Needed_info[i];
        }
        Upsilon_decayIDToTree = temp.Upsilon_decayID;
        Bsig_decayIDToTree = temp.Bsig_decayID;

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

    int temp_EventDataToTree[N_event_info];
    double temp_UpsilonDataToTree[N_Upsilon_info];
    double temp_BsigDataToTree[N_Bsig_info];
    double temp_BtagDataToTree[N_Btag_info];
    double temp_DecayDataToTree[N_decay];
    double temp_DataToTree[N_Needed_info];
    double temp_Upsilon_decayIDToTree;
    double temp_Bsig_decayIDToTree;

    /*================================================================*/
    // get event_info
    temp_tree_upsilon->Branch("__experiment__", &temp_EventDataToTree[0]);
    temp_tree_upsilon->Branch("__run__", &temp_EventDataToTree[1]);
    temp_tree_upsilon->Branch("__event__", &temp_EventDataToTree[2]);
    temp_tree_upsilon->Branch("__candidate__", &temp_EventDataToTree[3]);
    temp_tree_upsilon->Branch("__ncandidates__", &temp_EventDataToTree[4]);
    temp_tree_Bsig->Branch("__experiment__", &temp_EventDataToTree[5]);
    temp_tree_Bsig->Branch("__run__", &temp_EventDataToTree[6]);
    temp_tree_Bsig->Branch("__event__", &temp_EventDataToTree[7]);
    temp_tree_Bsig->Branch("__candidate__", &temp_EventDataToTree[8]);
    temp_tree_Bsig->Branch("__ncandidates__", &temp_EventDataToTree[9]);
    temp_tree_Btag->Branch("__experiment__", &temp_EventDataToTree[10]);
    temp_tree_Btag->Branch("__run__", &temp_EventDataToTree[11]);
    temp_tree_Btag->Branch("__event__", &temp_EventDataToTree[12]);
    temp_tree_Btag->Branch("__candidate__", &temp_EventDataToTree[13]);
    temp_tree_Btag->Branch("__ncandidates__", &temp_EventDataToTree[14]);

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

    // get Bsig_info
    temp_tree_Bsig->Branch("Bsig_E", &temp_BsigDataToTree[0]);
    temp_tree_Bsig->Branch("Bsig_useCMSFrame_E", &temp_BsigDataToTree[1]);
    temp_tree_upsilon->Branch("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc", &temp_BsigDataToTree[2]);
    temp_tree_Bsig->Branch("Bsig_p", &temp_BsigDataToTree[3]);
    temp_tree_Bsig->Branch("Bsig_useCMSFrame_p", &temp_BsigDataToTree[4]);
    temp_tree_upsilon->Branch("useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc", &temp_BsigDataToTree[5]);
    temp_tree_Bsig->Branch("Bsig_M", &temp_BsigDataToTree[6]);

    // get Btag_info
    temp_tree_Btag->Branch("Btag_extraInfo_decayModeID", &temp_BtagDataToTree[0]);
    temp_tree_Btag->Branch("Btag_Mbc", &temp_BtagDataToTree[1]);
    temp_tree_Btag->Branch("Btag_deltaE", &temp_BtagDataToTree[2]);
    temp_tree_Btag->Branch("Btag_E", &temp_BtagDataToTree[3]);
    temp_tree_Btag->Branch("Btag_useCMSFrame_E", &temp_BtagDataToTree[4]);
    temp_tree_Btag->Branch("Btag_extraInfo_SignalProbability", &temp_BtagDataToTree[5]);
    temp_tree_Btag->Branch("Btag_chiProb", &temp_BtagDataToTree[6]);

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
    }
    else {
        for (int i = 0; i < N_decay; i++)  temp_DecayDataToTree[i] = -1;
    }
    /*================================================================*/

    std::queue<Data> temp_queue;
    temp_queue.swap(TotalData);
    while (!temp_queue.empty()) {
        Data temp = temp_queue.front();
        temp_queue.pop();

        for (int i = 0; i < N_event_info; i++) {
            temp_EventDataToTree[i] = temp.event_info[i];
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
        for (int i = 0; i < N_Needed_info; i++) {
            temp_DataToTree[i] = temp.Needed_info[i];
        }
        temp_Upsilon_decayIDToTree = temp.Upsilon_decayID;
        temp_Bsig_decayIDToTree = temp.Bsig_decayID;

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

void Loader::ConvertIntoSeparateDataFile(std::string output_name, int flag = 0) {

    TFile* temp_file = new TFile(output_name.c_str(), "recreate");
    temp_file->cd();
    TTree* temp_tree = new TTree("data", "");

    int temp_EventDataToTree[N_event_info / 3];
    double temp_UpsilonDataToTree[N_Upsilon_info];
    double temp_BsigDataToTree[N_Bsig_info];
    double temp_BtagDataToTree[N_Btag_info];
    double temp_DataToTree[N_Needed_info];
    double temp_DecayDataToTree[N_decay];
    double temp_Upsilon_decayIDToTree;
    double temp_Bsig_decayIDToTree;
    int temp_flag;

    /*================================================================*/
    // get event_info
    temp_tree->Branch("__experiment__", &temp_EventDataToTree[0]);
    temp_tree->Branch("__run__", &temp_EventDataToTree[1]);
    temp_tree->Branch("__event__", &temp_EventDataToTree[2]);
    temp_tree->Branch("__candidate__", &temp_EventDataToTree[3]);
    temp_tree->Branch("__ncandidates__", &temp_EventDataToTree[4]);

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

    // get Bsig_info
    temp_tree->Branch("Bsig_E", &temp_BsigDataToTree[0]);
    temp_tree->Branch("Bsig_useCMSFrame_E", &temp_BsigDataToTree[1]);
    temp_tree->Branch("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc", &temp_BsigDataToTree[2]);
    temp_tree->Branch("Bsig_p", &temp_BsigDataToTree[3]);
    temp_tree->Branch("Bsig_useCMSFrame_p", &temp_BsigDataToTree[4]);
    temp_tree->Branch("useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc", &temp_BsigDataToTree[5]);
    temp_tree->Branch("Bsig_M", &temp_BsigDataToTree[6]);

    // get Btag_info
    temp_tree->Branch("Btag_extraInfo_decayModeID", &temp_BtagDataToTree[0]);
    temp_tree->Branch("Btag_Mbc", &temp_BtagDataToTree[1]);
    temp_tree->Branch("Btag_deltaE", &temp_BtagDataToTree[2]);
    temp_tree->Branch("Btag_E", &temp_BtagDataToTree[3]);
    temp_tree->Branch("Btag_useCMSFrame_E", &temp_BtagDataToTree[4]);
    temp_tree->Branch("Btag_extraInfo_SignalProbability", &temp_BtagDataToTree[5]);
    temp_tree->Branch("Btag_chiProb", &temp_BtagDataToTree[6]);

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
    }
    else {
        for (int i = 0; i < N_decay; i++)  temp_DecayDataToTree[i] = -1;
    }

    // flag
    temp_tree->Branch("flag", &temp_flag);
    /*================================================================*/

    std::queue<Data> temp_queue;
    temp_queue.swap(TotalData);
    while (!temp_queue.empty()) {
        Data temp = temp_queue.front();
        temp_queue.pop();

        for (int i = 0; i < N_event_info / 3; i++) {
            temp_EventDataToTree[i] = temp.event_info[i];
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
        temp_Upsilon_decayIDToTree = temp.Upsilon_decayID;
        temp_Bsig_decayIDToTree = temp.Bsig_decayID;
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


void Loader::PrintConfusionMatrix() {
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

        Confusion[decaymodeid][decaymodeid_MC]++;

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
        else if (decaymodeid_MC == 15)decaymodeid_MC_for_square = 11;
        else if (decaymodeid_MC == 16)decaymodeid_MC_for_square = 12;
        else if (decaymodeid_MC == 17)decaymodeid_MC_for_square = 13;
        else if (decaymodeid_MC == 18)decaymodeid_MC_for_square = 14;
        else if (decaymodeid_MC == 19)decaymodeid_MC_for_square = 15;
        else if (decaymodeid_MC == 20)decaymodeid_MC_for_square = 16;
        else if (decaymodeid_MC == 21)decaymodeid_MC_for_square = 17;
        else if (decaymodeid_MC == 22)decaymodeid_MC_for_square = 18;
        else if (decaymodeid_MC == 23)decaymodeid_MC_for_square = 19;
        else if (decaymodeid_MC == 24)decaymodeid_MC_for_square = 20;
        Confusion_square[decaymodeid][decaymodeid_MC_for_square]++;

        TotalData.push(temp);
    }

    Confusion_matrixIsOn = true;
    current_Confusion_matrix++;
}

void main(int argc, char* argv[]) { // argv[1]: Ntuple input with path, argv[2]: destination of output, 

    std::string filename_with_path = std::string(argv[1]); // /home/file.root
    std::string destination = std::string(argv[2]);

    std::string base_filename = filename_with_path.substr(filename_with_path.find_last_of("/") + 1); // file.root
    std::string::size_type const p(base_filename.find_last_of('.'));
    std::string file_without_extension = base_filename.substr(0, p); // file

    char endch = destination.back();
    if (endch != '/') destination = destination + std::string("/");

    Loader loader;

    loader.initialize();

    TFile* input_file = new TFile(filename_with_path.c_str(), "read");
    printf("%s\n", filename_with_path.c_str());
    loader.GetData(input_file);
    if (loader.event_info_is_valid() == false) { printf("error!\n"); return; }

    loader.PrintInformation(std::string("========== inital =========="));

    //loader.DrawTH1F("Btag_chiProb", "chiProb of Btag;chiProb;evt", 100, -1.2, 1.2, Loader::Btag, 6);
    loader.Cut(Loader::Btag, 6, Loader::larger_than, -0.5);
    loader.PrintInformation(std::string("========== chiProb_Btag > -0.5 =========="));

    loader.PrintSeparateRootFile(destination + file_without_extension + std::string("_before_nRawtrack_cut.root"));
    //loader.DrawTH1F("nROE_track_Upsilon_after_initial", "number of raw tracks in ROE of #Upsilon(4S);number of raw tracks;evt", 100, -0.5, 13.5, Loader::Upsilon, 10);
    loader.Cut(Loader::Upsilon, 10, Loader::smaller_than, 0.5);
    loader.PrintInformation(std::string("========== nRawtrack = 0 =========="));

    loader.PrintSeparateRootFile(destination + file_without_extension + std::string("_before_npi0_cut.root"));
    //loader.DrawTH1F("nROE_pi0_after_ntrack_cut", "number of #pi^{0} candidates in ROE of #Upsilon(4S);number of #pi^{0} candidates;evt", 100, -0.5, 13.5, Loader::Upsilon, 6);
    loader.Cut(Loader::Upsilon, 6, Loader::smaller_than, 0.5);
    loader.PrintInformation(std::string("========== npi0 = 0 =========="));

    loader.PrintSeparateRootFile(destination + file_without_extension + std::string("_before_missing_momentum_theta_cut.root"));
    //loader.DrawTH1F("missing_momentum_theta_after_npi0_cut", "#theta_{missing};#theta_{missing} [rad];evt", 100, 0, 3.2, Loader::Upsilon, 7);
    loader.Cut(Loader::Upsilon, 7, Loader::smaller_than, 2.618);
    loader.Cut(Loader::Upsilon, 7, Loader::larger_than, 0.297);
    loader.PrintInformation(std::string("========== 0.297 < missing momentum theta < 2.618 =========="));

    loader.PrintSeparateRootFile(destination + file_without_extension + std::string("_before_psig_cut.root"));
    //loader.DrawTH1F("momentum_Bsig_after_missing_theta_cut", "momentum of B_{sig} at CMS;momentum [GeV];evt", 100, 0, 3.2, Loader::Bsig, 4);
    loader.Cut(Loader::Bsig, 4, Loader::smaller_than, 2.96);
    loader.Cut(Loader::Bsig, 4, Loader::larger_than, 0.5);
    loader.PrintInformation(std::string("========== 0.5 < momentum of signal side < 2.96 =========="));

    loader.PrintSeparateRootFile(destination + file_without_extension + std::string("_before_delE_cut.root"));
    //loader.DrawTH2F("MbcVSdeltaE_after_psig_cut", ";Mbc of B_{tag} [GeV];#DeltaE of B_{tag} [GeV]", 100, 5.24, 5.3, 100, -0.2, 0.2, Loader::Btag, 1, Loader::Btag, 2);
    loader.Cut(Loader::Btag, 2, Loader::larger_than, -0.1);
    loader.Cut(Loader::Btag, 2, Loader::smaller_than, 0.1);
    loader.PrintInformation(std::string("========== abs(deltaE) < 0.1 =========="));
    loader.DrawTH2F("MbcVSdeltaE_after_deltaE_strict_cut", ";Mbc of B_{tag} [GeV];#DeltaE of B_{tag} [GeV]", 100, 5.24, 5.3, 100, -0.2, 0.2, Loader::Btag, 1, Loader::Btag, 2);

    //loader.PrintSeparateRootFile(destination + file_without_extension + std::string("_before_Mbc_cut.root"));
    //loader.Cut(Loader::Btag, 1, Loader::larger_than, 5.27);
    //loader.PrintInformation(std::string("========== Mbc > 5.27 =========="));
    //loader.DrawTH2F("MbcVSdeltaE_after_Mbc_strict_cut", ";Mbc of B_{tag} [GeV];#DeltaE of B_{tag} [GeV]", 100, 5.24, 5.3, 100, -0.2, 0.2, Loader::Btag, 1, Loader::Btag, 2);

    loader.PrintSeparateRootFile(destination + file_without_extension + std::string("_before_Eecl_cut.root"));
    //loader.DrawTH1F("ROE_Eecl_Upsilon_after_Mbc_strict_cut", "E_ecl in ROE of #Upsilon(4S);E_{ecl} [GeV];candidates", 100, -0.1, 8, Loader::Upsilon, 3);
    loader.Cut(Loader::Upsilon, 3, Loader::smaller_than, 2);
    loader.PrintInformation(std::string("========== E_ecl < 2 GeV =========="));

    //loader.DrawTH1F("SignalProbability_Btag_before_BCS", "SignalProbability of B_{tag};log_{10}(SignalProbability);Num of candidate", 100, -10, 0, Loader::Btag, 5, Loader::Log);
    loader.BCS(Loader::Btag, 5, Loader::Highest);
    if (loader.IsBCSValid() == false) {
        printf("ERROR! it is not valid\n");
        exit(1);
    }
    //loader.DrawTH1F("SignalProbability_Btag_after_BCS", "SignalProbability of B_{tag};log_{10}(SignalProbability);Num of candidate", 100, -10, 0, Loader::Btag, 5, Loader::Log);
    loader.PrintInformation(std::string("========== BCS =========="));

    //loader.DrawTH1F("nROE_ECLcluster_Upsilon", "number of ECL clusters in ROE of #Upsilon(4S);number of ECL clusters;evt", 14, -0.5, 13.5, Loader::Upsilon, 0);
    //loader.DrawTH1F("nROE_KLMcluster_Upsilon", "number of KLM clusters in ROE of #Upsilon(4S);number of KLM clusters;evt", 14, -0.5, 13.5, Loader::Upsilon, 1);
    //loader.DrawTH1F("Bsig_p_LAB", "momentum of B_{sig} at LAB frame;p [GeV];evt", 50, -0.5, 6, Loader::Bsig, 3);
    //loader.DrawTH1F("Bsig_p_CMS", "momentum of B_{sig} at CMS frame;p [GeV];evt", 50, -0.5, 6, Loader::Bsig, 4);
    //loader.DrawTH1F("Bsig_p_RecoilRest", "momentum of B_{sig} at rest frame of recoil system;p [GeV];evt", 50, -0.5, 6, Loader::Bsig, 5);
    //loader.DrawTH1F("Btag_dmID", "decay ID of B_{tag};decay ID;evt", 74, -0.5, 36.5, Loader::Btag, 0);
    //loader.DrawTH1F("nROE_K_S0", "number of K_S0 candidates in ROE of #Upsilon(4S);number of K_{S}^{0} candidates in ROE;evt", 100, -0.5, 5.5, Loader::Upsilon, 5);
    //loader.DrawTH1F("theta_missing_momentum", "#theta of missing momentum;#theta [rad];evt", 50, 0, 3.2, Loader::Upsilon, 7);
    //loader.DrawTH1F("M_Xs", "mass of X_{s};M_{Xs} [GeV];evt", 100, 0, 3.5, Loader::Bsig, 6);
    //loader.DrawTH1F("ROE_Eecl_Upsilon_final", "E_ecl in ROE of #Upsilon(4S);E_{ecl} [GeV];candidates", 100, -0.1, 8, Loader::Upsilon, 3);

    loader.PrintSeparateRootFile(destination + file_without_extension + std::string("_final_output.root"));
    //loader.PrintRootFile(file_without_extension + std::string("_final_output_merge.root"));
    //loader.ConvertIntoSeparateDataFile("data.root");
    //loader.PrintConfusionMatrix();

    loader.End();
}
