// last update: 2021-09-21-00
// for Belle2 data

# define N_Needed_info 29
# define N_event_info 15
# define N_Upsilon_info 11
# define N_Bsig_info 7
# define N_Btag_info 6

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

    double Bsig_info[N_Bsig_info];
    // 0: Bsig_E, 1: Bsig_E_CMS, 2: Bsig_E_Recoil
    // 3: Bsig_p, 4: Bsig_p_CMS, 5: Bsig_p_Recoil
    // 6: M

    double Btag_info[N_Btag_info];
    // 0: Btag_dmID, 1: Btag_Mbc, 2: Btag_deltaE
    // 3: Btag_E, 4: Btag_E_CMS, 5: Btag_signalprobability

    double Needed_info[N_Needed_info];
    // 0: R2, 1: thrustBm, 2: thrustOm. 3: cosTBTO
    // 4: cosTBz, 5: KSFWVariables(et), 6: KSFWVariables(mm2)
    // 7: KSFWVariables(hso00), 8: KSFWVariables(hso02), 9: KSFWVariables(hso04)
    // 10: KSFWVariables(hso10), 11: KSFWVariables(hso12), 12: KSFWVariables(hso14)
    // 13: KSFWVariables(hso20), 14: KSFWVariables(hso22), 15: KSFWVariables(hso24)
    // 16: CleoConeCS(1), 17: CleoConeCS(2), 18: CleoConeCS(3)
    // 19: CleoConeCS(4), 20: CleoConeCS(5), 21: CleoConeCS(6)
    // 22: CleoConeCS(7), 23: CleoConeCS(8), 24:CleoConeCS(9)
    // 25: missing mass^2, 26: visible energy, 27: theta_CMS_of_tag
    // 28: chiProb_tag

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
    enum DecayMode {
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
    int current_file;

    std::vector<THStack*> THStacks;
    std::vector<TH1F*> TH1Fs_THStack[Loader::MAX_NUM_DECAYMODE];
    int current_THStack;

    int EventDataToTree[N_event_info];
    double UpsilonDataToTree[N_Upsilon_info];
    double BsigDataToTree[N_Bsig_info];
    double BtagDataToTree[N_Btag_info];
    double DataToTree[N_Needed_info];
    double Upsilon_decayIDToTree;
    double Bsig_decayIDToTree;

    bool TrueIfDecayModeMatch(Data temp_data, Loader::DecayMode decaymode);

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
    void PrintDebugLogIf(Loader::Variable variable, int i, Loader::Inequality inq, double value);
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
}

void Loader::GetData(TFile* input_file) {
    TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
    TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
    TTree* tree_Btag = (TTree*)input_file->Get("Btag");

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

    // other information I need
    tree_Btag->SetBranchAddress("Btag_R2", &temp.Needed_info[0]);
    tree_Btag->SetBranchAddress("Btag_thrustBm", &temp.Needed_info[1]);
    tree_Btag->SetBranchAddress("Btag_thrustOm", &temp.Needed_info[2]);
    tree_Btag->SetBranchAddress("Btag_cosTBTO", &temp.Needed_info[3]);
    tree_Btag->SetBranchAddress("Btag_cosTBz", &temp.Needed_info[4]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_et", &temp.Needed_info[5]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_mm2", &temp.Needed_info[6]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hso00", &temp.Needed_info[7]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hso02", &temp.Needed_info[8]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hso04", &temp.Needed_info[9]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hso10", &temp.Needed_info[10]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hso12", &temp.Needed_info[11]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hso14", &temp.Needed_info[12]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hso20", &temp.Needed_info[13]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hso22", &temp.Needed_info[14]);
    tree_Btag->SetBranchAddress("Btag_KSFWVariables_hso24", &temp.Needed_info[15]);
    tree_Btag->SetBranchAddress("Btag_CleoConeCS_1", &temp.Needed_info[16]);
    tree_Btag->SetBranchAddress("Btag_CleoConeCS_2", &temp.Needed_info[17]);
    tree_Btag->SetBranchAddress("Btag_CleoConeCS_3", &temp.Needed_info[18]);
    tree_Btag->SetBranchAddress("Btag_CleoConeCS_4", &temp.Needed_info[19]);
    tree_Btag->SetBranchAddress("Btag_CleoConeCS_5", &temp.Needed_info[20]);
    tree_Btag->SetBranchAddress("Btag_CleoConeCS_6", &temp.Needed_info[21]);
    tree_Btag->SetBranchAddress("Btag_CleoConeCS_7", &temp.Needed_info[22]);
    tree_Btag->SetBranchAddress("Btag_CleoConeCS_8", &temp.Needed_info[23]);
    tree_Btag->SetBranchAddress("Btag_CleoConeCS_9", &temp.Needed_info[24]);
    tree_upsilon->SetBranchAddress("missingMass2OfEvent", &temp.Needed_info[25]);
    tree_upsilon->SetBranchAddress("visibleEnergyOfEventCMS", &temp.Needed_info[26]);
    tree_Btag->SetBranchAddress("Btag_useCMSFrame_theta", &temp.Needed_info[27]);
    tree_Btag->SetBranchAddress("Btag_chiProb", &temp.Needed_info[28]);

    printf("%lld entries...\n", tree_upsilon->GetEntries());
    for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
        tree_upsilon->GetEntry(j);
        tree_Bsig->GetEntry(j);
        tree_Btag->GetEntry(j);
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
    TotalData = temp_queue;
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
        printf("ERROR!\n");
        exit(1);
    }

    TH1F* temp_hist = TH1Fs.at(current_TH1F);
    std::queue<Data> temp_queue = TotalData;
    while (!temp_queue.empty()) {
        Data temp_data = temp_queue.front();
        temp_queue.pop();
        if (variable == Loader::Upsilon) {
            if(dr == Loader::Linear) temp_hist->Fill(temp_data.Upsilon_info[i]);
            else if (dr == Loader::Log) temp_hist->Fill(TMath::Log10(temp_data.Upsilon_info[i]));
            else { printf("ERROR!\n"); exit(1); }
        }
        else if(variable == Loader::Bsig) {
            if (dr == Loader::Linear) temp_hist->Fill(temp_data.Bsig_info[i]);
            else if (dr == Loader::Log) temp_hist->Fill(TMath::Log10(temp_data.Bsig_info[i]));
            else { printf("ERROR!\n"); exit(1); }
        }
        else if(variable == Loader::Btag) {
            if (dr == Loader::Linear) temp_hist->Fill(temp_data.Btag_info[i]);
            else if (dr == Loader::Log) temp_hist->Fill(TMath::Log10(temp_data.Btag_info[i]));
            else { printf("ERROR!\n"); exit(1); }
        }
        else {
            printf("ERROR!\n");
            exit(1);
        }
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
        printf("ERROR!\n");
        exit(1);
    }

    TH1F* temp_hist = TH1Fs.at(current_TH1F);
    std::queue<Data> temp_queue = TotalData;
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
            printf("ERROR!\n");
            exit(1);
        }

        if (variable == Loader::Upsilon) {
            if (dr == Loader::Linear) temp_hist->Fill(temp_data.Upsilon_info[i]);
            else if (dr == Loader::Log) temp_hist->Fill(TMath::Log10(temp_data.Upsilon_info[i]));
            else { printf("ERROR!\n"); exit(1); }
        }
        else if (variable == Loader::Bsig) {
            if (dr == Loader::Linear) temp_hist->Fill(temp_data.Bsig_info[i]);
            else if (dr == Loader::Log) temp_hist->Fill(TMath::Log10(temp_data.Bsig_info[i]));
            else { printf("ERROR!\n"); exit(1); }
        }
        else if (variable == Loader::Btag) {
            if (dr == Loader::Linear) temp_hist->Fill(temp_data.Btag_info[i]);
            else if (dr == Loader::Log) temp_hist->Fill(TMath::Log10(temp_data.Btag_info[i]));
            else { printf("ERROR!\n"); exit(1); }
        }
        else {
            printf("ERROR!\n");
            exit(1);
        }
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
        printf("ERROR!\n");
        exit(1);
    }

    TH2F* temp_hist = TH2Fs.at(current_TH2F);
    std::queue<Data> temp_queue = TotalData;
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
                printf("ERROR!\n");
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
                printf("ERROR!\n");
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
                printf("ERROR!\n");
                exit(1);
            }
        }
        else {
            printf("ERROR!\n");
            exit(1);
        }
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
        printf("ERROR!\n");
        exit(1);
    }

    TH2F* temp_hist = TH2Fs.at(current_TH2F);
    std::queue<Data> temp_queue = TotalData;
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
            printf("ERROR!\n");
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
                printf("ERROR!\n");
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
                printf("ERROR!\n");
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
                printf("ERROR!\n");
                exit(1);
            }
        }
        else {
            printf("ERROR!\n");
            exit(1);
        }
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
        printf("ERROR!\n");
        exit(1);
    }

    TH1F* temp_hist[Loader::MAX_NUM_DECAYMODE];
    for (int i = 0; i < MAX_NUM_DECAYMODE;i++) {
        temp_hist[i] = TH1Fs_THStack[i].at(current_THStack);
    }

    std::queue<Data> temp_queue = TotalData;
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
            printf("ERROR!\n");
            exit(1);
        }

        if (variable == Loader::Upsilon) {
            if (dr == Loader::Linear) temp_hist[decaymodeid]->Fill(temp_data.Upsilon_info[i]);
            else if (dr == Loader::Log) temp_hist[decaymodeid]->Fill(TMath::Log10(temp_data.Upsilon_info[i]));
            else { printf("ERROR!\n"); exit(1); }
        }
        else if (variable == Loader::Bsig) {
            if (dr == Loader::Linear) temp_hist[decaymodeid]->Fill(temp_data.Bsig_info[i]);
            else if (dr == Loader::Log) temp_hist[decaymodeid]->Fill(TMath::Log10(temp_data.Bsig_info[i]));
            else { printf("ERROR!\n"); exit(1); }
        }
        else if (variable == Loader::Btag) {
            if (dr == Loader::Linear) temp_hist[decaymodeid]->Fill(temp_data.Btag_info[i]);
            else if (dr == Loader::Log) temp_hist[decaymodeid]->Fill(TMath::Log10(temp_data.Btag_info[i]));
            else { printf("ERROR!\n"); exit(1); }
        }
        else {
            printf("ERROR!\n");
            exit(1);
        }
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
        printf("ERROR!\n");
        exit(1);
    }

    std::queue<Data> temp_queue = TotalData;
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
            printf("ERROR!\n");
            exit(1);
        }
        N_candidates_modes[decaymodeid].at(current_N_candidate)++;

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
            printf("ERROR!\n");
            exit(1);
        }
    }
    TotalData = temp_queue;
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
            printf("ERROR!\n");
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
            printf("ERROR!\n");
            exit(1);
        }
    }
    TotalData = temp_queue;
}

void Loader::PrintDebugLogIf(Loader::Variable variable, int i, Loader::Inequality inq, double value) {
    if (current_N_experiment_index > 0) { // allocate new int
        printf("The number of PrintDebugLogIf should not be larger than 1\n");
        printf("Only first PrintDebugLogIf is accepted\n");
        return;
    }

    std::queue<Data> temp_queue;
    temp_queue = TotalData;

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
            printf("ERROR!\n");
            exit(1);
        }
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
            printf("ERROR!\n");
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
    typedef struct labels {
        int __experiment__;
        int __run__;
        int __event__;
        int __ncandidates__;
    } Labels;

    std::vector<Labels> label_list;
    std::queue<Data> TotalData_;

    TotalData_ = TotalData;
    while (!TotalData_.empty()) {
        Data temp = TotalData_.front();
        TotalData_.pop();
        for (unsigned int i = 0; i < label_list.size(); i++) {
            if (label_list.at(i).__experiment__ == temp.__experiment__ && label_list.at(i).__run__ == temp.__run__ && label_list.at(i).__event__ == temp.__event__ && label_list.at(i).__ncandidates__ == temp.__ncandidates__) return false;
        }
        Labels temp_Labels;
        temp_Labels.__experiment__ = temp.__experiment__;
        temp_Labels.__run__ = temp.__run__;
        temp_Labels.__event__ = temp.__event__;
        temp_Labels.__ncandidates__ = temp.__ncandidates__;
        label_list.push_back(temp_Labels);
    }
    return true;
}

void Loader::End() {
    for (int i = 0; i < N_events.size();i++) {
        printf("%s\n", titles.at(i).c_str());
        printf("Number of event: %d\n", N_events.at(i));
        printf("Number of candidate: %d\n", N_candidates.at(i));
        for (int j = 0; j < Loader::MAX_NUM_DECAYMODE; j++) printf("Number of candidate of decayID %d: %d\n", j, N_candidates_modes[j].at(i));
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
        temp_file->cd();
        temp_tree_upsilon->Write();
        temp_tree_Bsig->Write();
        temp_tree_Btag->Write();
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
    if (files.size() == current_file && trees_upsilon.size() == current_file && trees_Bsig.size() == current_file && trees_Btag.size() == current_file) { // allocate new TFile and TTree
        TFile* file = new TFile(output_name.c_str(), "recreate");
        file->cd();
        TTree* tree_upsilon = new TTree("Upsilon", "");
        TTree* tree_Bsig = new TTree("Bsig", "");
        TTree* tree_Btag = new TTree("Btag", "");

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
        tree_upsilon->Branch("extraInfo__bodecayModeID__bc", &temp.Upsilon_decayIDToTree);
        tree_Bsig->Branch("Bsig_daughter_0_extraInfo_decayModeID", &temp.Bsig_decayIDToTree);

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

        // other information I need
        tree_Btag->Branch("Btag_R2", &DataToTree[0]);
        tree_Btag->Branch("Btag_thrustBm", &DataToTree[1]);
        tree_Btag->Branch("Btag_thrustOm", &DataToTree[2]);
        tree_Btag->Branch("Btag_cosTBTO", &DataToTree[3]);
        tree_Btag->Branch("Btag_cosTBz", &DataToTree[4]);
        tree_Btag->Branch("Btag_KSFWVariables_et", &DataToTree[5]);
        tree_Btag->Branch("Btag_KSFWVariables_mm2", &DataToTree[6]);
        tree_Btag->Branch("Btag_KSFWVariables_hso00", &DataToTree[7]);
        tree_Btag->Branch("Btag_KSFWVariables_hso02", &DataToTree[8]);
        tree_Btag->Branch("Btag_KSFWVariables_hso04", &DataToTree[9]);
        tree_Btag->Branch("Btag_KSFWVariables_hso10", &DataToTree[10]);
        tree_Btag->Branch("Btag_KSFWVariables_hso12", &DataToTree[11]);
        tree_Btag->Branch("Btag_KSFWVariables_hso14", &DataToTree[12]);
        tree_Btag->Branch("Btag_KSFWVariables_hso20", &DataToTree[13]);
        tree_Btag->Branch("Btag_KSFWVariables_hso22", &DataToTree[14]);
        tree_Btag->Branch("Btag_KSFWVariables_hso24", &DataToTree[15]);
        tree_Btag->Branch("Btag_CleoConeCS_1", &DataToTree[16]);
        tree_Btag->Branch("Btag_CleoConeCS_2", &DataToTree[17]);
        tree_Btag->Branch("Btag_CleoConeCS_3", &DataToTree[18]);
        tree_Btag->Branch("Btag_CleoConeCS_4", &DataToTree[19]);
        tree_Btag->Branch("Btag_CleoConeCS_5", &DataToTree[20]);
        tree_Btag->Branch("Btag_CleoConeCS_6", &DataToTree[21]);
        tree_Btag->Branch("Btag_CleoConeCS_7", &DataToTree[22]);
        tree_Btag->Branch("Btag_CleoConeCS_8", &DataToTree[23]);
        tree_Btag->Branch("Btag_CleoConeCS_9", &DataToTree[24]);
        tree_upsilon->Branch("missingMass2OfEvent", &DataToTree[25]);
        tree_upsilon->Branch("visibleEnergyOfEventCMS", &DataToTree[26]);
        tree_Btag->Branch("Btag_useCMSFrame_theta", &DataToTree[27]);
        tree_Btag->Branch("Btag_chiProb", &DataToTree[28]);
        /*================================================================*/
        files.push_back(file);
        trees_upsilon.push_back(tree_upsilon);
        trees_Bsig.push_back(tree_Bsig);
        trees_Btag.push_back(tree_Btag);
    }
    else if (files.size() > current_file && trees_upsilon.size() > current_file && files.size() == trees_upsilon.size() && files.size() == trees_Bsig.size() && files.size() == trees_Btag.size()) { // use what I have
    }
    else { // error
        printf("ERROR!\n");
        exit(1);
    }

    TFile* temp_file = files.at(current_file);
    temp_file->cd();
    TTree* temp_tree_upsilon = trees_upsilon.at(current_file);
    TTree* temp_tree_Bsig = trees_Bsig.at(current_file);
    TTree* temp_tree_Btag = trees_Btag.at(current_file);
    std::queue<Data> temp_queue = TotalData;
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
        for (int i = 0; i < N_Needed_info; i++) {
            DataToTree[i] = temp.Needed_info[i];
        }
        Upsilon_decayIDToTree = temp.Upsilon_decayID;
        Bsig_decayIDToTree = temp.Bsig_decayID;

        temp_tree_upsilon->Fill();
        temp_tree_Bsig->Fill();
        temp_tree_Btag->Fill();
    }

    current_file++;
}

bool Loader::TrueIfDecayModeMatch(Data temp_data, Loader::DecayMode decaymode) {
    /*
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
    */
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
        printf("ERROR!\n");
        exit(1);
        break;
    }

    printf("ERROR!\n");
    exit(1);
    return false;
}

void ReadRootFiles_r_sp(){

    std::vector<string> names;
    const char* dirname = "/home/jwpark/storage/B2Xsnunu_1/GEN_SIG/output/Ntuple";

    load_files(dirname, &names);

    Loader loader;

    for(unsigned int i = 0; i<names.size(); i++){
        loader.initialize();

        TFile *input_file = new TFile( (dirname+std::string("/")+names.at(i)).c_str(),"read");
        printf("%s (%d/%zu)\n",("Read "+names.at(i) + "... ").c_str(), i, names.size());
        loader.GetData(input_file);
        if (loader.event_info_is_valid() == false) { printf("error!\n"); return; }

        loader.PrintInformation(std::string("========== inital =========="));
        loader.DrawTH2F("MbcVSdeltaE_initial", ";Mbc of B_{tag} [GeV];#DeltaE of B_{tag} [GeV]", 100, 5.24, 5.3, 100, -0.2, 0.2, Loader::Btag, 1, Loader::Btag, 2);

        loader.Cut(Loader::Btag,1,Loader::larger_than,5.2);
        loader.PrintInformation(std::string("========== Mbc > 5.2 =========="));

        loader.Cut(Loader::Btag, 2, Loader::larger_than, -0.5);
        loader.Cut(Loader::Btag, 2, Loader::smaller_than, 0.5);
        loader.PrintInformation(std::string("========== abs(deltaE) < 0.5 =========="));
        loader.DrawTH2F("MbcVSdeltaE_after_loose_MbcDeltaE_cut", ";Mbc of B_{tag} [GeV];#DeltaE of B_{tag} [GeV]", 100, 5.24, 5.3, 100, -0.2, 0.2, Loader::Btag, 1, Loader::Btag, 2);

        loader.DrawTH1F("SignalProbability_Btag_after_loose_MbcDeltaE_cut", "SignalProbability of B_{tag};log_{10}(SignalProbability);Num of candidate", 100, -10, 0, Loader::Btag, 5, Loader::Log);
        loader.Cut(Loader::Btag, 5, Loader::larger_than, 0.01);
        loader.PrintInformation(std::string("========== SignalProbability > 0.01 =========="));
        loader.DrawTH2F("MbcVSdeltaE_after_SignalProbability_cut", ";Mbc of B_{tag} [GeV];#DeltaE of B_{tag} [GeV]", 100, 5.24, 5.3, 100, -0.2, 0.2, Loader::Btag, 1, Loader::Btag, 2);

        loader.DrawTH1F("ROE_Eecl_Upsilon_after_SignalProbability_cut", "E_ecl in ROE of #Upsilon(4S);E_{ecl} [GeV];candidates", 100, -0.1, 8, Loader::Upsilon, 3);
        loader.Cut(Loader::Upsilon, 3, Loader::smaller_than, 1.2);
        loader.PrintInformation(std::string("========== E_ecl < 1.2 GeV =========="));

        loader.DrawTH1F("nROE_track_Upsilon_after_E_ROE_cut", "number of raw tracks in ROE of #Upsilon(4S);number of raw tracks;evt", 100, -0.5, 13.5, Loader::Upsilon, 10);
        loader.Cut(Loader::Upsilon, 10, Loader::smaller_than, 0.5);
        loader.PrintInformation(std::string("========== nRawtrack = 0 =========="));

        loader.DrawTH1F("nROE_pi0_after_ntrack_cut", "number of #pi^{0} candidates in ROE of #Upsilon(4S);number of #pi^{0} candidates;evt", 100, -0.5, 13.5, Loader::Upsilon, 6);
        loader.Cut(Loader::Upsilon, 6, Loader::smaller_than, 0.5);
        loader.PrintInformation(std::string("========== npi0 = 0 =========="));

        loader.DrawTH1F("momentum_Bsig_after_missingmomentumtheta_cut", "momentum of B_{sig} at CMS;momentum [GeV];evt", 100, 0, 3.2, Loader::Bsig, 4);
        loader.Cut(Loader::Bsig, 4, Loader::smaller_than, 2.96);
        loader.Cut(Loader::Bsig, 4, Loader::larger_than, 1.6);
    	loader.PrintInformation(std::string("========== 1.6 < momentum of signal side < 2.96 =========="));

        loader.DrawTH1F("missing_momentum_theta_after_BCS", "#theta_{missing};#theta_{missing} [rad];evt", 100, 0, 3.2, Loader::Upsilon, 7);
        loader.Cut(Loader::Upsilon, 7, Loader::smaller_than, 2.618);
        loader.Cut(Loader::Upsilon, 7, Loader::larger_than, 0.297);
	    loader.PrintInformation(std::string("========== 0.297 < missing momentum theta < 2.618 =========="));

        loader.DrawTH1F("SignalProbability_Btag_before_BCS", "SignalProbability of B_{tag};log_{10}(SignalProbability);Num of candidate", 100, -10, 0, Loader::Btag, 5, Loader::Log);
        loader.BCS(Loader::Btag, 5, Loader::Highest);
        if (loader.IsBCSValid() == false) {
            printf("ERROR!\n");
            exit(1);
        }
        loader.DrawTH1F("SignalProbability_Btag_after_BCS", "SignalProbability of B_{tag};log_{10}(SignalProbability);Num of candidate", 100, -10, 0, Loader::Btag, 5, Loader::Log);
        loader.PrintInformation(std::string("========== BCS =========="));
        loader.DrawTH2F("MbcVSdeltaE_after_BCS", ";Mbc of B_{tag} [GeV];#DeltaE of B_{tag} [GeV]", 100, 5.24, 5.3, 100, -0.2, 0.2, Loader::Btag, 1, Loader::Btag, 2);

        loader.Cut(Loader::Btag, 1, Loader::larger_than, 5.27);
        loader.PrintInformation(std::string("========== Mbc > 5.27 =========="));
        loader.DrawTH2F("MbcVSdeltaE_after_Mbc_strict_cut", ";Mbc of B_{tag} [GeV];#DeltaE of B_{tag} [GeV]", 100, 5.24, 5.3, 100, -0.2, 0.2, Loader::Btag, 1, Loader::Btag, 2);

        loader.Cut(Loader::Btag, 2, Loader::larger_than, -0.1);
        loader.Cut(Loader::Btag, 2, Loader::smaller_than, 0.1);
        loader.PrintInformation(std::string("========== abs(deltaE) < 0.1 =========="));
        loader.DrawTH2F("MbcVSdeltaE_after_deltaE_strict_cut", ";Mbc of B_{tag} [GeV];#DeltaE of B_{tag} [GeV]", 100, 5.24, 5.3, 100, -0.2, 0.2, Loader::Btag, 1, Loader::Btag, 2);

        loader.DrawTH1F("nROE_ECLcluster_Upsilon", "number of ECL clusters in ROE of #Upsilon(4S);number of ECL clusters;evt", 14, -0.5, 13.5, Loader::Upsilon, 0);
        loader.DrawTH1F("nROE_KLMcluster_Upsilon", "number of KLM clusters in ROE of #Upsilon(4S);number of KLM clusters;evt", 14, -0.5, 13.5, Loader::Upsilon, 1);
        loader.DrawTH1F("Bsig_p_LAB", "momentum of B_{sig} at LAB frame;p [GeV];evt", 50, -0.5, 6, Loader::Bsig, 3);
        loader.DrawTH1F("Bsig_p_CMS", "momentum of B_{sig} at CMS frame;p [GeV];evt", 50, -0.5, 6, Loader::Bsig, 4);
        loader.DrawTH1F("Bsig_p_RecoilRest", "momentum of B_{sig} at rest frame of recoil system;p [GeV];evt", 50, -0.5, 6, Loader::Bsig, 5);
        loader.DrawTH1F("Btag_dmID", "decay ID of B_{tag};decay ID;evt", 74, -0.5, 36.5, Loader::Btag, 0);
        loader.DrawTH1F("nROE_K_S0", "number of K_S0 candidates in ROE of #Upsilon(4S);number of K_{S}^{0} candidates in ROE;evt", 100, -0.5, 5.5, Loader::Upsilon, 5);
        loader.DrawTH1F("theta_missing_momentum", "#theta of missing momentum;#theta [rad];evt", 50, 0, 3.2, Loader::Upsilon, 7);
        loader.DrawTH1F("M_Xs", "mass of X_{s};M_{Xs} [GeV];evt", 100, 0, 3.5, Loader::Bsig, 6);

        loader.PrintRootFile(std::string("output.root"));
    }
    loader.End();
}
