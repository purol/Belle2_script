# define N_Needed_info 25

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

    int event_info[15];
    // 0: upsilon_experiment, 1: upsilon_run, 2: upsilon_event, 3: upsilon_candidate, 4: upsilon_ncandidates
    // 5: Bsig_experiment, 6: Bsig_run, 7: Bsig_event, 8: Bsig_candidate, 9: Bsig_ncandidates
    // 10: Btag_experiment, 11: Btag_run, 12: Btag_event, 13: Btag_candidate, 14: Btag_ncandidates

    double Upsilon_info[13];
    // 0: Upsilon_isSignal; 1: number of ECL clusters in ROE(cleanMask), 2: number of KLM clusters in ROE
    // 3: energy in ROE(cleanMask), 4: number of tracks in ROE(cleanMask), 5: roeEextra(cleanMask)
    // 6: nROE_NeutralECLClusters(cleanMask), 7: roeNeextra(cleanMask), 8: energy in ROE(cleanMask) at CMS
    // 9: roeExtra(cleanMask) at CMS, 10: roeNeextra(cleanMask) at CMS, 11: nROE_K_S0
    // 12: nROE_pi0, 

    double Bsig_info[14];
    // 0: Bsig_isSignal, 1: Bsig_E, 2: Bsig_E_CMS, 3: Bsig_E_Recoil, 4: Bsig_dmID
    // 5: Bsig_first_daughter's_actPID(3,2) 6: Bsig_first_daughter's_mcPDG
    // 7: Bsig_p, 8: Bsig_p_CMS, 9: Bsig_p_Recoil
    // 10: Kaon dr, 11: Kaon dz, 12: Bsig_first_daughter's eIDBelle
    // 13: Bsig_first_daughter's_muIDBelle

    double Btag_info[7];
    // 0: Btag_isSignal, 1:Btag_dmID, 2: Btag_Mbc, 3: Btag_deltaE
    // 4: Btag_E, 5: Btag_E_CMS, 6: Btag_signalprobability

    double Needed_info[N_Needed_info];
    // 0: R2, 1: thrustBm, 2: thrustOm. 3: cosTBTO
    // 4: cosTBz, 5: KSFWVariables(et), 6: KSFWVariables(mm2)
    // 7: KSFWVariables(hso00), 8: KSFWVariables(hso02), 9: KSFWVariables(hso04)
    // 10: KSFWVariables(hso10), 11: KSFWVariables(hso12), 12: KSFWVariables(hso14)
    // 13: KSFWVariables(hso20), 14: KSFWVariables(hso22), 15: KSFWVariables(hso24)
    // 16: CleoConeCS(1), 17: CleoConeCS(2), 18: CleoConeCS(3)
    // 19: CleoConeCS(4), 20: CleoConeCS(5), 21: CleoConeCS(6)
    // 22: CleoConeCS(7), 23: CleoConeCS(8), 24:CleoConeCS(9)

} Data; 

class Loader {
private:
    std::queue<Data> TotalData;
    std::vector<TH1F*> TH1Fs;
    std::vector<TH2F*> TH2Fs;
    int current_TH1F;
    int current_TH2F;

    std::vector<int> N_events;
    std::vector<int> N_candidates;
    std::vector<std::string> titles;
    int current_N_event;
    int current_N_candidate;

    std::vector<TFile*> files;
    std::vector<TTree*> trees;
    int current_file;

    double DataToTree[N_Needed_info];

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

    Loader();
    void initialize();
    void GetData(TFile* input_file);
    bool event_info_is_valid();
    void DrawTH1F(TH1F* hist, int i, Loader::Variable variable, Loader::ValueOption dr = Loader::Linear);
    void DrawTH2F(TH2F* hist, int i, int j, Loader::Variable variable_1, Loader::Variable variable_2);
    void PrintInformation(std::string title);
    void Cut(Loader::Variable variable, int i, Loader::Inequality inq, double value);
    void BCS(Loader::Variable variable, int index, Loader::BCS_criterion crit);
    bool IsBCSValid();
    void End();
    void PrintRootFile(std::string output_name);
};

Loader::Loader() {
    std::queue<Data> empty;
    std::swap(TotalData, empty);
    current_TH1F = 0;
    current_TH2F = 0;
    current_N_event = 0;
    current_N_candidate = 0;
    current_file = 0;
}

void Loader::initialize() {
    std::queue<Data> empty;
    std::swap(TotalData, empty);
    current_TH1F = 0;
    current_TH2F = 0;
    current_N_event = 0;
    current_N_candidate = 0;
    current_file = 0;
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

    // get Upsilon_info
    tree_upsilon->SetBranchAddress("isSignal", &temp.Upsilon_info[0]);
    tree_upsilon->SetBranchAddress("nROE_ECLClusters__bocleanMask__bc", &temp.Upsilon_info[1]);
    tree_upsilon->SetBranchAddress("nROE_KLMClusters", &temp.Upsilon_info[2]);
    tree_upsilon->SetBranchAddress("roeE__bocleanMask__bc", &temp.Upsilon_info[3]);
    tree_upsilon->SetBranchAddress("nROE_Tracks__bocleanMask__bc", &temp.Upsilon_info[4]);
    tree_upsilon->SetBranchAddress("roeEextra__bocleanMask__bc", &temp.Upsilon_info[5]);
    tree_upsilon->SetBranchAddress("nROE_NeutralECLClusters__bocleanMask__bc", &temp.Upsilon_info[6]);
    tree_upsilon->SetBranchAddress("roeNeextra__bocleanMask__bc", &temp.Upsilon_info[7]);
    tree_upsilon->SetBranchAddress("useCMSFrame__boroeE__bocleanMask__bc__bc", &temp.Upsilon_info[8]);
    tree_upsilon->SetBranchAddress("useCMSFrame__boroeEextra__bocleanMask__bc__bc", &temp.Upsilon_info[9]);
    tree_upsilon->SetBranchAddress("useCMSFrame__boroeNeextra__bocleanMask__bc__bc", &temp.Upsilon_info[10]);
    tree_upsilon->SetBranchAddress("nROE_ParticlesInList__boK_S0__clgood__bc", &temp.Upsilon_info[11]);
    tree_upsilon->SetBranchAddress("nROE_ParticlesInList__bopi0__clgood__bc", &temp.Upsilon_info[12]);

    // get Bsig_info
    tree_Bsig->SetBranchAddress("Bsig_isSignal", &temp.Bsig_info[0]);
    tree_Bsig->SetBranchAddress("Bsig_E", &temp.Bsig_info[1]);
    tree_Bsig->SetBranchAddress("Bsig_useCMSFrame_E", &temp.Bsig_info[2]);
    tree_upsilon->SetBranchAddress("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc", &temp.Bsig_info[3]);
    tree_Bsig->SetBranchAddress("Bsig_extraInfo_decayModeID", &temp.Bsig_info[4]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_atcPIDBelle_3_2", &temp.Bsig_info[5]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_mcPDG", &temp.Bsig_info[6]);
    tree_Bsig->SetBranchAddress("Bsig_p", &temp.Bsig_info[7]);
    tree_Bsig->SetBranchAddress("Bsig_useCMSFrame_p", &temp.Bsig_info[8]);
    tree_upsilon->SetBranchAddress("useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc", &temp.Bsig_info[9]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_dr", &temp.Bsig_info[10]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_dz", &temp.Bsig_info[11]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_eIDBelle", &temp.Bsig_info[12]);
    tree_Bsig->SetBranchAddress("Bsig_daughter_0_muIDBelle", &temp.Bsig_info[13]);

    // get Btag_info
    tree_Btag->SetBranchAddress("Btag_isSignal", &temp.Btag_info[0]);
    tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &temp.Btag_info[1]);
    tree_Btag->SetBranchAddress("Btag_Mbc", &temp.Btag_info[2]);
    tree_Btag->SetBranchAddress("Btag_deltaE", &temp.Btag_info[3]);
    tree_Btag->SetBranchAddress("Btag_E", &temp.Btag_info[4]);
    tree_Btag->SetBranchAddress("Btag_useCMSFrame_E", &temp.Btag_info[5]);
    tree_Btag->SetBranchAddress("Btag_extraInfo_SignalProbability", &temp.Btag_info[6]);

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

void Loader::DrawTH1F(TH1F* hist, int i, Loader::Variable variable, Loader::ValueOption dr = Loader::Linear) {
    if (TH1Fs.size() == current_TH1F) { // allocate new hist
        TH1Fs.push_back(hist);
    }
    else if (TH1Fs.size() > current_TH1F) { // use what I have
        delete hist;
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

void Loader::DrawTH2F(TH2F* hist, int i, int j, Loader::Variable variable_1, Loader::Variable variable_2) {
    if (TH2Fs.size() == current_TH2F) { // allocate new hist
        TH2Fs.push_back(hist);
    }
    else if (TH2Fs.size() > current_TH2F) { // use what I have
        delete hist;
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
        TTree* temp_tree = trees.at(i);
        temp_file->cd();
        temp_tree->Write();
        temp_file->Close();
    }
}

void Loader::PrintRootFile(std::string output_name) {
    if (files.size() == current_file && trees.size() == current_file) { // allocate new TFile and TTree
        TFile* file = new TFile(output_name.c_str(), "recreate");
        file->cd();
        TTree* tree = new TTree("data", "data");
        tree->Branch("R2", &DataToTree[0]);
        tree->Branch("thrustBm", &DataToTree[1]);
        tree->Branch("thrustOm", &DataToTree[2]);
        tree->Branch("cosTBTO", &DataToTree[3]);
        tree->Branch("cosTBz", &DataToTree[4]);
        tree->Branch("KSFWVariables_et", &DataToTree[5]);
        tree->Branch("KSFWVariables_mm2", &DataToTree[6]);
        tree->Branch("KSFWVariables_hso00", &DataToTree[7]);
        tree->Branch("KSFWVariables_hso02", &DataToTree[8]);
        tree->Branch("KSFWVariables_hso04", &DataToTree[9]);
        tree->Branch("KSFWVariables_hso10", &DataToTree[10]);
        tree->Branch("KSFWVariables_hso12", &DataToTree[11]);
        tree->Branch("KSFWVariables_hso14", &DataToTree[12]);
        tree->Branch("KSFWVariables_hso20", &DataToTree[13]);
        tree->Branch("KSFWVariables_hso22", &DataToTree[14]);
        tree->Branch("KSFWVariables_hso24", &DataToTree[15]);
        tree->Branch("CleoConeCS_1", &DataToTree[16]);
        tree->Branch("CleoConeCS_2", &DataToTree[17]);
        tree->Branch("CleoConeCS_3", &DataToTree[18]);
        tree->Branch("CleoConeCS_4", &DataToTree[19]);
        tree->Branch("CleoConeCS_5", &DataToTree[20]);
        tree->Branch("CleoConeCS_6", &DataToTree[21]);
        tree->Branch("CleoConeCS_7", &DataToTree[22]);
        tree->Branch("CleoConeCS_8", &DataToTree[23]);
        tree->Branch("CleoConeCS_9", &DataToTree[24]);
        files.push_back(file);
        trees.push_back(tree);
    }
    else if (files.size() > current_file && trees.size() > current_file && files.size() == trees.size()) { // use what I have
    }
    else { // error
        printf("ERROR!\n");
        exit(1);
    }

    TFile* temp_file = files.at(current_file);
    temp_file->cd();
    TTree* temp_tree = trees.at(current_file);
    std::queue<Data> temp_queue = TotalData;
    while (!temp_queue.empty()) {
        Data temp = temp_queue.front();
        temp_queue.pop();

        for (int i = 0; i < N_Needed_info; i++) {
            DataToTree[i] = temp.Needed_info[i];
        }
        temp_tree->Fill();
    }

    current_file++;
}

void ReadRootFiles_r_sp(){

    std::vector<string> names;
    const char* dirname = "/home/jwpark/storage/Ntuple32";

    load_files(dirname, &names);

    Loader loader;

    for(unsigned int i = 0; i<names.size(); i++){
        loader.initialize();

        TFile *input_file = new TFile( (dirname+std::string("/")+names.at(i)).c_str(),"read");
        printf("%s (%d/%zu)\n",("Read "+names.at(i) + "... ").c_str(), i, names.size());
        loader.GetData(input_file);
        if (loader.event_info_is_valid() == false) { printf("error!\n"); return; }

        loader.PrintInformation(std::string("========== inital =========="));
        loader.DrawTH2F(new TH2F("MbcVSdeltaE_initial", ";Mbc of B_{tag} [GeV];#DeltaE of B_{tag} [GeV]", 100, 5.2, 5.3, 100, -0.5, 0.5), 2, 3, Loader::Btag, Loader::Btag);

        loader.Cut(Loader::Btag,2,Loader::larger_than,5.2);
        loader.PrintInformation(std::string("========== Mbc > 5.2 =========="));

        loader.Cut(Loader::Btag, 3, Loader::larger_than, -0.5);
        loader.Cut(Loader::Btag, 3, Loader::smaller_than, 0.5);
        loader.PrintInformation(std::string("========== abs(deltaE) < 0.5 =========="));
        loader.DrawTH2F(new TH2F("MbcVSdeltaE_after_loose_MbcDeltaE_cut", ";Mbc of B_{tag} [GeV];#DeltaE of B_{tag} [GeV]", 100, 5.2, 5.3, 100, -0.5, 0.5), 2, 3, Loader::Btag, Loader::Btag);

        loader.DrawTH1F(new TH1F("SignalProbability_Btag_after_loose_MbcDeltaE_cut", "SignalProbability of B_{tag};log_{10}(SignalProbability);Num of candidate", 100, -10, 0), 6, Loader::Btag, Loader::Log);
        loader.Cut(Loader::Btag, 6, Loader::larger_than, 0.01);
        loader.PrintInformation(std::string("========== SignalProbability > 0.01 =========="));
        loader.DrawTH2F(new TH2F("MbcVSdeltaE_after_SignalProbability_cut", ";Mbc of B_{tag} [GeV];#DeltaE of B_{tag} [GeV]", 100, 5.2, 5.3, 100, -0.5, 0.5), 2, 3, Loader::Btag, Loader::Btag);

        loader.DrawTH1F(new TH1F("atcPID(3,2)", "atcPID(3,2) of daughter of B_{sig};atcPID(3,2);Num of candidate", 100, -0.1, 1.1), 5, Loader::Bsig);
        loader.Cut(Loader::Bsig, 5, Loader::larger_than, 0.6);
        loader.PrintInformation(std::string("========== atcPID(3,2) > 0.6 =========="));
        loader.DrawTH2F(new TH2F("MbcVSdeltaE_after_atcPID_cut", ";Mbc of B_{tag} [GeV];#DeltaE of B_{tag} [GeV]", 100, 5.2, 5.3, 100, -0.5, 0.5), 2, 3, Loader::Btag, Loader::Btag);

        loader.DrawTH1F(new TH1F("dr_Kaon_from_Bsig", "dr of Kaon from B_{sig};dr [cm];candidates", 100, -0.1, 3.5), 10, Loader::Bsig);
        loader.Cut(Loader::Bsig, 10, Loader::smaller_than, 2);
        loader.PrintInformation(std::string("========== dr < 2 cm =========="));

        loader.DrawTH1F(new TH1F("dz_Kaon_from_Bsig", "dz of Kaon from B_{sig};dz [cm];candidates", 100, -6, 6), 11, Loader::Bsig);
        loader.Cut(Loader::Bsig, 11, Loader::smaller_than, 4);
        loader.Cut(Loader::Bsig, 11, Loader::larger_than, -4);
        loader.PrintInformation(std::string("========== abs(dz) < 4 cm =========="));

        loader.DrawTH1F(new TH1F("eIDBelle_Bsig_first_daughter", "eIDBelle of Kaon from B_{sig};eIDBelle;candidates", 100, 0, 1), 12, Loader::Bsig);
        loader.Cut(Loader::Bsig, 12, Loader::smaller_than, 0.9);
        loader.PrintInformation(std::string("========== eIDBelle < 0.9 =========="));

        loader.DrawTH1F(new TH1F("muIDBelle_Bsig_first_daughter", "muIDBelle of Kaon from B_{sig};muIDBelle;candidates", 100, 0, 1), 13, Loader::Bsig);
        loader.Cut(Loader::Bsig, 13, Loader::smaller_than, 0.9);
        loader.PrintInformation(std::string("========== muIDBelle < 0.9 =========="));

        loader.DrawTH1F(new TH1F("ROE_E_Upsilon_after_muIDBelle_cut", "Energy in ROE of #Upsilon(4S) at LAB;energy [GeV];candidates", 100, -0.1, 8), 3, Loader::Upsilon);
        loader.DrawTH1F(new TH1F("ROE_E_Upsilon_CMS_after_muIDBelle_cut", "Energy in ROE of #Upsilon(4S) at CMS;energy [GeV];candidates", 100, -0.1, 8), 8, Loader::Upsilon);
        loader.DrawTH1F(new TH1F("ROE_ECLC_Upsilon_after_muIDBelle_cut", "Energy in ECLClusters in ROE of #Upsilon(4S) at CMS;energy [GeV];candidates", 100, -0.1, 8), 9, Loader::Upsilon);
        loader.DrawTH1F(new TH1F("ROE_NECLC_Upsilon_after_muIDBelle_cut", "Energy in neutral ECLClusters in ROE of #Upsilon(4S) at CMS;energy [GeV];candidates", 100, -0.1, 8), 10, Loader::Upsilon);
        loader.Cut(Loader::Upsilon, 3, Loader::smaller_than, 1);
        loader.PrintInformation(std::string("========== E_ROE < 1 GeV =========="));

        loader.DrawTH1F(new TH1F("nROE_track_Upsilon_after_E_ROE_cut", "number of tracks in ROE of #Upsilon(4S);number of tracks;evt", 100, -0.5, 13.5), 4, Loader::Upsilon);
        loader.Cut(Loader::Upsilon, 4, Loader::smaller_than, 0.5);
        loader.PrintInformation(std::string("========== ntrack = 0 =========="));

        loader.DrawTH1F(new TH1F("SignalProbability_Btag_before_BCS", "SignalProbability of B_{tag};log_{10}(SignalProbability);Num of candidate", 100, -10, 0), 6, Loader::Btag, Loader::Log);
        loader.BCS(Loader::Btag, 6, Loader::Highest);
        if (loader.IsBCSValid() == false) {
            printf("ERROR!\n");
            exit(1);
        }
        loader.DrawTH1F(new TH1F("SignalProbability_Btag_after_BCS", "SignalProbability of B_{tag};log_{10}(SignalProbability);Num of candidate", 100, -10, 0), 6, Loader::Btag, Loader::Log);
        loader.PrintInformation(std::string("========== BCS =========="));
        loader.DrawTH2F(new TH2F("MbcVSdeltaE_after_BCS", ";Mbc of B_{tag} [GeV];#DeltaE of B_{tag} [GeV]", 100, 5.2, 5.3, 100, -0.5, 0.5), 2, 3, Loader::Btag, Loader::Btag);

        loader.Cut(Loader::Btag, 2, Loader::larger_than, 5.27);
        loader.PrintInformation(std::string("========== Mbc > 5.27 =========="));
        loader.DrawTH2F(new TH2F("MbcVSdeltaE_after_Mbc_strict_cut", ";Mbc of B_{tag} [GeV];#DeltaE of B_{tag} [GeV]", 100, 5.2, 5.3, 100, -0.5, 0.5), 2, 3, Loader::Btag, Loader::Btag);

        loader.Cut(Loader::Btag, 3, Loader::larger_than, -0.1);
        loader.Cut(Loader::Btag, 3, Loader::smaller_than, 0.1);
        loader.PrintInformation(std::string("========== abs(deltaE) < 0.1 =========="));
        loader.DrawTH2F(new TH2F("MbcVSdeltaE_after_deltaE_strict_cut", ";Mbc of B_{tag} [GeV];#DeltaE of B_{tag} [GeV]", 100, 5.2, 5.3, 100, -0.5, 0.5), 2, 3, Loader::Btag, Loader::Btag);

        loader.DrawTH1F(new TH1F("nROE_ECLcluster_Upsilon", "number of ECL clusters in ROE of #Upsilon(4S);number of ECL clusters;evt", 14, -0.5, 13.5), 1, Loader::Upsilon);
        loader.DrawTH1F(new TH1F("nROE_KLMcluster_Upsilon", "number of KLM clusters in ROE of #Upsilon(4S);number of KLM clusters;evt", 14, -0.5, 13.5), 2, Loader::Upsilon);
        loader.DrawTH1F(new TH1F("nROE_energy_Upsilon", "energy of ROE of #Upsilon(4S);energy of ROE [GeV];evt", 50, 0, 3), 3, Loader::Upsilon);
        loader.DrawTH1F(new TH1F("nROE_track_Upsilon", "number of tracks in ROE of #Upsilon(4S);number of tracks;evt", 14, -0.5, 13.5), 4, Loader::Upsilon);
        loader.DrawTH1F(new TH1F("ROE_Eextra", "ROE Eextra of #Upsilon(4S);number of tracks;evt", 100, -0.5, 5.5), 5, Loader::Upsilon);
        loader.DrawTH1F(new TH1F("ROE_NECLC_Upsilon", "Energy in neutral ECLClusters in ROE of #Upsilon(4S) at CMS;energy [GeV];candidates", 100, -0.1, 8), 10, Loader::Upsilon);
        loader.DrawTH1F(new TH1F("Bsig_p_LAB", "momentum of B_{sig} at LAB frame;p [GeV];evt", 50, -0.5, 6), 7, Loader::Bsig);
        loader.DrawTH1F(new TH1F("Bsig_p_CMS", "momentum of B_{sig} at CMS frame;p [GeV];evt", 50, -0.5, 6), 8, Loader::Bsig);
        loader.DrawTH1F(new TH1F("Bsig_p_RecoilRest", "momentum of B_{sig} at rest frame of recoil system;p [GeV];evt", 50, -0.5, 6), 9, Loader::Bsig);
        loader.DrawTH1F(new TH1F("Btag_dmID", "decay ID of B_{tag};decay ID;evt", 74, -0.5, 36.5), 1, Loader::Btag);
        loader.DrawTH1F(new TH1F("nROE_K_S0", "number of K_S0:good candidates in ROE of #Upsilon(4S);number of good K_{S}^{0} candidates in ROE;evt", 100, -0.5, 5.5), 11, Loader::Upsilon);
        loader.DrawTH1F(new TH1F("nROE_pi0", "number of #pi^{0} candidates in ROE of #Upsilon(4S);number of #pi^{0} candidates in ROE;evt", 100, -0.5, 8.5), 12, Loader::Upsilon);

        loader.PrintRootFile(std::string("output.root"));
    }
    loader.End();
}
