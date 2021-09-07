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

    double variables[3];
    // 0: M, 1: dM, 2: significanceOfDistance

} Data; 

class Loader {
public:
    enum Variable
    {
        variables = 0
    };
    enum Inequality
    {
        larger_than = 0,
        smaller_than
    };
    enum ValueOption {
        Linear = 0,
        Log
    };

private:
    std::queue<Data> TotalData;
    std::vector<TH1F*> TH1Fs;
    int current_TH1F;

    std::vector<int> N_events;
    std::vector<int> N_candidates;
    std::vector<std::string> titles;
    int current_N_event;
    int current_N_candidate;

public:
    Loader();
    void initialize();
    void GetData(TFile* input_file);
    bool event_info_is_valid();
    void DrawTH1F(const char* name, const char* title, int nbins, double x_low, double x_high, Loader::Variable variable, int i, Loader::ValueOption dr = Loader::Linear);
    void PrintInformation(std::string title);
    void Cut(Loader::Variable variable, int i, Loader::Inequality inq, double value);
    void End();
};

Loader::Loader() {
    std::queue<Data> empty;
    std::swap(TotalData, empty);
    current_TH1F = 0;
    current_N_event = 0;
    current_N_candidate = 0;
}

void Loader::initialize() {
    std::queue<Data> empty;
    std::swap(TotalData, empty);
    current_TH1F = 0;
    current_N_event = 0;
    current_N_candidate = 0;
}

void Loader::GetData(TFile* input_file) {
    TTree* tree_KaonS0 = (TTree*)input_file->Get("KaonS0");

    Data temp = { 0 };

    // get event_info
    tree_KaonS0->SetBranchAddress("__experiment__", &temp.__experiment__);
    tree_KaonS0->SetBranchAddress("__run__", &temp.__run__);
    tree_KaonS0->SetBranchAddress("__event__", &temp.__event__);
    tree_KaonS0->SetBranchAddress("__candidate__", &temp.__candidate__);
    tree_KaonS0->SetBranchAddress("__ncandidates__", &temp.__ncandidates__);

    // get Upsilon_info
    tree_KaonS0->SetBranchAddress("M", &temp.variables[0]);
    tree_KaonS0->SetBranchAddress("dM", &temp.variables[1]);
    tree_KaonS0->SetBranchAddress("significanceOfDistance", &temp.variables[2]);
    // 0: M, 1: dM, 2: significanceOfDistance

    printf("%lld entries...\n", tree_KaonS0->GetEntries());
    for (unsigned int j = 0; j < tree_KaonS0->GetEntries(); j++) { // Fill
        tree_Gamma->GetEntry(j);
        TotalData.push(temp);
    }
    input_file->Close();
}

bool Loader::event_info_is_valid() {
    std::queue<Data> temp_queue;
    while (!TotalData.empty()) {
        Data temp_data = TotalData.front();
        TotalData.pop();
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

        if (variable == Loader::variables) {
            if (dr == Loader::Linear) temp_hist->Fill(temp_data.variables[i]);
            else if (dr == Loader::Log) temp_hist->Fill(TMath::Log10(temp_data.variables[i]));
            else { printf("ERROR!\n"); exit(1); }
        }
        else {
            printf("ERROR!\n");
            exit(1);
        }

    }

    current_TH1F++;
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
        if (variable == Loader::variables) {
            if (inq == Loader::larger_than && temp_data.variables[i] > value) temp_queue.push(temp_data);
            else if (inq == Loader::smaller_than && temp_data.variables[i] < value) temp_queue.push(temp_data);
        }
        else {
            printf("ERROR!\n");
            exit(1);
        }
    }
    TotalData = temp_queue;
}

void Loader::End() {

    for (int i = 0; i < N_events.size();i++) {
        printf("%s\n", titles.at(i).c_str());
        printf("Number of event: %d\n", N_events.at(i));
        printf("Number of candidate: %d\n", N_candidates.at(i));
    }

    for (int i = 0; i < TH1Fs.size(); i++) {
        TCanvas* c_temp = new TCanvas("c", "", 1500, 1200); c_temp->cd();
        TH1Fs.at(i)->Draw("Hist"); c_temp->SaveAs((std::string(TH1Fs.at(i)->GetName()) + ".png").c_str());
        delete c_temp;
    }

}

void ReadKaonS0Files_r_sp(){

    std::vector<string> names;
    const char* dirname = "/home/jwpark/storage/B02Xsnunu_decay_modes/output";

    load_files(dirname, &names);

    Loader loader;

    for(unsigned int i = 0; i<names.size(); i++){
        loader.initialize();

        TFile *input_file = new TFile( (dirname+std::string("/")+names.at(i)).c_str(),"read");
        printf("%s (%d/%zu)\n",("Read "+names.at(i) + "... ").c_str(), i, names.size());
        loader.GetData(input_file);
        if (loader.event_info_is_valid() == false) { printf("error!\n"); return; }

        loader.PrintInformation(std::string("========== inital =========="));

        loader.DrawTH1F("significance_of_distance", "significance of distance of K_{S}^{0} candidates;significance of distance;evt", 100, 0, 100, Loader::variables, 2);
        loader.Cut(Loader::variables, 2, Loader::larger_than, 50);

        loader.DrawTH1F("M", "Mass of K_{S}^{0} candidates;mass [GeV];evt", 100, 0.45, 0.55, Loader::variables, 0);
        loader.DrawTH1F("dM", "dM of K_{S}^{0} candidates;dM [GeV];evt", 100, -0.05, 0.05, Loader::variables, 1);
        loader.Cut(Loader::variables, 1, Loader::larger_than, -0.01);
        loader.Cut(Loader::variables, 1, Loader::smaller_than, 0.01);
    }
    loader.End();
}
