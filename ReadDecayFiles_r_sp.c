# define N_decay 6 // five decay mode + others

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

    double Mxs;
    double N_B_ch_total;
    double N_B_ch_1;
    double N_B_ch_2;
    double N_B_ch_3;
    double N_B_ch_4;
    double N_B_ch_5;

} Data; 

class Loader {
private:
    std::queue<Data> TotalData;

    std::vector<int> N_events;
    std::vector<int> N_candidates;
    std::vector<std::string> titles;
    int current_N_event;
    int current_N_candidate;

    std::vector<THStack*> THStacks;
    std::vector<TH1F*> TH1Fs_THStack[N_decay];
    int current_THStack;

public:
    Loader();
    void initialize();
    void GetData(TFile* input_file);
    bool event_info_is_valid();
    void DrawTHStack(const char* name, const char* title, int nbins, double x_low, double x_high);
    void PrintInformation(std::string title);
    void End();
};

Loader::Loader() {
    std::queue<Data> empty;
    std::swap(TotalData, empty);
    current_N_event = 0;
    current_N_candidate = 0;
    current_THStack = 0;
}

void Loader::initialize() {
    std::queue<Data> empty;
    std::swap(TotalData, empty);
    current_N_event = 0;
    current_N_candidate = 0;
    current_THStack = 0;
}

void Loader::GetData(TFile* input_file) {
    TTree* tree_Xsd = (TTree*)input_file->Get("Xsd");

    Data temp = { 0 };

    // get event_info
    tree_Xsd->SetBranchAddress("__experiment__", &temp.__experiment__);
    tree_Xsd->SetBranchAddress("__run__", &temp.__run__);
    tree_Xsd->SetBranchAddress("__event__", &temp.__event__);
    tree_Xsd->SetBranchAddress("__candidate__", &temp.__candidate__);
    tree_Xsd->SetBranchAddress("__ncandidates__", &temp.__ncandidates__);

    // get Upsilon_info
    tree_Xsd->SetBranchAddress("daughter__bo0__cmM__bc", &temp.Mxs);
    tree_Xsd->SetBranchAddress("nParticlesInList__boB0__clch_total__bc", &temp.N_B_ch_total);
    tree_Xsd->SetBranchAddress("nParticlesInList__boB0__clch1_total__bc", &temp.N_B_ch_1);
    tree_Xsd->SetBranchAddress("nParticlesInList__boB0__clch2_total__bc", &temp.N_B_ch_2);
    tree_Xsd->SetBranchAddress("nParticlesInList__boB0__clch3_total__bc", &temp.N_B_ch_3);
    tree_Xsd->SetBranchAddress("nParticlesInList__boB0__clch4_total__bc", &temp.N_B_ch_4);
    tree_Xsd->SetBranchAddress("nParticlesInList__boB0__clch5_total__bc", &temp.N_B_ch_5);

    printf("%lld entries...\n", tree_Xsd->GetEntries());
    for (unsigned int j = 0; j < tree_Xsd->GetEntries(); j++) { // Fill
        tree_Xsd->GetEntry(j);
        TotalData.push(temp);
    }
    input_file->Close();
}

bool Loader::event_info_is_valid() {
    std::queue<Data> temp_queue;
    while (!TotalData.empty()) {
        Data temp_data = TotalData.front();
        TotalData.pop();

        if (temp_data.N_B_ch_total > 1.5) { // error
	printf("N_B_ch_total > 1: this event is rejected!\n");
	continue;
	//            return false;
        }

        temp_queue.push(temp_data);
    }
    TotalData = temp_queue;
    return true;
}

void Loader::DrawTHStack(const char* name, const char* title, int nbins, double x_low, double x_high) {
    if (THStacks.size() == current_THStack) { // allocate new thstacks
        THStack* stack = new THStack(name, title);
        THStacks.push_back(stack);
        for (int i = 0; i < N_decay; i++) {
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

    TH1F* temp_hist[N_decay];
    for (int i = 0; i < N_decay;i++) {
        temp_hist[i] = TH1Fs_THStack[i].at(current_THStack);
    }

    std::queue<Data> temp_queue = TotalData;
    while (!temp_queue.empty()) {
        Data temp_data = temp_queue.front();
        temp_queue.pop();

        if (temp_data.N_B_ch_total > 1.5) { // error
            printf("ERROR!\n");
            exit(1);
        }

        int decaymodeid = -1;

        if(temp_data.N_B_ch_1 > 0.5) decaymodeid = 1;
        else if (temp_data.N_B_ch_2 > 0.5) decaymodeid = 2;
        else if (temp_data.N_B_ch_3 > 0.5) decaymodeid = 3;
        else if (temp_data.N_B_ch_4 > 0.5) decaymodeid = 4;
        else if (temp_data.N_B_ch_5 > 0.5) decaymodeid = 5;
        else decaymodeid = 0;

        temp_hist[decaymodeid]->Fill(temp_data.Mxs);
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

void Loader::End() {

    for (int i = 0; i < N_events.size();i++) {
        printf("%s\n", titles.at(i).c_str());
        printf("Number of event: %d\n", N_events.at(i));
        printf("Number of candidate: %d\n", N_candidates.at(i));
    }

    for (int i = 0; i < THStacks.size(); i++) {
        TCanvas* c_temp = new TCanvas("c", "", 1500, 1200); c_temp->cd();
        gStyle->SetPalette(kOcean);

        for (int j = 0; j < N_decay; j++) {
		THStacks.at(i)->Add(TH1Fs_THStack[j].at(i));
                TH1Fs_THStack[j].at(i)->Draw("Hist");
		c_temp->SaveAs( (std::string(TH1Fs_THStack[j].at(i)->GetName()) + ".png").c_str() );
	}
        THStacks.at(i)->Draw("pfc"); 
	//gPad->BuildLegend();
	c_temp->SaveAs((std::string(THStacks.at(i)->GetName()) + ".png").c_str());
        delete c_temp;
    }
}

void ReadDecayFiles_r_sp(){

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
        loader.DrawTHStack("Mxs", ";M_{Xs} [GeV]; evt", 100, 1.05, 4);

    }
    loader.End();
}
