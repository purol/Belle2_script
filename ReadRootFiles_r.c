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
    int __event__;
    int __candidate__;
    int __ncandidates__;

    double Upsilon_isSignal;

    double Bsig_isSignal;
    double Bsig_E;
    double Bsig_E_CMS;
    double Bsig_E_Recoil;
    double Bsig_dmID;

    double Btag_isSignal;
    double Btag_dmID;
    double Btag_Mbc;
    double Btag_deltaE;
    double Btag_E;
    double Btag_E_CMS;

} Data; 

void ReadRootFiles_r(){

    std::vector<string> names;
    const char* dirname = "/home/jwpark/Ntuple18";

    load_files(dirname, &names);

    std::queue<Data> TotalData;
    Data temp = {0};

    int Bsig_experiment;
    int Bsig_event;
    int Bsig_candidate;
    int Bsig_ncandidates;
    int Btag_experiment;
    int Btag_event;
    int Btag_candidate;
    int Btag_ncandidates;

    unsigned int totalnum_entry = 0;
    for(unsigned int i = 0; i<names.size(); i++){
        TFile *input_file = new TFile( (dirname+std::string("/")+names.at(i)).c_str(),"read");
        TTree *tree_upsilon = (TTree *) input_file -> Get("Upsilon");
        TTree *tree_Bsig = (TTree *) input_file -> Get("Bsig");
        TTree *tree_Btag = (TTree *) input_file -> Get("Btag");

        tree_upsilon->SetBranchAddress("__experiment__", &temp.__experiment__);
        tree_upsilon->SetBranchAddress("__event__", &temp.__event__);
        tree_upsilon->SetBranchAddress("__candidate__", &temp.__candidate__);
        tree_upsilon->SetBranchAddress("__ncandidates__", &temp.__ncandidates__);
        tree_upsilon->SetBranchAddress("isSignal", &temp.Upsilon_isSignal);
        tree_Bsig->SetBranchAddress("Bsig_isSignal", &temp.Bsig_isSignal);
        tree_Bsig->SetBranchAddress("Bsig_E", &temp.Bsig_E);
        tree_Bsig->SetBranchAddress("Bsig_useCMSFrame_E", &temp.Bsig_E_CMS);
        tree_upsilon->SetBranchAddress("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc", &temp.Bsig_E_Recoil);
        tree_Bsig->SetBranchAddress("Bsig_extraInfo_decayModeID", &temp.Bsig_dmID);
        tree_Btag->SetBranchAddress("Btag_isSignal", &temp.Btag_isSignal);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &temp.Btag_dmID);
        tree_Btag->SetBranchAddress("Btag_Mbc", &temp.Btag_Mbc);
        tree_Btag->SetBranchAddress("Btag_deltaE", &temp.Btag_deltaE);
        tree_Btag->SetBranchAddress("Btag_E", &temp.Btag_E);
        tree_Btag->SetBranchAddress("Btag_useCMSFrame_E", &temp.Btag_E_CMS);

        tree_Bsig->SetBranchAddress("__experiment__", &Bsig_experiment);
        tree_Bsig->SetBranchAddress("__event__", &Bsig_event);
        tree_Bsig->SetBranchAddress("__candidate__", &Bsig_candidate);
        tree_Bsig->SetBranchAddress("__ncandidates__", &Bsig_ncandidates);
        tree_Btag->SetBranchAddress("__experiment__", &Btag_experiment);
        tree_Btag->SetBranchAddress("__event__", &Btag_event);
        tree_Btag->SetBranchAddress("__candidate__", &Btag_candidate);
        tree_Btag->SetBranchAddress("__ncandidates__", &Btag_ncandidates);

        printf("%s",("Read "+names.at(i) + "... ").c_str());
        printf("%lld entries...\n",tree_upsilon->GetEntries());
        totalnum_entry += tree_upsilon->GetEntries();
        for(unsigned int j = 0; j < tree_upsilon->GetEntries(); j++){ // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            if(temp.__experiment__ != Bsig_experiment || temp.__experiment__ != Btag_experiment || Btag_experiment != Bsig_experiment) {printf("error\n"); return;
            if(temp.__event__ != Bsig_event || temp.__event__ != Btag_event || Btag_event != Bsig_event) {printf("error\n"); return;}
            if(temp.__candidate__ != Bsig_candidate || temp.__candidate__ != Btag_candidate || Btag_candidate != Bsig_candidate) {printf("error\n"); return;}
            if(temp.__ncandidates__ != Bsig_ncandidates || temp.__ncandidates__ != Btag_ncandidates || Btag_ncandidates != Bsig_ncandidates) {printf("error\n"); return;}}

            TotalData.push(temp);
        }
    }
    printf("Total %d entries\n", totalnum_entry);

    while(!TotalData.empty()){
        Data temp_temp = TotalData.front();
        TotalData.pop();
        if(temp_temp.Btag_Mbc < 5.2) printf("dmID: %f\n", temp_temp.Btag_dmID);
    }
}
