// last update: 2021-10-13
// for Belle2 data

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

void ReadInfoFiles_r_sp(){

    std::vector<string> names;
    const char* dirname = "/home/jwpark/storage/BKG_gbasf2/BKG_total/CHG_after_FEISKIM_light_0_fix/train";

    load_files(dirname, &names);

    Float_t Lower[6];

    Lower[0] = 0;
    Lower[1] = 4;
    Lower[2] = 8;
    Lower[3] = 12;
    Lower[4] = 16;
    Lower[5] = 22.9;
    TH1F* Plot = new TF1F("q2", ";q^{2} [GeV^{2}]; evt", 5, Lower);

    for(unsigned int i = 0; i<names.size(); i++){
        loader.initialize();

        TFile *input_file = new TFile( (dirname+std::string("/")+names.at(i)).c_str(),"read");
        printf("%s (%d/%zu)\n",("Read "+names.at(i) + "... ").c_str(), i, names.size());

        // Get Data
        TTree* tree_info = (TTree*)input_file->Get("info");

        int __experiment__;
        int __run__;
        int __event__;
        int __candidate__;
        int __ncandidates__;
        double s;

        // get event_info
        tree_info->SetBranchAddress("__experiment__", &__experiment__);
        tree_info->SetBranchAddress("__run__", &__run__);
        tree_info->SetBranchAddress("__event__", &__event__);
        tree_info->SetBranchAddress("__candidate__", &__candidate__);
        tree_info->SetBranchAddress("__ncandidates__", &__ncandidates__);
        tree_info->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &s);

        printf("%lld entries...\n", tree_info->GetEntries());
        for (unsigned int j = 0; j < tree_info->GetEntries(); j++) { // Fill
            tree_info->GetEntry(j);
            Plot->Fill(s*s);
        }
        input_file->Close();
    }

    TCanvas* c_temp = new TCanvas("c", "", 1500, 1200); c_temp->cd();
    Plot->Draw("Hist");
}
