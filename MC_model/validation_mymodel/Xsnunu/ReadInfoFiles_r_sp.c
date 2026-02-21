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

void LetsFill(const char* dirname, TH1D* hist_q2, TH1D* hist_MXs) {

    double invM = -1;
    double MXs = -1;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_info = (TTree*)input_file->Get("info");

        tree_info->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM); 
        tree_info->SetBranchAddress("daughter__bo0__cmM__bc", &MXs);

        printf("%lld entries...\n", tree_info->GetEntries());
        for (unsigned int j = 0; j < tree_info->GetEntries(); j++) { // Fill
            tree_info->GetEntry(j);

            hist_q2->Fill(invM* invM);
            hist_MXs->Fill(MXs);

        }
        input_file->Close();

    }

}

void ReadInfoFiles_r_sp(){

    const char* Previous_dirname = "./previous";
    const char* CC_dirname = "./CC";

    TH1D* previous_q2 = new TH1D("previous_q2", ";;", 100, 0, 17.6);
    TH1D* previous_MXs = new TH1D("previous_MXs", ";;", 100, 1.1, 3.6);

    TH1D* CC_q2 = new TH1D("CC_q2", ";;", 100, 0, 17.6);
    TH1D* CC_MXs = new TH1D("CC_MXs", ";;", 100, 1.1, 3.6);

    LetsFill(Previous_dirname, previous_q2, previous_MXs);
    LetsFill(CC_dirname, CC_q2, CC_MXs);

    TCanvas* c_temp = new TCanvas("c", "", 1200, 1200); c_temp->cd();

    previous_q2->Draw("HIST");
    c_temp->SaveAs("previous_q2.png");

    previous_MXs->Draw("HIST");
    c_temp->SaveAs("previous_MXs.png");

    CC_q2->Draw("HIST");
    c_temp->SaveAs("CC_q2.png");

    CC_MXs->Draw("HIST");
    c_temp->SaveAs("CC_MXs.png");
}
