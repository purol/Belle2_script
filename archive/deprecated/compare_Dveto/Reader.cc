void load_files(const char* dirname, std::vector<std::string>* names) {
    TSystemDirectory dir(dirname, dirname);
    TList* files = dir.GetListOfFiles();
    if (files) {
        TSystemFile* file;
        TString fname;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(".root")) {
                names->push_back(fname.Data());
            }
        }
    }
}

int Reader() {

    TH1D* hist = new TH1D("Dc", ";Mass of D candidates [GeV];arbitrary unit", 100, 0.0, 2.1);

    std::vector<std::string> names;
    const char* dirname = "./";

    load_files(dirname, &names);

    double var = -1;

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Dc_M_3", &var);

        printf("%lld entries...\n", tree_Bsig->GetEntries());
        for (unsigned int j = 0; j < tree_Bsig->GetEntries(); j++) { // Fill
            tree_Bsig->GetEntry(j);

            hist->Fill(var);

        }

    }

    TCanvas* c_temp = new TCanvas("c", "", 800, 800); c_temp->cd();
    
    hist->SetStats(false);
    hist->Scale(1.0 / hist->Integral(), "width");
    hist->Draw("Hist");

    c_temp->SaveAs("mass.png");

	return 0;
}
