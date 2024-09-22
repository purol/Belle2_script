#define LT_number 51

void load_files(const char* dirname, std::vector<string>* names) {
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

void load_files(const char* dirname, std::vector<string>* names, const char* included_string) {
    TSystemDirectory dir(dirname, dirname);
    TList* files = dir.GetListOfFiles();
    if (files) {
        TSystemFile* file;
        TString fname;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(".root") && fname.Contains(included_string)) {
                names->push_back(fname.Data());
            }
        }
    }
}

void ReadLinearityRootFile_bayesian() {

    double step = 1.0;

    // define LT results
    std::vector<double> in_mus;
    std::vector<double> out_mus;

    // fill LT result
    std::vector<std::string> names;
    const char* dirname = "./";
    load_files(dirname, &names, "LT_result_");
    for (unsigned int k = 0; k < names.size(); k++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(k)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), k, names.size());
        TTree* temp_tree = (TTree*)input_file->Get("LT_result");

        double in_mu;
        double out_mu;
        unsigned int Nentry = temp_tree->GetEntries();

        temp_tree->SetBranchAddress("mu_MXs3_true", &in_mu);
        temp_tree->SetBranchAddress("mu_MXs3_value", &out_mu);

        for (unsigned int j = 0; j < Nentry; j++) { // Fill
            temp_tree->GetEntry(j);
            in_mus.push_back(in_mu);
            out_mus.push_back(out_mu);
        }
        input_file->Close();
    }

    // draw Linearity map
    std::vector<double>::iterator max_it = std::max_element(in_mus.begin(), in_mus.end());
    double in_mu_MAX = *max_it;
    TH2D* Linearity_map = new TH2D("Linearity_map", ";input #mu;output #mu", 100, 0.0, in_mu_MAX, 100, 0.0, in_mu_MAX);
    for (int i = 0; i < in_mus.size(); i++) {
        Linearity_map->Fill(in_mus.at(i), out_mus.at(i));
    }
    TCanvas* c = new TCanvas("Linearity test canvas", "", 800, 800);
    Linearity_map->Draw("BOX");
    c->SaveAs("Linearity_MXs3_test_map.png");

    // calculate fitter bias
    double specific_out_mu = 0.0 * step;
    for (int i = 0; specific_out_mu < in_mu_MAX * 0.8; i++) {
        specific_out_mu = i * step;

        std::vector<double> in_mus_for_specific_out_mu;

        for (int j = 0; j < in_mus.size(); j++) {
            if (((specific_out_mu + 0.1) > out_mus.at(j)) && ((specific_out_mu - 0.1) < out_mus.at(j))) in_mus_for_specific_out_mu.push_back(in_mus.at(j));
        }

        TH1D* temp_hist = new TH1D("temp_hist", ";input #mu;toys", 100, specific_out_mu - 10.0, specific_out_mu + 10.0);
        for (int j = 0; j < in_mus_for_specific_out_mu.size(); j++) temp_hist->Fill(in_mus_for_specific_out_mu.at(j));
        temp_hist->Draw("Hist");
        c->SaveAs(("Linearity_MXs3_test_bayesian_" + std::to_string(specific_out_mu) + ".png").c_str());

        delete temp_hist;

    }

    delete c;

}
