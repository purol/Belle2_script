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

void PrintToyResult(){

    std::vector<string> names;
    const char* dirname = "./";
    load_files(dirname, &names, "TOY_result");

    // read TOY files
    std::vector<double> mu1_fitting;
    std::vector<double> mu2_fitting;
    std::vector<double> mu3_fitting;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        TTree* temp_tree = (TTree*)input_file->Get("TOY_result");

        double temp_mu1_fitting = -1;
        double temp_mu2_fitting = -1;
        double temp_mu3_fitting = -1;

        int temp_covQual = -1;
        int temp_status = -1;
        double temp_edm = -1;

        temp_tree->SetBranchAddress("mu_MXs1_value", &temp_mu1_fitting);
        temp_tree->SetBranchAddress("mu_MXs2_value", &temp_mu2_fitting);
        temp_tree->SetBranchAddress("mu_MXs3_value", &temp_mu3_fitting);

        temp_tree->SetBranchAddress("covQual", &temp_covQual);
        temp_tree->SetBranchAddress("status", &temp_status);
        temp_tree->SetBranchAddress("edm", &temp_edm);

        for (unsigned int j = 0; j < temp_tree->GetEntries(); j++) { // Fill
            temp_tree->GetEntry(j);
            mu1_fitting.push_back(temp_mu1_fitting);
            mu2_fitting.push_back(temp_mu2_fitting);
            mu3_fitting.push_back(temp_mu3_fitting);
        }

        input_file->Close();
    }

    // write in txt file
    FILE* fp = fopen("mu_fit_values.csv", "w");
    fprintf(fp, "mu1,mu2,mu3\n");
    for (int i = 0; i < mu2_fitting.size(); i++) fprintf(fp, "%lf,%lf,%lf\n", mu1_fitting.at(i), mu2_fitting.at(i), mu3_fitting.at(i));
    fclose(fp);

}
