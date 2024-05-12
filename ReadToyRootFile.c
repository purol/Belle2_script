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

void ReadToyRootFile(){

    TH1F* ToyMCmu = new TH1F("ToyMCmu", ";#mu;Toys", 40, -20, 25);
    TH1F* ToyMCmuerror = new TH1F("ToyMCmuerror", ";error of #mu;Toys", 50, 1, 7);
    TH1F* ToyMCmuHIerror = new TH1F("ToyMCmuHIerror", ";error of #mu;Toys", 50, 5, 15);
    TH1F* ToyMCmuLOerror = new TH1F("ToyMCmuLOerror", ";error of #mu;Toys", 50, -15, -5);
    TH1F* ToyMCmupull = new TH1F("ToyMCmupull", ";pull of #mu;Toys", 40, -4, 4);

    ToyMCmu->SetMarkerStyle(kFullCircle);
    ToyMCmu->SetLineColor(kBlack);
    ToyMCmu->SetMarkerColor(kBlack);

    ToyMCmuerror->SetMarkerStyle(kFullCircle);
    ToyMCmuerror->SetLineColor(kBlack);
    ToyMCmuerror->SetMarkerColor(kBlack);

    ToyMCmuHIerror->SetMarkerStyle(kFullCircle);
    ToyMCmuHIerror->SetLineColor(kBlack);
    ToyMCmuHIerror->SetMarkerColor(kBlack);

    ToyMCmuLOerror->SetMarkerStyle(kFullCircle);
    ToyMCmuLOerror->SetLineColor(kBlack);
    ToyMCmuLOerror->SetMarkerColor(kBlack);

    ToyMCmupull->SetMarkerStyle(kFullCircle);
    ToyMCmupull->SetLineColor(kBlack);
    ToyMCmupull->SetMarkerColor(kBlack);

    std::vector<string> names;
    const char* dirname = "./";
    load_files(dirname, &names, "TOY_result");
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        TTree* temp_tree = (TTree*)input_file->Get("TOY_result");

        double temp_mu_true = -1;
        double temp_mu_fitting = -1;
        double temp_mu_error = -1;
        double temp_mu_HIerror = -1;
        double temp_mu_LOerror = -1;

        int temp_covQual = -1;
        int temp_status = -1;
        double temp_edm = -1;

        temp_tree->SetBranchAddress("mu_MXs3_true", &temp_mu_true);
        temp_tree->SetBranchAddress("mu_MXs3_value", &temp_mu_fitting);
        temp_tree->SetBranchAddress("mu_MXs3_error", &temp_mu_error);
        temp_tree->SetBranchAddress("mu_MXs3_HIerror", &temp_mu_HIerror);
        temp_tree->SetBranchAddress("mu_MXs3_LOerror", &temp_mu_LOerror);

        temp_tree->SetBranchAddress("covQual", &temp_covQual);
        temp_tree->SetBranchAddress("status", &temp_status);
        temp_tree->SetBranchAddress("edm", &temp_edm);

        for (unsigned int j = 0; j < temp_tree->GetEntries(); j++) { // Fill
            temp_tree->GetEntry(j);
            ToyMCmu->Fill(temp_mu_fitting);
            ToyMCmuerror->Fill(temp_mu_error);
            ToyMCmuHIerror->Fill(temp_mu_HIerror);
            ToyMCmuLOerror->Fill(temp_mu_LOerror);
            if (temp_mu_true >= temp_mu_fitting) ToyMCmupull->Fill((temp_mu_true - temp_mu_fitting) / temp_mu_HIerror);
            else ToyMCmupull->Fill((temp_mu_fitting - temp_mu_true) / temp_mu_LOerror);
        }
        input_file->Close();
    }

    gStyle->SetOptFit(11);

    TCanvas* c = new TCanvas("canvas_ToyMC_study", "", 800, 800);
    ToyMCmu->Draw("PE1");
    c->SaveAs("TOYMC_mu_MXs3.png");
    delete c;

    c = new TCanvas("canvas_ToyMC_study", "", 800, 800);
    ToyMCmuerror->Draw("PE1");
    c->SaveAs("TOYMC_mu_MXs3_error.png");
    delete c;

    c = new TCanvas("canvas_ToyMC_study", "", 800, 800);
    ToyMCmuHIerror->Draw("PE1");
    c->SaveAs("TOYMC_mu_MXs3_HIerror.png");
    delete c;

    c = new TCanvas("canvas_ToyMC_study", "", 800, 800);
    ToyMCmuLOerror->Draw("PE1");
    c->SaveAs("TOYMC_mu_MXs3_LOerror.png");
    delete c;

    c = new TCanvas("canvas_ToyMC_study", "", 800, 800);
    ToyMCmupull->Fit("gaus");
    ToyMCmupull->Draw("PE1");
    c->SaveAs("TOYMC_mu_MXs3_pull.png");
    delete c;

    TF1* fit_pull_gauss = ToyMCmupull->GetFunction("gaus");
    const double constant_gauss = fit_pull_gauss->GetParameter(0);
    const double mean_gauss = fit_pull_gauss->GetParameter(1);
    const double sigma_gauss = fit_pull_gauss->GetParameter(2);

    delete ToyMCmu;
    delete ToyMCmuerror;
    delete ToyMCmuHIerror;
    delete ToyMCmuLOerror;
    delete ToyMCmupull;

}
