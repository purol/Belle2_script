void ReadToyRootFile(){

    std::string fname = "TOY_result.root";

    TFile* input_file = new TFile(fname.c_str(), "read");
    TTree* temp_tree = (TTree*)input_file->Get("TOY_result");

    TH1F* ToyMCmu = new TH1F("ToyMCmu", ";#mu;Toys", 40, -10, 10);
    TH1F* ToyMCmuerror = new TH1F("ToyMCmuerror", ";error of #mu;Toys", 50, 1, 4);
    TH1F* ToyMCmupull = new TH1F("ToyMCmupull", ";pull of #mu;Toys", 40, -4, 4);

    ToyMCmu->SetMarkerStyle(kFullCircle);
    ToyMCmu->SetLineColor(kBlack);
    ToyMCmu->SetMarkerColor(kBlack);

    ToyMCmuerror->SetMarkerStyle(kFullCircle);
    ToyMCmuerror->SetLineColor(kBlack);
    ToyMCmuerror->SetMarkerColor(kBlack);

    ToyMCmupull->SetMarkerStyle(kFullCircle);
    ToyMCmupull->SetLineColor(kBlack);
    ToyMCmupull->SetMarkerColor(kBlack);

    {
        double temp_mu = -1;
        double temp_muerror = -1;
        double temp_mupull = -1;

        temp_tree->SetBranchAddress("mu", &temp_mu);
        temp_tree->SetBranchAddress("mu_error", &temp_muerror);
        temp_tree->SetBranchAddress("mu_pull", &temp_mupull);

        for (unsigned int j = 0; j < temp_tree->GetEntries(); j++) { // Fill
            temp_tree->GetEntry(j);
            ToyMCmu->Fill(temp_mu);
            ToyMCmuerror->Fill(temp_muerror);
            ToyMCmupull->Fill(temp_mupull);
        }
    }

    TCanvas* c = new TCanvas("canvas_ToyMC_study", "", 800, 800);
    ToyMCmu->Draw("PE1");
    c->SaveAs("TOYMC_mu.png");
    delete c;

    c = new TCanvas("canvas_ToyMC_study", "", 800, 800);
    ToyMCmuerror->Draw("PE1");
    c->SaveAs("TOYMC_muerror.png");
    delete c;

    c = new TCanvas("canvas_ToyMC_study", "", 800, 800);
    ToyMCmupull->Fit("gaus");
    ToyMCmupull->Draw("PE1");
    c->SaveAs("TOYMC_mupull.png");
    delete c;

    TF1* fit_pull_gauss = ToyMCmupull->GetFunction("gaus");
    const double constant_gauss = fit_pull_gauss->GetParameter(0);
    const double mean_gauss = fit_pull_gauss->GetParameter(1);
    const double sigma_gauss = fit_pull_gauss->GetParameter(2);

    delete ToyMCmu;
    delete ToyMCmuerror;
    delete ToyMCmupull;

}
