void ReadLinearityRootFile(){

    for (int i = 0; i < 51; i++) {
        double injected_mu = i * 0.2;
        std::string fname = "LT_result_" + std::to_string(injected_mu) + ".root";
        
        TFile* input_file = new TFile(fname.c_str(), "read");
        TTree* temp_tree = (TTree*)input_file->Get("LT_result");

        double out_mu;
        std::vector<double> out_mus;
        unsigned int Nentry = temp_tree->GetEntries();

        temp_tree->SetBranchAddress("mu", &out_mu);

        for (unsigned int j = 0; j < Nentry; j++) { // Fill
            temp_tree->GetEntry(j);
            out_mus.push_back(out_mu);
        }
        input_file->Close();



        RooRealVar  mu_roorealvar("mu_roorealvar", "mu", -100, 100);
        RooDataSet mu_RooDataSet("mu_RooDataSet", "mu_RooDataSet", RooArgSet(mu_roorealvar));
        for (int j = 0; j < Nentry; j++) {
            mu_roorealvar = out_mus.at(j);
            mu_RooDataSet.add(RooArgSet(mu_roorealvar));
        }
        RooRealVar gausmean("gausmean", "", SIGNAL_num, SIGNAL_num - ERR_nsig, SIGNAL_num + ERR_nsig);
        RooRealVar gauswidth("gauswidth", "", ERR_nsig, 0, 2 * ERR_nsig);
        RooGaussian gauss("gauss", "gauss", nsig_roorealvar, gausmean, gauswidth);
        RooRealVar nentry("nentry", "number of entries", 1000, 900, 1100);
        RooExtendPdf egauss("egauss", "extended gauss", gauss, nentry);
        RooFitResult* result_LT = egauss.fitTo(nsig_RooDataSet, Save());
    }

    /*
    gStyle->SetOptFit(11); //gStyle->SetStatH(0.05);
    TCanvas* c = new TCanvas("Linearity test canvas", "", 800, 800);
    TGraphErrors* gr = new TGraphErrors(LT_number, Inputnsig, outputnsig, Inputnsigerror, outputnsigerror);
    gr->SetMarkerStyle(21); gr->SetTitle(";input n_{sig};output n_{sig}");
    gr->Fit("pol1"); gr->Draw("AP");
    c->SaveAs("Linearity_test.png");
    */
}
