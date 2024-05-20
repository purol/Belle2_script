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

void ReadLinearityRootFile_median(){

    std::vector<double> mus;
    std::vector<double> mus_err;

    double Inputmu[LT_number] = { 0 };
    double Inputmuerror[LT_number] = { 0 };
    double outputmu[LT_number] = { 0 };
    double outputmuerror[LT_number] = { 0 };

    for (int i = 0; i < LT_number; i++) {
        double injected_mu = i * 0.2;
        Inputmu[i] = injected_mu;

        std::vector<double> out_mus;
        unsigned int Nentry = 0;

        std::vector<string> names;
        const char* dirname = "./";
        load_files(dirname, &names, ("LT_result_" + std::to_string(injected_mu)).c_str() );
        for (unsigned int k = 0; k < names.size(); k++) {

            TFile* input_file = new TFile((dirname + std::string("/") + names.at(k)).c_str(), "read");
            TTree* temp_tree = (TTree*)input_file->Get("LT_result");

            double out_mu;
            Nentry = temp_tree->GetEntries();

            temp_tree->SetBranchAddress("mu_MXs3_value", &out_mu);

            for (unsigned int j = 0; j < Nentry; j++) { // Fill
                temp_tree->GetEntry(j);
                out_mus.push_back(out_mu);
            }
            input_file->Close();
        }

        RooRealVar  mu_roorealvar("mu_roorealvar", "mu", injected_mu - 10.0, injected_mu + 10.0);
        RooDataSet mu_RooDataSet("mu_RooDataSet", "mu_RooDataSet", RooArgSet(mu_roorealvar));
        for (int j = 0; j < out_mus.size(); j++) {
            mu_roorealvar = out_mus.at(j);
            mu_RooDataSet.add(RooArgSet(mu_roorealvar));
        }

        // get median
        std::sort(out_mus.begin(), out_mus.end());
        size_t size = out_mus.size();

        outputmu[i] = out_mus[size / 2];
        outputmuerror[i] = 0.0; // it is median so it is zero

        RooPlot* muframe = mu_roorealvar.frame(RooFit::Bins(400), RooFit::Title(" "));
        mu_RooDataSet.plotOn(muframe);
        TCanvas* c = new TCanvas("mu_gauss_fit", "mu_gauss_fit", 600, 600);
        gPad->SetLeftMargin(0.15); muframe->GetYaxis()->SetTitleOffset(1.4); muframe->Draw(); c->SaveAs(("mu_" + to_string(injected_mu) + "_distribution.png").c_str());

    }

    gStyle->SetOptFit(11); //gStyle->SetStatH(0.05);
    gStyle->SetStatX(0.6);  gStyle->SetStatY(0.9);
    TCanvas* c = new TCanvas("Linearity test canvas", "", 800, 800);
    TGraphErrors* gr = new TGraphErrors(LT_number, Inputmu, outputmu, Inputmuerror, outputmuerror);
    gr->SetMarkerStyle(21); gr->SetTitle(";input #mu;median of output #mu");
    gr->Fit("pol1"); gr->Draw("AP");
    c->SaveAs("Linearity_MXs3_test_median.png");

    /*
    gr->SetMinimum(-1);
    gr->SetMaximum(12.0);
    gr->Draw("AP");
    c->Modified(); c->Update(); c->SaveAs("Linearity_test_magnify.png");
    */

    delete gr;
    delete c;

    gStyle->SetOptFit(0); //gStyle->SetStatH();

}
