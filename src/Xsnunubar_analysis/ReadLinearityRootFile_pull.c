#define LT_number 10

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

void ReadLinearityRootFile_pull(){

    std::vector<double> mus;
    std::vector<double> mus_err;

    double Inputmu[LT_number] = { 0 };
    double Inputmuerror[LT_number] = { 0 };
    double outputmu[LT_number] = { 0 };
    double outputmuerror[LT_number] = { 0 };

    for (int i = 0; i < LT_number; i++) {
        double injected_mu = i * 1.0;
        Inputmu[i] = injected_mu;

        std::vector<double> out_pulls;
        unsigned int Nentry = 0;

        std::vector<string> names;
        const char* dirname = "./";
        load_files(dirname, &names, ("LT_result_" + std::to_string(injected_mu)).c_str() );
        for (unsigned int k = 0; k < names.size(); k++) {

            TFile* input_file = new TFile((dirname + std::string("/") + names.at(k)).c_str(), "read");
            TTree* temp_tree = (TTree*)input_file->Get("LT_result");

            double temp_mu_true = -1;
            double temp_mu_fitting = -1;
            double temp_mu_error = -1;
            double temp_mu_HIerror = -1;
            double temp_mu_LOerror = -1;
            Nentry = temp_tree->GetEntries();

            temp_tree->SetBranchAddress("mu_MXs3_true", &temp_mu_true);
            temp_tree->SetBranchAddress("mu_MXs3_value", &temp_mu_fitting);
            temp_tree->SetBranchAddress("mu_MXs3_HIerror", &temp_mu_HIerror);
            temp_tree->SetBranchAddress("mu_MXs3_LOerror", &temp_mu_LOerror);

            for (unsigned int j = 0; j < Nentry; j++) { // Fill
                temp_tree->GetEntry(j);
                if (temp_mu_true >= temp_mu_fitting) out_pulls.push_back((temp_mu_true - temp_mu_fitting) / temp_mu_HIerror);
                else out_pulls.push_back((temp_mu_fitting - temp_mu_true) / temp_mu_LOerror);
            }
            input_file->Close();
	}

        RooRealVar  pull_roorealvar("pull_roorealvar", "mu", -4.0, 4.0);
        RooDataSet pull_RooDataSet("pull_RooDataSet", "pull_RooDataSet", RooArgSet(pull_roorealvar));
        for (int j = 0; j < out_pulls.size(); j++) {
            pull_roorealvar = out_pulls.at(j);
            pull_RooDataSet.add(RooArgSet(pull_roorealvar));
        }
        RooRealVar gausmean("mean", "", 0.0, -1.0, 1.0);
        RooRealVar gauswidth("width", "", 1.0, 0.8, 1.2);
        RooGaussian gauss("gauss", "gauss", pull_roorealvar, gausmean, gauswidth);
        RooRealVar nentry("nentry", "number of entries", 10000, 9500, 10500);
        RooExtendPdf egauss("egauss", "extended gauss", gauss, nentry);
        RooFitResult* result_LT = egauss.fitTo(pull_RooDataSet, RooFit::Save(), RooFit::Minimizer("Minuit2", "migrad"));

        RooArgSet fitargs_LT = result_LT->floatParsFinal();
        TIterator* iter_LT(fitargs_LT.createIterator());

        for (TObject* a_LT = iter_LT->Next(); a_LT != 0; a_LT = iter_LT->Next()) {
            RooRealVar* rrv_LT = dynamic_cast<RooRealVar*>(a_LT);
            std::string name_LT = rrv_LT->GetName();
            double val_LT = rrv_LT->getVal();
            double err_LT = rrv_LT->getError();
            if (name_LT == std::string("mean")) {
                outputmu[i] = val_LT;
                outputmuerror[i] = err_LT;
            }
        }

        RooPlot* pullframe = pull_roorealvar.frame(RooFit::Bins(20), RooFit::Title(" "));
        pull_RooDataSet.plotOn(pullframe);
        egauss.plotOn(pullframe, RooFit::LineColor(kBlue));
        egauss.paramOn(pullframe, RooFit::Layout(0.55, 0.9, 0.9)); pullframe->getAttText()->SetTextSize(0.03);
        TCanvas* c = new TCanvas("mu_gauss_fit", "mu_gauss_fit", 600, 600);
        gPad->SetLeftMargin(0.15); pullframe->GetYaxis()->SetTitleOffset(1.4); pullframe->Draw(); c->SaveAs(("mu_MXs3_" + to_string(injected_mu) + "_pull_distribution.png").c_str());

    }

    gStyle->SetOptFit(11); //gStyle->SetStatH(0.05);
    gStyle->SetStatX(0.6);  gStyle->SetStatY(0.9);
    TCanvas* c = new TCanvas("Linearity test canvas", "", 800, 800);
    TGraphErrors* gr = new TGraphErrors(LT_number, Inputmu, outputmu, Inputmuerror, outputmuerror);
    gr->SetMarkerStyle(21); gr->SetTitle(";input #mu;mean of pull");
    gr->Fit("pol1"); gr->Draw("AP");
    c->SaveAs("Linearity_MXs3_test_pull.png");

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
