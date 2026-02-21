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

TH1D* Kstar_delta_weight;
TH1D* K0star_delta_weight;
double ReadWeightHist(TH1D* hist, double value);

void ReadInfo(const char* dirname, TH1D* hist) {

    std::vector<string> names;

    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        // Get Data
        TTree* tree_info = (TTree*)input_file->Get("Xs");

        int __experiment__;
        int __run__;
        unsigned int __event__;
        int __candidate__;
        int __ncandidates__;
        double s;

        double MXs = -1.0;

        // get event_info
        tree_info->SetBranchAddress("__experiment__", &__experiment__);
        tree_info->SetBranchAddress("__run__", &__run__);
        tree_info->SetBranchAddress("__event__", &__event__);
        tree_info->SetBranchAddress("__candidate__", &__candidate__);
        tree_info->SetBranchAddress("__ncandidates__", &__ncandidates__);
        tree_info->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &s);
        tree_info->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &MXs);

        printf("%lld entries...\n", tree_info->GetEntries());
        for (unsigned int j = 0; j < tree_info->GetEntries(); j++) { // Fill
            tree_info->GetEntry(j);
            hist->Fill(s * s);
        }
        input_file->Close();
    }

}

void CalculatemKstar(){

    const int Nbin = 100;
    double MIN_q2 = 0.0;
    double MAX_q2 = 21.2;

    const char* dir_nominal = "B2Kstar";
    const char* dir_fix = "B2Kstar_fix";

    TH1D* Kstar_nominal = new TH1D("Kstar_nominal", "Kstar_nominal", Nbin, MIN_q2, MAX_q2);
    TH1D* Kstar_fix = new TH1D("Kstar_fix", "Kstar_fix", Nbin, MIN_q2, MAX_q2);
    TH1D* Kstar_ratio = new TH1D("Kstar_fix", "Kstar_fix", Nbin, MIN_q2, MAX_q2);

    ReadInfo(dir_nominal, Kstar_nominal);
    ReadInfo(dir_fix, Kstar_fix);
    Kstar_ratio->Divide(Kstar_fix, Kstar_nominal);

    double factor = 1.0;

    Kstar_nominal->Scale(factor / Kstar_nominal->Integral(), "width");
    Kstar_fix->Scale(factor / Kstar_fix->Integral(), "width");

    Kstar_nominal->SetFillStyle(3004);
    Kstar_nominal->SetLineColor(kBlue);
    Kstar_nominal->SetFillColor(kBlue);

    Kstar_fix->SetFillStyle(3005);
    Kstar_fix->SetLineColor(kRed);
    Kstar_fix->SetFillColor(kRed);

    gStyle->SetOptStat(0);

    TCanvas* c_temp = new TCanvas("c", "", 1200, 800); c_temp->cd();
    Kstar_nominal->Draw("Hist"); Kstar_fix->Draw("HistSAME");
    TLegend* legend = gPad->BuildLegend(0.9, 0.9, 0.6, 0.6); legend->SetFillStyle(0); legend->SetLineWidth(0);
    c_temp->SaveAs("mKstar_VS_nominal.png");

    FILE* fp = fopen("Kstar_mKstar_uncertainty.txt", "w");
    fprintf(fp, "%d %lf %lf\n", Nbin, MIN_q2, MAX_q2);
    for (int i = 0; i < Nbin; i++) fprintf(fp, "%lf\n", Kstar_ratio->GetBinContent(i + 1));
    fclose(fp);
}
