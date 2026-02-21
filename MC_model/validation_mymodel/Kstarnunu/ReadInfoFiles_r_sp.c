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

double ReadWeightHist(TH1D* hist, double value) {
    int Bin = hist->FindBin(value);
    if (Bin < 1) Bin = 1;
    else if (Bin > hist->GetNbinsX()) Bin = hist->GetNbinsX();
    return hist->GetBinContent(Bin);
}

void ReadSignalModelingFile() {

    FILE* fp;
    int Nentry = 0;
    double RangeMIN = -1;
    double RangeMAX = -1;

    // Kstar_delta_weight
    fp = fopen("/home/jwpark/storage/BKG_gbasf2/systematic/signal_modeling/MKstar_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &Nentry, &RangeMIN, &RangeMAX);
    Kstar_delta_weight = new TH1D("Kstar_delta_weight", ";;", Nentry, RangeMIN, RangeMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        Kstar_delta_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

    // K0star_delta_weight
    fp = fopen("/home/jwpark/storage/BKG_gbasf2/systematic/signal_modeling/MK0star_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &Nentry, &RangeMIN, &RangeMAX);
    K0star_delta_weight = new TH1D("K0star_delta_weight", ";;", Nentry, RangeMIN, RangeMAX);
    for (int i = 0; i < Nentry; i++) {
        double weight = 0;
        fscanf(fp, "%lf\n", &weight);
        K0star_delta_weight->SetBinContent(i + 1, weight);
    }
    fclose(fp);

}

void ReadInfoFiles_r_sp(){

    std::vector<string> names;
    const char* dirname = "./";

    load_files(dirname, &names);

    Float_t Lower[7];

    Lower[0] = 0;
    Lower[1] = 4;
    Lower[2] = 8;
    Lower[3] = 12;
    Lower[4] = 16;
    Lower[5] = 19.2;
    Lower[6] = 21.2;
    TH1F* Plot = new TH1F("q2", ";q^{2} [GeV^{2}]; evt", 6, Lower);

    for(unsigned int i = 0; i<names.size(); i++){

        TFile *input_file = new TFile( (dirname+std::string("/")+names.at(i)).c_str(),"read");
        printf("%s (%d/%zu)\n",("Read "+names.at(i) + "... ").c_str(), i, names.size());

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
            Plot->Fill(s*s, ReadWeightHist(Kstar_delta_weight, MXs));
        }
        input_file->Close();
    }

    Double_t factor = 9.19;
    Float_t ratio = factor / Plot->GetEntries();

    // draw my MC data
    double x_MC[6] = { 2, 6, 10, 14, 17.6, 20.2 };
    double y_MC[6] = { Plot->GetBinContent(1) * ratio, Plot->GetBinContent(2) * ratio, Plot->GetBinContent(3) * ratio, Plot->GetBinContent(4) * ratio, Plot->GetBinContent(5) * ratio, Plot->GetBinContent(6) * ratio };

    // create the statistical error for my MC data
    double ex_MC[6] = { 2, 2, 2, 2, 1.6, 1.0 };
    double ey_MC[6] = { std::sqrt(Plot->GetBinContent(1)) * ratio, std::sqrt(Plot->GetBinContent(2)) * ratio, std::sqrt(Plot->GetBinContent(3)) * ratio, std::sqrt(Plot->GetBinContent(4)) * ratio, std::sqrt(Plot->GetBinContent(5)) * ratio, std::sqrt(Plot->GetBinContent(6)) * ratio };

    TGraphErrors* gr_MC = new TGraphErrors(6, x_MC, y_MC, ex_MC, ey_MC);
    gr_MC->SetMarkerColor(4);
    gr_MC->SetMarkerStyle(21);

    Plot->Scale(factor / Plot->GetEntries());
    //Plot->SetMaximum(1.2);
    //Plot->SetMinimum(0.0);
    Plot->SetMarkerColor(4);
    Plot->SetMarkerStyle(21);
    //Plot->SetMarkerSize(2);

    // draw theoratical data
    Float_t x[5] = { 2, 6, 10, 14, 17.6 };
    Float_t y[5] = { 1.38, 1.88, 2.27, 2.36, 1.30 };

    // create the error arrays
    Float_t ex[5] = { 2, 2, 2, 2, 1.6 };
    Float_t ey[5] = { 0.28, 0.32, 0.34, 0.31, 0.17 };
    // uncertainty by form factor

    // create the TGraphErrors and draw it
    TGraphErrors* gr = new TGraphErrors(5, x, y, ex, ey);
    //gr->SetMarkerColor(4);
    //gr->SetMarkerStyle(21);
    gr->SetFillColor(2);
    gr->SetFillStyle(3001);
    gr->SetMinimum(0.0);
    gr->SetTitle(";q^{2} [GeV^{2}];arbitrary unit");

    TCanvas* c_temp = new TCanvas("c", "", 1200, 800); c_temp->cd();
    TLegend *legend = new TLegend(0.35, 0.35, 0.1, 0.1);
    legend->AddEntry(gr_MC,"MC","lpfe");
    legend->AddEntry(gr,"[arXiv:1409.4557]","f");
    gr->Draw("a2");
    //Plot->Draw("E1 SAME");
    gr_MC->Draw("SAMEP");
    legend->Draw("SAME");
    c_temp->SaveAs("q2distribution.png");
}
