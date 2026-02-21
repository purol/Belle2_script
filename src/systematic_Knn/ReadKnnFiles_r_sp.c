// last update: 2021-10-13
// for Belle2 data

# define Nentry 50
# define MIN_invM 1.876
# define MAX_invM 4.782

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

void ReadKnnFiles_r_sp(){
    /*
    get correction factor for B+ --> K+ n nbar MC
    MC:   B+ --> K+ n nbar
    data: B0 --> K0 p pbar
    */
    
    std::vector<string> names;
    const char* dirname = "./Knn";

    load_files(dirname, &names);

    TH1D* weight_hist = new TH1D("Knn weight", ";weight;", 50, 0.0, 10.0);
    TH1D* Plot = new TH1D("M", ";M [GeV/c^{2}]; evt", Nentry, MIN_invM, MAX_invM);

    for(unsigned int i = 0; i<names.size(); i++) {

        TFile *input_file = new TFile( (dirname+std::string("/")+names.at(i)).c_str(),"read");
        printf("%s (%d/%zu)\n",("Read "+names.at(i) + "... ").c_str(), i, names.size());

        // Get Data
        TTree* tree_info = (TTree*)input_file->Get("Knn");

        int __experiment__;
        int __run__;
        int __event__;
        int __candidate__;
        int __ncandidates__;
        double invM;

        // get event_info
        tree_info->SetBranchAddress("__experiment__", &__experiment__);
        tree_info->SetBranchAddress("__run__", &__run__);
        tree_info->SetBranchAddress("__event__", &__event__);
        tree_info->SetBranchAddress("__candidate__", &__candidate__);
        tree_info->SetBranchAddress("__ncandidates__", &__ncandidates__);
        tree_info->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM);

        printf("%lld entries...\n", tree_info->GetEntries());
        for (unsigned int j = 0; j < tree_info->GetEntries(); j++) { // Fill
            tree_info->GetEntry(j);
            Plot->Fill(invM);
        }
        input_file->Close();
    }

    Double_t factor = 1.;
    Plot->Scale(factor / Plot->Integral(), "width");

    //Plot->Scale(factor / Plot->GetEntries());
    //Plot->SetMaximum(1.2);
    //Plot->SetMinimum(0.0);
    //Plot->SetMarkerColor(4);
    //Plot->SetMarkerStyle(21);
    //Plot->SetMarkerSize(2);

    // draw experimental data [Phys. Rev. D 76, 092004]
    // B0 --> K0 p pbar data
    Float_t Knn_bd_x[12] = { 1.988, 2.225, 2.475, 2.725, 3.0, 3.275, 3.5, 3.75, 4.025, 4.2685, 4.4835, 4.681 };
    Float_t Knn_bd_y[12] = { 160, 151, 80, 60, 48, -2, 3, 13, 41, 24, 49, 25 };

    // create the error arrays
    Float_t ex[12] = { 0.112, 0.125, 0.125, 0.125, 0.15, 0.125, 0.1, 0.15, 0.125, 0.1185, 0.0965, 0.101 };
    Float_t Knn_bd_y_plus[12] = { 50, 50, 40, 31, 23, 13, 22, 27, 25, 25, 29, 34 };
    Float_t Knn_bd_y_minus[12] = { 50, 50, 40, 31, 23, 13, 22, 27, 25, 25, 29, 34 };
    Float_t Knn_bd_y_sym[12] = { 14, 13, 7, 6, 8, 5, 8, 5, 10 ,9, 10, 6 };
    Float_t Knn_bd_y_plus_tot[12] = { 0.0 };
    Float_t Knn_bd_y_minus_tot[12] = { 0.0 };
    for (int j = 0; j < 12; j++) {
        Knn_bd_y_plus_tot[j] = std::sqrt(Knn_bd_y_plus[j] * Knn_bd_y_plus[j] + Knn_bd_y_sym[j] * Knn_bd_y_sym[j]);
        Knn_bd_y_minus_tot[j] = std::sqrt(Knn_bd_y_minus[j] * Knn_bd_y_minus[j] + Knn_bd_y_sym[j] * Knn_bd_y_sym[j]);
    }

    // calculate dBR/dM
    for (int j = 0; j < 12; j++) {
        Knn_bd_y[j] = Knn_bd_y[j] / (2 * ex[j]);
        Knn_bd_y_plus_tot[j] = Knn_bd_y_plus_tot[j] / (2 * ex[j]);
        Knn_bd_y_minus_tot[j] = Knn_bd_y_minus_tot[j] / (2 * ex[j]);
    }

    // normalization
    for (int j = 0; j < 12; j++) {
        Knn_bd_y[j] = Knn_bd_y[j] / 831.0;
        Knn_bd_y_plus_tot[j] = Knn_bd_y_plus_tot[j] / 831.0;
        Knn_bd_y_minus_tot[j] = Knn_bd_y_minus_tot[j] / 831.0;
    }


    // create the TGraphErrors and draw it
    TGraphAsymmErrors* gr = new TGraphAsymmErrors(12, Knn_bd_x, Knn_bd_y, ex, ex, Knn_bd_y_plus_tot, Knn_bd_y_minus_tot);
    gr->SetMarkerColor(1);
    gr->SetMarkerStyle(21);
    gr->SetLineWidth(2);
    gr->SetMinimum(0.0);
    gr->SetTitle(";M [GeV/c^{2}];dBR/dM [arbitrary unit]");

    // set x-axis range
    auto axis = gr->GetXaxis();
    axis->SetLimits(MIN_invM, MAX_invM);
    gr->GetYaxis()->SetRangeUser(-0.1, 1.3);

    // fit by exponential
    TF1* f = new TF1("f", "[0]*(x-[1])*(x-[1])+[2]", MIN_invM, MAX_invM);
    f->SetParameter(0, 0.3);
    f->SetParameter(1, 3.5);
    f->SetParameter(2, 0.0);
    gr->Fit("f", "", "", MIN_invM, MAX_invM);
    double Area = f->Integral(MIN_invM, MAX_invM);
    double First = f->GetParameter(0);
    double Second = f->GetParameter(1);

    TCanvas* c_temp = new TCanvas("c", "", 1200, 1200); c_temp->cd();
    TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(Plot,"MC","f");
    legend->AddEntry(gr,"[Phys. Rev. D 76, 092004]","lpfe");
    legend->SetFillStyle(0);
    legend->SetLineWidth(0);
    gr->Draw("AP");
    //Plot->Draw("E1 SAME");
    Plot->Draw("SAME Hist");
    legend->Draw("SAME");
    c_temp->SaveAs("invMdistribution_Knn.png");

    FILE* fp;
    fp = fopen("Knn_weight.txt","w");
    fprintf(fp,"%d %lf %lf\n", Nentry, MIN_invM, MAX_invM);
    for (int j = 0; j < Nentry; j++) {
        double PHSP_bin_center = Plot->GetBinCenter(j + 1);
        double PHSP_value = Plot->GetBinContent(j + 1);
        double FUNC_value = f->Eval(PHSP_bin_center);
        FUNC_value = FUNC_value / Area;

        double weight = 0;
        if (PHSP_value > 0) weight = FUNC_value / PHSP_value;
        else weight = 0;

        fprintf(fp, "%lf\n", weight);
        weight_hist->Fill(weight);
    }

    fclose(fp);

    weight_hist->Draw("Hist");
    c_temp->SaveAs("weights_Knn.png");
}
