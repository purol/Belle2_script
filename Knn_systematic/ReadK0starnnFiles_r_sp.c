// last update: 2021-10-13
// for Belle2 data

# define Nentry 500
# define MIN_QSQ 1.876
# define MAX_QSQ 4.79

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

void ReadK0starnnFiles_r_sp(){
    
    std::vector<string> names;
    const char* dirname = "./";

    load_files(dirname, &names);

    TH1D* Plot = new TH1D("M_{p#bar{p}}", ";M_{p#bar{p}} [GeV/c^{2}]; evt", Nentry, MIN_QSQ, MAX_QSQ);

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

    // draw experimental data [Phys. Rev. Lett. 100, 251801]
    Float_t Knn_bd_x[7] = { 1.938, 2.1, 2.3, 2.5, 2.725, 3.2215, 4.0675 };
    Float_t Knn_bd_y[7] = { 0.3, 0.31, 0.26, 0.22, 0.02, 0.05, 0.02 };

    // create the error arrays
    Float_t ex[7] = { 0.062, 0.1, 0.1, 0.1, 0.125, 0.0935, 0.3325 };
    Float_t Knn_bd_y_plus_tot[7] = { 0.11, 0.12, 0.10, 0.11, 0.09, 0.09, 0.14 };
    Float_t Knn_bd_y_minus_tot[7] = { 0.10, 0.11, 0.09, 0.07, 0.06, 0.11 };

    // calculate dBR/dM
    for (int j = 0; j < 7; j++) {
        Knn_bd_y[j] = Knn_bd_y[j] / (2 * ex[j]);
        Knn_bd_y_plus_tot[j] = Knn_bd_y_plus_tot[j] / (2 * ex[j]);
        Knn_bd_y_minus_tot[j] = Knn_bd_y_minus_tot[j] / (2 * ex[j]);
    }

    // normalization
    for (int j = 0; j < 7; j++) {
        Knn_bd_y[j] = Knn_bd_y[j] / 1.18;
        Knn_bd_y_plus_tot[j] = Knn_bd_y_plus_tot[j] / 1.18;
        Knn_bd_y_minus_tot[j] = Knn_bd_y_minus_tot[j] / 1.18;
    }


    // create the TGraphErrors and draw it
    TGraphAsymmErrors* gr = new TGraphAsymmErrors(7, Knn_bd_x, Knn_bd_y, ex, ex, Knn_bd_y_plus_tot, Knn_bd_y_minus_tot);
    gr->SetMarkerColor(1);
    gr->SetMarkerStyle(21);
    gr->SetLineWidth(2);
    gr->SetMinimum(0.0);
    gr->SetTitle(";M_{p#bar{p}} [GeV/c^{2}];dBR/dM_{p#bar{p}} [arbitrary unit]");

    // set x-axis range
    auto axis = gr->GetXaxis();
    axis->SetLimits(MIN_QSQ, MAX_QSQ);

    // fit by exponential
    TF1* f = new TF1("f", "[0]*exp(x*[1])", MIN_QSQ, MAX_QSQ);
    f->SetParameter(0, 1);
    f->SetParameter(1, -0.4);
    gr->Fit("f", "", "", MIN_QSQ, MAX_QSQ);
    double Area = f->Integral(MIN_QSQ, MAX_QSQ);
    double First = f->GetParameter(0);
    double Second = f->GetParameter(1);

    TCanvas* c_temp = new TCanvas("c", "", 1200, 1200); c_temp->cd();
    TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(Plot,"Phase space","f");
    legend->AddEntry(gr,"[Phys. Rev. Lett. 100, 251801]","lpfe");
    legend->SetFillStyle(0);
    legend->SetLineWidth(0);
    gr->Draw("AP");
    //Plot->Draw("E1 SAME");
    Plot->Draw("SAME Hist");
    legend->Draw("SAME");
    c_temp->SaveAs("invMdistribution_K0starnn.png");

    FILE* fp;
    fp = fopen("K0starnn_weight.txt","w");
    fprintf(fp,"%d %lf %lf\n", Nentry, MIN_QSQ, MAX_QSQ);
    for (int j = 0; j < Nentry; j++) {
        double PHSP_bin_center = Plot->GetBinCenter(j + 1);
        double PHSP_value = Plot->GetBinContent(j + 1);
        double FUNC_value = f->Eval(PHSP_bin_center);
        FUNC_value = FUNC_value / Area;

        double weight = 0;
        if (PHSP_value > 0) weight = FUNC_value / PHSP_value;
        else weight = 0;

        fprintf(fp, "%lf\n", weight);
    }

    fclose(fp);
}
