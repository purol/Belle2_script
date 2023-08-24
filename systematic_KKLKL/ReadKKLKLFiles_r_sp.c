// last update: 2021-10-13
// for Belle2 data

# define Nentry 32

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

void ReadKKLKLFiles_r_sp(){
    /*
    get correction factor for B+ --> K+ KL0 KL0 MC
    MC:   B+ --> K+ KL0 KL0
    */
    
    std::vector<string> names;
    const char* dirname = "./B2KKLKL";

    load_files(dirname, &names);

    TH1D* Plot = new TH1D("M_{KS0KS0}", ";M_{K_{S}^{0}K_{S}^{0}} [GeV/c^{2}]; evt", Nentry, 4.77 - 0.06 * 63, 4.77 + 0.06);
    TH1D* PHSP = new TH1D("M_{KS0KS0}PHSP", ";M_{K_{S}^{0}K_{S}^{0}} [GeV/c^{2}]; evt", Nentry, 4.77 - 0.06 * 63, 4.77 + 0.06);

    // PhysRevD.85.112010
    // use Engauge dizitizer
    const double FitValue[Nentry] = {
        76.5,
        53.6,
        63.3,
        82.3,
        73.5,
        38.6,
        13.3,
        31.6,
        24.1,
        23.3,
        20.4,
        16.4,
        14.7,
        12.0,
        11.0,
        9.4,
        8.0,
        7.8,
        7.6,
        12.3,
        34.1,
        1.9,
        2.2,
        1.8,
        1.2,
        0.7,
        0.2,
        0.1,
        0.1,
        0.4,
        0.8,
        0.7
    };

    for (int i = 0; i < Nentry; i++) Plot->SetBinContent(i + 1, FitValue[i]);

    for(unsigned int i = 0; i<names.size(); i++) {

        TFile *input_file = new TFile( (dirname+std::string("/")+names.at(i)).c_str(),"read");
        printf("%s (%d/%zu)\n",("Read "+names.at(i) + "... ").c_str(), i, names.size());

        // Get Data
        TTree* tree_info = (TTree*)input_file->Get("KKLKL");

        int __experiment__;
        int __run__;
        unsigned int __event__;
        int __candidate__;
        int __ncandidates__;

        int NKKLKLNR = -1;
        double invM;

        // get event_info
        tree_info->SetBranchAddress("__experiment__", &__experiment__);
        tree_info->SetBranchAddress("__run__", &__run__);
        tree_info->SetBranchAddress("__event__", &__event__);
        tree_info->SetBranchAddress("__candidate__", &__candidate__);
        tree_info->SetBranchAddress("__ncandidates__", &__ncandidates__);
        tree_info->SetBranchAddress("nParticlesInList__boB__pl__clKKLKL_NR__bc", &NKKLKLNR);
        tree_info->SetBranchAddress("invMassInLists__boK_L0__clfrom_KKLKL_NR__bc", &invM);

        printf("%lld entries...\n", tree_info->GetEntries());
        for (unsigned int j = 0; j < tree_info->GetEntries(); j++) { // Fill
            tree_info->GetEntry(j);
            if(NKKLKLNR == 1) PHSP->Fill(invM);
        }
        input_file->Close();
    }

    Double_t factor = 1.;
    Plot->Scale(factor / Plot->Integral(), "width");
    PHSP->Scale(factor / PHSP->Integral(), "width");

    Plot->SetFillStyle(3004);
    Plot->SetLineColor(kBlue);
    Plot->SetFillColor(kBlue);

    PHSP->SetFillStyle(3005);
    PHSP->SetLineColor(kRed);
    PHSP->SetFillColor(kRed);

    Plot->SetStats(false);
    PHSP->SetStats(false);

    //Plot->Scale(factor / Plot->GetEntries());
    //Plot->SetMaximum(1.2);
    //Plot->SetMinimum(0.0);
    //Plot->SetMarkerColor(4);
    //Plot->SetMarkerStyle(21);
    //Plot->SetMarkerSize(2);

    TCanvas* c_temp = new TCanvas("c", "", 1200, 1200); c_temp->cd();
    double y_max = -1;
    if (Plot->GetMaximum() > PHSP->GetMaximum()) y_max = Plot->GetMaximum();
    else y_max = PHSP->GetMaximum();
    PHSP->SetMaximum(1.05 * y_max);
    TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(PHSP,"Phase space","f");
    legend->AddEntry(Plot,"[Phys. Rev. D 85, 112010] FIG. 12","f");
    legend->SetFillStyle(0);
    legend->SetLineWidth(0);
    PHSP->Draw("Hist");
    Plot->Draw("SAME Hist");
    legend->Draw("SAME");
    c_temp->SaveAs("invMdistribution_KKLKL.png");

    FILE* fp;
    fp = fopen("KKLKL_weight.txt","w");
    fprintf(fp,"%d %lf %lf\n", Nentry, 4.77 - 0.06 * 63, 4.77 + 0.06);
    for (int j = 0; j < Nentry; j++) {
        double PHSP_value = PHSP->GetBinContent(j + 1);
        double Fit_value = Plot->GetBinContent(j + 1);

        double weight = 0;
        if (PHSP_value > 0) weight = Fit_value / PHSP_value;
        else weight = 0;

        fprintf(fp, "%lf\n", weight);
    }

    fclose(fp);
}
