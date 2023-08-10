// last update: 2021-10-13
// for Belle2 data

# define Nentry 1500
# define MIN_Q2 0.6
# define MAX_Q2 1.65
# define Nfine 1000000

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

void LetsFillMK0star(const char* dirname, TH1D* hist) {

    double invM = 0.0;
    double MXs = 0.0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_Xs = (TTree*)input_file->Get("Xs");

        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);
        tree_Xs->SetBranchAddress("daughter__bo1__cm__spdaughter__bo0__cm__spM__bc__bc", &MXs);

        printf("%lld entries...\n", tree_Xs->GetEntries());
        for (unsigned int j = 0; j < tree_Xs->GetEntries(); j++) { // Fill
            tree_Xs->GetEntry(j);

            hist->Fill(MXs);

        }
        input_file->Close();

    }

}

void ReadK0starFiles_r_sp(){
    const double MK0star_evtpdl = 0.89555;
    const double MK0star_Gamma_evtpdl = 0.01;

    const char* K0star_dirname = "./K0starnunu";

    TH1D* hist = new TH1D("MK0star", ";M_{K^{*}} [GeV/c^{2}]; evt", Nentry, MIN_Q2, MAX_Q2);
    TH1D* hist_delta = new TH1D("MK0star_delta", ";M_{K^{*}} [GeV/c^{2}]; evt", Nentry, MIN_Q2, MAX_Q2);

    LetsFillMK0star(K0star_dirname, hist);
    for (int i = 0; i < Nfine; i++) {
        double temp = MK0star_evtpdl - (MK0star_Gamma_evtpdl / 2.0) + i * MK0star_Gamma_evtpdl / Nfine;
        hist_delta->Fill(temp);
    }

    Double_t factor = 1.;
    hist->Scale(factor / hist->Integral(), "width");
    hist_delta->Scale(factor / hist_delta->Integral(), "width");

    hist_delta->SetStats(false);

    hist_delta->SetLineColor(2);
    hist->SetLineColor(4);

    double ymax_delta = hist_delta->GetMaximum();
    double ymax = hist->GetMaximum();
    double real_ymax = 0;
    if (ymax_delta > ymax) real_ymax = ymax_delta;
    else real_ymax = ymax;

    TCanvas* c_temp = new TCanvas("c", "", 1200, 1200); c_temp->cd();
    hist_delta->SetMaximum(real_ymax * 1.1);
    TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(hist_delta,"#Gamma = 0.01 GeV/c^{2}","f");
    legend->AddEntry(hist, "evt.pdl", "f");
    legend->SetFillStyle(0);
    legend->SetLineWidth(0);
    hist_delta->Draw("Hist");
    hist->Draw("SAME Hist");
    legend->Draw("SAME");
    c_temp->SaveAs("MK0star.png");

    FILE* fp;
    fp = fopen("MK0star_weight.txt","w");
    fprintf(fp,"%d %lf %lf\n", Nentry, MIN_Q2, MAX_Q2);
    for (int j = 0; j < Nentry; j++) {
        double value = hist->GetBinContent(j + 1);
        double delta_value = hist_delta->GetBinContent(j + 1);

        double weight = 0;
        if (value > 0) weight = delta_value / value;
        else weight = 0;

        fprintf(fp, "%lf\n", weight);
    }
    fclose(fp);

}
