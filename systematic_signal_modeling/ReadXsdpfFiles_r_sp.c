// last update: 2021-10-13
// for Belle2 data

# define Nentry 200
# define MIN_MXs 1.1
# define MAX_MXs 3.52

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

void LetsFillMXs(const char* dirname, TH1D* hist_q2) {

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

            hist_q2->Fill(MXs);

        }
        input_file->Close();

    }

}

void ReadXsdpfFiles_r_sp(){
    const char* Xsd_Hpf_dirname = "./Xsd_Hpf";
    const char* Xsd_Lpf_dirname = "./Xsd_Lpf";
    const char* Xsd_dirname = "./Xsd";

    TH1D* hist_Hpf = new TH1D("MXs_Hpf", ";M_{Xs} [GeV/c^{2}]; evt", Nentry, MIN_MXs, MAX_MXs);
    TH1D* hist_Lpf = new TH1D("MXs_Lpf", ";M_{Xs} [GeV/c^{2}]; evt", Nentry, MIN_MXs, MAX_MXs);
    TH1D* hist = new TH1D("MXs", ";M_{Xs} [GeV/c^{2}]; evt", Nentry, MIN_MXs, MAX_MXs);

    LetsFillMXs(Xsd_Hpf_dirname, hist_Hpf);
    LetsFillMXs(Xsd_Lpf_dirname, hist_Lpf);
    LetsFillMXs(Xsd_dirname, hist);

    Double_t factor = 1.;
    hist_Hpf->Scale(factor / hist_Hpf->Integral(), "width");
    hist_Lpf->Scale(factor / hist_Lpf->Integral(), "width");
    hist->Scale(factor / hist->Integral(), "width");

    hist_Hpf->SetStats(false);

    hist_Hpf->SetLineColor(2);
    hist_Lpf->SetLineColor(3);
    hist->SetLineColor(4);

    double ymax_Hpf = hist_Hpf->GetMaximum();
    double ymax_Lpf = hist_Lpf->GetMaximum();
    double ymax = hist->GetMaximum();
    double real_ymax = 0;
    if (ymax_Hpf > ymax_Lpf) real_ymax = ymax_Hpf;
    else real_ymax = ymax_Lpf;
    if (real_ymax < ymax) real_ymax = ymax;

    TCanvas* c_temp = new TCanvas("c", "", 1200, 1200); c_temp->cd();
    hist_Hpf->SetMaximum(real_ymax * 1.1);
    TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(hist_Hpf,"p_{F} = 0.497 GeV/c","f");
    legend->AddEntry(hist_Lpf, "p_{F} = 0.422 GeV/c", "f");
    legend->AddEntry(hist, "p_{F} = 0.461 GeV/c", "f");
    legend->SetFillStyle(0);
    legend->SetLineWidth(0);
    hist_Hpf->Draw("Hist");
    hist_Lpf->Draw("SAME Hist");
    hist->Draw("SAME Hist");
    legend->Draw("SAME");
    c_temp->SaveAs("Xsd_pf.png");

    FILE* fp;
    fp = fopen("Xsd_Hpf_weight.txt","w");
    fprintf(fp,"%d %lf %lf\n", Nentry, MIN_MXs, MAX_MXs);
    for (int j = 0; j < Nentry; j++) {
        double value = hist->GetBinContent(j + 1);
        double Hpf_value = hist_Hpf->GetBinContent(j + 1);

        double weight = 0;
        if (value > 0) weight = Hpf_value / value;
        else weight = 0;

        fprintf(fp, "%lf\n", weight);
    }
    fclose(fp);

    fp = fopen("Xsd_Lpf_weight.txt", "w");
    fprintf(fp, "%d %lf %lf\n", Nentry, MIN_MXs, MAX_MXs);
    for (int j = 0; j < Nentry; j++) {
        double value = hist->GetBinContent(j + 1);
        double Lpf_value = hist_Lpf->GetBinContent(j + 1);

        double weight = 0;
        if (value > 0) weight = Lpf_value / value;
        else weight = 0;

        fprintf(fp, "%lf\n", weight);
    }
}
