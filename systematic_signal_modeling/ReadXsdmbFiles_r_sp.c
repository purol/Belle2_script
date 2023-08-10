// last update: 2021-10-13
// for Belle2 data

# define Nentry 200
# define MIN_Q2 0.0
# define MAX_Q2 17.64

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

void LetsFillq2(const char* dirname, TH1D* hist_q2) {

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

            hist_q2->Fill(invM * invM);

        }
        input_file->Close();

    }

}

void ReadXsdmbFiles_r_sp(){
    const char* Xsd_Hmb_dirname = "./Xsd_Hmb";
    const char* Xsd_Lmb_dirname = "./Xsd_Lmb";
    const char* Xsd_dirname = "./Xsd";

    TH1D* hist_Hmb = new TH1D("Q2_Hmb", ";q^{2} [GeV^{2}/c^{4}]; evt", Nentry, MIN_Q2, MAX_Q2);
    TH1D* hist_Lmb = new TH1D("Q2_Lmb", ";q^{2} [GeV^{2}/c^{4}]; evt", Nentry, MIN_Q2, MAX_Q2);
    TH1D* hist = new TH1D("Q2", ";q^{2} [GeV^{2}/c^{4}]; evt", Nentry, MIN_Q2, MAX_Q2);

    LetsFillq2(Xsd_Hmb_dirname, hist_Hmb);
    LetsFillq2(Xsd_Lmb_dirname, hist_Lmb);
    LetsFillq2(Xsd_dirname, hist);

    Double_t factor = 1.;
    hist_Hmb->Scale(factor / hist_Hmb->Integral(), "width");
    hist_Lmb->Scale(factor / hist_Lmb->Integral(), "width");
    hist->Scale(factor / hist->Integral(), "width");

    hist_Hmb->SetStats(false);

    hist_Hmb->SetLineColor(2);
    hist_Lmb->SetLineColor(3);
    hist->SetLineColor(4);

    double ymax_Hmb = hist_Hmb->GetMaximum();
    double ymax_Lmb = hist_Lmb->GetMaximum();
    double ymax = hist->GetMaximum();
    double real_ymax = 0;
    if (ymax_Hmb > ymax_Lmb) real_ymax = ymax_Hmb;
    else real_ymax = ymax_Lmb;
    if (real_ymax < ymax) real_ymax = ymax;

    TCanvas* c_temp = new TCanvas("c", "", 1200, 1200); c_temp->cd();
    hist_Hmb->SetMaximum(real_ymax * 1.1);
    TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(hist_Hmb,"m_{b} = 4.83 GeV/c^{2}","f");
    legend->AddEntry(hist_Lmb, "m_{b} = 4.53 GeV/c^{2}", "f");
    legend->AddEntry(hist, "m_{b} = 4.68 GeV/c^{2}", "f");
    legend->SetFillStyle(0);
    legend->SetLineWidth(0);
    hist_Hmb->Draw("Hist");
    hist_Lmb->Draw("SAME Hist");
    hist->Draw("SAME Hist");
    legend->Draw("SAME");
    c_temp->SaveAs("Xsd_mb.png");

    FILE* fp;
    fp = fopen("Xsd_Hmb_weight.txt","w");
    fprintf(fp,"%d %lf %lf\n", Nentry, MIN_Q2, MAX_Q2);
    for (int j = 0; j < Nentry; j++) {
        double value = hist->GetBinContent(j + 1);
        double Hmb_value = hist_Hmb->GetBinContent(j + 1);

        double weight = 0;
        if (value > 0) weight = Hmb_value / value;
        else weight = 0;

        fprintf(fp, "%lf\n", weight);
    }
    fclose(fp);

    fp = fopen("Xsd_Lmb_weight.txt", "w");
    fprintf(fp, "%d %lf %lf\n", Nentry, MIN_Q2, MAX_Q2);
    for (int j = 0; j < Nentry; j++) {
        double value = hist->GetBinContent(j + 1);
        double Lmb_value = hist_Lmb->GetBinContent(j + 1);

        double weight = 0;
        if (value > 0) weight = Lmb_value / value;
        else weight = 0;

        fprintf(fp, "%lf\n", weight);
    }
}
