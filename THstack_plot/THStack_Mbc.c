// last update: 2021-10-08
// for Belle2 data

/*
when you add new variables:
revise # define N_Needed_info ...
revise typedef struct data
revise void Loader::GetData(TFile* input_file)
revise void Loader::PrintRootFile(std::string output_name)
revise void Loader::PrintSeparateRootFile(std::string output_name)
revise void Loader::ConvertIntoSeparateDataFile(std::string output_name, double flag = 0)
*/

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

THStack* Stack = new THStack("thstack", ";M_{bc}^{tag} [GeV];arbitrary unit");
TH1F* SIGNAL_hist = new TH1F("SIGNAL_hist", ";M_{bc}^{tag} [GeV];arbitrary unit", 100, 5.24, 5.3);
TH1F* CHG_hist = new TH1F("CHG_hist", ";M_{bc}^{tag} [GeV];arbitrary unit", 100, 5.24, 5.3);
TH1F* MIX_hist = new TH1F("MIX_hist", ";M_{bc}^{tag} [GeV];arbitrary unit", 100, 5.24, 5.3);
TH1F* UUBAR_hist = new TH1F("UUBAR_hist", ";M_{bc}^{tag} [GeV];arbitrary unit", 100, 5.24, 5.3);
TH1F* DDBAR_hist = new TH1F("DDBAR_hist", ";M_{bc}^{tag} [GeV];arbitrary unit", 100, 5.24, 5.3);
TH1F* SSBAR_hist = new TH1F("SSBAR_hist", ";M_{bc}^{tag} [GeV];arbitrary unit", 100, 5.24, 5.3);
TH1F* CHARM_hist = new TH1F("CHARM_hist", ";M_{bc}^{tag} [GeV];arbitrary unit", 100, 5.24, 5.3);

void LetsFill(const char* dirname, TH1F* hist) {
    double var = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_Btag->SetBranchAddress("Btag_Mbc", &var);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            hist->Fill(var);
        }
        input_file->Close();

    }

}

void THStack_Mbc() {

    const char* SIGNAL_dirname = "/home/jwpark/storage/SIGNAL_Aqua/test_v000/before_Mbc_cut";
    const char* CHG_dirname = "/home/jwpark/storage/BKG_gbasf2/Aqua/CHG_analysis/test_v000/before_Mbc_cut";
    const char* MIX_dirname = "/home/jwpark/storage/BKG_gbasf2/Aqua/MIX_analysis/test_v000/before_Mbc_cut";
    const char* UUBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Aqua/UUBAR_analysis/test_v000/before_Mbc_cut";
    const char* DDBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Aqua/DDBAR_analysis/test_v000/before_Mbc_cut";
    const char* SSBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Aqua/SSBAR_analysis/test_v000/before_Mbc_cut";
    const char* CHARM_dirname = "/home/jwpark/storage/BKG_gbasf2/Aqua/CHARM_analysis/test_v000/before_Mbc_cut";

    LetsFill(SIGNAL_dirname, SIGNAL_hist);
    LetsFill(CHG_dirname, CHG_hist);
    LetsFill(MIX_dirname, MIX_hist);
    LetsFill(UUBAR_dirname, UUBAR_hist);
    LetsFill(DDBAR_dirname, DDBAR_hist);
    LetsFill(SSBAR_dirname, SSBAR_hist);
    LetsFill(CHARM_dirname, CHARM_hist);

    double CHG_int = CHG_hist->Integral();
    double MIX_int = MIX_hist->Integral();
    double UUBAR_int = UUBAR_hist->Integral();
    double DDBAR_int = DDBAR_hist->Integral();
    double SSBAR_int = SSBAR_hist->Integral();
    double CHARM_int = CHARM_hist->Integral();
    double SIGNAL_int = SIGNAL_hist->Integral();

    double BKG_int = CHG_int + MIX_int + UUBAR_int + DDBAR_int + SSBAR_int + CHARM_int;

    CHG_hist->Scale(1.0 / BKG_int, "width");
    MIX_hist->Scale(1.0 / BKG_int, "width");
    UUBAR_hist->Scale(1.0 / BKG_int, "width");
    DDBAR_hist->Scale(1.0 / BKG_int, "width");
    SSBAR_hist->Scale(1.0 / BKG_int, "width");
    CHARM_hist->Scale(1.0 / BKG_int, "width");

    Stack->Add(CHG_hist);
    Stack->Add(MIX_hist);
    Stack->Add(UUBAR_hist);
    Stack->Add(DDBAR_hist);
    Stack->Add(SSBAR_hist);
    Stack->Add(CHARM_hist);
    //Stack->SetMaximum(100.0);


    double norm_for_SIG = 0.003385;

    SIGNAL_hist->Scale(norm_for_SIG*600.0/BKG_int, "width");
    //SIGNAL_hist->Scale(1.0 / SIGNAL_int, "width");
    SIGNAL_hist->SetLineWidth(3);
    SIGNAL_hist->SetLineColor(2);
    SIGNAL_hist->SetFillStyle(0);

    TCanvas* c_temp = new TCanvas("c", "", 1500, 1200); c_temp->cd();
    gStyle->SetPalette(kGistEarth);

    Stack->Draw("pfc Hist"); SIGNAL_hist->Draw("HistSAME");
    c_temp->SaveAs("Plot.png");

    gPad->BuildLegend();
    c_temp->SaveAs("Plot_with_legend.png");
}
