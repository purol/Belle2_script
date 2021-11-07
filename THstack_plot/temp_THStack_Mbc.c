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
TH1F* SIGNAL_hist = new TH1F("SIGNAL_hist", ";#theta_{missing};",100, 5.24, 5.3);
TH1F* UUBAR_hist = new TH1F("UUBAR_hist", ";#theta_{missing};", 100, 5.24, 5.3);
TH1F* DDBAR_hist = new TH1F("DDBAR_hist", ";#theta_{missing};", 100, 5.24, 5.3);
TH1F* SSBAR_hist = new TH1F("SSBAR_hist", ";#theta_{missing};", 100, 5.24, 5.3);
TH1F* CHARM_hist = new TH1F("CHARM_hist", ";#theta_{missing};", 100, 5.24, 5.3);

void temp_THStack(){

    double var;

    std::vector<string> SIGNAL_names;
    const char* SIGNAL_dirname = "/home/jwpark/storage/SIGNAL/Eecl3_test_withchi/before_Mbc_cut";
    const char* UUBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/BKG_total/UUBAR_analysis/Eecl3_test_withchi/before_Mbc_cut";
    const char* DDBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/BKG_total/DDBAR_analysis/Eecl3_test_withchi_exactly_half/before_Mbc_cut";
    const char* SSBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/BKG_total/SSBAR_analysis/Eecl3_test_withchi_exactly_half/before_Mbc_cut";
    const char* CHARM_dirname = "/home/jwpark/storage/BKG_gbasf2/BKG_total/CHARM_analysis/Eecl3_test_withchi/before_Mbc_cut";

    load_files(SIGNAL_dirname, &SIGNAL_names);

    for(unsigned int i = 0; i< SIGNAL_names.size(); i++){

        TFile *input_file = new TFile( (SIGNAL_dirname+std::string("/")+ SIGNAL_names.at(i)).c_str(),"read");
        printf("%s (%d/%zu)\n",("Read "+ SIGNAL_names.at(i) + "... ").c_str(), i, SIGNAL_names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_Btag->SetBranchAddress("Btag_Mbc", &var);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            SIGNAL_hist->Fill(var);
        }
        input_file->Close();

    }

    std::vector<string> UUBAR_names;

    load_files(UUBAR_dirname, &UUBAR_names);

    for (unsigned int i = 0; i < UUBAR_names.size(); i++) {

        TFile* input_file = new TFile((UUBAR_dirname + std::string("/") + UUBAR_names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + UUBAR_names.at(i) + "... ").c_str(), i, UUBAR_names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_Btag->SetBranchAddress("Btag_Mbc", &var);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            UUBAR_hist->Fill(var);
        }
        input_file->Close();

    }

    std::vector<string> DDBAR_names;

    load_files(DDBAR_dirname, &DDBAR_names);

    for (unsigned int i = 0; i < DDBAR_names.size(); i++) {

        TFile* input_file = new TFile((DDBAR_dirname + std::string("/") + DDBAR_names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + DDBAR_names.at(i) + "... ").c_str(), i, DDBAR_names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_Btag->SetBranchAddress("Btag_Mbc", &var);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            DDBAR_hist->Fill(var);
        }
        input_file->Close();

    }

    std::vector<string> SSBAR_names;

    load_files(SSBAR_dirname, &SSBAR_names);

    for (unsigned int i = 0; i < SSBAR_names.size(); i++) {

        TFile* input_file = new TFile((SSBAR_dirname + std::string("/") + SSBAR_names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + SSBAR_names.at(i) + "... ").c_str(), i, SSBAR_names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_Btag->SetBranchAddress("Btag_Mbc", &var);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            SSBAR_hist->Fill(var);
        }
        input_file->Close();

    }

    std::vector<string> CHARM_names;

    load_files(CHARM_dirname, &CHARM_names);

    for (unsigned int i = 0; i < CHARM_names.size(); i++) {

        TFile* input_file = new TFile((CHARM_dirname + std::string("/") + CHARM_names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + CHARM_names.at(i) + "... ").c_str(), i, CHARM_names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_Btag->SetBranchAddress("Btag_Mbc", &var);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            CHARM_hist->Fill(var);
        }
        input_file->Close();

    }

    double UUBAR_int = UUBAR_hist->Integral();
    double DDBAR_int = DDBAR_hist->Integral();
    double SSBAR_int = SSBAR_hist->Integral();
    double CHARM_int = CHARM_hist->Integral();
    double SIGNAL_int = SIGNAL_hist->Integral();

    double BKG_int = UUBAR_int + DDBAR_int + SSBAR_int + CHARM_int;

    UUBAR_hist->Scale(1.0/BKG_int, "width");
    DDBAR_hist->Scale(1.0/BKG_int, "width");
    SSBAR_hist->Scale(1.0/BKG_int, "width");
    CHARM_hist->Scale(1.0/BKG_int, "width");    

    Stack->Add(UUBAR_hist);
    Stack->Add(DDBAR_hist);
    Stack->Add(SSBAR_hist);
    Stack->Add(CHARM_hist);
    //Stack->SetMaximum(100.0);


    double norm_for_SIG = 0.003385;

    SIGNAL_hist->Scale(norm_for_SIG*600.0/BKG_int, "width");
    SIGNAL_hist->SetLineWidth(2);
    SIGNAL_hist->SetLineColor(2);
    SIGNAL_hist->SetFillStyle(0);

    TCanvas* c_temp = new TCanvas("c", "", 1500, 1200); c_temp->cd();
    gStyle->SetPalette(kOcean);

    Stack->Draw("pfc Hist"); SIGNAL_hist->Draw("HistSAME");
    c_temp->SaveAs("Plot.png");

    gPad->BuildLegend();
    c_temp->SaveAs("Plot_with_legend.png");
}
