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

// arXiv:1409.4557v2
# define TB0 1.5195 // (Table. 1)
# define TBp 1.6384 // (Table. 1)
# define BR_Kplus_nunubar 0.00000398 // (eq. 10)
# define BR_K0star_nunubar 0.00000919 // (eq. 11)
# define BR_K0_nunubar (BR_Kplus_nunubar*TB0/TBp) // under (eq. 15)
# define BR_Kplusstar_nunubar (BR_K0star_nunubar*TBp/TB0) // under (eq. 15)
# define BR_Xs_nunubar 0.000029 // (eq. 23)
# define BR_Xsu_nonresonant_nunubar (BR_Xs_nunubar - BR_Kplus_nunubar - BR_Kplusstar_nunubar)
# define BR_Xsd_nonresonant_nunubar (BR_Xs_nunubar - BR_K0_nunubar - BR_K0star_nunubar)

// https://confluence.desy.de/pages/viewpage.action?pageId=107054222
# define N_BpBp_1invab 565400000.0
# define N_B0B0_1invab 534600000.0

# define N_Kplus_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Kplus_nunubar)
# define N_Kplusstar_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Kplusstar_nunubar)
# define N_Xsu_nonresonant_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Xsu_nonresonant_nunubar)
# define N_K0_nunubar_1invab (2.0 * N_B0B0_1invab * BR_K0_nunubar)
# define N_K0star_nunubar_1invab (2.0 * N_B0B0_1invab * BR_K0star_nunubar)
# define N_Xsd_nunubar_1invab (2.0 * N_B0B0_1invab * BR_Xsd_nonresonant_nunubar)

// my MC sample number
# define N_Kplus_nunubar 10000000.0
# define N_K0_nunubar 10000000.0
# define N_Kplusstar_nunubar 10000000.0
# define N_K0star_nunubar 10000000.0
# define N_Xsu_nonresonant_nunubar 50000000.0
# define N_Xsd_nonresonant_nunubar 50000000.0

// scale factor for each MC sample
# define Scale_Kplus (N_Kplus_nunubar_1invab/N_Kplus_nunubar)
# define Scale_Kplusstar (N_Kplusstar_nunubar_1invab/N_Kplusstar_nunubar)
# define Scale_Xsu_nonresonant (N_Xsu_nonresonant_nunubar_1invab/N_Xsu_nonresonant_nunubar)
# define Scale_K0 (N_K0_nunubar_1invab/N_K0_nunubar)
# define Scale_K0star (N_K0star_nunubar_1invab/N_K0star_nunubar)
# define Scale_Xsd_nonresonant (N_Xsd_nunubar_1invab/N_Xsd_nonresonant_nunubar)

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

void load_files(const char* dirname, std::vector<string>* names, const char* included_string) {
    TSystemDirectory dir(dirname, dirname);
    TList* files = dir.GetListOfFiles();
    if (files) {
        TSystemFile* file;
        TString fname;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(".root") && fname.Contains(included_string)) {
                names->push_back(fname.Data());
            }
        }
    }
}

THStack* Stack = new THStack("thstack", ";M_{bc}^{sig} [GeV];Arbitrary unit");
TH1F* SIGNAL_hist = new TH1F("signal", ";M_{bc}^{sig} [GeV];Arbitrary unit", 100, 5.25, 5.29);
TH1F* CHG_hist = new TH1F("charge", ";E_{ecl} [GeV];Arbitrary unit", 100, 5.25, 5.29);
TH1F* MIX_hist = new TH1F("mix", ";E_{ecl} [GeV];Arbitrary unit", 100, 5.25, 5.29);
TH1F* UUBAR_hist = new TH1F("uubar", ";E_{ecl} [GeV];Arbitrary unit", 100, 5.25, 5.29);
TH1F* DDBAR_hist = new TH1F("ddbar", ";E_{ecl} [GeV];Arbitrary unit", 100, 5.25, 5.29);
TH1F* SSBAR_hist = new TH1F("ssbar", ";E_{ecl} [GeV];Arbitrary unit", 100, 5.25, 5.29);
TH1F* CHARM_hist = new TH1F("ccbar", ";E_{ecl} [GeV];Arbitrary unit", 100, 5.25, 5.29);

void LetsFill(const char* dirname, TH1F* hist, double weight = 1) {
    double var = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Jpsi");

        tree_upsilon->SetBranchAddress("daughter__bo1__cm__spMbc__bc", &var);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            hist->Fill(var, weight);
        }
        input_file->Close();

    }

}

void LetsFill(const char* dirname, TH1F* hist, const char* included_string, double weight = 1) {
    double var = 0;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("roeEextra__bocleanMask__bc", &var);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            hist->Fill(var, weight);
        }
        input_file->Close();

    }

}

void THStack_Mbc_Jpsi() {

    const char* SIGNAL_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori_Jpsi/Jpsi/SIGNAL_after_FEISKIM_Nitori";
    const char* CHG_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori_Jpsi/Jpsi/CHG_after_FEISKIM_Nitori";
    const char* MIX_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori_Jpsi/Jpsi/MIX_after_FEISKIM_Nitori";
    const char* UUBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori_Jpsi/Jpsi/UUBAR_after_FEISKIM_Nitori";
    const char* DDBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori_Jpsi/Jpsi/DDBAR_after_FEISKIM_Nitori";
    const char* SSBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori_Jpsi/Jpsi/SSBAR_after_FEISKIM_Nitori";
    const char* CHARM_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori_Jpsi/Jpsi/CHARM_after_FEISKIM_Nitori";

    LetsFill(SIGNAL_dirname, SIGNAL_hist, (0.3/20.0));
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

    //CHG_hist->Scale(1.0 / BKG_int, "width");
    //MIX_hist->Scale(1.0 / BKG_int, "width");
    //UUBAR_hist->Scale(1.0 / BKG_int, "width");
    //DDBAR_hist->Scale(1.0 / BKG_int, "width");
    //SSBAR_hist->Scale(1.0 / BKG_int, "width");
    //CHARM_hist->Scale(1.0 / BKG_int, "width");

    Stack->Add(CHG_hist);
    Stack->Add(MIX_hist);
    Stack->Add(UUBAR_hist);
    Stack->Add(DDBAR_hist);
    Stack->Add(SSBAR_hist);
    Stack->Add(CHARM_hist);
    //Stack->SetMaximum(100.0);

    //SIGNAL_hist->Scale(10000.0/BKG_int, "width");
    //SIGNAL_hist->Scale(1.0 / SIGNAL_int, "width");
    SIGNAL_hist->SetLineWidth(3);
    SIGNAL_hist->SetLineColor(2);
    SIGNAL_hist->SetFillStyle(0);

    TCanvas* c_temp = new TCanvas("c", "", 1500, 1200); c_temp->cd();

    gStyle->SetPalette(kGistEarth);

    Float_t ymax = SIGNAL_hist->GetMaximum(); Stack->SetMaximum(ymax*1.1);
    Stack->Draw("pfc Hist"); SIGNAL_hist->Draw("HistSAME");

    TLegend* legend = gPad->BuildLegend(0.9, 0.9, 0.7, 0.7);
    legend->SetFillStyle(0);

    c_temp->SaveAs("Plot_Mbc_Jpsi.png");

//    free(line);
//    gPad->BuildLegend();
//    c_temp->SaveAs("Plot_with_legend.png");
}
