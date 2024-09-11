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

#include "constants.h"

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

THStack* Stack = new THStack("thstack", ";dM [GeV];Candidates");
TH1F* SIGNAL_hist = new TH1F("signal", ";dM [GeV];Candidates", 100, -0.05, 0.05);
TH1F* CHG_hist = new TH1F("charge", ";E_{ecl} [GeV];Candidates", 100, -0.05, 0.05);
TH1F* MIX_hist = new TH1F("mix", ";E_{ecl} [GeV];Candidates", 100, -0.05, 0.05);
TH1F* UUBAR_hist = new TH1F("uubar", ";E_{ecl} [GeV];Candidates", 100, -0.05, 0.05);
TH1F* DDBAR_hist = new TH1F("ddbar", ";E_{ecl} [GeV];Candidates", 100, -0.05, 0.05);
TH1F* SSBAR_hist = new TH1F("ssbar", ";E_{ecl} [GeV];Candidates", 100, -0.05, 0.05);
TH1F* CHARM_hist = new TH1F("ccbar", ";E_{ecl} [GeV];Candidates", 100, -0.05, 0.05);

void LetsFill(const char* dirname, TH1F* hist, double weight = 1) {
    double var = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon_Jpsi");

        tree_upsilon->SetBranchAddress("daughter__bo1__cm__spdaughter__bo1__cm__spdM__bc__bc", &var);

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

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon_Jpsi");

        tree_upsilon->SetBranchAddress("daughter__bo1__cm__spdaughter__bo1__cm__spdM__bc__bc", &var);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            hist->Fill(var, weight);
        }
        input_file->Close();

    }

}

void THStack_dm_Jpsi() {

    const char* SIGNAL_dirname = "/home/jwpark/storage/BKG_gbasf2/Nazrin_LS_MC_Jpsi/Jpsi/SIGNAL";
    const char* CHG_dirname = "/home/jwpark/storage/BKG_gbasf2/Nazrin_LS_MC_Jpsi/Jpsi/CHG";
    const char* MIX_dirname = "/home/jwpark/storage/BKG_gbasf2/Nazrin_LS_MC_Jpsi/Jpsi/MIX";
    const char* UUBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Nazrin_LS_MC_Jpsi/Jpsi/UUBAR";
    const char* DDBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Nazrin_LS_MC_Jpsi/Jpsi/DDBAR";
    const char* SSBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Nazrin_LS_MC_Jpsi/Jpsi/SSBAR";
    const char* CHARM_dirname = "/home/jwpark/storage/BKG_gbasf2/Nazrin_LS_MC_Jpsi/Jpsi/CHARM";

    LetsFill(SIGNAL_dirname, SIGNAL_hist, (0.3616 / 1.8));
    LetsFill(CHG_dirname, CHG_hist, (0.3616 / 0.8));
    LetsFill(MIX_dirname, MIX_hist, (0.3616 / 0.8));
    LetsFill(UUBAR_dirname, UUBAR_hist, (0.3616 / 1.0));
    LetsFill(DDBAR_dirname, DDBAR_hist, (0.3616 / 1.0));
    LetsFill(SSBAR_dirname, SSBAR_hist, (0.3616 / 1.0));
    LetsFill(CHARM_dirname, CHARM_hist, (0.3616 / 1.0));


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

    TCanvas* c_temp = new TCanvas("c", "", 1200, 1200); c_temp->cd();

    gStyle->SetPalette(kGistEarth);

    Float_t ymax = SIGNAL_hist->GetMaximum(); Stack->SetMaximum(ymax*1.1);
    Stack->Draw("pfc Hist"); SIGNAL_hist->Draw("HistSAME");

    TLegend* legend = gPad->BuildLegend(0.9, 0.9, 0.7, 0.7);
    legend->SetFillStyle(0);

    c_temp->SaveAs("Plot_dM_Jpsi.png");

//    free(line);
//    gPad->BuildLegend();
//    c_temp->SaveAs("Plot_with_legend.png");
}
