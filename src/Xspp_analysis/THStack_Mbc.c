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

THStack* Stack = new THStack("thstack", ";M_{bc} [GeV/c^{2}];number of candidates");
TH1D* SIGNAL_hist = new TH1D("signal", ";M_{bc} [GeV/c^{2}];number of candidates", 30, 5.20, 5.3);
TH1D* CHG_hist = new TH1D("charged", ";M_{bc} [GeV/c^{2}];number of candidates", 30, 5.20, 5.3);
TH1D* MIX_hist = new TH1D("mixed", ";M_{bc} [GeV/c^{2}];number of candidates", 30, 5.20, 5.3);
TH1D* UUBAR_hist = new TH1D("u#bar{u}", ";M_{bc} [GeV/c^{2}];number of candidates", 30, 5.20, 5.3);
TH1D* DDBAR_hist = new TH1D("d#bar{d}", ";M_{bc} [GeV/c^{2}];number of candidates", 30, 5.20, 5.3);
TH1D* SSBAR_hist = new TH1D("s#bar{s}", ";M_{bc} [GeV/c^{2}];number of candidates", 30, 5.20, 5.3);
TH1D* CHARM_hist = new TH1D("c#bar{c}", ";M_{bc} [GeV/c^{2}];number of candidates", 30, 5.20, 5.3);

void LetsFill(const char* dirname, TH1D* hist, double weight = 1, bool isSignal = false) {
    double var = 0;
    double nBplusToXpp_var = -1;
    double nBzeroToXpp_var = -1;
    double Mpp = -1;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("daughter__bo1__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &Mpp);
        tree_Bsig->SetBranchAddress("Bsig_Mbc", &var);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clXpp__bc", &nBplusToXpp_var);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clXpp__bc", &nBzeroToXpp_var);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            if(isSignal == true){
                if((names.at(i).find("CHG") == std::string::npos) && (names.at(i).find("MIX") == std::string::npos)) continue;
                else{
                    if((nBplusToXpp_var > 0.5) || (nBzeroToXpp_var > 0.5)) {}
                    else continue;
                }
            }
            else{
                if((names.at(i).find("CHG") == std::string::npos) && (names.at(i).find("MIX") == std::string::npos)) {}
                else{
                    if((nBplusToXpp_var > 0.5) || (nBzeroToXpp_var > 0.5)) continue;
                    else {}
                }
            }

            hist->Fill(var, weight);
        }
        input_file->Close();

    }

}

void LetsFill(const char* dirname, TH1D* hist, const char* included_string, double weight = 1, bool isSignal = false) {
    double var = 0;
    double nBplusToXpp_var = -1;
    double nBzeroToXpp_var = -1;
    double Mpp = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("daughter__bo1__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &Mpp);
        tree_Bsig->SetBranchAddress("Bsig_Mbc", &var);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clXpp__bc", &nBplusToXpp_var);
        tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clXpp__bc", &nBzeroToXpp_var);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            if(isSignal == true){
                if((names.at(i).find("CHG") == std::string::npos) && (names.at(i).find("MIX") == std::string::npos)) continue;
                else{
                    if((nBplusToXpp_var > 0.5) || (nBzeroToXpp_var > 0.5)) {}
                    else continue;
                }
            }
            else{
                if((names.at(i).find("CHG") == std::string::npos) && (names.at(i).find("MIX") == std::string::npos)) {}
                else{
                    if((nBplusToXpp_var > 0.5) || (nBzeroToXpp_var > 0.5)) continue;
                    else {}
                }
            }

            hist->Fill(var, weight);
        }
        input_file->Close();

    }

}

void THStack_Mbc() {

    const char* CHG_dirname = "/home/belle2/junewoo/storage_ghi/pp_SS_systematic_BCS_bug_fix/pp/pp_sample_output/CHG";
    const char* MIX_dirname = "/home/belle2/junewoo/storage_ghi/pp_SS_systematic_BCS_bug_fix/pp/pp_sample_output/MIX";
    const char* UUBAR_dirname = "/home/belle2/junewoo/storage_ghi/pp_SS_systematic_BCS_bug_fix/pp/pp_sample_output/UUBAR";
    const char* DDBAR_dirname = "/home/belle2/junewoo/storage_ghi/pp_SS_systematic_BCS_bug_fix/pp/pp_sample_output/DDBAR";
    const char* SSBAR_dirname = "/home/belle2/junewoo/storage_ghi/pp_SS_systematic_BCS_bug_fix/pp/pp_sample_output/SSBAR";
    const char* CHARM_dirname = "/home/belle2/junewoo/storage_ghi/pp_SS_systematic_BCS_bug_fix/pp/pp_sample_output/CHARM";

    LetsFill(CHG_dirname, SIGNAL_hist, "final_final", Scale_CHG_validation_MC15rd, true);
    LetsFill(MIX_dirname, SIGNAL_hist, "final_final", Scale_MIX_validation_MC15rd, true);

    LetsFill(CHG_dirname, CHG_hist, "final_final", Scale_CHG_validation_MC15rd, false);
    LetsFill(MIX_dirname, MIX_hist, "final_final", Scale_MIX_validation_MC15rd, false);
    LetsFill(UUBAR_dirname, UUBAR_hist, "final_final", Scale_UUBAR_validation_MC15rd, false);
    LetsFill(DDBAR_dirname, DDBAR_hist, "final_final", Scale_DDBAR_validation_MC15rd, false);
    LetsFill(SSBAR_dirname, SSBAR_hist, "final_final", Scale_SSBAR_validation_MC15rd, false);
    LetsFill(CHARM_dirname, CHARM_hist, "final_final", Scale_CHARM_validation_MC15rd, false);


    double CHG_int = CHG_hist->Integral();
    double MIX_int = MIX_hist->Integral();
    double UUBAR_int = UUBAR_hist->Integral();
    double DDBAR_int = DDBAR_hist->Integral();
    double SSBAR_int = SSBAR_hist->Integral();
    double CHARM_int = CHARM_hist->Integral();

    double BKG_int = CHG_int + MIX_int + UUBAR_int + DDBAR_int + SSBAR_int + CHARM_int;

    /*
    CHG_hist->Scale(1.0 / BKG_int, "width");
    MIX_hist->Scale(1.0 / BKG_int, "width");
    UUBAR_hist->Scale(1.0 / BKG_int, "width");
    DDBAR_hist->Scale(1.0 / BKG_int, "width");
    SSBAR_hist->Scale(1.0 / BKG_int, "width");
    CHARM_hist->Scale(1.0 / BKG_int, "width");
    */

    TGaxis::SetMaxDigits(3);

    Stack->Add(CHG_hist);
    Stack->Add(MIX_hist);
    Stack->Add(UUBAR_hist);
    Stack->Add(DDBAR_hist);
    Stack->Add(SSBAR_hist);
    Stack->Add(CHARM_hist);
    Stack->SetMaximum(25.0);

    //SIGNAL_hist->Scale(10000.0/BKG_int, "width");
    //SIGNAL_hist->Scale(1.0 / SIGNAL_int, "width");
    //SIGNAL_hist->Scale(20.0);
    SIGNAL_hist->SetLineWidth(3);
    SIGNAL_hist->SetLineColor(2);
    SIGNAL_hist->SetFillStyle(0);

    TCanvas* c_temp = new TCanvas("c", "", 1500, 1200);
    c_temp->SetCanvasSize(1500, 1200);
    c_temp->SetLeftMargin(0.17);
    c_temp->SetRightMargin(0.04);
    c_temp->SetBottomMargin(0.16);
    c_temp->SetTopMargin(0.05);
    c_temp->cd();

    gStyle->SetPalette(kGistEarth);
    gStyle->SetOptStat(0);

    Stack->SetMinimum(0.0);
    Stack->Draw("pfc Hist");
    Stack->GetXaxis()->SetTitle("M_{bc} [GeV/c^{2}]");
    Stack->GetYaxis()->SetTitle("Number of candidates");
    Stack->GetXaxis()->SetLabelFont(42);
    Stack->GetYaxis()->SetLabelFont(42);
    Stack->GetXaxis()->SetTitleFont(42);
    Stack->GetYaxis()->SetTitleFont(42);
    Stack->GetXaxis()->SetLabelSize(0.040);
    Stack->GetYaxis()->SetLabelSize(0.040);
    Stack->GetXaxis()->SetTitleSize(0.045);
    Stack->GetYaxis()->SetTitleSize(0.045);
    Stack->GetXaxis()->SetTitleOffset(1.25);
    Stack->GetYaxis()->SetTitleOffset(1.40);
    Stack->GetXaxis()->CenterTitle();
    Stack->GetYaxis()->CenterTitle();
    Stack->GetXaxis()->SetNdivisions(505);
    Stack->GetYaxis()->SetNdivisions(505);
    SIGNAL_hist->Draw("HistSAME");
    c_temp->Modified();
    c_temp->Update();

    TLegend* legend = new TLegend(0.23, 0.76, 0.94, 0.93);
    legend->SetNColumns(3);
    legend->SetTextFont(42);
    legend->SetTextSize(0.040);
    legend->SetMargin(0.18);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->AddEntry(CHG_hist, "charged", "f");
    legend->AddEntry(MIX_hist, "mixed", "f");
    legend->AddEntry(UUBAR_hist, "u#bar{u}", "f");
    legend->AddEntry(DDBAR_hist, "d#bar{d}", "f");
    legend->AddEntry(SSBAR_hist, "s#bar{s}", "f");
    legend->AddEntry(CHARM_hist, "c#bar{c}", "f");
    legend->AddEntry(SIGNAL_hist, "B #rightarrow X_{s}p#bar{p}", "l");
    legend->Draw();

    c_temp->RedrawAxis();
    c_temp->Modified();
    c_temp->Update();
    c_temp->SaveAs("Plot_Mbc.png");
    c_temp->Print("Plot_Mbc.pdf", "pdf Portrait");

//    free(line);
//    gPad->BuildLegend();
//    c_temp->SaveAs("Plot_with_legend.png");
}



