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

THStack* Stack = new THStack("thstack", ";#theta_{missing} [rad];number of candidates");
TH1D* SIGNAL_hist = new TH1D("signal #times 10000", ";#theta_{missing};number of candidates", 100, 0, 3.142);
TH1D* CHG_hist = new TH1D("charged", ";#theta_{missing};number of candidates", 100, 0, 3.142);
TH1D* MIX_hist = new TH1D("mixed", ";#theta_{missing};number of candidates", 100, 0, 3.142);
TH1D* UUBAR_hist = new TH1D("u#bar{u}", ";#theta_{missing};number of candidates", 100, 0, 3.142);
TH1D* DDBAR_hist = new TH1D("d#bar{d}", ";#theta_{missing};number of candidates", 100, 0, 3.142);
TH1D* SSBAR_hist = new TH1D("s#bar{s}", ";#theta_{missing};number of candidates", 100, 0, 3.142);
TH1D* CHARM_hist = new TH1D("c#bar{c}", ";#theta_{missing};number of candidates", 100, 0, 3.142);

void LetsFill(const char* dirname, TH1D* hist, double weight = 1) {
    double var = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("missingMomentumOfEvent_theta", &var);

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

void LetsFill(const char* dirname, TH1D* hist, const char* included_string, double weight = 1) {
    double var = 0;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("missingMomentumOfEvent_theta", &var);

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

void THStack_missing_theta() {

    const char* Knunu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/test_v000/before_missing_momentum_theta_cut";
    const char* Kstarnunu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/test_v000/before_missing_momentum_theta_cut";
    const char* Xsununu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/test_v000/before_missing_momentum_theta_cut";
    const char* K0nunu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/test_v000/before_missing_momentum_theta_cut";
    const char* K0starnunu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/test_v000/before_missing_momentum_theta_cut";
    const char* Xsdnunu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/test_v000/before_missing_momentum_theta_cut";
    const char* CHG_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/CHG_analysis/test_v000/before_missing_momentum_theta_cut";
    const char* MIX_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/MIX_analysis/test_v000/before_missing_momentum_theta_cut";
    const char* UUBAR_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/UUBAR_analysis/test_v000/before_before_missing_momentum_theta_cut";
    const char* DDBAR_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/DDBAR_analysis/test_v000/before_missing_momentum_theta_cut";
    const char* SSBAR_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SSBAR_analysis/test_v000/before_missing_momentum_theta_cut";
    const char* CHARM_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/CHARM_analysis/test_v000/before_missing_momentum_theta_cut";

    LetsFill(Knunu_dirname, SIGNAL_hist, "B2Knunu", Scale_Kplus_validation_MC15rd);
    LetsFill(Kstarnunu_dirname, SIGNAL_hist, "B2Kstarnunu", Scale_Kplusstar_validation_MC15rd);
    LetsFill(Xsununu_dirname, SIGNAL_hist, "B2Xsnunu", Scale_Xsu_nonresonant_validation_MC15rd);
    LetsFill(K0nunu_dirname, SIGNAL_hist, "B02K0nunu", Scale_K0_validation_MC15rd);
    LetsFill(K0starnunu_dirname, SIGNAL_hist, "B02Kstar0nunu", Scale_K0star_validation_MC15rd);
    LetsFill(Xsdnunu_dirname, SIGNAL_hist, "B02Xsnunu", Scale_Xsd_nonresonant_validation_MC15rd);
    LetsFill(CHG_dirname, CHG_hist, Scale_CHG_validation_MC15rd);
    LetsFill(MIX_dirname, MIX_hist, Scale_MIX_validation_MC15rd);
    LetsFill(UUBAR_dirname, UUBAR_hist, Scale_UUBAR_validation_MC15rd);
    LetsFill(DDBAR_dirname, DDBAR_hist, Scale_DDBAR_validation_MC15rd);
    LetsFill(SSBAR_dirname, SSBAR_hist, Scale_SSBAR_validation_MC15rd);
    LetsFill(CHARM_dirname, CHARM_hist, Scale_CHARM_validation_MC15rd);


    double CHG_int = CHG_hist->Integral();
    double MIX_int = MIX_hist->Integral();
    double UUBAR_int = UUBAR_hist->Integral();
    double DDBAR_int = DDBAR_hist->Integral();
    double SSBAR_int = SSBAR_hist->Integral();
    double CHARM_int = CHARM_hist->Integral();
    double SIGNAL_int = SIGNAL_hist->Integral();

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
    Stack->SetMaximum(Stack->GetMaximum() * 1.3);

    //SIGNAL_hist->Scale(10000.0/BKG_int, "width");
    //SIGNAL_hist->Scale(1.0 / SIGNAL_int, "width");
    SIGNAL_hist->Scale(10000.0);
    SIGNAL_hist->SetLineWidth(3);
    SIGNAL_hist->SetLineColor(2);
    SIGNAL_hist->SetFillStyle(0);

    TCanvas* c_temp = new TCanvas("c", "", 1500, 1200); c_temp->cd();

    gStyle->SetPalette(kGistEarth);

    Stack->Draw("pfc Hist"); Stack->GetXaxis()->SetTitleSize(0.042); Stack->GetYaxis()->SetTitleSize(0.042); Stack->GetYaxis()->SetMaxDigits(2);
    SIGNAL_hist->Draw("HistSAME");

    TPaveText* pt_lumi = new TPaveText(0.12, 0.83, 0.42, 0.85, "NDC NB");
    pt_lumi->SetTextSize(0.035); pt_lumi->SetFillStyle(0); pt_lumi->SetLineWidth(0); pt_lumi->SetTextAlign(11); pt_lumi->AddText("Belle II simulation #int L dt = 365.4 fb^{-1}"); pt_lumi->Draw();

    TLegend* legend = gPad->BuildLegend(0.95, 0.9, 0.68, 0.65);
    legend->SetFillStyle(0); legend->SetLineWidth(0);

    Float_t ymax = Stack->GetMaximum();
    TLine* lineA = new TLine(0.297, 0, 0.297, ymax);
    lineA->SetLineColor(kRed);
    lineA->SetLineStyle(5);
    lineA->SetLineWidth(3);
    lineA->Draw();
    TLine* lineB = new TLine(2.618, 0, 2.618, ymax);
    lineB->SetLineColor(kRed);
    lineB->SetLineStyle(5);
    lineB->SetLineWidth(3);
    lineB->Draw();
    c_temp->SaveAs("Plot_missing_theta.png");

//    free(line);
//    gPad->BuildLegend();
//    c_temp->SaveAs("Plot_with_legend.png");
}
