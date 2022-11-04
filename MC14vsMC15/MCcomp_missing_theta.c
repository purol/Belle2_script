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

// my MC sample number for MC14
# define N_Kplus_nunubar_MC14 10000000.0
# define N_K0_nunubar_MC14 10000000.0
# define N_Kplusstar_nunubar_MC14 10000000.0
# define N_K0star_nunubar_MC14 10000000.0
# define N_Xsu_nonresonant_nunubar_MC14 50000000.0
# define N_Xsd_nonresonant_nunubar_MC14 50000000.0

// scale factor for MC14 (364.436/fb)
# define Scale_Kplus_MC14 (0.364436 * N_Kplus_nunubar_1invab/N_Kplus_nunubar_MC14)
# define Scale_Kplusstar_MC14 (0.364436 * N_Kplusstar_nunubar_1invab/N_Kplusstar_nunubar_MC14)
# define Scale_Xsu_nonresonant_MC14 (0.364436 * N_Xsu_nonresonant_nunubar_1invab/N_Xsu_nonresonant_nunubar_MC14)
# define Scale_K0_MC14 (0.364436 * N_K0_nunubar_1invab/N_K0_nunubar_MC14)
# define Scale_K0star_MC14 (0.364436 * N_K0star_nunubar_1invab/N_K0star_nunubar_MC14)
# define Scale_Xsd_nonresonant_MC14 (0.364436 * N_Xsd_nunubar_1invab/N_Xsd_nonresonant_nunubar_MC14)

// SIGNAL MC sample number for MC15
# define N_Kplus_test_MC15 2961000.0
# define N_K0_test_MC15 2833376.0
# define N_Kplusstar_test_MC15 2961000.0
# define N_K0star_test_MC15 2833376.0
# define N_Xsu_nonresonant_test_MC15 14805000.0
# define N_Xsd_nonresonant_test_MC15 15059570.0

// scale factor for SIGNAL MC15 sample (364.436/fb)
# define Scale_Kplus_test_MC15 (0.364436 * N_Kplus_nunubar_1invab/N_Kplus_test_MC15)
# define Scale_Kplusstar_test_MC15 (0.364436 * N_Kplusstar_nunubar_1invab/N_Kplusstar_test_MC15)
# define Scale_Xsu_nonresonant_test_MC15 (0.364436 * N_Xsu_nonresonant_nunubar_1invab/N_Xsu_nonresonant_test_MC15)
# define Scale_K0_test_MC15 (0.364436 * N_K0_nunubar_1invab/N_K0_test_MC15)
# define Scale_K0star_test_MC15 (0.364436 * N_K0star_nunubar_1invab/N_K0star_test_MC15)
# define Scale_Xsd_nonresonant_test_MC15 (0.364436 * N_Xsd_nunubar_1invab/N_Xsd_nonresonant_test_MC15)

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

TH1D* SIGNAL_hist_MC14 = new TH1D("MC14ria", ";#theta_{missing};number of candidates", 100, 0, 3.142);
TH1D* SIGNAL_hist_MC15 = new TH1D("MC15rib", ";#theta_{missing};number of candidates", 100, 0, 3.142);
TH1D* Ratio_hist = new TH1D("MC15rib/MC14ria", ";;MC15rib/MC14ria", 100, 0, 3.142);

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

void MCcomp_missing_theta() {

    const char* MC15_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v000/before_missing_momentum_theta_cut";
    const char* MC14_dirname = "/home/jwpark/storage/BKG_gbasf2/Kasen_MC14/SIGNAL_analysis/validation_v000/before_missing_momentum_theta_cut";

    LetsFill(MC14_dirname, SIGNAL_hist_MC14, "B2Knunu", Scale_Kplus_MC14);
    LetsFill(MC14_dirname, SIGNAL_hist_MC14, "B2Kstarnunu", Scale_Kplusstar_MC14);
    LetsFill(MC14_dirname, SIGNAL_hist_MC14, "B2Xsnunu", Scale_Xsu_nonresonant_MC14);
    LetsFill(MC14_dirname, SIGNAL_hist_MC14, "B02K0nunu", Scale_K0_MC14);
    LetsFill(MC14_dirname, SIGNAL_hist_MC14, "B02Kstar0nunu", Scale_K0star_MC14);
    LetsFill(MC14_dirname, SIGNAL_hist_MC14, "B02Xsnunu", Scale_Xsd_nonresonant_MC14);

    LetsFill(MC15_dirname, SIGNAL_hist_MC15, "B2Knunu", Scale_Kplus_test_MC15);
    LetsFill(MC15_dirname, SIGNAL_hist_MC15, "B2Kstarnunu", Scale_Kplusstar_test_MC15);
    LetsFill(MC15_dirname, SIGNAL_hist_MC15, "B2Xsnunu", Scale_Xsu_nonresonant_test_MC15);
    LetsFill(MC15_dirname, SIGNAL_hist_MC15, "B02K0nunu", Scale_K0_test_MC15);
    LetsFill(MC15_dirname, SIGNAL_hist_MC15, "B02Kstar0nunu", Scale_K0star_test_MC15);
    LetsFill(MC15_dirname, SIGNAL_hist_MC15, "B02Xsnunu", Scale_Xsd_nonresonant_test_MC15);

    SIGNAL_hist_MC15->SetFillStyle(3004);
    SIGNAL_hist_MC15->SetLineColor(kBlue);
    SIGNAL_hist_MC15->SetFillColor(kBlue);

    SIGNAL_hist_MC14->SetFillStyle(3005);
    SIGNAL_hist_MC14->SetLineColor(kRed);
    SIGNAL_hist_MC14->SetFillColor(kRed);

    Ratio_hist->SetLineColor(kBlack); Ratio_hist->SetMarkerStyle(21); Ratio_hist->Sumw2(); Ratio_hist->SetStats(0);
    Ratio_hist->Divide(SIGNAL_hist_MC15, SIGNAL_hist_MC14);

    gStyle->SetOptStat(0);

    TCanvas* c_temp = new TCanvas("c", "", 800, 800); c_temp->cd();

    TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.35, 1.0, 1.0);
    pad1->SetBottomMargin(0.08); pad1->SetLeftMargin(0.15);
    pad1->SetGridx(); pad1->Draw(); pad1->cd();

    gStyle->SetPalette(kPastel);

    Float_t ymax_1 = SIGNAL_hist_MC15->GetMaximum();
    Float_t ymax_2 = SIGNAL_hist_MC14->GetMaximum();
    double real_max = 0;
    if (ymax_1 > ymax_2) real_max = ymax_1;
    else real_max = ymax_2;

    SIGNAL_hist_MC14->SetMaximum(real_max * 1.1);

    SIGNAL_hist_MC14->Draw("pfc Hist");
    SIGNAL_hist_MC15->Draw("SAME Hist");
    TLegend* legend = pad1->BuildLegend(0.9, 0.9, 0.7, 0.7);
    legend->SetFillStyle(0); legend->SetLineWidth(0);

    c_temp->cd();
    TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1, 0.3); pad2->SetBottomMargin(0.15); pad2->SetLeftMargin(0.15); pad2->SetGridx(); pad2->Draw(); pad2->cd();
    Ratio_hist->SetMinimum(0.5); Ratio_hist->SetMaximum(1.5); Ratio_hist->SetLineWidth(2);
    Ratio_hist->GetYaxis()->SetTitleSize(0.08); Ratio_hist->GetYaxis()->SetTitleOffset(0.5);
    Ratio_hist->GetXaxis()->SetLabelSize(0.08); Ratio_hist->GetYaxis()->SetLabelSize(0.08);
    Ratio_hist->Draw("e0p");
    TLine* line = new TLine(Ratio_hist->GetXaxis()->GetXmin(), 1.0, Ratio_hist->GetXaxis()->GetXmax(), 1.0);
    line->SetLineColor(kRed);
    line->SetLineStyle(1); line->SetLineWidth(3);
    line->Draw();

    c_temp->SetBottomMargin(0.0);
    c_temp->SaveAs("missing_MC14vsMC15.png");

    delete c_temp;
}
