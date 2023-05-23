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

/*
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
*/

// SIGNAL MC sample number
# define N_Kplus_train 7039000.0
# define N_K0_train 7166624.0
# define N_Kplusstar_train 7039000.0
# define N_K0star_train 7166624.0
# define N_Xsu_nonresonant_train 35195000.0
# define N_Xsd_nonresonant_train 34940430.0
# define N_Kplus_test 2961000.0
# define N_K0_test 2833376.0
# define N_Kplusstar_test 2961000.0
# define N_K0star_test 2833376.0
# define N_Xsu_nonresonant_test 14805000.0
# define N_Xsd_nonresonant_test 15059570.0

// scale factor for SIGNAL MC sample (364.436/fb)
# define Scale_Kplus_train (0.364436 * N_Kplus_nunubar_1invab/N_Kplus_train)
# define Scale_Kplusstar_train (0.364436 * N_Kplusstar_nunubar_1invab/N_Kplusstar_train)
# define Scale_Xsu_nonresonant_train (0.364436 * N_Xsu_nonresonant_nunubar_1invab/N_Xsu_nonresonant_train)
# define Scale_K0_train (0.364436 * N_K0_nunubar_1invab/N_K0_train)
# define Scale_K0star_train (0.364436 * N_K0star_nunubar_1invab/N_K0star_train)
# define Scale_Xsd_nonresonant_train (0.364436 * N_Xsd_nunubar_1invab/N_Xsd_nonresonant_train)
# define Scale_Kplus_test (0.364436 * N_Kplus_nunubar_1invab/N_Kplus_test)
# define Scale_Kplusstar_test (0.364436 * N_Kplusstar_nunubar_1invab/N_Kplusstar_test)
# define Scale_Xsu_nonresonant_test (0.364436 * N_Xsu_nonresonant_nunubar_1invab/N_Xsu_nonresonant_test)
# define Scale_K0_test (0.364436 * N_K0_nunubar_1invab/N_K0_test)
# define Scale_K0star_test (0.364436 * N_K0star_nunubar_1invab/N_K0star_test)
# define Scale_Xsd_nonresonant_test (0.364436 * N_Xsd_nunubar_1invab/N_Xsd_nonresonant_test)

// BKG MC sample number (2.8/ab for BB, 1.0/ab for qq)
# define N_CHG_test 139768443.0
# define N_MIX_test 104591690.0
# define N_UUBAR_test 107641351.0
# define N_DDBAR_test 25588403.0
# define N_SSBAR_test 20668794.0
# define N_CHARM_test 116108850.0
# define N_CHG_train 151915459.0
# define N_MIX_train 120136353.0
# define N_UUBAR_train 158509639.0
# define N_DDBAR_train 38644413.0
# define N_SSBAR_train 31001866.0
# define N_CHARM_train 174901296.0

// new scale factor for BKG MC sample with additional 1/ab (364.436 - 2.763 = 361.673/fb), until LS1
# define Scale_CHG_train ((N_BB_LS1* (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab * (N_CHG_train / (N_CHG_train + N_CHG_test)) ))
# define Scale_MIX_train ((N_BB_LS1* (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab * (N_MIX_train / (N_MIX_train + N_MIX_test)) ))
# define Scale_UUBAR_train (0.361673/((N_UUBAR_train/(N_UUBAR_train + N_UUBAR_test))*1.0))
# define Scale_DDBAR_train (0.361673/((N_DDBAR_train/(N_DDBAR_train + N_DDBAR_test))*1.0))
# define Scale_SSBAR_train (0.361673/((N_SSBAR_train/(N_SSBAR_train + N_SSBAR_test))*1.0))
# define Scale_CHARM_train (0.361673/((N_CHARM_train/(N_CHARM_train + N_CHARM_test))*1.0))
# define Scale_CHG_test ((N_BB_LS1* (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab * (N_CHG_test / (N_CHG_train + N_CHG_test)) ))
# define Scale_MIX_test ((N_BB_LS1* (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab * (N_MIX_test / (N_MIX_train + N_MIX_test)) ))
# define Scale_UUBAR_test (0.361673/((N_UUBAR_test/(N_UUBAR_train + N_UUBAR_test))*1.0))
# define Scale_DDBAR_test (0.361673/((N_DDBAR_test/(N_DDBAR_train + N_DDBAR_test))*1.0))
# define Scale_SSBAR_test (0.361673/((N_SSBAR_test/(N_SSBAR_train + N_SSBAR_test))*1.0))
# define Scale_CHARM_test (0.361673/((N_CHARM_test/(N_CHARM_train + N_CHARM_test))*1.0))

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

    const char* Knunu_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/SIGNAL_analysis/validation_v000/before_missing_momentum_theta_cut";
    const char* Kstarnunu_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/SIGNAL_analysis/validation_v000/before_missing_momentum_theta_cut";
    const char* Xsununu_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/SIGNAL_analysis/validation_v000/before_missing_momentum_theta_cut";
    const char* K0nunu_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/SIGNAL_analysis/validation_v000/before_missing_momentum_theta_cut";
    const char* K0starnunu_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/SIGNAL_analysis/validation_v000/before_missing_momentum_theta_cut";
    const char* Xsdnunu_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/SIGNAL_analysis/validation_v000/before_missing_momentum_theta_cut";
    const char* CHG_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/CHG_analysis/validation_v000/before_missing_momentum_theta_cut";
    const char* MIX_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/MIX_analysis/validation_v000/before_missing_momentum_theta_cut";
    const char* UUBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/UUBAR_analysis/validation_v000/before_before_missing_momentum_theta_cut";
    const char* DDBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/DDBAR_analysis/validation_v000/before_missing_momentum_theta_cut";
    const char* SSBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/SSBAR_analysis/validation_v000/before_missing_momentum_theta_cut";
    const char* CHARM_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/CHARM_analysis/validation_v000/before_missing_momentum_theta_cut";

    LetsFill(Knunu_dirname, SIGNAL_hist, "B2Knunu", Scale_Kplus_test);
    LetsFill(Kstarnunu_dirname, SIGNAL_hist, "B2Kstarnunu", Scale_Kplusstar_test);
    LetsFill(Xsununu_dirname, SIGNAL_hist, "B2Xsnunu", Scale_Xsu_nonresonant_test);
    LetsFill(K0nunu_dirname, SIGNAL_hist, "B02K0nunu", Scale_K0_test);
    LetsFill(K0starnunu_dirname, SIGNAL_hist, "B02Kstar0nunu", Scale_K0star_test);
    LetsFill(Xsdnunu_dirname, SIGNAL_hist, "B02Xsnunu", Scale_Xsd_nonresonant_test);
    LetsFill(CHG_dirname, CHG_hist, Scale_CHG_test);
    LetsFill(MIX_dirname, MIX_hist, Scale_MIX_test);
    LetsFill(UUBAR_dirname, UUBAR_hist, Scale_UUBAR_test);
    LetsFill(DDBAR_dirname, DDBAR_hist, Scale_DDBAR_test);
    LetsFill(SSBAR_dirname, SSBAR_hist, Scale_SSBAR_test);
    LetsFill(CHARM_dirname, CHARM_hist, Scale_CHARM_test);


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
    //Stack->SetMaximum(100.0);

    //SIGNAL_hist->Scale(10000.0/BKG_int, "width");
    //SIGNAL_hist->Scale(1.0 / SIGNAL_int, "width");
    SIGNAL_hist->Scale(10000.0);
    SIGNAL_hist->SetLineWidth(3);
    SIGNAL_hist->SetLineColor(2);
    SIGNAL_hist->SetFillStyle(0);

    TCanvas* c_temp = new TCanvas("c", "", 1500, 1200); c_temp->cd();

    gStyle->SetPalette(kGistEarth);

    Stack->Draw("pfc Hist"); SIGNAL_hist->Draw("HistSAME");

    TLegend* legend = gPad->BuildLegend(0.9, 0.9, 0.7, 0.7);
    legend->SetFillStyle(0);

    Float_t ymax = Stack->GetMaximum();
    TLine* lineA = new TLine(0.297, 0, 0.297, ymax);
    lineA->SetLineColor(kRed);
    lineA->SetLineStyle(5);
    lineA->Draw();
    TLine* lineB = new TLine(2.618, 0, 2.618, ymax);
    lineB->SetLineColor(kRed);
    lineB->SetLineStyle(5);
    lineB->Draw();
    c_temp->SaveAs("Plot_missing_theta.png");

//    free(line);
//    gPad->BuildLegend();
//    c_temp->SaveAs("Plot_with_legend.png");
}
