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

// arXiv:1409.4557v2, PhysRevD.107.014511
# define TB0 1.5195 // (Table. 1)
# define TBp 1.6384 // (Table. 1)
# define BR_Kplus_nunubar 0.000005044 // Table VI = (5.044 +- 0.402) * 10^{-6}
# define BR_K0star_nunubar 0.00000919 // (eq. 11)
# define BR_K0_nunubar (BR_Kplus_nunubar*TB0/TBp) // under (eq. 15). In Table VI = (4.6669 +- 0.3707) * 10^{-6}
# define BR_Kplusstar_nunubar (BR_K0star_nunubar*TBp/TB0) // under (eq. 15)
# define BR_Xs_nunubar 0.000029 // (eq. 23)
# define BR_Xsu_nonresonant_nunubar (BR_Xs_nunubar - BR_Kplus_nunubar - BR_Kplusstar_nunubar)
# define BR_Xsd_nonresonant_nunubar (BR_Xs_nunubar - BR_K0_nunubar - BR_K0star_nunubar)

// according to DIRAC
# define N_BpBp_1invab 540000000.0
# define N_B0B0_1invab 510000000.0

# define BR_BpBp 0.514
# define BR_B0B0 0.486

// https://confluence.desy.de/pages/viewpage.action?spaceKey=BI&title=Conference+readiness
# define N_BB_LS1 387100000.0 // NBB = (387.1 +/- 5.6) x 10^6

# define N_Kplus_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_BpBp/(BR_BpBp+BR_B0B0)) * BR_Kplus_nunubar)
# define N_Kplusstar_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_BpBp/(BR_BpBp+BR_B0B0)) * BR_Kplusstar_nunubar)
# define N_Xsu_nonresonant_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_BpBp/(BR_BpBp+BR_B0B0)) * BR_Xsu_nonresonant_nunubar)
# define N_K0_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_B0B0/(BR_BpBp+BR_B0B0)) * BR_K0_nunubar)
# define N_K0star_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_B0B0/(BR_BpBp+BR_B0B0)) * BR_K0star_nunubar)
# define N_Xsd_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_B0B0/(BR_BpBp+BR_B0B0)) * BR_Xsd_nonresonant_nunubar)

// SIGNAL MC sample number before skimming
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

// SIGNAL MC sample number before skimming for MC15rd
# define N_Kplus_validation_MC15rd 3962022.0
# define N_K0_validation_MC15rd 3962022.0
# define N_Kplusstar_validation_MC15rd 3962022.0
# define N_K0star_validation_MC15rd 3962022.0
# define N_Xsu_nonresonant_validation_MC15rd 15846594.0
# define N_Xsd_nonresonant_validation_MC15rd 15846594.0

// scale factor for SIGNAL MC sample until LS1
# define Scale_Kplus_train (N_Kplus_nunubar_LS1/N_Kplus_train)
# define Scale_Kplusstar_train (N_Kplusstar_nunubar_LS1/N_Kplusstar_train)
# define Scale_Xsu_nonresonant_train (N_Xsu_nonresonant_nunubar_LS1/N_Xsu_nonresonant_train)
# define Scale_K0_train (N_K0_nunubar_LS1/N_K0_train)
# define Scale_K0star_train (N_K0star_nunubar_LS1/N_K0star_train)
# define Scale_Xsd_nonresonant_train (N_Xsd_nunubar_LS1/N_Xsd_nonresonant_train)
# define Scale_Kplus_test (N_Kplus_nunubar_LS1/N_Kplus_test)
# define Scale_Kplusstar_test (N_Kplusstar_nunubar_LS1/N_Kplusstar_test)
# define Scale_Xsu_nonresonant_test (N_Xsu_nonresonant_nunubar_LS1/N_Xsu_nonresonant_test)
# define Scale_K0_test (N_K0_nunubar_LS1/N_K0_test)
# define Scale_K0star_test (N_K0star_nunubar_LS1/N_K0star_test)
# define Scale_Xsd_nonresonant_test (N_Xsd_nunubar_LS1/N_Xsd_nonresonant_test)

// scale factor for SIGNAL MC sample until LS1 for MC15rd
# define Scale_Kplus_validation_MC15rd (N_Kplus_nunubar_LS1/N_Kplus_validation_MC15rd)
# define Scale_Kplusstar_validation_MC15rd (N_Kplusstar_nunubar_LS1/N_Kplusstar_validation_MC15rd)
# define Scale_Xsu_nonresonant_validation_MC15rd (N_Xsu_nonresonant_nunubar_LS1/N_Xsu_nonresonant_validation_MC15rd)
# define Scale_K0_validation_MC15rd (N_K0_nunubar_LS1/N_K0_validation_MC15rd)
# define Scale_K0star_validation_MC15rd (N_K0star_nunubar_LS1/N_K0star_validation_MC15rd)
# define Scale_Xsd_nonresonant_validation_MC15rd (N_Xsd_nunubar_LS1/N_Xsd_nonresonant_validation_MC15rd)

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
# define Scale_CHG_validation ((N_BB_LS1* (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab))
# define Scale_MIX_validation ((N_BB_LS1* (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab))
# define Scale_UUBAR_validation (0.361673)
# define Scale_DDBAR_validation (0.361673)
# define Scale_SSBAR_validation (0.361673)
# define Scale_CHARM_validation (0.361673)

// BKG MC sample number for MC15rd
# define N_CHG_validation_MC15rd 785108449.0 // 1458.959/fb
# define N_MIX_validation_MC15rd 741492304.0 // 1458.959/fb
# define N_UUBAR_validation_MC15rd 2306265848.0 // 1458.959/fb
# define N_DDBAR_validation_MC15rd 576209482.0 // 1458.959/fb
# define N_SSBAR_validation_MC15rd 526874294.0 // 1458.959/fb
# define N_CHARM_validation_MC15rd 1889822323.0 // 1458.959/fb

# define Scale_CHG_validation_MC15rd (0.361673/1.458959)
# define Scale_MIX_validation_MC15rd (0.361673/1.458959)
# define Scale_UUBAR_validation_MC15rd (0.361673/1.458959)
# define Scale_DDBAR_validation_MC15rd (0.361673/1.458959)
# define Scale_SSBAR_validation_MC15rd (0.361673/1.458959)
# define Scale_CHARM_validation_MC15rd (0.361673/1.458959)

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

TH1D* SIGNAL_hist = new TH1D("signal", ";##Upsilon(4S) multiplicity;Fraction of events", 6, -0.5, 5.5);
TH1D* CHG_hist = new TH1D("charged", ";#Upsilon(4S) multiplicity;Fraction of events", 6, -0.5, 5.5);
TH1D* MIX_hist = new TH1D("mixed", ";#Upsilon(4S) multiplicity;Fraction of events", 6, -0.5, 5.5);
TH1D* UUBAR_hist = new TH1D("u#bar{u}", ";#Upsilon(4S) multiplicity;Fraction of events", 6, -0.5, 5.5);
TH1D* DDBAR_hist = new TH1D("d#bar{d}", ";#Upsilon(4S) multiplicity;Fraction of events", 6, -0.5, 5.5);
TH1D* SSBAR_hist = new TH1D("s#bar{s}", ";#Upsilon(4S) multiplicity;Fraction of events", 6, -0.5, 5.5);
TH1D* CHARM_hist = new TH1D("c#bar{c}", ";#Upsilon(4S) multiplicity;Fraction of events", 6, -0.5, 5.5);

void LetsFill(const char* dirname, TH1D* hist, double weight = 1) {

    double MXs = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        int __experiment__ = -1;
        int __run__ = -1;
        unsigned int __event__ = -1;
        int __ncandidates__ = -1;

        int PREV__experiment__ = -1;
        int PREV__run__ = -1;
        unsigned int PREV__event__ = -1;
        int PREV__ncandidates__ = -1;

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_Bsig->SetBranchAddress("Bsig_M", &MXs);
        tree_upsilon->SetBranchAddress("__experiment__", &__experiment__);
        tree_upsilon->SetBranchAddress("__run__", &__run__);
        tree_upsilon->SetBranchAddress("__event__", &__event__);
        tree_upsilon->SetBranchAddress("__ncandidates__", &__ncandidates__);

        int Nevt = 0;

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            // D veto
            if ((MXs > 1.84) && (MXs < 1.89)) continue;

            if ((-1 == PREV__experiment__) && (-1 == PREV__run__) && (-1 == PREV__event__) && (-1 == PREV__ncandidates__)) {
                Nevt++;
            }
            else if ((__experiment__ == PREV__experiment__) && (__run__ == PREV__run__) && (__event__ == PREV__event__) && (__ncandidates__ == PREV__ncandidates__)) {
                Nevt++;
            }
            else {
                hist->Fill(Nevt, weight);
                Nevt = 0;
            }

            PREV__experiment__ = __experiment__;
            PREV__run__ = __run__;
            PREV__event__ = __event__;
            PREV__ncandidates__ = __ncandidates__;

        }
        input_file->Close();

    }

}

void LetsFill(const char* dirname, TH1D* hist, const char* included_string, double weight = 1) {

    double MXs = 0;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    for (unsigned int i = 0; i < names.size(); i++) {

        int __experiment__ = -1;
        int __run__ = -1;
        unsigned int __event__ = -1;
        int __ncandidates__ = -1;

        int PREV__experiment__ = -1;
        int PREV__run__ = -1;
        unsigned int PREV__event__ = -1;
        int PREV__ncandidates__ = -1;

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_Bsig->SetBranchAddress("Bsig_M", &MXs);
        tree_upsilon->SetBranchAddress("__experiment__", &__experiment__);
        tree_upsilon->SetBranchAddress("__run__", &__run__);
        tree_upsilon->SetBranchAddress("__event__", &__event__);
        tree_upsilon->SetBranchAddress("__ncandidates__", &__ncandidates__);

        int Nevt = 0;

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            // D veto
            if ((MXs > 1.84) && (MXs < 1.89)) continue;

            if ((-1 == PREV__experiment__) && (-1 == PREV__run__) && (-1 == PREV__event__) && (-1 == PREV__ncandidates__)) {
                Nevt++;
            }
            else if ((__experiment__ == PREV__experiment__) && (__run__ == PREV__run__) && (__event__ == PREV__event__) && (__ncandidates__ == PREV__ncandidates__)) {
                Nevt++;
            }
            else {
                hist->Fill(Nevt, weight);
                Nevt = 0;
            }

            PREV__experiment__ = __experiment__;
            PREV__run__ = __run__;
            PREV__event__ = __event__;
            PREV__ncandidates__ = __ncandidates__;

        }
        input_file->Close();

    }

}

void THStack_multiplicity() {

    const char* Knunu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/validation_v004/before_Dveto_cut";
    const char* Kstarnunu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/validation_v004/before_Dveto_cut";
    const char* Xsununu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/validation_v004/before_Dveto_cut";
    const char* K0nunu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/validation_v004/before_Dveto_cut";
    const char* K0starnunu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/validation_v004/before_Dveto_cut";
    const char* Xsdnunu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/validation_v004/before_Dveto_cut";
    const char* CHG_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/SatoriRD/CHG_analysis/validation_v004/before_Dveto_cut/Merge";
    const char* MIX_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/MIX_analysis/validation_v004/before_Dveto_cut";
    const char* UUBAR_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/UUBAR_analysis/validation_v004/before_Dveto_cut";
    const char* DDBAR_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/DDBAR_analysis/validation_v004/before_Dveto_cut";
    const char* SSBAR_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SSBAR_analysis/validation_v004/before_Dveto_cut";
    const char* CHARM_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/CHARM_analysis/validation_v004/before_Dveto_cut";

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

    
    CHG_hist->Scale(1.0 / CHG_int, "width");
    MIX_hist->Scale(1.0 / MIX_int, "width");
    UUBAR_hist->Scale(1.0 / UUBAR_int, "width");
    DDBAR_hist->Scale(1.0 / DDBAR_int, "width");
    SSBAR_hist->Scale(1.0 / SSBAR_int, "width");
    CHARM_hist->Scale(1.0 / CHARM_int, "width");
    SIGNAL_hist->Scale(1.0 / SIGNAL_int, "width");

    TGaxis::SetMaxDigits(3);

    CHG_hist->SetLineWidth(3);
    CHG_hist->SetLineColor(6);
    CHG_hist->SetFillStyle(0);

    MIX_hist->SetLineWidth(3);
    MIX_hist->SetLineColor(3);
    MIX_hist->SetFillStyle(0);

    UUBAR_hist->SetLineWidth(3);
    UUBAR_hist->SetLineColor(4);
    UUBAR_hist->SetFillStyle(0);

    DDBAR_hist->SetLineWidth(3);
    DDBAR_hist->SetLineColor(1);
    DDBAR_hist->SetFillStyle(0);

    SSBAR_hist->SetLineWidth(3);
    SSBAR_hist->SetLineColor(7);
    SSBAR_hist->SetFillStyle(0);

    CHARM_hist->SetLineWidth(3);
    CHARM_hist->SetLineColor(8);
    CHARM_hist->SetFillStyle(0);

    SIGNAL_hist->SetLineWidth(3);
    SIGNAL_hist->SetLineColor(2);
    SIGNAL_hist->SetFillStyle(0);

    TCanvas* c_temp = new TCanvas("c", "", 1500, 1200); c_temp->cd();

    gStyle->SetPalette(kGistEarth);

    CHG_hist->SetStats(false);
    MIX_hist->SetStats(false);
    UUBAR_hist->SetStats(false);
    DDBAR_hist->SetStats(false);
    SSBAR_hist->SetStats(false);
    CHARM_hist->SetStats(false);
    SIGNAL_hist->SetStats(false);

    CHG_hist->Draw("Hist");
    MIX_hist->Draw("HistSAME");
    UUBAR_hist->Draw("HistSAME");
    DDBAR_hist->Draw("HistSAME");
    SSBAR_hist->Draw("HistSAME");
    CHARM_hist->Draw("HistSAME");
    SIGNAL_hist->Draw("HistSAME");

    TLegend* legend = gPad->BuildLegend(0.9, 0.9, 0.7, 0.7);
    legend->SetFillStyle(0); legend->SetLineWidth(0);

    c_temp->SaveAs("Plot_multiplicity.png");

//    free(line);
//    gPad->BuildLegend();
//    c_temp->SaveAs("Plot_with_legend.png");
}
