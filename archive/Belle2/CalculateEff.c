// last update: 2022-07-22
// for Belle2 data

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

TH1D* Nevt_initial = new TH1D("q2_initial", ";q^{2} [GeV^{2}]; evt", 3, 0, 25.1);
TH1D* Nevt_final = new TH1D("q2_final", ";q^{2} [GeV^{2}]; evt", 3, 0, 25.1);

TH1D* Ratio_hist new TH1D("eff", ";q^{2} [GeV^{2}]; efficiency", 3, 0, 25.1);

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

void LetsFillInitial(const char* dirname, TH1D* hist, double weight = 1) {
    double var = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_info = (TTree*)input_file->Get("info");

        tree_info->SetBranchAddress("invMassInLists__bonu_e__clPrimaryMC_signal__bc", &var);

        printf("%lld entries...\n", tree_info->GetEntries());
        for (unsigned int j = 0; j < tree_info->GetEntries(); j++) { // Fill
            tree_info->GetEntry(j);
            hist->Fill(var * var, weight);
        }
        input_file->Close();

    }

}

void LetsFillFinal(const char* dirname, TH1D* hist, double weight = 1) {
    double var = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_Xs = (TTree*)input_file->Get("Xs");

        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &var);

        printf("%lld entries...\n", tree_Xs->GetEntries());
        for (unsigned int j = 0; j < tree_Xs->GetEntries(); j++) { // Fill
            tree_Xs->GetEntry(j);
            hist->Fill(var * var, weight);
        }
        input_file->Close();

    }

}

void CalculateEff(){

    const char* Knunu_dirname_init = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/CHARM_analysis/validation_v000/before_Eecl_cut";
    const char* Kstarnunu_dirname_init = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/CHARM_analysis/validation_v000/before_Eecl_cut";
    const char* Xsununu_dirname_init = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/CHARM_analysis/validation_v000/before_Eecl_cut";
    const char* K0nunu_dirname_init = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/CHARM_analysis/validation_v000/before_Eecl_cut";
    const char* K0starnunu_dirname_init = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/CHARM_analysis/validation_v000/before_Eecl_cut";
    const char* Xsdnunu_dirname_init = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/CHARM_analysis/validation_v000/before_Eecl_cut";

    const char* Knunu_dirname_final = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/CHARM_analysis/validation_v000/before_Eecl_cut";
    const char* Kstarnunu_dirname_final = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/CHARM_analysis/validation_v000/before_Eecl_cut";
    const char* Xsununu_dirname_final = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/CHARM_analysis/validation_v000/before_Eecl_cut";
    const char* K0nunu_dirname_final = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/CHARM_analysis/validation_v000/before_Eecl_cut";
    const char* K0starnunu_dirname_final = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/CHARM_analysis/validation_v000/before_Eecl_cut";
    const char* Xsdnunu_dirname_final = "/home/jwpark/storage/BKG_gbasf2/Nitori_ad/CHARM_analysis/validation_v000/before_Eecl_cut";

    LetsFillInitial(Knunu_dirname_init, Nevt_initial, Scale_Kplus);
    LetsFillInitial(Kstarnunu_dirname_init, Nevt_initial, Scale_Kplusstar);
    LetsFillInitial(Xsununu_dirname_init, Nevt_initial, Scale_Xsu_nonresonant);
    LetsFillInitial(K0nunu_dirname_init, Nevt_initial, Scale_K0);
    LetsFillInitial(K0starnunu_dirname_init, Nevt_initial, Scale_K0star);
    LetsFillInitial(Xsdnunu_dirname_init, Nevt_initial, Scale_Xsd_nonresonant);

    LetsFillInitial(Knunu_dirname_final, Nevt_final, Scale_Kplus);
    LetsFillInitial(Kstarnunu_dirname_final, Nevt_final, Scale_Kplusstar);
    LetsFillInitial(Xsununu_dirname_final, Nevt_final, Scale_Xsu_nonresonant);
    LetsFillInitial(K0nunu_dirname_final, Nevt_final, Scale_K0);
    LetsFillInitial(K0starnunu_dirname_final, Nevt_final, Scale_K0star);
    LetsFillInitial(Xsdnunu_dirname_final, Nevt_final, Scale_Xsd_nonresonant);

    Ratio_hist->SetLineColor(kBlack); Ratio_hist->SetMarkerStyle(21); Ratio_hist->SetStats(0);
    Ratio_hist->Divide(Nevt_final, Nevt_initial);

    TCanvas* c_temp = new TCanvas("c", "", 1200, 1200); c_temp->cd();
    Ratio_hist->Draw("ep");

    TLine* line = new TLine(Ratio_hist[k]->GetXaxis()->GetXmin(), 0.00140, Ratio_hist[k]->GetXaxis()->GetXmax(), 0.00140);
    line->SetLineColor(kRed);
    line->SetLineStyle(1); line->SetLineWidth(3);
    line->Draw();

    c_temp->SaveAs("Efficiency.png");

    delete c_temp;
}
