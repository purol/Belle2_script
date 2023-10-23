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

TH2F* SIGNAL_hist = new TH2F("SIGNAL_hist", ";MBsigop0 [GeV^{2}];M_{bc} [GeV]", 15, -5, 35, 15, 5.27, 5.29);
TH2F* BKG_hist = new TH2F("BKG_hist", ";MBsigop0 [GeV^{2}];M_{bc} [GeV]", 15, -5, 35, 15, 5.27, 5.29);

void LetsFill(const char* dirname, TH2F* hist, double weight) {
    double var1 = 0;
    double var2 = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("MsquaredBsig_op0", &var1);
        tree_Btag->SetBranchAddress("Btag_Mbc", &var2);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            hist->Fill(var1,var2,weight);
        }
        input_file->Close();

    }

}

void twoDplot() {

    const char* B2Knunu_dirname = "./SIGNAL/B2Knunu";
    const char* B2Kstarnunu_dirname = "./SIGNAL/B2Kstarnunu";
    const char* B2Xsnunu_dirname = "./SIGNAL/B2Xsnunu";
    const char* B02K0nunu_dirname = "./SIGNAL/B02K0nunu";
    const char* B02K0starnunu_dirname = "./SIGNAL/B02K0starnunu";
    const char* B02Xsnunu_dirname = "./SIGNAL/B02Xsnunu";
    const char* BKG_dirname = "./BKG";

    LetsFill(B2Knunu_dirname, SIGNAL_hist, Scale_Kplus);
    LetsFill(B2Kstarnunu_dirname, SIGNAL_hist, Scale_Kplusstar);
    LetsFill(B2Xsnunu_dirname, SIGNAL_hist, Scale_Xsu_nonresonant);
    LetsFill(B02K0nunu_dirname, SIGNAL_hist, Scale_K0);
    LetsFill(B02K0starnunu_dirname, SIGNAL_hist, Scale_K0star);
    LetsFill(B02Xsnunu_dirname, SIGNAL_hist, Scale_Xsd_nonresonant);
    LetsFill(BKG_dirname, BKG_hist, 1);

    TCanvas* c_temp = new TCanvas("c", "", 1500, 1200); c_temp->cd();

    SIGNAL_hist->Draw("colz"); c_temp->SaveAs("SIGNAL.png");
    BKG_hist->Draw("colz"); c_temp->SaveAs("BKG.png");
}
