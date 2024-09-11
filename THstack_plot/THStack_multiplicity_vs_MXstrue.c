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

const int NBins = 3;
double edges[NBins + 1] = { 0.0, 0.6, 1.0, 5.3 };

TH1D* SIGNAL_Nevt = new TH1D("SIGNAL_Nevt", ";M_{X_{s}}^{true} [GeV/c^{2}];the number of events", NBins, edges);
TH1D* SIGNAL_Ncandidates = new TH1D("SIGNAL_Ncandidates", ";M_{X_{s}}^{true} [GeV/c^{2}];the number of candidates", NBins, edges);
TH1D* SIGNAL_multiplicity = new TH1D("SIGNAL_multiplicity", ";M_{X_{s}}^{true} [GeV/c^{2}];multiplicity", NBins, edges);

void LetsFill(const char* dirname, TH1D* hist_Nevt, TH1D* hist_Ncandidate, const char* included_string, bool IsItBplus, double weight = 1) {

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

        double Mxs_Bc_MC = -1;
        double Mxs_B0_MC = -1;
        double MC_MXs = -1;

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");
        TTree* tree_Xs = (TTree*)input_file->Get("Xs");

        tree_Bsig->SetBranchAddress("Bsig_M", &MXs);
        tree_upsilon->SetBranchAddress("__experiment__", &__experiment__);
        tree_upsilon->SetBranchAddress("__run__", &__run__);
        tree_upsilon->SetBranchAddress("__event__", &__event__);
        tree_upsilon->SetBranchAddress("__ncandidates__", &__ncandidates__);

        tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
        tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);

        int Ncandidates = 0;

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            tree_Xs->GetEntry(j);

            if (IsItBplus) MC_MXs = Mxs_Bc_MC;
            else MC_MXs = Mxs_B0_MC;

            // D veto
            if ((MXs > 1.84) && (MXs < 1.89)) continue;

            if ((-1 == PREV__experiment__) && (-1 == PREV__run__) && (-1 == PREV__event__) && (-1 == PREV__ncandidates__)) {
                hist_Nevt->Fill(MC_MXs, weight);
                hist_Ncandidate->Fill(MC_MXs, weight);
            }
            else if ((__experiment__ == PREV__experiment__) && (__run__ == PREV__run__) && (__event__ == PREV__event__) && (__ncandidates__ == PREV__ncandidates__)) {
                hist_Ncandidate->Fill(MC_MXs, weight);
            }
            else {
                hist_Nevt->Fill(MC_MXs, weight);
                hist_Ncandidate->Fill(MC_MXs, weight);
            }

            PREV__experiment__ = __experiment__;
            PREV__run__ = __run__;
            PREV__event__ = __event__;
            PREV__ncandidates__ = __ncandidates__;

        }
        input_file->Close();

    }

}

void THStack_multiplicity_vs_MXstrue() {

    const char* Knunu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/validation_v004/before_Dveto_cut";
    const char* Kstarnunu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/validation_v004/before_Dveto_cut";
    const char* Xsununu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/validation_v004/before_Dveto_cut";
    const char* K0nunu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/validation_v004/before_Dveto_cut";
    const char* K0starnunu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/validation_v004/before_Dveto_cut";
    const char* Xsdnunu_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/SatoriRD/SIGNAL_analysis/validation_v004/before_Dveto_cut";

    LetsFill(Knunu_dirname, SIGNAL_Nevt, SIGNAL_Ncandidates, "B2Knunu", true, Scale_Kplus_validation_MC15rd);
    LetsFill(Kstarnunu_dirname, SIGNAL_Nevt, SIGNAL_Ncandidates, "B2Kstarnunu", true, Scale_Kplusstar_validation_MC15rd);
    LetsFill(Xsununu_dirname, SIGNAL_Nevt, SIGNAL_Ncandidates, "B2Xsnunu", true, Scale_Xsu_nonresonant_validation_MC15rd);
    LetsFill(K0nunu_dirname, SIGNAL_Nevt, SIGNAL_Ncandidates, "B02K0nunu", false, Scale_K0_validation_MC15rd);
    LetsFill(K0starnunu_dirname, SIGNAL_Nevt, SIGNAL_Ncandidates, "B02Kstar0nunu", false, Scale_K0star_validation_MC15rd);
    LetsFill(Xsdnunu_dirname, SIGNAL_Nevt, SIGNAL_Ncandidates, "B02Xsnunu", false, Scale_Xsd_nonresonant_validation_MC15rd);

    SIGNAL_multiplicity->Divide(SIGNAL_Ncandidates, SIGNAL_Nevt);

    TCanvas* c_temp = new TCanvas("c", "", 1500, 1200); c_temp->cd();

    SIGNAL_multiplicity->SetMarkerSize(0.0);
    SIGNAL_multiplicity->SetLineWidth(2);
    SIGNAL_multiplicity->SetMarkerColor(4);
    SIGNAL_multiplicity->SetMarkerStyle(21);
    SIGNAL_multiplicity->SetStats(false);

    SIGNAL_multiplicity->Draw("P");

    c_temp->SaveAs("Plot_multiplicity_vs_MXstrue.png");

    for (int i = 0; i < NBins; i++) printf("%lf ", SIGNAL_multiplicity->GetBinContent(i + 1));

//    free(line);
//    gPad->BuildLegend();
//    c_temp->SaveAs("Plot_with_legend.png");
}
