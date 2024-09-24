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

TH1D* SIGNAL_hist = new TH1D("signal", ";##Upsilon(4S) multiplicity;Fraction of events", 5, 0.5, 5.5);
TH1D* CHG_hist = new TH1D("charged", ";#Upsilon(4S) multiplicity;Fraction of events", 5, 0.5, 5.5);
TH1D* MIX_hist = new TH1D("mixed", ";#Upsilon(4S) multiplicity;Fraction of events", 5, 0.5, 5.5);
TH1D* UUBAR_hist = new TH1D("u#bar{u}", ";#Upsilon(4S) multiplicity;Fraction of events", 5, 0.5, 5.5);
TH1D* DDBAR_hist = new TH1D("d#bar{d}", ";#Upsilon(4S) multiplicity;Fraction of events", 5, 0.5, 5.5);
TH1D* SSBAR_hist = new TH1D("s#bar{s}", ";#Upsilon(4S) multiplicity;Fraction of events", 5, 0.5, 5.5);
TH1D* CHARM_hist = new TH1D("c#bar{c}", ";#Upsilon(4S) multiplicity;Fraction of events", 5, 0.5, 5.5);

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

        int Ncandidates = 0;

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            // D veto
            if ((MXs > 1.84) && (MXs < 1.89)) continue;

            if ((-1 == PREV__experiment__) && (-1 == PREV__run__) && (-1 == PREV__event__) && (-1 == PREV__ncandidates__)) {
                Ncandidates = 1;
            }
            else if ((__experiment__ == PREV__experiment__) && (__run__ == PREV__run__) && (__event__ == PREV__event__) && (__ncandidates__ == PREV__ncandidates__)) {
                Ncandidates++;
            }
            else {
                hist->Fill(Ncandidates, weight);
                Ncandidates = 1;
            }

            PREV__experiment__ = __experiment__;
            PREV__run__ = __run__;
            PREV__event__ = __event__;
            PREV__ncandidates__ = __ncandidates__;

        }

        // read the final entry
        if(tree_upsilon->GetEntries() != 0) hist->Fill(Ncandidates, weight);

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

        int Ncandidates = 0;

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            // D veto
            if ((MXs > 1.84) && (MXs < 1.89)) continue;

            if ((-1 == PREV__experiment__) && (-1 == PREV__run__) && (-1 == PREV__event__) && (-1 == PREV__ncandidates__)) {
                Ncandidates = 1;
            }
            else if ((__experiment__ == PREV__experiment__) && (__run__ == PREV__run__) && (__event__ == PREV__event__) && (__ncandidates__ == PREV__ncandidates__)) {
                Ncandidates++;
            }
            else {
                hist->Fill(Ncandidates, weight);
                Ncandidates = 1;
            }

            PREV__experiment__ = __experiment__;
            PREV__run__ = __run__;
            PREV__event__ = __event__;
            PREV__ncandidates__ = __ncandidates__;

        }

        // read the final entry
        if (Ncandidates != 0) hist->Fill(Ncandidates, weight);

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
    const char* CHG_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/SatoriRD/CHG_analysis/validation_v004/before_Dveto_cut";
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
    legend->Draw();

    c_temp->SaveAs("Plot_multiplicity.png");

//    free(line);
//    gPad->BuildLegend();
//    c_temp->SaveAs("Plot_with_legend.png");
}
