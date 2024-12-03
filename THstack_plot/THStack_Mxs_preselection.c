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

#include "TLine.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "THStack.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"

#include <string>

#include "constants.h"

using namespace std;

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

THStack* Stack = new THStack("thstack", ";M_{Xs}^{reco} [GeV];number of candidates");
TH1D* SIGNAL_hist = new TH1D("signal #times 10^{5}", ";M_{Xs}^{reco} [GeV];number of candidates", 70, 0.4, 2.8);
TH1D* CHG_hist = new TH1D("charged", ";M_{Xs}^{reco} [GeV];number of candidates", 70, 0.4, 2.8);
TH1D* MIX_hist = new TH1D("mixed", ";M_{Xs}^{reco} [GeV];number of candidates", 70, 0.4, 2.8);
TH1D* UUBAR_hist = new TH1D("u#bar{u}", ";M_{Xs}^{reco} [GeV];number of candidates", 70, 0.4, 2.8);
TH1D* DDBAR_hist = new TH1D("d#bar{d}", ";M_{Xs}^{reco} [GeV];number of candidates", 70, 0.4, 2.8);
TH1D* SSBAR_hist = new TH1D("s#bar{s}", ";M_{Xs}^{reco} [GeV];number of candidates", 70, 0.4, 2.8);
TH1D* CHARM_hist = new TH1D("c#bar{c}", ";M_{Xs}^{reco} [GeV];number of candidates", 70, 0.4, 2.8);

void LetsFill(const char* dirname, TH1D* hist, double weight = 1) {
    double var = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_Xsu = (TTree*)input_file->Get("Xsu_preselection");
        TTree* tree_Xsd = (TTree*)input_file->Get("Xsd_preselection");

        tree_Xsu->SetBranchAddress("M", &var);
        printf("%lld entries...\n", tree_Xsu->GetEntries());
        for (unsigned int j = 0; j < tree_Xsu->GetEntries(); j++) { // Fill
            tree_Xsu->GetEntry(j);
            hist->Fill(var, weight);
        }

        tree_Xsd->SetBranchAddress("M", &var);
        printf("%lld entries...\n", tree_Xsd->GetEntries());
        for (unsigned int j = 0; j < tree_Xsd->GetEntries(); j++) { // Fill
            tree_Xsd->GetEntry(j);
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

        TTree* tree_Xsu = (TTree*)input_file->Get("Xsu_preselection");
        TTree* tree_Xsd = (TTree*)input_file->Get("Xsd_preselection");

        tree_Xsu->SetBranchAddress("M", &var);
        printf("%lld entries...\n", tree_Xsu->GetEntries());
        for (unsigned int j = 0; j < tree_Xsu->GetEntries(); j++) { // Fill
            tree_Xsu->GetEntry(j);
            hist->Fill(var, weight);
        }

        tree_Xsd->SetBranchAddress("M", &var);
        printf("%lld entries...\n", tree_Xsd->GetEntries());
        for (unsigned int j = 0; j < tree_Xsd->GetEntries(); j++) { // Fill
            tree_Xsd->GetEntry(j);
            hist->Fill(var, weight);
        }

        input_file->Close();

    }

}

int main() {

    const char* Knunu_dirname = "/home/belle2/junewoo/storage_ghi/Ntuple/KumoiRD/SIGNAL/validation";
    const char* Kstarnunu_dirname = "/home/belle2/junewoo/storage_ghi/Ntuple/KumoiRD/SIGNAL/validation";
    const char* Xsununu_dirname = "/home/belle2/junewoo/storage_ghi/Ntuple/KumoiRD/SIGNAL/validation";
    const char* K0nunu_dirname = "/home/belle2/junewoo/storage_ghi/Ntuple/KumoiRD/SIGNAL/validation";
    const char* K0starnunu_dirname = "/home/belle2/junewoo/storage_ghi/Ntuple/KumoiRD/SIGNAL/validation";
    const char* Xsdnunu_dirname = "/home/belle2/junewoo/storage_ghi/Ntuple/KumoiRD/SIGNAL/validation";
    const char* CHG_dirname = "/home/belle2/junewoo/storage_ghi/Ntuple/KumoiRD/CHG/validation";
    const char* MIX_dirname = "/home/belle2/junewoo/storage_ghi/Ntuple/KumoiRD/MIX/validation";
    const char* UUBAR_dirname = "/home/belle2/junewoo/storage_ghi/Ntuple/KumoiRD/UUBAR/validation";
    const char* DDBAR_dirname = "/home/belle2/junewoo/storage_ghi/Ntuple/KumoiRD/DDBAR/validation";
    const char* SSBAR_dirname = "/home/belle2/junewoo/storage_ghi/Ntuple/KumoiRD/SSBAR/validation";
    const char* CHARM_dirname = "/home/belle2/junewoo/storage_ghi/Ntuple/KumoiRD/CHARM/validation";

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

    double BKG_int = CHG_int + MIX_int;

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

    SIGNAL_hist->Scale(5000.0);
    //SIGNAL_hist->Scale(1.0 / SIGNAL_int, "width");
    SIGNAL_hist->SetLineWidth(3);
    SIGNAL_hist->SetLineColor(2);
    SIGNAL_hist->SetFillStyle(0);

    TCanvas* c_temp = new TCanvas("c", "", 1500, 1200); c_temp->cd();

    gStyle->SetPalette(kGistEarth);

    Stack->Draw("pfc Hist"); SIGNAL_hist->Draw("HistSAME");

    TLegend* legend = gPad->BuildLegend(0.9, 0.9, 0.7, 0.7);
    legend->SetFillStyle(0); legend->SetLineWidth(0);

    Float_t ymax = Stack->GetMaximum();
    TLine* lineA = new TLine(2.0, 0, 2.0, ymax);
    lineA->SetLineColor(kRed);
    lineA->SetLineStyle(5);
    lineA->Draw();
    c_temp->SaveAs("Plot_Mxs_preselection.png");

//    free(lineA);
//    free(lineB);
//    gPad->BuildLegend();
//    c_temp->SaveAs("Plot_with_legend.png");

    return 0;
}
