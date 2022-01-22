#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit ;

void load_files(const char* dirname, std::vector<string>* names) {
    TSystemDirectory dir(dirname, dirname);
    TList* files = dir.GetListOfFiles();
    if (files) {
        TSystemFile* file;
        TString fname;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(".root")) {
                names->push_back(fname.Data());
            }
        }
    }
}

THStack* Stack = new THStack("thstack", ";E_{ecl} [GeV];evt");
TH1F* Eecl_DATA = new TH1F("Eecl_DATA", ";E_{ecl} [GeV];", 12, 0, 1.2);
TH1F* Eecl_MC = new TH1F("Eecl_MC", ";E_{ecl} [GeV];", 12, 0, 1.2);
TH1F* Eecl_MC_SIGNAL = new TH1F("Eecl_MC_SIGNAL", ";E_{ecl} [GeV];", 12, 0, 1.2);
TH1F* Eecl_MC_CHG = new TH1F("Eecl_MC_CHG", ";E_{ecl} [GeV];", 12, 0, 1.2);
TH1F* Eecl_MC_MIX = new TH1F("Eecl_MC_MIX", ";E_{ecl} [GeV];", 12, 0, 1.2);
TH1F* Eecl_MC_UUBAR = new TH1F("Eecl_MC_UUBAR", ";E_{ecl} [GeV];", 12, 0, 1.2);
TH1F* Eecl_MC_DDBAR = new TH1F("Eecl_MC_DDBAR", ";E_{ecl} [GeV];", 12, 0, 1.2);
TH1F* Eecl_MC_SSBAR = new TH1F("Eecl_MC_SSBAR", ";E_{ecl} [GeV];", 12, 0, 1.2);
TH1F* Eecl_MC_CHARM = new TH1F("Eecl_MC_CHARM", ";E_{ecl} [GeV];", 12, 0, 1.2);

void LetsFill(const char* dirname, TH1F* hist, double weight_var = 1.0) {
    double Eecl_var = 0;
    double Mbc_var = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("roeEextra__bocleanMask__bc", &Eecl_var); // Eecl
        tree_Btag->SetBranchAddress("Btag_Mbc", &Mbc_var); // Mbc

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            
            hist->Fill(Eecl_var, weight_var);
        }
        input_file->Close();

    }

}

void KS_test()
{
    const double weight = 0.003385;

    // get data from root files
    const char* MC_dirname_signal = "./output_after_TMVA_test/SIGNAL";
    LetsFill(MC_dirname_signal, Eecl_MC_SIGNAL, weight);
    LetsFill(MC_dirname_signal, Eecl_MC, weight);

    const char* MC_dirname_CHG = "./output_after_TMVA_test/CHG";
    const char* MC_dirname_MIX = "./output_after_TMVA_test/MIX";
    const char* MC_dirname_UUBAR = "./output_after_TMVA_test/UUBAR";
    const char* MC_dirname_DDBAR = "./output_after_TMVA_test/DDBAR";
    const char* MC_dirname_SSBAR = "./output_after_TMVA_test/SSBAR";
    const char* MC_dirname_CHARM = "./output_after_TMVA_test/CHARM";
    LetsFill(MC_dirname_CHG, Eecl_MC_CHG);
    LetsFill(MC_dirname_MIX, Eecl_MC_MIX);
    LetsFill(MC_dirname_UUBAR, Eecl_MC_UUBAR);
    LetsFill(MC_dirname_DDBAR, Eecl_MC_DDBAR);
    LetsFill(MC_dirname_CHARM, Eecl_MC_CHARM);
    LetsFill(MC_dirname_CHG, Eecl_MC);
    LetsFill(MC_dirname_MIX, Eecl_MC);
    LetsFill(MC_dirname_UUBAR, Eecl_MC);
    LetsFill(MC_dirname_DDBAR, Eecl_MC);
    LetsFill(MC_dirname_CHARM, Eecl_MC);

    const char* DATA_dirname_BKG = "./output_after_TMVA_train/BKG";
    const char* DATA_dirname_SIGNAL = "./output_after_TMVA_train/SIGNAL";
    LetsFill(DATA_dirname_BKG, Eecl_DATA);
    LetsFill(DATA_dirname_SIGNAL, Eecl_DATA, weight);

    Eecl_DATA->SetMarkerStyle(8);
    Eecl_DATA->SetLineColor(1);
    Eecl_DATA->SetLineWidth(1);

    Stack->Add(Eecl_MC_SIGNAL);
    Stack->Add(Eecl_MC_CHG);
    Stack->Add(Eecl_MC_MIX);
    Stack->Add(Eecl_MC_UUBAR);
    Stack->Add(Eecl_MC_DDBAR);
    Stack->Add(Eecl_MC_SSBAR);
    Stack->Add(Eecl_MC_CHARM);
    Stack->SetMaximum(28.0);

    TCanvas* c_temp = new TCanvas("c", "", 600, 600); c_temp->cd();
    //gStyle->SetPalette(kBlackBody);
    Eecl_MC_SIGNAL->SetLineColor(kRed - 7);
    Eecl_MC_SIGNAL->SetFillColor(kRed - 7);
    Eecl_MC_CHG->SetLineColor(kYellow - 7);
    Eecl_MC_CHG->SetFillColor(kYellow - 7);
    Eecl_MC_MIX->SetLineColor(kGreen - 7);
    Eecl_MC_MIX->SetFillColor(kGreen - 7);
    Eecl_MC_UUBAR->SetLineColor(kBlue - 7);
    Eecl_MC_UUBAR->SetFillColor(kBlue - 7);
    Eecl_MC_DDBAR->SetLineColor(kMagenta - 7);
    Eecl_MC_DDBAR->SetFillColor(kMagenta - 7);
    Eecl_MC_SSBAR->SetLineColor(kOrange - 7);
    Eecl_MC_SSBAR->SetFillColor(kOrange - 7);
    Eecl_MC_CHARM->SetLineColor(kTeal - 7);
    Eecl_MC_CHARM->SetFillColor(kTeal - 7);

    double p_value = Eecl_DATA->KolmogorovTest(Eecl_MC);
    printf("p value: %lf\n", p_value);

    Stack->Draw("Hist"); Eecl_DATA->Draw("SAME");
    c_temp->SaveAs("Plot.png");

    gPad->BuildLegend();
    c_temp->SaveAs("Plot_with_legend.png");

}
