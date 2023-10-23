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
TH1F* Eecl_DATA = new TH1F("Eecl_DATA", ";E_{ecl} [GeV];", 15, 0, 4.0);
TH1F* Eecl_MC = new TH1F("Eecl_MC", ";E_{ecl} [GeV];", 15, 0, 4.0);
TH1F* Eecl_MC_SIGNAL = new TH1F("Eecl_MC_SIGNAL", ";E_{ecl} [GeV];", 15, 0, 4.0);
TH1F* Eecl_MC_CHG = new TH1F("Eecl_MC_CHG", ";E_{ecl} [GeV];", 15, 0, 4.0);
TH1F* Eecl_MC_MIX = new TH1F("Eecl_MC_MIX", ";E_{ecl} [GeV];", 15, 0, 4.0);
TH1F* Eecl_MC_UUBAR = new TH1F("Eecl_MC_UUBAR", ";E_{ecl} [GeV];", 15, 0, 4.0);
TH1F* Eecl_MC_DDBAR = new TH1F("Eecl_MC_DDBAR", ";E_{ecl} [GeV];", 15, 0, 4.0);
TH1F* Eecl_MC_SSBAR = new TH1F("Eecl_MC_SSBAR", ";E_{ecl} [GeV];", 15, 0, 4.0);
TH1F* Eecl_MC_CHARM = new TH1F("Eecl_MC_CHARM", ";E_{ecl} [GeV];", 15, 0, 4.0);

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

    // get data from root files
    const char* MC_dirname_Knunu = "./SIGNAL_analysis/test_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Knunu";
    const char* MC_dirname_Kstarnunu = "./SIGNAL_analysis/test_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Kstarnunu";
    const char* MC_dirname_Xsununu = "./SIGNAL_analysis/test_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Xsnunu";
    const char* MC_dirname_K0nunu = "./SIGNAL_analysis/test_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0nunu";
    const char* MC_dirname_K0starnunu = "./SIGNAL_analysis/test_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0starnunu";
    const char* MC_dirname_Xsdnunu = "./SIGNAL_analysis/test_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02Xsnunu";
    LetsFill(MC_dirname_Knunu, Eecl_MC_SIGNAL, Scale_Kplus);
    LetsFill(MC_dirname_Kstarnunu, Eecl_MC_SIGNAL, Scale_Kplusstar);
    LetsFill(MC_dirname_Xsununu, Eecl_MC_SIGNAL, Scale_Xsu_nonresonant);
    LetsFill(MC_dirname_K0nunu, Eecl_MC_SIGNAL, Scale_K0);
    LetsFill(MC_dirname_K0starnunu, Eecl_MC_SIGNAL, Scale_K0star);
    LetsFill(MC_dirname_Xsdnunu, Eecl_MC_SIGNAL, Scale_Xsd_nonresonant);
    LetsFill(MC_dirname_Knunu, Eecl_MC, Scale_Kplus);
    LetsFill(MC_dirname_Kstarnunu, Eecl_MC, Scale_Kplusstar);
    LetsFill(MC_dirname_Xsununu, Eecl_MC, Scale_Xsu_nonresonant);
    LetsFill(MC_dirname_K0nunu, Eecl_MC, Scale_K0);
    LetsFill(MC_dirname_K0starnunu, Eecl_MC, Scale_K0star);
    LetsFill(MC_dirname_Xsdnunu, Eecl_MC, Scale_Xsd_nonresonant);

    const char* MC_dirname_CHG = "./CHG_analysis/test_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_MIX = "./MIX_analysis/test_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_UUBAR = "./UUBAR_analysis/test_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_DDBAR = "./DDBAR_analysis/test_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_SSBAR = "./SSBAR_analysis/test_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_CHARM = "./CHARM_analysis/test_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
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

    const char* DATA_dirname_Knunu = "./SIGNAL_analysis/validation_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Knunu";
    const char* DATA_dirname_Kstarnunu = "./SIGNAL_analysis/validation_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Kstarnunu";
    const char* DATA_dirname_Xsununu = "./SIGNAL_analysis/validation_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Xsnunu";
    const char* DATA_dirname_K0nunu = "./SIGNAL_analysis/validation_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0nunu";
    const char* DATA_dirname_K0starnunu = "./SIGNAL_analysis/validation_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0starnunu";
    const char* DATA_dirname_Xsdnunu = "./SIGNAL_analysis/validation_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02Xsnunu";
    LetsFill(DATA_dirname_Knunu, Eecl_DATA, Scale_Kplus);
    LetsFill(DATA_dirname_Kstarnunu, Eecl_DATA, Scale_Kplusstar);
    LetsFill(DATA_dirname_Xsununu, Eecl_DATA, Scale_Xsu_nonresonant);
    LetsFill(DATA_dirname_K0nunu, Eecl_DATA, Scale_K0);
    LetsFill(DATA_dirname_K0starnunu, Eecl_DATA, Scale_K0star);
    LetsFill(DATA_dirname_Xsdnunu, Eecl_DATA, Scale_Xsd_nonresonant);

    const char* DATA_dirname_CHG = "./CHG_analysis/validation_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_MIX = "./MIX_analysis/validation_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_UUBAR = "./UUBAR_analysis/validation_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_DDBAR = "./DDBAR_analysis/validation_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_SSBAR = "./SSBAR_analysis/validation_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_CHARM = "./CHARM_analysis/validation_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    LetsFill(DATA_dirname_CHG, Eecl_DATA);
    LetsFill(DATA_dirname_MIX, Eecl_DATA);
    LetsFill(DATA_dirname_UUBAR, Eecl_DATA);
    LetsFill(DATA_dirname_DDBAR, Eecl_DATA);
    LetsFill(DATA_dirname_CHARM, Eecl_DATA);

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
