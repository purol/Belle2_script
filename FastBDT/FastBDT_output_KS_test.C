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

TH2F* OBB_Oqq_BKG_train = new TH2F("OBB_Oqq_BKG_train", ";O_{BB};O{qq}", 40, 0, 1.0, 40, 0, 1.0);
TH1F* OBB_BKG_train = new TH1F("OBB_BKG_train", ";O_{BB};", 40, 0, 1.0);
TH1F* Oqq_BKG_train = new TH1F("Oqq_BKG_train", ";O_{qq};", 40, 0, 1.0);
TH2F* OBB_Oqq_BKG_test = new TH2F("OBB_Oqq_BKG_test", ";O_{BB};O{qq}", 40, 0, 1.0, 40, 0, 1.0);
TH1F* OBB_BKG_test = new TH1F("OBB_BKG_test", ";O_{BB};", 40, 0, 1.0);
TH1F* Oqq_BKG_test = new TH1F("Oqq_BKG_test", ";O_{qq};", 40, 0, 1.0);

TH2F* OBB_Oqq_SIGNAL_train = new TH2F("OBB_Oqq_SIGNAL_train", ";O_{BB};O{qq}", 40, 0, 1.0, 40, 0, 1.0);
TH1F* OBB_SIGNAL_train = new TH1F("OBB_SIGNAL_train", ";O_{BB};", 40, 0, 1.0);
TH1F* Oqq_SIGNAL_train = new TH1F("Oqq_SIGNAL_train", ";O_{qq};", 40, 0, 1.0);
TH2F* OBB_Oqq_SIGNAL_test = new TH2F("OBB_Oqq_SIGNAL_test", ";O_{BB};O{qq}", 40, 0, 1.0, 40, 0, 1.0);
TH1F* OBB_SIGNAL_test = new TH1F("OBB_SIGNAL_test", ";O_{BB};", 40, 0, 1.0);
TH1F* Oqq_SIGNAL_test = new TH1F("Oqq_SIGNAL_test", ";O_{qq};", 40, 0, 1.0);

void LetsFill(const char* dirname, TH2F* OBB_Oqq_hist, TH1F* OBB_hist, TH1F* Oqq_hist, double weight_var = 1.0) {
    double OBB_var = 0;
    double Oqq_var = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("MVA_BB", &OBB_var);
        tree_Btag->SetBranchAddress("MVA_Continuum", &Oqq_var);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            
            OBB_Oqq_hist->Fill(OBB_var, Oqq_var, weight_var);
            OBB_hist->Fill(OBB_var, Oqq_var, weight_var);
            Oqq_hist->Fill(OBB_var, Oqq_var, weight_var);
        }
        input_file->Close();

    }

}

void FastBDT_output_KS_test()
{

    // get data from root files
    const char* train_dirname_Knunu = "./SIGNAL_analysis/train_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Knunu";
    const char* train_dirname_Kstarnunu = "./SIGNAL_analysis/train_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Kstarnunu";
    const char* train_dirname_Xsununu = "./SIGNAL_analysis/train_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Xsnunu";
    const char* train_dirname_K0nunu = "./SIGNAL_analysis/train_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0nunu";
    const char* train_dirname_K0starnunu = "./SIGNAL_analysis/train_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0starnunu";
    const char* train_dirname_Xsdnunu = "./SIGNAL_analysis/train_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02Xsnunu";
    LetsFill(train_dirname_Knunu, OBB_Oqq_SIGNAL_train, OBB_SIGNAL_train, Oqq_SIGNAL_train, Scale_Kplus);
    LetsFill(train_dirname_Kstarnunu, OBB_Oqq_SIGNAL_train, OBB_SIGNAL_train, Oqq_SIGNAL_train, Scale_Kplusstar);
    LetsFill(train_dirname_Xsununu, OBB_Oqq_SIGNAL_train, OBB_SIGNAL_train, Oqq_SIGNAL_train, Scale_Xsu_nonresonant);
    LetsFill(train_dirname_K0nunu, OBB_Oqq_SIGNAL_train, OBB_SIGNAL_train, Oqq_SIGNAL_train, Scale_K0);
    LetsFill(train_dirname_K0starnunu, OBB_Oqq_SIGNAL_train, OBB_SIGNAL_train, Oqq_SIGNAL_train, Scale_K0star);
    LetsFill(train_dirname_Xsdnunu, OBB_Oqq_SIGNAL_train, OBB_SIGNAL_train, Oqq_SIGNAL_train, Scale_Xsd_nonresonant);

    const char* train_dirname_CHG = "./CHG_analysis/train_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* train_dirname_MIX = "./MIX_analysis/train_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* train_dirname_UUBAR = "./UUBAR_analysis/train_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* train_dirname_DDBAR = "./DDBAR_analysis/train_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* train_dirname_SSBAR = "./SSBAR_analysis/train_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* train_dirname_CHARM = "./CHARM_analysis/train_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    LetsFill(train_dirname_CHG, OBB_Oqq_BKG_train, OBB_BKG_train, Oqq_BKG_train);
    LetsFill(train_dirname_MIX, OBB_Oqq_BKG_train, OBB_BKG_train, Oqq_BKG_train);
    LetsFill(train_dirname_UUBAR, OBB_Oqq_BKG_train, OBB_BKG_train, Oqq_BKG_train);
    LetsFill(train_dirname_DDBAR, OBB_Oqq_BKG_train, OBB_BKG_train, Oqq_BKG_train);
    LetsFill(train_dirname_CHARM, OBB_Oqq_BKG_train, OBB_BKG_train, Oqq_BKG_train);

    const char* test_dirname_Knunu = "./SIGNAL_analysis/test_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Knunu";
    const char* test_dirname_Kstarnunu = "./SIGNAL_analysis/test_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Kstarnunu";
    const char* test_dirname_Xsununu = "./SIGNAL_analysis/test_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Xsnunu";
    const char* test_dirname_K0nunu = "./SIGNAL_analysis/test_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0nunu";
    const char* test_dirname_K0starnunu = "./SIGNAL_analysis/test_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0starnunu";
    const char* test_dirname_Xsdnunu = "./SIGNAL_analysis/test_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02Xsnunu";
    LetsFill(test_dirname_Knunu, OBB_Oqq_SIGNAL_test, OBB_SIGNAL_test, Oqq_SIGNAL_test, Scale_Kplus);
    LetsFill(test_dirname_Kstarnunu, OBB_Oqq_SIGNAL_test, OBB_SIGNAL_test, Oqq_SIGNAL_test, Scale_Kplusstar);
    LetsFill(test_dirname_Xsununu, OBB_Oqq_SIGNAL_test, OBB_SIGNAL_test, Oqq_SIGNAL_test, Scale_Xsu_nonresonant);
    LetsFill(test_dirname_K0nunu, OBB_Oqq_SIGNAL_test, OBB_SIGNAL_test, Oqq_SIGNAL_test, Scale_K0);
    LetsFill(test_dirname_K0starnunu, OBB_Oqq_SIGNAL_test, OBB_SIGNAL_test, Oqq_SIGNAL_test, Scale_K0star);
    LetsFill(test_dirname_Xsdnunu, OBB_Oqq_SIGNAL_test, OBB_SIGNAL_test, Oqq_SIGNAL_test, Scale_Xsd_nonresonant);

    const char* test_dirname_CHG = "./CHG_analysis/validation_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* test_dirname_MIX = "./MIX_analysis/validation_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* test_dirname_UUBAR = "./UUBAR_analysis/validation_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* test_dirname_DDBAR = "./DDBAR_analysis/validation_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* test_dirname_SSBAR = "./SSBAR_analysis/validation_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* test_dirname_CHARM = "./CHARM_analysis/validation_v002/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    LetsFill(test_dirname_CHG, OBB_Oqq_BKG_test, OBB_BKG_test, Oqq_BKG_test);
    LetsFill(test_dirname_MIX, OBB_Oqq_BKG_test, OBB_BKG_test, Oqq_BKG_test);
    LetsFill(test_dirname_UUBAR, OBB_Oqq_BKG_test, OBB_BKG_test, Oqq_BKG_test);
    LetsFill(test_dirname_DDBAR, OBB_Oqq_BKG_test, OBB_BKG_test, Oqq_BKG_test);
    LetsFill(test_dirname_CHARM, OBB_Oqq_BKG_test, OBB_BKG_test, Oqq_BKG_test);

    // normalization
    OBB_Oqq_SIGNAL_train->Scale(factor / OBB_Oqq_SIGNAL_train->Integral(), "width");
    OBB_SIGNAL_train->Scale(factor / OBB_SIGNAL_train->Integral(), "width");
    Oqq_SIGNAL_train->Scale(factor / Oqq_SIGNAL_train->Integral(), "width");

    OBB_Oqq_SIGNAL_test->Scale(factor / OBB_Oqq_SIGNAL_test->Integral(), "width");
    OBB_SIGNAL_test->Scale(factor / OBB_SIGNAL_test->Integral(), "width");
    Oqq_SIGNAL_test->Scale(factor / Oqq_SIGNAL_test->Integral(), "width");

    OBB_Oqq_BKG_train->Scale(factor / OBB_Oqq_BKG_train->Integral(), "width");
    OBB_BKG_train->Scale(factor / OBB_BKG_train->Integral(), "width");
    Oqq_BKG_train->Scale(factor / Oqq_BKG_train->Integral(), "width");

    OBB_Oqq_BKG_test->Scale(factor / OBB_Oqq_BKG_test->Integral(), "width");
    OBB_BKG_test->Scale(factor / OBB_BKG_test->Integral(), "width");
    Oqq_BKG_test->Scale(factor / Oqq_BKG_test->Integral(), "width");

    // set color (BKG: kRed, SIGNAL: kBlue)
    // solid line: test, histogram: traing
    OBB_SIGNAL_test->SetMarkerStyle(8);
    OBB_SIGNAL_test->SetLineColor(kBlue);
    OBB_SIGNAL_test->SetLineWidth(1);

    Oqq_SIGNAL_test->SetMarkerStyle(8);
    Oqq_SIGNAL_test->SetLineColor(kBlue);
    Oqq_SIGNAL_test->SetLineWidth(1);

    OBB_BKG_test->SetMarkerStyle(8);
    OBB_BKG_test->SetLineColor(kRed);
    OBB_BKG_test->SetLineWidth(1);

    Oqq_BKG_test->SetMarkerStyle(8);
    Oqq_BKG_test->SetLineColor(kRed);
    Oqq_BKG_test->SetLineWidth(1);

    OBB_SIGNAL_train->SetFillStyle(3004);
    OBB_SIGNAL_train->SetLineColor(kBlue);
    OBB_SIGNAL_train->SetFillColor(kBlue);

    Oqq_SIGNAL_train->SetFillStyle(3004);
    Oqq_SIGNAL_train->SetLineColor(kBlue);
    Oqq_SIGNAL_train->SetFillColor(kBlue);

    OBB_BKG_train->SetFillStyle(3005);
    OBB_BKG_train->SetLineColor(kRed);
    OBB_BKG_train->SetFillColor(kRed);

    Oqq_BKG_train->SetFillStyle(3005);
    Oqq_BKG_train->SetLineColor(kRed);
    Oqq_BKG_train->SetFillColor(kRed);

    // draw plot and get p-value
    double p_value_SIGNAL = OBB_Oqq_SIGNAL_test->KolmogorovTest(OBB_Oqq_SIGNAL_train);
    double p_value_BKG = OBB_Oqq_BKG_test->KolmogorovTest(OBB_Oqq_BKG_train);
    printf("p value SIGNAL: %lf\n", p_value_SIGNAL);
    printf("p value BKG: %lf\n", p_value_BKG);

    TCanvas* c_temp = new TCanvas("c", "", 600, 600); c_temp->cd();
    OBB_SIGNAL_train->Draw("Hist"); OBB_BKG_train->Draw("HistSAME");
    OBB_SIGNAL_test->Draw("SAME"); OBB_BKG_test->Draw("SAME");
    c_temp->SaveAs("OBB_Plot.png");

    TCanvas* c_temp_2 = new TCanvas("c2", "", 600, 600); c_temp->cd();
    Oqq_SIGNAL_train->Draw("Hist"); Oqq_BKG_train->Draw("HistSAME");
    Oqq_SIGNAL_test->Draw("SAME"); Oqq_BKG_test->Draw("SAME");
    c_temp->SaveAs("Oqq_Plot.png");

}
