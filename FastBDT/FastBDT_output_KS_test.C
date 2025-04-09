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

#include "constants.h"

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
TH1F* OBB_BKG_train = new TH1F("FastBDT BKG train", ";FastBDT;Arbitrary unit", 40, 0, 1.0);
TH1F* Oqq_BKG_train = new TH1F("FastBDT BKG train", ";FastBDT;Arbitrary unit", 40, 0, 1.0);
TH2F* OBB_Oqq_BKG_test = new TH2F("OBB_Oqq_BKG_test", ";O_{BB};O{qq}", 40, 0, 1.0, 40, 0, 1.0);
TH1F* OBB_BKG_test = new TH1F("FastBDT BKG test", ";FastBDT;Arbitrary unit", 40, 0, 1.0);
TH1F* Oqq_BKG_test = new TH1F("FastBDT BKG test", ";FastBDT;Arbitrary unit", 40, 0, 1.0);

TH2F* OBB_Oqq_SIGNAL_train = new TH2F("OBB_Oqq_SIGNAL_train", ";O_{BB};O{qq}", 40, 0, 1.0, 40, 0, 1.0);
TH1F* OBB_SIGNAL_train = new TH1F("FastBDT SIGNAL train", ";FastBDT;Arbitrary unit", 40, 0, 1.0);
TH1F* Oqq_SIGNAL_train = new TH1F("FastBDT SIGNAL train", ";FastBDT;Arbitrary unit", 40, 0, 1.0);
TH2F* OBB_Oqq_SIGNAL_test = new TH2F("OBB_Oqq_SIGNAL_test", ";O_{BB};O{qq}", 40, 0, 1.0, 40, 0, 1.0);
TH1F* OBB_SIGNAL_test = new TH1F("FastBDT SIGNAL test", ";FastBDT;Arbitrary unit", 40, 0, 1.0);
TH1F* Oqq_SIGNAL_test = new TH1F("FastBDT SIGNAL test", ";FastBDT;Arbitrary unit", 40, 0, 1.0);

void LetsFill(const char* dirname, TH2F* OBB_Oqq_hist, TH1F* OBB_hist, TH1F* Oqq_hist, double weight_var = 1.0) {
    float OBB_var = 0;
    float Oqq_var = 0;
    double Mx = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("MVA_BB", &OBB_var);
        tree_upsilon->SetBranchAddress("MVA_Continuum", &Oqq_var);
        tree_Bsig->SetBranchAddress("Bsig_M", &Mx);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            
            OBB_Oqq_hist->Fill(OBB_var, Oqq_var, weight_var);
            OBB_hist->Fill(OBB_var, weight_var);
            Oqq_hist->Fill(Oqq_var, weight_var);
        }
        input_file->Close();

    }

}

void FastBDT_output_KS_test()
{

    // get data from root files
    const char* train_dirname_Knunu = "./temp_v008_train/B2Knunu";
    const char* train_dirname_Kstarnunu = "./temp_v008_train/B2Kstarnunu";
    const char* train_dirname_Xsununu = "./temp_v008_train/B2Xsnunu";
    const char* train_dirname_K0nunu = "./temp_v008_train/B02K0nunu";
    const char* train_dirname_K0starnunu = "./temp_v008_train/B02K0starnunu";
    const char* train_dirname_Xsdnunu = "./temp_v008_train/B02Xsnunu";
    LetsFill(train_dirname_Knunu, OBB_Oqq_SIGNAL_train, OBB_SIGNAL_train, Oqq_SIGNAL_train, Scale_Kplus_train);
    LetsFill(train_dirname_Kstarnunu, OBB_Oqq_SIGNAL_train, OBB_SIGNAL_train, Oqq_SIGNAL_train, Scale_Kplusstar_train);
    LetsFill(train_dirname_Xsununu, OBB_Oqq_SIGNAL_train, OBB_SIGNAL_train, Oqq_SIGNAL_train, Scale_Xsu_nonresonant_train);
    LetsFill(train_dirname_K0nunu, OBB_Oqq_SIGNAL_train, OBB_SIGNAL_train, Oqq_SIGNAL_train, Scale_K0_train);
    LetsFill(train_dirname_K0starnunu, OBB_Oqq_SIGNAL_train, OBB_SIGNAL_train, Oqq_SIGNAL_train, Scale_K0star_train);
    LetsFill(train_dirname_Xsdnunu, OBB_Oqq_SIGNAL_train, OBB_SIGNAL_train, Oqq_SIGNAL_train, Scale_Xsd_nonresonant_train);

    const char* train_dirname_CHG = "./temp_v008_train/CHG";
    const char* train_dirname_MIX = "./temp_v008_train/MIX";
    const char* train_dirname_UUBAR = "./temp_v008_train/UUBAR";
    const char* train_dirname_DDBAR = "./temp_v008_train/DDBAR";
    const char* train_dirname_SSBAR = "./temp_v008_train/SSBAR";
    const char* train_dirname_CHARM = "./temp_v008_train/CHARM";
    LetsFill(train_dirname_CHG, OBB_Oqq_BKG_train, OBB_BKG_train, Oqq_BKG_train, Scale_CHG_train);
    LetsFill(train_dirname_MIX, OBB_Oqq_BKG_train, OBB_BKG_train, Oqq_BKG_train, Scale_MIX_train);
    LetsFill(train_dirname_UUBAR, OBB_Oqq_BKG_train, OBB_BKG_train, Oqq_BKG_train, Scale_UUBAR_train);
    LetsFill(train_dirname_DDBAR, OBB_Oqq_BKG_train, OBB_BKG_train, Oqq_BKG_train, Scale_DDBAR_train);
    LetsFill(train_dirname_SSBAR, OBB_Oqq_BKG_train, OBB_BKG_train, Oqq_BKG_train, Scale_SSBAR_train);
    LetsFill(train_dirname_CHARM, OBB_Oqq_BKG_train, OBB_BKG_train, Oqq_BKG_train, Scale_CHARM_train);

    const char* test_dirname_Knunu = "./temp_v008/B2Knunu";
    const char* test_dirname_Kstarnunu = "./temp_v008/B2Kstarnunu";
    const char* test_dirname_Xsununu = "./temp_v008/B2Xsnunu";
    const char* test_dirname_K0nunu = "./temp_v008/B02K0nunu";
    const char* test_dirname_K0starnunu = "./temp_v008/B02K0starnunu";
    const char* test_dirname_Xsdnunu = "./temp_v008/B02Xsnunu";
    LetsFill(test_dirname_Knunu, OBB_Oqq_SIGNAL_test, OBB_SIGNAL_test, Oqq_SIGNAL_test, Scale_Kplus_test);
    LetsFill(test_dirname_Kstarnunu, OBB_Oqq_SIGNAL_test, OBB_SIGNAL_test, Oqq_SIGNAL_test, Scale_Kplusstar_test);
    LetsFill(test_dirname_Xsununu, OBB_Oqq_SIGNAL_test, OBB_SIGNAL_test, Oqq_SIGNAL_test, Scale_Xsu_nonresonant_test);
    LetsFill(test_dirname_K0nunu, OBB_Oqq_SIGNAL_test, OBB_SIGNAL_test, Oqq_SIGNAL_test, Scale_K0_test);
    LetsFill(test_dirname_K0starnunu, OBB_Oqq_SIGNAL_test, OBB_SIGNAL_test, Oqq_SIGNAL_test, Scale_K0star_test);
    LetsFill(test_dirname_Xsdnunu, OBB_Oqq_SIGNAL_test, OBB_SIGNAL_test, Oqq_SIGNAL_test, Scale_Xsd_nonresonant_test);

    const char* test_dirname_CHG = "./temp_v008/CHG";
    const char* test_dirname_MIX = "./temp_v008/MIX";
    const char* test_dirname_UUBAR = "./temp_v008/UUBAR";
    const char* test_dirname_DDBAR = "./temp_v008/DDBAR";
    const char* test_dirname_SSBAR = "./temp_v008/SSBAR";
    const char* test_dirname_CHARM = "./temp_v008/CHARM";
    LetsFill(test_dirname_CHG, OBB_Oqq_BKG_test, OBB_BKG_test, Oqq_BKG_test, Scale_CHG_test);
    LetsFill(test_dirname_MIX, OBB_Oqq_BKG_test, OBB_BKG_test, Oqq_BKG_test, Scale_MIX_test);
    LetsFill(test_dirname_UUBAR, OBB_Oqq_BKG_test, OBB_BKG_test, Oqq_BKG_test, Scale_UUBAR_test);
    LetsFill(test_dirname_DDBAR, OBB_Oqq_BKG_test, OBB_BKG_test, Oqq_BKG_test, Scale_DDBAR_test);
    LetsFill(test_dirname_SSBAR, OBB_Oqq_BKG_test, OBB_BKG_test, Oqq_BKG_test, Scale_SSBAR_test);
    LetsFill(test_dirname_CHARM, OBB_Oqq_BKG_test, OBB_BKG_test, Oqq_BKG_test, Scale_CHARM_test);

    double factor = 1.0;

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
    OBB_SIGNAL_test->SetMarkerStyle(kFullCircle);
    OBB_SIGNAL_test->SetLineColor(kBlue);
    OBB_SIGNAL_test->SetMarkerColor(kBlue);
    OBB_SIGNAL_test->SetLineWidth(1);

    Oqq_SIGNAL_test->SetMarkerStyle(kFullCircle);
    Oqq_SIGNAL_test->SetLineColor(kBlue);
    Oqq_SIGNAL_test->SetMarkerColor(kBlue);
    Oqq_SIGNAL_test->SetLineWidth(1);

    OBB_BKG_test->SetMarkerStyle(kFullCircle);
    OBB_BKG_test->SetLineColor(kRed);
    OBB_BKG_test->SetMarkerColor(kRed);
    OBB_BKG_test->SetLineWidth(1);

    Oqq_BKG_test->SetMarkerStyle(kFullCircle);
    Oqq_BKG_test->SetLineColor(kRed);
    Oqq_BKG_test->SetMarkerColor(kRed);
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

    gStyle->SetOptStat(0);

    TCanvas* c_temp = new TCanvas("c", "", 600, 600); c_temp->cd(); gPad->SetLogy();
    double OBB_BKG_train_max = OBB_BKG_train->GetMaximum(); double OBB_SIGNAL_train_max = OBB_SIGNAL_train->GetMaximum();
    if(OBB_BKG_train_max > OBB_SIGNAL_train_max) OBB_BKG_train->SetMaximum(1.05 * OBB_BKG_train_max);
    else OBB_BKG_train->SetMaximum(1.05 * OBB_SIGNAL_train_max);
    OBB_BKG_train->SetMinimum(0.03);
    OBB_BKG_train->Draw("Hist"); OBB_BKG_train->GetYaxis()->SetTitleOffset(1.3);
    OBB_SIGNAL_train->Draw("HistSAME");
    OBB_BKG_test->Draw("AP SAME"); OBB_SIGNAL_test->Draw("AP SAME");
    TLegend* legend = gPad->BuildLegend(0.8, 0.9, 0.3, 0.6); legend->SetFillStyle(0); legend->SetLineWidth(0);
    c_temp->SaveAs("OBB_Plot.png");

    TCanvas* c_temp_2 = new TCanvas("c2", "", 600, 600); c_temp_2->cd(); gPad->SetLogy();
    double Oqq_BKG_train_max = Oqq_BKG_train->GetMaximum(); double Oqq_SIGNAL_train_max = Oqq_SIGNAL_train->GetMaximum();
    if(Oqq_BKG_train_max > Oqq_SIGNAL_train_max) Oqq_BKG_train->SetMaximum(1.05 * Oqq_BKG_train_max);
    else Oqq_BKG_train->SetMaximum(1.05 * Oqq_SIGNAL_train_max);
    Oqq_BKG_train->Draw("Hist"); Oqq_BKG_train->GetYaxis()->SetTitleOffset(1.3);
    Oqq_SIGNAL_train->Draw("HistSAME");
    Oqq_BKG_test->Draw("AP SAME"); Oqq_SIGNAL_test->Draw("AP SAME");
    TLegend* legend_2 = gPad->BuildLegend(0.8, 0.9, 0.3, 0.6); legend_2->SetFillStyle(0); legend_2->SetLineWidth(0);
    c_temp_2->SaveAs("Oqq_Plot.png");

}
