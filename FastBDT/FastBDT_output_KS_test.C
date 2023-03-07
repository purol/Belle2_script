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

// SIGNAL MC sample number befor skimming
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

// BKG MC sample number (0.8/ab for BB, 1.0/ab for qq)
# define N_CHG_test 32042497.0
# define N_MIX_test 24693710.0
# define N_UUBAR_test 94447089.0
# define N_DDBAR_test 22664556.0
# define N_SSBAR_test 19244661.0
# define N_CHARM_test 107541168.0
# define N_CHG_train 48052238.0
# define N_MIX_train 37030486.0
# define N_UUBAR_train 141671998.0
# define N_DDBAR_train 34114182.0
# define N_SSBAR_train 28859338.0
# define N_CHARM_train 161280679.0

// new scale factor for BKG MC sample with additional 1/ab (364.436 - 2.763 = 361.673/fb), until LS1
# define Scale_CHG_train ((N_BB_LS1* (BR_BpBp / (BR_BpBp + BR_B0B0))) / (0.8 * N_BpBp_1invab * (N_CHG_train / (N_CHG_train + N_CHG_test)) + N_BpBp_1invab))
# define Scale_MIX_train ((N_BB_LS1* (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (0.8 * N_B0B0_1invab * (N_MIX_train / (N_MIX_train + N_MIX_test)) + N_B0B0_1invab))
# define Scale_UUBAR_train (0.361673/((N_UUBAR_train/(N_UUBAR_train + N_UUBAR_test))*1.0))
# define Scale_DDBAR_train (0.361673/((N_DDBAR_train/(N_DDBAR_train + N_DDBAR_test))*1.0))
# define Scale_SSBAR_train (0.361673/((N_SSBAR_train/(N_SSBAR_train + N_SSBAR_test))*1.0))
# define Scale_CHARM_train (0.361673/((N_CHARM_train/(N_CHARM_train + N_CHARM_test))*1.0))
# define Scale_CHG_test ((N_BB_LS1* (BR_BpBp / (BR_BpBp + BR_B0B0))) / (0.8 * N_BpBp_1invab * (N_CHG_test / (N_CHG_train + N_CHG_test)) + N_BpBp_1invab))
# define Scale_MIX_test ((N_BB_LS1* (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (0.8 * N_B0B0_1invab * (N_MIX_test / (N_MIX_train + N_MIX_test)) + N_B0B0_1invab))
# define Scale_UUBAR_test (0.361673/((N_UUBAR_test/(N_UUBAR_train + N_UUBAR_test))*1.0))
# define Scale_DDBAR_test (0.361673/((N_DDBAR_test/(N_DDBAR_train + N_DDBAR_test))*1.0))
# define Scale_SSBAR_test (0.361673/((N_SSBAR_test/(N_SSBAR_train + N_SSBAR_test))*1.0))
# define Scale_CHARM_test (0.361673/((N_CHARM_test/(N_CHARM_train + N_CHARM_test))*1.0))

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
TH1F* OBB_BKG_train = new TH1F("FastBDT_{2} BKG train", ";FastBDT_{2};", 40, 0, 1.0);
TH1F* Oqq_BKG_train = new TH1F("FastBDT_{1} BKG train", ";FastBDT_{1};", 40, 0, 1.0);
TH2F* OBB_Oqq_BKG_test = new TH2F("OBB_Oqq_BKG_test", ";O_{BB};O{qq}", 40, 0, 1.0, 40, 0, 1.0);
TH1F* OBB_BKG_test = new TH1F("FastBDT_{2} BKG test", ";FastBDT_{2};", 40, 0, 1.0);
TH1F* Oqq_BKG_test = new TH1F("FastBDT_{1} BKG test", ";FastBDT_{1};", 40, 0, 1.0);

TH2F* OBB_Oqq_SIGNAL_train = new TH2F("OBB_Oqq_SIGNAL_train", ";O_{BB};O{qq}", 40, 0, 1.0, 40, 0, 1.0);
TH1F* OBB_SIGNAL_train = new TH1F("FastBDT_{2} SIGNAL train", ";FastBDT_{2};", 40, 0, 1.0);
TH1F* Oqq_SIGNAL_train = new TH1F("FastBDT_{1} SIGNAL train", ";FastBDT_{1};", 40, 0, 1.0);
TH2F* OBB_Oqq_SIGNAL_test = new TH2F("OBB_Oqq_SIGNAL_test", ";O_{BB};O{qq}", 40, 0, 1.0, 40, 0, 1.0);
TH1F* OBB_SIGNAL_test = new TH1F("FastBDT_{2} SIGNAL test", ";FastBDT_{2};", 40, 0, 1.0);
TH1F* Oqq_SIGNAL_test = new TH1F("FastBDT_{1} SIGNAL test", ";FastBDT_{1};", 40, 0, 1.0);

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
            
            if(Mx < 1.1) continue;
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
    const char* train_dirname_Knunu = "./temp_v000_train/B2Knunu";
    const char* train_dirname_Kstarnunu = "./temp_v000_train/B2Kstarnunu";
    const char* train_dirname_Xsununu = "./temp_v000_train/B2Xsnunu";
    const char* train_dirname_K0nunu = "./temp_v000_train/B02K0nunu";
    const char* train_dirname_K0starnunu = "./temp_v000_train/B02K0starnunu";
    const char* train_dirname_Xsdnunu = "./temp_v000_train/B02Xsnunu";
    LetsFill(train_dirname_Knunu, OBB_Oqq_SIGNAL_train, OBB_SIGNAL_train, Oqq_SIGNAL_train, Scale_Kplus_train);
    LetsFill(train_dirname_Kstarnunu, OBB_Oqq_SIGNAL_train, OBB_SIGNAL_train, Oqq_SIGNAL_train, Scale_Kplusstar_train);
    LetsFill(train_dirname_Xsununu, OBB_Oqq_SIGNAL_train, OBB_SIGNAL_train, Oqq_SIGNAL_train, Scale_Xsu_nonresonant_train);
    LetsFill(train_dirname_K0nunu, OBB_Oqq_SIGNAL_train, OBB_SIGNAL_train, Oqq_SIGNAL_train, Scale_K0_train);
    LetsFill(train_dirname_K0starnunu, OBB_Oqq_SIGNAL_train, OBB_SIGNAL_train, Oqq_SIGNAL_train, Scale_K0star_train);
    LetsFill(train_dirname_Xsdnunu, OBB_Oqq_SIGNAL_train, OBB_SIGNAL_train, Oqq_SIGNAL_train, Scale_Xsd_nonresonant_train);

    const char* train_dirname_CHG = "./temp_v000_train/CHG";
    const char* train_dirname_MIX = "./temp_v000_train/MIX";
    const char* train_dirname_UUBAR = "./temp_v000_train/UUBAR";
    const char* train_dirname_DDBAR = "./temp_v000_train/DDBAR";
    const char* train_dirname_SSBAR = "./temp_v000_train/SSBAR";
    const char* train_dirname_CHARM = "./temp_v000_train/CHARM";
    LetsFill(train_dirname_CHG, OBB_Oqq_BKG_train, OBB_BKG_train, Oqq_BKG_train, Scale_CHG_train);
    LetsFill(train_dirname_MIX, OBB_Oqq_BKG_train, OBB_BKG_train, Oqq_BKG_train, Scale_MIX_train);
    LetsFill(train_dirname_UUBAR, OBB_Oqq_BKG_train, OBB_BKG_train, Oqq_BKG_train, Scale_UUBAR_train);
    LetsFill(train_dirname_DDBAR, OBB_Oqq_BKG_train, OBB_BKG_train, Oqq_BKG_train, Scale_DDBAR_train);
    LetsFill(train_dirname_SSBAR, OBB_Oqq_BKG_train, OBB_BKG_train, Oqq_BKG_train, Scale_SSBAR_train);
    LetsFill(train_dirname_CHARM, OBB_Oqq_BKG_train, OBB_BKG_train, Oqq_BKG_train, Scale_CHARM_train);

    const char* test_dirname_Knunu = "./temp_v000/B2Knunu";
    const char* test_dirname_Kstarnunu = "./temp_v000/B2Kstarnunu";
    const char* test_dirname_Xsununu = "./temp_v000/B2Xsnunu";
    const char* test_dirname_K0nunu = "./temp_v000/B02K0nunu";
    const char* test_dirname_K0starnunu = "./temp_v000/B02K0starnunu";
    const char* test_dirname_Xsdnunu = "./temp_v000/B02Xsnunu";
    LetsFill(test_dirname_Knunu, OBB_Oqq_SIGNAL_test, OBB_SIGNAL_test, Oqq_SIGNAL_test, Scale_Kplus_test);
    LetsFill(test_dirname_Kstarnunu, OBB_Oqq_SIGNAL_test, OBB_SIGNAL_test, Oqq_SIGNAL_test, Scale_Kplusstar_test);
    LetsFill(test_dirname_Xsununu, OBB_Oqq_SIGNAL_test, OBB_SIGNAL_test, Oqq_SIGNAL_test, Scale_Xsu_nonresonant_test);
    LetsFill(test_dirname_K0nunu, OBB_Oqq_SIGNAL_test, OBB_SIGNAL_test, Oqq_SIGNAL_test, Scale_K0_test);
    LetsFill(test_dirname_K0starnunu, OBB_Oqq_SIGNAL_test, OBB_SIGNAL_test, Oqq_SIGNAL_test, Scale_K0star_test);
    LetsFill(test_dirname_Xsdnunu, OBB_Oqq_SIGNAL_test, OBB_SIGNAL_test, Oqq_SIGNAL_test, Scale_Xsd_nonresonant_test);

    const char* test_dirname_CHG = "./temp_v000/CHG";
    const char* test_dirname_MIX = "./temp_v000/MIX";
    const char* test_dirname_UUBAR = "./temp_v000/UUBAR";
    const char* test_dirname_DDBAR = "./temp_v000/DDBAR";
    const char* test_dirname_SSBAR = "./temp_v000/SSBAR";
    const char* test_dirname_CHARM = "./temp_v000/CHARM";
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
    OBB_BKG_train->Draw("Hist"); OBB_SIGNAL_train->Draw("HistSAME");
    OBB_BKG_test->Draw("AP SAME"); OBB_SIGNAL_test->Draw("AP SAME");
    gPad->BuildLegend(0.9, 0.9, 0.6, 0.6);
    c_temp->SaveAs("OBB_Plot.png");

    TCanvas* c_temp_2 = new TCanvas("c2", "", 600, 600); c_temp_2->cd(); gPad->SetLogy();
    double Oqq_BKG_train_max = Oqq_BKG_train->GetMaximum(); double Oqq_SIGNAL_train_max = Oqq_SIGNAL_train->GetMaximum();
    if(Oqq_BKG_train_max > Oqq_SIGNAL_train_max) Oqq_BKG_train->SetMaximum(1.05 * Oqq_BKG_train_max);
    else Oqq_BKG_train->SetMaximum(1.05 * Oqq_SIGNAL_train_max);
    Oqq_BKG_train->Draw("Hist"); Oqq_SIGNAL_train->Draw("HistSAME");
    Oqq_BKG_test->Draw("AP SAME"); Oqq_SIGNAL_test->Draw("AP SAME");
    gPad->BuildLegend(0.9, 0.9, 0.6, 0.6);
    c_temp_2->SaveAs("Oqq_Plot.png");

}
