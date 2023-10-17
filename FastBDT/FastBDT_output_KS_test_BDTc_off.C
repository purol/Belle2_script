#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"

#include <string>
# include <vector>

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

// BKG MC sample number (2.8/ab for BB, 1.0/ab for qq)
# define N_CHG_test 135437652.0
# define N_MIX_test 104176964.0
# define N_UUBAR_test 106465323.0
# define N_DDBAR_test 25669347.0
# define N_SSBAR_test 20666484.0
# define N_CHARM_test 116381155.0
# define N_CHG_train 156250364.0
# define N_MIX_train 120553880.0
# define N_UUBAR_train 159685667.0
# define N_DDBAR_train 38563469.0
# define N_SSBAR_train 31004176.0
# define N_CHARM_train 174628991.0

// new scale factor for BKG MC sample with additional 1/ab (364.436 - 2.763 = 361.673/fb), until LS1
# define Scale_CHG_train ((N_BB_LS1* (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab * (N_CHG_train / (N_CHG_train + N_CHG_test)) ))
# define Scale_MIX_train ((N_BB_LS1* (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab * (N_MIX_train / (N_MIX_train + N_MIX_test)) ))
# define Scale_UUBAR_train (0.361673/((N_UUBAR_train/(N_UUBAR_train + N_UUBAR_test))*1.0))
# define Scale_DDBAR_train (0.361673/((N_DDBAR_train/(N_DDBAR_train + N_DDBAR_test))*1.0))
# define Scale_SSBAR_train (0.361673/((N_SSBAR_train/(N_SSBAR_train + N_SSBAR_test))*1.0))
# define Scale_CHARM_train (0.361673/((N_CHARM_train/(N_CHARM_train + N_CHARM_test))*1.0))
# define Scale_CHG_test ((N_BB_LS1* (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab * (N_CHG_test / (N_CHG_train + N_CHG_test)) ))
# define Scale_MIX_test ((N_BB_LS1* (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab * (N_MIX_test / (N_MIX_train + N_MIX_test)) ))
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

TH1F* FBDTc_MC_train = new TH1F("MC train", ";BDT_{c};", 40, 0.0, 1.0);
TH1F* FBDTc_MC_test = new TH1F("MC test", ";BDT_{c};", 40, 0.0, 1.0);

TH1F* FBDTc_MC_train_correction = new TH1F("MC train with correction", ";BDT_{c};", 40, 0.0, 1.0);
TH1F* FBDTc_MC_test_correction = new TH1F("MC test with correction", ";BDT_{c};", 40, 0.0, 1.0);

TH1F* FBDTc_data_train = new TH1F("data train", ";BDT_{c};", 40, 0.0, 1.0);
TH1F* FBDTc_data_test = new TH1F("data test", ";BDT_{c};", 40, 0.0, 1.0);

typedef struct _Nevt {
    double NevtwithoutCorrection;
    double NevtwithCorrection;
    double BDTcSum;
} Nevt;

double BDTcToWeight(double BDTc) {

    if ( BDTc > (5.0 / 6.0) ) return std::sqrt(5.0);
    else return std::sqrt( BDTc / (1.0 - BDTc) );

    //return 2*BDTc;

}

double LetsFill(const char* filename, TH1F* FBDTc_hist, Nevt* nevt, double weight_var = 1.0) {
    float FBDT_var = 0;
    double NormFactor = 0;

    TFile* input_file = new TFile(filename, "read");

    TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
    tree_upsilon->SetBranchAddress("MVA_Continuum", &FBDT_var);

    for (unsigned int i = 0; i < tree_upsilon->GetEntries(); i++) {
        tree_upsilon->GetEntry(i);
        FBDTc_hist->Fill(FBDT_var, weight_var);

        NormFactor = NormFactor + BDTcToWeight(FBDT_var) * weight_var;

        nevt->NevtwithoutCorrection = nevt->NevtwithoutCorrection + weight_var;
        nevt->NevtwithCorrection = nevt->NevtwithCorrection + BDTcToWeight(FBDT_var) * weight_var;
        nevt->BDTcSum = nevt->BDTcSum + FBDT_var * weight_var;

    }
    input_file->Close();

    return NormFactor;
}

double LetsFill(const char* filename, TH1F* FBDTc_hist, double weight_var = 1.0) {
    float FBDT_var = 0;
    double NormFactor = 0;

    TFile* input_file = new TFile(filename, "read");

    TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
    tree_upsilon->SetBranchAddress("MVA_Continuum", &FBDT_var);

    for (unsigned int i = 0; i < tree_upsilon->GetEntries(); i++) {
        tree_upsilon->GetEntry(i);
        FBDTc_hist->Fill(FBDT_var, weight_var);

        NormFactor = NormFactor + BDTcToWeight(FBDT_var) * weight_var;

    }
    input_file->Close();

    return NormFactor;
}

double LetsFillwithCorrection(const char* filename, TH1F* FBDTc_hist, double weight_var = 1.0, double NormFactor = 1.0) {
    float FBDT_var = 0;

    TFile* input_file = new TFile(filename, "read");

    TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
    tree_upsilon->SetBranchAddress("MVA_Continuum", &FBDT_var);

    for (unsigned int i = 0; i < tree_upsilon->GetEntries(); i++) {
        tree_upsilon->GetEntry(i);

        FBDTc_hist->Fill(FBDT_var, BDTcToWeight(FBDT_var) * weight_var / NormFactor);

    }
    input_file->Close();

    return NormFactor;
}

void FastBDT_output_KS_test_BDTc_off()
{
    Nevt nevt_train = { 0.0, 0.0, 0.0 };
    Nevt nevt_test = { 0.0, 0.0, 0.0 };

    Nevt data_all = { 0.0, 0.0, 0.0 };

    // dirnames
    const char* side_data = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_LS_data_off/temp_v000/SIGNAL";

    const char* side_MC_UUBAR_train = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_LS_MC_off/temp_v000_train/UUBAR";
    const char* side_MC_DDBAR_train = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_LS_MC_off/temp_v000_train/DDBAR";
    const char* side_MC_SSBAR_train = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_LS_MC_off/temp_v000_train/SSBAR";
    const char* side_MC_CHARM_train = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_LS_MC_off/temp_v000_train/CHARM";

    const char* side_MC_UUBAR_test = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_LS_MC_off/temp_v000/UUBAR";
    const char* side_MC_DDBAR_test = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_LS_MC_off/temp_v000/DDBAR";
    const char* side_MC_SSBAR_test = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_LS_MC_off/temp_v000/SSBAR";
    const char* side_MC_CHARM_test = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_LS_MC_off/temp_v000/CHARM";


    // train sample
    {
        std::vector<string> names;
        load_files(side_data, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 == 0) continue;
            LetsFill((side_data + std::string("/") + names.at(i)).c_str(), FBDTc_data_train, &data_all, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_UUBAR_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            LetsFill((side_MC_UUBAR_train + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train, &nevt_train, Scale_UUBAR_train * (0.0361697 / 0.361673));
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_DDBAR_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            LetsFill((side_MC_DDBAR_train + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train, &nevt_train, Scale_DDBAR_train * (0.0361697 / 0.361673));
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_SSBAR_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            LetsFill((side_MC_SSBAR_train + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train, &nevt_train, Scale_SSBAR_train * (0.0361697 / 0.361673));
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_CHARM_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            LetsFill((side_MC_CHARM_train + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train, &nevt_train, Scale_CHARM_train * (0.0361697 / 0.361673));
        }
    }

    // test sample
    {
        std::vector<string> names;
        load_files(side_data, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 != 0) continue;
            LetsFill((side_data + std::string("/") + names.at(i)).c_str(), FBDTc_data_test, &data_all, 1.0);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_UUBAR_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            LetsFill((side_MC_UUBAR_test + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test, &nevt_test, Scale_UUBAR_test * (0.0061593 / 0.361673));
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_DDBAR_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            LetsFill((side_MC_DDBAR_test + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test, &nevt_test, Scale_DDBAR_test * (0.0061593 / 0.361673));
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_SSBAR_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            LetsFill((side_MC_SSBAR_test + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test, &nevt_test, Scale_SSBAR_test * (0.0061593 / 0.361673));
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_CHARM_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            LetsFill((side_MC_CHARM_test + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test, &nevt_test, Scale_CHARM_test * (0.0061593 / 0.361673));
        }
    }


    double NormFactorForTrain = nevt_train.NevtwithoutCorrection / nevt_train.NevtwithCorrection;
    double NormFactorForTest = nevt_test.NevtwithoutCorrection / nevt_test.NevtwithCorrection;


    // Corrected MC train
    {
        std::vector<string> names;
        load_files(side_MC_UUBAR_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            LetsFillwithCorrection((side_MC_UUBAR_train + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train_correction, Scale_UUBAR_train * (0.0361697 / 0.361673), NormFactorForTrain);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_DDBAR_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            LetsFillwithCorrection((side_MC_DDBAR_train + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train_correction, Scale_DDBAR_train* (0.0361697 / 0.361673), NormFactorForTrain);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_SSBAR_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            LetsFillwithCorrection((side_MC_SSBAR_train + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train_correction, Scale_SSBAR_train* (0.0361697 / 0.361673), NormFactorForTrain);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_CHARM_train, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            LetsFillwithCorrection((side_MC_CHARM_train + std::string("/") + names.at(i)).c_str(), FBDTc_MC_train_correction, Scale_CHARM_train* (0.0361697 / 0.361673), NormFactorForTrain);
        }
    }


    // Corrected MC test
    {
        std::vector<string> names;
        load_files(side_MC_UUBAR_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            LetsFillwithCorrection((side_MC_UUBAR_test + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test_correction, Scale_UUBAR_test* (0.0061593 / 0.361673), NormFactorForTest);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_DDBAR_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            LetsFillwithCorrection((side_MC_DDBAR_test + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test_correction, Scale_DDBAR_test* (0.0061593 / 0.361673), NormFactorForTest);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_SSBAR_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            LetsFillwithCorrection((side_MC_SSBAR_test + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test_correction, Scale_SSBAR_test* (0.0061593 / 0.361673), NormFactorForTest);
        }
    }
    {
        std::vector<string> names;
        load_files(side_MC_CHARM_test, &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            LetsFillwithCorrection((side_MC_CHARM_test + std::string("/") + names.at(i)).c_str(), FBDTc_MC_test_correction, Scale_CHARM_test* (0.0061593 / 0.361673), NormFactorForTest);
        }
    }



    // normalization
    double factor = 1.0;
    FBDTc_MC_train->Scale(factor / FBDTc_MC_train->Integral(), "width");
    FBDTc_MC_test->Scale(factor / FBDTc_MC_test->Integral(), "width");
    FBDTc_MC_train_correction->Scale(factor / FBDTc_MC_train_correction->Integral(), "width");
    FBDTc_MC_test_correction->Scale(factor / FBDTc_MC_test_correction->Integral(), "width");
    FBDTc_data_train->Scale(factor / FBDTc_data_train->Integral(), "width");
    FBDTc_data_test->Scale(factor / FBDTc_data_test->Integral(), "width");


    // set color (MC: kRed, data: kBlue)
    // solid line: test, histogram: traing
    FBDTc_data_test->SetMarkerStyle(kFullCircle);
    FBDTc_data_test->SetLineColor(kBlue);
    FBDTc_data_test->SetMarkerColor(kBlue);
    FBDTc_data_test->SetLineWidth(1);

    FBDTc_MC_test->SetMarkerStyle(kFullCircle);
    FBDTc_MC_test->SetLineColor(kRed);
    FBDTc_MC_test->SetMarkerColor(kRed);
    FBDTc_MC_test->SetLineWidth(1);

    FBDTc_data_train->SetFillStyle(3004);
    FBDTc_data_train->SetLineColor(kBlue);
    FBDTc_data_train->SetFillColor(kBlue);

    FBDTc_MC_train->SetFillStyle(3005);
    FBDTc_MC_train->SetLineColor(kRed);
    FBDTc_MC_train->SetFillColor(kRed);

    FBDTc_MC_test_correction->SetMarkerStyle(kFullCircle);
    FBDTc_MC_test_correction->SetLineColor(kRed);
    FBDTc_MC_test_correction->SetMarkerColor(kRed);
    FBDTc_MC_test_correction->SetLineWidth(1);

    FBDTc_MC_train_correction->SetFillStyle(3005);
    FBDTc_MC_train_correction->SetLineColor(kRed);
    FBDTc_MC_train_correction->SetFillColor(kRed);

    gStyle->SetOptStat(0);

    TCanvas* c_temp = new TCanvas("c", "", 600, 600); c_temp->cd(); // gPad->SetLogy();

    double MC_train_max = FBDTc_MC_train->GetMaximum(); double data_train_max = FBDTc_data_train->GetMaximum();
    if(MC_train_max > data_train_max) FBDTc_MC_train->SetMaximum(1.1 * MC_train_max);
    else FBDTc_MC_train->SetMaximum(1.1 * data_train_max);

    FBDTc_MC_train->Draw("Hist"); FBDTc_data_train->Draw("HistSAME");
    FBDTc_MC_test->Draw("AP SAME"); FBDTc_data_test->Draw("AP SAME");
    gPad->BuildLegend(0.9, 0.9, 0.6, 0.6);
    c_temp->SaveAs("BDTc_BeforeCorrection.png");

    TCanvas* c_temp_2 = new TCanvas("c2", "", 600, 600); c_temp_2->cd(); // gPad->SetLogy();

    if (MC_train_max > data_train_max) FBDTc_MC_train_correction->SetMaximum(1.1 * MC_train_max);
    else FBDTc_MC_train_correction->SetMaximum(1.1 * data_train_max);

    FBDTc_MC_train_correction->Draw("Hist"); FBDTc_data_train->Draw("HistSAME");
    FBDTc_MC_test_correction->Draw("AP SAME"); FBDTc_data_test->Draw("AP SAME");
    gPad->BuildLegend(0.9, 0.9, 0.6, 0.6);
    c_temp_2->SaveAs("BDTc_AfterCorrection.png");

    printf("average BDTc of data: %lf\n", data_all.BDTcSum / data_all.NevtwithoutCorrection);
}
