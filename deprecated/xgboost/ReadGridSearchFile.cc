
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
# include <vector>
#include <fstream>

#include <TMath.h>
#include <TColor.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TString.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THStack.h>
#include <TPaveText.h>
#include <TKey.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <numeric>

#include <algorithm>

#include "Classifier.h"

# define N_Needed_info 37
# define N_event_info 15
# define N_Upsilon_info 47
# define N_Bsig_info 78
# define N_Btag_info 7
# define N_decay 38 // five decay mode + others

# define Nvar 34
# define DvetoNvar 6

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

using std::string;

std::vector<std::string> var_names;

typedef struct data {
    double eta = 0;
    double gamma = 0;
    int max_depth = 0;
    double min_child_weight = 0;
    double subsample = 0;
    double colsample_bytree = 0;
    int Niter = 0;
    double train_AUC = 0;
    double test_AUC = 0;
} Data;

bool data_sorter(Data const& lhs, Data const& rhs) {
    return lhs.test_AUC < rhs.test_AUC;
}

int ReadGridSearchFile()
{
    const char* fname = "";
    const int Nlist = 1500;
    std::vector<Data> Datas;

    FILE* fp;
    fp = fopen(fname, "r");

    for (int i = 0; i < Nlist; i++) {
        double eta = 0;
        double gamma = 0;
        int max_depth = 0;
        double min_child_weight = 0;
        double subsample = 0;
        double colsample_bytree = 0;
        int Niter = 0;
        double train_AUC = 0;
        double test_AUC = 0;
        int dummy = 0;

        fscanf(fp,"%lf_%lf_%d_%lf_%lf_%lf_%d %lf %lf %d\n", eta, gamma, max_depth, min_child_weight, subsample, colsample_bytree, Niter, train_AUC, test_AUC, dummy);
        Data temp_data = { eta, gamma, max_depth, min_child_weight, subsample, colsample_bytree, Niter, train_AUC, test_AUC };
        Datas.push_back(temp_data);
    }

    fclose(fp);

    std::sort(Datas.begin(), Datas.end(), &data_sorter);

    double Rank[Nlist];
    for (int i = 0; i < Nlist; i++) Rank[i] = i + 1;
    double train_AUCs[Nlist];
    for (int i = 0; i < Nlist; i++) train_AUCs[i] = Datas.at(i).train_AUC;
    double test_AUCs[Nlist];
    for (int i = 0; i < Nlist; i++) test_AUCs[i] = Datas.at(i).test_AUC;

    for (int i = 0; i < Nlist; i++) printf("%lf_%lf_%d_%lf_%lf_%lf_%d %lf %lf\n", Datas.at(i).eta, Datas.at(i).gamma, Datas.at(i).max_depth, Datas.at(i).min_child_weight, Datas.at(i).subsample, Datas.at(i).colsample_bytree, Datas.at(i).Niter, Datas.at(i).train_AUC, Datas.at(i).test_AUC);

    TGraph* gr_train = new TGraph(Nlist, Rank, train_AUCs);
    TGraph* gr_test = new TGraph(Nlist, Rank, test_AUCs);

    gr_train->SetMarkerStyle(8); gr_train->SetMarkerSize(0.8);
    gr_test->SetMarkerStyle(8); gr_test->SetMarkerSize(0.8);

    gr_train->SetMarkerColor(kRed + 1);
    gr_test->SetMarkerColor(kBlue + 1);

    gr_train->SetTitle(";AUC rank;AUC");

    TCanvas* c = new TCanvas("c1", "AUC", 200, 10, 600, 600);

    gr_train->Draw("AP");
    gr_test->Draw("P");

    TLegend* legend = new TLegend(0.15, 0.8, 0.35, 0.9); legend->SetFillStyle(0); legend->SetLineWidth(0);
    legend->AddEntry(gr_train, "train", "P"); legend->AddEntry(gr_test, "test", "P");
    legend->Draw();
    c->SaveAs("AUC.png");

    return 0;
}
