
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

#include "Classifier.h"

#include "correctors.h"
#include "constants.h"
#include "base.h"
#include "ObtainWeight.h"
#include "MyFBDT.h"

char* MCTYPE;

int main(int argc, char* argv[])
{
    // grid search
    // unsigned int nTrees[5] = { 100, 500, 1000, 1500, 2000 };  default is 100
    // unsigned int depth[3] = { 2, 3, 4 };  default is 3 
    // double shrinkage[4] = { 0.05, 0.1, 0.15, 0.2 };  default is 0.1
    // double subsample[5] = { 0.3, 0.4, 0.5, 0.6, 0.7 };  default is 0.5
    // unsigned int binning[4] = { 6, 7, 8, 9 };  default is 2^8 bins per feature

    /*
    * argv[1]: nTrees
    * argv[2]: depth
    * argv[3]: shrinkage path
    * argv[4]: subsample type
    * argv[5]: binning
    * argv[6]: version name (ex. Aqua, Kokoro, Satori, ...)
    * argb[7]: dirname (ex. v000, v001, ...)
    * argv[8]: MC version: {MC15ri|MC15rd}
    */

    unsigned int nTrees = (unsigned int)atoi(argv[1]);
    unsigned int depth = (unsigned int)atoi(argv[2]);
    double shrinkage = atof(argv[3]);
    double subsample = atof(argv[4]);
    unsigned int binning_num = (unsigned int)atoi(argv[5]);

    MCTYPE = argv[8];

    // set classifier option
    FastBDT::Classifier classifier;
    classifier.SetNTrees(nTrees);
    classifier.SetDepth(depth);
    classifier.SetShrinkage(shrinkage);
    classifier.SetSubsample(subsample);
    std::vector<unsigned int> binning(Nvar, binning_num); classifier.SetBinning(binning);



    // input file
    std::string SIGNAL_input_train = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/SIGNAL_analysis/train_" + std::string(argv[7]) + "/final_output_data";
    std::string CHG_input_train = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/CHG_analysis/train_" + std::string(argv[7]) + "/final_output_data";
    std::string MIX_input_train = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/MIX_analysis/train_" + std::string(argv[7]) + "/final_output_data";
    std::string UUBAR_input_train = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/UUBAR_analysis/train_" + std::string(argv[7]) + "/final_output_data";
    std::string DDBAR_input_train = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/DDBAR_analysis/train_" + std::string(argv[7]) + "/final_output_data";
    std::string SSBAR_input_train = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/SSBAR_analysis/train_" + std::string(argv[7]) + "/final_output_data";
    std::string CHARM_input_train = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/CHARM_analysis/train_" + std::string(argv[7]) + "/final_output_data";

    std::string SIGNAL_input_test = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/SIGNAL_analysis/test_" + std::string(argv[7]) + "/final_output_data";
    std::string CHG_input_test = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/CHG_analysis/test_" + std::string(argv[7]) + "/final_output_data";
    std::string MIX_input_test = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/MIX_analysis/test_" + std::string(argv[7]) + "/final_output_data";
    std::string UUBAR_input_test = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/UUBAR_analysis/test_" + std::string(argv[7]) + "/final_output_data";
    std::string DDBAR_input_test = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/DDBAR_analysis/test_" + std::string(argv[7]) + "/final_output_data";
    std::string SSBAR_input_test = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/SSBAR_analysis/test_" + std::string(argv[7]) + "/final_output_data";
    std::string CHARM_input_test = "/home/belle2/junewoo/storage_b1/bsub/Analysis/" + std::string(argv[6]) + "/CHARM_analysis/test_" + std::string(argv[7]) + "/final_output_data";



    // define input of the classifier
    std::vector<std::vector<float>> InputVariables;
    std::vector<bool> IsSignal;
    std::vector<float> weight;

    // define input variables
    std::vector<float> input_vars[Nvar];

    {
        std::vector<std::string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B2Knunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, MCTYPE, 1.0, "B2Knunu");
        }
    }
    {
        std::vector<std::string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B2Kstarnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, MCTYPE, 1.0, "otherwise");
        }
    }
    {
        std::vector<std::string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B2Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, MCTYPE, 1.0, "B2Xsnunu");
        }
    }
    {
        std::vector<std::string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B02K0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, MCTYPE, 1.0, "B02K0nunu");
        }
    }
    {
        std::vector<std::string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B02Kstar0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, MCTYPE, 1.0, "otherwise");
        }
    }
    {
        std::vector<std::string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B02Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, MCTYPE, 1.0, "B02Xsnunu");
        }
    }
    {
        std::vector<std::string> names;
        load_files(CHG_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHG_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, MCTYPE, ObtainWeight("CHG", MCTYPE, "train", std::string("")), "otherwise");
        }
    }
    {
        std::vector<std::string> names;
        load_files(MIX_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((MIX_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, MCTYPE, ObtainWeight("MIX", MCTYPE, "train", std::string("")), "otherwise");
        }
    }
    {
        std::vector<std::string> names;
        load_files(UUBAR_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((UUBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, MCTYPE, ObtainWeight("UUBAR", MCTYPE, "train", std::string("")), "otherwise");
        }
    }
    {
        std::vector<std::string> names;
        load_files(DDBAR_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((DDBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, MCTYPE, ObtainWeight("DDBAR", MCTYPE, "train", std::string("")), "otherwise");
        }
    }
    {
        std::vector<std::string> names;
        load_files(SSBAR_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SSBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, MCTYPE, MCTYPE, ObtainWeight("SSBAR", MCTYPE, "train", std::string("")), "otherwise");
        }
    }
    {
        std::vector<std::string> names;
        load_files(CHARM_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHARM_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, MCTYPE, ObtainWeight("CHARM", MCTYPE, "train", std::string("")), "otherwise");
        }
    }

    // convert to double vector
    for (int i = 0; i < Nvar; i++) {
        InputVariables.push_back(input_vars[i]);
    }



    // fit
    classifier.fit(InputVariables, IsSignal, weight);

    // clear vector to save memory
    for (unsigned int i = 0; i < InputVariables.size(); ++i) std::vector<float>().swap(InputVariables.at(i));
    std::vector<std::vector<float>>().swap(InputVariables);
    std::vector<bool>().swap(IsSignal);
    std::vector<float>().swap(weight);
    for (int i = 0; i < Nvar; i++) std::vector<float>().swap(input_vars[i]);



    // test sample
    std::vector<std::vector<float>> InputVariables2;
    std::vector<bool> IsSignal2;
    std::vector<float> weight2;

    std::vector<float> input_vars2[Nvar];

    {
        std::vector<std::string> names;
        load_files(SIGNAL_input_test.c_str(), &names, "B2Knunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, MCTYPE, ObtainWeight("SIGNAL", MCTYPE, "test", std::string("B2Knunu")), "B2Knunu");
        }
    }
    {
        std::vector<std::string> names;
        load_files(SIGNAL_input_test.c_str(), &names, "B2Kstarnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, MCTYPE, ObtainWeight("SIGNAL", MCTYPE, "test", std::string("B2Kstarnunu")), "otherwise");
        }
    }
    {
        std::vector<std::string> names;
        load_files(SIGNAL_input_test.c_str(), &names, "B2Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, MCTYPE, ObtainWeight("SIGNAL", MCTYPE, "test", std::string("B2Xsnunu")), "B2Xsnunu");
        }
    }
    {
        std::vector<std::string> names;
        load_files(SIGNAL_input_test.c_str(), &names, "B02K0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, MCTYPE, ObtainWeight("SIGNAL", MCTYPE, "test", std::string("B02K0nunu")), "B02K0nunu");
        }
    }
    {
        std::vector<std::string> names;
        load_files(SIGNAL_input_test.c_str(), &names, "B02Kstar0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, MCTYPE, ObtainWeight("SIGNAL", MCTYPE, "test", std::string("B02Kstar0nunu")), "otherwise");
        }
    }
    {
        std::vector<std::string> names;
        load_files(SIGNAL_input_test.c_str(), &names, "B02Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, MCTYPE, ObtainWeight("SIGNAL", MCTYPE, "test", std::string("B02Xsnunu")), "B02Xsnunu");
        }
    }
    {
        std::vector<std::string> names;
        load_files(CHG_input_test.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((CHG_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, MCTYPE, ObtainWeight("CHG", MCTYPE, "test", std::string("")), "otherwise");
        }
    }
    {
        std::vector<std::string> names;
        load_files(MIX_input_test.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((MIX_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, MCTYPE, ObtainWeight("MIX", MCTYPE, "test", std::string("")), "otherwise");
        }
    }
    {
        std::vector<std::string> names;
        load_files(UUBAR_input_test.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((UUBAR_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, MCTYPE, ObtainWeight("UUBAR", MCTYPE, "test", std::string("")), "otherwise");
        }
    }
    {
        std::vector<std::string> names;
        load_files(DDBAR_input_test.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((DDBAR_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, MCTYPE, ObtainWeight("DDBAR", MCTYPE, "test", std::string("")), "otherwise");
        }
    }
    {
        std::vector<std::string> names;
        load_files(SSBAR_input_test.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SSBAR_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, MCTYPE, ObtainWeight("SSBAR", MCTYPE, "test", std::string("")), "otherwise");
        }
    }
    {
        std::vector<std::string> names;
        load_files(CHARM_input_test.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((CHARM_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, MCTYPE, ObtainWeight("CHARM", MCTYPE, "test", std::string("")), "otherwise");
        }
    }

    // convert to double vector
    for (int i = 0; i < Nvar; i++) {
        InputVariables2.push_back(input_vars2[i]);
    }



    // get FOM for testing sample
    // double test_FOM = PrintMaximumFOM(classifier, InputVariables2, IsSignal2, weight2);

    // get AUC for testing sample
    double test_AUC = PrintAUC(classifier, InputVariables2, IsSignal2, weight2);

    // clear vector to save memory
    for (unsigned int i = 0; i < InputVariables2.size(); ++i) std::vector<float>().swap(InputVariables2.at(i));
    std::vector<std::vector<float>>().swap(InputVariables2);
    std::vector<bool>().swap(IsSignal2);
    std::vector<float>().swap(weight2);
    for (int i = 0; i < Nvar; i++) std::vector<float>().swap(input_vars2[i]);



    // train sample with proper weight
    std::vector<std::vector<float>> InputVariables3;
    std::vector<bool> IsSignal3;
    std::vector<float> weight3;

    std::vector<float> input_vars3[Nvar];

    {
        std::vector<std::string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B2Knunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, MCTYPE, ObtainWeight("SIGNAL", MCTYPE, "train", std::string("B2Knunu")), "B2Knunu");
        }
    }
    {
        std::vector<std::string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B2Kstarnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, MCTYPE, ObtainWeight("SIGNAL", MCTYPE, "train", std::string("B2Kstarnunu")), "otherwise");
        }
    }
    {
        std::vector<std::string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B2Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, MCTYPE, ObtainWeight("SIGNAL", MCTYPE, "train", std::string("B2Xsnunu")), "B2Xsnunu");
        }
    }
    {
        std::vector<std::string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B02K0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, MCTYPE, ObtainWeight("SIGNAL", MCTYPE, "train", std::string("B02K0nunu")), "B02K0nunu");
        }
    }
    {
        std::vector<std::string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B02Kstar0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, MCTYPE, ObtainWeight("SIGNAL", MCTYPE, "train", std::string("B02Kstar0nunu")), "otherwise");
        }
    }
    {
        std::vector<std::string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B02Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, MCTYPE, ObtainWeight("SIGNAL", MCTYPE, "train", std::string("B02Xsnunu")), "B02Xsnunu");
        }
    }
    {
        std::vector<std::string> names;
        load_files(CHG_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHG_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, MCTYPE, ObtainWeight("CHG", MCTYPE, "train", std::string("")), "otherwise");
        }
    }
    {
        std::vector<std::string> names;
        load_files(MIX_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((MIX_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, MCTYPE, ObtainWeight("MIX", MCTYPE, "train", std::string("")), "otherwise");
        }
    }
    {
        std::vector<std::string> names;
        load_files(UUBAR_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((UUBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, MCTYPE, ObtainWeight("UUBAR", MCTYPE, "train", std::string("")), "otherwise");
        }
    }
    {
        std::vector<std::string> names;
        load_files(DDBAR_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((DDBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, MCTYPE, ObtainWeight("DDBAR", MCTYPE, "train", std::string("")), "otherwise");
        }
    }
    {
        std::vector<std::string> names;
        load_files(SSBAR_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SSBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, MCTYPE, ObtainWeight("SSBAR", MCTYPE, "train", std::string("")), "otherwise");
        }
    }
    {
        std::vector<std::string> names;
        load_files(CHARM_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHARM_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, MCTYPE, ObtainWeight("CHARM", MCTYPE, "train", std::string("")), "otherwise");
        }
    }

    // convert to double vector
    for (int i = 0; i < Nvar; i++) {
        InputVariables3.push_back(input_vars3[i]);
    }



    // get FOM for training sample
    // double train_FOM = PrintMaximumFOM(classifier, InputVariables3, IsSignal3, weight3);

    // get AUC for training sample
    double train_AUC = PrintAUC(classifier, InputVariables3, IsSignal3, weight3);

    // clear vector to save memory
    for (unsigned int i = 0; i < InputVariables3.size(); ++i) std::vector<float>().swap(InputVariables3.at(i));
    std::vector<std::vector<float>>().swap(InputVariables3);
    std::vector<bool>().swap(IsSignal3);
    std::vector<float>().swap(weight3);
    for (int i = 0; i < Nvar; i++) std::vector<float>().swap(input_vars3[i]);



    printf("%u_%u_%lf_%lf_%u %lf %lf\n", nTrees, depth, shrinkage, subsample, binning_num, train_AUC, test_AUC);

    FILE* fp;
    fp = fopen(("/home/belle2/junewoo/storage_b1/bsub/Analysis/GridSearch/out/Result_" + std::string(argv[1]) + "_" + std::string(argv[2]) + "_" + std::string(argv[3]) + "_" + std::string(argv[4]) + "_" + std::string(argv[5])).c_str(), "w");
    fprintf(fp, "%u_%u_%lf_%lf_%u %lf %lf\n", nTrees, depth, shrinkage, subsample, binning_num, train_AUC, test_AUC);
    fclose(fp);



    // save model
    std::fstream out_stream(("/home/belle2/junewoo/storage_b1/bsub/Analysis/GridSearch/out/classifier_" + std::string(argv[1]) + "_" + std::string(argv[2]) + "_" + std::string(argv[3]) + "_" + std::string(argv[4]) + "_" + std::string(argv[5])+".weightfile").c_str(), std::ios_base::out | std::ios_base::trunc);
    out_stream << classifier << std::endl;
    out_stream.close();

    return 0;
}
