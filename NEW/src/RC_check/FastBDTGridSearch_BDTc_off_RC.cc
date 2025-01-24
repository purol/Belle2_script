
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

#include "constants.h"
#include "base.h"
#include "ObtainWeight.h"
#include "MyFBDT.h"
#include "correctors.h"

Corrector_Knn corrector_Knn;

int main(int argc, char* argv[]) // offres total: 42.329/fb
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

    // set classifier option
    FastBDT::Classifier classifier;
    classifier.SetNTrees(nTrees);
    classifier.SetDepth(depth);
    classifier.SetShrinkage(shrinkage);
    classifier.SetSubsample(subsample);
    std::vector<unsigned int> binning(Nvar, binning_num); classifier.SetBinning(binning);



    // input file
    std::string off_data = ("/home/belle2/junewoo/storage_ghi/Analysis/" + std::string(argv[6]) + "_LS_data_off/SIGNAL_analysis/validation_" + std::string("v009") + "/final_output_data").c_str();
   
    std::string off_MC_UUBAR_train = ("/home/belle2/junewoo/storage_ghi/Analysis/" + std::string(argv[6]) + "_LS_MC_off/UUBAR_analysis/train_" + std::string("v009") + "/final_output_data").c_str();
    std::string off_MC_DDBAR_train = ("/home/belle2/junewoo/storage_ghi/Analysis/" + std::string(argv[6]) + "_LS_MC_off/DDBAR_analysis/train_" + std::string("v009") + "/final_output_data").c_str();
    std::string off_MC_SSBAR_train = ("/home/belle2/junewoo/storage_ghi/Analysis/" + std::string(argv[6]) + "_LS_MC_off/SSBAR_analysis/train_" + std::string("v009") + "/final_output_data").c_str();
    std::string off_MC_CHARM_train = ("/home/belle2/junewoo/storage_ghi/Analysis/" + std::string(argv[6]) + "_LS_MC_off/CHARM_analysis/train_" + std::string("v009") + "/final_output_data").c_str();

    std::string off_MC_UUBAR_test = ("/home/belle2/junewoo/storage_ghi/Analysis/" + std::string(argv[6]) + "_LS_MC_off/UUBAR_analysis/test_" + std::string("v009") + "/final_output_data").c_str();
    std::string off_MC_DDBAR_test = ("/home/belle2/junewoo/storage_ghi/Analysis/" + std::string(argv[6]) + "_LS_MC_off/DDBAR_analysis/test_" + std::string("v009") + "/final_output_data").c_str();
    std::string off_MC_SSBAR_test = ("/home/belle2/junewoo/storage_ghi/Analysis/" + std::string(argv[6]) + "_LS_MC_off/SSBAR_analysis/test_" + std::string("v009") + "/final_output_data").c_str();
    std::string off_MC_CHARM_test = ("/home/belle2/junewoo/storage_ghi/Analysis/" + std::string(argv[6]) + "_LS_MC_off/CHARM_analysis/test_" + std::string("v009") + "/final_output_data").c_str();

    // define input of the classifier
    std::vector<std::vector<float>> InputVariables;
    std::vector<bool> IsSignal;
    std::vector<float> weight;

    // define input variables
    std::vector<float> input_vars[Nvar];

    {
        std::vector<std::string> names;
        load_files(off_data.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 == 0) continue; // offresdata: 7758 of 9012 (36.1697/fb)
            FillVariables_Off((off_data + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 1.0);
        }
    }
    {
        std::vector<std::string> names;
        load_files(off_MC_UUBAR_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables_Off((off_MC_UUBAR_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, ObtainWeight("UUBAR", argv[8], "train", std::string("")) * (0.0364390 / 0.36537) * corrector_Knn.GetCorrectionFactorCancelOutObtainWeight(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, names.at(i), argv[8], false));
        }
    }
    {
        std::vector<std::string> names;
        load_files(off_MC_DDBAR_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables_Off((off_MC_DDBAR_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, ObtainWeight("DDBAR", argv[8], "train", std::string("")) * (0.0364390 / 0.36537) * corrector_Knn.GetCorrectionFactorCancelOutObtainWeight(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, names.at(i), argv[8], false));
        }
    }
    {
        std::vector<std::string> names;
        load_files(off_MC_SSBAR_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables_Off((off_MC_SSBAR_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, ObtainWeight("SSBAR", argv[8], "train", std::string("")) * (0.0364390 / 0.36537) * corrector_Knn.GetCorrectionFactorCancelOutObtainWeight(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, names.at(i), argv[8], false));
        }
    }
    {
        std::vector<std::string> names;
        load_files(off_MC_CHARM_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables_Off((off_MC_CHARM_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, ObtainWeight("CHARM", argv[8], "train", std::string("")) * (0.0364390 / 0.36537) * corrector_Knn.GetCorrectionFactorCancelOutObtainWeight(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, names.at(i), argv[8], false));
        }
    }


    // convert to double vector
    for (int i = 0; i < Nvar; i++) {
        InputVariables.push_back(input_vars[i]);
    }



    // fit
    classifier.fit(InputVariables, IsSignal, weight);


    // get FOM for training sample
    // double train_FOM = PrintMaximumFOM(classifier, InputVariables3, IsSignal3, weight3);

    // get AUC for training sample
    double train_AUC = PrintAUC(classifier, InputVariables, IsSignal, weight);


    // get average value for BKG and SIG
    double train_AVG_SIG = PrintAVG(classifier, InputVariables, IsSignal, weight, true);
    double train_AVG_BKG = PrintAVG(classifier, InputVariables, IsSignal, weight, false);


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
        load_files(off_data.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            if (i % 7 != 0) continue; // offresdata: 1254 of 9012  (6.1593/fb)
            FillVariables_Off((off_data + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, 1.0);
        }
    }
    {
        std::vector<std::string> names;
        load_files(off_MC_UUBAR_test.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables_Off((off_MC_UUBAR_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, ObtainWeight("UUBAR", argv[8], "test", std::string("")) * (0.0058900 / 0.36537) * corrector_Knn.GetCorrectionFactorCancelOutObtainWeight(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, names.at(i), argv[8], false));
        }
    }
    {
        std::vector<std::string> names;
        load_files(off_MC_DDBAR_test.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables_Off((off_MC_DDBAR_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, ObtainWeight("DDBAR", argv[8], "test", std::string("")) * (0.0058900 / 0.36537) * corrector_Knn.GetCorrectionFactorCancelOutObtainWeight(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, names.at(i), argv[8], false));
        }
    }
    {
        std::vector<std::string> names;
        load_files(off_MC_SSBAR_test.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables_Off((off_MC_SSBAR_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, ObtainWeight("SSBAR", argv[8], "test", std::string(""))* (0.0058900 / 0.36537) * corrector_Knn.GetCorrectionFactorCancelOutObtainWeight(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, names.at(i), argv[8], false));
        }
    }
    {
        std::vector<std::string> names;
        load_files(off_MC_CHARM_test.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables_Off((off_MC_CHARM_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, ObtainWeight("CHARM", argv[8], "test", std::string(""))* (0.0058900 / 0.36537) * corrector_Knn.GetCorrectionFactorCancelOutObtainWeight(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, names.at(i), argv[8], false));
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

    // get average value for BKG and SIG
    double AVG_SIG = PrintAVG(classifier, InputVariables2, IsSignal2, weight2, true);
    double AVG_BKG = PrintAVG(classifier, InputVariables2, IsSignal2, weight2, false);

    // clear vector to save memory
    for (unsigned int i = 0; i < InputVariables2.size(); ++i) std::vector<float>().swap(InputVariables2.at(i));
    std::vector<std::vector<float>>().swap(InputVariables2);
    std::vector<bool>().swap(IsSignal2);
    std::vector<float>().swap(weight2);
    for (int i = 0; i < Nvar; i++) std::vector<float>().swap(input_vars2[i]);



    printf("%u_%u_%lf_%lf_%u %lf %lf %lf %lf %lf %lf\n", nTrees, depth, shrinkage, subsample, binning_num, train_AUC, test_AUC, train_AVG_SIG, train_AVG_BKG, AVG_SIG, AVG_BKG);

    FILE* fp;
    fp = fopen(("/home/belle2/junewoo/storage_b1/bsub/Analysis/GridSearch_BDTc/" + std::string(argv[6]) + "_" + std::string(argv[7]) + "/out/Result_" + std::string(argv[1]) + "_" + std::string(argv[2]) + "_" + std::string(argv[3]) + "_" + std::string(argv[4]) + "_" + std::string(argv[5])).c_str(), "w");
    fprintf(fp, "%u_%u_%lf_%lf_%u %lf %lf\n", nTrees, depth, shrinkage, subsample, binning_num, train_AUC, test_AUC);
    fclose(fp);



    // save model
    std::fstream out_stream(("/home/belle2/junewoo/storage_b1/bsub/Analysis/GridSearch_BDTc/" + std::string(argv[6]) + "_" + std::string(argv[7]) + "/out/classifier_" + std::string(argv[1]) + "_" + std::string(argv[2]) + "_" + std::string(argv[3]) + "_" + std::string(argv[4]) + "_" + std::string(argv[5])+".weightfile").c_str(), std::ios_base::out | std::ios_base::trunc);
    out_stream << classifier << std::endl;
    out_stream.close();

    // save model again
    std::fstream out_stream_2(("/home/belle2/junewoo/storage_b1/bsub/Analysis/GridSearch_BDTc/" + std::string(argv[6]) + "_" + std::string(argv[7]) + "/Continuum.weightfile").c_str(), std::ios_base::out | std::ios_base::trunc);
    out_stream_2 << classifier << std::endl;
    out_stream_2.close();

    return 0;
}
