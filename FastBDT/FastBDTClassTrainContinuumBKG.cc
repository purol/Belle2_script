
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

#include "Classifier.h"

# define N_Needed_info 37
# define N_event_info 15
# define N_Upsilon_info 45
# define N_Bsig_info 60
# define N_Btag_info 7
# define N_decay 38 // five decay mode + others

# define Nvar 15

void FillVariables(const char* filename, std::vector<float> input_vars[Nvar], std::vector<bool>* IsSignal, std::vector<float>* weight) {
    TFile* input_file = new TFile(filename, "read");

    TTree* tree_data = (TTree*)input_file->Get("data");

    double Vars[Nvar];
    double flag;

    tree_data->SetBranchAddress("MsquaredBsig_op1", &Vars[0]);
    tree_data->SetBranchAddress("nRemainingTracksInEvent", &Vars[1]);
    tree_data->SetBranchAddress("missingMomentumOfEvent_theta", &Vars[2]);
    tree_data->SetBranchAddress("Btag_extraInfo_SignalProbability", &Vars[3]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hoo1", &Vars[4]);
    tree_data->SetBranchAddress("harmonicMomentThrust1", &Vars[5]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hso00", &Vars[6]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hso22", &Vars[7]);
    tree_data->SetBranchAddress("Btag_useCMSFrame_theta", &Vars[8]);
    tree_data->SetBranchAddress("Btag_deltaE", &Vars[9]);
    tree_data->SetBranchAddress("Btag_cosTBTO", &Vars[10]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hso24", &Vars[11]);
    tree_data->SetBranchAddress("Btag_chiProb", &Vars[12]);
    tree_data->SetBranchAddress("Btag_thrustBm", &Vars[13]);
    tree_data->SetBranchAddress("cleoConeThrust0", &Vars[14]);
    tree_data->SetBranchAddress("flag", &flag);

    printf("%lld entries...\n", tree_data->GetEntries());
    for (unsigned int j = 0; j < tree_data->GetEntries(); j++) { // Fill
        tree_data->GetEntry(j);

        for (unsigned int k = 0; k < Nvar; k++) input_vars[k].push_back((float)Vars[k]);
        if (flag > 0.5) IsSignal->push_back(true);
        else IsSignal->push_back(false);

        weight->push_back(1.0f);

    }

    input_file->Close();
}

float GetScore(const FastBDT::Classifier& classifier, std::vector<std::vector<float>> InputVariables, std::vector<bool> IsSignal) {
    float sum = 0;
    for (unsigned int i = 0; i < IsSignal.size(); ++i) {
        std::vector<float> temp;
        for (int j = 0; j < Nvar; j++) temp.push_back(InputVariables.at(j).at(i));
        float p = classifier.predict(temp);
        sum += (static_cast<int>(IsSignal[i]) - p) * (static_cast<int>(IsSignal[i]) - p);
    }
    return sum / IsSignal.size();
}

int main()
{
    // define classifier and set options
    FastBDT::Classifier classifier;
    classifier.SetNTrees(300);

    // define input of the classifier
    std::vector<std::vector<float>> InputVariables;
    std::vector<bool> IsSignal;
    std::vector<float> weight;

    // define input variables
    std::vector<float> input_vars[Nvar];

    // input file
    const char* SIGNAL_input = "/media/sf_virtualbox_folder/20211109/For_TMVA_SIGNAL/train/Mxs_large/final_output_merge_Mxs_larger_SIGNAL_train_data.root";
    const char* UUBAR_input = "/media/sf_virtualbox_folder/20211109/For_TMVA_BKG/train/Mxs_large/final_output_merge_Mxs_larger_UUBAR_train_data.root";
    const char* DDBAR_input = "/media/sf_virtualbox_folder/20211109/For_TMVA_BKG/train/Mxs_large/final_output_merge_Mxs_larger_DDBAR_train_data.root";
    const char* SSBAR_input = "/media/sf_virtualbox_folder/20211109/For_TMVA_BKG/train/Mxs_large/final_output_merge_Mxs_larger_SSBAR_train_data.root";
    const char* CHARM_input = "/media/sf_virtualbox_folder/20211109/For_TMVA_BKG/train/Mxs_large/final_output_merge_Mxs_larger_CHARM_train_data.root";
    const char* SIGNAL_input_test = "/media/sf_virtualbox_folder/20211109/For_TMVA_SIGNAL/test/Mxs_large/final_output_merge_Mxs_larger_SIGNAL_test_data.root";
    const char* UUBAR_input_test = "/media/sf_virtualbox_folder/20211109/For_TMVA_BKG/test/Mxs_large/final_output_merge_Mxs_larger_UUBAR_test_data.root";
    const char* DDBAR_input_test = "/media/sf_virtualbox_folder/20211109/For_TMVA_BKG/test/Mxs_large/final_output_merge_Mxs_larger_DDBAR_test_data.root";
    const char* SSBAR_input_test = "/media/sf_virtualbox_folder/20211109/For_TMVA_BKG/test/Mxs_large/final_output_merge_Mxs_larger_SSBAR_test_data.root";
    const char* CHARM_input_test = "/media/sf_virtualbox_folder/20211109/For_TMVA_BKG/test/Mxs_large/final_output_merge_Mxs_larger_CHARM_test_data.root";

    FillVariables(SIGNAL_input, input_vars, &IsSignal, &weight);
    FillVariables(UUBAR_input, input_vars, &IsSignal, &weight);
    FillVariables(DDBAR_input, input_vars, &IsSignal, &weight);
    FillVariables(SSBAR_input, input_vars, &IsSignal, &weight);
    FillVariables(CHARM_input, input_vars, &IsSignal, &weight);

    // convert to double vector
    for (int i = 0; i < Nvar; i++) {
        InputVariables.push_back(input_vars[i]);
    }

    // fit
    classifier.fit(InputVariables, IsSignal, weight);

    // print result
    std::cout << "Train Score " << GetScore(classifier, InputVariables, IsSignal) << std::endl;

    std::fstream out_stream("unittest.weightfile", std::ios_base::out | std::ios_base::trunc);
    out_stream << classifier << std::endl;
    out_stream.close();

    classifier.Print();

    std::map<unsigned int, double> rank;
    rank = classifier.GetVariableRanking();
    printf("Variable importance:\n");
        for (auto iter = rank.begin(); iter != rank.end(); iter++)
        {
            std::cout << "[" << iter->first << ", " << iter->second << "]" << " ";
        }
    printf("\n\n");


    // print result of testing sample
    std::fstream in_stream("unittest.weightfile", std::ios_base::in);
    FastBDT::Classifier classifier2(in_stream);

    std::vector<std::vector<float>> InputVariables2;
    std::vector<bool> IsSignal2;
    std::vector<float> weight2;

    std::vector<float> input_vars2[Nvar];

    FillVariables(SIGNAL_input_test, input_vars2, &IsSignal2, &weight2);
    FillVariables(UUBAR_input_test, input_vars2, &IsSignal2, &weight2);
    FillVariables(DDBAR_input_test, input_vars2, &IsSignal2, &weight2);
    FillVariables(SSBAR_input_test, input_vars2, &IsSignal2, &weight2);
    FillVariables(CHARM_input_test, input_vars2, &IsSignal2, &weight2);

    // convert to double vector
    for (int i = 0; i < Nvar; i++) {
        InputVariables2.push_back(input_vars2[i]);
    }

    std::cout << "Test Score " << GetScore(classifier2, InputVariables2, IsSignal2) << std::endl;

    classifier2.Print();

    return 0;
}