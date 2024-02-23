
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

# define Nvar 33
# define DvetoNvar 4

char* MCTYPE;

using std::string;

std::vector<std::string> var_names;

void FillVariables(const char* filename, std::vector<float> input_vars[Nvar], std::vector<bool>* IsSignal, std::vector<float>* weight, bool tempissignal, float weight_N = 1.0, std::string CorrectionType = "otherwise") {
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */
    TFile* input_file = new TFile(filename, "read");

    TTree* tree_data = (TTree*)input_file->Get("data");

    double Vars[Nvar];
    int flag;

    double Dc_chiProb; // 0.0
    double Dc_pvalue_med;
    double Dc_pvalue_std; // 0.0
    double Dc_dr; // -1.0
    double Dc_dz; // -100.0
    double Dc_M; // 0.0
    double D0_chiProb;
    double D0_pvalue_med;
    double D0_pvalue_std;
    double D0_dr;
    double D0_dz;
    double D0_M;

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double Mxs = -1;
    double Pcms = -1;

    double invM = -1;

    tree_data->SetBranchAddress("Bsig_cosTBTO", &Vars[0]);
    tree_data->SetBranchAddress("Bsig_KSFWVariables_hso01", &Vars[1]);
    tree_data->SetBranchAddress("Bsig_KSFWVariables_hso04", &Vars[2]);
    tree_data->SetBranchAddress("Bsig_thrustBm", &Vars[3]);
    tree_data->SetBranchAddress("Bsig_useCMSFrame_p", &Vars[4]);
    tree_data->SetBranchAddress("Btag_CleoConeCS_1", &Vars[5]);
    tree_data->SetBranchAddress("Btag_CleoConeCS_2", &Vars[6]);
    tree_data->SetBranchAddress("Btag_CleoConeCS_3", &Vars[7]);
    tree_data->SetBranchAddress("Btag_cosTBTO", &Vars[8]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hoo1", &Vars[9]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hoo2", &Vars[10]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hoo3", &Vars[11]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hoo4", &Vars[12]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hso02", &Vars[13]);
    tree_data->SetBranchAddress("Btag_KSFWVariables_hso24", &Vars[14]);
    tree_data->SetBranchAddress("Btag_useCMSFrame_theta", &Vars[15]);
    tree_data->SetBranchAddress("extraInfo__boEeclv200__bc", &Vars[16]);
    tree_data->SetBranchAddress("extraInfo__boNgammav200__bc", &Vars[17]);
    tree_data->SetBranchAddress("foxWolframR1", &Vars[18]);
    tree_data->SetBranchAddress("foxWolframR3", &Vars[19]);
    tree_data->SetBranchAddress("foxWolframR4", &Vars[20]);
    tree_data->SetBranchAddress("harmonicMomentThrust1", &Vars[21]);
    tree_data->SetBranchAddress("harmonicMomentThrust2", &Vars[22]);
    tree_data->SetBranchAddress("missingEnergyOfEventCMS", &Vars[23]);
    tree_data->SetBranchAddress("missingMomentumOfEvent", &Vars[24]);
    tree_data->SetBranchAddress("missingMomentumOfEvent_theta", &Vars[25]);
    tree_data->SetBranchAddress("nRemainingTracksInEvent", &Vars[26]);
    tree_data->SetBranchAddress("roePTheta__bocleanMask__bc", &Vars[27]);
    tree_data->SetBranchAddress("useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc", &Vars[28]);

    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dc_pValue_med", &Dc_pvalue_med);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dc_pValue_std", &Dc_pvalue_std);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_chiProb", &Dc_chiProb);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_dz", &Dc_dz);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_M", &Dc_M);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_D0_pValue_med", &D0_pvalue_med);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_D0_pValue_std", &D0_pvalue_std);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_D0simpleveto_chiProb", &D0_chiProb);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_D0simpleveto_dz", &D0_dz);
    tree_data->SetBranchAddress("Bsig_daughter_0_extraInfo_D0simpleveto_M", &D0_M);

    tree_data->SetBranchAddress("flag", &flag);

    if (tempissignal) {
        tree_data->SetBranchAddress("nParticlesInList__boB__pl__clKcharge_total__bc", &Decay[0]);
        tree_data->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &Decay[1]);
        tree_data->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &Decay[2]);
        tree_data->SetBranchAddress("nParticlesInList__boXsu__clMCcomb__bc", &Decay[3]);
        tree_data->SetBranchAddress("nParticlesInList__boXsu__clMCch1__bc", &Decay[4]);
        tree_data->SetBranchAddress("nParticlesInList__boXsu__clMCch2__bc", &Decay[5]);
        tree_data->SetBranchAddress("nParticlesInList__boXsu__clMCch3__bc", &Decay[6]);
        tree_data->SetBranchAddress("nParticlesInList__boXsu__clMCch4__bc", &Decay[7]);
        tree_data->SetBranchAddress("nParticlesInList__boXsu__clMCch5__bc", &Decay[8]);
        tree_data->SetBranchAddress("nParticlesInList__boXsu__clMCch6__bc", &Decay[9]);
        tree_data->SetBranchAddress("nParticlesInList__boXsu__clMCch7__bc", &Decay[10]);
        tree_data->SetBranchAddress("nParticlesInList__boXsu__clMCch8__bc", &Decay[11]);
        tree_data->SetBranchAddress("nParticlesInList__boXsu__clMCch9__bc", &Decay[12]);
        tree_data->SetBranchAddress("nParticlesInList__boXsu__clMCch10__bc", &Decay[13]);
        tree_data->SetBranchAddress("nParticlesInList__boXsu__clMCch11__bc", &Decay[14]);
        tree_data->SetBranchAddress("nParticlesInList__boXsu__clMCch12__bc", &Decay[15]);
        tree_data->SetBranchAddress("nParticlesInList__boXsu__clMCch13__bc", &Decay[16]);
        tree_data->SetBranchAddress("nParticlesInList__boXsu__clMCch14__bc", &Decay[17]);
        tree_data->SetBranchAddress("nParticlesInList__boXsu__clMCch15__bc", &Decay[18]);
        tree_data->SetBranchAddress("nParticlesInList__boB0__clKneutral_total__bc", &Decay[19]);
        tree_data->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &Decay[20]);
        tree_data->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &Decay[21]);
        tree_data->SetBranchAddress("nParticlesInList__boXsd__clMCcomb__bc", &Decay[22]);
        tree_data->SetBranchAddress("nParticlesInList__boXsd__clMCch16__bc", &Decay[23]);
        tree_data->SetBranchAddress("nParticlesInList__boXsd__clMCch17__bc", &Decay[24]);
        tree_data->SetBranchAddress("nParticlesInList__boXsd__clMCch18__bc", &Decay[25]);
        tree_data->SetBranchAddress("nParticlesInList__boXsd__clMCch19__bc", &Decay[26]);
        tree_data->SetBranchAddress("nParticlesInList__boXsd__clMCch20__bc", &Decay[27]);
        tree_data->SetBranchAddress("nParticlesInList__boXsd__clMCch21__bc", &Decay[28]);
        tree_data->SetBranchAddress("nParticlesInList__boXsd__clMCch22__bc", &Decay[29]);
        tree_data->SetBranchAddress("nParticlesInList__boXsd__clMCch23__bc", &Decay[30]);
        tree_data->SetBranchAddress("nParticlesInList__boXsd__clMCch24__bc", &Decay[31]);
        tree_data->SetBranchAddress("nParticlesInList__boXsd__clMCch25__bc", &Decay[32]);
        tree_data->SetBranchAddress("nParticlesInList__boXsd__clMCch26__bc", &Decay[33]);
        tree_data->SetBranchAddress("nParticlesInList__boXsd__clMCch27__bc", &Decay[34]);
        tree_data->SetBranchAddress("nParticlesInList__boXsd__clMCch28__bc", &Decay[35]);
        tree_data->SetBranchAddress("nParticlesInList__boXsd__clMCch29__bc", &Decay[36]);
        tree_data->SetBranchAddress("nParticlesInList__boXsd__clMCch30__bc", &Decay[37]);

        tree_data->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

        tree_data->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
        tree_data->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);
    }
    tree_data->SetBranchAddress("Bsig_M", &Mxs);

    int Nevt = 0;
    //printf("%lld entries...\n", tree_data->GetEntries());
    for (unsigned int j = 0; j < tree_data->GetEntries(); j++) { // Fill
        tree_data->GetEntry(j);
        if (tempissignal == true && (Decay[0] > 0.5 || Decay[1] > 0.5 || Decay[2] > 0.5 || Decay[19] > 0.5 || Decay[20] > 0.5 || Decay[21] > 0.5) && Mxs > 1.1) continue;
        else if (tempissignal == true && (Decay[0] < 0.5 && Decay[1] < 0.5 && Decay[2] < 0.5 && Decay[19] < 0.5 && Decay[20] < 0.5 && Decay[21] < 0.5) && Mxs < 1.1) continue;
        Nevt++;

        for (unsigned int k = 0; k < Nvar - DvetoNvar; k++) input_vars[k].push_back((float)Vars[k]);

        if (Dc_chiProb > -0.5) {
            input_vars[Nvar - DvetoNvar + 0].push_back((float)Dc_pvalue_std);
            input_vars[Nvar - DvetoNvar + 1].push_back((float)Dc_M);
        }
        else {
            input_vars[Nvar - DvetoNvar + 0].push_back((float)0.0);
            input_vars[Nvar - DvetoNvar + 1].push_back((float)0.0);
        }
        if (D0_chiProb > -0.5) {
            input_vars[Nvar - DvetoNvar + 2].push_back((float)D0_pvalue_std);
            input_vars[Nvar - DvetoNvar + 3].push_back((float)D0_M);
        }
        else {
            input_vars[Nvar - DvetoNvar + 2].push_back((float)0.0);
            input_vars[Nvar - DvetoNvar + 3].push_back((float)0.0);
        }

        IsSignal->push_back(tempissignal);

        double FF_reweight = 1.0;
        double Fragmentation_reweight = 1.0;
        if (CorrectionType == "B2Knunu") FF_reweight = FF_reweight * corrector.GetCorrectionFactor(invM * invM, "Bplus");
        else if (CorrectionType == "B02K0nunu") FF_reweight = FF_reweight * corrector.GetCorrectionFactor(invM * invM, "Bzero");
        else if (CorrectionType == "B2Xsnunu") Fragmentation_reweight = Fragmentation_reweight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);
        else if (CorrectionType == "B02Xsnunu") Fragmentation_reweight = Fragmentation_reweight * corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_B0_MC, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MCTYPE);

        weight->push_back(weight_N * FF_reweight * Fragmentation_reweight);

    }

    input_file->Close();
    //printf("==> Total %d events survive...\n", Nevt);
}

double PrintMaximumFOM(const FastBDT::Classifier& classifier, std::vector<std::vector<float>> InputVariables, std::vector<bool> IsSignal, std::vector<float> weight) {
    const int step = 100;
    double FOM_max = -1;

    for (int i = 0; i < step; i++) {
        float value = ((float)i) / ((float)step);
        double NBKG = 0;
        double NSIG = 0;

        for (unsigned int i = 0; i < IsSignal.size(); ++i) {
            std::vector<float> temp;
            for (int j = 0; j < Nvar; j++) temp.push_back(InputVariables.at(j).at(i));
            float p = classifier.predict(temp);
            if (p > value) {
                if(IsSignal[i]) NSIG = NSIG + weight[i];
                else NBKG = NBKG + weight[i];
            }
        }

        double FOM = NSIG / std::sqrt(NBKG + NSIG);
        if (FOM > FOM_max) FOM_max = FOM;

    }

    return FOM_max;
}

double PrintAUC(const FastBDT::Classifier& classifier, std::vector<std::vector<float>> InputVariables, std::vector<bool> IsSignal, std::vector<float> weight) {
    const int step = 100;
    double AUC = 0;
    double NBKG_total = 0;
    double NSIG_total = 0;
    std::vector<double> TPRs;
    std::vector<double> FPRs;

    for (unsigned int i = 0; i < IsSignal.size(); ++i) {
        if (IsSignal[i]) NSIG_total = NSIG_total + weight[i];
        else NBKG_total = NBKG_total + weight[i];
    }

    for (int i = 0; i < step; i++) {
        float value = ((float)i) / ((float)step);
        double NBKG = 0;
        double NSIG = 0;

        for (unsigned int i = 0; i < IsSignal.size(); ++i) {
            std::vector<float> temp;
            for (int j = 0; j < Nvar; j++) temp.push_back(InputVariables.at(j).at(i));
            float p = classifier.predict(temp);
            if (p >= value) {
                if (IsSignal[i]) NSIG = NSIG + weight[i];
                else NBKG = NBKG + weight[i];
            }
        }

        double TPR = NSIG / NSIG_total;
        double FPR = NBKG / NBKG_total;

        TPRs.push_back(TPR);
        FPRs.push_back(FPR);
    }

    for (unsigned int i = 0; i < TPRs.size(); ++i) {
        if ( i != TPRs.size() - 1) {
            double del_FPR = FPRs.at(i) - FPRs.at(i + 1);
            double avg_TPR = (TPRs.at(i) + TPRs.at(i + 1)) / 2.0;
            AUC = AUC + del_FPR * avg_TPR;
        }
        else {
            double del_FPR = FPRs.at(i) - 0.0;
            double avg_TPR = (TPRs.at(i) + 0.0) / 2.0;
            AUC = AUC + del_FPR * avg_TPR;
        }
    }

    return AUC;
}

double PrintAVG(const FastBDT::Classifier& classifier, std::vector<std::vector<float>> InputVariables, std::vector<bool> IsSignal, std::vector<float> weight, bool SelectSignal) {
    double NBKG_total = 0;
    double NSIG_total = 0;

    double NBKG_AVG = 0;
    double NSIG_AVG = 0;

    for (unsigned int i = 0; i < IsSignal.size(); ++i) {
        if (IsSignal[i]) NSIG_total = NSIG_total + weight[i];
        else NBKG_total = NBKG_total + weight[i];
    }

    for (unsigned int i = 0; i < IsSignal.size(); ++i) {
        std::vector<float> temp;
        for (int j = 0; j < Nvar; j++) temp.push_back(InputVariables.at(j).at(i));
        float p = classifier.predict(temp);

        if (IsSignal[i]) NSIG_AVG = NSIG_AVG + p * weight[i];
        else NBKG_AVG = NBKG_AVG + p * weight[i];

    }

    NSIG_AVG = NSIG_AVG / NSIG_total;
    NBKG_AVG = NBKG_AVG / NBKG_total;

    if (SelectSignal) return NSIG_AVG;
    else return NBKG_AVG;
}

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
        std::vector<string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B2Knunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 1.0, "B2Knunu");
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B2Kstarnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 1.0, "otherwise");
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B2Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 1.0, "B2Xsnunu");
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B02K0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 1.0, "B02K0nunu");
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B02Kstar0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 1.0, "otherwise");
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B02Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, true, 1.0, "B02Xsnunu");
        }
    }
    {
        std::vector<string> names;
        load_files(CHG_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHG_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, ObtainWeight("CHG", MCTYPE, "train", std::string("")), "otherwise");
        }
    }
    {
        std::vector<string> names;
        load_files(MIX_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((MIX_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, ObtainWeight("MIX", MCTYPE, "train", std::string("")), "otherwise");
        }
    }
    {
        std::vector<string> names;
        load_files(UUBAR_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((UUBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, ObtainWeight("UUBAR", MCTYPE, "train", std::string("")), "otherwise");
        }
    }
    {
        std::vector<string> names;
        load_files(DDBAR_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((DDBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, ObtainWeight("DDBAR", MCTYPE, "train", std::string("")), "otherwise");
        }
    }
    {
        std::vector<string> names;
        load_files(SSBAR_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SSBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, ObtainWeight("SSBAR", MCTYPE, "train", std::string("")), "otherwise");
        }
    }
    {
        std::vector<string> names;
        load_files(CHARM_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHARM_input_train + std::string("/") + names.at(i)).c_str(), input_vars, &IsSignal, &weight, false, ObtainWeight("CHARM", MCTYPE, "train", std::string("")), "otherwise");
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
        std::vector<string> names;
        load_files(SIGNAL_input_test.c_str(), &names, "B2Knunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, ObtainWeight("SIGNAL", MCTYPE, "test", std::string("B2Knunu")), "B2Knunu");
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test.c_str(), &names, "B2Kstarnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, ObtainWeight("SIGNAL", MCTYPE, "test", std::string("B2Kstarnunu")), "otherwise");
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test.c_str(), &names, "B2Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, ObtainWeight("SIGNAL", MCTYPE, "test", std::string("B2Xsnunu")), "B2Xsnunu");
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test.c_str(), &names, "B02K0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, ObtainWeight("SIGNAL", MCTYPE, "test", std::string("B02K0nunu")), "B02K0nunu");
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test.c_str(), &names, "B02Kstar0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, ObtainWeight("SIGNAL", MCTYPE, "test", std::string("B02Kstar0nunu")), "otherwise");
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_test.c_str(), &names, "B02Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SIGNAL_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, true, ObtainWeight("SIGNAL", MCTYPE, "test", std::string("B02Xsnunu")), "B02Xsnunu");
        }
    }
    {
        std::vector<string> names;
        load_files(CHG_input_test.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((CHG_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, ObtainWeight("CHG", MCTYPE, "test", std::string("")), "otherwise");
        }
    }
    {
        std::vector<string> names;
        load_files(MIX_input_test.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((MIX_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, ObtainWeight("MIX", MCTYPE, "test", std::string("")), "otherwise");
        }
    }
    {
        std::vector<string> names;
        load_files(UUBAR_input_test.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((UUBAR_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, ObtainWeight("UUBAR", MCTYPE, "test", std::string("")), "otherwise");
        }
    }
    {
        std::vector<string> names;
        load_files(DDBAR_input_test.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((DDBAR_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, ObtainWeight("DDBAR", MCTYPE, "test", std::string("")), "otherwise");
        }
    }
    {
        std::vector<string> names;
        load_files(SSBAR_input_test.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((SSBAR_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, ObtainWeight("SSBAR", MCTYPE, "test", std::string("")), "otherwise");
        }
    }
    {
        std::vector<string> names;
        load_files(CHARM_input_test.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            // take even only
            if (i % 2 == 1) continue;
            FillVariables((CHARM_input_test + std::string("/") + names.at(i)).c_str(), input_vars2, &IsSignal2, &weight2, false, ObtainWeight("CHARM", MCTYPE, "test", std::string("")), "otherwise");
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
        std::vector<string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B2Knunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, ObtainWeight("SIGNAL", MCTYPE, "train", std::string("B2Knunu")), "B2Knunu");
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B2Kstarnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, ObtainWeight("SIGNAL", MCTYPE, "train", std::string("B2Kstarnunu")), "otherwise");
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B2Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, ObtainWeight("SIGNAL", MCTYPE, "train", std::string("B2Xsnunu")), "B2Xsnunu");
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B02K0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, ObtainWeight("SIGNAL", MCTYPE, "train", std::string("B02K0nunu")), "B02K0nunu");
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B02Kstar0nunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, ObtainWeight("SIGNAL", MCTYPE, "train", std::string("B02Kstar0nunu")), "otherwise");
        }
    }
    {
        std::vector<string> names;
        load_files(SIGNAL_input_train.c_str(), &names, "B02Xsnunu");
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SIGNAL_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, true, ObtainWeight("SIGNAL", MCTYPE, "train", std::string("B02Xsnunu")), "B02Xsnunu");
        }
    }
    {
        std::vector<string> names;
        load_files(CHG_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHG_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, ObtainWeight("CHG", MCTYPE, "train", std::string("")), "otherwise");
        }
    }
    {
        std::vector<string> names;
        load_files(MIX_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((MIX_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, ObtainWeight("MIX", MCTYPE, "train", std::string("")), "otherwise");
        }
    }
    {
        std::vector<string> names;
        load_files(UUBAR_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((UUBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, ObtainWeight("UUBAR", MCTYPE, "train", std::string("")), "otherwise");
        }
    }
    {
        std::vector<string> names;
        load_files(DDBAR_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((DDBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, ObtainWeight("DDBAR", MCTYPE, "train", std::string("")), "otherwise");
        }
    }
    {
        std::vector<string> names;
        load_files(SSBAR_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((SSBAR_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, ObtainWeight("SSBAR", MCTYPE, "train", std::string("")), "otherwise");
        }
    }
    {
        std::vector<string> names;
        load_files(CHARM_input_train.c_str(), &names);
        for (unsigned int i = 0; i < names.size(); ++i) {
            FillVariables((CHARM_input_train + std::string("/") + names.at(i)).c_str(), input_vars3, &IsSignal3, &weight3, false, ObtainWeight("CHARM", MCTYPE, "train", std::string("")), "otherwise");
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
