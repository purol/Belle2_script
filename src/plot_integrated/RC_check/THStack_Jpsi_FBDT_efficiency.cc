// last update: 2021-10-08
// for Belle2 data

/*
when you add new variables:
revise # define N_Needed_info ...
revise typedef struct data
revise void Loader::GetData(TFile* input_file)
revise void Loader::PrintRootFile(std::string output_name)
revise void Loader::PrintSeparateRootFile(std::string output_name)
revise void Loader::ConvertIntoSeparateDataFile(std::string output_name, double flag = 0)
*/

# include <algorithm>
# include <float.h>

#include "constants.h"
#include "base.h"
#include "THStack_plot.h"

int main(int argc, char* argv[]) {
    /*
    * argv[1]: input version (ex. v031)
    * argv[2]: output path
    * argv[3]: output file name
    */

    // dirnames
    const char* Jpsi_MC_SIGNAL_before_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_MC_Jpsi/SIGNAL_analysis/validation_v009/final_output_root_after_MVA_Application";
    const char* Jpsi_MC_CHG_before_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_MC_Jpsi/CHG_analysis/validation_v009/final_output_root_after_MVA_Application";
    const char* Jpsi_MC_MIX_before_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_MC_Jpsi/MIX_analysis/validation_v009/final_output_root_after_MVA_Application";
    const char* Jpsi_MC_UUBAR_before_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_MC_Jpsi/UUBAR_analysis/validation_v009/final_output_root_after_MVA_Application";
    const char* Jpsi_MC_DDBAR_before_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_MC_Jpsi/DDBAR_analysis/validation_v009/final_output_root_after_MVA_Application";
    const char* Jpsi_MC_SSBAR_before_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_MC_Jpsi/SSBAR_analysis/validation_v009/final_output_root_after_MVA_Application";
    const char* Jpsi_MC_CHARM_before_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_MC_Jpsi/CHARM_analysis/validation_v009/final_output_root_after_MVA_Application";
    const char* Jpsi_data_before_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_data_Jpsi/SIGNAL_analysis/validation_v009/final_output_root_after_MVA_Application";

    const char* Jpsi_MC_SIGNAL_after_dirname = strdup(("/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_MC_Jpsi/SIGNAL_analysis/validation_" + std::string(argv[1]) + "/final_output_root_after_MVA_Application_after_cut").c_str());
    const char* Jpsi_MC_CHG_after_dirname = strdup(("/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_MC_Jpsi/CHG_analysis/validation_" + std::string(argv[1]) + "/final_output_root_after_MVA_Application_after_cut").c_str());
    const char* Jpsi_MC_MIX_after_dirname = strdup(("/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_MC_Jpsi/MIX_analysis/validation_" + std::string(argv[1]) + "/final_output_root_after_MVA_Application_after_cut").c_str());
    const char* Jpsi_MC_UUBAR_after_dirname = strdup(("/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_MC_Jpsi/UUBAR_analysis/validation_" + std::string(argv[1]) + "/final_output_root_after_MVA_Application_after_cut").c_str());
    const char* Jpsi_MC_DDBAR_after_dirname = strdup(("/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_MC_Jpsi/DDBAR_analysis/validation_" + std::string(argv[1]) + "/final_output_root_after_MVA_Application_after_cut").c_str());
    const char* Jpsi_MC_SSBAR_after_dirname = strdup(("/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_MC_Jpsi/SSBAR_analysis/validation_" + std::string(argv[1]) + "/final_output_root_after_MVA_Application_after_cut").c_str());
    const char* Jpsi_MC_CHARM_after_dirname = strdup(("/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_MC_Jpsi/CHARM_analysis/validation_" + std::string(argv[1]) + "/final_output_root_after_MVA_Application_after_cut").c_str());
    const char* Jpsi_data_after_dirname = strdup(("/home/belle2/junewoo/storage_ghi/Analysis/KumoiRD_LS_data_Jpsi/SIGNAL_analysis/validation_" + std::string(argv[1]) + "/final_output_root_after_MVA_Application_after_cut").c_str());

    std::vector<std::string> variable_names;
    std::vector<std::string> branch_names;
    variable_names.push_back("MVA_BB"); branch_names.push_back("Upsilon");

    Nvar_num = static_cast<int>(variable_names.size());

    std::vector<double>* Jpsi_MC_values_before = new std::vector<double>[Nvar_num];
    std::vector<double>* Jpsi_MC_values_after = new std::vector<double>[Nvar_num];
    std::vector<int> Jpsi_MC_numbering_before;
    std::vector<int> Jpsi_MC_numbering_after;

    std::vector<double>* Jpsi_data_values_before = new std::vector<double>[Nvar_num];
    std::vector<double>* Jpsi_data_values_after = new std::vector<double>[Nvar_num];

    std::vector<double> weights_before;
    std::vector<double> weights_after;

    LetsFillMC(Jpsi_MC_SIGNAL_before_dirname, variable_names, branch_names, Jpsi_MC_values_before, &Jpsi_MC_numbering_before, &weights_before, "SIGNAL");
    LetsFillMC(Jpsi_MC_CHG_before_dirname, variable_names, branch_names, Jpsi_MC_values_before, &Jpsi_MC_numbering_before, &weights_before, "CHG");
    LetsFillMC(Jpsi_MC_MIX_before_dirname, variable_names, branch_names, Jpsi_MC_values_before, &Jpsi_MC_numbering_before, &weights_before, "MIX");
    LetsFillMC(Jpsi_MC_UUBAR_before_dirname, variable_names, branch_names, Jpsi_MC_values_before, &Jpsi_MC_numbering_before, &weights_before, "UUBAR");
    LetsFillMC(Jpsi_MC_DDBAR_before_dirname, variable_names, branch_names, Jpsi_MC_values_before, &Jpsi_MC_numbering_before, &weights_before, "DDBAR");
    LetsFillMC(Jpsi_MC_SSBAR_before_dirname, variable_names, branch_names, Jpsi_MC_values_before, &Jpsi_MC_numbering_before, &weights_before, "SSBAR");
    LetsFillMC(Jpsi_MC_CHARM_before_dirname, variable_names, branch_names, Jpsi_MC_values_before, &Jpsi_MC_numbering_before, &weights_before, "CHARM");
    LetsFillMC(Jpsi_MC_SIGNAL_after_dirname, variable_names, branch_names, Jpsi_MC_values_after, &Jpsi_MC_numbering_after, &weights_after, "SIGNAL");
    LetsFillMC(Jpsi_MC_CHG_after_dirname, variable_names, branch_names, Jpsi_MC_values_after, &Jpsi_MC_numbering_after, &weights_after, "CHG");
    LetsFillMC(Jpsi_MC_MIX_after_dirname, variable_names, branch_names, Jpsi_MC_values_after, &Jpsi_MC_numbering_after, &weights_after, "MIX");
    LetsFillMC(Jpsi_MC_UUBAR_after_dirname, variable_names, branch_names, Jpsi_MC_values_after, &Jpsi_MC_numbering_after, &weights_after, "UUBAR");
    LetsFillMC(Jpsi_MC_DDBAR_after_dirname, variable_names, branch_names, Jpsi_MC_values_after, &Jpsi_MC_numbering_after, &weights_after, "DDBAR");
    LetsFillMC(Jpsi_MC_SSBAR_after_dirname, variable_names, branch_names, Jpsi_MC_values_after, &Jpsi_MC_numbering_after, &weights_after, "SSBAR");
    LetsFillMC(Jpsi_MC_CHARM_after_dirname, variable_names, branch_names, Jpsi_MC_values_after, &Jpsi_MC_numbering_after, &weights_after, "CHARM");
    LetsFilldata(Jpsi_data_before_dirname, variable_names, branch_names, Jpsi_data_values_before);
    LetsFilldata(Jpsi_data_after_dirname, variable_names, branch_names, Jpsi_data_values_after);

    TH1D* MC_SIGNAL_before_one_bin = new TH1D("MC_SIGNAL_before_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* MC_BKG_before_one_bin = new TH1D("MC_BKG_before_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* data_before_one_bin = new TH1D("data_before_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);

    TH1D* MC_SIGNAL_after_one_bin = new TH1D("MC_SIGNAL_after_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* MC_BKG_after_one_bin = new TH1D("MC_BKG_after_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* data_after_one_bin = new TH1D("data_after_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);

    for (int k = 0; k < (int)Jpsi_MC_numbering_before.size(); k++) {
        if (Jpsi_MC_numbering_before.at(k) >= 0 && Jpsi_MC_numbering_before.at(k) <= 13) { // BKG
            for (int l = 0; l < (int)variable_names.size(); l++) MC_BKG_before_one_bin->Fill(Jpsi_MC_values_before[l].at(k), weights_before.at(k));
        }
        else if (Jpsi_MC_numbering_before.at(k) == 14) { // signal
            for (int l = 0; l < (int)variable_names.size(); l++) MC_SIGNAL_before_one_bin->Fill(Jpsi_MC_values_before[l].at(k), weights_before.at(k));
        }
        else {
            printf("undefined numbering!\n");
            exit(1);
        }
    }

    for (int k = 0; k < (int)Jpsi_data_values_before[0].size(); k++) {
        for (int l = 0; l < (int)variable_names.size(); l++) data_before_one_bin->Fill(Jpsi_data_values_before[l].at(k));
    }

    for (int k = 0; k < (int)Jpsi_MC_numbering_after.size(); k++) {
        if (Jpsi_MC_numbering_after.at(k) >= 0 && Jpsi_MC_numbering_after.at(k) < static_cast<int>(MCsample::SIGNAL)) { // BKG
            for (int l = 0; l < (int)variable_names.size(); l++) MC_BKG_after_one_bin->Fill(Jpsi_MC_values_after[l].at(k), weights_after.at(k));
        }
        else if (Jpsi_MC_numbering_after.at(k) == static_cast<int>(MCsample::SIGNAL)) { // signal
            for (int l = 0; l < (int)variable_names.size(); l++) MC_SIGNAL_after_one_bin->Fill(Jpsi_MC_values_after[l].at(k), weights_after.at(k));
        }
        else {
            printf("undefined numbering!\n");
            exit(1);
        }
    }

    for (int k = 0; k < (int)Jpsi_data_values_after[0].size(); k++) {
        for (int l = 0; l < (int)variable_names.size(); l++) data_after_one_bin->Fill(Jpsi_data_values_after[l].at(k));
    }

    double Nsig_before = MC_SIGNAL_before_one_bin->GetBinContent(1);
    double NBKG_before = MC_BKG_before_one_bin->GetBinContent(1);
    double Ndata_before = data_before_one_bin->GetBinContent(1);
    double Nsig_after = MC_SIGNAL_after_one_bin->GetBinContent(1);
    double NBKG_after = MC_BKG_after_one_bin->GetBinContent(1);
    double Ndata_after = data_after_one_bin->GetBinContent(1);

    double Nsig_before_uncer = MC_SIGNAL_before_one_bin->GetBinError(1);
    double NBKG_before_uncer = MC_BKG_before_one_bin->GetBinError(1);
    double Ndata_before_uncer = data_before_one_bin->GetBinError(1);
    double Nsig_after_uncer = MC_SIGNAL_after_one_bin->GetBinError(1);
    double NBKG_after_uncer = MC_BKG_after_one_bin->GetBinError(1);
    double Ndata_after_uncer = data_after_one_bin->GetBinError(1);

    // Print data-MC discrepancy
    printf("Nsig before: %lf +- %lf\n", Nsig_before, Nsig_before_uncer);
    printf("NBKG before: %lf +- %lf\n", NBKG_before, NBKG_before_uncer);
    printf("Ndata before: %lf +- %lf\n", Ndata_before, Ndata_before_uncer);
    printf("Nsig after: %lf +- %lf\n", Nsig_after, Nsig_after_uncer);
    printf("NBKG after: %lf +- %lf\n", NBKG_after, NBKG_after_uncer);
    printf("Ndata after: %lf +- %lf\n", Ndata_after, Ndata_after_uncer);

    double efficiency_data = (Ndata_after - NBKG_after) / (Ndata_before - NBKG_before);
    double efficiency_MC = Nsig_after / Nsig_before;

    // define TH1D to calculate efficiency
    TH1D* MC_SIGNAL_efficiency = new TH1D("MC_SIGNAL_efficiency", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* data_SIGNAL_before_one_bin = new TH1D("data_SIGNAL_before_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* data_SIGNAL_after_one_bin = new TH1D("data_SIGNAL_after_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* data_SIGNAL_efficiency = new TH1D("data_SIGNAL_efficiency", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* efficiency_ratio = new TH1D("efficiency_ratio", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* MC_after_one_bin = new TH1D("MC_after_one_bin", ";FBDT;Nevt", 1, 0.0, 1.0);
    TH1D* Nevt_ratio = new TH1D("Nevt_ratio", ";FBDT;Nevt", 1, 0.0, 1.0);

    // calculate efficiency!
    MC_SIGNAL_efficiency->Divide(MC_SIGNAL_after_one_bin, MC_SIGNAL_before_one_bin, 1.0, 1.0, "B");
    data_SIGNAL_before_one_bin->Add(data_before_one_bin, MC_BKG_before_one_bin, 1.0, -1.0);
    data_SIGNAL_after_one_bin->Add(data_after_one_bin, MC_BKG_after_one_bin, 1.0, -1.0);
    data_SIGNAL_efficiency->Divide(data_SIGNAL_after_one_bin, data_SIGNAL_before_one_bin, 1.0, 1.0, "B");
    efficiency_ratio->Divide(data_SIGNAL_efficiency, MC_SIGNAL_efficiency);
    MC_after_one_bin->Add(MC_SIGNAL_after_one_bin, MC_BKG_after_one_bin);
    Nevt_ratio->Divide(data_after_one_bin, MC_after_one_bin);

    printf("eps_data for FBDT = %lf +- %lf\n", data_SIGNAL_efficiency->GetBinContent(1), data_SIGNAL_efficiency->GetBinError(1));
    printf("eps_MC for FBDT = %lf +- %lf\n", MC_SIGNAL_efficiency->GetBinContent(1), MC_SIGNAL_efficiency->GetBinError(1));
    printf("eps_data/eps_MC for FBDT = %lf +- %lf\n", efficiency_ratio->GetBinContent(1), efficiency_ratio->GetBinError(1));
    printf("N_{data}/N_{MC}          = %lf +- %lf\n", Nevt_ratio->GetBinContent(1), Nevt_ratio->GetBinError(1));

    FILE* fp;
    fp = fopen((std::string(argv[2]) + "/" + std::string(argv[3])).c_str(), "w");
    fprintf(fp, "%lf %lf\n", efficiency_ratio->GetBinContent(1), efficiency_ratio->GetBinError(1));
    fclose(fp);

    delete[] Jpsi_MC_values_before;
    delete[] Jpsi_MC_values_after;

    delete[] Jpsi_data_values_before;
    delete[] Jpsi_data_values_after;

    delete MC_SIGNAL_before_one_bin;
    delete MC_BKG_before_one_bin;
    delete data_before_one_bin;

    delete MC_SIGNAL_after_one_bin;
    delete MC_BKG_after_one_bin;
    delete data_after_one_bin;
}
