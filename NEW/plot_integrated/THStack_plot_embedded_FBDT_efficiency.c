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

void THStack_plot_embedded_FBDT_efficiency() {

    // dirnames
    const char* Embedded_MC_CHG_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_embedded/SIGNAL_analysis/validation_v005/final_output_root_after_MVA_Application_after_cut/for_plot/CHG";
    const char* Embedded_MC_MIX_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_embedded/SIGNAL_analysis/validation_v005/final_output_root_after_MVA_Application_after_cut/for_plot/MIX";
    const char* Embedded_MC_UUBAR_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_embedded/SIGNAL_analysis/validation_v005/final_output_root_after_MVA_Application_after_cut/for_plot/UUBAR";
    const char* Embedded_MC_DDBAR_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_embedded/SIGNAL_analysis/validation_v005/final_output_root_after_MVA_Application_after_cut/for_plot/DDBAR";
    const char* Embedded_MC_SSBAR_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_embedded/SIGNAL_analysis/validation_v005/final_output_root_after_MVA_Application_after_cut/for_plot/SSBAR";
    const char* Embedded_MC_CHARM_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_embedded/SIGNAL_analysis/validation_v005/final_output_root_after_MVA_Application_after_cut/for_plot/CHARM";
    const char* Embedded_data_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_embedded/SIGNAL_analysis/validation_v005/final_output_root_after_MVA_Application_after_cut/for_plot/data";

    const char* Embedded_MC_CHG_before_FBDT_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_embedded/SIGNAL_analysis/validation_v005/final_output_root_after_MVA_Application/for_plot/CHG";
    const char* Embedded_MC_MIX_before_FBDT_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_embedded/SIGNAL_analysis/validation_v005/final_output_root_after_MVA_Application/for_plot/MIX";
    const char* Embedded_MC_UUBAR_before_FBDT_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_embedded/SIGNAL_analysis/validation_v005/final_output_root_after_MVA_Application/for_plot/UUBAR";
    const char* Embedded_MC_DDBAR_before_FBDT_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_embedded/SIGNAL_analysis/validation_v005/final_output_root_after_MVA_Application/for_plot/DDBAR";
    const char* Embedded_MC_SSBAR_before_FBDT_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_embedded/SIGNAL_analysis/validation_v005/final_output_root_after_MVA_Application/for_plot/SSBAR";
    const char* Embedded_MC_CHARM_before_FBDT_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_embedded/SIGNAL_analysis/validation_v005/final_output_root_after_MVA_Application/for_plot/CHARM";
    const char* Embedded_data_before_FBDT_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_embedded/SIGNAL_analysis/validation_v005/final_output_root_after_MVA_Application/for_plot/data";

    std::vector<std::string> variable_names;
    std::vector<std::string> branch_names;
    variable_names.push_back("MVA_BB"); branch_names.push_back("Upsilon");

    Nvar_num = static_cast<int>(variable_names.size());

    std::vector<double>* Jpsi_MC_values = new std::vector<double>[Nvar_num];
    std::vector<double>* Jpsi_MC_before_FBDT_values = new std::vector<double>[Nvar_num];
    std::vector<double>* Jpsi_data_values = new std::vector<double>[Nvar_num];
    std::vector<double>* Jpsi_data_before_FBDT_values = new std::vector<double>[Nvar_num];

    std::vector<double> weights;
    std::vector<double> weights_before_FBDT;

    LetsFillembeddedMC(Embedded_MC_CHG_dirname, variable_names, branch_names, Jpsi_MC_values, &weights, "CHG");
    LetsFillembeddedMC(Embedded_MC_MIX_dirname, variable_names, branch_names, Jpsi_MC_values, &weights, "MIX");
    LetsFillembeddedMC(Embedded_MC_UUBAR_dirname, variable_names, branch_names, Jpsi_MC_values, &weights, "UUBAR");
    LetsFillembeddedMC(Embedded_MC_DDBAR_dirname, variable_names, branch_names, Jpsi_MC_values, &weights, "DDBAR");
    LetsFillembeddedMC(Embedded_MC_SSBAR_dirname, variable_names, branch_names, Jpsi_MC_values, &weights, "SSBAR");
    LetsFillembeddedMC(Embedded_MC_CHARM_dirname, variable_names, branch_names, Jpsi_MC_values, &weights, "CHARM");
    LetsFilldata(Embedded_data_dirname, variable_names, branch_names, Jpsi_data_values);

    LetsFillembeddedMC(Embedded_MC_CHG_before_FBDT_dirname, variable_names, branch_names, Jpsi_MC_before_FBDT_values, &weights_before_FBDT, "CHG");
    LetsFillembeddedMC(Embedded_MC_MIX_before_FBDT_dirname, variable_names, branch_names, Jpsi_MC_before_FBDT_values, &weights_before_FBDT, "MIX");
    LetsFillembeddedMC(Embedded_MC_UUBAR_before_FBDT_dirname, variable_names, branch_names, Jpsi_MC_before_FBDT_values, &weights_before_FBDT, "UUBAR");
    LetsFillembeddedMC(Embedded_MC_DDBAR_before_FBDT_dirname, variable_names, branch_names, Jpsi_MC_before_FBDT_values, &weights_before_FBDT, "DDBAR");
    LetsFillembeddedMC(Embedded_MC_SSBAR_before_FBDT_dirname, variable_names, branch_names, Jpsi_MC_before_FBDT_values, &weights_before_FBDT, "SSBAR");
    LetsFillembeddedMC(Embedded_MC_CHARM_before_FBDT_dirname, variable_names, branch_names, Jpsi_MC_before_FBDT_values, &weights_before_FBDT, "CHARM");
    LetsFilldata(Embedded_data_before_FBDT_dirname, variable_names, branch_names, Jpsi_data_before_FBDT_values);

    TH1D** MC_one_bin = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** MC_before_one_bin = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** MC_efficiency = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** data_one_bin = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** data_before_one_bin = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** data_efficiency = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** Ratio_one_bin = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);

    for (int k = 0; k < (int)variable_names.size(); k++) { // malloc TH1D
        std::vector<double> temp_v;
        temp_v.insert(temp_v.end(), Jpsi_MC_values[k].begin(), Jpsi_MC_values[k].end());
        temp_v.insert(temp_v.end(), Jpsi_data_values[k].begin(), Jpsi_data_values[k].end());


        double min = *min_element(temp_v.begin(), temp_v.end());
        double max = *max_element(temp_v.begin(), temp_v.end());
        int bins = 30;
        min = 0.0;
        max = 1.0;

        MC_one_bin[k] = new TH1D((variable_names.at(k) + "_MC_one_bin").c_str(), ";number of candidates", 1, min, max);
        MC_before_one_bin[k] = new TH1D((variable_names.at(k) + "_MC_before_one_bin").c_str(), ";number of candidates", 1, min, max);
        MC_efficiency[k] = new TH1D((variable_names.at(k) + "_MC_efficiency").c_str(), ";number of candidates", 1, min, max);
        data_one_bin[k] = new TH1D((variable_names.at(k) + "_data_one_bin").c_str(), ";number of candidates", 1, min, max);
        data_before_one_bin[k] = new TH1D((variable_names.at(k) + "_data_before_one_bin").c_str(), ";number of candidates", 1, min, max);
        data_efficiency[k] = new TH1D((variable_names.at(k) + "_data_efficiency").c_str(), ";number of candidates", 1, min, max);
        Ratio_one_bin[k] = new TH1D((variable_names.at(k) + "_ratio_one_bin").c_str(), ";number of candidates", 1, min, max);
    }

    for (int k = 0; k < (int)variable_names.size(); k++) { // fill

        for (int i = 0; i < (int)Jpsi_MC_values[k].size(); i++) MC_one_bin[k]->Fill(Jpsi_MC_values[k].at(i), weights.at(i));
        for (int i = 0; i < (int)Jpsi_data_values[k].size(); i++) data_one_bin[k]->Fill(Jpsi_data_values[k].at(i));

        for (int i = 0; i < (int)Jpsi_MC_before_FBDT_values[k].size(); i++) MC_before_one_bin[k]->Fill(Jpsi_MC_before_FBDT_values[k].at(i), weights_before_FBDT.at(i));
        for (int i = 0; i < (int)Jpsi_data_before_FBDT_values[k].size(); i++) data_before_one_bin[k]->Fill(Jpsi_data_before_FBDT_values[k].at(i));
    }

    for (int k = 0; k < (int)variable_names.size(); k++) {
        MC_efficiency[k]->Divide(MC_one_bin[k], MC_before_one_bin[k]);
        data_efficiency[k]->Divide(data_one_bin[k], data_before_one_bin[k]);
        Ratio_one_bin[k]->Divide(data_efficiency[k], MC_efficiency[k]);
    }

    // Print data-MC discrepancy
    double Nevt_MC_before = MC_before_one_bin[0]->GetBinContent(1);
    double Nevt_MC_err_before = MC_before_one_bin[0]->GetBinError(1);
    double MC_relative_err_before = Nevt_MC_err_before / Nevt_MC_before;
    double Nevt_MC_after = MC_one_bin[0]->GetBinContent(1);
    double Nevt_MC_err_after = MC_one_bin[0]->GetBinError(1);
    double MC_relative_err_after = Nevt_MC_err_after / Nevt_MC_after;
    double MC_FBDT_efficiency = Nevt_MC_after / Nevt_MC_before;
    double MC_FBDT_efficiency_relative_error = sqrt(MC_relative_err_before * MC_relative_err_before + MC_relative_err_after * MC_relative_err_after - 2 * MC_relative_err_before * MC_relative_err_after);

    double Nevt_data_before = data_before_one_bin[0]->GetBinContent(1);
    double Nevt_data_err_before = data_before_one_bin[0]->GetBinError(1);
    double data_relative_err_before = Nevt_data_err_before / Nevt_data_before;
    double Nevt_data_after = data_one_bin[0]->GetBinContent(1);
    double Nevt_data_err_after = data_one_bin[0]->GetBinError(1);
    double data_relative_err_after = Nevt_data_err_after / Nevt_data_after;
    double data_FBDT_efficiency = Nevt_data_after / Nevt_data_before;
    double data_FBDT_efficiency_relative_error = sqrt(data_relative_err_before * data_relative_err_before + data_relative_err_after * data_relative_err_after - 2 * data_relative_err_before * data_relative_err_after);

    double eff_ratio = data_FBDT_efficiency / MC_FBDT_efficiency;
    double eff_ratio_relative_err = std::sqrt(MC_FBDT_efficiency_relative_error * MC_FBDT_efficiency_relative_error + data_FBDT_efficiency_relative_error * data_FBDT_efficiency_relative_error);

    printf("MC before FBDT: %lf +- %lf\n", Nevt_MC_before, Nevt_MC_err_before);
    printf("MC after FBDT: %lf +- %lf\n", Nevt_MC_after, Nevt_MC_err_after);
    printf("MC FBDT efficiency: %lf +- %lf\n", MC_FBDT_efficiency, MC_FBDT_efficiency * MC_FBDT_efficiency_relative_error);
    printf("data before FBDT: %lf +- %lf\n", Nevt_data_before, Nevt_data_err_before);
    printf("data after FBDT: %lf +- %lf\n", Nevt_data_after, Nevt_data_err_after);
    printf("data FBDT efficiency: %lf +- %lf\n", data_FBDT_efficiency, data_FBDT_efficiency * data_FBDT_efficiency_relative_error);
    printf("eff_{data}/eff_{MC}: %lf +- %lf\n", eff_ratio, eff_ratio * eff_ratio_relative_err);

    // free
    delete[] Jpsi_MC_values;
    delete[] Jpsi_MC_before_FBDT_values;
    delete[] Jpsi_data_values;
    delete[] Jpsi_data_before_FBDT_values;

    for (int k = 0; k < Nvar_num; k++) {
        delete MC_one_bin[k];
        delete MC_before_one_bin[k];
        delete MC_efficiency[k];
        delete data_one_bin[k];
        delete data_before_one_bin[k];
        delete data_efficiency[k];
        delete Ratio_one_bin[k];
    }

    delete[] MC_one_bin;
    delete[] MC_before_one_bin;
    delete[] MC_efficiency;
    delete[] data_one_bin;
    delete[] data_before_one_bin;
    delete[] data_efficiency;
    delete[] Ratio_one_bin;
}
