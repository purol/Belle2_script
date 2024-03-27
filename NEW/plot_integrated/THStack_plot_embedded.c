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

void THStack_plot_embedded() {

    // dirnames
    const char* Embedded_MC_CHG_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/SatoriRD_embedded/SIGNAL_analysis/validation_v004/final_output_root_after_MVA_Application/for_plot/CHG";
    const char* Embedded_MC_MIX_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/SatoriRD_embedded/SIGNAL_analysis/validation_v004/final_output_root_after_MVA_Application/for_plot/MIX";
    const char* Embedded_MC_UUBAR_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/SatoriRD_embedded/SIGNAL_analysis/validation_v004/final_output_root_after_MVA_Application/for_plot/UUBAR";
    const char* Embedded_MC_DDBAR_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/SatoriRD_embedded/SIGNAL_analysis/validation_v004/final_output_root_after_MVA_Application/for_plot/DDBAR";
    const char* Embedded_MC_SSBAR_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/SatoriRD_embedded/SIGNAL_analysis/validation_v004/final_output_root_after_MVA_Application/for_plot/SSBAR";
    const char* Embedded_MC_CHARM_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/SatoriRD_embedded/SIGNAL_analysis/validation_v004/final_output_root_after_MVA_Application/for_plot/CHARM";

    const char* Embedded_data_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/SatoriRD_embedded/SIGNAL_analysis/validation_v004/final_output_root_after_MVA_Application/for_plot/data";

    std::vector<std::string> variable_names;
    std::vector<std::string> branch_names;
    variable_names.push_back("nRemainingTracksInEvent"); branch_names.push_back("Upsilon");
    variable_names.push_back("Btag_chiProb"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_extraInfo_SignalProbability"); branch_names.push_back("Btag");
    variable_names.push_back("thrustAxisCosTheta"); branch_names.push_back("Upsilon");
    variable_names.push_back("missingMomentumOfEvent_theta"); branch_names.push_back("Upsilon");
    variable_names.push_back("missingEnergyOfEventCMS"); branch_names.push_back("Upsilon");
    variable_names.push_back("missingMomentumOfEvent"); branch_names.push_back("Upsilon");
    variable_names.push_back("Btag_deltaE"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_Mbc"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_useCMSFrame_theta"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_R2"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_cosTBTO"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_et"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_mm2"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso00"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso01"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso02"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso03"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso04"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso10"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso12"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso14"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso20"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso22"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso24"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hoo0"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hoo1"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hoo2"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hoo3"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hoo4"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_dr"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_dz"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_useCMSFrame_p"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_useCMSFrame_phi"); branch_names.push_back("Btag");
    variable_names.push_back("roeEextra__bocleanMask__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("extraInfo__boEeclv133__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("extraInfo__boNgammav133__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("extraInfo__boEeclv200__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("extraInfo__boNgammav200__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("Btag_thrustOm"); branch_names.push_back("Btag");
    variable_names.push_back("nParticlesInList__boe__pl__clElectronFBDT__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("nParticlesInList__bomu__pl__clMuonFBDT__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("nParticlesInList__boe__pl__clElectronFBDT_loose__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("nParticlesInList__bomu__pl__clMuonFBDT_loose__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("nParticlesInList__boe__pl__clElectronFBDT_tight__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("nParticlesInList__bomu__pl__clMuonFBDT_tight__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("Bsig_M"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_useCMSFrame_p"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dc_pValue_med"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dc_pValue_std"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dcsimpleveto_chiProb"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dcsimpleveto_dr"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dcsimpleveto_dz"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dcsimpleveto_M"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_D0_pValue_med"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_D0_pValue_std"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_D0simpleveto_chiProb"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_D0simpleveto_dr"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_D0simpleveto_dz"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_D0simpleveto_M"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_nDc_noDCS"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dc_pValue_med_noDCS"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dc_pValue_std_noDCS"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dcsimpleveto_chiProb_noDCS"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dcsimpleveto_dr_noDCS"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dcsimpleveto_dz_noDCS"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dcsimpleveto_M_noDCS"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_nDc_yespizero"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dc_pValue_med_yespizero"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dc_pValue_std_yespizero"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dcsimpleveto_chiProb_yespizero"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dcsimpleveto_dr_yespizero"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dcsimpleveto_dz_yespizero"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_Dcsimpleveto_M_yespizero"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_nD0_yespizero"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_D0_pValue_med_yespizero"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_D0_pValue_std_yespizero"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_D0simpleveto_chiProb_yespizero"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_D0simpleveto_dr_yespizero"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_D0simpleveto_dz_yespizero"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_D0simpleveto_M_yespizero"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_thrustBm"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_thrustOm"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_cosTBTO"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_cosTBz"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_et"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_mm2"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hso00"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hso01"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hso02"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hso03"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hso04"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hso10"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hso12"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hso14"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hso20"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hso22"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hso24"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hoo0"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hoo1"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hoo2"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hoo3"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_KSFWVariables_hoo4"); branch_names.push_back("Bsig");
    variable_names.push_back("foxWolframR1"); branch_names.push_back("Upsilon");
    variable_names.push_back("foxWolframR2"); branch_names.push_back("Upsilon");
    variable_names.push_back("foxWolframR3"); branch_names.push_back("Upsilon");
    variable_names.push_back("foxWolframR4"); branch_names.push_back("Upsilon");
    variable_names.push_back("cleoConeThrust0"); branch_names.push_back("Upsilon");
    variable_names.push_back("cleoConeThrust1"); branch_names.push_back("Upsilon");
    variable_names.push_back("cleoConeThrust2"); branch_names.push_back("Upsilon");
    variable_names.push_back("cleoConeThrust3"); branch_names.push_back("Upsilon");
    variable_names.push_back("cleoConeThrust4"); branch_names.push_back("Upsilon");
    variable_names.push_back("cleoConeThrust5"); branch_names.push_back("Upsilon");
    variable_names.push_back("cleoConeThrust6"); branch_names.push_back("Upsilon");
    variable_names.push_back("cleoConeThrust7"); branch_names.push_back("Upsilon");
    variable_names.push_back("cleoConeThrust8"); branch_names.push_back("Upsilon");
    variable_names.push_back("harmonicMomentThrust0"); branch_names.push_back("Upsilon");
    variable_names.push_back("harmonicMomentThrust1"); branch_names.push_back("Upsilon");
    variable_names.push_back("harmonicMomentThrust2"); branch_names.push_back("Upsilon");
    variable_names.push_back("harmonicMomentThrust3"); branch_names.push_back("Upsilon");
    variable_names.push_back("harmonicMomentThrust4"); branch_names.push_back("Upsilon");
    variable_names.push_back("sphericity"); branch_names.push_back("Upsilon");
    variable_names.push_back("aplanarity"); branch_names.push_back("Upsilon");
    variable_names.push_back("Btag_thrustBm"); branch_names.push_back("Btag");
    variable_names.push_back("roePTheta__bocleanMask__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("Btag_cosTBz"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_CleoConeCS_1"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_CleoConeCS_2"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_CleoConeCS_3"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_CleoConeCS_4"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_CleoConeCS_5"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_CleoConeCS_6"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_CleoConeCS_7"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_CleoConeCS_8"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_CleoConeCS_9"); branch_names.push_back("Btag");
    variable_names.push_back("Bsig_daughter_0_extraInfo_mychiProb"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_mydr"); branch_names.push_back("Bsig");
    variable_names.push_back("Bsig_daughter_0_extraInfo_mydz"); branch_names.push_back("Bsig");
    variable_names.push_back("nROE_KLMClusters"); branch_names.push_back("Upsilon");
    variable_names.push_back("nROE_ECLClusters__bocleanMask__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("nROE_ParticlesInList__bopi0__clmyneutralPion__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("nROE_ParticlesInList__bogamma__clmygamma__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("missingMass2OfEvent"); branch_names.push_back("Upsilon");
    variable_names.push_back("visibleEnergyOfEventCMS"); branch_names.push_back("Upsilon");
    variable_names.push_back("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc"); branch_names.push_back("Upsilon");
    variable_names.push_back("chiProb"); branch_names.push_back("Upsilon");
    variable_names.push_back("dr"); branch_names.push_back("Upsilon");
    variable_names.push_back("dz"); branch_names.push_back("Upsilon");
    variable_names.push_back("MsquaredBsig_op0"); branch_names.push_back("Upsilon");
    variable_names.push_back("MsquaredBsig_op1"); branch_names.push_back("Upsilon");
    variable_names.push_back("MsquaredBsig_op2"); branch_names.push_back("Upsilon");
    variable_names.push_back("MsquaredBsig_op3"); branch_names.push_back("Upsilon");
    variable_names.push_back("MsquaredBsig_op4"); branch_names.push_back("Upsilon");
    variable_names.push_back("MsquaredBsig_op7"); branch_names.push_back("Upsilon");
    variable_names.push_back("MVA_BB"); branch_names.push_back("Upsilon");

    Nvar_num = static_cast<int>(variable_names.size());

    std::vector<double>* Jpsi_MC_values = new std::vector<double>[Nvar_num];
    std::vector<double>* Jpsi_data_values = new std::vector<double>[Nvar_num];

    std::vector<double> weights;

    LetsFillembeddedMC(Embedded_MC_CHG_dirname, variable_names, branch_names, Jpsi_MC_values, &weights, "CHG");
    LetsFillembeddedMC(Embedded_MC_MIX_dirname, variable_names, branch_names, Jpsi_MC_values, &weights, "MIX");
    LetsFillembeddedMC(Embedded_MC_UUBAR_dirname, variable_names, branch_names, Jpsi_MC_values, &weights, "UUBAR");
    LetsFillembeddedMC(Embedded_MC_DDBAR_dirname, variable_names, branch_names, Jpsi_MC_values, &weights, "DDBAR");
    LetsFillembeddedMC(Embedded_MC_SSBAR_dirname, variable_names, branch_names, Jpsi_MC_values, &weights, "SSBAR");
    LetsFillembeddedMC(Embedded_MC_CHARM_dirname, variable_names, branch_names, Jpsi_MC_values, &weights, "CHARM");
    LetsFilldata(Embedded_data_dirname, variable_names, branch_names, Jpsi_data_values);

    TH1D** MC_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** stat_error_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** data_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** Ratio_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);

    TH1D* MC_one_bin = new TH1D("MC_one_bin", ";number of candidates", 1, -100, 100);
    TH1D* data_one_bin = new TH1D("data_one_bin", ";number of candidates", 1, -100, 100);
    TH1D* Ratio_one_bin = new TH1D("ratio_one_bin", ";number of candidates", 1, -100, 100);

    for (int k = 0; k < (int)variable_names.size(); k++) { // malloc TH1D
        std::vector<double> temp_v;
        temp_v.insert(temp_v.end(), Jpsi_MC_values[k].begin(), Jpsi_MC_values[k].end());
        temp_v.insert(temp_v.end(), Jpsi_data_values[k].begin(), Jpsi_data_values[k].end());


        double min = *min_element(temp_v.begin(), temp_v.end());
        double max = *max_element(temp_v.begin(), temp_v.end());
        int bins = 30;

        if (hasEnding(variable_names.at(k), std::string("dr"))) { // exceptions
            max = 0.2;
            min = 0.0;
        }
        else if (hasEnding(variable_names.at(k), std::string("dz"))) {
            max = 0.2;
            min = -0.2;
        }
        else if (hasEnding(variable_names.at(k), std::string("Bsig_M"))) {
            max = 2.0;
            min = 0.0;
        }
        else if (hasEnding(variable_names.at(k), std::string("chiProb"))) {
            max = 1.0;
            min = 0.0;
        }
        else if (hasEnding(variable_names.at(k), std::string("harmonicMomentThrust1"))) {
            max = 0.6;
            min = -0.6;
        }
        else if (hasEnding(variable_names.at(k), std::string("harmonicMomentThrust2"))) {
            max = 1.0;
            min = 0.0;
        }
        else if (hasEnding(variable_names.at(k), std::string("harmonicMomentThrust3"))) {
            max = 1.0;
            min = -1.0;
        }
        else if (hasEnding(variable_names.at(k), std::string("harmonicMomentThrust4"))) {
            max = 1.0;
            min = -0.5;
        }
        else if (hasEnding(variable_names.at(k), std::string("missingMomentumOfEvent"))) {
            max = 5.0;
        }
        else if (hasEnding(variable_names.at(k), std::string("missingEnergyOfEventCMS"))) {
            min = -1.5;
        }
        else if (hasEnding(variable_names.at(k), std::string("Btag_extraInfo_SignalProbability"))) {
            max = 0;
            min = -3;
            variable_names.at(k) = std::string("log_{10}SignalProbability");
        }
        else if (hasEnding(variable_names.at(k), std::string("Btag_thrustOm"))) {
            min = 0.5;
        }
        else if (hasEnding(variable_names.at(k), std::string("foxWolframR1"))) {
            max = 0.25;
        }
        else if (hasEnding(variable_names.at(k), std::string("foxWolframR3"))) {
            max = 0.4;
        }
        else if (hasEnding(variable_names.at(k), std::string("missingMass2OfEvent"))) {
            min = -20.0;
            max = 40.0;
        }
        else if (hasEnding(variable_names.at(k), std::string("visibleEnergyOfEventCMS"))) {
            max = 15.0;
        }
        else if (hasEnding(variable_names.at(k), std::string("Btag_KSFWVariables_hoo4"))) {
            max = 0.08;
        }
        else if (hasEnding(variable_names.at(k), std::string("Bsig_KSFWVariables_et"))) {
            max = 10;
        }
        else if (hasEnding(variable_names.at(k), std::string("Bsig_KSFWVariables_hso24"))) {
            max = 0.2;
        }
        else if (hasEnding(variable_names.at(k), std::string("MsquaredBsig_op0"))) {
            min = 7.0;
            max = 33.0;
        }
        else if (hasEnding(variable_names.at(k), std::string("MsquaredBsig_op1"))) {
            min = 0.0;
            max = 40.0;
        }
        else if (hasEnding(variable_names.at(k), std::string("MsquaredBsig_op2"))) {
            min = -1.5;
            max = 0.1;
        }
        else if (hasEnding(variable_names.at(k), std::string("MsquaredBsig_op3"))) {
            min = -0.4;
            max = 0.1;
        }
        else if (hasEnding(variable_names.at(k), std::string("MsquaredBsig_op4"))) {
            min = -0.6;
            max = 0.1;
        }
        else if (hasEnding(variable_names.at(k), std::string("MsquaredBsig_op7"))) {
            min = 0.0;
            max = 15.0;
        }
        else if ((variable_names.at(k).find("Bsig_daughter_0_extraInfo_D") != std::string::npos) && (variable_names.at(k).find("_M") != std::string::npos)) {
            min = 0.6;
            max = 2.1;
        }

        MC_hist[k] = new TH1D("embedded MC", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        stat_error_hist[k] = new TH1D("MC stat error", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        data_hist[k] = new TH1D("embedded data", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        Ratio_hist[k] = new TH1D((variable_names.at(k) + "_ratio").c_str(), ";;data/MC", bins, min, max);
    }

    int index = std::find(variable_names.begin(), variable_names.end(), std::string("log_{10}SignalProbability")) - variable_names.begin();
    for (int i = 0; i < (int)Jpsi_MC_values[index].size(); i++) Jpsi_MC_values[index].at(i) = log10l(Jpsi_MC_values[index].at(i));
    for (int i = 0; i < (int)Jpsi_data_values[index].size(); i++) Jpsi_data_values[index].at(i) = log10l(Jpsi_data_values[index].at(i));

    for (int k = 0; k < (int)variable_names.size(); k++) { // fill
        for (int i = 0; i < (int)Jpsi_MC_values[k].size(); i++) {
            MC_hist[k]->Fill(Jpsi_MC_values[k].at(i), weights.at(i));
            stat_error_hist[k]->Fill(Jpsi_MC_values[k].at(i), weights.at(i));
        }
        for (int i = 0; i < (int)Jpsi_data_values[k].size(); i++) data_hist[k]->Fill(Jpsi_data_values[k].at(i));

        if (variable_names.at(k) == "nROE_ParticlesInList__bopi0__clmyneutralPion__bc") {
            for (int i = 0; i < (int)Jpsi_MC_values[k].size(); i++) MC_one_bin->Fill(Jpsi_MC_values[k].at(i), weights.at(i));
            for (int i = 0; i < (int)Jpsi_data_values[k].size(); i++) data_one_bin->Fill(Jpsi_data_values[k].at(i));
        }
    }
    Ratio_one_bin->Divide(data_one_bin, MC_one_bin);

    printf("MC: %d\n", (int)Jpsi_MC_values[0].size());
    printf("data: %d\n", (int)Jpsi_data_values[0].size());

    for (int k = 0; k < (int)variable_names.size(); k++) { // draw

        Ratio_hist[k]->SetLineColor(kBlack); Ratio_hist[k]->SetMarkerStyle(21); Ratio_hist[k]->Sumw2(); Ratio_hist[k]->SetStats(0);
        Ratio_hist[k]->Divide(data_hist[k], MC_hist[k]);

        TCanvas* c_temp = new TCanvas("c", "", 800, 800); c_temp->cd();

        TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.35, 1.0, 1.0);
        pad1->SetBottomMargin(0.08); pad1->SetLeftMargin(0.15);
        pad1->SetGridx(); pad1->Draw(); pad1->cd();

        gStyle->SetPalette(kPastel);

        Float_t ymax_1 = MC_hist[k]->GetMaximum();
        Float_t ymax_2 = data_hist[k]->GetMaximum();
        double real_max = 0;
        if (ymax_1 > ymax_2) real_max = ymax_1;
        else real_max = ymax_2;

        MC_hist[k]->SetStats(0);
        MC_hist[k]->SetFillStyle(3001);
        MC_hist[k]->SetLineColor(33);
        MC_hist[k]->SetFillColor(33);
        MC_hist[k]->GetYaxis()->SetRangeUser(0.0, real_max * 1.1);
        MC_hist[k]->Draw("Hist");

        stat_error_hist[k]->SetFillColor(12); stat_error_hist[k]->SetLineWidth(0); stat_error_hist[k]->SetFillStyle(3004); stat_error_hist[k]->Draw("e2 SAME");
        data_hist[k]->SetLineWidth(2); data_hist[k]->SetLineColor(kBlack); data_hist[k]->SetMarkerStyle(8); data_hist[k]->Draw("SAME eP");
        TLegend* legend = pad1->BuildLegend(0.9, 0.9, 0.7, 0.7);
        legend->SetFillStyle(0); legend->SetLineWidth(0);
        TPaveText* pt = new TPaveText(0.135, 0.88, 0.5, 1.0, "NDC NB"); pt->SetFillStyle(0); pt->SetLineWidth(0); pt->AddText(("MC scaled to data, Data/MC= " + std::to_string(CAL)).c_str()); pt->Draw();

        c_temp->cd();
        TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1, 0.3); pad2->SetBottomMargin(0.15); pad2->SetLeftMargin(0.15); pad2->SetGridx(); pad2->Draw(); pad2->cd();
        Ratio_hist[k]->SetMinimum(0.5); Ratio_hist[k]->SetMaximum(1.5); Ratio_hist[k]->SetLineWidth(2);
        Ratio_hist[k]->GetYaxis()->SetTitleSize(0.08); Ratio_hist[k]->GetYaxis()->SetTitleOffset(0.5);
        Ratio_hist[k]->GetXaxis()->SetLabelSize(0.08); Ratio_hist[k]->GetYaxis()->SetLabelSize(0.08);
        Ratio_hist[k]->Draw("e0p");
        TLine* line = new TLine(Ratio_hist[k]->GetXaxis()->GetXmin(), 1.0, Ratio_hist[k]->GetXaxis()->GetXmax(), 1.0);
        line->SetLineColor(kRed);
        line->SetLineStyle(1); line->SetLineWidth(3);
        line->Draw();

        c_temp->SetBottomMargin(0.0);
        c_temp->SaveAs((variable_names.at(k) + "_embedded_Jpsi.png").c_str());

        delete c_temp;
    }

    // Print data-MC discrepancy
    double MC_sum = 0;
    for (int i = 0; i < (int)Jpsi_MC_values[0].size(); i++) MC_sum = MC_sum + weights.at(i);
    printf("data num: %ld\n", Jpsi_data_values[0].size());
    printf("MC num with calibration: %lf\n", MC_sum);
    printf("MC with calibration: %lf +- %lf\n", MC_one_bin->GetBinContent(1), MC_one_bin->GetBinError(1));
    printf("data with calibration: %lf +- %lf\n", data_one_bin->GetBinContent(1), data_one_bin->GetBinError(1));
    printf("data/MC with calibration: %lf +- %lf\n", Ratio_one_bin->GetBinContent(1), Ratio_one_bin->GetBinError(1));

    // free
    delete[] Jpsi_MC_values;
    delete[] Jpsi_data_values;

    for (int k = 0; k < Nvar_num; k++) {
        delete MC_hist[k];
        delete stat_error_hist[k];
        delete data_hist[k];
        delete Ratio_hist[k];
    }

    delete[] MC_hist;
    delete[] stat_error_hist;
    delete[] data_hist;
    delete[] Ratio_hist;

    delete MC_one_bin;
    delete data_one_bin;
    delete Ratio_one_bin;
}
