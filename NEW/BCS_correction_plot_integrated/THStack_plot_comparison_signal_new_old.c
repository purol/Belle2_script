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
# include "THStack_plot_MC15rdVSMC15ri.h"
#include "template.h"

/* ====================================== */

// for background
// (old sample + new FBDT) VS (new sample + new FBDT)

void THStack_plot_comparison_signal_new_old() {

    const char* old_signal_validation_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD/SIGNAL_analysis/validation_v010/final_output_root_after_MVA_Application/";

    const char* new_signal_validation_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD/SIGNAL_analysis/validation_v009/final_output_root_after_MVA_Application/";

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
    variable_names.push_back("bin index"); branch_names.push_back("Upsilon");

    Nvar_num = static_cast<int>(variable_names.size());

    std::vector<double>* MC15rd_old_values = new std::vector<double>[Nvar_num];
    std::vector<double>* MC15rd_values = new std::vector<double>[Nvar_num];
    std::vector<int> MC_numbering;

    std::vector<double>* data_values = new std::vector<double>[Nvar_num];

    std::vector<double> weights;
    std::vector<double> MC15rd_old_weights;
    std::vector<double> MC15rd_weights;

    double CAL_MC15rd_old = 1.0;
    double CAL_MC15rd = 1.0;

    LetsFill_ri(old_signal_validation_dirname, variable_names, branch_names, MC15rd_old_values, &MC_numbering, &MC15rd_old_weights, "SIGNAL", "MC15rd", "validation", CAL_MC15rd_old, 1.0, 1.0, 1.0, "B2Knunu", true);
    LetsFill_ri(old_signal_validation_dirname, variable_names, branch_names, MC15rd_old_values, &MC_numbering, &MC15rd_old_weights, "SIGNAL", "MC15rd", "validation", CAL_MC15rd_old, 1.0, 1.0, 1.0, "B2Kstarnunu", true);
    LetsFill_ri(old_signal_validation_dirname, variable_names, branch_names, MC15rd_old_values, &MC_numbering, &MC15rd_old_weights, "SIGNAL", "MC15rd", "validation", CAL_MC15rd_old, 1.0, 1.0, 1.0, "B2Xsnunu", true);
    LetsFill_ri(old_signal_validation_dirname, variable_names, branch_names, MC15rd_old_values, &MC_numbering, &MC15rd_old_weights, "SIGNAL", "MC15rd", "validation", CAL_MC15rd_old, 1.0, 1.0, 1.0, "B02K0nunu", true);
    LetsFill_ri(old_signal_validation_dirname, variable_names, branch_names, MC15rd_old_values, &MC_numbering, &MC15rd_old_weights, "SIGNAL", "MC15rd", "validation", CAL_MC15rd_old, 1.0, 1.0, 1.0, "B02Kstar0nunu", true);
    LetsFill_ri(old_signal_validation_dirname, variable_names, branch_names, MC15rd_old_values, &MC_numbering, &MC15rd_old_weights, "SIGNAL", "MC15rd", "validation", CAL_MC15rd_old, 1.0, 1.0, 1.0, "B02Xsnunu", true);
    LetsFill_ri(new_signal_validation_dirname, variable_names, branch_names, MC15rd_values, &MC_numbering, &MC15rd_weights, "SIGNAL", "MC15rd", "validation", CAL_MC15rd, 1.0, 1.0, 1.0, "B2Knunu", false);
    LetsFill_ri(new_signal_validation_dirname, variable_names, branch_names, MC15rd_values, &MC_numbering, &MC15rd_weights, "SIGNAL", "MC15rd", "validation", CAL_MC15rd, 1.0, 1.0, 1.0, "B2Kstarnunu", false);
    LetsFill_ri(new_signal_validation_dirname, variable_names, branch_names, MC15rd_values, &MC_numbering, &MC15rd_weights, "SIGNAL", "MC15rd", "validation", CAL_MC15rd, 1.0, 1.0, 1.0, "B2Xsnunu", false);
    LetsFill_ri(new_signal_validation_dirname, variable_names, branch_names, MC15rd_values, &MC_numbering, &MC15rd_weights, "SIGNAL", "MC15rd", "validation", CAL_MC15rd, 1.0, 1.0, 1.0, "B02K0nunu", false);
    LetsFill_ri(new_signal_validation_dirname, variable_names, branch_names, MC15rd_values, &MC_numbering, &MC15rd_weights, "SIGNAL", "MC15rd", "validation", CAL_MC15rd, 1.0, 1.0, 1.0, "B02Kstar0nunu", false);
    LetsFill_ri(new_signal_validation_dirname, variable_names, branch_names, MC15rd_values, &MC_numbering, &MC15rd_weights, "SIGNAL", "MC15rd", "validation", CAL_MC15rd, 1.0, 1.0, 1.0, "B02Xsnunu", false);

    THStack** Stack = (THStack**)malloc(sizeof(THStack*) * Nvar_num);
    TH1D** MC15rd_old_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** MC15rd_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** stat_error_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** data_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** Ratio_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);

    for (int k = 0; k < (int)variable_names.size(); k++) { // malloc TH1D
        std::vector<double> temp_v;
        temp_v.insert(temp_v.end(), MC15rd_old_values[k].begin(), MC15rd_old_values[k].end());
        temp_v.insert(temp_v.end(), MC15rd_values[k].begin(), MC15rd_values[k].end());
        temp_v.insert(temp_v.end(), data_values[k].begin(), data_values[k].end());


        double min = *min_element(temp_v.begin(), temp_v.end());
        double max = *max_element(temp_v.begin(), temp_v.end());
        int bins = 100;

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
        else if (hasEnding(variable_names.at(k), std::string("bin index"))) {
            min = 0.0;
            max = RarityBins;
            bins = RarityBins;
        }

        Stack[k] = new THStack(variable_names.at(k).c_str(), (";" + variable_names.at(k) + ";number of candidates").c_str());
        MC15rd_old_hist[k] = new TH1D("old signal", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        MC15rd_hist[k] = new TH1D("new signal", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        stat_error_hist[k] = new TH1D("MC stat error", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        data_hist[k] = new TH1D("data", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        Ratio_hist[k] = new TH1D((variable_names.at(k) + "_ratio").c_str(), ";;new/old", bins, min, max);
    }

    int index = std::find(variable_names.begin(), variable_names.end(), std::string("log_{10}SignalProbability")) - variable_names.begin();
    for (int i = 0; i < (int)MC15rd_old_values[index].size(); i++) MC15rd_old_values[index].at(i) = log10l(MC15rd_old_values[index].at(i));
    for (int i = 0; i < (int)MC15rd_values[index].size(); i++) MC15rd_values[index].at(i) = log10l(MC15rd_values[index].at(i));
    for (int i = 0; i < (int)data_values[index].size(); i++) data_values[index].at(i) = log10l(data_values[index].at(i));

    for (int k = 0; k < (int)variable_names.size(); k++) { // fill
        for (int i = 0; i < (int)MC15rd_old_values[k].size(); i++) MC15rd_old_hist[k]->Fill(MC15rd_old_values[k].at(i), MC15rd_old_weights.at(i));
        for (int i = 0; i < (int)MC15rd_values[k].size(); i++) MC15rd_hist[k]->Fill(MC15rd_values[k].at(i), MC15rd_weights.at(i));
        for (int i = 0; i < (int)data_values[k].size(); i++) data_hist[k]->Fill(data_values[k].at(i));
    }

    double MC15rd_old_sum = 0;
    double MC15rd_sum = 0;
    for (int i = 0; i < (int)MC15rd_old_values[0].size(); i++) MC15rd_old_sum = MC15rd_old_sum + MC15rd_old_weights.at(i);
    for (int i = 0; i < (int)MC15rd_values[0].size(); i++) MC15rd_sum = MC15rd_sum + MC15rd_weights.at(i);

    for (int k = 0; k < (int)variable_names.size(); k++) { // draw
        // Scale the histogram
        //CAL_MC15rd_old = MC15rd_sum / MC15rd_old_sum;
        //CAL_MC15rd = 1.0;
        //MC15rd_hist[k]->Scale(1.0);
        //MC15rd_old_hist[k]->Scale(CAL_MC15rd_old);

        //Stack[k]->Add(charged_hist[k]);
        //Stack[k]->Add(mixed_hist[k]);
        //Stack[k]->Add(uubar_hist[k]);
        //Stack[k]->Add(ddbar_hist[k]);
        //Stack[k]->Add(ssbar_hist[k]);
        //Stack[k]->Add(ccbar_hist[k]);
        //Stack[k]->Add(taupair_hist[k]);
        //Stack[k]->Add(mumu_hist[k]);
        //Stack[k]->Add(gg_hist[k]);
        //Stack[k]->Add(ee_hist[k]);
        //Stack[k]->Add(eeee_hist[k]);
        //Stack[k]->Add(eemumu_hist[k]);
        //Stack[k]->Add(llXX_hist[k]);
        //Stack[k]->Add(hhISR_hist[k]);

        Ratio_hist[k]->SetLineColor(kBlack); Ratio_hist[k]->SetMarkerStyle(21); Ratio_hist[k]->Sumw2(); Ratio_hist[k]->SetStats(0);
        Ratio_hist[k]->Divide(MC15rd_hist[k], MC15rd_old_hist[k]);

        TCanvas* c_temp = new TCanvas("c", "", 800, 800); c_temp->cd();

        TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.35, 1.0, 1.0);
        pad1->SetBottomMargin(0.08); pad1->SetLeftMargin(0.15);
        pad1->SetGridx(); pad1->Draw(); pad1->cd();

        gStyle->SetPalette(kPastel);

        Float_t ymax_1 = MC15rd_hist[k]->GetMaximum();
        Float_t ymax_2 = MC15rd_old_hist[k]->GetMaximum();
        double real_max = 0;
        if (ymax_1 > ymax_2) real_max = ymax_1;
        else real_max = ymax_2;

        MC15rd_hist[k]->SetMaximum(real_max * 1.1);

        MC15rd_old_hist[k]->SetStats(0);
        MC15rd_hist[k]->SetStats(0);

        MC15rd_hist[k]->SetFillColor(4); MC15rd_hist[k]->SetLineColor(4); MC15rd_hist[k]->SetLineWidth(1); MC15rd_hist[k]->SetFillStyle(3002); MC15rd_hist[k]->Draw("e1 Hist");
        MC15rd_old_hist[k]->SetFillColor(2); MC15rd_old_hist[k]->SetLineColor(2); MC15rd_old_hist[k]->SetLineWidth(1); MC15rd_old_hist[k]->SetFillStyle(3003); MC15rd_old_hist[k]->Draw("e1 Hist SAME");
        TLegend* legend = pad1->BuildLegend(0.35, 0.9, 0.15, 0.6);
        legend->SetFillStyle(0); legend->SetLineWidth(0);
        //TPaveText* pt = new TPaveText(0.135, 0.88, 0.6, 1.0, "NDC NB"); pt->SetFillStyle(0); pt->SetLineWidth(0); pt->AddText(("MC15rd_old scaled to MC15rd, MC15rd/MC15rd_old= " + std::to_string(CAL_MC15rd_old/ CAL_MC15rd)).c_str()); pt->Draw();

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
        c_temp->SaveAs((variable_names.at(k) + "_MC15rd_old_vs_MC15rd_signal.png").c_str());

        delete c_temp;
    }

    // Print new-old discrepancy
    printf("data num: %ld\n", data_values[0].size());
    printf("old num with calibration: %lf\n", MC15rd_old_sum);
    printf("new num with calibration: %lf\n", MC15rd_sum);

    // free
    delete[] MC15rd_old_values;
    delete[] MC15rd_values;

    delete[] data_values;

    for (int k = 0; k < Nvar_num; k++) {
        delete Stack[k];
        delete MC15rd_old_hist[k];
        delete MC15rd_hist[k];
        delete stat_error_hist[k];
        delete data_hist[k];
        delete Ratio_hist[k];
    }

    free(Stack);
    free(MC15rd_old_hist);
    free(MC15rd_hist);
    free(stat_error_hist);
    free(data_hist);
    free(Ratio_hist);
}
