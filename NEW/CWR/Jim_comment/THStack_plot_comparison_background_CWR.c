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
# include "THStack_plot_el7vsel9.h"
#include "template.h"

/* ====================================== */

void THStack_plot_comparison_background_CWR() {

    const char* el7_CHG_validation_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/CHG_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut/";
    const char* el7_MIX_validation_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/MIX_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut/";
    const char* el7_UUBAR_validation_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/UUBAR_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut/";
    const char* el7_DDBAR_validation_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/DDBAR_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut/";
    const char* el7_SSBAR_validation_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/SSBAR_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut/";
    const char* el7_CHARM_validation_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD/CHARM_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut/";

    const char* el9_CHG_validation_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/ByakurenRD/CHG_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut/";
    const char* el9_MIX_validation_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/ByakurenRD/MIX_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut/";
    const char* el9_UUBAR_validation_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/ByakurenRD/UUBAR_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut/";
    const char* el9_DDBAR_validation_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/ByakurenRD/DDBAR_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut/";
    const char* el9_SSBAR_validation_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/ByakurenRD/SSBAR_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut/";
    const char* el9_CHARM_validation_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/ByakurenRD/CHARM_analysis/validation_v009/final_output_root_after_MVA_Application_after_cut/";

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

    std::vector<double>* MC_values = new std::vector<double>[Nvar_num];
    std::vector<double>* charged_el7_values = new std::vector<double>[Nvar_num];
    std::vector<double>* mixed_el7_values = new std::vector<double>[Nvar_num];
    std::vector<double>* uubar_el7_values = new std::vector<double>[Nvar_num];
    std::vector<double>* ddbar_el7_values = new std::vector<double>[Nvar_num];
    std::vector<double>* ssbar_el7_values = new std::vector<double>[Nvar_num];
    std::vector<double>* ccbar_el7_values = new std::vector<double>[Nvar_num];
    std::vector<double>* taupair_el7_values = new std::vector<double>[Nvar_num];
    std::vector<double>* mumu_el7_values = new std::vector<double>[Nvar_num];
    std::vector<double>* gg_el7_values = new std::vector<double>[Nvar_num];
    std::vector<double>* ee_el7_values = new std::vector<double>[Nvar_num];
    std::vector<double>* eeee_el7_values = new std::vector<double>[Nvar_num];
    std::vector<double>* eemumu_el7_values = new std::vector<double>[Nvar_num];
    std::vector<double>* llXX_el7_values = new std::vector<double>[Nvar_num];
    std::vector<double>* hhISR_el7_values = new std::vector<double>[Nvar_num];
    std::vector<double>* charged_el9_values = new std::vector<double>[Nvar_num];
    std::vector<double>* mixed_el9_values = new std::vector<double>[Nvar_num];
    std::vector<double>* uubar_el9_values = new std::vector<double>[Nvar_num];
    std::vector<double>* ddbar_el9_values = new std::vector<double>[Nvar_num];
    std::vector<double>* ssbar_el9_values = new std::vector<double>[Nvar_num];
    std::vector<double>* ccbar_el9_values = new std::vector<double>[Nvar_num];
    std::vector<double>* el7_values = new std::vector<double>[Nvar_num];
    std::vector<double>* el9_values = new std::vector<double>[Nvar_num];
    std::vector<int> MC_numbering;

    std::vector<double>* data_values = new std::vector<double>[Nvar_num];

    std::vector<double> weights;
    std::vector<double> charged_el7_weights;
    std::vector<double> mixed_el7_weights;
    std::vector<double> uubar_el7_weights;
    std::vector<double> ddbar_el7_weights;
    std::vector<double> ssbar_el7_weights;
    std::vector<double> ccbar_el7_weights;
    std::vector<double> taupair_el7_weights;
    std::vector<double> mumu_el7_weights;
    std::vector<double> gg_el7_weights;
    std::vector<double> ee_el7_weights;
    std::vector<double> eeee_el7_weights;
    std::vector<double> eemumu_el7_weights;
    std::vector<double> llXX_el7_weights;
    std::vector<double> hhISR_el7_weights;
    std::vector<double> charged_el9_weights;
    std::vector<double> mixed_el9_weights;
    std::vector<double> uubar_el9_weights;
    std::vector<double> ddbar_el9_weights;
    std::vector<double> ssbar_el9_weights;
    std::vector<double> ccbar_el9_weights;
    std::vector<double> el7_weights;
    std::vector<double> el9_weights;

    double CAL_el7 = 1.0;
    double CAL_el9 = 1.0;

    LetsFill_ri(el7_CHG_validation_dirname, variable_names, branch_names, MC_values, &MC_numbering, &weights, "CHG", "el7", "validation", CAL_el7, 1.0, 1.0, 1.0, "root");
    LetsFill_ri(el7_MIX_validation_dirname, variable_names, branch_names, MC_values, &MC_numbering, &weights, "MIX", "el7", "validation", CAL_el7, 1.0, 1.0, 1.0, "root");
    LetsFill_ri(el7_UUBAR_validation_dirname, variable_names, branch_names, MC_values, &MC_numbering, &weights, "UUBAR", "el7", "validation", CAL_el7, 1.0, 1.0, 1.0, "root");
    LetsFill_ri(el7_DDBAR_validation_dirname, variable_names, branch_names, MC_values, &MC_numbering, &weights, "DDBAR", "el7", "validation", CAL_el7, 1.0, 1.0, 1.0, "root");
    LetsFill_ri(el7_SSBAR_validation_dirname, variable_names, branch_names, MC_values, &MC_numbering, &weights, "SSBAR", "el7", "validation", CAL_el7, 1.0, 1.0, 1.0, "root");
    LetsFill_ri(el7_CHARM_validation_dirname, variable_names, branch_names, MC_values, &MC_numbering, &weights, "CHARM", "el7", "validation", CAL_el7, 1.0, 1.0, 1.0, "root");
    LetsFill_ri(el9_CHG_validation_dirname, variable_names, branch_names, MC_values, &MC_numbering, &weights, "CHG", "el9", "validation", CAL_el9, 1.0, 1.0, 1.0, "root");
    LetsFill_ri(el9_MIX_validation_dirname, variable_names, branch_names, MC_values, &MC_numbering, &weights, "MIX", "el9", "validation", CAL_el9, 1.0, 1.0, 1.0, "root");
    LetsFill_ri(el9_UUBAR_validation_dirname, variable_names, branch_names, MC_values, &MC_numbering, &weights, "UUBAR", "el9", "validation", CAL_el9, 1.0, 1.0, 1.0, "root");
    LetsFill_ri(el9_DDBAR_validation_dirname, variable_names, branch_names, MC_values, &MC_numbering, &weights, "DDBAR", "el9", "validation", CAL_el9, 1.0, 1.0, 1.0, "root");
    LetsFill_ri(el9_SSBAR_validation_dirname, variable_names, branch_names, MC_values, &MC_numbering, &weights, "SSBAR", "el9", "validation", CAL_el9, 1.0, 1.0, 1.0, "root");
    LetsFill_ri(el9_CHARM_validation_dirname, variable_names, branch_names, MC_values, &MC_numbering, &weights, "CHARM", "el9", "validation", CAL_el9, 1.0, 1.0, 1.0, "root");

    // sort variables
    for (int k = 0; k < (int)MC_numbering.size(); k++) {
        if (MC_numbering.at(k) == 0) {
            for (int l = 0; l < (int)variable_names.size(); l++) charged_el7_values[l].push_back(MC_values[l].at(k));
            charged_el7_weights.push_back(weights.at(k));
        }
        else if (MC_numbering.at(k) == 1) {
            for (int l = 0; l < (int)variable_names.size(); l++) mixed_el7_values[l].push_back(MC_values[l].at(k));
            mixed_el7_weights.push_back(weights.at(k));
        }
        else if (MC_numbering.at(k) == 2) {
            for (int l = 0; l < (int)variable_names.size(); l++) uubar_el7_values[l].push_back(MC_values[l].at(k));
            uubar_el7_weights.push_back(weights.at(k));
        }
        else if (MC_numbering.at(k) == 3) {
            for (int l = 0; l < (int)variable_names.size(); l++) ddbar_el7_values[l].push_back(MC_values[l].at(k));
            ddbar_el7_weights.push_back(weights.at(k));
        }
        else if (MC_numbering.at(k) == 4) {
            for (int l = 0; l < (int)variable_names.size(); l++) ssbar_el7_values[l].push_back(MC_values[l].at(k));
            ssbar_el7_weights.push_back(weights.at(k));
        }
        else if (MC_numbering.at(k) == 5) {
            for (int l = 0; l < (int)variable_names.size(); l++) ccbar_el7_values[l].push_back(MC_values[l].at(k));
            ccbar_el7_weights.push_back(weights.at(k));
        }
        else if (MC_numbering.at(k) == 6) {
            for (int l = 0; l < (int)variable_names.size(); l++) taupair_el7_values[l].push_back(MC_values[l].at(k));
            taupair_el7_weights.push_back(weights.at(k));
        }
        else if (MC_numbering.at(k) == 7) {
            for (int l = 0; l < (int)variable_names.size(); l++) mumu_el7_values[l].push_back(MC_values[l].at(k));
            mumu_el7_weights.push_back(weights.at(k));
        }
        else if (MC_numbering.at(k) == 8) {
            for (int l = 0; l < (int)variable_names.size(); l++) gg_el7_values[l].push_back(MC_values[l].at(k));
            gg_el7_weights.push_back(weights.at(k));
        }
        else if (MC_numbering.at(k) == 9) {
            for (int l = 0; l < (int)variable_names.size(); l++) ee_el7_values[l].push_back(MC_values[l].at(k));
            ee_el7_weights.push_back(weights.at(k));
        }
        else if (MC_numbering.at(k) == 10) {
            for (int l = 0; l < (int)variable_names.size(); l++) eeee_el7_values[l].push_back(MC_values[l].at(k));
            eeee_el7_weights.push_back(weights.at(k));
        }
        else if (MC_numbering.at(k) == 11) {
            for (int l = 0; l < (int)variable_names.size(); l++) eemumu_el7_values[l].push_back(MC_values[l].at(k));
            eemumu_el7_weights.push_back(weights.at(k));
        }
        else if (MC_numbering.at(k) == 12) {
            for (int l = 0; l < (int)variable_names.size(); l++) llXX_el7_values[l].push_back(MC_values[l].at(k));
            llXX_el7_weights.push_back(weights.at(k));
        }
        else if (MC_numbering.at(k) == 13) {
            for (int l = 0; l < (int)variable_names.size(); l++) hhISR_el7_values[l].push_back(MC_values[l].at(k));
            hhISR_el7_weights.push_back(weights.at(k));
        }
        else if (MC_numbering.at(k) == (0 + IndexShift)) {
            for (int l = 0; l < (int)variable_names.size(); l++) charged_el9_values[l].push_back(MC_values[l].at(k));
            charged_el9_weights.push_back(weights.at(k));
        }
        else if (MC_numbering.at(k) == (1 + IndexShift)) {
            for (int l = 0; l < (int)variable_names.size(); l++) mixed_el9_values[l].push_back(MC_values[l].at(k));
            mixed_el9_weights.push_back(weights.at(k));
        }
        else if (MC_numbering.at(k) == (2 + IndexShift)) {
            for (int l = 0; l < (int)variable_names.size(); l++) uubar_el9_values[l].push_back(MC_values[l].at(k));
            uubar_el9_weights.push_back(weights.at(k));
        }
        else if (MC_numbering.at(k) == (3 + IndexShift)) {
            for (int l = 0; l < (int)variable_names.size(); l++) ddbar_el9_values[l].push_back(MC_values[l].at(k));
            ddbar_el9_weights.push_back(weights.at(k));
        }
        else if (MC_numbering.at(k) == (4 + IndexShift)) {
            for (int l = 0; l < (int)variable_names.size(); l++) ssbar_el9_values[l].push_back(MC_values[l].at(k));
            ssbar_el9_weights.push_back(weights.at(k));
        }
        else if (MC_numbering.at(k) == (5 + IndexShift)) {
            for (int l = 0; l < (int)variable_names.size(); l++) ccbar_el9_values[l].push_back(MC_values[l].at(k));
            ccbar_el9_weights.push_back(weights.at(k));
        }
        else {
            printf("undefined numbering!\n");
            exit(1);
        }

        // for el7 vs el9
        if (MC_numbering.at(k) < IndexShift) { // el7
            for (int l = 0; l < (int)variable_names.size(); l++) el7_values[l].push_back(MC_values[l].at(k));
            el7_weights.push_back(weights.at(k));
        }
        else { // el9
            for (int l = 0; l < (int)variable_names.size(); l++) el9_values[l].push_back(MC_values[l].at(k));
            el9_weights.push_back(weights.at(k));
        }
    }

    THStack** Stack = (THStack**)malloc(sizeof(THStack*) * Nvar_num);
    TH1D** charged_el7_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** mixed_el7_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** uubar_el7_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** ddbar_el7_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** ssbar_el7_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** ccbar_el7_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** taupair_el7_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** mumu_el7_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** gg_el7_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** ee_el7_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** eeee_el7_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** eemumu_el7_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** llXX_el7_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** hhISR_el7_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** charged_el9_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** mixed_el9_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** uubar_el9_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** ddbar_el9_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** ssbar_el9_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** ccbar_el9_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** el7_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** el9_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** stat_error_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** data_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** Ratio_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);

    for (int k = 0; k < (int)variable_names.size(); k++) { // malloc TH1D
        std::vector<double> temp_v;
        temp_v.insert(temp_v.end(), charged_el7_values[k].begin(), charged_el7_values[k].end());
        temp_v.insert(temp_v.end(), mixed_el7_values[k].begin(), mixed_el7_values[k].end());
        temp_v.insert(temp_v.end(), uubar_el7_values[k].begin(), uubar_el7_values[k].end());
        temp_v.insert(temp_v.end(), ddbar_el7_values[k].begin(), ddbar_el7_values[k].end());
        temp_v.insert(temp_v.end(), ssbar_el7_values[k].begin(), ssbar_el7_values[k].end());
        temp_v.insert(temp_v.end(), ccbar_el7_values[k].begin(), ccbar_el7_values[k].end());
        temp_v.insert(temp_v.end(), taupair_el7_values[k].begin(), taupair_el7_values[k].end());
        temp_v.insert(temp_v.end(), mumu_el7_values[k].begin(), mumu_el7_values[k].end());
        temp_v.insert(temp_v.end(), gg_el7_values[k].begin(), gg_el7_values[k].end());
        temp_v.insert(temp_v.end(), ee_el7_values[k].begin(), ee_el7_values[k].end());
        temp_v.insert(temp_v.end(), eeee_el7_values[k].begin(), eeee_el7_values[k].end());
        temp_v.insert(temp_v.end(), eemumu_el7_values[k].begin(), eemumu_el7_values[k].end());
        temp_v.insert(temp_v.end(), llXX_el7_values[k].begin(), llXX_el7_values[k].end());
        temp_v.insert(temp_v.end(), hhISR_el7_values[k].begin(), hhISR_el7_values[k].end());
        temp_v.insert(temp_v.end(), charged_el9_values[k].begin(), charged_el9_values[k].end());
        temp_v.insert(temp_v.end(), mixed_el9_values[k].begin(), mixed_el9_values[k].end());
        temp_v.insert(temp_v.end(), uubar_el9_values[k].begin(), uubar_el9_values[k].end());
        temp_v.insert(temp_v.end(), ddbar_el9_values[k].begin(), ddbar_el9_values[k].end());
        temp_v.insert(temp_v.end(), ssbar_el9_values[k].begin(), ssbar_el9_values[k].end());
        temp_v.insert(temp_v.end(), ccbar_el9_values[k].begin(), ccbar_el9_values[k].end());
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
        charged_el7_hist[k] = new TH1D("charged el7", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        mixed_el7_hist[k] = new TH1D("mixed el7", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        uubar_el7_hist[k] = new TH1D("u#bar{u} el7", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        ddbar_el7_hist[k] = new TH1D("d#bar{d} el7", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        ssbar_el7_hist[k] = new TH1D("s#bar{s} el7", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        ccbar_el7_hist[k] = new TH1D("c#bar{c} el7", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        taupair_el7_hist[k] = new TH1D("#tau#bar{#tau} el7", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        mumu_el7_hist[k] = new TH1D("#mu#bar{#mu} el7", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        gg_el7_hist[k] = new TH1D("gg el7", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        ee_el7_hist[k] = new TH1D("e#bar{e} el7", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        eeee_el7_hist[k] = new TH1D("ee#bar{e}#bar{e} el7", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        eemumu_el7_hist[k] = new TH1D("e#bar{e}#mu#bar{#mu} el7", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        llXX_el7_hist[k] = new TH1D("\ell#bar{\ell}XX el7", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        hhISR_el7_hist[k] = new TH1D("hhISR el7", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        charged_el9_hist[k] = new TH1D("charged el9", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        mixed_el9_hist[k] = new TH1D("mixed el9", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        uubar_el9_hist[k] = new TH1D("u#bar{u} el9", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        ddbar_el9_hist[k] = new TH1D("d#bar{d} el9", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        ssbar_el9_hist[k] = new TH1D("s#bar{s} el9", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        ccbar_el9_hist[k] = new TH1D("c#bar{c} el9", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        el7_hist[k] = new TH1D("el7 background", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        el9_hist[k] = new TH1D("el9 background", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        stat_error_hist[k] = new TH1D("MC stat error", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        data_hist[k] = new TH1D("data", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        Ratio_hist[k] = new TH1D((variable_names.at(k) + "_ratio").c_str(), ";;el9/el7", bins, min, max);
    }

    int index = std::find(variable_names.begin(), variable_names.end(), std::string("log_{10}SignalProbability")) - variable_names.begin();
    for (int i = 0; i < (int)charged_el7_values[index].size(); i++) charged_el7_values[index].at(i) = log10l(charged_el7_values[index].at(i));
    for (int i = 0; i < (int)mixed_el7_values[index].size(); i++) mixed_el7_values[index].at(i) = log10l(mixed_el7_values[index].at(i));
    for (int i = 0; i < (int)uubar_el7_values[index].size(); i++) uubar_el7_values[index].at(i) = log10l(uubar_el7_values[index].at(i));
    for (int i = 0; i < (int)ddbar_el7_values[index].size(); i++) ddbar_el7_values[index].at(i) = log10l(ddbar_el7_values[index].at(i));
    for (int i = 0; i < (int)ssbar_el7_values[index].size(); i++) ssbar_el7_values[index].at(i) = log10l(ssbar_el7_values[index].at(i));
    for (int i = 0; i < (int)ccbar_el7_values[index].size(); i++) ccbar_el7_values[index].at(i) = log10l(ccbar_el7_values[index].at(i));
    for (int i = 0; i < (int)taupair_el7_values[index].size(); i++) taupair_el7_values[index].at(i) = log10l(taupair_el7_values[index].at(i));
    for (int i = 0; i < (int)mumu_el7_values[index].size(); i++) mumu_el7_values[index].at(i) = log10l(mumu_el7_values[index].at(i));
    for (int i = 0; i < (int)gg_el7_values[index].size(); i++) gg_el7_values[index].at(i) = log10l(gg_el7_values[index].at(i));
    for (int i = 0; i < (int)ee_el7_values[index].size(); i++) ee_el7_values[index].at(i) = log10l(ee_el7_values[index].at(i));
    for (int i = 0; i < (int)eeee_el7_values[index].size(); i++) eeee_el7_values[index].at(i) = log10l(eeee_el7_values[index].at(i));
    for (int i = 0; i < (int)eemumu_el7_values[index].size(); i++) eemumu_el7_values[index].at(i) = log10l(eemumu_el7_values[index].at(i));
    for (int i = 0; i < (int)llXX_el7_values[index].size(); i++) llXX_el7_values[index].at(i) = log10l(llXX_el7_values[index].at(i));
    for (int i = 0; i < (int)hhISR_el7_values[index].size(); i++) hhISR_el7_values[index].at(i) = log10l(hhISR_el7_values[index].at(i));
    for (int i = 0; i < (int)charged_el9_values[index].size(); i++) charged_el9_values[index].at(i) = log10l(charged_el9_values[index].at(i));
    for (int i = 0; i < (int)mixed_el9_values[index].size(); i++) mixed_el9_values[index].at(i) = log10l(mixed_el9_values[index].at(i));
    for (int i = 0; i < (int)uubar_el9_values[index].size(); i++) uubar_el9_values[index].at(i) = log10l(uubar_el9_values[index].at(i));
    for (int i = 0; i < (int)ddbar_el9_values[index].size(); i++) ddbar_el9_values[index].at(i) = log10l(ddbar_el9_values[index].at(i));
    for (int i = 0; i < (int)ssbar_el9_values[index].size(); i++) ssbar_el9_values[index].at(i) = log10l(ssbar_el9_values[index].at(i));
    for (int i = 0; i < (int)ccbar_el9_values[index].size(); i++) ccbar_el9_values[index].at(i) = log10l(ccbar_el9_values[index].at(i));
    for (int i = 0; i < (int)el7_values[index].size(); i++) el7_values[index].at(i) = log10l(el7_values[index].at(i));
    for (int i = 0; i < (int)el9_values[index].size(); i++) el9_values[index].at(i) = log10l(el9_values[index].at(i));
    for (int i = 0; i < (int)MC_values[index].size(); i++) MC_values[index].at(i) = log10l(MC_values[index].at(i));
    for (int i = 0; i < (int)data_values[index].size(); i++) data_values[index].at(i) = log10l(data_values[index].at(i));

    for (int k = 0; k < (int)variable_names.size(); k++) { // fill
        for (int i = 0; i < (int)charged_el7_values[k].size(); i++) charged_el7_hist[k]->Fill(charged_el7_values[k].at(i), charged_el7_weights.at(i));
        for (int i = 0; i < (int)mixed_el7_values[k].size(); i++) mixed_el7_hist[k]->Fill(mixed_el7_values[k].at(i), mixed_el7_weights.at(i));
        for (int i = 0; i < (int)uubar_el7_values[k].size(); i++) uubar_el7_hist[k]->Fill(uubar_el7_values[k].at(i), uubar_el7_weights.at(i));
        for (int i = 0; i < (int)ddbar_el7_values[k].size(); i++) ddbar_el7_hist[k]->Fill(ddbar_el7_values[k].at(i), ddbar_el7_weights.at(i));
        for (int i = 0; i < (int)ssbar_el7_values[k].size(); i++) ssbar_el7_hist[k]->Fill(ssbar_el7_values[k].at(i), ssbar_el7_weights.at(i));
        for (int i = 0; i < (int)ccbar_el7_values[k].size(); i++) ccbar_el7_hist[k]->Fill(ccbar_el7_values[k].at(i), ccbar_el7_weights.at(i));
        for (int i = 0; i < (int)taupair_el7_values[k].size(); i++) taupair_el7_hist[k]->Fill(taupair_el7_values[k].at(i), taupair_el7_weights.at(i));
        for (int i = 0; i < (int)mumu_el7_values[k].size(); i++) mumu_el7_hist[k]->Fill(mumu_el7_values[k].at(i), mumu_el7_weights.at(i));
        for (int i = 0; i < (int)gg_el7_values[k].size(); i++) gg_el7_hist[k]->Fill(gg_el7_values[k].at(i), gg_el7_weights.at(i));
        for (int i = 0; i < (int)ee_el7_values[k].size(); i++) ee_el7_hist[k]->Fill(ee_el7_values[k].at(i), ee_el7_weights.at(i));
        for (int i = 0; i < (int)eeee_el7_values[k].size(); i++) eeee_el7_hist[k]->Fill(eeee_el7_values[k].at(i), eeee_el7_weights.at(i));
        for (int i = 0; i < (int)eemumu_el7_values[k].size(); i++) eemumu_el7_hist[k]->Fill(eemumu_el7_values[k].at(i), eemumu_el7_weights.at(i));
        for (int i = 0; i < (int)llXX_el7_values[k].size(); i++) llXX_el7_hist[k]->Fill(llXX_el7_values[k].at(i), llXX_el7_weights.at(i));
        for (int i = 0; i < (int)hhISR_el7_values[k].size(); i++) hhISR_el7_hist[k]->Fill(hhISR_el7_values[k].at(i), hhISR_el7_weights.at(i));
        for (int i = 0; i < (int)charged_el9_values[k].size(); i++) charged_el9_hist[k]->Fill(charged_el9_values[k].at(i), charged_el9_weights.at(i));
        for (int i = 0; i < (int)mixed_el9_values[k].size(); i++) mixed_el9_hist[k]->Fill(mixed_el9_values[k].at(i), mixed_el9_weights.at(i));
        for (int i = 0; i < (int)uubar_el9_values[k].size(); i++) uubar_el9_hist[k]->Fill(uubar_el9_values[k].at(i), uubar_el9_weights.at(i));
        for (int i = 0; i < (int)ddbar_el9_values[k].size(); i++) ddbar_el9_hist[k]->Fill(ddbar_el9_values[k].at(i), ddbar_el9_weights.at(i));
        for (int i = 0; i < (int)ssbar_el9_values[k].size(); i++) ssbar_el9_hist[k]->Fill(ssbar_el9_values[k].at(i), ssbar_el9_weights.at(i));
        for (int i = 0; i < (int)ccbar_el9_values[k].size(); i++) ccbar_el9_hist[k]->Fill(ccbar_el9_values[k].at(i), ccbar_el9_weights.at(i));
        for (int i = 0; i < (int)el7_values[k].size(); i++) el7_hist[k]->Fill(el7_values[k].at(i), el7_weights.at(i));
        for (int i = 0; i < (int)el9_values[k].size(); i++) el9_hist[k]->Fill(el9_values[k].at(i), el9_weights.at(i));
        for (int i = 0; i < (int)MC_values[k].size(); i++) stat_error_hist[k]->Fill(MC_values[k].at(i), weights.at(i));
        for (int i = 0; i < (int)data_values[k].size(); i++) data_hist[k]->Fill(data_values[k].at(i));
    }

    double el7_sum = 0;
    double el9_sum = 0;
    for (int i = 0; i < (int)el7_values[0].size(); i++) el7_sum = el7_sum + el7_weights.at(i);
    for (int i = 0; i < (int)el9_values[0].size(); i++) el9_sum = el9_sum + el9_weights.at(i);

    for (int k = 0; k < (int)variable_names.size(); k++) { // draw
        // Scale the histogram
        CAL_el7 = el9_sum / el7_sum;
        CAL_el9 = 1.0;
        el9_hist[k]->Scale(1.0);
        el7_hist[k]->Scale(CAL_el7);

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
        Ratio_hist[k]->Divide(el9_hist[k], el7_hist[k]);

        TCanvas* c_temp = new TCanvas("c", "", 800, 800); c_temp->cd();

        TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.35, 1.0, 1.0);
        pad1->SetBottomMargin(0.08); pad1->SetLeftMargin(0.15);
        pad1->SetGridx(); pad1->Draw(); pad1->cd();

        gStyle->SetPalette(kPastel);

        Float_t ymax_1 = el9_hist[k]->GetMaximum();
        Float_t ymax_2 = el7_hist[k]->GetMaximum();
        double real_max = 0;
        if (ymax_1 > ymax_2) real_max = ymax_1;
        else real_max = ymax_2;

        el9_hist[k]->SetMaximum(real_max * 1.1);

        el7_hist[k]->SetStats(0);
        el9_hist[k]->SetStats(0);

        el9_hist[k]->SetFillColor(4); el9_hist[k]->SetLineColor(4); el9_hist[k]->SetLineWidth(1); el9_hist[k]->SetFillStyle(3002); el9_hist[k]->Draw("e1 Hist");
        el7_hist[k]->SetFillColor(2); el7_hist[k]->SetLineColor(2); el7_hist[k]->SetLineWidth(1); el7_hist[k]->SetFillStyle(3003); el7_hist[k]->Draw("e1 Hist SAME");
        TLegend* legend = pad1->BuildLegend(0.95, 0.9, 0.75, 0.6);
        legend->SetFillStyle(0); legend->SetLineWidth(0);
        TPaveText* pt = new TPaveText(0.135, 0.88, 0.6, 1.0, "NDC NB"); pt->SetFillStyle(0); pt->SetLineWidth(0); pt->AddText(("el7 scaled to el9, el9/el7= " + std::to_string(CAL_el7/ CAL_el9)).c_str()); pt->Draw();

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
        c_temp->SaveAs((variable_names.at(k) + "_el7_vs_el9.png").c_str());

        delete c_temp;
    }

    // Print rd-ri discrepancy
    printf("data num: %ld\n", data_values[0].size());
    printf("el7 num with calibration: %lf\n", el7_sum);
    printf("el9 num with calibration: %lf\n", el9_sum);

    // free
    delete[] MC_values;
    delete[] charged_el7_values;
    delete[] mixed_el7_values;
    delete[] uubar_el7_values;
    delete[] ddbar_el7_values;
    delete[] ssbar_el7_values;
    delete[] ccbar_el7_values;
    delete[] taupair_el7_values;
    delete[] mumu_el7_values;
    delete[] gg_el7_values;
    delete[] ee_el7_values;
    delete[] eeee_el7_values;
    delete[] eemumu_el7_values;
    delete[] llXX_el7_values;
    delete[] hhISR_el7_values;
    delete[] charged_el9_values;
    delete[] mixed_el9_values;
    delete[] uubar_el9_values;
    delete[] ddbar_el9_values;
    delete[] ssbar_el9_values;
    delete[] ccbar_el9_values;
    delete[] el7_values;
    delete[] el9_values;

    delete[] data_values;

    for (int k = 0; k < Nvar_num; k++) {
        delete Stack[k];
        delete charged_el7_hist[k];
        delete mixed_el7_hist[k];
        delete uubar_el7_hist[k];
        delete ddbar_el7_hist[k];
        delete ssbar_el7_hist[k];
        delete ccbar_el7_hist[k];
        delete taupair_el7_hist[k];
        delete mumu_el7_hist[k];
        delete gg_el7_hist[k];
        delete ee_el7_hist[k];
        delete eeee_el7_hist[k];
        delete eemumu_el7_hist[k];
        delete llXX_el7_hist[k];
        delete hhISR_el7_hist[k];
        delete charged_el9_hist[k];
        delete mixed_el9_hist[k];
        delete uubar_el9_hist[k];
        delete ddbar_el9_hist[k];
        delete ssbar_el9_hist[k];
        delete ccbar_el9_hist[k];
        delete el7_hist[k];
        delete el9_hist[k];
        delete stat_error_hist[k];
        delete data_hist[k];
        delete Ratio_hist[k];
    }

    free(Stack);
    free(charged_el7_hist);
    free(mixed_el7_hist);
    free(uubar_el7_hist);
    free(ddbar_el7_hist);
    free(ssbar_el7_hist);
    free(ccbar_el7_hist);
    free(taupair_el7_hist);
    free(mumu_el7_hist);
    free(gg_el7_hist);
    free(ee_el7_hist);
    free(eeee_el7_hist);
    free(eemumu_el7_hist);
    free(llXX_el7_hist);
    free(hhISR_el7_hist);
    free(charged_el9_hist);
    free(mixed_el9_hist);
    free(uubar_el9_hist);
    free(ddbar_el9_hist);
    free(ssbar_el9_hist);
    free(ccbar_el9_hist);
    free(el7_hist);
    free(el9_hist);
    free(stat_error_hist);
    free(data_hist);
    free(Ratio_hist);
}
