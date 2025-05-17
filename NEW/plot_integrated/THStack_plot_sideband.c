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
#include "template.h"

void THStack_plot_sideband() {

    NormalizeAtEachMXs = true;

    const char* Sideband_MC_CHG_validation_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/SatoriRD_LS_MC_side/CHG_analysis/validation_v004/final_output";
    const char* Sideband_MC_MIX_validation_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/SatoriRD_LS_MC_side/MIX_analysis/validation_v004/final_output";
    const char* Sideband_MC_UUBAR_validation_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/SatoriRD_LS_MC_side/UUBAR_analysis/validation_v004/final_output";
    const char* Sideband_MC_DDBAR_validation_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/SatoriRD_LS_MC_side/DDBAR_analysis/validation_v004/final_output";
    const char* Sideband_MC_SSBAR_validation_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/SatoriRD_LS_MC_side/SSBAR_analysis/validation_v004/final_output";
    const char* Sideband_MC_CHARM_validation_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/SatoriRD_LS_MC_side/CHARM_analysis/validation_v004/final_output";

    const char* Sideband_data_dirname = "/home/belle2/junewoo/storage_ghi/Analysis/SatoriRD_LS_data_side/SIGNAL_analysis/validation_v004/final_output";

    // Count event to normalize at each MXs region
    LetsCountMC(Sideband_MC_CHG_validation_dirname, "CHG");
    LetsCountMC(Sideband_MC_MIX_validation_dirname, "MIX");
    LetsCountMC(Sideband_MC_UUBAR_validation_dirname, "UUBAR");
    LetsCountMC(Sideband_MC_DDBAR_validation_dirname, "DDBAR");
    LetsCountMC(Sideband_MC_SSBAR_validation_dirname, "SSBAR");
    LetsCountMC(Sideband_MC_CHARM_validation_dirname, "CHARM");

    LetsCountdata(Sideband_data_dirname);

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

    std::vector<double>* Sideband_MC_values = new std::vector<double>[Nvar_num];
    std::vector<double>* charged_values = new std::vector<double>[Nvar_num];
    std::vector<double>* mixed_values = new std::vector<double>[Nvar_num];
    std::vector<double>* uubar_values = new std::vector<double>[Nvar_num];
    std::vector<double>* ddbar_values = new std::vector<double>[Nvar_num];
    std::vector<double>* ssbar_values = new std::vector<double>[Nvar_num];
    std::vector<double>* ccbar_values = new std::vector<double>[Nvar_num];
    std::vector<double>* taupair_values = new std::vector<double>[Nvar_num];
    std::vector<double>* mumu_values = new std::vector<double>[Nvar_num];
    std::vector<double>* gg_values = new std::vector<double>[Nvar_num];
    std::vector<double>* ee_values = new std::vector<double>[Nvar_num];
    std::vector<double>* eeee_values = new std::vector<double>[Nvar_num];
    std::vector<double>* eemumu_values = new std::vector<double>[Nvar_num];
    std::vector<double>* llXX_values = new std::vector<double>[Nvar_num];
    std::vector<double>* hhISR_values = new std::vector<double>[Nvar_num];
    std::vector<int> Sideband_MC_numbering;

    std::vector<double>* Sideband_data_values = new std::vector<double>[Nvar_num];

    std::vector<double> weights;
    std::vector<double> charged_weights;
    std::vector<double> mixed_weights;
    std::vector<double> uubar_weights;
    std::vector<double> ddbar_weights;
    std::vector<double> ssbar_weights;
    std::vector<double> ccbar_weights;
    std::vector<double> taupair_weights;
    std::vector<double> mumu_weights;
    std::vector<double> gg_weights;
    std::vector<double> ee_weights;
    std::vector<double> eeee_weights;
    std::vector<double> eemumu_weights;
    std::vector<double> llXX_weights;
    std::vector<double> hhISR_weights;
    /*
    LetsFillMC(Sideband_MC_CHG_validation_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "CHG");
    LetsFillMC(Sideband_MC_MIX_validation_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "MIX");
    LetsFillMC_correction(Sideband_MC_UUBAR_validation_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "UUBAR", NormFactor_UUBAR);
    LetsFillMC_correction(Sideband_MC_DDBAR_validation_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "DDBAR", NormFactor_DDBAR);
    LetsFillMC_correction(Sideband_MC_SSBAR_validation_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "SSBAR", NormFactor_SSBAR);
    LetsFillMC_correction(Sideband_MC_CHARM_validation_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "CHARM", NormFactor_CHARM);
    */
    LetsFillMC(Sideband_MC_CHG_validation_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "CHG");
    LetsFillMC(Sideband_MC_MIX_validation_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "MIX");
    LetsFillMC(Sideband_MC_UUBAR_validation_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "UUBAR");
    LetsFillMC(Sideband_MC_DDBAR_validation_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "DDBAR");
    LetsFillMC(Sideband_MC_SSBAR_validation_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "SSBAR");
    LetsFillMC(Sideband_MC_CHARM_validation_dirname, variable_names, branch_names, Sideband_MC_values, &Sideband_MC_numbering, &weights, "CHARM");
    
    LetsFilldata(Sideband_data_dirname, variable_names, branch_names, Sideband_data_values);

    // sort variables
    for (int k = 0; k < (int)Sideband_MC_numbering.size(); k++) {
        if (Sideband_MC_numbering.at(k) == static_cast<int>(MCsample::CHG)) {
            for (int l = 0; l < (int)variable_names.size(); l++) charged_values[l].push_back(Sideband_MC_values[l].at(k));
            charged_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == static_cast<int>(MCsample::MIX)) {
            for (int l = 0; l < (int)variable_names.size(); l++) mixed_values[l].push_back(Sideband_MC_values[l].at(k));
            mixed_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == static_cast<int>(MCsample::UUBAR)) {
            for (int l = 0; l < (int)variable_names.size(); l++) uubar_values[l].push_back(Sideband_MC_values[l].at(k));
            uubar_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == static_cast<int>(MCsample::DDBAR)) {
            for (int l = 0; l < (int)variable_names.size(); l++) ddbar_values[l].push_back(Sideband_MC_values[l].at(k));
            ddbar_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == static_cast<int>(MCsample::SSBAR)) {
            for (int l = 0; l < (int)variable_names.size(); l++) ssbar_values[l].push_back(Sideband_MC_values[l].at(k));
            ssbar_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == static_cast<int>(MCsample::CHARM)) {
            for (int l = 0; l < (int)variable_names.size(); l++) ccbar_values[l].push_back(Sideband_MC_values[l].at(k));
            ccbar_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == static_cast<int>(MCsample::TAU)) {
            for (int l = 0; l < (int)variable_names.size(); l++) taupair_values[l].push_back(Sideband_MC_values[l].at(k));
            taupair_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == static_cast<int>(MCsample::MUMU)) {
            for (int l = 0; l < (int)variable_names.size(); l++) mumu_values[l].push_back(Sideband_MC_values[l].at(k));
            mumu_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == static_cast<int>(MCsample::GG)) {
            for (int l = 0; l < (int)variable_names.size(); l++) gg_values[l].push_back(Sideband_MC_values[l].at(k));
            gg_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == static_cast<int>(MCsample::EE)) {
            for (int l = 0; l < (int)variable_names.size(); l++) ee_values[l].push_back(Sideband_MC_values[l].at(k));
            ee_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == static_cast<int>(MCsample::EEEE)) {
            for (int l = 0; l < (int)variable_names.size(); l++) eeee_values[l].push_back(Sideband_MC_values[l].at(k));
            eeee_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == static_cast<int>(MCsample::EEMUMU)) {
            for (int l = 0; l < (int)variable_names.size(); l++) eemumu_values[l].push_back(Sideband_MC_values[l].at(k));
            eemumu_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == static_cast<int>(MCsample::LLXX)) {
            for (int l = 0; l < (int)variable_names.size(); l++) llXX_values[l].push_back(Sideband_MC_values[l].at(k));
            llXX_weights.push_back(weights.at(k));
        }
        else if (Sideband_MC_numbering.at(k) == static_cast<int>(MCsample::HHISR)) {
            for (int l = 0; l < (int)variable_names.size(); l++) hhISR_values[l].push_back(Sideband_MC_values[l].at(k));
            hhISR_weights.push_back(weights.at(k));
        }
        else {
            printf("undefined numbering!\n");
            exit(1);
        }
    }

    THStack** Stack = (THStack**)malloc(sizeof(THStack*) * Nvar_num);
    TH1D** charged_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** mixed_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** uubar_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** ddbar_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** ssbar_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** ccbar_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** taupair_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** mumu_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** gg_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** ee_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** eeee_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** eemumu_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** llXX_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** hhISR_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** stat_error_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** data_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);
    TH1D** Ratio_hist = (TH1D**)malloc(sizeof(TH1D*) * Nvar_num);

    for (int k = 0; k < (int)variable_names.size(); k++) { // malloc TH1D
        std::vector<double> temp_v;
        temp_v.insert(temp_v.end(), charged_values[k].begin(), charged_values[k].end());
        temp_v.insert(temp_v.end(), mixed_values[k].begin(), mixed_values[k].end());
        temp_v.insert(temp_v.end(), uubar_values[k].begin(), uubar_values[k].end());
        temp_v.insert(temp_v.end(), ddbar_values[k].begin(), ddbar_values[k].end());
        temp_v.insert(temp_v.end(), ssbar_values[k].begin(), ssbar_values[k].end());
        temp_v.insert(temp_v.end(), ccbar_values[k].begin(), ccbar_values[k].end());
        temp_v.insert(temp_v.end(), taupair_values[k].begin(), taupair_values[k].end());
        temp_v.insert(temp_v.end(), mumu_values[k].begin(), mumu_values[k].end());
        temp_v.insert(temp_v.end(), gg_values[k].begin(), gg_values[k].end());
        temp_v.insert(temp_v.end(), ee_values[k].begin(), ee_values[k].end());
        temp_v.insert(temp_v.end(), eeee_values[k].begin(), eeee_values[k].end());
        temp_v.insert(temp_v.end(), eemumu_values[k].begin(), eemumu_values[k].end());
        temp_v.insert(temp_v.end(), llXX_values[k].begin(), llXX_values[k].end());
        temp_v.insert(temp_v.end(), hhISR_values[k].begin(), hhISR_values[k].end());
        temp_v.insert(temp_v.end(), Sideband_data_values[k].begin(), Sideband_data_values[k].end());


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
        charged_hist[k] = new TH1D("charged", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        mixed_hist[k] = new TH1D("mixed", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        uubar_hist[k] = new TH1D("u#bar{u}", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        ddbar_hist[k] = new TH1D("d#bar{d}", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        ssbar_hist[k] = new TH1D("s#bar{s}", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        ccbar_hist[k] = new TH1D("c#bar{c}", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        taupair_hist[k] = new TH1D("#tau#bar{#tau}", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        mumu_hist[k] = new TH1D("#mu#bar{#mu}", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        gg_hist[k] = new TH1D("gg", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        ee_hist[k] = new TH1D("e#bar{e}", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        eeee_hist[k] = new TH1D("ee#bar{e}#bar{e}", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        eemumu_hist[k] = new TH1D("e#bar{e}#mu#bar{#mu}", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        llXX_hist[k] = new TH1D("\ell#bar{\ell}XX", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        hhISR_hist[k] = new TH1D("hhISR", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        stat_error_hist[k] = new TH1D("MC stat error", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        data_hist[k] = new TH1D("data", (";" + variable_names.at(k) + ";number of candidates").c_str(), bins, min, max);
        Ratio_hist[k] = new TH1D((variable_names.at(k) + "_ratio").c_str(), ";;data/MC", bins, min, max);
    }

    int index = std::find(variable_names.begin(), variable_names.end(), std::string("log_{10}SignalProbability")) - variable_names.begin();
    for (int i = 0; i < (int)charged_values[index].size(); i++) charged_values[index].at(i) = log10l(charged_values[index].at(i));
    for (int i = 0; i < (int)mixed_values[index].size(); i++) mixed_values[index].at(i) = log10l(mixed_values[index].at(i));
    for (int i = 0; i < (int)uubar_values[index].size(); i++) uubar_values[index].at(i) = log10l(uubar_values[index].at(i));
    for (int i = 0; i < (int)ddbar_values[index].size(); i++) ddbar_values[index].at(i) = log10l(ddbar_values[index].at(i));
    for (int i = 0; i < (int)ssbar_values[index].size(); i++) ssbar_values[index].at(i) = log10l(ssbar_values[index].at(i));
    for (int i = 0; i < (int)ccbar_values[index].size(); i++) ccbar_values[index].at(i) = log10l(ccbar_values[index].at(i));
    for (int i = 0; i < (int)taupair_values[index].size(); i++) taupair_values[index].at(i) = log10l(taupair_values[index].at(i));
    for (int i = 0; i < (int)mumu_values[index].size(); i++) mumu_values[index].at(i) = log10l(mumu_values[index].at(i));
    for (int i = 0; i < (int)gg_values[index].size(); i++) gg_values[index].at(i) = log10l(gg_values[index].at(i));
    for (int i = 0; i < (int)ee_values[index].size(); i++) ee_values[index].at(i) = log10l(ee_values[index].at(i));
    for (int i = 0; i < (int)eeee_values[index].size(); i++) eeee_values[index].at(i) = log10l(eeee_values[index].at(i));
    for (int i = 0; i < (int)eemumu_values[index].size(); i++) eemumu_values[index].at(i) = log10l(eemumu_values[index].at(i));
    for (int i = 0; i < (int)llXX_values[index].size(); i++) llXX_values[index].at(i) = log10l(llXX_values[index].at(i));
    for (int i = 0; i < (int)hhISR_values[index].size(); i++) hhISR_values[index].at(i) = log10l(hhISR_values[index].at(i));
    for (int i = 0; i < (int)Sideband_MC_values[index].size(); i++) Sideband_MC_values[index].at(i) = log10l(Sideband_MC_values[index].at(i));
    for (int i = 0; i < (int)Sideband_data_values[index].size(); i++) Sideband_data_values[index].at(i) = log10l(Sideband_data_values[index].at(i));

    for (int k = 0; k < (int)variable_names.size(); k++) { // fill
        for (int i = 0; i < (int)charged_values[k].size(); i++) charged_hist[k]->Fill(charged_values[k].at(i), charged_weights.at(i));
        for (int i = 0; i < (int)mixed_values[k].size(); i++) mixed_hist[k]->Fill(mixed_values[k].at(i), mixed_weights.at(i));
        for (int i = 0; i < (int)uubar_values[k].size(); i++) uubar_hist[k]->Fill(uubar_values[k].at(i), uubar_weights.at(i));
        for (int i = 0; i < (int)ddbar_values[k].size(); i++) ddbar_hist[k]->Fill(ddbar_values[k].at(i), ddbar_weights.at(i));
        for (int i = 0; i < (int)ssbar_values[k].size(); i++) ssbar_hist[k]->Fill(ssbar_values[k].at(i), ssbar_weights.at(i));
        for (int i = 0; i < (int)ccbar_values[k].size(); i++) ccbar_hist[k]->Fill(ccbar_values[k].at(i), ccbar_weights.at(i));
        for (int i = 0; i < (int)taupair_values[k].size(); i++) taupair_hist[k]->Fill(taupair_values[k].at(i), taupair_weights.at(i));
        for (int i = 0; i < (int)mumu_values[k].size(); i++) mumu_hist[k]->Fill(mumu_values[k].at(i), mumu_weights.at(i));
        for (int i = 0; i < (int)gg_values[k].size(); i++) gg_hist[k]->Fill(gg_values[k].at(i), gg_weights.at(i));
        for (int i = 0; i < (int)ee_values[k].size(); i++) ee_hist[k]->Fill(ee_values[k].at(i), ee_weights.at(i));
        for (int i = 0; i < (int)eeee_values[k].size(); i++) eeee_hist[k]->Fill(eeee_values[k].at(i), eeee_weights.at(i));
        for (int i = 0; i < (int)eemumu_values[k].size(); i++) eemumu_hist[k]->Fill(eemumu_values[k].at(i), eemumu_weights.at(i));
        for (int i = 0; i < (int)llXX_values[k].size(); i++) llXX_hist[k]->Fill(llXX_values[k].at(i), llXX_weights.at(i));
        for (int i = 0; i < (int)hhISR_values[k].size(); i++) hhISR_hist[k]->Fill(hhISR_values[k].at(i), hhISR_weights.at(i));
        for (int i = 0; i < (int)Sideband_MC_values[k].size(); i++) stat_error_hist[k]->Fill(Sideband_MC_values[k].at(i), weights.at(i));
        for (int i = 0; i < (int)Sideband_data_values[k].size(); i++) data_hist[k]->Fill(Sideband_data_values[k].at(i));
    }
    if (NormalizeAtEachMXs == false) {
        Ratio_one_bin->Divide(data_one_bin, MC_one_bin);
    }
    else {
        Ratio_Nevt_MXs1->Divide(data_Nevt_MXs1, MC_Nevt_MXs1);
        Ratio_Nevt_MXs2->Divide(data_Nevt_MXs2, MC_Nevt_MXs2);
        Ratio_Nevt_MXs3->Divide(data_Nevt_MXs3, MC_Nevt_MXs3);
    }

    printf("charged: %d\n", (int)charged_values[0].size());
    printf("mixed: %d\n", (int)mixed_values[0].size());
    printf("uubar: %d\n", (int)uubar_values[0].size());
    printf("ddbar: %d\n", (int)ddbar_values[0].size());
    printf("ssbar: %d\n", (int)ssbar_values[0].size());
    printf("ccbar: %d\n", (int)ccbar_values[0].size());
    printf("taupair: %d\n", (int)taupair_values[0].size());
    printf("mumu: %d\n", (int)mumu_values[0].size());
    printf("gg: %d\n", (int)gg_values[0].size());
    printf("ee: %d\n", (int)ee_values[0].size());
    printf("eeee: %d\n", (int)eeee_values[0].size());
    printf("eemuu: %d\n", (int)eemumu_values[0].size());
    printf("llXX: %d\n", (int)llXX_values[0].size());
    printf("hhISR: %d\n", (int)hhISR_values[0].size());
    printf("data: %d\n", (int)Sideband_data_values[0].size());

    for (int k = 0; k < (int)variable_names.size(); k++) { // draw
        // Scale the histogram if `NormalizeAtEachMXs` is not turned on
        if (NormalizeAtEachMXs == false) {
            CAL = Ratio_one_bin->GetBinContent(1);
            charged_hist[k]->Scale(CAL);
            mixed_hist[k]->Scale(CAL);
            uubar_hist[k]->Scale(CAL);
            ddbar_hist[k]->Scale(CAL);
            ssbar_hist[k]->Scale(CAL);
            ccbar_hist[k]->Scale(CAL);
            stat_error_hist[k]->Scale(CAL);
        }

        Stack[k]->Add(charged_hist[k]);
        Stack[k]->Add(mixed_hist[k]);
        Stack[k]->Add(uubar_hist[k]);
        Stack[k]->Add(ddbar_hist[k]);
        Stack[k]->Add(ssbar_hist[k]);
        Stack[k]->Add(ccbar_hist[k]);
        //Stack[k]->Add(taupair_hist[k]);
        //Stack[k]->Add(mumu_hist[k]);
        //Stack[k]->Add(gg_hist[k]);
        //Stack[k]->Add(ee_hist[k]);
        //Stack[k]->Add(eeee_hist[k]);
        //Stack[k]->Add(eemumu_hist[k]);
        //Stack[k]->Add(llXX_hist[k]);
        //Stack[k]->Add(hhISR_hist[k]);

        Ratio_hist[k]->SetLineColor(kBlack); Ratio_hist[k]->SetMarkerStyle(21); Ratio_hist[k]->Sumw2(); Ratio_hist[k]->SetStats(0);
        Ratio_hist[k]->Divide(data_hist[k], stat_error_hist[k]);

        TCanvas* c_temp = new TCanvas("c", "", 800, 800); c_temp->cd();

        TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
        pad1->SetBottomMargin(0.13); pad1->SetLeftMargin(0.15);
        pad1->SetGridx(); pad1->Draw(); pad1->cd();

        if((variable_names.at(k).find("MVA") != std::string::npos) && (Sideband_MC_values[k].size() > 10000)) pad1->SetLogy(1);
        else pad1->SetLogy(0);

        gStyle->SetPalette(kPastel);

        Float_t ymax_1 = Stack[k]->GetMaximum();
        Float_t ymax_2 = data_hist[k]->GetMaximum();
        double real_max = 0;
        if (ymax_1 > ymax_2) real_max = ymax_1;
        else real_max = ymax_2;

        Stack[k]->SetMaximum(real_max * 1.2);

        Stack[k]->Draw("pfc Hist");
        Stack[k]->GetXaxis()->SetTitleSize(0.05);
        Stack[k]->GetYaxis()->SetTitleSize(0.05);
        Stack[k]->GetXaxis()->SetLabelSize(0.05);
        Stack[k]->GetYaxis()->SetLabelSize(0.05);
        stat_error_hist[k]->SetFillColor(12); stat_error_hist[k]->SetLineWidth(0); stat_error_hist[k]->SetFillStyle(3004); stat_error_hist[k]->Draw("e2 SAME");
        data_hist[k]->SetLineWidth(2); data_hist[k]->SetLineColor(kBlack); data_hist[k]->SetMarkerStyle(8); data_hist[k]->Draw("SAME eP EX0");
        TPaveText* pt_lumi = new TPaveText(0.18, 0.83, 0.48, 0.85, "NDC NB");
        pt_lumi->SetTextSize(0.035); pt_lumi->SetFillStyle(0); pt_lumi->SetLineWidth(0); pt_lumi->SetTextAlign(11); pt_lumi->AddText("Belle II #int L dt = 365.4 fb^{-1}"); pt_lumi->Draw();
        TLegend* legend = pad1->BuildLegend(0.92, 0.9, 0.67, 0.5);
        legend->SetFillStyle(0); legend->SetLineWidth(0);
        if (NormalizeAtEachMXs == false) {
            TPaveText* pt = new TPaveText(0.25, 0.88, 0.65, 1.0, "NDC NB");
            pt->SetFillStyle(0);
            pt->SetLineWidth(0);
            pt->SetTextSize(0.06);
            pt->AddText(("MC scaled to data, Data/MC= " + to_string_precision(CAL)).c_str());
            pt->Draw();
        }
        else {
            TPaveText* pt = new TPaveText(0.15, 0.88, 0.915, 1.0, "NDC NB");
            pt->SetFillStyle(0);
            pt->SetLineWidth(0);
            pt->SetTextSize(0.035);
            pt->AddText(("MC scaled to data at each MXs region, Data/MC= " + to_string_precision(Ratio_Nevt_MXs1->GetBinContent(1)) + ", " + to_string_precision(Ratio_Nevt_MXs2->GetBinContent(1)) + ", " + to_string_precision(Ratio_Nevt_MXs3->GetBinContent(1)) ).c_str());
            pt->Draw();
        }

        c_temp->cd();
        TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1, 0.3); pad2->SetTopMargin(0.03); pad2->SetBottomMargin(0.15); pad2->SetLeftMargin(0.15); pad2->SetGridx(); pad2->Draw(); pad2->cd();
        Ratio_hist[k]->SetMinimum(0.5); Ratio_hist[k]->SetMaximum(1.5); Ratio_hist[k]->SetLineWidth(2);
        Ratio_hist[k]->GetYaxis()->SetTitleSize(0.1); Ratio_hist[k]->GetYaxis()->SetTitleOffset(0.5);
        Ratio_hist[k]->GetXaxis()->SetLabelSize(0.1); Ratio_hist[k]->GetYaxis()->SetLabelSize(0.1);
        Ratio_hist[k]->Draw("e0p");
        TLine* line = new TLine(Ratio_hist[k]->GetXaxis()->GetXmin(), 1.0, Ratio_hist[k]->GetXaxis()->GetXmax(), 1.0);
        line->SetLineColor(kRed);
        line->SetLineStyle(1); line->SetLineWidth(3);
        line->Draw();

        c_temp->SetBottomMargin(0.0);
        c_temp->SaveAs((variable_names.at(k) + "_sideband.png").c_str());

        // print data/MC if it is bin-index
        if (variable_names.at(k) == "bin index") {
            PrintDataMCRatio(Stack[k], data_hist[k], Ratio_hist[k], "dataMCratio_sideband.txt");
            PrintDataMC(Stack[k], data_hist[k]);
            Printchi2(Stack[k], data_hist[k]);
        }

        delete c_temp;
    }

    // Print data-MC discrepancy
    if (NormalizeAtEachMXs == false) {
        double MC_sum = 0;
        for (int i = 0; i < (int)Sideband_MC_values[0].size(); i++) MC_sum = MC_sum + weights.at(i);
        printf("data num: %ld\n", Sideband_data_values[0].size());
        printf("MC num with calibration: %lf\n", MC_sum);
        printf("MC with calibration: %lf +- %lf\n", MC_one_bin->GetBinContent(1), MC_one_bin->GetBinError(1));
        printf("data with calibration: %lf +- %lf\n", data_one_bin->GetBinContent(1), data_one_bin->GetBinError(1));
        printf("data/MC with calibration: %lf +- %lf\n", Ratio_one_bin->GetBinContent(1), Ratio_one_bin->GetBinError(1));
    }
    else {
        printf("data num in 1st region: %lf +- %lf\n", data_Nevt_MXs1->GetBinContent(1), data_Nevt_MXs1->GetBinError(1));
        printf("data num in 2nd region: %lf +- %lf\n", data_Nevt_MXs2->GetBinContent(1), data_Nevt_MXs2->GetBinError(1));
        printf("data num in 3rd region: %lf +- %lf\n", data_Nevt_MXs3->GetBinContent(1), data_Nevt_MXs3->GetBinError(1));

        printf("MC num in 1st region: %lf +- %lf\n", MC_Nevt_MXs1->GetBinContent(1), MC_Nevt_MXs1->GetBinError(1));
        printf("MC num in 2nd region: %lf +- %lf\n", MC_Nevt_MXs2->GetBinContent(1), MC_Nevt_MXs2->GetBinError(1));
        printf("MC num in 3rd region: %lf +- %lf\n", MC_Nevt_MXs3->GetBinContent(1), MC_Nevt_MXs3->GetBinError(1));

        printf("data/MC in 1st region: %lf +- %lf\n", Ratio_Nevt_MXs1->GetBinContent(1), Ratio_Nevt_MXs1->GetBinError(1));
        printf("data/MC in 2nd region: %lf +- %lf\n", Ratio_Nevt_MXs2->GetBinContent(1), Ratio_Nevt_MXs2->GetBinError(1));
        printf("data/MC in 3rd region: %lf +- %lf\n", Ratio_Nevt_MXs3->GetBinContent(1), Ratio_Nevt_MXs3->GetBinError(1));
    }

    // free
    delete[] Sideband_MC_values;
    delete[] charged_values;
    delete[] mixed_values;
    delete[] uubar_values;
    delete[] ddbar_values;
    delete[] ssbar_values;
    delete[] ccbar_values;
    delete[] taupair_values;
    delete[] mumu_values;
    delete[] gg_values;
    delete[] ee_values;
    delete[] eeee_values;
    delete[] eemumu_values;
    delete[] llXX_values;
    delete[] hhISR_values;

    delete[]  Sideband_data_values;

    for (int k = 0; k < Nvar_num; k++) {
        delete Stack[k];
        delete charged_hist[k];
        delete mixed_hist[k];
        delete uubar_hist[k];
        delete ddbar_hist[k];
        delete ssbar_hist[k];
        delete ccbar_hist[k];
        delete taupair_hist[k];
        delete mumu_hist[k];
        delete gg_hist[k];
        delete ee_hist[k];
        delete eeee_hist[k];
        delete eemumu_hist[k];
        delete llXX_hist[k];
        delete hhISR_hist[k];
        delete stat_error_hist[k];
        delete data_hist[k];
        delete Ratio_hist[k];
    }

    free(Stack);
    free(charged_hist);
    free(mixed_hist);
    free(uubar_hist);
    free(ddbar_hist);
    free(ssbar_hist);
    free(ccbar_hist);
    free(taupair_hist);
    free(mumu_hist);
    free(gg_hist);
    free(ee_hist);
    free(eeee_hist);
    free(eemumu_hist);
    free(llXX_hist);
    free(hhISR_hist);
    free(stat_error_hist);
    free(data_hist);
    free(Ratio_hist);

    delete MC_one_bin;
    delete data_one_bin;
    delete Ratio_one_bin;
}
