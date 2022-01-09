/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This macro provides a simple example on how to use the trained classifiers
/// within an analysis module
/// - Project   : TMVA - a Root-integrated toolkit for multivariate data analysis
/// - Package   : TMVA
/// - Exectuable: TMVAClassificationApplication
///
/// \macro_output
/// \macro_code
/// \author Andreas Hoecker

# define N_Needed_info 37
# define N_event_info 15
# define N_Upsilon_info 14
# define N_Bsig_info 28
# define N_Btag_info 7
# define N_decay 38 // five decay mode + others

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;

void TMVAClassificationApplication(const char* filename, const char* option)
{
    std::string OPTION = std::string(option);
    if (OPTION == std::string("large")) {}
    else if (OPTION == std::string("small")) {}
    else {
        printf("unknown option!\n");
        exit(0);
    }

    std::string string_filename(filename);
    std::string OnlyFileName = string_filename.substr(string_filename.find_last_of("\\/") + 1, string_filename.size() - string_filename.find_last_of("\\/"));

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;

   // Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   float temp_EventDataToTree_f[N_event_info / 3];
   float temp_UpsilonDataToTree_f[N_Upsilon_info];
   float temp_BsigDataToTree_f[N_Bsig_info];
   float temp_BtagDataToTree_f[N_Btag_info];
   float temp_DataToTree_f[N_Needed_info];
   //float temp_DecayDataToTree_f[N_decay];
   float temp_Upsilon_decayIDToTree_f;
   float temp_Bsig_decayIDToTree_f;
   float temp_flag_f;

   int temp_EventDataToTree[N_event_info / 3];
   double temp_UpsilonDataToTree[N_Upsilon_info];
   double temp_BsigDataToTree[N_Bsig_info];
   double temp_BtagDataToTree[N_Btag_info];
   double temp_DataToTree[N_Needed_info];
   double temp_DecayDataToTree[N_decay];
   double temp_Upsilon_decayIDToTree;
   double temp_Bsig_decayIDToTree;
   int temp_flag;
   float Output_BB;
   float Output_Continuum;

   reader->AddSpectator("Btag_R2", &temp_DataToTree_f[0]);
   reader->AddSpectator("Btag_thrustBm", &temp_DataToTree_f[1]);
   reader->AddVariable("Btag_thrustOm", &temp_DataToTree_f[2]);
   reader->AddVariable("Btag_cosTBTO", &temp_DataToTree_f[3]);
   reader->AddSpectator("Btag_cosTBz", &temp_DataToTree_f[4]);
   reader->AddSpectator("Btag_KSFWVariables_et", &temp_DataToTree_f[5]);
   reader->AddSpectator("Btag_KSFWVariables_mm2", &temp_DataToTree_f[6]);
   reader->AddVariable("Btag_KSFWVariables_hso00", &temp_DataToTree_f[7]);
   reader->AddSpectator("Btag_KSFWVariables_hso01", &temp_DataToTree_f[8]);
   reader->AddSpectator("Btag_KSFWVariables_hso02", &temp_DataToTree_f[9]);
   reader->AddSpectator("Btag_KSFWVariables_hso03", &temp_DataToTree_f[10]);
   reader->AddSpectator("Btag_KSFWVariables_hso04", &temp_DataToTree_f[11]);
   reader->AddVariable("Btag_KSFWVariables_hso10", &temp_DataToTree_f[12]);
   reader->AddSpectator("Btag_KSFWVariables_hso12", &temp_DataToTree_f[13]);
   reader->AddSpectator("Btag_KSFWVariables_hso14", &temp_DataToTree_f[14]);
   reader->AddVariable("Btag_KSFWVariables_hso20", &temp_DataToTree_f[15]);
   reader->AddSpectator("Btag_KSFWVariables_hso22", &temp_DataToTree_f[16]);
   reader->AddSpectator("Btag_KSFWVariables_hso24", &temp_DataToTree_f[17]);
   reader->AddSpectator("Btag_KSFWVariables_hoo0", &temp_DataToTree_f[18]);
   reader->AddSpectator("Btag_KSFWVariables_hoo1", &temp_DataToTree_f[19]);
   reader->AddSpectator("Btag_KSFWVariables_hoo2", &temp_DataToTree_f[20]);
   reader->AddSpectator("Btag_KSFWVariables_hoo3", &temp_DataToTree_f[21]);
   reader->AddSpectator("Btag_KSFWVariables_hoo4", &temp_DataToTree_f[22]);
   reader->AddSpectator("Btag_CleoConeCS_1", &temp_DataToTree_f[23]);
   reader->AddSpectator("Btag_CleoConeCS_2", &temp_DataToTree_f[24]);
   reader->AddSpectator("Btag_CleoConeCS_3", &temp_DataToTree_f[25]);
   reader->AddSpectator("Btag_CleoConeCS_4", &temp_DataToTree_f[26]);
   reader->AddSpectator("Btag_CleoConeCS_5", &temp_DataToTree_f[27]);
   reader->AddSpectator("Btag_CleoConeCS_6", &temp_DataToTree_f[28]);
   reader->AddSpectator("Btag_CleoConeCS_7", &temp_DataToTree_f[29]);
   reader->AddSpectator("Btag_CleoConeCS_8", &temp_DataToTree_f[30]);
   reader->AddSpectator("Btag_CleoConeCS_9", &temp_DataToTree_f[31]);
   reader->AddVariable("Btag_useCMSFrame_theta", &temp_DataToTree_f[34]);
   reader->AddVariable("Btag_extraInfo_SignalProbability", &temp_BtagDataToTree_f[5]);
   reader->AddVariable("Btag_chiProb", &temp_BtagDataToTree_f[6]);
   //reader->AddSpectator("extraInfo__boDecayHash__bc", &temp_DataToTree_f[35]);
   //reader->AddSpectator("extraInfo__boDecayHashExtended__bc", &temp_DataToTree_f[36]);

   reader->AddSpectator("__experiment__", &temp_EventDataToTree_f[0]);
   reader->AddSpectator("__run__", &temp_EventDataToTree_f[1]);
   reader->AddSpectator("__event__", &temp_EventDataToTree_f[2]);
   reader->AddSpectator("__candidate__", &temp_EventDataToTree_f[3]);
   reader->AddSpectator("__ncandidates__", &temp_EventDataToTree_f[4]);
   reader->AddSpectator("extraInfo__bodecayModeID__bc", &temp_Upsilon_decayIDToTree_f);
   reader->AddSpectator("Bsig_daughter_0_extraInfo_decayModeID", &temp_Bsig_decayIDToTree_f);
   reader->AddSpectator("nROE_ECLClusters__bocleanMask__bc", &temp_UpsilonDataToTree_f[0]);
   reader->AddSpectator("nROE_KLMClusters", &temp_UpsilonDataToTree_f[1]);
   reader->AddSpectator("nROE_Tracks__bocleanMask__bc", &temp_UpsilonDataToTree_f[2]);
   reader->AddSpectator("roeEextra__bocleanMask__bc", &temp_UpsilonDataToTree_f[3]);
   reader->AddSpectator("nROE_NeutralECLClusters__bocleanMask__bc", &temp_UpsilonDataToTree_f[4]);
   reader->AddSpectator("nROE_ParticlesInList__boK_S0__clmyKaonshort__bc", &temp_UpsilonDataToTree_f[5]);
   reader->AddSpectator("nROE_ParticlesInList__bopi0__clmyneutralPion__bc", &temp_UpsilonDataToTree_f[6]);
   reader->AddSpectator("missingMomentumOfEvent_theta", &temp_UpsilonDataToTree_f[7]);
   reader->AddSpectator("missingMomentumOfEvent", &temp_UpsilonDataToTree_f[8]);
   reader->AddSpectator("missingEnergyOfEventCMS", &temp_UpsilonDataToTree_f[9]);
   reader->AddSpectator("nRemainingTracksInEvent", &temp_UpsilonDataToTree_f[10]);
   reader->AddSpectator("roeNeextra__bocleanMask__bc", &temp_UpsilonDataToTree_f[11]);
   reader->AddSpectator("useCMSFrame__boroeNeextra__bocleanMask__bc__bc", &temp_UpsilonDataToTree_f[12]);
   reader->AddSpectator("nROE_ParticlesInList__bogamma__clmygamma__bc", &temp_UpsilonDataToTree_f[13]);
   reader->AddSpectator("Bsig_E", &temp_BsigDataToTree_f[0]);
   reader->AddSpectator("Bsig_useCMSFrame_E", &temp_BsigDataToTree_f[1]);
   reader->AddSpectator("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc", &temp_BsigDataToTree_f[2]);
   reader->AddSpectator("Bsig_p", &temp_BsigDataToTree_f[3]);
   reader->AddSpectator("Bsig_useCMSFrame_p", &temp_BsigDataToTree_f[4]);
   reader->AddSpectator("useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc", &temp_BsigDataToTree_f[5]);
   reader->AddSpectator("Bsig_M", &temp_BsigDataToTree_f[6]);
   reader->AddSpectator("Bsig_daughter_0_extraInfo_Dcvetomass", &temp_BsigDataToTree_f[7]);
   reader->AddSpectator("Bsig_daughter_0_extraInfo_DcvetodmID", &temp_BsigDataToTree_f[8]);
   reader->AddSpectator("Bsig_daughter_0_extraInfo_Dcvetoabsdm", &temp_BsigDataToTree_f[9]);
   reader->AddSpectator("Bsig_daughter_0_extraInfo_Dnvetomass", &temp_BsigDataToTree_f[10]);
   reader->AddSpectator("Bsig_daughter_0_extraInfo_DnvetodmID", &temp_BsigDataToTree_f[11]);
   reader->AddSpectator("Bsig_daughter_0_extraInfo_Dnvetoabsdm", &temp_BsigDataToTree_f[12]);
   reader->AddSpectator("Bsig_daughter_0_extraInfo_KS0_3D_distance", &temp_BsigDataToTree_f[13]);
   reader->AddSpectator("Bsig_daughter_0_extraInfo_Kaon_PID_correction", &temp_BsigDataToTree_f[14]);
   reader->AddSpectator("Bsig_daughter_0_extraInfo_Kaon_PID_rel_uncer_up", &temp_BsigDataToTree_f[15]);
   reader->AddSpectator("Bsig_daughter_0_extraInfo_Kaon_PID_rel_uncer_dn", &temp_BsigDataToTree_f[16]);
   reader->AddSpectator("Bsig_daughter_0_extraInfo_Pion_PID_correction", &temp_BsigDataToTree_f[17]);
   reader->AddSpectator("Bsig_daughter_0_extraInfo_Pion_PID_rel_uncer_up", &temp_BsigDataToTree_f[18]);
   reader->AddSpectator("Bsig_daughter_0_extraInfo_Pion_PID_rel_uncer_dn", &temp_BsigDataToTree_f[19]);
   reader->AddSpectator("Bsig_daughter_0_extraInfo_nKslow1", &temp_BsigDataToTree_f[20]);
   reader->AddSpectator("Bsig_daughter_0_extraInfo_nKslow2", &temp_BsigDataToTree_f[21]);
   reader->AddSpectator("Bsig_daughter_0_extraInfo_nKslow3", &temp_BsigDataToTree_f[22]);
   reader->AddSpectator("Bsig_daughter_0_extraInfo_nPislow1", &temp_BsigDataToTree_f[23]);
   reader->AddSpectator("Bsig_daughter_0_extraInfo_nPislow2", &temp_BsigDataToTree_f[24]);
   reader->AddSpectator("Bsig_daughter_0_extraInfo_nPislow3", &temp_BsigDataToTree_f[25]);
   reader->AddSpectator("Bsig_daughter_0_extraInfo_nKexcep", &temp_BsigDataToTree_f[26]);
   reader->AddSpectator("Bsig_daughter_0_extraInfo_nPiexcep", &temp_BsigDataToTree_f[27]);
   reader->AddSpectator("Btag_extraInfo_decayModeID", &temp_BtagDataToTree_f[0]);
   reader->AddSpectator("Btag_Mbc", &temp_BtagDataToTree_f[1]);
   reader->AddSpectator("Btag_deltaE", &temp_BtagDataToTree_f[2]);
   reader->AddSpectator("Btag_E", &temp_BtagDataToTree_f[3]);
   reader->AddSpectator("missingMass2OfEvent", &temp_DataToTree_f[32]);
   reader->AddSpectator("visibleEnergyOfEventCMS", &temp_DataToTree_f[33]);
   reader->AddSpectator("Btag_useCMSFrame_E", &temp_BtagDataToTree_f[4]);
   reader->AddSpectator("flag", &temp_flag_f);

   TString dir_continuum    = "";
   TString dir_BB = "";
   TString prefix = "";
   if (OPTION == std::string("large")) {
       dir_continuum = "dataset_Continuum_Mxs_large/weights/";
       dir_BB = "dataset_BB_Mxs_large/weights/";
       prefix = "TMVAClassification";
   }
   else if (OPTION == std::string("small")) {
       dir_continuum = "dataset_Continuum_Mxs_small/weights/";
       dir_BB = "dataset_BB_Mxs_small/weights/";
       prefix = "TMVAClassification";
   }

   TString methodName_continuum = TString("Continuum");
   TString weightfile_continuum = dir_continuum + prefix + TString("_") + TString("MLP") + TString(".weights.xml");
   reader->BookMVA(methodName_continuum, weightfile_continuum);

   TString methodName_BB = TString("BB");
   TString weightfile_BB = dir_BB + prefix + TString("_") + TString("MLP") + TString(".weights.xml");
   reader->BookMVA(methodName_BB, weightfile_BB);

   int nbin = 100;
   TH1F *histMLP_Continuum(0);
   TH1F* histMLP_BB(0);

   histMLP_Continuum = new TH1F("MVA_MLP_Continuum", "MVA_MLP_Continuum", nbin, -0.8, 0.8);
   histMLP_BB = new TH1F("MVA_MLP_BB", "MVA_MLP_BB", nbin, -0.8, 0.8);

   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //
   TFile *input(0);
   TString fname = std::string(filename);
   input = TFile::Open( fname );
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;

   // Event loop

   // Prepare the event tree
   // - Here the variable names have to corresponds to your tree
   // - You can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   std::cout << "--- Select signal sample" << std::endl;
   TTree* theTree = (TTree*)input->Get("data");
   // get event_info
   theTree->SetBranchAddress("__experiment__", &temp_EventDataToTree[0]);
   theTree->SetBranchAddress("__run__", &temp_EventDataToTree[1]);
   theTree->SetBranchAddress("__event__", &temp_EventDataToTree[2]);
   theTree->SetBranchAddress("__candidate__", &temp_EventDataToTree[3]);
   theTree->SetBranchAddress("__ncandidates__", &temp_EventDataToTree[4]);

   // get decaymodeID
   theTree->SetBranchAddress("extraInfo__bodecayModeID__bc", &temp_Upsilon_decayIDToTree);
   theTree->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &temp_Bsig_decayIDToTree);

   // get Upsilon_info
   theTree->SetBranchAddress("nROE_ECLClusters__bocleanMask__bc", &temp_UpsilonDataToTree[0]);
   theTree->SetBranchAddress("nROE_KLMClusters", &temp_UpsilonDataToTree[1]);
   theTree->SetBranchAddress("nROE_Tracks__bocleanMask__bc", &temp_UpsilonDataToTree[2]);
   theTree->SetBranchAddress("roeEextra__bocleanMask__bc", &temp_UpsilonDataToTree[3]);
   theTree->SetBranchAddress("nROE_NeutralECLClusters__bocleanMask__bc", &temp_UpsilonDataToTree[4]);
   theTree->SetBranchAddress("nROE_ParticlesInList__boK_S0__clmyKaonshort__bc", &temp_UpsilonDataToTree[5]);
   theTree->SetBranchAddress("nROE_ParticlesInList__bopi0__clmyneutralPion__bc", &temp_UpsilonDataToTree[6]);
   theTree->SetBranchAddress("missingMomentumOfEvent_theta", &temp_UpsilonDataToTree[7]);
   theTree->SetBranchAddress("missingMomentumOfEvent", &temp_UpsilonDataToTree[8]);
   theTree->SetBranchAddress("missingEnergyOfEventCMS", &temp_UpsilonDataToTree[9]);
   theTree->SetBranchAddress("nRemainingTracksInEvent", &temp_UpsilonDataToTree[10]);
   theTree->SetBranchAddress("roeNeextra__bocleanMask__bc", &temp_UpsilonDataToTree[11]);
   theTree->SetBranchAddress("useCMSFrame__boroeNeextra__bocleanMask__bc__bc", &temp_UpsilonDataToTree[12]);
   theTree->SetBranchAddress("nROE_ParticlesInList__bogamma__clmygamma__bc", &temp_UpsilonDataToTree[13]);

   // get Bsig_info
   theTree->SetBranchAddress("Bsig_E", &temp_BsigDataToTree[0]);
   theTree->SetBranchAddress("Bsig_useCMSFrame_E", &temp_BsigDataToTree[1]);
   theTree->SetBranchAddress("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc", &temp_BsigDataToTree[2]);
   theTree->SetBranchAddress("Bsig_p", &temp_BsigDataToTree[3]);
   theTree->SetBranchAddress("Bsig_useCMSFrame_p", &temp_BsigDataToTree[4]);
   theTree->SetBranchAddress("useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc", &temp_BsigDataToTree[5]);
   theTree->SetBranchAddress("Bsig_M", &temp_BsigDataToTree[6]);
   theTree->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcvetomass", &temp_BsigDataToTree[7]);
   theTree->SetBranchAddress("Bsig_daughter_0_extraInfo_DcvetodmID", &temp_BsigDataToTree[8]);
   theTree->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcvetoabsdm", &temp_BsigDataToTree[9]);
   theTree->SetBranchAddress("Bsig_daughter_0_extraInfo_Dnvetomass", &temp_BsigDataToTree[10]);
   theTree->SetBranchAddress("Bsig_daughter_0_extraInfo_DnvetodmID", &temp_BsigDataToTree[11]);
   theTree->SetBranchAddress("Bsig_daughter_0_extraInfo_Dnvetoabsdm", &temp_BsigDataToTree[12]);
   theTree->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_3D_distance", &temp_BsigDataToTree[13]);
   theTree->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_correction", &temp_BsigDataToTree[14]);
   theTree->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_rel_uncer_up", &temp_BsigDataToTree[15]);
   theTree->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_rel_uncer_dn", &temp_BsigDataToTree[16]);
   theTree->SetBranchAddress("Bsig_daughter_0_extraInfo_Pion_PID_correction", &temp_BsigDataToTree[17]);
   theTree->SetBranchAddress("Bsig_daughter_0_extraInfo_Pion_PID_rel_uncer_up", &temp_BsigDataToTree[18]);
   theTree->SetBranchAddress("Bsig_daughter_0_extraInfo_Pion_PID_rel_uncer_dn", &temp_BsigDataToTree[19]);
   theTree->SetBranchAddress("Bsig_daughter_0_extraInfo_nKslow1", &temp_BsigDataToTree[20]);
   theTree->SetBranchAddress("Bsig_daughter_0_extraInfo_nKslow2", &temp_BsigDataToTree[21]);
   theTree->SetBranchAddress("Bsig_daughter_0_extraInfo_nKslow3", &temp_BsigDataToTree[22]);
   theTree->SetBranchAddress("Bsig_daughter_0_extraInfo_nPislow1", &temp_BsigDataToTree[23]);
   theTree->SetBranchAddress("Bsig_daughter_0_extraInfo_nPislow2", &temp_BsigDataToTree[24]);
   theTree->SetBranchAddress("Bsig_daughter_0_extraInfo_nPislow3", &temp_BsigDataToTree[25]);
   theTree->SetBranchAddress("Bsig_daughter_0_extraInfo_nKexcep", &temp_BsigDataToTree[26]);
   theTree->SetBranchAddress("Bsig_daughter_0_extraInfo_nPiexcep", &temp_BsigDataToTree[27]);

   // get Btag_info
   theTree->SetBranchAddress("Btag_extraInfo_decayModeID", &temp_BtagDataToTree[0]);
   theTree->SetBranchAddress("Btag_Mbc", &temp_BtagDataToTree[1]);
   theTree->SetBranchAddress("Btag_deltaE", &temp_BtagDataToTree[2]);
   theTree->SetBranchAddress("Btag_E", &temp_BtagDataToTree[3]);
   theTree->SetBranchAddress("Btag_useCMSFrame_E", &temp_BtagDataToTree[4]);
   theTree->SetBranchAddress("Btag_extraInfo_SignalProbability", &temp_BtagDataToTree[5]);
   theTree->SetBranchAddress("Btag_chiProb", &temp_BtagDataToTree[6]);

   // other information I need
   theTree->SetBranchAddress("Btag_R2", &temp_DataToTree[0]);
   theTree->SetBranchAddress("Btag_thrustBm", &temp_DataToTree[1]);
   theTree->SetBranchAddress("Btag_thrustOm", &temp_DataToTree[2]);
   theTree->SetBranchAddress("Btag_cosTBTO", &temp_DataToTree[3]);
   theTree->SetBranchAddress("Btag_cosTBz", &temp_DataToTree[4]);
   theTree->SetBranchAddress("Btag_KSFWVariables_et", &temp_DataToTree[5]);
   theTree->SetBranchAddress("Btag_KSFWVariables_mm2", &temp_DataToTree[6]);
   theTree->SetBranchAddress("Btag_KSFWVariables_hso00", &temp_DataToTree[7]);
   theTree->SetBranchAddress("Btag_KSFWVariables_hso01", &temp_DataToTree[8]);
   theTree->SetBranchAddress("Btag_KSFWVariables_hso02", &temp_DataToTree[9]);
   theTree->SetBranchAddress("Btag_KSFWVariables_hso03", &temp_DataToTree[10]);
   theTree->SetBranchAddress("Btag_KSFWVariables_hso04", &temp_DataToTree[11]);
   theTree->SetBranchAddress("Btag_KSFWVariables_hso10", &temp_DataToTree[12]);
   theTree->SetBranchAddress("Btag_KSFWVariables_hso12", &temp_DataToTree[13]);
   theTree->SetBranchAddress("Btag_KSFWVariables_hso14", &temp_DataToTree[14]);
   theTree->SetBranchAddress("Btag_KSFWVariables_hso20", &temp_DataToTree[15]);
   theTree->SetBranchAddress("Btag_KSFWVariables_hso22", &temp_DataToTree[16]);
   theTree->SetBranchAddress("Btag_KSFWVariables_hso24", &temp_DataToTree[17]);
   theTree->SetBranchAddress("Btag_KSFWVariables_hoo0", &temp_DataToTree[18]);
   theTree->SetBranchAddress("Btag_KSFWVariables_hoo1", &temp_DataToTree[19]);
   theTree->SetBranchAddress("Btag_KSFWVariables_hoo2", &temp_DataToTree[20]);
   theTree->SetBranchAddress("Btag_KSFWVariables_hoo3", &temp_DataToTree[21]);
   theTree->SetBranchAddress("Btag_KSFWVariables_hoo4", &temp_DataToTree[22]);
   theTree->SetBranchAddress("Btag_CleoConeCS_1", &temp_DataToTree[23]);
   theTree->SetBranchAddress("Btag_CleoConeCS_2", &temp_DataToTree[24]);
   theTree->SetBranchAddress("Btag_CleoConeCS_3", &temp_DataToTree[25]);
   theTree->SetBranchAddress("Btag_CleoConeCS_4", &temp_DataToTree[26]);
   theTree->SetBranchAddress("Btag_CleoConeCS_5", &temp_DataToTree[27]);
   theTree->SetBranchAddress("Btag_CleoConeCS_6", &temp_DataToTree[28]);
   theTree->SetBranchAddress("Btag_CleoConeCS_7", &temp_DataToTree[29]);
   theTree->SetBranchAddress("Btag_CleoConeCS_8", &temp_DataToTree[30]);
   theTree->SetBranchAddress("Btag_CleoConeCS_9", &temp_DataToTree[31]);
   theTree->SetBranchAddress("missingMass2OfEvent", &temp_DataToTree[32]);
   theTree->SetBranchAddress("visibleEnergyOfEventCMS", &temp_DataToTree[33]);
   theTree->SetBranchAddress("Btag_useCMSFrame_theta", &temp_DataToTree[34]);
   theTree->SetBranchAddress("extraInfo__boDecayHash__bc", &temp_DataToTree[35]);
   theTree->SetBranchAddress("extraInfo__boDecayHashExtended__bc", &temp_DataToTree[36]);

   bool DoesItHaveXsBranch = false;
   for (int i = 0; i < theTree->GetListOfBranches()->LastIndex(); i++) {
       if (theTree->GetListOfBranches()->At(i)->GetName() == std::string("nParticlesInList__boB__pl__clKcharge_total__bc")) DoesItHaveXsBranch = true;
   }
   if (DoesItHaveXsBranch) {
       // decay mode (MC level)
       theTree->SetBranchAddress("nParticlesInList__boB__pl__clKcharge_total__bc", &temp_DecayDataToTree[0]);
       theTree->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &temp_DecayDataToTree[1]);
       theTree->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &temp_DecayDataToTree[2]);
       theTree->SetBranchAddress("nParticlesInList__boXsu__clMCcomb__bc", &temp_DecayDataToTree[3]);
       theTree->SetBranchAddress("nParticlesInList__boXsu__clMCch1__bc", &temp_DecayDataToTree[4]);
       theTree->SetBranchAddress("nParticlesInList__boXsu__clMCch2__bc", &temp_DecayDataToTree[5]);
       theTree->SetBranchAddress("nParticlesInList__boXsu__clMCch3__bc", &temp_DecayDataToTree[6]);
       theTree->SetBranchAddress("nParticlesInList__boXsu__clMCch4__bc", &temp_DecayDataToTree[7]);
       theTree->SetBranchAddress("nParticlesInList__boXsu__clMCch5__bc", &temp_DecayDataToTree[8]);
       theTree->SetBranchAddress("nParticlesInList__boXsu__clMCch6__bc", &temp_DecayDataToTree[9]);
       theTree->SetBranchAddress("nParticlesInList__boXsu__clMCch7__bc", &temp_DecayDataToTree[10]);
       theTree->SetBranchAddress("nParticlesInList__boXsu__clMCch8__bc", &temp_DecayDataToTree[11]);
       theTree->SetBranchAddress("nParticlesInList__boXsu__clMCch9__bc", &temp_DecayDataToTree[12]);
       theTree->SetBranchAddress("nParticlesInList__boXsu__clMCch10__bc", &temp_DecayDataToTree[13]);
       theTree->SetBranchAddress("nParticlesInList__boXsu__clMCch11__bc", &temp_DecayDataToTree[14]);
       theTree->SetBranchAddress("nParticlesInList__boXsu__clMCch12__bc", &temp_DecayDataToTree[15]);
       theTree->SetBranchAddress("nParticlesInList__boXsu__clMCch13__bc", &temp_DecayDataToTree[16]);
       theTree->SetBranchAddress("nParticlesInList__boXsu__clMCch14__bc", &temp_DecayDataToTree[17]);
       theTree->SetBranchAddress("nParticlesInList__boXsu__clMCch15__bc", &temp_DecayDataToTree[18]);
       theTree->SetBranchAddress("nParticlesInList__boB0__clKneutral_total__bc", &temp_DecayDataToTree[19]);
       theTree->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &temp_DecayDataToTree[20]);
       theTree->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &temp_DecayDataToTree[21]);
       theTree->SetBranchAddress("nParticlesInList__boXsd__clMCcomb__bc", &temp_DecayDataToTree[22]);
       theTree->SetBranchAddress("nParticlesInList__boXsd__clMCch16__bc", &temp_DecayDataToTree[23]);
       theTree->SetBranchAddress("nParticlesInList__boXsd__clMCch17__bc", &temp_DecayDataToTree[24]);
       theTree->SetBranchAddress("nParticlesInList__boXsd__clMCch18__bc", &temp_DecayDataToTree[25]);
       theTree->SetBranchAddress("nParticlesInList__boXsd__clMCch19__bc", &temp_DecayDataToTree[26]);
       theTree->SetBranchAddress("nParticlesInList__boXsd__clMCch20__bc", &temp_DecayDataToTree[27]);
       theTree->SetBranchAddress("nParticlesInList__boXsd__clMCch21__bc", &temp_DecayDataToTree[28]);
       theTree->SetBranchAddress("nParticlesInList__boXsd__clMCch22__bc", &temp_DecayDataToTree[29]);
       theTree->SetBranchAddress("nParticlesInList__boXsd__clMCch23__bc", &temp_DecayDataToTree[30]);
       theTree->SetBranchAddress("nParticlesInList__boXsd__clMCch24__bc", &temp_DecayDataToTree[31]);
       theTree->SetBranchAddress("nParticlesInList__boXsd__clMCch25__bc", &temp_DecayDataToTree[32]);
       theTree->SetBranchAddress("nParticlesInList__boXsd__clMCch26__bc", &temp_DecayDataToTree[33]);
       theTree->SetBranchAddress("nParticlesInList__boXsd__clMCch27__bc", &temp_DecayDataToTree[34]);
       theTree->SetBranchAddress("nParticlesInList__boXsd__clMCch28__bc", &temp_DecayDataToTree[35]);
       theTree->SetBranchAddress("nParticlesInList__boXsd__clMCch29__bc", &temp_DecayDataToTree[36]);
       theTree->SetBranchAddress("nParticlesInList__boXsd__clMCch30__bc", &temp_DecayDataToTree[37]);
   }

   theTree->SetBranchAddress("flag", &temp_flag);

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;

   TFile* temp_file = new TFile( ("TMVAoutput_"+ OnlyFileName).c_str(), "recreate");
   temp_file->cd();
   TTree* temp_tree = new TTree("data", "");

   /*================================================================*/
   // get event_info
   temp_tree->Branch("__experiment__", &temp_EventDataToTree[0]);
   temp_tree->Branch("__run__", &temp_EventDataToTree[1]);
   temp_tree->Branch("__event__", &temp_EventDataToTree[2]);
   temp_tree->Branch("__candidate__", &temp_EventDataToTree[3]);
   temp_tree->Branch("__ncandidates__", &temp_EventDataToTree[4]);

   // get decaymodeID
   temp_tree->Branch("extraInfo__bodecayModeID__bc", &temp_Upsilon_decayIDToTree);
   temp_tree->Branch("Bsig_daughter_0_extraInfo_decayModeID", &temp_Bsig_decayIDToTree);

   // get Upsilon_info
   temp_tree->Branch("nROE_ECLClusters__bocleanMask__bc", &temp_UpsilonDataToTree[0]);
   temp_tree->Branch("nROE_KLMClusters", &temp_UpsilonDataToTree[1]);
   temp_tree->Branch("nROE_Tracks__bocleanMask__bc", &temp_UpsilonDataToTree[2]);
   temp_tree->Branch("roeEextra__bocleanMask__bc", &temp_UpsilonDataToTree[3]);
   temp_tree->Branch("nROE_NeutralECLClusters__bocleanMask__bc", &temp_UpsilonDataToTree[4]);
   temp_tree->Branch("nROE_ParticlesInList__boK_S0__clmyKaonshort__bc", &temp_UpsilonDataToTree[5]);
   temp_tree->Branch("nROE_ParticlesInList__bopi0__clmyneutralPion__bc", &temp_UpsilonDataToTree[6]);
   temp_tree->Branch("missingMomentumOfEvent_theta", &temp_UpsilonDataToTree[7]);
   temp_tree->Branch("missingMomentumOfEvent", &temp_UpsilonDataToTree[8]);
   temp_tree->Branch("missingEnergyOfEventCMS", &temp_UpsilonDataToTree[9]);
   temp_tree->Branch("nRemainingTracksInEvent", &temp_UpsilonDataToTree[10]);
   temp_tree->Branch("roeNeextra__bocleanMask__bc", &temp_UpsilonDataToTree[11]);
   temp_tree->Branch("useCMSFrame__boroeNeextra__bocleanMask__bc__bc", &temp_UpsilonDataToTree[12]);
   temp_tree->Branch("nROE_ParticlesInList__bogamma__clmygamma__bc", &temp_UpsilonDataToTree[13]);

   // get Bsig_info
   temp_tree->Branch("Bsig_E", &temp_BsigDataToTree[0]);
   temp_tree->Branch("Bsig_useCMSFrame_E", &temp_BsigDataToTree[1]);
   temp_tree->Branch("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc", &temp_BsigDataToTree[2]);
   temp_tree->Branch("Bsig_p", &temp_BsigDataToTree[3]);
   temp_tree->Branch("Bsig_useCMSFrame_p", &temp_BsigDataToTree[4]);
   temp_tree->Branch("useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc", &temp_BsigDataToTree[5]);
   temp_tree->Branch("Bsig_M", &temp_BsigDataToTree[6]);
   temp_tree->Branch("Bsig_daughter_0_extraInfo_Dcvetomass", &temp_BsigDataToTree[7]);
   temp_tree->Branch("Bsig_daughter_0_extraInfo_DcvetodmID", &temp_BsigDataToTree[8]);
   temp_tree->Branch("Bsig_daughter_0_extraInfo_Dcvetoabsdm", &temp_BsigDataToTree[9]);
   temp_tree->Branch("Bsig_daughter_0_extraInfo_Dnvetomass", &temp_BsigDataToTree[10]);
   temp_tree->Branch("Bsig_daughter_0_extraInfo_DnvetodmID", &temp_BsigDataToTree[11]);
   temp_tree->Branch("Bsig_daughter_0_extraInfo_Dnvetoabsdm", &temp_BsigDataToTree[12]);
   temp_tree->Branch("Bsig_daughter_0_extraInfo_KS0_3D_distance", &temp_BsigDataToTree[13]);
   temp_tree->Branch("Bsig_daughter_0_extraInfo_Kaon_PID_correction", &temp_BsigDataToTree[14]);
   temp_tree->Branch("Bsig_daughter_0_extraInfo_Kaon_PID_rel_uncer_up", &temp_BsigDataToTree[15]);
   temp_tree->Branch("Bsig_daughter_0_extraInfo_Kaon_PID_rel_uncer_dn", &temp_BsigDataToTree[16]);
   temp_tree->Branch("Bsig_daughter_0_extraInfo_Pion_PID_correction", &temp_BsigDataToTree[17]);
   temp_tree->Branch("Bsig_daughter_0_extraInfo_Pion_PID_rel_uncer_up", &temp_BsigDataToTree[18]);
   temp_tree->Branch("Bsig_daughter_0_extraInfo_Pion_PID_rel_uncer_dn", &temp_BsigDataToTree[19]);
   temp_tree->Branch("Bsig_daughter_0_extraInfo_nKslow1", &temp_BsigDataToTree[20]);
   temp_tree->Branch("Bsig_daughter_0_extraInfo_nKslow2", &temp_BsigDataToTree[21]);
   temp_tree->Branch("Bsig_daughter_0_extraInfo_nKslow3", &temp_BsigDataToTree[22]);
   temp_tree->Branch("Bsig_daughter_0_extraInfo_nPislow1", &temp_BsigDataToTree[23]);
   temp_tree->Branch("Bsig_daughter_0_extraInfo_nPislow2", &temp_BsigDataToTree[24]);
   temp_tree->Branch("Bsig_daughter_0_extraInfo_nPislow3", &temp_BsigDataToTree[25]);
   temp_tree->Branch("Bsig_daughter_0_extraInfo_nKexcep", &temp_BsigDataToTree[26]);
   temp_tree->Branch("Bsig_daughter_0_extraInfo_nPiexcep", &temp_BsigDataToTree[27]);

   // get Btag_info
   temp_tree->Branch("Btag_extraInfo_decayModeID", &temp_BtagDataToTree[0]);
   temp_tree->Branch("Btag_Mbc", &temp_BtagDataToTree[1]);
   temp_tree->Branch("Btag_deltaE", &temp_BtagDataToTree[2]);
   temp_tree->Branch("Btag_E", &temp_BtagDataToTree[3]);
   temp_tree->Branch("Btag_useCMSFrame_E", &temp_BtagDataToTree[4]);
   temp_tree->Branch("Btag_extraInfo_SignalProbability", &temp_BtagDataToTree[5]);
   temp_tree->Branch("Btag_chiProb", &temp_BtagDataToTree[6]);

   // other information I need
   temp_tree->Branch("Btag_R2", &temp_DataToTree[0]);
   temp_tree->Branch("Btag_thrustBm", &temp_DataToTree[1]);
   temp_tree->Branch("Btag_thrustOm", &temp_DataToTree[2]);
   temp_tree->Branch("Btag_cosTBTO", &temp_DataToTree[3]);
   temp_tree->Branch("Btag_cosTBz", &temp_DataToTree[4]);
   temp_tree->Branch("Btag_KSFWVariables_et", &temp_DataToTree[5]);
   temp_tree->Branch("Btag_KSFWVariables_mm2", &temp_DataToTree[6]);
   temp_tree->Branch("Btag_KSFWVariables_hso00", &temp_DataToTree[7]);
   temp_tree->Branch("Btag_KSFWVariables_hso01", &temp_DataToTree[8]);
   temp_tree->Branch("Btag_KSFWVariables_hso02", &temp_DataToTree[9]);
   temp_tree->Branch("Btag_KSFWVariables_hso03", &temp_DataToTree[10]);
   temp_tree->Branch("Btag_KSFWVariables_hso04", &temp_DataToTree[11]);
   temp_tree->Branch("Btag_KSFWVariables_hso10", &temp_DataToTree[12]);
   temp_tree->Branch("Btag_KSFWVariables_hso12", &temp_DataToTree[13]);
   temp_tree->Branch("Btag_KSFWVariables_hso14", &temp_DataToTree[14]);
   temp_tree->Branch("Btag_KSFWVariables_hso20", &temp_DataToTree[15]);
   temp_tree->Branch("Btag_KSFWVariables_hso22", &temp_DataToTree[16]);
   temp_tree->Branch("Btag_KSFWVariables_hso24", &temp_DataToTree[17]);
   temp_tree->Branch("Btag_KSFWVariables_hoo0", &temp_DataToTree[18]);
   temp_tree->Branch("Btag_KSFWVariables_hoo1", &temp_DataToTree[19]);
   temp_tree->Branch("Btag_KSFWVariables_hoo2", &temp_DataToTree[20]);
   temp_tree->Branch("Btag_KSFWVariables_hoo3", &temp_DataToTree[21]);
   temp_tree->Branch("Btag_KSFWVariables_hoo4", &temp_DataToTree[22]);
   temp_tree->Branch("Btag_CleoConeCS_1", &temp_DataToTree[23]);
   temp_tree->Branch("Btag_CleoConeCS_2", &temp_DataToTree[24]);
   temp_tree->Branch("Btag_CleoConeCS_3", &temp_DataToTree[25]);
   temp_tree->Branch("Btag_CleoConeCS_4", &temp_DataToTree[26]);
   temp_tree->Branch("Btag_CleoConeCS_5", &temp_DataToTree[27]);
   temp_tree->Branch("Btag_CleoConeCS_6", &temp_DataToTree[28]);
   temp_tree->Branch("Btag_CleoConeCS_7", &temp_DataToTree[29]);
   temp_tree->Branch("Btag_CleoConeCS_8", &temp_DataToTree[30]);
   temp_tree->Branch("Btag_CleoConeCS_9", &temp_DataToTree[31]);
   temp_tree->Branch("missingMass2OfEvent", &temp_DataToTree[32]);
   temp_tree->Branch("visibleEnergyOfEventCMS", &temp_DataToTree[33]);
   temp_tree->Branch("Btag_useCMSFrame_theta", &temp_DataToTree[34]);
   temp_tree->Branch("extraInfo__boDecayHash__bc", &temp_DataToTree[35]);
   temp_tree->Branch("extraInfo__boDecayHashExtended__bc", &temp_DataToTree[36]);

   if (DoesItHaveXsBranch) {
       // decay mode (MC level)
       temp_tree->Branch("nParticlesInList__boB__pl__clKcharge_total__bc", &temp_DecayDataToTree[0]);
       temp_tree->Branch("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &temp_DecayDataToTree[1]);
       temp_tree->Branch("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &temp_DecayDataToTree[2]);
       temp_tree->Branch("nParticlesInList__boXsu__clMCcomb__bc", &temp_DecayDataToTree[3]);
       temp_tree->Branch("nParticlesInList__boXsu__clMCch1__bc", &temp_DecayDataToTree[4]);
       temp_tree->Branch("nParticlesInList__boXsu__clMCch2__bc", &temp_DecayDataToTree[5]);
       temp_tree->Branch("nParticlesInList__boXsu__clMCch3__bc", &temp_DecayDataToTree[6]);
       temp_tree->Branch("nParticlesInList__boXsu__clMCch4__bc", &temp_DecayDataToTree[7]);
       temp_tree->Branch("nParticlesInList__boXsu__clMCch5__bc", &temp_DecayDataToTree[8]);
       temp_tree->Branch("nParticlesInList__boXsu__clMCch6__bc", &temp_DecayDataToTree[9]);
       temp_tree->Branch("nParticlesInList__boXsu__clMCch7__bc", &temp_DecayDataToTree[10]);
       temp_tree->Branch("nParticlesInList__boXsu__clMCch8__bc", &temp_DecayDataToTree[11]);
       temp_tree->Branch("nParticlesInList__boXsu__clMCch9__bc", &temp_DecayDataToTree[12]);
       temp_tree->Branch("nParticlesInList__boXsu__clMCch10__bc", &temp_DecayDataToTree[13]);
       temp_tree->Branch("nParticlesInList__boXsu__clMCch11__bc", &temp_DecayDataToTree[14]);
       temp_tree->Branch("nParticlesInList__boXsu__clMCch12__bc", &temp_DecayDataToTree[15]);
       temp_tree->Branch("nParticlesInList__boXsu__clMCch13__bc", &temp_DecayDataToTree[16]);
       temp_tree->Branch("nParticlesInList__boXsu__clMCch14__bc", &temp_DecayDataToTree[17]);
       temp_tree->Branch("nParticlesInList__boXsu__clMCch15__bc", &temp_DecayDataToTree[18]);
       temp_tree->Branch("nParticlesInList__boB0__clKneutral_total__bc", &temp_DecayDataToTree[19]);
       temp_tree->Branch("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &temp_DecayDataToTree[20]);
       temp_tree->Branch("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &temp_DecayDataToTree[21]);
       temp_tree->Branch("nParticlesInList__boXsd__clMCcomb__bc", &temp_DecayDataToTree[22]);
       temp_tree->Branch("nParticlesInList__boXsd__clMCch16__bc", &temp_DecayDataToTree[23]);
       temp_tree->Branch("nParticlesInList__boXsd__clMCch17__bc", &temp_DecayDataToTree[24]);
       temp_tree->Branch("nParticlesInList__boXsd__clMCch18__bc", &temp_DecayDataToTree[25]);
       temp_tree->Branch("nParticlesInList__boXsd__clMCch19__bc", &temp_DecayDataToTree[26]);
       temp_tree->Branch("nParticlesInList__boXsd__clMCch20__bc", &temp_DecayDataToTree[27]);
       temp_tree->Branch("nParticlesInList__boXsd__clMCch21__bc", &temp_DecayDataToTree[28]);
       temp_tree->Branch("nParticlesInList__boXsd__clMCch22__bc", &temp_DecayDataToTree[29]);
       temp_tree->Branch("nParticlesInList__boXsd__clMCch23__bc", &temp_DecayDataToTree[30]);
       temp_tree->Branch("nParticlesInList__boXsd__clMCch24__bc", &temp_DecayDataToTree[31]);
       temp_tree->Branch("nParticlesInList__boXsd__clMCch25__bc", &temp_DecayDataToTree[32]);
       temp_tree->Branch("nParticlesInList__boXsd__clMCch26__bc", &temp_DecayDataToTree[33]);
       temp_tree->Branch("nParticlesInList__boXsd__clMCch27__bc", &temp_DecayDataToTree[34]);
       temp_tree->Branch("nParticlesInList__boXsd__clMCch28__bc", &temp_DecayDataToTree[35]);
       temp_tree->Branch("nParticlesInList__boXsd__clMCch29__bc", &temp_DecayDataToTree[36]);
       temp_tree->Branch("nParticlesInList__boXsd__clMCch30__bc", &temp_DecayDataToTree[37]);
   }

   // flag
   temp_tree->Branch("flag", &temp_flag);

   // TMVA output
   temp_tree->Branch("TMVA_BB", &Output_BB);
   temp_tree->Branch("TMVA_Continuum", &Output_Continuum);

   for (Long64_t ievt = 0; ievt < theTree->GetEntries(); ievt++) {

       if (ievt % 1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

       theTree->GetEntry(ievt);

       for (int i = 0; i < N_event_info / 3; i++) temp_EventDataToTree_f[i] = (float)temp_EventDataToTree[i];
       for (int i = 0; i < N_Upsilon_info; i++) temp_UpsilonDataToTree_f[i] = (float)temp_UpsilonDataToTree[i];
       for (int i = 0; i < N_Bsig_info; i++) temp_BsigDataToTree_f[i] = (float)temp_BsigDataToTree[i];
       for (int i = 0; i < N_Btag_info; i++) temp_BtagDataToTree_f[i] = (float)temp_BtagDataToTree[i];
       for (int i = 0; i < N_Needed_info ; i++) temp_DataToTree_f[i] = (float)temp_DataToTree[i];
       temp_Upsilon_decayIDToTree_f = (float)temp_Upsilon_decayIDToTree;
       temp_Bsig_decayIDToTree_f = (float)temp_Bsig_decayIDToTree;
       temp_flag_f = (float)temp_flag;

       Output_BB = reader->EvaluateMVA("BB");
       Output_Continuum = reader->EvaluateMVA("Continuum");

       temp_tree->Fill();
   }

   temp_file->cd();
   temp_tree->Write();
   temp_file->Close();

   delete reader;

   std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
   
}