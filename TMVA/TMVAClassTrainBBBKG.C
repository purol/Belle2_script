/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This macro provides a simple example for the training and testing of the TMVA
/// multiclass classification
/// - Project   : TMVA - a Root-integrated toolkit for multivariate data analysis
/// - Package   : TMVA
/// - Root Macro: TMVAMulticlass
///
/// \macro_output
/// \macro_code
/// \author Andreas Hoecker

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"


#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVAMultiClassGui.h"


using namespace TMVA;

void TMVAClassTrainBBBKG( TString myMethodList = "" )
{

   // This loads the library
   TMVA::Tools::Instance();

   // to get access to the GUI and all tmva macros
   //
   //     TString tmva_dir(TString(gRootDir) + "/tmva");
   //     if(gSystem->Getenv("TMVASYS"))
   //        tmva_dir = TString(gSystem->Getenv("TMVASYS"));
   //     gROOT->SetMacroPath(tmva_dir + "/test/:" + gROOT->GetMacroPath() );
   //     gROOT->ProcessLine(".L TMVAMultiClassGui.C");


   //---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAClass" << std::endl;

   // Create a new root output file.
   TString outfileName = "TMVAClass.root";
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                                  "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");

   dataloader->AddSpectator("Btag_R2", 'F');
   dataloader->AddSpectator("Btag_thrustBm", 'F');
   dataloader->AddVariable("Btag_thrustOm", 'F');
   dataloader->AddVariable("Btag_cosTBTO", 'F');
   dataloader->AddSpectator("Btag_cosTBz", 'F');
   dataloader->AddSpectator("Btag_KSFWVariables_et", 'F');
   dataloader->AddSpectator("Btag_KSFWVariables_mm2", 'F');
   dataloader->AddVariable("Btag_KSFWVariables_hso00", 'F');
   dataloader->AddSpectator("Btag_KSFWVariables_hso01", 'F');
   dataloader->AddSpectator("Btag_KSFWVariables_hso02", 'F');
   dataloader->AddSpectator("Btag_KSFWVariables_hso03", 'F');
   dataloader->AddSpectator("Btag_KSFWVariables_hso04", 'F');
   dataloader->AddVariable("Btag_KSFWVariables_hso10", 'F');
   dataloader->AddSpectator("Btag_KSFWVariables_hso12", 'F');
   dataloader->AddSpectator("Btag_KSFWVariables_hso14", 'F');
   dataloader->AddVariable("Btag_KSFWVariables_hso20", 'F');
   dataloader->AddSpectator("Btag_KSFWVariables_hso22", 'F');
   dataloader->AddSpectator("Btag_KSFWVariables_hso24", 'F');
   dataloader->AddSpectator("Btag_KSFWVariables_hoo0", 'F');
   dataloader->AddSpectator("Btag_KSFWVariables_hoo1", 'F');
   dataloader->AddSpectator("Btag_KSFWVariables_hoo2", 'F');
   dataloader->AddSpectator("Btag_KSFWVariables_hoo3", 'F');
   dataloader->AddSpectator("Btag_KSFWVariables_hoo4", 'F');
   dataloader->AddSpectator("Btag_CleoConeCS_1", 'F');
   dataloader->AddSpectator("Btag_CleoConeCS_2", 'F');
   dataloader->AddSpectator("Btag_CleoConeCS_3", 'F');
   dataloader->AddSpectator("Btag_CleoConeCS_4", 'F');
   dataloader->AddSpectator("Btag_CleoConeCS_5", 'F');
   dataloader->AddSpectator("Btag_CleoConeCS_6", 'F');
   dataloader->AddSpectator("Btag_CleoConeCS_7", 'F');
   dataloader->AddSpectator("Btag_CleoConeCS_8", 'F');
   dataloader->AddSpectator("Btag_CleoConeCS_9", 'F');
   dataloader->AddVariable("Btag_useCMSFrame_theta", 'F');
   dataloader->AddVariable("Btag_extraInfo_SignalProbability", 'F');
   dataloader->AddVariable("Btag_chiProb", 'F');
   //dataloader->AddSpectator("extraInfo__boDecayHash__bc", 'F');
   //dataloader->AddSpectator("extraInfo__boDecayHashExtended__bc", 'F');

   dataloader->AddSpectator("__experiment__", 'I');
   dataloader->AddSpectator("__run__", 'I');
   dataloader->AddSpectator("__event__", 'I');
   dataloader->AddSpectator("__candidate__", 'I');
   dataloader->AddSpectator("__ncandidates__", 'I');
   dataloader->AddSpectator("extraInfo__bodecayModeID__bc", 'F');
   dataloader->AddSpectator("Bsig_daughter_0_extraInfo_decayModeID", 'F');
   dataloader->AddSpectator("nROE_ECLClusters__bocleanMask__bc", 'F');
   dataloader->AddSpectator("nROE_KLMClusters", 'F');
   dataloader->AddSpectator("nROE_Tracks__bocleanMask__bc", 'F');
   dataloader->AddSpectator("roeEextra__bocleanMask__bc", 'F');
   dataloader->AddSpectator("nROE_NeutralECLClusters__bocleanMask__bc", 'F');
   dataloader->AddSpectator("nROE_ParticlesInList__boK_S0__clmyKaonshort__bc", 'F');
   dataloader->AddSpectator("nROE_ParticlesInList__bopi0__clmyneutralPion__bc", 'F');
   dataloader->AddSpectator("missingMomentumOfEvent_theta", 'F');
   dataloader->AddSpectator("missingMomentumOfEvent", 'F');
   dataloader->AddSpectator("missingEnergyOfEventCMS", 'F');
   dataloader->AddSpectator("nRemainingTracksInEvent", 'F');
   dataloader->AddSpectator("roeNeextra__bocleanMask__bc", 'F');
   dataloader->AddSpectator("useCMSFrame__boroeNeextra__bocleanMask__bc__bc", 'F');
   dataloader->AddSpectator("nROE_ParticlesInList__bogamma__clmygamma__bc", 'F');
   dataloader->AddSpectator("Bsig_E", 'F');
   dataloader->AddSpectator("Bsig_useCMSFrame_E", 'F');
   dataloader->AddSpectator("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc", 'F');
   dataloader->AddSpectator("Bsig_p", 'F');
   dataloader->AddSpectator("Bsig_useCMSFrame_p", 'F');
   dataloader->AddSpectator("useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc", 'F');
   dataloader->AddSpectator("Bsig_M", 'F');
   dataloader->AddSpectator("Bsig_daughter_0_extraInfo_Dcvetomass", 'F');
   dataloader->AddSpectator("Bsig_daughter_0_extraInfo_DcvetodmID", 'F');
   dataloader->AddSpectator("Bsig_daughter_0_extraInfo_Dcvetoabsdm", 'F');
   dataloader->AddSpectator("Bsig_daughter_0_extraInfo_Dnvetomass", 'F');
   dataloader->AddSpectator("Bsig_daughter_0_extraInfo_DnvetodmID", 'F');
   dataloader->AddSpectator("Bsig_daughter_0_extraInfo_Dnvetoabsdm", 'F');
   dataloader->AddSpectator("Bsig_daughter_0_extraInfo_KS0_3D_distance", 'F');
   dataloader->AddSpectator("Bsig_daughter_0_extraInfo_Kaon_PID_correction", 'F');
   dataloader->AddSpectator("Bsig_daughter_0_extraInfo_Kaon_PID_rel_uncer_up", 'F');
   dataloader->AddSpectator("Bsig_daughter_0_extraInfo_Kaon_PID_rel_uncer_dn", 'F');
   dataloader->AddSpectator("Bsig_daughter_0_extraInfo_Pion_PID_correction", 'F');
   dataloader->AddSpectator("Bsig_daughter_0_extraInfo_Pion_PID_rel_uncer_up", 'F');
   dataloader->AddSpectator("Bsig_daughter_0_extraInfo_Pion_PID_rel_uncer_dn", 'F');
   dataloader->AddSpectator("Bsig_daughter_0_extraInfo_nKslow1", 'F');
   dataloader->AddSpectator("Bsig_daughter_0_extraInfo_nKslow2", 'F');
   dataloader->AddSpectator("Bsig_daughter_0_extraInfo_nKslow3", 'F');
   dataloader->AddSpectator("Bsig_daughter_0_extraInfo_nPislow1", 'F');
   dataloader->AddSpectator("Bsig_daughter_0_extraInfo_nPislow2", 'F');
   dataloader->AddSpectator("Bsig_daughter_0_extraInfo_nPislow3", 'F');
   dataloader->AddSpectator("Bsig_daughter_0_extraInfo_nKexcep", 'F');
   dataloader->AddSpectator("Bsig_daughter_0_extraInfo_nPiexcep", 'F');
   dataloader->AddSpectator("Btag_extraInfo_decayModeID", 'F');
   dataloader->AddSpectator("Btag_Mbc", 'F');
   dataloader->AddSpectator("Btag_deltaE", 'F');
   dataloader->AddSpectator("Btag_E", 'F');
   dataloader->AddSpectator("missingMass2OfEvent", 'F');
   dataloader->AddSpectator("visibleEnergyOfEventCMS", 'F');
   dataloader->AddSpectator("Btag_useCMSFrame_E", 'F');
   dataloader->AddSpectator("flag", 'I');

   TFile *SIGNAL_input(0);
   TFile *CHG_input(0);
   TFile *MIX_input(0);
   TFile* SIGNAL_input_test(0);
   TFile* CHG_input_test(0);
   TFile* MIX_input_test(0);

   SIGNAL_input = TFile::Open( "/media/sf_virtualbox_folder/20211109/For_TMVA_SIGNAL/train/Mxs_large/final_output_merge_Mxs_larger_SIGNAL_train_data.root" );
   CHG_input = TFile::Open("/media/sf_virtualbox_folder/20211109/For_TMVA_BKG/train/Mxs_large/final_output_merge_Mxs_larger_CHG_train_data.root");
   MIX_input = TFile::Open("/media/sf_virtualbox_folder/20211109/For_TMVA_BKG/train/Mxs_large/final_output_merge_Mxs_larger_MIX_train_data.root");
   SIGNAL_input_test = TFile::Open("/media/sf_virtualbox_folder/20211109/For_TMVA_SIGNAL/test/Mxs_large/final_output_merge_Mxs_larger_SIGNAL_test_data.root");
   CHG_input_test = TFile::Open("/media/sf_virtualbox_folder/20211109/For_TMVA_BKG/test/Mxs_large/final_output_merge_Mxs_larger_CHG_test_data.root");
   MIX_input_test = TFile::Open("/media/sf_virtualbox_folder/20211109/For_TMVA_BKG/test/Mxs_large/final_output_merge_Mxs_larger_MIX_test_data.root");

   TTree* SIGNAL_Tree = (TTree*)SIGNAL_input->Get("data");
   TTree* CHG_Tree = (TTree*)CHG_input->Get("data");
   TTree* MIX_Tree = (TTree*)MIX_input->Get("data");
   TTree* SIGNAL_Tree_test = (TTree*)SIGNAL_input_test->Get("data");
   TTree* CHG_Tree_test = (TTree*)CHG_input_test->Get("data");
   TTree* MIX_Tree_test = (TTree*)MIX_input_test->Get("data");

   gROOT->cd( outfileName+TString(":/") );
   dataloader->AddSignalTree    (SIGNAL_Tree, 1.0, TMVA::Types::kTraining);
   dataloader->AddBackgroundTree    (CHG_Tree, 1.0, TMVA::Types::kTraining);
   dataloader->AddBackgroundTree    (MIX_Tree, 1.0, TMVA::Types::kTraining);
   dataloader->AddSignalTree(SIGNAL_Tree_test, 1.0, TMVA::Types::kTesting);
   dataloader->AddBackgroundTree(CHG_Tree_test, 1.0, TMVA::Types::kTesting);
   dataloader->AddBackgroundTree(MIX_Tree_test, 1.0, TMVA::Types::kTesting);

   dataloader->PrepareTrainingAndTestTree( "", "SplitMode=Random:NormMode=NumEvents:!V" );

   factory->BookMethod(dataloader, TMVA::Types::kMLP, "MLP", "H:!V:VarTransform=N");

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClass is done!" << std::endl;

   delete factory;
   delete dataloader;

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVAGui( outfileName );


}