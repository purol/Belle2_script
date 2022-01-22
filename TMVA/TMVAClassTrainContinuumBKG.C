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

void TMVAClassTrainContinuumBKG( TString myMethodList = "" )
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
   dataloader->AddVariable("nROE_Tracks__bocleanMask__bc", 'F');
   dataloader->AddSpectator("roeEextra__bocleanMask__bc", 'F');
   dataloader->AddSpectator("nROE_NeutralECLClusters__bocleanMask__bc", 'F');
   dataloader->AddSpectator("nROE_ParticlesInList__boK_S0__clmyKaonshort__bc", 'F');
   dataloader->AddSpectator("nROE_ParticlesInList__bopi0__clmyneutralPion__bc", 'F');
   dataloader->AddSpectator("missingMomentumOfEvent_theta", 'F');
   dataloader->AddSpectator("missingMomentumOfEvent", 'F');
   dataloader->AddSpectator("missingEnergyOfEventCMS", 'F');
   dataloader->AddVariable("nRemainingTracksInEvent", 'F');
   dataloader->AddSpectator("roeNeextra__bocleanMask__bc", 'F');
   dataloader->AddSpectator("useCMSFrame__boroeNeextra__bocleanMask__bc__bc", 'F');
   dataloader->AddSpectator("nROE_ParticlesInList__bogamma__clmygamma__bc", 'F');
   dataloader->AddSpectator("foxWolframR1", 'F');
   dataloader->AddSpectator("foxWolframR2", 'F');
   dataloader->AddSpectator("foxWolframR3", 'F');
   dataloader->AddSpectator("foxWolframR4", 'F');
   dataloader->AddSpectator("harmonicMomentThrust0", 'F');
   dataloader->AddSpectator("harmonicMomentThrust1", 'F');
   dataloader->AddSpectator("harmonicMomentThrust2", 'F');
   dataloader->AddSpectator("harmonicMomentThrust3", 'F');
   dataloader->AddSpectator("harmonicMomentThrust4", 'F');
   dataloader->AddSpectator("cleoConeThrust0", 'F');
   dataloader->AddSpectator("cleoConeThrust1", 'F');
   dataloader->AddSpectator("cleoConeThrust2", 'F');
   dataloader->AddSpectator("cleoConeThrust3", 'F');
   dataloader->AddSpectator("cleoConeThrust4", 'F');
   dataloader->AddSpectator("cleoConeThrust5", 'F');
   dataloader->AddSpectator("cleoConeThrust6", 'F');
   dataloader->AddSpectator("cleoConeThrust7", 'F');
   dataloader->AddSpectator("cleoConeThrust8", 'F');
   dataloader->AddVariable("sphericity", 'F');
   dataloader->AddVariable("aplanarity", 'F');
   dataloader->AddSpectator("thrust", 'F');
   dataloader->AddSpectator("thrustAxisCosTheta", 'F');
   dataloader->AddSpectator("MsquaredBsig_op0", 'F');
   dataloader->AddSpectator("MsquaredBsig_op1", 'F');
   dataloader->AddSpectator("MsquaredBsig_op2", 'F');
   dataloader->AddSpectator("MsquaredBsig_op3", 'F');
   dataloader->AddSpectator("MsquaredBsig_op4", 'F');
   dataloader->AddSpectator("MsquaredBsig_op7", 'F');
   dataloader->AddVariable("roeP__bocleanMask__bc", 'F');
   dataloader->AddVariable("roeM__bocleanMask__bc", 'F');
   dataloader->AddVariable("roePTheta__bocleanMask__bc", 'F');
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
   dataloader->AddSpectator("Bsig_R2", 'F');
   dataloader->AddSpectator("Bsig_thrustBm", 'F');
   dataloader->AddSpectator("Bsig_thrustOm", 'F');
   dataloader->AddSpectator("Bsig_cosTBTO", 'F');
   dataloader->AddSpectator("Bsig_cosTBz", 'F');
   dataloader->AddSpectator("Bsig_KSFWVariables_et", 'F');
   dataloader->AddSpectator("Bsig_KSFWVariables_mm2", 'F');
   dataloader->AddSpectator("Bsig_KSFWVariables_hso00", 'F');
   dataloader->AddSpectator("Bsig_KSFWVariables_hso01", 'F');
   dataloader->AddSpectator("Bsig_KSFWVariables_hso02", 'F');
   dataloader->AddSpectator("Bsig_KSFWVariables_hso03", 'F');
   dataloader->AddSpectator("Bsig_KSFWVariables_hso04", 'F');
   dataloader->AddSpectator("Bsig_KSFWVariables_hso10", 'F');
   dataloader->AddSpectator("Bsig_KSFWVariables_hso12", 'F');
   dataloader->AddSpectator("Bsig_KSFWVariables_hso14", 'F');
   dataloader->AddSpectator("Bsig_KSFWVariables_hso20", 'F');
   dataloader->AddSpectator("Bsig_KSFWVariables_hso22", 'F');
   dataloader->AddSpectator("Bsig_KSFWVariables_hso24", 'F');
   dataloader->AddSpectator("Bsig_KSFWVariables_hoo0", 'F');
   dataloader->AddSpectator("Bsig_KSFWVariables_hoo1", 'F');
   dataloader->AddSpectator("Bsig_KSFWVariables_hoo2", 'F');
   dataloader->AddSpectator("Bsig_KSFWVariables_hoo3", 'F');
   dataloader->AddSpectator("Bsig_KSFWVariables_hoo4", 'F');
   dataloader->AddSpectator("Bsig_CleoConeCS_1", 'F');
   dataloader->AddSpectator("Bsig_CleoConeCS_2", 'F');
   dataloader->AddSpectator("Bsig_CleoConeCS_3", 'F');
   dataloader->AddSpectator("Bsig_CleoConeCS_4", 'F');
   dataloader->AddSpectator("Bsig_CleoConeCS_5", 'F');
   dataloader->AddSpectator("Bsig_CleoConeCS_6", 'F');
   dataloader->AddSpectator("Bsig_CleoConeCS_7", 'F');
   dataloader->AddSpectator("Bsig_CleoConeCS_8", 'F');
   dataloader->AddSpectator("Bsig_CleoConeCS_9", 'F');
   dataloader->AddSpectator("Btag_extraInfo_decayModeID", 'F');
   dataloader->AddSpectator("Btag_Mbc", 'F');
   dataloader->AddVariable("Btag_deltaE", 'F');
   dataloader->AddSpectator("Btag_E", 'F');
   dataloader->AddSpectator("missingMass2OfEvent", 'F');
   dataloader->AddSpectator("visibleEnergyOfEventCMS", 'F');
   dataloader->AddSpectator("Btag_useCMSFrame_E", 'F');
   dataloader->AddSpectator("flag", 'I');

   TFile *SIGNAL_input(0);
   TFile *UUBAR_input(0);
   TFile *DDBAR_input(0);
   TFile *SSBAR_input(0);
   TFile* CHARM_input(0);
   TFile* SIGNAL_input_test(0);
   TFile* UUBAR_input_test(0);
   TFile* DDBAR_input_test(0);
   TFile* SSBAR_input_test(0);
   TFile* CHARM_input_test(0);

   SIGNAL_input = TFile::Open( "/media/sf_virtualbox_folder/20211109/For_TMVA_SIGNAL/train/Mxs_large/final_output_merge_Mxs_larger_SIGNAL_train_data.root" );
   UUBAR_input = TFile::Open("/media/sf_virtualbox_folder/20211109/For_TMVA_BKG/train/Mxs_large/final_output_merge_Mxs_larger_UUBAR_train_data.root");
   DDBAR_input = TFile::Open("/media/sf_virtualbox_folder/20211109/For_TMVA_BKG/train/Mxs_large/final_output_merge_Mxs_larger_DDBAR_train_data.root");
   SSBAR_input = TFile::Open("/media/sf_virtualbox_folder/20211109/For_TMVA_BKG/train/Mxs_large/final_output_merge_Mxs_larger_SSBAR_train_data.root");
   CHARM_input = TFile::Open("/media/sf_virtualbox_folder/20211109/For_TMVA_BKG/train/Mxs_large/final_output_merge_Mxs_larger_CHARM_train_data.root");
   SIGNAL_input_test = TFile::Open("/media/sf_virtualbox_folder/20211109/For_TMVA_SIGNAL/test/Mxs_large/final_output_merge_Mxs_larger_SIGNAL_test_data.root");
   UUBAR_input_test = TFile::Open("/media/sf_virtualbox_folder/20211109/For_TMVA_BKG/test/Mxs_large/final_output_merge_Mxs_larger_UUBAR_test_data.root");
   DDBAR_input_test = TFile::Open("/media/sf_virtualbox_folder/20211109/For_TMVA_BKG/test/Mxs_large/final_output_merge_Mxs_larger_DDBAR_test_data.root");
   SSBAR_input_test = TFile::Open("/media/sf_virtualbox_folder/20211109/For_TMVA_BKG/test/Mxs_large/final_output_merge_Mxs_larger_SSBAR_test_data.root");
   CHARM_input_test = TFile::Open("/media/sf_virtualbox_folder/20211109/For_TMVA_BKG/test/Mxs_large/final_output_merge_Mxs_larger_CHARM_test_data.root");

   TTree* SIGNAL_Tree = (TTree*)SIGNAL_input->Get("data");
   TTree* UUBAR_Tree = (TTree*)UUBAR_input->Get("data");
   TTree* DDBAR_Tree = (TTree*)DDBAR_input->Get("data");
   TTree* SSBAR_Tree = (TTree*)SSBAR_input->Get("data");
   TTree* CHARM_Tree = (TTree*)CHARM_input->Get("data");
   TTree* SIGNAL_Tree_test = (TTree*)SIGNAL_input_test->Get("data");
   TTree* UUBAR_Tree_test = (TTree*)UUBAR_input_test->Get("data");
   TTree* DDBAR_Tree_test = (TTree*)DDBAR_input_test->Get("data");
   TTree* SSBAR_Tree_test = (TTree*)SSBAR_input_test->Get("data");
   TTree* CHARM_Tree_test = (TTree*)CHARM_input_test->Get("data");

   gROOT->cd( outfileName+TString(":/") );
   dataloader->AddSignalTree    (SIGNAL_Tree, 1.0, TMVA::Types::kTraining);
   dataloader->AddBackgroundTree    (UUBAR_Tree, 1.0, TMVA::Types::kTraining);
   dataloader->AddBackgroundTree    (DDBAR_Tree, 1.0, TMVA::Types::kTraining);
   dataloader->AddBackgroundTree    (SSBAR_Tree, 1.0, TMVA::Types::kTraining);
   dataloader->AddBackgroundTree(CHARM_Tree, 1.0, TMVA::Types::kTraining);
   dataloader->AddSignalTree(SIGNAL_Tree_test, 1.0, TMVA::Types::kTesting);
   dataloader->AddBackgroundTree(UUBAR_Tree_test, 1.0, TMVA::Types::kTesting);
   dataloader->AddBackgroundTree(DDBAR_Tree_test, 1.0, TMVA::Types::kTesting);
   dataloader->AddBackgroundTree(SSBAR_Tree_test, 1.0, TMVA::Types::kTesting);
   dataloader->AddBackgroundTree(CHARM_Tree_test, 1.0, TMVA::Types::kTesting);

   dataloader->PrepareTrainingAndTestTree( "", "SplitMode=Random:NormMode=NumEvents:!V" );

   factory->BookMethod(dataloader, TMVA::Types::kBDT, "MLP", "H:!V");

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