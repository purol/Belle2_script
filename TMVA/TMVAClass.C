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

void TMVAClass( TString myMethodList = "" )
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

   dataloader->AddVariable("Btag_R2", 'F');
   dataloader->AddVariable("Btag_thrustBm", 'F');
   dataloader->AddVariable("Btag_thrustOm", 'F');
   dataloader->AddVariable("Btag_cosTBTO", 'F');
   dataloader->AddVariable("Btag_cosTBz", 'F');
   dataloader->AddVariable("Btag_KSFWVariables_et", 'F');
   dataloader->AddVariable("Btag_KSFWVariables_mm2", 'F');
   dataloader->AddVariable("Btag_KSFWVariables_hso00", 'F');
   dataloader->AddVariable("Btag_KSFWVariables_hso02", 'F');
   dataloader->AddVariable("Btag_KSFWVariables_hso04", 'F');
   dataloader->AddVariable("Btag_KSFWVariables_hso10", 'F');
   dataloader->AddVariable("Btag_KSFWVariables_hso12", 'F');
   dataloader->AddVariable("Btag_KSFWVariables_hso14", 'F');
   dataloader->AddVariable("Btag_KSFWVariables_hso20", 'F');
   dataloader->AddVariable("Btag_KSFWVariables_hso22", 'F');
   dataloader->AddVariable("Btag_KSFWVariables_hso24", 'F');
   dataloader->AddVariable("Btag_CleoConeCS_1", 'F');
   dataloader->AddVariable("Btag_CleoConeCS_2", 'F');
   dataloader->AddVariable("Btag_CleoConeCS_3", 'F');
   dataloader->AddVariable("Btag_CleoConeCS_4", 'F');
   dataloader->AddVariable("Btag_CleoConeCS_5", 'F');
   dataloader->AddVariable("Btag_CleoConeCS_6", 'F');
   dataloader->AddVariable("Btag_CleoConeCS_7", 'F');
   dataloader->AddVariable("Btag_CleoConeCS_8", 'F');
   dataloader->AddVariable("Btag_CleoConeCS_9", 'F');
   dataloader->AddVariable("missingMass2OfEvent", 'F');
   dataloader->AddVariable("visibleEnergyOfEventCMS", 'F');
   dataloader->AddVariable("Btag_useCMSFrame_theta", 'F');
   dataloader->AddVariable("Btag_extraInfo_SignalProbability", 'F');

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
   dataloader->AddSpectator("Bsig_E", 'F');
   dataloader->AddSpectator("Bsig_useCMSFrame_E", 'F');
   dataloader->AddSpectator("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc", 'F');
   dataloader->AddSpectator("Bsig_p", 'F');
   dataloader->AddSpectator("Bsig_useCMSFrame_p", 'F');
   dataloader->AddSpectator("useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc", 'F');
   dataloader->AddSpectator("Bsig_M", 'F');
   dataloader->AddSpectator("Btag_extraInfo_decayModeID", 'F');
   dataloader->AddSpectator("Btag_Mbc", 'F');
   dataloader->AddSpectator("Btag_deltaE", 'F');
   dataloader->AddSpectator("Btag_E", 'F');
   dataloader->AddSpectator("Btag_useCMSFrame_E", 'F');
   dataloader->AddSpectator("Btag_chiProb", 'F');
   dataloader->AddSpectator("flag", 'I');

   dataloader->AddSpectator("nParticlesInList__boB__pl__clKcharge_total__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsu__clMCcomb__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsu__clMCch1__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsu__clMCch2__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsu__clMCch3__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsu__clMCch4__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsu__clMCch5__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsu__clMCch6__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsu__clMCch7__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsu__clMCch8__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsu__clMCch9__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsu__clMCch10__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsu__clMCch11__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsu__clMCch12__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsu__clMCch13__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsu__clMCch14__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsu__clMCch15__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boB0__clKneutral_total__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsd__clMCcomb__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsd__clMCch16__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsd__clMCch17__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsd__clMCch18__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsd__clMCch19__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsd__clMCch20__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsd__clMCch21__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsd__clMCch22__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsd__clMCch23__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsd__clMCch24__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsd__clMCch25__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsd__clMCch26__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsd__clMCch27__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsd__clMCch28__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsd__clMCch29__bc", 'F');
   dataloader->AddSpectator("nParticlesInList__boXsd__clMCch30__bc", 'F');

   TFile *SIGNAL_input(0);
   TFile *CHG_input(0);
   TFile *MIX_input(0);
   TFile *UUBAR_input(0);
   TFile *DDBAR_input(0);
   TFile *SSBAR_input(0);
   TFile* CHARM_input(0);

   SIGNAL_input = TFile::Open( "SIGNAL_output_merge_data.root" );
   CHG_input = TFile::Open("CHG_output_merge_data.root");
   MIX_input = TFile::Open("MIX_output_merge_data.root");
   UUBAR_input = TFile::Open("UUBAR_output_merge_data.root");
   DDBAR_input = TFile::Open("DDBAR_output_merge_data.root");
   SSBAR_input = TFile::Open("SSBAR_output_merge_data.root");
   CHARM_input = TFile::Open("CHARM_output_merge_data.root");

   TTree* SIGNAL_Tree = (TTree*)SIGNAL_input->Get("data");
   TTree* CHG_Tree = (TTree*)CHG_input->Get("data");
   TTree* MIX_Tree = (TTree*)MIX_input->Get("data");
   TTree* UUBAR_Tree = (TTree*)UUBAR_input->Get("data");
   TTree* DDBAR_Tree = (TTree*)DDBAR_input->Get("data");
   TTree* SSBAR_Tree = (TTree*)SSBAR_input->Get("data");
   TTree* CHARM_Tree = (TTree*)CHARM_input->Get("data");

   gROOT->cd( outfileName+TString(":/") );
   dataloader->AddSignalTree    (SIGNAL_Tree, 1.0 );
   dataloader->AddBackgroundTree    (CHG_Tree, 1.0);
   dataloader->AddBackgroundTree    (MIX_Tree, 1.0);
   dataloader->AddBackgroundTree    (UUBAR_Tree, 1.0);
   dataloader->AddBackgroundTree    (DDBAR_Tree, 1.0);
   dataloader->AddBackgroundTree    (SSBAR_Tree, 1.0);
   dataloader->AddBackgroundTree(CHARM_Tree, 1.0);

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

int main( int argc, char** argv )
{
   // Select methods (don't look at this code - not of interest)
   TString methodList;
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   TMVAClass(methodList);
   return 0;
}

