// last update: 2023-07-14
// for Belle2 data

#include <iostream>
#include <stdio.h>
#include <string>
#include <algorithm>
#include <vector>
#include <numeric>
#include <time.h>
#include <random>
#include <sstream>
#include <queue>

#include "TH1.h"
#include "TH2.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TCollection.h"
#include "TPad.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TMath.h"

/*
when you add new variables:
revise # define N_Needed_info ...
revise typedef struct data
revise void Loader::GetData(TFile* input_file)
revise void Loader::PrintRootFile(std::string output_name)
revise void Loader::PrintSeperateRootFile(std::string output_name)
revise void Loader::ConvertIntoSeparateRootFile(std::string output_name, double flag = 0)
*/

# define N_Needed_info 37
//# define N_event_info 15
# define N_Upsilon_info 164
# define N_Bsig_info 738
# define N_Btag_info 11
# define N_decay 38 // five decay mode + others
# define N_decay_nparticles 5 // # of nu_e, B, B0
# define N_decay_syst_ff 7 // helicity angle + q2
# define N_PID_syst 73
# define N_fakeE_syst 37
# define N_fakeMU_syst 49
# define N_pi0_syst 8

// small: temp_BB_output > 0.6 && temp_Continuum_output > 0.94
// large: temp_BB_output > 0.6 && temp_Continuum_output > 0.84

void load_files(const char* dirname, std::vector<std::string>* names) {
    TSystemDirectory dir(dirname, dirname);
    TList* files = dir.GetListOfFiles();
    if (files) {
        TSystemFile* file;
        TString fname;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(".root")) {
                names->push_back(fname.Data());
            }
        }
    }
}

void main(int argc, char* argv[]){
        /*
        * argv[1]: dirname
        * argv[2]: OBB
        * argv[3]: OContinuum
        * argv[4]: output path
        */

        double OBB = std::stod(argv[2]);
        double OContinuum = std::stod(argv[3]);

        std::vector<std::string> names;
        load_files(argv[1], &names);

        for (unsigned int i = 0; i < names.size(); i++) {
            std::string string_filename = argv[1] + std::string("/") + names.at(i);
            const char* filename = string_filename.c_str();
            std::string OnlyFileName = string_filename.substr(string_filename.find_last_of("\\/") + 1, string_filename.size() - string_filename.find_last_of("\\/"));
            size_t lastindex = OnlyFileName.find_last_of(".");
            std::string rawname = OnlyFileName.substr(0, lastindex);

            TFile* input_file = new TFile(filename, "read");

            TTree* temp_tree = (TTree*)input_file->Get("data");

            int temp_ExperimentToTree;
            int temp_RunToTree;
            unsigned int temp_EventToTree;
            int temp_CandidateToTree;
            int temp_NcandidatesToTree;
            double temp_UpsilonDataToTree[N_Upsilon_info];
            double temp_BsigDataToTree[N_Bsig_info];
            double temp_BtagDataToTree[N_Btag_info];
            double temp_DataToTree[N_Needed_info];
            int temp_DecayDataToTree[N_decay];
            int temp_DecayNparticlesDataToTree[N_decay_nparticles];
            double temp_DecaySystFFDataToTree[N_decay_syst_ff];
            double temp_Upsilon_decayIDToTree;
            double temp_Bsig_decayIDToTree;
            int temp_flag;
            float temp_BB_output;
            float temp_Continuum_output;
            int temp_nROE_mu;
            int temp_nROE_Jpsi;
            int temp_nROE_Upsilon_BCS;
            int temp_nROE_Upsilon;

            bool DoesItHaveXsBranch = false;
            for (int i = 0; i < temp_tree->GetListOfBranches()->LastIndex(); i++) {
                if (temp_tree->GetListOfBranches()->At(i)->GetName() == std::string("nParticlesInList__boB__pl__clKcharge_total__bc")) DoesItHaveXsBranch = true;
            }

            bool DoesItHaveJpsiOutput = false;
            if (temp_tree->FindLeaf("nParticlesInList__bomu__pl__clfromUpsilonmychargedMuon__bc") == 0 || temp_tree->FindLeaf("nParticlesInList__boJ__slpsi__clfromUpsilontemp__bc") == 0 || temp_tree->FindLeaf("nParticlesInList__boUpsilon__bo4S__bc__cltemp_withoutneutrino_BCS__bc") == 0 || temp_tree->FindLeaf("nParticlesInList__boUpsilon__bo4S__bc__cltemp_withoutneutrino__bc") == 0) DoesItHaveJpsiOutput = false;
            else DoesItHaveJpsiOutput = true;

            // Get data
            /*================================================================*/
            // get event_info
            temp_tree->SetBranchAddress("__experiment__", &temp_ExperimentToTree);
            temp_tree->SetBranchAddress("__run__", &temp_RunToTree);
            temp_tree->SetBranchAddress("__event__", &temp_EventToTree);
            temp_tree->SetBranchAddress("__candidate__", &temp_CandidateToTree);
            temp_tree->SetBranchAddress("__ncandidates__", &temp_NcandidatesToTree);

            // get decaymodeID
            temp_tree->SetBranchAddress("extraInfo__bodecayModeID__bc", &temp_Upsilon_decayIDToTree);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &temp_Bsig_decayIDToTree);

            // get Upsilon_info
            temp_tree->SetBranchAddress("nROE_ECLClusters__bocleanMask__bc", &temp_UpsilonDataToTree[0]);
            temp_tree->SetBranchAddress("nROE_KLMClusters", &temp_UpsilonDataToTree[1]);
            temp_tree->SetBranchAddress("nROE_Tracks__bocleanMask__bc", &temp_UpsilonDataToTree[2]);
            temp_tree->SetBranchAddress("roeEextra__bocleanMask__bc", &temp_UpsilonDataToTree[3]);
            temp_tree->SetBranchAddress("nROE_NeutralECLClusters__bocleanMask__bc", &temp_UpsilonDataToTree[4]);
            temp_tree->SetBranchAddress("nROE_ParticlesInList__boK_S0__clmyKaonshort__bc", &temp_UpsilonDataToTree[5]);
            temp_tree->SetBranchAddress("nROE_ParticlesInList__bopi0__clmyneutralPion__bc", &temp_UpsilonDataToTree[6]);
            temp_tree->SetBranchAddress("missingMomentumOfEvent_theta", &temp_UpsilonDataToTree[7]);
            temp_tree->SetBranchAddress("missingMomentumOfEvent", &temp_UpsilonDataToTree[8]);
            temp_tree->SetBranchAddress("missingEnergyOfEventCMS", &temp_UpsilonDataToTree[9]);
            temp_tree->SetBranchAddress("nRemainingTracksInEvent", &temp_UpsilonDataToTree[10]);
            temp_tree->SetBranchAddress("roeNeextra__bocleanMask__bc", &temp_UpsilonDataToTree[11]);
            temp_tree->SetBranchAddress("useCMSFrame__boroeNeextra__bocleanMask__bc__bc", &temp_UpsilonDataToTree[12]);
            temp_tree->SetBranchAddress("nROE_ParticlesInList__bogamma__clmygamma__bc", &temp_UpsilonDataToTree[13]);
            temp_tree->SetBranchAddress("foxWolframR1", &temp_UpsilonDataToTree[14]);
            temp_tree->SetBranchAddress("foxWolframR2", &temp_UpsilonDataToTree[15]);
            temp_tree->SetBranchAddress("foxWolframR3", &temp_UpsilonDataToTree[16]);
            temp_tree->SetBranchAddress("foxWolframR4", &temp_UpsilonDataToTree[17]);
            temp_tree->SetBranchAddress("harmonicMomentThrust0", &temp_UpsilonDataToTree[18]);
            temp_tree->SetBranchAddress("harmonicMomentThrust1", &temp_UpsilonDataToTree[19]);
            temp_tree->SetBranchAddress("harmonicMomentThrust2", &temp_UpsilonDataToTree[20]);
            temp_tree->SetBranchAddress("harmonicMomentThrust3", &temp_UpsilonDataToTree[21]);
            temp_tree->SetBranchAddress("harmonicMomentThrust4", &temp_UpsilonDataToTree[22]);
            temp_tree->SetBranchAddress("cleoConeThrust0", &temp_UpsilonDataToTree[23]);
            temp_tree->SetBranchAddress("cleoConeThrust1", &temp_UpsilonDataToTree[24]);
            temp_tree->SetBranchAddress("cleoConeThrust2", &temp_UpsilonDataToTree[25]);
            temp_tree->SetBranchAddress("cleoConeThrust3", &temp_UpsilonDataToTree[26]);
            temp_tree->SetBranchAddress("cleoConeThrust4", &temp_UpsilonDataToTree[27]);
            temp_tree->SetBranchAddress("cleoConeThrust5", &temp_UpsilonDataToTree[28]);
            temp_tree->SetBranchAddress("cleoConeThrust6", &temp_UpsilonDataToTree[29]);
            temp_tree->SetBranchAddress("cleoConeThrust7", &temp_UpsilonDataToTree[30]);
            temp_tree->SetBranchAddress("cleoConeThrust8", &temp_UpsilonDataToTree[31]);
            temp_tree->SetBranchAddress("sphericity", &temp_UpsilonDataToTree[32]);
            temp_tree->SetBranchAddress("aplanarity", &temp_UpsilonDataToTree[33]);
            temp_tree->SetBranchAddress("thrust", &temp_UpsilonDataToTree[34]);
            temp_tree->SetBranchAddress("thrustAxisCosTheta", &temp_UpsilonDataToTree[35]);
            temp_tree->SetBranchAddress("MsquaredBsig_op0", &temp_UpsilonDataToTree[36]);
            temp_tree->SetBranchAddress("MsquaredBsig_op1", &temp_UpsilonDataToTree[37]);
            temp_tree->SetBranchAddress("MsquaredBsig_op2", &temp_UpsilonDataToTree[38]);
            temp_tree->SetBranchAddress("MsquaredBsig_op3", &temp_UpsilonDataToTree[39]);
            temp_tree->SetBranchAddress("MsquaredBsig_op4", &temp_UpsilonDataToTree[40]);
            temp_tree->SetBranchAddress("MsquaredBsig_op7", &temp_UpsilonDataToTree[41]);
            temp_tree->SetBranchAddress("roeP__bocleanMask__bc", &temp_UpsilonDataToTree[42]);
            temp_tree->SetBranchAddress("roeM__bocleanMask__bc", &temp_UpsilonDataToTree[43]);
            temp_tree->SetBranchAddress("roePTheta__bocleanMask__bc", &temp_UpsilonDataToTree[44]);
            temp_tree->SetBranchAddress("qsquared", &temp_UpsilonDataToTree[45]);
            temp_tree->SetBranchAddress("chiProb", &temp_UpsilonDataToTree[46]);
            temp_tree->SetBranchAddress("dr", &temp_UpsilonDataToTree[47]);
            temp_tree->SetBranchAddress("dz", &temp_UpsilonDataToTree[48]);
            temp_tree->SetBranchAddress("nParticlesInList__boe__pl__clElectronFBDT__bc", &temp_UpsilonDataToTree[49]);
            temp_tree->SetBranchAddress("nParticlesInList__bomu__pl__clMuonFBDT__bc", &temp_UpsilonDataToTree[50]);
            temp_tree->SetBranchAddress("nParticlesInList__boe__pl__clElectronFBDT_loose__bc", &temp_UpsilonDataToTree[51]);
            temp_tree->SetBranchAddress("nParticlesInList__bomu__pl__clMuonFBDT_loose__bc", &temp_UpsilonDataToTree[52]);
            temp_tree->SetBranchAddress("nParticlesInList__boe__pl__clElectronFBDT_tight__bc", &temp_UpsilonDataToTree[53]);
            temp_tree->SetBranchAddress("nParticlesInList__bomu__pl__clMuonFBDT_tight__bc", &temp_UpsilonDataToTree[54]);
            temp_tree->SetBranchAddress("beamE", &temp_UpsilonDataToTree[55]);
            temp_tree->SetBranchAddress("nROE_Tracks__bolooseMask__bc", &temp_UpsilonDataToTree[56]);
            temp_tree->SetBranchAddress("Ecms", &temp_UpsilonDataToTree[57]);
            temp_tree->SetBranchAddress("extraInfo__boNgamma__bc", &temp_UpsilonDataToTree[58]);
            temp_tree->SetBranchAddress("extraInfo__boEecl__bc", &temp_UpsilonDataToTree[59]);
            temp_tree->SetBranchAddress("extraInfo__boEecl_matched__bc", &temp_UpsilonDataToTree[60]);
            temp_tree->SetBranchAddress("extraInfo__boEecl_unmatched__bc", &temp_UpsilonDataToTree[61]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav133__bc", &temp_UpsilonDataToTree[62]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv133__bc", &temp_UpsilonDataToTree[63]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv133_matched__bc", &temp_UpsilonDataToTree[64]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv133_unmatched__bc", &temp_UpsilonDataToTree[65]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200__bc", &temp_UpsilonDataToTree[66]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_matched__bc", &temp_UpsilonDataToTree[67]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_unmatched__bc", &temp_UpsilonDataToTree[68]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200__bc", &temp_UpsilonDataToTree[69]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_matched__bc", &temp_UpsilonDataToTree[70]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_unmatched__bc", &temp_UpsilonDataToTree[71]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_800__bc", &temp_UpsilonDataToTree[72]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_800_matched__bc", &temp_UpsilonDataToTree[73]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_800_unmatched__bc", &temp_UpsilonDataToTree[74]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_800__bc", &temp_UpsilonDataToTree[75]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_800_matched__bc", &temp_UpsilonDataToTree[76]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_800_unmatched__bc", &temp_UpsilonDataToTree[77]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_825__bc", &temp_UpsilonDataToTree[78]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_825_matched__bc", &temp_UpsilonDataToTree[79]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_825_unmatched__bc", &temp_UpsilonDataToTree[80]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_825__bc", &temp_UpsilonDataToTree[81]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_825_matched__bc", &temp_UpsilonDataToTree[82]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_825_unmatched__bc", &temp_UpsilonDataToTree[83]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_850__bc", &temp_UpsilonDataToTree[84]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_850_matched__bc", &temp_UpsilonDataToTree[85]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_850_unmatched__bc", &temp_UpsilonDataToTree[86]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_850__bc", &temp_UpsilonDataToTree[87]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_850_matched__bc", &temp_UpsilonDataToTree[88]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_850_unmatched__bc", &temp_UpsilonDataToTree[89]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_875__bc", &temp_UpsilonDataToTree[90]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_875_matched__bc", &temp_UpsilonDataToTree[91]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_875_unmatched__bc", &temp_UpsilonDataToTree[92]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_875__bc", &temp_UpsilonDataToTree[93]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_875_matched__bc", &temp_UpsilonDataToTree[94]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_875_unmatched__bc", &temp_UpsilonDataToTree[95]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_900__bc", &temp_UpsilonDataToTree[96]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_900_matched__bc", &temp_UpsilonDataToTree[97]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_900_unmatched__bc", &temp_UpsilonDataToTree[98]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_900__bc", &temp_UpsilonDataToTree[99]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_900_matched__bc", &temp_UpsilonDataToTree[100]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_900_unmatched__bc", &temp_UpsilonDataToTree[101]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_925__bc", &temp_UpsilonDataToTree[102]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_925_matched__bc", &temp_UpsilonDataToTree[103]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_925_unmatched__bc", &temp_UpsilonDataToTree[104]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_925__bc", &temp_UpsilonDataToTree[105]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_925_matched__bc", &temp_UpsilonDataToTree[106]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_925_unmatched__bc", &temp_UpsilonDataToTree[107]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_950__bc", &temp_UpsilonDataToTree[108]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_950_matched__bc", &temp_UpsilonDataToTree[109]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_950_unmatched__bc", &temp_UpsilonDataToTree[110]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_950__bc", &temp_UpsilonDataToTree[111]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_950_matched__bc", &temp_UpsilonDataToTree[112]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_950_unmatched__bc", &temp_UpsilonDataToTree[113]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_975__bc", &temp_UpsilonDataToTree[114]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_975_matched__bc", &temp_UpsilonDataToTree[115]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_975_unmatched__bc", &temp_UpsilonDataToTree[116]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_975__bc", &temp_UpsilonDataToTree[117]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_975_matched__bc", &temp_UpsilonDataToTree[118]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_975_unmatched__bc", &temp_UpsilonDataToTree[119]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_025__bc", &temp_UpsilonDataToTree[120]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_025_matched__bc", &temp_UpsilonDataToTree[121]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_025_unmatched__bc", &temp_UpsilonDataToTree[122]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_025__bc", &temp_UpsilonDataToTree[123]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_025_matched__bc", &temp_UpsilonDataToTree[124]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_025_unmatched__bc", &temp_UpsilonDataToTree[125]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_050__bc", &temp_UpsilonDataToTree[126]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_050_matched__bc", &temp_UpsilonDataToTree[127]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_050_unmatched__bc", &temp_UpsilonDataToTree[128]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_050__bc", &temp_UpsilonDataToTree[129]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_050_matched__bc", &temp_UpsilonDataToTree[130]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_050_unmatched__bc", &temp_UpsilonDataToTree[131]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_075__bc", &temp_UpsilonDataToTree[132]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_075_matched__bc", &temp_UpsilonDataToTree[133]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_075_unmatched__bc", &temp_UpsilonDataToTree[134]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_075__bc", &temp_UpsilonDataToTree[135]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_075_matched__bc", &temp_UpsilonDataToTree[136]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_075_unmatched__bc", &temp_UpsilonDataToTree[137]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_100__bc", &temp_UpsilonDataToTree[138]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_100_matched__bc", &temp_UpsilonDataToTree[139]);
            temp_tree->SetBranchAddress("extraInfo__boNgammav200_100_unmatched__bc", &temp_UpsilonDataToTree[140]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_100__bc", &temp_UpsilonDataToTree[141]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_100_matched__bc", &temp_UpsilonDataToTree[142]);
            temp_tree->SetBranchAddress("extraInfo__boEeclv200_100_unmatched__bc", &temp_UpsilonDataToTree[143]);
            temp_tree->SetBranchAddress("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &temp_UpsilonDataToTree[144]);
            temp_tree->SetBranchAddress("nParticlesInList__boD0__clDecayIntoKL0__bc", &temp_UpsilonDataToTree[145]);
            temp_tree->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &temp_UpsilonDataToTree[146]);
            temp_tree->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &temp_UpsilonDataToTree[147]);
            temp_tree->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &temp_UpsilonDataToTree[148]);
            temp_tree->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &temp_UpsilonDataToTree[149]);
            temp_tree->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &temp_UpsilonDataToTree[150]);
            temp_tree->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &temp_UpsilonDataToTree[151]);
            temp_tree->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &temp_UpsilonDataToTree[152]);
            temp_tree->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &temp_UpsilonDataToTree[153]);
            temp_tree->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_all__bc", &temp_UpsilonDataToTree[154]);
            temp_tree->SetBranchAddress("nParticlesInList__boB__pl__clKpKLKL_NR__bc", &temp_UpsilonDataToTree[155]);
            temp_tree->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &temp_UpsilonDataToTree[156]);
            temp_tree->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &temp_UpsilonDataToTree[157]);
            temp_tree->SetBranchAddress("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &temp_UpsilonDataToTree[158]);
            temp_tree->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_all__bc", &temp_UpsilonDataToTree[159]);
            temp_tree->SetBranchAddress("nParticlesInList__boB0__clKSKLKL_NR__bc", &temp_UpsilonDataToTree[160]);
            temp_tree->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &temp_UpsilonDataToTree[161]);
            temp_tree->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &temp_UpsilonDataToTree[162]);
            temp_tree->SetBranchAddress("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &temp_UpsilonDataToTree[163]);

            // get Bsig_info
            temp_tree->SetBranchAddress("Bsig_E", &temp_BsigDataToTree[0]);
            temp_tree->SetBranchAddress("Bsig_useCMSFrame_E", &temp_BsigDataToTree[1]);
            temp_tree->SetBranchAddress("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc", &temp_BsigDataToTree[2]);
            temp_tree->SetBranchAddress("Bsig_p", &temp_BsigDataToTree[3]);
            temp_tree->SetBranchAddress("Bsig_useCMSFrame_p", &temp_BsigDataToTree[4]);
            temp_tree->SetBranchAddress("useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc", &temp_BsigDataToTree[5]);
            temp_tree->SetBranchAddress("Bsig_M", &temp_BsigDataToTree[6]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_3D_distance", &temp_BsigDataToTree[7]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_nKslow1", &temp_BsigDataToTree[8]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_nKslow2", &temp_BsigDataToTree[9]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_nKslow3", &temp_BsigDataToTree[10]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_nPislow1", &temp_BsigDataToTree[11]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_nPislow2", &temp_BsigDataToTree[12]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_nPislow3", &temp_BsigDataToTree[13]);
            temp_tree->SetBranchAddress("Bsig_R2", &temp_BsigDataToTree[14]);
            temp_tree->SetBranchAddress("Bsig_thrustBm", &temp_BsigDataToTree[15]);
            temp_tree->SetBranchAddress("Bsig_thrustOm", &temp_BsigDataToTree[16]);
            temp_tree->SetBranchAddress("Bsig_cosTBTO", &temp_BsigDataToTree[17]);
            temp_tree->SetBranchAddress("Bsig_cosTBz", &temp_BsigDataToTree[18]);
            temp_tree->SetBranchAddress("Bsig_KSFWVariables_et", &temp_BsigDataToTree[19]);
            temp_tree->SetBranchAddress("Bsig_KSFWVariables_mm2", &temp_BsigDataToTree[20]);
            temp_tree->SetBranchAddress("Bsig_KSFWVariables_hso00", &temp_BsigDataToTree[21]);
            temp_tree->SetBranchAddress("Bsig_KSFWVariables_hso01", &temp_BsigDataToTree[22]);
            temp_tree->SetBranchAddress("Bsig_KSFWVariables_hso02", &temp_BsigDataToTree[23]);
            temp_tree->SetBranchAddress("Bsig_KSFWVariables_hso03", &temp_BsigDataToTree[24]);
            temp_tree->SetBranchAddress("Bsig_KSFWVariables_hso04", &temp_BsigDataToTree[25]);
            temp_tree->SetBranchAddress("Bsig_KSFWVariables_hso10", &temp_BsigDataToTree[26]);
            temp_tree->SetBranchAddress("Bsig_KSFWVariables_hso12", &temp_BsigDataToTree[27]);
            temp_tree->SetBranchAddress("Bsig_KSFWVariables_hso14", &temp_BsigDataToTree[28]);
            temp_tree->SetBranchAddress("Bsig_KSFWVariables_hso20", &temp_BsigDataToTree[29]);
            temp_tree->SetBranchAddress("Bsig_KSFWVariables_hso22", &temp_BsigDataToTree[30]);
            temp_tree->SetBranchAddress("Bsig_KSFWVariables_hso24", &temp_BsigDataToTree[31]);
            temp_tree->SetBranchAddress("Bsig_KSFWVariables_hoo0", &temp_BsigDataToTree[32]);
            temp_tree->SetBranchAddress("Bsig_KSFWVariables_hoo1", &temp_BsigDataToTree[33]);
            temp_tree->SetBranchAddress("Bsig_KSFWVariables_hoo2", &temp_BsigDataToTree[34]);
            temp_tree->SetBranchAddress("Bsig_KSFWVariables_hoo3", &temp_BsigDataToTree[35]);
            temp_tree->SetBranchAddress("Bsig_KSFWVariables_hoo4", &temp_BsigDataToTree[36]);
            temp_tree->SetBranchAddress("Bsig_CleoConeCS_1", &temp_BsigDataToTree[37]);
            temp_tree->SetBranchAddress("Bsig_CleoConeCS_2", &temp_BsigDataToTree[38]);
            temp_tree->SetBranchAddress("Bsig_CleoConeCS_3", &temp_BsigDataToTree[39]);
            temp_tree->SetBranchAddress("Bsig_CleoConeCS_4", &temp_BsigDataToTree[40]);
            temp_tree->SetBranchAddress("Bsig_CleoConeCS_5", &temp_BsigDataToTree[41]);
            temp_tree->SetBranchAddress("Bsig_CleoConeCS_6", &temp_BsigDataToTree[42]);
            temp_tree->SetBranchAddress("Bsig_CleoConeCS_7", &temp_BsigDataToTree[43]);
            temp_tree->SetBranchAddress("Bsig_CleoConeCS_8", &temp_BsigDataToTree[44]);
            temp_tree->SetBranchAddress("Bsig_CleoConeCS_9", &temp_BsigDataToTree[45]);
            temp_tree->SetBranchAddress("Bsig_pt", &temp_BsigDataToTree[46]);
            temp_tree->SetBranchAddress("Bsig_useCMSFrame_pt", &temp_BsigDataToTree[47]);
            temp_tree->SetBranchAddress("Bsig_theta", &temp_BsigDataToTree[48]);
            temp_tree->SetBranchAddress("Bsig_useCMSFrame_theta", &temp_BsigDataToTree[49]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_nDc", &temp_BsigDataToTree[50]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_Dc_pValue_med", &temp_BsigDataToTree[51]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_Dc_pValue_std", &temp_BsigDataToTree[52]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_chiProb", &temp_BsigDataToTree[53]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_dr", &temp_BsigDataToTree[54]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_dz", &temp_BsigDataToTree[55]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcsimpleveto_M", &temp_BsigDataToTree[56]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_nD0", &temp_BsigDataToTree[57]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_D0_pValue_med", &temp_BsigDataToTree[58]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_D0_pValue_std", &temp_BsigDataToTree[59]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_D0simpleveto_chiProb", &temp_BsigDataToTree[60]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_D0simpleveto_dr", &temp_BsigDataToTree[61]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_D0simpleveto_dz", &temp_BsigDataToTree[62]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_D0simpleveto_M", &temp_BsigDataToTree[63]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_mychiProb", &temp_BsigDataToTree[64]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_mydr", &temp_BsigDataToTree[65]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_mydz", &temp_BsigDataToTree[66]);
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                temp_tree->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKtruebin" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[67 + 4 * i_PID]);
                temp_tree->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKmisbin" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[67 + 4 * i_PID + 1]);
                temp_tree->SetBranchAddress(("Bsig_daughter_0_extraInfo_npitruebin" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[67 + 4 * i_PID + 2]);
                temp_tree->SetBranchAddress(("Bsig_daughter_0_extraInfo_npimisbin" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[67 + 4 * i_PID + 3]);
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) temp_tree->SetBranchAddress(("Bsig_daughter_0_extraInfo_npi0bin" + std::to_string(i_pi0)).c_str(), &temp_BsigDataToTree[359 + i_pi0]);
            for (int i_PID = 0; i_PID < N_fakeE_syst; i_PID++) {
                temp_tree->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeEbin_n" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[367 + 4 * i_PID]);
                temp_tree->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeEbin_p" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[367 + 4 * i_PID + 1]);
                temp_tree->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeEbin_n" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[367 + 4 * i_PID + 2]);
                temp_tree->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeEbin_p" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[367 + 4 * i_PID + 3]);
            }
            for (int i_PID = 0; i_PID < N_fakeMU_syst; i_PID++) {
                temp_tree->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeMUbin_n" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[542 + 4 * i_PID]);
                temp_tree->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeMUbin_p" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[542 + 4 * i_PID + 1]);
                temp_tree->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeMUbin_n" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[542 + 4 * i_PID + 2]);
                temp_tree->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeMUbin_p" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[542 + 4 * i_PID + 3]);
            }

            // get Btag_info
            temp_tree->SetBranchAddress("Btag_extraInfo_decayModeID", &temp_BtagDataToTree[0]);
            temp_tree->SetBranchAddress("Btag_Mbc", &temp_BtagDataToTree[1]);
            temp_tree->SetBranchAddress("Btag_deltaE", &temp_BtagDataToTree[2]);
            temp_tree->SetBranchAddress("Btag_E", &temp_BtagDataToTree[3]);
            temp_tree->SetBranchAddress("Btag_useCMSFrame_E", &temp_BtagDataToTree[4]);
            temp_tree->SetBranchAddress("Btag_extraInfo_SignalProbability", &temp_BtagDataToTree[5]);
            temp_tree->SetBranchAddress("Btag_chiProb", &temp_BtagDataToTree[6]);
            temp_tree->SetBranchAddress("Btag_dr", &temp_BtagDataToTree[7]);
            temp_tree->SetBranchAddress("Btag_dz", &temp_BtagDataToTree[8]);
            temp_tree->SetBranchAddress("Btag_useCMSFrame_p", &temp_BtagDataToTree[9]);
            temp_tree->SetBranchAddress("Btag_useCMSFrame_phi", &temp_BtagDataToTree[10]);

            // other information I need
            temp_tree->SetBranchAddress("Btag_R2", &temp_DataToTree[0]);
            temp_tree->SetBranchAddress("Btag_thrustBm", &temp_DataToTree[1]);
            temp_tree->SetBranchAddress("Btag_thrustOm", &temp_DataToTree[2]);
            temp_tree->SetBranchAddress("Btag_cosTBTO", &temp_DataToTree[3]);
            temp_tree->SetBranchAddress("Btag_cosTBz", &temp_DataToTree[4]);
            temp_tree->SetBranchAddress("Btag_KSFWVariables_et", &temp_DataToTree[5]);
            temp_tree->SetBranchAddress("Btag_KSFWVariables_mm2", &temp_DataToTree[6]);
            temp_tree->SetBranchAddress("Btag_KSFWVariables_hso00", &temp_DataToTree[7]);
            temp_tree->SetBranchAddress("Btag_KSFWVariables_hso01", &temp_DataToTree[8]);
            temp_tree->SetBranchAddress("Btag_KSFWVariables_hso02", &temp_DataToTree[9]);
            temp_tree->SetBranchAddress("Btag_KSFWVariables_hso03", &temp_DataToTree[10]);
            temp_tree->SetBranchAddress("Btag_KSFWVariables_hso04", &temp_DataToTree[11]);
            temp_tree->SetBranchAddress("Btag_KSFWVariables_hso10", &temp_DataToTree[12]);
            temp_tree->SetBranchAddress("Btag_KSFWVariables_hso12", &temp_DataToTree[13]);
            temp_tree->SetBranchAddress("Btag_KSFWVariables_hso14", &temp_DataToTree[14]);
            temp_tree->SetBranchAddress("Btag_KSFWVariables_hso20", &temp_DataToTree[15]);
            temp_tree->SetBranchAddress("Btag_KSFWVariables_hso22", &temp_DataToTree[16]);
            temp_tree->SetBranchAddress("Btag_KSFWVariables_hso24", &temp_DataToTree[17]);
            temp_tree->SetBranchAddress("Btag_KSFWVariables_hoo0", &temp_DataToTree[18]);
            temp_tree->SetBranchAddress("Btag_KSFWVariables_hoo1", &temp_DataToTree[19]);
            temp_tree->SetBranchAddress("Btag_KSFWVariables_hoo2", &temp_DataToTree[20]);
            temp_tree->SetBranchAddress("Btag_KSFWVariables_hoo3", &temp_DataToTree[21]);
            temp_tree->SetBranchAddress("Btag_KSFWVariables_hoo4", &temp_DataToTree[22]);
            temp_tree->SetBranchAddress("Btag_CleoConeCS_1", &temp_DataToTree[23]);
            temp_tree->SetBranchAddress("Btag_CleoConeCS_2", &temp_DataToTree[24]);
            temp_tree->SetBranchAddress("Btag_CleoConeCS_3", &temp_DataToTree[25]);
            temp_tree->SetBranchAddress("Btag_CleoConeCS_4", &temp_DataToTree[26]);
            temp_tree->SetBranchAddress("Btag_CleoConeCS_5", &temp_DataToTree[27]);
            temp_tree->SetBranchAddress("Btag_CleoConeCS_6", &temp_DataToTree[28]);
            temp_tree->SetBranchAddress("Btag_CleoConeCS_7", &temp_DataToTree[29]);
            temp_tree->SetBranchAddress("Btag_CleoConeCS_8", &temp_DataToTree[30]);
            temp_tree->SetBranchAddress("Btag_CleoConeCS_9", &temp_DataToTree[31]);
            temp_tree->SetBranchAddress("missingMass2OfEvent", &temp_DataToTree[32]);
            temp_tree->SetBranchAddress("visibleEnergyOfEventCMS", &temp_DataToTree[33]);
            temp_tree->SetBranchAddress("Btag_useCMSFrame_theta", &temp_DataToTree[34]);
            temp_tree->SetBranchAddress("extraInfo__boDecayHash__bc", &temp_DataToTree[35]);
            temp_tree->SetBranchAddress("extraInfo__boDecayHashExtended__bc", &temp_DataToTree[36]);

            if (DoesItHaveXsBranch) {
                // decay mode (MC level)
                temp_tree->SetBranchAddress("nParticlesInList__boB__pl__clKcharge_total__bc", &temp_DecayDataToTree[0]);
                temp_tree->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &temp_DecayDataToTree[1]);
                temp_tree->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &temp_DecayDataToTree[2]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsu__clMCcomb__bc", &temp_DecayDataToTree[3]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsu__clMCch1__bc", &temp_DecayDataToTree[4]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsu__clMCch2__bc", &temp_DecayDataToTree[5]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsu__clMCch3__bc", &temp_DecayDataToTree[6]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsu__clMCch4__bc", &temp_DecayDataToTree[7]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsu__clMCch5__bc", &temp_DecayDataToTree[8]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsu__clMCch6__bc", &temp_DecayDataToTree[9]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsu__clMCch7__bc", &temp_DecayDataToTree[10]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsu__clMCch8__bc", &temp_DecayDataToTree[11]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsu__clMCch9__bc", &temp_DecayDataToTree[12]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsu__clMCch10__bc", &temp_DecayDataToTree[13]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsu__clMCch11__bc", &temp_DecayDataToTree[14]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsu__clMCch12__bc", &temp_DecayDataToTree[15]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsu__clMCch13__bc", &temp_DecayDataToTree[16]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsu__clMCch14__bc", &temp_DecayDataToTree[17]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsu__clMCch15__bc", &temp_DecayDataToTree[18]);
                temp_tree->SetBranchAddress("nParticlesInList__boB0__clKneutral_total__bc", &temp_DecayDataToTree[19]);
                temp_tree->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &temp_DecayDataToTree[20]);
                temp_tree->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &temp_DecayDataToTree[21]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsd__clMCcomb__bc", &temp_DecayDataToTree[22]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsd__clMCch16__bc", &temp_DecayDataToTree[23]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsd__clMCch17__bc", &temp_DecayDataToTree[24]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsd__clMCch18__bc", &temp_DecayDataToTree[25]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsd__clMCch19__bc", &temp_DecayDataToTree[26]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsd__clMCch20__bc", &temp_DecayDataToTree[27]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsd__clMCch21__bc", &temp_DecayDataToTree[28]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsd__clMCch22__bc", &temp_DecayDataToTree[29]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsd__clMCch23__bc", &temp_DecayDataToTree[30]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsd__clMCch24__bc", &temp_DecayDataToTree[31]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsd__clMCch25__bc", &temp_DecayDataToTree[32]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsd__clMCch26__bc", &temp_DecayDataToTree[33]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsd__clMCch27__bc", &temp_DecayDataToTree[34]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsd__clMCch28__bc", &temp_DecayDataToTree[35]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsd__clMCch29__bc", &temp_DecayDataToTree[36]);
                temp_tree->SetBranchAddress("nParticlesInList__boXsd__clMCch30__bc", &temp_DecayDataToTree[37]);
                temp_tree->SetBranchAddress("nParticlesInList__bonu_e__clMC_signal__bc", &temp_DecayNparticlesDataToTree[0]);
                temp_tree->SetBranchAddress("nParticlesInList__boB__pl__clMC_signal_total_e__bc", &temp_DecayNparticlesDataToTree[1]);
                temp_tree->SetBranchAddress("nParticlesInList__boB0__clMC_signal_total_e__bc", &temp_DecayNparticlesDataToTree[2]);
                temp_tree->SetBranchAddress("nParticlesInList__boB__pl__clPrimaryMC__bc", &temp_DecayNparticlesDataToTree[3]);
                temp_tree->SetBranchAddress("nParticlesInList__boB0__clPrimaryMC__bc", &temp_DecayNparticlesDataToTree[4]);
                temp_tree->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &temp_DecaySystFFDataToTree[0]);
                temp_tree->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", &temp_DecaySystFFDataToTree[1]);
                temp_tree->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", &temp_DecaySystFFDataToTree[2]);
                temp_tree->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &temp_DecaySystFFDataToTree[3]);
                temp_tree->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &temp_DecaySystFFDataToTree[4]);
                temp_tree->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spM__bc", &temp_DecaySystFFDataToTree[5]);
                temp_tree->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spM__bc", &temp_DecaySystFFDataToTree[6]);
            }

            // flag
            temp_tree->SetBranchAddress("flag", &temp_flag);

            // MLP
            temp_tree->SetBranchAddress("MVA_BB", &temp_BB_output);
            temp_tree->SetBranchAddress("MVA_Continuum", &temp_Continuum_output);

            if (DoesItHaveJpsiOutput) {
                temp_tree->SetBranchAddress("nParticlesInList__bomu__pl__clfromUpsilonmychargedMuon__bc", &temp_nROE_mu);
                temp_tree->SetBranchAddress("nParticlesInList__boJ__slpsi__clfromUpsilontemp__bc", &temp_nROE_Jpsi);
                temp_tree->SetBranchAddress("nParticlesInList__boUpsilon__bo4S__bc__cltemp_withoutneutrino_BCS__bc", &temp_nROE_Upsilon_BCS);
                temp_tree->SetBranchAddress("nParticlesInList__boUpsilon__bo4S__bc__cltemp_withoutneutrino__bc", &temp_nROE_Upsilon);
            }
            else {
                temp_nROE_mu = -1.0;
                temp_nROE_Jpsi = -1.0;
                temp_nROE_Upsilon_BCS = -1.0;
                temp_nROE_Upsilon = -1.0;
            }
            /*================================================================*/


            TFile* temp_file = new TFile((argv[4] + std::string("/") + rawname + "_after_MVA.root").c_str(), "recreate");
            temp_file->cd();
            TTree* temp_tree_upsilon = new TTree("Upsilon", "");
            TTree* temp_tree_Bsig = new TTree("Bsig", "");
            TTree* temp_tree_Btag = new TTree("Btag", "");
            TTree* temp_tree_Xs = nullptr;
            if (DoesItHaveXsBranch) temp_tree_Xs = new TTree("Xs", "");

            // print root file
                /*================================================================*/
            // get event_info
            temp_tree_upsilon->Branch("__experiment__", &temp_ExperimentToTree);
            temp_tree_upsilon->Branch("__run__", &temp_RunToTree);
            temp_tree_upsilon->Branch("__event__", &temp_EventToTree);
            temp_tree_upsilon->Branch("__candidate__", &temp_CandidateToTree);
            temp_tree_upsilon->Branch("__ncandidates__", &temp_NcandidatesToTree);
            temp_tree_Bsig->Branch("__experiment__", &temp_ExperimentToTree);
            temp_tree_Bsig->Branch("__run__", &temp_RunToTree);
            temp_tree_Bsig->Branch("__event__", &temp_EventToTree);
            temp_tree_Bsig->Branch("__candidate__", &temp_CandidateToTree);
            temp_tree_Bsig->Branch("__ncandidates__", &temp_NcandidatesToTree);
            temp_tree_Btag->Branch("__experiment__", &temp_ExperimentToTree);
            temp_tree_Btag->Branch("__run__", &temp_RunToTree);
            temp_tree_Btag->Branch("__event__", &temp_EventToTree);
            temp_tree_Btag->Branch("__candidate__", &temp_CandidateToTree);
            temp_tree_Btag->Branch("__ncandidates__", &temp_NcandidatesToTree);

            // get decaymodeID
            temp_tree_upsilon->Branch("extraInfo__bodecayModeID__bc", &temp_Upsilon_decayIDToTree);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_decayModeID", &temp_Bsig_decayIDToTree);

            // get Upsilon_info
            temp_tree_upsilon->Branch("nROE_ECLClusters__bocleanMask__bc", &temp_UpsilonDataToTree[0]);
            temp_tree_upsilon->Branch("nROE_KLMClusters", &temp_UpsilonDataToTree[1]);
            temp_tree_upsilon->Branch("nROE_Tracks__bocleanMask__bc", &temp_UpsilonDataToTree[2]);
            temp_tree_upsilon->Branch("roeEextra__bocleanMask__bc", &temp_UpsilonDataToTree[3]);
            temp_tree_upsilon->Branch("nROE_NeutralECLClusters__bocleanMask__bc", &temp_UpsilonDataToTree[4]);
            temp_tree_upsilon->Branch("nROE_ParticlesInList__boK_S0__clmyKaonshort__bc", &temp_UpsilonDataToTree[5]);
            temp_tree_upsilon->Branch("nROE_ParticlesInList__bopi0__clmyneutralPion__bc", &temp_UpsilonDataToTree[6]);
            temp_tree_upsilon->Branch("missingMomentumOfEvent_theta", &temp_UpsilonDataToTree[7]);
            temp_tree_upsilon->Branch("missingMomentumOfEvent", &temp_UpsilonDataToTree[8]);
            temp_tree_upsilon->Branch("missingEnergyOfEventCMS", &temp_UpsilonDataToTree[9]);
            temp_tree_upsilon->Branch("nRemainingTracksInEvent", &temp_UpsilonDataToTree[10]);
            temp_tree_upsilon->Branch("roeNeextra__bocleanMask__bc", &temp_UpsilonDataToTree[11]);
            temp_tree_upsilon->Branch("useCMSFrame__boroeNeextra__bocleanMask__bc__bc", &temp_UpsilonDataToTree[12]);
            temp_tree_upsilon->Branch("nROE_ParticlesInList__bogamma__clmygamma__bc", &temp_UpsilonDataToTree[13]);
            temp_tree_upsilon->Branch("foxWolframR1", &temp_UpsilonDataToTree[14]);
            temp_tree_upsilon->Branch("foxWolframR2", &temp_UpsilonDataToTree[15]);
            temp_tree_upsilon->Branch("foxWolframR3", &temp_UpsilonDataToTree[16]);
            temp_tree_upsilon->Branch("foxWolframR4", &temp_UpsilonDataToTree[17]);
            temp_tree_upsilon->Branch("harmonicMomentThrust0", &temp_UpsilonDataToTree[18]);
            temp_tree_upsilon->Branch("harmonicMomentThrust1", &temp_UpsilonDataToTree[19]);
            temp_tree_upsilon->Branch("harmonicMomentThrust2", &temp_UpsilonDataToTree[20]);
            temp_tree_upsilon->Branch("harmonicMomentThrust3", &temp_UpsilonDataToTree[21]);
            temp_tree_upsilon->Branch("harmonicMomentThrust4", &temp_UpsilonDataToTree[22]);
            temp_tree_upsilon->Branch("cleoConeThrust0", &temp_UpsilonDataToTree[23]);
            temp_tree_upsilon->Branch("cleoConeThrust1", &temp_UpsilonDataToTree[24]);
            temp_tree_upsilon->Branch("cleoConeThrust2", &temp_UpsilonDataToTree[25]);
            temp_tree_upsilon->Branch("cleoConeThrust3", &temp_UpsilonDataToTree[26]);
            temp_tree_upsilon->Branch("cleoConeThrust4", &temp_UpsilonDataToTree[27]);
            temp_tree_upsilon->Branch("cleoConeThrust5", &temp_UpsilonDataToTree[28]);
            temp_tree_upsilon->Branch("cleoConeThrust6", &temp_UpsilonDataToTree[29]);
            temp_tree_upsilon->Branch("cleoConeThrust7", &temp_UpsilonDataToTree[30]);
            temp_tree_upsilon->Branch("cleoConeThrust8", &temp_UpsilonDataToTree[31]);
            temp_tree_upsilon->Branch("sphericity", &temp_UpsilonDataToTree[32]);
            temp_tree_upsilon->Branch("aplanarity", &temp_UpsilonDataToTree[33]);
            temp_tree_upsilon->Branch("thrust", &temp_UpsilonDataToTree[34]);
            temp_tree_upsilon->Branch("thrustAxisCosTheta", &temp_UpsilonDataToTree[35]);
            temp_tree_upsilon->Branch("MsquaredBsig_op0", &temp_UpsilonDataToTree[36]);
            temp_tree_upsilon->Branch("MsquaredBsig_op1", &temp_UpsilonDataToTree[37]);
            temp_tree_upsilon->Branch("MsquaredBsig_op2", &temp_UpsilonDataToTree[38]);
            temp_tree_upsilon->Branch("MsquaredBsig_op3", &temp_UpsilonDataToTree[39]);
            temp_tree_upsilon->Branch("MsquaredBsig_op4", &temp_UpsilonDataToTree[40]);
            temp_tree_upsilon->Branch("MsquaredBsig_op7", &temp_UpsilonDataToTree[41]);
            temp_tree_upsilon->Branch("roeP__bocleanMask__bc", &temp_UpsilonDataToTree[42]);
            temp_tree_upsilon->Branch("roeM__bocleanMask__bc", &temp_UpsilonDataToTree[43]);
            temp_tree_upsilon->Branch("roePTheta__bocleanMask__bc", &temp_UpsilonDataToTree[44]);
            temp_tree_upsilon->Branch("qsquared", &temp_UpsilonDataToTree[45]);
            temp_tree_upsilon->Branch("chiProb", &temp_UpsilonDataToTree[46]);
            temp_tree_upsilon->Branch("dr", &temp_UpsilonDataToTree[47]);
            temp_tree_upsilon->Branch("dz", &temp_UpsilonDataToTree[48]);
            temp_tree_upsilon->Branch("nParticlesInList__boe__pl__clElectronFBDT__bc", &temp_UpsilonDataToTree[49]);
            temp_tree_upsilon->Branch("nParticlesInList__bomu__pl__clMuonFBDT__bc", &temp_UpsilonDataToTree[50]);
            temp_tree_upsilon->Branch("nParticlesInList__boe__pl__clElectronFBDT_loose__bc", &temp_UpsilonDataToTree[51]);
            temp_tree_upsilon->Branch("nParticlesInList__bomu__pl__clMuonFBDT_loose__bc", &temp_UpsilonDataToTree[52]);
            temp_tree_upsilon->Branch("nParticlesInList__boe__pl__clElectronFBDT_tight__bc", &temp_UpsilonDataToTree[53]);
            temp_tree_upsilon->Branch("nParticlesInList__bomu__pl__clMuonFBDT_tight__bc", &temp_UpsilonDataToTree[54]);
            temp_tree_upsilon->Branch("beamE", &temp_UpsilonDataToTree[55]);
            temp_tree_upsilon->Branch("nROE_Tracks__bolooseMask__bc", &temp_UpsilonDataToTree[56]);
            temp_tree_upsilon->Branch("Ecms", &temp_UpsilonDataToTree[57]);
            temp_tree_upsilon->Branch("extraInfo__boNgamma__bc", &temp_UpsilonDataToTree[58]);
            temp_tree_upsilon->Branch("extraInfo__boEecl__bc", &temp_UpsilonDataToTree[59]);
            temp_tree_upsilon->Branch("extraInfo__boEecl_matched__bc", &temp_UpsilonDataToTree[60]);
            temp_tree_upsilon->Branch("extraInfo__boEecl_unmatched__bc", &temp_UpsilonDataToTree[61]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav133__bc", &temp_UpsilonDataToTree[62]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv133__bc", &temp_UpsilonDataToTree[63]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv133_matched__bc", &temp_UpsilonDataToTree[64]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv133_unmatched__bc", &temp_UpsilonDataToTree[65]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200__bc", &temp_UpsilonDataToTree[66]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_matched__bc", &temp_UpsilonDataToTree[67]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_unmatched__bc", &temp_UpsilonDataToTree[68]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200__bc", &temp_UpsilonDataToTree[69]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_matched__bc", &temp_UpsilonDataToTree[70]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_unmatched__bc", &temp_UpsilonDataToTree[71]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_800__bc", &temp_UpsilonDataToTree[72]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_800_matched__bc", &temp_UpsilonDataToTree[73]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_800_unmatched__bc", &temp_UpsilonDataToTree[74]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_800__bc", &temp_UpsilonDataToTree[75]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_800_matched__bc", &temp_UpsilonDataToTree[76]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_800_unmatched__bc", &temp_UpsilonDataToTree[77]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_825__bc", &temp_UpsilonDataToTree[78]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_825_matched__bc", &temp_UpsilonDataToTree[79]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_825_unmatched__bc", &temp_UpsilonDataToTree[80]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_825__bc", &temp_UpsilonDataToTree[81]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_825_matched__bc", &temp_UpsilonDataToTree[82]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_825_unmatched__bc", &temp_UpsilonDataToTree[83]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_850__bc", &temp_UpsilonDataToTree[84]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_850_matched__bc", &temp_UpsilonDataToTree[85]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_850_unmatched__bc", &temp_UpsilonDataToTree[86]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_850__bc", &temp_UpsilonDataToTree[87]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_850_matched__bc", &temp_UpsilonDataToTree[88]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_850_unmatched__bc", &temp_UpsilonDataToTree[89]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_875__bc", &temp_UpsilonDataToTree[90]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_875_matched__bc", &temp_UpsilonDataToTree[91]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_875_unmatched__bc", &temp_UpsilonDataToTree[92]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_875__bc", &temp_UpsilonDataToTree[93]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_875_matched__bc", &temp_UpsilonDataToTree[94]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_875_unmatched__bc", &temp_UpsilonDataToTree[95]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_900__bc", &temp_UpsilonDataToTree[96]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_900_matched__bc", &temp_UpsilonDataToTree[97]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_900_unmatched__bc", &temp_UpsilonDataToTree[98]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_900__bc", &temp_UpsilonDataToTree[99]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_900_matched__bc", &temp_UpsilonDataToTree[100]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_900_unmatched__bc", &temp_UpsilonDataToTree[101]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_925__bc", &temp_UpsilonDataToTree[102]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_925_matched__bc", &temp_UpsilonDataToTree[103]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_925_unmatched__bc", &temp_UpsilonDataToTree[104]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_925__bc", &temp_UpsilonDataToTree[105]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_925_matched__bc", &temp_UpsilonDataToTree[106]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_925_unmatched__bc", &temp_UpsilonDataToTree[107]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_950__bc", &temp_UpsilonDataToTree[108]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_950_matched__bc", &temp_UpsilonDataToTree[109]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_950_unmatched__bc", &temp_UpsilonDataToTree[110]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_950__bc", &temp_UpsilonDataToTree[111]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_950_matched__bc", &temp_UpsilonDataToTree[112]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_950_unmatched__bc", &temp_UpsilonDataToTree[113]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_975__bc", &temp_UpsilonDataToTree[114]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_975_matched__bc", &temp_UpsilonDataToTree[115]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_975_unmatched__bc", &temp_UpsilonDataToTree[116]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_975__bc", &temp_UpsilonDataToTree[117]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_975_matched__bc", &temp_UpsilonDataToTree[118]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_975_unmatched__bc", &temp_UpsilonDataToTree[119]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_025__bc", &temp_UpsilonDataToTree[120]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_025_matched__bc", &temp_UpsilonDataToTree[121]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_025_unmatched__bc", &temp_UpsilonDataToTree[122]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_025__bc", &temp_UpsilonDataToTree[123]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_025_matched__bc", &temp_UpsilonDataToTree[124]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_025_unmatched__bc", &temp_UpsilonDataToTree[125]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_050__bc", &temp_UpsilonDataToTree[126]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_050_matched__bc", &temp_UpsilonDataToTree[127]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_050_unmatched__bc", &temp_UpsilonDataToTree[128]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_050__bc", &temp_UpsilonDataToTree[129]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_050_matched__bc", &temp_UpsilonDataToTree[130]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_050_unmatched__bc", &temp_UpsilonDataToTree[131]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_075__bc", &temp_UpsilonDataToTree[132]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_075_matched__bc", &temp_UpsilonDataToTree[133]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_075_unmatched__bc", &temp_UpsilonDataToTree[134]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_075__bc", &temp_UpsilonDataToTree[135]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_075_matched__bc", &temp_UpsilonDataToTree[136]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_075_unmatched__bc", &temp_UpsilonDataToTree[137]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_100__bc", &temp_UpsilonDataToTree[138]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_100_matched__bc", &temp_UpsilonDataToTree[139]);
            temp_tree_upsilon->Branch("extraInfo__boNgammav200_100_unmatched__bc", &temp_UpsilonDataToTree[140]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_100__bc", &temp_UpsilonDataToTree[141]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_100_matched__bc", &temp_UpsilonDataToTree[142]);
            temp_tree_upsilon->Branch("extraInfo__boEeclv200_100_unmatched__bc", &temp_UpsilonDataToTree[143]);
            temp_tree_upsilon->Branch("nParticlesInList__boD__pl__clDecayIntoKL0__bc", &temp_UpsilonDataToTree[144]);
            temp_tree_upsilon->Branch("nParticlesInList__boD0__clDecayIntoKL0__bc", &temp_UpsilonDataToTree[145]);
            temp_tree_upsilon->Branch("nParticlesInList__boB__pl__clKnn__bc", &temp_UpsilonDataToTree[146]);
            temp_tree_upsilon->Branch("invMassInLists__bon0__clKnn__bc", &temp_UpsilonDataToTree[147]);
            temp_tree_upsilon->Branch("nParticlesInList__boB__pl__clKstarnn__bc", &temp_UpsilonDataToTree[148]);
            temp_tree_upsilon->Branch("invMassInLists__bon0__clKstarnn__bc", &temp_UpsilonDataToTree[149]);
            temp_tree_upsilon->Branch("nParticlesInList__boB0__clK0nn__bc", &temp_UpsilonDataToTree[150]);
            temp_tree_upsilon->Branch("invMassInLists__bon0__clK0nn__bc", &temp_UpsilonDataToTree[151]);
            temp_tree_upsilon->Branch("nParticlesInList__boB0__clKstar0nn__bc", &temp_UpsilonDataToTree[152]);
            temp_tree_upsilon->Branch("invMassInLists__bon0__clKstar0nn__bc", &temp_UpsilonDataToTree[153]);
            temp_tree_upsilon->Branch("nParticlesInList__boB__pl__clKpKLKL_all__bc", &temp_UpsilonDataToTree[154]);
            temp_tree_upsilon->Branch("nParticlesInList__boB__pl__clKpKLKL_NR__bc", &temp_UpsilonDataToTree[155]);
            temp_tree_upsilon->Branch("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &temp_UpsilonDataToTree[156]);
            temp_tree_upsilon->Branch("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &temp_UpsilonDataToTree[157]);
            temp_tree_upsilon->Branch("averageValueInList__boB__pl__clKpKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &temp_UpsilonDataToTree[158]);
            temp_tree_upsilon->Branch("nParticlesInList__boB0__clKSKLKL_all__bc", &temp_UpsilonDataToTree[159]);
            temp_tree_upsilon->Branch("nParticlesInList__boB0__clKSKLKL_NR__bc", &temp_UpsilonDataToTree[160]);
            temp_tree_upsilon->Branch("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp1__bc__bc", &temp_UpsilonDataToTree[161]);
            temp_tree_upsilon->Branch("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo0__cm__sp2__bc__bc", &temp_UpsilonDataToTree[162]);
            temp_tree_upsilon->Branch("averageValueInList__boB0__clKSKLKL_NR__cm__spdaughterInvariantMass__bo1__cm__sp2__bc__bc", &temp_UpsilonDataToTree[163]);

            // get Bsig_info
            temp_tree_Bsig->Branch("Bsig_E", &temp_BsigDataToTree[0]);
            temp_tree_Bsig->Branch("Bsig_useCMSFrame_E", &temp_BsigDataToTree[1]);
            temp_tree_upsilon->Branch("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc", &temp_BsigDataToTree[2]);
            temp_tree_Bsig->Branch("Bsig_p", &temp_BsigDataToTree[3]);
            temp_tree_Bsig->Branch("Bsig_useCMSFrame_p", &temp_BsigDataToTree[4]);
            temp_tree_upsilon->Branch("useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc", &temp_BsigDataToTree[5]);
            temp_tree_Bsig->Branch("Bsig_M", &temp_BsigDataToTree[6]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_KS0_3D_distance", &temp_BsigDataToTree[7]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nKslow1", &temp_BsigDataToTree[8]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nKslow2", &temp_BsigDataToTree[9]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nKslow3", &temp_BsigDataToTree[10]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nPislow1", &temp_BsigDataToTree[11]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nPislow2", &temp_BsigDataToTree[12]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nPislow3", &temp_BsigDataToTree[13]);
            temp_tree_Bsig->Branch("Bsig_R2", &temp_BsigDataToTree[14]);
            temp_tree_Bsig->Branch("Bsig_thrustBm", &temp_BsigDataToTree[15]);
            temp_tree_Bsig->Branch("Bsig_thrustOm", &temp_BsigDataToTree[16]);
            temp_tree_Bsig->Branch("Bsig_cosTBTO", &temp_BsigDataToTree[17]);
            temp_tree_Bsig->Branch("Bsig_cosTBz", &temp_BsigDataToTree[18]);
            temp_tree_Bsig->Branch("Bsig_KSFWVariables_et", &temp_BsigDataToTree[19]);
            temp_tree_Bsig->Branch("Bsig_KSFWVariables_mm2", &temp_BsigDataToTree[20]);
            temp_tree_Bsig->Branch("Bsig_KSFWVariables_hso00", &temp_BsigDataToTree[21]);
            temp_tree_Bsig->Branch("Bsig_KSFWVariables_hso01", &temp_BsigDataToTree[22]);
            temp_tree_Bsig->Branch("Bsig_KSFWVariables_hso02", &temp_BsigDataToTree[23]);
            temp_tree_Bsig->Branch("Bsig_KSFWVariables_hso03", &temp_BsigDataToTree[24]);
            temp_tree_Bsig->Branch("Bsig_KSFWVariables_hso04", &temp_BsigDataToTree[25]);
            temp_tree_Bsig->Branch("Bsig_KSFWVariables_hso10", &temp_BsigDataToTree[26]);
            temp_tree_Bsig->Branch("Bsig_KSFWVariables_hso12", &temp_BsigDataToTree[27]);
            temp_tree_Bsig->Branch("Bsig_KSFWVariables_hso14", &temp_BsigDataToTree[28]);
            temp_tree_Bsig->Branch("Bsig_KSFWVariables_hso20", &temp_BsigDataToTree[29]);
            temp_tree_Bsig->Branch("Bsig_KSFWVariables_hso22", &temp_BsigDataToTree[30]);
            temp_tree_Bsig->Branch("Bsig_KSFWVariables_hso24", &temp_BsigDataToTree[31]);
            temp_tree_Bsig->Branch("Bsig_KSFWVariables_hoo0", &temp_BsigDataToTree[32]);
            temp_tree_Bsig->Branch("Bsig_KSFWVariables_hoo1", &temp_BsigDataToTree[33]);
            temp_tree_Bsig->Branch("Bsig_KSFWVariables_hoo2", &temp_BsigDataToTree[34]);
            temp_tree_Bsig->Branch("Bsig_KSFWVariables_hoo3", &temp_BsigDataToTree[35]);
            temp_tree_Bsig->Branch("Bsig_KSFWVariables_hoo4", &temp_BsigDataToTree[36]);
            temp_tree_Bsig->Branch("Bsig_CleoConeCS_1", &temp_BsigDataToTree[37]);
            temp_tree_Bsig->Branch("Bsig_CleoConeCS_2", &temp_BsigDataToTree[38]);
            temp_tree_Bsig->Branch("Bsig_CleoConeCS_3", &temp_BsigDataToTree[39]);
            temp_tree_Bsig->Branch("Bsig_CleoConeCS_4", &temp_BsigDataToTree[40]);
            temp_tree_Bsig->Branch("Bsig_CleoConeCS_5", &temp_BsigDataToTree[41]);
            temp_tree_Bsig->Branch("Bsig_CleoConeCS_6", &temp_BsigDataToTree[42]);
            temp_tree_Bsig->Branch("Bsig_CleoConeCS_7", &temp_BsigDataToTree[43]);
            temp_tree_Bsig->Branch("Bsig_CleoConeCS_8", &temp_BsigDataToTree[44]);
            temp_tree_Bsig->Branch("Bsig_CleoConeCS_9", &temp_BsigDataToTree[45]);
            temp_tree_Bsig->Branch("Bsig_pt", &temp_BsigDataToTree[46]);
            temp_tree_Bsig->Branch("Bsig_useCMSFrame_pt", &temp_BsigDataToTree[47]);
            temp_tree_Bsig->Branch("Bsig_theta", &temp_BsigDataToTree[48]);
            temp_tree_Bsig->Branch("Bsig_useCMSFrame_theta", &temp_BsigDataToTree[49]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nDc", &temp_BsigDataToTree[50]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Dc_pValue_med", &temp_BsigDataToTree[51]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Dc_pValue_std", &temp_BsigDataToTree[52]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Dcsimpleveto_chiProb", &temp_BsigDataToTree[53]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Dcsimpleveto_dr", &temp_BsigDataToTree[54]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Dcsimpleveto_dz", &temp_BsigDataToTree[55]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Dcsimpleveto_M", &temp_BsigDataToTree[56]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nD0", &temp_BsigDataToTree[57]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_D0_pValue_med", &temp_BsigDataToTree[58]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_D0_pValue_std", &temp_BsigDataToTree[59]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_D0simpleveto_chiProb", &temp_BsigDataToTree[60]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_D0simpleveto_dr", &temp_BsigDataToTree[61]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_D0simpleveto_dz", &temp_BsigDataToTree[62]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_D0simpleveto_M", &temp_BsigDataToTree[63]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_mychiProb", &temp_BsigDataToTree[64]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_mydr", &temp_BsigDataToTree[65]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_mydz", &temp_BsigDataToTree[66]);
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_nKtruebin" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[67 + 4 * i_PID]);
                temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_nKmisbin" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[67 + 4 * i_PID + 1]);
                temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_npitruebin" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[67 + 4 * i_PID + 2]);
                temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_npimisbin" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[67 + 4 * i_PID + 3]);
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_npi0bin" + std::to_string(i_pi0)).c_str(), &temp_BsigDataToTree[359 + i_pi0]);
            for (int i_PID = 0; i_PID < N_fakeE_syst; i_PID++) {
                temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_nKfakeEbin_n" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[367 + 4 * i_PID]);
                temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_nKfakeEbin_p" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[367 + 4 * i_PID + 1]);
                temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_npifakeEbin_n" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[367 + 4 * i_PID + 2]);
                temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_npifakeEbin_p" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[367 + 4 * i_PID + 3]);
            }
            for (int i_PID = 0; i_PID < N_fakeMU_syst; i_PID++) {
                temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_nKfakeMUbin_n" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[542 + 4 * i_PID]);
                temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_nKfakeMUbin_p" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[542 + 4 * i_PID + 1]);
                temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_npifakeMUbin_n" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[542 + 4 * i_PID + 2]);
                temp_tree_Bsig->Branch(("Bsig_daughter_0_extraInfo_npifakeMUbin_p" + std::to_string(i_PID)).c_str(), &temp_BsigDataToTree[542 + 4 * i_PID + 3]);
            }

            // get Btag_info
            temp_tree_Btag->Branch("Btag_extraInfo_decayModeID", &temp_BtagDataToTree[0]);
            temp_tree_Btag->Branch("Btag_Mbc", &temp_BtagDataToTree[1]);
            temp_tree_Btag->Branch("Btag_deltaE", &temp_BtagDataToTree[2]);
            temp_tree_Btag->Branch("Btag_E", &temp_BtagDataToTree[3]);
            temp_tree_Btag->Branch("Btag_useCMSFrame_E", &temp_BtagDataToTree[4]);
            temp_tree_Btag->Branch("Btag_extraInfo_SignalProbability", &temp_BtagDataToTree[5]);
            temp_tree_Btag->Branch("Btag_chiProb", &temp_BtagDataToTree[6]);
            temp_tree_Btag->Branch("Btag_dr", &temp_BtagDataToTree[7]);
            temp_tree_Btag->Branch("Btag_dz", &temp_BtagDataToTree[8]);
            temp_tree_Btag->Branch("Btag_useCMSFrame_p", &temp_BtagDataToTree[9]);
            temp_tree_Btag->Branch("Btag_useCMSFrame_phi", &temp_BtagDataToTree[10]);

            // other information I need
            temp_tree_Btag->Branch("Btag_R2", &temp_DataToTree[0]);
            temp_tree_Btag->Branch("Btag_thrustBm", &temp_DataToTree[1]);
            temp_tree_Btag->Branch("Btag_thrustOm", &temp_DataToTree[2]);
            temp_tree_Btag->Branch("Btag_cosTBTO", &temp_DataToTree[3]);
            temp_tree_Btag->Branch("Btag_cosTBz", &temp_DataToTree[4]);
            temp_tree_Btag->Branch("Btag_KSFWVariables_et", &temp_DataToTree[5]);
            temp_tree_Btag->Branch("Btag_KSFWVariables_mm2", &temp_DataToTree[6]);
            temp_tree_Btag->Branch("Btag_KSFWVariables_hso00", &temp_DataToTree[7]);
            temp_tree_Btag->Branch("Btag_KSFWVariables_hso01", &temp_DataToTree[8]);
            temp_tree_Btag->Branch("Btag_KSFWVariables_hso02", &temp_DataToTree[9]);
            temp_tree_Btag->Branch("Btag_KSFWVariables_hso03", &temp_DataToTree[10]);
            temp_tree_Btag->Branch("Btag_KSFWVariables_hso04", &temp_DataToTree[11]);
            temp_tree_Btag->Branch("Btag_KSFWVariables_hso10", &temp_DataToTree[12]);
            temp_tree_Btag->Branch("Btag_KSFWVariables_hso12", &temp_DataToTree[13]);
            temp_tree_Btag->Branch("Btag_KSFWVariables_hso14", &temp_DataToTree[14]);
            temp_tree_Btag->Branch("Btag_KSFWVariables_hso20", &temp_DataToTree[15]);
            temp_tree_Btag->Branch("Btag_KSFWVariables_hso22", &temp_DataToTree[16]);
            temp_tree_Btag->Branch("Btag_KSFWVariables_hso24", &temp_DataToTree[17]);
            temp_tree_Btag->Branch("Btag_KSFWVariables_hoo0", &temp_DataToTree[18]);
            temp_tree_Btag->Branch("Btag_KSFWVariables_hoo1", &temp_DataToTree[19]);
            temp_tree_Btag->Branch("Btag_KSFWVariables_hoo2", &temp_DataToTree[20]);
            temp_tree_Btag->Branch("Btag_KSFWVariables_hoo3", &temp_DataToTree[21]);
            temp_tree_Btag->Branch("Btag_KSFWVariables_hoo4", &temp_DataToTree[22]);
            temp_tree_Btag->Branch("Btag_CleoConeCS_1", &temp_DataToTree[23]);
            temp_tree_Btag->Branch("Btag_CleoConeCS_2", &temp_DataToTree[24]);
            temp_tree_Btag->Branch("Btag_CleoConeCS_3", &temp_DataToTree[25]);
            temp_tree_Btag->Branch("Btag_CleoConeCS_4", &temp_DataToTree[26]);
            temp_tree_Btag->Branch("Btag_CleoConeCS_5", &temp_DataToTree[27]);
            temp_tree_Btag->Branch("Btag_CleoConeCS_6", &temp_DataToTree[28]);
            temp_tree_Btag->Branch("Btag_CleoConeCS_7", &temp_DataToTree[29]);
            temp_tree_Btag->Branch("Btag_CleoConeCS_8", &temp_DataToTree[30]);
            temp_tree_Btag->Branch("Btag_CleoConeCS_9", &temp_DataToTree[31]);
            temp_tree_upsilon->Branch("missingMass2OfEvent", &temp_DataToTree[32]);
            temp_tree_upsilon->Branch("visibleEnergyOfEventCMS", &temp_DataToTree[33]);
            temp_tree_Btag->Branch("Btag_useCMSFrame_theta", &temp_DataToTree[34]);
            temp_tree_upsilon->Branch("extraInfo__boDecayHash__bc", &temp_DataToTree[35]);
            temp_tree_upsilon->Branch("extraInfo__boDecayHashExtended__bc", &temp_DataToTree[36]);

            if (DoesItHaveXsBranch) {
                // decay mode (MC level)
                temp_tree_Xs->Branch("nParticlesInList__boB__pl__clKcharge_total__bc", &temp_DecayDataToTree[0]);
                temp_tree_Xs->Branch("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &temp_DecayDataToTree[1]);
                temp_tree_Xs->Branch("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &temp_DecayDataToTree[2]);
                temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCcomb__bc", &temp_DecayDataToTree[3]);
                temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch1__bc", &temp_DecayDataToTree[4]);
                temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch2__bc", &temp_DecayDataToTree[5]);
                temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch3__bc", &temp_DecayDataToTree[6]);
                temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch4__bc", &temp_DecayDataToTree[7]);
                temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch5__bc", &temp_DecayDataToTree[8]);
                temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch6__bc", &temp_DecayDataToTree[9]);
                temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch7__bc", &temp_DecayDataToTree[10]);
                temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch8__bc", &temp_DecayDataToTree[11]);
                temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch9__bc", &temp_DecayDataToTree[12]);
                temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch10__bc", &temp_DecayDataToTree[13]);
                temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch11__bc", &temp_DecayDataToTree[14]);
                temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch12__bc", &temp_DecayDataToTree[15]);
                temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch13__bc", &temp_DecayDataToTree[16]);
                temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch14__bc", &temp_DecayDataToTree[17]);
                temp_tree_Xs->Branch("nParticlesInList__boXsu__clMCch15__bc", &temp_DecayDataToTree[18]);
                temp_tree_Xs->Branch("nParticlesInList__boB0__clKneutral_total__bc", &temp_DecayDataToTree[19]);
                temp_tree_Xs->Branch("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &temp_DecayDataToTree[20]);
                temp_tree_Xs->Branch("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &temp_DecayDataToTree[21]);
                temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCcomb__bc", &temp_DecayDataToTree[22]);
                temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch16__bc", &temp_DecayDataToTree[23]);
                temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch17__bc", &temp_DecayDataToTree[24]);
                temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch18__bc", &temp_DecayDataToTree[25]);
                temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch19__bc", &temp_DecayDataToTree[26]);
                temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch20__bc", &temp_DecayDataToTree[27]);
                temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch21__bc", &temp_DecayDataToTree[28]);
                temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch22__bc", &temp_DecayDataToTree[29]);
                temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch23__bc", &temp_DecayDataToTree[30]);
                temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch24__bc", &temp_DecayDataToTree[31]);
                temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch25__bc", &temp_DecayDataToTree[32]);
                temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch26__bc", &temp_DecayDataToTree[33]);
                temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch27__bc", &temp_DecayDataToTree[34]);
                temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch28__bc", &temp_DecayDataToTree[35]);
                temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch29__bc", &temp_DecayDataToTree[36]);
                temp_tree_Xs->Branch("nParticlesInList__boXsd__clMCch30__bc", &temp_DecayDataToTree[37]);
                temp_tree_Xs->Branch("nParticlesInList__bonu_e__clMC_signal__bc", &temp_DecayNparticlesDataToTree[0]);
                temp_tree_Xs->Branch("nParticlesInList__boB__pl__clMC_signal_total_e__bc", &temp_DecayNparticlesDataToTree[1]);
                temp_tree_Xs->Branch("nParticlesInList__boB0__clMC_signal_total_e__bc", &temp_DecayNparticlesDataToTree[2]);
                temp_tree_Xs->Branch("nParticlesInList__boB__pl__clPrimaryMC__bc", &temp_DecayNparticlesDataToTree[3]);
                temp_tree_Xs->Branch("nParticlesInList__boB0__clPrimaryMC__bc", &temp_DecayNparticlesDataToTree[4]);
                temp_tree_Xs->Branch("invMassInLists__bonu_e__clMC_signal__bc", &temp_DecaySystFFDataToTree[0]);
                temp_tree_Xs->Branch("averageValueInList__boB__pl__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", &temp_DecaySystFFDataToTree[1]);
                temp_tree_Xs->Branch("averageValueInList__boB0__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", &temp_DecaySystFFDataToTree[2]);
                temp_tree_Xs->Branch("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &temp_DecaySystFFDataToTree[3]);
                temp_tree_Xs->Branch("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &temp_DecaySystFFDataToTree[4]);
                temp_tree_Xs->Branch("averageValueInList__boB__pl__clMC_signal_total_e__cm__spM__bc", &temp_DecaySystFFDataToTree[5]);
                temp_tree_Xs->Branch("averageValueInList__boB0__clMC_signal_total_e__cm__spM__bc", &temp_DecaySystFFDataToTree[6]);
            }

            temp_tree_upsilon->Branch("MVA_BB", &temp_BB_output);
            temp_tree_upsilon->Branch("MVA_Continuum", &temp_Continuum_output);

            if (DoesItHaveJpsiOutput) {
                temp_tree_upsilon->Branch("nParticlesInList__bomu__pl__clfromUpsilonmychargedMuon__bc", &temp_nROE_mu);
                temp_tree_upsilon->Branch("nParticlesInList__boJ__slpsi__clfromUpsilontemp__bc", &temp_nROE_Jpsi);
                temp_tree_upsilon->Branch("nParticlesInList__boUpsilon__bo4S__bc__cltemp_withoutneutrino_BCS__bc", &temp_nROE_Upsilon_BCS);
                temp_tree_upsilon->Branch("nParticlesInList__boUpsilon__bo4S__bc__cltemp_withoutneutrino__bc", &temp_nROE_Upsilon);
            }

            /*================================================================*/

            for (unsigned int j = 0; j < temp_tree->GetEntries(); j++) { // Fill
                temp_tree->GetEntry(j);

                if (temp_Continuum_output > OContinuum && temp_BB_output > OBB) {
                    temp_tree_upsilon->Fill();
                    temp_tree_Bsig->Fill();
                    temp_tree_Btag->Fill();
                    if (DoesItHaveXsBranch) temp_tree_Xs->Fill();
                }

            }

            input_file->Close();

            temp_file->cd();
            temp_tree_upsilon->Write();
            temp_tree_Bsig->Write();
            temp_tree_Btag->Write();
            if (DoesItHaveXsBranch) temp_tree_Xs->Write();
            temp_file->Close();
        }

}
