// last update: 2022-01-13-00
// for Belle2 data

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
# define N_event_info 15
# define N_Upsilon_info 14
# define N_Bsig_info 28
# define N_Btag_info 7
# define N_decay 38 // five decay mode + others

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

void Deconvertor(const char* dirname, double OBB, double OContinuum){

        std::vector<std::string> names;
        load_files(dirname, &names);

        for (unsigned int i = 0; i < names.size(); i++) {
            const char* filename = (dirname + std::string("/") + names.at(i)).c_str();
            std::string string_filename(filename);
            std::string OnlyFileName = string_filename.substr(string_filename.find_last_of("\\/") + 1, string_filename.size() - string_filename.find_last_of("\\/"));
            size_t lastindex = OnlyFileName.find_last_of(".");
            std::string rawname = OnlyFileName.substr(0, lastindex);

            TFile* input_file = new TFile(filename, "read");

            TTree* temp_tree = (TTree*)input_file->Get("data");

            int temp_EventDataToTree[N_event_info / 3];
            double temp_UpsilonDataToTree[N_Upsilon_info];
            double temp_BsigDataToTree[N_Bsig_info];
            double temp_BtagDataToTree[N_Btag_info];
            double temp_DataToTree[N_Needed_info];
            double temp_DecayDataToTree[N_decay];
            double temp_Upsilon_decayIDToTree;
            double temp_Bsig_decayIDToTree;
            int temp_flag;
            float temp_BB_output;
            float temp_Continuum_output;

            bool DoesItHaveXsBranch = false;
            for (int i = 0; i < temp_tree->GetListOfBranches()->LastIndex(); i++) {
                if (temp_tree->GetListOfBranches()->At(i)->GetName() == std::string("nParticlesInList__boB__pl__clKcharge_total__bc")) DoesItHaveXsBranch = true;
            }

            // Get data
            /*================================================================*/
            // get event_info
            temp_tree->SetBranchAddress("__experiment__", &temp_EventDataToTree[0]);
            temp_tree->SetBranchAddress("__run__", &temp_EventDataToTree[1]);
            temp_tree->SetBranchAddress("__event__", &temp_EventDataToTree[2]);
            temp_tree->SetBranchAddress("__candidate__", &temp_EventDataToTree[3]);
            temp_tree->SetBranchAddress("__ncandidates__", &temp_EventDataToTree[4]);

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

            // get Bsig_info
            temp_tree->SetBranchAddress("Bsig_E", &temp_BsigDataToTree[0]);
            temp_tree->SetBranchAddress("Bsig_useCMSFrame_E", &temp_BsigDataToTree[1]);
            temp_tree->SetBranchAddress("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc", &temp_BsigDataToTree[2]);
            temp_tree->SetBranchAddress("Bsig_p", &temp_BsigDataToTree[3]);
            temp_tree->SetBranchAddress("Bsig_useCMSFrame_p", &temp_BsigDataToTree[4]);
            temp_tree->SetBranchAddress("useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc", &temp_BsigDataToTree[5]);
            temp_tree->SetBranchAddress("Bsig_M", &temp_BsigDataToTree[6]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcvetomass", &temp_BsigDataToTree[7]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_DcvetodmID", &temp_BsigDataToTree[8]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_Dcvetoabsdm", &temp_BsigDataToTree[9]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_Dnvetomass", &temp_BsigDataToTree[10]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_DnvetodmID", &temp_BsigDataToTree[11]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_Dnvetoabsdm", &temp_BsigDataToTree[12]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_3D_distance", &temp_BsigDataToTree[13]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_correction", &temp_BsigDataToTree[14]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_rel_uncer_up", &temp_BsigDataToTree[15]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_rel_uncer_dn", &temp_BsigDataToTree[16]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_Pion_PID_correction", &temp_BsigDataToTree[17]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_Pion_PID_rel_uncer_up", &temp_BsigDataToTree[18]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_Pion_PID_rel_uncer_dn", &temp_BsigDataToTree[19]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_nKslow1", &temp_BsigDataToTree[20]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_nKslow2", &temp_BsigDataToTree[21]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_nKslow3", &temp_BsigDataToTree[22]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_nPislow1", &temp_BsigDataToTree[23]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_nPislow2", &temp_BsigDataToTree[24]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_nPislow3", &temp_BsigDataToTree[25]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_nKexcep", &temp_BsigDataToTree[26]);
            temp_tree->SetBranchAddress("Bsig_daughter_0_extraInfo_nPiexcep", &temp_BsigDataToTree[27]);

            // get Btag_info
            temp_tree->SetBranchAddress("Btag_extraInfo_decayModeID", &temp_BtagDataToTree[0]);
            temp_tree->SetBranchAddress("Btag_Mbc", &temp_BtagDataToTree[1]);
            temp_tree->SetBranchAddress("Btag_deltaE", &temp_BtagDataToTree[2]);
            temp_tree->SetBranchAddress("Btag_E", &temp_BtagDataToTree[3]);
            temp_tree->SetBranchAddress("Btag_useCMSFrame_E", &temp_BtagDataToTree[4]);
            temp_tree->SetBranchAddress("Btag_extraInfo_SignalProbability", &temp_BtagDataToTree[5]);
            temp_tree->SetBranchAddress("Btag_chiProb", &temp_BtagDataToTree[6]);

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
            }

            // flag
            temp_tree->SetBranchAddress("flag", &temp_flag);

            // MLP
            temp_tree->SetBranchAddress("TMVA_BB", &temp_BB_output);
            temp_tree->SetBranchAddress("TMVA_Continuum", &temp_Continuum_output);
            /*================================================================*/


            TFile* temp_file = new TFile((rawname + "_after_TMVA.root").c_str(), "recreate");
            temp_file->cd();
            TTree* temp_tree_upsilon = new TTree("Upsilon", "");
            TTree* temp_tree_Bsig = new TTree("Bsig", "");
            TTree* temp_tree_Btag = new TTree("Btag", "");
            TTree* temp_tree_Xs = nullptr;
            if (DoesItHaveXsBranch) temp_tree_Xs = new TTree("Xs", "");

            // print root file
                /*================================================================*/
            // get event_info
            temp_tree_upsilon->Branch("__experiment__", &temp_EventDataToTree[0]);
            temp_tree_upsilon->Branch("__run__", &temp_EventDataToTree[1]);
            temp_tree_upsilon->Branch("__event__", &temp_EventDataToTree[2]);
            temp_tree_upsilon->Branch("__candidate__", &temp_EventDataToTree[3]);
            temp_tree_upsilon->Branch("__ncandidates__", &temp_EventDataToTree[4]);
            temp_tree_Bsig->Branch("__experiment__", &temp_EventDataToTree[0]);
            temp_tree_Bsig->Branch("__run__", &temp_EventDataToTree[1]);
            temp_tree_Bsig->Branch("__event__", &temp_EventDataToTree[2]);
            temp_tree_Bsig->Branch("__candidate__", &temp_EventDataToTree[3]);
            temp_tree_Bsig->Branch("__ncandidates__", &temp_EventDataToTree[4]);
            temp_tree_Btag->Branch("__experiment__", &temp_EventDataToTree[0]);
            temp_tree_Btag->Branch("__run__", &temp_EventDataToTree[1]);
            temp_tree_Btag->Branch("__event__", &temp_EventDataToTree[2]);
            temp_tree_Btag->Branch("__candidate__", &temp_EventDataToTree[3]);
            temp_tree_Btag->Branch("__ncandidates__", &temp_EventDataToTree[4]);

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

            // get Bsig_info
            temp_tree_Bsig->Branch("Bsig_E", &temp_BsigDataToTree[0]);
            temp_tree_Bsig->Branch("Bsig_useCMSFrame_E", &temp_BsigDataToTree[1]);
            temp_tree_upsilon->Branch("useTagSideRecoilRestFrame__bodaughter__bo1__cmE__bc__cm0__bc", &temp_BsigDataToTree[2]);
            temp_tree_Bsig->Branch("Bsig_p", &temp_BsigDataToTree[3]);
            temp_tree_Bsig->Branch("Bsig_useCMSFrame_p", &temp_BsigDataToTree[4]);
            temp_tree_upsilon->Branch("useTagSideRecoilRestFrame__bodaughter__bo1__cmp__bc__cm0__bc", &temp_BsigDataToTree[5]);
            temp_tree_Bsig->Branch("Bsig_M", &temp_BsigDataToTree[6]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Dcvetomass", &temp_BsigDataToTree[7]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_DcvetodmID", &temp_BsigDataToTree[8]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Dcvetoabsdm", &temp_BsigDataToTree[9]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Dnvetomass", &temp_BsigDataToTree[10]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_DnvetodmID", &temp_BsigDataToTree[11]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Dnvetoabsdm", &temp_BsigDataToTree[12]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_KS0_3D_distance", &temp_BsigDataToTree[13]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Kaon_PID_correction", &temp_BsigDataToTree[14]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Kaon_PID_rel_uncer_up", &temp_BsigDataToTree[15]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Kaon_PID_rel_uncer_dn", &temp_BsigDataToTree[16]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Pion_PID_correction", &temp_BsigDataToTree[17]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Pion_PID_rel_uncer_up", &temp_BsigDataToTree[18]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_Pion_PID_rel_uncer_dn", &temp_BsigDataToTree[19]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nKslow1", &temp_BsigDataToTree[20]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nKslow2", &temp_BsigDataToTree[21]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nKslow3", &temp_BsigDataToTree[22]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nPislow1", &temp_BsigDataToTree[23]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nPislow2", &temp_BsigDataToTree[24]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nPislow3", &temp_BsigDataToTree[25]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nKexcep", &temp_BsigDataToTree[26]);
            temp_tree_Bsig->Branch("Bsig_daughter_0_extraInfo_nPiexcep", &temp_BsigDataToTree[27]);

            // get Btag_info
            temp_tree_Btag->Branch("Btag_extraInfo_decayModeID", &temp_BtagDataToTree[0]);
            temp_tree_Btag->Branch("Btag_Mbc", &temp_BtagDataToTree[1]);
            temp_tree_Btag->Branch("Btag_deltaE", &temp_BtagDataToTree[2]);
            temp_tree_Btag->Branch("Btag_E", &temp_BtagDataToTree[3]);
            temp_tree_Btag->Branch("Btag_useCMSFrame_E", &temp_BtagDataToTree[4]);
            temp_tree_Btag->Branch("Btag_extraInfo_SignalProbability", &temp_BtagDataToTree[5]);
            temp_tree_Btag->Branch("Btag_chiProb", &temp_BtagDataToTree[6]);

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
            }

            temp_tree_upsilon->Branch("TMVA_BB", &temp_BB_output);
            temp_tree_upsilon->Branch("TMVA_Continuum", &temp_Continuum_output);

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
