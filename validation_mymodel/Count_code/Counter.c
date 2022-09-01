// last update: 2022-07-16
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

# define N_Needed_info 37
# define N_event_info 15
# define N_Upsilon_info 55
# define N_Bsig_info 81
# define N_Btag_info 9
# define N_decay 48 // five decay mode + others + 10 variables for systematics

# define Nstep 20
# define start 0.8
# define end 1.0

// arXiv:1409.4557v2
# define TB0 1.5195 // (Table. 1)
# define TBp 1.6384 // (Table. 1)
# define BR_Kplus_nunubar 0.00000398 // (eq. 10)
# define BR_K0star_nunubar 0.00000919 // (eq. 11)
# define BR_K0_nunubar (BR_Kplus_nunubar*TB0/TBp) // under (eq. 15)
# define BR_Kplusstar_nunubar (BR_K0star_nunubar*TBp/TB0) // under (eq. 15)
# define BR_Xs_nunubar 0.000029 // (eq. 23)
# define BR_Xsu_nonresonant_nunubar (BR_Xs_nunubar - BR_Kplus_nunubar - BR_Kplusstar_nunubar)
# define BR_Xsd_nonresonant_nunubar (BR_Xs_nunubar - BR_K0_nunubar - BR_K0star_nunubar)

// https://confluence.desy.de/pages/viewpage.action?pageId=107054222
# define N_BpBp_1invab 565400000.0
# define N_B0B0_1invab 534600000.0

# define N_Kplus_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Kplus_nunubar)
# define N_Kplusstar_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Kplusstar_nunubar)
# define N_Xsu_nonresonant_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Xsu_nonresonant_nunubar)
# define N_K0_nunubar_1invab (2.0 * N_B0B0_1invab * BR_K0_nunubar)
# define N_K0star_nunubar_1invab (2.0 * N_B0B0_1invab * BR_K0star_nunubar)
# define N_Xsd_nunubar_1invab (2.0 * N_B0B0_1invab * BR_Xsd_nonresonant_nunubar)

// my MC sample number
# define N_Kplus_nunubar 10000000.0
# define N_K0_nunubar 10000000.0
# define N_Kplusstar_nunubar 10000000.0
# define N_K0star_nunubar 10000000.0
# define N_Xsu_nonresonant_nunubar 50000000.0
# define N_Xsd_nonresonant_nunubar 50000000.0

// scale factor for each MC sample
# define Scale_Kplus (N_Kplus_nunubar_1invab/N_Kplus_nunubar)
# define Scale_Kplusstar (N_Kplusstar_nunubar_1invab/N_Kplusstar_nunubar)
# define Scale_Xsu_nonresonant (N_Xsu_nonresonant_nunubar_1invab/N_Xsu_nonresonant_nunubar)
# define Scale_K0 (N_K0_nunubar_1invab/N_K0_nunubar)
# define Scale_K0star (N_K0star_nunubar_1invab/N_K0star_nunubar)
# define Scale_Xsd_nonresonant (N_Xsd_nunubar_1invab/N_Xsd_nonresonant_nunubar)

void load_files(const char *dirname, std::vector<std::string>* names){
   TSystemDirectory dir(dirname, dirname);
   TList *files = dir.GetListOfFiles();
   if (files) {
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next())) {
         fname = file->GetName();
         if (!file->IsDirectory() && fname.EndsWith(".root")) {
            names->push_back(fname.Data());
         }
      }
   }
}

void load_files(const char* dirname, std::vector<string>* names, const char* included_string) {
    TSystemDirectory dir(dirname, dirname);
    TList* files = dir.GetListOfFiles();
    if (files) {
        TSystemFile* file;
        TString fname;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(".root") && fname.Contains(included_string)) {
                names->push_back(fname.Data());
            }
        }
    }
}

enum DecayMode { // reco level
    B2Kc = 0,
    B2KcPi0,
    B2Ks0Pic,
    B2KcPicPic,
    B2Ks0PicPi0,
    B2KcPicPicPi0,
    B2Ks0PicPicPic,
    B2KcPicPicPicPic,
    B2Ks0PicPicPicPi0,
    B2KcPi0Pi0,
    B2Ks0PicPi0Pi0,
    B2KcPicPicPi0Pi0,
    B2KcKcKc,
    B2KcKcKs0Pic,
    B2KcKcKcPi0,
    B02Ks0,
    B02KcPic,
    B02Ks0Pi0,
    B02KcPicPi0,
    B02Ks0PicPic,
    B02KcPicPicPic,
    B02Ks0PicPicPi0,
    B02KcPicPicPicPi0,
    B02Ks0PicPicPicPic,
    B02Ks0Pi0Pi0,
    B02KcPicPi0Pi0,
    B02Ks0PicPicPi0Pi0,
    B02KcKcKs0,
    B02KcKcKcPic,
    B02KcKcKs0Pi0,
    MAX_NUM_DECAYMODE
};
enum DecayModeMC { // MC level
    Xsu2Kc_MC = 0,
    Xsu2Kcstar2KcPi0_MC,
    Xsu2Kcstar2K0Pic_MC,
    Xsu2KcPi0_MC,
    Xsu2K0Pic_MC,
    Xsu2KcPicPic_MC,
    Xsu2K0PicPi0_MC,
    Xsu2KcPicPicPi0_MC,
    Xsu2K0PicPicPic_MC,
    Xsu2KcPicPicPicPic_MC,
    Xsu2K0PicPicPicPi0_MC,
    Xsu2KcPi0Pi0_MC,
    Xsu2K0PicPi0Pi0_MC,
    Xsu2KcPicPicPi0Pi0_MC,
    Xsu2KcKcKc_MC,
    Xsu2KcKcK0Pic_MC,
    Xsu2KcKcKcPi0_MC,
    Xsd2K0_MC,
    Xsd2K0star2KcPic_MC,
    Xsd2K0star2K0Pi0_MC,
    Xsd2KcPic_MC,
    Xsd2K0Pi0_MC,
    Xsd2KcPicPi0_MC,
    Xsd2K0PicPic_MC,
    Xsd2KcPicPicPic_MC,
    Xsd2K0PicPicPi0_MC,
    Xsd2KcPicPicPicPi0_MC,
    Xsd2K0PicPicPicPic_MC,
    Xsd2K0Pi0Pi0_MC,
    Xsd2KcPicPi0Pi0_MC,
    Xsd2K0PicPicPi0Pi0_MC,
    Xsd2KcKcK0_MC,
    Xsd2KcKcKcPic_MC,
    Xsd2KcKcK0Pi0_MC,
    other,
    MAX_NUM_DECAYMODE_MC
};

bool TrueIfDecayModeMatch_MC(double Decay[N_decay], DecayModeMC decaymodeMC) {

    switch (decaymodeMC) {
    case Xsu2Kc_MC:
        if (Decay[0] > 0) return true;
        return false;
        break;
    case Xsu2Kcstar2KcPi0_MC:
        if (Decay[2] > 0) return true;
        return false;
        break;
    case Xsu2Kcstar2K0Pic_MC:
        if (Decay[1] > 0) return true;
        return false;
        break;
    case Xsu2KcPi0_MC:
        if (Decay[5] > 0) return true;
        return false;
        break;
    case Xsu2K0Pic_MC:
        if (Decay[6] > 0) return true;
        return false;
        break;
    case Xsu2KcPicPic_MC:
        if (Decay[7] > 0) return true;
        return false;
        break;
    case Xsu2K0PicPi0_MC:
        if (Decay[8] > 0) return true;
        return false;
        break;
    case Xsu2KcPicPicPi0_MC:
        if (Decay[9] > 0) return true;
        return false;
        break;
    case Xsu2K0PicPicPic_MC:
        if (Decay[10] > 0) return true;
        return false;
        break;
    case Xsu2KcPicPicPicPic_MC:
        if (Decay[11] > 0) return true;
        return false;
        break;
    case Xsu2K0PicPicPicPi0_MC:
        if (Decay[12] > 0) return true;
        return false;
        break;
    case Xsu2KcPi0Pi0_MC:
        if (Decay[13] > 0) return true;
        return false;
        break;
    case Xsu2K0PicPi0Pi0_MC:
        if (Decay[14] > 0) return true;
        return false;
        break;
    case Xsu2KcPicPicPi0Pi0_MC:
        if (Decay[15] > 0) return true;
        return false;
        break;
    case Xsu2KcKcKc_MC:
        if (Decay[16] > 0) return true;
        return false;
        break;
    case Xsu2KcKcK0Pic_MC:
        if (Decay[17] > 0) return true;
        return false;
        break;
    case Xsu2KcKcKcPi0_MC:
        if (Decay[18] > 0) return true;
        return false;
        break;
    case Xsd2K0_MC:
        if (Decay[19] > 0) return true;
        return false;
        break;
    case Xsd2K0star2KcPic_MC:
        if (Decay[20] > 0) return true;
        return false;
        break;
    case Xsd2K0star2K0Pi0_MC:
        if (Decay[21] > 0) return true;
        return false;
        break;
    case Xsd2KcPic_MC:
        if (Decay[24] > 0) return true;
        return false;
        break;
    case Xsd2K0Pi0_MC:
        if (Decay[25] > 0) return true;
        return false;
        break;
    case Xsd2KcPicPi0_MC:
        if (Decay[26] > 0) return true;
        return false;
        break;
    case Xsd2K0PicPic_MC:
        if (Decay[27] > 0) return true;
        return false;
        break;
    case Xsd2KcPicPicPic_MC:
        if (Decay[28] > 0) return true;
        return false;
        break;
    case Xsd2K0PicPicPi0_MC:
        if (Decay[29] > 0) return true;
        return false;
        break;
    case Xsd2KcPicPicPicPi0_MC:
        if (Decay[30] > 0) return true;
        return false;
        break;
    case Xsd2K0PicPicPicPic_MC:
        if (Decay[31] > 0) return true;
        return false;
        break;
    case Xsd2K0Pi0Pi0_MC:
        if (Decay[32] > 0) return true;
        return false;
        break;
    case Xsd2KcPicPi0Pi0_MC:
        if (Decay[33] > 0) return true;
        return false;
        break;
    case Xsd2K0PicPicPi0Pi0_MC:
        if (Decay[34] > 0) return true;
        return false;
        break;
    case Xsd2KcKcK0_MC:
        if (Decay[35] > 0) return true;
        return false;
        break;
    case Xsd2KcKcKcPic_MC:
        if (Decay[36] > 0) return true;
        return false;
        break;
    case Xsd2KcKcK0Pi0_MC:
        if (Decay[37] > 0) return true;
        return false;
        break;
    case other:
        return true;
        break;
    default:
        printf("ERROR! Input value of TrueIfDecayModeMatch_MC is not appropriate\n");
        exit(1);
        break;
    }

    printf("ERROR! Input value of TrueIfDecayModeMatch_MC is not appropriate\n");
    exit(1);
    return false;
}

void LetsFill(const char* dirname, double Decay_num[MAX_NUM_DECAYMODE + 1], double weight = 1) {
    double var = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    double Decay[N_decay] = { 0.0 };

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_Xs = (TTree*)input_file->Get("Upsilon");

        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKcharge_total__bc", &Decay[0]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &Decay[1]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &Decay[2]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCcomb__bc", &Decay[3]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch1__bc", &Decay[4]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch2__bc", &Decay[5]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch3__bc", &Decay[6]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch4__bc", &Decay[7]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch5__bc", &Decay[8]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch6__bc", &Decay[9]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch7__bc", &Decay[10]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch8__bc", &Decay[11]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch9__bc", &Decay[12]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch10__bc", &Decay[13]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch11__bc", &Decay[14]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch12__bc", &Decay[15]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch13__bc", &Decay[16]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch14__bc", &Decay[17]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch15__bc", &Decay[18]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKneutral_total__bc", &Decay[19]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &Decay[20]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &Decay[21]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCcomb__bc", &Decay[22]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch16__bc", &Decay[23]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch17__bc", &Decay[24]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch18__bc", &Decay[25]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch19__bc", &Decay[26]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch20__bc", &Decay[27]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch21__bc", &Decay[28]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch22__bc", &Decay[29]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch23__bc", &Decay[30]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch24__bc", &Decay[31]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch25__bc", &Decay[32]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch26__bc", &Decay[33]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch27__bc", &Decay[34]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch28__bc", &Decay[35]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch29__bc", &Decay[36]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch30__bc", &Decay[37]);
        tree_Xs->SetBranchAddress("nParticlesInList__bonu_e__clMC_signal__bc", &Decay[38]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clMC_signal_total_e__bc", &Decay[39]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clMC_signal_total_e__bc", &Decay[40]);
        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &Decay[41]);
        tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", &Decay[42]);
        tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", &Decay[43]);
        tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Decay[44]);
        tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Decay[45]);
        tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spM__bc", &Decay[46]);
        tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spM__bc", &Decay[47]);

        printf("%lld entries...\n", tree_Xs->GetEntries());
        for (unsigned int j = 0; j < tree_Xs->GetEntries(); j++) { // Fill
            tree_Xs->GetEntry(j);

            int decaymodeid_MC = -1;

            for (int j = 0; j < MAX_NUM_DECAYMODE_MC; j++) { // find MC decay mode
                if (TrueIfDecayModeMatch_MC(Decay, static_cast<DecayModeMC>(j))) {
                    decaymodeid_MC = j;
                    break;
                }
            }
            if (decaymodeid_MC == MAX_NUM_DECAYMODE_MC) {
                printf("ERROR! MC decay id cannot be found\n");
                exit(1);
            }

            int decaymodeid_MC_for_square = -1;
            if (decaymodeid_MC == 0) decaymodeid_MC_for_square = 0;
            else if (decaymodeid_MC == 1)decaymodeid_MC_for_square = 1;
            else if (decaymodeid_MC == 2)decaymodeid_MC_for_square = 2;
            else if (decaymodeid_MC == 3)decaymodeid_MC_for_square = 1;
            else if (decaymodeid_MC == 4)decaymodeid_MC_for_square = 2;
            else if (decaymodeid_MC == 5)decaymodeid_MC_for_square = 3;
            else if (decaymodeid_MC == 6)decaymodeid_MC_for_square = 4;
            else if (decaymodeid_MC == 7)decaymodeid_MC_for_square = 5;
            else if (decaymodeid_MC == 8)decaymodeid_MC_for_square = 6;
            else if (decaymodeid_MC == 9)decaymodeid_MC_for_square = 7;
            else if (decaymodeid_MC == 10)decaymodeid_MC_for_square = 8;
            else if (decaymodeid_MC == 11)decaymodeid_MC_for_square = 9;
            else if (decaymodeid_MC == 12)decaymodeid_MC_for_square = 10;
            else if (decaymodeid_MC == 13)decaymodeid_MC_for_square = 11;
            else if (decaymodeid_MC == 14)decaymodeid_MC_for_square = 12;
            else if (decaymodeid_MC == 15)decaymodeid_MC_for_square = 13;
            else if (decaymodeid_MC == 16)decaymodeid_MC_for_square = 14;
            else if (decaymodeid_MC == 17)decaymodeid_MC_for_square = 15;
            else if (decaymodeid_MC == 18)decaymodeid_MC_for_square = 16;
            else if (decaymodeid_MC == 19)decaymodeid_MC_for_square = 17;
            else if (decaymodeid_MC == 20)decaymodeid_MC_for_square = 16;
            else if (decaymodeid_MC == 21)decaymodeid_MC_for_square = 17;
            else if (decaymodeid_MC == 22)decaymodeid_MC_for_square = 18;
            else if (decaymodeid_MC == 23)decaymodeid_MC_for_square = 19;
            else if (decaymodeid_MC == 24)decaymodeid_MC_for_square = 20;
            else if (decaymodeid_MC == 25)decaymodeid_MC_for_square = 21;
            else if (decaymodeid_MC == 26)decaymodeid_MC_for_square = 22;
            else if (decaymodeid_MC == 27)decaymodeid_MC_for_square = 23;
            else if (decaymodeid_MC == 28)decaymodeid_MC_for_square = 24;
            else if (decaymodeid_MC == 29)decaymodeid_MC_for_square = 25;
            else if (decaymodeid_MC == 30)decaymodeid_MC_for_square = 26;
            else if (decaymodeid_MC == 31)decaymodeid_MC_for_square = 27;
            else if (decaymodeid_MC == 32)decaymodeid_MC_for_square = 28;
            else if (decaymodeid_MC == 33)decaymodeid_MC_for_square = 29;
            else if (decaymodeid_MC == 34)decaymodeid_MC_for_square = 30;

            Decay_num[decaymodeid_MC_for_square] = Decay_num[decaymodeid_MC_for_square] + weight;

        }
        input_file->Close();

    }

}

int Counter(){

    std::vector<string> names;
    const char* Knunu_dirname = "/home/jwpark/storage/BKG_gbasf2/BKG_total/CHG_after_FEISKIM_light_0_fix/train";
    const char* Kstarnunu_dirname = "/home/jwpark/storage/BKG_gbasf2/BKG_total/CHG_after_FEISKIM_light_0_fix/train";
    const char* Xsununu_dirname = "/home/jwpark/storage/BKG_gbasf2/BKG_total/CHG_after_FEISKIM_light_0_fix/train";
    const char* K0nunu_dirname = "/home/jwpark/storage/BKG_gbasf2/BKG_total/CHG_after_FEISKIM_light_0_fix/train";
    const char* K0starnunu_dirname = "/home/jwpark/storage/BKG_gbasf2/BKG_total/CHG_after_FEISKIM_light_0_fix/train";
    const char* Xsdnunu_dirname = "/home/jwpark/storage/BKG_gbasf2/BKG_total/CHG_after_FEISKIM_light_0_fix/train";

    double Decay_num[MAX_NUM_DECAYMODE + 1] = { 0.0 };

    LetsFill(Knunu_dirname, Decay_num, Scale_Kplus);
    LetsFill(Kstarnunu_dirname, Decay_num, Scale_Kplusstar);
    LetsFill(Xsununu_dirname, Decay_num, Scale_Xsu_nonresonant);
    LetsFill(K0nunu_dirname, Decay_num, Scale_K0);
    LetsFill(K0starnunu_dirname, Decay_num, Scale_K0star);
    LetsFill(Xsdnunu_dirname, Decay_num, Scale_Xsd_nonresonant);

    for (int i = 0; i < MAX_NUM_DECAYMODE + 1;i++) {
        printf("Num of event with MC decay %d: %lf\n", Decay_num[i]);
    }

    return 0;
}
