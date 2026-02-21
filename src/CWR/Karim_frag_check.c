// last update: 2021-10-13
// for Belle2 data

typedef struct {
    double N_Xsu_ch4;
    double N_Xsu_ch5;
    double N_Xsu_ch6;
    double N_Xsu_ch7;
    double N_Xsu_ch11;
    double N_Xsd_ch19;
    double N_Xsd_ch20;
    double N_Xsd_ch21;
    double N_Xsd_ch22;
    double N_Xsd_ch26;

    double N_Xsu_Kpipi_ch[8];
    double N_Xsu_Kpipipi_ch[25];
    double N_Xsd_Kpipi_ch[8];
    double N_Xsd_Kpipipi_ch[23];
} DecayInfo;

void load_files(const char *dirname, std::vector<string>* names){
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

void LetsFillMXs(const char* dirname, DecayInfo* decayinfo) {

    int Decay[74] = { 0 };

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_decay = (TTree*)input_file->Get("decay");

        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clch4__bc", &Decay[0]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clch5__bc", &Decay[1]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clch6__bc", &Decay[2]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clch7__bc", &Decay[3]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clch11__bc", &Decay[4]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clch19__bc", &Decay[5]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clch20__bc", &Decay[6]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clch21__bc", &Decay[7]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clch22__bc", &Decay[8]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clch26__bc", &Decay[9]);

        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipi_ch1__bc", &Decay[10]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipi_ch2__bc", &Decay[11]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipi_ch3__bc", &Decay[12]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipi_ch4__bc", &Decay[13]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipi_ch5__bc", &Decay[14]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipi_ch6__bc", &Decay[15]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipi_ch7__bc", &Decay[16]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipi_ch8__bc", &Decay[17]);

        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch1__bc", &Decay[18]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch2__bc", &Decay[19]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch3__bc", &Decay[20]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch4__bc", &Decay[21]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch5__bc", &Decay[22]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch6__bc", &Decay[23]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch7__bc", &Decay[24]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch8__bc", &Decay[25]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch9__bc", &Decay[26]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch10__bc", &Decay[27]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch11__bc", &Decay[28]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch12__bc", &Decay[29]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch13__bc", &Decay[30]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch14__bc", &Decay[31]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch15__bc", &Decay[32]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch16__bc", &Decay[33]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch17__bc", &Decay[34]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch18__bc", &Decay[35]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch19__bc", &Decay[36]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch20__bc", &Decay[37]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch21__bc", &Decay[38]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch22__bc", &Decay[39]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch23__bc", &Decay[40]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch24__bc", &Decay[41]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsu__clKpipipi_ch25__bc", &Decay[42]);

        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipi_ch1__bc", &Decay[43]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipi_ch2__bc", &Decay[44]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipi_ch3__bc", &Decay[45]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipi_ch4__bc", &Decay[46]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipi_ch5__bc", &Decay[47]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipi_ch6__bc", &Decay[48]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipi_ch7__bc", &Decay[49]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipi_ch8__bc", &Decay[50]);

        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch1__bc", &Decay[51]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch2__bc", &Decay[52]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch3__bc", &Decay[53]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch4__bc", &Decay[54]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch5__bc", &Decay[55]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch6__bc", &Decay[56]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch7__bc", &Decay[57]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch8__bc", &Decay[58]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch9__bc", &Decay[59]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch10__bc", &Decay[60]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch11__bc", &Decay[61]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch12__bc", &Decay[62]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch13__bc", &Decay[63]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch14__bc", &Decay[64]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch15__bc", &Decay[65]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch16__bc", &Decay[66]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch17__bc", &Decay[67]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch18__bc", &Decay[68]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch19__bc", &Decay[69]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch20__bc", &Decay[70]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch21__bc", &Decay[71]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch22__bc", &Decay[72]);
        tree_decay->SetBranchAddress("nParticlesInList__boXsd__clKpipipi_ch23__bc", &Decay[73]);


        printf("%lld entries...\n", tree_decay->GetEntries());
        for (unsigned int j = 0; j < tree_decay->GetEntries(); j++) { // Fill
            tree_decay->GetEntry(j);

            decayinfo->N_Xsu_ch4 = decayinfo->N_Xsu_ch4 + Decay[0];
            decayinfo->N_Xsu_ch5 = decayinfo->N_Xsu_ch5 + Decay[1];
            decayinfo->N_Xsu_ch6 = decayinfo->N_Xsu_ch6 + Decay[2];
            decayinfo->N_Xsu_ch7 = decayinfo->N_Xsu_ch7 + Decay[3];
            decayinfo->N_Xsu_ch11 = decayinfo->N_Xsu_ch11 + Decay[4];
            decayinfo->N_Xsd_ch19 = decayinfo->N_Xsd_ch19 + Decay[5];
            decayinfo->N_Xsd_ch20 = decayinfo->N_Xsd_ch20 + Decay[6];
            decayinfo->N_Xsd_ch21 = decayinfo->N_Xsd_ch21 + Decay[7];
            decayinfo->N_Xsd_ch22 = decayinfo->N_Xsd_ch22 + Decay[8];
            decayinfo->N_Xsd_ch26 = decayinfo->N_Xsd_ch26 + Decay[9];

            if (Decay[0] > 0.5) {
                for (int i = 0; i < 4; i++) decayinfo->N_Xsu_Kpipi_ch[i] = decayinfo->N_Xsu_Kpipi_ch[i] + Decay[i + 10];
            }

            if (Decay[1] > 0.5) {
                for (int i = 4; i < 8; i++) decayinfo->N_Xsu_Kpipi_ch[i] = decayinfo->N_Xsu_Kpipi_ch[i] + Decay[i + 10];
            }

            if (Decay[2] > 0.5) {
                for (int i = 0; i < 13; i++) decayinfo->N_Xsu_Kpipipi_ch[i] = decayinfo->N_Xsu_Kpipipi_ch[i] + Decay[i + 18];
            }

            if (Decay[3] > 0.5) {
                for (int i = 13; i < 20; i++) decayinfo->N_Xsu_Kpipipi_ch[i] = decayinfo->N_Xsu_Kpipipi_ch[i] + Decay[i + 18];
            }

            if (Decay[4] > 0.5) {
                for (int i = 20; i < 25; i++) decayinfo->N_Xsu_Kpipipi_ch[i] = decayinfo->N_Xsu_Kpipipi_ch[i] + Decay[i + 18];
            }

            if (Decay[5] > 0.5) {
                for (int i = 0; i < 4; i++) decayinfo->N_Xsd_Kpipi_ch[i] = decayinfo->N_Xsd_Kpipi_ch[i] + Decay[i + 43];
            }

            if (Decay[6] > 0.5) {
                for (int i = 4; i < 8; i++) decayinfo->N_Xsd_Kpipi_ch[i] = decayinfo->N_Xsd_Kpipi_ch[i] + Decay[i + 43];
            }

            if (Decay[7] > 0.5) {
                for (int i = 0; i < 6; i++) decayinfo->N_Xsd_Kpipipi_ch[i] = decayinfo->N_Xsd_Kpipipi_ch[i] + Decay[i + 51];
            }

            if (Decay[8] > 0.5) {
                for (int i = 6; i < 18; i++) decayinfo->N_Xsd_Kpipipi_ch[i] = decayinfo->N_Xsd_Kpipipi_ch[i] + Decay[i + 51];
            }

            if (Decay[9] > 0.5) {
                for (int i = 18; i < 23; i++) decayinfo->N_Xsd_Kpipipi_ch[i] = decayinfo->N_Xsd_Kpipipi_ch[i] + Decay[i + 51];
            }

        }
        input_file->Close();

    }

}

void Karim_frag_check(){

    const char* Xsu_dirname = "./Xsu";
    const char* Xsd_dirname = "./Xsd";

    DecayInfo decayinfo = { 0.0 };

    LetsFillMXs(Xsu_dirname, &decayinfo);
    LetsFillMXs(Xsd_dirname, &decayinfo);

    printf("======================   B+ -> K+ pi- pi+   ======================\n");
    printf("%lf\n", decayinfo.N_Xsu_Kpipi_ch[0] / decayinfo.N_Xsu_ch4);
    printf("%lf\n", decayinfo.N_Xsu_Kpipi_ch[1] / decayinfo.N_Xsu_ch4);
    printf("%lf\n", decayinfo.N_Xsu_Kpipi_ch[2] / decayinfo.N_Xsu_ch4);
    printf("%lf\n", decayinfo.N_Xsu_Kpipi_ch[3] / decayinfo.N_Xsu_ch4);

    printf("======================   B+ -> K0 pi+ pi0   ======================\n");
    printf("%lf\n", decayinfo.N_Xsu_Kpipi_ch[4] / decayinfo.N_Xsu_ch5);
    printf("%lf\n", decayinfo.N_Xsu_Kpipi_ch[5] / decayinfo.N_Xsu_ch5);
    printf("%lf\n", decayinfo.N_Xsu_Kpipi_ch[6] / decayinfo.N_Xsu_ch5);
    printf("%lf\n", decayinfo.N_Xsu_Kpipi_ch[7] / decayinfo.N_Xsu_ch5);

    printf("======================   B+ -> K+ pi- pi+ pi0   ======================\n");
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[0] / decayinfo.N_Xsu_ch6);
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[1] / decayinfo.N_Xsu_ch6);
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[2] / decayinfo.N_Xsu_ch6);
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[3] / decayinfo.N_Xsu_ch6);
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[4] / decayinfo.N_Xsu_ch6);
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[5] / decayinfo.N_Xsu_ch6);
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[6] / decayinfo.N_Xsu_ch6);
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[7] / decayinfo.N_Xsu_ch6);
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[8] / decayinfo.N_Xsu_ch6);
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[9] / decayinfo.N_Xsu_ch6);
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[10] / decayinfo.N_Xsu_ch6);
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[11] / decayinfo.N_Xsu_ch6);
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[12] / decayinfo.N_Xsu_ch6);

    printf("======================   B+ -> K0 pi+ pi- pi+   ======================\n");
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[13] / decayinfo.N_Xsu_ch7);
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[14] / decayinfo.N_Xsu_ch7);
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[15] / decayinfo.N_Xsu_ch7);
    // printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[16] / decayinfo.N_Xsu_ch7);
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[17] / decayinfo.N_Xsu_ch7);
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[18] / decayinfo.N_Xsu_ch7);
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[19] / decayinfo.N_Xsu_ch7);

    printf("======================   B+ -> K0 pi+ pi0 pi0   ======================\n");
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[20] / decayinfo.N_Xsu_ch11);
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[21] / decayinfo.N_Xsu_ch11);
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[22] / decayinfo.N_Xsu_ch11);
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[23] / decayinfo.N_Xsu_ch11);
    printf("%lf\n", decayinfo.N_Xsu_Kpipipi_ch[24] / decayinfo.N_Xsu_ch11);

    printf("======================   B0 -> K+ pi- pi0   ======================\n");
    printf("%lf\n", decayinfo.N_Xsd_Kpipi_ch[0] / decayinfo.N_Xsd_ch19);
    printf("%lf\n", decayinfo.N_Xsd_Kpipi_ch[1] / decayinfo.N_Xsd_ch19);
    printf("%lf\n", decayinfo.N_Xsd_Kpipi_ch[2] / decayinfo.N_Xsd_ch19);
    printf("%lf\n", decayinfo.N_Xsd_Kpipi_ch[3] / decayinfo.N_Xsd_ch19);

    printf("======================   B0 -> K0 pi+ pi-   ======================\n");
    printf("%lf\n", decayinfo.N_Xsd_Kpipi_ch[4] / decayinfo.N_Xsd_ch20);
    printf("%lf\n", decayinfo.N_Xsd_Kpipi_ch[5] / decayinfo.N_Xsd_ch20);
    printf("%lf\n", decayinfo.N_Xsd_Kpipi_ch[6] / decayinfo.N_Xsd_ch20);
    printf("%lf\n", decayinfo.N_Xsd_Kpipi_ch[7] / decayinfo.N_Xsd_ch20);

    printf("======================   B0 -> K+ pi- pi+ pi-   ======================\n");
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[0] / decayinfo.N_Xsd_ch21);
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[1] / decayinfo.N_Xsd_ch21);
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[2] / decayinfo.N_Xsd_ch21);
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[3] / decayinfo.N_Xsd_ch21);
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[4] / decayinfo.N_Xsd_ch21);
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[5] / decayinfo.N_Xsd_ch21);

    printf("======================   B0 -> K+ pi- pi+ pi-   ======================\n");
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[6] / decayinfo.N_Xsd_ch22);
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[7] / decayinfo.N_Xsd_ch22);
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[8] / decayinfo.N_Xsd_ch22);
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[9] / decayinfo.N_Xsd_ch22);
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[10] / decayinfo.N_Xsd_ch22);
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[11] / decayinfo.N_Xsd_ch22);
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[12] / decayinfo.N_Xsd_ch22);
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[13] / decayinfo.N_Xsd_ch22);
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[14] / decayinfo.N_Xsd_ch22);
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[15] / decayinfo.N_Xsd_ch22);
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[16] / decayinfo.N_Xsd_ch22);
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[17] / decayinfo.N_Xsd_ch22);

    printf("======================   B0 -> K+ pi- pi+ pi-   ======================\n");
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[18] / decayinfo.N_Xsd_ch26);
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[19] / decayinfo.N_Xsd_ch26);
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[20] / decayinfo.N_Xsd_ch26);
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[21] / decayinfo.N_Xsd_ch26);
    printf("%lf\n", decayinfo.N_Xsd_Kpipipi_ch[22] / decayinfo.N_Xsd_ch26);

}
