// last update: 2023-07-14
// for Belle2 data

#include <stdio.h>
#include <string>
#include <vector>

#include "TFile.h"

#include "constants.h"
#include "base.h"
#include "correctors.h"

# define MCTYPE "MC15rd"

// global variables for the simplicity
double Ntotal = 0;
double totalTrueKaon = 0;
double totalTruePion = 0;
double totalTrueMyChargedKaon = 0;
double totalTrueMyChargedPion = 0;
double totalTrueKaonID = 0;
double totalTruePionID = 0;

Corrector corrector;
Corrector_Knn corrector_Knn;
Corrector_Fragmentation corrector_Fragmentation;

/* ====================================== */

void GetNominalNevt(const char* dirname, const char* included_string, const char* type, const char* MC_version, const char* category, double weight_var = 1.0) { // get nominal PDF with appropriate correction
    /*
    CorrectionType for new form factors
    B2Knunu
    B02K0nunu
    B2Xsnunu
    B02Xsnunu
    otherwise
    */

    int Decay[N_decay] = { 0 };

    int NTrueKaon = 0;
    int NTruePion = 0;
    int NTrueMyChargedKaon = 0;
    int NTrueMyChargedPion = 0;
    int NTrueKaonID = 0;
    int NTruePionID = 0;

    double MC_MXs = -1;
    double qsquared = -1;

    std::vector<std::string> names;
    load_files(dirname, &names, included_string);

    for (unsigned int i = 0; i < names.size(); i++) {

        double Nevt = 0;
        std::string filename = names.at(i);

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");

        if (strcmp(type, "SIGNAL") == 0) {
            tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKcharge_total__bc", &Decay[0]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &Decay[1]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &Decay[2]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsu__clMCcomb__bc", &Decay[3]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsu__clMCch1__bc", &Decay[4]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsu__clMCch2__bc", &Decay[5]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsu__clMCch3__bc", &Decay[6]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsu__clMCch4__bc", &Decay[7]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsu__clMCch5__bc", &Decay[8]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsu__clMCch6__bc", &Decay[9]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsu__clMCch7__bc", &Decay[10]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsu__clMCch8__bc", &Decay[11]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsu__clMCch9__bc", &Decay[12]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsu__clMCch10__bc", &Decay[13]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsu__clMCch11__bc", &Decay[14]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsu__clMCch12__bc", &Decay[15]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsu__clMCch13__bc", &Decay[16]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsu__clMCch14__bc", &Decay[17]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsu__clMCch15__bc", &Decay[18]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKneutral_total__bc", &Decay[19]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &Decay[20]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &Decay[21]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsd__clMCcomb__bc", &Decay[22]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsd__clMCch16__bc", &Decay[23]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsd__clMCch17__bc", &Decay[24]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsd__clMCch18__bc", &Decay[25]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsd__clMCch19__bc", &Decay[26]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsd__clMCch20__bc", &Decay[27]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsd__clMCch21__bc", &Decay[28]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsd__clMCch22__bc", &Decay[29]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsd__clMCch23__bc", &Decay[30]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsd__clMCch24__bc", &Decay[31]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsd__clMCch25__bc", &Decay[32]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsd__clMCch26__bc", &Decay[33]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsd__clMCch27__bc", &Decay[34]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsd__clMCch28__bc", &Decay[35]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsd__clMCch29__bc", &Decay[36]);
            tree_upsilon->SetBranchAddress("nParticlesInList__boXsd__clMCch30__bc", &Decay[37]);

            tree_upsilon->SetBranchAddress("mcDaughter__bo0__cm__spmcDaughter__bo0__cm__spM__bc__bc", &MC_MXs);
            tree_upsilon->SetBranchAddress("qsquared", &qsquared);
        }

        tree_upsilon->SetBranchAddress("nParticlesInList__boK__pl__clKaon_true__bc", &NTrueKaon);
        tree_upsilon->SetBranchAddress("nParticlesInList__bopi__pl__clPion_true__bc", &NTruePion);
        tree_upsilon->SetBranchAddress("nParticlesInList__boK__pl__clmychargedKaon_true__bc", &NTrueMyChargedKaon);
        tree_upsilon->SetBranchAddress("nParticlesInList__bopi__pl__clmychargedPion_true__bc", &NTrueMyChargedPion);
        tree_upsilon->SetBranchAddress("nParticlesInList__boK__pl__clKaonID_true__bc", &NTrueKaonID);
        tree_upsilon->SetBranchAddress("nParticlesInList__bopi__pl__clPionID_true__bc", &NTruePionID);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);

            double total_weight = weight_var;

            if (strcmp(type, "SIGNAL") == 0) {
                if (filename.find("B2Knunu") != std::string::npos) {
                    double correction_weight = corrector.GetCorrectionFactor(qsquared, "Bplus");
                    total_weight = total_weight * ObtainWeight(type, MC_version, category, filename) * correction_weight;
                }
                else if (filename.find("B2Kstarnunu") != std::string::npos) total_weight = total_weight * ObtainWeight(type, MC_version, category, filename);
                else if (filename.find("B2Xsnunu") != std::string::npos) {
                    double correction_fragmentation = corrector_Fragmentation.GetCorrectionFactor(Decay, MC_MXs, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MC_version);
                    total_weight = total_weight * ObtainWeight(type, MC_version, category, filename) * correction_fragmentation;
                }
                else if (filename.find("B02K0nunu") != std::string::npos) {
                    double correction_weight = corrector.GetCorrectionFactor(qsquared, "Bzero");
                    total_weight = total_weight * ObtainWeight(type, MC_version, category, filename) * correction_weight;
                }
                else if (filename.find("B02Kstar0nunu") != std::string::npos) total_weight = total_weight * ObtainWeight(type, MC_version, category, filename);
                else if (filename.find("B02Xsnunu") != std::string::npos) {
                    double correction_fragmentation = corrector_Fragmentation.GetCorrectionFactor(Decay, MC_MXs, Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, MC_version);
                    total_weight = total_weight * ObtainWeight(type, MC_version, category, filename) * correction_fragmentation;
                }
                else { total_weight = total_weight * ObtainWeight(type, MC_version, category, filename); }
            }
            else total_weight = total_weight * ObtainWeight(type, MC_version, category, filename) * corrector_Knn.GetCorrectionFactorCancelOutObtainWeight(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, filename, MC_version, false);

            Nevt = Nevt + total_weight;

            totalTrueKaon = totalTrueKaon + NTrueKaon * total_weight;
            totalTruePion = totalTruePion + NTruePion * total_weight;
            totalTrueMyChargedKaon = totalTrueMyChargedKaon + NTrueMyChargedKaon * total_weight;
            totalTrueMyChargedPion = totalTrueMyChargedPion + NTrueMyChargedPion * total_weight;
            totalTrueKaonID = totalTrueKaonID + NTrueKaonID * total_weight;
            totalTruePionID = totalTruePionID + NTruePionID * total_weight;

        }
        input_file->Close();

        printf("%s has %lf events (with correction)\n", filename.c_str(), Nevt);

        Ntotal = Ntotal + Nevt;

    }

    return;
}

int Calculate_PID_efficiency() {

    const char* dirname = "";

    GetNominalNevt(dirname, "B2Knunu", "SIGNAL", MCTYPE, "validation", 1.0);
    GetNominalNevt(dirname, "B2Kstarnunu", "SIGNAL", MCTYPE, "validation", 1.0);
    GetNominalNevt(dirname, "B2Xsnunu", "SIGNAL", MCTYPE, "validation", 1.0);
    GetNominalNevt(dirname, "B02K0nunu", "SIGNAL", MCTYPE, "validation", 1.0);
    GetNominalNevt(dirname, "B02Kstar0nunu", "SIGNAL", MCTYPE, "validation", 1.0);
    GetNominalNevt(dirname, "B02Xsnunu", "SIGNAL", MCTYPE, "validation", 1.0);

    printf("Ntotal: %lf\n", Ntotal);
    printf("totalTrueKaon: %lf\n", totalTrueKaon);
    printf("totalTruePion: %lf\n", totalTruePion);
    printf("totalTrueMyChargedKaon: %lf\n", totalTrueMyChargedKaon);
    printf("totalTrueMyChargedPion: %lf\n", totalTrueMyChargedPion);
    printf("totalTrueKaonID: %lf\n", totalTrueKaonID);
    printf("totalTruePionID: %lf\n", totalTruePionID);

    return 0;
}
