#ifndef SYSTEMATIC_H
#define SYSTEMATIC_H

#include "constants.h"
#include <vector>
#include <algorithm>

// define struct and sorting function for a binary search
typedef struct EventInformation {
    int experiment;
    int run;
    unsigned int event;
    int candidate;
    int ncandidates;
    int DMID1;
    int DMID2;
} EvtInfo;

bool compare(EvtInfo first, EvtInfo second) {
    if (first.experiment > second.experiment) return true;
    else if (first.experiment < second.experiment) return false;
    else {
        if (first.run > second.run) return true;
        else if (first.run < second.run) return false;
        else {
            if (first.event > second.event) return true;
            else if (first.event < second.event) return false;
            else {
                if (first.candidate > second.candidate) return true;
                else if (first.candidate < second.candidate) return false;
                else {
                    if (first.ncandidates > second.ncandidates) return true;
                    else return false;
                }
            }
        }
    }
}

std::vector<EvtInfo> EvtInfos;

// hard coded uncertainty
struct _BRuncertainty {
    std::vector<int> DMID;
    std::vector<double> RelativeUncertainty;
} BRuncertainty = {
    {
        -1,

        0, 1, 2, 3, 4,
        5, 6, 7, 8, 9,
        10, 11, 12, 13, 14,
        15, 16, 17, 18, 19,
        20, 21, 22, 23, 24,
        25, 26,

        100, 101, 102, 103, 104,
        105, 106, 107, 108, 109,
        110, 111, 112, 113, 114,
        115, 116, 117, 118, 119,
        120, 121, 122
    },
    {
        0,

        0.03943, 0.03913, 0.03913, 0.13433, 0.03943,
        0.02778, 0.26316, 0.46667, 0.02896, 0.80000,
        0.71429, 0.17347, 0.32468, 0.10638, 0.22222,
        //1.00000, 0.10000, 0.03809, 0.22018, 0.04336, care about BR (B --> K+ KL KL) seperately
        1.00000, 0.10000, 0.00000, 0.22018, 0.04336,
        0.11650, 0.29091, 0.21053, 0.07381, 0.16296,
        0.01863, 0.20732,

        0.02414, 0.02414, 0.04018, 0.33333, 0.04018,
        0.20769, 0.15789, 0.13235, 0.15341, 0.56140,
        0.04745, 0.21905, 0.03187, 0.05696, 0.48718,
        //0.55000, 0.13750, 0.08333, 0.07910, 0.21622, care about BR (B0 --> KS KL KL) seperately
        0.55000, 0.13750, 0.00000, 0.07910, 0.21622,
        0.04545, 0.11111, 0.07910
    },
};

void ReadEvtFile() {
    printf("[ReadEvtFile] try to read Evt files...\n");

    const char* CHG_Evt_file = "CHG_Evt";
    const char* MIX_Evt_file = "MIX_Evt";
    const char* SIGNAL_Evt_file = "SIGNAL_Evt";

    FILE* fp_CHG_Evt = fopen(CHG_Evt_file, "r");
    FILE* fp_MIX_Evt = fopen(MIX_Evt_file, "r");
    FILE* fp_SIGNAL_Evt = fopen(SIGNAL_Evt_file, "r");

    double temp_experiment = -1;
    double temp_run = -1;
    double temp_event = -1;
    double temp_candidate = -1;
    double temp_ncandidates = -1;
    int temp_DMID1 = -1;
    int temp_DMID2 = -1;

    while (true) {
        if (fscanf(fp_CHG_Evt, "%lf\n", &temp_experiment) == EOF) break;
        if (fscanf(fp_CHG_Evt, "%lf\n", &temp_run) == EOF) break;
        if (fscanf(fp_CHG_Evt, "%lf\n", &temp_event) == EOF) break;
        if (fscanf(fp_CHG_Evt, "%lf\n", &temp_candidate) == EOF) break;
        if (fscanf(fp_CHG_Evt, "%lf\n", &temp_ncandidates) == EOF) break;
        if (fscanf(fp_CHG_Evt, "%d\n", &temp_DMID1) == EOF) break;
        if (fscanf(fp_CHG_Evt, "%d\n", &temp_DMID2) == EOF) break;

        if ((temp_DMID1 > 99) || (temp_DMID2 > 99)) {
            printf("[ReadEvtFile] DMID for CHG is larger than 99!\n");
            exit(1);
        }

        EvtInfo temp_EvtInfo = { (int)lround(temp_experiment), (int)lround(temp_run), (unsigned int)lround(temp_event), (int)lround(temp_candidate), (int)lround(temp_ncandidates), temp_DMID1, temp_DMID2 };
        EvtInfos.push_back(temp_EvtInfo);

    }

    while (true) {
        if (fscanf(fp_MIX_Evt, "%lf\n", &temp_experiment) == EOF) break;
        if (fscanf(fp_MIX_Evt, "%lf\n", &temp_run) == EOF) break;
        if (fscanf(fp_MIX_Evt, "%lf\n", &temp_event) == EOF) break;
        if (fscanf(fp_MIX_Evt, "%lf\n", &temp_candidate) == EOF) break;
        if (fscanf(fp_MIX_Evt, "%lf\n", &temp_ncandidates) == EOF) break;
        if (fscanf(fp_MIX_Evt, "%d\n", &temp_DMID1) == EOF) break;
        if (fscanf(fp_MIX_Evt, "%d\n", &temp_DMID2) == EOF) break;

        if (((temp_DMID1 < 100) && (temp_DMID1 >= 0)) || ((temp_DMID2 < 100) && (temp_DMID2 >= 0))) {
            printf("[ReadEvtFile] DMID for MIX is smaller than 100!\n");
            exit(1);
        }

        EvtInfo temp_EvtInfo = { (int)lround(temp_experiment), (int)lround(temp_run), (unsigned int)lround(temp_event), (int)lround(temp_candidate), (int)lround(temp_ncandidates), temp_DMID1, temp_DMID2 };
        EvtInfos.push_back(temp_EvtInfo);

    }

    while (true) {
        if (fscanf(fp_SIGNAL_Evt, "%lf\n", &temp_experiment) == EOF) break;
        if (fscanf(fp_SIGNAL_Evt, "%lf\n", &temp_run) == EOF) break;
        if (fscanf(fp_SIGNAL_Evt, "%lf\n", &temp_event) == EOF) break;
        if (fscanf(fp_SIGNAL_Evt, "%lf\n", &temp_candidate) == EOF) break;
        if (fscanf(fp_SIGNAL_Evt, "%lf\n", &temp_ncandidates) == EOF) break;
        if (fscanf(fp_SIGNAL_Evt, "%d\n", &temp_DMID1) == EOF) break;
        if (fscanf(fp_SIGNAL_Evt, "%d\n", &temp_DMID2) == EOF) break;

        if ((temp_DMID1 > 0) && (temp_DMID2 > 0)) {
            printf("[ReadEvtFile] both DMID for SIGNAL is larger than 0!\n");
            exit(1);
        }

        EvtInfo temp_EvtInfo = { (int)lround(temp_experiment), (int)lround(temp_run), (unsigned int)lround(temp_event), (int)lround(temp_candidate), (int)lround(temp_ncandidates), temp_DMID1, temp_DMID2 };
        EvtInfos.push_back(temp_EvtInfo);

    }

    sort(EvtInfos.begin(), EvtInfos.end(), compare);

    fclose(fp_CHG_Evt);
    fclose(fp_MIX_Evt);
    fclose(fp_SIGNAL_Evt);
}

double GetBRFluctuation(int experiment, int run, unsigned int event, int candidate, int ncandidates, std::string fname, int dmID, bool IsItUp) {

    if (fname.find("B2Knn_flat") != std::string::npos) return 1.0; // it is B+ --> K+ n nbar special MC. We do not fluctuate it
    else if (fname.find("B2Kstarnn_flat") != std::string::npos) return 1.0; // it is B+ --> K*+ n nbar special MC. We do not fluctuate it
    else if (fname.find("B02K0nn_flat") != std::string::npos) return 1.0; // it is B0 --> K0 n nbar special MC. We do not fluctuate it
    else if (fname.find("B02K0starnn_flat") != std::string::npos) return 1.0; // it is B0 --> K0 n nbar special MC. We do not fluctuate it

    int temp_Evt_DMID1 = -100;
    int temp_Evt_DMID2 = -100;

    EvtInfo temp_EvtInfo = { experiment, run, event, candidate, ncandidates, -100, -100 };
    int temp_index = lower_bound(EvtInfos.begin(), EvtInfos.end(), temp_EvtInfo, compare) - EvtInfos.begin();
    temp_Evt_DMID1 = EvtInfos.at(temp_index).DMID1;
    temp_Evt_DMID2 = EvtInfos.at(temp_index).DMID2;

    if (!((experiment == EvtInfos.at(temp_index).experiment) && (run == EvtInfos.at(temp_index).run) && (event == EvtInfos.at(temp_index).event) && (candidate == EvtInfos.at(temp_index).candidate) && (ncandidates == EvtInfos.at(temp_index).ncandidates))) {
        printf("[GetBRRelativeUncertainty] Cannot find!\n");
        exit(1);
    }

    double Correction = 1.0;

    if (temp_Evt_DMID1 == dmID) {
        if (IsItUp) Correction = Correction * (1.0 + BRuncertainty.RelativeUncertainty.at(std::find(BRuncertainty.DMID.begin(), BRuncertainty.DMID.end(), temp_Evt_DMID1) - BRuncertainty.DMID.begin()));
        else Correction = Correction * (1.0 - BRuncertainty.RelativeUncertainty.at(std::find(BRuncertainty.DMID.begin(), BRuncertainty.DMID.end(), temp_Evt_DMID1) - BRuncertainty.DMID.begin()));
    }
    if (temp_Evt_DMID2 == dmID) {
        if (IsItUp) Correction = Correction * (1.0 + BRuncertainty.RelativeUncertainty.at(std::find(BRuncertainty.DMID.begin(), BRuncertainty.DMID.end(), temp_Evt_DMID2) - BRuncertainty.DMID.begin()));
        else Correction = Correction * (1.0 - BRuncertainty.RelativeUncertainty.at(std::find(BRuncertainty.DMID.begin(), BRuncertainty.DMID.end(), temp_Evt_DMID2) - BRuncertainty.DMID.begin()));
    }

    return Correction;

}

int NBRdmID() {

    // I want to return the number of dmID whose relative uncertainty is not equal to 0
    int number_dmID = 0;

    for (int i = 0; i < BRuncertainty.DMID.size(); i++) {
        if (BRuncertainty.RelativeUncertainty.at(i) > MyEPSILON) number_dmID++;
    }

    return number_dmID;
}

int GetBRdmID(int index) {

    // I want to get dmID whose relative uncertainty is not equal to 0
    int NdmID = NBRdmID();

    if (index >= NdmID) {
        printf("[GetBRdmID] index is larger than or equal to the number of dmID: %d >= %d\n", index, NdmID);
        exit(1);
    }

    for (int i = 0; i < BRuncertainty.DMID.size(); i++) {
        int dmID = BRuncertainty.DMID.at(i);
        double relative_uncertainty = BRuncertainty.RelativeUncertainty.at(i);

        if (relative_uncertainty <= MyEPSILON) continue;

        if (index == 0) return dmID;
        index = index - 1;
    }

    printf("[GetBRdmID] cannot find dmID. Something wrong.\n");
    return -1;

}


#endif 