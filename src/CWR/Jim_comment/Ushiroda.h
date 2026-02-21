#ifndef USHIRODA_H
#define USHIRODA_H

#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include "TSystemDirectory.h"
#include "TList.h"
#include "TSystemFile.h"
#include "TString.h"
#include "TCollection.h"

typedef struct metainfo {
    int __experiment__;
    int __run__;
    unsigned int __event__;
} MetaInfo;

std::vector<MetaInfo> SIGNAL_uncommon = { 
    {17, 265, 1912}, {22, 520, 27}, {12, 5418, 1539}, {22, 187, 314}, 
    {24, 1912, 484}, {18, 1789, 295}, {16, 1002, 6895}, {22, 171, 5203}, 
    {12, 5901, 4716}, {24, 1641, 2096}, {8, 1288, 88}, {18, 203, 584},
    {18, 2908, 8660}, {18, 1314, 774} };

std::vector<MetaInfo> CHG_uncommon = {
    {17,192,387070}, {12,3500,50874}, {12,6315,82426}, {12,3795,83021},
    {22,461,286853} };

std::vector<MetaInfo> MIX_uncommon = {
    {22,600,141873}, {12,3274,44556} };

std::vector<MetaInfo> UUBAR_uncommon = {
    {18,1415,112399}, {24,2148,1037930}, {18,2983,1779445} };

std::vector<MetaInfo> DDBAR_uncommon = {};

std::vector<MetaInfo> SSBAR_uncommon = {
    {18,2392,258143}, {24,1248,193135}, {24,829,36752} };

std::vector<MetaInfo> CHARM_uncommon = {
    {26,1485,1565119}, {22,210,833992}, {26,775,351095}, {12,2841,4258} };

bool IsIsUnCommon(int experiment, int run, unsigned int event, const char* sample) {
    if (strcmp(sample, "SIGNAL") == 0) {
        for (int i = 0; i < SIGNAL_uncommon.size(); i++) {
            if ((SIGNAL_uncommon.at(i).__experiment__ == experiment) && (SIGNAL_uncommon.at(i).__run__ == run) && (SIGNAL_uncommon.at(i).__event__ == event)) return true;
        }
    }
    else if (strcmp(sample, "CHG") == 0) {
        for (int i = 0; i < CHG_uncommon.size(); i++) {
            if ((CHG_uncommon.at(i).__experiment__ == experiment) && (CHG_uncommon.at(i).__run__ == run) && (CHG_uncommon.at(i).__event__ == event)) return true;
        }
    }
    else if (strcmp(sample, "MIX") == 0) {
        for (int i = 0; i < MIX_uncommon.size(); i++) {
            if ((MIX_uncommon.at(i).__experiment__ == experiment) && (MIX_uncommon.at(i).__run__ == run) && (MIX_uncommon.at(i).__event__ == event)) return true;
        }
    }
    else if (strcmp(sample, "UUBAR") == 0) {
        for (int i = 0; i < UUBAR_uncommon.size(); i++) {
            if ((UUBAR_uncommon.at(i).__experiment__ == experiment) && (UUBAR_uncommon.at(i).__run__ == run) && (UUBAR_uncommon.at(i).__event__ == event)) return true;
        }
    }
    else if (strcmp(sample, "DDBAR") == 0) {
        for (int i = 0; i < DDBAR_uncommon.size(); i++) {
            if ((DDBAR_uncommon.at(i).__experiment__ == experiment) && (DDBAR_uncommon.at(i).__run__ == run) && (DDBAR_uncommon.at(i).__event__ == event)) return true;
        }
    }
    else if (strcmp(sample, "SSBAR") == 0) {
        for (int i = 0; i < SSBAR_uncommon.size(); i++) {
            if ((SSBAR_uncommon.at(i).__experiment__ == experiment) && (SSBAR_uncommon.at(i).__run__ == run) && (SSBAR_uncommon.at(i).__event__ == event)) return true;
        }
    }
    else if (strcmp(sample, "CHARM") == 0) {
        for (int i = 0; i < CHARM_uncommon.size(); i++) {
            if ((CHARM_uncommon.at(i).__experiment__ == experiment) && (CHARM_uncommon.at(i).__run__ == run) && (CHARM_uncommon.at(i).__event__ == event)) return true;
        }
    }

    return false;
}

#endif 