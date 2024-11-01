// last update: 2023-07-14
// for Belle2 data

#include <stdio.h>
#include <string>
#include <vector>

#include "TFile.h"

int ReadnuisanceparameterPullTypeone(int argc, char* argv[]) {

    std::vector<double> pulls;

    FILE* fp = fopen("number.txt", "r");

    double pull = -1;
    while (true) {
        if (fscanf(fp, "%lf\n", &pull) == EOF) break;
        pulls.push_back(pull);
        printf("%lf\n", pull);
    }

    fclose(fp);

    return 0;
}
