#include <stdio.h>
#include <string>
#include <vector>
#include <map>

#include "TFile.h"

#include "Loader.h"
#include "constants.h"
#include "MyObtainWeight.h"
#include "MyModule.h"



int main(int argc, char* argv[]) {

    Loader loader("pi0");

    loader.Load("./", "root", "label");

    loader.PrintInformation("========== initial ==========");
    loader.Cut("(0.5 < isSignal)");
    loader.PrintInformation("========== isSignal > 0.5 ==========");
    loader.DrawTH1D("M", "M", 100, 0.11, 0.155, "M.png");
    loader.Cut("(0.120 < M) && (M < 0.145)");
    loader.PrintInformation("========== 0.12 < M < 0.145 ==========");

    loader.end();

    return 0;
}
