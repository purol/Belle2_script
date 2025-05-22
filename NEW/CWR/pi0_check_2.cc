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
    loader.DrawTH1D("EE", "EE", 100, 0.02, 0.05, "E_daughter.png");

    loader.end();

    return 0;
}
