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

#include "constants.h"
#include "base.h"
#include "ObtainWeight.h"
#include "Loader.h"

/*
when you add new variables:
revise # define N_Needed_info ...
revise typedef struct dataLoader::ConvertIntoSeparateDataFile(std::string output_name, double flag = 0)
revise void Loader::GetData(TFile* input_file)
check MakeShiftDoubleToInt if it affects the code
revise void Loader::PrintRootFile(std::string output_name)
revise void Loader::PrintSeparateRootFile(std::string output_name)
revise void Loader::ConvertIntoSeparateDataFile(std::string output_name, double flag = 0)
*/

char* MCTYPE;

/* ====================================== */

int main(int argc, char* argv[]) {
    /*
    * argv[1]: dirname
    * argv[2]: output path
    * argv[3]: sample type: {SIGNAL|CHG|MIX|UUBAR|DDBAR|SSBAR|CHARM}
    * argv[4]: MC version: {data|MC15ri|MC15rd}
    * argv[5]: sample category: {test|train|validation}
    */

    MCTYPE = argv[4];

    std::vector<std::string> names;

    load_files(argv[1], &names);

    Loader loader;

    for (unsigned int i = 0; i < names.size(); i++) {
        loader.initialize();

        TFile* input_file = new TFile((argv[1] + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());
        loader.GetData(input_file);
        if (loader.event_info_is_valid() == false) { printf("error!\n"); return 0; }

        std::string::size_type const p(names.at(i).find_last_of('.'));
        std::string file_without_extension = names.at(i).substr(0, p);

        loader.PrintInformation(std::string("========== inital =========="), names.at(i), argv[3], MCTYPE, argv[5], true);
        loader.MVACut(0.86, -1.0, Loader::KaonMass);
        loader.MVACut(0.86, -1.0, Loader::KstarMass);
        loader.MVACut(0.86, -1.0, Loader::XsMass);
        loader.PrintInformation(std::string("========== MVA CUT =========="), names.at(i), argv[3], MCTYPE, argv[5], true);

        loader.PrintSeparateRootFile(std::string(argv[2]) + "/final_output_root_after_MVA_Application_after_cut/" + file_without_extension + std::string("_after_MVA_cut.root"));
    }
    loader.End();

    return 0;
}
