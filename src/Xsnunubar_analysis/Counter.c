// last update: 2023-07-14
// for Belle2 data

#include <stdio.h>
#include <string>
#include <vector>

#include "TFile.h"

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

/* ====================================== */

int Counter() {
    /*
    * argv[1]: dirname
    * argv[2]: outputname
    * argv[3]: output path
    * argv[4]: sample type: {SIGNAL|CHG|MIX|UUBAR|DDBAR|SSBAR|CHARM}
    * argv[5]: MC version: {data|MC15ri|MC15rd}
    * argv[6]: sample category: {test|train|validation}
    */

    std::vector<std::string> names;

    load_files("./final_output_root_after_MVA_Application_after_cut/Merge", &names);

    Loader loader;

    for (unsigned int i = 0; i < names.size(); i++) {
        loader.initialize();

        TFile* input_file = new TFile(("./final_output_root_after_MVA_Application_after_cut/Merge" + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());
        loader.GetData(input_file);
        if (loader.event_info_is_valid() == false) { printf("error!\n"); return 0; }

        std::string::size_type const p(names.at(i).find_last_of('.'));
        std::string file_without_extension = names.at(i).substr(0, p);

        loader.PrintInformation(std::string("========== inital =========="), names.at(i), "SIGNAL", "MC15rd", "validation", true);

        loader.PrintConfusionMatrix(names.at(i), "SIGNAL", "MC15rd", "validation", true);

    }
    loader.End();

    return 0;
}
