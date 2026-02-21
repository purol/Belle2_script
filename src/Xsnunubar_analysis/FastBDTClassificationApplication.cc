#include <string>
#include <vector>

#include "Classifier.h"

#include "constants.h"
#include "base.h"
#include "MyFBDT.h"

int main(int argc, char* argv[]) { 
    /*
    * argv[1]: dirname (ex. /home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_LS_MC_side/CHG_analysis/test_v000)
    * argv[2]: weightfile version name (ex. Aqua, Kokoro, Satori, ...)
    * argv[3]: simple dirname (ex. v000, v001, ...)
    */

    std::string file_path = std::string(argv[1]) + "/final_output_data";
    std::string BB_weightfile_path = "/home/belle2/junewoo/storage_b1/bsub/Analysis/GridSearch/" + std::string(argv[2]) + "_" + std::string(argv[3]);
    std::string Continuum_weightfile_path = "/home/belle2/junewoo/storage_b1/bsub/Analysis/GridSearch_BDTc/" + std::string(argv[2]) + "_" + std::string(argv[3]);
    std::string output_path = std::string(argv[1]) + "/final_output_data_after_MVA_Application";

    std::vector<std::string> names;
    load_files(file_path.c_str(), &names);

    for (unsigned int i = 0; i < names.size(); i++) {
        std::string string_filename = file_path + std::string("/") + names.at(i);
        const char* filename = string_filename.c_str();

        ApplicationEachFile(filename, BB_weightfile_path.c_str(), Continuum_weightfile_path.c_str(), output_path.c_str());
    }

    return 0;
}
