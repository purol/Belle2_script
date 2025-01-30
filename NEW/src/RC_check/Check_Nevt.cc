// last update: 2023-07-14
// for Belle2 data

#include <stdio.h>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"

#include "constants.h"
#include "base.h"

/*
when you add new variables:
revise # define N_Needed_info ...
revise typedef struct data
revise void Loader::GetData(TFile* input_file)
check MakeShiftDoubleToInt if it affects the code
revise void Loader::PrintRootFile(std::string output_name)
revise void Loader::PrintSeparateRootFile(std::string output_name)
revise void Loader::ConvertIntoSeparateDataFile(std::string output_name, double flag = 0)
*/

/* ====================================== */

int main(int argc, char* argv[]) {
    /*
    * argv[1]: dirname
    */

    std::vector<std::string> names;

    load_files(argv[1], &names);

    FILE* fp = fopen("empty.txt", "w");

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((argv[1] + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* temp_tree = (TTree*)input_file->Get("data");

        if (temp_tree->GetEntries() == 0) fprintf(fp, "%s\n", names.at(i).c_str());

    }

    fclose(fp);

    return 0;
}
