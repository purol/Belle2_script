// last update: 2023-07-14
// for Belle2 data

#include <iostream>
#include <vector>
#include <numeric>  // for std::accumulate
#include <cmath>    // for std::sqrt

/*
when you add new variables:
revise # define N_Needed_info ...
revise typedef struct
revise void Loader::GetData(TFile* input_file)
check MakeShiftDoubleToInt if it affects the code
revise void Loader::PrintRootFile(std::string output_name)
revise void Loader::PrintSeparateRootFile(std::string output_name)
revise void Loader::ConvertIntoSeparateDataFile(std::string output_name, double flag = 0)
*/

/* ====================================== */

void Load(const char* filename, std::vector<double>* variable) {
    double temp = 0.0;

    TFile* input_file = new TFile(filename, "read");

    TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
    TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
    TTree* tree_Btag = (TTree*)input_file->Get("Btag");

    // Bsig_useCMSFrame_p
    // Bsig_M
    tree_Bsig->SetBranchAddress("Bsig_useCMSFrame_E", &temp);

    printf("%lld entries...\n", tree_upsilon->GetEntries());
    for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
        tree_upsilon->GetEntry(j);
        tree_Bsig->GetEntry(j);
        tree_Btag->GetEntry(j);

        variable->push_back(temp);
    }
    input_file->Close();
}

int Get_mean_std() {
    std::vector<double> variable;

    const char* filename = "./CHG_Merge.root";

    Load(filename, &variable);

    // Compute mean
    double sum = std::accumulate(variable.begin(), variable.end(), 0.0);
    double mean = sum / variable.size();

    // Compute standard deviation (population)
    double sq_sum = 0.0;
    for (double x : variable) {
        sq_sum += (x - mean) * (x - mean);
    }
    double stdev = std::sqrt(sq_sum / variable.size());

    std::cout << "Mean = " << mean << "\n";
    std::cout << "Standard deviation = " << stdev << "\n";

    return 0;
}
