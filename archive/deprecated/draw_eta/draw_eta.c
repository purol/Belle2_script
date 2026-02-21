// last update: 2021-10-08
// for Belle2 data

/*
when you add new variables:
revise # define N_Needed_info ...
revise typedef struct data
revise void Loader::GetData(TFile* input_file)
revise void Loader::PrintRootFile(std::string output_name)
revise void Loader::PrintSeparateRootFile(std::string output_name)
revise void Loader::ConvertIntoSeparateDataFile(std::string output_name, double flag = 0)
*/

# include <algorithm>
# include <float.h>

/* ====================================== */

bool hasEnding(std::string const& fullString, std::string const& ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
    }
    else {
        return false;
    }
}


void load_files(const char *dirname, std::vector<string>* names){
   TSystemDirectory dir(dirname, dirname);
   TList *files = dir.GetListOfFiles();
   if (files) {
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next())) {
         fname = file->GetName();
         if (!file->IsDirectory() && fname.EndsWith(".root")) {
            names->push_back(fname.Data());
         }
      }
   }
}

void load_files(const char* dirname, std::vector<string>* names, const char* included_string) {
    TSystemDirectory dir(dirname, dirname);
    TList* files = dir.GetListOfFiles();
    if (files) {
        TSystemFile* file;
        TString fname;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(".root") && fname.Contains(included_string)) {
                names->push_back(fname.Data());
            }
        }
    }
}

void LetsFill(const char* dirname, const char* branch_name, const char* var_name, TH1D* hist ) {
    double var = 0.0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree = (TTree*)input_file->Get(branch_name);
        tree->SetBranchAddress(var_name, &var);

        printf("%lld entries...\n", tree->GetEntries());
        for (unsigned int j = 0; j < tree->GetEntries(); j++) { // Fill
            tree->GetEntry(j);

            if ((var < 0.45) || (var > 0.6)) continue;
            hist->Fill(var);
        }
        input_file->Close();

    }

    hist->Scale(1.0 / hist->Integral(), "width");

}


void draw_eta() {

    const char* dirname = "./";

    TH1D* eta_1_M = new TH1D("eta_1_M", ";M [GeV];arbitrary unit", 30, 0.45, 0.6);
    TH1D* eta_2_M = new TH1D("eta_2_M", ";M [GeV];arbitrary unit", 30, 0.45, 0.6);
    TH1D* eta_3_M = new TH1D("eta_3_M", ";M [GeV];arbitrary unit", 30, 0.45, 0.6);

    LetsFill(dirname, "eta_1", "M", eta_1_M);
    LetsFill(dirname, "eta_2", "M", eta_2_M);
    LetsFill(dirname, "eta_3", "M", eta_3_M);

    TCanvas* c_temp = new TCanvas("c", "", 800, 800);

    eta_1_M->SetStats(false);
    eta_2_M->SetStats(false);
    eta_3_M->SetStats(false);

    eta_1_M->Draw("hist"); c_temp->SaveAs("eta_1_M.png");
    eta_2_M->Draw("hist"); c_temp->SaveAs("eta_2_M.png");
    eta_3_M->Draw("hist"); c_temp->SaveAs("eta_3_M.png");

    delete c_temp;

}
