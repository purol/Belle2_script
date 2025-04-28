// last update: 2021-10-13
// for Belle2 data

# define N_decay 38 // five decay mode + others

# define Nentry 200
# define MIN_MXs 1.1
# define MAX_MXs 5.3

# define Nvar 23
# define DvetoNvar 4

# define MyEPSILON 0.000001
# define MCTYPE "MC15ri"

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

void LetsFillMXs(const char* dirname, std::vector<double>* list) {

    std::vector<string> names;
    load_files(dirname, &names);

    double MXs = -1;

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");

        tree_Bsig->SetBranchAddress("Bsig_M", &MXs);

        printf("%lld entries...\n", tree_Bsig->GetEntries());
        for (unsigned int j = 0; j < tree_Bsig->GetEntries(); j++) { // Fill
            tree_Bsig->GetEntry(j);

            list->push_back(MXs);

        }
        input_file->Close();

    }

}

void ReadMXs(){

    const char* dirname = "./";

    std::vector<double> MXs_list;

    LetsFillMXs(dirname, &MXs_list);

    FILE* fp;
    fp = fopen("MXs_data.csv","w");
    for (int i = 0; i < MXs_list.size(); i++) {
        fprintf(fp, "%lf\n", MXs_list.at(i));
    }
    fclose(fp);


}
