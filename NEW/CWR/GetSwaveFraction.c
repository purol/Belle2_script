// last update: 2021-10-13
// for Belle2 data

#include <vector>

# define Nentry 2000
# define MIN_Q2 0.0
# define MAX_Q2 20.0
# define Nfine 1000000

typedef struct {
    double q_square;
    double MKstar;
} Data;

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

void load_files(const char* dirname, std::vector<std::string>* names, const char* included_string) {
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

void LetsFillMK0star(const char* dirname, std::vector<Data>* list const char* included_string) {

    double invM = 0.0;
    double MXs = 0.0;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_Xs = (TTree*)input_file->Get("Xs");

        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);
        tree_Xs->SetBranchAddress("daughter__bo1__cm__spdaughter__bo0__cm__spM__bc__bc", &MXs);

        printf("%lld entries...\n", tree_Xs->GetEntries());
        for (unsigned int j = 0; j < tree_Xs->GetEntries(); j++) { // Fill
            tree_Xs->GetEntry(j);

            Data temp_data = { invM * invM , MXs };
            list->push_back(temp_data);

        }
        input_file->Close();

    }

}

void GetSwaveFraction(){

    const char* K0star_dirname = "/home/belle2/junewoo/storage_ghi/20220929_SIGNAL_decayInfo_again/small";

    std::vector<Data> list;

    LetsFillMK0star(K0star_dirname, &list, "B02Kstar0nunu");

    double N_total_796_996 = 0;
    double N_total_644_1200 = 0;

    double avg_FS_796_996 = 0;
    double avg_FS_644_1200 = 0;

    for (int i = 0; i < list.size(); i++) {
        if ((0.796 < list.at(i).MKstar) && (list.at(i).MKstar < 0.996)) {
            N_total_796_996 = N_total_796_996 + 1.0;

            if ((list.at(i).q_square < 0.98)) avg_FS_796_996 = avg_FS_796_996 + 1.0 * 0.021; // assume
            else if ((0.98 < list.at(i).q_square) && (list.at(i).q_square < 1.1)) avg_FS_796_996 = avg_FS_796_996 + 1.0 * (0.021 + 0.144) / 2.0; // interpolate
            else if ((1.1 < list.at(i).q_square) && (list.at(i).q_square < 2.5)) avg_FS_796_996 = avg_FS_796_996 + 1.0 * 0.144;
            else if ((2.5 < list.at(i).q_square) && (list.at(i).q_square < 4.0)) avg_FS_796_996 = avg_FS_796_996 + 1.0 * 0.029;
            else if ((4.0 < list.at(i).q_square) && (list.at(i).q_square < 6.0)) avg_FS_796_996 = avg_FS_796_996 + 1.0 * 0.117;
            else if ((6.0 < list.at(i).q_square) && (list.at(i).q_square < 8.0)) avg_FS_796_996 = avg_FS_796_996 + 1.0 * 0.033;
            else if ((8.0 < list.at(i).q_square) && (list.at(i).q_square < 11.0)) avg_FS_796_996 = avg_FS_796_996 + 1.0 * (0.033 + 0.021) / 2.0; // interpolate
            else if ((11.0 < list.at(i).q_square) && (list.at(i).q_square < 12.5)) avg_FS_796_996 = avg_FS_796_996 + 1.0 * 0.021;
            else if ((12.5 < list.at(i).q_square) && (list.at(i).q_square < 15.0)) avg_FS_796_996 = avg_FS_796_996 + 1.0 * (0.021 - 0.008) / 2.0; // interpolate
            else if ((15.0 < list.at(i).q_square) && (list.at(i).q_square < 17.0)) avg_FS_796_996 = avg_FS_796_996 + 1.0 * (-0.008);
            else if (17.0 < list.at(i).q_square) avg_FS_796_996 = avg_FS_796_996 + 1.0 * 0.144; // assume

        }



        if ((0.644 < list.at(i).MKstar) && (list.at(i).MKstar < 1.200)) {
            N_total_644_1200 = N_total_644_1200 + 1.0;

            if ((list.at(i).q_square < 0.98)) avg_FS_644_1200 = avg_FS_644_1200 + 1.0 * 0.021; // assume
            else if ((0.98 < list.at(i).q_square) && (list.at(i).q_square < 1.1)) avg_FS_644_1200 = avg_FS_644_1200 + 1.0 * (0.021 + 0.144) / 2.0; // interpolate
            else if ((1.1 < list.at(i).q_square) && (list.at(i).q_square < 2.5)) avg_FS_644_1200 = avg_FS_644_1200 + 1.0 * 0.144;
            else if ((2.5 < list.at(i).q_square) && (list.at(i).q_square < 4.0)) avg_FS_644_1200 = avg_FS_644_1200 + 1.0 * 0.029;
            else if ((4.0 < list.at(i).q_square) && (list.at(i).q_square < 6.0)) avg_FS_644_1200 = avg_FS_644_1200 + 1.0 * 0.117;
            else if ((6.0 < list.at(i).q_square) && (list.at(i).q_square < 8.0)) avg_FS_644_1200 = avg_FS_644_1200 + 1.0 * 0.033;
            else if ((8.0 < list.at(i).q_square) && (list.at(i).q_square < 11.0)) avg_FS_644_1200 = avg_FS_644_1200 + 1.0 * (0.033 + 0.021) / 2.0; // interpolate
            else if ((11.0 < list.at(i).q_square) && (list.at(i).q_square < 12.5)) avg_FS_644_1200 = avg_FS_644_1200 + 1.0 * 0.021;
            else if ((12.5 < list.at(i).q_square) && (list.at(i).q_square < 15.0)) avg_FS_644_1200 = avg_FS_644_1200 + 1.0 * (0.021 - 0.008) / 2.0; // interpolate
            else if ((15.0 < list.at(i).q_square) && (list.at(i).q_square < 17.0)) avg_FS_644_1200 = avg_FS_644_1200 + 1.0 * (-0.008);
            else if (17.0 < list.at(i).q_square) avg_FS_644_1200 = avg_FS_644_1200 + 1.0 * 0.144; // assume
        }

    }

    avg_FS_796_996 = avg_FS_796_996 / N_total_796_996;
    avg_FS_644_1200 = avg_FS_644_1200 / N_total_644_1200;

    printf("%lf %lf\n", avg_FS_796_996, avg_FS_644_1200);


}
