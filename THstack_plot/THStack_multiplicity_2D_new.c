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

#include "constants.h"
#include <string>

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

typedef struct Event_ {
    int __experiment__ = -1;
    int __run__ = -1;
    unsigned int __event__ = -1;
    int __ncandidates__ = -1;
    double MXs_reco;
    double weight;
    std::string original_filename; // ex) B2Xsnunu_11180500_01354_job362894809_00.udst
} Event;

void LetsFill(const char* dirname, std::vector<Event>* event_vector, const char* included_string, bool IsItBplus, double weight = 1) {

    double MXs = 0;
    double SignalProbability = -1;

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    for (unsigned int i = 0; i < names.size(); i++) {

        int __experiment__ = -1;
        int __run__ = -1;
        unsigned int __event__ = -1;
        int __ncandidates__ = -1;

        double Mxs_Bc_MC = -1;
        double Mxs_B0_MC = -1;
        double MC_MXs = -1;

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");
        TTree* tree_Xs = (TTree*)input_file->Get("Xs");

        tree_Bsig->SetBranchAddress("Bsig_M", &MXs);
        tree_Btag->SetBranchAddress("Btag_extraInfo_SignalProbability", &SignalProbability);
        tree_upsilon->SetBranchAddress("__experiment__", &__experiment__);
        tree_upsilon->SetBranchAddress("__run__", &__run__);
        tree_upsilon->SetBranchAddress("__event__", &__event__);
        tree_upsilon->SetBranchAddress("__ncandidates__", &__ncandidates__);

        tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
        tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);

        // define vector to save MXs^reco
        std::vector<double> MXs_recos;
        std::vector<double> Signal_Probabilities;

        // extract original filename
        std::string original_filename;
        size_t pos = names.at(i).find(".udst");
        if (pos != std::string::npos) {
            // Extract the substring until ".udst"
            original_filename = names.at(i).substr(0, pos + 5); // Include ".udst"
        }
        else {
            printf("'.udst' not found in the string.\n");
            exit(1);
        }

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            tree_Xs->GetEntry(j);

            if (IsItBplus) MC_MXs = Mxs_Bc_MC;
            else MC_MXs = Mxs_B0_MC;

            // D veto
            if ((MXs > 1.84) && (MXs < 1.89)) continue;

            Event temp_ = { __experiment__, __run__, __event__, __ncandidates__, MXs, original_filename, weight };
            event_vector->push_back(temp_);

        }

        input_file->Close();

    }

}

void FillMatrix(double Ncandidates[3][3], std::vector<Event> before_BCS_Event, std::vector<Event> after_BCS_Event) {
    for (int i = 0; i < after_BCS_Event.size(); i++) {
        for (int j = 0; j < before_BCS_Event.size(); j++) {
            if ((after_BCS_Event.at(i).__experiment__ == before_BCS_Event.at(j).__experiment__) && (after_BCS_Event.at(i).__run__ == before_BCS_Event.at(j).__run__) && (after_BCS_Event.at(i).__event__ == before_BCS_Event.at(j).__event__) && (after_BCS_Event.at(i).__ncandidates__ == before_BCS_Event.at(j).__ncandidates__) && (after_BCS_Event.at(i).original_filename == before_BCS_Event.at(j).original_filename)) {
                int before_MXs_bin = -1;
                int after_MXs_bin = -1;

                if ((before_BCS_Event.at(j).MXs_reco > 0) && (before_BCS_Event.at(j).MXs_reco < 0.6)) before_MXs_bin = 0;
                else if ((before_BCS_Event.at(j).MXs_reco > 0.6) && (before_BCS_Event.at(j).MXs_reco < 1.0)) before_MXs_bin = 1;
                else before_MXs_bin = 2;

                if ((after_BCS_Event.at(i).MXs_reco > 0) && (after_BCS_Event.at(i).MXs_reco < 0.6)) after_MXs_bin = 0;
                else if ((after_BCS_Event.at(i).MXs_reco > 0.6) && (after_BCS_Event.at(i).MXs_reco < 1.0)) after_MXs_bin = 1;
                else after_MXs_bin = 2;

                Ncandidates[before_MXs_bin][after_MXs_bin] = Ncandidates[before_MXs_bin][after_MXs_bin] + before_BCS_Event.at(j).weight;

            }
        }
    }
}

void THStack_multiplicity_2D_new() {

    double Ncandidates[3][3] = { 0.0 }; // [MXs^reco before BCS][MXs^reco after BCS]
    std::vector<Event> before_BCS_Event;
    std::vector<Event> after_BCS_Event;

    const char* Knunu_before_BCS_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD/SIGNAL_analysis/validation_v004/before_Dveto_cut";
    const char* Kstarnunu_before_BCS_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD/SIGNAL_analysis/validation_v004/before_Dveto_cut";
    const char* Xsununu_before_BCS_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD/SIGNAL_analysis/validation_v004/before_Dveto_cut";
    const char* K0nunu_before_BCS_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD/SIGNAL_analysis/validation_v004/before_Dveto_cut";
    const char* K0starnunu_before_BCS_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD/SIGNAL_analysis/validation_v004/before_Dveto_cut";
    const char* Xsdnunu_before_BCS_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD/SIGNAL_analysis/validation_v004/before_Dveto_cut";

    const char* Knunu_after_BCS_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD/SIGNAL_analysis/validation_v004/final_output";
    const char* Kstarnunu_after_BCS_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD/SIGNAL_analysis/validation_v004/final_output";
    const char* Xsununu_after_BCS_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD/SIGNAL_analysis/validation_v004/final_output";
    const char* K0nunu_after_BCS_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD/SIGNAL_analysis/validation_v004/final_output";
    const char* K0starnunu_after_BCS_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD/SIGNAL_analysis/validation_v004/final_output";
    const char* Xsdnunu_after_BCS_dirname = "/home/belle2/junewoo/storage_b1/bsub/Analysis/KumoiRD/SIGNAL_analysis/validation_v004/final_output";

    LetsFill(Knunu_before_BCS_dirname, &before_BCS_Event, "B2Knunu", true, Scale_Kplus_validation_MC15rd);
    LetsFill(Kstarnunu_before_BCS_dirname, &before_BCS_Event, "B2Kstarnunu", true, Scale_Kplusstar_validation_MC15rd);
    LetsFill(Xsununu_before_BCS_dirname, &before_BCS_Event, "B2Xsnunu", true, Scale_Xsu_nonresonant_validation_MC15rd);
    LetsFill(K0nunu_before_BCS_dirname, &before_BCS_Event, "B02K0nunu", false, Scale_K0_validation_MC15rd);
    LetsFill(K0starnunu_before_BCS_dirname, &before_BCS_Event, "B02Kstar0nunu", false, Scale_K0star_validation_MC15rd);
    LetsFill(Xsdnunu_before_BCS_dirname, &before_BCS_Event, "B02Xsnunu", false, Scale_Xsd_nonresonant_validation_MC15rd);

    LetsFill(Knunu_after_BCS_dirname, &after_BCS_Event, "B2Knunu", true, Scale_Kplus_validation_MC15rd);
    LetsFill(Kstarnunu_after_BCS_dirname, &after_BCS_Event, "B2Kstarnunu", true, Scale_Kplusstar_validation_MC15rd);
    LetsFill(Xsununu_after_BCS_dirname, &after_BCS_Event, "B2Xsnunu", true, Scale_Xsu_nonresonant_validation_MC15rd);
    LetsFill(K0nunu_after_BCS_dirname, &after_BCS_Event, "B02K0nunu", false, Scale_K0_validation_MC15rd);
    LetsFill(K0starnunu_after_BCS_dirname, &after_BCS_Event, "B02Kstar0nunu", false, Scale_K0star_validation_MC15rd);
    LetsFill(Xsdnunu_after_BCS_dirname, &after_BCS_Event, "B02Xsnunu", false, Scale_Xsd_nonresonant_validation_MC15rd);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%lf ", Ncandidates[j][i]);
        }
        printf("\n");
    }

    //    free(line);
    //    gPad->BuildLegend();
    //    c_temp->SaveAs("Plot_with_legend.png");
}
