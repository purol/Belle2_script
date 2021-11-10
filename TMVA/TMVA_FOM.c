// last update: 2021-11-10-00
// for Belle2 data
# define Nstep 20
# define start 0.5
# define end 1.0
# define NormFactorForSignal (31481.0/9300000.0)

std::vector<double> SIGNAL_BB_outputs;
std::vector<double> SIGNAL_Continuum_outputs;
std::vector<double> BKG_BB_outputs;
std::vector<double> BKG_Continuum_outputs;

float BB_output = 0;
float Continuum_output = 0;

void Fill_TMVA_outputs(const char * name, bool IsItSignal) {

    TFile* input_file = new TFile(name, "read");
    TTree* temp_tree = (TTree*)input_file->Get("data");
    temp_tree->SetBranchAddress("TMVA_BB", &BB_output);
    temp_tree->SetBranchAddress("TMVA_Continuum", &Continuum_output);
    for (unsigned int j = 0; j < temp_tree->GetEntries(); j++) { // Fill
        temp_tree->GetEntry(j);
        if (IsItSignal) {
            SIGNAL_BB_outputs.push_back(BB_output);
            SIGNAL_Continuum_outputs.push_back(Continuum_output);
        }
        else {
            BKG_BB_outputs.push_back(BB_output);
            BKG_Continuum_outputs.push_back(Continuum_output);
        }
    }
    input_file->Close();

}

void TMVA_FOM() {

    //Fill_TMVA_outputs("/media/sf_virtualbox_folder/20211109/TMVA_output/TMVAoutput_SIGNAL_final_output_merge_Mxs_smaller_test.root", true);
    Fill_TMVA_outputs("/media/sf_virtualbox_folder/20211109/TMVA_output/TMVAoutput_SIGNAL_final_output_merge_Mxs_larger_test.root", true);
    //Fill_TMVA_outputs("/media/sf_virtualbox_folder/20211109/TMVA_output/TMVAoutput_CHG_final_output_merge_Mxs_smaller_test.root", false);
    Fill_TMVA_outputs("/media/sf_virtualbox_folder/20211109/TMVA_output/TMVAoutput_CHG_final_output_merge_Mxs_larger_test.root", false);
    //Fill_TMVA_outputs("/media/sf_virtualbox_folder/20211109/TMVA_output/TMVAoutput_MIX_final_output_merge_Mxs_smaller_test.root", false);
    Fill_TMVA_outputs("/media/sf_virtualbox_folder/20211109/TMVA_output/TMVAoutput_MIX_final_output_merge_Mxs_larger_test.root", false);
    //Fill_TMVA_outputs("/media/sf_virtualbox_folder/20211109/TMVA_output/TMVAoutput_UUBAR_final_output_merge_Mxs_smaller_test.root", false);
    Fill_TMVA_outputs("/media/sf_virtualbox_folder/20211109/TMVA_output/TMVAoutput_UUBAR_final_output_merge_Mxs_larger_test.root", false);
    //Fill_TMVA_outputs("/media/sf_virtualbox_folder/20211109/TMVA_output/TMVAoutput_DDBAR_final_output_merge_Mxs_smaller_test.root", false);
    Fill_TMVA_outputs("/media/sf_virtualbox_folder/20211109/TMVA_output/TMVAoutput_DDBAR_final_output_merge_Mxs_larger_test.root", false);
    //Fill_TMVA_outputs("/media/sf_virtualbox_folder/20211109/TMVA_output/TMVAoutput_SSBAR_final_output_merge_Mxs_smaller_test.root", false);
    Fill_TMVA_outputs("/media/sf_virtualbox_folder/20211109/TMVA_output/TMVAoutput_SSBAR_final_output_merge_Mxs_larger_test.root", false);
    //Fill_TMVA_outputs("/media/sf_virtualbox_folder/20211109/TMVA_output/TMVAoutput_CHARM_final_output_merge_Mxs_smaller_test.root", false);
    Fill_TMVA_outputs("/media/sf_virtualbox_folder/20211109/TMVA_output/TMVAoutput_CHARM_final_output_merge_Mxs_larger_test.root", false);

    // output step: 0.05 0.1 0.15 ... 0.95
    double Matrix[Nstep][Nstep];
    for(int i = 0; i < Nstep; i++) for(int j = 0; i < Nstep; i++) Matrix[i][j] = 0.0;

    for (int i = 0; i < Nstep; i++) {
        for (int j = 0; j < Nstep; j++) {
            double BB_output = start + (end - start) * i / Nstep;
            double Continuum_output = start + (end - start) * j / Nstep;

            double SIGNAL_num = 0.0;
            double BKG_num = 0.0;

            for (size_t k = 0; k < SIGNAL_BB_outputs.size(); ++k) if (SIGNAL_BB_outputs.at(k) > BB_output && SIGNAL_Continuum_outputs.at(k) > Continuum_output) SIGNAL_num = SIGNAL_num + 1.0;
            for (size_t k = 0; k < BKG_BB_outputs.size(); ++k) if (BKG_BB_outputs.at(k) > BB_output && BKG_Continuum_outputs.at(k) > Continuum_output) BKG_num = BKG_num + 1.0;

            double Norm_SIGNAL_num = NormFactorForSignal * SIGNAL_num;

            Matrix[i][j] = Norm_SIGNAL_num / std::sqrt(Norm_SIGNAL_num + BKG_num);
        }
    }

    for (int i = 0; i < Nstep; i++) {
        for (int j = 0; j < Nstep; j++) {
            printf("%lf ", Matrix[i][j]);
        }
        printf("\n");
    }

}
