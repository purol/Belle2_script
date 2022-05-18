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

// arXiv:1409.4557v2
# define TB0 1.5195 // (Table. 1)
# define TBp 1.6384 // (Table. 1)
# define BR_Kplus_nunubar 0.00000398 // (eq. 10)
# define BR_K0star_nunubar 0.00000919 // (eq. 11)
# define BR_K0_nunubar (BR_Kplus_nunubar*TB0/TBp) // under (eq. 15)
# define BR_Kplusstar_nunubar (BR_K0star_nunubar*TBp/TB0) // under (eq. 15)
# define BR_Xs_nunubar 0.000029 // (eq. 23)
# define BR_Xsu_nonresonant_nunubar (BR_Xs_nunubar - BR_Kplus_nunubar - BR_Kplusstar_nunubar)
# define BR_Xsd_nonresonant_nunubar (BR_Xs_nunubar - BR_K0_nunubar - BR_K0star_nunubar)

// https://confluence.desy.de/pages/viewpage.action?pageId=107054222
# define N_BpBp_1invab 565400000.0
# define N_B0B0_1invab 534600000.0

# define N_Kplus_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Kplus_nunubar)
# define N_Kplusstar_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Kplusstar_nunubar)
# define N_Xsu_nonresonant_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Xsu_nonresonant_nunubar)
# define N_K0_nunubar_1invab (2.0 * N_B0B0_1invab * BR_K0_nunubar)
# define N_K0star_nunubar_1invab (2.0 * N_B0B0_1invab * BR_K0star_nunubar)
# define N_Xsd_nunubar_1invab (2.0 * N_B0B0_1invab * BR_Xsd_nonresonant_nunubar)

// my MC sample number
# define N_Kplus_nunubar 10000000.0
# define N_K0_nunubar 10000000.0
# define N_Kplusstar_nunubar 10000000.0
# define N_K0star_nunubar 10000000.0
# define N_Xsu_nonresonant_nunubar 50000000.0
# define N_Xsd_nonresonant_nunubar 50000000.0

// scale factor for each MC sample
# define Scale_Kplus (N_Kplus_nunubar_1invab/N_Kplus_nunubar)
# define Scale_Kplusstar (N_Kplusstar_nunubar_1invab/N_Kplusstar_nunubar)
# define Scale_Xsu_nonresonant (N_Xsu_nonresonant_nunubar_1invab/N_Xsu_nonresonant_nunubar)
# define Scale_K0 (N_K0_nunubar_1invab/N_K0_nunubar)
# define Scale_K0star (N_K0star_nunubar_1invab/N_K0star_nunubar)
# define Scale_Xsd_nonresonant (N_Xsd_nunubar_1invab/N_Xsd_nonresonant_nunubar)

# define Nvar_num 17

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

void LetsFill(const char* dirname, std::vector<std::string> variable_names, std::vector<std::string> branch_names, std::vector<double> variable_values[Nvar_num]) {
    double var[Nvar_num] = { 0.0 };

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        for (int k = 0; k < (int)variable_names.size(); k++) {
            if(branch_names.at(k) == std::string("Upsilon")) tree_upsilon->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else if (branch_names.at(k) == std::string("Bsig")) tree_Bsig->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else if (branch_names.at(k) == std::string("Btag")) tree_Btag->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else {
                printf("ERROR! \n");
                exit(1);
            }
        }

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            for (int k = 0; k < (int)variable_names.size(); k++) variable_values[k].push_back(var[k]);
        }
        input_file->Close();

    }

}

void LetsFill(const char* dirname, std::vector<std::string> variable_names, std::vector<std::string> branch_names, std::vector<double> variable_values[Nvar_num], const char* included_string) {
    double var[Nvar_num] = { 0.0 };

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        for (int k = 0; k < (int)variable_names.size(); k++) {
            if (branch_names.at(k) == std::string("Upsilon")) tree_upsilon->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else if (branch_names.at(k) == std::string("Bsig")) tree_Bsig->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else if (branch_names.at(k) == std::string("Btag")) tree_Btag->SetBranchAddress(variable_names.at(k).c_str(), &var[k]);
            else {
                printf("ERROR! \n");
                exit(1);
            }
        }

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            for (int k = 0; k < (int)variable_names.size(); k++) variable_values[k].push_back(var[k]);
        }
        input_file->Close();

    }

}

void THStack_plot() {

    const char* Knunu_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori/SIGNAL_analysis/train_v000/final_output";
    const char* Kstarnunu_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori/SIGNAL_analysis/train_v000/final_output";
    const char* Xsununu_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori/SIGNAL_analysis/train_v000/final_output";
    const char* K0nunu_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori/SIGNAL_analysis/train_v000/final_output";
    const char* K0starnunu_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori/SIGNAL_analysis/train_v000/final_output";
    const char* Xsdnunu_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori/SIGNAL_analysis/train_v000/final_output";
    const char* CHG_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori/CHG_analysis/train_v000/final_output";
    const char* MIX_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori/MIX_analysis/train_v000/final_output";
    const char* UUBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori/UUBAR_analysis/train_v000/final_output";
    const char* DDBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori/DDBAR_analysis/train_v000/final_output";
    const char* SSBAR_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori/SSBAR_analysis/train_v000/final_output";
    const char* CHARM_dirname = "/home/jwpark/storage/BKG_gbasf2/Nitori/CHARM_analysis/train_v000/final_output";

    std::vector<std::string> variable_names;
    std::vector<std::string> branch_names;
    variable_names.push_back("nRemainingTracksInEvent"); branch_names.push_back("Upsilon");
    variable_names.push_back("Btag_chiProb"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_extraInfo_SignalProbability"); branch_names.push_back("Btag");
    variable_names.push_back("thrustAxisCosTheta"); branch_names.push_back("Upsilon");
    variable_names.push_back("missingMomentumOfEvent_theta"); branch_names.push_back("Upsilon");
    variable_names.push_back("Btag_deltaE"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_useCMSFrame_theta"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_cosTBTO"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso00"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso01"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso03"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso04"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso10"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso14"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hso24"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hoo1"); branch_names.push_back("Btag");
    variable_names.push_back("Btag_KSFWVariables_hoo3"); branch_names.push_back("Btag");

    int Nvar = static_cast<int>(variable_names.size());
    if (Nvar != Nvar_num) exit(1);
    std::vector<double> Knunu_values[Nvar_num];
    std::vector<double> Kstarnunu_values[Nvar_num];
    std::vector<double> Xsununu_values[Nvar_num];
    std::vector<double> K0nunu_values[Nvar_num];
    std::vector<double> K0starnunu_values[Nvar_num];
    std::vector<double> Xsdnunu_values[Nvar_num];
    std::vector<double> CHG_values[Nvar_num];
    std::vector<double> MIX_values[Nvar_num];
    std::vector<double> UUBAR_values[Nvar_num];
    std::vector<double> DDBAR_values[Nvar_num];
    std::vector<double> SSBAR_values[Nvar_num];
    std::vector<double> CHARM_values[Nvar_num];

    LetsFill(Knunu_dirname, variable_names, branch_names, Knunu_values, "B2Knunu");
    LetsFill(Kstarnunu_dirname, variable_names, branch_names, Kstarnunu_values, "B2Kstarnunu");
    LetsFill(Xsununu_dirname, variable_names, branch_names, Xsununu_values, "B2Xsnunu");
    LetsFill(K0nunu_dirname, variable_names, branch_names, K0nunu_values, "B02K0nunu");
    LetsFill(K0starnunu_dirname, variable_names, branch_names, K0starnunu_values, "B02Kstar0nunu");
    LetsFill(Xsdnunu_dirname, variable_names, branch_names, Xsdnunu_values, "B02Xsnunu");
    LetsFill(CHG_dirname, variable_names, branch_names, CHG_values);
    LetsFill(MIX_dirname, variable_names, branch_names, MIX_values);
    LetsFill(UUBAR_dirname, variable_names, branch_names, UUBAR_values);
    LetsFill(DDBAR_dirname, variable_names, branch_names, DDBAR_values);
    LetsFill(SSBAR_dirname, variable_names, branch_names, SSBAR_values);
    LetsFill(CHARM_dirname, variable_names, branch_names, CHARM_values);

    THStack* Stack[Nvar_num];
    TH1F* SIGNAL_hist[Nvar_num];
    TH1F* CHG_hist[Nvar_num];
    TH1F* MIX_hist[Nvar_num];
    TH1F* UUBAR_hist[Nvar_num];
    TH1F* DDBAR_hist[Nvar_num];
    TH1F* SSBAR_hist[Nvar_num];
    TH1F* CHARM_hist[Nvar_num];

    for (int k = 0; k < (int)variable_names.size(); k++) { // malloc th1f
        std::vector<double> temp_v;
        temp_v.insert(temp_v.end(), Knunu_values[k].begin(), Knunu_values[k].end());
        temp_v.insert(temp_v.end(), Kstarnunu_values[k].begin(), Kstarnunu_values[k].end());
        temp_v.insert(temp_v.end(), Xsununu_values[k].begin(), Xsununu_values[k].end());
        temp_v.insert(temp_v.end(), K0nunu_values[k].begin(), K0nunu_values[k].end());
        temp_v.insert(temp_v.end(), K0starnunu_values[k].begin(), K0starnunu_values[k].end());
        temp_v.insert(temp_v.end(), Xsdnunu_values[k].begin(), Xsdnunu_values[k].end());
        temp_v.insert(temp_v.end(), CHG_values[k].begin(), CHG_values[k].end());
        temp_v.insert(temp_v.end(), MIX_values[k].begin(), MIX_values[k].end());
        temp_v.insert(temp_v.end(), UUBAR_values[k].begin(), UUBAR_values[k].end());
        temp_v.insert(temp_v.end(), DDBAR_values[k].begin(), DDBAR_values[k].end());
        temp_v.insert(temp_v.end(), SSBAR_values[k].begin(), SSBAR_values[k].end());
        temp_v.insert(temp_v.end(), CHARM_values[k].begin(), CHARM_values[k].end());

        double min = *min_element(temp_v.begin(), temp_v.end());
        double max = *max_element(temp_v.begin(), temp_v.end());

        Stack[k] = new THStack(variable_names.at(k).c_str(), (";"+ variable_names.at(k) + ";arbitrary unit").c_str());
        SIGNAL_hist[k] = new TH1F("signal", (";" + variable_names.at(k) + ";arbitrary unit").c_str(), 100, min, max);
        CHG_hist[k] = new TH1F("charge", (";" + variable_names.at(k) + ";arbitrary unit").c_str(), 100, min, max);
        MIX_hist[k] = new TH1F("mix", (";" + variable_names.at(k) + ";arbitrary unit").c_str(), 100, min, max);
        UUBAR_hist[k] = new TH1F("uubar", (";" + variable_names.at(k) + ";arbitrary unit").c_str(), 100, min, max);
        DDBAR_hist[k] = new TH1F("ddbar", (";" + variable_names.at(k) + ";arbitrary unit").c_str(), 100, min, max);
        SSBAR_hist[k] = new TH1F("ssbar", (";" + variable_names.at(k) + ";arbitrary unit").c_str(), 100, min, max);
        CHARM_hist[k] = new TH1F("ccbar", (";" + variable_names.at(k) + ";arbitrary unit").c_str(), 100, min, max);
    }

    for (int k = 0; k < (int)variable_names.size(); k++) { // fill
        for (int i = 0; i < (int)Knunu_values[k].size(); i++) SIGNAL_hist[k]->Fill(Knunu_values[k].at(i), Scale_Kplus);
        for (int i = 0; i < (int)Kstarnunu_values[k].size(); i++) SIGNAL_hist[k]->Fill(Kstarnunu_values[k].at(i), Scale_Kplusstar);
        for (int i = 0; i < (int)Xsununu_values[k].size(); i++) SIGNAL_hist[k]->Fill(Xsununu_values[k].at(i), Scale_Xsu_nonresonant);
        for (int i = 0; i < (int)K0nunu_values[k].size(); i++) SIGNAL_hist[k]->Fill(K0nunu_values[k].at(i), Scale_K0);
        for (int i = 0; i < (int)K0starnunu_values[k].size(); i++) SIGNAL_hist[k]->Fill(K0starnunu_values[k].at(i), Scale_K0star);
        for (int i = 0; i < (int)Xsdnunu_values[k].size(); i++) SIGNAL_hist[k]->Fill(Xsdnunu_values[k].at(i), Scale_Xsd_nonresonant);

        for (int i = 0; i < (int)CHG_values[k].size(); i++) CHG_hist[k]->Fill(CHG_values[k].at(i));
        for (int i = 0; i < (int)MIX_values[k].size(); i++) MIX_hist[k]->Fill(MIX_values[k].at(i));
        for (int i = 0; i < (int)UUBAR_values[k].size(); i++) UUBAR_hist[k]->Fill(UUBAR_values[k].at(i));
        for (int i = 0; i < (int)DDBAR_values[k].size(); i++) DDBAR_hist[k]->Fill(DDBAR_values[k].at(i));
        for (int i = 0; i < (int)SSBAR_values[k].size(); i++) SSBAR_hist[k]->Fill(SSBAR_values[k].at(i));
        for (int i = 0; i < (int)CHARM_values[k].size(); i++) CHARM_hist[k]->Fill(CHARM_values[k].at(i));
    }

    for (int k = 0; k < (int)variable_names.size(); k++) { // draw
        double CHG_int = CHG_hist[k]->Integral();
        double MIX_int = MIX_hist[k]->Integral();
        double UUBAR_int = UUBAR_hist[k]->Integral();
        double DDBAR_int = DDBAR_hist[k]->Integral();
        double SSBAR_int = SSBAR_hist[k]->Integral();
        double CHARM_int = CHARM_hist[k]->Integral();
        double SIGNAL_int = SIGNAL_hist[k]->Integral();

        double BKG_int = CHG_int + MIX_int + UUBAR_int + DDBAR_int + SSBAR_int + CHARM_int;

        CHG_hist[k]->Scale(1.0 / BKG_int, "width");
        MIX_hist[k]->Scale(1.0 / BKG_int, "width");
        UUBAR_hist[k]->Scale(1.0 / BKG_int, "width");
        DDBAR_hist[k]->Scale(1.0 / BKG_int, "width");
        SSBAR_hist[k]->Scale(1.0 / BKG_int, "width");
        CHARM_hist[k]->Scale(1.0 / BKG_int, "width");

        Stack[k]->Add(CHG_hist[k]);
        Stack[k]->Add(MIX_hist[k]);
        Stack[k]->Add(UUBAR_hist[k]);
        Stack[k]->Add(DDBAR_hist[k]);
        Stack[k]->Add(SSBAR_hist[k]);
        Stack[k]->Add(CHARM_hist[k]);

        SIGNAL_hist[k]->Scale(1.0 / SIGNAL_int, "width");
        SIGNAL_hist[k]->SetLineWidth(3);
        SIGNAL_hist[k]->SetLineColor(2);
        SIGNAL_hist[k]->SetFillStyle(0);

        TCanvas* c_temp = new TCanvas("c", "", 1200, 1200); c_temp->cd();

        gStyle->SetPalette(kGistEarth);

        Float_t ymax_1 = Stack[k]->GetMaximum();
        Float_t ymax_2 = SIGNAL_hist[k]->GetMaximum();
        double real_max = 0;
        if (ymax_1 > ymax_2) real_max = ymax_1;
        else real_max = ymax_2;

        Stack[k]->GetYaxis()->SetRangeUser(0, real_max * 1.1);

        gPad->BuildLegend(0.9, 0.9, 0.7, 0.7);
        //gPad->BuildLegend();
        legend.SetFillStyle(0);
        Stack[k]->Draw("pfc Hist"); SIGNAL_hist[k]->Draw("HistSAME");
        c_temp->SaveAs((variable_names+".png").c_str());

        delete c_temp;
    }
}
