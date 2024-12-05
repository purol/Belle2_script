#include "constants.h"
#include "correctors.h"
#include "base.h"

# define MCTYPE "MC15ri"

// my MC sample number
# define N_Kplus_nunubar 1000000.0
# define N_K0_nunubar 1000000.0
# define N_Kplusstar_nunubar 1000000.0
# define N_K0star_nunubar 1000000.0
# define N_Xsu_nonresonant_nunubar 5000000.0
# define N_Xsd_nonresonant_nunubar 5000000.0

// scale factor for each MC sample until LS1
# define Scale_Kplus (N_Kplus_nunubar_LS1/N_Kplus_nunubar)
# define Scale_Kplusstar (N_Kplusstar_nunubar_LS1/N_Kplusstar_nunubar)
# define Scale_Xsu_nonresonant (N_Xsu_nonresonant_nunubar_LS1/N_Xsu_nonresonant_nunubar)
# define Scale_K0 (N_K0_nunubar_LS1/N_K0_nunubar)
# define Scale_K0star (N_K0star_nunubar_LS1/N_K0star_nunubar)
# define Scale_Xsd_nonresonant (N_Xsd_nunubar_LS1/N_Xsd_nonresonant_nunubar)

Corrector corrector;
Corrector_Fragmentation corrector_Fragmentation;

/* ==========================================================================================*/

const char* decay_names[12] = {
"K",
"K1#pi w/o #pi^{0}",
"K1#pi w/ 1#pi^{0}",
"K2#pi w/o #pi^{0}",
"K2#pi w/ 1#pi^{0}",
"K3#pi w/o #pi^{0}",
"K3#pi w/ 1#pi^{0}",
"K4#pi w/o #pi^{0}",
"K4#pi w/ 1#pi^{0}",
"K w/ 2#pi^{0}",
"3K w/ atmost 1#pi^{0}",
"others"
};

typedef struct data{
    int __experiment__;
    int __run__;
    unsigned int __event__;
    int __candidate__;
    int __ncandidates__;

    double Tree[12];
    int Decay[N_decay];

    double invM;

    double MXs_c;
    double MXs_0;

} Data; 

class Loader {
public:
    enum DecayMode {
        K = 0,
        K1Pi_wo_Pi0,
        K1Pi_w_1Pi0,
        K2Pi_wo_Pi0,
        K2Pi_w_1Pi0,
        K3Pi_wo_Pi0,
        K3Pi_w_1Pi0,
        K4Pi_wo_Pi0,
        K4Pi_w_1Pi0,
        K_w_2pi0_w_atmost_4Pi,
        _3K_w_atmost_1Pi,
        others,
        MAX_NUM_DECAYMODE
    };
    enum DecayModeMC { // MC level
        Xsu2Kc_MC = 0,
        Xsu2Kcstar2KcPi0_MC,
        Xsu2Kcstar2K0Pic_MC,
        Xsu2KcPi0_MC,
        Xsu2K0Pic_MC,
        Xsu2KcPicPic_MC,
        Xsu2K0PicPi0_MC,
        Xsu2KcPicPicPi0_MC,
        Xsu2K0PicPicPic_MC,
        Xsu2KcPicPicPicPic_MC,
        Xsu2K0PicPicPicPi0_MC,
        Xsu2KcPi0Pi0_MC,
        Xsu2K0PicPi0Pi0_MC,
        Xsu2KcPicPicPi0Pi0_MC,
        Xsu2KcKcKc_MC,
        Xsu2KcKcK0Pic_MC,
        Xsu2KcKcKcPi0_MC,
        Xsd2K0_MC,
        Xsd2K0star2KcPic_MC,
        Xsd2K0star2K0Pi0_MC,
        Xsd2KcPic_MC,
        Xsd2K0Pi0_MC,
        Xsd2KcPicPi0_MC,
        Xsd2K0PicPic_MC,
        Xsd2KcPicPicPic_MC,
        Xsd2K0PicPicPi0_MC,
        Xsd2KcPicPicPicPi0_MC,
        Xsd2K0PicPicPicPic_MC,
        Xsd2K0Pi0Pi0_MC,
        Xsd2KcPicPi0Pi0_MC,
        Xsd2K0PicPicPi0Pi0_MC,
        Xsd2KcKcK0_MC,
        Xsd2KcKcKcPic_MC,
        Xsd2KcKcK0Pi0_MC,
        other,
        MAX_NUM_DECAYMODE_MC
    };
private:
    std::queue<Data> TotalData;

    std::vector<int> N_events;
    std::vector<int> N_candidates;
    std::vector<std::string> titles;
    int current_N_event;
    int current_N_candidate;

    std::vector<THStack*> THStacks;
    std::vector<TH1D*> TH1Ds_THStack[Loader::MAX_NUM_DECAYMODE];
    int current_THStack;

    std::vector<double> N_MC_modes[Loader::MAX_NUM_DECAYMODE_MC];

    bool TrueIfDecayModeMatch_MC(Data temp_data, Loader::DecayModeMC decaymodeMC);

public:
    Loader();
    void initialize();
    void GetData(TFile* input_file);
    bool event_info_is_valid();
    void DrawTHStack(const char* name, const char* title, int nbins, double x_low, double x_high, std::string filename = std::string(""), bool smart_mode = true);
    void PrintInformation(std::string title, std::string filename);
    double Mxs(Data data);
    bool AreTheyNeutrinosAndConj(double pdg1, double pdg2);
    bool AreTheyBmesonAndXs(double Bpdg, double Xspdg);
    Loader::DecayMode PrintDecayClassification(Data data);
    void End();
};

Loader::Loader() {
    std::queue<Data> empty;
    std::swap(TotalData, empty);
    current_N_event = 0;
    current_N_candidate = 0;
    current_THStack = 0;
}

void Loader::initialize() {
    std::queue<Data> empty;
    std::swap(TotalData, empty);
    current_N_event = 0;
    current_N_candidate = 0;
    current_THStack = 0;
}

void Loader::GetData(TFile* input_file) {
    TTree* tree_Xs = (TTree*)input_file->Get("Xs");

    Data temp = { 0 };

    // get event_info
    tree_Xs->SetBranchAddress("__experiment__", &temp.__experiment__);
    tree_Xs->SetBranchAddress("__run__", &temp.__run__);
    tree_Xs->SetBranchAddress("__event__", &temp.__event__);
    tree_Xs->SetBranchAddress("__candidate__", &temp.__candidate__);
    tree_Xs->SetBranchAddress("__ncandidates__", &temp.__ncandidates__);

    // get Upsilon_info
    tree_Xs->SetBranchAddress("daughter__bo0__cm__spnDaughters__bc", &temp.Tree[0]);
    tree_Xs->SetBranchAddress("daughter__bo0__cm__spPDG__bc", &temp.Tree[1]);
    tree_Xs->SetBranchAddress("daughter__bo0__cm__spdaughter__bo0__cmPDG__bc__bc", &temp.Tree[2]);
    tree_Xs->SetBranchAddress("daughter__bo0__cm__spdaughter__bo1__cmPDG__bc__bc", &temp.Tree[3]);
    tree_Xs->SetBranchAddress("daughter__bo0__cm__spdaughter__bo2__cmPDG__bc__bc", &temp.Tree[4]);
    tree_Xs->SetBranchAddress("daughter__bo0__cm__spdaughter__bo0__cmM__bc__bc", &temp.Tree[5]);
    tree_Xs->SetBranchAddress("daughter__bo1__cm__spnDaughters__bc", &temp.Tree[6]);
    tree_Xs->SetBranchAddress("daughter__bo1__cm__spPDG__bc", &temp.Tree[7]);
    tree_Xs->SetBranchAddress("daughter__bo1__cm__spdaughter__bo0__cmPDG__bc__bc", &temp.Tree[8]);
    tree_Xs->SetBranchAddress("daughter__bo1__cm__spdaughter__bo1__cmPDG__bc__bc", &temp.Tree[9]);
    tree_Xs->SetBranchAddress("daughter__bo1__cm__spdaughter__bo2__cmPDG__bc__bc", &temp.Tree[10]);
    tree_Xs->SetBranchAddress("daughter__bo1__cm__spdaughter__bo0__cmM__bc__bc", &temp.Tree[11]);

    tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKcharge_total__bc", &temp.Decay[0]);
    tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &temp.Decay[1]);
    tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &temp.Decay[2]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clcomb__bc", &temp.Decay[3]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch1__bc", &temp.Decay[4]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch2__bc", &temp.Decay[5]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch3__bc", &temp.Decay[6]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch4__bc", &temp.Decay[7]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch5__bc", &temp.Decay[8]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch6__bc", &temp.Decay[9]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch7__bc", &temp.Decay[10]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch8__bc", &temp.Decay[11]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch9__bc", &temp.Decay[12]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch10__bc", &temp.Decay[13]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch11__bc", &temp.Decay[14]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch12__bc", &temp.Decay[15]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch13__bc", &temp.Decay[16]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch14__bc", &temp.Decay[17]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clch15__bc", &temp.Decay[18]);
    tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKneutral_total__bc", &temp.Decay[19]);
    tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &temp.Decay[20]);
    tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &temp.Decay[21]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clcomb__bc", &temp.Decay[22]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch16__bc", &temp.Decay[23]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch17__bc", &temp.Decay[24]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch18__bc", &temp.Decay[25]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch19__bc", &temp.Decay[26]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch20__bc", &temp.Decay[27]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch21__bc", &temp.Decay[28]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch22__bc", &temp.Decay[29]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch23__bc", &temp.Decay[30]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch24__bc", &temp.Decay[31]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch25__bc", &temp.Decay[32]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch26__bc", &temp.Decay[33]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch27__bc", &temp.Decay[34]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch28__bc", &temp.Decay[35]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch29__bc", &temp.Decay[36]);
    tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clch30__bc", &temp.Decay[37]);

    tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &temp.invM);

    tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &temp.MXs_c);
    tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &temp.MXs_0);

    printf("%lld entries...\n", tree_Xs->GetEntries());
    for (unsigned int j = 0; j < tree_Xs->GetEntries(); j++) { // Fill
        tree_Xs->GetEntry(j);
        TotalData.push(temp);
    }
    input_file->Close();
}

bool Loader::event_info_is_valid() {
    std::queue<Data> temp_queue;
    while (!TotalData.empty()) {
        Data temp_data = TotalData.front();
        TotalData.pop();

        if ((temp_data.Decay[0] + temp_data.Decay[1] + temp_data.Decay[2] + temp_data.Decay[3] + temp_data.Decay[19] + temp_data.Decay[20] + temp_data.Decay[21] + temp_data.Decay[22]) > 1.5) { // error
            printf("Decay num > 1: this event is rejected!\n");
            continue;
            //            return false;
        }

        temp_queue.push(temp_data);
    }
    TotalData = temp_queue;
    return true;
}

void Loader::DrawTHStack(const char* name, const char* title, int nbins, double x_low, double x_high, std::string filename, bool smart_mode) {
    if (THStacks.size() == current_THStack) { // allocate new thstacks
        THStack* stack = new THStack(name, title);
        THStacks.push_back(stack);
        for (int i = 0; i < Loader::MAX_NUM_DECAYMODE; i++) {
            TH1D* hist = new TH1D(decay_names[i], title, nbins, x_low, x_high);
            TH1Ds_THStack[i].push_back(hist);
        }
    }
    else if (THStacks.size() > current_THStack) { // use what I have
    }
    else { // error
        printf("ERROR!\n");
        exit(1);
    }

    TH1D* temp_hist[Loader::MAX_NUM_DECAYMODE];
    for (int i = 0; i < Loader::MAX_NUM_DECAYMODE;i++) {
        temp_hist[i] = TH1Ds_THStack[i].at(current_THStack);
    }

    std::queue<Data> temp_queue = TotalData;
    while (!temp_queue.empty()) {
        Data temp_data = temp_queue.front();
        temp_queue.pop();

        if (temp_data.Decay[0] + temp_data.Decay[1] + temp_data.Decay[2] + temp_data.Decay[3] + temp_data.Decay[19] + temp_data.Decay[20] + temp_data.Decay[21] + temp_data.Decay[22] > 1.5) { // error
            printf("ERROR!\n");
            exit(1);
        }

        int decaymodeid = PrintDecayClassification(temp_data);

        if (smart_mode == false) temp_hist[decaymodeid]->Fill(Mxs(temp_data));
        else {
            if (filename.find("B2Knunu") != string::npos) {
                double correction_weight = corrector.GetCorrectionFactor(temp_data.invM * temp_data.invM, "Bplus");
                // double correction_weight = 1.0;
                temp_hist[decaymodeid]->Fill(Mxs(temp_data), Scale_Kplus * correction_weight);
            }
            else if (filename.find("B2Kstarnunu") != string::npos) temp_hist[decaymodeid]->Fill(Mxs(temp_data), Scale_Kplusstar);
            else if (filename.find("B2Xsnunu") != string::npos) {
                double Correction_Fragmentation = corrector_Fragmentation.GetCorrectionFactor(temp_data.Decay, Mxs(temp_data), Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, "MC15ri");
                temp_hist[decaymodeid]->Fill(Mxs(temp_data), Scale_Xsu_nonresonant * Correction_Fragmentation);
                //temp_hist[decaymodeid]->Fill(Mxs(temp_data), Scale_Xsu_nonresonant);
            }
            else if (filename.find("B02K0nunu") != string::npos) {
                double correction_weight = corrector.GetCorrectionFactor(temp_data.invM * temp_data.invM, "Bzero");
                // double correction_weight = 1.0;
                temp_hist[decaymodeid]->Fill(Mxs(temp_data), Scale_K0 * correction_weight);
            }
            else if (filename.find("B02Kstar0nunu") != string::npos) temp_hist[decaymodeid]->Fill(Mxs(temp_data), Scale_K0star);
            else if (filename.find("B02Xsnunu") != string::npos) {
                double Correction_Fragmentation = corrector_Fragmentation.GetCorrectionFactor(temp_data.Decay, Mxs(temp_data), Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, "MC15ri");
                temp_hist[decaymodeid]->Fill(Mxs(temp_data), Correction_Fragmentation * Scale_Xsd_nonresonant);
                //temp_hist[decaymodeid]->Fill(Mxs(temp_data), Scale_Xsd_nonresonant);
            }
            else { temp_hist[decaymodeid]->Fill(Mxs(temp_data)); }
        }
    }

    current_THStack++;
}

void Loader::PrintInformation(std::string title, std::string filename) {
    typedef struct labels {
        int __experiment__;
        int __run__;
        unsigned int __event__;
        int __ncandidates__;
    } Labels;
    std::vector<Labels> label_list;

    if (N_events.size() == current_N_event && N_candidates.size() == current_N_candidate && N_events.size() == N_candidates.size()) { // allocate new int
        N_events.push_back(0);
        N_candidates.push_back(0);
        titles.push_back(title);
        for (int i = 0; i < Loader::MAX_NUM_DECAYMODE_MC; i++) N_MC_modes[i].push_back(0);
    }
    else if (N_events.size() > current_N_event && N_candidates.size() > current_N_candidate && N_events.size() == N_candidates.size() && current_N_event == current_N_candidate) { // use what I have
    }
    else { // error
        printf("ERROR!\n");
        exit(1);
    }

    std::queue<Data> temp_queue = TotalData;
    while (!temp_queue.empty()) {
        Data temp = temp_queue.front();
        temp_queue.pop();

        bool overlap = false;
        for (unsigned int i = 0; i < label_list.size(); i++) {
            if (label_list.at(i).__experiment__ == temp.__experiment__ && label_list.at(i).__run__ == temp.__run__ && label_list.at(i).__event__ == temp.__event__ && label_list.at(i).__ncandidates__ == temp.__ncandidates__) { overlap = true; }
        }
        if (overlap == false) {
            N_events.at(current_N_event) = N_events.at(current_N_event) + 1;
            Labels temp_Labels;
            temp_Labels.__experiment__ = temp.__experiment__;
            temp_Labels.__run__ = temp.__run__;
            temp_Labels.__event__ = temp.__event__;
            temp_Labels.__ncandidates__ = temp.__ncandidates__;
            label_list.push_back(temp_Labels);

            // get detail MC decay information
            Loader::DecayModeMC decaymodeid_MC = Loader::MAX_NUM_DECAYMODE_MC;
            for (int i = 0; i < Loader::MAX_NUM_DECAYMODE_MC; i++) { // find MC decay mode
                if (TrueIfDecayModeMatch_MC(temp, static_cast<Loader::DecayModeMC>(i))) {
                    decaymodeid_MC = static_cast<Loader::DecayModeMC>(i);
                    if ((filename.find("B2Xsnunu") != string::npos) && (decaymodeid_MC == Loader::Xsu2Kc_MC)) {} // something wrong, try to find another decay mode
                    else if ((filename.find("B2Xsnunu") != string::npos) && (decaymodeid_MC == Loader::Xsu2Kcstar2KcPi0_MC || decaymodeid_MC == Loader::Xsu2Kcstar2K0Pic_MC)) {} // something wrong, try to find another decay mode
                    else if ((filename.find("B02Xsnunu") != string::npos) && (decaymodeid_MC == Loader::Xsd2K0_MC)) {} // something wrong, try to find another decay mode
                    else if ((filename.find("B02Xsnunu") != string::npos) && (decaymodeid_MC == Loader::Xsd2K0star2KcPic_MC || decaymodeid_MC == Loader::Xsd2K0star2K0Pi0_MC)) {} // something wrong, try to find another decay mode
                    else break;
                }
            }
            if (decaymodeid_MC == Loader::MAX_NUM_DECAYMODE_MC) {
                printf("ERROR! MC decay id cannot be found\n");
                exit(1);
            }

            // Number of event with MC decayID (scaled)
            double temp_N = -1;
            if (decaymodeid_MC == Loader::Xsu2Kc_MC) {
                double correction_weight = corrector.GetCorrectionFactor(temp.invM * temp.invM, "Bplus");
                temp_N = Scale_Kplus * correction_weight;
            }
            else if (decaymodeid_MC == Loader::Xsu2Kcstar2KcPi0_MC || decaymodeid_MC == Loader::Xsu2Kcstar2K0Pic_MC) temp_N = Scale_Kplusstar;
            else if (static_cast<int>(Xsu2KcPi0_MC) <= static_cast<int>(decaymodeid_MC) && static_cast<int>(decaymodeid_MC) <= static_cast<int>(Xsu2KcKcKcPi0_MC)) {
                //double correction_fragmentation = corrector_Fragmentation.GetCorrectionFactor(temp.Decay, Mxs(temp), Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, "MC15ri");
                double correction_fragmentation = corrector_Fragmentation.GetCorrectionFactor(temp.Decay, Mxs(temp), Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, "MC15ri");
                temp_N = Scale_Xsu_nonresonant * correction_fragmentation;
            }
            else if (decaymodeid_MC == Loader::Xsd2K0_MC) {
                double correction_weight = corrector.GetCorrectionFactor(temp.invM * temp.invM, "Bzero");
                temp_N = Scale_K0 * correction_weight;
            }
            else if (decaymodeid_MC == Loader::Xsd2K0star2KcPic_MC || decaymodeid_MC == Loader::Xsd2K0star2K0Pi0_MC) temp_N = Scale_K0star;
            else if (static_cast<int>(Xsd2KcPic_MC) <= static_cast<int>(decaymodeid_MC) && static_cast<int>(decaymodeid_MC) <= static_cast<int>(other)) {
                //double correction_fragmentation = corrector_Fragmentation.GetCorrectionFactor(temp.Decay, Mxs(temp), Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, "MC15ri");
                double correction_fragmentation = corrector_Fragmentation.GetCorrectionFactor(temp.Decay, Mxs(temp), Corrector_Fragmentation::SystType::Nominal, Corrector_Fragmentation::Sample::gamma, "MC15ri");
                temp_N = Scale_Xsd_nonresonant * correction_fragmentation;
            }
            else {
                printf("ERROR 265\n");
                exit(1);
            }
            N_MC_modes[decaymodeid_MC].at(current_N_event) = N_MC_modes[decaymodeid_MC].at(current_N_event) + temp_N;
        }
    }
    N_candidates.at(current_N_candidate) = N_candidates.at(current_N_candidate) + TotalData.size();

    current_N_event++;
    current_N_candidate++;
}

void Loader::End() {

    for (int i = 0; i < N_events.size();i++) {
        printf("%s\n", titles.at(i).c_str());
        printf("Number of event: %d\n", N_events.at(i));
        printf("Number of candidate: %d\n", N_candidates.at(i));
        for (int j = 0; j < Loader::MAX_NUM_DECAYMODE_MC; j++) printf("Number of event with MC decayID %d(scaled): %lf\n", j, N_MC_modes[j].at(i));
    }

    for (int i = 0; i < THStacks.size(); i++) {
        TCanvas* c_temp = new TCanvas("c", "", 1500, 1200); c_temp->cd();
        gStyle->SetPalette(kColorPrintableOnGrey);

        for (int j = 0; j < Loader::MAX_NUM_DECAYMODE; j++) {
            THStacks.at(i)->Add(TH1Ds_THStack[j].at(i));
            // TH1Ds_THStack[j].at(i)->Draw("Hist");
            // c_temp->SaveAs( (std::string(TH1Ds_THStack[j].at(i)->GetName()) + ".png").c_str() );
        }
        THStacks.at(i)->Draw("pfc Hist"); 
        c_temp->SaveAs((std::string(THStacks.at(i)->GetName()) + ".png").c_str());
        TLegend* legend = gPad->BuildLegend(0.9, 0.9, 0.65, 0.45);
        legend->SetFillStyle(0); legend->SetLineWidth(0);
        // c_temp->SetLogy();
        c_temp->SaveAs((std::string(THStacks.at(i)->GetName()) + "_legend.png").c_str());
        delete c_temp;
    }
}

double Loader::Mxs(Data data) {
    if (data.Tree[0] > 2.5 && data.Tree[0] < 3.5) { // Ndaughter of B == 3 
        if (AreTheyNeutrinosAndConj(data.Tree[2], data.Tree[3]) || AreTheyNeutrinosAndConj(data.Tree[2], data.Tree[4]) || AreTheyNeutrinosAndConj(data.Tree[4], data.Tree[3])) { // check B-> ? nu nubar
            if (AreTheyBmesonAndXs(data.Tree[1], data.Tree[2])) { // check B->Xs ? ?
                return data.Tree[5];
            }
        }
    }

    if (data.Tree[6] > 2.5 && data.Tree[6] < 3.5) { // Ndaughter of B == 3 
        if (AreTheyNeutrinosAndConj(data.Tree[8], data.Tree[9]) || AreTheyNeutrinosAndConj(data.Tree[8], data.Tree[10]) || AreTheyNeutrinosAndConj(data.Tree[10], data.Tree[9])) { // check B-> ? nu nubar
            if (AreTheyBmesonAndXs(data.Tree[7], data.Tree[8])) { // check B->Xs ? ?
                return data.Tree[11];
            }
        }
    }

    printf("ERROR! Cannot find B->Xs nu nubar decay!\n");
    printf("This may originate from the number of daughter! Condition wrt the num of daughter is removed and search again!\n");

    if (AreTheyNeutrinosAndConj(data.Tree[8], data.Tree[9]) || AreTheyNeutrinosAndConj(data.Tree[8], data.Tree[10]) || AreTheyNeutrinosAndConj(data.Tree[10], data.Tree[9])) { // check B-> ? nu nubar
        if (AreTheyBmesonAndXs(data.Tree[7], data.Tree[8])) { // check B->Xs ? ?
            return data.Tree[11];
        }
    }
    if (AreTheyNeutrinosAndConj(data.Tree[2], data.Tree[3]) || AreTheyNeutrinosAndConj(data.Tree[2], data.Tree[4]) || AreTheyNeutrinosAndConj(data.Tree[4], data.Tree[3])) { // check B-> ? nu nubar
        if (AreTheyBmesonAndXs(data.Tree[1], data.Tree[2])) { // check B->Xs ? ?
            return data.Tree[5];
        }
    }
    printf("ERROR! Again, cannot find B->Xs nu nubar decay!\n");
    exit(1);
    return -1;
}

bool Loader::AreTheyNeutrinosAndConj(double pdg1, double pdg2) {
    if (pdg1 * pdg2 > 0) return false;

    if (std::fabs(pdg1) > 11.5 && std::fabs(pdg1) < 12.5 && std::fabs(pdg2) > 11.5 && std::fabs(pdg2) < 12.5) return true; // nu_e
    else if (std::fabs(pdg1) > 13.5 && std::fabs(pdg1) < 14.5 && std::fabs(pdg2) > 13.5 && std::fabs(pdg2) < 14.5) return true; // nu_mu
    else if (std::fabs(pdg1) > 15.5 && std::fabs(pdg1) < 16.5 && std::fabs(pdg2) > 15.5 && std::fabs(pdg2) < 16.5) return true; // nu_tau

    return false;
}

bool Loader::AreTheyBmesonAndXs(double Bpdg, double Xspdg) {
    if (std::fabs(Xspdg) > 310.5 && std::fabs(Xspdg) < 311.5 && std::fabs(Bpdg) > 510.5 && std::fabs(Bpdg) < 511.5) return true; // B0/B0bar -> K0 nu nubar
    if (Bpdg * Xspdg < 0) return false; // it is not proper for B->K0 nu nubar

    if (std::fabs(Bpdg) > 520.5 && std::fabs(Bpdg) < 521.5) { // B+/B- -> Xsu/Xsubar
        if (std::fabs(Xspdg) > 30352.5 && std::fabs(Xspdg) < 30353.5) { // Xs is non-resonant
            return true;
        }
        else if (std::fabs(Xspdg) > 320.5 && std::fabs(Xspdg) < 321.5) { // Xs is K
            return true;
        }
        else if (std::fabs(Xspdg) > 322.5 && std::fabs(Xspdg) < 323.5) { // Xs is K*
            return true;
        }
    }

    if (std::fabs(Bpdg) > 510.5 && std::fabs(Bpdg) < 511.5) { // B0/B0bar -> Xsd/Xsdbar
        if (std::fabs(Xspdg) > 30342.5 && std::fabs(Xspdg) < 30343.5) { // Xs is non-resonant
            return true;
        }
        else if (std::fabs(Xspdg) > 310.5 && std::fabs(Xspdg) < 311.5) { // Xs K
            return true;
        }
        else if (std::fabs(Xspdg) > 312.5 && std::fabs(Xspdg) < 313.5) { // Xs K*
            return true;
        }
    }

    return false;
}

Loader::DecayMode Loader::PrintDecayClassification(Data data) {
    if (data.Decay[0] > 0 || data.Decay[4] > 0 || data.Decay[19] > 0 || data.Decay[23] > 0) return Loader::K; // B+->K+ / Xsu->K+ / B0->K0 / Xsd->K0
    else if (data.Decay[1] > 0 || data.Decay[6] > 0 || data.Decay[20] > 0 || data.Decay[24] > 0) return Loader::K1Pi_wo_Pi0; // B+->K*+->K0 pi+ / Xsu->K0 pi+ / B0->K*0->K+ pi- / Xsd->K+ pi-
    else if (data.Decay[2] > 0 || data.Decay[5] > 0 || data.Decay[21] > 0 || data.Decay[25] > 0) return Loader::K1Pi_w_1Pi0; // B+->K*+->K+ pi0 / Xsu->K+ pi0 / B0->K*0->K0 pi0 / Xsd->K0 pi0
    else if (data.Decay[7] > 0 || data.Decay[27] > 0) return Loader::K2Pi_wo_Pi0; // Xsu->K+ pi- pi+ / Xsd->K0 pi+ pi-
    else if (data.Decay[8] > 0 || data.Decay[26] > 0) return Loader::K2Pi_w_1Pi0; // Xsu->K0 pi+ K0 / Xsd->K+ pi- pi0
    else if (data.Decay[10] > 0 || data.Decay[28] > 0) return Loader::K3Pi_wo_Pi0; // Xsu->K0 pi+ pi- pi+ / Xsd->K+ pi- pi+ pi-
    else if (data.Decay[9] > 0 || data.Decay[29] > 0) return Loader::K3Pi_w_1Pi0; // Xsu->K+ pi- pi+ pi0 / Xsd->K0 pi+ pi- pi0
    else if (data.Decay[11] > 0 || data.Decay[31] > 0) return Loader::K4Pi_wo_Pi0; // Xsu->K+ pi- pi+ pi- pi+ / Xsd->K0 pi+ pi- pi+ pi-
    else if (data.Decay[12] > 0 || data.Decay[30] > 0) return Loader::K4Pi_w_1Pi0; // Xsu->K0 pi+ pi- pi+ pi0 / Xsd->K+ pi- pi+ pi- pi0
    else if (data.Decay[13] > 0 || data.Decay[14] > 0 || data.Decay[15] > 0 || data.Decay[32] > 0 || data.Decay[33] > 0 || data.Decay[34] > 0) return Loader::K_w_2pi0_w_atmost_4Pi;
    else if (data.Decay[16] > 0 || data.Decay[17] > 0 || data.Decay[18] > 0 || data.Decay[35] > 0 || data.Decay[36] > 0 || data.Decay[37] > 0) return Loader::_3K_w_atmost_1Pi;
    return Loader::others;
}

bool Loader::TrueIfDecayModeMatch_MC(Data temp_data, Loader::DecayModeMC decaymodeMC) {

    switch (decaymodeMC) {
    case Loader::Xsu2Kc_MC:
        if (temp_data.Decay[0] > 0) return true;
        return false;
        break;
    case Loader::Xsu2Kcstar2KcPi0_MC:
        if (temp_data.Decay[2] > 0) return true;
        return false;
        break;
    case Loader::Xsu2Kcstar2K0Pic_MC:
        if (temp_data.Decay[1] > 0) return true;
        return false;
        break;
    case Loader::Xsu2KcPi0_MC:
        if (temp_data.Decay[5] > 0) return true;
        return false;
        break;
    case Loader::Xsu2K0Pic_MC:
        if (temp_data.Decay[6] > 0) return true;
        return false;
        break;
    case Loader::Xsu2KcPicPic_MC:
        if (temp_data.Decay[7] > 0) return true;
        return false;
        break;
    case Loader::Xsu2K0PicPi0_MC:
        if (temp_data.Decay[8] > 0) return true;
        return false;
        break;
    case Loader::Xsu2KcPicPicPi0_MC:
        if (temp_data.Decay[9] > 0) return true;
        return false;
        break;
    case Loader::Xsu2K0PicPicPic_MC:
        if (temp_data.Decay[10] > 0) return true;
        return false;
        break;
    case Loader::Xsu2KcPicPicPicPic_MC:
        if (temp_data.Decay[11] > 0) return true;
        return false;
        break;
    case Loader::Xsu2K0PicPicPicPi0_MC:
        if (temp_data.Decay[12] > 0) return true;
        return false;
        break;
    case Loader::Xsu2KcPi0Pi0_MC:
        if (temp_data.Decay[13] > 0) return true;
        return false;
        break;
    case Loader::Xsu2K0PicPi0Pi0_MC:
        if (temp_data.Decay[14] > 0) return true;
        return false;
        break;
    case Loader::Xsu2KcPicPicPi0Pi0_MC:
        if (temp_data.Decay[15] > 0) return true;
        return false;
        break;
    case Loader::Xsu2KcKcKc_MC:
        if (temp_data.Decay[16] > 0) return true;
        return false;
        break;
    case Loader::Xsu2KcKcK0Pic_MC:
        if (temp_data.Decay[17] > 0) return true;
        return false;
        break;
    case Loader::Xsu2KcKcKcPi0_MC:
        if (temp_data.Decay[18] > 0) return true;
        return false;
        break;
    case Loader::Xsd2K0_MC:
        if (temp_data.Decay[19] > 0) return true;
        return false;
        break;
    case Loader::Xsd2K0star2KcPic_MC:
        if (temp_data.Decay[20] > 0) return true;
        return false;
        break;
    case Loader::Xsd2K0star2K0Pi0_MC:
        if (temp_data.Decay[21] > 0) return true;
        return false;
        break;
    case Loader::Xsd2KcPic_MC:
        if (temp_data.Decay[24] > 0) return true;
        return false;
        break;
    case Loader::Xsd2K0Pi0_MC:
        if (temp_data.Decay[25] > 0) return true;
        return false;
        break;
    case Loader::Xsd2KcPicPi0_MC:
        if (temp_data.Decay[26] > 0) return true;
        return false;
        break;
    case Loader::Xsd2K0PicPic_MC:
        if (temp_data.Decay[27] > 0) return true;
        return false;
        break;
    case Loader::Xsd2KcPicPicPic_MC:
        if (temp_data.Decay[28] > 0) return true;
        return false;
        break;
    case Loader::Xsd2K0PicPicPi0_MC:
        if (temp_data.Decay[29] > 0) return true;
        return false;
        break;
    case Loader::Xsd2KcPicPicPicPi0_MC:
        if (temp_data.Decay[30] > 0) return true;
        return false;
        break;
    case Loader::Xsd2K0PicPicPicPic_MC:
        if (temp_data.Decay[31] > 0) return true;
        return false;
        break;
    case Loader::Xsd2K0Pi0Pi0_MC:
        if (temp_data.Decay[32] > 0) return true;
        return false;
        break;
    case Loader::Xsd2KcPicPi0Pi0_MC:
        if (temp_data.Decay[33] > 0) return true;
        return false;
        break;
    case Loader::Xsd2K0PicPicPi0Pi0_MC:
        if (temp_data.Decay[34] > 0) return true;
        return false;
        break;
    case Loader::Xsd2KcKcK0_MC:
        if (temp_data.Decay[35] > 0) return true;
        return false;
        break;
    case Loader::Xsd2KcKcKcPic_MC:
        if (temp_data.Decay[36] > 0) return true;
        return false;
        break;
    case Loader::Xsd2KcKcK0Pi0_MC:
        if (temp_data.Decay[37] > 0) return true;
        return false;
        break;
    case Loader::other:
        return true;
        break;
    default:
        printf("ERROR! Input value of TrueIfDecayModeMatch_MC is not appropriate\n");
        exit(1);
        break;
    }

    printf("ERROR! Input value of TrueIfDecayModeMatch_MC is not appropriate\n");
    exit(1);
    return false;
}

void ReadDecayFiles_r_sp(){

    std::vector<string> names;
    const char* dirname = "/home/belle2/junewoo/storage_b2/MXs_distribution/small";

    load_files(dirname, &names);

    Loader loader;

    for(unsigned int i = 0; i<names.size(); i++){
        loader.initialize();

        TFile *input_file = new TFile( (dirname+std::string("/")+names.at(i)).c_str(),"read");
        printf("%s (%d/%zu)\n",("Read "+names.at(i) + "... ").c_str(), i, names.size());
        loader.GetData(input_file);
        if (loader.event_info_is_valid() == false) { printf("error!\n"); return; }

        loader.PrintInformation(std::string("========== inital =========="), names.at(i));
        loader.DrawTHStack("Mxs", ";M_{Xs}^{gen} [GeV]; arbitrary unit", 100, 0.45, 3.5, names.at(i), true);

    }
    loader.End();
}
