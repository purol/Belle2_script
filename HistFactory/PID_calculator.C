#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include <cmath>
#include <float.h>
#include <string>
#include <string.h>
#include "TFile.h"
#include "TTree.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "RooGenericPdf.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "TH1.h"
#include "TF1.h"
#include "RooDataHist.h"
#include <time.h>
#include <random>
using namespace RooFit;
using std::string;
using std::to_string;

# define MyEPSILON 0.000001

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

// uncertainty of BR from parametric reason
# define Sigma_BR_Kplus_nunubar 0.00000019
# define Sigma_BR_K0star_nunubar 0.00000050
# define Sigma_BR_K0_nunubar (Sigma_BR_Kplus_nunubar*TB0/TBp)
# define Sigma_BR_Kplusstar_nunubar (Sigma_BR_K0star_nunubar*TBp/TB0)
# define Sigma_BR_Xs_nunubar 0.000003

// https://confluence.desy.de/pages/viewpage.action?pageId=107054222
# define N_BpBp_1invab 565400000.0
# define N_B0B0_1invab 534600000.0

# define N_Kplus_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Kplus_nunubar)
# define N_Kplusstar_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Kplusstar_nunubar)
# define N_Xsu_nonresonant_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Xsu_nonresonant_nunubar)
# define N_K0_nunubar_1invab (2.0 * N_B0B0_1invab * BR_K0_nunubar)
# define N_K0star_nunubar_1invab (2.0 * N_B0B0_1invab * BR_K0star_nunubar)
# define N_Xsd_nunubar_1invab (2.0 * N_B0B0_1invab * BR_Xsd_nonresonant_nunubar)

# define Lumi_validation_MC (0.364436 - 0.002763) // ab-1

/*
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
*/

// SIGNAL MC sample number
# define N_Kplus_train 7039000.0
# define N_K0_train 7166624.0
# define N_Kplusstar_train 7039000.0
# define N_K0star_train 7166624.0
# define N_Xsu_nonresonant_train 35195000.0
# define N_Xsd_nonresonant_train 34940430.0
# define N_Kplus_test 2961000.0
# define N_K0_test 2833376.0
# define N_Kplusstar_test 2961000.0
# define N_K0star_test 2833376.0
# define N_Xsu_nonresonant_test 14805000.0
# define N_Xsd_nonresonant_test 15059570.0

// scale factor for SIGNAL MC sample (364.436 - 2.763 = 361.673/fb)
# define Scale_Kplus_train ((0.364436 - 0.002763) * N_Kplus_nunubar_1invab/N_Kplus_train)
# define Scale_Kplusstar_train ((0.364436 - 0.002763) * N_Kplusstar_nunubar_1invab/N_Kplusstar_train)
# define Scale_Xsu_nonresonant_train ((0.364436 - 0.002763) * N_Xsu_nonresonant_nunubar_1invab/N_Xsu_nonresonant_train)
# define Scale_K0_train ((0.364436 - 0.002763) * N_K0_nunubar_1invab/N_K0_train)
# define Scale_K0star_train ((0.364436 - 0.002763) * N_K0star_nunubar_1invab/N_K0star_train)
# define Scale_Xsd_nonresonant_train ((0.364436 - 0.002763) * N_Xsd_nunubar_1invab/N_Xsd_nonresonant_train)
# define Scale_Kplus_test ((0.364436 - 0.002763) * N_Kplus_nunubar_1invab/N_Kplus_test)
# define Scale_Kplusstar_test ((0.364436 - 0.002763) * N_Kplusstar_nunubar_1invab/N_Kplusstar_test)
# define Scale_Xsu_nonresonant_test ((0.364436 - 0.002763) * N_Xsu_nonresonant_nunubar_1invab/N_Xsu_nonresonant_test)
# define Scale_K0_test ((0.364436 - 0.002763) * N_K0_nunubar_1invab/N_K0_test)
# define Scale_K0star_test ((0.364436 - 0.002763) * N_K0star_nunubar_1invab/N_K0star_test)
# define Scale_Xsd_nonresonant_test ((0.364436 - 0.002763) * N_Xsd_nunubar_1invab/N_Xsd_nonresonant_test)

// BKG MC sample number (0.8/ab for BB, 1.0/ab for qq)
# define N_CHG_train 32042497.0
# define N_MIX_train 24693710.0
# define N_UUBAR_train 94447089.0
# define N_DDBAR_train 22664556.0
# define N_SSBAR_train 19244661.0
# define N_CHARM_train 107541168.0
# define N_CHG_test 48052238.0
# define N_MIX_test 37030486.0
# define N_UUBAR_test 141671998.0
# define N_DDBAR_test 34114182.0
# define N_SSBAR_test 28859338.0
# define N_CHARM_test 161280679.0

// scale factor for BKG MC sample (364.436 - 2.763 = 361.673/fb)
# define Scale_CHG_train (((0.364436 - 0.002763)/0.8)*(N_CHG_train + N_CHG_test)/ N_CHG_train)
# define Scale_MIX_train (((0.364436 - 0.002763)/0.8)*(N_MIX_train + N_MIX_test)/ N_MIX_train)
# define Scale_UUBAR_train (((0.364436 - 0.002763)/1.0)*(N_UUBAR_train + N_UUBAR_test)/ N_UUBAR_train)
# define Scale_DDBAR_train (((0.364436 - 0.002763)/1.0)*(N_DDBAR_train + N_DDBAR_test)/ N_DDBAR_train)
# define Scale_SSBAR_train (((0.364436 - 0.002763)/1.0)*(N_SSBAR_train + N_SSBAR_test)/ N_SSBAR_train)
# define Scale_CHARM_train (((0.364436 - 0.002763)/1.0)*(N_CHARM_train + N_CHARM_test)/ N_CHARM_train)
# define Scale_CHG_test (((0.364436 - 0.002763)/0.8)*(N_CHG_train + N_CHG_test)/ N_CHG_test)
# define Scale_MIX_test (((0.364436 - 0.002763)/0.8)*(N_MIX_train + N_MIX_test)/ N_MIX_test)
# define Scale_UUBAR_test (((0.364436 - 0.002763)/1.0)*(N_UUBAR_train + N_UUBAR_test)/ N_UUBAR_test)
# define Scale_DDBAR_test (((0.364436 - 0.002763)/1.0)*(N_DDBAR_train + N_DDBAR_test)/ N_DDBAR_test)
# define Scale_SSBAR_test (((0.364436 - 0.002763)/1.0)*(N_SSBAR_train + N_SSBAR_test)/ N_SSBAR_test)
# define Scale_CHARM_test (((0.364436 - 0.002763)/1.0)*(N_CHARM_train + N_CHARM_test)/ N_CHARM_test)

// scale facto for systematic MC sample
# define N_K0star_nunubar_syst 10000000.0
# define N_Kplusstar_nunubar_syst 10000000.0
# define N_Xsu_nonresonant_nunubar_syst 10000000.0
# define N_Xsd_nonresonant_nunubar_syst 10000000.0

// scale factor for each systematic MC sample (364.436 - 2.763 = 361.673/fb)
# define Scale_Kplusstar_syst ((0.364436 - 0.002763) * N_Kplusstar_nunubar_1invab/N_Kplusstar_nunubar_syst)
# define Scale_Xsu_nonresonant_syst ((0.364436 - 0.002763) * N_Xsu_nonresonant_nunubar_1invab/N_Xsu_nonresonant_nunubar_syst)
# define Scale_K0star_syst ((0.364436 - 0.002763) * N_K0star_nunubar_1invab/N_K0star_nunubar_syst)
# define Scale_Xsd_nonresonant_syst ((0.364436 - 0.002763) * N_Xsd_nunubar_1invab/N_Xsd_nonresonant_nunubar_syst)

# define KS0_rel_uncertainty 0.6 // %/cm
# define track_rel_uncertainty 0.69 // %
# define pi0_correction 0.932
# define pi0_rel_uncertainty ((0.0369 / 0.932) * 100.0) // %
# define Kaon_PID_max_uncertainty 0.1 // not percentage. relative uncertainty
// https://indico.belle2.org/event/6872/contributions/37447/attachments/17127/25504/FEIperformance_B2GM.pdf
# define FEI_cal_Bc 0.679
# define FEI_cal_Bc_uncertainty (0.017/FEI_cal_Bc) // not percentage. relative uncertainty
# define FEI_cal_B0 0.713
# define FEI_cal_B0_uncertainty (0.019/FEI_cal_B0) // not percentage. relative uncertainty

# define Xsu_frag_decay1 531.4
# define Xsu_frag_decay2 1062.1
# define Xsu_frag_decay3 673.9
# define Xsu_frag_decay4 825.3
# define Xsu_frag_decay5 473.3
# define Xsu_frag_decay6 119.1
# define Xsu_frag_decay7 24.6
# define Xsu_frag_decay8 148.2
# define Xsu_frag_decay9 127.6
# define Xsu_frag_decay10 102.9
# define Xsu_frag_decay11 91.8
# define Xsu_frag_decay12 54.8
# define Xsu_frag_decay13 21.7
# define Xsu_frag_decay14 11.0
# define Xsu_frag_decay30 858.4

# define Xsd_frag_decay16 1072.6
# define Xsd_frag_decay17 536.2
# define Xsd_frag_decay18 835.5
# define Xsd_frag_decay19 677.4
# define Xsd_frag_decay20 120.5
# define Xsd_frag_decay21 475.9
# define Xsd_frag_decay22 149.8
# define Xsd_frag_decay23 25.0
# define Xsd_frag_decay24 129.2
# define Xsd_frag_decay25 104.4
# define Xsd_frag_decay26 92.4
# define Xsd_frag_decay27 54.9
# define Xsd_frag_decay28 14.5
# define Xsd_frag_decay29 10.9
# define Xsd_frag_decay30 871.0

std::random_device rd;
std::default_random_engine generator(rd());

# define N_PID_syst 73
double PID_correction[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
double PID_correction_stat_uncer[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
double PID_correction_sys_uncer[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
double PID_correction_uncer[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
double PID_correction_fluctuated[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss

enum DecayMode { // reco level
    B2Kc = 0,
    B2KcPi0,
    B2Ks0Pic,
    B2KcPicPic,
    B2Ks0PicPi0,
    B2KcPicPicPi0,
    B2Ks0PicPicPic,
    B2KcPicPicPicPic,
    B2Ks0PicPicPicPi0,
    B2KcPi0Pi0,
    B2Ks0PicPi0Pi0,
    B2KcPicPicPi0Pi0,
    B2KcKcKc,
    B2KcKcKs0Pic,
    B2KcKcKcPi0,
    B02Ks0,
    B02KcPic,
    B02Ks0Pi0,
    B02KcPicPi0,
    B02Ks0PicPic,
    B02KcPicPicPic,
    B02Ks0PicPicPi0,
    B02KcPicPicPicPi0,
    B02Ks0PicPicPicPic,
    B02Ks0Pi0Pi0,
    B02KcPicPi0Pi0,
    B02Ks0PicPicPi0Pi0,
    B02KcKcKs0,
    B02KcKcKcPic,
    B02KcKcKs0Pi0,
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

# define N_decay 38 // five decay mode + others

//# define Htransition_Xsu_change (-0.142)
//# define Htransition_Xsd_change (-0.099)
//# define Ltransition_Xsu_change (-0.002)
//# define Ltransition_Xsd_change (0.067)
//# define Hmb_Xsu_change (-0.067)
//# define Hmb_Xsd_change (-0.058)
//# define Lmb_Xsu_change (-0.062)
//# define Lmb_Xsd_change (-0.046)
//# define Hpf_Xsu_change (-0.120)
//# define Hpf_Xsd_change (-0.073)
//# define Lpf_Xsu_change (-0.014)
//# define Lpf_Xsd_change (0.019)

//# define FEI_cal_Bc 1.0
//# define FEI_cal_Bc_uncertainty 0.0
//# define FEI_cal_B0 1.0
//# define FEI_cal_B0_uncertainty 0.0

//# define Htransition_Xsu_change 0.0
//# define Htransition_Xsd_change 0.0
//# define Ltransition_Xsu_change 0.0
//# define Ltransition_Xsd_change 0.0
//# define Hmb_Xsu_change 0.0
//# define Hmb_Xsd_change 0.0
//# define Lmb_Xsu_change 0.0
//# define Lmb_Xsd_change 0.0
//# define Hpf_Xsu_change 0.0
//# define Hpf_Xsd_change 0.0
//# define Lpf_Xsu_change 0.0
//# define Lpf_Xsd_change 0.0

# define RarityBins 15
# define NToys 500

/* ====================================== */
// Rarity module
bool sorting(pair<double, double> a, pair<double, double> b) {

    return a.first < b.first;

}

class CDF {
private:
    std::vector<double> m_values;
    std::vector<double> m_weights;
    double m_weights_sum;

    std::vector<double> m_signal_output;
    std::vector<double> m_signal_weight;
    std::vector<double> m_background_output;
    std::vector<double> m_background_weight;

    const double m_max = 1.0;
    const double m_min = 0.92;

    template <typename A, typename B>
    void zip(const std::vector<A>& a, const std::vector<B>& b, std::vector<std::pair<A, B>>& zipped)
    {
        for (size_t i = 0; i < a.size(); ++i)
        {
            zipped.push_back(std::make_pair(a[i], b[i]));
        }
    }

    template <typename A, typename B>
    void unzip(const std::vector<std::pair<A, B>>& zipped, std::vector<A>& a, std::vector<B>& b)
    {
        for (size_t i = 0; i < a.size(); i++)
        {
            a[i] = zipped[i].first;
            b[i] = zipped[i].second;
        }
    }

    void m_load_files(const char* dirname, std::vector<std::string>* names) {
        TSystemDirectory dir(dirname, dirname);
        TList* files = dir.GetListOfFiles();
        if (files) {
            TSystemFile* file;
            TString fname;
            TIter next(files);
            while ((file = (TSystemFile*)next())) {
                fname = file->GetName();
                if (!file->IsDirectory() && fname.EndsWith(".root")) {
                    names->push_back(fname.Data());
                }
            }
        }
    }


public:
    void init(std::vector<double> values, std::vector<double> weights) {
        // Zip the vectors together
        std::vector<std::pair<double, double>> zipped;
        zip(values, weights, zipped);

        sort(zipped.begin(), zipped.end(), sorting);

        // Write the sorted pairs back to the vectors
        unzip(zipped, m_values, m_weights);
        m_weights_sum = std::accumulate(m_weights.begin(), m_weights.end(), 0.0);
    }
    void initbypath(const char* dirname, double weight) {
        std::vector<string> names;

       m_load_files(dirname, &names);

        for (unsigned int i = 0; i < names.size(); i++) {
            float MVA_BB = 0;

            TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");

            TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
            tree_upsilon->SetBranchAddress("MVA_BB", &MVA_BB);

            for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
                tree_upsilon->GetEntry(j);

                m_values.push_back(MVA_BB);
                m_weights.push_back(weight);

            }
            input_file->Close();
        }

        // Zip the vectors together
        std::vector<std::pair<double, double>> zipped;
        zip(m_values, m_weights, zipped);

        sort(zipped.begin(), zipped.end(), sorting);

        // Write the sorted pairs back to the vectors
        unzip(zipped, m_values, m_weights);
        m_weights_sum = std::accumulate(m_weights.begin(), m_weights.end(), 0.0);
    }
    double GetCDFvalue(float value) {
        const double max = 1.0;
        const double min = 0.85;

        auto larger_or_equal = std::lower_bound(m_values.begin(), m_values.end(), value);
        auto larger = std::upper_bound(m_values.begin(), m_values.end(), value);

        int until_index = -1;
        if (larger_or_equal - m_values.begin() == larger - m_values.begin()) { // input value is not same with one of value which is saved
            until_index = larger_or_equal - m_values.begin() - 1;
        }
        else { // input value is same with one of value which is saved
            until_index = larger_or_equal - m_values.begin();
            double weights_partial_sum = std::accumulate(m_weights.begin(), until_index + m_weights.begin(), 0.0);
            return weights_partial_sum / m_weights_sum;
        }

        if (value < m_values.at(0)) {
            int temp_index = std::upper_bound(m_values.begin(), m_values.end(), m_values.at(0)) - m_values.begin() - 1;
            double weights_partial_sum = (value - m_min) * std::accumulate(m_weights.begin(), temp_index + m_weights.begin(), 0.0) / (m_values.at(0) - m_min);
            return weights_partial_sum / m_weights_sum;
        }
        else if (value >= *(m_values.end() - 1)) {
            return 1.0;
        }

        double weights_partial_sum_back = std::accumulate(m_weights.begin(), until_index + m_weights.begin(), 0.0);
        double weights_partial_sum_front = std::accumulate(m_weights.begin(), until_index + m_weights.begin() + 1, 0.0);

        int index_back = until_index;
        int index_front = until_index + 1;

        double weights_partial_sum = weights_partial_sum_back + (value - m_values.at(index_back)) * (weights_partial_sum_front - weights_partial_sum_back) / (m_values.at(index_front) - m_values.at(index_back));

        return weights_partial_sum / m_weights_sum;

    }
    void print() {
        for (unsigned int i = 0; i < m_values.size(); i++) {
            printf("[%lf, %lf] ", m_values.at(i), m_weights.at(i));
        }
        printf("\n");
    }
    void Evaluate(const char* dirname, double weight, bool IsSignal) {
        std::vector<string> names;

        m_load_files(dirname, &names);

        for (unsigned int i = 0; i < names.size(); i++) {
            float MVA_BB = 0;

            TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");

            TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
            tree_upsilon->SetBranchAddress("MVA_BB", &MVA_BB);

            for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
                tree_upsilon->GetEntry(j);

                if (IsSignal) {
                    m_signal_output.push_back(GetCDFvalue(MVA_BB));
                    m_signal_weight.push_back(weight);
                }
                else {
                    m_background_output.push_back(GetCDFvalue(MVA_BB));
                    m_background_weight.push_back(weight);
                }
            }
            input_file->Close();

        }
    }
    void PrintPNG() {
        TH1D* temp_signal_hist = new TH1D("CDF_signal", ";CDF;arbitrary unit", 100, 0.0, 1.0);
        TH1D* temp_background_hist = new TH1D("CDF_background", ";CDF;arbitrary unit", 100, 0.0, 1.0);

        for (unsigned int i = 0; i < m_signal_output.size(); i++) temp_signal_hist->Fill(m_signal_output.at(i), m_signal_weight.at(i));
        for (unsigned int i = 0; i < m_background_output.size(); i++) temp_background_hist->Fill(m_background_output.at(i), m_background_weight.at(i));

        temp_signal_hist->Scale(1.0 / temp_signal_hist->Integral(), "width");
        temp_background_hist->Scale(1.0 / temp_background_hist->Integral(), "width");

        TCanvas* c_temp = new TCanvas("c", "", 800, 800);

        temp_signal_hist->Draw("hist"); c_temp->SaveAs("signal_cdf.png");
        temp_background_hist->Draw("hist"); c_temp->SaveAs("background_cdf.png");

        delete c_temp;
    }
};

CDF cdf;

/* ====================================== */

// global variables to calculate uncertainties
std::vector<double> Ns;
std::vector<int> ntracks;
std::vector<int> npi0s;
std::vector<int> dmID;
std::vector<double> KS0_3D_distance;
std::vector<double> KaonID_correction;
std::vector<double> KaonID_rel_up;
std::vector<double> KaonID_rel_dn;

void load_files(const char* dirname, std::vector<string>* names);

void load_files(const char* dirname, std::vector<string>* names) {
    TSystemDirectory dir(dirname, dirname);
    TList* files = dir.GetListOfFiles();
    if (files) {
        TSystemFile* file;
        TString fname;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(".root")) {
                names->push_back(fname.Data());
            }
        }
    }
}

double GetEvtNum(const char* dirname, double weight_var = 1.0) {

    double Total_Evt_Num = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        Total_Evt_Num = Total_Evt_Num + weight_var * tree_upsilon->GetEntries();

        input_file->Close();

    }

    return Total_Evt_Num;

}

double GetNpi0(double Upsilon_ID, double Bsig_ID) {

    if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) return 0.0; // B2Kc
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) return 1.0; // B2KcPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) return 0.0; // B2Ks0Pic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) return 0.0; // B2KcPicPic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) return 1.0; // B2Ks0PicPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) return 1.0; // B2KcPicPicPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) return 0.0; // B2Ks0PicPicPic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) return 0.0; // B2KcPicPicPicPic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) return 1.0; // B2Ks0PicPicPicPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) return 2.0; // B2KcPi0Pi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) return 2.0; // B2Ks0PicPi0Pi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) return 2.0; // B2KcPicPicPi0Pi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) return 0.0; // B2KcKcKc
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) return 0.0; // B2KcKcKs0Pic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) return 1.0; // B2KcKcKcPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) return 0.0; // B02Ks0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) return 0.0; // B02KcPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) return 1.0; // B02Ks0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) return 1.0; // B02KcPicPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) return 0.0; // B02Ks0PicPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) return 0.0; // B02KcPicPicPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) return 1.0; // B02Ks0PicPicPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) return 1.0; // B02KcPicPicPicPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) return 0.0; // B02Ks0PicPicPicPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) return 2.0; // B02Ks0Pi0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) return 2.0; // B02KcPicPi0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) return 2.0; // B02Ks0PicPicPi0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) return 0.0; // B02KcKcKs0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) return 0.0; // B02KcKcKcPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) return 1.0; // B02KcKcKs0Pi0
    else {
        printf("[ERROR] unexpected decay ID\n");
        exit(1);
    }

    return -1.0;
}

double GetNtracks(double Upsilon_ID, double Bsig_ID) {

    if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) return 1.0; // B2Kc
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) return 1.0; // B2KcPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) return 1.0; // B2Ks0Pic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) return 3.0; // B2KcPicPic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) return 1.0; // B2Ks0PicPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) return 3.0; // B2KcPicPicPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) return 3.0; // B2Ks0PicPicPic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) return 5.0; // B2KcPicPicPicPic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) return 3.0; // B2Ks0PicPicPicPi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) return 1.0; // B2KcPi0Pi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) return 1.0; // B2Ks0PicPi0Pi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) return 3.0; // B2KcPicPicPi0Pi0
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) return 3.0; // B2KcKcKc
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) return 3.0; // B2KcKcKs0Pic
    else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) return 3.0; // B2KcKcKcPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) return 0.0; // B02Ks0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) return 2.0; // B02KcPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) return 0.0; // B02Ks0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) return 2.0; // B02KcPicPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) return 2.0; // B02Ks0PicPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) return 4.0; // B02KcPicPicPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) return 2.0; // B02Ks0PicPicPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) return 4.0; // B02KcPicPicPicPi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) return 4.0; // B02Ks0PicPicPicPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) return 0.0; // B02Ks0Pi0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) return 2.0; // B02KcPicPi0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) return 2.0; // B02Ks0PicPicPi0Pi0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) return 2.0; // B02KcKcKs0
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) return 4.0; // B02KcKcKcPic
    else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) return 2.0; // B02KcKcKs0Pi0
    else {
        printf("[ERROR] unexpected decay ID\n");
        exit(1);
    }

    return -1.0;
}

void GetNominalNevt(const char* dirname, TH1D* hist, const char* type, const char* sample, double Nevt_nominal[RarityBins * 7], double weight_var = 1.0) { // get nominal PDF with appropriate correction
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    if (strcmp(sample, "CHG") == 0) {}
    else if (strcmp(sample, "MIX") == 0) {}
    else if (strcmp(sample, "UUBAR") == 0) {}
    else if (strcmp(sample, "DDBAR") == 0) {}
    else if (strcmp(sample, "SSBAR") == 0) {}
    else if (strcmp(sample, "CHARM") == 0) {}
    else if (strcmp(sample, "SIGNAL") == 0) {}
    else {
        printf("[ERROR] unexpected sample name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss

    std::vector<string> names;
    load_files(dirname, &names);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKtruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[0][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKmisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[1][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npitruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[2][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npimisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[3][i_PID]);
        }

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

            double Correction_pi0 = std::pow(pi0_correction, Npi0);
            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = FEI_cal_Bc;
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = FEI_cal_B0;
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(PID_correction[0][i_PID], temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(PID_correction[1][i_PID], temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(PID_correction[2][i_PID], temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(PID_correction[3][i_PID], temp_N_bin_PID[3][i_PID]); // mis PID
            }

            double total_weight = weight_var * Correction_pi0 * Correction_FEI * Correction_KID * Correction_PID;

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

        printf("%s has %lf events (with correction)\n", dirname, Nevt);

    }

    int ArrayBinID = -1;
    
    if (strcmp(sample, "CHG") == 0) ArrayBinID = 0;
    else if (strcmp(sample, "MIX") == 0) ArrayBinID = 1;
    else if (strcmp(sample, "UUBAR") == 0) ArrayBinID = 2;
    else if (strcmp(sample, "DDBAR") == 0) ArrayBinID = 3;
    else if (strcmp(sample, "SSBAR") == 0) ArrayBinID = 4;
    else if (strcmp(sample, "CHARM") == 0) ArrayBinID = 5;
    else if (strcmp(sample, "SIGNAL") == 0) ArrayBinID = 6;

    for (int i = 0; i < RarityBins; i++) {
        Nevt_nominal[ArrayBinID * RarityBins + i] = Nevt_nominal[ArrayBinID * RarityBins + i] + hist->GetBinContent(i + 1);
    }

    return;
}

void GetFlucNevt(const char* dirname, TH1D* hist, const char* type, const char* sample, double Nevt_fluc[NToys][RarityBins * 7], int ToyNum, bool IsItKID, double weight_var = 1.0) { // get nominal PDF with appropriate correction
    if (strcmp(type, "Bplus") == 0) {}
    else if (strcmp(type, "Bzero") == 0) {}
    else if (strcmp(type, "Continuum") == 0) {}
    else {
        printf("[ERROR] unexpected type name\n");
        exit(1);
    }

    if (strcmp(sample, "CHG") == 0) {}
    else if (strcmp(sample, "MIX") == 0) {}
    else if (strcmp(sample, "UUBAR") == 0) {}
    else if (strcmp(sample, "DDBAR") == 0) {}
    else if (strcmp(sample, "SSBAR") == 0) {}
    else if (strcmp(sample, "CHARM") == 0) {}
    else if (strcmp(sample, "SIGNAL") == 0) {}
    else {
        printf("[ERROR] unexpected sample name\n");
        exit(1);
    }

    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss

    std::vector<string> names;
    load_files(dirname, &names);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("MVA_BB", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKtruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[0][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKmisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[1][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npitruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[2][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npimisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[3][i_PID]);
        }

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            double Npi0 = GetNpi0(Upsilon_ID, Bsig_ID);

            double Correction_pi0 = std::pow(pi0_correction, Npi0);
            double Correction_FEI = 1.0;
            if (strcmp(type, "Bplus") == 0) Correction_FEI = FEI_cal_Bc;
            else if (strcmp(type, "Bzero") == 0) Correction_FEI = FEI_cal_B0;
            else if (strcmp(type, "Continuum") == 0) Correction_FEI = 1.0;
            double Correction_KID = 1;
            double Correction_PID = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {

                if (IsItKID) {
                    Correction_KID = Correction_KID * std::pow(PID_correction[0][i_PID] * PID_correction_fluctuated[0][i_PID], temp_N_bin_PID[0][i_PID]); // true KID
                    Correction_KID = Correction_KID * std::pow(PID_correction[1][i_PID] * PID_correction_fluctuated[1][i_PID], temp_N_bin_PID[1][i_PID]); // mis KID
                    Correction_PID = Correction_PID * std::pow(PID_correction[2][i_PID], temp_N_bin_PID[2][i_PID]); // true PID
                    Correction_PID = Correction_PID * std::pow(PID_correction[3][i_PID], temp_N_bin_PID[3][i_PID]); // mis PID
                }
                else {
                    Correction_KID = Correction_KID * std::pow(PID_correction[0][i_PID], temp_N_bin_PID[0][i_PID]); // true KID
                    Correction_KID = Correction_KID * std::pow(PID_correction[1][i_PID], temp_N_bin_PID[1][i_PID]); // mis KID
                    Correction_PID = Correction_PID * std::pow(PID_correction[2][i_PID] * PID_correction_fluctuated[2][i_PID], temp_N_bin_PID[2][i_PID]); // true PID
                    Correction_PID = Correction_PID * std::pow(PID_correction[3][i_PID] * PID_correction_fluctuated[3][i_PID], temp_N_bin_PID[3][i_PID]); // mis PID
                }

            }

            double total_weight = weight_var * Correction_pi0 * Correction_FEI * Correction_KID * Correction_PID;

            Nevt = Nevt + total_weight;

            hist->Fill(MVA_var, total_weight);
        }
        input_file->Close();

        printf("%s has %lf events (with correction)\n", dirname, Nevt);

    }

    int ArrayBinID = -1;

    if (strcmp(sample, "CHG") == 0) ArrayBinID = 0;
    else if (strcmp(sample, "MIX") == 0) ArrayBinID = 1;
    else if (strcmp(sample, "UUBAR") == 0) ArrayBinID = 2;
    else if (strcmp(sample, "DDBAR") == 0) ArrayBinID = 3;
    else if (strcmp(sample, "SSBAR") == 0) ArrayBinID = 4;
    else if (strcmp(sample, "CHARM") == 0) ArrayBinID = 5;
    else if (strcmp(sample, "SIGNAL") == 0) ArrayBinID = 6;

    for (int i = 0; i < RarityBins; i++) {
        Nevt_fluc[ToyNum][ArrayBinID * RarityBins + i] = Nevt_fluc[ToyNum][ArrayBinID * RarityBins + i] + hist->GetBinContent(i + 1);
    }

    return;
}

void ReadPIDFile() {
    const char* KID_true_file = "KaonEff.csv";
    const char* KID_mis_file = "Kaonmis.csv";
    const char* PID_true_file = "PionEff.csv";
    const char* PID_mis_file = "Pionmis.csv";

    FILE* fp_KID_true = fopen(KID_true_file, "r");
    FILE* fp_KID_mis = fopen(KID_mis_file, "r");
    FILE* fp_PID_true = fopen(PID_true_file, "r");
    FILE* fp_PID_mis = fopen(PID_mis_file, "r");

    fscanf(fp_KID_true, "p_min,p_max,cosTheta_min,cosTheta_max,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,threshold,variable\n");
    fscanf(fp_KID_mis, "p_min,p_max,cosTheta_min,cosTheta_max,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,threshold,variable\n");
    fscanf(fp_PID_true, "p_min,p_max,cosTheta_min,cosTheta_max,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,threshold,variable\n");
    fscanf(fp_PID_mis, "p_min,p_max,cosTheta_min,cosTheta_max,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,threshold,variable\n");

    double temp_p_min;
    double temp_p_max;
    double temp_cosTheta_min;
    double temp_cosTheta_max;
    double temp_data_MC_ratio;
    double temp_data_MC_uncertainty_stat_up;
    double temp_data_MC_uncertainty_stat_dn;
    double temp_data_MC_uncertainty_sys_up;
    double temp_data_MC_uncertainty_sys_dn;
    double temp_data_efficiency;
    double temp_data_uncertainty_stat_up;
    double temp_data_uncertainty_stat_dn;
    double temp_data_uncertainty_sys_up;
    double temp_data_uncertainty_sys_dn;
    double temp_MC_efficiency;
    double temp_MC_uncertainty_stat_up;
    double temp_MC_uncertainty_stat_dn;
    double temp_MC_uncertainty_sys_up;
    double temp_MC_uncertainty_sys_dn;
    double temp_threshold;

    for (int i = 0; i < N_PID_syst - 1; i++) {
        fscanf(fp_KID_true, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,kaonID\n", &temp_p_min, &temp_p_max, &temp_cosTheta_min, &temp_cosTheta_max, &PID_correction[0][i], &PID_correction_stat_uncer[0][i], &temp_data_MC_uncertainty_stat_dn, &PID_correction_sys_uncer[0][i], &temp_data_MC_uncertainty_sys_dn, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_threshold);
        fscanf(fp_KID_mis, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,kaonID\n", &temp_p_min, &temp_p_max, &temp_cosTheta_min, &temp_cosTheta_max, &PID_correction[1][i], &PID_correction_stat_uncer[1][i], &temp_data_MC_uncertainty_stat_dn, &PID_correction_sys_uncer[1][i], &temp_data_MC_uncertainty_sys_dn, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_threshold);
        fscanf(fp_PID_true, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,pionID\n", &temp_p_min, &temp_p_max, &temp_cosTheta_min, &temp_cosTheta_max, &PID_correction[2][i], &PID_correction_stat_uncer[2][i], &temp_data_MC_uncertainty_stat_dn, &PID_correction_sys_uncer[2][i], &temp_data_MC_uncertainty_sys_dn, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_threshold);
        fscanf(fp_PID_mis, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,pionID\n", &temp_p_min, &temp_p_max, &temp_cosTheta_min, &temp_cosTheta_max, &PID_correction[3][i], &PID_correction_stat_uncer[3][i], &temp_data_MC_uncertainty_stat_dn, &PID_correction_sys_uncer[3][i], &temp_data_MC_uncertainty_sys_dn, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_threshold);
 
        if ( (std::abs(PID_correction[0][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer[0][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer[0][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction[0][i]) > 10000.0) ) {
            PID_correction[0][i] = 1.0;
            PID_correction_stat_uncer[0][i] = 0.0;
            PID_correction_sys_uncer[0][i] = 0.0;
        }
        if ( (std::abs(PID_correction[1][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer[1][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer[1][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction[1][i]) > 10000.0) ) {
            PID_correction[1][i] = 1.0;
            PID_correction_stat_uncer[1][i] = 0.0;
            PID_correction_sys_uncer[1][i] = 0.0;
        }
        if ( (std::abs(PID_correction[2][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer[2][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer[2][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction[2][i]) > 10000.0)) {
            PID_correction[2][i] = 1.0;
            PID_correction_stat_uncer[2][i] = 0.0;
            PID_correction_sys_uncer[2][i] = 0.0;
        }
        if ( (std::abs(PID_correction[3][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer[3][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer[3][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction[3][i]) > 10000.0)) {
            PID_correction[3][i] = 1.0;
            PID_correction_stat_uncer[3][i] = 0.0;
            PID_correction_sys_uncer[3][i] = 0.0;
        }

        PID_correction_uncer[0][i] = std::sqrt(PID_correction_stat_uncer[0][i] * PID_correction_stat_uncer[0][i] + PID_correction_sys_uncer[0][i] * PID_correction_sys_uncer[0][i]);
        PID_correction_uncer[1][i] = std::sqrt(PID_correction_stat_uncer[1][i] * PID_correction_stat_uncer[1][i] + PID_correction_sys_uncer[1][i] * PID_correction_sys_uncer[1][i]);
        PID_correction_uncer[2][i] = std::sqrt(PID_correction_stat_uncer[2][i] * PID_correction_stat_uncer[2][i] + PID_correction_sys_uncer[2][i] * PID_correction_sys_uncer[2][i]);
        PID_correction_uncer[3][i] = std::sqrt(PID_correction_stat_uncer[3][i] * PID_correction_stat_uncer[3][i] + PID_correction_sys_uncer[3][i] * PID_correction_sys_uncer[3][i]);
    }

    PID_correction[0][N_PID_syst - 1] = 1.0;
    PID_correction_stat_uncer[0][N_PID_syst - 1] = 0.0;
    PID_correction_sys_uncer[0][N_PID_syst - 1] = 0.0;
    PID_correction_uncer[0][N_PID_syst - 1] = 0.0;

    PID_correction[1][N_PID_syst - 1] = 1.0;
    PID_correction_stat_uncer[1][N_PID_syst - 1] = 0.0;
    PID_correction_sys_uncer[1][N_PID_syst - 1] = 0.0;
    PID_correction_uncer[1][N_PID_syst - 1] = 0.0;

    PID_correction[2][N_PID_syst - 1] = 1.0;
    PID_correction_stat_uncer[2][N_PID_syst - 1] = 0.0;
    PID_correction_sys_uncer[2][N_PID_syst - 1] = 0.0;
    PID_correction_uncer[2][N_PID_syst - 1] = 0.0;

    PID_correction[3][N_PID_syst - 1] = 1.0;
    PID_correction_stat_uncer[3][N_PID_syst - 1] = 0.0;
    PID_correction_sys_uncer[3][N_PID_syst - 1] = 0.0;
    PID_correction_uncer[3][N_PID_syst - 1] = 0.0;

    fclose(fp_KID_true);
    fclose(fp_KID_mis);
    fclose(fp_PID_true);
    fclose(fp_PID_mis);
}

void FluctuatePIDCorrection(bool IsItKID) {

    for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {

        std::lognormal_distribution<double> KID_true_distribution(0.0, PID_correction_uncer[0][i_PID] / PID_correction[0][i_PID]);
        std::lognormal_distribution<double> KID_mis_distribution(0.0, PID_correction_uncer[1][i_PID] / PID_correction[1][i_PID]);
        std::lognormal_distribution<double> PID_true_distribution(0.0, PID_correction_uncer[2][i_PID] / PID_correction[2][i_PID]);
        std::lognormal_distribution<double> PID_mis_distribution(0.0, PID_correction_uncer[3][i_PID] / PID_correction[3][i_PID]);

        if (IsItKID) {
            if (std::abs(PID_correction_uncer[0][i_PID] / PID_correction[0][i_PID]) < MyEPSILON) PID_correction_fluctuated[0][i_PID] = 1.0; // true KID
            else PID_correction_fluctuated[0][i_PID] = KID_true_distribution(generator);
            if (std::abs(PID_correction_uncer[1][i_PID] / PID_correction[1][i_PID]) < MyEPSILON) PID_correction_fluctuated[1][i_PID] = 1.0; // mis KID
            else PID_correction_fluctuated[1][i_PID] = KID_mis_distribution(generator);
            PID_correction_fluctuated[2][i_PID] = 1.0;
            PID_correction_fluctuated[3][i_PID] = 1.0;
        }
        else {
            PID_correction_fluctuated[0][i_PID] = 1.0;
            PID_correction_fluctuated[1][i_PID] = 1.0;
            if (std::abs(PID_correction_uncer[2][i_PID] / PID_correction[2][i_PID]) < MyEPSILON) PID_correction_fluctuated[2][i_PID] = 1.0; // true PID
            else PID_correction_fluctuated[2][i_PID] = PID_true_distribution(generator);
            if (std::abs(PID_correction_uncer[3][i_PID] / PID_correction[3][i_PID]) < MyEPSILON) PID_correction_fluctuated[3][i_PID] = 1.0; // mis PID
            else PID_correction_fluctuated[3][i_PID] = PID_mis_distribution(generator);
        }

    }
}

void PID_calculator()
{
    ReadPIDFile();

    RooRandom::randomGenerator()->SetSeed(rd());

    double Nevt_nominal[RarityBins * 7] = { 0.0 }; // CHG MIX UUBAR DDBAR SSBAR CHARM SIGNAL
    double Nevt_fluc_KID[NToys][RarityBins * 7] = { 0.0 }; // CHG MIX UUBAR DDBAR SSBAR CHARM SIGNAL
    double Nevt_fluc_PID[NToys][RarityBins * 7] = { 0.0 }; // CHG MIX UUBAR DDBAR SSBAR CHARM SIGNAL
    double Relative_Uncertainty_KID[NToys][RarityBins * 7] = { 0.0 };
    double Relative_Uncertainty_PID[NToys][RarityBins * 7] = { 0.0 };
    double Covariance_KID[RarityBins * 7][RarityBins * 7] = { 0.0 };
    double Covariance_PID[RarityBins * 7][RarityBins * 7] = { 0.0 };

    /* ====================================== */
    // Seting CDF module
    cdf.initbypath("/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Knunu", Scale_Kplus_test);
    cdf.initbypath("/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Kstarnunu", Scale_Kplusstar_test);
    cdf.initbypath("/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Xsnunu", Scale_Xsu_nonresonant_test);
    cdf.initbypath("/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0nunu", Scale_K0_test);
    cdf.initbypath("/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0starnunu", Scale_K0star_test);
    cdf.initbypath("/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02Xsnunu", Scale_Xsd_nonresonant_test);
    /* ====================================== */



    /* ====================================== */
    // define TH1D for temporary usage
    TH1D* temp_hist = new TH1D("temp_hist", "temp_hist", RarityBins, 0.7, 1.0);
    temp_hist->Reset();
    /* ====================================== */



    /* ====================================== */
    // define path for Ntuple
    const char* MC_dirname_Knunu = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Knunu";
    const char* MC_dirname_Kstarnunu = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Kstarnunu";
    const char* MC_dirname_Xsununu = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Xsnunu";
    const char* MC_dirname_K0nunu = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0nunu";
    const char* MC_dirname_K0starnunu = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0starnunu";
    const char* MC_dirname_Xsdnunu = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/SIGNAL_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02Xsnunu";

    const char* MC_dirname_CHG = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/CHG_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_MIX = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/MIX_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_UUBAR = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/UUBAR_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_DDBAR = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/DDBAR_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_SSBAR = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/SSBAR_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_CHARM = "/home/jwpark/storage/BKG_gbasf2/Kasen_again/CHARM_analysis/test_v005/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    /* ====================================== */



    /* ====================================== */
    // get nominal Nevt
    GetNominalNevt(MC_dirname_Knunu, temp_hist, "Bplus", "SIGNAL", Nevt_nominal, Scale_Kplus_test);
    temp_hist->Reset();
    GetNominalNevt(MC_dirname_Kstarnunu, temp_hist, "Bplus", "SIGNAL", Nevt_nominal, Scale_Kplusstar_test);
    temp_hist->Reset();
    GetNominalNevt(MC_dirname_Xsununu, temp_hist, "Bplus", "SIGNAL", Nevt_nominal, Scale_Xsu_nonresonant_test);
    temp_hist->Reset();
    GetNominalNevt(MC_dirname_K0nunu, temp_hist, "Bzero", "SIGNAL", Nevt_nominal, Scale_K0_test);
    temp_hist->Reset();
    GetNominalNevt(MC_dirname_K0starnunu, temp_hist, "Bzero", "SIGNAL", Nevt_nominal, Scale_K0star_test);
    temp_hist->Reset();
    GetNominalNevt(MC_dirname_Xsdnunu, temp_hist, "Bzero", "SIGNAL", Nevt_nominal, Scale_Xsd_nonresonant_test);
    temp_hist->Reset();

    GetNominalNevt(MC_dirname_CHG, temp_hist, "Bplus", "CHG", Nevt_nominal, Scale_CHG_test);
    temp_hist->Reset();
    GetNominalNevt(MC_dirname_MIX, temp_hist, "Bzero", "MIX", Nevt_nominal, Scale_MIX_test);
    temp_hist->Reset();
    GetNominalNevt(MC_dirname_UUBAR, temp_hist, "Continuum", "UUBAR", Nevt_nominal, Scale_UUBAR_test);
    temp_hist->Reset();
    GetNominalNevt(MC_dirname_DDBAR, temp_hist, "Continuum", "DDBAR", Nevt_nominal, Scale_DDBAR_test);
    temp_hist->Reset();
    GetNominalNevt(MC_dirname_SSBAR, temp_hist, "Continuum", "SSBAR", Nevt_nominal, Scale_SSBAR_test);
    temp_hist->Reset();
    GetNominalNevt(MC_dirname_CHARM, temp_hist, "Continuum", "CHARM", Nevt_nominal, Scale_CHARM_test);
    temp_hist->Reset();
    /* ====================================== */



    /* ====================================== */
    // get fluctuated Nevt for KID
    for (int i = 0; i < NToys; i++) {
        FluctuatePIDCorrection(true);

        GetFlucNevt(MC_dirname_Knunu, temp_hist, "Bplus", "SIGNAL", Nevt_fluc_KID, i, true, Scale_Kplus_test);
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_Kstarnunu, temp_hist, "Bplus", "SIGNAL", Nevt_fluc_KID, i, true, Scale_Kplusstar_test);
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_Xsununu, temp_hist, "Bplus", "SIGNAL", Nevt_fluc_KID, i, true, Scale_Xsu_nonresonant_test);
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_K0nunu, temp_hist, "Bzero", "SIGNAL", Nevt_fluc_KID, i, true, Scale_K0_test);
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_K0starnunu, temp_hist, "Bzero", "SIGNAL", Nevt_fluc_KID, i, true, Scale_K0star_test);
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_Xsdnunu, temp_hist, "Bzero", "SIGNAL", Nevt_fluc_KID, i, true, Scale_Xsd_nonresonant_test);
        temp_hist->Reset();

        GetFlucNevt(MC_dirname_CHG, temp_hist, "Bplus", "CHG", Nevt_fluc_KID, i, true, Scale_CHG_test);
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_MIX, temp_hist, "Bzero", "MIX", Nevt_fluc_KID, i, true, Scale_MIX_test);
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_UUBAR, temp_hist, "Continuum", "UUBAR", Nevt_fluc_KID, i, true, Scale_UUBAR_test);
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_DDBAR, temp_hist, "Continuum", "DDBAR", Nevt_fluc_KID, i, true, Scale_DDBAR_test);
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_SSBAR, temp_hist, "Continuum", "SSBAR", Nevt_fluc_KID, i, true, Scale_SSBAR_test);
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_CHARM, temp_hist, "Continuum", "CHARM", Nevt_fluc_KID, i, true, Scale_CHARM_test);
        temp_hist->Reset();
    }
    /* ====================================== */



    /* ====================================== */
    // get fluctuated Nevt for PID
    for (int i = 0; i < NToys; i++) {
        FluctuatePIDCorrection(false);

        GetFlucNevt(MC_dirname_Knunu, temp_hist, "Bplus", "SIGNAL", Nevt_fluc_PID, i, false, Scale_Kplus_test);
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_Kstarnunu, temp_hist, "Bplus", "SIGNAL", Nevt_fluc_PID, i, false, Scale_Kplusstar_test);
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_Xsununu, temp_hist, "Bplus", "SIGNAL", Nevt_fluc_PID, i, false, Scale_Xsu_nonresonant_test);
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_K0nunu, temp_hist, "Bzero", "SIGNAL", Nevt_fluc_PID, i, false, Scale_K0_test);
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_K0starnunu, temp_hist, "Bzero", "SIGNAL", Nevt_fluc_PID, i, false, Scale_K0star_test);
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_Xsdnunu, temp_hist, "Bzero", "SIGNAL", Nevt_fluc_PID, i, false, Scale_Xsd_nonresonant_test);
        temp_hist->Reset();

        GetFlucNevt(MC_dirname_CHG, temp_hist, "Bplus", "CHG", Nevt_fluc_PID, i, false, Scale_CHG_test);
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_MIX, temp_hist, "Bzero", "MIX", Nevt_fluc_PID, i, false, Scale_MIX_test);
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_UUBAR, temp_hist, "Continuum", "UUBAR", Nevt_fluc_PID, i, false, Scale_UUBAR_test);
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_DDBAR, temp_hist, "Continuum", "DDBAR", Nevt_fluc_PID, i, false, Scale_DDBAR_test);
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_SSBAR, temp_hist, "Continuum", "SSBAR", Nevt_fluc_PID, i, false, Scale_SSBAR_test);
        temp_hist->Reset();
        GetFlucNevt(MC_dirname_CHARM, temp_hist, "Continuum", "CHARM", Nevt_fluc_PID, i, false, Scale_CHARM_test);
        temp_hist->Reset();
    }
    /* ====================================== */


    /* ====================================== */
    // get relative uncertainty
    for (int i = 0; i < NToys; i++) {
        for (int j = 0; j < RarityBins * 7; j++) {
            if (std::abs(Nevt_nominal[j]) < MyEPSILON) Relative_Uncertainty_KID[i][j] = 1.0;
            else Relative_Uncertainty_KID[i][j] = Nevt_fluc_KID[i][j] / Nevt_nominal[j];
        }
    }
    for (int i = 0; i < NToys; i++) {
        for (int j = 0; j < RarityBins * 7; j++) {
            if (std::abs(Nevt_nominal[j]) < MyEPSILON) Relative_Uncertainty_PID[i][j] = 1.0;
            else Relative_Uncertainty_PID[i][j] = Nevt_fluc_PID[i][j] / Nevt_nominal[j];
        }
    }
    /* ====================================== */



    /* ====================================== */
    // file output
    FILE* fp;
    
    fp = fopen("KID_toys.txt","w");
    for (int i = 0; i < NToys; i++) {
        for (int j = 0; j < RarityBins * 7; j++) {
            fprintf(fp, "%lf ", Relative_Uncertainty_KID[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    fp = fopen("PID_toys.txt", "w");
    for (int i = 0; i < NToys; i++) {
        for (int j = 0; j < RarityBins * 7; j++) {
            fprintf(fp, "%lf ", Relative_Uncertainty_PID[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    /* ====================================== */

    delete temp_hist;
}
