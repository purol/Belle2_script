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
using namespace RooFit;
using std::string;
using std::to_string;

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

# define Lumi_validation_MC 0.3 // ab-1

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

# define KS0_rel_uncertainty 0.6 // %/cm
# define track_rel_uncertainty 0.69 // %
# define pi0_correction 0.932
# define pi0_rel_uncertainty ((0.0369 / 0.932) * 100.0) // %
# define Kaon_PID_max_uncertainty 0.1 // not percentage. relative uncertainty
// https://indico.belle2.org/event/6872/contributions/37447/attachments/17127/25504/FEIperformance_B2GM.pdf
//# define FEI_cal_Bc 0.679
//# define FEI_cal_Bc_uncertainty (0.017/FEI_cal_Bc)
//# define FEI_cal_B0 0.713
//# define FEI_cal_B0_uncertainty (0.019/FEI_cal_B0)
//
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

# define FEI_cal_Bc 1.0
# define FEI_cal_Bc_uncertainty 0.0
# define FEI_cal_B0 1.0
# define FEI_cal_B0_uncertainty 0.0

# define Htransition_Xsu_change 0.0
# define Htransition_Xsd_change 0.0
# define Ltransition_Xsu_change 0.0
# define Ltransition_Xsd_change 0.0
# define Hmb_Xsu_change 0.0
# define Hmb_Xsd_change 0.0
# define Lmb_Xsu_change 0.0
# define Lmb_Xsd_change 0.0
# define Hpf_Xsu_change 0.0
# define Hpf_Xsd_change 0.0
# define Lpf_Xsu_change 0.0
# define Lpf_Xsd_change 0.0

# define RarityBins 100

# define nsig_min (-280)
# define nsig_max 200

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
            float MVA_Continuum = 0;

            TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");

            TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
            tree_upsilon->SetBranchAddress("MVA_Continuum", &MVA_Continuum);

            for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
                tree_upsilon->GetEntry(j);

                m_values.push_back(MVA_Continuum);
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
            float MVA_Continuum = 0;

            TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");

            TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
            tree_upsilon->SetBranchAddress("MVA_Continuum", &MVA_Continuum);

            for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
                tree_upsilon->GetEntry(j);

                if (IsSignal) {
                    m_signal_output.push_back(GetCDFvalue(MVA_Continuum));
                    m_signal_weight.push_back(weight);
                }
                else {
                    m_background_output.push_back(GetCDFvalue(MVA_Continuum));
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

double transition_point_uncertainty(double N_K, double N_Kstar, double N_Xsu, double N_K0, double N_K0star, double N_Xsd) {

    // calculate efficiency
    double efficiency[3] = { 0.0 };
    const double total_Xs_Num = Lumi_validation_MC * (2.0 * N_BpBp_1invab * BR_Xs_nunubar + 2.0 * N_B0B0_1invab * BR_Xs_nunubar);
    efficiency[0] = (N_K + N_Kstar + N_Xsu + N_K0 + N_K0star + N_Xsd) / total_Xs_Num;
    efficiency[1] = (N_K + N_Kstar + N_Xsu * (1 + Htransition_Xsu_change) + N_K0 + N_K0star + N_Xsd * (1 + Htransition_Xsd_change)) / total_Xs_Num; // High transition point
    efficiency[2] = (N_K + N_Kstar + N_Xsu * (1 + Ltransition_Xsu_change) + N_K0 + N_K0star + N_Xsd * (1 + Ltransition_Xsd_change)) / total_Xs_Num; // Low transition point

    // get relative uncertainty of efficiency
    double efficiency_relative_uncertainty = std::max(std::abs((efficiency[0] - efficiency[1]) / efficiency[0]), std::abs((efficiency[0] - efficiency[2]) / efficiency[0]));

    // show summary
    printf("relative uncertainty from K*-Xs transition point: +-%lf\n", efficiency_relative_uncertainty);

    return efficiency_relative_uncertainty;
}

double mb_uncertainty(double N_K, double N_Kstar, double N_Xsu, double N_K0, double N_K0star, double N_Xsd) {

    // calculate efficiency
    double efficiency[3] = { 0.0 };
    const double total_Xs_Num = Lumi_validation_MC * (2.0 * N_BpBp_1invab * BR_Xs_nunubar + 2.0 * N_B0B0_1invab * BR_Xs_nunubar);
    efficiency[0] = (N_K + N_Kstar + N_Xsu + N_K0 + N_K0star + N_Xsd) / total_Xs_Num;
    efficiency[1] = (N_K + N_Kstar + N_Xsu * (1 + Hmb_Xsu_change) + N_K0 + N_K0star + N_Xsd * (1 + Hmb_Xsd_change)) / total_Xs_Num; // High transition point
    efficiency[2] = (N_K + N_Kstar + N_Xsu * (1 + Lmb_Xsu_change) + N_K0 + N_K0star + N_Xsd * (1 + Lmb_Xsd_change)) / total_Xs_Num; // Low transition point

    // get relative uncertainty of efficiency
    double efficiency_relative_uncertainty = std::max(std::abs((efficiency[0] - efficiency[1]) / efficiency[0]), std::abs((efficiency[0] - efficiency[2]) / efficiency[0]));

    // show summary
    printf("relative uncertainty from mb: +-%lf\n", efficiency_relative_uncertainty);

    return efficiency_relative_uncertainty;
}

double pf_uncertainty(double N_K, double N_Kstar, double N_Xsu, double N_K0, double N_K0star, double N_Xsd) {

    // calculate efficiency
    double efficiency[3] = { 0.0 };
    const double total_Xs_Num = Lumi_validation_MC * (2.0 * N_BpBp_1invab * BR_Xs_nunubar + 2.0 * N_B0B0_1invab * BR_Xs_nunubar);
    efficiency[0] = (N_K + N_Kstar + N_Xsu + N_K0 + N_K0star + N_Xsd) / total_Xs_Num;
    efficiency[1] = (N_K + N_Kstar + N_Xsu * (1 + Hpf_Xsu_change) + N_K0 + N_K0star + N_Xsd * (1 + Hpf_Xsd_change)) / total_Xs_Num; // High transition point
    efficiency[2] = (N_K + N_Kstar + N_Xsu * (1 + Lpf_Xsu_change) + N_K0 + N_K0star + N_Xsd * (1 + Lpf_Xsd_change)) / total_Xs_Num; // Low transition point

    // get relative uncertainty of efficiency
    double efficiency_relative_uncertainty = std::max(std::abs((efficiency[0] - efficiency[1]) / efficiency[0]), std::abs((efficiency[0] - efficiency[2]) / efficiency[0]));

    // show summary
    printf("relative uncertainty from pf: +-%lf\n", efficiency_relative_uncertainty);

    return efficiency_relative_uncertainty;
}

double fraction_uncertainty(double N_K, double N_Kstar, double N_Xsu, double N_K0, double N_K0star, double N_Xsd) {

    // calculate efficiency
    double efficiency[7] = { 0.0 };
    const double total_Xs_Num = Lumi_validation_MC * (2.0 * N_BpBp_1invab * BR_Xs_nunubar + 2.0 * N_B0B0_1invab * BR_Xs_nunubar);
    efficiency[0] = (N_K + N_Kstar + N_Xsu + N_K0 + N_K0star + N_Xsd) / total_Xs_Num;
    efficiency[1] = (N_K * ((BR_Kplus_nunubar + Sigma_BR_Kplus_nunubar) / BR_Kplus_nunubar) + N_Kstar + N_Xsu * ((BR_Xsu_nonresonant_nunubar - Sigma_BR_Kplus_nunubar) / BR_Xsu_nonresonant_nunubar) + N_K0 * ((BR_K0_nunubar + Sigma_BR_K0_nunubar) / BR_K0_nunubar) + N_K0star + N_Xsd * ((BR_Xsd_nonresonant_nunubar - Sigma_BR_K0_nunubar) / BR_Xsd_nonresonant_nunubar)) / total_Xs_Num;
    efficiency[2] = (N_K * ((BR_Kplus_nunubar - Sigma_BR_Kplus_nunubar) / BR_Kplus_nunubar) + N_Kstar + N_Xsu * ((BR_Xsu_nonresonant_nunubar + Sigma_BR_Kplus_nunubar) / BR_Xsu_nonresonant_nunubar) + N_K0 * ((BR_K0_nunubar - Sigma_BR_K0_nunubar) / BR_K0_nunubar) + N_K0star + N_Xsd * ((BR_Xsd_nonresonant_nunubar + Sigma_BR_K0_nunubar) / BR_Xsd_nonresonant_nunubar)) / total_Xs_Num;
    efficiency[3] = (N_K + N_Kstar * ((BR_Kplusstar_nunubar + Sigma_BR_Kplusstar_nunubar) / BR_Kplusstar_nunubar) + N_Xsu * ((BR_Xsu_nonresonant_nunubar - Sigma_BR_Kplusstar_nunubar) / BR_Xsu_nonresonant_nunubar) + N_K0 + N_K0star * ((BR_K0star_nunubar + Sigma_BR_K0star_nunubar) / BR_K0star_nunubar) + N_Xsd * ((BR_Xsd_nonresonant_nunubar - Sigma_BR_K0star_nunubar) / BR_Xsd_nonresonant_nunubar)) / total_Xs_Num;
    efficiency[4] = (N_K + N_Kstar * ((BR_Kplusstar_nunubar - Sigma_BR_Kplusstar_nunubar) / BR_Kplusstar_nunubar) + N_Xsu * ((BR_Xsu_nonresonant_nunubar + Sigma_BR_Kplusstar_nunubar) / BR_Xsu_nonresonant_nunubar) + N_K0 + N_K0star * ((BR_K0star_nunubar - Sigma_BR_K0star_nunubar) / BR_K0star_nunubar) + N_Xsd * ((BR_Xsd_nonresonant_nunubar + Sigma_BR_K0star_nunubar) / BR_Xsd_nonresonant_nunubar)) / total_Xs_Num;
    const double total_Xs_Num_plus = Lumi_validation_MC * (2.0 * N_BpBp_1invab * (BR_Xs_nunubar + Sigma_BR_Xs_nunubar) + 2.0 * N_B0B0_1invab * (BR_Xs_nunubar + Sigma_BR_Xs_nunubar));
    const double total_Xs_Num_minus = Lumi_validation_MC * (2.0 * N_BpBp_1invab * (BR_Xs_nunubar - Sigma_BR_Xs_nunubar) + 2.0 * N_B0B0_1invab * (BR_Xs_nunubar - Sigma_BR_Xs_nunubar));
    efficiency[5] = (N_K + N_Kstar + N_Xsu * ((BR_Xsu_nonresonant_nunubar + Sigma_BR_Xs_nunubar) / BR_Xsu_nonresonant_nunubar) + N_K0 + N_K0star + N_Xsd * ((BR_Xsd_nonresonant_nunubar + Sigma_BR_Xs_nunubar) / BR_Xsd_nonresonant_nunubar)) / total_Xs_Num_plus;
    efficiency[6] = (N_K + N_Kstar + N_Xsu * ((BR_Xsu_nonresonant_nunubar - Sigma_BR_Xs_nunubar) / BR_Xsu_nonresonant_nunubar) + N_K0 + N_K0star + N_Xsd * ((BR_Xsd_nonresonant_nunubar - Sigma_BR_Xs_nunubar) / BR_Xsd_nonresonant_nunubar)) / total_Xs_Num_minus;

    // get relative uncertainty of efficiency
    double efficiency_correction[3] = { 0.0 };
    efficiency_correction[0] = std::max(std::abs((efficiency[0] - efficiency[1]) / efficiency[0]), std::abs((efficiency[0] - efficiency[2]) / efficiency[0]));
    efficiency_correction[1] = std::max(std::abs((efficiency[0] - efficiency[3]) / efficiency[0]), std::abs((efficiency[0] - efficiency[4]) / efficiency[0]));
    efficiency_correction[2] = std::max(std::abs((efficiency[0] - efficiency[5]) / efficiency[0]), std::abs((efficiency[0] - efficiency[6]) / efficiency[0]));
    const double efficiency_relative_uncertainty = std::sqrt(efficiency_correction[0] * efficiency_correction[0] + efficiency_correction[1] * efficiency_correction[1] + efficiency_correction[2] * efficiency_correction[2]);

    // show summary
    printf("relative uncertainty from fraction: %lf (+-%lf +-%lf +-%lf)\n", efficiency_relative_uncertainty, efficiency_correction[0], efficiency_correction[1], efficiency_correction[2]);

    return efficiency_relative_uncertainty;
}

double K_formfactor_uncertainty(const char* dirname, int charge, double weight, double N_K, double N_Kstar, double N_Xsu, double N_K0, double N_K0star, double N_Xsd) {
    if (charge == 0 || charge == 1 || charge == -1) {}
    else {
        printf("charge should be 0 or +-1 for B->K nu nubar decay\n");
        exit(1);
    }

    double m_b = -1; // B meson mass
    double m_k = -1; // Kaon mass

    const double alpha0 = 0.432; // +-0.011
    const double alpha1 = -0.664; // +-0.096
    const double alpha2 = -1.2; // +-0.69
    const double alpha0_sigma = 0.011;
    const double alpha1_sigma = -0.664;
    const double alpha2_sigma = -1.2;
    /*
    <PCA>
    e0 = -0.00586 alpha0 + 0.03681 alpha1 + 0.99931 alpha2
    e1 = 0.05181 alpha0 + 0.99799 alpha1 - 0.03646 alpha2
    e2 = 0.99864 alpha0 -0.05156 alpha2 + 0.00776 alpha2
    lambda0 = 0.69047
    lambda1 = 0.09274
    lambda2 = 0.01007
    */
    const double LinearCoefficients[3][3] = {
        {-0.005860659921244,   0.051812539908648,   0.998639631385367},
        {0.036811794399807,   0.997991063981209, -0.051562854912436},
        {0.999305030738039, -0.036459524432847,   0.007756198796883}
    };
    const double Lambdas[7][3] = {
        {0.0, 0.0, 0.0},
        {0.690465567619637, 0.0, 0.0},
        {-0.690465567619637, 0.0, 0.0},
        {0.0, 0.092737087984295, 0.0},
        {0.0, -0.092737087984295, 0.0},
        {0.0, 0.0, 0.010065727287981},
        {0.0, 0.0, -0.010065727287981}
    };
    double fluctuations[7][3] = {0.0};
    for (int i = 0; i < 7; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) fluctuations[i][j] = fluctuations[i][j] + Lambdas[i][k] * LinearCoefficients[k][j];
        }
    }
    /*
    fluctuations[0][0] = 0.0; // alpha0 change
    fluctuations[0][1] = 0.0; // alpha1 change
    fluctuations[0][2] = 0.0; // alpha2 change
    fluctuations[1][0] = 0.69047 * (-0.00586); // alpha0 change Lambdas[1][0]*LinearCoefficients[0][0] + Lambdas[1][1]*LinearCoefficients[1][0] + Lambdas[1][2]*LinearCoefficients[2][0]
    fluctuations[1][1] = 0.69047 * 0.03681; // alpha1 change Lambdas[1][0]*LinearCoefficients[0][1] + Lambdas[1][1]*LinearCoefficients[1][1] + Lambdas[1][2]*LinearCoefficients[2][1]
    fluctuations[1][2] = 0.69047 * 0.99931; // alpha2 change
    fluctuations[2][0] = (-0.69047) * (-0.00586); // alpha0 change
    fluctuations[2][1] = (-0.69047) * 0.03681; // alpha1 change
    fluctuations[2][2] = (-0.69047) * 0.99931; // alpha2 change
    fluctuations[3][0] = 0.09274 * 0.05181; // alpha0 change
    fluctuations[3][1] = 0.09274 * 0.99799; // alpha1 change
    fluctuations[3][2] = 0.09274 * (-0.03646); // alpha2 change
    fluctuations[4][0] = (-0.09274) * 0.05181; // alpha0 change
    fluctuations[4][1] = (-0.09274) * 0.99799; // alpha1 change
    fluctuations[4][2] = (-0.09274) * (-0.03646); // alpha2 change
    fluctuations[5][0] = 0.01007 * 0.99864; // alpha0 change
    fluctuations[5][1] = 0.01007 * (-0.05156); // alpha1 change
    fluctuations[5][2] = 0.01007 * 0.00775; // alpha2 change
    fluctuations[6][0] = 0.01007 * 0.99864; // alpha0 change
    fluctuations[6][1] = 0.01007 * (-0.05156); // alpha1 change
    fluctuations[6][2] = 0.01007 * 0.00775; // alpha2 change
    */
    double value[7] = { 0.0 }; // value of lambda^1.5 * fp*fp
    double Nevts[7] = { 0.0 }; // number of events at each fluctuations
    double tot_value[7] = { 0.0 }; // value of integral lambda^1.5 * fp*fp
    double Correction_factor_BR[7] = { 0.0 }; // correction factor of BR

    double q2 = -1;

    // load files
    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        //printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_Xs = (TTree*)input_file->Get("Xs");

        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &q2);
        if (charge == 0) {
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spM__bc", &m_b);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &m_k);
        }
        else {
            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spM__bc", &m_b);
            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &m_k);
        }

        //printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_Xs->GetEntries(); j++) { // Fill
            tree_Xs->GetEntry(j);

            q2 = q2 * q2;

            for (int k = 0; k < 7; k++) {
                const double alpha0_fluc = alpha0 + fluctuations[k][0];
                const double alpha1_fluc = alpha1 + fluctuations[k][1];
                const double alpha2_fluc = alpha2 + fluctuations[k][2];

                double mp = m_b + 0.046;
                double tp = (m_b + m_k) * (m_b + m_k);
                double tm = (m_b - m_k) * (m_b - m_k);
                double t0 = tp * (1 - sqrt(1 - tm / tp));
                double z = (sqrt(tp - q2) - sqrt(tp - t0)) / (sqrt(tp - q2) + sqrt(tp - t0));
                double fp = (1 / (1 - q2 / (mp * mp))) * (alpha0_fluc + alpha1_fluc * z + alpha2_fluc * z * z + (-alpha1_fluc + 2 * alpha2_fluc) * z * z * z / 3);
                double lambda = (m_b * m_b * m_b * m_b) + (m_k * m_k * m_k * m_k) + (q2 * q2) - 2 * (m_b * m_b * m_k * m_k + m_b * m_b * q2 + m_k * m_k * q2);

                value[k] = std::pow(lambda, 1.5) * fp * fp;
                Nevts[k] = Nevts[k] + (value[k] / value[0]) * weight;
            }

        }
        input_file->Close();

    }

    // calculate fluctuation of total BR
    for (int k = 0; k < 7; k++) {
        const int step = 2000;
        double m_b_fix = -1;
        double m_k_fix = -1;
        if (charge == 0) {
            m_b_fix = 5.27965;
            m_k_fix = 0.497611;
        }
        else if (charge == 1 || charge == -1) {
            m_b_fix = 5.27934;
            m_k_fix = 0.493677;
        }
        const double q2_min = 0;
        const double q2_max = (m_b_fix - m_k_fix) * (m_b_fix - m_k_fix);
        const double delq2 = (q2_max - q2_min) / step; // max: (mb-mk)^2

        const double alpha0_fluc = alpha0 + fluctuations[k][0];
        const double alpha1_fluc = alpha1 + fluctuations[k][1];
        const double alpha2_fluc = alpha2 + fluctuations[k][2];

        for (int i = 0; i < step; i++) {
            const double q2_integral = q2_min + delq2 * i;

            double mp = m_b_fix + 0.046;
            double tp = (m_b_fix + m_k_fix) * (m_b_fix + m_k_fix);
            double tm = (m_b_fix - m_k_fix) * (m_b_fix - m_k_fix);
            double t0 = tp * (1 - sqrt(1 - tm / tp));
            double z = (sqrt(tp - q2_integral) - sqrt(tp - t0)) / (sqrt(tp - q2_integral) + sqrt(tp - t0));
            double fp = (1 / (1 - q2_integral / (mp * mp))) * (alpha0_fluc + alpha1_fluc * z + alpha2_fluc * z * z + (-alpha1_fluc + 2 * alpha2_fluc) * z * z * z / 3);
            double lambda = (m_b_fix * m_b_fix * m_b_fix * m_b_fix) + (m_k_fix * m_k_fix * m_k_fix * m_k_fix) + (q2_integral * q2_integral) - 2 * (m_b_fix * m_b_fix * m_k_fix * m_k_fix + m_b_fix * m_b_fix * q2_integral + m_k_fix * m_k_fix * q2_integral);

            tot_value[k] = tot_value[k] + std::pow(lambda, 1.5) * fp * fp * delq2;
        }
        if (charge == 0) {
            Correction_factor_BR[k] = tot_value[k] / tot_value[0];
        }
        else if (charge == 1 || charge == -1) {
            Correction_factor_BR[k] = tot_value[k] / tot_value[0];
        }
    }

    /*
    double Min = DBL_MAX;
    double Max = -1;
    for (int i = 1; i < 7; i++) {
        if (Nevts[i] > Max) Max = Nevts[i];
        if (Nevts[i] < Min) Min = Nevts[i];
    }
    if (charge == 0) {
        printf("B0->K0 nu nubar num evt: %lf + %lf -%lf\n", Nevts[0], Max - Nevts[0], Nevts[0] - Min);
    }
    else if (charge == 1 || charge == -1) {
        printf("B+->K+ nu nubar num evt: %lf + %lf -%lf\n", Nevts[0], Max - Nevts[0], Nevts[0] - Min);
    }
    */

    // calculate efficiency
    double efficiency[7] = { 0.0 };
    for (int i = 0; i < 7; i++) {
        const double total_Xs_Num = Lumi_validation_MC * (2.0 * N_BpBp_1invab * BR_Xs_nunubar + 2.0 * N_B0B0_1invab * BR_Xs_nunubar);
        if (charge == 0) {
            const double non_res_Xsd_correction_factor = (BR_Xsd_nonresonant_nunubar + BR_K0_nunubar - BR_K0_nunubar * Correction_factor_BR[i]) / BR_Xsd_nonresonant_nunubar;
            efficiency[i] = (N_K + N_Kstar + N_Xsu + Nevts[i] + N_K0star + N_Xsd * non_res_Xsd_correction_factor) / total_Xs_Num;
        }
        else if (charge == 1 || charge == -1) {
            const double non_res_Xsu_correction_factor = (BR_Xsu_nonresonant_nunubar + BR_Kplus_nunubar - BR_Kplus_nunubar * Correction_factor_BR[i]) / BR_Xsu_nonresonant_nunubar;
            efficiency[i] = (Nevts[i] + N_Kstar + N_Xsu * non_res_Xsu_correction_factor + N_K0 + N_K0star + N_Xsd) / total_Xs_Num;
        }
    }

    // get relative uncertainty of efficiency
    double efficiency_correction[3] = { 0.0 };
    efficiency_correction[0] = std::max(std::abs((efficiency[0] - efficiency[1]) / efficiency[0]), std::abs((efficiency[0] - efficiency[2]) / efficiency[0]));
    efficiency_correction[1] = std::max(std::abs((efficiency[0] - efficiency[3]) / efficiency[0]), std::abs((efficiency[0] - efficiency[4]) / efficiency[0]));
    efficiency_correction[2] = std::max(std::abs((efficiency[0] - efficiency[5]) / efficiency[0]), std::abs((efficiency[0] - efficiency[6]) / efficiency[0]));
    const double efficiency_relative_uncertainty = std::sqrt(efficiency_correction[0]* efficiency_correction[0]+ efficiency_correction[1] * efficiency_correction[1]+ efficiency_correction[2] * efficiency_correction[2]);

    // show summary
    if (charge == 0) {
        printf("B0->K0 nu nubar relative uncertainty: %lf (+-%lf +-%lf +-%lf)\n", efficiency_relative_uncertainty, efficiency_correction[0], efficiency_correction[1], efficiency_correction[2]);
    }
    else if (charge == 1 || charge == -1) {
        printf("B+->K+ nu nubar relative uncertainty: %lf (+-%lf +-%lf +-%lf)\n", efficiency_relative_uncertainty, efficiency_correction[0], efficiency_correction[1], efficiency_correction[2]);
    }

    return efficiency_relative_uncertainty;
}

double Kstar_formfactor_uncertainty(const char* dirname, int charge, double weight, double N_K, double N_Kstar, double N_Xsu, double N_K0, double N_K0star, double N_Xsd) {
    if (charge == 0 || charge == 1 || charge == -1) {}
    else {
        printf("charge should be 0 or +-1 for B->K nu nubar decay\n");
        exit(1);
    }

    double m_b = -1; // B meson mass
    double m_k = -1; // Kaon star mass
    double costheta = -100; // costheta

    const double alpha0_A1 = 0.3;
    const double alpha1_A1 = 0.39;
    const double alpha2_A1 = 1.19;
    const double alpha0_A12 = 0.27;
    const double alpha1_A12 = 0.53;
    const double alpha2_A12 = 0.48;
    const double alpha0_v0 = 0.38;
    const double alpha1_v0 = -1.17;
    const double alpha2_v0 = 2.42;
    const double mR_A1 = 5.829;
    const double mR_A12 = 5.829;
    const double mR_v0 = 5.415;
    /*
    <PCA>
    e0 = - 0.00477 A1a0 + 0.01377 A1a1 + 0.29659 A1a2 + 0.00314 A12a0 + 0.02574 A12a1 + 0.11862 A12a2 - 0.00662 Va0 + 0.02683 Va1 + 0.94674 Va2
    e1 = 0.00955 A1a0 + 0.14336 A1a1 + 0.83303 A1a2 + 0.00253 A12a0 + 0.04200 A12a1 + 0.40190 A12a2 + 0.01359 Va0 + 0.14322 Va1 - 0.31847 Va2
    e2 = 0.00490 A1a0 - 0.02719 A1a1 - 0.43464 A1a2 - 0.00490 A12a0 + 0.14067 A12a1 + 0.88604 A12a2 + 0.00321 Va0 + 0.07099 Va1 + 0.01977 Va2
    e = 0.A1a0 + A1a1 + A1a2 + A12a0 + A12a1 + A12a2 + Va0 + Va1 + Va2
    e = A1a0 + A1a1 + A1a2 + A12a0 + A12a1 + A12a2 + Va0 + Va1 + Va2
    e = A1a0 + A1a1 + A1a2 + A12a0 + A12a1 + A12a2 + Va0 + Va1 + Va2
    e = A1a0 + A1a1 + A1a2 + A12a0 + A12a1 + A12a2 + Va0 + Va1 + Va2
    e = A1a0 + A1a1 + A1a2 + A12a0 + A12a1 + A12a2 + Va0 + Va1 + Va2
    e = A1a0 + A1a1 + A1a2 + A12a0 + A12a1 + A12a2 + Va0 + Va1 + Va2
    lambda0 = 1.57728
    lambda1 = 1.06016
    lambda2 = 0.52088
    lambda3 = 0.23284
    lambda4 = 0.08996
    lambda5 = 0.05088
    lambda6 = 0.01842
    lambda7 = 0.00740
    lambda8 = 0.00177
    */
    const double LinearCoefficients[9][9] = {
    {-0.004767833929074,   0.013772238042338,   0.296589616507393,   0.003139670856264,   0.025743553369441,   0.118620792973715, -0.006615160985686,   0.026825078669141,   0.946739334135628},
    {0.009545743410425,   0.143364495729883,   0.833028724155009,   0.002528027950777,   0.041995533293814,   0.401904805723855,   0.013592144331223,   0.143215823688685, -0.318473920017488},
    {0.004903939104258, -0.027192276054885, -0.434636497785463, -0.004895872506074,   0.140668610588867,   0.886044604252808,   0.003208697626475,   0.070992606469561,   0.019766944252755},
    {0.081726107544427,   0.418943528790113, - 0.148580744618415,   0.008936060075137,   0.016121321335068, - 0.134221136343898,   0.101230277583689,   0.875232126303976,   0.033121258010682},
    {-0.047393076509321, -0.020252979416420,   0.020465460529877,   0.194484738990744,   0.966871985468078, -0.140890698347323, -0.062525020565583, -0.015984271025904, -0.015622549984120},
    {0.130012725270740,   0.880660653289362, -0.078326566909849, -0.023912474389986,   0.027560190412870,   0.018608679240451,   0.008992181636485, -0.446264603521217,   0.022087208404295},
    {0.512122358575286, -0.109791833791786,   0.012308677320526,   0.445594862213020, -0.020341462354386,   0.009165569689839,   0.721186808802633, -0.079569219353840,   0.005540910930872},
    {-0.519183885001133,   0.030474942310095,   0.001530220137484, -0.500080935622905,   0.116801392077100, -0.015739303946084,   0.680871786694491, -0.044168374846196,   0.002926015416780},
    {-0.664990460261152,   0.118522184960495, -0.020121877756987,   0.716134125194759, -0.166939357181738,   0.027983481968456,   0.042738598277210, -0.002952732689141,   0.000271198372049}
    };
    const double Lambdas[19][9] = {
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {1.577276362747197, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {-1.577276362747197, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 1.060164618950133, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, -1.060164618950133, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.520878893100167, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, -0.520878893100167, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.232840750162789, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, -0.232840750162789, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.089960782159514, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, -0.089960782159514, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.050880752840377, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, -0.050880752840377, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.018416666017616, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.018416666017616, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.007399540696482, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.007399540696482, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001770848832452},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.001770848832452}
    };
    double fluctuations[19][9] = {0.0};
    for (int i = 0; i < 19; i++) {
        for (int j = 0; j < 9; j++) {
            for (int k = 0; k < 9; k++) fluctuations[i][j] = fluctuations[i][j] + Lambdas[i][k] * LinearCoefficients[k][j];
        }
    }

    double value[19] = { 0.0 }; // 
    double Nevts[19] = { 0.0 }; // number of events at each fluctuations
    double tot_value[19] = { 0.0 }; // value of integral dBR/dsdcostheta
    double Correction_factor_BR[19] = { 0.0 }; // correction factor of BR

    double q2 = -1;

    // load files
    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        //printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_Xs = (TTree*)input_file->Get("Xs");

        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &q2);
        if (charge == 0) {
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spM__bc", &m_b);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &m_k);
            tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", &costheta);
        }
        else {
            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spM__bc", &m_b);
            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &m_k);
            tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spextraInfo__bohelicityangle__bc__bc", &costheta);
        }

        //printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_Xs->GetEntries(); j++) { // Fill
            tree_Xs->GetEntry(j);

            q2 = q2 * q2;

            for (int k = 0; k < 19; k++) {
                const double alpha0_A1_fluc = alpha0_A1 + fluctuations[k][0];
                const double alpha1_A1_fluc = alpha1_A1 + fluctuations[k][1];
                const double alpha2_A1_fluc = alpha2_A1 + fluctuations[k][2];
                const double alpha0_A12_fluc = alpha0_A12 + fluctuations[k][3];
                const double alpha1_A12_fluc = alpha1_A12 + fluctuations[k][4];
                const double alpha2_A12_fluc = alpha2_A12 + fluctuations[k][5];
                const double alpha0_v0_fluc = alpha0_v0 + fluctuations[k][6];
                const double alpha1_v0_fluc = alpha1_v0 + fluctuations[k][7];
                const double alpha2_v0_fluc = alpha2_v0 + fluctuations[k][8];

                double tp = (m_b + m_k) * (m_b + m_k);
                double tm = (m_b - m_k) * (m_b - m_k);
                double t0 = tp * (1 - sqrt(1 - tm / tp));
                double z = (sqrt(tp - q2) - sqrt(tp - t0)) / (sqrt(tp - q2) + sqrt(tp - t0));
                double z0 = (sqrt(tp) - sqrt(tp - t0)) / (sqrt(tp) + sqrt(tp - t0));

                double v0 = (1 / (1 - q2 / (mR_v0 * mR_v0))) * (alpha0_v0_fluc + alpha1_v0_fluc * (z - z0) + alpha2_v0_fluc * (z - z0) * (z - z0));
                double A1 = (1 / (1 - q2 / (mR_A1 * mR_A1))) * (alpha0_A1_fluc + alpha1_A1_fluc * (z - z0) + alpha2_A1_fluc * (z - z0) * (z - z0));
                double A12 = (1 / (1 - q2 / (mR_A12 * mR_A12))) * (alpha0_A12_fluc + alpha1_A12_fluc * (z - z0) + alpha2_A12_fluc * (z - z0) * (z - z0));
                double lambda = (tp - q2) * (tm - q2);
                double A2 = ((m_b + m_k) * (m_b + m_k) * (m_b * m_b - m_k * m_k - q2) * A1 - A12 * 16 * m_b * m_k * m_k * (m_b + m_k)) / lambda;

                double sB = q2 / (m_b * m_b);
                double m_k_tilda = m_k / m_b;
                double Lambda = 1 + std::pow(m_k_tilda, 4) + sB * sB - 2 * (m_k_tilda * m_k_tilda + sB + sB * m_k_tilda * m_k_tilda);

                double Amp_parallel = -2 * (std::sqrt(sB)) * (std::pow(Lambda, 1.0 / 4.0)) * (std::sqrt(2)) * (1 + m_k_tilda) * A1;
                double Amp_vertical = 2 * (std::sqrt(sB)) * (std::pow(Lambda, 1.0 / 4.0)) * (std::sqrt(2)) * (std::sqrt(Lambda)) * v0 / (1 + m_k_tilda);
                double Amp_0 = -1 * (std::sqrt(sB)) * (std::pow(Lambda, 1.0 / 4.0)) * (1.0 / m_k_tilda) * (1.0 / std::pow(sB, 0.5)) * ((1 - m_k_tilda * m_k_tilda - sB) * (1 + m_k_tilda) * A1 - Lambda * A2 / (1 + m_k_tilda));

                value[k] = (3.0 / 4.0) * (Amp_vertical * Amp_vertical + Amp_parallel * Amp_parallel) * (1 - costheta * costheta) + (3.0 / 2.0) * Amp_0 * Amp_0 * costheta * costheta;
                Nevts[k] = Nevts[k] + (value[k] / value[0]) * weight;
            }

        }
        input_file->Close();

    }

    // calculate fluctuation of total BR
    for (int k = 0; k < 19; k++) {
        const int step = 100;
        double m_b_fix = -1;
        double m_k_fix = -1;
        if (charge == 0) {
            m_b_fix = 5.27965;
            m_k_fix = 0.89555;
        }
        else if (charge == 1 || charge == -1) {
            m_b_fix = 5.27934;
            m_k_fix = 0.89166;
        }
        const double q2_min = 0;
        const double q2_max = (m_b_fix - m_k_fix) * (m_b_fix - m_k_fix);
        const double delq2 = (q2_max - q2_min) / step; // max: (mb-mk)^2

        const double costheta_min = -1;
        const double costheta_max = 1;
        const double delcostheta = (costheta_max - costheta_min) / step;

        const double alpha0_A1_fluc = alpha0_A1 + fluctuations[k][0];
        const double alpha1_A1_fluc = alpha1_A1 + fluctuations[k][1];
        const double alpha2_A1_fluc = alpha2_A1 + fluctuations[k][2];
        const double alpha0_A12_fluc = alpha0_A12 + fluctuations[k][3];
        const double alpha1_A12_fluc = alpha1_A12 + fluctuations[k][4];
        const double alpha2_A12_fluc = alpha2_A12 + fluctuations[k][5];
        const double alpha0_v0_fluc = alpha0_v0 + fluctuations[k][6];
        const double alpha1_v0_fluc = alpha1_v0 + fluctuations[k][7];
        const double alpha2_v0_fluc = alpha2_v0 + fluctuations[k][8];

        for (int i = 0; i < step; i++) {
            const double q2_integral = q2_min + delq2 * i;

            for (int j = 0; j < step; j++) {
                const double costheta_integral = costheta_min + delcostheta * j;

                double tp = (m_b_fix + m_k_fix) * (m_b_fix + m_k_fix);
                double tm = (m_b_fix - m_k_fix) * (m_b_fix - m_k_fix);
                double t0 = tp * (1 - sqrt(1 - tm / tp));
                double z = (sqrt(tp - q2_integral) - sqrt(tp - t0)) / (sqrt(tp - q2_integral) + sqrt(tp - t0));
                double z0 = (sqrt(tp) - sqrt(tp - t0)) / (sqrt(tp) + sqrt(tp - t0));

                double v0 = (1 / (1 - q2_integral / (mR_v0 * mR_v0))) * (alpha0_v0_fluc + alpha1_v0_fluc * (z - z0) + alpha2_v0_fluc * (z - z0) * (z - z0));
                double A1 = (1 / (1 - q2_integral / (mR_A1 * mR_A1))) * (alpha0_A1_fluc + alpha1_A1_fluc * (z - z0) + alpha2_A1_fluc * (z - z0) * (z - z0));
                double A12 = (1 / (1 - q2_integral / (mR_A12 * mR_A12))) * (alpha0_A12_fluc + alpha1_A12_fluc * (z - z0) + alpha2_A12_fluc * (z - z0) * (z - z0));
                double lambda = (tp - q2_integral) * (tm - q2_integral);
                double A2 = ((m_b_fix + m_k_fix) * (m_b_fix + m_k_fix) * (m_b_fix * m_b_fix - m_k_fix * m_k_fix - q2_integral) * A1 - A12 * 16 * m_b_fix * m_k_fix * m_k_fix * (m_b_fix + m_k_fix)) / lambda;

                double sB = q2_integral / (m_b_fix * m_b_fix);
                double m_k_fix_tilda = m_k_fix / m_b_fix;
                double Lambda = 1 + std::pow(m_k_fix_tilda, 4) + sB * sB - 2 * (m_k_fix_tilda * m_k_fix_tilda + sB + sB * m_k_fix_tilda * m_k_fix_tilda);

                double Amp_parallel = -2 * (std::sqrt(sB)) * (std::pow(Lambda, 1.0 / 4.0)) * (std::sqrt(2)) * (1 + m_k_fix_tilda) * A1;
                double Amp_vertical = 2 * (std::sqrt(sB)) * (std::pow(Lambda, 1.0 / 4.0)) * (std::sqrt(2)) * (std::sqrt(Lambda)) * v0 / (1 + m_k_fix_tilda);
                double Amp_0 = -1 * (std::pow(Lambda, 1.0 / 4.0)) * (1.0 / m_k_fix_tilda) * ((1 - m_k_fix_tilda * m_k_fix_tilda - sB) * (1 + m_k_fix_tilda) * A1 - Lambda * A2 / (1 + m_k_fix_tilda));

                tot_value[k] = tot_value[k] + ((3.0 / 4.0) * (Amp_vertical * Amp_vertical + Amp_parallel * Amp_parallel) * (1 - costheta_integral * costheta_integral) + (3.0 / 2.0) * Amp_0 * Amp_0 * costheta_integral * costheta_integral) * delcostheta * (delq2 / (m_b_fix * m_b_fix));
            }
        }
        if (charge == 0) {
            Correction_factor_BR[k] = tot_value[k] / tot_value[0];
        }
        else if (charge == 1 || charge == -1) {
            Correction_factor_BR[k] = tot_value[k] / tot_value[0];
        }
    }

    /*
    double Min = DBL_MAX;
    double Max = -1;
    for (int i = 1; i < 7; i++) {
        if (Nevts[i] > Max) Max = Nevts[i];
        if (Nevts[i] < Min) Min = Nevts[i];
    }
    if (charge == 0) {
        printf("B*0->K*0 nu nubar num evt: %lf + %lf -%lf\n", Nevts[0], Max - Nevts[0], Nevts[0] - Min);
    }
    else if (charge == 1 || charge == -1) {
        printf("B*+->K*+ nu nubar num evt: %lf + %lf -%lf\n", Nevts[0], Max - Nevts[0], Nevts[0] - Min);
    }
    */

    // calculate efficiency
    double efficiency[19] = { 0.0 };
    for (int i = 0; i < 19; i++) {
        const double total_Xs_Num = Lumi_validation_MC * (2.0 * N_BpBp_1invab * BR_Xs_nunubar + 2.0 * N_B0B0_1invab * BR_Xs_nunubar);
        if (charge == 0) {
            const double non_res_Xsd_correction_factor = (BR_Xsd_nonresonant_nunubar + BR_K0star_nunubar - BR_K0star_nunubar * Correction_factor_BR[i]) / BR_Xsd_nonresonant_nunubar;
            efficiency[i] = (N_K + N_Kstar + N_Xsu + N_K0 + Nevts[i] + N_Xsd * non_res_Xsd_correction_factor) / total_Xs_Num;
        }
        else if (charge == 1 || charge == -1) {
            const double non_res_Xsu_correction_factor = (BR_Xsu_nonresonant_nunubar + BR_Kplusstar_nunubar - BR_Kplusstar_nunubar * Correction_factor_BR[i]) / BR_Xsu_nonresonant_nunubar;
            efficiency[i] = (N_K + Nevts[i] + N_Xsu * non_res_Xsu_correction_factor + N_K0 + N_K0star + N_Xsd) / total_Xs_Num;
        }
    }

    // get relative uncertainty of efficiency
    double efficiency_correction[9] = { 0.0 };
    efficiency_correction[0] = std::max(std::abs((efficiency[0] - efficiency[1]) / efficiency[0]), std::abs((efficiency[0] - efficiency[2]) / efficiency[0]));
    efficiency_correction[1] = std::max(std::abs((efficiency[0] - efficiency[3]) / efficiency[0]), std::abs((efficiency[0] - efficiency[4]) / efficiency[0]));
    efficiency_correction[2] = std::max(std::abs((efficiency[0] - efficiency[5]) / efficiency[0]), std::abs((efficiency[0] - efficiency[6]) / efficiency[0]));
    efficiency_correction[3] = std::max(std::abs((efficiency[0] - efficiency[7]) / efficiency[0]), std::abs((efficiency[0] - efficiency[8]) / efficiency[0]));
    efficiency_correction[4] = std::max(std::abs((efficiency[0] - efficiency[9]) / efficiency[0]), std::abs((efficiency[0] - efficiency[10]) / efficiency[0]));
    efficiency_correction[5] = std::max(std::abs((efficiency[0] - efficiency[11]) / efficiency[0]), std::abs((efficiency[0] - efficiency[12]) / efficiency[0]));
    efficiency_correction[6] = std::max(std::abs((efficiency[0] - efficiency[13]) / efficiency[0]), std::abs((efficiency[0] - efficiency[14]) / efficiency[0]));
    efficiency_correction[7] = std::max(std::abs((efficiency[0] - efficiency[15]) / efficiency[0]), std::abs((efficiency[0] - efficiency[16]) / efficiency[0]));
    efficiency_correction[8] = std::max(std::abs((efficiency[0] - efficiency[17]) / efficiency[0]), std::abs((efficiency[0] - efficiency[18]) / efficiency[0]));
    const double efficiency_relative_uncertainty = std::sqrt(
        efficiency_correction[0] * efficiency_correction[0] +
        efficiency_correction[1] * efficiency_correction[1] +
        efficiency_correction[2] * efficiency_correction[2] +
        efficiency_correction[3] * efficiency_correction[3] +
        efficiency_correction[4] * efficiency_correction[4] +
        efficiency_correction[5] * efficiency_correction[5] +
        efficiency_correction[6] * efficiency_correction[6] +
        efficiency_correction[7] * efficiency_correction[7] +
        efficiency_correction[8] * efficiency_correction[8]
    );

    // show summary
    if (charge == 0) {
        printf("B*0->K*0 nu nubar relative uncertainty: %lf", efficiency_relative_uncertainty);
        printf(" (");
        for (int i = 0; i < 9; i++) printf("+-%lf ", efficiency_correction[i]);
        printf(" )");
        printf("\n");
    }
    else if (charge == 1 || charge == -1) {
        printf("B*+->K*+ nu nubar relative uncertainty: %lf", efficiency_relative_uncertainty);
        printf(" (");
        for (int i = 0; i < 9; i++) printf("+-%lf ", efficiency_correction[i]);
        printf(" )");
        printf("\n");
    }

    return efficiency_relative_uncertainty;
}

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

RooRealVar  Rarity_DATA("Rarity", "Rarity_DATA", 0.0, 1.0);
RooRealVar weight_DATA("weight", "weight_DATA", 0.0, 1.0);
RooDataSet info_DATA("info", "info_DATA", RooArgSet(Rarity_DATA, weight_DATA), WeightVar("weight"));

RooRealVar  Rarity_MC_signal("Rarity", "Rarity_MC_signal", 0.0, 1.0);
RooRealVar weight_MC_signal("weight", "weight_MC_signal", 0.0, 1.0);
RooDataSet info_MC_signal("info", "info_MC_signal", RooArgSet(Rarity_MC_signal, weight_MC_signal), WeightVar("weight"));

RooRealVar  Rarity_MC_background("Rarity", "Rarity_MC_background", 0.0, 1.0);
RooRealVar weight_MC_background("weight", "weight_MC_background", 0.0, 1.0);
RooDataSet info_MC_background("info", "info_MC_background", RooArgSet(Rarity_MC_background, weight_MC_background), WeightVar("weight"));

double LetsAdd_Rarity(const char* dirname, RooRealVar* Rarity_, RooRealVar* weight_, RooDataSet* info_, double weight_var = 1.0) {
    float MVA_var = 0;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;

    std::vector<string> names;
    load_files(dirname, &names);

    double Nevt = 0;
    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("MVA_Continuum", &MVA_var); // MVA

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            //            if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) {} // B2Kc
            //            else continue;
            //            if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5 && Mxs < 1.1) {} // B2KcPi0
            //            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5 && Mxs < 1.1) {} // B2Ks0Pic
            //            else continue;
            Nevt = Nevt + weight_var;

            *Rarity_ = cdf.GetCDFvalue(MVA_var);
            *weight_ = weight_var;
            info_->add(RooArgSet(*Rarity_), weight_->getVal());
        }
        input_file->Close();

        printf("%s has %lf events\n", dirname, Nevt);

    }
    return Nevt;
}

void LetsCalculateUncertainties(const char* dirname, double weight) {

    double Upsilon_ID = -1;
    double Bsig_ID = -1;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        double temp_KS0_3D_distance = -1;
        double temp_KaonID_correction = -1;
        double temp_KaonID_rel_up = -1;
        double temp_KaonID_rel_dn = -1;
        double temp_nKaon_excep = -1;

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID); // charged: 0, mixed: 1
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_3D_distance", &temp_KS0_3D_distance);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_correction", &temp_KaonID_correction);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_rel_uncer_up", &temp_KaonID_rel_up);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_Kaon_PID_rel_uncer_dn", &temp_KaonID_rel_dn);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_nKexcep", &temp_nKaon_excep);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) { // B2Kc
                ntracks.push_back(1);
                npi0s.push_back(0);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) { // B2KcPi0
                ntracks.push_back(1);
                npi0s.push_back(1);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) { // B2Ks0Pic
                ntracks.push_back(1);
                npi0s.push_back(0);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) { // B2KcPicPic
                ntracks.push_back(3);
                npi0s.push_back(0);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) { // B2Ks0PicPi0
                ntracks.push_back(1);
                npi0s.push_back(1);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) { // B2KcPicPicPi0
                ntracks.push_back(3);
                npi0s.push_back(1);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) { // B2Ks0PicPicPic
                ntracks.push_back(3);
                npi0s.push_back(0);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) { // B2KcPicPicPicPic
                ntracks.push_back(5);
                npi0s.push_back(0);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) { // B2Ks0PicPicPicPi0
                ntracks.push_back(3);
                npi0s.push_back(1);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) { // B2KcPi0Pi0
                ntracks.push_back(1);
                npi0s.push_back(2);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) { // B2Ks0PicPi0Pi0
                ntracks.push_back(1);
                npi0s.push_back(2);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) { // B2KcPicPicPi0Pi0
                ntracks.push_back(3);
                npi0s.push_back(2);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) { // B2KcKcKc
                ntracks.push_back(3);
                npi0s.push_back(0);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) { // B2KcKcKs0Pic
                ntracks.push_back(3);
                npi0s.push_back(0);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) { // B2KcKcKcPi0
                ntracks.push_back(3);
                npi0s.push_back(1);
                dmID.push_back(0);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) { // B02Ks0
                ntracks.push_back(0);
                npi0s.push_back(0);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) { // B02KcPic
                ntracks.push_back(2);
                npi0s.push_back(0);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) { // B02Ks0Pi0
                ntracks.push_back(0);
                npi0s.push_back(1);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) { // B02KcPicPi0
                ntracks.push_back(2);
                npi0s.push_back(1);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) { // B02Ks0PicPic
                ntracks.push_back(2);
                npi0s.push_back(0);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) { // B02KcPicPicPic
                ntracks.push_back(4);
                npi0s.push_back(0);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) { // B02Ks0PicPicPi0
                ntracks.push_back(2);
                npi0s.push_back(1);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 6.5 && Bsig_ID < 7.5) { // B02KcPicPicPicPi0
                ntracks.push_back(4);
                npi0s.push_back(1);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 7.5 && Bsig_ID < 8.5) { // B02Ks0PicPicPicPic
                ntracks.push_back(4);
                npi0s.push_back(0);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 8.5 && Bsig_ID < 9.5) { // B02Ks0Pi0Pi0
                ntracks.push_back(0);
                npi0s.push_back(2);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 9.5 && Bsig_ID < 10.5) { // B02KcPicPi0Pi0
                ntracks.push_back(2);
                npi0s.push_back(2);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 10.5 && Bsig_ID < 11.5) { // B02Ks0PicPicPi0Pi0
                ntracks.push_back(2);
                npi0s.push_back(2);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) { // B02KcKcKs0
                ntracks.push_back(2);
                npi0s.push_back(0);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) { // B02KcKcKcPic
                ntracks.push_back(4);
                npi0s.push_back(0);
                dmID.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) { // B02KcKcKs0Pi0
                ntracks.push_back(2);
                npi0s.push_back(1);
                dmID.push_back(1);
            }
            else {
                printf("[ERROR] unexpected decay ID\n");
                exit(1);
            }
            KS0_3D_distance.push_back(temp_KS0_3D_distance);

            // correct low momentum region for KaonPID
            KaonID_correction.push_back(temp_KaonID_correction * std::pow(-1, temp_nKaon_excep));
            KaonID_rel_up.push_back(temp_KaonID_rel_up + std::sqrt(2) * temp_nKaon_excep + temp_nKaon_excep * Kaon_PID_max_uncertainty);
            KaonID_rel_dn.push_back(temp_KaonID_rel_dn + std::sqrt(2) * temp_nKaon_excep + temp_nKaon_excep * Kaon_PID_max_uncertainty);

            Ns.push_back(weight);
        }
        input_file->Close();

    }
}

std::vector<double> PrintUncertainties() {
    // start to calculate correction factor/uncertainties
    std::vector<double> corrected_Ns;
    std::vector<double> KS0_rel_uncertainties;
    std::vector<double> track_rel_uncertainties;
    std::vector<double> pi0_rel_uncertainties;
    std::vector<double> FEI_rel_uncertainties;
    for (unsigned int j = 0; j < Ns.size(); j++) {
        double FEI_rel_uncertainty = -1;
        double FEI_correction_factor = -1;
        if (dmID.at(j) == 0) {
            FEI_rel_uncertainty = FEI_cal_Bc_uncertainty;
            FEI_correction_factor = FEI_cal_Bc;
        }
        else if (dmID.at(j) == 1) {
            FEI_rel_uncertainty = FEI_cal_B0_uncertainty;
            FEI_correction_factor = FEI_cal_B0;
        }
        else {
            printf("[ERROR] unexpected decay ID\n");
            exit(1);
        }

        corrected_Ns.push_back(pow(pi0_correction, npi0s.at(j)) * KaonID_correction.at(j) * FEI_correction_factor * Ns.at(j) );
        track_rel_uncertainties.push_back(track_rel_uncertainty * ntracks.at(j));
        pi0_rel_uncertainties.push_back(pi0_rel_uncertainty * npi0s.at(j));
        KS0_rel_uncertainties.push_back(KS0_rel_uncertainty * KS0_3D_distance.at(j));
        FEI_rel_uncertainties.push_back(FEI_rel_uncertainty);
    }
    double corrected_N = 0;
    double avg_track_rel_uncertainty = 0;
    double avg_pi0_rel_uncertainty = 0;
    double avg_KS0_rel_uncertainty = 0;
    double avg_KaonID_rel_uncertainty_up = 0;
    double avg_KaonID_rel_uncertainty_dn = 0;
    double avg_FEI_rel_uncertainty = 0;

    for (unsigned int j = 0; j < Ns.size(); j++) {
        corrected_N = corrected_N + corrected_Ns.at(j);
        avg_track_rel_uncertainty = avg_track_rel_uncertainty + corrected_Ns.at(j) * track_rel_uncertainties.at(j);
        avg_pi0_rel_uncertainty = avg_pi0_rel_uncertainty + corrected_Ns.at(j) * pi0_rel_uncertainties.at(j);
        avg_KS0_rel_uncertainty = avg_KS0_rel_uncertainty + corrected_Ns.at(j) * KS0_rel_uncertainties.at(j);
        avg_KaonID_rel_uncertainty_up = avg_KaonID_rel_uncertainty_up + corrected_Ns.at(j) * KaonID_rel_up.at(j);
        avg_KaonID_rel_uncertainty_dn = avg_KaonID_rel_uncertainty_dn + corrected_Ns.at(j) * KaonID_rel_dn.at(j);
        avg_FEI_rel_uncertainty = avg_FEI_rel_uncertainty + corrected_Ns.at(j) * FEI_rel_uncertainties.at(j);
    }
    avg_track_rel_uncertainty = avg_track_rel_uncertainty / corrected_N;
    avg_pi0_rel_uncertainty = avg_pi0_rel_uncertainty / corrected_N;
    avg_KS0_rel_uncertainty = avg_KS0_rel_uncertainty / corrected_N;
    avg_KaonID_rel_uncertainty_up = avg_KaonID_rel_uncertainty_up / corrected_N;
    avg_KaonID_rel_uncertainty_dn = avg_KaonID_rel_uncertainty_dn / corrected_N;
    avg_FEI_rel_uncertainty = avg_FEI_rel_uncertainty / corrected_N;

    double total_N = 0;
    for (unsigned int j = 0; j < Ns.size(); j++) total_N = total_N + Ns.at(j);

    printf("total uncorrected signal num: %lf\n", total_N);
    printf("Average correction factor: %lf\n", corrected_N / total_N);
    printf("Average relative uncertainty from track: %lf%%\n", avg_track_rel_uncertainty);
    printf("Average relative uncertainty from pi0: %lf%%\n", avg_pi0_rel_uncertainty);
    printf("Average relative uncertainty from KS0: %lf%%\n", avg_KS0_rel_uncertainty);
    printf("Average relative uncertainty from Kaon ID: %lf\n", std::max(avg_KaonID_rel_uncertainty_up, avg_KaonID_rel_uncertainty_dn));
    printf("Average relative uncertainty from FEI: %lf\n", avg_FEI_rel_uncertainty);

    std::vector<double> outputs;
    outputs.push_back(total_N);
    outputs.push_back(corrected_N / total_N);
    outputs.push_back(avg_track_rel_uncertainty);
    outputs.push_back(avg_pi0_rel_uncertainty);
    outputs.push_back(avg_KS0_rel_uncertainty);
    outputs.push_back(std::max(avg_KaonID_rel_uncertainty_up, avg_KaonID_rel_uncertainty_dn));
    outputs.push_back(avg_FEI_rel_uncertainty);
    return outputs;
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

void LinearityTest(RooFitResult* r, RooRealVar RarityFit, double BKG_num) {
    const int LT_number = 35;
    const int LT_iterate_number = 1000;

    RooArgSet fitargs = r->floatParsFinal();
    TIterator* iter(fitargs.createIterator());

    double N_nbkg = -1;
    double N_nsig = -1;
    double ERR_nbkg = -1;
    double ERR_nsig = -1;

    for (TObject* a = iter->Next(); a != 0; a = iter->Next()) {
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(a);
        std::string name = rrv->GetName();
        double val = rrv->getVal();
        double err = rrv->getError();
        if (name == std::string("nbkg")) {
            N_nbkg = val;
            ERR_nbkg = err;
        }
        else if (name == std::string("nsig")) {
            N_nsig = val;
            ERR_nsig = err;
        }
    }

    std::vector<double> n_sigs[LT_number];
    std::vector<double> n_sigs_err[LT_number];
    for (int i = 0; i < LT_number; i++) {

        for (int j = 0; j < LT_iterate_number; j++) {
            // produce toy MC sample

            RooPolynomial histpdf_Rarity_signal_fix("pol0_fix", "pol0_fix", RarityFit, RooArgList());
            RooDataSet* d1 = histpdf_Rarity_signal_fix.generate(RarityFit, i);

            RooRealVar pa_fix("pa_fix", "pa_fix", -1.883); // -2.286 -5 0
            RooRealVar pb_fix("pb_fix", "pb_fix", 1.09); // 1.33 -5 5
            RooGenericPdf histpdf_Rarity_background_fix("g_fix", "(1-@0)*exp(@0*@1+@0*@0*@0*@2)", RooArgList(RarityFit, pa_fix, pb_fix));
            int N_nbkg_int = round(BKG_num);
            RooDataSet* genData = histpdf_Rarity_background_fix.generate(RarityFit, N_nbkg_int);

            genData->append(*d1);
//            RooDataHist gen_binned_data_Rarity("gen binned Rarity data", "gen binned Rarity data", RarityFit, *genData);

            // construct fitting function
            RooPolynomial histpdf_Rarity_signal_LT("pol0_LT", "pol0_LT", RarityFit, RooArgList());
            RooRealVar nsig("nsig", "n_{sig}", 26, nsig_min, nsig_max);
            RooExtendPdf esig("esignal", "extended signal p.d.f", histpdf_Rarity_signal_LT, nsig);

            RooRealVar pa_LT("pa_LT", "pa_LT", -1.883); // -2.286 -5 0
            RooRealVar pb_LT("pb_LT", "pb_LT", 1.09); // 1.33 -5 5
            RooGenericPdf histpdf_Rarity_background_LT("g_LT", "(1-@0)*exp(@0*@1+@0*@0*@0*@2)", RooArgList(RarityFit, pa_LT, pb_LT));

            RooRealVar nbkg("nbkg", "number of background events", 1900, 1700, 2100);
            RooExtendPdf ebkg("ebkg", "extended background p.d.f", histpdf_Rarity_background_LT, nbkg);

            RooAddPdf  totalpdf("model", "b+n", RooArgList(ebkg, esig));

            // fit
            RooFitResult* result = totalpdf.fitTo(*genData, Verbose(false), PrintLevel(-1),  Save());

            RooArgSet fitargs_LT = result->floatParsFinal();
            TIterator* iter_LT(fitargs_LT.createIterator());

            for (TObject* a_LT = iter_LT->Next(); a_LT != 0; a_LT = iter_LT->Next()) {
                RooRealVar* rrv_LT = dynamic_cast<RooRealVar*>(a_LT);
                std::string name_LT = rrv_LT->GetName();
                double val_LT = rrv_LT->getVal();
                double err_LT = rrv_LT->getError();
                if (name_LT == std::string("nsig")) {
                    n_sigs[i].push_back(val_LT);
                    n_sigs_err[i].push_back(err_LT);
                }
            }

        }

    }

    // create root file
    TFile* temp_file = new TFile("Linearity.root", "recreate");
    temp_file->cd();
    TTree* temp_tree = new TTree("Linearity", "");
    double nsig_fit[LT_number];
    double nsig_fit_err[LT_number];

    for (int i = 0; i < LT_number;i++) temp_tree->Branch(("nsig_"+ to_string(i)).c_str(), &nsig_fit[i]);
    for (int i = 0; i < LT_number; i++) temp_tree->Branch(("nsig_error_" + to_string(i)).c_str(), &nsig_fit_err[i]);

    for (int i = 0; i < LT_iterate_number; i++) {
        for (int j = 0; j < LT_number; j++) {
            nsig_fit[j] = n_sigs[j].at(i);
            nsig_fit_err[j] = n_sigs_err[j].at(i);
        }
        temp_tree->Fill();
    }
    temp_tree->Write();
    temp_file->Close();

    // print png file
    double Inputnsig[LT_number] = { 0 };
    double Inputnsigerror[LT_number] = { 0 };
    double outputnsig[LT_number] = { 0 };
    double outputnsigerror[LT_number] = { 0 };
    for (int i = 0; i < LT_number; i++) {
        RooRealVar  nsig_roorealvar("nsig_roorealvar", "n_{sig}", nsig_min, nsig_max);
        RooDataSet nsig_RooDataSet("nsig_RooDataSet", "nsig_RooDataSet", RooArgSet(nsig_roorealvar));
        for (int j = 0; j < LT_iterate_number; j++) {
            nsig_roorealvar = n_sigs[i].at(j);
            nsig_RooDataSet.add(RooArgSet(nsig_roorealvar));
        }
        RooRealVar gausmean("gausmean","",i, i - ERR_nsig, i + ERR_nsig);
        RooRealVar gauswidth("gauswidth","", ERR_nsig, 0, 2* ERR_nsig);
        RooGaussian gauss("gauss","gauss", nsig_roorealvar,gausmean,gauswidth);
        RooRealVar nentry("nentry", "number of entries", 1000, 900, 1100);
        RooExtendPdf egauss("egauss", "extended gauss", gauss, nentry);
        RooFitResult* result_LT = egauss.fitTo(nsig_RooDataSet, Save());

        RooPlot* nsigframe = nsig_roorealvar.frame(Bins(400), Title(" "));
        nsig_RooDataSet.plotOn(nsigframe);
        egauss.plotOn(nsigframe, LineColor(kBlue));
        egauss.paramOn(nsigframe);
        TCanvas* c = new TCanvas("nsig_gauss_fit", "nsig_gauss_fit", 600, 600);
        gPad->SetLeftMargin(0.15); nsigframe->GetYaxis()->SetTitleOffset(1.4); nsigframe->Draw(); c->SaveAs(("nsig_"+ to_string(i) +"_distribution.png").c_str());

        Inputnsig[i] = i;
        Inputnsigerror[i] = 0;
        outputnsig[i] = gausmean.getValV();
        outputnsigerror[i] = gauswidth.getValV();

        delete c;
    }
    gStyle->SetOptFit(11); //gStyle->SetStatH(0.05);
    TCanvas* c = new TCanvas("Linearity test canvas", "", 800, 800);
    TGraphErrors* gr = new TGraphErrors(LT_number, Inputnsig, outputnsig, Inputnsigerror, outputnsigerror);
    gr->SetMarkerStyle(21); gr->SetTitle(";input n_{sig};output n_{sig}");
    gr->Fit("pol1"); gr->Draw("AP");
    c->SaveAs("Linearity_test.png");

    gr->SetMinimum(-1);
    gr->SetMaximum(12.0);
    gr->Draw("AP");
    c->Modified(); c->Update(); c->SaveAs("Linearity_test_magnify.png");

    delete gr;
    delete c;

    gStyle->SetOptFit(0); //gStyle->SetStatH();

}

void LinearityTestSeparately(RooFitResult* r, RooRealVar RarityFit, int SIGNAL_num, double BKG_num) {
    const int LT_iterate_number = 1000;

    RooArgSet fitargs = r->floatParsFinal();
    TIterator* iter(fitargs.createIterator());

    double N_nbkg = -1;
    double N_nsig = -1;
    double ERR_nbkg = -1;
    double ERR_nsig = -1;

    for (TObject* a = iter->Next(); a != 0; a = iter->Next()) {
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(a);
        std::string name = rrv->GetName();
        double val = rrv->getVal();
        double err = rrv->getError();
        if (name == std::string("nbkg")) {
            N_nbkg = val;
            ERR_nbkg = err;
        }
        else if (name == std::string("nsig")) {
            N_nsig = val;
            ERR_nsig = err;
        }
    }

    std::vector<double> n_sigs;
    std::vector<double> n_sigs_err;

    for (int j = 0; j < LT_iterate_number; j++) {
        RooRealVar pa_fix("pa_fix", "pa_fix", 1.450);
        RooRealVar pb_fix("pb_fix", "pb_fix", 753);
        RooRealVar pc_fix("pc_fix", "pc_fix", 0.339);
        RooRealVar pd_fix("pd_fix", "pd_fix", 10.0);
        RooGenericPdf histpdf_Rarity_signal_fix("g_fix", "((1-exp(@2*(@0-1)))^@1)*((1-@0)^(@3-1))*exp(@0^@4)", RooArgList(RarityFit, pa_fix, pb_fix, pc_fix, pd_fix));
        RooDataSet* d1 = histpdf_Rarity_signal_fix.generate(RarityFit, SIGNAL_num);

        RooRealVar m0b_fix("m0b_fix", "m0b_fix", 0.99637);
        RooRealVar sigmab_fix("sigmab_fix", "sigmab_fix", 0.00143);
        RooRealVar alphab_fix("alphab_fix", "alphab_fix", 0.0087);
        RooRealVar nb_fix("nb_fix", "nb_fix", 0.35);
        RooCBShape histpdf_Rarity_background_fix("genpdfb_fix", "genpdfb_fix", RarityFit, m0b_fix, sigmab_fix, alphab_fix, nb_fix);
        int N_nbkg_int = round(BKG_num);
        RooDataSet* genData = histpdf_Rarity_background_fix.generate(RarityFit, N_nbkg_int);

        genData->append(*d1);

        RooRealVar pa_LT("pa_LT", "pa_LT", 1.450);
        RooRealVar pb_LT("pb_LT", "pb_LT", 753);
        RooRealVar pc_LT("pc_LT", "pc_LT", 0.339);
        RooRealVar pd_LT("pd_LT", "pd_LT", 10.0);
        RooGenericPdf histpdf_Rarity_signal_LT("g_LT", "((1-exp(@2*(@0-1)))^@1)*((1-@0)^(@3-1))*exp(@0^@4)", RooArgList(RarityFit, pa_LT, pb_LT, pc_LT, pd_LT));
        RooRealVar nsig("nsig", "n_{sig}", 26, nsig_min, nsig_max);
        RooExtendPdf esig("esignal", "extended signal p.d.f", histpdf_Rarity_signal_LT, nsig);

        RooRealVar m0b_LT("m0b_LT", "m0b_LT", 0.99637);
        RooRealVar sigmab_LT("sigmab_LT", "sigmab_LT", 0.00143);
        RooRealVar alphab_LT("alphab_LT", "alphab_LT", 0.0087);
        RooRealVar nb_LT("nb_LT", "nb_LT", 0.35);
        RooCBShape histpdf_Rarity_background_LT("genpdfb_LT", "genpdfb_LT", RarityFit, m0b_LT, sigmab_LT, alphab_LT, nb_LT);

        RooRealVar nbkg("nbkg", "number of background events", 1900, 1700, 2100);
        RooExtendPdf ebkg("ebkg", "extended background p.d.f", histpdf_Rarity_background_LT, nbkg);

        RooAddPdf  totalpdf("model", "b+n", RooArgList(ebkg, esig));

        // fit
        RooFitResult* result = totalpdf.fitTo(*genData, Verbose(false), PrintLevel(-1), Save());

        RooArgSet fitargs_LT = result->floatParsFinal();
        TIterator* iter_LT(fitargs_LT.createIterator());

        for (TObject* a_LT = iter_LT->Next(); a_LT != 0; a_LT = iter_LT->Next()) {
            RooRealVar* rrv_LT = dynamic_cast<RooRealVar*>(a_LT);
            std::string name_LT = rrv_LT->GetName();
            double val_LT = rrv_LT->getVal();
            double err_LT = rrv_LT->getError();
            if (name_LT == std::string("nsig")) {
                n_sigs.push_back(val_LT);
                n_sigs_err.push_back(err_LT);
            }
        }

    }

    // create root file
    TFile* temp_file = new TFile(("Linearity_"+ std::to_string(SIGNAL_num) +".root").c_str(), "recreate");
    temp_file->cd();
    TTree* temp_tree = new TTree("Linearity", "");
    double nsig_fit;
    double nsig_fit_err;

    temp_tree->Branch(("nsig_" + to_string(SIGNAL_num)).c_str(), &nsig_fit);
    temp_tree->Branch(("nsig_error_" + to_string(SIGNAL_num)).c_str(), &nsig_fit_err);

    for (int i = 0; i < LT_iterate_number; i++) {
        nsig_fit = n_sigs.at(i);
        nsig_fit_err = n_sigs_err.at(i);
        temp_tree->Fill();
    }
    temp_tree->Write();
    temp_file->Close();

    // print png file
    double Inputnsig =  0;
    double Inputnsigerror = 0;
    double outputnsig = 0;
    double outputnsigerror = 0;

    RooRealVar  nsig_roorealvar("nsig_roorealvar", "n_{sig}", nsig_min, nsig_max);
    RooDataSet nsig_RooDataSet("nsig_RooDataSet", "nsig_RooDataSet", RooArgSet(nsig_roorealvar));
    for (int j = 0; j < LT_iterate_number; j++) {
        nsig_roorealvar = n_sigs.at(j);
        nsig_RooDataSet.add(RooArgSet(nsig_roorealvar));
    }
    RooRealVar gausmean("gausmean", "", SIGNAL_num, SIGNAL_num - ERR_nsig, SIGNAL_num + ERR_nsig);
    RooRealVar gauswidth("gauswidth", "", ERR_nsig, 0, 2 * ERR_nsig);
    RooGaussian gauss("gauss", "gauss", nsig_roorealvar, gausmean, gauswidth);
    RooRealVar nentry("nentry", "number of entries", 1000, 900, 1100);
    RooExtendPdf egauss("egauss", "extended gauss", gauss, nentry);
    RooFitResult* result_LT = egauss.fitTo(nsig_RooDataSet, Save());

    RooPlot* nsigframe = nsig_roorealvar.frame(Bins(400), Title(" "));
    nsig_RooDataSet.plotOn(nsigframe);
    egauss.plotOn(nsigframe, LineColor(kBlue));
    egauss.paramOn(nsigframe);
    TCanvas* c = new TCanvas("nsig_gauss_fit", "nsig_gauss_fit", 600, 600);
    gPad->SetLeftMargin(0.15); nsigframe->GetYaxis()->SetTitleOffset(1.4); nsigframe->Draw(); c->SaveAs(("nsig_" + to_string(SIGNAL_num) + "_distribution.png").c_str());

    Inputnsig = SIGNAL_num;
    Inputnsigerror = 0;
    outputnsig = gausmean.getValV();
    outputnsigerror = gauswidth.getValV();

    printf("input signal: %d\n", SIGNAL_num);
    printf("output signal: %f\n", gausmean.getValV());
    printf("output signal error: %f\n", gauswidth.getValV());

    delete c;

}

std::vector<double> ToyMCstudy(RooExtendPdf ExtendedSIGNALPDF, RooExtendPdf ExtendedBKGPDF, RooRealVar RarityFit, double SIGNAL_num, double BKG_num) {
    const int TOY_iterate_number = 1000;

    std::vector<double> n_sigs;
    std::vector<double> n_sigs_err;

    for (int j = 0; j < TOY_iterate_number; j++) {
        // produce toy MC sample
        RooDataSet* d1 = ExtendedSIGNALPDF.generate(RarityFit, SIGNAL_num, false, true, "", false, true);
        RooDataSet* genData = ExtendedBKGPDF.generate(RarityFit, BKG_num, false, true, "", false, true);

        genData->append(*d1);
        //            RooDataHist gen_binned_data_Rarity("gen binned Rarity data", "gen binned Rarity data", RarityFit, *genData);

                    // construct fitting function
        RooPolynomial histpdf_Rarity_signal_TOY("pol0_TOY", "pol0_TOY", RarityFit, RooArgList());
        RooRealVar nsig("nsig", "n_{sig}", 26.0, nsig_min, nsig_max);
        RooExtendPdf esig("esignal", "extended signal p.d.f", histpdf_Rarity_signal_TOY, nsig);

        RooRealVar pa_TOY("pa_TOY", "pa_TOY", -1.883); // -2.286 -5 0
        RooRealVar pb_TOY("pb_TOY", "pb_TOY", 1.09); // 1.33 -5 5
        RooGenericPdf histpdf_Rarity_background_TOY("g_TOY", "(1-@0)*exp(@0*@1+@0*@0*@0*@2)", RooArgList(RarityFit, pa_TOY, pb_TOY));

        RooRealVar nbkg("nbkg", "number of background events", 1900, 1700, 2100); // 1900 1700 2100
        RooExtendPdf ebkg("ebkg", "extended background p.d.f", histpdf_Rarity_background_TOY, nbkg);

        RooAddPdf  totalpdf("model", "b+n", RooArgList(ebkg, esig));

        // fit
        RooFitResult* result = totalpdf.fitTo(*genData, Verbose(false), PrintLevel(-1), Save());

        RooArgSet fitargs_TOY = result->floatParsFinal();
        TIterator* iter_TOY(fitargs_TOY.createIterator());

        for (TObject* a_TOY = iter_TOY->Next(); a_TOY != 0; a_TOY = iter_TOY->Next()) {
            RooRealVar* rrv_TOY = dynamic_cast<RooRealVar*>(a_TOY);
            std::string name_TOY = rrv_TOY->GetName();
            double val_TOY = rrv_TOY->getVal();
            double err_TOY = rrv_TOY->getError();
            if (name_TOY == std::string("nsig")) {
                n_sigs.push_back(val_TOY);
                n_sigs_err.push_back(err_TOY);
            }
        }

    }

    // print png
    TH1F* ToyMCnsig = new TH1F("ToyMCnsig", ";n_{sig};evt", 40, nsig_min, nsig_max);
    TH1F* ToyMCnsigerror = new TH1F("ToyMCnsigerror", ";error of n_{sig};evt", 50, 10, 60);
    TH1F* ToyMCnsigpull = new TH1F("ToyMCnsigpull", ";pull of n_{sig};evt", 40, -5, 5);

    ToyMCnsig->SetMarkerStyle(kFullCircle);
    ToyMCnsig->SetLineColor(kBlack);
    ToyMCnsig->SetMarkerColor(kBlack);

    ToyMCnsigerror->SetMarkerStyle(kFullCircle);
    ToyMCnsigerror->SetLineColor(kBlack);
    ToyMCnsigerror->SetMarkerColor(kBlack);

    ToyMCnsigpull->SetMarkerStyle(kFullCircle);
    ToyMCnsigpull->SetLineColor(kBlack);
    ToyMCnsigpull->SetMarkerColor(kBlack);

    for (int i = 0; i < TOY_iterate_number; i++) {
        ToyMCnsig->Fill(n_sigs.at(i));
        ToyMCnsigerror->Fill(n_sigs_err.at(i));
        ToyMCnsigpull->Fill((n_sigs.at(i) - SIGNAL_num)/ n_sigs_err.at(i));
    }

//    gStyle->SetOptFit(11); gStyle->SetStatH(0.05);
    TCanvas* c = new TCanvas("canvas_ToyMC_study", "", 800, 800);
    ToyMCnsig->Draw("PE1");
    c->SaveAs("TOYMC_nsig.png");
    delete c;

    c = new TCanvas("canvas_ToyMC_study", "", 800, 800);
    ToyMCnsigerror->Draw("PE1");
    c->SaveAs("TOYMC_nsigerror.png");
    delete c;

    c = new TCanvas("canvas_ToyMC_study", "", 800, 800);
    ToyMCnsigpull->Fit("gaus");
    ToyMCnsigpull->Draw("PE1");
    c->SaveAs("TOYMC_nsigpull.png");
    delete c;

    TF1* fit_pull_gauss = ToyMCnsigpull->GetFunction("gaus");
    const double constant_gauss = fit_pull_gauss->GetParameter(0);
    const double mean_gauss = fit_pull_gauss->GetParameter(1);
    const double sigma_gauss = fit_pull_gauss->GetParameter(2);

    delete ToyMCnsig;
    delete ToyMCnsigerror;
    delete ToyMCnsigpull;

    gStyle->SetOptFit(0); //gStyle->SetStatH();

    std::vector<double> outputs;
    outputs.push_back(constant_gauss);
    outputs.push_back(mean_gauss);
    outputs.push_back(sigma_gauss);

    printf("TOY MC study pull mean: %lf\n", mean_gauss);
    printf("TOY MC study pull sigma: %lf\n", sigma_gauss);
    return outputs;
}

// functions for TF1. it will be used for likelihood
TGraph *TGraph_Likelihood;
TF1* TF1_Likelihood;
TGraph *TGraph_smeared_Likelihood;
TF1* TF1_smeared_Likelihood;
double ConvertToTF1(double *xx, double *){
    return TGraph_Likelihood->Eval(xx[0]);
}

double ConvertToTF1_smeared(double *xx, double *){
    return TGraph_smeared_Likelihood->Eval(xx[0]);
}

void GetLikelihood(RooRealVar* nsig, RooDataSet* d_Rarity, RooAddPdf totalpdf, double uncertainty_additive, double uncertainty_multiplicative){
    RooPlot* nllframe = nsig->frame();
    RooAbsReal* nll = totalpdf.createNLL(*d_Rarity, NumCPU(40), Verbose(false));
    RooAbsReal* pll_nsig = nll->createProfile(*nsig);
    pll_nsig->plotOn(nllframe);

    TF1 *pll_tf;
    pll_tf = pll_nsig->asTF(RooArgList(*nsig));

    const int Nbins = 3000;
    const double xmin = nsig_min;
    const double xmax = nsig_max;
    const double delx = (xmax-xmin)/Nbins;
    double x_pll[Nbins] = {0};
    double y_pll[Nbins] = {0};

    for(int i=0; i<Nbins; i++){
      x_pll[i] = xmin + ((xmax - xmin)*i)/Nbins;
      if(isnan(pll_tf->Eval(x_pll[i])) == 0) { y_pll[i] = 1.0/std::exp(pll_tf->Eval(x_pll[i])); printf("%lf %lf %lf\n",x_pll[i], pll_tf->Eval(x_pll[i]), y_pll[i]); }
      else y_pll[i] = 0.0;
    }
    TGraph_Likelihood = new TGraph(Nbins, x_pll, y_pll);
    TF1_Likelihood = new TF1("TF1_Likelihood", ConvertToTF1, nsig_min, nsig_max);
    //TF1* Gaussian_smear = new TF1("Gauss_smear", GaussianDis, -100, 130);
    //TF1Convolution *smeared_likelihood_conv = new TF1Convolution("TF1_Likelihood","Gauss_smear", -100, 130);
    //TF1* smeared_likelihood = new TF1("smeared_likelihood", *smeared_likelihood_conv, -100, 130, smeared_likelihood_conv->GetNpar());

    // obtain smeared likelihood
    double unc_add = uncertainty_additive;
    double unc_mul = uncertainty_multiplicative;
    double x_pll_smeared[Nbins] = {0};
    double y_pll_smeared[Nbins] = {0};
    for(int i=0; i<Nbins; i++){
        x_pll_smeared[i] = xmin + ((xmax - xmin)*i)/Nbins;
        double integral_value = 0;
        for(int j=0; j<Nbins; j++){ // convolution
            double x_conv = xmin + j * delx;
            double unc_tot = std::sqrt(unc_add*unc_add + unc_mul*x_conv*unc_mul*x_conv);
            double Gaussian_value = std::exp( -0.5*(x_pll_smeared[i]-x_conv)*(x_pll_smeared[i]-x_conv)/( unc_tot*unc_tot ) )/( unc_tot*std::sqrt(2*TMath::Pi()) );
            integral_value = integral_value + TF1_Likelihood->Eval(x_conv) * Gaussian_value * delx;
        }
        y_pll_smeared[i] = integral_value;
    }
    TGraph_smeared_Likelihood = new TGraph(Nbins, x_pll_smeared, y_pll_smeared);
    TF1_smeared_Likelihood = new TF1("TF1_smeared_Likelihood", ConvertToTF1_smeared, nsig_min, nsig_max);;

    // plot -Log profile likelihood
    TCanvas* c_pll = new TCanvas("pll", "pll", 600, 600);
    gPad->SetLeftMargin(0.15); nllframe->GetYaxis()->SetTitleOffset(1.4); nllframe->Draw(); c_pll->SaveAs("PLL.png");
    delete c_pll;

    // calculate CL90 point
    const double total_Area_only_stat = TF1_Likelihood->Integral(0, nsig_max);
    const double total_Area = TF1_smeared_Likelihood->Integral(0, nsig_max);
    double nsig_limit_only_stat = 40;
    double nsig_limit = 20;
    while( total_Area_only_stat*0.9 > TF1_Likelihood->Integral(0, nsig_limit_only_stat) ){
        nsig_limit_only_stat = nsig_limit_only_stat + 0.002;
    }
    nsig_limit = nsig_limit_only_stat;
    while( total_Area*0.9 > TF1_smeared_Likelihood->Integral(0, nsig_limit) ){
        nsig_limit = nsig_limit + 0.002;
    }
    printf("nsig limit only stat: %lf\n", nsig_limit_only_stat);
    printf("nsig limit: %lf\n", nsig_limit);

    // plot profile likelihood
    TCanvas* c_pll_again = new TCanvas("pll_again", "pll_again", 600, 600);
    TF1_Likelihood->SetLineColor(1); TF1_Likelihood->GetHistogram()->GetYaxis()->SetTitle("profile likelihood"); TF1_Likelihood->GetHistogram()->GetXaxis()->SetTitle("n_{sig}");
    TF1_smeared_Likelihood->SetLineColor(2); TF1_smeared_Likelihood->SetLineStyle(7);
    TF1_Likelihood->Draw(""); TF1_smeared_Likelihood->Draw("SAME");
    TLine* line_only_stat = new TLine(nsig_limit_only_stat, 0, nsig_limit_only_stat, 0.6); line_only_stat->SetLineColor(1); line_only_stat->Draw("");
    TLine* line = new TLine(nsig_limit, 0, nsig_limit, 0.6); line->SetLineColor(2); line->SetLineStyle(7); line->Draw("");
    c_pll_again->SaveAs("Likelihood.png");
    delete c_pll_again;
}

void Signal_yield_fit_BDT_Rarity_CLs()
{
    /* ====================================== */
    // Seting CDF module
    cdf.initbypath("./SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Knunu", Scale_Kplus);
    cdf.initbypath("./SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Kstarnunu", Scale_Kplusstar);
    cdf.initbypath("./SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Xsnunu", Scale_Xsu_nonresonant);
    cdf.initbypath("./SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0nunu", Scale_K0);
    cdf.initbypath("./SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0starnunu", Scale_K0star);
    cdf.initbypath("./SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02Xsnunu", Scale_Xsd_nonresonant);
    /* ====================================== */

    // to extract signal yield
    RooRealVar RarityFit("Rarity", "Transformed BDT_{1}", 0.0, 1.0, "");
    RarityFit.setBins(RarityBins);
    RooPlot* Rarityframe = RarityFit.frame(Bins(RarityBins), Title(" "));
//    RooPlot* Rarityframe = RarityFit.frame(Title(" "));

    // get data from root files
    double BKG_total_Num = 0;
    double SIGNAL_total_Num = 0;

    const char* MC_dirname_Knunu = "./SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Knunu";
    const char* MC_dirname_Kstarnunu = "./SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Kstarnunu";
    const char* MC_dirname_Xsununu = "./SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Xsnunu";
    const char* MC_dirname_K0nunu = "./SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0nunu";
    const char* MC_dirname_K0starnunu = "./SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0starnunu";
    const char* MC_dirname_Xsdnunu = "./SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02Xsnunu";
    LetsAdd_Rarity(MC_dirname_Knunu, &Rarity_MC_signal, &weight_MC_signal, &info_MC_signal, Scale_Kplus);
    LetsAdd_Rarity(MC_dirname_Kstarnunu, &Rarity_MC_signal, &weight_MC_signal, &info_MC_signal, Scale_Kplusstar);
    LetsAdd_Rarity(MC_dirname_Xsununu, &Rarity_MC_signal, &weight_MC_signal, &info_MC_signal, Scale_Xsu_nonresonant);
    LetsAdd_Rarity(MC_dirname_K0nunu, &Rarity_MC_signal, &weight_MC_signal, &info_MC_signal, Scale_K0);
    LetsAdd_Rarity(MC_dirname_K0starnunu, &Rarity_MC_signal, &weight_MC_signal, &info_MC_signal, Scale_K0star);
    LetsAdd_Rarity(MC_dirname_Xsdnunu, &Rarity_MC_signal, &weight_MC_signal, &info_MC_signal, Scale_Xsd_nonresonant);

    const char* MC_dirname_CHG = "./CHG_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_MIX = "./MIX_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_UUBAR = "./UUBAR_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_DDBAR = "./DDBAR_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_SSBAR = "./SSBAR_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_CHARM = "./CHARM_analysis/test_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    LetsAdd_Rarity(MC_dirname_CHG, &Rarity_MC_background, &weight_MC_background, &info_MC_background);
    LetsAdd_Rarity(MC_dirname_MIX, &Rarity_MC_background, &weight_MC_background, &info_MC_background);
    LetsAdd_Rarity(MC_dirname_UUBAR, &Rarity_MC_background, &weight_MC_background, &info_MC_background);
    LetsAdd_Rarity(MC_dirname_DDBAR, &Rarity_MC_background, &weight_MC_background, &info_MC_background);
    LetsAdd_Rarity(MC_dirname_SSBAR, &Rarity_MC_background, &weight_MC_background, &info_MC_background);
    LetsAdd_Rarity(MC_dirname_CHARM, &Rarity_MC_background, &weight_MC_background, &info_MC_background);

    const char* DATA_dirname_Knunu = "./SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Knunu";
    const char* DATA_dirname_Kstarnunu = "./SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Kstarnunu";
    const char* DATA_dirname_Xsununu = "./SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Xsnunu";
    const char* DATA_dirname_K0nunu = "./SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0nunu";
    const char* DATA_dirname_K0starnunu = "./SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0starnunu";
    const char* DATA_dirname_Xsdnunu = "./SIGNAL_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02Xsnunu";
    const double Knunu_total_num = LetsAdd_Rarity(DATA_dirname_Knunu, &Rarity_DATA, &weight_DATA, &info_DATA, Scale_Kplus);
    const double Kstarnunu_total_num = LetsAdd_Rarity(DATA_dirname_Kstarnunu, &Rarity_DATA, &weight_DATA, &info_DATA, Scale_Kplusstar);
    const double Xsununu_total_num = LetsAdd_Rarity(DATA_dirname_Xsununu, &Rarity_DATA, &weight_DATA, &info_DATA, Scale_Xsu_nonresonant);
    const double K0nunu_total_num = LetsAdd_Rarity(DATA_dirname_K0nunu, &Rarity_DATA, &weight_DATA, &info_DATA, Scale_K0);
    const double K0starnunu_total_num = LetsAdd_Rarity(DATA_dirname_K0starnunu, &Rarity_DATA, &weight_DATA, &info_DATA, Scale_K0star);
    const double Xsdnunu_total_num = LetsAdd_Rarity(DATA_dirname_Xsdnunu, &Rarity_DATA, &weight_DATA, &info_DATA, Scale_Xsd_nonresonant);
    SIGNAL_total_Num = Knunu_total_num + Kstarnunu_total_num + Xsununu_total_num + K0nunu_total_num + K0starnunu_total_num + Xsdnunu_total_num;

    const char* DATA_dirname_CHG = "./CHG_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_MIX = "./MIX_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_UUBAR = "./UUBAR_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_DDBAR = "./DDBAR_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_SSBAR = "./SSBAR_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_CHARM = "./CHARM_analysis/validation_v000/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    BKG_total_Num += LetsAdd_Rarity(DATA_dirname_CHG, &Rarity_DATA, &weight_DATA, &info_DATA);
    BKG_total_Num += LetsAdd_Rarity(DATA_dirname_MIX, &Rarity_DATA, &weight_DATA, &info_DATA);
    BKG_total_Num += LetsAdd_Rarity(DATA_dirname_UUBAR, &Rarity_DATA, &weight_DATA, &info_DATA);
    BKG_total_Num += LetsAdd_Rarity(DATA_dirname_DDBAR, &Rarity_DATA, &weight_DATA, &info_DATA);
    BKG_total_Num += LetsAdd_Rarity(DATA_dirname_SSBAR, &Rarity_DATA, &weight_DATA, &info_DATA);
    BKG_total_Num += LetsAdd_Rarity(DATA_dirname_CHARM, &Rarity_DATA, &weight_DATA, &info_DATA);


    // define frame and get ready to make pdfs
    Rarity_DATA.setBins(RarityBins);
    RooDataSet* d_Rarity = (RooDataSet*)info_DATA.reduce(RooArgSet(Rarity_DATA));
    // RooDataHist binned_data_Rarity("binned Rarity data", "binned Rarity data", RarityFit, *d_Rarity);

    Rarity_MC_signal.setBins(RarityBins);
    Rarity_MC_background.setBins(RarityBins);

    RooPolynomial histpdf_Rarity_signal("pol0", "pol0", RarityFit, RooArgList());

    RooRealVar nsig("nsig", "n_{sig}", 26.0, nsig_min, nsig_max);
    RooExtendPdf esig("esignal", "extended signal p.d.f", histpdf_Rarity_signal, nsig);

//    RooRealVar pa("pa", "pa", -5.56);
//    RooRealVar pb("pb", "pb", -3.079);
//    RooGenericPdf histpdf_Rarity_background("g", "exp(@0*@1)+exp(@0*@0*@2)", RooArgList(RarityFit, pa, pb));

    RooRealVar pa("pa", "pa", -1.883); // -2.286 -5 0
    RooRealVar pb("pb", "pb", 1.09); // 1.33 -5 5
    RooGenericPdf histpdf_Rarity_background("g", "(1-@0)*exp(@0*@1+@0*@0*@0*@2)", RooArgList(RarityFit, pa, pb));

    RooRealVar nbkg("nbkg", "number of background events", 1900, 1700, 2100); // 1900 1700 2100
    RooExtendPdf ebkg("ebkg", "extended background p.d.f", histpdf_Rarity_background, nbkg);

    RooAddPdf  totalpdf("model", "b+n", RooArgList(ebkg, esig));

    // fit
    RooFitResult* r = totalpdf.fitTo(*d_Rarity, Save(), SumW2Error(true));
    RooArgSet fitargs = r->floatParsFinal();
    TIterator* iter(fitargs.createIterator());

    double nsig_val = -1;
    double nsig_err = -1;
    for (TObject* a = iter->Next(); a != 0; a = iter->Next()) {
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(a);
        std::string name = rrv->GetName();
        double val = rrv->getVal();
        double err = rrv->getError();
        if (name == std::string("nsig")) {
            printf("nsig: %lf +- %lf\n",val, err);
            nsig_val = val;
            nsig_err = err;
        }
    }

    // Draw result
//    d_Rarity->plotOn(Rarityframe, DataError(RooAbsData::Poisson));
    d_Rarity->plotOn(Rarityframe, DataError(RooAbsData::SumW2));
    totalpdf.plotOn(Rarityframe, LineColor(kRed));
    totalpdf.plotOn(Rarityframe, Components(esig), LineColor(kBlue), LineStyle(kDashed));
    totalpdf.plotOn(Rarityframe, Components(ebkg), LineColor(kViolet), LineStyle(kDashed));
    totalpdf.paramOn(Rarityframe);
    
    TCanvas* c = new TCanvas("Rarity", "Rarity", 900, 600);
    gPad->SetLeftMargin(0.15); Rarityframe->GetYaxis()->SetTitleOffset(1.4); Rarityframe->Draw(); c->SaveAs("Rarity_distribution.png");
    delete c;

    
    /* ============== Linearity test ============== */
    //LinearityTest(r, RarityFit, BKG_total_Num);

    /* ============== toy MC study ============== */
    RooPolynomial histpdf_Rarity_signal_TOY_basis("pol0_TOY_basis", "pol0_TOY_basis", RarityFit, RooArgList());
    RooRealVar nsig_TOY_basis("nsig_TOY_basis", "n_{sig}_TOY_basis", SIGNAL_total_Num);
    RooExtendPdf esig_TOY_basis("esignal_TOY_basis", "extended signal p.d.f_TOY_basis", histpdf_Rarity_signal_TOY_basis, nsig_TOY_basis);

    RooRealVar pa_TOY_basis("pa_TOY_basis", "pa_TOY_basis", -1.883); // -2.29 -5 0
    RooRealVar pb_TOY_basis("pb_TOY_basis", "pb_TOY_basis", 1.09); // 1.33 -5 5
    RooGenericPdf histpdf_Rarity_background_TOY_basis("g_TOY_basis", "(1-@0)*exp(@0*@1+@0*@0*@0*@2)", RooArgList(RarityFit, pa_TOY_basis, pb_TOY_basis));
    RooRealVar nbkg_TOY_basis("nbkg_TOY_basis", "number of background events_TOY_basis", BKG_total_Num); // 1900 1700 2100
    RooExtendPdf ebkg_TOY_basis("ebkg_TOY_basis", "extended background p.d.f_TOY_basis", histpdf_Rarity_background_TOY_basis, nbkg_TOY_basis);
    std::vector<double> pull_result = ToyMCstudy(esig_TOY_basis, ebkg_TOY_basis, RarityFit, SIGNAL_total_Num, BKG_total_Num); // const, mean, sigma

    /* ============== Uncertainty Calculation ============== */
    printf("=================== start to calculate uncertainties ===================\n");
    LetsCalculateUncertainties(MC_dirname_Knunu, Scale_Kplus);
    LetsCalculateUncertainties(MC_dirname_Kstarnunu, Scale_Kplusstar);
    LetsCalculateUncertainties(MC_dirname_Xsununu, Scale_Xsu_nonresonant);
    LetsCalculateUncertainties(MC_dirname_K0nunu, Scale_K0);
    LetsCalculateUncertainties(MC_dirname_K0starnunu, Scale_K0star);
    LetsCalculateUncertainties(MC_dirname_Xsdnunu, Scale_Xsd_nonresonant);
    std::vector<double> uncertainties_basic = PrintUncertainties();
    const double K_ff_uncertainty = K_formfactor_uncertainty(MC_dirname_Knunu, 1, Scale_Kplus, Knunu_total_num, Kstarnunu_total_num, Xsununu_total_num, K0nunu_total_num, K0starnunu_total_num, Xsdnunu_total_num);
    const double K0_ff_uncertainty = K_formfactor_uncertainty(MC_dirname_K0nunu, 0, Scale_K0, Knunu_total_num, Kstarnunu_total_num, Xsununu_total_num, K0nunu_total_num, K0starnunu_total_num, Xsdnunu_total_num);
    const double Kstar_ff_uncertainty = Kstar_formfactor_uncertainty(MC_dirname_Kstarnunu, 1, Scale_Kplusstar, Knunu_total_num, Kstarnunu_total_num, Xsununu_total_num, K0nunu_total_num, K0starnunu_total_num, Xsdnunu_total_num);
    const double K0star_ff_uncertainty = Kstar_formfactor_uncertainty(MC_dirname_K0starnunu, 0, Scale_K0star, Knunu_total_num, Kstarnunu_total_num, Xsununu_total_num, K0nunu_total_num, K0starnunu_total_num, Xsdnunu_total_num);
    const double Fraction_uncertainty = fraction_uncertainty(Knunu_total_num, Kstarnunu_total_num, Xsununu_total_num, K0nunu_total_num, K0starnunu_total_num, Xsdnunu_total_num);

    const double transition_point_uncertainty_value = transition_point_uncertainty(Knunu_total_num, Kstarnunu_total_num, Xsununu_total_num, K0nunu_total_num, K0starnunu_total_num, Xsdnunu_total_num);
    const double mb_uncertainty_value = mb_uncertainty(Knunu_total_num, Kstarnunu_total_num, Xsununu_total_num, K0nunu_total_num, K0starnunu_total_num, Xsdnunu_total_num);
    const double pf_uncertainty_value = pf_uncertainty(Knunu_total_num, Kstarnunu_total_num, Xsununu_total_num, K0nunu_total_num, K0starnunu_total_num, Xsdnunu_total_num);

    const double total_uncertainty_mul = std::sqrt(
        (uncertainties_basic.at(2) * 0.01) * (uncertainties_basic.at(2) * 0.01) +
        (uncertainties_basic.at(3) * 0.01) * (uncertainties_basic.at(3) * 0.01) +
        (uncertainties_basic.at(4) * 0.01) * (uncertainties_basic.at(4) * 0.01) +
        uncertainties_basic.at(5) * uncertainties_basic.at(5) +
        uncertainties_basic.at(6) * uncertainties_basic.at(6) +
        (K_ff_uncertainty + K0_ff_uncertainty) * (K_ff_uncertainty + K0_ff_uncertainty) +
        (Kstar_ff_uncertainty + K0star_ff_uncertainty) * (Kstar_ff_uncertainty + K0star_ff_uncertainty) +
        Fraction_uncertainty * Fraction_uncertainty + 
        transition_point_uncertainty_value * transition_point_uncertainty_value +
        mb_uncertainty_value * mb_uncertainty_value +
        pf_uncertainty_value * pf_uncertainty_value
    );
    const double total_uncertainty_add = std::sqrt(
        (nsig_err * pull_result.at(1)) * (nsig_err * pull_result.at(1)) +
        (nsig_err * (pull_result.at(2) - 1)) * (nsig_err * (pull_result.at(2) - 1))
    );

    /* ============== Print profile likelihood  ============== */
    GetLikelihood(&nsig, d_Rarity, totalpdf, total_uncertainty_add, total_uncertainty_mul);

    /* ============== End!  ============== */
    printf("nsig: %lf +- %lf\n",nsig_val, nsig_err);
    printf("total additive %lf\n", total_uncertainty_add);
    printf("total multiplicative factor %lf\n", total_uncertainty_mul);

    /* ============== Work space  ============== */
    RooWorkspace* w = new RooWorkspace("w", "workspace");

    // obs
    w->factory("Rarity[0,1]");

    // pdfs
    w->factory("EXPR::bkg('(1-Rarity)*exp(Rarity*a+Rarity*Rarity*b)',Rarity,a[-1.883],b[1.09])");
    w->factory("Polynomial::sig(Rarity)");

    // systematics
    w->factory("Gaussian::add_sub(0,alpha[-10.0,10.0],1.39)");
    w->factory("Gaussian::mul_sub(0,beta[-0.5,0.5],0.0555)");

    w->factory("sig_eff[0.0014]");
    w->factory("NBB[330000000]");
    w->factory("sigBR[0.00008, 0.0, 0.0005]");
    w->factory("prod::sig_yield(2,NBB,sig_eff,sigBR)");
    w->factory("expr::sig_yield_with_unc('sig_yield*(1+beta)+alpha', sig_yield, alpha, beta)");
    w->factory("bkg_yield[1900, 1700, 2100]");
    w->factory("SUM::model(sig_yield*sig, bkg_yield*bkg)");
    w->factory("PROD::model_with_unc(model,add_sub,mul_sub)");

    RooStats::ModelConfig mc("ModelConfig", w);
    mc.SetPdf(*w->pdf("model_with_unc"));
    mc.SetParametersOfInterest(*w->var("sigBR"));
    mc.SetNuisanceParameters(RooArgSet(*w->var("alpha"), *w->var("beta"), *w->var("bkg_yield")));

    w->import(*d_Rarity, RooFit::Rename("observed_data"));
    w->import(mc);
    w->Print();
    w->writeToFile("CLs_test.root");
}
