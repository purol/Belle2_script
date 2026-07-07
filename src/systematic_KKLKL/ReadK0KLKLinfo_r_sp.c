// last update: 2021-10-13
// for Belle2 data

#include <algorithm>
#include <complex>
#include <string.h>
#include <float.h>
#include <cmath>
#include <string>

#include "TCanvas.h"
#include "RooHist.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TFile.h"
#include "TTree.h"

# define mB0 5.27965 // GeV
# define mKp 0.493677 // GeV
# define mKL0 0.497611 // GeV
# define mKS0 0.497611 // GeV
# define mK0 0.497611 // GeV

# define smax_min 0.0
# define smax_max 23.2
# define smin_min 0.0
# define smin_max 23.2
# define NBin 200
# define NBin_cal 10000

# define MyEPSILON 0.000001

# define spin_f980 0
# define m0_f980 0.965
# define gpi_f980 0.165
# define gK_f980 0.695
# define c_f980 1.0
# define phi_f980 0.09

# define spin_f1710 0
# define m0_f1710 1.724
# define Gamma0_f1710 0.137
double q0_f1710 = 0;
double pstar0_f1710 = 0;
# define c_f1710 1.0
# define phi_f1710 1.11

# define spin_f2010 2
# define m0_f2010 2.011
# define Gamma0_f2010 0.202
double q0_f2010 = 0;
double pstar0_f2010 = 0;
# define c_f2010 1.0
# define phi_f2010 2.5

# define spin_chic0 0
# define m0_chic0 3.41475
# define Gamma0_chic0 0.0102
double q0_chic0 = 0;
double pstar0_chic0 = 0;
# define c_chic0 1.0
# define phi_chic0 0.63

# define alpha_NR (-0.14)
# define c_NR 1.0
# define phi_NR 0.0

# define mK 0.495644 // avg between Kp and K0

# define r 1.5
# define rprime 1.5

# define mpi0 0.1349768
# define mpic 0.13957039

using namespace std::complex_literals;

typedef struct info {
    double smax;
    double smin;
    double smed;
} Information;

// global variable

std::vector<Information> Infos;

TH2D* N_evt = new TH2D("N_evt", ";s_{max} [(GeV/c^{2})^{2}];s_{min} [(GeV/c^{2})^{2}];arbitrary unit", NBin, smax_min, smax_max, NBin, smin_min, smin_max);

TH2D* Prob_PHSP = new TH2D("Prob_PHSP", ";s_{max} [(GeV/c^{2})^{2}];s_{min} [(GeV/c^{2})^{2}];arbitrary unit", NBin, smax_min, smax_max, NBin, smin_min, smin_max);

TH2D* Prob = new TH2D("Prob", ";s_{max} [(GeV/c^{2})^{2}];s_{min} [(GeV/c^{2})^{2}];arbitrary unit", NBin, smax_min, smax_max, NBin, smin_min, smin_max);

TH1D* Prob_1D_sqrt_smin = new TH1D("Prob_1D_sqrt_smin", ";#sqrt{s_{min}};", 20, 1.0, 3.3);
TH1D* Prob_1D_sqrt_smax = new TH1D("Prob_1D_sqrt_smax", ";#sqrt{s_{max}};", 20, 3.1, 4.8);

void load_files(const char *dirname, std::vector<std::string>* names){
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

double Calculate_m(double s13, double s23) { // resonance: 1+2
    return std::sqrt(mB0 * mB0 + mKS0 * mKS0 + mKL0 * mKL0 + mKL0 * mKL0 - s13 - s23);
}

double Calculate_q_mag(double s13, double s23) { // resonance: 1+2
    double m = Calculate_m(s13, s23);
    return std::sqrt(m * m / 4.0 - mKL0 * mKL0);
}

double Calculate_pstar_mag(double s13, double s23) { // resonance: 1+2
    double m = Calculate_m(s13, s23);
    return std::sqrt(std::pow(mB0 * mB0 - m * m - mKS0 * mKS0, 2) / 4.0 - m * m * mKS0 * mKS0) / mB0;
}

double Calculate_p_mag(double s13, double s23) { // resonance: 1+2
    double m = Calculate_m(s13, s23);
    double s12 = m * m;
    return std::sqrt(std::pow(mB0 * mB0 - s12 - mKS0 * mKS0, 2) / (4 * s12) - mKS0 * mKS0);
}

double Calculate_q_dot_p_mag(double s13, double s23) { // resonance: 1+2
    return std::abs((s13 - s23) / 4.0);
}

void GetZeros() {

    // get q0 and pstar0 for f1710
    q0_f1710 = std::sqrt((m0_f1710 * m0_f1710) / 4.0 - mKL0 * mKL0);
    pstar0_f1710 = std::sqrt(std::pow(mB0 * mB0 - m0_f1710 * m0_f1710 - mKS0 * mKS0, 2) / 4.0 - m0_f1710 * m0_f1710 * mKS0 * mKS0) / mB0;

    // get q0 and pstar0 for f2010
    q0_f2010 = std::sqrt((m0_f2010 * m0_f2010) / 4.0 - mKL0 * mKL0);
    pstar0_f2010 = std::sqrt(std::pow(mB0 * mB0 - m0_f2010 * m0_f2010 - mKS0 * mKS0, 2) / 4.0 - m0_f2010 * m0_f2010 * mKS0 * mKS0) / mB0;

    // get q0 and pstar0 for chic0
    q0_chic0 = std::sqrt((m0_chic0 * m0_chic0) / 4.0 - mKL0 * mKL0);
    pstar0_chic0 = std::sqrt(std::pow(mB0 * mB0 - m0_chic0 * m0_chic0 - mKS0 * mKS0, 2) / 4.0 - m0_chic0 * m0_chic0 * mKS0 * mKS0) / mB0;

}

double BlattWeisskopf_qr(double s13, double s23, const char* resonance) { // resonance: 1+2
    // X_L(|q|r) when s12 and s23
    int spin = -1;
    double z0 = -1;
    double z = Calculate_q_mag(s13, s23) * r;

    if (strcmp(resonance, "f980") == 0) {
        spin = spin_f980;
        z0 = -1;
    }
    else if (strcmp(resonance, "f1710") == 0) {
        spin = spin_f1710;
        z0 = q0_f1710 * r;
    }
    else if (strcmp(resonance, "f2010") == 0) {
        spin = spin_f2010;
        z0 = q0_f2010 * r;
    }
    else if (strcmp(resonance, "chic0") == 0) {
        spin = spin_chic0;
        z0 = q0_chic0 * r;
    }
    else {
        printf("[BlattWeisskopf_qr] unsupported resonance\n");
        exit(1);
    }

    if (spin == 0) return 1;
    else if (spin == 1) return std::sqrt((1 + z0 * z0) / (1 + z * z));
    else if (spin == 2) return std::sqrt((9 + 3 * z0 * z0 + z0 * z0 * z0 * z0) / (9 + 3 * z * z + z * z * z * z));
    else {
        printf("[BlattWeisskopf_qr] unsupported spin\n");
        exit(1);
    }

}

double BlattWeisskopf_pstarrprime(double s13, double s23, const char* resonance) { // resonance: 1+2
    // X_L(|p*|r') when s12 and s23
    int spin = -1;
    double z0 = -1;
    double z = Calculate_pstar_mag(s13, s23) * rprime;

    if (strcmp(resonance, "f980") == 0) {
        spin = spin_f980;
        z0 = -1;
    }
    else if (strcmp(resonance, "f1710") == 0) {
        spin = spin_f1710;
        z0 = pstar0_f1710 * rprime;
    }
    else if (strcmp(resonance, "f2010") == 0) {
        spin = spin_f2010;
        z0 = pstar0_f2010 * rprime;
    }
    else if (strcmp(resonance, "chic0") == 0) {
        spin = spin_chic0;
        z0 = pstar0_chic0 * rprime;
    }
    else {
        printf("[BlattWeisskopf_pstarrprime] unsupported resonance\n");
        exit(1);
    }

    if (spin == 0) return 1;
    else if (spin == 1) return std::sqrt((1 + z0 * z0) / (1 + z * z));
    else if (spin == 2) return std::sqrt((9 + 3 * z0 * z0 + z0 * z0 * z0 * z0) / (9 + 3 * z * z + z * z * z * z));
    else {
        printf("[BlattWeisskopf_pstarrprime] unsupported spin\n");
        exit(1);
    }

}

double MassDepWidth(double s13, double s23, const char* resonance) { // resonance: 1+2
    // Gamma(m) when s12 and s23
    int spin = -1;
    double m0;
    double q0;
    double Gamma0;

    double q_mag = Calculate_q_mag(s13, s23);
    double m = Calculate_m(s13, s23);

    if (strcmp(resonance, "f1710") == 0) {
        spin = spin_f1710;
        m0 = m0_f1710;
        q0 = q0_f1710;
        Gamma0 = Gamma0_f1710;
    }
    else if (strcmp(resonance, "f2010") == 0) {
        spin = spin_f2010;
        m0 = m0_f2010;
        q0 = q0_f2010;
        Gamma0 = Gamma0_f2010;
    }
    else if (strcmp(resonance, "chic0") == 0) {
        spin = spin_chic0;
        m0 = m0_chic0;
        q0 = q0_chic0;
        Gamma0 = Gamma0_chic0;
    }
    else {
        printf("[MassDepWidth] unsupported resonance\n");
        exit(1);
    }

    return Gamma0 * std::pow(q_mag / q0, 2 * spin + 1) * (m0 / m) * std::pow(BlattWeisskopf_qr(s13, s23, resonance), 2);
}

std::complex<double> RBW(double s13, double s23, const char* resonance) { // resonance: 1+2
    double m0;
    double Gamma0;
    int spin = -1;

    double m = Calculate_m(s13, s23);

    if (strcmp(resonance, "f1710") == 0) {
        spin = spin_f1710;
        m0 = m0_f1710;
        Gamma0 = Gamma0_f1710;
    }
    else if (strcmp(resonance, "f2010") == 0) {
        spin = spin_f2010;
        m0 = m0_f2010;
        Gamma0 = Gamma0_f2010;
    }
    else if (strcmp(resonance, "chic0") == 0) {
        spin = spin_chic0;
        m0 = m0_chic0;
        Gamma0 = Gamma0_chic0;
    }
    else {
        printf("[RBW] unsupported resonance\n");
        exit(1);
    }

    return 1.0 / (m0 * m0 - m * m - 1i * m0 * MassDepWidth(s13, s23, resonance));
}

std::complex<double> Flatte(double s13, double s23, const char* resonance) { // resonance: 1+2
    double m0;
    double gpi;
    double gK;

    if (strcmp(resonance, "f980") == 0) {
        m0 = m0_f980;
        gpi = gpi_f980;
        gK = gK_f980;
    }
    else {
        printf("[Flatte] unsupported resonance\n");
        exit(1);
    }

    double m = Calculate_m(s13, s23);

    double Gammapipi = gpi * ((1.0 / 3.0) * std::sqrt(1 - 4.0 * mpi0 * mpi0 / (m * m)) + (2.0 / 3.0) * std::sqrt(1 - 4.0 * mpic * mpic / (m * m)));
    double GammaKK = gK * (0.5 * std::sqrt(1 - 4.0 * mKp * mKp / (m * m)) + 0.5 * std::sqrt(1 - 4.0 * mK0 * mK0 / (m * m)));

    return 1.0 / ((m0 * m0 - m * m) - 1i * (Gammapipi + GammaKK));
}

double Zemach(double s13, double s23, const char* resonance) { // resonance: 1+2

    int spin;
    double p_dot_q = Calculate_q_dot_p_mag(s13, s23);
    double p_mag = Calculate_p_mag(s13, s23);
    double q_mag = Calculate_q_mag(s13, s23);

    if (strcmp(resonance, "f980") == 0) {
        spin = spin_f980;
    }
    else if (strcmp(resonance, "f1710") == 0) {
        spin = spin_f1710;
    }
    else if (strcmp(resonance, "f2010") == 0) {
        spin = spin_f2010;
    }
    else if (strcmp(resonance, "chic0") == 0) {
        spin = spin_chic0;
    }
    else {
        printf("[Zemach] unsupported resonance\n");
        exit(1);
    }

    if (spin == 0) return 1.0;
    else if (spin == 2) return (8.0 / 3.0) * (3 * p_dot_q * p_dot_q - p_mag * p_mag * q_mag * q_mag);
    else {
        printf("[Zemach] unsupported spin\n");
        exit(1);
    }

}

std::complex<double> DynamicalAmplitude(double s13, double s23, const char* resonance) { // resonance: 1+2

    if (strcmp(resonance, "f980") == 0) {
        return Flatte(s13, s23, resonance) * BlattWeisskopf_pstarrprime(s13, s23, resonance) * BlattWeisskopf_qr(s13, s23, resonance) * Zemach(s13, s23, resonance);
    }
    if (strcmp(resonance, "f1710") == 0) {
        return RBW(s13, s23, resonance) * BlattWeisskopf_pstarrprime(s13, s23, resonance) * BlattWeisskopf_qr(s13, s23, resonance) * Zemach(s13, s23, resonance);
    }
    else if (strcmp(resonance, "f2010") == 0) {
        return RBW(s13, s23, resonance) * BlattWeisskopf_pstarrprime(s13, s23, resonance) * BlattWeisskopf_qr(s13, s23, resonance) * Zemach(s13, s23, resonance);
    }
    else if (strcmp(resonance, "chic0") == 0) {
        return RBW(s13, s23, resonance) * BlattWeisskopf_pstarrprime(s13, s23, resonance) * BlattWeisskopf_qr(s13, s23, resonance) * Zemach(s13, s23, resonance);
    }
    else {
        printf("[Amplitude] unsupported resonance\n");
        exit(1);
    }

}

std::complex<double> Amplitude(double s13, double s23, const char* resonance) { // resonance: 1+2

    /*
     * my FF
     * f980:  30.222805945
     * f1710: 136.555962030
     * f2010: 310762.546712675
     * chic0: 677.282967269
     * NR:    51.262701105
     *
     * paper's FF
     * f980:  0.44
     * f1710: 0.07
     * f2010: 0.09
     * chic0: 0.07
     * NR:    2.16
     */

    double MagicNumber_f980 = std::sqrt(0.44 / 30.222805945);
    double MagicNumber_f1710 = std::sqrt(0.07 / 136.555962030);
    double MagicNumber_f2010 = std::sqrt(0.09 / 310762.546712675);
    double MagicNumber_chic0 = std::sqrt(0.07 / 677.282967269);
    double MagicNumber_NR = std::sqrt(2.16 / 51.262701105);

    double s12 = mB0 * mB0 + mKS0 * mKS0 + mKL0 * mKL0 + mKL0 * mKL0 - s13 - s23;

    if (strcmp(resonance, "f980") == 0) {
        std::complex<double> a;
        a = c_f980 * std::exp(1i * phi_f980);
        return MagicNumber_f980 * a * (DynamicalAmplitude(s13, s23, resonance) + DynamicalAmplitude(s12, s23, resonance) + DynamicalAmplitude(s12, s13, resonance));
    }
    else if (strcmp(resonance, "f1710") == 0) {
        std::complex<double> a;
        a = c_f1710 * std::exp(1i * phi_f1710);
        return MagicNumber_f1710 * a * (DynamicalAmplitude(s13, s23, resonance) + DynamicalAmplitude(s12, s23, resonance) + DynamicalAmplitude(s12, s13, resonance));
    }
    else if (strcmp(resonance, "f2010") == 0) {
        std::complex<double> a;
        a = c_f2010 * std::exp(1i * phi_f2010);
        return MagicNumber_f2010 * a * (DynamicalAmplitude(s13, s23, resonance) + DynamicalAmplitude(s12, s23, resonance) + DynamicalAmplitude(s12, s13, resonance));
    }
    else if (strcmp(resonance, "chic0") == 0) {
        std::complex<double> a;
        a = c_chic0 * std::exp(1i * phi_chic0);
        return MagicNumber_chic0 * a * (DynamicalAmplitude(s13, s23, resonance) + DynamicalAmplitude(s12, s23, resonance) + DynamicalAmplitude(s12, s13, resonance));
    }
    else if (strcmp(resonance, "NR") == 0) {

        std::complex<double> a;
        a = c_NR * std::exp(1i * phi_NR);

        return MagicNumber_NR * a * (std::exp(alpha_NR * s13) + std::exp(alpha_NR * s12) + std::exp(alpha_NR * s23));
    }
    else {
        printf("[Amplitude] unsupported resonance\n");
        exit(1);
    }

}

double Probability(double s13, double s23) { // resonance: 1+2
    std::complex<double> total_amplitude = Amplitude(s13, s23, "f980") + Amplitude(s13, s23, "f1710") + Amplitude(s13, s23, "f2010") + Amplitude(s13, s23, "chic0") + Amplitude(s13, s23, "NR");
    return std::abs(total_amplitude) * std::abs(total_amplitude);
}

void GetProbabilityHist() {

    for (int i = 0; i < NBin_cal; i++) {
        for (int j = 0; j < NBin_cal; j++) {

            double smax = (i + 0.5) * (smax_max - smax_min) / NBin_cal + smax_min;
            double smin = (j + 0.5) * (smin_max - smin_min) / NBin_cal + smin_min;
            double smed = mB0 * mB0 + mKS0 * mKS0 + mKL0 * mKL0 + mKL0 * mKL0 - smax - smin;
            double Probability_value = Probability(smax, smin) * ((smax_max - smax_min) / NBin_cal) * ((smin_max - smin_min) / NBin_cal);

            // https://halldweb.jlab.org/DocDB/0033/003345/002/dalitz.pdf.
            double DalitzValue = smax * smin * smin + (smax * smax - smax * (mB0 * mB0 + mKS0 * mKS0 + mKL0 * mKL0 + mKL0 * mKL0)) * smin + (mKL0 * mKL0 * mB0 * mB0 - mKS0 * mKS0 * mKL0 * mKL0) * (mB0 * mB0 - mKS0 * mKS0 + mKL0 * mKL0 - mKL0 * mKL0);

            if (DalitzValue < 0 && smax > smin && smax > smed && smin < smed) {
                if (isnan(Probability_value)) continue;
                Prob->Fill(smax, smin, Probability_value);
            }

        }
    }

    for (int i = 0; i < NBin_cal; i++) {
        for (int j = 0; j < NBin_cal; j++) {

            double mmax = (i + 0.5) * (std::sqrt(smax_max) - std::sqrt(smax_min)) / NBin_cal + std::sqrt(smax_min);
            double mmin = (j + 0.5) * (std::sqrt(smin_max) - std::sqrt(smin_min)) / NBin_cal + std::sqrt(smin_min);

            double smax = mmax * mmax;
            double smin = mmin * mmin;
            double smed = mB0 * mB0 + mKS0 * mKS0 + mKL0 * mKL0 + mKL0 * mKL0 - smax - smin;
            double Probability_value = Probability(smax, smin) * ((smax_max - smax_min) / NBin_cal) * ((smin_max - smin_min) / NBin_cal);

            // https://halldweb.jlab.org/DocDB/0033/003345/002/dalitz.pdf.
            double DalitzValue = smax * smin * smin + (smax * smax - smax * (mB0 * mB0 + mKS0 * mKS0 + mKL0 * mKL0 + mKL0 * mKL0)) * smin + (mKL0 * mKL0 * mB0 * mB0 - mKS0 * mKS0 * mKL0 * mKL0) * (mB0 * mB0 - mKS0 * mKS0 + mKL0 * mKL0 - mKL0 * mKL0);

            if (DalitzValue < 0 && smax > smin && smax > smed && smin < smed) {
                if (isnan(Probability_value)) continue;
                Prob_1D_sqrt_smin->Fill(mmin, 4 * mmax * mmin * Probability_value); // multiply Jacobian
                Prob_1D_sqrt_smax->Fill(mmax, 4 * mmax * mmin * Probability_value); // multiply Jacobian
            }

        }
    }

    double norm = 1.0 / Prob->Integral();
    printf("Initial area of Prob: %.9lf\nnormalized to one\n", 1.0 / norm);
    Prob->Scale(norm);

    double norm_1D;

    norm_1D = 1.0 / Prob_1D_sqrt_smin->Integral();
    Prob_1D_sqrt_smin->Scale(norm_1D);

    norm_1D = 1.0 / Prob_1D_sqrt_smax->Integral();
    Prob_1D_sqrt_smax->Scale(norm_1D);
}

void Draw1DPlots() {

    TCanvas* c = new TCanvas("c", "", 700, 700);

    Prob_1D_sqrt_smin->SetStats(0);
    Prob_1D_sqrt_smin->Draw("Hist");
    c->SaveAs("1D_sqrt_smin_K0KLKL.png");

    Prob_1D_sqrt_smax->SetStats(0);
    Prob_1D_sqrt_smax->Draw("Hist");
    c->SaveAs("1D_sqrt_smax_K0KLKL.png");

    delete c;
}

void DrawDalitz(TH2D* Dalitz, const char* name = "Dalitz.png") {

    Dalitz->SetStats(0);

    TCanvas* c = new TCanvas("c", "", 700, 700);
    Dalitz->Draw("COLZ");

    c->SaveAs(name);

    delete c;

}

void FillHist() {

    // initialization
    for (int i = 0; i < NBin; i++) {
        for (int j = 0; j < NBin; j++) {
            N_evt->SetBinContent(i + 1, j + 1, 0.0);
            Prob->SetBinContent(i + 1, j + 1, 0.0);
            Prob_PHSP->SetBinContent(i + 1, j + 1, 0.0);
        }
    }

    for (size_t i = 0; i < Infos.size(); i++) {
        double smax_ = Infos.at(i).smax;
        double smin_ = Infos.at(i).smin;
        double smed_ = Infos.at(i).smed;

        int GLobalBin = 0;

        // find Bin index
        GLobalBin = N_evt->FindBin(smax_, smin_);
        N_evt->SetBinContent(GLobalBin, N_evt->GetBinContent(GLobalBin) + 1.0);

        // find Bin index
        GLobalBin = Prob_PHSP->FindBin(smax_, smin_);
        Prob_PHSP->SetBinContent(GLobalBin, Prob_PHSP->GetBinContent(GLobalBin) + 1.0);

    }

    Prob_PHSP->Scale(1.0 / Prob_PHSP->Integral());
}

void FillInfo(const char* dirname) {
    std::vector<std::string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        // Get Data
        TTree* tree_info = (TTree*)input_file->Get("info");

        int __experiment__;
        int __run__;
        unsigned int __event__;
        int __candidate__;
        int __ncandidates__;

        double m = -1;

        double m12;
        double m13;
        double m23;

        // get event_info
        tree_info->SetBranchAddress("__experiment__", &__experiment__);
        tree_info->SetBranchAddress("__run__", &__run__);
        tree_info->SetBranchAddress("__event__", &__event__);
        tree_info->SetBranchAddress("__candidate__", &__candidate__);
        tree_info->SetBranchAddress("__ncandidates__", &__ncandidates__);

        tree_info->SetBranchAddress("daughterInvariantMass__bo0__cm__sp1__bc", &m12);
        tree_info->SetBranchAddress("daughterInvariantMass__bo0__cm__sp2__bc", &m13);
        tree_info->SetBranchAddress("daughterInvariantMass__bo1__cm__sp2__bc", &m23);

        printf("%lld entries...\n", tree_info->GetEntries());
        for (unsigned int j = 0; j < tree_info->GetEntries(); j++) { // Fill

            tree_info->GetEntry(j);

            double s12_ = m12 * m12;
            double s13_ = m13 * m13;
            double s23_ = m23 * m23;

            double smax = std::max(std::max(s12_, s13_), s23_);
            double smin = std::min(std::min(s12_, s13_), s23_);
            double smed = s12_ + s13_ + s23_ - smax - smin;

            // Makeshift to remove photon production decay
            if (std::abs(s12_ - (mB0 * mB0 + mKS0 * mKS0 + mKL0 * mKL0 + mKL0 * mKL0 - s13_ - s23_)) > MyEPSILON) continue;

            Information temp_info = { smax, smin, smed };

            Infos.push_back(temp_info);

        }
        input_file->Close();
    }
}

void GetWeights() {
    FILE* fp;

    fp = fopen("K0KLKL_weight.txt","w");

    fprintf(fp, "smax: %d %lf %lf\n", NBin, smax_min, smax_max);
    fprintf(fp, "smin: %d %lf %lf\n", NBin, smin_min, smin_max);
    for (int i = 0; i < NBin; i++) {
        for (int j = 0; j < NBin; j++) {
            double smax = (i + 0.5) * (smax_max - smax_min) / NBin + smax_min;
            double smin = (j + 0.5) * (smin_max - smin_min) / NBin + smin_min;

            int GLobalBin_PHSP = Prob_PHSP->FindBin(smax, smin);
            double PHSP_val = Prob_PHSP->GetBinContent(GLobalBin_PHSP);

            int GLobalBin_model = Prob->FindBin(smax, smin);
            double model_val = Prob->GetBinContent(GLobalBin_model);

            int GLobalBin_Nevt = N_evt->FindBin(smax, smin);
            double Nevt = N_evt->GetBinContent(GLobalBin_Nevt);

            if (Nevt < MyEPSILON) fprintf(fp, "%lf %lf %lf\n", smax, smin, 0.0);
            else fprintf(fp, "%lf %lf %lf\n", smax, smin, model_val / PHSP_val);
        }
    }

    fclose(fp);
}

int main(){
    /*
    get correction factor for B+ --> K+ KL0 KL0 MC
    MC:   B+ --> K+ KL0 KL0
    */
    
    const char* dirname = "./files";

    FillInfo(dirname);
    FillHist();
    GetZeros();

    DrawDalitz(N_evt, "PHSP_K0KLKL.png");

    GetProbabilityHist();

    DrawDalitz(Prob, "model_K0KLKL.png");

    Draw1DPlots();

    GetWeights();

    return 0;
}
