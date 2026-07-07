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

# define mB 5.27934 // GeV
# define mKp 0.493677 // GeV
# define mKL0 0.497611 // GeV

# define s12_min 0.0
# define s12_max 23.2
# define s13_min 0.0
# define s13_max 23.2
# define s23_min 0.0
# define s23_max 23.2
# define NBin 200
# define NBin_cal 10000

# define MyEPSILON 0.000001

# define spin_f980 0
# define m0_f980 0.965
# define gpi_f980 0.165
# define gKgpi_f980 4.21
# define c_f980 3.35
# define phi_f980 31.0

# define spin_f1500 0
# define m0_f1500 1.505
# define Gamma0_f1500 0.109
double q0_f1500 = 0;
double pstar0_f1500 = 0;
# define c_f1500 0.2
# define phi_f1500 (-83.0)

# define spin_f1525 2
# define m0_f1525 1.525
# define Gamma0_f1525 0.073
double q0_f1525 = 0;
double pstar0_f1525 = 0;
# define c_f1525 0.00179
# define phi_f1525 (-58.0)

# define spin_f1710 0
# define m0_f1710 1.720
# define Gamma0_f1710 0.135
double q0_f1710 = 0;
double pstar0_f1710 = 0;
# define c_f1710 0.24
# define phi_f1710 (-22.0)

# define spin_chic0 0
# define m0_chic0 3.41475
# define Gamma0_chic0 0.0103
double q0_chic0 = 0;
double pstar0_chic0 = 0;
# define c_chic0 0.113
# define phi_chic0 45.0
# define delta_chic0 (-12.0)

# define b_NR (-0.018)
# define aS0_c_NR 1.0
# define aS0_phi_NR 0.0
# define aS1_c_NR 1.0
# define aS1_phi_NR 129.0
# define aS2_c_NR 0.51
# define aS2_phi_NR (-85.0)

# define mK 0.495644 // avg between Kp and K0
# define mpic 0.13957039

# define r 4.0
# define rprime 0.0

using namespace std::complex_literals;

typedef struct info {
    double s12;
    double s13;
    double s23;
    double m;
    double pstar_mag;
    double p_mag;
    double q_mag;
    double p_dot_q;
} Information;

// global variable

std::vector<Information> Infos;

TH2D* N_evt = new TH2D("N_evt", ";s13 [(GeV/c^{2})^{2}];s23 [(GeV/c^{2})^{2}];arbitrary unit", NBin, s13_min, s13_max, NBin, s23_min, s23_max);
TH2D* N_evt_conv = new TH2D("N_evt_conv", ";s12;s13", NBin, s12_min, s12_max, NBin, s13_min, s13_max);

TH2D* Prob_PHSP = new TH2D("Prob_PHSP", ";s13;s23", NBin, s13_min, s13_max, NBin, s23_min, s23_max);
TH2D* Prob_PHSP_conv = new TH2D("Prob_PHSP_conv", ";s12;s13", NBin, s12_min, s12_max, NBin, s13_min, s13_max);

TH2D* Prob = new TH2D("Prob", ";s13 [(GeV/c^{2})^{2}];s23 [(GeV/c^{2})^{2}];arbitrary unit", NBin, s13_min, s13_max, NBin, s23_min, s23_max);
TH2D* Prob_conv = new TH2D("Prob_conv", ";s12 [(GeV/c^{2})^{2}];s13 [(GeV/c^{2})^{2}];arbitrary unit", NBin, s12_min, s12_max, NBin, s13_min, s13_max);

TH1D* Prob_1D_s12 = new TH1D("Prob_1D_s12", ";m_{K_{S}K_{S}};", 32, 4.77 - 0.06 * 63, 4.77 + 0.06);
TH1D* Prob_1D_s13 = new TH1D("Prob_1D_s13", ";m_{K^{+}K_{S},low};", 29, 0.99, 3.74);
TH1D* Prob_1D_s23 = new TH1D("Prob_1D_s23", ";m_{K^{+}K_{S},high};", 31, 1.685, 4.784);

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

double DegreeToRadian(double degree) {
    return (3.141592 * degree) / 180.0;
}

double Calculate_m(double s13, double s23) {
    return std::sqrt(mB * mB + mKp * mKp + mKL0 * mKL0 + mKL0 * mKL0 - s13 - s23);
}

double Calculate_q_mag(double s13, double s23) {
    double m = Calculate_m(s13, s23);
    return std::sqrt(m * m / 4.0 - mKL0 * mKL0);
}

double Calculate_pstar_mag(double s13, double s23) {
    double m = Calculate_m(s13, s23);
    return std::sqrt(std::pow(mB * mB - m * m - mKp * mKp, 2) / 4.0 - m * m * mKp * mKp) / mB;
}

double Calculate_p_mag(double s13, double s23) {
    double m = Calculate_m(s13, s23);
    double s12 = m * m;
    return std::sqrt(std::pow(mB * mB - s12 - mKp * mKp, 2) / (4 * s12) - mKp * mKp);
}

double Calculate_q_dot_p_mag(double s13, double s23) {
    return std::abs((s13 - s23) / 4.0);
}

void GetZeros() {

    // get q0 and pstar0 for f1500
    q0_f1500 = std::sqrt((m0_f1500 * m0_f1500) / 4.0 - mKL0 * mKL0);
    pstar0_f1500 = std::sqrt(std::pow(mB * mB - m0_f1500 * m0_f1500 - mKp * mKp, 2) / 4.0 - m0_f1500 * m0_f1500 * mKp * mKp) / mB;

    // get q0 and pstar0 for f1525
    q0_f1525 = std::sqrt((m0_f1525 * m0_f1525) / 4.0 - mKL0 * mKL0);
    pstar0_f1525 = std::sqrt(std::pow(mB * mB - m0_f1525 * m0_f1525 - mKp * mKp, 2) / 4.0 - m0_f1525 * m0_f1525 * mKp * mKp) / mB;

    // get q0 and pstar0 for f1710
    q0_f1710 = std::sqrt((m0_f1710 * m0_f1710) / 4.0 - mKL0 * mKL0);
    pstar0_f1710 = std::sqrt(std::pow(mB * mB - m0_f1710 * m0_f1710 - mKp * mKp, 2) / 4.0 - m0_f1710 * m0_f1710 * mKp * mKp) / mB;

    // get q0 and pstar0 for chic0
    q0_chic0 = std::sqrt((m0_chic0 * m0_chic0) / 4.0 - mKL0 * mKL0);
    pstar0_chic0 = std::sqrt(std::pow(mB * mB - m0_chic0 * m0_chic0 - mKp * mKp, 2) / 4.0 - m0_chic0 * m0_chic0 * mKp * mKp) / mB;

}

double BlattWeisskopf_qr(double s13, double s23, const char* resonance) {
    // X_L(|q|r) when s12 and s23
    int spin = -1;
    double z0 = -1;
    double z = Calculate_q_mag(s13, s23) * r;

    if (strcmp(resonance, "f980") == 0) {
        spin = spin_f980;
        z0 = -1;
    }
    else if (strcmp(resonance, "f1500") == 0) {
        spin = spin_f1500;
        z0 = q0_f1500 * r;
    }
    else if (strcmp(resonance, "f1525") == 0) {
        spin = spin_f1525;
        z0 = q0_f1525 * r;
    }
    else if (strcmp(resonance, "f1710") == 0) {
        spin = spin_f1710;
        z0 = q0_f1710 * r;
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

double BlattWeisskopf_pstarrprime(double s13, double s23, const char* resonance) {
    // X_L(|p*|r') when s12 and s23
    int spin = -1;
    double z0 = -1;
    double z = Calculate_pstar_mag(s13, s23) * rprime;

    if (strcmp(resonance, "f980") == 0) {
        spin = spin_f980;
        z0 = -1;
    }
    else if (strcmp(resonance, "f1500") == 0) {
        spin = spin_f1500;
        z0 = pstar0_f1500 * rprime;
    }
    else if (strcmp(resonance, "f1525") == 0) {
        spin = spin_f1525;
        z0 = pstar0_f1525 * rprime;
    }
    else if (strcmp(resonance, "f1710") == 0) {
        spin = spin_f1710;
        z0 = pstar0_f1710 * rprime;
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

double MassDepWidth(double s13, double s23, const char* resonance) {
    // Gamma(m) when s12 and s23
    int spin = -1;
    double m0;
    double q0;
    double Gamma0;

    double q_mag = Calculate_q_mag(s13, s23);
    double m = Calculate_m(s13, s23);

    if (strcmp(resonance, "f1500") == 0) {
        spin = spin_f1500;
        m0 = m0_f1500;
        q0 = q0_f1500;
        Gamma0 = Gamma0_f1500;
    }
    else if (strcmp(resonance, "f1525") == 0) {
        spin = spin_f1525;
        m0 = m0_f1525;
        q0 = q0_f1525;
        Gamma0 = Gamma0_f1525;
    }
    else if (strcmp(resonance, "f1710") == 0) {
        spin = spin_f1710;
        m0 = m0_f1710;
        q0 = q0_f1710;
        Gamma0 = Gamma0_f1710;
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

std::complex<double> RBW(double s13, double s23, const char* resonance) {
    double m0;
    double Gamma0;
    int spin = -1;

    double m = Calculate_m(s13, s23);

    if (strcmp(resonance, "f1500") == 0) {
        spin = spin_f1500;
        m0 = m0_f1500;
        Gamma0 = Gamma0_f1500;
    }
    else if (strcmp(resonance, "f1525") == 0) {
        spin = spin_f1525;
        m0 = m0_f1525;
        Gamma0 = Gamma0_f1525;
    }
    else if (strcmp(resonance, "f1710") == 0) {
        spin = spin_f1710;
        m0 = m0_f1710;
        Gamma0 = Gamma0_f1710;
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

std::complex<double> Flatte(double s13, double s23, const char* resonance) {
    double m0;
    double gpi;
    double gK;

    if (strcmp(resonance, "f980") == 0) {
        m0 = m0_f980;
        gpi = gpi_f980;
        gK = gKgpi_f980 * gpi_f980;
    }
    else {
        printf("[Flatte] unsupported resonance\n");
        exit(1);
    }

    double m = Calculate_m(s13, s23);

    double rho_pipi = std::sqrt(1 - 4 * mpic * mpic / (m * m));
    double rho_KK = std::sqrt(1 - 4 * mK * mK / (m * m));

    return 1.0 / ((m0 * m0 - m * m) - 1i * (gpi * rho_pipi + gK * rho_KK));
}

double Zemach(double s13, double s23, const char* resonance) {

    int spin;
    double p_dot_q = Calculate_q_dot_p_mag(s13, s23);
    double p_mag = Calculate_p_mag(s13, s23);
    double q_mag = Calculate_q_mag(s13, s23);

    if (strcmp(resonance, "f980") == 0) {
        spin = spin_f980;
    }
    else if (strcmp(resonance, "f1500") == 0) {
        spin = spin_f1500;
    }
    else if (strcmp(resonance, "f1525") == 0) {
        spin = spin_f1525;
    }
    else if (strcmp(resonance, "f1710") == 0) {
        spin = spin_f1710;
    }
    else if (strcmp(resonance, "chic0") == 0) {
        spin = spin_chic0;
    }
    else {
        printf("[Zemach] unsupported resonance\n");
        exit(1);
    }

    if (spin == 0) return 1.0;
    else if (spin == 1) return 4 * p_dot_q;
    else if (spin == 2) return (16.0 / 3.0) * (3 * p_dot_q * p_dot_q - p_mag * p_mag * q_mag * q_mag);
    else {
        printf("[Zemach] unsupported spin\n");
        exit(1);
    }

}

std::complex<double> DynamicalAmplitude(double s13, double s23, const char* resonance) {

    if (strcmp(resonance, "f980") == 0) {
        return Flatte(s13, s23, resonance) * BlattWeisskopf_pstarrprime(s13, s23, resonance) * BlattWeisskopf_qr(s13, s23, resonance) * Zemach(s13, s23, resonance);
    }
    else if (strcmp(resonance, "f1500") == 0) {
        return RBW(s13, s23, resonance) * BlattWeisskopf_pstarrprime(s13, s23, resonance) * BlattWeisskopf_qr(s13, s23, resonance) * Zemach(s13, s23, resonance);
    }
    else if (strcmp(resonance, "f1525") == 0) {
        return RBW(s13, s23, resonance) * BlattWeisskopf_pstarrprime(s13, s23, resonance) * BlattWeisskopf_qr(s13, s23, resonance) * Zemach(s13, s23, resonance);
    }
    else if (strcmp(resonance, "f1710") == 0) {
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

std::complex<double> Amplitude(double s13, double s23, const char* resonance, bool isobar = false) {

    /*
    * my FF
    * total: 1368.050340092
    * f980: 1871.741841813 | 1.368181993717517
    * f1500: 54.952240702  | 0.0401682884697829
    * f1525: 76.873254637  | 0.0561918318238424
    * f1710: 59.078317296  | 0.0431843153461934
    * chic0: 67.950681349  | 0.0496697229317091
    * NR: 2515.144152802   | 1.838488013995785
    * 
    * paper's FF
    * f980  | 1.386792452830189
    * f1500 | 0.039622641509434
    * f1525 | 0.0575471698113208
    * f1710 | 0.0452830188679245
    * chic0 | 0.05
    * NR    | 1.867924528301887
    */

    if (strcmp(resonance, "f980") == 0) {
        std::complex<double> a;
        if(isobar == false) a = c_f980 * std::exp(1i * DegreeToRadian(phi_f980));
        else a = c_f980 * std::exp(1i * DegreeToRadian(phi_f980));
        return a * (DynamicalAmplitude(s13, s23, resonance) + DynamicalAmplitude(s23, s13, resonance));
    }
    else if (strcmp(resonance, "f1500") == 0) {
        std::complex<double> a;
        if(isobar == false) a = c_f1500 * std::exp(1i * DegreeToRadian(phi_f1500));
        else a = c_f1500 * std::exp(1i * DegreeToRadian(phi_f1500));
        return a * (DynamicalAmplitude(s13, s23, resonance) + DynamicalAmplitude(s23, s13, resonance));
    }
    else if (strcmp(resonance, "f1525") == 0) {
        std::complex<double> a;
        if(isobar == false) a = c_f1525 * std::exp(1i * DegreeToRadian(phi_f1525));
        else a = c_f1525 * std::exp(1i * DegreeToRadian(phi_f1525));
        return a * (DynamicalAmplitude(s13, s23, resonance) + DynamicalAmplitude(s23, s13, resonance));
    }
    else if (strcmp(resonance, "f1710") == 0) {
        std::complex<double> a;
        if(isobar == false) a = c_f1710 * std::exp(1i * DegreeToRadian(phi_f1710));
        else a = c_f1710 * std::exp(1i * DegreeToRadian(phi_f1710));
        return a * (DynamicalAmplitude(s13, s23, resonance) + DynamicalAmplitude(s23, s13, resonance));
    }
    else if (strcmp(resonance, "chic0") == 0) {
        std::complex<double> a;
        if(isobar == false) a = c_chic0 * std::exp(1i * (DegreeToRadian(phi_chic0) + DegreeToRadian(delta_chic0)));
        else a = c_chic0 * std::exp(1i * (DegreeToRadian(phi_chic0) - DegreeToRadian(delta_chic0)));
        return a * (DynamicalAmplitude(s13, s23, resonance) + DynamicalAmplitude(s23, s13, resonance));
    }
    else if (strcmp(resonance, "NR") == 0) {

        double Omega = 0.5 * (mB + (1.0 / 3.0) * (mKp + mKL0 + mKL0));
        double m12 = std::sqrt(mB * mB + mKp * mKp + mKL0 * mKL0 + mKL0 * mKL0 - s13 - s23); // sqrt(s12)
        double x = m12 - Omega;

        std::complex<double> aS0;
        std::complex<double> aS1;
        std::complex<double> aS2;

        if (isobar == false) {
            aS0 = aS0_c_NR * (1 + b_NR) * std::exp(1i * DegreeToRadian(aS0_phi_NR));
            aS1 = aS1_c_NR * (1 + b_NR) * std::exp(1i * DegreeToRadian(aS1_phi_NR));
            aS2 = aS2_c_NR * (1 + b_NR) * std::exp(1i * DegreeToRadian(aS2_phi_NR));
        }
        else {
            aS0 = aS0_c_NR * (1 - b_NR) * std::exp(1i * DegreeToRadian(aS0_phi_NR));
            aS1 = aS1_c_NR * (1 - b_NR) * std::exp(1i * DegreeToRadian(aS1_phi_NR));
            aS2 = aS2_c_NR * (1 - b_NR) * std::exp(1i * DegreeToRadian(aS2_phi_NR));
        }

        return 2.0 * (aS0 + aS1 * x + aS2 * x * x);
    }
    else {
        printf("[Amplitude] unsupported resonance\n");
        exit(1);
    }

}

double Probability(double s13, double s23) {
    std::complex<double> total_amplitude = Amplitude(s13, s23, "f980", false) + Amplitude(s13, s23, "f1500", false) + Amplitude(s13, s23, "f1525", false) + Amplitude(s13, s23, "f1710", false) + Amplitude(s13, s23, "chic0", false) + Amplitude(s13, s23, "NR", false);
    std::complex<double> total_amplitude_isobar = Amplitude(s13, s23, "f980", true) + Amplitude(s13, s23, "f1500", true) + Amplitude(s13, s23, "f1525", true) + Amplitude(s13, s23, "f1710", true) + Amplitude(s13, s23, "chic0", true) + Amplitude(s13, s23, "NR", true);
    return std::abs(total_amplitude) * std::abs(total_amplitude) + std::abs(total_amplitude_isobar) * std::abs(total_amplitude_isobar);
}

void GetProbabilityHist() {

    for (int i = 0; i < NBin_cal; i++) {
        for (int j = 0; j < NBin_cal; j++) {

            double s13 = (i + 0.5) * (s13_max - s13_min) / NBin_cal + s13_min;
            double s23 = (j + 0.5) * (s23_max - s23_min) / NBin_cal + s23_min;
            double s12 = mB * mB + mKp * mKp + mKL0 * mKL0 + mKL0 * mKL0 - s13 - s23;
            double Probability_value = Probability(s13, s23) * ((s13_max - s13_min) / NBin_cal) * ((s23_max - s23_min) / NBin_cal);

            // https://halldweb.jlab.org/DocDB/0033/003345/002/dalitz.pdf.
            double DalitzValue = s23 * s12 * s12 + (s23 * s23 + (mKp * mKp - mKL0 * mKL0) * (mB * mB - mKL0 * mKL0) - s23 * (mB * mB + mKL0 * mKL0 + mKL0 * mKL0 + mKp * mKp)) * s12 + ((mKL0 * mKL0 * mB * mB - mKL0 * mKL0 * mKp * mKp) * (mB * mB - mKL0 * mKL0 + mKL0 * mKL0 - mKp * mKp));

            if (DalitzValue < 0 && s23 >= s13) {
                if (isnan(Probability_value)) continue;
                Prob->Fill(s13, s23, Probability_value);
            }

        }
    }

    for (int i = 0; i < NBin_cal; i++) {
        for (int j = 0; j < NBin_cal; j++) {

            double m13 = (i + 0.5) * (std::sqrt(s13_max) - std::sqrt(s13_min)) / NBin_cal + std::sqrt(s13_min);
            double m23 = (j + 0.5) * (std::sqrt(s23_max) - std::sqrt(s23_min)) / NBin_cal + std::sqrt(s23_min);

            double s13 = m13 * m13;
            double s23 = m23 * m23;
            double s12 = mB * mB + mKp * mKp + mKL0 * mKL0 + mKL0 * mKL0 - s13 - s23;
            double Probability_value = Probability(s13, s23) * ((std::sqrt(s13_max) - std::sqrt(s13_min)) / NBin_cal) * ((std::sqrt(s23_max) - std::sqrt(s23_min)) / NBin_cal);

            // https://halldweb.jlab.org/DocDB/0033/003345/002/dalitz.pdf.
            double DalitzValue = s23 * s12 * s12 + (s23 * s23 + (mKp * mKp - mKL0 * mKL0) * (mB * mB - mKL0 * mKL0) - s23 * (mB * mB + mKL0 * mKL0 + mKL0 * mKL0 + mKp * mKp)) * s12 + ((mKL0 * mKL0 * mB * mB - mKL0 * mKL0 * mKp * mKp) * (mB * mB - mKL0 * mKL0 + mKL0 * mKL0 - mKp * mKp));

            if (DalitzValue < 0 && s23 >= s13) {
                if (isnan(Probability_value)) continue;
                Prob_1D_s23->Fill(m23, 4 * m13 * m23 * Probability_value); // multiply Jacobian
            }

        }
    }

    for (int i = 0; i < NBin_cal; i++) {
        for (int j = 0; j < NBin_cal; j++) {

            double s12 = (i + 0.5) * (s12_max - s12_min) / NBin_cal + s12_min;
            double s13 = (j + 0.5) * (s13_max - s13_min) / NBin_cal + s13_min;
            double s23 = mB * mB + mKp * mKp + mKL0 * mKL0 + mKL0 * mKL0 - s12 - s13;
            double Probability_value = Probability(s13, s23) * ((s12_max - s12_min) / NBin_cal) * ((s13_max - s13_min) / NBin_cal);

            // https://halldweb.jlab.org/DocDB/0033/003345/002/dalitz.pdf.
            double DalitzValue = s23 * s12 * s12 + (s23 * s23 + (mKp * mKp - mKL0 * mKL0) * (mB * mB - mKL0 * mKL0) - s23 * (mB * mB + mKL0 * mKL0 + mKL0 * mKL0 + mKp * mKp)) * s12 + ((mKL0 * mKL0 * mB * mB - mKL0 * mKL0 * mKp * mKp) * (mB * mB - mKL0 * mKL0 + mKL0 * mKL0 - mKp * mKp));

            if (DalitzValue < 0 && s23 >= s13) {
                if (isnan(Probability_value)) continue;
                Prob_conv->Fill(s12, s13, Probability_value);
            }

        }
    }

    for (int i = 0; i < NBin_cal; i++) {
        for (int j = 0; j < NBin_cal; j++) {

            double m12 = (i + 0.5) * (std::sqrt(s12_max) - std::sqrt(s12_min)) / NBin_cal + std::sqrt(s12_min);
            double m13 = (j + 0.5) * (std::sqrt(s13_max) - std::sqrt(s13_min)) / NBin_cal + std::sqrt(s13_min);

            double s12 = m12 * m12;
            double s13 = m13 * m13;
            double s23 = mB * mB + mKp * mKp + mKL0 * mKL0 + mKL0 * mKL0 - s12 - s13;
            double Probability_value = Probability(s13, s23) * ((std::sqrt(s12_max) - std::sqrt(s12_min)) / NBin_cal) * ((std::sqrt(s13_max) - std::sqrt(s13_min)) / NBin_cal);

            // https://halldweb.jlab.org/DocDB/0033/003345/002/dalitz.pdf.
            double DalitzValue = s23 * s12 * s12 + (s23 * s23 + (mKp * mKp - mKL0 * mKL0) * (mB * mB - mKL0 * mKL0) - s23 * (mB * mB + mKL0 * mKL0 + mKL0 * mKL0 + mKp * mKp)) * s12 + ((mKL0 * mKL0 * mB * mB - mKL0 * mKL0 * mKp * mKp) * (mB * mB - mKL0 * mKL0 + mKL0 * mKL0 - mKp * mKp));

            if (DalitzValue < 0 && s23 >= s13) {
                if (isnan(Probability_value)) continue;
                Prob_1D_s12->Fill(m12, 4 * m13 * m12 * Probability_value); // multiply Jacobian
                Prob_1D_s13->Fill(m13, 4 * m13 * m12 * Probability_value); // multiply Jacobian
            }

        }
    }

    double norm = 1.0 / Prob->Integral();
    printf("Initial area of Prob: %.9lf\nnormalized to one\n", 1.0 / norm);
    Prob->Scale(norm);

    double norm_conv = 1.0 / Prob_conv->Integral();
    Prob_conv->Scale(norm_conv);

    double norm_1D;

    norm_1D = 1.0 / Prob_1D_s12->Integral();
    Prob_1D_s12->Scale(norm_1D);

    norm_1D = 1.0 / Prob_1D_s13->Integral();
    Prob_1D_s13->Scale(norm_1D);

    norm_1D = 1.0 / Prob_1D_s23->Integral();
    Prob_1D_s23->Scale(norm_1D);
}

void Draw1DPlots() {

    TCanvas* c = new TCanvas("c", "", 700, 700);

    Prob_1D_s12->SetStats(0);
    Prob_1D_s12->Draw("Hist");
    c->SaveAs("1D_sqrt_s12.png");

    Prob_1D_s13->SetStats(0);
    Prob_1D_s13->Draw("Hist");
    c->SaveAs("1D_sqrt_s13.png");

    Prob_1D_s23->SetStats(0);
    Prob_1D_s23->Draw("Hist");
    c->SaveAs("1D_sqrt_s23.png");

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

            N_evt_conv->SetBinContent(i + 1, j + 1, 0.0);
            Prob_conv->SetBinContent(i + 1, j + 1, 0.0);
            Prob_PHSP_conv->SetBinContent(i + 1, j + 1, 0.0);
        }
    }

    for (size_t i = 0; i < Infos.size(); i++) {
        double s12_ = Infos.at(i).s12;
        double s13_ = Infos.at(i).s13;
        double s23_ = Infos.at(i).s23;
        double m_ = Infos.at(i).m;
        double pstar_mag_ = Infos.at(i).pstar_mag;
        double p_mag_ = Infos.at(i).p_mag;
        double q_mag_ = Infos.at(i).q_mag;
        double p_dot_q_ = Infos.at(i).p_dot_q;

        int GLobalBin = 0;

        // find Bin index
        GLobalBin = N_evt->FindBin(s13_, s23_);
        N_evt->SetBinContent(GLobalBin, N_evt->GetBinContent(GLobalBin) + 1.0);

        // find Bin index
        GLobalBin = Prob_PHSP->FindBin(s13_, s23_);
        Prob_PHSP->SetBinContent(GLobalBin, Prob_PHSP->GetBinContent(GLobalBin) + 1.0);

        // find Bin index in s12-s13 plane
        GLobalBin = N_evt_conv->FindBin(s12_, s13_);
        N_evt_conv->SetBinContent(GLobalBin, N_evt_conv->GetBinContent(GLobalBin) + 1.0);

        // find Bin index in s12-s13 plane
        GLobalBin = Prob_PHSP_conv->FindBin(s12_, s13_);
        Prob_PHSP_conv->SetBinContent(GLobalBin, Prob_PHSP_conv->GetBinContent(GLobalBin) + 1.0);
    }

    Prob_PHSP->Scale(1.0 / Prob_PHSP->Integral());
    Prob_PHSP_conv->Scale(1.0 / Prob_PHSP_conv->Integral());
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

        double p_mag = -1;
        double Ep = -1;
        double px = -1;
        double py = -1;
        double pz = -1;

        double q_mag = -1;
        double Eq = -1;
        double qx = -1;
        double qy = -1;
        double qz = -1;

        double s12 = -1;
        double s_first = -1;
        double s_second = -1;

        double pstar_mag = -1;

        // get event_info
        tree_info->SetBranchAddress("__experiment__", &__experiment__);
        tree_info->SetBranchAddress("__run__", &__run__);
        tree_info->SetBranchAddress("__event__", &__event__);
        tree_info->SetBranchAddress("__candidate__", &__candidate__);
        tree_info->SetBranchAddress("__ncandidates__", &__ncandidates__);

        tree_info->SetBranchAddress("invMassInLists__boK_L0__clKpKLKL__bc", &m);

        tree_info->SetBranchAddress("useDaughterRestFrame__bodaughter__bo0__cm__spmcP__bc__cm__sp1__cm__sp2__bc", &p_mag);
        tree_info->SetBranchAddress("useDaughterRestFrame__bodaughter__bo0__cm__spmcE__bc__cm__sp1__cm__sp2__bc", &Ep);
        tree_info->SetBranchAddress("useDaughterRestFrame__bodaughter__bo0__cm__spmcPX__bc__cm__sp1__cm__sp2__bc", &px);
        tree_info->SetBranchAddress("useDaughterRestFrame__bodaughter__bo0__cm__spmcPY__bc__cm__sp1__cm__sp2__bc", &py);
        tree_info->SetBranchAddress("useDaughterRestFrame__bodaughter__bo0__cm__spmcPZ__bc__cm__sp1__cm__sp2__bc", &pz);

        tree_info->SetBranchAddress("useDaughterRestFrame__bodaughter__bo1__cm__spmcP__bc__cm__sp1__cm__sp2__bc", &q_mag);
        tree_info->SetBranchAddress("useDaughterRestFrame__bodaughter__bo1__cm__spmcE__bc__cm__sp1__cm__sp2__bc", &Eq);
        tree_info->SetBranchAddress("useDaughterRestFrame__bodaughter__bo1__cm__spmcPX__bc__cm__sp1__cm__sp2__bc", &qx);
        tree_info->SetBranchAddress("useDaughterRestFrame__bodaughter__bo1__cm__spmcPY__bc__cm__sp1__cm__sp2__bc", &qy);
        tree_info->SetBranchAddress("useDaughterRestFrame__bodaughter__bo1__cm__spmcPZ__bc__cm__sp1__cm__sp2__bc", &qz);

        tree_info->SetBranchAddress("daughterInvariantMass__bo1__cm__sp2__bc", &s12);
        tree_info->SetBranchAddress("daughterInvariantMass__bo0__cm__sp1__bc", &s_first);
        tree_info->SetBranchAddress("daughterInvariantMass__bo0__cm__sp2__bc", &s_second);

        tree_info->SetBranchAddress("daughter__bo0__cm__spuseParticleRestFrame__bomcP__cm__spB__pl__clKpKLKL_BCS__bc__bc", &pstar_mag);

        printf("%lld entries...\n", tree_info->GetEntries());
        for (unsigned int j = 0; j < tree_info->GetEntries(); j++) { // Fill

            tree_info->GetEntry(j);

            double s12_ = s12 * s12;
            double s13_ = std::min(s_first, s_second) * std::min(s_first, s_second);
            double s23_ = std::max(s_first, s_second) * std::max(s_first, s_second);
            double m_ = m;
            double pstar_mag_ = pstar_mag;
            double p_mag_ = p_mag;
            double q_mag_ = q_mag;
            double p_dot_q_ = px * qx + py * qy + pz * qz;

            // Makeshift to remove photon production decay
            if (std::abs(s12_ -(mB * mB + mKp * mKp + mKL0 * mKL0 + mKL0 * mKL0 - s13_ - s23_)) > MyEPSILON) continue;

            Information temp_info = { s12_, s13_, s23_, m_, pstar_mag_, p_mag_, q_mag_, p_dot_q_ };

            Infos.push_back(temp_info);

        }
        input_file->Close();
    }
}

void GetWeights() {
    FILE* fp;

    fp = fopen("KpKLKL_weight.txt","w");

    fprintf(fp, "s13: %d %lf %lf\n", NBin, s13_min, s13_max);
    fprintf(fp, "s23: %d %lf %lf\n", NBin, s23_min, s23_max);
    for (int i = 0; i < NBin; i++) {
        for (int j = 0; j < NBin; j++) {
            double s13 = (i + 0.5) * (s13_max - s13_min) / NBin + s13_min;
            double s23 = (j + 0.5) * (s23_max - s23_min) / NBin + s23_min;

            int GLobalBin_PHSP = Prob_PHSP->FindBin(s13, s23);
            double PHSP_val = Prob_PHSP->GetBinContent(GLobalBin_PHSP);

            int GLobalBin_model = Prob->FindBin(s13, s23);
            double model_val = Prob->GetBinContent(GLobalBin_model);

            int GLobalBin_Nevt = N_evt->FindBin(s13, s23);
            double Nevt = N_evt->GetBinContent(GLobalBin_Nevt);

            if (Nevt < MyEPSILON) fprintf(fp, "%lf %lf %lf\n", s13, s23, 0.0);
            else fprintf(fp, "%lf %lf %lf\n", s13, s23, model_val / PHSP_val);
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

    DrawDalitz(N_evt, "PHSP.png");

    GetProbabilityHist();

    DrawDalitz(Prob, "model.png");
    DrawDalitz(Prob_conv, "model_conv.png");

    Draw1DPlots();

    GetWeights();

    return 0;
}
