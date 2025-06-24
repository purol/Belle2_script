// last update: 2023-07-14
// for Belle2 data

#include <stdio.h>
#include <vector>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <algorithm>

#include <TH1D.h>
#include <TCanvas.h>

# define NToys 10000

int main(int argc, char* argv[]) {

    std::random_device rd{};
    std::mt19937 gen{ rd() };
    std::normal_distribution normal{ 0.0, 1.0 };

    double q2 = 0.5;
    std::vector<double> fp_values;

    const double m_b = 5.27934; // B meson mass
    const double m_k = 0.493677; // Kaon mass

    const double alpha0 = 0.2545; // 0.2545 +- 0.0090
    const double alpha1 = -0.71; // -0.71 +- 0.14
    const double alpha2 = 0.32; // 0.32 +- 0.59

    /*
    <PCA>
    e0 = 0.00090 alpha0 + 0.19987 alpha1 + 0.97982 alpha2
    e1 = 0.06182 alpha0 + 0.97794 alpha1 - 0.19954 alpha2
    e2 = 0.99809 alpha0 - 0.06075 alpha2 + 0.01148 alpha2
    lambda0 = 0.60197
    lambda1 = 0.07314
    lambda2 = 0.00778
    */
    const double LinearCoefficients[3][3] = {
        {0.000898287544607,   0.199871550874459,   0.979821696141969},
        {0.061817247836003,   0.977937066768286,  -0.199543782941545},
        {0.998087080844888,  -0.060749128320321,   0.011477040509285}
    };
    const double Lambdas[7][3] = {
        { 0.0,                0.0,                0.0              },
        { 0.601972208713257,  0.0,                0.0              },
        {-0.601972208713257,  0.0,                0.0              },
        { 0.0,                0.073143483537865,  0.0              },
        { 0.0,               -0.073143483537865,  0.0              },
        { 0.0,                0.0,                0.007777700923050},
        { 0.0,                0.0,               -0.007777700923050}
    };
    double fluctuations[7][3] = { 0.0 };
    for (int i = 0; i < 7; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) fluctuations[i][j] = fluctuations[i][j] + Lambdas[i][k] * LinearCoefficients[k][j];
        }
    }
    /*
    fluctuations[0][0] = 0.0; // alpha0 change
    fluctuations[0][1] = 0.0; // alpha1 change
    fluctuations[0][2] = 0.0; // alpha2 change
    fluctuations[1][0] = 0.60197 * 0.00090; // alpha0 change Lambdas[1][0]*LinearCoefficients[0][0] + Lambdas[1][1]*LinearCoefficients[1][0] + Lambdas[1][2]*LinearCoefficients[2][0]
    fluctuations[1][1] = 0.60197 * 0.19987; // alpha1 change Lambdas[1][0]*LinearCoefficients[0][1] + Lambdas[1][1]*LinearCoefficients[1][1] + Lambdas[1][2]*LinearCoefficients[2][1]
    fluctuations[1][2] = 0.60197 * 0.97982; // alpha2 change
    fluctuations[2][0] = (-0.60197) * 0.00090; // alpha0 change
    fluctuations[2][1] = (-0.60197) * 0.19987; // alpha1 change
    fluctuations[2][2] = (-0.60197) * 0.97982; // alpha2 change
    fluctuations[3][0] = 0.07314 * 0.06182; // alpha0 change
    fluctuations[3][1] = 0.07314 * 0.97794; // alpha1 change
    fluctuations[3][2] = 0.07314 * (-0.19954); // alpha2 change
    fluctuations[4][0] // alpha0 change
    fluctuations[4][1] // alpha1 change
    fluctuations[4][2] // alpha2 change
    fluctuations[5][0] // alpha0 change
    fluctuations[5][1] // alpha1 change
    fluctuations[5][2] // alpha2 change
    fluctuations[6][0] // alpha0 change
    fluctuations[6][1] // alpha1 change
    fluctuations[6][2] // alpha2 change
    */

    for (int i = 0; i < NToys; i++) {

        double values[3] = { 0.0 };
        for(int j = 0; j < 3; j++) values[j] = normal(gen);

        double alpha0_fluc = alpha0;
        double alpha1_fluc = alpha1;
        double alpha2_fluc = alpha2;

        for (int j = 0; j < 3; j++) {
            alpha0_fluc = alpha0_fluc + values[j] * fluctuations[2 * j + 1][0];
            alpha1_fluc = alpha1_fluc + values[j] * fluctuations[2 * j + 1][1];
            alpha2_fluc = alpha2_fluc + values[j] * fluctuations[2 * j + 1][2];
        }

        double mp = 5.4158;
        double tp = (m_b + m_k) * (m_b + m_k);
        double tm = (m_b - m_k) * (m_b - m_k);
        double t0 = 0.0;
        double z = (sqrt(tp - q2) - sqrt(tp - t0)) / (sqrt(tp - q2) + sqrt(tp - t0));
        double fp = (1 / (1 - q2 / (mp * mp))) * (alpha0_fluc + alpha1_fluc * z + alpha2_fluc * z * z + (-alpha1_fluc + 2 * alpha2_fluc) * z * z * z / 3);

        fp_values.push_back(fp);
    }

    // get distribution
    auto min_it = std::min_element(fp_values.begin(), fp_values.end());
    auto max_it = std::max_element(fp_values.begin(), fp_values.end());

    TH1D* dist_fp = new TH1D("dist_fp", ";fp;arbitrary unit", 100, *min_it, *max_it);
    for (int i = 0; i < fp_values.size(); i++) dist_fp->Fill(fp_values.at(i));

    TCanvas* c_temp = new TCanvas("c", "", 800, 800);
    dist_fp->Draw("Hist");
    c_temp->SaveAs("fp_dist.png");

    return 0;
}
