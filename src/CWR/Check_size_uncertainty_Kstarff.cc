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
    std::vector<double> A1_values;
    std::vector<double> A12_values;
    std::vector<double> v0_values;

    const double m_b = 5.27934; // B meson mass
    const double m_k = 0.89167; // Kaon mass

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
    {0.081726107544427,   0.418943528790113, -0.148580744618415,   0.008936060075137,   0.016121321335068, -0.134221136343898,   0.101230277583689,   0.875232126303976,   0.033121258010682},
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
    double fluctuations[19][9] = { 0.0 };
    for (int i = 0; i < 19; i++) {
        for (int j = 0; j < 9; j++) {
            for (int k = 0; k < 9; k++) fluctuations[i][j] = fluctuations[i][j] + Lambdas[i][k] * LinearCoefficients[k][j];
        }
    }


    for (int i = 0; i < NToys; i++) {

        double values[9] = { 0.0 };
        for (int j = 0; j < 9; j++) values[j] = normal(gen);

        double alpha0_A1_fluc = alpha0_A1;
        double alpha1_A1_fluc = alpha1_A1;
        double alpha2_A1_fluc = alpha2_A1;
        double alpha0_A12_fluc = alpha0_A12;
        double alpha1_A12_fluc = alpha1_A12;
        double alpha2_A12_fluc = alpha2_A12;
        double alpha0_v0_fluc = alpha0_v0;
        double alpha1_v0_fluc = alpha1_v0;
        double alpha2_v0_fluc = alpha2_v0;

        for (int j = 0; j < 9; j++) {
            alpha0_A1_fluc = alpha0_A1_fluc + values[j] * fluctuations[2 * j + 1][0];
            alpha1_A1_fluc = alpha1_A1_fluc + values[j] * fluctuations[2 * j + 1][1];
            alpha2_A1_fluc = alpha2_A1_fluc + values[j] * fluctuations[2 * j + 1][2];
            alpha0_A12_fluc = alpha0_A12_fluc + values[j] * fluctuations[2 * j + 1][3];
            alpha1_A12_fluc = alpha1_A12_fluc + values[j] * fluctuations[2 * j + 1][4];
            alpha2_A12_fluc = alpha2_A12_fluc + values[j] * fluctuations[2 * j + 1][5];
            alpha0_v0_fluc = alpha0_v0_fluc + values[j] * fluctuations[2 * j + 1][6];
            alpha1_v0_fluc = alpha1_v0_fluc + values[j] * fluctuations[2 * j + 1][7];
            alpha2_v0_fluc = alpha2_v0_fluc + values[j] * fluctuations[2 * j + 1][8];
        }

        double tp = (m_b + m_k) * (m_b + m_k);
        double tm = (m_b - m_k) * (m_b - m_k);
        double t0 = tp * (1 - sqrt(1 - tm / tp));
        double z = (sqrt(tp - q2) - sqrt(tp - t0)) / (sqrt(tp - q2) + sqrt(tp - t0));
        double z0 = (sqrt(tp) - sqrt(tp - t0)) / (sqrt(tp) + sqrt(tp - t0));

        double v0 = (1 / (1 - q2 / (mR_v0 * mR_v0))) * (alpha0_v0_fluc + alpha1_v0_fluc * (z - z0) + alpha2_v0_fluc * (z - z0) * (z - z0));
        double A1 = (1 / (1 - q2 / (mR_A1 * mR_A1))) * (alpha0_A1_fluc + alpha1_A1_fluc * (z - z0) + alpha2_A1_fluc * (z - z0) * (z - z0));
        double A12 = (1 / (1 - q2 / (mR_A12 * mR_A12))) * (alpha0_A12_fluc + alpha1_A12_fluc * (z - z0) + alpha2_A12_fluc * (z - z0) * (z - z0));

        A1_values.push_back(A1);
        A12_values.push_back(A12);
        v0_values.push_back(v0);
    }

    // get distribution
    auto min_it_A1 = std::min_element(A1_values.begin(), A1_values.end());
    auto max_it_A1 = std::max_element(A1_values.begin(), A1_values.end());

    auto min_it_A12 = std::min_element(A12_values.begin(), A12_values.end());
    auto max_it_A12 = std::max_element(A12_values.begin(), A12_values.end());

    auto min_it_v0 = std::min_element(v0_values.begin(), v0_values.end());
    auto max_it_v0 = std::max_element(v0_values.begin(), v0_values.end());

    TH1D* dist_A1 = new TH1D("dist_A1", ";A1;arbitrary unit", 100, *min_it_A1, *max_it_A1);
    for (int i = 0; i < A1_values.size(); i++) dist_A1->Fill(A1_values.at(i));

    TH1D* dist_A12 = new TH1D("dist_A12", ";A12;arbitrary unit", 100, *min_it_A12, *max_it_A12);
    for (int i = 0; i < A12_values.size(); i++) dist_A12->Fill(A12_values.at(i));

    TH1D* dist_v0 = new TH1D("dist_v0", ";v0;arbitrary unit", 100, *min_it_v0, *max_it_v0);
    for (int i = 0; i < v0_values.size(); i++) dist_v0->Fill(v0_values.at(i));

    TCanvas* c_temp = new TCanvas("c", "", 800, 800);
    dist_A1->Draw("Hist");
    c_temp->SaveAs("A1_dist.png");

    dist_A12->Draw("Hist");
    c_temp->SaveAs("A12_dist.png");

    dist_v0->Draw("Hist");
    c_temp->SaveAs("v0_dist.png");

    return 0;
}
