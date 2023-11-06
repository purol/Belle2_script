// last update: 2021-10-13
// for Belle2 data

# define N_decay 38 // five decay mode + others

# define Nentry 200
# define MIN_Q2 0.0
# define MAX_Q2 17.64

# define Nvar 23
# define DvetoNvar 4

# define MyEPSILON 0.000001
# define MCTYPE "MC15ri"

class Corrector_Fragmentation {
private:
    static constexpr int N_Category_gamma = 9;
    static constexpr int N_Bin_gamma = 5;
    static constexpr double Bins_gamma[N_Bin_gamma - 1] = { 1.15, 1.5, 2.0, 2.4 };
    static constexpr double Fragmentation_Xsgamma[N_Bin_gamma][N_Category_gamma] = { // unit: percentage with respect to "un-missing" mode
        {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00},
        {10.6, 5.88, 23.2, 44.5, 0.46, 9.99, 0.52, 4.78, 0.00},
        {3.12, 1.13, 15.7, 20.6, 9.48, 26.9, 5.29, 15.5, 2.20},
        {1.65, 0.51, 9.02, 7.98, 16.8, 26.7, 10.3, 21.8, 5.38},
        {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}
    };
    static constexpr double Fragmentation_Uncertainty_Xsgamma[N_Bin_gamma][N_Category_gamma] = { // relative uncertainty
        {0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000},
        {0.0651, 0.0990, 0.0554, 0.0400, 0.6304, 0.1175, 0.7115, 0.2395, 0.0000},
        {0.1053, 0.2727, 0.0686, 0.0741, 0.0942, 0.0733, 0.2777, 0.1715, 0.1168},
        {0.4250, 1.4490, 0.2646, 0.4003, 0.2509, 0.2671, 0.8409, 0.6682, 0.2898},
        {0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000}
    };

    static constexpr double Nevt_Nominal_before_Xsgamma_MC15[N_Bin_gamma][N_Category_gamma] = {
        {228.010067, 113.889034, 97.106455, 99.597443, 2.350444, 17.694690, 0.265738, 25.726503, 0.000112}, // [-inf, 1.15]
        {1217.857056, 608.996474, 683.434991, 821.802920, 54.319919, 320.575746, 29.362053, 194.574254, 0.829998}, // [1.15, 1.5]
        {745.156145, 372.620424, 554.609584, 717.744130, 132.646240, 507.837684, 190.020473, 312.058874, 120.779580}, // [1.5, 2.0]
        {125.738141, 62.841816, 119.699133, 149.610098, 58.094970, 157.286832, 123.180903, 139.932048, 49.973636}, // [2.0, 2.4]
        {22.820347, 11.451874, 24.859499, 30.438579, 15.176837, 37.050655, 37.906534, 38.685844, 11.884932} // [2.4, inf]
    };
    static constexpr double Nevt_Hmb_before_Xsgamma_MC15[N_Bin_gamma][N_Category_gamma] = {
        {234.717370, 116.931851, 100.045961, 102.811438, 2.443726, 18.213362, 0.283348, 26.352459, 0.000000}, // [-inf, 1.15]
        {1237.919103, 618.966253, 694.727432, 833.619261, 54.873905, 324.784821, 29.598600, 196.698641, 0.862629}, // [1.15, 1.5]
        {739.765455, 368.897914, 550.204743, 712.129243, 131.851781, 504.040211, 187.911583, 309.873035, 119.610335}, // [1.5, 2.0]
        {122.975491, 61.761023, 116.640122, 146.325654, 56.851972, 153.571830, 120.538083, 137.029154, 49.297575}, // [2.0, 2.4]
        {22.077625, 11.336836, 24.298764, 29.398639, 14.646905, 35.886317, 36.466259, 37.697931, 11.619447} // [2.4, inf]
    };
    static constexpr double Nevt_Hpf_before_Xsgamma_MC15[N_Bin_gamma][N_Category_gamma] = {
        {211.414635, 105.941646, 89.984692, 92.845350, 2.132717, 16.210070, 0.250627, 23.805551, 0.000569}, // [-inf, 1.15]
        {1145.424694, 573.572296, 643.905532, 775.309971, 51.161825, 303.573493, 27.884826, 183.335197, 0.805736}, // [1.15, 1.5]
        {736.485305, 369.574161, 550.673698, 711.747889, 132.674662, 505.825855, 190.843751, 312.993403, 121.016612}, // [1.5, 2.0]
        {137.460490, 68.603354, 131.880622, 163.688524, 63.757702, 172.478521, 135.887581, 153.778110, 54.689402}, // [2.0, 2.4]
        {29.109281, 14.773702, 32.144158, 39.024643, 19.707985, 47.791272, 49.197184, 50.159756, 15.393167} // [2.4, inf]
    };
    static constexpr double Nevt_Htransition_before_Xsgamma_MC15[N_Bin_gamma][N_Category_gamma] = {
        {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, // [-inf, 1.15]
        {1126.011312, 562.651975, 658.733369, 809.726662, 55.294719, 337.131227, 31.789385, 187.359070, 0.920546}, // [1.15, 1.5]
        {838.951829, 419.735045, 624.498178, 809.052632, 149.312766, 571.803869, 214.105973, 351.790571, 136.066782}, // [1.5, 2.0]
        {141.771711, 71.022628, 135.144794, 168.801789, 65.504295, 176.864279, 139.131693, 157.662450, 56.155419}, // [2.0, 2.4]
        {25.757474, 12.884694, 27.851586, 34.061168, 16.946366, 41.803329, 42.925015, 43.435334, 13.297957} // [2.4, inf]
    };
    static constexpr double Nevt_Lmb_before_Xsgamma_MC15[N_Bin_gamma][N_Category_gamma] = {
        {219.785163, 110.184748, 93.899905, 96.479643, 2.274163, 17.081852, 0.246506, 24.692561, 0.000000}, // [-inf, 1.15]
        {1195.125885, 599.131943, 672.175942, 810.139571, 53.494062, 316.373461, 29.121914, 191.567680, 0.804038}, // [1.15, 1.5]
        {750.803012, 374.158691, 558.976621, 723.091791, 134.075592, 512.438695, 192.030343, 314.970543, 122.127927}, // [1.5, 2.0]
        {128.850611, 64.181819, 121.971830, 153.567926, 59.265244, 160.822175, 126.599064, 143.704766, 51.333129}, // [2.0, 2.4]
        {23.591953, 11.740537, 25.576218, 31.342805, 15.603211, 38.154400, 39.437907, 39.746047, 12.199234} // [2.4, inf]
    };
    static constexpr double Nevt_Lpf_before_Xsgamma_MC15[N_Bin_gamma][N_Category_gamma] = {
        {250.065120, 125.068681, 106.236894, 109.408332, 2.586983, 19.616770, 0.297757, 28.249527, 0.000000}, // [-inf, 1.15]
        {1307.816445, 653.786860, 732.097197, 879.672426, 58.091454, 342.326888, 31.063100, 208.314601, 0.883302}, // [1.15, 1.5]
        {747.709095, 373.918779, 554.600603, 718.457434, 131.129261, 504.948261, 186.036870, 308.211200, 119.650476}, // [1.5, 2.0]
        {109.711774, 55.152876, 104.872451, 130.564821, 50.850965, 137.527349, 107.684824, 121.713044, 43.870882}, // [2.0, 2.4]
        {16.301694, 8.096028, 17.764354, 21.644725, 10.719722, 26.442162, 26.669031, 27.403075, 8.479496} // [2.4, inf]
    };
    static constexpr double Nevt_Ltransition_before_Xsgamma_MC15[N_Bin_gamma][N_Category_gamma] = {
        {687.014043, 343.100287, 227.441995, 234.500737, 5.644878, 35.769345, 0.466491, 61.879282, 0.000000}, // [-inf, 1.15]
        {1095.635246, 549.051012, 615.383772, 741.741024, 48.963482, 288.728398, 26.536985, 175.580393, 0.736949}, // [1.15, 1.5]
        {671.592637, 336.051687, 499.387200, 646.831974, 119.008094, 456.176166, 171.088734, 281.417112, 108.985261}, // [1.5, 2.0]
        {113.378578, 56.429701, 107.697531, 134.644214, 52.336623, 141.679228, 110.495689, 126.585786, 45.203748}, // [2.0, 2.4]
        {20.750444, 10.246957, 22.414418, 27.261420, 13.458927, 33.490579, 34.139560, 34.988253, 10.795822} // [2.4, inf]
    };

    static constexpr double Nevt_Nominal_missing_before_Xsgamma_MC15[N_Bin_gamma] = { 48.015884, 360.363091, 676.491505, 518.892840, 292.939844 }; // number of missing mode (MC)
    static constexpr double Nevt_Hmb_missing_before_Xsgamma_MC15[N_Bin_gamma] = { 49.450546, 365.252264, 669.873825, 507.353413, 284.711591 }; // number of missing mode (MC)
    static constexpr double Nevt_Hpf_missing_before_Xsgamma_MC15[N_Bin_gamma] = { 44.686813, 340.083522, 676.515003, 574.044827, 388.939028 }; // number of missing mode (MC)
    static constexpr double Nevt_Htransition_missing_before_Xsgamma_MC15[N_Bin_gamma] = { 0.000000, 348.935210, 762.920563, 585.338309, 330.039376 }; // number of missing mode (MC)
    static constexpr double Nevt_Lmb_missing_before_Xsgamma_MC15[N_Bin_gamma] = { 46.502106, 354.423602, 683.769325, 533.204382, 302.346798 }; // number of missing mode (MC)
    static constexpr double Nevt_Lpf_missing_before_Xsgamma_MC15[N_Bin_gamma] = { 52.743632, 386.345000, 669.316781, 448.947976, 200.139315 }; // number of missing mode (MC)
    static constexpr double Nevt_Ltransition_missing_before_Xsgamma_MC15[N_Bin_gamma] = { 93.606503, 323.683107, 609.133190, 468.259061, 263.808865 }; // number of missing mode (MC)

    double Total_Nevt_Nominal_before_Xsgamma_MC15[N_Bin_gamma]; // without missing mode
    double Total_Nevt_Hmb_before_Xsgamma_MC15[N_Bin_gamma];
    double Total_Nevt_Hpf_before_Xsgamma_MC15[N_Bin_gamma];
    double Total_Nevt_Htransition_before_Xsgamma_MC15[N_Bin_gamma];
    double Total_Nevt_Lmb_before_Xsgamma_MC15[N_Bin_gamma];
    double Total_Nevt_Lpf_before_Xsgamma_MC15[N_Bin_gamma];
    double Total_Nevt_Ltransition_before_Xsgamma_MC15[N_Bin_gamma];
public:
    enum class SystType { // what kind of Xs sample is corrected?
        Nominal,
        Hmb,
        Hpf,
        Htransition,
        Lmb,
        Lpf,
        Ltransition
    };
    enum class Sample {
        gamma // B -> Xs gamma is used
    };

    Corrector_Fragmentation();
    int Classify(int Decay[N_decay], Sample sample);
    int GetMxBin(double MXs, Sample sample);
    double GetCorrectionFactor(int Decay[N_decay], double MXs, SystType systtype, Sample sample, std::string type);
    double FluctuateCorrection(int Decay[N_decay], double MXs, SystType systtype, int TargetMxsBin, int TargetCategory, bool IsTargetCategoryUp, Sample sample, std::string type);
    int GetNMxsBin(Sample sample);
    int GetNCategory(Sample sample);
};

Corrector_Fragmentation corrector_Fragmentation;

Corrector_Fragmentation::Corrector_Fragmentation() {
    for (int i = 0; i < N_Bin_gamma; i++) {
        Total_Nevt_Nominal_before_Xsgamma_MC15[i] = 0.0;
        Total_Nevt_Hmb_before_Xsgamma_MC15[i] = 0.0;
        Total_Nevt_Hpf_before_Xsgamma_MC15[i] = 0.0;
        Total_Nevt_Htransition_before_Xsgamma_MC15[i] = 0.0;
        Total_Nevt_Lmb_before_Xsgamma_MC15[i] = 0.0;
        Total_Nevt_Lpf_before_Xsgamma_MC15[i] = 0.0;
        Total_Nevt_Ltransition_before_Xsgamma_MC15[i] = 0.0;
    }

    for (int i = 0; i < N_Bin_gamma; i++) {
        for (int j = 0; j < N_Category_gamma; j++) {
            Total_Nevt_Nominal_before_Xsgamma_MC15[i] = Total_Nevt_Nominal_before_Xsgamma_MC15[i] + Nevt_Nominal_before_Xsgamma_MC15[i][j];
            Total_Nevt_Hmb_before_Xsgamma_MC15[i] = Total_Nevt_Hmb_before_Xsgamma_MC15[i] + Nevt_Hmb_before_Xsgamma_MC15[i][j];
            Total_Nevt_Hpf_before_Xsgamma_MC15[i] = Total_Nevt_Hpf_before_Xsgamma_MC15[i] + Nevt_Hpf_before_Xsgamma_MC15[i][j];
            Total_Nevt_Htransition_before_Xsgamma_MC15[i] = Total_Nevt_Htransition_before_Xsgamma_MC15[i] + Nevt_Htransition_before_Xsgamma_MC15[i][j];
            Total_Nevt_Lmb_before_Xsgamma_MC15[i] = Total_Nevt_Lmb_before_Xsgamma_MC15[i] + Nevt_Lmb_before_Xsgamma_MC15[i][j];
            Total_Nevt_Lpf_before_Xsgamma_MC15[i] = Total_Nevt_Lpf_before_Xsgamma_MC15[i] + Nevt_Lpf_before_Xsgamma_MC15[i][j];
            Total_Nevt_Ltransition_before_Xsgamma_MC15[i] = Total_Nevt_Ltransition_before_Xsgamma_MC15[i] + Nevt_Ltransition_before_Xsgamma_MC15[i][j];
        }
    }
}

int Corrector_Fragmentation::Classify(int Decay[N_decay], Sample sample) {
    if (sample == Corrector_Fragmentation::Sample::gamma) {
        /*
            0: Xs2Kpi_wopi0_MC // Xs -> K pi w/o pi0
            1: Xs2Kpi_wpi0_MC, // Xs -> K pi w/ 1pi0
            2: Xs2Kpipi_wopi0_MC, // Xs -> K pi pi w/o pi0
            3: Xs2Kpipi_wpi0_MC, // Xs -> K pi pi w/ 1pi0
            4: Xs2Kpipipi_wopi0_MC, // Xs -> K pi pi pi w/o pi0
            5: Xs2Kpipipi_wpi0_MC, // Xs -> K pi pi pi w/ 1pi0
            6: Xs2Kpipipipi_atmost1pi0_MC, // Xs -> K pi pi pi pi at most 1pi0
            7: Xs2pi0pi0_MC, // including 2pi0
            8: Xs2KKK_MC, // including 3K
            9: MC_other, // other non-resonant Xs decay
        */
        if (Decay[6] > 0 || Decay[24] > 0) return 0;
        else if (Decay[5] > 0 || Decay[25] > 0) return 1;
        else if (Decay[7] > 0 || Decay[27] > 0) return 2;
        else if (Decay[8] > 0 || Decay[26] > 0) return 3;
        else if (Decay[10] > 0 || Decay[28] > 0) return 4;
        else if (Decay[9] > 0 || Decay[29] > 0) return 5;
        else if (Decay[11] > 0 || Decay[12] > 0 || Decay[30] > 0 || Decay[31] > 0) return 6;
        else if (Decay[13] > 0 || Decay[14] > 0 || Decay[15] > 0 || Decay[32] > 0 || Decay[33] > 0 || Decay[34] > 0) return 7;
        else if (Decay[16] > 0 || Decay[17] > 0 || Decay[18] > 0 || Decay[35] > 0 || Decay[36] > 0 || Decay[37] > 0) return 8;
        else return N_Category_gamma;
    }
    else {
        printf("[Corrector_Fragmentation::Classify] cannot find sample");
        exit(1);
        return 1;
    }
}

int Corrector_Fragmentation::GetMxBin(double MXs, Sample sample) {
    if (sample == Corrector_Fragmentation::Sample::gamma) { // [-inf, 1.15, 1.5, 2.0, 2.4, inf]
        for (int i = 0; i < N_Bin_gamma - 1; i++) {
            if (MXs < Bins_gamma[i]) return i;
        }
        return N_Bin_gamma - 1;
    }
    else {
        printf("[Corrector_Fragmentation::GetMxBin] cannot find sample");
        exit(1);
        return 1;
    }
}

double Corrector_Fragmentation::GetCorrectionFactor(int Decay[N_decay], double MXs, SystType systtype, Sample sample, std::string type) {
    int Category = Classify(Decay, sample);
    int MxsBin = GetMxBin(MXs, sample);

    if (sample == Corrector_Fragmentation::Sample::gamma) {

        if (Category == N_Category_gamma) return 1.0; // no correction if it is missing mode

        if ((MxsBin == 0) || (MxsBin == N_Bin_gamma - 1)) return 1.0; // no correction if Mxs bin is out of range

        if (type == "MC15ri") {
            if (systtype == Corrector_Fragmentation::SystType::Nominal) {
                if (Nevt_Nominal_before_Xsgamma_MC15[MxsBin][Category] < MyEPSILON) return 1.0; // no correction if Nevt is 0 at MC
                return (0.01 * Fragmentation_Xsgamma[MxsBin][Category]) / (Nevt_Nominal_before_Xsgamma_MC15[MxsBin][Category] / Total_Nevt_Nominal_before_Xsgamma_MC15[MxsBin]);
            }
            else if (systtype == Corrector_Fragmentation::SystType::Hmb) {
                if (Nevt_Hmb_before_Xsgamma_MC15[MxsBin][Category] < MyEPSILON) return 1.0; // no correction if Nevt is 0 at MC
                return (0.01 * Fragmentation_Xsgamma[MxsBin][Category]) / (Nevt_Hmb_before_Xsgamma_MC15[MxsBin][Category] / Total_Nevt_Hmb_before_Xsgamma_MC15[MxsBin]);
            }
            else if (systtype == Corrector_Fragmentation::SystType::Hpf) {
                if (Nevt_Hpf_before_Xsgamma_MC15[MxsBin][Category] < MyEPSILON) return 1.0; // no correction if Nevt is 0 at MC
                return (0.01 * Fragmentation_Xsgamma[MxsBin][Category]) / (Nevt_Hpf_before_Xsgamma_MC15[MxsBin][Category] / Total_Nevt_Hpf_before_Xsgamma_MC15[MxsBin]);
            }
            else if (systtype == Corrector_Fragmentation::SystType::Htransition) {
                if (Nevt_Htransition_before_Xsgamma_MC15[MxsBin][Category] < MyEPSILON) return 1.0; // no correction if Nevt is 0 at MC
                return (0.01 * Fragmentation_Xsgamma[MxsBin][Category]) / (Nevt_Htransition_before_Xsgamma_MC15[MxsBin][Category] / Total_Nevt_Htransition_before_Xsgamma_MC15[MxsBin]);
            }
            else if (systtype == Corrector_Fragmentation::SystType::Lmb) {
                if (Nevt_Lmb_before_Xsgamma_MC15[MxsBin][Category] < MyEPSILON) return 1.0; // no correction if Nevt is 0 at MC
                return (0.01 * Fragmentation_Xsgamma[MxsBin][Category]) / (Nevt_Lmb_before_Xsgamma_MC15[MxsBin][Category] / Total_Nevt_Lmb_before_Xsgamma_MC15[MxsBin]);
            }
            else if (systtype == Corrector_Fragmentation::SystType::Lpf) {
                if (Nevt_Lpf_before_Xsgamma_MC15[MxsBin][Category] < MyEPSILON) return 1.0; // no correction if Nevt is 0 at MC
                return (0.01 * Fragmentation_Xsgamma[MxsBin][Category]) / (Nevt_Lpf_before_Xsgamma_MC15[MxsBin][Category] / Total_Nevt_Lpf_before_Xsgamma_MC15[MxsBin]);
            }
            else if (systtype == Corrector_Fragmentation::SystType::Ltransition) {
                if (Nevt_Ltransition_before_Xsgamma_MC15[MxsBin][Category] < MyEPSILON) return 1.0; // no correction if Nevt is 0 at MC
                return (0.01 * Fragmentation_Xsgamma[MxsBin][Category]) / (Nevt_Ltransition_before_Xsgamma_MC15[MxsBin][Category] / Total_Nevt_Ltransition_before_Xsgamma_MC15[MxsBin]);
            }
        }
        else {
            printf("[Corrector_Fragmentation::GetCorrectionFactor] Invalid type!\n");
            exit(1);
            return 0;
        }
    }
    else {
        printf("[Corrector_Fragmentation::GetCorrectionFactor] cannot find sample");
        exit(1);
        return 1;
    }
}

double Corrector_Fragmentation::FluctuateCorrection(int Decay[N_decay], double MXs, SystType systtype, int TargetMxsBin, int TargetCategory, bool IsTargetCategoryUp, Sample sample, std::string type) {
    /*
    Nevt of TargetCategory in TargetMxsBin is fluctuated.
    we calculate the change of Decay[N_decay] (Category)
    */
    int Category = Classify(Decay, sample);
    int MxsBin = GetMxBin(MXs, sample);

    if (TargetMxsBin != MxsBin) return 1.0; // no fluctuation if it is not in target MXs region

    const double RelativeUncertainty = 1.0;

    if (sample == Corrector_Fragmentation::Sample::gamma) {

        double TotalNevtAtMxsBinWithMissing = Total_Nevt_Nominal_before_Xsgamma_MC15[TargetMxsBin] + Nevt_Nominal_missing_before_Xsgamma_MC15[TargetMxsBin];
        double TargetNevtAtMxsBin;
        if (TargetCategory == N_Category_gamma) TargetNevtAtMxsBin = Nevt_Nominal_missing_before_Xsgamma_MC15[TargetMxsBin]; // if it is missing mode
        else TargetNevtAtMxsBin = Total_Nevt_Nominal_before_Xsgamma_MC15[TargetMxsBin] * Fragmentation_Xsgamma[TargetMxsBin][TargetCategory];
        double TotalNevtAtMxsBinWithMissingWithoutTargetCategory = TotalNevtAtMxsBinWithMissing - TargetNevtAtMxsBin;

        if (Category == TargetCategory) {
            if (IsTargetCategoryUp) return (1.0 + RelativeUncertainty);
            else return (1.0 - RelativeUncertainty);
        }
        else {
            if (IsTargetCategoryUp) {
                double output = (TotalNevtAtMxsBinWithMissingWithoutTargetCategory - TargetNevtAtMxsBin * RelativeUncertainty) / TotalNevtAtMxsBinWithMissingWithoutTargetCategory;
                if (output < 0) {
                    printf("[Corrector_Fragmentation::FluctuateCorrection] minus Nevt!\n");
                    exit(1);
                }
                return output;
            }
            else {
                double output = (TotalNevtAtMxsBinWithMissingWithoutTargetCategory + TargetNevtAtMxsBin * RelativeUncertainty) / TotalNevtAtMxsBinWithMissingWithoutTargetCategory;
                if (output < 0) {
                    printf("[Corrector_Fragmentation::FluctuateCorrection] minus Nevt!\n");
                    exit(1);
                }
                return output;
            }
        }
    }
    else {
        printf("[Corrector_Fragmentation::GetCorrectionFactor] cannot find sample");
        exit(1);
        return 1;
    }
}

int Corrector_Fragmentation::GetNMxsBin(Sample sample) {
    if (sample == Corrector_Fragmentation::Sample::gamma) return N_Bin_gamma;
    else {
        printf("[Corrector_Fragmentation::GetMxBin] cannot find sample");
        exit(1);
        return 1;
    }
}

int Corrector_Fragmentation::GetNCategory(Sample sample) {
    if (sample == Corrector_Fragmentation::Sample::gamma) return (N_Category_gamma + 1); // visible + missing
    else {
        printf("[Corrector_Fragmentation::GetMxBin] cannot find sample");
        exit(1);
        return 1;
    }
}

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

void LetsFillq2(const char* dirname, TH1D* hist_q2, Corrector_Fragmentation::SystType systtype) {

    int Decay[N_decay] = { 0 };
    double Mxs_Bc_MC = -1;
    double Mxs_B0_MC = -1;

    double invM = 0.0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_Xs = (TTree*)input_file->Get("Xs");

        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKcharge_total__bc", &Decay[0]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &Decay[1]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &Decay[2]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCcomb__bc", &Decay[3]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch1__bc", &Decay[4]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch2__bc", &Decay[5]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch3__bc", &Decay[6]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch4__bc", &Decay[7]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch5__bc", &Decay[8]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch6__bc", &Decay[9]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch7__bc", &Decay[10]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch8__bc", &Decay[11]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch9__bc", &Decay[12]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch10__bc", &Decay[13]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch11__bc", &Decay[14]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch12__bc", &Decay[15]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch13__bc", &Decay[16]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch14__bc", &Decay[17]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsu__clMCch15__bc", &Decay[18]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKneutral_total__bc", &Decay[19]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &Decay[20]);
        tree_Xs->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &Decay[21]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCcomb__bc", &Decay[22]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch16__bc", &Decay[23]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch17__bc", &Decay[24]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch18__bc", &Decay[25]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch19__bc", &Decay[26]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch20__bc", &Decay[27]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch21__bc", &Decay[28]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch22__bc", &Decay[29]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch23__bc", &Decay[30]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch24__bc", &Decay[31]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch25__bc", &Decay[32]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch26__bc", &Decay[33]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch27__bc", &Decay[34]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch28__bc", &Decay[35]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch29__bc", &Decay[36]);
        tree_Xs->SetBranchAddress("nParticlesInList__boXsd__clMCch30__bc", &Decay[37]);

        tree_Xs->SetBranchAddress("invMassInLists__bonu_e__clMC_signal__bc", &invM);

        tree_Xs->SetBranchAddress("averageValueInList__boB__pl__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_Bc_MC);
        tree_Xs->SetBranchAddress("averageValueInList__boB0__clMC_signal_total_e__cm__spdaughter__bo0__cm__spM__bc__bc", &Mxs_B0_MC);

        printf("%lld entries...\n", tree_Xs->GetEntries());
        for (unsigned int j = 0; j < tree_Xs->GetEntries(); j++) { // Fill
            tree_Xs->GetEntry(j);

            double Fragmentation_reweight = corrector_Fragmentation.GetCorrectionFactor(Decay, Mxs_Bc_MC, systtype, Corrector_Fragmentation::Sample::gamma, MCTYPE);

            hist_q2->Fill(invM * invM, Fragmentation_reweight);

        }
        input_file->Close();

    }

}

void ReadXsumbFiles_r_sp(){
    const char* Xsu_Hmb_dirname = "./Xsu_Hmb";
    const char* Xsu_Lmb_dirname = "./Xsu_Lmb";
    const char* Xsu_dirname = "./Xsu";

    TH1D* hist_Hmb = new TH1D("Q2_Hmb", ";q^{2} [GeV^{2}/c^{4}]; evt", Nentry, MIN_Q2, MAX_Q2);
    TH1D* hist_Lmb = new TH1D("Q2_Lmb", ";q^{2} [GeV^{2}/c^{4}]; evt", Nentry, MIN_Q2, MAX_Q2);
    TH1D* hist = new TH1D("Q2", ";q^{2} [GeV^{2}/c^{4}]; evt", Nentry, MIN_Q2, MAX_Q2);

    LetsFillq2(Xsu_Hmb_dirname, hist_Hmb, Corrector_Fragmentation::SystType::Hmb);
    LetsFillq2(Xsu_Lmb_dirname, hist_Lmb, Corrector_Fragmentation::SystType::Lmb);
    LetsFillq2(Xsu_dirname, hist, Corrector_Fragmentation::SystType::Nominal);

    Double_t factor = 1.;
    hist_Hmb->Scale(factor / hist_Hmb->Integral(), "width");
    hist_Lmb->Scale(factor / hist_Lmb->Integral(), "width");
    hist->Scale(factor / hist->Integral(), "width");

    hist_Hmb->SetStats(false);

    hist_Hmb->SetLineColor(2);
    hist_Lmb->SetLineColor(3);
    hist->SetLineColor(4);

    double ymax_Hmb = hist_Hmb->GetMaximum();
    double ymax_Lmb = hist_Lmb->GetMaximum();
    double ymax = hist->GetMaximum();
    double real_ymax = 0;
    if (ymax_Hmb > ymax_Lmb) real_ymax = ymax_Hmb;
    else real_ymax = ymax_Lmb;
    if (real_ymax < ymax) real_ymax = ymax;

    TCanvas* c_temp = new TCanvas("c", "", 1200, 1200); c_temp->cd();
    hist_Hmb->SetMaximum(real_ymax * 1.1);
    TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(hist_Hmb,"m_{b} = 4.83 GeV/c^{2}","f");
    legend->AddEntry(hist_Lmb, "m_{b} = 4.53 GeV/c^{2}", "f");
    legend->AddEntry(hist, "m_{b} = 4.68 GeV/c^{2}", "f");
    legend->SetFillStyle(0);
    legend->SetLineWidth(0);
    hist_Hmb->Draw("Hist");
    hist_Lmb->Draw("SAME Hist");
    hist->Draw("SAME Hist");
    legend->Draw("SAME");
    c_temp->SaveAs("Xsu_mb.png");

    FILE* fp;
    fp = fopen("Xsu_Hmb_weight.txt","w");
    fprintf(fp,"%d %lf %lf\n", Nentry, MIN_Q2, MAX_Q2);
    for (int j = 0; j < Nentry; j++) {
        double value = hist->GetBinContent(j + 1);
        double Hmb_value = hist_Hmb->GetBinContent(j + 1);

        double weight = 0;
        if (value > 0) weight = Hmb_value / value;
        else weight = 0;

        fprintf(fp, "%lf\n", weight);
    }
    fclose(fp);

    fp = fopen("Xsu_Lmb_weight.txt", "w");
    fprintf(fp, "%d %lf %lf\n", Nentry, MIN_Q2, MAX_Q2);
    for (int j = 0; j < Nentry; j++) {
        double value = hist->GetBinContent(j + 1);
        double Lmb_value = hist_Lmb->GetBinContent(j + 1);

        double weight = 0;
        if (value > 0) weight = Lmb_value / value;
        else weight = 0;

        fprintf(fp, "%lf\n", weight);
    }
}
