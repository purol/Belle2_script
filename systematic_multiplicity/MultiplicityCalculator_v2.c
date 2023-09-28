// last update: 2021-10-13
// for Belle2 data

# define NgammaMAX 20

# define MyEPSILON 0.000001

// arXiv:1409.4557v2, PhysRevD.107.014511
# define TB0 1.5195 // (Table. 1)
# define TBp 1.6384 // (Table. 1)
# define BR_Kplus_nunubar 0.000005044 // Table VI = (5.044 +- 0.402) * 10^{-6}
# define BR_K0star_nunubar 0.00000919 // (eq. 11)
# define BR_K0_nunubar (BR_Kplus_nunubar*TB0/TBp) // under (eq. 15). In Table VI = (4.6669 +- 0.3707) * 10^{-6}
# define BR_Kplusstar_nunubar (BR_K0star_nunubar*TBp/TB0) // under (eq. 15)
# define BR_Xs_nunubar 0.000029 // (eq. 23)
# define BR_Xsu_nonresonant_nunubar (BR_Xs_nunubar - BR_Kplus_nunubar - BR_Kplusstar_nunubar)
# define BR_Xsd_nonresonant_nunubar (BR_Xs_nunubar - BR_K0_nunubar - BR_K0star_nunubar)

// according to DIRAC
# define N_BpBp_1invab 540000000.0
# define N_B0B0_1invab 510000000.0

# define BR_BpBp 0.514
# define BR_B0B0 0.486

// https://confluence.desy.de/pages/viewpage.action?spaceKey=BI&title=Conference+readiness
# define N_BB_LS1 387100000.0 // NBB = (387.1 +/- 5.6) x 10^6

# define N_Kplus_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_BpBp/(BR_BpBp+BR_B0B0)) * BR_Kplus_nunubar)
# define N_Kplusstar_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_BpBp/(BR_BpBp+BR_B0B0)) * BR_Kplusstar_nunubar)
# define N_Xsu_nonresonant_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_BpBp/(BR_BpBp+BR_B0B0)) * BR_Xsu_nonresonant_nunubar)
# define N_K0_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_B0B0/(BR_BpBp+BR_B0B0)) * BR_K0_nunubar)
# define N_K0star_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_B0B0/(BR_BpBp+BR_B0B0)) * BR_K0star_nunubar)
# define N_Xsd_nunubar_LS1 (2.0 * N_BB_LS1 * (BR_B0B0/(BR_BpBp+BR_B0B0)) * BR_Xsd_nonresonant_nunubar)

// SIGNAL MC sample number befor skimming
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

// scale factor for SIGNAL MC sample until LS1
# define Scale_Kplus_train (N_Kplus_nunubar_LS1/N_Kplus_train)
# define Scale_Kplusstar_train (N_Kplusstar_nunubar_LS1/N_Kplusstar_train)
# define Scale_Xsu_nonresonant_train (N_Xsu_nonresonant_nunubar_LS1/N_Xsu_nonresonant_train)
# define Scale_K0_train (N_K0_nunubar_LS1/N_K0_train)
# define Scale_K0star_train (N_K0star_nunubar_LS1/N_K0star_train)
# define Scale_Xsd_nonresonant_train (N_Xsd_nunubar_LS1/N_Xsd_nonresonant_train)
# define Scale_Kplus_test (N_Kplus_nunubar_LS1/N_Kplus_test)
# define Scale_Kplusstar_test (N_Kplusstar_nunubar_LS1/N_Kplusstar_test)
# define Scale_Xsu_nonresonant_test (N_Xsu_nonresonant_nunubar_LS1/N_Xsu_nonresonant_test)
# define Scale_K0_test (N_K0_nunubar_LS1/N_K0_test)
# define Scale_K0star_test (N_K0star_nunubar_LS1/N_K0star_test)
# define Scale_Xsd_nonresonant_test (N_Xsd_nunubar_LS1/N_Xsd_nonresonant_test)

// BKG MC sample number (2.8/ab for BB, 1.0/ab for qq)
# define N_CHG_test 139768443.0
# define N_MIX_test 104591690.0
# define N_UUBAR_test 107641351.0
# define N_DDBAR_test 25588403.0
# define N_SSBAR_test 20668794.0
# define N_CHARM_test 116108850.0
# define N_CHG_train 151915459.0
# define N_MIX_train 120136353.0
# define N_UUBAR_train 158509639.0
# define N_DDBAR_train 38644413.0
# define N_SSBAR_train 31001866.0
# define N_CHARM_train 174901296.0

// new scale factor for BKG MC sample with additional 1/ab (364.436 - 2.763 = 361.673/fb), until LS1
# define Scale_CHG_train ((N_BB_LS1* (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab * (N_CHG_train / (N_CHG_train + N_CHG_test)) ))
# define Scale_MIX_train ((N_BB_LS1* (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab * (N_MIX_train / (N_MIX_train + N_MIX_test)) ))
# define Scale_UUBAR_train (0.361673/((N_UUBAR_train/(N_UUBAR_train + N_UUBAR_test))*1.0))
# define Scale_DDBAR_train (0.361673/((N_DDBAR_train/(N_DDBAR_train + N_DDBAR_test))*1.0))
# define Scale_SSBAR_train (0.361673/((N_SSBAR_train/(N_SSBAR_train + N_SSBAR_test))*1.0))
# define Scale_CHARM_train (0.361673/((N_CHARM_train/(N_CHARM_train + N_CHARM_test))*1.0))
# define Scale_CHG_test ((N_BB_LS1* (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab * (N_CHG_test / (N_CHG_train + N_CHG_test)) ))
# define Scale_MIX_test ((N_BB_LS1* (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab * (N_MIX_test / (N_MIX_train + N_MIX_test)) ))
# define Scale_UUBAR_test (0.361673/((N_UUBAR_test/(N_UUBAR_train + N_UUBAR_test))*1.0))
# define Scale_DDBAR_test (0.361673/((N_DDBAR_test/(N_DDBAR_train + N_DDBAR_test))*1.0))
# define Scale_SSBAR_test (0.361673/((N_SSBAR_test/(N_SSBAR_train + N_SSBAR_test))*1.0))
# define Scale_CHARM_test (0.361673/((N_CHARM_test/(N_CHARM_train + N_CHARM_test))*1.0))

# define KS0_rel_uncertainty 0.6 // %/cm
# define track_rel_uncertainty 0.69 // %
// # define pi0_correction 0.932
# define pi0_rel_uncertainty ((0.0369 / 0.932) * 100.0) // %
# define Kaon_PID_max_uncertainty 0.1 // not percentage. relative uncertainty

# define FEI_cal_Bc_num 12
# define FEI_cal_B0_num 11
double FEI_cal_Bc[FEI_cal_Bc_num] = { 1.04, 0.79, 0.69, 0.56, 0.97, 0.95, 0.74, 0.57, 0.91, 0.51, 0.34, 0.59 };
double FEI_cal_Bc_uncertainty[FEI_cal_Bc_num] = { 0.03, 0.03, 0.05, 0.11, 0.03, 0.03, 0.02, 0.06, 0.1, 0.13, 0.07, 0.02 }; // not relative uncertainty. absolute uncertainty
double FEI_cal_Bc_modeID[FEI_cal_Bc_num] = { 0.0, 1.0, 3.0, 4.0, 15.0, 16.0, 18.0, 19.0, 23.0, 24.0, 30.0, -1.0 };
double FEI_cal_B0[FEI_cal_B0_num] = { 1.16, 0.95, 0.84, 0.78, 0.99, 1.01, 0.67, 0.65, 0.69, 0.58, 0.81 };
double FEI_cal_B0_uncertainty[FEI_cal_B0_num] = { 0.04, 0.03, 0.02, 0.02, 0.03, 0.03, 0.02, 0.02, 0.02, 0.16, 0.13 }; // not relative uncertainty. absolute uncertainty
double FEI_cal_B0_modeID[FEI_cal_B0_num] = { 0.0, 1.0, 3.0, 4.0, 5.0, 15.0, 16.0, 18.0, 19.0, 26.0, -1.0 };

# define CAL 1.1728
# define CAL_qq 1.0
# define Stream 0.25

# define N_PID_syst 73
# define N_fakeE_syst 37
# define N_fakeMU_syst 49
double PID_correction[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
double PID_correction_stat_uncer[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
double PID_correction_sys_uncer[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
double PID_correction_uncer[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss

double PID_fakeE_correction[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
double PID_fakeE_uncer[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+

double PID_fakeMU_correction[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+
double PID_fakeMU_uncer[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

# define N_pi0_syst 8
const double pi0_correction[N_pi0_syst] = {
    0.917, 0.965, 0.988, 1.013, 1.042, 1.044, 1.011, 1.0
};
const double pi0_stat_uncer[N_pi0_syst] = {
    0.004, 0.004, 0.004, 0.005, 0.004, 0.005, 0.005, 0.0
};
const double pi0_sys_uncer1[N_pi0_syst] = {
    0.049, 0.036, 0.079, 0.058, 0.045, 0.041, 0.040, 0.0
};
const double pi0_sys_uncer2[N_pi0_syst] = {
    0.0, 0.0, 0.0, 0.0, 0.039, 0.051, 0.030, 0.0
};

class Corrector_Knn {
private:

    // K+nn
    int STEP_Knn;
    double mininvM_Knn;
    double maxinvM_Knn;
    TH1D* weights_Knn;
    const double DECAY_DEC_BR_Knn;
    const double new_BR_K0pp;
    double new_BR_Knn;
    const double Nraw_initial_Knn;
    double Nscale_initial_Knn;

    // K*+nn
    int STEP_Kstarnn;
    double mininvM_Kstarnn;
    double maxinvM_Kstarnn;
    TH1D* weights_Kstarnn;
    const double DECAY_DEC_BR_Kstarnn;
    const double new_BR_K0starpp;
    double new_BR_Kstarnn;
    const double Nraw_initial_Kstarnn;
    double Nscale_initial_Kstarnn;

    // K0nn
    int STEP_K0nn;
    double mininvM_K0nn;
    double maxinvM_K0nn;
    TH1D* weights_K0nn;
    const double DECAY_DEC_BR_K0nn;
    const double new_BR_Kpp;
    double new_BR_K0nn;
    const double Nraw_initial_K0nn;
    double Nscale_initial_K0nn;

    // K0*nn
    int STEP_K0starnn;
    double mininvM_K0starnn;
    double maxinvM_K0starnn;
    TH1D* weights_K0starnn;
    const double DECAY_DEC_BR_K0starnn;
    const double new_BR_Kstarpp;
    double new_BR_K0starnn;
    const double Nraw_initial_K0starnn;
    double Nscale_initial_K0starnn;

    const double N_EPSILON;
    const double CUTOFF;

    const double tau_Bp;
    const double tau_B0;

public:
    Corrector_Knn();
    double GetCorrectionFactor(double invM_Knn, double invM_Kstarnn, double invM_K0nn, double invM_K0starnn, double N_Knn, double N_Kstarnn, double N_K0nn, double N_K0starnn);
    double GetCorrectionFactorAtGeneric(double invM_Knn, double invM_Kstarnn, double invM_K0nn, double invM_K0starnn, double N_Knn, double N_Kstarnn, double N_K0nn, double N_K0starnn);
};

Corrector_Knn corrector_Knn;

Corrector_Knn::Corrector_Knn() :
    DECAY_DEC_BR_Knn(0.0000057),
    new_BR_K0pp(0.00000266),
    Nraw_initial_Knn(1000000.0),
    DECAY_DEC_BR_Kstarnn(0.0000057),
    new_BR_K0starpp(0.00000124),
    Nraw_initial_Kstarnn(1000000.0),
    DECAY_DEC_BR_K0nn(0.000002),
    new_BR_Kpp(0.0000059),
    Nraw_initial_K0nn(1000000.0),
    DECAY_DEC_BR_K0starnn(0.0000056),
    new_BR_Kstarpp(0.0000036),
    Nraw_initial_K0starnn(1000000.0),
    N_EPSILON(0.01),
    CUTOFF(50.0),
    tau_Bp(1.6384), // ps
    tau_B0(1.5195) // ps
{
    FILE* fp;

    // read Knn weights
    fp = fopen("/home/jwpark/storage/BKG_gbasf2/systematic/Knn_weight/Knn_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &STEP_Knn, &mininvM_Knn, &maxinvM_Knn);
    weights_Knn = new TH1D("Knn_weights", ";;", STEP_Knn, mininvM_Knn, maxinvM_Knn);
    for (int i = 0; i < STEP_Knn; i++) {
        double temp;
        fscanf(fp, "%lf\n", &temp);
        if (temp < CUTOFF) weights_Knn->SetBinContent(i + 1, temp);
        else weights_Knn->SetBinContent(i + 1, CUTOFF);
    }
    fclose(fp);

    // read Kstarnn weights
    fp = fopen("/home/jwpark/storage/BKG_gbasf2/systematic/Knn_weight/Kstarnn_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &STEP_Kstarnn, &mininvM_Kstarnn, &maxinvM_Kstarnn);
    weights_Kstarnn = new TH1D("Kstarnn_weights", ";;", STEP_Kstarnn, mininvM_Kstarnn, maxinvM_Kstarnn);
    for (int i = 0; i < STEP_Kstarnn; i++) {
        double temp;
        fscanf(fp, "%lf\n", &temp);
        if (temp < CUTOFF) weights_Kstarnn->SetBinContent(i + 1, temp);
        else weights_Kstarnn->SetBinContent(i + 1, CUTOFF);
    }
    fclose(fp);

    // read K0nn weights
    fp = fopen("/home/jwpark/storage/BKG_gbasf2/systematic/Knn_weight/K0nn_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &STEP_K0nn, &mininvM_K0nn, &maxinvM_K0nn);
    weights_K0nn = new TH1D("K0nn_weights", ";;", STEP_K0nn, mininvM_K0nn, maxinvM_K0nn);
    for (int i = 0; i < STEP_K0nn; i++) {
        double temp;
        fscanf(fp, "%lf\n", &temp);
        if (temp < CUTOFF) weights_K0nn->SetBinContent(i + 1, temp);
        else weights_K0nn->SetBinContent(i + 1, CUTOFF);
    }
    fclose(fp);

    // read K0starnn weights
    fp = fopen("/home/jwpark/storage/BKG_gbasf2/systematic/Knn_weight/K0starnn_weight.txt", "r");
    fscanf(fp, "%d %lf %lf\n", &STEP_K0starnn, &mininvM_K0starnn, &maxinvM_K0starnn);
    weights_K0starnn = new TH1D("K0starnn_weights", ";;", STEP_K0starnn, mininvM_K0starnn, maxinvM_K0starnn);
    for (int i = 0; i < STEP_K0starnn; i++) {
        double temp;
        fscanf(fp, "%lf\n", &temp);
        if (temp < CUTOFF) weights_K0starnn->SetBinContent(i + 1, temp);
        else weights_K0starnn->SetBinContent(i + 1, CUTOFF);
    }
    fclose(fp);

    // calculate the expected number of event
    // use proper isospin
    new_BR_Knn = new_BR_K0pp * (tau_Bp / tau_B0);
    new_BR_Kstarnn = new_BR_K0starpp * (tau_Bp / tau_B0);
    new_BR_K0nn = new_BR_Kpp * (tau_B0 / tau_Bp);
    new_BR_K0starnn = new_BR_Kstarpp * (tau_B0 / tau_Bp);

    Nscale_initial_Knn = (2.0 * N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0)) * new_BR_Knn);
    Nscale_initial_Kstarnn = (2.0 * N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0)) * new_BR_Kstarnn);
    Nscale_initial_K0nn = (2.0 * N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0)) * new_BR_K0nn);
    Nscale_initial_K0starnn = (2.0 * N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0)) * new_BR_K0starnn);
}

double Corrector_Knn::GetCorrectionFactor(double invM_Knn, double invM_Kstarnn, double invM_K0nn, double invM_K0starnn, double N_Knn, double N_Kstarnn, double N_K0nn, double N_K0starnn) {

    double Correction_Knn = 1;
    double Correction_Kstarnn = 1;
    double Correction_K0nn = 1;
    double Correction_K0starnn = 1;

    if (N_Knn < N_EPSILON) Correction_Knn = 1;
    else {
        int Bin = weights_Knn->FindBin(invM_Knn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_Knn) Bin = STEP_Knn;
        Correction_Knn = std::pow((Nscale_initial_Knn / Nraw_initial_Knn) * weights_Knn->GetBinContent(Bin), N_Knn); // BR correction * invM correction
    }

    if (N_Kstarnn < N_EPSILON) Correction_Kstarnn = 1;
    else {
        int Bin = weights_Kstarnn->FindBin(invM_Kstarnn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_Kstarnn) Bin = STEP_Kstarnn;
        Correction_Kstarnn = std::pow((Nscale_initial_Kstarnn / Nraw_initial_Kstarnn) * weights_Kstarnn->GetBinContent(Bin), N_Kstarnn); // BR correction * invM correction
    }

    if (N_K0nn < N_EPSILON) Correction_K0nn = 1;
    else {
        int Bin = weights_K0nn->FindBin(invM_K0nn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_K0nn) Bin = STEP_K0nn;
        Correction_K0nn = std::pow((Nscale_initial_K0nn / Nraw_initial_K0nn) * weights_K0nn->GetBinContent(Bin), N_K0nn); // BR correction * invM correction
    }

    if (N_K0starnn < N_EPSILON) Correction_K0starnn = 1;
    else {
        int Bin = weights_K0starnn->FindBin(invM_K0starnn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_K0starnn) Bin = STEP_K0starnn;
        Correction_K0starnn = std::pow((Nscale_initial_K0starnn / Nraw_initial_K0starnn) * weights_K0starnn->GetBinContent(Bin), N_K0starnn); // BR correction * invM correction
    }

    return Correction_Knn * Correction_Kstarnn * Correction_K0nn * Correction_K0starnn;
}

double Corrector_Knn::GetCorrectionFactorAtGeneric(double invM_Knn, double invM_Kstarnn, double invM_K0nn, double invM_K0starnn, double N_Knn, double N_Kstarnn, double N_K0nn, double N_K0starnn) {

    double Correction_Knn = 1;
    double Correction_Kstarnn = 1;
    double Correction_K0nn = 1;
    double Correction_K0starnn = 1;

    if (N_Knn < N_EPSILON) Correction_Knn = 1;
    else {
        int Bin = weights_Knn->FindBin(invM_Knn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_Knn) Bin = STEP_Knn;
        Correction_Knn = std::pow((new_BR_Knn / DECAY_DEC_BR_Knn) * weights_Knn->GetBinContent(Bin), N_Knn); // BR correction * invM correction
    }

    if (N_Kstarnn < N_EPSILON) Correction_Kstarnn = 1;
    else {
        int Bin = weights_Kstarnn->FindBin(invM_Kstarnn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_Kstarnn) Bin = STEP_Kstarnn;
        Correction_Kstarnn = std::pow((new_BR_Kstarnn / DECAY_DEC_BR_Kstarnn) * weights_Kstarnn->GetBinContent(Bin), N_Kstarnn); // BR correction * invM correction
    }

    if (N_K0nn < N_EPSILON) Correction_K0nn = 1;
    else {
        int Bin = weights_K0nn->FindBin(invM_K0nn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_K0nn) Bin = STEP_K0nn;
        Correction_K0nn = std::pow((new_BR_K0nn / DECAY_DEC_BR_K0nn) * weights_K0nn->GetBinContent(Bin), N_K0nn); // BR correction * invM correction
    }

    if (N_K0starnn < N_EPSILON) Correction_K0starnn = 1;
    else {
        int Bin = weights_K0starnn->FindBin(invM_K0starnn);
        if (Bin < 1) Bin = 1;
        else if (Bin > STEP_K0starnn) Bin = STEP_K0starnn;
        Correction_K0starnn = std::pow((new_BR_K0starnn / DECAY_DEC_BR_K0starnn) * weights_K0starnn->GetBinContent(Bin), N_K0starnn); // BR correction * invM correction
    }

    return Correction_Knn * Correction_Kstarnn * Correction_K0nn * Correction_K0starnn;
}

// Corrector_Multiplicity is fixed for MultiplicityCalculator! Do not copy-paste it!

class Corrector_Multiplicity {
private:

    int NgammaMAX_;
    TH1D* weights_Ngamma_zero;
    TH1D* weights_Ngamma_one;
    TH1D* weights_Ngamma_two;
    const double CUTOFF;

public:
    Corrector_Multiplicity();
    double GetCorrectionFactor(double Ngamma, int pi0_multiplicity);
};

Corrector_Multiplicity *corrector_Multiplicity;

Corrector_Multiplicity::Corrector_Multiplicity() :
    CUTOFF(50.0)
{
    FILE* fp;

    // read Knn weights
    fp = fopen("./multiplicity_weight_zero.txt", "r");
    fscanf(fp, "%d\n", &NgammaMAX_);
    weights_Ngamma_zero = new TH1D("weights_Ngamma_zero", ";;", NgammaMAX_ + 1, -0.5, NgammaMAX_ + 0.5);
    for (int i = 0; i < NgammaMAX_ + 1; i++) {
        double temp1;
        double temp2;
        double temp3;
        fscanf(fp, "%lf %lf %lf\n", &temp1, &temp2, &temp3);
        if (temp3 < CUTOFF) weights_Ngamma_zero->SetBinContent(i + 1, temp3);
        else weights_Ngamma_zero->SetBinContent(i + 1, CUTOFF);
    }
    fclose(fp);

    fp = fopen("./multiplicity_weight_one.txt", "r");
    fscanf(fp, "%d\n", &NgammaMAX_);
    weights_Ngamma_one = new TH1D("weights_Ngamma_one", ";;", NgammaMAX_ + 1, -0.5, NgammaMAX_ + 0.5);
    for (int i = 0; i < NgammaMAX_ + 1; i++) {
        double temp1;
        double temp2;
        double temp3;
        fscanf(fp, "%lf %lf %lf\n", &temp1, &temp2, &temp3);
        if (temp3 < CUTOFF) weights_Ngamma_one->SetBinContent(i + 1, temp3);
        else weights_Ngamma_one->SetBinContent(i + 1, CUTOFF);
    }
    fclose(fp);

    fp = fopen("./multiplicity_weight_two.txt", "r");
    fscanf(fp, "%d\n", &NgammaMAX_);
    weights_Ngamma_two = new TH1D("weights_Ngamma_two", ";;", NgammaMAX_ + 1, -0.5, NgammaMAX_ + 0.5);
    for (int i = 0; i < NgammaMAX_ + 1; i++) {
        double temp1;
        double temp2;
        double temp3;
        fscanf(fp, "%lf %lf %lf\n", &temp1, &temp2, &temp3);
        if (temp3 < CUTOFF) weights_Ngamma_two->SetBinContent(i + 1, temp3);
        else weights_Ngamma_two->SetBinContent(i + 1, CUTOFF);
    }
    fclose(fp);

}

double Corrector_Multiplicity::GetCorrectionFactor(double Ngamma, int pi0_multiplicity) {
    if (pi0_multiplicity == 0) {
        int Bin = weights_Ngamma_zero->FindBin(Ngamma);
        if (Bin < 1) {
            printf("[ERROR] Ngamma is smaller than 0!\n");
            exit(1);
        }
        else if (Bin > NgammaMAX_ + 1) return 1.0;

        return weights_Ngamma_zero->GetBinContent(Bin);
    }
    else if (pi0_multiplicity == 1) {
        int Bin = weights_Ngamma_one->FindBin(Ngamma);
        if (Bin < 1) {
            printf("[ERROR] Ngamma is smaller than 0!\n");
            exit(1);
        }
        else if (Bin > NgammaMAX_ + 1) return 1.0;

        return weights_Ngamma_one->GetBinContent(Bin);
    }
    else if (pi0_multiplicity == 2) {
        int Bin = weights_Ngamma_two->FindBin(Ngamma);
        if (Bin < 1) {
            printf("[ERROR] Ngamma is smaller than 0!\n");
            exit(1);
        }
        else if (Bin > NgammaMAX_ + 1) return 1.0;

        return weights_Ngamma_two->GetBinContent(Bin);
    }
    else {
        exit(1);
        return 1;
    }

}

/* ====================================== */

void ReadPIDFile() {
    const char* KID_true_file = "/home/jwpark/storage/BKG_gbasf2/systematic/MC15ri_PID/KaonEff.csv";
    const char* KID_mis_file = "/home/jwpark/storage/BKG_gbasf2/systematic/MC15ri_PID/Kaonmis.csv";
    const char* PID_true_file = "/home/jwpark/storage/BKG_gbasf2/systematic/MC15ri_PID/PionEff.csv";
    const char* PID_mis_file = "/home/jwpark/storage/BKG_gbasf2/systematic/MC15ri_PID/Pionmis.csv";

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

        if ((std::abs(PID_correction[0][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer[0][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer[0][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction[0][i]) > 10000.0)) {
            PID_correction[0][i] = 1.0;
            PID_correction_stat_uncer[0][i] = 0.0;
            PID_correction_sys_uncer[0][i] = 0.0;
        }
        if ((std::abs(PID_correction[1][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer[1][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer[1][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction[1][i]) > 10000.0)) {
            PID_correction[1][i] = 1.0;
            PID_correction_stat_uncer[1][i] = 0.0;
            PID_correction_sys_uncer[1][i] = 0.0;
        }
        if ((std::abs(PID_correction[2][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer[2][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer[2][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction[2][i]) > 10000.0)) {
            PID_correction[2][i] = 1.0;
            PID_correction_stat_uncer[2][i] = 0.0;
            PID_correction_sys_uncer[2][i] = 0.0;
        }
        if ((std::abs(PID_correction[3][i]) < MyEPSILON && std::abs(PID_correction_stat_uncer[3][i] - 1.0) < MyEPSILON && std::abs(PID_correction_sys_uncer[3][i] - 1.0) < MyEPSILON) || (std::abs(PID_correction[3][i]) > 10000.0)) {
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

void ReadFakePIDFile() {
    // initialization
    for (int i = 0; i < N_fakeE_syst; i++) {
        PID_fakeE_correction[0][i] = 1.0; //  K-, K+, pi-, pi+
        PID_fakeE_correction[1][i] = 1.0;
        PID_fakeE_correction[2][i] = 1.0;
        PID_fakeE_correction[3][i] = 1.0;

        PID_fakeE_uncer[0][i] = 0.0;
        PID_fakeE_uncer[1][i] = 0.0;
        PID_fakeE_uncer[2][i] = 0.0;
        PID_fakeE_uncer[3][i] = 0.0;
    }

    for (int i = 0; i < N_fakeMU_syst; i++) {
        PID_fakeMU_correction[0][i] = 1.0;
        PID_fakeMU_correction[1][i] = 1.0;
        PID_fakeMU_correction[2][i] = 1.0;
        PID_fakeMU_correction[3][i] = 1.0;

        PID_fakeMU_uncer[0][i] = 0.0;
        PID_fakeMU_uncer[1][i] = 0.0;
        PID_fakeMU_uncer[2][i] = 0.0;
        PID_fakeMU_uncer[3][i] = 0.0;
    }

    const char* K_fromE_file = "/home/jwpark/storage/BKG_gbasf2/systematic/MC15ri_PID/kaonID_efficiency_electron.csv";
    const char* K_fromMU_file = "/home/jwpark/storage/BKG_gbasf2/systematic/MC15ri_PID/kaonID_efficiency_muon.csv";
    const char* pi_fromE_file = "/home/jwpark/storage/BKG_gbasf2/systematic/MC15ri_PID/pionID_efficiency_electron.csv";
    const char* pi_fromMU_file = "/home/jwpark/storage/BKG_gbasf2/systematic/MC15ri_PID/pionID_efficiency_muon.csv";

    FILE* fp_K_fromE = fopen(K_fromE_file, "r");
    FILE* fp_K_fromMU = fopen(K_fromMU_file, "r");
    FILE* fp_pi_fromE = fopen(pi_fromE_file, "r");
    FILE* fp_pi_fromMU = fopen(pi_fromMU_file, "r");

    fscanf(fp_K_fromE, "variable,charge,p_min,p_max,theta_min,theta_max,iso_score_min,iso_score_max,working_point,threshold,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn\n");
    fscanf(fp_K_fromMU, "variable,charge,p_min,p_max,theta_min,theta_max,iso_score_min,iso_score_max,working_point,threshold,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn\n");
    fscanf(fp_pi_fromE, "variable,charge,p_min,p_max,theta_min,theta_max,iso_score_min,iso_score_max,working_point,threshold,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn\n");
    fscanf(fp_pi_fromMU, "variable,charge,p_min,p_max,theta_min,theta_max,iso_score_min,iso_score_max,working_point,threshold,data_efficiency,data_uncertainty_stat_up,data_uncertainty_stat_dn,data_uncertainty_sys_up,data_uncertainty_sys_dn,MC_efficiency,MC_uncertainty_stat_up,MC_uncertainty_stat_dn,MC_uncertainty_sys_up,MC_uncertainty_sys_dn,data_MC_ratio,data_MC_uncertainty_stat_up,data_MC_uncertainty_stat_dn,data_MC_uncertainty_sys_up,data_MC_uncertainty_sys_dn\n");

    char temp_charge;
    double temp_p_min;
    double temp_p_max;
    double temp_theta_min;
    double temp_theta_max;
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

    // Kaon from fake electron
    while (fscanf(fp_K_fromE, "kaonID,%c,%lf,%lf,%lf,%lf,0.0,1.0,FixedThresh06,0.6,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &temp_charge, &temp_p_min, &temp_p_max, &temp_theta_min, &temp_theta_max, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_data_MC_ratio, &temp_data_MC_uncertainty_stat_up, &temp_data_MC_uncertainty_stat_dn, &temp_data_MC_uncertainty_sys_up, &temp_data_MC_uncertainty_sys_dn) != EOF) {
        int p_bin = -1;
        int theta_bin = -1;

        if (std::abs(temp_p_min - 0.4) < MyEPSILON && std::abs(temp_p_max - 0.5) < MyEPSILON) p_bin = 0;
        else if (std::abs(temp_p_min - 0.5) < MyEPSILON && std::abs(temp_p_max - 1.0) < MyEPSILON) p_bin = 1;
        else if (std::abs(temp_p_min - 1.0) < MyEPSILON && std::abs(temp_p_max - 1.5) < MyEPSILON) p_bin = 2;
        else if (std::abs(temp_p_min - 1.5) < MyEPSILON && std::abs(temp_p_max - 2.0) < MyEPSILON) p_bin = 3;
        else if (std::abs(temp_p_min - 2.0) < MyEPSILON && std::abs(temp_p_max - 2.5) < MyEPSILON) p_bin = 4;
        else if (std::abs(temp_p_min - 2.5) < MyEPSILON && std::abs(temp_p_max - 3.0) < MyEPSILON) p_bin = 5;
        else {
            printf("[ERROR] unknown p bin!\n");
            exit(1);
        }

        if (std::abs(temp_theta_min - 0.22) < MyEPSILON && std::abs(temp_theta_max - 0.56) < MyEPSILON) theta_bin = 0;
        else if (std::abs(temp_theta_min - 0.56) < MyEPSILON && std::abs(temp_theta_max - 1.13) < MyEPSILON) theta_bin = 1;
        else if (std::abs(temp_theta_min - 1.13) < MyEPSILON && std::abs(temp_theta_max - 1.57) < MyEPSILON) theta_bin = 2;
        else if (std::abs(temp_theta_min - 1.57) < MyEPSILON && std::abs(temp_theta_max - 1.88) < MyEPSILON) theta_bin = 3;
        else if (std::abs(temp_theta_min - 1.88) < MyEPSILON && std::abs(temp_theta_max - 2.23) < MyEPSILON) theta_bin = 4;
        else if (std::abs(temp_theta_min - 2.23) < MyEPSILON && std::abs(temp_theta_max - 2.71) < MyEPSILON) theta_bin = 5;
        else {
            printf("[ERROR] unknown theta bin!\n");
            exit(1);
        }

        int bin = theta_bin + 6 * p_bin;

        if (temp_charge == '+') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeE_correction[1][bin] = temp_data_MC_ratio;
            PID_fakeE_uncer[1][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeE_correction[1][bin]) < MyEPSILON) {
                PID_fakeE_correction[1][bin] = 1.0;
                PID_fakeE_uncer[1][bin] = 0.0;
            }

        }
        else if (temp_charge == '-') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeE_correction[0][bin] = temp_data_MC_ratio;
            PID_fakeE_uncer[0][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeE_correction[0][bin]) < MyEPSILON) {
                PID_fakeE_correction[0][bin] = 1.0;
                PID_fakeE_uncer[0][bin] = 0.0;
            }

        }
        else {
            printf("[ERROR] unknown charge!\n");
            exit(1);
        }
    }
    fclose(fp_K_fromE);

    // Pion from fake electron
    while (fscanf(fp_pi_fromE, "pionID,%c,%lf,%lf,%lf,%lf,0.0,1.0,FixedThresh06,0.6,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &temp_charge, &temp_p_min, &temp_p_max, &temp_theta_min, &temp_theta_max, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_data_MC_ratio, &temp_data_MC_uncertainty_stat_up, &temp_data_MC_uncertainty_stat_dn, &temp_data_MC_uncertainty_sys_up, &temp_data_MC_uncertainty_sys_dn) != EOF) {
        int p_bin = -1;
        int theta_bin = -1;

        if (std::abs(temp_p_min - 0.4) < MyEPSILON && std::abs(temp_p_max - 0.5) < MyEPSILON) p_bin = 0;
        else if (std::abs(temp_p_min - 0.5) < MyEPSILON && std::abs(temp_p_max - 1.0) < MyEPSILON) p_bin = 1;
        else if (std::abs(temp_p_min - 1.0) < MyEPSILON && std::abs(temp_p_max - 1.5) < MyEPSILON) p_bin = 2;
        else if (std::abs(temp_p_min - 1.5) < MyEPSILON && std::abs(temp_p_max - 2.0) < MyEPSILON) p_bin = 3;
        else if (std::abs(temp_p_min - 2.0) < MyEPSILON && std::abs(temp_p_max - 2.5) < MyEPSILON) p_bin = 4;
        else if (std::abs(temp_p_min - 2.5) < MyEPSILON && std::abs(temp_p_max - 3.0) < MyEPSILON) p_bin = 5;
        else {
            printf("[ERROR] unknown p bin!\n");
            exit(1);
        }

        if (std::abs(temp_theta_min - 0.22) < MyEPSILON && std::abs(temp_theta_max - 0.56) < MyEPSILON) theta_bin = 0;
        else if (std::abs(temp_theta_min - 0.56) < MyEPSILON && std::abs(temp_theta_max - 1.13) < MyEPSILON) theta_bin = 1;
        else if (std::abs(temp_theta_min - 1.13) < MyEPSILON && std::abs(temp_theta_max - 1.57) < MyEPSILON) theta_bin = 2;
        else if (std::abs(temp_theta_min - 1.57) < MyEPSILON && std::abs(temp_theta_max - 1.88) < MyEPSILON) theta_bin = 3;
        else if (std::abs(temp_theta_min - 1.88) < MyEPSILON && std::abs(temp_theta_max - 2.23) < MyEPSILON) theta_bin = 4;
        else if (std::abs(temp_theta_min - 2.23) < MyEPSILON && std::abs(temp_theta_max - 2.71) < MyEPSILON) theta_bin = 5;
        else {
            printf("[ERROR] unknown theta bin!\n");
            exit(1);
        }

        int bin = theta_bin + 6 * p_bin;

        if (temp_charge == '+') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeE_correction[3][bin] = temp_data_MC_ratio;
            PID_fakeE_uncer[3][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeE_correction[3][bin]) < MyEPSILON) {
                PID_fakeE_correction[3][bin] = 1.0;
                PID_fakeE_uncer[3][bin] = 0.0;
            }

        }
        else if (temp_charge == '-') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeE_correction[2][bin] = temp_data_MC_ratio;
            PID_fakeE_uncer[2][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeE_correction[2][bin]) < MyEPSILON) {
                PID_fakeE_correction[2][bin] = 1.0;
                PID_fakeE_uncer[2][bin] = 0.0;
            }

        }
        else {
            printf("[ERROR] unknown charge!\n");
            exit(1);
        }
    }
    fclose(fp_pi_fromE);

    // Kaon from fake muon
    while (fscanf(fp_K_fromMU, "kaonID,%c,%lf,%lf,%lf,%lf,0.0,1.0,FixedThresh06,0.6,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &temp_charge, &temp_p_min, &temp_p_max, &temp_theta_min, &temp_theta_max, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_data_MC_ratio, &temp_data_MC_uncertainty_stat_up, &temp_data_MC_uncertainty_stat_dn, &temp_data_MC_uncertainty_sys_up, &temp_data_MC_uncertainty_sys_dn) != EOF) {
        int p_bin = -1;
        int theta_bin = -1;

        if (std::abs(temp_p_min - 0.4) < MyEPSILON && std::abs(temp_p_max - 0.5) < MyEPSILON) p_bin = 0;
        else if (std::abs(temp_p_min - 0.5) < MyEPSILON && std::abs(temp_p_max - 0.7) < MyEPSILON) p_bin = 1;
        else if (std::abs(temp_p_min - 0.7) < MyEPSILON && std::abs(temp_p_max - 1.0) < MyEPSILON) p_bin = 2;
        else if (std::abs(temp_p_min - 1.0) < MyEPSILON && std::abs(temp_p_max - 1.5) < MyEPSILON) p_bin = 3;
        else if (std::abs(temp_p_min - 1.5) < MyEPSILON && std::abs(temp_p_max - 2.0) < MyEPSILON) p_bin = 4;
        else if (std::abs(temp_p_min - 2.0) < MyEPSILON && std::abs(temp_p_max - 2.5) < MyEPSILON) p_bin = 5;
        else {
            printf("[ERROR] unknown p bin!\n");
            exit(1);
        }

        if (std::abs(temp_theta_min - 0.4) < MyEPSILON && std::abs(temp_theta_max - 0.64) < MyEPSILON) theta_bin = 0;
        else if (std::abs(temp_theta_min - 0.64) < MyEPSILON && std::abs(temp_theta_max - 0.82) < MyEPSILON) theta_bin = 1;
        else if (std::abs(temp_theta_min - 0.82) < MyEPSILON && std::abs(temp_theta_max - 1.16) < MyEPSILON) theta_bin = 2;
        else if (std::abs(temp_theta_min - 1.16) < MyEPSILON && std::abs(temp_theta_max - 1.46) < MyEPSILON) theta_bin = 3;
        else if (std::abs(temp_theta_min - 1.46) < MyEPSILON && std::abs(temp_theta_max - 1.78) < MyEPSILON) theta_bin = 4;
        else if (std::abs(temp_theta_min - 1.78) < MyEPSILON && std::abs(temp_theta_max - 2.13) < MyEPSILON) theta_bin = 5;
        else if (std::abs(temp_theta_min - 2.13) < MyEPSILON && std::abs(temp_theta_max - 2.22) < MyEPSILON) theta_bin = 6;
        else if (std::abs(temp_theta_min - 2.22) < MyEPSILON && std::abs(temp_theta_max - 2.6) < MyEPSILON) theta_bin = 7;
        else {
            printf("[ERROR] unknown theta bin!\n");
            exit(1);
        }

        int bin = theta_bin + 8 * p_bin;

        if (temp_charge == '+') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeMU_correction[1][bin] = temp_data_MC_ratio;
            PID_fakeMU_uncer[1][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeMU_correction[1][bin]) < MyEPSILON) {
                PID_fakeMU_correction[1][bin] = 1.0;
                PID_fakeMU_uncer[1][bin] = 0.0;
            }

        }
        else if (temp_charge == '-') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeMU_correction[0][bin] = temp_data_MC_ratio;
            PID_fakeMU_uncer[0][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeMU_correction[0][bin]) < MyEPSILON) {
                PID_fakeMU_correction[0][bin] = 1.0;
                PID_fakeMU_uncer[0][bin] = 0.0;
            }

        }
        else {
            printf("[ERROR] unknown charge!\n");
            exit(1);
        }
    }
    fclose(fp_K_fromMU);

    // Pion from fake muon
    while (fscanf(fp_pi_fromMU, "pionID,%c,%lf,%lf,%lf,%lf,0.0,1.0,FixedThresh06,0.6,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &temp_charge, &temp_p_min, &temp_p_max, &temp_theta_min, &temp_theta_max, &temp_data_efficiency, &temp_data_uncertainty_stat_up, &temp_data_uncertainty_stat_dn, &temp_data_uncertainty_sys_up, &temp_data_uncertainty_sys_dn, &temp_MC_efficiency, &temp_MC_uncertainty_stat_up, &temp_MC_uncertainty_stat_dn, &temp_MC_uncertainty_sys_up, &temp_MC_uncertainty_sys_dn, &temp_data_MC_ratio, &temp_data_MC_uncertainty_stat_up, &temp_data_MC_uncertainty_stat_dn, &temp_data_MC_uncertainty_sys_up, &temp_data_MC_uncertainty_sys_dn) != EOF) {
        int p_bin = -1;
        int theta_bin = -1;

        if (std::abs(temp_p_min - 0.4) < MyEPSILON && std::abs(temp_p_max - 0.5) < MyEPSILON) p_bin = 0;
        else if (std::abs(temp_p_min - 0.5) < MyEPSILON && std::abs(temp_p_max - 0.7) < MyEPSILON) p_bin = 1;
        else if (std::abs(temp_p_min - 0.7) < MyEPSILON && std::abs(temp_p_max - 1.0) < MyEPSILON) p_bin = 2;
        else if (std::abs(temp_p_min - 1.0) < MyEPSILON && std::abs(temp_p_max - 1.5) < MyEPSILON) p_bin = 3;
        else if (std::abs(temp_p_min - 1.5) < MyEPSILON && std::abs(temp_p_max - 2.0) < MyEPSILON) p_bin = 4;
        else if (std::abs(temp_p_min - 2.0) < MyEPSILON && std::abs(temp_p_max - 2.5) < MyEPSILON) p_bin = 5;
        else {
            printf("[ERROR] unknown p bin!\n");
            exit(1);
        }

        if (std::abs(temp_theta_min - 0.4) < MyEPSILON && std::abs(temp_theta_max - 0.64) < MyEPSILON) theta_bin = 0;
        else if (std::abs(temp_theta_min - 0.64) < MyEPSILON && std::abs(temp_theta_max - 0.82) < MyEPSILON) theta_bin = 1;
        else if (std::abs(temp_theta_min - 0.82) < MyEPSILON && std::abs(temp_theta_max - 1.16) < MyEPSILON) theta_bin = 2;
        else if (std::abs(temp_theta_min - 1.16) < MyEPSILON && std::abs(temp_theta_max - 1.46) < MyEPSILON) theta_bin = 3;
        else if (std::abs(temp_theta_min - 1.46) < MyEPSILON && std::abs(temp_theta_max - 1.78) < MyEPSILON) theta_bin = 4;
        else if (std::abs(temp_theta_min - 1.78) < MyEPSILON && std::abs(temp_theta_max - 2.13) < MyEPSILON) theta_bin = 5;
        else if (std::abs(temp_theta_min - 2.13) < MyEPSILON && std::abs(temp_theta_max - 2.22) < MyEPSILON) theta_bin = 6;
        else if (std::abs(temp_theta_min - 2.22) < MyEPSILON && std::abs(temp_theta_max - 2.6) < MyEPSILON) theta_bin = 7;
        else {
            printf("[ERROR] unknown theta bin!\n");
            exit(1);
        }

        int bin = theta_bin + 8 * p_bin;

        if (temp_charge == '+') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeMU_correction[3][bin] = temp_data_MC_ratio;
            PID_fakeMU_uncer[3][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeMU_correction[3][bin]) < MyEPSILON) {
                PID_fakeMU_correction[3][bin] = 1.0;
                PID_fakeMU_uncer[3][bin] = 0.0;
            }

        }
        else if (temp_charge == '-') {
            double temp_data_MC_uncertainty_up = std::sqrt(temp_data_MC_uncertainty_stat_up * temp_data_MC_uncertainty_stat_up + temp_data_MC_uncertainty_sys_up * temp_data_MC_uncertainty_sys_up);
            double temp_data_MC_uncertainty_dn = std::sqrt(temp_data_MC_uncertainty_stat_dn * temp_data_MC_uncertainty_stat_dn + temp_data_MC_uncertainty_sys_dn * temp_data_MC_uncertainty_sys_dn);

            PID_fakeMU_correction[2][bin] = temp_data_MC_ratio;
            PID_fakeMU_uncer[2][bin] = (temp_data_MC_uncertainty_up + temp_data_MC_uncertainty_dn) / 2.0;

            if (std::abs(PID_fakeMU_correction[2][bin]) < MyEPSILON) {
                PID_fakeMU_correction[2][bin] = 1.0;
                PID_fakeMU_uncer[2][bin] = 0.0;
            }

        }
        else {
            printf("[ERROR] unknown charge!\n");
            exit(1);
        }
    }
    fclose(fp_pi_fromMU);
}

double GetFEICalFactor(double UpsilonID, double BtagID) {
    // UpsilonID => charged: 0, mixed: 1

    if (UpsilonID > -0.5 && UpsilonID < 0.5) { // charged
        for (int i = 0; i < FEI_cal_Bc_num - 1; i++) {
            if (BtagID > FEI_cal_Bc_modeID[i] - 0.5 && BtagID < FEI_cal_Bc_modeID[i] + 0.5) return FEI_cal_Bc[i];
        }
        return FEI_cal_Bc[FEI_cal_Bc_num - 1];
    }
    else if (UpsilonID > 0.5 && UpsilonID < 1.5) { // mixed
        for (int i = 0; i < FEI_cal_B0_num - 1; i++) {
            if (BtagID > FEI_cal_B0_modeID[i] - 0.5 && BtagID < FEI_cal_B0_modeID[i] + 0.5) return FEI_cal_B0[i];
        }
        return FEI_cal_B0[FEI_cal_B0_num - 1];
    }

    printf("[GetFEICalFactor] error! unexpected decay ID\n");
    exit(1);
    return 0;

}

double GetFEICalFactorUncer(double UpsilonID, double BtagID) {
    // UpsilonID => charged: 0, mixed: 1

    if (UpsilonID > -0.5 && UpsilonID < 0.5) { // charged
        for (int i = 0; i < FEI_cal_Bc_num - 1; i++) {
            if (BtagID > FEI_cal_Bc_modeID[i] - 0.5 && BtagID < FEI_cal_Bc_modeID[i] + 0.5) return FEI_cal_Bc_uncertainty[i];
        }
        return FEI_cal_Bc_uncertainty[FEI_cal_Bc_num - 1];
    }
    else if (UpsilonID > 0.5 && UpsilonID < 1.5) { // mixed
        for (int i = 0; i < FEI_cal_B0_num - 1; i++) {
            if (BtagID > FEI_cal_B0_modeID[i] - 0.5 && BtagID < FEI_cal_B0_modeID[i] + 0.5) return FEI_cal_B0[i];
        }
        return FEI_cal_B0_uncertainty[FEI_cal_B0_num - 1];
    }

    printf("[GetFEICalFactor] error! unexpected decay ID\n");
    exit(1);
    return 0;

}

bool hasEnding(std::string const& fullString, std::string const& ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
    }
    else {
        return false;
    }
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

int Getpi0Multiplicity(double Upsilon_ID, double Bsig_ID) {
    if (-0.5 < Upsilon_ID && 0.5 > Upsilon_ID) { // charged
        if (-0.5 < Bsig_ID && 0.5 > Bsig_ID) return 0;
        else if (0.5 < Bsig_ID && 1.5 > Bsig_ID) return 1;
        else if (1.5 < Bsig_ID && 2.5 > Bsig_ID) return 0;
        else if (2.5 < Bsig_ID && 3.5 > Bsig_ID) return 0;
        else if (3.5 < Bsig_ID && 4.5 > Bsig_ID) return 1;
        else if (4.5 < Bsig_ID && 5.5 > Bsig_ID) return 1;
        else if (5.5 < Bsig_ID && 6.5 > Bsig_ID) return 0;
        else if (6.5 < Bsig_ID && 7.5 > Bsig_ID) return 0;
        else if (7.5 < Bsig_ID && 8.5 > Bsig_ID) return 1;
        else if (8.5 < Bsig_ID && 9.5 > Bsig_ID) return 2;
        else if (9.5 < Bsig_ID && 10.5 > Bsig_ID) return 2;
        else if (10.5 < Bsig_ID && 11.5 > Bsig_ID) return 2;
        else if (11.5 < Bsig_ID && 12.5 > Bsig_ID) return 0;
        else if (12.5 < Bsig_ID && 13.5 > Bsig_ID) return 0;
        else if (13.5 < Bsig_ID && 14.5 > Bsig_ID) return 1;
        else exit(1);
    }
    else if (0.5 < Upsilon_ID && 1.5 > Upsilon_ID) { // mixed
        if (-0.5 < Bsig_ID && 0.5 > Bsig_ID) return 0;
        else if (0.5 < Bsig_ID && 1.5 > Bsig_ID) return 0;
        else if (1.5 < Bsig_ID && 2.5 > Bsig_ID) return 1;
        else if (2.5 < Bsig_ID && 3.5 > Bsig_ID) return 1;
        else if (3.5 < Bsig_ID && 4.5 > Bsig_ID) return 0;
        else if (4.5 < Bsig_ID && 5.5 > Bsig_ID) return 0;
        else if (5.5 < Bsig_ID && 6.5 > Bsig_ID) return 1;
        else if (6.5 < Bsig_ID && 7.5 > Bsig_ID) return 1;
        else if (7.5 < Bsig_ID && 8.5 > Bsig_ID) return 0;
        else if (8.5 < Bsig_ID && 9.5 > Bsig_ID) return 2;
        else if (9.5 < Bsig_ID && 10.5 > Bsig_ID) return 2;
        else if (10.5 < Bsig_ID && 11.5 > Bsig_ID) return 2;
        else if (11.5 < Bsig_ID && 12.5 > Bsig_ID) return 0;
        else if (12.5 < Bsig_ID && 13.5 > Bsig_ID) return 0;
        else if (13.5 < Bsig_ID && 14.5 > Bsig_ID) return 1;
        else exit(1);
    }

    exit(1);
    return -1;
}

void LetsFillNgamma(const char* dirname, TH1D* hist_Ngamma_zero, TH1D* hist_Ngamma_one, TH1D* hist_Ngamma_two, std::string SampleName, int option, bool IsMultiplicityCorrectionApplied) {
    /*
    SampleName for Knn
    CHG
    MIX
    UUBAR
    DDBAR
    SSBAR
    Knn
    Kstarnn
    K0nn
    K0starnn
    */
    /*
    option 0: select all Btag
    option 1: select Btag+
    option 2: select Btag0
    */
    /*
    0: charged
    1: mixed
    2: uubar
    3: ddbar
    4: ssbar
    5: ccbar
    6: tautau
    7: mumu
    8: gg
    9: ee
    10: eeee
    11: eemumu
    12: llXX
    13: hhISR
    */

    double var = -1;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;
    double Btag_ID = -1;
    double temp_N_bin_PID[4][N_PID_syst] = { 0.0 }; // K-true, K-mis, pi-true, pi-miss
    double temp_N_bin_pi0[N_pi0_syst] = { 0.0 };
    double temp_N_bin_fakeE[4][N_fakeE_syst] = { 0.0 }; //  K-, K+, pi-, pi+
    double temp_N_bin_fakeMU[4][N_fakeMU_syst] = { 0.0 }; //  K-, K+, pi-, pi+

    double FEI_calibration_factor = -1;

    double invM_Knn = 0;
    double invM_Kstarnn = 0;
    double invM_K0nn = 0;
    double invM_K0starnn = 0;
    double N_Knn = 0;
    double N_Kstarnn = 0;
    double N_K0nn = 0;
    double N_K0starnn = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &var);

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID); // charged: 0, mixed: 1
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Btag->SetBranchAddress("Btag_extraInfo_decayModeID", &Btag_ID);
        for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKtruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[0][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKmisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[1][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npitruebin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[2][i_PID]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npimisbin" + std::to_string(i_PID)).c_str(), &temp_N_bin_PID[3][i_PID]);
        }
        for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npi0bin" + std::to_string(i_pi0)).c_str(), &temp_N_bin_pi0[i_pi0]);
        for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeEbin_n" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeE[0][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeEbin_p" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeE[1][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeEbin_n" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeE[2][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeEbin_p" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeE[3][i_fake]);
        }
        for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeMUbin_n" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeMU[0][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_nKfakeMUbin_p" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeMU[1][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeMUbin_n" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeMU[2][i_fake]);
            tree_Bsig->SetBranchAddress(("Bsig_daughter_0_extraInfo_npifakeMUbin_p" + std::to_string(i_fake)).c_str(), &temp_N_bin_fakeMU[3][i_fake]);
        }
        if ((SampleName == "CHG") || (SampleName == "MIX")) {
            tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKnn__bc", &N_Knn);
            tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKnn__bc", &invM_Knn);
            tree_upsilon->SetBranchAddress("nParticlesInList__boB__pl__clKstarnn__bc", &N_Kstarnn);
            tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstarnn__bc", &invM_Kstarnn);
            tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clK0nn__bc", &N_K0nn);
            tree_upsilon->SetBranchAddress("invMassInLists__bon0__clK0nn__bc", &invM_K0nn);
            tree_upsilon->SetBranchAddress("nParticlesInList__boB0__clKstar0nn__bc", &N_K0starnn);
            tree_upsilon->SetBranchAddress("invMassInLists__bon0__clKstar0nn__bc", &invM_K0starnn);
        }

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            if (option == 1 && Upsilon_ID != 0) continue;
            else if (option == 2 && Upsilon_ID != 1) continue;

            // Fill numberings
            double weight_ri = 0.0;
            if (SampleName == "CHG") {
                FEI_calibration_factor = GetFEICalFactor(Upsilon_ID, Btag_ID);
                weight_ri = ((N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab)); // total 2.8/ab for BB
            }
            else if (SampleName == "MIX") {
                FEI_calibration_factor = GetFEICalFactor(Upsilon_ID, Btag_ID);
                weight_ri = ((N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab)); // total 2.8/ab for BB
            }
            else if (SampleName == "UUBAR") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = ((0.364436 - 0.002763) / 1.0); // total 1.0/ab for qq
            }
            else if (SampleName == "DDBAR") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = ((0.364436 - 0.002763) / 1.0); // total 1.0/ab for qq
            }
            else if (SampleName == "SSBAR") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = ((0.364436 - 0.002763) / 1.0); // total 1.0/ab for qq
            }
            else if (SampleName == "CHARM") {
                FEI_calibration_factor = CAL_qq;
                weight_ri = ((0.364436 - 0.002763) / 1.0); // total 1.0/ab for qq
            }
            //else if (job_id >= 256846858 && job_id <= 256847295) numberings->push_back(6);
            //else if (job_id >= 256847296 && job_id <= 256847807) numberings->push_back(7);
            //else if (job_id >= 256847808 && job_id <= 256848291) numberings->push_back(8);
            //else if (job_id >= 256848292 && job_id <= 256848743) numberings->push_back(9);
            //else if (job_id >= 256848744 && job_id <= 256849128) numberings->push_back(10);
            //else if (job_id >= 256849129 && job_id <= 256849396) numberings->push_back(11);
            else {
                printf("undefined job id!\n");
                exit(1);
            }

            // Fill calibration factors
            double Correction_KID = 1;
            double Correction_PID = 1;
            double Correction_pi0 = 1;
            double Correction_fake = 1;
            for (int i_PID = 0; i_PID < N_PID_syst; i_PID++) {
                Correction_KID = Correction_KID * std::pow(PID_correction[0][i_PID], temp_N_bin_PID[0][i_PID]); // true KID
                Correction_KID = Correction_KID * std::pow(PID_correction[1][i_PID], temp_N_bin_PID[1][i_PID]); // mis KID
                Correction_PID = Correction_PID * std::pow(PID_correction[2][i_PID], temp_N_bin_PID[2][i_PID]); // true PID
                Correction_PID = Correction_PID * std::pow(PID_correction[3][i_PID], temp_N_bin_PID[3][i_PID]); // mis PID
            }
            for (int i_pi0 = 0; i_pi0 < N_pi0_syst; i_pi0++) Correction_pi0 = Correction_pi0 * std::pow(pi0_correction[i_pi0], temp_N_bin_pi0[i_pi0]);
            for (int i_fake = 0; i_fake < N_fakeE_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(PID_fakeE_correction[0][i_fake], temp_N_bin_fakeE[0][i_fake]); // K- from e
                Correction_fake = Correction_fake * std::pow(PID_fakeE_correction[1][i_fake], temp_N_bin_fakeE[1][i_fake]); // K+ from e
                Correction_fake = Correction_fake * std::pow(PID_fakeE_correction[2][i_fake], temp_N_bin_fakeE[2][i_fake]); // pi- from e
                Correction_fake = Correction_fake * std::pow(PID_fakeE_correction[3][i_fake], temp_N_bin_fakeE[3][i_fake]); // pi+ from e
            }
            for (int i_fake = 0; i_fake < N_fakeMU_syst; i_fake++) {
                Correction_fake = Correction_fake * std::pow(PID_fakeMU_correction[0][i_fake], temp_N_bin_fakeMU[0][i_fake]); // K- from mu
                Correction_fake = Correction_fake * std::pow(PID_fakeMU_correction[1][i_fake], temp_N_bin_fakeMU[1][i_fake]); // K+ from mu
                Correction_fake = Correction_fake * std::pow(PID_fakeMU_correction[2][i_fake], temp_N_bin_fakeMU[2][i_fake]); // pi- from mu
                Correction_fake = Correction_fake * std::pow(PID_fakeMU_correction[3][i_fake], temp_N_bin_fakeMU[3][i_fake]); // pi+ from mu
            }

            // Knn correction factor
            double Correction_Knn = corrector_Knn.GetCorrectionFactorAtGeneric(invM_Knn, invM_Kstarnn, invM_K0nn, invM_K0starnn, N_Knn, N_Kstarnn, N_K0nn, N_K0starnn);

            // Multiplicity correction factor, it is not applied now. it is for a systematic uncertainty
            double Ngamma_v200 = var;
            int pi0_multiplicity = Getpi0Multiplicity(Upsilon_ID, Bsig_ID);
            double Correction_multiplicity = 1.0;
            if (IsMultiplicityCorrectionApplied) Correction_multiplicity = corrector_Multiplicity->GetCorrectionFactor(Ngamma_v200, pi0_multiplicity);

            double weight = 1.0;
            if (IsMultiplicityCorrectionApplied) weight = FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn * Correction_multiplicity;
            else weight = FEI_calibration_factor * CAL * weight_ri * Correction_pi0 * Correction_KID * Correction_PID * Correction_fake * Correction_Knn;

            if (pi0_multiplicity == 0) hist_Ngamma_zero->Fill(var, weight);
            else if (pi0_multiplicity == 1) hist_Ngamma_one->Fill(var, weight);
            else if (pi0_multiplicity == 2) hist_Ngamma_two->Fill(var, weight);


        }
        input_file->Close();

    }

}

void LetsFillNgamma(const char* dirname, TH1D* hist_Ngamma_zero, TH1D* hist_Ngamma_one, TH1D* hist_Ngamma_two, int option) {
    /*
    option 0: select all Btag
    option 1: select Btag+
    option 2: select Btag0
    */
    /*
    0: charged
    1: mixed
    2: uubar
    3: ddbar
    4: ssbar
    5: ccbar
    6: tautau
    7: mumu
    8: gg
    9: ee
    10: eeee
    11: eemumu
    12: llXX
    13: hhISR
    */

    double var = -1;

    double Upsilon_ID = -1;
    double Bsig_ID = -1;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("extraInfo__boNgammav200__bc", &var);

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID); // charged: 0, mixed: 1
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            if (option == 1 && Upsilon_ID != 0) continue;
            else if (option == 2 && Upsilon_ID != 1) continue;

            hist_Ngamma->Fill(var);

        }
        input_file->Close();

    }

}

void MultiplicityCalculator_v2(){

    ReadPIDFile();
    ReadFakePIDFile();
    
    const char* ChargeSideband_MC_CHG_train_dirname = "/home/jwpark/storage/BKG_gbasf2/Kokoro_LS_MC_cside/CHG_analysis/train_v000/final_output";
    const char* ChargeSideband_MC_MIX_train_dirname = "/home/jwpark/storage/BKG_gbasf2/Kokoro_LS_MC_cside/MIX_analysis/train_v000/final_output";
    const char* ChargeSideband_MC_UUBAR_train_dirname = "/home/jwpark/storage/BKG_gbasf2/Kokoro_LS_MC_cside/UUBAR_analysis/train_v000/final_output";
    const char* ChargeSideband_MC_DDBAR_train_dirname = "/home/jwpark/storage/BKG_gbasf2/Kokoro_LS_MC_cside/DDBAR_analysis/train_v000/final_output";
    const char* ChargeSideband_MC_SSBAR_train_dirname = "/home/jwpark/storage/BKG_gbasf2/Kokoro_LS_MC_cside/SSBAR_analysis/train_v000/final_output";
    const char* ChargeSideband_MC_CHARM_train_dirname = "/home/jwpark/storage/BKG_gbasf2/Kokoro_LS_MC_cside/CHARM_analysis/train_v000/final_output";

    const char* ChargeSideband_MC_CHG_test_dirname = "/home/jwpark/storage/BKG_gbasf2/Kokoro_LS_MC_cside/CHG_analysis/test_v000/final_output";
    const char* ChargeSideband_MC_MIX_test_dirname = "/home/jwpark/storage/BKG_gbasf2/Kokoro_LS_MC_cside/MIX_analysis/test_v000/final_output";
    const char* ChargeSideband_MC_UUBAR_test_dirname = "/home/jwpark/storage/BKG_gbasf2/Kokoro_LS_MC_cside/UUBAR_analysis/test_v000/final_output";
    const char* ChargeSideband_MC_DDBAR_test_dirname = "/home/jwpark/storage/BKG_gbasf2/Kokoro_LS_MC_cside/DDBAR_analysis/test_v000/final_output";
    const char* ChargeSideband_MC_SSBAR_test_dirname = "/home/jwpark/storage/BKG_gbasf2/Kokoro_LS_MC_cside/SSBAR_analysis/test_v000/final_output";
    const char* ChargeSideband_MC_CHARM_test_dirname = "/home/jwpark/storage/BKG_gbasf2/Kokoro_LS_MC_cside/CHARM_analysis/test_v000/final_output";

    const char* ChargeSideband_data_dirname = "/home/jwpark/storage/BKG_gbasf2/Kokoro_LS_data_cside/SIGNAL_analysis/validation_v000/final_output";

    // Lets fill!
    TH1D* Ngamma_v200_MC = new TH1D("Ngamma_v200_MC", ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);
    TH1D* Ngamma_v200_data = new TH1D("Ngamma_v200_data", ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);
    TH1D* Ngamma_zero_v200_MC = new TH1D("Ngamma_zero_v200_MC", ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);
    TH1D* Ngamma_zero_v200_data = new TH1D("Ngamma_zero_v200_data", ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);
    TH1D* Ngamma_one_v200_MC = new TH1D("Ngamma_one_v200_MC", ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);
    TH1D* Ngamma_one_v200_data = new TH1D("Ngamma_one_v200_data", ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);
    TH1D* Ngamma_two_v200_MC = new TH1D("Ngamma_two_v200_MC", ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);
    TH1D* Ngamma_two_v200_data = new TH1D("Ngamma_two_v200_data", ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);

    LetsFillNgamma(ChargeSideband_MC_CHG_train_dirname, Ngamma_zero_v200_MC, Ngamma_one_v200_MC, Ngamma_two_v200_MC, "CHG", 1, false);
    LetsFillNgamma(ChargeSideband_MC_MIX_train_dirname, Ngamma_zero_v200_MC, Ngamma_one_v200_MC, Ngamma_two_v200_MC, "MIX", 1, false);
    LetsFillNgamma(ChargeSideband_MC_UUBAR_train_dirname, Ngamma_zero_v200_MC, Ngamma_one_v200_MC, Ngamma_two_v200_MC, "UUBAR", 1, false);
    LetsFillNgamma(ChargeSideband_MC_DDBAR_train_dirname, Ngamma_zero_v200_MC, Ngamma_one_v200_MC, Ngamma_two_v200_MC, "DDBAR", 1, false);
    LetsFillNgamma(ChargeSideband_MC_SSBAR_train_dirname, Ngamma_zero_v200_MC, Ngamma_one_v200_MC, Ngamma_two_v200_MC, "SSBAR", 1, false);
    LetsFillNgamma(ChargeSideband_MC_CHARM_train_dirname, Ngamma_zero_v200_MC, Ngamma_one_v200_MC, Ngamma_two_v200_MC, "CHARM", 1, false);
    LetsFillNgamma(ChargeSideband_MC_CHG_test_dirname, Ngamma_zero_v200_MC, Ngamma_one_v200_MC, Ngamma_two_v200_MC, "CHG", 1, false);
    LetsFillNgamma(ChargeSideband_MC_MIX_test_dirname, Ngamma_zero_v200_MC, Ngamma_one_v200_MC, Ngamma_two_v200_MC, "MIX", 1, false);
    LetsFillNgamma(ChargeSideband_MC_UUBAR_test_dirname, Ngamma_zero_v200_MC, Ngamma_one_v200_MC, Ngamma_two_v200_MC, "UUBAR", 1, false);
    LetsFillNgamma(ChargeSideband_MC_DDBAR_test_dirname, Ngamma_zero_v200_MC, Ngamma_one_v200_MC, Ngamma_two_v200_MC, "DDBAR", 1, false);
    LetsFillNgamma(ChargeSideband_MC_SSBAR_test_dirname, Ngamma_zero_v200_MC, Ngamma_one_v200_MC, Ngamma_two_v200_MC, "SSBAR", 1, false);
    LetsFillNgamma(ChargeSideband_MC_CHARM_test_dirname, Ngamma_zero_v200_MC, Ngamma_one_v200_MC, Ngamma_two_v200_MC, "CHARM", 1, false);

    LetsFillNgamma(ChargeSideband_data_dirname, Ngamma_zero_v200_data, Ngamma_one_v200_data, Ngamma_two_v200_data, 1);

    // print weight file
    FILE* fp = fopen("multiplicity_weight_zero.txt","w");
    fprintf(fp, "%d\n", NgammaMAX);
    for (int i = 0; i < NgammaMAX + 1; i++) {
        double MC_num_bin = Ngamma_zero_v200_MC->GetBinContent(i + 1);
        double data_num_bin = Ngamma_zero_v200_data->GetBinContent(i + 1);
        if (MC_num_bin > MyEPSILON) fprintf(fp, "%lf %lf %lf\n", data_num_bin, MC_num_bin, data_num_bin / MC_num_bin);
        else fprintf(fp, "%lf %lf %lf\n", data_num_bin, MC_num_bin, 1.0);
    }
    fclose(fp);

    fp = fopen("multiplicity_weight_one.txt", "w");
    fprintf(fp, "%d\n", NgammaMAX);
    for (int i = 0; i < NgammaMAX + 1; i++) {
        double MC_num_bin = Ngamma_one_v200_MC->GetBinContent(i + 1);
        double data_num_bin = Ngamma_one_v200_data->GetBinContent(i + 1);
        if (MC_num_bin > MyEPSILON) fprintf(fp, "%lf %lf %lf\n", data_num_bin, MC_num_bin, data_num_bin / MC_num_bin);
        else fprintf(fp, "%lf %lf %lf\n", data_num_bin, MC_num_bin, 1.0);
    }
    fclose(fp);

    fp = fopen("multiplicity_weight_two.txt", "w");
    fprintf(fp, "%d\n", NgammaMAX);
    for (int i = 0; i < NgammaMAX + 1; i++) {
        double MC_num_bin = Ngamma_two_v200_MC->GetBinContent(i + 1);
        double data_num_bin = Ngamma_two_v200_data->GetBinContent(i + 1);
        if (MC_num_bin > MyEPSILON) fprintf(fp, "%lf %lf %lf\n", data_num_bin, MC_num_bin, data_num_bin / MC_num_bin);
        else fprintf(fp, "%lf %lf %lf\n", data_num_bin, MC_num_bin, 1.0);
    }
    fclose(fp);

    double MC_num = 0;
    double data_num = 0;

    for (int i = 0; i < NgammaMAX + 1; i++) {
        Ngamma_v200_MC->SetBinContent(i + 1, Ngamma_zero_v200_MC->GetBinContent(i + 1) + Ngamma_one_v200_MC->GetBinContent(i + 1) + Ngamma_two_v200_MC->GetBinContent(i + 1));
        Ngamma_v200_data->SetBinContent(i + 1, Ngamma_zero_v200_data->GetBinContent(i + 1) + Ngamma_one_v200_data->GetBinContent(i + 1) + Ngamma_two_v200_data->GetBinContent(i + 1));
        MC_num = MC_num + Ngamma_v200_MC->GetBinContent(i + 1);
        data_num = data_num + Ngamma_v200_data->GetBinContent(i + 1);
    }

    printf("option1 data num: %lf\n", data_num);
    printf("option1 MC num with calibration: %lf\n", MC_num);




    corrector_Multiplicity = new Corrector_Multiplicity();

    // Lets fill!
    TH1D* Ngamma_v200_MC_test = new TH1D("Ngamma_v200_MC_test", ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);
    TH1D* Ngamma_v200_data_test = new TH1D("Ngamma_v200_data_test", ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);
    TH1D* Ngamma_zero_v200_MC_test = new TH1D("Ngamma_zero_v200_MC_test", ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);
    TH1D* Ngamma_zero_v200_data_test = new TH1D("Ngamma_zero_v200_data_test", ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);
    TH1D* Ngamma_one_v200_MC_test = new TH1D("Ngamma_one_v200_MC_test", ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);
    TH1D* Ngamma_one_v200_data_test = new TH1D("Ngamma_one_v200_data_test", ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);
    TH1D* Ngamma_two_v200_MC_test = new TH1D("Ngamma_two_v200_MC_test", ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);
    TH1D* Ngamma_two_v200_data_test = new TH1D("Ngamma_two_v200_data_test", ";;", NgammaMAX + 1, -0.5, NgammaMAX + 0.5);

    LetsFillNgamma(ChargeSideband_MC_CHG_train_dirname, Ngamma_zero_v200_MC_test, Ngamma_one_v200_MC_test, Ngamma_two_v200_MC_test, "CHG", 2, true);
    LetsFillNgamma(ChargeSideband_MC_MIX_train_dirname, Ngamma_zero_v200_MC_test, Ngamma_one_v200_MC_test, Ngamma_two_v200_MC_test, "MIX", 2, true);
    LetsFillNgamma(ChargeSideband_MC_UUBAR_train_dirname, Ngamma_zero_v200_MC_test, Ngamma_one_v200_MC_test, Ngamma_two_v200_MC_test, "UUBAR", 2, true);
    LetsFillNgamma(ChargeSideband_MC_DDBAR_train_dirname, Ngamma_zero_v200_MC_test, Ngamma_one_v200_MC_test, Ngamma_two_v200_MC_test, "DDBAR", 2, true);
    LetsFillNgamma(ChargeSideband_MC_SSBAR_train_dirname, Ngamma_zero_v200_MC_test, Ngamma_one_v200_MC_test, Ngamma_two_v200_MC_test, "SSBAR", 2, true);
    LetsFillNgamma(ChargeSideband_MC_CHARM_train_dirname, Ngamma_zero_v200_MC_test, Ngamma_one_v200_MC_test, Ngamma_two_v200_MC_test, "CHARM", 2, true);
    LetsFillNgamma(ChargeSideband_MC_CHG_test_dirname, Ngamma_zero_v200_MC_test, Ngamma_one_v200_MC_test, Ngamma_two_v200_MC_test, "CHG", 2, true);
    LetsFillNgamma(ChargeSideband_MC_MIX_test_dirname, Ngamma_zero_v200_MC_test, Ngamma_one_v200_MC_test, Ngamma_two_v200_MC_test, "MIX", 2, true);
    LetsFillNgamma(ChargeSideband_MC_UUBAR_test_dirname, Ngamma_zero_v200_MC_test, Ngamma_one_v200_MC_test, Ngamma_two_v200_MC_test, "UUBAR", 2, true);
    LetsFillNgamma(ChargeSideband_MC_DDBAR_test_dirname, Ngamma_zero_v200_MC_test, Ngamma_one_v200_MC_test, Ngamma_two_v200_MC_test, "DDBAR", 2, true);
    LetsFillNgamma(ChargeSideband_MC_SSBAR_test_dirname, Ngamma_zero_v200_MC_test, Ngamma_one_v200_MC_test, Ngamma_two_v200_MC_test, "SSBAR", 2, true);
    LetsFillNgamma(ChargeSideband_MC_CHARM_test_dirname, Ngamma_zero_v200_MC_test, Ngamma_one_v200_MC_test, Ngamma_two_v200_MC_test, "CHARM", 2, true);

    LetsFillNgamma(ChargeSideband_data_dirname, Ngamma_zero_v200_data_test, Ngamma_one_v200_data_test, Ngamma_two_v200_data_test, 2);

    // print weight file
    fp = fopen("multiplicity_weight_uncertainty_zero.txt", "w");
    fprintf(fp, "%d\n", NgammaMAX);
    for (int i = 0; i < NgammaMAX + 1; i++) {
        double MC_num_bin = Ngamma_zero_v200_MC_test->GetBinContent(i + 1);
        double data_num_bin = Ngamma_zero_v200_data_test->GetBinContent(i + 1);
        if (MC_num_bin > MyEPSILON) fprintf(fp, "%lf %lf %lf\n", data_num_bin, MC_num_bin, data_num_bin / MC_num_bin);
        else fprintf(fp, "%lf %lf %lf\n", data_num_bin, MC_num_bin, 1.0);
    }
    fclose(fp);

    fp = fopen("multiplicity_weight_uncertainty_one.txt", "w");
    fprintf(fp, "%d\n", NgammaMAX);
    for (int i = 0; i < NgammaMAX + 1; i++) {
        double MC_num_bin = Ngamma_one_v200_MC_test->GetBinContent(i + 1);
        double data_num_bin = Ngamma_one_v200_data_test->GetBinContent(i + 1);
        if (MC_num_bin > MyEPSILON) fprintf(fp, "%lf %lf %lf\n", data_num_bin, MC_num_bin, data_num_bin / MC_num_bin);
        else fprintf(fp, "%lf %lf %lf\n", data_num_bin, MC_num_bin, 1.0);
    }
    fclose(fp);

    fp = fopen("multiplicity_weight_uncertainty_two.txt", "w");
    fprintf(fp, "%d\n", NgammaMAX);
    for (int i = 0; i < NgammaMAX + 1; i++) {
        double MC_num_bin = Ngamma_two_v200_MC_test->GetBinContent(i + 1);
        double data_num_bin = Ngamma_two_v200_data_test->GetBinContent(i + 1);
        if (MC_num_bin > MyEPSILON) fprintf(fp, "%lf %lf %lf\n", data_num_bin, MC_num_bin, data_num_bin / MC_num_bin);
        else fprintf(fp, "%lf %lf %lf\n", data_num_bin, MC_num_bin, 1.0);
    }
    fclose(fp);

    MC_num = 0;
    data_num = 0;

    for (int i = 0; i < NgammaMAX + 1; i++) {
        Ngamma_v200_MC_test->SetBinContent(i + 1, Ngamma_zero_v200_MC_test->GetBinContent(i + 1) + Ngamma_one_v200_MC_test->GetBinContent(i + 1) + Ngamma_two_v200_MC_test->GetBinContent(i + 1));
        Ngamma_v200_data_test->SetBinContent(i + 1, Ngamma_zero_v200_data_test->GetBinContent(i + 1) + Ngamma_one_v200_data_test->GetBinContent(i + 1) + Ngamma_two_v200_data_test->GetBinContent(i + 1));
        MC_num = MC_num + Ngamma_v200_MC_test->GetBinContent(i + 1);
        data_num = data_num + Ngamma_v200_data_test->GetBinContent(i + 1);
    }

    printf("option2 data num: %lf\n", data_num);
    printf("option2 MC num with calibration: %lf\n", MC_num);

}
