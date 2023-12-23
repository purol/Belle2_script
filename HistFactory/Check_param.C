#include <iostream>
#include <stdio.h>

#include "TRandom3.h"
#include "TCanvas.h"
#include "TDatime.h"
#include "TStopwatch.h"
#include "TLegend.h"
#include "TIterator.h"
#include "TH3.h"
#include "TLatex.h"

#include "RooChi2Var.h"
#include "RooAbsData.h"
#include "RooRealSumPdf.h"
#include "RooPoisson.h"
#include "RooGaussian.h"
#include "RooRealVar.h"
#include "RooMCStudy.h"
#include "RooMinuit.h"
#include "RooCategory.h"
#include "RooHistPdf.h"
#include "RooSimultaneous.h"
#include "RooExtendPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "RooParamHistFunc.h"
#include "RooHist.h"
#include "RooRandom.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/MinNLLTestStat.h"

#include "RooStats/HistFactory/FlexibleInterpVar.h"
#include "RooStats/HistFactory/PiecewiseInterpolation.h"
#include "RooStats/HistFactory/HistFactorySimultaneous.h"
#include "RooStats/HistFactory/Channel.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/ParamHistFunc.h"
#include "RooStats/HistFactory/HistFactoryModelUtils.h"
#include "RooStats/HistFactory/RooBarlowBeestonLL.h"

using namespace RooFit;
using namespace RooStats;
using namespace HistFactory;

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

# define FEI_cal_Bc_num 12
# define FEI_cal_B0_num 11
double FEI_cal_Bc[FEI_cal_Bc_num] = { 1.04, 0.79, 0.69, 0.56, 0.97, 0.95, 0.74, 0.57, 0.91, 0.51, 0.34, 0.59 };
double FEI_cal_Bc_uncertainty[FEI_cal_Bc_num] = { 0.03, 0.03, 0.05, 0.11, 0.03, 0.03, 0.02, 0.06, 0.1, 0.13, 0.07, 0.02 }; // not relative uncertainty. absolute uncertainty
double FEI_cal_Bc_modeID[FEI_cal_Bc_num] = { 0.0, 1.0, 3.0, 4.0, 15.0, 16.0, 18.0, 19.0, 23.0, 24.0, 30.0, -1.0 };
double FEI_cal_B0[FEI_cal_B0_num] = { 1.16, 0.94, 0.81, 0.79, 0.99, 1.03, 0.67, 0.66, 0.69, 0.49, 0.79 };
double FEI_cal_B0_uncertainty[FEI_cal_B0_num] = { 0.04, 0.05, 0.06, 0.02, 0.03, 0.06, 0.02, 0.03, 0.02, 0.02, 0.12 }; // not relative uncertainty. absolute uncertainty
double FEI_cal_B0_modeID[FEI_cal_B0_num] = { 0.0, 1.0, 3.0, 4.0, 5.0, 15.0, 16.0, 18.0, 19.0, 26.0, -1.0 };

std::vector<std::string> Sample_names = {
    "L_x_Signal_nominal_channel_overallSyst_x_StatUncert_x_channel_Signal_all_uncorr_uncer_ShapeSys",
    "L_x_CHG_nominal_channel_overallSyst_x_StatUncert_x_channel_CHG_all_uncorr_uncer_ShapeSys",
    "L_x_MIX_nominal_channel_overallSyst_x_StatUncert_x_channel_MIX_all_uncorr_uncer_ShapeSys",
    "L_x_UUBAR_nominal_channel_overallSyst_x_StatUncert_x_channel_UUBAR_all_uncorr_uncer_ShapeSys",
    "L_x_DDBAR_nominal_channel_overallSyst_x_StatUncert_x_channel_DDBAR_all_uncorr_uncer_ShapeSys",
    "L_x_SSBAR_nominal_channel_overallSyst_x_StatUncert_x_channel_SSBAR_all_uncorr_uncer_ShapeSys",
    "L_x_CHARM_nominal_channel_overallSyst_x_StatUncert_x_channel_CHARM_all_uncorr_uncer_ShapeSys"
};

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

# define RarityBins 6

double weight_KIDsys[RarityBins * 7] = { 0.0 };
double weight_PIDsys[RarityBins * 7] = { 0.0 };
double weight_BRsys[RarityBins * 3] = { 0.0 };
double weight_pi0sys[RarityBins * 7] = { 0.0 };
double weight_FEIsys[RarityBins * 3] = { 0.0 };

std::vector<std::string> split(std::string str, char Delimiter) {
    std::istringstream iss(str);
    std::string buffer;

    std::vector<std::string> result;

    while (getline(iss, buffer, Delimiter)) {
        result.push_back(buffer);
    }

    return result;
}

typedef struct Options
{
    // turn on it?
    bool track;
    bool PID;
    bool KID;
    bool KS0;
    bool pi0;
    bool FBDT;
    bool FEI;
    bool qqbar;
    bool multiplicity;
    bool Kff;
    bool Kstarff;
    bool pf;
    bool Transition;
    bool mb;
    bool fraction;
    bool MCstat;
    bool Fragmentation;
    bool mKstar;
    bool BKGNorm;
    bool BDTc;
    bool BBcounting;
    bool BBBR;
    bool BRKnn;
    bool BRDKL0;
    bool uncorrelated;

    int NEntryFEI;
    int NEntryKID;
    int NEntryPID;
    int NEntryBR;
    int NEntrypi0;
    int NEntryMultiplicity;

    int Nsyst = 25;
} OPTIONS;


int ReadNFEIEigenVector(const char* dirname) {
    int Nentry = 0; // number of eigen values/vectors
    double eigen_value = 0; // eigen value
    double weight_sys[RarityBins * 3] = { 0.0 }; // eigen vector

    FILE* fp;
    fp = fopen(dirname, "r");
    while (true) {
        if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
        for (int i = 0; i < RarityBins * 3; i++) {
            if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
        }
        Nentry++;
    }
    fclose(fp);

    return Nentry;
}

int ReadNPIDEigenVector(const char* dirname) {
    int Nentry = 0; // number of eigen values/vectors
    double eigen_value = 0; // eigen value
    double weight_sys[RarityBins * 7] = { 0.0 }; // eigen vector

    FILE* fp;
    fp = fopen(dirname, "r");
    while (true) {
        if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
        for (int i = 0; i < RarityBins * 7; i++) {
            if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
        }
        Nentry++;
    }
    fclose(fp);

    return Nentry;
}

int ReadNBREigenVector(const char* dirname) {
    int Nentry = 0; // number of eigen values/vectors
    double eigen_value = 0; // eigen value
    double weight_sys[RarityBins * 3] = { 0.0 }; // eigen vector

    FILE* fp;
    fp = fopen(dirname, "r");
    while (true) {
        if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
        for (int i = 0; i < RarityBins * 3; i++) {
            if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
        }
        Nentry++;
    }
    fclose(fp);

    return Nentry;
}

int ReadNpi0EigenVector(const char* dirname) {
    int Nentry = 0; // number of eigen values/vectors
    double eigen_value = 0; // eigen value
    double weight_sys[RarityBins * 7] = { 0.0 }; // eigen vector

    FILE* fp;
    fp = fopen(dirname, "r");
    while (true) {
        if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
        for (int i = 0; i < RarityBins * 7; i++) {
            if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
        }
        Nentry++;
    }
    fclose(fp);

    return Nentry;
}

int ReadMultiplicityInfo(const char* dirname) {

    int N_hist = -1;

    FILE* fp;
    fp = fopen(dirname, "r");
    fscanf(fp, "%d", &N_hist);
    fclose(fp);

    return N_hist;
}

void Initialize_options(OPTIONS* options_, const char* tested_param) {
    // initialize
    options_->track = false;
    options_->PID = false;
    options_->KID = false;
    options_->KS0 = false;
    options_->pi0 = false;
    options_->FBDT = false;
    options_->FEI = false;
    options_->qqbar = false;
    options_->multiplicity = false;
    options_->Kff = false;
    options_->Kstarff = false;
    options_->pf = false;
    options_->Transition = false;
    options_->mb = false;
    options_->fraction = false;
    options_->MCstat = false;
    options_->Fragmentation = false;
    options_->mKstar = false;
    options_->BKGNorm = false;
    options_->BDTc = false;
    options_->BBcounting = false;
    options_->BBBR = false;
    options_->BRKnn = false;
    options_->BRDKL0 = false;
    options_->uncorrelated = false;

    if (std::string(tested_param) == std::string("all")) {
        options_->track = true;
        options_->PID = true;
        options_->KID = true;
        options_->KS0 = true;
        options_->pi0 = true;
        options_->FBDT = true;
        options_->FEI = true;
        options_->qqbar = true;
        options_->multiplicity = true;
        options_->Kff = true;
        options_->Kstarff = true;
        options_->pf = true;
        options_->Transition = true;
        options_->mb = true;
        options_->fraction = true;
        options_->MCstat = true;
        options_->Fragmentation = true;
        options_->mKstar = true;
        options_->BKGNorm = true;
        options_->BDTc = true;
        options_->BBcounting = true;
        options_->BBBR = true;
        options_->BRKnn = true;
        options_->BRDKL0 = true;
        options_->uncorrelated = true;
    }
    else if (std::string(tested_param) == std::string("none")) {}
    else if (std::string(tested_param) == std::string("track")) options_->track = true;
    else if (std::string(tested_param) == std::string("PID")) options_->PID = true;
    else if (std::string(tested_param) == std::string("KID")) options_->KID = true;
    else if (std::string(tested_param) == std::string("KS0")) options_->KS0 = true;
    else if (std::string(tested_param) == std::string("pi0")) options_->pi0 = true;
    else if (std::string(tested_param) == std::string("FBDT_efficiency")) options_->FBDT = true;
    else if (std::string(tested_param) == std::string("FEI")) options_->FEI = true;
    else if (std::string(tested_param) == std::string("qqbar_CAL")) options_->qqbar = true;
    else if (std::string(tested_param) == std::string("photon_multiplicity")) options_->multiplicity = true;
    else if (std::string(tested_param) == std::string("Kff")) options_->Kff = true;
    else if (std::string(tested_param) == std::string("Kstarff")) options_->Kstarff = true;
    else if (std::string(tested_param) == std::string("pf")) options_->pf = true;
    else if (std::string(tested_param) == std::string("Transition")) options_->Transition = true;
    else if (std::string(tested_param) == std::string("mb")) options_->mb = true;
    else if (std::string(tested_param) == std::string("fraction")) options_->fraction = true;
    else if (std::string(tested_param) == std::string("MCstat")) options_->MCstat = true;
    else if (std::string(tested_param) == std::string("Fragmentation")) options_->Fragmentation = true;
    else if (std::string(tested_param) == std::string("mKstar")) options_->mKstar = true;
    else if (std::string(tested_param) == std::string("BKGNorm")) options_->BKGNorm = true;
    else if (std::string(tested_param) == std::string("BDTc")) options_->BDTc = true;
    else if (std::string(tested_param) == std::string("BBcounting")) options_->BBcounting = true;
    else if (std::string(tested_param) == std::string("BBBR")) options_->BBBR = true;
    else if (std::string(tested_param) == std::string("BRKnn")) options_->BRKnn = true;
    else if (std::string(tested_param) == std::string("BtoDtoXKL")) options_->BRDKL0 = true;
    else if (std::string(tested_param) == std::string("uncorrelated")) options_->uncorrelated = true;
    else {
        printf("inappropriate parameter name\n");
        exit(1);
    }

    // read entry for nuisance parameters
    options_->NEntryFEI = ReadNFEIEigenVector("./FEI_selected.txt");
    options_->NEntryKID = ReadNPIDEigenVector("./KID_selected.txt");
    options_->NEntryPID = ReadNPIDEigenVector("./PID_selected.txt");
    options_->NEntryBR = ReadNBREigenVector("./BR_selected.txt");
    options_->NEntrypi0 = ReadNpi0EigenVector("./pi0_selected.txt");
    options_->NEntryMultiplicity = ReadMultiplicityInfo("./Multiplicity_info.txt");
}

void FixParameters(RooWorkspace* w, OPTIONS* options_) {

    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    w->loadSnapshot("NominalParamValues");

    RooRealVar* x_val = w->var("obs_x_channel");
    std::unique_ptr<RooArgSet> params{model->getParameters(*x_val)};

    // track
    w->var("alpha_track_eff_uncer")->setConstant(options_->track);

    // PID
    for (int i = 0; i < options_->NEntryPID; i++) w->var(("alpha_PID" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->PID);

    // KID
    for (int i = 0; i < options_->NEntryKID; i++) w->var(("alpha_KID" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->KID);

    // KS0
    w->var("alpha_KS0_reco_uncer")->setConstant(options_->KS0);

    // pi0
    for (int i = 0; i < options_->NEntrypi0; i++) w->var(("alpha_pi0" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->pi0);

    // FBDT efficiency
    w->var("alpha_FBDT_efficiency_uncer")->setConstant(options_->FBDT);

    // FEI
    for (int i = 0; i < options_->NEntryFEI; i++) w->var(("alpha_FEI" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->FEI);

    // qqbar normalization
    w->var("alpha_qq_CAL_UUBAR_uncer")->setConstant(options_->qqbar);
    w->var("alpha_qq_CAL_DDBAR_uncer")->setConstant(options_->qqbar);
    w->var("alpha_qq_CAL_SSBAR_uncer")->setConstant(options_->qqbar);
    w->var("alpha_qq_CAL_CHARM_uncer")->setConstant(options_->qqbar);

    // photon multiplicity correction
    for (int i = 0; i < options_->NEntryMultiplicity; i++) w->var(("alpha_multiplicity" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->multiplicity);

    // B->K form factor
    w->var("alpha_Kff1_uncer")->setConstant(options_->Kff);
    w->var("alpha_Kff2_uncer")->setConstant(options_->Kff);
    w->var("alpha_Kff3_uncer")->setConstant(options_->Kff);

    // B->Kstar form factor
    w->var("alpha_Kstarff1_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff2_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff3_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff4_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff5_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff6_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff7_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff8_uncer")->setConstant(options_->Kstarff);
    w->var("alpha_Kstarff9_uncer")->setConstant(options_->Kstarff);

    // fermi motion moment
    w->var("alpha_pf_uncer")->setConstant(options_->pf);

    // K*-Xs transition point
    w->var("alpha_transition_uncer")->setConstant(options_->Transition);

    // b-quark mass
    w->var("alpha_mb_uncer")->setConstant(options_->mb);

    // relative fraction
    w->var("alpha_Kfrac_uncer")->setConstant(options_->fraction);
    w->var("alpha_Kstarfrac_uncer")->setConstant(options_->fraction);

    // MC statistics
    for (int i = 0; i < RarityBins; i++) w->var(("gamma_stat_channel_bin_" + std::to_string(i)).c_str())->setConstant(options_->MCstat);

    // Fragmentation
    for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
        for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
            int temp_index = MxsBin * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma) + Category;
            w->var(("alpha_Xs_fragmentation" + std::to_string(temp_index) + "_uncer").c_str())->setConstant(options_->Fragmentation);
        }
    }

    // mKstar
    w->var("alpha_mKstar_uncer")->setConstant(options_->mKstar);

    // background normalization
    w->var("alpha_mu_CHG")->setConstant(options_->BKGNorm);
    w->var("alpha_mu_MIX")->setConstant(options_->BKGNorm);
    w->var("alpha_mu_UUBAR")->setConstant(options_->BKGNorm);
    w->var("alpha_mu_DDBAR")->setConstant(options_->BKGNorm);
    w->var("alpha_mu_SSBAR")->setConstant(options_->BKGNorm);
    w->var("alpha_mu_CHARM")->setConstant(options_->BKGNorm);

    // BDTc
    w->var("alpha_BDTc_shape_qq")->setConstant(options_->BDTc);

    // BB counting
    w->var("alpha_BB_counting_uncer")->setConstant(options_->BBcounting);

    // BB BR
    for (int i = 0; i < options_->NEntryBR; i++) w->var(("alpha_BBBR" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->BBBR);

    // B->Knn BR
    w->var("alpha_Knn_BR_uncer")->setConstant(options_->BRKnn);
    w->var("alpha_Kstarnn_BR_uncer")->setConstant(options_->BRKnn);
    w->var("alpha_K0nn_BR_uncer")->setConstant(options_->BRKnn);
    w->var("alpha_K0starnn_BR_uncer")->setConstant(options_->BRKnn);

    // B->(D->X KL0) + anything
    w->var("alpha_BtoDtoXKL_uncer")->setConstant(options_->BRDKL0);

    // uncorrelated
    for (int i = 0; i < RarityBins; i++) {
        w->var(("gamma_Signal_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
        w->var(("gamma_CHG_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
        w->var(("gamma_MIX_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
        w->var(("gamma_UUBAR_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
        w->var(("gamma_DDBAR_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
        w->var(("gamma_SSBAR_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
        w->var(("gamma_CHARM_all_uncorr_uncer_bin_" + std::to_string(i)).c_str())->setConstant(options_->uncorrelated);
    }

    // save snapshot
    w->saveSnapshot("ParamValues", *params, true);
}

RooFitResult* MyMinimizeNLL(RooWorkspace* w, RooDataSet* data, RooAbsReal** nll, double tolerance = -1.0) {
    // what we have done
    w->loadSnapshot("ParamValues");
    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    // get nll
    RooArgSet* allParams = model->getParameters(*data);
    RooStats::RemoveConstantParameters(allParams);
    RooArgSet fGlobalObs = *mc->GetGlobalObservables();
    RooArgSet fConditionalObs;
    Bool_t fLOffset = RooStats::IsNLLOffset();
    (*nll) = model->createNLL(*data, RooFit::CloneData(kFALSE), RooFit::Constrain(*allParams), RooFit::GlobalObservables(fGlobalObs), RooFit::ConditionalObservables(fConditionalObs), RooFit::Offset(fLOffset));

    // minimizer option
    TString fMinimizer = ::ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str();
    TString minimizer = fMinimizer;

    TString algorithm = ::ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo();

    Int_t fStrategy = ::ROOT::Math::MinimizerOptions::DefaultStrategy();

    Double_t fTolerance;
    if (tolerance < 0) fTolerance = TMath::Max(1., ::ROOT::Math::MinimizerOptions::DefaultTolerance());
    else fTolerance = tolerance;

    Int_t fPrintLevel = ::ROOT::Math::MinimizerOptions::DefaultPrintLevel();
    //LM: RooMinimizer.setPrintLevel has +1 offset - so subtract  here -1 + an extra -1
    int level = (fPrintLevel == 0) ? -1 : fPrintLevel - 2;


    // follow what ProfileLikelihoodTestStat.cxx does
    const auto& config = RooStats::GetGlobalRooStatsConfig();
    RooMinimizer minim(*(*nll));
    minim.setStrategy(fStrategy);
    minim.setEvalErrorWall(config.useEvalErrorWall);
    minim.setEps(fTolerance);
    minim.setPrintLevel(level);
    // this causes a memory leak
    minim.optimizeConst(2);
    minim.migrad();
    minim.minos(RooArgSet(*w->var("mu")));

    // fit!
    int status;
    status = minim.minimize(minimizer, algorithm);

    return minim.save();
}

void Drawpull(RooWorkspace* w, TIterator* iter) {
    std::vector<double> pulls;
    std::vector<double> pull_errors;
    std::vector<std::string> names;

    for (TObject* Tobj = iter->Next(); Tobj != 0; Tobj = iter->Next()) {
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(Tobj);
        std::string name = rrv->GetName();
        double val = rrv->getVal();
        double err = rrv->getError();

        std::cout.width(25);
        std::cout << name;
        std::cout.width(15);
        std::cout << val;
        std::cout.width(5);
        std::cout << "+-" << err << std::endl;

        if (name.find("alpha") != std::string::npos) {
            pulls.push_back((val - 0.0) / 1.0);
            pull_errors.push_back(err / 1.0);
            names.push_back(name);
        }
        else if (name.find("gamma_stat") != std::string::npos) {
            RooRealVar* norm = w->var(("nom_" + name).c_str());
            double width = std::pow(norm->getValV(), -0.5);
            pulls.push_back((val - 1.0) / width);
            pull_errors.push_back(err / width);
            names.push_back(name);
        }
        else if ((name.find("gamma") != std::string::npos) && (name.find("uncorr") != std::string::npos)) {
            int sample_index = -1;
            int bin_index = -1;

            if (name.find("CHG") != std::string::npos) sample_index = 0;
            else if (name.find("MIX") != std::string::npos) sample_index = 1;
            else if (name.find("UUBAR") != std::string::npos) sample_index = 2;
            else if (name.find("DDBAR") != std::string::npos) sample_index = 3;
            else if (name.find("SSBAR") != std::string::npos) sample_index = 4;
            else if (name.find("CHARM") != std::string::npos) sample_index = 5;
            else if (name.find("Signal") != std::string::npos) sample_index = 6;

            std::vector<std::string> temp_strings = split(name, '_');
            bin_index = stoi(temp_strings.back()); // from 0

            if (name.find("all") != std::string::npos) {
                double KID_uncertainty = weight_KIDsys[RarityBins * sample_index + bin_index];
                double PID_uncertainty = weight_PIDsys[RarityBins * sample_index + bin_index];
                double BR_uncertainty = 0.0;
                double pi0_uncertainty = weight_pi0sys[RarityBins * sample_index + bin_index];
                double FEI_uncertainty = 0.0;
                if ((name.find("CHG") != std::string::npos) || (name.find("MIX") != std::string::npos)) {
                    BR_uncertainty = weight_BRsys[RarityBins * sample_index + bin_index];
                    FEI_uncertainty = weight_FEIsys[RarityBins * sample_index + bin_index];
                }
                else if (name.find("Signal") != std::string::npos) {
                    BR_uncertainty = weight_BRsys[RarityBins * 2 + bin_index]; // exception for signal BB BR uncorrelated uncertainty!
                    FEI_uncertainty = weight_FEIsys[RarityBins * 2 + bin_index]; // exception for signal FEI uncorrelated uncertainty!
                }

                double total_uncertainty = std::sqrt(KID_uncertainty * KID_uncertainty + PID_uncertainty * PID_uncertainty + BR_uncertainty * BR_uncertainty + pi0_uncertainty * pi0_uncertainty + FEI_uncertainty * FEI_uncertainty);
                pulls.push_back((val - 1.0) / total_uncertainty);
                pull_errors.push_back(err / total_uncertainty);
                names.push_back(name);
            }
        }
    }

    for (unsigned int i = 0; i < names.size(); i++) {
        std::cout << names.at(i) << " " << pulls.at(i) << "+-" << pull_errors.at(i) << std::endl;
    }

    // draw pull
    int size_pull = pulls.size();

    TH1D* pull_ht = new TH1D("pull data hist", "pull of parameters;;", size_pull, 0, size_pull);
    for (int i = 0; i < size_pull; i++) {
        pull_ht->SetBinContent(i + 1, pulls.at(i));
        pull_ht->SetBinError(i + 1, pull_errors.at(i));
    }
    pull_ht->SetLineWidth(2.0);
    pull_ht->SetMarkerColor(1);
    pull_ht->SetMarkerStyle(21);
    pull_ht->SetLineColor(1);
    char** label_name = (char**)malloc(sizeof(char*) * size_pull);
    for (int i = 0; i < size_pull; i++) {
        label_name[i] = (char*)malloc(sizeof(char) * names.at(i).size() + 1);
        memcpy(label_name[i], names.at(i).c_str(), names.at(i).size() + 1);
    }

    TH1D* pull_one_sigma = new TH1D("1sig hist", "1sig;;", size_pull, 0.0, size_pull);
    for (int i = 0; i < size_pull; i++) {
        pull_one_sigma->SetBinContent(i + 1, 0.0);
        pull_one_sigma->SetBinError(i + 1, 1.0);
    }
    pull_one_sigma->SetFillColor(kGreen);
    pull_one_sigma->SetFillStyle(1001);

    TH1D* pull_two_sigma = new TH1D("2sig hist", "2sig;;", size_pull, 0.0, size_pull);
    for (int i = 0; i < size_pull; i++) {
        pull_two_sigma->SetBinContent(i + 1, 0.0);
        pull_two_sigma->SetBinError(i + 1, 2.0);
    }
    pull_two_sigma->SetFillColor(kYellow);
    pull_two_sigma->SetFillStyle(1001);
    for (int i = 0; i < size_pull; i++) {
        pull_two_sigma->GetXaxis()->SetBinLabel(i + 1, names.at(i).c_str());
    }
    pull_two_sigma->SetStats(false);
    pull_two_sigma->GetYaxis()->SetTitle("(#hat{#theta}-#theta)/#Delta#theta");
    pull_two_sigma->GetYaxis()->SetTitleOffset(1.4);
    pull_two_sigma->GetXaxis()->LabelsOption("v");

    TLine* line = new TLine(0.0, 0.0, size_pull, 0.0);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2); line->SetLineWidth(1);

    TCanvas* cpull = new TCanvas("pull_Plot", "pull Plot", 1800, 800); cpull->SetBottomMargin(0.3);
    pull_two_sigma->Draw("E2");
    pull_one_sigma->Draw("E2 same");
    pull_ht->Draw("e1 same");
    line->Draw();

    cpull->SaveAs("param_pull.png");

    delete cpull;
}

void ReadPIDuncorrsysFile(const char* dirname_KID, const char* dirname_PID) {
    FILE* fp;

    fp = fopen(dirname_KID, "r");
    for (int i = 0; i < RarityBins * 7; i++) fscanf(fp, "%lf\n", &weight_KIDsys[i]);
    fclose(fp);
    for (int i = 0; i < RarityBins * 7; i++) weight_KIDsys[i] = std::sqrt(weight_KIDsys[i]);

    fp = fopen(dirname_PID, "r");
    for (int i = 0; i < RarityBins * 7; i++) fscanf(fp, "%lf\n", &weight_PIDsys[i]);
    fclose(fp);
    for (int i = 0; i < RarityBins * 7; i++) weight_PIDsys[i] = std::sqrt(weight_PIDsys[i]);
}

void ReadBRuncorrsysFile(const char* dirname_BR) {
    FILE* fp;

    fp = fopen(dirname_BR, "r");
    for (int i = 0; i < RarityBins * 3; i++) fscanf(fp, "%lf\n", &weight_BRsys[i]);
    fclose(fp);
    for (int i = 0; i < RarityBins * 3; i++) weight_BRsys[i] = std::sqrt(weight_BRsys[i]);

}

void Readpi0uncorrsysFile(const char* dirname_pi0) {
    FILE* fp;

    fp = fopen(dirname_pi0, "r");
    for (int i = 0; i < RarityBins * 7; i++) fscanf(fp, "%lf\n", &weight_pi0sys[i]);
    fclose(fp);
    for (int i = 0; i < RarityBins * 7; i++) weight_pi0sys[i] = std::sqrt(weight_pi0sys[i]);

}

void ReadFEIuncorrsysFile(const char* dirname_FEI) {
    FILE* fp;

    fp = fopen(dirname_FEI, "r");
    for (int i = 0; i < RarityBins * 3; i++) fscanf(fp, "%lf\n", &weight_FEIsys[i]);
    fclose(fp);
    for (int i = 0; i < RarityBins * 3; i++) weight_FEIsys[i] = std::sqrt(weight_FEIsys[i]);

}

double GetNumEvts(RooWorkspace* w, const char* sample_type) {

    double Nevt = 0.0;

    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();
    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x_val = w->var("obs_x_channel");

    int index = -1;
    if (strcmp(sample_type, "Signal") == 0) index = 0;
    else if (strcmp(sample_type, "CHG") == 0) index = 1;
    else if (strcmp(sample_type, "MIX") == 0) index = 2;
    else if (strcmp(sample_type, "UUBAR") == 0) index = 3;
    else if (strcmp(sample_type, "DDBAR") == 0) index = 4;
    else if (strcmp(sample_type, "SSBAR") == 0) index = 5;
    else if (strcmp(sample_type, "CHARM") == 0) index = 6;
    else {
        printf("[ERROR] unexpected sample type!\n");
        exit(1);
    }

    /* ================================ cal Nexpected ================================*/
    RooAbsBinning const& binning = x_val->getBinning();
    const double oldVal = x_val->getVal();

    for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
        double binCenter = binning.binCenter(iBin);
        double binWidth = binning.binWidth(iBin);

        *x_val = binCenter; // set x value

        RooAbsReal* temp_func = w->function(Sample_names.at(index).c_str());
        Nevt = Nevt + temp_func->getValV();
        if (temp_func->getValV() < 0) {
            printf("[ERROR] negative count!\n");
            exit(1);
        }

    }

    *x_val = oldVal;

    return Nevt;

}

int Check_param() {

    ::ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit"); // default: Minuit Migrad
    ::ROOT::Math::MinimizerOptions::SetDefaultStrategy(1); // default 1

    ReadPIDuncorrsysFile("./KID_cov_remain_truncated.txt", "./PID_cov_remain_truncated.txt");
    ReadBRuncorrsysFile("./BR_cov_remain_truncated.txt");
    Readpi0uncorrsysFile("./pi0_cov_remain_truncated.txt");
    ReadFEIuncorrsysFile("./FEI_cov_remain_truncated.txt");

    const char* fname = "./PDFandDATA_workspace.root";

    TFile* f = TFile::Open(fname);

    RooWorkspace* w = (RooWorkspace*)f->Get("combined");

    w->Print();
    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    OPTIONS* options = (OPTIONS*)malloc(sizeof(OPTIONS));
    Initialize_options(options, "none");
    FixParameters(w, options);

    // test
    //RooRealVar* alpha = w->var("nom_gamma_stat_channel_bin_0");
    //printf("%lf", alpha->getValV());

    // Lets tell roofit the right names for our histogram variables //
    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x = (RooRealVar*)obs->find("obs_x_channel");
    x->SetTitle("FBDT output");
    x->setUnit("");

    // get Category and data
    RooCategory* idx = (RooCategory*)obs->find("channelCat");
    //RooAbsData* data = (RooAbsData*)w->data("obsData");
    RooDataSet* data = (RooDataSet*)w->data("asimovData");

    // fit
    double eps = 0.001;
    RooAbsReal* nll;
    RooFitResult* fitres = MyMinimizeNLL(w, data, &nll, eps);

    RooArgSet fitargs = fitres->floatParsFinal();
    TIterator* iter(fitargs.createIterator());

    for (TObject* a = iter->Next(); a != 0; a = iter->Next()) {
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(a);
        std::string name = rrv->GetName();
        double val = rrv->getVal();
        double err = rrv->getError();
        double HIerr = rrv->getAsymErrorHi();
        double LOerr = rrv->getAsymErrorLo();

        if (name == "mu") {
            printf("fit result mu = %lf +- %lf\n", val, err);
            printf("MINOS error: %lf %lf\n", HIerr, LOerr);
        }

    }

    // draw pull
    Drawpull(w, iter);

    // define frame
    RooPlot* x_frame = x->frame(Title("fit result"));

    // get expected num of evts for PDFs
    double Signal_Nevts = GetNumEvts(w, "Signal");
    double CHG_Nevts = GetNumEvts(w, "CHG");
    double MIX_Nevts = GetNumEvts(w, "MIX");
    double UUBAR_Nevts = GetNumEvts(w, "UUBAR");
    double DDBAR_Nevts = GetNumEvts(w, "DDBAR");
    double SSBAR_Nevts = GetNumEvts(w, "SSBAR");
    double CHARM_Nevts = GetNumEvts(w, "CHARM");

    // draw
    model->plotOn(x_frame, Name("signal_name"), Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kRed - 6), LineWidth(0), Components("L_x_*Signal*_ShapeSys,L_x_*CHG*_ShapeSys,L_x_*MIX*_ShapeSys,L_x_*UUBAR*_ShapeSys,L_x_*DDBAR*_ShapeSys,L_x_*SSBAR*_ShapeSys,L_x_*CHARM*_ShapeSys"), Normalization(Signal_Nevts + CHG_Nevts + MIX_Nevts + UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent));
    model->plotOn(x_frame, Name("CHG_name"), Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kBlue - 6), LineWidth(0), Components("L_x_*CHG*_ShapeSys,L_x_*MIX*_ShapeSys,L_x_*UUBAR*_ShapeSys,L_x_*DDBAR*_ShapeSys,L_x_*SSBAR*_ShapeSys,L_x_*CHARM*_ShapeSys"), Normalization(CHG_Nevts + MIX_Nevts + UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent));
    model->plotOn(x_frame, Name("MIX_name"), Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kCyan - 6), LineWidth(0), Components("L_x_*MIX*_ShapeSys,L_x_*UUBAR*_ShapeSys,L_x_*DDBAR*_ShapeSys,L_x_*SSBAR*_ShapeSys,L_x_*CHARM*_ShapeSys"), Normalization(MIX_Nevts + UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent));
    model->plotOn(x_frame, Name("UUBAR_name"), Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kOrange - 6), LineWidth(0), Components("L_x_*UUBAR*_ShapeSys,L_x_*DDBAR*_ShapeSys,L_x_*SSBAR*_ShapeSys,L_x_*CHARM*_ShapeSys"), Normalization(UUBAR_Nevts + DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent));
    model->plotOn(x_frame, Name("DDBAR_name"), Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kViolet - 6), LineWidth(0), Components("L_x_*DDBAR*_ShapeSys,L_x_*SSBAR*_ShapeSys,L_x_*CHARM*_ShapeSys"), Normalization(DDBAR_Nevts + SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent));
    model->plotOn(x_frame, Name("SSBAR_name"), Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kGreen - 6), LineWidth(0), Components("L_x_*SSBAR*_ShapeSys,L_x_*CHARM*_ShapeSys"), Normalization(SSBAR_Nevts + CHARM_Nevts, RooAbsReal::ScaleType::NumEvent));
    model->plotOn(x_frame, Name("CHARM_name"), Slice(*idx), ProjWData(*idx, *data), DrawOption("F"), FillColor(kYellow - 6), LineWidth(0), Components("L_x_*CHARM*_ShapeSys"), Normalization(CHARM_Nevts, RooAbsReal::ScaleType::NumEvent));
    data->plotOn(x_frame, Name("data_name"), DataError(RooAbsData::Poisson), Cut("channelCat==0"), RooFit::MarkerSize(1.0), RooFit::LineWidth(1.0), DrawOption("ZP"));
    //model->paramOn(x_frame);
    TCanvas* cdata = new TCanvas("sPlot", "sPlot demo", 700, 700);
    x_frame->Draw();

    TLegend* leg = new TLegend(0.7, 0.6, 0.9, 0.9);
    leg->SetFillStyle(0);
    leg->SetLineWidth(0);
    leg->AddEntry("data_name", "Data", "LP");
    leg->AddEntry("CHG_name", "charged", "F");
    leg->AddEntry("MIX_name", "mixed", "F");
    leg->AddEntry("UUBAR_name", "u#bar{u}", "F");
    leg->AddEntry("DDBAR_name", "d#bar{d}", "F");
    leg->AddEntry("SSBAR_name", "s#bar{s}", "F");
    leg->AddEntry("CHARM_name", "c#bar{c}", "F");
    leg->AddEntry("signal_name", "signal", "F");
    leg->Draw();

    cdata->SaveAs("FitResult.png");
    delete cdata;

    // draw profile likelihood
    RooRealVar* mu = w->var("mu");
    mu->setRange(mu->getValV() - 4, mu->getValV() + 4);
    RooPlot* mu_frame = mu->frame();
    RooAbsReal* pll = nll->createProfile(*mu);
    pll->plotOn(mu_frame);

    TCanvas* cmu = new TCanvas("pllPlot", "pllPlot", 700, 700);
    mu_frame->Draw();
    cmu->SaveAs("profile_likelihood_mu.png");
    delete cmu;

    return 0;
}
