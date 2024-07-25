#ifndef FITTER_H
#define FITTER_H

# include <vector>
# include <deque>
# include <string>
# include "RooFitResult.h"
# include "template.h"
# include "correctors.h"
# include <cstring>
# include "THStack.h"
# include "TStyle.h"
# include "TLine.h"
# include "TColor.h"
# include "TPaveText.h"
#include "systematic.h"

using namespace RooFit;
using namespace RooStats;
using namespace HistFactory;

Corrector_Fragmentation corrector_Fragmentation;

std::vector<std::string> scaleFactors_pdf_names = {
    "Signal_MXs1_nominal_MXs1_channel_MXs1_scaleFactors",
    "Signal_MXs2_nominal_MXs1_channel_MXs1_scaleFactors",
    "Signal_MXs3_nominal_MXs1_channel_MXs1_scaleFactors",
    "Signal_MXs1_nominal_MXs2_channel_MXs2_scaleFactors",
    "Signal_MXs2_nominal_MXs2_channel_MXs2_scaleFactors",
    "Signal_MXs3_nominal_MXs2_channel_MXs2_scaleFactors",
    "Signal_MXs1_nominal_MXs3_channel_MXs3_scaleFactors",
    "Signal_MXs2_nominal_MXs3_channel_MXs3_scaleFactors",
    "Signal_MXs3_nominal_MXs3_channel_MXs3_scaleFactors",
    "CHG_nominal_MXs1_channel_MXs1_scaleFactors",
    "CHG_nominal_MXs2_channel_MXs2_scaleFactors",
    "CHG_nominal_MXs3_channel_MXs3_scaleFactors",
    "MIX_nominal_MXs1_channel_MXs1_scaleFactors",
    "MIX_nominal_MXs2_channel_MXs2_scaleFactors",
    "MIX_nominal_MXs3_channel_MXs3_scaleFactors",
    "UUBAR_nominal_MXs1_channel_MXs1_scaleFactors",
    "UUBAR_nominal_MXs2_channel_MXs2_scaleFactors",
    "UUBAR_nominal_MXs3_channel_MXs3_scaleFactors",
    "DDBAR_nominal_MXs1_channel_MXs1_scaleFactors",
    "DDBAR_nominal_MXs2_channel_MXs2_scaleFactors",
    "DDBAR_nominal_MXs3_channel_MXs3_scaleFactors",
    "SSBAR_nominal_MXs1_channel_MXs1_scaleFactors",
    "SSBAR_nominal_MXs2_channel_MXs2_scaleFactors",
    "SSBAR_nominal_MXs3_channel_MXs3_scaleFactors",
    "CHARM_nominal_MXs1_channel_MXs1_scaleFactors",
    "CHARM_nominal_MXs2_channel_MXs2_scaleFactors",
    "CHARM_nominal_MXs3_channel_MXs3_scaleFactors"
};

std::vector<std::string> shapes_pdf_names = {
    "Signal_MXs1_nominal_MXs1_channel_MXs1_shapes",
    "Signal_MXs2_nominal_MXs1_channel_MXs1_shapes",
    "Signal_MXs3_nominal_MXs1_channel_MXs1_shapes",
    "Signal_MXs1_nominal_MXs2_channel_MXs2_shapes",
    "Signal_MXs2_nominal_MXs2_channel_MXs2_shapes",
    "Signal_MXs3_nominal_MXs2_channel_MXs2_shapes",
    "Signal_MXs1_nominal_MXs3_channel_MXs3_shapes",
    "Signal_MXs2_nominal_MXs3_channel_MXs3_shapes",
    "Signal_MXs3_nominal_MXs3_channel_MXs3_shapes",
    "CHG_nominal_MXs1_channel_MXs1_shapes",
    "CHG_nominal_MXs2_channel_MXs2_shapes",
    "CHG_nominal_MXs3_channel_MXs3_shapes",
    "MIX_nominal_MXs1_channel_MXs1_shapes",
    "MIX_nominal_MXs2_channel_MXs2_shapes",
    "MIX_nominal_MXs3_channel_MXs3_shapes",
    "UUBAR_nominal_MXs1_channel_MXs1_shapes",
    "UUBAR_nominal_MXs2_channel_MXs2_shapes",
    "UUBAR_nominal_MXs3_channel_MXs3_shapes",
    "DDBAR_nominal_MXs1_channel_MXs1_shapes",
    "DDBAR_nominal_MXs2_channel_MXs2_shapes",
    "DDBAR_nominal_MXs3_channel_MXs3_shapes",
    "SSBAR_nominal_MXs1_channel_MXs1_shapes",
    "SSBAR_nominal_MXs2_channel_MXs2_shapes",
    "SSBAR_nominal_MXs3_channel_MXs3_shapes",
    "CHARM_nominal_MXs1_channel_MXs1_shapes",
    "CHARM_nominal_MXs2_channel_MXs2_shapes",
    "CHARM_nominal_MXs3_channel_MXs3_shapes"
};

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
    bool BRBtoXKLKL;
    bool EffECLKL;
    bool NEWFEICAL;
    bool BRXnn;
    bool BRDKL0;
    bool dataMC;
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
    double weight_sys[RarityBins * 5] = { 0.0 }; // eigen vector

    FILE* fp;
    fp = fopen(dirname, "r");
    while (true) {
        if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
        for (int i = 0; i < RarityBins * 5; i++) {
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
    double weight_sys[RarityBins * 9] = { 0.0 }; // eigen vector

    FILE* fp;
    fp = fopen(dirname, "r");
    while (true) {
        if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
        for (int i = 0; i < RarityBins * 9; i++) {
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
    double weight_sys[RarityBins * 9] = { 0.0 }; // eigen vector

    FILE* fp;
    fp = fopen(dirname, "r");
    while (true) {
        if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
        for (int i = 0; i < RarityBins * 9; i++) {
            if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
        }
        Nentry++;
    }
    fclose(fp);

    return Nentry;
}

int ReadMultiplicityInfo(const char* dirname) {
    int Nentry = 0; // number of eigen values/vectors
    double eigen_value = 0; // eigen value
    double weight_sys[RarityBins * 9] = { 0.0 }; // eigen vector

    FILE* fp;
    fp = fopen(dirname, "r");
    while (true) {
        if (fscanf(fp, "%lf\n", &eigen_value) == EOF) break;
        for (int i = 0; i < RarityBins * 9; i++) {
            if (fscanf(fp, "%lf\n", &weight_sys[i]) == EOF) break;
        }
        Nentry++;
    }
    fclose(fp);

    return Nentry;
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
    options_->BRBtoXKLKL = false;
    options_->EffECLKL = false;
    options_->NEWFEICAL = false;
    options_->BRXnn = false;
    options_->BRDKL0 = false;
    options_->dataMC = false;
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
        options_->BRBtoXKLKL = true;
        options_->EffECLKL = true;
        options_->NEWFEICAL = true;
        options_->BRXnn = true;
        options_->BRDKL0 = true;
        options_->dataMC = true;
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
    else if (std::string(tested_param) == std::string("BRBtoXKLKL")) options_->BRBtoXKLKL = true;
    else if (std::string(tested_param) == std::string("EffECLKL")) options_->EffECLKL = true;
    else if (std::string(tested_param) == std::string("NEWFEICAL")) options_->NEWFEICAL = true;
    else if (std::string(tested_param) == std::string("BRXnn")) options_->BRXnn = true;
    else if (std::string(tested_param) == std::string("BtoDtoXKL")) options_->BRDKL0 = true;
    else if (std::string(tested_param) == std::string("dataMC")) options_->dataMC = true;
    else if (std::string(tested_param) == std::string("uncorrelated")) options_->uncorrelated = true;
    else {
        printf("inappropriate parameter name\n");
        exit(1);
    }

    // read entry for nuisance parameters
    options_->NEntryFEI = ReadNFEIEigenVector("./FEI_selected.txt");
    options_->NEntryKID = ReadNPIDEigenVector("./KID_selected.txt");
    options_->NEntryPID = ReadNPIDEigenVector("./PID_selected.txt");
    options_->NEntrypi0 = ReadNpi0EigenVector("./pi0_selected.txt");
    options_->NEntryMultiplicity = ReadMultiplicityInfo("./multiplicity_selected.txt");
}

void FixParameters(RooWorkspace* w, OPTIONS* options_) {

    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    w->loadSnapshot("NominalParamValues");

    RooRealVar* x_val_MXs1 = w->var("obs_x_channel_MXs1");
    RooRealVar* x_val_MXs2 = w->var("obs_x_channel_MXs2");
    RooRealVar* x_val_MXs3 = w->var("obs_x_channel_MXs3");
    std::unique_ptr<RooArgSet> params{model->getParameters(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3))};

    // track
    if (options_->track) w->var("alpha_track_eff_uncer")->setConstant(options_->track);

    // PID
    if (options_->PID) for (int i = 0; i < options_->NEntryPID; i++) w->var(("alpha_PID" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->PID);

    // KID
    if (options_->KID) for (int i = 0; i < options_->NEntryKID; i++) w->var(("alpha_KID" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->KID);

    // KS0
    if (options_->KS0) w->var("alpha_KS0_reco_uncer")->setConstant(options_->KS0);

    // pi0
    if (options_->pi0) for (int i = 0; i < options_->NEntrypi0; i++) w->var(("alpha_pi0" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->pi0);

    // FBDT efficiency
    if (options_->FBDT) {
        w->var("alpha_FBDT_efficiency_uncer_MXs1")->setConstant(options_->FBDT);
        w->var("alpha_FBDT_efficiency_uncer_MXs2")->setConstant(options_->FBDT);
        w->var("alpha_FBDT_efficiency_uncer_MXs3")->setConstant(options_->FBDT);
    }

    // FEI
    if (options_->FEI) for (int i = 0; i < options_->NEntryFEI; i++) w->var(("alpha_FEI" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->FEI);

    // qqbar normalization
    if (options_->qqbar) {
        w->var("alpha_qq_CAL_UUBAR_uncer")->setConstant(options_->qqbar);
        w->var("alpha_qq_CAL_DDBAR_uncer")->setConstant(options_->qqbar);
        w->var("alpha_qq_CAL_SSBAR_uncer")->setConstant(options_->qqbar);
        w->var("alpha_qq_CAL_CHARM_uncer")->setConstant(options_->qqbar);
    }

    // photon multiplicity correction
    if (options_->multiplicity) for (int i = 0; i < options_->NEntryMultiplicity; i++) w->var(("alpha_multiplicity" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->multiplicity);

    // B->K form factor
    if (options_->Kff) {
        w->var("alpha_Kff1_uncer")->setConstant(options_->Kff);
        w->var("alpha_Kff2_uncer")->setConstant(options_->Kff);
        w->var("alpha_Kff3_uncer")->setConstant(options_->Kff);
    }

    // B->Kstar form factor
    if (options_->Kstarff) {
        w->var("alpha_Kstarff1_uncer")->setConstant(options_->Kstarff);
        w->var("alpha_Kstarff2_uncer")->setConstant(options_->Kstarff);
        w->var("alpha_Kstarff3_uncer")->setConstant(options_->Kstarff);
        w->var("alpha_Kstarff4_uncer")->setConstant(options_->Kstarff);
        w->var("alpha_Kstarff5_uncer")->setConstant(options_->Kstarff);
        w->var("alpha_Kstarff6_uncer")->setConstant(options_->Kstarff);
        w->var("alpha_Kstarff7_uncer")->setConstant(options_->Kstarff);
        w->var("alpha_Kstarff8_uncer")->setConstant(options_->Kstarff);
        w->var("alpha_Kstarff9_uncer")->setConstant(options_->Kstarff);
    }

    // fermi motion moment
    if (options_->pf) w->var("alpha_pf_uncer")->setConstant(options_->pf);

    // K*-Xs transition point
    if (options_->Transition) w->var("alpha_transition_uncer")->setConstant(options_->Transition);

    // b-quark mass
    if (options_->mb) w->var("alpha_mb_uncer")->setConstant(options_->mb);

    // relative fraction
    if (options_->fraction) {
        w->var("alpha_Kfrac_uncer")->setConstant(options_->fraction);
        w->var("alpha_Kstarfrac_uncer")->setConstant(options_->fraction);
    }

    // MC statistics
    if (options_->MCstat) {
        for (int i = 0; i < RarityBins_MX1; i++) w->var(("gamma_stat_channel_MXs1_bin_" + std::to_string(i)).c_str())->setConstant(options_->MCstat);
        for (int i = 0; i < RarityBins_MX2; i++) w->var(("gamma_stat_channel_MXs2_bin_" + std::to_string(i)).c_str())->setConstant(options_->MCstat);
        for (int i = 0; i < RarityBins_MX3; i++) w->var(("gamma_stat_channel_MXs3_bin_" + std::to_string(i)).c_str())->setConstant(options_->MCstat);
    }

    // Fragmentation
    if (options_->Fragmentation) {
        for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
            for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
                int temp_index = MxsBin * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma) + Category;
                std::string uncertainty_name = "Xs_fragmentation" + std::to_string(temp_index) + "_uncer";
                w->var(uncertainty_name.c_str())->setConstant(options_->Fragmentation);
            }
        }
    }

    // mKstar
    //if (options_->mKstar) w->var("alpha_mKstar_uncer")->setConstant(options_->mKstar);

    // background normalization
    if (options_->BKGNorm) {
        w->var("alpha_mu_CHG_MXs1")->setConstant(options_->BKGNorm);
        w->var("alpha_mu_CHG_MXs2")->setConstant(options_->BKGNorm);
        w->var("alpha_mu_CHG_MXs3")->setConstant(options_->BKGNorm);
        w->var("alpha_mu_MIX_MXs1")->setConstant(options_->BKGNorm);
        w->var("alpha_mu_MIX_MXs2")->setConstant(options_->BKGNorm);
        w->var("alpha_mu_MIX_MXs3")->setConstant(options_->BKGNorm);
        w->var("alpha_mu_UUBAR_MXs1")->setConstant(options_->BKGNorm);
        w->var("alpha_mu_UUBAR_MXs2")->setConstant(options_->BKGNorm);
        w->var("alpha_mu_UUBAR_MXs3")->setConstant(options_->BKGNorm);
        w->var("alpha_mu_DDBAR_MXs1")->setConstant(options_->BKGNorm);
        w->var("alpha_mu_DDBAR_MXs2")->setConstant(options_->BKGNorm);
        w->var("alpha_mu_DDBAR_MXs3")->setConstant(options_->BKGNorm);
        w->var("alpha_mu_SSBAR_MXs1")->setConstant(options_->BKGNorm);
        w->var("alpha_mu_SSBAR_MXs2")->setConstant(options_->BKGNorm);
        w->var("alpha_mu_SSBAR_MXs3")->setConstant(options_->BKGNorm);
        w->var("alpha_mu_CHARM_MXs1")->setConstant(options_->BKGNorm);
        w->var("alpha_mu_CHARM_MXs2")->setConstant(options_->BKGNorm);
        w->var("alpha_mu_CHARM_MXs3")->setConstant(options_->BKGNorm);
    }

    // BDTc
    if (options_->BDTc) w->var("alpha_BDTc_shape_qq")->setConstant(options_->BDTc);

    // BB counting
    if (options_->BBcounting) w->var("alpha_BB_counting_uncer")->setConstant(options_->BBcounting);

    // BB BR
    if (options_->BBBR) for (int i = 0; i < NBRdmID(); i++) w->var(("alpha_BBBR" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->BBBR);

    // B->K KL KL BR
    if (options_->BRBtoXKLKL) w->var("alpha_BRBtoXKLKL_uncer")->setConstant(options_->BRBtoXKLKL);

    // efficiency for ECL cluster from KL0
    if (options_->EffECLKL) w->var("alpha_EffECLKL_uncer")->setConstant(options_->EffECLKL);

    // New FEI CAL
    if (options_->NEWFEICAL) w->var("alpha_NEWFEICAL_uncer")->setConstant(options_->NEWFEICAL);

    // B->Xnn BR
    if (options_->BRXnn) w->var("alpha_Xnn_BR_uncer")->setConstant(options_->BRXnn);

    // B->(D->X KL0) + anything
    if (options_->BRDKL0) w->var("alpha_BtoDtoXKL_uncer")->setConstant(options_->BRDKL0);

    // data MC discrepancy
    if (options_->dataMC) {
        /*
        for (int i = 0; i < RarityBins_MX1; i++) {
            w->var(("gamma_CHG_additional_uncorr_uncerMXs1_bin_" + std::to_string(i)).c_str())->setConstant(options_->dataMC);
            w->var(("gamma_MIX_additional_uncorr_uncerMXs1_bin_" + std::to_string(i)).c_str())->setConstant(options_->dataMC);
            w->var(("gamma_UUBAR_additional_uncorr_uncerMXs1_bin_" + std::to_string(i)).c_str())->setConstant(options_->dataMC);
            w->var(("gamma_DDBAR_additional_uncorr_uncerMXs1_bin_" + std::to_string(i)).c_str())->setConstant(options_->dataMC);
            w->var(("gamma_SSBAR_additional_uncorr_uncerMXs1_bin_" + std::to_string(i)).c_str())->setConstant(options_->dataMC);
            w->var(("gamma_CHARM_additional_uncorr_uncerMXs1_bin_" + std::to_string(i)).c_str())->setConstant(options_->dataMC);
        }
        for (int i = 0; i < RarityBins_MX2; i++) {
            w->var(("gamma_CHG_additional_uncorr_uncerMXs2_bin_" + std::to_string(i)).c_str())->setConstant(options_->dataMC);
            w->var(("gamma_MIX_additional_uncorr_uncerMXs2_bin_" + std::to_string(i)).c_str())->setConstant(options_->dataMC);
            w->var(("gamma_UUBAR_additional_uncorr_uncerMXs2_bin_" + std::to_string(i)).c_str())->setConstant(options_->dataMC);
            w->var(("gamma_DDBAR_additional_uncorr_uncerMXs2_bin_" + std::to_string(i)).c_str())->setConstant(options_->dataMC);
            w->var(("gamma_SSBAR_additional_uncorr_uncerMXs2_bin_" + std::to_string(i)).c_str())->setConstant(options_->dataMC);
            w->var(("gamma_CHARM_additional_uncorr_uncerMXs2_bin_" + std::to_string(i)).c_str())->setConstant(options_->dataMC);
        }
        for (int i = 0; i < RarityBins_MX3; i++) {
            w->var(("gamma_CHG_additional_uncorr_uncerMXs3_bin_" + std::to_string(i)).c_str())->setConstant(options_->dataMC);
            w->var(("gamma_MIX_additional_uncorr_uncerMXs3_bin_" + std::to_string(i)).c_str())->setConstant(options_->dataMC);
            w->var(("gamma_UUBAR_additional_uncorr_uncerMXs3_bin_" + std::to_string(i)).c_str())->setConstant(options_->dataMC);
            w->var(("gamma_DDBAR_additional_uncorr_uncerMXs3_bin_" + std::to_string(i)).c_str())->setConstant(options_->dataMC);
            w->var(("gamma_SSBAR_additional_uncorr_uncerMXs3_bin_" + std::to_string(i)).c_str())->setConstant(options_->dataMC);
            w->var(("gamma_CHARM_additional_uncorr_uncerMXs3_bin_" + std::to_string(i)).c_str())->setConstant(options_->dataMC);
        }
        */
    }

    // save snapshot
    w->saveSnapshot("ParamValues", *params, true);
}

double GetNumEvts(RooWorkspace* w, const char* sample_type) {

    double Nevt = 0.0;

    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();
    RooArgSet* obs = (RooArgSet*)mc->GetObservables();
    RooRealVar* x_val;
    if (std::strstr(sample_type, "MX1") != nullptr) x_val = w->var("obs_x_channel_MXs1");
    else if (std::strstr(sample_type, "MX2") != nullptr) x_val = w->var("obs_x_channel_MXs2");
    else if (std::strstr(sample_type, "MX3") != nullptr) x_val = w->var("obs_x_channel_MXs3");
    else {
        printf("[ERROR] unexpected sample type!\n");
        exit(1);
    }

    std::vector<int> indices;
    if (strcmp(sample_type, "Signal_MX1") == 0) {
        indices.push_back(0);
        indices.push_back(1);
        indices.push_back(2);
    }
    else if (strcmp(sample_type, "Signal_MX2") == 0) {
        indices.push_back(3);
        indices.push_back(4);
        indices.push_back(5);
    }
    else if (strcmp(sample_type, "Signal_MX3") == 0) {
        indices.push_back(6);
        indices.push_back(7);
        indices.push_back(8);
    }
    else if (strcmp(sample_type, "CHG_MX1") == 0) indices.push_back(9);
    else if (strcmp(sample_type, "CHG_MX2") == 0) indices.push_back(10);
    else if (strcmp(sample_type, "CHG_MX3") == 0) indices.push_back(11);
    else if (strcmp(sample_type, "MIX_MX1") == 0) indices.push_back(12);
    else if (strcmp(sample_type, "MIX_MX2") == 0) indices.push_back(13);
    else if (strcmp(sample_type, "MIX_MX3") == 0) indices.push_back(14);
    else if (strcmp(sample_type, "UUBAR_MX1") == 0) indices.push_back(15);
    else if (strcmp(sample_type, "UUBAR_MX2") == 0) indices.push_back(16);
    else if (strcmp(sample_type, "UUBAR_MX3") == 0) indices.push_back(17);
    else if (strcmp(sample_type, "DDBAR_MX1") == 0) indices.push_back(18);
    else if (strcmp(sample_type, "DDBAR_MX2") == 0) indices.push_back(19);
    else if (strcmp(sample_type, "DDBAR_MX3") == 0) indices.push_back(20);
    else if (strcmp(sample_type, "SSBAR_MX1") == 0) indices.push_back(21);
    else if (strcmp(sample_type, "SSBAR_MX2") == 0) indices.push_back(22);
    else if (strcmp(sample_type, "SSBAR_MX3") == 0) indices.push_back(23);
    else if (strcmp(sample_type, "CHARM_MX1") == 0) indices.push_back(24);
    else if (strcmp(sample_type, "CHARM_MX2") == 0) indices.push_back(25);
    else if (strcmp(sample_type, "CHARM_MX3") == 0) indices.push_back(26);
    else {
        printf("[ERROR] unexpected sample type!\n");
        exit(1);
    }


    /* ================================ cal Nexpected ================================*/
    for (int i = 0; i < indices.size(); i++) {
        int index = indices.at(i);

        RooAbsBinning const& binning = x_val->getBinning();
        const double oldVal = x_val->getVal();

        for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
            double binCenter = binning.binCenter(iBin);
            double binWidth = binning.binWidth(iBin);

            *x_val = binCenter; // set x value

            RooAbsReal* temp_func_scaleFactors = w->function(scaleFactors_pdf_names.at(index).c_str());
            RooAbsReal* temp_func_shapes = w->function(shapes_pdf_names.at(index).c_str());
            Nevt = Nevt + (temp_func_scaleFactors->getValV() * temp_func_shapes->getValV());
            //if (temp_func->getValV() < 0) {
            //    printf("[ERROR] negative count!\n");
            //    exit(1);
            //}

        }

        *x_val = oldVal;
    }

    return Nevt;

}

RooFitResult* MyMinimizeNLL(RooWorkspace* w, RooDataSet* data, RooAbsReal** nll, double tolerance = -1.0, bool Minos = true) {
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
    if (Minos) {
        minim.minos(RooArgSet(*w->var("mu_MXs1")));
        minim.minos(RooArgSet(*w->var("mu_MXs2")));
        minim.minos(RooArgSet(*w->var("mu_MXs3")));
    }

    // fit!
    int status;
    status = minim.minimize(minimizer, algorithm);

    return minim.save();
}

RooFitResult* MyMinimizeNLLReuse(RooWorkspace* w, RooDataSet* data, RooAbsReal** nll, double tolerance = -1.0, bool Minos = true) {
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
    if (Minos) {
        minim.minos(RooArgSet(*w->var("mu_MXs1")));
        minim.minos(RooArgSet(*w->var("mu_MXs2")));
        minim.minos(RooArgSet(*w->var("mu_MXs3")));
    }

    // fit!
    int status;
    status = minim.minimize(minimizer, algorithm);

    return minim.save();
}

RooFitResult* MyMinimizeNLLWithAsymError(RooWorkspace* w, RooDataSet* data, RooAbsReal** nll, double tolerance = -1.0, bool Minos = true) {
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
    if (Minos) {
        minim.minos(RooArgSet(*w->var("mu_MXs1")));
        minim.minos(RooArgSet(*w->var("mu_MXs2")));
        minim.minos(RooArgSet(*w->var("mu_MXs3")));
    }

    // fit!
    int status;
    status = minim.minimize(minimizer, algorithm);

    // save snapshot at global minimum
    RooRealVar* x_val_MXs1 = w->var("obs_x_channel_MXs1");
    RooRealVar* x_val_MXs2 = w->var("obs_x_channel_MXs2");
    RooRealVar* x_val_MXs3 = w->var("obs_x_channel_MXs3");
    std::unique_ptr<RooArgSet> params{model->getParameters(RooArgSet(*x_val_MXs1, *x_val_MXs2, *x_val_MXs3))};
    w->saveSnapshot("GlobalMinimumParamValues", *params, true);
    w->saveSnapshot("CurrentMinimumParamValues_MXs1", *params, true);
    w->saveSnapshot("CurrentMinimumParamValues_MXs2", *params, true);
    w->saveSnapshot("CurrentMinimumParamValues_MXs3", *params, true);

    // get global minimum
    double Global_mu_MXs1 = -100;
    double Global_mu_MXs2 = -100;
    double Global_mu_MXs3 = -100;
    double mu_MXs1_HI_error = -100;
    double mu_MXs2_HI_error = -100;
    double mu_MXs3_HI_error = -100;
    double mu_MXs1_LO_error = -100;
    double mu_MXs2_LO_error = -100;
    double mu_MXs3_LO_error = -100;
    double Global_MinusLogLikelihood = (*nll)->getVal();
    RooFitResult* fitres = minim.save();
    RooArgSet fitargs = fitres->floatParsFinal();
    TIterator* iter(fitargs.createIterator());
    for (TObject* a = iter->Next(); a != 0; a = iter->Next()) {
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(a);
        std::string name = rrv->GetName();
        if (name == std::string("mu_MXs1")) {
            Global_mu_MXs1 = rrv->getVal();
            mu_MXs1_HI_error = rrv->getAsymErrorHi();
            mu_MXs1_LO_error = rrv->getAsymErrorLo();
        }
        else if (name == std::string("mu_MXs2")) {
            Global_mu_MXs2 = rrv->getVal();
            mu_MXs2_HI_error = rrv->getAsymErrorHi();
            mu_MXs2_LO_error = rrv->getAsymErrorLo();
        }
        else if (name == std::string("mu_MXs3")) {
            Global_mu_MXs3 = rrv->getVal();
            mu_MXs3_HI_error = rrv->getAsymErrorHi();
            mu_MXs3_LO_error = rrv->getAsymErrorLo();
        }
    }

    // vector to save mu
    std::deque<double> mu_MXs1; mu_MXs1.push_back(Global_mu_MXs1);
    std::deque<double> mu_MXs2; mu_MXs2.push_back(Global_mu_MXs2);
    std::deque<double> mu_MXs3; mu_MXs3.push_back(Global_mu_MXs3);

    // vector to save -log(L)
    std::deque<double> MinusLogLikelihood_MXs1; MinusLogLikelihood_MXs1.push_back(Global_MinusLogLikelihood);
    std::deque<double> MinusLogLikelihood_MXs2; MinusLogLikelihood_MXs2.push_back(Global_MinusLogLikelihood);
    std::deque<double> MinusLogLikelihood_MXs3; MinusLogLikelihood_MXs3.push_back(Global_MinusLogLikelihood);

    double delta = 0.1;

    // re-define minimizer to remove MINOS option
    RooMinimizer MyMinim(*(*nll));
    MyMinim.setStrategy(fStrategy);
    MyMinim.setEvalErrorWall(config.useEvalErrorWall);
    MyMinim.setEps(fTolerance);
    MyMinim.setPrintLevel(level - 1);
    // this causes a memory leak
    MyMinim.optimizeConst(2);
    MyMinim.migrad();

    // get asym error for mu_MXs1
    for (int i = 1; i < 15; i++) {
        w->loadSnapshot("CurrentMinimumParamValues_MXs1");
        w->var("mu_MXs1")->setVal(Global_mu_MXs1 + i * delta * mu_MXs1_HI_error);
        w->var("mu_MXs1")->setConstant(true);

        // fit with fixed mu
        MyMinim.minimize(minimizer, algorithm);

        // save snapshot
        w->saveSnapshot("CurrentMinimumParamValues_MXs1", *params, true);

        // save result
        mu_MXs1.push_back(Global_mu_MXs1 + i * delta * mu_MXs1_HI_error);
        MinusLogLikelihood_MXs1.push_back((*nll)->getVal());
    }
    w->loadSnapshot("GlobalMinimumParamValues");
    w->saveSnapshot("CurrentMinimumParamValues_MXs1", *params, true);

    for (int i = 1; i < 15; i++) {
        w->loadSnapshot("CurrentMinimumParamValues_MXs1");
        w->var("mu_MXs1")->setVal(Global_mu_MXs1 + i * delta * mu_MXs1_LO_error);
        w->var("mu_MXs1")->setConstant(true);

        // fit with fixed mu
        MyMinim.minimize(minimizer, algorithm);

        // save snapshot
        w->saveSnapshot("CurrentMinimumParamValues_MXs1", *params, true);

        // save result
        mu_MXs1.push_front(Global_mu_MXs1 + i * delta * mu_MXs1_LO_error);
        MinusLogLikelihood_MXs1.push_front((*nll)->getVal());
    }
    w->loadSnapshot("GlobalMinimumParamValues");
    w->saveSnapshot("CurrentMinimumParamValues_MXs1", *params, true);

    // get asym error for mu_MXs2
    for (int i = 1; i < 15; i++) {
        w->loadSnapshot("CurrentMinimumParamValues_MXs2");
        w->var("mu_MXs2")->setVal(Global_mu_MXs2 + i * delta * mu_MXs2_HI_error);
        w->var("mu_MXs2")->setConstant(true);

        // fit with fixed mu
        MyMinim.minimize(minimizer, algorithm);

        // save snapshot
        w->saveSnapshot("CurrentMinimumParamValues_MXs2", *params, true);

        // save result
        mu_MXs2.push_back(Global_mu_MXs2 + i * delta * mu_MXs2_HI_error);
        MinusLogLikelihood_MXs2.push_back((*nll)->getVal());
    }
    w->loadSnapshot("GlobalMinimumParamValues");
    w->saveSnapshot("CurrentMinimumParamValues_MXs2", *params, true);

    for (int i = 1; i < 15; i++) {
        w->loadSnapshot("CurrentMinimumParamValues_MXs2");
        w->var("mu_MXs2")->setVal(Global_mu_MXs2 + i * delta * mu_MXs2_LO_error);
        w->var("mu_MXs2")->setConstant(true);

        // fit with fixed mu
        MyMinim.minimize(minimizer, algorithm);

        // save snapshot
        w->saveSnapshot("CurrentMinimumParamValues_MXs2", *params, true);

        // save result
        mu_MXs2.push_front(Global_mu_MXs2 + i * delta * mu_MXs2_LO_error);
        MinusLogLikelihood_MXs2.push_front((*nll)->getVal());
    }
    w->loadSnapshot("GlobalMinimumParamValues");
    w->saveSnapshot("CurrentMinimumParamValues_MXs2", *params, true);

    // get asym error for mu_MXs3
    for (int i = 1; i < 15; i++) {
        w->loadSnapshot("CurrentMinimumParamValues_MXs3");
        w->var("mu_MXs3")->setVal(Global_mu_MXs3 + i * delta * mu_MXs3_HI_error);
        w->var("mu_MXs3")->setConstant(true);

        // fit with fixed mu
        MyMinim.minimize(minimizer, algorithm);

        // save snapshot
        w->saveSnapshot("CurrentMinimumParamValues_MXs3", *params, true);

        // save result
        mu_MXs3.push_back(Global_mu_MXs3 + i * delta * mu_MXs3_HI_error);
        MinusLogLikelihood_MXs3.push_back((*nll)->getVal());
    }
    w->loadSnapshot("GlobalMinimumParamValues");
    w->saveSnapshot("CurrentMinimumParamValues_MXs3", *params, true);

    for (int i = 1; i < 15; i++) {
        w->loadSnapshot("CurrentMinimumParamValues_MXs3");
        w->var("mu_MXs3")->setVal(Global_mu_MXs3 + i * delta * mu_MXs3_LO_error);
        w->var("mu_MXs3")->setConstant(true);

        // fit with fixed mu
        MyMinim.minimize(minimizer, algorithm);

        // save snapshot
        w->saveSnapshot("CurrentMinimumParamValues_MXs3", *params, true);

        // save result
        mu_MXs3.push_front(Global_mu_MXs3 + i * delta * mu_MXs3_LO_error);
        MinusLogLikelihood_MXs3.push_front((*nll)->getVal());
    }
    w->loadSnapshot("GlobalMinimumParamValues");
    w->saveSnapshot("CurrentMinimumParamValues_MXs3", *params, true);

    // calculate asym error
    double My_mu_MXs1_HI_error = -100;
    double My_mu_MXs2_HI_error = -100;
    double My_mu_MXs3_HI_error = -100;
    double My_mu_MXs1_LO_error = -100;
    double My_mu_MXs2_LO_error = -100;
    double My_mu_MXs3_LO_error = -100;

    for (int i = 0; i < mu_MXs1.size() - 1; i++) {
        double previous_profile_likelihood = MinusLogLikelihood_MXs1.at(i) - Global_MinusLogLikelihood;
        double current_profile_likelihood = MinusLogLikelihood_MXs1.at(i + 1) - Global_MinusLogLikelihood;
        double previous_mu = mu_MXs1.at(i);
        double current_mu = mu_MXs1.at(i + 1);
        if ((previous_profile_likelihood > 0.5) && (current_profile_likelihood < 0.5)) { // we just passed -1 sigma point
            double mu_nominal_minus_sigma = current_mu - (0.5 - current_profile_likelihood) * (current_mu - previous_mu) / (previous_profile_likelihood - current_profile_likelihood);
            double minus_sigma = mu_nominal_minus_sigma - Global_mu_MXs1; // it is negative value
            My_mu_MXs1_LO_error = minus_sigma;
        }
        else if ((previous_profile_likelihood < 0.5) && (current_profile_likelihood > 0.5)) { // we just passed +1 sigma point
            double mu_nominal_plus_sigma = previous_mu + (0.5 - previous_profile_likelihood) * (current_mu - previous_mu) / (current_profile_likelihood - previous_profile_likelihood);
            double plus_sigma = mu_nominal_plus_sigma - Global_mu_MXs1; // it is positive value
            My_mu_MXs1_HI_error = plus_sigma;
        }
    }

    for (int i = 0; i < mu_MXs2.size() - 1; i++) {
        double previous_profile_likelihood = MinusLogLikelihood_MXs2.at(i) - Global_MinusLogLikelihood;
        double current_profile_likelihood = MinusLogLikelihood_MXs2.at(i + 1) - Global_MinusLogLikelihood;
        double previous_mu = mu_MXs2.at(i);
        double current_mu = mu_MXs2.at(i + 1);
        if ((previous_profile_likelihood > 0.5) && (current_profile_likelihood < 0.5)) { // we just passed -1 sigma point
            double mu_nominal_minus_sigma = current_mu - (0.5 - current_profile_likelihood) * (current_mu - previous_mu) / (previous_profile_likelihood - current_profile_likelihood);
            double minus_sigma = mu_nominal_minus_sigma - Global_mu_MXs2; // it is negative value
            My_mu_MXs2_LO_error = minus_sigma;
        }
        else if ((previous_profile_likelihood < 0.5) && (current_profile_likelihood > 0.5)) { // we just passed +1 sigma point
            double mu_nominal_plus_sigma = previous_mu + (0.5 - previous_profile_likelihood) * (current_mu - previous_mu) / (current_profile_likelihood - previous_profile_likelihood);
            double plus_sigma = mu_nominal_plus_sigma - Global_mu_MXs2; // it is positive value
            My_mu_MXs2_HI_error = plus_sigma;
        }
    }

    for (int i = 0; i < mu_MXs3.size() - 1; i++) {
        double previous_profile_likelihood = MinusLogLikelihood_MXs3.at(i) - Global_MinusLogLikelihood;
        double current_profile_likelihood = MinusLogLikelihood_MXs3.at(i + 1) - Global_MinusLogLikelihood;
        double previous_mu = mu_MXs3.at(i);
        double current_mu = mu_MXs3.at(i + 1);
        if ((previous_profile_likelihood > 0.5) && (current_profile_likelihood < 0.5)) { // we just passed -1 sigma point
            double mu_nominal_minus_sigma = current_mu - (0.5 - current_profile_likelihood) * (current_mu - previous_mu) / (previous_profile_likelihood - current_profile_likelihood);
            double minus_sigma = mu_nominal_minus_sigma - Global_mu_MXs3; // it is negative value
            My_mu_MXs3_LO_error = minus_sigma;
        }
        else if ((previous_profile_likelihood < 0.5) && (current_profile_likelihood > 0.5)) { // we just passed +1 sigma point
            double mu_nominal_plus_sigma = previous_mu + (0.5 - previous_profile_likelihood) * (current_mu - previous_mu) / (current_profile_likelihood - previous_profile_likelihood);
            double plus_sigma = mu_nominal_plus_sigma - Global_mu_MXs3; // it is positive value
            My_mu_MXs3_HI_error = plus_sigma;
        }
    }

    printf("=======================================\n");
    printf("Let's compare MINOS and homemade error\n");
    printf("MINOS:\n");
    printf("mu MXs1: %lf %lf\n", mu_MXs1_HI_error, mu_MXs1_LO_error);
    printf("mu MXs2: %lf %lf\n", mu_MXs2_HI_error, mu_MXs2_LO_error);
    printf("mu MXs3: %lf %lf\n", mu_MXs3_HI_error, mu_MXs3_LO_error);
    printf("homemade:\n");
    printf("mu MXs1: %lf %lf\n", My_mu_MXs1_HI_error, My_mu_MXs1_LO_error);
    printf("mu MXs2: %lf %lf\n", My_mu_MXs2_HI_error, My_mu_MXs2_LO_error);
    printf("mu MXs3: %lf %lf\n", My_mu_MXs3_HI_error, My_mu_MXs3_LO_error);
    printf("-Log(L) scan:\n");
    printf("mu MXs1: ");
    for(int i = 0; i < MinusLogLikelihood_MXs1.size() - 1; i++) printf("%lf ", MinusLogLikelihood_MXs1.at(i));
    printf("\n");
    printf("mu MXs2: ");
    for (int i = 0; i < MinusLogLikelihood_MXs2.size() - 1; i++) printf("%lf ", MinusLogLikelihood_MXs2.at(i));
    printf("\n");
    printf("mu MXs3: ");
    for (int i = 0; i < MinusLogLikelihood_MXs3.size() - 1; i++) printf("%lf ", MinusLogLikelihood_MXs3.at(i));
    printf("=======================================\n");

    return minim.save();
}

void GetPlotTemplate(RooWorkspace* w, RooDataSet* data = nullptr) {

    bool Allchargednull = true;
    bool Allmixednull = true;
    bool Alluubarnull = true;
    bool Allddbarnull = true;
    bool Allssbarnull = true;
    bool Allccbarnull = true;
    bool AllSIGANLnull = true;

    THStack* Stack = new THStack("Stack", ";bin index;number of event");
    TH1D* charged_hist = new TH1D("charged", ";bin index;number of event", RarityBins, BinMIN, BinMAX);
    TH1D* mixed_hist = new TH1D("mixed", ";bin index;number of event", RarityBins, BinMIN, BinMAX);
    TH1D* uubar_hist = new TH1D("u#bar{u}", ";bin index;number of event", RarityBins, BinMIN, BinMAX);
    TH1D* ddbar_hist = new TH1D("d#bar{d}", ";bin index;number of event", RarityBins, BinMIN, BinMAX);
    TH1D* ssbar_hist = new TH1D("s#bar{s}", ";bin index;number of event", RarityBins, BinMIN, BinMAX);
    TH1D* ccbar_hist = new TH1D("c#bar{c}", ";bin index;number of event", RarityBins, BinMIN, BinMAX);
    TH1D* SIGNAL_hist = new TH1D("SIGNAL", ";bin index;number of event", RarityBins, BinMIN, BinMAX);
    TH1D* all_hist = new TH1D("all", ";bin index;number of event", RarityBins, BinMIN, BinMAX);
    TH1D* data_hist = nullptr;
    if (data != nullptr) {
        data_hist = new TH1D("data", ";bin index;number of event", RarityBins, BinMIN, BinMAX);
        data_hist->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    }
    TH1D* Ratio_hist = new TH1D("Ratio", ";bin index;data/MC", RarityBins, BinMIN, BinMAX);

    // fill histogram
    for (int i = 0; i < scaleFactors_pdf_names.size(); i++) {

        RooRealVar* x_val;
        if (std::strstr(scaleFactors_pdf_names.at(i).c_str(), "channel_MXs1") != nullptr) x_val = w->var("obs_x_channel_MXs1");
        else if (std::strstr(scaleFactors_pdf_names.at(i).c_str(), "channel_MXs2") != nullptr) x_val = w->var("obs_x_channel_MXs2");
        else if (std::strstr(scaleFactors_pdf_names.at(i).c_str(), "channel_MXs3") != nullptr) x_val = w->var("obs_x_channel_MXs3");
        else {
            printf("[ERROR] unexpected sample type!\n");
            exit(1);
        }

        RooAbsBinning const& binning = x_val->getBinning();
        const double oldVal = x_val->getVal();

        TH1D* temp_hist;
        if (std::strstr(scaleFactors_pdf_names.at(i).c_str(), "Signal") != nullptr) temp_hist = SIGNAL_hist;
        else if (std::strstr(scaleFactors_pdf_names.at(i).c_str(), "CHG") != nullptr) temp_hist = charged_hist;
        else if (std::strstr(scaleFactors_pdf_names.at(i).c_str(), "MIX") != nullptr) temp_hist = mixed_hist;
        else if (std::strstr(scaleFactors_pdf_names.at(i).c_str(), "UUBAR") != nullptr) temp_hist = uubar_hist;
        else if (std::strstr(scaleFactors_pdf_names.at(i).c_str(), "DDBAR") != nullptr) temp_hist = ddbar_hist;
        else if (std::strstr(scaleFactors_pdf_names.at(i).c_str(), "SSBAR") != nullptr) temp_hist = ssbar_hist;
        else if (std::strstr(scaleFactors_pdf_names.at(i).c_str(), "CHARM") != nullptr) temp_hist = ccbar_hist;
        else {
            printf("[ERROR] unexpected sample type!\n");
            exit(1);
        }

        for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
            double binCenter = binning.binCenter(iBin);
            double binWidth = binning.binWidth(iBin);

            *x_val = binCenter; // set x value

            RooAbsReal* temp_func_scaleFactors = w->function(scaleFactors_pdf_names.at(i).c_str());
            RooAbsReal* temp_func_shapes = w->function(shapes_pdf_names.at(i).c_str());
            if ((temp_func_scaleFactors == nullptr) || (temp_func_shapes == nullptr)) {
                printf("[WARNING] cannot find %s or %s. Just skip.\n", scaleFactors_pdf_names.at(i).c_str(), shapes_pdf_names.at(i).c_str());
                break;
            }
            else {
                if (std::strstr(scaleFactors_pdf_names.at(i).c_str(), "Signal") != nullptr) AllSIGANLnull = false;
                else if (std::strstr(scaleFactors_pdf_names.at(i).c_str(), "CHG") != nullptr) Allchargednull = false;
                else if (std::strstr(scaleFactors_pdf_names.at(i).c_str(), "MIX") != nullptr) Allmixednull = false;
                else if (std::strstr(scaleFactors_pdf_names.at(i).c_str(), "UUBAR") != nullptr) Alluubarnull = false;
                else if (std::strstr(scaleFactors_pdf_names.at(i).c_str(), "DDBAR") != nullptr) Allddbarnull = false;
                else if (std::strstr(scaleFactors_pdf_names.at(i).c_str(), "SSBAR") != nullptr) Allssbarnull = false;
                else if (std::strstr(scaleFactors_pdf_names.at(i).c_str(), "CHARM") != nullptr) Allccbarnull = false;
            }
            double Nevt = (temp_func_scaleFactors->getValV() * temp_func_shapes->getValV());

            int index = -1;
            if (std::strstr(scaleFactors_pdf_names.at(i).c_str(), "channel_MXs1") != nullptr) index = iBin + 1;
            else if (std::strstr(scaleFactors_pdf_names.at(i).c_str(), "channel_MXs2") != nullptr) index = iBin + RarityBins_MX1 + 1;
            else if (std::strstr(scaleFactors_pdf_names.at(i).c_str(), "channel_MXs3") != nullptr) index = iBin + RarityBins_MX1 + RarityBins_MX2 + 1;
            else {
                printf("[ERROR] unexpected sample type!\n");
                exit(1);
            }

            temp_hist->SetBinContent(index, temp_hist->GetBinContent(index) + Nevt);
            all_hist->Fill(((double)index) - 0.5, Nevt);
            all_hist->SetBinError(index, 0);
        }

        *x_val = oldVal;

    }

    // Fill data
    if (data != nullptr) {
        for (int i = 0; i < RarityBins; i++) {
            const RooArgSet* argSet = data->get(i);
            data_hist->SetBinContent(i + 1, data->weight());
        }
    }

    // fill stack
    if (Allchargednull == false) Stack->Add(charged_hist);
    if (Allmixednull == false) Stack->Add(mixed_hist);
    if (Alluubarnull == false) Stack->Add(uubar_hist);
    if (Allddbarnull == false) Stack->Add(ddbar_hist);
    if (Allssbarnull == false) Stack->Add(ssbar_hist);
    if (Allccbarnull == false) Stack->Add(ccbar_hist);
    if (AllSIGANLnull == false) Stack->Add(SIGNAL_hist);

    // fill ratio
    Ratio_hist->SetLineColor(kBlack); Ratio_hist->SetMarkerStyle(21); Ratio_hist->Sumw2(); Ratio_hist->SetStats(0);
    Ratio_hist->Divide(data_hist, all_hist);

    // draw plot
    TCanvas* c_temp = new TCanvas("c", "", 800, 800); c_temp->cd();

    TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.35, 1.0, 1.0);
    pad1->SetBottomMargin(0.08); pad1->SetLeftMargin(0.15);
    pad1->Draw(); pad1->cd();

    gStyle->SetPalette(kPastel);

    Stack->Draw("pfc Hist");
    if (data != nullptr) {
        data_hist->SetLineWidth(2);
        data_hist->SetLineColor(kBlack);
        data_hist->SetMarkerStyle(8);
        data_hist->Draw("SAME eP EX0");
    }
    TLegend* legend = pad1->BuildLegend(0.95, 0.9, 0.75, 0.6);
    legend->SetFillStyle(0); legend->SetLineWidth(0);

    // vertical line to separate MXs region
    c_temp->Update();
    TLine* line_12 = new TLine((double)RarityBins_MX1, pad1->GetUymin(), (double)RarityBins_MX1, pad1->GetUymax());
    line_12->SetLineColor(kBlack); line_12->SetLineStyle(2); line_12->SetLineWidth(3);
    TLine* line_23 = new TLine((double)RarityBins_MX1 + RarityBins_MX2, pad1->GetUymin(), (double)RarityBins_MX1 + RarityBins_MX2, pad1->GetUymax());
    line_23->SetLineColor(kBlack); line_23->SetLineStyle(2); line_23->SetLineWidth(3);
    line_12->Draw(); line_23->Draw();

    // write MXs bin text
    TPaveText* pt_1 = new TPaveText(0.14, 0.9, 0.8 / 3.0 + 0.14, 1.0, "NDC NB");
    pt_1->SetTextSize(0.04); pt_1->SetFillStyle(0); pt_1->SetLineWidth(0); pt_1->AddText("0.0 < M_{X_{s}}^{reco} < 0.6 GeV/c^{2}"); pt_1->Draw();
    TPaveText* pt_2 = new TPaveText(0.8 / 3.0 + 0.14, 0.9, 2.0 * 0.8 / 3.0 + 0.14, 1.0, "NDC NB");
    pt_2->SetTextSize(0.04); pt_2->SetFillStyle(0); pt_2->SetLineWidth(0); pt_2->AddText("0.6 < M_{X_{s}}^{reco} < 1.0 GeV/c^{2}"); pt_2->Draw();
    TPaveText* pt_3 = new TPaveText(2.0 * 0.8 / 3.0 + 0.1, 0.9, 3.0 * 0.8 / 3.0 + 0.1, 1.0, "NDC NB");
    pt_3->SetTextSize(0.04); pt_3->SetFillStyle(0); pt_3->SetLineWidth(0); pt_3->AddText("1.0 GeV/c^{2} < M_{X_{s}}^{reco}"); pt_3->Draw();

    c_temp->cd();
    TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1, 0.3); pad2->SetBottomMargin(0.15); pad2->SetLeftMargin(0.15); pad2->Draw(); pad2->cd();
    Ratio_hist->SetMinimum(0.5); Ratio_hist->SetMaximum(1.5); Ratio_hist->SetLineWidth(2);
    Ratio_hist->GetYaxis()->SetTitleSize(0.08); Ratio_hist->GetYaxis()->SetTitleOffset(0.5);
    Ratio_hist->GetXaxis()->SetLabelSize(0.08); Ratio_hist->GetYaxis()->SetLabelSize(0.08);
    Ratio_hist->Draw("e0p");
    TLine* line = new TLine(Ratio_hist->GetXaxis()->GetXmin(), 1.0, Ratio_hist->GetXaxis()->GetXmax(), 1.0);
    line->SetLineColor(kRed);
    line->SetLineStyle(1); line->SetLineWidth(3);
    line->Draw();

    // vertical line to separate MXs region
    c_temp->Update();
    TLine* line_12_pad2 = new TLine((double)RarityBins_MX1, pad2->GetUymin(), (double)RarityBins_MX1, pad2->GetUymax());
    line_12_pad2->SetLineColor(kBlack); line_12_pad2->SetLineStyle(2); line_12_pad2->SetLineWidth(3);
    TLine* line_23_pad2 = new TLine((double)RarityBins_MX1 + RarityBins_MX2, pad2->GetUymin(), (double)RarityBins_MX1 + RarityBins_MX2, pad2->GetUymax());
    line_23_pad2->SetLineColor(kBlack); line_23_pad2->SetLineStyle(2); line_23_pad2->SetLineWidth(3);
    line_12_pad2->Draw(); line_23_pad2->Draw();

    c_temp->SetBottomMargin(0.0);
    c_temp->SaveAs("hist_data_plot.png");

    // print values
    printf("data:\n");
    for (int i = 0; i < RarityBins; i++) printf("%lf +- %lf\n", data_hist->GetBinContent(i + 1), data_hist->GetBinError(i + 1));
    printf("\nMC:\n");
    for (int i = 0; i < RarityBins; i++) printf("%lf +- %lf\n", all_hist->GetBinContent(i + 1), all_hist->GetBinError(i + 1));
    printf("\n");

    // delete
    delete c_temp;

    delete Stack;
    delete charged_hist;
    delete mixed_hist;
    delete uubar_hist;
    delete ddbar_hist;
    delete ssbar_hist;
    delete ccbar_hist;
    delete SIGNAL_hist;
    delete all_hist;
    if (data != nullptr) {
        delete data_hist;
    }
    delete Ratio_hist;

    delete line;

    delete line_12;
    delete line_23;

    delete line_12_pad2;
    delete line_23_pad2;
}

void ObtainNLL(RooWorkspace* w, RooDataSet* data, RooAbsReal** nll) {
    // what we have done
    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    // get nll
    RooArgSet* allParams = model->getParameters(*data);
    RooStats::RemoveConstantParameters(allParams);
    RooArgSet fGlobalObs = *mc->GetGlobalObservables();
    RooArgSet fConditionalObs;
    Bool_t fLOffset = RooStats::IsNLLOffset();
    (*nll) = model->createNLL(*data, RooFit::CloneData(kFALSE), RooFit::Constrain(*allParams), RooFit::GlobalObservables(fGlobalObs), RooFit::ConditionalObservables(fConditionalObs), RooFit::Offset(fLOffset));
}

#endif 