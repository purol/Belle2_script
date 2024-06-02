#ifndef FITTER_H
#define FITTER_H

# include <vector>
# include <string>
# include "RooFitResult.h"
# include "template.h"
# include "correctors.h"
# include <cstring>
# include "THStack.h"
# include "TStyle.h"
# include "TLine.h"
# include "TColor.h"

using namespace RooFit;
using namespace RooStats;
using namespace HistFactory;

Corrector_Fragmentation corrector_Fragmentation;

std::vector<std::string> Sample_names = {
    "L_x_Signal_MXs1_nominal_MXs1_channel_MXs1_overallSyst_x_StatUncert",
    "L_x_Signal_MXs2_nominal_MXs1_channel_MXs1_overallSyst_x_StatUncert",
    "L_x_Signal_MXs3_nominal_MXs1_channel_MXs1_overallSyst_x_StatUncert",
    "L_x_Signal_MXs1_nominal_MXs2_channel_MXs2_overallSyst_x_StatUncert",
    "L_x_Signal_MXs2_nominal_MXs2_channel_MXs2_overallSyst_x_StatUncert",
    "L_x_Signal_MXs3_nominal_MXs2_channel_MXs2_overallSyst_x_StatUncert",
    "L_x_Signal_MXs1_nominal_MXs3_channel_MXs3_overallSyst_x_StatUncert",
    "L_x_Signal_MXs2_nominal_MXs3_channel_MXs3_overallSyst_x_StatUncert",
    "L_x_Signal_MXs3_nominal_MXs3_channel_MXs3_overallSyst_x_StatUncert",
    "L_x_CHG_nominal_MXs1_channel_MXs1_overallSyst_x_StatUncert",
    "L_x_CHG_nominal_MXs2_channel_MXs2_overallSyst_x_StatUncert",
    "L_x_CHG_nominal_MXs3_channel_MXs3_overallSyst_x_StatUncert",
    "L_x_MIX_nominal_MXs1_channel_MXs1_overallSyst_x_StatUncert",
    "L_x_MIX_nominal_MXs2_channel_MXs2_overallSyst_x_StatUncert",
    "L_x_MIX_nominal_MXs3_channel_MXs3_overallSyst_x_StatUncert",
    "L_x_UUBAR_nominal_MXs1_channel_MXs1_overallSyst_x_StatUncert",
    "L_x_UUBAR_nominal_MXs2_channel_MXs2_overallSyst_x_StatUncert",
    "L_x_UUBAR_nominal_MXs3_channel_MXs3_overallSyst_x_StatUncert",
    "L_x_DDBAR_nominal_MXs1_channel_MXs1_overallSyst_x_StatUncert",
    "L_x_DDBAR_nominal_MXs2_channel_MXs2_overallSyst_x_StatUncert",
    "L_x_DDBAR_nominal_MXs3_channel_MXs3_overallSyst_x_StatUncert",
    "L_x_SSBAR_nominal_MXs1_channel_MXs1_overallSyst_x_StatUncert",
    "L_x_SSBAR_nominal_MXs2_channel_MXs2_overallSyst_x_StatUncert",
    "L_x_SSBAR_nominal_MXs3_channel_MXs3_overallSyst_x_StatUncert",
    "L_x_CHARM_nominal_MXs1_channel_MXs1_overallSyst_x_StatUncert",
    "L_x_CHARM_nominal_MXs2_channel_MXs2_overallSyst_x_StatUncert",
    "L_x_CHARM_nominal_MXs3_channel_MXs3_overallSyst_x_StatUncert"
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
    bool BRXnn;
    bool BRDKL0;
    bool uncorrelated;

    int NEntryFEI;
    int NEntryKID;
    int NEntryPID;
    int NEntryBR;
    int NEntrypi0;
    int NEntryMultiplicity;
    int NEntryFragmentation;

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

int ReadNBREigenVector(const char* dirname) {
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

int ReadNFragmentationEigenVector(const char* dirname) {
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
    options_->BRXnn = false;
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
        options_->BRBtoXKLKL = true;
        options_->BRXnn = true;
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
    else if (std::string(tested_param) == std::string("BRBtoXKLKL")) options_->BRBtoXKLKL = true;
    else if (std::string(tested_param) == std::string("BRXnn")) options_->BRXnn = true;
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
    options_->NEntryMultiplicity = ReadMultiplicityInfo("./multiplicity_selected.txt");
    options_->NEntryFragmentation = ReadNFragmentationEigenVector("./Fragmentation_selected.txt");
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
    //if (options_->qqbar) {
    //    w->var("alpha_qq_CAL_UUBAR_uncer")->setConstant(options_->qqbar);
    //    w->var("alpha_qq_CAL_DDBAR_uncer")->setConstant(options_->qqbar);
    //    w->var("alpha_qq_CAL_SSBAR_uncer")->setConstant(options_->qqbar);
    //    w->var("alpha_qq_CAL_CHARM_uncer")->setConstant(options_->qqbar);
    //}

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
    if (options_->Fragmentation) for (int i = 0; i < options_->NEntryFragmentation; i++) w->var(("alpha_Xs_fragmentation" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->Fragmentation);

    // mKstar
    if (options_->mKstar) w->var("alpha_mKstar_uncer")->setConstant(options_->mKstar);

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
    if (options_->BBBR) for (int i = 0; i < options_->NEntryBR; i++) w->var(("alpha_BBBR" + std::to_string(i) + "_uncer").c_str())->setConstant(options_->BBBR);

    // B->K KL KL BR
    if (options_->BRBtoXKLKL) w->var("alpha_BRBtoXKLKL_uncer")->setConstant(options_->BRBtoXKLKL);

    // B->Xnn BR
    if (options_->BRXnn) w->var("alpha_Xnn_BR_uncer")->setConstant(options_->BRXnn);

    // B->(D->X KL0) + anything
    if (options_->BRDKL0) w->var("alpha_BtoDtoXKL_uncer")->setConstant(options_->BRDKL0);

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

            RooAbsReal* temp_func = w->function(Sample_names.at(index).c_str());
            Nevt = Nevt + temp_func->getValV();
            if (temp_func->getValV() < 0) {
                printf("[ERROR] negative count!\n");
                exit(1);
            }

        }

        *x_val = oldVal;
    }

    return Nevt;

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
    minim.minos(RooArgSet(*w->var("mu_MXs1")));
    minim.minos(RooArgSet(*w->var("mu_MXs2")));
    minim.minos(RooArgSet(*w->var("mu_MXs3")));

    // fit!
    int status;
    status = minim.minimize(minimizer, algorithm);

    return minim.save();
}

void GetPlotTemplate(RooWorkspace* w, RooDataSet* data = nullptr) {

    THStack* Stack = new THStack("Stack", ";FBDT index;number of event");
    TH1D* charged_hist = new TH1D("charged", ";FBDT index;number of event", RarityBins, BinMIN, BinMAX);
    TH1D* mixed_hist = new TH1D("mixed", ";FBDT index;number of event", RarityBins, BinMIN, BinMAX);
    TH1D* uubar_hist = new TH1D("u#bar{u}", ";FBDT index;number of event", RarityBins, BinMIN, BinMAX);
    TH1D* ddbar_hist = new TH1D("d#bar{d}", ";FBDT index;number of event", RarityBins, BinMIN, BinMAX);
    TH1D* ssbar_hist = new TH1D("s#bar{s}", ";FBDT index;number of event", RarityBins, BinMIN, BinMAX);
    TH1D* ccbar_hist = new TH1D("c#bar{c}", ";FBDT index;number of event", RarityBins, BinMIN, BinMAX);
    TH1D* SIGNAL_hist = new TH1D("SIGNAL", ";FBDT index;number of event", RarityBins, BinMIN, BinMAX);
    TH1D* all_hist = new TH1D("all", ";FBDT index;number of event", RarityBins, BinMIN, BinMAX);
    TH1D* data_hist = nullptr;
    if (data != nullptr) {
        data_hist = new TH1D("data", ";FBDT index;number of event", RarityBins, BinMIN, BinMAX);
        data_hist->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    }
    TH1D* Ratio_hist = new TH1D("Ratio", ";FBDT index;data/MC", RarityBins, BinMIN, BinMAX);

    // fill histogram
    for (int i = 0; i < Sample_names.size(); i++) {

        RooRealVar* x_val;
        if (std::strstr(Sample_names.at(i).c_str(), "channel_MXs1") != nullptr) x_val = w->var("obs_x_channel_MXs1");
        else if (std::strstr(Sample_names.at(i).c_str(), "channel_MXs2") != nullptr) x_val = w->var("obs_x_channel_MXs2");
        else if (std::strstr(Sample_names.at(i).c_str(), "channel_MXs3") != nullptr) x_val = w->var("obs_x_channel_MXs3");
        else {
            printf("[ERROR] unexpected sample type!\n");
            exit(1);
        }

        RooAbsBinning const& binning = x_val->getBinning();
        const double oldVal = x_val->getVal();

        TH1D* temp_hist;
        if (std::strstr(Sample_names.at(i).c_str(), "Signal") != nullptr) temp_hist = SIGNAL_hist;
        else if (std::strstr(Sample_names.at(i).c_str(), "CHG") != nullptr) temp_hist = charged_hist;
        else if (std::strstr(Sample_names.at(i).c_str(), "MIX") != nullptr) temp_hist = mixed_hist;
        else if (std::strstr(Sample_names.at(i).c_str(), "UUBAR") != nullptr) temp_hist = uubar_hist;
        else if (std::strstr(Sample_names.at(i).c_str(), "DDBAR") != nullptr) temp_hist = ddbar_hist;
        else if (std::strstr(Sample_names.at(i).c_str(), "SSBAR") != nullptr) temp_hist = ssbar_hist;
        else if (std::strstr(Sample_names.at(i).c_str(), "CHARM") != nullptr) temp_hist = ccbar_hist;
        else {
            printf("[ERROR] unexpected sample type!\n");
            exit(1);
        }

        for (std::size_t iBin = 0; iBin < binning.numBins(); ++iBin) {
            double binCenter = binning.binCenter(iBin);
            double binWidth = binning.binWidth(iBin);

            *x_val = binCenter; // set x value

            RooAbsReal* temp_func = w->function(Sample_names.at(i).c_str());
            double Nevt = temp_func->getValV();

            int index = -1;
            if (std::strstr(Sample_names.at(i).c_str(), "channel_MXs1") != nullptr) index = iBin + 1;
            else if (std::strstr(Sample_names.at(i).c_str(), "channel_MXs2") != nullptr) index = iBin + RarityBins_MX1 + 1;
            else if (std::strstr(Sample_names.at(i).c_str(), "channel_MXs3") != nullptr) index = iBin + RarityBins_MX1 + RarityBins_MX2 + 1;
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
    Stack->Add(charged_hist);
    Stack->Add(mixed_hist);
    Stack->Add(uubar_hist);
    Stack->Add(ddbar_hist);
    Stack->Add(ssbar_hist);
    Stack->Add(ccbar_hist);
    Stack->Add(SIGNAL_hist);

    // fill ratio
    Ratio_hist->SetLineColor(kBlack); Ratio_hist->SetMarkerStyle(21); Ratio_hist->Sumw2(); Ratio_hist->SetStats(0);
    Ratio_hist->Divide(data_hist, all_hist);

    // draw plot
    TCanvas* c_temp = new TCanvas("c", "", 800, 800); c_temp->cd();

    TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.35, 1.0, 1.0);
    pad1->SetBottomMargin(0.08); pad1->SetLeftMargin(0.15);
    pad1->SetGridx(); pad1->Draw(); pad1->cd();

    gStyle->SetErrorX(0.0);
    gStyle->SetEndErrorSize(0.0);
    gStyle->SetPalette(kPastel);

    Stack->Draw("pfc Hist");
    if (data != nullptr) {
        data_hist->SetLineWidth(2);
        data_hist->SetLineColor(kBlack);
        data_hist->SetMarkerStyle(8);
        data_hist->Draw("SAME eP");
    }
    TLegend* legend = pad1->BuildLegend(0.95, 0.9, 0.75, 0.6);
    legend->SetFillStyle(0); legend->SetLineWidth(0);

    c_temp->cd();
    TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1, 0.3); pad2->SetBottomMargin(0.15); pad2->SetLeftMargin(0.15); pad2->SetGridx(); pad2->Draw(); pad2->cd();
    Ratio_hist->SetMinimum(0.5); Ratio_hist->SetMaximum(1.5); Ratio_hist->SetLineWidth(2);
    Ratio_hist->GetYaxis()->SetTitleSize(0.08); Ratio_hist->GetYaxis()->SetTitleOffset(0.5);
    Ratio_hist->GetXaxis()->SetLabelSize(0.08); Ratio_hist->GetYaxis()->SetLabelSize(0.08);
    Ratio_hist->Draw("e0p");
    TLine* line = new TLine(Ratio_hist->GetXaxis()->GetXmin(), 1.0, Ratio_hist->GetXaxis()->GetXmax(), 1.0);
    line->SetLineColor(kRed);
    line->SetLineStyle(1); line->SetLineWidth(3);
    line->Draw();

    c_temp->SetBottomMargin(0.0);
    c_temp->SaveAs("hist_data_plot.png");

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