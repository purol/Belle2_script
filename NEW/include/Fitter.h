#ifndef FITTER_H
#define FITTER_H

# include <vector>
# include <string>
# include "RooFitResult.h"
# include "template.h"
#include "correctors.h"

using namespace RooFit;
using namespace RooStats;
using namespace HistFactory;

Corrector_Fragmentation corrector_Fragmentation;

std::vector<std::string> Sample_names = {
    "L_x_Signal_nominal_MXs1_channel_overallSyst_x_StatUncert",
    "L_x_Signal_nominal_MXs2_channel_overallSyst_x_StatUncert",
    "L_x_Signal_nominal_MXs3_channel_overallSyst_x_StatUncert",
    "L_x_CHG_nominal_channel_overallSyst_x_StatUncert",
    "L_x_MIX_nominal_channel_overallSyst_x_StatUncert",
    "L_x_UUBAR_nominal_channel_overallSyst_x_StatUncert",
    "L_x_DDBAR_nominal_channel_overallSyst_x_StatUncert",
    "L_x_SSBAR_nominal_channel_overallSyst_x_StatUncert",
    "L_x_CHARM_nominal_channel_overallSyst_x_StatUncert"
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
    bool OLD_Kff;
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
    options_->OLD_Kff = false;
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
        options_->OLD_Kff = true;
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
    else if (std::string(tested_param) == std::string("OLDKff")) options_->OLD_Kff = true;
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
}

void FixParameters(RooWorkspace* w, OPTIONS* options_) {

    ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig"); // Get model manually
    RooSimultaneous* model = (RooSimultaneous*)mc->GetPdf();

    w->loadSnapshot("NominalParamValues");

    RooRealVar* x_val = w->var("obs_x_channel");
    std::unique_ptr<RooArgSet> params{model->getParameters(*x_val)};

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
    if (options_->FBDT) w->var("alpha_FBDT_efficiency_uncer")->setConstant(options_->FBDT);

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

    // new B->K form factor
    if (options_->OLD_Kff) w->var("alpha_Kff_OLD_uncer")->setConstant(options_->OLD_Kff);

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
    if (options_->MCstat) for (int i = 0; i < RarityBins; i++) w->var(("gamma_stat_channel_bin_" + std::to_string(i)).c_str())->setConstant(options_->MCstat);

    // Fragmentation
    if (options_->Fragmentation) {
        for (int MxsBin = 0; MxsBin < corrector_Fragmentation.GetNMxsBin(Corrector_Fragmentation::Sample::gamma); MxsBin++) {
            for (int Category = 0; Category < corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma); Category++) {
                int temp_index = MxsBin * corrector_Fragmentation.GetNCategory(Corrector_Fragmentation::Sample::gamma) + Category;
                w->var(("alpha_Xs_fragmentation" + std::to_string(temp_index) + "_uncer").c_str())->setConstant(options_->Fragmentation);
            }
        }
    }

    // mKstar
    if (options_->mKstar) w->var("alpha_mKstar_uncer")->setConstant(options_->mKstar);

    // background normalization
    if (options_->BKGNorm) {
        w->var("alpha_mu_CHG")->setConstant(options_->BKGNorm);
        w->var("alpha_mu_MIX")->setConstant(options_->BKGNorm);
        w->var("alpha_mu_UUBAR")->setConstant(options_->BKGNorm);
        w->var("alpha_mu_DDBAR")->setConstant(options_->BKGNorm);
        w->var("alpha_mu_SSBAR")->setConstant(options_->BKGNorm);
        w->var("alpha_mu_CHARM")->setConstant(options_->BKGNorm);
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
    RooRealVar* x_val = w->var("obs_x_channel");

    int index = -1;
    if (strcmp(sample_type, "Signal_MX1") == 0) index = 0;
    else if (strcmp(sample_type, "Signal_MX2") == 0) index = 1;
    else if (strcmp(sample_type, "Signal_MX3") == 0) index = 2;
    else if (strcmp(sample_type, "CHG") == 0) index = 3;
    else if (strcmp(sample_type, "MIX") == 0) index = 4;
    else if (strcmp(sample_type, "UUBAR") == 0) index = 5;
    else if (strcmp(sample_type, "DDBAR") == 0) index = 6;
    else if (strcmp(sample_type, "SSBAR") == 0) index = 7;
    else if (strcmp(sample_type, "CHARM") == 0) index = 8;
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