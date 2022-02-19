#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include <cmath>
using namespace RooFit ;

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

// https://confluence.desy.de/pages/viewpage.action?pageId=107054222
# define N_BpBp_1invab 565400000.0
# define N_B0B0_1invab 534600000.0

# define N_Kplus_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Kplus_nunubar)
# define N_Kplusstar_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Kplusstar_nunubar)
# define N_Xsu_nonresonant_nunubar_1invab (2.0 * N_BpBp_1invab * BR_Xsu_nonresonant_nunubar)
# define N_K0_nunubar_1invab (2.0 * N_B0B0_1invab * BR_K0_nunubar)
# define N_K0star_nunubar_1invab (2.0 * N_B0B0_1invab * BR_K0star_nunubar)
# define N_Xsd_nunubar_1invab (2.0 * N_B0B0_1invab * BR_Xsd_nonresonant_nunubar)

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

# define EeclBins 15

std::vector<double> Ns;
std::vector<int> ntracks;
std::vector<int> npi0s;
std::vector<double> KS0_3D_distance;

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

RooRealVar  Mbc_DATA("Mbc", "Mbc_DATA", 5.27, 5.29);
RooRealVar  Eecl_DATA("Eecl", "Eecl_DATA", 0, 4);
RooRealVar weight_DATA("weight", "weight_DATA", 0.0, 1.0);
RooDataSet info_DATA("2Dinfo", "2Dinfo_DATA", RooArgSet(Mbc_DATA, Eecl_DATA, weight_DATA), WeightVar("weight"));

RooRealVar  Mbc_MC_signal("Mbc", "Mbc_MC_signal", 5.27, 5.29);
RooRealVar  Eecl_MC_signal("Eecl", "Eecl_MC_signal", 0, 4);
RooRealVar weight_MC_signal("weight", "weight_MC_signal", 0.0, 1.0);
RooDataSet info_MC_signal("2Dinfo", "2Dinfo_MC_signal", RooArgSet(Mbc_MC_signal, Eecl_MC_signal, weight_MC_signal), WeightVar("weight"));

RooRealVar  Mbc_MC_background("Mbc", "Mbc_MC_background", 5.27, 5.29);
RooRealVar  Eecl_MC_background("Eecl", "Eecl_MC_background", 0, 4);
RooRealVar weight_MC_background("weight", "weight_MC_background", 0.0, 1.0);
RooDataSet info_MC_background("2Dinfo", "2Dinfo_MC_background", RooArgSet(Mbc_MC_background, Eecl_MC_background, weight_MC_background), WeightVar("weight"));

void LetsAdd(const char* dirname, RooRealVar* Mbc_, RooRealVar*  Eecl_, RooRealVar* weight_, RooDataSet* info_, double weight_var = 1.0) {
    double Eecl_var = 0;
    double Mbc_var = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("roeEextra__bocleanMask__bc", &Eecl_var); // Eecl
        tree_Btag->SetBranchAddress("Btag_Mbc", &Mbc_var); // Mbc

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);
            
            *Eecl_ = Eecl_var;
            *Mbc_ = Mbc_var;
            *weight_ = weight_var;
            info_->add(RooArgSet(*Mbc_, *Eecl_), weight_->getVal());
        }
        input_file->Close();

    }

}

void LetsCalculateUncertainties(const char* dirname, double weight) {

    double Upsilon_ID = -1;
    double Bsig_ID = -1;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        double temp_KS0_3D_distance = -1;

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        tree_upsilon->SetBranchAddress("extraInfo__bodecayModeID__bc", &Upsilon_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_decayModeID", &Bsig_ID);
        tree_Bsig->SetBranchAddress("Bsig_daughter_0_extraInfo_KS0_3D_distance", &temp_KS0_3D_distance);

        printf("%lld entries...\n", tree_upsilon->GetEntries());
        for (unsigned int j = 0; j < tree_upsilon->GetEntries(); j++) { // Fill
            tree_upsilon->GetEntry(j);
            tree_Bsig->GetEntry(j);
            tree_Btag->GetEntry(j);

            if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) { // B2Kc
                ntracks.push_back(1);
                npi0s.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) { // B2KcPi0
                ntracks.push_back(1);
                npi0s.push_back(1);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) { // B2Ks0Pic
                ntracks.push_back(1);
                npi0s.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) { // B2KcPicPic
                ntracks.push_back(3);
                npi0s.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) { // B2Ks0PicPi0
                ntracks.push_back(1);
                npi0s.push_back(1);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) { // B2KcPicPicPi0
                ntracks.push_back(3);
                npi0s.push_back(1);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) { // B2Ks0PicPicPic
                ntracks.push_back(3);
                npi0s.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) { // B2KcKcKc
                ntracks.push_back(3);
                npi0s.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) { // B2KcKcKs0Pic
                ntracks.push_back(3);
                npi0s.push_back(0);
            }
            else if (Upsilon_ID > -0.5 && Upsilon_ID < 0.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) { // B2KcKcKcPi0
                ntracks.push_back(3);
                npi0s.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > -0.5 && Bsig_ID < 0.5) { // B02Ks0
                ntracks.push_back(0);
                npi0s.push_back(0);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 0.5 && Bsig_ID < 1.5) { // B02KcPic
                ntracks.push_back(2);
                npi0s.push_back(0);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 1.5 && Bsig_ID < 2.5) { // B02Ks0Pi0
                ntracks.push_back(0);
                npi0s.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 2.5 && Bsig_ID < 3.5) { // B02KcPicPi0
                ntracks.push_back(2);
                npi0s.push_back(1);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 3.5 && Bsig_ID < 4.5) { // B02Ks0PicPic
                ntracks.push_back(2);
                npi0s.push_back(0);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 4.5 && Bsig_ID < 5.5) { // B02KcPicPicPic
                ntracks.push_back(4);
                npi0s.push_back(0);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 5.5 && Bsig_ID < 6.5) { // B02Ks0PicPicPi0
                ntracks.push_back(2);
                npi0s.push_back(0);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 11.5 && Bsig_ID < 12.5) { // B02KcKcKs0
                ntracks.push_back(2);
                npi0s.push_back(0);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 12.5 && Bsig_ID < 13.5) { // B02KcKcKcPic
                ntracks.push_back(4);
                npi0s.push_back(0);
            }
            else if (Upsilon_ID > 0.5 && Upsilon_ID < 1.5 && Bsig_ID > 13.5 && Bsig_ID < 14.5) { // B02KcKcKs0Pi0
                ntracks.push_back(2);
                npi0s.push_back(1);
            }
            else {
                printf("[ERROR] unexpected decay ID\n");
                exit(1);
            }
            KS0_3D_distance.push_back(temp_KS0_3D_distance);

            Ns.push_back(weight);
        }
        input_file->Close();

    }
}

void PrintUncertainties() {
    // start to calculate correction factor/uncertainties
    std::vector<double> corrected_Ns;
    std::vector<double> KS0_rel_uncertainties;
    std::vector<double> track_rel_uncertainties;
    std::vector<double> pi0_rel_uncertainties;
    for (unsigned int j = 0; j < Ns.size(); j++) {
        corrected_Ns.push_back(pow(pi0_correction, npi0s.at(j)) * Ns.at(j) );
        track_rel_uncertainties.push_back(track_rel_uncertainty * ntracks.at(j));
        pi0_rel_uncertainties.push_back(pi0_rel_uncertainty * npi0s.at(j));
        KS0_rel_uncertainties.push_back(KS0_rel_uncertainty * KS0_3D_distance.at(j));
    }
    double corrected_N = 0;
    double avg_track_rel_uncertainty = 0;
    double avg_pi0_rel_uncertainty = 0;
    double avg_KS0_rel_uncertainty = 0;

    for (unsigned int j = 0; j < Ns.size(); j++) {
        corrected_N = corrected_N + corrected_Ns.at(j);
        avg_track_rel_uncertainty = avg_track_rel_uncertainty + corrected_Ns.at(j) * track_rel_uncertainties.at(j);
        avg_pi0_rel_uncertainty = avg_pi0_rel_uncertainty + corrected_Ns.at(j) * pi0_rel_uncertainties.at(j);
        avg_KS0_rel_uncertainty = avg_KS0_rel_uncertainty + corrected_Ns.at(j) * KS0_rel_uncertainties.at(j);
    }
    avg_track_rel_uncertainty = avg_track_rel_uncertainty / corrected_N;
    avg_pi0_rel_uncertainty = avg_pi0_rel_uncertainty / corrected_N;
    avg_KS0_rel_uncertainty = avg_KS0_rel_uncertainty / corrected_N;

    double total_N = 0;
    for (unsigned int j = 0; j < Ns.size(); j++) total_N = total_N + Ns.at(j);

    printf("total uncorrected signal num: %lf\n", total_N);
    printf("Average correction factor: %lf\n", corrected_N / total_N);
    printf("Average relative uncertainty from track: %lf%%\n", avg_track_rel_uncertainty);
    printf("Average relative uncertainty from pi0: %lf%%\n", avg_pi0_rel_uncertainty);
    printf("Average relative uncertainty from KS0: %lf%%\n", avg_KS0_rel_uncertainty);
}

double GetEvtNum(const char* dirname, double weight_var = 1.0) {

    double Total_Evt_Num = 0;

    std::vector<string> names;
    load_files(dirname, &names);

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");

        TTree* tree_upsilon = (TTree*)input_file->Get("Upsilon");
        TTree* tree_Bsig = (TTree*)input_file->Get("Bsig");
        TTree* tree_Btag = (TTree*)input_file->Get("Btag");

        Total_Evt_Num = Total_Evt_Num + weight_var * tree_upsilon->GetEntries();

        input_file->Close();

    }

    return Total_Evt_Num;

}

void LinearityTest(RooFitResult* r, RooRealVar EeclFit, double BKG_num) {
    const int LT_number = 50;
    const int LT_iterate_number = 1000;

    RooArgSet fitargs = r->floatParsFinal();
    TIterator* iter(fitargs.createIterator());

    double N_nbkg = -1;
    double N_nsig = -1;
    double N_p2 = -1;
    double ERR_nbkg = -1;
    double ERR_nsig = -1;
    double ERR_p2 = -1;

    for (TObject* a = iter->Next(); a != 0; a = iter->Next()) {
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(a);
        std::string name = rrv->GetName();
        double val = rrv->getVal();
        double err = rrv->getError();
        if (name == std::string("nbkg")) {
            N_nbkg = val;
            ERR_nbkg = err;
        }
        else if (name == std::string("nsig")) {
            N_nsig = val;
            ERR_nsig = err;
        }
        else if (name == std::string("p2")) {
            N_p2 = val;
            ERR_p2 = err;
        }
    }

    std::vector<double> n_sigs[LT_number];
    std::vector<double> n_sigs_err[LT_number];
    for (int i = 0; i < LT_number; i++) {

        for (int j = 0; j < LT_iterate_number; j++) {
            // produce toy MC sample
            RooRealVar m0_fix("m0_fix", "m0_fix", 0.65);
            RooRealVar sigma_fix("sigma_fix", "sigma_fix", 0.53);
            RooRealVar alpha_fix("alpha_fix", "alpha_fix", -0.85);
            RooRealVar n_fix("n_fix", "n_fix", 60);
            RooCBShape histpdf_Eecl_signal_fix("genpdfs_fix", "genpdfs_fix", EeclFit, m0_fix, sigma_fix, alpha_fix, n_fix);
            RooDataSet* d1 = histpdf_Eecl_signal_fix.generate(EeclFit, i);

            RooRealVar m0b_fix("m0b_fix", "m0b_fix", 1.6);
            RooRealVar sigmab_fix("sigmab_fix", "sigmab_fix", 0.67);
            RooRealVar alphab_fix("alphab_fix", "alphab_fix", -0.368);
            RooRealVar nb_fix("nb_fix", "nb_fix", 20);
            RooCBShape histpdf_Eecl_background_fix("genpdfb_fix", "genpdfb_fix", EeclFit, m0b_fix, sigmab_fix, alphab_fix, nb_fix);
            int N_nbkg_int = round(BKG_num);
            RooDataSet* genData = histpdf_Eecl_background_fix.generate(EeclFit, N_nbkg_int);

            genData->append(*d1);
//            RooDataHist gen_binned_data_Eecl("gen binned Eecl data", "gen binned Eecl data", EeclFit, *genData);

            // construct fitting function
            RooRealVar m0_LT("m0_LT", "m0_LT", 0.65);
            RooRealVar sigma_LT("sigma_LT", "sigma_LT", 0.53);
            RooRealVar alpha_LT("alpha_LT", "alpha_LT", -0.85);
            RooRealVar n_LT("n_LT", "n_LT", 60);
            RooCBShape histpdf_Eecl_signal_LT("genpdfs_LT", "genpdfs_LT", EeclFit, m0_LT, sigma_LT, alpha_LT, n_LT);
            RooRealVar nsig("nsig", "n_{sig}", 7.5, -100, 150);
            RooExtendPdf esig("esignal", "extended signal p.d.f", histpdf_Eecl_signal_LT, nsig);

            RooRealVar m0b_LT("m0b_LT", "m0b_LT", 1.6);
            RooRealVar sigmab_LT("sigmab_LT", "sigmab_LT", 0.67);
            RooRealVar alphab_LT("alphab_LT", "alphab_LT", -0.368);
            RooRealVar nb_LT("nb_LT", "nb_LT", 20);
            RooCBShape histpdf_Eecl_background_LT("genpdfb_LT", "genpdfb_LT", EeclFit, m0b_LT, sigmab_LT, alphab_LT, nb_LT);

            RooRealVar nbkg("nbkg", "number of background events", 410, 200, 600);
            RooExtendPdf ebkg("ebkg", "extended background p.d.f", histpdf_Eecl_background_LT, nbkg);

            RooAddPdf  totalpdf("model", "b+n", RooArgList(ebkg, esig));

            // fit
            RooFitResult* result = totalpdf.fitTo(*genData, Verbose(false), PrintLevel(-1),  Save());

            RooArgSet fitargs_LT = result->floatParsFinal();
            TIterator* iter_LT(fitargs_LT.createIterator());

            for (TObject* a_LT = iter_LT->Next(); a_LT != 0; a_LT = iter_LT->Next()) {
                RooRealVar* rrv_LT = dynamic_cast<RooRealVar*>(a_LT);
                std::string name_LT = rrv_LT->GetName();
                double val_LT = rrv_LT->getVal();
                double err_LT = rrv_LT->getError();
                if (name_LT == std::string("nsig")) {
                    n_sigs[i].push_back(val_LT);
                    n_sigs_err[i].push_back(err_LT);
                }
            }

        }

    }

    // create root file
    TFile* temp_file = new TFile("Linearity.root", "recreate");
    temp_file->cd();
    TTree* temp_tree = new TTree("Linearity", "");
    double nsig_fit[LT_number];
    double nsig_fit_err[LT_number];

    for (int i = 0; i < LT_number;i++) temp_tree->Branch(("nsig_"+ to_string(i)).c_str(), &nsig_fit[i]);
    for (int i = 0; i < LT_number; i++) temp_tree->Branch(("nsig_error_" + to_string(i)).c_str(), &nsig_fit_err[i]);

    for (int i = 0; i < LT_iterate_number; i++) {
        for (int j = 0; j < LT_number; j++) {
            nsig_fit[j] = n_sigs[j].at(i);
            nsig_fit_err[j] = n_sigs_err[j].at(i);
        }
        temp_tree->Fill();
    }
    temp_tree->Write();
    temp_file->Close();

    // print png file
    double Inputnsig[LT_number] = { 0 };
    double Inputnsigerror[LT_number] = { 0 };
    double outputnsig[LT_number] = { 0 };
    double outputnsigerror[LT_number] = { 0 };
    for (int i = 0; i < LT_number; i++) {
        RooRealVar  nsig_roorealvar("nsig_roorealvar", "n_{sig}", -100, 100);
        RooDataSet nsig_RooDataSet("nsig_RooDataSet", "nsig_RooDataSet", RooArgSet(nsig_roorealvar));
        for (int j = 0; j < LT_iterate_number; j++) {
            nsig_roorealvar = n_sigs[i].at(j);
            nsig_RooDataSet.add(RooArgSet(nsig_roorealvar));
        }
        RooRealVar gausmean("gausmean","",i, i - ERR_nsig, i + ERR_nsig);
        RooRealVar gauswidth("gauswidth","", ERR_nsig, 0, 2* ERR_nsig);
        RooGaussian gauss("gauss","gauss", nsig_roorealvar,gausmean,gauswidth);
        RooRealVar nentry("nentry", "number of entries", 1000, 900, 1100);
        RooExtendPdf egauss("egauss", "extended gauss", gauss, nentry);
        RooFitResult* result_LT = egauss.fitTo(nsig_RooDataSet, Save());

        RooPlot* nsigframe = nsig_roorealvar.frame(Bins(400), Title(" "));
        nsig_RooDataSet.plotOn(nsigframe);
        egauss.plotOn(nsigframe, LineColor(kBlue));
        egauss.paramOn(nsigframe);
        TCanvas* c = new TCanvas("nsig_gauss_fit", "nsig_gauss_fit", 600, 600);
        gPad->SetLeftMargin(0.15); nsigframe->GetYaxis()->SetTitleOffset(1.4); nsigframe->Draw(); c->SaveAs(("nsig_"+ to_string(i) +"_distribution.png").c_str());

        Inputnsig[i] = i;
        Inputnsigerror[i] = 0;
        outputnsig[i] = gausmean.getValV();
        outputnsigerror[i] = gauswidth.getValV();

        delete c;
    }
    gStyle->SetOptFit(010); gStyle->SetStatH(0.05);
    TCanvas* c = new TCanvas("Linearity test canvas", "", 800, 800);
    TGraphErrors* gr = new TGraphErrors(LT_number, Inputnsig, outputnsig, Inputnsigerror, outputnsigerror);
    gr->SetMarkerStyle(21); gr->SetTitle(";input n_{sig};output n_{sig}");
    gr->Fit("pol1"); gr->Draw("AP");
    c->SaveAs("Linearity_test.png");

    gr->SetMinimum(-1);
    gr->SetMaximum(12.0);
    gr->Draw("AP");
    c->Modified(); c->Update(); c->SaveAs("Linearity_test_magnify.png");

    delete gr;
    delete c;

    gStyle->SetOptFit(0); gStyle->SetStatH();

}

void Signal_yield_fit()
{
    // to extract signal yield
    RooRealVar EeclFit("Eecl", "Eecl", 0, 4, "GeV");
//    EeclFit.setBins(EeclBins);
//    RooPlot* Eeclframe = EeclFit.frame(Bins(EeclBins), Title(" "));
    RooPlot* Eeclframe = EeclFit.frame(Title(" "));

    // get data from root files
    const char* MC_dirname_Knunu = "./SIGNAL_analysis/test_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Knunu";
    const char* MC_dirname_Kstarnunu = "./SIGNAL_analysis/test_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Kstarnunu";
    const char* MC_dirname_Xsununu = "./SIGNAL_analysis/test_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Xsnunu";
    const char* MC_dirname_K0nunu = "./SIGNAL_analysis/test_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0nunu";
    const char* MC_dirname_K0starnunu = "./SIGNAL_analysis/test_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0starnunu";
    const char* MC_dirname_Xsdnunu = "./SIGNAL_analysis/test_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02Xsnunu";
    LetsAdd(MC_dirname_Knunu, &Mbc_MC_signal, &Eecl_MC_signal, &weight_MC_signal, &info_MC_signal, Scale_Kplus);
    LetsAdd(MC_dirname_Kstarnunu, &Mbc_MC_signal, &Eecl_MC_signal, &weight_MC_signal, &info_MC_signal, Scale_Kplusstar);
    LetsAdd(MC_dirname_Xsununu, &Mbc_MC_signal, &Eecl_MC_signal, &weight_MC_signal, &info_MC_signal, Scale_Xsu_nonresonant);
    LetsAdd(MC_dirname_K0nunu, &Mbc_MC_signal, &Eecl_MC_signal, &weight_MC_signal, &info_MC_signal, Scale_K0);
    LetsAdd(MC_dirname_K0starnunu, &Mbc_MC_signal, &Eecl_MC_signal, &weight_MC_signal, &info_MC_signal, Scale_K0star);
    LetsAdd(MC_dirname_Xsdnunu, &Mbc_MC_signal, &Eecl_MC_signal, &weight_MC_signal, &info_MC_signal, Scale_Xsd_nonresonant);

    const char* MC_dirname_CHG = "./CHG_analysis/test_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_MIX = "./MIX_analysis/test_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_UUBAR = "./UUBAR_analysis/test_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_DDBAR = "./DDBAR_analysis/test_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_SSBAR = "./SSBAR_analysis/test_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* MC_dirname_CHARM = "./CHARM_analysis/test_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    LetsAdd(MC_dirname_CHG, &Mbc_MC_background, &Eecl_MC_background, &weight_MC_background, &info_MC_background);
    LetsAdd(MC_dirname_MIX, &Mbc_MC_background, &Eecl_MC_background, &weight_MC_background, &info_MC_background);
    LetsAdd(MC_dirname_UUBAR, &Mbc_MC_background, &Eecl_MC_background, &weight_MC_background, &info_MC_background);
    LetsAdd(MC_dirname_DDBAR, &Mbc_MC_background, &Eecl_MC_background, &weight_MC_background, &info_MC_background);
    LetsAdd(MC_dirname_SSBAR, &Mbc_MC_background, &Eecl_MC_background, &weight_MC_background, &info_MC_background);
    LetsAdd(MC_dirname_CHARM, &Mbc_MC_background, &Eecl_MC_background, &weight_MC_background, &info_MC_background);

    const char* DATA_dirname_Knunu = "./SIGNAL_analysis/validation_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Knunu";
    const char* DATA_dirname_Kstarnunu = "./SIGNAL_analysis/validation_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Kstarnunu";
    const char* DATA_dirname_Xsununu = "./SIGNAL_analysis/validation_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B2Xsnunu";
    const char* DATA_dirname_K0nunu = "./SIGNAL_analysis/validation_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0nunu";
    const char* DATA_dirname_K0starnunu = "./SIGNAL_analysis/validation_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02K0starnunu";
    const char* DATA_dirname_Xsdnunu = "./SIGNAL_analysis/validation_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge/B02Xsnunu";
    LetsAdd(DATA_dirname_Knunu, &Mbc_DATA, &Eecl_DATA, &weight_DATA, &info_DATA, Scale_Kplus);
    LetsAdd(DATA_dirname_Kstarnunu, &Mbc_DATA, &Eecl_DATA, &weight_DATA, &info_DATA, Scale_Kplusstar);
    LetsAdd(DATA_dirname_Xsununu, &Mbc_DATA, &Eecl_DATA, &weight_DATA, &info_DATA, Scale_Xsu_nonresonant);
    LetsAdd(DATA_dirname_K0nunu, &Mbc_DATA, &Eecl_DATA, &weight_DATA, &info_DATA, Scale_K0);
    LetsAdd(DATA_dirname_K0starnunu, &Mbc_DATA, &Eecl_DATA, &weight_DATA, &info_DATA, Scale_K0star);
    LetsAdd(DATA_dirname_Xsdnunu, &Mbc_DATA, &Eecl_DATA, &weight_DATA, &info_DATA, Scale_Xsd_nonresonant);

    const char* DATA_dirname_CHG = "./CHG_analysis/validation_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_MIX = "./MIX_analysis/validation_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_UUBAR = "./UUBAR_analysis/validation_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_DDBAR = "./DDBAR_analysis/validation_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_SSBAR = "./SSBAR_analysis/validation_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    const char* DATA_dirname_CHARM = "./CHARM_analysis/validation_v003/final_output_root_after_MVA_Application_after_cut/BCS_only/Merge";
    LetsAdd(DATA_dirname_CHG, &Mbc_DATA, &Eecl_DATA, &weight_DATA, &info_DATA);
    LetsAdd(DATA_dirname_MIX, &Mbc_DATA, &Eecl_DATA, &weight_DATA, &info_DATA);
    LetsAdd(DATA_dirname_UUBAR, &Mbc_DATA, &Eecl_DATA, &weight_DATA, &info_DATA);
    LetsAdd(DATA_dirname_DDBAR, &Mbc_DATA, &Eecl_DATA, &weight_DATA, &info_DATA);
    LetsAdd(DATA_dirname_SSBAR, &Mbc_DATA, &Eecl_DATA, &weight_DATA, &info_DATA);
    LetsAdd(DATA_dirname_CHARM, &Mbc_DATA, &Eecl_DATA, &weight_DATA, &info_DATA);


    // define frame and get ready to make pdfs
//    Eecl_DATA.setBins(EeclBins);
    RooDataSet* d_Eecl = (RooDataSet*)info_DATA.reduce(RooArgSet(Eecl_DATA));
//    RooDataHist binned_data_Eecl("binned Eecl data", "binned Eecl data", EeclFit, *d_Eecl);

//    Eecl_MC_signal.setBins(EeclBins);
//    Eecl_MC_background.setBins(EeclBins);

    // define pdf and extended pdf
    RooRealVar m0("m0", "m0", 0.65);
    RooRealVar sigma("sigma", "sigma", 0.53);
    RooRealVar alpha("alpha", "alpha", -0.85);
    RooRealVar n("n", "n", 60);
    RooCBShape histpdf_Eecl_signal("genpdfs", "genpdfs", EeclFit, m0, sigma, alpha, n);

//    RooDataSet* dataset_Eecl_MC_signal = (RooDataSet*)info_MC_signal.reduce(RooArgSet(Eecl_MC_signal));
//    RooDataHist hist_Eecl_MC_signal("hist_Eecl_MC_signal", "histogram for Eecl of MC signal samples", EeclFit, *dataset_Eecl_MC_signal);
//    RooHistPdf histpdf_Eecl_signal("histpdf_Eecl_signal", "histpdf_Eecl_signal", EeclFit, hist_Eecl_MC_signal, 0);

    RooRealVar nsig("nsig", "n_{sig}", 7.5, -100, 100);
    RooExtendPdf esig("esignal", "extended signal p.d.f", histpdf_Eecl_signal, nsig);

    //    RooRealVar p1("p1","coeff #1", 1.65, -10, 15.0);
    //    RooPolynomial bkg_linear("bkg_linear","bkg_linear", EeclFit, RooArgList(p1));
//        RooRealVar gausmean("gausmean","",2.0, 1.7, 2.3);
//        RooRealVar gauswidth("gauswidth","", 0.547, 0.2, 1.2);
//        RooGaussian histpdf_Eecl_background("bkg_gauss","bkg_gauss",EeclFit,gausmean,gauswidth);
    //    RooRealVar bkgfrac("bkgfrac","fraction in bkg",0.589,0.,1.);   
    //    RooAddPdf histpdf_Eecl_background("int","int",RooArgList(bkg_gauss,bkg_linear),bkgfrac);
//        RooRealVar p1("p1","coeff #1", 153);
//        RooRealVar p2("p2","coeff #2", -48);
//        RooRealVar p3("p3","coeff #3", 2.6);
//        RooRealVar p1("p1","coeff #1", 153, 40.0, 190);
//        RooRealVar p2("p2","coeff #2", -48, -80, -30);
//        RooRealVar p3("p3","coeff #3", 2.6, -4.0, 10.0);

//    RooRealVar p1("p1", "coeff #1", 225);
//       RooRealVar p2("p2","coeff #2", -69);
//    RooRealVar p3("p3", "coeff #3", 3.69);
//        RooRealVar p1("p1","coeff #1", 90, 20, 1000.0);
//    RooRealVar p2("p2", "coeff #2", -10, -70, 10);
//        RooRealVar p3("p3","coeff #3", 1.5, -20, 50);
//     RooPolynomial histpdf_Eecl_background("histpdf_Eecl_background", "histpdf_Eecl_background", EeclFit, RooArgList(p1, p2, p3));

    RooRealVar m0b("m0b", "m0b", 1.6);
    RooRealVar sigmab("sigmab", "sigmab", 0.67);
    RooRealVar alphab("alphab", "alphab", -0.368);
    RooRealVar nb("nb", "nb", 20);
//    RooRealVar m0b("m0b", "m0b", 1.6, 1.5, 1.8);
//    RooRealVar sigmab("sigmab", "sigmab", 1, 0.4, 2.0);
//    RooRealVar alphab("alphab", "alphab", -1, -2.0, 0.0);
//    RooRealVar nb("nb", "nb", 50, 0, 60);
    RooCBShape histpdf_Eecl_background("genpdfb", "genpdfb", EeclFit, m0b, sigmab, alphab, nb);

//        RooRealVar p1("p1","coeff #1", -3.6e-03);
//        RooRealVar p2("p2","coeff #2", -6.37e-01);
//        RooRealVar p3("p3","coeff #3", 2.78e-01);
//       RooRealVar p1("p1", "coeff #1", 6.65e-02, -2.0e-01, 2.0e-01);
//        RooRealVar p2("p2", "coeff #2", -8.67e-01, -9.0e-01, -3.0e-01);
//        RooRealVar p3("p3", "coeff #3", 8.73e-02, -2.5e-01, 8.5e-01);
//        RooChebychev histpdf_Eecl_background("histpdf_Eecl_background", "histpdf_Eecl_background", EeclFit, RooArgList(p1, p2, p3));

//    RooDataSet* dataset_Eecl_MC_background = (RooDataSet*)info_MC_background.reduce(RooArgSet(Eecl_MC_background));
//    RooDataHist hist_Eecl_MC_background("hist_Eecl_MC_background", "histogram for Eecl of MC background samples", EeclFit, *dataset_Eecl_MC_background);
//    RooHistPdf histpdf_Eecl_background("histpdf_Eecl_background", "histpdf_Eecl_background", EeclFit, hist_Eecl_MC_background, 0);

//    RooRealVar p1("p1", "coeff #1", 39.89);
//       RooRealVar p2("p2","coeff #2", -26.03);
//        RooRealVar p1("p1","coeff #1", 40, 0.0, 4000.0);
//    RooRealVar p2("p2", "coeff #2", -10, -1000, -1);
//    RooPolynomial histpdf_Eecl_background("histpdf_Eecl_background", "histpdf_Eecl_background", EeclFit, RooArgList(p1, p2));

    RooRealVar nbkg("nbkg", "number of background events", 450, 300, 600);
    RooExtendPdf ebkg("ebkg", "extended background p.d.f", histpdf_Eecl_background, nbkg);

    RooAddPdf  totalpdf("model", "b+n", RooArgList(ebkg, esig));

    // fit
    RooFitResult* r = totalpdf.fitTo(*d_Eecl, Save(), SumW2Error(false));

    // Draw result
    d_Eecl->plotOn(Eeclframe, DataError(RooAbsData::Poisson));
    totalpdf.plotOn(Eeclframe, LineColor(kRed));
    totalpdf.plotOn(Eeclframe, Components(esig), LineColor(kBlue), LineStyle(kDashed));
    totalpdf.plotOn(Eeclframe, Components(ebkg), LineColor(kViolet), LineStyle(kDashed));

    TCanvas* c = new TCanvas("Eecl", "Eecl", 600, 600);
    gPad->SetLeftMargin(0.15); Eeclframe->GetYaxis()->SetTitleOffset(1.4); Eeclframe->Draw(); c->SaveAs("Eecl_distribution.png");
    delete c;

    
    /* ============== Linearity test ============== */
//    double BKG_total_Num = 0;
//    BKG_total_Num = BKG_total_Num + GetEvtNum(DATA_dirname_CHG);
//    BKG_total_Num = BKG_total_Num + GetEvtNum(DATA_dirname_MIX);
//    BKG_total_Num = BKG_total_Num + GetEvtNum(DATA_dirname_UUBAR);
//    BKG_total_Num = BKG_total_Num + GetEvtNum(DATA_dirname_DDBAR);
//    BKG_total_Num = BKG_total_Num + GetEvtNum(DATA_dirname_SSBAR);
//    BKG_total_Num = BKG_total_Num + GetEvtNum(DATA_dirname_CHARM);
//    LinearityTest(r, EeclFit, BKG_total_Num);

    /* ============== Print profile likelihood  ============== */
//    RooPlot* nllframe = nsig.frame();
//    RooAbsReal* nll = totalpdf.createNLL(*d_Eecl, NumCPU(12), Verbose(false));
//    RooAbsReal* pll_nsig = nll->createProfile(nsig);
//    pll_nsig->plotOn(nllframe); //nllframe->SetMaximum(10.0);

//    TCanvas* c_pll = new TCanvas("pll", "pll", 600, 600);
//    gPad->SetLeftMargin(0.15); nllframe->GetYaxis()->SetTitleOffset(1.4); nllframe->Draw(); c_pll->SaveAs("PLL.png");
//    delete c_pll;

    /* ============== toy MC study ============== */
    RooRealVar  Eecl_TOY("Eecl", "Eecl_TOY", 0, 4);
//    Eecl_TOY.setBins(EeclBins);

    RooMCStudy* mcstudy = new RooMCStudy(totalpdf, Eecl_TOY, Binned(kTRUE), Silence(), Extended(),FitOptions(Save(kTRUE), PrintEvalErrors(0)));
    mcstudy->generateAndFit(1000);

    // Make plots of the distributions of mean, the error on mean and the pull of mean
    RooPlot* frame1 = mcstudy->plotParam(nsig, Bins(40));
    RooPlot* frame2 = mcstudy->plotError(nsig, Bins(40));
    RooPlot* frame3 = mcstudy->plotPull(nsig, Bins(40), FitGauss(kTRUE));

    // Draw all plots on a canvas
    gStyle->SetOptStat(0);
    TCanvas* cf = new TCanvas("rf801_mcstudy", "rf801_mcstudy", 1200, 400);
    cf->Divide(3, 1);
    cf->cd(1); gPad->SetLeftMargin(0.15); frame1->GetYaxis()->SetTitleOffset(1.4); frame1->Draw();
    cf->cd(2); gPad->SetLeftMargin(0.15); frame2->GetYaxis()->SetTitleOffset(1.4); frame2->Draw();
    cf->cd(3); gPad->SetLeftMargin(0.15); frame3->GetYaxis()->SetTitleOffset(1.4); frame3->Draw();
    cf->SaveAs("ToyStudy.png");

    /* ============== Uncertainty Calculation ============== */
    LetsCalculateUncertainties(MC_dirname_Knunu, Scale_Kplus);
    LetsCalculateUncertainties(MC_dirname_Kstarnunu, Scale_Kplusstar);
    LetsCalculateUncertainties(MC_dirname_Xsununu, Scale_Xsu_nonresonant);
    LetsCalculateUncertainties(MC_dirname_K0nunu, Scale_K0);
    LetsCalculateUncertainties(MC_dirname_K0starnunu, Scale_K0star);
    LetsCalculateUncertainties(MC_dirname_Xsdnunu, Scale_Xsd_nonresonant);
    PrintUncertainties();
}
