// last update: 2023-07-14
// for Belle2 data

#include <stdio.h>
#include <string>
#include <vector>

#include "TFile.h"

#include "constants.h"
#include "base.h"
#include "ObtainWeight.h"
#include "Loader.h"

/*
when you add new variables:
revise # define N_Needed_info ...
revise typedef struct dataLoader::ConvertIntoSeparateDataFile(std::string output_name, double flag = 0)
revise void Loader::GetData(TFile* input_file)
check MakeShiftDoubleToInt if it affects the code
revise void Loader::PrintRootFile(std::string output_name)
revise void Loader::PrintSeparateRootFile(std::string output_name)
revise void Loader::ConvertIntoSeparateDataFile(std::string output_name, double flag = 0)
*/

/* ====================================== */

int main(int argc, char* argv[]) {
    /*
    * argv[1]: dirname
    * argv[2]: outputname
    * argv[3]: output path
    * argv[4]: sample type: {SIGNAL|CHG|MIX|UUBAR|DDBAR|SSBAR|CHARM}
    * argv[5]: MC version: {data|MC15ri|MC15rd}
    * argv[6]: sample category: {test|train|validation}
    */

    std::vector<std::string> names;

    load_files(argv[1], &names);

    Loader loader;

    for (unsigned int i = 0; i < names.size(); i++) {
        loader.initialize();

        TFile* input_file = new TFile((argv[1] + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());
        loader.GetData(input_file);
        if (loader.event_info_is_valid() == false) { printf("error!\n"); return 0; }

        std::string::size_type const p(names.at(i).find_last_of('.'));
        std::string file_without_extension = names.at(i).substr(0, p);

        loader.Cut(Loader::Upsilon, 56, Loader::smaller_than, 1.5); // nloose track
        loader.Cut(Loader::Upsilon, 6, Loader::smaller_than, 0.5); // N pi0
        loader.Cut(Loader::Upsilon, 5, Loader::smaller_than, 0.5); // N KS0
        loader.Cut(Loader::Bsig, 64, Loader::larger_than, -0.5); // Bsig vertex fit
        loader.Cut(Loader::Btag, 6, Loader::larger_than, -0.5); // Btag vertex fit
        loader.PrintInformation(std::string("========== inital =========="), names.at(i), argv[4], argv[5], argv[6], true);

        loader.PrintSeparateRootFile("./" + file_without_extension + std::string("_before_Mbcsig_") + std::string(argv[2]) + ".root");
        //loader.DrawTH1F("Bsig_Mbc", "Bsig_Mbc;Bsig_Mbc [GeV];candidates", 100, 5.24, 5.3, Loader::Bsig, 738);
        loader.Cut(Loader::Bsig, 738, Loader::larger_than, 5.2);
        loader.PrintInformation(std::string("========== sig Mbc > 5.2 =========="), names.at(i), argv[4], argv[5], argv[6], true);
        //loader.DrawTH2F("MbcVSdeltaE_after_Mbc_strict_cut", ";Mbc of B_{tag} [GeV];#DeltaE of B_{tag} [GeV]", 100, 5.24, 5.3, 100, -0.2, 0.2, Loader::Btag, 1, Loader::Btag, 2);

        loader.PrintSeparateRootFile("./" + file_without_extension + std::string("_before_deltaE_sig_") + std::string(argv[2]) + ".root");
        //loader.DrawTH1F("Bsig_deltaE", "Bsig_deltaE;Bsig_deltaE [GeV];candidates", 100, -0.3, 0.3, Loader::Bsig, 739);
        loader.Cut(Loader::Bsig, 739, Loader::larger_than, -0.05);
        loader.Cut(Loader::Bsig, 739, Loader::smaller_than, 0.05);
        loader.PrintInformation(std::string("========== sig abs(deltaE) < 0.05 =========="), names.at(i), argv[4], argv[5], argv[6], true);
        //loader.DrawTH1F("Bsig_Mbc", "Bsig_Mbc;Bsig_Mbc [GeV];candidates", 100, 5.24, 5.3, Loader::Bsig, 738);
        //loader.DrawTH2F("MbcVSdeltaE_after_deltaE_strict_cut", ";Mbc of B_{tag} [GeV];#DeltaE of B_{tag} [GeV]", 100, 5.24, 5.3, 100, -0.2, 0.2, Loader::Btag, 1, Loader::Btag, 2);

        //loader.DrawTH1F("Btag_Mbc", "Btag_Mbc;Btag_Mbc [GeV];candidates", 100, 5.24, 5.3, Loader::Btag, 1);
//        loader.Cut(Loader::Btag, 1, Loader::larger_than, 5.27);
//        loader.PrintInformation(std::string("========== Mbc > 5.27 =========="), names.at(i), argv[4], argv[5], argv[6], true);
        //loader.DrawTH2F("MbcVSdeltaE_after_Mbc_strict_cut", ";Mbc of B_{tag} [GeV];#DeltaE of B_{tag} [GeV]", 100, 5.24, 5.3, 100, -0.2, 0.2, Loader::Btag, 1, Loader::Btag, 2);

        //loader.DrawTH1F("Btag_deltaE", "Btag_deltaE;Btag_deltaE [GeV];candidates", 100, -0.3, 0.3, Loader::Btag, 2);
        loader.PrintSeparateRootFile("./" + file_without_extension + std::string("_before_deltaE_tag_") + std::string(argv[2]) + ".root");
        loader.Cut(Loader::Btag, 2, Loader::larger_than, -0.2);
        loader.Cut(Loader::Btag, 2, Loader::smaller_than, 0.2);
        loader.PrintInformation(std::string("========== abs(deltaE) < 0.2 =========="), names.at(i), argv[4], argv[5], argv[6], true);
        //loader.DrawTH2F("MbcVSdeltaE_after_deltaE_strict_cut", ";Mbc of B_{tag} [GeV];#DeltaE of B_{tag} [GeV]", 100, 5.24, 5.3, 100, -0.2, 0.2, Loader::Btag, 1, Loader::Btag, 2);

        //loader.DrawTH1F("missing_momentum_theta_after_npi0_cut", "#theta_{missing};#theta_{missing} [rad];evt", 100, 0, 3.2, Loader::Upsilon, 7);
        loader.PrintSeparateRootFile("./" + file_without_extension + std::string("_before_Eecl_") + std::string(argv[2]) + ".root");
        loader.Cut(Loader::Upsilon, 69, Loader::smaller_than, 1.3);
        loader.PrintInformation(std::string("========== Eecl_v200 < 1.3 =========="), names.at(i), argv[4], argv[5], argv[6], true);
/*
        //loader.DrawTH1F("missing_momentum_theta_after_npi0_cut", "#theta_{missing};#theta_{missing} [rad];evt", 100, 0, 3.2, Loader::Upsilon, 7);
        loader.Cut(Loader::Upsilon, 7, Loader::smaller_than, 2.618);
        loader.Cut(Loader::Upsilon, 7, Loader::larger_than, 0.297);
        loader.PrintInformation(std::string("========== 0.297 < missing momentum theta < 2.618 =========="), names.at(i), argv[4], argv[5], argv[6], true);

        //loader.DrawTH1F("momentum_Bsig_after_missing_theta_cut", "momentum of B_{sig} at CMS;momentum [GeV];evt", 100, 0, 3.2, Loader::Bsig, 4);
        loader.Cut(Loader::Bsig, 4, Loader::smaller_than, 2.96);
        loader.Cut(Loader::Bsig, 4, Loader::larger_than, 0.5);
        loader.PrintInformation(std::string("========== 0.5 < momentum of signal side < 2.96 =========="), names.at(i), argv[4], argv[5], argv[6], true);

        //loader.DrawTH1F("Bsig_M_Xs", "mass of X_{s};M_{Xs} [GeV];evt", 100, 0, 3.5, Loader::Bsig, 6);
        loader.DvetoFor(Loader::Bsig, 6, 1.84, 1.89);
        loader.PrintInformation(std::string("========== D veto =========="), names.at(i), argv[4], argv[5], argv[6], true);
*/

        loader.PrintSeparateRootFile("./" + file_without_extension + std::string("_before_sphericity_") + std::string(argv[2]) + ".root");
        loader.Cut(Loader::Upsilon, 32, Loader::larger_than, 0.2);
        loader.PrintInformation(std::string("========== sphericity > 0.2 =========="), names.at(i), argv[4], argv[5], argv[6], true);

        loader.PrintSeparateRootFile("./" + file_without_extension + std::string("_before_R2_") + std::string(argv[2]) + ".root");
        loader.Cut(Loader::Bsig, 14, Loader::smaller_than, 0.4);
        loader.PrintInformation(std::string("========== Bsig_R2 < 0.4 =========="), names.at(i), argv[4], argv[5], argv[6], true);

        loader.PrintSeparateRootFile("./" + file_without_extension + std::string("_before_cosTBTO_") + std::string(argv[2]) + ".root");
        loader.Cut(Loader::Bsig, 17, Loader::smaller_than, 0.9);
        loader.PrintInformation(std::string("========== Bsig_cosTBTO < 0.9 =========="), names.at(i), argv[4], argv[5], argv[6], true);

        loader.Cut(Loader::Bsig, 740, Loader::larger_than, 1.0);
        loader.PrintInformation(std::string("========== sig MXs > 1.0 =========="), names.at(i), argv[4], argv[5], argv[6], true);

        loader.BCS(Loader::Btag, 5, Loader::Highest);
        if (loader.IsBCSValid() == false) {
            printf("ERROR! it is not valid\n");
            exit(1);
        }
        loader.PrintInformation(std::string("========== BCS =========="), names.at(i), argv[4], argv[5], argv[6], true);
        //loader.DrawTH1F("Mppbar", "Mppbar;Mppbar [GeV];candidates", 10, 1.8, 3.8, Loader::Upsilon, 164);

//        loader.Cut(Loader::Upsilon, 164, Loader::larger_than, 1.8);
//        loader.Cut(Loader::Upsilon, 164, Loader::smaller_than, 2.46666666667);

        loader.PrintSeparateRootFile("./" + file_without_extension + std::string("_final_") + std::string(argv[2]) + ".root");
    }
    loader.End();

    return 0;
}
