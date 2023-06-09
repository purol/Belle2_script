# define N_decay 38 // five decay mode + others

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

enum DecayMode { // reco level
    B2Kc = 0,
    B2KcPi0,
    B2Ks0Pic,
    B2KcPicPic,
    B2Ks0PicPi0,
    B2KcPicPicPi0,
    B2Ks0PicPicPic,
    B2KcPicPicPicPic,
    B2Ks0PicPicPicPi0,
    B2KcPi0Pi0,
    B2Ks0PicPi0Pi0,
    B2KcPicPicPi0Pi0,
    B2KcKcKc,
    B2KcKcKs0Pic,
    B2KcKcKcPi0,
    B02Ks0,
    B02KcPic,
    B02Ks0Pi0,
    B02KcPicPi0,
    B02Ks0PicPic,
    B02KcPicPicPic,
    B02Ks0PicPicPi0,
    B02KcPicPicPicPi0,
    B02Ks0PicPicPicPic,
    B02Ks0Pi0Pi0,
    B02KcPicPi0Pi0,
    B02Ks0PicPicPi0Pi0,
    B02KcKcKs0,
    B02KcKcKcPic,
    B02KcKcKs0Pi0,
    MAX_NUM_DECAYMODE
};
enum DecayModeMC { // MC level
    Xsu2Kc_MC = 0,
    Xsu2Kcstar2KcPi0_MC,
    Xsu2Kcstar2K0Pic_MC,
    Xsu2KcPi0_MC,
    Xsu2K0Pic_MC,
    Xsu2KcPicPic_MC,
    Xsu2K0PicPi0_MC,
    Xsu2KcPicPicPi0_MC,
    Xsu2K0PicPicPic_MC,
    Xsu2KcPicPicPicPic_MC,
    Xsu2K0PicPicPicPi0_MC,
    Xsu2KcPi0Pi0_MC,
    Xsu2K0PicPi0Pi0_MC,
    Xsu2KcPicPicPi0Pi0_MC,
    Xsu2KcKcKc_MC,
    Xsu2KcKcK0Pic_MC,
    Xsu2KcKcKcPi0_MC,
    Xsd2K0_MC,
    Xsd2K0star2KcPic_MC,
    Xsd2K0star2K0Pi0_MC,
    Xsd2KcPic_MC,
    Xsd2K0Pi0_MC,
    Xsd2KcPicPi0_MC,
    Xsd2K0PicPic_MC,
    Xsd2KcPicPicPic_MC,
    Xsd2K0PicPicPi0_MC,
    Xsd2KcPicPicPicPi0_MC,
    Xsd2K0PicPicPicPic_MC,
    Xsd2K0Pi0Pi0_MC,
    Xsd2KcPicPi0Pi0_MC,
    Xsd2K0PicPicPi0Pi0_MC,
    Xsd2KcKcK0_MC,
    Xsd2KcKcKcPic_MC,
    Xsd2KcKcK0Pi0_MC,
    other,
    MAX_NUM_DECAYMODE_MC
};
enum SimpleDecayModeMC { // MC level
    Xs2KKstar = 0,
    Xs2Kpi,
    Xs2Kpipi,
    simple_other,
    MAX_NUM_SIMPLE_DECAYMODE_MC
};

bool TrueIfDecayModeMatch(double Upsilon_decayID, double Bsig_decayID, DecayMode decaymode) {
    switch (decaymode) {
    case DecayMode::B2Kc:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > -0.5 && Bsig_decayID < 0.5) return true;
        return false;
        break;
    case DecayMode::B2KcPi0:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 0.5 && Bsig_decayID < 1.5) return true;
        return false;
        break;
    case DecayMode::B2Ks0Pic:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 1.5 && Bsig_decayID < 2.5) return true;
        return false;
        break;
    case DecayMode::B2KcPicPic:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 2.5 && Bsig_decayID < 3.5) return true;
        return false;
        break;
    case DecayMode::B2Ks0PicPi0:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 3.5 && Bsig_decayID < 4.5) return true;
        return false;
        break;
    case DecayMode::B2KcPicPicPi0:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 4.5 && Bsig_decayID < 5.5) return true;
        return false;
        break;
    case DecayMode::B2Ks0PicPicPic:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 5.5 && Bsig_decayID < 6.5) return true;
        return false;
        break;
    case DecayMode::B2KcPicPicPicPic:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 6.5 && Bsig_decayID < 7.5) return true;
        return false;
        break;
    case DecayMode::B2Ks0PicPicPicPi0:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 7.5 && Bsig_decayID < 8.5) return true;
        return false;
        break;
    case DecayMode::B2KcPi0Pi0:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 8.5 && Bsig_decayID < 9.5) return true;
        return false;
        break;
    case DecayMode::B2Ks0PicPi0Pi0:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 9.5 && Bsig_decayID < 10.5) return true;
        return false;
        break;
    case DecayMode::B2KcPicPicPi0Pi0:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 10.5 && Bsig_decayID < 11.5) return true;
        return false;
        break;
    case DecayMode::B2KcKcKc:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 11.5 && Bsig_decayID < 12.5) return true;
        return false;
        break;
    case DecayMode::B2KcKcKs0Pic:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 12.5 && Bsig_decayID < 13.5) return true;
        return false;
        break;
    case DecayMode::B2KcKcKcPi0:
        if (Upsilon_decayID > -0.5 && Upsilon_decayID < 0.5 && Bsig_decayID > 13.5 && Bsig_decayID < 14.5) return true;
        return false;
        break;
    case DecayMode::B02Ks0:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > -0.5 && Bsig_decayID < 0.5) return true;
        return false;
        break;
    case DecayMode::B02KcPic:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 0.5 && Bsig_decayID < 1.5) return true;
        return false;
        break;
    case DecayMode::B02Ks0Pi0:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 1.5 && Bsig_decayID < 2.5) return true;
        return false;
        break;
    case DecayMode::B02KcPicPi0:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 2.5 && Bsig_decayID < 3.5) return true;
        return false;
        break;
    case DecayMode::B02Ks0PicPic:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 3.5 && Bsig_decayID < 4.5) return true;
        return false;
        break;
    case DecayMode::B02KcPicPicPic:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 4.5 && Bsig_decayID < 5.5) return true;
        return false;
        break;
    case DecayMode::B02Ks0PicPicPi0:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 5.5 && Bsig_decayID < 6.5) return true;
        return false;
        break;
    case DecayMode::B02KcPicPicPicPi0:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 6.5 && Bsig_decayID < 7.5) return true;
        return false;
        break;
    case DecayMode::B02Ks0PicPicPicPic:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 7.5 && Bsig_decayID < 8.5) return true;
        return false;
        break;
    case DecayMode::B02Ks0Pi0Pi0:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 8.5 && Bsig_decayID < 9.5) return true;
        return false;
        break;
    case DecayMode::B02KcPicPi0Pi0:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 9.5 && Bsig_decayID < 10.5) return true;
        return false;
        break;
    case DecayMode::B02Ks0PicPicPi0Pi0:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 10.5 && Bsig_decayID < 11.5) return true;
        return false;
        break;
    case DecayMode::B02KcKcKs0:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 11.5 && Bsig_decayID < 12.5) return true;
        return false;
        break;
    case DecayMode::B02KcKcKcPic:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 12.5 && Bsig_decayID < 13.5) return true;
        return false;
        break;
    case DecayMode::B02KcKcKs0Pi0:
        if (Upsilon_decayID > 0.5 && Upsilon_decayID < 1.5 && Bsig_decayID > 13.5 && Bsig_decayID < 14.5) return true;
        return false;
        break;
    default:
        printf("ERROR! 036\n");
        exit(1);
        break;
    }

    printf("ERROR! 037\n");
    exit(1);
    return false;
}

bool TrueIfDecayModeMatch_MC(int Decay[N_decay], DecayModeMC decaymodeMC) {

    switch (decaymodeMC) {
    case DecayModeMC::Xsu2Kc_MC:
        if (Decay[0] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2Kcstar2KcPi0_MC:
        if (Decay[2] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2Kcstar2K0Pic_MC:
        if (Decay[1] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2KcPi0_MC:
        if (Decay[5] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2K0Pic_MC:
        if (Decay[6] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2KcPicPic_MC:
        if (Decay[7] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2K0PicPi0_MC:
        if (Decay[8] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2KcPicPicPi0_MC:
        if (Decay[9] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2K0PicPicPic_MC:
        if (Decay[10] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2KcPicPicPicPic_MC:
        if (Decay[11] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2K0PicPicPicPi0_MC:
        if (Decay[12] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2KcPi0Pi0_MC:
        if (Decay[13] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2K0PicPi0Pi0_MC:
        if (Decay[14] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2KcPicPicPi0Pi0_MC:
        if (Decay[15] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2KcKcKc_MC:
        if (Decay[16] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2KcKcK0Pic_MC:
        if (Decay[17] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsu2KcKcKcPi0_MC:
        if (Decay[18] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2K0_MC:
        if (Decay[19] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2K0star2KcPic_MC:
        if (Decay[20] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2K0star2K0Pi0_MC:
        if (Decay[21] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2KcPic_MC:
        if (Decay[24] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2K0Pi0_MC:
        if (Decay[25] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2KcPicPi0_MC:
        if (Decay[26] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2K0PicPic_MC:
        if (Decay[27] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2KcPicPicPic_MC:
        if (Decay[28] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2K0PicPicPi0_MC:
        if (Decay[29] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2KcPicPicPicPi0_MC:
        if (Decay[30] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2K0PicPicPicPic_MC:
        if (Decay[31] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2K0Pi0Pi0_MC:
        if (Decay[32] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2KcPicPi0Pi0_MC:
        if (Decay[33] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2K0PicPicPi0Pi0_MC:
        if (Decay[34] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2KcKcK0_MC:
        if (Decay[35] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2KcKcKcPic_MC:
        if (Decay[36] > 0) return true;
        return false;
        break;
    case DecayModeMC::Xsd2KcKcK0Pi0_MC:
        if (Decay[37] > 0) return true;
        return false;
        break;
    case DecayModeMC::other:
        return true;
        break;
    default:
        printf("ERROR! Input value of TrueIfDecayModeMatch_MC is not appropriate\n");
        exit(1);
        break;
    }

    printf("ERROR! Input value of TrueIfDecayModeMatch_MC is not appropriate\n");
    exit(1);
    return false;
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

void LetsFill(const char* dirname, double OneDConfusion[DecayModeMC::MAX_NUM_DECAYMODE_MC], double weight = 1.0) {

    std::vector<string> names;
    load_files(dirname, &names);

    int Decay[N_decay] = { 0 }; // MC level info

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_Jpsi = (TTree*)input_file->Get("Jpsi");

        tree_Jpsi->SetBranchAddress("nParticlesInList__boB__pl__clKcharge_total__bc", &Decay[0]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &Decay[1]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &Decay[2]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCcomb__bc", &Decay[3]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch1__bc", &Decay[4]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch2__bc", &Decay[5]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch3__bc", &Decay[6]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch4__bc", &Decay[7]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch5__bc", &Decay[8]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch6__bc", &Decay[9]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch7__bc", &Decay[10]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch8__bc", &Decay[11]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch9__bc", &Decay[12]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch10__bc", &Decay[13]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch11__bc", &Decay[14]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch12__bc", &Decay[15]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch13__bc", &Decay[16]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch14__bc", &Decay[17]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch15__bc", &Decay[18]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boB0__clKneutral_total__bc", &Decay[19]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &Decay[20]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &Decay[21]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCcomb__bc", &Decay[22]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch16__bc", &Decay[23]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch17__bc", &Decay[24]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch18__bc", &Decay[25]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch19__bc", &Decay[26]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch20__bc", &Decay[27]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch21__bc", &Decay[28]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch22__bc", &Decay[29]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch23__bc", &Decay[30]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch24__bc", &Decay[31]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch25__bc", &Decay[32]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch26__bc", &Decay[33]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch27__bc", &Decay[34]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch28__bc", &Decay[35]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch29__bc", &Decay[36]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch30__bc", &Decay[37]);

        printf("%lld entries...\n", tree_Jpsi->GetEntries());
        for (unsigned int j = 0; j < tree_Jpsi->GetEntries(); j++) { // Fill
            tree_Jpsi->GetEntry(j);

            int decaymodeid_MC = -1;
            for (int i = 0; i < DecayModeMC::MAX_NUM_DECAYMODE_MC; i++) { // find MC decay mode
                if (TrueIfDecayModeMatch_MC(Decay, static_cast<DecayModeMC>(i))) {
                    decaymodeid_MC = i;
                    break;
                }
            }
            if (decaymodeid_MC == DecayModeMC::MAX_NUM_DECAYMODE_MC) {
                printf("ERROR! MC decay id cannot be found\n");
                exit(1);
            }

            OneDConfusion[decaymodeid_MC] = OneDConfusion[decaymodeid_MC] + weight;

        }
        input_file->Close();

    }

}

void LetsFill(const char* dirname, double OneDConfusion[DecayModeMC::MAX_NUM_DECAYMODE_MC], const char* included_string, double weight = 1.0) {

    std::vector<string> names;
    load_files(dirname, &names, included_string);

    int Decay[N_decay] = { 0 }; // MC level info

    for (unsigned int i = 0; i < names.size(); i++) {

        TFile* input_file = new TFile((dirname + std::string("/") + names.at(i)).c_str(), "read");
        printf("%s (%d/%zu)\n", ("Read " + names.at(i) + "... ").c_str(), i, names.size());

        TTree* tree_Jpsi = (TTree*)input_file->Get("Jpsi");

        tree_Jpsi->SetBranchAddress("nParticlesInList__boB__pl__clKcharge_total__bc", &Decay[0]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch1_total__bc", &Decay[1]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boB__pl__clKstarcharge_ch2_total__bc", &Decay[2]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCcomb__bc", &Decay[3]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch1__bc", &Decay[4]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch2__bc", &Decay[5]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch3__bc", &Decay[6]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch4__bc", &Decay[7]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch5__bc", &Decay[8]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch6__bc", &Decay[9]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch7__bc", &Decay[10]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch8__bc", &Decay[11]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch9__bc", &Decay[12]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch10__bc", &Decay[13]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch11__bc", &Decay[14]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch12__bc", &Decay[15]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch13__bc", &Decay[16]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch14__bc", &Decay[17]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsu__clMCch15__bc", &Decay[18]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boB0__clKneutral_total__bc", &Decay[19]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch1_total__bc", &Decay[20]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boB0__clKstarneutral_ch2_total__bc", &Decay[21]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCcomb__bc", &Decay[22]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch16__bc", &Decay[23]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch17__bc", &Decay[24]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch18__bc", &Decay[25]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch19__bc", &Decay[26]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch20__bc", &Decay[27]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch21__bc", &Decay[28]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch22__bc", &Decay[29]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch23__bc", &Decay[30]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch24__bc", &Decay[31]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch25__bc", &Decay[32]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch26__bc", &Decay[33]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch27__bc", &Decay[34]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch28__bc", &Decay[35]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch29__bc", &Decay[36]);
        tree_Jpsi->SetBranchAddress("nParticlesInList__boXsd__clMCch30__bc", &Decay[37]);

        printf("%lld entries...\n", tree_Jpsi->GetEntries());
        for (unsigned int j = 0; j < tree_Jpsi->GetEntries(); j++) { // Fill
            tree_Jpsi->GetEntry(j);

            int decaymodeid_MC = -1;
            for (int i = 0; i < DecayModeMC::MAX_NUM_DECAYMODE_MC; i++) { // find MC decay mode
                if (TrueIfDecayModeMatch_MC(Decay, static_cast<DecayModeMC>(i))) {
                    decaymodeid_MC = i;
                    break;
                }
            }
            if (decaymodeid_MC == DecayModeMC::MAX_NUM_DECAYMODE_MC) {
                printf("ERROR! MC decay id cannot be found\n");
                exit(1);
            }

            OneDConfusion[decaymodeid_MC] = OneDConfusion[decaymodeid_MC] + weight;

        }
        input_file->Close();

    }

}

void DecayMatrixToSimpleDecayMatrix(double Confusion[DecayMode::MAX_NUM_DECAYMODE][DecayModeMC::MAX_NUM_DECAYMODE_MC], double SimpleConfusion[DecayMode::MAX_NUM_DECAYMODE][SimpleDecayModeMC::MAX_NUM_SIMPLE_DECAYMODE_MC]) {
    for (int i = 0; i < DecayMode::MAX_NUM_DECAYMODE; i++) {
        SimpleConfusion[i][SimpleDecayModeMC::Xs2KKstar]
            = Confusion[i][DecayModeMC::Xsu2Kc_MC]
            + Confusion[i][DecayModeMC::Xsu2Kcstar2KcPi0_MC]
            + Confusion[i][DecayModeMC::Xsu2Kcstar2K0Pic_MC]
            + Confusion[i][DecayModeMC::Xsd2K0_MC]
            + Confusion[i][DecayModeMC::Xsd2K0star2KcPic_MC]
            + Confusion[i][DecayModeMC::Xsd2K0star2K0Pi0_MC];

        SimpleConfusion[i][SimpleDecayModeMC::Xs2Kpi]
            = Confusion[i][DecayModeMC::Xsu2KcPi0_MC]
            + Confusion[i][DecayModeMC::Xsu2K0Pic_MC]
            + Confusion[i][DecayModeMC::Xsd2KcPic_MC]
            + Confusion[i][DecayModeMC::Xsd2K0Pi0_MC];

        SimpleConfusion[i][SimpleDecayModeMC::Xs2Kpipi]
            = Confusion[i][DecayModeMC::Xsu2KcPicPic_MC]
            + Confusion[i][DecayModeMC::Xsu2K0PicPi0_MC]
            + Confusion[i][DecayModeMC::Xsu2KcPi0Pi0_MC]
            + Confusion[i][DecayModeMC::Xsd2KcPicPi0_MC]
            + Confusion[i][DecayModeMC::Xsd2K0PicPic_MC]
            + Confusion[i][DecayModeMC::Xsd2K0Pi0Pi0_MC];

        SimpleConfusion[i][SimpleDecayModeMC::simple_other]
            = Confusion[i][DecayModeMC::Xsu2KcPicPicPi0_MC]
            + Confusion[i][DecayModeMC::Xsu2K0PicPicPic_MC]
            + Confusion[i][DecayModeMC::Xsu2KcPicPicPicPic_MC]
            + Confusion[i][DecayModeMC::Xsu2K0PicPicPicPi0_MC]
            + Confusion[i][DecayModeMC::Xsu2K0PicPi0Pi0_MC]
            + Confusion[i][DecayModeMC::Xsu2KcPicPicPi0Pi0_MC]
            + Confusion[i][DecayModeMC::Xsu2KcKcKc_MC]
            + Confusion[i][DecayModeMC::Xsu2KcKcK0Pic_MC]
            + Confusion[i][DecayModeMC::Xsu2KcKcKcPi0_MC]
            + Confusion[i][DecayModeMC::Xsd2KcPicPicPic_MC]
            + Confusion[i][DecayModeMC::Xsd2K0PicPicPi0_MC]
            + Confusion[i][DecayModeMC::Xsd2KcPicPicPicPi0_MC]
            + Confusion[i][DecayModeMC::Xsd2K0PicPicPicPic_MC]
            + Confusion[i][DecayModeMC::Xsd2KcPicPi0Pi0_MC]
            + Confusion[i][DecayModeMC::Xsd2K0PicPicPi0Pi0_MC]
            + Confusion[i][DecayModeMC::Xsd2KcKcK0_MC]
            + Confusion[i][DecayModeMC::Xsd2KcKcKcPic_MC]
            + Confusion[i][DecayModeMC::Xsd2KcKcK0Pi0_MC]
            + Confusion[i][DecayModeMC::other]
    }
}

void ReadJpsiFiles_r_sp(){

    double OneDEvt[DecayModeMC::MAX_NUM_DECAYMODE_MC] = { 0.0 }; // [reco]
    double SimpleEvt[DecayMode::MAX_NUM_DECAYMODE][SimpleDecayModeMC::MAX_NUM_SIMPLE_DECAYMODE_MC] = { 0.0 }; // [reco][MC truth]

    std::vector<string> names;
    const char* dirname = "/home/jwpark/storage/DecayInfo/small";

    LetsFill(dirname, OneDEvt, "CHG", ((N_BB_LS1 * (BR_BpBp / (BR_BpBp + BR_B0B0))) / (2.8 * N_BpBp_1invab))); // total 2.8/ab
    LetsFill(dirname, OneDEvt, "MIX", ((N_BB_LS1 * (BR_B0B0 / (BR_BpBp + BR_B0B0))) / (2.8 * N_B0B0_1invab))); // total 2.8/ab

    DecayMatrixToSimpleDecayMatrix(OneDEvt, SimpleEvt);

    printf("--------------- Number of events ---------------\n");
    for (int i = 0; i < DecayModeMC::MAX_NUM_DECAYMODE_MC; i++) {
        printf("%f ", OneDEvt[i]);
        printf("\n");
    }
    printf("--------------- Number of events ---------------\n");
    printf("--------------- Number of simplified events ---------------\n");
    for (int i = 0; i < SimpleDecayModeMC::MAX_NUM_SIMPLE_DECAYMODE_MC; i++) {
        printf("%f ", SimpleEvt[i]);
        printf("\n");
    }
    printf("--------------- Number of simplified events ---------------\n");
}
