#ifndef CONSTANTS_H
#define CONSTANTS_H

# define N_Needed_info 37
# define N_Upsilon_info 189
# define N_Bsig_info 869
# define N_Btag_info 12
# define N_decay 38 // five decay mode + others
# define N_decay_nparticles 5 // # of nu_e, B->Xs nu_e nu_e_bar, B0->Xs nu_e nu_e_bar, B+-, B0
# define N_decay_syst_ff 7 // helicity angle + q2
# define N_PID_syst 73
# define N_fakeE_syst 37
# define N_fakeMU_syst 49
# define N_pi0_syst 8
# define N_pi0_syst_MC15ri 56
# define N_pi0_syst_MC15rd 49
# define index_q2 0
# define index_MXs_Bc 3
# define index_MXs_B0 4

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

// f+-/f0 [https://arxiv.org/abs/2411.18639]
# define fpm_f0 1.052
# define f00 0.4861
# define fpm 0.5113

// https://confluence.desy.de/pages/viewpage.action?spaceKey=BI&title=Conference+readiness
# define N_Upsilon_LS1 387100000.0 // NUpsilon(4S) = (387.1 +/- 5.6) x 10^6
# define N_Bp (2 * fpm * N_Upsilon_LS1)
# define N_B0 (2 * f00 * N_Upsilon_LS1)

# define N_Kplus_nunubar_LS1 (N_Bp * BR_Kplus_nunubar)
# define N_Kplusstar_nunubar_LS1 (N_Bp * BR_Kplusstar_nunubar)
# define N_Xsu_nonresonant_nunubar_LS1 (N_Bp * BR_Xsu_nonresonant_nunubar)
# define N_K0_nunubar_LS1 (N_B0 * BR_K0_nunubar)
# define N_K0star_nunubar_LS1 (N_B0 * BR_K0star_nunubar)
# define N_Xsd_nunubar_LS1 (N_B0 * BR_Xsd_nonresonant_nunubar)

// SIGNAL MC sample number before skimming
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

// SIGNAL MC sample number before skimming for MC15rd
# define N_Kplus_validation_MC15rd 3962022.0
# define N_K0_validation_MC15rd 3962022.0
# define N_Kplusstar_validation_MC15rd 3962022.0
# define N_K0star_validation_MC15rd 3962022.0
# define N_Xsu_nonresonant_validation_MC15rd 15846594.0
# define N_Xsd_nonresonant_validation_MC15rd 15846594.0

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

// scale factor for SIGNAL MC sample until LS1 for MC15rd
# define Scale_Kplus_validation_MC15rd (N_Kplus_nunubar_LS1/N_Kplus_validation_MC15rd)
# define Scale_Kplusstar_validation_MC15rd (N_Kplusstar_nunubar_LS1/N_Kplusstar_validation_MC15rd)
# define Scale_Xsu_nonresonant_validation_MC15rd (N_Xsu_nonresonant_nunubar_LS1/N_Xsu_nonresonant_validation_MC15rd)
# define Scale_K0_validation_MC15rd (N_K0_nunubar_LS1/N_K0_validation_MC15rd)
# define Scale_K0star_validation_MC15rd (N_K0star_nunubar_LS1/N_K0star_validation_MC15rd)
# define Scale_Xsd_nonresonant_validation_MC15rd (N_Xsd_nunubar_LS1/N_Xsd_nonresonant_validation_MC15rd)

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

// untile LS1, integrated luminosity is 365.37/fb
// https://arxiv.org/abs/2407.00965
# define Scale_CHG_train (N_Bp / (2 * 2.8 * N_BpBp_1invab * (N_CHG_train / (N_CHG_train + N_CHG_test)) ))
# define Scale_MIX_train (N_B0 / (2 * 2.8 * N_B0B0_1invab * (N_MIX_train / (N_MIX_train + N_MIX_test)) ))
# define Scale_UUBAR_train (0.36537/((N_UUBAR_train/(N_UUBAR_train + N_UUBAR_test))*1.0))
# define Scale_DDBAR_train (0.36537/((N_DDBAR_train/(N_DDBAR_train + N_DDBAR_test))*1.0))
# define Scale_SSBAR_train (0.36537/((N_SSBAR_train/(N_SSBAR_train + N_SSBAR_test))*1.0))
# define Scale_CHARM_train (0.36537/((N_CHARM_train/(N_CHARM_train + N_CHARM_test))*1.0))
# define Scale_CHG_test (N_Bp / (2 * 2.8 * N_BpBp_1invab * (N_CHG_test / (N_CHG_train + N_CHG_test)) ))
# define Scale_MIX_test (N_B0 / (2 * 2.8 * N_B0B0_1invab * (N_MIX_test / (N_MIX_train + N_MIX_test)) ))
# define Scale_UUBAR_test (0.36537/((N_UUBAR_test/(N_UUBAR_train + N_UUBAR_test))*1.0))
# define Scale_DDBAR_test (0.36537/((N_DDBAR_test/(N_DDBAR_train + N_DDBAR_test))*1.0))
# define Scale_SSBAR_test (0.36537/((N_SSBAR_test/(N_SSBAR_train + N_SSBAR_test))*1.0))
# define Scale_CHARM_test (0.36537/((N_CHARM_test/(N_CHARM_train + N_CHARM_test))*1.0))
# define Scale_CHG_validation (N_Bp / (2 * 2.8 * N_BpBp_1invab))
# define Scale_MIX_validation (N_B0 / (2 * 2.8 * N_B0B0_1invab))
# define Scale_UUBAR_validation (0.36537)
# define Scale_DDBAR_validation (0.36537)
# define Scale_SSBAR_validation (0.36537)
# define Scale_CHARM_validation (0.36537)

// BKG MC sample number for MC15rd
# define N_CHG_validation_MC15rd 785108449.0 // 1458.959/fb
# define N_MIX_validation_MC15rd 741492304.0 // 1458.959/fb
# define N_UUBAR_validation_MC15rd 2306265848.0 // 1458.959/fb
# define N_DDBAR_validation_MC15rd 576209482.0 // 1458.959/fb
# define N_SSBAR_validation_MC15rd 526874294.0 // 1458.959/fb
# define N_CHARM_validation_MC15rd 1889822323.0 // 1458.959/fb

# define Scale_CHG_validation_MC15rd (0.36537/1.458959)
# define Scale_MIX_validation_MC15rd (0.36537/1.458959)
# define Scale_UUBAR_validation_MC15rd (0.36537/1.458959)
# define Scale_DDBAR_validation_MC15rd (0.36537/1.458959)
# define Scale_SSBAR_validation_MC15rd (0.36537/1.458959)
# define Scale_CHARM_validation_MC15rd (0.36537/1.458959)

// off-resonance sample for MC15rd
# define Scale_UUBAR_offres_MC15rd (0.04274/0.169328) // 169.328/fb
# define Scale_DDBAR_offres_MC15rd (0.04274/0.169328) // 169.328/fb
# define Scale_SSBAR_offres_MC15rd (0.04274/0.169328) // 169.328/fb
# define Scale_CHARM_offres_MC15rd (0.04274/0.169328) // 169.328/fb

// uncertainty of BR from parametric reason
# define Sigma_BR_Kplus_nunubar 0.000000402
# define Sigma_BR_K0star_nunubar 0.00000050
# define Sigma_BR_K0_nunubar (Sigma_BR_Kplus_nunubar*TB0/TBp)
# define Sigma_BR_Kplusstar_nunubar (Sigma_BR_K0star_nunubar*TBp/TB0)
# define Sigma_BR_Xs_nunubar 0.000003

// scale factor for systematic MC sample
# define N_K0star_nunubar_syst 10000000.0
# define N_Kplusstar_nunubar_syst 10000000.0
# define N_Xsu_nonresonant_nunubar_syst 10000000.0
# define N_Xsd_nonresonant_nunubar_syst 10000000.0

//# define KS0_rel_uncertainty 0.6 // %/cm
# define track_rel_uncertainty_MC15ri 0.24 // %
# define track_rel_uncertainty_MC15rd 0.27 // %
# define Kaon_PID_max_uncertainty 0.1 // not percentage. relative uncertainty

# define BR_KpKLKL_uncertainty (0.04/1.05) // not percentage. relative uncertainty
# define BR_KSKLKL_uncertainty (0.5/6.0) // not percentage. relative uncertainty

# define N_Knn_type 4 //  B2Knn B2Kstarnn B02K0nn B02K0starnn
const double B2Knn_up_uncer[N_Knn_type] = { // relative uncertainty
    0.32 / 2.66, 0.28 / 1.24, 0.5 / 5.9, 0.8 / 3.6
};
const double B2Knn_dn_uncer[N_Knn_type] = { // relative uncertainty
    0.32 / 2.66, 0.25 / 1.24, 0.5 / 5.9, 0.7 / 3.6
};

// scale factor for each systematic MC sample 
# define Scale_Xsu_nonresonant_syst (N_Xsu_nonresonant_nunubar_LS1/N_Xsu_nonresonant_nunubar_syst)
# define Scale_Xsd_nonresonant_syst (N_Xsd_nunubar_LS1/N_Xsd_nonresonant_nunubar_syst)

/* for the projection of analysis
* This factor is included at `corrector`, `NEWFEI_calculator`, and `Signal_yield_fit_BDT_Rarity_HistFactory`
* Because we are not interested in plots, we do not touch anything on `plot_integrated`.
* Therefore, `NEWFEIcal.txt` value does not change. The projection is applied on `NEWFEI_calculator`.
*/
# define projection_multiplication 1.0 // If it is 1.0, the result is based on MC before LS1. If it is 2.0, then we get the projection if we have 2 times of the MC before LS1.

#endif 