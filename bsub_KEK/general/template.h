#ifndef TEMPLATE_H
#define TEMPLATE_H

#include "TH1.h"
#include <cmath>

// 0.0 < MXs < 0.6
# define RarityBins_MX1 6
# define BinMIN_MX1 0.96
# define BinMAX_MX1 1.0

// 0.6 < MXs < 1.0
# define RarityBins_MX2 6
# define BinMIN_MX2 0.96
# define BinMAX_MX2 1.0

// 1.0 < MXs < 2.0
# define RarityBins_MX3 4
# define BinMIN_MX3 0.96
# define BinMAX_MX3 1.0

const int RarityBins = (RarityBins_MX1 + RarityBins_MX2 + RarityBins_MX3);
const int BinMIN = 0;
const int BinMAX = RarityBins;

// TH1D* hist = new TH1D("histogram", ";;", RarityBins, BinMIN, BinMAX);

void FillTemplate(TH1D* hist, double FBDToutput, double total_weight, double MXs) {
	if ((MXs > 0.0) && (MXs < 0.6)) { // MXs1
		if ((FBDToutput >= BinMIN_MX1) && (FBDToutput <= BinMAX_MX1)) {
			double binwidth = (BinMAX_MX1 - BinMIN_MX1) / RarityBins_MX1;
			double binindex = std::ceil((FBDToutput - BinMIN_MX1) / binwidth); // start from 0. global bin index
			if (binindex >= RarityBins_MX1) {
				printf("Bin index is larger than expected. Just put maximum value.\n");
				binindex = RarityBins_MX1 - 1.0;
			}
			else if (binindex < 0) {
				printf("Bin index is smaller than expected. Just put minimum value.\n");
				binindex = 0.0;
			}
			hist->Fill(binindex + 0.5, total_weight); // add 0.5 to avoid unexpected rounding error
		}
		else {
			printf("FBDT output is out of range!\n");
			exit(1);
		}
	}
	else if ((MXs >= 0.6) && (MXs < 1.0)) { // MX2
		if ((FBDToutput >= BinMIN_MX2) && (FBDToutput <= BinMAX_MX2)) {
			double binwidth = (BinMAX_MX2 - BinMIN_MX2) / RarityBins_MX2;
			double binindex = RarityBins_MX1 + std::ceil((FBDToutput - BinMIN_MX2) / binwidth); // start from RarityBins_MX1. global bin index
			if (binindex >= (RarityBins_MX1 + RarityBins_MX2)) {
				printf("Bin index is larger than expected. Just put maximum value.\n");
				binindex = RarityBins_MX1 + RarityBins_MX2 - 1.0;
			}
			else if (binindex < 0) {
				printf("Bin index is smaller than expected. Just put minimum value.\n");
				binindex = RarityBins_MX1;
			}
			hist->Fill(binindex + 0.5, total_weight); // add 0.5 to avoid unexpected rounding error
		}
		else {
			printf("FBDT output is out of range!\n");
			exit(1);
		}
	}
	else if ((MXs >= 1.0) && (MXs < 2.0)) { // MX3
		if ((FBDToutput >= BinMIN_MX3) && (FBDToutput <= BinMAX_MX3)) {
			double binwidth = (BinMAX_MX3 - BinMIN_MX3) / RarityBins_MX3;
			double binindex = RarityBins_MX1 + RarityBins_MX2 + std::ceil((FBDToutput - BinMIN_MX3) / binwidth); // start from RarityBins_MX1 + RarityBins_MX2. global bin index
			if (binindex >= (RarityBins_MX1 + RarityBins_MX2 + RarityBins_MX3)) {
				printf("Bin index is larger than expected. Just put maximum value.\n");
				binindex = RarityBins_MX1 + RarityBins_MX2 + RarityBins_MX3 - 1.0;
			}
			else if (binindex < 0) {
				printf("Bin index is smaller than expected. Just put minimum value.\n");
				binindex = RarityBins_MX1 + RarityBins_MX2;
			}
			hist->Fill(binindex + 0.5, total_weight); // add 0.5 to avoid unexpected rounding error
		}
		else {
			printf("FBDT output is out of range!\n");
			exit(1);
		}
	}
	else {
		printf("MXs is out of range!\n");
		exit(1);
	}
}

#endif 