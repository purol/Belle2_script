#ifndef TEMPLATE_H
#define TEMPLATE_H

#include "TH1.h"
#include <cmath>

// 0.0 < MXs < 0.6
# define RarityBins_MX1 5
# define BinMIN_MX1 0.86
# define BinMAX_MX1 1.0
const double boundary_MXs1[RarityBins_MX1 + 1] = { BinMIN_MX1, 0.9195, 0.979, 0.986, 0.993, BinMAX_MX1 };

// 0.6 < MXs < 1.0
# define RarityBins_MX2 5
# define BinMIN_MX2 0.86
# define BinMAX_MX2 1.0
const double boundary_MXs2[RarityBins_MX2 + 1] = { BinMIN_MX2, 0.8925, 0.925, 0.95, 0.975, BinMAX_MX2 };

// 1.0 < MXs < 2.0
# define RarityBins_MX3 5
# define BinMIN_MX3 0.86
# define BinMAX_MX3 1.0
const double boundary_MXs3[RarityBins_MX3 + 1] = { BinMIN_MX3, 0.9, 0.94, 0.96, 0.98, BinMAX_MX3 };

const int RarityBins = (RarityBins_MX1 + RarityBins_MX2 + RarityBins_MX3);
const int BinMIN = 0;
const int BinMAX = RarityBins;

// TH1D* hist = new TH1D("histogram", ";;", RarityBins, BinMIN, BinMAX);

double ReturnBinIndex(double FBDToutput, double MXs) {
	if ((MXs > 0.0) && (MXs < 0.6)) { // MXs1
		if ((FBDToutput > BinMIN_MX1) && (FBDToutput <= BinMAX_MX1)) {
			for (int i = 0; i < RarityBins_MX1; i++) {
				if ((FBDToutput > boundary_MXs1[i]) && (FBDToutput <= boundary_MXs1[i + 1])) return ((double)(i) + 0.5);
			}
		}
		else {
			printf("FBDT output is out of range!\n");
			exit(1);
		}
	}
	else if ((MXs >= 0.6) && (MXs < 1.0)) { // MX2
		if ((FBDToutput >= BinMIN_MX2) && (FBDToutput <= BinMAX_MX2)) {
			for (int i = 0; i < RarityBins_MX2; i++) {
				if ((FBDToutput > boundary_MXs2[i]) && (FBDToutput <= boundary_MXs2[i + 1])) return ((double)(i + RarityBins_MX1) + 0.5);
			}
		}
		else {
			printf("FBDT output is out of range!\n");
			exit(1);
		}
	}
	else if ((MXs >= 1.0) && (MXs < 2.0)) { // MX3
		if ((FBDToutput >= BinMIN_MX3) && (FBDToutput <= BinMAX_MX3)) {
			for (int i = 0; i < RarityBins_MX3; i++) {
				if ((FBDToutput > boundary_MXs3[i]) && (FBDToutput <= boundary_MXs3[i + 1])) return ((double)(i + RarityBins_MX1 + RarityBins_MX2) + 0.5);
			}
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
	exit(1);
	return -1;
}

double FillTemplate(TH1D* hist, double FBDToutput, double total_weight, double MXs) {
	double binindex = ReturnBinIndex(FBDToutput, MXs);
	hist->Fill(binindex, total_weight); // add 0.5 to avoid unexpected rounding error
	return total_weight;
}

double DoesNotFillTemplate(double FBDToutput, double total_weight, double MXs) {
	double binindex = ReturnBinIndex(FBDToutput, MXs);
	return total_weight;
}

#endif 