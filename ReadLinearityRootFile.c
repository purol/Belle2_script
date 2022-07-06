#define LT_number 36

void ReadLinearityRootFile(){

    double Inputnsig[LT_number] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0 };
    double Inputnsigerror[LT_number] = { 0 };
    double outputnsig[LT_number] = { -0.657, -0.288, 1.81, 2.79, 1.78, 4.52, 5.8, 6.99, 6.26, 10.04, 7.41, 10.40, 10.48, 12.06, 13.6, 13.47, 18.0, 15.34, 17.74, 17.6, 18.48, 20.4, 22.14, 22.0, 23.43, 24.74, 27.4, 27.52, 27.6, 28.60, 31.1, 30.9, 30.3, 33.2, 35.1, 36.6 };
    double outputnsigerror[LT_number] = { 0.96, 0.97, 0.96, 0.99, 0.93, 0.92, 0.91, 0.95, 0.97, 0.94, 0.94, 0.96, 0.97, 0.97, 1.0, 0.97, 1.0, 0.97, 0.95, 1.0, 0.93, 1.0, 0.96, 1.0, 0.98, 0.99, 1.0, 0.99, 1.0, 1.00, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

    gStyle->SetOptFit(11); //gStyle->SetStatH(0.05);
    TCanvas* c = new TCanvas("Linearity test canvas", "", 800, 800);
    TGraphErrors* gr = new TGraphErrors(LT_number, Inputnsig, outputnsig, Inputnsigerror, outputnsigerror);
    gr->SetMarkerStyle(21); gr->SetTitle(";input n_{sig};output n_{sig}");
    gr->Fit("pol1"); gr->Draw("AP");
    c->SaveAs("Linearity_test.png");
}
