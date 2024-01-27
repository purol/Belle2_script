// last update: 2021-10-13
// for Belle2 data

void DrawEfficiency(){
    
    // draw my MC data
    double x_Eff[9] = { 0.0 };
    double y_Eff[9] = { 0.00225531482, 0.00074224593, 0.00056996206, 0.00038427068, 0.00023724124, 0.00011063264, 0.00004409573, 0.00002578873, 0.00002414053 };

    // create the statistical error for my MC data
    double ex_Eff[9] = { 0.14166666666, 0.14166666666, 0.14166666666, 0.14166666666, 0.14166666666, 0.14166666666, 0.14166666666, 0.14166666666, 0.14166666666 };
    double ey_Eff[9] = { 0.00693455447, 0.01259696666, 0.01650860019, 0.0180980196, 0.02547758934, 0.0458877967, 0.0982922851, 0.18784792295, 0.34000834833 }; // relative uncertainty

    for (int i = 0; i < 9; i++) {
        x_Eff[i] = 0.45 + 0.14166666666 + i * (3.0 - 0.45) / 9;
        ey_Eff[i] = ey_Eff[i] * y_Eff[i];
    }

    TGraphErrors* gr_Eff = new TGraphErrors(9, x_Eff, y_Eff, ex_Eff, ey_Eff);
    gr_Eff->SetMarkerSize(0.0);
    gr_Eff->SetLineWidth(2);
    gr_Eff->SetMarkerColor(4);
    gr_Eff->SetMarkerStyle(21);

    TCanvas* c_temp = new TCanvas("c", "", 1200, 1000); c_temp->cd();
    gr_Eff->SetTitle(";M_{Xs};Efficiency");
    gr_Eff->Draw("AP");
    c_temp->SaveAs("Efficiency.png");
    
    // corrector.Drawq2Plot(false, "Bplus");
    // corrector.Drawq2Plot(true, "Bplus");
}
