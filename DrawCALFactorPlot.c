int DrawCALFactorPlot(){

    double CAL[4] = { 0.661, 0.963355, 2.087818, 1.216868 };
    double CAL_error[4] = { 0.447531, 0.157872, 0.666075, 0.083132 };
    std::string sample_names[5] = {"proc12 chunk1 off-resonance", "proc12 chunk1 on-resonance sideband", "proc12 chunk2 off-resonance", "proc12 chunk2 on-resonance sideband", "combined" };

    double CAL_mean;
    double CAL_mean_error;

    for (int i = 0; i < 4; i++) {
        CAL_mean_error = CAL_mean_error + 1 / (CAL_error[i] * CAL_error[i]);
    }
    CAL_mean_error = 1 / std::sqrt(CAL_mean_error);

    for (int i = 0; i < 4; i++) {
        CAL_mean = CAL_mean + CAL[i] / (CAL_error[i] * CAL_error[i]);
    }
    CAL_mean = CAL_mean * CAL_mean_error * CAL_mean_error;

    TH1D* MC_dividedby_data = new TH1D("MC/data", "MC/data", 5, 0.0, 5);

    MC_dividedby_data->SetBinContent(1, CAL[0]);
    MC_dividedby_data->SetBinError(1, CAL_error[0]);

    MC_dividedby_data->SetBinContent(2, CAL[1]);
    MC_dividedby_data->SetBinError(2, CAL_error[1]);

    MC_dividedby_data->SetBinContent(3, CAL[2]);
    MC_dividedby_data->SetBinError(3, CAL_error[2]);

    MC_dividedby_data->SetBinContent(4, CAL[3]);
    MC_dividedby_data->SetBinError(4, CAL_error[3]);

    MC_dividedby_data->SetBinContent(5, CAL_mean);
    MC_dividedby_data->SetBinError(5, CAL_mean_error);

    for (int i = 0; i < 5; i++) {
        MC_dividedby_data->GetXaxis()->SetBinLabel(i + 1, sample_names[i].c_str());
    }

    TH1D* two_sigma = new TH1D("2sig hist", ";;", 5, 0.0, 5);
    for (int i = 0; i < 5; i++) {
        two_sigma->SetBinContent(i + 1, CAL_mean);
        two_sigma->SetBinError(i + 1, CAL_mean_error * 2);
    }
    two_sigma->SetFillColor(kYellow);
    two_sigma->SetFillStyle(1001);
    for (int i = 0; i < 5; i++) {
        two_sigma->GetXaxis()->SetBinLabel(i + 1, sample_names[i].c_str());
    }
    two_sigma->SetStats(false);
    two_sigma->GetYaxis()->SetTitleOffset(1.4);
    two_sigma->GetYaxis()->SetTitle("MC/data");
    //two_sigma->GetXaxis()->LabelsOption("v");



    TH1D* one_sigma = new TH1D("1sig hist", "1sig;;", 5, 0.0, 5);
    for (int i = 0; i < 5; i++) {
        one_sigma->SetBinContent(i + 1, CAL_mean);
        one_sigma->SetBinError(i + 1, CAL_mean_error);
    }
    one_sigma->SetFillColor(kGreen);
    one_sigma->SetFillStyle(1001);


    TCanvas* c = new TCanvas("c", "c", 600, 600); c->SetBottomMargin(0.3);

    MC_dividedby_data->SetLineWidth(2.0);
    MC_dividedby_data->SetMarkerColor(1);
    MC_dividedby_data->SetMarkerStyle(21);
    MC_dividedby_data->SetLineColor(1);

    two_sigma->GetYaxis()->SetRangeUser(0., 3.);
    two_sigma->Draw("E2");
    one_sigma->Draw("E2 same");
    MC_dividedby_data->Draw("e1 same");

    TLine* line = new TLine(0.0, CAL_mean, 5.0, CAL_mean);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2); line->SetLineWidth(1);
    line->Draw();

    TPaveText* pt = new TPaveText(0.129, 0.86, 0.494, 0.98, "NDC NB"); pt->SetFillStyle(0); pt->SetLineWidth(0); pt->AddText(("combined= " + std::to_string(CAL_mean) + " #pm " + std::to_string(CAL_mean_error)).c_str()); pt->Draw();

    return 0;
}
