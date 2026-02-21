// last update: 2021-10-13
// for Belle2 data

class Corrector {
private:
    const int STEP;

    const double mKp;
    const double mK0;

    const double mBp;
    const double mB0;

    const double alpha0_old;
    const double alpha1_old;
    const double alpha2_old;
    const double mp_Bp_old;
    const double mp_B0_old;

    const double alpha0_new; // 0.2545 +- 0.0090
    const double alpha1_new; // -0.71 +- 0.14
    const double alpha2_new; // 0.32 +- 0.59
    const double mp_B_new; //  5.4158 +- 0.0015

    double Total_Bp_old;
    double Total_B0_old;
    double Total_Bp_new;
    double Total_B0_new;
public:
    Corrector();
    double GetCorrectionFactor(double q2, std::string type);
    void Drawq2Plot(bool IsItNew, std::string type);
};

Corrector corrector;

Corrector::Corrector() :
    STEP(100),
    mKp(0.493677),
    mK0(0.497611),
    mBp(5.27934),
    mB0(5.27965),
    alpha0_old(0.432),
    alpha1_old(-0.664),
    alpha2_old(-1.2),
    mp_Bp_old(5.27934 + 0.046),
    mp_B0_old(5.27965 + 0.046),
    alpha0_new(0.2545),
    alpha1_new(-0.71),
    alpha2_new(0.32),
    mp_B_new(5.4158)
{
    Total_Bp_old = 0;
    Total_B0_old = 0;
    Total_Bp_new = 0;
    Total_B0_new = 0;

    // Bp old
    for (int i = 0; i < STEP; i++) {
        double q2 = ((double)i) * (mBp - mKp) * (mBp - mKp) / ((double)STEP);
        double delq2 = (mBp - mKp) * (mBp - mKp) / ((double)STEP);

        double tp = (mBp + mKp) * (mBp + mKp);
        double tm = (mBp - mKp) * (mBp - mKp);
        double t0 = tp * (1 - std::sqrt(1 - tm / tp));
        double z = (std::sqrt(tp - q2) - std::sqrt(tp - t0)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0));
        double fp = (1 / (1 - q2 / (mp_Bp_old * mp_Bp_old))) * (alpha0_old + alpha1_old * z + alpha2_old * z * z + (-alpha1_old + 2 * alpha2_old) * z * z * z / 3);
        double lambda = std::pow(mBp, 4) + std::pow(mKp, 4) + std::pow(q2, 2) - 2 * (std::pow(mBp, 2) * std::pow(mKp, 2) + std::pow(mKp, 2) * q2 + std::pow(mBp, 2) * q2);

        Total_Bp_old = Total_Bp_old + std::pow(lambda, 1.5) * std::pow(fp, 2) * delq2;
    }

    // B0 old
    for (int i = 0; i < STEP; i++) {
        double q2 = ((double)i) * (mB0 - mK0) * (mB0 - mK0) / ((double)STEP);
        double delq2 = (mB0 - mK0) * (mB0 - mK0) / ((double)STEP);

        double tp = (mB0 + mK0) * (mB0 + mK0);
        double tm = (mB0 - mK0) * (mB0 - mK0);
        double t0 = tp * (1 - std::sqrt(1 - tm / tp));
        double z = (std::sqrt(tp - q2) - std::sqrt(tp - t0)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0));
        double fp = (1 / (1 - q2 / (mp_B0_old * mp_B0_old))) * (alpha0_old + alpha1_old * z + alpha2_old * z * z + (-alpha1_old + 2 * alpha2_old) * z * z * z / 3);
        double lambda = std::pow(mB0, 4) + std::pow(mK0, 4) + std::pow(q2, 2) - 2 * (std::pow(mB0, 2) * std::pow(mK0, 2) + std::pow(mK0, 2) * q2 + std::pow(mB0, 2) * q2);

        Total_B0_old = Total_B0_old + std::pow(lambda, 1.5) * std::pow(fp, 2) * delq2;
    }

    // Bp new
    for (int i = 0; i < STEP; i++) {
        double q2 = ((double)i) * (mBp - mKp) * (mBp - mKp) / ((double)STEP);
        double delq2 = (mBp - mKp) * (mBp - mKp) / ((double)STEP);

        double tp = (mBp + mKp) * (mBp + mKp);
        double tm = (mBp - mKp) * (mBp - mKp);
        double t0 = 0.0;
        double z = (std::sqrt(tp - q2) - std::sqrt(tp - t0)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0));
        double fp = (1 / (1 - q2 / (mp_B_new * mp_B_new))) * (alpha0_new + alpha1_new * z + alpha2_new * z * z + (-alpha1_new + 2 * alpha2_new) * z * z * z / 3);
        double lambda = std::pow(mBp, 4) + std::pow(mKp, 4) + std::pow(q2, 2) - 2 * (std::pow(mBp, 2) * std::pow(mKp, 2) + std::pow(mKp, 2) * q2 + std::pow(mBp, 2) * q2);

        Total_Bp_new = Total_Bp_new + std::pow(lambda, 1.5) * std::pow(fp, 2) * delq2;
    }

    // B0 new
    for (int i = 0; i < STEP; i++) {
        double q2 = ((double)i) * (mB0 - mK0) * (mB0 - mK0) / ((double)STEP);
        double delq2 = (mB0 - mK0) * (mB0 - mK0) / ((double)STEP);

        double tp = (mB0 + mK0) * (mB0 + mK0);
        double tm = (mB0 - mK0) * (mB0 - mK0);
        double t0 = 0.0;
        double z = (std::sqrt(tp - q2) - std::sqrt(tp - t0)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0));
        double fp = (1 / (1 - q2 / (mp_B_new * mp_B_new))) * (alpha0_new + alpha1_new * z + alpha2_new * z * z + (-alpha1_new + 2 * alpha2_new) * z * z * z / 3);
        double lambda = std::pow(mB0, 4) + std::pow(mK0, 4) + std::pow(q2, 2) - 2 * (std::pow(mB0, 2) * std::pow(mK0, 2) + std::pow(mK0, 2) * q2 + std::pow(mB0, 2) * q2);

        Total_B0_new = Total_B0_new + std::pow(lambda, 1.5) * std::pow(fp, 2) * delq2;
    }
}

double Corrector::GetCorrectionFactor(double q2, std::string type) {
    if (type == "Bplus") {
        double tp = (mBp + mKp) * (mBp + mKp);
        double tm = (mBp - mKp) * (mBp - mKp);
        double t0_old = tp * (1 - std::sqrt(1 - tm / tp));
        double t0_new = 0.0;
        double z_old = (std::sqrt(tp - q2) - std::sqrt(tp - t0_old)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0_old));
        double z_new = (std::sqrt(tp - q2) - std::sqrt(tp - t0_new)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0_new));

        double fp_old = (1 / (1 - q2 / (mp_Bp_old * mp_Bp_old))) * (alpha0_old + alpha1_old * z_old + alpha2_old * z_old * z_old + (-alpha1_old + 2 * alpha2_old) * z_old * z_old * z_old / 3);
        double fp_new = (1 / (1 - q2 / (mp_B_new * mp_B_new))) * (alpha0_new + alpha1_new * z_new + alpha2_new * z_new * z_new + (-alpha1_new + 2 * alpha2_new) * z_new * z_new * z_new / 3);
        double lambda = std::pow(mBp, 4) + std::pow(mKp, 4) + std::pow(q2, 2) - 2 * (std::pow(mBp, 2) * std::pow(mKp, 2) + std::pow(mKp, 2) * q2 + std::pow(mBp, 2) * q2);

        if ((std::pow(lambda, 1.5) * std::pow(fp_old, 2) / Total_Bp_old) == 0) {
            printf("[ERROR] try to divide by 0!\n");
            exit(1);
        }
        return (std::pow(lambda, 1.5) * std::pow(fp_new, 2) / Total_Bp_new) / (std::pow(lambda, 1.5) * std::pow(fp_old, 2) / Total_Bp_old);
    }
    else if (type == "Bzero") {
        double tp = (mB0 + mK0) * (mB0 + mK0);
        double tm = (mB0 - mK0) * (mB0 - mK0);
        double t0_old = tp * (1 - std::sqrt(1 - tm / tp));
        double t0_new = 0.0;
        double z_old = (std::sqrt(tp - q2) - std::sqrt(tp - t0_old)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0_old));
        double z_new = (std::sqrt(tp - q2) - std::sqrt(tp - t0_new)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0_new));

        double fp_old = (1 / (1 - q2 / (mp_B0_old * mp_B0_old))) * (alpha0_old + alpha1_old * z_old + alpha2_old * z_old * z_old + (-alpha1_old + 2 * alpha2_old) * z_old * z_old * z_old / 3);
        double fp_new = (1 / (1 - q2 / (mp_B_new * mp_B_new))) * (alpha0_new + alpha1_new * z_new + alpha2_new * z_new * z_new + (-alpha1_new + 2 * alpha2_new) * z_new * z_new * z_new / 3);
        double lambda = std::pow(mB0, 4) + std::pow(mK0, 4) + std::pow(q2, 2) - 2 * (std::pow(mB0, 2) * std::pow(mK0, 2) + std::pow(mK0, 2) * q2 + std::pow(mB0, 2) * q2);

        if ((std::pow(lambda, 1.5) * std::pow(fp_old, 2) / Total_B0_old) == 0) {
            printf("[ERROR] try to divide by 0!\n");
            exit(1);
        }
        return (std::pow(lambda, 1.5) * std::pow(fp_new, 2) / Total_B0_new) / (std::pow(lambda, 1.5) * std::pow(fp_old, 2) / Total_B0_old);
    }

    printf("[ERROR] undefined type name");
    exit(1);

    return 0;
}

void Corrector::Drawq2Plot(bool IsItNew, std::string type) {

    Float_t Lower[8];

    Lower[0] = 0;
    Lower[1] = 4;
    Lower[2] = 8;
    Lower[3] = 12;
    Lower[4] = 16;
    Lower[5] = 20;
    Lower[6] = 22.9;
    Lower[7] = 25.1;
    TH1D* Plot = new TH1D("q2", ";q^{2} [GeV^{2}]; evt", 7, Lower);

    if (IsItNew == false && type == "Bplus") { // Bp old
        for (int i = 0; i < STEP; i++) {
            double q2 = ((double)i) * (mBp - mKp) * (mBp - mKp) / ((double)STEP);
            double delq2 = (mBp - mKp) * (mBp - mKp) / ((double)STEP);

            double tp = (mBp + mKp) * (mBp + mKp);
            double tm = (mBp - mKp) * (mBp - mKp);
            double t0 = tp * (1 - std::sqrt(1 - tm / tp));
            double z = (std::sqrt(tp - q2) - std::sqrt(tp - t0)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0));
            double fp = (1 / (1 - q2 / (mp_Bp_old * mp_Bp_old))) * (alpha0_old + alpha1_old * z + alpha2_old * z * z + (-alpha1_old + 2 * alpha2_old) * z * z * z / 3);
            double lambda = std::pow(mBp, 4) + std::pow(mKp, 4) + std::pow(q2, 2) - 2 * (std::pow(mBp, 2) * std::pow(mKp, 2) + std::pow(mKp, 2) * q2 + std::pow(mBp, 2) * q2);

            Plot->Fill(q2, std::pow(lambda, 1.5) * std::pow(fp, 2) * delq2);
        }
    }
    else if (IsItNew == false && type == "Bzero") { // B0 old
        for (int i = 0; i < STEP; i++) {
            double q2 = ((double)i) * (mB0 - mK0) * (mB0 - mK0) / ((double)STEP);
            double delq2 = (mB0 - mK0) * (mB0 - mK0) / ((double)STEP);

            double tp = (mB0 + mK0) * (mB0 + mK0);
            double tm = (mB0 - mK0) * (mB0 - mK0);
            double t0 = tp * (1 - std::sqrt(1 - tm / tp));
            double z = (std::sqrt(tp - q2) - std::sqrt(tp - t0)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0));
            double fp = (1 / (1 - q2 / (mp_B0_old * mp_B0_old))) * (alpha0_old + alpha1_old * z + alpha2_old * z * z + (-alpha1_old + 2 * alpha2_old) * z * z * z / 3);
            double lambda = std::pow(mB0, 4) + std::pow(mK0, 4) + std::pow(q2, 2) - 2 * (std::pow(mB0, 2) * std::pow(mK0, 2) + std::pow(mK0, 2) * q2 + std::pow(mB0, 2) * q2);

            Plot->Fill(q2, std::pow(lambda, 1.5) * std::pow(fp, 2) * delq2);
        }
    }
    else if (IsItNew == true && type == "Bplus") { // Bp new
        for (int i = 0; i < STEP; i++) {
            double q2 = ((double)i) * (mBp - mKp) * (mBp - mKp) / ((double)STEP);
            double delq2 = (mBp - mKp) * (mBp - mKp) / ((double)STEP);

            double tp = (mBp + mKp) * (mBp + mKp);
            double tm = (mBp - mKp) * (mBp - mKp);
            double t0 = 0.0;
            double z = (std::sqrt(tp - q2) - std::sqrt(tp - t0)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0));
            double fp = (1 / (1 - q2 / (mp_B_new * mp_B_new))) * (alpha0_new + alpha1_new * z + alpha2_new * z * z + (-alpha1_new + 2 * alpha2_new) * z * z * z / 3);
            double lambda = std::pow(mBp, 4) + std::pow(mKp, 4) + std::pow(q2, 2) - 2 * (std::pow(mBp, 2) * std::pow(mKp, 2) + std::pow(mKp, 2) * q2 + std::pow(mBp, 2) * q2);

            Plot->Fill(q2, std::pow(lambda, 1.5) * std::pow(fp, 2) * delq2);
        }
    }
    else if (IsItNew == true && type == "Bzero") { // B0 new
        for (int i = 0; i < STEP; i++) {
            double q2 = ((double)i) * (mB0 - mK0) * (mB0 - mK0) / ((double)STEP);
            double delq2 = (mB0 - mK0) * (mB0 - mK0) / ((double)STEP);

            double tp = (mB0 + mK0) * (mB0 + mK0);
            double tm = (mB0 - mK0) * (mB0 - mK0);
            double t0 = 0.0;
            double z = (std::sqrt(tp - q2) - std::sqrt(tp - t0)) / (std::sqrt(tp - q2) + std::sqrt(tp - t0));
            double fp = (1 / (1 - q2 / (mp_B_new * mp_B_new))) * (alpha0_new + alpha1_new * z + alpha2_new * z * z + (-alpha1_new + 2 * alpha2_new) * z * z * z / 3);
            double lambda = std::pow(mB0, 4) + std::pow(mK0, 4) + std::pow(q2, 2) - 2 * (std::pow(mB0, 2) * std::pow(mK0, 2) + std::pow(mK0, 2) * q2 + std::pow(mB0, 2) * q2);

            Plot->Fill(q2, std::pow(lambda, 1.5) * std::pow(fp, 2) * delq2);
        }
    }

    TCanvas* c_temp = new TCanvas("c", "", 1200, 1200); c_temp->cd();
    Plot->Draw();
    if(IsItNew == true) c_temp->SaveAs( ("q2distribution_new_" + type + ".png").c_str() );
    else c_temp->SaveAs(("q2distribution_old_" + type + ".png").c_str());

    delete Plot;
    delete c_temp;
}

void load_files(const char *dirname, std::vector<string>* names){
   TSystemDirectory dir(dirname, dirname);
   TList *files = dir.GetListOfFiles();
   if (files) {
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next())) {
         fname = file->GetName();
         if (!file->IsDirectory() && fname.EndsWith(".root")) {
            names->push_back(fname.Data());
         }
      }
   }
}

void ReadInfoFiles_r_sp(){
    
    std::vector<string> names;
    const char* dirname = "./";

    load_files(dirname, &names);

    Float_t Lower[8];

    Lower[0] = 0;
    Lower[1] = 4;
    Lower[2] = 8;
    Lower[3] = 12;
    Lower[4] = 16;
    Lower[5] = 20;
    Lower[6] = 22.9;
    Lower[7] = 25.1;
    TH1D* Plot = new TH1D("q2", ";q^{2} [GeV^{2}]; evt", 7, Lower);

    for(unsigned int i = 0; i<names.size(); i++){

        TFile *input_file = new TFile( (dirname+std::string("/")+names.at(i)).c_str(),"read");
        printf("%s (%d/%zu)\n",("Read "+names.at(i) + "... ").c_str(), i, names.size());

        // Get Data
        TTree* tree_info = (TTree*)input_file->Get("info");

        int __experiment__;
        int __run__;
        int __event__;
        int __candidate__;
        int __ncandidates__;
        double s;

        // get event_info
        tree_info->SetBranchAddress("__experiment__", &__experiment__);
        tree_info->SetBranchAddress("__run__", &__run__);
        tree_info->SetBranchAddress("__event__", &__event__);
        tree_info->SetBranchAddress("__candidate__", &__candidate__);
        tree_info->SetBranchAddress("__ncandidates__", &__ncandidates__);
        tree_info->SetBranchAddress("invMassInLists__bonu_e__clPrimaryMC_signal__bc", &s);

        printf("%lld entries...\n", tree_info->GetEntries());
        for (unsigned int j = 0; j < tree_info->GetEntries(); j++) { // Fill
            tree_info->GetEntry(j);
            double correction_weight = corrector.GetCorrectionFactor(s * s, "Bplus");
            Plot->Fill(s*s, correction_weight);
        }
        input_file->Close();
    }

    Double_t factor = 5.044;
    Float_t ratio = factor / Plot->GetEntries();
    
    // draw my MC data
    double x_MC[7] = { 2, 6, 10, 14, 18, 21.45, 24 };
    double y_MC[7] = { Plot->GetBinContent(1) * ratio, Plot->GetBinContent(2) * ratio, Plot->GetBinContent(3) * ratio, Plot->GetBinContent(4) * ratio, Plot->GetBinContent(5) * ratio, Plot->GetBinContent(6) * ratio, Plot->GetBinContent(7) * ratio };

    // create the statistical error for my MC data
    double ex_MC[7] = { 2, 2, 2, 2, 2, 1.45, 1.1 };
    double ey_MC[7] = { std::sqrt(Plot->GetBinContent(1)) * ratio, std::sqrt(Plot->GetBinContent(2)) * ratio, std::sqrt(Plot->GetBinContent(3)) * ratio, std::sqrt(Plot->GetBinContent(4)) * ratio, std::sqrt(Plot->GetBinContent(5)) * ratio, std::sqrt(Plot->GetBinContent(6)) * ratio, std::sqrt(Plot->GetBinContent(7)) * ratio };

    TGraphErrors* gr_MC = new TGraphErrors(7, x_MC, y_MC, ex_MC, ey_MC);
    gr_MC->SetMarkerColor(4);
    gr_MC->SetMarkerStyle(21);

    //Plot->Scale(factor / Plot->GetEntries());
    //Plot->SetMaximum(1.2);
    //Plot->SetMinimum(0.0);
    Plot->SetMarkerColor(4);
    Plot->SetMarkerStyle(21);
    //Plot->SetMarkerSize(2);

    // draw theoratical data
    Float_t x[6] = { 2, 6, 10, 14, 18, 21.45 };
    Float_t y[6] = { 1.189, 1.155, 1.071, 0.905, 0.597, 0.127 };

    // create the error arrays
    Float_t ex[6] = { 2, 2, 2, 2, 2, 1.45 };
    Float_t ey[6] = { 0.097, 0.090, 0.084, 0.072, 0.048, 0.011 };
    // uncertainty by form factor

    // create the TGraphErrors and draw it
    TGraphErrors* gr = new TGraphErrors(6, x, y, ex, ey);
    //gr->SetMarkerColor(4);
    //gr->SetMarkerStyle(21);
    gr->SetFillColor(2);
    gr->SetFillStyle(3001);
    gr->SetMinimum(0.0);
    gr->SetTitle(";q^{2} [GeV^{2}];arbitrary unit");

    TCanvas* c_temp = new TCanvas("c", "", 1200, 1200); c_temp->cd();
    TLegend *legend = new TLegend(0.35, 0.35, 0.1, 0.1);
    legend->AddEntry(gr_MC,"MC","lpfe");
    legend->AddEntry(gr,"[PHYS. REV. D 107, 014511]","f");
    gr->Draw("a2");
    //Plot->Draw("E1 SAME");
    gr_MC->Draw("SAMEP");
    legend->Draw("SAME");
    c_temp->SaveAs("q2distribution.png");
    
    // corrector.Drawq2Plot(false, "Bplus");
    // corrector.Drawq2Plot(true, "Bplus");
}
