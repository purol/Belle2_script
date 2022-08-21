# define Nstep 20
# define start 0.9
# define end 1.0

void ReadFOMFiles1D(){

    const char* Kplus_large = "./FOM_plot_Kplus_large";
    const char* Kplusstar_large = "./FOM_plot_Kplusstar_large";
    const char* Xsu_nonres_large = "./FOM_plot_Xsu_large";
    const char* K0_large = "./FOM_plot_K0_large";
    const char* K0star_large = "./FOM_plot_K0star_large";
    const char* Xsd_nonres_large = "./FOM_plot_Xsd_large";
    const char* BKG_large = "./FOM_plot_BKG_large";

    const char* Kplus_small = "./FOM_plot_Kplus_small";
    const char* Kplusstar_small = "./FOM_plot_Kplusstar_small";
    const char* Xsu_nonres_small = "./FOM_plot_Xsu_small";
    const char* K0_small = "./FOM_plot_K0_small";
    const char* K0star_small = "./FOM_plot_K0star_small";
    const char* Xsd_nonres_small = "./FOM_plot_Xsd_small";
    const char* BKG_small = "./FOM_plot_BKG_small";

    double Matrix_SIGNAL_small[Nstep];
    for (int i = 0; i < Nstep; i++) Matrix_SIGNAL_small[i] = 0;

    double Matrix_SIGNAL_large[Nstep];
    for (int i = 0; i < Nstep; i++) Matrix_SIGNAL_large[i] = 0;

    double Matrix_BKG_small[Nstep];
    for (int i = 0; i < Nstep; i++) Matrix_BKG_small[i] = 0;

    double Matrix_BKG_large[Nstep];
    for (int i = 0; i < Nstep; i++) Matrix_BKG_large[i] = 0;

    double FOM_Matrix_small[Nstep];
    for (int i = 0; i < Nstep; i++) FOM_Matrix_small[i] = 0;

    double FOM_Matrix_large[Nstep];
    for (int i = 0; i < Nstep; i++) FOM_Matrix_large[i] = 0;

    double Matrix_SIGNAL[Nstep];
    for (int i = 0; i < Nstep; i++) Matrix_SIGNAL[i] = 0;

    double Matrix_BKG[Nstep];
    for (int i = 0; i < Nstep; i++) Matrix_BKG[i] = 0;

    double FOM_Matrix[Nstep];
    for (int i = 0; i < Nstep; i++) FOM_Matrix[i] = 0;
    
    /* ========================== read FOM files ========================== */
    FILE* pf;
    double temp_var = 0;
    pf = fopen(Kplus_small, "r");
    for (int i = 0; i < Nstep; i++) {
        fscanf(pf, "%lf ", &temp_var);
        Matrix_SIGNAL_small[i] = Matrix_SIGNAL_small[i] + temp_var;
        Matrix_SIGNAL[i] = Matrix_SIGNAL[i] + temp_var;
    }
    fclose(pf);

    pf = fopen(Kplusstar_small, "r");
    for (int i = 0; i < Nstep; i++) {
        fscanf(pf, "%lf ", &temp_var);
        Matrix_SIGNAL_small[i] = Matrix_SIGNAL_small[i] + temp_var;
        Matrix_SIGNAL[i] = Matrix_SIGNAL[i] + temp_var;
    }
    fclose(pf);

    pf = fopen(Xsu_nonres_small, "r");
    for (int i = 0; i < Nstep; i++) {
        fscanf(pf, "%lf ", &temp_var);
        Matrix_SIGNAL_small[i] = Matrix_SIGNAL_small[i] + temp_var;
        Matrix_SIGNAL[i] = Matrix_SIGNAL[i] + temp_var;
    }
    fclose(pf);

    pf = fopen(K0_small, "r");
    for (int i = 0; i < Nstep; i++) {
        fscanf(pf, "%lf ", &temp_var);
        Matrix_SIGNAL_small[i] = Matrix_SIGNAL_small[i] + temp_var;
        Matrix_SIGNAL[i] = Matrix_SIGNAL[i] + temp_var;
    }
    fclose(pf);

    pf = fopen(K0star_small, "r");
    for (int i = 0; i < Nstep; i++) {
        fscanf(pf, "%lf ", &temp_var);
        Matrix_SIGNAL_small[i] = Matrix_SIGNAL_small[i] + temp_var;
        Matrix_SIGNAL[i] = Matrix_SIGNAL[i] + temp_var;
    }
    fclose(pf);

    pf = fopen(Xsd_nonres_small, "r");
    for (int i = 0; i < Nstep; i++) {
        fscanf(pf, "%lf ", &temp_var);
        Matrix_SIGNAL_small[i] = Matrix_SIGNAL_small[i] + temp_var;
        Matrix_SIGNAL[i] = Matrix_SIGNAL[i] + temp_var;
    }
    fclose(pf);

    pf = fopen(BKG_small, "r");
    for (int i = 0; i < Nstep; i++) {
        fscanf(pf, "%lf ", &temp_var);
        Matrix_BKG_small[i] = Matrix_BKG_small[i] + temp_var;
        Matrix_BKG[i] = Matrix_BKG[i] + temp_var;
    }
    fclose(pf);


    pf = fopen(Kplus_large, "r");
    for (int i = 0; i < Nstep; i++) {
        fscanf(pf, "%lf ", &temp_var);
        Matrix_SIGNAL_large[i] = Matrix_SIGNAL_large[i] + temp_var;
        Matrix_SIGNAL[i] = Matrix_SIGNAL[i] + temp_var;
    }
    fclose(pf);

    pf = fopen(Kplusstar_large, "r");
    for (int i = 0; i < Nstep; i++) {
        fscanf(pf, "%lf ", &temp_var);
        Matrix_SIGNAL_large[i] = Matrix_SIGNAL_large[i] + temp_var;
        Matrix_SIGNAL[i] = Matrix_SIGNAL[i] + temp_var;
    }
    fclose(pf);

    pf = fopen(Xsu_nonres_large, "r");
    for (int i = 0; i < Nstep; i++) {
        fscanf(pf, "%lf ", &temp_var);
        Matrix_SIGNAL_large[i] = Matrix_SIGNAL_large[i] + temp_var;
        Matrix_SIGNAL[i] = Matrix_SIGNAL[i] + temp_var;
    }
    fclose(pf);

    pf = fopen(K0_large, "r");
    for (int i = 0; i < Nstep; i++) {
        fscanf(pf, "%lf ", &temp_var);
        Matrix_SIGNAL_large[i] = Matrix_SIGNAL_large[i] + temp_var;
        Matrix_SIGNAL[i] = Matrix_SIGNAL[i] + temp_var;
    }
    fclose(pf);

    pf = fopen(K0star_large, "r");
    for (int i = 0; i < Nstep; i++) {
        fscanf(pf, "%lf ", &temp_var);
        Matrix_SIGNAL_large[i] = Matrix_SIGNAL_large[i] + temp_var;
        Matrix_SIGNAL[i] = Matrix_SIGNAL[i] + temp_var;
    }
    fclose(pf);

    pf = fopen(Xsd_nonres_large, "r");
    for (int i = 0; i < Nstep; i++) {
        fscanf(pf, "%lf ", &temp_var);
        Matrix_SIGNAL_large[i] = Matrix_SIGNAL_large[i] + temp_var;
        Matrix_SIGNAL[i] = Matrix_SIGNAL[i] + temp_var;
    }
    fclose(pf);

    pf = fopen(BKG_large, "r");
    for (int i = 0; i < Nstep; i++) {
        fscanf(pf, "%lf ", &temp_var);
        Matrix_BKG_large[i] = Matrix_BKG_large[i] + temp_var;
        Matrix_BKG[i] = Matrix_BKG[i] + temp_var;
    }
    fclose(pf);
    /* ========================== read FOM files ========================== */
    
    for (int i = 0; i < Nstep; i++) FOM_Matrix_small[i] = Matrix_SIGNAL_small[i] / sqrt(Matrix_SIGNAL_small[i] + Matrix_BKG_small[i]);
    for (int i = 0; i < Nstep; i++) FOM_Matrix_large[i] = Matrix_SIGNAL_large[i] / sqrt(Matrix_SIGNAL_large[i] + Matrix_BKG_large[i]);
    for (int i = 0; i < Nstep; i++) FOM_Matrix[i] = Matrix_SIGNAL[i] / sqrt(Matrix_SIGNAL[i] + Matrix_BKG[i]);

    for (int i = 0; i < Nstep; i++) {
        printf("%f ", FOM_Matrix_small[i]);
        printf("\n");
    }

    printf("\n\n\n");

    for (int i = 0; i < Nstep; i++) {
        printf("%f ", FOM_Matrix_large[i]);
        printf("\n");
    }

    printf("\n\n\n");

    for (int i = 0; i < Nstep; i++) {
        printf("%f ", FOM_Matrix[i]);
        printf("\n");
    }

    double O_Continuum[Nstep];
    for (int i = 0; i < Nstep; i++) O_Continuum[i] = start + (end - start) * i / Nstep;

    TCanvas* c1 = new TCanvas("c1", "Graph Draw Options", 200, 10, 600, 600);

    TGraph* gr = new TGraph(Nstep, O_Continuum, FOM_Matrix_small);
    gr->SetTitle("M_{Xs} < 1.1 [GeV];O_{Continuum};#frac{S}{#sqrt{S+B}}");
    gr->SetMarkerStyle(0);
    gr->Draw("");
    c1->SaveAs("FOM_small.png");


    TCanvas* c2 = new TCanvas("c2", "Graph Draw Options", 200, 10, 600, 600);

    TGraph* gr2 = new TGraph(Nstep, O_Continuum, FOM_Matrix_large);
    gr2->SetTitle("M_{Xs} > 1.1 [GeV];O_{Continuum};#frac{S}{#sqrt{S+B}}");
    gr2->SetMarkerStyle(0);
    gr2->Draw("");
    c2->SaveAs("FOM_large.png");


    TCanvas* c3 = new TCanvas("c3", "Graph Draw Options", 200, 10, 600, 600);

    TGraph* gr3 = new TGraph(Nstep, O_Continuum, FOM_Matrix);
    gr3->SetTitle(";FBDT;#frac{S}{#sqrt{S+B}}");
    gr3->SetMarkerStyle(0);
    gr3->Draw("");
    c3->SaveAs("FOM_total.png");
}
