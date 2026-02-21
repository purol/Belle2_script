# define Nstep 100
# define start 0.5
# define end 1.0

void ReadFOMFiles1D(){

    const char* Kplus = "./FOM_plot_Kplus";
    const char* Kplusstar = "./FOM_plot_Kplusstar";
    const char* Xsu_nonres = "./FOM_plot_Xsu";
    const char* K0 = "./FOM_plot_K0";
    const char* K0star = "./FOM_plot_K0star";
    const char* Xsd_nonres = "./FOM_plot_Xsd";
    const char* CHG = "./FOM_plot_CHG";
    const char* MIX = "./FOM_plot_MIX";
    const char* UUBAR = "./FOM_plot_UUBAR";
    const char* DDBAR = "./FOM_plot_DDBAR";
    const char* SSBAR = "./FOM_plot_SSBAR";
    const char* CHARM = "./FOM_plot_CHARM";

    const char *dir_names[12] = { Kplus, Kplusstar, Xsu_nonres, K0, K0star, Xsd_nonres, CHG, MIX, UUBAR, DDBAR, SSBAR, CHARM };

    double Matrix_SIGNAL[Nstep];
    for (int i = 0; i < Nstep; i++) Matrix_SIGNAL[i] = 0;

    double Matrix_BKG[Nstep];
    for (int i = 0; i < Nstep; i++) Matrix_BKG[i] = 0;

    double FOM_Matrix[Nstep];
    for (int i = 0; i < Nstep; i++) FOM_Matrix[i] = 0;
    
    /* ========================== read FOM files ========================== */
    for (int i = 0; i < 6; i++) {
        FILE* pf;
        double temp_var = 0;
        pf = fopen(dir_names[i], "r");
        for (int i = 0; i < Nstep; i++) {
            fscanf(pf, "%lf ", &temp_var);
            Matrix_SIGNAL[i] = Matrix_SIGNAL[i] + temp_var;
        }
        fclose(pf);
    }

    for (int i = 0; i < 6; i++) {
        FILE* pf;
        double temp_var = 0;
        pf = fopen(dir_names[i + 6], "r");
        for (int i = 0; i < Nstep; i++) {
            fscanf(pf, "%lf ", &temp_var);
            Matrix_BKG[i] = Matrix_BKG[i] + temp_var;
        }
        fclose(pf);
    }

    /* ========================== read FOM files ========================== */
    
    for (int i = 0; i < Nstep; i++) FOM_Matrix[i] = Matrix_SIGNAL[i] / sqrt(Matrix_SIGNAL[i] + Matrix_BKG[i]);

    for (int i = 0; i < Nstep; i++) {
        printf("%f ", FOM_Matrix[i]);
        printf("\n");
    }

    double O_Continuum[Nstep];
    for (int i = 0; i < Nstep; i++) O_Continuum[i] = start + (end - start) * i / Nstep;

    TCanvas* c3 = new TCanvas("c3", "Graph Draw Options", 200, 10, 600, 600);

    TGraph* gr3 = new TGraph(Nstep, O_Continuum, FOM_Matrix);
    gr3->SetTitle(";FBDT;#frac{S}{#sqrt{S+B}}");
    gr3->SetMarkerStyle(0);
    gr3->SetMinimum(0.30);
    gr3->Draw("");

    TLine* line = new TLine(0.7, 0.30, 0.7, 0.60);
    line->SetLineColor(kRed);
    line->SetLineStyle(5);
    line->Draw();

    c3->SaveAs("FOM.png");
}
