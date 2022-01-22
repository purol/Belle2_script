# define NormFactorForSignal (31481.0/9300000.0)

# define Nstep 20
# define start 0.6
# define end 1.0

void ReadFOMFiles1D(){

    const char* SIGNAL_small = "./FOM_SIGNAL_small";
    const char* SIGNAL_large = "./FOM_SIGNAL_large";
    const char* BKG_small = "./FOM_BKG_small";
    const char* BKG_large = "./FOM_BKG_large";

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
    
    FILE* pf;
    pf = fopen(SIGNAL_small, "r");
    for (int i = 0; i < Nstep; i++) {
        fscanf(pf, "%lf ", &Matrix_SIGNAL_small[i]);
    }
    fclose(pf);
    
    pf = fopen(SIGNAL_large, "r");
    for (int i = 0; i < Nstep; i++) {
        fscanf(pf, "%lf ", &Matrix_SIGNAL_large[i]);
    }
    fclose(pf);

    pf = fopen(BKG_small, "r");
    for (int i = 0; i < Nstep; i++) {
        fscanf(pf, "%lf ", &Matrix_BKG_small[i]);
    }
    fclose(pf);

    pf = fopen(BKG_large, "r");
    for (int i = 0; i < Nstep; i++) {
        fscanf(pf, "%lf ", &Matrix_BKG_large[i]);
    }
    fclose(pf);
    
    for (int i = 0; i < Nstep; i++) Matrix_SIGNAL_small[i] = Matrix_SIGNAL_small[i] * NormFactorForSignal;
    for (int i = 0; i < Nstep; i++) Matrix_SIGNAL_large[i] = Matrix_SIGNAL_large[i] * NormFactorForSignal;

    for (int i = 0; i < Nstep; i++) FOM_Matrix_small[i] = Matrix_SIGNAL_small[i] / sqrt(Matrix_SIGNAL_small[i] + Matrix_BKG_small[i]);
    for (int i = 0; i < Nstep; i++) FOM_Matrix_large[i] = Matrix_SIGNAL_large[i] / sqrt(Matrix_SIGNAL_large[i] + Matrix_BKG_large[i]);

    for (int i = 0; i < Nstep; i++) {
        printf("%f ", FOM_Matrix_small[i]);
        printf("\n");
    }

    printf("\n\n\n");

    for (int i = 0; i < Nstep; i++) {
        printf("%f ", FOM_Matrix_large[i]);
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
    c1->SaveAs("FOM_large.png");
}
