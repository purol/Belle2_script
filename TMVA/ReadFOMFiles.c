# define NormFactorForSignal (31481.0/9300000.0)

# define Nstep 20
# define start 0.8
# define end 1.0

void ReadFOMFiles(){

    const char* SIGNAL_small = "./FOM_plot_SIGNAL_smaller";
    const char* SIGNAL_large = "./FOM_plot_SIGNAL_larger";
    const char* BKG_small = "./FOM_plot_BKG_smaller";
    const char* BKG_large = "./FOM_plot_BKG_larger";

    double Matrix_SIGNAL_small[Nstep][Nstep];
    for (int i = 0; i < Nstep; i++) for (int j = 0; j < Nstep; j++) Matrix_SIGNAL_small[i][j] = 0;

    double Matrix_SIGNAL_large[Nstep][Nstep];
    for (int i = 0; i < Nstep; i++) for (int j = 0; j < Nstep; j++) Matrix_SIGNAL_large[i][j] = 0;

    double Matrix_BKG_small[Nstep][Nstep];
    for (int i = 0; i < Nstep; i++) for (int j = 0; j < Nstep; j++) Matrix_BKG_small[i][j] = 0;

    double Matrix_BKG_large[Nstep][Nstep];
    for (int i = 0; i < Nstep; i++) for (int j = 0; j < Nstep; j++) Matrix_BKG_large[i][j] = 0;

    double FOM_Matrix_small[Nstep][Nstep];
    for (int i = 0; i < Nstep; i++) for (int j = 0; j < Nstep; j++) FOM_Matrix_small[i][j] = 0;

    double FOM_Matrix_large[Nstep][Nstep];
    for (int i = 0; i < Nstep; i++) for (int j = 0; j < Nstep; j++) FOM_Matrix_large[i][j] = 0;

    FILE* pf;
    pf = fopen(SIGNAL_small, "r");
    for (int i = 0; i < Nstep; i++) {
        for (int j = 0; j < Nstep; j++) {
            fscanf(pf, "%lf ", &Matrix_SIGNAL_small[i][j]);
        }
        fscanf(pf, "\n");
    }
    fclose(pf);

    pf = fopen(SIGNAL_large, "r");
    for (int i = 0; i < Nstep; i++) {
        for (int j = 0; j < Nstep; j++) {
            fscanf(pf, "%lf ", &Matrix_SIGNAL_large[i][j]);
        }
        fscanf(pf, "\n");
    }
    fclose(pf);

    pf = fopen(BKG_small, "r");
    for (int i = 0; i < Nstep; i++) {
        for (int j = 0; j < Nstep; j++) {
            fscanf(pf, "%lf ", &Matrix_BKG_small[i][j]);
        }
        fscanf(pf, "\n");
    }
    fclose(pf);

    pf = fopen(BKG_large, "r");
    for (int i = 0; i < Nstep; i++) {
        for (int j = 0; j < Nstep; j++) {
            fscanf(pf, "%lf ", &Matrix_BKG_large[i][j]);
        }
        fscanf(pf, "\n");
    }
    fclose(pf);

    for (int i = 0; i < Nstep; i++) for (int j = 0; j < Nstep; j++) Matrix_SIGNAL_small[i][j] = Matrix_SIGNAL_small[i][j] * NormFactorForSignal;
    for (int i = 0; i < Nstep; i++) for (int j = 0; j < Nstep; j++) Matrix_SIGNAL_large[i][j] = Matrix_SIGNAL_large[i][j] * NormFactorForSignal;

    for (int i = 0; i < Nstep; i++) for (int j = 0; j < Nstep; j++) FOM_Matrix_small[i][j] = Matrix_SIGNAL_small[i][j] / sqrt(Matrix_SIGNAL_small[i][j] + Matrix_BKG_small[i][j]);
    for (int i = 0; i < Nstep; i++) for (int j = 0; j < Nstep; j++) FOM_Matrix_large[i][j] = Matrix_SIGNAL_large[i][j] / sqrt(Matrix_SIGNAL_large[i][j] + Matrix_BKG_large[i][j]);

    for (int i = 0; i < Nstep; i++) {
        for (int j = 0; j < Nstep; j++) {
            printf("%f ", FOM_Matrix_small[i][j]);
        }
        printf("\n");
    }

    printf("\n\n\n");

    for (int i = 0; i < Nstep; i++) {
        for (int j = 0; j < Nstep; j++) {
            printf("%f ", FOM_Matrix_large[i][j]);
        }
        printf("\n");
    }
}
