# define Nstep 25
# define start 0.25
# define end 0.45

void FillMatrix(const char* filename, double Matrix[Nstep][Nstep]) {
    FILE* pf;
    pf = fopen(filename, "r");
    for (int i = 0; i < Nstep; i++) {
        for (int j = 0; j < Nstep; j++) {
            double temp = 0;
            fscanf(pf, "%lf ", &temp);
            Matrix[i][j] = Matrix[i][j] + temp;
        }
        fscanf(pf, "\n");
    }
    fclose(pf);
}

void ReadFOMFiles(){

    const char* Kplus_small = "./FOM_plot_Kplus_small";
    const char* Kplusstar_small = "./FOM_plot_Kplusstar_small";
    const char* Xsu_nonres_small = "./FOM_plot_Xsu_small";
    const char* K0_small = "./FOM_plot_K0_small";
    const char* K0star_small = "./FOM_plot_K0star_small";
    const char* Xsd_nonres_small = "./FOM_plot_Xsd_small";
    const char* BKG_small = "./FOM_plot_BKG_small";

    const char* Kplus_large = "./FOM_plot_Kplus_large";
    const char* Kplusstar_large = "./FOM_plot_Kplusstar_large";
    const char* Xsu_nonres_large = "./FOM_plot_Xsu_large";
    const char* K0_large = "./FOM_plot_K0_large";
    const char* K0star_large = "./FOM_plot_K0star_large";
    const char* Xsd_nonres_large = "./FOM_plot_Xsd_large";
    const char* BKG_large = "./FOM_plot_BKG_large";

    double Matrix_SIGNAL[Nstep][Nstep];
    for (int i = 0; i < Nstep; i++) for (int j = 0; j < Nstep; j++) Matrix_SIGNAL[i][j] = 0;

    double Matrix_BKG[Nstep][Nstep];
    for (int i = 0; i < Nstep; i++) for (int j = 0; j < Nstep; j++) Matrix_BKG[i][j] = 0;

    double FOM_Matrix[Nstep][Nstep];
    for (int i = 0; i < Nstep; i++) for (int j = 0; j < Nstep; j++) FOM_Matrix[i][j] = 0;


    /* ========================== read FOM files ========================== */
    FillMatrix(Kplus_small, Matrix_SIGNAL);
    FillMatrix(Kplusstar_small, Matrix_SIGNAL);
    FillMatrix(Xsu_nonres_small, Matrix_SIGNAL);
    FillMatrix(K0_small, Matrix_SIGNAL);
    FillMatrix(K0star_small, Matrix_SIGNAL);
    FillMatrix(Xsd_nonres_small, Matrix_SIGNAL);
    FillMatrix(BKG_small, Matrix_BKG);

    FillMatrix(Kplus_large, Matrix_SIGNAL);
    FillMatrix(Kplusstar_large, Matrix_SIGNAL);
    FillMatrix(Xsu_nonres_large, Matrix_SIGNAL);
    FillMatrix(K0_large, Matrix_SIGNAL);
    FillMatrix(K0star_large, Matrix_SIGNAL);
    FillMatrix(Xsd_nonres_large, Matrix_SIGNAL);
    FillMatrix(BKG_large, Matrix_BKG);

    for (int i = 0; i < Nstep; i++) for (int j = 0; j < Nstep; j++) { 
        if(Matrix_SIGNAL[i][j] + Matrix_BKG[i][j] != 0 ) FOM_Matrix[i][j] = Matrix_SIGNAL[i][j] / sqrt(Matrix_SIGNAL[i][j] + Matrix_BKG[i][j]);
        else FOM_Matrix[i][j] = 0;
    }

    for (int i = 0; i < Nstep; i++) {
        for (int j = 0; j < Nstep; j++) {
            printf("%f ", Matrix_SIGNAL[i][j]);
        }
        printf("\n");
    }

    printf("\n\n\n");

    for (int i = 0; i < Nstep; i++) {
        for (int j = 0; j < Nstep; j++) {
            printf("%f ", Matrix_BKG[i][j]);
        }
        printf("\n");
    }

    printf("\n\n\n");

    for (int i = 0; i < Nstep; i++) {
        for (int j = 0; j < Nstep; j++) {
            printf("%f ", FOM_Matrix[i][j]);
        }
        printf("\n");
    }
}
