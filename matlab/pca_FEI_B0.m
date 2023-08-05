clear all
format long

N_total_bins = 11;
N_selected_vectors = 11;

% FEI
FEIfileID = fopen('FEI_cal_B0_cov_all.txt','r');
formatSpec = strjoin(repelem(["%f "], [N_total_bins]), '');
sizeFEI = [N_total_bins Inf];
CovFEI = fscanf(FEIfileID,formatSpec,sizeFEI);
CorrFEI = corrcov(CovFEI);
fclose(FEIfileID);

[FEIcoeff,FEIlatent] = pcacov(CovFEI);

% print selected KID info into file
FEIfileID_w = fopen('FEI_cal_B0_eigen.txt','w');
for j = 1:N_selected_vectors
    fprintf(FEIfileID_w,"%1.15f\n", sqrt(FEIlatent(j)));
    for i = 1:N_total_bins
        temp_eigen_vector = FEIcoeff(:,j);
        fprintf(FEIfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
end
fclose(FEIfileID_w);