clear all
format long

N_total_bins = 20 * 5;
N_selected_vectors = 3;

% FEI
FEIfileID = fopen('FEI_toys.txt','r');
formatSpec = strjoin(repelem(["%f "], [N_total_bins]), '');
sizeFEI = [N_total_bins Inf];
FEIMatrix = fscanf(FEIfileID,formatSpec,sizeFEI);
FEIMatrix = FEIMatrix';
CovFEI = cov(FEIMatrix);
CorrFEI = corrcov(CovFEI);
fclose(FEIfileID);

[FEIcoeff,FEIlatent] = pcacov(CovFEI);

% print covariance all
FEIfileID_w = fopen('FEI_cov_all.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(FEIfileID_w,"%1.15f ",CovFEI(i,j));
    end
    fprintf(FEIfileID_w,"\n");
end
fclose(FEIfileID_w);

% print truncated covariance matrix all
FEIfileID_w = fopen('FEI_cov_all_truncated.txt','w');
for i = 1:N_total_bins
    fprintf(FEIfileID_w,"%1.15f\n",CovFEI(i,i));
end
fclose(FEIfileID_w);

% print correlation all
FEIfileID_w = fopen('FEI_corr_all.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(FEIfileID_w,"%1.15f ",CorrFEI(i,j));
    end
    fprintf(FEIfileID_w,"\n");
end
fclose(FEIfileID_w);

% print selected KID info into file
FEIfileID_w = fopen('FEI_selected.txt','w');
for j = 1:N_selected_vectors
    fprintf(FEIfileID_w,"%1.15f\n", sqrt(FEIlatent(j)));
    for i = 1:N_total_bins
        temp_eigen_vector = FEIcoeff(:,j);
        fprintf(FEIfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
end
fclose(FEIfileID_w);

FEIpartialcov = zeros(N_total_bins, N_total_bins);
for j = 1:N_selected_vectors
    FEIpartialcov = FEIpartialcov + FEIcoeff(:,j) * FEIlatent(j) * FEIcoeff(:,j)';
end

FEIpartialcorr = corrcov(FEIpartialcov);

for j = 1:N_selected_vectors
    FEIcoeff(:,j) = 0;
    FEIlatent(j) = 0;
end

% print covariance partial
FEIfileID_w = fopen('FEI_cov_partial.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(FEIfileID_w,"%1.15f ",FEIpartialcov(i,j));
    end
    fprintf(FEIfileID_w,"\n");
end
fclose(FEIfileID_w);

% print correlation partial
FEIfileID_w = fopen('FEI_corr_partial.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(FEIfileID_w,"%1.15f ",FEIpartialcorr(i,j));
    end
    fprintf(FEIfileID_w,"\n");
end
fclose(FEIfileID_w);

% get covariance for remaining
FEIremaincov = FEIcoeff*diag(FEIlatent)*FEIcoeff';

% print truncated remaining covariance matrix
FEIfileID_w = fopen('FEI_cov_remain_truncated.txt','w');
for i = 1:N_total_bins
    fprintf(FEIfileID_w,"%1.15f\n",FEIremaincov(i,i));
end
fclose(FEIfileID_w);

% latent: eigen value
% coeff: eigen vector