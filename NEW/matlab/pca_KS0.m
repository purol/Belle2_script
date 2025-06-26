clear all
format long

N_total_bins = 15 * 9;
N_selected_vectors = 3;

% KS0 correction
KS0fileID = fopen('KS0_toys.txt','r');
formatSpec = strjoin(repelem(["%f "], [N_total_bins]), '');
sizeKS0 = [N_total_bins Inf];
KS0Matrix = fscanf(KS0fileID,formatSpec,sizeKS0);
KS0Matrix = KS0Matrix';
CovKS0 = cov(KS0Matrix);
CorrKS0 = corrcov(CovKS0);
fclose(KS0fileID);

[KS0coeff,KS0latent] = pcacov(CovKS0);

total_variance_sum = 0.0;
selected_variance_sum = 0.0;
for j = 1:N_total_bins
    total_variance_sum = total_variance_sum + KS0latent(j);
end
for j = 1:N_selected_vectors
    selected_variance_sum = selected_variance_sum + KS0latent(j);
end
selected_variance_sum / total_variance_sum

% print covariance all
KS0fileID_w = fopen('KS0_cov_all.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(KS0fileID_w,"%1.15f ",CovKS0(i,j));
    end
    fprintf(KS0fileID_w,"\n");
end
fclose(KS0fileID_w);

% print truncated covariance matrix all
KS0fileID_w = fopen('KS0_cov_all_truncated.txt','w');
for i = 1:N_total_bins
    fprintf(KS0fileID_w,"%1.15f\n",CovKS0(i,i));
end
fclose(KS0fileID_w);

% print correlation all
KS0fileID_w = fopen('KS0_corr_all.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(KS0fileID_w,"%1.15f ",CorrKS0(i,j));
    end
    fprintf(KS0fileID_w,"\n");
end
fclose(KS0fileID_w);

% print selected KS0 info into file
KS0fileID_w = fopen('KS0_selected.txt','w');
for j = 1:N_selected_vectors
    fprintf(KS0fileID_w,"%1.15f\n", sqrt(KS0latent(j)));
    for i = 1:N_total_bins
        temp_eigen_vector = KS0coeff(:,j);
        fprintf(KS0fileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
end
fclose(KS0fileID_w);

KS0partialcov = zeros(N_total_bins, N_total_bins);
for j = 1:N_selected_vectors
    KS0partialcov = KS0partialcov + KS0coeff(:,j) * KS0latent(j) * KS0coeff(:,j)';
end

KS0partialcorr = corrcov(KS0partialcov);

for j = 1:N_selected_vectors
    KS0coeff(:,j) = 0;
    KS0latent(j) = 0;
end

% print covariance partial
KS0fileID_w = fopen('KS0_cov_partial.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(KS0fileID_w,"%1.15f ",KS0partialcov(i,j));
    end
    fprintf(KS0fileID_w,"\n");
end
fclose(KS0fileID_w);

% print correlation partial
KS0fileID_w = fopen('KS0_corr_partial.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(KS0fileID_w,"%1.15f ",KS0partialcorr(i,j));
    end
    fprintf(KS0fileID_w,"\n");
end
fclose(KS0fileID_w);

% get covariance for remaining
KS0remaincov = KS0coeff*diag(KS0latent)*KS0coeff';

% print truncated remaining covariance matrix
KS0fileID_w = fopen('KS0_cov_remain_truncated.txt','w');
for i = 1:N_total_bins
    fprintf(KS0fileID_w,"%1.15f\n",KS0remaincov(i,i));
end
fclose(KS0fileID_w);

% latent: eigen value
% coeff: eigen vector