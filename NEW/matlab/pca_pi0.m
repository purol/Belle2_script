clear all
format long

N_total_bins = 22 * 7;
N_selected_vectors = 1;

% pi0 correction
pi0fileID = fopen('pi0_toys.txt','r');
formatSpec = strjoin(repelem(["%f "], [N_total_bins]), '');
sizepi0 = [N_total_bins Inf];
pi0Matrix = fscanf(pi0fileID,formatSpec,sizepi0);
pi0Matrix = pi0Matrix';
Covpi0 = cov(pi0Matrix);
Corrpi0 = corrcov(Covpi0);
fclose(pi0fileID);

[pi0coeff,pi0latent] = pcacov(Covpi0);

% print covariance all
pi0fileID_w = fopen('pi0_cov_all.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(pi0fileID_w,"%1.15f ",Covpi0(i,j));
    end
    fprintf(pi0fileID_w,"\n");
end
fclose(pi0fileID_w);

% print truncated covariance matrix all
pi0fileID_w = fopen('pi0_cov_all_truncated.txt','w');
for i = 1:N_total_bins
    fprintf(pi0fileID_w,"%1.15f\n",Covpi0(i,i));
end
fclose(pi0fileID_w);

% print correlation all
pi0fileID_w = fopen('pi0_corr_all.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(pi0fileID_w,"%1.15f ",Corrpi0(i,j));
    end
    fprintf(pi0fileID_w,"\n");
end
fclose(pi0fileID_w);

% print selected pi0 info into file
pi0fileID_w = fopen('pi0_selected.txt','w');
for j = 1:N_selected_vectors
    fprintf(pi0fileID_w,"%1.15f\n", sqrt(pi0latent(j)));
    for i = 1:N_total_bins
        temp_eigen_vector = pi0coeff(:,j);
        fprintf(pi0fileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
end
fclose(pi0fileID_w);

pi0partialcov = zeros(N_total_bins, N_total_bins);
for j = 1:N_selected_vectors
    pi0partialcov = pi0partialcov + pi0coeff(:,j) * pi0latent(j) * pi0coeff(:,j)';
end

pi0partialcorr = corrcov(pi0partialcov);

for j = 1:N_selected_vectors
    pi0coeff(:,j) = 0;
    pi0latent(j) = 0;
end

% print covariance partial
pi0fileID_w = fopen('pi0_cov_partial.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(pi0fileID_w,"%1.15f ",pi0partialcov(i,j));
    end
    fprintf(pi0fileID_w,"\n");
end
fclose(pi0fileID_w);

% print correlation partial
pi0fileID_w = fopen('pi0_corr_partial.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(pi0fileID_w,"%1.15f ",pi0partialcorr(i,j));
    end
    fprintf(pi0fileID_w,"\n");
end
fclose(pi0fileID_w);

% get covariance for remaining
pi0remaincov = pi0coeff*diag(pi0latent)*pi0coeff';

% print truncated remaining covariance matrix
pi0fileID_w = fopen('pi0_cov_remain_truncated.txt','w');
for i = 1:N_total_bins
    fprintf(pi0fileID_w,"%1.15f\n",pi0remaincov(i,i));
end
fclose(pi0fileID_w);

% latent: eigen value
% coeff: eigen vector