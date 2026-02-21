clear all
format long

N_total_bins = 42;
N_selected_vectors = 3;

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

pi0fileID_w = fopen('pi0_MX1_selected.txt','w');
for j = 1:N_selected_vectors
    fprintf(pi0fileID_w,"%1.15f\n", sqrt(pi0latent(j)));
    for i = 1:6
        temp_eigen_vector = pi0coeff(:,j);
        fprintf(pi0fileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 19:24
        temp_eigen_vector = pi0coeff(:,j);
        fprintf(pi0fileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 37:42
        temp_eigen_vector = pi0coeff(:,j);
        fprintf(pi0fileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 55:60
        temp_eigen_vector = pi0coeff(:,j);
        fprintf(pi0fileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 73:78
        temp_eigen_vector = pi0coeff(:,j);
        fprintf(pi0fileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 91:96
        temp_eigen_vector = pi0coeff(:,j);
        fprintf(pi0fileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 109:114
        temp_eigen_vector = pi0coeff(:,j);
        fprintf(pi0fileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
end
fclose(pi0fileID_w);

pi0fileID_w = fopen('pi0_MX2_selected.txt','w');
for j = 1:N_selected_vectors
    fprintf(pi0fileID_w,"%1.15f\n", sqrt(pi0latent(j)));
    for i = 7:12
        temp_eigen_vector = pi0coeff(:,j);
        fprintf(pi0fileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 25:30
        temp_eigen_vector = pi0coeff(:,j);
        fprintf(pi0fileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 43:48
        temp_eigen_vector = pi0coeff(:,j);
        fprintf(pi0fileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 61:66
        temp_eigen_vector = pi0coeff(:,j);
        fprintf(pi0fileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 79:84
        temp_eigen_vector = pi0coeff(:,j);
        fprintf(pi0fileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 97:102
        temp_eigen_vector = pi0coeff(:,j);
        fprintf(pi0fileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 115:120
        temp_eigen_vector = pi0coeff(:,j);
        fprintf(pi0fileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
end
fclose(pi0fileID_w);

pi0fileID_w = fopen('pi0_MX3_selected.txt','w');
for j = 1:N_selected_vectors
    fprintf(pi0fileID_w,"%1.15f\n", sqrt(pi0latent(j)));
    for i = 13:18
        temp_eigen_vector = pi0coeff(:,j);
        fprintf(pi0fileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 31:36
        temp_eigen_vector = pi0coeff(:,j);
        fprintf(pi0fileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 49:54
        temp_eigen_vector = pi0coeff(:,j);
        fprintf(pi0fileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 67:72
        temp_eigen_vector = pi0coeff(:,j);
        fprintf(pi0fileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 85:90
        temp_eigen_vector = pi0coeff(:,j);
        fprintf(pi0fileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 103:108
        temp_eigen_vector = pi0coeff(:,j);
        fprintf(pi0fileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 121:126
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

pi0fileID_w = fopen('pi0_MX1_cov_remain_truncated.txt','w');
for i = 1:6
    fprintf(pi0fileID_w,"%1.15f\n",pi0remaincov(i,i));
end
for i = 19:24
    fprintf(pi0fileID_w,"%1.15f\n",pi0remaincov(i,i));
end
for i = 37:42
    fprintf(pi0fileID_w,"%1.15f\n",pi0remaincov(i,i));
end
for i = 55:60
    fprintf(pi0fileID_w,"%1.15f\n",pi0remaincov(i,i));
end
for i = 73:78
    fprintf(pi0fileID_w,"%1.15f\n",pi0remaincov(i,i));
end
for i = 91:96
    fprintf(pi0fileID_w,"%1.15f\n",pi0remaincov(i,i));
end
for i = 109:114
    fprintf(pi0fileID_w,"%1.15f\n",pi0remaincov(i,i));
end
fclose(pi0fileID_w);

pi0fileID_w = fopen('pi0_MX2_cov_remain_truncated.txt','w');
for i = 7:12
    fprintf(pi0fileID_w,"%1.15f\n",pi0remaincov(i,i));
end
for i = 25:30
    fprintf(pi0fileID_w,"%1.15f\n",pi0remaincov(i,i));
end
for i = 43:48
    fprintf(pi0fileID_w,"%1.15f\n",pi0remaincov(i,i));
end
for i = 61:66
    fprintf(pi0fileID_w,"%1.15f\n",pi0remaincov(i,i));
end
for i = 79:84
    fprintf(pi0fileID_w,"%1.15f\n",pi0remaincov(i,i));
end
for i = 97:102
    fprintf(pi0fileID_w,"%1.15f\n",pi0remaincov(i,i));
end
for i = 115:120
    fprintf(pi0fileID_w,"%1.15f\n",pi0remaincov(i,i));
end
fclose(pi0fileID_w);

pi0fileID_w = fopen('pi0_MX3_cov_remain_truncated.txt','w');
for i = 13:18
    fprintf(pi0fileID_w,"%1.15f\n",pi0remaincov(i,i));
end
for i = 31:36
    fprintf(pi0fileID_w,"%1.15f\n",pi0remaincov(i,i));
end
for i = 49:54
    fprintf(pi0fileID_w,"%1.15f\n",pi0remaincov(i,i));
end
for i = 67:72
    fprintf(pi0fileID_w,"%1.15f\n",pi0remaincov(i,i));
end
for i = 85:90
    fprintf(pi0fileID_w,"%1.15f\n",pi0remaincov(i,i));
end
for i = 103:108
    fprintf(pi0fileID_w,"%1.15f\n",pi0remaincov(i,i));
end
for i = 121:126
    fprintf(pi0fileID_w,"%1.15f\n",pi0remaincov(i,i));
end
fclose(pi0fileID_w);

% latent: eigen value
% coeff: eigen vector