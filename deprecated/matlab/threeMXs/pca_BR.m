clear all
format long

N_total_bins = 18;
N_selected_vectors = 5;

% BR
BRfileID = fopen('BR_toys.txt','r');
formatSpec = strjoin(repelem(["%f "], [N_total_bins]), '');
sizeBR = [N_total_bins Inf];
BRMatrix = fscanf(BRfileID,formatSpec,sizeBR);
BRMatrix = BRMatrix';
CovBR = cov(BRMatrix);
CorrBR = corrcov(CovBR);
fclose(BRfileID);

[BRcoeff,BRlatent] = pcacov(CovBR);

% print covariance all
BRfileID_w = fopen('BR_cov_all.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(BRfileID_w,"%1.15f ",CovBR(i,j));
    end
    fprintf(BRfileID_w,"\n");
end
fclose(BRfileID_w);

% print truncated covariance matrix all
BRfileID_w = fopen('BR_cov_all_truncated.txt','w');
for i = 1:N_total_bins
    fprintf(BRfileID_w,"%1.15f\n",CovBR(i,i));
end
fclose(BRfileID_w);

% print correlation all
BRfileID_w = fopen('BR_corr_all.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(BRfileID_w,"%1.15f ",CorrBR(i,j));
    end
    fprintf(BRfileID_w,"\n");
end
fclose(BRfileID_w);

% print selected KID info into file
BRfileID_w = fopen('BR_selected.txt','w');
for j = 1:N_selected_vectors
    fprintf(BRfileID_w,"%1.15f\n", sqrt(BRlatent(j)));
    for i = 1:N_total_bins
        temp_eigen_vector = BRcoeff(:,j);
        fprintf(BRfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
end
fclose(BRfileID_w);

BRfileID_w = fopen('BR_MX1_selected.txt','w');
for j = 1:N_selected_vectors
    fprintf(BRfileID_w,"%1.15f\n", sqrt(BRlatent(j)));
    for i = 1:6
        temp_eigen_vector = BRcoeff(:,j);
        fprintf(BRfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 19:24
        temp_eigen_vector = BRcoeff(:,j);
        fprintf(BRfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 37:42
        temp_eigen_vector = BRcoeff(:,j);
        fprintf(BRfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
end
fclose(BRfileID_w);

BRfileID_w = fopen('BR_MX2_selected.txt','w');
for j = 1:N_selected_vectors
    fprintf(BRfileID_w,"%1.15f\n", sqrt(BRlatent(j)));
    for i = 7:12
        temp_eigen_vector = BRcoeff(:,j);
        fprintf(BRfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 25:30
        temp_eigen_vector = BRcoeff(:,j);
        fprintf(BRfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 43:48
        temp_eigen_vector = BRcoeff(:,j);
        fprintf(BRfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
end
fclose(BRfileID_w);

BRfileID_w = fopen('BR_MX3_selected.txt','w');
for j = 1:N_selected_vectors
    fprintf(BRfileID_w,"%1.15f\n", sqrt(BRlatent(j)));
    for i = 13:18
        temp_eigen_vector = BRcoeff(:,j);
        fprintf(BRfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 31:36
        temp_eigen_vector = BRcoeff(:,j);
        fprintf(BRfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 49:54
        temp_eigen_vector = BRcoeff(:,j);
        fprintf(BRfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
end
fclose(BRfileID_w);

BRpartialcov = zeros(N_total_bins, N_total_bins);
for j = 1:N_selected_vectors
    BRpartialcov = BRpartialcov + BRcoeff(:,j) * BRlatent(j) * BRcoeff(:,j)';
end

BRpartialcorr = corrcov(BRpartialcov);

for j = 1:N_selected_vectors
    BRcoeff(:,j) = 0;
    BRlatent(j) = 0;
end

% print covariance partial
BRfileID_w = fopen('BR_cov_partial.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(BRfileID_w,"%1.15f ",BRpartialcov(i,j));
    end
    fprintf(BRfileID_w,"\n");
end
fclose(BRfileID_w);

% print correlation partial
BRfileID_w = fopen('BR_corr_partial.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(BRfileID_w,"%1.15f ",BRpartialcorr(i,j));
    end
    fprintf(BRfileID_w,"\n");
end
fclose(BRfileID_w);

% get covariance for remaining
BRremaincov = BRcoeff*diag(BRlatent)*BRcoeff';

% print truncated remaining covariance matrix
BRfileID_w = fopen('BR_cov_remain_truncated.txt','w');
for i = 1:N_total_bins
    fprintf(BRfileID_w,"%1.15f\n",BRremaincov(i,i));
end
fclose(BRfileID_w);

BRfileID_w = fopen('BR_MX1_cov_remain_truncated.txt','w');
for i = 1:6
    fprintf(BRfileID_w,"%1.15f\n",BRremaincov(i,i));
end
for i = 19:14
    fprintf(BRfileID_w,"%1.15f\n",BRremaincov(i,i));
end
for i = 37:42
    fprintf(BRfileID_w,"%1.15f\n",BRremaincov(i,i));
end
fclose(BRfileID_w);

BRfileID_w = fopen('BR_MX2_cov_remain_truncated.txt','w');
for i = 7:12
    fprintf(BRfileID_w,"%1.15f\n",BRremaincov(i,i));
end
for i = 25:30
    fprintf(BRfileID_w,"%1.15f\n",BRremaincov(i,i));
end
for i = 43:48
    fprintf(BRfileID_w,"%1.15f\n",BRremaincov(i,i));
end
fclose(BRfileID_w);

BRfileID_w = fopen('BR_MX3_cov_remain_truncated.txt','w');
for i = 13:18
    fprintf(BRfileID_w,"%1.15f\n",BRremaincov(i,i));
end
for i = 31:36
    fprintf(BRfileID_w,"%1.15f\n",BRremaincov(i,i));
end
for i = 49:54
    fprintf(BRfileID_w,"%1.15f\n",BRremaincov(i,i));
end
fclose(BRfileID_w);

% latent: eigen value
% coeff: eigen vector