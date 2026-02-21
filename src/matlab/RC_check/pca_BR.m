clear all
format long

N_total_bins = 15 * 5;
N_selected_vectors = 7;

% Generate base directory names
base_dirs = arrayfun(@(x) sprintf('KumoiRD_v%03d', x), 11:110, 'UniformOutput', false);

for d = 1:length(base_dirs)
    base_dir = base_dirs{d}; 
% BR
BRfileID = fopen(fullfile(base_dir, 'BR_toys.txt'),'r');
formatSpec = strjoin(repelem(["%f "], [N_total_bins]), '');
sizeBR = [N_total_bins Inf];
BRMatrix = fscanf(BRfileID,formatSpec,sizeBR);
BRMatrix = BRMatrix';
CovBR = cov(BRMatrix);
CorrBR = corrcov(CovBR);
fclose(BRfileID);

[BRcoeff,BRlatent] = pcacov(CovBR);

total_variance_sum = 0.0;
selected_variance_sum = 0.0;
for j = 1:N_total_bins
    total_variance_sum = total_variance_sum + BRlatent(j);
end
for j = 1:N_selected_vectors
    selected_variance_sum = selected_variance_sum + BRlatent(j);
end
selected_variance_sum / total_variance_sum

% print covariance all
BRfileID_w = fopen(fullfile(base_dir, 'BR_cov_all.txt'),'w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(BRfileID_w,"%1.15f ",CovBR(i,j));
    end
    fprintf(BRfileID_w,"\n");
end
fclose(BRfileID_w);

% print truncated covariance matrix all
BRfileID_w = fopen(fullfile(base_dir, 'BR_cov_all_truncated.txt'),'w');
for i = 1:N_total_bins
    fprintf(BRfileID_w,"%1.15f\n",CovBR(i,i));
end
fclose(BRfileID_w);

% print correlation all
BRfileID_w = fopen(fullfile(base_dir, 'BR_corr_all.txt'),'w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(BRfileID_w,"%1.15f ",CorrBR(i,j));
    end
    fprintf(BRfileID_w,"\n");
end
fclose(BRfileID_w);

% print selected KID info into file
BRfileID_w = fopen(fullfile(base_dir, 'BR_selected.txt'),'w');
for j = 1:N_selected_vectors
    fprintf(BRfileID_w,"%1.15f\n", sqrt(BRlatent(j)));
    for i = 1:N_total_bins
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
BRfileID_w = fopen(fullfile(base_dir, 'BR_cov_partial.txt'),'w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(BRfileID_w,"%1.15f ",BRpartialcov(i,j));
    end
    fprintf(BRfileID_w,"\n");
end
fclose(BRfileID_w);

% print correlation partial
BRfileID_w = fopen(fullfile(base_dir, 'BR_corr_partial.txt'),'w');
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
BRfileID_w = fopen(fullfile(base_dir, 'BR_cov_remain_truncated.txt'),'w');
for i = 1:N_total_bins
    fprintf(BRfileID_w,"%1.15f\n",BRremaincov(i,i));
end
fclose(BRfileID_w);
end
% latent: eigen value
% coeff: eigen vector