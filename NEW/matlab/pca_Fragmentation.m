clear all
format long

N_total_bins = 15 * 5;
N_selected_vectors = 15;

% Fragmentation
FragmentationfileID = fopen('Fragmentation_toys.txt','r');
formatSpec = strjoin(repelem(["%f "], [N_total_bins]), '');
sizeFragmentation = [N_total_bins Inf];
FragmentationMatrix = fscanf(FragmentationfileID,formatSpec,sizeFragmentation);
FragmentationMatrix = FragmentationMatrix';
CovFragmentation = cov(FragmentationMatrix);
CovFragmentation(isnan(CovFragmentation))=0;
CorrFragmentation = corrcov(CovFragmentation);
CorrFragmentation(isnan(CorrFragmentation))=0;
fclose(FragmentationfileID);

[Fragmentationcoeff,Fragmentationlatent] = pcacov(CovFragmentation);

total_variance_sum = 0.0;
selected_variance_sum = 0.0;
for j = 1:N_total_bins
    total_variance_sum = total_variance_sum + Fragmentationlatent(j);
end
for j = 1:N_selected_vectors
    selected_variance_sum = selected_variance_sum + Fragmentationlatent(j);
end
selected_variance_sum / total_variance_sum

% print covariance all
FragmentationfileID_w = fopen('Fragmentation_cov_all.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(FragmentationfileID_w,"%1.15f ",CovFragmentation(i,j));
    end
    fprintf(FragmentationfileID_w,"\n");
end
fclose(FragmentationfileID_w);

% print truncated covariance matrix all
FragmentationfileID_w = fopen('Fragmentation_cov_all_truncated.txt','w');
for i = 1:N_total_bins
    fprintf(FragmentationfileID_w,"%1.15f\n",CovFragmentation(i,i));
end
fclose(FragmentationfileID_w);

% print correlation all
FragmentationfileID_w = fopen('Fragmentation_corr_all.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(FragmentationfileID_w,"%1.15f ",CorrFragmentation(i,j));
    end
    fprintf(FragmentationfileID_w,"\n");
end
fclose(FragmentationfileID_w);

% print selected KID info into file
FragmentationfileID_w = fopen('Fragmentation_selected.txt','w');
for j = 1:N_selected_vectors
    fprintf(FragmentationfileID_w,"%1.15f\n", sqrt(Fragmentationlatent(j)));
    for i = 1:N_total_bins
        temp_eigen_vector = Fragmentationcoeff(:,j);
        fprintf(FragmentationfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
end
fclose(FragmentationfileID_w);

Fragmentationpartialcov = zeros(N_total_bins, N_total_bins);
for j = 1:N_selected_vectors
    Fragmentationpartialcov = Fragmentationpartialcov + Fragmentationcoeff(:,j) * Fragmentationlatent(j) * Fragmentationcoeff(:,j)';
end

Fragmentationpartialcorr = corrcov(Fragmentationpartialcov);
Fragmentationpartialcorr(isnan(Fragmentationpartialcorr))=0;

for j = 1:N_selected_vectors
    Fragmentationcoeff(:,j) = 0;
    Fragmentationlatent(j) = 0;
end

% print covariance partial
FragmentationfileID_w = fopen('Fragmentation_cov_partial.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(FragmentationfileID_w,"%1.15f ",Fragmentationpartialcov(i,j));
    end
    fprintf(FragmentationfileID_w,"\n");
end
fclose(FragmentationfileID_w);

% print correlation partial
FragmentationfileID_w = fopen('Fragmentation_corr_partial.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(FragmentationfileID_w,"%1.15f ",Fragmentationpartialcorr(i,j));
    end
    fprintf(FragmentationfileID_w,"\n");
end
fclose(FragmentationfileID_w);

% get covariance for remaining
Fragmentationremaincov = Fragmentationcoeff*diag(Fragmentationlatent)*Fragmentationcoeff';

% print truncated remaining covariance matrix
FragmentationfileID_w = fopen('Fragmentation_cov_remain_truncated.txt','w');
for i = 1:N_total_bins
    fprintf(FragmentationfileID_w,"%1.15f\n",Fragmentationremaincov(i,i));
end
fclose(FragmentationfileID_w);

% latent: eigen value
% coeff: eigen vector