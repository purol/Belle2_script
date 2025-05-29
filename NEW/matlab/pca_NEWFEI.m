clear all
format long

N_total_bins = 15 * 3;
N_selected_vectors = 4;

% NEWFEICAL
NEWFEICALfileID = fopen('NEWFEICAL_toys.txt','r');
formatSpec = strjoin(repelem(["%f "], [N_total_bins]), '');
sizeNEWFEICAL = [N_total_bins Inf];
NEWFEICALMatrix = fscanf(NEWFEICALfileID,formatSpec,sizeNEWFEICAL);
NEWFEICALMatrix = NEWFEICALMatrix';
CovNEWFEICAL = cov(NEWFEICALMatrix);
CovNEWFEICAL(isnan(CovNEWFEICAL))=0;
CorrNEWFEICAL = corrcov(CovNEWFEICAL);
CorrNEWFEICAL(isnan(CorrNEWFEICAL))=0;
fclose(NEWFEICALfileID);

[NEWFEICALcoeff,NEWFEICALlatent] = pcacov(CovNEWFEICAL);

total_variance_sum = 0.0;
selected_variance_sum = 0.0;
for j = 1:N_total_bins
    total_variance_sum = total_variance_sum + NEWFEICALlatent(j);
end
for j = 1:N_selected_vectors
    selected_variance_sum = selected_variance_sum + NEWFEICALlatent(j);
end
selected_variance_sum / total_variance_sum

% print covariance all
NEWFEICALfileID_w = fopen('NEWFEICAL_cov_all.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(NEWFEICALfileID_w,"%1.15f ",CovNEWFEICAL(i,j));
    end
    fprintf(NEWFEICALfileID_w,"\n");
end
fclose(NEWFEICALfileID_w);

% print truncated covariance matrix all
NEWFEICALfileID_w = fopen('NEWFEICAL_cov_all_truncated.txt','w');
for i = 1:N_total_bins
    fprintf(NEWFEICALfileID_w,"%1.15f\n",CovNEWFEICAL(i,i));
end
fclose(NEWFEICALfileID_w);

% print correlation all
NEWFEICALfileID_w = fopen('NEWFEICAL_corr_all.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(NEWFEICALfileID_w,"%1.15f ",CorrNEWFEICAL(i,j));
    end
    fprintf(NEWFEICALfileID_w,"\n");
end
fclose(NEWFEICALfileID_w);

% print selected KID info into file
NEWFEICALfileID_w = fopen('NEWFEICAL_selected.txt','w');
for j = 1:N_selected_vectors
    fprintf(NEWFEICALfileID_w,"%1.15f\n", sqrt(NEWFEICALlatent(j)));
    for i = 1:N_total_bins
        temp_eigen_vector = NEWFEICALcoeff(:,j);
        fprintf(NEWFEICALfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
end
fclose(NEWFEICALfileID_w);

NEWFEICALpartialcov = zeros(N_total_bins, N_total_bins);
for j = 1:N_selected_vectors
    NEWFEICALpartialcov = NEWFEICALpartialcov + NEWFEICALcoeff(:,j) * NEWFEICALlatent(j) * NEWFEICALcoeff(:,j)';
end

NEWFEICALpartialcorr = corrcov(NEWFEICALpartialcov);
NEWFEICALpartialcorr(isnan(NEWFEICALpartialcorr))=0;

for j = 1:N_selected_vectors
    NEWFEICALcoeff(:,j) = 0;
    NEWFEICALlatent(j) = 0;
end

% print covariance partial
NEWFEICALfileID_w = fopen('NEWFEICAL_cov_partial.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(NEWFEICALfileID_w,"%1.15f ",NEWFEICALpartialcov(i,j));
    end
    fprintf(NEWFEICALfileID_w,"\n");
end
fclose(NEWFEICALfileID_w);

% print correlation partial
NEWFEICALfileID_w = fopen('NEWFEICAL_corr_partial.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(NEWFEICALfileID_w,"%1.15f ",NEWFEICALpartialcorr(i,j));
    end
    fprintf(NEWFEICALfileID_w,"\n");
end
fclose(NEWFEICALfileID_w);

% get covariance for remaining
NEWFEICALremaincov = NEWFEICALcoeff*diag(NEWFEICALlatent)*NEWFEICALcoeff';

% print truncated remaining covariance matrix
NEWFEICALfileID_w = fopen('NEWFEICAL_cov_remain_truncated.txt','w');
for i = 1:N_total_bins
    fprintf(NEWFEICALfileID_w,"%1.15f\n",NEWFEICALremaincov(i,i));
end
fclose(NEWFEICALfileID_w);

% latent: eigen value
% coeff: eigen vector