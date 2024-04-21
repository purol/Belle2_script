clear all
format long

N_total_bins = 22 * 3;
N_selected_vectors = 2;

% Kstarff
KstarfffileID = fopen('Kstarff_toys.txt','r');
formatSpec = strjoin(repelem(["%f "], [N_total_bins]), '');
sizeKstarff = [N_total_bins Inf];
KstarffMatrix = fscanf(KstarfffileID,formatSpec,sizeKstarff);
KstarffMatrix = KstarffMatrix';
CovKstarff = cov(KstarffMatrix);
CovKstarff(isnan(CovKstarff))=0;
CorrKstarff = corrcov(CovKstarff);
CorrKstarff(isnan(CorrKstarff))=0;
fclose(KstarfffileID);

[Kstarffcoeff,Kstarfflatent] = pcacov(CovKstarff);

% print covariance all
KstarfffileID_w = fopen('Kstarff_cov_all.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(KstarfffileID_w,"%1.15f ",CovKstarff(i,j));
    end
    fprintf(KstarfffileID_w,"\n");
end
fclose(KstarfffileID_w);

% print truncated covariance matrix all
KstarfffileID_w = fopen('Kstarff_cov_all_truncated.txt','w');
for i = 1:N_total_bins
    fprintf(KstarfffileID_w,"%1.15f\n",CovKstarff(i,i));
end
fclose(KstarfffileID_w);

% print correlation all
KstarfffileID_w = fopen('Kstarff_corr_all.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(KstarfffileID_w,"%1.15f ",CorrKstarff(i,j));
    end
    fprintf(KstarfffileID_w,"\n");
end
fclose(KstarfffileID_w);

% print selected KID info into file
KstarfffileID_w = fopen('Kstarff_selected.txt','w');
for j = 1:N_selected_vectors
    fprintf(KstarfffileID_w,"%1.15f\n", sqrt(Kstarfflatent(j)));
    for i = 1:N_total_bins
        temp_eigen_vector = Kstarffcoeff(:,j);
        fprintf(KstarfffileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
end
fclose(KstarfffileID_w);

Kstarffpartialcov = zeros(N_total_bins, N_total_bins);
for j = 1:N_selected_vectors
    Kstarffpartialcov = Kstarffpartialcov + Kstarffcoeff(:,j) * Kstarfflatent(j) * Kstarffcoeff(:,j)';
end

Kstarffpartialcorr = corrcov(Kstarffpartialcov);
Kstarffpartialcov(isnan(Kstarffpartialcov))=0;

for j = 1:N_selected_vectors
    Kstarffcoeff(:,j) = 0;
    Kstarfflatent(j) = 0;
end

% print covariance partial
KstarfffileID_w = fopen('Kstarff_cov_partial.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(KstarfffileID_w,"%1.15f ",Kstarffpartialcov(i,j));
    end
    fprintf(KstarfffileID_w,"\n");
end
fclose(KstarfffileID_w);

% print correlation partial
KstarfffileID_w = fopen('Kstarff_corr_partial.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(KstarfffileID_w,"%1.15f ",Kstarffpartialcorr(i,j));
    end
    fprintf(KstarfffileID_w,"\n");
end
fclose(KstarfffileID_w);

% get covariance for remaining
Kstarffremaincov = Kstarffcoeff*diag(Kstarfflatent)*Kstarffcoeff';

% print truncated remaining covariance matrix
KstarfffileID_w = fopen('Kstarff_cov_remain_truncated.txt','w');
for i = 1:N_total_bins
    fprintf(KstarfffileID_w,"%1.15f\n",Kstarffremaincov(i,i));
end
fclose(KstarfffileID_w);

% latent: eigen value
% coeff: eigen vector