clear all
format long

N_total_bins = 22 * 7;
N_selected_vectors = 15;

% Kaon ID
KIDfileID = fopen('KID_toys.txt','r');
formatSpec = strjoin(repelem(["%f "], [N_total_bins]), '');
sizeK = [N_total_bins Inf];
KaonMatrix = fscanf(KIDfileID,formatSpec,sizeK);
KaonMatrix = KaonMatrix';
CovKaonID = cov(KaonMatrix);
CorrKaonID = corrcov(CovKaonID);
fclose(KIDfileID);

[Kaoncoeff,Kaonlatent] = pcacov(CovKaonID);

% print covariance all
KIDfileID_w = fopen('KID_cov_all.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(KIDfileID_w,"%1.15f ",CovKaonID(i,j));
    end
    fprintf(KIDfileID_w,"\n");
end
fclose(KIDfileID_w);

% print truncated covariance matrix all
KIDfileID_w = fopen('KID_cov_all_truncated.txt','w');
for i = 1:N_total_bins
    fprintf(KIDfileID_w,"%1.15f\n",CovKaonID(i,i));
end
fclose(KIDfileID_w);

% print correlation all
KIDfileID_w = fopen('KID_corr_all.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(KIDfileID_w,"%1.15f ",CorrKaonID(i,j));
    end
    fprintf(KIDfileID_w,"\n");
end
fclose(KIDfileID_w);

% print selected KID info into file
KIDfileID_w = fopen('KID_selected.txt','w');
for j = 1:N_selected_vectors
    fprintf(KIDfileID_w,"%1.15f\n", sqrt(Kaonlatent(j)));
    for i = 1:N_total_bins
        temp_eigen_vector = Kaoncoeff(:,j);
        fprintf(KIDfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
end
fclose(KIDfileID_w);

Kaonpartialcov = zeros(N_total_bins, N_total_bins);
for j = 1:N_selected_vectors
    Kaonpartialcov = Kaonpartialcov + Kaoncoeff(:,j) * Kaonlatent(j) * Kaoncoeff(:,j)';
end

Kaonpartialcorr = corrcov(Kaonpartialcov);

for j = 1:N_selected_vectors
    Kaoncoeff(:,j) = 0;
    Kaonlatent(j) = 0;
end

% print covariance partial
KIDfileID_w = fopen('KID_cov_partial.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(KIDfileID_w,"%1.15f ",Kaonpartialcov(i,j));
    end
    fprintf(KIDfileID_w,"\n");
end
fclose(KIDfileID_w);

% print correlation partial
KIDfileID_w = fopen('KID_corr_partial.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(KIDfileID_w,"%1.15f ",Kaonpartialcorr(i,j));
    end
    fprintf(KIDfileID_w,"\n");
end
fclose(KIDfileID_w);

% get covariance for remaining
Kaonremaincov = Kaoncoeff*diag(Kaonlatent)*Kaoncoeff';

% print truncated remaining covariance matrix
KIDfileID_w = fopen('KID_cov_remain_truncated.txt','w');
for i = 1:N_total_bins
    fprintf(KIDfileID_w,"%1.15f\n",Kaonremaincov(i,i));
end
fclose(KIDfileID_w);


% Pion ID
PIDfileID = fopen('PID_toys.txt','r');
formatSpec = strjoin(repelem(["%f "], [N_total_bins]), '');
sizeK = [N_total_bins Inf];
PionMatrix = fscanf(PIDfileID,formatSpec,sizeK);
PionMatrix = PionMatrix';
CovPionID = cov(PionMatrix);
CorrPionID = corrcov(CovPionID);
fclose(PIDfileID);

[Pioncoeff,Pionlatent] = pcacov(CovPionID);

% print covariance all
PIDfileID_w = fopen('PID_cov_all.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(PIDfileID_w,"%1.15f ",CovPionID(i,j));
    end
    fprintf(PIDfileID_w,"\n");
end
fclose(PIDfileID_w);

% print truncated covariance all
PIDfileID_w = fopen('PID_cov_all_truncated.txt','w');
for i = 1:N_total_bins
    fprintf(PIDfileID_w,"%1.15f\n",CovPionID(i,i));
end
fclose(PIDfileID_w);

% print correlation all
PIDfileID_w = fopen('PID_corr_all.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(PIDfileID_w,"%1.15f ",CorrPionID(i,j));
    end
    fprintf(PIDfileID_w,"\n");
end
fclose(PIDfileID_w);

% print selected PID info into file
PIDfileID_w = fopen('PID_selected.txt','w');
for j = 1:N_selected_vectors
    fprintf(PIDfileID_w,"%1.15f\n", sqrt(Pionlatent(j)));
    for i = 1:N_total_bins
        temp_eigen_vector = Pioncoeff(:,j);
        fprintf(PIDfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
end
fclose(PIDfileID_w);

Pionpartialcov = zeros(N_total_bins, N_total_bins);
for j = 1:N_selected_vectors
    Pionpartialcov = Pionpartialcov + Pioncoeff(:,j) * Pionlatent(j) * Pioncoeff(:,j)';
end

Pionpartialcorr = corrcov(Pionpartialcov);

for j = 1:N_selected_vectors
    Pioncoeff(:,j) = 0;
    Pionlatent(j) = 0;
end

% print covariance partial
PIDfileID_w = fopen('PID_cov_partial.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(PIDfileID_w,"%1.15f ",Pionpartialcov(i,j));
    end
    fprintf(PIDfileID_w,"\n");
end
fclose(PIDfileID_w);

% print correlation partial
PIDfileID_w = fopen('PID_corr_partial.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(PIDfileID_w,"%1.15f ",Pionpartialcorr(i,j));
    end
    fprintf(PIDfileID_w,"\n");
end
fclose(PIDfileID_w);

% get covariance for remaining
Pionremaincov = Pioncoeff*diag(Pionlatent)*Pioncoeff';

% print truncated remaining covariance matrix
PIDfileID_w = fopen('PID_cov_remain_truncated.txt','w');
for i = 1:N_total_bins
    fprintf(PIDfileID_w,"%1.15f\n",Pionremaincov(i,i));
end
fclose(PIDfileID_w);

% latent: eigen value
% coeff: eigen vector