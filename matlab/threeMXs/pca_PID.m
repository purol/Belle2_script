clear all
format long

N_total_bins = 126;
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

KIDfileID_w = fopen('KID_MX1_selected.txt','w');
for j = 1:N_selected_vectors
    fprintf(KIDfileID_w,"%1.15f\n", sqrt(Kaonlatent(j)));
    for i = 1:6
        temp_eigen_vector = Kaoncoeff(:,j);
        fprintf(KIDfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 19:24
        temp_eigen_vector = Kaoncoeff(:,j);
        fprintf(KIDfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 37:42
        temp_eigen_vector = Kaoncoeff(:,j);
        fprintf(KIDfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 55:60
        temp_eigen_vector = Kaoncoeff(:,j);
        fprintf(KIDfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 73:78
        temp_eigen_vector = Kaoncoeff(:,j);
        fprintf(KIDfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 91:96
        temp_eigen_vector = Kaoncoeff(:,j);
        fprintf(KIDfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 109:114
        temp_eigen_vector = Kaoncoeff(:,j);
        fprintf(KIDfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
end
fclose(KIDfileID_w);

KIDfileID_w = fopen('KID_MX2_selected.txt','w');
for j = 1:N_selected_vectors
    fprintf(KIDfileID_w,"%1.15f\n", sqrt(Kaonlatent(j)));
    for i = 7:12
        temp_eigen_vector = Kaoncoeff(:,j);
        fprintf(KIDfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 25:30
        temp_eigen_vector = Kaoncoeff(:,j);
        fprintf(KIDfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 43:48
        temp_eigen_vector = Kaoncoeff(:,j);
        fprintf(KIDfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 61:66
        temp_eigen_vector = Kaoncoeff(:,j);
        fprintf(KIDfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 79:84
        temp_eigen_vector = Kaoncoeff(:,j);
        fprintf(KIDfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 97:102
        temp_eigen_vector = Kaoncoeff(:,j);
        fprintf(KIDfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 115:120
        temp_eigen_vector = Kaoncoeff(:,j);
        fprintf(KIDfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
end
fclose(KIDfileID_w);

KIDfileID_w = fopen('KID_MX3_selected.txt','w');
for j = 1:N_selected_vectors
    fprintf(KIDfileID_w,"%1.15f\n", sqrt(Kaonlatent(j)));
    for i = 13:18
        temp_eigen_vector = Kaoncoeff(:,j);
        fprintf(KIDfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 31:36
        temp_eigen_vector = Kaoncoeff(:,j);
        fprintf(KIDfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 49:54
        temp_eigen_vector = Kaoncoeff(:,j);
        fprintf(KIDfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 67:72
        temp_eigen_vector = Kaoncoeff(:,j);
        fprintf(KIDfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 85:90
        temp_eigen_vector = Kaoncoeff(:,j);
        fprintf(KIDfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 103:108
        temp_eigen_vector = Kaoncoeff(:,j);
        fprintf(KIDfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
    for i = 121:126
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

PIDfileID_w = fopen('PID_MX1_cov_remain_truncated.txt','w');
for i = 1:6
    fprintf(PIDfileID_w,"%1.15f\n",Pionremaincov(i,i));
end
for i = 19:24
    fprintf(PIDfileID_w,"%1.15f\n",Pionremaincov(i,i));
end
for i = 37:42
    fprintf(PIDfileID_w,"%1.15f\n",Pionremaincov(i,i));
end
for i = 55:60
    fprintf(PIDfileID_w,"%1.15f\n",Pionremaincov(i,i));
end
for i = 73:78
    fprintf(PIDfileID_w,"%1.15f\n",Pionremaincov(i,i));
end
for i = 91:96
    fprintf(PIDfileID_w,"%1.15f\n",Pionremaincov(i,i));
end
for i = 109:114
    fprintf(PIDfileID_w,"%1.15f\n",Pionremaincov(i,i));
end
fclose(PIDfileID_w);

PIDfileID_w = fopen('PID_MX2_cov_remain_truncated.txt','w');
for i = 7:12
    fprintf(PIDfileID_w,"%1.15f\n",Pionremaincov(i,i));
end
for i = 24:30
    fprintf(PIDfileID_w,"%1.15f\n",Pionremaincov(i,i));
end
for i = 43:48
    fprintf(PIDfileID_w,"%1.15f\n",Pionremaincov(i,i));
end
for i = 61:66
    fprintf(PIDfileID_w,"%1.15f\n",Pionremaincov(i,i));
end
for i = 79:84
    fprintf(PIDfileID_w,"%1.15f\n",Pionremaincov(i,i));
end
for i = 97:102
    fprintf(PIDfileID_w,"%1.15f\n",Pionremaincov(i,i));
end
for i = 115:120
    fprintf(PIDfileID_w,"%1.15f\n",Pionremaincov(i,i));
end
fclose(PIDfileID_w);

PIDfileID_w = fopen('PID_MX3_cov_remain_truncated.txt','w');
for i = 13:18
    fprintf(PIDfileID_w,"%1.15f\n",Pionremaincov(i,i));
end
for i = 30:36
    fprintf(PIDfileID_w,"%1.15f\n",Pionremaincov(i,i));
end
for i = 49:54
    fprintf(PIDfileID_w,"%1.15f\n",Pionremaincov(i,i));
end
for i = 67:72
    fprintf(PIDfileID_w,"%1.15f\n",Pionremaincov(i,i));
end
for i = 85:90
    fprintf(PIDfileID_w,"%1.15f\n",Pionremaincov(i,i));
end
for i = 91:96
    fprintf(PIDfileID_w,"%1.15f\n",Pionremaincov(i,i));
end
for i = 121:126
    fprintf(PIDfileID_w,"%1.15f\n",Pionremaincov(i,i));
end
fclose(PIDfileID_w);

% latent: eigen value
% coeff: eigen vector