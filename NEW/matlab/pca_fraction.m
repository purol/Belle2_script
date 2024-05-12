clear all
format long

N_total_bins = 22 * 3;
N_selected_vectors = 2;

% fraction
fractionfileID = fopen('fraction_toys.txt','r');
formatSpec = strjoin(repelem(["%f "], [N_total_bins]), '');
sizefraction = [N_total_bins Inf];
fractionMatrix = fscanf(fractionfileID,formatSpec,sizefraction);
fractionMatrix = fractionMatrix';
Covfraction = cov(fractionMatrix);
Covfraction(isnan(Covfraction))=0;
Corrfraction = corrcov(Covfraction);
Corrfraction(isnan(Corrfraction))=0;
fclose(fractionfileID);

[fractioncoeff,fractionlatent] = pcacov(Covfraction);

% print covariance all
fractionfileID_w = fopen('fraction_cov_all.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(fractionfileID_w,"%1.15f ",Covfraction(i,j));
    end
    fprintf(fractionfileID_w,"\n");
end
fclose(fractionfileID_w);

% print truncated covariance matrix all
fractionfileID_w = fopen('fraction_cov_all_truncated.txt','w');
for i = 1:N_total_bins
    fprintf(fractionfileID_w,"%1.15f\n",Covfraction(i,i));
end
fclose(fractionfileID_w);

% print correlation all
fractionfileID_w = fopen('fraction_corr_all.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(fractionfileID_w,"%1.15f ",Corrfraction(i,j));
    end
    fprintf(fractionfileID_w,"\n");
end
fclose(fractionfileID_w);

% print selected KID info into file
fractionfileID_w = fopen('fraction_selected.txt','w');
for j = 1:N_selected_vectors
    fprintf(fractionfileID_w,"%1.15f\n", sqrt(fractionlatent(j)));
    for i = 1:N_total_bins
        temp_eigen_vector = fractioncoeff(:,j);
        fprintf(fractionfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
end
fclose(fractionfileID_w);

fractionpartialcov = zeros(N_total_bins, N_total_bins);
for j = 1:N_selected_vectors
    fractionpartialcov = fractionpartialcov + fractioncoeff(:,j) * fractionlatent(j) * fractioncoeff(:,j)';
end

fractionpartialcorr = corrcov(fractionpartialcov);
fractionpartialcorr(isnan(fractionpartialcorr))=0;

for j = 1:N_selected_vectors
    fractioncoeff(:,j) = 0;
    fractionlatent(j) = 0;
end

% print covariance partial
fractionfileID_w = fopen('fraction_cov_partial.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(fractionfileID_w,"%1.15f ",fractionpartialcov(i,j));
    end
    fprintf(fractionfileID_w,"\n");
end
fclose(fractionfileID_w);

% print correlation partial
fractionfileID_w = fopen('fraction_corr_partial.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(fractionfileID_w,"%1.15f ",fractionpartialcorr(i,j));
    end
    fprintf(fractionfileID_w,"\n");
end
fclose(fractionfileID_w);

% get covariance for remaining
fractionremaincov = fractioncoeff*diag(fractionlatent)*fractioncoeff';

% print truncated remaining covariance matrix
fractionfileID_w = fopen('fraction_cov_remain_truncated.txt','w');
for i = 1:N_total_bins
    fprintf(fractionfileID_w,"%1.15f\n",fractionremaincov(i,i));
end
fclose(fractionfileID_w);

% latent: eigen value
% coeff: eigen vector