clear all
format long

N_total_bins = 18 * 7;
N_selected_vectors = 2;

% multiplicity correction
multiplicityfileID = fopen('multiplicity_toys.txt','r');
formatSpec = strjoin(repelem(["%f "], [N_total_bins]), '');
sizemultiplicity = [N_total_bins Inf];
multiplicityMatrix = fscanf(multiplicityfileID,formatSpec,sizemultiplicity);
multiplicityMatrix = multiplicityMatrix';
Covmultiplicity = cov(multiplicityMatrix);
Corrmultiplicity = corrcov(Covmultiplicity);
fclose(multiplicityfileID);

[multiplicitycoeff,multiplicitylatent] = pcacov(Covmultiplicity);

% print covariance all
multiplicityfileID_w = fopen('multiplicity_cov_all.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(multiplicityfileID_w,"%1.15f ",Covmultiplicity(i,j));
    end
    fprintf(multiplicityfileID_w,"\n");
end
fclose(multiplicityfileID_w);

% print truncated covariance matrix all
multiplicityfileID_w = fopen('multiplicity_cov_all_truncated.txt','w');
for i = 1:N_total_bins
    fprintf(multiplicityfileID_w,"%1.15f\n",Covmultiplicity(i,i));
end
fclose(multiplicityfileID_w);

% print correlation all
multiplicityfileID_w = fopen('multiplicity_corr_all.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(multiplicityfileID_w,"%1.15f ",Corrmultiplicity(i,j));
    end
    fprintf(multiplicityfileID_w,"\n");
end
fclose(multiplicityfileID_w);

% print selected multiplicity info into file
multiplicityfileID_w = fopen('multiplicity_selected.txt','w');
for j = 1:N_selected_vectors
    fprintf(multiplicityfileID_w,"%1.15f\n", sqrt(multiplicitylatent(j)));
    for i = 1:N_total_bins
        temp_eigen_vector = multiplicitycoeff(:,j);
        fprintf(multiplicityfileID_w,"%1.15f\n", temp_eigen_vector(i));
    end
end
fclose(multiplicityfileID_w);

multiplicitypartialcov = zeros(N_total_bins, N_total_bins);
for j = 1:N_selected_vectors
    multiplicitypartialcov = multiplicitypartialcov + multiplicitycoeff(:,j) * multiplicitylatent(j) * multiplicitycoeff(:,j)';
end

multiplicitypartialcorr = corrcov(multiplicitypartialcov);

for j = 1:N_selected_vectors
    multiplicitycoeff(:,j) = 0;
    multiplicitylatent(j) = 0;
end

% print covariance partial
multiplicityfileID_w = fopen('multiplicity_cov_partial.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(multiplicityfileID_w,"%1.15f ",multiplicitypartialcov(i,j));
    end
    fprintf(multiplicityfileID_w,"\n");
end
fclose(multiplicityfileID_w);

% print correlation partial
multiplicityfileID_w = fopen('multiplicity_corr_partial.txt','w');
for i = 1:N_total_bins
    for j = 1:N_total_bins
        fprintf(multiplicityfileID_w,"%1.15f ",multiplicitypartialcorr(i,j));
    end
    fprintf(multiplicityfileID_w,"\n");
end
fclose(multiplicityfileID_w);

% get covariance for remaining
multiplicityremaincov = multiplicitycoeff*diag(multiplicitylatent)*multiplicitycoeff';

% print truncated remaining covariance matrix
multiplicityfileID_w = fopen('multiplicity_cov_remain_truncated.txt','w');
for i = 1:N_total_bins
    fprintf(multiplicityfileID_w,"%1.15f\n",multiplicityremaincov(i,i));
end
fclose(multiplicityfileID_w);

% latent: eigen value
% coeff: eigen vector