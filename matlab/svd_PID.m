clear all
format long

% Kaon ID
KIDfileID = fopen('KID_toys.txt','r');
formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ';
sizeK = [105 Inf];
KaonMatrix = fscanf(KIDfileID,formatSpec,sizeK);
KaonMatrix = KaonMatrix';
CovKaonID = cov(KaonMatrix);
fclose(KIDfileID);

[U,Kaonlatent,Kaoncoeff] = svd(CovKaonID);
Kaonlatent
Kaon_eigen_vector1 = Kaoncoeff(:,1);
Kaon_eigen_vector2 = Kaoncoeff(:,2);
Kaon_eigen_vector3 = Kaoncoeff(:,3);
Kaon_eigen_value1 = sqrt(Kaonlatent(1,1));
Kaon_eigen_value2 = sqrt(Kaonlatent(2,2));
Kaon_eigen_value3 = sqrt(Kaonlatent(3,3));
Kaoncoeff(:,1) = 0;
Kaoncoeff(:,2) = 0;
Kaoncoeff(:,3) = 0;

Kaonpartialcov = Kaon_eigen_vector1*Kaon_eigen_value1*Kaon_eigen_value1*Kaon_eigen_vector1'+...
    Kaon_eigen_vector2*Kaon_eigen_value2*Kaon_eigen_value2*Kaon_eigen_vector2'+...
    Kaon_eigen_vector3*Kaon_eigen_value3*Kaon_eigen_value3*Kaon_eigen_vector3';

CorrKaonID = corrcov(CovKaonID);
Kaonpartialcorr = corrcov(Kaonpartialcov);

% print KID info into file
KIDfileID_w = fopen('KID_selected.txt','w');
fprintf(KIDfileID_w,"%1.15f\n",Kaon_eigen_value1);
for i = 1:105
        fprintf(KIDfileID_w,"%1.15f ",Kaon_eigen_vector1(i));
end
fprintf(KIDfileID_w,"\n");

fprintf(KIDfileID_w,"%1.15f\n",Kaon_eigen_value2);
for i = 1:105
        fprintf(KIDfileID_w,"%1.15f ",Kaon_eigen_vector2(i));
end
fprintf(KIDfileID_w,"\n");

fprintf(KIDfileID_w,"%1.15f\n",Kaon_eigen_value3);
for i = 1:105
        fprintf(KIDfileID_w,"%1.15f ",Kaon_eigen_vector3(i));
end
fprintf(KIDfileID_w,"\n");

fclose(KIDfileID_w);

KIDfileID_w = fopen('KID_cov_all.txt','w');
for i = 1:105
    for j = 1:105
        fprintf(KIDfileID_w,"%1.15f ",CovKaonID(i,j));
    end
    fprintf(KIDfileID_w,"\n");
end
fclose(KIDfileID_w);

KIDfileID_w = fopen('KID_cov_partial.txt','w');
for i = 1:105
    for j = 1:105
        fprintf(KIDfileID_w,"%1.15f ",Kaonpartialcov(i,j));
    end
    fprintf(KIDfileID_w,"\n");
end
fclose(KIDfileID_w);

KIDfileID_w = fopen('KID_corr_all.txt','w');
for i = 1:105
    for j = 1:105
        fprintf(KIDfileID_w,"%1.15f ",CorrKaonID(i,j));
    end
    fprintf(KIDfileID_w,"\n");
end
fclose(KIDfileID_w);

KIDfileID_w = fopen('KID_corr_partial.txt','w');
for i = 1:105
    for j = 1:105
        fprintf(KIDfileID_w,"%1.15f ",Kaonpartialcorr(i,j));
    end
    fprintf(KIDfileID_w,"\n");
end
fclose(KIDfileID_w);


% latent: eigen value
% coeff: eigen vector