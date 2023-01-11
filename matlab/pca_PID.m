clear all
format long

% Kaon ID
KIDfileID = fopen('KID_toys.txt','r');
formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ';
sizeK = [105 Inf];
KaonMatrix = fscanf(KIDfileID,formatSpec,sizeK);
KaonMatrix = KaonMatrix';
CovKaonID = cov(KaonMatrix);

[Kaoncoeff,Kaonlatent] = pcacov(CovKaonID);
Kaoncoeff
Kaonlatent


% Pion ID
PIDfileID = fopen('PID_toys.txt','r');
formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ';
sizeP = [105 Inf];
PionMatrix = fscanf(KIDfileID,formatSpec,sizeP);
PionMatrix = PionMatrix';
CovPionID = cov(PionMatrix);

[Pioncoeff,Pionlatent] = pcacov(CovPionID);
Pioncoeff
Pionlatent

% latent: eigen value
% coeff: eigen vector