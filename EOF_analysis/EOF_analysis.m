function [EOFs,lambda,contribution,PCs]=EOF_analysis(anomaly)
% the programe contaims eof decompositon and significance test
% 1)the data must be anomaly without 'NaN', not actual value
% 2)mXn,expressing m grids and n times
% 3)[anomalies should be area-weighted£¬sqrt(cos(lat)),not cos(lat)]
% Syntax:
% 1)[EOFs,lambda,contribution,PCs]=EOF_analysis(anomaly);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2015 Chunl¨¹e Zhou, from SYSU and BNU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0.0 size of anomaly
[m,n]=size(anomaly);
% 1.0 calculate covariance matrix
C=anomaly*anomaly'/n;

% 2.0 Do eigenanalysis
[EOFs,D]=eig(C);
PCs=EOFs'*anomaly;

% 3.0 reverse the order
EOFs=fliplr(EOFs);
PCs=flipud(PCs);
D=fliplr(flipud(D));
lambda=diag(D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.0 North et al., (1982) significant test at 95% confidence
% 1) effective freedoms:
for j=1:m
    %alpha(j,:)=autocorr(anomaly(j,:),1);
    alpha(:,j)=acf(anomaly(j,:)',1);
end
ralph=(mean(abs(alpha(:,2)))+max(alpha(:,2)))/2.;
tau=-1.0/log(ralph);
Nstar=n/(2*tau);
factor=sqrt(2./Nstar);
% 3£©Plot Eigenvalues
L=lambda;
L_error=L*factor;
contribution=[100*L/sum(L),cumsum(100*L/sum(L))];
figure
errorbar(1:20,L(1:20),L_error(1:20));
title(['Eigenvalue: 1-20']);
xlabel('x mode');