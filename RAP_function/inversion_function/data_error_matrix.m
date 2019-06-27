function [Cd] = data_error_matrix(x_err,y_err,z_err,tx_lat,tx_lon,tx_heading,theta_l,theta_r,SNR,cl,cr,hyd_lat,hyd_lon,hyd_depth)
% generate measurement uncertainty matrix          
% diagonal elements are sum of time measurement uncertainty and the source location uncertianty 
% transmission time uncertainty is mcuh greater than reception 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. RMS time measurement error
wrms = 795.661;     % Hz
sig_t = 1./(wrms.*sqrt(SNR));
var_t = sig_t.^2;

% 2.source location uncertainty
%  Phi
azmth = azimuth(ones(1,length(tx_lat))*hyd_lat,ones(1,length(tx_lat))*hyd_lon,tx_lat,tx_lon);
phi = (360-azmth)-90;
phi(find(phi<0)) = phi(find(phi<0))+360;
% phi = phi/180*pi;

cl = 1537.77;   % october
cr = 1539.026;      % october

x_anomaly=(((sin(theta_l)./cl).*cosd(phi)).^2).*(x_err).^2;
y_anomaly=(((sin(theta_l)./cl).*cosd(phi)).^2).*(y_err).^2;
z_anomaly=(cos(theta_l)./cl).*(z_err).^2;

% transducer location variance
% var_p = ((sin(theta_r).*cosd(phi)./cr).^2).*(x_err.^2)+((sin(theta_r).*sind(phi)./cr).^2).*(y_err.^2)+((sin(theta_r)./cr).^2).*(z_err.^2);
var_p2 = x_anomaly+y_anomaly+z_anomaly;


% 3. transmission time RMS uncertainty
rms_tx_t = 0.05*1e-3;       % 0.05 ms
var_tx_t = rms_tx_t^2;

% var_p(find(var_p > 1e9*median(var_p))) = median(var_p);
% 4. Form the martix
Cd = var_t+var_tx_t+var_p2;
ind_outlier = find(Cd >= 3*median(Cd));
Cd(ind_outlier) = median(Cd);



end