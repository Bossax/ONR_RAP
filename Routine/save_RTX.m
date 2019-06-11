% save RTX data
clear
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/RTX/before
load rtx10Hz.dat
A = rtx10Hz;
%%
% Timestamps
% Modified Julian Date
mjd = A(:,1);
jd = mjd +2400000.5;
rtx_t = [];
% Convert to Gregorian Date
rtx_t = datenum(datetime(jd,'ConvertFrom','juliandate'));

%%
RTX.rtx_t = rtx_t';
% other parameters
RTX.rtx_lat = A(:,2);
RTX.rtx_lon= A(:,3);
RTX.rtx_altitude= A(:,4);
RTX.rtx_lat_err = A(:,5);
RTX.rtx_lon_err = A(:,6);
RTX.rtx_altitude_err = A(:,7);
RTX.geoid_height = A(:,8);
%% save
save('RTX_beforecruise.mat','RTX')