% RTX data exploration
% Plot samples versus time
% interval between 2 contiguous samples
% Load
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/RTX/
load kilo_moana_1819_RTX.dat
A = kilo_moana_1819_RTX;

% Timestamps
% Modified Julian Date
mjd = A(:,1);
jd = mjd +2400000.5;
rtx_t = [];
% Convert to Gregorian Date
for i = 1:length(jd)
   rtx_t(end+1) = datenum(datetime(jd(i),'ConvertFrom','juliandate'));
end

% other parameters
rtx_lat_ant = A(:,2)';
rtx_lon_ant = A(:,3)';
rtx_altitude_ant = A(:,4)';
geoid_height = A(:,8)';
%%
int_sample =[];
for k = 2:length(rtx_t)
   int_sample(end+1) = (rtx_t(k)-rtx_t(k-1))*3600*24;   % sec
end

%% Plot
figure(1)
n = 1:length(rtx_t);
scatter(n,rtx_t);
xlabel('Sample')
ylabel('Timestamp')
grid on
datetick('y')
title('RTX Sample')

figure(2)
n = 1:length(rtx_t);
plot(rtx_t(2:end),int_sample);
xlabel('Time')
ylabel('Interval (sec)')
grid on
datetick('x')
headline = sprintf('RTX Sample intervals: Median = %.3f s',median(int_sample))
title(headline)

%%
figure(3)
plot(rtx_t,rtx_altitude_ant-28.47)
datetick('x')
grid on
axis tight