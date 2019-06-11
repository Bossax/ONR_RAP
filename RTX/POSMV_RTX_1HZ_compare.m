% POS MV and RTX 1Hz comparison
% 1. Download POSMV 1Hz data and RTX
% 2. plot 2 datasets on top of each other
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/posmv/1Hz_4_RTX
d_posmv = dir('pos*');
posmv_ta = [];
posmv_lata = [];
posmv_lona = [];
posmv_altitudea = [];
posmv_rolla = [];
posmv_pitcha = [];
posmv_headinga = [];
for k = 1:length(d_posmv)
    posmv_fname = d_posmv(k).name
    load(posmv_fname)
    posmv_t = posmv.t;
    posmv_lat = posmv.lat;
    posmv_lon = posmv.lon;
    posmv_altitude = posmv.altitude';
    posmv_roll = posmv.roll;
    posmv_pitch = posmv.pitch;
    posmv_heading = posmv.heading;
    
    posmv_ta = horzcat(posmv_ta,posmv_t');
    posmv_lata = horzcat(posmv_lata,posmv_lat');
    posmv_lona = horzcat(posmv_lona,posmv_lon');
    posmv_altitudea = horzcat(posmv_altitudea,posmv_altitude);
    posmv_rolla = horzcat(posmv_rolla,posmv_roll);
   posmv_pitcha =  horzcat(posmv_pitcha,posmv_pitch);
   posmv_headinga = horzcat(posmv_headinga,posmv_heading);
end
%% RTX
% Load
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/RTX/
load kilo_moana_1819_RTX.dat
A = kilo_moana_1819_RTX;

% Timestamps
% Modified Julian Date
mjd = A(:,1);
jd = mjd +2400000.5;
rtx_ta = [];
% Convert to Gregorian Date
for i = 1:length(jd)
   rtx_ta(end+1) = datenum(datetime(jd(i),'ConvertFrom','juliandate'));
end

% other parameters
rtx_lat_ant = A(:,2)';
rtx_lon_ant = A(:,3)';
rtx_altitude_ant = A(:,4)';
geoid_height = A(:,8)';

% Transformation
[rtx_lata,rtx_lona,rtx_altitudea] = transducer_pos_rtx(rtx_lat_ant,rtx_lon_ant,rtx_altitude_ant,rtx_ta,posmv_rolla,posmv_pitcha,posmv_headinga,posmv_ta);


%% Plot
f = figure(1);
clf
f.Units = 'normalized';
f.Position = [0.2 0.2 0.6 0.7];
subplot(3,1,1)
scatter(posmv_ta,posmv_lona,12,'k','filled')
hold on
scatter(rtx_ta,rtx_lona,'x')
datetick('x')
legend('POSMV','RTX')
grid on
title('Absolute Transducer Positions of POSMV and RTX')
axis tight
ylabel('Lon')

subplot(3,1,2)
scatter(posmv_ta,posmv_lata,12,'k','filled')
hold on
scatter(rtx_ta,rtx_lata,'x')
datetick('x')
legend('POSMV','RTX')
grid on
ylabel('Lat')
axis tight

subplot(3,1,3)
scatter(posmv_ta,posmv_altitudea,12,'k','filled')
hold on
scatter(rtx_ta,rtx_altitudea,'x')
datetick('x')
legend('POSMV','RTX')
grid on
ylabel('Z (m)')
axis tight






