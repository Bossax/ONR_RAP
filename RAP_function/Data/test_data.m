close all
clear
load('RTX_beforecruise_1Hz.mat')
t1Hz = RTX.rtx_t;
lat1Hz= RTX.rtx_lat;
lon1Hz= RTX.rtx_lon;
altitude1Hz= RTX.rtx_altitude;

load('RTX_beforecruise_10Hz.mat')
t10Hz = RTX.rtx_t;
lat10Hz= RTX.rtx_lat;
lon10Hz= RTX.rtx_lon;
altitude10Hz= RTX.rtx_altitude;

load('posmv_20181024.mat')
posmvs2_t = posmv.t;
posmvs2_lat = posmv.lat;
posmvs2_lon = posmv.lon;
posmvs2_altitude = posmv.altitude;

posmvs1_t = posmv.t;
posmvs1_lat = posmv.lat_s;
posmvs1_lon = posmv.lon_s;
posmvs1_altitude = posmv.altitude_s;

figure(1)
clf
subplot(2,1,1)
plot((posmvs2_t-t1Hz(1))*3600*24,posmvs2_altitude)
hold on
plot((t1Hz-t1Hz(1))*3600*24,altitude1Hz)
plot((t10Hz-t1Hz(1))*(3600*24),altitude10Hz)
ylim([37.2 38])
grid on
% xlim([20000 20500])

subplot(2,1,2)
plot((posmvs2_t-t1Hz(1))*3600*24,posmvs2_altitude)
hold on
plot((t1Hz-t1Hz(1))*3600*24,altitude1Hz)
plot((t10Hz-t1Hz(1))*(3600*24)-35,altitude10Hz)
ylim([37.2 38])
grid on
% xlim([20000 20500])