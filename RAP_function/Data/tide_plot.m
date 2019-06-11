% Tidal changes
% after cruise
clear
cd /Users/testuser/Documents/MATLAB/Script/Data/
load('rtx_posmvGGA_aftercruise.mat')

cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/posmv/after
d_posmv = dir('posmv*');
bin_t = [];
bin_altitude  = [];

for iii =1:length(d_posmv)
    fname = d_posmv(iii).name;
    load (fname)
    bin_t = vertcat(bin_t,posmv.t);
    bin_altitude = vertcat(bin_altitude,posmv.altitude);
end
[bin_t,I] = unique(bin_t);
bin_altitude = bin_altitude(I);

tide = csvread('ioc_stations_31Oct2018.csv',1,4);
t_tide = datenum('20181031 00:00:00','yyyymmdd HH:MM:SS') + [0:60/(3600*24):(length(tide)-1)*60/(3600*24)];


figure(1)
clf
plot(rtx_t,rtx_altitude-median(rtx_altitude),'g');
hold on
plot(posmvgga_t,posmvgga_altitude-median(posmvgga_altitude),'r');
plot(bin_t,bin_altitude-median(bin_altitude),'b')
plot(t_tide,tide-median(tide),'k')
datetick('x')
grid on
ylabel(' Median Offset (m)')
xlabel('Time')
title('Tidal signal of the GPS systems')
legend('RTX','POSMV GGA','POSMV Binary','Observed Tide')