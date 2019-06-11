%  geoidheight function test 
% geoid height statisti of the transmission data
clear
close all
day = 27:30;
start_hour = 3;
end_hour = 14;
hydrophone= 'icListen';
[tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,x_err,y_err,z_err,act_arrival,est_arrival,SNR] = tx_rx_extraction_Oct(day,start_hour,end_hour,hydrophone);
ttp = (act_arrival-est_arrival)*3600*24*1000;
% geoid height

gh = geoidheight(tx_lat,tx_lon+360,'EGM96');

GH_med = median(gh);
GH_rms = rms(gh-GH_med);
%%
figure(1)
clf
scatter(tx_t,gh,20,ttp,'filled')

grid on
datetick('x')

colormap jet
cbar = colorbar;
cbar.Label.String = 'TTP(ms)';
headline = sprintf('Median = %f m; RMS = %f m',GH_med,GH_rms);
title(headline)
axis tight
%%
figure(2)
clf
scatter(range,gh,20,ttp,'filled')
grid on
colormap jet
cbar = colorbar;
cbar.Label.String = 'TTP(ms)';
ylabel('Geiod Height (m)')
xlabel('Range (km)')
% title('Geiod Height vs Range')
set(gca,'fontsize',13)
%% 
figure(3)
scatter(tx_lon,tx_lat,20,gh,'filled');
grid on
colormap jet
cbar = colorbar;
cbar.Label.String = 'Geiod Height (m)';
axis tight
set(gca,'fontsize',13)
title('Geoid Height Map')
%%
figure(4)
tx_altitude2 = tx_altitude+2.31;
tx_altitude2 = tx_altitude2-gh;

scatter(range,tx_altitude - tx_altitude2,20,ttp,'filled')
grid on

colormap jet
cbar = colorbar;
cbar.Label.String = 'TTP(ms)';
set(gca,'fontsize',13)
axis tight


