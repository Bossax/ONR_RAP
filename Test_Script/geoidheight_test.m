% geoid height function test 
% geoid height statistics of the transmission data
% visualize geoid heights of all acoustic tranmission locations
clear
close all

ACO_lat= 22.738772;                  % June 2017
ACO_lon= -158.006186;                % June 2017
ACO_depth = -4727.6;                 % Original depth

day = 27:30;
start_hour = 3;
end_hour = 14;
hydrophone= 'HEM';
% download tx rx data from October 2018
[tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,x_err,y_err,z_err,act_arrival,est_arrival,SNR] = tx_rx_extraction_Oct(day,start_hour,end_hour,hydrophone);
ttp = (act_arrival-est_arrival)*3600*24*1000;
% geoid height

gh = geoidheight(tx_lat,tx_lon+360,'EGM96');

GH_med = median(gh);
GH_rms = rms(gh-GH_med);

% Calculate Azimute
azmth = ones(length(tx_lat),1);
% 0-360

for i=1:length(tx_lat)
    azmth(i) = azimuth(ACO_lat,ACO_lon,tx_lat(i),tx_lon(i));
end
theta = tx_heading';
theta = tx_heading' - azmth;

for i = 1:length(theta)
   if theta(i) < 0 
       theta(i) = theta(i)+360;
   end
    
end
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
scatter(tx_lon,tx_lat,20,gh-2.3,'filled');
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

%% 
figure(5)
scatter(azmth,gh,[],range,'filled')
grid on
xticks(0:60:360)
xticks(0:30:360)
xlim([0 360])
title('Geoid Height vs Azimuth')
xlabel('Azimuth')
ylabel('Geoid HEight (ms)')
set(gca,'fontsize',15)
c = colorbar;
c.Label.String = 'Range (km)';
c.Ticks = 0:5:30;
caxis([0 30])
colormap jet

