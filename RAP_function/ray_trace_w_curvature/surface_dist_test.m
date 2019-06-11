% comparoison between with and without cv correction
% x_dist = random walk until reah 25 km
% z_offset = random between -7 to -5
% z positive downward
clear 
close all

ACO_lat=22.738772;                  % June 2018
ACO_lon=-158.006186;                % June2018  

icListen_lat=22.73911569;                  % March 2019 #3
icListen_lon=-158.006106601;                % March 2019

slant_range1 = [];
day =27:30;
start_hour = 3;
end_hour = 14;
[tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,~,x_err,y_err,z_err,act_arrival,est_arrival,SNR] = tx_rx_extraction_Oct(day,start_hour,end_hour,'HEM');

%% dist.m ellipsoid
for ii=1:length(tx_lat)
    dist_dist(ii)=dist([ACO_lat tx_lat(ii)],[ACO_lon tx_lon(ii)],'wgs84');
end

% deg2utmUTM
[x,y,utmzone] = deg2utm(tx_lat,tx_lon);
[ACO_x,ACO_y,utmzone] = deg2utm(ACO_lat,ACO_lon);
utm_dist = sqrt((x-ACO_x).^2+(y-ACO_y).^2);
delta_s_deg = dist_dist' -utm_dist;

% pos2dist (plane estimate on spherical surface)
for iii = 1:length(tx_lat)
    pos_dist(iii) = pos2dist(ACO_lat,ACO_lon,tx_lat(iii),tx_lon(iii),1);
    gsw_dist(iii) = gsw_distance([tx_lon(iii) ACO_lon],[tx_lat(iii) ACO_lat],0);
end
del_s_pos = dist_dist - pos_dist*1000;

% gsw_distance (spherical)
del_s_gsw =  dist_dist - gsw_dist;

% Matlab distance
[dist_matlab,~] = distance(ACO_lat,ACO_lon,tx_lat,tx_lon,wgs84Ellipsoid);
del_s_matlab = dist_dist - dist_matlab;

% lldiskm
for iii =1:length(tx_lat)
    [d1km(iii) d2km(iii)]=lldistkm([ACO_lat ACO_lon],[tx_lat(iii) tx_lon(iii)]);
end
del_s_lldiskm = dist_dist-d1km*1000;

%% Calculate Azimuth and Heading relative to the ACO
azmth = ones(length(tx_lat),1);
% 0-360
for ii=1:length(tx_lat)
    azmth(ii) = azimuth(icListen_lat,icListen_lon,tx_lat(ii),tx_lon(ii));
end

%%
figure(1)
scatter(tx_lon,tx_lat,20,del_s_gsw,'filled')
grid on
% caxis([2 8])
cbar = colorbar;
title('Surface Distance Difference (dist.m - deg2utm.m)')
cbar.Label.String = '\Delta m';
colormap jet
axis tight
set(gca,'fontsize',13)
%%
figure(22)
scatter(dist_dist,delta_s_deg,20,azmth,'filled')
grid on
ylabel('Distance Difference (m)')
xlabel('Range (km)')
colormap jet
cbar = colorbar;
cbar.Label.String = 'Azimuth (degrees)';
set(gca,'fontsize',13)
title('Surface Distance Difference (dist.m - deg2utm.m) vs Range')
