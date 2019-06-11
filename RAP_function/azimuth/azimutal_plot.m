% Travel time perturbation plotted agaisnt Azimutal angle of the circular path
%
clearvars
close all

%% Load Tx data 
cd /Volumes/LaCie_RAP_Backup/RAP/June2018Cruise/Tx_Rx_Output/script/azimuth
cd ../..
cd tx_file_path/CCW    % Edit

d_tx=dir('tx_data*.mat');
tx_t = [];
tx_lat = [];
tx_lon = [];
tx_altitude = [];
tx_heading = [];
for p = 1:length(d_tx)
    fname = d_tx(p).name;
    load(fname)
    tx_t = horzcat(tx_t,tx_data.t);      %TX time                                  
    tx_lat = horzcat(tx_lat,tx_data.lat);   %TX lat
    tx_lon = horzcat(tx_lon,tx_data.lon);
    tx_altitude = horzcat(tx_altitude,tx_data.altitude);
    %tx_heading = horzcat(tx_heading,tx_data.heading);
end

tx_altitude=tx_altitude-2.31;           % WSG84 Geoid height relative to MSL

% Correct t_tx for Scarlett offset
tx_t = tx_t +3/(3600*24);       % 3 second June 2018

%ACO LAT/LON
ACO_lat=22.738772;                  % Aug 2018
ACO_lon=-158.006186; 
cntr = [ACO_lon ACO_lat];
%Boat distance from ACO
for i=1:length(tx_lat)
    x_dist(i)=dist([ACO_lat tx_lat(i)],[ACO_lon tx_lon(i)]);
end

%Remove data outside theoretical range (29.5 km)
remove_pts=find(x_dist>29400);
x_dist(remove_pts)=[];
tx_t(remove_pts)=[];
tx_lat(remove_pts)=[];
tx_lon(remove_pts)=[];

%% Load Rx files
cd ../..
cd rx_files_info
load('rx_data_CCW.mat')
est_arrival = rx_data.est_arrival;
act_arrival = rx_data.act_arrival;
tt_perb = act_arrival-est_arrival;

%% Take out tx data points
rm_ind = find((tx_t - act_arrival(1))*3600*24 < -20);
tx_lat(rm_ind) = [];
tx_lon(rm_ind) = [];
tx_t(rm_ind) = [];
rm_ind = find((tx_t - act_arrival(end))*3600*24 >20);
tx_lat(rm_ind) = [];
tx_lon(rm_ind) = [];
tx_t(rm_ind) = [];
%% Calculate azimuthal angles
azmth = [];
addpath('/Volumes/LaCie_RAP_Backup/RAP/June2018Cruise/Tx_Rx_Output/script/azimuth');
for i=1:length(tx_lat)
    azmth(end+1) = azimuth(ACO_lat,ACO_lon,tx_lat(i),tx_lon(i));
end
ins = find(199<azmth & azmth<202);

%% Plot
f = figure(1);
f.Units = 'normalized';
f.Position = [0.1 0.7 0.6 0.4];
c = 1:1:length(tt_perb);    % colormap scale
scatter(azmth,tt_perb*3600*24*1000,[],c,'fill')
colormap jet
cbar = colorbar;
cbar.Label.String = 'Transmission Number';
cbar.Ticks = 0:50:350;
grid on
xlabel('Azimuthal Angle')
ylabel('Travel Time Perturbation (ms)')
xlim([0 360])
title('Counterclockwise Circle')


f = figure(2);
clf
f.Units = 'normalized';
f.Position = [0.1 0.7 0.4 0.5];
c2 = tt_perb*3600*24*1000;
scatter(tx_lon,tx_lat,[],c2);

colormap jet
cbar = colorbar;
cbar.Label.String = 'Travel Time Perturbation (ms)';
caxis([4 9])

axis tight
grid on
title('Plane View CounterClockwise')
hold on
scatter(ACO_lon,ACO_lat,200,'kx')



