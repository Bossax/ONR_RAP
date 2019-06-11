% Calculate ray arc lengths for corresponding transmission parameters
% Save output file "ray_data contains layerwise surface distance, layerwise ray arc length, and total ray arc length"

clearvars
close all
addpath('/Volumes/LaCie_RAP_Backup/RAP/June2018Cruise/Tx_Rx_Output/script/ray_trace_w_curvature')
addpath('/Volumes/LaCie_RAP_Backup/RAP/June2018Cruise/Tx_Rx_Output/script')
%ACO LAT/LON
ACO_lat=22.738772;                  % Aug 2018
ACO_lon=-158.006186;    

cd ..
cd Path/Radial/Westward    % Edit

%% Load TX data
d_tx=dir('tx_data*.mat');
tx_t = [];
tx_lat = [];
tx_lon = [];
tx_altitude = [];
for p = 1:length(d_tx)
    fname = d_tx(p).name;
    load(fname)
    tx_t = horzcat(tx_t,tx_data.t);      %TX time                                  
    tx_lat = horzcat(tx_lat,tx_data.lat);   %TX lat
    tx_lon = horzcat(tx_lon,tx_data.lon);
    tx_altitude = horzcat(tx_altitude,tx_data.altitude);
end

tx_altitude=tx_altitude-2.31;           %WSG84 Geoid height relative to MSL
cd ../../..
cd script
% Correct t_tx for Scarlett offset
tx_t = tx_t +3/(3600*24);



%% Boat distance from ACO
for i=1:length(tx_lat)
    x_dist(i)=dist([ACO_lat tx_lat(i)],[ACO_lon tx_lon(i)]);
end

%Remove data outside theoretical range (29.5 km)
remove_pts=find(x_dist>29400);
x_dist(remove_pts)=[];
tx_t(remove_pts)=[];
tx_lat(remove_pts)=[];
tx_lon(remove_pts)=[];

%% Main
arc_lengths = [];
tot_dist = [];
surface_dist = [];
for i=1:length(x_dist)
    [arc_handle,tot_handle,~,~,~,~,surface_dist_handle,~,~]=ray_trace_w_curvature(x_dist(i),tx_altitude(i));
    arc_lengths(end+1,:) = arc_handle;
    tot_dist(end+1) = tot_handle;
    surface_dist(end+1,:) = surface_dist_handle;

end

 %% Save

 cd ../variable/Westward
 
 ray_data.ray_arc_dist = tot_dist;
 ray_data.arc_length = arc_lengths;
 ray_data.surface_dist = surface_dist;
 save('ray_data_westward.mat','ray_data')
