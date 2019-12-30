clear
close all
%Load TX data
cd /Users/testuser/Documents/ONR_RAP/Data/Tx_Rx_Output/October2018/tx_file
fname = dir('*.mat');
tx_t = [];
tx_lat = [];
tx_lon = [];
tx_altitude = [];
tx_heading = [];
tx_xvel = [];
x_err = [];
y_err = [];
z_err = [];
for p = 1:length(fname)
    load(fname(p).name)
    tx_t = horzcat(tx_t,tx_data.t);      %TX time
    tx_lat = horzcat(tx_lat,tx_data.lat);   %TX lat
    tx_lon = horzcat(tx_lon,tx_data.lon);
    tx_altitude = horzcat(tx_altitude,tx_data.altitude);
    tx_heading = horzcat(tx_heading,tx_data.heading);
    tx_xvel = horzcat(tx_xvel,tx_data.x_vel);
    x_err = horzcat(x_err,tx_data.lon_err);
    y_err = horzcat(y_err,tx_data.lat_err);
    z_err = horzcat(z_err,tx_data.altitude_err);
end

HEM_lat=22.738772;                  % Aug 2018
HEM_lon=-158.006186;                % Aug 2018

lat2 = HEM_lat*ones(size(tx_lat));
lon2 = HEM_lon*ones(size(tx_lon));

[x1,~] = distance(tx_lat,tx_lon,lat2,lon2,referenceEllipsoid('WGS84'));
x2 = [];
for ii = 1:length(x1)
    x2(ii) = vdist(tx_lat(ii),tx_lon(ii),HEM_lat,HEM_lon);
    
end


