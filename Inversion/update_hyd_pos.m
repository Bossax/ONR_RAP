% update hydrophone lat/lon
close all
clear
% 
HEM_lat= 22.738772;                  % June 2017
HEM_lon= -158.006186;                % June 2017
% HEM_depth = 4736.266;              % at 4,736.266 m June 2017

HEM_depth = 4729.92;                 % original depth

% HEM_lat= 22.738783;                    % June 2018
% HEM_lon= -158.00619898;                % June 2018

HEM_lat= 22.7387643;                   % Oct 2018
HEM_lon= -158.00617623;                % Oct 2018
HEM_depth = 4734.58;                   % Oct 2018

icListen_lat=22.739153;                  % June 2018
icListen_lon=-158.0061254;               % June 2018
icListen_depth = 4728.17;                % original depth

% icListen_lat=22.73912734;                  % March 2019#2 
% icListen_lon=-158.006111197;                % March 2019
% icListen_depth = 4735.71;                  %  March 2019

% icListen_lat = 22.7391066;                  % Oct 2018
% icListen_lon = -158.00610724;              % Oct 2018
% icListen_depth = 4733.24;                  %  Oct 2018

% offset in m z down positive
% HEM

% x = 1.01;
% y = -0.84;
% z = 4.66;


% icListen
x = 1.9;
y = -5.12;
z = 5.07;

[X,Y,utmzone] = deg2utm(icListen_lat,icListen_lon);
X_new = X+x ; Y_new = Y+y;
[new_Lat,new_Lon] = utm2deg(X_new,Y_new,utmzone)
 new_depth = icListen_depth+z
 pressure = gsw_p_from_z(-new_depth,icListen_lat)
 % new depth in SS inversion/ ray trace
 % new lat/lon in find_reception/ SS_inversion