% conversion of earth-centered earth-fixed coordiate to local horizon coordinate
clear 
close all

% Geodetic Position of the ACO (Ellipsoid)
ACO_lat = 22.738772;
ACO_lon = -158.006186;
ACO_depth = -4736.226;  
ACO = [ACO_lat ACO_lon ACO_depth];
% ship coordinates
[tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,~,~,~,~,~,~,~] = tx_rx_extraction_Oct(27:30,3,14,'icListen');

% shift the vertical reference to the ellipsoid
tx_altitude = tx_altitude+2.31;

ship_coordinate = [tx_lat' tx_lon'];

%{
% flat earth build in function (Aerospace Tool Box)
N_mean = 2.32;  % mean geoid height
for ii =1:length(tx_lat)
    flatearth_pos(ii,:) = lla2flat(ACO, ship_coordinate(ii,:),90, 0,'WGS84');
end

local_range = sqrt(flatearth_pos(:,1).^2+flatearth_pos(:,2).^2);
%}

% ellipsoidal arc
for ii = 1:length(tx_lat)
    elps_arc(ii) = dist([ACO_lat tx_lat(ii)],[ACO_lon tx_lon(ii)],'wgs84');
end


%% Chadwell 2010
% geodetic to ENU
[e,n,u,X_sc,Y_sc,Z_sc,X_rv,Y_rv,Z_rv,N] = geodetic2enu_chadwell(tx_lat,tx_lon,tx_altitude,ACO_lat,ACO_lon,ACO_depth+N_mean);

local_range2 = sqrt(e.^2+n.^2);
length_diff2 = elps_arc-local_range2;
 
%% Matlab Mapping tool
%  geodetic to ECEF
[X_sc2,Y_sc2,Z_sc2] = geodetic2ecef(tx_lat,tx_lon,tx_altitude,referenceEllipsoid('WGS84'));
[X_rv2,Y_rv2,Z_rv2] = geodetic2ecef(ACO_lat,ACO_lon,ACO_depth+N_mean,referenceEllipsoid('WGS84'));

% geodetic to ENU
xEast = [];
yNorth =[];
zUp = [];
for ii = 1:length(tx_lat)
    [xEast(ii),yNorth(ii),zUp(ii)] = geodetic2enu(ACO_lat,ACO_lon,ACO_depth+N_mean,tx_lat(ii),tx_lon(ii)+360,tx_altitude(ii),referenceEllipsoid('WGS84'));
end
local_range3 = sqrt(xEast.^2+yNorth.^2);
length_diff3 = elps_arc-local_range3;
%% Calculate Azimuth and Heading relative to the ACO
azmth = ones(length(tx_lat),1);
% 0-360
for ii=1:length(tx_lat)
    azmth(ii) = azimuth(ACO_lat,ACO_lon,tx_lat(ii),tx_lon(ii));
end
%% Plot
figure(1)
clf
scatter(tx_lon,tx_lat,20,length_diff3,'filled')
grid on
cbar = colorbar;
cbar.Label.String = '\Deltam';
colormap jet
axis tight
title('Horizntal Length Difference (Ellipsoidal Arc - Local Horizon Length)')

figure(2)
scatter(local_range3/1000,length_diff3,20,azmth,'filled')
ylabel('Distance Difference (m)')
xlabel('Range (km)')
colormap jet
cbar = colorbar;
cbar.Label.String = 'Azimuth (degrees)';
set(gca,'fontsize',13)
grid on
% ylim([-10 10])
title('Horizntal Length Difference (Ellipsoidal Arc - Local Horizon Length)')

%% 
function [e,n,u,X_sc,Y_sc,Z_sc,X_rv,Y_rv,Z_rv,N] = geodetic2enu_chadwell(lat_sc,lon_sc,h_sc,lat_rv,lon_rv,h_rv)
%%%%%%%% Input %%%%%%%%%%%
% Geodetic positions of the source and the receiver (does not support vector)
%%%%%%%%% Output %%%%%%%%
% ECEF and local horizon coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geodetic to ENU Chadwell
% reference ellipsoid
ref_ellipsoid =  referenceEllipsoid('WGS84');
a = ref_ellipsoid.SemimajorAxis;
b = ref_ellipsoid.SemiminorAxis;
ecc = ref_ellipsoid.Eccentricity;
% source ECEF
if lon_sc<0
    lambda_sc =(lon_sc+360)/180*pi;
else
    lambda_sc =(lon_sc)/180*pi;
end

phi_sc = (lat_sc)/180*pi;
N = a./sqrt(1-ecc^2*sin(phi_sc).^2);
X_sc = (N+h_sc).*cos(phi_sc).*cos(lambda_sc);
Y_sc= (N+h_sc).*cos(phi_sc).*sin(lambda_sc);
Z_sc= (N.*(1-ecc^2)+h_sc).*sin(phi_sc);
% receiver ECEF
if lon_rv<0
    lambda_rv =(lon_rv+360)/180*pi;
else
    lambda_rv =(lon_rv)/180*pi;
end
phi_rv= (lat_rv)/180*pi;

Rc_rv = a./sqrt(1-ecc^2*sin(phi_rv).^2);
X_rv = (Rc_rv+h_rv).*cos(phi_rv).*cos(lambda_rv);
Y_rv = (Rc_rv+h_rv).*cos(phi_rv).*sin(lambda_rv);
Z_rv = (Rc_rv.*(1-ecc^2)+h_rv).*sin(phi_rv);

e = [];
n = [];
u = [];
% Local Horizon Coordinate

for iii =1:length(lat_sc)
    tf_mat = [      -sin(lambda_sc(iii))               cos(lambda_sc(iii))                     0;
        -sin(phi_sc(iii))*cos(lambda_sc(iii))       -sin(phi_sc(iii))*sin(lambda_sc(iii)) cos(phi_sc(iii));
        cos(phi_sc(iii))*cos(lambda_sc(iii))        cos(phi_sc(iii))*sin(lambda_sc(iii))  sin(phi_sc(iii))];
    
    lpos=  tf_mat*[X_rv-X_sc(iii);Y_rv-Y_sc(iii);Z_rv-Z_sc(iii)];
    e(iii) = lpos(1);
    n(iii) = lpos(2);
    u(iii) = lpos(3);
end
end