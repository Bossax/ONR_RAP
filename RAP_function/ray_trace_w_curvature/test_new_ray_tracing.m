% new ray tracing code test
clear
close all
%% Hydrophone Position on the Ellipsoid 
N_mean = 2.32;
HEM_lat = 22.738772;
HEM_lon = -158.006186;
HEM_depth = -4736.226+N_mean;  

icListen_lat=22.73911569;                  % March 2019 #3
icListen_lon=-158.006106601;                % March 2019
icListen_depth = -4734.646+N_mean;
%% load sample tx points

[tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,x_err,y_err,z_err,act_arrival,est_arrival,SNR] = tx_rx_extraction_Oct(27:30,3,14,'HEM');

%% Geodetic to Local Horizon ENU (receiver)
for ii = 1:length(tx_lat)
    [e(ii),n(ii),u(ii)] = geodetic2enu(HEM_lat,HEM_lon,HEM_depth,tx_lat(ii),tx_lon(ii)+360,tx_altitude(ii),referenceEllipsoid('WGS84'));
end

%% compute estimated travel times of the new code and the old code
surface_dis_layer = zeros(length(tx_lat),2415);
for ii =1:length(tx_lat)
    ii
    geo_S(ii)=dist([HEM_lat tx_lat(ii)],[HEM_lon tx_lon(ii)]);   % Geodesic Length (S)
    azmth(ii) = azimuth(tx_lat(ii),tx_lon(ii),HEM_lat,HEM_lon);
    % new ray tracing
    [~,tot_dist(ii),theta0(ii),SS1,z1,~,surf_handle,est_tt1(ii),~ ] = ray_trace_w_earth_flattening(geo_S(ii),tx_altitude(ii),tx_lat(ii),azmth(ii),HEM_lat,HEM_lon,HEM_depth);
    surface_dis_layer(ii,1:length(surf_handle)) = surf_handle;
    
    [~,tot_dist2(ii),theta02(ii),~,z2,~,surf_handle2,est_tt2(ii),~ ] = ray_trace_planar_coordinate(geo_S(ii),u(ii),tx_altitude(ii),tx_lat(ii),azmth(ii),HEM_lat,HEM_lon);
    
    % old ray tracing
    [~,~,~,SS3,z3,~,~,est_tt3(ii),~]=ray_trace_w_curvature_v3(geo_S(ii),tx_altitude(ii)-2.31);    
    
    hyd_vert_dist1(ii) = z1(end)-z1(1);
    hyd_vert_dist2(ii) = z2(end)-z2(1);
    hyd_vert_dist3(ii) = z3(end)-z3(1);
end
%% azimuth of the ship rel to the ACO
azmth = [];
for ii = 1:length(tx_lat)
   azmth(ii) = azimuth(HEM_lat,HEM_lon,tx_lat(ii),tx_lon(ii)); 
end
range = sum(surface_dis_layer,2);
%% Plot
figure(1)
clf
scatter(range(1:4710)/1000,(est_tt1-est_tt3)*1000,10,azmth(1:4710),'filled')
grid on
xlabel('Range (km)')
ylabel('(msec)')
title('Estimated Time Difference (Planar Coordinate - Old)')
set(gca,'fontsize',12)
colormap jet

figure(2)
scatter(geo_S/1000,range-geo_S','.')
grid on
xlabel('Range (km)')
ylabel('m')

set(gca,'fontsize',12)
title('Surface Range Differecne (New - Old)')

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