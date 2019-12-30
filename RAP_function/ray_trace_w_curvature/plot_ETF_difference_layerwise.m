% new ray tracing code test (Chadwell 2010)
clear
close all
%% Hydrophone Position
N_mean = 2.32;
HEM_lat = 22.738772;
HEM_lon = -158.006186;
HEM_depth = -4729.92;                   % original depth MSL local u/e/n frame
% HEM_depth = -2000;
%% load sample tx points
[tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,~,x_err,y_err,z_err,~,~,~] = tx_rx_extraction_Oct(27:30,3,14,'HEM');


%% compute estimated travel times of the new code and the old code
est_tt1 = [];
est_tt2 = [];
est_tt3 = [];
geo_S = [10 1000 5000 10000 15000 26000];
% surface_dis_layer = zeros(length(tx_lat),2415);
arclength1 = [];
arclength2 = [];
tx_altitude = -4.2 ; % ellipsoidal height
tx_lat = 22.7;
tx_lon = -158.1;
azmth = 30;

figure(1001)
clf
plot_xoffset = 0.005;
margin = 0.06;
plotwidth = 0.1;
set(gcf,'Units','normalized','Position',[0 .7 0.9 0.6])
for ii =1:length(geo_S)
   
    
%     geo_S(ii) = distance(tx_lat(ii),tx_lon(ii),HEM_lat,HEM_lon,referenceEllipsoid('WGS84')); % geodesic
%     azmth(ii) = azimuth(tx_lat(ii),tx_lon(ii),HEM_lat,HEM_lon);
    disp(geo_S(ii))
    %ray tracing with EFT
    [arc1,arclength1(ii),~,SS1,z1,~,~,~,est_tt1(ii),tt1,~,~] = ray_trace_w_earth_flattening(geo_S(ii),tx_altitude,tx_lon,tx_lat,azmth,HEM_lat,HEM_lon,HEM_depth,'Oct','2018');
    
    %ray tracing without EFT
    [arc2,arclength2(ii),~,SS2,z2,~,~,~,est_tt2(ii),tt2,~,~] = ray_trace_w_earth_flattening(geo_S(ii),tx_altitude,tx_lon,tx_lat,azmth,HEM_lat,HEM_lon,HEM_depth,'Oct','2018',false);

%{
     surface_dis_layer(ii,1:length(surf_handle)) = surf_handle;
    % new ray tracing Chadwell and Sweeney 2010 with respect to geoid
%     [~,~,~,SS2,z2,~,~,~,est_tt2(ii),~,~] = ray_trace_w_earth_flattening_noVshift(geo_S(ii),tx_altitude,tx_lon(1),tx_lat(1),azmth(1),HEM_lat,HEM_lon,HEM_depth,'Oct','2018');
    
    % old ray tracing
%     [~,~,~,SS3,z3,~,~,est_tt3(ii),~]=ray_trace_w_curvature_v3(geo_S(ii),tx_altitude(ii)-2.31,HEM_lat,HEM_lon,HEM_depth,'Oct','2018');
    %}
    
    hyd_vert_dist1(ii) = z1(end)-z1(1);
    hyd_vert_dist2(ii) = z2(end)-z2(1);
%     hyd_vert_dist3(ii) = z3(end)-z3(1);

    ax = axes;
    set(ax,'Units','normalized','Position',[plot_xoffset+(ii-1)*(plotwidth+margin)+0.05 0.1 plotwidth 0.75]);
    plot(z2(2:end),cumsum((tt2-tt1)*1000),'k')
    ylabel('msec')
    
    grid on
  
    axis tight
    xlabel('Depth (m)')
    title(sprintf('%i m, %.2f msec, %.2f m',geo_S(ii),sum(tt2-tt1)*1000,sum(arc2-arc1)))


    
    yyaxis right
    plot(z2(2:end),cumsum((arc2(1:end)-arc1(1:end))))
    ylabel('meter')
%      ylim([-10 10])
%     set(gca,'View',[90 90])
end
t = annotation('textbox',[.4 .88 .1 .1],'String','Cumulative Sum: Travel Time and Arc Length','Fontsize',14,'BackgroundColor','white');
%% azimuth of the ship rel to the ACO
azmth = [];
for ii = 1:length(tx_lat)
   azmth(ii) = azimuth(HEM_lat,HEM_lon,tx_lat(ii),tx_lon(ii)); 
end
% range = sum(surface_dis_layer,2);
%% Plot
figure(1)
clf
scatter(geo_S/1000,(est_tt2-est_tt1)*1000,20,'filled')
grid on
xlabel('Range (km)')
ylabel('msec')
title({'Estimated Travel Time Difference', 'No EFT -  with EFT'})
% title({'Estimated Travel Time Difference','(vertical reference at MSL - ellipsoidal surface)'})
set(gca,'fontsize',14)
% yticks(-2:2:16)
% ylim([-0.007 0])

%% 
figure(2)
scatter(geo_S/1000,range-geo_S','.')
grid on
xlabel('Range (km)')
ylabel('m')

set(gca,'fontsize',12)
title('Surface Range Differecne (New - Old)')

%% SSP and depth comparison
figure(3)
clf
set(gcf,'Units','normalized','Position',[0.2 0.6 .6 0.5]);
subplot(1,2,1)
plot(z1(1:length(z2)-1)-z2(1:end-1))
grid on
ylabel('Depth difference (m)')
xlabel('Layer')
xlim([0 length(z2)-1])
% ylim([-6 0])
% yticks(horzcat(horzcat(-6:1:-3,-2.32),-2:1:0))
% text(40,z1(1) - z3(1)+0.1,'Mean Geoid Height')
set(gca,'fontsize',14)

subplot(1,2,2)
plot(SS1(1:length(SS2))-SS2)
grid on
ylabel('Sound speed difference (m/s)')
xlabel('Layer')
xlim([0 length(z2)-1])
set(gca,'fontsize',14)

t = annotation('textbox',[.35 .9 .1 .1],'String','Flattened Sound Speed Profile Difference','Fontsize',15,'EdgeColor','none','FontWeight','bold');

%%
figure
arc_diff = arclength2-arclength1;
plot(arc_diff)
yyaxis right
plot(geo_S/1000)
grid on


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