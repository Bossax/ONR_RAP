% testing algorithm to convert geodetic position to local east north depth
close all
clear

hyd_lat = 22.738772; 
hyd_lon =  -158.006186; 
hyd_height = -4730;
% generate geodetic positions
lat = hyd_lat:-0.01:22.558772;    
lon =  hyd_lon:0.01:-157.826186;
height = -4*ones(1,length(lat));
% convert geodetic positions to ECEFG
source = geodetic2ECEFG(lat,lon,height);
hyd = geodetic2ECEFG(hyd_lat,hyd_lon,hyd_height);
% convert to surface ranges and depths in ENU
[enu_range,enu_depth] = XYZ2ENU(source.X,source.Y,source.Z,source.lambda,source.phi,hyd.X,hyd.Y,hyd.Z);
enu_range(1) = [];
enu_depth(1) = [];
height(1) = [];
hyd_depth_enu = enu_depth-height-2.32;
slant_range_enu = sqrt(enu_range.^2+enu_depth.^2);
% Vincenty Geodesic
geodesic = distance(lat,lon,hyd_lat,hyd_lon,referenceEllipsoid('WGS84'));
geodesic(1) = [];
slant_range_geo = sqrt(geodesic.^2+(hyd_height-height).^2);
%% Ray tracing

for ii = 10:length(enu_range)
    [~,S_enu(ii),~,~,z_enu,~,~,~,tt_enu(ii),~,~] = ray_tracing_planar_const_SS(...
                                                enu_range(ii),height(ii),hyd_depth_enu(ii),3);
    [~,S_geo(ii),~,~,z_geo,~,~,~,tt_geodesic(ii),~,~] = ray_tracing_planar_const_SS(...
                                                geodesic(ii),height(ii),-hyd_height,3);
end
%%
figure(21)
clf
subplot(3,1,1)
scatter(enu_range/1000,(tt_geodesic-tt_enu)*1000)
grid on
xlabel('Range (km)')
ylabel('msec')
title('Travel Time Difference (Geodesic - ENU)')
subplot(3,1,2)
scatter(enu_range/1000,(S_geo-S_enu))
grid on
xlabel('Range (km)')
ylabel('meter')
title('Ray Arc Difference (Geodesic - ENU)')

subplot(3,1,3)
scatter(enu_range/1000,(slant_range_geo-slant_range_enu))
grid on
xlabel('Range (km)')
ylabel('meter')
title('Slant Range Difference (Geodesic - ENU)')

figure(22)
subplot(2,1,1)
scatter(enu_range/1000,geodesic-enu_range)
grid on
xlabel('Range (km)')
ylabel('meter')
title('Horizontal Range Difference (Geodesic - ENU)')

subplot(2,1,2)
scatter(enu_range/1000,(-hyd_height)-hyd_depth_enu)
grid on
xlabel('Range (km)')
ylabel('meter')
title('Depth Difference (Geodesic - ENU)')

%%
function ECEFG = geodetic2ECEFG(lat_sc,lon_sc,h_sc)
%%%%%%%% Input %%%%%%%%%%%
% Geodetic position
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

ECEFG.X = X_sc;
ECEFG.Y = Y_sc;
ECEFG.Z = Z_sc;
ECEFG.lambda = lambda_sc;
ECEFG.phi = phi_sc;

end

function [enu_range,enu_depth] = XYZ2ENU(X_sc,Y_sc,Z_sc,lambda_sc,phi_sc,X_rv,Y_rv,Z_rv)

for iii =1:length(X_sc)
    
    tf_mat = [      -sin(lambda_sc(iii))               cos(lambda_sc(iii))                     0;
        -sin(phi_sc(iii))*cos(lambda_sc(iii))       -sin(phi_sc(iii))*sin(lambda_sc(iii)) cos(phi_sc(iii));
        cos(phi_sc(iii))*cos(lambda_sc(iii))        cos(phi_sc(iii))*sin(lambda_sc(iii))  sin(phi_sc(iii))];
    
    lpos=  tf_mat*[X_rv-X_sc(iii);Y_rv-Y_sc(iii);Z_rv-Z_sc(iii)];
    e(iii) = lpos(1);
    n(iii) = lpos(2);
    u(iii) = lpos(3);
end

enu_range = sqrt(e.^2+n.^2);
enu_depth = -u;
end

function [arc_lengths_ly,rayarc_dist_tot,theta0,SS_flt,z_flt,SS_flt_HOT_avg,horz_dist,inc_r_ly,est_tt,tt,ray_angles_ly] = ray_tracing_planar_const_SS(x_dist,td_h,hyd_depth,EFT,frame)
% arc_length= arc length in each layer
% tot_dist  = total surface distance calculated from numerical solution
% theta0    = calculated launch angle
% SS        = Sound speed profile (with depth z)
% z         = Depth (negative/ below the transducer)
% SS_HOT    = depth averaged sound speed
% surface_dist = surface distance in each layer
% est_tt    = estimated travel time
% ray angle = incident ray angle in each layer
% ETF: 1 = general formula, 2 = Chadwell and Swiney
% Reference Ellipsoid
ref_ellipsoid =  referenceEllipsoid('WGS84');
major_ax = ref_ellipsoid.SemimajorAxis;
minor_ax = ref_ellipsoid.SemiminorAxis;
ecc = ref_ellipsoid.Eccentricity;

% remove bad data
if td_h >=0
    td_h = 0;
end

ACO_lat = 22.738772;
ACO_lon = -158.006186;

%% 1. CTD data (MSL)
% CTD file directory
% cd '/Users/testuser/Documents/MATLAB/Script/Data'
fid=fopen('h306a0202.ctd');                       % 12 October 2018

D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f','headerlines',6));
fclose(fid);

pres_MSL = D(:,1);        %pressure (dbars)
temp = D(:,2);        %temperature
sal = D(:,3);         %salinity (Sp)
sal = gsw_SA_from_SP(sal,pres_MSL,ACO_lon,ACO_lat);

% Depths in the ellipsoidal coordinate; downward positive (-h)
ACO_gh = geoidheight(ACO_lat,ACO_lon+360,'EGM96');
z_td = -td_h;
z_hyd = hyd_depth;

% disp(sprintf('transducer depth = %f',z_td))
%% 2 truncate CTD data to fit the range between the receiver depth and the source depth
% on the ellipsoidal surface Orthometrc height

N_mean = 2.32;                                      % average geiod height over the area of coverage
z_CTD = -(gsw_z_from_p(pres_MSL,ACO_lat) );         % relative to MSL
z_elps = z_CTD-N_mean;                              % from MSL to ellipsoid


% 2.1 truncate the upper bound using transducer depth
rm_ind = find(z_elps<=z_td);
z_elps(rm_ind(1:end-1)) = [];
z_CTD(rm_ind(1:end-1)) = [];
temp(rm_ind(1:end-1)) = [];
sal(rm_ind(1:end-1)) = [];

% set the first depth to be the transducer depth
z_diff = z_td - z_elps(1);
temp(1) = temp(1)+(temp(2)-temp(1))*(z_diff/(z_elps(2)-z_elps(1)));   % interpolate the first temp and salinity
sal(1)=sal(1)+(sal(2)-sal(1))*(z_diff/(z_elps(2)-z_elps(1)));
z_elps(1) = z_td;
z_CTD(1) = z_elps(1)+N_mean;    % transducer depth + mean geoid undulation = transducer depth relative to MSL


% 2.2 truncate the lower bound
if z_elps(end) > z_hyd
    rm_ind = find(z_elps>z_hyd);
    z_elps(rm_ind) = [];
    z_CTD(rm_ind) = [];
    temp(rm_ind) = [];
    sal(rm_ind) = [];
else
    % extrapolate the lower bound
    ly_thk = -mean(z_elps(end-5:end-1) - z_elps(end-4:end));
    add_layer = floor((z_hyd -z_elps(end))/ly_thk);
    z_elps = vertcat(z_elps,z_elps(end)+ly_thk*(1:add_layer)');
    z_elps(end+1) = z_hyd;
    
    z_CTD = vertcat(z_CTD,z_elps(end-add_layer:end)+N_mean);
    
    temp_gradient = temp(end-1)-temp(end);
    temp = vertcat(temp,temp(end)+temp_gradient*(1:add_layer+1)');
    sal_gradient = sal(end-1)-sal(end);
    sal = vertcat(sal,sal(end)+sal_gradient*(1:add_layer+1)');
end

pres_CTD = gsw_p_from_z(-z_CTD,ACO_lat);
%% 3. Creat Sound Speed Profile
SS_elps = gsw_sound_speed(sal,temp,pres_CTD);


%% Earth Flattening Transformation
R_e = 6371000;
if EFT == 1
    % Dushaw and Colosi
    E = z_elps./R_e;
    z_flt = z_elps.*(1+(E./2)+((E.^2)./3));
    SS_flt = SS_elps.*(1+E+E.^2);
    
elseif EFT == 2
    % Chadwell 2010
    z_flt = -R_e.*log(R_e./(R_e-(-z_elps)));
    SS_flt = SS_elps.*(R_e./(R_e-(-z_elps)));
elseif EFT == 3
    % no curvature correction
    z_flt = z_elps;
    SS_flt = SS_elps;
else
    z_flt = -R_e.*log((R_e-(z_elps))./R_e);
    SS_flt = SS_elps.*(R_e./(R_e-(z_elps)));
end

r = [];

%Find theta0 of acoustic transmission
% incident angle
% slant anle
slant_angle = atan(x_dist/z_hyd);
theta0_all= slant_angle-10/180*pi:1*(pi/180):slant_angle+10/180*pi;   % all possible launch angles
theta0_all(find(theta0_all <= 0.0001)) = [];
theta0_all(find(theta0_all >= 89.5/180*pi)) = [];


% loop over all possible thetas
% calculate a small horizontal distance within an interval betweenf two depths
for ii = 1:length(theta0_all)
    
    % Ray parameter or each launch angle
    a=sin(theta0_all(ii))/SS_flt(1);
    
    % incident angles of each layer
    theta(1) = theta0_all(ii);
    for i = 2:length(SS_flt)-1
        theta(i)=asin(a*SS_flt(i));
    end
    
    % horizontal distance
    
    r = diff(z_flt).*a.*SS_flt(1:end-1)./sqrt(1-a^2*SS_flt(1:end-1).^2);
    
    
    %Total x distance
    x_diff(ii)=sum(r);  % total horizontal distance for one SS profile launching theta
    
    
end


%Find min and max launch angle for Netwon-Raph method

t_low = find(x_diff<x_dist);      % find all launch angles which yield x distance within the bound
theta_low = theta0_all(t_low(end));   % label the lowest poSS_fltible launch angle (correspond to the steepest angle)
x_low = x_diff(t_low(end));           % The shortest ray arc length, lower bound of the sur dist

t_high = find(x_diff>x_dist);
theta_high = theta0_all(t_high(1));   % The upper bound of sur dist (the least steep angle)
x_high = x_diff(t_high(1));

% Intrapolation the launch angle to find the optimum ray arc length of the direct path
theta_new = theta_low+((x_dist)/((x_low-x_high)/(theta_low-theta_high)))-((x_low)/((x_low-x_high)/(theta_low-theta_high)));

% We now have a range of launch angles which gives the nearest number of horizontal distance
theta1 = theta_low;
x1 = x_low;           % lower bound of the horizontal distance

x_diff = 1000000;   % arbitary number

% Loop until x_diff is close to x_dist to 1 mm

while abs(x_diff-x_dist) > 0.0001
    arc_distance = [];
    tt = [];
    r = [];
    % Ray parameter
    a=sin(theta_new)/SS_flt(1);
    
    % incident angle
    theta=asin(a*SS_flt(2:end-1));
    theta = vertcat(theta_new,theta);
    
    
    % horizontal range in each layer
    for i=1:length(SS_flt)-1
        r(i)=(z_flt(i+1)-z_flt(i))*tan(theta(i));
    end
    
    %Total x distance
    x_diff = sum(r);
    
    %Total distance and time traveled
    for i=1:length(SS_flt)-1
        arc_distance(i)=(z_flt(i+1)-z_flt(i))/cos(theta(i));
        %         tt(i)=(z_flt(i+1)-z_flt(i))*(1/SS_flt(i)^2)/sqrt((1/SS_flt(i)^2)-a^2);
        tt(i) = arc_distance(i)/SS_flt(i);
    end
    tot_dist_all=sum(arc_distance);
    
    arc_lengths_hold = arc_distance;
    surface_dist_hold = r;
    ray_angles_hold=theta(2:end);
    
    
    if abs(x_diff-x_dist)<0.0001
        break;
    end
    
    % Newton-Raphson method
    % Interpolate again to find a new launch angle which is closer to that
    % of the eigenray
    
    theta2 = theta_new;
    
    
    x2=x_diff;
    
    theta_new=theta2+((x_dist)/((x2-x1)/(theta2-theta1)))-((x2)/((x2-x1)/(theta2-theta1)));
    
    % Set theta_low to theta 2 (new angle)
    theta1=theta2;
    % update nearest horizontal distance
    x1=x2;
end


% packaging outputs
theta0 = theta_new;                        % Launch angle
rayarc_dist_tot = tot_dist_all;                   % Ray arc length
arc_lengths_ly = arc_lengths_hold;            % Arc length in each layer
horz_dist = x_diff;          % horizontal distance
ray_angles_ly = ray_angles_hold;                % Incident ray angle at each interval

% incremental horizontal distance
inc_r_ly = zeros(1,length(surface_dist_hold)+1);

for i= 2:length(inc_r_ly)
    inc_r_ly(i) = inc_r_ly(i-1)+surface_dist_hold(i-1);
end



est_tt=sum(tt);       %Estimated TT

% Vertically averaged SS_flt for HOTS CTD
SS_flt_HOT_avg=mean(SS_flt);

end