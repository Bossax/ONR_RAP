% comparoison between with and without cv correction
% x_dist = random walk until reah 25 km
% z_offset = random between -5 to `-3 meters (mean = -4 m)
% z positive downward
clear
close all

z_offset0 = -4;
td_offset = -4;
% x_dist = [500];
x_dist = [10 500 1000 7000 10000 15000 20000 25000 26000];
% x_dist = [5000 10000 20000];
% while x_dist(end) <= 26000
%     x_dist = cat(1,x_dist,x_dist(end)+rand*5000+200);
%     td_offset = cat(1,td_offset,z_offset0+rand*-2);
% end


%% Ray tracing looping
% figure(99)
% figure(100)
% set(gcf,'Unit','normalized','Position',[0 0.5 0.4 0.4])
for iii = 1:length(x_dist)
    
    % APL
    [arc1,tot_dist1(iii),theta01(iii),SS1,z1,~,surface_dist1(iii),inc_r1,est_tt1(iii),tt1,ray_angles1] = ray_tracing_planar_const_SS(x_dist(iii),td_offset(1),1);
    theta01(iii) = pi-theta01(iii);
%     % Chadwell
%     [arc2,tot_dist2(iii),theta02(iii),SS2,z2,~,surface_dist2(iii),inc_r2,est_tt2(iii),tt2,ray_angles2] = ray_tracing_planar_const_SS(x_dist(iii),td_offset(1),2);
%     theta02(iii) = pi-theta02(iii);
%     % No correction
    [arc3,tot_dist3(iii),theta03(iii),SS3,z3,~,surface_dist3(iii),inc_r3,est_tt3(iii),tt3,ray_angles3] = ray_tracing_planar_const_SS(x_dist(iii),td_offset(1),3);
    theta03(iii) = pi-theta03(iii);
%     % spherical coordinate
%     [arc4,tot_dist4(iii),theta04(iii),SS4,z4,~,surface_dist4(iii),inc_r4,est_tt4(iii),tt4,ray_angles4] = ray_tracing_sperhical(x_dist(iii),td_offset(1));
%     theta04(iii) = pi-theta04(iii);
%      With ETF and ss gradient
    [arc5,tot_dist5(iii),theta05(iii),SS5,z5,~,surface_dist5(iii),inc_r5,est_tt5(iii),tt5,ray_angles5] = ray_tracing_planar_vary_SS(x_dist(iii),td_offset(1),3);
	
    % APl integration equation
%     [est_tt7(iii),~,theta07(iii)] = ray_tracing_planar_const_SS_APL(x_dist(iii),td_offset(1));
%     theta07(iii) = pi-theta07(iii);
    % spherical coordinate with incremental horizontal range
    %   [arc6,tot_dist6(iii),theta06(iii),SS6,z6,~,surface_dist6(iii),inc_r6,est_tt6(iii),tt6,ray_angles6] = ray_tracing_sperhical_dr(x_dist(iii),td_offset(iii));
   

end
%% plot
figure(1)
clf
ax1 = axes();
scatter(x_dist(1:iii)/1000,(est_tt3-est_tt1)*1000,'*');
xlabel('Range (km)')
ylabel('msec')
yyaxis right 
scatter(x_dist(1:iii)/1000,(theta03-theta01)*180/pi,'*')
% title('The planar ray tracing with EFT and the spherical coordinate ray tracing (m)')
ylabel('degree')

ax1.FontSize = 14;
grid on

%%
%{
figure(2)
scatter(x_dist(1:iii)/1000,(tot_dist5-tot_dist4),'*')
title('Arc Length Difference (inaccurate) (m)')
grid on
%}
%% Function
function [arc_lengths_ly,rayarc_dist_tot,theta0,SS_flt,z_flt,SS_flt_HOT_avg,horz_dist,inc_r_ly,est_tt,tt,ray_angles_ly] = ray_tracing_planar_const_SS(x_dist,td_h,EFT)
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
ACO_depth = -4729.92;                   % original depth MSL local u/e/n frame
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
z_hyd = -(ACO_depth+ACO_gh);

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

function [arc_lengths_ly,rayarc_dist_tot,theta0,SS_flt,z_flt,SS_flt_HOT_avg,horz_dist,inc_r_ly,est_tt,dt,ray_angles_ly] = ray_tracing_planar_vary_SS(x_dist,td_h,EFT)
% linearly changing sound speed
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
ACO_depth = -4729.92;                   % original depth MSL local u/e/n frame
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
z_hyd = -(ACO_depth+ACO_gh);

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
% grazing angle
% slant anle
slant_angle = acot(x_dist/z_hyd);
upper_ang = slant_angle+15/180*pi;
if upper_ang >=89.9/pi*180
    upper_ang = 89.9/pi*180;
end

lower_ang = slant_angle-15/180*pi;
if lower_ang<= 1/180*pi
    lower_ang = 1/180*pi;
end
theta0_all= lower_ang:0.5*(pi/180):upper_ang;   % all possible launch angles


% loop over all possible thetas
% calculate a small horizontal distance within an interval betweenf two depths
% sound speed gradient
g = diff(SS_flt)./diff(z_flt);

for ii = 1:length(theta0_all)
    
    % Ray parameter of each launch angle
    a=cos(theta0_all(ii))/SS_flt(1);
    
    %  grazing angles of each layer
    theta(1) = theta0_all(ii);
    for i = 2:length(SS_flt)-1
        theta(i)=acos(a*SS_flt(i));
    end
    
    % horizontal distance
    
    r = (sqrt(1-a^2*SS_flt(1:end-1).^2) - sqrt(1-a^2*SS_flt(2:end).^2))./(a*g);
    
    
    %Total x distance
    x_diff(ii)=sum(r);  % total horizontal distance for one SS profile launching theta
    
    
end
rm_ind = find(imag(x_diff) ~= 0 );
x_diff(rm_ind) = [];
theta0_all(rm_ind) = [];

%Find min and max launch angle for Netwon-Raph method

% farther range bound
t_high = find(x_diff >= x_dist ,1 ,'last');      
theta_low = theta0_all(t_high);   
x_high = x_diff(t_high);           

t_low = find(x_diff <= x_dist,1,'first');
theta_high = theta0_all(t_low);   % The upper bound of sur dist (the least steep angle)
x_low = x_diff(t_low);

% Intrapolation the launch angle to find the optimum ray arc length of the direct path
theta_new = theta_low+((x_dist)/((x_low-x_high)/(theta_low-theta_high)))-((x_low)/((x_low-x_high)/(theta_low-theta_high)));

% We now have a range of launch angles which gives the nearest number of horizontal distance
theta1 = theta_low;
x1 = x_low;           % lower bound of the horizontal distance

geodesic = 100000;   % arbitary number

% Loop until x_diff is close to x_dist to 1 mm

while abs(geodesic-x_dist) > 0.0001
    
    % Ray parameter
    a=cos(theta_new)/SS_flt(1);
    
    % incident angle
    theta=asin(a*SS_flt(2:end-1));
    theta = vertcat(theta_new,theta);
    
    dr = (sqrt(1-a^2*SS_flt(1:end-1).^2) - sqrt(1-a^2*SS_flt(2:end).^2))./(a*g);
    dt = log((SS_flt(2:end)./SS_flt(1:end-1)).*(1+sqrt(1-a^2*SS_flt(1:end-1).^2))./(1+sqrt(1-a^2*SS_flt(2:end).^2)))./g;
    ds = (acos(a*SS_flt(1:end-1)) - acos(a*SS_flt(2:end)))./(a*g);
    
    %Total x distance
    geodesic = sum(dr);
    
    if abs(geodesic-x_dist)<0.0001
        break;
    end
    
    % Newton-Raphson method
    % Interpolate again to find a new launch angle which is closer to that
    % of the eigenray
    
    theta2 = theta_new;
    
    
    x2 = geodesic;
    
    theta_new = theta2+((x_dist)/((x2-x1)/(theta2-theta1)))-((x2)/((x2-x1)/(theta2-theta1)));
    
    % Set theta_low to theta 2 (new angle)
    theta1=theta2;
    % update nearest horizontal distance
    x1=x2;
end



% packaging outputs
theta0 = theta_new;                     % Launch angle
rayarc_dist_tot = sum(ds);               % Ray arc length
arc_lengths_ly = ds;                   % Arc length in each layer
horz_dist = geodesic;                  % horizontal distance
ray_angles_ly = theta;                 % Incident ray angle at each interval
est_tt=sum(dt);                        % Estimated T
    
% cumulative horizontal distance
inc_r_ly = cumsum(dr);





% Vertically averaged SS_flt for HOTS CTD
SS_flt_HOT_avg=mean(SS_flt);

end

function [est_tt,depth,theta0] = ray_tracing_planar_const_SS_APL(x_dist,td_h)
% employs APL ray tracing equaitons and algorithm with EFT
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
ACO_depth = -4729.92;                   % original depth MSL local u/e/n frame
%% 1. CTD data (MSL)
% CTD file directory
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
z_hyd = -(ACO_depth+ACO_gh);

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

% Dushaw and Colosi
E = z_elps./R_e;
z_flt = z_elps.*(1+(E./2)+((E.^2)./3));
SS_flt = SS_elps.*(1+E+E.^2);

% Launch rays
% incident launch angles
theta0_all = [0.01:0.1:89]/180*pi;

r = []; % horizontal range

% loop over all possible theta0
% calculate a small horizontal distance within an interval betweenf two depths

nray = length(theta0_all);
nlayer = length(z_elps)-1;
% generate dr vector
interval = 500;
dr = round(x_dist/interval);

if dr >= 1 % 1 m
    dr = 1;
elseif dr == 0
    dr = 0.1;
end

ddr = dr:dr:x_dist;
nr = length(ddr)
% SS gradient 
ssg = diff(SS_flt)./diff(z_flt);

for ii = 1:nray
    
    z_n = z_elps(1);
    theta_n = theta0_all(ii);
    a = sin(theta_n)/SS_flt(1);
    dz = [];
    % loop over horizontal inervals
    for k = 1:nr
        n = find(z_elps <= z_n,1,'last');
        if n >= length(ssg)
            n = length(ssg);
        end
        c_n = SS_flt(n)+ssg(n)*(z_n- z_elps(n));
        theta_n = asin(a*c_n);
        dz(k) = cot(theta_n)*dr;
        z_n = z_n +  dz(k);
        
    end
    
    depth(ii) = sum(dz);
end
    

    %Find min and max launch angle for Netwon-Raph method
    
    t_low = find(depth > z_hyd,1,'last');     
    theta_low = theta0_all(t_low);   
    d_low = depth(t_low);           
    
    t_high = find(depth < z_hyd,1,'first');
    theta_high = theta0_all(t_high);   % The upper bound of sur dist (the least steep angle)
    d_high = depth(t_high);
    
    % Intrapolation the launch angle to find the optimum ray arc length of the direct path
    theta_new = theta_low+((z_hyd)/((d_low-d_high)/(theta_low-theta_high)))-((d_low)/((d_low-d_high)/(theta_low-theta_high)));
    
    % We now have a range of launch angles which gives the nearest number of horizontal distance
    theta1 = theta_low;
    d1 = d_low;           % lower bound of the horizontal distance
    
    depth = 1000000;   % arbitary number
    
    % Loop until  depth is within 0.1 mm
    
    while abs(depth - z_hyd) > 0.0001
        z_n = z_elps(1);
        a = sin(theta_new)/SS_flt(1);
        dz = [];
        theta_n = [] ;
        c_n = [];
        % loop over horizontal inervals
        for k = 1:nr
            n = find(z_n >= z_elps,1,'last');
            if n >= length(ssg)
            n = length(ssg);
        end
            
            c_n(k) = SS_flt(n)+ssg(n)*(z_n- z_elps(n));
            theta_n(k) = asin(a*c_n(k));
            dz(k) = cot(theta_n(k))*dr;
            z_n = z_n +  dz(k);
            
        end
       
        depth = sum(dz);
        
        % Newton-Raphson method
        % Interpolate again to find a new launch angle which is closer to that
        % of the eigenray
        
        theta2 = theta_new;
        d2=depth;
        
        theta_new=theta2+((z_hyd)/((d2-d1)/(theta2-theta1)))-((d2)/((d2-d1)/(theta2-theta1)));
        
        % Set theta_low to theta 2 (new angle)
        theta1 = theta2;
        % update nearest horizontal distance
        d1 = d2;
        
    end
    
    % calculate travel time and horizontal distance
    dt = dr./c_n./sin(theta_n);
    
    %% packaging outputs
    theta0 = theta_new;                        % Launch angle 
    est_tt=sum(dt);       %Estimated TT
    
    
    figure(99)
    clf
    scatter(SS_flt,z_flt,'*')
    hold on
    scatter(c_n,cumsum(dz),'*')
    set(gca,'YDir','reverse')
    grid on
    pause
    
end

function [arc_lengths_ly,rayarc_dist_tot,theta0,SS_elps,z_elps,SS_flt_HOT_avg,horz_dist,inc_r_ly,est_tt,tt,ray_angles_ly ] = ray_tracing_sperhical(x_dist,td_h)
        
        % remove bad data
        if td_h >=0
            td_h = 0;
        end
        
        ACO_lat = 22.738772;
        ACO_lon = -158.006186;
        ACO_depth = -4729.92;                   % original depth MSL local u/e/n frame
        R_e = 6371000;          % Earth Radius
        %% 1. CTD data (MSL)
        fid=fopen('h306a0202.ctd');                       % 12 October 2018
        
        D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f','headerlines',6));
        fclose(fid);
        
        pres_MSL = D(:,1);    % pressure (dbars)
        temp = D(:,2);        % temperature
        sal = D(:,3);         % salinity (Sp)
        sal = gsw_SA_from_SP(sal,pres_MSL,ACO_lon,ACO_lat);
        
        % Depths in the ellipsoidal coordinate; downward positive (-h)
        ACO_gh = geoidheight(ACO_lat,ACO_lon+360,'EGM96');
        z_td = -td_h;
        z_hyd = -(ACO_depth+ACO_gh);
        
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
        
        % 3. Creat Sound Speed Profile
        SS_elps = gsw_sound_speed(sal,temp,pres_CTD);
        
        
        
        
        %Find theta0 of acoustic transmission
        theta0_all = [0.01:0.1:89.5]/180*pi;
        
        
        % loop over all possible thetas
        % calculate a small horizontal distance within an interval betweenf two depths
        ang_dis = [];             % layer-wise angular displacement
        R0 = R_e - z_td;           % Launch radial distance
        R_hyd = R_e - z_hyd;      % Hydrophone radial distance
        Ri = R_e - z_elps;
        for ii=1:length(theta0_all)
            
            % 1st incident angle at interface of layer 2
            incident_angle = asin(Ri(1)/Ri(2)*sin(theta0_all(ii)));
            
            %Ray parameter
            p = R0*sin(incident_angle)/SS_elps(1);
            
            
            % leaving angles of each layer
            theta = asin(p*SS_elps(2:end)./Ri(2:end));
            theta = vertcat(incident_angle,theta);
            
            % angular displacement in each layer
            % circumferential displace ment in each layer
            for gg = 1:length(SS_elps)-1
                
                %         Rmean = (Ri(gg)+Ri(gg+1))/2;
                %         ang_dis(gg)= p*SS_elps(gg)*(Ri(gg) - Ri(gg+1))/(Rmean^2*sqrt(1-(p*SS_elps(gg)/Rmean)^2));
                %         c_dis(gg) = ang_dis(gg)*Rmean;
                
                % APL with curavture correction factor
                fe = (R_e - z_elps(gg))/R_e;
                c_dis(gg) = (z_elps(gg+1)-z_elps(gg))*tan(theta(gg))/fe;
                
            end
            
            % Total geodetic distance
            geo_dis(ii)=sum(c_dis);  % total horizontal distance
            
        end
        
        
        %Find min and max launch angle for Netwon-Raph method
        
        t_low=find(geo_dis<x_dist);      % find all launch angles which yield x distance within the bound
        theta_low=theta0_all(t_low(end));   % label the lowest poSS_elpsible launch angle (correspond to the steepest angle)
        x_low=geo_dis(t_low(end));           % The shortest ray arc length, lower bound of the sur dist
        
        t_high=find(geo_dis>x_dist);
        theta_high=theta0_all(t_high(1));   % The upper bound of sur dist (the least steep angle)
        x_high=geo_dis(t_high(1));
        
        % Intrapolation to find the optimum ray arc length of the direct path
        theta_new = theta_low+((x_dist)/((x_low-x_high)/(theta_low-theta_high)))-((x_low)/((x_low-x_high)/(theta_low-theta_high)));
        
        % We now have a range of launch angles which gives the nearest number of
        % horizontal distance
        theta1=theta_low;
        x1=x_low;           % lower bound of the horizontal distance
        
        geo_dis=1000000;
        c_dis = [];
        ang_dis = [];
        % Loop until x_diff is close to x_dist to 0.1 mm
        while abs(geo_dis-x_dist)>0.0001
            
            % 1st incident angle at interface of layer 2
            incident_angle = asin(Ri(1)/Ri(2)*sin(theta_new));
            
            %Ray parameter
            p = R0*sin(incident_angle)/SS_elps(1);
            
            % incident angles of each layer
            theta = asin(p*SS_elps(2:end)./Ri(2:end));
            theta = vertcat(incident_angle,theta);
            
            % Horizontal displacement
            for gg = 1:length(Ri)-1
                % APL
                fe = (R_e - z_elps(gg))/R_e;
                c_dis(gg) = (z_elps(gg+1)-z_elps(gg))*tan(theta(gg))/fe;
                
            end
            % Total geodetic distance
            geo_dis = sum(c_dis);
            
            if abs(geo_dis-x_dist)<0.0001
                break;
            end
            
            % Newton-Raphson method
            % Interpolate again to find a new launch angle which is closer to that of the eigenray
            
            theta2 = theta_new;
            
            
            x2 = geo_dis;
            
            theta_new = theta2+((x_dist)/((x2-x1)/(theta2-theta1)))-((x2)/((x2-x1)/(theta2-theta1)));
            
            % Set theta_low to theta 2 (new angle)
            theta1 = theta2;
            
            % update nearest horizontal distance
            x1 = x2;
            
        end
        
        incident_angle = asin(Ri(1)/Ri(2)*sin(theta_new));
        % incident angles of each layer
        theta = asin(p*SS_elps(2:end)./Ri(2:end));
        theta = vertcat(incident_angle,theta);
        % calculate travel time
        fe = (R_e - z_elps(1:end-1))/R_e;
        tt = c_dis'.*fe.*csc(theta(1:end-1))./SS_elps(1:end-1);
        arc_distance = tt.*SS_elps(1:end-1);  % approximate from travel time and sound speed
        
        
        surface_dist_hold = c_dis;
        ray_angles_hold = theta(2:end);
        
        %% packaging outputs
        theta0 = theta_new;                        % Launch angle
        rayarc_dist_tot = sum(arc_distance);                   % Ray arc length
        arc_lengths_ly = arc_distance;            % Arc length in each layer
        horz_dist = geo_dis;          % horizontal distance
        ray_angles_ly = ray_angles_hold';                % Incident ray angle at each interval
        
        % incremental horizontal distance
        inc_r_ly = zeros(1,length(surface_dist_hold)+1);
        
        for i= 2:length(inc_r_ly)
            inc_r_ly(i) = inc_r_ly(i-1)+surface_dist_hold(i-1);
        end
        
        est_tt=sum(tt);       %Estimated TT
        
        
        SS_flt_HOT_avg = mean(SS_elps);
        
    end


