% Test sound speed profile error
clear 
close all

x_dist = [10];

z_offset0 = -4;
td_offset = -4;


while x_dist(end) <= 26000
    x_dist = cat(1,x_dist,x_dist(end)+rand*300+200);
    td_offset = cat(1,td_offset,z_offset0+rand*-2);
end


%% Ray tracing looping
for iii = 1:length(x_dist)
    iii
    [arc1,tot_dist1(iii),theta01(iii),SS1,z1,~,surface_dist1(iii),inc_r1,est_tt1(iii),tt1,ray_angles1] = ray_tracing_sperhical(x_dist(iii),td_offset(iii),0);
    [arc2,tot_dist2(iii),theta02(iii),SS2,z2,~,surface_dist2(iii),inc_r2,est_tt2(iii),tt2,ray_angles2] = ray_tracing_sperhical(x_dist(iii),td_offset(iii),1);
    
end
%% Plot
scatter(x_dist/1000,(est_tt2-est_tt1)*1000,'.')
grid on
xlabel('km')
ylabel('msec')
%%


function [arc_lengths_ly,rayarc_dist_tot,theta0,SS_elps,z_elps,SS_flt_HOT_avg,horz_dist,inc_r_ly,est_tt,tt,ray_angles_ly ] = ray_tracing_sperhical(x_dist,td_h,SS_error)

% remove bad data
if td_h >=0
    td_h = 0;
end

ACO_lat = 22.738772;
ACO_lon = -158.006186;
ACO_depth = -4729.92;                   % original depth MSL local u/e/n frame
R_e = 6371000;          % Earth Radius
%% 1. CTD data (MSL)
% if SS_error == 1
%     fname = 'h306a0202.ctd';
% else
%     fname = 'h306a0215.ctd';
% end

fname = 'h306a0202.ctd';
fid=fopen(fname);                       % 12 October 2018

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


if SS_error == 1
%     ss_error = -std(SS_elps)*0.05*rand(length(SS_elps),1);
%     ss_error = -std(SS_elps)*0.05;
% specific inetrval
    start_l = 100;
    end_l = 2200;
    ss_error = zeros(length(SS_elps),1);
    ss_error(start_l:end_l) = -(SS_elps(start_l:end_l))*0.001;

else
    ss_error = 0;
end

SS_elps = SS_elps + ss_error;
%}

%Find theta0 of acoustic transmission
% incident angle
% slant angle
slant_angle = atan(x_dist/z_hyd);
theta0_all= slant_angle-10/180*pi:.1*(pi/180):slant_angle+10/180*pi;   
theta0_all(find(theta0_all <= 0.0001)) = [];
theta0_all(find(theta0_all >= 89.5/180*pi)) = [];


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
