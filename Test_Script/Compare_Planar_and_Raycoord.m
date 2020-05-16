% Use ray equations and integration technique from "Computationl Ocean 
% Acoustics" textbook (page 201).
%
% Equaitons are expressed in Ray coordinates (along s).
%
% Use both depth and horizontal range to determine eigenrays.
%
% Arc length (ds) is the differential variable
%
% Compare the result of this method with the planar ray tarcing used in RAP (Dushaw & Colosi)

clear
close all
 
% surface distances used for computation
x_dist = [500 1000 5000 10000 15000 20000 25000];
% hydrophone depth
z_hyd = 4730;
% transducer depth
z_td = 5;

% load sound speed profile (including horizontal gradient dc/dr(z)
SS = soundspeedprofile(z_hyd);

%% Calculating loop over all distances
for ii = 1:length(x_dist)
    disp(ii)
    % Planar ray tracing with Earth Flattening Transformation
    [est_time2(ii),theta2(ii),arclength2(ii)]= ray_tracing_planar_const_SS(x_dist(ii),z_td,SS);
    % Ray tacing in ray coordinate system 
    [est_time1(ii),theta1(ii),arclength1(ii)] = ray_trace_raycoor(x_dist(ii),z_td,SS);

end
theta2 = theta2/pi*180;
theta1 = 90-theta1;
%% Plot the results
subplot(1,2,1)
line(x_dist/1000,((est_time2-est_time1)*1000),'Linewidth',2,'Marker','o')
grid on
xlabel('Range (km)')
ylabel('\delta t')
ylabel('\delta t (ms)')
title({'Travel time difference (TT_{Snell} - TT_{Ekinal})'})
set(gca,'fontsize',14)
xlim([0 25])
ylim([-5 5])

subplot(1,2,2)
% line(SS.dcdr,SS.z,'Linewidth',2)
line(zeros(1,length(SS.z)),SS.z,'Linewidth',2);
grid on
set(gca,'YDir','reverse')
xlabel('dc/dr')
ylabel('Depth (m)')
title('Horizontal Sound Speed Gradient')
set(gca,'fontsize',14)
axis tight
%% ray tracing in ray coordinate
function [est_time,theta0,arclength] = ray_trace_raycoor(geodesic,td_z,SS)
% input 1. geodesic
%       2. hyd depth
%       3. transducer depth
%       4. sound speed profile
rcv_depth = SS.z(end);

slant_angle = acot(geodesic/rcv_depth)/pi*180;
upper_ang = slant_angle+10;
if upper_ang >=89.9
    upper_ang = 89.9;
end



lower_ang = slant_angle-15;
if lower_ang<= 1
    lower_ang = 1;
end
theta0= upper_ang:-0.5:lower_ang;   % all possible launch angles

sinc = .5;     % arc length increment 1m 

r0 = 0; % source r0
z0 = td_z; % source depth z0


% store depth and range of each launch angle
Range = [];
Depth = [];
break_trig = 0;
for n = 1:length(theta0)
    
% initial state parameters  
  idt = 1; idz = 2; idr = 3; idsh = 4; idsv = 5;
    [c0,~] = find_soundspeed(z0,SS);
    X0(idt) = 0;                                % time      
    X0(idz)= z0;                             % depth
    X0(idr)= r0;                            % range
    X0(idsh)= cos(pi*theta0(n)/180)./c0;    % horizontal slowness
    X0(idsv)= sin(pi*theta0(n)/180)./c0;    % vertical slowness
   
    % initial value
    snow = 0; 
    Xnow = X0;
    znow = z0;
    rnow = r0;
    X = Xnow;
    % loop until the horizonal range exceed rend
    while (rnow < geodesic)
        % calculate the changes of parameters within this arc length interval
        f = diff_X(snow,Xnow,SS);
    
        % update the state vector 
        Xnow = Xnow + sinc*f;
        snow = snow+sinc;
        znow = Xnow(end,idz);
        rnow = Xnow(end,idr);
        
        % concatenate the state vector
        X = vertcat(X,Xnow);
        
        % check if the ray hits the bottom
        if znow >= rcv_depth 
            break;
        end        
    end
    Range(n) = rnow;
    Depth(n) = znow;
    
    % break the loop
    if break_trig == 1
             
             fprintf('Stop scanning at launhc angle =  %f\n',theta0(n-1))
             break_trig = 0;
             break;
    end
    
     % check if the depth starts to decrease
     if n >= 2
         if ((Depth(n-1) - znow) >= 5)
             % proceed for one more angle
             break_trig = 1;
         end
         
     end
    
end
%% find angle
ind = n-1;
theta1= theta0(ind);      % longer range/ smaller angle
theta2 = theta0(ind-1);
r1= Range(ind);            % longer range
r2 = Range(ind-1);         
depth1 = Depth(ind);  % longer range/ shallower
depth2 = Depth(ind-1);    
theta_new = theta2+(rcv_depth-depth2)/(depth2 - depth1)*(theta2-theta1);

depth_new = 5000;  
r_new = 100000;
%%
sinc = 0.1;
error = 0.05;
error_ang = 0.001;
weight = 1.4;
add_ang = 1;
op = 2;
while true 
    if ((abs(depth_new- rcv_depth) <= error) & (abs(r_new-geodesic) <= error)) & (add_ang <= error_ang) 
       break; 
    end
    % initial state parameters  
    [c0,~] = find_soundspeed(z0,SS);
    X0(idt) = 0;                                % time      
    X0(idz)= z0;                             % depth
    X0(idr)= r0;                            % range
    X0(idsh)= cos(pi*theta_new/180)./c0;    % horizontal slowness
    X0(idsv)= sin(pi*theta_new/180)./c0;    % vertical slowness
   
    % initial value
    snow = 0; 
    Xnow = X0;
    znow = z0;
    rnow = r0;
    X = Xnow;
    % loop until the horizonal range exceed geodesic
     op_d = op;
    while (rnow <= geodesic)
        % calculate the changes of parameters within this arc length interval
        f = diff_X(snow,Xnow,SS);
    
        % update the state vector 
        Xnow = Xnow + sinc*f;
        snow = snow+sinc;
        znow = Xnow(idz);
        rnow = Xnow(idr);
        % check if the ray hits the bottom
      
        if znow >= rcv_depth    
            op = 1;     % range optimize
            break;
        else
            op = 2;     % depth optimize
        end
         
    end
    
    if op_d ~= op
        fprintf(' op = %i\n',op)
        weight = weight+0.2 ;
        fprintf(' Weight = %f\n',weight)
        
    end
%     sinc = sinc/2;
    depth_new = znow;
    r_new = rnow;
    
    theta2 = theta_new;
    
    depth2 = depth_new;
    r2 = r_new;
    if op == 2
        add_ang = (rcv_depth-depth2)/(depth2 - depth1)*(theta2-theta1);
        theta_new = theta2+add_ang/weight;
    elseif op == 1
        add_ang = (geodesic-r2)/(r1 - r2)*(theta2-theta1);
        theta_new = theta2-add_ang/weight;
    end
    
    fprintf(' Depth Diff = %f\n Range Diff = %f\n Angle Adjust = %f \n',...
        (depth_new - rcv_depth),r_new-geodesic,theta_new-theta2)
    theta1 = theta2;
     
    depth1 = depth2;
    r1 = r2;
    
    if sinc >= 0.02
        sinc = sinc/1.5;
        fprintf('Arc Interval = %f\n',sinc)
    end
    disp('----------')
end

% return
est_time = Xnow(idt);
arclength = snow;
theta0 = theta_new;


end


%% ray tracing in planar coordinate
function [est_tt,theta0,rayarc_dist_tot] = ray_tracing_planar_const_SS(x_dist,td_z,SS)
% load sound speed profile
z_flt = SS.z;
SS_flt = SS.c;

% adjust the top layer according to the souce depth
ind = find(z_flt < td_z);
z_flt(ind) = [];
SS_flt(ind) = [];
z_flt = vertcat(td_z,z_flt);
[c0,dcdz] = find_soundspeed(td_z,SS);
SS_flt = vertcat(c0,SS_flt);

rcv_depth= z_flt(end);

%Find theta0 of acoustic transmission
% incident angle
% slant anle
slant_angle = atan(x_dist/rcv_depth);
theta0_all= slant_angle-10/180*pi:1*(pi/180):slant_angle+10/180*pi;   % all possible launch angles
theta0_all(find(theta0_all <= 0.0001)) = [];
theta0_all(find(theta0_all >= 89.5/180*pi)) = [];

r = [];    % horizontal range
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
 %% Supplementary function
 function SS = soundspeedprofile(ACO_depth)
ACO_lat = 22.738772;
ACO_lon = -158.006186;
ACO_depth = -ACO_depth;            % original depth MSL local u/e/n frame

fid=fopen('h306a0202.ctd');                       % 12 October 2018
D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f','headerlines',6));
fclose(fid);

pres_MSL = D(:,1);        %pressure (dbars)
temp = D(:,2);        %temperature
sal = D(:,3);         %salinity (Sp)
sal = gsw_SA_from_SP(sal,pres_MSL,ACO_lon,ACO_lat);

% Depths in the ellipsoidal coordinate; downward positive (-h)
% ACO_gh = geoidheight(ACO_lat,ACO_lon+360,'EGM96');
ACO_gh = 0;
z_hyd = -(ACO_depth+ACO_gh);

% create SS profile with the final depth at the hyd depth
% N_mean = 2.32;                                      % average geiod height over the area of coverage
N_mean = 0;
z_CTD = -(gsw_z_from_p(pres_MSL,ACO_lat) );         % relative to MSL
z_elps = z_CTD-N_mean;                              % from MSL to ellipsoid


% Truncate the lower bound incase the hyd depth is shallower
if z_elps(end) > z_hyd
    rm_ind = find(z_elps>z_hyd);
    z_elps(rm_ind) = [];
    z_CTD(rm_ind) = [];
    temp(rm_ind) = [];
    sal(rm_ind) = [];
else
%     Interpolate the lower bound in case the hyd depth is deeper
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

% calculate pressure
pres_CTD = gsw_p_from_z(-z_CTD,ACO_lat);
% 3. Creat Sound Speed Profile
SS_elps = gsw_sound_speed(sal,temp,pres_CTD);


%% Earth Flattening Transformation
R_e = 6371000;

% Dushaw and Colosi
% E = z_elps./R_e;
% z_flt = z_elps.*(1+(E./2)+((E.^2)./3));
% SS_flt = SS_elps.*(1+E+E.^2);

% Horizontal Sound Speed Gradient
a = 3.8e-4;
% dcdr = a*exp(-z_elps/600);
dcdr = 0;
% return Sound speed profile in the ellipsoidal coordinates
SS.z = z_elps;
SS.c = SS_elps;
SS.dcdr = dcdr;
 end
 
% function to find sound speed, vertical sound speed gradient, and horizontal
% sound speed gradient (c_z, dcdz,dcdr) at a given depth (znow) and
% a given sound speed profile (SS).
 function [c_z,dcdz,dcdr] = find_soundspeed(znow,SS)
 c = SS.c;
 z = SS.z;
 Hgrad = SS.dcdr;
 
 dcdz_all = diff(c)./diff(z);
 Hgrad_dz = diff(Hgrad)./diff(z);
 
 % find the interval where znow falls 
 ind = find(znow >= z,1,'last');
 c_lower = c(ind);
 Hgrad_lower = Hgrad(ind);
 
 del_z = znow - z(ind);
 Hgrad_dz = Hgrad_dz(ind);
 
 dcdz = dcdz_all(ind);
 dcdr = Hgrad_lower+del_z*Hgrad_dz;
 c_z = c_lower+del_z*dcdz;
 
 
 end
    
% differentials 
function dx =  diff_X(S,x,SS)
   idt = 1; idz = 2; idr = 3; idsh = 4; idsv = 5;
   % sound speed
   [c_z,dcdz,dcdr] = find_soundspeed(x(idz),SS);
%    dcdr = 0.00005;        % horizontal ss gradient
   dcdr = 0;
   
   dx(idt) = 1/c_z;
   dx(idz) = c_z*x(idsv);
   dx(idr) = c_z*x(idsh);
   dx(idsh) = -1/(c_z.^2)*dcdr;
   dx(idsv) = -1/(c_z.^2)*dcdz;
   
end
    
   