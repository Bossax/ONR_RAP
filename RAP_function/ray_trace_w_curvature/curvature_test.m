% comparoison between with and without cv correction
% x_dist = random walk until reah 25 km
% z_offset = random between -7 to -5
% z positive downward
clear 
close all
dist_now = 0;
z_offset0 = -4;
z_offset = -4;
% iii = 1;
slant_range1 = [];
% while dist_now <= 25000
% for iii = 1:2
dist_now = dist_now+ rand*300+200;
x_dist = 500:500:24550;
offset = (1:50)*.45;
% x_dist3 = x_dist2+offset;
for iii = 1:length(x_dist)
        %     x_dist(iii) = dist_now;
%     x_dist(1) = 10;
%     z_offset(iii) = z_offset0+(1.5*rand-7)
%     [arc_lengths1,tot_dist1(iii),theta01,SS1,z1,SS_HOT_avg1,surface_dist1,inc_r1,est_tt1,ray_angles1 ] = ray_tracing_const_SS_no_cv( x_dist(iii),z_offset(iii) );
  [arc_lengths2,tot_dist2(1),theta02(iii),SS2,z2,SS_HOT_avg2,surface_dist2,inc_r2,est_tt2(iii),ray_angles2] = ray_tracing_const_SS(x_dist(iii),z_offset(1) -2.32);
%   [arc_lengths3,tot_dist3(1),theta03(iii),SS3,z3,SS_HOT_avg3,surface_dist3,inc_r3,est_tt3(iii),ray_angles3] = ray_tracing_const_SS_new(x_dist(iii),z_offset(1) );
%   [arc_lengths4,tot_dist4(iii),theta04,SS4,z4,SS_HOT_avg4,surface_dist4,inc_r4,est_tt4(iii),ray_angles4,theta_loop] = ray_tracing_const_SS_gradient(x_dist(iii),z_offset(1) );
  
%     slant_range1(iii) = sqrt(x_dist(iii)^2 + (z3(end)+z_offset(iii))^2);
%     slant_range2(iii) = sqrt(dist_now^2 + (z4(end)+z_offset0)^2);
%     diff_slant_arc1(iii) = slant_range1(iii) - tot_dist4;
%     diff_slant_arc2(iii) = slant_range1(iii) - tot_dist1;
%     diff_slant_arc3(iii) = slant_range2(iii) - tot_dist4;
%     diff_lang(iii) = (theta02-theta03)/pi*180;
   
%      iii = iii+1;
   
end
%%
diff_t = est_tt2 - est_tt4;
diff_arc =  tot_dist2-tot_dist4;
%% plot
% % p = polyfit(x_dist/1000,diff_t*1000,1);
figure(1)
clf
scatter(x_dist/1000,diff_t*1000,'.')
grid on
xlabel('Range (km)')
ylabel('Time difference (ms)')
title('(Geoid - Ellipsoid)')
% xlim([0 25])

hold on
% yyaxis right
% scatter(x_dist/1000,diff_arc,'.')
% ylabel('Arc Length Difference (meter)')
% set(gca,'Yscale','log')
% txt = sprintf('%fRange  %f',p(1),p(2));
% % 
% figure(2)
% plot(x_dist/1000,diff_lang)
% grid on
% title('Launch Angle difference (degrees): (Model2 - Model3 with CV)')
% %%
% figure(3)
% clf
% plot(x_dist/1000,diff_slant_arc2)
% grid on
% title('Slant Range - Total Arc Length (no curvature correction)')
% axis tight
% xlabel('Range (km)')
% ylabel('meter')

%% 
figure(3)
clf
scatter(x_dist2/1000,(est_tt3-est_tt2)*1000,20,offset,'filled')
cbar = colorbar;
cbar.Label.String = 'Additional Distance (m)';
colormap jet
grid on
yticks([1:15])
ylabel('Travel Time Difference (ms)')
xlabel('Range (km)')
title('Simulation: Additional Surface Range as a function of range')
set(gca,'fontsize',13)
%%
function [arc_lengths,tot_dist,theta0,SS,z,SS_HOT_avg,surface_dist,inc_r,est_tt,ray_angles ] = ray_tracing_const_SS_no_cv( x_dist,z_offset )
% arc_length=
% tot_dist  =
% theta0    = launch angle
% SS        = Sound speed profile (with depth z)
% z         = Depth (negative/ below the transducer)
% SS_HOT    = 
% surface_dist = 
% est_tt    = estimated travel time
% ray angle = 
%fid=fopen('h294a0202.ctd');                       % June 2018
fid=fopen('h306a0202.ctd');                       % October 2018

D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f','headerlines',6));
fclose(fid);
        
%z_dist=4812.7+.669;      %Paroscientific pressure sensor (dbar) + hyd depth 4729.92 m 
z_dist = 4819.897;     % at 4,736.266 m June 2017

pres=D(:,1);        %pressure (dbars)
temp=D(:,2);        %temperature
sal=D(:,3);         %salinity


%% Fix z offset for transducer location (m=dbar at this level)
% Find depths above the onboard transducer depth, and discard them and
% associated data points
pres_offset = gsw_p_from_z(z_offset,22.738772);
A=find(pres_offset>=pres);
pres(A(1:end-1))=[];
pres_diff= pres_offset-pres(1);
pres(1)= pres_offset;      % set the first depth to be the transducer depth
temp(A(1:end-1))=[];
temp(1)=(temp(1)*(2-pres_diff)+temp(2)*pres_diff)/2;   % interpolate the first temp and salinity
sal(A(1:end-1))=[];
sal(1)=(sal(1)*(2-pres_diff)+sal(2)*pres_diff)/2;

%Remove/Add deeper depths
if pres(end)>z_dist
remove_points=find(pres>z_dist);
pres(remove_points)=[];
temp(remove_points)=[];
sal(remove_points)=[];
else
    
    add_len = floor((z_dist - pres(end))/2);
    pres(end+1:end+1+add_len)=[pres(end)+transpose(2:2:add_len*2) ;z_dist];
    temp(end+1:end+1+add_len)=ones(add_len+1,1)*temp(end);
    sal(end+1:end+1+add_len)=ones(add_len+1,1)*sal(end);
    %sal(end+1)=sal(end)+(sal(end)-sal(end-1));
end

%A priori SS field
SS=gsw_sound_speed(sal,temp,pres);

z=-gsw_z_from_p(pres,22.738772);

% %Correct for Earth Curvature
E_radius=6371000;        % epsilonray tracing APL code
% coordinate transformation
% E=z./E_radius;
% z=z.*(1+(E./2)+((E.^2)./3));
% %z(end)=hyd_depth;
% SS=SS.*(1+E+E.^2);

%Find theta0 of acoustic transmission
% theta0_all=0.1:0.01*(pi/180):89*(pi/180);
% incident angle
theta0_all=0.0001:1*(pi/180):89/180*pi;   % all possible launch angles
theta0_all(end)=89/180*pi;

% loop over all possible thetas 
% calculate a small horizontal distance within an interval betweenf two depths
for ii=1:length(theta0_all)
    
    %Ray parameter
    a=sin(theta0_all(ii))/SS(1);
    
    % incident angles of each layer
    theta(1) = theta0_all(ii);
    for i=2:length(SS)-1
        theta(i)=asin(a*SS(i));
    end    
    % SS gradient in each layer for the slowly-varying ss profile
%     for i=1:length(SS)-1
%         b(i)=(SS(i+1)-SS(i))/(z(i+1)-z(i));
%     end
    
    % horizontl distance
    for i=1:length(SS)-1
        r(i)=(z(i+1)-z(i))*tan(theta(i));
    end
    
    %Total x distance
    x_diff(ii)=sum(r);  % total horizontal distance for one possible launching theta
    
    
end


%Find min and max launch angle for Netwon-Raph method 

t_low=find(x_diff<x_dist);      % find all launch angles which yield x distance within the bound
theta_low=theta0_all(t_low(end));   % label the lowest possible launch angle (correspond to the steepest angle)
x_low=x_diff(t_low(end));           % The shortest ray arc length, lower bound of the sur dist

t_high=find(x_diff>x_dist);
theta_high=theta0_all(t_high(1));   % The upper bound of sur dist (the least steep angle)
x_high=x_diff(t_high(1));

% Intrapolation to find the optimum ray arc length of the direct path
theta_new=theta_low+((x_dist)/((x_low-x_high)/(theta_low-theta_high)))-((x_low)/((x_low-x_high)/(theta_low-theta_high)));

% We now have a range of launch angles which gives the nearest number of
% horizontal distance 
theta1=theta_low;
x1=x_low;           % lower bound of the horizontal distance

x_diff=1000000;

% Loop until x_diff is close to x_dist to 0.1 mm 
while abs(x_diff-x_dist)>0.0001
    
    %Ray parameter
    a=sin(theta_new)/SS(1);
    
    % incident angle
    theta(1) = theta_new;
    for i=2:length(SS)-1
        theta(i)=asin(a*SS(i));
    end
    
    r = [];
    % horizontal range in each layer
    for i=1:length(SS)-1
        r(i)=(z(i+1)-z(i))*tan(theta(i));
    end
    
    %Total x distance
    x_diff=sum(r);
    
    %Total distance and time traveled
    for i=1:length(SS)-1
        arc_distance(i)=(z(i+1)-z(i))/cos(theta(i));
        tt(i)=(z(i+1)-z(i))/(SS(i)*cos(theta(i)));
    end
    tot_dist_all=sum(arc_distance);
    
    arc_lengths_hold=arc_distance;
    surface_dist_hold=r;
    ray_angles_hold=theta(2:end);
    
    %Newton-Raphson method
    % Interpolate again to find a new launch angle which is closer to that
    % of the eigenray
    
    theta2=theta_new;
    
    
    x2=x_diff;
    
    theta_new=theta2+((x_dist)/((x2-x1)/(theta2-theta1)))-((x2)/((x2-x1)/(theta2-theta1)));
    
    % Set theta_low to theta 2 (new angle)
    theta1=theta2;
    % update nearest horizontal distance
    x1=x2;
    
end


%% packaging outputs
theta0=theta_new;                        %Launch angle
tot_dist=tot_dist_all;                 %Ray arc length
arc_lengths=arc_lengths_hold;           %arc length at each interval
surface_dist=surface_dist_hold;         %Surface distance
ray_angles=ray_angles_hold;             %ray angle at each interval

% incremental horizontal distance
inc_r = zeros(1,length(surface_dist)+1);

for i= 2:length(inc_r)
    inc_r(i) = inc_r(i-1)+surface_dist(i-1);
end

clear x_diff tot_dist_all arc_lengths_hold

est_tt=sum(tt);       %Estimated TT

%Vertically averaged SS for HOTS CTD
SS_HOT_avg=mean(SS);


%Make variables same size
% if length(A)>1
%     arc_lengths=horzcat(zeros(1,length(A)-1),arc_lengths);
%     surface_dist=horzcat(zeros(1,length(A)-1),surface_dist);
%     ray_angles=horzcat(ones(1,length(A)-1).*ray_angles(1,1),ray_angles);
% end

end

function [arc_lengths,tot_dist,theta0,SS,z,SS_HOT_avg,surface_dist,inc_r,est_tt,ray_angles ] = ray_tracing_const_SS( x_dist,z_offset )
% arc_length=
% tot_dist  =
% theta0    = launch angle
% SS        = Sound speed profile (with depth z)
% z         = Depth (negative/ below the transducer)
% SS_HOT    = 
% surface_dist = 
% est_tt    = estimated travel time
% ray angle = 
%fid=fopen('h294a0202.ctd');                       % June 2018
fid=fopen('h306a0202.ctd');                       % October 2018

D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f','headerlines',6));
fclose(fid);
icListen_lat=22.73911569;
icListen_lon=-158.006106601; 
%z_dist=4812.7+.669;      %Paroscientific pressure sensor (dbar) + hyd depth 4729.92 m 
z_dist = 4819.897;     % at 4,736.266 m June 2017

pres=D(:,1);        %pressure (dbars)
temp=D(:,2);        %temperature
sal=D(:,3);         %salinity
sal = gsw_SA_from_SP(sal,pres,icListen_lon,icListen_lat);

% Fix z offset for transducer location (m=dbar at this level)
% Find depths above the onboard transducer depth, and discard them and
% associated data points
pres_offset = gsw_p_from_z(z_offset,22.738772);
A=find(pres_offset>=pres);
pres(A(1:end-1))=[];
pres_diff= pres_offset-pres(1);
pres(1)= pres_offset;      % set the first depth to be the transducer depth
temp(A(1:end-1))=[];
temp(1)=(temp(1)*(2-pres_diff)+temp(2)*pres_diff)/2;   % interpolate the first temp and salinity
sal(A(1:end-1))=[];
sal(1)=(sal(1)*(2-pres_diff)+sal(2)*pres_diff)/2;

%Remove/Add deeper depths
if pres(end)>z_dist
    remove_points=find(pres>z_dist);
    pres(remove_points)=[];
    temp(remove_points)=[];
    sal(remove_points)=[];
else
    
    add_len = floor((z_dist - pres(end))/2);
    pres(end+1:end+1+add_len)=[pres(end)+transpose(2:2:add_len*2) ;z_dist];
    temp(end+1:end+1+add_len)=ones(add_len+1,1)*temp(end);
    sal(end+1:end+1+add_len)=ones(add_len+1,1)*sal(end);
    %sal(end+1)=sal(end)+(sal(end)-sal(end-1));
end

%A priori SS field
SS=gsw_sound_speed(sal,temp,pres);
z=-gsw_z_from_p(pres,22.738772);

% %Correct for Earth Curvature
% coordinate transformation
E_radius=6371000;        % epsilonray tracing APL code
E=z./E_radius;
z=z.*(1+(E./2)+((E.^2)./3));
SS=SS.*(1+E+E.^2);

% Find theta0 of acoustic transmission
% theta0_all=0.1:0.01*(pi/180):89*(pi/180);
% incident angle
theta0_all=0.0001:1*(pi/180):89/180*pi;   % all possible launch angles
theta0_all(end)=89/180*pi;

% loop over all possible thetas 
% calculate a small horizontal distance within an interval betweenf two depths
for ii=1:length(theta0_all)
    
    %Ray parameter
    a=sin(theta0_all(ii))/SS(1);
    
    % incident angles of each layer
    theta(1) = theta0_all(ii);
    for i=2:length(SS)-1
        theta(i)=asin(a*SS(i));
    end    
    % SS gradient in each layer for the slowly-varying ss profile
%     for i=1:length(SS)-1
%         b(i)=(SS(i+1)-SS(i))/(z(i+1)-z(i));
%     end
    
    % horizontl distance
    for i=1:length(SS)-1
        r(i)=(z(i+1)-z(i))*tan(theta(i));
    end
    
    %Total x distance
    x_diff(ii)=sum(r);  % total horizontal distance for one possible launching theta
    
    
end


%Find min and max launch angle for Netwon-Raph method 

t_low=find(x_diff<x_dist);      % find all launch angles which yield x distance within the bound
theta_low=theta0_all(t_low(end));   % label the lowest possible launch angle (correspond to the steepest angle)
x_low=x_diff(t_low(end));           % The shortest ray arc length, lower bound of the sur dist

t_high=find(x_diff>x_dist);
theta_high=theta0_all(t_high(1));   % The upper bound of sur dist (the least steep angle)
x_high=x_diff(t_high(1));

% Intrapolation to find the optimum ray arc length of the direct path
theta_new=theta_low+((x_dist)/((x_low-x_high)/(theta_low-theta_high)))-((x_low)/((x_low-x_high)/(theta_low-theta_high)));

% We now have a range of launch angles which gives the nearest number of
% horizontal distance 
theta1=theta_low;
x1=x_low;           % lower bound of the horizontal distance

x_diff=1000000;

% Loop until x_diff is close to x_dist to 0.1 mm 
while abs(x_diff-x_dist)>0.0001
    
    %Ray parameter
    a=sin(theta_new)/SS(1);
    
    % incident angle
    theta(1) = theta_new;
    for i=2:length(SS)-1
        theta(i)=asin(a*SS(i));
    end
    
    r = [];
    % horizontal range in each layer
    for i=1:length(SS)-1
        r(i)=(z(i+1)-z(i))*tan(theta(i));
    end
    
    %Total x distance
    x_diff=sum(r);
    
    %Total distance and time traveled
    for i=1:length(SS)-1
        arc_distance(i)=(z(i+1)-z(i))/cos(theta(i));
        tt(i)=(z(i+1)-z(i))/(SS(i)*cos(theta(i)));
    end
    tot_dist_all=sum(arc_distance);
    
    arc_lengths_hold=arc_distance;
    surface_dist_hold=r;
    ray_angles_hold=theta(2:end);
    
    %Newton-Raphson method
    % Interpolate again to find a new launch angle which is closer to that
    % of the eigenray
    
    theta2=theta_new;
    
    
    x2=x_diff;
    
    theta_new=theta2+((x_dist)/((x2-x1)/(theta2-theta1)))-((x2)/((x2-x1)/(theta2-theta1)));
    
    % Set theta_low to theta 2 (new angle)
    theta1=theta2;
    % update nearest horizontal distance
    x1=x2;
    
end


%% packaging outputs
theta0=theta_new;                        %Launch angle
tot_dist=tot_dist_all;                 %Ray arc length
arc_lengths=arc_lengths_hold;           %arc length at each interval
surface_dist=surface_dist_hold;         %Surface distance
ray_angles=ray_angles_hold;             %ray angle at each interval

% incremental horizontal distance
inc_r = zeros(1,length(surface_dist)+1);

for i= 2:length(inc_r)
    inc_r(i) = inc_r(i-1)+surface_dist(i-1);
end

clear x_diff tot_dist_all arc_lengths_hold

est_tt=sum(tt);       %Estimated TT

%Vertically averaged SS for HOTS CTD
SS_HOT_avg=mean(SS);


%Make variables same size
% if length(A)>1
%     arc_lengths=horzcat(zeros(1,length(A)-1),arc_lengths);
%     surface_dist=horzcat(zeros(1,length(A)-1),surface_dist);
%     ray_angles=horzcat(ones(1,length(A)-1).*ray_angles(1,1),ray_angles);
% end

end

function [arc_lengths,tot_dist,theta0,SS_fl,z_fl,SS_HOT_avg,surface_dist,inc_r,est_tt,ray_angles ] = ray_tracing_const_SS_new( x_dist,z_offset )
% arc_length=
% tot_dist  =
% theta0    = launch angle
% SS        = Sound speed profile (with depth z)
% z         = Depth (negative/ below the transducer)
% SS_HOT    = 
% surface_dist = 
% est_tt    = estimated travel time
% ray angle = 
%fid=fopen('h294a0202.ctd');                       % June 2018
fid=fopen('h306a0202.ctd');                       % October 2018
icListen_lat=22.73911569;
icListen_lon=-158.006106601;   
D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f','headerlines',6));
fclose(fid);
        
%z_dist=4812.7+.669;      %Paroscientific pressure sensor (dbar) + hyd depth 4729.92 m 
% z_dist = 4819.897;     % at 4,736.266 m June 2017 % Orthometric Height
z_hyd_MSL = 4736.266;

pres_MSL=D(:,1);        %pressure (dbars)
temp=D(:,2);        %temperature
sal=D(:,3);         %salinity (Sp)
sal = gsw_SA_from_SP(sal,pres_MSL,icListen_lon,icListen_lat);

%% Fix z offset for transducer location (m=dbar at this level)
% Find depths above the onboard transducer depth, and discard them and associated data points
% transducer depth in ellipsoidal height
% pressure is equivalenth to depth relative to MSL

N_mean = 2.32;  % median geoid height 
z_hyd_elps = z_hyd_MSL-N_mean;
pres_hyd_elps = gsw_p_from_z(-z_hyd_elps,icListen_lat);

z_CTD_MSL = -gsw_z_from_p(pres_MSL,22.738772);
z_CTD_elps = z_CTD_MSL-N_mean;
pres_elps = gsw_p_from_z(-z_CTD_elps,icListen_lat);

A=find(z_CTD_elps<=-z_offset);
pres_elps(A(1:end-1))=[];

z_diff= z_offset-z_CTD_elps(1);

% set the first depth to be the transducer depth
pres_elps(1) = gsw_p_from_z(z_offset,icListen_lat);

temp(A(1:end-1))=[];
temp(1)=(temp(1)*(2-z_diff)+temp(2)*z_diff)/2;   % interpolate the first temp and salinity
sal(A(1:end-1))=[];
sal(1)=(sal(1)*(2-z_diff)+sal(2)*z_diff)/2;

%Remove/Add deeper depths
if pres_elps(end)>pres_hyd_elps
    
    remove_ind=find(pres_elps>pres_hyd_elps);
    pres_elps(remove_ind)=[];
    temp(remove_ind)=[];
    sal(remove_ind)=[];
else
    
    add_len = floor((pres_hyd_elps - pres_elps(end))/2);
    pres_elps(end+1:end+1+add_len)=[pres_elps(end)+transpose(2:2:add_len*2) ;pres_hyd_elps];
    temp(end+1:end+1+add_len)=ones(add_len+1,1)*temp(end);
    sal(end+1:end+1+add_len)=ones(add_len+1,1)*sal(end);
    
end

%A priori SS field in Cartesian coordinate relative to reference ellipsoid
SS_elps=gsw_sound_speed(sal,temp,pres_elps);
z_elps = -gsw_z_from_p(pres_elps,icListen_lat);


% %Correct for Earth Curvature
Re=6371000;        % epsilonray tracing APL code
Rc = Re;
% coordinate transformation
epsilon = z_elps./Re;
% z_r =z.*(1-epsilon./2+epsilon.*epsilon/6);
% SS_r = SS.*(1-z_r./Re);
z_r = z_elps.*(1+epsilon/2+(epsilon.^2)/3);
SS_r = SS_elps.*(1+epsilon+epsilon.^2);
% SS_r = SS;

% Find theta0 of acoustic transmission
% theta0_all=0.1:0.01*(pi/180):89*(pi/180);
% incident angle
theta0_all=[0.0001:1*(pi/180):86]/180*pi;   % all possible launch angles
theta0_all(end)=87/180*pi;

% From Peter M Shearer
r = Rc-z_elps;
z_fl = -Rc*log(r./Rc);
SS_fl = Rc./r.*SS_elps;

% loop over all possible thetas 
% calculate a small horizontal distance within an interval betweenf two depths
for ii=1:length(theta0_all)
    
    %Ray parameter
    a=sin(theta0_all(ii))/SS_fl(1);
    
 
    % Snell's Law incident angles of each layer
    theta(1) = theta0_all(ii);
    theta2(1) =theta0_all(ii);
    
    for iii=2:length(SS_fl)-1
        theta(iii)=asin(a*SS_fl(iii));
        theta2(iii) = asin(a*SS_elps(iii));
    end    
    
    % horizontl distance
    for iii=1:length(SS_fl)-1
        x(iii)=(z_fl(iii+1)-z_fl(iii))*a*SS_fl(iii)/sqrt(1-(a*SS_fl(iii))^2);
        x2(iii)=(z_elps(iii+1)-z_elps(iii))*a*SS_elps(iii)/sqrt(1-(a*SS_elps(iii))^2);
    end
  
    %{
    ang_dis = [];
    x_inc = [];
    % angular displacement
    for iii=1:length(SS_r)-1
        ang_dis(iii)=log(r(iii)/r(iii+1))*a*SS_r2(iii)/sqrt(1-(a*SS_r2(iii)));
        x_inc(iii) = real(ang_dis(iii))*r(iii);
    end
      %}
    
    %Total x distance
    X_dist(ii)=sum(x);  % total horizontal distance for one possible launching theta
    X_dist2(ii)=sum(x2);  
    %{
    % Total angular distance
    ang_dist_tot(ii) = real(sum(ang_dis));
    %}
    
   
end


%Find min and max launch angle for Netwon-Raph method 

t_low=find(X_dist<x_dist);      % find all launch angles which yield x distance within the bound
theta_low=theta0_all(t_low(end));   % label the lowest possible launch angle (correspond to the steepest angle)
x_low=X_dist(t_low(end));           % The shortest ray arc length, lower bound of the sur dist

t_high=find(X_dist>x_dist);
theta_high=theta0_all(t_high(1));   % The upper bound of sur dist (the least steep angle)
x_high=X_dist(t_high(1));

% Intrapolation to find the optimum ray arc length of the direct path
theta_new=theta_low+((x_dist)/((x_low-x_high)/(theta_low-theta_high)))-((x_low)/((x_low-x_high)/(theta_low-theta_high)));

% We now have a range of launch angles which gives the nearest number of
% horizontal distance 
theta1=theta_low;
x1=x_low;           % lower bound of the horizontal distance

X_dist=1000000;

% Loop until x_diff is close to x_dist to 0.1 mm 
while abs(X_dist-x_dist)>0.0001
    
    %Ray parameter
    a=sin(theta_new)/SS_fl(1);
    
    % incident angle
    theta(1) = theta_new;
    for iii=2:length(SS_fl)-1
        theta(iii)=asin(a*SS_fl(iii));
    end
   
    %{
     ang_dis = [];
     tt = [];
     arc_distance = [];
     x_inc = [];
    % angular displacement
    for iii=1:length(SS_r)-1
        ang_dis(iii)=log(r(iii)/r(iii+1))*a*SS_r2(iii)/sqrt(1-(a*SS_r2(iii)));
        arc_distance(iii) = Rc*log(r(iii)/r(iii+1))/sqrt(1-(a*SS_r2(iii))^2);
        tt(iii)=log(r(iii)/r(iii+1))*r(iii)/(SS_r2(iii)*sqrt(1-(a*SS_r2(iii))^2));
        x_inc(iii) = real(ang_dis(iii))*r(iii);
    end
    %Total x distance
    x_diff=sum(x_inc);
    
    %}
    x = [];
    tt = [];
    arc_distance = [];
    for iii=1:length(SS_fl)-1
        
        arc_distance(iii)=(z_fl(iii+1)-z_fl(iii))/cos(theta(iii));
        tt(iii)= (z_fl(iii+1)-z_fl(iii))/(SS_fl(iii)*sqrt(1-(a*SS_fl(iii))^2));
        x(iii)=(z_fl(iii+1)-z_fl(iii))*a*SS_fl(iii)/sqrt(1-(a*SS_fl(iii))^2);
    end
    
    L_length =sum(arc_distance);
    X_dist = sum(x);
    travel_time = sum(tt);
    %Newton-Raphson method
    % Interpolate again to find a new launch angle which is closer to that
    % of the eigenray
    
    theta2=theta_new;
    
    
    x2=X_dist;
    
    theta_new=theta2+((x_dist)/((x2-x1)/(theta2-theta1)))-((x2)/((x2-x1)/(theta2-theta1)));
    
    % Set theta_low to theta 2 (new angle)
    theta1=theta2;
    % update nearest horizontal distance
    x1=x2;
    
end


%% packaging outputs
theta0=theta_new;                        %Launch angle
tot_dist =L_length;                 %Ray arc length
arc_lengths= arc_distance;           %arc length at each interval
surface_dist=x;         %Surface distance
% ray_angles=ray_angles_hold;             %ray angle at each interval
ray_angles = [];

% incremental horizontal distance
inc_r = zeros(1,length(surface_dist)+1);

for iii= 2:length(inc_r)
    inc_r(iii) = inc_r(iii-1)+surface_dist(iii-1);
end

clear x_diff tot_dist_all arc_lengths_hold

est_tt=real(sum(tt));       %Estimated TT

%Vertically averaged SS for HOTS CTD
SS_HOT_avg=mean(SS_fl);


%Make variables same size
% if length(A)>1
%     arc_lengths=horzcat(zeros(1,length(A)-1),arc_lengths);
%     surface_dist=horzcat(zeros(1,length(A)-1),surface_dist);
%     ray_angles=horzcat(ones(1,length(A)-1).*ray_angles(1,1),ray_angles);
% end

end

function [arc_lengths,tot_dist,theta0,SS_fl,z_fl,SS_HOT_avg,surface_dist,inc_r,est_tt,ray_angles,theta_loop ] = ray_tracing_const_SS_gradient( x_dist,z_offset )
% arc_length= layerwise arc length
% tot_dist  = total arc length of the eigen ray
% theta0    = launch angle
% SS        = Sound speed profile (with depth z)
% z         = Depth (negative/ below the transducer)
% SS_HOT    = Average sound speed
% surface_dist = Total horizontal distance traveled by the ray
% est_tt    = estimated travel time
% ray angle = eigen angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fid=fopen('h294a0202.ctd');                       % June 2018
fid=fopen('h306a0202.ctd');                       % October 2018
icListen_lat=22.73911569;
icListen_lon=-158.006106601;   
D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f','headerlines',6));
fclose(fid);
        
%z_dist=4812.7+.669;      %Paroscientific pressure sensor (dbar) + hyd depth 4729.92 m 
z_dist = 4819.897;     % at 4,736.266 m June 2017
z_hyd_MSL = 4736.266;

pres_MSL=D(:,1);        %pressure (dbars)
temp=D(:,2);        %temperature
sal=D(:,3);         %salinity
sal = gsw_SA_from_SP(sal,pres_MSL,icListen_lon,icListen_lat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sound speed profile preparation

%% Fix z offset for transducer location (m=dbar at this level)
% Find depths above the onboard transducer depth, and discard them and associated data points
% transducer depth in ellipsoidal height
% pressure is equivalenth to depth relative to MSL

N_mean = 2.32;  % median geoid height 
z_hyd_elps = z_hyd_MSL-N_mean;
pres_hyd_elps = gsw_p_from_z(-z_hyd_elps,icListen_lat);

z_CTD_MSL = -gsw_z_from_p(pres_MSL,22.738772);
z_CTD_elps = z_CTD_MSL-N_mean;
pres_elps = gsw_p_from_z(-z_CTD_elps,icListen_lat);

A=find(z_CTD_elps<=-z_offset);
pres_elps(A(1:end-1))=[];

z_diff= z_offset-z_CTD_elps(1);

% set the first depth to be the transducer depth
pres_elps(1) = gsw_p_from_z(z_offset,icListen_lat);

temp(A(1:end-1))=[];
temp(1)=(temp(1)*(2-z_diff)+temp(2)*z_diff)/2;   % interpolate the first temp and salinity
sal(A(1:end-1))=[];
sal(1)=(sal(1)*(2-z_diff)+sal(2)*z_diff)/2;

%Remove/Add deeper depths
if pres_elps(end)>pres_hyd_elps
    
    remove_ind=find(pres_elps>pres_hyd_elps);
    pres_elps(remove_ind)=[];
    temp(remove_ind)=[];
    sal(remove_ind)=[];
else
    
    add_len = floor((pres_hyd_elps - pres_elps(end))/2);
    pres_elps(end+1:end+1+add_len)=[pres_elps(end)+transpose(2:2:add_len*2) ;pres_hyd_elps];
    temp(end+1:end+1+add_len)=ones(add_len+1,1)*temp(end);
    sal(end+1:end+1+add_len)=ones(add_len+1,1)*sal(end);
    
end

%A priori SS field in Cartesian coordinate relative to reference ellipsoid
SS_elps=gsw_sound_speed(sal,temp,pres_elps);
z_elps = -gsw_z_from_p(pres_elps,icListen_lat);


% %Correct for Earth Curvature
Re=6371000;        % epsilonray tracing APL code
Rc = Re;
% coordinate transformation
epsilon = z_elps./Re;
% z_r =z.*(1-epsilon./2+epsilon.*epsilon/6);
% SS_r = SS.*(1-z_r./Re);
z_r = z_elps.*(1+epsilon/2+(epsilon.^2)/3);
SS_r = SS_elps.*(1+epsilon+epsilon.^2);


% Find theta0 of acoustic transmission
% theta0_all=0.1:0.01*(pi/180):89*(pi/180);
% incident angle
theta0_all=0.0001:1*(pi/180):85/180*pi;   % all possible launch angles
theta0_all(end)=85/180*pi;

% From Peter M Shearer
r = Re-z_elps;
z_fl = Re*log(Re./r);
SS_fl = Re./r.*SS_elps;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the eigen ray
%Find theta0 of acoustic transmission

% incident angle
theta0_all=0.0001:1*(pi/180):89/180*pi;   % all possible launch angles
theta0_all(end)=88/180*pi;

% loop over all possible thetas 
% calculate a small horizontal distance within an interval betweenf two depths

for ii=1:length(theta0_all)
    r = [];
    b = [];
    
    %Ray parameter
    a=sin(theta0_all(ii))/SS_fl(1);
    
    % SS_fl gradient in each layer for the slowly-varying SS_fl profile
    for i= 1:length(SS_fl) - 1
        b(i)=(SS_fl(i+1)-SS_fl(i))/(z_fl(i+1)-z_fl(i));
     end

    % horiz_flontl distance
    for j= 1:length(SS_fl) - 1
        v = SS_fl(j+1)/SS_fl(j);
        h = (sqrt(1-(a*SS_fl(j+1))^2)+1)/(sqrt(1-(a*SS_fl(j))^2)+1);
        r(j) = (1/(a*b(j)))*(sqrt(1-(a*SS_fl(j))^2) - sqrt(1-(a*SS_fl(j+1))^2));
       
    end   

    
    %Total horiz_flontal distance to reach the hydrophone
    x_diff(ii)=sum(r);  
    
    
end


%Find min and max launch angle for Netwon-Raph method 

t_low=find(x_diff<x_dist);      % find all launch angles which yield x distance within the bound
theta_low=theta0_all(t_low(end));   % label the lowest poSS_flible launch angle (correspond to the steepest angle)
x_low=x_diff(t_low(end));           % The shortest ray arc length, lower bound of the sur dist

t_high=find(x_diff>x_dist);
theta_high=theta0_all(t_high(1));   % The upper bound of sur dist (the least steep angle)
x_high=x_diff(t_high(1));

% Intrapolation to find the optimum ray arc length of the direct path
theta_new=theta_low+((x_dist)/((x_low-x_high)/(theta_low-theta_high)))-((x_low)/((x_low-x_high)/(theta_low-theta_high)));

% We now have a range of launch angles which gives the nearest number of
% horiz_flontal distance 
theta1=theta_low;
x1=x_low;           % lower bound of the horiz_flontal distance

x_diff=1000000;
theta_loop = theta_new;
% Loop until x_diff is close to x_dist to 0.1 mm 
while abs(x_diff-x_dist)>0.0001
    r = [];         % horiz_flontal distance in each layer
    t = [];         % elapsed time in each layer
    b = [];         % SS_fl gradient in each layer
    
    %Ray parameter
    a=sin(theta_new)/SS_fl(1);
    
    % incident angle at each interface
    theta(1) = theta_new;
    for i= 2:length(SS_fl)
        theta(i)=asin(a*SS_fl(i));
    end
    
    % SS_fl gradient in each layer for the slowly-varying SS_fl profile
    for i= 1:length(SS_fl)-1
        b(i)=(SS_fl(i+1)-SS_fl(i))/(z_fl(i+1)-z_fl(i));
    end

% horiz_flontl distance, arc length and travel time
    for j= 1:length(SS_fl)-1
        v = SS_fl(j+1)/SS_fl(j);
        h = (sqrt(1-(a*SS_fl(j+1))^2)+1)/(sqrt(1-(a*SS_fl(j))^2)+1);
       t(j) = (1/b(j))*(log(v)-log(h));
       r(j) = (1/(a*b(j)))*(sqrt(1-(a*SS_fl(j))^2) - sqrt(1-(a*SS_fl(j+1))^2));
       arc_distance(j)= 1/abs(a*b(j))*abs(theta(j+1) - theta(j));
    end   
   
    % Total ditances
    x_diff=sum(r);                  % horiz_flontal distance
    tot_dist_all=sum(arc_distance); % arc length
    
    % Layerwise distance/ angles
    arc_lengths_hold=arc_distance;
    surface_dist_hold=r;
    ray_angles_hold=theta(2:end);
    
    % %%%% Newton-Raphson method
    % Interpolate again to find a new launch angle which is closer to that
    % of the eigenray
    
    theta2=theta_new;
    
    
    x2=x_diff;
    
    theta_new=theta2+((x_dist)/((x2-x1)/(theta2-theta1)))-((x2)/((x2-x1)/(theta2-theta1)));
    
    theta_loop(end+1) = theta_new;
    
    % Set theta_low to theta 2 (new angle)
    theta1=theta2;
    % update nearest horiz_flontal distance
    x1=x2;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% packaging outputs
theta0=theta_new;                        %Launch angle
tot_dist=tot_dist_all;                 %Ray arc length
arc_lengths=arc_lengths_hold;           %arc length at each interval
surface_dist=surface_dist_hold;         %Surface distance
ray_angles=ray_angles_hold;             %ray angle at each interval

% incremental horiz_flontal distance
inc_r = zeros(1,length(surface_dist)+1);

for i= 2:length(inc_r)
    inc_r(i) = inc_r(i-1)+surface_dist(i-1);
end

clear x_diff tot_dist_all arc_lengths_hold

est_tt=sum(t);       %Estimated TT

%Vertically averaged SS_fl for HOTS CTD
SS_HOT_avg=mean(SS_fl);
end

