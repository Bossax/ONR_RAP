clear
close all

ACO_lat= 22.738772;                  % June 2017
ACO_lon= -158.006186;                % June 2017
ACO_depth = -4733.946;            % at 4,736.266 m June 2017
td_h = -6;
x_dist = 300:500:25800;

% Radial
n_test = length(x_dist);
angle = 90;
L = 0;
tx_lon = [];
tx_lat = [];

for iii = 1:n_test
    L = x_dist(iii);
    [tx_lon(iii),tx_lat(iii),~] = m_fdist(ACO_lon,ACO_lat,angle,L);
    tx_lon(iii) = tx_lon(iii)-360;
    
end
%% Oct
diff_arc = [];
diff_tt = [];
alpha = 90;
for iii =1:length(tx_lat)
    [arc_lengths_o,SS_o,z_o,est_tto(iii),to] = ray_trace(x_dist(iii),td_h,tx_lat(iii),alpha,ACO_lat,ACO_lon,ACO_depth,'Oct');
    [arc_lengths_n,SS_n,z_n,est_ttn(iii),tn]= ray_trace(x_dist(iii),td_h,tx_lat(iii),alpha,ACO_lat,ACO_lon,ACO_depth,'Nov');
%     [arc_lengths,tot_dist,theta0,SS_flt,z_flt,SS_HOT_avg,surface_dist,est_tt,tt,ray_angles,R_alpha]
    diff_arc(iii,:) = arc_lengths_n-arc_lengths_o;
    diff_tt(iii,:) = real(tn)-real(to);
    
end

%% Plot
load('EOF_SS.mat')
mode1 = EOF_SS.mode1;
mode2 = EOF_SS.mode2;
mode3 = EOF_SS.mode3;
mode4 = EOF_SS.mode4;
%%
figure(1)
clf
subplot(1,4,1)
plot(SS_n-SS_o,z_o,'LineWidth',3)
hold on
plot(mode1*5,z_o(1:length(mode1)),'-r')
set(gca,'YDir','reverse')
grid on
xlim([-10 10])
ylim([0 4700])
subplot(1,4,2)
plot(SS_n-SS_o,z_o,'LineWidth',3)
hold on
plot(mode2*5,z_o(1:length(mode1)),'-k')
set(gca,'YDir','reverse')
grid on
xlim([-10 10])
ylim([0 4700])

subplot(1,4,3)
plot(SS_n-SS_o,z_o,'LineWidth',3)
hold on
plot(mode3*5,z_o(1:length(mode1)),'-g')
set(gca,'YDir','reverse')
grid on
xlim([-10 10])
ylim([0 4700])

subplot(1,4,4)
plot(SS_n-SS_o,z_o,'LineWidth',3)
hold on
plot(mode4*5,z_o(1:length(mode1)),'-m')
set(gca,'YDir','reverse')
grid on
xlim([-10 10])
ylim([0 4700])
%%
figure(2)
clf
plot(diff_tt(10,:)*10000000,z_o(1:end-1))
set(gca,'YDir','reverse')
grid on
hold on
plot(SS_n-SS_o,z_o,'LineWidth',2)
function [arc_lengths,SS_flt,z_flt,est_tt,tt] = ray_trace(S,td_h,tx_lat,alphas,ACO_lat,ACO_lon,ACO_depth,month)
%%%%% Functions %%%%%%%%%%
% 1. ctd file
% 2. sound speed calculation function gsw_sound_speed
% 3. Earth Flatenning Transformation
% 4. enthalpy
%%%%% inputs %%%%%%%%%%%%
% 1. Geodesic distance S
% 2. up in local horizon coordinate
% 3. transducer ellipsoidal height = td_elps_h 
% 4. azimuth of the receiver relative to the transducer
% 5. Hydrophone geodetic positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference Ellipsoid
ref_ellipsoid =  referenceEllipsoid('WGS84');
major_ax = ref_ellipsoid.SemimajorAxis;
minor_ax = ref_ellipsoid.SemiminorAxis;
ecc = ref_ellipsoid.Eccentricity;

% remove bad data
if td_h >=0
    td_h = 0;
end
%% 1. CTD data (MSL)
if month == 'Oct'
fid=fopen('h306a0202.ctd');                       % 12 October 2018
elseif month == 'Nov'
fid=fopen('h307a0202.ctd');                       % 16 November 2018
end

D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f','headerlines',6));
fclose(fid);

pres_MSL = D(:,1);        %pressure (dbars)
temp = D(:,2);        %temperature
sal = D(:,3);         %salinity (Sp)
sal = gsw_SA_from_SP(sal,pres_MSL,ACO_lon,ACO_lat);

% Depths in the ellipsoidal coordinate 
z_td = -td_h;
z_hyd = -ACO_depth;      


%% 2 truncate CTD data to fit the range between the receiver depth and the source depth
% on the ellipsoidal surface Orthometrc height

N_mean = 2.32;      % average geiod height over the area of coverage 
z_CTD = -(gsw_z_from_p(pres_MSL,ACO_lat) + N_mean); % from MSL to ellipsoid

% 2.1 truncate the upper bound using transducer depth
rm_ind = find(z_CTD<=z_td);
z_CTD(rm_ind(1:end-1)) = [];
temp(rm_ind(1:end-1)) = [];
sal(rm_ind(1:end-1)) = [];

% set the first depth to be the transducer depth
z_diff = z_td - z_CTD(1);
temp(1) = temp(1)+(temp(2)-temp(1))*(z_diff/(z_CTD(2)-z_CTD(1)));   % interpolate the first temp and salinity
sal(1)=sal(1)+(sal(2)-sal(1))*(z_diff/(z_CTD(2)-z_CTD(1)));   
z_CTD(1) = z_td;

% 2.2 truncate the lower bound
if z_CTD(end) > z_hyd
    rm_ind = find(z_CTD>z_hyd);
    z_CTD(rm_ind) = [];
    temp(rm_ind) = [];
    sal(rm_ind) = [];
else
    ly_thk = -mean(z_CTD(end-5:end-1) - z_CTD(end-4:end));
    add_layer = floor((z_hyd -z_CTD(end))/ly_thk);
    z_CTD = vertcat(z_CTD,z_CTD(end)+ly_thk*(1:add_layer)');
    z_CTD(end+1) = z_hyd;
    
    temp_gradient = temp(end-1)-temp(end);
    temp = vertcat(temp,temp(end)+temp_gradient*(1:add_layer+1)');
    sal_gradient = sal(end-1)-sal(end);
    sal = vertcat(sal,sal(end)+sal_gradient*(1:add_layer+1)');
end

pres_CTD = gsw_p_from_z(-z_CTD,ACO_lat);
%% 3. Creat Sound Speed Profile 
SS_elps = gsw_sound_speed(sal,temp,pres_CTD);
z_elps = z_CTD;

%% Earth Flattening Transformation
% Chadwell 2010
M_r = major_ax*(1-ecc^2)/(1-ecc^2*sin(tx_lat/180*pi)^2)^1.5;
N_r = major_ax/sqrt(1-ecc^2*sin(tx_lat/180*pi)^2);
alphas = alphas/180*pi;
R_alpha = (cos(alphas)^2/M_r+sin(alphas)^2/N_r)^-1;

z_flt = -R_alpha.*log(R_alpha./(R_alpha-(-z_elps)));
SS_flt = SS_elps.*(R_alpha./(R_alpha-(-z_elps)));


r = [];
%Find theta0 of acoustic transmission
theta0_all=0.0001:1*(pi/180):89/180*pi;   % all possible launch angles
theta0_all(end)=89/180*pi;

% loop over all possible thetas 
% calculate a small horizontal distance within an interval betweenf two depths
for jj=1:length(theta0_all)
    
    %Ray parameter
    a=sin(theta0_all(jj))/SS_flt(1);
    
    % incident angles of each layer
    theta(1) = theta0_all(jj);
    
    for ii=2:length(SS_flt)-1
        theta(ii)=asin(a*SS_flt(ii));
    end    
    % SS gradient in each layer for the slowly-varying ss profile
%     for i=1:length(SS)-1
%         b(i)=(SS(i+1)-SS(i))/(z(i+1)-z(i));
%     end
    
    % horizontl distance
    for ii=1:length(SS_flt)-1
        r(ii)=(z_flt(ii+1)-z_flt(ii))*tan(theta(ii));
    end
    
    %Total x distance
    x_diff(jj)=sum(r);  % total horizontal distance for one possible launching theta
    
    
end


%Find min and max launch angle for Netwon-Raph method 

t_low=find(x_diff<S);      % find all launch angles which yield x distance within the bound
theta_low=theta0_all(t_low(end));   % label the lowest possible launch angle (correspond to the steepest angle)
x_low=x_diff(t_low(end));           % The shortest ray arc length, lower bound of the sur dist

t_high=find(x_diff>S);
theta_high=theta0_all(t_high(1));   % The upper bound of sur dist (the least steep angle)
x_high=x_diff(t_high(1));

% Intrapolation to find the optimum ray arc length of the direct path
theta_new=theta_low+((S)/((x_low-x_high)/(theta_low-theta_high)))-((x_low)/((x_low-x_high)/(theta_low-theta_high)));

% We now have a range of launch angles which gives the nearest number of
% horizontal distance 
theta1=theta_low;
x1=x_low;           % lower bound of the horizontal distance

x_diff=1000000;

% Loop until x_diff is close to x_dist to 0.1 mm 
while abs(x_diff-S)>0.0001
    
    %Ray parameter
    a=sin(theta_new)/SS_flt(1);
    
    % incident angle
    theta(1) = theta_new;
    for ii=2:length(SS_flt)-1
        theta(ii)=asin(a*SS_flt(ii));
    end
    
    r = [];
    % horizontal range in each layer
    for ii=1:length(SS_flt)-1
        r(ii)=(z_flt(ii+1)-z_flt(ii))*tan(theta(ii));
    end
    
    %Total x distance
    x_diff=sum(r);
    
    %Total distance and time traveled
    for ii=1:length(SS_flt)-1
        arc_distance(ii)=(z_flt(ii+1)-z_flt(ii))/cos(theta(ii));
        tt(ii)=(z_flt(ii+1)-z_flt(ii))/(SS_flt(ii)*cos(theta(ii)));
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
    
    theta_new=theta2+((S)/((x2-x1)/(theta2-theta1)))-((x2)/((x2-x1)/(theta2-theta1)));
    
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

for ii= 2:length(inc_r)
    inc_r(ii) = inc_r(ii-1)+surface_dist(ii-1);
end


est_tt=sum(tt);       %Estimated TT

%Vertically averaged SS for HOTS CTD
SS_HOT_avg=mean(SS_flt);



% %Make variables same size
% if length(A)>1
%     arc_lengths=horzcat(zeros(1,length(A)-1),arc_lengths);
%     surface_dist=horzcat(zeros(1,length(A)-1),surface_dist);
%     ray_angles=horzcat(ones(1,length(A)-1).*ray_angles(1,1),ray_angles);
% end




end

