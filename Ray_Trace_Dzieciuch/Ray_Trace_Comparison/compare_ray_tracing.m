%% Compare Ray Tracing (Varamo and Dzieciuch)
clearvars
close all

seafloor_z = 4730;
lat=22.738894;

z_c=(0:seafloor_z);                          %DEPTH PROFILE (m) 
SS_c=cssprofile(z_c,lat);
 
% load('SS_profile_Oct.mat')
% SS_c = SS;
% z_c = z;

% smoothing by MA order 5
% order = 10;
% SS_sm(1:floor(order/2)) = SS_c(1:floor(order/2));
% for ii = floor(order/2)+1:length(z_c)-floor(order/2)-1
%     SS_sm(ii) = mv_avg_SS(order,SS_c,ii);
% end
% 
% SS_sm = vertcat(SS_sm',SS_c(end-floor(order/2):end));
%% Compare tt's for different distances
z_offset= -5;
% x_dist = 5000;
% x_dist = logspace(3,4.1,5);
%x_dist=(100:100:1000);
x_dist=(1000:2000:20000);
for i=1:length(x_dist)

[theta0_RAP(i),est_tt_RAP(i),SS_RAP,z_RAP] = ray_trace_RAP(z_c,SS_c,z_offset,-seafloor_z,x_dist(i),1);
[theta0_V(i),est_tt_V(i),SS_EC,z_EC]=ray_trace_Varamo_no_RAP(x_dist(i),z_offset,SS_c,z_c);
[est_tt_D(i),theta0_D(i)]=ray_trace_Dzieciuch2(z_RAP,SS_RAP,-z_offset,x_dist(i));
% pause
end

theta0_D = 90+theta0_D; %    Incident angle

%% PLOT
figure(1)
scatter(x_dist./1000,(est_tt_V -est_tt_D).*1000)
hold on
grid on
xlabel('Surface Distance (km)')
ylabel('Time Difference (ms)')
title('Varamo vs. Dzieciuch Ray Trace Estimation')

figure(2)
clf
plot(x_dist/1000,theta0_RAP/pi*180-theta0_D)
grid on
xlabel('Surface Distance (km)')
ylabel('Angle Difference in degree')
title('Varamo vs. Dzieciuch Ray Trace Estimation')

%% Function
function [theta0,est_tt,SS_flt,z_flt] = ray_trace_RAP(z_elps,SS_elps,td_h,ACO_depth,S,EFT)
ref_ellipsoid =  referenceEllipsoid('WGS84');
major_ax = ref_ellipsoid.SemimajorAxis;
minor_ax = ref_ellipsoid.SemiminorAxis;
ecc = ref_ellipsoid.Eccentricity;

if nargin < 11
    if ~exist('EFT','var')
        EFT = 1;
    end
end

switch EFT
    case 1
        EFT =1;
    case 'spherical'
        EFT = 2;
    case 'NoEFT'
        EFT = 3;
end

% remove bad data
if td_h >=0
    td_h = 0;
end

% Depths in the ellipsoidal coordinate; downward positive (-h)
% ACO_gh = geoidheight(22.738772,-158.006186+360,'EGM96');
ACO_gh = 0;
z_td = -td_h;
z_hyd = -(ACO_depth+ACO_gh);

% disp(sprintf('transducer depth = %f',z_td))
%% 2 truncate CTD data to fit the range between the receiver depth and the source depth
% on the ellipsoidal surface Orthometrc height
% N_mean = 2.32;                                      % average geiod height over the area of coverage
N_mean = 0;
z_elps = z_elps-N_mean;                              % from MSL to ellipsoid


% 2.1 truncate the upper bound using transducer depth
rm_ind = find(z_elps<=z_td);
z_elps(rm_ind(1:end-1)) = [];
SS_elps(rm_ind(1:end-1)) = [];
% set the first depth to be the transducer depth
z_diff = z_td - z_elps(1);
SS_elps(1) = SS_elps(1)+(SS_elps(2)-SS_elps(1))*(z_diff/(z_elps(2)-z_elps(1)));   % interpolate the first temp and salinity
z_elps(1) = z_td;

% 2.2 truncate the lower bound
if z_elps(end) > z_hyd
    rm_ind = find(z_elps>z_hyd);
    z_elps(rm_ind) = [];
    SS_elps(rm_ind) = [];
else
    ly_thk = -mean(z_elps(end-5:end-1) - z_elps(end-4:end));
    add_layer = floor((z_hyd -z_elps(end))/ly_thk);
    if add_layer ~= 0 
        z_elps = vertcat(z_elps,z_elps(end)+ly_thk*(1:add_layer+1)');
        SS_gradient = (SS_elps(end-1)-SS_elps(end))/ly_thk;
        SS_elps = vertcat(SS_elps,SS_elps(end)+SS_gradient*ly_thk*(1:add_layer+1)');
    end
   
end

%% Earth Flattening Transformation
% Chadwell 2010
% M_r = major_ax*(1-ecc^2)/(1-ecc^2*sin(tx_lat/180*pi)^2)^1.5;
% N_r = major_ax/sqrt(1-ecc^2*sin(tx_lat/180*pi)^2);
% alpha = alpha/180*pi;
% R_alpha = (cos(alpha)^2/M_r+sin(alpha)^2/N_r)^-1;
R_alpha = 6371000;
% hs = -z_elps;
% hf_AR = R_alpha.*log((R_alpha+hs)./R_alpha);
% hf_CS =  R_alpha.*log(R_alpha./(R_alpha-hs));





if EFT == 1
    % Chadwell 2010
    %     z_flt = -R_alpha.*log(R_alpha./(R_alpha-(-z_elps)));
    %     SS_flt = SS_elps.*(R_alpha./(R_alpha-(-z_elps)));
    % Dushaw and Colosi
    E_radius=6371000;
    E = z_elps'./E_radius;
    z_flt= z_elps'.*(1+(E./2)+((E.^2)./3));
%     z(end)=hyd_depth;
    SS_flt=SS_elps'.*(1+E+E.^2);
    
elseif EFT == 2
    z_flt = z_elps';
    SS_flt = SS_elps';
elseif EFT == 3
    z_flt = z_elps';
    SS_flt = SS_elps';
end

r = [];

%Find theta0 of acoustic transmission
theta0_all=0.00001:1*(pi/180):89.5/180*pi;   % all possible launch angles


R0 = R_alpha - z_td;           % Launch radial distance
R_hyd = R_alpha - z_hyd;      % Hydrophone radial distance
Ri = R_alpha - z_elps;

% loop over all possible thetas
% calculate a small horizontal distance within an interval betweenf two depths
for jj=1:length(theta0_all)
    if EFT == 2 % Spherical ray tracing
        % 1st incident angle at interface of layer 2
        incident_angle = asin(Ri(1)/Ri(2)*sin(theta0_all(jj)));
        
        %Ray parameter
        p = R0*sin(incident_angle)/SS_elps(1);
        
        
        % leaving angles of each layer
        theta = asin(p*SS_elps(2:end)./Ri(2:end));
        theta = vertcat(incident_angle,theta);
        
        
        % circumferential displacement in each layer
        for gg = 1:length(SS_elps)-1
            
            % APL with curavture correction factor
            fe = (R_alpha - z_elps(gg))/R_alpha;
            r(gg) = (z_elps(gg+1)-z_elps(gg))*tan(theta(gg))/fe;
            
        end
        
    else
        % Planar Ray Tracing
        %Ray parameter
        a= sin(theta0_all(jj))/SS_flt(1);
        
        % incident angles of each layer
        theta(1) = theta0_all(jj);
        
        for ii=2:length(SS_flt)-1
            theta(ii)=asin(a*SS_flt(ii));
        end
        
        % horizontl distance
        for ii=1:length(SS_flt)-1
            r(ii)=(z_flt(ii+1)-z_flt(ii))*tan(theta(ii));
        end
        
    end
    
    %Total x distance
    x_diff(jj)=sum(r);  % total horizontal distance for one possible launching theta
    
end


%Find min and max launch angle for Netwon-Raph method

t_low=find(x_diff<S);      % find all launch angles which yield x distance within the bound
theta_low=theta0_all(t_low(end));   % label the lowest possible launch angle (correspond to the steepest angle)
x_low=x_diff(t_low(end));           % The shortest ray arc length, lower bound of the sur dist

t_high=find(x_diff>S);


if abs(imag(x_low)) < 1e-12
    
    try
        theta_high=theta0_all(t_high(1));   % The upper bound of sur dist (the least steep angle)
        
        x_high=x_diff(t_high(1));
        
        % check if there is a real solution
        assert(isreal(x_high))
        
        % Intrapolation to find the optimum ray arc length of the direct path
        theta_new=theta_low+((S)/((x_low-x_high)/(theta_low-theta_high)))-((x_low)/((x_low-x_high)/(theta_low-theta_high)));
        
        % We now have a range of launch angles which gives the nearest number of
        % horizontal distance
        theta1=theta_low;
        x1=x_low;           % lower bound of the horizontal distance
        
        x_diff=1000000;
        
        % Loop until x_diff is close to x_dist to 0.1 mm
        while abs(x_diff-S)>0.0001
            
            r = [];
            tt = [];
            arc_distance = [];
            
            if EFT == 2 
                % Spherical Ray Tracing
                % 1st incident angle at interface of layer 2
                incident_angle = asin(Ri(1)/Ri(2)*sin(theta_new));
                
                %Ray parameter
                p = R0*sin(incident_angle)/SS_elps(1);
                
                % incident angles of each layer
                theta = asin(p*SS_elps(2:end)./Ri(2:end));
                theta = vertcat(incident_angle,theta);
                
                
                % circumferential displace ment in each layer
                % ray arc in each layer
                % travel time in each layer
                
                for gg = 1:length(Ri)-1
                    % APL
                    fe = (R_alpha - z_elps(gg))/R_alpha;
                    r(gg) = (z_elps(gg+1)-z_elps(gg))*tan(theta(gg))/fe;
                    tt(gg) = r(gg) *fe*csc(theta(gg))/SS_elps(gg);
                    arc_distance(gg) = tt(gg)*SS_elps(gg);  % approximate from travel time and sound speed
                    
                end
                
                
                
            else
                
                % Planar Ray Tracing
                %Ray parameter
                a=sin(theta_new)/SS_flt(1);
                
                % incident angle
                theta(1) = theta_new;
                for ii=2:length(SS_flt)-1
                    theta(ii)=asin(a*SS_flt(ii));
                end
                
                for ii=1:length(SS_flt)-1
                    r(ii)=(z_flt(ii+1)-z_flt(ii))*tan(theta(ii));
                    arc_distance(ii)=(z_flt(ii+1)-z_flt(ii))/cos(theta(ii));
                    tt(ii)=(z_flt(ii+1)-z_flt(ii))/(SS_flt(ii)*cos(theta(ii)));
                end
                
                
            end
            
            %Total x distance
            x_diff=sum(r);
            
            arc_lengths_hold = arc_distance;
            surface_dist_hold = r;
            ray_angles_hold = theta(2:end);
            
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
        theta0 = theta_new;                       % Launch angle
        tot_dist = sum(arc_lengths_hold);         % Ray arc length
        arc_lengths = arc_lengths_hold;           % arc length at each interval
        surface_dist = surface_dist_hold;         % Surface distance
        ray_angles = ray_angles_hold;             % ray angle at each interval
        t_layer = tt;                             % time in each interval
        
        % incremental horizontal distance
        inc_r = zeros(1,length(surface_dist)+1);
        
        for ii= 2:length(inc_r)
            inc_r(ii) = inc_r(ii-1)+surface_dist(ii-1);
        end
        
        
        
        est_tt=sum(tt);       %Estimated TT
        
        
    catch
        
        theta0= 0;                          % Launch angle
        tot_dist= 0 ;                       % Ray arc length
        arc_lengths = 0;                    % arc length at each interval
        surface_dist= 0;                    % Surface distance
        ray_angles= 0 ;                     % ray angle at each interval
        r = 0;
        est_tt= 0;                         %Estimated TT
        t_layer = 0;
    end
else
    theta0= 0;                          % Launch angle
    tot_dist= 0 ;                       % Ray arc length
    arc_lengths = 0;                    % arc length at each interval
    surface_dist= 0;                    % Surface distance
    ray_angles= 0 ;                     % ray angle at each interval
    r = 0;
    est_tt= 0;                         %Estimated TT
    t_layer = 0;
end

% disp(max(z_elps)-min(z_elps))


end

function SS_avg = mv_avg_SS(order,SS,ind)
    side = floor(order/2);
    handle_SS = 0;
    for r = 0:order-1
        handle_SS = handle_SS+SS(ind-side+r);
        
    end
    SS_avg = handle_SS/order;
       
end




