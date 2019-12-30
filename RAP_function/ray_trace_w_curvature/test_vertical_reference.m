% Test effects of the vertical reference 
% test cases
% 1. 
clear
close all
%% Hydrophone Position
N_mean = 2.32;
HEM_lat = 22.738772;
HEM_lon = -158.006186;
HEM_depth = -4729.92;                   % original depth MSL local u/e/n frame

%% compute estimated travel times of 2 cases
% 1. all relative to MSL (geoid)
% 2. all relative to ellipsoid
est_tt_MSL = [];            
est_tt_elps = [];
est_tt_cs = [];
detph_inv_MSL = [];
detph_inv_elps = [];
detph_inv_cs = [];

ray_arc_length_MSL = [];
ray_arc_length_elps= [];
ray_arc_length_cs = [];

med_SS_diff_MSL_elps = [];
med_SS_diff_CS_DC_EFT = [];

% SSP for each surface range
% geodesics
geo_S = [10 50 100:100:400 500:100:27000];

% fixed transducer depth
td_h = -4 ;  % ellipsoidal height

for ii =1:length(geo_S)
   
    % MSL 
    [SS_MSL,z_MSL,SS_MSL_flt,z_MSL_flt,est_tt_MSL(ii),ray_arc_length_MSL(ii)] = ray_tracing(geo_S(ii),td_h,HEM_lat,HEM_lon,HEM_depth,'MSL','MSL',2);
    
    % Ellipsoid
    [SS_elps,z_elps,SS_elps_flt,z_elps_flt,est_tt_elps(ii),ray_arc_length_elps(ii)] = ray_tracing(geo_S(ii),td_h,HEM_lat,HEM_lon,HEM_depth,'Ellipsoid','Ellipsoid',2); 
    
    % Chadwell Sweeney EFT
    [SS_cs,z_cs,SS_cs_flt,z_cs_flt,est_tt_cs(ii),ray_arc_length_cs(ii)] = ray_tracing(geo_S(ii),td_h,HEM_lat,HEM_lon,HEM_depth,'Ellipsoid','Ellipsoid',3); 
   
   med_SS_diff_MSL_elps(ii) =  median(SS_MSL_flt - SS_elps_flt);
   med_SS_diff_CS_DC_EFT(ii) = median(SS_cs_flt - SS_elps_flt);
   detph_inv_MSL(ii) = range(z_MSL_flt);
   detph_inv_elps(ii) = range(z_elps_flt);
   detph_inv_cs(ii) = range(z_cs_flt);
end

%% Plot
%% tt different
figure(1)
scatter(geo_S/1000,(est_tt_elps-est_tt_MSL)*1000,'.')
grid on
xlabel('Range (km)')
ylabel('msec')
% title({'Estimated Travel Time Difference','(Chadwell and Sweeney EFT - Aki and Richards EFT)'})
title({'Estimated Travel Time Difference','(vertical reference at ellipsoidal surface - MSL)'})
set(gca,'fontsize',14)
% ylim([-0.014 -0.04])

%% ray arc length difference
figure(2)
scatter(geo_S/1000,(ray_arc_length_MSL-ray_arc_length_elps)*100,'.')
ylim([-2 0])
ylabel('Arc Length Difference (cm)')
xlabel('Range (km)')
title('Ray Arc Length Difference (MSL - Ellipsoid/ APL EFT)')
grid on
set(gca,'fontsize',14)
%% Functions
function [SS,z_adjusted,SS_flt,z_flt,est_tt,tot_dist] = ray_tracing(S,td_h,ACO_lat,ACO_lon,ACO_depth,RT_ref,SSP_ref,EFT)
% RT_ref = 'MSL', 'Ellipsoid'
% Reference Ellipsoid
ref_ellipsoid =  referenceEllipsoid('WGS84');
major_ax = ref_ellipsoid.SemimajorAxis;
minor_ax = ref_ellipsoid.SemiminorAxis;
ecc = ref_ellipsoid.Eccentricity;
N_mean = 2.32;                                      % average geiod height over the area of coverage 
%% 1. CTD data (MSL)
% CTD file directory
fid=fopen('h306a0202.ctd');                       % 12 October 2018
       
D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f','headerlines',6));
fclose(fid);

pres_MSL = D(:,1);        %pressure (dbars)
temp = D(:,2);        %temperature
sal = D(:,3);         %salinity (Sp)
sal = gsw_SA_from_SP(sal,pres_MSL,ACO_lon,ACO_lat);


% Depths of the transducer and the hydrophone
ACO_gh = geoidheight(ACO_lat,ACO_lon+360,'EGM96');
switch RT_ref 
    case'Ellipsoid'
        z_td = -td_h;
        z_hyd = -(ACO_depth+ACO_gh);    % h = (H+N) dF = -(H+N)
    case 'MSL'
        z_td = -(td_h-N_mean);              % d_MSL = -(h-N)
        z_hyd = -ACO_depth;
end


%% 2 truncate CTD data to fit the range between the receiver depth and the source depth
z_CTD = -(gsw_z_from_p(pres_MSL,ACO_lat));         % relative to MSL
switch SSP_ref
    case 'Ellipsoid'
        z_adjusted = z_CTD-N_mean;                              % from MSL to ellipsoid
    case 'MSL'
        z_adjusted = z_CTD;                                     % remain at MSL
end

% 2.1 truncate the upper bound using transducer depth
rm_ind = find(z_adjusted<=z_td);
z_adjusted(rm_ind(1:end-1)) = [];
z_CTD(rm_ind(1:end-1)) = [];
temp(rm_ind(1:end-1)) = [];
sal(rm_ind(1:end-1)) = [];

% set the first depth to be the transducer depth
z_diff = z_td - z_adjusted(1);
temp(1) = temp(1)+(temp(2)-temp(1))*(z_diff/(z_adjusted(2)-z_adjusted(1)));   % interpolate the first temp and salinity
sal(1)=sal(1)+(sal(2)-sal(1))*(z_diff/(z_adjusted(2)-z_adjusted(1)));
z_adjusted(1) = z_td;
switch SSP_ref
    case 'Ellipsoid'
        z_CTD(1) = z_adjusted(1)+N_mean;    % transducer depth + mean geoid undulation move back to MSL
    case 'MSL'
        z_CTD(1) = z_adjusted(1);            % MSL
end



% 2.2 truncate the lower bound
if z_adjusted(end) > z_hyd
    rm_ind = find(z_adjusted>z_hyd);
    z_adjusted(rm_ind) = [];
    z_CTD(rm_ind) = [];
    temp(rm_ind) = [];
    sal(rm_ind) = [];
else
    ly_thk = -mean(z_adjusted(end-5:end-1) - z_adjusted(end-4:end));
    add_layer = floor((z_hyd -z_adjusted(end))/ly_thk);
    z_adjusted = vertcat(z_adjusted,z_adjusted(end)+ly_thk*(1:add_layer)');
    z_adjusted(end+1) = z_hyd;
    
    
    switch SSP_ref
    case 'Ellipsoid'
       z_CTD = vertcat(z_CTD,z_adjusted(end-add_layer:end)+N_mean);% transducer depth + mean geoid undulation move back to MSL
    case 'MSL'
       z_CTD = vertcat(z_CTD,z_adjusted(end-add_layer:end));          % MSL
    end
    
    temp_gradient = temp(end-1)-temp(end);
    temp = vertcat(temp,temp(end)+temp_gradient*(1:add_layer+1)');
    sal_gradient = sal(end-1)-sal(end);
    sal = vertcat(sal,sal(end)+sal_gradient*(1:add_layer+1)');
end

pres_CTD = gsw_p_from_z(-z_CTD,ACO_lat);
%% 3. Creat Sound Speed Profile 
SS = gsw_sound_speed(sal,temp,pres_CTD);


%% Earth Flattening Transformation
% formula 1 Dushaw and Colosi
E_radius=6370229;     

if EFT == 1
    
    E = z_adjusted./E_radius;
    z_flt = z_adjusted.*(1+(E./2)+((E.^2)./3));
    SS_flt = SS.*(1+E+E.^2);

% formula 2 Aki and Richard
elseif EFT == 2
    z_flt = E_radius.*log(E_radius./(E_radius-z_adjusted));
    SS_flt= SS.*(E_radius./(E_radius-(z_adjusted)));

elseif EFT == 3
    z_flt = -E_radius.*log(E_radius./(E_radius-(-z_adjusted)));    
    SS_flt = SS.*(E_radius./(E_radius-(-z_adjusted)));

end
r = [];

%Find theta0 of acoustic transmission
theta0_all=0.00001:1*(pi/180):89.5/180*pi;   % all possible launch angles



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
    
    
catch

    theta0= 0;                          % Launch angle
    tot_dist= 0 ;                       % Ray arc length
    arc_lengths = 0;                    % arc length at each interval
    surface_dist= 0;                    % Surface distance
    ray_angles= 0 ;                     % ray angle at each interval
    r = 0;
    est_tt= 0;                         %Estimated TT
    
end
else
    theta0= 0;                          % Launch angle
    tot_dist= 0 ;                       % Ray arc length
    arc_lengths = 0;                    % arc length at each interval
    surface_dist= 0;                    % Surface distance
    ray_angles= 0 ;                     % ray angle at each interval
    r = 0;
    est_tt= 0;                         %Estimated TT
    
end


    %Vertically averaged SS for HOTS CTD
    SS_HOT_avg=mean(SS_flt);


end


