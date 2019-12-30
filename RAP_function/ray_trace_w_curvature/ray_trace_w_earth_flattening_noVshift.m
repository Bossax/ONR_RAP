function [arc_lengths,tot_dist,theta0,SS_flt,z_flt,SS_HOT_avg,surface_dist,r,est_tt,ray_angles,R_alpha] = ray_trace_w_earth_flattening_no_Vshift(S,td_h,tx_lon,tx_lat,alpha,ACO_lat,ACO_lon,ACO_depth,month,year)
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
% Ray tracing using esitmated spherical earth model
% z positive downward (d)
% ellipsoidal heights positive upward
% Reference Ellipsoid
ref_ellipsoid =  referenceEllipsoid('WGS84');
major_ax = ref_ellipsoid.SemimajorAxis;
minor_ax = ref_ellipsoid.SemiminorAxis;
ecc = ref_ellipsoid.Eccentricity;

% remove bad data
if td_h >=0
    td_h = 0;
end

% shift transducer depth to be relative to geoid
% td_h = td_h - geoidheight(tx_lat,tx_lon+360,'EGM96');
td_h = td_h - 2.31;

% disp(sprintf('Geoid Height = %f',geoidheight(tx_lat,tx_lon+360,'EGM96')))
%% 1. CTD data (MSL)
% CTD file directory
% cd '/Users/testuser/Documents/MATLAB/Script/Data'
if year == '2017'
    fid=fopen('h294a0202.ctd');                       % June 2017
elseif year == '2018'
    switch month
        case 'Jun'
            fid = fopen('h302a0201.ctd');                    % June 2018
        case 'Sep'
            fid = fopen('h305a0202.ctd');                       % 10 Sep 2018
        case 'Oct'
            fid=fopen('h306a0202.ctd');                       % 12 October 2018
        case 'Nov'
            fid=fopen('h307a0202.ctd');                       % 16 November 2018
    end
end

D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f','headerlines',6));
fclose(fid);

pres_MSL = D(:,1);        %pressure (dbars)
temp = D(:,2);        %temperature
sal = D(:,3);         %salinity (Sp)
sal = gsw_SA_from_SP(sal,pres_MSL,ACO_lon,ACO_lat);

% Depths relative to geoid
z_td = -td_h;
z_hyd = -ACO_depth;

% disp(sprintf('transducer depth = %f',z_td))
%% 2 truncate CTD data to fit the range between the receiver depth and the source depth

z_CTD = -(gsw_z_from_p(pres_MSL,ACO_lat) );         % relative to MSL

% 2.1 truncate the upper bound using transducer depth
rm_ind = find(z_CTD<=z_td); 
z_CTD(rm_ind(1:end-1)) = [];
temp(rm_ind(1:end-1)) = [];
sal(rm_ind(1:end-1)) = [];
pres_MSL(rm_ind(1:end-1)) = [];

% set the first depth to be the transducer depth
z_diff = z_td - z_CTD(1);
temp(1) = temp(1)+(temp(2)-temp(1))*(z_diff/(z_CTD(2)-z_CTD(1)));   % interpolate the first temp and salinity
sal(1)=sal(1)+(sal(2)-sal(1))*(z_diff/(z_CTD(2)-z_CTD(1))); 
pres_MSL(1) = pres_MSL(1)+(pres_MSL(2)-pres_MSL(1))*(z_diff/(z_CTD(2)-z_CTD(1))); 
z_CTD(1) = z_td;

% 2.2 truncate the lower bound
if z_CTD(end) > z_hyd
    rm_ind = find(z_CTD>z_hyd);
    z_CTD(rm_ind) = [];
    pres_MSL(rm_ind) = [];
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
    pres_MSL_gradient = pres_MSL(end-1)-pres_MSL(end);
    pres_MSL = vertcat(pres_MSL,pres_MSL(end)+pres_MSL_gradient*(1:add_layer+1)');
end

%% 3. Creat Sound Speed Profile 
SS_elps = gsw_sound_speed(sal,temp,pres_MSL);


%% Earth Flattening Transformation
% Chadwell 2010
M_r = major_ax*(1-ecc^2)/(1-ecc^2*sin(tx_lat/180*pi)^2)^1.5;
N_r = major_ax/sqrt(1-ecc^2*sin(tx_lat/180*pi)^2);
alpha = alpha/180*pi;
R_alpha = (cos(alpha)^2/M_r+sin(alpha)^2/N_r)^-1;

hs = -z_CTD;
hf_AR = R_alpha.*log((R_alpha+hs)./R_alpha);
hf_CS =  R_alpha.*log(R_alpha./(R_alpha-hs));

z_flt = R_alpha.*log(R_alpha./(R_alpha-z_CTD));    
SS_flt= SS_elps.*(R_alpha./(R_alpha-(z_CTD)));

% z_flt = -R_alpha.*log(R_alpha./(R_alpha-(-z_CTD)));    
% SS_flt = SS_elps.*(R_alpha./(R_alpha-(-z_CTD)));

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

