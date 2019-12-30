function [arc_lengths,tot_dist,theta0,SS_flt,z_flt,SS_HOT_avg,surface_dist,r,est_tt,t_layer,ray_angles,R_alpha] = ray_trace_w_earth_flattening(S,td_h,tx_lon,tx_lat,alpha,ACO_lat,ACO_lon,ACO_depth,month,year,EFT)
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
% Chadwell 2010
M_r = major_ax*(1-ecc^2)/(1-ecc^2*sin(tx_lat/180*pi)^2)^1.5;
N_r = major_ax/sqrt(1-ecc^2*sin(tx_lat/180*pi)^2);
alpha = alpha/180*pi;
R_alpha = (cos(alpha)^2/M_r+sin(alpha)^2/N_r)^-1;

% hs = -z_elps;
% hf_AR = R_alpha.*log((R_alpha+hs)./R_alpha);
% hf_CS =  R_alpha.*log(R_alpha./(R_alpha-hs));





if EFT == 1
    % Chadwell 2010
%     z_flt = -R_alpha.*log(R_alpha./(R_alpha-(-z_elps)));
%     SS_flt = SS_elps.*(R_alpha./(R_alpha-(-z_elps)));
% Dushaw and Colosi
    E = z_elps./R_alpha;
    z_flt = z_elps.*(1+(E./2)+((E.^2)./3));
    SS_flt = SS_elps.*(1+E+E.^2);
    
elseif EFT == 2
    z_flt = z_elps;
    SS_flt = SS_elps;
elseif EFT == 3
    z_flt = z_elps;
    SS_flt = SS_elps;
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
    if EFT == 2
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
            %Ray parameter
            a=sin(theta_new)/SS_flt(1);
            
            % incident angle
            theta(1) = theta_new;
            for ii=2:length(SS_flt)-1
                theta(ii)=asin(a*SS_flt(ii));
            end
            
            
            % horizontal range in each layer
            for ii=1:length(SS_flt)-1
                r(ii)=(z_flt(ii+1)-z_flt(ii))*tan(theta(ii));
            end
            
            %Total distance and time traveled
            for ii=1:length(SS_flt)-1
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


    %Vertically averaged SS for HOTS CTD
    SS_HOT_avg=mean(SS_flt);


end

