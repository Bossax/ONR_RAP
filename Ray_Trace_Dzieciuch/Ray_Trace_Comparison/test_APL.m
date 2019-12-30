
[est_tt,depth,theta0] = ray_tracing_planar_const_SS_APL(1000,-5)

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
plot(SS_flt,z_flt,'--')
hold on
plot(c_n,cumsum(dz),'--')
set(gca,'YDir','reverse')
end