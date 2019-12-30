function [arc_lengths,tot_dist,theta0,SS_new,z_new,SS_HOT_avg,surface_dist,est_tt,ray_angles ] = ray_trace_w_curvature_v3(x_dist,z_offset,ACO_lat,ACO_lon,ACO_depth,month,year)
%READ_CTD_HOTS Summary of this function goes here
%%%%%%%%%NOTE: THEORETICAL MAX RANGE ALLOWABLE IS 29.7 km%%%%%%%%%%%%
%%%%%%%%%NOTE: THEORETICAL MAX RANGE ALLOWABLE IS 29.50 km%%%%%%%%%%%%
%%%%% parameters %%%%%%%%%%
% 1. ctd file
% 2. sound speed calculation function gsw_sound_speed
% 3. height calculation from presuure function gsw_z_from_p
% 4. enthalpy
%%%%% inputs %%%%%%%%%%%%
% 1. surface distances = x_dist (single value)
% 2. transducer altitudes = z_offset  ( single value)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%{
% HEMs        
% z_dist=4811.4;      %Hydrophone pressure at 4728m depth
% z_dist = 4812.7+.669;      %Paroscientific pressure sensor (dbar) + hyd depth (4729.92 m) 
% z_dist = 4819.897;     % at 4,736.266 m June 2017
%z_dist = 4822.107;   % at 4,738.415 m June 2018
% z_dist = 4818.461;  % at 4734.870 m March 2019 #2


% iclisten
% z_dist =  4818.0339   ;    % 4728+6.455 = 4,734.46 m icListen Sep 2018
% z_dist = 4819.4318;      % 4,735.71 m  March 2019 #2 0.556 shallower from 4,736.266
% z_dist = 4818.231;      % 4,734.646 m  March 2019 #3

%hyd_depth=4728+6.346;      % Hydrophone depth June 2017
%hyd_depth = 4736.494;      % Sep 2018
%}
%%%%%%%%%%%%

z_dist = gsw_p_from_z(ACO_depth,ACO_lat);

pres=D(:,1);        %pressure (dbars)
temp=D(:,2);        %temperature
sal=D(:,3);         %salinity
sal = gsw_SA_from_SP(sal,pres,ACO_lon,ACO_lat);

% Fix z offset for transducer location (m=dbar at this level)
% Find depths above the onboard transducer depth, and discard them and
% associated data points
pres_offset = gsw_p_from_z(z_offset,ACO_lat);
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

%Correct for Earth Curvature
% E_radius=6371000;        % epsilonray tracing APL code
E_radius=6370229;        

% coordinate transformation
E=z./E_radius;
z_new=z.*(1+(E./2)+((E.^2)./3));
SS_new=SS.*(1+E+E.^2);


% From Peter M Shearer
% r = Rc-z_elps;
% z_new = -Rc*log(r./Rc);
% SS_new = Rc./r.*SS_elps;

r = [];
%Find theta0 of acoustic transmission
% theta0_all=0.1:0.01*(pi/180):89*(pi/180);
% incident angle
theta0_all=0.0001:1*(pi/180):87.4/180*pi;   % all possible launch angles
theta0_all(end)=87.4/180*pi;

% loop over all possible thetas 
% calculate a small horizontal distance within an interval betweenf two depths
for jj=1:length(theta0_all)
    
    %Ray parameter
    a=sin(theta0_all(jj))/SS_new(1);
    
    % incident angles of each layer
    theta(1) = theta0_all(jj);
    
    for ii=2:length(SS_new)-1
        theta(ii)=asin(a*SS_new(ii));
    end    
    % SS gradient in each layer for the slowly-varying ss profile
%     for i=1:length(SS)-1
%         b(i)=(SS(i+1)-SS(i))/(z(i+1)-z(i));
%     end
    
    % horizontl distance
    for ii=1:length(SS_new)-1
        r(ii)=(z_new(ii+1)-z_new(ii))*tan(theta(ii));
    end
    
    %Total x distance
    x_diff(jj)=sum(r);  % total horizontal distance for one possible launching theta
    
    
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
    a=sin(theta_new)/SS_new(1);
    
    % incident angle
    theta(1) = theta_new;
    for ii=2:length(SS_new)-1
        theta(ii)=asin(a*SS_new(ii));
    end
    
    r = [];
    % horizontal range in each layer
    for ii=1:length(SS_new)-1
        r(ii)=(z_new(ii+1)-z_new(ii))*tan(theta(ii));
    end
    
    %Total x distance
    x_diff=sum(r);
    
    %Total distance and time traveled
    for ii=1:length(SS_new)-1
        arc_distance(ii)=(z_new(ii+1)-z_new(ii))/cos(theta(ii));
        tt(ii)=(z_new(ii+1)-z_new(ii))/(SS_new(ii)*cos(theta(ii)));
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

for ii= 2:length(inc_r)
    inc_r(ii) = inc_r(ii-1)+surface_dist(ii-1);
end

clear x_diff tot_dist_all arc_lengths_hold

est_tt=sum(tt);       %Estimated TT

%Vertically averaged SS for HOTS CTD
SS_HOT_avg=mean(SS_new);



% %Make variables same size
% if length(A)>1
%     arc_lengths=horzcat(zeros(1,length(A)-1),arc_lengths);
%     surface_dist=horzcat(zeros(1,length(A)-1),surface_dist);
%     ray_angles=horzcat(ones(1,length(A)-1).*ray_angles(1,1),ray_angles);
% end




end

