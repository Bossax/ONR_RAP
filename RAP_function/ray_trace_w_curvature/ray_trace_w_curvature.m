function [arc_lengths,tot_dist,theta0,SS,z,SS_HOT_avg,surface_dist,est_tt,ray_angles ] = ray_trace_w_curvature(x_dist,z_offset )
%READ_CTD_HOTS Summary of this function goes here
%%%%%%%%%NOTE: THEORETICAL MAX RANGE ALLOWABLE IS 29.50 km%%%%%%%%%%%%
%%%%% parameters %%%%%%%%%%
% 1. ctd file
% 2. sound speed calculation function gsw_sound_speed
% 3. height calculation from presuure function gsw_z_from_p
% 4. enthalpy
%%%%% inputs %%%%%%%%%%%%
% 1. surface distances = x_dist (single value)
% 2. ship altitudes = z_offset  ( single value)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cd '/Volumes/LaCie_RAP_Backup/RAP/June2018Cruise/Tx_Rx_Output/script/Data'

%fid=fopen('h269a0215.ctd');       %February 2015
%fid=fopen('h281a0202.ctd');       %February 2016

%fid=fopen('h289a0201.ctd');       %January 2017
fid = fopen('h294a0202.ctd');      %June 2017
fid = fopen('h306a0202.ctd'); % October

D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f','headerlines',6));
fclose(fid);

% Grouping CTD parameters
% z_dist=4728;        %Hydrophone depth
z_dist=4812.4;        %Hydrophone pressure at 4728m depth 
z_dist = 4819.897;     % at 4,736.266 m June 2017

pres=D(:,1);         %pressure (dbars)
temp=D(:,2);         %temperature
sal=D(:,3);          %salinity
cd ../ray_trace_w_curvature
%% Fix z offset for transducer location (m=dbar at this level)
% Find depths above the onboard transducer depth, and discard them and
% associated data points
A=find(z_offset>pres);
pres(A(1:end-1))=[];
pres(1)=z_offset;       % set the first depth to be the transducer depth
temp(A(1:end-1))=[];
temp(1)=(temp(1)+temp(2))/2;    % interpolate the first temp and salinity
sal(A(1:end-1))=[];
sal(1)=(sal(1)+sal(2))/2;

%% Remove/Add deeper depths
if pres(end)>z_dist
remove_points=find(pres>z_dist);
pres(remove_points)=[];
temp(remove_points)=[];
sal(remove_points)=[];
else
    pres(end+1)=z_dist;
    temp(end+1)=temp(end);
    sal(end+1)=sal(end)+(sal(end)-sal(end-1));
end

%% A priori SS field
% get the priori sound speed field as a function of depth
SS= gsw_sound_speed(sal,temp,pres);

z=-gsw_z_from_p(pres,22.738894);

%% Correct for Earth Curvature
E_radius=6371000;
E=z./E_radius;                      % dummy variable
z=z.*(1+(E./2)+((E.^2)./3));        % depth
SS=SS.*(1+E+E.^2);                  % sound spped

%% Find theta0 of acoustic transmission
% theta0_all=0.1:0.01*(pi/180):89*(pi/180);
theta0_all=0.01:1*(pi/180):1.474139396298460;   % 0 to  170 degree
theta0_all(end)=1.474139396298460;

% loop over all possible thetas 
% calculate a small distance between an interval of two depths
for ii=1:length(theta0_all)
    
    %Ray calculations
    a=sin(theta0_all(ii))/SS(1);
    
    for i=1:length(SS)
        theta(i)=asin(a*SS(i));
    end    
    for i=1:length(SS)-1
        b(i)=(SS(i+1)-SS(i))/(z(i+1)-z(i));
    end
    for i=1:length(b)
        r(i)=(1/(a*b(i)))*(cos(theta(i))-cos(theta(i+1)));
    end
    
    %Total x distance
    x_diff(ii)=sum(r);      % total distance for one possible launching theta
    
%     %Total distance traveled
%     for i=1:length(b)
%         arc_distance(i)=(1/(a*b(i)))*(theta(i+1)-theta(i));
%     end
%     tot_dist_all(ii)=sum(arc_distance);
%     
%     arc_lengths_hold(:,ii)=arc_distance;
    
end


%% Filtering the result

t_low = find(x_diff<x_dist);          % find the ray arc lengths shorter than the surface distance
theta_low = theta0_all(t_low(end));   % label the lowest possible launch angle (correspond to the greatest distance close to the surface distance )
x_low = x_diff(t_low(end));           % The shortest ray arc length, lower bound of the sur dist

t_high = find(x_diff>x_dist);         % The upper bound of sur dist
theta_high=theta0_all(t_high(1));
x_high=x_diff(t_high(1));

%% Intrapolation to find the optimum ray arc length of the direct path
theta_new=theta_low+((x_dist)/((x_low-x_high)/(theta_low-theta_high)))-((x_low)/((x_low-x_high)/(theta_low-theta_high)));

theta1=theta_low;
x1=x_low;

x_diff=1000000;
% narrowen down the range of interest to between theta_low to theta_new
while abs(x_diff-x_dist)>0.0001
    
    %Ray calculations
    a=sin(theta_new)/SS(1);
    
    for i=1:length(SS)
        theta(i)=asin(a*SS(i));
    end    
    for i=1:length(SS)-1
        b(i)=(SS(i+1)-SS(i))/(z(i+1)-z(i));
    end
    for i=1:length(b)
        r(i)=(1/(a*b(i)))*(cos(theta(i))-cos(theta(i+1)));
    end
    
    %Total x distance
    x_diff=sum(r);
    
    %Total distance traveled
    for i=1:length(b)
        arc_distance(i)=(1/(a*b(i)))*(theta(i+1)-theta(i));
    end
    tot_dist_all=sum(arc_distance);
    
    arc_lengths_hold=arc_distance;
    surface_dist_hold=r;
    ray_angles_hold=theta(2:end);
    
    theta2=theta_new;
    x2=x_diff;
    
    % Udate theta_new by interpolation (Newton Rhapson)
    theta_new=theta2+((x_dist)/((x2-x1)/(theta2-theta1)))-((x2)/((x2-x1)/(theta2-theta1)));
    
    % Set theta_low to theta 2
    theta1=theta2;
    
    % update ray arc length
    x1=x2;
    
end


%% packaging outputs
theta0=theta_new;                           % Launch angle
tot_dist=tot_dist_all;                      % ray arc length of the direct path

arc_lengths=arc_lengths_hold;               % segments of ray arc length
surface_dist=surface_dist_hold;             % segments of surface distance
ray_angles=ray_angles_hold;                 % segments of ray angles at each depth interval

clear x_diff tot_dist_all arc_lengths_hold


est_tt=sum(arc_distance'./SS(2:end));

%Vertically averaged SS for HOTS CTD
SS_HOT_avg=mean(SS);

end

