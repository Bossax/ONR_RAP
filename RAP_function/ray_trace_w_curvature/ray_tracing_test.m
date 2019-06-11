% ray_tracing_const_SS_gradient
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
clear
close all
x_dist = 1680;    % meter
z_offset = -6.2; % meter
cd '/Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Script/Data'
%fid=fopen('h294a0202.ctd');                       % June 2018
fid=fopen('h306a0202.ctd');                       % October 2018

D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f','headerlines',6));
fclose(fid);
        
%z_dist=4812.7+.669;      %Paroscientific pressure sensor (dbar) + hyd depth 4729.92 m 
z_dist = 4819.897;     % at 4,736.266 m June 2017

pres=D(:,1);        %pressure (dbars)
temp=D(:,2);        %temperature
sal=D(:,3);         %salinity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sound speed profile preparation

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

%Correct for Earth Curvature
E_radius=6371000;        % epsilonray tracing APL code
% coordinate transformation
E=z./E_radius;
z=z.*(1+(E./2)+((E.^2)./3));
%z(end)=hyd_depth;
SS=SS.*(1+E+E.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants
theta0= 60;             % degree
theta0=theta0/180*pi;   % radian      
layer1 = 1;
layer2 = length(SS)-1;             % max layers length(SS)-1


    %Ray parameter
    a=sin(theta0)/SS(1);
    
    % incident angles of each layer
    theta(1) = theta0;
    for i=2:length(SS)
        theta(i)=asin(a*SS(i));
    end    
    
%% Constant SS calculation  
inc_r = 0;
r = [];
tt = [];
arc_distance = [];
    % horizontl distance
    for i=layer1:layer2
        r(i)=(z(i+1)-z(i))*tan(theta(i));
        inc_r(end+1) = inc_r(end)+r(i);
    end
    
    % arc distance and travel time
    for i=layer1:layer2
        arc_distance(i)=(z(i+1)-z(i))/cos(theta(i));
        tt(i)=(z(i+1)-z(i))/(SS(i)*cos(theta(i)));
    end
    
    %Total r distance and 
    x_range=sum(r);  % total horizontal distance for one possible launching theta
    tot_acr_length = sum(arc_distance);
    est_tt = sum(tt);

f1 = figure(1);
f1.Units = 'normalized';
f1.Position = [0.7 0.5 0.3 0.6];
plot(inc_r,-z(layer1:layer2+1))
grid on
ylim([-z(layer2+1) 0])
xlabel('Horizontal Distance (m)')
ylabel('Depth')
title('Constant SS Model')
distances = sprintf('Distance from the hydrophone = %.4f m\nTotal Arc Length = %.4f m \nTotal Horizontal Distance = %.6f m \nTravel Time = %.6f s\n\n.......\n',x_dist,tot_acr_length,x_range,est_tt);
fprintf(distances)
%% Constant SS gradient calculation  
inc_r2 = 0;
r2 = [];
t = [];
% SS gradient in each layer for the slowly-varying ss profile
    for i=layer1:layer2
        b(i)=(SS(i+1)-SS(i))/(z(i+1)-z(i));
    end

% horizontl distance, arc length and travel time
    for j=layer1:layer2
        v = SS(j+1)/SS(j);
        h = (sqrt(1-(a*SS(j+1))^2)+1)/(sqrt(1-(a*SS(j))^2)+1);
       t(j) = (1/b(j))*(log(v)-log(h));
       r2(j) = (1/(a*b(j)))*(sqrt(1-(a*SS(j))^2) - sqrt(1-(a*SS(j+1))^2));
       arc_distance2(j)= 1/abs(a*b(j))*abs(theta(j+1) - theta(j));
       inc_r2(end+1) = inc_r2(end)+r2(j); 
    end   
   
 
    %Total r distance and 
    x_range2=sum(r2);  % total horizontal distance for one possible launching theta
    tot_acr_length2 = sum(arc_distance2);
    est_tt2 = sum(t);

f2 = figure(2);
f2.Units = 'normalized';
f2.Position = [0.4 0.5 0.3 0.6];
plot(inc_r2,-z(layer1:layer2+1))
grid on
ylim([-z(layer2+1) 0])
distances = sprintf('Distance from the hydrophone = %.4f m\nTotal Arc Length = %.4f m \nTotal Horizontal Distance = %.6f m \nTravel Time = %.6f s\n\n.......\n',x_dist,tot_acr_length2,x_range2,est_tt2);
fprintf(distances)
xlabel('Horizontal Distance (m)')
ylabel('Depth')
title('Constant Gradient SS Model')

f3 = figure(3);
f3.Units = 'normalized';
f3.Position = [0 0.5 0.3 0.6];
plot(theta/pi*180,-z)
grid on