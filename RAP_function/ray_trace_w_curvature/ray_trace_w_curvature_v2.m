function [ arc_lengths,tot_dist,theta0,SS,z,SS_HOT_avg,surface_dist,est_tt,ray_angles ] = ray_trace_w_curvature_v2( x_dist,z_offset )
%READ_CTD_HOTS Summary of this function goes here
%%%%%%%%%NOTE: THEORETICAL MAX RANGE ALLOWABLE IS 29.7 km%%%%%%%%%%%%

% fid=fopen('h296a0214.ctd');       %October 2017
fid=fopen('h306a0202.ctd');                       % October 2018

D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f','headerlines',6));
fclose(fid);
        
%z_dist=4811.4;      %Hydrophone pressure at 4728m depth
z_dist=4812.7+.669;      %Paroscientific pressure sensor (dbar) + hyd depth 
z_dist = 4819.897;     % at 4,736.266 m June 2017
z_dist = 4818.231;      % 4,734.646 m  March 2019 #3
% hyd_depth=4728;     %Hydrophone depth
pres=D(:,1);        %pressure (dbars)
temp=D(:,2);        %temperature
sal=D(:,3);         %salinity

%Fix z offset for transducer location (m=dbar at this level)
A=find(-z_offset>=pres);
pres(A(1:end-1))=[];
pres_diff=-z_offset-pres(1);
pres(1)=-z_offset;
temp(A(1:end-1))=[];
temp(1)=(temp(1)*(2-pres_diff)+temp(2)*pres_diff)/2;
sal(A(1:end-1))=[];
sal(1)=(sal(1)*(2-pres_diff)+sal(2)*pres_diff)/2;

%Remove/Add deeper depths
if pres(end)>z_dist
remove_points=find(pres>z_dist);
pres(remove_points)=[];
temp(remove_points)=[];
sal(remove_points)=[];
else
    pres(end+1)=z_dist;
    temp(end+1)=temp(end);
    sal(end+1)=sal(end);
    %sal(end+1)=sal(end)+(sal(end)-sal(end-1));
end

%Apriori SS field
SS=gsw_sound_speed(sal,temp,pres);

z=-gsw_z_from_p(pres,22.738894);

%Correct for Earth Curvature
E_radius=6371000;

E=z./E_radius;
z=z.*(1+(E./2)+((E.^2)./3));
%z(end)=hyd_depth;
SS=SS.*(1+E+E.^2);

%Find theta0 of acoustic transmission
theta0_all=0.1:0.01*(pi/180):89*(pi/180);
% theta0_all=0.0001:1*(pi/180):1.51534;
% theta0_all(end)=1.51534;
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
        r(i)=(z(i+1)-z(i))*tan(theta(i));
    end
    
    %Total x distance
    x_diff(ii)=sum(r);
    
    
end



%Find min and max launch angle for Netwon-Raph method 
t_low=find(x_diff<x_dist);
theta_low=theta0_all(t_low(end));
x_low=x_diff(t_low(end));
t_high=find(x_diff>x_dist);
theta_high=theta0_all(t_high(1));
x_high=x_diff(t_high(1));

theta_new=theta_low+((x_dist)/((x_low-x_high)/(theta_low-theta_high)))-((x_low)/((x_low-x_high)/(theta_low-theta_high)));

theta1=theta_low;
x1=x_low;

x_diff=1000000;
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
        r(i)=(z(i+1)-z(i))*tan(theta(i));
    end
    
    %Total x distance
    x_diff=sum(r);
    
    %Total distance and time traveled
    for i=1:length(b)
        arc_distance(i)=(z(i+1)-z(i))/cos(theta(i));
        tt(i)=(z(i+1)-z(i))/(SS(i)*cos(theta(i)));
    end
    tot_dist_all=sum(arc_distance);
    
    arc_lengths_hold=arc_distance;
    surface_dist_hold=r;
    ray_angles_hold=theta(2:end);
    
    %Newton-Raphson method
    theta2=theta_new;
    x2=x_diff;
    
    theta_new=theta2+((x_dist)/((x2-x1)/(theta2-theta1)))-((x2)/((x2-x1)/(theta2-theta1)));
    
    theta1=theta2;
    x1=x2;
    
end



theta0=theta_new;               %Launch angle
tot_dist=tot_dist_all;          %Ray arc length

arc_lengths=arc_lengths_hold;           %arc length at each interval
surface_dist=surface_dist_hold;         %Surface distance
ray_angles=ray_angles_hold;             %ray angle at each interval

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

