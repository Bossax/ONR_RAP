function [ theta0,est_tt,SS,z ] = ray_trace_Varamo_no_RAP( x_dist,z_offset,SS,z )
%READ_CTD_HOTS Summary of this function goes here
%%%%%%%%%NOTE: THEORETICAL MAX RANGE ALLOWABLE IS 29.7 km%%%%%%%%%%%%

%Fix z offset for transducer location (m=dbar at this level)
A=find(-z_offset>=z);
z(A(1:end-1))=[];
z_diff=-z_offset-z(1);
z(1)=-z_offset;
SS(A(1:end-1))=[];
SS(1)=(SS(1)*(2-z_diff)+SS(2)*z_diff)/2;


%Correct for Earth Curvature
E_radius=6371000;
E=z./E_radius;
z=z.*(1+(E./2)+((E.^2)./3));
%z(end)=hyd_depth;
SS=SS.*(1+E+E.^2);

%Find theta0 of acoustic transmission
% theta0_all=0.1:0.01*(pi/180):89*(pi/180);
theta0_all=0.0001:1*(pi/180):1.51534;
theta0_all(end)=1.51534;
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
if length(A)>1
    arc_lengths=horzcat(zeros(1,length(A)-1),arc_lengths);
    surface_dist=horzcat(zeros(1,length(A)-1),surface_dist);
    ray_angles=horzcat(ones(1,length(A)-1).*ray_angles(1,1),ray_angles);
end

end

