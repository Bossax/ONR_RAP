function [ est_tt,theta0 ] = ray_trace_Dzieciuch2( z,SS,z_offset,x_dist )
%RAY_TRACE_DZIECIUCH Summary of this function goes here
%   Detailed explanation goes here


seafloor_z=z(end);
theta0_lower= 1;
theta0_upper=89.9;

zw=z;                          %DEPTH PROFILE (m)
cw=SS;                         %SS PROFILE (m/s)

cw=cw(:);
zw=zw(:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%source depth
srcdep=z_offset;                                 %SOURCE DEPTH (m)

% generate launch angles
% grazing angle
% upper negative
% lower positive
theta0=linspace(theta0_lower,theta0_upper,100);
% reciever range (x_dist)

% make a regularly-spaced set of ranges with interval rinc
rinc=5;     % 5 m
rgsv=linspace( 0, x_dist, ceil(x_dist/rinc)+1);
rgsv=rgsv(:);
clear rinc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ssr=[0:x_dist:x_dist];              % 0 range and max range
ssc=repmat( cw, 1, length(ssr));    % make sound speed profile a row vector

bthr = ssr;     % range
zb = seafloor_z;
bthz = zb*ones(size(bthr)); % depth

% polynomial structure
global zenv;
% calling function
% Cubic Spline Interpolation
zenv = zrayenv( zw, ssr, ssc, bthr, bthz);  % z_environemental input


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%raytrace accuracy
rtol = 1e-6;

% ray tracing function
% pass launch angles, source depth, incremental range vector, tolerable error
% upper-shot rays, find upper limit (surface bounces, limitingn upper turning rays
zry = zray2( theta0, srcdep, rgsv, rtol);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Find theta0 bottom limits%%%%%
% fprintf('Bottom Bounce\n')
% Find the index of the first Bottom Bounce of the first turning ray (in
% RAP case this would be the bottom bounce)
for i=1:length(zry)
    if isempty(zry(i).zltp) ~= 1
        break
    end
end

% reverse the sign of the theta0
theta_low = - theta0(i-1);
theta_high = -theta0(i);
z_low=zry(i-1).finalz;
z_high=seafloor_z+(seafloor_z-zry(i).finalz);


%Newton-Raphson Method
theta_new=theta_low+((seafloor_z)/((z_low-z_high)/(theta_low-theta_high)))-((z_low)/((z_low-z_high)/(theta_low-theta_high)));

theta1=theta_low;
z1=z_low;

z2=1000000;

disp(theta1)
while abs(z2-seafloor_z)>0.0001
    
    
zry=zray2( theta_new, srcdep, rgsv, rtol);

% check ray turning 
if isempty(zry.zltp)==1
    z2=zry.finalz;
elseif isempty(zry.zltp)~=1
    z2=seafloor_z+(seafloor_z-zry.finalz);
end

theta2=theta_new;

theta_new=theta2+((seafloor_z)/((z2-z1)/(theta2-theta1)))-((z2)/((z2-z1)/(theta2-theta1)));

theta1=theta2;
z1=z2;
fprintf('depth Diff = %f\n',(z2-seafloor_z))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
est_tt=zry.finalt;

theta0=theta_new;



end

