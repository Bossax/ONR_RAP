% use ray equations from Computationl Ocean Acoustics Textbook
% equaitons are expressed in Ray coordinates
% use both depth and horizontal range to determine eigenrays
% arc length is the differential variable

clear
close all
 
figure(8)
global c_v z_v
c_v = [];
z_v = [];

geodesic = 20000;   % test 20 km
hyd_depth = 4730;
%% load sound speed profile
SS = soundspeedprofile(hyd_depth);

slant_angle = acot(geodesic/hyd_depth)/pi*180;
upper_ang = slant_angle+10;
if upper_ang >=89.9 % min range ~ 9m
    upper_ang = 89.9;
end

lower_ang = slant_angle-15;
if lower_ang<= 1
    lower_ang = 1;
end
theta0= upper_ang:-0.5:lower_ang;   % all possible launch angles

sinc = .5;     % arc length increment 1m 

r0 = 0; % source r0
z0 = 5; % source depth z0

rcv_depth = SS.z(end);
% store depth and range of each launch angle
Range = [];
Depth = [];
break_trig = 0;
for n = 1:length(theta0)
    n
% initial state parameters  
  idt = 1; idz = 2; idr = 3; idsh = 4; idsv = 5;
    [c0,~] = find_soundspeed(z0,SS);
    X0(idt) = 0;                             % time      
    X0(idz)= z0;                            % depth
    X0(idr)= r0;                            % range
    X0(idsh)= cos(pi*theta0(n)/180)./c0;    % horizontal slowness
    X0(idsv)= sin(pi*theta0(n)/180)./c0;    % vertical slowness
   
    % initial value
    snow = 0; 
    Xnow = X0;
    znow = z0;
    rnow = r0;
    X = Xnow;
    % loop until the horizonal range exceed rend
    while (rnow < geodesic)
        % calculate the changes of parameters within this arc length interval
        f = diff_X(snow,Xnow,SS);
    
        % update the state vector 
        Xnow = Xnow + sinc*f;
        snow = snow+sinc;
        znow = Xnow(end,idz);
        rnow = Xnow(end,idr);
        
        % concatenate the state vector
        X = vertcat(X,Xnow);
        
        % check if the ray hits the bottom
        if znow >= rcv_depth 
            break;
        end
        
       
    end
    
   figure(8)
   clf
   scatter(c_v,z_v,'o')
   grid on
   hold on
   scatter(SS.c,SS.z,'.')
   set(gca,'YDir','reverse')
   axis tight
   c_v = [];
   z_v = [];
    

    Range(n) = rnow;
    Depth(n) = znow;
    % break the loop
    if break_trig == 1
             
             fprintf('Stop scanning at launhc angle =  %f\n',theta0(n-1))
             break_trig = 0;
             break;
    end
    
     % check if the depth starts to decrease
     if n >= 2
         if ((Depth(n-1) - znow) >= 5)
             % proceed for one more angle
             break_trig = 1;
         end
         
     end
    Range(n) = rnow;
    Depth(n) = znow;
    
end
%% find angle
ind = n-1;
theta1= theta0(ind);      % longer range/ smaller angle
theta2 = theta0(ind-1);
x1= Range(ind);            % longer range
x2 = Range(ind-1);         

depth1 = Depth(ind);       % longer range/ shallower
depth2 = Depth(ind-1);    
theta_new = theta2+(rcv_depth-depth2)/(depth2 - depth1)*(theta2-theta1);
depth_new = 5000;  
%%
sinc = 0.05;
d = 1;
depht2_v = [];
ang_ad = [];
d_err = 0.04; % 3 cm
r_err = 0.04; % 3 cm
while abs(depth_new- rcv_depth) >= d_err 
    
    % initial state parameters  
    [c0,~] = find_soundspeed(z0,SS);
    X0(idt) = 0;                              % time      
    X0(idz) = z0;                             % depth
    X0(idr) = r0;                            % range
    X0(idsh) = cos(pi*theta_new/180)./c0;    % horizontal slowness
    X0(idsv) = sin(pi*theta_new/180)./c0;    % vertical slowness
   
    % initial value
    snow = 0; 
    Xnow = X0;
    znow = z0;
    rnow = r0;
    X = Xnow;
    % loop until the horizonal range exceed geodesic
    
    while (rnow <= geodesic)
        % calculate the changes of parameters within this arc length interval
        f = diff_X(snow,Xnow,SS);
    
        % update the state vector 
        Xnow = Xnow + sinc*f;
        snow = snow+sinc;
        znow = Xnow(idz);
        rnow = Xnow(idr);
%         X = vertcat(X,Xnow);
    end
%     sinc = sinc/2;
    depth_new = znow;
    
    if abs(depth_new- rcv_depth) <= d_err
        break;
    end
    
    theta2 = theta_new;
    depth2 = depth_new;
    depht2_v(end+1) = depth2;
    theta_new = theta2+(rcv_depth-depth2)/(depth2 - depth1)*(theta2-theta1);
    ang_ad(end+1) = (rcv_depth-depth2)/(depth2 - depth1)*(theta2-theta1);
    fprintf(' Depth Diff = %f\n Angle Adjust = %f \n',(rcv_depth - depth2),theta_new-theta2)
    theta1 = theta2;
    depth1 = depth2;
    
%     if sinc>=0.02
    sinc = sinc/1.5;
%     end
d = d+1;
end


 %% Supplementary function
 function SS = soundspeedprofile(ACO_depth)
ACO_lat = 22.738772;
ACO_lon = -158.006186;            % original depth MSL local u/e/n frame

fid=fopen('h306a0202.ctd');                       % 12 October 2018
D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f','headerlines',6));
fclose(fid);

pres_MSL = D(:,1);        %pressure (dbars)
temp = D(:,2);        %temperature
sal = D(:,3);         %salinity (Sp)
sal = gsw_SA_from_SP(sal,pres_MSL,ACO_lon,ACO_lat);

% Depths in the ellipsoidal coordinate; downward positive (-h)
% ACO_gh = geoidheight(ACO_lat,ACO_lon+360,'EGM96');
ACO_gh = 0;
z_hyd = ACO_depth-ACO_gh;
 

% create SS profile with the final depth at the hyd depth
% N_mean = 2.32;                                      % average geiod height over the area of coverage
N_mean = 0;
z_CTD = -(gsw_z_from_p(pres_MSL,ACO_lat) );         % relative to MSL
z_elps = z_CTD-N_mean;                              % from MSL to ellipsoid


% Truncate the lower bound incase the hyd depth is shallower
if z_elps(end) > z_hyd
    rm_ind = find(z_elps>z_hyd);
    z_elps(rm_ind) = [];
    z_CTD(rm_ind) = [];
    temp(rm_ind) = [];
    sal(rm_ind) = [];
else
%     Interpolate the lower bound in case the hyd depth is deeper
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

% calculate pressure
pres_CTD = gsw_p_from_z(-z_CTD,ACO_lat);
% 3. Creat Sound Speed Profile
SS_elps = gsw_sound_speed(sal,temp,pres_CTD);


%% Earth Flattening Transformation
R_e = 6371000;

% Dushaw and Colosi
E = z_elps./R_e;
z_flt = z_elps.*(1+(E./2)+((E.^2)./3));
SS_flt = SS_elps.*(1+E+E.^2);


% return Sound speed profile in the ellipsoidal coordinates
SS.z = z_flt;
SS.c = SS_flt;
 end
 
 function [c_z,dcdz] = find_soundspeed(znow,SS)
 c = SS.c;
 z = SS.z;
 
 dcdz_all = diff(c)./diff(z);
 % find the interval where znow falls 
 ind = find(znow >= z,1,'last');
 if ind > length(dcdz_all)
     ind = length(dcdz_all);
 end
 try
     
     c_lower = c(ind);
     del_z = znow - z(ind);
     dcdz = dcdz_all(ind);
     c_z = c_lower+del_z*dcdz;
     
 catch
     disp("error");
     disp(sprinf("Index of z = %.1f, Size dcdz_all = %.1f",ind(end),size(dcdz_all)));
     
 end
 
 end
    
function dx =  diff_X(S,x,SS)
   idt = 1; idz = 2; idr = 3; idsh = 4; idsv = 5;
   % sound speed
   [c_z,dcdz] = find_soundspeed(x(idz),SS);
   dcdr = 0;        % no horizontal gradient
   
   global c_v z_v
   c_v(end+1) = c_z;
   z_v(end+1) = x(idz);
   
   
   dx(idt) = 1/c_z;
   dx(idz) = c_z*x(idsv);
   dx(idr) = c_z*x(idsh);
   dx(idsh) = -1/(c_z.^2)*dcdr;
   dx(idsv) = -1/(c_z.^2)*dcdz;
   
end
    
   