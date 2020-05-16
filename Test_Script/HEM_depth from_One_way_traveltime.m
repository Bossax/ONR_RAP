clear 
%%%%
close all

%ACO LAT/LON
HEM_lat=22.738772;                  % Aug 2018
HEM_lon=-158.006186;                % Aug 2018

% icListen LAT/LON
icListen_lat = 22.739153;
icListen_lon = -158.0061254;

% day = 27:30 ;               %  Edit
% start_hour = 3;             % Edit
% end_hour = 14;              % EDIT
% hydrophone = "HEM";    % EDIT
% % extract tx rx H
% [tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,x_err,y_err,z_err,...
%     act_arrival,est_arrival,SNR]  = tx_rx_extraction_Oct(day,start_hour,end_hour,hydrophone);


day = 7:12 ;               %  Edit
start_hour = 13;             % Edit
end_hour = 5;              % EDIT
[tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,x_err,y_err,z_err,act_arrival,est_arrival,SNR]  = tx_rx_extraction_June2017(day,start_hour,end_hour);


lim_range = 0.05;
lim_td = 2.3;
ind_kp = (find(range <= lim_range));
tvt = (act_arrival(ind_kp)-tx_t(ind_kp))*3600*24;
td_depth = -tx_altitude(ind_kp);     % ellipsoid
range = range(ind_kp);

ind_kp2 = (find(td_depth >= lim_td));
tvt = tvt(ind_kp2);
td_depth = td_depth(ind_kp2);     % ellipsoid
range = range(ind_kp2);

bm_tvt = min(tvt);
%%

% get sound speed profile
SS = soundspeedprofile(4750);

% z = integrate(v dt)
% dz = v dt
dt = 1e-1;
total_depth = [];

for ii = 1:length(tvt)
   t_now = 0;
   z_now = td_depth(ii);
    while(1)
        t_now = t_now+dt;
        if  t_now < tvt(ii)
            [SS_now,~,~] = find_soundspeed(z_now,SS);
            dz = SS_now*dt;
            z_now = z_now+dz;
        elseif t_now >= tvt(ii)
           t_now = t_now-dt;
           dt_last =  tvt(ii)-t_now;
           [SS_now,~,~] = find_soundspeed(z_now,SS);
            dz = SS_now*dt_last;
            z_now = z_now+dz;
            total_depth(ii) = z_now;
            break;
        end
        
    end
    
end

%%
scatter(td_depth-2.32,total_depth-2.32,30,'filled')
grid on
xlabel('Transducer Depth (m)')
ylabel('HEM hydrophone Depth (m)')
title({'HEM hydrophone depth relative to Mean Sea Level','from one-way acoustic travel time data'})
set(gca,'fontsize',14)
line([0 3.5],[mean(total_depth)-2.32 mean(total_depth)-2.32],'COlor','r')
text(2.8,mean(total_depth)-2.23,sprintf('Mean = %.2f m',mean(total_depth)-2.32),'fontsize',14);
%%
%
% funtions
% 
% 
% 
%%%%%%%%%%%%%
function SS = soundspeedprofile(ACO_depth)
ACO_lat = 22.73873;
ACO_lon = -158.62;
ACO_depth = -ACO_depth;            % original depth MSL local u/e/n frame

% fid=fopen('h306a0202.ctd');                       % 12 October 2018
fid=fopen('h294a0202.ctd');                       % June 2017
D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f','headerlines',6));
fclose(fid);

pres_MSL = D(:,1);        %pressure (dbars)
temp = D(:,2);        %temperature
sal = D(:,3);         %salinity (Sp)
sal = gsw_SA_from_SP(sal,pres_MSL,ACO_lon,ACO_lat);

% Depths in the ellipsoidal coordinate; downward positive (-h)
% ACO_gh = geoidheight(ACO_lat,ACO_lon+360,'EGM96');
ACO_gh = 0;
z_hyd = -(ACO_depth+ACO_gh);

% create SS profile with the final depth at the hyd depth
N_mean = 2.32;                                      % average geiod height over the area of coverage
% N_mean = 0;
z_CTD = -(gsw_z_from_p(pres_MSL,ACO_lat) );         % relative to MSL
z_elps = z_CTD-N_mean;                              % from MSL to ellipsoid


% Truncate the lower bound in case the hyd depth is shallower
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
% E = z_elps./R_e;
% z_flt = z_elps.*(1+(E./2)+((E.^2)./3));
% SS_flt = SS_elps.*(1+E+E.^2);

% Horizontal Sound Speed Gradient
dcdr = 0;

% return Sound speed profile in the ellipsoidal coordinates
SS.z = z_elps;
SS.c = SS_elps;
SS.dcdr = dcdr;
 end


function [c_z,dcdz,dcdr] = find_soundspeed(znow,SS)
 c = SS.c;
 z = SS.z;
 Hgrad = SS.dcdr;
 
 dcdz_all = diff(c)./diff(z);
 
 % find the interval where znow falls 
 ind = find(znow >= z,1,'last');
 c_lower = c(ind);
 
 
 del_z = znow - z(ind);
 
 dcdr = 0;
 dcdz = dcdz_all(ind);
 c_z = c_lower+del_z*dcdz;
 
 
 end


