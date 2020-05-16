close all
clear


%ACO LAT/LON
HEM_lat=22.738772;                  % Aug 2018
HEM_lon=-158.006186;                % Aug 2018

day = 27:30 ;               %  Edit
start_hour = 3;             % Edit
end_hour = 14;              % EDIT
%% 2 Load Tx data
cd /Users/testuser/Documents/ONR_RAP/Data/Tx_Rx_Output/October2018/tx_file

% create a set of file names
fname = [];

now_hour = start_hour;
now_day = day(1);
while true
    while now_hour <= 23
        if now_hour < 10
            fname_d = "tx_data_2018_10_"+string(now_day)+"_0"+string(now_hour);
        else
            fname_d = "tx_data_2018_10_"+string(now_day)+"_"+string(now_hour);
        end
       fname = vertcat(fname,fname_d);
       now_hour = now_hour+1;
       if (now_hour > end_hour)& (now_day == day(end))
            break;
       end
    end
        if (now_hour > end_hour)& (now_day == day(end))
            break;
       end
        now_hour = 0;
        now_day = now_day+1;
end

%Load TX data
tx_t = [];
tx_lat = [];
tx_lon = [];
tx_altitude = [];
tx_heading = [];
tx_xvel = [];
x_err = [];
y_err = [];
z_err = [];
for p = 1:length(fname)
    load(fname(p))
    tx_t = horzcat(tx_t,tx_data.t);      %TX time
    tx_lat = horzcat(tx_lat,tx_data.lat);   %TX lat
    tx_lon = horzcat(tx_lon,tx_data.lon);
    tx_altitude = horzcat(tx_altitude,tx_data.altitude);
    tx_heading = horzcat(tx_heading,tx_data.heading);
    tx_xvel = horzcat(tx_xvel,tx_data.x_vel);
    x_err = horzcat(x_err,tx_data.lon_err);
    y_err = horzcat(y_err,tx_data.lat_err);
    z_err = horzcat(z_err,tx_data.altitude_err);
end


%% Check surface distance/arc length
cd /Users/testuser/Documents/ONR_RAP/Data/Tx_Rx_Output/October2018/rx_file/APL_ray_tracing/HEM
fname = dir('*.mat');

x_dist = [];
arclength = [];
hyd_depth = [];
ray_trace_x_dist = [];
est_arrival =[];
act_arrival= [];
SNR= [];

for i = 1:length(fname)
    load(fname(i).name)
    est_arrival  = horzcat(est_arrival ,rx_data.est_arrival);
    act_arrival  = horzcat(act_arrival  ,rx_data.act_arrival);
    SNR  = horzcat(SNR,rx_data.SNR);
    x_dist = horzcat(x_dist,rx_data.x_dist);
    arclength = horzcat(arclength,rx_data.arclength);
    hyd_depth = horzcat(hyd_depth,rx_data.hyd_depth_kept);
    ray_trace_x_dist = horzcat(ray_trace_x_dist,rx_data.ray_trace_x_dist_kept);
    
end
% travel time perturbation with EFT
ttp = (act_arrival-est_arrival)*3600*24*1000;

%% With out EFT
cd /Users/testuser/Documents/ONR_RAP/Data/Tx_Rx_Output/October2018/rx_file/No_EFT/HEM
fname = dir('*.mat');

x_dist2 = [];
arclength2 = [];
hyd_depth2 = [];
ray_trace_x_dist2 = [];
est_arrival2 =[];
act_arrival2 = [];
SNR2= [];

for i = 1:length(fname)
    load(fname(i).name)
    est_arrival2  = horzcat(est_arrival2 ,rx_data.est_arrival);
    act_arrival2  = horzcat(act_arrival2  ,rx_data.act_arrival);
    SNR2  = horzcat(SNR2,rx_data.SNR);
    x_dist2 = horzcat(x_dist2,rx_data.x_dist);
    arclength2 = horzcat(arclength2,rx_data.arclength);
    hyd_depth2 = horzcat(hyd_depth2,rx_data.hyd_depth_kept);
    ray_trace_x_dist2 = horzcat(ray_trace_x_dist2,rx_data.ray_trace_x_dist_kept);
    
end
ttp2 = (act_arrival2-est_arrival2)*3600*24*1000;

%% Clean data 

%%%%% HEM %%%%%%
keep_ind = [];
% find the nearest arrival time (within 20 sec)
for l = 1:length(tx_t)
    dum_arrival = act_arrival;
    dum_arrival(find(dum_arrival < tx_t(l))) = [];
    h_ind = find((dum_arrival -tx_t(l))*3600*24 < 20);
    if ~isempty(h_ind)
        keep_ind(end+1) = l;
        
    end
end


tx_t =tx_t(keep_ind);
tx_lat = tx_lat(keep_ind);
tx_lon = tx_lon(keep_ind);
tx_altitude  = tx_altitude(keep_ind);
tx_heading = tx_heading(keep_ind);
tx_xvel = tx_xvel(keep_ind);
x_err = x_err(keep_ind);
y_err = y_err(keep_ind);
z_err = z_err(keep_ind);

azmth = ones(length(tx_lat),1);
% 0-360
for i=1:length(tx_lat)
    azmth(i) = azimuth(HEM_lat,HEM_lon,tx_lat(i),tx_lon(i));
end

%Boat distance from the hydrophone
for i=1:length(tx_lat)
    range(i)=dist([HEM_lat tx_lat(i)],[HEM_lon tx_lon(i)]);
end
range = range/1000;

%% Figures
figure(1)
scatter(x_dist/1000,arclength2 - arclength,[],(est_arrival2-est_arrival)*3600*24*1000,'.')
grid on
xlabel('Range (km)')
ylabel('Arc Length Difference (m)')
title('No EFT - With EFT')
set(gca,'fontsize',14)
colormap jet
%%
figure(2)
clf
set(gcf,'Units','normalized','Position',[0.1 0.6 0.5 0.5])
scatter(x_dist,(est_arrival2-est_arrival)*3600*24*1000,10,'filled')
hold on

colormap jet 
% cbar = colorbar;
caxis([0 360])

xlabel('Range (km)')
ylabel('Travel Time Difference (ms)')
title('No EFT - With EFT ')
grid on
set(gca,'fontsize',15)