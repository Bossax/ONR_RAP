% Plot data points of the entire dataset of HEM and icListen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
% HEM LAT/LON
HEM_lat=22.738772;                  % Aug 2018
HEM_lon=-158.006186;                % Aug 2018

% icListen LAT/LON
icListen_lat = 22.739153;
icListen_lon = -158.0061254;

day =[27:30]' ;            %  Edit
start_hour = 3;         % Edit
end_hour = 14;

% extract tx rx information
[tx_t_HEM,tx_lon_HEM,tx_lat_HEM,tx_heading_HEM,tx_altitude_HEM,tx_xvel_HEM,range_HEM,~,~,~,act_arrival_HEM,est_arrival_HEM,SNR_HEM] = tx_rx_extraction_Oct(day,start_hour,end_hour,'HEM');
[tx_t_icListen,tx_lon_icListen,tx_lat_icListen,tx_heading_icListen,tx_altitude_icListen,tx_xvel_icListen,range_icListen,~,~,~,act_arrival_icListen,est_arrival_icListen,SNR_icListen] = tx_rx_extraction_Oct(day,start_hour,end_hour,'icListen');

% [tx_t_HEM,tx_lon_HEM,tx_lat_HEM,tx_heading_HEM,tx_altitude_HEM,tx_xvel_HEM,range_HEM,act_arrival_HEM,est_arrival_HEM] = tx_rx_extraction_June(day,start_hour,end_hour,'HEM');
% [tx_t_icListen,tx_lon_icListen,tx_lat_icListen,tx_heading_icListen,tx_altitude_icListen,tx_xvel_icListen,range_icListen,act_arrival_icListen,est_arrival_icListen] = tx_rx_extraction_June(day,start_hour+3,end_hour+4,'HEM');
%% calculate azimuth
% HEM
HEM_azmth = [];
for i=1:length(tx_lon_HEM)
    HEM_azmth(end+1) = azimuth(HEM_lat,HEM_lon,tx_lat_HEM(i),tx_lon_HEM(i));
end


% icListen
icListen_azmth = [];
for i=1:length(tx_lon_icListen)
    icListen_azmth(end+1) = azimuth(icListen_lat,icListen_lon,tx_lat_icListen(i),tx_lon_icListen(i));
end

%% match HEM and icListen data points based on mutual transmission times
tx_t_share = intersect(tx_t_HEM,tx_t_icListen);

ind_HEM = [];
ind_ic = [];
for i = 1:length(tx_t_share)
    ind_HEM(end+1) = find(tx_t_HEM == tx_t_share(i));
    ind_ic(end+1) = find(tx_t_icListen == tx_t_share(i));
end

ttp_HEM_share = (act_arrival_HEM(ind_HEM)-est_arrival_HEM(ind_HEM))*3600*24*1000;
ttp_icListen_share = (act_arrival_icListen(ind_ic)-est_arrival_icListen(ind_ic))*3600*24*1000;
ttp_diff = ttp_HEM_share-ttp_icListen_share;


ttp_HEM = (act_arrival_HEM-est_arrival_HEM)*3600*24*1000;
ttp_ic = (act_arrival_icListen-est_arrival_icListen)*3600*24*1000;


act_trtime_HEM = (act_arrival_HEM-tx_t_HEM)*3600*24;
est_trtime_HEM = (est_arrival_HEM-tx_t_HEM)*3600*24;
act_trtime_ic = (act_arrival_icListen-tx_t_icListen)*3600*24;
est_trtime_ic = (est_arrival_icListen-tx_t_icListen)*3600*24;


lat_s = tx_lat_HEM(ind_HEM);
lon_s = tx_lon_HEM(ind_HEM);
azmth_s = HEM_azmth(ind_HEM);
act_arrival_diff = act_trtime_ic(ind_ic) - act_trtime_HEM(ind_HEM);
est_arrival_diff = est_trtime_ic(ind_ic) - est_trtime_HEM(ind_HEM);

%% time framing data for plotting
start_time = "3:05";             %EDIT
end_time = "15:30";              % EDIT
start_time = sprintf('2018-10-%s %s',string(day(1)),start_time);
end_time = sprintf('2018-10-%s %s',string(day(end)),end_time);
start_time = datenum(start_time,'yyyy-mm-dd HH:MM');
end_time = datenum(end_time,'yyyy-mm-dd HH:MM');


% crop timestamps
HEM_ind = (tx_t_HEM >= start_time)&(tx_t_HEM <= end_time);
tx_t_HEM_crp = tx_t_HEM(HEM_ind);
tx_lat_HEM_crp = tx_lat_HEM(HEM_ind);
tx_lon_HEM_crp = tx_lon_HEM(HEM_ind);
tx_altitude_HEM_crp = tx_altitude_HEM(HEM_ind);
tx_heading_HEM_crp = tx_heading_HEM(HEM_ind);
ttp_HEM_crp =ttp_HEM(HEM_ind);
range_HEM_crp = range_HEM(HEM_ind);
azmth_HEM_crp = HEM_azmth(HEM_ind); 


ic_ind = (tx_t_icListen >= start_time)&(tx_t_icListen <= end_time);
tx_t_ic_crp = tx_t_icListen(ic_ind);
tx_lat_ic_crp = tx_lat_icListen(ic_ind);
tx_lon_ic_crp = tx_lon_icListen(ic_ind);
tx_altitude_ic_crp = tx_altitude_icListen(ic_ind);
tx_heading_ic_crp = tx_heading_icListen(ic_ind);
ttp_ic_crp =ttp_ic(ic_ind);
range_icListen_crp = range_icListen(ic_ind);
azmth_ic_crp = icListen_azmth(ic_ind);

% % mutual points
keep_ind = (tx_t_share>= start_time)&(tx_t_share<= end_time);
tx_t_share_crp = tx_t_share(keep_ind);
tx_lat_mu_ic_crp = tx_lat_icListen(keep_ind);
tx_lon_mu_crp = tx_lon_icListen(keep_ind);
tx_altitude_mu_crp = tx_altitude_icListen(keep_ind);
tx_heading_mu8_crp = tx_heading_icListen(keep_ind);
ttp_diff_mu_crp =ttp_diff(keep_ind);
azmth_mu_crp = azmth_s(keep_ind);
range_share_crp = range_icListen_crp(keep_ind);

% Y limits
 ylower_lim = min([min(ttp_HEM_crp) min(ttp_ic_crp)])-1;
 yupper_lim = max([max(ttp_HEM_crp) max(ttp_ic_crp)])+1;
ylower_lim = min(ttp_ic_crp) -1;
yupper_lim = max(ttp_ic_crp)+1;

yrange = [ ylower_lim yupper_lim ];

%% Plot tx map
f1 = figure(1);
f1.Units = 'normalized';
f1.Position = [0.5 0.5 0.5 0.6];
scatter(lon_s,lat_s,20,ttp_diff,'filled')
hold on
scatter(icListen_lon,icListen_lat,200,'rh','filled')
scatter(HEM_lon,HEM_lat,200,'kp','filled')
grid on
cbar = colorbar;
cbar.Label.String = 'TTP_{HEM} - TTP_{icListen} (ms)';
colormap jet
% caxis([-4 4])
title('Trave Time Perturbation Difference (HEM - icListen)')
axis tight
xlabel('Lon')
ylabel('Lat')
set(gca,'fontsize',12)
%% icListen TTP Data
f2 = figure(2);
f2.Units = 'normalized';
f2.Position = [0.01 0.5 0.5 0.6];
scatter(tx_lon_icListen,tx_lat_icListen,20,ttp_ic,'filled')
grid on
cbar = colorbar;
cbar.Label.String = 'ms';
colormap jet
% caxis([-5 18])
title(' Travel Time Pertrubation (icListen)')
axis tight
set(gca,'fontsize',14)
%% HEM TTP Data
f3 = figure(3);
f3.Units = 'normalized';
f3.Position = [0.5 0.05 0.5 0.6];
scatter(tx_lon_HEM,tx_lat_HEM,20,ttp_HEM,'filled')
grid on
cbar = colorbar;
cbar.Label.String = 'ms';
colormap jet
title(' Travel Time Pertrubation(HEM)')
set(gca,'fontsize',14)
%% Actual Arrival
f4 = figure(4);
f4.Units = 'normalized';
f4.Position = [0.5 0.05 0.5 0.6];
scatter(lon_s,lat_s,20,act_arrival_diff*1000,'filled')
grid on
cbar = colorbar;
cbar.Label.String = 'msecond';
colormap jet
title('Actual Arrival Time Difference (icListen - HEM)')

%% Estimate Arrival
f4 = figure(5);
f4.Units = 'normalized';
f4.Position = [0.1 0.05 0.5 0.6];
scatter(tx_lon_icListen,tx_lat_icListen,20,est_trtime_ic,'filled')
grid on
cbar = colorbar;
cbar.Label.String = 'second';
colormap jet
title('Estimated Arrival Time (icListen)')
axis tight
caxis([2 19])
%}
%% Plot
f = figure(1);
clf
f.Units = 'normalized';
f.Position = [0.1 0.1 0.55 0.9];
subplot(3,1,1)
scatter(azmth_HEM_crp,ttp_HEM_crp,20,tx_t_HEM_crp,'filled')
grid on
% datetick('x')
xlim([0 360])
yticks(-20:2:20)
colormap jet
ylabel('TTP (ms)')
xlabel('Azimuth (degree)')

headline = sprintf('square CW circle\nHEM Data  Median = %.2f ms ',median(ttp_HEM_crp));
title(headline)
xticks(0:30:360)

subplot(3,1,2)
scatter(azmth_ic_crp,ttp_ic_crp,20,tx_t_ic_crp,'filled')
grid on
% datetick('x')
yticks(-20:2:20)
colormap jet
ylabel('TTP (ms)')
xlabel('Azimuth (degree)')

headline = sprintf('square CW circle\nicListenData  Median = %.2f ms ',median(ttp_ic_crp));
title(headline)
xlim([0 360])
xticks(0:30:360)

subplot(3,1,3)
scatter(azmth_mu_crp,ttp_diff_mu_crp,20,tx_t_share_crp,'filled')
grid on
% datetick('x')
yticks(-20:2:20)
% ylim([-10 10])
colormap jet
ylabel('TTP Difference (ms)')
xlabel('Azimuth (degree)')
headline = sprintf('Diffference (HEM-icListen)  Median = %.2f ms ',median(ttp_diff));
title(headline)
xlim([0 360])
xticks(0:30:360)

% subplot(3,1,3)
% scatter([azmth_HEM_crp azmth_ic_crp],[ttp_HEM_crp ttp_ic_crp],20,1:733,'filled')
% hold on
% grid on
% % datetick('x')
% yticks(-20:2:20)
% % ylim([-10 10])
% colormap jet
% ylabel('TTP Difference (ms)')
% xlabel('Azimtuh')
% headline = sprintf('CE then CCW  Median = %.2f ms ',median([ttp_HEM_crp ttp_ic_crp]));
% title(headline)
% xlim([0 360])
% xticks(0:30:360)

% subplot(4,1,4)
% scatter(tx_t_crp,range_crp,'x')
% ylabel('Surface Distance (km)')
% hold on
% yyaxis right
% scatter(tx_t_crp,tx_heading_crp,'x')
% ylim([0 360])
% yticks([0:45:360])
% grid on
% datetick('x')
% xlabel('Time')
% ylabel('heading (degree)')

%% time series
%% Plot
%{
f = figure(222);
clf
f.Units = 'normalized';
f.Position = [0.1 0.1 0.55 0.9];
subplot(3,1,1)
scatter(tx_t_HEM_crp,ttp_HEM_crp,20,range_HEM_crp,'filled')
grid on
datetick('x')
colormap jet
ylabel('TTP (ms)')
xlabel('Time')
title('HEM TTP square course without curvature correction')
cbar = colorbar ;
cbar.Label.String = 'Range (km)';
ylim([-5 15])

subplot(3,1,2)
scatter(tx_t_ic_crp,ttp_ic_crp,20,range_icListen_crp,'filled')
grid on
datetick('x')
colormap jet
ylabel('TTP (ms)')
xlabel('Time')
title('HEM TTP square course with curvature correction')
cbar = colorbar ;
cbar.Label.String = 'Range (km)';
ylim([-5 15])

subplot(3,1,3)
scatter(tx_t_share_crp,ttp_diff_mu_crp,20,range_icListen_crp,'filled')
grid on
datetick('x')
colormap jet
ylabel('TTP (ms)')
xlabel('Time')
title('Difference')
yticks(-10:2:5)
% ylim([-1 4])
cbar = colorbar ;
cbar.Label.String = 'Range (km)';
%}
