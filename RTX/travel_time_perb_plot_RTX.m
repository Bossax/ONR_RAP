 % Plot Travel time perturbations RTX and POS MV
clear
close all 

%% 1 Load RTX Tx/Rx files
%ACO LAT/LON
ACO_lat=22.738772;                  % June 2017
ACO_lon=-158.006186;                %June 2017
day = 27:30 ;            %  Edit
start_hour = 3;         % Edit
end_hour = 14;
% Load TX RX file
[rtx_tx_t,rtx_tx_lon,rtx_tx_lat,rtx_tx_heading,rtx_tx_altitude,rtx_tx_xvel,rtx_range,rtx_act_arrival,rtx_est_arrival] = tx_rx_extraction_RTX(day,start_hour,end_hour,'HEM');
[tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,act_arrival,est_arrival] = tx_rx_extraction_Oct(day,start_hour,end_hour,'HEM');
%% Calculate azimuthal angles
rtx_azmth = [];
azmth = [];
for i=1:length(rtx_tx_lat)
    rtx_azmth(end+1) = azimuth(ACO_lat,ACO_lon,rtx_tx_lat(i),rtx_tx_lon(i));
end

% travel time perturbation
rtx_ttp = (rtx_act_arrival - rtx_est_arrival)*3600*24*1000;
 
% posmv
for i=1:length(tx_lat)
    azmth(end+1) = azimuth(ACO_lat,ACO_lon,tx_lat(i),tx_lon(i));
end
%
ttp = (act_arrival - est_arrival)*3600*24*1000;

%% Crop data for plotting
start_time = "20181027 00:00";       %EDIT
end_time = "20181030 15:00";         % EDIT
start_time = datenum(start_time,'yyyymmdd HH:MM');
end_time = datenum(end_time,'yyyymmdd HH:MM');

% crop timestamps
posmv_ind = (tx_t >= start_time)&(tx_t <= end_time);
tx_t_crp = tx_t(posmv_ind);
tx_lat_crp = tx_lat(posmv_ind);
tx_lon_crp = tx_lon(posmv_ind);
tx_altitude_crp = tx_altitude(posmv_ind);
tx_heading_crp = tx_heading(posmv_ind);
ttp_crp =ttp(posmv_ind);
range_crp = range(posmv_ind); 
azmth_crp = azmth(posmv_ind); 


rtx_ind = (rtx_tx_t >= start_time)&(rtx_tx_t <= end_time);
rtx_tx_t_crp = rtx_tx_t(rtx_ind);
rtx_tx_lat_crp = rtx_tx_lat(rtx_ind);
rtx_tx_lon_crp = rtx_tx_lon(rtx_ind);
rtx_tx_altitude_crp = rtx_tx_altitude(rtx_ind);
rtx_tx_heading_crp = rtx_tx_heading(rtx_ind);
rtx_ttp_crp =rtx_ttp(rtx_ind);
rtx_range_crp = rtx_range(rtx_ind);
rtx_azmth_crp = rtx_azmth(rtx_ind);


%% differencing
rm_indposmv = [];
rm_indrtx = [];
t_diffa = [];
for i = 1:length(tx_t_crp)
    [t_diff,~] = min(abs(tx_t_crp(i) - rtx_tx_t_crp)); 
    if t_diff ~= 0 
        rm_indposmv(end+1) = i;
%         t_diffa(end+1) = t_diff*3600*24;
    end
end
tx_t_4diff = tx_t_crp;
tx_t_4diff(rm_indposmv) = [];
ttp_4diff = ttp_crp;
ttp_4diff(rm_indposmv) = [];

for i = 1:length(rtx_tx_t_crp)
    [t_diff,~] = min(abs(rtx_tx_t_crp(i) - tx_t_4diff)); 
    if t_diff ~= 0 
        rm_indrtx(end+1) = i;
%         t_diffa(end+1) = t_diff*3600*24;
    end
end
rtx_tx_t_4diff = rtx_tx_t_crp;
rtx_tx_t_4diff(rm_indrtx) = [];
rtx_ttp_4diff = rtx_ttp_crp;
rtx_ttp_4diff(rm_indrtx) = [];

lon_diff = rtx_tx_lon_crp;
lon_diff(rm_indrtx) = [];
lat_diff = rtx_tx_lat_crp;
lat_diff(rm_indrtx) = [];
ttp_diff = rtx_ttp_4diff- ttp_4diff;

% Y limits
 ylower_lim = min([min(rtx_ttp_crp) min(ttp_crp)])-1;
 yupper_lim = max([max(rtx_ttp_crp) max(ttp_crp)])+1;
ylower_lim = min(rtx_ttp_crp) -1;
yupper_lim = max(rtx_ttp_crp)+1;


yrange = [ ylower_lim yupper_lim ];

%% Offset tx positions to get a better view
%{
% 2nd westaward
lat_off = 0.0065;
lon_off = 0;
start_time = '10/28/2018 18:30';
end_time = '10/28/2018 23:30';
[rtx_tx_lat,rtx_tx_lon] = tx_pos_ofset(rtx_tx_lat,rtx_tx_lon,rtx_tx_t,lat_off,lon_off,start_time,end_time);

% 3rd eastward
lat_off = 0.0130;
lon_off = 0;
start_time = '10/28/2018 23:30';
end_time = '10/29/2018 1:22';
[rtx_tx_lat,rtx_tx_lon] = tx_pos_ofset(rtx_tx_lat,rtx_tx_lon,rtx_tx_t,lat_off,lon_off,start_time,end_time);

% 1 North
lat_off = 0;
lon_off = 0.003;
start_time = '10/29/2018 01:20';
end_time = '10/29/2018 3:00';
[rtx_tx_lat,rtx_tx_lon] = tx_pos_ofset(rtx_tx_lat,rtx_tx_lon,rtx_tx_t,lat_off,lon_off,start_time,end_time);

% 2nd south
lat_off = 0;
lon_off = 0.006;
start_time = '10/29/2018 03:00';
end_time = '10/29/2018 5:30';
[rtx_tx_lat,rtx_tx_lon] = tx_pos_ofset(rtx_tx_lat,rtx_tx_lon,rtx_tx_t,lat_off,lon_off,start_time,end_time);

% 1st square
lat_off = 0;
lon_off = 0.006;
start_time = '10/29/2018 08:00';
end_time = '10/29/2018 08:45';
[rtx_tx_lat,rtx_tx_lon] = tx_pos_ofset(rtx_tx_lat,rtx_tx_lon,rtx_tx_t,lat_off,lon_off,start_time,end_time);

lat_off = -0.004;
lon_off = 0;
start_time = '10/29/2018 08:40';
end_time = '10/29/2018 09:45';
[rtx_tx_lat,rtx_tx_lon] = tx_pos_ofset(rtx_tx_lat,rtx_tx_lon,rtx_tx_t,lat_off,lon_off,start_time,end_time);


lat_off = 0;
lon_off = -0.006;
start_time = '10/29/2018 09:45';
end_time = '10/29/2018 10:45';
[rtx_tx_lat,rtx_tx_lon] = tx_pos_ofset(rtx_tx_lat,rtx_tx_lon,rtx_tx_t,lat_off,lon_off,start_time,end_time);


lat_off = 0.006;
lon_off = 0;
start_time = '10/29/2018 10:45';
end_time = '10/29/2018 11:55';
[rtx_tx_lat,rtx_tx_lon] = tx_pos_ofset(rtx_tx_lat,rtx_tx_lon,rtx_tx_t,lat_off,lon_off,start_time,end_time);

lat_off = 0;
lon_off = 0.006;
start_time = '10/29/2018 11:55';
end_time = '10/29/2018 12:30';
[rtx_tx_lat,rtx_tx_lon] = tx_pos_ofset(rtx_tx_lat,rtx_tx_lon,rtx_tx_t,lat_off,lon_off,start_time,end_time);
%}

%% circle Plot
%{
% Plot AZ
f = figure(1);
f.Units = 'normalized';
f.Position = [0.1 0.7 0.6 0.4];
clf
subplot(2,1,1)
scatter(rtx_azmth_crp,rtx_ttp_crp,[],rtx_tx_t_crp,'fill')
colormap jet
cbar = colorbar;
cbar.Label.String = 'Time';
cbdate(rtx_tx_t_crp(1:40:end),'HH:MM')
xlabel('Azimuth')
ylabel('Travel Time Perturbation (ms)')
set(gca,'Xtick',[0:30:360])
xlim([0 360])
headline = sprintf('RTX Data Median = %.2f ms',median(rtx_ttp_crp));
title(headline)
grid on
ylim([yrange])

subplot(2,1,2)
scatter(azmth_crp,ttp_crp,[],tx_t_crp,'fill')
colormap jet
cbar = colorbar;
cbar.Label.String = 'Time';
cbdate(tx_t_crp(1:60:end),'HH:MM')
grid on
xlabel('Azimuth')
ylabel('Travel Time Perturbation (ms)')
set(gca,'Xtick',[0:30:360])
xlim([0 360])
headline = sprintf('POS MV Data  Median = %.2f ms ',median(ttp_crp));
title(headline)
ylim(yrange)
%}
%% Map
f = figure(2);
clf
f.Units = 'normalized';
f.Position = [0.01 0.7 0.6 0.7];
scatter(rtx_tx_lon_crp,rtx_tx_lat_crp,[],rtx_ttp_crp,'filled')
grid on
colormap jet
 caxis([-5 20])
cbar = colorbar;
cbar.Label.String = 'Travel Time Perturbation (ms)';
set(gca,'Fontsize',11)
title('Transmission Map RTX')
hold on
scatter(ACO_lon,ACO_lat,400,'rp','filled')
%% ttp diff
f = figure(22);
clf
f.Units = 'normalized';
f.Position = [0.01 0.7 0.6 0.7];
scatter(lon_diff,lat_diff,30,ttp_diff,'filled')
grid on
colormap jet
caxis([-4 4])
cbar = colorbar;
cbar.Label.String = 'Travel Time Perturbation Difference(ms)';
title('Transmission Map Travel Time Perturbation Difference (RTX-POSMV)')
hold on
scatter(ACO_lon,ACO_lat,400,'rp','filled')
%% TTP Heading and Surface Distance

% RTX
f = figure(3);
clf
f.Units = 'normalized';
f.Position = [0.0 0.7 0.55 0.4];
subplot(2,1,1)
scatter(rtx_tx_t_crp,rtx_ttp_crp,20,rtx_tx_t_crp,'filled')
grid on
datetick('x')
yticks(-20:2:20)
ylim(yrange)
colormap jet
ylabel('Travel Time Perturbation (ms)')
xlabel('Time')
title(' RTX data')


subplot(2,1,2)
scatter(rtx_tx_t_crp,rtx_range_crp,'x')
ylabel('Surface Distance (km)')
hold on
yyaxis right
scatter(rtx_tx_t_crp,rtx_tx_heading_crp,'x')
ylim([0 360])
yticks([0:45:360])
grid on
datetick('x')
xlabel('Time')
ylabel('heading (degree)')

%{
% POSMV
f = figure(4);
clf
f.Units = 'normalized';
f.Position = [0.0 0.0 0.55 0.4];
subplot(2,1,1)
scatter(tx_t_crp,ttp_crp,20,tx_t_crp,'filled')
grid on
datetick('x')
yticks(-10:5:20)
colormap jet
ylabel('Travel Time Perturbation (ms)')
xlabel('Time')
title('POSMV data')
ylim(yrange)

subplot(2,1,2)
scatter(tx_t_crp,x_dist_crp/1000,'x')
ylabel('Surface Distance (km)')
hold on
yyaxis right
scatter(tx_t_crp,tx_heading_crp,'x')
ylim([0 360])
yticks([0:45:360])
grid on
datetick('x')
xlabel('Time')
ylabel('heading (degree)')
%}
%% ttp difference between RTX and POSMV

f = figure(5);
clf
f.Units = 'normalized';
f.Position = [0.45 0.1 0.55 0.9];
subplot(4,1,1)
scatter(tx_t_crp,ttp_crp,20,tx_t_crp,'filled')
grid on
datetick('x')
yticks(-20:5:20)
colormap jet
ylabel('TTP (ms)')
xlabel('Time')
title('POSMV data')
ylim(yrange)
headline = sprintf('POS MV Data  Median = %.2f ms ',median(ttp_crp));
title(headline)

subplot(4,1,2)
scatter(rtx_tx_t_crp,rtx_ttp_crp,20,rtx_tx_t_crp,'filled')
grid on
datetick('x')
yticks(-20:5:20)
colormap jet
ylabel('TTP (ms)')
xlabel('Time')

ylim(yrange)
headline = sprintf('RTX Data  Median = %.2f ms ',median(rtx_ttp_crp));
title(headline)

subplot(4,1,3)
scatter(tx_t_4diff,ttp_diff,20,tx_t_4diff,'filled')
grid on
datetick('x')
yticks(-20:5:20)
ylim([-10 10])
colormap jet
ylabel('TTP Difference (ms)')
xlabel('Time')
headline = sprintf('Diffference (RTX - POSMV)  Median = %.2f ms ',median(ttp_diff));
title(headline)

subplot(4,1,4)
scatter(tx_t_crp,range_crp,'x')
ylabel('Surface Distance (km)')
hold on
yyaxis right
scatter(tx_t_crp,tx_heading_crp,'x')
ylim([0 360])
yticks([0:45:360])
grid on
datetick('x')
xlabel('Time')
ylabel('heading (degree)')


cd /Users/testuser/Documents/Oct2018Cruise/Figure/RTX_calculation/traveltimeperturbation
%}
%% Tranducer map
scatter(rtx_tx_lon_crp,rtx_tx_lat_crp,[],rtx_tx_t_crp,'*')
grid on
hold on
scatter(tx_lon_crp,tx_lat_crp,[],tx_t_crp,'x')
colormap jet
legend('RTX','POSMV')

%%%%%%
function [tx_lat_out,tx_lon_out] = tx_pos_ofset(tx_lat,tx_lon,tx_t,lat_off,lon_off,start_time,end_time)
% offsets are in deg
% time of date 'mm/dd/yyyy HH:MM'
% create offseting matrices
date_format =  'mm/dd/yyyy HH:MM';
l_tx = length(tx_lat);
start_time = datenum(start_time,date_format);
end_time = datenum(end_time,date_format);
[~,ind_start] = min(abs( tx_t - start_time));
[~,ind_end] = min(abs( tx_t - end_time));
off_mat_lat = [zeros(1,ind_start-1) ones(1,ind_end-ind_start+1)*lat_off zeros(1,l_tx - ind_end)];
off_mat_lon = [zeros(1,ind_start-1) ones(1,ind_end-ind_start+1)*lon_off zeros(1,l_tx - ind_end)]; 
tx_lat_out = tx_lat + off_mat_lat;
tx_lon_out = tx_lon + off_mat_lon;

end