% Check tx files
% examine and plot data in a tx file
clear
close all

cd /Volumes/ACO_RAP_2/RAP/June2017Cruise/Tx_Rx_Output/tx_file


%ACO LAT/LON
ACO_lat = 22.738772;                  % Aug 2018
ACO_lon = -158.006186;                % Aug 2018

d = dir('tx*');
tx_t =[];
lat = [];
lon = [];
altitude = [];
heading = [];
x_vel = [];
y_vel = [];
z_vel = [];

for i = 10:length(d)
   load(d(i).name);
   tx_t = horzcat(tx_t,tx_data.t);
   lat = horzcat(lat,tx_data.lat);
   lon = horzcat(lon,tx_data.lon);
   altitude = horzcat(altitude,tx_data.altitude);
   heading = horzcat(heading,tx_data.heading);
   x_vel = horzcat(x_vel,tx_data.x_vel);
   y_vel = horzcat(y_vel,tx_data.y_vel);
   z_vel = horzcat(z_vel,tx_data.z_vel);
   
end

% remove sample
pos1 = 1;
pos2 = 0;

tx_t(1:pos1) =[];
lat(1:pos1) = [];
lon(1:pos1) = [];
altitude(1:pos1) = [];
heading(1:pos1) = [];
x_vel(1:pos1) = [];
y_vel(1:pos1) = [];
z_vel(1:pos1) = [];

tx_t(end-pos2:end) =[];
lat(end-pos2:end) = [];
lon(end-pos2:end) = [];
altitude(end-pos2:end) = [];
heading(end-pos2:end) = [];
x_vel(end-pos2:end) = [];
y_vel(end-pos2:end) = [];
z_vel(end-pos2:end) = [];

% check tx times
% today = datenum('28-10-2018','dd-mm-yyyy');
today = datenum('07-06-2017','dd-mm-yyyy');
diff = tx_t-today;
hr = floor(diff*24);
min = floor((diff*24-hr)*60);
sec = round(((diff*24-hr)*60-min)*60);
mmm = ((((diff*24-hr)*60-min)*60)-sec)*1000; % msec

% Calculate Surface Distance
x_dist = zeros(1,length(tx_t));
for i=1:length(tx_t)
    x_dist(i)=dist([ACO_lat lat(i)],[ACO_lon lon(i)]);
end

%% plot
scatter(tx_t,sec,[],mmm)
c = colorbar;
c.Label.String = 'msec';
colormap jet;
datetick('x')
grid on
sec_label = 0:1:60;
set(gca,'YTick',sec_label)
set(gca,'YTickLabel',sec_label)
%%

figure
scatter(tx_t,mmm,[],'filled')
datetick('x')
grid on
ylabel('Time offset from integer second (msec)')

%% Plot coordinate
f = figure(5);
f.Units = 'normalized';
f.Position = [0.1 0.2 0.8 0.5];
subplot(1,2,1)
scatter(lon,lat,[],tx_t)

colormap jet

grid on
title('Tx Lat/LON')
subplot(1,2,2)
scatter(tx_t,x_dist/1000,[],tx_t,'*')
grid on
ylabel('Surface Distance (km)')
xlabel('Time')
datetick('x')
title('Surface Distance')

%%
f = figure(6);
f.Units = 'normalized';
f.Position = [0.4 0.2 0.5 0.4];
subplot(2,1,1)
scatter(tx_t,heading,[],tx_t,'*')
colormap jet
ylabel('Heading')
xlabel('Time')
datetick('x')
grid on
title('Angular Motion')
subplot(2,1,2)
scatter(tx_t,z_vel,[],tx_t,'*')
colormap jet
ylabel('Yaw speed (deg/s)')
xlabel('Time')
datetick('x')
grid on



