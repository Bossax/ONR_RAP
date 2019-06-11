% pre rtx rotation
clear
close all

cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/posmv/
load posmv2911
load posmvb2911

ACO_lat=22.738772;                  % Aug 2018
ACO_lon=-158.006186;                % Aug 2018
x_dist = [];

% Transducer parameters
tr_lat = posmv.lat;
tr_lon = posmv.lon;
tr_altitude = posmv.altitude;
roll = posmv.roll;
pitch = posmv.pitch;
tr_heading = posmv.heading;
tr_heave = posmv.heave;
tr_t = posmv.t;

% granite block parameters
gr_lat = posmvb.lat;
gr_lon = posmvb.lon;
gr_altitude = posmvb.altitude;
gr_roll = posmvb.roll;
gr_pitch = posmvb.pitch;
gr_heading = posmvb.heading;
gr_heave = posmvb.heave;
gr_t = posmvb.t;

%Tx distance from ACO
for i=1:length(tr_lat)
    x_dist(i)=dist([ACO_lat tr_lat(i)],[ACO_lon tr_lon(i)]);
end
%% Plot lat/lon/altitude
figure(1)
subplot(3,1,1)
scatter(tr_lon,tr_lat,'*')
grid on
xlabel('lon')
ylabel('lat')
axis tight
title('POS MV')
subplot(3,1,2)
plot(tr_t,tr_altitude)
datetick('x')
grid on
xlabel('Time')
ylabel('Altitude (m)')
axis tight

subplot(3,1,3)
plot(tr_t,x_dist)
datetick('x')
grid on
xlabel('Time')
ylabel('Surface Distance (m)')
axis tight
%% RTX 
% examine RTX data
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/RTX/
load kilo_moana_1819_RTX.dat
A = kilo_moana_1819_RTX;
%% Timestamps
% Modified Julian Date
mjd = A(:,1);
jd = mjd +2400000.5;
rtx_t = [];
% Convert to Gregorian Date
for i = 1:length(jd)
   rtx_t(end+1) = datenum(datetime(jd(i),'ConvertFrom','juliandate'));
end

%% other parameters
rtx_lat = A(:,2)';
rtx_lon = A(:,3)';
rtx_altitude = A(:,4)';
geoid_height = A(:,8)';
%%
diff1 = [];
diff2 = [];
for i = 1:length(rtx_t)
    diff1(end+1) = abs(rtx_t(i) - tr_t(1));
    diff2(end+1) = abs(rtx_t(i) - tr_t(end));
end

[~,I1] = min(diff1);
[~,I2] = min(diff2);

% truncate
rtx_t = rtx_t(I1:I2);
rtx_lon = rtx_lon(I1:I2);
rtx_lat = rtx_lat(I1:I2);
rtx_altitude = rtx_altitude(I1:I2);


%RTX antenna distance from ACO
for i=1:length(rtx_lat)
    rtx_dist(i)=dist([ACO_lat rtx_lat(i)],[ACO_lon rtx_lon(i)]);
end


%% Plot RTX
figure(2)

subplot(3,1,1)
scatter(rtx_lon,rtx_lat,'*')
grid on
xlabel('lon')
ylabel('lat')
axis tight
title('RTX')
subplot(3,1,2)
plot(rtx_t,rtx_altitude)
datetick('x')
grid on
xlabel('Time')
ylabel('Altitude (m)')
axis tight

subplot(3,1,3)
plot(rtx_t,rtx_dist)
datetick('x')
grid on
xlabel('Time')
ylabel('Surface Distance (m)')
axis tight


%% Transformation
[rtx_tx_lat,rtx_tx_lon,rtx_tx_altitude] = transducer_pos_rtx(rtx_lat,rtx_lon,rtx_altitude,rtx_t,roll,pitch,tr_heading,tr_t);


%% 1Hz POSMV
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/posmv/1Hz_4_RTX
load('posmv_2018102911_1Hz')
    
    posmv_t = posmv.t;
    posmv_lat = posmv.lat;
    posmv_lon = posmv.lon;
    posmv_altitude = posmv.altitude';
    posmv_roll = posmv.roll;
    posmv_pitch = posmv.pitch;
    posmv_heading = posmv.heading;


%% Plot altitudes/ roll/pitch/yaw
figure(3)
subplot(3,1,1)
plot(posmv_t,posmv_altitude)
datetick('x')
grid on
xlabel('Time')
ylabel('Altitude (m)')
axis tight
title('POS MV (1 Hz)')

subplot(3,1,2)
plot(rtx_t,rtx_tx_altitude-2.31)
datetick('x')
grid on
xlabel('Time')
ylabel('Altitude (m)')
axis tight
title('RTX (1 Hz)')
ylim([-7 0])

subplot(3,1,3)
plot(posmv_t,posmv_altitude)
hold on
plot(rtx_t,rtx_tx_altitude-2.31)
datetick('x')
grid on
xlabel('Time')
ylabel('Altitude (m)')
axis tight
legend('POSMV','RTX')
ylim([-8 0])
yticks(-8:1:0)
%% lat lon plot
figure(4)
scatter(rtx_tx_lon,rtx_tx_lat)
grid on
hold on
scatter(tr_lon,tr_lat)
hold on
scatter(rtx_lon,rtx_lat)
legend('transducer rtx','transducer posmv','antenna rtx')
%% Displacement different)
for k =1:length(rtx_t)
    [~,I] = min(abs(tr_t - rtx_t(k)));
    [rtx_x,rtx_y,~] = deg2utm(rtx_tx_lat(k),rtx_tx_lon(k));
    [tr_x,tr_y,~] = deg2utm(tr_lat(I),tr_lon(I));
    dis_diff(k) = sqrt((rtx_x-tr_x)^2+(rtx_y-tr_y)^2+(rtx_tx_altitude(k)-tr_altitude(I))^2);
end
mean_dis = mean(dis_diff);
std_dis = std(dis_diff);
min_x = mean_dis - 4*std_dis;
max_x = mean_dis +4*std_dis;
x = min_x:(max_x-min_x)/100:max_x;
f = exp(-(x-mean(dis_diff)).^2./(2*std(dis_diff)));
% Distribution of distance differences between posmv and rtx
figure(5)
plot(x,f)
line([mean_dis mean_dis],[0 1],'Color','r')
grid on
xlabel('Distance difference (m)')
ylabel('f(x)')
title('Distribution')
text(mean_dis-std_dis/2,0.5,'mean = 12.28 m')
text(mean_dis+std_dis ,0.3,'std = 1.068 m')
line([mean_dis+std_dis mean_dis+std_dis],[0 0.58],'Color','k')
