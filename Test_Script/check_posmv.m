% check posmv 
% examine and plot data in a posmv file
close all
clear
%ACO LAT/LON
ACO_lat=22.738772;                  % Aug 2018
ACO_lon=-158.006186;                % Aug 201
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/posmv/all
fname = ['posmv_2018102908.mat';'posmv_2018102909.mat';'posmv_2018102910.mat']
%% Read in all var
t=[];
lat=[];
lon=[];
altitude=[];
heading=[];
accel_lon =[];
accel_tran = [];
accel_down = [];
rate_lon = [];
rate_tran = [];
rate_down = [];
x_vel = [];
y_vel = [];
z_vel = [];
roll = [];
pitch = [];
heave = [];
altitude_err = [];
lat_err = [];
lon_err = [];
roll_err = [];
pitch_err = [];
heading_err = [];
t_err = [];

 for i = 1:3
    load(fname(i,:));
    t=vertcat(t,posmv.t);
    lat=vertcat(lat,posmv.lat);
    lon=vertcat(lon,posmv.lon);
    altitude=vertcat(altitude,posmv.altitude);
    heading=vertcat(heading,posmv.heading);
    accel_lon =vertcat(accel_lon,posmv.accel_lon);
    accel_tran =vertcat(accel_tran,posmv.accel_tran);
    accel_down =vertcat(accel_down,posmv.accel_down);
    rate_lon =vertcat(rate_lon,posmv.rate_lon);
    rate_tran =vertcat(rate_tran,posmv.rate_tran);
    rate_down =vertcat(rate_down,posmv.rate_down);
    x_vel = vertcat(x_vel,posmv.x_vel);
    y_vel = vertcat(y_vel,posmv.y_vel);
    z_vel = vertcat(z_vel,posmv.z_vel);
    roll = vertcat(roll,posmv.roll);
    pitch = vertcat(pitch,posmv.pitch);
    heave =  vertcat(heave,posmv.heave);
    altitude_err =  vertcat(altitude,posmv.altitude);
    lat_err =  vertcat(lat_err,posmv.lat_err);
    lon_err =  vertcat(lon_err,posmv.lon_err);
    roll_err =  vertcat(roll_err,posmv.roll_err);
    pitch_err =  vertcat(pitch_err,posmv.pitch_err);
    heading_err =  vertcat(heading_err,posmv.heading_err);
    t_err =  vertcat(t_err,posmv.t_err);
    
 end
altitude=altitude-2.31;           %WSG84 Geoid height relative to MSL (at ACO)

%% massage data
%{ 
get rid off out-of-bound data points

date_ref = mean(t);
datestr(date_ref);
one_hr = 1/24;
ind = find(abs(t-date_ref) > one_hr);

t(ind) = [];
lat(ind) = [];
lon(ind) = [];
altitude(ind) = [];
heading(ind) = [];
accel_lon(ind) = [];
accel_tran(ind) = [];
accel_down(ind) = [];
x_vel(ind) = [];
y_vel(ind) = [];
z_vel(ind) = [];
rate_down(ind) = [];
rate_lon(ind) = [];
rate_tran(ind) = [];
roll(ind) = [];
pitch(ind) = [];
%}

% Calculate Surface Distance
x_dist = zeros(1,length(t));
for i=1:length(t)
    x_dist(i)=dist([ACO_lat lat(i)],[ACO_lon lon(i)]);
end

%% Plot
f1 = figure(1);
f1.Units = 'normalized';
f1.Position = [0.01 0.6 0.4 0.3];
clf
% Scoping var
start_t = datenum('29-10-2018 17:30:00','dd-mm-yyyy HH:MM:SS');
end_t = datenum('29-10-2018 19:25:00','dd-mm-yyyy HH:MM:SS');
ind = find(t >= start_t & t<=end_t);

% Coordinates
subplot(2,1,1)
scatter(lon(ind),lat(ind),[],t(ind),'filled');
grid on
ylabel('Lat')
xlabel('Lon')
title('Ship Trajectory')

subplot(2,1,2)
scatter(t(ind),x_dist(ind)/1000,[],t(ind))
grid on
datetick('x')
ylabel('surface distance (km) ')
colormap jet

% Angular Motion
f2 = figure(2);
f2.Units = 'normalized';
f2.Position = [0.5 0.6 0.4 0.3];
clf

subplot(2,1,1)
scatter(t(ind),heading(ind),[],t(ind),'filled');
grid on
ylabel('Heading')
xlabel('Time')
title('Yaw Motion')
datetick('x')
yticks([0:45:360])

subplot(2,1,2)
scatter(t(ind),rate_down(ind),[],t(ind))
grid on
datetick('x')
ylabel('Yaw Velocity (deg/sec) ')
colormap jet

% Horizontal Motion
horz_vel=x_vel.^2+y_vel.^2;
horz_accel = accel_lon.^2 + accel_tran.^2;

f3 = figure(3);
f3.Units = 'normalized';
f3.Position = [0.01 0.2 0.4 0.4];
clf

subplot(2,1,1)
scatter(t(ind),horz_vel(ind),[],t(ind),'filled');
grid on
ylabel('Horizontal Velocity (m/s)')
xlabel('Time')
title('Horizontal Motion')
datetick('x')

subplot(2,1,2)
scatter(t(ind),horz_accel(ind),[],t(ind))
grid on
datetick('x')
ylabel('Horizontal Acceleration (m/s^{2}) ')
colormap jet

% Roll Pitch
f4 = figure(4);
f4.Units = 'normalized';
f4.Position = [0.5 0.2 0.4 0.3];
clf

subplot(2,1,1)
scatter(t(ind),roll(ind),[],t(ind),'filled');
grid on
ylabel('Roll (deg)')
title('Roll/Pitch Motions')
datetick('x')

subplot(2,1,2)
scatter(t(ind),pitch(ind),[],t(ind))
grid on
datetick('x')
ylabel('Pitch (deg) ')
colormap jet

% altitude
f5 = figure(5);
f5.Units = 'normalized';
f5.Position = [0.4 0.5 0.4 0.3];
clf

subplot(2,1,1)
scatter(t(ind),altitude(ind),[],t(ind),'filled');
grid on
ylabel('Altitude (m)')
title('Z Motions')
datetick('x')

subplot(2,1,2)
scatter(t(ind),heave(ind),[],t(ind))
grid on
datetick('x')
ylabel('Heave (m) ')
colormap jet
%%
%Rx data
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/rx_file/spin/10km
load('rx_data_spin_10km')
est_arrival = rx_data.est_arrival;
act_arrival =rx_data.act_arrival;
indc = [];
for i =1:length(est_arrival)
   
    set = find(abs(est_arrival(i)-act_arrival)*3600*24*1000 < 50);
    if isempty(set)
        indc(end+1) = i;
    end
    
end
est_arrival(indc) = [];
tt_perb = (act_arrival - est_arrival)*3600*24;
ind = find(est_arrival >= start_t & est_arrival<=end_t);

f6 = figure(6);
f6.Units = 'normalized';
f6.Position = [0.4 0.5 0.5 0.3];
clf

scatter(est_arrival(ind),tt_perb(ind)*1000,[],est_arrival(ind),'filled')
datetick('x')
grid on
colormap jet
ylabel('TT perturbation (ms)')
title('Square CW')

%% Save figures
%{
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Document/Fig/Square/CW/5th_turn
figure(1)
saveas(f1,'Ship Trajectory.jpg')

figure(2)
saveas(f2,'Yaw Motion.jpg')

figure(3)
saveas(f3,'Horizontal Velocity.jpg')

figure(4)
saveas(f4,'Roll and Pitch.jpg')

figure(5)
saveas(f5,'Z Motions.jpg')

figure(6)
saveas(f6,'TT perturbation.jpg')
%}

