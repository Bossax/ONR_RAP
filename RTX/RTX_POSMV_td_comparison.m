% Calculate offsets of the transducer position dfference between RTX (raw)
% and POS MV and plot several plots
% 1. Load RTX files and store parameters
% 2. Load POS MV file
% 3. Coordinate Transformation
% 4. Plot
% 5. Calculate time lags between POSMV and RTX and the corresponding
% displacement due to the time lags (and velocities)
clear
close all
%% Download RTX file
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/RTX/
load kilo_moana_1819_RTX.dat
A = kilo_moana_1819_RTX;
rtx_lat = transpose(A(:,2));
rtx_lon = transpose(A(:,3));
rtx_altitude = transpose(A(:,4));
rtx_altitude = rtx_altitude - 2.31;     % subtract geoid height
%% Format timestamps
% Modified Julian Date 
mjd = A(:,1);
jd = mjd +2400000.5;
rtx_t_all = [];
% Convert to Gregorian Date
rtx_t_all = datenum(datetime(jd,'ConvertFrom','juliandate'));
rtx_t_all = rtx_t_all';

% Adjust GPS time to UTC time
rtx_t_all  = rtx_t_all  - 2/(3600*24); % subtract 2 second from the GPS time

%% Download POS MV files
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/posmv/all
d_posmv = dir('posmv*');
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/posmv/granite_block/all
d_posmv_gn = dir('gb*');

%% Loop over all files and read the tx times
close all
% find nearest corresponding rtx times
% calculate transducer position
% pack data and save
t_1sec = 0;
ind_1sec_rtx = 0;
t_diff = [];
% create space for transducer position vectors
% 1Hz
% RTX
rtx_tx_lat = [];
rtx_tx_lon = [];
rtx_tx_altitude = [];
rtx_x_tx_gb = [];
rtx_y_tx_gb = [];
rtx_t = [];
rtx_x_ant_gb = [];
rtx_y_ant_gb = [];
rtx_z_ant_gb = [];

% POSMV 
posmv_t = [];
posmv_lat = [];
posmv_lon = [];
posmv_altitude = [];
posmv_heading = [];
posmv_x_gb = [];
posmv_y_gb = [];

posmv_xvel = [];
posmv_yvel = [];
posmv_zvel = [];
posmv_xacc = [];
posmv_yacc = [];
posmv_zacc = [];

% POSMV antenna
posmv_x_ant_gb = [];
posmv_y_ant_gb = [];
posmv_z_ant_gb = [];

% RTX td offset from POSMV's
offset_rtx_posmv_tx_s3 = [];
offset_rtx_tx_gn_s3 = [];
offset_posmv_tx_gn_s3 = [];
offset_rtx_tx_ant_s3 = [];
offset_rtx_posmv_ant_s3=[];
% 10 Hz
% POSMV
posmv_t10hz = [];
posmv_lat10hz  = [];
posmv_lon10hz  = [];
posmv_altitude10hz  = [];
posmv_heading10hz  = [];

posmv_xvel10hz  = [];
posmv_yvel10hz  = [];
posmv_zvel10hz  = [];
posmv_xaccel10hz  = [];
posmv_yaccel10hz  = [];
posmv_zaccel10hz  = [];


% coordinates of devices with respect to the granite block
ant = [6.439;6.505;-27.761]; % from KM coordinate system May 2018 spreasheet
tx5 = [0.568;19.599;0.715];
ant2tx = tx5-ant;    

% rotation matrix: Global, UTM coordiante to unrotated ship coordinate
thetax_gb = 180/180*pi;
thetaz_gb = -90/180*pi;
Rxgb = [1 0 0;0 cos(thetax_gb) -sin(thetax_gb);0 sin(thetax_gb) cos(thetax_gb)];
Rzgb = [cos(thetaz_gb) -sin(thetaz_gb) 0 ; sin(thetaz_gb) cos(thetaz_gb) 0; 0 0 1];
Rgb_s = Rxgb*Rzgb;
Rs_gb = transpose(Rgb_s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop over posmv files
% EDIT start file and end file
for ii = 1:length(d_posmv)
    ii
    cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/posmv/all
    % extract posmv data
    % transducer
    fname = d_posmv(ii).name;
    load (fname)
    posmv_latd = posmv.lat;
    posmv_lond = posmv.lon;
    posmv_altituded =posmv.altitude;
    posmv_rolld = posmv.roll;
    posmv_pitchd = posmv.pitch;
    posmv_yawd = posmv.heading;
    posmv_td = posmv.t;
    
    ind_l = [];
    % eliminate outliers
    for l = 2:length(posmv_td)
        if (posmv_td(l) - posmv_td(l-1)) < 0 
            ind_l = l;
            break;
        end
    end
    posmv_latd(ind_l) = [];
    posmv_lond(ind_l) = [];
    posmv_altituded(ind_l) = [];
    posmv_rolld(ind_l) = [];
    posmv_pitchd(ind_l) = [];
    posmv_yawd(ind_l) = [];
    posmv_td(ind_l) = [];
    
    
    % store raw 10Hz data
    posmv_t10hz = horzcat(posmv_t10hz,posmv_td');
    posmv_lat10hz  = horzcat(posmv_lat10hz,posmv_latd');
    posmv_lon10hz  = horzcat(posmv_lon10hz,posmv_lond');;
    posmv_altitude10hz  = horzcat(posmv_altitude10hz,posmv_altituded');;
    
    
    
    
    % vessel dynamics
    posmv_xveld = posmv.x_vel;
    posmv_yveld = posmv.y_vel;
    posmv_zveld = posmv.z_vel;
    posmv_xacceld = posmv.accel_lon;
    posmv_yacceld = posmv.accel_tran;
    posmv_zacceld = posmv.accel_down;
    
    
    posmv_xvel10hz  = horzcat(posmv_xvel10hz,posmv_xveld');
    posmv_yvel10hz  = horzcat(posmv_yvel10hz,posmv_yveld');
    posmv_zvel10hz  = horzcat(posmv_yvel10hz,posmv_yveld');
    posmv_xacc10hz  = horzcat(posmv_xaccel10hz,posmv_xacceld');
    posmv_yacc10hz  = horzcat(posmv_yaccel10hz,posmv_yacceld');
    posmv_zacc10hz  = horzcat(posmv_zaccel10hz,posmv_zacceld');
    
    
    
    % granite block
    cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/posmv/granite_block/all
    fname2 = d_posmv_gn(ii).name;
    load(fname2)
    posmv_lat_gnd = posmv_gb.lat;
    posmv_lon_gnd = posmv_gb.lon;
    posmv_altitude_gnd =posmv_gb.altitude;
    posmv_t_gnd = posmv_gb.t;
    
    
    % match rtx and posmv data points
    ind = [];           % posmv indices to subtract data
    rtx_t_w = [];       % times
    ind_1sec_rtx = [];   % eliminated indicies
    % window rtx timeseries to a specific interval
    start_t = posmv_td(1);
    end_t = start_t+1/24;
    [~,start_ind] = min(abs(rtx_t_all-start_t));
    [~,end_ind] = min(abs(rtx_t_all-end_t));
    rtx_t_w = rtx_t_all(start_ind:end_ind);
    
    % antenna position
    rtx_lat_w = rtx_lat(start_ind:end_ind);
    rtx_lon_w = rtx_lon(start_ind:end_ind);
    rtx_altitude_w = rtx_altitude(start_ind:end_ind);
    
    % find the neast neighbor posmv data points of the corresponding RTX data points
    for k = 1:length(rtx_t_w)
       [t_difff,I] =  min(abs(rtx_t_w(k) - posmv_td));
       if t_difff*3600*24 <= 0.5    % if time difference is larger than 0.5 sec
            ind(end+1) = I;         % extract ship attitudes from POS-MV 
            t_diff(end+1) = t_difff*3600*24;
       else
           ind_1sec_rtx(end+1) = k;   % discard a point whose nearest neighbor POSMV data point is farther than 0.5 sec
            t_1sec = t_1sec+1;
        end
    end
    
    % eliminate RTX data points without POS MV match
    if length(ind_1sec_rtx) > 0
        rtx_t_w(ind_1sec_rtx) = [];
    end
    rtx_roll_w = posmv_rolld(ind);
    rtx_pitch_w = posmv_pitchd(ind);
    rtx_heading_w = posmv_yawd(ind);
    
    % calculate transducer position (global coordinate)/ antenna position
    [tx_lat_rtx,tx_lon_rtx,tx_altitude_rtx] = transducer_pos_rtx(rtx_lat_w,rtx_lon_w,rtx_altitude_w,rtx_t_w,posmv_rolld,posmv_pitchd,posmv_yawd,posmv_td);
%     [pos_ant_lat,pos_ant_lon,pos_ant_altitude] = ant_pos_posmv(posmv_latd,posmv_lond,posmv_altituded,posmv_td,posmv_rolld,posmv_pitchd,posmv_yawd);
    
    % convert lat lon to UTM x,y,z POSMV
    [posmv_x_d,posmv_y_d,utmzone] = deg2utm(posmv_latd(ind),posmv_lond(ind));
    posmv_z_d = posmv_altituded(ind);
    posmv_t_d = posmv_td(ind);
    posmv_tx5_gb = [posmv_x_d'; posmv_y_d' ; posmv_z_d'];
    
%     [posmv_x_ant_d,posmv_y_ant_d,utmzone] = deg2utm(pos_ant_lat(ind),pos_ant_lon(ind));
%     posmv_z_ant_d = pos_ant_altitude(ind);
%     posmv_ant_gb = [posmv_x_ant_d'; posmv_y_ant_d' ; posmv_z_ant_d'];
%     
    
    [posmv_x_gnd,posmv_y_gnd,utmzone] = deg2utm(posmv_lat_gnd(ind),posmv_lon_gnd(ind));
    posmv_z_gnd = posmv_altitude_gnd(ind);
    posmv_t_d = posmv_t_gnd(ind);
    posmv_gn_gb = [posmv_x_gnd'; posmv_y_gnd' ; posmv_z_gnd'];
    
    
    % RTX transducer/antenna positions in UTM
    [rtx_x_w,rtx_y_w,utmzone] = deg2utm(tx_lat_rtx,tx_lon_rtx);
    rtx_z_w = tx_altitude_rtx;
    rtx_tx5_gb = [rtx_x_w';rtx_y_w';rtx_z_w'];
    
    % antenna position
    [rtx_x_ant,rtx_y_ant,utmzone] = deg2utm(rtx_lat_w,rtx_lon_w);
     rtx_ant_gb = [rtx_x_ant';rtx_y_ant'; rtx_altitude_w];
     
    % transform lat/lon in the global UTM coordinate back to ship locally level frame
    tx5_rtx_s3_refant = [];
    tx5_rtx_s3_reftx = [];
    tx5_rtx_s3_refgn = [];
    tx5_posmv_s3_refgn = [];
    ant_posmv_s3_refgn = [];
    ant_rtx_s3_refgn = [];
    
    rtx_roll_w2 = rtx_roll_w/180*pi; 
    rtx_pitch_w2 = rtx_pitch_w/180*pi; 
    rtx_heading_w2 = rtx_heading_w/180*pi; 
    
    % Perform coordinate transformation from UTM to ship locally level frame
    for k = 1:length(rtx_roll_w)
        Rz = [cos(rtx_heading_w2(k)) -sin(rtx_heading_w2(k)) 0; sin(rtx_heading_w2(k)) cos(rtx_heading_w2(k)) 0; 0 0 1];    % yaw
        Ry = [cos(rtx_pitch_w2(k)) 0 sin(rtx_pitch_w2(k)); 0 1 0; -sin(rtx_pitch_w2(k)) 0 cos(rtx_pitch_w2(k))]; %pitch
        Rx = [1 0 0; 0 cos(rtx_roll_w2(k)) -sin(rtx_roll_w2(k)); 0 sin(rtx_roll_w2(k)) cos(rtx_roll_w2(k))]; % roll
        R = Rz*Ry*Rx;     
        % the offset of RTX's transducer from POSMV's
        tx5_rtx_s3_handle3 = R^-1*Rs_gb*(rtx_tx5_gb(:,k) - rtx_ant_gb(:,k));
        tx5_rtx_s3_handle = R^-1*Rs_gb*(rtx_tx5_gb(:,k) - posmv_tx5_gb(:,k));
        tx5_rtx_s3_handle2 = R^-1*Rs_gb*(rtx_tx5_gb(:,k) - posmv_gn_gb(:,k));
        tx5_posmv_s3_handle = R^-1*Rs_gb*(posmv_tx5_gb(:,k) - posmv_gn_gb(:,k));
%         ant_posmv_s3_handle =  R^-1*Rs_gb*(posmv_ant_gb(:,k) - posmv_gn_gb(:,k));
        ant_rtx_s3_handle = R^-1*Rs_gb*(rtx_ant_gb(:,k) - posmv_gn_gb(:,k));
        
        tx5_rtx_s3_refant = horzcat(tx5_rtx_s3_refant,tx5_rtx_s3_handle3);;
        tx5_rtx_s3_reftx = horzcat(tx5_rtx_s3_reftx,tx5_rtx_s3_handle);
        tx5_rtx_s3_refgn = horzcat(tx5_rtx_s3_refgn,tx5_rtx_s3_handle2);
        tx5_posmv_s3_refgn = horzcat(tx5_posmv_s3_refgn,tx5_posmv_s3_handle);
%         ant_posmv_s3_refgn = horzcat(ant_posmv_s3_refgn,ant_posmv_s3_handle);;
        ant_rtx_s3_refgn = horzcat(ant_rtx_s3_refgn,ant_rtx_s3_handle);;
    end

   % store variables
   % RTX
   % transducer
    rtx_tx_lat = horzcat(rtx_tx_lat,tx_lat_rtx');
    rtx_tx_lon = horzcat(rtx_tx_lon,tx_lon_rtx');
    rtx_tx_altitude = horzcat(rtx_tx_altitude,tx_altitude_rtx');
    rtx_x_tx_gb = horzcat(rtx_x_tx_gb,rtx_x_w');
    rtx_y_tx_gb = horzcat(rtx_y_tx_gb,rtx_y_w');
    rtx_t = horzcat(rtx_t,rtx_t_w);
    % antenna
    rtx_x_ant_gb = horzcat(rtx_x_ant_gb,rtx_x_ant');
    rtx_y_ant_gb = horzcat(rtx_y_ant_gb,rtx_y_ant');
    rtx_z_ant_gb = horzcat(rtx_z_ant_gb,rtx_altitude_w);
    
%     posmv_x_ant_gb = horzcat(posmv_x_ant_gb,posmv_x_ant_d');
%     posmv_y_ant_gb = horzcat(posmv_y_ant_gb,posmv_y_ant_d');
%     posmv_z_ant_gb = horzcat(posmv_z_ant_gb,posmv_z_ant_d');
   % POSMV concatinate
    posmv_t = horzcat(posmv_t,posmv_td(ind)');
    posmv_lat = horzcat(posmv_lat ,posmv_latd(ind)');
    posmv_lon = horzcat(posmv_lon,posmv_lond(ind)');
    posmv_altitude = horzcat(posmv_altitude,posmv_altituded(ind)');
    posmv_x_gb = horzcat(posmv_x_gb,posmv_x_d');
    posmv_y_gb = horzcat(posmv_y_gb,posmv_y_d');
    posmv_heading = horzcat(posmv_heading,posmv_yawd(ind)');
    
    posmv_xvel = horzcat(posmv_xvel,posmv_xveld(ind)');
    posmv_yvel = horzcat(posmv_yvel,posmv_yveld(ind)');
    posmv_zvel = horzcat(posmv_zvel,posmv_zveld(ind)');
    posmv_xacc = horzcat(posmv_xacc,posmv_xacceld(ind)');
    posmv_yacc = horzcat(posmv_yacc,posmv_yacceld(ind)');
    posmv_zacc = horzcat(posmv_zacc,posmv_zacceld(ind)');
    
    offset_rtx_tx_ant_s3 = horzcat(offset_rtx_tx_ant_s3, tx5_rtx_s3_refant) ;
    offset_rtx_posmv_tx_s3 = horzcat(offset_rtx_posmv_tx_s3,tx5_rtx_s3_reftx);
    offset_rtx_tx_gn_s3 = horzcat(offset_rtx_tx_gn_s3,tx5_rtx_s3_refgn);
    offset_posmv_tx_gn_s3 = horzcat(offset_posmv_tx_gn_s3,tx5_posmv_s3_refgn);
%     offset_rtx_posmv_ant_s3= horzcat(offset_rtx_posmv_ant_s3,ant_rtx_s3_refgn-ant_posmv_s3_refgn);
    
end
%% time frame
tstart = datenum('20181027 03:00:00','yyyymmdd HH:MM:SS');
tend= datenum('20181029 03:00:00','yyyymmdd HH:MM:SS');
[~,indstart] = min(abs(posmv_t - tstart));
[~,indend] = min(abs(posmv_t - tend));


%% transducer offsets timeseries
%{
% global frame
horz_dis_vector = [rtx_x_tx_gb - posmv_x_gb ; rtx_y_tx_gb - posmv_y_gb];
horz_dis = zeros(1,length(horz_dis_vector));
for jj = 1:length(horz_dis_vector)
    horz_dis(jj) = sqrt(horz_dis_vector(1,jj)^2 +horz_dis_vector(2,jj)^2);
end
figure(1)
plot(rtx_t,horz_dis)
grid on
ylabel('meter')
datetick('x')
xlabel('Time')
axis tight
headline1 = "Horizontal Distance Offset " +string(datestr(rtx_t(1),'mm/dd HH:MM'))+" to "+string(datestr(rtx_t(end),'HH:MM'));
headline = headline1 + sprintf('\n RMS offset = %.3f m ',rms(horz_dis));
title(headline)

% locally level ship frame
horz_dis_s3 = zeros(1,length(offset_rtx_posmv_tx_s3));
for jj = 1:length(horz_dis_s3)
    horz_dis_s3(jj) = sqrt(offset_rtx_posmv_tx_s3(1,jj)^2 +offset_rtx_posmv_tx_s3(2,jj)^2);
end
%}
%% offsets RTX td /POSMV td

f = figure(2);
f.Units = 'normalized';
f.Position = [0.01 0.9 0.5 0.5];
subplot(3,1,1)
scatter(rtx_t(indstart:indend),offset_rtx_posmv_tx_s3(1,indstart:indend),1,rtx_t(indstart:indend),'.')
colormap jet
grid on
ylabel('Bow-Stern Offset (m)')
datetick('x')
axis tight
ylim([-1 1])
yticks(-1:.2:1)
headline = sprintf('RTX transducer offset from POSMV transducer\nMedian = %.3f m',median(offset_rtx_posmv_tx_s3(1,indend:indend)));
title(headline)

subplot(3,1,2)
scatter(rtx_t(indstart:indend),offset_rtx_posmv_tx_s3(2,indstart:indend),1,rtx_t(indstart:indend),'.')
colormap jet
grid on
ylabel('Startboard-Port Offset (m)')
datetick('x')
axis tight
ylim([-1 1])
headline = sprintf('Median = %.3f m',median(offset_rtx_posmv_tx_s3(2,indstart:indend)));
title(headline)

subplot(3,1,3)
scatter(rtx_t(indstart:indend),offset_rtx_posmv_tx_s3(3,indstart:indend),1,rtx_t(indstart:indend),'.')
colormap jet
grid on
ylabel('Z Offset (m)')
datetick('x')
axis tight
ylim([0 3])
headline = sprintf('Median = %.3f m',median(offset_rtx_posmv_tx_s3(3,indend:indend)));
title(headline)
%}
%% offset RTX td/POS MV GN
%{
f = figure(22);
f.Units = 'normalized';
f.Position = [0.5 0.9 0.5 0.5];
subplot(3,1,1)
scatter(rtx_t,offset_rtx_posmv_gn_s3(1,:),[],rtx_t,'.')
colormap jet
grid on
ylabel('Bow-Stern Offset (m)')
datetick('x')
axis tight
ylim([-15 5])
headline = sprintf('RTX transducer offset from POSMV granite block\nMedian = %.3f m',median(offset_rtx_posmv_gn_s3(1,:)));
title(headline)

subplot(3,1,2)
scatter(rtx_t,offset_rtx_posmv_gn_s3(2,:),[],rtx_t,'.')
colormap jet
grid on
ylabel('Startboard-Port Offset (m)')
datetick('x')
axis tight
ylim([10 30])
headline = sprintf('Median = %.3f m',median(offset_rtx_posmv_gn_s3(2,:)));
title(headline)


subplot(3,1,3)
scatter(rtx_t,offset_rtx_posmv_gn_s3(3,:),[],rtx_t,'.')
colormap jet
grid on
ylabel('Z Offset (m)')
datetick('x')
axis tight
ylim([0 5])
headline = sprintf('Median = %.3f m',median(offset_rtx_posmv_gn_s3(3,:)));
title(headline)
%}

%% POS Mv td to GN
%{
f = figure(23);
f.Units = 'normalized';
f.Position = [0.01 0.01 0.5 0.5];
subplot(4,1,1)
scatter(rtx_t,offset_posmv_tx_gn_s3(1,:),[],rtx_t,'.')
colormap jet
grid on
ylabel('Bow-Stern Offset (m)')
datetick('x')
axis tight
ylim([0 1.2])
headline = sprintf('POSMV transducer offset from POSMV granite block \n Median = %.3f m',median(offset_posmv_tx_gn_s3(1,:)));
title(headline)


subplot(4,1,2)
scatter(rtx_t,offset_posmv_tx_gn_s3(2,:),[],rtx_t,'.')
colormap jet
grid on
ylabel('Startboard-Port Offset (m)')
datetick('x')
axis tight
ylim([19 20])
headline = sprintf('Median = %.3f m',median(offset_posmv_tx_gn_s3(2,:)));
title(headline)


subplot(4,1,3)
scatter(rtx_t,offset_posmv_tx_gn_s3(3,:),[],rtx_t,'.')
colormap jet
grid on
ylabel('Z Offset (m)')
datetick('x')
axis tight
ylim([-10 10])
headline = sprintf('Median = %.3f m',median(offset_posmv_tx_gn_s3(3,:)));
title(headline)

tx_posmv_dis = sqrt(offset_posmv_tx_gn_s3(1,:).^2+offset_posmv_tx_gn_s3(2,:).^2+offset_posmv_tx_gn_s3(3,:).^2);
subplot(4,1,4)
scatter(rtx_t,tx_posmv_dis,[],rtx_t,'.')
colormap jet
grid on
ylabel('Absolute Offset (m)')
datetick('x')
axis tight
ylim([19.62 19.68])
headline = sprintf('Median = %.3f m',median(tx_posmv_dis ));
title(headline)
%}
%% Map/Velocity/ heading
%{
aco_loc = [22.7388,-158.006186];
f = figure(3);
f.Units = 'normalized';
f.Position = [0.55 0.9 0.4 0.5];
scatter(rtx_tx_lon,rtx_tx_lat,[],rtx_t)
colormap jet
grid on
c = colorbar;
cbdate('HH:SS')
c.Label.String = 'Transmission Time ';
hold on
scatter(aco_loc(2),aco_loc(1),300,'kp','filled')

%
f = figure(4);
f.Units = 'normalized';
f.Position = [0.5 0.1 0.55 0.25];
scatter(rtx_t,posmv_heading,[],rtx_t,'.')
colormap jet
xlabel('Time')
ylabel('Heading (degree)')
datetick('x','HH:MM')
grid on
yticks(0:30:360)
ylim([0 360])


f = figure(5);
f.Units = 'normalized';
f.Position = [0.01 0.1 0.55 0.4];
subplot(3,1,1)
scatter(rtx_t,posmv_xvel,[],rtx_t,'.')
colormap jet
xlabel('Time')
ylabel('X Speed (m/s)')
datetick('x','HH:MM')
grid on
axis tight

subplot(3,1,2)
scatter(rtx_t,posmv_yvel,[],rtx_t,'.')
colormap jet
xlabel('Time')
ylabel('Y Speed (m/s)')
datetick('x','HH:MM')
grid on
axis tight

subplot(3,1,3)
scatter(rtx_t,posmv_zvel,[],rtx_t,'.')
colormap jet
xlabel('Time')
ylabel('Z Speed (m/s)')
datetick('x','HH:MM')
grid on
axis tight
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% scatter velocity vs offset
%{
figure(56)
scatter(posmv_xvel,offset_rtx_posmv_tx_s3(1,:),[],rtx_t,'.')
ylim([-15 5])
grid on
colormap jet
xlabel('Longitudinal Speed (m/s)')
ylabel('Offset in Bow-Stern direction (m)')
colorbar
cbdate()

figure(76)
scatter(posmv_yvel,offset_rtx_posmv_tx_s3(2,:),[],rtx_t,'.')
ylim([-10 10])
grid on
colormap jet
xlabel('Transverse Speed (m/s)')
ylabel('Offset in Starboard-Port direction (m)')
colorbar
cbdate()

figure(86)
scatter(posmv_zvel,offset_rtx_posmv_tx_s3(3,:),[],rtx_t,'.')
ylim([-15 15])
grid on
colormap jet
xlabel('Down Speed (m/s)')
ylabel('Offset in vertical direction (m)')
colorbar
cbdate()
%}
%% scatter velocity in X,Y and horizontal offset
%{
figure(66)
clf
scatter(posmv_xvel,posmv_yvel,[],offset_rtx_posmv_tx_s3(2,:),'.')
grid on
c = colorbar;
caxis([-5 5])
xlabel('Longitudinal Speed (m/s)')
ylabel('Transvers Speed (m/s)')
c.Label.String = 'Port-Starboard offset (m)';
title('Port-Starboard Offset vs Velocities in X and Y')
colormap(flipud(jet))
xlim([-7 7])
ylim([-7 7])

%}
%% RTX and POSMV transducer points in the UTM coordinate system
%{
f = figure(5);
f.Units = 'normalized';
f.Position = [0.3 0.5 0.5 0.55];
clf
for l = 40:50%length(rtx_x_ant_gb)
    drawnow
    grid on
    xlabel('Easting (m)')
    ylabel('Northing (m)')
    scatter(rtx_x_ant_gb(l),rtx_y_ant_gb(l),'kx')
    hold on
    scatter(rtx_x_tx_gb(l),rtx_y_tx_gb(l),'r.')
    scatter(posmv_x_gb(l),posmv_y_gb(l),'b.')
    legend('RTX antenna','RTX transducer','POS MV transducer','Location','northwest'	)
end
%}
%% Vertical
%{
f = figure(7);
f.Units = 'normalized';
f.Position = [0.3 0.5 0.5 0.55];
clf
for l = 10:50%length(rtx_x_ant_gb)
    drawnow
    grid on
    xlabel('Easting (m)')
    ylabel('Northing (m)')
    scatter(rtx_t(l),rtx_z_ant_gb(l),'kx')
    hold on
    scatter(rtx_t(l),rtx_tx_altitude(l),'r.')
    scatter(posmv_t(l),posmv_altitude(l),'b.')
end
%}
%% verify transformed coordinate
%{
f = figure(777);
f.Units = 'normalized';
f.Position = [0.3 0.5 0.4 0.6];
clf
subplot(3,1,1)
plot(rtx_t,offset_rtx_tx_ant_s3(1,:) )
datetick('x')
grid on
ylabel('Bow-Stern (m)')
headline = sprintf('Median = %.4f m',median(offset_rtx_tx_ant_s3(1,:)));
title(headline)
ylim([-5.8725 -5.8695])

subplot(3,1,2)
plot(rtx_t,offset_rtx_tx_ant_s3(2,:) )
datetick('x')
grid on
ylabel('Port-Starboard (m)')
headline = sprintf('Median = %.4f m',median(offset_rtx_tx_ant_s3(2,:)));
title(headline)
ylim([13.092 13.096])

subplot(3,1,3)
plot(rtx_t,offset_rtx_tx_ant_s3(3,:) )
datetick('x')
grid on
ylabel('Vertical (m)')
headline = sprintf('Median = %.4f m',median(offset_rtx_tx_ant_s3(3,:)));
title(headline)
ylim([28.475 28.477])
%}
%% offset RTX POSMV transducer
med_xoff = median(offset_rtx_posmv_tx_s3(1,indstart:indend))
med_yoff = median(offset_rtx_posmv_tx_s3(2,indstart:indend))
med_zoff = median(offset_rtx_posmv_tx_s3(3,indstart:indend))
rmsx = offset_rtx_posmv_tx_s3(1,indstart:indend) - med_xoff;
rmsy = offset_rtx_posmv_tx_s3(2,indstart:indend) - med_yoff;
rmsz = offset_rtx_posmv_tx_s3(3,indstart:indend) - med_zoff;

figure(23)
subplot(3,1,1)
scatter(posmv_t(indstart:indend),offset_rtx_posmv_tx_s3(1,indstart:indend),10,posmv_xvel(indstart:indend),'filled')
cbar = colorbar;
cbar.Label.String = 'X velocity (m/s)';
colormap jet;
grid on
datetick('x')
ylabel('Bow-Stern Offset (m)')
ylim([-1 1])
yticks([-1:.2:1])
title('RTX - POSMV (bin) transducer offset ')

subplot(3,1,2)
scatter(posmv_t(indstart:indend),offset_rtx_posmv_tx_s3(2,indstart:indend),10,posmv_yvel(indstart:indend),'filled')
cbar = colorbar;
cbar.Label.String = 'Y velocity (m/s)';
colormap jet;
grid on
datetick('x')
ylabel('Port-Starboard Offset (m)')
ylim([-1 1])
yticks([-1:.2:1])

subplot(3,1,3)
scatter(posmv_t(indstart:indend),offset_rtx_posmv_tx_s3(3,indstart:indend),10,posmv_zvel(indstart:indend),'filled')
cbar = colorbar;
cbar.Label.String = 'Z velocity (m/s)';
colormap jet;
grid on
datetick('x')
ylabel('Vertical Offset (m)')
ylim([0 4])
yticks([-1:.2:1])
%%






%% Cal Lag between RTX and POSMV
% EDIT manually pick time frame of interest
hour_start =  '10:10:00';              % start time
date = '10-27-2018';                    % date

duration_sec = 60/(3600*24);           % sub-window length
t_start = datenum([date hour_start],'mm-dd-yyyyHH:MM:SS');
avg_min = 5;                           % edit window length
avg_duration = avg_min/(60*24);              % window length
t_step = 1/(60*24);                     % time step
t_plot = [];
avg_xspeed = zeros(1,avg_min);
avg_xl = zeros(1,avg_min);
avg_yl = zeros(1,avg_min);
cal_offset = zeros(1,avg_min);
avg_xoffset = zeros(1,avg_min);
x_speed = [];
x_offset = [];
for p = 0:avg_min
    t_start = datenum([date hour_start],'mm-dd-yyyyHH:MM:SS')+ t_step*p;
    t_end = t_start + duration_sec;
    datestr(t_start)
    ind_rtx = find(t_start-0.06/(3600*24) <= rtx_t & rtx_t<= t_end+6.06/(3600*24));
    ind_posmv = find(t_start-0.06/(3600*24) <= posmv_t & posmv_t <= t_end+0.06/(3600*24));
    ind_posmv10Hz = find(t_start <= posmv_t10hz & posmv_t10hz <= t_end);

    % UTM coordinate transformation
    %[posmv_x_gb10hz,posmv_y_gb10hz,utmzone] = deg2utm(posmv_lat10hz,posmv_lon10hz);
    
    % velocity
    posmv_xvelw = posmv_xvel(ind_posmv);
    avg_xspeed(p+1) = mean(posmv_xvelw);
    x_speed = horzcat(x_speed,posmv_xvelw);
    % Heading
    yaww =  posmv_heading(ind_posmv)/180*pi;
    
   
    rtx_x_tx_gbw = rtx_x_tx_gb(ind_rtx);
    rtx_y_tx_gbw = rtx_y_tx_gb(ind_rtx);
    rtx_z_tx_gbw = rtx_tx_altitude(ind_rtx);

    posmv_x_gbw = posmv_x_gb(ind_posmv);
    posmv_y_gbw = posmv_y_gb(ind_posmv);
    posmv_z_gbw = posmv_altitude10hz(ind_posmv10Hz);

    % Calculate Lags by least sum squared- difference
    [lsd_x,x_shift]= lag_rtx_posmv(rtx_x_tx_gbw,posmv_x_gbw);
    [lsd_y,y_shift]= lag_rtx_posmv(rtx_y_tx_gbw,posmv_y_gbw);
    
    avg_xl(p+1) = x_shift;
    avg_yl(p+1) = y_shift;
    
    % for plot
    t_plot= horzcat(t_plot,posmv_t(ind_posmv));
    x_offset = horzcat(x_offset,offset_rtx_posmv_tx_s3(1,ind_rtx));
    % calculated offset
    cal_offset(p+1) =  mean(posmv_xvelw.*sqrt(cos(yaww).^2*y_shift^2+sin(yaww).^2*x_shift^2));
     % avg offset
    avg_xoffset(p+1) = mean(offset_rtx_posmv_tx_s3(1,ind_rtx));
    

end

fprintf(sprintf('Average Longitudinal Speed = %.3f m/s\nCalculated X Offset in ship coordinate = %.3f m \nAveraged X Offset = %.3f m\n',mean(avg_xspeed),-mean(cal_offset),mean(avg_xoffset)))


% Plot lags
%%%%%%%%%%%%%%%%%%
% X
figure(99)
clf
subplot(2,1,1)
scatter(posmv_t(ind_posmv),posmv_x_gbw,'.');
hold on
scatter(rtx_t(ind_rtx),rtx_x_tx_gbw,'.')
datetick('x')
grid on
legend('POS MV', 'RTX')
ylabel('Easting (m)')
title('Original timeseries')

subplot(2,1,2)
scatter(posmv_t(ind_posmv)+(x_shift)/(3600*24),posmv_x_gbw,'.');
hold on
scatter(rtx_t(ind_rtx),rtx_x_tx_gbw,'.')
datetick('x')
grid on
ylabel('Easting (m)')
%%%%%%%%%%%%%%%
% Y
figure(109)
clf
subplot(2,1,1)
scatter(posmv_t(ind_posmv),posmv_y_gbw,'.');
hold on
scatter(rtx_t(ind_rtx),rtx_y_tx_gbw,'.')
datetick('x')
grid on
legend('POS MV', 'RTX')
ylabel('Northing (m)')
title('Original timeseries')

subplot(2,1,2)
scatter(posmv_t(ind_posmv)+(y_shift)/(3600*24),posmv_y_gbw,'.');
hold on
scatter(rtx_t(ind_rtx),rtx_y_tx_gbw,'.')
datetick('x')
grid on
legend('POS MV', 'RTX')
ylabel('Northing (m)')
title('Shifted timeseries')

figure(129)
plot(posmv_t(ind_posmv),posmv_xvelw)
grid on
datetick('x')
headline = sprintf('Mean = %.3f',mean(posmv_xvelw));
title(headline)
ylabel('Ship Speed (m/s)')
%}



