% Find transmission time and transducer position based on POSMV GGA and interpolation
% 1. RTX file truncation based on Time of Day and add transducer positions
% 2. extract transmission times from tx_file, match RTX data points 
% 3. find corresponding POSMV GGA transducer position
% 4. read in 10 Hz POS MV data
% 5. find 2 nearest 10 Hz POS MV data to extract roll/pitch/yaw/velocity x z
% 6. Use time differences between the POSMV GGA data point and those 2 POS MV
% data points to interpolate for the antenna position
% 7. convert back to lat/lon altitude
clear
close all
%% Download RTX file

load posmvGGA_1819
gga_lat = posmv_gga.lat;
gga_lon = posmv_gga.lon;
gga_altitude = posmv_gga.altitude;
gga_t = posmv_gga.t;
gga_x_err = posmv_gga.lat_err;
gga_y_err = posmv_gga.lon_er;
gga_z_err = posmv_gga.altitude_err;


%% Download POS MV files
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/tx_file/all
d = dir('tx*');

%% Loop over all files and read the tx times
% find nearest corresponding rtx times
% calculate transducer position
% pack data and save

% outputs 
% rtx_lat_nc = [];
% rtx_lon_nc = [];

x_err= [];
y_err = [];
z_err = [];
tx_lat_tot = [];
tx_lon_tot = [];
tx_altitude_tot = [];

% stat
x_offset_tot = [];
y_offset_tot = [];
z_offset_tot = []; 
t_offset_tot = [];  % offset between chosen RTX time and transmission time
xvel_tot = [];
zvel_tot = [];
yaw_tot = [];
 
t_diff = [];        % offset between all RTX times and transmission times
t_1sec = 0;
gga_t_f = [];
tx_t_f  = [];
for ii = 1:length(d)
    cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/tx_file/all
    % extract tx data
     fname = d(ii).name;
    load (fname)
    posmv_tx_roll = tx_data.roll;
    posmv_tx_pitch = tx_data.pitch;
    posmv_tx_yaw = tx_data.heading;
    posmv_heave = tx_data.heave;
    posmv_tx_t = tx_data.t;
    
    clearvars tx_data
    ind_pos = [];
    rm_ind_t = [];
    
    tx_lat_gga = [];
    tx_lon_gga = [];
    tx_altitude_gga = [];
    
    t_offset_hrly = [];
    xvel_hrly = [];
    zvel_hrly = [];
     heading_hrly = [];
    
    % match transmission times with gga data point
    for k = 1:length(posmv_tx_t)
        % find corresponding gga data points to use its antenna position
        [t_difff,I] = min(abs(posmv_tx_t(k) - gga_t));
        t_diff(end+1) = t_difff*3600*24;
        if t_difff*3600*24 > 0.5          % discard transmissions with a gap larger than 0.5 sec
            rm_ind_t(end+1) = k;
            t_1sec = t_1sec+1;
        else
            ind_pos(end+1) = I;     % index to extract lon/lat/alt from the gga
        end
    end
    gga_t_w = posmv_tx_t;       % chosen tx times
    gga_t_w(rm_ind_t) = [];
    
    % antenna position
    gga_t_pick = gga_t(ind_pos);    % gga timestamps
    
    % lat lon altitude fro gga
    gga_lat_w = gga_lat(ind_pos);
    gga_lon_w = gga_lon(ind_pos);
    gga_altitude_w = gga_altitude(ind_pos);
    gga_x_err_w = gga_x_err(ind_pos);
    gga_y_err_w = gga_y_err(ind_pos);
    gga_z_err_w = gga_z_err(ind_pos);
     
     % time diff anal
%      gga_t_f = horzcat(gga_t_f,gga_t_pick);
%      tx_t_f = horzcat(tx_t_f,gga_t_w);
%      
     
     % 4. Read in 10Hz POSZ MV binary file
     % search the file by hour
     
     posmvbin_t = [];
     posmvbin_xvel = [];
     posmvbin_yvel = [];
     posmvbin_zvel = [];
     posmvbin_roll = [];
     posmvbin_pitch = [];
     posmvbin_heading = [];
     
     cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/posmv/all
     hour = str2num(fname(end-5:end-4));
     day =  str2num(fname(end-8:end-7));
     % load 2 posmv files
     for p = 0:1
        hour = hour+p;
         if hour <= 23
          if hour < 10
              shour = "0"+string(hour);
              sday = string(day);
          else
              shour = string(hour);
              sday = string(day);
          end
         elseif hour>=24
             shour = "00"
             sday = string(day+1); 
         end
         
     posmv_10Hz_fname = "posmv_201810"+ sday + shour;
     load(posmv_10Hz_fname)
     posmv_t_d = posmv.t;
     posmv_xvel_d = posmv.x_vel;
     posmv_yvel_d = posmv.y_vel;
     posmv_zvel_d = posmv.z_vel;
     posmv_roll_d = posmv.roll;
     posmv_pitch_d = posmv.pitch;
     posmv_yaw_d = posmv.heading;
     
     % eliminate outlier
     posmv_t_d(find(posmv_t_d < posmv_t_d(1))) = [];
     posmv_xvel_d(find(posmv_t_d < posmv_t_d(1))) = [];
     posmv_yvel_d(find(posmv_t_d < posmv_t_d(1))) = [];
     posmv_zvel_d(find(posmv_t_d < posmv_t_d(1))) = [];
     posmv_roll_d(find(posmv_t_d < posmv_t_d(1))) = [];
     posmv_pitch_d(find(posmv_t_d < posmv_t_d(1))) = [];
     posmv_yaw_d(find(posmv_t_d < posmv_t_d(1))) = [];
     
     % concatenate data
     posmvbin_t = horzcat(posmvbin_t,posmv_t_d');
     posmvbin_xvel = horzcat(posmvbin_xvel,posmv_xvel_d');
     posmvbin_yvel = horzcat(posmvbin_yvel,posmv_yvel_d');
     posmvbin_zvel = horzcat(posmvbin_zvel,posmv_zvel_d');
     posmvbin_roll = horzcat(posmvbin_roll,posmv_roll_d');
     posmvbin_pitch = horzcat(posmvbin_pitch,posmv_pitch_d');
     posmvbin_heading = horzcat(posmvbin_heading,posmv_yaw_d');
     
     end
     
     % cut off out-of-frame data points
     len = length(gga_t_w(gga_t_w<posmvbin_t(end)));
     
     % getting lat/lon/alt attitudes
     % loop over tx times
     % interpolation
     for r = 1:len
         
         % time offset of the actually gga times from the transmission times
         t_offset = 3600*24*(gga_t_pick(r)  - gga_t_w(r));
         
         % find attitudes (POSMV / Tx file)
         ind_lower_posmv_att = find((posmvbin_t - gga_t_w(r))<0);
         ind_lower_posmv_att = ind_lower_posmv_att(end);
         ind_upper_posmv_att = find((posmvbin_t - gga_t_w(r))>0);
         ind_upper_posmv_att = ind_upper_posmv_att(1);
         t_lower = posmvbin_t(ind_lower_posmv_att);
         t_upper = posmvbin_t(ind_upper_posmv_att);
         t_offset_flower = (gga_t_w(r) - t_lower)*3600*24;
         
             % find roll pitch yaw
             lower_roll =  posmvbin_roll(ind_lower_posmv_att);
             lower_pitch =  posmvbin_pitch(ind_lower_posmv_att);
             lower_yaw =  posmvbin_heading(ind_lower_posmv_att);
             upper_roll =  posmvbin_roll(ind_upper_posmv_att);
             upper_pitch =  posmvbin_pitch(ind_upper_posmv_att);
             upper_yaw =  posmvbin_heading(ind_upper_posmv_att);
        
        %instantenous attitudes
        roll_now = (upper_roll - lower_roll)*(t_offset_flower/((t_upper - t_lower)*3600*24))+lower_roll;
        pitch_now = (upper_pitch - lower_pitch)*(t_offset_flower/((t_upper - t_lower)*3600*24))+lower_pitch;
        yaw_now = (upper_yaw - lower_yaw)*(t_offset_flower/((t_upper - t_lower)*3600*24))+lower_yaw;
        
        % extract x/z vel 
        if t_offset > 0
            ind_upper2 = min(find((posmvbin_t - gga_t_pick(r))>0));
            ind_lower2 = max(find((posmvbin_t - gga_t_w(r))<0));
        elseif t_offset < 0 
            ind_upper2 = min(find((posmvbin_t - gga_t_w(r))>0));
            ind_lower2 = max(find((posmvbin_t - gga_t_pick(r))<0));
        end
        
        t_upper2 =  posmvbin_t(ind_upper2);
        t_lower2 = posmvbin_t(ind_lower2);
        
        xvel_upper2 = posmvbin_xvel(ind_upper2);
        xvel_lower2 = posmvbin_xvel(ind_lower2);
        yvel_upper2 = posmvbin_yvel(ind_upper2);
        yvel_lower2 = posmvbin_yvel(ind_lower2);
        zvel_upper2 = posmvbin_zvel(ind_upper2);
        zvel_lower2 = posmvbin_zvel(ind_lower2);
        
        mean_xvel =  mean([xvel_upper2 xvel_lower2]);
        mean_yvel =  mean([yvel_upper2 yvel_lower2]);
        mean_zvel = mean([zvel_upper2 zvel_lower2]);
        
        % x velocity
        yaw_tot = horzcat(yaw_tot,yaw_now);
        
        Vx = mean_xvel*sin(yaw_now/180*pi);
        Vy = mean_xvel*cos(yaw_now/180*pi);
        X_offset = -Vx*(t_offset);      % meter
        Y_offset = -Vy*(t_offset);      % meter
        Z_offset = -mean_zvel*(t_offset);
        
        % To UTM of the antenna
%         gga_lat_nc = horzcat(gga_lat_nc ,gga_lat_w(r));
%         gga_lon_nc = horzcat(gga_lon_nc,gga_lon_w(r) );

        [td_x,td_y,utmzone] =  deg2utm(gga_lat_w(r),gga_lon_w(r));
        
        % offset the antenna to the time of transmission
        td_x2 = td_x+X_offset;
        td_y2 = td_y+Y_offset;
        td_z2 = gga_altitude_w(r)+Z_offset;
        
        % To degree
        [gga_lat_now,gga_lon_now] = utm2deg(td_x2,td_y2,utmzone);
%         gga_lat_td(end+1) = gga_lat_now;
%         gga_lon_td(end+1) = gga_lon_now;
        
        % calculate transducer position (global coordinate) of the actual gga time
        
        [td_lat,td_lon,td_altitude] = transducer_pos_rtx(gga_lat_now,gga_lon_now,td_z2,1,roll_now,pitch_now,yaw_now,1);
        
        
        tx_lat_gga = horzcat(tx_lat_gga,td_lat);
        tx_lon_gga = horzcat(tx_lon_gga,td_lon);
        tx_altitude_gga = horzcat(tx_altitude_gga,td_altitude);
        
        % stat
        x_offset_tot = horzcat(x_offset_tot ,X_offset);
        y_offset_tot = horzcat(y_offset_tot ,Y_offset);
        z_offset_tot = horzcat(z_offset_tot ,Z_offset); 
        t_offset_tot = horzcat(t_offset_tot ,t_offset);
        xvel_tot = horzcat(xvel_tot ,mean_xvel);
        zvel_tot = horzcat(zvel_tot ,mean_zvel);
        
        
        heading_hrly = horzcat(heading_hrly,yaw_now);
        xvel_hrly = horzcat(xvel_hrly,mean_xvel);
        zvel_hrly = horzcat(zvel_hrly,mean_zvel);
        t_offset_hrly = horzcat(t_offset_hrly ,t_offset);
    
     
     end 
     
     
     
   
     % 7. pack data and save
    tx_lat_tot = horzcat(tx_lat_tot,tx_lat_gga) ;
    tx_lon_tot = horzcat(tx_lon_tot,tx_lon_gga) ;
    tx_altitude_tot = horzcat(tx_altitude_tot,tx_altitude_gga) ;

%     size(gga_t_w)
%     size(tx_lat_gga)
    
    tx_data.t = gga_t_w;
    tx_data.lat = tx_lat_gga;
    tx_data.lon = tx_lon_gga;
    tx_data.altitude = tx_altitude_gga; 
    tx_data.t_offset = t_offset;
    tx_data.x_err = gga_x_err_w ; 
    tx_data.y_err = gga_y_err_w ; 
    tx_data.z_err =  gga_z_err_w; 
%     gga_tx_data.xvel = xvel_hrly;
%     gga_tx_data.zvel = zvel_hrly;
%     gga_tx_data.heading = heading_hrly;
    cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/POSMV_GGA/tx_file/all
    savename = sprintf("posmvGGA_"+string(fname(1:end-4))+".mat")
    save(savename,'tx_data');
 
end
% %% Plot stat
% figure(1)
% subplot(4,1,1)
% plot(tx_t_f,x_offset_tot)
% grid on
% axis tight
% ylabel('X Offset (m)')
% datetick('x')
% 
% subplot(4,1,2)
% plot(tx_t_f,y_offset_tot)
% grid on
% axis tight
% ylabel('Y Offset (m)')
% datetick('x')
% 
% subplot(4,1,3)
% plot(tx_t_f,t_offset_tot)
% grid on
% axis tight
% ylabel('Time Offset (s)')
% datetick('x')
% 
% subplot(4,1,4)
% plot(xvel_tot)
% grid on
% axis tight
% ylabel('X Velocity (m/s)')
datetick('x')
% %% Plot coordinates
% figure(2)
% clf
% index = 3677:3681;
% scatter(gga_lon_nc(index),gga_lat_nc(index),[],xvel_tot(index),'filled')
% hold on
% scatter(gga_lon_c(index),gga_lat_c(index),[],xvel_tot(index),'x')
% hold on 
% scatter(tx_lon_tot(index),tx_lat_tot(index),[],xvel_tot(index),'o')
% colormap jet
% grid on
% colorbar
% 
% figure(3)
% subplot(3,1,1)
% plot(tx_t_f(index),x_offset_tot(index))
% grid on
% datetick('x','dd HH:MM')
% axis tight
% 
% subplot(3,1,2)
% plot(tx_t_f(index),y_offset_tot(index))
% grid on
% datetick('x')
% axis tight
% 
% subplot(3,1,3)
% plot(tx_t_f(index),yaw_tot(index))
% grid on
% datetick('x')
% axis tight