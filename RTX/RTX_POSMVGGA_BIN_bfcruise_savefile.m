% compare RTX and POSMV GGA ans POSMV binary antenna positions
% transform to the ship position
% save the output file
close all
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load RTX and POSMV GGA data
cd /Users/testuser/Documents/MATLAB/Script/Data
load('rtx10Hz_posmvGGA_beforecruise.mat')
rtx_t = rtx_t - 35/(3600*24);        % adjust rtx time
[rtx_t,~,~] = unique(rtx_t,'first');

%% main Loop
% create space for antenna position vectors
% RTX
rtx_t_all = [];
rtx_x_ant_gb = [];
rtx_y_ant_gb = [];
rtx_z_ant_gb = [];

% POSMV GGA
posmv_t_all = [];
posmv_x_ant_gb = [];
posmv_y_ant_gb = [];
posmv_z_ant_gb =[];

% POS MV binary output
% Sensor 2
posmvbin_x_ant_gb = [];
posmvbin_y_ant_gb = [];
posmvbin_z_ant_gb =[];

% sensor 1
posmvbin_x_ant_gb2 = [];
posmvbin_y_ant_gb2 = [];
posmvbin_z_ant_gb2 = [];

posmv_td_all = [];
posmv_gn_all = [];


% % raw data of sensor 1 and 2
% % Sensor 2
% posmvbin_x_td_gb = [];
% posmvbin_y_td_gb = [];
% posmvbin_z_td_gb =[];
% 
% % sensor 1
% posmvbin_x_gn_gb2 = [];
% posmvbin_y_gn_gb2 = [];
% posmvbin_z_gn_gb2 = [];

% parameters
posmv_heading = [];
posmv_xvel = [];
posmv_yvel = [];
posmv_zvel = [];
posmv_xacc = [];
posmv_yacc = [];
posmv_zacc = [];

% Interpolation residual vectors
rtx_X_dis = [];
rtx_Y_dis = [];
rtx_Z_dis = [];
rtx_toffset = [];

gga_X_dis = [];
gga_Y_dis = [];
gga_Z_dis = [];
gga_toffset = [];

% RTX/GGA/Bin antenna offset from POSMV's granite block (ship locally level frame)
offset_rtx_ref_gn_s3 = [];          % RTX
offset_posmvgga_ref_gn_s3 = [];     % GGA
offset_posmvbins2_ref_gn_s3 = [];   % sensor 2
offset_posmvbins1_ref_gn_s3 = [];   % sensor 1

posmvbins2_ref_gn_s3 = [];          % sensor 2 (transducer)

posmvbins2_ref_rtx_s3 = [];          % sensor 2 (transducer)
posmvbins1_ref_rtx_s3 = [];           % sensor 1 (granite block)
posmvbins2_ref_gga_s3 = [];
posmvbins1_ref_gga_s3 = [];

% coordinates of devices with respect to the granite block
ant = [6.439;6.505;-27.761]; % from KM coordinate system May 2018 spreasheet
tx5 = [0.568;19.599;0.715];
ant2tx = tx5-ant;            % vector pointing from the transducer to the antenna

% Rotation matrix: Global, UTM coordiante to unrotated ship coordinate
thetax_gb = 180/180*pi;
thetaz_gb = -90/180*pi;
Rxgb = [1 0 0;0 cos(thetax_gb) -sin(thetax_gb);0 sin(thetax_gb) cos(thetax_gb)];
Rzgb = [cos(thetaz_gb) -sin(thetaz_gb) 0 ; sin(thetaz_gb) cos(thetaz_gb) 0; 0 0 1];
Rgb_s = Rxgb*Rzgb;
Rs_gb = transpose(Rgb_s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop over posmv binary files
% EDIT start file and end file

    
    cd /Users/testuser/Documents/MATLAB/Script/Data
    % Load POSMV Binary file
    % transducer
    
    load ('posmv_20181024.mat')
    posmv_latd = posmv.lat;
    posmv_lond = posmv.lon;
    posmv_altituded =posmv.altitude;
    posmv_td = posmv.t;
    
    % vessel dynamics
    posmv_xveld = posmv.x_vel;
    posmv_yveld = posmv.y_vel;
    posmv_zveld = posmv.z_vel;
    posmv_rolld = posmv.roll;
    posmv_pitchd = posmv.pitch;
    posmv_headingd = posmv.heading;
    
    % load granite block
    posmv_lat_gnd = posmv.lat_s;
    posmv_lon_gnd = posmv.lon_s;
    posmv_altitude_gnd =posmv.altitude_s;
    posmv_t_gnd = posmv.t;
    
    ind_l = [];
    % eliminate POSMV binary outliers
    for l = 2:length(posmv_td)
        if (posmv_td(l) - posmv_td(l-1)) < 0 
            ind_l = l;
            break;
        end
    end
    if ~isempty(ind_l)
        posmv_latd(ind_l) = [];
        posmv_lond(ind_l) = [];
        posmv_altituded(ind_l) = [];
        
        posmv_rolld(ind_l) = [];
        posmv_pitchd(ind_l) = [];
        posmv_headingd(ind_l) = [];
        posmv_td(ind_l) = [];
        posmv_xveld(ind_l) = [];
        posmv_yveld(ind_l) = [];
        posmv_zveld(ind_l) = [];
        
  
    end
    
    %%%%%% Truncate RTX and GGA data to the duration within POSMV binary file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rtx_t_w = [];       % rtx times
    gga_t_w = [];       % POSMV GGA times
    ind_1sec_rtx = [];   % eliminated indicies rtx
    ind_1sec_gga = [];   % eliminated indicies GGA
    
    % window RTX/GGA timeseries to a specific interval
    start_t = posmv_td(find(posmv_td >= datenum('20181024 00:00:00','yyyymmdd HH:MM:SS'),1,'first'));
    end_t = start_t+20/24;
    [~,start_ind] = min(abs(rtx_t-start_t));
    [~,end_ind] = min(abs(rtx_t-end_t));
    
    [~,start_ind2] = min(abs(posmvgga_t-start_t));
    [~,end_ind2] = min(abs(posmvgga_t-end_t));
    rtx_t_w = rtx_t(start_ind:end_ind);
    gga_t_w = posmvgga_t(start_ind2:end_ind2);
    
    % window antenna positions
    rtx_lat_w = rtx_lat(start_ind:end_ind);
    rtx_lon_w = rtx_lon(start_ind:end_ind);
    rtx_altitude_w = rtx_altitude(start_ind:end_ind);
    
    gga_lat_w = posmvgga_lat(start_ind2:end_ind2);
    gga_lon_w = posmvgga_lon(start_ind2:end_ind2);
    gga_altitude_w = posmvgga_altitude(start_ind2:end_ind2);
    
     % match RTX with POS MV GGA
    rm_ind_rtx = [];
    rm_ind_gga = [];
    for k = 1:length(rtx_t_w)
       [t_difff,I] =  min(abs(rtx_t_w(k) - gga_t_w));
       if t_difff*3600*24 >= 0.5    % if time difference is larger than 0.5 sec
            rm_ind_rtx(end+1) = k;         % remove data point                    
       end
    end
    rtx_t_w(rm_ind_rtx) = [];
    rtx_lat_w(rm_ind_rtx) = [];
    rtx_lon_w(rm_ind_rtx) = [];
    rtx_altitude_w(rm_ind_rtx) = [];
    
    
    for k = 1:length(gga_t_w)
       [t_difff,I] =  min(abs(gga_t_w(k) - rtx_t_w));
       if t_difff*3600*24 >= 0.5    % if time difference is larger than 0.5 sec
            rm_ind_gga(end+1) = k;         % extract ship attitudes from POS-MV                    
       end
    end
    gga_t_w(rm_ind_gga) = [];
    gga_lat_w(rm_ind_gga) = [];
    gga_lon_w(rm_ind_gga) = [];
    gga_altitude_w(rm_ind_gga) = [];
   
    % Match POSMV binary data points
    % find the nearest neighbor POSMV binary data points of the corresponding RTX data points
    ind = [];           % posmv bin indices  to keep data
    for k = 1:length(rtx_t_w)
       [t_difff,I] =  min(abs(rtx_t_w(k) - posmv_td));
       if t_difff*3600*24 <= 0.5    % if time difference is larger than 0.5 sec
            ind(end+1) = I;         % extract ship attitudes from POS-MV 
            %t_diff(end+1) = t_difff*3600*24;
       else
           ind_1sec_rtx(end+1) = k;   % discard a point whose nearest neighbor POSMV data point is farther than 0.5 sec
           % t_1sec = t_1sec+1;
        end
    end
    
    % eliminate RTX data points without POS MV match
    if length(ind_1sec_rtx) > 0
        rtx_t_w(ind_1sec_rtx) = [];       
        rtx_lat_w(ind_1sec_rtx) = [];
        rtx_lon_w(ind_1sec_rtx) = [];
        rtx_altitude_w (ind_1sec_rtx) = [];
        gga_t_w(ind_1sec_rtx) = [];       
        gga_lat_w(ind_1sec_rtx) = [];
        gga_lon_w(ind_1sec_rtx) = [];
        gga_altitude_w(ind_1sec_rtx) = [];
    
    end
    
    % attitudes for coordinate transformation (at POS MV bin time)
    rtx_roll_w = posmv_rolld(ind);
    rtx_pitch_w = posmv_pitchd(ind);
    rtx_heading_w = posmv_headingd(ind);
         
    % Interpolate RTX and GGA position to the nearest POSMV bin data point
%     [rtx_lat_w,rtx_lon_w,rtx_altitude_w,rtxx,rtxy,rtxz,rtxt] = position_interpolation(rtx_lat_w,rtx_lon_w,rtx_altitude_w,rtx_t_w,posmv_xveld,posmv_yveld,posmv_zveld,posmv_rolld,posmv_pitchd,posmv_headingd,posmv_td,posmv_td(ind));
%     [gga_lat_w,gga_lon_w,gga_altitude_w,ggax,ggay,ggaz,ggat] = position_interpolation(gga_lat_w,gga_lon_w,gga_altitude_w,gga_t_w,posmv_xveld,posmv_yveld,posmv_zveld,posmv_rolld,posmv_pitchd,posmv_headingd,posmv_td,posmv_td(ind));
    rtxx = [];
    rtxy = [];
    rtxz = [];
    rtxt = [];
    ggax = [];
    ggay = [];
    ggaz = [];
    ggat = [];


    rtx_X_dis = vertcat(rtx_X_dis,rtxx);
    rtx_Y_dis = vertcat(rtx_Y_dis,rtxy);
    rtx_Z_dis = vertcat(rtx_Z_dis,rtxz);
    rtx_toffset = vertcat(rtx_toffset,rtxt);
    
    gga_X_dis = vertcat(gga_X_dis,ggax);
    gga_Y_dis = vertcat(gga_Y_dis,ggay);
    gga_Z_dis = vertcat(gga_Z_dis,ggaz);
    gga_toffset  = vertcat(gga_toffset,ggat);
    
    % calculate POSMV Bin antenna position from Sensor 2 (global coordinate)
    [pos_ant_lat,pos_ant_lon,pos_ant_altitude] = ant_pos_posmv(posmv_latd(ind),posmv_lond(ind),posmv_altituded(ind),posmv_td(ind),posmv_rolld(ind),posmv_pitchd(ind),posmv_headingd(ind),2);
    % calculate POSMV Bin antenna position from Sensor 1 (global coordinate)
    [pos_ant_lat2,pos_ant_lon2,pos_ant_altitude2] = ant_pos_posmv(posmv_lat_gnd(ind),posmv_lon_gnd(ind),posmv_altitude_gnd(ind),posmv_t_gnd(ind),posmv_rolld(ind),posmv_pitchd(ind),posmv_headingd(ind),1);
    
    % convert lat lon to UTM x,y,z  transducer/ granite block POSMV Bin
    [posmv_x_td,posmv_y_td,utmzone] = deg2utm(posmv_latd(ind),posmv_lond(ind));
    posmv_z_td = posmv_altituded(ind);
    posmv_td_gb = [posmv_x_td'; posmv_y_td' ; posmv_z_td'];  
    
    [posmv_x_gnd,posmv_y_gnd,utmzone] = deg2utm(posmv_lat_gnd(ind),posmv_lon_gnd(ind));
    posmv_z_gnd = posmv_altitude_gnd(ind);
    posmv_gn_gb = [posmv_x_gnd'; posmv_y_gnd' ; posmv_z_gnd'];  
    
    % RTX/GGA position
    [rtx_x_ant,rtx_y_ant,~] = deg2utm(rtx_lat_w,rtx_lon_w);
     rtx_ant_gb = [rtx_x_ant';rtx_y_ant'; rtx_altitude_w'];
     
     [gga_x_ant,gga_y_ant,utmzone] = deg2utm(gga_lat_w,gga_lon_w);
     gga_ant_gb = [gga_x_ant';gga_y_ant'; gga_altitude_w'];
     
     % Sensor 2 antenna
    [posmvbin_x_d,posmvbin_y_d,utmzone] = deg2utm(pos_ant_lat,pos_ant_lon);
    posmvbin_z_d = pos_ant_altitude;
    posmvbin_ant_gb = [posmvbin_x_d'; posmvbin_y_d' ; posmvbin_z_d'];
    % Sensor 1 antenna
    [posmvbin_x_d2,posmvbin_y_d2,utmzone] = deg2utm(pos_ant_lat2,pos_ant_lon2);
    posmvbin_z_d2 = pos_ant_altitude2;
    posmvbin_ant_gb2 = [posmvbin_x_d2'; posmvbin_y_d2' ; posmvbin_z_d2'];
     
    % transform antenna positions the global UTM coordinate back to ship locally level frame
    rtx_s3_refgn = [];
    gga_s3_refgn = [];
    bin_s3_refgn = [];
    bin2_s3_refgn = [];
    % transducer
    td_s3_refgn = [];
    % granite block and transducer
    gn_s3_refrtx = [];
    td_s3_refrtx = [];
    gn_s3_refgga = [];
    td_s3_refgga = [];
    
    rtx_roll_w2 = rtx_roll_w/180*pi; 
    rtx_pitch_w2 = rtx_pitch_w/180*pi; 
    rtx_heading_w2 = rtx_heading_w/180*pi; 
    
    % Perform coordinate transformation from UTM to ship locally level frame
    for k = 1:length(rtx_roll_w)
        Rz = [cos(rtx_heading_w2(k)) -sin(rtx_heading_w2(k)) 0; sin(rtx_heading_w2(k)) cos(rtx_heading_w2(k)) 0; 0 0 1];    % yaw
        Ry = [cos(rtx_pitch_w2(k)) 0 sin(rtx_pitch_w2(k)); 0 1 0; -sin(rtx_pitch_w2(k)) 0 cos(rtx_pitch_w2(k))]; %pitch
        Rx = [1 0 0; 0 cos(rtx_roll_w2(k)) -sin(rtx_roll_w2(k)); 0 sin(rtx_roll_w2(k)) cos(rtx_roll_w2(k))]; % roll
        R = Rz*Ry*Rx;     
        
       % with respect to the granite block
        rtx_s3_handle2 = R^-1*Rs_gb*(rtx_ant_gb(:,k) - posmv_gn_gb(:,k));
        gga_s3_handle2 = R^-1*Rs_gb*(gga_ant_gb(:,k) - posmv_gn_gb(:,k));
        bin_s3_handle2 =  R^-1*Rs_gb*(posmvbin_ant_gb(:,k) - posmv_gn_gb(:,k)); % Sensor 2
        bin2_s3_handle2 =  R^-1*Rs_gb*(posmvbin_ant_gb2(:,k) - posmv_gn_gb(:,k)); % Sensor1
        
        td_s3_handle2 =  R^-1*Rs_gb*(posmv_td_gb(:,k) - posmv_gn_gb(:,k)); % Sensor2 with respect to the granite block
        
        gn_s3_handle3 =  ant+R^-1*Rs_gb*(posmv_gn_gb(:,k) - rtx_ant_gb(:,k)); % Sensor1 with respect to the granite block
        td_s3_handle3 =  ant+R^-1*Rs_gb*(posmv_td_gb(:,k) - rtx_ant_gb(:,k)); % Sensor2 with respect to the granite block
        gn_s3_handle4 =  ant+R^-1*Rs_gb*(posmv_gn_gb(:,k) - gga_ant_gb(:,k)); % Sensor1 with respect to the granite block
        td_s3_handle4 =  ant+R^-1*Rs_gb*(posmv_td_gb(:,k) - gga_ant_gb(:,k)); % Sensor2 with respect to the granite block
        
        
        rtx_s3_refgn = horzcat(rtx_s3_refgn,rtx_s3_handle2);
        gga_s3_refgn = horzcat(gga_s3_refgn,gga_s3_handle2);
        bin_s3_refgn = horzcat(bin_s3_refgn,bin_s3_handle2);
        bin2_s3_refgn = horzcat(bin2_s3_refgn,bin2_s3_handle2);
        
        td_s3_refgn = horzcat(td_s3_refgn,td_s3_handle2);
        
        gn_s3_refrtx = horzcat(gn_s3_refrtx,gn_s3_handle3);
        td_s3_refrtx = horzcat(td_s3_refrtx,td_s3_handle3);
        gn_s3_refgga= horzcat(gn_s3_refgga,gn_s3_handle4);
        td_s3_refgga= horzcat(td_s3_refgga,td_s3_handle4);
        
    end

   % store variables
   % RTX
    rtx_t_all =vertcat(rtx_t_all,rtx_t_w);
    rtx_x_ant_gb = vertcat(rtx_x_ant_gb,rtx_x_ant);
    rtx_y_ant_gb = vertcat(rtx_y_ant_gb,rtx_y_ant);
    rtx_z_ant_gb = vertcat(rtx_z_ant_gb,rtx_altitude_w);
    
   % GGA
    posmv_t_all = vertcat(posmv_t_all,gga_t_w);
    posmv_x_ant_gb = vertcat(posmv_x_ant_gb,gga_x_ant);
    posmv_y_ant_gb = vertcat(posmv_y_ant_gb,gga_y_ant);
    posmv_z_ant_gb = vertcat(posmv_z_ant_gb,gga_altitude_w);
    
    
   % POSMV concatinate
    posmvbin_x_ant_gb = vertcat(posmvbin_x_ant_gb,posmvbin_x_d);
    posmvbin_y_ant_gb = vertcat(posmvbin_y_ant_gb,posmvbin_y_d);
    posmvbin_z_ant_gb =vertcat(posmvbin_z_ant_gb,posmvbin_z_d);
    
    posmvbin_x_ant_gb2 = vertcat(posmvbin_x_ant_gb2,posmvbin_x_d2);
    posmvbin_y_ant_gb2 = vertcat(posmvbin_y_ant_gb2,posmvbin_y_d2);
    posmvbin_z_ant_gb2 =vertcat(posmvbin_z_ant_gb2,posmvbin_z_d2);
    
    posmv_td_all = horzcat(posmv_td_all,posmv_td_gb);
    posmv_gn_all = horzcat(posmv_gn_all,posmv_gn_gb);
    
    posmv_heading = vertcat(posmv_heading,posmv_headingd(ind));   
    posmv_xvel = vertcat(posmv_xvel,posmv_xveld(ind));
    posmv_yvel = vertcat(posmv_yvel,posmv_yveld(ind));
    posmv_zvel = vertcat(posmv_zvel,posmv_zveld(ind));
   
    
    offset_rtx_ref_gn_s3 = horzcat(offset_rtx_ref_gn_s3,rtx_s3_refgn);
    offset_posmvgga_ref_gn_s3 = horzcat(offset_posmvgga_ref_gn_s3,gga_s3_refgn);
    offset_posmvbins2_ref_gn_s3 = horzcat(offset_posmvbins2_ref_gn_s3,bin_s3_refgn);
    offset_posmvbins1_ref_gn_s3 = horzcat(offset_posmvbins1_ref_gn_s3,bin2_s3_refgn);
    
    posmvbins2_ref_gn_s3 = horzcat(posmvbins2_ref_gn_s3,td_s3_refgn);

    posmvbins2_ref_rtx_s3 = horzcat(posmvbins2_ref_rtx_s3,td_s3_refrtx);
    posmvbins1_ref_rtx_s3 = horzcat(posmvbins1_ref_rtx_s3,gn_s3_refrtx);
    posmvbins2_ref_gga_s3 = horzcat(posmvbins2_ref_gga_s3,td_s3_refgga);
    posmvbins1_ref_gga_s3 = horzcat(posmvbins1_ref_gga_s3,gn_s3_refgga);

    

%% save file

 cd /Users/testuser/Documents/MATLAB/Script/Data/
 
 clearvars -except posmvbins1_ref_gga_s3 posmvbins2_ref_gga_s3 posmv_gn_all posmv_td_all posmvbins2_ref_gn_s3 posmvbins2_ref_rtx_s3 posmvbins1_ref_rtx_s3 offset_posmvgga_ref_gn_s3 offset_posmvbins1_ref_gn_s3 offset_posmvbins2_ref_gn_s3 posmvbin_x_ant_gb posmvbin_y_ant_gb posmvbin_z_ant_gb rtx_t_all rtx_x_ant_gb rtx_y_ant_gb rtx_z_ant_gb posmv_t_all posmv_x_ant_gb posmv_y_ant_gb posmv_z_ant_gb posmv_heading posmv_xvel posmv_yvel posmv_zvel posmv_xacc posmv_yacc posmv_zacc offset_rtx_ref_gn_s3 ant rtx_X_dis rtx_Y_dis rtx_Z_dis rtx_toffset gga_X_dis gga_Y_dis gga_Z_dis gga_toffset

save RTX_GGA_BIN_compare_wo_interpolation_beforecruise2


