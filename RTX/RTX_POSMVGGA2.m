% compare RTX and POSMV GGA ans POSMV binary antenna positions
% transform to the ship position
% plot the position offsets in 3 directions and as a function of velocity
clear
%close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load RTX and POSMV GGA data
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Script/Data
load('rtx_posmvGGA.mat')
rtx_t = rtx_t - 2/(3600*24);        % adjust rtx time
% Load POSMV to extract heading and velo city
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/posmv/all
d_posmv = dir('posmv*');
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/posmv/granite_block/all
d_posmv_gn = dir('gb*');
%%
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
posmvbin_x_ant_gb = [];
posmvbin_y_ant_gb = [];
posmvbin_z_ant_gb =[];

% parameters
posmv_heading = [];
posmv_xvel = [];
posmv_yvel = [];
posmv_zvel = [];
posmv_xacc = [];
posmv_yacc = [];
posmv_zacc = [];

% RTX/GGA/Bin offset from POSMV's granite block (ship locally level frame)
offset_rtx_ref_gn_s3 = [];
offset_posmvgga_ref_gn_s3 = [];
offset_posmvbin_ref_gn_s3 = [];

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
% Loop over posmv binary files
% EDIT start file and end file
for ii = 5:length(d_posmv)-1
    ii
    cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/posmv/all
    % Load POSMV Binary file
    % transducer
    fname = d_posmv(ii).name;
    load (fname)
    posmv_latd = posmv.lat;
    posmv_lond = posmv.lon;
    posmv_altituded =posmv.altitude;
    posmv_td = posmv.t;
    
    % vessel dynamics
    posmv_xveld = posmv.x_vel;
    posmv_yveld = posmv.y_vel;
    posmv_zveld = posmv.z_vel;
    posmv_xacceld = posmv.accel_lon;
    posmv_yacceld = posmv.accel_tran;
    posmv_zacceld = posmv.accel_down;
    posmv_rolld = posmv.roll;
    posmv_pitchd = posmv.pitch;
    posmv_headingd = posmv.heading;
    
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
        posmv_xacceld(ind_l) = [];
        posmv_yacceld(ind_l) = [];
        posmv_zacceld(ind_l) = [];
  
    end
    % load granite block
    cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/posmv/granite_block/all
    fname2 = d_posmv_gn(ii).name;
    load(fname2)
    posmv_lat_gnd = posmv_gb.lat;
    posmv_lon_gnd = posmv_gb.lon;
    posmv_altitude_gnd =posmv_gb.altitude;
    posmv_t_gnd = posmv_gb.t;
    
    
    
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
    
   
    %%%%%Interpolation RTX/GGA to POS-MV binnary data points%%%%%%%%%%%%%%%%%%%
    %%%%% Output: rtx_lat/lon/altitude gga_lat/lon/altitude  %%%%%%
    
    
    
    
    
    % calculate POSMV Bin antenna position (global coordinate)
    [pos_ant_lat,pos_ant_lon,pos_ant_altitude] = ant_pos_posmv(posmv_latd(ind),posmv_lond(ind),posmv_altituded(ind),posmv_td(ind),posmv_rolld(ind),posmv_pitchd(ind),posmv_headingd(ind));

    % convert lat lon to UTM x,y,z  granite block POSMV 
    [posmv_x_gnd,posmv_y_gnd,utmzone] = deg2utm(posmv_lat_gnd(ind),posmv_lon_gnd(ind));
    posmv_z_gnd = posmv_altitude_gnd(ind);
    posmv_gn_gb = [posmv_x_gnd'; posmv_y_gnd' ; posmv_z_gnd'];  
    
    % RTX/GGA position
    [rtx_x_ant,rtx_y_ant,utmzone] = deg2utm(rtx_lat_w,rtx_lon_w);
     rtx_ant_gb = [rtx_x_ant';rtx_y_ant'; rtx_altitude_w'];
     
     [gga_x_ant,gga_y_ant,utmzone] = deg2utm(gga_lat_w,gga_lon_w);
     gga_ant_gb = [gga_x_ant';gga_y_ant'; gga_altitude_w'];
     
    [posmvbin_x_d,posmvbin_y_d,utmzone] = deg2utm(pos_ant_lat,pos_ant_lon);
    posmvbin_z_d = pos_ant_altitude;
    posmvbin_ant_gb = [posmvbin_x_d'; posmvbin_y_d' ; posmvbin_z_d'];
    
     
    % transform positionin the global UTM coordinate back to ship locally level frame
    rtx_s3_refgn = [];
    gga_s3_refgn = [];
    bin_s3_refgn = [];
   
    rtx_roll_w2 = rtx_roll_w/180*pi; 
    rtx_pitch_w2 = rtx_pitch_w/180*pi; 
    rtx_heading_w2 = rtx_heading_w/180*pi; 
    
    % Perform coordinate transformation from UTM to ship locally level frame
    for k = 1:length(rtx_roll_w)
        Rz = [cos(rtx_heading_w2(k)) -sin(rtx_heading_w2(k)) 0; sin(rtx_heading_w2(k)) cos(rtx_heading_w2(k)) 0; 0 0 1];    % yaw
        Ry = [cos(rtx_pitch_w2(k)) 0 -sin(rtx_pitch_w2(k)); 0 1 0; sin(rtx_pitch_w2(k)) 0 cos(rtx_pitch_w2(k))]; %pitch
        Rx = [1 0 0; 0 cos(rtx_roll_w2(k)) sin(rtx_roll_w2(k)); 0 -sin(rtx_roll_w2(k)) cos(rtx_roll_w2(k))]; % roll
        R = Rz*Ry*Rx;     
        % the offset of RTX's transducer from POSMV's
        
        rtx_s3_handle2 = R^-1*Rs_gb*(rtx_ant_gb(:,k) - posmv_gn_gb(:,k));
        gga_s3_handle2 = R^-1*Rs_gb*(gga_ant_gb(:,k) - posmv_gn_gb(:,k));
        bin_s3_handle2 =  R^-1*Rs_gb*(posmvbin_ant_gb(:,k) - posmv_gn_gb(:,k));
        
        rtx_s3_refgn = horzcat(rtx_s3_refgn,rtx_s3_handle2);
        gga_s3_refgn = horzcat(gga_s3_refgn,gga_s3_handle2);
        bin_s3_refgn = horzcat(bin_s3_refgn,bin_s3_handle2);
        
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
    posmv_heading = vertcat(posmv_heading,posmv_headingd(ind));   
    posmv_xvel = vertcat(posmv_xvel,posmv_xveld(ind));
    posmv_yvel = vertcat(posmv_yvel,posmv_yveld(ind));
    posmv_zvel = vertcat(posmv_zvel,posmv_zveld(ind));
    posmv_xacc = vertcat(posmv_xacc,posmv_xacceld(ind));
    posmv_yacc = vertcat(posmv_yacc,posmv_yacceld(ind));
    posmv_zacc = vertcat(posmv_zacc,posmv_zacceld(ind));
    
    
    offset_rtx_ref_gn_s3 = horzcat(offset_rtx_ref_gn_s3,rtx_s3_refgn);
    offset_posmvgga_ref_gn_s3 = horzcat(offset_posmvgga_ref_gn_s3,gga_s3_refgn);
   offset_posmvbin_ref_gn_s3 = horzcat(offset_posmvbin_ref_gn_s3,bin_s3_refgn);
    
end
%%
clearvars -except offset_posmvgga_ref_gn_s3 offset_posmvbin_ref_gn_s3 posmvbin_x_ant_gb posmvbin_y_ant_gb posmvbin_z_ant_gb rtx_t_all rtx_x_ant_gb rtx_y_ant_gb rtx_z_ant_gb posmv_t_all posmv_x_ant_gb posmv_y_ant_gb posmv_z_ant_gb posmv_heading posmv_xvel posmv_yvel posmv_zvel posmv_xacc posmv_yacc posmv_zacc offset_rtx_ref_gn_s3  ant 
%%
offset_rtx_posmv_sn3 = offset_rtx_ref_gn_s3-offset_posmvbin_ref_gn_s3;

%% plot

% time frame
tstart = datenum('20181028 12:00:01','yyyymmdd HH:MM:SS');
tend= datenum('20181028 18:00:00','yyyymmdd HH:MM:SS');
[~,indstart] = min(abs(rtx_t_all - tstart));
[~,indend] = min(abs(rtx_t_all - tend));



med_xoff = median(offset_rtx_posmv_sn3(1,indstart:indend))
med_yoff = median(offset_rtx_posmv_sn3(2,indstart:indend))
med_zoff = median(offset_rtx_posmv_sn3(3,indstart:indend))
med_xvel = median(posmv_xvel(indstart:indend))
med_yvel = median(posmv_yvel(indstart:indend))
med_zvel = median(posmv_zvel(indstart:indend))

% RTX and GGA offset (RTX - POS MV)
figure(1)
clf
subplot(3,1,1)
scatter(rtx_t_all(indstart:indend),offset_rtx_posmv_sn3(1,indstart:indend),abs(posmv_xvel(indstart:indend)*7)+2,posmv_xacc(indstart:indend),'filled')
grid on
datetick('x')
ylim([ -6 6])
yticks([-16:2:6])
ylabel('Bow-Stern offset (m)')
headline = sprintf('RTX and POSMV Binary Difference (Antenna 2 seconds subtracted from the RTX)\n Median = %.3f m Median x Vel = %.3f m/s',med_xoff, med_xvel);
title(headline)
cbar = colorbar;
cbar.Label.String = 'X Acc (m/s^{2})';
colormap jet

subplot(3,1,2)
scatter(rtx_t_all(indstart:indend),offset_rtx_posmv_sn3(2,indstart:indend),abs(posmv_yvel(indstart:indend)*6)+5,posmv_yacc(indstart:indend),'filled')
grid on
datetick('x')
ylim([-6 6])
yticks([-9:2:8])
ylabel('Startboard-Portside offset (m)')
colorbar
colormap jet
cbar = colorbar;
cbar.Label.String = 'Y Acc (m/s^{2})';
headline = sprintf('Median = %.3f m Median y Vel = %.3f m/s',med_yoff,med_yvel);
title(headline)
caxis([-0.8 0.8])

subplot(3,1,3)
scatter(rtx_t_all(indstart:indend),offset_rtx_posmv_sn3(3,indstart:indend),abs(posmv_zvel(indstart:indend)*10)+5,posmv_zacc(indstart:indend),'filled')
grid on
datetick('x')
ylim([ -6 8])
yticks([-6:2:6])
ylabel('Vertical offset (m)')
cbar = colorbar;
cbar.Label.String = 'Z Acc (m/s^{2})';
colormap jet
headline = sprintf('Median = %.3f m Median z Vel = %.3f m/s',med_zoff,med_zvel);
title(headline)
caxis([-0.8 0.8])

