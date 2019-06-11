% Compare RTX and POSMV tx transducer positions in tx files
% 1. Loop over each hourly file
% 2. Load data
% 3. transform lat/lon to x,y
% 4. Do comparisons
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
close all
% POS MV and RTX tx file names
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/tx_file/all
d_posmv = dir('tx*');
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/RTX/tx_file/corrected
d_rtx = dir('rtx*');
% 1. Loop
utmzone = '4Q';
posmv_t = [];
rtx_t = [];
posmv_toffset = [];
rtx_toffset = [];
posmv_lat = [];
rtx_lat = [];
posmv_lon = [];
rtx_lon = [];
posmv_altitude = [];
rtx_altitude = [];
posmv_xvel = [];
rtx_xvel = [];
posmv_zvel = [];
rtx_zvel = [];
posmv_heading = [];
rtx_heading = [];
posmv_roll = [];
posmv_pitch = [];
td = [];

for k = 35:37%length(d_posmv)
   % 2. Load data
   cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/tx_file/all
   posmv_fname = d_posmv(k).name
   load(posmv_fname)
   posmv_ta = tx_data.t;
   posmv_lata = tx_data.lat;
   posmv_lona = tx_data.lon;
   posmv_altitudea = tx_data.altitude;
   %posmv_toffseta = tx_data.t_offset;
   posmv_xvela = tx_data.x_vel;
   posmv_zvela = tx_data.z_vel;
   posmv_headinga = tx_data.heading;
   posmv_rolla = tx_data.roll;
   posmv_pitcha = tx_data.pitch;
   
   cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/RTX/tx_file/corrected
   rtx_fname = d_rtx(k).name
   load(rtx_fname)
   rtx_ta = rtx_tx_data.t;
   rtx_lata = rtx_tx_data.lat;
   rtx_lona = rtx_tx_data.lon;
   rtx_altitudea = rtx_tx_data.altitude;
   rtx_xvela = rtx_tx_data.xvel;
   rtx_zvela = rtx_tx_data.zvel;
   rtx_toffseta = rtx_tx_data.t_offset;
   rtx_headinga = rtx_tx_data.heading;
   
   
   % match POS MV and RTX tx times
   rm_ind = [];
   for i = 1:length(posmv_ta)
      [t_diff,I] = min(abs(posmv_ta(i) - rtx_ta)); 
      td(end+1) = t_diff*3600*24;
      if abs(t_diff*3600*24) > 1 
          rm_ind(end+1) = i;
      end
       
   end
  
   posmv_ta(rm_ind) = [];
   posmv_lata(rm_ind) = [];
   posmv_lona(rm_ind) = [];
   posmv_altitudea(rm_ind) = [];
   %posmv_toffseta(rm_ind) = [];
   posmv_xvela(rm_ind) = [];
   posmv_zvela(rm_ind) = [];
   posmv_headinga(rm_ind) = [];
   posmv_rolla(rm_ind) = [];
   posmv_pitcha(rm_ind) = [];
   
   
   %%% concatenate the data
    posmv_t = horzcat(posmv_t,posmv_ta);
    rtx_t = horzcat(rtx_t,rtx_ta);
%     posmv_toffset = horzcat(posmv_toffset,posmv_toffseta);
    rtx_toffset = horzcat(rtx_toffset,rtx_toffseta);
    posmv_lat = horzcat(posmv_lat,posmv_lata);
    rtx_lat = horzcat(rtx_lat,rtx_lata);
    posmv_lon = horzcat(posmv_lon,posmv_lona);
    rtx_lon = horzcat(rtx_lon,rtx_lona);
    posmv_altitude = horzcat(posmv_altitude,posmv_altitudea);
    rtx_altitude = horzcat(rtx_altitude,rtx_altitudea);
    posmv_xvel = horzcat(posmv_xvel,posmv_xvela);
    rtx_xvel = horzcat(rtx_xvel,rtx_xvela);
    posmv_zvel = horzcat(posmv_zvel,posmv_zvela);
    rtx_zvel = horzcat(rtx_zvel,rtx_zvela); 
    posmv_heading = horzcat(posmv_heading,posmv_headinga);
    rtx_heading = horzcat(rtx_heading,rtx_headinga);
    posmv_roll = horzcat(posmv_roll,posmv_rolla);
    posmv_pitch = horzcat(posmv_pitch,posmv_pitcha);
     
end
%%
% RTX transducer offset from POS MV
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

% to UTM coordinate system
[posmv_x,posmv_y,utmzone] = deg2utm(posmv_lat,posmv_lon);
[rtx_x,rtx_y,utmzone] = deg2utm(rtx_lat,rtx_lon);
posmv_tx5_gb = [posmv_x';posmv_y';posmv_altitude];
rtx_tx5_gb = [rtx_x';rtx_y';rtx_altitude];

% transform lat/lon in the global UTM coordinate back to ship locally level frame
tx5_rtx_s3_reftx = [];
posmv_roll = posmv_roll/180*pi;
posmv_pitch = posmv_pitch/180*pi;
posmv_heading = posmv_heading/180*pi;
     
    for k = 1:length(posmv_roll)
        Rz = [cos(posmv_heading(k)) -sin(posmv_heading(k)) 0; sin(posmv_heading(k)) cos(posmv_heading(k)) 0; 0 0 1];    % yaw
        Ry = [cos(posmv_pitch(k)) 0 -sin(posmv_pitch(k)); 0 1 0; sin(posmv_pitch(k)) 0 cos(posmv_pitch(k))]; %pitch
        Rx = [1 0 0; 0 cos(posmv_roll(k)) sin(posmv_roll(k)); 0 -sin(posmv_roll(k)) cos(posmv_roll(k))]; % roll
        R = Rz*Ry*Rx;     
        % the offset of RTX's transducer from POSMV's
        tx5_rtx_s3_handle = R^-1*Rs_gb*(rtx_tx5_gb(:,k) - posmv_tx5_gb(:,k));
        tx5_rtx_s3_reftx = horzcat(tx5_rtx_s3_reftx,tx5_rtx_s3_handle);
       
    end



% % Time offset
% figure(1)
% subplot(2,1,1)
% scatter(posmv_t,-posmv_toffset*3600*24*1000,20,'filled')
% grid on
% ylabel('msec')
% xlabel('Time')
% title('POSMV time offset from the actual tx time')
% datetick('x','dd hh:MM')
% 
% subplot(2,1,2)
% scatter(rtx_t,rtx_toffset*1000,20,'filled')
% grid on
% ylabel('msec')
% xlabel('Time')
% title('RTX time offset from the actual tx time')
% datetick('x','dd hh:MM')

%% Transducer Offset
figure(2)
subplot(2,1,1)
scatter(rtx_t,tx5_rtx_s3_reftx(1,:),20,posmv_xvel,'filled')
datetick('x')
grid on
c = colorbar;
c.Label.String = 'X speed (m/s)';
colormap jet
title('RTX transducer offset from the POS MV ')
ylabel('Bow/Stern (m)')

subplot(2,1,2)
scatter(rtx_t,tx5_rtx_s3_reftx(3,:),20,posmv_zvel,'filled')
datetick('x')
grid on
c = colorbar;
c.Label.String = 'Z speed (m/s)';
colormap jet
ylim([-4 4])
ylabel('Vertical (m)')





