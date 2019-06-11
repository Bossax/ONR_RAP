% Compare statistics of RTX and POSMV tx transducer positions
% 1. Loop over each hourly file
% 2. Load data
% 3. transform lat/lon to x,y
% 4. Do stat cal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
close all
% POS MV and RTX file names
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/tx_file/all
d_posmv = dir('tx*');
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/RTX/tx_file/2
d_rtx = dir('rtx*');
% 1. Loop
utmzone = '4Q'
diff_x = [];
diff_y = [];
diff_z = [];
diff_t = [];
posmv_ta = [];
rtx_ta = [];
for k = 1:length(d_posmv)
   % 2. Load data
   cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/tx_file/all
   posmv_fname = d_posmv(k).name
   load(posmv_fname)
   posmv_t = tx_data.t;
   posmv_lat = tx_data.lat;
   posmv_lon = tx_data.lon;
   posmv_altitude = tx_data.altitude';
   cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/RTX/tx_file
   rtx_fname = d_rtx(k).name
   load(rtx_fname)
   rtx_t = rtx_tx_data.t;
   rtx_lat = rtx_tx_data.lat;
   rtx_lon = rtx_tx_data.lon;
   rtx_altitude = rtx_tx_data.altitude';
   
   posmv_ta = horzcat(posmv_ta,posmv_t);
   rtx_ta = horzcat(rtx_ta,posmv_t);
   % 3. Coordinate transformation
   [posmv_x,posmv_y,utmzone] = deg2utm(posmv_lat,posmv_lon);
   [rtx_x,rtx_y,utmzone] = deg2utm(rtx_lat,rtx_lon);
   % 4. Calculation
   diff_x_d = abs(posmv_x - rtx_x);
   diff_y_d = abs(posmv_y - rtx_y);
   diff_z_d = (posmv_altitude - rtx_altitude');
   diff_t_d = (posmv_t' - rtx_t')*3600*24;
   
   diff_x = vertcat(diff_x,diff_x_d);
   diff_y = vertcat(diff_y,diff_y_d);
   diff_z = vertcat(diff_z,diff_z_d);
   diff_t = vertcat(diff_t,diff_t_d);
   
end

%% Plot

figure(1)
subplot(2,2,1)
plot(posmv_ta,diff_t)
grid on
datetick('x','mmm/dd')
axis tight
grid minor
xlabel('Time')
ylabel('Time Difference (s)')
headline = sprintf('Median = %.3f s SD = %.3f s',median(diff_t),std(diff_t));
title(headline)

subplot(2,2,2)
plot(posmv_ta,diff_z)
grid on
datetick('x','mmm/dd')
axis tight
grid minor
xlabel('Time')
ylabel('Z Difference (m)')
headline = sprintf('Median = %.3f m SD = %.3f m',median(diff_z),std(diff_z));
title(headline)

subplot(2,2,3)
plot(posmv_ta,diff_x)
grid on
datetick('x','mmm/dd')
grid minor
xlabel('Time')
ylabel('X Difference (m)')
axis tight
headline = sprintf('Median = %.3f m SD = %.3f m',median(diff_x),std(diff_x));
title(headline)

subplot(2,2,4)
plot(posmv_ta,diff_y)
grid on
datetick('x','mmm/dd')
axis tight
grid minor
xlabel('Time')
ylabel('Y Difference (m)')
headline = sprintf('Median = %.3f m SD = %.3f m',median(diff_y),std(diff_y));
title(headline)

