% Compare corrected and uncorrected RTX tx files
clear
close all
%% uncorrected
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/RTX/tx_file/uncorrected
day = 28;            %  Edit
hour = 13:15;         % Edit
% create a set of file names
fname = [];
for i = 1:length(hour)
    if hour(i) < 10
        fname_d = "rtx_tx_data_2018_10_"+string(day)+"_0"+string(hour(i));
    else
        fname_d = "rtx_tx_data_2018_10_"+string(day)+"_"+string(hour(i));
    end
   fname = vertcat(fname,fname_d);
end


%Load Tx data
rtx_tx_t_uc = [];
rtx_tx_lat_uc = [];
rtx_tx_lon_uc = [];
rtx_tx_altitude_uc = [];

for p = 1:length(fname)
    lname = fname(p);
    load(lname)
    rtx_tx_t_uc = horzcat(rtx_tx_t_uc,rtx_tx_data.t);      %TX time                                  
    rtx_tx_lat_uc = horzcat(rtx_tx_lat_uc,rtx_tx_data.lat);   %TX lat
    rtx_tx_lon_uc = horzcat(rtx_tx_lon_uc,rtx_tx_data.lon);
    rtx_tx_altitude_uc = horzcat(rtx_tx_altitude_uc,rtx_tx_data.altitude);
    
end

% to UTM
[x_uc,y_uc,utmzone] = deg2utm(rtx_tx_lat_uc,rtx_tx_lon_uc);



%% corrected
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/RTX/tx_file/corrected
% create a set of file names
fname = [];
for i = 1:length(hour)
    if hour(i) < 10
        fname_d = "rtx_tx_data_2018_10_"+string(day)+"_0"+string(hour(i));
    else
        fname_d = "rtx_tx_data_2018_10_"+string(day)+"_"+string(hour(i));
    end
   fname = vertcat(fname,fname_d);
end

%Load Tx data
rtx_tx_t_c = [];
rtx_tx_lat_c = [];
rtx_tx_lon_c = [];
rtx_tx_altitude_c = [];
t_offset = [];
xvel = [];
zvel = [];
heading = [];
t_offset = [];
for p = 1:length(fname)
    lname = fname(p);
    load(lname)
    rtx_tx_t_c = horzcat(rtx_tx_t_c,rtx_tx_data.t);      %TX time                                  
    rtx_tx_lat_c = horzcat(rtx_tx_lat_c,rtx_tx_data.lat);   %TX lat
    rtx_tx_lon_c = horzcat(rtx_tx_lon_c,rtx_tx_data.lon);
    rtx_tx_altitude_c = horzcat(rtx_tx_altitude_c,rtx_tx_data.altitude);
    t_offset = horzcat(t_offset,rtx_tx_data.t_offset);  
    xvel = horzcat(xvel,rtx_tx_data.xvel);
    zvel = horzcat(zvel,rtx_tx_data.zvel);
    heading = horzcat(heading,rtx_tx_data.heading);
   
end

% to UTM
[x_c,y_c,utmzone] = deg2utm(rtx_tx_lat_c,rtx_tx_lon_c);

%% Plot
diff_x = x_c - x_uc;
diff_y = y_c - y_uc;
diff_z = rtx_tx_altitude_c - rtx_tx_altitude_uc';

figure(1)
clf
subplot(3,1,1)
scatter(rtx_tx_t_uc,diff_x,[],heading,'filled')
c = colorbar;
c.Label.String = 'Heading';
colormap jet
datetick('x')
grid on
axis tight
title('X offset')
ylabel('m')
yticks(-2:0.5:2 )
ylim([-2 2])

subplot(3,1,2)
scatter(rtx_tx_t_uc,diff_y,[],heading,'filled')
datetick('x')
grid on
axis tight
c = colorbar;
c.Label.String = 'Heading';
colormap jet
title('Y offset')
ylabel('m')
yticks(-2:1:2 )
ylim([-3 3])

subplot(3,1,3)
scatter(rtx_tx_t_uc,xvel,[],heading,'filled')
datetick('x')
grid on
axis tight
c = colorbar;
c.Label.String = 'Heading';
colormap jet
title('Logitudinal Speed (m/s)')
ylabel('m')

%% time offset and distance offset
figure(2)
scatter(t_offset,diff_x,[],xvel,'filled')
grid on
xlabel('Time Offset (s)')
ylabel('X Offset (m)')
c = colorbar;
c.Label.String = 'Ship Velocity (Longitudinal (m/s)';
colormap jet
caxis([-1 3.5])
title('Transducer Position Correction: Westward Transect ')
ylim([-2 2])

figure(3)
scatter(t_offset,diff_y,[],xvel,'filled')
grid on
xlabel('Time Offset (s)')
ylabel('Y Offset (m)')
c = colorbar;
c.Label.String = 'Ship Velocity (Longitudinal (m/s)';
colormap jet
caxis([-1 3.5])
title('Transducer Position Correction: Westward Transect ')
ylim([-2 2])

figure(33)
scatter(diff_x,diff_y,20,xvel,'filled')
xlim([-2 2])
ylim([-1 1])
grid on
xlabel('X offset(m)')
ylabel('Y offset(m)')
c = colorbar;
c.Label.String  = 'Longitudinal Speed (m/s)';
colormap jet
%% time offset
figure(4)
scatter(rtx_tx_t_c,t_offset,'filled')
grid on
xlabel('Time')
ylabel('second')
datetick('x')
title('RTX time offset from the actual tx time')

%% Map plot
figure(5)
scatter(rtx_tx_lon_uc,rtx_tx_lat_uc,15,rtx_tx_t_uc,'filled')
grid on
colormap jet


