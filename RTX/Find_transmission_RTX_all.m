% Find transmission time and transducer position based on RTX 
% 1. RTX file truncation based on Time of Day and add transducer positions
% 2. extract transmission times from tx_file, match RTX data points 
% 3. find corresponding RTX transducer position
clear
close all
%% Download RTX file
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/RTX/
load kilo_moana_1819_RTX.dat
A = kilo_moana_1819_RTX;
rtx_lat = transpose(A(:,2));
rtx_lon = transpose(A(:,3));
rtx_altitude = transpose(A(:,4));
%% Format timestamps
% Modified Julian Date 
mjd = A(:,1);
jd = mjd +2400000.5;
rtx_t = [];
% Convert to Gregorian Date
for i = 1:length(jd)
   rtx_t(end+1) = datenum(datetime(jd(i),'ConvertFrom','juliandate'));
end
%% Download POS MV files
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/tx_file/all
d = dir('tx*');

%% Loop over all files and read the tx times
% find nearest corresponding rtx times
% calculate transducer position
% pack data and save
tx_lat_tot = [];
tx_lon_tot = [];


t_diff = [];        % offset between all RTX times and transmission times
t_1sec = 0;
rtx_t_f = [];
tx_t_f  = [];
for ii = 1:length(d)
    cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/tx_file/all
    % extract tx posmv data
     fname = d(ii).name;
    load (fname)
    posmv_tx_roll = tx_data.roll;
    posmv_tx_pitch = tx_data.pitch;
    posmv_tx_yaw = tx_data.heading;
    posmv_heave = tx_data.heave;
    posmv_tx_t = tx_data.t;
    ind_pos = [];
    rm_ind_t = [];
    tx_lat_rtx = [];
    tx_lon_rtx = [];
    tx_altitude_rtx = [];
    tx_heading_rtx = [];
    % 2 match transmission times with RTX data point
    for k = 1:length(posmv_tx_t)
        % find corresponding rtx data points to use its antenna position
        [t_difff,I] = min(abs(posmv_tx_t(k) - rtx_t));
        t_diff(end+1) = t_difff*3600*24;
        if t_difff*3600*24 > 0.5          % discard a sample with a gap larger than 0.5 sec
            rm_ind_t(end+1) = k;
            t_1sec = t_1sec+1;
        else
            ind_pos(end+1) = I;     % index to extract lon/lat/alt
        end
    end
    rtx_t_w = posmv_tx_t;   
    rtx_t_w(rm_ind_t) = [];
    
    tx_heading_rtx = posmv_tx_yaw;
    tx_heading_rtx(rm_ind_t) = [];
    % antenna position
    rtx_t_pick = rtx_t(ind_pos);
    
    % lat lon altitude fro RTX
    rtx_lat_w = rtx_lat(ind_pos);
    rtx_lon_w = rtx_lon(ind_pos);
    rtx_altitude_w = rtx_altitude(ind_pos);
    
     % time diff anal
     rtx_t_f = horzcat(rtx_t_f,rtx_t_pick);
     tx_t_f = horzcat(tx_t_f,rtx_t_w);   
    
    % calculate transducer position (global coordinate) of the actual rtx time
    [td_lat,td_lon,td_altitude] = transducer_pos_rtx(rtx_lat_w,rtx_lon_w,rtx_altitude_w,rtx_t_w,posmv_tx_roll,posmv_tx_pitch,posmv_tx_yaw,posmv_tx_t);
        
     tx_lat_rtx = transpose(vertcat(tx_lat_rtx,td_lat));
     tx_lon_rtx = transpose(vertcat(tx_lon_rtx,td_lon));
     tx_altitude_rtx = transpose(vertcat(tx_altitude_rtx,td_altitude));
        
       
    % 7. pack data and save
    tx_lat_tot = horzcat(tx_lat_tot,tx_lat_rtx) ;
    tx_lon_tot = horzcat(tx_lon_tot,tx_lon_rtx) ;
    
    size(rtx_t_w)
    size(tx_lat_rtx)
    rtx_tx_data.t = rtx_t_w;
    rtx_tx_data.lat = tx_lat_rtx;
    rtx_tx_data.lon = tx_lon_rtx;
    rtx_tx_data.altitude = tx_altitude_rtx - 2.31;     % subtract 2.31 Geiod Height above MSL
    rtx_tx_data.heading = tx_heading_rtx;
    
    cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/RTX/tx_file/
    savename = sprintf("rtx_"+string(fname(1:end-4))+".mat")
    save(savename,'rtx_tx_data');
    
 
end
