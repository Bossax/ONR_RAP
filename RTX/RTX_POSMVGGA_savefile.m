% compared RTX and POS-MV GGA
clear
close all

% load RTX
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/RTX/before
load('RTX_beforecruise_10Hz.mat')

rtx_t = transpose(RTX.rtx_t);
rtx_lat = RTX.rtx_lat;
rtx_lon = RTX.rtx_lon;
rtx_altitude = RTX.rtx_altitude;

% load POSMV GGA
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/posmv/GGA/before
load('posmv_gga_beforecruise')
posmv_t = posmv_gga.t;
posmv_lat=posmv_gga.lat;
posmv_lon = posmv_gga.lon;
posmv_altitude = posmv_gga.altitude;

clear RTX posmv_gga
%% truncate data
tstart = datenum('2018-10-24 00:00:00','yyyy-mm-dd HH:MM:SS');
tend = datenum('2018-10-24 19:30:00','yyyy-mm-dd HH:MM:SS');
pos_ind = find((posmv_t >=tstart)&(posmv_t <= tend));
rtx_ind = find((rtx_t >= tstart)&(rtx_t <= tend));

rtx_t = rtx_t(rtx_ind); 
rtx_lat = rtx_lat(rtx_ind);
rtx_lon = rtx_lon(rtx_ind);
rtx_altitude = rtx_altitude(rtx_ind);

posmv_t = posmv_t(pos_ind);
posmv_lat=posmv_lat(pos_ind);
posmv_lon = posmv_lon(pos_ind);
posmv_altitude = posmv_altitude(pos_ind);
%% match RTX with POSMV
keep_ind = [];
for iii =1:length(rtx_t)
    if mod(iii,1000) == 0
        iii/length(rtx_t)*100;
    end
    [tdiff,~] = min(abs(rtx_t(iii)- posmv_t));
    if tdiff*3600*24 <=0.5
        keep_ind(end+1) = iii;
    end
end
rtx_t = rtx_t(keep_ind);
rtx_lat= rtx_lat(keep_ind);
rtx_lon = rtx_lon(keep_ind);
rtx_altitude = rtx_altitude(keep_ind);

%% match POSMV with RTX
keep_ind = [];
for iii =1:length(posmv_t)
    [tdiff,~] = min(abs(posmv_t(iii)- rtx_t));
    if tdiff*3600*24 <=0.5
        keep_ind(end+1) = iii;
    end
end

posmvgga_t = posmv_t(keep_ind);
posmvgga_lat= posmv_lat(keep_ind);
posmvgga_lon = posmv_lon(keep_ind);
posmvgga_altitude = posmv_altitude(keep_ind);

%% save
rtx_t = rtx_t';
cd /Users/testuser/Documents/MATLAB/Script/Data
save rtx10Hz_posmvGGA_beforecruise posmvgga_t posmvgga_lat posmvgga_lon posmvgga_altitude rtx_t rtx_lat rtx_lon rtx_altitude 



