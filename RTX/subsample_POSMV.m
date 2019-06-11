% POS-MV
clear 
close all
% sub-sample pos mv
    cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/posmv/all
    d_posmv=dir('pos*.mat');  
   
    % load posmv files
    for i=1:length(d_posmv)
        cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/posmv/all
        fname=d_posmv(i).name;
        load(fname)

        t_pos = posmv.t;
        lat = posmv.lat;
        lon = posmv.lon;
        altitude = posmv.altitude;
        heading = posmv.heading;
        roll = posmv.roll;
        pitch = posmv.pitch;
        heave = posmv.heave;
        altitude = altitude-2.31;           %WSG84 Geoid height relative to MSL (at ACO)
        
        % subsample by d=10    
        t_pos_1Hz=t_pos(1:10:end);
        lat_1Hz=lat(1:10:end);
        lon_1Hz=lon(1:10:end);
        altitude_1Hz = altitude(1:10:end);
        heading_1Hz = heading(1:10:end);
        roll_1Hz = roll(1:10:end);
        pitch_1Hz = pitch(1:10:end);
        heave_1Hz = heave(1:10:end);   

        posmv.t = t_pos_1Hz;
        posmv.lat = lat_1Hz;
        posmv.lon = lon_1Hz;
        posmv.altitude = altitude_1Hz;
        posmv.heading = heading_1Hz;
        posmv.roll = roll_1Hz;
        posmv.pitch = pitch_1Hz;
    cd ../1Hz_4_RTX
        savename = [fname(1:end-4) '_1Hz'];
        save(savename,'posmv');
    end

   

