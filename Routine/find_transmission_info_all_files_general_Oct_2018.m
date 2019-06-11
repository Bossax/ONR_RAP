%% Find Transmission time
% 1. Group Scarlette files created in the same hour (using file name)
% 2. Match POS MV file of the corresponding hour of each Scarlette file group
% no 1-second subtraction when calculating tx times
% compensate for 0.304 ms electronic delay
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read Scarlette file names EDIT
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Scarlett_Data
d=dir('LFM_*.mat');     

mod_date = [];
mod_hr = [];
% group files of the same hour based on start time
for i = 1:length(d)
    s = d(i).name;
    mod_date(end+1) = str2num(s(7:8))-69;     % EDIT to map Scarlette file name to actual date (27-30 Oct) ex. file name LFM_1396030114 96-69 = 27 start date of Oct cruise 
    mod_hr(end+1) = str2num(s(9:10));         % EDIT 
end

% extract pos_mv files date hr EDIT
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/posmv/all
d_pos=dir('pos*.mat');          

pos_date = [];
pos_hr = [];
for i = 1:length(d_pos)
    pos_name = d_pos(i).name;
    pos_date(end+1) = str2num(pos_name(end-7:end-6));   % edit
    pos_hr(end+1) = str2num(pos_name(end-5:end-4));     % edit
end

% Group Scarlett files created in the same hour
group = [];         % array whose rows contain groups and coulumns contain index of start file and end file of each group
counter = 1;        % counter of the Scarlett file
trig = 1;           % trigger creation of a new group
start_point = 1;
end_point = 1;

while true
    % first file of a group
    if trig == 1
        
        present_date = mod_date(counter);
        present_hr = mod_hr(counter);
        trig = 0;
        start_point = counter;
        counter = counter+1;
        
    end
    % check last file
    if counter > length(d)
        group(end+1,:) = [start_point counter-1];
        break
    end
    
 
    % consecutive files of a group
    % if still created in the same hour
   if trig == 0 && present_date == mod_date(counter) && present_hr == mod_hr(counter) 
       counter = counter+1
      
   else
        % end the group
        end_point = counter-1;
        group(end+1,:) = [start_point end_point];     
        trig =1;   % trigger a new group
        
   end
   
     % last Scarlett file is reached
        if counter == length(d)+1
            group(end+1,:) = [start_point counter];
        end

end

fs=44100;


%% Loop over all groups
for k = 1:length(group)     % 142 10/29 00:00
    
    cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Scarlett_Data  % EDIT
    tx_time=[];
    start_c = group(k,1);
    end_c = group(k,2);
    
    % loop files in the same group
    % find tx time of all files in the group
    for i=start_c:end_c
        y_tx=[];
        t_tx=[];
        v_time=[];
        pps_time=[];
        num_offset=[];
        
        fname=d(i).name;
        load(fname)


        pps=txrx.sam(2,:);              %PPS
        V=txrx.sam(4,:);                %Voltage
        t=txrx.taxis;
        start_t=txrx.start_time;

        t_tx=horzcat(t_tx,start_t+(t./(3600*24)));
        y_tx=horzcat(y_tx,V);


        peak_num=max(V);
        med = median(abs(V));
        if peak_num <= 4*med
            threshold = 4*med;
        else
            threshold = peak_num/2;
        end
        
        A=find(abs(V)>=threshold);          %Find voltage spike

        %Find start of voltage spikes
        for i=1:length(A)
            if i==1
                v_time(end+1)=A(i);
            else
                if (A(i)-A(i-1))>fs
                    v_time(end+1)=A(i);
                end
            end
        end
        v_time=v_time-1;


        A=find(abs(pps)>=.008);              %Find PPS spike

        %Find start of PPS spikes
        for i=1:length(A)
            if i==1
                pps_time(end+1)=A(i);
            else
                if (A(i)-A(i-1))>fs/2
                    pps_time(end+1)=A(i);
                end
            end
        end


        %Find offset between PPS spike and integer second
        for i=1:length(v_time)

            [~,B]=min(abs(v_time(i)-pps_time));
            num_offset(i)=t(pps_time(B))-round(t(pps_time(B)));

        end



        %Transmission time
        tx_time_sec=t(v_time)-num_offset;       % correct for time drift
        tx_time_sec = tx_time_sec+0.000304;     % cable delay 0.304 ms
        tx_time_date=start_t+(tx_time_sec./(3600*24));

        tx_time=horzcat(tx_time,tx_time_date);
    end

    %% POS-MV
    % search for corresponding posmv files
    present_date = mod_date(start_c);   % date of Scarlett file
    present_hr = mod_hr(start_c);       % hr of Scarlett file
    f_counter = 1;
    pos_file_index = [];
    % loop to find the corresponding posmv files
    % include the file of the next hour just in case
    for i = 1:length(d_pos)
       if present_date == pos_date(f_counter) && present_hr == pos_hr(f_counter) 
           pos_file_index(end+1,:) = [f_counter f_counter+1];
       end
       f_counter = f_counter +1;
    end
 
    cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/posmv/all/  % EDIT
 
    t_pos=[];
    lat=[];
    lon=[];
    altitude=[];
    heading=[];
    accel_lon =[];
    accel_tran = [];
    accel_down = [];
    rate_lon = [];
    rate_tran = [];
    rate_down = [];
    x_vel = [];
    y_vel = [];
    z_vel = [];
    roll = [];
    pitch = [];
    heave = [];
    altitude_err = [];
    lat_err = [];
    lon_err = [];
    roll_err = [];
    pitch_err = [];
    heading_err = [];
    t_err = [];

    % load posmv files
    % load an additional file which is of the next hour
    for i=1:2
        
        fname=d_pos(pos_file_index(i)).name;
        load(fname)

        t_pos=vertcat(t_pos,posmv.t);
        lat=vertcat(lat,posmv.lat);
        lon=vertcat(lon,posmv.lon);
        altitude=vertcat(altitude,posmv.altitude);
        heading=vertcat(heading,posmv.heading);
        accel_lon =vertcat(accel_lon,posmv.accel_lon);
        accel_tran =vertcat(accel_tran,posmv.accel_tran);
        accel_down =vertcat(accel_down,posmv.accel_down);
        rate_lon =vertcat(rate_lon,posmv.rate_lon);
        rate_tran =vertcat(rate_tran,posmv.rate_tran);
        rate_down =vertcat(rate_down,posmv.rate_down);
        x_vel = vertcat(x_vel,posmv.x_vel);
        y_vel = vertcat(y_vel,posmv.y_vel);
        z_vel = vertcat(z_vel,posmv.z_vel);
        roll = vertcat(roll,posmv.roll);
        pitch = vertcat(pitch,posmv.pitch);
        heave =  vertcat(heave,posmv.heave);
        altitude_err =  vertcat(altitude_err,posmv.altitude_err);
        lat_err =  vertcat(lat_err,posmv.lat_err);
        lon_err =  vertcat(lon_err,posmv.lon_err);
        roll_err =  vertcat(roll_err,posmv.roll_err);
        pitch_err =  vertcat(pitch_err,posmv.pitch_err);
        heading_err =  vertcat(heading_err,posmv.heading_err);
        t_err =  vertcat(t_err,posmv.t_err);
    end

    

    tx_lat=[];
    tx_lon=[];
    tx_altitude=[];
    tx_heading=[];
    tx_accel_lon =[];
    tx_accel_tran =[];
    tx_accel_down =[];
    tx_rate_lon =[];
    tx_rate_tran =[];
    tx_rate_down =[];
    tx_x_vel = [];
    tx_y_vel = [];
    tx_z_vel = [];
    tx_roll = [];
    tx_pitch = [];
    tx_heave = [];
    tx_altitude_err =   [];
    tx_lat_err =   [];
    tx_lon_err =  [];
    tx_roll_err =  [];
    tx_pitch_err =  [];
    tx_heading_err =  [];
    tx_t_err =   [];
    

    % Find the nearest neighbor
    t_offset = [];       % time offset between the nearest POS MV and the tx time
    for i=1:length(tx_time)

        [t_diff,B]=min(abs(tx_time(i)-t_pos));
        t_offset = horzcat(t_offset,t_diff);
        tx_lat=horzcat(tx_lat,lat(B));
        tx_lon=horzcat(tx_lon,lon(B));
        tx_altitude=horzcat(tx_altitude,altitude(B));
        tx_heading=horzcat(tx_heading,heading(B));
        tx_accel_lon=horzcat(tx_accel_lon,accel_lon(B));
        tx_accel_tran=horzcat(tx_accel_tran,accel_tran(B));
        tx_accel_down=horzcat(tx_accel_down,accel_down(B));
        tx_rate_lon=horzcat(tx_rate_lon,rate_lon(B));
        tx_rate_tran=horzcat(tx_rate_tran,rate_tran(B));
        tx_rate_down=horzcat(tx_rate_down,rate_down(B));
        tx_x_vel = horzcat(tx_x_vel,x_vel(B));
        tx_y_vel = horzcat(tx_y_vel,y_vel(B));
        tx_z_vel = horzcat(tx_z_vel,z_vel(B));
        tx_roll = horzcat(tx_roll,roll(B));
        tx_pitch = horzcat(tx_pitch,pitch(B));
        tx_heave = horzcat(tx_heave,heave(B));
        % errors are reported at 1 Hz unlike other variables
        tx_altitude_err =   horzcat(tx_altitude_err,altitude_err(round(B/10)));
        tx_lat_err =  horzcat(tx_lat_err,lat_err((round(B/10))));
        tx_lon_err = horzcat(tx_lon_err,lon_err((round(B/10))));
        tx_roll_err = horzcat(tx_roll_err,roll_err((round(B/10))));
        tx_pitch_err = horzcat(tx_pitch_err,pitch_err((round(B/10))));
        tx_heading_err = horzcat(tx_heading_err ,heading_err((round(B/10))));
        tx_t_err =  horzcat(tx_t_err,t_err((round(B/10))));

    end


    %Sort by time of TX
    [a,b]=sort(tx_time);
    tx_time=tx_time(b);
    tx_lat=tx_lat(b);
    tx_lon=tx_lon(b);
    tx_altitude=tx_altitude(b);
    tx_heading=tx_heading(b);
    tx_accel_lon =tx_accel_lon(b); 
    tx_accel_tran =tx_accel_tran(b);
    tx_accel_down =tx_accel_down(b);
    tx_rate_lon =tx_rate_lon(b);
    tx_rate_tran = tx_rate_tran(b);
    tx_rate_down = tx_rate_down(b);
    tx_x_vel = tx_x_vel(b);
    tx_y_vel = tx_y_vel(b);
    tx_z_vel =  tx_z_vel(b);
    tx_roll = tx_roll(b);
    tx_pitch = tx_pitch(b);
    tx_heave = tx_heave(b);
    tx_altitude_err = tx_altitude_err(b);
    tx_lat_err = tx_lat_err(b);
    tx_lon_err =  tx_lon_err(b);
    tx_roll_err = tx_roll_err(b);
    tx_pitch_err = tx_pitch_err(b);
    tx_heading_err = tx_heading_err(b);
    tx_t_err = tx_t_err(b);
    
    
    
    

    tx_data.t=tx_time;
    tx_data.lat=tx_lat;
    tx_data.lon=tx_lon;
    tx_data.altitude=tx_altitude;
    tx_data.heading=tx_heading;
    tx_data.accel_lon = tx_accel_lon;
    tx_data.accel_tran = tx_accel_tran;
    tx_data.accel_down = tx_accel_down;
    tx_data.rate_lon = tx_rate_lon;
    tx_data.rate_tran = tx_rate_tran;
    tx_data.rate_down = tx_rate_down;
    tx_data.rate_lon = tx_rate_lon;
    tx_data.rate_tran = tx_rate_tran;
    tx_data.rate_down = tx_rate_down;
    tx_data.x_vel = tx_x_vel;
    tx_data.y_vel = tx_y_vel;
    tx_data.z_vel = tx_z_vel;
    tx_data.roll = tx_roll;
    tx_data.pitch = tx_pitch;
    tx_data.heave = tx_heave;
    tx_data.t_offset  = t_offset;
    tx_data.altitude_err =  tx_altitude_err;
    tx_data.lat_err = tx_lat_err ;
    tx_data.lon_err = tx_lon_err ;
    tx_data.roll_err = tx_roll_err;
    tx_data.pitch_err = tx_pitch_err ;
    tx_data.heading_err = tx_heading_err; 
    tx_data.t_err = tx_t_err;
    
    %%% EDIT save directory
    cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/tx_file/all/3
    save_fname = sprintf("tx_data_2018_10_%.0f_%.0f",pos_date(pos_file_index(1)),pos_hr(pos_file_index(1)));
    if pos_hr(pos_file_index(1)) <10
       save_fname = sprintf("tx_data_2018_10_%.0f_0%.0f",pos_date(pos_file_index(1)),pos_hr(pos_file_index(1))) ;
    end
    save_fname
    save(save_fname,'tx_data')
end



