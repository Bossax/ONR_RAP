function [tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,x_err,y_err,z_err,act_arrival,est_arrival,SNR] = tx_rx_extraction_June(day,start_hour,end_hour,hydrophone)
% return time-corrected HEM or icListen transmission and recpetion
% information
% Input: 1. month vector 2.day vector 3.start hour/ end hour 4. HEM or icListen
% PS: month should be in 'mm' format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% HEM LAT/LON
HEM_lat=22.738772;                  % Aug 2018
HEM_lon=-158.006186;                % Aug 2018

% icListen LAT/LON
icListen_lat = 22.739153;
icListen_lon = -158.0061254;

% create a set of file names
fname = [];
now_hour = start_hour;
now_day = day(1);



%% 2 Load POS MV Tx data
cd /Volumes/ACO_RAP_2/RAP/June2018Cruise/Tx_Rx_Output/tx_file/all
% create a set of file names
fname = [];

now_hour = start_hour;
now_day = day(1);
while true
    while now_hour <= 23
        if now_hour < 10
            fname_d = "tx_data_2018_06_"+string(now_day)+"_0"+string(now_hour);
        else
            fname_d = "tx_data_2018_06_"+string(now_day)+"_"+string(now_hour);
        end
       fname = vertcat(fname,fname_d);
       now_hour = now_hour+1;
       if (now_hour > end_hour)& (now_day == day(end))
            break;
       end
    end
        if (now_hour > end_hour)& (now_day == day(end))
            break;
       end
        now_hour = 0;
        now_day = now_day+1;
end

%Load TX data
tx_t = [];
tx_lat = [];
tx_lon = [];
tx_altitude = [];
tx_heading = [];
tx_xvel = [];
x_err = [];
y_err = [];
z_err = [];
for p = 1:length(fname)
    load(fname(p))
    tx_t = horzcat(tx_t,tx_data.t);      %TX time                                  
    tx_lat = horzcat(tx_lat,tx_data.lat);   %TX lat
    tx_lon = horzcat(tx_lon,tx_data.lon);
    tx_altitude = horzcat(tx_altitude,tx_data.altitude);
    tx_heading = horzcat(tx_heading,tx_data.heading);
    tx_xvel = horzcat(tx_xvel,tx_data.x_vel);
    x_err = horzcat(x_err,tx_data.lon_err);
    y_err = horzcat(y_err,tx_data.lat_err);
    z_err = horzcat(z_err,tx_data.altitude_err);
end

switch string(hydrophone)
%%%%%%%%%%%%%%%%%%%% HEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case  "HEM"
        % Load  HEM rx files
        % create a set of file names
        cd /Users/testuser/Documents/ONR_RAP/Data/Tx_Rx_Output/June2018/rx_file/original_depth
        fname = [];
        now_hour = start_hour;
        now_day = day(1);
        while true
            while now_hour <= 23
                if now_hour < 10
                    fname_d = "rx_data_2018_06_"+string(now_day)+"_0"+string(now_hour);
                else
                    fname_d = "rx_data_2018_06_"+string(now_day)+"_"+string(now_hour);
                end
               fname = vertcat(fname,fname_d);
               now_hour = now_hour+1;
               if (now_hour > end_hour)& (now_day == day(end))
                    break;
               end
            end
                if (now_hour > end_hour)& (now_day == day(end))
                    break;
               end
                now_hour = 0;
                now_day = now_day+1;
        end

        est_arrival =[];
        act_arrival = [];
        SNR = [];
        for i = 1:length(fname)
            try
            load(fname(i)+".mat")
            est_arrival = horzcat(est_arrival,rx_data.est_arrival);
            act_arrival = horzcat(act_arrival ,rx_data.act_arrival);
            SNR = horzcat(SNR ,rx_data.SNR);
            catch
                
            end

        end
        %%%%%%%%%%%%%%%%% Clean data %%%%%%%%%%%%%%%%%%%%%%%%
       
        %%%%% HEM %%%%%%
        rm_ind = [];

        for l = 1:length(tx_t)
            dum_arrival = act_arrival;
            dum_arrival(find(dum_arrival < tx_t(l))) = [];          
            h_ind = find((tx_t(l) - dum_arrival)*3600*24 > -18);   
            if isempty(h_ind)
                rm_ind(end+1) = l;
            end
        end
        keep_ind = setdiff(1:length(tx_t),rm_ind);

        tx_t =tx_t(keep_ind);
        tx_lat = tx_lat(keep_ind);
        tx_lon = tx_lon(keep_ind);
        tx_altitude  = tx_altitude(keep_ind);
        tx_heading = tx_heading(keep_ind);
        tx_xvel = tx_xvel(keep_ind);
        x_err = x_err(keep_ind);
        y_err = y_err(keep_ind);
        z_err = z_err(keep_ind);

        rm_ind = [];
        for l = 1:length(act_arrival)
            dum_tx_t = tx_t;
             dum_tx_t(find( dum_tx_t > act_arrival(l))) = [];
            % if greater than 20 sec = no reception was picked up 
            h_ind = find((act_arrival(l) - dum_tx_t)*3600*24 < 20);
            if isempty(h_ind)
                rm_ind(end+1) = l;
            end
        end
        act_arrival(rm_ind) = [];
        est_arrival(rm_ind) = [];
        SNR(rm_ind) = [];

        %Boat distance from the hydrophone
        for i=1:length(tx_lat)
            range(i)=dist([HEM_lat tx_lat(i)],[HEM_lon tx_lon(i)]);
        end
        range = range/1000;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    case "icListen"
        % Load  icListen rx files
        % create a set of file names    
        cd /Volumes/ACO_RAP_2/RAP/June2018Cruise/Tx_Rx_Output/rx_file/icListen/all/
        fname = [];
        now_hour = start_hour;
        now_day = day(1);

        while true
            while now_hour <= 23
                if now_hour < 10
                    fname_d = "rx_data_2018_06_"+string(now_day)+"_0"+string(now_hour)+"_icListen";
                else
                    fname_d = "rx_data_2018_06_"+string(now_day)+"_"+string(now_hour)+"_icListen";
                end
               fname = vertcat(fname,fname_d);
               now_hour = now_hour+1;
               if (now_hour > end_hour)& (now_day == day(end))
                    break;
               end
            end
                if (now_hour > end_hour)& (now_day == day(end))
                    break;
               end
                now_hour = 0;
                now_day = now_day+1;
        end

        est_arrival =[];
        act_arrival= [];
        SNR= [];
        for i = 1:length(fname)
            load(fname(i)+".mat")
            est_arrival  = horzcat(est_arrival ,rx_data.est_arrival);
            act_arrival  = horzcat(act_arrival  ,rx_data.act_arrival);
            SNR  = horzcat(SNR  ,rx_data.SNR);

        end
        
    
            %%%%%%%%%% Clean Data %%%%%%%%%%%%%%%%%%%%
            %%%%% icListen %%%%%%%
            rm_ind = [];

            for l = 1:length(tx_t)
                dum_arrival = act_arrival;
                dum_arrival(find(dum_arrival < tx_t(l))) = [];
                % if greater than 20 sec = no reception was picked up 
                h_ind = find((tx_t(l) - dum_arrival)*3600*24 > -20);
                if isempty(h_ind)
                    rm_ind(end+1) = l;
                end
            end

            keep_ind = setdiff(1:length(tx_t),rm_ind);
            tx_t=tx_t(keep_ind);
            tx_lat= tx_lat(keep_ind);
            tx_lon= tx_lon(keep_ind);
            tx_altitude= tx_altitude(keep_ind);
            tx_heading= tx_heading(keep_ind);
            tx_xvel= tx_xvel(keep_ind);
            x_err = x_err(keep_ind);
            y_err = y_err(keep_ind);
            z_err = z_err(keep_ind);

            rm_ind = [];
            for l = 1:length(act_arrival)
                dum_tx_t = tx_t;
                 dum_tx_t(find( dum_tx_t > act_arrival(l))) = [];
                % if greater than 20 sec = no reception was picked up 
                h_ind = find((act_arrival(l) - dum_tx_t)*3600*24 < 20);
                if isempty(h_ind)
                    rm_ind(end+1) = l;
                end
            end

            act_arrival(rm_ind) = [];
            est_arrival(rm_ind) = [];
            SNR(rm_ind) = [];

            act_trtime = (act_arrival-tx_t)*3600*24;
            est_trtime= (est_arrival-tx_t)*3600*24;
            
            % Vessel Surface Distance
            
            for i=1:length(tx_lat)
                range(i)=dist([icListen_lat tx_lat(i)],[icListen_lon tx_lon(i)]);
            end
            range = range/1000;
    otherwise
        fprintf("Hydrophone input must be either 'HEM' or 'icListen'\n ")
    end





end