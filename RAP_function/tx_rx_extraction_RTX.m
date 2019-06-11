function [tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,act_arrival,est_arrival] = tx_rx_extraction_RTX(day,start_hour,end_hour,hydrophone)
% return time-corrected HEM or icListen transmission and recpetion
% information
% Input: 1. month vector 2.day vector 3.start hour/ end hour 4. HEM or icListen
% PS: month should be in 'mm' format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% HEM LAT/LON
HEM_lat=22.738772;                  % June 2017
HEM_lon=-158.006186;                % June 2017

% icListen LAT/LON
icListen_lat = 22.739153;
icListen_lon = -158.0061254;

% create a set of file names
fname = [];
now_hour = start_hour;
now_day = day(1);



%% 2 Load POS MV Tx data
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/RTX/tx_file/corrected/t_offset
% create a set of file names
fname = [];

now_hour = start_hour;
now_day = day(1);
while true
    while now_hour <= 23
        if now_hour < 10
            fname_d = "rtx_tx_data_2018_10_"+string(now_day)+"_0"+string(now_hour);
        else
            fname_d = "rtx_tx_data_2018_10_"+string(now_day)+"_"+string(now_hour);
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
for p = 1:length(fname)
    load(fname(p))
    tx_t = horzcat(tx_t,rtx_tx_data.t);      %TX time                                  
    tx_lat = horzcat(tx_lat,rtx_tx_data.lat);   %TX lat
    tx_lon = horzcat(tx_lon,rtx_tx_data.lon);
    tx_altitude = horzcat(tx_altitude,rtx_tx_data.altitude);
    tx_heading = horzcat(tx_heading,rtx_tx_data.heading);
    tx_xvel = horzcat(tx_xvel,rtx_tx_data.xvel);
end

switch string(hydrophone)
%%%%%%%%%%%%%%%%%%%% HEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case  "HEM"
        % Load  HEM rx files
        % create a set of file names
        cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/RTX/rx_file/corrected/ideal_pulse/t_offset
        fname = [];
        now_hour = start_hour;
        now_day = day(1);
        while true
            while now_hour <= 23
                if now_hour < 10
                    fname_d = "rtx_rx_data_2018_10_"+string(now_day)+"_0"+string(now_hour);
                else
                    fname_d = "rtx_rx_data_2018_10_"+string(now_day)+"_"+string(now_hour);
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
            load(fname(i)+".mat")
            est_arrival = horzcat(est_arrival,rtx_rx_data.est_arrival);
            act_arrival = horzcat(act_arrival ,rtx_rx_data.act_arrival);
            SNR = horzcat(SNR ,rtx_rx_data.SNR);

        end
        %%%%%%%%%%%%%%%%% Clean data %%%%%%%%%%%%%%%%%%%%%%%%
       
        %%%%% HEM %%%%%%
        rm_ind = [];


        % correct for the reception time
        date_mark1= "20181028 01:00";
        date_mark2= "20181029 01:00";
        date_mark3= "20181030 01:00";
        date_mark1 = datenum(date_mark1,'yyyymmdd HH:MM');
        date_mark2 = datenum(date_mark2,'yyyymmdd HH:MM');
        date_mark3 = datenum(date_mark3,'yyyymmdd HH:MM');

        for p = 1:length(act_arrival)
            if act_arrival(p) <= date_mark1
                act_arrival(p) = act_arrival(p) -5/(3600*24);
                est_arrival(p) = est_arrival(p) -5/(3600*24);
            elseif (date_mark1 <= act_arrival(p))&(act_arrival(p) <= date_mark2)
                act_arrival(p) = act_arrival(p)-6/(3600*24);
                est_arrival(p) = est_arrival(p) -6/(3600*24);
            elseif (date_mark2 <= act_arrival(p))&(act_arrival(p) <= date_mark3)
                act_arrival(p) = act_arrival(p)-7/(3600*24);
                est_arrival(p) = est_arrival(p) -7/(3600*24);
            elseif  date_mark3 <= act_arrival(p)
                act_arrival(p) = act_arrival(p)-8/(3600*24);
                est_arrival(p) = est_arrival(p) -8/(3600*24);
            end
        end



        for l = 1:length(tx_t)
            dum_arrival = act_arrival;
            dum_arrival(find(dum_arrival < tx_t(l))) = [];          % shift forward in time 6 seconds
            h_ind = find((tx_t(l) - dum_arrival)*3600*24 > -20);   % HEM is behind UTC for ~ 5-8sec
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
        cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/rx_file/icListen/all/
        fname = [];
        now_hour = start_hour;
        now_day = day(1);

        while true
            while now_hour <= 23
                if now_hour < 10
                    fname_d = "rx_data_2018_10_"+string(now_day)+"_0"+string(now_hour)+"_icListen";
                else
                    fname_d = "rx_data_2018_10_"+string(now_day)+"_"+string(now_hour)+"_icListen";
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

            % correct for the reception time
            date_mark= "20181030 05:55";
            date_mark = datenum(date_mark,'yyyymmdd HH:MM');

            for p = 1:length(act_arrival)
                if act_arrival(p) >= date_mark
                    act_arrival(p) = act_arrival(p) -2/(3600*24);
                    est_arrival(p) = est_arrival(p) -2/(3600*24);
                else
                    act_arrival(p) = act_arrival(p)+6/(3600*24);
                    est_arrival(p) = est_arrival(p) +6/(3600*24);
                end
            end

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