function [tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,x_err,y_err,z_err,act_arrival,est_arrival,SNR] = tx_rx_extraction_Oct(day,start_hour,end_hour,hydrophone)
% return time-corrected HEM or icListen transmission and recpetion
% information
% Hydrophone times are corrected for their offsets from the UTC time
% Input: 1. month vector 2.day vector 3.start hour/ end hour 4. HEM or icListen
% PS: month must be in 'mm' format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% HEM LAT/LON
HEM_lat=22.738772;                  % Aug 2018
HEM_lon=-158.006186;                % Aug 2018


% icListen LAT/LON
icListen_lat = 22.739153;
icListen_lon = -158.0061254;

icListen_lat=22.73912734;                  % March 2019
icListen_lon=-158.006111197;                % March 2019

% create a set of file names
fname = [];
now_hour = start_hour;
now_day = day(1);



%% 2 Load Tx data
% EDIT  tx file directory
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/tx_file/all/3

% create a set of file names
fname = [];

now_hour = start_hour;
now_day = day(1);
while true
    while now_hour <= 23
        if now_hour < 10
            fname_d = "tx_data_2018_10_"+string(now_day)+"_0"+string(now_hour);
        else
            fname_d = "tx_data_2018_10_"+string(now_day)+"_"+string(now_hour);
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
        
      % rx file directory EDIT  
      cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/rx_file/HEM/final/Oct/original_depth

        fname = [];
        now_hour = start_hour;
        now_day = day(1);
        while true
            while now_hour <= 23
                if now_hour < 10
                    fname_d = "rx_data_2018_10_"+string(now_day)+"_0"+string(now_hour);
                else
                    fname_d = "rx_data_2018_10_"+string(now_day)+"_"+string(now_hour);
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
            est_arrival = horzcat(est_arrival,rx_data.est_arrival);
            act_arrival = horzcat(act_arrival ,rx_data.act_arrival);
            SNR = horzcat(SNR ,rx_data.SNR);

        end
        
         % correct for the reception time
         %{
            date_mark1= "20181028 01:00";
            date_mark2= "20181029 01:00";
            date_mark3= "20181030 01:00";
            date_mark1 = datenum(date_mark1,'yyyymmdd HH:MM');
            date_mark2 = datenum(date_mark2,'yyyymmdd HH:MM');
            date_mark3 = datenum(date_mark3,'yyyymmdd HH:MM');

            for p = 1:length(act_arrival)
                if act_arrival(p) <= date_mark1
                    act_arrival(p) = act_arrival(p) -4/(3600*24);
                    est_arrival(p) = est_arrival(p) -4/(3600*24);
                elseif (date_mark1 <= act_arrival(p))&(act_arrival(p) <= date_mark2)
                    act_arrival(p) = act_arrival(p)-5/(3600*24);
                    est_arrival(p) = est_arrival(p) -5/(3600*24);
                elseif (date_mark2 <= act_arrival(p))&(act_arrival(p) <= date_mark3)
                    act_arrival(p) = act_arrival(p)-6/(3600*24);
                    est_arrival(p) = est_arrival(p) -6/(3600*24);
                elseif  date_mark3 <= act_arrival(p)
                    act_arrival(p) = act_arrival(p)-7/(3600*24);
                    est_arrival(p) = est_arrival(p) -7/(3600*24);
                end
            end
         %}
        %%%%%%%%%%%%%%%%% Clean data %%%%%%%%%%%%%%%%%%%%%%%%
        
            %%%%% HEM %%%%%%
            keep_ind = [];
            % find the nearest arrival time (within 20 sec)
            for l = 1:length(tx_t)
                dum_arrival = act_arrival;
                dum_arrival(find(dum_arrival < tx_t(l))) = [];        
                h_ind = find((dum_arrival -tx_t(l))*3600*24 < 20);   
                if ~isempty(h_ind)
                    keep_ind(end+1) = l;
                    
                end
            end
            
            
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
        
        % adjust the geiod heights
        tx_altitude = tx_altitude+2.31;
%         gh = geoidheight(tx_lat,tx_lon+360);
%         tx_altitude = tx_altitude - gh;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    case "icListen"
        % Load  icListen rx files  EDIT        
        cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/rx_file/icListen/final/Oct/original_depth
        
        % create a set of file names 
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
            SNR  = horzcat(SNR,rx_data.SNR);

        end
        
    
            %%%%%%%%%% Clean Data %%%%%%%%%%%%%%%%%%%%
            %%%%% icListen %%%%%%%
            rm_ind = [];
            
           % correct for the reception time
%                         date_mark= "20181030 05:55";
%                         date_mark = datenum(date_mark,'yyyymmdd HH:MM');
%             
%                         for p = 1:length(act_arrival)
%                             if act_arrival(p) >= date_mark
%                                 act_arrival(p) = act_arrival(p) -1/(3600*24);
%                                 est_arrival(p) = est_arrival(p) -1/(3600*24);
%                             else
%                                 act_arrival(p) = act_arrival(p)+7/(3600*24);
%                                 est_arrival(p) = est_arrival(p) +7/(3600*24);
%                             end
%                         end

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
            
            % adjust the geiod heights
        tx_altitude = tx_altitude+2.31;
%         gh = geoidheight(tx_lat,tx_lon+360);
%         tx_altitude = tx_altitude - gh;
            
    otherwise
        fprintf("Hydrophone input must be either 'HEM' or 'icListen'\n ")
    end





end