% Find Estimated Reception Time with associated information
% 1. estimated arrival time of individual transmissions
% 2. actual arrival time of individual transmissions
% 3. SNR of individual transmission
% 4.Output hourly data
% RTX data/ lat/lon/altitude of the tx transducer
%%%%%%%%%%%%%%%%%
clearvars
close all

%ACO LAT/LON
ACO_lat=22.738772;                  % Aug 2018
ACO_lon=-158.006186;                % Aug 2018

%ACO_lat= 22.738783;                  % Sep 2018
%ACO_lon= -158.00619898;                % Sep 2018

%ACO_lat=22.7387938;                    % Sep 2018 z suppress
%ACO_lon=-158.00620893;                % Sep 2018 z suppress

day = 27:30;            %  Edit
start_hour = 3;         % Edit
end_hour = 14;          % Edit
%%  Bandpass filter for HEM 2000 Hz to 6000 Hz
fl = 2000;
fh = 6000;
Fs = 24000;
trans_band = 500;  %Hz
fcut = [fl fl fh fh] + [-1 1 -1 1]*trans_band;
mags = [0 1 0];
devs = [0.01 0.05 0.01];
[n,Wn,beta,ftype] = kaiserord(fcut,mags,devs,Fs);
b= fir1(n,Wn);
[h,freq] = freqz(b,1,16000,Fs);


%% Actaul Arrival
% create a set of file names
au_fname = [];
now_hour = start_hour
now_day = day(1)
while true
    cd("/Volumes/ACO_RAP_2/RAP/Oct2018Cruise/wav_data/HEM/"+string(now_day)+"_Oct")
    while now_hour <= 23
        if now_hour < 10
            fname_d = "1810"+string(now_day)+"-0"+string(now_hour);
        else
            fname_d =  "1810"+string(now_day)+"-"+string(now_hour);
        end
        d = dir(fname_d+"*");
        au_fname = vertcat(au_fname,d.name);
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
now_hour = start_hour;
now_day = day(1);
counter = 1;

% Find Actual Arrival Times
while true
    
    % load tx file and calculate estimated arrival times
    [tx_t,tx_lat,tx_lon,tx_heading,x_dist,est_arrival] = posmv_tx_load(now_day,now_hour);
    % minutes
    m = 1;
    % hourly file
    act_arrival=[];
    SNR=[];
    estimate = [];
    for m = 1:12
        cd("/Volumes/ACO_RAP_2/RAP/Oct2018Cruise/wav_data/HEM/"+string(now_day)+"_Oct")
        %Load Audio Data
        mat_name = au_fname(counter,:);
        wav_name = au_fname(counter+1,:);
        % Timing audio data 
        [y,t_date] = hyd_audio_prep(mat_name,wav_name);
        datestr(t_date(1))
        % Filtering
        yf = filter(b,1,y); 
        delay = round(mean(grpdelay(b,1,6000,Fs)));      % group delay
        yf = yf(delay+1:end);                       % shift output signal forward in time
        yf(end+1:end+1+delay-1) = zeros(delay,1); 


        
        % Cut off est_arrival out of the time frame
        if m == 1
            est = est_arrival(find(est_arrival > t_date(1)));
        end
        
        %%%%%%Cross Correlation%%%%%%

        [demod,demod_pos,demod_max,demod_snr]=ACO_cross_correlation_sample(yf);
         % Arrival timestamps
         arrivals=t_date(demod_pos);
         
        
         
         % Match the calculated actual arrival time with the estimated time calculated in the previous section
        for ii=1:length(arrivals)
            [arrival_diff,arrival_pos]=min(abs(est-arrivals(ii)));
            arrival_diff*(3600*24)
          % if the discrepency is less than 50 ms 
            if abs(arrival_diff*(3600*24))<0.05
                estimate(end+1) = est(arrival_pos);
                act_arrival(end+1)=arrivals(ii);   % concatenate that timestamp and corresponding snr
                SNR(end+1)=demod_snr(ii);         
            end


        end

        %fprintf('Results \n len of demod_snr = %d \n len of arrivals = %d \n ',length(demod_snr),length(arrivals));
        %
        %Remove est_arrival without a match of act_arrival (no receptions made)
        k  = 1; % counter for the estimates
        while true
            
           if k >= length(estimate)+1
               break;
           end
           if min(abs(estimate(k) - act_arrival)*3600*24) > 0.05
               estimate(k) = [];
           else
               k = k+1;
           end

        end
       %}
          m = m+1;
         counter = counter+2;
        
    end
   
    % Save
    cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/RTX/rx_file/corrected/ideal_pulse/t_offset
    %Save Variables
    rtx_rx_data.est_arrival=estimate;
    rtx_rx_data.act_arrival=act_arrival;
    rtx_rx_data.SNR=SNR;
    sname = mat_name(1:end-8);
    sname = ['rtx_rx_data' '_20' sname(1:2) '_' sname(3:4) '_' sname(5:6) '_' sname(8:9)]
    save(sname,'rtx_rx_data')

    
   if now_hour < 23  
        now_hour = now_hour+1
   else
       now_hour = 0;
       now_day = now_day+1
   end
   
     if (now_hour > end_hour)& (now_day == day(end))
        break;
     end
end


%%
function [tx_t,tx_lat,tx_lon,tx_heading,x_dist,est_arrival] = posmv_tx_load(day,hour)
    %% Load POS MV TX/RX files
    cd('/Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/RTX/tx_file/corrected/t_offset')
    % create a set of file names
    %ACO LAT/LON
    ACO_lat=22.738772;                  % Aug 2018
    ACO_lon=-158.006186;                % Aug 2018
    fname = [];
    now_hour = hour;
    now_day = day;
    % load 2 files
    if now_day ==27 & now_hour == 3
        fname = "rtx_tx_data_2018_10_27_03";
    else
       for p = 1:2
            if p ==1 & now_hour == 0
                fname_d = "rtx_tx_data_2018_10_"+string(now_day-1)+"_"+string(23);
            elseif p==1 & now_hour == 10
                fname_d = "rtx_tx_data_2018_10_"+string(now_day)+"_0"+string(9);
            elseif now_hour < 10
                fname_d = "rtx_tx_data_2018_10_"+string(now_day)+"_0"+string(now_hour-2+p);
            else
                fname_d = "rtx_tx_data_2018_10_"+string(now_day)+"_"+string(now_hour-2+p);
            end
            fname = vertcat(fname,fname_d );
       end
    end
   fname
    %Load TX data
    tx_t = [];
    tx_lat = [];
    tx_lon = [];
    tx_altitude = [];
    tx_heading = [];

   for p = 1:length(fname)
        load(fname(p))
        tx_t = horzcat(tx_t,rtx_tx_data.t);      %TX time                                  
        tx_lat = horzcat(tx_lat,rtx_tx_data.lat);   %TX lat
        tx_lon = horzcat(tx_lon,rtx_tx_data.lon);   
        tx_altitude = horzcat(tx_altitude,rtx_tx_data.altitude);   
        tx_heading = horzcat(tx_heading,rtx_tx_data.heading);
    
   end

    % offset t_tx for HEM time offset
    % correct for the reception time
        date_mark1= "20181028 00:50";
        date_mark2= "20181029 00:50";
        date_mark3= "20181030 00:50";
        date_mark1 = datenum(date_mark1,'yyyymmdd HH:MM');
        date_mark2 = datenum(date_mark2,'yyyymmdd HH:MM');
        date_mark3 = datenum(date_mark3,'yyyymmdd HH:MM');

        for p = 1:length(tx_t)
            if tx_t(p) <= date_mark1
                tx_t(p) = tx_t(p) +5/(3600*24);
                
            elseif (date_mark1 <= tx_t(p))&(tx_t(p) <= date_mark2)
                tx_t(p) = tx_t(p)+6/(3600*24);
                
            elseif (date_mark2 <= tx_t(p))&(tx_t(p) <= date_mark3)
                tx_t(p) = tx_t(p)+7/(3600*24);
                
            elseif  date_mark3 <= tx_t(p)
                tx_t(p) = tx_t(p)+8/(3600*24);
                
            end
        end

    %Boat distance from ACO
    for i=1:length(tx_lat)
        x_dist(i)=dist([ACO_lat tx_lat(i)],[ACO_lon tx_lon(i)]);
    end

    
    %Estimate travel time based on CTD cast

    for i=1:length(x_dist)
    [~,~,~,~,~,~,~,est_tt(i),~]=ray_trace_w_curvature_v3(x_dist(i),tx_altitude(i));
    end

    %Estimate arrival time
    est_arrival=tx_t+(est_tt./(3600*24));

    tx_t = tx_t';
    tx_lat = tx_lat';
    tx_lon = tx_lon';
    tx_altitude = tx_altitude';
    tx_heading = tx_heading';
    x_dist = x_dist';
    est_arrival =est_arrival';

    
    
end










