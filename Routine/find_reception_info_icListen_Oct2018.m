% Produce Rx files
% Setup
% 1. Days and Hours of files on October Cruise
% 2. Tx File Directory
% 3. Audio File Directpry
% 4. Rx File Directory
% OUTPUT
% 1. estimated arrival times of individual transmissions
% 2. actual arrival time of individual transmissions
% 3. SNR of individual transmission
%%%%%%%%%%%%%%%%%
clearvars
close all
% Days and Hours
% ex. 14:00 UTC of 27 Oct 2018 to 05:00 UTC of 28 Oct 2018
% day = 27:28
% hour = 14:5

day = 27:30;            %  Edit
hour = 3:14;            %  Edit

%%  Bandpass filter for icListen 2000 Hz to 6000 Hz
fl = 2000;
fh = 6000;
Fs = 32000;
trans_band = 500;  %Hz
fcut = [fl fl fh fh] + [-1 1 -1 1]*trans_band;
mags = [0 1 0];
devs = [0.01 0.05 0.01];
[n,Wn,beta,ftype] = kaiserord(fcut,mags,devs,Fs);
b= fir1(n,Wn);
[h,freq] = freqz(b,1,16000,Fs);

%% Actaul Arrival
au_fname = [];
now_hour = hour(1);
now_day = day(1);   
while true
    cd("/Volumes/ACO_RAP_2/RAP/Oct2018Cruise/wav_data/icListen/"+string(now_day)+"_Oct")
    while now_hour <= 23
        if now_hour < 10
            fname_d = "SBW1391_201810"+string(now_day)+"_0"+string(now_hour);
        else
            fname_d =  "SBW1391_201810"+string(now_day)+"_"+string(now_hour);
        end
        d = dir(fname_d+"*");
        au_fname = vertcat(au_fname,d.name);
       now_hour = now_hour+1;
       if (now_hour > hour(end))& (now_day == day(end))
            break;
       end
    end
        if (now_hour > hour(end))& (now_day == day(end))
            break;
       end
        now_hour = 0;
        now_day = now_day+1;
end

%% Find Actual Arrival Times
now_day = day(1);
now_hour = hour(1)
counter = 1;
[len,~] = size(au_fname);

while true
    
    % Load Tx File and Calculate Estimated Travel Time of this 1 hour
    [tx_t,tx_lat,tx_lon,tx_heading,x_dist,est_arrival] = posmv_tx_load(now_day,now_hour);
    
    % hourly file
    act_arrival=[];
    SNR=[];
    estimate = [];
    
    % store the current hour
    m = 1;
    wav_name = au_fname(counter,:);
    current_hour = str2num(wav_name(18:19))
    
    % Loop over minutes in the hour
    while true
    
        % Audio File Directory
        cd("/Volumes/ACO_RAP_2/RAP/Oct2018Cruise/wav_data/icListen/"+string(now_day)+"_Oct")    % EDIT
        %Load Audio Data
        wav_name = au_fname(counter,:)
        
        % Timing audio data 
        [y,t_date] = icListen_audio_prep(wav_name);
        time_offset = time_correction(t_date(1));       % load time offset
        t_date = t_date - time_offset;                  % icListen time offset from UTC (ahead of UTC)
        datestr(t_date(1))
        
        % Filtering
        yf = filter(b,1,y); 
        delay = round(mean(grpdelay(b,1,6000,Fs)));      % group delay
        yf = yf(delay+1:end);                            % shift output signal forward in time
        yf(end+1:end+1+delay-1) = zeros(delay,1); 

        % drop out est_arrival outside this 5-min time frame
        if m == 1
            est = est_arrival(find(est_arrival > t_date(1)));
            datestr(est(1))
        end
        
        %%%%%%Cross Correlation%%%%%%

        [demod,demod_pos,demod_max,demod_snr]=ACO_cross_correlation_ideal_icListen(yf);
         % Arrival timestamps
         arrivals=t_date(demod_pos);
         
         % for debugging
         %{
         figure(1)
         clf
         plot(t_date,demod)
         datetick('x')
         grid on
         hold on
         scatter(t_date(demod_pos),demod_max,'r')
         axis tight
         
         figure(2)
         clf
         plot(t_date,y)
          hold on
         plot(t_date,yf)
         grid on
         datetick('x')
         axis tight
         %}
         
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
        
        %Remove est_arrival without an act_arrival match (no receptions at the hydrophone end)
        
        k  = 1; % counter for the estimates
        if ~isempty(arrivals) & ~isempty(estimate)
            while true
               if min(abs(estimate(k) - act_arrival)*3600*24) > 0.05
                   estimate(k) = [];
               else
                   k = k+1;
               end

               if k >= length(estimate)+1
                   break;
               end
            end
        end
        
        counter = counter+1;
        m = m+1
        % check if the file jumps to the next hour
        if counter > len
            break;
        end
        nx_fname = au_fname(counter,:);
        chk_hour =  str2double(nx_fname(18:19))
       
        if chk_hour ~= current_hour
            break;
        end
        
    end
    
%% Save hourly file
% 3. Rx File Directory
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/rx_file/icListen/final/Oct/original_position % EDIT
% Save Variables
rx_data.est_arrival=estimate;
rx_data.act_arrival=act_arrival;
rx_data.SNR=SNR;
sname =wav_name(1:end-6);
sname = ['rx_data' '_20' sname(11:12) '_' sname(13:14) '_' sname(15:16) '_' sname(18:19) '_icListen']
save(sname,'rx_data')

% update day and hour for the next hour file
if now_hour < 23
    now_hour = now_hour+1
else
    now_hour = 0;
    now_day = now_day+1
end

if (now_hour > hour(end))& (now_day == day(end))
    break;
end


end



%% Functions
%%%%
 function [tx_t,tx_lat,tx_lon,tx_heading,x_dist,est_arrival] = posmv_tx_load(now_day,now_hour)
 % setup
 % 1. Hydrophone Position and Depth
 %%% icListen LAT/LON
% icListen_lat = 22.738894;                  % original (HEM)
% icListen_lon = -158.006009;                % original (HEM)

% Grant Numbers: 
% icListen_lat = 22.7390046;          % 26 m  north of HEM (June 2017 Position)   
% icListen_lon =  -158.00582401;      % 37 m east of HEM (June 2017 Position)   

% Vincent's inversion result
icListen_lat = 22.739153;                   % June 2018 north of HEM  42.22 m nort of HEM(June 2017 Position)   
icListen_lon = -158.0061254;                % June 2018 east of HEM  5.94 m east of HEM(June 2017 Position)   

% icListen_lat = 22.7391066;                   % Oct 2018  from first iteraion
% icListen_lon = -158.00610724;                % Oct 2018 from first iteraion

% icListen_lat=22.73912734;                  % March 2019 #2 
% icListen_lon=-158.006111197;                % March 2019 

% icListen_lat = 22.73911569;                  % March 2019 #3
% icListen_lon = -158.006106601;                % March 2019

%%% icListen Depth
% icListen_depth = -4736.226+2.32;            % June 2018 (from HEM)
% icListen_depth = -4734.646+2.32;            % March 2019 # 2
icListen_depth = -4729.92+2.32 +1.75;          % original depth ellipsoid height+ 1.75 m higher than HEM
% icListen_depth = -4733.24+2.32;                    %  Oct 2018 from first iteraion
 
% 2. Tx File Directory
cd("/Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/tx_file/all/3")  %% EDIT 

% 3. edit CTD file in the ray traing code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load TX data
tx_t = [];
tx_lat = [];
tx_lon = [];
tx_altitude = [];
 tx_heading  = [];
fname =[];
% load 2 files
    if now_day ==27 & now_hour == 3
        fname = "tx_data_2018_10_27_03";
    else
       for p = 1:2
            if p ==1 & now_hour == 0
                fname_d = "tx_data_2018_10_"+string(now_day-1)+"_"+string(23);
            elseif p==1 & now_hour == 10
                fname_d = "tx_data_2018_10_"+string(now_day)+"_0"+string(9);
            elseif now_hour < 10
                fname_d = "tx_data_2018_10_"+string(now_day)+"_0"+string(now_hour-2+p);
            else
                fname_d = "tx_data_2018_10_"+string(now_day)+"_"+string(now_hour-2+p);
            end
            fname = vertcat(fname,fname_d );
       end
    end
   fname
   
for p = 1:length(fname)
    load(fname(p))
    tx_t = horzcat(tx_t,tx_data.t);      %TX time                                  
    tx_lat = horzcat(tx_lat,tx_data.lat);   %TX lat
    tx_lon = horzcat(tx_lon,tx_data.lon);   
    tx_altitude = horzcat(tx_altitude,tx_data.altitude);   
    tx_heading = horzcat(tx_heading,tx_data.heading);
    
end


x_dist = [];
est_tt = [];

% Geodic Length
for ii=1:length(tx_lat)
    x_dist(ii) = distance(tx_lat(ii),tx_lon(ii),icListen_lat,icListen_lon,referenceEllipsoid('WGS84'));
end

% Estimated travel time based on CTD cast
for ii=1:length(x_dist)
    azmth(ii) = azimuth(tx_lat(ii),tx_lon(ii),icListen_lat,icListen_lon);
    [~,~,~,~,~,~,~,~,est_tt(ii),~,~] = ray_trace_w_earth_flattening(x_dist(ii),tx_altitude(ii),tx_lat(ii),azmth(ii),icListen_lat,icListen_lon,icListen_depth);
end

% Estimated arrival time
est_arrival = tx_t+(est_tt./(3600*24));

 end
 
 function time_offset = time_correction(time)
 % return time offset for timestamp correction
 date_mark1=  "20181030 05:55";
 
 date_mark1 = datenum(date_mark1,'yyyymmdd HH:MM');
 
 if time <= date_mark1
     time_offset =  -7/(3600*24); % 7 secbehind
     
 else
     time_offset = 1/(3600*24);
     
 end
 
 end