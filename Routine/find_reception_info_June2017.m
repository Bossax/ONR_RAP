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

day = 7:12;            %  Edit
start_hour = 15            %  Edit
end_hour = 5;           % EDIT



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
now_hour = start_hour;
now_day = day(1);

while true
    cd("/Volumes/ACO_RAP_2/RAP/June2017Cruise/wav_data/" + string(now_day)+"_June")
    if now_day < 10
        str_now_day = "0"+string(now_day);
    else
        str_now_day = string(now_day);
    end

    if now_hour < 10
        fname_d = "1706"+str_now_day+"-0"+string(now_hour);
    else
        fname_d = "1706"+str_now_day+"-"+string(now_hour);
    end
    d = dir(fname_d+"*");
    au_fname = vertcat(au_fname,d.name);
     if now_hour >= 23
           now_hour = 0;
           now_day = now_day+1;
       else
           now_hour = now_hour+1;
       end

       if now_hour > end_hour && now_day == day(end)
           break;
       end
end

%% Find Actual Arrival Times
now_hour = start_hour;
now_day = day(1)
counter = 1;
[len,~] = size(au_fname);
while true

    % Load Tx File and Calculate Estimated Travel Time of this 1 hour
    [tx_t,tx_lat,tx_lon,tx_heading,x_dist,est_arrival] = posmv_tx_load(now_day,now_hour);

    % hourly file
    act_arrival=[];
    SNR=[];
    estimate = [];
    surface_range = [];
    % store the current hour
    m = 1;
    wav_name = au_fname(counter,:);
    current_hour = str2num(wav_name(8:9))

    % Loop over minutes in the hour
    while true

        % Audio File Directory
        cd("/Volumes/ACO_RAP_2/RAP/June2017Cruise/wav_data/"+string(now_day)+"_June")

        %Load Audio Data
        mat_name = au_fname(counter,:);
        wav_name = au_fname(counter+1,:);

        % Timing audio data
        [y,t_date] = hyd_audio_prep(mat_name,wav_name);
        time_offset =  -1/(3600*24);
        t_date = t_date - time_offset;      % HEM time is ahead UTC time
        datestr(t_date(1))

        % Filtering
        yf = filter(b,1,y);
        delay = round(mean(grpdelay(b,1,6000,Fs)));      % group delay
        yf = yf(delay+1:end);                       % shift output signal forward in time
        yf(end+1:end+1+delay-1) = zeros(delay,1);


        % drop out est_arrival outside this 5-min time frame
        if m == 1
            est = est_arrival(find(est_arrival > t_date(1)));
            datestr(est(1))
        end

        %%%%%%Cross Correlation%%%%%%

        [demod,demod_pos,demod_max,demod_snr] = ACO_cross_correlation_ideal(yf);
         % Arrival timestamps
         arrivals=t_date(demod_pos);
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
            disp(arrival_diff*(3600*24)*1000)
          % if the discrepency is less than 50 ms
            if abs(arrival_diff*(3600*24))<0.05
                estimate(end+1) = est(arrival_pos);
                act_arrival(end+1)=arrivals(ii);   % concatenate that timestamp and corresponding snr
                SNR(end+1)=demod_snr(ii);
                surface_range(end+1) = x_dist(ii);
            end


        end

        %fprintf('Results \n len of demod_snr = %d \n len of arrivals = %d \n ',length(demod_snr),length(arrivals));
        %
        %Remove est_arrival without a match of act_arrival (no receptions made)
        k  = 1; % counter for the estimates
        if length(estimate) ~= 0
            while true
                if min(abs(estimate(k) - act_arrival)*3600*24) > 0.05
                    estimate(k) = [];
                    surface_range(k) = [];
                    SNR(k) = [];
                else
                    k = k+1;
                end

                if k >= length(estimate)+1
                    break;
                end
            end
        end


        counter = counter+2;
        cm = m+1;
        % check if the file jumps to the next hour
        if counter > len
            break;
        end
        nx_fname = au_fname(counter,:);
        chk_hour =  str2double(nx_fname(8:9));

        if chk_hour ~= current_hour
            break;
        end

    end
%% Save hourly file
% 3. Rx File Directory
cd /Users/testuser/Documents/ONR_RAP/Data/Tx_Rx_Output/June2017/rx_file/No_EFT

%Save Variables
rx_data.est_arrival = estimate;
rx_data.act_arrival = act_arrival;
rx_data.SNR = SNR;
rx_data.x_dist = surface_range;
sname = mat_name(1:end-8);
sname = ['rx_data' '_20' sname(1:2) '_' sname(3:4) '_' sname(5:6) '_' sname(8:9)]
save(sname,'rx_data')

 % update day and hour for the next hour file
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
function [tx_t,tx_lat,tx_lon,tx_heading,x_dist,est_arrival] = posmv_tx_load(now_day,now_hour);
% setup
% 1. Hydrophone Position and Depth
%%% ACO LAT/LON
% ACO_lat = 22.738894;                  % original
% ACO_lon = -158.006009;                % original

% ACO_lat = 22.738772;                  % June 2017
% ACO_lon = -158.006186;                % June2017

ACO_lat= 22.738764;                  % June 2017 1st iteration
ACO_lon= -158.0061781;               % June 2017



%%% ACO Depth
ACO_depth = -4729.92;                         % original depth MSL
% ACO_depth = -4735.29;                          % June 2017 1st iteration

 % 2. Tx File Directory
cd /Volumes/ACO_RAP_2/RAP/June2017Cruise/Tx_Rx_Output/tx_file %% EDIT

% 3. edit CTD file in the ray traing code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create a set of file names
% Load TX data
tx_t = [];
tx_lat = [];
tx_lon = [];
tx_altitude = [];
tx_heading  = [];
fname =[];


% load 2 file
if now_day == 7 & now_hour == 2
    fname = "tx_data_2017_06_07_02";

else
    for p = 1:2
        if p ==1 & now_hour == 0
            if now_day-1 < 10
                str_now_day = "0"+string(now_day-1);
            else
                str_now_day = string(now_day-1);
            end
            fname_d = "tx_data_2017_06_"+str_now_day+"_"+string(23);

        elseif p==1 & now_hour == 10
            if now_day < 10
                str_now_day = "0"+string(now_day);
            else
                str_now_day = string(now_day);
            end
            fname_d = "tx_data_2017_06_"+str_now_day+"_0"+string(9);

        elseif now_hour < 10
            if now_day < 10
                str_now_day = "0"+string(now_day);
            else
                str_now_day = string(now_day);
            end
            fname_d = "tx_data_2017_06_"+str_now_day+"_0"+string(now_hour-2+p);

        else
            if now_day < 10
                str_now_day = "0"+string(now_day);
            else
                str_now_day = string(now_day);
            end
            fname_d = "tx_data_2017_06_"+str_now_day+"_"+string(now_hour-2+p);

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
    x_dist(ii) = distance(tx_lat(ii),tx_lon(ii),ACO_lat,ACO_lon,referenceEllipsoid('WGS84'));
end


%Estimate travel time based on CTD cast
for ii=1:length(x_dist)
    azmth(ii) = azimuth(tx_lat(ii),tx_lon(ii),ACO_lat,ACO_lon);
    [~,~,~,~,~,~,~,~,est_tt(ii),~,~] = ray_trace_w_earth_flattening(x_dist(ii),tx_altitude(ii),tx_lon(ii),tx_lat(ii),azmth(ii),ACO_lat,ACO_lon,ACO_depth,'June','2017','NoEFT');
end

%Estimate arrival time
est_arrival = tx_t+(est_tt./(3600*24));

end
