% Return signal trains matrix containing data points of each reception
% starting from its transmitted time
% Radial Eastward 21-June-2018 11:20 - 12:20
% Radial Westward 21-June-2018 14:40 - 16:20
% Read Hydrophone data
%

clearvars
close all

Arrival_signal = [];               % Output matrix
Envelope = [];
y_scatter = [];             % y displacement of markers in demod

lat = [];
lon = [];
t_tx = [];
sum_ind = 0;
Fs = 24000;

% Bandpass filter 3000-5000 Hz
fc1 = 3000/(Fs/2);
fc2 = 5000/(Fs/2);
Wp = [fc1 fc2];
Ws = Wp + [-500/(Fs/2) +500/(Fs/2)];
Rp = 3; % 3dB
Rs = 30; % -40 dB stopband attenuation  
[n,Wn] = buttord(Wp,Ws,Rp,Rs);
[b,a] = butter(n,Wn,'bandpass');
[h,w] = freqz(b,a,300);

%% Load reception times of an interest interval

cd /Volumes/ACO_RAP_2/RAP/June2018Cruise/Tx_Rx_Output/rx_files_arrival

load('rx_data_ccw.mat')       % Edit
act_arrival = rx_data.act_arrival;
demod_num = 1;


%% Load transmission times of an interest interval

cd /Volumes/ACO_RAP_2/RAP/June2018Cruise/Tx_Rx_Output/tx_file_path/CCW
d_tx=dir('tx_data_*.mat');
lat = [];
lon = [];
t_tx = [];

% Cut out unusable data points 

for p = 1:length(d_tx)
    rm_ind = [];
    fname = d_tx(p).name;
    load(fname)
    t_dum = tx_data.t;
    lat_dum = tx_data.lat;
    lon_dum = tx_data.lon;
    for ii = 1:length(t_dum)
        gap = (act_arrival-t_dum(ii))*3600*24;
        if isempty(find(0 < gap & gap < 28))
           rm_ind(end+1) = ii; 
        end
    end
    t_dum(rm_ind) = [];
    lat_dum(rm_ind) = [];
    lon_dum(rm_ind) = [];
    t_tx = horzcat(t_tx,t_dum);                      % Transmission timestamps
    % store lat and lon data
    lat = horzcat(lat,lat_dum);
    lon = horzcat(lon,lon_dum);
end
lat = transpose(lat);
lon = transpose(lon);

%% find matching timestamps in the hydrophone file
add_nx = 0;
add_nx2 = 0;
y_len = [];
% Deal with one audio file at a time
cd /Volumes/ACO_RAP_2/RAP/June2018Cruise/wav_data/Hydrophone_Raw/Hydrophone_Data/CCW% Edit
hyd = dir('18*');
addpath('/Volumes/ACO_RAP_2/RAP/June2018Cruise/Script/Signal_Arrival')

for jj = 1:2:length(hyd)-1
    (jj+1)/2
    mat_name = hyd(jj).name;
    wav_name = hyd(jj+1).name; 
    
    %Load Audio Data
    [y,t_date] = hyd_audio_prep(mat_name,wav_name);
    y = y';
    
    % Filtering
    y = filter(b,a,y);
    t_rx = t_date;

    % find indices in rx times that are close to tx times
    ind = zeros(1,length(t_tx));                % indices of transmission timestamps in rx times
    C = zeros(1,length(t_tx));                  % values of timestamps in rx times
    
    for i = 1:length(t_tx)
        [c,I] = min(abs(t_rx - t_tx(i)));           % the index whose element is the minimum in t_rx
        C(i) = c;
        ind(i) = I;   
        if i ==1 && I == 1 && c*3600*24 >1 && jj == 1
            % the first tx time is of the previous time frame --> zero-padding 
            time_padding = c*3600*24;
            sample_padding = time_padding*Fs;   % the number of samples need padding back in time from the first rx timestamp
            ind(1) = -sample_padding;
        end
    end

    % eliminate  of which indices are at the first and the last samples as they belong to other time frames
    one_loc = (ind == 1);
    one_ind = find(one_loc == 1);
    end_loc = (ind == length(t_rx));
    end_ind = find(end_loc == 1);

    ind([one_ind end_ind]) = [];    % deletion


    % Generate Envelope signal
    load LFM_pulse_sample.mat                         % Replica pulse for cross correlation
    [xc,lag]=xcorr(y,replica_pulse);
    xc(1:find(lag==0)-1)=[];    
    demod=abs(hilbert(xc)); 
     
    % find indices in rx times that are close to demod peak time (arrival time) 
        ind_act = zeros(1,length(act_arrival));         % indices of demod peak timestamps in the audio file           
        C = zeros(1,length(act_arrival));                     
        for i = 1:length(act_arrival)
            [c,I] = min(abs(t_rx - act_arrival(i)));   
            C(i) = c;
            ind_act(i) = I;  
        end

        % eliminate  of which indices are at the first and the last samples as they belong to other time frames
        one_loc = (ind_act == 1);
        one_ind = find(one_loc == 1);
        end_loc = (ind_act == length(t_rx));
        end_ind = find(end_loc == 1);

        ind_act([one_ind end_ind]) = [];   % deletion
        
        % markers in demod plot
        y_scatter(end+1:end+length(ind_act)) = demod(ind_act);
    
    %%%%%%%%Cropping%%%%%%%%%%%
    
    % crop signal and demod
    len = 18 *Fs;                         % 18 sec after tx times (sample)
    sig_train = zeros(length(ind),len);  % pre-allocate space for time series
    demod_train = zeros(length(ind),len); 
    start_row = 1;

    % Deal with a ping in between 2 files
        if add_nx == 1                                         % get the first reception if there is a residual from the previous
            if l_zeropad_sam > 0                               % leading edge is in the last file -> zero-padding
                sig_train(start_row,:) = [zeros(1,l_zeropad_sam ) y(1:t_zeropad_sam)] ; % crop the first reception
                demod_train(start_row,:) = [zeros(1,l_zeropad_sam ) demod(1:t_zeropad_sam)] ; % crop the first reception
            else
                                                                % leading edge is in the same file                                                     
               sig_train(start_row,:) = y(abs(l_zeropad_sam):t_zeropad_sam) ; % crop the first reception 
               demod_train(start_row,:) = demod(abs(l_zeropad_sam):t_zeropad_sam) ; % crop the first reception 
            end
            add_nx = 0 ;                                         % set add_nx to be inactive
            start_row = 2;                                      % increment the row number 
            sig_train = vertcat(sig_train,zeros(1,len));        % add an additional row
            demod_train = vertcat(demod_train,zeros(1,len));    % add an additional row
        end
        %%%%%%%% Rountine 
        
        row = start_row;                                % start row of sig_train
        for j = 1:length(ind)                               % main routine loop over indices x

            if ind(j) < 0
            % Zero-padding
            fprintf('Zero-padding \n')
                leading_edge = ceil(abs(ind(j)))                     % check for the fit of time frame of each index 
                trailing_edge = ind(j)+len-1
                
                sig_train(row,:) = [zeros(1,leading_edge) y(1:ceil(trailing_edge))] ;
                demod_train(row,:) = [zeros(1,leading_edge) demod(1:ceil(trailing_edge))] ;
            
            else

                leading_edge = ind(j);                       % check for the fit of time frame of each index 
                trailing_edge = ind(j)+len-1;

                if trailing_edge < length(y)                               %window in the same file ALL GOOD!
                            sig_train(row,:) = y(leading_edge:trailing_edge) ;
                            demod_train(row,:) = demod(leading_edge:trailing_edge) ;
                else
                                                                             % This can happen to the last index   
                    handle = min(find(y(leading_edge:end) >= 1.5e5)) ;  
                    % thresdhole 1.5e5

                    if isempty(handle)                                   % reception is in the next file
                        l_zeropad_sam = length(y) - leading_edge+1     % zero pad next file
                        t_zeropad_sam = trailing_edge - length(y)
                        add_nx = 1                                  % set add_nx to be active
                        sig_train(row,:) = [];                            % delete the row
                         demod_train(row,:) = [];     
                     else
                        l_zeropad_sam = 0;                               % reception is in the same file 
                        t_zeropad_sam = trailing_edge - length(y);       % zero pad the traling edge  
                        sig_train(row,:) = [y(leading_edge:end) zeros(1,t_zeropad_sam)] ;
                         demod_train(row,:) = [demod(leading_edge:end) zeros(1,t_zeropad_sam)] ;
                    end
                 end 
            end 
            row = row+1;
        end
        % Timestamps for arrival signals
        taxis = 0:1/Fs:1/Fs*(len-1); 
        
         %%%%%%%%Plotting%%%%%%%%
         %{     
        figure(1)
        clf
        subplot(2,1,1)
        plot(y)
        grid on
        axis tight
        hold on
        scatter(ind,zeros(1,length(ind)))
         
        subplot(2,1,2)
        plot(demod)
        grid on
        axis tight
        hold on
        scatter(ind_act,demod(ind_act))
         
       pause
       
            figure(3)
            clf
            plot(y)
            grid on
            axis tight
            hold on
            scatter(ind,zeros(1,length(ind)))
            pause
     
         
          %%% Debug demodulation
            figure(4)
           
            s = size(demod_train);
            
            for k =1:s(1)
                plot(taxis,demod_train(k,:))
                hold on
                %axis tight
                grid on
                title(num2str(demod_num))
                xlim([6 14])
                pause
                demod_num = demod_num+1;
            end
         
             %}
         if jj == 0
          
          s = size(demod_train);
             figure(5)
           f_sam = 0.01*Fs;  % 10msec
           s_sam = 0.020*Fs; % 20 msec
            
            for k =1:s(1)
                clf
                plot(t_rx,demod)
                hold on
                scatter(t_rx(ind_act(k)),demod(ind_act(k)),'o')
                %axis tight
                grid on
                title(num2str(demod_num))
                xlim([t_rx(ind_act(k)-f_sam) t_rx(ind_act(k)+s_sam)])
                pause
                demod_num = demod_num+1;
            end
         end

        %%%%%% Plot to adjust time windowing length

%{
        figure(2)
        clf
        s = size(sig_train);

        for k =1:s(1)
            subplot(2,1,1)

            plot(taxis,sig_train(k,:))
            hold on
            axis tight
            grid on
            title(num2str(k))

            subplot(2,1,2)

            plot(taxis,demod_train(k,:))
            hold on
            axis tight
            grid on
            title(num2str(k))
            pause

        end
%} 
    %%%%% Stack time-series
    Arrival_signal = vertcat(Arrival_signal,sig_train);
    Envelope = vertcat(Envelope,demod_train);
    sum_ind = sum_ind+length(ind); 
end


%% Mark rx time corresond to cropped signals
time_offset_tx_act = (act_arrival - t_tx)*3600*24;  % sec from tx timestamp

%% Calculate Surface Distance
ACO_lat = 22.738772*ones(length(lat),1) ;       % North Hydrophone position updated by Vincent March 2018
ACO_lon = -158.006186*ones(length(lon),1);      % West
geodesic_point = 10000;
LAT = [lat ACO_lat];
LON = [lon ACO_lon];
range = [];
lat_cur = [];
lon_cur = [];
for p = 1:length(lat)
[range(p),lat_cur(p),lon_cur(p)]=dist(LAT(p,:),LON(p,:)); % range in meter
end
range = range/1000 ; % range in km

%% Prepare parameters for plotting figures
start_tx_num = 210;
tx_num_max = 230;               %Edit
              %Edit
tx_numbering = start_tx_num:tx_num_max;
% Time Framing
start_time = 7.31;                 % Edit
end_time = 7.4;                 % Edit
start_time_indx = (start_time)*Fs +1;
end_time_indx = (end_time)*Fs +1;

scale_signal = 7/10^7;               % Edit
signal = bsxfun(@plus, Arrival_signal(start_tx_num:tx_num_max,start_time_indx:end_time_indx )*scale_signal, tx_numbering(:));  

scale_snr = 2/10^15;               % Edit
SNR = bsxfun(@plus, Envelope(start_tx_num:tx_num_max,start_time_indx:end_time_indx )*scale_snr, tx_numbering(:));  
%% Plot SNR on Raw signal

% Surface Distance Plot
f = figure(3);
f.Units = 'normalized';
f.Position = [0.1 0.7 0.8 0.7];

subplot(1,2,1, 'Position',[0.05 0.1 0.2 0.8])
scatter(range(start_tx_num:tx_num_max),tx_numbering(:),'*');

axis tight
grid on
grid minor;
xlabel('Surface DIstance (km)')
set(gca,'Xdir','normal')
set(gca,'YAxisLocation','right');

% Axis Ticks
yvalue = [start_tx_num:5:tx_num_max-1 tx_num_max];                    %Edit
y_label = yvalue;                  
xvalue = [start_time:0.01:end_time];                 %Edit                  
x_label = xvalue  ;               
set(gca,'YTick',yvalue);
set(gca,'Yticklabel',y_label);

title('CCW Hyd')
                       

% Label Timestamps
t_step = 5; %Edit
label_offset = -0.0085;
date_label = datestr([t_tx(start_tx_num:t_step:tx_num_max-1) t_tx(tx_num_max)],'HH:MM');        
text([range(start_tx_num:t_step:tx_num_max-1) range(tx_num_max)]+label_offset,[start_tx_num:t_step :tx_num_max-1 tx_num_max],date_label,'FontSize',13)

% Arrival Pattern + SNR Plot
subplot(1,2,2, 'Position',[0.33 0.1 0.65 0.8])
plot(taxis(start_time_indx :end_time_indx ),signal,'k')                   % Edit
hold on 
plot(taxis(start_time_indx :end_time_indx ),SNR,'Color',[0 0.6 0.8])                  

axis tight 
title('Envelopes on Raw signal ')


%Rx times plot
hold on

rx_plot = scatter(time_offset_tx_act(start_tx_num:tx_num_max),tx_numbering+y_scatter(start_tx_num:tx_num_max)*scale_snr);
rx_plot.Marker = 'o';
rx_plot.MarkerEdgeColor = 'red';
rx_plot.LineWidth = 0.8;


%Rx times plot on raw signal
hold on
rx_plot = scatter(time_offset_tx_act(start_tx_num:tx_num_max-1),tx_numbering(1:end-1),100);
rx_plot.Marker = 'x';
rx_plot.MarkerEdgeColor = 'red';
rx_plot.LineWidth = 1.5;

grid on
grid minor
set(gca,'YTick',yvalue);
set(gca,'Yticklabel',y_label);
set(gca,'XTick',xvalue);
set(gca,'Xticklabel',x_label);
ylim([start_tx_num tx_num_max])                                  
xlabel('Time from transmission')
ylabel('Transmission order Number')

%% Plot arrivals patterns
% Surface Distance Plot
f = figure(1);
f.Units = 'normalized';
f.Position = [0.1 0.7 0.8 0.7];

subplot(1,2,1, 'Position',[0.05 0.1 0.2 0.8])
scatter(range(start_tx_num:tx_num_max),tx_numbering(:),'*');

axis tight
grid on
grid minor;
xlabel('Surface DIstance (km)')
set(gca,'Xdir','reverse')
set(gca,'YAxisLocation','right');

% Axis Ticks
yvalue = [start_tx_num:10:tx_num_max-1 tx_num_max];                    %Edit
y_label = yvalue;                  
xvalue = [start_time:0.5:end_time];                 %Edit                  
x_label = xvalue  ;               
set(gca,'YTick',yvalue);
set(gca,'Yticklabel',y_label);

title('Departure')
                       

% Label Timestamps
t_step = 10; %Edit
label_offset = 1.4;
date_label = datestr([t_tx(start_tx_num:t_step:tx_num_max-1) t_tx(tx_num_max)],'HH:MM');        
text([range(start_tx_num:t_step:tx_num_max-1) range(tx_num_max)]+label_offset,[start_tx_num:t_step :tx_num_max-1 tx_num_max],date_label,'FontSize',13)

% Arrival Pattern Plot
subplot(1,2,2, 'Position',[0.33 0.1 0.65 0.8])
plot(taxis(start_time_indx :end_time_indx ),signal,'k')                   % Edit
axis tight 

%{
%Rx times plot
hold on
rx_plot = scatter(time_offset_tx_act(start_tx_num:tx_num_max),tx_numbering);
rx_plot.Marker = 'o';
rx_plot.MarkerEdgeColor = 'red';
rx_plot.LineWidth = 0.8;
%}

title('Raw Signal Arrival Pattern')
grid on

set(gca,'YTick',yvalue);
set(gca,'Yticklabel',y_label);
set(gca,'XTick',xvalue);
set(gca,'Xticklabel',x_label);
ylim([start_tx_num tx_num_max])                                  
xlabel('Time from transmission')
ylabel('Transmission order Number')


%% Plot SNR
% Surface Distance Plot
f = figure(2);
f.Units = 'normalized';
f.Position = [0.1 0.7 0.8 0.7];
title('Eastward Radial')

subplot(1,2,1, 'Position',[0.05 0.1 0.2 0.8])
scatter(range(start_tx_num:tx_num_max),tx_numbering(:),'*');

axis tight
grid on
grid minor;
xlabel('Surface DIstance (km)')
set(gca,'Xdir','reverse')
set(gca,'YAxisLocation','right');

% Axis Ticks
yvalue = [start_tx_num:5:tx_num_max-1 tx_num_max];                    %Edit
y_label = yvalue;                  
xvalue = [start_time:0.5:end_time];                 %Edit                  
x_label = xvalue  ;               
set(gca,'YTick',yvalue);
set(gca,'Yticklabel',y_label);

title('Eastward Radial')
                       

% Label Timestamps
date_label = datestr([t_tx(start_tx_num:10:tx_num_max-1) t_tx(tx_num_max)],'HH:MM');        %Edit
text([range(start_tx_num:10:tx_num_max-1) range(tx_num_max)]+3,[start_tx_num:10:tx_num_max-1 tx_num_max],date_label,'FontSize',13)

% SNR Plot
subplot(1,2,2, 'Position',[0.33 0.1 0.65 0.8])
plot(taxis(start_time_indx :end_time_indx ),SNR,'Color',[0 0.6 0.8])                  
axis tight 
title('Envelope Signal from Demodulation')

%Rx times plot
hold on
y_scatter_db = max(Envelope,[],2)';
rx_plot = scatter(time_offset_tx_act(start_tx_num:tx_num_max),tx_numbering+y_scatter(start_tx_num:tx_num_max)*scale_snr);
rx_plot.Marker = 'o';
rx_plot.MarkerEdgeColor = 'red';
rx_plot.LineWidth = 0.8;


grid on

set(gca,'YTick',yvalue);
set(gca,'Yticklabel',y_label);
set(gca,'XTick',xvalue);
set(gca,'Xticklabel',x_label);
ylim([start_tx_num tx_num_max])                                  
xlabel('Time from transmission')
ylabel('Transmission order Number')

