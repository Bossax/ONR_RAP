% Bottom Bouncing detection
% 1. Cross_correlate the audio signal using the ideal replica
% 2. Detect all peaks in the cluster of audio reception and pick only the first 2 peaks
% 3. Calculate the time lapse between the first peak and the second peak
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

day = 28 ;            %  Edit
start_hour = 16;         % Edit
end_hour = 16;          % Edit
Fs = 32000;
%icListen LAT/LON
icListen_lat=22.739153;                  % Aug 2018
icListen_lon=-158.0061254;                % Aug 2018

%% Load POS MV TX/RX files

cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/tx_file/all
% create a set of file names
fname = [];
now_hour = start_hour-1;
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
for p = 1:length(fname)
    load(fname(p))
    tx_t = horzcat(tx_t,tx_data.t);      %TX time                                  
    tx_lat = horzcat(tx_lat,tx_data.lat);   %TX lat
    tx_lon = horzcat(tx_lon,tx_data.lon);
    tx_altitude = horzcat(tx_altitude,tx_data.altitude);
    tx_heading = horzcat(tx_heading,tx_data.heading);
end


% offset t_tx for icListen time offset
tx_t = tx_t - 6/(3600*24);      

%Boat distance from ACO
for i=1:length(tx_lat)
    x_dist(i)=dist([icListen_lat tx_lat(i)],[icListen_lon tx_lon(i)]);
end



% Load RX file
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

act_arrival = [];
SNR = [];
for i = 1:length(fname)
    load(fname(i)+".mat")
    act_arrival = horzcat(act_arrival ,rx_data.act_arrival);
    SNR = horzcat(SNR ,rx_data.SNR);
    
end

%% Match a TX point with an RX point POSMV
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
tx_lat(rm_ind) = [];
tx_lon(rm_ind) = [];
tx_t(rm_ind) = [];
tx_altitude(rm_ind) = [];
x_dist(rm_ind) = [];
tx_heading(rm_ind) = [];

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
SNR(rm_ind) = [];


%% Load audio file name
% create a set of file names
au_fname = [];
now_hour = start_hour;
now_day = day(1);
while true
    cd("/Volumes/ACO_RAP_2/RAP/Oct2018Cruise/wav_data/icListen/"+string(now_day)+"_Oct")
    while now_hour <= 23
        if now_hour < 10
           fname_d = "SBW1391_201810"+string(now_day)+"_0"+string(now_hour);
        else
           fname_d = "SBW1391_201810"+string(now_day)+"_"+string(now_hour);
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

%%  Bandpass filter for HEM 2000 Hz to 6000 Hz
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

%% Cross Correlation
Arrival_signal = [];               % Output matrix
Envelope = [];
t_sep_bottbounce = [];
rx_t = [];
t_reception = [];
env_peak = [];
y_scatter = [];
envpk_pos = [];
env_pk = [];
add_nx = 0;
add_nx2 = 0;
sum_ind = 0;
cd("/Volumes/ACO_RAP_2/RAP/Oct2018Cruise/wav_data/icListen/"+string(now_day)+"_Oct")
[r,~] = size(au_fname);
for i =1:r
        %Load Audio Data
        wav_name = au_fname(i,:)
        % Timing audio data 
        [y,t_date] = icListen_audio_prep(wav_name);

        % Filtering
        yf = filter(b,1,y); 
        delay = round(mean(grpdelay(b,1,6000,Fs)));      % group delay
        yf = yf(delay+1:end);                       % shift output signal forward in time
        yf(end+1:end+1+delay-1) = zeros(delay,1); 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Envelope Signal
        % Generate Envelope signal and detect the peaks
        [demod,demod_pos_pks,demod_pos_1st_pk,demod_snr_pks,t_sep] = ACO_cross_correlation_ideal_bottom_bounce_icListen(y);
     
    if length(demod_pos_pks) ~= 0
        t_reception = t_date(demod_pos_1st_pk);
        rx_t = horzcat(rx_t,t_reception);
        
        % eliminate peaks without tx points
        rm_ind = [];
        ind_acta = [];
        for k =1:length(t_reception)
            [t_diff,I] = min(abs(t_reception(k) - act_arrival));
             ind_acta(end+1) = I;
            if t_diff*3600*24>1
                rm_ind(end+1) = k;
            end
        end
        % 
        demod_pos_1st_pk(rm_ind) = [];
        t_sep(rm_ind) = [];
         t_reception = t_date(demod_pos_1st_pk);
        
        rm_ind_2pos = [];
        for p = 1:length(rm_ind)
            rm_ind_2pos = horzcat(rm_ind_2pos,2*rm_ind(p)-1:(2*rm_ind(p)));
        end
            demod_pos_pks(rm_ind_2pos) = [];
        
        % t seperation storage
        t_sep_bottbounce = vertcat(t_sep_bottbounce,t_sep);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Audiosignal
      % find indices in time in the audio file that are closest tx times
        ind = zeros(1,length(tx_t));                % indices of transmission timestamps in rx times
        C = zeros(1,length(tx_t));                  % timestamps in t_date
    
        for jj = 1:length(tx_t)
            [c,I] = min(abs(t_date - tx_t(jj)));           % the index whose element is the minimum in t_date
            C(jj) = c;
            ind(jj) = I;
            
            if jj ==1 && I == 1 && c*3600*24 >1 && i == 1
                % the first tx time is of the previous time frame --> zero-padding 
                time_padding = c*3600*24;
                sample_padding = time_padding*Fs;   % the number of samples need padding back in time from the first rx timestamp
                ind(1) = -sample_padding;
            end
            
        end
        % eliminate of which indices are at the first and the last samples as they belong to other time frames
        one_loc = (ind == 1);
        one_ind = find(one_loc == 1);
        end_loc = (ind == length(t_date));
        end_ind = find(end_loc == 1);

        ind([one_ind end_ind]) = [];    % deletion
        
        
         %%%%%% eliminate tx times without peaks
        rm_ind = [];
 
        for k =1:length(ind)
            ind_trep = min(find((t_date(ind(k)) - t_reception)<0));
            t_diff = t_date(ind(k)) - t_reception(ind_trep) ;
            if (t_diff*3600*24 < -20)
                rm_ind(end+1) = k;
            elseif  isempty(t_diff)
                rm_ind(end+1) = k;
            end
        end
        ind(rm_ind) = []; 
      
    %%%%%%%%Cropping%%%%%%%%%%%
    % crop audio signal and envelope
    len = 18 *Fs;                         % 18 sec after tx times (sample)
    ausig_train = zeros(length(ind),len);  % pre-allocate space for time series
    env_train = zeros(length(ind),len); 
    start_row = 1;

    % Deal with a ping in between 2 files
        if add_nx == 1                                         % get the first reception if there is a residual from the previous
            if l_zeropad_sam > 0                               % leading edge is in the last file -> zero-padding
                ausig_train(start_row,:) = [zeros(1,l_zeropad_sam ) yf(1:t_zeropad_sam)'] ; % crop the first reception
                env_train(start_row,:) = [zeros(1,l_zeropad_sam ) demod(1:t_zeropad_sam)'] ; % crop the first reception
            else
                                                                % leading edge is in the same file                                                     
               ausig_train(start_row,:) = yf(abs(l_zeropad_sam):t_zeropad_sam)' ; % crop the first reception 
               env_train(start_row,:) = demod(abs(l_zeropad_sam):t_zeropad_sam)' ; % crop the first reception 
            end
            add_nx = 0 ;                                         % set add_nx to be inactive
            start_row = 2;                                      % increment the row number 
            ausig_train = vertcat(ausig_train,zeros(1,len));        % add an additional row
            env_train = vertcat(env_train,zeros(1,len));    % add an additional row
        end
        %%%%%%%% Rountine 
        
        row = start_row;                         % start row of sig_train
        for j = 1:length(ind)                    % main routine loop over indices tx_time

            if ind(j) < 0
            % Zero-padding
            fprintf('Zero-padding \n')
                leading_edge = ceil(abs(ind(j)))                     % check for the fit of time frame of each index 
                trailing_edge = ind(j)+len-1
                
                ausig_train(row,:) = [zeros(1,leading_edge) yf(1:ceil(trailing_edge))'] ;
                env_train(row,:) = [zeros(1,leading_edge) demod(1:ceil(trailing_edge))'] ;
            
            else

                leading_edge = ind(j);                       % check for the fit of time frame of each index 
                trailing_edge = ind(j)+len-1;

                if trailing_edge < length(yf)                               %window in the same file ALL GOOD!
                            ausig_train(row,:) = yf(leading_edge:trailing_edge)' ;
                            env_train(row,:) = demod(leading_edge:trailing_edge)' ;
                else
                                                                             % This can happen to the last index   
                    handle = min(find(yf(leading_edge:end) >= 1.5e5)) ;  
                    % thredshold 1.5e5

                    if isempty(handle)                                   % reception is in the next file
                        l_zeropad_sam = length(yf) - leading_edge+1;     % zero pad next file
                        t_zeropad_sam = trailing_edge - length(yf);
                        add_nx = 1;                                  % set add_nx to be active
                        ausig_train(row,:) = [];                            % delete the row
                         env_train(row,:) = [];     
                     else
                        l_zeropad_sam = 0;                               % reception is in the same file 
                        t_zeropad_sam = trailing_edge - length(yf);       % zero pad the traling edge  
                        ausig_train(row,:) = [yf(leading_edge:end)' zeros(1,t_zeropad_sam)] ;
                         env_train(row,:) = [demod(leading_edge:end)' zeros(1,t_zeropad_sam)] ;
                    end
                 end 
            end 
            row = row+1;
        end
        
        %%% outputs
        % storage envelope peaks and locations (in the taxis)
        
            for w = 1:length(ind)
                env_pk = vertcat(env_pk,demod(demod_pos_pks((2*w-1):2*w)));
                envpk_pos = vertcat(envpk_pos ,demod_pos_pks((2*w-1):2*w)-ind(w));
            end
       
        % Timestamps for arrival signals
        taxis = 0:1/Fs:1/Fs*(len-1); 
       
        %%%%% Stack time-series
        Arrival_signal = vertcat(Arrival_signal,ausig_train);
        Envelope = vertcat(Envelope,env_train);
        sum_ind = sum_ind+length(ind); 
     end 
%{
        t_plot = (t_date - t_date(1))*24*3600;
        figure(1)
        clf
        plot(t_date,demod)
        hold on
        scatter(t_date(demod_pos_pks),demod(demod_pos_pks))
        plot(t_date,y,'k')
        datetick('x')
        axis tight
        grid minor
        pause

%}
        
end


%% Plot SNR + Audio Signal
stx = 1;                      % Edit
ltx = length(tx_t);               %Edit
range = x_dist(stx:ltx)/1000;
scale_signal = 2/10^5;               % Edit
scale_snr = 1/10^7;               % Edit 

% parameters for display
tx_disp = 3;                   % edit
overlap = 0;                    % edit
move = tx_disp - overlap;
disp_round = ceil(((ltx - stx+1)- move)/move+1);
% Loop to display 
for j = 1:disp_round
    j
    % select signals to display
    if j == 1 
        ssig = stx;
        lsig = ssig+tx_disp-1;
        env_pos= envpk_pos(ssig:2*lsig);
        env_pk_d = env_pk(ssig:2*lsig);
    elseif j == disp_round
        ssig = ltx-tx_disp+1;
        lsig = ltx;
        env_pos= envpk_pos(ssig:2*lsig);
         env_pk_d = env_pk(ssig:2*lsig);
    else
        ssig = (j-1)*move+1;
        lsig = ssig+tx_disp-1;
        env_pos= envpk_pos(2*ssig-1:2*ssig+2*(lsig-ssig));
        env_pk_d = env_pk(2*ssig-1:2*ssig+2*(lsig-ssig));
    end
    % crop audio signal and envelope signals based on arrival time
    %+10 ms and -30 ms after the arrival time
    h_time = 0.01*Fs;
    t_time = 0.03*Fs;
    duration = 0:1/Fs:0.04;
    env_pos_1pk = env_pos(1:2:end);           % fisrt peak of each signal
    % start time and end time of each signal
    start_time_indx = env_pos_1pk - h_time;
    end_time_indx = env_pos_1pk +t_time;
    
    
    % truncate signals
    trun_Arrival_signal = [];
    trun_Envelope = [];
    pk_ind = [];            % index in duration
    env_peak = [];
    % Loop over all receptions
    % truncate signals for plotting
    for k = 1:length(env_pos_1pk)
        trun_Arrival_signal(k,:) = Arrival_signal(ssig+k-1, start_time_indx(k):end_time_indx(k) )-mean(Arrival_signal(ssig+k-1, start_time_indx(k):end_time_indx(k) ));
        trun_Envelope(k,:) = Envelope(ssig+k-1, start_time_indx(k):end_time_indx(k) )-mean( Envelope(ssig+k-1, start_time_indx(k):end_time_indx(k) ));
        pk_ind = vertcat(pk_ind,[env_pos(2*k-1) env_pos(2*k)] -[1 1]*start_time_indx(k) );
        env_peak = vertcat(env_peak,[env_pk_d(2*k-1) env_pk_d(2*k)]-mean( Envelope(ssig+k-1, start_time_indx(k):end_time_indx(k) )));
    end
    
    tx_num  = transpose(ssig:lsig);
    signal = bsxfun(@plus, trun_Arrival_signal*scale_signal, tx_num);  
    SNR = bsxfun(@plus, trun_Envelope*scale_snr, tx_num);  
    
    %%%%%%%%%%%%%%%
    % Plot
    %%%% Surface distance
    f = figure(2);
    f.Units = 'normalized';
    f.Position = [0.1 0.7 0.5 0.7];
    clf
    subplot(1,2,1, 'Position',[0.05 0.1 0.2 0.8])
    scatter(range(ssig:lsig),tx_num,'*');
 
    axis tight
    grid on
    grid minor;
    xlabel('Surface DIstance (km)')
    set(gca,'Xdir','normal')
    set(gca,'YAxisLocation','right');

    % Axis Ticks
    yvalue = [stx:5:ltx-1 ltx];                    %Edit
    y_label = yvalue;                               
    set(gca,'YTick',yvalue);
    set(gca,'Yticklabel',y_label);
    % Label Timestamps
    t_step = 2;                      %Edit
    label_offset = 0.001;             % edit
    date_label = datestr([act_arrival(ssig:t_step:lsig-1) act_arrival(lsig)],'HH:MM');        
%     text([range(ssig:t_step:lsig-1) range(lsig)]+label_offset,[ssig:t_step :lsig-1 lsig],date_label,'FontSize',13)
    
    %%% Arrival Pattern + SNR Plot
    subplot(1,2,2, 'Position',[0.33 0.1 0.65 0.8])
    plot(duration,signal,'k')                   % Edit
    hold on 
    plot(duration,SNR,'Color',[0 0.6 0.8],'LineWidth',1)                  
    grid minor
    axis tight 
    title('Complex Envelopes on Raw signal: icListen')
    set(gca,'fontsize',15)
    xlabel('msec')
    % peak locations
    hold on
    for k  =1:length(env_pos_1pk)
        pk_plot = scatter(duration(pk_ind(k,:)),env_peak(k,:)*scale_snr+tx_num(k),'r');
        
    end
 
    pause
end


%% plot time seperation
figure(1)
scatter(x_dist/1000,t_sep_bottbounce*1000,20,tx_t,'filled')
grid on
xlabel('surface Distance (km)')
ylabel('Time seperation (ms)')
c = colorbar;
c.Label.String = 'Tx Time';
cbdate('HH:MM')
colormap jet
title('Time seperation between the direct arrival and the bottom bounce: Ideal LFM Pulse')
