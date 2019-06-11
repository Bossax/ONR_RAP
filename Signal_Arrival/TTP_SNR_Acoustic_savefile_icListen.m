% Travel Time Perturbation Plot combined with signal arrival patterns
% 1. Cross_correlate the audio signal using the ideal replica
% 2. Detect all peaks in the cluster of audio reception and pick only the first 2 peaks
% 3. Calculate the time lapse between the first peak and the second peak
% 4. Calculate Travel time perturbation
% 5. Plot ttp and signal arrival pattern on the same figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

day = 27:27 ;            %  Edit
start_hour = 3;         % Edit
end_hour = 9;          % Edit
Fs = 32000;
    
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
rx_t = [];
est_reception = [];
act_reception = [];               % arrival time
env_pk = [];                      % peak magnitude of envelope signals
envpk_pos = [];
tx_time =[];
arrival_mark = [];                  % markers of arrival points of time in the audio file
tx_range = [];
heading = [];
sum_ind = 0;

now_hour = start_hour;
now_day = day(1);
counter = 1;
[len_fname,~] = size(au_fname);

add_nx = 0;
add_nx2 = 0;
tx_from_pre  = [];
tx_sd2nex  = [];
t_zeropad_sam = [];
l_zeropad_sam = [];
while true
    % day and hour
    % load tx file and calculate estimated arrival times
    [tx_t,tx_lat,tx_lon,tx_heading,x_dist,est_arrival] = posmv_tx_load(now_day,now_hour);
    
    % minutes
    m = 1;
    wav_name = au_fname(counter,:);
    current_hour = str2num(wav_name(18:19))
       while true  
             
             cd("/Volumes/ACO_RAP_2/RAP/Oct2018Cruise/wav_data/icListen/"+string(now_day)+"_Oct")
              %Load Audio Data
              wav_name = au_fname(counter,:);
              % Timing audio data 
             % Timing audio data 
             [y,t_date] = icListen_audio_prep(wav_name);
             time_offset = time_correction(t_date(1));       % load time offset
             t_date = t_date - time_offset;                  % icListen time offset from UTC
             sprintf("file end time %s",datestr(t_date(end)))
             % Filtering
            r = 0.1;
            window = tukeywin(length(y),r)';
            y_ic4_f =  y.*window';
            invwin_ic = [0 0 0 0 0 0 0 1./window(8:end-7) 0 0 0 0 0 0 0];
             
             yf = filter(b,1,y_ic4_f); 
             delay = round(mean(grpdelay(b,1,6000,Fs)));      % group delay
             yf = yf(delay+1:end);                       % shift output signal forward in time
             yf(end+1:end+1+delay-1) = zeros(delay,1); 
             yf = yf.*invwin_ic';
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             % Envelope Signal
             % Generate Envelope signal and detect the peaks
             [demod,demod_pos_pks,demod_pos_1st_pk,demod_snr_pks,t_sep] = ACO_cross_correlation_ideal_bottom_bounce_icListen(yf);
             
             
             t_reception = t_date(demod_pos_1st_pk);
             
             %%% Eliminate peaks without tx points
                rm_ind = [];                
                for k =1:length(t_reception)

                    [~,I] = max(find((t_reception(k) -tx_t)>0));
                     t_diff =t_reception(k)- tx_t(I);
                  
                    if t_diff*3600*24>20
                        rm_ind(end+1) = k;
                    end
                end

                demod_pos_1st_pk(rm_ind) = [];
                t_sep(rm_ind) = [];
                t_reception = t_date(demod_pos_1st_pk);
               
                rm_ind_2pos = [];
                for p = 1:length(rm_ind)
                    rm_ind_2pos = horzcat(rm_ind_2pos,2*rm_ind(p)-1:(2*rm_ind(p)));
                end
                demod_pos_pks(rm_ind_2pos) = [];
                
            %%%%%%%%Cropping%%%%%%%%%%%
            % crop audio signal and envelope
            len = 20 *Fs;                         % 18 sec after tx times (sample)
%             sprintf("The added tx time %s",datestr(tx_sd2nex))
            [trun_ausig_train,trun_env_train,pk_ind,env_peak,tx_td,range,heading_d,est_reception_d,t_reception,add_nx,add_nx2,l_zeropad_sam,t_zeropad_sam,tx_from_pre,tx_sd2nex] = truncate_timeseries(yf,t_date,demod,demod_pos_1st_pk,demod_pos_pks,tx_t,x_dist,tx_heading,est_arrival,t_reception,len,add_nx,add_nx2,l_zeropad_sam,t_zeropad_sam,tx_from_pre,tx_sd2nex);
                

            %%% outputs
            % storage envelope peaks and locations (in the ttp axis)
            env_pk = vertcat(env_pk,env_peak);
            envpk_pos = vertcat(envpk_pos,pk_ind);


            %%%%% Stack time-series
            Arrival_signal = vertcat(Arrival_signal,trun_ausig_train);
            Envelope = vertcat(Envelope,trun_env_train);
            act_reception = vertcat(act_reception,t_reception');             
            est_reception = vertcat(est_reception,est_reception_d);
            tx_time = vertcat(tx_time,tx_td);
            tx_range = vertcat(tx_range,range);
            heading = vertcat(heading,heading_d);

            if m == 12
            %%% Save data
             ttp_snr_sig_data.sig = Arrival_signal;
             ttp_snr_sig_data.snr = Envelope;
             ttp_snr_sig_data.est_arrival = est_reception;
             ttp_snr_sig_data.act_arrival = act_reception;
             ttp_snr_sig_data.tx_t = tx_time;
             ttp_snr_sig_data.env_pk = env_pk;
             ttp_snr_sig_data.envpk_pos = envpk_pos;
             ttp_snr_sig_data.range = tx_range;
             ttp_snr_sig_data.heading = heading;
             sname = "ttp_snr_sig_data"+ "_2018_10_"+string(now_day)+ "_"+ string(current_hour)+"_icListen"
             cd '/Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/ttp_snr_plot/icListen'
             save(sname,'ttp_snr_sig_data' );

              Arrival_signal = [];               % Output matrix
              Envelope = [];
              est_reception = [];
              act_reception = [];                   % arrival time
              env_pk = [];                      % peak magnitude of envelope signals
              envpk_pos = [];
              tx_time = [];
              tx_range = [];
              heading = [];
            end
            counter = counter+1;
            m = m+1;
            % check if the file jumps to the next hour
        
            if counter > len_fname
                break;
            end
            nx_fname = au_fname(counter,:);
            chk_hour =  str2double(nx_fname(18:19));

            if chk_hour ~= current_hour
                break;
            end
            
        end
        
   
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tx_t,tx_lat,tx_lon,tx_heading,x_dist,est_arrival] = posmv_tx_load(day,hour)
    %% Load POS MV TX/RX files
    cd('/Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/tx_file/all/3')
    % create a set of file names
    %ACO LAT/LON
    icListen_lat = 22.739153;                   % June 2018 north of HEM  42.22 m nort of HEM(June 2017 Position)   
    icListen_lon = -158.0061254;                % June 2018 east of HEM  5.94 m east of HEM(June 2017 Position)   
    icListen_depth = -4729.92+2.32 +1.75;          % original depth ellipsoid height+ 1.75 m higher
    fname = [];
    now_hour = hour;
    now_day = day;
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

   % Geodic Length
   for ii=1:length(tx_lat)
       x_dist(ii) = distance(tx_lat(ii),tx_lon(ii),icListen_lat,icListen_lon,referenceEllipsoid('WGS84'));
   end
   
   % Estimated travel time based on CTD cast
   for ii=1:length(x_dist)
       azmth(ii) = azimuth(tx_lat(ii),tx_lon(ii),icListen_lat,icListen_lon);
       [~,~,~,~,~,~,~,~,est_tt(ii),~,~] = ray_trace_w_earth_flattening(x_dist(ii),tx_altitude(ii),tx_lat(ii),azmth(ii),icListen_lat,icListen_lon,icListen_depth);
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

function time_offset = time_correction(time)
 % return time offset for timestamp correction
 date_mark1=  "20181030 05:55";
 
 date_mark1 = datenum(date_mark1,'yyyymmdd HH:MM');
 
 if time <= date_mark1
     time_offset =  -7/(3600*24);
     
 else
     time_offset = +1/(3600*24);
     
 end
 
 end

%%
function [trun_ausig_train,trun_env_train,pk_ind,env_peak,tx_td,range,heading,est_reception_d,t_reception,add_nx,add_nx2,l_zeropad_sam,t_zeropad_sam,tx_from_pre,tx_sd2nex] = truncate_timeseries(y,t_date,demod,demod_pos_1st_pk,demod_pos_pks,tx_t,x_dist,tx_heading,est_arrival,t_reception,len,add_nx,add_nx2,l_zeropad_sam,t_zeropad_sam,tx_from_pre,tx_sd2nex);
    Fs = 32000; 
    %%%%% find indices in time in the audio file that are closest tx times
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
    % liminate of which indices are at the first and the last samples as they belong to other time frames
    one_loc = (ind == 1);
    one_ind = find(one_loc == 1);
    end_loc = (ind == length(t_date));
    end_ind = find(end_loc == 1);

    ind([one_ind end_ind]) = [];    % deletion
    %%%%%%%%%%%%%%%%%%
        if length(ind) ~=0
    
    
        ausig_train = zeros(length(ind),len);  % pre-allocate space for time series
        env_train = zeros(length(ind),len); 
        start_row = 1;
        % check if the first sample is zero-padded
        if add_nx == 1
            n_offset = l_zeropad_sam;
        else
            n_offset = 0;
        end
    
    % Deal with a ping in between 2 files
        if add_nx == 1                                         % get the first reception if there is a residual from the previous
            if l_zeropad_sam > 0                               % leading edge is in the last file -> zero-padding
                ausig_train(start_row,:) = [zeros(1,l_zeropad_sam ) y(1:t_zeropad_sam)'] ; % crop the first reception
                env_train(start_row,:) = [zeros(1,l_zeropad_sam ) demod(1:t_zeropad_sam)'] ; % crop the first reception
            else
                                                                % leading edge is in the same file                                                     
               ausig_train(start_row,:) = y(abs(l_zeropad_sam):t_zeropad_sam)' ; % crop the first reception 
               env_train(start_row,:) = demod(abs(l_zeropad_sam):t_zeropad_sam)' ; % crop the first reception 
            end
            add_nx = 0 ;                                         % set add_nx to be inactive
            tx_from_pre = tx_sd2nex;    
            tx_sd2nex = [];     % clear the send next
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
                
                ausig_train(row,:) = [zeros(1,leading_edge) y(1:ceil(trailing_edge))'] ;
                env_train(row,:) = [zeros(1,leading_edge) demod(1:ceil(trailing_edge))'] ;
            
            else

                leading_edge = ind(j);                       % check for the fit of time frame of each index 
                trailing_edge = ind(j)+len-1;

                if trailing_edge < length(y)                               %window in the same file ALL GOOD!
                            ausig_train(row,:) = y(leading_edge:trailing_edge)' ;
                            env_train(row,:) = demod(leading_edge:trailing_edge)' ;
                else
                    
                    %  handle = min(find(y(leading_edge:end) >= 0.25e4));   % thredshold 0.5e5
                    
                    % last reception
                    if ~isempty(demod_pos_1st_pk)
                        if  demod_pos_1st_pk(end) > length(y)              % reception is in the next file
                            l_zeropad_sam = length(y) - leading_edge+1;     % zero pad next file
                            t_zeropad_sam = trailing_edge - length(y);
                            add_nx = 1                                       % set add_nx to be active
                            tx_sd2nex= t_date(ind(j));                      % send the tx t to the next file
                            ausig_train(row,:) = [];                            % delete the row
                            env_train(row,:) = [];
                        else
                            l_zeropad_sam = 0;                               % reception is in the same file
                            t_zeropad_sam = trailing_edge - length(y);       % zero pad the traling edge
                            ausig_train(row,:) = [y(leading_edge:end)' zeros(1,t_zeropad_sam)] ;
                            env_train(row,:) = [demod(leading_edge:end)' zeros(1,t_zeropad_sam)] ;
                        end
                        
                    end
                end
            end 
            row = row+1;
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Truncation for ttp plot
    %%% shift the axis from zero at tx time to zero at zero ttp %%%%%%%%%%%%%%%%%%%        
    
    % estimated/actual arrivals and ttp and tx times
     if n_offset ~= 0
        tx_td = [tx_from_pre t_date(ind)];
     else
         tx_td = t_date(ind);
     end
         
    est_reception_d = [];
    
    for j = 1:length(tx_td)
        
        [I,~] = min(find((est_arrival - tx_td(j))>0));
        t_diff = est_arrival(I) - tx_td(j);
        if t_diff*3600*24 < 20
        est_reception_d(end+1) = est_arrival(I);
        
        end
    end
                
    
    
    % crop audio signal and envelope signals based on ttp
    % -20 ms and +20 ms about zero travel time perturbation
    h_time = -0.03;
    t_time = 0.03;
    

    %%%%%%%%%%%%%%%%%% prepare peak positions %%%%%%%%%%%%%%%%%%%%%
    %%% in case of no reception is detected 
    add_ind = [];
    for j = 1:length(est_reception_d)
        % check for missing times
        t_diff = min(abs((est_reception_d(j) - t_date(demod_pos_1st_pk))*3600*24));
        if t_diff >0.05 
            % add time index in the demod
            [~,I] = min(abs(est_reception_d(j)-t_date));
            add_ind(end+1) = I;
       elseif isempty(t_diff)
            % add time index in the demod
            [~,I] = min(abs(est_reception_d(j)-t_date));
            add_ind(end+1) = I;
       end
        
    end
    % belongs to the next file
    nx_f = find(add_ind == length(t_date));
    if length(nx_f)~=0
        add_ind(nx_f) = [];
        tx_td(end) = [];
        est_reception_d(end) = [];
        ind(end) = [];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%% add missing reception times %%%%%%%%%%%%%%%%%%%%
    demod_pos_1st_pk = sort(vertcat(demod_pos_1st_pk,add_ind'));
    for k = 1:length(add_ind)
        demod_pos_pks = sort(vertcat(demod_pos_pks,[1;1]*add_ind(k)));
    end
   
   t_reception = t_date(demod_pos_1st_pk);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   %%%%%%%%%%%%%% remove t_reception without est times %%%%%%%%%%%%%%%%
   %%%% in case of error in peak detection
    rm_ind = [];         
    for k =1:length(t_reception)

        [~,I] = max(find((t_reception(k) -tx_td)>0));
        t_diff =t_reception(k)- tx_td(I);
                  
        if t_diff*3600*24>20
            rm_ind(end+1) = k;
        elseif isempty(I)
             rm_ind(end+1) = k;
        end
    end
       
    demod_pos_1st_pk(rm_ind) = [];
    t_reception = t_date(demod_pos_1st_pk);

    rm_ind_2pos = [];
    for w = 1:length(rm_ind)
        rm_ind_2pos = horzcat(rm_ind_2pos,2*rm_ind(w)-1:(2*rm_ind(w)));
    end
    demod_pos_pks(rm_ind_2pos) = [];
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
   
    %%%%%% eliminate act arrival without est arrival
    %%%% in case of error during peak detection
    rm_ind = [];         
    for k =1:length(t_reception)

        [t_diff,I] = min(abs(t_reception(k) -est_reception_d));
                  
        if t_diff*3600*24>0.05
            rm_ind(end+1) = k;
        elseif isempty(I)
             rm_ind(end+1) = k;
        end
    end
       
    demod_pos_1st_pk(rm_ind) = [];
    t_reception = t_date(demod_pos_1st_pk);

    rm_ind_2pos = [];
    for w = 1:length(rm_ind)
        rm_ind_2pos = horzcat(rm_ind_2pos,2*rm_ind(w)-1:(2*rm_ind(w)));
    end
    demod_pos_pks(rm_ind_2pos) = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   %%%%%%%%% shift peak locations to the tx axis (zero point at tx time)
   envpk_pos_hd = [];           % offset samples from the tx time to the peaks
   env_pk_hd  = [];
   u = 1;
    for w = 1:length(t_reception)
        if (w==1)& (n_offset ~= 0)
            env_pk_hd = vertcat(env_pk_hd,demod(demod_pos_pks(1:2))');
            envpk_pos_hd = vertcat(envpk_pos_hd,demod_pos_pks(1:2)'+n_offset);
        else
            env_pk_hd = vertcat(env_pk_hd,demod(demod_pos_pks((2*w-1):2*w))');
            envpk_pos_hd = vertcat(envpk_pos_hd,demod_pos_pks((2*w-1):2*w)'-ind(u));
            u = u+1;
        end
    end 
    env_pos_1pk = transpose(envpk_pos_hd(1:end,1));           % fisrt peak of each signal (in tx time domain)
  
   
    % find zero ttp indicies
    % start time and end time of each signal
    shift_sam = zeros(length(tx_td),1);        % in tx axis
    start_time_indx = zeros(length(tx_td),1);  % in tx axis
    end_time_indx = zeros(length(tx_td),1);    % in tx axis
    for n= 1:length(tx_td)
        shift_sam(n) = round((est_reception_d(n) - tx_td(n))*3600*24*Fs);    % shift 0 from tx time to zerp ttp
        start_time_indx(n) = shift_sam(n) + h_time*Fs;
        end_time_indx(n) = shift_sam(n) +t_time*Fs;
     end


     % truncate signals
     trun_ausig_train = [];
     trun_env_train = [];
     pk_ind = [];            % index in duration
     env_peak = [];


     % Loop over all receptions
     % truncate signals for plotting
     % shift peak locs to zero at est arrival (ttp)
     for k = 1:length(tx_td)
        trun_ausig_train(k,:) = ausig_train(k, start_time_indx(k):end_time_indx(k) )-mean(ausig_train(k, start_time_indx(k):end_time_indx(k) ));
        trun_env_train(k,:) = env_train(k, start_time_indx(k):end_time_indx(k) )-mean(env_train(k, start_time_indx(k):end_time_indx(k) ));
        pk_ind = vertcat(pk_ind,[envpk_pos_hd(k,1) envpk_pos_hd(k,2)] -[1 1]*start_time_indx(k) );
        env_peak = vertcat(env_peak,([env_pk_hd(k,1) env_pk_hd(k,2)]-mean(env_train(k, start_time_indx(k):end_time_indx(k))))/max(trun_env_train(k,:)));
        trun_ausig_train(k,:)=trun_ausig_train(k,:)/max(trun_ausig_train(k,:));
        trun_env_train(k,:) = trun_env_train(k,:)/max(trun_env_train(k,:));

     end
     tx_td = tx_td';
     range = [];
     heading = [];
     est_reception_d = est_reception_d';
     for u = 1:length(tx_td);
         [~,I]= min(abs(tx_t - tx_td(u)));
         range = vertcat(range,x_dist(I));
         heading = vertcat(heading,tx_heading(I));
     end
     
     pk_ind
    else
     fprintf('No ind')   
     tx_from_pre = [];   
     tx_sd2nex = [];
     tx_td = [];
     t_reception = [];
     trun_ausig_train = [];
     est_reception_d = [];
     trun_env_train = [];
     pk_ind = [];            % index in duration
     env_peak = [];
     range = [];
     heading = [];
    end

end