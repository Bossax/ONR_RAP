% Plot Envelope signals and Audio signals on the travel time perturbation axis
% 1. Edit day, start_hour, and end-hour to disply only the duration of interest
% ex day = 27:28; start_hour = 3 end_hour = 2 will display receptions
% during 27 Oct 3:00 to 28 Oct 2:00 (UTC)
% 2. Edit tx_disp to change the number of receptions displayed at once
% 3. scale_signal and scale_snr control the size of the audio signal and the envelope respectively
% 4. All complex Demodulations (Envelope signals) are computed using the
% ideal LFM 22.5-ms pulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
%%%%%%%%
day = 27:27 ;              %  Edit
start_hour = 3;            % Edit
end_hour = 12;             % Edit
Fs = 24000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a file name array
now_hour = start_hour;
now_day = day(1);
fname = [];
fname_tx = [];
while true
    while now_hour <= 23
        if now_hour < 10
            fname_d = "ttp_snr_sig_data_2018_10_"+string(now_day)+"_0"+string(now_hour);
            fname_txd =  "tx_data_2018_10_"+string(now_day)+"_0"+string(now_hour);
        else
            fname_d = "ttp_snr_sig_data_2018_10_"+string(now_day)+"_"+string(now_hour);
            fname_txd =  "tx_data_2018_10_"+string(now_day)+"_"+string(now_hour);
        end
       fname = vertcat(fname,fname_d);
       fname_tx = vertcat(fname_tx,fname_txd);
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

%% Load inputs
% stack timeseries on top of each other
scale_signal = 0.5;
scale_snr = 0.5;


for i = 1:length(fname)
    f_count = 1;
    tx_disp = 5;                   % Edit
    load(fname(i))
    load(fname_tx(i));
    % store data into variaples
    audi_train = ttp_snr_sig_data.sig;
    env_train = ttp_snr_sig_data.snr;
    act_rep_t = ttp_snr_sig_data.act_arrival;
    est_rep_t = ttp_snr_sig_data.est_arrival;
    tx_t = ttp_snr_sig_data.tx_t;
    env_pk_mag = ttp_snr_sig_data.env_pk;
    env_pk_pos = ttp_snr_sig_data.envpk_pos;
    range_surface =  ttp_snr_sig_data.range/1000;  % km
    heading =  ttp_snr_sig_data.heading;
    [tx_len,~] = size(audi_train);
    lat = tx_data.lat;lon = tx_data.lon;
    if tx_len < tx_disp
        tx_disp = tx_len;
    end
    % DELETE NEGATIVE POSITION
    [rm_ind,~] = find(env_pk_pos(:,1) < 0);
    env_pk_pos(rm_ind,:) =ones(length(rm_ind),2);
    env_pk_mag_plot(rm_ind,:) = NaN;
    % start plotting
    while true
        
        %check if already reached the last index
        if f_count+tx_disp-1 >= tx_len
            tx_num  = (tx_len-tx_disp+1:tx_len);
            draw_plot(audi_train,env_train,tx_num,env_pk_mag,env_pk_pos,range_surface,heading,act_rep_t,lat,lon,scale_signal,scale_snr,tx_len,tx_disp);
            break;
        else
            tx_num = (f_count:f_count+tx_disp-1);
            draw_plot(audi_train,env_train,tx_num,env_pk_mag,env_pk_pos,range_surface,heading,act_rep_t,lat,lon,scale_signal,scale_snr,tx_len,tx_disp);
           f_count = f_count+tx_disp;
        end
           pause
            
    end
    
    pause
end

%% %%%%%%%%%
function draw_plot(audi_train,env_train,tx_num,env_pk_mag,env_pk_pos,range_surface,heading,act_arrival,lat,lon,scale_signal,scale_snr,tx_len,tx_disp)
    ACO_lat=22.738772;                  % Aug 2018
    ACO_lon=-158.006186;                % Aug 2018
    %%% create variables stacking each signal on top of each other
    Audio = bsxfun(@plus,audi_train(tx_num,:)*scale_signal, tx_num');  
    Envelope = bsxfun(@plus, env_train(tx_num,:)*scale_snr, tx_num');
    ssig = tx_num(1);
    lsig = tx_num(end);
    Fs = 24000;
    h_time = -0.03;
    t_time = 0.03;
    ttp_duration = h_time:1/Fs:t_time;
    %%%%%%%%%%%%%%%
    % Clean data points without detected peaks
    env_pk_mag_plot = env_pk_mag(tx_num,:);
    env_pk_pos_plot = env_pk_pos(tx_num,:);
    rm_ind = find(((env_pk_pos_plot(:,1) <=723) & (721 <= env_pk_pos_plot(:,1))) & (env_pk_pos_plot(:,1)==env_pk_pos_plot(:,2)));
    env_pk_mag_plot(rm_ind,:) = NaN;
    env_pk_pos_plot(rm_ind,:) = ones(length(rm_ind),2);
    
    
    %%%%%%%%%%%%%%
    % Plot
    %%%% Surface distance
    
    f = figure(2);
    clf
    f.Units = 'normalized';
    f.Position = [0.05 0.9 0.95 0.95];
    
    subplot('Position',[0.02 0.1 0.2 0.5])
    scatter(range_surface(ssig:lsig),tx_num,'*');
    
    axis tight
    grid on
    grid minor;
    xlabel('Surface DIstance (km)');
    ax1 = gca;
    set(gca,'Xdir','normal');
    set(gca,'YAxisLocation','right');

    % Axis Ticks
    yvalue = [1:5:tx_len-1 tx_len];                    %Edit
    y_label = yvalue;                               
    set(gca,'YTick',yvalue);
    set(gca,'Yticklabel',y_label);
    % Label Timestamps
    t_step = tx_disp/5;                      %Edit
    label_offset = 0.001;             % edit
    date_label = datestr([act_arrival(ssig:t_step:lsig-1)' act_arrival(lsig)'],'HH:MM');        
    text([range_surface(ssig:t_step:lsig-1); range_surface(lsig)]'+label_offset,[ssig:t_step:lsig-1 lsig],date_label,'FontSize',13)
    
    % heading plot
    ax1_pos = ax1.Position;             % position of first axes
    ax2 = axes('Position',ax1_pos,'Color','none');
    hold on
    scatter(ax2,heading(ssig:lsig),tx_num,'rx')
    grid on
    set(ax2,'XAxisLocation','top','YAxisLocation','right')
    set(ax2,'YTick',[]);
    set(ax2,'Yticklabel',y_label);
    x2value = 0:90:360;
    x2_label = x2value;
    set(ax2,'XTick',x2value);
    set(ax2,'Xticklabel',x2_label);
    ax2.XLabel.String = 'Heading(degrees)';
    ax2.XColor = 'r';
    xlim(ax2,[0 360]);
    
    % Tx map
    subplot('Position',[0.03 0.67 0.2 0.3])
    if tx_num(end) <= length(lon)
        scatter(lon(tx_num),lat(tx_num),20,1:length(tx_num),'filled');
        hold on
        scatter(ACO_lon,ACO_lat,200,'pk','filled');
        grid on
        colormap jet
    else
        scatter(lon(tx_num(1):end),lat(tx_num(1):end),20,1:length(lat(tx_num(1):end)),'filled');
        hold on
        scatter(ACO_lon,ACO_lat,200,'pk','filled');
        grid on
        colormap jet
    end
       ylim([22.5 23])
    xlim([-158.25 -157.75])  
    
    %%% Arrival Pattern + Envelope Plot
    subplot('Position',[0.26 0.1 0.7 0.8])
    plot(ttp_duration*1000,Audio,'k')                   % Edit
    hold on 
    plot(ttp_duration*1000,Envelope,'Color',[0 0.6 0.8])                  
    grid on
    axis tight 
    headline = sprintf('Envelopes on Raw Audio Signal\n Date: %s',datestr(act_arrival(ssig),'mm/dd/yyyy'))
    title(headline)
    
    % peak locations
    hold on
    for k  =1:length(tx_num)
        pk_plot = scatter(ttp_duration(env_pk_pos_plot(k,:))*1000,env_pk_mag_plot(k,:)*scale_snr+tx_num(k),'r');
        
    end
    xticks([-30:2:30])
    xlabel('Travel Time Perturbation (ms)')
    line([0 0],[tx_num(1) tx_num(end)],'Color','k')    
    clc
end
