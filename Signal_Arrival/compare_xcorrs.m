% Look at cross correlations of given audio files with corresponding tx
% file inputs
clearvars
close all

ACO_lat=22.738894;               %22 44.324'
ACO_lon=-158.006009;              %158 00.372'

%% Load Tx times
cd /Volumes/LaCie_RAP_Backup/RAP/June2018Cruise/Tx_Rx_Output/Path/Radial/Eastward_2

d=dir('tx*.mat');
t_tx =[];
lat_tx = [];
lon_tx = [];
z_offset = [];

% Concatenate tx variables
for i=1:length(d)
    load(d(i).name)
    t_tx=horzcat(t_tx,tx_data.t);                           % EDIT
    lat_tx=horzcat(lat_tx,tx_data.lat);
    lon_tx=horzcat(lon_tx,tx_data.lon);
    z_offset=horzcat(z_offset,tx_data.altitude);
end


%Boat distance from ACO
for i=1:length(lat_tx)
    x_dist(i)=dist([ACO_lat lat_tx(i)],[ACO_lon lon_tx(i)]);
end

%Remove data outside theoretical range (29.5 km)
remove_pts=find(x_dist>29400);
lat_tx(remove_pts)=[];
lon_tx(remove_pts)=[];
x_dist(remove_pts)=[];


%% Find Cross correlation + demodulation+ peak of the envelope
cd /Volumes/LaCie_RAP_Backup/RAP/June2018Cruise/wav_data/Hydrophone_Raw/Hydrophone_Data/Radial/Eastward_2
d=dir('180621-*');
SNR=[];
addpath('/Volumes/LaCie_RAP_Backup/RAP/June2018Cruise/Script/Signal_Arrival')
for i=1:2:(length(d)-1)
    mat_name = d(i).name;
    wav_name = d(i+1).name;
    %Load Audio Data
    [y,t_date] = hyd_audio_prep(mat_name,wav_name);
    
    fname = mat_name;
    % time subtracted starting date/time to the scale of minute
    file_start_time = datenum(2018,str2double(fname(3:4)),str2double(fname(5:6)),str2double(fname(8:9)),str2double(fname(10:11)),0);
    t_num=(t_date-file_start_time).*(3600*24);
  
    %%%%%%Cross Correlation%%%%%%
    fprintf('Cross Correlation \n')
 
    [demod,demod_pos,demod_max]=ACO_cross_correlation_sample(y);
    
    % SNR criterion to filter out bouncers
    SNR=10*log10(demod_max/median(demod));
    ind_remv = find(SNR <= 13);
    demod_pos(ind_remv) = [];
    demod_max(ind_remv) = [];
    
    %Arrivals
    arrivals=t_date(demod_pos);    % actual arrival timestamps
    
    
    % Arrival time criterion to filter out bouncers
    % eliminate bottom bouncer / loop over arrivals / 
    
    start_time = datenum('21-Jun-2018', 'dd-mmm-yyyy');     % EDIT
    start_hour = 11;                                        % EDIT 11 am
    start_time = start_time + start_hour/24;                % fraction of day
    dum_arrivals = arrivals - start_time;                   % sec from start hour
    
if ~isempty(arrivals )  
        
    not_end = 1;
    ar_ind_to_remove = [];
    step = 1;
    
    while not_end == 1
        % cal the interval of 2 contiguous arrivals
        diff_int = (dum_arrivals(step+1) - dum_arrivals(step))*(3600*24);
        % check wheter to keep it
        if diff_int < 20 || diff_int > 32
            %remove that arrival
            ar_ind_to_remove(end+1) = step+1;
        end
        % check if the end of arrival matrix is reached
        if step == length(dum_arrivals)-1
            not_end = 0        
        end
        step =step+1;
        
    end
    
    % Remove bottom bouncer arrivals
    arrivals(ar_ind_to_remove) = [];
    demod_pos(ar_ind_to_remove) = [];
    demod_max(ar_ind_to_remove) = [];
    
    %}
    
    % Remove rx time without tx time
    ar_ind_to_remove = [];
    tt_time = []; 
     for o =1:length(t_tx)
        tt_time(o,:) = abs(arrivals - t_tx(o))*(3600*24);
     end
     
     size_tt = size(tt_time);
     
     % Loop along the column of dum_arrivals to find which index does not have corresponding tx time
     
     for   kk =1:size_tt(2) 
         if isempty(find(tt_time(:,kk) < 16))      % Remove if there is no neighbpring tx time within 16secs apart
            ar_ind_to_remove(end+1) = kk;
           
         end
     end    

     arrivals(ar_ind_to_remove) = [];
     demod_pos(ar_ind_to_remove) = [];
     demod_max(ar_ind_to_remove) = [];
    
    % Plot
    % loop over arrivals
    % window demod signal for each arrival
    
    
    for ii=1:length(arrivals)

        figure(4)
        clf
         subplot(2,1,1)
            
            plot((t_num(demod_pos(ii)-200:demod_pos(ii)+500)-0),y(demod_pos(ii)-200:demod_pos(ii)+500))
            title(sprintf('%s; Range: %d; ',datestr(file_start_time),x_dist(ii)))
            grid on
            hold on
            scatter(t_num(demod_pos(ii)),y(demod_pos(ii)),'r')
            %set(gca,'XLim',[t_num(demod_pos(ii)-200) t_num(demod_pos(ii)+500)])
            ylabel('Raw Signal')
            timetick = get(gca,'Xtick');
            l_ttick = min(timetick);
            r_ttick = max(timetick);
            timetick = l_ttick:0.002:r_ttick;
            set(gca,'XTick',timetick)
            set(gca,'Xticklabel',timetick)
            axis tight
            
          subplot(2,1,2)
            
          
            plot((t_num(demod_pos(ii)-200:demod_pos(ii)+500)-0),demod(demod_pos(ii)-200:demod_pos(ii)+500))
            title(sprintf('SNR: %3f',SNR(ii)))
            grid on
            hold on
            ylabel('Envelope')
            scatter(t_num(demod_pos(ii)),demod(demod_pos(ii)),'r')
            %set(gca,'XLim',[t_num(demod_pos(ii)-200) t_num(demod_pos(ii)+500)])
            set(gcf,'Position',[100 150 900 600])
            timetick = get(gca,'Xtick');
            l_ttick = min(timetick);
            r_ttick = max(timetick);
            timetick = l_ttick:0.002:r_ttick;
            set(gca,'XTick',timetick)
            set(gca,'Xticklabel',timetick)
            axis tight
            %saveas(gcf,[pwd sprintf('/Xcorr_figures/Xcorr%d.jpg',arrival_pos)])
            pause
        end
        
    end
    
end










