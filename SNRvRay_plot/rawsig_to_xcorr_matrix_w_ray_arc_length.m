% Extract hydrophone data for each transmission cycle for SNR plot % Boris Aug2018
% 1. seperate individual cycle of transmissions and package all cycles of transmissions in a matrix y_tx
% whose column contains acoustic waveform data of each transmission cycle. Each cycle is 30-second long (2018)
% 2. produce ray arc parameters (total ray length, layerwise surface and ray arc distances and)and save the output file
% 3. produce a croxs-correlation matrix of signal of each transmission cycle
% 4. produce surface distance and ship altitude output files
clearvars
close all

ACO_lat=22.738772;               % Update Aug 2018
ACO_lon=-158.006186;              

%% Pack tx parameters
cd /Volumes/LaCie_RAP_Backup/RAP/June2018Cruise/Tx_Rx_Output/tx_files

t_tx = [];
lat_tx = [];
lon_tx = [];
z_offset = [];
% Edit
file_num = 65:67;

for i=1:length(file_num)
    load(['tx_data',int2str(file_num(i)),'.mat'])
    t_tx = horzcat(t_tx,tx_data.t);
    lat_tx = horzcat(lat_tx,tx_data.lat);
    lon_tx = horzcat(lon_tx,tx_data.lon);
    z_offset =horzcat(z_offset,tx_data.altitude);
end

z_offset=z_offset-2.31;     % geoid compensation

% Correct t_tx for Scarlett offset
t_tx = t_tx +3/(3600*24);           % Aug 2018

%% Calculate ship distance
cd ..

%Boat distance from ACO
for i=1:length(lat_tx)
    x_dist(i)=dist([ACO_lat lat_tx(i)],[ACO_lon lon_tx(i)]);
end

%Remove data outside theoretical range (29.5 km)
remove_pts=find(x_dist>29400);
t_tx(remove_pts)=[];
lat_tx(remove_pts)=[];
lon_tx(remove_pts)=[];
x_dist(remove_pts)=[];

z_offset(remove_pts)=[];

t_tx_all=t_tx;

%% Loop over tx times and split them into group of 300 transmissions

%split_num=300;
%num_hold=0;

%num_hold=4;
%for k=1170:split_num:length(t_tx_all)
    %num_hold=num_hold+1;
    
    clearvars -except t_tx_all split_num k num_hold t_pulse x_dist z_offset

%{   
if k+split_num>length(t_tx_all) % the remaining section in the array
t_tx=t_tx_all(k:end);    
else
t_tx=t_tx_all(k:k+split_num-1); % the next 300 timestamps
end
%}
    
%% Load audio files

cd /Volumes/LaCie_RAP_Backup/RAP/June2018Cruise/wav_data/Hydrophone_Raw/Hydrophone_Data/Radial/Westward

d=dir('18*');

y_all=[];
t_date_all=[];

for i=1:2:length(d)

    %Load Audio Data
    [y,fs]=audioread(d(i+1).name);  
    M1=matfile(d(i).name);

    %% Timing Variables
    time_exact=M1.Date_num_count-(.00005/(3600*24));        %50 us lag between 96 and 26kHz data
    time_interval=diff(time_exact).*(3600*24);              %Find time intervals between headers
    time_interval(time_interval<0)=[];
    mean_time_interval=mean(time_interval);     %Mean time interval
    single_time_interval=mean_time_interval/4096;      %Single time interval
    mean_date_interval=mean_time_interval/(3600*24);
    single_date_interval=single_time_interval/(3600*24);
    time_diff=(time_exact(end)-time_exact(1))/(1/(3600*24))+mean_time_interval;

    t_date=[];
    t_date(1:512)=linspace(time_exact(1)-(mean_date_interval/8),time_exact(1),512);
    for jj=2:length(time_exact)
        t_date(((jj-2)*4096+513):((jj-1)*4096+512))=linspace((time_exact(jj-1)+single_date_interval),time_exact(jj),4096);
    end
    t_date(end+1:end+3584)=linspace((time_exact(end)+single_date_interval),(time_exact(end)+(7*mean_date_interval/8)),3584);
    y=y(1:length(t_date));              %Fix y (if needed/should no longer be necessary)
    
    y_all=vertcat(y_all,y);                 % Pack audio data 
    t_date_all=vertcat(t_date_all,t_date'); % pack timing data of the corresponding set of audio files

end

%% Loop over tx times and find the closet timestamp in the packed audio file coressponding to individual tx time
t_int = 30;          % interval of each transmission cycle
num_pts=t_int*fs;   % Number of samples for 30 sec interval

% Remove tx timestamps which donot belong to the audio file time frame
rm_ind = find(t_tx_all < t_date_all(1));
t_tx_all(rm_ind) = [];
x_dist(rm_ind) = [];
z_offset(rm_ind) = [];
rm_ind = find(t_tx_all > t_date_all(end));
t_tx_all(rm_ind) = [];
x_dist(rm_ind) = [];
z_offset(rm_ind) = [];

y_tx=zeros(num_pts,length(t_tx_all)); % pre-allocate space


for i=1:length(t_tx_all)
    [~,B]=min(abs(t_tx_all(i)-t_date_all)); % 
    if B+num_pts-1 <= length(y_all)
        y_tx(:,i)=y_all(B:B+num_pts-1);
    else
        zeropad =  num_pts - (length(y_all)-B+1);
        y_tx(:,i)=vertcat(y_all(B:end),zeros(zeropad,1));
    end
end
%% Calculate and save ray larc length
cd ../../../../..
cd Tx_Rx_Output/variable/Westward

ray_arc_length(t_tx_all,x_dist,z_offset);



%% Calculate xcorr for truncated signal matrix y_tx

 load('LFM_pulse_sample.mat')
    
    for ii=1:size(y_tx,2)
       
        y_tx(:,ii)=y_tx(:,ii)./(max(abs(y_tx(:,ii))));
        
        %%Cross Correlation%%
        [xc,lag]=xcorr(y_tx(:,ii),replica_pulse);
        xc(1:find(lag==0)-1)=[];
        xcorr_tx(:,ii)=xc;
        
    end    
    
%% SAVE surface distance and altitude
save('y_tx_westward','y_tx')        % EDIT                                    
save('x_dist_westward','x_dist')
save('z_offset_westward','z_offset')
 save('xcorr_tx_westward','xcorr_tx')
%% 
function ray_arc_length(t_tx_all,x_dist,z_offset)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate ray arc length
%%%% Inputs %%%%
% 1. tx files
% 2. x dist
% 3. z offset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main
arc_lengths = [];
tot_dist = [];
surface_dist = [];

for i=1:length(x_dist)
    [arc_handle,tot_handle,~,~,~,~,surface_dist_handle,~,~]=ray_trace_w_curvature(x_dist(i),z_offset(i));
    arc_lengths(end+1,:) = arc_handle;
    tot_dist(end+1) = tot_handle;
    surface_dist(end+1,:) = surface_dist_handle;

end

 % Save

 ray_data.ray_arc_dist = tot_dist;
 ray_data.arc_length = arc_lengths;
 ray_data.surface_dist = surface_dist;
 save('ray_data_westward.mat','ray_data')


end



