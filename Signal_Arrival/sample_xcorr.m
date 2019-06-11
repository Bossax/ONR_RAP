%% Sample to read in audio files and compute cross correlation and plot the result
clearvars
close all

cd '/Volumes/LaCie_RAP_Backup/RAP/June2018Cruise/wav_data/Hydrophone_Raw/Hydrophone_Data/Radial/Eastward_1'
%%%%Desired Audio File (EDIT)%%%%
fname1='180621-111500.wav';
fname2='180621-111500.mat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Load Audio Data
[y,fs]=audioread(fname1);
M1=matfile(fname2);


%Timing Variables
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


%%%%%%Cross Correlation%%%%%%
addpath('/Volumes/LaCie_RAP_Backup/RAP/June2018Cruise/Script/Signal_Arrival');
[demod,demod_pos,demod_max]=ACO_cross_correlation_sample(y);

%Arrivals
arrivals=t_date(demod_pos);


%Plot audio data with arrivals circled
figure(1)
plot(t_date,y)
hold on
grid on
datetick
scatter(t_date(demod_pos),y(demod_pos),'r')
ylabel('Amplitude (uPa)')
xlabel('Time (HH:MM)')
set(gca,'Xlim',[min(t_date) max(t_date)])
title('Arrival points')







