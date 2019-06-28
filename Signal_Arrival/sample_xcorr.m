%% Sample to read in audio files and compute cross correlation and plot the result
clearvars
close all

cd '/Volumes/ACO_RAP_2/RAP/Oct2018Cruise/wav_data/HEM/27_Oct'
%%%%Desired Audio File (EDIT)%%%%
wav_name='181027-042500.wav';
mat_name='181027-042500.mat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Timing audio data
[y,t_date] = hyd_audio_prep(mat_name,wav_name);



%%%%%%Cross Correlation%%%%%%
 [demod,demod_pos_pks,demod_max,demod_snr] = ACO_cross_correlation_ideal(y);
SNR=10*log10(demod/median(demod));
%Arrivals
arrivals=t_date(demod_pos_pks);


%% Plot audio data with arrivals circled
figure(1)
clf
set(gcf,'Units','normalized','Position',[0 0.6 0.6 0.5])
subplot(2,1,1)
plot(t_date,y)
hold on
grid on
datetick
scatter(t_date(demod_pos_pks),y(demod_pos_pks),'r')
ylabel('Amplitude (uPa)')
xlabel('Time (HH:MM)')
set(gca,'Xlim',[min(t_date) max(t_date)])
title('Arrival points')

subplot(2,1,2)
plot(t_date,SNR)
hold on
grid on
datetick
scatter(t_date(demod_pos_pks),SNR(demod_pos_pks),'r')
ylim([0 20])
ylabel('SNR (dB)')
xlabel('Time (HH:MM)')
set(gca,'Xlim',[min(t_date) max(t_date)])
title('Complex Envelope')
line([t_date(1) t_date(end)],[10 10],'color','r')








