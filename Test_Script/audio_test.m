clear
close all
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/wav_data/HEM/28_Oct
fs = 24000;
[y,t_date] = hyd_audio_prep('181028-124000.mat','181028-124000.wav');

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

% Filtering
yf = filter(b,1,y);
delay = round(mean(grpdelay(b,1,6000,Fs)));      % group delay
yf = yf(delay+1:end);                            % shift output signal forward in time
yf(end+1:end+1+delay-1) = zeros(delay,1);
%% Plot a signal
figure(1);
clf
ax1 = axes();
t = (t_date-t_date(1))*3600*24;
line1 = line(ax1,t,yf/max(yf),'LineWidth',2);
line([0 t(end)],[0 0],'Color','k','LineWidth',1)
ax1.XAxisLocation = 'bottom';
grid on
axis tight
xlabel('Second')

% Cross Correlation
rep = load('LFM_pulse_ideal.mat');
rep = rep.replica_pulse;


[xc,lag]=xcorr(yf,rep);
xc(1:find(lag==0)-1)=[];    % set correlation values to zero (no negative time lag)
demod=abs(hilbert(xc));     %Complex Demod of signal (Magnitude of a complex signal)
lag = lag(find(lag ==0):end);   % keep only positive lag


figure(1)
hold on
plot(t,demod/max(demod),'m','LineWidth',0.2)

%%
plot(t,xc/max(xc),'g--','LineWidth',0.1)




