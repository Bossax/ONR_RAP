clear
close all
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/wav_data/HEM/28_Oct
fname = '181028-102000';
[y,t_date] = hyd_audio_prep([fname '.mat'],[fname '.wav']);

%%  Bandpass filter for HEM 2000 Hz to 6000 Hz
fl = 2000;
fh = 10000;
Fs = 24000;
trans_band = 500;  %Hz
fcut = [fl fl fh fh] + [-1 1 -1 1]*trans_band;
mags = [0 1 0];
devs = [0.01 0.05 0.01];
[n,Wn,beta,ftype] = kaiserord(fcut,mags,devs,Fs);
b= fir1(n,Wn);
[h,freq] = freqz(b,1,16000,Fs);

yf = filter(b,1,y); 
delay = round(mean(grpdelay(b,1,6000,Fs)));      % group delay
yf = yf(delay+1:end);                       % shift output signal forward in time
yf(end+1:end+1+delay-1) = zeros(delay,1);
figure(22)
clf
plot(yf)
for iii =1:2
    [x,y] = ginput(2);
    % HEM
    ind1 = x(1);
    ind2 = x(2);
    yf = yf(ind1:ind2);
    clf 
    plot(yf)
end
%% periodogram
nfft = 2^ceil(log2(length(yf)));
w = hann(length(yf));
[pxx,f] = periodogram(yf.*w,[],nfft,24000);
f(2)-f(1)
%% 
figure(1)
plot(f,pxx)
grid on
set(gca,'Xscale','log')
set(gca,'Yscale','log')
xlim([3000 6000])
ylabel('PSD')
xlabel('Frequency')
date = datestr(t_date(round(x(1))),'mm/dd HH:MM')
title({'HEM PSD: ', date})