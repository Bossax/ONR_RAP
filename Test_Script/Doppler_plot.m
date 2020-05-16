% Doppler plot
% 1. read in 6 audio files of 3 different velocities (0, spin, 3in/out 5in/out m/s)
% 2. build a starndard LFM function
% 3. specify f shift for each iteration and the range of doppler test
% 4. loop over each audio file pick the 5th transmission (arbritarily) to present its corss-correlation function
% 5. plot 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
close all

au_fname1 = '181028-110000'; % 0m/s 5km circle
au_fname2 = '181028-140000'; % spin 5 km
au_fname3 = '181028-172000'; % 5m/s outward
au_fname4 = '181028-200000'; % 3m/s inward
au_fname5 = '181028-220000'; % 3m/s outward
au_fname6 = '181029-004500';  % 5m/s inward

% original LFM
fs=24000;       % HEM
dt=1/fs;
fc=4134.375;    %carrier frequency   default
B=1378.125;     %Bandwidth
Tp=0.0225;
length_sig=Tp*fs;
sigo=chirp(0:Tp/length_sig:Tp,(fc-B/2),Tp,(fc+B/2),'linear',-90)';

% F range
frange = 50;
f_step = -50:5:50;
fc_test = fc +f_step;

% initiate variables for plots
duration = 0.01;    % 20 msec
demod_f = zeros(length(f_step),duration*fs*2+1);    % Matrix to contain demod for each  f

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
%%
%%%%%%%%%%%%%%%%%%%%
% choose 
au_fname = string(au_fname5);
%%%%%%%%%%%%%%%%%%%%
% load audio file
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/wav_data/HEM/28_Oct

mat_name = au_fname+".mat";
wav_name = au_fname+".wav";
% Timing audio data 
[y,t_date] = hyd_audio_prep(mat_name,wav_name);
datestr(t_date(end))
 % Filtering
y = filter(b,1,y); 
delay = round(mean(grpdelay(b,1,6000,Fs)));      % group delay
y = y(delay+1:end);                       % shift output signal forward in time
y(end+1:end+1+delay-1) = zeros(delay,1); 


% find the original arrival time
[demod,demod_pos,demod_max,demod_snr] = ACO_cross_correlation_sample(y);

ind0 = demod_pos(2);
t0 = t_date(ind0);
%% Loop over f range 
for ii = 1:length(f_step)
    % create replica pulse
    replica_pulse=chirp(0:Tp/length_sig:Tp,(fc_test(ii)-B/2),Tp,(fc_test(ii)+B/2),'linear',-90)';
    
    % cross correlation
    [demod,demod_pos,~,~] = ACO_cross_correlation(y,replica_pulse);
    
    % truncate demod
    indstart = ind0 - duration*fs;
    indend= ind0 + duration*fs;
    demod_trunc = demod(indstart:indend);
    
    % store
    demod_f(ii,:) = demod_trunc;
end

%% plot
scale = 0.75e-7
scale_offset = 5;
demod_f_plot = bsxfun(@plus,demod_f*scale, transpose(-(length(f_step)-1)/2:(length(f_step)-1)/2)*scale_offset);  
%%
fontsize = 12;
t_sec = (-duration:1/fs:duration)*1000; % msec
f1 = figure(1);
f1.Units = 'normalized';
f1.Position = [0.3 0.5 0.5 0.7];
clf
plot(t_sec,demod_f_plot(1:10,:),'b')
hold on
plot(t_sec,demod_f_plot(11,:),'r')
plot(t_sec,demod_f_plot(12:21,:),'b')
grid on

% adjust y axis labl
set(gca,'Ytick',f_step)
set(gca,'YtickLabel',f_step)
set(gca,'fontsize',fontsize)
ylim([-80 80])
xticks(-10:1:10)

line([0 0],[-80 80],'linewidth',2,'Color','k')
xlabel('\tau (msec)')
ylabel('\deltaf (Hz)')
title('Doppler check (5 m/s inward)')

figure(3)
plot(t_date,y)
datetick('x')
axis tight
grid on
%%
 function [demod,demod_pos,demod_max,demod_snr] = ACO_cross_correlation(y,replica_pulse)
%Create pulse replica for time series
% return
fs=24000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%Cross Correlation%%
[xc,lag]=xcorr(y,replica_pulse);
xc(1:find(lag==0)-1)=[];    % set correlation values to zero (no negative time lag)
demod=abs(hilbert(xc));     %Complex Demod of signal (Magnitude of a complex signal)
lag = lag(find(lag ==0):end);   % keep only positive lag

%{
figure(1)
plot(lag,xc)
axis tight
grid on
title(sprintf('Cross-Correlation: Signal length = %f \n Lag = %f',length(y),length(lag)))
xlabel('Lag')
%}

%threshold=mean(demod)+10*std(demod);
threshold=9;
%threshold=6;                %dB threshold
SNR=10*log10(demod/median(demod));  % / noise
threshold2=max(SNR);
threshold2=threshold2*(3/5);  % scaling??

if threshold2>threshold
    threshold=threshold2;
end
n=find(SNR>threshold);            %Find large values for demod xcorr
%%%%


%% Get each xcorr arrival %%
demod_separation=[];

%%Separate each pulse cross correlation%%

for i=2:length(n)
    if n(i)>(n(i-1)+(fs*2))     %check if it's within 2 seconds
        demod_separation(end+1)=i;
    end
end
%%%%


if isempty(demod_separation)
    demod_max=[];
    demod_pos=[];
    demod_snr = [];
   
else
    for i=1:length(demod_separation)+1
        if i==1
            if n(1)-200>0
                [pks,locs]=findpeaks(demod(n(1)-200:n(demod_separation(i)-1)));
            else
                [pks,locs]=findpeaks(demod(1:n(demod_separation(i)-1)+5));
            end
            if isempty(pks)==1
                demod_max(i)=demod(1);
                demod_snr(i) = 10*log10(demod(1)/median(demod));
                demod_pos(i)=1;
            else
                
                [~,B]=max(pks);
                
                A=find((pks.*3/2)>=pks(B)); % find the neighbor peaks which potentially belongs to the first pulse
                
                demod_max(i)=pks(A(1));
                demod_snr(i) = 10*log10(pks(A(1))/median(demod));
                if n(1)-200>0
                    demod_pos(i)=n(1)-200+locs(A(1))-1;
                else
                    demod_pos(i)=locs(A(1));
                end
            end
            
        elseif i==length(demod_separation)+1
            [pks,locs]=findpeaks(demod(n(demod_separation(i-1))-200:n(end)));
            [~,B]=max(pks);
            
            A=find((pks.*3/2)>=pks(B));
            
            demod_max(i)=pks(A(1));
            demod_snr(i) = 10*log10(pks(A(1))/median(demod));
            demod_pos(i)=n(demod_separation(i-1))-200+locs(A(1))-1;
            
        else
            [pks,locs]=findpeaks(demod(n(demod_separation(i-1))-200:n(demod_separation(i)-1)));
            [~,B]=max(pks);
            
            A=find((pks.*3/2)>=pks(B));
            
            demod_max(i)=pks(A(1));
            demod_snr(i) = 10*log10(pks(A(1))/median(demod));
            demod_pos(i)=n(demod_separation(i-1))-200+locs(A(1))-1;
        end
        
    end
end
%{
figure(2)
plot(demod)
axis tight
grid on
title('Envelope Signal')
hold on
scatter(demod_pos,demod(demod_pos))
  %}  
    


end

