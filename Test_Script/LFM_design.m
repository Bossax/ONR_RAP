% LFM pulse design

close all
clear

% Sampling period
fs=44100;
dt=1/fs;

% Carrier frequency   default
fc=4134.375;  

% Bandwidth
B=1378.125;     

% Period
Tp=0.0225;

% signal samples
length_sig=Tp*fs;

% chirp() generates cosine wave
% sine LFM
t = 0:Tp/length_sig:Tp;
fi_0 = fc-B/2;
fi_1 = fc+B/2;
del_f = B/length_sig; 
fi = fi_0:del_f:fi_1;
sig = chirp(t,fi_0,Tp,fi_1,'linear',-90); 

% zero pad
z_sig =[zeros(1,ceil(length_sig/3)) sig zeros(1,floor(length_sig/3))];

figure(1)
plot(z_sig)
grid on

S = fft(z_sig);
S = abs(S)./length(z_sig);
S(2:end) = 2*S(2:end);
S = S(1:end/2);
f_step = fs/length(z_sig);
bin_f = 0:f_step:fs/2-f_step;

% plot amplitude spectrum

[B,I] = min(abs(fc-bin_f));
amp_c = S(I+1);
ind_fcut1 = find(10*log10(S(1:end/5)./amp_c) < -3,1,'last');
cutoff_f1 = bin_f(ind_fcut1)
ind_fcut2 =floor(length(S)/5)+find(10*log10(S(end/5:end/4)./amp_c) < -3,1,'first');
cutoff_f2 = bin_f(ind_fcut2)

figure(2)
plot(bin_f,10*log10(S./amp_c))
grid on
hold on
line([fc fc],[-20 2],'Color','k')
set(gca,'XScale','log')
xlim([3000 6000])
ylim([-20 2])
ylabel('Amplitude')
scatter(cutoff_f1,S(ind_fcut1))
