% LFM pulse design

close all
clear

% Sampling period
fs=44100;
dt=1/fs;

% Carrier frequency   default
fc=4000;  


% Period
Tp1=0.1;
Tp2=0.08;
Tp3=0.01;

% signal samples
length_sig1=Tp1*fs;
length_sig2=Tp2*fs;
length_sig3 =Tp3*fs;

t1 = 0:Tp1/length_sig1:Tp1;
t2 = 0:Tp2/length_sig2:Tp2;
t3 = 0:Tp3/length_sig3:Tp3;
sig1 = cos(2*pi*fc*t1);
sig2 = cos(2*pi*fc*t2);
sig3 = cos(2*pi*fc*t3);

% zero pad
z_sig1 =[zeros(1,ceil(length_sig1/3)) sig1 zeros(1,floor(length_sig1/3))];
z_sig2 =[zeros(1,ceil(length_sig2/3)) sig2 zeros(1,floor(length_sig2/3))];
z_sig3 =[zeros(1,ceil(length_sig3/3)) sig3 zeros(1,floor(length_sig3/3))];


S1 = fft(z_sig1);
S1 = abs(S1)./length(z_sig1);
S1(2:end) = 2*S1(2:end);
S1 = S1(1:end/2);
f_step1 = fs/length(z_sig1);
bin_f1 = 0:f_step1:fs/2-f_step1;

S2 = fft(z_sig2);
S2 = abs(S2)./length(z_sig2);
S2(2:end) = 2*S2(2:end);
S2 = S2(1:end/2);
f_step2 = fs/length(z_sig2);
bin_f2 = 0:f_step2:fs/2-f_step2;

S3 = fft(z_sig3);
S3 = abs(S3)./length(z_sig3);
S3(2:end) = 2*S3(2:end);
S3 = S3(1:end/2);
f_step3 = fs/length(z_sig3);
bin_f3 = 0:f_step3:fs/2-f_step3;
%% plot amplitude spectrum

[B,I] = min(abs(fc-bin_f1));
amp_c = S1(I+1);
ind_fcut1 = find(10*log10(S1(1:end/5)./amp_c) < -3,1,'last');
cutoff_f1 = bin_f1(ind_fcut1)
ind_fcut2 =floor(length(S1)/5)+find(10*log10(S1(end/5:end/4)./amp_c) < -3,1,'first');
cutoff_f2 = bin_f1(ind_fcut2)

figure(2)
% plot(bin_f1,10*log10(S1./amp_c))
stem(bin_f1,S1/amp_c)
grid on
hold on
% line([fc fc],[-20 2],'Color','k')
set(gca,'XScale','log')
xlim([3000 5000])
ylabel('Amplitude')
%%
hold on

stem(bin_f2,S2/amp_c)
stem(bin_f3,S3/amp_c)

