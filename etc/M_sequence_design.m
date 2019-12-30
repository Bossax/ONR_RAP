% M sequence signal design

close all
clear

% Sampling period
fs=44100;
dt=1/fs;

% Carrier frequency   default
fc = 4134.375;  



law = [ 1 0 0 0 0 0 0 1 0 0 1];    % 2011 oct (2018 RAP)
pp = 2;                            % (GF)
nr = length(law)-1;                % degree
L = pp^nr-1;                       % number of digits

Q = 3;                             % cyc_dig = Q;
B = fc/Q;                           % chiprate

fsx = fs;                   % xducer sample rate 
sam_cyc = fsx/fc;
Tp= sam_cyc*Q*L/fsx;

[sig, law]=msqgen_edit( Q, nr, sam_cyc);      

length_sig = length(sig);
% zero pad
z_sig =[zeros(1,ceil(length_sig/3)) sig' zeros(1,floor(length_sig/3))];

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


amp_c = max(S)
ind_fcut1 = find(10*log10(S(1:end/5)./amp_c) < -3,1,'last');
cutoff_f1 = bin_f(ind_fcut1)
ind_fcut2 =floor(length(S)/5)+find(10*log10(S(end/5:end/4)./amp_c) < -3,1,'first');
cutoff_f2 = bin_f(ind_fcut2)
%%
figure(2)
clf
plot(bin_f,10*log10(S./amp_c)+3)
grid on
hold on

set(gca,'XScale','log')
xlim([3000 6000])
ylim([-20 5])
ylabel('Amplitude')
