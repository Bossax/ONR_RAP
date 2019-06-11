%% Find time offset between PPS and computer time for each transmission
% slightly modifed from Vincent's by Bruce 6June2017
% check Scarlette output file
clearvars
clear
close all

fs=44100;                       %Sample Rate (Hz)

%oldfolder = cd('/Volumes/Seagate Backup Plus Drive/RAP_20170601_ROV_Cruise')
%oldfolder = cd('/Users/testuser/Documents/Vincent_RAP/Scarlett_Code/data/')
cd '/Volumes/ACO_RAP_3/RAP/Oct2018Cruise/Scarlett_Data'
%%%%%%%%%%EDIT%%%%%%%%%%
filename = 'LFM_1397081757.mat'
load(filename)%
%load LFM_888024629.mat         %Load Data Structure
%%%%%%%%%%%%%%%%%%%%%%%%

ref_tx=txrx.sam(1,:);           %Ref Transducer
pps=txrx.sam(2,:);              %PPS
current=txrx.sam(3,:);          %Current
voltage=txrx.sam(4,:);          %Voltage
t_axis=txrx.taxis;              %Time (s)


%Gain/Scaling
current=current.*7;             %7 for gain
voltage=voltage.*7.*100;        %7 for gain and 100 for scale
pps=pps.*2.66;                  %2.66 for gain
ref_tx=ref_tx.*2.66;            %2.66 for gain

%Plots
close all
fname = filename

figure(1)
plot(t_axis,ref_tx)
xlabel('Time (s)')
ylabel('Amplitude (V)')
tit = ['Reference Hydrophone  ',fname]
title(tit)
grid on

figure(2)
plot(t_axis,pps)
xlabel('Time (s)')
ylabel('Amplitude (V)')
tit = ['PPS  ',fname]
title(tit)
grid on

figure(3)
plot(t_axis,current)
xlabel('Time (s)')
ylabel('Amplitude (A)')
tit = ['Current Measurement  ',fname]
title(tit)
grid on

figure(4)
plot(t_axis,voltage)
xlabel('Time (s)')
ylabel('Amplitude (V)')
tit = ['Voltage Measurement  ',fname]
title(tit)
grid on
%{
figure(5)
max1 = max(abs(ref_tx))
plot(t_axis,tx_ref/max1)
hold on
max2 = max(abs(pps))
plot(t_axis,pps/max2)
max3 = max(abs(current))
plot(t_axis,current/max3)
max4 = max(abs(voltage))
plot(t_axis,voltage/max4)
xlabel('Time (s)')
ylabel('Amplitude (V)')
tit = ['All Measurement(scaled)  ',fname]
title(tit)
%}





