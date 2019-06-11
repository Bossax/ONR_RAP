%% Plot distance vs time with log(xcorr)
clearvars
close all

%% PART 1

cd ../../..
cd Tx_Rx_Output/variable/Westward

% Load ray parameters
for i =1:1
   
     load ray_data_westward
end

%% Load cross correlation function and concatenate them together into a
% single matrix
tot_dist = ray_data.ray_arc_dist;
arc_length = ray_data.arc_length;
x_dist = ray_data.surface_dist;
d=dir('xcorr*.mat');

%num_pos=0;
xcorr_all=[];
for i=1:1length(d)

    load(d(i).name);        %Load xcorr
    
    %num_pos(i+1)=num_pos(i)+size(xcorr_tx,2);
    
    %get_pos=find(B>num_pos(i) & B<=num_pos(i+1));
    
    xcorr_all=horzcat(xcorr_all,xcorr_tx(1:2:end,:));           % sample every 2 
    %xcorr_all=horzcat(xcorr_all,xcorr_tx(:,B(get_pos)-num_pos(i)));
    
    % clear xcorr_tx
end

%Remove every (50) redundant transmission from the same distance
% A=horzcat((1:130),(131:50:1169),(1170:1199),(1200:50:1864),(1865:1905),(1906:50:2233),(2234:2253),(2254:50:2400),(2401:2514),(2515:50:end));
% tot_dist=tot_dist(A);
% xcorr_all=xcorr_all(:,A);

%% Sorting from mow to high

[~,I] = sort(tot_dist);
tot_dist = tot_dist(I);
xcorr_all = xcorr_all(:,I);

%% Get SNR to plot for xcorr
noise=median(median(abs(xcorr_all)));
SNR=abs(xcorr_all)./noise;
SNR_dB=10*log10(SNR./noise);
a=find(SNR_dB<-20);
b = find(SNR_dB > 40);
SNR_dB(a)=NaN;      % discard SNR less than -20 dB 
SNR_dB(b)=NaN;      % discard SNR greater than 40 dB
clear SNR
%Remove distances that out of bound
s = size(xcorr_all);
tot_dist(s(2)+1:end) = [];

%% Plotting
%surf(x_dist./1000,linspace(0,30,length(xcorr_all)),SNR_dB,'EdgeColor','None')
%imagesc(linspace(0,16,size(xcorr_all,2)),linspace(0,30,length(xcorr_all)),SNR_dB)
surf(tot_dist./1000,linspace(0,30,length(xcorr_all)),SNR_dB,'EdgeColor','None')
view(2)

set(gca,'Xlim',[min(tot_dist./1000) max(tot_dist./1000)])
xlabel('Ray Arc Length (km)')
%xlabel('Surface Distance (km)')
ylabel('Travel Time (s)')
c = colorbar;
c.Label.String = 'SNR (dB)';
caxis([-20 20]);  % Tune
colormap jet
set(gcf,'Position',[100 100 1000 750])
ytick = yticks;
max_tick = max(ytick);
yvalue = 0:1:max(ytick);
set(gca,'YTick',yvalue);
set(gca,'Yticklabel',yvalue);
set(gca,'YLim',[0 20]) ;
figname = 'Westward radial SNR'
title(figname)
axis tight
saveas(gcf,[figname,'.jpg'])










