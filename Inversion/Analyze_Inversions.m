% This script is executed section-wise
% The first section is to compare 2 inversion solution datasets
% The second section is to plot reconstructed measurement d^ from the observation matrix and the model vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1st section
%% save first month data
m_recov1_o = m_recov_o(1:tot_grid_num);
m_recov2_o = m_recov_o(tot_grid_num+1:2*tot_grid_num);
m_recov3_o = m_recov_o(2*tot_grid_num+1:3*tot_grid_num);

recov_SS1o=reshape(m_recov1_o(1:end),grid_num,grid_num)'  ;
recov_SS2o=reshape(m_recov2_o(1:end),grid_num,grid_num)'  ;
recov_SS3o=reshape(m_recov3_o(1:end),grid_num,grid_num)'  ;


%% save second month data
m_recov1_n = m_recov_n(1:tot_grid_num);
m_recov2_n = m_recov_n(tot_grid_num+1:2*tot_grid_num);
m_recov3_n = m_recov_n(2*tot_grid_num+1:3*tot_grid_num);

recov_SS1n=reshape(m_recov1_n(1:end),grid_num,grid_num)'  ;
recov_SS2n=reshape(m_recov2_n(1:end),grid_num,grid_num)'  ;
recov_SS3n=reshape(m_recov3_n(1:end),grid_num,grid_num)'  ;
%% Plot SS perturbation map
figure(2)
clf
set(gcf,'Units','normalized','Position',[0. 0.9 .7 0.3]);
ax1 = axes;
set(ax1,'Units','normalized','Position',[0.03 0.07 .26 0.68]);
imagesc(x_cen,y_cen,recov_SS1o-recov_SS1n)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
caxis([-4 4])
% cbar.Ticks = -.5:.1:5;
title('Recovered SS Perturbation Field (1st mode)')
hold on
scatter(icListen_lon,icListen_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal')
colormap jet

ax2 = axes;
set(ax2,'Units','normalized','Position',[0.37 0.07 .26 0.68]);
imagesc(x_cen,y_cen,recov_SS2o-recov_SS2n)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title('Recovered SS Perturbation Field (2nd mode) ')
hold on
scatter(icListen_lon,icListen_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal')
caxis([-2 2])
colormap jet

icListen_pos = sprintf('Lat = %.5f\nLon = %.5f\n Depth =%.2f m',icListen_lat,icListen_lon,icListen_depth);
% hyd_off = sprintf('HEM: Lat = %.5f Lon = %.5f Depth =%.2f m\nHydrophone Offset: x = %.3f m, y = %.3f m, z = %.3f m',icListen_lat,icListen_lon,icListen_depth,m_p(1),m_p(2),m_p(3));
t = annotation('textbox',[.3 .88 .1 .1],'String','Sound Speed Perturbation Difference 2nd iteration - 1st iteration (icListen)','Fontsize',12,'BackgroundColor','white');

ax3 = axes;
set(ax3,'Units','normalized','Position',[0.7 0.07 .26 0.68]);
imagesc(x_cen,y_cen,recov_SS3o-recov_SS3n)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title('Recovered SS Perturbation Field (3rd mode) ')
hold on
scatter(icListen_lon,icListen_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal')
caxis([-2 2])
colormap jet

% Depth average diff
figure(21)
clf
set(gcf,'name','2-D sound speed perturbation field','Units','normalized','Position',[0 0.1 0.35 .45])
imagesc(x_cen,y_cen,recov_SS_d_avgo - recov_SS_d_avgn)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Difference (m/s)';
% title('Difference between the origianl and the spatially-filtered measurements')
title({'Difference Sound Speed_{Nov} - Sound Speed_{Oct}','HEM'})
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',13)
caxis([-.05 .05])
xlabel('Long')
ylabel('Lat')
%% 2nd Section
%% Use with SS_inversion to plot reconstructed measurements
% Calculate Azimute
azmth = ones(length(tx_lat),1);
% 0-360

for i=1:length(tx_lat)
    azmth(i) = azimuth(icListen_lat,icListen_lon,tx_lat(i),tx_lon(i));
end
theta = tx_heading';
theta = tx_heading' - azmth;

for i = 1:length(theta)
   if theta(i) < 0 
       theta(i) = theta(i)+360;
   end
    
end

% plot reconstructed measurement from ocean perturbations
figure(3)
clf
set(gcf,'Units','normalized','Position',[0 0.5 0.4 0.5])
grid on
scatter(tx_lon,tx_lat,20,d_recov_ocean*1000,'filled')
axis tight
grid on
colormap jet
cbar = colorbar;
cbar.Label.String = 'ms';
title('Reconstructed TTP Map from Ocean Perturbations (icListen)')
xlabel('Long')
ylabel('Lat')
set(gca,'fontsize',13)
caxis([-3  3])

% Plot hydrophone offset and ocean perturbation 
d_recov_ocean = G(:,1:end-3)*m_recov(1:end-3);
d_recov_hyd = G(:,end-2:end)*m_recov(end-2:end);
d_recov = G*m_recov;
figure(4)
set(gcf,'Units','normalized','Position',[0.4 0.5 0.4 0.8])
subplot(3,1,1)

scatter(range,d_recov_ocean*1000,10,azmth,'filled')
colormap jet
ylim([-5 5])
grid on
yticks(-5:5)
xlim([0 26])
xlabel('Range (km)')
ylabel('TTP (ms)')
title('Reconstructed Measurements: Ocean Features (icListen)')
set(gca,'fontsize',12)
cbar = colorbar;
cbar.Label.String = 'Azimuth (degrees)' ;


subplot(3,1,2)
scatter(range,d_recov_hyd*1000,10,azmth,'filled')
colormap jet
ylim([-5 5])
grid on
yticks(-5:5)
xlim([0 26])
xlabel('Range (km)')
ylabel('TTP (ms)')
title('Reconstructed Measurements: Hydrophone Offsets (icListen)')
set(gca,'fontsize',12)
cbar = colorbar;
cbar.Label.String = 'Azimuth (degrees)' ;

subplot(3,1,3)
scatter(range,d_recov*1000,10,azmth,'filled')
colormap jet
ylim([-5 5])
grid on
yticks(-5:5)
xlim([0 26])
xlabel('Range (km)')
ylabel('TTP (ms)')
title('Reconstructed Measurements (icListen)')
set(gca,'fontsize',12)
cbar = colorbar;
cbar.Label.String = 'Azimuth (degrees)' ;
