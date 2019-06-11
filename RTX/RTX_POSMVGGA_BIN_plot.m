% Plot RTX/GGA?Bin comparison
close all
clear
cd /Users/testuser/Documents/MATLAB/Script/Data
load('RTX_GGA_BIN_compare_wo_interpolation_aftercruise.mat')

ant = [6.439;6.505;-27.761]; % from KM coordinate system May 2018 spreasheet
tx5 = [0.568;19.599;0.715];
ant2tx = tx5-ant;            % vector pointing from the transducer to the antenna

%%
offset_rtx_posmv_sn3 = offset_posmvgga_ref_gn_s3-offset_posmvbins2_ref_gn_s3;%offset_rtx_ref_gn_s3-offset_posmvbins2_ref_gn_s3; %offset_rtx_ref_gn_s3 - offset_posmvgga_ref_gn_s3 ; 
%% plot

% time frame
tstart1 = datenum('20181031 03:00:02','yyyymmdd HH:MM:SS');
tend1= datenum('20181031 11:00:00','yyyymmdd HH:MM:SS');
[~,indstart1] = min(abs(rtx_t_all - tstart1));
[~,indend1] = min(abs(rtx_t_all - tend1));
tstart2 = datenum('20181031 12:00:02','yyyymmdd HH:MM:SS');
tend2= datenum('20181031 24:00:00','yyyymmdd HH:MM:SS');
[~,indstart2] = min(abs(rtx_t_all - tstart2));
[~,indend2] = min(abs(rtx_t_all - tend2));


fsize = 14;
step_size = 0.2;
datetick1 = datenum('2018-10-31 00:00','yyyy-mm-dd HH:MM');
datetick2 = datenum('2018-10-31 24:00','yyyy-mm-dd HH:MM');
int_tick =3/24;
d_xtick = datetick1:int_tick:datetick2;

dateformat ='HH:MM';
mn_xoff = mean(offset_rtx_posmv_sn3(1,[indstart1:indend1 indstart2:indend2]))
mn_yoff = mean(offset_rtx_posmv_sn3(2,[indstart1:indend1 indstart2:indend2]))
mn_zoff = mean(offset_rtx_posmv_sn3(3,[indstart1:indend1 indstart2:indend2]))
med_xoff = median(offset_rtx_posmv_sn3(1,[indstart1:indend1 indstart2:indend2]))
med_yoff = median(offset_rtx_posmv_sn3(2,[indstart1:indend1 indstart2:indend2]))
med_zoff = median(offset_rtx_posmv_sn3(3,[indstart1:indend1 indstart2:indend2]))
rms_xoff = rms(offset_rtx_posmv_sn3(1,[indstart1:indend1 indstart2:indend2])-med_xoff)
rms_yoff = rms(offset_rtx_posmv_sn3(2,[indstart1:indend1 indstart2:indend2])-med_yoff)
rms_zoff = rms(offset_rtx_posmv_sn3(3,[indstart1:indend1 indstart2:indend2])-med_zoff)
med_xvel = median(posmv_xvel([indstart1:indend1 indstart2:indend2]))
med_yvel = median(posmv_yvel([indstart1:indend1 indstart2:indend2]))
med_zvel = median(posmv_zvel([indstart1:indend1 indstart2:indend2]))

% Standarad deviation
sdx = std(offset_rtx_posmv_sn3(1,[indstart1:indend1 indstart2:indend2]));
sdy = std(offset_rtx_posmv_sn3(2,[indstart1:indend1 indstart2:indend2]));
sdz = std(offset_rtx_posmv_sn3(3,[indstart1:indend1 indstart2:indend2]));

% RTX and GGA offset (RTX - POS MV)
f = figure(1);
f.Units = 'normalized';
f.Position = [0 1 0.6 0.7];
clf
subplot(3,1,1)
scatter(rtx_t_all([indstart1:indend1 indstart2:indend2]),offset_rtx_posmv_sn3(1,[indstart1:indend1 indstart2:indend2]),15,posmv_xvel([indstart1:indend1 indstart2:indend2]),'filled')
grid on
xticks(d_xtick)
datetick('x',dateformat,'keepticks')
ylim([-1 1])
yticks(-1:0.2:1)
ylabel('Bow-Stern offset (m)')
% headline = sprintf('Transducer offset from the RTX antenna \n Mean = %.3f m,Median = %.3f m, RMS offset = %.3f m',mn_xoff,med_xoff,rms_xoff);
headline = sprintf('Antenna Position Difference at-dock after cruise (POSMV GGA - POSMV Binary Sensor 2)  \nMean = %.3f m, Median = %.3f m, RMS = %.4f m',mn_xoff,med_xoff,rms_xoff);
title(headline)
cbar = colorbar;
cbar.Label.String = 'X Vel (m/s)';
colormap jet
line([datetick1 datetick2],[2*sdx+med_xoff 2*sdx+med_xoff],'Color','k')
line([datetick1 datetick2],[-2*sdx+med_xoff -2*sdx+med_xoff],'Color','k')
text(d_xtick(round(length(d_xtick)*0.5)),6*sdx+med_xoff,sprintf('%.4f m',4*sdx),'fontsize',fsize);
set(gca,'fontsize',fsize)

subplot(3,1,2)
scatter(rtx_t_all([indstart1:indend1 indstart2:indend2]),offset_rtx_posmv_sn3(2,[indstart1:indend1 indstart2:indend2]),15,posmv_yvel([indstart1:indend1 indstart2:indend2]),'filled')
grid on
datetick('x')
ylim([-1 1])
yticks(-1:0.2:1)
ylabel('Startboard-Portside offset (m)')
colorbar
colormap jet
cbar = colorbar;
cbar.Label.String = 'Y Vel (m/s)';
headline = sprintf('Mean = %.3f m, Median = %.3f m, RMS offset = %.3f m',mn_yoff,med_yoff,rms_yoff);
title(headline)
caxis([-0.8 0.8])
set(gca,'fontsize',fsize)
xticks(d_xtick)
datetick('x',dateformat,'keepticks')
line([datetick1 datetick2],[2*sdy+med_yoff 2*sdy+med_yoff],'Color','k')
line([datetick1 datetick2],[-2*sdy+med_yoff -2*sdy+med_yoff],'Color','k')
text(d_xtick(round(length(d_xtick)*0.5)),5*sdy+med_yoff,sprintf('%.4f m',4*sdy),'fontsize',fsize);

subplot(3,1,3)
scatter(rtx_t_all([indstart1:indend1 indstart2:indend2]),offset_rtx_posmv_sn3(3,[indstart1:indend1 indstart2:indend2]),15,posmv_zvel([indstart1:indend1 indstart2:indend2]),'filled')
grid on
datetick('x')
ylim([-1 1])
yticks(-1:0.2:1)
ylabel('Vertical offset (m)')
cbar = colorbar;
cbar.Label.String = 'Z Vel (m/s)';
colormap jet
headline = sprintf('Mean = %.3f m,Median = %.3f m, RMS offset = %.3f m',mn_zoff,med_zoff,rms_zoff);
title(headline)
caxis([-0.8 0.8])
set(gca,'fontsize',fsize)
xticks(d_xtick)
datetick('x',dateformat,'keepticks')
line([datetick1 datetick2],[2*sdz+med_zoff 2*sdz+med_zoff],'Color','k')
line([datetick1 datetick2],[-2*sdz+med_zoff -2*sdz+med_zoff],'Color','k')
text(d_xtick(round(length(d_xtick)*0.5)),8*sdz+med_zoff,sprintf('%.4f m',4*sdz),'fontsize',fsize);
%% Plot
mn_x = mean(posmvbins2_ref_gn_s3(1,[indstart1:indend1 indstart2:indend2]));
mn_y = mean(posmvbins2_ref_gn_s3(2,[indstart1:indend1 indstart2:indend2]));
mn_z = mean(posmvbins2_ref_gn_s3(3,[indstart1:indend1 indstart2:indend2]));

med_x = median(posmvbins2_ref_gn_s3(1,[indstart1:indend1 indstart2:indend2]));
med_y = median(posmvbins2_ref_gn_s3(2,[indstart1:indend1 indstart2:indend2]));
med_z = median(posmvbins2_ref_gn_s3(3,[indstart1:indend1 indstart2:indend2]));
rmsx =rms( posmvbins2_ref_gn_s3(1,[indstart1:indend1 indstart2:indend2]) - med_x);
rmsy =rms( posmvbins2_ref_gn_s3(2,[indstart1:indend1 indstart2:indend2]) - med_y);
rmsz =rms( posmvbins2_ref_gn_s3(3,[indstart1:indend1 indstart2:indend2]) - med_z);
fsize = 11;
figure(2)
clf
subplot(4,1,1)
scatter(posmv_t_all([indstart1:indend1 indstart2:indend2]),posmvbins2_ref_gn_s3(1,[indstart1:indend1 indstart2:indend2]),15,posmv_xvel([indstart1:indend1 indstart2:indend2]),'filled')
grid on
datetick('x',dateformat)
ylim([6.48 6.488])
yticks(6.478:0.002:6.490)
% ylim([0.69 .72])
% yticks([-1:.01:1])
ylabel('Bow-Stern offset (m)')
headline = sprintf('Antenna position from the granite block (POSMV Binary) \n Mean = %.3f m, Median = %.3f m, RMS Offset = %.3f m',mn_x,med_x,rmsx);
title(headline)
cbar = colorbar;
cbar.Label.String = 'X Vel (m/s)';
colormap jet
set(gca,'fontsize',fsize)
xticks(d_xtick)
datetick('x',dateformat,'keepticks')

subplot(4,1,2)
scatter(posmv_t_all([indstart1:indend1 indstart2:indend2]),posmvbins2_ref_gn_s3(2,[indstart1:indend1 indstart2:indend2]),15,posmv_yvel([indstart1:indend1 indstart2:indend2]),'filled')
grid on
datetick('x')
% ylim([-1 1])
% yticks([18:.3:19.8])
ylim([6.458 6.462])
% yticks(19:0.01:20)
ylabel('Startboard-Portside offset (m)')
colorbar
colormap jet
cbar = colorbar;
cbar.Label.String = 'Y Vel (m/s)';
headline = sprintf('Mean = %.3f m, Median = %.4f m, RMS Offset = %.3f m',mn_y,med_y,rmsy);
title(headline)
caxis([-0.8 0.8])
set(gca,'fontsize',fsize)
xticks(d_xtick)
datetick('x',dateformat,'keepticks')

subplot(4,1,3)
scatter(posmv_t_all([indstart1:indend1 indstart2:indend2]),posmvbins2_ref_gn_s3(3,[indstart1:indend1 indstart2:indend2]),15,posmv_zvel([indstart1:indend1 indstart2:indend2]),'filled')
grid on
datetick('x')
% ylim([ -6 6])
ylim([-27.762 -27.758])
ylabel('Vertical offset (m)')
cbar = colorbar;
cbar.Label.String = 'Z Vel (m/s)';
colormap jet
headline = sprintf('Mean = %.3f m, Median = %.4f m, RMS Offset = %.3f m',mn_z,med_z,rmsz);
title(headline)
caxis([-0.8 0.8])
set(gca,'fontsize',fsize)
xticks(d_xtick)
datetick('x',dateformat,'keepticks')

subplot(4,1,4)
displacement = sqrt(posmvbins2_ref_gn_s3(1,[indstart1:indend1 indstart2:indend2]).^2+posmvbins2_ref_gn_s3(2,[indstart1:indend1 indstart2:indend2]).^2+posmvbins2_ref_gn_s3(3,[indstart1:indend1 indstart2:indend2]).^2);
scatter(posmv_t_all([indstart1:indend1 indstart2:indend2]),displacement,15,'filled')
grid on
datetick('x',dateformat)
%  ylim([19.66 19.68])
% yticks([19.66:.005:19.68])
ylabel('Absolute Distance(m)')
headline= sprintf('Median = %.4f RMS Offset = %.4f m',median(displacement),rms(displacement-median(displacement)));
title(headline)
set(gca,'fontsize',fsize)
xticks(d_xtick)
datetick('x',dateformat,'keepticks')
%%
plot(posmv_t_all,sqrt(rtx_X_dis.^2+rtx_Y_dis.^2+gga_Z_dis.^2));
datetick('x')
grid on
axis tight
title(sprintf('Additional Dispalcement of the POS MV GGA: RMS = %.4f',rms(sqrt(gga_X_dis.^2+gga_Y_dis.^2+gga_Z_dis.^2))))
ylabel('m')
%% at deck
%{
n1 = 1;
n2 = 29000;

n3 = 32000;
n4 = length(posmv_t_all);

mean_heading =(posmv_heading([n1:n2 n3:n4]));
star_ax = mean_heading+90-270;
td_x = posmv_td_all(1,[n1:n2 n3:n4])-rtx_x_ant_gb([n1:n2 n3:n4])';
td_y = posmv_td_all(2,[n1:n2 n3:n4])-rtx_y_ant_gb([n1:n2 n3:n4])';
td_z = -posmv_td_all(3,[n1:n2 n3:n4]) + rtx_z_ant_gb([n1:n2 n3:n4])';
td_ab = sqrt(td_x.^2+td_y.^2);
theta = atan(abs(td_y./td_x))/pi*180;
diff_ang = theta-star_ax';
td_x_ship = -td_ab.*sin(diff_ang/180*pi);
td_y_ship = td_ab.*cos(diff_ang/180*pi);

med_x = median(td_x_ship);
med_y = median(td_y_ship);
med_z = median(td_z);
rmsx = rms(td_x_ship - med_x);
rmsy = rms(td_y_ship - med_y);
rmsz = rms(td_z - med_z);

figure(8)
clf
subplot(3,1,1)
plot(posmv_t_all([n1:n2 n3:n4]),td_x_ship)
grid on
datetick('x')
axis tight
headline = sprintf('Median = %.4f m, rms offset from median = %.4f',med_x,rmsx);
title({'Transducer offset from the RTX antenna',headline})
hold on
line([posmv_t_all(n1) posmv_t_all(n4) ],[-5.87 -5.87],'Color','r')
ylabel('Bow-Stern (m)')

subplot(3,1,2)
plot(posmv_t_all([n1:n2 n3:n4]),td_y_ship)
datetick('x')
grid on
axis tight
headline = sprintf('Median = %.4f m, rms offset from median = %.4f',med_y,rmsy);
title(headline)
hold on
line([posmv_t_all(n1) posmv_t_all(n4) ],[13.1 13.1],'Color','r')
ylabel('STarboard-Portside (m)')

subplot(3,1,3)
plot(posmv_t_all([n1:n2 n3:n4]),td_z)
datetick('x')
grid on
axis tight
headline = sprintf('Median = %.4f m, rms offset from median = %.4f',med_z,rmsz);
title(headline)
hold on
line([posmv_t_all(n1) posmv_t_all(n4)] ,[28.52 28.52],'Color','r')
ylabel('Vertical (m)')
%}
%%  Interpolation analysis
figure(7)
subplot(2,1,1)
scatter(rtx_t_all([indstart1:indend1 indstart2:indend2]),offset_rtx_posmv_sn3(1,[indstart1:indend1 indstart2:indend2]),15,gga_X_dis([indstart1:indend1 indstart2:indend2]),'filled')
grid on
datetick('x',dateformat)
ylim([-.4 .4])
% ylim([-6 -5])
% yticks([-16:.2:6])
yticks([-4:.1:5])
ylabel('Bow-Stern offset (m)')
% headline = sprintf('Transducer offset from the RTX antenna \n Mean = %.3f m,Median = %.3f m, RMS offset = %.3f m',mn_xoff,med_xoff,rms_xoff);
headline = sprintf('Antenna position difference at-sea (POSMV GGA - POSMV Binary Sensor 1 data)\n with interpolation \nMean = %.3f m, Median = %.3f m RMS offset = %.3f m',mn_xoff,med_xoff,rms_xoff);
title(headline)
cbar = colorbar;
cbar.Label.String = 'X Displacement (m)';
caxis([-0.3 0.3])
colormap jet

subplot(2,1,2)
scatter(rtx_t_all([indstart1:indend1 indstart2:indend2]),offset_rtx_posmv_sn3(1,[indstart1:indend1 indstart2:indend2]),15,gga_toffset([indstart1:indend1 indstart2:indend2]),'filled')
grid on
datetick('x',dateformat)
ylim([-.4 .4])
% ylim([-6 -5])
% yticks([-16:.2:6])
yticks([-4:.1:5])
ylabel('Bow-Stern offset (m)')
% headline = sprintf('Transducer offset from the RTX antenna \n Mean = %.3f m,Median = %.3f m, RMS offset = %.3f m',mn_xoff,med_xoff,rms_xoff);
title('Eastward')
cbar = colorbar;
cbar.Label.String = 'Time Offset (sec)';
caxis([-0.05 0.05])
colormap jet
