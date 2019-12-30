clc% Sound Speed Profile Test
clear
close all

ACO_lat = 22.738772;
ACO_lon = -158.006186;
ACO_depth = -4729.92;                   % original depth MSL local u/e/n frame
R_e = 6371000;          % Earth Radius
%% 1. CTD data (MSL)
 fid=fopen('h304a0215.ctd');                       % 12 October 2018 cast2
D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f','headerlines',6));
fclose(fid);

pres = D(:,1);    % pressure (dbars)
temp = D(:,2);        % temperature
sal = D(:,3);         % salinity (Sp)
sal = gsw_SA_from_SP(sal,pres,ACO_lon,ACO_lat);

SS1 = gsw_sound_speed(sal,temp,pres);
z1 = -(gsw_z_from_p(pres,ACO_lat) );   

fid=fopen('h304a0202.ctd');                       % 12 October 2018 cast2


D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f','headerlines',6));
fclose(fid);

pres = D(:,1);    % pressure (dbars)
temp = D(:,2);        % temperature
sal = D(:,3);         % salinity (Sp)
sal = gsw_SA_from_SP(sal,pres,ACO_lon,ACO_lat);

gsw_infunnel(sal,temp,pres)


SS2 = gsw_sound_speed(sal,temp,pres);
z2 = -(gsw_z_from_p(pres,ACO_lat) );   
%% plot
%{
figure(1)
depth = 4000;
ax1 = axes();
line(sal,z1,'Color','b')
ax1.XLabel.String = 'Salinity';
ax1.XLabel.FontSize= 14;
ax1.YDir = 'reverse';
ax1.XColor = 'b';
ylim([0 depth])

ax2 = axes('Color', 'none','XAxisLocation','top','YColor','none','Position',ax1.Position);
line(temp,z1,'Color','r')
ax2.XLabel.String = 'Temperature';
ax2.XLabel.FontSize= 14;
ax2.YDir = 'reverse';
ax2.XColor = 'r';   
ylim([0 depth])

ax3 = axes('Color','none','XColor','none','YColor','none','Position',ax1.Position);
line(SS1,z1,'Color','k','LineStyle','--')
ax3.YDir = 'reverse';
ylim([0 depth])
%}
%%
figure(2)
plot(SS1,z1)
hold on
plot(SS2,z2)
set(gca,'YDir','reverse')
xlabel('Sound Speed (m/s)')
ylabel('Depth (m)')
set(gca,'fontsize',14)
legend('Cast 15','Cast 2','location','southwest')
grid on

figure(3)
plot(SS1(1:end-2)-SS2,z2)

set(gca,'YDir','reverse')
xlabel('Sound Speed (m/s)')
ylabel('Depth (m)')
title('Sound Speed Difference')
set(gca,'fontsize',14)
grid on
