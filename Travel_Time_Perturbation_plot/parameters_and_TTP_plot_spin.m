% Plot data points of the entire dataset
% Spins 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
%ACO LAT/LON
HEM_lat=22.738772;                  % Aug 2018
HEM_lon=-158.006186;                % Aug 2018

hydrophone = "HEM";    
% extract tx rx information
tx_t = [];
tx_lon = [];
tx_lat = [];
tx_heading = [];
tx_altitude = [];
tx_xvel = [];
range = [];
act_arrival = [];
est_arrival = [];

% June


[tx_td,tx_lond,tx_latd,tx_headingd,tx_altituded,tx_xveld,ranged,act_arrivald,est_arrivald] = tx_rx_extraction_June(20,7,8,hydrophone);
[tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,act_arrival,est_arrival] = concat_tx(tx_td,tx_lond,tx_latd,tx_headingd,tx_altituded,tx_xveld,ranged,act_arrivald,est_arrivald,tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,act_arrival,est_arrival);

[tx_td,tx_lond,tx_latd,tx_headingd,tx_altituded,tx_xveld,ranged,act_arrivald,est_arrivald] = tx_rx_extraction_June(21,12,14,hydrophone);
[tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,act_arrival,est_arrival] = concat_tx(tx_td,tx_lond,tx_latd,tx_headingd,tx_altituded,tx_xveld,ranged,act_arrivald,est_arrivald,tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,act_arrival,est_arrival);

[tx_td,tx_lond,tx_latd,tx_headingd,tx_altituded,tx_xveld,ranged,act_arrivald,est_arrivald] = tx_rx_extraction_June(21,20,22,hydrophone);
[tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,act_arrival,est_arrival] = concat_tx(tx_td,tx_lond,tx_latd,tx_headingd,tx_altituded,tx_xveld,ranged,act_arrivald,est_arrivald,tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,act_arrival,est_arrival);

% Oct
[tx_td,tx_lond,tx_latd,tx_headingd,tx_altituded,tx_xveld,ranged,x_err,y_err,z_err,act_arrivald,est_arrivald,SNR] = tx_rx_extraction_Oct(28,13,15,hydrophone);
[tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,act_arrival,est_arrival] = concat_tx(tx_td,tx_lond,tx_latd,tx_headingd,tx_altituded,tx_xveld,ranged,act_arrivald,est_arrivald,tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,act_arrival,est_arrival);
% 
[tx_td,tx_lond,tx_latd,tx_headingd,tx_altituded,tx_xveld,ranged,x_err,y_err,z_err,act_arrivald,est_arrivald,SNR] = tx_rx_extraction_Oct(29,17,19,hydrophone);
[tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,act_arrival,est_arrival] = concat_tx(tx_td,tx_lond,tx_latd,tx_headingd,tx_altituded,tx_xveld,ranged,act_arrivald,est_arrivald,tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,act_arrival,est_arrival);
% 
% [tx_td,tx_lond,tx_latd,tx_headingd,tx_altituded,tx_xveld,ranged,act_arrivald,est_arrivald] = tx_rx_extraction_Oct(30,9,12,hydrophone);
% [tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,act_arrival,est_arrival] = concat_tx(tx_td,tx_lond,tx_latd,tx_headingd,tx_altituded,tx_xveld,ranged,act_arrivald,est_arrivald,tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,act_arrival,est_arrival);
% 
% [tx_td,tx_lond,tx_latd,tx_headingd,tx_altituded,tx_xveld,ranged,act_arrivald,est_arrivald] = tx_rx_extraction_Oct(30,12,14,hydrophone);
% [tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,act_arrival,est_arrival] = concat_tx(tx_td,tx_lond,tx_latd,tx_headingd,tx_altituded,tx_xveld,ranged,act_arrivald,est_arrivald,tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,act_arrival,est_arrival);
% %}
%% window time
frame1_start= "20180620 07:40";
frame1_end= "20180620 08:40";
frame2_start= "20180621 12:30";
frame2_end= "20180621 14:40";
frame3_start= "20180621 20:28";
frame3_end= "20180621 22:28";
frame4_start= "20181028 13:25";
frame4_end= "20181028 15:40";
frame5_start= "20181029 17:00";
frame5_end= "20181029 19:15";
frame6_start= "20181030 09:40";
frame6_end= "20181030 11:10";


frame1_start = datenum(frame1_start,'yyyymmdd HH:MM');
frame1_end = datenum(frame1_end,'yyyymmdd HH:MM');
frame2_start = datenum(frame2_start,'yyyymmdd HH:MM');
frame2_end = datenum(frame2_end,'yyyymmdd HH:MM');
frame3_start = datenum(frame3_start,'yyyymmdd HH:MM');
frame3_end = datenum(frame3_end,'yyyymmdd HH:MM');
frame4_start = datenum(frame4_start,'yyyymmdd HH:MM');
frame4_end = datenum(frame4_end,'yyyymmdd HH:MM');
frame5_start = datenum(frame5_start,'yyyymmdd HH:MM');
frame5_end = datenum(frame5_end,'yyyymmdd HH:MM');
frame6_start = datenum(frame6_start,'yyyymmdd HH:MM');
frame6_end = datenum(frame6_end,'yyyymmdd HH:MM');



rm_ind = [];
ind1 = find(tx_t <= frame1_start);
ind2 = find(tx_t >= frame1_end);
ind3 = find(tx_t <= frame2_start);
ind4 = find(tx_t >= frame2_end);
ind5 = find(tx_t <= frame3_start);
ind6 = find(tx_t >= frame3_end);
ind7 = find(tx_t <= frame4_start);
ind8 = find(tx_t >= frame4_end);
ind9 = find(tx_t <= frame5_start);
ind10 = find(tx_t >= frame5_end);
ind11 = find(tx_t <= frame6_start);
ind12 = find(tx_t >= frame6_end);

C1 = intersect(ind2,ind3);
C2 = intersect(ind4,ind5);
C3 = intersect(ind6,ind7);
C4 = intersect(ind8,ind9);
C5 = intersect(ind10,ind11);
rm_ind = union(rm_ind,C1);
rm_ind = union(rm_ind,C2);
rm_ind = union(rm_ind,C3);
rm_ind = union(rm_ind,C4);
rm_ind = union(rm_ind,C5);
rm_ind = union(rm_ind,ind1);
rm_ind = union(rm_ind,ind12);
%%
tx_t(rm_ind) = [];
tx_lon(rm_ind)  = [];
tx_lat(rm_ind)  = [];
tx_heading(rm_ind) = [];
tx_altitude(rm_ind) = [];
tx_xvel(rm_ind) = [];
range(rm_ind) = [];
act_arrival(rm_ind) = [];
est_arrival(rm_ind) = [];




%% Calculate Azimuth and Heading relative to the ACO

azmth = ones(length(tx_lat),1);
for i=1:length(tx_lat)
    azmth(i) = azimuth(HEM_lat,HEM_lon,tx_lat(i),tx_lon(i));
end
theta = tx_heading;
theta = tx_heading - azmth';

for i = 1:length(theta)
   if theta(i) < 0 
       theta(i) = theta(i)+360;
   end
    
end

% relative radial velocity t.o the ACO 
% with markers representing orientations
vmarkers = zeros(1,length(theta));  % 1 = dot (right-hand) 2 = * (left-hand)
rel_v =  zeros(1,length(theta));  
for i =1:length(theta)
        if (0 <= theta(i))& (theta(i) <90)
            vmarkers(i) = 1;
            rel_v(i) = tx_xvel(i)*cos(theta(i)/180*pi); 
        elseif (90 <= theta(i))& (theta(i) < 180)
            vmarkers(i) = 2;
            rel_v(i) = tx_xvel(i)*cos(theta(i)/180*pi);         
        elseif (180 <= theta(i))& (theta(i) <270)
            vmarkers(i) = 1;
            rel_v(i) = tx_xvel(i)*cos(theta(i)/180*pi); 
         elseif (270 <= theta(i))& (theta(i) <= 360)
             vmarkers(i) = 2;
            rel_v(i) = tx_xvel(i)*cos(theta(i)/180*pi); 
        end
end


%% Calculate parameters
ttp = (act_arrival-est_arrival)*1000*3600*24;

%{
% unit vector of the azimuthal vector anf the ship heading vector
v_azmth = [sin(azmth'/180*pi);cos(azmth'/180*pi)];
v_heading = [sin(tx_heading/180*pi);cos(tx_heading/180*pi)];
theta = ones(length(v_azmth),1);
for i = 1:length(v_azmth)
    dot_pro =dot(v_azmth(:,i),v_heading(:,i));
    theta_no_s = acos(dot_pro);
    R = [cos(theta_no_s) -sin(theta_no_s);sin(theta_no_s) cos(theta_no_s)] ;
    test_v_heading = R* v_heading(:,i);
    %test_v_azmth = R'* v_azmth(:,i);
    % assume cw rotation positive
    if norm((test_v_heading - v_azmth(:,i))) < (1^-10)
        theta(i) = theta_no_s/pi*180;
    else
        theta(i) = -theta_no_s/pi*180;
    end
    
    % 90 -> 0 , -90 -> 0
        theta(i) = -abs(theta(i))+90;
    
end
%}

%%
% Plot
% Plot
figure(1)
clf
dot = find(vmarkers == 1);
star = find(vmarkers == 2);
scatter(rel_v(dot),range(dot),40,ttp(dot),'s')
hold on
scatter(rel_v(star),range(star),40,ttp(star),'*')
grid on
xlabel('Velocity relative to the ACO (Radially) (m/s)')
ylabel('Range (km)')
colormap jet
cbar = colorbar;
cbar.Label.String = 'Travel Time Perturbation (ms)';
cbar.FontSize = 12;
xlim([-1 1])
set(gca,'FontSize',14)
title('Spins')
label = sprintf('Squre = right-hand\nStar = left-hand');
yticks(5:20)
text(0.62,13,label,'Fontsize',12)
%% funtional form data fitting
r = range';
amp = (1+cos(theta/180*pi))/2;
f_sim = 10/25*r;
f_1 = (r*15/25+5).*amp + (r*10/25).*(1-amp);
f_2 = 10/25*r+(r/50).*abs(cos(theta/180*pi)-pi/2);
f_3 = (10/25)*r+r/3.*cos(theta);
%% Polar plot
f2 = figure(2);
f2.Units = 'normalized';
f2.Position = [0.4 0.5 0.6 0.7 ];
clf
axp = polaraxes(f2);
polarscatter(theta/180*pi,range,20,ttp,'filled')

axp.ThetaZeroLocation = 'top';
axp.ThetaDir = 'clockwise';
axp.RLim = [0 30];
axp.RTick = 0:5:30;
axp.RTickLabelRotation = 45;
cbar = colorbar;
cbar.Label.String = 'ms';
colormap jet;
title('Travel Time Perturbations of Spin Courses')
%% TTP vs range Conparison function and real data
figure(5)
clf
scatter(range,ttp,[],theta,'filled')
hold on
%scatter(range,f_3,[],'r','filled')
colormap jet
cbar = colorbar;
cbar.Label.String = 'Heading Rel ACO';
%cbdate('dd HH')
xlabel(' Range (km)')
ylabel('ttp (ms)')
title('Spin Course ')
grid on

%% 4D
f4 = figure(4);
f4.Units = 'normalized';
f4.Position = [0.0 0.5 0.5 0.6];
scatter3(tx_xvel,theta,range,20,ttp,'filled');
cbar = colorbar;
cbar.Label.String = 'ms';
xlabel('Ship Speed (m/s)')
ylabel('Heading(degrees)')
zlabel('Range (km)')
ylim([0 360])
colormap jet;


%% 
function [tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,act_arrival,est_arrival] = concat_tx(tx_td,tx_lond,tx_latd,tx_headingd,tx_altituded,tx_xveld,ranged,act_arrivald,est_arrivald,tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,act_arrival,est_arrival)
tx_t = horzcat(tx_t,tx_td);
tx_lon = horzcat(tx_lon,tx_lond);
tx_lat = horzcat(tx_lat,tx_latd);
tx_heading = horzcat(tx_heading,tx_headingd);
tx_altitude = horzcat(tx_altitude,tx_altituded);
tx_xvel = horzcat(tx_xvel,tx_xveld);
range = horzcat(range,ranged);
act_arrival = horzcat(act_arrival,act_arrivald);
est_arrival = horzcat(est_arrival,est_arrivald);


end