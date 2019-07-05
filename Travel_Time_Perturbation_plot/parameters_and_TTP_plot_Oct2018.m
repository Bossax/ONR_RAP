% Plot data points of the entire dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
%ACO LAT/LON
HEM_lat=22.738772;                  % Aug 2018
HEM_lon=-158.006186;                % Aug 2018

% icListen LAT/LON
icListen_lat = 22.739153;
icListen_lon = -158.0061254;

day = 27:30 ;               %  Edit
start_hour = 3;             % Edit
end_hour = 14;              % EDIT
hydrophone = "icListen";    % EDIT
% extract tx rx H
[tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,x_err,y_err,z_err,act_arrival,est_arrival,SNR]  = tx_rx_extraction_Oct(day,start_hour,end_hour,hydrophone);


%% Calculate Azimuth and Heading relative to the ACO
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



% -90 90
%{
% unit vector of the azimuthal vector anf the ship heading vector
% v_azmth = [sin(azmth'/180*pi);cos(azmth'/180*pi)];
% v_heading = [sin(tx_heading/180*pi);cos(tx_heading/180*pi)];
% theta = ones(length(v_azmth),1);
% for i = 1:length(v_azmth)
%     dot_pro =dot(v_azmth(:,i),v_heading(:,i));
%     theta_no_s = acos(dot_pro);
%     R = [cos(theta_no_s) -sin(theta_no_s);sin(theta_no_s) cos(theta_no_s)] ;
%     test_v_heading = R* v_heading(:,i);
%     %test_v_azmth = R'* v_azmth(:,i);
%     % assume cw rotation positive
%     if norm((test_v_heading - v_azmth(:,i))) < (1^-10)
%         theta(i) = theta_no_s/pi*180;
%     else
%         theta(i) = -theta_no_s/pi*180;
%     end
%     
%     % 90 -> 0 , -90 -> 0
%         theta(i) = -abs(theta(i))+90;
%     
% end
%}
%% Calculate parameters
ttp = (act_arrival-real(est_arrival))*1000*3600*24;

%% funtional form data fitting
r = range';
amp = (1+cos(theta/180*pi))/2;
f_sim = 10/25*r;
f_1 = (r*15/25+5).*amp + (r*10/25).*(1-amp);
f_2 = 10/25*r+(r/50).*abs(cos(theta/180*pi)-pi/2);
f_3 = (10/25)*r+r/3.*cos(theta/180*pi);
% function from Lasso Regression
f_lasso = 0.4585*r+0.2620*tx_xvel'.*cos(theta/180*pi)+0.0838*r.*cos(theta/180*pi);
f_lasso2 = 0.4855*r+0.3948*tx_xvel'.*cos(theta/180*pi)+0.0596*r.*cos(theta/180*pi);
f_lasso3 = 0.4867*r+0.0614*r.*cos(theta/180*pi)+0.3992*tx_xvel'.*cos(theta/180*pi);
f_lasso_wcable_delay = 0.4619*r;
f_lasso_wcable_delay_SA = 0.5666*r;
f_lasso_wcable_delay2 = 0.4458*r;

f_lasso_wcable_delay_nocv = 0.1978*r;
f_lasso_ic_ncable_delay = 0.4263*r+0.3486*tx_xvel'.*cos(theta/180*pi);
f_lasso_ic = 0.4276*r;
f_lasso_ic2 = 0.3961*r;
f_lasso_ic3 = 0.4082*r;
f_lasso_ic4 = -0.0616*r;
%% scope range
rmin = 0;
rmax = 30;
indl = (find(range >= rmin));
indh = (find(range <= rmax));
keep_ind = intersect(indl,indh);

range_p = range(keep_ind);
ttp_p = ttp(keep_ind);
theta_p = theta(keep_ind);
rel_v_p =rel_v(keep_ind);
f_lasso_p = f_lasso(keep_ind);
f_lasso2_p = f_lasso2(keep_ind);
f_lasso3_p = f_lasso3(keep_ind);
f_lasso_wcable_delay_p = f_lasso_wcable_delay(keep_ind);
f_lasso_wcable_delay_SA_p = f_lasso_wcable_delay_SA(keep_ind);
f_lasso_wcable_delay2_p = f_lasso_wcable_delay2(keep_ind);

f_lasso_wcable_delay_nocv_p = f_lasso_wcable_delay_nocv(keep_ind);
f_lasso_ic_ncable_delay_p = f_lasso_ic_ncable_delay(keep_ind);
f_lasso_ic_p = f_lasso_ic(keep_ind);
f_lasso_ic2_p = f_lasso_ic2(keep_ind);
f_lasso_ic3_p = f_lasso_ic3(keep_ind);
f_lasso_ic4_p = f_lasso_ic4(keep_ind);
% EDIT

function2plot =f_lasso_ic4_p;

%{
%% With respect to rel v
% figure(1)
% clf
% dot = find(vmarkers == 1);
% star = find(vmarkers == 2);
% % scatter(theta,range,40,ttp,'s')
% scatter(theta(dot),range(dot),2.8.^(abs(rel_v(dot))/1.1),ttp(dot),'*')
% hold on
% scatter(theta(star),range(star),2.8.^(abs(rel_v(star))/1.1),ttp(star),'*')
% grid on
% xlabel('Velocity relative to the ACO (Radially) (m/s)')
% ylabel('Range (km)')
% colormap jet
% cbar = colorbar;
% cbar.Label.String = 'Travel Time Perturbation (ms)';
% cbar.FontSize = 12;
% set(gca,'FontSize',14)
% label = sprintf('Squre = right-hand\nStar = left-hand');
% % yticks(0:2:26)
% % text(4,26,label,'Fontsize',12)
%  xticks(0:30:360)
%  xlim([0 360])
% 
% title('October 2018 Cruise (icListen)')
%% Polar plot
% with rel v embedded
% f2 = figure(2);
% f2.Units = 'normalized';
% f2.Position = [0.01 0.5 0.6 0.7 ];
% clf
% axp = polaraxes(f2);
% polarscatter(theta_p/180*pi,range_p,7.^(abs(rel_v_p)/2.5)+1.4,ttp_p,'filled') %
% 
% axp.ThetaZeroLocation = 'top';
% axp.ThetaDir = 'clockwise';
% axp.RLim = [0 25.5];
% axp.RTick = 0:5:30;
% axp.RTickLabelRotation = 45;
% cbar = colorbar;
% cbar.Label.String = 'ms';
% colormap jet;
% title('Travel Time Perturbation: icListen (marker size = relative v)')
%}
%% Polar Plot
% with out rel v
f3 = figure(3);
f3.Units = 'normalized';
f3.Position = [0.0 0.5 0.5 0.6];
clf
axp = polaraxes(f3);
polarscatter(theta_p/180*pi,range_p,15,ttp_p,'filled') %
axp.ThetaZeroLocation = 'top';
axp.ThetaDir = 'clockwise';
axp.RLim = [0 30];
 axp.RTick = 0:5:30;
axp.RTickLabelRotation = 45;
cbar = colorbar;
cbar.Label.String = 'ms';
colormap jet;
title('Travel Time Perturbation: icListen')
set(gca,'fontsize',18)
% caxis([-4 16]);


% Lasso regression result
%{
f7 = figure(7);
f7.Units = 'normalized';
f7.Position = [0.5 0.5 0.5 0.6];
clf
axp = polaraxes(f7);
polarscatter(theta_p/180*pi,range_p,15,function2plot,'filled') %
axp.ThetaZeroLocation = 'top';
axp.ThetaDir = 'clockwise';
axp.RLim = [0 30];
axp.RTick = 0:5:30;
axp.RTickLabelRotation = 45;
% caxis([-5 16]);
cbar = colorbar;
cbar.Label.String = 'ms';
colormap jet;
title('Travel Time Perturbation: icListen (Lasso Regression)')
%}
%% Residual
%{
f8 = figure(8);
f8.Units = 'normalized';
f8.Position = [0.5 0.1 0.5 0.6];
clf
axp = polaraxes(f8);
axp.RLim = [0 30];
axp.RTick = 0:5:30;
axp.RTickLabelRotation = 45;
cbar = colorbar;
cbar.Label.String = 'ms';
colormap jet;
caxis([-8 5]);

f4 = figure(6);
f4.Units = 'normalized';
f4.Position = [0.0 0.1 0.5 0.6];
% truncate data
for ii = length(ttp_p)
    figure(8);
    drawnow
    axp.ThetaZeroLocation = 'top';
    axp.ThetaDir = 'clockwise';
    drawnow
    polarscatter(theta_p(1:ii)/180*pi,range_p(1:ii),15,ttp_p(1:ii)'-function2plot(1:ii),'filled') 
    axp.ThetaZeroLocation = 'top';
    axp.ThetaDir = 'clockwise';
    cbar = colorbar;
    cbar.Label.String = 'ms';
    colormap jet;
    caxis([-8 5]);
    
    MSE = (sum((ttp_p(1:ii)'-function2plot(1:ii)).^2))/length(ttp);
    title(sprintf('Residual (Measurement - Funtion): MSE = %.2f',MSE))
    
    
    figure(6);
    scatter(tx_lon(1:ii),tx_lat(1:ii),20,ttp_p(1:ii)'-function2plot(1:ii),'filled')
    c = colorbar;
    c.Label.String =  'TTP (ms)';
    colormap jet
    grid on
    caxis([-8 5]);
    title(' Travel Time Perturbation Residual (Measurement - Function): Map' )

    pause(1)
    
end

%}
%% TTP vs range Comparison function and real data
figure(9)
clf
set(gcf,'Units','normalized','Position',[0.4 0.4 0.5 0.5])
scatter(range,ttp,10,rel_v,'filled')
hold on
% scatter(range_p,function2plot,10,'k','filled')
colormap jet 
cbar = colorbar;
% ylim([-5 5])
% caxis([0 360])
% cbar.Label.String = 'Azimuth (Degrees)';
cbar.Label.String = 'Relative Velocity (m/s)';
xlabel('Range (km)')
ylabel('TTP (ms)')
title('HEM Hydrophone: Travel time perturbation vs Range')
grid on
% yticks(-5:5)
% legend('Actual Data','Function','location','northwest')
%xlim([4 6])
% ylim([-6 6])
set(gca,'fontsize',15)

%% 4D
%{
f4 = figure(4);
f4.Units = 'normalized';
f4.Position = [0.0 0.5 0.5 0.6];
scatter3(rel_v,theta,range,20,ttp,'filled');
cbar = colorbar;
cbar.Label.String = 'ms';
xlabel('Relative Velocity (m/s)')
ylabel('Heading(degrees)')
zlabel('range (km)')

colormap jet;
%% TX Map
f6 = figure(6);
f6.Units = 'normalized';
f6.Position = [0.5 0.1 0.5 0.6];
scatter(icListen_lon,icListen_lat,200,'kp','filled')
grid on
colorbar;
caxis([0 8]);
colormap jet
for i = 1:length(tx_lon)
    drawnow
    hold on
    scatter(tx_lon(i),tx_lat(i),[],ttp(i),'filled')
    
end
%}
%% Map
f4 = figure(6);
f4.Units = 'normalized';
f4.Position = [0.0 0.5 0.5 0.6];
scatter(tx_lon,tx_lat,20,ttp_p,'filled')%'-function2plot
c  = colorbar;
c.Label.String =  'TTP (ms)';
colormap jet
grid on
caxis([-5 5])
% caxis([-4 16])
xlabel('Longitutde')
ylabel('Latitude')
axis tight
title(' Travel Time Perturbation Map (icListen)')
set(gca,'fontsize',16)
% title(' Travel Time Perturbation Residual (Measurement - Function): Map' )
%% Histogram
figure(7777)
histogram(ttp_p,30); %'-function2plot
xlabel('travel time perturbation (ms)')
% title('Histrogram: resdiual (Measurement - Function)')
grid on
grid minor
xticks(-8:8)
xlim([-6 6])
residual =ttp_p'-function2plot;
RMS = rms(ttp_p - mean(ttp_p))
med = median(ttp_p);
set(gca,'fontsize',13)
% title(sprintf('Histrogram: TTP Resdiual (Measurement - Function)\n Median = %.2f RMS = %.2f ms',med,RMS))
title(sprintf('Histrogram (icListen): Travel Time Perturbation \n Median = %.2f RMS = %.2f ms',med,RMS))
%% for PSD ob
% figure(10)
% scatter(tx_t,theta_p,[],range_p,'.')
% datetick('x','HH:MM')
% colormap jet
% grid on
% colorbar
% yticks([0:30:360])