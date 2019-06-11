%% Calculate SS field with hydrophone position
clearvars
close all
%%
load EOF_mode_1
%EOF_mode=ones(2405,1);
%EOF_mode(EOF_mode<0)=0;
load EOF_lambda_1
%EOF_lambda=1;
load RAP_test_IC_cable_delay
EOF_mode = vertcat(EOF_mode,EOF_mode(end)*ones(5,1));
% load RMS_SS
% EOF_mode=RMS_SS./(max(RMS_SS));
% EOF_lambda=max(abs(RMS_SS))^2;


% ACO_lat=22.738772;              %Updated (18/6/18)
% ACO_lon=-158.006186;            %Updated (18/6/18)

icListen_lat=22.73911569;                  % March 2019 #3
icListen_lon=-158.006106601;                % March 2019

%%%%%%%%Transmissions%%%%%%%%%%
t_act=t_act;
t_tx=t_tx;
t_diff=t_diff;
boat_lon=lon;
boat_lat=lat;
z_tx=z;
heading=tx_heading;
x_err=lon_err;
y_err=lat_err;
z_err=z_err;

%Boat distance from hydrophone
for i=1:length(boat_lat)
    x_dist(i)=dist([icListen_lat boat_lat(i)],[icListen_lon boat_lon(i)]);
end






%%%%%Initialize variables%%%%%%
num_tx=length(boat_lon);
% g_bound=28450;
g_bound=50000;
%Set grid bounds to 30km
[~,max_lat,~]=m_fdist(icListen_lon,icListen_lat,0,g_bound);
[max_lon,~,~]=m_fdist(icListen_lon,icListen_lat,90,g_bound);
max_lon=max_lon-360;
[~,min_lat,~]=m_fdist(icListen_lon,icListen_lat,180,g_bound);
[min_lon,~,~]=m_fdist(icListen_lon,icListen_lat,270,g_bound);
min_lon=min_lon-360;
grid_num = 21;
n_cols=grid_num;
n_rows=grid_num;
x_tick=min_lon:(max_lon-min_lon)/n_cols:max_lon;
y_tick=min_lat:(max_lat-min_lat)/n_rows:max_lat;
tot_grid_num = n_rows^2;
for i=1:n_cols
x_cen(i)=(x_tick(i)+x_tick(i+1))/2;
end
for i=1:n_rows
y_cen(i)=(y_tick(i)+y_tick(i+1))/2;
end
%%


G_SS=zeros(length(y_cen)*length(x_cen),num_tx);      %Observation matrix

%Boat Bearing
for i=1:length(boat_lat)
    
    [r,b1,b2]=dist([boat_lat(i) icListen_lat], [boat_lon(i) icListen_lon]);
    
    boat_bearing(i)=-b1+90;
    
end

% Transmission Map
% draw circle
R = 25000;
num_point = 2000;
inc_ang =360/num_point;
x_cir = zeros(1,num_point);
y_cir = zeros(1,num_point);
for ii = 1:num_point
    angle = inc_ang*ii;
    [x_cir(ii),y_cir(ii),~]=m_fdist(icListen_lon,icListen_lat,angle,R); 
end
x_cir = x_cir-360;

figure(1)
clf
imagesc(x_cen,y_cen,rand(length(y_cen),length(y_cen)))
colormap gray
hold on
scatter(icListen_lon,icListen_lat,200,'pr','filled')
set(gca,'YDir','normal')
plot(x_cir,y_cir,'r')

%% Estimate travel time and total distance traveled in each pixel (HOTS CTD)
total_pixel_distance=zeros(num_tx,max([n_cols;n_rows]));
total_pixel_num=zeros(num_tx,max([n_cols;n_rows]));
arc_dist_all=zeros(num_tx,max([n_cols;n_rows]));
arc_lengths = zeros(length(EOF_mode),length(boat_lat));
ray_angles = zeros(length(EOF_mode),length(boat_lat));;
M = zeros(length(boat_lat),1);
for i=1:length(boat_lat)
    i
    %Surface distance for each pixel
    [arc_length_hd,~,~,SS_z_apriori,z,SS_avg_apriori,surface_dist,~,ray_angles_hd]=ray_trace_w_curvature_v3(x_dist(i),z_tx(i));
    [tot_pix_dist,tot_pix_num]=lat_lon_dist(boat_lat(i),boat_lon(i),icListen_lat,icListen_lon,x_tick,y_tick);
    
    total_pixel_distance(i,1:length(tot_pix_dist))=tot_pix_dist;
    total_pixel_num(i,1:length(tot_pix_num))=tot_pix_num;
    M(i) = length(tot_pix_num)-1;
    
    len = length(arc_length_hd);
    arc_lengths(1:len,i) = arc_length_hd;
    len2 = length(ray_angles_hd);
    if len2 <size(ray_angles,1)
        add_len = size(ray_angles,1) - len2;
        add_ones = ones(1,add_len);
        ray_angles(:,i) = horzcat(ray_angles_hd,ray_angles_hd(end).*add_ones);
    else
        ray_angles(:,i) = ray_angles_hd;
    end
    
    %Total Arc Lengths for each pixel
    start_pos=1;
    arc_depth(1)=1;
    for ii=1:length(tot_pix_dist)-1
        A=find(cumsum(surface_dist(start_pos:end))>=tot_pix_dist(ii));
        arc_dist(ii)=sum(arc_lengths(start_pos:A(1)+start_pos-1,i));
        arc_depth(ii+1)=A(1)+start_pos-1;
        start_pos=start_pos+A(1)-1;
    end
    arc_dist(length(tot_pix_dist))=sum(arc_lengths(start_pos:end,i));
    arc_depth(length(tot_pix_dist)+1)=length(SS_z_apriori);
    
    arc_deptharc_dist_all(i,1:length(arc_dist))=arc_dist;
    
    
    %G(tot_pix_num,i)=arc_dist;
    for ii=1:length(arc_dist)
        if arc_depth(ii+1)>arc_depth(ii)
            G_SS(tot_pix_num(ii),i)=sum(EOF_mode(arc_depth(ii):arc_depth(ii+1)-1).*(-arc_lengths(arc_depth(ii):arc_depth(ii+1)-1,i)./(((SS_z_apriori(arc_depth(ii):arc_depth(ii+1)-1))).^2)));
        end
    end
    % fprintf('%d / %d \n',i,num_tx)
    
    
    % plot
%     figure(1)
%     
%     hold on
%     plot([lon(i) icListen_lon],[lat(i) icListen_lat],'Linewidth',1,'Color','m')
%     tot_pix_dist
%     tot_pix_num;
    clear arc_dist arc_length_hd arc_depth
    
end
arc_lengths=real(arc_lengths);
arc_dist_all=real(arc_dist_all);
G_SS=real(G_SS);
ray_angles=real(ray_angles);


%Apriori SS field
SS_xy_apriori=SS_avg_apriori.*ones(n_cols,n_rows);
m_prior=SS_avg_apriori*ones(n_cols*n_rows,1);


%Uncertainty in apriori SS field (Cm) and hydrophone position (x,y,z)
%SS_uncertainty=1;                                   %Averaged SS field uncertainty
SS_uncertainty=sqrt(EOF_lambda);
pos_uncertainty=[10,10,5];                          %Hydrophone position uncertainty
%pos_uncertainty=[100,100,0.0001];
%pos_uncertainty=[0.0001,0.0001,0.0001];


%% Uncertainty in travel time (Cd)
t_uncertainty=0.0004;               %Timing uncertainty
% t_uncertainty=0.003;
x_anomaly=(((sin(ray_angles(1,:))./SS_avg_apriori).*cosd(boat_bearing)).^2).*(x_err).^2;
y_anomaly=(((sin(ray_angles(1,:))./SS_avg_apriori).*sind(boat_bearing)).^2).*(y_err).^2;
z_anomaly=((cos(ray_angles(1,:))./SS_avg_apriori).^2).*(z_err).^2;
% x_anomaly=(((sin(ray_angles(1,:))./SS_avg_apriori).*cosd(boat_bearing)).*(x_err')).^2;
% y_anomaly=(((sin(ray_angles(1,:))./SS_avg_apriori).*sind(boat_bearing)).*(y_err')).^2;
% z_anomaly=(((cos(ray_angles(1,:))./SS_avg_apriori)).*(z_err')).^2;
SS_anomaly=x_anomaly+y_anomaly+z_anomaly;
% Cd=(t_uncertainty)^2.*diag(ones(num_tx,1));
Cd=diag(SS_anomaly+(t_uncertainty)^2);

%%
%Create Model error covariance (P)
[~,P,~]=gaussian_covariance_SS_w_position(SS_uncertainty,pos_uncertainty,n_cols,n_rows,x_tick,y_tick,SS_avg_apriori);
% SS_actual_reshape=reshape(SS_xy_actual',[n_cols*n_rows,1]);
% m0_reshape=reshape(m0(1:end-3),[n_cols,n_rows]);
% m0_reshape=rot90(fliplr(m0_reshape));

%% Priori cov matrix
figure(1)
clf
set(gcf,'name','Covariance Matrix','Units','normalized','Position',[0.3 .5 0.35 0.45])
imagesc(P)
axis tight
grid off
title('Gaussian Covariance Matrix')
colorbar
caxis([0 14])
colormap jet



%% Inverse techniques to find new SS field

%MATRIX EQUATION [G][m]=[d]
%G=>observation
%m=>sound speed
%d=>travel time


G_SS=G_SS';


arc_d=sum(arc_dist_all,2);
%%
% for i=1:size(arc_lengths,2)
% G_x(i)=sum(cos(ray_angles(:,i)).*((arc_lengths(:,i)./arc_d(i))./SS_z_apriori(2:end)))*cosd(boat_bearing(i));
% G_y(i)=sum(cos(ray_angles(:,i)).*((arc_lengths(:,i)./arc_d(i))./SS_z_apriori(2:end)))*sind(boat_bearing(i));
% G_z(i)=sum(sin(ray_angles(:,i)).*((arc_lengths(:,i)./arc_d(i))./SS_z_apriori(2:end)));
% end
G_x=cosd(boat_bearing).*sin(ray_angles(end,:))./SS_z_apriori(end);
G_y=sind(boat_bearing).*sin(ray_angles(end,:))./SS_z_apriori(end);
G_z=cos(ray_angles(end,:))./SS_z_apriori(end);
G_pos=horzcat(G_x',G_y',G_z');

G=horzcat(G_SS,G_pos);
% G=G_pos;

% act_tt=(t_act-t_tx).*3600*24;
% est_tt=act_tt-t_diff;

Function = 0.4082*x_dist/1000;
t_diff_new = t_diff*1000 - Function;
d=t_diff_new'/1000;
%%
% Cm_all=P(end-2:end,end-2:end);
Cm_all = P;

% GI2=Cm*G'*pinv(G*Cm*G'+Cd);
% GI=Cm_all*G'*pinv(G*Cm_all*G'+Cd);
% G = G_SS;
GI=Cm_all*transpose(G)*inv((G*Cm_all*transpose(G)+diag(Cd)));


m_recov= GI*d;
% m2=GI2*d;

% m_SS=m(1:end-3);                    %Modeled SS field
% m_pos=m(end-2:end);                 %Hydrophone position corrections (x,y,z)

P_post=(eye(length(Cm_all))-(GI*G))*Cm_all;
model_std=sqrt(diag(P_post));

% R (Correlation Matrix)
for i=1:size(P_post,1)
    for j=1:size(P_post,2)
        
        R(i,j)=P_post(i,j)/(model_std(i)*model_std(j));
        
    end
end

Data_Res_Matrix=G*GI;
Model_Res_Matrix=GI*G;

%{
% SS_xy_post=SS_avg_apriori+m_SS;
% SS_xy_post=reshape(SS_xy_post,[n_rows n_cols])';
% SS_change=reshape(m_SS,[n_rows n_cols])';
% %SS_change2=reshape(m2,[n_rows n_cols])';
% 
% 
% 
% 
% 
% % max_lim=max([max(max(SS_xy_actual)) max(max(SS_xy_post))]);
% % min_lim=min([min(min(SS_xy_actual)) min(min(SS_xy_post))]);
% max_lim=max(max(SS_change));
% min_lim=min(min(SS_change));
%}
%% Plot Simulation
%{
% % figure(1)
% % imagesc(center_x,center_y,m0_reshape)
% % set(gca,'ydir','normal')
% % colormap jet
% % hold on
% % plot(boat_lon(137:340),boat_lat(137:340),'k')
% % % scatter(boat_lon,boat_lat,'k')
% % colorbar
% % caxis([min_lim max_lim])
% % title('Actual Sound Speed Field')
% 
% figure(20)
% imagesc(center_x,center_y,SS_change)
% set(gca,'ydir','normal')
% colormap jet
% hold on
% for i=1:length(boat_lat)
%     plot([boat_lon(i) ACO_lon],[boat_lat(i) ACO_lat],'k')
% end
% colorbar
% caxis([min_lim max_lim])
% title('Estimated Sound Speed Field')
% 
% % figure(3)
% % imagesc(center_x,center_y,m0_reshape-SS_change)
% % set(gca,'ydir','normal')
% % colormap jet
% % hold on
% % plot(boat_lon(137:340),boat_lat(137:340),'k')
% % % scatter(boat_lon,boat_lat,'k')
% % colorbar
% % title('Difference between actual and modeled SS')
% 
% figure(40)
% imagesc(center_x,center_y,reshape(model_std(1:end-3),[n_rows n_cols])')
% set(gca,'ydir','normal')
% colormap jet
% hold on
% for i=1:length(boat_lat)
%     plot([boat_lon(i) ACO_lon],[boat_lat(i) ACO_lat],'k')
% end
% % scatter(boat_lon,boat_lat,'k')
% colorbar
% title('A Posteriori Error')
% 

% %Expected Travel Time Perturbation
% pointsize=10;
% rms_tt_w_pos=diag(G*Cm_all*G');
% rms_tt_no_pos=diag(G(:,1:end-3)*Cm_all(1:end-3,1:end-3)*G(:,1:end-3)');
% figure(234)
% %scatter(x_dist./1000,sqrt(rms_tt_w_pos).*1000,'b')
% hold on
% scatter(x_dist./1000,sqrt(rms_tt_no_pos).*1000,'b')
% scatter(x_dist./1000,-sqrt(rms_tt_no_pos).*1000,'b')
% scatter(x_dist./1000,d.*1000,pointsize,linspace(1,length(x_dist),length(x_dist)))
% colormap jet
% %scatter(x_dist./1000,-sqrt(rms_tt_w_pos).*1000,'b')
% xlabel('Range (km)')
% ylabel('Perturbation (ms)')
% title('Travel Time Perturbation')
% %legend('w/ position','w/o position')
% legend('Expected')
% grid on


%%%%%% Scatter Plots of Data

% t_tx_all=t_tx;
% boat_heading=heading;
% %Lat/Lon of residual (offset)
% figure(234234)
% pointsize=10;
% hold on
% grid on
% colormap jet
% title('Time residual (w/ pos and w/o SS correction)')
% c=colorbar;
% c.Label.String='residual (ms)';
% 
% 
% %25 km Paths
% A=find(t_tx_all<datenum(2017,6,7,23,0,0));      %First
% scatter(boat_lon(A),boat_lat(A),pointsize,(d(A)-(G(A,:)*m)).*1000)
% 
% A=find(t_tx_all>datenum(2017,6,11,4,45,0) & t_tx_all<datenum(2017,6,11,13,40,0));   %Second
% clear plot_lon plot_lat quiv_lon quiv_lat
% for i=1:length(A)
%     if boat_lon(A(i))>=ACO_lon
%         plot_lon(i)=boat_lon(A(i))+0.02;
%         quiv_lon(i)=boat_lon(A(i))+0.01;
%     elseif boat_lon(A(i))<ACO_lon
%         plot_lon(i)=boat_lon(A(i))-0.02;
%         quiv_lon(i)=boat_lon(A(i))-0.01;
%     end
%     if boat_lat(A(i))>=ACO_lat
%         plot_lat(i)=boat_lat(A(i))+0.02;
%         quiv_lat(i)=boat_lat(A(i))+0.01;
%     elseif boat_lat(A(i))<ACO_lat
%         plot_lat(i)=boat_lat(A(i))-0.02;
%         quiv_lat(i)=boat_lat(A(i))-0.01;
%     end
% end
% scatter(plot_lon,plot_lat,pointsize,(d(A)-(G(A,:)*m)).*1000)
% quiver(quiv_lon(1:10:end),quiv_lat(1:10:end),sind(boat_heading(A(1:10:end))'),cosd(boat_heading(A(1:10:end))'),0.25,'k')
% 
% 
% 
% %15 km Paths
% A=find(t_tx_all>datenum(2017,6,8,16,37,0) & t_tx_all<datenum(2017,6,8,21,50,0));    %First
% scatter(boat_lon(A),boat_lat(A),pointsize,(d(A)-(G(A,:)*m)).*1000)
% 
% A=find(t_tx_all>datenum(2017,6,11,13,50,0) & t_tx_all<datenum(2017,6,11,18,50,0));  %Second
% clear plot_lon plot_lat quiv_lon quiv_lat
% for i=1:length(A)
%     if boat_lon(A(i))>=ACO_lon
%         plot_lon(i)=boat_lon(A(i))+0.02;
%         quiv_lon(i)=boat_lon(A(i))+0.01;
%     elseif boat_lon(A(i))<ACO_lon
%         plot_lon(i)=boat_lon(A(i))-0.02;
%         quiv_lon(i)=boat_lon(A(i))-0.01;
%     end
%     if boat_lat(A(i))>=ACO_lat
%         plot_lat(i)=boat_lat(A(i))+0.02;
%         quiv_lat(i)=boat_lat(A(i))+0.01;
%     elseif boat_lat(A(i))<ACO_lat
%         plot_lat(i)=boat_lat(A(i))-0.02;
%         quiv_lat(i)=boat_lat(A(i))-0.01;
%     end
% end
% scatter(plot_lon,plot_lat,pointsize,(d(A)-(G(A,:)*m)).*1000)
% quiver(quiv_lon(1:10:end),quiv_lat(1:10:end),sind(boat_heading(A(1:10:end))'),cosd(boat_heading(A(1:10:end))'),0.3,'k')
% 
% 
% %5 km Paths
% A=find(t_tx_all>datenum(2017,6,8,22,30,0) & t_tx_all<datenum(2017,6,8,23,57,0));    %First
% scatter(boat_lon(A),boat_lat(A),pointsize,(d(A)-(G(A,:)*m)).*1000)
% 
% A=find(t_tx_all>datenum(2017,6,11,19,20,0) & t_tx_all<datenum(2017,6,11,21,3,0));   %Second
% clear plot_lon plot_lat quiv_lon quiv_lat
% for i=1:length(A)
%     if boat_lon(A(i))>=ACO_lon
%         plot_lon(i)=boat_lon(A(i))+0.02;
%         quiv_lon(i)=boat_lon(A(i))+0.01;
%     elseif boat_lon(A(i))<ACO_lon
%         plot_lon(i)=boat_lon(A(i))-0.02;
%         quiv_lon(i)=boat_lon(A(i))-0.01;
%     end
%     if boat_lat(A(i))>=ACO_lat
%         plot_lat(i)=boat_lat(A(i))+0.02;
%         quiv_lat(i)=boat_lat(A(i))+0.01;
%     elseif boat_lat(A(i))<ACO_lat
%         plot_lat(i)=boat_lat(A(i))-0.02;
%         quiv_lat(i)=boat_lat(A(i))-0.01;
%     end
% end
% scatter(plot_lon,plot_lat,pointsize,(d(A)-(G(A,:)*m)).*1000)
% quiver(quiv_lon(1:10:end),quiv_lat(1:10:end),sind(boat_heading(A(1:10:end))'),cosd(boat_heading(A(1:10:end))'),0.35,'k')
% 
% 
% %Radials
% A=find(t_tx_all>datenum(2017,6,7,23,0,0) & t_tx_all<datenum(2017,6,8,10,0,0));      %First
% clear plot_lon plot_lat quiv_lon quiv_lat
% for i=1:length(A)
% if boat_heading(A(i))<=(0+45) || boat_heading(A(i))>=(0+315)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<=(180+45) && boat_heading(A(i))>=(180-45)
%     plot_lon(i)=boat_lon(A(i))-0.01;
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<(90+45) && boat_heading(A(i))>(90-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<(270+45) && boat_heading(A(i))>(270-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i))+0.01;    
% end
% end
% scatter(plot_lon,plot_lat,pointsize,(d(A)-(G(A,:)*m)).*1000)
% quiver(plot_lon(1:10:end),plot_lat(1:10:end),sind(boat_heading(A(1:10:end))'),cosd(boat_heading(A(1:10:end))'),0.35,'k')
% 
% A=find(t_tx_all>datenum(2017,6,10,12,2,0) & t_tx_all<datenum(2017,6,11,4,45,0));    %Second
% clear plot_lon plot_lat quiv_lon quiv_lat
% for i=1:length(A)
% if boat_heading(A(i))<=(0+45) || boat_heading(A(i))>=(0+315)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<=(180+45) && boat_heading(A(i))>=(180-45)
%     plot_lon(i)=boat_lon(A(i))-0.01;
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<(90+45) && boat_heading(A(i))>(90-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<(270+45) && boat_heading(A(i))>(270-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i))+0.01;    
% end
% end
% plot_lon=plot_lon+0.02;
% plot_lat=plot_lat-0.02;
% scatter(plot_lon,plot_lat,pointsize,(d(A)-(G(A,:)*m)).*1000)
% quiver(plot_lon(1:10:end),plot_lat(1:10:end),sind(boat_heading(A(1:10:end))'),cosd(boat_heading(A(1:10:end))'),0.35,'k')
% 
% 
% A=find(t_tx_all>datenum(2017,6,11,21,3,0));             %Third
% clear plot_lon plot_lat quiv_lon quiv_lat
% for i=1:length(A)
% if boat_heading(A(i))<=(0+45) || boat_heading(A(i))>=(0+315)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<=(180+45) && boat_heading(A(i))>=(180-45)
%     plot_lon(i)=boat_lon(A(i))-0.01;
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<(90+45) && boat_heading(A(i))>(90-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<(270+45) && boat_heading(A(i))>(270-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i))+0.01;    
% end
% end
% plot_lon=plot_lon-0.02;
% plot_lat=plot_lat+0.02;
% scatter(plot_lon,plot_lat,pointsize,(d(A)-(G(A,:)*m)).*1000)
% quiver(plot_lon(1:10:end),plot_lat(1:10:end),sind(boat_heading(A(1:10:end))'),cosd(boat_heading(A(1:10:end))'),0.35,'k')
% 
% 
% 
% %Stationary
% A=find(t_tx_all>datenum(2017,6,8,10,0,0) & t_tx_all<datenum(2017,6,8,15,0,0));
% scatter(boat_lon(A),boat_lat(A),pointsize,(d(A)-(G(A,:)*m)).*1000)
% A=find(t_tx_all>datenum(2017,6,9,0,10,0) & t_tx_all<datenum(2017,6,10,12,2,0));
% scatter(boat_lon(A),boat_lat(A),pointsize,(d(A)-(G(A,:)*m)).*1000)
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %Lat/Lon of actual v estimated (offset)
% figure(234235)
% pointsize=10;
% hold on
% grid on
% colormap jet
% title('Actual vs. Estimated TT (w/o pos and w/o SS correction)')
% c=colorbar;
% c.Label.String='time diff (ms)';
% 
% %25 km Paths
% A=find(t_tx_all<datenum(2017,6,7,23,0,0));      %First
% scatter(boat_lon(A),boat_lat(A),pointsize,d(A).*1000)
% 
% A=find(t_tx_all>datenum(2017,6,11,4,45,0) & t_tx_all<datenum(2017,6,11,13,40,0));   %Second
% clear plot_lon plot_lat
% for i=1:length(A)
%     if boat_lon(A(i))>=ACO_lon
%         plot_lon(i)=boat_lon(A(i))+0.02;
%     elseif boat_lon(A(i))<ACO_lon
%         plot_lon(i)=boat_lon(A(i))-0.02;
%     end
%     if boat_lat(A(i))>=ACO_lat
%         plot_lat(i)=boat_lat(A(i))+0.02;
%     elseif boat_lat(A(i))<ACO_lat
%         plot_lat(i)=boat_lat(A(i))-0.02;
%     end
% end
% scatter(plot_lon,plot_lat,pointsize,d(A).*1000)
% 
% 
% %15 km Paths
% A=find(t_tx_all>datenum(2017,6,8,16,37,0) & t_tx_all<datenum(2017,6,8,21,50,0));    %First
% scatter(boat_lon(A),boat_lat(A),pointsize,d(A).*1000)
% 
% A=find(t_tx_all>datenum(2017,6,11,13,50,0) & t_tx_all<datenum(2017,6,11,18,50,0));  %Second
% clear plot_lon plot_lat
% for i=1:length(A)
%     if boat_lon(A(i))>=ACO_lon
%         plot_lon(i)=boat_lon(A(i))+0.02;
%     elseif boat_lon(A(i))<ACO_lon
%         plot_lon(i)=boat_lon(A(i))-0.02;
%     end
%     if boat_lat(A(i))>=ACO_lat
%         plot_lat(i)=boat_lat(A(i))+0.02;
%     elseif boat_lat(A(i))<ACO_lat
%         plot_lat(i)=boat_lat(A(i))-0.02;
%     end
% end
% scatter(plot_lon,plot_lat,pointsize,d(A).*1000)
% 
% 
% %5 km Paths
% A=find(t_tx_all>datenum(2017,6,8,22,30,0) & t_tx_all<datenum(2017,6,8,23,57,0));    %First
% scatter(boat_lon(A),boat_lat(A),pointsize,d(A).*1000)
% 
% A=find(t_tx_all>datenum(2017,6,11,19,20,0) & t_tx_all<datenum(2017,6,11,21,3,0));   %Second
% clear plot_lon plot_lat
% for i=1:length(A)
%     if boat_lon(A(i))>=ACO_lon
%         plot_lon(i)=boat_lon(A(i))+0.02;
%     elseif boat_lon(A(i))<ACO_lon
%         plot_lon(i)=boat_lon(A(i))-0.02;
%     end
%     if boat_lat(A(i))>=ACO_lat
%         plot_lat(i)=boat_lat(A(i))+0.02;
%     elseif boat_lat(A(i))<ACO_lat
%         plot_lat(i)=boat_lat(A(i))-0.02;
%     end
% end
% scatter(plot_lon,plot_lat,pointsize,d(A).*1000)
% 
% 
% %Radials
% A=find(t_tx_all>datenum(2017,6,7,23,0,0) & t_tx_all<datenum(2017,6,8,10,0,0));      %First
% clear plot_lon plot_lat
% for i=1:length(A)
% if boat_heading(A(i))<=(0+45) || boat_heading(A(i))>=(0+315)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<=(180+45) && boat_heading(A(i))>=(180-45)
%     plot_lon(i)=boat_lon(A(i))-0.01;
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<(90+45) && boat_heading(A(i))>(90-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<(270+45) && boat_heading(A(i))>(270-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i))+0.01;    
% end
% end
% scatter(plot_lon,plot_lat,pointsize,d(A).*1000)
% 
% A=find(t_tx_all>datenum(2017,6,10,12,2,0) & t_tx_all<datenum(2017,6,11,4,45,0));    %Second
% clear plot_lon plot_lat
% for i=1:length(A)
% if boat_heading(A(i))<=(0+45) || boat_heading(A(i))>=(0+315)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<=(180+45) && boat_heading(A(i))>=(180-45)
%     plot_lon(i)=boat_lon(A(i))-0.01;
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<(90+45) && boat_heading(A(i))>(90-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<(270+45) && boat_heading(A(i))>(270-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i))+0.01;    
% end
% end
% plot_lon=plot_lon+0.02;
% plot_lat=plot_lat-0.02;
% scatter(plot_lon,plot_lat,pointsize,d(A).*1000)
% 
% A=find(t_tx_all>datenum(2017,6,11,21,3,0));             %Third
% clear plot_lon plot_lat
% for i=1:length(A)
% if boat_heading(A(i))<=(0+45) || boat_heading(A(i))>=(0+315)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<=(180+45) && boat_heading(A(i))>=(180-45)
%     plot_lon(i)=boat_lon(A(i))-0.01;
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<(90+45) && boat_heading(A(i))>(90-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<(270+45) && boat_heading(A(i))>(270-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i))+0.01;    
% end
% end
% plot_lon=plot_lon-0.02;
% plot_lat=plot_lat+0.02;
% scatter(plot_lon,plot_lat,pointsize,d(A).*1000)
% 
% 
% %Stationary
% A=find(t_tx_all>datenum(2017,6,8,10,0,0) & t_tx_all<datenum(2017,6,8,15,0,0));
% scatter(boat_lon(A),boat_lat(A),pointsize,d(A).*1000)
% A=find(t_tx_all>datenum(2017,6,9,0,10,0) & t_tx_all<datenum(2017,6,10,12,2,0));
% scatter(boat_lon(A),boat_lat(A),pointsize,d(A).*1000)
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %Lat/Lon of vessel heading (offset)
% figure(234236)
% pointsize=10;
% hold on
% grid on
% colormap jet
% title('Vessel Heading')
% c=colorbar;
% c.Label.String='degree';
% 
% %25 km Paths
% A=find(t_tx_all<datenum(2017,6,7,23,0,0));      %First
% scatter(boat_lon(A),boat_lat(A),pointsize,boat_heading(A))
% 
% A=find(t_tx_all>datenum(2017,6,11,4,45,0) & t_tx_all<datenum(2017,6,11,13,40,0));   %Second
% clear plot_lon plot_lat
% for i=1:length(A)
%     if boat_lon(A(i))>=ACO_lon
%         plot_lon(i)=boat_lon(A(i))+0.02;
%     elseif boat_lon(A(i))<ACO_lon
%         plot_lon(i)=boat_lon(A(i))-0.02;
%     end
%     if boat_lat(A(i))>=ACO_lat
%         plot_lat(i)=boat_lat(A(i))+0.02;
%     elseif boat_lat(A(i))<ACO_lat
%         plot_lat(i)=boat_lat(A(i))-0.02;
%     end
% end
% scatter(plot_lon,plot_lat,pointsize,boat_heading(A))
% 
% 
% %15 km Paths
% A=find(t_tx_all>datenum(2017,6,8,16,37,0) & t_tx_all<datenum(2017,6,8,21,50,0));    %First
% scatter(boat_lon(A),boat_lat(A),pointsize,boat_heading(A))
% 
% A=find(t_tx_all>datenum(2017,6,11,13,50,0) & t_tx_all<datenum(2017,6,11,18,50,0));  %Second
% clear plot_lon plot_lat
% for i=1:length(A)
%     if boat_lon(A(i))>=ACO_lon
%         plot_lon(i)=boat_lon(A(i))+0.02;
%     elseif boat_lon(A(i))<ACO_lon
%         plot_lon(i)=boat_lon(A(i))-0.02;
%     end
%     if boat_lat(A(i))>=ACO_lat
%         plot_lat(i)=boat_lat(A(i))+0.02;
%     elseif boat_lat(A(i))<ACO_lat
%         plot_lat(i)=boat_lat(A(i))-0.02;
%     end
% end
% scatter(plot_lon,plot_lat,pointsize,boat_heading(A))
% 
% 
% %5 km Paths
% A=find(t_tx_all>datenum(2017,6,8,22,30,0) & t_tx_all<datenum(2017,6,8,23,57,0));    %First
% scatter(boat_lon(A),boat_lat(A),pointsize,boat_heading(A))
% 
% A=find(t_tx_all>datenum(2017,6,11,19,20,0) & t_tx_all<datenum(2017,6,11,21,3,0));   %Second
% clear plot_lon plot_lat
% for i=1:length(A)
%     if boat_lon(A(i))>=ACO_lon
%         plot_lon(i)=boat_lon(A(i))+0.02;
%     elseif boat_lon(A(i))<ACO_lon
%         plot_lon(i)=boat_lon(A(i))-0.02;
%     end
%     if boat_lat(A(i))>=ACO_lat
%         plot_lat(i)=boat_lat(A(i))+0.02;
%     elseif boat_lat(A(i))<ACO_lat
%         plot_lat(i)=boat_lat(A(i))-0.02;
%     end
% end
% scatter(plot_lon,plot_lat,pointsize,boat_heading(A))
% 
% 
% %Radials
% A=find(t_tx_all>datenum(2017,6,7,23,0,0) & t_tx_all<datenum(2017,6,8,10,0,0));      %First
% clear plot_lon plot_lat
% for i=1:length(A)
% if boat_heading(A(i))<=(0+45) || boat_heading(A(i))>=(0+315)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<=(180+45) && boat_heading(A(i))>=(180-45)
%     plot_lon(i)=boat_lon(A(i))-0.01;
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<(90+45) && boat_heading(A(i))>(90-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<(270+45) && boat_heading(A(i))>(270-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i))+0.01;    
% end
% end
% scatter(plot_lon,plot_lat,pointsize,boat_heading(A))
% 
% A=find(t_tx_all>datenum(2017,6,10,12,2,0) & t_tx_all<datenum(2017,6,11,4,45,0));    %Second
% clear plot_lon plot_lat
% for i=1:length(A)
% if boat_heading(A(i))<=(0+45) || boat_heading(A(i))>=(0+315)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<=(180+45) && boat_heading(A(i))>=(180-45)
%     plot_lon(i)=boat_lon(A(i))-0.01;
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<(90+45) && boat_heading(A(i))>(90-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<(270+45) && boat_heading(A(i))>(270-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i))+0.01;    
% end
% end
% plot_lon=plot_lon+0.02;
% plot_lat=plot_lat-0.02;
% scatter(plot_lon,plot_lat,pointsize,boat_heading(A))
% 
% A=find(t_tx_all>datenum(2017,6,11,21,3,0));             %Third
% clear plot_lon plot_lat
% for i=1:length(A)
% if boat_heading(A(i))<=(0+45) || boat_heading(A(i))>=(0+315)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<=(180+45) && boat_heading(A(i))>=(180-45)
%     plot_lon(i)=boat_lon(A(i))-0.01;
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<(90+45) && boat_heading(A(i))>(90-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
% elseif boat_heading(A(i))<(270+45) && boat_heading(A(i))>(270-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i))+0.01;    
% end
% end
% plot_lon=plot_lon-0.02;
% plot_lat=plot_lat+0.02;
% scatter(plot_lon,plot_lat,pointsize,boat_heading(A))
% 
% 
% %Stationary
% % A=find(t_tx_all>datenum(2017,6,8,10,0,0) & t_tx_all<datenum(2017,6,8,15,0,0));
% % scatter(boat_lon(A),boat_lat(A),pointsize,boat_heading(A))
% % A=find(t_tx_all>datenum(2017,6,9,0,10,0) & t_tx_all<datenum(2017,6,10,12,2,0));
% % scatter(boat_lon(A),boat_lat(A),pointsize,boat_heading(A))
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %Separate plots for each radial path
% %Lat/Lon of residual (offset)
% figure(234908)
% subplot(2,3,1)
% hold on
% grid on
% colormap jet
% set(gca,'Xlim',[-158.25 -157.725])
% set(gca,'Ylim',[22.55 22.95])
% 
% %Radials
% A=find(t_tx_all>datenum(2017,6,7,23,0,0) & t_tx_all<datenum(2017,6,8,10,0,0));      %First
% clear plot_lon plot_lat quiv_lon quiv_lat head_col
% for i=1:length(A)
% if boat_heading(A(i))<=(0+45) || boat_heading(A(i))>=(0+315)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
%     head_col(i)=1;
% elseif boat_heading(A(i))<=(180+45) && boat_heading(A(i))>=(180-45)
%     plot_lon(i)=boat_lon(A(i))-0.01;
%     plot_lat(i)=boat_lat(A(i));
%     head_col(i)=2;
% elseif boat_heading(A(i))<(90+45) && boat_heading(A(i))>(90-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
%     head_col(i)=3;
% elseif boat_heading(A(i))<(270+45) && boat_heading(A(i))>(270-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i))+0.01;
%     head_col(i)=4;
% end
% end
% scatter(plot_lon,plot_lat,pointsize,(d(A)-(G(A,:)*m)).*1000)
% quiver(plot_lon(1:10:end),plot_lat(1:10:end),sind(boat_heading(A(1:10:end))'),cosd(boat_heading(A(1:10:end))'),0.35,'k')
% 
% subplot(2,3,4)
% hold on
% grid on
% colormap jet
% set(gca,'Xlim',[min(t_tx_all(A))-0.01 max(t_tx_all(A))+0.01])
% set(gca,'Ylim',[-5 15])
% ylabel('Residual (ms)')
% for i=1:length(head_col)
% if head_col(i)==1   
%     scatter(t_tx_all(A(i)),(d(A(i))-(G(A(i),:)*m)).*1000,'.b')       %North
% elseif head_col(i)==2
%     scatter(t_tx_all(A(i)),(d(A(i))-(G(A(i),:)*m)).*1000,'.r')       %South
% elseif head_col(i)==3
%     scatter(t_tx_all(A(i)),(d(A(i))-(G(A(i),:)*m)).*1000,'.g')       %East
% elseif head_col(i)==4
%     scatter(t_tx_all(A(i)),(d(A(i))-(G(A(i),:)*m)).*1000,'.k')       %West
% end
% end
% datetick(gca,'keeplimits')
% 
% subplot(2,3,2)
% hold on
% grid on
% colormap jet
% set(gca,'Xlim',[-158.25 -157.725])
% set(gca,'Ylim',[22.55 22.95])
% A=find(t_tx_all>datenum(2017,6,10,12,2,0) & t_tx_all<datenum(2017,6,11,4,45,0));    %Second
% clear plot_lon plot_lat quiv_lon quiv_lat head_col
% for i=1:length(A)
% if boat_heading(A(i))<=(0+45) || boat_heading(A(i))>=(0+315)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
%     head_col(i)=1;
% elseif boat_heading(A(i))<=(180+45) && boat_heading(A(i))>=(180-45)
%     plot_lon(i)=boat_lon(A(i))-0.01;
%     plot_lat(i)=boat_lat(A(i));
%     head_col(i)=2;
% elseif boat_heading(A(i))<(90+45) && boat_heading(A(i))>(90-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
%     head_col(i)=3;
% elseif boat_heading(A(i))<(270+45) && boat_heading(A(i))>(270-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i))+0.01;
%     head_col(i)=4;
% end
% end
% scatter(plot_lon,plot_lat,pointsize,(d(A)-(G(A,:)*m)).*1000)
% quiver(plot_lon(1:10:end),plot_lat(1:10:end),sind(boat_heading(A(1:10:end))'),cosd(boat_heading(A(1:10:end))'),0.35,'k')
% 
% subplot(2,3,5)
% hold on
% grid on
% colormap jet
% set(gca,'Xlim',[min(t_tx_all(A))-0.01 max(t_tx_all(A))+0.01])
% set(gca,'Ylim',[-5 15])
% for i=1:length(head_col)
% if head_col(i)==1   
%     scatter(t_tx_all(A(i)),(d(A(i))-(G(A(i),:)*m)).*1000,'.b')       %North
% elseif head_col(i)==2
%     scatter(t_tx_all(A(i)),(d(A(i))-(G(A(i),:)*m)).*1000,'.r')       %South
% elseif head_col(i)==3
%     scatter(t_tx_all(A(i)),(d(A(i))-(G(A(i),:)*m)).*1000,'.g')       %East
% elseif head_col(i)==4
%     scatter(t_tx_all(A(i)),(d(A(i))-(G(A(i),:)*m)).*1000,'.k')       %West
% end
% end
% datetick(gca,'keeplimits')
% 
% subplot(2,3,3)
% hold on
% grid on
% colormap jet
% c=colorbar;
% c.Label.String='residual (ms)';
% set(gca,'Xlim',[-158.25 -157.725])
% set(gca,'Ylim',[22.55 22.95])
% A=find(t_tx_all>datenum(2017,6,11,21,3,0));             %Third
% clear plot_lon plot_lat quiv_lon quiv_lat head_col
% for i=1:length(A)
% if boat_heading(A(i))<=(0+45) || boat_heading(A(i))>=(0+315)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
%     head_col(i)=1;
% elseif boat_heading(A(i))<=(180+45) && boat_heading(A(i))>=(180-45)
%     plot_lon(i)=boat_lon(A(i))-0.01;
%     plot_lat(i)=boat_lat(A(i));
%     head_col(i)=2;
% elseif boat_heading(A(i))<(90+45) && boat_heading(A(i))>(90-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
%     head_col(i)=3;
% elseif boat_heading(A(i))<(270+45) && boat_heading(A(i))>(270-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i))+0.01;
%     head_col(i)=4;
% end
% end
% scatter(plot_lon,plot_lat,pointsize,(d(A)-(G(A,:)*m)).*1000)
% quiver(plot_lon(1:10:end),plot_lat(1:10:end),sind(boat_heading(A(1:10:end))'),cosd(boat_heading(A(1:10:end))'),0.35,'k')
% 
% subplot(2,3,6)
% hold on
% grid on
% colormap jet
% set(gca,'Xlim',[min(t_tx_all(A))-0.01 max(t_tx_all(A))+0.01])
% set(gca,'Ylim',[-5 15])
% for i=1:length(head_col)
% if head_col(i)==1   
%     l_1=scatter(t_tx_all(A(i)),(d(A(i))-(G(A(i),:)*m)).*1000,'.b');       %North
% elseif head_col(i)==2
%     l_2=scatter(t_tx_all(A(i)),(d(A(i))-(G(A(i),:)*m)).*1000,'.r');       %South
% elseif head_col(i)==3
%     l_3=scatter(t_tx_all(A(i)),(d(A(i))-(G(A(i),:)*m)).*1000,'.g');       %East
% elseif head_col(i)==4
%     l_4=scatter(t_tx_all(A(i)),(d(A(i))-(G(A(i),:)*m)).*1000,'.k');       %West
% end
% end
% datetick(gca,'keeplimits')
% legend([l_1 l_2 l_3 l_4],'N','S','E','W')
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %Separate plots for each radial path
% %Lat/Lon of residual (offset) vs. Range
% figure(234909)
% subplot(2,3,1)
% hold on
% grid on
% colormap jet
% set(gca,'Xlim',[-158.25 -157.725])
% set(gca,'Ylim',[22.55 22.95])
% 
% %Radials
% A=find(t_tx_all>datenum(2017,6,7,23,0,0) & t_tx_all<datenum(2017,6,8,10,0,0));      %First
% clear plot_lon plot_lat quiv_lon quiv_lat head_col
% for i=1:length(A)
% if boat_heading(A(i))<=(0+45) || boat_heading(A(i))>=(0+315)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
%     head_col(i)=1;
% elseif boat_heading(A(i))<=(180+45) && boat_heading(A(i))>=(180-45)
%     plot_lon(i)=boat_lon(A(i))-0.01;
%     plot_lat(i)=boat_lat(A(i));
%     head_col(i)=2;
% elseif boat_heading(A(i))<(90+45) && boat_heading(A(i))>(90-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
%     head_col(i)=3;
% elseif boat_heading(A(i))<(270+45) && boat_heading(A(i))>(270-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i))+0.01;
%     head_col(i)=4;
% end
% end
% scatter(plot_lon,plot_lat,pointsize,(d(A)-(G(A,:)*m)).*1000)
% quiver(plot_lon(1:10:end),plot_lat(1:10:end),sind(boat_heading(A(1:10:end))'),cosd(boat_heading(A(1:10:end))'),0.35,'k')
% 
% subplot(2,3,4)
% hold on
% grid on
% colormap jet
% %set(gca,'Xlim',[min(t_tx_all(A))-0.01 max(t_tx_all(A))+0.01])
% set(gca,'Ylim',[-5 15])
% ylabel('Residual (ms)')
% for i=1:length(head_col)
% if head_col(i)==1   
%     scatter(x_dist(A(i))./1000,(d(A(i))-(G(A(i),:)*m)).*1000,'.b')       %North
% elseif head_col(i)==2
%     scatter(x_dist(A(i))./1000,(d(A(i))-(G(A(i),:)*m)).*1000,'.r')       %South
% elseif head_col(i)==3
%     scatter(x_dist(A(i))./1000,(d(A(i))-(G(A(i),:)*m)).*1000,'.g')       %East
% elseif head_col(i)==4
%     scatter(x_dist(A(i))./1000,(d(A(i))-(G(A(i),:)*m)).*1000,'.k')       %West
% end
% end
% 
% subplot(2,3,2)
% hold on
% grid on
% colormap jet
% set(gca,'Xlim',[-158.25 -157.725])
% set(gca,'Ylim',[22.55 22.95])
% A=find(t_tx_all>datenum(2017,6,10,12,2,0) & t_tx_all<datenum(2017,6,11,4,45,0));    %Second
% clear plot_lon plot_lat quiv_lon quiv_lat head_col
% for i=1:length(A)
% if boat_heading(A(i))<=(0+45) || boat_heading(A(i))>=(0+315)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
%     head_col(i)=1;
% elseif boat_heading(A(i))<=(180+45) && boat_heading(A(i))>=(180-45)
%     plot_lon(i)=boat_lon(A(i))-0.01;
%     plot_lat(i)=boat_lat(A(i));
%     head_col(i)=2;
% elseif boat_heading(A(i))<(90+45) && boat_heading(A(i))>(90-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
%     head_col(i)=3;
% elseif boat_heading(A(i))<(270+45) && boat_heading(A(i))>(270-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i))+0.01;
%     head_col(i)=4;
% end
% end
% scatter(plot_lon,plot_lat,pointsize,(d(A)-(G(A,:)*m)).*1000)
% quiver(plot_lon(1:10:end),plot_lat(1:10:end),sind(boat_heading(A(1:10:end))'),cosd(boat_heading(A(1:10:end))'),0.35,'k')
% 
% subplot(2,3,5)
% hold on
% grid on
% colormap jet
% %set(gca,'Xlim',[min(t_tx_all(A))-0.01 max(t_tx_all(A))+0.01])
% set(gca,'Ylim',[-5 15])
% for i=1:length(head_col)
% if head_col(i)==1   
%     scatter(x_dist(A(i))./1000,(d(A(i))-(G(A(i),:)*m)).*1000,'.b')       %North
% elseif head_col(i)==2
%     scatter(x_dist(A(i))./1000,(d(A(i))-(G(A(i),:)*m)).*1000,'.r')       %South
% elseif head_col(i)==3
%     scatter(x_dist(A(i))./1000,(d(A(i))-(G(A(i),:)*m)).*1000,'.g')       %East
% elseif head_col(i)==4
%     scatter(x_dist(A(i))./1000,(d(A(i))-(G(A(i),:)*m)).*1000,'.k')       %West
% end
% end
% xlabel('Range (km)')
% 
% subplot(2,3,3)
% hold on
% grid on
% colormap jet
% c=colorbar;
% c.Label.String='residual (ms)';
% set(gca,'Xlim',[-158.25 -157.725])
% set(gca,'Ylim',[22.55 22.95])
% A=find(t_tx_all>datenum(2017,6,11,21,3,0));             %Third
% clear plot_lon plot_lat quiv_lon quiv_lat head_col
% for i=1:length(A)
% if boat_heading(A(i))<=(0+45) || boat_heading(A(i))>=(0+315)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
%     head_col(i)=1;
% elseif boat_heading(A(i))<=(180+45) && boat_heading(A(i))>=(180-45)
%     plot_lon(i)=boat_lon(A(i))-0.01;
%     plot_lat(i)=boat_lat(A(i));
%     head_col(i)=2;
% elseif boat_heading(A(i))<(90+45) && boat_heading(A(i))>(90-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i));
%     head_col(i)=3;
% elseif boat_heading(A(i))<(270+45) && boat_heading(A(i))>(270-45)
%     plot_lon(i)=boat_lon(A(i));
%     plot_lat(i)=boat_lat(A(i))+0.01;
%     head_col(i)=4;
% end
% end
% scatter(plot_lon,plot_lat,pointsize,(d(A)-(G(A,:)*m)).*1000)
% quiver(plot_lon(1:10:end),plot_lat(1:10:end),sind(boat_heading(A(1:10:end))'),cosd(boat_heading(A(1:10:end))'),0.35,'k')
% 
% subplot(2,3,6)
% hold on
% grid on
% colormap jet
% %set(gca,'Xlim',[min(t_tx_all(A))-0.01 max(t_tx_all(A))+0.01])
% set(gca,'Ylim',[-5 15])
% for i=1:length(head_col)
% if head_col(i)==1   
%     l_1=scatter(x_dist(A(i))./1000,(d(A(i))-(G(A(i),:)*m)).*1000,'.b');       %North
% elseif head_col(i)==2
%     l_2=scatter(x_dist(A(i))./1000,(d(A(i))-(G(A(i),:)*m)).*1000,'.r');       %South
% elseif head_col(i)==3
%     l_3=scatter(x_dist(A(i))./1000,(d(A(i))-(G(A(i),:)*m)).*1000,'.g');       %East
% elseif head_col(i)==4
%     l_4=scatter(x_dist(A(i))./1000,(d(A(i))-(G(A(i),:)*m)).*1000,'.k');       %West
% end
% end
% legend([l_1 l_2 l_3 l_4],'N','S','E','W')


%}

%% Massage data
m_recov1 = m_recov(1:tot_grid_num);
% m_p = m_recov(tot_grid_num+1:tot_grid_num+3);

% 6.3. Resolution Matrix
Res_mat = GI*G;
Res_mat1 = Res_mat(1:tot_grid_num,1:tot_grid_num);
% Res_mat_p  = Res_mat(tot_grid_num+1:tot_grid_num+3,tot_grid_num+1:tot_grid_num+3);

% 6.4. Posteriori covariance
P_post = (eye(tot_grid_num+3)-Res_mat)*P;
P_post1 = P_post(1:tot_grid_num,1:tot_grid_num);
% P_post_p = P_post(tot_grid_num+1:tot_grid_num+3,tot_grid_num+1:tot_grid_num+3);

% 6.5 reshape the SS field
recov_SS1=reshape(m_recov1(1:end),grid_num,grid_num)'  ;

% 6.6 residual
residual = d-G*m_recov;
norm_res = residual./sqrt(diag(Cd));

% 6.7 model resolution
model_res1 = reshape(diag(Res_mat1),grid_num,grid_num);
% model_res_p = reshape(diag(Res_mat_p),3,1);

% 6.8 variance reduction
prior_SD1 = reshape(sqrt(diag(P(1:tot_grid_num,1:tot_grid_num))),grid_num,grid_num);
% prior_p = reshape(sqrt(diag(P(tot_grid_num+1:tot_grid_num+3,tot_grid_num+1:tot_grid_num+3))),3,1);

post_SD1 = reshape(sqrt(diag(P_post1)),grid_num,grid_num);
% post_p = reshape(sqrt(diag(P_post_p)),3,1);

sd_reduction1 = real(post_SD1)./prior_SD1*100;
% sd_reduction_p = post_p./prior_p*100;

%% 7. plot the recovered ss pertrubation field

% draw circle
R = 25000;
num_point = 2000;
inc_ang =360/num_point;
x_cir = zeros(1,num_point);
y_cir = zeros(1,num_point);
for ii = 1:num_point
    angle = inc_ang*ii;
    [x_cir(ii),y_cir(ii),~]=m_fdist(icListen_lon,icListen_lat,angle,R); 
end
x_cir = x_cir-360;

figure(2)
clf
set(gcf,'Units','normalized','Position',[0.2 0.5 0.3 0.4]);
imagesc(x_cen,y_cen,recov_SS1)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title('Recovered SS Perturbation Field (1st mode)')
hold on
scatter(icListen_lon,icListen_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal')
% caxis([-7 7])



%% Solution Quality

%%% Plot
figure(3)
set(gcf,'Units','normalized','Position',[.3 1 0.5 .5])
subplot(2,2,1)
imagesc(x_cen,y_cen,real(sd_reduction1))
colormap jet
cbar = colorbar;
cbar.Label.String = '%';
title('RMS Error Reduction 1st mode')
hold on
scatter(icListen_lon,icListen_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal')
% caxis([0 100])
cax_lim = cbar.Limits;


subplot(2,2,2)
imagesc(3,1,sd_reduction_p)
colormap jet
cbar = colorbar;
cbar.Label.String = '%';
title('RMS Error Reduction Hydrophone Position')
hold on
scatter(icListen_lon,icListen_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal')
%%
% subplot(2,2,3)
imagesc(x_cen,y_cen,model_res1)
axis tight
grid off
title('Model Resolution 1st mode')
colorbar
% caxis([0 1])
% caxis([min(P_post(:)) max(P_post(:))])
colormap jet
%%

subplot(2,2,4)
imagesc(1,3,model_res_p)
axis tight
grid off
title('Model Resolution Hydrophone Offset')
colorbar
% caxis([0 1])
% caxis([min(P_post(:)) max(P_post(:))])
colormap jet


%% G matrix
figure(3)
plot(1:tot_grid_num,sum(G_SS,1))
xticks(1:20:tot_grid_num)
grid on
ylabel('meter')
xlabel('pixel')
xlim([1 tot_grid_num])
title('Obserbvation Matrix of the first mode')


%% Travel Time Errors
d_diff=d-(G*m_recov);
d_error=sqrt(diag(Cd));
d_normalized=d_diff./diag(d_error);
figure(324)
scatter(d_normalized,x_dist./1000,'.')
xlabel('Data Error Normalized')
ylabel('Surface Distance (km)')
grid on



%Plot actual and corrected d
figure(8544)
hold on
grid on
plot(d.*1000)
plot((d-(G*m_recov)).*1000,'r')
ylabel('(ms)')
xlabel('Transmission #')
title('Actual vs Estimated Travel Times')
legend('Original','After Inversion')







%%

save inverse_soluiton_21x21_20km_V_code G_SS G m_recov Cd d P P_post




