%{
% % Check ray tracing
% % RTX/POSMV on HEM
% close all
% clear
% addpath('/Users/testuser/Documents/MATLAB/Script/Ray_Trace_Dzieciuch/Ray_Trace_Comparison')
% ACO_lat=22.738772;                  % June 2017
% ACO_lon=-158.006186;                % June 2017
% ACO_altitude = 4736.266;         % at 4,736.266 m June 2017
% 
% day = 28:28;            %  Edit
% hour = 17:23;              % Edit
% %% Load tx files
% [tx_t,tx_lat,tx_lon,tx_altitude,tx_heading,range,tx_xvel] = posmv_tx_load(day,hour);
% 
% 
% %% Load rx files
% % POS MV
% cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/rx_file/HEM/all/ideal_pulse
% fname = [];
% start_hour = hour(1);
% end_hour = hour(end);
% now_hour = start_hour;
% now_day = day(1);
% while true
%     while now_hour <= 23
%         if now_hour < 10
%             fname_d = "rx_data_2018_10_"+string(now_day)+"_0"+string(now_hour)+"_i";
%         else
%             fname_d = "rx_data_2018_10_"+string(now_day)+"_"+string(now_hour)+"_i";
%         end
%         fname = vertcat(fname,fname_d);
%         now_hour = now_hour+1;
%         
%         if (now_hour > end_hour)& (now_day == day(end))
%             break;
%         end
%     end
%     
%     if (now_hour > end_hour)& (now_day == day(end))
%         break;
%     end
%     now_hour = 0;
%     now_day = now_day+1;
% end
% 
% act_arrival = [];
% SNR = [];
% for ii = 1:length(fname)
%     load(fname(ii)+".mat")
%     act_arrival = horzcat(act_arrival ,rx_data.act_arrival);
%     SNR = horzcat(SNR ,rx_data.SNR);
% end
% 
% % RTX
% cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/RTX/rx_file/corrected/ideal_pulse
% fname = [];
% start_hour = hour(1);
% end_hour = hour(end);
% now_hour = start_hour;
% now_day = day(1);
% 
% while true
%     while now_hour <= 23
%         if now_hour < 10
%             fname_d = "rtx_rx_data_2018_10_"+string(now_day)+"_0"+string(now_hour);
%         else
%             fname_d = "rtx_rx_data_2018_10_"+string(now_day)+"_"+string(now_hour);
% 
%         end
%        fname = vertcat(fname,fname_d);
%         now_hour = now_hour+1;
% 
%         if (now_hour > end_hour)& (now_day == day(end))
%             break;
%         end
%     end
%     
%     if (now_hour > end_hour)& (now_day == day(end))
%         break;
%     end
%     now_hour = 0;
%     now_day = now_day+1;
% end
% 
% 
% rtx_act_arrival = [];
% rtx_SNR = [];
% for ii = 1:length(fname)
%     load(fname(ii)+".mat")
%     rtx_act_arrival = horzcat(rtx_act_arrival ,rtx_rx_data.act_arrival);
%     rtx_SNR = horzcat(rtx_SNR ,rtx_rx_data.SNR);
% 
% end
% 
% %% Clean Data
% % POS MV
% rm_ind = [];
% for l = 1:length(tx_t)
%     dum_arrival = act_arrival;
%     dum_arrival(find(dum_arrival < tx_t(l))) = [];         
%     h_ind = find((tx_t(l) - dum_arrival)*3600*24 > -20);   
%     if isempty(h_ind)
%         rm_ind(end+1) = l;
%     end
% end
% keep_ind = setdiff(1:length(tx_t),rm_ind);
% 
% tx_t =tx_t(keep_ind);
% tx_lat = tx_lat(keep_ind);
% tx_lon = tx_lon(keep_ind);
% tx_altitude  = tx_altitude(keep_ind);
% tx_heading = tx_heading(keep_ind);
% tx_xvel = tx_xvel(keep_ind);
% range = range(keep_ind);
% 
% 
% rm_ind = [];
% for l = 1:length(act_arrival)
%     dum_tx_t = tx_t;
%     dum_tx_t(find( dum_tx_t > act_arrival(l))) = [];
%     % if greater than 20 sec = no reception was picked up 
%     h_ind = find((act_arrival(l) - dum_tx_t)*3600*24 < 20);
%     if isempty(h_ind)
%         rm_ind(end+1) = l;
%     end
% end
% act_arrival(rm_ind) = [];
% SNR(rm_ind) = [];
% 
% 
% %{ 
% % 
% % % RTX
% % for l = 1:length(rtx_tx_t)
% %     dum_arrival = rtx_act_arrival;
% %     dum_arrival(find(dum_arrival < rtx_tx_t(l))) = [];         
% %     h_ind = find((rtx_tx_t(l) - dum_arrival)*3600*24 > -20);   
% %     if isempty(h_ind)
% %         rm_ind(end+1) = l;
% %     end
% % end
% % keep_ind = setdiff(1:length(rtx_tx_t),rm_ind);
% % 
% % rtx_tx_t = rtx_tx_t(keep_ind);
% % rtx_tx_lat = rtx_tx_lat(keep_ind);
% % rtx_tx_lon = rtx_tx_lon(keep_ind);
% % rtx_tx_altitude  = rtx_tx_altitude(keep_ind);
% % rtx_tx_heading = rtx_tx_heading(keep_ind);
% % rtx_tx_xvel = rtx_tx_xvel(keep_ind);
% % rtx_range = rtx_range(keep_ind);
% % 
% % 
% % rm_ind = [];
% % for l = 1:length(rtx_act_arrival)
% %     dum_tx_t = rtx_tx_t;
% %     dum_tx_t(find( dum_tx_t > rtx_act_arrival(l))) = [];
% %     % if greater than 20 sec = no reception was picked up 
% %     h_ind = find((rtx_act_arrival(l) - dum_tx_t)*3600*24 < 20);
% %     if isempty(h_ind)
% %         rm_ind(end+1) = l;
% %     end
% % end
% % rtx_act_arrival(rm_ind) = [];
% % rtx_SNR(rm_ind) = [];
% %}
% %% Match POS MV/ RTX
% %{
% rm_indposmv = [];
% rm_indrtx = [];
% t_diffa = [];
% for i = 1:length(tx_t)
%     [t_diff,~] = min(abs(tx_t(i) - rtx_tx_t)); 
%     if t_diff ~= 0 
%         rm_indposmv(end+1) = i;
% %         t_diffa(end+1) = t_diff*3600*24;
%     end
% end
% tx_t(rm_indposmv) = [];
% tx_lat(rm_indposmv) = []; 
% tx_lon(rm_indposmv) = [];
% tx_altitude(rm_indposmv) = [];
% range(rm_indposmv) = [];
% tx_heading(rm_indposmv) = [];
% tx_xvel(rm_indposmv) = [];
% act_arrival(rm_indposmv) = [];
% SNR(rm_indposmv) = [];
% %{
% for i = 1:length(rtx_tx_t)
%     [t_diff,~] = min(abs(rtx_tx_t(i) - tx_t)); 
%     if t_diff ~= 0 
%         rm_indrtx(end+1) = i;
% %         t_diffa(end+1) = t_diff*3600*24;
%     end
% end
% rtx_tx_t(rm_indrtx) = [];
% rtx_tx_lat(rm_indrtx) = [];
% rtx_tx_lon(rm_indrtx) = [];
% rtx_tx_altitude(rm_indrtx) = [];
% rtx_range(rm_indrtx) = [];
% rtx_tx_heading(rm_indrtx) = [];
% rtx_tx_xvel(rm_indrtx) = [];
% rtx_act_arrival(rm_indrtx) = [];
% rtx_SNR(rm_indrtx) = [];
% %}
% %}
% %% UTM coordinate
% [ACO_x,ACO_y,utmzone] = deg2utm(ACO_lat,ACO_lon);
% [tx_x,tx_y,~] = deg2utm(tx_lat',tx_lon');
% %[rtx_tx_x,rtx_tx_y,~] = deg2utm(rtx_tx_lat',rtx_tx_lon');
% 
% %% Correct for Earth Curvature (Hydrophone)
% 
% E_radius=6371000;               % epsilonray tracing APL code
% % coordinate transformation
% E = tx_altitude./E_radius;
% %rtx_E = rtx_tx_altitude./E_radius;
% ACO_E = ACO_altitude./E_radius;
% tx_z= tx_altitude.*(1+(E./2)+((E.^2)./3));
% %rtx_tx_z = rtx_tx_altitude.*(1+(rtx_E./2)+((rtx_E.^2)./3));
% ACO_z = ACO_altitude.*(1+(ACO_E./2)+((ACO_E.^2)./3));
% 
% 
% 
% %% Calculate estimated Time
% range = [0.5 1 5 10 20];
% tx_altitude = ones(1,length(range))*-6.6;
% % POS MV
%  for ii= 1:length(range)
%      
%     [arc_lengths,tot_dist(ii),theta0_1,SS,z,SS_HOT_avg(ii),surface_dist,inc_r,est_tt(ii),~] = ray_tracing_const_SS(range(ii)*1000,tx_altitude(ii));
%     [arc_lengths2,tot_dist2(ii),theta0_2,SS2,z2,SS_HOT_avg2(ii),surface_dist2,inc_r2,est_tt2(ii),~] = ray_tracing_const_SS_gradient(range(ii)*1000,tx_altitude(ii));
%     
%     save SS_profile SS z
%     % z with spherical earth correction/ SS with corrected z 
%     %[est_tt_D(i),theta0_D(i)]=ray_trace_Dzieciuch(z,SS,-tx_altitude(i),range(i)*1000);
%     
%     %[rtx_arc_length,rtx_tot_dist(i),~,rtx_SS,rtx_z,rtx_SS_HOT_avg(i),rtx_surface_dist,rtx_inc_r,rtx_est_tt(i),~]=ray_tracing_const_SS(rtx_range(i)*1000,rtx_tx_altitude(i));
%  
%     % tot_arc_diff(i) = tot_dist2(i)-tot_dist(i);
%     
%     %{
% %     figure(1)
% %     clf
% %     if tx_x(i) >ACO_x
% %         scatter(0,ACO_z,200,'kp','filled')
% %         hold on
% %         plot(range(i)-inc_r/1000,z2,'b')
% %         plot(range(i)-inc_r2/1000,z2,'g')
% %        % plot(rtx_range(i)-rtx_inc_r/1000,rtx_z,'r')
% %         set(gca,'YDir','reverse')
% %         xlim([-max(range) max(range)])
% %         ylim([0 ACO_z])
% %         grid on
% %         xlabel('Range (km)')
% %         ylabel('Depth (m)')
% %     else
% %         scatter(0,ACO_z,200,'kp','filled')
% %         hold on
% %         plot(-range(i)+inc_r/1000,z2,'b')
% %         plot(-range(i)+inc_r2/1000,z2,'g')
% %         %plot(-rtx_range(i)+rtx_inc_r/1000,rtx_z,'r')
% %         set(gca,'YDir','reverse')
% %         xlim([-max(range) max(range)])
% %         ylim([0 ACO_z])
% %         grid on
% %         xlabel('Range (km)')
% %         ylabel('Depth (m)')
% %     end
%     
%     %}
%     %}
%  end
%  %%
% % est_arrival = est_tt/3600/24 +tx_t;
% % est_arrival2 = est_tt2/3600/24 +tx_t;
% est_tt_diff = est_tt2-est_tt ;
% tot_arc_diff = tot_dist2-tot_dist;
% %est_tt_diff_D = est_tt2 - est_tt_D;
% %% Travel time perturbation 
% ttp = (act_arrival - est_arrival)*3600*24*1000;
% ttp2 = (act_arrival - est_arrival2)*3600*24*1000;
% %% Estimated travel time difference
% figure(3)
% scatter(range,est_tt_diff*1000,10,'filled')
% ylabel('millisecond')
% title('Estimated Travel Time Differecne (New model - Old Model)')
% xlabel('Range (km)')
% grid on
%  %% Arc Length Difference
% figure(2)
% scatter(range,tot_arc_diff*100,10,'filled')
% ylabel('cm')
% grid on
% xlabel('Range (km)')
% title('Total Arc Length Difference (New Model - Old Model)')
% axis tight
%}
%%%%%
% x_dist = random walk until reah 25 km
% z_offset = random between -7 to -5

%%%%%%%%%%%%

clear 
dist_now = 10;
z_offset0 = -6;
iii = 1;
while dist_now <= 25000
    dist_now = dist_now+ rand*100+200;
    x_dist(iii) = dist_now;
    z_offset(iii)  = 0;
%     z_offset(iii) = z_offset0+(2*rand-7);
    [arc_lengths1,tot_dist1,theta01,SS1,z1,SS_HOT_avg1,surface_dist1,inc_r1,est_tt1,ray_angles1 ] = ray_tracing_const_SS( x_dist(iii),z_offset(iii) );
    [arc_lengths2,tot_dist2,theta02,SS2,z2,SS_HOT_avg2,surface_dist2,inc_r2,est_tt2,ray_angles2] = ray_tracing_const_SS_gradient( x_dist(iii),z_offset(iii) );
    diff_t(iii) = est_tt1 - est_tt2;
    diff_arc(iii) = tot_dist1 - tot_dist2;
     iii = iii+1;
   
end
%%
figure(1)
clf
scatter(x_dist/1000,diff_t*1000,'.')
grid on
xlabel('Range (km)')
ylabel('Time difference (ms)')
title('Comparison (Constant SS - Varying SS)')
xlim([0 25])
hold on
yyaxis right
scatter(x_dist/1000,diff_arc*100,'.')
ylabel('Arc Length Difference (cm)')

%% Ray tracing
function [arc_lengths,tot_dist,theta0,SS,z,SS_HOT_avg,surface_dist,inc_r,est_tt,ray_angles ] = ray_tracing_const_SS( x_dist,z_offset )
% arc_length=
% tot_dist  =
% theta0    = launch angle
% SS        = Sound speed profile (with depth z)
% z         = Depth (negative/ below the transducer)
% SS_HOT    = 
% surface_dist = 
% est_tt    = estimated travel time
% ray angle = 
cd '/Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Script/Data'
%fid=fopen('h294a0202.ctd');                       % June 2018
fid=fopen('h306a0202.ctd');                       % October 2018

D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f','headerlines',6));
fclose(fid);
        
%z_dist=4812.7+.669;      %Paroscientific pressure sensor (dbar) + hyd depth 4729.92 m 
z_dist = 4819.897;     % at 4,736.266 m June 2017

pres=D(:,1);        %pressure (dbars)
temp=D(:,2);        %temperature
sal=D(:,3);         %salinity


%% Fix z offset for transducer location (m=dbar at this level)
% Find depths above the onboard transducer depth, and discard them and
% associated data points
pres_offset = gsw_p_from_z(z_offset,22.738772);
A=find(pres_offset>=pres);
pres(A(1:end-1))=[];
pres_diff= pres_offset-pres(1);
pres(1)= pres_offset;      % set the first depth to be the transducer depth
temp(A(1:end-1))=[];
temp(1)=(temp(1)*(2-pres_diff)+temp(2)*pres_diff)/2;   % interpolate the first temp and salinity
sal(A(1:end-1))=[];
sal(1)=(sal(1)*(2-pres_diff)+sal(2)*pres_diff)/2;

%Remove/Add deeper depths
if pres(end)>z_dist
remove_points=find(pres>z_dist);
pres(remove_points)=[];
temp(remove_points)=[];
sal(remove_points)=[];
else
    
    add_len = floor((z_dist - pres(end))/2);
    pres(end+1:end+1+add_len)=[pres(end)+transpose(2:2:add_len*2) ;z_dist];
    temp(end+1:end+1+add_len)=ones(add_len+1,1)*temp(end);
    sal(end+1:end+1+add_len)=ones(add_len+1,1)*sal(end);
    %sal(end+1)=sal(end)+(sal(end)-sal(end-1));
end

%A priori SS field
SS=gsw_sound_speed(sal,temp,pres);

z=-gsw_z_from_p(pres,22.738772);

% %Correct for Earth Curvature
E_radius=6371000;        % epsilonray tracing APL code
% coordinate transformation
E=z./E_radius;
z=z.*(1+(E./2)+((E.^2)./3));
%z(end)=hyd_depth;
SS=SS.*(1+E+E.^2);

%Find theta0 of acoustic transmission
% theta0_all=0.1:0.01*(pi/180):89*(pi/180);
% incident angle
theta0_all=0.0001:1*(pi/180):89/180*pi;   % all possible launch angles
theta0_all(end)=89/180*pi;

% loop over all possible thetas 
% calculate a small horizontal distance within an interval betweenf two depths
for ii=1:length(theta0_all)
    
    %Ray parameter
    a=sin(theta0_all(ii))/SS(1);
    
    % incident angles of each layer
    theta(1) = theta0_all(ii);
    for i=2:length(SS)-1
        theta(i)=asin(a*SS(i));
    end    
    % SS gradient in each layer for the slowly-varying ss profile
%     for i=1:length(SS)-1
%         b(i)=(SS(i+1)-SS(i))/(z(i+1)-z(i));
%     end
    
    % horizontl distance
    for i=1:length(SS)-1
        r(i)=(z(i+1)-z(i))*tan(theta(i));
    end
    
    %Total x distance
    x_diff(ii)=sum(r);  % total horizontal distance for one possible launching theta
    
    
end


%Find min and max launch angle for Netwon-Raph method 

t_low=find(x_diff<x_dist);      % find all launch angles which yield x distance within the bound
theta_low=theta0_all(t_low(end));   % label the lowest possible launch angle (correspond to the steepest angle)
x_low=x_diff(t_low(end));           % The shortest ray arc length, lower bound of the sur dist

t_high=find(x_diff>x_dist);
theta_high=theta0_all(t_high(1));   % The upper bound of sur dist (the least steep angle)
x_high=x_diff(t_high(1));

% Intrapolation to find the optimum ray arc length of the direct path
theta_new=theta_low+((x_dist)/((x_low-x_high)/(theta_low-theta_high)))-((x_low)/((x_low-x_high)/(theta_low-theta_high)));

% We now have a range of launch angles which gives the nearest number of
% horizontal distance 
theta1=theta_low;
x1=x_low;           % lower bound of the horizontal distance

x_diff=1000000;

% Loop until x_diff is close to x_dist to 0.1 mm 
while abs(x_diff-x_dist)>0.0001
    
    %Ray parameter
    a=sin(theta_new)/SS(1);
    
    % incident angle
    theta(1) = theta_new;
    for i=2:length(SS)-1
        theta(i)=asin(a*SS(i));
    end
    
    r = [];
    % horizontal range in each layer
    for i=1:length(SS)-1
        r(i)=(z(i+1)-z(i))*tan(theta(i));
    end
    
    %Total x distance
    x_diff=sum(r);
    
    %Total distance and time traveled
    for i=1:length(SS)-1
        arc_distance(i)=(z(i+1)-z(i))/cos(theta(i));
        tt(i)=(z(i+1)-z(i))/(SS(i)*cos(theta(i)));
    end
    tot_dist_all=sum(arc_distance);
    
    arc_lengths_hold=arc_distance;
    surface_dist_hold=r;
    ray_angles_hold=theta(2:end);
    
    %Newton-Raphson method
    % Interpolate again to find a new launch angle which is closer to that
    % of the eigenray
    
    theta2=theta_new;
    
    
    x2=x_diff;
    
    theta_new=theta2+((x_dist)/((x2-x1)/(theta2-theta1)))-((x2)/((x2-x1)/(theta2-theta1)));
    
    % Set theta_low to theta 2 (new angle)
    theta1=theta2;
    % update nearest horizontal distance
    x1=x2;
    
end


%% packaging outputs
theta0=theta_new;                        %Launch angle
tot_dist=tot_dist_all;                 %Ray arc length
arc_lengths=arc_lengths_hold;           %arc length at each interval
surface_dist=surface_dist_hold;         %Surface distance
ray_angles=ray_angles_hold;             %ray angle at each interval

% incremental horizontal distance
inc_r = zeros(1,length(surface_dist)+1);

for i= 2:length(inc_r)
    inc_r(i) = inc_r(i-1)+surface_dist(i-1);
end

clear x_diff tot_dist_all arc_lengths_hold

est_tt=sum(tt);       %Estimated TT

%Vertically averaged SS for HOTS CTD
SS_HOT_avg=mean(SS);


%Make variables same size
% if length(A)>1
%     arc_lengths=horzcat(zeros(1,length(A)-1),arc_lengths);
%     surface_dist=horzcat(zeros(1,length(A)-1),surface_dist);
%     ray_angles=horzcat(ones(1,length(A)-1).*ray_angles(1,1),ray_angles);
% end

end

function [arc_lengths,tot_dist,theta0,SS,z,SS_HOT_avg,surface_dist,inc_r,est_tt,ray_angles ] = ray_tracing_const_SS_gradient( x_dist,z_offset )
% arc_length= layerwise arc length
% tot_dist  = total arc length of the eigen ray
% theta0    = launch angle
% SS        = Sound speed profile (with depth z)
% z         = Depth (negative/ below the transducer)
% SS_HOT    = Average sound speed
% surface_dist = Total horizontal distance traveled by the ray
% est_tt    = estimated travel time
% ray angle = eigen angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd '/Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Script/Data'
%fid=fopen('h294a0202.ctd');                       % June 2018
fid=fopen('h306a0202.ctd');                       % October 2018

D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f','headerlines',6));
fclose(fid);
        
%z_dist=4812.7+.669;      %Paroscientific pressure sensor (dbar) + hyd depth 4729.92 m 
z_dist = 4819.897;     % at 4,736.266 m June 2017

pres=D(:,1);        %pressure (dbars)
temp=D(:,2);        %temperature
sal=D(:,3);         %salinity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sound speed profile preparation

% Fix z offset for transducer location (m=dbar at this level)
% Find depths above the onboard transducer depth, and discard them and
% associated data points
pres_offset = gsw_p_from_z(z_offset,22.738772);
A=find(pres_offset>=pres);
pres(A(1:end-1))=[];
pres_diff= pres_offset-pres(1);
pres(1)= pres_offset;      % set the first depth to be the transducer depth
temp(A(1:end-1))=[];
temp(1)=(temp(1)*(2-pres_diff)+temp(2)*pres_diff)/2;   % interpolate the first temp and salinity
sal(A(1:end-1))=[];
sal(1)=(sal(1)*(2-pres_diff)+sal(2)*pres_diff)/2;

%Remove/Add deeper depths
if pres(end)>z_dist
remove_points=find(pres>z_dist);
pres(remove_points)=[];
temp(remove_points)=[];
sal(remove_points)=[];
else
    
    add_len = floor((z_dist - pres(end))/2);
    pres(end+1:end+1+add_len)=[pres(end)+transpose(2:2:add_len*2) ;z_dist];
    temp(end+1:end+1+add_len)=ones(add_len+1,1)*temp(end);
    sal(end+1:end+1+add_len)=ones(add_len+1,1)*sal(end);
    %sal(end+1)=sal(end)+(sal(end)-sal(end-1));
end

%A priori SS field
SS=gsw_sound_speed(sal,temp,pres);

z=-gsw_z_from_p(pres,22.738772);

%Correct for Earth Curvature
E_radius=6371000;        % epsilonray tracing APL code
% coordinate transformation
E=z./E_radius;
z=z.*(1+(E./2)+((E.^2)./3));
SS=SS.*(1+E+E.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the eigen ray
%Find theta0 of acoustic transmission

% incident angle
theta0_all=0.0001:1*(pi/180):89/180*pi;   % all possible launch angles
theta0_all(end)=89/180*pi;

% loop over all possible thetas 
% calculate a small horizontal distance within an interval betweenf two depths

for ii=1:length(theta0_all)
    r = [];
    b = [];
    
    %Ray parameter
    a=sin(theta0_all(ii))/SS(1);
    
    % SS gradient in each layer for the slowly-varying ss profile
    for i= 1:length(SS) - 1
        b(i)=(SS(i+1)-SS(i))/(z(i+1)-z(i));
     end

    % horizontl distance
    for j= 1:length(SS) - 1
        v = SS(j+1)/SS(j);
        h = (sqrt(1-(a*SS(j+1))^2)+1)/(sqrt(1-(a*SS(j))^2)+1);
        r(j) = (1/(a*b(j)))*(sqrt(1-(a*SS(j))^2) - sqrt(1-(a*SS(j+1))^2));
       
    end   

    
    %Total horizontal distance to reach the hydrophone
    x_diff(ii)=sum(r);  
    
    
end


%Find min and max launch angle for Netwon-Raph method 

t_low=find(x_diff<x_dist);      % find all launch angles which yield x distance within the bound
theta_low=theta0_all(t_low(end));   % label the lowest possible launch angle (correspond to the steepest angle)
x_low=x_diff(t_low(end));           % The shortest ray arc length, lower bound of the sur dist

t_high=find(x_diff>x_dist);
theta_high=theta0_all(t_high(1));   % The upper bound of sur dist (the least steep angle)
x_high=x_diff(t_high(1));

% Intrapolation to find the optimum ray arc length of the direct path
theta_new=theta_low+((x_dist)/((x_low-x_high)/(theta_low-theta_high)))-((x_low)/((x_low-x_high)/(theta_low-theta_high)));

% We now have a range of launch angles which gives the nearest number of
% horizontal distance 
theta1=theta_low;
x1=x_low;           % lower bound of the horizontal distance

x_diff=1000000;

% Loop until x_diff is close to x_dist to 0.1 mm 
while abs(x_diff-x_dist)>0.0001
    r = [];         % horizontal distance in each layer
    t = [];         % elapsed time in each layer
    b = [];         % ss gradient in each layer
    
    %Ray parameter
    a=sin(theta_new)/SS(1);
    
    % incident angle at each interface
    theta(1) = theta_new;
    for i= 2:length(SS)
        theta(i)=asin(a*SS(i));
    end
    
    % SS gradient in each layer for the slowly-varying ss profile
    for i= 1:length(SS)-1
        b(i)=(SS(i+1)-SS(i))/(z(i+1)-z(i));
    end

% horizontl distance, arc length and travel time
    for j= 1:length(SS)-1
        v = SS(j+1)/SS(j);
        h = (sqrt(1-(a*SS(j+1))^2)+1)/(sqrt(1-(a*SS(j))^2)+1);
       t(j) = (1/b(j))*(log(v)-log(h));
       r(j) = (1/(a*b(j)))*(sqrt(1-(a*SS(j))^2) - sqrt(1-(a*SS(j+1))^2));
       arc_distance(j)= 1/abs(a*b(j))*abs(theta(j+1) - theta(j));
    end   
   
    % Total ditances
    x_diff=sum(r);                  % horizontal distance
    tot_dist_all=sum(arc_distance); % arc length
    
    % Layerwise distance/ angles
    arc_lengths_hold=arc_distance;
    surface_dist_hold=r;
    ray_angles_hold=theta(2:end);
    
    % %%%% Newton-Raphson method
    % Interpolate again to find a new launch angle which is closer to that
    % of the eigenray
    
    theta2=theta_new;
    
    
    x2=x_diff;
    
    theta_new=theta2+((x_dist)/((x2-x1)/(theta2-theta1)))-((x2)/((x2-x1)/(theta2-theta1)));
    
    % Set theta_low to theta 2 (new angle)
    theta1=theta2;
    % update nearest horizontal distance
    x1=x2;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% packaging outputs
theta0=theta_new;                        %Launch angle
tot_dist=tot_dist_all;                 %Ray arc length
arc_lengths=arc_lengths_hold;           %arc length at each interval
surface_dist=surface_dist_hold;         %Surface distance
ray_angles=ray_angles_hold;             %ray angle at each interval

% incremental horizontal distance
inc_r = zeros(1,length(surface_dist)+1);

for i= 2:length(inc_r)
    inc_r(i) = inc_r(i-1)+surface_dist(i-1);
end

clear x_diff tot_dist_all arc_lengths_hold

est_tt=sum(t);       %Estimated TT

%Vertically averaged SS for HOTS CTD
SS_HOT_avg=mean(SS);


%Make variables same size
% if length(A)>1
%     arc_lengths=horzcat(zeros(1,length(A)-1),arc_lengths);
%     surface_dist=horzcat(zeros(1,length(A)-1),surface_dist);
%     ray_angles=horzcat(ones(1,length(A)-1).*ray_angles(1,1),ray_angles);
% end

end




%% TX file Load
function [tx_t,tx_lat,tx_lon,tx_altitude,tx_heading,x_dist,tx_xvel] = posmv_tx_load(day,hour)
    %% Load Tx files
    cd('/Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/tx_file/all')
    % create a set of file names
    %ACO LAT/LON
    ACO_lat=22.738772;                  % Aug 2018
    ACO_lon=-158.006186;                % Aug 2018
    
    fname = [];
    start_hour = hour(1);
    end_hour = hour(end);
    now_hour = hour(1);
    now_day = day(1);
    while true
        while now_hour <= 23
            if now_hour < 10
                fname_d = "tx_data_2018_10_"+string(now_day)+"_0"+string(now_hour);
            else
                fname_d = "tx_data_2018_10_"+string(now_day)+"_"+string(now_hour);
            end
           fname = vertcat(fname,fname_d);
           now_hour = now_hour+1;
           if (now_hour > end_hour)& (now_day == day(end))
                break;
           end
        end
            if (now_hour > end_hour)& (now_day == day(end))
                break;
           end
            now_hour = 0;
            now_day = now_day+1;
    end
   fname
    %Load TX data
    tx_t = [];
    tx_lat = [];
    tx_lon = [];
    tx_altitude = [];
    tx_heading = [];
    tx_xvel = [];
   for p = 1:length(fname)
        load(fname(p))
        tx_t = horzcat(tx_t,tx_data.t);      %TX time                                  
        tx_lat = horzcat(tx_lat,tx_data.lat);   %TX lat
        tx_lon = horzcat(tx_lon,tx_data.lon);   
        tx_altitude = horzcat(tx_altitude,tx_data.altitude);   
        tx_heading = horzcat(tx_heading,tx_data.heading);
        tx_xvel = horzcat(tx_xvel,tx_data.x_vel);
   end

    % offset t_tx for HEM time offset
    % correct for the reception time
        date_mark1= "20181028 00:50";
        date_mark2= "20181029 00:50";
        date_mark3= "20181030 00:50";
        date_mark1 = datenum(date_mark1,'yyyymmdd HH:MM');
        date_mark2 = datenum(date_mark2,'yyyymmdd HH:MM');
        date_mark3 = datenum(date_mark3,'yyyymmdd HH:MM');

        for p = 1:length(tx_t)
            if tx_t(p) <= date_mark1
                tx_t(p) = tx_t(p) +5/(3600*24);
                
            elseif (date_mark1 <= tx_t(p))&(tx_t(p) <= date_mark2)
                tx_t(p) = tx_t(p)+6/(3600*24);
                
            elseif (date_mark2 <= tx_t(p))&(tx_t(p) <= date_mark3)
                tx_t(p) = tx_t(p)+7/(3600*24);
                
            elseif  date_mark3 <= tx_t(p)
                tx_t(p) = tx_t(p)+8/(3600*24);
                
            end
        end

    %Boat distance from ACO
    for i=1:length(tx_lat)
        x_dist(i)=dist([ACO_lat tx_lat(i)],[ACO_lon tx_lon(i)])/1000;
    end

    
    
end

function [tx_t,tx_lat,tx_lon,tx_altitude,tx_heading,x_dist,tx_xvel] = posmv_rtx_tx_load(day,hour)
    %% Load TX files
    cd('/Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/RTX/tx_file/corrected')
    % create a set of file names
    %ACO LAT/LON
    ACO_lat=22.738772;                  % Aug 2018
    ACO_lon=-158.006186;                % Aug 2018
    fname = [];
    start_hour = hour(1);
    end_hour = hour(end);
    now_hour = start_hour;
    now_day = day(1);
    while true
        while now_hour <= 23
            if now_hour < 10
                fname_d = "rtx_tx_data_2018_10_"+string(now_day)+"_0"+string(now_hour);
            else
                fname_d = "rtx_tx_data_2018_10_"+string(now_day)+"_"+string(now_hour);
            end
           fname = vertcat(fname,fname_d);
           now_hour = now_hour+1;
           if (now_hour > end_hour)& (now_day == day(end))
                break;
           end
        end
            if (now_hour > end_hour)& (now_day == day(end))
                break;
           end
            now_hour = 0;
            now_day = now_day+1;
    end

    %Load TX data
    tx_t = [];
    tx_lat = [];
    tx_lon = [];
    tx_altitude = [];
    tx_heading = [];
    tx_xvel = [];
   for p = 1:length(fname)
        load(fname(p))
        tx_t = horzcat(tx_t,rtx_tx_data.t);      %TX time                                  
        tx_lat = horzcat(tx_lat,rtx_tx_data.lat);   %TX lat
        tx_lon = horzcat(tx_lon,rtx_tx_data.lon);   
        tx_altitude = horzcat(tx_altitude,rtx_tx_data.altitude);   
        tx_heading = horzcat(tx_heading,rtx_tx_data.heading);
        tx_xvel = horzcat(tx_xvel,rtx_tx_data.xvel);
   end

    % offset t_tx for HEM time offset
    % correct for the reception time
        date_mark1= "20181028 00:50";
        date_mark2= "20181029 00:50";
        date_mark3= "20181030 00:50";
        date_mark1 = datenum(date_mark1,'yyyymmdd HH:MM');
        date_mark2 = datenum(date_mark2,'yyyymmdd HH:MM');
        date_mark3 = datenum(date_mark3,'yyyymmdd HH:MM');

        for p = 1:length(tx_t)
            if tx_t(p) <= date_mark1
                tx_t(p) = tx_t(p) +5/(3600*24);
                
            elseif (date_mark1 <= tx_t(p))&(tx_t(p) <= date_mark2)
                tx_t(p) = tx_t(p)+6/(3600*24);
                
            elseif (date_mark2 <= tx_t(p))&(tx_t(p) <= date_mark3)
                tx_t(p) = tx_t(p)+7/(3600*24);
                
            elseif  date_mark3 <= tx_t(p)
                tx_t(p) = tx_t(p)+8/(3600*24);
                
            end
        end

    %Boat distance from ACO
    for i=1:length(tx_lat)
        x_dist(i)=dist([ACO_lat tx_lat(i)],[ACO_lon tx_lon(i)])/1000;
    end

    
    
end