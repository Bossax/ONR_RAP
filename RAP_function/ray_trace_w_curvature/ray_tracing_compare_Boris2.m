%% ray trace calculation

clear 
close all
% Sound speed profile
cd /Users/testuser/Documents/MATLAB/Script/ray_trace_w_curvature
load('SS_profile_Oct.mat')
plot(SS,z)
title('Sound Speed Profile')
xlabel('Sound Speed (m/s)')
ylabel('depth')
set(gca,'YDir','reverse')
grid on
figure(2)
clf
%% Genrate surface distances
dist_now = 10;
z_offset0 = -6;
iii = 1;
while dist_now <= 23000
    dist_now = dist_now+ rand*100+200;
    x_dist(iii) = dist_now;
    z_offset(iii)  = 0;
%     z_offset(iii) = z_offset0+(2*rand-7);  
     iii = iii+1;
   
end
%%
%%%%
%%% constant sound speed in each layer
theta_all = [2:0.01:89]./(180).*pi;  % critical grazing angle
n_layer_ant = length(z)-1;
horz_range = [];
r_inc =[];
for k = 1:length(theta_all)
    theta_layer_n = [];
    r = [];
    t = [];
    a = cos(theta_all(k))/SS(1);    % ray parameter at the top-most layer

    % grazing angle of each layer
    theta_layer_n(1) = theta_all(k);
    for j =2:n_layer_ant
        theta_layer_n(j) = acos(a*SS(j)); 
    end
    r = (z(2)-z(1))./tan(theta_layer_n);    % horizontal distance in all layers
    % incremental horizontal distance array
    for l =1:length(r);
        r_inc(l) = sum(r(1:l));
    end
%     if mod(k,20) == 0
%         hold on
%         plot(flip(r_inc)/1000,z(2:end))
%     end
    % horizontal range 
    horz_range1(k) = sum(r);
    
end
%%
keep_ind =[];
horz_dist = [];
for iii = 1:length(x_dist)
    keep_ind(iii) = find(horz_range1<x_dist(iii),1,'first');
end
theta_kept = [theta_all(keep_ind)' theta_all(keep_ind-1)'];
horz_range1_kept = [horz_range1(keep_ind)' horz_range1(keep_ind-1)'];
%%
max_x_err = 0.001;
theta_graz_launch = [];
est_tt = [];
horz_dist = [];
arc_length = [];
for n = 1:length(x_dist)
    n
    theta_new = -(theta_kept(n,1)-theta_kept(n,2))*((x_dist(n) - horz_range1_kept(n,1))/(horz_range1_kept(n,2)-horz_range1_kept(n,1)))+theta_kept(n,1);
    theta_high = theta_kept(n,2);
%     theta_high/pi*180
    theta_low = theta_new;
    x_high = horz_range1_kept(n,2);
    x_err = 10;
    while abs(x_err) > max_x_err
        theta_layer_n = [];
        r = [];
        t = [];
        a = cos(theta_new)/SS(1);    % ray parameter at the top-most layer
        
        % grazing angle of each layer
        theta_layer_n(1) = theta_all(k);
        for j =2:n_layer_ant
            theta_layer_n(j) = acos(a*SS(j));
        end
        r = (z(2)-z(1))./tan(theta_layer_n);    % horizontal distance in all layers
        % incremental horizontal distance array
        for l =1:length(r);
            r_inc(l) = sum(r(1:l));
        end
        if mod(k,20) == 0
            hold on
            plot(flip(r_inc)/1000,z(2:end))
        end
        % horizontal range
        horz_range_hold = sum(r);
        % sound speed gradient of each layer
        
        % ray arc length of each layer and elasped time
        for j = 1:n_layer_ant
            s(j) = -(z(j+1)-z(j))/sin(theta_layer_n(j));
            t(j) = s(j)/SS(j);
        end
        
        % tot ray arc length
        ray_arc_length_hold = sum(s);
        
        % tot time
        elapsed_time_hold = sum(t);
        
        x_low  = horz_range_hold;
        x_err = x_dist(n) -x_low;
        
        del_theta = (theta_low-theta_high)
        if del_theta/pi*180 <= 0.01
            del_theta = 0.01/180*pi
        end
        
        theta_new = theta_low-x_err/(x_high - x_low)*del_theta;
        
        
        theta_low = theta_new;
        pause
     
    end
    theta_graz_launch(n) = theta_new;
    horz_dist(n) = horz_range_hold;
    est_tt(n) = elapsed_time_hold;
    arc_length(n) = ray_arc_length_hold;
    
end
%%  plot ray loop at launch angle of 11.89
figure(2)
r_inc = [r_inc flip(2*r_inc(end)-r_inc(1:end))];
z_plot = [z(1:end-1) flip(z(1:end-1))];
plot(r_inc/1000,-z_plot)
set(gca,'YDir','reverse')

grid on
xlabel('Range (km)')
ylabel('Depth (m)')
titlename = sprintf('Bottom Limiting Ray \n Layer ThickneSS = 1m');
title(titlename)
set(gca,'XTick',[0:2:38])
set(gca,'YTick',[0:200:2000])
set(gca,'XAxisLocation','top')
xlim([0 38])
fprintf('Layer ThickneSS = %d m\n',z(1)-z(2))  
fprintf('Constant Sound Speed \n')  
fprintf('Ray Arc Length = %f \n',ray_arc_length1*2)
fprintf('Horizontal Length = %f\n',horz_range1*2)
fprintf('Elapsed Time = %f\n',elapsed_time1*2)
%% constant sound speed gradient 
% theta_all = pi/2 - 11.89.*pi/180;   % critical incident angle
figure(3)
for k = 1:length(theta_all)
    theta_layer_n = [];
    r = [];
    t = [];
    b = [];
   
    a = sin(theta_all(k))/SS(1);    % ray parameter at the top-most layer

    % incident angle of each layer
    theta_layer_n(1) = theta_all(k);
    for j =2:length(SS)
        theta_layer_n(j) = asin(a*SS(j)); 
    end
     
    % sound speed gradient
    for j =1:n_layer_ant
       b(j) = -(SS(j+1)-SS(j))/(z(j+1)-z(j));    % SS profile slope
    end
    
%     % horizontal range and elapsed time
%     for j = 1:n_layer_ant
%        trans_coeff = SS(j)/b(j) - z(j); % translation coefficient
%        w1 = z(j) + trans_coeff; % upper bound
%        w2 = z(j+1) + trans_coeff;   % lower bound
%        r(j) = abs(1/(a*b(j)))*(sqrt(1-a^2*b(j)^2*w1^2)-sqrt(1-a^2*b(j)^2*w2^2)); 
%        t(j) = abs(1/b(j))*log((w2*(1+sqrt(1-a^2*b(j)^2*w1^2)))/(w1*(1+sqrt(1-a^2*b(j)^2*w2^2)))); 
%     end
    for j= 1:n_layer_ant
        v = SS(j+1)/SS(j);
        h = (sqrt(1-(a*SS(j+1))^2)+1)/(sqrt(1-(a*SS(j))^2)+1);
        t(j) = (1/b(j))*(log(v)-log(h));
        r(j) = (1/(a*b(j)))*(sqrt(1-(a*SS(j))^2) - sqrt(1-(a*SS(j+1))^2));
        arc_distance(j)= 1/abs(a*b(j))*abs(theta_layer_n(j+1) - theta_layer_n(j));
    end   
       



    % incremental horizontal distance array
    for l =1:length(r);
        r_inc(l) = sum(r(1:l));
    end
    if mod(k,20) == 0
        hold on
        plot(flip(r_inc)/1000,z(2:end))
    end
    
%     % ray arc length of each layer and elasped time
%     for j = 1:n_layer_ant
%         s(j) = -(z(j+1)-z(j))/cos(theta_layer_n(j));
%     end
%     
    % tot ray arc length
    ray_arc_length(k) = sum(arc_distance);
    % total horizontal range 
    horz_range(k) = sum(r);
    % elapsed time
    elapsed_time(k) = sum(t);
end
%%
    
fprintf('-------------------------- \n')   
fprintf('Slowly Changing Sound Speed \n')  
fprintf('Ray Arc Length = %f \n',ray_arc_length*2)
fprintf('Horizontal Length = %f\n',horz_range*2)
fprintf('Elapsed Time = %f\n',elapsed_time*2)
fprintf('-------------------------- \n')  
fprintf('Arc Length Difference (SS constant - SS Gradient = %f\n',ray_arc_length1 -ray_arc_length)
fprintf('Travel Time Difference (SS constant - SS Gradient = %f\n',elapsed_time1 -elapsed_time)

