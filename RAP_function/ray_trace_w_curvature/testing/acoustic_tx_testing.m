% test the ray tracing code by varying lat,lon,altitude,geometry
% use a tx sample file to represent the transducer location
% test hypotheses
% 1. show the ttp sensitivity to the position error with diff ranges
% 2. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear
range = [100,1000,5000,10000];%,15000,25000];
offstep = 0:0.1:10;
est_tt = zeros(1,length(offstep));
figure(1)
z_offset = -6.5;
for jjj = 1:length(range)
    est = [];
    for iii = 1:length(offstep)
        x_dist = range(jjj)+offstep(iii);
        
        [~,tot_dist(iii),theta0(iii),~,~,~,~,est_tt(iii),~] = ray_trace_w_curvature_v3(x_dist,z_offset);
    end
    y = abs(est_tt-(est_tt(1)))*1000;
    plot(offstep,y)
    text(10.2,y(end),num2str(y(end)));
    hold on
    grid on
end
title('TTP error vs position error with range')
ylabel('TTP error (ms)')
xlabel('position error (m)')
legend('0.1 km','1 km','5 km','10 km')
%% read CTD data
z_offset = 0;
[~,~,~,SS,z,~,~,~,~] = ray_trace_w_curvature_v3(x_dist,z_offset);
figure(2)
plot(SS(1:6),z(1:6))
set(gca,'YDir','reverse')
ylabel('depth')
xlabel('SS (m/s)')
title('Sound Speed in shallow layer')
grid on
xlim([1532 1533])
%% travel time of different path
% 1 direct
x_dist =9000;      % 5km
z_offset = 7.5;
[~,tot_dist_direct,theta0_direct,SS,z,~,~,est_tt_direct,~] = ray_trace_w_curvature_v3(x_dist,z_offset);

% 2 hull-reflected 1
del_s = 15;
del_t1 = del_s/SS(1);
x_dist1 = x_dist - del_s;
[~,tot_dist_hrf1,theta0_hrf1,~,~,~,~,est_tt_hrf1,~] = ray_trace_w_curvature_v3(x_dist1,z_offset);
est_tt_hrf1 = est_tt_hrf1+del_t1;


del_s2 = sqrt(del_s^2+6.6^2);
del_t2 = del_s2/SS(1);
x_dist2 = x_dist - del_s2;
[~,tot_dist_hrf2,theta0_hrf2,~,~,~,~,est_tt_hrf2,~] = ray_trace_w_curvature_v3(x_dist2,z_offset);
est_tt_hrf2 = est_tt_hrf2+del_t2;


vertdis = 0;
del_s3 = del_s;
del_t3 = del_s3/SS(1);
z_offset3 = z_offset-vertdis;
theta = asin(vertdis/del_s3);
x_dist3 = x_dist - del_s3*cos(theta);
[~,tot_dist_hrf3,theta0_hrf3,~,~,~,~,est_tt_hrf3,~] = ray_trace_w_curvature_v3(x_dist3,z_offset3);
est_tt_hrf3 = est_tt_hrf3+del_t3;

fprintf('T diff(Transvers hull-reflecetd path1 - Direct path)\n = %.2f ms\n',(est_tt_hrf1-est_tt_direct)*1000)
fprintf('T diff(Transvers hull-reflecetd path2 - Direct path)\n = %.2f ms\n',(est_tt_hrf2-est_tt_direct)*1000)
fprintf('T diff(Transvers hull-reflecetd path3 - Direct path)\n = %.2f ms\n',(est_tt_hrf3-est_tt_direct)*1000)

