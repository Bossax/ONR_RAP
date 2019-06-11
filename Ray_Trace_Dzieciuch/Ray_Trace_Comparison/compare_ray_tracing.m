%% Compare Ray Tracing (Varamo and Dzieciuch)
clearvars
close all
% 
% seafloor_z=4730;
% lat=22.738894;
% 
% z_c=(0:seafloor_z);                          %DEPTH PROFILE (m) 
% SS_c=cssprofile(z_c,lat);

load('SS_profile.mat')
SS_c = SS;
z_c = z;

%% Compare tt's for different distances
z_offset= -6;
x_dist = 15000;
%x_dist = logspace(3,4.1,5);
%x_dist=(100:100:1000);
% x_dist=(1000:1000:25000);
for i=1:length(x_dist)

[theta0_V(i),est_tt_V(i),SS_EC,z_EC]=ray_trace_Varamo_no_RAP(x_dist(i),z_offset,SS_c,z_c);

[est_tt_D(i),theta0_D(i)]=ray_trace_Dzieciuch(z_EC,SS_EC,-z_offset,x_dist(i));

end

%%
figure(1)
scatter(x_dist./1000,(est_tt_V-est_tt_D).*1000)
hold on
grid on
xlabel('Surface Distance (km)')
ylabel('Time Difference (ms)')
title('Varamo vs. Dzieciuch Ray Trace Estimation')
% 
figure(2)
plot(x_dist,est_tt_V)
hold on
plot(x_dist,est_tt_D)
% plot(SS,-z)
% hold on
% plot(SS_EC,-z_EC)
% plot(SS,-z)
% grid on
% ylabel('Depth (m)')
% xlabel('SS (m/s)')
% % 
% figure(3)
% scatter(x_dist./1000,est_tt_D)
% hold on
% grid on
% xlabel('Surface Distance (km)')
% ylabel('Travel Time (s)')
% title('Ray Trace TT Estimation')




