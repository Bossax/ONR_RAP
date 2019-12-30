close all
clear

ACO_lat = 22.73877;
ACO_lon = -158.00619;

% Read CTD casts and create sound speed profiles
ctd_name ={'h305a0202.ctd','h306a0202.ctd','h307a0202.ctd'};


for k = 1:length(ctd_name)
    fid = fopen(ctd_name{k});
    D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f','headerlines',6));
    fclose(fid);
    
    pres = D(:,1);        %pressure (dbars)
    temp = D(:,2);        %temperature
    sal = D(:,3);         %salinity (Sp)
    sal = gsw_SA_from_SP(sal,pres,ACO_lon,ACO_lat);
    
    % Creat depth from the surface
    depth_profile{k} = -gsw_z_from_p(pres,ACO_lat); 
    
    % Creat Sound Speed Profile
    SS_profile{k} = gsw_sound_speed(sal,temp,pres);
   
end

%% plot SS
f = figure(1);
clf
f.Units = 'normalized';
f.Position = [0 .5 0.5 0.8];
color = {'red','black','blue'};
linewidth = [1,2,1];
subplot(1,2,1)
for k = 1:length(ctd_name)
    hold on
    p = plot(SS_profile{k},depth_profile{k});
    p.Color = color{k};
    p.LineWidth = linewidth(k);
end
set(gca,'YDir','reverse','fontsize',14)
xlabel('Sound Speed (m/s)')
ylabel('Depth (m)')
legend('Sep 2018','Oct 2018','November 2018','Location','southwest')
title('Sound Speed Profiles')
grid on
ylim([10 4700])
mid_depth = depth_profile{2}(end);


subplot(1,2,2)
for k = 1:length(ctd_name)
    hold on
    SS_len = length(SS_profile{k});
    % check depths
    if mid_depth >= depth_profile{k}
        rel_SS = SS_profile{k}-SS_profile{2}(1:SS_len);
    else
        rel_SS = SS_profile{k}(1:SS_len)-SS_profile{2};
    end
    p = plot(rel_SS,depth_profile{k});
    p.Color = color{k};
    p.LineWidth = linewidth(k);
end
depth_mark = [10 100 1000 2000 3000 4000];
set(gca,'YDir','reverse','YScale','log','YTick',depth_mark,'YTickLabel',depth_mark,'fontsize',14)
xlabel('Sound Speed Anomaly (m/s)')
ylabel('Depth (m)')
legend('Sep 2018','Oct 2018','November 2018','Location','southeast')
title('Relative Sound Speed Profiles')
ylim([10 4700])
grid minor

