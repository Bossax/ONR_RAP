close all
clearvars -except arc1 arc2 t2 t1 z1 SS1

ACO_lat = 22.738772;
ACO_lon = -158.006186;
ACO_depth = -4479.92;                   % original depth MSL local u/e/n frame

x_dist = 21000;
td_h = -4;

%% 1. CTD data (MSL)
% CTD file directory
% cd '/Users/testuser/Documents/MATLAB/Script/Data'
fid=fopen('h306a0202.ctd');                       % 12 October 2018

D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f','headerlines',6));
fclose(fid);

pres_MSL = D(:,1);        %pressure (dbars)
temp = D(:,2);        %temperature
sal = D(:,3);         %salinity (Sp)
sal = gsw_SA_from_SP(sal,pres_MSL,ACO_lon,ACO_lat);

% Depths in the ellipsoidal coordinate; downward positive (-h)
ACO_gh = geoidheight(ACO_lat,ACO_lon+360,'EGM96');
z_td = -td_h;
z_hyd = -(ACO_depth+ACO_gh);      

% disp(sprintf('transducer depth = %f',z_td))
%% 2 truncate CTD data to fit the range between the receiver depth and the source depth
% on the ellipsoidal surface Orthometrc height

N_mean = 2.32;                                      % average geiod height over the area of coverage 
z_CTD = -(gsw_z_from_p(pres_MSL,ACO_lat) );         % relative to MSL
z_elps = z_CTD-N_mean;                              % from MSL to ellipsoid


% 2.1 truncate the upper bound using transducer depth
rm_ind = find(z_elps <= z_td); 
z_elps(rm_ind(1:end-1)) = [];
z_CTD(rm_ind(1:end-1)) = [];
temp(rm_ind(1:end-1)) = [];
sal(rm_ind(1:end-1)) = [];

% set the first depth to be the transducer depth
z_diff = z_td - z_elps(1);
temp(1) = temp(1)+(temp(2)-temp(1))*(z_diff/(z_elps(2)-z_elps(1)));   % interpolate the first temp and salinity
sal(1)=sal(1)+(sal(2)-sal(1))*(z_diff/(z_elps(2)-z_elps(1)));   
z_elps(1) = z_td;
z_CTD(1) = z_elps(1)+N_mean;    % transducer depth + mean geoid undulation = transducer depth relative to MSL


% 2.2 truncate the lower bound
if z_elps(end) > z_hyd
    rm_ind = find(z_elps>z_hyd);
    z_elps(rm_ind) = [];
    z_CTD(rm_ind) = [];
    temp(rm_ind) = [];
    sal(rm_ind) = [];
else
    ly_thk = -mean(z_elps(end-5:end-1) - z_elps(end-4:end));
    add_layer = floor((z_hyd -z_elps(end))/ly_thk);
    z_elps = vertcat(z_elps,z_elps(end)+ly_thk*(1:add_layer)');
    z_elps(end+1) = z_hyd;
    
    z_CTD = vertcat(z_CTD,z_elps(end-add_layer:end)+N_mean);
    
    temp_gradient = temp(end-1)-temp(end);
    temp = vertcat(temp,temp(end)+temp_gradient*(1:add_layer+1)');
    sal_gradient = sal(end-1)-sal(end);
    sal = vertcat(sal,sal(end)+sal_gradient*(1:add_layer+1)');
end

pres_CTD = gsw_p_from_z(-z_CTD,ACO_lat);

%% 3. Creat Sound Speed Profile 
SS_elps = gsw_sound_speed(sal,temp,pres_CTD);

% sound speed profile is accurate enough (Nov27)
%% Earth Flattening Transformation
R_e = 6371000;
z1 = -R_e.*log((R_e-(z_elps))./R_e);
SS1 = SS_elps.*(R_e./(R_e-(z_elps)));

z2 = z_elps;
SS2 = SS_elps;
%%
ax1 = subplot(1,2,1);
plot(ax1,z1-z2,z1);
ax1.YDir = 'reverse';
grid on
axis tight

ax2 = subplot(1,2,2);
plot(ax2,SS1-SS2,z1)
grid on
ax2.YDir = 'reverse';
axis tight

%%
r1 = [];
r2 = [];
%Find theta0 of acoustic transmission
% incident angle
% slant anle
slant_angle = atan(x_dist/z_hyd);
theta0_all= slant_angle-10/180*pi:1*(pi/180):slant_angle+10/180*pi;   % all possible launch angles
theta0_all(find(theta0_all <= 0.0001)) = [];
theta0_all(find(theta0_all >= 89.5/180*pi)) = [];


%%
% loop over all possible thetas 
% calculate a small horizontal distance within an interval betweenf two depths
ang_dis = [];             % layer-wise angular displacement
R0 = R_e - z_td;           % Launch radial distance
R_hyd = R_e - z_hyd;      % Hydrophone radial distance
Ri = R_e - z_elps;

figure(88)
for ii = 1:length(theta0_all)
    
    % 1st incident angle at interface of layer 2
    incident_angle = asin(Ri(1)/Ri(2)*sin(theta0_all(ii)));
    
    % Ray parameter or each launch angle
    a1 = sin(theta0_all(ii))/SS1(1);
    a2 = sin(theta0_all(ii))/SS2(1);
    p = R0*sin(incident_angle)/SS2(1);
    
    % incident angles of each layer
    theta_k1(1) = theta0_all(ii);
    theta_d1(1) = theta0_all(ii);
    for i = 2:length(SS1)-1
        theta_k1(i)=asin(a1*SS1(i));
        theta_d1(i)=asin(a2*SS2(i));
    end    
    theta = asin(p*SS2(2:end)./Ri(2:end));
    theta = vertcat(incident_angle,theta);

    % horizontal distance
    for i = 1:length(SS1) - 1
        r1(i)=(z1(i+1)-z1(i))*tan(theta_k1(i));
    end
    
    for i = 1:length(SS2) - 1
        r2(i)=(z2(i+1)-z2(i))*tan(theta_d1(i));
    end
    
    % spherical coordinates
    for gg = 1:length(SS2)-1
        Rmean = (Ri(gg)+Ri(gg+1))/2;
        ang_dis(gg)= p*(Ri(gg) - Ri(gg+1))/(sqrt((Rmean/SS2(gg))^2-p^2)*Rmean);
%         ang_dis(gg) = p*SS2(gg)*(Ri(gg) - Ri(gg+1))/(Rmean^2*sqrt(1-(p*SS2(gg)/Rmean)^2));
%         r3(gg) = ang_dis(gg)*Rmean;
        fe = (R_e - z_elps(gg))/R_e;
        r3(gg) = (z_elps(gg+1)-z_elps(gg))/(fe*cot(theta(gg))); % APL
    end
    
    %Total x distance
    x_diff1(ii)=sum(r1);  % total horizontal distance for one SS profile launching theta
    x_diff2(ii)=sum(r2);  
    x_diff3(ii)=sum(r3);  
    
    figure(88)
    hold on
    plot(cumsum(r1)/1000,z1(1:end-1),'k')
    plot(cumsum(r2)/1000,z2(1:end-1),'r')
    plot(cumsum(r3)/1000,z2(1:end-1),'g')
    axis tight
    grid on
    set(gca,'YDir','reverse')
    
end
legend('EFT','Planar','Spherical')

%%
%Find min and max launch angle for Newton-Raph method 

t_low1 = find(x_diff1<x_dist);      % find all launch angles which yield x distance within the bound
theta_low1 = theta0_all(t_low1(end));   % label the lowest poSS_fltible launch angle (correspond to the steepest angle)
x_low1 = x_diff1(t_low1(end));           % The shortest ray arc length, lower bound of the sur dist

t_high1 = find(x_diff1>x_dist);
theta_high1 = theta0_all(t_high1(1));   % The upper bound of sur dist (the least steep angle)
x_high1 = x_diff1(t_high1(1));

% Intrapolation the launch angle to find the optimum ray arc length of the direct path
theta_new1 = theta_low1+((x_dist)/((x_low1-x_high1)/(theta_low1-theta_high1)))-((x_low1)/((x_low1-x_high1)/(theta_low1-theta_high1)));


t_low2 = find(x_diff2<x_dist);      % find all launch angles which yield x distance within the bound
theta_low2 = theta0_all(t_low2(end));   % label the lowest poSS_fltible launch angle (correspond to the steepest angle)
x_low2 = x_diff2(t_low2(end));           % The shortest ray arc length, lower bound of the sur dist

t_high2 = find(x_diff2>x_dist);
theta_high2 = theta0_all(t_high2(1));   % The upper bound of sur dist (the least steep angle)
x_high2 = x_diff2(t_high2(1));

% Intrapolation the launch angle to find the optimum ray arc length of the direct path
theta_new2 = theta_low2+((x_dist)/((x_low2-x_high2)/(theta_low2-theta_high2)))-((x_low2)/((x_low2-x_high2)/(theta_low2-theta_high2)));


t_low3 = find(x_diff3<x_dist);      % find all launch angles which yield x distance within the bound
theta_low3 = theta0_all(t_low3(end));   % label the lowest poSS_fltible launch angle (correspond to the steepest angle)
x_low3 = x_diff3(t_low3(end));           % The shortest ray arc length, lower bound of the sur dist

t_high3 = find(x_diff3>x_dist);
theta_high3 = theta0_all(t_high3(1));   % The upper bound of sur dist (the least steep angle)
x_high3 = x_diff3(t_high3(1));

% Intrapolation the launch angle to find the optimum ray arc length of the direct path
theta_new3 = theta_low3+((x_dist)/((x_low3-x_high3)/(theta_low3-theta_high3)))-((x_low3)/((x_low3-x_high3)/(theta_low3-theta_high3)));


figure(88)
scatter(x_high1/1000,z1(end-1),'kx')
scatter(x_low1/1000,z1(end-1),'kx')

scatter(x_high2/1000,z2(end-1),'rx')
scatter(x_low2/1000,z2(end-1),'rx')


scatter(x_high3/1000,z2(end-1),'gx')
scatter(x_low3/1000,z2(end-1),'gx')

disp(sprintf('Range of angles: %f - %f',theta_low1/pi*180,theta_high1/pi*180))
disp(sprintf('delta theta = %f',(theta_new1 -theta_low1)/pi*180) )
%%
% figure(3)
% plot(theta0_all/pi*180,x_diff1/1000)
% hold on
% plot(slant_angle/pi*180,x_diff1(find(slant_angle <= theta0_all,1,'first'))/1000,'ro')
% plot(theta_low1/pi*180,x_diff1(find(theta_low1 <= theta0_all,1,'first'))/1000,'kx')
% plot(theta_high1/pi*180,x_diff1(find(theta_high1 <= theta0_all,1,'first'))/1000,'kx')
% 
% line([theta0_all(1) theta0_all(end)]/pi*180,[x_dist x_dist]/1000,'Color','k')
% grid on
% ylabel('Range (km)')
% xlabel('Angle (degree)')
% axis tight

%%

% We now have a range of launch angles which gives the nearest number of horizontal distance 
 theta_k1 = theta_low1;
x_old1 = x_low1;           % lower bound of the horizontal distance

x_diff1 = 1000000;   % arbitary number

theta_k2 = theta_low2;
x_old2 = x_low2;           % lower bound of the horizontal distance

x_diff2 = 1000000;   % arbitary number

theta_k3 = theta_low3;
x_old3 = x_low3;           % lower bound of the horizontal distance

x_diff3 = 1000000;   % arbitary number


count = 1;

% Loop until x_diff is close to x_dist to 1 cm 
while abs(x_diff1-x_dist) > 0.0001
    
    arc_distance1 = [];
    tt1 = [];
    r1 = [];
    %Ray parameter
    a1=sin(theta_new1)/SS1(1);
    
    % incident angle
    theta1(1) = theta_new1;
    for i=2:length(SS1)-1
        theta1(i)=asin(a1*SS1(i));
    end
    
    
    % horizontal range in each layer
    for i=1:length(SS1)-1
        r1(i)=(z1(i+1)-z1(i))*tan(theta1(i));
    end
    
    %Total x distance
    x_diff1=sum(r1);
    
    %Total distance and time traveled
    for i=1:length(SS1)-1
        arc_distance1(i)=(z1(i+1)-z1(i))/cos(theta1(i));
%         tt(i)=(z_flt(i+1)-z_flt(i))*(1/SS_flt(i)^2)/sqrt((1/SS_flt(i)^2)-a^2);    
        tt1(i) = arc_distance1(i)/SS1(i);
    end
    tot_dist_all1 = sum(arc_distance1);
    
    arc_lengths_hold1 = arc_distance1;
    surface_dist_hold1 = r1;
    ray_angles_hold1 = theta_k1(2:end);
    
    % Newton-Raphson method
    % Interpolate again to find a new launch angle which is closer to that
    % of the eigenray
    
    theta_d1 = theta_new1;
    
    
    x_new1 = x_diff1;
    
    theta_new1 = theta_d1+((x_dist)/((x_new1-x_old1)/(theta_d1-theta_k1)))-((x_new1)/((x_new1-x_old1)/(theta_d1-theta_k1)));
    
    % Set theta_low to theta 2 (new angle)
    theta_k1 = theta_d1;
    % update nearest horizontal distance
    x_old1 = x_new1;
    
    disp(sprintf('delta theta = %f',(theta_new1 -theta_d1)/pi*180) )
    disp(sprintf('difference between calculated range and actual range =  %f m',x_diff1-x_dist))
    
    figure(89)
    hold on
    plot(cumsum(r1)/1000,z1(1:end-1),'k','LineWidth',count*0.5)
    axis tight
    grid on
    set(gca,'YDir','reverse')
    
    if abs(x_diff1-x_dist)<0.0001
        break;
    end
    count = count+1;
    pause
end
count = 1;
while abs(x_diff2-x_dist) > 0.0001
    
    arc_distance2 = [];
    tt2 = [];
    r2 = [];
    %Ray parameter
    a2=sin(theta_new2)/SS2(1);
    
    % incident angle
    theta2(1) = theta_new2;
    for i = 2:length(SS2)-1
        theta2(i)=asin(a2*SS2(i));
    end 
    
    
    % horizontal range in each layer
    for i = 1:length(SS2)-1
        r2(i)=(z2(i+1)-z2(i))*tan(theta2(i));
    end
    
    %Total x distance
    x_diff2 = sum(r2);
    
    %Total distance and time traveled
    for i=1:length(SS2)-1
        arc_distance2(i)=(z2(i+1)-z2(i))/cos(theta2(i));
%         tt(i)=(z_flt(i+1)-z_flt(i))*(1/SS_flt(i)^2)/sqrt((1/SS_flt(i)^2)-a^2);    
        tt2(i) = arc_distance2(i)/SS2(i);
    end
    tot_dist_all2 = sum(arc_distance2);
    
    arc_lengths_hold2 = arc_distance2;
    surface_dist_hold2 = r2;
    ray_angles_hold2 = theta_k2(2:end);
    
    % Newton-Raphson method
    % Interpolate again to find a new launch angle which is closer to that
    % of the eigenray
    
    theta_d2 = theta_new2;
    
    
    x_new2 = x_diff2;
    
    theta_new2 = theta_d2+((x_dist)/((x_new2-x_old2)/(theta_d2-theta_k2)))-((x_new2)/((x_new2-x_old2)/(theta_d2-theta_k2)));
    
    % Set theta_low to theta 2 (new angle)
    theta_k2 = theta_d2;
    % update nearest horizontal distance
    x_old2 = x_new2;
    
    disp(sprintf('delta theta = %f',(theta_new2 -theta_d2)/pi*180) )
    disp(sprintf('difference between calculated range and actual range =  %f m',x_diff2-x_dist))
    
    
    figure(89)
    hold on
    plot(cumsum(r2)/1000,z2(1:end-1),'b','LineWidth',count*0.5)
    axis tight
    grid on
    set(gca,'YDir','reverse')
    
    if abs(x_diff2-x_dist)<0.0001
        break;
    end
    count = count+1;
    
    pause
end

count = 1;
while abs(x_diff3-x_dist) > 0.0001
    
    arc_distance3 = [];
    tt3 = [];
    r3 = [];
    
    incident_angle = asin(Ri(1)/Ri(2)*sin(theta_new3));
    %Ray parameter
    p = R0*sin(incident_angle)/SS2(1);
    % incident angle
    theta3 = asin(p*SS2(2:end)./Ri(2:end));
    theta3 = vertcat(incident_angle,theta3);
    
    
    % horizontal range in each layer
    for gg = 1:length(SS2)-1
        
        
        R_mean = (Ri(gg)+Ri(gg+1))/2;
        %       ang_dis(gg)= a3*(Ri(gg) - Ri(gg+1))/(sqrt((R_mean/SS2(gg))^2-a3^2)*R_mean);
        ang_dis(gg) = tan(theta3(gg))*(Ri(gg) - Ri(gg+1))/R_mean;
        arc_distance3(gg) = (Ri(gg)-Ri(gg+1)*cos(ang_dis(gg)))/cos(theta3(gg));
        %         r3(gg) = ang_dis(gg)*Ri(gg);
        %         tt3(gg) = arc_distance3(gg)/SS2(gg);
        % APL
        fe = (R_e - z_elps(gg))/R_e;        
        r3(gg) = (z_elps(gg+1)-z_elps(gg))/(fe*cot(theta3(gg))); 
        tt3(gg) = r3(gg)*fe*csc(theta3(gg))/SS2(gg);
    end
    
    %Total x distance
    x_diff3 = sum(r3);
   
    tot_dist_all3 = sum(arc_distance3);
    
    arc_lengths_hold3 = arc_distance3;
    surface_dist_hold3 = r3;
    ray_angles_hold3 = theta3(2:end);
    
    % Newton-Raphson method
    % Interpolate again to find a new launch angle which is closer to that
    % of the eigenray
    
    theta_d3 = theta_new3;
    
    
    x_new3 = x_diff3;
    
    theta_new3 = theta_d3+((x_dist)/((x_new3-x_old3)/(theta_d3-theta_k3)))-((x_new3)/((x_new3-x_old3)/(theta_d3-theta_k3)));
    
    % Set theta_low to theta 2 (new angle)
    theta_k3 = theta_d3;
    % update nearest horizontal distance
    x_old3 = x_new3;
    
    disp(sprintf('delta theta = %f',(theta_new3 -theta_d3)/pi*180) )
    disp(sprintf('difference between calculated range and actual range =  %f m',x_diff3-x_dist))
    
    
    figure(89)
    hold on
    plot(cumsum(r3)/1000,z2(1:end-1),'g','LineWidth',count*0.5)
    axis tight
    grid on
    set(gca,'YDir','reverse')
    
    if abs(x_diff3-x_dist)<0.0001
        break;
    end
    count = count+1;
    
    pause
end

%%
figure(90)
% plot(cumsum(arc_distance1)/1000,z1(1:end-1),'k','LineWidth',1)
% hold on
% plot(cumsum(arc_distance2)/1000,z2(1:end-1),'r','LineWidth',1)
% plot(cumsum(arc_distance3)/1000,z2(1:end-1),'--g','LineWidth',1)
% axis tight
% grid on
% set(gca,'YDir','reverse')

plot(cumsum(tt1),z1(1:end-1),'k','LineWidth',1)
hold on
plot(cumsum(tt2),z2(1:end-1),'r','LineWidth',1)
plot(cumsum(tt3),z2(1:end-1),'--g','LineWidth',1)
axis tight
grid on
set(gca,'YDir','reverse')
%% packaging outputs
theta01 = theta_new1;                        % Launch angle
rayarc_dist_tot1 = tot_dist_all1;                   % Ray arc length
arc_lengths_ly1 = arc_lengths_hold1;            % Arc length in each layer
horz_dist_ly1 = surface_dist_hold1;          % horizontal distance in each layer
ray_angles_ly1 = ray_angles_hold1;                % Incident ray angle at each interval

% incremental horizontal distance
inc_r_ly1 = zeros(1,length(horz_dist_ly1)+1);

for i= 2:length(inc_r_ly1)
    inc_r_ly1(i) = inc_r_ly1(i-1)+horz_dist_ly1(i-1);
end

theta02 = theta_new2;                        % Launch angle
rayarc_dist_tot2 = tot_dist_all2;                   % Ray arc length
arc_lengths_ly2 = arc_lengths_hold2;            % Arc length in each layer
horz_dist_ly2 = surface_dist_hold2;          % horizontal distance in each layer
ray_angles_ly2 = ray_angles_hold2;                % Incident ray angle at each interval

% incremental horizontal distance
inc_r_ly2 = zeros(1,length(horz_dist_ly2)+1);

for i= 2:length(inc_r_ly2)
    inc_r_ly2(i) = inc_r_ly2(i-1)+horz_dist_ly2(i-1);
end

theta03 = theta_new3;                        % Launch angle
rayarc_dist_tot3 = tot_dist_all3;                   % Ray arc length
arc_lengths_ly3 = arc_lengths_hold3;            % Arc length in each layer
horz_dist_ly3 = surface_dist_hold3;          % horizontal distance in each layer
ray_angles_ly3 = ray_angles_hold3;                % Incident ray angle at each interval

% incremental horizontal distance
inc_r_ly3 = zeros(1,length(horz_dist_ly3)+1);

for i= 2:length(inc_r_ly3)
    inc_r_ly3(i) = inc_r_ly3(i-1)+horz_dist_ly3(i-1);
end



est_tt1 = sum(tt1);       %Estimated TT
est_tt2 = sum(tt2);       %Estimated TT
est_tt3 = sum(tt3);       %Estimated TT

%Vertically averaged SS_flt for HOTS CTD
SS_flt_HOT_avg2 = mean(SS2);
