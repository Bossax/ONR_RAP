% Verfify Rotation procedure
% function [pos_ant_lat,pos_ant_lon,pos_ant_altitude] = ant_pos_posmv(posmv_lat,posmv_lon,posmv_altitude,posmv_t,posmv_roll,posmv_pitch,posmv_yaw,sensor)
% 1. Calculate the antenna position by using the transmit transducer position provied by posmv binary and the attitude data 
% 2. find the nearest posmv data point to the rtx antenna data point and use the posmv attitude data from that point of time
% 3. report the transducer position in lat and lon
% 4. Column vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input section
posmv_latgn = 22.2;
posmv_longn = -158.6;
posmv_altitudegn = 8;
posmv_t = 1;
posmv_roll = 0;
posmv_pitch = 0;
posmv_yaw = 0;


% Transformation 
% coordinates of devices with respect to the granite block
% ant = [6.439;6.505;-27.761]; % from KM coordinate system May 2018 spreasheet
% tx5 = [0.568;19.599;0.715];
tx5 = [1;2;3];
ant = [4;4;4];
gn = [0;0;0];

ant2tx = tx5-ant;    
ant2gn = -ant;
lever_arm = ant2tx;
tx_gb = [];                % vector for antenna position in the global UTM coordinate


% reference point
% convert lat lon of gn to northing and easting (y,x) UTM
for k =1:length(posmv_longn)
    [ant_x,ant_y,utmzone] = deg2utm(posmv_latgn(k),posmv_longn(k));
    tx_gb(k,:) = [ant_x,ant_y]; % easting, northing
end
tx_gb = horzcat(tx_gb,posmv_altitudegn)';

% Create the transducer vector
td_gb = tx_gb +tx5;
% Visualiztion
figure(1)
clf
% UTM
line([0 10],[0 0],[0 0],'Color','r','Linewidth',1)
hold on
line([0 0],[0 10],[0 0],'Color','b','Linewidth',1)
line([0 0],[0 0],[0 10],'Color','g','Linewidth',1)
scatter3(10,0,0,40,'r','filled')
scatter3(0,10,0,40,'b','filled')
scatter3(00,0,10,40,'g','filled')
% plot3([tx_gb(1)-tx_gb(1) td_gb(1)-tx_gb(1)],[tx_gb(2)-tx_gb(2) td_gb(2)-tx_gb(2)],[tx_gb(3)-tx_gb(3) td_gb(3)-tx_gb(3)],'--r')
% scatter3(td_gb(1)-tx_gb(1),td_gb(2)-tx_gb(2),td_gb(3)-tx_gb(3),40,'r','filled')
grid on
xlim([-10 10])
ylim([-10 10])
zlim([-10 10])
xlabel('X')
ylabel('Y')
zlabel('Z')
view([1,1,1])
% rotation matrix: Global, UTM coordiante to unrotated ship coordinate
thetax_gb = 180/180*pi;
thetaz_gb = -90/180*pi;
Rxgb = [1 0 0;0 cos(thetax_gb) -sin(thetax_gb);0 sin(thetax_gb) cos(thetax_gb)];
Rzgb = [cos(thetaz_gb) -sin(thetaz_gb) 0 ; sin(thetaz_gb) cos(thetaz_gb) 0; 0 0 1];
Rgb_s = Rxgb*Rzgb;  % X then Z

% Rotation matrix of the ship attitudes
% retreive data from the POS-MV
utmx = [10;0;0];
utmy = [0;10;0];
utmz = [0;0;10];
roll_posmv = posmv_roll/180*pi;
pitch_posmv = posmv_pitch/180*pi;
yaw_posmv = posmv_yaw/180*pi;


% ship unrotated axes
ship1x = Rgb_s*utmx;
ship1y = Rgb_s*utmy;
ship1z = Rgb_s*utmz;

% update axes
pause
figure(1)
clf
line([0 ship1x(1)],[0 ship1x(2)],[0 ship1x(3)],'Color','r','Linewidth',2)
hold on
line([0 ship1y(1)],[0 ship1y(2)],[0 ship1y(3)],'Color','b','Linewidth',2)
line([0 ship1z(1)],[0 ship1z(2)],[0 ship1z(3)],'Color','g','Linewidth',2)
scatter3(ship1x(1),ship1x(2),ship1x(3),40,'r','filled')
scatter3(ship1y(1),ship1y(2),ship1y(3),40,'b','filled')
scatter3(ship1z(1),ship1z(2),ship1z(3),40,'g','filled')
grid on
xlim([-10 10])
ylim([-10 10])
zlim([-10 10])
xlabel('X')
ylabel('Y')
zlabel('Z')
view([1,1,1])
title('Unrotated ship coordinate frame')



% Calculation to UTM
ant_gb = zeros(3,length(roll_posmv));   % for antenna position

for k = 1:length(roll_posmv)
    Rz = [cos(yaw_posmv(k)) -sin(yaw_posmv(k)) 0; sin(yaw_posmv(k)) cos(yaw_posmv(k)) 0; 0 0 1];    % yaw
    Ry = [cos(pitch_posmv(k)) 0 sin(pitch_posmv(k)); 0 1 0; -sin(pitch_posmv(k)) 0 cos(pitch_posmv(k))]; %pitch
    Rx = [1 0 0; 0 cos(roll_posmv(k)) -sin(roll_posmv(k)); 0 sin(roll_posmv(k)) cos(roll_posmv(k))]; % roll
    R = Rz*Ry*Rx;
    ant_gb(:,k) = tx_gb(:,k)-Rgb_s*R*lever_arm;
end
ship3ax =  Rxgb*Rzgb*Rz*Ry*Rx*[utmx utmy utmz];
ship3x = ship3ax(:,1);
ship3y = ship3ax(:,2);
ship3z = ship3ax(:,3);
% update axes
pause
figure(1)
clf
line([0 ship3x(1)],[0 ship3x(2)],[0 ship3x(3)],'Color','r','Linewidth',3)
hold on
line([0 ship3y(1)],[0 ship3y(2)],[0 ship3y(3)],'Color','b','Linewidth',3)
line([0 ship3z(1)],[0 ship3z(2)],[0 ship3z(3)],'Color','g','Linewidth',3)
scatter3(ship3x(1),ship3x(2),ship3x(3),40,'r','filled')
scatter3(ship3y(1),ship3y(2),ship3y(3),40,'b','filled')
scatter3(ship3z(1),ship3z(2),ship3z(3),40,'g','filled')
grid on
xlim([-10 10])
ylim([-10 10])
zlim([-10 10])
xlabel('X')
ylabel('Y')
zlabel('Z')
view([1,1,1])
title('Rotated ship coordinate frame')

% % convert back to lat lon
% col = size(ant_gb);
% col = col(2);
% utm = cell(col,1);
% utm(:) = {utmzone};
% utm = cell2mat(utm);
% [pos_ant_lat,pos_ant_lon] = utm2deg(ant_gb(1,:)',ant_gb(2,:)',utm);
% pos_ant_altitude = ant_gb(3,:)';


