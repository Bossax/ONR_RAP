function [pos_ant_lat,pos_ant_lon,pos_ant_altitude] = ant_pos_posmv(posmv_lat,posmv_lon,posmv_altitude,posmv_t,posmv_roll,posmv_pitch,posmv_yaw,sensor)
% 1. Calculate the antenna position using the POSMV binary data
% 2. find the nearest posmv data point to the corresponding data point and use the posmv attitude data from that point of time
% 3. report the transducer position in lat and lon
% 4. column vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if the inputs are column vectors
[r,c]= size(posmv_lat);
if r == 1
   posmv_lat =posmv_lat'; 
end
[r,c]= size(posmv_lon);
if r == 1
   posmv_lon =posmv_lon'; 
end
[r,c] = size(posmv_altitude);
if r == 1
   posmv_altitude =posmv_altitude'; 
end
[r,c] = size(posmv_t);
if r == 1
   posmv_t =posmv_t'; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transformation 
% coordinates of devices with respect to the granite block
ant = [6.439;6.505;-27.761]; % from KM coordinate system May 2018 spreasheet
tx5 = [0.568;19.599;0.715];
ant2tx = tx5-ant;    
ant2gn = -ant;

if sensor == 1
    lever_arm = ant2gn;
elseif sensor == 2
    lever_arm = ant2tx;
else
    fprintf('Sensor input invalid\n')
    return;
end
lever_arm
tx_gb = [];                % vector for antenna position in the global UTM coordinate
% convert lat lon of tx_pos to northing and easting (y,x) UTM
for k =1:length(posmv_lon)
    [ant_x,ant_y,utmzone1] = deg2utm(posmv_lat(k),posmv_lon(k));
    tx_gb(k,:) = [ant_x,ant_y]; % easting, northing
end
tx_gb = horzcat(tx_gb,posmv_altitude)';

% rotation matrix: Global, UTM coordiante to unrotated ship coordinate
thetax_gb = 180/180*pi;
thetaz_gb = -90/180*pi;
Rxgb = [1 0 0;0 cos(thetax_gb) -sin(thetax_gb);0 sin(thetax_gb) cos(thetax_gb)];
Rzgb = [cos(thetaz_gb) -sin(thetaz_gb) 0 ; sin(thetaz_gb) cos(thetaz_gb) 0; 0 0 1];
Rgb_s = Rxgb*Rzgb;

% Rotation matrix of the ship attitudes
% retreive data from the POS-MV

roll_posmv = posmv_roll/180*pi;
pitch_posmv = posmv_pitch/180*pi;
yaw_posmv = posmv_yaw/180*pi;

% Calculation
ant_gb = zeros(3,length(roll_posmv));   % for antenna position

for k = 1:length(roll_posmv)
    Rz = [cos(yaw_posmv(k)) -sin(yaw_posmv(k)) 0; sin(yaw_posmv(k)) cos(yaw_posmv(k)) 0; 0 0 1];    % yaw
    Ry = [cos(pitch_posmv(k)) 0 sin(pitch_posmv(k)); 0 1 0; -sin(pitch_posmv(k)) 0 cos(pitch_posmv(k))]; %pitch
    Rx = [1 0 0; 0 cos(roll_posmv(k)) -sin(roll_posmv(k)); 0 sin(roll_posmv(k)) cos(roll_posmv(k))]; % roll
    R = Rz*Ry*Rx;
    ant_gb(:,k) = tx_gb(:,k)-Rgb_s*R*lever_arm;
end

% convert back to lat lon
col = size(ant_gb);
col = col(2);
utm = cell(col,1);
utm(:) = {utmzone1};
utm = cell2mat(utm);
[pos_ant_lat,pos_ant_lon] = utm2deg(ant_gb(1,:)',ant_gb(2,:)',utm);
pos_ant_altitude = ant_gb(3,:)';



end