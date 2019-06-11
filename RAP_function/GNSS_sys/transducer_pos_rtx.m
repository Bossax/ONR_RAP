function [rtx_tx_lat,rtx_tx_lon,rtx_tx_altitude] = transducer_pos_rtx(rtx_lat,rtx_lon,rtx_altitude,rtx_t,posmv_roll,posmv_pitch,posmv_yaw,posmv_t)
% 1. Calculate the transmit transducer position using the antenna positiondata provied by RTX and the attitude data from posmv
% 2. find the nearest posmv data point to the rtx antenna data point and use the posmv attitude data from that point of time
% 3. report the transducer position in lat and lon
% 4. Column vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if the inputs are column vectors
[r,c]= size(rtx_lat);
if r == 1
   rtx_lat =rtx_lat'; 
end
[r,c]= size(rtx_lon);
if r == 1
   rtx_lon =rtx_lon'; 
end
[r,c] = size(rtx_altitude);
if r == 1
   rtx_altitude =rtx_altitude'; 
end
[r,c] = size(rtx_t);
if r == 1
   rtx_t =rtx_t'; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transformation 
% coordinates of devices with respect to the granite block
ant = [6.439;6.505;-27.761]; % from KM coordinate system May 2018 spreasheet
tx5 = [0.568;19.599;0.715];
ant2tx = tx5-ant;    

ant_gb = [];                % vector for antenna position in the global UTM coordinate
% convert lat lon of ant_pos to northing and easting (y,x) UTM
for k =1:length(rtx_lon)
    [ant_x,ant_y,utmzone] = deg2utm(rtx_lat(k),rtx_lon(k));
    ant_gb(k,:) = [ant_x,ant_y]; % easting, northing
end
ant_gb = horzcat(ant_gb,rtx_altitude)';

% rotation matrix: Global, UTM coordiante to unrotated ship coordinate
thetax_gb = 180/180*pi;
thetaz_gb = -90/180*pi;
Rxgb = [1 0 0;0 cos(thetax_gb) -sin(thetax_gb);0 sin(thetax_gb) cos(thetax_gb)];
Rzgb = [cos(thetaz_gb) -sin(thetaz_gb) 0 ; sin(thetaz_gb) cos(thetaz_gb) 0; 0 0 1];
Rgb_s = Rxgb*Rzgb;

% Rotation matrix of the ship attitudes
% retreive data from the POS-MV
% pick nearest neighbor based on timestamps

roll_rtx = [];
pitch_rtx = [];
yaw_rtx = [];
for k =1:length(rtx_t)
   [t_diff,I] = min(abs(posmv_t - rtx_t(k)));
   if t_diff <= 0.5
        roll_rtx(end+1) = posmv_roll(I);
        pitch_rtx(end+1) = posmv_pitch(I);
        yaw_rtx(end+1) = posmv_yaw(I);
   end
end
roll_rtx = roll_rtx/180*pi;
pitch_rtx = pitch_rtx/180*pi;
yaw_rtx = yaw_rtx/180*pi;

% Calculation
tx5_gb = zeros(3,length(roll_rtx));   % for transducer position
for k = 1:length(roll_rtx)
    Rz = [cos(yaw_rtx(k)) -sin(yaw_rtx(k)) 0; sin(yaw_rtx(k)) cos(yaw_rtx(k)) 0; 0 0 1];    % yaw
    Ry = [cos(pitch_rtx(k)) 0 sin(pitch_rtx(k)); 0 1 0; -sin(pitch_rtx(k)) 0 cos(pitch_rtx(k))]; %pitch
    Rx = [1 0 0; 0 cos(roll_rtx(k)) -sin(roll_rtx(k)); 0 sin(roll_rtx(k)) cos(roll_rtx(k))]; % roll
    R = Rz*Ry*Rx;
    tx5_gb(:,k) = Rgb_s*R*ant2tx+ant_gb(:,k);
end

% convert back to lat lon
col = size(tx5_gb);
col = col(2);
utm = cell(col,1);
utm(:) = {utmzone};
utm = cell2mat(utm);
[rtx_tx_lat,rtx_tx_lon] = utm2deg(tx5_gb(1,:)',tx5_gb(2,:)',utm);
rtx_tx_altitude = tx5_gb(3,:)';



end