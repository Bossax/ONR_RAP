function [lat_all,lon_all,altitude_all,X_offset,Y_offset,Z_offset,t_offset] = position_interpolation(rtx_lat,rtx_lon,rtx_altitude,rtx_t,posmvbin_xvel,posmvbin_yvel,posmvbin_zvel,posmvbin_roll,posmvbin_pitch,posmvbin_heading,posmvbin_t,target_t)
% Interpolation positions based on nearest neigbor POS-MV data points and ship dynamics
% Input: 1. Lat 2. Lon 3. Altitude 4. rtxtimestamps 
% target time = the time we want to interpolate for
% Output: Updated 1. Lat 2. Lon 3. Altitude 4. timestamps of the position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
len = length(rtx_t);
% interpolation
lat_all = [];
lon_all = [];
altitude_all = [];
X_offset = [];
Y_offset = [];
Z_offset = [];
t_offset = [];
ind_test = [];
% target time length must be equal to rtx time length
% loop over rtx time
for r = 1:len
      
    % 1. time offset of the actually rtx times from the transmission times
    t_offsetd = 3600*24*(rtx_t(r)  - target_t(r));
         
    % find attitudes indices (at the target time)
    ind_lower_posmv_att = find((posmvbin_t - target_t(r))<=0); % no equal?
    ind_upper_posmv_att = find((posmvbin_t - target_t(r)) >= 0); % no equal?
    
    t_lower = posmvbin_t(ind_lower_posmv_att);
    t_upper = posmvbin_t(ind_upper_posmv_att);
    t_offset_flower = (target_t(r) - t_lower)*3600*24;
         
    if isempty(ind_lower_posmv_att)
        ind_lower_posmv_att = 1;
    else
        ind_lower_posmv_att = ind_lower_posmv_att(end);
    end
    
    if isempty(ind_upper_posmv_att)
        ind_upper_posmv_att = length(posmvbin_t);
    else
        ind_upper_posmv_att = ind_upper_posmv_att(1);
    end
    
    
    % find roll pitch yaw
    lower_roll =  posmvbin_roll(ind_lower_posmv_att);
    lower_pitch =  posmvbin_pitch(ind_lower_posmv_att);
    lower_yaw =  posmvbin_heading(ind_lower_posmv_att);
    upper_roll =  posmvbin_roll(ind_upper_posmv_att);
    upper_pitch =  posmvbin_pitch(ind_upper_posmv_att);
    upper_yaw =  posmvbin_heading(ind_upper_posmv_att);
        
        %instantenous attitudes
        if ind_upper_posmv_att == ind_lower_posmv_att
            roll_now = (upper_roll );
            pitch_now = (upper_pitch );
            yaw_now = (upper_yaw);
            
        else
            roll_now = (upper_roll - lower_roll)*(t_offset_flower/((t_upper - t_lower)*3600*24))+lower_roll;
            pitch_now = (upper_pitch - lower_pitch)*(t_offset_flower/((t_upper - t_lower)*3600*24))+lower_pitch;
            yaw_now = (upper_yaw - lower_yaw)*(t_offset_flower/((t_upper - t_lower)*3600*24))+lower_yaw;
        
        end
        % extract x/y/z vel (of the posmv bin points sandwiching rtx_t and the target t)
        if t_offsetd > 0
            
            ind_upper2 = (find((posmvbin_t - rtx_t(r)) >= 0));
            if isempty(ind_upper2)
                ind_upper2 = length(posmvbin_t);
            else
                ind_upper2 = min(ind_upper2);
            end
            
            ind_lower2 = max(find((posmvbin_t - target_t(r)) <= 0));
            
        elseif t_offsetd <= 0 
            ind_upper2 = min(find((posmvbin_t - target_t(r)) >= 0));
            ind_lower2 = (find((posmvbin_t - rtx_t(r)) <= 0));
            if isempty(ind_lower2)
                ind_lower2 = 1;
            else
                ind_lower2 = max(ind_lower2 );
            end
        end
        
        t_upper2 =  posmvbin_t(ind_upper2);
        t_lower2 = posmvbin_t(ind_lower2);
        
        xvel_upper2 = posmvbin_xvel(ind_upper2);
        xvel_lower2 = posmvbin_xvel(ind_lower2);
        
        yvel_upper2 = posmvbin_yvel(ind_upper2);
        yvel_lower2 = posmvbin_yvel(ind_lower2);
        
        zvel_upper2 = posmvbin_zvel(ind_upper2);
        zvel_lower2 = posmvbin_zvel(ind_lower2);
        
        mean_xvel =  mean([xvel_upper2 xvel_lower2]);
        mean_yvel =  mean([yvel_upper2 yvel_lower2]);
        mean_zvel = mean([zvel_upper2 zvel_lower2]);
        
        % x velocity
        % yaw_tot = horzcat(yaw_tot,yaw_now);
        
        Vx = mean_xvel*sin(yaw_now/180*pi)+mean_yvel*sin(yaw_now/180*pi+pi/2);
        Vy = mean_xvel*cos(yaw_now/180*pi)+mean_yvel*cos(yaw_now/180*pi+pi/2);
        X_offsetd = -Vx*(t_offsetd);      % meter
        Y_offsetd = -Vy*(t_offsetd);      % meter
        Z_offsetd = mean_zvel*(t_offsetd);  % down +
        
        % to UTM
        [td_x,td_y,utmzone] =  deg2utm(rtx_lat(r),rtx_lon(r));
        
        % translate the position to that of the target time
        td_x2 = td_x+X_offsetd;
        td_y2 = td_y+Y_offsetd;
        td_z2 = rtx_altitude(r)+Z_offsetd;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Outputs
        [lat_now,lon_now] = utm2deg(td_x2,td_y2,utmzone);
        
        lat_all = vertcat(lat_all,lat_now);
        lon_all = vertcat(lon_all,lon_now);
        altitude_all = vertcat(altitude_all,td_z2);
        X_offset = vertcat(X_offset,X_offsetd);
        Y_offset = vertcat(Y_offset,Y_offsetd);
        Z_offset = vertcat(Z_offset,Z_offsetd);
        t_offset = vertcat(t_offset,t_offsetd);
        
     
end 
end