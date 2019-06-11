function obs_G =obs_matrix_hyd(tx_lat,tx_lon,tx_altitude,theta_r,hyd_lat,hyd_lon,hyd_depth)
    % return observation matrix for hydrophone position uncertainty    
    azmth = azimuth(hyd_lat,hyd_lon,tx_lat,tx_lon);
    
    % phi
    phi = (360-azmth)-90;
    if phi<0
        phi = phi+360;
    end
    
    % x_dist =dist([ACO_lat tx_lat],[ACO_lon tx_lon]);
    % theta
    cr = 1539.057;
    
    obs_G = [cos(phi/180*pi)*sin(theta_r)/cr sin(phi/180*pi)*sin(theta_r)/cr cos(theta_r)/cr];
    
end