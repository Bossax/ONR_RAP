% interpolation of the antenna position in x,y,z directions in the ship coordinate
% transform to the UTM coordinate
% 1. Find RTX data point matches of transmission times
% 2. read in 10 Hz POS MV data
% 3. find 2 nearest 10 Hz POS MV data to extract roll/pitch/yaw/velocity
% 4. Use time differences between the RTX data point and those 2 POS MV
% data points to interpolate for the antenna position
% 5. convert back to lat/lon altitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



















