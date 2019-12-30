function [enu_range,enu_depth] = geodetic2enu(lat_sc,lon_sc,h_sc,lat_rv,lon_rv,h_rv)
%%%%%%%% Input %%%%%%%%%%%
% Geodetic position
%%%%%%%%% Output %%%%%%%%
% ECEF and local horizon coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geodetic to ENU Chadwell
% reference ellipsoid
ref_ellipsoid =  referenceEllipsoid('WGS84');
a = ref_ellipsoid.SemimajorAxis;
b = ref_ellipsoid.SemiminorAxis;
ecc = ref_ellipsoid.Eccentricity;
% source ECEF
if lon_sc<0
    lambda_sc =(lon_sc+360)/180*pi;
else
    lambda_sc =(lon_sc)/180*pi;
end

phi_sc = (lat_sc)/180*pi;
N = a./sqrt(1-ecc^2*sin(phi_sc).^2);
X_sc = (N+h_sc).*cos(phi_sc).*cos(lambda_sc);
Y_sc= (N+h_sc).*cos(phi_sc).*sin(lambda_sc);
Z_sc= (N.*(1-ecc^2)+h_sc).*sin(phi_sc);

% receiver ECEF
if lon_sc<0
    lambda_rv =(lon_rv+360)/180*pi;
else
    lambda_rv =(lon_rv)/180*pi;
end
phi_rv= (lat_rv)/180*pi;
N_rv = a./sqrt(1-ecc^2*sin(phi_rv).^2);
X_rv = (N_rv+h_rv).*cos(phi_rv).*cos(lambda_rv);
Y_rv= (N_rv+h_rv).*cos(phi_rv).*sin(lambda_rv);
Z_rv= (N_rv.*(1-ecc^2)+h_rv).*sin(phi_rv);


for iii =1:length(X_sc)
    
    tf_mat = [      -sin(lambda_sc(iii))               cos(lambda_sc(iii))                     0;
        -sin(phi_sc(iii))*cos(lambda_sc(iii))       -sin(phi_sc(iii))*sin(lambda_sc(iii)) cos(phi_sc(iii));
        cos(phi_sc(iii))*cos(lambda_sc(iii))        cos(phi_sc(iii))*sin(lambda_sc(iii))  sin(phi_sc(iii))];
    
    lpos=  tf_mat*[X_rv-X_sc(iii);Y_rv-Y_sc(iii);Z_rv-Z_sc(iii)];
    e(iii) = lpos(1);
    n(iii) = lpos(2);
    u(iii) = lpos(3);
end

enu_range = sqrt(e.^2+n.^2)';
enu_depth = -u';
end
