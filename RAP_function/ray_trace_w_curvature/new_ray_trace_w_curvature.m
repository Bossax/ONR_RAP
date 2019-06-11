
clear
close all
%%%%%%%%%
td_depth = 3.454;
icListen_lat=22.73911569;
icListen_lon=-158.006106601;                
x_dist = 25123;

%%%%%%%%%
fid=fopen('h306a0202.ctd');                       % October 2018

D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f','headerlines',6));
fclose(fid);
        
z_hyd =  4736.266;

pres=D(:,1);        %pressure (dbars)
temp=D(:,2);        %temperature
sal=D(:,3);         %salinity
sal_A = gsw_SA_from_SP(sal,pres,icListen_lon,icListen_lat);

%A priori SS field in Cartesian coordinate
SS1 = gsw_sound_speed(sal,temp,pres);
SS2 = gsw_sound_speed(sal_A,temp,pres);
z = gsw_z_from_p(-pres,icListen_lat);


%%
% truncate sound speed profile
% upper end
rm_ind = find(td_depth >= z,1,'last');
z(1:rm_ind) = [];
SS1(1:rm_ind) = [];
% lower end
add_z = z_hyd  - z(end);
layer_thickness = z(end)-z(end-1);
layer_num = floor(add_z/layer_thickness);
z = vertcat(z,z(end)+transpose([1:layer_num].*layer_thickness));
z(end+1) = z_hyd;

p_end = gsw_p_from_z(icListen_lat,z_hyd);
%%
% %Correct for Earth Curvature
Re=6371000;        % epsilonray tracing APL code
Rc = Re;
% coordinate transformation
epsilon = z./Re;
% z_r =z.*(1-epsilon./2+epsilon.*epsilon/6);
% SS_r = SS.*(1-z_r./Re);
z_r = z.*(1+epsilon/2+(epsilon.^2)/3);
SS_r = SS1.*(1+epsilon+epsilon.^2);

