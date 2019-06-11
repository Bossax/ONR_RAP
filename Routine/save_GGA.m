% save posmv GGA dat file to matlab structure file
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/RTX/before
A = load('rtx1Hz.dat');
% store variable
year = A(:,1);
mon = A(:,2);
day = A(:,3);
hour = A(:,4);
min = A(:,5);
sec = A(:,6);
lat = A(:,7);
lon = A(:,8);
altitude = A(:,9);
lat_err = A(:,10);
lon_err = A(:,11);
altitude_err = A(:,12);

%% datenum
keep_ind = find(day <= 31);
year =year(keep_ind);
mon =  mon(keep_ind);
day = day(keep_ind);
hour = hour(keep_ind);
min = min(keep_ind);
sec = sec(keep_ind);
lat = lat(keep_ind);
lon = lon(keep_ind);
altitude = altitude(keep_ind);
lat_err = lat_err(keep_ind);
lon_err = lon_err(keep_ind);
altitude_err = altitude_err(keep_ind);

t_format = 'yyyy-mm-dd HH:MM:SS';

for iii = 1:length(sec)
   sec_hd = sec(iii);
   if sec_hd <=9
       sec_hd = "0"+string(sec_hd);
   else
       sec_hd = string(sec_hd);
   end
  
   min_hd = min(iii);
   if min_hd <=9
       min_hd = "0"+string(min_hd);
   else
       min_hd = string(min_hd);
   end
  
   hr_hd = hour(iii);
   if hr_hd <=9
       hr_hd = "0"+string(hr_hd);
   else
       hr_hd = string(hr_hd);
   end
   
   day_hd = string(day(iii));
   mon_hd = string(mon(iii));
   year_hd = string(year(iii));
  
   t_str(iii) = year_hd+"-"+mon_hd+"-"+day_hd+" "+hr_hd+":"+min_hd+":"+sec_hd;
   
end
t_str = t_str';
t = datenum(t_str,t_format);
%% save
RTX.rtx_t = t;
RTX.rtx_lat = lat;
RTX.rtx_lon = lon;
RTX.rtx_altitude = altitude;
RTX.rtx_lat_err = lat_err;
RTX.rtx_lon_er = lon_err;
RTX.rtx_altitude_err = altitude_err;
cd /Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Tx_Rx_Output/posmv/GGA

save('RTX_beforecruise_1Hz','RTX')