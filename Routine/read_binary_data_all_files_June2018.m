%% Read Binary data from Sensor 2 (Group 103) 
% Read "all" POS-MV binary files encoded with transducer position/granite block data and performance matrics
% save the output files in .mat extension with associated labels

clearvars
close all
bin_file_dir = '/Volumes/ACO_RAP_2/RAP/June2018Cruise/POS-MV/after';
save_file_dir = '/Volumes/ACO_RAP_2/RAP/June2018Cruise/Tx_Rx_Output/posmv/after';
cd(bin_file_dir)
%fname = 'km1809.004';
directory=dir('km1809*');

for j =1:length(directory)
fname = directory(j).name;

label = string(fname(end-9:end));


num_bytes=0;
i=1;
cd(bin_file_dir)
fid=fopen(fname,'r');

%for i=1:length(start_pos)
while feof(fid)<1                       %check if the file ends     0 = don't end yet

%Group Start
try
A(i,:)=fread(fid,4,'char','l')';          %$GRP
catch
    break
end
num_bytes=num_bytes+4;

%Group ID and byte count
B(i,:)=fread(fid,2,'ushort','l')';        %Group # and byte size (2 bytes)
num_bytes=num_bytes+4;
num_bytes=num_bytes+B(i,2);

if (B(i,1)~=103 && B(i,1)~=102 && B(i,1)~=105) || (B(i,1)==105 && B(i,2)~=68)
    %remove_data=fread(fid,start_pos(i+1)-(num_bytes+1),'char','l');
    remove_data=fread(fid,B(i,2),'char','l');
else

%Time/Distance Fields
C(i,:)=fread(fid,3,'double','l')';        %Time 1, Time 2, and distance tag (8 bytes)
D(i,:)=fread(fid,2,'ubit4','l')';         %Time types for 1&2 (4 bits)
E(i,:)=fread(fid,1,'uint8','l');         %Distance type (1 byte)

if B(i,1)==103 || B(i,1)==102
    

    %Latitude/Longitude/Altitude
    F(i,:)=fread(fid,3,'double','l')';        %Lat, Lon, and Altitude (8 bytes)
    
    %Velocity
    G(i,:)=fread(fid,3,'float','l')';         %Along, Across, and Down Velocity (4 bytes)
    
    %Roll, Pitch, Heading, Wander Angle
    H(i,:)=fread(fid,4,'double','l')';        %Roll, Pitch, Heading, and Angle (8 bytes)
    
    %Heave, Angular rate about long, Angular rate about transverse, Angular
    %rate about down, Long accel, Transverse accel, Down accel  (4 bytes)
    I(i,:)=fread(fid,7,'float','l')';   
    
    %Padding (2 bytes)
    J(i,:)=fread(fid,1,'uint16','l');
    
    %Checksum (2 bytes)
    K(i,:)=fread(fid,1,'ushort','l');
    
    %Group End (read characters)
    L(i,:)=fread(fid,2,'char','l')';

elseif B(i,1)==105                      % group 105 reports errors
    %Latitude/Longitude/Altitude
    F1(i,:)=fread(fid,3,'float','l')';        %Lat, Lon, and Altitude
    
    %Velocity
    G1(i,:)=fread(fid,3,'float','l')';         %Along, Across, and Down Velocity
    
    %Roll, Pitch, Heading
    H1(i,:)=fread(fid,3,'float','l')';        %Roll, Pitch, and Heading
    
    %Padding
    J1(i,:)=fread(fid,1,'uint16','l');
    
    %Checksum
    K1(i,:)=fread(fid,1,'ushort','l');
    
    %Group End
    L1(i,:)=fread(fid,2,'char','l')';
end

i=i+1;
end

end

fclose(fid)
fid = 0
keep_rows1=find(B(:,1)==103);

%Extract data from vectors for transducer
t1=C(keep_rows1,1);              %seconds of the week
t2=C(keep_rows1,2);              %seconds since pos-mv powered on
dist_tag=C(keep_rows1,3);        %meters
lat=F(keep_rows1,1);             %degrees
lon=F(keep_rows1,2);             %degrees
altitude=F(keep_rows1,3);        %meters
x_vel=G(keep_rows1,1);           %m/s
y_vel=G(keep_rows1,2);           %m/s
z_vel=G(keep_rows1,3);           %m/s
roll=H(keep_rows1,1);            %degrees (-180,180)
pitch=H(keep_rows1,2);           %degrees (-90,90)
heading=H(keep_rows1,3);         %degrees(0,360)
wander_angle=H(keep_rows1,4);    %degrees (-180,180)
heave=-I(keep_rows1,1);           %meters
rate_lon=I(keep_rows1,2);        %deg/s
rate_tran=I(keep_rows1,3);       %deg/s
rate_down=I(keep_rows1,4);       %deg/s
accel_lon=I(keep_rows1,5);       %m/s^2
accel_tran=I(keep_rows1,6);      %m/s^2
accel_down=I(keep_rows1,7);      %m/s^2


keep_rows2=find(B(:,1)==102);
%Extract data for granite block
lat_b=F(keep_rows2,1);             %degrees
lon_b=F(keep_rows2,2);             %degrees
altitude_b=F(keep_rows2,3);        %meters


% Performance
keep_rows1_1=find(B(:,1)==105);
t1_1=C(keep_rows1_1,1);                   %seconds of the week
lat_err=F1(keep_rows1_1,1);             %meter
lon_err=F1(keep_rows1_1,2);             %meter
altitude_err=F1(keep_rows1_1,3);        %meter
x_vel_err=G1(keep_rows1_1,1);           %m/s
y_vel_err=G1(keep_rows1_1,2);           %m/s
z_vel_err=G1(keep_rows1_1,3);           %m/s (down)
roll_err=H1(keep_rows1_1,1);            %degrees (-180,180)
pitch_err=H1(keep_rows1_1,2);           %degrees (-90,90)
heading_err=H1(keep_rows1_1,3);         %degrees(0,360)

% check sunday (Oct)
%{
date =char(label);
date = str2num(string(date(7:8)));
sunday = 21
if date >=28
   sunday = 28
end
%}
sunday = 17;        % June

%Convert time to datenum       %%%%%  (NEEDS TO BE EDITED)  %%%%%
t=datenum(2018,6,sunday)+t1./(3600*24);
t_err=datenum(2018,6,sunday)+t1_1./(3600*24);



%% Save variables in .mat file

posmv.altitude=altitude;
posmv.lat=lat;
posmv.lon=lon;
posmv.x_vel=x_vel;
posmv.y_vel=y_vel;
posmv.z_vel=z_vel;
posmv.roll=roll;
posmv.pitch=pitch;
posmv.heading=heading;
posmv.heave=heave;
posmv.t=t;
posmv.altitude_err=altitude_err;
posmv.lat_err=lat_err;
posmv.lon_err=lon_err;
posmv.x_vel_err=x_vel_err;
posmv.y_vel_err=y_vel_err;
posmv.z_vel_err=z_vel_err;
posmv.roll_err=roll_err;
posmv.pitch_err=pitch_err;
posmv.heading_err=heading_err;
posmv.t_err=t_err;
posmv.rate_lon = rate_lon;
posmv.rate_tran = rate_tran;
posmv.rate_down = rate_down;
posmv.accel_lon=accel_lon;
posmv.accel_tran=accel_tran;
posmv.accel_down=accel_down;


cd(save_file_dir)
hour = datestr(t(1),'HH');
day = datestr(t(1),'dd');

save_fname = "posmv_201806"+day+hour;
fprintf(save_fname +"\n.\n.\n")
save(save_fname+".mat",'posmv')

posmv_gb.t = t;
posmv_gb.lat = lat_b;
posmv_gb.lon = lon_b;
posmv_gb.altitude = altitude_b;

cd ../granite_block/all/
save("gb_"+save_fname+".mat",'posmv_gb');
end

%{
%% Save in delimiter-seprated-valued excel file
% separate date and time
for i=1:length(t)
    r = datestr(t(i));
    time(i) = string(r(end-7:end));
    d(i) = string(r(1:end-8));
end


filename = "km1809-"+label+".csv";
variable_name = {'Date','Time','Altitude','Latitude' ,'Longitude','x_velocity' ,'y_velocity' ,'z_velocity','Roll','Pitch','Heading','Heave','Altitude_Error','Latitude_Error','Longitude_Error'};
Package= table(d',time',altitude,lat,lon,x_vel,y_vel,z_vel,roll,pitch,heading,heave,altitude_err,lat_err,lon_err,'VariableNames',variable_name);%,'VariableNames',varname);
writetable(Package,filename);

%}

cd /Users/testuser/Documents/MATLAB/Script/RTX/

RTX_POSMVGGA_savefile


