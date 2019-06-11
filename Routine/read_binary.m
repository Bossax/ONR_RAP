%% Read Binary data from Sensor 2 (Group 103)
clearvars
close all

% direct to POSMV binary folder
cd ../..
cd POS-MV

fname='pos-bin2018102911';          %EDIT !!!

num_bytes=0;
i=1;
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

elseif B(i,1)==105
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
fclose(fid);

%%
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
x_vel_b=G(keep_rows2,1);           %m/s
y_vel_b=G(keep_rows2,2);           %m/s
z_vel_b=G(keep_rows2,3);           %m/s
roll_b=H(keep_rows2,1);            %degrees (-180,180)
pitch_b=H(keep_rows2,2);           %degrees (-90,90)
heading_b=H(keep_rows2,3);         %degrees(0,360)
wander_angle_b=H(keep_rows2,4);    %degrees (-180,180)
heave_b=-I(keep_rows2,1);           %meters
rate_lon_b=I(keep_rows2,2);        %deg/s
rate_tran_b=I(keep_rows2,3);       %deg/s
rate_down_b=I(keep_rows2,4);       %deg/s
accel_lon_b=I(keep_rows2,5);       %m/s^2
accel_tran_b=I(keep_rows2,6);      %m/s^2
accel_down_b=I(keep_rows2,7);      %m/s^2


% Performance Matrics
keep_rows1_1=find(B(:,1)==105);
t1_1=C(keep_rows1_1,1);                   %seconds of the week
lat_err=F1(keep_rows1_1,1);             %degrees
lon_err=F1(keep_rows1_1,2);             %degrees
altitude_err=F1(keep_rows1_1,3);        %meters
x_vel_err=G1(keep_rows1_1,1);           %m/s
y_vel_err=G1(keep_rows1_1,2);           %m/s
z_vel_err=G1(keep_rows1_1,3);           %m/s
roll_err=H1(keep_rows1_1,1);            %degrees (-180,180)
pitch_err=H1(keep_rows1_1,2);           %degrees (-90,90)
heading_err=H1(keep_rows1_1,3);         %degrees(0,360)


%Convert time to datenum       %%%%%  (NEEDS TO BE EDITED)  %%%%%
t=datenum(2018,10,28)+t1./(3600*24);
t_err=datenum(2018,10,28)+t1_1./(3600*24);



%Save variables in .mat file
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

posmvb.altitude=altitude_b;
posmvb.lat=lat_b;
posmvb.lon=lon_b;
posmvb.x_vel=x_vel_b;
posmvb.y_vel=y_vel_b;
posmvb.z_vel=z_vel_b;
posmvb.roll=roll_b;
posmvb.pitch=pitch_b;
posmvb.heading=heading_b;
posmvb.heave=heave_b;
posmvb.t=t;


cd ../..
cd Tx_Rx_Output/posmv/
save("posmv"+string(fname(end-3:end)),'posmv')
save("posmvb"+string(fname(end-3:end)),'posmvb')

%% Save in excel
%{
Data = [];
Data(:,1) = altitude;
Data(:,2) = lat;
Data(:,3) = lon;
Data(:,4) = x_vel;
Data(:,5) = y_vel;
Data(:,6) = z_vel;
Data(:,7) = roll;
Data(:,8) = pitch;
Data(:,9) = heading;
Data(:,10) = heave;
Data(:,11) = altitude_err;
Data(:,12) = lat_err;
Data(:,13) = lon_err;
Data(:,14) =  x_vel_err;
Data(:,15) = y_vel_err;
Data(:,16) = z_vel_err;
Data(:,17) = roll_err;
Data(:,18) = pitch_err;
Data(:,19) = heading_err;
dlmwrite('km1809-46.csv',Data,'precision', '%f20');
%string(datestr(t)),string(datestr(t_err)));
%}




