%% Create Empirical Orthogonal Functions for HOTS Data
clearvars
close all
cd /Users/testuser/Documents/Vincent_RAP/HOT_289_Data/EOF_analysis
cd Data

ACO_lat=22.73911569;                  % March 2019 #3
ACO_lon=-158.006106601;               

d=dir('*.mn');

for iii=1:length(d)
    
fid=fopen(d(iii).name);

B=fgets(fid);
f_mo=regexp(B,' ');
month(iii)=str2double(B(f_mo(3)+1:f_mo(3)+2));

D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f%f'));
fclose(fid);

z_dist=4728;        %Hydrophone depth

pres_MSL=D(:,1);        %pressure (dbars)+pres_MSL(end)
temp=D(:,2);        %temperature
sal=D(:,3);         %salinity
sal = gsw_SA_from_SP(sal,pres_MSL,ACO_lon,ACO_lat);

%Remove shallow entries
if length(pres_MSL)<2405
else

%Remove deeper depths
% remove_points=find(pres>z_dist);
pres_MSL(2406:end)=[];
temp(2406:end)=[];
sal(2406:end)=[];

%Get SS field
SS_all(:,iii)=gsw_sound_speed(sal,temp,pres_MSL);
end

end
% depth rel to ellipsoid
z = -gsw_z_from_p(pres_MSL,22.738772)-2.32;;
E_radius=6371000;                  % epsilonray tracing APL code
% E_radius=6370229;        
% coordinate transformation
% E=z./E_radius;
% z=z.*(1+(E./2)+((E.^2)./3));

cd ..

%Remove Zero entries
month(:,sum(abs(SS_all))==0)=[];
SS_all(:,sum(abs(SS_all))==0)=[];

%Remove SS with bad entries
remove_col=[];
for iii=1:size(SS_all,1)
    for j=1:size(SS_all,2)
    
    check_1=find(abs(SS_all(iii,j)-SS_all(iii,:))>20);
    
    if length(check_1)>size(SS_all,2)/2
        remove_col(end+1)=j;
    end
    
    end
end
SS_all(:,unique(remove_col))=[];
month(unique(remove_col))=[];


% mean_SS=mean(SS_all,2);
% 
% for i=1:size(SS_all,2)
% SS_anomaly(:,i)=SS_all(:,i)-mean_SS;
% end


%Create SS anomaly for each month
for iii=1:12
    
    mo_pos=find(iii==month);
    
    mean_SS(:,iii)=mean(SS_all(:,mo_pos),2);
    
    for ii=1:length(mo_pos)
    SS_anomaly(:,mo_pos(ii))=SS_all(:,mo_pos(ii))-mean_SS(:,iii);
    end
    
end

RMS_SS = rms(SS_anomaly,2);
[EOFs,lambda,contribution,PCs]=EOF_analysis(SS_anomaly);

%% scale EOFs and singular value
mode1 = EOFs(:,1)./max( EOFs(:,1));
mode2 = EOFs(:,2)./max( EOFs(:,2));
mode3 = EOFs(:,3)./max( EOFs(:,3));
mode4 = EOFs(:,4)./max( EOFs(:,4));
singv1 = sqrt(lambda(1))*max( EOFs(:,1));
singv2 = sqrt(lambda(2))*max( EOFs(:,2));
singv3 = sqrt(lambda(3))*max( EOFs(:,3));
singv4 = sqrt(lambda(4))*max( EOFs(:,4));
EOF_SS.mode1 = mode1;
EOF_SS.mode2 = mode2;
EOF_SS.mode3 = mode3;
EOF_SS.mode4 = mode4;
EOF_SS.singv1 = singv1;
EOF_SS.singv2 = singv2;
EOF_SS.singv3 = singv3;
EOF_SS.singv4 = singv4;
EOF_SS.z = z;
save('EOF_SS','EOF_SS')

%% Plot

subplot(1,4,1)
plot(pres_MSL,mode1)
view([90 -90])
set(gca,'Xdir','reverse')
xlabel('Pressure (dbar)')
ylabel('Amplitude')
title(sprintf('Mode 1 (%.3f%%)',contribution(1,1)))
line([0 5000],[0 0],'color','k')
% title(sprintf('Mode 1 (%.3f%%), SS variance = %.3f (m/s)^{2}',contribution(1,1),singv1^2))
grid on
yticks(-1:0.2:1)
set(gca,'fontsize',13)

subplot(1,4,2)
plot(pres_MSL,mode2)
view([90 -90])
set(gca,'Xdir','reverse')
xlabel('Pressure (dbar)')
ylabel('Amplitude')
title(sprintf('Mode 2 (%.3f%%)',contribution(2,1)))
line([0 5000],[0 0],'color','k')
% title(sprintf('Mode 2 (%.3f%%), SS variance = %.3f (m/s)^{2}',contribution(2,1),singv2^2))
grid on
yticks(-1:0.5:1)
set(gca,'fontsize',13)

subplot(1,4,3)
plot(linspace(0,z_dist,length(EOFs)),mode3)
view([90 -90])
set(gca,'Xdir','reverse')
xlabel('Pressure (dbar)')
ylabel('Amplitude')
line([0 5000],[0 0],'color','k')
title(sprintf('Mode 3 (%.3f%%)',contribution(3,1)))
% title(sprintf('Mode 3 (%.3f%%),SS variance = %.3f (m/s)^{2}',contribution(3,1),singv3^2))
grid on
yticks(-1:0.5:1)
set(gca,'fontsize',13)

subplot(1,4,4)
plot(pres_MSL,mode4)
view([90 -90])
set(gca,'Xdir','reverse')
xlabel('Pressure (dbar)')
ylabel('Amplitude')
title(sprintf('Mode 4 (%.3f%%)',contribution(4,1)))
line([0 5000],[0 0],'color','k')
% title(sprintf('Mode 1 (%.3f%%), SS variance = %.3f (m/s)^{2}',contribution(1,1),singv1^2))
grid on
yticks(-1:0.5:1)
set(gca,'fontsize',13)
%%
figure(4)
plot(linspace(0,z_dist,length(EOFs)),mode4)
view([90 -90])
set(gca,'Xdir','reverse')
xlabel('Pressure (dbar)')
ylabel('Amplitude')
title(sprintf('Mode 4 (%.3f%%),SS variance = %.3f (m/s)^{2}',contribution(4,1),singv4^2))
grid on

figure(1000)
plot(linspace(0,z_dist,length(EOFs)),SS_all)
view([90 -90])
set(gca,'Xdir','reverse')
xlabel('Pressure (dbar)')
ylabel('SS (m/s)')
grid on

