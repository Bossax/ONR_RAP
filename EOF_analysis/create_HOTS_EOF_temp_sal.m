%% Create Empirical Orthogonal Functions for HOTS Data (Temp and Sal)
clearvars
close all

cd Data

d=dir('*.mn');

for i=1:length(d)
    
fid=fopen(d(i).name);

B=fgets(fid);               % read the header section
f_mo=regexp(B,' ');          % specifiy delimiter indices
month(i)=str2double(B(f_mo(3)+1:f_mo(3)+2));    % extract the month of the file

D=cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f%f')); % convert the text to table
fclose(fid);

z_dist=4728;        %Hydrophone depth

pres=D(:,1);        %pressure (dbars)
%pres=gsw_z_from_p(pres,22.5);       %pressure (m)
temp=D(:,2);        %temperature
sal=D(:,3);         %salinity

%keep shallow entries
if length(pres)<500
else

%Remove deeper depths
% remove_points=find(pres>z_dist);
pres(501:end)=[];
temp(501:end)=[];
sal(501:end)=[];

% stack temp and sal in the arrays
temp_all(:,i)=temp;
sal_all(:,i)=sal;
end

end

cd ..

%Remove Zero entries
month(:,sum(abs(temp_all))==0)=[];
sal_all(:,sum(abs(temp_all))==0)=[];
temp_all(:,sum(abs(temp_all))==0)=[];



%Remove SS with bad entries
remove_col=[];
for i=1:size(sal_all,1)
    for j=1:size(sal_all,2)
    
    check_1=find(abs(sal_all(i,j)-sal_all(i,:))>0.25);
    check_2=find(abs(temp_all(i,j)-temp_all(i,:))>2);
    
    if length(check_1)>size(sal_all,2)/2 || length(check_2)>size(temp_all,2)/2
        remove_col(end+1)=j;
    end
    
    end
end
sal_all(:,unique(remove_col))=[];
temp_all(:,unique(remove_col))=[];
month(unique(remove_col))=[];
%%

% mean_SS=mean(SS_all,2);
% 
% for i=1:size(SS_all,2)
% SS_anomaly(:,i)=SS_all(:,i)-mean_SS;
% end


%Create temp/sal anomaly for each month
for i=1:12
    
    mo_pos=find(i==month);      % find columns with the correspinding month
    
    mean_sal(:,i)=mean(sal_all(:,mo_pos),2);
    mean_temp(:,i)=mean(temp_all(:,mo_pos),2);
    
    for ii=1:length(mo_pos)
    sal_anomaly(:,mo_pos(ii))=sal_all(:,mo_pos(ii))-mean_sal(:,i);
    temp_anomaly(:,mo_pos(ii))=temp_all(:,mo_pos(ii))-mean_temp(:,i);
    end
    
end


[EOFs_sal,lambda_sal,contribution_sal,PCs_sal]=EOF_analysis(sal_anomaly);
[EOFs_temp,lambda_temp,contribution_temp,PCs_temp]=EOF_analysis(temp_anomaly);





figure(1)
plot(pres(:,1),EOFs_temp(:,1),'r')
hold on
plot(pres(:,1),EOFs_sal(:,1),'b')
view([90 -90])
set(gca,'Xdir','reverse')
xlabel('Pressure (dbar)')
ylabel('Amplitude')
title(sprintf('Mode 1 - Temperature (%f%%), Salinity (%f%%)',contribution_temp(1,1),contribution_sal(1,1)))
legend('Temperature','Salinity')
figure(2)
plot(pres(:,1),EOFs_temp(:,2),'r')
hold on
plot(pres(:,1),EOFs_sal(:,2),'b')
view([90 -90])
set(gca,'Xdir','reverse')
xlabel('Pressure (dbar)')
ylabel('Amplitude')
title(sprintf('Mode 2 - Temperature (%f%%), Salinity (%f%%)',contribution_temp(2,1),contribution_sal(2,1)))
legend('Temperature','Salinity')


figure(1000)
plot(pres(:,1),temp_all)
view([90 -90])
set(gca,'Xdir','reverse')
xlabel('Pressure (dbar)')
ylabel('Temperature')

figure(1001)
plot(pres(:,1),sal_all)
view([90 -90])
set(gca,'Xdir','reverse')
xlabel('Pressure (dbar)')
ylabel('Salinity')


