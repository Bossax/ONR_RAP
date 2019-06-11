close all
clear
%%
aco_loc = [22.7388,-158.006186]; % updated June 2017
[aco_x,aco_y,utmzone] = deg2utm(aco_loc(1),aco_loc(2));
waypoint = [];
waypoint(1,:) = aco_loc;
%% Circular 25 km
R = 25; % km
Re = 6371; %km
% UTM easting northing
E = aco_x+R*1000;
W = aco_x-R*1000;
X = linspace(W,E,100); %geo
Y1 = sqrt((R*1000)^2-(X-aco_x).^2)+ aco_y;
Y2 = -sqrt((R*1000)^2-(X-aco_x).^2)+ aco_y;
% Degree
utm = cell(length(X),1);
utm(:) = {utmzone};
utm = cell2mat(utm);
[lat,lon] = utm2deg(X',Y1',utm);
Y1 = lat;
[lat,lon] = utm2deg(X',Y2',utm);
Y2 = lat;
X = lon;
X = [X flip(X)];
Y = [Y1 Y2];
ind = [];
for i = 1:length(Y)
   if ~isreal(Y(i))
       ind(end+1) = i;
   end
end
Y(ind) = [];
X(ind) = [];

% ACO Location
f = figure(1);
f.Units = 'normalized';
f.Position = [0.2 0.2 0.7 0.8];
aco = scatter(aco_loc(2),aco_loc(1),200,'kx','LineWidth',5);
text(aco_loc(2)+0.0042,aco_loc(1)-0.022,'ACO','FontSize',13);
hold on

% plot r 25
r25 = plot(X,real(Y),'Color',[229 83 .0]/255);
grid on

% ACO label
% lat = (aco_loc(1));
% lon = aco_loc(2);
% label = sprintf("lat = %.5f \nlon = %.5f",lat,lon);
% text(lon+0.03,lat-0.02,label);


% label 1
x = aco_x;
y = aco_y - R*1000;
[lat,lon] = utm2deg(x,y,utm(1,:));
waypoint(end+1,:) = [lat lon];
text(lon,lat+.005,'25km','FontSize',14,'FontWeight','bold');

%%
% 15 km
R = 15; % km
Re = 6371; %km
% UTM easting northing
E = aco_x+R*1000;
W = aco_x-R*1000;
X = linspace(W,E,100); %geo
Y1 = sqrt((R*1000)^2-(X-aco_x).^2)+ aco_y;
Y2 = -sqrt((R*1000)^2-(X-aco_x).^2)+ aco_y;
% Degree
utm = cell(length(X),1);
utm(:) = {utmzone};
utm = cell2mat(utm);
[lat,lon] = utm2deg(X',Y1',utm);
Y1 = lat;
[lat,lon] = utm2deg(X',Y2',utm);
Y2 = lat;
X = lon;
X = [X flip(X)];
Y = [Y1 Y2];

% plot r 15
r15 = plot(X,real(Y),'g');
grid on

% label 2
x = aco_x;
y = aco_y - R*1000;
[lat,lon] = utm2deg(x,y,utm(1,:));
waypoint(end+1,:) = [lat lon];
text(lon+0.005,lat-0.005,'15km','FontSize',14,'FontWeight','bold');



%5km
R = 5; % km
Re = 6371; %km
% UTM easting northing
E = aco_x+R*1000;
W = aco_x-R*1000;
X = linspace(W,E,100); %geo
Y1 = sqrt((R*1000)^2-(X-aco_x).^2)+ aco_y;
Y2 = -sqrt((R*1000)^2-(X-aco_x).^2)+ aco_y;
% Degree
utm = cell(length(X),1);
utm(:) = {utmzone};
utm = cell2mat(utm);
[lat,lon] = utm2deg(X',Y1',utm);
Y1 = lat;
[lat,lon] = utm2deg(X',Y2',utm);
Y2 = lat;
X = lon;
X = [X flip(X)];
Y = [Y1 Y2];
hold on
% plot r 5
r5 = plot(X,real(Y),'r');
grid on

% label 3
x = aco_x;
y = aco_y - R*1000;
[lat,lon] = utm2deg(x,y,utm(1,:));
waypoint(end+1,:) = [lat lon];
text(lon+0.005,lat-0.005,'5km','FontSize',14,'FontWeight','bold');


% 10km
R = 10; % km
Re = 6371; %km
% UTM easting northing
E = aco_x+R*1000;
W = aco_x-R*1000;
X = linspace(W,E,100); %geo
Y1 = sqrt((R*1000)^2-(X-aco_x).^2)+ aco_y;
Y2 = -sqrt((R*1000)^2-(X-aco_x).^2)+ aco_y;
% Degree
utm = cell(length(X),1);
utm(:) = {utmzone};
utm = cell2mat(utm);
[lat,lon] = utm2deg(X',Y1',utm);
Y1 = lat;
[lat,lon] = utm2deg(X',Y2',utm);
Y2 = lat;
X = lon;
X = [X flip(X)];
Y = [Y1 Y2];
hold on
% plot r 5
r10 = plot(X,real(Y),'m');
grid on

% label 3
x = aco_x;
y = aco_y - R*1000;
[lat,lon] = utm2deg(x,y,utm(1,:));
waypoint(end+1,:) = [lat lon];
text(lon+0.005,lat-0.005,'10km','FontSize',14,'FontWeight','bold');



%{
%% Radial angled
R = 17.67;
N =  aco_y+R*1000;
S =  aco_y-R*1000;
E =  aco_x+R*1000;
W =  aco_x-R*1000;
X2 = linspace(W,E,100);
Y2 = linspace(S,N,100);
% Degree
utm = cell(length(X),1);
utm(:) = {utmzone};
utm = cell2mat(utm);
[Y2,X2] = utm2deg(X2',Y2',utm);

hold on
radial1 = plot(X2,Y2,'b--');
hold on
plot(flip(X2),Y2,'b--')

% label 4
x = W;
y = N;
[lat,lon] = utm2deg(x,y,utm(1,:));
waypoint(end+1,:) = [lat lon];
text(lon,lat,'4','FontSize',14,'FontWeight','bold');

% label 5
x = E;
y = S;
[lat,lon] = utm2deg(x,y,utm(1,:));
waypoint(end+1,:) = [lat lon];
text(lon,lat,'5','FontSize',14,'FontWeight','bold');

% label 6
x = W;
y = S;
[lat,lon] = utm2deg(x,y,utm(1,:));
waypoint(end+1,:) = [lat lon];
text(lon,lat,'6','FontSize',14,'FontWeight','bold');

% label 7
x = E;
y = N;
[lat,lon] = utm2deg(x,y,utm(1,:));
waypoint(end+1,:) = [lat lon];
text(lon,lat,'7','FontSize',14,'FontWeight','bold');
%}
%% Radial 
N =  aco_y+25*1000;
S =  aco_y-25*1000;
E =  aco_x+25*1000;
W =  aco_x-25*1000;
X3 = aco_loc(2).*ones(1,50);
[N_lat,lon] = utm2deg(aco_x,N,utm(1,:));
[S_lat,lon] = utm2deg(aco_x,S,utm(1,:));
Y3 = linspace(S_lat,N_lat,50);
hold on
radial2 = plot(X3,Y3,'c');
%%
Y4= aco_loc(1).*ones(1,50);
[lat,W_lon] = utm2deg(W,aco_y,utm(1,:));
[lat,E_lon] = utm2deg(E,aco_y,utm(1,:));
X4 = linspace(W_lon,E_lon,50);
hold on
plot(X4,Y4,'c')
%{
% label 8
x = E;
y = aco_y;
[lat,lon] = utm2deg(x,y,utm(1,:));
waypoint(end+1,:) = [lat lon];
text(lon,lat,'8','FontSize',14,'FontWeight','bold');

% label 9
x = W;
y = aco_y;
[lat,lon] = utm2deg(x,y,utm(1,:));
waypoint(end+1,:) = [lat lon];
text(lon,lat,'9','FontSize',14,'FontWeight','bold');

%% label 10
x = aco_x;
y = N;
[lat,lon] = utm2deg(x,y,utm(1,:));
waypoint(end+1,:) = [lat lon];
text(lon,lat,'10','FontSize',14,'FontWeight','bold');

%% label 11
x = aco_x;
y = S;
[lat,lon] = utm2deg(x,y,utm(1,:));
waypoint(end+1,:) = [lat lon];
text(lon-0.01,lat+0.02,'11','FontSize',14,'FontWeight','bold');

%% Spin
S1 = [aco_x,aco_y-15000];
S2 = [aco_x,aco_y];
S3 = [aco_x+2000,aco_y];
S4 = [aco_x+4000,aco_y];
S5 = [aco_x+6000,aco_y];
S = [S1;S2 ;S3 ;S4 ;S5];
% Degree
utm = cell(length(S),1);
utm(:) = {utmzone};
utm = cell2mat(utm);
[S_lat,S_lon] = utm2deg(S(:,1),S(:,2),utm);

hold on
spin = scatter(S_lon,S_lat,'ro');

% label 12
x = S3(1);
y = S3(2);
[lat,lon] = utm2deg(x,y,utm(1,:));
waypoint(end+1,:) = [lat lon];
text(lon,lat,'12','FontSize',14,'FontWeight','bold');

% label 13
x = S4(1);
y = S4(2);
[lat,lon] = utm2deg(x,y,utm(1,:));
waypoint(end+1,:) = [lat lon];
text(lon,lat,'13','FontSize',14,'FontWeight','bold');

% label 14
x = S5(1);
y = S5(2);
[lat,lon] = utm2deg(x,y,utm(1,:));
waypoint(end+1,:) = [lat lon];
text(lon,lat,'14','FontSize',14,'FontWeight','bold');
%}

%% Square
R = 10;
N= aco_y+R*1000; 
S= aco_y-R*1000;
E = aco_x+R*1000;
W = aco_x-R*1000;
% 1
X = linspace(W,E,100); 
Y = S*ones(1,100);
% Degree
utm = cell(length(X),1);
utm(:) = {utmzone};
utm = cell2mat(utm);
[lat,lon] = utm2deg(X,Y,utm);
hold on
sq = plot(lon,lat,'k');

% 2
X = linspace(W,E,100); 
Y = N*ones(1,100);
% Degree
utm = cell(length(X),1);
utm(:) = {utmzone};
utm = cell2mat(utm);
[lat,lon] = utm2deg(X,Y,utm);
hold on
plot(lon,lat,'k')

%3
X = W*ones(1,100);
Y = linspace(S,N,100);
% Degree
utm = cell(length(X),1);
utm(:) = {utmzone};
utm = cell2mat(utm);
[lat,lon] = utm2deg(X,Y,utm);
hold on
plot(lon,lat,'k')

%4
X = E*ones(1,100);
Y = linspace(S,N,100);
% Degree
utm = cell(length(X),1);
utm(:) = {utmzone};
utm = cell2mat(utm);
[lat,lon] = utm2deg(X,Y,utm);
hold on
plot(lon,lat,'k')

%{
% label 15
x = E;
y = aco_y;
[lat,lon] = utm2deg(x,y,utm(1,:));
waypoint(end+1,:) = [lat lon];
text(lon,lat,'15','FontSize',14,'FontWeight','bold');


% label 16
x = E;
y = S;
[lat,lon] = utm2deg(x,y,utm(1,:));
waypoint(end+1,:) = [lat lon];
text(lon,lat,'16','FontSize',14,'FontWeight','bold');


% label 17
x = W;
y = S;
[lat,lon] = utm2deg(x,y,utm(1,:));
waypoint(end+1,:) = [lat lon];
text(lon,lat,'17','FontSize',14,'FontWeight','bold');

% label 18
x = W;
y = N;
[lat,lon] = utm2deg(x,y,utm(1,:));
waypoint(end+1,:) = [lat lon];
text(lon,lat,'18','FontSize',14,'FontWeight','bold');

% label 19
x = E;
y = N;
[lat,lon] = utm2deg(x,y,utm(1,:));
waypoint(end+1,:) = [lat lon];
text(lon,lat,'19','FontSize',14,'FontWeight','bold');

legend([aco r25 r15(1) r5(1) radial1 radial2 spin sq],{'ACO','25km','15km','5km','Angled Radial','Radial','Spin','10km Square'})

%}
%% lawn mow
U = [aco_x,aco_y].*ones(6,2) + 5000*[2:-1:-3 ; 3.*ones(1,6)]';
L = [aco_x,aco_y].*ones(7,2) + 5000*[3:-1:-3 ; -3.*ones(1,7)]';
p = [aco_x+15000 ,aco_y];
utm = cell(length(U),1);
utm(:) = {utmzone};
utm = cell2mat(utm);
[U_lat,U_lon] =  utm2deg(U(:,1),U(:,2),utm);

utm = cell(length(L),1);
utm(:) = {utmzone};
utm = cell2mat(utm);
[L_lat,L_lon] =  utm2deg(L(:,1),L(:,2),utm);

[p_lat,p_lon] = utm2deg(p(1),p(2),utm(1,:));

for i = 1:2:length(U)
    hold on
   plot([U_lon(i+1) U_lon(i)],[mean(U_lat) mean(U_lat)],'b');
    
end

for i = 1:2:length(L)-1
    hold on
   plot([L_lon(i+1) L_lon(i)],[mean(L_lat) mean(L_lat)],'b');
    
end

for i = 1:6
    hold on
    lon = (L_lon(i+1)+U_lon(i))/2;
   plot([lon lon],[L_lat(i+1) U_lat(i)],'b');
    
end
hold on
plot([p_lon L_lon(1)],[p_lat L_lat(1)],'b');

%% Spin
S1 = [aco_x,aco_y-5000];
S2 = [aco_x+10000,aco_y];
S3 = [aco_x-10000,aco_y-10000];
S4 = [aco_x+10000,aco_y-10000];
S5 = [aco_x,aco_y+10000];
S6 = [aco_x+15000,aco_y];
S = [S1;S2 ;S3 ;S4;S5;S6];

% Degree
utm = cell(length(S),1);
utm(:) = {utmzone};
utm = cell2mat(utm);
[S_lat,S_lon] = utm2deg(S(:,1),S(:,2),utm);

hold on
spin = scatter(S_lon,S_lat,'r','filled');

title('RAP Tomography Ship Paths')

legend([spin],{'Spin'})

 axis tight
xlabel('Lon')
ylabel('Lat')
set(gca,'fontsize',13)
