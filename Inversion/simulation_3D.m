% 3-D  simulation
% simulate ocean sound speed perturbation field (m) with Gaussian model uncertainty
% Impose empirical orthogonal functions (EOF) on the vertical structure of ocean sound speed variation
% Exclude hydrophone position offsets

close all
clearvars
% Work on the spherical surface
%%%%%%%%%%%%%%%%
% 1. Ocean domain
ACO_lat= 22.738772;                  % June 2017
ACO_lon= -158.006186;                % June 2017
ACO_depth = -4736.266+2.32;         % at 4,736.266 m June 2017


% Initialize the domain
L = 60000;      % meter EDIT
grid_num = 25;  % grid_num x grid_num  pixels EDIT
tot_grid_num = grid_num^2;

[~,lat_u,~] = m_fdist(ACO_lon,ACO_lat,0,L);
[lon_u,~,~] = m_fdist(ACO_lon,ACO_lat,90,L);
[~,lat_l,~] = m_fdist(ACO_lon,ACO_lat,180,L);
[lon_l,~,~] = m_fdist(ACO_lon,ACO_lat,270,L);
n_wall = grid_num+1;
lon_l = lon_l-360;
lon_u = lon_u-360;
x_step = (lon_u - lon_l)/grid_num;
y_step = (lat_u - lat_l)/grid_num;

% x/y boundary coordinates
x_node = lon_l:x_step:lon_u;
y_node = flip(lat_l:y_step:lat_u);

% find center points of the pixels
x_cen = zeros(1,grid_num);
y_cen = zeros(1,grid_num);
for ii = 1:grid_num
    x_cen(ii)= (x_node(ii)+x_node(ii+1))/2;
    y_cen(ii)= (y_node(ii)+y_node(ii+1))/2; 
end

% Domain 
[X,Y] = meshgrid(x_cen,y_cen);
%%%% pixels do not have equal sides !! 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Generate Gaussian Covariance Matrix 
[P_mode1,var_c1] = gaussian_cov_mx(x_cen,y_cen,1); % priori uncerrtainty of mode 1
[P_mode2,var_c2] = gaussian_cov_mx(x_cen,y_cen,2); % priori uncerrtainty of mode 1
[P_mode3,var_c3] = gaussian_cov_mx(x_cen,y_cen,3); % priori uncerrtainty of mode 1
[P_mode4,var_c4] = gaussian_cov_mx(x_cen,y_cen,4); % priori uncerrtainty of mode 1

% bundle up P matrices
dimP = size(P_mode1);
% P = P_mode1;
P = [P_mode1 zeros(dimP) zeros(dimP) zeros(dimP); zeros(dimP) P_mode2 zeros(dimP) zeros(dimP);zeros(dimP) zeros(dimP) P_mode3 zeros(dimP); zeros(dimP) zeros(dimP) zeros(dimP) P_mode4];
 
%% 2.1 For SS simulation
[P_mode1_sim,var_c1] = gaussian_cov_sim(x_cen,y_cen,1); % priori uncerrtainty of mode 1
[P_mode2_sim,var_c2] = gaussian_cov_sim(x_cen,y_cen,2); % priori uncerrtainty of mode 1
[P_mode3_sim,var_c3] = gaussian_cov_sim(x_cen,y_cen,3); % priori uncerrtainty of mode 1
[P_mode4_sim,var_c4] = gaussian_cov_sim(x_cen,y_cen,4); % priori uncerrtainty of mode 1

% bundle up P matrices
dimP = size(P_mode1_sim);
% P = P_mode1;
P_sim = [P_mode1_sim zeros(dimP) zeros(dimP) zeros(dimP); zeros(dimP) P_mode2_sim zeros(dimP) zeros(dimP);zeros(dimP) zeros(dimP) P_mode3_sim zeros(dimP); zeros(dimP) zeros(dimP) zeros(dimP) P_mode4_sim];
%% Priori cov matrix
figure(1)
set(gcf,'name','Covariance Matrix','Units','normalized','Position',[0.3 .5 0.35 0.45])
imagesc(P)
axis tight
grid off
title('Gaussian Covariance Matrix')
colorbar
caxis([0 14])
colormap jet

%% 3 Simulate SS perturbation field
%Cholesky Decomposition to create a simulated SS field
% P is not eactly symmetric, need to asjust by nearestSPD
P1=P_sim;
h =1;
clearvars D
while exist('D')==0
try
D=chol(P1,'lower');
end
P1=nearestSPD(P1);
end
rng(13)
rand_sst=normrnd(0,1,length(P_sim),1);      %Random sound speed perturbation vector with normal distribution
m0 = D*rand_sst;                                %Simulated SS field for each pixel
m01 = m0(1:tot_grid_num);
m02 = m0(tot_grid_num+1:2*tot_grid_num);
m03 = m0(2*tot_grid_num+1:3*tot_grid_num);
m04 = m0(3*tot_grid_num+1:4*tot_grid_num);
m0 = [m01;m02;m03;m04];
%% delta function at the center
mid_p = (tot_grid_num-1)/2;
m01 = [zeros(mid_p,1);0;zeros(mid_p,1)];
m02 = [zeros(mid_p,1);0;zeros(mid_p,1)];
m03 = [zeros(mid_p,1);0;zeros(mid_p,1)];
m04 = [zeros(mid_p,1);0;zeros(mid_p,1)];
m0 = [m01;m02;m03;m04];

%% Gaussian function at the center
x_distance = [];

for iii = 1:grid_num % y
    for jjj = 1:grid_num % x
 [x_distance(end+1),~] = distance(y_cen(11),x_cen(15),y_cen(iii),x_cen(jjj),referenceEllipsoid('WGS84'));
    end
end
G_function = exp(-(x_distance.^2)/(20000*20000));
m0 = [4*G_function';3*G_function';2*G_function';1*G_function'];
m01 = m0(1:tot_grid_num);
m02 = m0(tot_grid_num+1:2*tot_grid_num);
m03 = m0(2*tot_grid_num+1:3*tot_grid_num);
m04 = m0(3*tot_grid_num+1:4*tot_grid_num);

%% on the circle
pixel_num = 260;
m01 = [zeros(pixel_num,1);4;zeros(tot_grid_num-pixel_num-1,1)];
m02 = [zeros(pixel_num,1);3;zeros(tot_grid_num-pixel_num-1,1)];
m03 = [zeros(pixel_num,1);2;zeros(tot_grid_num-pixel_num-1,1)];
m04 = [zeros(pixel_num,1);1;zeros(tot_grid_num-pixel_num-1,1)];
m0 = [m01;m02;m03;m04];

%% in the middle point of a radial line
pixel_num = 315;
m01 = [zeros(pixel_num,1);3;zeros(tot_grid_num-pixel_num-1,1)];
m02 = [zeros(pixel_num,1);2;zeros(tot_grid_num-pixel_num-1,1)];
m03 = [zeros(pixel_num,1);1;zeros(tot_grid_num-pixel_num-1,1)];
m04 = [zeros(pixel_num,1);.5;zeros(tot_grid_num-pixel_num-1,1)];
m0 = [m01;m02;m03;m04];
%%
sim_SS1=reshape(m01,grid_num,grid_num)'  ;% reshape works on a column basis need to swap row and column 
sim_SS2=reshape(m02,grid_num,grid_num)'  ;% reshape works on a column basis need to swap row and column 
sim_SS3 = reshape(m03,grid_num,grid_num)'  ;% reshape works on a column basis need to swap row and column 
sim_SS4 = reshape(m04,grid_num,grid_num)'  ;% reshape works on a column basis need to swap row and column 

% 4 Plot the Sound Speed Field
% draw circle
R = 26000;
num_point = 2000;
xpoint = linspace(lon_l,lon_u,num_point);
inc_ang =360/num_point;
x_cir = zeros(1,num_point);
y_cir = zeros(1,num_point);
for ii = 1:num_point
    angle = inc_ang*ii;
    [x_cir(ii),y_cir(ii),~]=m_fdist(ACO_lon,ACO_lat,angle,R); 
end
x_cir = x_cir-360;

f2 = figure(2);
clf
set(gcf,'name','2-D sound speed perturbation field simulation','Units','normalized','Position',[0 1 0.6 .75])
subplot(2,2,1)
imagesc(x_cen,y_cen,sim_SS1)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title(sprintf('1st mode SS perturbation field'))
% title(sprintf('1st mode SS perturbation field \n(ss variance = %.3f (m/s)^{2})',var_c1))
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
set(gca,'YDir','normal','fontsize',15)
xlabel('Long')
ylabel('Lat')
plot(x_cir,y_cir,'k')
caxis_lim_priori = cbar.Limits;
caxis([-4 4])

subplot(2,2,2)
imagesc(x_cen,y_cen,sim_SS2)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title(sprintf('2nd mode SS perturbation field'));
% title(sprintf('2nd mode SS perturbation field \n(ss variance = %.3f (m/s)^{2})',var_c2))
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
set(gca,'YDir','normal','fontsize',15)
xlabel('Long')
ylabel('Lat')
plot(x_cir,y_cir,'k')
caxis_lim_priori2 = cbar.Limits;
caxis([-4 4])

subplot(2,2,3)
imagesc(x_cen,y_cen,sim_SS3)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title(sprintf('3rd mode SS perturbation field'))
% title(sprintf('3rd mode SS perturbation field \n(ss variance = %.3f (m/s)^{2})',var_c3))
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
set(gca,'YDir','normal','fontsize',15)
xlabel('Long')
ylabel('Lat')
plot(x_cir,y_cir,'k')
caxis_lim_priori3 = cbar.Limits;
caxis([-4 4])

subplot(2,2,4)
imagesc(x_cen,y_cen,sim_SS4)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title(sprintf('4th mode SS perturbation field'))
% title(sprintf('3rd mode SS perturbation field \n(ss variance = %.3f (m/s)^{2})',var_c3))
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
set(gca,'YDir','normal','fontsize',15)
xlabel('Long')
ylabel('Lat')
plot(x_cir,y_cir,'k')
caxis_lim_priori3 = cbar.Limits;
caxis([-4 4])
%% 5 Create Ray Paths
% generate tx points
tx_lon = [];
tx_lat = [];
tx_ind = 1;
%{

% first circle
L =25000;
angle = 0;
angle_all = [];
n1= 300;     % EDIT
figure(2)

subplot(1,3,1)
for iii = 1:n1
    angle = 360/n1*(iii-1);
    [tx_lon(iii),tx_lat(iii),~] = m_fdist(ACO_lon,ACO_lat,angle,L);
    tx_lon(iii) = tx_lon(iii)-360;
    hold on
%     scatter(tx_lon(iii),tx_lat(iii),'po')
     plot([tx_lon(iii) ACO_lon],[tx_lat(iii) ACO_lat],'Linewidth',1,'Color','m')
    
end

subplot(1,3,2)
for iii = 1:n1
    angle = 360/n1*(iii-1);
    [tx_lon(iii),tx_lat(iii),~] = m_fdist(ACO_lon,ACO_lat,angle,L);
    tx_lon(iii) = tx_lon(iii)-360;
    hold on
%     scatter(tx_lon(iii),tx_lat(iii),'po')
     plot([tx_lon(iii) ACO_lon],[tx_lat(iii) ACO_lat],'Linewidth',1,'Color','m')
    
end

subplot(1,3,3)
for iii = 1:n1
    angle = 360/n1*(iii-1);
    [tx_lon(iii),tx_lat(iii),~] = m_fdist(ACO_lon,ACO_lat,angle,L);
    tx_lon(iii) = tx_lon(iii)-360;
    hold on
%     scatter(tx_lon(iii),tx_lat(iii),'po')
     plot([tx_lon(iii) ACO_lon],[tx_lat(iii) ACO_lat],'Linewidth',1,'Color','m')
    
end

L =15000;
angle = 0;
angle_all = [];
n2 = 200;

subplot(1,3,1)
% second circle
for iii = n1+1:n2+n1
    angle = 360/n2*(iii-1);
    [tx_lon(iii),tx_lat(iii),~] = m_fdist(ACO_lon,ACO_lat,angle,L);
    tx_lon(iii) = tx_lon(iii)-360;
    hold on
    %scatter(tx_lon(iii),tx_lat(iii),'po')
     plot([tx_lon(iii) ACO_lon],[tx_lat(iii) ACO_lat],'Linewidth',1,'Color','k')
    
end
subplot(1,3,2)
% second circle
for iii = n1+1:n2+n1
    angle = 360/n2*(iii-1);
    [tx_lon(iii),tx_lat(iii),~] = m_fdist(ACO_lon,ACO_lat,angle,L);
    tx_lon(iii) = tx_lon(iii)-360;
    hold on
    %scatter(tx_lon(iii),tx_lat(iii),'po')
     plot([tx_lon(iii) ACO_lon],[tx_lat(iii) ACO_lat],'Linewidth',1,'Color','k')
    
end

subplot(1,3,3)
% second circle
for iii = n1+1:n2+n1
    angle = 360/n2*(iii-1);
    [tx_lon(iii),tx_lat(iii),~] = m_fdist(ACO_lon,ACO_lat,angle,L);
    tx_lon(iii) = tx_lon(iii)-360;
    hold on
    %scatter(tx_lon(iii),tx_lat(iii),'po')
     plot([tx_lon(iii) ACO_lon],[tx_lat(iii) ACO_lat],'Linewidth',1,'Color','k')
    
end


L =5000;
angle = 0;
angle_all = [];
n3= 100;
subplot(1,3,1)
% 3rd circle
for iii = n1+n2+1:n2+n1+n3
    angle = 360/n3*(iii-1);
    [tx_lon(iii),tx_lat(iii),~] = m_fdist(ACO_lon,ACO_lat,angle,L);
    tx_lon(iii) = tx_lon(iii)-360;
    hold on
    %scatter(tx_lon(iii),tx_lat(iii),'po')
     plot([tx_lon(iii) ACO_lon],[tx_lat(iii) ACO_lat],'Linewidth',1,'Color','b')
    
end
subplot(1,3,2)
% 3rd circle
for iii = n1+n2+1:n2+n1+n3
    angle = 360/n3*(iii-1);
    [tx_lon(iii),tx_lat(iii),~] = m_fdist(ACO_lon,ACO_lat,angle,L);
    tx_lon(iii) = tx_lon(iii)-360;
    hold on
    %scatter(tx_lon(iii),tx_lat(iii),'po')
     plot([tx_lon(iii) ACO_lon],[tx_lat(iii) ACO_lat],'Linewidth',1,'Color','b')
    
end

subplot(1,3,3)
% 3rd circle
for iii = n1+n2+1:n2+n1+n3
    angle = 360/n3*(iii-1);
    [tx_lon(iii),tx_lat(iii),~] = m_fdist(ACO_lon,ACO_lat,angle,L);
    tx_lon(iii) = tx_lon(iii)-360;
    hold on
    %scatter(tx_lon(iii),tx_lat(iii),'po')
     plot([tx_lon(iii) ACO_lon],[tx_lat(iii) ACO_lat],'Linewidth',1,'Color','b')
    
end
%}
%% radial
n4 = 50;
% east
angle = 0;
L = 0;
for pp =1:n4
    L = L+26000/n4
    for iii = 1:4
        subplot(2,2,iii)
        [tx_lon(tx_ind),tx_lat(tx_ind),~] = m_fdist(ACO_lon,ACO_lat,angle,L);
        tx_lon(tx_ind) = tx_lon(tx_ind)-360;
        hold on 
        plot([tx_lon(tx_ind) ACO_lon],[tx_lat(tx_ind) ACO_lat],'Linewidth',1,'Color','r')

    end
    tx_ind = tx_ind+1; 
end
% west

angle = 270;
L = 0;
for pp =1:n4
    L = L+26000/n4;
    for iii = 1:4
        subplot(2,2,iii)
        [tx_lon(tx_ind),tx_lat(tx_ind),~] = m_fdist(ACO_lon,ACO_lat,angle,L);
        tx_lon(tx_ind) = tx_lon(tx_ind)-360;
        hold on 
        plot([tx_lon(tx_ind) ACO_lon],[tx_lat(tx_ind) ACO_lat],'Linewidth',1,'Color','r')

    end
    tx_ind = tx_ind+1;
end 


% north
angle = 0;
L = 0;
for pp =1:n4
    L = L+26000/n4;
    for iii = 1:4
        subplot(2,2,iii)
        [tx_lon(tx_ind),tx_lat(tx_ind),~] = m_fdist(ACO_lon,ACO_lat,angle,L);
        tx_lon(tx_ind) = tx_lon(tx_ind)-360;
        hold on 
        plot([tx_lon(tx_ind) ACO_lon],[tx_lat(tx_ind) ACO_lat],'Linewidth',1,'Color','r')

    end
    tx_ind = tx_ind+1;
end 

%west
angle = 180;
L = 0;
for pp =1:n4
    L = L+26000/n4;
    for iii = 1:4
        subplot(2,2,iii)
        [tx_lon(tx_ind),tx_lat(tx_ind),~] = m_fdist(ACO_lon,ACO_lat,angle,L);
        tx_lon(tx_ind) = tx_lon(tx_ind)-360;
        hold on 
        plot([tx_lon(tx_ind) ACO_lon],[tx_lat(tx_ind) ACO_lat],'Linewidth',1,'Color','r')

    end
    tx_ind = tx_ind+1;
end 

%% concentric circles
R1 = 26000;
R2 = 20000;
R3 = 15000;
R4 = 10000;
R5 = 5000;
R = [R1 R2 R3 R4 R5];
spacing = 500;
n1 = floor(R1/500*2*pi);
n2 = floor(R2/500*2*pi);
n3 = floor(R3/500*2*pi);
n4 = floor(R4/500*2*pi);
n5 = floor(R5/500*2*pi);
n = [n1 n2 n3 n4 n5];
n_tot = sum(n);

figure(1)
angle = 0;
tx_ind = 1;
tx_lon = [];
tx_lat = [];

for iii = 1:length(n)
    for jjj = 1:n(iii)
        for kkk = 1:4
            subplot(2,2,kkk)
            hold on
            angle = 360/n(iii)*(jjj-1);
            [tx_lon(tx_ind),tx_lat(tx_ind),~] = m_fdist(ACO_lon,ACO_lat,angle,R(iii));
            tx_lon(tx_ind) = tx_lon(tx_ind)-360;
            hold on
            %scatter(tx_lon(iii),tx_lat(iii),'po')
            plot([tx_lon(tx_ind) ACO_lon],[tx_lat(tx_ind) ACO_lat],'Linewidth',1,'Color','m')
            
        end
        tx_ind =tx_ind+1;
    end
end



%{
n6 = 50;
ang1 = 42;
ang2 = 48;
angle = ang1:(ang2-ang1)/(n6-1):ang2;
for jjj = 1:n6
    for kkk = 1:4
        subplot(2,2,kkk)
        hold on
        [tx_lon(tx_ind),tx_lat(tx_ind),~] = m_fdist(ACO_lon,ACO_lat,angle(jjj),26000);
        tx_lon(tx_ind) = tx_lon(tx_ind)-360;
        hold on
        %scatter(tx_lon(iii),tx_lat(iii),'po')
        plot([tx_lon(tx_ind) ACO_lon],[tx_lat(tx_ind) ACO_lat],'Linewidth',1,'Color','k')
        
    end
    tx_ind =tx_ind+1;
end

n6 = 50;
ang1 = 26;
ang2 = 35;
angle = ang1:(ang2-ang1)/(n6-1):ang2;
for jjj = 1:n6
    for kkk = 1:4
        subplot(2,2,kkk)
        hold on
        [tx_lon(tx_ind),tx_lat(tx_ind),~] = m_fdist(ACO_lon,ACO_lat,angle(jjj),27000);
        tx_lon(tx_ind) = tx_lon(tx_ind)-360;
        hold on
        %scatter(tx_lon(iii),tx_lat(iii),'po')
        plot([tx_lon(tx_ind) ACO_lon],[tx_lat(tx_ind) ACO_lat],'Linewidth',1,'Color','k')
        
    end
    tx_ind =tx_ind+1;
end


n6 = 50;
ang1 = 52;
ang2 = 60;
angle = ang1:(ang2-ang1)/(n6-1):ang2;
for jjj = 1:n6
    for kkk = 1:4
        subplot(2,2,kkk)
        hold on
        [tx_lon(tx_ind),tx_lat(tx_ind),~] = m_fdist(ACO_lon,ACO_lat,angle(jjj),27000);
        tx_lon(tx_ind) = tx_lon(tx_ind)-360;
        hold on
        %scatter(tx_lon(iii),tx_lat(iii),'po')
        plot([tx_lon(tx_ind) ACO_lon],[tx_lat(tx_ind) ACO_lat],'Linewidth',1,'Color','k')
        
    end
    tx_ind =tx_ind+1;
end
%}

%% range 
range = [];
for iii =1:length(tx_lon)
   range(iii) = dist([tx_lon(iii) ACO_lon],[tx_lat(iii) ACO_lat]); 
end
%}
%% For testing G integrand
%{
% Radial
n_test = 20;
angle = 90;
L = 0;
tx_lon = [];
tx_lat = [];
for pp =1:3
    subplot(1,3,pp)
    L = 1;
    for iii = 1:n_test
        [tx_lon(iii),tx_lat(iii),~] = m_fdist(ACO_lon,ACO_lat,angle,L);
        tx_lon(iii) = tx_lon(iii)-360;
        hold on 
        plot([tx_lon(iii) ACO_lon],[tx_lat(iii) ACO_lat],'Linewidth',1,'Color','r')
        L = L+25000/(n_test-1);
    end
end 

% Angled Radial
n_test2 = 10;
angle = 45;
L = 0;
for pp =1:3
    subplot(1,3,pp)
    L = 1;
    for iii = n_test+1:n_test+n_test2
        [tx_lon(iii),tx_lat(iii),~] = m_fdist(ACO_lon,ACO_lat,angle,L);
        tx_lon(iii) = tx_lon(iii)-360;
        hold on 
        plot([tx_lon(iii) ACO_lon],[tx_lat(iii) ACO_lat],'Linewidth',1,'Color','r')
        L = L+25000/(n_test2-1);
    end
end 

% circle
n_test3 = 11;
L = 25000;
for pp =1:3
    subplot(1,3,pp)
    angle = 30;
    for iii = n_test+n_test2+1:n_test+n_test2+n_test3
        [tx_lon(iii),tx_lat(iii),~] = m_fdist(ACO_lon,ACO_lat,angle,L);
        tx_lon(iii) = tx_lon(iii)-360;
        hold on 
        plot([tx_lon(iii) ACO_lon],[tx_lat(iii) ACO_lat],'Linewidth',1,'Color','r')
        angle = angle+30;
        
    end
end 
% circle
n_test4 = 11;
L = 15000;
for pp =1:3
    subplot(1,3,pp)
    angle = 30;
    for iii = n_test+n_test2+n_test3+1:n_test+n_test2+n_test3+n_test4
        [tx_lon(iii),tx_lat(iii),~] = m_fdist(ACO_lon,ACO_lat,angle,L);
        tx_lon(iii) = tx_lon(iii)-360;
        hold on 
        plot([tx_lon(iii) ACO_lon],[tx_lat(iii) ACO_lat],'Linewidth',1,'Color','r')
        angle = angle+30;
        
    end
end 

% circle
n_test5 = 11;
L = 5000;
for pp =1:3
    subplot(1,3,pp)
    angle = 30;
    for iii = n_test+n_test2+n_test3+n_test4+1:n_test+n_test2+n_test3+n_test4+n_test5
        [tx_lon(iii),tx_lat(iii),~] = m_fdist(ACO_lon,ACO_lat,angle,L);
        tx_lon(iii) = tx_lon(iii)-360;
        hold on 
        plot([tx_lon(iii) ACO_lon],[tx_lat(iii) ACO_lat],'Linewidth',1,'Color','r')
        angle = angle+30;
        
    end
end 


%}
for iii = 1:length(tx_lon)
            hold on
            plot([tx_lon(iii) ACO_lon],[tx_lat(iii) ACO_lat],'Linewidth',1,'Color','m')
            
end
%% Obervation Matrix
G1 = zeros(length(tx_lon),grid_num*grid_num);
G2 = zeros(length(tx_lon),grid_num*grid_num);
G3 = zeros(length(tx_lon),grid_num*grid_num);
G4 = zeros(length(tx_lon),grid_num*grid_num);
M = zeros(length(tx_lon),1);
tx_altitude = -6;

% figure(3)
% clf
% set(gcf,'name','Integrad Value','Units','normalized','Position',[0 0.5 0.4 0.6]);
% figure(44)
% clf
% set(gcf,'Units','normalized','Position',[0.5 0.5 0.4 0.6]);
% indx_srt = ceil(grid_num*0.1);
% indx_stp = grid_num - ceil(grid_num*0.1);
% loop over each ray
for iii=1:length(tx_lon)
    iii
%     xdist  = distance(tx_lat(iii),tx_lon(iii),ACO_lat, ACO_lon,referenceEllipsoid('wgs84'))
    % find pixels and the distance in each pixel
    [G_of_nray1,total_pixel_num,~,~,intd_G1,z_cross,intd_G_all1,z,~,~]=obs_matrix3D_v2(tx_lat(iii),tx_lon(iii),tx_altitude,ACO_lat,ACO_lon,ACO_depth,x_node,y_node,1,"simulation");
    [G_of_nray2,~,~,~,intd_G2,~,intd_G_all2,~,~,~]=obs_matrix3D_v2(tx_lat(iii),tx_lon(iii),tx_altitude,ACO_lat,ACO_lon,ACO_depth,x_node,y_node,2,"simulation");
    [G_of_nray3,~,~,~,intd_G3,~,intd_G_all3,~,~,~]=obs_matrix3D_v2(tx_lat(iii),tx_lon(iii),tx_altitude,ACO_lat,ACO_lon,ACO_depth,x_node,y_node,3,"simulation");
    [G_of_nray4,~,~,~,intd_G4,~,intd_G_all4,~,~,~]=obs_matrix3D_v2(tx_lat(iii),tx_lon(iii),tx_altitude,ACO_lat,ACO_lon,ACO_depth,x_node,y_node,4,"simulation");
    % form matrix G
    G1(iii,:) =real(G_of_nray1);
    G2(iii,:) = real(G_of_nray2);
    G3(iii,:) = real(G_of_nray3);
    G4(iii,:) = real(G_of_nray4);
    M(iii) = length(total_pixel_num)-1;
    
    %{
    intd_G_remat1 = repmat(intd_G1',4,1);
    intd_G_remat2 = repmat(intd_G2',4,1);
    intd_G_remat3 = repmat(intd_G3',4,1);
    intd_G_remat4 = repmat(intd_G4',4,1);
    
    
    figure(3)
    clf
    subplot(1,4,1)
    cla
    plot(intd_G_remat1(:),[ones(1,2)*z_cross(1) sort(repmat(z_cross(2:end-1),1,4)) ones(1,2)*z_cross(end)],'Color','b')
    hold on
    plot(intd_G_all1*120,linspace(z_cross(1),z_cross(end),length(intd_G_all1)),'--k')
    grid on
    axis tight
    set(gca,'YDir','reverse','fontsize',13)
    xlim([-1e-3 1e-3])
    title('Mode1')
    ylabel('Depth (m)')
    
    subplot(1,4,2)
    cla
    plot(intd_G_remat2(:),[ones(1,2)*z_cross(1) sort(repmat(z_cross(2:end-1),1,4)) ones(1,2)*z_cross(end)],'Color','b')
    hold on
    plot(intd_G_all2*120,linspace(z_cross(1),z_cross(end),length(intd_G_all2)),'--k')
    grid on
    axis tight
set(gca,'YDir','reverse','fontsize',13)
    xlim([-1e-3 1e-3])
    title('Mode2')
    
    subplot(1,4,3)
    cla
    plot(intd_G_remat3(:),[ones(1,2)*z_cross(1) sort(repmat(z_cross(2:end-1),1,4)) ones(1,2)*z_cross(end)],'Color','b')
    hold on
    plot(intd_G_all3*120,linspace(z_cross(1),z_cross(end),length(intd_G_all3)),'--k')
    grid on
    axis tight
set(gca,'YDir','reverse','fontsize',13)
    xlim([-1e-3 1e-3])
    title('Mode3')
    
    subplot(1,4,4)
    cla
    plot(intd_G_remat4(:),[ones(1,2)*z_cross(1) sort(repmat(z_cross(2:end-1),1,4)) ones(1,2)*z_cross(end)],'Color','b')
    hold on
    plot(intd_G_all4*120,linspace(z_cross(1),z_cross(end),length(intd_G_all4)),'--k')
    grid on
    axis tight
  set(gca,'YDir','reverse','fontsize',13)
    xlim([-1e-3 1e-3])
    title('Mode4')
    
    t = annotation('textbox',[.48 .9 .1 .1],'String',sprintf('%.0f',xdist)+" m",'Fontsize',12,'BackgroundColor','white');
    
    figure(44)
    rng(1)
    imagesc(x_cen(indx_srt:indx_stp),y_cen(indx_srt:indx_stp),0.5*rand(length(x_cen(indx_srt:indx_stp)),length(x_cen(indx_srt:indx_stp))))
    colormap gray
    hold on
    scatter(ACO_lon,ACO_lat,200,'pr','filled')
    scatter(tx_lon(iii),tx_lat(iii),'or','filled')
    set(gca,'YDir','normal')
    plot(x_cir,y_cir,'r')
    caxis([0 0.3])
    pause
    %}
end
% close(gcf)
% close(gcf)

G = [G1 G2 G3 G4];

%% plot elements of the observation matrix

figure(4)
subplot(4,1,1)
plot(1:length(m01),sum((G1),1))
% xticks(1:1:length(m01))
grid on
ylabel('meter')
xlabel('pixel')
xlim([1 length(m01)])
title('Obserbvation Matrix of the first mode')

subplot(4,1,2)
plot(1:length(m02),sum((G2),1))
% xticks(1:10:length(m02))
grid on
ylabel('meter')
xlabel('pixel')
xlim([1 length(m02)])
title('Obserbvation Matrix of the second mode')

subplot(4,1,3)
plot(1:length(m02),sum((G3),1))
% xticks(1:10:length(m03))
grid on
ylabel('meter')
xlabel('pixel')
xlim([1 length(m03)])
title('Obserbvation Matrix of the third mode')
subplot(4,1,4)
plot(1:length(m02),sum((G4),1))
% xticks(1:10:length(m03))
grid on
ylabel('meter')
xlabel('pixel')
xlim([1 length(m03)])
title('Obserbvation Matrix of the fourth mode')

figure(9191)
plot(1:length(m02),sum((G1),1)+sum((G2),1)+sum((G3),1)+sum((G4),1))
grid on
ylabel('meter')
xlabel('pixel')
xlim([1 length(m03)])
title('Obserbvation Matrix')

sum_G_pixelwise = real(sum((G1),1)+sum((G2),1)+sum((G3),1)+sum((G4),1))';

reshape_sum_G = reshape(sum_G_pixelwise,grid_num,grid_num)';
reshape_sum_Gu = triu(reshape_sum_G);
reshape_sum_Gl = tril(reshape_sum_G);
for iii = 1:grid_num
    for jjj = 1:grid_num
        reshape_sum_Gl_flip(iii,jjj) = reshape_sum_Gl(jjj,iii);
        
    end
end
figure(19)
imagesc(x_cen,y_cen,reshape_sum_G)
colormap jet
cbar = colorbar;
cbar.Label.String = 'G';
title('Pixel-wise Sum of integrand G ')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
set(gca,'YDir','normal','fontsize',12)
plot(x_cir,y_cir,'k')
% caxis([-0.03 0])
%}

%% Sensitivity vs pixel crossed

[M_sorted,ind_M] = sort(M);
M_ind = unique(M_sorted);
sumG1= sum(G1,2);
sumG2= sum(G2,2);
sumG3= sum(G3,2);
sumG4= sum(G4,2);
sumG1_sorted = sumG1(ind_M);
sumG2_sorted = sumG2(ind_M);
sumG3_sorted = sumG3(ind_M);
sumG4_sorted = sumG4(ind_M);
sens_m1 = zeros(length(M_ind),1);
sens_m2 = zeros(length(M_ind),1);
sens_m3 = zeros(length(M_ind),1);
sens_m4 = zeros(length(M_ind),1);
% Sensitivity
counter =1;
n = 1;
kp_ind = 1;
while true
    if M_sorted(counter+1) ~= M_sorted(counter)
        kp_ind(end+1) = counter;
        sens_m1(n) = mean(sumG1_sorted(kp_ind(n):kp_ind(n+1)));
        sens_m2(n) = mean(sumG2_sorted(kp_ind(n):kp_ind(n+1)));
        sens_m3(n) = mean(sumG3_sorted(kp_ind(n):kp_ind(n+1)));
        sens_m4(n) = mean(sumG4_sorted(kp_ind(n):kp_ind(n+1)));
        n = n+1;
    end
    counter = counter +1;
    if counter == length(M)
        kp_ind(end+1) = counter;
        sens_m1(n) = mean(sumG1_sorted(kp_ind(n):kp_ind(n+1)));
        sens_m2(n) = mean(sumG2_sorted(kp_ind(n):kp_ind(n+1)));
        sens_m3(n) = mean(sumG3_sorted(kp_ind(n):kp_ind(n+1)));
        sens_m4(n) = mean(sumG4_sorted(kp_ind(n):kp_ind(n+1)));
        break
    end
    
end

figure(10)
clf
plot(M_ind,abs(sens_m1),'k')
grid on
xlabel('Crossed Pixel')
hold on
plot(M_ind,abs(sens_m2),'r')
plot(M_ind,abs(sens_m3),'b')
plot(M_ind,abs(sens_m4),'g')
% scatter(M,abs(sumG1),'xk')
% scatter(M,abs(sumG2),'xr')
% scatter(M,abs(sumG3),'xb')
% scatter(M,abs(sumG4),'xg')
ylabel('|Sensitivity|')
title('Sensitivity of the TTP to the Number of Crossed Pixels')
legend('Mode1','Mode2','Mode3','Mode4','location','northwest')
set(gca,'fontsize',14)
%}
%% Row of G matrix 
%{
indx_srt = ceil(grid_num*0.1);
indx_stp = grid_num - ceil(grid_num*0.1);
figure(11)
clf
figure(121)
clf
for iii = 1:size(G,1)
    figure(11)
    plot(G(iii,:))
    ylim([-6 4]*1e-4)
    figure(121)
    rng(1)
    imagesc(x_cen(indx_srt:indx_stp),y_cen(indx_srt:indx_stp),0.5*rand(length(x_cen(indx_srt:indx_stp)),length(x_cen(indx_srt:indx_stp))))
    colormap gray
    hold on
    scatter(ACO_lon,ACO_lat,200,'pr','filled')
    scatter(tx_lon(iii),tx_lat(iii),'or','filled')
    set(gca,'YDir','normal')
    plot(x_cir,y_cir,'r')
    caxis([0 0.3])
    pause
end
%}
%% FInd the solution
% 6. Create measurement matrix from the observation matrix (TTP matrix)
d = G*m0;   % ms

% 7 Inversion process to recover the ss aperturbation field
% travel time perturbation error
Cd =diag(0.0002^2*ones(1,length(tx_lon)))     ; % ms
% 1. generalized inverse matrig G^-g
inv_U = (G*P*G'+Cd)^-1;
G_geninv = P*G'*(G*P*G'+Cd)^-1;
% recovered m
m_recov = G_geninv*d;
m_recov1 = m_recov(1:tot_grid_num);
m_recov2 = m_recov(tot_grid_num+1:2*tot_grid_num);
m_recov3 = m_recov(2*tot_grid_num+1:3*tot_grid_num);
m_recov4 = m_recov(3*tot_grid_num+1:4*tot_grid_num);

alpha0 = m_recov;
alpha_1 = alpha0(1:tot_grid_num);
alpha_2 = alpha0(tot_grid_num+1:2*tot_grid_num);
alpha_3 = alpha0(2*tot_grid_num+1:3*tot_grid_num);
alpha_4 = alpha0(3*tot_grid_num+1:4*tot_grid_num);

% Posteriori covariance
P_post = round((eye(tot_grid_num*4)-G_geninv*G)*P,10,'significant');
P_post_triu = triu(P_post);
P_post_tril = tril(P_post);
for iii = 1:tot_grid_num*4
    for jjj = 1:tot_grid_num*4
        P_post_trilflip(iii,jjj) = P_post_tril(jjj,iii);
        
    end
end
P_post_newu = (P_post_triu+P_post_triu')-diag(diag(P_post_triu));
P_post_newl = (P_post_tril+P_post_tril')-diag(diag(P_post_tril));
error_P = P_post_newu-P_post_newl;
error_Pu = P_post-P_post_newu;
P_post_new = (P_post_newu+P_post_newl)/2;
P_post_alpha = P_post_new;
P_post1 = P_post_new(1:tot_grid_num,1:tot_grid_num);
P_post2 = P_post_new(tot_grid_num+1:2*tot_grid_num,tot_grid_num+1:2*tot_grid_num);
P_post3 = P_post_new(2*tot_grid_num+1:3*tot_grid_num,2*tot_grid_num+1:3*tot_grid_num);
P_post4 = P_post_new(3*tot_grid_num+1:4*tot_grid_num,3*tot_grid_num+1:4*tot_grid_num);


%  Resolution Matrix
inv_PP0 = P_post*P^-1;
inv_PP_u = triu(inv_PP0);
inv_PP_l = tril(inv_PP0);
inv_PP_u = (inv_PP_u+inv_PP_u')-diag(diag(inv_PP_u));
inv_PP_l = (inv_PP_l+inv_PP_l')-diag(diag(inv_PP_l));
inv_PP = (inv_PP_u+inv_PP_l)/2;

Res_mat = eye(tot_grid_num*4)-inv_PP ;
% Res_mat = G_geninv*G;
Res_mat1 = Res_mat(1:tot_grid_num,1:tot_grid_num);
Res_mat2 = Res_mat(tot_grid_num+1:2*tot_grid_num,tot_grid_num+1:2*tot_grid_num);
Res_mat3 = Res_mat(2*tot_grid_num+1:3*tot_grid_num,2*tot_grid_num+1:3*tot_grid_num);
Res_mat4 = Res_mat(3*tot_grid_num+1:4*tot_grid_num,3*tot_grid_num+1:4*tot_grid_num);

% reshape the SS field
recov_SS1=reshape(m_recov1(1:end),grid_num,grid_num)'  ;
recov_SS2=reshape(m_recov2(1:end),grid_num,grid_num)'  ;
recov_SS3=reshape(m_recov3(1:end),grid_num,grid_num)'  ;
recov_SS4=reshape(m_recov4(1:end),grid_num,grid_num)'  ;

% residual
d_recov = G*m_recov;
residual = d-d_recov;

% model resolution
model_res1 = reshape(diag(Res_mat1),grid_num,grid_num)';
model_res2 = reshape(diag(Res_mat2),grid_num,grid_num)';
model_res3 = reshape(diag(Res_mat3),grid_num,grid_num)';
model_res4 = reshape(diag(Res_mat4),grid_num,grid_num)';


% variance reduction
prior_SD1 = reshape(sqrt(diag(P_mode1)),grid_num,grid_num);
prior_SD2 = reshape(sqrt(diag(P_mode2)),grid_num,grid_num);
prior_SD3 = reshape(sqrt(diag(P_mode3)),grid_num,grid_num);
prior_SD4 = reshape(sqrt(diag(P_mode4)),grid_num,grid_num);
prior_SD = sqrt(prior_SD1.^2+prior_SD2.^2+prior_SD3.^2+prior_SD4.^2);

alpha_SD1 = reshape(sqrt(diag(P_post1)),grid_num,grid_num);
alpha_SD2 = reshape(sqrt(diag(P_post2)),grid_num,grid_num);
alpha_SD3 = reshape(sqrt(diag(P_post3)),grid_num,grid_num);
alpha_SD4 = reshape(sqrt(diag(P_post4)),grid_num,grid_num);
alpha_SD = sqrt(alpha_SD1.^2+alpha_SD2.^2+alpha_SD3.^2+alpha_SD4.^2);
sd_reduction_alpha = alpha_SD./prior_SD*100;

post_SD1 = reshape(sqrt(diag(P_post1)),grid_num,grid_num);
post_SD2 = reshape(sqrt(diag(P_post2)),grid_num,grid_num);
post_SD3 = reshape(sqrt(diag(P_post3)),grid_num,grid_num);
post_SD4 = reshape(sqrt(diag(P_post4)),grid_num,grid_num);

sd_reduction1 = (post_SD1)./prior_SD1*100;
sd_reduction2 = (post_SD2)./prior_SD2*100;
sd_reduction3 = (post_SD3)./prior_SD3*100;
sd_reduction4 = (post_SD4)./prior_SD4*100;

%Combine 4 modes
load('EOF_SS.mat')
f1 = EOF_SS.mode1;
f2 = EOF_SS.mode2;
f3 = EOF_SS.mode3;
f4 = EOF_SS.mode4;
f1 = [f1;f1(end)*ones(2,1)];
f2 = [f2;f2(end)*ones(2,1)];
f3 = [f3;f3(end)*ones(2,1)];
f4 = [f4;f4(end)*ones(2,1)];
F = [f1 f2 f3 f4];
alpha = [alpha_1';alpha_2';alpha_3';alpha_4'];
SSP = F*alpha;

% depth averaged modes
f1_avg = sum(f1(1:end-1).*(z(2:end)-z(1:end-1)))/z(end);
f2_avg = sum(f2(1:end-1).*(z(2:end)-z(1:end-1)))/z(end);
f3_avg = sum(f3(1:end-1).*(z(2:end)-z(1:end-1)))/z(end);
f4_avg = sum(f4(1:end-1).*(z(2:end)-z(1:end-1)))/z(end);

% depth averaging and mode combining operator
F_operator = [diag(ones(1,625)*f1_avg) diag(ones(1,625)*f2_avg) diag(ones(1,625)*f3_avg) diag(ones(1,625)*f4_avg)];
SSP_d_avg2 = F_operator*alpha0;
SSP_d_avg1 = F_operator*m0;
P_prior_d_avg = F_operator*P*F_operator';
P_SSP_d_avg = F_operator*P_post_alpha*F_operator';


% recalculate the SSP field
% reshape the SS field
recov_SS_d_avg=reshape(SSP_d_avg2(1:end),grid_num,grid_num)'  ;
initial_SS_d_avg=reshape(SSP_d_avg1(1:end),grid_num,grid_num)'  ;

% Resolution Matrix

Res_mat_d_avg = eye(tot_grid_num)-P_SSP_d_avg*P_prior_d_avg^-1;


% 5. Posteriori covariance
prior_SD_d_avg = reshape(sqrt(diag(P_prior_d_avg)),grid_num,grid_num);

% residual
% d_recov = G*alpha0;
% residual = d-d_recov;

% model resolution
model_res_d_avg = reshape(diag(Res_mat_d_avg),grid_num,grid_num)';


% variance reduction
post_SD_d_avg = reshape(sqrt(diag(P_SSP_d_avg)),grid_num,grid_num)';
sd_reduction_d_avg = abs((post_SD_d_avg)./prior_SD_d_avg*100);%./prior_SD1*100;


%% 8 plot the recovered ss pertrubation field 4 modes
figure(31)
clf
set(gcf,'name','2-D sound speed perturbation field','Units','normalized','Position',[0 0 0.55 .6])
subplot(2,2,1)
imagesc(x_cen,y_cen,recov_SS1)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title('Recovered SSP Field Mode1')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',14)
xlabel('Long')
ylabel('Lat')
caxis([-4 4])

subplot(2,2,2)
imagesc(x_cen,y_cen,recov_SS2)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title('Recovered SSP Field: Mode2')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',14)
xlabel('Long')
ylabel('Lat')
caxis([-4 4])

subplot(2,2,3)
imagesc(x_cen,y_cen,recov_SS3)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title('Recovered SSP Field: Mode3')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',14)
xlabel('Long')
ylabel('Lat')
caxis([-4 4])

subplot(2,2,4)
imagesc(x_cen,y_cen,recov_SS4)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title('Recovered SSP Field: Mode4')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',14)
xlabel('Long')
ylabel('Lat')
caxis([-4 4])

%% 9 Depth-Averaged solution
figure(2)
clf
set(gcf,'name','2-D sound speed perturbation field','Units','normalized','Position',[0.5 0.1 0.35 .4])
imagesc(x_cen,y_cen,initial_SS_d_avg)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title('Ideal Depth-Averaged SSP Field')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',12)
xlabel('Long')
ylabel('Lat')
caxis([-.8 .8])

figure(3)
clf
set(gcf,'name','2-D sound speed perturbation field','Units','normalized','Position',[0 0.1 0.35 .4])
imagesc(x_cen,y_cen,recov_SS_d_avg)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title('Recovered Depth-Averaged SSP Field')
title('Ideal Sampling SSP Field')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',12)
xlabel('Long')
ylabel('Lat')
caxis([-0.8 .8])


%% Difference

figure(4)
clf
set(gcf,'name','2-D sound speed perturbation field Difference','Units','normalized','Position',[0 1 0.55 .6])
subplot(2,2,1)

imagesc(x_cen,y_cen,sim_SS1 - recov_SS1)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title('Difference Mode 1')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',12)
xlabel('Long')
ylabel('Lat')
caxis([-4 4])
 
subplot(2,2,2)
imagesc(x_cen,y_cen,sim_SS2 - recov_SS2)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title('Difference Mode 2')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',12)
xlabel('Long')
ylabel('Lat')
  caxis([-4 4])
  
subplot(2,2,3)
imagesc(x_cen,y_cen,sim_SS3 - recov_SS3)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title('Difference Mode 3')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',12)
xlabel('Long')
ylabel('Lat')
caxis([-4 4])

  
subplot(2,2,4)
imagesc(x_cen,y_cen,sim_SS4 - recov_SS4)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title('Difference Mode 3')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',12)
xlabel('Long')
ylabel('Lat')
caxis([-4 4])
%%
figure(42)
clf
set(gcf,'name','2-D sound speed perturbation field','Units','normalized','Position',[0 0.1 0.35 .4])
imagesc(x_cen,y_cen,abs((recov_SS_d_avg-initial_SS_d_avg))./post_SD_d_avg*100)
% imagesc(x_cen,y_cen,(recov_SS_d_avg-initial_SS_d_avg))
colormap jet
cbar = colorbar;
cbar.Label.String = '% (Difference/Posterior Error) ';
% cbar.Label.String = 'Sound Speed Perturbation (m/s) ';
title('Ideal Sampling Difference SSP Field')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',14)
xlabel('Long')
ylabel('Lat')
caxis([0 100])
% caxis([-0.08 0.08])
%}
%% Solution Quality
% RMS error reduction
%%% Plot
figure(5)
clf
set(gcf,'Units','normalized','Position',[.4 1 0.5 .55])
subplot(2,2,1)
imagesc(x_cen,y_cen,sd_reduction1)
colormap jet
cbar = colorbar;
cbar.Label.String = '%';
title('RMS Error Reduction 1st mode')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal')
set(gca,'YDir','normal','fontsize',13)
xlabel('Long')
ylabel('Lat')
caxis([0 100])

subplot(2,2,2)
imagesc(x_cen,y_cen,sd_reduction2)
colormap jet
cbar = colorbar;
cbar.Label.String = '%';
title('RMS Error Reduction 2nd mode')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',13)
xlabel('Long')
ylabel('Lat')
caxis([0 100])

subplot(2,2,3)
imagesc(x_cen,y_cen,sd_reduction3)
colormap jet
cbar = colorbar;
cbar.Label.String = '%';
title('RMS Error Reduction  3rd mode')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',13)
xlabel('Long')
ylabel('Lat')
caxis([0 100])

subplot(2,2,4)
imagesc(x_cen,y_cen,sd_reduction4)
colormap jet
cbar = colorbar;
cbar.Label.String = '%';
title('RMS Error Reduction  4th mode')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',13)
xlabel('Long')
ylabel('Lat')
caxis([0 100])


figure(33)
imagesc(x_cen,y_cen,sd_reduction_d_avg)
colormap jet
cbar = colorbar;
cbar.Label.String = '%';
title('RMS Error: Depth-Averaged Solution (icListen)')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',13)
xlabel('Long')
ylabel('Lat')
caxis([0 100])

figure(32)
imagesc(x_cen,y_cen,sd_reduction_alpha)
colormap jet
cbar = colorbar;
cbar.Label.String = '%';
title('RMS Error: Separate-mode Solution')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',13)
xlabel('Long')
ylabel('Lat')
caxis([0 100])
%% Resolution Matrix
figure(6)
clf
set(gcf,'Units','normalized','Position',[.4 .3 0.55 .6])
subplot(2,2,1)
imagesc(x_cen,y_cen,model_res1)
hold on
plot(x_cir,y_cir,'k')
axis tight
grid off
title('Model Resolution 1st mode')
colorbar
set(gca,'YDir','normal','fontsize',14)
xlabel('Long')
ylabel('Lat')
caxis([0 1])
% caxis([min(P_post(:)) max(P_post(:))])
colormap jet

subplot(2,2,2)
imagesc(x_cen,y_cen,model_res2)
hold on
plot(x_cir,y_cir,'k')
axis tight
grid off
title('Model Resolution 2nd mode')
colorbar
caxis([0 1])
% caxis([min(P_post(:)) max(P_post(:))])
colormap jet
set(gca,'YDir','normal','fontsize',14)
xlabel('Long')
ylabel('Lat')

subplot(2,2,3)
imagesc(x_cen,y_cen,model_res3)
hold on
plot(x_cir,y_cir,'k')
axis tight
grid off
title('Model Resolution 3rd mode')
colorbar
xlabel('Long')
ylabel('Lat')
caxis([0 1])
% caxis([min(P_post(:)) max(P_post(:))])
colormap jet
set(gca,'YDir','normal','fontsize',14)

subplot(2,2,4)
imagesc(x_cen,y_cen,model_res4)
hold on
plot(x_cir,y_cir,'k')
axis tight
grid off
title('Model Resolution 4th mode')
colorbar
caxis([0 1])
% caxis([min(P_post(:)) max(P_post(:))])
colormap jet
set(gca,'YDir','normal','fontsize',14)
xlabel('Long')
ylabel('Lat')

figure(62)
clf
set(gcf,'Units','normalized','Position',[.4 .3 0.4 .5])
imagesc(x_cen,y_cen,model_res_d_avg)
hold on
plot(x_cir,y_cir,'k')
axis tight
grid off
title('Model Resolution: Combined Modes')
colorbar
caxis([0 1])
xlabel('Long')
ylabel('Lat')
% caxis([min(P_post(:)) max(P_post(:))])
colormap jet
set(gca,'YDir','normal','fontsize',12)
%% Resolution vs Range
kk = 1;
sur_dist = [];
azmth = [];
for ii = 1:length(y_cen)
    for jj = 1:length(x_cen)
            sur_dist(kk) = dist([x_cen(jj) ACO_lon],[y_cen(ii) ACO_lat]);
            azmth(kk) = azimuth(ACO_lat,ACO_lon,y_cen(ii),x_cen(jj));
            kk = kk+1;
    end
end      
[sur_dist,I] = sort(sur_dist);
res_m1 = model_res1(:);
res_m2 = model_res2(:);
res_m3 = model_res3(:);
res_m4 = model_res4(:);
res_d_avg = model_res_d_avg(:);
res_m1 = res_m1(I);
res_m2 = res_m2(I);
res_m3 = res_m3(I);
res_m4 = res_m4(I);
res_d_avg = res_d_avg(I);
figure(89)
clf
plot(sur_dist/1000,res_m1)
hold on
plot(sur_dist/1000,res_m2,'--g')
plot(sur_dist/1000,res_m3,'-r')
plot(sur_dist/1000,res_m3,'Color',[1 89/255 .1])
plot(sur_dist/1000,res_m4,'--','Color',[.5 0 .5])
line([0 50],[0 0],'color','k')
grid minor
xlim([0 30])
xlabel('Range (km)')
ylabel('Resolution')
legend('Mode1','Mode2','Mode3','Mode4')
title('Resolution versus Range (Covariance Length = 20km)')
set(gca,'fontsize',14)


figure(819)
clf
plot(sur_dist/1000,res_d_avg)
grid minor
hold on
scatter(sur_dist/1000,res_d_avg,'b*')
xlim([0 30])
xlabel('Range (km)')
ylabel('Resolution')
title({'Resolution versus Range: Depth-averaged solution',' (Correlation Length = 8.6km)'})
set(gca,'fontsize',14)
line([0 50],[0 0],'color','k')
ylim([0 1])
% subplot(3,2,5)
% imagesc(x_cen,y_cen,recov_SS1+recov_SS2)
% colormap jet
% cbar = colorbar;
% cbar.Label.String = 'Sound Speed Perturbation (m/s)';
% title('Recovered SS Perturbation Field (2modes)')
% hold on
% scatter(ACO_lon,ACO_lat,200,'pk','filled')
% plot(x_cir,y_cir,'k')
% set(gca,'YDir','normal')


% cd /Users/testuser/Documents/Oct2018Cruise/Figure/Inversion/2D
% s1 = sprintf('px_%i_ray_%i_P',grid_num,ntot);
% s2 = sprintf('px_%i_ray_%i',grid_num,ntot);
% saveas(f1,[s1 '.png'])
% saveas(f2,[s2 '.png'])

%%
figure(90)
for iii =400:475
    plot(SSP(:,iii),z.*ones(2407,1))
    grid on
    set(gca,'YDir','reverse')
    xlim([-1 6])
    pause
end
%% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P,var_c] = gaussian_cov_mx(x,y,mode)
% load singular vaues of each mode
cd /Users/testuser/Documents/RAP_Script/Data
load EOF_SS.mat
var_c = (eval(['EOF_SS.singv' num2str(mode)]))^2;
% calculate distance between pixel
len = length(x);
range = zeros(len,len);
P = zeros(len^2,len^2);

%X and Y matrix
x_matrix=repmat(x,1,len);
y_matrix=reshape(repmat(y,len,1),1,len*len);

% Gaussian cov length
sigma = 20000;

for ii=1:size(x_matrix,2) 
    for jj=ii:size(x_matrix,2) 
        
            
        %Distance between pixel centers
        pixel_distance=dist([y_matrix(ii) y_matrix(jj)],[x_matrix(ii) x_matrix(jj)]);
        %Gaussian covariance matrix
        P(ii,jj)=(var_c)*exp(-(pixel_distance^2)/(sigma*sigma));
        
    end
end
P = P+P'-diag(diag(P));
end



function [P,var_c] = gaussian_cov_sim(x,y,mode)
% load singular vaues of each mode
cd /Users/testuser/Documents/RAP_Script/Data
load EOF_SS.mat
var_c = (eval(['EOF_SS.singv' num2str(mode)]))^2;
% calculate distance between pixel
len = length(x);
range = zeros(len,len);
P = zeros(len^2,len^2);

%X and Y matrix
x_matrix=repmat(x,1,len);
y_matrix=reshape(repmat(y,len,1),1,len*len);

% Gaussian cov length
sigma = 25000;

for ii=1:size(x_matrix,2) 
    for jj=ii:size(x_matrix,2) 
        
            
        %Distance between pixel centers
        pixel_distance=dist([y_matrix(ii) y_matrix(jj)],[x_matrix(ii) x_matrix(jj)]);
        %Gaussian covariance matrix
        P(ii,jj)=(var_c)*exp(-(pixel_distance^2)/(sigma*sigma));
        
    end
end
P = P+P'-diag(diag(P));
end
