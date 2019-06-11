% 3-D  simulation
% simulate ocean sound speed perturbation field (m) with Gaussian model uncertainty
% Impose empirical orthogonal functions (EOF) on the vertical structure of ocean sound speed variation
% Exclude hydrophone position offsets

close all
clearvars
%%
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
%  P = P_mode1;
 P = [P_mode1 zeros(dimP) zeros(dimP) zeros(dimP); zeros(dimP) P_mode2 zeros(dimP) zeros(dimP);zeros(dimP) zeros(dimP) P_mode3 zeros(dimP); zeros(dimP) zeros(dimP) zeros(dimP) P_mode4];
 
%% Priori cov matrix
figure(1)
set(gcf,'name','Covariance Matrix','Units','normalized','Position',[0.3 .5 0.35 0.45])
imagesc(P)
axis tight
grid off
title('Gaussian Covariance Matrix')
colorbar
caxis([0 1])
colormap jet

%% 3 Simulate SS perturbation field
%Cholesky Decomposition to create a simulated SS field
% P is not eactly symmetric, need to asjust by nearestSPD
P1=P;
while exist('D')==0
try
D=chol(P1,'lower');
end
P1=nearestSPD(P1);
end
rng(12)
rand_sst=normrnd(0,1,length(P),1);      %Random sound speed perturbation vector with normal distribution
m0 = D*rand_sst;                                %Simulated SS field for each pixel
m01 = m0(1:tot_grid_num);
m02 = m0(tot_grid_num+1:2*tot_grid_num);
m03 = m0(2*tot_grid_num+1:3*tot_grid_num);
m04 = m0(3*tot_grid_num+1:4*tot_grid_num);
%% delta function at the center
mid_p = (tot_grid_num-1)/2;
m01 = [zeros(mid_p,1);3;zeros(mid_p,1)];
m02 = [zeros(mid_p,1);2;zeros(mid_p,1)];
m03 = [zeros(mid_p,1);1;zeros(mid_p,1)];
m04 = [zeros(mid_p,1);0.5;zeros(mid_p,1)];
m0 = [m01;m02;m03;m04];

%% on the circle
pixel_num = 317;
m01 = [zeros(pixel_num,1);3;zeros(tot_grid_num-pixel_num-1,1)];
m02 = [zeros(pixel_num,1);2;zeros(tot_grid_num-pixel_num-1,1)];
m03 = [zeros(pixel_num,1);1;zeros(tot_grid_num-pixel_num-1,1)];
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

%% 4 Plot the Sound Speed Field
% draw circle
R = 25000;
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
set(gcf,'name','2-D sound speed perturbation field simulation','Units','normalized','Position',[0 1 0.5 .6])
subplot(2,2,1)
imagesc(x_cen,y_cen,sim_SS1)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title(sprintf('1st mode SS perturbation field'))
% title(sprintf('1st mode SS perturbation field \n(ss variance = %.3f (m/s)^{2})',var_c1))
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
set(gca,'YDir','normal')
plot(x_cir,y_cir,'k')
caxis_lim_priori = cbar.Limits;
caxis([-6 6])

subplot(2,2,2)
imagesc(x_cen,y_cen,sim_SS2)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title(sprintf('2nd mode SS perturbation field'));
% title(sprintf('2nd mode SS perturbation field \n(ss variance = %.3f (m/s)^{2})',var_c2))
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
set(gca,'YDir','normal')
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
set(gca,'YDir','normal')
plot(x_cir,y_cir,'k')
caxis_lim_priori3 = cbar.Limits;
caxis([-4 4])


subplot(2,2,4)
imagesc(x_cen,y_cen,sim_SS4)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Spexed Perturbation (m/s)';
title(sprintf('4th mode SS perturbation field'))
% title(sprintf('3rd mode SS perturbation field \n(ss variance = %.3f (m/s)^{2})',var_c3))
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
set(gca,'YDir','normal')
plot(x_cir,y_cir,'k')
caxis_lim_priori3 = cbar.Limits;
caxis([-4 4])

%% 5 Create Ray Paths
% from centers of pixels
% generate tx points
tx_lon = x_cen;
tx_lat = y_cen;

%range 
range = [];
cc = 1;
for iii =1:length(tx_lat)
    for jjj =1:length(tx_lon)
         range(iii,jjj) = dist([tx_lon(jjj) ACO_lon],[tx_lat(iii) ACO_lat]); 
    end
end

[ind_y,ind_x] = find((range<26500));
tx_lon = tx_lon(ind_x);
tx_lat= tx_lat(ind_y);

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

figure(2)
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
%% Obervation Matrix
G1 = zeros(length(tx_lon),grid_num*grid_num);
G2 = zeros(length(tx_lon),grid_num*grid_num);
G3 = zeros(length(tx_lon),grid_num*grid_num);
G4 = zeros(length(tx_lon),grid_num*grid_num);
M = zeros(length(tx_lon),1);
tx_altitude = -6;
% figure(3)
% clf
% set(gcf,'name','Integrad Value','Units','normalized','Position',[0 0.5 0.6 0.6]);
% figure(44)
% clf
% set(gcf,'Units','normalized','Position',[0.65 0.5 0.4 0.6]);
% loop over each ray
for iii=1:length(tx_lon)
    iii
%     xdist  = dist([tx_lon(iii) ACO_lon],[tx_lat(iii) ACO_lat])
    % find pixels and the distance in each pixel
    [G_of_nray1,total_pixel_num,~,~,~,z_cross,~,z,SS,~]=obs_matrix3D(tx_lat(iii),tx_lon(iii),tx_altitude,ACO_lat,ACO_lon,ACO_depth,x_node,y_node,1);
    [G_of_nray2,~,~,~,~,~,~]=obs_matrix3D(tx_lat(iii),tx_lon(iii),tx_altitude,ACO_lat,ACO_lon,ACO_depth,x_node,y_node,2);
    [G_of_nray3,~,~,~,~,~,~]=obs_matrix3D(tx_lat(iii),tx_lon(iii),tx_altitude,ACO_lat,ACO_lon,ACO_depth,x_node,y_node,3);
    [G_of_nray4,~,~,~,~,~,~]=obs_matrix3D(tx_lat(iii),tx_lon(iii),tx_altitude,ACO_lat,ACO_lon,ACO_depth,x_node,y_node,4);
    % form matrix G
    G1(iii,:) = G_of_nray1;
    G2(iii,:) = G_of_nray2;
    G3(iii,:) = G_of_nray3;
    G4(iii,:) = G_of_nray4;
    M(iii) = length(total_pixel_num)-1;
    
%     intd_G_remat1 = real(repmat(intd_G1',4,1));
%     intd_G_remat2 = real(repmat(intd_G2',4,1));
%     intd_G_remat3 = real(repmat(intd_G3',4,1));
%     intd_G_remat4 = real(repmat(intd_G4',4,1));
    
%     plot_integrand_G(intd_G_all1,intd_G_all2,intd_G_all3,intd_G_all4,intd_G_remat1,intd_G_remat2,intd_G_remat3,intd_G_remat4,z_cross,x_cen,y_cen,tx_lat(iii),tx_lon(iii),r_layer,xdist,grid_num,ACO_lon,ACO_lat,x_cir,y_cir)
end
% close(gcf)
% close(gcf)
G = real([G1 G2 G3 G4]);
% G = G1;
%%  plot elements of the observation matrix
%{
figure(4)
subplot(3,1,1)
plot(1:length(m01),sum((G1),1))
% xticks(1:1:length(m01))
grid on
ylabel('meter')
xlabel('pixel')
xlim([1 length(m01)])
title('Obserbvation Matrix of the first mode')

subplot(3,1,2)
plot(1:length(m02),sum((G2),1))
% xticks(1:10:length(m02))
grid on
ylabel('meter')
xlabel('pixel')
xlim([1 length(m02)])
title('Obserbvation Matrix of the second mode')

subplot(3,1,3)
plot(1:length(m02),sum((G3),1))
% xticks(1:10:length(m03))
grid on
ylabel('meter')
xlabel('pixel')
xlim([1 length(m03)])
title('Obserbvation Matrix of the third mode')
%}

%% Row of G matrix 

indx_srt = 1;% ceil(grid_num*0.1);
indx_stp = 25;%grid_num - ceil(grid_num*0.1);
figure(11)
clf
figure(121)
clf
for iii = 1:size(G,1)
    figure(11)
    subplot(4,1,1)
    plot(G(iii,1:625))
    ylim([-6 4]*1e-4)
    
    grid on
    axis tight
    xlim([130 450])
    
    subplot(4,1,2)
    plot(G(iii,626:1250))
    ylim([-6 4]*1e-4)
    
    grid on
    axis tight
    xlim([130 450])
    
    subplot(4,1,3)
    plot(G(iii,1251:1900))
    ylim([-6 4]*1e-4)
    
    grid on
    axis tight
    xlim([130 450])
    
    subplot(4,1,4)
    plot(G(iii,1901:2500))
    ylim([-6 4]*1e-4)
    xlim([180 450])
    grid on
    axis tight
    xlim([130 450])
    
    figure(121)
    rng(1)
    imagesc(x_cen(indx_srt:indx_stp),y_cen(indx_srt:indx_stp),0.5*rand(length(x_cen(indx_srt:indx_stp)),length(x_cen(indx_srt:indx_stp))))
    colormap gray
    hold on
    scatter(ACO_lon,ACO_lat,200,'pr','filled')
    scatter(tx_lon(iii),tx_lat(iii),'or','filled')
    line([ACO_lon tx_lon(iii)],[ACO_lat tx_lat(iii)],'color','r')
    set(gca,'YDir','normal')
    plot(x_cir,y_cir,'r')
    caxis([0 0.3])
    pause
end
%}
%%
% 6. Creat measurement matrix from the observation matrix (TTP matrix)
d = G*m0;   % ms

% 7 Inversion process to recover the ss aperturbation field
% travel time perturbation error
ntot =length(tx_lon);          % number of rays
Cd =diag(0.0002^2*ones(1,ntot))     ; % ms
% 1. generalized inverse matrig G^-g
G_geninv = P*transpose(G)*inv((G*P*transpose(G)+Cd));
% 2. recovered m
alpha0 = G_geninv*d;
alpha_1 = alpha0(1:tot_grid_num);
alpha_2 = alpha0(tot_grid_num+1:2*tot_grid_num);
alpha_3 = alpha0(2*tot_grid_num+1:3*tot_grid_num);
alpha_4 = alpha0(3*tot_grid_num+1:4*tot_grid_num);

m_recov= alpha0;
m_recov1 = alpha_1;
m_recov2= alpha_2;
m_recov3= alpha_3;
m_recov4 = alpha_4;

% Uncertianties 
prior_SD1 = reshape(sqrt(diag(P_mode1)),grid_num,grid_num);
prior_SD2 = reshape(sqrt(diag(P_mode2)),grid_num,grid_num);
prior_SD3 = reshape(sqrt(diag(P_mode3)),grid_num,grid_num);
prior_SD4 = reshape(sqrt(diag(P_mode4)),grid_num,grid_num);
prior_SD = sqrt(prior_SD1.^2+prior_SD2.^2+prior_SD3.^2+prior_SD4.^2);

P_post_alpha = (eye(tot_grid_num*4)-G_geninv*G)*P;
P_post_alpha1 = P_post_alpha(1:tot_grid_num,1:tot_grid_num);
P_post_alpha2 = P_post_alpha(tot_grid_num+1:2*tot_grid_num,tot_grid_num+1:2*tot_grid_num);
P_post_alpha3 = P_post_alpha(2*tot_grid_num+1:3*tot_grid_num,2*tot_grid_num+1:3*tot_grid_num);
P_post_alpha4 = P_post_alpha(3*tot_grid_num+1:4*tot_grid_num,3*tot_grid_num+1:4*tot_grid_num);
P_post_combined = P_post_alpha1.^2+P_post_alpha2.^2+P_post_alpha3.^2+P_post_alpha4.^2;

Res_mat = G_geninv*G;
Res_mat1 = Res_mat(1:tot_grid_num,1:tot_grid_num);
Res_mat2 = Res_mat(1:tot_grid_num,1:tot_grid_num);
model_res1 = reshape(diag(Res_mat1),grid_num,grid_num)';
% RMS Error
alpha_SD1 = reshape(sqrt(diag(P_post_alpha1)),grid_num,grid_num);
alpha_SD2 = reshape(sqrt(diag(P_post_alpha2)),grid_num,grid_num);
alpha_SD3 = reshape(sqrt(diag(P_post_alpha3)),grid_num,grid_num);
alpha_SD4 = reshape(sqrt(diag(P_post_alpha3)),grid_num,grid_num);
alpha_SD = sqrt(alpha_SD1.^2+alpha_SD2.^2+alpha_SD3.^2+alpha_SD4.^2);
sd_reduction1 = alpha_SD./prior_SD*100;

% Combine 4 modes
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
% 3 reshape the SS field
recov_SS1=reshape(SSP_d_avg2(1:end),grid_num,grid_num)'  ;
initial_SS=reshape(SSP_d_avg1(1:end),grid_num,grid_num)'  ;

% 4. Resolution Matrix

Res_mat_d_avg = eye(tot_grid_num)-P_SSP_d_avg*P_prior_d_avg^-1;


% 5. Posteriori covariance
prior_SD = reshape(sqrt(diag(P_prior_d_avg)),grid_num,grid_num);


% 6 residual
% d_recov = G*alpha0;
% residual = d-d_recov;

% 7 model resolution
model_res_d_avg = reshape(diag(Res_mat_d_avg),grid_num,grid_num)';


% variance reduction
post_SD_d_avg = reshape(sqrt(diag(P_SSP_d_avg)),grid_num,grid_num)';

sd_reduction2 = abs((post_SD_d_avg)./prior_SD*100);%./prior_SD1*100;


%% 8 plot the recovered ss pertrubation field
figure(2)
clf
set(gcf,'name','2-D sound speed perturbation field','Units','normalized','Position',[0.5 0.1 0.35 .4])
imagesc(x_cen,y_cen,initial_SS)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title('Initial Depth-Averaged SSP Field')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',12)
caxis([-.4 .4])

figure(3)
clf
set(gcf,'name','2-D sound speed perturbation field','Units','normalized','Position',[0 0.1 0.35 .4])
imagesc(x_cen,y_cen,recov_SS1)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title('Recovered Depth-Averaged SSP Field')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',12)
caxis([-.4 .4])

figure(4)
clf
set(gcf,'name','2-D sound speed perturbation field','Units','normalized','Position',[0 0.1 0.35 .4])
imagesc(x_cen,y_cen,recov_SS1-initial_SS)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title('Difference SSP Field')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',12)
caxis([-.4 .4])
% 

%% Solution Quality
ind1 = 1;
ind2 = 25;

%%% Plot
figure(5)
clf
set(gcf,'Units','normalized','Position',[.5 1 0.35 .4])
imagesc(x_cen(ind1:ind2),y_cen(ind1:ind2),real(sd_reduction1(ind1:ind2,ind1:ind2)))
cbar = colorbar;
% cbar.Label.String = 'RMS Error (m/s)';
cbar.Label.String = '%';
% caxis([0 .2])
title('RMS Error Reduction: Separate Modes')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',12)
colormap jet

figure(52)
clf
set(gcf,'Units','normalized','Position',[.1 1 0.35 .4])
imagesc(x_cen(ind1:ind2),y_cen(ind1:ind2),sd_reduction2(ind1:ind2,ind1:ind2))
cbar = colorbar;
% cbar.Label.String = 'RMS Error (m/s)';
cbar.Label.String = '%';
% caxis([0 .2])
title('RMS Error Reduction: Combined Modes')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',12)
colormap jet
%% Resolution Matrix
figure(61)
clf
set(gcf,'Units','normalized','Position',[.4 .3 0.4 .5])
imagesc(x_cen,y_cen,model_res1)
hold on
plot(x_cir,y_cir,'k')
axis tight
grid off
title('Model Resolution: Separate Modes (Mode 1)')
colorbar
caxis([0 .7])
% caxis([min(P_post(:)) max(P_post(:))])
colormap jet
set(gca,'YDir','normal','fontsize',12)

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
caxis([0 .7])
% caxis([min(P_post(:)) max(P_post(:))])
colormap jet
set(gca,'YDir','normal','fontsize',12)

% cd /Users/testuser/Documents/Oct2018Cruise/Figure/Inversion/2D
% s1 = sprintf('px_%i_ray_%i_P',grid_num,ntot);
% s2 = sprintf('px_%i_ray_%i',grid_num,ntot);
% saveas(f1,[s1 '.png'])
% saveas(f2,[s2 '.png'])
%% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P,var_c] = gaussian_cov_mx(x,y,mode)
% load singular vaues of each mode
cd /Users/testuser/Documents/MATLAB/Script/Data
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

function plot_integrand_G(intd_G_all1,intd_G_all2,intd_G_all3,intd_G_all4,intd_G_remat1,intd_G_remat2,intd_G_remat3,intd_G_remat4,z_cross,x_cen,y_cen,tx_lat,tx_lon,r_layer,xdist,grid_num,ACO_lon,ACO_lat,x_cir,y_cir)
indx_srt = ceil(grid_num*0.1);
indx_stp = grid_num - ceil(grid_num*0.1);

    figure(3)
    clf
    subplot(1,5,1)
    cla
    plot(intd_G_remat1(:),[ones(1,2)*z_cross(1) sort(repmat(z_cross(2:end-1),1,4)) ones(1,2)*z_cross(end)],'Color','b')
    hold on
    plot(intd_G_all1*120,linspace(z_cross(1),z_cross(end),length(intd_G_all1)),'--k')
    for k =2:length(z_cross)-1
        line([-1 1],ones(1,2)*z_cross(k),'Color','r')
    end
    grid on
    axis tight
    set(gca,'YDir','reverse')
    xlim([-1e-3 1e-3])
    title('Mode1')
    
    subplot(1,5,2)
    cla
    plot(intd_G_remat2(:),[ones(1,2)*z_cross(1) sort(repmat(z_cross(2:end-1),1,4)) ones(1,2)*z_cross(end)],'Color','b')
    hold on
    plot(intd_G_all2*120,linspace(z_cross(1),z_cross(end),length(intd_G_all2)),'--k')
    for k =2:length(z_cross)-1
        line([-1 1],ones(1,2)*z_cross(k),'Color','r')
    end
    grid on
    axis tight
    set(gca,'YDir','reverse')
    xlim([-1e-3 1e-3])
    title('Mode2')
    
    subplot(1,5,3)
    cla
    plot(intd_G_remat3(:),[ones(1,2)*z_cross(1) sort(repmat(z_cross(2:end-1),1,4)) ones(1,2)*z_cross(end)],'Color','b')
    hold on
    plot(intd_G_all3*120,linspace(z_cross(1),z_cross(end),length(intd_G_all3)),'--k')
    for k =2:length(z_cross)-1
        line([-1 1],ones(1,2)*z_cross(k),'Color','r')
    end
    grid on
    axis tight
    set(gca,'YDir','reverse')
    xlim([-1e-3 1e-3])
    title('Mode3')
    
    subplot(1,5,4)
    cla
    plot(intd_G_remat4(:),[ones(1,2)*z_cross(1) sort(repmat(z_cross(2:end-1),1,4)) ones(1,2)*z_cross(end)],'Color','b')
    hold on
    plot(intd_G_all4*120,linspace(z_cross(1),z_cross(end),length(intd_G_all4)),'--k')
    for k =2:length(z_cross)-1
        line([-1 1],ones(1,2)*z_cross(k),'Color','r')
    end
    grid on
    axis tight
    set(gca,'YDir','reverse')
    xlim([-1e-3 1e-3])
    title('Mode4')
    
    r_layer = real(repmat(cumsum(r_layer)',4,1));
    subplot(1,5,5)
    plot(r_layer(:)/1000,[ones(1,2)*z_cross(1) sort(repmat(z_cross(2:end-1),1,4)) ones(1,2)*z_cross(end)],'Color','b')
    grid on
    axis tight
    xlim([0 26])
    set(gca,'YDir','reverse')
    xlabel('r (km)')
    
    t = annotation('textbox',[.45 .9 .1 .1],'String',sprintf('%.0f',xdist)+" m",'Fontsize',12,'BackgroundColor','white');
    
    figure(44)
    rng(1)
    imagesc(x_cen(indx_srt:indx_stp),y_cen(indx_srt:indx_stp),0.5*rand(length(x_cen(indx_srt:indx_stp)),length(x_cen(indx_srt:indx_stp))))
    colormap gray
    hold on
    scatter(ACO_lon,ACO_lat,200,'pr','filled')
    scatter(tx_lon,tx_lat,'or','filled')
    line([ACO_lon tx_lon],[ACO_lat tx_lat],'color','r')
    set(gca,'YDir','normal')
    plot(x_cir,y_cir,'r')
    caxis([0 0.3])
    pause
end