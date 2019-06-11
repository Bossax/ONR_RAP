% 2-D  simulation
% simulate ocean sound speed perturbation field (m) with Gaussian model uncertainty
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
%%%% pixels do not have equal sizes
%%%%%%%%%%%%%%%%%%%%
%% 2. Generate Gaussian Covariance Matrix 
P = gaussian_cov_mx(x_cen,y_cen,1); % variance of sst is in m/s
%%
figure(1)
set(gcf,'name','Covariance Matrix','Units','normalized','Position',[0.6 .5 0.35 0.4])
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
P1=P;
while exist('D')==0
try
D=chol(P1,'lower');
end
P1=nearestSPD(P1);
end
rng(1)
rand_sst=normrnd(0,1,grid_num*grid_num,1);      %Random sound speed perturbation vector with normal distribution
m0=D*rand_sst;                                %Simulated SS field for each pixel

%% delta function at the center
mid_p = (tot_grid_num-1)/2;
m0 = [zeros(mid_p,1);3;zeros(mid_p,1)];
%% in the middle point of a radial line
pixel_num = 315;
m01 = [zeros(pixel_num,1);3;zeros(tot_grid_num-pixel_num-1,1)];

m0 = [m01];
%% Reshape
sim_SS=reshape(m0(1:end),grid_num,grid_num)'  ;% reshape works on a column basis need to swap row and column 
%% 4 Plot the Sound Speed Field
f2 = figure(2);
clf
set(gcf,'name','2-D sound speed perturbation field simulation','Units','normalized','Position',[0 1 0.35 0.4])
imagesc(x_cen,y_cen,sim_SS)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
caxis([-3 3])
title('Simulated Sound Speed Perturbation Field (30 km circle)')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
set(gca,'YDir','normal')
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

figure(2)
plot(x_cir,y_cir,'k')
%% 5 Create Ray Paths
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

figure(2)
for iii = 1:length(tx_lat)
    
    hold on
    plot([tx_lon(iii) ACO_lon],[tx_lat(iii) ACO_lat],'Linewidth',1,'Color','m')
    
end

%% concentric circles
R1 = 26000;
R2 = 20000;
% R3 = 15000;
% R4 = 10000;
% R5 = 5000;
R = [R1 R2];% R3 R4 R5];
spacing = 500;
n1 = floor(R1/500*2*pi);
n2 = floor(R2/500*2*pi);
% n3 = floor(R3/500*2*pi);
% n4 = floor(R4/500*2*pi);
% n5 = floor(R5/500*2*pi);
n = [n1 n2];% n3 n4 n5];
n_tot = sum(n);

figure(1)
angle = 0;
tx_ind = 1;
tx_lon = [];
tx_lat = [];
for iii = 1:length(n)
    for jjj = 1:n(iii)
            hold on
            angle = 360/n(iii)*(jjj-1);
            [tx_lon(tx_ind),tx_lat(tx_ind),~] = m_fdist(ACO_lon,ACO_lat,angle,R(iii));
            tx_lon(tx_ind) = tx_lon(tx_ind)-360;
            hold on
            %scatter(tx_lon(iii),tx_lat(iii),'po')
            plot([tx_lon(tx_ind) ACO_lon],[tx_lat(tx_ind) ACO_lat],'Linewidth',1,'Color','m')
            
        
        tx_ind =tx_ind+1;
    end
end


%% Obervation Matrix
G_prep = zeros(length(tx_lon),grid_num*grid_num);
for iii=1:length(tx_lon)
    iii
    % find pixels and the distance in each pixel
    [total_pixel_distance,total_pixel_num]=obs_matrix(tx_lat(iii),tx_lon(iii),ACO_lat,ACO_lon,x_node,y_node);
    
    % form matrix G
    g = zeros(1,grid_num*grid_num);
    for k = 1:length(total_pixel_num)
        G_prep(iii,total_pixel_num(k)) = total_pixel_distance(k);
     
    end
    total_pixel_num = total_pixel_num;
end 
priori_ss = 1506^2*ones(size(G_prep));
% matrix G
G = -G_prep.*(1./priori_ss);

%% plot elements og the observation matrix
figure(3)
plot(1:length(m0),sum((G),1))
xticks(1:100:length(m0))
grid on
% ylabel('meter')
xlabel('pixel')
xlim([1 length(m0)])
title('Obserbvation Matrix')

%% 6. Creat measurement matrix from the observation matrix (TTP matrix)
d = G*m0;   % ms

% 7 Inversion process to recover the ss aperturbation field
% travel time perturbation error
ntot =length(tx_lon);          % number of rays
Cd =diag(0.0002^2*ones(1,ntot))     ; % ms
% 1. generalized inverse matrig G^-g
G_geninv = P*transpose(G)*inv((G*P*transpose(G)+Cd));
% 2. recovered m
m_recov = G_geninv*d;
% 3. Resolution Matrix
Res_mat = G_geninv*G;
% 4. Posteriori covariance
P_post = (eye(grid_num*grid_num)-Res_mat)*P;

% 5 reshape the SS field
recov_SS=reshape(m_recov(1:end),grid_num,grid_num)'  ;
% 6 residual
residual = d-G*m_recov;
% 7 model resolution
prior_SD = reshape(sqrt(diag(P)),grid_num,grid_num);
model_res = reshape(diag(Res_mat),grid_num,grid_num);
post_SD = reshape(sqrt(diag(P_post)),grid_num,grid_num);

sd_reduction = post_SD./prior_SD*100;
% 8 plot the recovered ss pertrubation field

figure(3)
clf
imagesc(x_cen,y_cen,recov_SS)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title('Recovered Sound Speed Perturbation Field')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',12)
caxis([-6 6])
%% Difference
figure(4)
imagesc(x_cen,y_cen,sim_SS - recov_SS)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title('Difference between the actual SS field and the recovered SS field')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal')
 caxis([-1 1])
%%
figure(5)
imagesc(x_cen,y_cen,sd_reduction)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Reduction %';
title('RMS error reduction Uncertianty')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal')
set(gca,'fontsize',12)

figure(6)
imagesc(x_cen,y_cen,model_res)
axis tight
grid off
title('Model Resolution')
hold on
plot(x_cir,y_cir,'k')
colorbar
set(gca,'fontsize',12)
colormap jet

% cd /Users/testuser/Documents/Oct2018Cruise/Figure/Inversion/2D
% s1 = sprintf('px_%i_ray_%i_P',grid_num,ntot);
% s2 = sprintf('px_%i_ray_%i',grid_num,ntot);
% saveas(f1,[s1 '.png'])
% saveas(f2,[s2 '.png'])
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P,var_c] = gaussian_cov_mx(x,y,mode)
% load singular vaues of each mode
cd /Users/testuser/Documents/MATLAB/Script/Data
load EOF_SS.mat
var_c = (eval(['EOF_SS.singv' num2str(mode)]))^2;
% calculate distance between pixel
len = length(x);
range = zeros(len,len);


%X and Y matrix
x_matrix=repmat(x,1,len);
y_matrix=reshape(repmat(y,len,1),1,len*len);

% Gaussian cov length
sigma = 20000;

for ii=1:size(x_matrix,2) % pixel number
    for jj=1:size(x_matrix,2) % the second pixel
        
        %Distance between pixel centers
        pixel_distance=dist([y_matrix(ii) y_matrix(jj)],[x_matrix(ii) x_matrix(jj)]);
        
        %Gaussian covariance matrix
        P(ii,jj)=(var_c)*exp(-(pixel_distance^2)/(sigma*sigma));
        
    end
end
end