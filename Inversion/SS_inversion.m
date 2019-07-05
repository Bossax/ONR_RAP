% Inversion to calculate hydrophone position offsets and SS perturbation field
% Impose empirical orthogonal functions (EOF) on the vertical structure of ocean sound speed variation
% Specify the number of EOF modes included in the model
% SET-UP
% 1. hydrophone position
% 2. ocean domain
% 3. tx/rx files directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clearvars

% debug matrix G option
debug_obs = false; 

% spatial filter
spatial_filter = false;

%% 1. Initialize ocean domain
% Work on the spherical surface
%%%%%%%%%%%%%%%%
% 1. Ocean Domain
% ACO_lat = 22.738894;                  % original
% ACO_lon = -158.006009;                % original
 
ACO_lat= 22.738772;                  % June 2017
ACO_lon= -158.006186;                % June 2017

% ACO_lat= 22.738764;                  % June 2017 1st iteration
% ACO_lon= -158.0061781;               % June 2017 

ACO_depth =  -4729.92;                 % Original depth MSL
% ACO_depth = -4733.906;            % at 4,736.266 m June 2017
% ACO_depth = -4735.29;                 % June 2017 1st iteration


% Initialize the domain
L = 60000;      % meter EDIT
grid_num = 49;  % grid_num x grid_num  pixels EDIT
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
% 2.1 sound speed uncertainty
[P_mode1,var_c1] = gaussian_cov_mx(x_cen,y_cen,1); % priori uncerrtainty of mode 1
[P_mode2,var_c2] = gaussian_cov_mx(x_cen,y_cen,2); % priori uncerrtainty of mode 2
[P_mode3,var_c3] = gaussian_cov_mx(x_cen,y_cen,3); % priori uncerrtainty of mode 3
[P_mode4,var_c4] = gaussian_cov_mx(x_cen,y_cen,4); % priori uncerrtainty of mode 4

% 2.2 hydrophone uncertainty
var_x = 10; var_y = 10; var_z = 5;         % first round inverse
% var_x = 0.207; var_y = 0.2751; var_z = 0.0538;  % second round inverse
P_p = diag([var_x var_y var_z]);

% 2.3 bundle up P matrices
dimP = size(P_mode1);
zeroPc =zeros(dimP);
zeroPp_col = zeros(dimP(1),3);
zeroPp_row = zeros(3,dimP(1));
% P = [P_mode1 zeroPc zeroPc ; zeroPc P_mode2 zeroPc ;zeroPc zeroPc P_mode3];
P = [P_mode1 zeroPc zeroPc zeroPc zeroPp_col; zeroPc P_mode2 zeroPc zeroPc zeroPp_col;zeroPc zeroPc P_mode3 zeroPc zeroPp_col;zeroPc zeroPc zeroPc P_mode4 zeroPp_col;zeroPp_row zeroPp_row zeroPp_row zeroPp_row P_p];
 
%% Priori cov matrix
figure(1)
set(gcf,'name','Covariance Matrix','Units','normalized','Position',[0.3 .5 0.35 0.45])
imagesc(P)
axis tight
grid off
title('Gaussian Covariance Matrix')
colorbar
colormap jet


%% 3. Download Tx data points
% 3.1 download data files
month = 'Jun'
year = '2018'

switch year
    case '2017'
        %%% June 2017 %%% 
        day = 7:12;
        start_hour = 13;
        end_hour = 5;
        [tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,x_err,y_err,z_err,act_arrival,est_arrival,SNR] = tx_rx_extraction_June2017(day,start_hour,end_hour);

    case '2018'
        switch month
            case 'Jun'
                %%% June 2018
                day = 19:22;
                start_hour = 23;
                end_hour = 23;
                [tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,x_err,y_err,z_err,act_arrival,est_arrival,SNR] = tx_rx_extraction_June(day,start_hour,end_hour,'HEM');

                % remove outliers
                ttp = (act_arrival-real(est_arrival))*1000*3600*24;
                rm_ind = find(abs(ttp) > 4*(rms(ttp-median(ttp))));
                ttp(rm_ind) = [];
                tx_t(rm_ind) = [];
                tx_lon(rm_ind) = [];
                tx_lat(rm_ind) = [];
                tx_heading(rm_ind) = [];
                tx_altitude(rm_ind) = [];
                tx_xvel(rm_ind) = [];
                range(rm_ind) = [];
                x_err(rm_ind) = [];
                y_err(rm_ind) = [];
                z_err(rm_ind) = [];
                act_arrival(rm_ind) = [];
                est_arrival(rm_ind) = [];
                SNR(rm_ind) = [];

            case 'Oct'
                %%% October 2018
                day = 27:30;
                start_hour = 3;
                end_hour = 14;
                [tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,x_err,y_err,z_err,act_arrival,est_arrival,SNR] = tx_rx_extraction_Oct(day,start_hour,end_hour,'HEM');
                
        end
end

% limit range
keep_ind = find(range <15.2);

% 3.2 travel time perturbation
ttp_origin = (act_arrival - est_arrival)*3600*24*1000; 


ttp = ttp_origin(keep_ind);
tx_lat = tx_lat(keep_ind);
tx_lon = tx_lon(keep_ind);
tx_altitude = tx_altitude(keep_ind);
tx_heading = tx_heading(keep_ind);
x_err = x_err(keep_ind);
y_err = y_err(keep_ind);
z_err = z_err(keep_ind);
SNR = SNR(keep_ind);
range = range(keep_ind);

% 3.3 subtract ship factors 
%{
% relative heading
azmth = ones(1,length(tx_lat));

for ii=1:length(tx_lat)
    azmth(ii) = azimuth(ACO_lat,ACO_lon,tx_lat(ii),tx_lon(ii));
end
theta = tx_heading - azmth;

for ii = 1:length(theta)
   if theta(ii) < 0 
       theta(ii) = theta(ii)+360;
   end
end
% Empirical function (from Lasso regression)
% Function = 0.4619*range;
Function = 0.4458*range;

% TTP residueal
ttp = ttp_origin - Function;
%}

%% 4. draw circle
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

% plot
figure(1)
clf
imagesc(x_cen,y_cen,rand(length(x_cen),length(x_cen)))
colormap gray
hold on
scatter(ACO_lon,ACO_lat,200,'pr','filled')
set(gca,'YDir','normal')
plot(x_cir,y_cir,'r')

%% 5. Obervation Matrix
G1 = zeros(length(tx_lon),grid_num*grid_num);
G2 = zeros(length(tx_lon),grid_num*grid_num);
G3 = zeros(length(tx_lon),grid_num*grid_num);
G4 = zeros(length(tx_lon),grid_num*grid_num);
Ghyd =zeros(length(tx_lon),3);
M = zeros(length(tx_lon),1);
theta_r = zeros(1,length(tx_lon));  % received ray angle
theta_l = zeros(1,length(tx_lon));  % launched ray angle

% loop over each ray
for iii=1:length(tx_lon)
    iii
    % find pixels and the distance in each pixel
    
    [G_of_nray1,total_pixel_num,theta_l(iii),theta_r(iii),intd_G1,z_cross,intd_G_all1,z,SS,~]=obs_matrix3D_v2(tx_lat(iii),tx_lon(iii),tx_altitude(iii),ACO_lat,ACO_lon,ACO_depth,x_node,y_node,1,month,year);
    [G_of_nray2,~,~,~,intd_G2,~,intd_G_all2,~,~,~]=obs_matrix3D_v2(tx_lat(iii),tx_lon(iii),tx_altitude(iii),ACO_lat,ACO_lon,ACO_depth,x_node,y_node,2,month,year);
    [G_of_nray3,~,~,~,intd_G3,~,intd_G_all3,~,~,~]=obs_matrix3D_v2(tx_lat(iii),tx_lon(iii),tx_altitude(iii),ACO_lat,ACO_lon,ACO_depth,x_node,y_node,3,month,year);
    [G_of_nray4,~,~,~,intd_G4,~,intd_G_all4,~,~,~]=obs_matrix3D_v2(tx_lat(iii),tx_lon(iii),tx_altitude(iii),ACO_lat,ACO_lon,ACO_depth,x_node,y_node,4,month,year);
    [G_of_hyd] = obs_matrix_hyd(tx_lat(iii),tx_lon(iii),tx_altitude(iii),theta_r(iii),ACO_lat,ACO_lon,ACO_depth);
   
    % form matrix G
    G1(iii,:) = G_of_nray1;
    G2(iii,:) = G_of_nray2;
    G3(iii,:) = G_of_nray3;
    G4(iii,:) = G_of_nray4;
    Ghyd(iii,:) = G_of_hyd;
    M(iii) = length(total_pixel_num)-1;
    
end 

cl = SS(1);  cr = SS(end);  
G1 = real(G1);
G2 = real(G2);
G3 = real(G3);
G4 = real(G4);
Ghyd = real(Ghyd);

G = [G1 G2 G3 G4 Ghyd];

%% plot elements of the observation matrix
if debug_obs
    
    figure(4)
    subplot(4,1,1)
    plot(1:tot_grid_num,sum((G1),1))
    xticks(1:50:tot_grid_num)
    grid on
    ylabel('meter')
    xlabel('pixel')
    xlim([1 tot_grid_num])
    title('Obserbvation Matrix of the first mode')
    
    subplot(4,1,2)
    plot(1:tot_grid_num,sum((G2),1))
    xticks(1:50:tot_grid_num)
    grid on
    ylabel('meter')
    xlabel('pixel')
    xlim([1 tot_grid_num])
    title('Obserbvation Matrix of the second mode')
    
    subplot(4,1,3)
    plot(1:tot_grid_num,sum((G3),1))
    xticks(1:50:tot_grid_num)
    grid on
    ylabel('meter')
    xlabel('pixel')
    xlim([1 tot_grid_num])
    title('Obserbvation Matrix of the third mode')
    
    subplot(4,1,4)
    plot(1:tot_grid_num,sum((G4),1))
    xticks(1:50:tot_grid_num)
    grid on
    ylabel('meter')
    xlabel('pixel')
    xlim([1 tot_grid_num])
    title('Obserbvation Matrix of the third mode')
    
    
    figure
    plot(1:3,sum((Ghyd),1))
    % xticks(1:10:tot_grid_num)
    grid on
    ylabel('meter')
    xlabel('pixel')
    % xlim([1 tot_grid_num])
    title('Obserbvation Matrix of the hydrophne position offsets')


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inversion %%
% 6. Inversion process to recover the ss aperturbation field
% measurement matrix (TTP matrix)
d = real(ttp'/1000);   

% travel time perturbation error
Cd = data_error_matrix(x_err,y_err,z_err,tx_lat,tx_lon,tx_heading,theta_l,theta_r,SNR,cl,cr,ACO_lat,ACO_lon,ACO_depth);
Cd = real(diag(Cd));



% % % % % % form a filter matrix % % % % % % % 
if spatial_filter
    W = zeros(length(d),length(d));
    for iii = 1:length(d)
        iii
        
        tx_dist = distance(tx_lat(iii)*ones(1,length(d)-iii+1),tx_lon(iii)*ones(1,length(d)-iii+1),tx_lat(iii:end),tx_lon(iii:end),referenceEllipsoid('WGS84'));
        W(iii,iii:end) = exp(-(tx_dist/20000).^2);
        
    end
    
    W = W+W'-diag(diag(W));
end
% % % % % % %  pre-spatially filter d (before inversion)% % % % % % % 
% for spatial average
%{
d_f = (W*d)./sum(W,2);
figure(76)
scatter(tx_lon,tx_lat,20,d_f*1000,'filled')
grid on
xlabel('Long')
ylabel('Lat')
grid on
colormap jet
cbar = colorbar;
cbar.Label.String = 'Travel Time Perturbation (ms)';
title('Spatially-filtered TTP')
% caxis([.2 1.4])
colormap jet
set(gca,'fontsize',13)
axis tight

d = d_f;
%}
%% Calculation
% 6.1. generalized inverse matrig G^-g
inv_U = (G*P*G'+Cd)^-1;
G_geninv = P*G'*(G*P*G'+Cd)^-1;

% 6.2. recovered m
m_recov = G_geninv*d;

% performance metrics
Res_mat = G_geninv*G;                                                       % resolution 
P_post = round((eye(tot_grid_num*4+3)-G_geninv*G)*P,10,'significant');      % Posteriori uncertainty
P_post_triu = triu(P_post);
P_post_tril = tril(P_post);
P_post_newu = (P_post_triu+P_post_triu')-diag(diag(P_post_triu));
P_post_newl = (P_post_tril+P_post_tril')-diag(diag(P_post_tril));
P_post_new = (P_post_newu+P_post_newl)/2;


%% Arrange matrices
% recovered m
m_recov1 = m_recov(1:tot_grid_num);
m_recov2 = m_recov(tot_grid_num+1:2*tot_grid_num);
m_recov3 = m_recov(2*tot_grid_num+1:3*tot_grid_num);
m_recov4 = m_recov(3*tot_grid_num+1:4*tot_grid_num);
m_p = m_recov(4*tot_grid_num+1:4*tot_grid_num+3);
% m_p = [0;0;0];

% 6.3 for depth-average : alpha = m^
alpha0 = m_recov(1:end-3);
alpha_1 = alpha0(1:tot_grid_num);
alpha_2 = alpha0(tot_grid_num+1:2*tot_grid_num);
alpha_3 = alpha0(2*tot_grid_num+1:3*tot_grid_num);
alpha_4 = alpha0(3*tot_grid_num+1:4*tot_grid_num);

% 6.4 Resolution Matrix
Res_mat1 = Res_mat(1:tot_grid_num,1:tot_grid_num);
Res_mat2 = Res_mat(tot_grid_num+1:2*tot_grid_num,tot_grid_num+1:2*tot_grid_num);
Res_mat3 = Res_mat(2*tot_grid_num+1:3*tot_grid_num,2*tot_grid_num+1:3*tot_grid_num);
Res_mat4 = Res_mat(3*tot_grid_num+1:4*tot_grid_num,3*tot_grid_num+1:4*tot_grid_num);
Res_mat_p  = Res_mat(4*tot_grid_num+1:4*tot_grid_num+3,4*tot_grid_num+1:4*tot_grid_num+3);

% 6.5 Posteriori covariance
% Posteriori covariance
P_post1 = P_post_new(1:tot_grid_num,1:tot_grid_num);
P_post2 = P_post_new(tot_grid_num+1:2*tot_grid_num,tot_grid_num+1:2*tot_grid_num);
P_post3 = P_post_new(2*tot_grid_num+1:3*tot_grid_num,2*tot_grid_num+1:3*tot_grid_num);
P_post4 = P_post_new(3*tot_grid_num+1:4*tot_grid_num,3*tot_grid_num+1:4*tot_grid_num);
P_post_p = P_post_new(4*tot_grid_num+1:4*tot_grid_num+3,4*tot_grid_num+1:4*tot_grid_num+3);

% 6.6 reshape the SS field
recov_SS1=reshape(m_recov1(1:end),grid_num,grid_num)'  ;
recov_SS2=reshape(m_recov2(1:end),grid_num,grid_num)'  ;
recov_SS3=reshape(m_recov3(1:end),grid_num,grid_num)'  ;
recov_SS4=reshape(m_recov4(1:end),grid_num,grid_num)'  ;

% 6.7 residual
residual = d-G*m_recov;
norm_res = residual./sqrt(diag(Cd));
d_recov_ocean = G(:,1:end-3)*m_recov(1:end-3);

% 6.8 model resolution
model_res1 = reshape(diag(Res_mat1),grid_num,grid_num)';
model_res2 = reshape(diag(Res_mat2),grid_num,grid_num)';
model_res3 = reshape(diag(Res_mat3),grid_num,grid_num)';
model_res4 = reshape(diag(Res_mat4),grid_num,grid_num)';
model_res_p = reshape(diag(Res_mat_p),3,1);

% 6.9 variance reduction
prior_SD1 = reshape(sqrt(diag(P_mode1)),grid_num,grid_num);
prior_SD2 = reshape(sqrt(diag(P_mode2)),grid_num,grid_num);
prior_SD3 = reshape(sqrt(diag(P_mode3)),grid_num,grid_num);
prior_SD4 = reshape(sqrt(diag(P_mode4)),grid_num,grid_num);
prior_SD = sqrt(prior_SD1.^2+prior_SD2.^2+prior_SD3.^2+prior_SD4.^2);
prior_p = reshape(sqrt(diag(P_p)),3,1);

% 6.10 model uncertainty
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
post_p = reshape(sqrt(diag(P_post_p)),3,1);

sd_reduction1 = (post_SD1)./prior_SD1*100;
sd_reduction2 = (post_SD2)./prior_SD2*100;
sd_reduction3 = (post_SD3)./prior_SD3*100;
sd_reduction4 = (post_SD4)./prior_SD4*100;
sd_reduction_p = post_p./prior_p*100;

%%%%% Depth Averaing %%%%%%%%%
% Combine 4 modes
load('EOF_SS.mat')
f1 = EOF_SS.mode1;
f2 = EOF_SS.mode2;
f3 = EOF_SS.mode3;
f4 = EOF_SS.mode4;
ll=length(z)-length(f1);
% ll = 2408 - length(f1);
f1 = [f1;f1(end)*ones(ll,1)];
f2 = [f2;f2(end)*ones(ll,1)];
f3 = [f3;f3(end)*ones(ll,1)];
f4 = [f4;f4(end)*ones(ll,1)];
F = [f1 f2 f3 f4];
alpha = [alpha_1';alpha_2';alpha_3';alpha_4'];
SSP = F*alpha;

% depth averaged modes
f1_avg = sum(f1(1:end-1).*(z(2:end)-z(1:end-1)))/z(end);
f2_avg = sum(f2(1:end-1).*(z(2:end)-z(1:end-1)))/z(end);
f3_avg = sum(f3(1:end-1).*(z(2:end)-z(1:end-1)))/z(end);
f4_avg = sum(f4(1:end-1).*(z(2:end)-z(1:end-1)))/z(end);

% depth averaging and mode combining operator
P_post_alpha = P_post_new;

F_operator = [diag(ones(1,grid_num^2)*f1_avg) diag(ones(1,grid_num^2)*f2_avg) diag(ones(1,grid_num^2)*f3_avg) diag(ones(1,grid_num^2)*f4_avg)];

SSP_d_avg2 = F_operator*alpha0;
P_prior_d_avg = F_operator*P(1:end-3,1:end-3)*F_operator';
P_SSP_d_avg = F_operator*P_post_alpha(1:end-3,1:end-3)*F_operator';


% recalculate the SSP field
% reshape the SS field
recov_SS_d_avg=reshape(SSP_d_avg2(1:end),grid_num,grid_num)'  ;
% initial_SS_d_avg=reshape(SSP_d_avg1(1:end),grid_num,grid_num)'  ;

% Resolution Matrix
Res_mat_d_avg = eye(tot_grid_num)-P_SSP_d_avg*P_prior_d_avg^-1;


% 5. Posteriori covariance
prior_SD_d_avg = reshape(sqrt(diag(P_prior_d_avg)),grid_num,grid_num);

% model resolution
model_res_d_avg = reshape(diag(Res_mat_d_avg),grid_num,grid_num)';


% variance reduction
post_SD_d_avg = reshape(sqrt(diag(P_SSP_d_avg)),grid_num,grid_num)';
sd_reduction_d_avg = abs((post_SD_d_avg)./prior_SD_d_avg*100);%./prior_SD1*100;


%%%%%%%%%%%%%%%%  Spatially filtering the measurements  %%%%%%%%%%%%%%%%
%% for spatial averaging
%{
d_new = (W*d_recov_ocean)./sum(W,2);
%
figure(21)
scatter(tx_lon,tx_lat,20,d_new*1000,'filled')
grid on
xlabel('Long')
ylabel('Lat')
grid on
colormap jet
cbar = colorbar;
cbar.Label.String = 'Travel Time Perturbation (ms)';
title('Spatially-filtered TTP')
% caxis([1 2])
colormap jet
set(gca,'fontsize',13)
axis tight

%%%%% Repeat the calculation
G_m = G(:,1:end-3);
P_m = P(1:end-3,1:end-3);
G_geninv = P_m*G_m'*inv((G_m*P_m*G_m'+Cd));

% recovered m
m_recov2 = G_geninv*d_new;

% Arrange matrices
% recovered m
m_recov21 = m_recov2(1:tot_grid_num);
m_recov22 = m_recov2(tot_grid_num+1:2*tot_grid_num);
m_recov23 = m_recov2(2*tot_grid_num+1:3*tot_grid_num);
m_recov24 = m_recov2(3*tot_grid_num+1:4*tot_grid_num);


% for depth-average
alpha20 = m_recov2;
alpha_21 = alpha20(1:tot_grid_num);
alpha_22 = alpha20(tot_grid_num+1:2*tot_grid_num);
alpha_23 = alpha20(2*tot_grid_num+1:3*tot_grid_num);
alpha_24 = alpha20(3*tot_grid_num+1:4*tot_grid_num);

% reshape the SS field
recov_SS21=reshape(m_recov21(1:end),grid_num,grid_num)'  ;
recov_SS22=reshape(m_recov22(1:end),grid_num,grid_num)'  ;
recov_SS23=reshape(m_recov23(1:end),grid_num,grid_num)'  ;
recov_SS24=reshape(m_recov24(1:end),grid_num,grid_num)'  ;

%  residual
residual2 = d_new-G_m*m_recov2;
norm_res2 = residual2./sqrt(diag(Cd));
d_recov_ocean2 = G_m*m_recov2;

% Combine 4 modes
SSP_d_avg22 = F_operator*alpha20;

% recalculate the SSP field
% reshape the SS field
recov_SS_d_avg2 = reshape(SSP_d_avg22(1:end),grid_num,grid_num)'  ;
% initial_SS_d_avg=reshape(SSP_d_avg1(1:end),grid_num,grid_num)'  ;
%}
%% 7. plot the recovered ss pertrubation field
% draw circle
R = 15000;

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

% Mode-separate solution
figure(2)
clf
set(gcf,'Units','normalized','Position',[0. 0.9 .5 0.7]);
ax1 = axes;
set(ax1,'Units','normalized','Position',[0.05 0.5 .4 0.35]);
imagesc(x_cen,y_cen,recov_SS1)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
caxis([-5 5])
title('Recovered SS Perturbation Field (1st mode)')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal')
colormap jet

ax2 = axes;
set(ax2,'Units','normalized','Position',[0.55 0.5 .4 0.35]);
imagesc(x_cen,y_cen,recov_SS2)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title('Recovered SS Perturbation Field (2nd mode) ')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal')
caxis([-5 5])
colormap jet

ACO_pos = sprintf('Lat = %.5f\nLon = %.5f\n Depth =%.2f m',ACO_lat,ACO_lon,ACO_depth);
hyd_off = sprintf('HEM: Lat = %.5f Lon = %.5f Depth =%.2f m\nHydrophone Offset: x = %.3f m, y = %.3f m, z = %.3f m',ACO_lat,ACO_lon,ACO_depth,m_p(1),m_p(2),m_p(3));
t = annotation('textbox',[.25 .88 .1 .1],'String',hyd_off,'Fontsize',12,'BackgroundColor','white');

ax3 = axes;
set(ax3,'Units','normalized','Position',[0.05 0.05 .4 0.35]);
imagesc(x_cen,y_cen,recov_SS3)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title('Recovered SS Perturbation Field (3rd mode) ')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal')
caxis([-5 5])
colormap jet


ax4 = axes;
set(ax4,'Units','normalized','Position',[0.55 0.05 .4 0.35]);
imagesc(x_cen,y_cen,recov_SS4)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title('Recovered SS Perturbation Field (4th mode) ')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal')
caxis([-5 5])
colormap jet
%% Depth-Averaged solution
figure(23)
clf
set(gcf,'name','2-D sound speed perturbation field','Units','normalized','Position',[0 0.1 0.34 .4])
imagesc(x_cen,y_cen,recov_SS_d_avg)
colormap jet
cbar = colorbar;
cbar.Label.String = 'Sound Speed Perturbation (m/s)';
title('Depth-Averaged SSP Field: No Filtering (HEM)')
hold on
scatter(ACO_lon,ACO_lat,200,'pk','filled')
plot(x_cir,y_cir,'k')
set(gca,'YDir','normal','fontsize',12)
peak = 0.6
% caxis([-0.3-peak -0.3+peak])
caxis([-.5 .5])
xlabel('Long')
ylabel('Lat')


%% 8. Solution Quality
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
title('RMS Error Reduction: Depth-Averaged Solution')
hold on
% scatter(ACO_lon,ACO_lat,200,'pk','filled')
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
set(gcf,'Units','normalized','Position',[.4 .3 0.5 .55])
subplot(2,2,1)
imagesc(x_cen,y_cen,model_res1)
hold on
plot(x_cir,y_cir,'k')
axis tight
grid off
title('Model Resolution 1st mode')
colorbar
set(gca,'YDir','normal','fontsize',12)
xlabel('Long')
ylabel('Lat')
caxis([0 .7])
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
caxis([0 .7])
% caxis([min(P_post(:)) max(P_post(:))])
colormap jet
set(gca,'YDir','normal','fontsize',12)
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
caxis([0 .7])
% caxis([min(P_post(:)) max(P_post(:))])
colormap jet
set(gca,'YDir','normal','fontsize',12)

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
set(gca,'YDir','normal','fontsize',12)
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
title('Model Resolution')
colorbar
caxis([0 1])
xlabel('Lon')
ylabel('Lat')
% caxis([min(P_post(:)) max(P_post(:))])
colormap jet
set(gca,'YDir','normal','fontsize',12)

%% residual 
figure(5)
clf
% scatter(range,norm_res,'.')
histogram(norm_res,40)
grid on
ylabel('Frequency')
xlabel('Normalized Residual')
rms_res = rms(norm_res-median(norm_res));
% yticks(-4:4)
xlim([-3 3])
good_data = length(find(abs(norm_res)<=1));
percent_g_data =  good_data/length(norm_res)*100;
title(sprintf('Normalized Residual Median: %.3f, RMS: %.3f',median(residual),rms_res))
set(gca,'fontsize',12)
annotation('textbox',[.6 .8 0.1 .1],'String',[sprintf('%.2f',percent_g_data) '% of data is within +/-1'],'fontsize',11,'background','white');
line([1 1],[0 1500],'color','red')
line([-1 -1],[0 1500],'color','red')
ylim([0 1500])
%% Histogram
edges = -10.1:.2:10.1;
counts = -10:0.2:10;
figure(6)
subplot(2,1,1)
histogram(d*1000,'BinEdges',edges)
xlim([-6 6])
ylim([0 1500])
grid on
med_d = median(d*1000);
rms_d = rms(d*1000-med_d);
title(sprintf('Measurement: Median = %.2f ms, RMS = %.2f ms',med_d,rms_d))
ylabel('Frequency')
xlabel('msec')
subplot(2,1,2)
histogram(d_recov_ocean*1000,'BinEdges',edges)
grid on
xlim([-6 6])
ylim([0 1500])
med_drecov = median(d_recov_ocean*1000);
rms_drecov = rms(d_recov_ocean*1000-med_drecov);
title(sprintf('Recovered Measurement: Median = %.2f ms, RMS = %.2f ms',med_drecov,rms_drecov))
ylabel('Frequency')
xlabel('msec')

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
title('Resolution versus Range (Covariance Length = 8km)')
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
title({'Resolution versus Range: Depth-averaged solution',' (Correlation Length = 8 km)'})
set(gca,'fontsize',14)
line([0 50],[0 0],'color','k')
ylim([0 1])
 %% 9. Save file
 cd /Users/testuser/Documents/ONR_RAP/Data/Inversion_file/June2018
%  save HEM_inverse_solution_June2018 ACO_lat ACO_lon ACO_depth G G_geninv d Cd P P_mode1 P_mode2 P_mode3 P_mode4 P_p P_prior_d_avg P_post P_SSP_d_avg Res_mat Res_mat_d_avg SSP_d_avg2 m_recov x_cen y_cen z W 
 save HEM_inverse_solution_June2018_originaldepth_L20km_sizing_2_txuncer_newrx  ACO_lat ACO_lon ACO_depth G G_geninv d Cd P P_mode1 P_mode2 P_mode3 P_mode4 P_p P_prior_d_avg P_post_p P_post_new P_SSP_d_avg Res_mat Res_mat_d_avg SSP_d_avg2 m_recov x_cen y_cen z tx_lat tx_lon tx_heading range tot_grid_num
%% 10. write hydrophone position offsets to spreadsheet
% filename = 'Hydrophone_Position.xlsx';
% 
% % sheet
% switch year
%     case '2017'
%         %%% June 2017 %%% 
%         sheet = 'June2017_HEM'
% end
%% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P,var_c] = gaussian_cov_mx(x,y,mode)
% load singular vaues of each mode
cd /Users/testuser/Documents/ONR_RAP/Data/inversion_file
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
% sigma = 8000;   % correlation length from Jun3 2017 data (Variogram result)

for ii=1:size(x_matrix,2) % pixel number
    for jj=ii:size(x_matrix,2) % the second pixel
        
        %Distance between pixel centers
        pixel_distance=dist([y_matrix(ii) y_matrix(jj)],[x_matrix(ii) x_matrix(jj)]);
        
        %Gaussian covariance matrix
        P(ii,jj)=(var_c)*exp(-(pixel_distance^2)/(sigma*sigma));
        
    end
end
P = P+P'-diag(diag(P));
end

