function [obs_ray,total_pixel_num,theta_l,theta_r,intd_G,z_cross,intd_G_all,z,SS_priori,r_layer]=obs_matrix(tx_lat,tx_lon,tx_altitude,ACO_lat,ACO_lon,ACO_depth,x_tick,y_tick,mode,type,EOF_plot)
%  oberservation_matrix construction
% 1. Find the surface distance of an acoustic ray travels within one pixel
% 2. match the surface distance with the ray tracing calcluation to divide the ray path into pixel-wise arc lengths
% 3. y_tick = y coordinates of pixel boundaries from higher latitudes to  lower latitudes
% 4. the last section is for model representation analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if the function is called for a real inversion or a simulation case
if ~exist('type','var')
    % default is a real inversion
      type = "inversion";
end

if ~exist('EOF_plot','var')
    % default is a real inversion
      EOF_plot = 0;
end
type = string(type);

% Load EOF 
load('EOF_SS.mat')
EOF = eval(['EOF_SS.mode' num2str(mode)]);
z_EOF = EOF_SS.z;

domain_size = length(x_tick)-1;
stop_y=ACO_lat;
stop_x=ACO_lon;
start_y=tx_lat;
start_x=tx_lon;

% creates 10,000 points in the line between tx and ACO
res_dist = .2;  % horizontal distance resolution (larger than the minimum number a ray tarvels in one 2-m thick layer)
[tot_sep,az] = distance(start_y,start_x,stop_y,stop_x,referenceEllipsoid('WGS84'));
sep_point = floor(tot_sep/res_dist);
range = res_dist*(0:sep_point+1);
range(end) = tot_sep;
% lon_cur/lat_cur increases from tx point to ACO
[lon_cur,lat_cur,~] = m_fdist(start_x,start_y,az,range);
lon_cur = lon_cur-360;

%%% ray crosses longitude lines %%%%
%%%%%% X-INT %%%%%%%
% find xpos, boundary grids/ x_int, crossed lines
% tx point is on the west
if start_x<stop_x
    % find the range of X
    x_pos1=min(x_tick(start_x <= x_tick));    % lower limit
    x_pos2=max(x_tick(x_tick <= stop_x));     %  upper limit
    
    % points in between
    if size(x_pos2,2)~=0 && size(x_pos1,2)~=0
        x_int=x_tick(x_tick>=x_pos1 & x_tick<=x_pos2);
    elseif size(x_pos1,2)~=0
        x_int=x_pos1;
    elseif size(x_pos2,2)~=0
        x_int=x_pos2;
    end
    
% tx point is on the east
elseif start_x>stop_x
    x_pos1=max(x_tick(x_tick < start_x));
    x_pos2=min(x_tick(stop_x <x_tick ));
    if size(x_pos2,2)~=0 && size(x_pos1,2)~=0
        x_int=fliplr(x_tick(x_tick<=x_pos1 & x_tick>=x_pos2));
    elseif size(x_pos1,2)~=0
        x_int=x_pos1;
    elseif size(x_pos2,2)~=0
        x_int=x_pos2;
        
    end
elseif start_x == stop_x
        x_int = [];
end
    
% calculate surface distances of points in the interval
pixel_num_x =[];
pixel_distance_x = [];
if size(x_int,2)>=1
    for k=1:size(x_int,2)
        % find the nearest lon_cur to the x_interval point from the tx location side
        if start_x<stop_x
            dir = 'west';
            dist_pos = find(lon_cur >= x_int(k),1,'first'); % the first point on the right of the x_int line
        elseif start_x>stop_x
            dir = 'east';
            dist_pos = find(lon_cur >= x_int(k),1,'last'); % the first point on the right    
        end
        % Assign lat/lon to pixel boundary handles
        lon_pos(k)=lon_cur(dist_pos);
        lat_pos(k)=lat_cur(dist_pos);
        pixel_distance_x(k)=range(dist_pos);
        x_x = [];
        x_y = [];
        % Loop over all pixels to assign pixel number to the segment
        % find the crossed lines
        for ii=1:length(y_tick)-1
            for jj=1:length(x_tick)-1
                % check if the lon/lat is within that pixel
                if (x_tick(jj)<lon_pos(k))&&(lon_pos(k)<x_tick(jj+1)) && (y_tick(ii+1)<=lat_pos(k))&&(lat_pos(k)<=y_tick(ii))
                    % pixel number
%                     fprintf(sprintf('x = %i \n',jj)) 
                    x_x = jj;
                    x_y = ii;
                    
                end
                
            end
        end
        switch dir
            case 'west'
                pixel=(x_x-1)+(x_y-1)*(domain_size);
            case 'east'
                pixel=(x_x)+(x_y-1)*(domain_size);
        end
        pixel_num_x(k)=pixel;
        
    end
end

%%% ray crosses latitude lines %%%%
%%%%%%%Y-INT%%%%%%%%%%%%5
% tx lat is in the south
if start_y<stop_y
    % find the range of Y
        y_pos1=min(y_tick(y_tick > start_y));
        y_pos2=max(y_tick(y_tick < stop_y));
        
        % find points in between
        if size(y_pos2,2)~=0 && size(y_pos1,2)~=0
            y_int=fliplr(y_tick(y_tick>=y_pos1 & y_tick<=y_pos2));
        elseif size(y_pos1,2)~=0
            y_int=y_pos1;
        elseif size(x_pos2,2)~=0
            y_int=y_pos2;
        end
        
    % tx lat is in the north
    elseif start_y>stop_y
        % find the range of Y
        y_pos1=max(y_tick(y_tick < start_y));
        y_pos2=min(y_tick(y_tick > stop_y));
        % find points in between
        if size(y_pos2,2)~=0 && size(y_pos1,2)~=0
            y_int=(y_tick(y_tick<=y_pos1 & y_tick>=y_pos2));
        elseif size(y_pos1,2)~=0
            y_int=y_pos1;
        elseif size(x_pos2,2)~=0
            y_int=y_pos2;
        end
        
        % tx lat is at the same lat of the ACO
    elseif start_y==stop_y
        y_int=[];
end
    
    
    lat_pos = [];
    lon_pos = [];
    pixel_num_y =[];
    pixel_distance_y = [];
    if size(y_int,2)>=1
        for k=1:length(y_int)
            % find the nearest lat_cur to the y_interval point from the tx location side
            % in the south
           if start_y<stop_y
            dir = 'south';
            dist_pos = find(lat_cur <= y_int(k),1,'last'); % the first point in the south
        elseif start_y>stop_y
            dir = 'north';
            dist_pos = find(lat_cur <= y_int(k),1,'first'); % the first point in the south    
        end
        % Assign lat/lon to pixel boundary handles
        lon_pos(k)=lon_cur(dist_pos);
        lat_pos(k)=lat_cur(dist_pos);
        pixel_distance_y(k)=range(dist_pos);
        x_x = [];
        x_y = [];
        % Loop over all pixels to assign pixel number to the segment
        % find the crossed lines
        for ii=1:length(y_tick)-1
            for jj=1:length(x_tick)-1
                % check if the lon/lat is within that pixel
                if (x_tick(jj)<=lon_pos(k))&&(lon_pos(k)<=x_tick(jj+1)) && (y_tick(ii+1)<lat_pos(k))&&(lat_pos(k)<y_tick(ii))
                    % pixel number
%                     fprintf(sprintf('y = %i \n',ii)) 
                    x_x = jj;
                    x_y = ii;
                    
                end
                
            end
        end
        switch dir
            case 'south'
                pixel=(x_x)+(x_y-1)*(domain_size);
            case 'north'
                pixel=(x_x)+(x_y-2)*(domain_size);
        end
        pixel_num_y(k)=pixel;
        end
    end
    
    
    
    %Final Pixel
    pixel_num_last = [];
    pixel_distance_final = [];
    % ACO is not located at any node
    if isempty(intersect(x_tick,ACO_lon)) && isempty(intersect(y_tick,ACO_lat))
        pixel_distance_final=range(end);
        for ii=1:length(y_tick)-1
            for jj=1:length(x_tick)-1
                
                if (stop_x<x_tick(jj+1)) && (stop_x>x_tick(jj)) && (stop_y<y_tick(ii)) && (stop_y>y_tick(ii+1))
%                     fprintf(sprintf('last =  %i \n',ii))
                    pixel=(jj)+(ii-1)*domain_size;
                    
                end
                
            end
        end
        pixel_num_last=pixel(end);
    end
    
%%% sort pixel nunmbers and surface distances
if length(pixel_num_x)~=0 && length(pixel_num_y)~=0
        cum_pixel_distance=horzcat(pixel_distance_x,pixel_distance_y,pixel_distance_final);
        [cum_pixel_distance,b]=sort(cum_pixel_distance);
        total_pixel_num=horzcat(pixel_num_x,pixel_num_y,pixel_num_last);
        total_pixel_num=total_pixel_num(b);
    elseif length(pixel_num_x)~=0
        cum_pixel_distance=horzcat(pixel_distance_x,pixel_distance_final);
        [cum_pixel_distance,b]=sort(cum_pixel_distance);
        total_pixel_num=horzcat(pixel_num_x,pixel_num_last);
        total_pixel_num=total_pixel_num(b);
    elseif length(pixel_num_y)~=0
        cum_pixel_distance=horzcat(pixel_distance_y,pixel_distance_final);
        [cum_pixel_distance,b]=sort(cum_pixel_distance);
        total_pixel_num=horzcat(pixel_num_y,pixel_num_last);
        total_pixel_num=total_pixel_num(b);
    else
        cum_pixel_distance=pixel_distance_final;
        total_pixel_num=pixel_num_last;
    end
    
% truncate the total distance into sections
    if length(pixel_num_x)~=0 || length(pixel_num_y)~=0
        surface_pixel_distance(1)=cum_pixel_distance(1);
        for k=2:size(cum_pixel_distance,2)
            surface_pixel_distance(1,k)=cum_pixel_distance(k)-cum_pixel_distance(k-1);
        end
    end
    

%%%%%%%%%%%%%%%%%%%%%%%% summary of variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. cumulative surface distances from the tx point to the crossing points (cum_to)
% 2. surface distances in of each pixel (total_pixel_distance)
% 3. pixel numbers traveled by the ray

    
%%%% Calculate arc length of a ray in each layer
% x_dist =dist([ACO_lat tx_lat],[ACO_lon tx_lon]);
x_dist = tot_sep;
% azmth = azimuth(tx_lat,tx_lon,ACO_lat,ACO_lon);
azmth = az;
[arc_lengths,tot_arclength,theta0,SS_priori,z,SS_HOT_avg,surface_dist,r,est_tt,ray_angles,R_alpha] = ray_trace_w_earth_flattening(x_dist,tx_altitude,tx_lat,azmth,ACO_lat,ACO_lon,ACO_depth);
cl = SS_priori(1);
cr = SS_priori(end);

% Flatten EOF
z_EOF = -R_alpha.*log(R_alpha./(R_alpha-(-z_EOF)));
EOF = EOF.*(R_alpha./(R_alpha-(-z_EOF)));

cum_surface_dist = cumsum(surface_dist);
    
% match total surface distance from ray calculation with that of grid seperation
x_point = length(total_pixel_num)-1;
x_surface_dist = zeros(1,x_point);
z_ind = zeros(1,x_point);

for iii =1:x_point
    [~,z_ind(iii)] = min(abs(cum_surface_dist -cum_pixel_distance(iii)));
    x_surface_dist(iii) = cum_surface_dist(z_ind(iii));
    
end
x_surface_dist(end+1) = cum_surface_dist(end);    
z_ind = [1 z_ind length(z)];
z_cross =  z(z_ind) ;

% calculate arc in each pixel based on the z_cross
arc_section = zeros(1,length(total_pixel_num));
for iii = 1:length(total_pixel_num)
    arc_section(iii) = sum(arc_lengths(z_ind(iii):z_ind(iii+1)-1));
    
end

% adjust the length of EOF
[~,I]= min(abs(z_EOF-z(1)));
if I ~= 1
    z_EOF(1:I-1) = [];
    EOF(1:I-1) = [];
end

int = length(z)-length(z_EOF);
EOF =vertcat(EOF,ones(int,1)*EOF(end));

%%%%%%%%%%%%%% form a row of matrix G %%%%%%%%%%%%%%
G = zeros(1,(length(x_tick)-1)^2);
intd_G = zeros(length(total_pixel_num),1);      % integrand of G in each pixel
r_layer = zeros(length(total_pixel_num),1);     % horizontal distance in each pixel
intd_G_all = zeros(length(arc_lengths),1);
for iii = 1:length(total_pixel_num)    
    pixel_no = total_pixel_num(iii);
    indx = z_ind(iii):z_ind(iii+1)-1;
    G(pixel_no) = sum(-EOF(indx).*arc_lengths(indx)'./(SS_priori(indx).^2));
    
end
obs_ray = G;  
theta_l = theta0;
theta_r = ray_angles(end);

%%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% 
%% Analyze modal sensitivity/ representation
if type == "simulation"
    
    for iii = 1:length(total_pixel_num)    
        indx = z_ind(iii):z_ind(iii+1)-1;
        r_layer(iii) = sum(r(indx));
        intd_G(iii) = sum(-EOF(indx).*arc_lengths(indx)'./(SS_priori(indx).^2));
        
    end
    
    z_cross = z_cross';
    % all layers integrand values
    for iii = 1:length(arc_lengths)
        
        intd_G_all(iii) = sum(-EOF(iii).*arc_lengths(iii)'./(SS_priori(iii).^2));
    end
    
    
    intd_G_remat = repmat(intd_G',4,1);
%     figure(22)
%     plot(intd_G_remat(:),[ones(1,2)*z_cross(1) sort(repmat(z_cross(2:end-1),1,4)) ones(1,2)*z_cross(end)],'Color','b')
%     grid on
%     axis tight
%     set(gca,'YDir','reverse')
%     xlim([-1e-3 1e-3])
%     title('')
%     pause
    %% validate the calculation
    % plot S vs depth and EOF
    if EOF_plot
        mid_z = (z(2:end)+z(1:end-1))/2;
        figure(99)
        set(gcf,'Units','normalized','Position',[0.4 0.8 0.4 0.6])
        clf
        subplot(1,3,1)
        plot(EOF,[z(1); mid_z])
        set(gca,'YDir','reverse')
        grid on
        xlabel('Amplitude')
        ylabel('Depth(m)')
        axis tight
        set(gca,'Yscale','log')
        for k =1:length(z_cross)
            line([-1 1],ones(1,2)*z_cross(k),'Color','k')
        end
        z_cross = z_cross';
        xlim([min(EOF) max(EOF)])
        ylim([z(1) z(end)])
        title('EOF')
        
        subplot(1,3,2)
        new_arc_section = repmat(arc_section,4,1);
        plot(new_arc_section(:),[ones(1,2)*z_cross(1) sort(repmat(z_cross(2:end-1),1,4)) ones(1,2)*z_cross(end)],'Color','r')
        set(gca,'YDir','reverse')
        ylim([z(1) z(end)])
        ylabel('Depth (m)')
        xlabel('meter')
        grid on
        xlim([min(new_arc_section(:))-200 max(new_arc_section(:))+200])
        title('Arc Length (Pixel-Wise')
        
        subplot(1,3,3)
        plot(ray_angles*180/pi,z(2:end-1))
        set(gca,'YDir','reverse')
        ylim([z(1) z(end)])
        grid on
        xlabel('Degrees')
        ylabel('Depth (m)')
        title('Incident Ray Angle')
        pause
        
    end
end

end