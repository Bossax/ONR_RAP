function [total_pixel_distance,total_pixel_num]=lat_lon_dist(tx_lat,tx_lon,ACO_lat,ACO_lon,x_tick,y_tick)
%  oberservation_matrix construction
%   Find the surface distance an acoustic ray travels over X and Y
%  y_tick = y coordinates of pixel boundaries from higher to lower

domain_size = length(x_tick)-1;
stop_y=ACO_lat;
stop_x=ACO_lon;
start_y=tx_lat;
start_x=tx_lon;

% creates 10,000 points in the line between tx and ACO
res_dist = 20;
sep_point = 10000;
while res_dist > 0.5
    [range,lat_cur,lon_cur]=dist([start_y,stop_y],[start_x,stop_x],sep_point);
    res_dist = (range(end)-range(1))/(sep_point-1);
    sep_point = sep_point*2;
end

%%%%%%X-INT%%%%%%%
% tx point is on the west
if start_x<stop_x
    % find the range of X
    x_pos1=min(x_tick(x_tick>=start_x));    % lower limit
    x_pos2=max(x_tick(x_tick<=stop_x));     %  upper limit
    
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
    x_pos1=max(x_tick(x_tick<=start_x));
    x_pos2=min(x_tick(x_tick>=stop_x));
    if size(x_pos2,2)~=0 && size(x_pos1,2)~=0
        x_int=fliplr(x_tick(x_tick<=x_pos1 & x_tick>=x_pos2));
    elseif size(x_pos1,2)~=0
        x_int=x_pos1;
    elseif size(x_pos2,2)~=0
        x_int=x_pos2;
        
    elseif start_x==stop_x
        x_int=[];
    end
end

% calculate distances of points in the interval
if size(x_int,2)>=1
    for k=1:size(x_int,2)
        % find the nearest lon_cur to the x_interval point from the tx location side
        if start_x<stop_x
            dist_pos = find(lon_cur <= x_int(k),1,'last');
            % tx point is on the east
        elseif start_x>stop_x
            dist_pos = find(lon_cur >= x_int(k),1,'first');
        end
        % Assign lat/lon to pixel boundary handles
        lon_pos(k)=lon_cur(dist_pos);
        lat_pos(k)=lat_cur(dist_pos);
        pixel_distance_x(k)=range(dist_pos);
        
        % Loop over all pixels to assign pixel number to the segment
        % find nearest x,y nodes to the lat/lon above
        % counting direction: top to bot/left to right
        % count the pixel in which the ray crosses the boundary
        for ii=2:length(y_tick)
            for jj=2:length(x_tick)
                % check if the lon/lat is within that pixel
                if (lon_pos(k)<x_tick(jj)) && (lon_pos(k)>x_tick(jj-1)) && (lat_pos(k)<y_tick(ii)) && (lat_pos(k)>y_tick(ii+1))
                    % pixel number
                    
                    pixel=(jj-1)+(ii-1)*domain_size;    %  column + roll
                    
                end
                
            end
        end
        pixel_num_x(k)=pixel(end);
        
    end
end


%%%%%%%Y-INT%%%%%%%%%%%%5
% tx lat is in the south
if start_y<stop_y
    % find the range of Y
    y_pos1=min(y_tick(y_tick>=start_y));
    y_pos2=max(y_tick(y_tick<=stop_y));
    % find points in between
    if size(y_pos2,2)~=0 && size(y_pos1,2)~=0
        y_int=y_tick(y_tick>=y_pos1 & y_tick<=y_pos2);
    elseif size(y_pos1,2)~=0
        y_int=y_pos1;
    elseif size(x_pos2,2)~=0
        y_int=y_pos2;
    end
    
    % tx lat is in the north
elseif start_y>stop_y
    % find the range of Y
    y_pos1=max(y_tick(y_tick<=start_y));
    y_pos2=min(y_tick(y_tick>=stop_y));
    % find points in between
    if size(y_pos2,2)~=0 && size(y_pos1,2)~=0
        y_int=fliplr(y_tick(y_tick<=y_pos1 & y_tick>=y_pos2));
    elseif size(y_pos1,2)~=0
        y_int=y_pos1;
    elseif size(x_pos2,2)~=0
        y_int=y_pos2;
    end
elseif start_y==stop_y
    y_int=[];
end

y_int_handle = fliplr(y_int);
lat_pos = [];
lon_pos = [];
if size(y_int,2)>=1
    for k=1:length(y_int_handle)
        % find the nearest lat_cur to the y_interval point from the tx location side
        % in the south
        if start_y<stop_y
            dist_pos = find(lat_cur <= y_int_handle(k),1,'last');
            % tx point is on the north
        elseif start_y>stop_y
            dist_pos = find(lat_cur >= y_int_handle(k),1,'first');
        end
        lon_pos(k)=lon_cur(dist_pos);
        lat_pos(k)=lat_cur(dist_pos);
        pixel_distance_y(k)=range(dist_pos);
        
        for ii=2:size(y_tick,2)
            for jj=2:size(x_tick,2)
                
                if lon_pos(k)<x_tick(jj) && lon_pos(k)>x_tick(jj-1) && lat_pos(k)<y_tick(ii) && lat_pos(k)>y_tick(ii+1)
                    
                    pixel=(jj-1)+(ii-1)*(domain_size);
                    
                end
                
            end
        end
        pixel_num_y(k)=pixel(end);
        
    end
end



%Final Pixel
pixel_num_last = [];
pixel_distance_final = [];
if isempty(intersect(x_tick,ACO_lon)) && isempty(intersect(y_tick,ACO_lat))
pixel_distance_final=range(end);
for ii=2:length(y_tick)
    for jj=2:size(x_tick,2)
        
        if (stop_x<x_tick(jj)) && (stop_x>x_tick(jj-1)) && (stop_y<y_tick(ii)) && (stop_y>y_tick(ii+1))
            pixel=(jj-1)+(ii-1)*domain_size;
            
        end
        
    end
end
pixel_num_last=pixel(end);
end

if exist('pixel_num_x')~=0 && exist('pixel_num_y')~=0
    total_pixel_distance_temp=horzcat(pixel_distance_x,pixel_distance_y,pixel_distance_final);
    [total_pixel_distance_temp,b]=sort(total_pixel_distance_temp);
    total_pixel_num=horzcat(pixel_num_x,pixel_num_y,pixel_num_last);
    total_pixel_num=total_pixel_num(b);
elseif exist('pixel_num_x')~=0
    total_pixel_distance_temp=horzcat(pixel_distance_x,pixel_distance_final);
    [total_pixel_distance_temp,b]=sort(total_pixel_distance_temp);
    total_pixel_num=horzcat(pixel_num_x,pixel_num_last);
    total_pixel_num=total_pixel_num(b);
elseif exist('pixel_num_y')~=0
    total_pixel_distance_temp=horzcat(pixel_distance_y,pixel_distance_final);
    [total_pixel_distance_temp,b]=sort(total_pixel_distance_temp);
    total_pixel_num=horzcat(pixel_num_y,pixel_num_last);
    total_pixel_num=total_pixel_num(b);
else
    total_pixel_distance=pixel_distance_final;
    total_pixel_num=pixel_num_last;
end

if exist('pixel_num_x')~=0 || exist('pixel_num_y')~=0
    total_pixel_distance(1)=total_pixel_distance_temp(1);
    for k=2:size(total_pixel_distance_temp,2)
        total_pixel_distance(1,k)=total_pixel_distance_temp(k)-total_pixel_distance_temp(k-1);
    end
end




end