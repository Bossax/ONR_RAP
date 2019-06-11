function time_lock()
today = datenum('28-10-2018','dd-mm-yyyy');
% Extract the present timestamp when the program starts running
    p_time = datenum(datetime('now'));
    diff = p_time-today;
    hr = floor(diff*24)
    min = floor((diff*24-hr)*60)
    sec = floor(((diff*24-hr)*60-min)*60)
    mmm = ((((diff*24-hr)*60-min)*60)-sec)*1000
    
 % Degsinate the next integer minute
 if sec <= 53
    sec_next = 54;
    mmm_next = 800;
    min_next = min;
    hr_next = hr;
 elseif sec > 53
     % wait for the next minute
     sec_next = 54;
     min_next = min+1;
     hr_next = hr;
     mmm_next = 800;
     if min_next >=60
         hr_next = hr+1;
         min_next = 0;
     end
 end
 
% loop until the time is within one number
threshold = 1; % 1 ms
err = 100;  % initial error 100 ms
while err>threshold
    % read the real-time timestamp
    p_time = datenum(datetime('now'));
    diff = p_time-today;
    hr = floor(diff*24);
    min = floor((diff*24-hr)*60);
    sec = floor(((diff*24-hr)*60-min)*60);
    mmm = ((((diff*24-hr)*60-min)*60)-sec)*1000;
    if sec == sec_next && min == min_next && hr == hr_next && abs(mmm-mmm_next) <= threshold
        err = abs(mmm-mmm_next)
        mmm
        mmm_next
    end
    
end

time = sprintf("Run time = %.0f:%.0f:%.0f%.6f \n",hr,min,sec,mmm/1000);
fprintf(time)
end