function [y,t_date] = icListen_audio_prep(fname)
% Prepare/ assign timing variables for the icListen audio file given its name

    %Load Audio Data
    [y,fs] = audioread(fname,'double');                          % wav acoustic file
    y = y*2^23;
    % Timing Variables
    year = fname(9:12);
    mm = fname(13:14);
    dd = fname(15:16);
    hr = fname(18:19);
    min = fname(20:21);
    sec = fname(22:23);
    formatIn = 'yyyy-mm-dd HH-MM-SS';
    date_string = string(year)+string("-"+mm)+string("-"+dd)+string(" "+hr)+string("-"+min)+string("-"+sec);
    start_time = datenum(date_string,formatIn);
    % time increment
    T = 1/fs;
    time_inc = T /(3600*24); % Fraction of day
    timestamp_num = length(y);     % number of timestamps in the files
    t_date = zeros(1,timestamp_num);    % pre-allocate space for timestamps
    for jj = 1:timestamp_num
        t_date(jj) = start_time+(jj-1)*time_inc;
    end
    
end