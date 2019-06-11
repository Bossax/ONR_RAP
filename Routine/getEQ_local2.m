global per_cent_noise

% ACO location assumed for now
    aco_loc = [22.7387 -158.0062];
    cdMATLAB;
    cd('/Users/rhett/Desktop/Seismo/Hawaii_local')
    
%  zeroth, make sure the HEM directory is mounted 
%   (though not necessary for archived events)
    if exist('/Volumes/HEM_data' , 'dir') ~= 7
        fprintf('*** /Volumes/HEM_data is not mounted\n')
    end

% first, find the .xlsx file name
    disp('Please select a .xlsx event file...')
	[filename,pathname] = uigetfile('*.xlsx');
	source = fullfile(pathname,filename);
    if ~ischar(source) || ~exist(source,'file')
        error('File %s does not exist.',source);
    else
        fprintf('...accessing file: %s\n', source);
    end
    
% read the file...
    [k_events, k_text] = xlsread(source);
    n_events = length(k_events(:,1));
    if n_events == 1
        fprintf('...N = %d event\n', n_events)
    else
        fprintf('...N = %d events\n', n_events)
    end
    
% assign variables 
    latlon   =[k_events(:,5) k_events(:,6)];
    depth    = k_events(:,7);
    mag      = k_events(:,8);
    err_NE   = k_events(:,19);
    err_Z    = k_events(:,20);
    location = k_text(:,16);
    matlab_time = datenum(datetime(k_text(:,1),'Format', 'yyyy/MM/dd'))...
        + datenum(datetime(k_text(:,2),'Format', 'HH:mm:ss.SS')...
        - floor(datenum(datetime(k_text(:,2),'Format', 'HH:mm:ss.SS'))));
    [delta, az] = distance([latlon(:,1)	latlon(:,2)], aco_loc);
    km = delta * 111.19;
    [Year,Month,Day,Hour,Minute,Second] = datevec(matlab_time(:));
    EventTag = 'Kilauea_';
    
 %     event before sensor start date ?
%     icL_commence = datenum(2018,6,22);
    hydro_commence = datenum(2011,6,9);
%     n_events = 1;
    for is = 1:n_events
        if matlab_time(is) < hydro_commence
            fprintf('*** event %s is before hydro installation\n',datestr(matlab_time(is))) 
%         elseif matlab_time(is) < icL_commence
%             fprintf('*** event %s is before IcListen installation\n',datestr(matlab_time(is))) 
        end
    end
    
%   loop over events
%      n_events = 1;      % one for now
     N_5_minutes = 3;   % three 5-minute files
    for iq = 1:n_events
      % de-archive if available
        EQname = sprintf('%s%04d%02d%02d_%02d%02d%02d,%1.0f.mat', EventTag, Year(iq), Month(iq), Day(iq), Hour(iq), Minute(iq), floor(Second(iq)), 10*(Second(iq) - floor(Second(iq))  ))
        load_flag = 0;
        lastCD = cd('/Users/rhett/Documents/MATLAB/ACO_archive');
        if exist(EQname, 'file') == 2
            load(EQname);
            load_flag = 1;
            origin_time_arch = matlab_time(iq);
            Toffset_1 = (matlab_time(iq) - EQ.Time(1))*86400; % ok
            plota(EQ.Time, EQ.Data, origin_time_arch, Toffset_1)
            continue
        end
        
        % not archived, so download...
        [YearNow, MonthNow, DayNow] = datevec(now);
            % determine Hour and Minute of file start
            if YearNow == Year(iq) && MonthNow == Month(iq)
                % current month, not yet in main archive
                [File_0, Toffset, ACOdirectory] = read5minutes('interim', Year(iq), Month(iq), Day(iq), Hour(iq) , Minute(iq), Second(iq));
            else
                % archived          
               [File_0, Toffset, ACOdirectory] = read5minutes('archive', Year(iq), Month(iq), Day(iq), Hour(iq) , Minute(iq), Second(iq));
            end             
        [channel_24K, time_24K] = get_hydro(N_5_minutes, File_0, Toffset, ACOdirectory,icLdirectory);
        Toffset = (matlab_time(iq) - File_0)*86400
        h_400 = ACO_decimate24k(channel_24K);
        t_400 = downsample(time_24K, 3*4*5)';
        sz_min = min(length(h_400),length(t_400));
        h_400 = h_400(1:sz_min);
        t_400 = t_400(1:sz_min);
        if load_flag == 0
            archive_event(EQname, h_400, t_400);
        end
    
        origin_time = matlab_time(iq);
        plota(t_400, h_400, origin_time, Toffset)
        EQfigure = sprintf('%s%04d%02d%02d_%02d%02d%02d,%1.0f.fig', EventTag, Year(iq), Month(iq), Day(iq), Hour(iq), Minute(iq), floor(Second(iq)), 10*(Second(iq) - floor(Second(iq))  ))
        saveas(74, EQfigure);
        
    end  % n_events loop
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save the event
function archive_event(EQname, h_400, t_400);
    EQ = timeseries(h_400, t_400);
    lastCD = cd('/Users/rhett/Documents/MATLAB/ACO_archive');
    save(EQname, 'EQ');
    cd(lastCD);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plota(t_400, h_400, origin_time, Toffset)
global per_cent_noise
    
    monikker = datestr(origin_time, 'yyyy.mm.dd');
%   various plots du jour
    figure(71); clf;
        plot(t_400, h_400, 'Color', [0 0.5 0]);
        datetick('x','HH:MM');
        title(monikker);


    figure(72); clf;
        ts_400 = ( t_400 - origin_time )*846400;
        plot(ts_400, h_400, 'Color', [0 0.5 0.5]);
        xlabel('Time, sec');
        title(monikker);

    figure(73); clf;
        per_cent_noise = 1/3;
        nfft = 256;  % 1.6 sec
        rb_spectrogram_H5(h_400, -Toffset, nfft, ceil(nfft/2), nfft, 400, 'yaxis');
        title(monikker);
        colormap jet

                
     figure(74); clf;
        per_cent_noise = 1/3;
        nfft = 256;  % 1.6 sec
        rb_spectrogram_H5_hydroresp(h_400, -Toffset, nfft, ceil(nfft/2), nfft, 400, 'yaxis');
        title(monikker);          
        colormap jet
        
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [File_0, Toffset, ACOdir] = read5minutes(data_location, Year, Month, Day, Hour, Minute, Second)
        % Toffset is the time sec from file start to origin time
        origin_time = datenum(Year, Month, Day, Hour, Minute, Second);
        File_0 = datenum(Year, Month, Day, Hour, Minute - mod(Minute,5) -5 , 0);
        Toffset = (origin_time - File_0)*86400;
        [Year_0, Month_0, Day_0, ~, ~, ~] = datevec(File_0);
        if contains(data_location,'archive')
           ACOdir = sprintf('/Volumes/HEM_data/ARCHIVE_ACO/HYDR024K/%04d/%02d/%02d',[Year_0, Month_0, Day_0]');
%            icLdir = sprintf('/Volumes/iclisten_data/icListen_wav/%04d/%02d/%02d/SBW1391_%04d%02d%02d_%02d%02d00.wav', Year_0, Month_0, Day_0, Year_0, Month_0, Day_0, Hour_0, Minute_0);  
        else
           ACOdir = sprintf('/Volumes/HEM_data/DATA_ACO/HYDR024K/%02d/%02d',[Month_0, Day_0]');
%            icLdir = sprintf('/Volumes/iclisten_data/icListen_wav/%04d/%02d/%02d/SBW1391_%04d%02d%02d_%02d%02d00.wav', Year_0, Month_0, Day_0, Year_0, Month_0, Day_0, Hour_0, Minute_0);  
        end 
        if exist(ACOdir, 'dir') ~= 7
            warning('*** HEM directory not found: %s\n', ACOdir);
        end
%         if exist(icLdir, 'file') ~= 2
%             warning('*** icL file not found: %s\n', icLdir);
%         end

    end
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    
    function [channel_24K ,time_24K] = get_hydro(N_5_minutes, File_0, Toffset, ACOdirectory,icLdirectory);
            % confirm folder exists
            if exist(ACOdirectory, 'dir') ~= 7
                    warning('*** HEM directory not found: %s\n', ACOdirectory);
            end
            
%           ACO HEM data 24K data
            channel_24K = [];
            time_24K = [];
%           channel_24K =   zeros(N_5_minutes*5*60*24000,1); % zero out data vector 7200000
%           time_24K    =   zeros(N_5_minutes*5*60*24000,1);
            dt = 1/24000;
            
            record_number = 0;  % record number for nfile = 1:N_5_minutes
%           [length = number of files*5 minutes/file * 60 sec/min * 24000 samp/sec]           
            [Year_0, Month_0, Day_0, Hour_0, Minute_0, Second_0] = datevec(File_0);
            
            
        for nfile = 1:N_5_minutes
            % shift by 5 minutes per file
            File_1 = datenum(Year_0, Month_0, Day_0, Hour_0, Minute_0 + 5 * (nfile-1), Second_0);
            [Year_1, Month_1, Day_1, Hour_1, Minute_1, ~] = datevec(File_1);
            File_n = sprintf('%s/%04d-%02d-%02d--%02d.%02d.HYD24BBpk',ACOdirectory,Year_1, Month_1, Day_1, Hour_1, Minute_1);
            if exist(File_n, 'file') ~= 2
                warning('*** HEM file not found: %s\n', File_n);
            end      
            fprintf('%s\n',File_n)
            fid = fopen(File_n, 'r');           

%{
             record_number = 0;  % record number
%             while feof(fid) < 1  % straight from Fred Duennebier !
%                 Record          = fread(fid,1,'uint32','ieee-le');%1
%                 Decimation      = fread(fid,1,'uint8','ieee-le');%2
%                 StartofFile     = fread(fid,1,'uint8','ieee-le');%3
%                 Sync1           = fread(fid,1,'uint8','ieee-le');%4
%                 Sync2           = fread(fid,1,'uint8','ieee-le');%5
%                 Statusbyte1     = fread(fid,1,'uint8','ieee-le');%6
%                 Statusbyte2     = fread(fid,1,'uint8','ieee-le');%7
%                 pad1            = fread(fid,1,'uint8','ieee-le');%8
%                 LeftRightFlag   = fread(fid,1,'uint8','ieee-le');%9
%                 tSec            = fread(fid,1,'uint32','ieee-le');%10
%                 tuSec           = fread(fid,1,'uint32','ieee-le');%11
%                 timecount       = fread(fid,1,'uint32','ieee-le');%12
%                 Year            = fread(fid,1,'int16','ieee-le');%13
%                 yDay            = fread(fid,1,'int16','ieee-le');%14
%                 Hour            = fread(fid,1,'char','ieee-le');%15
%                 Min             = fread(fid,1,'char','ieee-le');%16
%                 Sec             = fread(fid,1,'char','ieee-le');%17
%                 Allignment      = fread(fid,1,'char','ieee-le');%18
%                 sSec            = fread(fid,1,'int16','ieee-le');%19
%                 pad2            = fread(fid,1,'uint8','ieee-le');%20
%                 bits            = fread(fid,1,'uint8','ieee-le');%21 in packed format, size of sample in bits 
            %}
            % data format is ,'ieee-le' (that's the letter L)
            tSec={};
            tuSec={};
            timecount={};
            AtSec = {};
            AtnSec= {};
            CorrectionOffset_HnSec= {};
            ID = {};
            LeftRightFlag = {};
            bits=[];
            y = [];         % data row vector
                    
                while feof(fid) < 1
                    Record                  = fread(fid,1,'uint32','ieee-le'); % bytes 1-4
                    Record_duplicate        = fread(fid,1,'uint32','ieee-le'); % bytes 5-8
                    Sync1                   = fread(fid,1,'uint8','ieee-le'); % 9
                    Sync2                   = fread(fid,1,'uint8','ieee-le'); % 10
                    Statusbyte1             = fread(fid,1,'uint8','ieee-le'); % 11
                    Statusbyte2             = fread(fid,1,'uint8','ieee-le'); % 12
                    tSec{end+1}             = fread(fid,1,'uint32','ieee-le'); % 13-16
                    tuSec{end+1}            = fread(fid,1,'uint32','ieee-le'); % 17-20
                    timecount{end+1}        = fread(fid,1,'uint32','ieee-le'); % 21-24
                    AtSec{end+1}            = fread(fid,1,'uint32','ieee-le'); % 25-28
                    AtnSec{end+1}           = fread(fid,1,'uint32','ieee-le'); % 29-32
                    CorrectionOffset_HnSec{end+1}  = fread(fid,1,'uint8','ieee-le'); % 33
                    ID{end+1}               = fread(fid,1,'uint8','ieee-le'); % 34
                    LeftRightFlag{end+1}    = fread(fid,1,'uint8','ieee-le'); % 35
                    bits                    = fread(fid,1,'uint8','ieee-le');  % 36 - SampleSize
                    
    
                    record_number=record_number+1;      % add 1 to records read

                    if feof(fid) < 1
                        bitn=['bit' num2str(bits)];      % this is the call in fread for sample bit length
                        cha=fread(fid, 4096, bitn,'ieee-le');
                        sizecha=size(cha);
                        if sizecha < 4096
                            cha(sizecha+1:4096) = 0;    % complete 4096 values  not read in
                        end
                        if record_number == 1           % if this is the first record
                            y(1:4096) = cha;
                        else
                            cstart = 4096*(record_number-1)+1;
                            cstop = (record_number)*4096;
                            y(cstart:cstop) = cha;       % add new data to data vector
                        end
                    end
                end             % while loop
                
                fclose(fid);    % end of 1 hydrophone file
                
                %%% Saving Timing Information of the hydrophone file %%%%
                %{
                AtSec = cell2mat(AtSec(~cellfun('isempty',AtSec)))
                AtnSec = cell2mat(AtnSec(~cellfun('isempty',AtnSec)))
                CorrectionOffset_HnSec=cell2mat(CorrectionOffset_HnSec(~cellfun('isempty',CorrectionOffset_HnSec)))
                %}
                timecount = cell2mat(timecount(~cellfun('isempty',timecount)))
                
                Date_num=datenum(1970,1,1,0,0,cell2mat(tSec)+(cell2mat(tuSec)./1000000));
                datestr(Date_num,'yyyy/mm/dd HH:MM:SS.FFF')
                
                Date_num_count=datenum(1970,1,1,0,0,cell2mat(tSec)+(timecount./10000000));
                datestr(Date_num_count,'yyyy/mm/dd HH:MM:SS.FFF')
                
                change_sec=Date_num-Date_num_count;
                Date_num_count(change_sec<0)=Date_num_count(change_sec<0)-(1/(3600*24));
                
                %%% align timestamps
                time_exact = Date_num_count-(.00005/(3600*24));           % 50 us lag between 96 and 24kHz data
                time_interval = diff(time_exact).*(3600*24);              % Find time intervals between headers
                time_interval(time_interval<0) = [];
                mean_time_interval = mean(time_interval);                  % Mean time interval
                single_time_interval = mean_time_interval/4096;            % Single time interval
                mean_date_interval = mean_time_interval/(3600*24);
                single_date_interval = single_time_interval/(3600*24);
                time_diff=(time_exact(end)-time_exact(1))/(1/(3600*24))+mean_time_interval;
                
                t_date=[];
                t_date(1:512)=linspace(time_exact(1)-(mean_date_interval/8),time_exact(1),512);
                
                for jj=2:length(time_exact)
                    t_date(((jj-2)*4096+513):((jj-1)*4096+512))=linspace((time_exact(jj-1)+single_date_interval),time_exact(jj),4096);
                end
                
                t_date(end+1:end+3584)=linspace((time_exact(end)+single_date_interval),(time_exact(end)+(7*mean_date_interval/8)),3584);
                y = y(1:length(t_date));                     %Fix y (if needed/should no longer be necessary)
                
                % time offset correction
                time_offset = time_correction(t_date(1));
                
                t_date = t_date - time_offset;
                
                %%% concatenate data %%%
                time_24K = horzcat(time_24K,t_date);
                channel_24K = horzcat(channel_24K,y);
                
        end  % for nfile loop
        
        %         Time_0   = datenum(Year_0, Month_0, Day_0, Hour_0, Minute_0, 0);
        %         time_24K = linspace(Time_0, Time_0 + N_5_minutes*5*60/86400, N_5_minutes*5*60/dt);
    end % function
    
    
    
function time_offset = time_correction(time)
% June 2018
date_mark1= "20180620 23:25";
date_mark1 = datenum(date_mark1,'yyyymmdd HH:MM');
     
 % Unknow event
 date_mark2= "20180831 00:00";
 date_mark2 = datenum(date_mark2,'yyyymmdd HH:MM');
 
% October
 date_mark3= "20181028 01:00";
 date_mark4= "20181029 01:00";
 date_mark5= "20181030 01:00";
 date_mark3 = datenum(date_mark3,'yyyymmdd HH:MM');
 date_mark4 = datenum(date_mark4,'yyyymmdd HH:MM');
 date_mark5 = datenum(date_mark5,'yyyymmdd HH:MM');

% November 
 date_mark6= "20181110 01:00";
 date_mark6 = datenum(date_mark6,'yyyymmdd HH:MM');
 
 if time <= date_mark1
     time_offset =  -1/(3600*24);
     
 elseif(date_mark1 <= time)&(time <= date_mark2)
     time_offset = 2/(3600*24);
     
 elseif (date_mark2 <= time)&(time <= date_mark3)
     time_offset =  4/(3600*24);
     
 elseif (date_mark3 <= time)&(time <= date_mark4)
     time_offset = 5/(3600*24);
     
 elseif (date_mark4 <= time)&(time<= date_mark5)
     time_offset = 6/(3600*24);
     
 elseif (date_mark5 <= time)&(time<= date_mark6)
     time_offset = 7/(3600*24);

 elseif (date_mark6 <= time)
     time_offset = 0;
     
 end
 
 end
    