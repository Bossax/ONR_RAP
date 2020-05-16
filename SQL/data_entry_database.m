% Insert data into database
% read .mat files from data folder
% 1. EDIT table name
% 2. EDIT data to insert to the table
%%%%%%%%%%%%%%%%%%%%%%
% field names in database
% 1. tx_time Primary Key
% 2. longitide
% 3. latitude
% 4. altitude
% 5. heading
% 6. x_velocities
% 7. surface_range
% 8. x_err
% 9. y_err
% 10. z_err
% 11. travel_time
% 12. SNR
%%%%%%%%%%%%%%%%%%%%%%
clear
close all
%%%%%%%%%%%%%%%%%%%%
% Octoebr 2018
day = 27:30;
start_hour = 3;
end_hour = 14;

% [tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,x_err,y_err,z_err,act_arrival,est_arrival,SNR] = tx_rx_extraction_Oct(day,start_hour,end_hour,'HEM');
[tx_t,tx_lon,tx_lat,tx_heading,tx_altitude,tx_xvel,range,x_err,y_err,z_err,act_arrival,est_arrival,SNR] = tx_rx_extraction_Oct(day,start_hour,end_hour,'icListen');

% create tabulated dataset
ttp = (act_arrival - est_arrival)*3600*24*2400; % millisecond
% Define field names
% refer to database table
date_format = 'yyyy-mm-dd HH:MM:SS'; % SQL timestamp YYYY-MM-DD HH:MI:SS
field_name = {'tx_time';'longitude';'latitude';'altitude';...
    'heading';'x_velocity';'surface_range';'x_err';'y_err';'z_err';...
    'traveltime_perturbation';'SNR'};

% data_table = table(string(datestr(tx_t,date_format)),tx_lon',tx_lat',tx_altitude',tx_heading',tx_xvel',range',x_err',y_err',z_err',ttp',SNR',...
%             'VariableNames',field_name);

data_table = table(string(datestr(tx_t,date_format)),tx_lon',tx_lat',tx_altitude',tx_heading',tx_xvel',range',x_err',y_err',z_err',ttp',SNR', ...
    'VariableNames',field_name);

disp(sprintf('Connecting to database ... \n'))
% create a connection to the database
datasource = 'RAP';     % schema
conn = database(datasource,'root','Kea2c7ra$',...
    'Vendor','MySQL', ...
    'Server','localhost', ...
    'PortNumber',3306);

% check connection status
if isempty(conn.Message)
    fprintf('Connection is successful \n')
end

% check existing table
% tablename = 'HEM_original_position';       % EDIT
tablename = 'icListen_original_position';       % EDIT

ex_data = sqlread(conn,tablename);
[row_no,~] = size(ex_data);
fprintf('Existing Rows = %i \nPreview..\n',row_no)
tail(ex_data,4)

% ignore duplicate rows
% append data to the table
while true
    try
        sqlwrite(conn,tablename,data_table)
        updated_data = sqlread(conn,tablename);
        [new_row_no,~] = size(ex_data);
        fprintf('Updated Table = %i \nPreview..\n',new_row_no)
        tail(updated_data,4)
        break
    catch ME
        err_msg = ME.message;
        if contains(err_msg,'Duplicate entry')
            % check for duplicate data
            [dup_row,rm_ind,~] = intersect(datenum(data_table.tx_time),datenum(ex_data.tx_time),'stable');
            data_table(rm_ind,:) = [];
            if isempty(data_table)
                disp('Redundant data: No new data has been added to the database')
                break
            end
        else
            rethrow(ME)
        end
    end
end
close(conn)



