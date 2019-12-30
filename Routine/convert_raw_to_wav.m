% Unpack a directory of packed HEM data files
% adapted from Ethan Roth's code
% Bruce Howe 
% 20141111
% hardwired to only reading from one directory specified, and writing wav
% files to local directory with output name format changed (I think to
% conform to Hildbrand's Triton format)
% EDIT the hydrophone file directory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format compact
format long
close all
clearvars

% get directory and info on HYD files
ii = 1;
%HYD.ddir = 'Z:\Data_ACO\HYDR024K\';     % default directory%
%%%%%%%%% EDIT %%%%%%%%%%%%%
% raw file directoty
HYD.ddir = '/Volumes/ACO_RAP_2/RAP/Oct2018Cruise/wav_data/HEM/27_Oct';
% save directory 
save_dir = '/Volumes/ACO_RAP_2/RAP/Oct2018Cruise/wav_data/HEM/27_Oct/2';


%HYD.idir{ii} = uigetdir(HYD.ddir,'Select directory with only HYDRO pk files');
% if the cancel button is pushed, then no file is loaded so exit this script
%if strcmp(num2str(HYD.idir{ii}),'0')
%    disp('canceled - no directory chosen')
%    return
%else
%    disp('HYDRO pk file directory : ')
%    disp([HYD.idir{ii}])
%end
%HYD.idir{ii} = HYD.ddir

% **** fix line
d = dir(fullfile(HYD.ddir,'/','*.HYD24BBpk'));    % directory info
%d = dir(fullfile(HYD.ddir,'/','*.HYD24HFpk'));    % directory info
HYD.fname{ii} = char(d.name);             % file names
fnsz = size(HYD.fname{ii})                 % file size
HYD.nfiles{ii} = fnsz(1)   % number of files in directory
disp(['number of files in directory is ',num2str(HYD.nfiles{ii})])
HYD.inpath = [HYD.ddir,'/'];
HYD.outpath = HYD.inpath
filenames=HYD.fname{ii}
numfiles=HYD.nfiles{ii}


% Convert files
newfile = 1;
ii = 1;
for jj = 1:HYD.nfiles{ii}
    
    HYD.infile = deblank(HYD.fname{ii}(jj,:));
    tstr = [];
    tstr = regexp(HYD.infile,'\d\d\d\d-\d\d-\d\d--\d\d.\d\d','match');
    HYD.dnum = datenum(tstr,'yyyy-mm-dd--HH.MM');
    HYD.dstr = datestr(HYD.dnum,'yymmdd-HHMMSS')
    filename2 = [HYD.dstr];
    filename1 = [HYD.inpath,HYD.infile];
    %[y,fs,n] = ACO3HYD_Compress(filename1);
    [y,fs,n] = ACO3HYD_Compress_Timing(filename1,filename2,save_dir);
    filename2 = sprintf('%s.wav',filename2)
    
    % Convert y to uPa and write .wav file
    % y = y/max(abs(y));
    % y = y.*(20.4);
    n = 32;
    %!pwd
    
    cd(save_dir)
    audiowrite(filename2,y,fs,'BitsPerSample',n)
    cd(HYD.ddir)
end




