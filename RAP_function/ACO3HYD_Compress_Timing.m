function [y, fs, n] = ACO3HYD_Compress_Timing(filename1,filename2,save_dir)
% read packed ALOHA hydrophone data
% Modified to comply with a new header format implemented on December 2016
% Outputs:
%   y = signal
%   fs = sample rate
%   n = date number

 fs = 24000; % decimated sample rate
%fs = 96000; % original sample rate

%channela=zeros(fs*5*60,1);
nn=0;
fid = fopen(filename1,'r');
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

    % bytes to read = bits*4096/8;
    % samples must be parsed into 4096 values
    if feof(fid)<1
        bitn=['bit' num2str(bits)]; % this is the call in fread for sample bit length
        cha=fread(fid,4096,bitn,'ieee-le');
        if size(cha) <4096
            cha(1:4096)=0;
        end
        if nn==0
            channela(1:4096)=cha;
        elseif nn>=1
            cstart=4096*nn+1;
            cstop=(nn+1)*4096;
            channela(cstart:cstop)=cha;
        end
    end
    nn=nn+1;
end
y = channela;
fclose(fid);

% ******************* Filename data ***********************
% filenames are all in format:
% 'YYYY-MM-DD--hh.mm.HYD'
YY = str2num(filename1(1:4));
MM = str2num(filename1(6:7));
DD = str2num(filename1(9:10));
hh = str2num(filename1(13:14));
mm = str2num(filename1(16:17));
ss = 00;
% Convert date/time to a number
n = datenum(YY,MM,DD,hh,mm,ss);





%%%%%%%Save Time Information%%%%%%%%%
AtSec=cell2mat(AtSec(~cellfun('isempty',AtSec)));
AtnSec=cell2mat(AtnSec(~cellfun('isempty',AtnSec)));
CorrectionOffset_HnSec=cell2mat(CorrectionOffset_HnSec(~cellfun('isempty',CorrectionOffset_HnSec)));
timecount=cell2mat(timecount(~cellfun('isempty',timecount)));

Date_num=datenum(1970,1,1,0,0,cell2mat(tSec)+(cell2mat(tuSec)./1000000));
Date_num_count=datenum(1970,1,1,0,0,cell2mat(tSec)+(timecount./10000000));

change_sec=Date_num-Date_num_count;
Date_num_count(change_sec<0)=Date_num_count(change_sec<0)-(1/(3600*24));

mat_filename=sprintf('%s.mat',filename2);
cd(save_dir) 
save(mat_filename,'AtSec','AtnSec','CorrectionOffset_HnSec','Date_num','Date_num_count')



