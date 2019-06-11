% IRIG-B decoder
close all
clear

cd '/Volumes/ACO_RAP_2/RAP/Oct2018Cruise/Scarlett_Data/'
fname = 'LFM_1398082157.mat';
load(fname);
fs = 44100;

%% Find TX times
y_tx=[];
t_tx=[];
v_time=[];
pps_ind=[];
num_offset=[];
tx_time = [];

irigb = txrx.sam(1,:);                       
pps=txrx.sam(2,:);                       % PPS     
V = txrx.sam(4,:);                       % Voltage
t =txrx.taxis;                            % timestamps from start
start_t=txrx.start_time;                 % start time

t_tx = horzcat(t_tx,start_t+(t./(3600*24)));    % datenum time

% t_tx = t_tx - 1/(3600*24);      % 1 sec ahead of UTC
peak_num=max(V);
med = median(abs(V));
if peak_num <= 4*med
    threshold = 4*med;
else
    threshold = peak_num/2;
end
        
v_spike=find(abs(V)>=threshold);          %Find voltage spike
 
%Find start of voltage spikes
for i=1:length(v_spike)
    if i==1
        v_time(end+1)=v_spike(i);
    else
    if (v_spike(i)-v_spike(i-1))>fs
        v_time(end+1)=v_spike(i);
    end
    end
end


v_time=v_time-1;

pps_spike=find(abs(pps)>=.008);              %Find PPS spike
pps_ind = [];                               % PPS integer sec index
%Find start of PPS spikes
for i=1:length(pps_spike)
    if i==1
        pps_ind(end+1)=pps_spike(i);
    else
        if (pps_spike(i)-pps_spike(i-1))>fs/2
            pps_ind(end+1)=pps_spike(i);
        end
    end
end


%Find offset between PPS spike and integer second
for i=1:length(v_time)
    [~,B]=min(abs(v_time(i)-pps_ind));  % match transmission with the nearest PPS (shift the tx time to PPS)
    v_offset_pps(i) = (v_time(i) - pps_ind(B))/fs;      % time offset of the tx from PPS
    num_offset(i)=t(pps_ind(B))-round(t(pps_ind(B)));   % time offset of the PPS to the integer second in the computer

end



%Transmission time
% tx_time_sec=t(v_time)-num_offset;

tx_time_sec=t(v_time)-num_offset;     % correct for time drift

tx_time_date=start_t+(tx_time_sec./(3600*24));

tx_time=horzcat(tx_time,tx_time_date);

pps_time = t_tx(pps_ind);           % PPS integer second (in t_tx)
%% Decoding
% signal containing edge locations
p_diff = [];

% identify step changes
stepup = [];
for i =2:length(irigb)
    p_diff(end+1) = irigb(i) - irigb(i-1);
    
end
p_diff = [0 p_diff];

threshold = max(p_diff)/(1.5);
uppk = find(p_diff >= threshold);
downpk = find(p_diff <= -threshold);

rm_ind = [];
for i =1:length(downpk)-1
    if (downpk(i+1)-downpk(i)) <= 2
        rm_ind(end+1) = i;
    end
end

stepdown = downpk;
stepdown(rm_ind) = [];
stepdown = stepdown+1;

rm_ind = [];
for i =1:length(uppk)-1
    if (uppk(i+1)-uppk(i)) <= 2
        rm_ind(end+1) = i;
    end
end
stepup= uppk;
stepup(rm_ind) = [];
if stepup(1) > stepdown(1)
    stepdown(1) = [];
end
if stepup(end) > stepdown(end)
    stepup(end) = [];
end
%% group pulse width
% pulse_width = unique(stepdown - stepup);
% for i = 1:length(pulse_width)-1
%    if pulse_width(i+1)-pulse_width(i) <=2
%        pulse_width(i) = pulse_width(i+1);
%    end
%     
% end
% pulse_width = unique(pulse_width);
% 

%% identify the start of each second

pulse = stepdown - stepup;
pulse = round(pulse/fs*1000);
sec_marker = [];               % position of the second markers among the spikes
for i=2:length(pulse)
    if (pulse(i) == 8)&(pulse(i-1) == 8)
        sec_marker(end+1) = i;
    end
end


%% identifiers
iden_marker = find(pulse == 8);
iden_marker = setdiff(iden_marker,sec_marker);
% get rid of the identifiers which come ahead of the first sec
iden_marker(find(iden_marker <sec_marker(1)))= [];

%% decode second/minute/hour/day/year
count = 1;
no_sec = [];
for jj = 1:length(sec_marker)-1
   p = iden_marker(find(sec_marker(jj) <iden_marker));
   p = p(1:10);
   if p(end) <sec_marker(jj+1)
       % second
       sec_code = pulse(sec_marker(jj)+1:p(1)-1);
       sec = 0;
       sec_wieght = [1 2 4 8 0 10 20 40];
       for ii = 1:length(sec_code)
           if sec_code(ii) == 2
               bin = 0;
           elseif sec_code(ii) == 5
               bin = 1;
           end
           sec = sec+sec_wieght(ii)*bin;
       end
       
       % minute
       min_code = pulse(p(1)+1:p(2)-1);
       minute = 0;
       min_wieght = [1 2 4 8 0 10 20 0 0];
       for ii = 1:length(min_code)
           if min_code(ii) == 2
               bin = 0;
           elseif min_code(ii) == 5
               bin = 1;
           end
           minute = minute+min_wieght(ii)*bin;
       end
       
       % hour
       hr_code = pulse(p(2)+1:p(3)-1);
       hr = 0;
       hr_wieght = [1 2 4 8 0 10 20 40 0];
       for ii = 1:length(hr_code)
           if hr_code(ii) == 2
               bin = 0;
           elseif hr_code(ii) == 5
               bin = 1;
           end
           hr = hr+hr_wieght(ii)*bin;
       end
       
       % day of year
       day_code = pulse(p(3)+1:p(5)-1);
       day = 0;
       day_wieght = [1 2 4 8 0 10 20 40 80 0 100 200 0 0 0 0 0 0 0];
       for ii = 1:length(day_code)
           if day_code(ii) == 2
               bin = 0;
           elseif day_code(ii) == 5
               bin = 1;
           else
               bin = 0;
           end
           day = day+day_wieght(ii)*bin;
       end
       
       % year
       year_code = pulse(p(5)+1:p(6)-1);
       year = 0;
       year_wieght = [1 2 4 8 0 10 20 40 80];
       for ii = 1:length(year_code)
           if year_code(ii) == 2
               bin = 0;
           elseif year_code(ii) == 5
               bin = 1;
           else
               bin = 0;
           end
           year = year+year_wieght(ii)*bin;
       end
       
       date = datenum(datetime("20"+num2str(year) + num2str(day),'InputFormat','uuuuDDD'));
       irigb_time(count) = date+hr/24+minute/(24*60)+sec/(24*3600);
       count = count +1;    
   else
       no_sec(end+1) = jj;
       
   end
end

%% subtract only integer seconds from p_diff
sec_marker(no_sec) = [];
sec_marker(end) = [];
irigb_ind = stepup(sec_marker);    % Irigb integer second index (in t_tx)   
filter = zeros(1,length(p_diff));
for iii = 1:length(irigb_ind)
    filter(irigb_ind(iii)-50:irigb_ind(iii)+50) = ones(1,101);    
end

p_diff_n = p_diff.*filter;


%% Match PPS and IRIG B
% irigb_ind/ pps_ind
% pps_time/ irigb_time
% check the difference between pps_ind and irigb_ind (both are with respect to t_tx)
% match indices of the corresponding PPS and Second Marker together/ eliminate those without pairs
rm_ind = [];
 for ii  = 1:length(pps_ind)
    if min(abs(pps_ind(ii)-irigb_ind)) >5      % if more than 5 samples apart
        rm_ind(end+1) = ii;
    end
 end
pps_ind(rm_ind) = [];
pps_time = t_tx(pps_ind);
%% Plot spikes and pulses
figure(1)
clf
plot(t_tx,pps)
hold on
plot(t_tx,p_diff_n)
scatter(t_tx(irigb_ind),p_diff_n(irigb_ind))
scatter(t_tx(pps_ind),pps(pps_ind))
datetick('x')
axis tight

%% linear regression
ind1 = 20;
ind2 = length(pps_time);
clock_drift = (pps_time(ind1:ind2)-irigb_time(ind1:ind2))*3600*24*1000; % ms
b = polyfit(t(pps_ind(ind1:ind2)),clock_drift,1);
slope = b(1);   % ms/sec
abs_drift = (t(pps_ind(ind2))-t(pps_ind(ind1)))*slope


% Plot time difference
figure(2)
clf
plot(t(pps_ind(ind1:ind2)),clock_drift)
title(sprintf('Time Difference (PPS - IRIG-B) \n drift rate = %.6f ms/sec',slope))
ylabel('msec','fontsize',15)
grid on
xlabel('Time of the start of the file (second)')