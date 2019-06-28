function [ demod,demod_pos,demod_max ] = ACO_cross_correlation_plots(y,t_date)
%Create pulse replica for time series
%{
%Find windows for arrival times using estimated travel times
f_low=find(min(t_date) <= est_arrival);          % indices of est_arrival fall within the time window of interest
f_low=f_low(1);                                  % lowest index of est_arrival
f_high=find(max(t_date)>=est_arrival);      
f_high=f_high(end);                              % last index of est_arrival
%}
t_window=0.025;                 %Find transmissions that arrive +/- 25 ms

demod_max=[];                   % Output
demod_pos=[];                   % Output
    
    %%%%%%Replica pulse for cross correlation%%%%%%%
    load LFM_pulse_ideal.mat
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%Cross Correlation%%
    [xc,lag]=xcorr(y,replica_pulse);
    xc(1:find(lag==0)-1)=[];
    lag = lag(find(lag ==0):end);   % keep only positive lag
    demod=abs(hilbert(xc));         %Complex Demod of signal
    
    %threshold=mean(demod)+10*std(demod);
    threshold=median(demod)+4*std(demod);       % why??
    %%%%
    
    
    %%Windows of XCORR to look through%%
    % Loop starting from the 1st est_arrival timestamp to the last one
    for i=f_low:f_high
        f_window=find((abs(est_arrival(i)-t_date).*(3600*24))<=t_window);
        
        %Find max in window
        [d_max,d_pos]=max(demod(f_window));
        demod_max(end+1)=d_max;
        demod_pos(end+1)=f_window(1)+d_pos-1;
        
        %Ensure max is larger than threshold
        if demod_max(end)<threshold
            demod_max=[];
            demod_pos=[];
        end
        
    end
    
    
    
    
end




