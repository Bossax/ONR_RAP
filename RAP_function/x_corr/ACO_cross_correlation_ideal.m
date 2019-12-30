function [demod,demod_pos,demod_max,demod_snr] = ACO_cross_correlation_ideal(y)
%Create pulse replica for time series
% return
fs=24000;

%%%%%%Replica pulse for cross correlation%%%%%%%
load LFM_pulse_ideal.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%Cross Correlation%%
[xc,lag]=xcorr(y,replica_pulse);
xc(1:find(lag==0)-1)=[];    % set correlation values to zero (no negative time lag)
demod=abs(hilbert(xc));     %Complex Demod of signal (Magnitude of a complex signal)
lag = lag(find(lag ==0):end);   % keep only positive lag

%{
figure(1)
plot(lag,xc)
axis tight
grid on
title(sprintf('Cross-Correlation: Signal length = %f \n Lag = %f',length(y),length(lag)))
xlabel('Lag')
%}

%threshold=mean(demod)+10*std(demod);
threshold=9;
%threshold=6;                %dB threshold
SNR=10*log10(demod/median(demod));  % / noise
threshold2=max(SNR);
threshold2=threshold2*(4/5);  % scaling??

if threshold2>threshold
    threshold=threshold2;
end
n=find(SNR>threshold);            %Find large values for demod xcorr
%%%%
% disp(threshold)

%% Get each xcorr arrival %%
demod_separation=[];

%%Separate each pulse cross correlation%%

for i=2:length(n)
    if n(i)>(n(i-1)+(fs*2))     %check if it's within 2 seconds
        demod_separation(end+1)=i;
    end
end
%%%%


if isempty(demod_separation)
    demod_max=[];
    demod_pos=[];
    demod_snr = [];
   
else
    for i=1:length(demod_separation)+1
        if i==1
            if n(1)-200>0
                [pks,locs]=findpeaks(demod(n(1)-200:n(demod_separation(i)-1)));
            else
                [pks,locs]=findpeaks(demod(1:n(demod_separation(i)-1)+5));
            end
            if isempty(pks)==1
                demod_max(i)=demod(1);
                demod_snr(i) = 10*log10(demod(1)/median(demod));
                demod_pos(i)=1;
            else
                
                [~,B]=max(pks);
                
                A=find((pks.*3/2)>=pks(B)); % find the neighbor peaks which potentially belongs to the first pulse
                
                demod_max(i)=pks(A(1));
                demod_snr(i) = 10*log10(pks(A(1))/median(demod));
                if n(1)-200>0
                    demod_pos(i)=n(1)-200+locs(A(1))-1;
                else
                    demod_pos(i)=locs(A(1));
                end
            end
            
        elseif i==length(demod_separation)+1
            [pks,locs]=findpeaks(demod(n(demod_separation(i-1))-200:n(end)));
            [~,B]=max(pks);
            
            A=find((pks.*3/2)>=pks(B));
            
            demod_max(i)=pks(A(1));
            demod_snr(i) = 10*log10(pks(A(1))/median(demod));
            demod_pos(i)=n(demod_separation(i-1))-200+locs(A(1))-1;
            
        else
            [pks,locs]=findpeaks(demod(n(demod_separation(i-1))-200:n(demod_separation(i)-1)));
            [~,B]=max(pks);
            
            A=find((pks.*3/2)>=pks(B));
            
            demod_max(i)=pks(A(1));
            demod_snr(i) = 10*log10(pks(A(1))/median(demod));
            demod_pos(i)=n(demod_separation(i-1))-200+locs(A(1))-1;
        end
        
    end
end
%{
figure(2)
plot(demod)
axis tight
grid on
title('Envelope Signal')
hold on
scatter(demod_pos,demod(demod_pos))
  %}  
    


end

