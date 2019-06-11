function  [demod,demod_pos_dr,demod_snr_dr,demod_pos_indr,demod_snr_indr] = ACO_cross_correlation_ideal_all(y)
% find all arrivals in a HEM audio using ideal replica
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs=24000;    
%%%%%%Replica pulse for cross correlation%%%%%%%
load LFM_pulse_ideal.mat

%%Cross Correlation%%
[xc,lag]=xcorr(y,replica_pulse);
xc(1:find(lag==0)-1)=[];    % set correlation values to zero (no negative time lag)
demod=abs(hilbert(xc));     %Complex Demod of signal (Magnitude of a complex signal)
lag = lag(find(lag ==0):end);   % keep only positive lag

    
threshold= 9.5 ;    % for indirect arrival
threshold_indr = 6.8;   % for the direct arrival
SNR=10*log10(demod/median(demod));  % / noise
% threshold2=max(SNR);
% threshold2=threshold2*(3/5);  
% 
% if threshold2>threshold
%     threshold=threshold2;
% end
n=find(SNR>threshold);            % Find all points exceeding threshold
n2=find(SNR>threshold_indr);            % Find all points exceeding threshold

%% Get each direct arrival %%%%%%%%%%%%%%
demod_separation=[];

%%Separate each pulse cross correlation%%

for i=2:length(n)
    if n(i)>(n(i-1)+(fs*2))     %check if it's within 2 seconds
        demod_separation(end+1)=i;
    end
end
%%%%


if isempty(demod_separation)
    demod_max_dr=[];
    demod_pos_dr=[];
    demod_snr_dr = [];
   else
    %%%% process each cluster
    for i=1:length(demod_separation)+1
        
        if i==1
            if n(1)-200>0
                [pks,locs]=findpeaks(demod(n(1)-200:n(demod_separation(i)-1)));
            else
                [pks,locs]=findpeaks(demod(1:n(demod_separation(i)-1)+5));
            end
            if isempty(pks)==1
                demod_max_dr(i)=demod(1);
                demod_snr_dr(i) = 10*log10(demod(1)/median(demod));
                demod_pos_dr(i)=1;
            else
                
                [~,B]=max(pks);
                
                A=find((pks.*3/2)>=pks(B)); % find the neighbor peaks which potentially belongs to the first pulse
                
                demod_max_dr(i)=pks(A(1));
                demod_snr_dr(i) = 10*log10(pks(A(1))/median(demod));
                if n(1)-200>0
                    demod_pos_dr(i)=n(1)-200+locs(A(1))-1;
                else
                    demod_pos_dr(i)=locs(A(1));
                end
            end
            
        elseif i==length(demod_separation)+1
            [pks,locs]=findpeaks(demod(n(demod_separation(i-1))-200:n(end)));
            [~,B]=max(pks);
            
            A=find((pks.*3/2)>=pks(B));
            
            demod_max_dr(i)=pks(A(1));
            demod_snr_dr(i) = 10*log10(pks(A(1))/median(demod));
            demod_pos_dr(i)=n(demod_separation(i-1))-200+locs(A(1))-1;
            
        else
            [pks,locs]=findpeaks(demod(n(demod_separation(i-1))-200:n(demod_separation(i)-1)));
            [~,B]=max(pks);
            
            A=find((pks.*3/2)>=pks(B));
            
            demod_max_dr(i)=pks(A(1));
            demod_snr_dr(i) = 10*log10(pks(A(1))/median(demod));
            demod_pos_dr(i)=n(demod_separation(i-1))-200+locs(A(1))-1;
        end
        
        
    end
end

%% Get all arrivals
demod_separation=[];

%%Separate each pulse cross correlation%%

for i=2:length(n2)
    if n2(i)>(n2(i-1)+(fs*2))     %check if it's within 2 seconds
        demod_separation(end+1)=i;
    end
end
%%%%


if isempty(demod_separation)
    demod_max_indr=[];
    demod_pos_indr=[];
    demod_snr_indr = [];
   
else
    %%%% process each cluster
    for i=1:length(demod_separation)+1
        if i==1
            if n2(1)-200>0
                [pks,locs]=findpeaks(demod(n2(1)-200:n2(demod_separation(i)-1)));
            else
                [pks,locs]=findpeaks(demod(1:n2(demod_separation(i)-1)+5));
            end
            if isempty(pks)==1
                demod_max_indr(i)=demod(1);
                demod_snr_indr(i) = 10*log10(demod(1)/median(demod));
                demod_pos_indr(i)=1;
            else
                
                [~,B]=max(pks);
                
                A=find((pks.*3/2)>=pks(B)); % find the neighbor peaks which potentially belongs to the first pulse
                
                demod_max_indr(i)=pks(A(1));
                demod_snr_indr(i) = 10*log10(pks(A(1))/median(demod));
                if n2(1)-200>0
                    demod_pos_indr(i)=n2(1)-200+locs(A(1))-1;
                else
                    demod_pos_indr(i)=locs(A(1));
                end
            end
            
        elseif i==length(demod_separation)+1
            [pks,locs]=findpeaks(demod(n2(demod_separation(i-1))-200:n2(end)));
            [~,B]=max(pks);
            
            A=find((pks.*3/2)>=pks(B));
            
            demod_max_indr(i)=pks(A(1));
            demod_snr_indr(i) = 10*log10(pks(A(1))/median(demod));
            demod_pos_indr(i)=n2(demod_separation(i-1))-200+locs(A(1))-1;
            
        else
            [pks,locs]=findpeaks(demod(n2(demod_separation(i-1))-200:n2(demod_separation(i)-1)));
            [~,B]=max(pks);
            
            A=find((pks.*3/2)>=pks(B));
            
            demod_max_indr(i)=pks(A(1));
            demod_snr_indr(i) = 10*log10(pks(A(1))/median(demod));
            demod_pos_indr(i)=n2(demod_separation(i-1))-200+locs(A(1))-1;
        end
        
    end
end
%% Screen out the direct arrival and the artifacts

rm_ind= [];

for i = 1:length(demod_pos_indr)
     [t_from_dr,I] = min(abs(demod_pos_dr-demod_pos_indr(i))); % t offset from the clear
      t_from_dr = t_from_dr/fs; % in second
      
      if t_from_dr < 2  % closer than 2 sec
          rm_ind(end+1) = i;
      end
end
demod_pos_indr(rm_ind) = [];
demod_snr_indr(rm_ind) = [];
%{
plot(SNR)
hold on
scatter(demod_pos_dr,SNR(demod_pos_dr),'r')
scatter(demod_pos_indr,SNR(demod_pos_indr),'*k')
%}
end