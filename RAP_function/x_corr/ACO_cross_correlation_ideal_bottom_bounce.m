function  [demod,demod_pos_pks,demod_pos_1st_pk,demod_snr_pks,t_sep] = ACO_cross_correlation_ideal_bottom_bounce(y)
   
fs=24000;    
%%%%%%Replica pulse for cross correlation%%%%%%%
load LFM_pulse_ideal.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%Cross Correlation%%
[xc,lag]=xcorr(y,replica_pulse);
xc(1:find(lag==0)-1)=[];    % set correlation values to zero (no negative time lag)
demod=abs(hilbert(xc));     %Complex Demod of signal (Magnitude of a complex signal)
lag = lag(find(lag ==0):end);   % keep only positive lag

    
threshold= 10;    % 9 is the original number
SNR=10*log10(demod/median(demod));  % / noise
threshold2=max(SNR);
threshold2=threshold2*(4/5)+2;  


if threshold2>threshold
    threshold=threshold2;
end
n=find(SNR>threshold);            % Find all points exceeding threshold

threshold
%% Get each xcorr arrival %%
demod_separation=[];            % indicies of the first points of the clusters

demod_max = [];
demod_snr_pks = [];
demod_pos_pks = [];
demod_pos_1st_pk = [];
t_sep     = [];
%%Separate each cluster of demod multi-peak hills
% mark the first point of each demod cluster
for i=2:length(n)
    if n(i)>(n(i-1)+(fs*2))             % if > 2 sec ahead of the last point
        demod_separation(end+1)=i;      % mark as the first point of the cluster
    end
end
%%%%


if isempty(demod_separation)
    demod_max=[];
    demod_pos_pks=[];
    demod_snr_pks = [];
   
else
    % process each cluster
    for i=1:length(demod_separation)+1
        if i==1                 % the first cluster
            
            if n(1)-200>0       % the 200th index
                % find local maxima
                % cropping the demod from 200 points ahead of the first point of the first
                % cluster to the next first point of the next cluster
                [pks,locs]=findpeaks(demod(n(1)-200:n(demod_separation(i)-1))); 
            else
                
                [pks,locs]=findpeaks(demod(1:n(demod_separation(i)-1)+5));
            end
            % no peak detected
            % just pick the 1st demod point
            if isempty(pks)==1
                demod_max= vertcat(demod_max,demod([1 ;1]));
                demod_snr_pks = vertcat(demod_snr_pks,10*log10(demod(1:2)'/median(demod)));
                demod_pos_pks = vertcat(demod_pos_pks,[1; 2]);
                demod_pos_1st_pk = vertcat(demod_pos_1st_pk,1)
                t_sep = vertcat(t_sep,0);
            else
            % if peaks are detected, output all demod_pos, demod_max,demod_snr
                 [~,B]=max(pks);
                 A=find((pks.*2.6/2)>=pks(B));    % A(1) = the RAP
                 l = length(pks);
                 
                 if length(A) > 1
                     % time seperation between RAP peak and the bottom bounce
                     % must be less than 1 ms (theoretically)
                     t_sep_d = (locs(A(2))-locs(A(1)))/24000; 
                     if t_sep_d < 0.0015
                         demod_snr_pks= vertcat(demod_snr_pks,10*log10(pks(A(1:2))/median(demod)));
                         demod_max = vertcat(demod_max,pks(A(1:2)));

                        if n(1)-200>0
                            demod_pos_pks = vertcat(demod_pos_pks,n(1)-200+locs(A(1:2)-1));
                            demod_pos_1st_pk = vertcat(demod_pos_1st_pk,n(1)-200+locs(A(1)-1));
                        else
                            demod_pos_pks = vertcat(demod_pos_pks,locs(A(1:2)));
                            demod_pos_1st_pk = vertcat(demod_pos_1st_pk,locs(A(1)));
                        end
                        t_sep =  vertcat(t_sep,(locs(A(2))-locs(A(1)))/24000);
                     else
                         demod_snr_pks= vertcat(demod_snr_pks,10*log10([pks(A(1)) ; pks(A(1))]./median(demod)));
                         demod_max = vertcat(demod_max,[pks(A(1)); pks(A(1))]);

                         if n(1)-200>0
                            demod_pos_pks = vertcat(demod_pos_pks,[(n(1)-200+locs(A(1))-1);(n(1)-200+locs(A(1))-1)]);
                            demod_pos_1st_pk = vertcat(demod_pos_1st_pk,n(1)-200+locs(A(1))-1);
                         else
                            demod_pos_pks = vertcat(demod_pos_pks,[locs(A(1));locs(A(1))]);
                            demod_pos_1st_pk = vertcat(demod_pos_1st_pk,locs(A(1)));
                         end

                         t_sep =  vertcat(t_sep,0); 
                     end
                
                else
                    if A(1) == length(locs)
                         demod_snr_pks= vertcat(demod_snr_pks,10*log10([pks(A(1)) ; pks(A(1))]./median(demod)));
                         demod_max = vertcat(demod_max,[pks(A(1)); pks(A(1))]);

                         if n(1)-200>0
                            demod_pos_pks = vertcat(demod_pos_pks,[(n(1)-200+locs(A(1))-1);(n(1)-200+locs(A(1))-1)]);
                            demod_pos_1st_pk = vertcat(demod_pos_1st_pk,n(1)-200+locs(A(1))-1);
                         else
                            demod_pos_pks = vertcat(demod_pos_pks,[locs(A(1));locs(A(1))]);
                            demod_pos_1st_pk = vertcat(demod_pos_1st_pk,locs(A(1)));
                         end

                         t_sep =  vertcat(t_sep,0);
                    else
                        t_sep_d = (locs(A(1)+1)-locs(A(1)))/24000;
                        if t_sep_d < 0.0015
                             demod_snr_pks= vertcat(demod_snr_pks,10*log10([pks(A(1)) ; pks(A(1)+1)]./median(demod)));
                             demod_max = vertcat(demod_max,[pks(A(1)); pks(A(1)+1)]);

                             if n(1)-200>0
                                demod_pos_pks = vertcat(demod_pos_pks,[(n(1)-200+locs(A(1))-1);(n(1)-200+locs(A(1)+1)-1)]);
                                demod_pos_1st_pk = vertcat(demod_pos_1st_pk,n(1)-200+locs(A(1))-1);
                             else
                                demod_pos_pks = vertcat(demod_pos_pks,[locs(A(1));locs(A(1)+1)]);
                                demod_pos_1st_pk = vertcat(demod_pos_1st_pk,locs(A(1)));
                             end
                             t_sep =  vertcat(t_sep,(locs(A(1)+1)-locs(A(1)))/24000);
                        else
                             demod_snr_pks= vertcat(demod_snr_pks,10*log10([pks(A(1)) ; pks(A(1))]./median(demod)));
                             demod_max = vertcat(demod_max,[pks(A(1)); pks(A(1))]);

                             if n(1)-200>0
                                demod_pos_pks = vertcat(demod_pos_pks,[(n(1)-200+locs(A(1))-1);(n(1)-200+locs(A(1))-1)]);
                                demod_pos_1st_pk = vertcat(demod_pos_1st_pk,n(1)-200+locs(A(1))-1);
                             else
                                demod_pos_pks = vertcat(demod_pos_pks,[locs(A(1));locs(A(1))]);
                                demod_pos_1st_pk = vertcat(demod_pos_1st_pk,locs(A(1)));
                             end

                             t_sep =  vertcat(t_sep,0); 
                        end
                    end

                         
                  end
               end 
              
          
            
        % the last cluster
        elseif i==length(demod_separation)+1
                     
            [pks,locs]=findpeaks(demod(n(demod_separation(i-1))-200:n(end)));
            [~,B]=max(pks);
            A=find((pks.*3/2)>=pks(B));
            if length(A) > 1
                % time seperation between RAP peak and the bottom bounce
                % must be less than 1 ms (theoretically)
                t_sep_d = (locs(A(2))-locs(A(1)))/24000; 
                if t_sep_d < 0.0015
                    demod_snr_pks= vertcat(demod_snr_pks,10*log10(pks(A(1:2))/median(demod)));
                    demod_max = vertcat(demod_max,pks(A(1:2)));
                    demod_pos_pks = vertcat(demod_pos_pks,n(demod_separation(i-1))-200+locs(A(1:2))-1);
                    demod_pos_1st_pk = vertcat(demod_pos_1st_pk,n(demod_separation(i-1))-200+locs(A(1))-1);
                    t_sep =  vertcat(t_sep,(locs(A(2))-locs(A(1)))/24000);
                else
                    demod_snr_pks= vertcat(demod_snr_pks,10*log10([pks(A(1));pks(A(1))]/median(demod)));
                    demod_max = vertcat(demod_max,[pks(A(1));pks(A(1))]);
                    demod_pos_pks = vertcat(demod_pos_pks,n(demod_separation(i-1))-200+[locs(A(1)); locs(A(1))]-1);
                    demod_pos_1st_pk = vertcat(demod_pos_1st_pk,n(demod_separation(i-1))-200+locs(A(1))-1);
                    t_sep =  vertcat(t_sep,0);
                end
            else
                if A(1) == length(locs)
                        demod_snr_pks= vertcat(demod_snr_pks,10*log10([pks(A(1));pks(A(1))]/median(demod)));
                        demod_max = vertcat(demod_max,[pks(A(1));pks(A(1))]);
                        demod_pos_pks = vertcat(demod_pos_pks,n(demod_separation(i-1))-200+[locs(A(1)); locs(A(1))]-1);
                        demod_pos_1st_pk = vertcat(demod_pos_1st_pk,n(demod_separation(i-1))-200+locs(A(1))-1);
                        t_sep =  vertcat(t_sep,0);
                else
                        t_sep_d = (locs(A(1)+1)-locs(A(1)))/24000;
                        if t_sep_d < 0.0015
                            demod_snr_pks= vertcat(demod_snr_pks,10*log10([pks(A(1));pks(A(1)+1)]/median(demod)));
                            demod_max = vertcat(demod_max,[pks(A(1));pks(A(1)+1)]);
                            demod_pos_pks = vertcat(demod_pos_pks,n(demod_separation(i-1))-200+[locs(A(1)); locs(A(1)+1)]-1);
                            demod_pos_1st_pk = vertcat(demod_pos_1st_pk,n(demod_separation(i-1))-200+locs(A(1))-1);
                            t_sep =  vertcat(t_sep,(locs(A(1)+1)-locs(A(1)))/24000);
                        else
                            demod_snr_pks= vertcat(demod_snr_pks,10*log10([pks(A(1));pks(A(1))]/median(demod)));
                            demod_max = vertcat(demod_max,[pks(A(1));pks(A(1))]);
                            demod_pos_pks = vertcat(demod_pos_pks,n(demod_separation(i-1))-200+[locs(A(1)); locs(A(1))]-1);
                            demod_pos_1st_pk = vertcat(demod_pos_1st_pk,n(demod_separation(i-1))-200+locs(A(1))-1);
                            t_sep =  vertcat(t_sep,0);    
                        end
                        
                end
            end
        % the clusters in between    
        else
            
                [pks,locs]=findpeaks(demod(n(demod_separation(i-1))-200:n(demod_separation(i)-1)));
                [~,B]=max(pks);
                A=find((pks.*3/2)>=pks(B));
                
                if length(A) > 1
                    % time seperation between RAP peak and the bottom bounce
                    % must be less than 1 ms (theoretically)
                    t_sep_d = (locs(A(2))-locs(A(1)))/24000; 
                    if t_sep_d < 0.0015
                        demod_snr_pks= vertcat(demod_snr_pks,10*log10(pks(A(1:2))/median(demod)));
                        demod_max = vertcat(demod_max,pks(A(1:2)));
                        demod_pos_pks = vertcat(demod_pos_pks,n(demod_separation(i-1))-200+locs(A(1:2))-1);
                        demod_pos_1st_pk = vertcat(demod_pos_1st_pk,n(demod_separation(i-1))-200+locs(A(1))-1);
                        t_sep =  vertcat(t_sep,(locs(A(2))-locs(A(1)))/24000);
                    else
                        demod_snr_pks= vertcat(demod_snr_pks,10*log10([pks(A(1));pks(A(1))]/median(demod)));
                        demod_max = vertcat(demod_max,[pks(A(1));pks(A(1))]);
                        demod_pos_pks = vertcat(demod_pos_pks,n(demod_separation(i-1))-200+[locs(A(1)); locs(A(1))]-1);
                        demod_pos_1st_pk = vertcat(demod_pos_1st_pk,n(demod_separation(i-1))-200+locs(A(1))-1);
                        t_sep =  vertcat(t_sep,0);
                    end
                else
                    if A(1) == length(locs)
                        demod_snr_pks= vertcat(demod_snr_pks,10*log10([pks(A(1));pks(A(1))]/median(demod)));
                        demod_max = vertcat(demod_max,[pks(A(1));pks(A(1))]);
                        demod_pos_pks = vertcat(demod_pos_pks,n(demod_separation(i-1))-200+[locs(A(1)); locs(A(1))]-1);
                        demod_pos_1st_pk = vertcat(demod_pos_1st_pk,n(demod_separation(i-1))-200+locs(A(1))-1);
                        t_sep =  vertcat(t_sep,0);
                    else
                        t_sep_d = (locs(A(1)+1)-locs(A(1)))/24000;
                        if t_sep_d < 0.0015
                            demod_snr_pks= vertcat(demod_snr_pks,10*log10([pks(A(1));pks(A(1)+1)]/median(demod)));
                            demod_max = vertcat(demod_max,[pks(A(1));pks(A(1)+1)]);
                            demod_pos_pks = vertcat(demod_pos_pks,n(demod_separation(i-1))-200+[locs(A(1)); locs(A(1)+1)]-1);
                            demod_pos_1st_pk = vertcat(demod_pos_1st_pk,n(demod_separation(i-1))-200+locs(A(1))-1);
                            t_sep =  vertcat(t_sep,(locs(A(1)+1)-locs(A(1)))/24000);
                        else
                            demod_snr_pks= vertcat(demod_snr_pks,10*log10([pks(A(1));pks(A(1))]/median(demod)));
                            demod_max = vertcat(demod_max,[pks(A(1));pks(A(1))]);
                            demod_pos_pks = vertcat(demod_pos_pks,n(demod_separation(i-1))-200+[locs(A(1)); locs(A(1))]-1);
                            demod_pos_1st_pk = vertcat(demod_pos_1st_pk,n(demod_separation(i-1))-200+locs(A(1))-1);
                            t_sep =  vertcat(t_sep,0);    
                        end
                    end
                end
        
        end
    end

end
end