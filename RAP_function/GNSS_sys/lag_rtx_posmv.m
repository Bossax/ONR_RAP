function [lsd,y_shift]= lag_rtx_posmv(x,y)
% x is RTX
% y is POSMV
% x is longer than y
% shift y and calculate a corresponding sum square difference
% find the least square difference
xl = length(x);
yl = length(y);
sx = size(x);
sy = size(y);
if sx(1) >1
    x = x';
end
if sy(1) > 1
    y=y';
end
lag = 0:(xl -yl);   % shift forward in time
ssd = zeros(1,length(lag));
for k = 1:length(lag)
   % truncate x
   x_truncated  = x(lag(k)+1:lag(k)+yl);
   ssd(k) = sum((x_truncated - y).^2);
end

[lsd,ind] =  min(ssd);

y_shift = lag(ind);

end