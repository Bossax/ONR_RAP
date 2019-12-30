function zenv=zrayenv( ssz, ssr, ssc, bthr, bthz);
%function zenv=zrayenv( ssz, ssr, ssc, bthr, bthz);
% inputs    ssz = depth profiles
%           ssr = start range and stop range
%           ssc = sound speed profile
%           bthr = max range
%           bthz = max depth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  nr=length(ssr);
  nd=length(ssz);   %numebr of depth intervals
  ssz=ssz(:);
  ssr=ssr(:);

    
  ppc = zeros((nd-1)*4,nr);
  dppc = zeros((nd-1)*4,nr);
  
  % loop over sound speed profiles
  for ir=1:nr
    cpp = spline(ssz,ssc(:,ir));

    [brk,cf,nbrk,ncf]=unmkpp(cpp);
    cfd=zeros(nbrk,4);
    
    % loop over the intervals
    for ib = 1:nbrk
      xcfd = polyder(cf(ib,:)); % slope function of the polynomial
      
      % store slope values
      cfd(ib,:)=[ zeros(1,4-length(xcfd)) xcfd];
    end

    ppc(:,ir)=cf(:); % store coefficients a polynomial describing a sound speed profile in a column vector
    dppc(:,ir)=cfd(:); % store slope values the polynomial                
  end
  clear ir cpp cf ncf brk nbrk cfd ib xcfd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ppr=spline(bthr, bthz);
  [brk,cf,nbrk,ncf]=unmkpp(ppr);
  cfd = zeros(nbrk,4);
  for ib=1:nbrk
    xcfd=polyder(cf(ib,:));
    cfd(ib,:)=[ zeros(1,4-length(xcfd)) xcfd];
  end
  dppr=mkpp(brk, cfd);
  clear brk cf nbrk ncf cfd ib xcfd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % return polynomial structures      
  zenv.ssz = ssz;
  zenv.ssr = ssr;
  zenv.ppc = ppc;
  zenv.dppc = dppc;
  zenv.bthr = bthr;
  zenv.ppr = ppr;
  zenv.dppr = dppr;

return

