function zenv=zrayenv( ssz, ssr, ssc, bthr, bthz);
%function zenv=zrayenv( ssz, ssr, ssc, bthr, bthz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  nr=length(ssr);
  nd=length(ssz);
  ssz=ssz(:);
  ssr=ssr(:);

  ppc=zeros((nd-1)*4,nr);
  dppc=zeros((nd-1)*4,nr);
  for ir=1:nr
    cpp=spline(ssz,ssc(:,ir));

    [brk,cf,nbrk,ncf]=unmkpp(cpp);
    cfd=zeros(nbrk,4);
    for ib=1:nbrk
      xcfd=polyder(cf(ib,:));
      cfd(ib,:)=[ zeros(1,4-length(xcfd)) xcfd];
    end

    ppc(:,ir)=cf(:); 
    dppc(:,ir)=cfd(:); 
  end
  clear ir cpp cf ncf brk nbrk cfd ib xcfd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ppr=spline(bthr, bthz);
  [brk,cf,nbrk,ncf]=unmkpp(ppr);
  cfd=zeros(nbrk,4);
  for ib=1:nbrk
    xcfd=polyder(cf(ib,:));
    cfd(ib,:)=[ zeros(1,4-length(xcfd)) xcfd];
  end
  dppr=mkpp(brk, cfd);
  clear brk cf nbrk ncf cfd ib xcfd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  zenv.ssz=ssz;
  zenv.ssr=ssr;
  zenv.ppc=ppc;
  zenv.dppc=dppc;
  zenv.bthr=bthr;
  zenv.ppr=ppr;
  zenv.dppr=dppr;

return

