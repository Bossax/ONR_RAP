function [v, vz]=zrayc( z, x)
% takes in depth and horizontal range
% return sound speed and gradient
  global zenv

  nr=length(zenv.ssr);
  nd=length(zenv.ssz);
  ndm1=nd-1;

  u = 1+(nr-1)*x/zenv.ssr(nr);
  ir = floor(u);
  if(ir<1)
    ir=1;
  elseif(ir>(nr-1))
    ir=nr-1;
  end
  fhi=u-ir;
  flo=1-fhi;

  z=real(z);
  iz=floor(1+(nd-1)*z/zenv.ssz(nd));
  if(iz<1)
    iz=1;
  elseif(iz>ndm1)
    iz=ndm1;
  end

  ind=[iz,iz+ndm1,iz+2*ndm1,iz+3*ndm1];
  clo=zenv.ppc(ind,ir);
  czlo=zenv.dppc(ind,ir);

  if fhi
    chi=zenv.ppc(ind,ir+1);
    czhi=zenv.dppc(ind,ir+1);
    cf=flo*clo+fhi*chi;
    czf=flo*czlo+fhi*czhi;
  else
    cf=flo*clo;
    czf=flo*czlo;
  end

  z=z-zenv.ssz(iz);
  v = cf(1);
  vz = czf(1);
  for i=2:4
    v = z*v + cf(i);
    vz = z*vz + czf(i);
  end

return


