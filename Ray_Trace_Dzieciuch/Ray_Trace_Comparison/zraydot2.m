function [ dx] = zraydot( r, x)
    % planar ray tracing / commented by Boris 
  [c, dcdz]=zrayc( x(2), r);    % sound speed and gradient
  cinv=1/c;     % slowness
  cinv2=cinv*cinv; % 1/c^2
  cinv3=cinv*cinv2;% 1/c^3

  %x33=x(3)*x(3);
  %if (x33<cinv2)
  %  pxinv=1/sqrt(cinv2-x33);
  %else
  %  pxinv=0.0;
  %end
  pxinv=1/sqrt(cinv2-x(3)*x(3)); % px = horizontal ray momentum (horizontal slowness)
  pxinv=real(pxinv);

  dx(1)=cinv2*pxinv;          %time  
  dx(2)=pxinv*x(3);           %depth
  dx(3)=-dcdz*pxinv*cinv3;    % angle

  dx(4)=pxinv*cinv3;          %ray weighting
  dx(5)=pxinv*cinv;           %arclength

  dx=dx(:);
return


