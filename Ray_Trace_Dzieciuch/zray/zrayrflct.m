function [value,isterminal,direction] = zrayrflct( r, x)

  zb=ppvalmad( r);
  value(1) =  x(2);               %ray hits surface
  isterminal(1) =  1;
  direction(1) =  -1;

  value(2) =  zb-x(2);            %ray hits bottom
  isterminal(2) =  1;
  direction(2) =  -1;

  value(3) =  x(3);               %ray upper turning point
  isterminal(3) =  1;
  direction(3) =   1;

  value(4) =  x(3);               %ray lower turning point
  isterminal(4) =  1;
  direction(4) =  -1;

return

function v=ppvalmad(xx)
  global zenv
  l = zenv.ppr.pieces;
  x = zenv.ppr.breaks;
  c = zenv.ppr.coefs;
  k = zenv.ppr.order;

  %%% for each data point, compute its breakpoint interval
  %[ignored,index] = sort([x(1:l) xx]);
  %indexa = max([find(index>l)-1; 1])

  if xx<x(1)
     xx=x(1);
  elseif xx>x(end)
    xx=x(end);
  end

  for index=1:l
    if  xx<=x(index);
      break;
    end
  end

  % now go to local coordinates ...
  xx = xx-x(index);

  % ... and apply nested multiplication:
  v = c(index,1).';
  for i=2:k
    v = xx.*v + c(index,i).';
  end
return


