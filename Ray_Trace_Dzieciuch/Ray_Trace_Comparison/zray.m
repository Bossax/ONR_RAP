function zry=zray( theta0, srcdep, rgsv, reltol);
%function zry=zray( theta0, srcdep, rgsv, reltol);
  global zenv
  
  if nargin<4
    reltol=1e-6;
  end
  nbblimit=1000;
  nbblimit=20;
  nsblimit=1000;
  swrgsv=0;

  nray=length(theta0);

  rgsv=rgsv(:);
  rstart=rgsv(1);
  rend=rgsv(end);

  %stop tracing any ray whose absolute angle becomes greater than this
  thetasteep=50;

  clear zrayc zraydot
  [c0, cz0]=zrayc(srcdep, rstart);

  idt=1; idz=2; idp=3; idw=4; ids=5;

  for ir=1:nray
    % write out angle being launched.
    %fprintf('theta0: %7.3f', theta0(ir));

    %global stphst nstpe nstps
    %stphst=[]; nstpe=0; nstps=0;

    X0(idt)=0;
    X0(idz)=srcdep;
    X0(idp)=-sin(pi*theta0(ir)/180)./c0;
    X0(idw)=0;
    X0(ids)=0;

    xerr(idt)=0.001;
    xerr(idz)=0.1;
    xerr(idp)=0.050*(sin(pi*15/180)/1475);
    xerr(idw)=1.000;
    xerr(ids)=1.000;

    nsb=0;    %number of surface bounces
    nbb=0;    %number of bottom bounces
    nutp=0;   %number of upper turning points
    nltp=0;   %number of lower turning points
    zutp=[];  %upper turning point depths
    zltp=[];  %lower turning point depths
    %zmltp=[]; %mixed-layer turning point depths

    %options.reltol=1e-4;
    options.reltol=reltol;
    options.abstol=1e-2*xerr/(rend-rstart);
    options.initialstep=100.0;
    options.events=@zrayrflct;

    opts=odeset( 'reltol', options.reltol);
    opts=odeset( opts, 'abstol', options.abstol);
    opts=odeset( opts, 'initialstep', options.initialstep);
    opts=odeset( opts, 'events', options.events);

    rnow=rstart; Xnow=X0;
    rx=rnow; X=Xnow;
    rev=[];
    while (rnow<rend)
      warning off
      rgo=[rnow, rend];
      [rint,Xint,rv,Xv,Ip]=ode45( @zraydot, rgo, Xnow, opts);
      %[rint,Xint,rv,Xv,Ip]=ode45x( rgo, Xnow, options);
      warning on

      rnow=rint(end); Xnow=Xint(end,:);
      rx=[rx; rint(2:end)]; X=[X; Xint(2:end,:)];

      if isempty(Ip) 
        Ip=0;
      else
        Ip=Ip(end);
        rev=[rev;rv(end)];
      end
    
      if Ip==1
        %surface bounce
        Xnow(idp)=-Xnow(idp);
        nsb=nsb+1;
	zutp=[zutp; 0.0];
      elseif Ip==2
        %bottom bounce
        zbslpx=ppval(zenv.dppr, rnow);
        [cnow, cznow]=zrayc( Xnow(idz), rnow);
        anow=asin(-cnow*Xnow(idp));
        anow=-anow-2*atan(zbslpx);
        Xnow(idp)=-sin(anow)./cnow;
        nbb=nbb+1;
        zbx=ppval(zenv.ppr, rnow);
	zltp=[zltp; zbx];
      elseif Ip==3
        %upper turning point
        nutp=nutp+1;
	zutp=[zutp; Xnow(idz)];
      elseif Ip==4
        %lower turning point
        nltp=nltp+1;
	zltp=[zltp; Xnow(idz)];
      end

      %Check if too many bottom bounces
      if (nbb>=nbblimit) break; end;
      %Check if too many surface bounces
      if (nsb>=nsblimit) break; end;
      %Check if too ray is too steep
      [cnow, cznow]=zrayc( Xnow(idz), rnow);
      anow=asin(-cnow*Xnow(idp));
      if (abs(anow)>(pi*thetasteep)/180) break; end;

      options.initialstep=rx(end)-rx(end-1);
      opts=odeset( opts, 'initialstep', options.initialstep);
    end

    %find mixed-layer turning points
    mxldep=400;
    rz=X(:,idz);
    ind=find( rz(2:end-1)>rz(1:end-2) & rz(2:end-1)>rz(3:end)); ind=ind+1;
    zind=find(rz(ind)<mxldep);
    nmltp=length(zind);
    clear rz ind zind

    %find rayid
    ryidu=nutp+nsb-nmltp;
    ryidl=nltp+nbb-nmltp;
    rayid=sign(theta0(ir)).*(ryidu+ryidl);
    if (theta0(ir)==0)
      rayid=ryidu+ryidl;
    end
    clear ryidu ryidl

    nrx=length(rx);
    zry(ir).theta0=theta0(ir);   %launch angle
    zry(ir).rayid=rayid;         %ray id
    zry(ir).nutp=nutp;           %number of upper turning points
    zry(ir).nltp=nltp;           %number of lower turning points
    zry(ir).nmltp=nmltp;         %number of mixed-layer turning points
    zry(ir).nsb=nsb;             %number of surface bounces
    zry(ir).nbb=nbb;             %number of bottom bounces

    zry(ir).finalx=rx(nrx);      %finalx
    zry(ir).finalz=X(nrx,idz);   %finalz
    zry(ir).finalt=X(nrx,idt);   %finalt
    zry(ir).finalp=X(nrx,idp);   %finalp
    zry(ir).finalw=X(nrx,idw);   %finalw
    zry(ir).finals=X(nrx,ids);   %finals

    zry(ir).finalc=zrayc( zry(ir).finalz, zry(ir).finalx);
    zry(ir).finala=180*asin(-zry(ir).finalc*zry(ir).finalp)/pi;

    zry(ir).zutp=zutp;
    zry(ir).zltp=zltp;
    %zry(ir).zmltp=zmltp;

    if swrgsv==1
      rxi=rgsv;
    elseif swrgsv==-1
      rxi=rx;
    else
      rxi=[rgsv; rev];
      [rxi, rind]=sort(rxi);
      rind=find(diff(rxi)==0);
      rxi(rind+1)=[];
      clear rind
    end

    %keyboard

    zry(ir).nrx=length(rxi);
    zry(ir).rx=rxi;                                     %ray range
    zry(ir).rz=interp1( rx, X(:,idz), rxi, 'linear');   %ray depth
    zry(ir).rt=interp1( rx, X(:,idt), rxi, 'linear');   %ray travel-time
    zry(ir).rp=interp1( rx, X(:,idp), rxi, 'linear');   %ray momentum
    zry(ir).rw=-interp1( rx, X(:,idw), rxi, 'linear');  %ray weighting
    zry(ir).rs=interp1( rx, X(:,ids), rxi, 'linear');   %ray arc length
    zry(ir).dwdx=gradient(zry(ir).rw,rxi);              %ray weighting
   
    %clf
    %plot( rx(2:end), diff(-X(:,idw))./diff(rx)); 
    %hold on
    %plot( rxi,  zry(ir).dwdx, 'r')

    %zry(ir).nrx=length(rxi);
    %zry(ir).rx=rxi;                                     %ray range
    %zry(ir).rz=interp1( rx, X(:,idz), rxi, 'spline');   %ray depth
    %zry(ir).rt=interp1( rx, X(:,idt), rxi, 'spline');   %ray travel-time
    %zry(ir).rp=interp1( rx, X(:,idp), rxi, 'spline');   %ray momentum
    %zry(ir).rw=-interp1( rx, X(:,idw), rxi, 'spline');  %ray weighting
    %zry(ir).rs=interp1( rx, X(:,ids), rxi, 'spline');   %ray arc length
    %zry(ir).dwdx=gradient(zry(ir).rw,rxi);              %ray weighting
    
    %plot( rxi,  zry(ir).dwdx, 'k')

%     fprintf(' rayid: %4d', zry(ir).rayid);
%     fprintf(' nutp: %2d nltp: %2d nmltp: %2d',...
%               zry(ir).nutp, zry(ir).nltp, zry(ir).nmltp);
%     fprintf(' sb: %2d bb: %2d', zry(ir).nsb, zry(ir).nbb);
%     fprintf(' t: %7.3f a: %6.2f z: %6.1f\n',...
%               zry(ir).finalt, zry(ir).finala, zry(ir).finalz);

    if (ir==1) zry=repmat(zry,1,nray); end
  end

return

function varargout = ode45x(tspan,y0,options)
%ODE45  Solve non-stiff differential equations, medium order method.
%   [TOUT,YOUT] = ODE45(ODEFUN,TSPAN,Y0) with TSPAN = [T0 TFINAL] integrates 
%   the system of differential equations y' = f(t,y) from time T0 to TFINAL 
%   with initial conditions Y0. ODEFUN is a function handle. For a scalar T
%   and a vector Y, ODEFUN(T,Y) must return a column vector corresponding 
%   to f(t,y). Each row in the solution array YOUT corresponds to a time 
%   returned in the column vector TOUT.  To obtain solutions at specific 
%   times T0,T1,...,TFINAL (all increasing or all decreasing), use TSPAN = 
%   [T0 T1 ... TFINAL].     
%   
%   ODE45 is an implementation of the explicit Runge-Kutta (4,5) pair of
%   Dormand and Prince called variously RK5(4)7FM, DOPRI5, DP(4,5) and DP54.
%   It uses a "free" interpolant of order 4 communicated privately by
%   Dormand and Prince.  Local extrapolation is done.

%   Details are to be found in The MATLAB ODE Suite, L. F. Shampine and
%   M. W. Reichelt, SIAM Journal on Scientific Computing, 18-1, 1997.

% Handle solver arguments
%htspan = abs(tspan(2) - tspan(1));  
tspan = tspan(:);
ntspan = length(tspan);
t0 = tspan(1);  
tfinal = tspan(end);     
tdir = sign(tfinal - t0);

%global stphst nstpe nstps

y0 = y0(:);
neq = length(y0);
f0 = zraydot(t0,y0);

dataType='double';
% Get the error control options, and set defaults.
rtol = options.reltol;
if rtol < 100*eps(dataType) 
  rtol = 100*eps(dataType);
  warning('ode45:reltol','RelTol has been increased to %g.',rtol)
end

atol = options.abstol;
atol = atol(:);
threshold = atol / rtol;

% By default, hmax is 1/10 of the interval.
hmax = abs(0.1*(tfinal-t0));
htry = options.initialstep;

% Handle the output
refine = 4;
S = (1:refine-1) / refine;

% the event function 
teout = []; yeout = []; ieout = [];
valt=zrayrflct(t0,y0);   

t = t0;
y = y0;

% Allocate memory if we're generating output.
nout = 0;
tout = []; yout = [];
if nargout > 0
  if ntspan > 2                         % output only at tspan points
    tout = zeros(1,ntspan,dataType);
    yout = zeros(neq,ntspan,dataType);
  else                                  % alloc in chunks
    chunk = min(max(100,50*refine), refine+floor((2^13)/neq));
    tout = zeros(1,chunk,dataType);
    yout = zeros(neq,chunk,dataType);
  end
  nout = 1;
  tout(nout) = t;
  yout(:,nout) = y;  
end

% Initialize method parameters.
pow = 1/5;
A = [1/5, 3/10, 4/5, 8/9, 1, 1];
B = [
    1/5         3/40    44/45   19372/6561      9017/3168       35/384
    0           9/40    -56/15  -25360/2187     -355/33         0
    0           0       32/9    64448/6561      46732/5247      500/1113
    0           0       0       -212/729        49/176          125/192
    0           0       0       0               -5103/18656     -2187/6784
    0           0       0       0               0               11/84
    0           0       0       0               0               0
    ];
E = [71/57600; 0; -71/16695; 71/1920; -17253/339200; 22/525; -1/40];
f = zeros(neq,7,dataType);
hmin = 16*eps(t);
%hmin = 5;
absh = min(hmax, max(hmin, htry));
f(:,1) = f0;

% THE MAIN LOOP
done = false;
while ~done
  
  % By default, hmin is a small number such that t+hmin is only slightly
  % different than t.  It might be 0 if t is 0.
  hmin = 16*eps(t);
  %hmin = 5;
  absh = min(hmax, max(hmin, absh));    % couldn't limit absh until new hmin
  h = tdir*absh;
  
  % Stretch the step if within 10% of tfinal-t.
  if 1.1*absh >= abs(tfinal - t)
    h = tfinal - t;
    absh = abs(h);
    done = true;
  end
  
  % LOOP FOR ADVANCING ONE STEP.
  nofailed = true;                      % no failed attempts
  while true
    hA = h*A;
    hB = h*B;
    f(:,2) = zraydot(t+hA(1),y+f*hB(:,1));
    f(:,3) = zraydot(t+hA(2),y+f*hB(:,2));
    f(:,4) = zraydot(t+hA(3),y+f*hB(:,3));
    f(:,5) = zraydot(t+hA(4),y+f*hB(:,4));
    f(:,6) = zraydot(t+hA(5),y+f*hB(:,5));

    tnew = t + hA(6);
    if done
      tnew = tfinal;   % Hit end point exactly.
    end
    h = tnew - t;      % Purify h.     
    
    ynew = y + f*hB(:,6);
    f(:,7) = zraydot(tnew,ynew);
    
    % Estimate the error.
    err = absh * norm((f*E) ./ max(max(abs(y),abs(ynew)),threshold),inf);
    
    % Accept the solution only if the weighted error is no more than the
    % tolerance rtol.  Estimate an h that will yield an error of rtol on
    % the next step or the next try at taking this step, as the case may be,
    % and use 0.8 of this value to avoid failures.

    if err > rtol   % Failed step
      %nstpe=nstpe+1;
      if absh <= hmin
        warning('MATLAB:ode45:IntegrationTolNotMet',['Failure at t=%e.  ' ...
                  'Unable to meet integration tolerances without reducing ' ...
                  'the step size below the smallest value allowed (%e) ' ...
                  'at time t.'],t,hmin);
        solver_output = {};
        if (nout > 0) % produce output
          solver_output{1} = tout(1:nout).';
          solver_output{2} = yout(:,1:nout).';
          solver_output{3} = teout.';
          solver_output{4} = yeout.';
          solver_output{5} = ieout.';
        end    
        if nargout > 0
          varargout = solver_output;
        end  
        return;
      end
      
      if nofailed
        nofailed = false;
        absh = max(hmin, absh * max(0.1, 0.8*(rtol/err)^pow));
      else
        absh = max(hmin, 0.5 * absh);
      end
      h = tdir * absh;
      done = false;
    else           % Successful step
      %nstps=nstps+1;
      %stphst=[stphst; absh;];
      break;
    end
  end
  
  % Successful step
  [te,ye,ie,valt,stop] = odezerox(valt,t,y,tnew,ynew,t0,h,f);
  if ~isempty(te)
    if (nargout > 2)
      teout = [teout, te];
      yeout = [yeout, ye];
      ieout = [ieout, ie];
    end

    if stop % Stop on a terminal event.               
      % Adjust the interpolation data to [t te(end)].   
      % Update the derivatives using the interpolating polynomial.
      taux = t + (te(end) - t)*A;        
      [ignore,f(:,2:7)] = ntrp45x(taux,t,y,h,f);        
        
      tnew = te(end);
      ynew = ye(:,end);
      h = tnew - t;
      done = true;
    end
  end

  tref = t + (tnew-t)*S;
  nout_new = refine;
  tout_new = [tref, tnew];
  yout_new = [ntrp45x(tref,t,y,h,f), ynew];
  
  if nout_new > 0
    oldnout = nout;
    nout = nout + nout_new;
    if nout > length(tout)
      tout = [tout, zeros(1,chunk,dataType)];  % requires chunk >= refine
      yout = [yout, zeros(neq,chunk,dataType)];
    end
    idx = oldnout+1:nout;        
    tout(idx) = tout_new;
    yout(:,idx) = yout_new;
  end  
  
  if done
    break
  end

  % If there were no failures compute a new h.
  if nofailed
    % Note that absh may shrink by 0.8, and that err may be 0.
    temp = 1.25*(err/rtol)^pow;
    if temp > 0.2
      absh = absh / temp;
    else
      absh = 5.0*absh;
    end
  end
  
  % Advance the integration one step.
  t = tnew;
  y = ynew;
  f(:,1) = f(:,7);  % Already have f(tnew,ynew)
end

solver_output = {};
if (nout > 0) % produce output
  solver_output{1} = tout(1:nout).';
  solver_output{2} = yout(:,1:nout).';
  solver_output{3} = teout.';
  solver_output{4} = yeout.';
  solver_output{5} = ieout.';
end    

if nargout > 0
  varargout = solver_output;
end  



function [tout,yout,iout,vnew,stop] = ...
    odezerox(v,t,y,tnew,ynew,t0,varargin)
%ODEZERO Locate any zero-crossings of event functions in a time step.
%   ODEZERO is an event location helper function for the ODE Suite.  ODEZERO
%   uses Regula Falsi and information passed from the ODE solver to locate
%   any zeros in the half open time interval (T,TNEW] of the event functions
%   coded in eventfun.
%   

% Initialize.
tol = 128*max(eps(t),eps(tnew));
tol = min(tol, abs(tnew - t));
tout = [];
yout = [];
iout = [];
tdir = sign(tnew - t);
stop = 0;
rmin = realmin;

% Set up tL, tR, yL, yR, vL, vR, isterminal and direction.
tL = t;
yL = y;
vL = v;
[vnew,isterminal,direction] = zrayrflct(tnew,ynew);
if isempty(direction)
  direction = zeros(size(vnew));   % zeros crossings in any direction
end
tR = tnew;
yR = ynew;
vR = vnew;

% Initialize ttry so that we won't extrapolate if vL or vR is zero.
ttry = tR;

% Find all events before tnew or the first terminal event.
while true
  
  lastmoved = 0;
  while true
    % Events of interest shouldn't have disappeared, but new ones might
    % be found in other elements of the v vector.
    indzc = find((sign(vL) ~= sign(vR)) & (direction .* (vR - vL) >= 0));
    if isempty(indzc)
      if lastmoved ~= 0
        error('MATLAB:odezerox:LostEvent',...
              'odezerox: an event disappeared (internal error)');
      end
      return;
    end
    
    % Check if the time interval is too short to continue looking.
    delta = tR - tL;
    if abs(delta) <= tol
      break;
    end
    
    if (tL == t) && any(vL(indzc) == 0 & vR(indzc) ~= 0)
      ttry = tL + tdir*0.5*tol;
      
    else
      % Compute Regula Falsi change, using leftmost possibility.
      change = 1;
      for j = indzc(:)'
        % If vL or vR is zero, try using old ttry to extrapolate.
        if vL(j) == 0
          if (tdir*ttry > tdir*tR) && (vtry(j) ~= vR(j))
            maybe = 1.0 - vR(j) * (ttry-tR) / ((vtry(j)-vR(j)) * delta);
            if (maybe < 0) || (maybe > 1)
              maybe = 0.5;
            end
          else
            maybe = 0.5;
          end
        elseif vR(j) == 0.0
          if (tdir*ttry < tdir*tL) && (vtry(j) ~= vL(j))
            maybe = vL(j) * (tL-ttry) / ((vtry(j)-vL(j)) * delta);
            if (maybe < 0) || (maybe > 1)
              maybe = 0.5;
            end
          else
            maybe = 0.5;
          end
        else
          maybe = -vL(j) / (vR(j) - vL(j)); % Note vR(j) ~= vL(j).
        end
        if maybe < change
          change = maybe;
        end
      end
      change = change * abs(delta);

      % Enforce minimum and maximum change.
      change = max(0.5*tol, min(change, abs(delta) - 0.5*tol));

      ttry = tL + tdir * change;
    end
    
    % Compute vtry.
    ytry = ntrp45x(ttry,t,y,varargin{:});
    vtry = zrayrflct(ttry,ytry);

    % Check for any crossings between tL and ttry.
    indzc = find((sign(vL) ~= sign(vtry)) & (direction .* (vtry - vL) >= 0));
    if ~isempty(indzc)
      % Move right end of bracket leftward, remembering the old value.
      tswap = tR; tR = ttry; ttry = tswap;
      yswap = yR; yR = ytry; ytry = yswap;
      vswap = vR; vR = vtry; vtry = vswap;
      % Illinois method.  If we've moved leftward twice, halve
      % vL so we'll move closer next time.
      if lastmoved == 2
        % Watch out for underflow and signs disappearing.
        maybe = 0.5 * vL;
        i = find(abs(maybe) >= rmin);
        vL(i) = maybe(i);
      end
      lastmoved = 2;
    else
      % Move left end of bracket rightward, remembering the old value.
      tswap = tL; tL = ttry; ttry = tswap;
      yswap = yL; yL = ytry; ytry = yswap;
      vswap = vL; vL = vtry; vtry = vswap;
      % Illinois method.  If we've moved rightward twice, halve
      % vR so we'll move closer next time.
      if lastmoved == 1
        % Watch out for underflow and signs disappearing.
        maybe = 0.5 * vR;
        i = find(abs(maybe) >= rmin);
        vR(i) = maybe(i);
      end
      lastmoved = 1;
    end
  end

  j = ones(1,length(indzc));
  tout = [tout, tR(j)];
  yout = [yout, yR(:,j)];
  iout = [iout, indzc'];
  if any(isterminal(indzc))
    if tL ~= t0
      stop = 1;
    end
    break;
  elseif abs(tnew - tR) <= tol
    %  We're not going to find events closer than tol.
    break;
  else
    % Shift bracket rightward from [tL tR] to [tR+0.5*tol tnew].
    ttry = tR; ytry = yR; vtry = vR;
    tL = tR + tdir*0.5*tol;
    yL = ntrp45x(tL,t,y,varargin{:});
    vL = zrayrflct(tL,yL);
    tR = tnew; yR = ynew; vR = vnew;
  end
end



function [yinterp,ypinterp] = ntrp45x(tinterp,t,y,h,f)
%NTRP45  Interpolation helper function for ODE45.
%   YINTERP = NTRP45(TINTERP,T,Y,TNEW,YNEW,H,F,IDX) uses data computed in ODE45
%   to approximate the solution at time TINTERP.  TINTERP may be a scalar 
%   or a row vector. 
%   [YINTERP,YPINTERP] = NTRP45(TINTERP,T,Y,TNEW,YNEW,H,F,IDX) returns also the
%   derivative of the polynomial approximating the solution. 
%
%   See also ODE45, DEVAL.

%   Mark W. Reichelt and Lawrence F. Shampine, 6-13-94
%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 1.13.4.6 $  $Date: 2005/04/18 22:11:53 $

BI = [
    1       -183/64     37/12       -145/128
    0       0       0       0
    0       1500/371    -1000/159   1000/371
    0       -125/32     125/12      -375/64 
    0       9477/3392   -729/106    25515/6784
    0       -11/7       11/3        -55/28
    0       3/2     -4      5/2
    ];

s = (tinterp - t)/h;  
yinterp = y(:,ones(size(tinterp))) + f*(h*BI)*cumprod([s;s;s;s]);

ypinterp = [];  
if nargout > 1
  ypinterp = f*BI*[ ones(size(s)); cumprod([2*s;3/2*s;4/3*s])];
end

