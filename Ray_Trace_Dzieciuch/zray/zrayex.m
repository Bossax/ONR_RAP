% zrayex.m
%
%compile the programs
%mex -O zrayc.c
%mex -O zraydot.c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zw=[0:5000];
cw=cssprofile(zw, 0);

clf
plot( cw, -zw);
xlabel( 'sound speed (m/s)')
ylabel( 'depth (m)');
grid on
drawnow

cw=cw(:);
zw=zw(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%source depth
srcdep=200;

% choose launch angles
theta0=-10:0.5:10;
% number of rays to trace
nray=length(theta0);
fprintf('nray: %d\n', nray);

% reciever range
rend=500000;
rend=250000;

%make a regularly-spaced set of ranges with interval rinc
rinc=1000;
rgsv=linspace( 0, rend, ceil(rend/rinc)+1);
rgsv=rgsv(:);
clear rinc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ssr=[0:rend:rend];
ssc=repmat( cw, 1, length(ssr));

bthr=ssr;
zb=5000;
bthz=zb*ones(size(bthr));

global zenv;
zenv=zrayenv( zw, ssr, ssc, bthr, bthz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%raytrace accuracy
rtol=1e-6;

cpu=cputime;
zry=zray( theta0, srcdep, rgsv, rtol);
cpu=cputime-cpu;
fprintf(1,'%d rays took %.1fs, %.5f rays/sec, %.4f ray-Mm/sec.\n',...
nray, cpu, (nray/cpu), (nray/cpu)*((rgsv(end)-rgsv(1))/1000000));
clear cpu

figure(1)
zrayplt(zry)
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find eigenrays
%receiver depth
rcvdep=3000;

% max depth error to accept
dacc=50;
% max number of iterations to allow for finding eigenrays.
itrmax=5;

[zrye, zry]=zrayeig( zry, srcdep, rcvdep, rgsv, rtol, dacc, itrmax);
pause
nraye=length(zrye);

swpplt=1;
figure(1)
zrayplt(zrye, swpplt)
pause

%kernel sampling
xinc=500; zinc=50;
%kernel smoothing
lfx=0; lfz=0; % no smoothing
figure(2)
for ir=1:nraye
  zrayk( zrye, ir, xinc, zinc, lfx, lfz);
  pause
end

