% zrayex.m
%
%compile the programs
%mex -O zrayc.c
%mex -O zraydot.c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
close all

theta0_lower=0;
theta0_upper=89.9;

seafloor_z=4728;

zw=(0:seafloor_z);                          %DEPTH PROFILE (m)
cw=cssprofile(zw, 0);                       %SS PROFILE (m/s)

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
srcdep=6;                                 %SOURCE DEPTH (m)

% choose launch angles
theta0=linspace(theta0_lower,theta0_upper,100);
% number of rays to trace
nray=length(theta0);
% fprintf('nray: %d\n', nray);

% reciever range
x_dist=15000;

%make a regularly-spaced set of ranges with interval rinc
rinc=10;
rgsv=linspace( 0, x_dist, ceil(x_dist/rinc)+1);
rgsv=rgsv(:);
clear rinc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ssr=[0:x_dist:x_dist];
ssc=repmat( cw, 1, length(ssr));

bthr=ssr;
zb=seafloor_z;
bthz=zb*ones(size(bthr));

global zenv;
zenv=zrayenv( zw, ssr, ssc, bthr, bthz);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%raytrace accuracy
rtol=1e-6;

% cpu=cputime;
zry=zray( theta0, srcdep, rgsv, rtol);
% cpu=cputime-cpu;
% fprintf(1,'%d rays took %.1fs, %.5f rays/sec, %.4f ray-Mm/sec.\n',...
% nray, cpu, (nray/cpu), (nray/cpu)*((rgsv(end)-rgsv(1))/1000000));
% clear cpu







%%%%%Find theta0 limits%%%%%

%First Bottom Bounce
for i=1:length(zry)
    if isempty(zry(i).zltp)~=1
        break
    end
end
theta_low=-theta0(i-1);
theta_high=-theta0(i);
z_low=zry(i-1).finalz;
z_high=seafloor_z+(seafloor_z-zry(i).finalz);


%Newton-Raphson Method
theta_new=theta_low+((seafloor_z)/((z_low-z_high)/(theta_low-theta_high)))-((z_low)/((z_low-z_high)/(theta_low-theta_high)));

theta1=theta_low;
z1=z_low;

z2=1000000;
while abs(z2-seafloor_z)>0.0001
    
zry=zray( theta_new, srcdep, rgsv, rtol);

if isempty(zry.zltp)==1
    z2=zry.finalz;
elseif isempty(zry.zltp)~=1
    z2=seafloor_z+(seafloor_z-zry.finalz);
end

theta2=theta_new;

theta_new=theta2+((seafloor_z)/((z2-z1)/(theta2-theta1)))-((z2)/((z2-z1)/(theta2-theta1)));

theta1=theta2;
z1=z2;

end

tt=zry.finalt;


% figure(1)
% zrayplt(zry)
% pause


