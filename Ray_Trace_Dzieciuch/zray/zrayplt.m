function zrayplt( zry, swplt, zbfr);
%function zrayplt( zry, swplt, zbfr);
global zenv

if ~exist( 'swplt', 'var')
  swplt=0;   %all plots
  %swplt=1;  %raypath plot only
  %swplt=2;  %ray kernel
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nray=length(zry);
srcdep=zry(1).rz(1);
rstart=zry(1).rx(1);
rend=zry(1).finalx;
rgsv=zry(1).rx;

%sound speed at rcvr
[c0, cz0]=zrayc(srcdep, rstart);
zzr=zenv.ssz;
for iz=1:length(zzr)
  czr(iz)=zrayc(zzr(iz),rend);
end
[caxr, cindr]=min( czr); zaxr=zzr(cindr);
clear cindr

%bathymetry
bthz=ppval(zenv.ppr, zenv.bthr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist( 'zbfr', 'var')
  zbfr=zaxr;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if swplt==0

  %angle-depth diagram
  clf
  plot( [zry.finala], -[zry.finalz], 'bo');
  hold on
  plot( [zry.finala], -[zry.finalz], 'b');
  plot( [zry.theta0], -srcdep*ones(size([zry.theta0])), 'r');
  hold off
  grid on
  xlabel('arrival-angle (degrees)');
  ylabel('depth (m)');
  pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %timefront (depth v. travel-time)
  clf
  plot( [zry.finalt], -[zry.finalz], 'b');
  hold on
  plot( [zry.finalt], -[zry.finalz], 'bo');
  hold off
  grid on
  xlabel('travel-time (s)');
  ylabel('depth (m)');
  pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %angle v. travel-time
  clf
  plot( [zry.finalt], [zry.finala], 'b');
  hold on
  plot( [zry.finalt], [zry.finala], 'bo');
  hold off
  grid on
  xlabel('travel-time (s)');
  ylabel('arrival-angle (degrees)');
  pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %depth v. launch-angle
  clf
  plot( [zry.theta0], -[zry.finalz], 'b');
  hold on
  plot( [zry.theta0], -[zry.finalz], 'bo');
  hold off
  grid on
  xlabel('launch-angle (degrees)');
  ylabel('depth (m)');
  pause
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %travel-time v. launch-angle
  clf
  plot( [zry.finalt], [zry.theta0], 'b');
  hold on
  plot( [zry.finalt], [zry.theta0], 'bo');
  hold off
  grid on
  xlabel('travel-time (s)');
  ylabel('launch-angle (degrees)');
  pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %arrival-angle v. launch-angle
  clf
  plot( [zry.theta0], [zry.finala], 'b');
  hold on
  plot( [zry.theta0], [zry.finala], 'bo');
  hold off
  grid on
  xlabel('launch-angle (degrees)');
  ylabel('arrival-angle (degrees)');
  pause
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %turning-points v. launch-angle
  clf
  plot( [zry.theta0], [zry.nutp], 'bo');
  hold on
  plot( [zry.theta0], [zry.nltp], 'rx');
  plot( [zry.theta0], [zry.nmltp], 'k*');
  plot( [zry.theta0], [zry.nutp], 'b');
  plot( [zry.theta0], [zry.nltp], 'r');
  plot( [zry.theta0], [zry.nmltp], 'k');
  hold off
  grid on
  xlabel('launch-angle (degrees)');
  ylabel('number of turning-points');
  legend('upper', 'lower', 'mixed layer')
  pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %surface/bottom bounces v. launch-angle
  clf
  plot( [zry.theta0], [zry.nsb], 'bo');
  hold on
  plot( [zry.theta0], [zry.nbb], 'rx');
  plot( [zry.theta0], [zry.nsb], 'b');
  plot( [zry.theta0], [zry.nbb], 'r');
  hold off
  grid on
  xlabel('launch-angle (degrees)');
  ylabel('number of bounces');
  legend('surface', 'bottom')
  pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %linear beamformer
  bft=[zry.finalt] + ([zry.finalz]-zbfr).*sin(pi*[zry.finala]/180)/caxr;
  bfa=180*acos(caxr*cos(pi*[zry.finala]/180)./[zry.finalc])/pi;
  bfa=bfa.*[zry.finala]./abs([zry.finala]);
  %fprintf(1,'bft: %10.5f bfa: %9.5f\n', [bft; bfa]); 

  clf
  ind=find([zry.theta0]>0);
  plot( bft(ind), bfa(ind), 'bo');
  hold on
  ind=find([zry.theta0]<0);
  plot( bft(ind), bfa(ind), 'ro');
  grid on; zoom on
  xlabel('axial travel-time (s)');
  ylabel('axial arrival-angle (degrees)');
  title('linear beamformer')
  clear ind bft bfa
  pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %turning-point filter
  tpfa=180*acos(caxr*cos(pi*[zry.finala]/180)./[zry.finalc])/pi;
  tpfa=tpfa.*[zry.finala]./abs([zry.finala]);
  tpft=zeros(size(tpfa));
  for ir=1:nray
    sv2=1./czr.^2 - (cos(pi*tpfa(ir)/180)/caxr).^2;
    sv=zeros(size(sv2));
    ind=find( sv2>=0.0);
    sv(ind)=sqrt(sv2(ind));
    if( tpfa(ir)<0.0)  sv=-sv;  end
    tau=cumtrapz(zzr,sv);
    td0=interp1(zzr,tau,zbfr,'linear');
    tdz=interp1(zzr,tau,zry(ir).finalz,'linear')-td0;
    tpft(ir)=zry(ir).finalt+tdz;
  end
  clear ir sv2 sv ind tau td0 tdz
  %fprintf(1,'tpft: %10.5f tpfa: %9.5f\n', [tpft; tpfa]); 

  clf
  ind=find([zry.theta0]>0);
  plot( tpft(ind), tpfa(ind), 'bo');
  hold on; grid on; zoom on;
  ind=find([zry.theta0]<0);
  plot( tpft(ind), tpfa(ind), 'ro');
  xlabel('axial travel-time (s)');
  ylabel('axial arrival-angle (degrees)');
  title('turning-point filter')
  clear ind tpft tpfa
  pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %lower turning-point depth v. launch-angle
  iltp=find( [zry.nltp]>0);
  if ~isempty(iltp)
    cltp=[zry.finalc]./cos(pi*[zry.finala]/180);
    zind=find( zzr>=zaxr);
    zltp=interp1( czr(zind), zzr(zind), cltp, 'spline');

    clf
    plot( [zry.theta0], zltp, 'b')
    hold on
    plot( [zry.theta0], zltp, 'bo')
    grid on; zoom on
    xlabel('launch-angle (degrees)');
    ylabel('lower turning-point depth (m)');
    clear zind cltp zltp
    pause
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%ray momentum
  %clf
  %hold on
  %warning off
  %for ir=1:1:nray
  %  plot( zry(ir).rx/1000, zry(ir).rp)
  %end
  %hold off;
  %grid on; zoom on
  %xlabel('range (km)');
  %ylabel('ray momentum');
  %warning on
  %clear ir
  %pause

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot raypaths
ntoshow=50;
if nray>ntoshow
  inc=floor(nray/(ntoshow-1));
  if (inc==0) inc=1; end
  irs=1:inc:nray;
else
  irs=1:1:nray;
end
clear inc ntoshow

clf
hold on
for ir=irs
  rx=zry(ir).rx;
  rz=zry(ir).rz;
  plot( rx/1000, -rz)

  %mark lower turning point
  ind=find( rz(2:end-1)>rz(1:end-2) & rz(2:end-1)>rz(3:end)); ind=ind+1;
  plot( rx(ind)/1000, -rz(ind), 'r*')
  %mark upper turning point
  ind=find( rz(2:end-1)<rz(1:end-2) & rz(2:end-1)<rz(3:end)); ind=ind+1;
  plot( rx(ind)/1000, -rz(ind), 'g*')

  if isfield( zry, 'eryt')
    eryt=zry(ir).eryt;
  else
    eryt=zry(ir).finalt;
    eryt=[];
  end

  zutp=mean(zry(ir).zutp);
  zltp=mean(zry(ir).zltp);

  nsb=zry(ir).nsb;
  nbb=zry(ir).nbb;

  fprintf('rayid: %4d', zry(ir).rayid);
  %fprintf(' t: %7.3f eryt: %7.3f zltp: %6.1f zutp: %6.1f',...
  %          zry(ir).finalt, eryt, zltp, zutp);
  fprintf(' eryt: %7.3f zltp: %6.1f zutp: %6.1f', eryt, zltp, zutp);
  fprintf(' nsb: %4d nbb: %4d ', nsb, nbb);
  fprintf(' theta0: %7.3f\n', zry(ir).theta0);
end
clear ir ind rx rz
plot( zenv.bthr/1000, -bthz, 'k')
plot( rgsv/1000, zeros(size(rgsv)), 'r')
grid on; zoom on
axis([rgsv(1)/1000 rgsv(end)/1000 -max(bthz)-200 500])
xlabel('range (km)');
ylabel('depth (m)');
drawnow

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%list rays
for ir=1:nray
  fprintf('rayid: %4d', zry(ir).rayid);
  fprintf(' theta0: %7.3f utp: %2d ltp: %2d',...
            zry(ir).theta0, zry(ir).nutp, zry(ir).nltp);
  fprintf(' sb: %2d bb: %2d mltp: %2d',...
            zry(ir).nsb, zry(ir).nbb, zry(ir).nmltp);
  fprintf(' t: %7.3f a: %6.2f z: %6.1f\n',...
            zry(ir).finalt, zry(ir).finala, zry(ir).finalz);
end
clear ir 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


