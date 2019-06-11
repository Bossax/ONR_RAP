function [tskr, xaxk, zaxk]=zrayk(zry, iry, xinc, zinc, lfxx, lfzz)
%function [tskr, xaxk, zaxk]=zrayk(zry, iry, xinc, zinc, lfxx, lfzz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global zenv

ncol=64;
cmappm=clrscl('bwr',ncol);

%zb=max(zenv.ssz);
zb=1000*ceil(max(zry(iry).rz)/1000);
zb=5400;

zb=0;
for ir=1:length(zry)
  zbg=max(zry(ir).rz);
  if (zbg>zb) zb=zbg; end
end
zb=1000*ceil(zb/1000);

rend=zry(iry).finalx;

%zinc=40;
%xinc=200;
zaxk=0:zinc:zb;
xaxk=0:xinc:rend;

rx=zry(iry).rx;
rz=zry(iry).rz;
rw=zry(iry).rw;
rs=zry(iry).rs;
dwdx=zry(iry).dwdx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dwdxk=interp1( rx, dwdx, xaxk, 'linear');
rzk=interp1( rx, rz, xaxk, 'linear');
%dwdxk=interp1( rx, dwdx, xaxk, 'spline');
%rzk=interp1( rx, rz, xaxk, 'spline');

tskr=zeros(length(zaxk), length(xaxk));
for ir=1:length(xaxk)
  delc=tripuls( zaxk-rzk(ir), 2*zinc);
  nrm=trapz(zaxk,delc);
  tskr(:,ir)=dwdxk(ir)*delc/nrm;
end

%rw(end)
%trapz( rx, dwdx)
%trapz( xaxk, dwdxk)
%trapz( xaxk, trapz( zaxk,tskr))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if lfxx~=0 | lfzz~=0

  nro=length(xaxk);
  dax=median(diff(xaxk));
  R=nro*dax;
  dk=2*pi/R;
  kxaxis=[0:nro-1]*dk;
  ks=2*pi/dax; kn=ks/2;
  ind=find(kxaxis>kn);
  kxaxis(ind)=kxaxis(ind)-ks;
  kxaxis=fftshift(kxaxis);

  nzo=length(zaxk);
  daz=median(diff(zaxk));
  H=nzo*daz;
  dk=2*pi/H;
  kzaxis=[0:nzo-1]*dk;
  ks=2*pi/daz; kn=ks/2;
  ind=find(kzaxis>kn);
  kzaxis(ind)=kzaxis(ind)-nzo*dk;
  kzaxis=fftshift(kzaxis);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  kwz=2*pi/lfzz;
  kwx=2*pi/lfxx;
  kwinz=exp(-(kzaxis/2/kwz).^2);
  kwinz=kwinz/max(kwinz(:));
  kwinx=exp(-(kxaxis/2/kwx).^2);
  kwinx=kwinx/max(kwinx(:));
  kwinz=ifftshift(kwinz);
  kwinx=ifftshift(kwinx);
  kwin=kwinz'*kwinx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  tskkk=fft( fft( tskr), [], 2);
  tskkk=kwin.*tskkk;
  tskr=ifft( ifft( tskkk, [], 2));
  tskr=real(tskr);

end

if (nargout>0)
  return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkdm=3*std(tskr(:));

clf
%subplot( 221)
imagesc( xaxk/1000, zaxk/1000, tskr, [-pkdm pkdm]);
hold on
%plot( rx/1000, rz/1000, 'k')
xlabel('x (km)');
ylabel('z (km)');
%set( gca, 'ylim', [zb 0]/1000);
set( gca, 'ylim', [0 zb]/1000);
%stt=sprintf( 'tskr.iry%d.ix%d', iry,ix);
%title(stt);
title('tskr');
grid on;
colormap(cmappm);
hd=colorbar;
set( get( hd, 'ylabel'), 'string', 'sensitivity (s^2/km^3)');

return

subplot( 222)
plot( trapz(xaxk, tskr.'), -zaxk/1000);
grid on
subplot( 223)
plot( xaxk/1000, trapz(zaxk, tskr));
hd=colorbar;
set( hd, 'visible', 'off');
set( gca, 'xlim', [0 max(xaxk)]/1000);
grid on
drawnow;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

