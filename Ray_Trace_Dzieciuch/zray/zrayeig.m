function [zrye, zry]=zrayeig( zry, srcdep, rcvdep, rgsv, rtol, dacc, itrmax);
%function [zrye, zry]=zrayeig( zry, srcdep, rcvdep, rgsv, rtol, dacc, itrmax);
%find eigenrays

  if nargin<5
    rtol=1e-6;
    dacc=25;
    itrmax=20;
  elseif nargin<6
    dacc=25;
    itrmax=5;
  elseif nargin<7
    itrmax=5;
  end

  zry(1).eryt=[];
  zry(1).erydep=rcvdep;

  itr=0;
  while 1
    itr=itr+1;
    theta0=[zry.theta0];
    finalz=[zry.finalz];
    nray=length(zry);

    %bracket root
    ddep=finalz-rcvdep;   %depth difference from receiver
    dthdz=gradient( theta0, ddep);
    zcrs=ddep(1:end-1).*ddep(2:end);
    ind=find( zcrs<=0);

    [ig, ii]=min([abs(ddep(ind)); abs(ddep(ind+1));]);
    indx=ind+ii-1;
    neig=length(indx);

    %candidates for next iteration 
    eigang=theta0(indx)-ddep(indx).*dthdz(indx);
    for ir=1:neig %check that they are still bracketed
      if ( eigang(ir)<theta0(ind(ir)) | eigang(ir)>theta0(ind(ir)+1))
        eigang(ir)=mean([theta0(ind(ir)) theta0(ind(ir)+1)]);
      end
    end
 
    % plot final depth vs launch-angle
    clf
    plot( theta0, -finalz, 'b');
    hold on
    plot( theta0, -finalz, 'bo');
    plot( theta0, -rcvdep*ones(size(theta0)), 'k');
    if (neig > 0)
      plot( theta0(indx), -finalz(indx), 'rx');
    end
    hold off
    grid on; zoom on;
    xlabel('launch-angle (degrees)');
    ylabel('depth (m)');
    drawnow

    for ir=1:neig
      %linear beamformer
      zutp=mean(zry(indx(ir)).zutp);
      zltp=mean(zry(indx(ir)).zltp);
      finalt=zry(indx(ir)).finalt;
      finala=zry(indx(ir)).finala;
      finalc=zry(indx(ir)).finalc;
      nsb=zry(indx(ir)).nsb;
      nbb=zry(indx(ir)).nbb;
      bft=finalt + (finalz(indx(ir))-rcvdep).*sin(pi*finala/180)/finalc;

      fprintf('rayid: %4d', zry(indx(ir)).rayid);
      %fprintf(' t: %7.3f bft: %7.3f zltp: %6.1f zutp: %6.1f',...
      %          finalt, bft, zltp, zutp);
      fprintf(' bft: %7.3f zltp: %6.1f zutp: %6.1f', bft, zltp, zutp);
      fprintf(' nsb: %4d nbb: %4d ', nsb, nbb);
      fprintf(' z: %6.1f\n', finalz(indx(ir)));

      zry(indx(ir)).eryt=bft;
      zry(indx(ir)).erydep=rcvdep;
    end
    
    %successful rays 
    eind=find( abs(finalz(indx)-rcvdep)<dacc);
    zrye=zry(indx(eind));

    %check if done
    if all( abs(finalz(indx)-rcvdep)<dacc) break; end
    if itr>itrmax break; end;

    eind=find( abs(finalz(indx)-rcvdep)>=dacc);
    eigang=eigang(eind);
    if itr>itrmax/2 & rem(itr,4)==0
      eigang=eigang + randn(size(eigang))/10000;
    end

    % trace more rays 
    zryx=zray( eigang, srcdep, rgsv, rtol);
    zryx(1).eryt=[];
    zryx(1).erydep=rcvdep;
    zrye=[zrye zryx];

    %append to the existing set of rays
    nray=length(zry);
    neig=length(zryx);
    zry(nray+1:nray+neig)=zryx;
    %sort in order of increasing angle
    theta0=[zry.theta0];
    [theta0,ind]=sort(theta0);
    zry=zry(ind);

    %remove duplicate rays
    ind=find( diff(theta0)==0);
    if (length(ind) > 0)
      disp(['remove ',num2str(length(ind)),' duplicate rays']);
      zry(ind+1)=[];
    end
  end

return



