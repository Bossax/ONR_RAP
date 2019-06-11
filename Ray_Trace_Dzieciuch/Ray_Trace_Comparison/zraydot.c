#include "mex.h"
#include "matrix.h"
#include "math.h"

void zraycf( double z, double x, double *c, double *dcdz)
{
  static int initialize=1;
  static int nr, nd, ndm1;
  static double *ssr, *ssz, *ppc, *dppc;
  const mxArray *zenv;
  mxArray *pa;

  int ir, iz;
  int ii, id, ik;
  double u, fhi, flo;
  double cf, czf;
  double clo, czlo, chi, czhi;

  if ( initialize) {
    zenv = mexGetVariablePtr("global", "zenv");

    pa=mxGetField( zenv, 0, "ssz");
    ssz=mxGetPr( pa);
    nd = mxGetNumberOfElements( pa);
    pa=mxGetField( zenv, 0, "ssr");
    ssr=mxGetPr( pa);
    nr = mxGetNumberOfElements( pa);
    pa=mxGetField( zenv, 0, "ppc");
    ppc=mxGetPr( pa);
    pa=mxGetField( zenv, 0, "dppc");
    dppc=mxGetPr( pa);

    ndm1=nd-1;
    initialize=0;
  }

  iz=(int)floor((nd-1)*z/ssz[nd-1]);
  if(iz<0)
    iz=0;
  else if(iz>ndm1-1)
    iz=ndm1-1;

  u=(nr-1.0)*x/ssr[nr-1];
  ir=(int)floor(u);
  if(ir<0) {
    ir=0; u=0;
  } else if (ir>nr) {
    ir=nr-1; u=nr-1;
  }
  fhi=u-ir;
  flo=1-fhi;

  z=z-ssz[iz];
  *c=0, *dcdz=0;
  if (fhi>0) {
    for( ii=0; ii<4; ii++) {
      id=iz+ii*ndm1;
      ik=4*ndm1*ir + id;
      clo=ppc[ik];
      czlo=dppc[ik];
      ik=4*ndm1*(ir+1) + id;
      chi=ppc[ik];
      czhi=dppc[ik];
      cf=flo*clo+fhi*chi;
      czf=flo*czlo+fhi*czhi;
      *c = *c*z + cf;
      *dcdz = *dcdz*z + czf;
    }
  } else {
    for( ii=0; ii<4; ii++) {
      id=iz+ii*ndm1;
      ik=4*ndm1*ir + id;
      clo=ppc[ik];
      czlo=dppc[ik];
      cf=flo*clo;
      czf=flo*czlo;
      *c = *c*z + cf;
      *dcdz = *dcdz*z + czf;
    }
  }

}

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[])
{
  double r, c, dcdz;
  double cinv, cinv2, cinv3, x22, pxinv;
  double *x, *dx;

  r = mxGetScalar(prhs[0]);
  x=mxGetPr(prhs[1]);

  plhs[0] = mxCreateDoubleMatrix(5,1,mxREAL);
  dx = mxGetPr(plhs[0]);

  zraycf( x[1], r, &c, &dcdz);
 
  cinv=1/c;
  cinv2=cinv*cinv;
  cinv3=cinv*cinv2;

  x22=x[2]*x[2];
  if (x22<cinv2)
    pxinv=1/sqrt(cinv2-x[2]*x[2]);
  else
    pxinv=0.0;

  dx[0]=cinv2*pxinv;
  dx[1]=pxinv*x[2];
  dx[2]=-dcdz*pxinv*cinv3;
  dx[3]=pxinv*cinv3;
  dx[4]=pxinv*cinv;
}




