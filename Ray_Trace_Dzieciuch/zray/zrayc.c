#include "mex.h"
#include "matrix.h"
#include "math.h"

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[])
{
  static int initialize=1;
  static int nr, nd, ndm1;
  static double *ssr, *ssz, *ppc, *dppc;
  const mxArray *zenv;
  mxArray *pa;

  double z, x;

  int ir, iz;
  int ii, id, ik;
  double u, fhi, flo;
  double *vp, *vzp, cf, czf;
  double clo, czlo, chi, czhi;
  double v=0, vz=0;

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

  z = mxGetScalar(prhs[0]);
  x = mxGetScalar(prhs[1]);

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
      v = v*z + cf;
      vz = vz*z + czf;
    }
  } else {
    for( ii=0; ii<4; ii++) {
      id=iz+ii*ndm1;
      ik=4*ndm1*ir + id;
      clo=ppc[ik];
      czlo=dppc[ik];
      cf=flo*clo;
      czf=flo*czlo;
      v = v*z + cf;
      vz = vz*z + czf;
    }
  }

  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  vp = mxGetPr(plhs[0]);
  *vp=v;
  if (nlhs>1) {
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    vzp = mxGetPr(plhs[1]);
    *vzp=vz;
  }

}

