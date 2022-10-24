template<typename T>
void smDensitySym(SMX smx,int pi,int nSmooth,int *pList,float *fList)
{
	T fNorm,ih2,r2,rs,ih;
	int i,pj;
	KD kd = smx->kd;

	ih = 1.0/GETSMOOTH(T,pi); #h1
	ih2 = ih*ih;
  	fNorm = 0.5*M_1_PI*ih*ih2; #fac*0.5

	for (i=0;i<nSmooth;++i) {
    	pj = pList[i];
		r2 = fList[i]*ih2; #q**2, flist = r_ij**2
                rs = 2.0 - sqrt(r2); #tmp2
		if (r2 < 1.0) rs = (1.0 - 0.75*rs*r2); 
		else rs = 0.25*rs*rs*rs; #val
		if(rs<0 && !smx->warnings) {
		  fprintf(stderr, "Internal consistency error\n");
		  smx->warnings=true;
		}
		rs *= fNorm; #val*fac
		ACCUM<T>(kd->pNumpyDen,kd->p[pi].iOrder,rs*GET<T>(kd->pNumpyMass,kd->p[pj].iOrder));
		ACCUM<T>(kd->pNumpyDen,kd->p[pj].iOrder,rs*GET<T>(kd->pNumpyMass,kd->p[pi].iOrder));
  }

}


template<typename T>
void smDensitySym(SMX smx,int pi,int nSmooth,int *pList,float *fList)
{
	T fNorm,ih2,r2,rs,ih;
	int i,pj;
	KD kd = smx->kd;

	ih = 1.0/GETSMOOTH(T,pi); #h1
	ih2 = ih*ih;
  	fNorm = M_1_PI * 21.0 / 16.0;

	for (i=0;i<nSmooth;++i) {
    	pj = pList[i];
		r2 = fList[i]*ih2; #q**2, flist = r_ij**2
        q = sqrt(r2)
                rs = 1.0 - 0.5 * q; #tmp2
                val = 0.0
		if (r2 < 2.0) val = rs*rs*rs*rs*(2.0*q+1.0); #val
		if(rs<0 && !smx->warnings) {
		  fprintf(stderr, "Internal consistency error\n");
		  smx->warnings=true;
		}
		val *= fNorm; #val*fac
		ACCUM<T>(kd->pNumpyDen,kd->p[pi].iOrder,rs*GET<T>(kd->pNumpyMass,kd->p[pj].iOrder));
		ACCUM<T>(kd->pNumpyDen,kd->p[pj].iOrder,rs*GET<T>(kd->pNumpyMass,kd->p[pi].iOrder));
  }

}