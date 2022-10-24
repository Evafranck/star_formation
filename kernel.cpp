template<typename T>
void smDensitySym(SMX smx,int pi,int nSmooth,int *pList,float *fList)
{
	T fNorm,ih2,r2,rs,ih;
	int i,pj;
	KD kd = smx->kd;

	ih = 1.0/GETSMOOTH(T,pi);; //h1
	ih2 = ih*ih;
  	fNorm = M_1_PI * 21.0 / 16.0;

	for (i=0;i<nSmooth;++i) {
    	pj = pList[i];
		r2 = fList[i]*ih2; //q**2, flist = r_ij**2
        q = sqrt(r2)
        rs = 1.0 - 0.5 * q; //tmp2
        val = 0.0
		if (r2 < 2.0) rs = rs*rs*rs*rs*(2.0*q+1.0); //val
		if(rs<0 && !smx->warnings) {
		  fprintf(stderr, "Internal consistency error\n");
		  smx->warnings=true;
		}
		val *= fNorm; //val*fac
		ACCUM<T>(kd->pNumpyDen,kd->p[pi].iOrder,rs*GET<T>(kd->pNumpyMass,kd->p[pj].iOrder));
		ACCUM<T>(kd->pNumpyDen,kd->p[pj].iOrder,rs*GET<T>(kd->pNumpyMass,kd->p[pi].iOrder));
  }

}

template<typename T>
void smDensity(SMX smx,int pi,int nSmooth,int *pList,float *fList)
{
	T fNorm,ih2,r2,rs,ih;
	int j,pj,pi_iord ;
	KD kd = smx->kd;

	pi_iord = kd->p[pi].iOrder;
	ih = 1.0/GET<T>(kd->pNumpySmooth, pi_iord);
	ih2 = ih*ih;
	fNorm = M_1_PI*21.0 / 16.0;
	SET<T>(kd->pNumpyDen,pi_iord,0.0);
	for (j=0;j<nSmooth;++j) {
		pj = pList[j];
		r2 = fList[j]*ih2;
        q = sqrt(r2);
		rs = 1.0 - 0.5 * q;
		if (q < 2.0) rs = rs*rs*rs*rs*(2.0*q+1.0);
        if(rs<0 && !smx->warnings) {
		  fprintf(stderr, "Internal consistency error\n");
		  smx->warnings=true;
		}
		rs *= fNorm;
		ACCUM<T>(kd->pNumpyDen,pi_iord,rs*GET<T>(kd->pNumpyMass,kd->p[pj].iOrder));
	}

}

template<typename Tf, typename Tq>
void smMeanQty1D(SMX smx,int pi,int nSmooth,int *pList,float *fList)
{
	Tf fNorm,ih2,r2,rs,ih,mass,rho;
	int j,pj,pi_iord ;
	KD kd = smx->kd;

	pi_iord = kd->p[pi].iOrder;
	ih = 1.0/GET<Tf>(kd->pNumpySmooth, pi_iord);
	ih2 = ih*ih;
	fNorm = M_1_PI*21.0 / 16.0;

	SET<Tq>(kd->pNumpyQtySmoothed,pi_iord,0.0);

	for (j=0;j<nSmooth;++j) {
		pj = pList[j];
		r2 = fList[j]*ih2;
		q = sqrt(r2)
        rs = 1.0 - 0.5 * q;
		if (q < 2.0) rs = rs*rs*rs*rs*(2.0*q+1.0);
        if(rs<0 && !smx->warnings) {
		  fprintf(stderr, "Internal consistency error\n");
		  smx->warnings=true;
		}
		rs *= fNorm;
		mass=GET<Tf>(kd->pNumpyMass,kd->p[pj].iOrder);
		rho=GET<Tf>(kd->pNumpyDen,kd->p[pj].iOrder);
		ACCUM<Tq>(kd->pNumpyQtySmoothed,pi_iord,
			  rs*mass*GET<Tq>(kd->pNumpyQty,kd->p[pj].iOrder)/rho);
	}

}

template<typename Tf, typename Tq>
void smMeanQtyND(SMX smx,int pi,int nSmooth,int *pList,float *fList)
{
	Tf fNorm,ih2,r2,rs,ih,mass,rho;
	int j,k,pj,pi_iord ;
	KD kd = smx->kd;

	pi_iord = kd->p[pi].iOrder;
	ih = 1.0/GET<Tf>(kd->pNumpySmooth, pi_iord);
	ih2 = ih*ih;
	fNorm = M_1_PI*21.0 / 16.0;

	for(k=0;k<3;++k)
		SET2<Tq>(kd->pNumpyQtySmoothed,pi_iord,k,0.0);

	for (j=0;j<nSmooth;++j) {
		pj = pList[j];
		r2 = fList[j]*ih2;
		q = sqrt(r2)
        rs = 1.0 - 0.5 * q;
		if (q < 2.0) rs = rs*rs*rs*rs*(2.0*q+1.0);
        if(rs<0 && !smx->warnings) {
		  fprintf(stderr, "Internal consistency error\n");
		  smx->warnings=true;
		}
		rs *= fNorm;
		mass=GET<Tf>(kd->pNumpyMass,kd->p[pj].iOrder);
		rho=GET<Tf>(kd->pNumpyDen,kd->p[pj].iOrder);
		for(k=0;k<3;++k) {
			ACCUM2<Tq>(kd->pNumpyQtySmoothed,pi_iord,k,
			    rs*mass*GET2<Tq>(kd->pNumpyQty,kd->p[pj].iOrder,k)/rho);
		}
	}

}

template<typename Tf, typename Tq>
void smCurlQty(SMX smx,int pi, int nSmooth,int *pList,float *fList)
{
	Tf fNorm,ih2,r2,r,rs,q2,q,ih,mass,rho, dqty[3], qty_i[3];
	int j,k,pj,pi_iord, pj_iord;
	KD kd = smx->kd;
	Tf curl[3], x,y,z,dx,dy,dz;

	pi_iord = kd->p[pi].iOrder;
	ih = 1.0/GET<Tf>(kd->pNumpySmooth, pi_iord);
	ih2 = ih*ih;
	fNorm = M_1_PI*21.0 / 16.0;

	for(k=0;k<3;++k) {
		SET2<Tq>(kd->pNumpyQtySmoothed, pi_iord, k, 0.0);
		qty_i[k] = GET2<Tq>(kd->pNumpyQty, pi_iord, k);
	}

	x = GET2<Tf>(kd->pNumpyPos, pi_iord, 0);
	y = GET2<Tf>(kd->pNumpyPos, pi_iord, 1);
	z = GET2<Tf>(kd->pNumpyPos, pi_iord, 2);

	for (j=0;j<nSmooth;++j) {
		pj = pList[j];
		pj_iord = kd->p[pj].iOrder;
		dx = x - GET2<Tf>(kd->pNumpyPos, pj_iord, 0);
		dy = y - GET2<Tf>(kd->pNumpyPos, pj_iord, 1);
		dz = z - GET2<Tf>(kd->pNumpyPos, pj_iord, 2);

		r2 = fList[j];
		q2 = r2*ih2;
		r = sqrt(r2);
		q = sqrt(q2);
		// Kernel gradient
		if (q < 1.0) rs = -3.0*ih + 2.25*r*ih2;
		else rs = -0.75*(2-q)*(2-q)/r;
		rs *= fNorm;

		mass=GET<Tf>(kd->pNumpyMass, pj_iord);
		rho=GET<Tf>(kd->pNumpyDen, pj_iord);

		for(k=0;k<3;++k)
			dqty[k] = GET2<Tq>(kd->pNumpyQty, pj_iord, k) - qty_i[k];

		curl[0] = dy * dqty[2] - dz * dqty[1];
		curl[1] = dz * dqty[0] - dx * dqty[2];
		curl[2] = dx * dqty[1] - dy * dqty[0];

		for(k=0;k<3;++k) {
			ACCUM2<Tq>(kd->pNumpyQtySmoothed, pi_iord, k, rs*curl[k]*mass/rho);
		}
	}
}

template<typename Tf, typename Tq>
void smDivQty(SMX smx,int pi, int nSmooth,int *pList,float *fList)
{
	Tf fNorm,ih2,r2,r,rs,q2,q,ih,mass,rho, div, dqty[3], qty_i[3];
	int j,k,pj,pi_iord, pj_iord;
	KD kd = smx->kd;
	Tf x,y,z,dx,dy,dz;

	pi_iord = kd->p[pi].iOrder;
	ih = 1.0/GET<Tf>(kd->pNumpySmooth, pi_iord);
	ih2 = ih*ih;
	fNorm = M_1_PI*21.0 / 16.0;

	SET<Tq>(kd->pNumpyQtySmoothed, pi_iord, 0.0);

	x = GET2<Tf>(kd->pNumpyPos, pi_iord, 0);
	y = GET2<Tf>(kd->pNumpyPos, pi_iord, 1);
	z = GET2<Tf>(kd->pNumpyPos, pi_iord, 2);

	for(k=0;k<3;++k)
		qty_i[k] = GET2<Tq>(kd->pNumpyQty, pi_iord, k);

	for (j=0;j<nSmooth;++j) {
		pj = pList[j];
		pj_iord = kd->p[pj].iOrder;
		dx = x - GET2<Tf>(kd->pNumpyPos, pj_iord, 0);
		dy = y - GET2<Tf>(kd->pNumpyPos, pj_iord, 1);
		dz = z - GET2<Tf>(kd->pNumpyPos, pj_iord, 2);

		r2 = fList[j];
		q2 = r2*ih2;
		r = sqrt(r2);
		q = sqrt(q2);
		// Kernel gradient
		if (q < 1.0) rs = -3.0*ih + 2.25*r*ih2;
		else rs = -0.75*(2-q)*(2-q)/r;
		rs *= fNorm;

		mass=GET<Tf>(kd->pNumpyMass, pj_iord);
		rho=GET<Tf>(kd->pNumpyDen, pj_iord);

		for(k=0;k<3;++k)
			dqty[k] = GET2<Tq>(kd->pNumpyQty, pj_iord, k) - qty_i[k];

		div = dx * dqty[0] + dy * dqty[1] + dz * dqty[2];

		ACCUM<Tq>(kd->pNumpyQtySmoothed, pi_iord, rs*div*mass/rho);
	}
}

template<typename Tf, typename Tq>
void smDispQtyND(SMX smx,int pi,int nSmooth,int *pList,float *fList)
{
	float fNorm,ih2,r2,rs,ih,mass,rho;
	int j,k,pj,pi_iord ;
	KD kd = smx->kd;
	float mean[3], tdiff;

	pi_iord = kd->p[pi].iOrder;
	ih = 1.0/GET<Tf>(kd->pNumpySmooth, pi_iord);
	ih2 = ih*ih;
	fNorm = M_1_PI*21.0 / 16.0;



	SET<Tq>(kd->pNumpyQtySmoothed,pi_iord,0.0);

	for(k=0;k<3;++k) {

		mean[k]=0;
	}

	// pass 1: find mean

	for (j=0;j<nSmooth;++j) {
		pj = pList[j];
		r2 = fList[j]*ih2;
		q = sqrt(r2)
        rs = 1.0 - 0.5 * q;
		if (q < 2.0) rs = rs*rs*rs*rs*(2.0*q+1.0);
        if(rs<0 && !smx->warnings) {
		  fprintf(stderr, "Internal consistency error\n");
		  smx->warnings=true;
		}
		rs *= fNorm;
		mass=GET<Tf>(kd->pNumpyMass,kd->p[pj].iOrder);
		rho=GET<Tf>(kd->pNumpyDen,kd->p[pj].iOrder);
		for(k=0;k<3;++k)
			mean[k]+=rs*mass*GET2<Tq>(kd->pNumpyQty,kd->p[pj].iOrder,k)/rho;
	}

	// pass 2: get variance

	for (j=0;j<nSmooth;++j) {
		pj = pList[j];
		r2 = fList[j]*ih2;
		q = sqrt(r2)
        rs = 1.0 - 0.5 * q;
		if (q < 2.0) rs = rs*rs*rs*rs*(2.0*q+1.0);
        if(rs<0 && !smx->warnings) {
		  fprintf(stderr, "Internal consistency error\n");
		  smx->warnings=true;
		}
		rs *= fNorm;
		mass=GET<Tf>(kd->pNumpyMass,kd->p[pj].iOrder);
		rho=GET<Tf>(kd->pNumpyDen,kd->p[pj].iOrder);
		for(k=0;k<3;++k) {
			tdiff = mean[k]-GET2<Tq>(kd->pNumpyQty,kd->p[pj].iOrder,k);
			ACCUM<Tq>(kd->pNumpyQtySmoothed,pi_iord,
				rs*mass*tdiff*tdiff/rho);
		}
	}

	// finally: take square root to get dispersion

	SET<Tq>(kd->pNumpyQtySmoothed,pi_iord,sqrt(GET<Tq>(kd->pNumpyQtySmoothed,pi_iord)));

}

template<typename Tf, typename Tq>
void smDispQty1D(SMX smx,int pi,int nSmooth,int *pList,float *fList)
{
	float fNorm,ih2,r2,rs,ih,mass,rho;
	int j,pj,pi_iord ;
	KD kd = smx->kd;
	Tq mean, tdiff;

	pi_iord = kd->p[pi].iOrder;
	ih = 1.0/GET<Tf>(kd->pNumpySmooth, pi_iord);
	ih2 = ih*ih;
	fNorm = M_1_PI*21.0 / 16.0;



	SET<Tq>(kd->pNumpyQtySmoothed,pi_iord,0.0);

	mean=0;

	// pass 1: find mean

	for (j=0;j<nSmooth;++j) {
		pj = pList[j];
		r2 = fList[j]*ih2;
		q = sqrt(r2)
        rs = 1.0 - 0.5 * q;
		if (q < 2.0) rs = rs*rs*rs*rs*(2.0*q+1.0);
        if(rs<0 && !smx->warnings) {
		  fprintf(stderr, "Internal consistency error\n");
		  smx->warnings=true;
		}
		rs *= fNorm;
		mass=GET<Tf>(kd->pNumpyMass,kd->p[pj].iOrder);
		rho=GET<Tf>(kd->pNumpyDen,kd->p[pj].iOrder);
		mean+=rs*mass*GET<Tq>(kd->pNumpyQty,kd->p[pj].iOrder)/rho;
	}

	// pass 2: get variance

	for (j=0;j<nSmooth;++j) {
		pj = pList[j];
		r2 = fList[j]*ih2;
		q = sqrt(r2)
        rs = 1.0 - 0.5 * q;
		if (q < 2.0) rs = rs*rs*rs*rs*(2.0*q+1.0);
        if(rs<0 && !smx->warnings) {
		  fprintf(stderr, "Internal consistency error\n");
		  smx->warnings=true;
		}
		rs *= fNorm;
		mass=GET<Tf>(kd->pNumpyMass,kd->p[pj].iOrder);
		rho=GET<Tf>(kd->pNumpyDen,kd->p[pj].iOrder);
		tdiff = mean-GET<Tq>(kd->pNumpyQty,kd->p[pj].iOrder);
		ACCUM<Tq>(kd->pNumpyQtySmoothed,pi_iord,rs*mass*tdiff*tdiff/rho);
	}

	// finally: take square root to get dispersion

	SET<Tq>(kd->pNumpyQtySmoothed,pi_iord,sqrt(GET<Tq>(kd->pNumpyQtySmoothed,pi_iord)));

}

