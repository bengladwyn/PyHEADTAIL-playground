#include <complex>
#include <cmath>
#include <stdlib.h>
#include <ablas.h>
#include <amp.h>
#include <mpfr.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_errno.h>

#include <IW2D.h>

struct params {unsigned int m; unsigned int n;
  	unsigned int N; unsigned int M;
	ap::template_1d_array< amp::campf<PRECISION> > eps1;
  	ap::template_1d_array< amp::campf<PRECISION> > mu1;
  	ap::template_1d_array< amp::campf<PRECISION> > nu2;
	ap::template_1d_array< amp::campf<PRECISION> > eps1ratio;
  	ap::template_1d_array< amp::campf<PRECISION> > mu1ratio;
  	ap::template_1d_array< amp::campf<PRECISION> > nu2ratio;
	ap::template_1d_array< amp::ampf<PRECISION> > b;
	ap::template_1d_array< amp::campf<PRECISION> > eps1m;
  	ap::template_1d_array< amp::campf<PRECISION> > mu1m;
  	ap::template_1d_array< amp::campf<PRECISION> > nu2m;
	ap::template_1d_array< amp::campf<PRECISION> > eps1mratio;
  	ap::template_1d_array< amp::campf<PRECISION> > mu1mratio;
  	ap::template_1d_array< amp::campf<PRECISION> > nu2mratio;
	ap::template_1d_array< amp::ampf<PRECISION> > bm;
	amp::ampf<PRECISION> beta;
	amp::ampf<PRECISION> kovergamma;
	memorycontainer* memory;};

// NB: eps1, eps1m, mu1 and mu1m are actually the inverse of the 'real' eps1, etc.

// // global arrays with memory of etas and chis
// extern ap::template_1d_array< amp::ampf<PRECISION> > kxmem;
// extern ap::template_1d_array< amp::campf<PRECISION> > eta1mem,eta2mem,chi1mem,chi2mem;
// extern unsigned long mem; // current number of elements of kxmem, eta1mem, etc.


/**************************************************
*** multilayerm_flat: computes matrix for field
*** matching (multiprecision) - flat chamber
***************************************************/

  void multilayerm_flat(ap::template_2d_array< amp::campf<PRECISION> >& m,
  	amp::campf<PRECISION>& kyn,
  	unsigned int N,
	ap::template_1d_array< amp::campf<PRECISION> >& eps1, 
  	ap::template_1d_array< amp::campf<PRECISION> >& mu1, 
  	ap::template_1d_array< amp::campf<PRECISION> >& nu2, 
	ap::template_1d_array< amp::campf<PRECISION> >& eps1ratio, 
  	ap::template_1d_array< amp::campf<PRECISION> >& mu1ratio, 
  	ap::template_1d_array< amp::campf<PRECISION> >& nu2ratio, 
	ap::template_1d_array< amp::ampf<PRECISION> >& b,
	amp::ampf<PRECISION> kx, amp::ampf<PRECISION> beta){
	
    /* computes the matrix M for the multilayer field matching, from mu1, eps1, b of each of the N layers
    and the horizontal wave number kx, angular frequency omega, relativistic velocity factor beta
    and wave number k.
    It also gives in output kyn, which is ky for the last layer. */

    // NB: eps1, eps1m, mu1 and mu1m are actually the inverse of the 'real' eps1, etc.
    
    amp::campf<PRECISION> kyp,kyp1,exp1,exp2,exp3,exp4;
    amp::campf<PRECISION> tmp1,tmp2,tmpplus,tmpminus,fac,fac2;
    amp::ampf<PRECISION> kx2,fac3;
    ap::template_2d_array< amp::campf<PRECISION> > mp1p;
    ap::template_2d_array< amp::campf<PRECISION> > mold;
   
    mp1p.setbounds(1,4,1,4);
    mold.setbounds(1,4,1,4);
  
    for (int i=1; i<=4; i++) {
      for (int j=1; j<=4; j++) {
        if (i==j) mold(i,j)=1;
	else mold(i,j)=0;
      }
    }
    
    m=mold;

    kx2=amp::sqr(kx);
    kyp=csqrtMP(kx2+nu2(1));
    fac3=kx/(2*beta);
    
    for (unsigned int p=1; p<=N-1; p++) {

      kyp1=csqrtMP(kx2+nu2(p+1));
      
      tmp1=kyp*nu2ratio(p)/kyp1;
      exp1=cexpMP((kyp-kyp1)*b(p));
      exp2=cexpMP(-(kyp+kyp1)*b(p));
      exp3=1/exp2;
      exp4=1/exp1;
      
      tmp2=tmp1*eps1ratio(p);
      tmpplus=(1+tmp2)/2;
      tmpminus=(1-tmp2)/2;
      mp1p(1,1)=tmpplus*exp1;
      mp1p(1,2)=tmpminus*exp2;
      mp1p(2,1)=tmpminus*exp3;
      mp1p(2,2)=tmpplus*exp4;
     
      fac=fac3*(nu2ratio(p)-1)/kyp1;
      fac2=fac*eps1(p+1);
      mp1p(1,3)=-fac2*exp1;
      mp1p(1,4)=-fac2*exp2;
      mp1p(2,3)=fac2*exp3;
      mp1p(2,4)=fac2*exp4;
      
      tmp2=tmp1*mu1ratio(p);
      tmpplus=(1+tmp2)/2;
      tmpminus=(1-tmp2)/2;
      mp1p(3,3)=tmpplus*exp1;
      mp1p(3,4)=tmpminus*exp2;
      mp1p(4,3)=tmpminus*exp3;
      mp1p(4,4)=tmpplus*exp4;
     
      fac2=fac*mu1(p+1);
      mp1p(3,1)=-fac2*exp1;
      mp1p(3,2)=-fac2*exp2;
      mp1p(4,1)=fac2*exp3;
      mp1p(4,2)=fac2*exp4;
      
      ablas::cmatrixgemm<PRECISION>(4,4,4,1,mp1p,1,1,0,mold,1,1,0,0,m,1,1);
      
      kyp=kyp1;
      mold=m;
    }
    
    kyn=kyp1;
  
    return;
  }
  

/**************************************************
*** etachi: computes eta1, eta2, chi1 and chi2 (multiprecision)
***
***************************************************/

  void etachi(amp::campf<PRECISION>& chi1,
  	amp::campf<PRECISION>& chi2,
  	amp::campf<PRECISION>& eta1,
  	amp::campf<PRECISION>& eta2,
  	unsigned int N, unsigned int M,
	ap::template_1d_array< amp::campf<PRECISION> >& eps1, 
  	ap::template_1d_array< amp::campf<PRECISION> >& mu1, 
  	ap::template_1d_array< amp::campf<PRECISION> >& nu2, 
	ap::template_1d_array< amp::campf<PRECISION> >& eps1ratio, 
  	ap::template_1d_array< amp::campf<PRECISION> >& mu1ratio, 
  	ap::template_1d_array< amp::campf<PRECISION> >& nu2ratio, 
	ap::template_1d_array< amp::ampf<PRECISION> >& b,
	ap::template_1d_array< amp::campf<PRECISION> >& eps1m, 
  	ap::template_1d_array< amp::campf<PRECISION> >& mu1m, 
  	ap::template_1d_array< amp::campf<PRECISION> >& nu2m, 
	ap::template_1d_array< amp::campf<PRECISION> >& eps1mratio, 
  	ap::template_1d_array< amp::campf<PRECISION> >& mu1mratio, 
  	ap::template_1d_array< amp::campf<PRECISION> >& nu2mratio, 
	ap::template_1d_array< amp::ampf<PRECISION> >& bm,
	amp::ampf<PRECISION> kx, amp::ampf<PRECISION> beta,
	memorycontainer* memory){
	
    /* function that computes chi1, chi2, eta1 and eta2 for a given horizontal wave number kx, from mu1, eps1, b 
    of each of the N upper layers and mu1m, eps1m, bm of each of the M lower layers, 
    and from the angular frequency omega, relativistic velocity factor beta and wave number k. */
    
    // NB: eps1, eps1m, mu1 and mu1m are actually the inverse of the 'real' eps1, etc.

    ap::template_2d_array< amp::campf<PRECISION> > mat,matprime; // 4*4 field matching matrices for upper and lower layers
    ap::template_2d_array< amp::campf<PRECISION> > pcal; // 4*4 final matrix to be inverted, and then its inverse
    amp::campf<PRECISION> kyn,kym; // ky for the last layers
    int i,l; // flags for last layers
    int info;
    unsigned int lprov;
    //matinv::matinvreport<Precision> repi;
    timeval c1,c2; // clock ticks


    // try to find kx in the table kxmem
    if (memory->mem==0) {
      lprov=0;
      memory->has_printed_warning = false; // If mem==0, then memory is either new or reset. In either case, an overflow warning has not been printed.
    }
    else lprov=locateMP(memory->kxmem,kx,memory->mem-1);
    if ( (memory->mem!=0)&&(kx==memory->kxmem(lprov)) ) {
      eta1=memory->eta1mem(lprov);eta2=memory->eta2mem(lprov);
      chi1=memory->chi1mem(lprov);chi2=memory->chi2mem(lprov);
    }
    else if ( ( (memory->mem!=0)&&(lprov>0)) && (kx==memory->kxmem(lprov-1)) ) {
      eta1=memory->eta1mem(lprov-1);eta2=memory->eta2mem(lprov-1);
      chi1=memory->chi1mem(lprov-1);chi2=memory->chi2mem(lprov-1);
    }
    else {

      /* setting bounds for matrices (beware, pcal first index has to be set to 0 instead of 1, because of the
      function matinv::cmatrixinverse */
      mat.setbounds(1,4,1,4);
      matprime.setbounds(1,4,1,4);
      pcal.setbounds(0,3,0,3);

      // compute the field matching 4*4 matrices (upper and lower layers)
      //gettimeofday(&c1,0);
      multilayerm_flat(mat, kyn, N, eps1, mu1, nu2, eps1ratio, mu1ratio, nu2ratio, b, kx, beta);
      multilayerm_flat(matprime, kym, M, eps1m, mu1m, nu2m, eps1mratio, mu1mratio, nu2mratio, bm, kx, beta);
      //gettimeofday(&c2,0);
      //std::cout << "multilayer: time= " << (c2.tv_sec-c1.tv_sec)*1.e6+(c2.tv_usec-c1.tv_usec) << " usec\n";

      i=1;l=2;

      // compute the final 4*4 matrix P and invert it
      for (int j=0; j<=3; j++) {
	pcal(0,j)=mat(i,j+1);
	pcal(1,j)=mat(i+2,j+1);
	pcal(2,j)=matprime(l,j+1);
	pcal(3,j)=matprime(l+2,j+1);
      }
      //pcal2=pcal;

      /*for (int k=1; k<=4; k++) {
	for (int j=1; j<=4; j++) {
          printf("%d %d : %s %s\n", k,j,mat(k,j).x.toDec().c_str(),mat(k,j).y.toDec().c_str());
	}
	printf("\n");
      }*/

      //matinv::cmatrixinverse(pcal,4,info,repi);
      matinv4(pcal); // ~50 times quicker but less accurate -> increase Precision
      /*amp::ampf<Precision> sum=0;
      for (unsigned int p=0;p<=1;p++) {
        for (unsigned int q=0;q<=3;q++) sum=amp::maximum(sum,amp::abscomplex((pcal(p,q)-pcal2(p,q))/pcal(p,q)));
      }
      std::cout << "difference: " << sum.toDec().c_str() << "\n";*/

      /*for (int k=0; k<=3; k++) {
	for (int j=0; j<=3; j++) {
          printf("%d %d : %s %s\n", k+1,j+1,pcal(k,j).x.toDec().c_str(),pcal(k,j).y.toDec().c_str());
	}
	printf("\n");
      }*/

      // computes chi1, chi2, eta1 and eta2 at kx
      // 4 next lines: matinv version
      chi1=pcal(0,0)*mat(i,2)+pcal(0,1)*mat(i+2,2);
      chi2=pcal(1,0)*mat(i,2)+pcal(1,1)*mat(i+2,2);
      eta1=pcal(0,2)*matprime(l,1)+pcal(0,3)*matprime(l+2,1);
      eta2=pcal(1,2)*matprime(l,1)+pcal(1,3)*matprime(l+2,1);
      //chi1=sol(0,0);chi2=sol(1,0);eta1=sol(0,1);eta2=sol(1,1);

      /*printf("kx : %s\n", kx.toDec().c_str());
      printf("chi1 : %s %s\n", chi1.x.toDec().c_str(),chi1.y.toDec().c_str());
      printf("chi2 : %s %s\n", chi2.x.toDec().c_str(),chi2.y.toDec().c_str());
      printf("eta1 : %s %s\n", eta1.x.toDec().c_str(),eta1.y.toDec().c_str());
      printf("eta2 : %s %s\n", eta2.x.toDec().c_str(),eta2.y.toDec().c_str());*/
     
      if (memory->mem <= memory->maxmem) {
        for(unsigned int k=memory->mem; k>=lprov+1; k--) {
          memory->kxmem(k)=memory->kxmem(k-1);
          memory->eta1mem(k)=memory->eta1mem(k-1);memory->eta2mem(k)=memory->eta2mem(k-1);
          memory->chi1mem(k)=memory->chi1mem(k-1);memory->chi2mem(k)=memory->chi2mem(k-1);
        }
        memory->kxmem(lprov)=kx;
        memory->eta1mem(lprov)=eta1;memory->eta2mem(lprov)=eta2;
        memory->chi1mem(lprov)=chi1;memory->chi2mem(lprov)=chi2;
        memory->mem++;
       
    } else if (!memory->has_printed_warning) {
      printf("Maximum memory capacity for eta's and chi's reached. Not saving new elements until memory is reset.\n");
      memory->has_printed_warning = true;
    }
    }

    return;
    
  }
  
/**************************************************
*** integrand: computes the complex integrand in u (std precision)
***
***************************************************/

  std::complex<double> integrand(unsigned int m, unsigned int n,
  	unsigned int N, unsigned int M,
	ap::template_1d_array< amp::campf<PRECISION> >& eps1, 
  	ap::template_1d_array< amp::campf<PRECISION> >& mu1, 
  	ap::template_1d_array< amp::campf<PRECISION> >& nu2, 
	ap::template_1d_array< amp::campf<PRECISION> >& eps1ratio, 
  	ap::template_1d_array< amp::campf<PRECISION> >& mu1ratio, 
  	ap::template_1d_array< amp::campf<PRECISION> >& nu2ratio, 
	ap::template_1d_array< amp::ampf<PRECISION> >& b,
	ap::template_1d_array< amp::campf<PRECISION> >& eps1m, 
  	ap::template_1d_array< amp::campf<PRECISION> >& mu1m, 
  	ap::template_1d_array< amp::campf<PRECISION> >& nu2m, 
	ap::template_1d_array< amp::campf<PRECISION> >& eps1mratio, 
  	ap::template_1d_array< amp::campf<PRECISION> >& mu1mratio, 
  	ap::template_1d_array< amp::campf<PRECISION> >& nu2mratio, 
	ap::template_1d_array< amp::ampf<PRECISION> >& bm,
	amp::ampf<PRECISION> u, amp::ampf<PRECISION> beta,
	amp::ampf<PRECISION> kovergamma,
	memorycontainer* memory){
    /* function that computes the integrand in alphamn for a given u (kx=k sinh(u)/gamma) 
    and given azimuthal mode numbers m and n, from mu1, eps1, b 
    of each of the N upper layers and mu1m, eps1m, bm of each of the M lower layers, 
    and from the angular frequency omega, relativistic velocity factor beta, wave number k and k/gamma. */
    
    // NB: eps1, eps1m, mu1 and mu1m are actually the inverse of the 'real' eps1, etc.
    
    amp::campf<PRECISION> chi1,chi2,eta1,eta2; // chi and eta functions of kx
    amp::ampf<PRECISION> kx;
    amp::campf<PRECISION> inte; // result in multiprecision
    int m1powm,m1pown,m1powmn;
    double x,y;
    std::complex<double> result; // result in standard double precision

    if (m % 2 ==0) m1powm=1;
    else m1powm=-1;
    if (n % 2 ==0) m1pown=1;
    else m1pown=-1;
    
    m1powmn=m1powm*m1pown;
    
    // computes kx
    kx=kovergamma*amp::sinh(u);
    
    // computes eta and chi functions at kx
    etachi(chi1, chi2, eta1, eta2, N, M, eps1, mu1, nu2, eps1ratio, mu1ratio, nu2ratio, b, eps1m, mu1m, nu2m, eps1mratio, mu1mratio, nu2mratio, bm, kx, beta, memory);
    
    // computes integrand
    inte=amp::cosh(m*u)*amp::cosh(n*u)*( chi1+m1powm*eta1+m1pown*chi2+m1powmn*eta2 );
    
    // check if real or imaginary part of inte is NaN. In that case, replace it by zero.
    x=double(amp::ampf<PRECISION>(inte.x).toDouble());
    y=double(amp::ampf<PRECISION>(inte.y).toDouble());
    if (x != x) x=0; // this is the way to check if it's a NaN (all comparisons give false with a NaN)
    if (y != y) y=0; // this is the way to check if it's a NaN (all comparisons give false with a NaN)

    result=std::complex<double>(x,y);
    //printf("%13.8e %13.8e %13.8e\n",double(amp::ampf<Precision>(u).toDouble()),x,y);
    return result;
    
  }
  
/**************************************************
*** integrand_real: computes the real part of the integrand in u (std precision)
***
***************************************************/

  double integrand_real(double x, void *p){
	
    /* function encapsulating the computation of the real part of the integrand of alphamn, 
    for gsl integration, see also function integrand */
    
    struct params *param=(struct params *)p;
    unsigned int N,M; // number of upper and lower layers
    unsigned int m,n; // indices of alphamn (azimuthal mode numbers)
    ap::template_1d_array< amp::campf<PRECISION> > eps1,eps1m,eps1ratio,eps1mratio,mu1,mu1m,mu1ratio,mu1mratio;
    ap::template_1d_array< amp::campf<PRECISION> > nu2,nu2m,nu2ratio,nu2mratio;
    ap::template_1d_array< amp::ampf<PRECISION> > b,bm;
    amp::ampf<PRECISION> beta,gamma,kovergamma,u;
    std::complex<double> inte;
    memorycontainer* memory;

    m=(param->m);
    n=(param->n);
    M=(param->M);
    N=(param->N);
    eps1=(param->eps1);
    eps1m=(param->eps1m);
    eps1ratio=(param->eps1ratio);
    eps1mratio=(param->eps1mratio);
    mu1=(param->mu1);
    mu1m=(param->mu1m);
    mu1ratio=(param->mu1ratio);
    mu1mratio=(param->mu1mratio);
    nu2=(param->nu2);
    nu2m=(param->nu2m);
    nu2ratio=(param->nu2ratio);
    nu2mratio=(param->nu2mratio);
    b=(param->b);
    bm=(param->bm);
    beta=(param->beta);
    kovergamma=(param->kovergamma);
    memory = (param->memory);
    u=x;

    inte=integrand(m, n, N, M, eps1, mu1, nu2, eps1ratio, mu1ratio, nu2ratio, b, eps1m, mu1m, nu2m, eps1mratio, mu1mratio, nu2mratio, bm, u, beta, kovergamma, memory);

    return inte.real();
    
  }

/**************************************************
*** integrand_imag: computes the imag. part of the integrand in u (std precision)
***
***************************************************/

  double integrand_imag(double x, void *p){
	
    /* function encapsulating the computation of the imaginary part of the integrand of alphamn, 
    for gsl integration, see also function integrand */
    
    struct params *param=(struct params *)p;
    unsigned int N,M; // number of upper and lower layers
    unsigned int m,n; // indices of alphamn (azimuthal mode numbers)
    ap::template_1d_array< amp::campf<PRECISION> > eps1,eps1m,eps1ratio,eps1mratio,mu1,mu1m,mu1ratio,mu1mratio;
    ap::template_1d_array< amp::campf<PRECISION> > nu2,nu2m,nu2ratio,nu2mratio;
    ap::template_1d_array< amp::ampf<PRECISION> > b,bm;
    amp::ampf<PRECISION> beta,gamma,kovergamma,u;
    std::complex<double> inte;
    memorycontainer* memory;

    m=(param->m);
    n=(param->n);
    M=(param->M);
    N=(param->N);
    eps1=(param->eps1);
    eps1m=(param->eps1m);
    eps1ratio=(param->eps1ratio);
    eps1mratio=(param->eps1mratio);
    mu1=(param->mu1);
    mu1m=(param->mu1m);
    mu1ratio=(param->mu1ratio);
    mu1mratio=(param->mu1mratio);
    nu2=(param->nu2);
    nu2m=(param->nu2m);
    nu2ratio=(param->nu2ratio);
    nu2mratio=(param->nu2mratio);
    b=(param->b);
    bm=(param->bm);
    beta=(param->beta);
    kovergamma=(param->kovergamma);
    memory = (param->memory);
    
    
    u=x;

    inte=integrand(m, n, N, M, eps1, mu1, nu2, eps1ratio, mu1ratio, nu2ratio, b, eps1m, mu1m, nu2m, eps1mratio, mu1mratio, nu2mratio, bm, u, beta, kovergamma, memory);

    return inte.imag();
    
  }


/**************************************************
*** integrand_real_modif: computes the real part of the integrand in t (std precision)
*** instead of u (change of variable u=(1-t)/t )
***************************************************/

  double integrand_real_modif(double t, void *p){
	
    /* function encapsulating the computation of the real part of the integrand of alphamn, 
    for gsl integration, see also function integrand */
    
    struct params *param=(struct params *)p;
    unsigned int N,M; // number of upper and lower layers
    unsigned int m,n; // indices of alphamn (azimuthal mode numbers)
    ap::template_1d_array< amp::campf<PRECISION> > eps1,eps1m,eps1ratio,eps1mratio,mu1,mu1m,mu1ratio,mu1mratio;
    ap::template_1d_array< amp::campf<PRECISION> > nu2,nu2m,nu2ratio,nu2mratio;
    ap::template_1d_array< amp::ampf<PRECISION> > b,bm;
    amp::ampf<PRECISION> beta,gamma,kovergamma,u;
    std::complex<double> inte;
    memorycontainer* memory;

    m=(param->m);
    n=(param->n);
    M=(param->M);
    N=(param->N);
    eps1=(param->eps1);
    eps1m=(param->eps1m);
    eps1ratio=(param->eps1ratio);
    eps1mratio=(param->eps1mratio);
    mu1=(param->mu1);
    mu1m=(param->mu1m);
    mu1ratio=(param->mu1ratio);
    mu1mratio=(param->mu1mratio);
    nu2=(param->nu2);
    nu2m=(param->nu2m);
    nu2ratio=(param->nu2ratio);
    nu2mratio=(param->nu2mratio);
    b=(param->b);
    bm=(param->bm);
    beta=(param->beta);
    kovergamma=(param->kovergamma);
    memory = (param->memory);
    
    
    if (t==0) u=1.e4; // integrand should be very close to zero anyway
    else u=(1.-t)/t;

    inte=integrand(m, n, N, M, eps1, mu1, nu2, eps1ratio, mu1ratio, nu2ratio, b, eps1m, mu1m, nu2m, eps1mratio, mu1mratio, nu2mratio, bm, u, beta, kovergamma, memory);

    if ((t==0)&&(std::abs(inte)>1.e-10)) printf ("Warning: approx for t=0 too rough, %13.8e\n", std::abs(inte));

    return inte.real()/(t*t);
    
  }


/**************************************************
*** integrand_imag_modif: computes the imag part of the integrand in t (std precision)
*** instead of u (change of variable u=(1-t)/t )
***************************************************/

  double integrand_imag_modif(double t, void *p){
	
    /* function encapsulating the computation of the real part of the integrand of alphamn, 
    for gsl integration, see also function integrand */
    
    struct params *param=(struct params *)p;
    unsigned int N,M; // number of upper and lower layers
    unsigned int m,n; // indices of alphamn (azimuthal mode numbers)
    ap::template_1d_array< amp::campf<PRECISION> > eps1,eps1m,eps1ratio,eps1mratio,mu1,mu1m,mu1ratio,mu1mratio;
    ap::template_1d_array< amp::campf<PRECISION> > nu2,nu2m,nu2ratio,nu2mratio;
    ap::template_1d_array< amp::ampf<PRECISION> > b,bm;
    amp::ampf<PRECISION> beta,gamma,kovergamma,u;
    std::complex<double> inte;
    memorycontainer* memory;

    m=(param->m);
    n=(param->n);
    M=(param->M);
    N=(param->N);
    eps1=(param->eps1);
    eps1m=(param->eps1m);
    eps1ratio=(param->eps1ratio);
    eps1mratio=(param->eps1mratio);
    mu1=(param->mu1);
    mu1m=(param->mu1m);
    mu1ratio=(param->mu1ratio);
    mu1mratio=(param->mu1mratio);
    nu2=(param->nu2);
    nu2m=(param->nu2m);
    nu2ratio=(param->nu2ratio);
    nu2mratio=(param->nu2mratio);
    b=(param->b);
    bm=(param->bm);
    beta=(param->beta);
    kovergamma=(param->kovergamma);
    memory = (param->memory);
    
    
    if (t==0) u=1.e4; // integrand should be very close to zero anyway
    else u=(1.-t)/t;

    inte=integrand(m, n, N, M, eps1, mu1, nu2, eps1ratio, mu1ratio, nu2ratio, b, eps1m, mu1m, nu2m, eps1mratio, mu1mratio, nu2mratio, bm, u, beta, kovergamma, memory);
    
    if ((t==0)&&(std::abs(inte)>1.e-10)) printf ("Warning: approx for t=0 too rough, %13.8e\n", std::abs(inte));

    return inte.imag()/(t*t);
    
  }


/**************************************************
*** integrate: performs the integration (using GSL)
*** of integrand_real or integrand_imag
***************************************************/

  double integrate(int flagreal, unsigned int M, unsigned int N, 
  	ap::template_1d_array< amp::ampf<PRECISION> > b, 
	ap::template_1d_array< amp::ampf<PRECISION> > bm, 
	amp::ampf<PRECISION> beta, ap::template_1d_array< amp::campf<PRECISION> > eps1,
	ap::template_1d_array< amp::campf<PRECISION> > eps1m,
	ap::template_1d_array< amp::campf<PRECISION> > mu1,
	ap::template_1d_array< amp::campf<PRECISION> > mu1m,
	amp::ampf<PRECISION> omega, amp::ampf<PRECISION> k, amp::ampf<PRECISION> kovergamma,
	unsigned int m, unsigned int n, size_t limit, gsl_integration_workspace *w,
	memorycontainer* memory){
	
    /* In input: flagreal : 1 to integrate integrand_real, otherwise integrate integrand_imag
       The final m and n are the azimuthal mode numbers (e.g. m=0, n=0 will give alpha_00 )
       The rest are the parameters of the multilayer computation (see etachi and integrand
       functions)*/
    
    struct params param; // input parameters for the integrand functions
    ap::template_1d_array< amp::campf<PRECISION> > eps1ratio,eps1mratio,mu1ratio,mu1mratio,nu2,nu2m,nu2ratio,nu2mratio;
    amp::ampf<PRECISION> beta2,k2;
    gsl_function F;
    double tolint=1.e-6; // relative error permitted for gsl adaptative integration
    double x,err; // result and error
    int status;
   
    if ((m==0)&&(n==0)) tolint*=1e-2; // higher precision for alpha00
 
    // precompute various ratio
    eps1ratio.setbounds(1,N);mu1ratio.setbounds(1,N);nu2.setbounds(1,N+1);nu2ratio.setbounds(1,N);
    eps1mratio.setbounds(1,M);mu1mratio.setbounds(1,M);nu2m.setbounds(1,M+1);nu2mratio.setbounds(1,M);
    beta2=amp::sqr(beta);
    k2=amp::sqr(k);
    //upper layers
    nu2(1)=k2*(1-beta2*eps1(1)*mu1(1));
    set_minus0jimag_to_positive(nu2(1));
    for (unsigned int p=1; p<=N; p++) {
      nu2(p+1)=k2*(1-beta2*eps1(p+1)*mu1(p+1));
      set_minus0jimag_to_positive(nu2(p+1));
      eps1ratio(p)=eps1(p)/eps1(p+1);
      mu1ratio(p)=mu1(p)/mu1(p+1);
      nu2ratio(p)=nu2(p+1)/nu2(p);
      eps1(p)=1/eps1(p);
      mu1(p)=1/mu1(p);
    }
    eps1(N+1)=1/eps1(N+1);
    mu1(N+1)=1/mu1(N+1);
    //lower layers
    nu2m(1)=k2*(1-beta2*eps1m(1)*mu1m(1));
    set_minus0jimag_to_positive(nu2m(1));
    for (unsigned int p=1; p<=M; p++) {
      nu2m(p+1)=k2*(1-beta2*eps1m(p+1)*mu1m(p+1));
      set_minus0jimag_to_positive(nu2m(p+1));
      eps1mratio(p)=eps1m(p)/eps1m(p+1);
      mu1mratio(p)=mu1m(p)/mu1m(p+1);
      nu2mratio(p)=nu2m(p+1)/nu2m(p);
      eps1m(p)=1/eps1m(p);
      mu1m(p)=1/mu1m(p);
    }
    eps1m(M+1)=1/eps1m(M+1);
    mu1m(M+1)=1/mu1m(M+1);
    
        
    // parameters
    param.M=M+1;
    param.N=N+1;
    param.b=b;
    param.bm=bm;
    param.beta=beta;
    param.eps1=eps1;
    param.eps1m=eps1m;
    param.eps1ratio=eps1ratio;
    param.eps1mratio=eps1mratio;
    param.mu1=mu1;
    param.mu1m=mu1m;
    param.mu1ratio=mu1ratio;
    param.mu1mratio=mu1mratio;
    param.nu2=nu2;
    param.nu2m=nu2m;
    param.nu2ratio=nu2ratio;
    param.nu2mratio=nu2mratio;
    param.kovergamma=kovergamma;
    param.m=m;
    param.n=n;
    param.memory = memory;

    // integration with adaptative integration on infinite interval (QAGIU GSL algorithm)
    F.params=&param;

    if (flagreal) {
      // real part
      F.function=&integrand_real;
    } else {
      // imaginary part
      F.function=&integrand_imag;
    }

    gsl_set_error_handler_off();

    status=gsl_integration_qagiu(&F, 0, 0., tolint, 15, w, &x, &err);

    // deal with GSL errors
    if (status) {
      //printf("alpha%d%d, omega= %s, flagreal= %d, result= %13.8e, rel. error= %13.8e\n",m,n,omega.toDec().c_str(),flagreal,x,std::abs(err/x));

      if ( ( (status==GSL_EMAXITER)||(status==GSL_EDIVERGE) )||(status==GSL_EROUND) ) {
	if (flagreal) F.function=&integrand_real_modif;
	else F.function=&integrand_imag_modif;
        status=gsl_integration_qag(&F, 0., 1., 0., tolint, limit, 1, w, &x, &err);
        //printf("GSL_EDIVERGE: alpha%d%d, omega= %s, flagreal= %d, result= %13.8e, rel. error= %13.8e\n",m,n,omega.toDec().c_str(),flagreal,x,std::abs(err/x));
	if (status) printf("Warning: alpha%d%d, omega= %s, flagreal= %d, result= %13.8e, rel. error= %13.8e\n",m,n,omega.toDec().c_str(),flagreal,x,std::abs(err/x));
      }
      /*else if (status==GSL_EROUND) {
	while ((std::abs(err/x)>2*tolint)&&((status==GSL_EROUND)||(status==GSL_EDIVERGE))) {
          tolint*=2;tolintabs*=2;
	  status=gsl_integration_qagiu(&F, 0, tolintabs, tolint, limit, w, &x, &err);
	}
	printf("GSL_EROUND: alpha%d%d, omega= %s, flagreal= %d, result= %13.8e, rel. error= %13.8e\n",m,n,omega.toDec().c_str(),flagreal,x,std::abs(err/x));
      }*/
      else {
        printf("Warning: alpha%d%d, omega= %s, flagreal= %d, result= %13.8e, rel. error= %13.8e\n",m,n,omega.toDec().c_str(),flagreal,x,std::abs(err/x));
      }
    } 
    
    /*if (std::abs(x)<tolintabs) {
      status=gsl_integration_qagiu(&F, 0, std::abs(x), tolint, limit, w, &x, &err);
      if (status) {
        printf("Warning: alpha%d%d, omega= %s, flagreal= %d, result= %13.8e, rel. error= %13.8e\n",m,n,omega.toDec().c_str(),flagreal,x,std::abs(err/x));
      }
    }*/

    return x;
     
  }


/**************************************************
*** alphamn: compute the alpha_mn coefficient
*** for given m and n
***************************************************/

  std::complex<double> alphamn(unsigned int flag_topbotsym, unsigned int M, unsigned int N, 
  	ap::template_1d_array< amp::ampf<PRECISION> > b, 
	ap::template_1d_array< amp::ampf<PRECISION> > bm, 
	amp::ampf<PRECISION> beta, ap::template_1d_array< amp::campf<PRECISION> > eps1,
	ap::template_1d_array< amp::campf<PRECISION> > eps1m,
	ap::template_1d_array< amp::campf<PRECISION> > mu1,
	ap::template_1d_array< amp::campf<PRECISION> > mu1m,
	amp::ampf<PRECISION> omega, amp::ampf<PRECISION> k, amp::ampf<PRECISION> kovergamma,
	unsigned int m, unsigned int n, size_t limit, gsl_integration_workspace *w,
	memorycontainer* memory){
	
    /* In input: flag_topbotsym: flag for top-bottom symmetry (1 if such a symmetry)
       -> this allows to put some terms to zero a priori (when m+n is not even).
       The final m and n are the azimuthal mode numbers (e.g. m=0, n=0 will give alpha_00 )
       The rest are the parameters of the multilayer computation (see etachi and integrand
       functions)*/

    std::complex<double> alpha=std::complex<double>(0.,0.); // alpha_mn (output)
    double x,y; // real and imaginary parts
    
    if ( flag_topbotsym==0 || (m+n)%2==0 ) {
      
      x=integrate(1,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,m,n,limit,w,memory);
      y=integrate(0,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,m,n,limit,w,memory);
      alpha=std::complex<double>(x,y);
    }

    return alpha;
    
  }

/***********************************************************************
*** function to compute the factorial of a positive integer		     ***
***********************************************************************/

  long factorial(int n){
    long fact = 1;

    while (n>1) {
      fact*=n;
      n--;
    }

    return fact;
  }


/*****************************************************
*** function to compute (-1)^n					   ***
*****************************************************/

  int minus1pow(int n){
    int res;
 
    if (n%2==0) res=1;
    else res= -1;
 
    return res;
  }

/**************************************************
*** coefmn: compute the coefficient in front of 
*** alpha_mn, for a given non-linear term of 
*** the impedance.
***************************************************/

  double coefmn(int a, int b, int c, int d, int m, int n, int p, int q){
	
    /* In input: the exponents involved in the non-linear term x1^a*y1^b*x^c*y^d,
       the m and n of the alpha_mn involved, and the integer indices 
       p, q of the sums*/

    double coef=0.;
    int aplusc2,dm0=1,dn0=1;

    aplusc2 = (a+c)/2;
    if (m==0) dm0 = 2;
    if (n==0) dn0 = 2;

    for (int k=0; k<=aplusc2; k++) {

      if ( (aplusc2 <= (q+k)) && ((q+k) <= aplusc2+d/2) ) {
        coef += double(minus1pow(c+k)*factorial(a+c)*factorial(n))/ \
          double(dm0*dn0*factorial(2*k)*factorial(n-2*k)*factorial(p)*factorial(b-p)* \
           factorial(aplusc2-k)*factorial(q+k-aplusc2)*factorial(n+q)*factorial(a)*factorial(c));
      }
    }

    return coef;
  }
