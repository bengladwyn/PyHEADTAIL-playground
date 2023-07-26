#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <string.h>

#include <complex>
#include <cmath>
#include <stdlib.h>

#include <IW2D.h>

using std::complex;
using std::log;
using std::exp;
using std::abs;
using std::max;
using std::min;
using std::cout;

/******************************************************************************
 *** locate: search a table of doubles ordered in ascending order	    ***
 *** Effect         : Function that gives the position lprov (integer) in   ***
 ***                  in table, such that table[lprov-1]<z<table[lprov]	    ***
 *** Parameters     : table, z, n (table is indexed from 0 to n)            ***
 ******************************************************************************/

unsigned long locate (double *table, double z, unsigned long n)

{ unsigned long il, iu, im, lprov;

 il=0;
 iu=n+1;
 while (iu-il >1) {
   im=(iu+il)/2; // midpoint
   if (z >= table[im])
     il=im;
   else
     iu=im;
   }
 if (z==table[0]) lprov=0;
 else if (z==table[n]) lprov=n; 
 else if (z<table[0]) lprov=0; 
 else if (z>table[n]) lprov=n+1; 
 else lprov=iu;

 return lprov;

}

/******************************************************************************
 *** locateMP: search a multiprecision table ordered in ascending order	    ***
 *** (same as above but multiprecision version)                             ***
 *** Effect         : Function that gives the position lprov (integer) in   ***
 ***                  in table, such that table[lprov-1]<z<table[lprov]	    ***
 *** Parameters     : table, z, n (table is indexed from 0 to n)            ***
 ******************************************************************************/

unsigned long locateMP (ap::template_1d_array< amp::ampf<PRECISION> >& table, 
	amp::ampf<PRECISION> z, unsigned long n)

{ unsigned long il, iu, im, lprov;

 il=0;
 iu=n+1;
 while (iu-il >1) {
   im=(iu+il)/2; // midpoint
   if (z >= table(im))
     il=im;
   else
     iu=im;
   }
 if (z==table(0)) lprov=0;
 else if (z==table(n)) lprov=n; 
 else if (z<table(0)) lprov=0; 
 else if (z>table(n)) lprov=n+1; 
 else lprov=iu;

 return lprov;

}

/**************************************************
*** pchip: interpolated value at z, given the function 
*** values fi at the points zi and the derivatives of the 
*** interpolated polynomial di at zi. Tables have size nz.
***************************************************/

  complex<double> pchip(double z, double *zi, 
  	complex<double> *fi, complex<double> *di, unsigned long nz) {
	
    double delta,t,t2,t3,tm,tm2,tm3,h1,h2,h3,h4;
    complex<double> value;
    unsigned long l;
    
    l=locate(zi, z, nz-1);
    if (z==zi[0]) l=1; /* change value from 0 to 1 in this case (because of the slightly modifed 
    		version of locate we use here) */
		
    delta=zi[l]-zi[l-1];
    //cout << l << " " << delta << "\n";
    t=(zi[l]-z)/delta;
    tm=(z-zi[l-1])/delta;
    t2=t*t;t3=t*t2;
    tm2=tm*tm;tm3=tm*tm2;
    
    // cubic Hermite basis functions
    h1=3.*t2-2.*t3;
    h2=3.*tm2-2.*tm3;
    h3=-delta*(t3-t2);
    h4=delta*(tm3-tm2);
    
    // interpolated value
    value=fi[l-1]*h1+fi[l]*h2+di[l-1]*h3+di[l]*h4;
	
    return value;

  }


/**************************************************
*** interp: interpolated value at z, given the function 
*** values fi at the points zi and the derivatives of the 
*** interpolated polynomial di at zi (in case of pchip interpolation).
*** Tables have size nz.
*** interp_type is the type of interpolation: 
*** 0 for pchip, 1 for linear, 2 for exponential, 3 for power law
***************************************************/

  complex<double> interp(double z, double *zi, complex<double> *fi, 
  	complex<double> *di, unsigned long nz, unsigned int interp_type) {
	
    complex<double> value;
    unsigned long l;
    
    switch (interp_type) {
      case 0:
        value=pchip(z, zi, fi, di, nz);
	return value;
    
      case 1:
        l=locate(zi, z, nz-1);
        if (z==zi[0]) l=1; /* change value from 0 to 1 in this case (because of the slightly modifed 
    		version of locate we use here) */
        // linearly interpolated value
        value=fi[l-1]+(fi[l]-fi[l-1])*(z-zi[l-1])/(zi[l]-zi[l-1]);
	return value;
	
      /*case 2:
        l=locate(zi, z, nz-1);
        if (z==zi[0]) l=1; // change value from 0 to 1 in this case (because of the slightly modifed 
    			   // version of locate we use here)
        // log of the function is linearly interpolated
        value=exp(log(fi[l-1])+(log(fi[l])-log(fi[l-1]))*(z-zi[l-1])/(zi[l]-zi[l-1]));
	return value;

      case 3:
        l=locate(zi, z, nz-1);
        if (z==zi[0]) l=1; // change value from 0 to 1 in this case (because of the slightly modifed 
    			   //version of locate we use here)
        // log of the function w.r.t to log of the variable is linearly interpolated
        value=exp(log(fi[l-1])+(log(fi[l])-log(fi[l-1]))*(log(z)-log(zi[l-1]))/(log(zi[l])-log(zi[l-1])));
	return value;
      */
    }
    
  }


/**************************************************
*** pchip_derivatives: derivatives of the pchip interpolating
*** polynomial (Sources: Monotone Piecewise Cubic Interpolation, 
*** F. N. Fritsch & R. E. Carlson, SIAM Journal on Numerical Analysis,
*** Vol. 17, No. 2, 1980, pp. 238-246, and F. N. Fritsch & J. Butland,
*** A method for constructing local monotone piecewise cubic interpolants,
*** SIAM J. Sci. Comput., 5(2), 1984, pp. 300-304)
***************************************************/

  void pchip_derivatives(complex<double>* d, double* x, complex<double> * y, 
	  unsigned long length){
    /* Derivative values for monotonic piecewise cubic Hermite interpolation
      d = pchip_derivatives(x,y,length) gives first derivatives,
      d[k] = P'(x[k]), of interpolating polynomial P */

    // derivatives from finite differences of y vs. x
    complex<double>* delta=new complex<double>[length-1];
    for (unsigned long i=0; i<length-1; i++) {
      delta[i] = (y[i+1] - y[i])/(x[i+1]-x[i]);
    }

    //  particular case n=2 -> linear interp.
    if (length==2) {
      //complex<double>* d = new complex<double>[2];
      d[0]=delta[0];
      d[1]=delta[0];
      return;
    }

    //  Derivatives for 0 < k < length-1:
    //      d[k] = weighted average of delta[k-1] and delta[k]
    //                 when they are of same sign.
    //      d[k] = 0 otherwise (different signs or one of them is zero).

    //complex<double>* d = new complex<double>[length];
    for (unsigned long i=0; i<length; i++) {
      d[i]=0.;
    }

    std::vector<unsigned long> k;
    for (unsigned long i=0; i<length-2; i++) {
      if (abs(delta[i])!=0. || abs(delta[i+1])!=0.) {
	k.push_back(i);
      }
    }

    double h[length-1];
    for (unsigned long i=0; i<length-1; i++){
      h[i]=x[i+1]-x[i];
    }

    double h_sum[k.size()];
    double w1[k.size()];
    double w2[k.size()];
    for (unsigned long i=0; i<k.size(); i++){
      h_sum[i] = h[k[i]] + h[k[i]+1];
      w1[i] = ( h[ k[i] ] + h_sum[i]) / (3.*h_sum[i]);
      w2[i] = ( h[ k[i]+1 ] + h_sum[i]) / (3.*h_sum[i]);
    }

    double dmax[k.size()];
    double dmin[k.size()];

    for (unsigned long i=0; i<k.size(); i++){
      dmax[i] = max ( abs(delta[k[i]]), abs(delta[k[i]+1] )); 
      dmin[i] = min ( abs(delta[k[i]]), abs(delta[k[i]+1] )); 
      d[k[i]+1] = dmin[i] / std::conj( (w1[i]*delta[k[i]] + w2[i]*delta[k[i]+1])/dmax[i] );
      /*cout << "i= " << i << " ; dmax= " << dmax[i] << " ; dmin= " << dmin[i] << '\n' ;
      cout << "k[i]= " << k[i] << " ; w1= " << w1[i] << " ; w2= " << w2[i] << '\n' ;
      cout << "delta[k[i]]= " << delta[k[i]] << " ; delta[k[i]+1]= " << delta[k[i]+1] << " ; d[k[i]+1]= " << d[k[i]+1] << '\n' ;
      cout << '\n' ;*/
    }

    //  Derivatives at the edges (one-sided, 3-point formula - see e.g. C. Moler,
    // 2004, DOI:10.1137/1.9780898717952).

    d[0] = (((double)2.*h[0]+h[1])*delta[0] - h[0]*delta[1])/(h[0]+h[1]);

    //if isreal(d) && (sign(d[1]) ~= sign(delta[1]))
    //d[1] = 0;
    if (delta[0]/abs(delta[0]) != delta[1]/abs(delta[1]) && (abs(d[0]) > (double)3.*abs(delta[0]) )) {
      d[0] = (double)3.*delta[0];
    }

    d[length-1] = (((double)2.*h[length-2]+h[length-3])*delta[length-2] - h[length-2]*delta[length-3])/(h[length-2]+h[length-3]);

    if (delta[length-2]/abs(delta[length-2]) != delta[length-3]/abs(delta[length-3]) && abs(d[length-1]) > abs((double)3.*delta[length-2])) {
      d[length-1] = (double)3.*delta[length-2];
    }
    
    delete[] delta;
    return;
  }
