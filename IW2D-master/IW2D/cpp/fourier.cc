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


/***********************************************************
*** Function to compute Phi(x)
***********************************************************/

  complex<long double> Phi(complex<long double> x,double eps){
    
    complex<long double> Phi(0.L,0.L),power(1.L,0.L);
    int condition;
    int i=0;
    long double fact=1.L,fact2,abspower,il,expoabs;
    
    if (abs(x)<1.L) {
      expoabs=exp(abs(x));
      do {  
            // fact is 1/factorial(i) and power is x^i
	    il=(long double)i;
	    fact2 = (il+6.L)/((il+3.L)*(il+4.L));
            Phi+= power * fact * fact2;
	    power *= x;
	    abspower = abs(power);
	    fact = fact/(il+1.L);
	    //cout << abs(std::real(Phi)) << " " << abs(std::imag(Phi)) <<
	    //		" " << min(abs(std::real(Phi)),abs(std::imag(Phi))) << '\n';
      	    condition =  (2.L* abspower * expoabs * fact / (il+5.L) > 
	    	(long double)eps*min(abs(std::real(Phi)),abs(std::imag(Phi))));
      	    i++;} while ( condition ) ;
    }
    else {
      complex<long double> expo = exp(x);
      Phi= (std::pow(x,3)*expo - 6.L*x*(expo+1.L) + 12.L*(expo-1.L))/(std::pow(x,4));
    }
    return Phi;
  }
  

/***********************************************************
*** Function to compute Psi(x)
***********************************************************/

  complex<long double> Psi(complex<long double> x,double eps){
    
    complex<long double> Psi(0.L,0.L),power(1.L,0.L);
    int condition;
    long double fact=1.L,fact2,abspower,il,expoabs;
    int i=0;
     
    if (abs(x)<1.L) {
      expoabs=exp(abs(x));
      do {  
            // fact is 1/factorial(i) and power is x^i
	    il=(long double)i;
	    fact2 = 1.L/((il+3.L)*(il+4.L));
	    Psi+= - power * fact * fact2;
	    power *= x;
	    fact = fact/(il+1.L);
	    abspower = abs(power);
      	    condition =  (abspower * expoabs * fact / ((il+4.L)*(il+5.L)) >
	    	(long double)eps*min(abs(std::real(Psi)),abs(std::imag(Psi))));
      	    i++;} while ( condition ) ;
    }
    else {
      complex<long double> expo = exp(x);
      Psi = (-std::pow(x,2)*expo + 2.L*x*(2.L*expo+1.L) - 6.L*(expo-1.L))/(std::pow(x,4));
    }
    return Psi;
  }
  
/***********************************************************
*** Function to compute Lambda(x)
***********************************************************/

  complex<long double> Lambda(complex<long double> x,double eps){
    
    complex<long double> Lambda(0.L,0.L),power(1.L,0.L);
    int condition;
    long double fact=1.L,fact2,abspower,il,expoabs;
    int i=0;
     
    if (abs(x)<1.L) {
      expoabs=exp(abs(x));
      do {  
            // fact is 1/factorial(i) and power is x^i
	    il=(long double)i;
	    fact2 = 1.L/(il+2.L);
	    Lambda+= power * fact * fact2;
	    power *= x;
	    abspower = abs(power);
      	    fact = fact/(il+1.L);
	    condition =  (abspower * expoabs * fact/(il+3.L) >
	    	(long double)eps*min(abs(std::real(Lambda)),abs(std::imag(Lambda))));
      	    i++;} while ( condition ) ;
    }
    else {
      complex<long double> expo = exp(x);
      Lambda = (x*expo - expo+1.L)/(x*x);
    }
    return Lambda;
  }

/***********************************************************
*** Function to compute Fourier integral on semi-infinite domain
***********************************************************/

complex<double> fourier_integral_inf(complex<double>* fi,complex<double>* d,complex<double>
	df, long double t, long double* omegai, long double* delta, unsigned long length, double eps,
	unsigned int* interp_type, int flaginf) {

// Computes the Fourier integral at t from omegai(1) to omegai(end) of f(omega)*exp(j*omega*t)
// where f is a function. Add a correcting term to take into account the
// rest of the integral between omegai(end) and infinity (or -infinity and
// -omegai(1) if omegai(end)<0)

// In input (omegai,fi) are the abscissae where f is interpolated to perform the 
// calculation, and the corresponding values of f. df is the derivative of f at
// the freqi which is maximum in absolute value.
// We use a cubic interpolation between the points omegai (monotonic piecewise cubic Hermite interpolation - pchip). 
// The derivatives of the interpolation
// are in d (of size length). delta (of size length-1) are the differences
// between sucessive freqi. eps is the relative precision in the computation of the auxiliary
// functions Phi and Psi.
// The type of interpolation on each interval is in interp_type (0 for pchip, 1 for linear)
 
// We use Filon's type method, plus an asymptotic method for the correcting term toward infinity if flaginf=1.

  long double s; 
  complex<long double> x,expo,fint(0.L,0.L),Lambdat,Lambdamt,Phit,Phimt,Psit,Psimt,dummy1,dummy2,aprime,bprime;

  // computes the fourier integral at t between omegai(1) and omegai(end)
  for (unsigned long i=0; i<length-1; i++) {
      s=delta[i]*t;
      x=complex<long double>(0.L,s);
      expo=exp(x);
      switch (interp_type[i]){
        case 0:
          // pchip interpolation
	  // computes some auxiliary variables
	  Phit=Phi(x,eps);
	  Phimt=Phi(-x,eps);
	  Psit=Psi(x,eps);
	  Psimt=Psi(-x,eps);      
	  //complex<double> dummy1 = exp(complex<double>(0.,omegai[i+1]*t) ) * (fi[i]*Phimt - d[i] * delta[i] * Psimt);
	  //complex<double> dummy2 = exp(complex<double>(0.,omegai[i]*t) ) * (fi[i+1]*Phit + d[i+1] * delta[i] * Psit);
	  //fint += delta[i]*(dummy1+dummy2);
	  dummy1 = delta[i]* (-expo * (complex<long double>)d[i]*Psimt +
	  	(complex<long double>)d[i+1]* Psit);
	  dummy2 = (complex<long double>)fi[i+1]*Phit + expo * (complex<long double>)fi[i]*Phimt;
	  fint += exp(complex<long double>(0.L,omegai[i]*t) )*delta[i]*(dummy1+dummy2);
	  break;
	  
        case 1:
          // linear interpolation
	  Lambdat=Lambda(x,eps);
	  Lambdamt=Lambda(-x,eps);
	  dummy1 = (complex<long double>)fi[i+1]*Lambdat + expo * (complex<long double>)fi[i]*Lambdamt;
	  fint += exp(complex<long double>(0.L,omegai[i]*t) )*delta[i]*dummy1;
	  break;
	  
        /*case 2:
          // exponential interpolation
	  break;
	*/
      }
    }
  //cout << std::real(d[length-1]) << "  " << std::imag(d[length-1]) << '\n'; // to check interpolation 
  //fi(1:n-1).*Phimt - d(1:n-1).*delta.*Psimt)

  // add the correcting term for the rest of the integral
  if (flaginf==1) {
    if (omegai[length-1]>=0) {
      fint+= exp(complex<long double>(0.L,t*omegai[length-1]))*
    	  (-(complex<long double>)fi[length-1]/complex<long double>(0.L,t) - 
	  (complex<long double>)df*(1.L/(t*t)));
      //cout << t << " " << exp(complex<double>(0.,t*omegai[length-1]))*(-fi[length-1]/complex<double>(0.,t) - df*(1./(t*t))) << "\n";
    }
    else {
      fint+=exp(complex<long double>(0.L,t*omegai[0]))*
    	  ((complex<long double>)fi[0]/complex<long double>(0.L,t) + 
	  (complex<long double>)df*(1.L/(t*t)));
    }
  }
  return (complex<double>)fint;
}
