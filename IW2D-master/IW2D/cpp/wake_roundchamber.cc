/*
 *  wake_roundchamber.cc 
 *  
 *  by Nicolas Mounet (Nicolas.Mounet@cern.ch)
 *
 *  computes the wakes in a flat chamber (see CERN note by N. Mounet and E. Metral, 
 "Electromagnetic fields and beam coupling impedances in a multilayer flat chamber", 2010)
 
In input : typical input file is

Machine:	LHC
Relativistic Gamma:	479.6
Impedance Length in m:	1
Number of layers:	3
Layer 1 inner radius in mm:	4
Layer 1 DC resistivity (Ohm.m):	4.3e-7
Layer 1 relaxation time for resistivity (ps):	0
Layer 1 real part of dielectric constant:	1
Layer 1 magnetic susceptibility:	0
Layer 1 relaxation frequency of permeability (MHz):	Infinity
Layer 1 thickness in mm:	0.003
Layer 2 DC resistivity (Ohm.m):	4e+12
Layer 2 relaxation time for resistivity (ps):	0
Layer 2 real part of dielectric constant:	4
Layer 2 magnetic susceptibility:	0
Layer 2 relaxation frequency of permeability (MHz):	Infinity
Layer 2 thickness in mm:	54
Layer 3 DC resistivity (Ohm.m):	7.2e-7
Layer 3 relaxation time for resistivity (ps):	0
Layer 3 real part of dielectric constant:	1
Layer 3 magnetic susceptibility:	0
Layer 3 relaxation frequency of permeability (MHz):	Infinity
Layer 3 thickness in mm:	Infinity
linear (1) or logarithmic (0) or both (2) scan in z for the wake:	0
sampling distance in m for the linear sampling:	0.5e-5
zmin in m of the linear sampling:	0.5e-5
zmax in m of the linear sampling:	0.01
Number of points per decade for the logarithmic sampling:	100
exponent (10^) of zmin (in m) of the logarithmic sampling:	-1.99
exponent (10^) of zmax (in m) of the logarithmic sampling:	6
added z [m]:	5e6 1e7
Yokoya factors long, xdip, ydip, xquad, yquad:	1 1 1 0 0
factor weighting the longitudinal impedance error:	100.
tolerance (in wake units) to achieve:	1.e9
frequency above which the mesh bisecting is linear [Hz]:	1.e11
Comments for the output files names:	_some_element

The order of the lines can be whatever, but the exact sentences and the TAB before the parameter
indicated, are necessary. If top-bottom symmetry is set (with "yes" or "y" or "1") the lower layers (with a
minus sign) are ignored. Also if there are more layers than indicated by the number of upper (lower) layers,
the additional one are ignored. The last layer is always assumed to go to infinity.

In output one gives five files with the wakes (longitudinal, x dipolar, y dipolar,
x quadrupolar, y quadrupolar and y constant term). Each has 2 columns : distance (behind source) and
wake (SI units).
In output one also gives six files with the impedances computed on the final mesh used to 
compute the wakes (files are for longitudinal, x dipolar, y dipolar,
x quadrupolar, y quadrupolar and y constant, impedances). Each file has 3 columns : frequency, real part
and imaginary part of the impedance (SI units).

 */

#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include <complex>
#include <cmath>
#include <stdlib.h>
#include <ablas.h>
#include <amp.h>
#include <mpfr.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sort_double.h>

#include "IW2D.h"
#include <globals.h>

/*
#define  Precision	160  // Number of bits for the mantissa of all numbers (classic doubles have 53 bits)
#define  MAXLINES	100  // Maximum number of lines in the input file (correspond to 5 layers in top and
			    // bottom parts)
#define  MAXCHAR	200  // Maximum number of characters in one line in the input file
#define  MAXCHARFILE	200  // Maximum number of characters in the output files name extension
*/

using std::complex;
using std::log;
using std::exp;
using std::abs;
using std::max;
using std::min;
using std::cout;


struct params_diff {unsigned int N;
        ap::template_1d_array< amp::ampf<PRECISION> > rho;
        ap::template_1d_array< amp::ampf<PRECISION> > tau;
        ap::template_1d_array< amp::ampf<PRECISION> > epsb;
        ap::template_1d_array< amp::ampf<PRECISION> > chi;
        ap::template_1d_array< amp::ampf<PRECISION> > fmu;
        ap::template_1d_array< amp::ampf<PRECISION> > b;
        amp::ampf<PRECISION> beta; amp::ampf<PRECISION> gamma;
	double L; double *freqi; unsigned int interp_type; unsigned long nf;
	complex<double> *Zxdipfi; complex<double> *Zxdipdi;
	complex<double> *Zxquadfi; complex<double> *Zxquaddi;
        complex<double> *Zlongfi; complex<double> *Zlongdi;
	};


// global arrays with memory of impedances
double *freqmem;
complex<double> *Zxdipmem,*Zlongmem,*Zxquadmem;
unsigned long impmem; // current number of elements of freqmem, Zxdipmem, etc.

double factlong=100.; // multiplication factor for longitudinal impedance (precision should be better for long. wake) (default value)

const unsigned long MAXMEM=150000;  // Maximum number of elements of arrays with impedances

/**************************************************
*** impedance: computes the impedances (std precision)
***
***************************************************/

  void impedance(complex<double>& Zxdip, complex<double>& Zxquad,
	complex<double>& Zlong, unsigned int N,
	ap::template_1d_array< amp::ampf<PRECISION> >& rho, 
  	ap::template_1d_array< amp::ampf<PRECISION> >& tau, 
  	ap::template_1d_array< amp::ampf<PRECISION> >& epsb, 
  	ap::template_1d_array< amp::ampf<PRECISION> >& chi, 
  	ap::template_1d_array< amp::ampf<PRECISION> >& fmu, 
	ap::template_1d_array< amp::ampf<PRECISION> >& b,
	amp::ampf<PRECISION> beta, amp::ampf<PRECISION> gamma, 
	double L, double freq) {
	
    /* function that computes the impedances at a given frequency freq, from rho, tau, epsb, chi, fmu, b 
    of each of the N layers, and from the relativistic velocity factor beta. */

    amp::ampf<PRECISION> omega, k;
    ap::template_1d_array< amp::campf<PRECISION> > eps1,mu1; /* eps1 and mu1 for the layers */
    amp::campf<PRECISION> jimagMP; // imaginary constant in multiprecision
    complex<double> cst,alphaTM0,alphaTM1,jimag; // some constant, alphamn constants and imaginary constant
    unsigned long lprov;
    
    // try to find freq in the table freqmem
    //cout.precision(18);
    //cout << "impedance: Frequency " << freq << "\n";
    //cout.flush();
    if (impmem==0) lprov=0;
    else lprov=locate(freqmem,freq,impmem-1);
 
    if ( (impmem!=0)&&(freq==freqmem[lprov]) ) {
      Zxdip=Zxdipmem[lprov];
      Zxquad=Zxquadmem[lprov];Zlong=Zlongmem[lprov];
      //cout << "Frequency found in impedance memory " << lprov << "\n";
    }
    else if ( ( (impmem!=0)&&(lprov>0)) && (freq==freqmem[lprov-1]) ) {
      Zxdip=Zxdipmem[lprov-1];
      Zxquad=Zxquadmem[lprov-1];Zlong=Zlongmem[lprov-1];
      //cout << "Frequency found in impedance memory " << lprov-1 << "\n";
    }
    
    
    else {

      //cout << "Frequency not found in impedance memory\n";
      //cout.flush();

      //cout << "lprov=" << lprov << "\n";
      //cout.flush();

      // constants
      jimagMP.x=0;jimagMP.y=1;
      jimag=complex<double>(0.,1.);

      omega=amp::twopi<PRECISION>()*freq;
      k=omega/(beta*C);
   
      // allocation for eps1, mu1
      eps1.setbounds(1,N+1);mu1.setbounds(1,N+1);

      // first layer (inside the chamber) is always vacuum
      eps1(1)=1;mu1(1)=1;
      
      // computes the layer properties for the angular freq. omega
      for (unsigned int p=2;p<=N+1; p++) {
        if (rho(p).isFiniteNumber()) {
	  eps1(p)=epsb(p)+1/(jimagMP*eps0*rho(p)*omega*(1+jimagMP*omega*tau(p)));
	} else {
	  eps1(p)=epsb(p);
	}
	mu1(p)=1+chi(p)/(1+jimagMP*omega/(amp::twopi<PRECISION>()*fmu(p)));
	//cout << p << "\n";
	//printf("%s %s %s %s\n%s %s\n",rho(p).toDec().c_str(),b(p).toDec().c_str(),eps1(p).x.toDec().c_str(),eps1(p).y.toDec().c_str(),
	//		mu1(p).x.toDec().c_str(),mu1(p).y.toDec().c_str());
      }

      //printf("%s \n",b(1).toDec().c_str());
      //printf("%s %s\n",omega.toDec().c_str(),k.toDec().c_str());
      //for (unsigned int i=0;i<500;i++) cout << "tata ";cout << "\n";cout.flush();

      // computes alphaTM0
      alphaTM0=alphaTM(N+1,0,eps1,mu1,b,k,beta);
      //printf("alphaTM0: %13.8e %13.8e\n",alphaTM0.real(),alphaTM0.imag());

      // computes alphaTM1
      alphaTM1=alphaTM(N+1,1,eps1,mu1,b,k,beta);
      //printf("alphaTM1: %13.8e %13.8e\n",alphaTM1.real(),alphaTM1.imag());


      // computes the impedances (without Yokoya factors, only axisymmetric)
      cst=jimag*L*double(amp::ampf<PRECISION>(k*Z0/(beta*amp::sqr(gamma)*amp::twopi<PRECISION>())).toDouble());
      Zlong=cst*alphaTM0;
      cst=cst*double(amp::ampf<PRECISION>(k/(2*amp::sqr(gamma))).toDouble());
      Zxdip=cst*alphaTM1;
      Zxquad=cst*alphaTM0;

      //cout.precision(18);
      //cout << freq << " " << Zlong << " " << Zxdip << " " << Zxquad << '\n';
      //cout << freq << " " << Zxdip.real() << " " << Zxdip.imag() << '\n';
      //cout.flush();

      // add to the memory
      for(unsigned int k=impmem; k>=lprov+1; k--) {
	freqmem[k]=freqmem[k-1];
	Zxdipmem[k]=Zxdipmem[k-1];
	Zxquadmem[k]=Zxquadmem[k-1];
	Zlongmem[k]=Zlongmem[k-1];
      }
      freqmem[lprov]=freq;
      Zxdipmem[lprov]=Zxdip;
      Zxquadmem[lprov]=Zxquad;
      Zlongmem[lprov]=Zlong;
      impmem++;
       
    }

    return;

  }


/***********************************************************
*** integrand_diff: function computing
*** the difference in norm between
*** the impedances at freq=exp(u) and their interpolation
***********************************************************/

  double integrand_diff(double u, void *p){

    /* function encapsulating the computation of the difference in norm between the 
    impedances and their interpolation */

    struct params_diff *param=(struct params_diff *)p;
    unsigned int N,interp_type; // number of layers, and interpolation type
    unsigned long nf; // number of frequencies in interpolation 
    // multilayer parameters
    ap::template_1d_array< amp::ampf<PRECISION> > rho,tau,epsb,chi,fmu,b;
    amp::ampf<PRECISION> beta,gamma; // relativistic velocity factor
    double *freqi; // frequencies of the interpolation
    complex<double> *Zxdipfi,*Zxdipdi; // values and derivatives of the interpolation at the frequencies freqi (Zxdip)
    complex<double> *Zxquadfi,*Zxquaddi; // values and derivatives of the interpolation at the frequencies freqi (Zyquad)
    complex<double> *Zlongfi,*Zlongdi; // values and derivatives of the interpolation at the frequencies freqi (Zlong)
    complex<double> Zxdip,Zxquad,Zlong; // exact impedances
    complex<double> pZxdip,pZxquad,pZlong; // interpolated impedances
    double L,result,freq;
    int flag_topbotsym;
    
    freq=exp(u);

    N=(param->N);
    rho=(param->rho);
    tau=(param->tau);
    epsb=(param->epsb);
    chi=(param->chi);
    fmu=(param->fmu);
    b=(param->b);
    beta=(param->beta);
    gamma=(param->gamma);
    L=(param->L);

    freqi=(param->freqi);
    interp_type=(param->interp_type);
    nf=(param->nf);
    Zxdipfi=(param->Zxdipfi);Zxdipdi=(param->Zxdipdi);
    Zxquadfi=(param->Zxquadfi);Zxquaddi=(param->Zxquaddi);
    Zlongfi=(param->Zlongfi);Zlongdi=(param->Zlongdi);
    

    // computes the impedances at freq
    impedance(Zxdip, Zxquad, Zlong, N, rho, tau, epsb, chi, fmu, b, 
	beta, gamma, L, freq);

    pZxdip=interp(freq,freqi,Zxdipfi,Zxdipdi,nf,interp_type);
    pZxquad=interp(freq,freqi,Zxquadfi,Zxquaddi,nf,interp_type);
    pZlong=interp(freq,freqi,Zlongfi,Zlongdi,nf,interp_type);

    /*result=max(abs(Zxdip-pZxdip),abs(Zydip-pZydip));
    result=max(result,max(abs(Zxquad-pZxquad),abs(Zyquad-pZyquad)));
    result=max(result,max(factlong*abs(Zlong-pZlong),abs(Zycst-pZycst)));*/
    result=abs(Zxdip-pZxdip)*abs(Zxdip-pZxdip);
    result+=factlong*factlong*abs(Zxquad-pZxquad)*abs(Zxquad-pZxquad);
    result+=factlong*factlong*abs(Zlong-pZlong)*abs(Zlong-pZlong);
    /*result=abs(Zxdip.real()-pZxdip.real())+abs(Zxdip.imag()-pZxdip.imag());
    result=max(result,abs(Zydip.real()-pZydip.real())+abs(Zydip.imag()-pZydip.imag()));
    result=max(result,abs(Zyquad.real()-pZyquad.real())+abs(Zyquad.imag()-pZyquad.imag()));
    result=max(result,abs(Zycst.real()-pZycst.real())+abs(Zycst.imag()-pZycst.imag()));
    result=max(result,factlong*(abs(Zlong.real()-pZlong.real())+abs(Zlong.imag()-pZlong.imag())));*/

    //cout.precision(18);
    //cout << u << "\t" << std::sqrt(result)*freq << '\n';
    //cout << freq << "\t" << Zlong.real() << " " << Zlong.imag() << " " << Zxdip.real() << " "<< Zxdip.imag() << " " << Zydip.real() << " "<< Zydip.imag() << " " << Zyquad.real()<< " " << Zyquad.imag() << " " << Zycst.real() << " " << Zycst.imag() << " " << '\n';
    //cout << freq << "\t" << pZlong.real() << " " << pZlong.imag() << " " << pZxdip.real() << " "<< pZxdip.imag() << " " << pZydip.real() << " "<< pZydip.imag() << " " << pZyquad.real()<< " " << pZyquad.imag() << " " << pZycst.real() << " " << pZycst.imag() << " " << '\n';
    //cout.flush();
    
    return std::sqrt(result)*freq;
    
  }

/***********************************************************
*** integrand_diff_freq: function computing
*** the difference in norm between
*** the impedances at freq and their interpolation.
*** variable is not u=log(freq), it is here freq directly
***********************************************************/

  double integrand_diff_freq(double freq, void *p){

    /* function encapsulating the computation of the difference in norm between the 
    impedances and their interpolation */
    double u;
    
    u=log(freq);
    return integrand_diff(u,p)/freq;
    
  }


/**************************************************
 *** 			main program		***
 ***						***
 ***************************************************/

main ()

{
 
 char *endline;
 char output[MAXCHARFILE],Zxdipoutput[MAXCHARFILE+10],Zydipoutput[MAXCHARFILE+10],Zlongoutput[MAXCHARFILE+10];
 char Zxquadoutput[MAXCHARFILE+10],Zyquadoutput[MAXCHARFILE+10],Input[MAXCHARFILE+10];
 char Zxdipoutput2[MAXCHARFILE+10],Zydipoutput2[MAXCHARFILE+10],Zlongoutput2[MAXCHARFILE+10];
 char Zxquadoutput2[MAXCHARFILE+10],Zyquadoutput2[MAXCHARFILE+10];
 char Wxdipoutput[MAXCHARFILE+10],Wydipoutput[MAXCHARFILE+10],Wlongoutput[MAXCHARFILE+10];
 char Wxquadoutput[MAXCHARFILE+10],Wyquadoutput[MAXCHARFILE+10];
 char Wxdipoutput2[MAXCHARFILE+10],Wydipoutput2[MAXCHARFILE+10],Wlongoutput2[MAXCHARFILE+10];
 char Wxquadoutput2[MAXCHARFILE+10],Wyquadoutput2[MAXCHARFILE+10];
 std::string data[MAXLINES],machine,commentoutput,dummy2;
 FILE *filZxdip, *filZydip, *filZxquad, *filZyquad, *filZlong, *filInput;
 FILE *filZxdip2, *filZydip2, *filZxquad2, *filZyquad2, *filZlong2;
 FILE *filWxdip, *filWydip, *filWxquad, *filWyquad,*filWlong;
 FILE *filWxdip2, *filWydip2, *filWxquad2, *filWyquad2, *filWlong2;
 unsigned int N,dummy0; // number of layers, then dummy parameter
 unsigned int n_input,n_added; /* number of input lines, number of individually added frequencies */
 double dzlin; // dz in linear scan
 unsigned int typescan,nzlog; /* flag for top-bottom symmetry (1 if such a symmetry), type of frequency scan,
 		number of z per decade in log. scan */
 unsigned long kmain,ind[20],lprov,lprov2,j,imax=ULONG_MAX,nz,nzdup,nf; // some temporary indices, total number of z in the scan, number of duplicate z, number of frequencies;
 ap::template_1d_array< amp::ampf<PRECISION> > b,thick; // position of boundaries, and thickness of the layers
 ap::template_1d_array< amp::ampf<PRECISION> > rho,tau,epsb,chi,fmu; /* layers
 						properties (DC resistivity, resistivity relaxation time,
						dielectric constant, magnetic susceptibility=mu_r-1,
						relaxation frequency of permeability)*/
 amp::ampf<PRECISION> beta,gamma,dummy3; // parameters
 // The next three are for gsl adaptative integration algorithm
 gsl_integration_workspace *w;
 gsl_function F;
 double tolintabs=1.e-2,tolintrel=1.e-2; // absolute and relative error permitted for gsl adaptative integration
 double freqlin=1.e11; // frequency limit above which we switch from a log mesh to a linear mesh (default value)
 double eps=1.e-20; // precision for the calculation of Phi, Psi and Lambda
 double x,y,xlin,maxi,err,L,zminlog,zmaxlog,zminlin,zmaxlin,zadded[15],dummy1,yokoya[5];
 double *freq,*z,*t,dif,sum,*newfreq,*inte;
 long double *omegai,*delta;
 double freqmin0=1.e-5,freqmin,freqmax,freqmax0;
 complex<double> *Zxdipfi,*Zxdipdi,*Zxquadfi,*Zxquaddi,*Zlongfi,*Zlongdi; //impedances
 complex<double> *newZxdipdi,*newZxquaddi,*newZlongdi; //impedances derivatives on subintervals
 complex<double> *newZxdipfi,*newZxquadfi,*newZlongfi; //impedances on subintervals
 double *Wakexdip,*Wakexquad,*Wakelong; //wakes
 double *Wakexdipold,*Wakexquadold,*Wakelongold; //wakes (previous loop)
 complex<double> jimag; // imaginary constant
 unsigned int *interp_type; // kind of interpolation on each interval (0 for pchip, 1 for linear)
 struct params_diff param; // input parameters for the integrand functions (gsl integration)
 size_t limit=1000,found,neval; // limit (number of intervals) for gsl integration algorithm
 time_t start,end,time1,time2; // times
 clock_t c1,c2; // clock ticks
 bool condition_int,condition_freqmin,condition_freqmax;
 //std::vector<double> newfreq;
 static const long double pi = 3.141592653589793238462643383279502884197;
 double tol=1.e9; /* absolute error permitted on the integral of the difference (in norm)
 	betweens impedances and their interpolation (default value) */
 
 
					
 // start time
 time(&start);
 
 // allocation for memory of impedances
 freqmem=new double[MAXMEM];
 Zxdipmem=new complex<double>[MAXMEM];
 Zxquadmem=new complex<double>[MAXMEM];Zlongmem=new complex<double>[MAXMEM];
 impmem=0;

 // allocation of freq, impedances and derivatives for the full interpolation
 freq=new double[MAXMEM];interp_type=new unsigned int[MAXMEM];inte=new double[MAXMEM];
 Zxdipfi=new complex<double>[MAXMEM];
 Zxquadfi=new complex<double>[MAXMEM];Zlongfi=new complex<double>[MAXMEM];
 Zxdipdi=new complex<double>[MAXMEM];
 Zxquaddi=new complex<double>[MAXMEM];Zlongdi=new complex<double>[MAXMEM];
 // allocation of freq, impedances and derivatives for interpolation on a subinterval
 newfreq=new double[5];
 newZxdipfi=new complex<double>[5];
 newZxquadfi=new complex<double>[5];newZlongfi=new complex<double>[5];
 newZxdipdi=new complex<double>[5];
 newZxquaddi=new complex<double>[5];newZlongdi=new complex<double>[5];

 // default values of the parameters (in case)
 typescan=0;
 zminlog=-2;zmaxlog=6;nzlog=10;n_added=0;
 zminlin=1.e-2;zmaxlin=1;dzlin=1.e-2;nz=0;
 N=2;L=1.;gamma="479.6";
 
 // read input file
 // first read everything to identify the strings in front of each parameters
 n_input=0;
 while (std::cin.eof()==0) {
   std::getline (std::cin,data[n_input+1]);
   /* next line is when the input file comes from windows or else and has some ^M characters in the
   end of each line */
   //data[n_input+1]=data[n_input+1].substr(0,data[n_input+1].length()-1);
   //cout << data[n_input+1] << '\n';
   n_input++;
 }
 n_input--;
 //printf("n_input: %d\n",n_input);
 // identify each argument
 for (unsigned int i=1; i<=n_input; i++) {
   read_input(data[i],"Machine",dummy0,dummy1,machine,dummy3,2);
   read_input(data[i],"Relativistic Gamma",dummy0,dummy1,dummy2,gamma,3);
   read_input(data[i],"Impedance Length in m",dummy0,L,dummy2,dummy3,1);
   read_input(data[i],"Number of layers",N,dummy1,dummy2,dummy3,0);
   //printf("yoyo %d %s %s %s %13.8e %d\n",i,data[i].c_str(),machine.c_str(),gamma.toDec().c_str(),L,N);
   read_input(data[i],"exponent (10^) of zmin (in m) of the logarithmic sampling",dummy0,zminlog,dummy2,dummy3,1);
   read_input(data[i],"exponent (10^) of zmax (in m) of the logarithmic sampling",dummy0,zmaxlog,dummy2,dummy3,1);
   read_input(data[i],"linear (1) or logarithmic (0) or both (2) scan in z for the wake",typescan,dummy1,dummy2,dummy3,0);
   read_input(data[i],"Number of points per decade for the logarithmic sampling",nzlog,dummy1,dummy2,dummy3,0);
   read_input(data[i],"zmin in m of the linear sampling",dummy0,zminlin,dummy2,dummy3,1);
   read_input(data[i],"zmax in m of the linear sampling",dummy0,zmaxlin,dummy2,dummy3,1);
   read_input(data[i],"sampling distance in m for the linear sampling",dummy0,dzlin,dummy2,dummy3,1);
   read_input(data[i],"factor weighting the longitudinal impedance error",dummy0,factlong,dummy2,dummy3,1);
   read_input(data[i],"tolerance (in wake units) to achieve",dummy0,tol,dummy2,dummy3,1);
   read_input(data[i],"frequency above which the mesh bisecting is linear",dummy0,freqlin,dummy2,dummy3,1);
   read_input(data[i],"Comments for the output files names",dummy0,dummy1,commentoutput,dummy3,2);
   read_input_double_list(data[i], "added z", zadded, n_added);
   read_input_yokoya_factors(data[i], yokoya);
 }
 
 tol/=(2.*(double)pi); // from angular frequancy to frequency units
 
 //printf("%s %d \n",machine,N);
 //printf("%13.8e %13.8e %d\n",zadded[1],zadded[2],n_added);
 //printf("%13.8e %13.8e %d %ld\n",zminlog,zmaxlog,typescan,nzlog);
 //printf("%13.8e %13.8e %13.8e\n",freqlin,tol,factlong);

 b.setbounds(1,N+1);thick.setbounds(1,N+1);rho.setbounds(1,N+1);tau.setbounds(1,N+1);
 epsb.setbounds(1,N+1);chi.setbounds(1,N+1);fmu.setbounds(1,N+1);

 // default values of the layers properties (in case)
 rho(2)="1e5";tau(2)=0;epsb(2)=1;chi(2)=0;fmu(2)="Infinity";b(1)=2;b(2)="Infinity";

 // find inner half gap(s) of the chamber
 for (unsigned int i=1; i<=n_input; i++) {
   read_input(data[i],"Layer 1 inner radius in mm",dummy0,dummy1,dummy2,b(1),3);
 }
 //printf("%s %s \n",b(1).toDec().c_str(),gamma.toDec().c_str());
 
 // find all the layers properties
 for (unsigned int i=1; i<=n_input; i++) {
   for (unsigned int p=1; p<=N; p++){
     read_input_layer(data[i],"DC resistivity",p,rho(p+1));
     read_input_layer(data[i],"relaxation time for resistivity (ps)",p,tau(p+1));
     read_input_layer(data[i],"real part of dielectric constant",p,epsb(p+1));
     read_input_layer(data[i],"magnetic susceptibility",p,chi(p+1));
     read_input_layer(data[i],"relaxation frequency of permeability",p,fmu(p+1));
     read_input_layer(data[i],"thickness in mm",p,thick(p));
   }
 }

 // units conversion to SI (tau was in ps, b in mm and fmu in MHz)
 b(1)=b(1)*1e-3;
 for (unsigned int p=2; p<=N+1; p++){
   tau(p)=tau(p)*1e-12;
   b(p)=thick(p-1)*1e-3+b(p-1);
   fmu(p)=fmu(p)*1e6;
   //printf("%s %s %s %s %s %s \n",rho(p).toDec().c_str(),tau(p).toDec().c_str(),epsb(p).toDec().c_str(),
   //	chi(p).toDec().c_str(),fmu(p).toDec().c_str(),b(p).toDec().c_str());
 }
  
 // relativistic velocity factor beta
 beta=amp::sqrt(1-1/amp::sqr(gamma));
 
 // construct the z scan
 // first estimation of the number of z (for memory allocation)
 switch(typescan) {
   case 0:
     nz=(int)ceil((zmaxlog-zminlog)*(double)nzlog)+1+n_added;
     break;
   case 1:
     nz=(int)ceil((zmaxlin-zminlin)/dzlin)+1+n_added;
     break;
   case 2:
     nz=(int)ceil((zmaxlog-zminlog)*(double)nzlog)+1+(int)ceil((zmaxlin-zminlin)/dzlin)+1+n_added;
     break;
   }
 z=new double[nz];
 
 // constructs unsorted version of the array z
 nz=0;z[0]=-1.;

 if (typescan==1) {
   do {
     z[nz]=zminlin+(double)nz * dzlin;
     nz++;
     } while(z[nz-1]<zmaxlin);
   }
 else {
   do {
     z[nz]=pow(10.,zminlog+(double)nz/(double)nzlog);
     nz++;
     } while(z[nz-1]<pow(10.,zmaxlog));
   if (typescan==2) {
     j=0;
     do {
       z[nz]=zminlin+(double)(j) * dzlin;
       nz++;j++;
       } while(z[nz-1]<zmaxlin);
     }
   }
 for (unsigned int i=0;i<n_added; i++) {
   z[nz]=zadded[i];
   nz++;
   }
 // nz is now the real number of z
 
 // sort the distances z
 gsl_sort(z, 1, nz);
 
 // remove duplicate distances
 nzdup=0;
 for (long i=nz-2;i>=0; i--) {
   if ( (z[i+1]-z[i]) ==0 ) {
       z[i]=DBL_MAX;
       nzdup++;
     }
   }
 gsl_sort(z, 1, nz);
 nz=nz-nzdup;
 // computes times
 t=new double[nz];
 for (unsigned long i=0;i<nz;i++) t[i]=z[i]/double(amp::ampf<PRECISION>(beta*C).toDouble());
 //printf("%d\n",nz);
 //for (unsigned long i=0;i<nz;i++) printf("%d %13.8e\n",i,z[i]);
 Wakexdip=new double[nz];
 Wakexquad=new double[nz];
 Wakelong=new double[nz];
 Wakexdipold=new double[nz];
 Wakexquadold=new double[nz];
 Wakelongold=new double[nz];
 
 
 // set the output files names
 //sprintf(output,"W%s_%dlayersup_%dlayersdown%.2lfmm%s",machine,N,M,
 //	1e3*double(amp::ampf<Precision>(b(1)).toDouble()),commentoutput);
 sprintf(output,"W%s_%dlayers%.2lfmm%s",machine.c_str(),N,
	1e3*double(amp::ampf<PRECISION>(b(1)).toDouble()),commentoutput.c_str());
 sprintf(Zxdipoutput,"Zxdip%s.dat",output);
 sprintf(Zydipoutput,"Zydip%s.dat",output);
 sprintf(Zxquadoutput,"Zxquad%s.dat",output);
 sprintf(Zyquadoutput,"Zyquad%s.dat",output);
 sprintf(Zlongoutput,"Zlong%s.dat",output);
 sprintf(Zxdipoutput2,"Zxdip%s_precise.dat",output);
 sprintf(Zydipoutput2,"Zydip%s_precise.dat",output);
 sprintf(Zxquadoutput2,"Zxquad%s_precise.dat",output);
 sprintf(Zyquadoutput2,"Zyquad%s_precise.dat",output);
 sprintf(Zlongoutput2,"Zlong%s_precise.dat",output);
 sprintf(Input,"InputData%s.dat",output);
 sprintf(Wxdipoutput,"Wxdip%s.dat",output);
 sprintf(Wydipoutput,"Wydip%s.dat",output);
 sprintf(Wxquadoutput,"Wxquad%s.dat",output);
 sprintf(Wyquadoutput,"Wyquad%s.dat",output);
 sprintf(Wlongoutput,"Wlong%s.dat",output);
 sprintf(Wxdipoutput2,"Wxdip%s_precise.dat",output);
 sprintf(Wydipoutput2,"Wydip%s_precise.dat",output);
 sprintf(Wxquadoutput2,"Wxquad%s_precise.dat",output);
 sprintf(Wyquadoutput2,"Wyquad%s_precise.dat",output);
 sprintf(Wlongoutput2,"Wlong%s_precise.dat",output);
 
 // writes the InputData file (copy of input file)
 std::ofstream InputFile (Input);
 for (unsigned int i=1; i<=n_input; i++) {
   InputFile << data[i] << '\n';
 }
 InputFile.close();

 // open the output files and write the first line (column description)
 filZxdip=fopen(Zxdipoutput,"w");
 filZydip=fopen(Zydipoutput,"w");
 filZxquad=fopen(Zxquadoutput,"w");
 filZyquad=fopen(Zyquadoutput,"w");
 filZlong=fopen(Zlongoutput,"w");
 fprintf(filZxdip,"Frequency [Hz]\tRe(Zxdip) [Ohm/m]\tIm(Zxdip) [Ohm/m]\n");
 fprintf(filZydip,"Frequency [Hz]\tRe(Zydip) [Ohm/m]\tIm(Zydip) [Ohm/m]\n");
 fprintf(filZxquad,"Frequency [Hz]\tRe(Zxquad) [Ohm/m]\tIm(Zxquad) [Ohm/m]\n");
 fprintf(filZyquad,"Frequency [Hz]\tRe(Zyquad) [Ohm/m]\tIm(Zyquad) [Ohm/m]\n");
 fprintf(filZlong,"Frequency [Hz]\tRe(Zlong) [Ohm]\tIm(Zlong) [Ohm]\n");
 filZxdip2=fopen(Zxdipoutput2,"w");
 filZydip2=fopen(Zydipoutput2,"w");
 filZxquad2=fopen(Zxquadoutput2,"w");
 filZyquad2=fopen(Zyquadoutput2,"w");
 filZlong2=fopen(Zlongoutput2,"w");
 fprintf(filZxdip2,"Frequency [Hz]\tRe(Zxdip) [Ohm/m]\tIm(Zxdip) [Ohm/m]\n");
 fprintf(filZydip2,"Frequency [Hz]\tRe(Zydip) [Ohm/m]\tIm(Zydip) [Ohm/m]\n");
 fprintf(filZxquad2,"Frequency [Hz]\tRe(Zxquad) [Ohm/m]\tIm(Zxquad) [Ohm/m]\n");
 fprintf(filZyquad2,"Frequency [Hz]\tRe(Zyquad) [Ohm/m]\tIm(Zyquad) [Ohm/m]\n");
 fprintf(filZlong2,"Frequency [Hz]\tRe(Zlong) [Ohm]\tIm(Zlong) [Ohm]\n");
 filWxdip=fopen(Wxdipoutput,"w");
 filWydip=fopen(Wydipoutput,"w");
 filWxquad=fopen(Wxquadoutput,"w");
 filWyquad=fopen(Wyquadoutput,"w");
 filWlong=fopen(Wlongoutput,"w");
 fprintf(filWxdip,"Distance [m]\tWake x dip [V/(C.m)]\n");
 fprintf(filWydip,"Distance [m]\tWake y dip [V/(C.m)]\n");
 fprintf(filWxquad,"Distance [m]\tWake x quad [V/(C.m)]\n");
 fprintf(filWyquad,"Distance [m]\tWake y quad [V/(C.m)]\n");
 fprintf(filWlong,"Distance [m]\tWake long [V/C]\n");
 filWxdip2=fopen(Wxdipoutput2,"w");
 filWydip2=fopen(Wydipoutput2,"w");
 filWxquad2=fopen(Wxquadoutput2,"w");
 filWyquad2=fopen(Wyquadoutput2,"w");
 filWlong2=fopen(Wlongoutput2,"w");
 fprintf(filWxdip2,"Distance [m]\tWake x dip [V/(C.m)]\n");
 fprintf(filWydip2,"Distance [m]\tWake y dip [V/(C.m)]\n");
 fprintf(filWxquad2,"Distance [m]\tWake x quad [V/(C.m)]\n");
 fprintf(filWyquad2,"Distance [m]\tWake y quad [V/(C.m)]\n");
 fprintf(filWlong2,"Distance [m]\tWake long [V/C]\n");
 
 
 // workspace allocation for gsl adaptative integration
 w=gsl_integration_workspace_alloc(limit);
 // deactivate gsl errors
 gsl_set_error_handler_off();
 
 // some parameters initialization
 param.N=N;
 param.rho=rho;
 param.tau=tau;
 param.epsb=epsb;
 param.chi=chi;
 param.fmu=fmu;
 param.b=b;
 param.beta=beta;
 param.gamma=gamma;
 param.L=L;
 
 // first guess for freqmin and freqmax (minimum and maximum frequencies of the interpolation)
 freqmin=freqmin0;
 freqmax0=10.*double(amp::ampf<PRECISION>(beta*gamma*C/(b(1)*amp::twopi<PRECISION>())).toDouble());
 freqmax=freqmax0;
	
 condition_freqmin=true;
 condition_freqmax=true;
 srand(time(NULL));
 
 nf=2; // begins with two frequencies
 freq[0]=freqmin;freq[1]=freqmax;
 condition_int=true;
 
 kmain=0;

 // loop to get accurate enough mesh of the impedances
 while ( condition_int ) {

   printf("Number of frequencies: %ld, iteration nb %ld\n",nf,kmain);
   cout.flush();

   if (kmain==0) {
     for (unsigned long i=0; i<=nf-1; i++) {
       //time(&time1);
       impedance(Zxdipfi[i], Zxquadfi[i], Zlongfi[i], N, rho, tau, epsb, chi, fmu, b,
		   beta, gamma, L, freq[i]);
       /*time(&time2);
       dif=difftime(time2,time1);
       printf("Elapsed time during calculation: %.5lf seconds\n",dif);
       cout.flush();*/
     }
     pchip_derivatives(Zxdipdi,freq,Zxdipfi,nf);
     pchip_derivatives(Zxquaddi,freq,Zxquadfi,nf);
     pchip_derivatives(Zlongdi,freq,Zlongfi,nf);

     // writing the rest of the input parameters structure
     param.freqi=freq;param.nf=nf;
     param.interp_type=0; // pchip interpolation chosen for this step
     param.Zxdipfi=Zxdipfi;param.Zxdipdi=Zxdipdi;
     param.Zxquadfi=Zxquadfi;param.Zxquaddi=Zxquaddi;
     param.Zlongfi=Zlongfi;param.Zlongdi=Zlongdi;

     // parameters for gsl integration
     F.params=&param;

     condition_int=false;
     //newfreq.clear();

     sum=0.;
      
     // make an adaptative integration to get a first trial mesh
     F.function=&integrand_diff; // function to integrate with gsl
     //F.function=&integrand_diff_freq; // function to integrate with gsl
     //time(&time1);
     gsl_integration_qag(&F, log(freq[0]), log(freq[1]), tolintabs, tolintrel, limit,1, w, &sum, &err);
     //gsl_integration_qag(&F, freq[0], freq[1], tolintabs, tolintrel, limit,1, w, &sum, &err);
     // replace freq by all the frequencies in the impedance memory
     //for (unsigned long i=0; i<=impmem-1; i++) newfreq.push_back(freqmem[i]);
     for (unsigned long i=0; i<=impmem-1; i++) {
       if (i<impmem-1) interp_type[i]=0; // pchip is chosen
       freq[i]=freqmem[i];Zxdipfi[i]=Zxdipmem[i];
       Zxquadfi[i]=Zxquadmem[i];Zlongfi[i]=Zlongmem[i];
     }
     nf=impmem;
     condition_int=true;
   }
   else if (kmain==1) {
     
     // compute the pchip derivatives
     pchip_derivatives(Zxdipdi,freq,Zxdipfi,nf);
     pchip_derivatives(Zxquaddi,freq,Zxquadfi,nf);
     pchip_derivatives(Zlongdi,freq,Zlongfi,nf);

     // rewrite the input parameters structure
     param.freqi=freq;param.nf=nf;
     param.Zxdipfi=Zxdipfi;param.Zxdipdi=Zxdipdi;
     param.Zxquadfi=Zxquadfi;param.Zxquaddi=Zxquaddi;
     param.Zlongfi=Zlongfi;param.Zlongdi=Zlongdi;
     
     err=DBL_MIN;sum=0.;
     for (unsigned long i=0; i<nf-1; i++) {
       if (freq[i]>=freqlin) y=(freq[i]+freq[i+1])/2.;	 
       else y=(log(freq[i])+log(freq[i+1]))/2.;
       // first try pchip (monotone piecewise cubic) interpolation
       param.interp_type=0;F.params=&param;interp_type[i]=0;
       if (freq[i]>=freqlin) inte[i]=integrand_diff_freq(freq[i],F.params)+4.*integrand_diff_freq(y,F.params)+integrand_diff_freq(freq[i+1],F.params);
       else inte[i]=integrand_diff(log(freq[i]),F.params)+4.*integrand_diff(y,F.params)+integrand_diff(log(freq[i+1]),F.params);
       // then try linear
       param.interp_type=1;F.params=&param;
       if (freq[i]>=freqlin) xlin=integrand_diff_freq(freq[i],F.params)+4.*integrand_diff_freq(y,F.params)+integrand_diff_freq(freq[i+1],F.params);
       else xlin=integrand_diff(log(freq[i]),F.params)+4.*integrand_diff(y,F.params)+integrand_diff(log(freq[i+1]),F.params);
       if (xlin<inte[i]) {
         inte[i]=xlin;interp_type[i]=1;
       }
       if (freq[i]>=freqlin) inte[i]*=(freq[i+1]-freq[i])/6.;
       else inte[i]*=(log(freq[i+1])-log(freq[i]))/6.;
       sum+=inte[i];
       if (inte[i]>err) {
	 err=inte[i];imax=i; // interval with largest error
       }
     }
     cout << "err=" << err << ", imax=" << imax << ", f[imax]=" << freq[imax] << ", f[imax+1]=" << freq[imax+1] << ", sum=" << sum << "\n";
     condition_int=(sum>tol);
     if ( (!condition_int)&&!(sum<=tol) ) condition_int=true; // case when sum=nan -> will still try another loop
     if (condition_int) nf++;
     /*if (condition_int) {
       newfreq.push_back(freq[0]);
       for (unsigned long i=0; i<nf-1; i++) {
         if (i==imax) {
           // refine the worst interval found above with the points in the impedance memory
	   lprov=locate(freqmem,freq[i],impmem-1);lprov2=locate(freqmem,freq[i+1],impmem-1);
	   if ( (lprov>0)&&(freq[i]==freqmem[lprov-1]) ) lprov=lprov-1;
	   if ( (lprov2>0)&&(freq[i+1]==freqmem[lprov2-1]) ) lprov2=lprov2-1;
           for (unsigned long k=lprov+1; k<lprov2; k++) {
	     //cout << freqmem[k] << "\n";
	     newfreq.push_back(freqmem[k]);
	   }
	 }
	 newfreq.push_back(freq[i+1]);
       }
     }*/
   } else {
   
     // first subtract the interval imax from the integral (in 'sum') previously computed
     sum-=inte[imax];

     // bisect the previously found subinterval imax (the one with maximum error)

     // construct the five points array to compute the new pchip derivatives on the 3 points freq[imax],
     // freq[imax+1] and the new frequency (half way in log between the two)
     for (unsigned long i=0; i<5; i++) {
       if (i!=2) {
         newfreq[i]=freq[imax-1+i-(i>2)];newZxdipfi[i]=Zxdipfi[imax-1+i-(i>2)];
	 newZxquadfi[i]=Zxquadfi[imax-1+i-(i>2)];
	 newZlongfi[i]=Zlongfi[imax-1+i-(i>2)];
       } else {
         if (freq[imax]>=freqlin) newfreq[i]=(freq[imax]+freq[imax+1])/2.;
         else newfreq[i]=exp((log(freq[imax])+log(freq[imax+1]))/2.);
	 impedance(newZxdipfi[i], newZxquadfi[i], newZlongfi[i], N, rho, tau, epsb, chi, fmu, b,
		   beta, gamma, L, newfreq[i]);
       }
     }
     // pchip derivatives on those five points
     pchip_derivatives(newZxdipdi,newfreq,newZxdipfi,5);
     pchip_derivatives(newZxquaddi,newfreq,newZxquadfi,5);
     pchip_derivatives(newZlongdi,newfreq,newZlongfi,5);
     // construct the new mesh with impedances and derivatives
     for (unsigned long i=nf-1; i>imax+1; i--) {
       if (i<(nf-1)) {
         interp_type[i]=interp_type[i-1];inte[i]=inte[i-1];
       }
       freq[i]=freq[i-1]; 
       Zxdipfi[i]=Zxdipfi[i-1];
       Zxquadfi[i]=Zxquadfi[i-1]; Zlongfi[i]=Zlongfi[i-1];
       Zxdipdi[i]=Zxdipdi[i-1];
       Zxquaddi[i]=Zxquaddi[i-1]; Zlongdi[i]=Zlongdi[i-1];
     }
     freq[imax+1]=newfreq[2];
     Zxdipfi[imax+1]=newZxdipfi[2];
     Zxquadfi[imax+1]=newZxquadfi[2]; Zlongfi[imax+1]=newZlongfi[2];
     interp_type[imax]=0; interp_type[imax+1]=0;// choose by default pchip interpolation on the two subintervals
     for (unsigned long i=imax; i<=imax+2; i++) {
       Zxdipdi[i]=newZxdipdi[1+i-imax];
       Zxquaddi[i]=newZxquaddi[1+i-imax]; Zlongdi[i]=newZlongdi[1+i-imax];
     }
     
     // rewrite the input parameters structure
     param.freqi=freq;param.nf=nf;
     param.Zxdipfi=Zxdipfi;param.Zxdipdi=Zxdipdi;
     param.Zxquadfi=Zxquadfi;param.Zxquaddi=Zxquaddi;
     param.Zlongfi=Zlongfi;param.Zlongdi=Zlongdi;
     
     // compute the integral (with Simpson's rule) on the two subintervals created
     for (unsigned long i=imax;i<=imax+1;i++) {
       if (freq[i]>=freqlin) y=(freq[i]+freq[i+1])/2.;	 
       else y=(log(freq[i])+log(freq[i+1]))/2.;
       // first try pchip
       param.interp_type=0;F.params=&param;interp_type[i]=0;
       if (freq[i]>=freqlin) inte[i]=integrand_diff_freq(freq[i],F.params)+4.*integrand_diff_freq(y,F.params)+integrand_diff_freq(freq[i+1],F.params);
       else inte[i]=integrand_diff(log(freq[i]),F.params)+4.*integrand_diff(y,F.params)+integrand_diff(log(freq[i+1]),F.params);
       // then try linear
       param.interp_type=1;F.params=&param;
       if (freq[i]>=freqlin) xlin=integrand_diff_freq(freq[i],F.params)+4.*integrand_diff_freq(y,F.params)+integrand_diff_freq(freq[i+1],F.params);
       else xlin=integrand_diff(log(freq[i]),F.params)+4.*integrand_diff(y,F.params)+integrand_diff(log(freq[i+1]),F.params);
       if (xlin<inte[i]) {
         inte[i]=xlin;interp_type[i]=1;
       }
       if (freq[i]>=freqlin) inte[i]*=(freq[i+1]-freq[i])/6.;
       else inte[i]*=(log(freq[i+1])-log(freq[i]))/6.;
       sum+=inte[i];
     }
     // find largest error 'err' and its index 'imax'
     err=DBL_MIN;
     for (unsigned long i=0;i<nf-1;i++) {
       if (inte[i]>err) {
         err=inte[i];imax=i;
       }
     }
     cout << "err=" << err << ", imax=" << imax << ", f[imax]=" << freq[imax] << ", f[imax+1]=" << freq[imax+1] << ", sum=" << sum << "\n";

     condition_int=(sum>tol); // condition to continue
     if ( (!condition_int)&&!(sum<=tol) ) condition_int=true; // case when sum=nan -> will still try another loop
     if (condition_int) nf++;
   
   }
   
   
   kmain++;

 }
 
 /*cout << "Final mesh and interpolation type: " << "\n";
 for (unsigned long i=0;i<nf-1;i++) {
   cout << freq[i] << " " << interp_type[i] << "\n";
 }*/

 // for the derivatives, convert frequencies to angular frequencies
 for (unsigned long i=0; i<=nf-1; i++) {
   Zxdipdi[i]=Zxdipdi[i]/(2.*(double)pi);
   Zxquaddi[i]=Zxquaddi[i]/(2.*(double)pi);Zlongdi[i]=Zlongdi[i]/(2.*(double)pi);
 }
 
 // loop to get low enough minimum frequency
 while (condition_freqmin) {
   
   printf("freq. min= %13.8e\n",freqmin);cout.flush();
   
   // computes omega and delta between successive omegas
   omegai=new long double[nf];delta=new long double[nf-1];
   for (unsigned long i=0; i<=nf-1; i++) omegai[i]=2.L*pi*(long double)freq[i];
   for (unsigned long i=0; i<nf-1; i++) delta[i]=omegai[i+1]-omegai[i];
   
   printf("Wake computation\n");cout.flush();
   // computes the wakes
   if (freqmin==freqmin0) {
     for (unsigned long i=0; i<=nz-1; i++) {
       Wakexdip[i]=std::imag(fourier_integral_inf(Zxdipfi,Zxdipdi,0.,(long double)t[i],omegai,delta,nf,eps,interp_type,1))/pi;
       Wakexquad[i]=std::imag(fourier_integral_inf(Zxquadfi,Zxquaddi,0.,(long double)t[i],omegai,delta,nf,eps,interp_type,1))/pi;
       Wakelong[i]=std::real(fourier_integral_inf(Zlongfi,Zlongdi,0.,(long double)t[i],omegai,delta,nf,eps,interp_type,1))/pi;
     }
   } else {
     for (unsigned long i=0; i<=nz-1; i++) {
       Wakexdip[i]=Wakexdipold[i]+std::imag(fourier_integral_inf(newZxdipfi,newZxdipdi,0.,(long double)t[i],omegai,delta,
     		  3,eps,interp_type,0))/pi - std::imag(fourier_integral_inf(&Zxdipfi[1],&Zxdipdi[1],0.,
		  (long double)t[i],&omegai[1],&delta[1],2,eps,&interp_type[1],0))/pi;
       Wakexquad[i]=Wakexquadold[i]+std::imag(fourier_integral_inf(newZxquadfi,newZxquaddi,0.,(long double)t[i],omegai,delta,
     		  3,eps,interp_type,0))/pi - std::imag(fourier_integral_inf(&Zxquadfi[1],&Zxquaddi[1],0.,
		  (long double)t[i],&omegai[1],&delta[1],2,eps,&interp_type[1],0))/pi;
       Wakelong[i]=Wakelongold[i]+std::real(fourier_integral_inf(newZlongfi,newZlongdi,0.,(long double)t[i],omegai,delta,
     		  3,eps,interp_type,0))/pi - std::real(fourier_integral_inf(&Zlongfi[1],&Zlongdi[1],0.,
		  (long double)t[i],&omegai[1],&delta[1],2,eps,&interp_type[1],0))/pi;
     }
     for (unsigned int i=0;i<=1;i++) {
       Zxdipdi[i]=newZxdipdi[i];
       Zxquaddi[i]=newZxquaddi[i];
       Zlongdi[i]=newZlongdi[i];
     }

   }
     
   
   if (freqmin!=freqmin0) {
     err=DBL_MIN;
     for (unsigned long i=0; i<=nz-1; i++) {
       err=max(abs(Wakexdip[i]-Wakexdipold[i]),err);
       err=max(abs(Wakexquad[i]-Wakexquadold[i]),err);
       err=max(factlong*abs(Wakelong[i]-Wakelongold[i]),err);
     }
     cout << "Max. error between two last minimum frequencies chosen : " << err << "\n";cout.flush();
     condition_freqmin=(err>=tol);
   }
   
   if ( condition_freqmin ) {
     freqmin/=10.;nf++;
     for (unsigned long i=nf-1; i>0; i--) {
       if (i<(nf-1)) interp_type[i]=interp_type[i-1];
       freq[i]=freq[i-1];
       Zxdipfi[i]=Zxdipfi[i-1];
       Zxquadfi[i]=Zxquadfi[i-1];Zlongfi[i]=Zlongfi[i-1];
       Zxdipdi[i]=Zxdipdi[i-1];
       Zxquaddi[i]=Zxquaddi[i-1];Zlongdi[i]=Zlongdi[i-1];
     }
     freq[0]=freqmin;interp_type[0]=0; // use pchip only
     impedance(Zxdipfi[0], Zxquadfi[0], Zlongfi[0], N, rho, tau, epsb, chi, fmu, b,
		   beta, gamma, L, freqmin);
     for (unsigned int i=0;i<=2;i++) {
       newfreq[i]=freq[i];newZxdipfi[i]=Zxdipfi[i];
       newZxquadfi[i]=Zxquadfi[i];
       newZlongfi[i]=Zlongfi[i];
     }
     pchip_derivatives(newZxdipdi,newfreq,newZxdipfi,3);newZxdipdi[2]=Zxdipdi[2]; // last derivative remains unchanged
     pchip_derivatives(newZxquaddi,newfreq,newZxquadfi,3);newZxquaddi[2]=Zxquaddi[2]; // last derivative remains unchanged
     pchip_derivatives(newZlongdi,newfreq,newZlongfi,3);newZlongdi[2]=Zlongdi[2]; // last derivative remains unchanged
     // convert to angular frequencies
     for (unsigned int i=0;i<=1;i++) {
       newZxdipdi[i]=newZxdipdi[i]/(2.*(double)pi);
       newZxquaddi[i]=newZxquaddi[i]/(2.*(double)pi);
       newZlongdi[i]=newZlongdi[i]/(2.*(double)pi);
     }     
     delete[] omegai;delete[] delta;
   }
   
   for (unsigned long i=0; i<=nz-1; i++) {
     Wakexdipold[i]=Wakexdip[i];
     Wakexquadold[i]=Wakexquad[i];
     Wakelongold[i]=Wakelong[i];
   }

 }

 // loop to get high enough maximum frequency
 while (condition_freqmax) {
   
   freqmax*=2.;nf++;
   printf("freq. max= %13.8e\n",freqmax);cout.flush();
   freq[nf-1]=freqmax;interp_type[nf-2]=0; // use pchip only
   impedance(Zxdipfi[nf-1], Zxquadfi[nf-1], Zlongfi[nf-1], N, rho, tau, epsb, chi, fmu, b,
		 beta, gamma, L, freqmax);
   for (unsigned long i=0;i<=2;i++) {
     newfreq[i]=freq[nf-3+i];newZxdipfi[i]=Zxdipfi[nf-3+i];
     newZxquadfi[i]=Zxquadfi[nf-3+i];
     newZlongfi[i]=Zlongfi[nf-3+i];
   }
   pchip_derivatives(newZxdipdi,newfreq,newZxdipfi,3);newZxdipdi[0]=Zxdipdi[nf-3]; // first derivative remains unchanged
   pchip_derivatives(newZxquaddi,newfreq,newZxquadfi,3);newZxquaddi[0]=Zxquaddi[nf-3]; // first derivative remains unchanged
   pchip_derivatives(newZlongdi,newfreq,newZlongfi,3);newZlongdi[0]=Zlongdi[nf-3]; // first derivative remains unchanged
   // convert to angular frequencies
   for (unsigned int i=1;i<=2;i++) {
     newZxdipdi[i]=newZxdipdi[i]/(2.*(double)pi);
     newZxquaddi[i]=newZxquaddi[i]/(2.*(double)pi);
     newZlongdi[i]=newZlongdi[i]/(2.*(double)pi);
   }     
   
   // computes omega and delta between successive omegas
   delete[] omegai;delete[] delta;
   omegai=new long double[3];delta=new long double[2];
   for (unsigned long i=0; i<=2; i++) omegai[i]=2.L*pi*(long double)freq[nf-3+i];
   for (unsigned long i=0; i<=1; i++) delta[i]=omegai[i+1]-omegai[i];
   
   printf("Wake computation\n");cout.flush();
   // computes the wakes
   for (unsigned long i=0; i<=nz-1; i++) {
     Wakexdip[i]=Wakexdipold[i]+std::imag(fourier_integral_inf(newZxdipfi,newZxdipdi,0.,(long double)t[i],omegai,delta,
     		3,eps,&interp_type[nf-3],1))/pi - std::imag(fourier_integral_inf(&Zxdipfi[nf-3],&Zxdipdi[nf-3],0.,
		(long double)t[i],omegai,delta,2,eps,&interp_type[nf-3],1))/pi;
     Wakexquad[i]=Wakexquadold[i]+std::imag(fourier_integral_inf(newZxquadfi,newZxquaddi,0.,(long double)t[i],omegai,delta,
     		3,eps,&interp_type[nf-3],1))/pi - std::imag(fourier_integral_inf(&Zxquadfi[nf-3],&Zxquaddi[nf-3],0.,
		(long double)t[i],omegai,delta,2,eps,&interp_type[nf-3],1))/pi;
     Wakelong[i]=Wakelongold[i]+std::real(fourier_integral_inf(newZlongfi,newZlongdi,0.,(long double)t[i],omegai,delta,
     		3,eps,&interp_type[nf-3],1))/pi - std::real(fourier_integral_inf(&Zlongfi[nf-3],&Zlongdi[nf-3],0.,
		(long double)t[i],omegai,delta,2,eps,&interp_type[nf-3],1))/pi;
   }
   for (unsigned long i=0;i<=1;i++) {
     Zxdipdi[nf-2+i]=newZxdipdi[i+1];
     Zxquaddi[nf-2+i]=newZxquaddi[i+1];
     Zlongdi[nf-2+i]=newZlongdi[i+1];
   }
   
   err=DBL_MIN;
   for (unsigned long i=0; i<=nz-1; i++) {
     err=max(abs(Wakexdip[i]-Wakexdipold[i]),err);
     err=max(abs(Wakexquad[i]-Wakexquadold[i]),err);
     err=max(factlong*abs(Wakelong[i]-Wakelongold[i]),err);
   }
   cout << "Max. error between two last maximum frequencies chosen : " << err << "\n";cout.flush();
   condition_freqmax=(err>=tol);
   
   for (unsigned long i=0; i<=nz-1; i++) {
     Wakexdipold[i]=Wakexdip[i];
     Wakexquadold[i]=Wakexquad[i];
     Wakelongold[i]=Wakelong[i];
   }

 }
 
 printf("Final wake computation with %ld frequencies\n",impmem);cout.flush();
 // compute the pchip derivatives on the full frequency range in memory
 pchip_derivatives(Zxdipdi,freqmem,Zxdipmem,impmem);
 pchip_derivatives(Zxquaddi,freqmem,Zxquadmem,impmem);
 pchip_derivatives(Zlongdi,freqmem,Zlongmem,impmem);
 // computes omega and delta between successive omegas, and choose interpolation type 
 // on each interval (pchip)
 delete[] omegai;delete[] delta;
 omegai=new long double[impmem];delta=new long double[impmem-1];
 for (unsigned long i=0; i<=impmem-1; i++) {
   omegai[i]=2.L*pi*(long double)freqmem[i];Zxdipdi[i]=Zxdipdi[i]/(2.*(double)pi);
   Zxquaddi[i]=Zxquaddi[i]/(2.*(double)pi);Zlongdi[i]=Zlongdi[i]/(2.*(double)pi);
 }
 for (unsigned long i=0; i<impmem-1; i++) {
   delta[i]=omegai[i+1]-omegai[i];interp_type[i]=0;
 }
 // computes the final wakes using all frequencies in memory
 for (unsigned long i=0; i<=nz-1; i++) {
   Wakexdip[i]=std::imag(fourier_integral_inf(Zxdipmem,Zxdipdi,0.,(long double)t[i],omegai,delta,impmem,eps,interp_type,1))/pi;
   Wakexquad[i]=std::imag(fourier_integral_inf(Zxquadmem,Zxquaddi,0.,(long double)t[i],omegai,delta,impmem,eps,interp_type,1))/pi;
   Wakelong[i]=std::real(fourier_integral_inf(Zlongmem,Zlongdi,0.,(long double)t[i],omegai,delta,impmem,eps,interp_type,1))/pi;
 }
 // compare with previous version on the converged mesh
 err=DBL_MIN;
 for (unsigned long i=0; i<=nz-1; i++) {
   err=max(abs(Wakexdip[i]-Wakexdipold[i]),err);
   err=max(abs(Wakexquad[i]-Wakexquadold[i]),err);
   err=max(abs(Wakelong[i]-Wakelongold[i]),err);
 }
 cout << "Max. error between two last frequency meshes : " << err << "\n";
  
   
   
 // writes the final impedances on the final mesh chosen
 for (unsigned long i=0; i<=nf-1; i++) {
   fprintf(filZxdip,"%13.8e %13.8e %13.8e\n",freq[i],yokoya[1]*Zxdipfi[i].real(),yokoya[1]*Zxdipfi[i].imag());
   fprintf(filZydip,"%13.8e %13.8e %13.8e\n",freq[i],yokoya[2]*Zxdipfi[i].real(),yokoya[2]*Zxdipfi[i].imag());
   if ( (yokoya[0]==1)&& ( (yokoya[1]==1)&&(yokoya[2]==1) ) ) {
     // case of axisymmetric geometry -> small quadrupolar impedance (see ref. cited at the beginning)
     fprintf(filZxquad,"%13.8e %13.8e %13.8e\n",freq[i],Zxquadfi[i].real(),Zxquadfi[i].imag());
     fprintf(filZyquad,"%13.8e %13.8e %13.8e\n",freq[i],Zxquadfi[i].real(),Zxquadfi[i].imag());
   } else {
     fprintf(filZxquad,"%13.8e %13.8e %13.8e\n",freq[i],yokoya[3]*Zxdipfi[i].real(),yokoya[3]*Zxdipfi[i].imag());
     fprintf(filZyquad,"%13.8e %13.8e %13.8e\n",freq[i],yokoya[4]*Zxdipfi[i].real(),yokoya[4]*Zxdipfi[i].imag());
   }
   fprintf(filZlong,"%13.8e %13.8e %13.8e\n",freq[i],yokoya[0]*Zlongfi[i].real(),yokoya[0]*Zlongfi[i].imag());
 }

 // writes the impedances with the finest possible mesh
 for (unsigned long i=0; i<=impmem-1; i++) {
   fprintf(filZxdip2,"%13.8e %13.8e %13.8e\n",freqmem[i],yokoya[1]*Zxdipmem[i].real(),yokoya[1]*Zxdipmem[i].imag());
   fprintf(filZydip2,"%13.8e %13.8e %13.8e\n",freqmem[i],yokoya[2]*Zxdipmem[i].real(),yokoya[2]*Zxdipmem[i].imag());
   if ( (yokoya[0]==1)&& ( (yokoya[1]==1)&&(yokoya[2]==1) ) ) {
     // case of axisymmetric geometry -> small quadrupolar impedance (see ref. cited at the beginning)
     fprintf(filZxquad2,"%13.8e %13.8e %13.8e\n",freqmem[i],Zxquadmem[i].real(),Zxquadmem[i].imag());
     fprintf(filZyquad2,"%13.8e %13.8e %13.8e\n",freqmem[i],Zxquadmem[i].real(),Zxquadmem[i].imag());
   } else {
     fprintf(filZxquad2,"%13.8e %13.8e %13.8e\n",freqmem[i],yokoya[3]*Zxdipmem[i].real(),yokoya[3]*Zxdipmem[i].imag());
     fprintf(filZyquad2,"%13.8e %13.8e %13.8e\n",freqmem[i],yokoya[4]*Zxdipmem[i].real(),yokoya[4]*Zxdipmem[i].imag());
   }
   fprintf(filZlong2,"%13.8e %13.8e %13.8e\n",freqmem[i],yokoya[0]*Zlongmem[i].real(),yokoya[0]*Zlongmem[i].imag());
 }

 // writes the final wakes 
 for (unsigned long i=0; i<=nz-1; i++) {
   fprintf(filWxdip,"%13.8e %13.8e\n",z[i],yokoya[1]*Wakexdipold[i]);
   fprintf(filWydip,"%13.8e %13.8e\n",z[i],yokoya[2]*Wakexdipold[i]);
   if ( (yokoya[0]==1)&& ( (yokoya[1]==1)&&(yokoya[2]==1) ) ) {
     // case of axisymmetric geometry -> small quadrupolar impedance (see ref. cited at the beginning)
     fprintf(filWxquad,"%13.8e %13.8e\n",z[i],Wakexquadold[i]);
     fprintf(filWyquad,"%13.8e %13.8e\n",z[i],Wakexquadold[i]);
   } else {
     fprintf(filWxquad,"%13.8e %13.8e\n",z[i],yokoya[3]*Wakexdipold[i]);
     fprintf(filWyquad,"%13.8e %13.8e\n",z[i],yokoya[4]*Wakexdipold[i]);
   }
   fprintf(filWlong,"%13.8e %13.8e\n",z[i],yokoya[0]*Wakelongold[i]);
 } 
 
 // writes the wakes with the finest possible mesh
 for (unsigned long i=0; i<=nz-1; i++) {
   fprintf(filWxdip2,"%13.8e %13.8e\n",z[i],yokoya[1]*Wakexdip[i]);
   fprintf(filWydip2,"%13.8e %13.8e\n",z[i],yokoya[2]*Wakexdip[i]);
   if ( (yokoya[0]==1)&& ( (yokoya[1]==1)&&(yokoya[2]==1) ) ) {
     // case of axisymmetric geometry -> small quadrupolar impedance (see ref. cited at the beginning)
     fprintf(filWxquad2,"%13.8e %13.8e\n",z[i],Wakexquad[i]);
     fprintf(filWyquad2,"%13.8e %13.8e\n",z[i],Wakexquad[i]);
   } else {
     fprintf(filWxquad2,"%13.8e %13.8e\n",z[i],yokoya[3]*Wakexdip[i]);
     fprintf(filWyquad2,"%13.8e %13.8e\n",z[i],yokoya[4]*Wakexdip[i]);
   }
   fprintf(filWlong2,"%13.8e %13.8e\n",z[i],yokoya[0]*Wakelong[i]);
 } 
 
 //finalization
 gsl_integration_workspace_free(w);
 
 delete[] freqmem;delete[] Zxdipmem;
 delete[] Zxquadmem;delete[] Zlongmem;
 
 delete[] freq;delete[] interp_type; delete[] inte;
 delete[] Zxdipfi;delete[] Zxdipdi;
 delete[] Zxquadfi;delete[] Zxquaddi;
 delete[] Zlongfi;delete[] Zlongdi;
 
 delete[] newfreq;
 delete[] newZxdipfi;delete[] newZxdipdi;
 delete[] newZxquadfi;delete[] newZxquaddi;
 delete[] newZlongfi;delete[] newZlongdi;
 
 delete[] omegai;delete[] delta;
 delete[] z;delete[] t;
 delete[] Wakexdip;
 delete[] Wakexquad;
 delete[] Wakelong;
 delete[] Wakexdipold;
 delete[] Wakexquadold;
 delete[] Wakelongold;

 fclose(filZxdip);
 fclose(filZydip);
 fclose(filZxquad);
 fclose(filZyquad);
 fclose(filZlong);
 fclose(filZxdip2);
 fclose(filZydip2);
 fclose(filZxquad2);
 fclose(filZyquad2);
 fclose(filZlong2);
 fclose(filWxdip);
 fclose(filWydip);
 fclose(filWxquad);
 fclose(filWyquad);
 fclose(filWlong);
 fclose(filWxdip2);
 fclose(filWydip2);
 fclose(filWxquad2);
 fclose(filWyquad2);
 fclose(filWlong2);
 
 time(&end);
 dif=difftime(end,start);
 printf("Elapsed time during calculation: %.2lf seconds\n",dif);

}
