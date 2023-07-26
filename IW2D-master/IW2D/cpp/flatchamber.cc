/*
 *  flatchamber.cc 
 *  
 *  by Nicolas Mounet (Nicolas.Mounet@cern.ch)
 *
 *  computes the impedance in a flat chamber (see CERN note by N. Mounet and E. Metral, 
 "Electromagnetic fields and beam coupling impedances in a multilayer flat chamber", 2010)
 
In input : typical input file is

Machine:	LHC
Relativistic Gamma:	479.6
Impedance Length in m:	1
Maximum order of non-linear terms:  0
Number of upper layers in the chamber wall:	2
Layer 1 inner half gap in mm:	18.375
Layer 1 DC resistivity (Ohm.m):	2e-10
Layer 1 relaxation time for resistivity (ps):	2.1
Layer 1 real part of dielectric constant:	1
Layer 1 magnetic susceptibility:	0
Layer 1 relaxation frequency of permeability (MHz):	Infinity
Layer 1 thickness in mm:	0.05
Layer 2 DC resistivity (Ohm.m):	7.2e-7
Layer 2 relaxation time for resistivity (ps):	0
Layer 2 real part of dielectric constant:	1.43
Layer 2 magnetic susceptibility:	0.02
Layer 2 relaxation frequency of permeability (MHz):	Infinity
Layer 2 thickness in mm:	Infinity
Top bottom symmetry (yes or no):	yes
Number of lower layers in the chamber wall:	2
Layer -1 inner half gap in mm:	18.375
Layer -1 DC resistivity (Ohm.m):	2e-10
Layer -1 relaxation time for resistivity (ps):	2.1
Layer -1 real part of dielectric constant:	1
Layer -1 magnetic susceptibility:	0
Layer -1 relaxation frequency of permeability (MHz):	Infinity
Layer -1 thickness in mm:	0.05
Layer -2 DC resistivity (Ohm.m):	7.2e-7
Layer -2 relaxation time for resistivity (ps)
:	0
Layer -2 real part of dielectric constant:	1.43
Layer -2 magnetic susceptibility:	0.02
Layer -2 relaxation frequency of permeability (MHz):	Infinity
Layer -2 thickness in mm:	Infinity
start frequency exponent (10^) in Hz:	0
stop frequency exponent (10^) in Hz:	14
linear (1) or logarithmic (0) or both (2) frequency scan:	2
sampling frequency exponent (10^) in Hz (for linear):	8
Number of points per decade (for log):	20
when both, fmin of the refinement (in THz):	5.5
when both, fmax of the refinement (in THz):	7
when both, number of points in the refinement:	100
added frequencies [Hz]:	1e-3 1e-1
Comments for the output files names:	_some_element

The order of the lines can be whatever, but the exact sentences and the TAB before the parameter
indicated, are necessary. If top-bottom symmetry is set (with "yes" or "y" or "1") the lower layers (with a
minus sign) are ignored. Also if there are more layers than indicated by the number of upper (lower) layers,
the additional ones are ignored. The last layer is always assumed to go to infinity.

Other options for layer properties (for any layer):

Layer 1 dielectric tan(delta):  0.1
Layer 1 frequency dependent conductivity:       conductivity_vs_freq_hBN_example.dat
Layer 1 frequency dependent relative complex permittivity:      relative_permittivity_vs_freq_hBN_example.dat
Layer 1 frequency dependent relative complex permeability:      ferrite_8C11_muprime_musecond.dat

These are mutually exclusive: one cannot use a tan(delta) with a resistivity or relaxation time
(it will raise an error), nor a frequency dependent conductivity - resp. permittivity - 
with any other permittivity related quantity (resistivity, relaxation time, dielectric constant, 
tan(delta), freq. dependent permittivity - resp. freq. dependent conductivity).
The same goes for the permeability, but permeability/permittivity properties
are independent from each other so not mutually exclusive.

Freq. dependent quantities are defined in files looking like

Frequency [Hz]	eps_prime	eps_second
0.	    3.12385	        -0.203868386315
106.082	3.12385	        -0.203868386315
106.289	3.12356849182	-0.203850014556
142.51	3.07431	        -0.183599648141
148.657	3.07207572348	-0.180577225441
203.092	3.05229	        -0.163971130987
207.913	3.05058721854	-0.162512712659

(here for the relative complex permittivity). The first line is always
ignored (headers). For the conductivity, a complex sigma should be given,
but one can set to zero the imaginary part if needed (still, there should always be 3 columns).
Units are S.I (relatice complex permittivity/permeability are dimensionless,
conductivity is in S/m), and beware of the sign convension - eps'' and mu'' 
should be negative for positive freq. (see N. Mounet's PhD thesis and references therehin).

NOTE: freq. dependent quantities are linearly extrapolated ouside the frequency
range defined in the file, using the slope from the edges. So to avoid any
unphysical value, better extend the frequency range artificially with well
chosen values.

In output one gives six files with the impedances (longitudinal, x dipolar, y dipolar,
x quadrupolar, y quadrupolar and y constant term). Each have 3 columns :frequency, real part and
imaginary part.
If "Maximum order of non-linear terms" is strictly higher than 0, it will also
give additional files of the same format, of name Z<plane><a><b><c><d>*.dat
with the general impedance term in front of x1^a y1^b x^c y^d , for all
orders such that a+b+c+d <= Max. order

 */

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sstream>
#include <sys/time.h>

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
#include <globals.h>

const unsigned long MAXMEM=50000;  // Maximum number of elements of arrays with etas and chis

/*
#define  Precision	200  // Number of bits for the mantissa of all numbers (classic doubles have 53 bits)
#define  MAXLINES	100  // Maximum number of lines in the input file (correspond to 5 layers in top and
			    // bottom parts)
#define  MAXCHAR	200  // Maximum number of characters in one line in the input file
#define  MAXCHARFILE	200  // Maximum number of characters in the output files name extension
#define  MAXMEM		50000 // Maximum number of elements of arrays with etas and chis

const  amp::ampf<Precision> C=299792458;    // velocity of light [m/s]
const amp::ampf<Precision> mu0=4e-7*amp::pi<Precision>(); // mu0
const amp::ampf<Precision> eps0=1/(mu0*amp::sqr(C)); // eps0

struct params {unsigned int m; unsigned int n;
  	unsigned int N; unsigned int M;
	ap::template_1d_array< amp::campf<Precision> > eps1;
  	ap::template_1d_array< amp::campf<Precision> > mu1;
  	ap::template_1d_array< amp::campf<Precision> > nu2;
	ap::template_1d_array< amp::campf<Precision> > eps1ratio;
  	ap::template_1d_array< amp::campf<Precision> > mu1ratio;
  	ap::template_1d_array< amp::campf<Precision> > nu2ratio;
	ap::template_1d_array< amp::ampf<Precision> > b;
	ap::template_1d_array< amp::campf<Precision> > eps1m;
  	ap::template_1d_array< amp::campf<Precision> > mu1m;
  	ap::template_1d_array< amp::campf<Precision> > nu2m;
	ap::template_1d_array< amp::campf<Precision> > eps1mratio;
  	ap::template_1d_array< amp::campf<Precision> > mu1mratio;
  	ap::template_1d_array< amp::campf<Precision> > nu2mratio;
	ap::template_1d_array< amp::ampf<Precision> > bm;
	amp::ampf<Precision> beta;
	amp::ampf<Precision> kovergamma;};
// NB: eps1, eps1m, mu1 and mu1m are actually the inverse of the 'real' eps1, etc.

// global arrays with memory of etas and chis
ap::template_1d_array< amp::ampf<Precision> > kxmem;
ap::template_1d_array< amp::campf<Precision> > eta1mem,eta2mem,chi1mem,chi2mem;
unsigned long mem; // current number of elements of kxmem, eta1mem, etc.
*/

/**************************************************
 *** 			main program		***
 ***				              	***
 **************************************************/

int main ()

{
 
 char *endline;
 char output[MAXCHARFILE],Zxdipoutput[MAXCHARFILE+10],Zydipoutput[MAXCHARFILE+10],Zlongoutput[MAXCHARFILE+10];
 char Zxquadoutput[MAXCHARFILE+10],Zyquadoutput[MAXCHARFILE+10],Zycstoutput[MAXCHARFILE+10],Input[MAXCHARFILE+10];
 char Zoutput[MAXCHARFILE+10],plane_char,buffer[10];
 std::string data[MAXLINES],machine,topbot,commentoutput,dummy2,unit,epsfile,sigmafile,mufile;
 FILE *filZxdip, *filZydip, *filZxquad, *filZyquad, *filZycst, *filZlong, *filInput, *filZ;
 unsigned int N,M,dummy0,fact,aplusc2; // number of upper and lower layers, then dummy parameter, then some factor, and (a+c)/2
 unsigned int m,n,n_input,n_added,nf,nfdup; /* indices of alphamn (azimuthal mode numbers), number of input lines,
 					number of individually added frequencies, total number of
					frequencies in the scan, number of duplicate frequencies */
 unsigned int flag_topbotsym,typescan,nflog,nflin,order; /* flag for top-bottom symmetry (1 if such a symmetry), type of frequency scan
 		number of freq. per decade, number of freq. in a lin. scan inside the log scan,
		maximum order for non-linear terms*/
 ap::template_1d_array< amp::campf<PRECISION> > eps1,eps1m,mu1,mu1m; /* eps1 and mu1 for upper (without m at
 					the end) and lower layers (with m at the end) */
 ap::template_1d_array< amp::ampf<PRECISION> > b,bm,thick,thickm; // position of upper and lower boundaries; thickness of the layers 
 ap::template_1d_array< amp::ampf<PRECISION> > rho,tau,epsb,tandelta,chi,fmu,rhom,taum,epsbm,tandeltam,chim,fmum; /* layers
 						properties (DC resistivity, resistivity relaxation time,
						dielectric constant, tan(delta), magnetic susceptibility=mu_r-1,
						relaxation frequency of permeability)*/
 amp::ampf<PRECISION> omega,beta,k,gamma,kovergamma,u,dummy3; // parameters
 amp::campf<PRECISION> jimagMP; // imaginary constant in multiprecision
 gsl_integration_workspace *w;
 size_t limit=1000; // limit (number of intervals) for gsl integration algorithm
 double x,y,L,fminlog,fmaxlog,fminlin,fmaxlin,fsamplin,fadded[15],*freq,dif,dummy1,coef;
 std::complex<double> *Z,Zxdip,Zydip,Zxquad,Zyquad,Zlong,Zycst,cst;
 std::complex<double> alpha,alpha00,alpha01,alpha02,alpha11; // alphamn constants
 double **eps_frequencies, **mu_frequencies, **epsm_frequencies, **mum_frequencies;
 std::complex<double> **eps_values, **mu_values, **epsm_values, **mum_values;
 int *eps_vs_freq_kind, *mu_vs_freq_kind, *epsm_vs_freq_kind, *mum_vs_freq_kind;
 long *n_eps_frequencies, *n_mu_frequencies, *n_epsm_frequencies, *n_mum_frequencies;
 size_t found;
 time_t start,end; // times
 memorycontainer memory(MAXMEM);	
 // start time
 time(&start);

 // imaginary constant in multiprecision
 jimagMP.x=0;jimagMP.y=1;

 // default values of the parameters (in case)
 flag_topbotsym=1;
 order=0;
 typescan=0;
 fminlog=2;fmaxlog=13;nflog=10;n_added=0;
 fsamplin=8;fminlin=1;fmaxlin=2;nflin=100;nf=0;
 N=2;M=2;L=1.;gamma="479.6";
 
 // read input file
 // first read everything to identify the strings in front of each parameters
 n_input=0;
 while (std::cin.eof()==0) {
   std::getline (std::cin,data[n_input+1]);
   /* next line is when the input file comes from windows or else and has some ^M characters in the
   end of each line */
   //data[n_input+1]=data[n_input+1].substr(0,data[n_input+1].length()-1);
   //std::cout << data[n_input+1] << '\n';
   n_input++;
 }
 n_input--;
 //printf("n_input: %d\n",n_input);
 // identify each argument
 for (unsigned int i=1; i<=n_input; i++) {
   read_input(data[i],"Machine",dummy0,dummy1,machine,dummy3,2);
   read_input(data[i],"Relativistic Gamma",dummy0,dummy1,dummy2,gamma,3);
   read_input(data[i],"Impedance Length in m",dummy0,L,dummy2,dummy3,1);
   read_input(data[i],"Number of upper layers",N,dummy1,dummy2,dummy3,0);
   read_input(data[i],"Number of lower layers",M,dummy1,dummy2,dummy3,0);
   read_input(data[i],"Top bottom symmetry (yes or no)",dummy0,dummy1,topbot,dummy3,2);
   read_input(data[i],"Maximum order of non-linear terms",order,dummy1,dummy2,dummy3,0);
   read_input(data[i],"start frequency exponent (10^)",dummy0,fminlog,dummy2,dummy3,1);
   read_input(data[i],"stop frequency exponent (10^)",dummy0,fmaxlog,dummy2,dummy3,1);
   read_input(data[i],"linear (1) or logarithmic (0) or both (2) frequency scan",typescan,dummy1,dummy2,dummy3,0);
   read_input(data[i],"sampling frequency exponent (10^) in Hz (for linear)",dummy0,fsamplin,dummy2,dummy3,1);
   read_input(data[i],"Number of points per decade (for log)",nflog,dummy1,dummy2,dummy3,0);
   read_input(data[i],"when both, fmin of the refinement",dummy0,fminlin,dummy2,dummy3,1);
   read_input(data[i],"when both, fmax of the refinement",dummy0,fmaxlin,dummy2,dummy3,1);
   read_input(data[i],"when both, number of points in the refinement",nflin,dummy1,dummy2,dummy3,0);
   read_input(data[i],"Comments for the output files names",dummy0,dummy1,commentoutput,dummy3,2);
   read_input_double_list(data[i], "added frequencies", fadded, n_added);
 }
 //printf("%s %s %d %d \n",machine.c_str(),topbot.c_str(),N,flag_topbotsym);
 //printf("%13.8e %13.8e %d\n",fadded[1],fadded[2],n_added);
 //printf("%13.8e %13.8e %d %d\n",fminlog,fmaxlog,nflog,typescan);

 // flag for top bottom symmetry (1 if there is such a symmetry)
 flag_topbotsym= ((strcmp(topbot.c_str(),"yes")==0 || strcmp(topbot.c_str(),"y")==0) || strcmp(topbot.c_str(),"1")==0);

 
 eps1.setbounds(1,N+1);mu1.setbounds(1,N+1);b.setbounds(1,N+1);thick.setbounds(1,N);
 rho.setbounds(1,N+1);tau.setbounds(1,N+1);epsb.setbounds(1,N+1);tandelta.setbounds(1,N+1);
 chi.setbounds(1,N+1);fmu.setbounds(1,N+1);
 eps_frequencies=new double*[N+2];
 eps_values=new std::complex<double>*[N+2];
 mu_frequencies=new double*[N+2];
 mu_values=new std::complex<double>*[N+2];
 eps_vs_freq_kind=new int[N+2];
 n_eps_frequencies=new long[N+2];
 mu_vs_freq_kind=new int[N+2];
 n_mu_frequencies=new long[N+2];
 if (flag_topbotsym) M=N;
 eps1m.setbounds(1,M+1);mu1m.setbounds(1,M+1);bm.setbounds(1,M+1);thickm.setbounds(1,M);
 rhom.setbounds(1,M+1);taum.setbounds(1,M+1);epsbm.setbounds(1,M+1);tandeltam.setbounds(1,M+1);
 chim.setbounds(1,M+1);fmum.setbounds(1,M+1);
 epsm_frequencies= new double*[M+2];
 epsm_values= new std::complex<double>*[M+2];
 mum_frequencies= new double*[M+2];
 mum_values= new std::complex<double>*[M+2];
 epsm_vs_freq_kind=new int[M+2];
 n_epsm_frequencies=new long[M+2];
 mum_vs_freq_kind=new int[M+2];
 n_mum_frequencies=new long[M+2];
 
 // default values of the layers properties
 for (unsigned int p=1; p<=N+1; p++){
   // upper layers
   rho(p)="Infinity";
   tau(p)=0;
   epsb(p)=1;
   tandelta(p)=0;
   chi(p)=0;
   fmu(p)="Infinity";
   eps_vs_freq_kind[p]=0;
   mu_vs_freq_kind[p]=0;
 }
 for (unsigned int p=1; p<=M+1; p++){
   // lower layers
   rhom(p)="Infinity";
   taum(p)=0;
   epsbm(p)=1;
   tandeltam(p)=0;
   chim(p)=0;
   fmum(p)="Infinity";
   epsm_vs_freq_kind[p]=0;
   mum_vs_freq_kind[p]=0;
 }
   
 // in case nothing is defined, we still compute a single layer of 2mm-separated graphite plates
 //rho(2)="1e5";tau(2)=0;epsb(2)=1;chi(2)=0;fmu(2)="Infinity";b(1)=2;thick(1)="Infinity";

 // find inner half gap(s) of the chamber
 for (unsigned int i=1; i<=n_input; i++) {
   read_input(data[i],"Layer 1 inner half gap in mm",dummy0,dummy1,dummy2,b(1),3);
   if (flag_topbotsym) bm(1)=b(1);
   else {
     read_input(data[i],"Layer -1 inner half gap in mm",dummy0,dummy1,dummy2,bm(1),3);
   }
 }
 bm(1)=-bm(1);
 //printf("%s %s \n",b(1).toDec().c_str(),bm(1).toDec().c_str());
 
 // find all the layers properties
 for (unsigned int i=1; i<=n_input; i++) {
   // upper layers
   for (unsigned int p=1; p<=N; p++){
     epsfile="";
     sigmafile="";
     mufile="";
     read_input_layer(data[i],"DC resistivity",p,rho(p+1));
     read_input_layer(data[i],"relaxation time for resistivity (ps)",p,tau(p+1));
     read_input_layer(data[i],"real part of dielectric constant",p,epsb(p+1));
     read_input_layer(data[i],"dielectric tan(delta)",p,tandelta(p+1));
     read_input_layer(data[i],"magnetic susceptibility",p,chi(p+1));
     read_input_layer(data[i],"relaxation frequency of permeability",p,fmu(p+1));
     read_input_layer_filename(data[i],"frequency dependent relative complex permittivity",p,epsfile);
     read_input_layer_filename(data[i],"frequency dependent conductivity",p,sigmafile);
     read_input_layer_filename(data[i],"frequency dependent relative complex permeability",p,mufile);
     read_input_layer(data[i],"thickness in mm",p,thick(p));
     if (!epsfile.empty()) {
       if (eps_vs_freq_kind[p+1]==2) {
         printf("Cannot use both a freq. dependent conductivity and a freq. dependent permittivity!\n");
         exit(2);
       }
       eps_vs_freq_kind[p+1]=1;
       read_complex_file_length(epsfile,n_eps_frequencies[p+1]);
       //eps_frequencies[p+1] = new double[n_eps_frequencies[p+1]+1];
       //eps_values[p+1] = new std::complex<double>[n_eps_frequencies[p+1]+1];
       read_complex_file(epsfile,eps_frequencies[p+1],eps_values[p+1],n_eps_frequencies[p+1]);
     }
     if (!sigmafile.empty()) {
       if (eps_vs_freq_kind[p+1]==1) {
         printf("Cannot use both a freq. dependent permittivity and a freq. dependent conductivity!\n");
         exit(2);
       }
       eps_vs_freq_kind[p+1]=2;
       read_complex_file_length(sigmafile,n_eps_frequencies[p+1]);
       eps_frequencies[p+1] = new double[n_eps_frequencies[p+1]+1];
       eps_values[p+1] = new std::complex<double>[n_eps_frequencies[p+1]+1];
       read_complex_file(sigmafile,eps_frequencies[p+1],eps_values[p+1],n_eps_frequencies[p+1]);
       //for (long j=0;j<=n_eps_frequencies[p+1]-1;j++) printf("%ld %13.8e %13.8e %13.8e\n",j,eps_frequencies[p+1][j],eps_values[p+1][j].real(),eps_values[p+1][j].imag());
     }
     if (!mufile.empty()) {
       mu_vs_freq_kind[p+1]=1;
       read_complex_file_length(mufile,n_mu_frequencies[p+1]);
       mu_frequencies[p+1] = new double[n_mu_frequencies[p+1]+1];
       mu_values[p+1] = new std::complex<double>[n_mu_frequencies[p+1]+1];
       read_complex_file(mufile,mu_frequencies[p+1],mu_values[p+1],n_mu_frequencies[p+1]);
     }
   }
   // lower layers
   if (flag_topbotsym) {
     for (unsigned int p=1; p<=M; p++){
       rhom(p+1)=rho(p+1);
       taum(p+1)=tau(p+1);
       epsbm(p+1)=epsb(p+1);
       tandeltam(p+1)=tandelta(p+1);
       chim(p+1)=chi(p+1);
       fmum(p+1)=fmu(p+1);
       epsm_vs_freq_kind[p+1]=eps_vs_freq_kind[p+1];
       n_epsm_frequencies[p+1]=n_eps_frequencies[p+1];
       if (epsm_vs_freq_kind[p+1]>=1) {
         //epsm_frequencies[p+1] = new double[n_epsm_frequencies[p+1]+1];
         //epsm_values[p+1] = new std::complex<double>[n_epsm_frequencies[p+1]+1];
         epsm_frequencies[p+1]=eps_frequencies[p+1];
         epsm_values[p+1]=eps_values[p+1];
       }
       mum_vs_freq_kind[p+1]=mu_vs_freq_kind[p+1];
       n_mum_frequencies[p+1]=n_mu_frequencies[p+1];
       if (mum_vs_freq_kind[p+1]>=1) {
         //mum_frequencies[p+1] = new double[n_mum_frequencies[p+1]+1];
         //mum_values[p+1] = new std::complex<double>[n_mu_frequencies[p+1]+1];
         mum_frequencies[p+1]=mu_frequencies[p+1];
         mum_values[p+1]=mu_values[p+1];
       }
       thickm(p)=thick(p);
     }
   }
   else {
     for (unsigned int p=1; p<=M; p++){
       epsfile="";
       sigmafile="";
       mufile="";
       read_input_layer(data[i],"DC resistivity",-p,rhom(p+1));
       read_input_layer(data[i],"relaxation time for resistivity (ps)",-p,taum(p+1));
       read_input_layer(data[i],"real part of dielectric constant",-p,epsbm(p+1));
       read_input_layer(data[i],"dielectric tan(delta)",-p,tandeltam(p+1));
       read_input_layer(data[i],"magnetic susceptibility",-p,chim(p+1));
       read_input_layer(data[i],"relaxation frequency of permeability",-p,fmum(p+1));
       read_input_layer_filename(data[i],"frequency dependent relative complex permittivity",-p,epsfile);
       read_input_layer_filename(data[i],"frequency dependent conductivity",-p,sigmafile);
       read_input_layer_filename(data[i],"frequency dependent relative complex permeability",-p,mufile);
       read_input_layer(data[i],"thickness in mm",-p,thickm(p));
       if (!epsfile.empty()) {
         if (epsm_vs_freq_kind[p+1]==2) {
           printf("Cannot use both a freq. dependent conductivity and a freq. dependent permittivity!\n");
           exit(2);
         }
         epsm_vs_freq_kind[p+1]=1;
         read_complex_file_length(epsfile,n_epsm_frequencies[p+1]);
         epsm_frequencies[p+1] = new double[n_epsm_frequencies[p+1]+1];
         epsm_values[p+1] = new std::complex<double>[n_epsm_frequencies[p+1]+1];
         read_complex_file(epsfile,epsm_frequencies[p+1],epsm_values[p+1],n_epsm_frequencies[p+1]);
       }
       if (!sigmafile.empty()) {
         if (epsm_vs_freq_kind[p+1]==1) {
           printf("Cannot use both a freq. dependent permittivity and a freq. dependent conductivity!\n");
           exit(2);
         }
         epsm_vs_freq_kind[p+1]=2;
         read_complex_file_length(sigmafile,n_epsm_frequencies[p+1]);
         epsm_frequencies[p+1] = new double[n_epsm_frequencies[p+1]+1];
         epsm_values[p+1] = new std::complex<double>[n_epsm_frequencies[p+1]+1];
         read_complex_file(sigmafile,epsm_frequencies[p+1],epsm_values[p+1],n_epsm_frequencies[p+1]);
       }
       if (!mufile.empty()) {
         mum_vs_freq_kind[p+1]=1;
         read_complex_file_length(mufile,n_mum_frequencies[p+1]);
         mum_frequencies[p+1] = new double[n_mum_frequencies[p+1]+1];
         mum_values[p+1] = new std::complex<double>[n_mum_frequencies[p+1]+1];
         read_complex_file(mufile,mum_frequencies[p+1],mum_values[p+1],n_mum_frequencies[p+1]);
       }
     }
   }
   //printf("%s %s %s %d %d %d %d %d\n",epsfile.c_str(),sigmafile.c_str(),mufile.c_str(),epsfile.empty(),sigmafile.empty(),mufile.empty(),i,N);
 }

 // units conversion to SI (fminlin and fmaxlin were in THz, tau was in ps, b in mm and fmu in MHz)
 fminlin*=1.e12;fmaxlin*=1.e12;
 b(1)=b(1)*1e-3;bm(1)=bm(1)*1e-3;
 for (unsigned int p=2; p<=N+1; p++){
   tau(p)=tau(p)*1e-12;
   b(p)=thick(p-1)*1e-3+b(p-1);
   fmu(p)=fmu(p)*1e6;
   //printf("%lg %lg %lg %lg %lg %lg \n",double(rho(p).toDouble()),double(tau(p).toDouble()),double(epsb(p).toDouble()),
   //	double(chi(p).toDouble()),double(fmu(p).toDouble()),double(b(p).toDouble()));
 }
 for (unsigned int p=2; p<=M+1; p++){
   taum(p)=taum(p)*1e-12;
   bm(p)=bm(p-1)-thickm(p-1)*1e-3;
   fmum(p)=fmum(p)*1e6;
   //printf("%lg %lg %lg %lg %lg %lg \n",double(rhom(p).toDouble()),double(taum(p).toDouble()),double(epsbm(p).toDouble()),
   //	double(chim(p).toDouble()),double(fmum(p).toDouble()),double(bm(p).toDouble()));
 }
  
 // first layer (inside the chamber) is always vacuum
 eps1(1)=1;mu1(1)=1;
 eps1m(1)=1;mu1m(1)=1;
 // relativistic velocity factor beta
 beta=amp::sqrt(1-1/amp::sqr(gamma));
 //printf("%lg %lg\n",double(eps1(1).x.toDouble()),double(eps1(1).y.toDouble()));
 
 // construct the frequency scan
 // first estimation of the number of frequencies (for memory allocation)
 switch(typescan) {
   case 0:
     nf=(int)ceil((fmaxlog-fminlog)*(double)nflog)+1+n_added;
     break;
   case 1:
     nf=(int)ceil((pow(10.,fmaxlog)-pow(10.,fminlog))/pow(10.,fsamplin))+1+n_added;
     break;
   case 2:
     nf=(int)ceil((fmaxlog-fminlog)*(double)nflog)+1+nflin+1+n_added;
     break;
   }
 freq=new double[nf];
 
 // constructs unsorted version of the array freq
 nf=0;freq[0]=-1.;
 if (typescan==1) {
   do {
     freq[nf]=pow(10.,fminlog)+(double)nf * pow(10.,fsamplin);
     nf++;
     } while(freq[nf-1]<pow(10.,fmaxlog));
   }
 else {
   do {
     freq[nf]=pow(10.,fminlog+(double)nf/(double)nflog);
     nf++;
     } while(freq[nf-1]<pow(10.,fmaxlog));
   if (typescan==2) {
     for (unsigned int i=0; i<=nflin; i++) {
       freq[nf]=fminlin+(double)(i) * (fmaxlin-fminlin)/(double)nflin;
       nf++;
       }
     }
   }
 for (unsigned int i=0;i<n_added; i++) {
   freq[nf]=fadded[i];
   nf++;
   }
 // nf is now the real number of frequencies
 
 // sort the frequencies
 gsl_sort(freq, 1, nf);
 
 // remove duplicate frequencies
 nfdup=0;
 for (int i=nf-2;i>=0; i--) {
   if ( (freq[i+1]-freq[i]) ==0 ) {
       freq[i]=DBL_MAX;
       nfdup++;
     }
   }
 gsl_sort(freq, 1, nf);
 nf=nf-nfdup;
 //printf("%d\n",nf);
 //for (unsigned int i=0;i<nf;i++) printf("%d %13.8e\n",i,freq[i]);
 
 
 // set the output files names
 //sprintf(output,"W%s_%dlayersup_%dlayersdown%.2lfmm%s",machine,N,M,
 //	1e3*double(amp::ampf<Precision>(b(1)).toDouble()),commentoutput);
 sprintf(output,"W%s_%dlayersup_%dlayersdown%.2lfmm%s",machine.c_str(),N,M,
	1e3*double(amp::ampf<PRECISION>(b(1)).toDouble()),commentoutput.c_str());
 sprintf(Zxdipoutput,"Zxdip%s.dat",output);
 sprintf(Zydipoutput,"Zydip%s.dat",output);
 sprintf(Zxquadoutput,"Zxquad%s.dat",output);
 sprintf(Zyquadoutput,"Zyquad%s.dat",output);
 sprintf(Zlongoutput,"Zlong%s.dat",output);
 sprintf(Zycstoutput,"Zycst%s.dat",output);
 sprintf(Input,"InputData%s.dat",output);
 
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
 filZycst=fopen(Zycstoutput,"w");
 fprintf(filZxdip,"Frequency [Hz]\tRe(Zxdip) [Ohm/m]\tIm(Zxdip) [Ohm/m]\n");
 fprintf(filZydip,"Frequency [Hz]\tRe(Zydip) [Ohm/m]\tIm(Zydip) [Ohm/m]\n");
 fprintf(filZxquad,"Frequency [Hz]\tRe(Zxquad) [Ohm/m]\tIm(Zxquad) [Ohm/m]\n");
 fprintf(filZyquad,"Frequency [Hz]\tRe(Zyquad) [Ohm/m]\tIm(Zyquad) [Ohm/m]\n");
 fprintf(filZlong,"Frequency [Hz]\tRe(Zlong) [Ohm]\tIm(Zlong) [Ohm]\n");
 fprintf(filZycst,"Frequency [Hz]\tRe(Zycst) [Ohm]\tIm(Zycst) [Ohm]\n");
 
 
 // workspace allocation for gsl adaptative integration
 w=gsl_integration_workspace_alloc(limit);
 gsl_set_error_handler_off();
 

 if (order<=0) {
   // only constant (for long.) and constant+linear (for trans.) terms are computed
   alpha00=std::complex<double>(0.,0.);
   alpha01=std::complex<double>(0.,0.);
   alpha02=std::complex<double>(0.,0.);
   alpha11=std::complex<double>(0.,0.);
   
   // impedance computation at each frequency: beginning of the loop
   for (unsigned int i=0; i<nf; i++) {
      
     memory.mem=0;  // initialize memory at each frequency
     omega=amp::twopi<PRECISION>()*freq[i];
     k=omega/(beta*C);
     kovergamma=k/gamma;
     
     // computes the layer properties for the angular freq. omega
     for (unsigned int p=2;p<=N+1; p++) {
       get_eps_and_mu_from_layer_properties(freq[i], rho(p), tau(p), epsb(p),
          tandelta(p), chi(p), fmu(p), eps_frequencies[p], eps_values[p],
          n_eps_frequencies[p], mu_frequencies[p], mu_values[p], n_mu_frequencies[p],
          eps_vs_freq_kind[p], mu_vs_freq_kind[p], eps1(p), mu1(p));
     
       //printf("%lg Hz: layer %d -> %lg m, eps1= %lg + %lg j, mu1= %lg + %lgj\n",freq[i],p,double(b(p).toDouble()),
       //   double(eps1(p).x.toDouble()),double(eps1(p).y.toDouble()),
       //   double(mu1(p).x.toDouble()),double(mu1(p).y.toDouble()));
       }
     for (unsigned int p=2;p<=M+1; p++) {
       get_eps_and_mu_from_layer_properties(freq[i], rhom(p), taum(p), epsbm(p),
          tandeltam(p), chim(p), fmum(p), epsm_frequencies[p], epsm_values[p],
          n_epsm_frequencies[p], mum_frequencies[p], mum_values[p], n_mum_frequencies[p],
          epsm_vs_freq_kind[p], mum_vs_freq_kind[p], eps1m(p), mu1m(p));

       //printf("%lg Hz: layer %d -> %lg m, eps1= %lg + %lg j, mu1= %lg + %lgj\n",freq[i],-p,double(bm(p).toDouble()),
       //   double(eps1m(p).x.toDouble()),double(eps1m(p).y.toDouble()),
       // double(mu1m(p).x.toDouble()),double(mu1m(p).y.toDouble()));
       }

     //printf("%s %s\n",b(1).toDec().c_str(),bm(1).toDec().c_str());
     //printf("%s %s %s\n",omega.toDec().c_str(),k.toDec().c_str(),kovergamma.toDec().c_str());


     // computes alpha00
     alpha00=alphamn(flag_topbotsym,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,0,0,limit,w, &memory);
     //printf("alpha00: %13.8e %13.8e\n",alpha00.real(),alpha00.imag());
 
     // computes alpha01
     alpha01=alphamn(flag_topbotsym,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,0,1,limit,w, &memory);
     //printf("alpha01: %13.8e %13.8e\n",alpha01.real(),alpha01.imag());
 
     // computes alpha02
     alpha02=alphamn(flag_topbotsym,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,0,2,limit,w, &memory);
     //printf("alpha02: %13.8e %13.8e\n",alpha02.real(),alpha02.imag());

     // computes alpha11
     alpha11=alphamn(flag_topbotsym,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,1,1,limit,w, &memory);
     //printf("alpha11: %13.8e %13.8e\n",alpha11.real(),alpha11.imag());
     
     // computes and writes the impedances
     cst=jimag*L*double(amp::ampf<PRECISION>(k*Z0/(beta*amp::sqr(gamma)*amp::twopi<PRECISION>())).toDouble());
     Zlong=cst*alpha00;
     Zycst=cst*alpha01/double(amp::ampf<PRECISION>(gamma).toDouble());
     cst=cst*double(amp::ampf<PRECISION>(k/(2*amp::sqr(gamma))).toDouble());
     Zxdip=cst*(alpha02-alpha00);
     Zydip=2.*cst*alpha11;
     Zxquad=-Zxdip;
     Zyquad=cst*(alpha02+alpha00);
     fprintf(filZlong,"%13.8e %13.8e %13.8e\n",freq[i],Zlong.real(),Zlong.imag());
     fprintf(filZycst,"%13.8e %13.8e %13.8e\n",freq[i],Zycst.real(),Zycst.imag());
     fprintf(filZxdip,"%13.8e %13.8e %13.8e\n",freq[i],Zxdip.real(),Zxdip.imag());
     fprintf(filZydip,"%13.8e %13.8e %13.8e\n",freq[i],Zydip.real(),Zydip.imag());
     fprintf(filZxquad,"%13.8e %13.8e %13.8e\n",freq[i],Zxquad.real(),Zxquad.imag());
     fprintf(filZyquad,"%13.8e %13.8e %13.8e\n",freq[i],Zyquad.real(),Zyquad.imag());

     }
     // end of loop on frequencies
   
   fclose(filZxdip);
   fclose(filZydip);
   fclose(filZxquad);
   fclose(filZyquad);
   fclose(filZlong);
   fclose(filZycst);
 
   /*for (int i=1;i<=10000; i++){
     u=i*0.01;
     alpha00=integrand(0, 0, N+1, M+1, eps1, mu1, b, eps1m, mu1m, bm, u, omega, beta, k, kovergamma);
     //alpha00=integrand(0, 0, N+1, M+1, eps1, mu1, b, eps1m, mu1m, bm, (1-u)/u, omega, beta, k,kovergamma)/amp::sqr(u);
     alpha02=integrand(0, 2, N+1, M+1, eps1, mu1, b, eps1m, mu1m, bm, u, omega, beta, k, kovergamma);
     //alpha02=integrand(0, 2, N+1, M+1, eps1, mu1, b, eps1m, mu1m, bm, (1-u)/u, omega, beta, k,kovergamma)/amp::sqr(u);
     printf("%13.8e %13.8e %13.8e\n",double(amp::ampf<Precision>(u).toDouble()),alpha02.real(),alpha02.imag());
   }*/

 }
 else {
   // we compute all non-linear terms up to total order 'order'
   Z = new std::complex<double>[nf];
   
   for (int plane=0; plane<=2; plane++) {
     // plane=0 -> long.; plane=1 -> hor.; plane=2 -> ver.
     if (plane==0) plane_char = 'l';
     if (plane==1) plane_char = 'x';
     if (plane==2) plane_char = 'y';
     
     for (int power_a=0; power_a<=order; power_a++) {
       for (int power_b=0; power_b<=order; power_b++) {
         for (int power_c=0; power_c<=order; power_c++) {
           for (int power_d=0; power_d<=order; power_d++) {

             //printf("plane=%c, a=%d, b=%d, c=%d, d=%d\n",plane_char,power_a,power_b,power_c,power_d);
             //printf("  -> (a+c+/-1) mod 2 =%d, a+b+c+d=%d, (a+b+c+d+/-1) mod 2=%d\n",(power_a+power_c+(plane==1))%2, \
                  power_a+power_b+power_c+power_d,((power_a+power_b+power_c+power_d+(plane>0))%2)*flag_topbotsym);

             if ( ((power_a+power_c+(plane==1))%2 == 0) && ((power_a+power_b+power_c+power_d)<=order) \
                  && (((power_a+power_b+power_c+power_d+(plane>0))%2)*flag_topbotsym==0) ) {

               if ( (power_a+power_b+power_c+power_d) == 0) unit="";
               else {
                 sprintf(buffer,"/m^%d",power_a+power_b+power_c+power_d);
                 unit = buffer;
               }
               
               //printf("selected, plane=%c, a=%d, b=%d, c=%d, d=%d, unit=%s\n",plane_char,power_a,power_b,power_c,power_d,unit.c_str());

               if ( (power_a+power_b+power_c+power_d)==0 ) {
                 if ( plane==0 ) filZ=filZlong;
                 else filZ=filZycst; // necessarily, plane==2 (it cannot be hor. because of the condition above (a+c+(plane==1))%2 == 0)
               }
               else if ( (power_a+power_b+power_c+power_d)==1 ) {
                 if ( power_a==1 && plane==1 ) filZ=filZxdip;
                 else if ( power_b==1 && plane==2 ) filZ=filZydip;
                 else if ( power_c==1 && plane==1 ) filZ=filZxquad;
                 else if ( power_d==1 && plane==2 ) filZ=filZyquad;
                 else {
                   sprintf(Zoutput,"Z%c%d%d%d%d%s.dat",plane_char,power_a,power_b,power_c,power_d,output);
                   filZ=fopen(Zoutput,"w");
                   fprintf(filZ,"Frequency [Hz]\tRe(Z) [Ohm%s]\tIm(Z) [Ohm%s]\n",unit.c_str(),unit.c_str());
                 }
               }
               else {
                 sprintf(Zoutput,"Z%c%d%d%d%d%s.dat",plane_char,power_a,power_b,power_c,power_d,output);
                 filZ=fopen(Zoutput,"w");
                 fprintf(filZ,"Frequency [Hz]\tRe(Z) [Ohm%s]\tIm(Z) [Ohm%s]\n",unit.c_str(),unit.c_str());
               }

               // initialize impedance table
               for (unsigned int i=0; i<nf; i++) Z[i]=std::complex<double>(0.,0.);

               aplusc2=(power_a+power_c+(plane==1))/2;
               
               for (int r=0; r<=power_b/2; r++) {
                 for (int q=0; q<=aplusc2+(power_d+(plane==2))/2; q++) {
                   
                   m = power_b-2*r;
                   n = power_a+power_c+power_d+(plane>0)-2*q;
                   if (plane == 0) fact=1;
                   else {
                     if (plane == 1) fact=power_c+1;
                     else fact=power_d+1;
                   }
                   coef = fact*coefmn(power_a,power_b,power_c+(plane==1),power_d+(plane==2),m,n,r,q);
                   
                   // impedance computation at each frequency: beginning of the loop
                   for (unsigned int i=0; i<nf; i++) {
                        
                     memory.mem=0;  // initialize memory at each frequency
                     omega=amp::twopi<PRECISION>()*freq[i];
                     k=omega/(beta*C);
                     kovergamma=k/gamma;
                     if (plane==0) cst=4.*jimag*L*double(amp::ampf<PRECISION>(k*Z0/(beta*amp::sqr(gamma)*amp::twopi<PRECISION>())).toDouble());
                     else cst=4.*jimag*L*double(amp::ampf<PRECISION>(Z0/(beta*amp::sqr(gamma)*amp::twopi<PRECISION>())).toDouble());

                     // computes the layer properties for the angular freq. omega
                     for (unsigned int p=2;p<=N+1; p++) {
                       get_eps_and_mu_from_layer_properties(freq[i], rho(p), tau(p), epsb(p),
                          tandelta(p), chi(p), fmu(p), eps_frequencies[p], eps_values[p],
                          n_eps_frequencies[p], mu_frequencies[p], mu_values[p], n_mu_frequencies[p],
                          eps_vs_freq_kind[p], mu_vs_freq_kind[p], eps1(p), mu1(p));
                     
                       //printf("%lg Hz: layer %d -> %lg m, eps1= %lg + %lg j, mu1= %lg + %lgj\n",freq[i],p,double(b(p).toDouble()),
                       //   double(eps1(p).x.toDouble()),double(eps1(p).y.toDouble()),
                       //   double(mu1(p).x.toDouble()),double(mu1(p).y.toDouble()));
                       }
                     for (unsigned int p=2;p<=M+1; p++) {
                       get_eps_and_mu_from_layer_properties(freq[i], rhom(p), taum(p), epsbm(p),
                          tandeltam(p), chim(p), fmum(p), epsm_frequencies[p], epsm_values[p],
                          n_epsm_frequencies[p], mum_frequencies[p], mum_values[p], n_mum_frequencies[p],
                          epsm_vs_freq_kind[p], mum_vs_freq_kind[p], eps1m(p), mu1m(p));

                       //printf("%lg Hz: layer %d -> %lg m, eps1= %lg + %lg j, mu1= %lg + %lgj\n",freq[i],-p,double(bm(p).toDouble()),
                       //   double(eps1m(p).x.toDouble()),double(eps1m(p).y.toDouble()),
                       // double(mu1m(p).x.toDouble()),double(mu1m(p).y.toDouble()));
                       }

                     // computes alpha
                     alpha = alphamn(flag_topbotsym,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,m,n,limit,w,&memory);

                     Z[i] = Z[i] + alpha*cst*coef*std::pow(double(kovergamma.toDouble())/2.,power_a+power_b+power_c+power_d+(plane>0));
                   }
                 }
               }
               // print all values of Z
               for (unsigned int i=0; i<nf; i++) fprintf(filZ,"%13.8e %13.8e %13.8e\n",freq[i],Z[i].real(),Z[i].imag());
               fclose(filZ);
             }
             
             // in top/bottom symmetric case, still fill (with zeros) the Zycst term (for back compatibility)
             if ( flag_topbotsym && (power_a+power_b+power_c+power_d)==0 && plane==2 ) {
               for (unsigned int i=0; i<nf; i++) fprintf(filZycst,"%13.8e %13.8e %13.8e\n",freq[i],0.,0.);
               fclose(filZycst);
             }
           }
         }
       }
     }
   }
 }
 
 gsl_integration_workspace_free(w);
  
 // deallocate all frequency-dependent tables
 delete[] n_eps_frequencies;
 delete[] n_epsm_frequencies;
 delete[] n_mu_frequencies;
 delete[] n_mum_frequencies;

 for (unsigned int p=2; p<=N+1; p++){
   if ((eps_vs_freq_kind[p]==1) || (eps_vs_freq_kind[p]==2)) {
     delete[] eps_frequencies[p];
     delete[] eps_values[p];
   }
   if (mu_vs_freq_kind[p]==1) {
     delete[] mu_frequencies[p];
     delete[] mu_values[p];
   }
 }
 if (!flag_topbotsym) {
   for (unsigned int p=2; p<=M+1; p++){
     if ((epsm_vs_freq_kind[p]==1) || (epsm_vs_freq_kind[p]==2)) {
       delete[] epsm_frequencies[p];
       delete[] epsm_values[p];
     }
     if (mum_vs_freq_kind[p]==1) {
       delete[] mum_frequencies[p];
       delete[] mum_values[p];
     }
   }
 }
 
 delete[] eps_frequencies;
 delete[] eps_values;
 delete[] mu_frequencies;
 delete[] mu_values;
 delete[] epsm_frequencies;
 delete[] epsm_values;
 delete[] mum_frequencies;
 delete[] mum_values;

 delete[] eps_vs_freq_kind;
 delete[] epsm_vs_freq_kind;
 delete[] mu_vs_freq_kind;
 delete[] mum_vs_freq_kind;
 
 time(&end);
 dif=difftime(end,start);
 printf("Elapsed time during calculation: %.2lf seconds\n",dif);

}
