/*
 *  roundchamber.cc 
 *  
 *  by Nicolas Mounet (Nicolas.Mounet@cern.ch)
 *
 *  computes the impedance in a round chamber (see CERN note by N. Mounet and E. Metral, 
 "Electromagnetic fields created by a macroparticle in an infinitely long
and axisymmetric multilayer beam pipe", 2009, CERN-BE-2009-039)
 
In input : typical input file is

Machine:	LHC
Relativistic Gamma:	479.6
Impedance Length in m:	2.8
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
start frequency exponent (10^) in Hz:	-1
stop frequency exponent (10^) in Hz:	15
linear (1) or logarithmic (0) or both (2) frequency scan:	2
sampling frequency exponent (10^) in Hz (for linear):	8
Number of points per decade (for log):	100
when both, fmin of the refinement (in THz):	0.0008
when both, fmax of the refinement (in THz):	0.1
when both, number of points in the refinement:	40000
added frequencies [Hz]:	1e-9 1e-8 1e-7 1e-6 1e-5 0.0001 0.001 0.01 1e+16 
Yokoya factors long, xdip, ydip, xquad, yquad:	1 1 1 0 0
Comments for the output files names:	_some_element

The order of the lines can be whatever, but the exact sentences and the TAB before the parameter
indicated, are necessary. If there are more layers than indicated by the number of layers,
the additional one are ignored. The last layer is always assumed to go to infinity.

In output one gives five files with the impedances (longitudinal, x dipolar, y dipolar,
x quadrupolar, y quadrupolar). Each have 3 columns :frequency, real part and
imaginary part.

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
#include <amp.h>
#include <ablas.h>
#include <mpfr.h>
#include <gsl/gsl_sort_double.h>

#include <IW2D.h>
#include <globals.h>

/*
#define  Precision	160  // Number of bits for the mantissa of all numbers (classic doubles have 53 bits)
#define  MAXLINES	100  // Maximum number of lines in the input file (correspond to 5 layers in top and
			    // bottom parts)
#define  MAXCHAR	200  // Maximum number of characters in one line in the input file
#define  MAXCHARFILE	200  // Maximum number of characters in the output files name extension
*/

/**************************************************
 *** 			main program		***
 ***						***
 **************************************************/

int main ()

{
 
 char *endline;
 char output[MAXCHARFILE],Zxdipoutput[MAXCHARFILE+10],Zydipoutput[MAXCHARFILE+10],Zlongoutput[MAXCHARFILE+10];
 char Zxquadoutput[MAXCHARFILE+10],Zyquadoutput[MAXCHARFILE+10],Input[MAXCHARFILE+10];
 std::string data[MAXLINES],machine,commentoutput,dummy2;
 FILE *filZxdip, *filZydip, *filZxquad, *filZyquad, *filZlong, *filInput;
 unsigned int N,dummy0; // number of layers, then dummy parameter
 unsigned int m,n_input,n_added,nf,nfdup; /* azimuthal mode number, number of input lines,
 					number of individually added frequencies, total number of
					frequencies in the scan, number of duplicate frequencies */
 unsigned int typescan,nflog,nflin; /* type of frequency scan, number of freq. per decade, number of freq. in a lin. scan inside the log scan*/
 ap::template_1d_array< amp::campf<PRECISION> > eps1,mu1; /* eps1 and mu1 for the layers */
 ap::template_1d_array< amp::ampf<PRECISION> > b,thick; // position of boundaries; thickness of the layers 
 ap::template_1d_array< amp::ampf<PRECISION> > rho,tau,epsb,chi,fmu; /* layers
 						properties (DC resistivity, resistivity relaxation time,
						dielectric constant, magnetic susceptibility=mu_r-1,
						relaxation frequency of permeability)*/
 amp::ampf<PRECISION> omega,beta,k,gamma,dummy3; // parameters
 amp::campf<PRECISION> jimagMP; // imaginary constant in multiprecision
 size_t limit=1000; // limit (number of intervals) for gsl integration algorithm
 double x,y,L,fminlog,fmaxlog,fminlin,fmaxlin,fsamplin,fadded[15],*freq,dif,dummy1,yokoya[5];
 std::complex<double> Zxdip,Zydip,Zxquad,Zyquad,Zlong,cst;
 std::complex<double> alphaTM0,alphaTM1,jimag; // alphaTM constants and imaginary constant
 size_t found;
 time_t start,end; // times
					
 // start time
 time(&start);
 
 // to test bessel functions
 /*amp::campf<Precision> zz,bes1,bes2,bes3,bes4;
 zz.x=1;zz.y=1;
 bessel(0,zz,bes1,bes2,bes3,bes4);bessel(0,zz*100000,bes1,bes2,bes3,bes4);
 bessel(2,zz,bes1,bes2,bes3,bes4);bessel(2,zz*100000,bes1,bes2,bes3,bes4);
 bessel(1,zz,bes1,bes2,bes3,bes4);bessel(1,zz*10,bes1,bes2,bes3,bes4);
 bessel(1,zz*100,bes1,bes2,bes3,bes4);
 zz.x="2.3";zz.y=0;
 bessel(0,zz,bes1,bes2,bes3,bes4);bessel(0,zz*100000,bes1,bes2,bes3,bes4);
 bessel(2,zz,bes1,bes2,bes3,bes4);bessel(2,zz*100000,bes1,bes2,bes3,bes4);
 bessel(1,zz,bes1,bes2,bes3,bes4);bessel(1,zz*10,bes1,bes2,bes3,bes4);
 bessel(1,zz*100,bes1,bes2,bes3,bes4);*/
 
 // constants
 jimagMP.x=0;jimagMP.y=1;
 jimag=std::complex<double>(0.,1.);
 
 // default values of the parameters (in case)
 typescan=0;
 fminlog=2;fmaxlog=13;nflog=10;n_added=0;
 fsamplin=8;fminlin=1;fmaxlin=2;nflin=100;nf=0;
 N=2;L=1.;gamma="479.6";
 
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
   read_input(data[i],"Number of layers",N,dummy1,dummy2,dummy3,0);
   //printf("yoyo %d %s %s %s %13.8e %d\n",i,data[i].c_str(),machine.c_str(),gamma.toDec().c_str(),L,N);
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
   read_input_yokoya_factors(data[i], yokoya);
 }

 //printf("%s %d \n",machine.c_str(),N);
 //printf("%13.8e %13.8e %d\n",fadded[1],fadded[2],n_added);
 //printf("%13.8e %13.8e %13.8e %13.8e %13.8e\n",yokoya[0],yokoya[1],yokoya[2],yokoya[3],yokoya[4]);
 //printf("%13.8e %13.8e %d %d\n",fminlog,fmaxlog,nflog,typescan);

 
 eps1.setbounds(1,N+1);mu1.setbounds(1,N+1);b.setbounds(1,N+1);thick.setbounds(1,N);
 rho.setbounds(1,N+1);tau.setbounds(1,N+1);epsb.setbounds(1,N+1);chi.setbounds(1,N+1);fmu.setbounds(1,N+1);

 // default values of the layers properties (in case)
 rho(2)="1e5";tau(2)=0;epsb(2)=1;chi(2)=0;fmu(2)="Infinity";b(1)=2;thick(1)="Infinity";

 // find inner radius(ii) of the chamber layers
 for (unsigned int i=1; i<=n_input; i++) {
   read_input(data[i],"Layer 1 inner radius in mm",dummy0,dummy1,dummy2,b(1),3);
 }
 //printf("%s \n",b(1).toDec().c_str());
 
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

 // units conversion to SI (fminlin and fmaxlin were in THz, tau was in ps, b in mm and fmu in MHz)
 fminlin*=1.e12;fmaxlin*=1.e12;
 b(1)=b(1)*1e-3;
 for (unsigned int p=2; p<=N+1; p++){
   tau(p)=tau(p)*1e-12;
   b(p)=thick(p-1)*1e-3+b(p-1);
   fmu(p)=fmu(p)*1e6;
   //printf("%s %s %s %s %s %s \n",rho(p).toDec().c_str(),tau(p).toDec().c_str(),epsb(p).toDec().c_str(),
   //	chi(p).toDec().c_str(),fmu(p).toDec().c_str(),b(p).toDec().c_str());
 }
  
 // first layer (inside the chamber) is always vacuum
 eps1(1)=1;mu1(1)=1;
 // relativistic velocity factor beta
 beta=amp::sqrt(1-1/amp::sqr(gamma));
 //printf("%s %s\n",eps1(1).x.toDec().c_str(),eps1(1).y.toDec().c_str());
 
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
 sprintf(output,"W%s_%dlayers%.2lfmm%s",machine.c_str(),N,
	1e3*double(amp::ampf<PRECISION>(b(1)).toDouble()),commentoutput.c_str());
 sprintf(Zxdipoutput,"Zxdip%s.dat",output);
 sprintf(Zydipoutput,"Zydip%s.dat",output);
 sprintf(Zxquadoutput,"Zxquad%s.dat",output);
 sprintf(Zyquadoutput,"Zyquad%s.dat",output);
 sprintf(Zlongoutput,"Zlong%s.dat",output);
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
 fprintf(filZxdip,"Frequency [Hz]\tRe(Zxdip) [Ohm/m]\tIm(Zxdip) [Ohm/m]\n");
 fprintf(filZydip,"Frequency [Hz]\tRe(Zydip) [Ohm/m]\tIm(Zydip) [Ohm/m]\n");
 fprintf(filZxquad,"Frequency [Hz]\tRe(Zxquad) [Ohm/m]\tIm(Zxquad) [Ohm/m]\n");
 fprintf(filZyquad,"Frequency [Hz]\tRe(Zyquad) [Ohm/m]\tIm(Zyquad) [Ohm/m]\n");
 fprintf(filZlong,"Frequency [Hz]\tRe(Zlong) [Ohm]\tIm(Zlong) [Ohm]\n");
 
 
 // impedance computation at each frequency: beginning of the loop
 for (unsigned int i=0; i<nf; i++) {
    
   omega=amp::twopi<PRECISION>()*freq[i];
   k=omega/(beta*C);
   
   // computes the layer properties for the angular freq. omega
   for (unsigned int p=2;p<=N+1; p++) {
     if (rho(p).isFiniteNumber()) {
       eps1(p)=epsb(p)+1/(jimagMP*eps0*rho(p)*omega*(1+jimagMP*omega*tau(p)));
     } else {
       eps1(p)=epsb(p);
     }
     mu1(p)=1+chi(p)/(1+jimagMP*omega/(amp::twopi<PRECISION>()*fmu(p)));
     //printf("%s %s %s\n%s %s\n",b(p).toDec().c_str(),eps1(p).x.toDec().c_str(),eps1(p).y.toDec().c_str(),
     //		mu1(p).x.toDec().c_str(),mu1(p).y.toDec().c_str());
     }
   //printf("%s \n",b(1).toDec().c_str());
   //printf("%s %s \n",omega.toDec().c_str(),k.toDec().c_str());
 
   // computes alphaTM0
   alphaTM0=alphaTM(N+1,0,eps1,mu1,b,k,beta);
   //printf("alphaTM0: %13.8e %13.8e\n",alphaTM0.real(),alphaTM0.imag());
 
   // computes alphaTM1
   alphaTM1=alphaTM(N+1,1,eps1,mu1,b,k,beta);
   //printf("alphaTM1: %13.8e %13.8e\n",alphaTM1.real(),alphaTM1.imag());
 
   
   // computes and writes the impedances (applying yokoya factors)
   cst=jimag*L*double(amp::ampf<PRECISION>(k*Z0/(beta*amp::sqr(gamma)*amp::twopi<PRECISION>())).toDouble());
   Zlong=cst*alphaTM0*yokoya[0];
   cst=cst*double(amp::ampf<PRECISION>(k/(2*amp::sqr(gamma))).toDouble());
   Zxdip=cst*alphaTM1*yokoya[1];
   Zydip=cst*alphaTM1*yokoya[2];
   if ( (yokoya[0]==1)&& ( (yokoya[1]==1)&&(yokoya[2]==1) ) ) {
     // case of axisymmetric geometry -> small quadrupolar impedance (see ref. cited at the beginning)
     Zxquad=cst*alphaTM0;
     Zyquad=Zxquad;
   } else {     
     Zxquad=cst*alphaTM1*yokoya[3];
     Zyquad=cst*alphaTM1*yokoya[4];
   }
   fprintf(filZlong,"%13.8e %13.8e %13.8e\n",freq[i],Zlong.real(),Zlong.imag());
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
 
 time(&end);
 dif=difftime(end,start);
 printf("Elapsed time during calculation: %.2lf seconds\n",dif);

}
