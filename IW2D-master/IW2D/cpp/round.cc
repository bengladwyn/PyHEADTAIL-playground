#include <iostream>
#include <fstream>
#include <stdio.h>

#include <complex>
#include <cmath>
#include <stdlib.h>
#include <amp.h>
#include <ablas.h>
#include <mpfr.h>

#include <IW2D.h>

/*
#define  Precision	160  // Number of bits for the mantissa of all numbers (classic doubles have 53 bits)
#define  MAXLINES	100  // Maximum number of lines in the input file (correspond to 5 layers in top and
			    // bottom parts)
#define  MAXCHAR	200  // Maximum number of characters in one line in the input file
#define  MAXCHARFILE	200  // Maximum number of characters in the output files name extension
*/


/**************************************************
*** multilayerm_round: computes matrix for field
*** matching (multiprecision) - round chamber
***************************************************/

  void multilayerm_round(ap::template_2d_array< amp::campf<PRECISION> >& mat,
  	unsigned int N, unsigned int m,
	ap::template_1d_array< amp::campf<PRECISION> >& eps1, 
  	ap::template_1d_array< amp::campf<PRECISION> >& mu1, 
	ap::template_1d_array< amp::ampf<PRECISION> >& b,
	amp::ampf<PRECISION> k, amp::ampf<PRECISION> beta){
	
    /* computes the matrix M for the multilayer field matching, from mu1, eps1, b of each of the N layers
    and the azimuthal mode number m, relativistic velocity factor beta
    and wave number k. */

    amp::campf<PRECISION> nup,nu2p,nup1,nu2p1,epsp1overnup1,epspovernup,mup1overnup1,mupovernup,xpp,xp1p;
    amp::campf<PRECISION> Impp,Imp1p,Kmpp,Kmp1p,Imppnorm,Imp1pnorm,Kmppnorm,Kmp1pnorm; // bessel functions of order m
    amp::campf<PRECISION> Im1pp,Im1p1p,Km1pp,Km1p1p,Im1ppnorm,Im1p1pnorm,Km1ppnorm,Km1p1pnorm; // bessel functions of order m-1
    amp::campf<PRECISION> Imppratio,Kmppratio,Imp1pratio,Kmp1pratio; // ratio of bessel functions
    amp::campf<PRECISION> ImIm,ImKm,KmIm,KmKm; // products of bessel functions
    amp::campf<PRECISION> tmp1;
    amp::ampf<PRECISION> mMP,beta2,k2;
    ap::template_2d_array< amp::campf<PRECISION> > mp1p;
    ap::template_2d_array< amp::campf<PRECISION> > matold;
   
    mp1p.setbounds(1,4,1,4);
    matold.setbounds(1,4,1,4);
    
    mMP=amp::ampf<PRECISION>(m);
  
    for (int i=1; i<=4; i++) {
      for (int j=1; j<=4; j++) {
        if (i==j) matold(i,j)=1;
	else matold(i,j)=0;
      }
    }
    
    mat=matold;

    beta2=amp::sqr(beta);
    k2=amp::sqr(k);
    // first layer
    nu2p=k2*(1-beta2*eps1(1)*mu1(1));
    set_minus0jimag_to_positive(nu2p);
    nup=csqrtMP(nu2p);
    epspovernup=eps1(1)/nup;
    mupovernup=mu1(1)/nup;

    for (unsigned int p=1; p<=N-1; p++) {

      nu2p1=k2*(1-beta2*eps1(p+1)*mu1(p+1));
      set_minus0jimag_to_positive(nu2p1);
      nup1=csqrtMP(nu2p1);
      epsp1overnup1=eps1(p+1)/nup1;
      mup1overnup1=mu1(p+1)/nup1;
      
      xpp=nup*b(p);
      xp1p=nup1*b(p);
      bessel(m,xpp,Impp,Kmpp);
      bessel(m,xp1p,Imp1p,Kmp1p);
      if (m==0) {
        bessel(1,xpp,Im1pp,Km1pp);
        bessel(1,xp1p,Im1p1p,Km1p1p);
	Imppratio = Im1pp/Impp;
	Kmppratio =-Km1pp/Kmpp;
	Imp1pratio= Im1p1p/Imp1p;
	Kmp1pratio=-Km1p1p/Kmp1p;
      } else {
        bessel(m-1,xpp,Im1pp,Km1pp);
        bessel(m-1,xp1p,Im1p1p,Km1p1p);
	Imppratio = Im1pp/Impp - mMP/xpp;
	Kmppratio =-Km1pp/Kmpp - mMP/xpp;
	Imp1pratio= Im1p1p/Imp1p - mMP/xp1p;
	Kmp1pratio=-Km1p1p/Kmp1p - mMP/xp1p;
      }
      
      ImIm=Impp*Imp1p;
      ImKm=Impp*Kmp1p;
      KmIm=Kmpp*Imp1p;
      KmKm=Kmpp*Kmp1p;
      
      // submatrix P
      tmp1=-xp1p/epsp1overnup1;
      mp1p(1,1)=tmp1*ImKm*(epsp1overnup1*Kmp1pratio-epspovernup*Imppratio);
      mp1p(1,2)=tmp1*KmKm*(epsp1overnup1*Kmp1pratio-epspovernup*Kmppratio);
      mp1p(2,1)=tmp1*ImIm*(-epsp1overnup1*Imp1pratio+epspovernup*Imppratio);
      mp1p(2,2)=tmp1*KmIm*(-epsp1overnup1*Imp1pratio+epspovernup*Kmppratio);
      
      // submatrix Q
      if (m==0) {
	mp1p(1,3)=0;
	mp1p(1,4)=0;
	mp1p(2,3)=0;
	mp1p(2,4)=0;
      } else {
        tmp1=-(nu2p1/nu2p-1)*mMP/(beta*eps1(p+1));
	mp1p(1,3)=-tmp1*ImKm;
	mp1p(1,4)=-tmp1*KmKm;
	mp1p(2,3)=tmp1*ImIm;
	mp1p(2,4)=tmp1*KmIm;
      }
      
      // submatrix R
      tmp1=-xp1p/mup1overnup1;
      mp1p(3,3)=tmp1*ImKm*(mup1overnup1*Kmp1pratio-mupovernup*Imppratio);
      mp1p(3,4)=tmp1*KmKm*(mup1overnup1*Kmp1pratio-mupovernup*Kmppratio);
      mp1p(4,3)=tmp1*ImIm*(-mup1overnup1*Imp1pratio+mupovernup*Imppratio);
      mp1p(4,4)=tmp1*KmIm*(-mup1overnup1*Imp1pratio+mupovernup*Kmppratio);
      
      // submatrix Q
      if (m==0) {
	mp1p(3,1)=0;
	mp1p(3,2)=0;
	mp1p(4,1)=0;
	mp1p(4,2)=0;
      } else {
	tmp1=eps1(p+1)/mu1(p+1);
	mp1p(3,1)=tmp1*mp1p(1,3);
	mp1p(3,2)=tmp1*mp1p(1,4);
	mp1p(4,1)=tmp1*mp1p(2,3);
	mp1p(4,2)=tmp1*mp1p(2,4);
      }
      
      // matrix multiplication
      ablas::cmatrixgemm<PRECISION>(4,4,4,1,mp1p,1,1,0,matold,1,1,0,0,mat,1,1);
      
      nu2p=nu2p1;
      nup=nup1;
      epspovernup=epsp1overnup1;
      mupovernup=mup1overnup1;

      matold=mat;
    }
    
  
    return;
  }
  
  
/**************************************************
*** alphaTM: computes alphaTM (multiprecision)
***
***************************************************/

  std::complex<double> alphaTM(unsigned int N, unsigned int m,
	ap::template_1d_array< amp::campf<PRECISION> >& eps1, 
  	ap::template_1d_array< amp::campf<PRECISION> >& mu1, 
	ap::template_1d_array< amp::ampf<PRECISION> >& b,
	amp::ampf<PRECISION> k, amp::ampf<PRECISION> beta){
	
    /* function that computes alphaTM for a given azimuthal mode number m, from mu1, eps1, b 
    of each of the N layers, and from the relativistic velocity factor beta and wave number k. */
    
    amp::campf<PRECISION> alphaTM;
    std::complex <double> result;
    ap::template_2d_array< amp::campf<PRECISION> > mat; // 4*4 field matching matrix
    //timeval c1,c2; // clock ticks


    /* setting bounds for matrix mat */
    mat.setbounds(1,4,1,4);

    // compute the field matching 4*4 matrices (upper and lower layers)
    //gettimeofday(&c1,0);
    multilayerm_round(mat, N, m, eps1, mu1, b, k, beta);
    //gettimeofday(&c2,0);
    //std::cout << "multilayer: time= " << (c2.tv_sec-c1.tv_sec)*1.e6+(c2.tv_usec-c1.tv_usec) << " usec\n";

    // compute alphaTM
    alphaTM=(mat(1,2)*mat(3,3)-mat(3,2)*mat(1,3))/(mat(1,1)*mat(3,3)-mat(1,3)*mat(3,1));
    
    // conversion to double complex
    result=std::complex<double>(amp::ampf<PRECISION>(alphaTM.x).toDouble(),amp::ampf<PRECISION>(alphaTM.y).toDouble());

    return result;
    
  }
  
