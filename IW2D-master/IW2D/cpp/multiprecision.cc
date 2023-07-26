//#include <matinv.h>
//#include <densesolver.h>
//#include <trfac.h>
#include <amp.h>
#include <mpfr.h>
#include <iostream>
#include "arb.h"
#include "acb.h"
#include "arf.h"
#include "acb_hypgeom.h"

#include <IW2D.h>


/**************************************************
*** acb_from_campf: Set the midpoint value of res_acb to be the same value as z_campf.
*** Sets radius of res_acb to 0.
*** Should produce deep copy, i.e. no two pointers to same memory afterwards.
***************************************************/

void acb_from_campf(acb_t res_acb, const amp::campf<PRECISION>& z_campf) {
    acb_zero(res_acb);
    arf_set_mpfr(arb_midref(acb_realref(res_acb)), z_campf.x.getReadPtr());
    arf_set_mpfr(arb_midref(acb_imagref(res_acb)), z_campf.y.getReadPtr());
}


/**************************************************
*** campf_from_acb: set res_campf to the midpoint value of z_acb.
*** Radius of z_acb is ignored.
*** Should produce deep copy, i.e. no two pointers to same memory afterwards.
***************************************************/

void campf_from_acb(amp::campf<PRECISION>& res_campf, acb_t z_acb) {
    
    // From AMP docs: "A pointer which has been received by getWritePtr()
    // can be used in operations which can change its contents,
    // but these operations should not change the precision of a pointer contents."
    // Therefore necessary to (re)setting the precision of the acb to PRECISION with acb_set_round.
    acb_set_round(z_acb, z_acb, PRECISION);
    arf_get_mpfr(res_campf.x.getWritePtr(), arb_midref(acb_realref(z_acb)), MPFR_RNDN);
    arf_get_mpfr(res_campf.y.getWritePtr(), arb_midref(acb_imagref(z_acb)), MPFR_RNDN);
}

/**************************************************
*** bessel: Compute the modified bessel functions I and K of order m and 
*** complex input z, and set the arguments besi and besk to the results, respectively for I and K.
***
***************************************************/

void bessel(unsigned int m, const amp::campf<PRECISION>& z, amp::campf<PRECISION>& besi, amp::campf<PRECISION>& besk) {
    acb_t z_acb;
    acb_init(z_acb);
    acb_from_campf(z_acb, z);

    acb_t nu;
    acb_init(nu);
    acb_set_ui(nu, m);

    acb_t res_besi;
    acb_init(res_besi);
    acb_t res_besk;
    acb_init(res_besk);

    // Calculate bessel functions with doubling precision until ouput is accurate to PRECISION bits
    // Breaks loop if BESSEL_PRECISION_SCALING_LIMIT is reached
    unsigned int precision_scaling = 2;
    acb_hypgeom_bessel_i(res_besi, nu, z_acb, precision_scaling*PRECISION);
    acb_hypgeom_bessel_k(res_besk, nu, z_acb, precision_scaling*PRECISION);

    while (
      (acb_rel_accuracy_bits(res_besi) < PRECISION || acb_rel_accuracy_bits(res_besk) < PRECISION)
       && precision_scaling <= BESSEL_PRECISION_SCALING_LIMIT) // 
    {
      precision_scaling *= 2;
      acb_hypgeom_bessel_i(res_besi, nu, z_acb, precision_scaling*PRECISION);
      acb_hypgeom_bessel_k(res_besk, nu, z_acb, precision_scaling*PRECISION);
    } ;

    if (precision_scaling == BESSEL_PRECISION_SCALING_LIMIT) {
      printf("Warning: Terminated bessel calculation after not reaching convergence using %i * %i precision bits.\n", precision_scaling, PRECISION);
      printf("Input z     =  "); acb_printd(z_acb, 15); printf(" (precision truncated)\n");
      printf("Input m     =  "); acb_printd(nu, 15); printf(" (precision truncated)\n");
      printf("Output I(z) =  "); acb_printd(res_besi, 15); printf(" (precision truncated)\n");
      printf("Output K(z) =  "); acb_printd(res_besk, 15); printf(" (precision truncated)\n");

    }

    campf_from_acb(besi, res_besi);
    campf_from_acb(besk, res_besk);
    acb_clear(z_acb);
    acb_clear(nu);
    
    acb_clear(res_besi);
    acb_clear(res_besk);
}


/* to use cmatrixgemm (ALGLIB routine -- 16.12.2009 Bochkanov Sergey)
INPUT PARAMETERS
  M - matrix size, M>0
  N - matrix size, N>0
  K - matrix size, K>0
  Alpha - coefficient 
  A - matrix 
  IA - submatrix offset 
  JA - submatrix offset 
  OpTypeA - transformation type: * 0 - no transformation * 1 - transposition * 2 - conjugate transposition 
  B - matrix 
  IB - submatrix offset 
  JB - submatrix offset 
  OpTypeB - transformation type: * 0 - no transformation * 1 - transposition * 2 - conjugate transposition 
  Beta - coefficient 
  C - matrix 
  IC - submatrix offset 
  JC - submatrix offset 

template<unsigned int Precision> void cmatrixgemm(int m, int n, int k, 
	amp::campf<Precision> alpha, const ap::template_2d_array< amp::campf<Precision> >& a,
	int ia, int ja, int optypea, const ap::template_2d_array< amp::campf<Precision> >& b, 
	int ib, int jb, int optypeb, amp::campf<Precision> beta, 
	ap::template_2d_array< amp::campf<Precision> >& c, int ic, int jc);
*/

  void set_minus0jimag_to_positive(amp::campf<PRECISION>& z) {
    if (z.y.toDouble()==0. && std::signbit(z.y.toDouble())) z.y = -z.y;
  }

/**************************************************
*** csqrtMP: Complex square root in multiprecision
***
***************************************************/

  amp::campf<PRECISION> csqrtMP(amp::campf<PRECISION> z){
  
    /* Complex square root of z in multiprecision. Convention specified in CERN note mentioned at the
    beginning. */
  
    amp::campf<PRECISION> result;
    amp::ampf<PRECISION> rho,phi;

    // if (z.y.toDouble() == 0.0 && std::signbit(z.y.toDouble())) printf("%lf\n",z.y.toDouble());
    
    rho=amp::sqrt(amp::abscomplex(z));
    phi=amp::atan2(z.y,z.x)/2;
    
    result.x=rho*amp::cos(phi);
    result.y=rho*amp::sin(phi);
    
    return result;
  
  }
  

/**************************************************
*** cexpMP: Complex exponential in multiprecision
***
***************************************************/

  amp::campf<PRECISION> cexpMP(amp::campf<PRECISION> z){
  
    /* Complex exponential of z in multiprecision. */
  
    amp::campf<PRECISION> result;
    amp::ampf<PRECISION> rho,phi;
    
    rho=amp::exp(z.x);
    phi=z.y;
    
    result.x=rho*amp::cos(phi);
    result.y=rho*amp::sin(phi);
    
    return result;
  
  }
  
  
/**************************************************
*** clogMP: Complex natural logarithm in multiprecision
***
***************************************************/

  amp::campf<PRECISION> clogMP(amp::campf<PRECISION> z){
  
    /* Complex natural logarithm of z in multiprecision. */
  
    amp::campf<PRECISION> result;
    amp::ampf<PRECISION> rho,phi;
    
    rho=amp::log(amp::abscomplex(z));
    phi=amp::atan2(z.y,z.x);
    
    result.x=rho;
    result.y=phi;
    
    return result;
  
  }


/**************************************************
*** matinv4: computes the 2 first lines of coefficients
*** of the inverse of a 4x4 matrix (multiprecision).
*** This is optimized.
***************************************************/

  void matinv4(ap::template_2d_array< amp::campf<PRECISION> >& mat) {
  
    /* input= matrix mat (bounds 0..3,0..3). output is also mat, with the first two lines
       filled with the coefficients of its inverse */
    
    // We need the two first columns of the cofactor matrix, and the total determinant
    ap::template_2d_array< amp::campf<PRECISION> > cof; // cofactor matrix (two first columns only)
    ap::template_1d_array< amp::campf<PRECISION> > det2; // array with the 2x2 determinants of the 2 last columns
    amp::campf<PRECISION> det4; //determinant of the initial 4x4 matrix mat
    
    cof.setbounds(0,3,0,1);
    det2.setbounds(0,5);
    
    // det2 is ordered from top to bottom
    det2(0)=mat(0,2)*mat(1,3)-mat(0,3)*mat(1,2);
    det2(1)=mat(0,2)*mat(2,3)-mat(0,3)*mat(2,2);
    det2(2)=mat(0,2)*mat(3,3)-mat(0,3)*mat(3,2);
    det2(3)=mat(1,2)*mat(2,3)-mat(1,3)*mat(2,2);
    det2(4)=mat(1,2)*mat(3,3)-mat(1,3)*mat(3,2);
    det2(5)=mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2);
    
    // cofactors of the first column
    cof(0,0)=  mat(1,1)*det2(5)-mat(2,1)*det2(4)+mat(3,1)*det2(3);
    cof(1,0)=-(mat(0,1)*det2(5)-mat(2,1)*det2(2)+mat(3,1)*det2(1));
    cof(2,0)=  mat(0,1)*det2(4)-mat(1,1)*det2(2)+mat(3,1)*det2(0);
    cof(3,0)=-(mat(0,1)*det2(3)-mat(1,1)*det2(1)+mat(2,1)*det2(0));
    
    // total determinant
    det4=mat(0,0)*cof(0,0)+mat(1,0)*cof(1,0)+mat(2,0)*cof(2,0)+mat(3,0)*cof(3,0);
    
    // cofactors of the second column
    cof(0,1)=-(mat(1,0)*det2(5)-mat(2,0)*det2(4)+mat(3,0)*det2(3));
    cof(1,1)=  mat(0,0)*det2(5)-mat(2,0)*det2(2)+mat(3,0)*det2(1);
    cof(2,1)=-(mat(0,0)*det2(4)-mat(1,0)*det2(2)+mat(3,0)*det2(0));
    cof(3,1)=  mat(0,0)*det2(3)-mat(1,0)*det2(1)+mat(2,0)*det2(0);
    
    // final coefficients sought for (we transpose the cofactors and divide by determinant)
    mat(0,0)=cof(0,0)/det4;
    mat(0,1)=cof(1,0)/det4;
    mat(0,2)=cof(2,0)/det4;
    mat(0,3)=cof(3,0)/det4;
    mat(1,0)=cof(0,1)/det4;
    mat(1,1)=cof(1,1)/det4;
    mat(1,2)=cof(2,1)/det4;
    mat(1,3)=cof(3,1)/det4;
    
    
  }
