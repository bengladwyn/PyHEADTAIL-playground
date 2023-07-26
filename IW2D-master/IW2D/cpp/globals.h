#ifndef GLOBALS_H
#define GLOBALS_H

#include <amp.h>


#define  Precision	160  // Number of bits for the mantissa of all numbers (classic doubles have 53 bits)

ap::template_1d_array< amp::ampf<PRECISION> > kxmem;
ap::template_1d_array< amp::campf<PRECISION> > eta1mem,eta2mem,chi1mem,chi2mem;
// current number of elements of kxmem, eta1mem, etc.
unsigned long mem;


#endif
