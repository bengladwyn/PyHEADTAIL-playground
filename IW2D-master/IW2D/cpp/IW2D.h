#ifndef IW2D_H
#define IW2D_H

#include <complex>
#include <amp.h>
#include <mpfr.h>
#include <gsl/gsl_integration.h>


#define  PRECISION	160  // Number of bits for the mantissa of all numbers (classic doubles have 53 bits)
#define  MAXLINES	100  // Maximum number of lines in the input file (correspond to 5 layers in top and
                         // bottom parts)
#define  MAXCHAR	200  // Maximum number of characters in one line in the input file
#define  MAXCHARFILE	200  // Maximum number of characters in the output files name extension


const unsigned int precision = PRECISION; // Precision as a variable, needed when binding to python

const unsigned int BESSEL_PRECISION_SCALING_LIMIT = 16; // Max precision used in adaptive bessel function calculation

// struct containing lookup tables to be passed between functions
struct memorycontainer{
  ap::template_1d_array< amp::ampf<PRECISION> > kxmem;
  ap::template_1d_array< amp::campf<PRECISION> > eta1mem;
  ap::template_1d_array< amp::campf<PRECISION> > eta2mem;
  ap::template_1d_array< amp::campf<PRECISION> > chi1mem;
  ap::template_1d_array< amp::campf<PRECISION> > chi2mem;
  unsigned long mem;
  unsigned long maxmem; // not strictly needed, since the arrays have a gethighbound() method, but this allows more readable code
  bool has_printed_warning; // used to prevent multiple warnings being printed for full memory between each memory reset

  memorycontainer(unsigned long maxmem)
  : mem{0}, maxmem{maxmem}, has_printed_warning{false}
  {
	kxmem.setbounds(0, maxmem);
	eta1mem.setbounds(0, maxmem);
	eta2mem.setbounds(0, maxmem);
	chi1mem.setbounds(0, maxmem);
	chi2mem.setbounds(0, maxmem);
  }
};

extern "C"
{

  const  amp::ampf<PRECISION> C=299792458;    // velocity of light [m/s]
  const amp::ampf<PRECISION> mu0=amp::ampf<PRECISION>("1.25663706212e-6"); // mu0
  const amp::ampf<PRECISION> eps0=1/(mu0*amp::sqr(C)); // eps0
  const amp::ampf<PRECISION> Z0=mu0*C; // free space impedance
  const std::complex<double> jimag=std::complex<double>(0.,1.); // imaginary constant
  const  amp::ampf<PRECISION> euler="0.577215664901532860606512090082402431042159335939923598805767234884867726777664670936947063"; // Euler's gamma constant
  const long double pi = 3.141592653589793238462643383279502884197;
  
  struct params;

  unsigned long locate (double *table, double z, unsigned long n);

  unsigned long locateMP (ap::template_1d_array< amp::ampf<PRECISION> >& table, amp::ampf<PRECISION> z, unsigned long n);
  
  std::complex<double> pchip(double z, double *zi, std::complex<double> *fi, std::complex<double> *di, unsigned long nz);
	
  std::complex<double> interp(double z, double *zi, std::complex<double> *fi, std::complex<double> *di, unsigned long nz, unsigned int interp_type);
	
  void pchip_derivatives(std::complex<double>* d, double* x, std::complex<double> * y, unsigned long length);

  std::complex<long double> Phi(std::complex<long double> x,double eps);
    
  std::complex<long double> Psi(std::complex<long double> x,double eps);
    
  std::complex<long double> Lambda(std::complex<long double> x,double eps);
    
  std::complex<double> fourier_integral_inf(std::complex<double>* fi,std::complex<double>* d,std::complex<double> df, long double t, long double* omegai, long double* delta, unsigned long length, double eps, unsigned int* interp_type, int flaginf);

  void set_minus0jimag_to_positive(amp::campf<PRECISION>& z);

  amp::campf<PRECISION> csqrtMP(amp::campf<PRECISION> z);
  
  amp::campf<PRECISION> cexpMP(amp::campf<PRECISION> z);
  
  amp::campf<PRECISION> clogMP(amp::campf<PRECISION> z);
  
void bessel(unsigned int m, const amp::campf<PRECISION>& z, amp::campf<PRECISION>& besi, amp::campf<PRECISION>& besk);
  
  void matinv4(ap::template_2d_array< amp::campf<PRECISION> >& mat);
  
  void read_input(std::string line, std::string description, unsigned int& param0, double& param1, std::string& param2, amp::ampf<PRECISION>& param3, int type);
  
  void read_input_layer(std::string line, std::string description, int p, amp::ampf<PRECISION>& param);
  
  void read_input_layer_filename(std::string line, std::string description, int p, std::string& param);
  
  void read_input_double_list(std::string line, std::string description, double *output_arr, unsigned int& n_added);

  void read_input_yokoya_factors(std::string line, double *output_arr);

  void read_complex_file_length(std::string filename, long& n);
  
  void read_complex_file(std::string filename,double *& frequencies, std::complex<double> *& values, long n);
  
  void read_real_file_length(std::string filename, long& n);
  
  void read_real_file(std::string filename,double *frequencies, std::complex<double> *values, long n);
  
  amp::campf<PRECISION> get_value_from_interp(double frequency, double *frequencies, std::complex<double> *values, long n);
  
  void get_eps_and_mu_from_layer_properties(double frequency,
    amp::ampf<PRECISION> rho, amp::ampf<PRECISION> tau, 
    amp::ampf<PRECISION> epsb, amp::ampf<PRECISION> tandelta, 
    amp::ampf<PRECISION> chi, amp::ampf<PRECISION> fmu,
    double *eps_frequencies, std::complex<double> *eps_values, long n_eps_frequencies, 
    double *mu_frequencies, std::complex<double> *mu_values, long n_mu_frequencies,
    int eps_vs_freq_kind, int mu_vs_freq_kind,
    amp::campf<PRECISION>& eps1, amp::campf<PRECISION>& mu1);
    
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
	amp::ampf<PRECISION> kx, amp::ampf<PRECISION> beta);

  void multilayerm_round(ap::template_2d_array< amp::campf<PRECISION> >& mat,
  	unsigned int N, unsigned int m,
	ap::template_1d_array< amp::campf<PRECISION> >& eps1, 
  	ap::template_1d_array< amp::campf<PRECISION> >& mu1, 
	ap::template_1d_array< amp::ampf<PRECISION> >& b,
	amp::ampf<PRECISION> k, amp::ampf<PRECISION> beta);
  
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
	amp::ampf<PRECISION> kx, amp::ampf<PRECISION> beta);
   
   std::complex<double> alphaTM(unsigned int N, unsigned int m,
	ap::template_1d_array< amp::campf<PRECISION> >& eps1, 
  	ap::template_1d_array< amp::campf<PRECISION> >& mu1, 
	ap::template_1d_array< amp::ampf<PRECISION> >& b,
	amp::ampf<PRECISION> k, amp::ampf<PRECISION> beta);
  
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
	amp::ampf<PRECISION> kovergamma);
    
  double integrand_real(double x, void *p);
  
  double integrand_imag(double x, void *p);
  
  double integrand_real_modif(double t, void *p);
  
  double integrand_imag_modif(double t, void *p);
  
  double integrate(int flagreal, unsigned int M, unsigned int N, 
  	ap::template_1d_array< amp::ampf<PRECISION> > b, 
	ap::template_1d_array< amp::ampf<PRECISION> > bm, 
	amp::ampf<PRECISION> beta, ap::template_1d_array< amp::campf<PRECISION> > eps1,
	ap::template_1d_array< amp::campf<PRECISION> > eps1m,
	ap::template_1d_array< amp::campf<PRECISION> > mu1,
	ap::template_1d_array< amp::campf<PRECISION> > mu1m,
	amp::ampf<PRECISION> omega, amp::ampf<PRECISION> k, amp::ampf<PRECISION> kovergamma,
	unsigned int m, unsigned int n, size_t limit, gsl_integration_workspace *w,
	memorycontainer* memory);

  std::complex<double> alphamn(unsigned int flag_topbotsym, unsigned int M, unsigned int N, 
  	ap::template_1d_array< amp::ampf<PRECISION> > b, 
	ap::template_1d_array< amp::ampf<PRECISION> > bm, 
	amp::ampf<PRECISION> beta, ap::template_1d_array< amp::campf<PRECISION> > eps1,
	ap::template_1d_array< amp::campf<PRECISION> > eps1m,
	ap::template_1d_array< amp::campf<PRECISION> > mu1,
	ap::template_1d_array< amp::campf<PRECISION> > mu1m,
	amp::ampf<PRECISION> omega, amp::ampf<PRECISION> k, amp::ampf<PRECISION> kovergamma,
	unsigned int m, unsigned int n, size_t limit, gsl_integration_workspace *w,
	memorycontainer* memory);

  long factorial(int n);

  int minus1pow(int n);

  double coefmn(int a, int b, int c, int d, int m, int n, int p, int q);

}

#endif
