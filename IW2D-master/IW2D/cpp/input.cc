#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <sys/time.h>

#include <complex>
#include <cmath>
#include <stdlib.h>
#include <amp.h>
#include <mpfr.h>

#include <IW2D.h>


/**************************************************
*** read_input: read an input file line
***
***************************************************/

  void read_input(std::string line, std::string description, unsigned int& param0, double& param1, 
  	std::string& param2, amp::ampf<PRECISION>& param3, int type){
	
    /* function reading the number or string at the end of an input line, identifying the description string.
    type=0 if output is an unsigned int, 1 if output is a double, 2 if output is a char array, and 3 if output is an
    amp high precision number.
    Output is in "param[type_number]". */
    
    size_t found;

    if (line.find(description) != std::string::npos){
      found=line.find("\t");

      if (found != std::string::npos){
         
        switch(type){
	  
	  case 0:
    
	    param0=atoi(line.substr(found+1,std::string::npos).c_str());
	    break;
	    
	  case 1:
	  
	    param1=strtod(line.substr(found+1,std::string::npos).c_str(),NULL);
	    break;
	    
	  case 2:
	  
	    param2=line.substr(found+1,std::string::npos);
	    break;
	    
	  case 3:
	  
	    param3=line.substr(found+1,std::string::npos).c_str();
	    break;
	
	}

      }

    }
    
    return;

  }

/**************************************************
*** read_input_layer: read an input file line
*** for a specific "layer property" line
***************************************************/

  void read_input_layer(std::string line, std::string description, int p, amp::ampf<PRECISION>& param){
	
    /* function reading the number or string at the end of an input line for a line containing a property of layer p.
    Output is an amp high precision number, in param. */
    
    unsigned int dummy0;
    double dummy1;
    std::string dummy2;
    
    std::ostringstream des;  
    des << "Layer " << p << " " << description;
    //std::cout << des.str() << "\n";
    //std::cout << line << "\n";
    read_input(line,des.str(),dummy0,dummy1,dummy2,param,3);
    
    return;
    
  }

/**************************************************
*** read_input_layer_filename: read an input file 
*** line for a specific "layer property" line, 
*** containing the name of a file containing a 
*** frequency dependent quantity (e.g. conductivity)
***************************************************/

  void read_input_layer_filename(std::string line, std::string description, int p, std::string& param){
	
    /* function reading the string at the end of an input line for a line containing a property of layer p.
    Output is a string, in param. */
    
    unsigned int dummy0;
    double dummy1;
    amp::ampf<PRECISION> dummy2;
    
    std::ostringstream des;  
    des << "Layer " << p << " " << description;
    //std::cout << des.str() << "\n";
    //std::cout << line << "\n";
    read_input(line,des.str(),dummy0,dummy1,param,dummy2,2);
    
    return;
    
  }

/**************************************************
*** read_input_double_list: Read an input file line
*** which contains a whitespace separated list of
*** arbitrarily many doubles.
***************************************************/
  void read_input_double_list(std::string line, std::string description, double *output_arr, unsigned int& n_added) {
    /*Args:
    line - Any line from the input file. If the start of the line matches description, the double list is read from after the \t character
    description - Line description to match against
    output_arr - array in which to store the values
    n_added - is set to the number of elements added to the array
    */

    // Check if line start matches description
    if (line.find(description) == std::string::npos) return;

    // Check if \t in line, find index
    size_t found = line.find("\t");
    if (found == std::string::npos) {
      printf("Warning: Could not find \\t character in input file, line:\n  %s", line.c_str());
      return;
    }

    n_added = 0;
    std::istringstream input_stream{line.substr(found+1, std::string::npos)};

    // Add all whitespace separated from line after \t character to output_arr
    // Ends on problems or end-of-line
    while (input_stream >> output_arr[n_added]) {
      n_added++;
    }

    // If end of line was not reached, print a warning
    if (!input_stream.eof()) {
      printf("Warning: Problem reading element %i (0-indexed) of input file line:\n  %s\n", n_added, line.c_str());
    }
  }

/**************************************************
*** read_input_yokoya_factors: Read an input file line
*** which contains a whitespace separated list of
*** yokoya factors. Throws an exception if it doesn't
*** find exactly 5 doubles in a line with the correct
*** description at the start.
***************************************************/
  void read_input_yokoya_factors(std::string line, double *output_arr) {
    /*
    Args:
    line - Any line from the input file. If the start of the line matches the Yokoya description, 5 doubles is read from after the \t character
    output_arr - array in which to store the values
    */

    // Check if line start matches Yokoya description
    if (line.find("Yokoya factors long, xdip, ydip, xquad, yquad") == std::string::npos) return;

    // Check if \t in line, find index
    size_t found = line.find("\t");
    if (found == std::string::npos) {
      printf("Warning: Could not find \\t character in input file, line:\n  %s", line.c_str());
      return;
    }
    
    // Read the yokoya factors, throw error if not exactly 5
    std::istringstream input_stream{line.substr(found+1,std::string::npos)};
    for (int k=0; k<5; k++) {
      if (!(input_stream >> output_arr[k])) {
        printf("Error: Problem reading Yokoya factor %i (0-indexed) from input file, line:\n  %s\n", k, line.c_str());
        throw "Yokoya error";
      }
    }
  }

/**************************************************
*** read_complex_file_length: find the length of a 
*** file containing complex numbers as a function 
*** of frequency
***************************************************/

  void read_complex_file_length(std::string filename, long& n){
    
    /* function getting the length of a file with freq. dependent complex numbers
    (3-column file). Note: first line ignored (headers). */
    
    FILE *f;
    double dummy1,dummy2,dummy3;
    char dummy_string[MAXCHAR];
    int condition;
    long j=0;
 
    // determines the number of lines n
    f = fopen(filename.c_str(),"r");
    fgets (dummy_string,MAXCHAR,f); // skip first line
    if (f == NULL)
      { printf ("\n\n    File %s does not exist.\n",filename.c_str()) ;
        exit(2); }
    n = 0;
    do {condition = fscanf (f, "%lg %lg %lg\n",&dummy1,&dummy2,&dummy3);
      ++n;
    } while (condition != EOF);
    n--;
    fclose(f);
    
  return;
  }

/**************************************************
*** read_complex_file: read a file containing complex 
*** numbers as a function of frequency (to be 
*** interpolated later), knowing its length
***************************************************/

  void read_complex_file(std::string filename,double *& frequencies, 
    std::complex<double> *& values, long n){
    
    /* function reading a file with freq. dependent complex numbers
    (3-column file), and put them in tables. Length of file is given by n.
    Note: first line ignored (headers). */
    
    FILE *f;
    double dummy1,dummy2,dummy3;
    char dummy_string[MAXCHAR];
    int condition;
    long j=0;
 
    frequencies = new double[n+1];
    values = new std::complex<double>[n+1];
    
    f = fopen(filename.c_str(),"r");
    fgets (dummy_string,MAXCHAR,f); // skip first line
    do {condition = fscanf (f, "%lg %lg %lg\n",&frequencies[j],&dummy2,&dummy3);
      values[j] = std::complex<double>(dummy2,dummy3);
 	  if ( (j>0) && (frequencies[j] == frequencies[j-1]) ) {
	    // remove duplicate frequencies
	    j--;
	    n--;
	  }
 	  j++;
     } while (condition != EOF);
     fclose(f);
 
  //for (j=0;j<=n-1;j++) printf("%ld %ld %13.8e %13.8e %13.8e\n",n,j,frequencies[j],values[j].real(),values[j].imag());
 
  return;
  }

/**************************************************
*** read_real_file_length: find the length of a 
*** file containing real numbers as a function 
*** of frequency
***************************************************/

  void read_real_file_length(std::string filename, long& n){
    
    /* function getting the length of a file with freq. dependent numbers
    (2-column file). Note: first line ignored (headers). */
    
    FILE *f;
    double dummy1,dummy2;
    char dummy_string[MAXCHAR];
    int condition;
    long j=0;
 
    // determines the number of lines n
    f = fopen(filename.c_str(),"r");
    fgets (dummy_string,MAXCHAR,f); // skip first line
    if (f == NULL)
      { printf ("\n\n    File %s does not exist.\n",filename.c_str()) ;
        exit(2); }
    n = 0;
    do {condition = fscanf (f, "%lg %lg\n",&dummy1,&dummy2);
      ++n;	 
    } while (condition != EOF);
    n--;
    fclose(f);

  return;
  }

/**************************************************
*** read_real_file: read a file containing real 
*** numbers as a function of frequency (to be 
*** interpolated later), knowing its length
***************************************************/

  void read_real_file(std::string filename,double *frequencies, 
    std::complex<double> *values, long n){
    
    /* function reading a file with freq. dependent numbers
    (2-column file), and put them in tables. Length of file is given by n.
    Note: first line ignored (headers). */
    
    FILE *f;
    double dummy1,dummy2;
    char dummy_string[MAXCHAR];
    int condition;
    long j=0;
 
    //values = new std::complex<double>[n];
    //frequencies = new double[n];
    
    f = fopen(filename.c_str(),"r");
    fgets (dummy_string,MAXCHAR,f); // skip first line
    do {condition = fscanf (f, "%lg %lg\n",&frequencies[j],&dummy2);
      values[j] = std::complex<double>(dummy2,0.);
      //printf("%s %ld %lg %lg %lg %lg %d %d\n",filename.c_str(),j,dummy2,frequencies[j],values[j].real(),values[j].imag(),condition,EOF);	 
 	  if ( (j>0) && (frequencies[j] == frequencies[j-1]) ) {
	    // remove duplicate frequencies
	    j--;
	    n--;
	  }
 	  j++;
     } while (condition != EOF);
     fclose(f);
 
  //for (j=0;j<=n-1;j++) printf("%ld %ld %13.8e %13.8e\n",n,j,frequencies[j],values[j].real());
 
  return;
  }

/**************************************************
*** get_value_from_interp: interpolate a frequency-
*** dependent quantity at a given frequency, and
*** convert it into a multiprecision complex
***************************************************/

  amp::campf<PRECISION> get_value_from_interp(double frequency, double *frequencies,
    std::complex<double> *values, long n){
    
    /* function interpolating at a given frequency the table defining the
    function f values=f(frequencies), and returns the corresponding
    multiprecision complex number. */
    
    double freq;
    std::complex<double> slopes[n],value;
    int condition;
    long j=0;
    amp::campf<PRECISION> valueMP;
    
    for (j=0;j<=n-1;j++) {slopes[j] = std::complex<double>(0.,0.);}
    value = interp(frequency, frequencies, values, slopes, n, 1); // linear interpolation
    valueMP.x = value.real();
    valueMP.y = value.imag();
    //printf("%13.8e %13.8e %13.8e\n",frequency,value.real(),value.imag());
 
    return valueMP;
  }

/**************************************************
*** get_eps_and_mu_from_layer_properties: compute
*** eps and mu for a given layer, from the layer
*** material properties
***************************************************/

  void get_eps_and_mu_from_layer_properties(double frequency,
    amp::ampf<PRECISION> rho, amp::ampf<PRECISION> tau, 
    amp::ampf<PRECISION> epsb, amp::ampf<PRECISION> tandelta, 
    amp::ampf<PRECISION> chi, amp::ampf<PRECISION> fmu,
    double *eps_frequencies, std::complex<double> *eps_values, long n_eps_frequencies, 
    double *mu_frequencies, std::complex<double> *mu_values, long n_mu_frequencies,
    int eps_vs_freq_kind, int mu_vs_freq_kind,
    amp::campf<PRECISION>& eps1, amp::campf<PRECISION>& mu1){
	
    /* function computing the relative complex permittivity and permeability of a layer,
       given its properties.
       
       eps_vs_freq_kind is:
       - 0 if no freq. dependent relative permittivity file is used,
       - 1 if a relative permittivity vs. frequency is given in a file,
       - 2 if a conductivity vs. frequency is given in a file.
       For 1 (resp. 2), frequencies are in eps_frequencies, and eps1 (resp. sigma)
       are in eps_values (complex numbers).
       The number of frequencies is n_eps_frequencies.
       
       mu_vs_freq_kind is:
       - 0 if no freq. dependent relative permeability file is used,
       - 1 if a relative permeability vs. frequency is given in a file
         (frequencies are in mu_frequencies, mu1 in mu_values, and the 
         number of frequencies in n_mu_frequencies). */
    
    amp::campf<PRECISION> jimagMP; // imaginary constant in multiprecision
    amp::ampf<PRECISION> omega; // angular frequency
    amp::campf<PRECISION> sigma;
    long n;
    
    omega=amp::twopi<PRECISION>()*frequency;
    
    jimagMP.x=0;jimagMP.y=1;
    
    // compute the relative complex permittivity
    if (eps_vs_freq_kind==1) {
      if ((rho.isFiniteNumber()) || (tau!=0) || (epsb!=1) || (tandelta!=0)) {
        printf("Cannot use rho, tau, epsb or tan(delta), when using a freq. dependent permittivity!\n");
        exit(2);
      }
      eps1=get_value_from_interp(frequency, eps_frequencies, eps_values, n_eps_frequencies);
    }
    else if (eps_vs_freq_kind==2) {
      if ((rho.isFiniteNumber()) || (tau!=0) || (tandelta!=0)) {
        printf("Cannot use rho, tau or tan(delta) when using a freq. dependent conductivity!\n");
        exit(2);
      }
      sigma=get_value_from_interp(frequency, eps_frequencies, eps_values, n_eps_frequencies);
      eps1=epsb+sigma/(jimagMP*eps0*omega);
    }
    else if (rho.isFiniteNumber()) {
      if (tandelta!=0) {
        printf("Cannot use both a resistivity and tan(delta)!\n");
        exit(2);
      }
      eps1=epsb+1/(jimagMP*eps0*rho*omega*(1+jimagMP*omega*tau));
    }
    else {
      if ((tandelta!=0) && (tau!=0)) {
        printf("Cannot use both a relaxation time and tan(delta)!\n");
        exit(2);
      }
      eps1=epsb*(1-jimagMP*tandelta);
    }

    // compute the relative complex permeability
    if (mu_vs_freq_kind==1) {
      if ((chi!=0) || (fmu.isFiniteNumber())) {
        printf("Cannot use chi or fmu, when using a freq. dependent permeability!\n");
        exit(2);
      }
      mu1=get_value_from_interp(frequency, mu_frequencies, mu_values, n_mu_frequencies);
    }
    else {
      mu1=1+chi/(1+jimagMP*omega/(amp::twopi<PRECISION>()*fmu));
    }
    
    return;
    
  }

/*************************************************
 *** read_general_input: read from std input the
 *** general input parameters - UNFINISHED & UNUSED
 *************************************************/

  void read_general_input(std::string *data, unsigned int n_input){
    
    /* Read from std input the general input parameters (machine, gamma,
    number of layers, frequency range, wake distance range, etc.).
    Chamber geometry (flat or round) is found from the lines:
     - "Number of upper/lower layers" (indicates flat)
     - "Top bottom symmetry (yes or no)" (indicates flat)
     - "Number of layers" (indicates round)
     - "Yokoya factors long, xdip, ydip, xquad, yquad" (indicates round)
     - "Geometry: flat"
     - "Geometry: round"
    */
     
     char *endline;
     std::string machine,topbot,commentoutput,dummy2;
     unsigned int N,M,dummy0; // number of upper and lower layers, then dummy parameter
     unsigned int nf_added,nz_added; /* number of individually added frequencies,
            number of individually added wake distances*/
     unsigned int flag_round,flag_topbotsym,ftypescan,nflog,nflin,ztypescan,nzlog,nzlin; /* flag for round geometry (1 if round),
            flag for top-bottom symmetry for flat geometries (1 if such a symmetry),
            type of frequency scan, number of freq. per decade, number of freq. in a lin. scan inside the log scan
            type of z scan, number of z per decade, number of z in a lin. scan inside the log scan*/
     amp::ampf<PRECISION> gamma,dummy3; // relativistic mass factor, dummy parameter
     double L,dummy1,yokoya[5]; // impedance length, dummy parameter, Yokoya factors
     double fminlog,fmaxlog,fminlin,fmaxlin,fsamplin,fadded[15];
     double zminlog,zmaxlog,zminlin,zmaxlin,dzlin,zadded[15];
     double tol; /* absolute error permitted on the integral of the difference (in norm)
            betweens impedances and their interpolation */
     double freqlin; // frequency limit above which we switch from a log mesh to a linear mesh
     size_t found;
     double factlong=100.; // multiplication factor for longitudinal impedance (precision should be better for long. wake) (default value)

     flag_round=1; // by default, assume round geometry
     flag_topbotsym=1; // by default, assume top/bottom symmetry for flat geometries

     // default values of the parameters - only for convergence and freq/z range
     // (material properties always have to be defined)
     ftypescan=0;
     fminlog=2;fmaxlog=13;nflog=10;
     fsamplin=8;fminlin=1;fmaxlin=2;nflin=100;
     ztypescan=0;
     zminlog=-2;zmaxlog=6;nzlog=10;nf_added=0;
     zminlin=1.e-2;zmaxlin=1;dzlin=1.e-2;nz_added=0;
     tol=1.e9;freqlin=1.e11;
     //N=2;M=2;L=1.;gamma="479.6";
     
     // identify each argument
     for (unsigned int i=1; i<=n_input; i++) {
       // general input parameters
       read_input(data[i],"Machine",dummy0,dummy1,machine,dummy3,2);
       read_input(data[i],"Relativistic Gamma",dummy0,dummy1,dummy2,gamma,3);
       read_input(data[i],"Impedance Length in m",dummy0,L,dummy2,dummy3,1);
       read_input(data[i],"Number of upper layers",N,dummy1,dummy2,dummy3,0);
       read_input(data[i],"Number of lower layers",M,dummy1,dummy2,dummy3,0);
       read_input(data[i],"Top bottom symmetry (yes or no)",dummy0,dummy1,topbot,dummy3,2);
       read_input(data[i],"Comments for the output files names",dummy0,dummy1,commentoutput,dummy3,2);
       if (data[i].find("Yokoya factors long, xdip, ydip, xquad, yquad") != std::string::npos){
         found=data[i].find("\t");
         if (found != std::string::npos){
           yokoya[0]=strtod(data[i].substr(found+1,std::string::npos).c_str(),&endline);
           for (int k=1; k<=4; k++) {
             yokoya[k]=strtod(endline,&endline);
           }
         }
       }
       // parameters related to frequency range explored
       read_input(data[i],"start frequency exponent (10^)",dummy0,fminlog,dummy2,dummy3,1);
       read_input(data[i],"stop frequency exponent (10^)",dummy0,fmaxlog,dummy2,dummy3,1);
       read_input(data[i],"linear (1) or logarithmic (0) or both (2) frequency scan",ftypescan,dummy1,dummy2,dummy3,0);
       read_input(data[i],"sampling frequency exponent (10^) in Hz (for linear)",dummy0,fsamplin,dummy2,dummy3,1);
       read_input(data[i],"Number of points per decade (for log)",nflog,dummy1,dummy2,dummy3,0);
       read_input(data[i],"when both, fmin of the refinement",dummy0,fminlin,dummy2,dummy3,1);
       read_input(data[i],"when both, fmax of the refinement",dummy0,fmaxlin,dummy2,dummy3,1);
       read_input(data[i],"when both, number of points in the refinement",nflin,dummy1,dummy2,dummy3,0);
       if (data[i].find("added frequencies") != std::string::npos){
         found=data[i].find("\t");
         if (found != std::string::npos){
           nf_added=1;fadded[nf_added]=strtod(data[i].substr(found+1,std::string::npos).c_str(),&endline);
           while (fadded[nf_added] != 0){
             nf_added++;fadded[nf_added]=strtod(endline,&endline);
           }
           nf_added--;
         }
       }
       // parameters related to wake distance range explored
       read_input(data[i],"exponent (10^) of zmin (in m) of the logarithmic sampling",dummy0,zminlog,dummy2,dummy3,1);
       read_input(data[i],"exponent (10^) of zmax (in m) of the logarithmic sampling",dummy0,zmaxlog,dummy2,dummy3,1);
       read_input(data[i],"linear (1) or logarithmic (0) or both (2) scan in z for the wake",ztypescan,dummy1,dummy2,dummy3,0);
       read_input(data[i],"sampling distance in m for the linear sampling",dummy0,dzlin,dummy2,dummy3,1);
       read_input(data[i],"Number of points per decade for the logarithmic sampling",nzlog,dummy1,dummy2,dummy3,0);
       read_input(data[i],"zmin in m of the linear sampling",dummy0,zminlin,dummy2,dummy3,1);
       read_input(data[i],"zmax in m of the linear sampling",dummy0,zmaxlin,dummy2,dummy3,1);
       read_input(data[i],"factor weighting the longitudinal impedance error",dummy0,factlong,dummy2,dummy3,1);
       read_input(data[i],"tolerance (in wake units) to achieve",dummy0,tol,dummy2,dummy3,1);
       read_input(data[i],"frequency above which the mesh bisecting is linear",dummy0,freqlin,dummy2,dummy3,1);
       if (data[i].find("added z") != std::string::npos){
         found=data[i].find("\t");
         if (found != std::string::npos){
           nz_added=1;zadded[nz_added]=strtod(data[i].substr(found+1,std::string::npos).c_str(),&endline);
           while (zadded[nz_added] != 0){
             nz_added++;zadded[nz_added]=strtod(endline,&endline);
           }
           nz_added--;
         }
       }
     }
     //printf("%s %s %d %d \n",machine.c_str(),topbot.c_str(),N,flag_topbotsym);
     //printf("%13.8e %13.8e %d\n",fadded[1],fadded[2],n_added);
     //printf("%13.8e %13.8e %d %d\n",fminlog,fmaxlog,nflog,ftypescan);
     //printf("%13.8e %13.8e %d\n",zadded[1],zadded[2],n_added);
     //printf("%13.8e %13.8e %d %ld\n",zminlog,zmaxlog,ztypescan,nzlog);
     //printf("%13.8e %13.8e %13.8e\n",freqlin,tol,factlong);

     // flag for top bottom symmetry (1 if there is such a symmetry)
     flag_topbotsym= ((strcmp(topbot.c_str(),"yes")==0 || strcmp(topbot.c_str(),"y")==0) || strcmp(topbot.c_str(),"1")==0);

  }

/**************************************************
 *** 			main program		            ***
 ***						                    ***
 **************************************************/
/*
main ()

{
 
 char *endline;
 char output[MAXCHARFILE],Zxdipoutput[MAXCHARFILE+10],Zydipoutput[MAXCHARFILE+10],Zlongoutput[MAXCHARFILE+10];
 char Zxquadoutput[MAXCHARFILE+10],Zyquadoutput[MAXCHARFILE+10],Zycstoutput[MAXCHARFILE+10],Input[MAXCHARFILE+10];
 std::string data[MAXLINES],machine,topbot,commentoutput,dummy2,epsfile,sigmafile,mufile;
 FILE *filZxdip, *filZydip, *filZxquad, *filZyquad, *filZycst, *filZlong, *filInput;
 unsigned int N,M,dummy0; // number of upper and lower layers, then dummy parameter
 unsigned int m,n,n_input,n_added,nf,nfdup; // indices of alphamn (azimuthal mode numbers), number of input lines,
 					//number of individually added frequencies, total number of
					//frequencies in the scan, number of duplicate frequencies
 unsigned int flag_topbotsym,typescan,nflog,nflin; // flag for top-bottom symmetry (1 if such a symmetry), type of frequency scan
 		//number of freq. per decade, number of freq. in a lin. scan inside the log scan
 ap::template_1d_array< amp::campf<Precision> > eps1,eps1m,mu1,mu1m; // eps1 and mu1 for upper (without m at
 					//the end) and lower layers (with m at the end) 
 ap::template_1d_array< amp::ampf<Precision> > b,bm,thick,thickm; // position of upper and lower boundaries; thickness of the layers 
 ap::template_1d_array< amp::ampf<Precision> > rho,tau,epsb,tandelta,chi,fmu,rhom,taum,epsbm,tandeltam,chim,fmum; // layers
 						//properties (DC resistivity, resistivity relaxation time,
						//dielectric constant, tan(delta), magnetic susceptibility=mu_r-1,
						//relaxation frequency of permeability)
 amp::ampf<Precision> omega,beta,k,gamma,kovergamma,u,dummy3; // parameters
 amp::campf<Precision> jimagMP; // imaginary constant in multiprecision
 gsl_integration_workspace *w;
 size_t limit=1000; // limit (number of intervals) for gsl integration algorithm
 double tolintrel=1.e-15; // relative error permitted (w.r.t. to value at the previous frequency) for integration
 double x,y,L,fminlog,fmaxlog,fminlin,fmaxlin,fsamplin,fadded[15],*freq,dif,dummy1;
 std::complex<double> Zxdip,Zydip,Zxquad,Zyquad,Zlong,Zycst,cst;
 std::complex<double> alpha00,alpha01,alpha02,alpha11; // alphamn constants and imaginary constant
 std::complex<double> tolintabs00,tolintabs01,tolintabs02,tolintabs11;
 double **eps_frequencies, **mu_frequencies, **epsm_frequencies, **mum_frequencies;
 std::complex<double> **eps_values, **mu_values, **epsm_values, **mum_values;
 int *eps_vs_freq_kind, *mu_vs_freq_kind, *epsm_vs_freq_kind, *mum_vs_freq_kind;
 long *n_eps_frequencies, *n_mu_frequencies, *n_epsm_frequencies, *n_mum_frequencies;
 size_t found;
 time_t start,end; // times
					
 // start time
 time(&start);
 
 // memory of etas and chis allocation
 kxmem.setbounds(0,MAXMEM);
 eta1mem.setbounds(0,MAXMEM);eta2mem.setbounds(0,MAXMEM);
 chi1mem.setbounds(0,MAXMEM);chi2mem.setbounds(0,MAXMEM);
 mem=0;

 // imaginary constant
 jimagMP.x=0;jimagMP.y=1;
 
 // default values of the parameters (in case)
 flag_topbotsym=1;
 typescan=0;
 fminlog=2;fmaxlog=13;nflog=10;n_added=0;
 fsamplin=8;fminlin=1;fmaxlin=2;nflin=100;nf=0;
 N=2;M=2;L=1.;gamma="479.6";
 
 // read input file
 // first read everything to identify the strings in front of each parameters
 n_input=0;
 while (std::cin.eof()==0) {
   std::getline (std::cin,data[n_input+1]);
   // next line is when the input file comes from windows or else and has some ^M characters in the
   //end of each line
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
   read_input(data[i],"start frequency exponent (10^)",dummy0,fminlog,dummy2,dummy3,1);
   read_input(data[i],"stop frequency exponent (10^)",dummy0,fmaxlog,dummy2,dummy3,1);
   read_input(data[i],"linear (1) or logarithmic (0) or both (2) frequency scan",typescan,dummy1,dummy2,dummy3,0);
   read_input(data[i],"sampling frequency exponent (10^) in Hz (for linear)",dummy0,fsamplin,dummy2,dummy3,1);
   read_input(data[i],"Number of points per decade (for log)",nflog,dummy1,dummy2,dummy3,0);
   read_input(data[i],"when both, fmin of the refinement",dummy0,fminlin,dummy2,dummy3,1);
   read_input(data[i],"when both, fmax of the refinement",dummy0,fmaxlin,dummy2,dummy3,1);
   read_input(data[i],"when both, number of points in the refinement",nflin,dummy1,dummy2,dummy3,0);
   read_input(data[i],"Comments for the output files names",dummy0,dummy1,commentoutput,dummy3,2);
   if (data[i].find("added frequencies") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos){
       n_added=1;fadded[n_added]=strtod(data[i].substr(found+1,std::string::npos).c_str(),&endline);
       while (fadded[n_added] != 0){
         n_added++;fadded[n_added]=strtod(endline,&endline);
       }
       n_added--;
     }
   }
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
       eps_frequencies[p+1] = new double[n_eps_frequencies[p+1]+1];
       eps_values[p+1] = new std::complex<double>[n_eps_frequencies[p+1]+1];
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
         epsm_values[p+1] = new complex<double>[n_epsm_frequencies[p+1]+1];
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
         epsm_values[p+1] = new complex<double>[n_epsm_frequencies[p+1]+1];
         read_complex_file(sigmafile,epsm_frequencies[p+1],epsm_values[p+1],n_epsm_frequencies[p+1]);
       }
       if (!mufile.empty()) {
         mum_vs_freq_kind[p+1]=1;
         read_complex_file_length(mufile,n_mum_frequencies[p+1]);
         mum_frequencies[p+1] = new double[n_mum_frequencies[p+1]+1];
         mum_values[p+1] = new complex<double>[n_mum_frequencies[p+1]+1];
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
   printf("%lg %lg %lg %lg %lg %lg \n",double(rho(p).toDouble()),double(tau(p).toDouble()),double(epsb(p).toDouble()),
   	double(chi(p).toDouble()),double(fmu(p).toDouble()),double(b(p).toDouble()));
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
 for (unsigned int i=1;i<=n_added; i++) {
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
	1e3*double(amp::ampf<Precision>(b(1)).toDouble()),commentoutput.c_str());
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
 
 alpha00=std::complex<double>(0.,0.);
 alpha01=std::complex<double>(0.,0.);
 alpha02=std::complex<double>(0.,0.);
 alpha11=std::complex<double>(0.,0.);
 tolintabs00=std::complex<double>(0.,0.);
 tolintabs01=std::complex<double>(0.,0.);
 tolintabs02=std::complex<double>(0.,0.);
 tolintabs11=std::complex<double>(0.,0.);
 
 // impedance computation at each frequency: beginning of the loop
 for (unsigned int i=0; i<nf; i++) {
    
   mem=0;  // initialize memory at eahc frequency
   omega=amp::twopi<Precision>()*freq[i];
   k=omega/(beta*C);
   kovergamma=k/gamma;
   
   // computes the layer properties for the angular freq. omega
   for (unsigned int p=2;p<=N+1; p++) {
     get_eps_and_mu_from_layer_properties(freq[i], rho(p), tau(p), epsb(p),
        tandelta(p), chi(p), fmu(p), eps_frequencies[p], eps_values[p],
        n_eps_frequencies[p], mu_frequencies[p], mu_values[p], n_mu_frequencies[p],
        eps_vs_freq_kind[p], mu_vs_freq_kind[p], eps1(p), mu1(p));
     
     printf("%lg Hz: layer %d -> %lg m, eps1= %lg + %lg j, mu1= %lg + %lgj\n",freq[i],p,double(b(p).toDouble()),
        double(eps1(p).x.toDouble()),double(eps1(p).y.toDouble()),
        double(mu1(p).x.toDouble()),double(mu1(p).y.toDouble()));
     }
   for (unsigned int p=2;p<=M+1; p++) {
     get_eps_and_mu_from_layer_properties(freq[i], rhom(p), taum(p), epsbm(p),
        tandeltam(p), chim(p), fmum(p), epsm_frequencies[p], epsm_values[p],
        n_epsm_frequencies[p], mum_frequencies[p], mum_values[p], n_mum_frequencies[p],
        epsm_vs_freq_kind[p], mum_vs_freq_kind[p], eps1m(p), mu1m(p));

     printf("%lg Hz: layer %d -> %lg m, eps1= %lg + %lg j, mu1= %lg + %lgj\n",freq[i],-p,double(bm(p).toDouble()),
        double(eps1m(p).x.toDouble()),double(eps1m(p).y.toDouble()),
     	double(mu1m(p).x.toDouble()),double(mu1m(p).y.toDouble()));
     }

   //printf("%s %s\n",b(1).toDec().c_str(),bm(1).toDec().c_str());
   //printf("%s %s %s\n",omega.toDec().c_str(),k.toDec().c_str(),kovergamma.toDec().c_str());
 
   // computes alpha00
   x=integrate(1,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,0,0,tolintabs00.real(),limit,w);
   y=integrate(0,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,0,0,tolintabs00.imag(),limit,w);
   alpha00=std::complex<double>(x,y);
   //printf("alpha00: %13.8e %13.8e\n",alpha00.real(),alpha00.imag());
 
   // computes alpha01
   if (flag_topbotsym==0) {
     // computes alpha01
     x=integrate(1,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,0,1,tolintabs01.real(),limit,w);
     y=integrate(0,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,0,1,tolintabs01.imag(),limit,w);
     alpha01=std::complex<double>(x,y);
     }
   else alpha01=std::complex<double>(0.,0.);
   //printf("alpha01: %13.8e %13.8e\n",alpha01.real(),alpha01.imag());
 
   // computes alpha02
   x=integrate(1,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega, k,kovergamma,0,2,tolintabs02.real(),limit,w);
   y=integrate(0,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,0,2,tolintabs02.imag(),limit,w);
   alpha02=std::complex<double>(x,y);
   //printf("alpha02: %13.8e %13.8e\n",alpha02.real(),alpha02.imag());

   // computes alpha11
   x=integrate(1,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,1,1,tolintabs11.real(),limit,w);
   y=integrate(0,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,1,1,tolintabs11.imag(),limit,w);
   alpha11=std::complex<double>(x,y);
   //printf("alpha11: %13.8e %13.8e\n",alpha11.real(),alpha11.imag());
   
   // computes and writes the impedances
   cst=jimag*L*double(amp::ampf<Precision>(k*Z0/(beta*amp::sqr(gamma)*amp::twopi<Precision>())).toDouble());
   Zlong=cst*alpha00;
   Zycst=cst*alpha01/double(amp::ampf<Precision>(gamma).toDouble());
   cst=cst*double(amp::ampf<Precision>(k/(2*amp::sqr(gamma))).toDouble());
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
 
  
   tolintabs00=tolintrel*std::complex<double>(std::abs(alpha00.real()),std::abs(alpha00.imag()));  
   tolintabs01=tolintrel*std::complex<double>(std::abs(alpha01.real()),std::abs(alpha01.imag()));  
   tolintabs02=tolintrel*std::complex<double>(std::abs(alpha02.real()),std::abs(alpha02.imag()));  
   tolintabs11=tolintrel*std::complex<double>(std::abs(alpha11.real()),std::abs(alpha11.imag()));
   //std::cout << "tolintabs02: " << tolintabs02.real() << " " << tolintabs02.imag() << "\n";

   }
   // end of loop on frequencies
 
 
 gsl_integration_workspace_free(w);
  
 fclose(filZxdip);
 fclose(filZydip);
 fclose(filZxquad);
 fclose(filZyquad);
 fclose(filZlong);
 fclose(filZycst);
 
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
*/
