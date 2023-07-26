#!/usr/bin/python

# library to compute Fourier integrals of smooth functions
# with a non-equidistant mesh

import sys
import subprocess
pymod=subprocess.check_output("echo $PYMOD",shell=True).strip().decode()
if pymod.startswith('local'):
    py_numpy=subprocess.check_output("echo $PY_NUMPY",shell=True).strip().decode()
    sys.path.insert(1,py_numpy)
    py_scipy=subprocess.check_output("echo $PY_SCIPY",shell=True).strip().decode()
    sys.path.insert(1,py_scipy)
from string import *
import numpy as np
from numpy import fft
import pylab,re,os
from tables_lib import *
from C_complex import *
from ctypes import *

def fourier_integral_asympt_correction(flaginf,tscan,omegai_end,fi_end,df_end=0.):

    ''' correcting term for the semi-infinite Fourier integral at tscan of the
    function f, above (if flaginf>0) or below (if flaginf<0) omegai_end.
    fi_end is the value of f at omegai_end, df_end it's derivative (can be zero).
    tscan can be a scalar or a list (or array). Result is always an array (with a single
    element if tscan is a scalar).
    '''

    tscan=create_list(tscan,n=1);
    res=np.zeros(len(tscan),dtype=complex);

    for it,t in enumerate(tscan):

        res[it] = -np.sign(flaginf)*np.exp(1j*t*omegai_end)*(fi_end/(1j*t) + df_end/(t**2));

    return res;



def fourier_integral_any_sampling(tscan,omegai,fi,df=0.,interp=0,flaginf=1,eps=1e-20,di=None):

    '''compute semi-infinite Fourier integral from omegai[0] to omegai[-1]
    (or +/- infinity if flaginf is True) of a function f, i.e. the integral
    of f(omega)*exp(j*omega*t) for any t in tscan (table) and for f given by
    its interpolation table fi of type 'interp' at the points omegai
    (0 for cubic 'pchip' interpolation, 1 for linear interpolation). In general
    f is complex.

    if flaginf is non-zero, a correcting term is added (at each t) to take into
    account the semi-infinite integral above (if flaginf>0) or below (if flaginf<0)
    the angular frequency range omegai. Then df gives the derivative of the function
    at the point of the range closest to this semi-infinite integral. It can be zero.

    interp can be a single integer or a tables of them, of the length of omegai -1,
    to specify a different kind of interpolation on each interval, if needed.
    eps is the relative precision on the computation of the functions used in the algorithm.

    If di is None, pchip derivatives are computed here, otherwise it should
    provide the pchip derivatives of the interpolation.

    We use Filon's type method on each interval, plus an asymptotic method for the correcting term toward infinity if flaginf=1.
    '''

    libIW2D=CDLL("libIW2D.so");

    n=len(omegai);
    # some test on arrays lengths
    interp=create_list(interp,n=n-1);
    if (len(interp)!=(n-1)): print("  Pb in fourier_integral_any_sampling: incorrect length of interp !");sys.exit();
    if (len(fi)!=n): print("  Pb in fourier_integral_any_sampling: fi and omegai not of same length !");sys.exit();

    # build C-compatible array types and functions
    DArray=c_double * n;
    LDArray=c_longdouble * n;
    LDArraym1=c_longdouble * (n-1);
    IArraym1=c_uint * (n-1);
    CArray=Complex * n;

    prototypepchip = CFUNCTYPE(c_void_p, CArray, DArray, CArray, c_ulong);
    pchip_derivatives = prototypepchip(('pchip_derivatives', libIW2D));
    prototypefourier = CFUNCTYPE(Complex, CArray, CArray, Complex, c_longdouble, LDArray,
        LDArraym1, c_ulong, c_double, IArraym1, c_int); # 1st one is the result
    fourier = prototypefourier(('fourier_integral_inf', libIW2D));

    # build C-compatible arrays
    # angular frequencies (double and long double versions)
    omegai_list=omegai.tolist();omegai_arr=DArray(*omegai_list);
    omegai_ld=np.longdouble(omegai);omegai_ldlist=omegai_ld.tolist();
    omegai_ldarr=LDArray(*omegai_ldlist);
    # differences between successive omegai
    delta=np.diff(omegai_ld);delta_list=delta.tolist();delta_arr=LDArraym1(*delta_list);
    # interpolating values of the function
    fi_list=[complex_to_Complex(elem) for elem in fi.tolist()];fi_arr=CArray(*fi_list);
    # interpolation types
    interp_arr=IArraym1(*interp);
    # array for pchip derivatives
    flagdi=False; # case when pchip derivatives are already computed
    if (di is None):
        di=np.zeros(n,dtype=complex);flagdi=True;
    di_list=[complex_to_Complex(elem) for elem in di.tolist()];
    di_arr=CArray(*di_list);

    if (flagdi)and(any(np.array(interp)==0)):
        # compute pchip derivatives
        pchip_derivatives(di_arr,omegai_arr,fi_arr,c_ulong(n));

    # table with resulting Fourier integrals for each t in tscan
    fint=np.zeros(len(tscan),dtype=complex);
    flaginf_int=0; # deactivate semi-infinite asymptotic correction in the C function
    # (add this correction here when flaginf is non-zero)

    for it,t in enumerate(tscan):

        res=fourier(fi_arr,di_arr,complex_to_Complex(df),c_longdouble(t),omegai_ldarr,
                delta_arr,c_ulong(n),c_double(eps),interp_arr, c_int(flaginf_int));

        fint[it]=res.real+1j*res.imag;

    if (flaginf!=0):
        # add correcting term for semi-infinite integral above or below the highest
        # omegai in absolute value
        if (flaginf>0): omegai_end=omegai[-1];fi_end=fi[-1];
        else: omegai_end=omegai[0];fi_end=fi[0];

        fint += fourier_integral_asympt_correction(flaginf,tscan,omegai_end,fi_end,df_end=df);


    return fint;


def reinterpolate_linear_sampling(tscan,omegai,fi,interp_type='linear'):

    ''' Make tscan and omegai be equidistant meshes with some correspondance
    (for DFT purposes).
    deltaomega is based on the difference between the 2 first points of the
    initial mesh (then keep the same beginning and end of the initial omegai range),
    while deltat is fixed to allow DFT computation:
    deltat= 2*pi/((n-1)*deltaomega) )
    where n is defined as 4*len(omegainew) (~ arbitrary, from "Numerical Recipes in C",
    2nd ed, chap. 13, p 586).
    tscan is then changed according to this deltat and such that the last
    point of the scan must be tscan[0]+2.*np.pi/deltaomega.

    fi is re-interpolated on the new mesh, based on interpolation scheme in
    interp_type ('linear', 'cubic', etc.). NOTE: 'cubic' is very slow.

    Then compute the interpolation table of the function g defined by:
    g(i)=f(omegai)*exp(1j*deltaomega*i*tscan[0])
    '''

    from scipy import interpolate as itp
    from copy import deepcopy

    # replace t and omega tables by equidistant tables (if they already were
    # and if deltat and deltaomega fit with the DFT, it won't change anything)
    deltaomega=omegai[1]-omegai[0];
    omegainew=np.arange(omegai[0],omegai[-1]+deltaomega/100.,deltaomega);
    m=len(omegainew);
    if (m!=len(omegai))or( np.abs(omegainew[-1]-omegai[-1]) > 1e-10*np.abs(omegai[-1]) ):
        print(" Warning in reinterpolate_linear_sampling: omegai has been changed for DFT, old length=",len(omegai),", new length=",n,", old last elem.=",omegai[-1],", new last elem.=",omegainew[-1]);
        flag_reinterp=True;
    else:
        flag_reinterp=(np.max(np.abs(omegai-omegainew))>0.);#print flag
        print(" reinterpolate_linear_sampling: max difference between old and new omegai sampling:",np.max(np.abs(omegai-omegainew)))

    # choose new tscan according to deltaomega
    tscanend=tscan[0]+2.*np.pi/deltaomega;
    # number of t (n) (see "Numerical Recipes in C", 2nd ed, chap. 13, p 586)
    n=4*m;
    deltat=2.*np.pi/(deltaomega*(n-1.));
    if ( np.abs(deltat-(tscan[1]-tscan[0])) > 1e-10*np.abs(deltat) ):
        print(" Warning in reinterpolate_linear_sampling: deltat has been changed for DFT, old=",tscan[1]-tscan[0],", new=",deltat);
    tscannew=np.arange(tscan[0],tscanend+deltat/100.,deltat)
    print(" tscanold[0]=",tscan[0],", tscanold[-1]=",tscan[-1]);
    print(" tscannew[0]=",tscannew[0],", tscannew[-1]=",tscannew[-1]);
    print(" len(tscannew)=",len(tscannew),", len(omegainew)=",m);

    if flag_reinterp:
        # re-interpolate f on the new equidistant angular frequencies
        f=itp.interp1d(omegai,fi,kind=interp_type,bounds_error=False,fill_value=fi[-1]); # interpolating function
        finew=f(omegainew); # interpolated values
    else:
        finew=deepcopy(fi);omegainew=deepcopy(omegai);

    # function on which to do the FFT (different from f to take into account
    # a possible non-zero tscan[0])
    gi=finew*np.exp(1j*deltaomega*np.arange(m)*tscannew[0]);

    if (len(tscannew)>m):
        # zero-pad gi (it's only for the FFT)
        ginew=np.hstack((gi,np.zeros(len(tscannew)-m,dtype=complex)));


    return tscannew,omegainew,finew,ginew,deltaomega;


def fourier_integral_linear_sampling(tscan,omegai,fi,df=0.,interp=0,flaginf=1):

    '''compute semi-infinite Fourier integral from omegai[0] to omegai[-1]
    (or +/- infinity if flaginf is True) of a function f, i.e. the integral
    of f(omega)*exp(j*omega*t) for any t in tscan (table) and for f given by
    its interpolation table fi of type 'interp' at the points omegai
    (0 for cubic spline interpolation, 1 for linear interpolation). In general
    f is complex.

    if flaginf is non-zero, a correcting term is added (at each t) to take into
    account the semi-infinite integral above (if flaginf>0) or below (if flaginf<0)
    the angular frequency range omegai. Then df gives the derivative of the function
    at the point of the range closest to this semi-infinite integral. It can be zero.

    For this function interp cannot be a table (only one kind of interpolation
    over the full range).
    In this version tscan and omegai should be made of equidistant points and
    with a certain relation between them to allow the computation by a DFT.
    If they are not, such meshes are anyway constructed by the function
    'reinterpolate_linear_sampling' just above.

    We use Filon's type method on each interval; computation based on the algorithm
    found in Numerical Recipes in C (2nd ed, chap. 13, p. 586) plus an
    asymptotic method for the correcting term toward infinity if flaginf=1.
    '''

    if interp==0: interp_type='cubic';
    else: interp_type='linear';

    # reinterpolate to get equidistant samplings, DFT compatible (if needed)
    tscannew,omegainew,finew,gi,deltaomega=reinterpolate_linear_sampling(tscan,omegai,fi,interp_type='linear');

    # FFT of g
    gfft=fft.fft(gi);

    fint=np.zeros(len(tscannew),dtype=complex);

    for jt,t in enumerate(tscannew):

        w,alpha0,alpha1,alpha2,alpha3=correction_functions(t*deltaomega,interp=interp_type);

        fint[jt]=deltaomega*np.exp(1j*omegainew[0]*t)*(w*gfft[jt] +
                alpha0*finew[0]+alpha1*finew[1]+alpha2*finew[2]+alpha3*finew[3] +
                np.conj(alpha0)*finew[-1]+np.conj(alpha1)*finew[-2]+np.conj(alpha2)*finew[-3]+np.conj(alpha3)*finew[-4]);

    if (flaginf!=0):
        # add correcting term for semi-infinite integral above or below the highest
        # omegai in absolute value
        if (flaginf>0): omegai_end=omegainew[-1];fi_end=finew[-1];
        else: omegai_end=omegainew[0];fi_end=finew[0];

        fint += fourier_integral_asympt_correction(flaginf,tscannew,omegai_end,fi_end,df_end=df);


    return tscannew,fint;


def fourier_integral_linear_sampling_naive(tscan,omegai,fi,df=0.,flaginf=1):

    '''compute semi-infinite Fourier integral from omegai[0] to omegai[-1]
    (or +/- infinity if flaginf is True) of a function f, i.e. the integral
    of f(omega)*exp(j*omega*t) for any t in tscan (table) and for f given by
    its interpolation table fi at the points omegai. In general f is complex.

    if flaginf is non-zero, a correcting term is added (at each t) to take into
    account the semi-infinite integral above (if flaginf>0) or below (if flaginf<0)
    the angular frequency range omegai. Then df gives the derivative of the function
    at the point of the range closest to this semi-infinite integral. It can be zero.

    In this version tscan and omegai should be made of equidistant points and
    with a certain relation between them to allow the computation by a DFT.
    If they are not, such meshes are anyway constructed by the function
    'reinterpolate_linear_sampling' just above.

    Here we use a "naive" FFT algorithm, that is generally BAD - see
    Numerical Recipes in C (2nd ed, chap. 13, p. 585). It is described in this
    reference and also in N. Mounet PhD thesis (EPFL 5305), p. 46.
    '''

    # reinterpolate to get equidistant samplings, DFT compatible (if needed)
    tscannew,omegainew,finew,gi,deltaomega=reinterpolate_linear_sampling(tscan,omegai,fi,interp_type='linear');

    # FFT of g
    gfft=fft.fft(gi);

    fint=np.zeros(len(tscannew),dtype=complex);

    for jt,t in enumerate(tscannew):

        fint[jt]=deltaomega*np.exp(1j*omegainew[0]*t)*gfft[jt];

    if (flaginf!=0):
        # add correcting term for semi-infinite integral above or below the highest
        # omegai in absolute value
        if (flaginf>0): omegai_end=omegainew[-1];fi_end=finew[-1];
        else: omegai_end=omegainew[0];fi_end=finew[0];

        fint += fourier_integral_asympt_correction(flaginf,tscannew,omegai_end,fi_end,df_end=df);


    return tscannew,fint;


def correction_functions(theta,interp='cubic',limit=1e-3):

    '''' computes correction functions W and alpha_j for Fourier integrals, as
    defined in Numerical Recipes in C (2nd ed, chap. 13, p. 586-587)
    - theta=t*deltaomega
    - interp is the interpolation type: 'linear' or 'cubic'
    - limit: limit in theta below which we use the truncated Taylor series
    instead of the exact formula
    '''

    if interp.startswith('linear'):

        # linear interpolation
        alpha1=0.;alpha2=0.;alpha3=0.;
        theta2=theta**2;

        if (theta<limit):
            # small parameter approximations
            theta4=theta2**2;theta6=theta2**3;
            w=1-theta2/12.+theta4/360.-theta6/20160.;
            alpha0=-w/2. + 1j*theta*(1./6.-theta2/120.+theta4/5040.-theta6/362880.);

        else:
            # general case
            w=2*(1.-np.cos(theta))/theta2;

            alpha0=-w/2. + 1j*(theta-np.sin(theta))/theta2;

    elif interp.startswith('cubic'):

        # cubic interpolation
        theta2=theta**2;theta4=theta2**2;

        if (theta<limit):
            # small parameter approximations
            theta6=theta2**3;

            w=1-11.*theta4/720.+23.*theta6/15120.;

            alpha0=-2./3.+theta2/45.+103.*theta4/15120.-169.*theta6/226800 + 1j*theta*(2./45.+2.*theta2/105.-8.*theta4/2835.+86.*theta6/467775.);

            alpha1=7./24.-7.*theta2/180.+5.*theta4/3456.-7.*theta6/259200. + 1j*theta*(7./72.-theta2/168.+11.*theta4/72576.-13.*theta6/5987520.);

            alpha2=-1./6.+theta2/45.-5.*theta4/6048.+theta6/64800. + 1j*theta*(-7./90.+theta2/210.-11.*theta4/90720.+13.*theta6/7484400.);

            alpha3=1./24.-theta2/180.+5.*theta4/24192.-theta6/259200. + 1j*theta*(7./360.-theta2/840.+11.*theta4/362880.-13.*theta6/29937600.);

        else:
            # general case
            costheta=np.cos(theta);
            sintheta=np.sin(theta);
            cos2theta=np.cos(2.*theta);
            sin2theta=np.sin(2.*theta);

            w=(6.+theta2)*(3.-4.*costheta+cos2theta)/(3.*theta4);

            alpha0=((-42.+5*theta2)+(6.+theta2)*(8.*costheta-cos2theta))/(6.*theta4) + 1j*(theta*(-12.+6.*theta2)+(6+theta2)*sin2theta)/(6.*theta4);

            alpha1=(14.*(3.-theta2)-7.*(6.+theta2)*costheta)/(6.*theta4) + 1j*(30.*theta-5*(6.+theta2)*sintheta)/(6.*theta4);

            alpha2=(-4.*(3.-theta2)+2.*(6.+theta2)*costheta)/(3.*theta4) + 1j*(-12.*theta+2.*(6+theta2)*sintheta)/(3.*theta4);

            alpha3=(2.*(3.-theta2)-(6.+theta2)*costheta)/(6.*theta4) + 1j*(6.*theta-(6.+theta2)*sintheta)/(6.*theta4);

    else:
        print(" Pb in correction_functions: interp not recognized !");sys.exit();

    return w,alpha0,alpha1,alpha2,alpha3;


def int_Simpson(a,b,fi):

    '''integrate a function given by its values in 3 points (in table fi),
    between a and b, with a simple Simpson's rule (source: Wikipedia).
    The 3 points where f is given must be a, (a+b)/2 and b.
    '''

    res = (b-a)*(fi[0]+4.*fi[1]+fi[2])/6.;

    return res;


def compute_intercalate_pchip_derivatives(i,omegai,fi,di):

    ''' compute pchip_derivatives from omegai and fi, and insert them at the places
    i, i+1 and i+2 in table di (shifting di to the right above i+1).
    '''

    libIW2D=CDLL("libIW2D.so");

    # build C-compatible array types and functions
    n=min(5,len(omegai));#print "n=",n;
    DArray=c_double * n;
    CArray=Complex * n;
    prototypepchip = CFUNCTYPE(c_void_p, CArray, DArray, CArray, c_ulong);
    pchip_derivatives = prototypepchip(('pchip_derivatives', libIW2D));

    # build C-compatible arrays
    # angular frequencies
    omegai_list=omegai[max(0,i-1):min(i+4,len(omegai))].tolist();
    omegai_arr=DArray(*omegai_list);
    # interpolating values of the function
    fi_list=[complex_to_Complex(elem) for elem in fi[max(0,i-1):min(i+4,len(omegai))].tolist()];
    fi_arr=CArray(*fi_list);
    # array for pchip derivatives
    di_list=[complex_to_Complex(elem) for elem in np.zeros(n,dtype=complex)];
    di_arr=CArray(*di_list);

    # compute pchip derivatives
    pchip_derivatives(di_arr,omegai_arr,fi_arr,c_ulong(n));
    # replace the values at the boundaries of the interval [omegai[i],omegai[i+1]],
    # and intercalate the derivative at the middle point.
    offset=min(0,i-1); # offset of -1 in case i=0 (because then we cannot
    # fetch an additional point on the left, to compute the derivatives)
    if (i>=0): di[i]=Complex_to_complex(di_arr[1+offset]);
    if (i<len(di)-1): di[i+1]=Complex_to_complex(di_arr[3+offset]);
    dinew=add_element_in_array(di,i,Complex_to_complex(di_arr[2+offset]));

    return dinew;


def bisect_and_calculate(f,i,omegai,fi,di,bisect_type='log',pchip_flag=True):

    ''' bisect the interval at the right of omegai[i],
    and recalculate and redefine accordingly omegai, fi (values of
    the function f at omegai) and the pchip derivatives di if
    pchip_flag is True.
    bisect_type can be 'log' (middle point in log scale) or 'lin'
    (middle point in linear scale)
    It also computes the interpolated value (before bisection) at the middle point.
    '''

    libIW2D=CDLL("libIW2D.so");

    # to compute the interpolated value
    DArray2=c_double * 2;
    CArray2=Complex * 2;
    prototypeinterp = CFUNCTYPE(Complex, c_double, DArray2, CArray2, CArray2, c_ulong, c_uint);
    interp = prototypeinterp(('interp', libIW2D));
    if pchip_flag: interp_int=0;
    else: interp_int=1;

    if bisect_type.startswith('lin'):
        omegai_added=(omegai[i]+omegai[i+1])/2.;
    else:
        omegai_added=np.exp((np.log(omegai[i])+np.log(omegai[i+1]))/2.);

    fi_added=f(omegai_added);
    # interpolated value at the same point (before bisection)
    omegai2_list=[omegai[i],omegai[i+1]];omegai2_arr=DArray2(*omegai2_list);
    fi2_list=[complex_to_Complex(fi[i]),complex_to_Complex(fi[i+1])];fi2_arr=CArray2(*fi2_list);
    try:
        di2_list=[complex_to_Complex(di[i]),complex_to_Complex(di[i+1])];
    except TypeError or IndexError:
        di2_list=[complex_to_Complex(0.+1j*0.),complex_to_Complex(0.+1j*0.)];

    di2_arr=CArray2(*di2_list);

    interp_value=Complex_to_complex(interp(c_double(omegai_added),omegai2_arr,fi2_arr,di2_arr,c_ulong(2),c_uint(interp_int)));


    # construct new tables
    omegainew=add_element_in_array(omegai,i,omegai_added);
    finew=add_element_in_array(fi,i,fi_added);

    if pchip_flag:

        dinew=compute_intercalate_pchip_derivatives(i,omegainew,finew,di)

    else:

        dinew=add_element_in_array(di,i,0.);

    return omegainew,finew,dinew,interp_value;


def fourier_int_convergence_boundary(f,df,omegabound,tscan,omegai,fi,di,omegainew,finew,dinew,min_or_max,
        pchip_flag=True,tolrel=0.001,criterion='abs'):

    ''' convergence loop for Fourier integral, vs minimum (if min_or_max=0) or maximum
    (if min_or_max !=0) angular frequency (given in omegabound) of the scan omegai.
    functions values at the omegai's are in fi, pchip derivatives (used when pchip_flag is True) in di.
    f and df are the function and its derivative.
    Criterion can be:
        - abs: difference in absolute value matters,
        - real: difference of real part matters,
        - imag: difference of imaginary part matters,
        - realimag: difference of both real & imaginary part matters.
    '''

    from copy import deepcopy

    # indices where to insert new point (at the beginning or the end of the freq. scan)
    if (min_or_max!=0): i=len(omegai)-1;inew=len(omegainew)-1;omegab=omegai[-1];ind=np.arange(-2,0);
    else: i=-1;inew=-1;omegab=omegai[0];ind=np.arange(2);

    if pchip_flag: interp_int=0;
    else: interp_int=1;

    if ~np.isinf(omegabound): flaginf=0;
    else: flaginf=np.sign(omegabound);

    # first trial
    fint=fourier_integral_any_sampling(tscan,omegai,fi,df=df(omegab),interp=interp_int,flaginf=flaginf,di=di);

    flag_continue=(omegab!=omegabound);


    while flag_continue:

        # here, convergence is on the Fourier integral itself
        # add one point
        if ~np.isinf(omegabound): omegai_added=(omegabound+omegab)/2.;
        else: omegai_added=omegab*2.; # note: we assume omegabound and omegab of same sign...


        fi_added=f(omegai_added);df_added=df(omegai_added);
        omegai=add_element_in_array(omegai,i,omegai_added);
        fi=add_element_in_array(fi,i,fi_added);
        if pchip_flag: di=compute_intercalate_pchip_derivatives(i,omegai,fi,di);
        else: di=add_element_in_array(di,i,0.);

        # also add this point to the refined sampling
        omegainew=add_element_in_array(omegainew,inew,omegai_added);
        finew=add_element_in_array(finew,inew,fi_added);
        if pchip_flag: dinew=compute_intercalate_pchip_derivatives(inew,omegainew,finew,dinew)
        else: dinew=add_element_in_array(dinew,inew,0.);

        # subtract old asymptotic part if needed
        if (flaginf==0):fintnew=deepcopy(fint);
        else: fintnew=fint - fourier_integral_asympt_correction(np.sign(omegab),tscan,omegab,f(omegab),df_end=df(omegab));
        # add contribution of added interval (plus asymptotic part when flaginf==1)
        fintnew += fourier_integral_any_sampling(tscan,omegai[ind],fi[ind],df=df_added,interp=interp_int,flaginf=flaginf,di=di[ind]);

        # max relative error
        maxerr,iterr=max_diff_rel_complex_array(fintnew,fint,criterion=criterion);

        flag_continue=(maxerr>tolrel);
        if (min_or_max!=0): i+=1;inew+=1;

        fint=deepcopy(fintnew);omegab=omegai_added;

    return omegai,fi,di,omegainew,finew,dinew;


def fourier_int_converged(f,df,tscan,omegamin,omegamax,tolabs=1e10,tolrel=1e-3,
        omegamin_ini=None,omegamax_ini=None,criterion='absabs',bisect='log',interp='pchip'):

    '''compute semi-infinite Fourier integral from omegamin to omegamax of a
    function f, i.e. the integral of f(omega)*exp(j*omega*t) for any t in
    tscan (table) and for f given in input.
    df is another function giving the derivative of f. It can be
    identically zero.
    We use the function fourier_integral_any_sampling, and make an automatic
    convergence vs sampling (see N. Mounet PhD thesis, EPFL 5305)

    omegamin can be -np.inf, omegamax can be np.inf. They cannot be both infinite
    (do 2 separate semi-infinite integrals in such as case).
    If omegamin_ini is None, begins interval at omegamin (it assumes no singularity
    of f there), otherwise the interval will initially begin at omegamin_ini
    and convergence will be checked with respect to this value (getting it
    closer and closer to omegamin). Same for omegamax_ini.

    2 possible convergence criterions:
    - criterion begins by 'abs': make
    int_omegamin^omegamax |current_interpolation-interpolation_when_bisecting_all_intervals|<=tolabs,
    Then for the convergence vs minimum and maximum frequencies (when needed), a criterion
    in relative is still used (see just below). Criterion can then be:
        * absabs: difference in absolute value matters,
        * absreal: difference of real part matters,
        * absimag: difference of imaginary part matters,
        * absrealimag: difference of both real & imaginary part matters.
    - criterion begins by 'rel': make Fourier integral difference between the current interpolation
    and the interpolation when bisecting all intervals, smaller than tolrel for
    any t in tscan, in relative. Then, criterion can be:
        * relabs: difference in absolute value matters,
        * relreal: difference of real part matters,
        * relimag: difference of imaginary part matters,
        * relrealimag: difference of both real & imaginary part matters.

    2 possible kinds of interval bisection to get the next interpolation of f (i.e. the
    next step in the convergence process):
    - bisect=='log': bisection at middle point in log scale,
    - bisect=='lin': bisection at middle point in linear scale.

    The kind of interpolation ('linear' or 'pchip') is chosen to be the same on all
    intervals and given by optional argument interp (default='pchip').
    '''

    pchip_flag=interp.startswith('pchip');

    if (np.isinf(omegamin))or(np.isinf(omegamax)):
        if np.isinf(omegamin): iinf=0;flaginf=np.sign(omegamin);
        else: iinf=-1;flaginf=np.sign(omegamax);
    else: flaginf=0;

    flag_continue=True;

    if omegamin_ini==None: omegamin_ini=omegamin;
    if omegamax_ini==None: omegamax_ini=omegamax;
    omegai=np.array([omegamin_ini,omegamax_ini]);
    fi=np.array([f(omegamin_ini),f(omegamax_ini)],dtype=complex);

    if pchip_flag: di=np.ones(2,dtype=complex)*(fi[1]-fi[0])/(omegai[1]-omegai[0]);interp_int=0;
    else: di=np.zeros(2,dtype=complex);interp_int=1;

    # "refined" sampling (to check convergence)
    omegainew,finew,dinew,interp_value=bisect_and_calculate(f,0,omegai,fi,di,bisect_type=bisect,pchip_flag=pchip_flag);

    if criterion.startswith('abs'):
        # initialize table with integrated absolute error between one interpolation
        # and its refinement
        integrand=np.abs(finew-np.array([fi[0],interp_value,fi[1]]));
        if bisect.startswith('lin'): inte=np.array([int_Simpson(omegai[0],omegai[1],integrand)]);
        else: inte=np.array([int_Simpson(np.log(omegai[0]),np.log(omegai[1]),omegainew*integrand)]);
        #print integrand,inte;
        sumerr=np.sum(inte); # total error

    else:
        fint_i=np.zeros((len(omegai)-1,len(tscan)),dtype=complex);
        fintnew_i=np.zeros((len(omegai)-1,len(tscan)),dtype=complex);
        fint_i[0,:]=fourier_integral_any_sampling(tscan,omegai,fi,interp=interp_int,flaginf=0,di=di);
        fintnew_i[0,:]=fourier_integral_any_sampling(tscan,omegainew,finew,interp=interp_int,flaginf=0,di=dinew);
        fint=np.sum(fint_i,axis=0);
        fintnew=np.sum(fintnew_i,axis=0);

    imax=0; # index of lower frequency of interval to bisect

    # convergence loop over sampling between omegamin_ini and omegamax_ini
    print("Convergence on full frequency mesh");
    while flag_continue:

        #print "bisection between omega=",omegai[imax],"and omega=",omegai[imax+1];

        # replace interval beginning at omegai[imax] by its bisection (2 subintervals)
        omegai,fi,di,interp_value=bisect_and_calculate(f,imax,omegai,fi,di,bisect_type=bisect,pchip_flag=pchip_flag);

        # bisect the 2 new subintervals to get the refined interpolation
        i=pylab.mlab.find(omegainew==omegai[imax])[0];
        omegainew,finew,dinew,interp_value1=bisect_and_calculate(f,i,omegainew,finew,dinew,bisect_type=bisect,pchip_flag=pchip_flag);
        omegainew,finew,dinew,interp_value2=bisect_and_calculate(f,i+2,omegainew,finew,dinew,bisect_type=bisect,pchip_flag=pchip_flag);

        if criterion.startswith('abs'):

            # here, convergence is on the integral of absolute difference between f and
            # its interpolation
            sumerr -= inte[imax]; # first substract from the total error  the integral of the previous interval beginning at omegai[imax] (now splitted in 2)

            # some temporary tests
            #print i,imax,finew[i:i+3],np.array([fi[imax],interp_value1,fi[imax+1]]),finew[i+2:i+5],np.array([fi[imax+1],interp_value2,fi[imax+2]]);
            #tmp1=np.abs(finew[i:i+3]-np.array([fi[imax],interp_value1,fi[imax+1]]));
            #tmp2=np.abs(finew[i+2:i+5]-np.array([fi[imax+1],interp_value2,fi[imax+2]]));
            #print "tmps",tmp1,tmp2;
            #print "integrals Simpson=",int_Simpson(omegai[imax],omegai[imax+1],tmp1),int_Simpson(omegai[imax+1],omegai[imax+2],tmp2);
            #print "integrals from middle points only =",2/3.*np.abs(interp_value1-finew[i+1])*(omegai[imax+1]-omegai[imax]),2/3.*np.abs(interp_value2-finew[i+3])*(omegai[imax+2]-omegai[imax+1]);

            inte[imax]=2/3.*np.abs(interp_value1-finew[i+1])*(omegai[imax+1]-omegai[imax]);
            inte=add_element_in_array(inte,imax,2/3.*np.abs(interp_value2-finew[i+3])*(omegai[imax+2]-omegai[imax+1]));

            # add to total error
            sumerr += inte[imax]+inte[imax+1];
            # find place with largest error
            imax=np.argmax(inte);
            # flag for continuation of the convergence loop
            flag_continue=(sumerr>tolabs);

        elif criterion.startswith('rel'):

            # here, convergence is on the Fourier integral itself

            # add and recompute some 'partial' Fourier integrals for tables fint_i and fintnew_i
            if (imax>0)and(interp_int==0): fint_i[imax-1,:]=fourier_integral_any_sampling(tscan,omegai[imax-1:imax+1],fi[imax-1:imax+1],interp=interp_int,flaginf=0,di=di[imax-1:imax+1]);
            fint_i[imax,:]=fourier_integral_any_sampling(tscan,omegai[imax:imax+2],fi[imax:imax+2],interp=interp_int,flaginf=0,di=di[imax:imax+2]);
            fint_i=add_line_in_2d_array(fint_i,imax,fourier_integral_any_sampling(tscan,omegai[imax+1:imax+3],fi[imax+1:imax+3],interp=interp_int,flaginf=0,di=di[imax+1:imax+3]));
            if (imax+3<len(omegai))and(interp_int==0): fint_i[imax+2,:]=fourier_integral_any_sampling(tscan,omegai[imax+2:imax+4],fi[imax+2:imax+4],interp=interp_int,flaginf=0,di=di[imax+2:imax+4]);

            if (imax>0)and(interp_int==0): fintnew_i[imax-1,:]=fourier_integral_any_sampling(tscan,omegainew[i-2:i+1],finew[i-2:i+1],interp=interp_int,flaginf=0,di=dinew[i-2:i+1]);
            fintnew_i[imax,:]=fourier_integral_any_sampling(tscan,omegainew[i:i+3],finew[i:i+3],interp=interp_int,flaginf=0,di=dinew[i:i+3]);
            fintnew_i=add_line_in_2d_array(fintnew_i,imax,fourier_integral_any_sampling(tscan,omegainew[i+2:i+5],finew[i+2:i+5],interp=interp_int,flaginf=0,di=dinew[i+2:i+5]));
            if (imax+3<len(omegai))and(interp_int==0): fintnew_i[imax+2,:]=fourier_integral_any_sampling(tscan,omegainew[i+4:i+7],finew[i+4:i+7],interp=interp_int,flaginf=0,di=dinew[i+4:i+7]);

            # if needed, add correction for semi-infinite part
            if flaginf!=0: fasympt=fourier_integral_asympt_correction(flaginf,tscan,omegai[iinf],fi[iinf],df_end=df(omegai[iinf]));
            else: fasympt=0;

            fint=np.sum(fint_i,axis=0)+fasympt;
            fintnew=np.sum(fintnew_i,axis=0)+fasympt;

            # some temporary checks
            #fint2=fourier_integral_any_sampling(tscan,omegai,fi,interp=interp_int,flaginf=0,di=di);
            #fintnew2=fourier_integral_any_sampling(tscan,omegainew,finew,interp=interp_int,flaginf=0,di=dinew);
            #print np.max(np.abs(fint2-fint)/np.abs(fint2)),np.max(np.abs(fintnew2-fintnew)/np.abs(fintnew2));
            #print omegainew[i:i+3],omegainew[i+2:i+5];
            #print omegai[imax:imax+2],omegai[imax+1:imax+3];
            #print finew[i:i+3],finew[i+2:i+5];
            #print fi[imax:imax+2],fi[imax+1:imax+3];
            #print dinew[i:i+3],dinew[i+2:i+5];
            #print di[imax:imax+2],di[imax+1:imax+3];
            #it=np.argmax(np.abs(fint2-fint)/np.abs(fint2));itnew=np.argmax(np.abs(fintnew2-fintnew)/np.abs(fintnew2));
            #print it,itnew;
            #print np.abs(fint2-fint)/np.abs(fint2),np.abs(fintnew2-fintnew)/np.abs(fintnew2);
            #for iomega,omega in enumerate(omegai[:-1]): print fourier_integral_any_sampling([tscan[it]],omegai[iomega:iomega+2],fi[iomega:iomega+2],interp=interp_int,flaginf=0,di=di[iomega:iomega+2]);
            #print fint_i[:,it];

            # max relative error
            maxerr,iterr=max_diff_rel_complex_array(fintnew,fint,criterion=criterion[3:]);
            # look for place in omega scan where error (for this t) is max
            maxi,imax=max_diff_abs_complex_array(fintnew_i[:,iterr],fint_i[:,iterr],criterion=criterion[3:]);

            # flag for continuation of the convergence loop
            flag_continue=(maxerr>=tolrel);


        flag_continue=((flag_continue)or(len(omegai)<20)); # do at least 20 iterations


    if flaginf:
        # convergence loops over minimum and maximum angular frequency

        # minimum
        print("Convergence on minimum frequency");
        omegai,fi,di,omegainew,finew,dinew=fourier_int_convergence_boundary(f,df,omegamin,tscan,omegai,fi,
            di,omegainew,finew,dinew,0,pchip_flag=pchip_flag,tolrel=tolrel,criterion=criterion[3:])

        # maximum
        print("Convergence on maximum frequency");
        omegai,fi,di,omegainew,finew,dinew=fourier_int_convergence_boundary(f,df,omegamax,tscan,omegai,fi,
            di,omegainew,finew,dinew,1,pchip_flag=pchip_flag,tolrel=tolrel,criterion=criterion[3:])

        # final computations (with both converged and refined sampling)
        fint=fourier_integral_any_sampling(tscan,omegai,fi,df=df(omegai[iinf]),interp=interp_int,flaginf=flaginf,di=di);
        fintnew=fourier_integral_any_sampling(tscan,omegainew,finew,df=df(omegainew[iinf]),interp=interp_int,flaginf=flaginf,di=dinew);
    
    else:
        # final computations (with both converged and refined sampling)
        fint=fourier_integral_any_sampling(tscan,omegai,fi,df=0.,interp=interp_int,flaginf=flaginf,di=di);
        fintnew=fourier_integral_any_sampling(tscan,omegainew,finew,df=0.,interp=interp_int,flaginf=flaginf,di=dinew);
        

    maxerr,iterr=max_diff_rel_complex_array(fintnew,fint,criterion=criterion[3:]);
    print("Final max. relative difference between converged and refined meshes:",maxerr);#,tscan[iterr],tolrel;

    # final temporary checks (on the pchip derivatives and fi)
    #fint2=fourier_integral_any_sampling(tscan,omegai,f(omegai),df=df(omegai[iinf]),interp=interp_int,flaginf=flaginf,di=None);
    #fintnew2=fourier_integral_any_sampling(tscan,omegainew,f(omegainew),df=df(omegainew[iinf]),interp=interp_int,flaginf=flaginf,di=None);
    #print np.max(np.abs(fint2-fint)/np.abs(fint)),np.max(np.abs(fintnew2-fintnew)/np.abs(fintnew));
    #print np.max(np.abs(np.real(fintnew2-fint2)/np.real(fintnew2))),np.max(np.abs(np.imag(fintnew2-fint2)/np.imag(fintnew2)));

    return omegai,fi,di,fint,omegainew,finew,dinew,fintnew;


if __name__ == "__main__":

    # some tests
    import time as ti;
    from plot_lib import *
    from scipy.special import *

    # test correction functions for DFT method
    #t=10**np.arange(-8,4,0.001);
    #w=np.zeros(len(t));
    #alpha0=np.zeros(len(t));alpha1=np.zeros(len(t));
    #alpha2=np.zeros(len(t));alpha3=np.zeros(len(t));

    #t1=ti.clock();
    #for it,theta in enumerate(t):
    #   w[it],alpha0[it],alpha1[it],alpha2[it],alpha3[it]=correction_functions(theta,interp='linear',limit=1e-3);

    #t2=ti.clock();
    #print "Elapsed time:",t2-t1,"seconds";

    #pylab.semilogx(t,w,t,alpha0,t,alpha1,t,alpha2,t,alpha3);

    if False:
        # test Fourier integrals with 1/sqrt(omega)
        tscan=10.**np.arange(-10,-2.8,0.2);
        omegai=10.**np.arange(0,12.2,0.2);fi=1./np.sqrt(np.abs(omegai));

        tbeg1=ti.clock();
        fint1=fourier_integral_any_sampling(tscan,omegai,fi,interp=0,flaginf=1);
        tend1=ti.clock();
        print("Elapsed time for new Fourier integral (",len(omegai)," frequencies):",tend1-tbeg1,"seconds");

        # for DFT
        N=2e6;omegamin=2*np.pi*1e2;deltaomega=2*np.pi*1e3;
        omegai_lin=np.arange(omegamin,omegamin+N*deltaomega,deltaomega);
        fi_lin=1./np.sqrt(np.abs(omegai_lin));

        tbeg2=ti.clock();
        tscan_lin,fint2=fourier_integral_linear_sampling(tscan,omegai_lin,fi_lin,interp=0,flaginf=1);
        tend2=ti.clock();
        print("Elapsed time for DFT Fourier integral (",N,"frequencies):",tend2-tbeg2,"seconds");

        tbeg3=ti.clock();
        tscan_lin,fint3=fourier_integral_linear_sampling_naive(tscan,omegai_lin,fi_lin,flaginf=1);
        tend3=ti.clock();
        print("Elapsed time for naive DFT Fourier integral (",N,"frequencies):",tend3-tbeg3,"seconds");

        fig,ax=init_figure();
        plot(tscan,np.real(fint1),"New approach with %g frequencies (computation time=" % len(omegai) +str(tend1-tbeg1)+' s)','xk','Real part of the Fourier integral',ax,3,xlab='Time [s]');
        plot(tscan_lin,np.real(fint2),"DFT approach with %g frequencies (computation time=" % N +str(tend2-tbeg2)+' s)','--b','Real part of the Fourier integral',ax,3,xlab='Time [s]',plotevery=1);
        plot(tscan_lin,np.real(fint3),"Naive DFT approach with %g frequencies (computation time=" % N +str(tend3-tbeg3)+' s)',':g','Real part of the Fourier integral',ax,3,xlab='Time [s]',plotevery=1);
        plot(tscan,np.sqrt(np.pi/(2.*np.abs(tscan))),"Analytical result",'-r','Real part of the Fourier integral',ax,3,xlab='Time [s]');

        end_figure(fig,ax);


    #tscan=10.**np.arange(-10,-2.8,0.2);
    #f=(lambda omega: 1./np.sqrt(np.abs(omega)));
    #tolabs=1e2;tolrel=1e-3

    tscan=10.**np.arange(-13,-2.99,0.01);
    L=1;mu0=4e-7*np.pi;c=299792458;Z0=mu0*c;gamma=479.6;beta=np.sqrt(1.-1/gamma**2);
    mu1=1.;b=2e-3;eps0=1./(Z0*c);rho=10e-6;tau=0.8e-12;v=beta*c;
    eps1=(lambda omega:1-1j/(eps0*omega*rho*(1.+1j*omega*tau)));
    nu=(lambda omega: np.abs(omega)*np.sqrt(1-beta**2*eps1(omega)*mu1)/v);
    x1=(lambda omega: omega*b/(v*gamma));
    x2=(lambda omega: nu(omega)*b);
    kprime1e=(lambda z: -kve(0,z)-kve(1,z)/z);
    iprime1=(lambda z: iv(0,z)-iv(1,z)/z);
    diff1=(lambda omega: iprime1(x1(omega))/i1(x1(omega)) - kprime1e(x1(omega))/k1e(x1(omega)) );
    diff2=(lambda omega: gamma*nu(omega)*iprime1(x1(omega))/i1(x1(omega))
        - (mu1*omega/v)*kprime1e(x2(omega))/kve(1,x2(omega)) );
    diff3=(lambda omega: gamma*nu(omega)*iprime1(x1(omega))/i1(x1(omega))
        - (eps1(omega)*omega/v)*kprime1e(x2(omega))/kve(1,x2(omega)) );

    f=(lambda omega: 1j*L*Z0*beta*omega**2*k1(x1(omega))*x1(omega)**2*x2(omega)**2*nu(omega)*diff1(omega)*diff2(omega)
        /(4.*np.pi*gamma**3*v**2*i1(x1(omega))*((gamma*nu(omega)*x2(omega)-omega*x1(omega)/v)**2
        -(beta*x1(omega)*x2(omega))**2*diff2(omega)*diff3(omega))));
    tolabs=1e12;tolrel=1e-2;
    df=(lambda omega: 0.);

    omegai=10.**np.arange(0,14,0.02);
    fi=np.array([f(omega) for omega in omegai]);
    fig,ax=init_figure();
    plot(omegai,np.real(fi),'Re','g','Z',ax,3,xlab=" $ \omega $ [rad/s]");
    plot(omegai,np.imag(fi),'Im','--g','Z',ax,3,xlab=" $ \omega $ [rad/s]");
    end_figure(fig,ax);

    tbeg1=ti.clock();
    omegai1,fi1,di1,fint1,omegainew1,finew1,dinew1,fintnew1=fourier_int_converged(f,df,tscan,0.,np.inf,tolabs=tolabs,tolrel=tolrel,
        omegamin_ini=1e-2,omegamax_ini=1e13,criterion='absabs',bisect='log',interp='pchip');
    tend1=ti.clock();
    print("Elapsed time for converged (abs) Fourier integral (",len(omegai1)," frequencies):",tend1-tbeg1,"seconds");

    tbeg2=ti.clock();
    omegai2,fi2,di2,fint2,omegainew2,finew2,dinew2,fintnew2=fourier_int_converged(f,df,tscan,0.,np.inf,tolabs=tolabs,tolrel=tolrel,
        omegamin_ini=1e-2,omegamax_ini=1e13,criterion='relimag',bisect='log',interp='pchip');
    tend2=ti.clock();
    print("Elapsed time for converged (rel) Fourier integral (",len(omegai2)," frequencies):",tend2-tbeg2,"seconds");

    # frequency domain plots
    fig,ax=init_figure();
    plot(omegai1,fi1,"Abs. convergence, %g frequencies (computation time=" % len(omegai1) +str(tend1-tbeg1)+' s)','xb','f',ax,3,xlab=" $ \omega $ [rad/s]");
    plot(omegainew1,finew1,"Abs. convergence (refined mesh), %g frequencies (computation time=" % len(omegainew1) +str(tend1-tbeg1)+' s)','ob','f',ax,3,xlab=" $ \omega $ [rad/s]");
    plot(omegai2,fi2,"Rel. convergence, %g frequencies (computation time=" % len(omegai2) +str(tend2-tbeg2)+' s)','xk','f',ax,3,xlab=" $ \omega $ [rad/s]");
    plot(omegainew2,finew2,"Rel. convergence (refined mesh), %g frequencies (computation time=" % len(omegainew2) +str(tend2-tbeg2)+' s)','ok','f',ax,3,xlab=" $ \omega $ [rad/s]");
    end_figure(fig,ax);

    if True:
        # histograms (distribution of frequencies)
        xbar=10.**np.hstack(([-5.],np.arange(-4,14,0.1),[15.]));
        ybar1=np.zeros(len(xbar));ybar2=np.zeros(len(xbar));
        ybarnew1=np.zeros(len(xbar));ybarnew2=np.zeros(len(xbar));
        for ix,x in enumerate(xbar[:-1]):
            ybar1[ix]=count_between(omegai1,np.zeros(len(omegai1)),0.,x,xbar[ix+1]);
            ybarnew1[ix]=count_between(omegainew1,np.zeros(len(omegainew1)),0.,x,xbar[ix+1]);
            ybar2[ix]=count_between(omegai2,np.zeros(len(omegai2)),0.,x,xbar[ix+1]);
            ybarnew2[ix]=count_between(omegainew2,np.zeros(len(omegainew2)),0.,x,xbar[ix+1]);

        fig,ax=init_figure();
        ax.bar(xbar,ybar1,facecolor='b',edgecolor='b',label="Abs. convergence, %g frequencies (computation time=" % len(omegai1) +str(tend1-tbeg1)+' s)');
        ax.set_xscale('log');ax.set_xlabel(" $ \omega $ [rad/s]");
        end_figure(fig,ax);

        fig,ax=init_figure();
        ax.bar(xbar,ybarnew1,facecolor='b',edgecolor='b',label="Abs. convergence (refined mesh), %g frequencies (computation time=" % len(omegainew1) +str(tend1-tbeg1)+' s)');
        ax.set_xscale('log');ax.set_xlabel(" $ \omega $ [rad/s]");
        end_figure(fig,ax);

        fig,ax=init_figure();
        ax.bar(xbar,ybar2,facecolor='b',edgecolor='b',label="Rel. convergence, %g frequencies (computation time=" % len(omegai2) +str(tend2-tbeg2)+' s)');
        ax.set_xscale('log');ax.set_xlabel(" $ \omega $ [rad/s]");
        end_figure(fig,ax);

        fig,ax=init_figure();
        ax.bar(xbar,ybarnew2,facecolor='b',edgecolor='b',label="Rel. convergence (refined mesh), %g frequencies (computation time=" % len(omegainew2) +str(tend2-tbeg2)+' s)');
        ax.set_xscale('log');ax.set_xlabel(" $ \omega $ [rad/s]");
        end_figure(fig,ax);

    # time domain plots
    plottype=3; # loglog plots
    plottype=1; # semilogx plots
    fig,ax=init_figure();
    #plot(tscan,np.real(fint1),"Abs. convergence, %g frequencies (computation time=" % len(omegai1) +str(tend1-tbeg1)+' s)','xb','Real part of the Fourier integral',ax,plottype,xlab='Time [s]');
    #plot(tscan,np.real(fintnew1),"Abs. convergence (refined mesh), %g frequencies (computation time=" % len(omegainew1) +str(tend1-tbeg1)+' s)','-ob','Real part of the Fourier integral',ax,plottype,xlab='Time [s]');
    #plot(tscan,np.real(fint2),"Rel. convergence, %g frequencies (computation time=" % len(omegai2) +str(tend2-tbeg2)+' s)','xk','Real part of the Fourier integral',ax,plottype,xlab='Time [s]');
    #plot(tscan,np.real(fintnew2),"Rel. convergence (refined mesh), %g frequencies (computation time=" % len(omegainew2) +str(tend2-tbeg2)+' s)','ok','Real part of the Fourier integral',ax,plottype,xlab='Time [s]');
    #plot(tscan,np.sqrt(np.pi/(2.*np.abs(tscan))),"Analytical result",'-r','Real part of the Fourier integral',ax,3,xlab='Time [s]');
    plot(tscan,np.abs(np.imag(fint1))/np.pi,"Abs. convergence, %g frequencies (computation time=" % len(omegai1) +str(tend1-tbeg1)+' s)','-xb','Real part of the Fourier integral',ax,plottype,xlab='Time [s]');
    plot(tscan,np.abs(np.imag(fintnew1))/np.pi,"Abs. convergence (refined mesh), %g frequencies (computation time=" % len(omegainew1) +str(tend1-tbeg1)+' s)','-ob','Real part of the Fourier integral',ax,plottype,xlab='Time [s]');
    plot(tscan,np.abs(np.imag(fint2))/np.pi,"Rel. convergence, %g frequencies (computation time=" % len(omegai2) +str(tend2-tbeg2)+' s)','-xk','Real part of the Fourier integral',ax,plottype,xlab='Time [s]');
    plot(tscan,np.abs(np.imag(fintnew2))/np.pi,"Rel. convergence (refined mesh), %g frequencies (computation time=" % len(omegainew2) +str(tend2-tbeg2)+' s)','-ok','Real part of the Fourier integral',ax,plottype,xlab='Time [s]');
    end_figure(fig,ax);


    pylab.show()
