#!/usr/bin/python

# library for generalized elliptic functions and integrals

import sys
import subprocess
pymod=subprocess.check_output("echo $PYMOD",shell=True).strip().decode()
if pymod.startswith('local'):
    py_numpy=subprocess.check_output("echo $PY_NUMPY",shell=True).strip().decode()
    sys.path.insert(1,py_numpy)
    py_scipy=subprocess.check_output("echo $PY_SCIPY",shell=True).strip().decode()
    sys.path.insert(1,py_scipy)
import numpy as np
from scipy import special as sp
from scipy import sqrt,arcsin,sin,cos # generalized functions for complex numbers


def ellipkinc_gen(phi,m):
    ''' extension of scipy incomplete elliptic integral of the first kind for:
    - any complex phi,
    - any real m.
    Uses formulas from Abramowitz-Stegun (pp. 592-593)'''

    if (m>=0):

        if (m<=1):
            # standard case 0<=m<=1
            return ellipkinc_comp(phi,m);

        else:
            # case m>1
            theta=arcsin(np.sqrt(m)*sin(phi));
            #if (np.imag(theta)>0):
            #    # NOTE: choose always negative imaginary part (not in Abramowitz -
            #    # it's to be the same as Mathematica)
            #    theta=np.pi-theta;
            return m**(-0.5)*ellipkinc_comp(theta,1./m)

    else:
        # case m<0
        m2=-m/(1.-m);
        return (1.-m)**(-0.5)*(sp.ellipk(m2)-ellipkinc_comp(np.pi/2.-phi,m2));


def ellipeinc_gen(phi,m):
    '''extension of scipy incomplete elliptic integral of the second kind for:
    - any complex phi,
    - any real m.
    Uses formulas from Abramowitz-Stegun (pp. 592-593)'''

    if (m>=0):

        if (m<=1):
            # standard case 0<=m<=1
            return ellipeinc_comp(phi,m);

        else:
            # case m>1
            u=ellipkinc_gen(phi,m);
            m2=np.sqrt(m);
            sn,cn,dn,ph=ellipj_comp(u*m2,1./m);
            return m2*ellipeinc_comp(ph,1./m)-(m-1.)*u;

    else:
        # case m<0
        m2=-m/(1.-m);
        m3=np.sqrt(1.-m);
        u=ellipkinc_gen(phi,m);
        sn,cn,dn,ph=ellipj_comp(u*m3,m2);
        cd=cn/dn;
        # NOTE: mistake (according to Mathematica) in Abramowitz: last "/m3" not there
        return m3*ellipeinc_comp(ph,m2)+m*sn*cd/m3;


def ellipj_comp(u,m):
    '''extension of scipy Jacobian elliptic functions for:
    - any complex u
    - still 0<=m<=1.
    Uses formulas from Abramowitz-Stegun (pp. 575 - section 16.21)'''

    x=np.real(u);y=np.imag(u);

    if (m<0)or(m>1):
        print(" ellipj_comp does not support parameter m outside [0,1] !");sys.exit()

    if (y==0):
        # standard case (real u)
        return ellipj_correct(u,m);

    else:
        # case of complex u
        m1=1.-m;
        s,c,d,ph=ellipj_correct(x,m);
        s1,c1,d1,ph1=ellipj_correct(y,m1);

        denom=c1**2+m*s**2*s1**2;
        sn=(s*d1+1j*c*d*s1*c1)/denom;
        cn=(c*c1-1j*s*d*s1*d1)/denom;
        dn=(d*c1*d1-1j*m*s*c*s1)/denom;
        ph=arcsin(sn);

        if np.abs((cos(ph)+cn)/cn)<1e-10: ph=np.pi-ph;
        if np.abs((cos(ph)-cn)/cn)>1e-10: print(" Warning in ellipj_comp: cn different from cos(ph)",cn,cos(ph));

        return sn,cn,dn,ph;


def ellipj_correct(u,m):
    '''correction of scipy Jacobian elliptic functions for real u
    and 0<=m<=1.
    Uses formulas from Abramowitz-Stegun (pp. 569 - section 16.1.5)'''

    if (m<0)or(m>1):
        print(" ellipj_correct does not support parameter m outside [0,1] !");sys.exit()

    if (np.imag(u)!=0):
        print(" ellipj_correct does not support complex u ! Use ellipj_comp instead.");sys.exit()

    sn,cn,dn,ph=sp.ellipj(u,m);
    # for some reason dn is sometimes wrong ! so we recompute it:
    dn=sqrt(1.-m*sin(ph)**2);

    return sn,cn,dn,ph;


def ellipkinc_comp(phi,m):
    '''extension of scipy incomplete elliptic integral of the first kind for:
     - any complex phi,
     - still 0<=m<=1.
    Uses formulas from Abramowitz-Stegun (pp. 592-593)'''

    if (m<0)or(m>1):
        print(" ellipkinc_comp does not support parameter m outside [0,1]; use ellipkinc_gen instead.");sys.exit()


    if (np.imag(phi)==0):
        # standard case (real phi)
        return sp.ellipkinc(phi,m);

    elif (np.real(phi)==0):
        # case of purely imaginary phi
        m1=1.-m;
        theta=np.arctan(np.sinh(np.imag(phi)));
        return 1j*sp.ellipkinc(theta,m1);

    else:
        # case of complex phi
        m1=1.-m;
        lambd,mu=solve_cot2lambda(np.real(phi),np.imag(phi),m);

        return sp.ellipkinc(lambd,m)+1j*sp.ellipkinc(mu,m1);


def ellipeinc_comp(phi,m):
    '''extension of scipy incomplete elliptic integral of the second kind for:
     - any complex phi,
     - still 0<=m<=1.
    Uses formulas from Abramowitz-Stegun (pp. 592-593)'''

    if (m<0)or(m>1):
        print(" ellipeinc_comp does not support parameter m outside [0,1]; use ellipeinc_gen instead.");sys.exit()


    if (np.imag(phi)==0):
        # standard case (real phi)
        return sp.ellipeinc(phi,m);

    elif (np.real(phi)==0):
        # case of purely imaginary phi
        m1=1.-m;
        sinhphi=np.sinh(np.imag(phi));
        theta=np.arctan(sinhphi);
        return -1j*sp.ellipeinc(theta,m1)+1j*sp.ellipkinc(theta,m1)+1j*sinhphi*np.sqrt(1.-m1*np.sin(theta)**2);

    else:
        # case of complex phi
        m1=1.-m;
        lambd,mu=solve_cot2lambda(np.real(phi),np.imag(phi),m);

        coslambd=np.cos(lambd);sinlambd=np.sin(lambd);
        cos2lambd=coslambd**2;sin2lambd=sinlambd**2;
        cosmu=np.cos(mu);sinmu=np.sin(mu);
        cos2mu=cosmu**2;sin2mu=sinmu**2;

        b1=m*sinlambd*coslambd*sin2mu*np.sqrt(1.-m*sin2lambd);
        b2=(1.-m*sin2lambd)*np.sqrt(1.-m1*sin2mu)*sinmu*cosmu;
        b3=cos2mu+m*sin2lambd*sin2mu;

        return sp.ellipeinc(lambd,m)-1j*sp.ellipeinc(mu,m1)+1j*sp.ellipkinc(mu,m1)+(b1+1j*b2)/b3;


def solve_cot2lambda(phi,psi,m):
    '''solve second order equation to get cot^2(lambda), then lambda and mu
    useful to compute complex amplitude case of incomplete elliptic integrals
    (Abram-Stegun top of p. 593, 17.4.11 and 17.4.12)

    NOTE: some modifications w.r.t. Abramowitz, to get same as Mathematica:
     - sign of mu is the same sign as psi,
     - case phi=pi/2 developped separately (I did a Taylor expansion to solve the equations)'''

    if (phi==0): print("Pb in solve_cot2lambda: phi cannot be zero (purely imaginary amplitude of elliptic integral)!");sys.exit();

    if (np.abs(phi-np.pi/2)<1e-10):
        # specific case phi=pi/2 (some infinities appear); see side calculations by N. Mounet
        x=m*np.cosh(psi)**2-1;
        if (x>0):
            # positive root of polynomial is nonzero
            lambd=np.arctan(1./np.sqrt(x));
            mu=np.sign(psi)*np.pi/2.;
        elif (x<0):
            # positive root of polynomial is zero
            lambd=np.pi/2.;
            mu=np.arctan(np.sinh(psi)/np.sqrt(-x));
        else:
            # case x=0 (roots both equal to 0) (not sure about this case...)
            lambd=np.pi/2.;
            mu=np.sign(psi)*np.pi/2.;

    else:
        # general case
        tan2phi=np.tan(phi)**2;
        cot2phi=1./tan2phi;
        csc2phi=1./(np.sin(phi))**2;
        sinh2psi=np.sinh(psi)**2;
        m1=1.-m;
        p=[1.,-(cot2phi+m*sinh2psi*csc2phi-m1),-m1*cot2phi];

        x=np.roots(p);cot2lambda=np.max(x);

        if not(cot2lambda>=0): print("Pb in solve_cot2lambda: no positive real root !",p);sys.exit();

        lambd=np.arctan((1./np.sqrt(cot2lambda)));
        mu=np.sign(psi)*np.arctan(np.sqrt((tan2phi*cot2lambda-1.)/m));

    return lambd,mu;
