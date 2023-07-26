#!/usr/bin/python

# library to try to find all roots of complex functions in the complex plane
# (in a given region)

import sys
import subprocess
pymod=subprocess.check_output("echo $PYMOD",shell=True).strip().decode()
if pymod.startswith('local'):
    py_numpy=subprocess.check_output("echo $PY_NUMPY",shell=True).strip().decode()
    sys.path.insert(1,py_numpy)
    py_matpl=subprocess.check_output("echo $PY_MATPL",shell=True).strip().decode()
    sys.path.insert(1,py_matpl)
    py_scipy=subprocess.check_output("echo $PY_SCIPY",shell=True).strip().decode()
    sys.path.insert(1,py_scipy)
from string import *
import numpy as np
from numpy import fft
import re,os
from tables_lib import *


def eqplane(A,B,C):

    ''' Find the equation of a plane from the coordinates of 3 points A, B and C
    in 3D. Outputs are the 4 coefficients a,b,c,d such that the equation of
    the plane is ax+by+cz=d for points of coordinates (x,y,z)'''

    a=np.cross(B-A,C-A);
    d=np.sum(a*A);

    return a[0],a[1],a[2],d;


def intersect_two_planes_z0(p1,p2):

    ''' Given 2 equations of planes in 3D p1 and p2 (each being an array of 4
    coefficents a,b,c,d such that eq. is ax+by+cz=d), find their intersection
    with the plane z=0 (if it's a single point, otherwise return inf if no
    such intersection can exist, or nan if it's more than a single point)
    '''

    det=p2[1]*p1[0]-p2[0]*p1[1];

    if (det==0.):

        if (p2[1]*p1[3]-p2[3]*p1[1])==0.: return np.nan,np.nan; # a line or a plane is solution

        else: return np.inf,np.inf; # no solution


    else:

        x= (p2[1]*p1[3]-p2[3]*p1[1])/det;
        y= (-p2[0]*p1[3]+p2[3]*p1[0])/det;

        return x,y;


def potential_roots_complex(f,real_min,real_max,imag_min,imag_max,npts=1000,meshkind='lin',offset_log=0.):

    ''' Find the potential roots in the complex plane of a function f (from C to C),
    inside the square of the complex plane delimited by real_min / real_max for the
    real part, and imag_min / imag_max for the imaginary part.
    "potential roots" mean roots of the linear interpolation of the function
    over a grid defined on this region, with npts points in each dimension.
    meshkind: indicates the kind of mesh, which can be can be either 'lin'
    (equidistant mesh) or 'log' (logarithmic mesh w.r.t. a constant offset given
    in offset_log - complex number).
    f must take a complex argument and returns also a complex.
    Outputs a list of potential roots.
    '''


    if meshkind.startswith('lin'):
        real_range=np.linspace(real_min,real_max,npts+1);
        imag_range=np.linspace(imag_min,imag_max,npts+1);

    elif meshkind.startswith('log'):

        if (real_min-np.real(offset_log))>0:
            # case when all range is above offset_log
            real_min_log=np.log10(real_min-np.real(offset_log));
            real_max_log=np.log10(real_max-np.real(offset_log));
            real_range=np.logspace(real_min_log,real_max_log,npts+1)+np.real(offset_log);
        elif (real_max-np.real(offset_log))<0:
            # case when all range is below offset_log
            real_min_log=np.log10(np.real(offset_log)-real_min);
            real_max_log=np.log10(np.real(offset_log)-real_max);
            real_range=-np.logspace(real_max_log,real_min_log,npts+1)+np.real(offset_log);
            real_range=real_range[::-1];
        else:
            # case when part of range is below offset_log and the rest is above offset_log
            real_min_log=np.log10(np.real(offset_log)-real_min);
            real_max_log=np.log10(real_max-np.real(offset_log));
            real_range1=-np.logspace(real_min_log-min(npts/4,6),real_min_log,npts/2+1);
            real_range2=np.logspace(real_max_log-min(npts/4,6),real_max_log,npts/2+1);
            real_range=np.hstack((real_range1[::-1],real_range2))+np.real(offset_log);


        #print "param:",imag_max,imag_min,offset_log;
        if (imag_min-np.imag(offset_log))>0:
            # case when all range is above offset_log
            imag_min_log=np.log10(imag_min-np.imag(offset_log));
            imag_max_log=np.log10(imag_max-np.imag(offset_log));
            imag_range=np.logspace(imag_min_log,imag_max_log,npts+1)+np.imag(offset_log);
        elif (imag_max-np.imag(offset_log))<0:
            # case when all range is below offset_log
            imag_min_log=np.log10(np.imag(offset_log)-imag_min);
            imag_max_log=np.log10(np.imag(offset_log)-imag_max);
            imag_range=-np.logspace(imag_max_log,imag_min_log,npts+1)+np.imag(offset_log);
            imag_range=imag_range[::-1];
        else:
            # case when part of range is below offset_log and the rest is above offset_log
            imag_min_log=np.log10(np.imag(offset_log)-imag_min);
            imag_max_log=np.log10(imag_max-np.imag(offset_log));
            imag_range1=-np.logspace(imag_min_log-min(npts/4,6),imag_min_log,npts/2+1);
            imag_range2=np.logspace(imag_max_log-min(npts/4,6),imag_max_log,npts/2+1);
            imag_range=np.hstack((imag_range1[::-1],imag_range2))+np.imag(offset_log);
            #print "lists:",imag_max_log,imag_min_log,imag_range1,imag_range2,imag_range;

    #print "real range:",real_range;
    #print "imag range:",imag_range;
    list_roots=[];

    ftable=np.zeros((len(real_range),len(imag_range)),dtype=complex);

    for i_im,im in enumerate(imag_range):

        for i_re,re in enumerate(real_range):

            ftable[i_re,i_im]=f(re+1j*im);

    # find roots using linear interpolation
    for i_im,im in enumerate(imag_range[:-1]):

        for i_re,re in enumerate(real_range[:-1]):

            rep1=real_range[i_re+1];imp1=imag_range[i_im+1];

            # points of the mesh cell on the surface defined by the real part of f
            Areal=np.array([re,im,ftable[i_re,i_im].real]);
            Breal=np.array([rep1,im,ftable[i_re+1,i_im].real])
            Creal=np.array([re,imp1,ftable[i_re,i_im+1].real])
            Dreal=np.array([rep1,imp1,ftable[i_re+1,i_im+1].real])
            # points of the mesh cell on the surface defined by the imag part of f
            Aimag=np.array([re,im,ftable[i_re,i_im].imag]);
            Bimag=np.array([rep1,im,ftable[i_re+1,i_im].imag])
            Cimag=np.array([re,imp1,ftable[i_re,i_im+1].imag])
            Dimag=np.array([rep1,imp1,ftable[i_re+1,i_im+1].imag])

            # generate the plane containing the lower-left triangle interpolating
            # the function f on half of the mesh cell
            # first real part
            p_real=eqplane(Areal,Breal,Creal);
            # then imaginary part
            p_imag=eqplane(Aimag,Bimag,Cimag);
            # intersect the two planes with the z=0 plane and check if
            # the resulting point is inside the mesh cell (note: it can be
            # in the other triangle of the cell)
            x,y=intersect_two_planes_z0(p_real,p_imag);
            # x & y are real, but their type is complex -> we extract the real part
            if np.isnan(x): list_roots.append((re.real,im.real));
            elif ((x>=re)and(x<=rep1))and((y>=im)and(y<=imp1)): list_roots.append((x.real,y.real));

            # same for upper-right triangle
            p_real=eqplane(Breal,Creal,Dreal);
            p_imag=eqplane(Bimag,Cimag,Dimag);
            x,y=intersect_two_planes_z0(p_real,p_imag);
            # x & y are real, but their type is complex -> we extract the real part
            if np.isnan(x): list_roots.append((rep1.real,imp1.real));
            elif ((x>=re)and(x<=rep1))and((y>=im)and(y<=imp1)): list_roots.append((x.real,y.real));

    # simpler version, that only uses sign of each mesh points (but find less roots)
#    for i_im,im in enumerate(imag_range[:-1]):
#
#       for i_re,re in enumerate(real_range[:-1]):
#
#           # points of the mesh cell on the surface defined by the real part of f
#           Areal=ftable[i_re,i_im].real;
#           Breal=ftable[i_re+1,i_im].real;
#           Creal=ftable[i_re,i_im+1].real;
#           Dreal=ftable[i_re+1,i_im+1].real;
#           # points of the mesh cell on the surface defined by the imag part of f
#           Aimag=ftable[i_re,i_im].imag;
#           Bimag=ftable[i_re+1,i_im].imag;
#           Cimag=ftable[i_re,i_im+1].imag;
#           Dimag=ftable[i_re+1,i_im+1].imag;
#
#           flag=True;
#           if np.abs(np.sign(Areal)+np.sign(Breal)+np.sign(Creal)+np.sign(Dreal))==4:
#               flag=False;
#           if np.abs(np.sign(Aimag)+np.sign(Bimag)+np.sign(Cimag)+np.sign(Dimag))==4:
#               flag=False;
#           # flag is True if some signs are not equal, both for real and imag. parts.
#           # then the center of the square is chosen as potential root
#           if flag: list_roots.append(((re+real_range[i_re+1])/2.,(im+real_range[i_im+1])/2.));

    #return list_roots,real_range,imag_range;
    return list_roots;


def roots_complex_list(f,list_pot_roots,tolf=1e-12,tolx=1e-14):

    ''' Find the roots in the complex plane of a function f (from C to C),
    taking as initial estimates the 'potential roots' in list_pot_roots.
    f must take a complex argument and returns also a complex.
    Outputs a list of roots.
    tolf and tolx are the tolerances for resp. the function value and
    the minimum acceptable difference between 2 roots.
    '''

    from scipy.optimize import fsolve,fmin_powell;#,root

    # function giving absolute value of f, and written as a 2D function (of
    # the real and imag. parts)
    #f_2d=(lambda x: (f(x[0]+1j*x[1]).real,f(x[0]+1j*x[1]).imag));

    list_roots=[];
    
    for x0 in list_pot_roots:
        # function giving absolute value of f, and written as a 2D function (of
        # the real and imag. parts)
        #f_2d=(lambda x: (f((1+x[0])*x0[0]+1j*(1+x[1])*x0[1]).real,f((1+x[0])*x0[0]+1j*(1+x[1])*x0[1]).imag));
        f_2d=(lambda x: (f(x[0]+1j*x[1]).real,f(x[0]+1j*x[1]).imag));

        xsol=fsolve(f_2d, x0);#, warning=False);
        #sol=root(f_2d, (0.,0.), method='hybr');xsol=sol.x;
        #xsol_comp=(1+xsol[0])*x0[0]+1j*(1+xsol[1])*x0[1];
        xsol_comp=xsol[0]+1j*xsol[1];#print "fsolve:",x0,xsol_comp,np.abs(f(xsol_comp));

        #if (len(list_roots)>0): print np.min(np.abs(xsol_comp-np.array(list_roots)));
        if (np.abs(f(xsol_comp))<tolf)and((len(list_roots)==0)or(np.min(np.abs(xsol_comp-np.array(list_roots)))>tolx)):
            list_roots.append(xsol_comp);

        # make another trial by maximizing first 1./abs(f) (actually,
        # minimizing -1/abs(f) ) (USELESS apparently)
        #invf=(lambda x: -1./np.abs(f(x[0]+1j*x[1])));
        #x1=fmin_powell(invf,x0);
        #xsol=fsolve(f_2d, x1);
        #xsol_comp=xsol[0]+1j*xsol[1];print "fmin:",x0,x1,xsol_comp;
        #if (np.abs(f(xsol_comp))<tolf)and((len(list_roots)==0)or(np.min(np.abs(xsol_comp-np.array(list_roots)))>tolx)):
        #    list_roots.append(xsol_comp);

    return np.array(list_roots);


def roots_complex(f,real_min,real_max,imag_min,imag_max,npts=1000,
        meshkind='lin',tolf=1e-12,tolx=1e-14):

    ''' function wrapping the two previous functions
    '''

    list_pot=potential_roots_complex(f,real_min,real_max,imag_min,imag_max,npts=npts,meshkind=meshkind);

    list_res=roots_complex_list(f,list_pot,tolf=tolf,tolx=tolx);

    return list_res;



if __name__ == "__main__":


    from scipy.special import j0,ive,kve,iv,kv
    import pylab
    from plot_lib import *

    if False:
        # some tests
        A=np.array([4.,0.,0.]);
        B=np.array([4.,2.,3.]);
        C=np.array([4.,0.,1.]);
        a,b,c,d=eqplane(A,B,C);
        print(a,b,c,d);

        A=np.array([2.,3.,0.]);
        B=np.array([4.,2.,3.]);
        C=np.array([5.,0.,1.]);
        D=np.array([1.,2.,3.]);
        E=np.array([-4.,5.,1.]);
        p1=eqplane(A,B,C);
        p2=eqplane(A,D,E);
        print(p1,p2,intersect_two_planes_z0(p1,p2));
        #OK (tested different combinations)

    if True:
        f=(lambda z: np.cos(z))
        f=(lambda z: np.cos(z.real)+1j*j0(z.imag))
        if True:
            list1=potential_roots_complex(f,-10,10,-10,10,npts=20,meshkind='lin');
            #print list1;
            #print f(np.array(list1)[:,0]+1j*np.array(list1)[:,1]);
            list2=roots_complex_list(f,list1,tolf=1e-10);
            #print list2;
            #print f(list2);
            # seems to work ! warning: tolf should not be too small

            # 2D plot
            npts=100;
            real_range=np.linspace(-10.,10.,npts+1);
            imag_range=np.linspace(-10.,10.,npts+1);
            ftable_real=np.zeros((npts+1,npts+1));
            ftable_imag=np.zeros((npts+1,npts+1));
            invftable=np.zeros((npts+1,npts+1));

            for ir,r in enumerate(real_range):
                for im,m in enumerate(imag_range):
                    ftable_real[im,ir]=f(r+1j*m).real;
                    ftable_imag[im,ir]=f(r+1j*m).imag;
                    invftable[im,ir]=1./np.abs(f(r+1j*m));

            fig,ax=init_figure();
            plot2D(ftable_real,-10,10,-10,10,'Re','Im','',ax,colorlabel='');
            plot(np.array(list1)[:,0],np.array(list1)[:,1],'potential zeros','ow','Im',ax,0,xlab='Re',ms=15);
            plot(np.real(list2),np.imag(list2),'zeros','xk','Im',ax,0,xlab='Re',ms=15);
            ax.set_xlim([-10,10]);ax.set_ylim([-10,10]);
            end_figure(fig,ax);

            fig,ax=init_figure();
            plot2D(ftable_imag,-10,10,-10,10,'Re','Im','',ax,colorlabel='');
            plot(np.array(list1)[:,0],np.array(list1)[:,1],'potential zeros','ow','Im',ax,0,xlab='Re',ms=15);
            plot(np.real(list2),np.imag(list2),'zeros','xk','Im',ax,0,xlab='Re',ms=15);
            ax.set_xlim([-10,10]);ax.set_ylim([-10,10]);
            end_figure(fig,ax);

            fig,ax=init_figure();
            plot2D(invftable,-10,10,-10,10,'Re','Im','',ax,colorlabel='');
            plot(np.array(list1)[:,0],np.array(list1)[:,1],'potential zeros','ow','Im',ax,0,xlab='Re',ms=15);
            plot(np.real(list2),np.imag(list2),'zeros','xk','Im',ax,0,xlab='Re',ms=15);
            ax.set_xlim([-10,10]);ax.set_ylim([-10,10]);
            end_figure(fig,ax);

            pylab.show()

        list3=roots_complex(f,-10,10,-10,10,npts=100,tolf=1e-10);
        ind=np.argsort(list3.imag);
        #print list3[ind];
        #print f(list3[ind])
        list_fin=np.conjugate(-1j*np.sort_complex(1j*np.conjugate(list3)));
        print(list_fin);


    if False:
        gamma=479.6;sigmaDC=2e5;b=2e-3;tauAC=0.;#tauAC=4.2e-12;
        beta=np.sqrt(1.-1./gamma**2);
        c=299792458;mu0=4e-7*np.pi;eps0=1./(mu0*c**2);
        sigma=(lambda omega: sigmaDC/(1+1j*omega*tauAC));
        v0=(lambda omega: omega/(beta*gamma*c));
        v1=(lambda omega: np.sqrt(v0(omega)**2+1j*sigma(omega)*mu0*omega));
        # equation to obtain resonance at high freq. for long. wall impedance
        # (condition of synchronicity between k of longitudinal waveguide mode &
        # beam wave number omega/v, for m=0)
        #print "m=0";res_eq=(lambda omega: ive(1,v0(omega)*b)/ive(0,v0(omega)*b)+v0(omega)*(1.-1j*sigma(omega)/(omega*eps0))/v1(omega) *(kve(1,v1(omega)*b)/kve(0,v1(omega)*b)));
        # equation to obtain resonance at high freq. for transverse wall impedance
        # (condition of synchronicity between k of transverse waveguide mode &
        # beam wave number omega/v, for m=1)
        ip1overi1=(lambda z: ive(0,z)/ive(1,z)-1./z);
        kp1overk1=(lambda z: -kve(0,z)/kve(1,z)-1./z);
        m11=(lambda omega: (ip1overi1(v0(omega)*b)/v0(omega)-kp1overk1(v1(omega)*b)/v1(omega)));
        m22=(lambda omega: (ip1overi1(v0(omega)*b)/v0(omega)-(1.-1j*sigma(omega)/(omega*eps0))*kp1overk1(v1(omega)*b)/v1(omega)));
        m21=(lambda omega: (1/v0(omega)**2-(1.-1j*sigma(omega)/(omega*eps0))/(v1(omega)**2))/b);
        m12=(lambda omega: (1/v0(omega)**2-1/v1(omega)**2)/b);
        print("m=1");res_eq=(lambda omega: m11(omega)*m22(omega) - m12(omega)*m21(omega));
        # after normalization
        f=(lambda fTHz:res_eq(fTHz*1e12*2*np.pi)/res_eq(1e12*2*np.pi));

        maximag=2;
        list1=potential_roots_complex(f,0.1,10,-maximag,maximag,npts=100,meshkind='lin');
        print(list1);
        if (len(list1)==0): list1.append((1.,0.));
        list2=roots_complex_list(f,list1,tolf=1e-8,tolx=1e-8);
        print(list2)
        # 2D plot
        npts=200;
        real_range=np.linspace(0.1,10.,npts+1);
        imag_range=np.linspace(-maximag,maximag,npts+1);
        ftable_real=np.zeros((npts+1,npts+1));
        ftable_imag=np.zeros((npts+1,npts+1));
        invftable=np.zeros((npts+1,npts+1));

        for ir,r in enumerate(real_range):
            for im,m in enumerate(imag_range):
                ftable_real[im,ir]=f(r+1j*m).real;
                ftable_imag[im,ir]=f(r+1j*m).imag;
                invftable[im,ir]=1./np.abs(f(r+1j*m));

        fig,ax=init_figure();
        plot2D(ftable_real,0.1,10,-maximag,maximag,'Re','Im','',ax,colorlabel='');
        plot(np.array(list1)[:,0],np.array(list1)[:,1],'potential zeros','ow','Im',ax,0,xlab='Re',ms=15);
        plot(np.real(list2),np.imag(list2),'zeros','xk','Im',ax,0,xlab='Re',ms=15);
        ax.set_xlim([0.1,10]);ax.set_ylim([-maximag,maximag]);
        end_figure(fig,ax);
        fig,ax=init_figure();
        plot2D(ftable_imag,0.1,10,-maximag,maximag,'Re','Im','',ax,colorlabel='');
        plot(np.array(list1)[:,0],np.array(list1)[:,1],'potential zeros','ow','Im',ax,0,xlab='Re',ms=15);
        plot(np.real(list2),np.imag(list2),'zeros','xk','Im',ax,0,xlab='Re',ms=15);
        ax.set_xlim([0.1,10]);ax.set_ylim([-maximag,maximag]);
        end_figure(fig,ax);
        fig,ax=init_figure();
        plot2D(invftable,0.1,10,-maximag,maximag,'Re','Im','',ax,colorlabel='');
        plot(np.array(list1)[:,0],np.array(list1)[:,1],'potential zeros','ow','Im',ax,0,xlab='Re',ms=15);
        plot(np.real(list2),np.imag(list2),'zeros','xk','Im',ax,0,xlab='Re',ms=15);
        ax.set_xlim([0.1,10]);ax.set_ylim([-maximag,maximag]);
        end_figure(fig,ax);

        pylab.show()

        # results: with tauAC=0
        # - m=0 (long.): f[THz]=1.74474324+1.01311935j
        # - m=1 (trans.): f[THz]=1.74441343+1.01504579j
        # results: with tauAC=4.2 ps
        # - m=0 (long.): f[THz]= 0.74513117+0.00950448j
        # - m=1 (trans.): f[THz]=7.44436724e-01 +9.52548055e-03j

    if False:
        # same as above but change a bit the perspective: solve in k and use omega=Re[k]*beta*c
        # BUT IT DOES NOT WORK...
        gamma=479.6;sigmaDC=2e5;b=2e-3;tauAC=0.;tauAC=4.2e-12;
        beta=np.sqrt(1.-1./gamma**2);
        c=299792458;mu0=4e-7*np.pi;eps0=1./(mu0*c**2);v=beta*c;
        sigma=(lambda omega: sigmaDC/(1+1j*omega*tauAC));
        v0=(lambda k: np.sqrt(k**2-(k.real*beta)**2));
        v1=(lambda k: np.sqrt(v0(k)**2+1j*sigma(k.real*v)*mu0*k.real*v));
        # equation to obtain resonance at high freq. for long. wall impedance
        # (condition of synchronicity between Re[k] of longitudinal waveguide mode &
        # beam wave number omega/v, for m=0)
        #print "m=0";res_eq=(lambda k: ive(1,v0(k)*b)/ive(0,v0(k)*b)+v0(k)*(1.-1j*sigma(k.real*v)/(k.real*v*eps0))/v1(k) *(kve(1,v1(k)*b)/kve(0,v1(k)*b)));
        # equation to obtain resonance at high freq. for transverse wall impedance
        # (condition of synchronicity between Re[k] of transverse waveguide mode &
        # beam wave number omega/v, for m=1)
        ip1overi1=(lambda z: ive(0,z)/ive(1,z)-1./z);
        kp1overk1=(lambda z: -kve(0,z)/kve(1,z)-1./z);
        m11=(lambda k: (ip1overi1(v0(k)*b)/v0(k)-kp1overk1(v1(k)*b)/v1(k)));
        m22=(lambda k: (ip1overi1(v0(k)*b)/v0(k)-(1.-1j*sigma(k.real*v)/(k.real*v*eps0))*kp1overk1(v1(k)*b)/v1(k)));
        m21=(lambda k: (1/v0(k)**2-(1.-1j*sigma(k.real*v)/(k.real*v*eps0))/(v1(k)**2))/b);
        m12=(lambda k: (1/v0(k)**2-1/v1(k)**2)/b);
        print("m=1");res_eq=(lambda k: m11(k)*m22(k) - m12(k)*m21(k));
        # after normalization
        f=(lambda fTHz:res_eq(fTHz*1e12*2*np.pi/v)/res_eq(1e12*2*np.pi/v));

        # 2D root solving: does not work here (maybe no zero ?)
        maximag=0.1;maxreal=2.;npts=200;
        list1=potential_roots_complex(f,0.1,maxreal,-maximag,maximag,npts=npts,meshkind='lin');
        #list1,r_ra,i_ra=potential_roots_complex(f,0.1,maxreal,-maximag,maximag,npts=200,meshkind='lin');
        #gridR,gridI=np.meshgrid(r_ra,i_ra)
        print(list1);
        real_range=np.linspace(0.1,maxreal,npts+1);
        imag_range=np.linspace(-maximag,maximag,npts+1);
        if (len(list1)==0):
            for r in real_range:
                for i in imag_range: list1.append((r,i));
        list2=roots_complex_list(f,list1,tolf=1e-6,tolx=1e-8);
        print(list2)
        # 2D plot
        npts=200;
        real_range=np.linspace(0.1,maxreal,npts+1);
        imag_range=np.linspace(-maximag,maximag,npts+1);
        gridR,gridI=np.meshgrid(real_range,imag_range);
        ftable_real=np.zeros((npts+1,npts+1));
        ftable_imag=np.zeros((npts+1,npts+1));
        invftable=np.zeros((npts+1,npts+1));

        for ir,r in enumerate(real_range):
            for im,m in enumerate(imag_range):
                ftable_real[im,ir]=f(r+1j*m).real;
                ftable_imag[im,ir]=f(r+1j*m).imag;
                invftable[im,ir]=1./np.abs(f(r+1j*m));

        fig,ax=init_figure();
        plot2D(ftable_real,0.1,maxreal,-maximag,maximag,'Re','Im','',ax,colorlabel='',colorlim=[-1,1]);
        #plot(gridR,gridI,'','.k','Im',ax,0,xlab='Re',ms=2);
        plot(np.array(list1)[:,0],np.array(list1)[:,1],'potential zeros','ow','Im',ax,0,xlab='Re',ms=5);
        plot(np.real(list2),np.imag(list2),'zeros','xk','Im',ax,0,xlab='Re',ms=15);
        ax.set_xlim([0.1,maxreal]);ax.set_ylim([-maximag,maximag]);
        end_figure(fig,ax);
        fig,ax=init_figure();
        plot2D(ftable_imag,0.1,maxreal,-maximag,maximag,'Re','Im','',ax,colorlabel='',colorlim=[-1,1]);
        #plot(gridR,gridI,'','.k','Im',ax,0,xlab='Re',ms=2);
        plot(np.array(list1)[:,0],np.array(list1)[:,1],'potential zeros','ow','Im',ax,0,xlab='Re',ms=5);
        plot(np.real(list2),np.imag(list2),'zeros','xk','Im',ax,0,xlab='Re',ms=15);
        ax.set_xlim([0.1,maxreal]);ax.set_ylim([-maximag,maximag]);
        end_figure(fig,ax);
        fig,ax=init_figure();
        plot2D(invftable,0.1,maxreal,-maximag,maximag,'Re','Im','',ax,colorlabel='');
        #plot(gridR,gridI,'','.k','Im',ax,0,xlab='Re',ms=2);
        plot(np.array(list1)[:,0],np.array(list1)[:,1],'potential zeros','ow','Im',ax,0,xlab='Re',ms=5);
        plot(np.real(list2),np.imag(list2),'zeros','xk','Im',ax,0,xlab='Re',ms=15);
        ax.set_xlim([0.1,maxreal]);ax.set_ylim([-maximag,maximag]);
        end_figure(fig,ax);

        if False:
            #3D plots
            from mpl_toolkits.mplot3d import Axes3D
            fig = pylab.figure();ax3D = Axes3D(fig);
            ax3D.plot(gridR, gridI, 0.);
            ax3D.plot_wireframe(gridR, gridI, ftable_real);
            fig = pylab.figure();ax3D = Axes3D(fig);
            ax3D.plot(gridR, gridI, 0.);
            ax3D.plot_wireframe(gridR, gridI, ftable_imag);
            fig = pylab.figure();ax3D = Axes3D(fig);
            ax3D.plot(gridR, gridI, 0.);
            ax3D.plot_wireframe(gridR, gridI, np.abs(ftable_real**2+ftable_imag**2));


        pylab.show()

        # results: with tauAC=0
        # - m=0 (long.): f[THz]=
        # - m=1 (trans.): f[THz]=
        # results: with tauAC=4.2 ps
        # - m=0 (long.): f[THz]= 0.74527449 -7.82611265e-05j
        # - m=1 (trans.): f[THz]=
