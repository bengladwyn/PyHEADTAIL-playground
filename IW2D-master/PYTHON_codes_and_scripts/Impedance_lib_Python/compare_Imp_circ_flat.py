#!/usr/bin/python

import sys
import subprocess
pymod=subprocess.check_output("echo $PYMOD",shell=True).strip().decode()
if pymod.startswith('local'):
    py_numpy=subprocess.check_output("echo $PY_NUMPY",shell=True).strip().decode()
    sys.path.insert(1,py_numpy)
    py_matpl=subprocess.check_output("echo $PY_MATPL",shell=True).strip().decode()
    sys.path.insert(1,py_matpl)
import pylab
import numpy as np
import math,cmath
from optparse import OptionParser


def parsse():
    parser = OptionParser()
    parser.add_option("-c", "--filec",
                      help="Specify the file name (longitudinal impedance, circular case)",
                      metavar="FILEC", default=None,dest="FILEC")
    parser.add_option("-f", "--filef",
                      help="Specify the file name (longitudinal impedance, flat case)",
                      metavar="FILEF", default=None,dest="FILEF")
    parser.add_option("-l", "--loglog",action="store_true",
                      help="Specify if loglog plot (semilogx by default)",
                      metavar="LOG", default=False,dest="LOG")
    parser.add_option("-m", "--meter",action="store_true",
                      help="Specify if the impedance are per meter length",
                      metavar="METER", default=False,dest="METER")
    parser.add_option("-v", "--vertical",action="store_true",
                      help="Specify if vertical structure (Zydip twice as large as Zxdip)",
                      metavar="VERT", default=False,dest="VERT")
    (opt, args) = parser.parse_args()
    print("Selected Files:", opt.FILEF, ", ", opt.FILEC)
    return opt, args

def logplot(freq1,freq2,Z1,Z2,lab,leg1,leg2):
    # plot two impedances (real and imag parts) and put legend and label
    leg1r="Re["+leg1+"]"
    leg1i="Im["+leg1+"]"
    leg2r="Re["+leg2+"]"
    leg2i="Im["+leg2+"]"
    Z1r=np.array([math.fabs(j) for j in Z1.real]);
    Z1i=np.array([math.fabs(j) for j in Z1.imag]);
    Z2r=np.array([math.fabs(j) for j in Z2.real]);
    Z2i=np.array([math.fabs(j) for j in Z2.imag]);
    pylab.figure()
    pylab.loglog(freq1,Z1r,'b-',label="|"+leg1r+"|")
    pylab.loglog(freq1,Z1i,'b--',label="|"+leg1i+"|")
    pylab.loglog(freq2,Z2r,'r-',label="|"+leg2r+"|")
    pylab.loglog(freq2,Z2i,'r--',label="|"+leg2i+"|")
    pylab.xlabel("Frequency [Hz]")
    pylab.ylabel("Z "+lab);
    pylab.legend(loc=0);

def semilogplot(freq1,freq2,Z1,Z2,lab,leg1,leg2):
    # plot two impedances (real and imag parts) and put legend and label
    leg1r="Re["+leg1+"]"
    leg1i="Im["+leg1+"]"
    leg2r="Re["+leg2+"]"
    leg2i="Im["+leg2+"]"
    pylab.figure()
    pylab.semilogx(freq1,Z1.real,'b-',label=leg1r)
    pylab.semilogx(freq1,Z1.imag,'b--',label=leg1i)
    pylab.semilogx(freq2,Z2.real,'r-',label=leg2r)
    pylab.semilogx(freq2,Z2.imag,'r--',label=leg2i)
    pylab.xlabel("Frequency [Hz]")
    pylab.ylabel("Z "+lab);
    pylab.legend(loc=0);


def read(filename):
    # read impedance file (3 columns: frequency, real and imag. parts)
    freq=[];Z1=[];Z2=[]
    fh=open(filename,"r")
    for l in fh.readlines():
        if not l.startswith('Freq'):
            ll=l.strip().split();
            freq.append(ll[0]);
            Z1.append(ll[1]);
            Z2.append(ll[2]);
    fh.close()
    f=np.array([float(j) for j in freq])
    Z=np.array([float(j) for j in Z1])+1j*np.array([float(j) for j in Z2])
    return f,Z


if __name__ == "__main__":
    opt,args=parsse();

    # constructs all the files names (circular, then flat impedances)
    filecdip=opt.FILEC.replace("Zlong","Ztransdip",1);
    filecquad=opt.FILEC.replace("Zlong","Ztransquad",1);
    filefxdip=opt.FILEF.replace("Zlong","Zxdip",1);
    filefxquad=opt.FILEF.replace("Zlong","Zxquad",1);
    filefydip=opt.FILEF.replace("Zlong","Zydip",1);
    filefyquad=opt.FILEF.replace("Zlong","Zyquad",1);

    # read longitudinal impedances and plot them
    freq1,Z1=read(opt.FILEC)
    freq2,Z2=read(opt.FILEF)
    if not opt.METER: lab=" [$\Omega$]"
    else: lab=" [$\Omega$ /m]"
    if opt.LOG: logplot(freq1,freq2,Z1,Z2,lab,"$Z_{\|\|}$ (circular)","$Z_{\|\|}$ (flat)");
    else: semilogplot(freq1,freq2,Z1,Z2,lab,"$Z_{\|\|}$ (circular)","$Z_{\|\|}$ (flat)");

    if not opt.METER: lab=" [$\Omega$ /m]"
    else: lab=" [$\Omega$ /$m^2$]"
    # Yokoya factors
    if opt.VERT:
        Yokxdip=math.pi*math.pi/24.;
        Yokydip=math.pi*math.pi/12.;
        Yokxquad=-math.pi*math.pi/24.;
    else:
        Yokxdip=math.pi*math.pi/12.;
        Yokydip=math.pi*math.pi/24.;
        Yokxquad=math.pi*math.pi/24.;
    Yokyquad=-Yokxquad;

    # read transverse dipolar impedances and plot them
    freq1,Z1=read(filecdip)
    freq2,Z2=read(filefxdip)
    freq3,Z3=read(filefydip)
    if opt.LOG:
        logplot(freq1,freq2,Z1*Yokxdip,Z2,lab,"$Z^{dip}$ (circular) $*$ Yok. factor %.2f" % Yokxdip,"$Z_x^{dip}$ (flat)");
        logplot(freq1,freq3,Z1*Yokydip,Z3,lab,"$Z^{dip}$ (circular) $*$ Yok. factor %.2f" % Yokydip,"$Z_y^{dip}$ (flat)");
    else:
        semilogplot(freq1,freq2,Z1*Yokxdip,Z2,lab,"$Z^{dip}$ (circular) $*$ Yok. factor %.2f" % Yokxdip,"$Z_x^{dip}$ (flat)");
        semilogplot(freq1,freq3,Z1*Yokydip,Z3,lab,"$Z^{dip}$ (circular) $*$ Yok. factor %.2f" % Yokydip,"$Z_y^{dip}$ (flat)");

    # read transverse quadrupolar impedances and plot them
    #freq1,Z1=read(filecquad)
    freq2,Z2=read(filefxquad)
    freq3,Z3=read(filefyquad)
    if opt.LOG:
        logplot(freq1,freq2,Z1*Yokxquad,Z2,lab,"$Z^{dip}$ (circular) $*$ Yok. factor %.2f" % Yokxquad,"$Z_x^{quad}$ (flat)");
        logplot(freq1,freq3,Z1*Yokyquad,Z3,lab,"$Z^{dip}$ (circular) $*$ Yok. factor %.2f" % Yokyquad,"$Z_y^{quad}$ (flat)");
    else:
        semilogplot(freq1,freq2,Z1*Yokxquad,Z2,lab,"$Z^{dip}$ (circular) $*$ Yok. factor %.2f" % Yokxquad,"$Z_x^{quad}$ (flat)");
        semilogplot(freq1,freq3,Z1*Yokyquad,Z3,lab,"$Z^{dip}$ (circular) $*$ Yok. factor %.2f" % Yokyquad,"$Z_y^{quad}$ (flat)");

    pylab.show();sys.exit()
