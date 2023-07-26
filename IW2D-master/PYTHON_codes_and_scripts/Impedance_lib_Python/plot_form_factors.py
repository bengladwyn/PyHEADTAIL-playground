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
from matplotlib.font_manager import FontProperties
from optparse import OptionParser
from plot_lib import init_figure,set_fontsize


def parsse():
    parser = OptionParser()
    parser.add_option("-a", "--all",action="store_true",
                      help="Specify if we plot all the form factors on the same plot",
                      metavar="ALL", default=False,dest="ALL")
    parser.add_option("-b", "--burovdanilov",action="store_true",
                      help="Specify if we plot the Burov-Danilov form factors for a single-plate flat chamber",
                      metavar="BUROV", default=False,dest="BUROV")
    parser.add_option("-c", "--circfile",action="append",
                      help="Specify the circular chamber file name (longitudinal impedance or wake)",
                      metavar="CFILE", default=None,dest="CFILE")
    parser.add_option("-f", "--flatfile",action="append",
                      help="Specify the flat chamber file name (longitudinal impedance or wake)",
                      metavar="FFILE", default=None,dest="FFILE")
    parser.add_option("-g", "--legend",action="append",
                      help="Specify the legend for each form factor",
                      metavar="LEG", default=None,dest="LEG")
    parser.add_option("-l", "--longi",action="store_true",
                      help="Specify if we plot also the longitudinal form factor",
                      metavar="LONG", default=False,dest="LONG")
    parser.add_option("-m", "--ylim",type=float,nargs=2,
                      help="Specify the limits on the y axis (with -a option)",
                      metavar="YLIM", default=[-1,3],dest="YLIM")
    parser.add_option("-o", "--loglog",action="store_true",
                      help="Specify if we do a loglog plot (instead of semilogx)",
                      metavar="LOG", default=False,dest="LOG")
    parser.add_option("-q", "--anchor_quad",type=float,nargs=2,
                      help="Specify the place to anchor the upper-right corner of the legend for quadrupolar form factors (only with -a option)",
                      metavar="ANC", default=(0.95,1.1),dest="ANCQ")
    parser.add_option("-r", "--anchor",type=float,nargs=2,
                      help="Specify the place to anchor the upper-right corner of the legend",
                      metavar="ANC", default=(0.45,1.1),dest="ANC")
    parser.add_option("-x", "--xlim",type=float,nargs=2,
                      help="Specify the limits on the x axis",
                      metavar="XLIM", default=None,dest="XLIM")
    parser.add_option("-y", "--yokoya",action="store_true",
                      help="Specify if we plot the Yokoya form factors for a flat chamber",
                      metavar="YOK", default=False,dest="YOK")
    (opt, args) = parser.parse_args()
    #print "Selected Files:", opt.FILE
    return opt, args

def semilogplot_comp(x,y,lab,leg,col,flaglog=False,lw=2.5,ms=12,ax=None):
    # y complex: plot real and imag parts, with labels, legend and color
    legr=leg+", real part"
    legi=leg+", imag. part"
    if ax is None: ax=pylab.gca()
    if flaglog:
        ax.loglog(x,np.abs(y.real),col,label=legr,lw=lw,ms=ms,mew=2.5)
        ax.loglog(x,np.abs(y.imag),col+'-',label=legi,lw=lw,ms=ms,mew=2.5)
    else:
        ax.semilogx(x,y.real,col,label=legr,lw=lw,ms=ms,mew=2.5)
        ax.semilogx(x,y.imag,col+'-',label=legi,lw=lw,ms=ms,mew=2.5)
    ax.set_xlabel(lab)
    ax.set_ylabel("Form factor");


def semilogplot(x,y,lab,leg,col,flaglog=False,lw=2.5,ms=12,ax=None):
    # y real: plot  with labels, legend and color
    if ax is None: ax=pylab.gca()
    if flaglog: ax.loglog(x,np.abs(y.real),col,label=leg,lw=lw,ms=ms,mew=2.5,mfc='w',mec=col[1])
    else: ax.semilogx(x,y.real,col,label=leg,lw=lw,ms=ms,mew=2.5,mfc='w',mec=col[1])
    ax.set_xlabel(lab)
    ax.set_ylabel("Form factor");


def read(filename):
    # read impedance or wake file
    # 3 columns for imp.: frequency, real and imag. parts
    # 2 columns for wake: z, wake
    # flagimp indicates if it's an impedance file (otherwise wake file)
    freq=[];Z1=[];Z2=[];z=[];W=[];
    flagimp=0;
    fh=open(filename,"r")
    for l in fh.readlines():
        if l.startswith('Freq'): flagimp=1;
        elif (not l.startswith('Dist'))and(flagimp==1):
            ll=l.strip().split();
            if (len(freq)==0)or(ll[0]!=freq[-1]): # to avoid duplicate frequencies
                freq.append(ll[0]);
                Z1.append(ll[1]);
                Z2.append(ll[2]);
        elif (not l.startswith('Dist')):
            ll=l.strip().split();
            if (len(z)==0)or(ll[0]!=z[-1]): # to avoid duplicate frequencies
                z.append(ll[0]);
                W.append(ll[1]);
    fh.close()

    if (flagimp==1):
        x=np.array([float(j) for j in freq])
        y=np.array([float(j) for j in Z1])+1j*np.array([float(j) for j in Z2])
    else:
        x=np.array([float(j) for j in z])
        y=np.array([float(j) for j in W])

    return x,y,flagimp


def form(fil1,fil2,i,leg,col,flaglog=False,lw=2.5, ax=None):
    # compute the form factor between two files
    # i is the figure number, leg the legend, col the color
    x1,y1,flagimp1=read(fil1);
    x2,y2,flagimp2=read(fil2);
    #print x1.size,y1.size,x2.size,y2.size
    if (flagimp1 != flagimp2): print("Mixing impedance and wake !");sys.exit();

    pylab.figure(i);
    if (flagimp1==1):
        y2rinterp=np.interp(x1,x2,y2.real);
        y2iinterp=np.interp(x1,x2,y2.imag);
        fact=(y2rinterp+1j*y2iinterp)/y1;
        lab="Frequency [Hz]"
        semilogplot_comp(x1,fact,lab,leg,col+'-',flaglog=flaglog,lw=lw,ax=ax);
    else:
        y2interp=np.interp(x1,x2,y2);
        fact=y2interp/y1;
        lab="Distance behind the source [m]";
        semilogplot(x1,fact,lab,leg,'-'+col,flaglog=flaglog,lw=lw,ax=ax);

    return x1,fact,lab;


if __name__ == "__main__":
    opt,args=parsse();


    if (opt.ALL):
        fig,ax=init_figure();ax=pylab.subplot(111);ax2=pylab.twinx(ax=ax)
    else:
        for i in range(5): init_figure();

    color=['b','r','g','m','c','y','k'];

    for j,fileflong in enumerate(opt.FFILE):

        # string for legend
        if opt.LEG==None: leg=fileflong.replace("ZlongW","",1).replace("_"," ").replace("WlongW","",1);
        else: leg=opt.LEG[j];

        # constructs all the files names (flat chamber)
        filefxdip=fileflong.replace("long","xdip",1);
        filefxquad=fileflong.replace("long","xquad",1);
        filefydip=fileflong.replace("long","ydip",1);
        filefyquad=fileflong.replace("long","yquad",1);

        # constructs all the files names (circular chamber)
        fileclong=opt.CFILE[j];
        filecxdip=fileclong.replace("long","xdip",1);

        i=1;kwargs={}
        npts=20;
        if opt.ALL: kwargs={'ax': ax}
        x,y,flagimp=read(fileclong);
        if (flagimp==1): xcst=10.**np.linspace(-3,15,npts); # absissae for constant form factors (impedances)
        else: xcst=10.**np.linspace(-6,6,12); # absissae for constant form factors (wakes)
        # longitudinal form factor
        if opt.LONG:
            x,y,lab=form(fileclong,fileflong,i,"Longitudinal "+leg,color[j],flaglog=opt.LOG,**kwargs)
            if (opt.YOK): semilogplot(xcst,np.ones(xcst.size),lab,"Longitudinal Yokoya factor=1",'s'+color[j],flaglog=opt.LOG,**kwargs);
            pylab.ylim(-0.1,1.5)

        # x dipolar form factor
        if not opt.ALL: i=2;
        else: j=j+1;kwargs={'ax': ax}
        x,y,lab=form(filecxdip,filefxdip,i,"Dipolar (x) "+leg,color[j],flaglog=opt.LOG,lw=8.,**kwargs)
        if (opt.YOK): semilogplot(xcst,np.ones(xcst.size)*np.pi*np.pi/24.,lab,"Dipolar (x) Yokoya factor=$\pi^2/24$",'o'+color[j],flaglog=opt.LOG,**kwargs);
        elif (opt.BUROV): semilogplot(xcst,np.ones(xcst.size)/4.,lab,"Dipolar (x) Burov-Danilov factor=$1/4$",'o'+color[j],flaglog=opt.LOG,**kwargs);
        pylab.ylim(-0.1,0.8)

        # y dipolar form factor
        if not opt.ALL: i=3;
        else: j=j+1;kwargs={'ax': ax}
        x,y,lab=form(filecxdip,filefydip,i,"Dipolar (y) "+leg,color[j],flaglog=opt.LOG,**kwargs)
        if (opt.YOK): semilogplot(xcst,np.ones(xcst.size)*np.pi*np.pi/12.,lab,"Dipolar (y) Yokoya factor=$\pi^2/12$",'D'+color[j],flaglog=opt.LOG,**kwargs);
        elif (opt.BUROV): semilogplot(xcst*np.sqrt(10.),np.ones(xcst.size)/4.,lab,"Dipolar (y) Burov-Danilov factor=$1/4$",'D'+color[j],flaglog=opt.LOG,**kwargs);
        pylab.ylim(-0.1,1.2)

        # x quadrupolar form factor
        if not opt.ALL:
            i=4;
        else:
            j=j+1;kwargs={'ax': ax2}
        if not opt.LOG:
            x,y,lab=form(filecxdip,filefxquad,i,"Quadrupolar (x) "+leg,color[j],flaglog=opt.LOG,**kwargs)
            if (opt.YOK): semilogplot(xcst,-np.ones(xcst.size)*np.pi*np.pi/24.,lab,"Quadrupolar (x) Yokoya factor=$-\pi^2/24$",'+'+color[j],flaglog=opt.LOG,**kwargs);
            elif (opt.BUROV): semilogplot(xcst,-np.ones(xcst.size)/4.,lab,"Quadrupolar (x) Burov-Danilov factor=$-1/4$",'+'+color[j],flaglog=opt.LOG,**kwargs);
            pylab.ylim(-0.8,0.1)

        # y quadrupolar form factor
        if not opt.ALL:
            i=5;
        else:
            j=j+1;kwargs={'ax': ax2}
        if not opt.BUROV:
            x,y,lab=form(filecxdip,filefyquad,i,"Quadrupolar (y) "+leg,color[j],flaglog=opt.LOG,**kwargs)
            if (opt.YOK): semilogplot(xcst*np.sqrt(10.),np.ones(xcst.size)*np.pi*np.pi/24.,lab,"Quadrupolar (y) Yokoya factor=$\pi^2/24$",'x'+color[j],flaglog=opt.LOG,**kwargs);
            elif (opt.BUROV): semilogplot(xcst*4.,np.ones(xcst.size)/4.,lab,"Quadrupolar (y) Burov-Danilov factor=$1/4$",'x'+color[j],flaglog=opt.LOG,**kwargs);
            pylab.ylim(-0.1,0.8)

        if (opt.ALL):
            ax.set_ylim(opt.YLIM[0],opt.YLIM[1]);
            if opt.XLIM is not None: ax.set_xlim(opt.XLIM[0],opt.XLIM[1])
            ax2.set_ylim(opt.YLIM[0],opt.YLIM[1]);
            set_fontsize(pylab.figure(1),'xx-large');
            if opt.LONG:
                fontP = FontProperties();fontP.set_size('x-large');
                ax.legend(bbox_to_anchor=opt.ANC,prop = fontP);
                ax2.legend(bbox_to_anchor=opt.ANCQ,prop = fontP);
            else:
                fontP = FontProperties();fontP.set_size('xx-large');
                ax.legend(bbox_to_anchor=opt.ANC,prop = fontP);
                ax2.legend(bbox_to_anchor=opt.ANCQ,prop = fontP);
        else:
            for i in range(5):
                pylab.figure(i+1);
                if opt.XLIM is not None: pylab.xlim(opt.XLIM[0],opt.XLIM[1])
                pylab.legend(loc=2);
                set_fontsize(pylab.figure(i+1),'xx-large');

    pylab.show();
    sys.exit()
