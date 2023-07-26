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
from matplotlib.patches import Rectangle
from optparse import OptionParser
from plot_lib import *


def parsse():
    parser = OptionParser()
    parser.add_option("-f", "--file",action="append",
                      help="Specify the curve file name (impedance or wake) (several -f options possible)",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-g", "--legend",action="append",
                      help="Specify the legend for each curve",
                      metavar="LEG", default=None,dest="LEG")
    parser.add_option("-r", "--anchor",type=float,nargs=2,action="append",
                      help="Specify the place to anchor the upper-right corner of the legend (2 -r options possible: real and imaginary parts))",
                      metavar="ANC", default=[],dest="ANC")
    parser.add_option("-t", "--totalfile",
                      help="Specify the file name for the full model (impedance or wake)",
                      metavar="TFILE", default=None,dest="TFILE")
    (opt, args) = parser.parse_args()
    #print "Selected Files:", opt.FILE
    return opt, args


def read(filename):
    # read impedance or wake file
    # 3 columns for imp.: frequency, real and imag. parts
    # 2 columns for wake: z, wake
    # flagimp indicates if it's an impedance file (otherwise wake file)
    freq=[];Z1=[];Z2=[];z=[];W=[];
    flagimp=0;
    fh=open(filename,"r")
    for l in fh.readlines():
        ll=l.strip().split();
        if (len(ll)==3):
            flagimp=1;
            if (len(freq)==0)or(ll[0]!=freq[-1]): # to avoid duplicate frequencies
                freq.append(ll[0]);
                Z1.append(ll[1]);
                Z2.append(ll[2]);
        elif (len(ll)==2):
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


if __name__ == "__main__":
    opt,args=parsse();

    if len(opt.ANC)==0: opt.ANC.append((0.5,1.1));

    #color=['b','r','g','m','c','y','k'];
    color=[(0.3,0.3,1,0),(0.7,0.3,0.3,0),(0.3,0.7,0.3,0),(0.7,0.3,0.8,0),(0,1,1,0),(1,1,0,0),(0,0,0,0)];

    # read the file with the total model
    xt,yt,flagimpt=read(opt.TFILE);

    if (flagimpt==1):
        figr,axr=init_figure();axr=pylab.subplot(111);
        figi,axi=init_figure();axi=pylab.subplot(111);
        pr=[];pi=[];
        if len(opt.ANC)==1: opt.ANC.append(opt.ANC[0])
    else:
        fig,ax=init_figure();ax=pylab.subplot(111);
        p=[];

    xsum=xt;
    ysum=np.zeros(xt.size);
    leg=[];

    for j,fil in enumerate(opt.FILE):

        # string for legend
        if opt.LEG==None: leg.append(fil.replace("_"," "));
        else: leg.append(opt.LEG[j]);

        # read the file
        x,y,flagimp=read(fil);

        if (flagimp != flagimpt): print("Mixing impedance and wake !");sys.exit();

        if (flagimp==1):
            # interpolate total model at the absissae x, and calculate the ratio
            ytrinterp=np.interp(x,xt,yt.real);
            ytiinterp=np.interp(x,xt,yt.imag);
            #ratior=np.abs(y.real/ytrinterp);
            #ratioi=np.abs(y.imag/ytiinterp);
            ratior=y.real/ytrinterp;# uncomment (and comment the 2 previous lines) when not taking absolute values
            ratioi=y.imag/ytiinterp;
            # calculate the ratio for the previous sum
            #sumratior=np.interp(x,xsum,np.abs(ysum.real/yt.real));
            #sumratioi=np.interp(x,xsum,np.abs(ysum.imag/yt.imag));
            sumratior=np.interp(x,xsum,ysum.real/yt.real);# uncomment (and comment the 2 previous lines) when not taking absolute values
            sumratioi=np.interp(x,xsum,ysum.imag/yt.imag);
            lab="Frequency [Hz]";
            # plot ratio with filling between curves (real and imag part)
            fillplot_percent(x,sumratior+ratior,sumratior,lab,leg,color[j],axr)
            fillplot_percent(x,sumratioi+ratioi,sumratioi,lab,leg,color[j],axi)
            # with fill_between we have to create the legend box ourselves
            pr.append(Rectangle((0, 0), 1, 1,axes=axr,color=color[j]));
            pi.append(Rectangle((0, 0), 1, 1,axes=axi,color=color[j]));
            #update the sum
            ysum=ysum+np.interp(xsum,x,y.real)+1j*np.interp(xsum,x,y.imag);
        else:
            # interpolate total model at the absissae x, and calculate the ratio
            ytinterp=np.interp(x,xt,yt);
            #ratio=np.abs(y/ytinterp);
            ratio=y/ytinterp;# uncomment (and comment the previous line) when not taking absolute values
            # calculate the ratio for the previous sum
            #sumratio=np.interp(x,xsum,np.abs(ysum/yt));
            sumratio=np.interp(x,xsum,ysum/yt);# uncomment (and comment the previous line) when not taking absolute values
            lab="Distance behind the source [m]";
            # plot ratio with filling between curves (real and imag part)
            fillplot_percent(x,sumratio+ratio,sumratio,lab,leg,color[j],ax)
            # with fill_between we have to create the legend box ourselves
            p.append(Rectangle((0, 0), 1, 1,axes=ax,color=color[j]));
            #update the sum
            ysum=ysum+np.interp(xsum,x,y);




    fontP = FontProperties();fontP.set_size('xx-large');
    if (flagimp==1):
        axr.set_xscale('log');axi.set_xscale('log');
        axr.set_ylim([0,110]);axi.set_ylim([0,110]);
        #axr.set_ylim([-100,300]);axi.set_ylim([-100,300]); # uncomment (and comment the previous line) when not taking absolute values
        set_fontsize(figr,'xx-large');set_fontsize(figi,'xx-large');
        axr.legend(pr,leg,bbox_to_anchor=opt.ANC[0],prop = fontP);
        axi.legend(pi,leg,bbox_to_anchor=opt.ANC[1],prop = fontP);
    else:
        ax.set_xscale('log');
        ax.set_ylim([0,110]);
        #ax.set_ylim([-100,300]);# uncomment (and comment the previous line) when not taking absolute values
        set_fontsize(fig,'xx-large');
        ax.legend(p,leg,bbox_to_anchor=opt.ANC[0],prop = fontP);



    pylab.show();
    sys.exit()
