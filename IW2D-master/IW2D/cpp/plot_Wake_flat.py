#!/usr/bin/python

import sys
import commands
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);
import pylab
import numpy as np
import math,cmath
from matplotlib.font_manager import FontProperties
import matplotlib.text as text
from string import split, replace
from optparse import OptionParser


def parsse():
    parser = OptionParser()
    parser.add_option("-c", "--constant",action="store_true",
                      help="Specify if we should plot the constant term Wycst",
                      metavar="CST", default=False,dest="CST")
    parser.add_option("-f", "--file",action="append",
                      help="Specify the file name (beginning with Wlong - the other wake files Wxdip, Wydip, etc. are selected automatically)",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-g", "--legend",action="append",
                      help="Specify the legend for each file",
                      metavar="LEG", default=None,dest="LEG")
    parser.add_option("-l", "--loglog",action="store_true",
                      help="Specify if loglog plot (semilogy by default)",
                      metavar="LOG", default=False,dest="LOG")
    parser.add_option("-m", "--meter",action="store_true",
                      help="Specify if the wakes are per meter length",
                      metavar="METER", default=False,dest="METER")
    parser.add_option("-o", "--outfile",
                      help="Specify the output file name if you want to output the plots in files (will create Wlong[this_option].png, Wxdip[this_option].png, etc. and the same in .eps format) (by default: print plots on screen -> CANNOT WORK IF YOU USE MATPLOTLIB WITH GUI DISABLED)",
                      metavar="OUT", default=None,dest="OUT")
    (opt, args) = parser.parse_args()
    #print "Selected Files:", opt.FILE
    return opt, args


def set_fontsize(fig,size):

    # set all fontsize to 'size' in the figure 'fig'
    for o in fig.findobj(text.Text): o.set_fontsize(size);


def init_figure():

    fig=pylab.figure(facecolor='w', edgecolor=None)
    pylab.axes([0.12,0.12,0.85,0.85])
    ax=fig.gca()
    
    return fig,ax
    

def logplot(z,W,lab,leg1,leg2,col):
    # plot wake and put legend and label
    pylab.loglog(z,np.abs(W),col+'-',label="|"+leg1+"|"+leg2,lw=2.5,ms=10.,mew=2.5)
    pylab.xlabel("Distance behind the source [m]")
    pylab.ylabel("W "+lab);


def semilogplot(z,W,lab,leg,col):
    # plot wake with labels, legend and color
    pylab.semilogy(z,W,col+'-',label=leg,lw=2.5,ms=10.,mew=2.5)
    pylab.xlabel("Distance behind the source [m]")
    pylab.ylabel("W "+lab);


def read(filename):
    # read wake file (2 columns: z, wake)
    z=[];W=[];
    fh=open(filename,"r")
    for l in fh.readlines():
    	if not l.startswith('Dist'):
		ll=l.strip().split();
		z.append(ll[0]);
		W.append(ll[1]);
    fh.close()
    z=np.array([float(j) for j in z])
    W=np.array([float(j) for j in W])
    return z,W
    
                  
if __name__ == "__main__":
    opt,args=parsse();
    
    if (opt.CST): nfig=6;
    else: nfig=5;
 
    for i in range(nfig): init_figure();
    
    color=['b','r','g','m','c','y','k'];
    comp=["Wlong","Wxdip","Wydip","Wxquad","Wyquad","Wycst"];

    for j,fil in enumerate(opt.FILE):

	# string for legend
	if opt.LEG==None: leg=fil.replace("WlongW","",1).replace("_"," ");
	else: leg=opt.LEG[j];
	
	# constructs all the files names (flat chamber impedances)
	filefxdip=fil.replace("Wlong","Wxdip",1);
	filefxquad=fil.replace("Wlong","Wxquad",1);
	filefydip=fil.replace("Wlong","Wydip",1);
	filefyquad=fil.replace("Wlong","Wyquad",1);
	filefycst=fil.replace("Wlong","Wycst",1);

	# read longitudinal impedance and plot
	z,W=read(fil)
	if not opt.METER: lab=" [V / C]"
	else: lab=" [V / (C.m)]"
	pylab.figure(1);
	if opt.LOG: logplot(z,W,lab,"$W_{\|\|}$",", "+leg,color[j]);
	else: semilogplot(z,W,lab,"$W_{\|\|}$, "+leg,color[j]);

	if opt.CST:
    	    # read constant vertical impedance and plot
    	    pylab.figure(6);
	    z,W=read(filefycst)
    	    if opt.LOG: logplot(z,W,lab,"$W_y^{cst}$",", "+leg,color[j]);
    	    else: semilogplot(z,W,lab,"$W_y^{cst}$, "+leg,color[j]);


	if not opt.METER: lab=" [V / (C.m)]"
	else: lab=" [V / (C.m$^2$)]"
	# read transverse dipolar impedances and plot them    
	z1,W1=read(filefxdip)
	z2,W2=read(filefydip)
	pylab.figure(2);
	if opt.LOG:
    	    logplot(z1,W1,lab,"$W_x^{dip}$",", "+leg,color[j]);
    	    pylab.figure(3);logplot(z2,W2,lab,"$W_y^{dip}$",", "+leg,color[j]);
	else:
    	    semilogplot(z1,W1,lab,"$W_x^{dip}$, "+leg,color[j]);
    	    pylab.figure(3);semilogplot(z2,W2,lab,"$W_y^{dip}$, "+leg,color[j]);

	# read transverse quadrupolar impedances and plot them
	z1,W1=read(filefxquad)
	z2,W2=read(filefyquad)
	pylab.figure(4);
	if opt.LOG:
    	    logplot(z1,W1,lab,"$W_x^{quad}$",", "+leg,color[j]);
    	    pylab.figure(5);logplot(z2,W2,lab,"$W_y^{quad}$",", "+leg,color[j]);
	else:
    	    semilogplot(z1,W1,lab,"$W_x^{quad}$, "+leg,color[j]);
    	    pylab.figure(5);semilogplot(z2,W2,lab,"$W_y^{quad}$, "+leg,color[j]);
    
    for i in range(nfig):
    	fig=pylab.figure(i+1);ax=fig.gca();
	set_fontsize(pylab.figure(i+1),'x-large');
    	ax.legend(loc=0);
    	l=ax.get_legend_handles_labels();
	set_fontsize(ax.get_legend(),max(18-2*len(l[0]),10));
	ax.grid();
    
	if opt.OUT!=None:
	    fig.savefig(comp[i]+opt.OUT+".png")
	    fig.savefig(comp[i]+opt.OUT+".eps",format='eps')
	    pylab.close(fig);

	    	
    if opt.OUT==None: pylab.show();
    
    sys.exit()
