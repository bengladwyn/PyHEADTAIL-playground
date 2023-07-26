#!/usr/bin/python

import sys
import commands
pymod=commands.getoutput("echo $PYMOD")
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy)
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl)
import pylab
import numpy as np
import math,cmath
from matplotlib.font_manager import FontProperties
import matplotlib.text as text
from string import split, replace
from optparse import OptionParser
from glob import glob
from copy import deepcopy


def parsse():
    parser = OptionParser()
    parser.add_option("-f", "--file",action="append",
                      help="Specify the file name (beginning with Zlong - the other impedance files Zxdip, Zydip, Zl1000, etc. are selected automatically)",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-g", "--legend",action="append",
                      help="Specify the legend for each file",
                      metavar="LEG", default=None,dest="LEG")
    parser.add_option("-l", "--loglog",action="store_true",
                      help="Specify if loglog plot (semilogx by default)",
                      metavar="LOG", default=False,dest="LOG")
    parser.add_option("-x", "--beamsizex",default=None,dest="SIZEX",
                      help="Specify a (geometric) horizontal beam size, in m (used to put all impedance terms on the same ground)",
                      metavar="SIZEX",type=float)
    parser.add_option("-y", "--beamsizey",default=None,dest="SIZEY",
                      help="Specify a (geometric) vertical beam size, in m (used to put all impedance terms on the same ground)",
                      metavar="SIZEY",type=float)
    parser.add_option("-o", "--outfile",
                      help="Specify the output file name if you want to output the plots in files (will create Zall[this_option].png and the same in .eps format) (by default: print plots on screen -> CANNOT WORK IF YOU USE MATPLOTLIB WITH GUI DISABLED)",
                      metavar="OUT", default=None,dest="OUT")
    (opt, args) = parser.parse_args()
    #print "Selected Files:", opt.FILE
    return opt, args


def set_fontsize(fig,size):

    # set all fontsize to 'size' in the figure 'fig'
    for o in fig.findobj(text.Text):
        o.set_fontsize(size)


def init_figure():

    fig=pylab.figure(facecolor='w', edgecolor=None)
    pylab.axes([0.15,0.15,0.8,0.8])
    ax=fig.gca()
    
    return fig,ax
    

def logplot(freq,Z,unit,leg1,leg2,col):
    # plot two impedances (real and imag parts) and put legend and label
    legr="|Re["+leg1+"]|"
    legi="|Im["+leg1+"]|"
    Zr=np.array([math.fabs(j) for j in Z.real])
    Zi=np.array([math.fabs(j) for j in Z.imag])
    pylab.loglog(freq,Zr,col+'-',label=legr+leg2,lw=2.5,ms=10.,mew=2.5)
    pylab.loglog(freq,Zi,col+'--',label=legi+leg2,lw=2.5,ms=10.,mew=2.5)
    pylab.xlabel("Frequency [Hz]")
    pylab.ylabel("Z "+unit)


def semilogplot(freq,Z,unit,leg1,leg2,col):
    # plot impedance (real and imag parts) with labels, legend and color
    legr="Re["+leg1+"]"
    legi="Im["+leg1+"]"
    pylab.semilogx(freq,Z.real,col+'-',label=legr+leg2,lw=2.5,ms=10.,mew=2.5)
    pylab.semilogx(freq,Z.imag,col+'--',label=legi+leg2,lw=2.5,ms=10.,mew=2.5)
    pylab.xlabel("Frequency [Hz]")
    pylab.ylabel("Z "+unit)


def read(filename):
    # read impedance file (3 columns: frequency, real and imag. parts)
    freq=[];Z1=[];Z2=[]
    fh=open(filename,"r")
    for l in fh.readlines():
        if not l.startswith('Freq'):
            ll=l.strip().split()
            freq.append(ll[0])
            Z1.append(ll[1])
            Z2.append(ll[2])
    fh.close()
    f=np.array([float(j) for j in freq])
    Z=np.array([float(j) for j in Z1])+1j*np.array([float(j) for j in Z2])
    return f,Z
    
                  
if __name__ == "__main__":
    opt,args=parsse()
    
    # one figure per plane
    for i in range(3):
        init_figure()
    
    color = ['b','r','g','m','c','y','k']
    colors = deepcopy(color)
    patterns = ['','x','.','+','d','^','v','<','>']
    for pattern in patterns:
        colors += [pattern+c for c in color]
    
    ind_color = 0
    for j,fil in enumerate(opt.FILE):

        # string for legend
        if opt.LEG==None: leg=fil.replace("ZlongW","",1).replace("_"," ")
        else: leg=opt.LEG[j]
        
        for iplane,(plane,Zlabel) in enumerate(zip(['l','x','y'],['$Z_{\|\|}$','$Z_x$','$Z_y$'])):
        
            pylab.figure(iplane+1)
            # constructs all the files names (flat chamber impedances)
            list_files = glob(fil.replace("Zlong","Z{}*".format(plane)))
            
            for iZ,fileZ in enumerate(list_files):
            
                # read impedance and plot
                term = fileZ.split('W')[0].split('Z')[1] # term is something like "y0100"
                term_kind = term.split(plane)[1]
                if term_kind in ['ong','cst']:
                    sum_power_x = 0
                    sum_power_y = 0
                elif term_kind in ['dip','quad']:
                    sum_power_x = 1 if plane=='x' else 0
                    sum_power_y = 1 if plane=='y' else 0
                else:
                    sum_power_x = int(term_kind[0])+int(term_kind[2]) # a+c
                    sum_power_y = int(term_kind[1])+int(term_kind[3]) # b+d
                
                term_kind = term_kind.replace('ong','cst')
                freq,Z = read(fileZ)
                unit="[$\Omega$]"

                if opt.LOG:
                    logplot(freq,Z*(opt.SIZEX**(sum_power_x)*opt.SIZEY**(sum_power_y)),unit,Zlabel+' '+term_kind,", "+leg,colors[ind_color])
                else:
                    semilogplot(freq,Z*(opt.SIZEX**(sum_power_x)*opt.SIZEY**(sum_power_y)),unit,Zlabel+' '+term_kind,", "+leg,colors[ind_color])
           
                ind_color += 1
        
        
    for i,plane in zip(range(3),['l','x','y']):
        fig=pylab.figure(i+1)
        ax=fig.gca()
        set_fontsize(pylab.figure(i+1),'x-large')
        ax.legend(loc=0)
        l=ax.get_legend_handles_labels()
        set_fontsize(ax.get_legend(),max(20-2*len(l[0]),4))
        ax.grid()
    
        if opt.OUT!=None:
            #ax.set_xlim([1e3,1e12])
            #ax.set_ylim([1e-2,1e4])
            fig.savefig("Z"+plane+opt.OUT+".png")
            fig.savefig("Z"+plane+opt.OUT+".eps",format='eps')
            pylab.close(fig)


    if opt.OUT==None: pylab.show()
    
