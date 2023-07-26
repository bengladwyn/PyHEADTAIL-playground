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
    # Note: if a range is given for several coordinates (e.g. dipx and quadx)
    # then they are both displaced together.
    parser.add_option("-f", "--file",action="append",
                      help="Specify the file name (beginning with Zlong - the other impedance files Zxdip, Zydip, Zl1000, etc. are selected automatically)",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-g", "--legend",action="append",
                      help="Specify the legend for each file",
                      metavar="LEG", default=None,dest="LEG")
    parser.add_option( "--xdip",default=(0.,0.),dest="DIPX",
                      help="Specify a (geometric) dipolar horizontal displacement, in mm (two values -*> provide a range)",
                      metavar="DIPX",type=float,nargs=2)
    parser.add_option( "--xquad",default=(0.,0.),dest="QUADX",
                      help="Specify a (geometric) quadrupolar horizontal displacement, in mm (two values -> provide a range)",
                      metavar="QUADX",type=float,nargs=2)
    parser.add_option( "--ydip",default=(0.,0.),dest="DIPY",
                      help="Specify a (geometric) dipolar vertical displacement, in mm (two values -*> provide a range)",
                      metavar="DIPY",type=float,nargs=2)
    parser.add_option( "--yquad",default=(0.,0.),dest="QUADY",
                      help="Specify a (geometric) quadrupolar vertical displacement, in mm (two values -> provide a range)",
                      metavar="QUADY",type=float,nargs=2)
    parser.add_option( "--freq",default=1.,dest="FREQ",
                      help="Specify a single frequency where to compute the impedance (in MHz)",
                      metavar="FREQ",type=float)
    parser.add_option( "-n","--norm",default=False,dest="NORM",
                      help="Specify if you want to normalize transverse impedances by the transverse displacement (default=False). WARNING: this gives problems at zero displacement.",
                      metavar="NORM",action="store_true")
    parser.add_option( "-l","--length",default=None,dest="LEN",
                      help="Specify if you want to normalize the impedances by the length of the device (default=None, i.e. one does not normalize). You should specify then the length.",
                      metavar="LEN",type=float)
    parser.add_option("-o", "--outfile",
                      help="Specify the output file name if you want to output the plots in files (will create Zvsdispl[this_option].png and the same in .eps format) (by default: print plots on screen -> CANNOT WORK IF YOU USE MATPLOTLIB WITH GUI DISABLED)",
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


def get_halfgap_from_input_(filename):
    # read impedance input file and get the half-gap
    with open(filename) as f:
        lines = f.readlines()
    hg_input_string = "Layer 1 inner half gap in mm:"
    return [float(l.split(hg_input_string)[1])*1e-3 for l in lines if hg_input_string in l][0]    


def plot(displ,Z,xlabel,unit,leg1,leg2,col):
    # plot impedance (real and imag parts) with labels, legend and color
    legr="Re["+leg1+"]"
    legi="Im["+leg1+"]"
    pylab.plot(displ,Z.real,col if 'x' in col else col+'-',label=legr+leg2,lw=2.5,ms=10.,mew=2.5)
    pylab.plot(displ,Z.imag,col.replace('x','+') if 'x' in col else col+'--',label=legi+leg2,lw=2.5,ms=10.,mew=2.5)
    pylab.xlabel("{} displacement [mm]".format(xlabel))
    pylab.ylabel("Z "+unit)


def get_power_from_term(term,plane,norm=False):
    # get a list of powers [a,b,c,d] (xdip, ydip, xquad, yquad) from a string
    # of the form 0120 or 'cst' or 'ong' or 'dip' or 'quad', and the plane under
    # consideration.
    # If norm is True, transverse terms are normalized by dipolar displacement
    # along the same plane (WARNING: this gives problems at zero displacement).
    term_kind = term.split(plane)[1]
    if term_kind in ['ong','cst']:
        power = [0, 0, 0, 0]
    elif term_kind=='dip':
        power = [1, 0 , 0, 0] if plane=='x' else [0, 1 , 0, 0]
    elif term_kind=='quad':
        power = [0, 0 , 1, 0] if plane=='x' else [0, 0 , 0, 1]
    else:
        power = [int(e) for e in term_kind]
        
    if norm and plane=='x':
        power[0] -= 1
    elif norm and plane=='y':
        power[1] -= 1
    return power


if __name__ == "__main__":
    opt,args=parsse()
       
    displ = []
    ndispl = 100
    indices = []
    terms = ['Dip. x', 'Dip. y', 'Quad. x', 'Quad. y']
    for i,(bottom,top) in enumerate([opt.DIPX,opt.DIPY,opt.QUADX,opt.QUADY]):
        displ.append(np.linspace(bottom,top,ndispl+1))
        if len(set(displ[-1]))!=1:
            indices.append(i)
    
    # one figure per plane
    for i in range(3):
        init_figure()
    
    
    color = ['b','xb','r','xr','g','xg','m','xm','k','xk','c','xc','y','xy']
    ind_color_old = 0
    
    for j,fil in enumerate(opt.FILE):

        # string for legend
        if opt.LEG==None: leg=fil.replace("ZlongW","",1).replace("_"," ")
        else: leg=opt.LEG[j]
        
        for iplane,(plane,Zlabel) in enumerate(zip(['l','x','y'],['$Z_{\|\|}$','$Z_x$','$Z_y$'])):
        
            pylab.figure(iplane+1)
            # constructs all the files names (flat chamber impedances)
            list_files = glob(fil.replace("Zlong","Z{}*".format(plane)))
            
            # first get the maximum order of all the terms
            max_order = 0
            for iZ,fileZ in enumerate(list_files):
                term = fileZ.split('W')[0].split('Z')[1] # term is something like "y0100"
                power = get_power_from_term(term,plane,norm=opt.NORM)
                if sum(power)>max_order:
                    max_order = sum(power)
            
            ind_color = ind_color_old
            for order in range(max_order+1):
                
                Zsum = np.zeros(ndispl+1,dtype=complex)
                for iZ,fileZ in enumerate(list_files):
                
                    # read impedance and plot
                    term = fileZ.split('W')[0].split('Z')[1] # term is something like "y0100"
                    power = get_power_from_term(term,plane,norm=opt.NORM)
                    
                    if sum(power)<=order:
                        all_freq,all_Z = read(fileZ)
                        if len(all_freq)>0:
                            Z = np.interp(1e6*opt.FREQ,all_freq,all_Z)
                            Zsum += Z*np.prod([(d*1e-3)**p for d,p in zip(displ,power)],axis=0)
                        else:
                            print "WARNING! file {} is empty (computation is not finished, most probably)".format(fileZ)
                
                if plane in ['x','y']:
                    unit="[$\Omega$]" if (not opt.NORM and not opt.LEN) else "[$\Omega/$m$^2$]" if (opt.NORM and opt.LEN) else "[$\Omega/$m]"
                else:
                    unit="[$\Omega$]" if not opt.LEN else "[$\Omega/$m]"
                plot(displ[indices[0]],Zsum,terms[indices[0]],unit,Zlabel,", up to order {} {}".format(order,leg),color[ind_color])
            
                ind_color += 1 
        ind_color_old = ind_color
        
    for i,plane in zip(range(3),['l','x','y']):
        fig=pylab.figure(i+1)
        ax=fig.gca()
        set_fontsize(pylab.figure(i+1),'x-large')
        ax.legend(loc=0)
        l=ax.get_legend_handles_labels()
        set_fontsize(ax.get_legend(),max(20-2*len(l[0]),4))
        ax.grid()
    
        if opt.OUT!=None:
            fig.savefig('Z'+plane+opt.OUT+".png")
            fig.savefig('Z'+plane+opt.OUT+".eps",format='eps')
            pylab.close(fig)


    if opt.OUT==None: pylab.show()
    
