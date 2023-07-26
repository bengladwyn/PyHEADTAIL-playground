#!/usr/bin/python

# library with plotting routines

import sys
import subprocess
pymod=subprocess.check_output("echo $PYMOD",shell=True).strip().decode()
if pymod.startswith('local'):
    py_numpy=subprocess.check_output("echo $PY_NUMPY",shell=True).strip().decode()
    sys.path.insert(1,py_numpy)
    py_matpl=subprocess.check_output("echo $PY_MATPL",shell=True).strip().decode()
    sys.path.insert(1,py_matpl)
import numpy as np
import pylab,os
import matplotlib
import matplotlib.text as text
import time as ti


def set_fontsize(obj,size):

    ''' set all fontsize to 'size' in the (pylab) object 'obj' (can be figure, axes, legend labels...) '''
    for o in obj.findobj(text.Text): o.set_fontsize(size);


def init_figure(axes=None,figsize=(16, 12)):
    ''' initialize a figure and its axes
    optional arguments:
    - figsize: total figure size
    - axes: for instance [0.15,0.15,0.8,0.8]: values between 0 and 1,
    in units of the total figure size: x & y coordinates of bottom left coorner,
    then horizontal and vertical sizes '''

    pylab.rcParams['xtick.major.pad']='10';
    pylab.rcParams['ytick.major.pad']='10';
    fig=pylab.figure(figsize=figsize,facecolor='w', edgecolor=None)
    if (axes!=None): pylab.axes(axes)
    ax=fig.gca()

    return fig,ax


def sciy(ax):
    ''' transform vertical tick labels from axes ax, into scientific notation (e.g. 1.5 10^-4)
    unless they are between 0.01 & 100 (in abs. value) '''
    try:
        ax.ticklabel_format(style='sci', scilimits=(-2,2),axis='y')
    except AttributeError:
        pass;

def scix(ax):
    ''' transform horizontal tick labels from axes ax, into scientific notation (e.g. 1.5 10^-4)
    unless they are between 0.01 & 100 (in abs. value) '''
    try:
        ax.ticklabel_format(style='sci', scilimits=(-2,2),axis='x')
    except AttributeError:
        pass;


def end_figure(fig,ax,save=None,legpos=0,fontsize=30,grid=True,legfontsize=None):
    ''' finalize a figure fig and its axes ax, to make it look 'nice'.
    Optional arguments:
    - save: if not None, should contain the name of the file where to save the plots
    (then it does not print on the screen but save in files [save].png and [save].eps).
    Otherwise plots are simply printed on screen,
    - legpos: position of the legend: from matplotlib help:
        best position   0
        upper right     1
        upper left      2
        lower left      3
        lower right     4
        right           5
        center left     6
        center right    7
        lower center    8
        upper center    9
        center          10
    It can also be a tuple with the coordinates of the up-left corner of the legend box.
    - fontsize: general fontsize of all texts in the plots. If legfontsize is None,
    legend text is with a fontsize equal to 'fontsize-2*(number of legends)' (but not less than 10),
    - grid: True or False to switch on/off the grid,
    - legfontsize: fontsize of legend text. Overrides the rule above if not None. '''

    font = {#'family' : 'normal',
            #'weight' : 'bold',
            'size'   : fontsize}
    matplotlib.rcdefaults()
    matplotlib.rc('font', **font)
    if (ax.xaxis.get_scale()=='linear'): scix(ax);
    if (ax.yaxis.get_scale()=='linear'): sciy(ax);

    l=ax.get_legend_handles_labels();
    if (len(l[0])>0):
        ax.legend(loc=legpos, prop=font);
        if (legfontsize==None): set_fontsize(ax.get_legend(),max(fontsize-2*len(l[0]),10));
        else: set_fontsize(ax.get_legend(),legfontsize);

    #ax.ticklabel_format(style='sci', scilimits=(0,0),axis='x')
    #ax.ticklabel_format(style='sci', scilimits=(0,0),axis='y')

    ax.xaxis.labelpad=12;
    ax.yaxis.labelpad=12;

    #set_fontsize(fig,'xx-large');
    #set_fontsize(fig,fontsize);

    if grid: ax.grid();

    if (save!=None)and(len(save)>0):
        # save figure (.png and .eps files) and close it.
        # "save" contains the root of the file names
        fig.savefig(save+".png")
        fig.savefig(save+".eps",format='eps')
        pylab.close(fig);

    return;


def plot(x,data,leg,patcol,ylab,ax,logflag,lw=3.,xlab="Number of turns",plotevery=1,colr=None,ms=10.):

    ''' plotting function.
    Arguments (the first 7 ones are mandatory, the others are optional):
    - x: x data,
    - data: y data,
    - leg: label for the legend,
    - patcol: pattern and/or color for the curve (e.g. '--xb'),
    - ylab: label for the y-axis,
    - ax: axes on which to plot (from e.g. init_figure),
    - logflag: 0 for linear plot, 1 for logarithmic plot in x / linear in y,
    2 for logarithmic plot in y / linear in x, 3 for log-log plots,
    - lw: linewidth of the curve,
    - xlab: label for the x-axis,
    - plotevery: plot every "plotevery" points,
    - colr: re-define the color of the curve (rgb format for instance)  ,
    - ms: size of the markers (crosses, circles, or other), if present. '''

    if (len(data)>0):
        if (colr is None):
            if (logflag==0): ax.plot(x[::plotevery],data[::plotevery],patcol,label=leg,linewidth=lw,ms=ms,mew=2.5);
            elif (logflag==1): ax.semilogx(x[::plotevery],data[::plotevery],patcol,label=leg,linewidth=lw,ms=ms,mew=2.5);
            elif (logflag==2): ax.semilogy(x[::plotevery],data[::plotevery],patcol,label=leg,linewidth=lw,ms=ms,mew=2.5);
            elif (logflag==3): ax.loglog(x[::plotevery],data[::plotevery],patcol,label=leg,linewidth=lw,ms=ms,mew=2.5);
        else:
            if (logflag==0): ax.plot(x[::plotevery],data[::plotevery],patcol,label=leg,linewidth=lw,ms=ms,mew=2.5,color=colr);
            elif (logflag==1): ax.semilogx(x[::plotevery],data[::plotevery],patcol,label=leg,linewidth=lw,ms=ms,mew=2.5,color=colr);
            elif (logflag==2): ax.semilogy(x[::plotevery],data[::plotevery],patcol,label=leg,linewidth=lw,ms=ms,mew=2.5,color=colr);
            elif (logflag==3): ax.loglog(x[::plotevery],data[::plotevery],patcol,label=leg,linewidth=lw,ms=10.,mew=2.5,color=colr);

    ax.set_xlabel(xlab);
    ax.set_ylabel(ylab);


def cmap(kys):
    ''' define a color map'''
    maps=[]
    for u in range(kys):
        maps.append((round(random.uniform(0,1),1),
                     round(random.uniform(0,1),1),
                     round(random.uniform(0,1),1)))
    return maps


def set_legend_fontsize(fig,ax,font='xx-large'):

    ''' puts a legend in axes 'ax' of figure 'fig', with fontsize 'font' '''
    ax.legend(loc=0);
    set_fontsize(fig,font);

    return;


def plot2D(data,nxmin,nxmax,nymin,nymax,xlab,ylab,tit,ax,colorlabel=None,colorlim=None,fig=None):

    ''' 2D color plot of data between nxmin, nxmax, nymin, nymax.
    xlab and ylab are the axes labels, tit the title,
    ax the axes to do the plot on.

    NOTE: data is of the format [iy,ix] (i.e. y axis in first dimension)

    if colorlabel or colorlim is not None, the color scale is shown.
    colorlabel is a label to be put next to the color scale, colorlim a list
    with 2 floats giving the limits of the scale (if None it is automatic). '''

    # NOTE 2: imshow plots the matrix "as it looks like when you read it" ->
    # this means that the lines with higher indices are at the BOTTOM of the graph
    # (contrary to what would be more logical, i.e. higher indices at the top). Columns
    # follow a more logical behaviour, i.e. higher indices on the right.
    # Therefore, one needs to flip "upside-down" the data in the vertical direction,
    # to reverse the orders of the lines -> this is done by the flipud command applied below.

    import matplotlib.pyplot as plt;
    cmap=plt.get_cmap('jet');


    im=ax.imshow(np.flipud(data), aspect='auto', extent=[nxmin,nxmax,nymin,nymax],cmap=cmap)
    # Plotting contour lines
    #ax.contour(data, aspect='auto', extent=[nxmin,nxmax,nymin,nymax])
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(tit)

    if (colorlabel!=None)or(colorlim!=None):
        # color scale
        if (fig==None): fig=plt.gcf();
        bar=fig.colorbar(im);
        if (colorlabel!=None): bar.set_label(colorlabel);
        if (colorlim!=None): bar.set_clim(colorlim);


def build_colors(ncol,randomize=False):

    ''' define ncol different colors (rgb format, for plots).
    If randomize==True, they are randomized (always in the same way, random permutation is hard-coded)'''

    if ncol==1: n=1;
    else: n=int(np.ceil((ncol+1)**(1./3.)));#print "build_colors: ncol=",ncol,", n=",n
    cl=np.linspace(0,0.9,n);col=[];
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if ((cl[i]!=0.9)or(cl[-j-1]!=0.9))or(cl[k]!=0.9): col.append((cl[i],cl[-j-1],cl[k]));

    c=np.array(col);

    if randomize:
        # apply some (fixed) random permutation (from np.random.permutation(len(col)) )
        #if n==1: ind=np.array([0]);
        #elif n==2: ind=np.array([5, 2, 0, 6, 3, 1, 4, 7]);
        #elif n==3: ind=np.array([24,13,0,4,16,11,19,26,18,14,22,21,15,23,12,17,2,6,5,7,10,8,3,9,25,1,20])
        #elif n==4: ind=np.array([49,27,14,26,5,45,51,62,20,52,57,12,31,61,25,63,22,43,4,21,24,42,0,39,11,23,56,7,35,28,41,3,40,32,33,8,37,50,53,55,58,10,60,9,6,44,29,54,17,46,18,1,15,59,47,2,16,30,38,36,48,19,13,34])
        if n==1: ind=np.array([0]);
        elif n==2: ind=np.array([5, 3, 1, 4, 2, 6, 0]);
        elif n==3: ind=np.array([4,23,12,22,20,13,1,16,11,0,19,3,21,8,5,18,25,6,17,15,9,7,14,24,10,2])
        elif n==4: ind=np.array([28,46,9,23,49,41,0,4,43,56,10,16,45,60,40,31,36,48,42,12,30,32,34,53,57,8,37,26,1,6,50,38,20,39,24,22,5,15,55,61,7,2,59,33,47,29,17,21,27,3,13,14,52,54,51,35,62,18,11,19,25,58,44])
        else: ind=np.arange(len(c));

        c=c[ind];

    return c


def make_movie(filename,rootname,flagrm=True):

    ''' make a movie (animated gif) from files rootname*.png
    The movie is entitled filename
    If flgrm==True, rm all files rootname*.png and rootname*.eps afterwards '''

    print('Making movie ',filename,' - this may take a while')
    os.system("convert "+rootname+"*.png "+filename)
    # then erase all initial files
    if flagrm:
        os.system("rm -f "+rootname+"*.png")
        os.system("rm -f "+rootname+"*.eps")


def plot_save_hist(x,data,leg,ax,fig,i,xlim=None,ylim=None,tit='',ylab="Bunch position",legposition='lower right',xlab="25ns-slot number"):

    ''' plot and save the histogram of data vs. x (save in '_tmp'+str(i)+'.png' and '_tmp'+str(i)+'.eps')
    - leg: label for the legend,
    - ax: axes where to plot,
    - fig: figure objectto which ax belongs,
    - i: integer to be put in the filename,
    - xlim: 2 numbers with x-axis limits (default=automatic),
    - ylim: 2 numbers with y-axis limits (default=automatic),
    - tit: title of the graph,
    - ylab: label for the y-axis,
    - legposition: position of the legend (see function end_figure),
    - xlab: label for the x-axis.
    '''

    ax.cla() # clear the plot
    ax.bar(x,data,facecolor='k',label=leg,edgecolor='k')
    ax.set_xlabel(xlab);
    ax.set_ylabel(ylab);
    ax.set_title(tit);
    if (xlim!=None): ax.set_xlim(xlim);
    if (ylim!=None): ax.set_ylim(ylim);
    fname = '_tmp%d'%i
    end_figure(fig,ax,save=fname,legpos=legposition);


def fillplot_percent(x,y1,y2,lab,leg,col,ax):
    ''' fill plot of percentages between y1 and y2, with x-label in 'lab', legend
    label in 'leg', and color 'col'.
    ax are the axes to plot on. '''

    ax.fill_between(x,y1*100,y2*100,where=(y2<=y1),color=col,label=leg,lw=2.5)
    ax.set_xlabel(lab);
    ax.set_ylabel("Percent of the total");
