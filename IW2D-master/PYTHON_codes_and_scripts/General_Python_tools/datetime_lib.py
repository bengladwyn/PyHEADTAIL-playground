#!/usr/bin/python

# library to use dates & times

import sys
import subprocess
pymod=subprocess.check_output("echo $PYMOD",shell=True).strip().decode()
if pymod.startswith('local'):
    py_numpy=subprocess.check_output("echo $PY_NUMPY",shell=True).strip().decode()
    sys.path.insert(1,py_numpy)
    py_matpl=subprocess.check_output("echo $PY_MATPL",shell=True).strip().decode()
    sys.path.insert(1,py_matpl)
import pylab,dateutil,pytz
from datetime import time,datetime,date
import matplotlib
import matplotlib.dates


def tt(f1): return datetime.strptime(f1,"%Y-%m-%d %H:%M:%S")


def tb(f1): return time.strftime("%Y-%m-%d %H:%M:%S",time.gmtime(f1))


def tb_(f1):

    ''' return a string in the format %Y-%m-%d %H:%M:%S corresponding to the pylab datetime (numeric format) f1'''
    tim=str(pylab.num2date(f1/86400.)); # string of the format 2012-10-09 13:40:00+00:00
    return tim.split('+')[0];


def set_axisdate(axi,axisname,timeinterval,tz=None):
    ''' set date on the axis 'axisname' ('x' or 'y') of the axes axi, using the timeinterval given.'''
    timeinterval=int(timeinterval);
    eval('axi.'+axisname+'axis_date(tz=tz)');
    eval('axi.'+axisname+'axis.set_major_locator(matplotlib.dates.SecondLocator(interval=timeinterval))');
    eval('axi.'+axisname+"axis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M:%S'))");

    return;


def plotdate(time,data,leg,dat,patcol,ylab,ax,logflag,lw=3.,plotevery=1,colr=None,ms=10.):

    ''' function plotting w.r.t time.
    Arguments (the first 8 ones are mandatory, the others are optional):
    - time: time data (datetime format),
    - data: y data,
    - leg: label for the legend,
    - dat: date (datetime format) at which the curve was taken,
    - patcol: pattern and/or color for the curve (e.g. '--xb'),
    - ylab: label for the y-axis,
    - ax: axes on which to plot (from e.g. init_figure),
    - logflag: 0 for linear plot, 1 for logarithmic plot in x / linear in y,
    2 for logarithmic plot in y / linear in x, 3 for log-log plots,
    - lw: linewidth of the curve,
    - plotevery: plot every "plotevery" points,
    - colr: re-define the color of the curve (rgb format for instance),
    - ms: size of the markers (crosses, circles, or other), if present. '''

    gmt=pytz.timezone('Europe/Amsterdam');

    if (colr==None):
        ax.plot_date(time[::plotevery],data[::plotevery],patcol,label=leg,linewidth=lw,ms=ms,mew=2.5,tz=gmt,color=colr);
    else:
        ax.plot_date(time[::plotevery],data[::plotevery],patcol,label=leg,linewidth=lw,ms=ms,mew=2.5,tz=gmt,color=colr);

    if (logflag>=2): ax.set_yscale('log');
    if (logflag==1)or(logflag==3): ax.set_xscale('log');

    ax.set_xlabel("Time on "+dat.strftime("%Y-%m-%d"));
    ax.set_ylabel(ylab);

    return;
