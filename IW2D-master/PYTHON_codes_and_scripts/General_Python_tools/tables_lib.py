#!/usr/bin/python

# library with tables and lists manipulation routines

import sys
import subprocess
pymod=subprocess.check_output("echo $PYMOD",shell=True).strip().decode()
if pymod.startswith('local'):
    py_numpy=subprocess.check_output("echo $PY_NUMPY",shell=True).strip().decode()
    sys.path.insert(1,py_numpy)
import numpy as np


def select_in_table(selection,table):
    ''' find indices ind of elements in table such that
    table[ind]=selection
    '''
    ind=-np.ones(len(selection),dtype=int);
    for isel,sel in enumerate(selection):
        for itab,tab in enumerate(table):
            if tab==sel:
                if (ind[isel]!=-1)or(any(ind-itab==0)): print("Pb in select_in_table: twice same selection!")
                else: ind[isel]=itab;

    if any(ind==-1): print("Pb in select_in_table: some non affected in selection!")

    return ind;


def intersect(tab1,tab2):
    ''' build array that is the intersection between two 1D arrays tab1 and tab2
    '''
    import pylab

    inter=[];ind1=[];ind2=[];
    # choose smallest of the 2 arrays
    if (len(tab1)<len(tab2)): i=1;
    else: i=2;

    for it,t in eval('enumerate(tab'+str(i)+')'):
        if (t in eval('tab'+str(3-i))):
            inter.append(t);
            exec('ind'+str(i)+'.append(it)');
            exec('ind'+str(3-i)+'.append(pylab.mlab.find(tab'+str(3-i)+'==t))');

    return np.array(inter),np.array(ind1).reshape(-1),np.array(ind2).reshape(-1);


def count_between(xdata,ydata,ythreshold,xmin,xmax):
    ''' count the number of data points in array 'ydata' above 'ythreshold' and such
    that their corresponding absissae in 'xdata' are between 'xmin' and 'xmax'.
    '''
    import pylab

    indy=pylab.mlab.find(ydata>=ythreshold);
    indx=pylab.mlab.find((xdata-xmax)*(xdata-xmin)<=0);
    inter,a,b=intersect(indx,indy);
    return len(inter);


def diff(data1,data2):
    ''' compare two sets of data'''
    if (len(data1)!=0):
        d1=np.max(np.abs((data1[data1!=0]-data2[data1!=0])/data1[data1!=0]));
        if (len(data1[data1==0])>0): d2=np.max(np.abs((data1[data1==0]-data2[data1==0])));
        else: d2=0.;
        d3=np.max(np.abs(data1-data2))/np.abs(np.average(data1));
        #d3=np.max(np.max(np.abs(data1-data2),axis=0)/np.abs(np.average(data1,axis=0))); # does not work (nan when one mean is zero))
    else: d1=0.;d2=0.;d3=0.;
    return d1,d2,d3;


def diffshape(data1,data2):
    ''' compare shape of two sets of data'''
    return (data1.shape==data2.shape);


def complementary(tab1,tab2):
    '''select the complementary of tab1 in tab2 (i.e. those in tab2 but not in tab1)
    '''
    tab=[];
    for t in tab2:
        if not (t in tab1): tab.append(t);

    return np.array(tab);


def sort_and_delete_duplicates(tab,tolerance=1e-8):

    ''' sort a table in ascending order and delete its duplicates
    tolerance is the relative tolerance accepted between 2 numbers (if their
    difference is lower than that, in relative, they are considered equal)
    '''

    tab2=np.sort(tab);
    #print len(tab),tab2[-10:]
    # delete duplicates
    indnotdup=np.where(np.diff(tab2)/np.abs(tab2[1:]+tolerance)>tolerance);
    indnotdup=np.concatenate((indnotdup[0],np.array([len(tab2)-1])));
    #print len(indnotdup),indnotdup[-10:],indnotdup[:10],tab2[indnotdup[-10:]],tab2[indnotdup[:10]]

    return tab2[indnotdup];


def create_list(a,n=1):
    ''' if a is a scalar, return a list containing n times the element a
    otherwise return a
    '''
    from copy import deepcopy

    #try:
    #   l=len(a); # to check if it's already a list
    if not(np.isscalar(a)):
        pass;
    else:
    #except TypeError:
        b=[];
        for i in range(n): b.append(a);
        a=deepcopy(b);

    return a;


def create_list_for_figure(a,n=1):
    ''' same as create_list but done differently, with a test on the length
    (because the one above fails for objects like axes and figures)
    also, no "deepcopy" (means that if all elements are duplicate,
    any change affecting one of them will affect the others)'''

    try:
        l=len(a); # to check if it's already a list
        b=a
    except TypeError:
        b=[];
        for i in range(n): b.append(a);

    return b;


def select_and_average(data,select_table,select_value,flagmax=False,flagvar=False):

    ''' select in list 'data' the values according to indices such that select_table(indices)=select_value
    then take out nan and compute average along the data lines.
    number of elements in list 'data' should be the same as number of elements of array 'select_table'

    if flagmax=True, also compute the maximum of the absolute value (along the lines)
    if flagvar=True, also compute the standard deviation (along the lines)'''

    # select according to select_value
    ind=pylab.mlab.find(select_table==select_value);

    if (len(ind)>0):

        data1=[data[i] for i in ind];
        # take out nan and infinities
        data1=np.array([el for el in data1 if isfinite(el).any()]);
        # average
        data_av=np.average(data1,axis=0);

        if flagmax:
            # maximum of absolute values along the lines of data
            data_max=np.max(np.abs(data1),axis=0);
        else: data_max=[];

        if flagvar:
            # standard deviation along the lines of data
            data_var=np.sqrt(np.var(data1,axis=0));
        else: data_var=[];

        return data_av,data_max,data_var;

    else: return nan,nan,nan;


def add_element_in_array(table,i,elem):

    ''' intercalate one element after place i in the array table.
    It means that elem takes the place at the index i+1 (between the previously i and
    i+1 elements of table) and that the right part of table is shifted to
    the right by one unit.
    Note that it also works with i=-1 (elem inserted at the beginning of the table)
    or with i=len(table)-1 (insertion at the end).
    '''

    return np.hstack((table[:i+1],elem,table[i+1:]));


def add_line_in_2d_array(table,i,line):

    ''' intercalate one line after line i in the 2D array table.
    It means that e'line' takes the place at the index i+1 (between the previously i and
    i+1 lines of table) and that the bottom part of table is shifted to
    the down by one unit.
    Note that it also works with i=-1 (eline inserted at the beginning of the table)
    or with i=len(table[:,0])-1 (insertion at the end).
    '''

    return np.concatenate((table[:i+1,:],line.reshape((1,-1)),table[i+1:,:]));


def max_diff_rel_complex_array(tab1,tab2,criterion='abs'):

    ''' maximum and index of the maximum of the difference between 2 complex
    arrays 'tab1' and 'tab2', in relative w.r.t. 'tab1', according to
    'criterion' that can be:
    - 'abs' -> maximum in absolute value,
    - 'real' -> maximum for the real part only,
    - 'imag' -> maximum for the imaginary part only,
    - 'realimag' -> maximum for both real and imaginary parts.
    '''

    if criterion=='abs':
        maxerr=np.max(np.abs((tab1-tab2)/tab1));
        ierr=np.argmax(np.abs((tab1-tab2)/tab1));

    elif criterion=='real':
        maxerr=np.max(np.abs(np.real(tab1-tab2)/np.real(tab1)));
        ierr=np.argmax(np.abs(np.real(tab1-tab2)/np.real(tab1)));

    elif criterion=='imag':
        maxerr=np.max(np.abs(np.imag(tab1-tab2)/np.imag(tab1)));
        ierr=np.argmax(np.abs(np.imag(tab1-tab2)/np.imag(tab1)));

    elif criterion=='realimag':
        maxerr_real=np.max(np.abs(np.real(tab1-tab2)/np.real(tab1)));
        maxerr_imag=np.max(np.abs(np.imag(tab1-tab2)/np.imag(tab1)));
        maxerr=max(maxerr_real,maxerr_imag);
        if (maxerr_real>maxerr_imag): ierr=np.argmax(np.abs(np.real(tab1-tab2)/np.real(tab1)));
        else: ierr=np.argmax(np.abs(np.imag(tab1-tab2)/np.imag(tab1)));

    return maxerr,ierr;


def max_diff_abs_complex_array(tab1,tab2,criterion='abs'):

    ''' maximum and index of the maximum of the difference between 2 complex
    arrays 'tab1' and 'tab2', in absolute, according to
    'criterion' that can be:
    - 'abs' -> maximum in absolute value,
    - 'real' -> maximum for the real part only,
    - 'imag' -> maximum for the imaginary part only,
    - 'realimag' -> maximum for both real and imaginary parts.
    '''

    if criterion=='abs':
        maxerr=np.max(np.abs(tab1-tab2));
        ierr=np.argmax(np.abs(tab1-tab2));

    elif criterion=='real':
        maxerr=np.max(np.abs(np.real(tab1-tab2)));
        ierr=np.argmax(np.abs(np.real(tab1-tab2)));

    elif criterion=='imag':
        maxerr=np.max(np.abs(np.imag(tab1-tab2)));
        ierr=np.argmax(np.abs(np.imag(tab1-tab2)));

    elif criterion=='realimag':
        maxerr_real=np.max(np.abs(np.real(tab1-tab2)));
        maxerr_imag=np.max(np.abs(np.imag(tab1-tab2)));
        maxerr=max(maxerr_real,maxerr_imag);
        if (maxerr_real>maxerr_imag): ierr=np.argmax(np.abs(np.real(tab1-tab2)));
        else: ierr=np.argmax(np.abs(np.imag(tab1-tab2)));

    return maxerr,ierr;
