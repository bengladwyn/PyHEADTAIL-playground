#!/usr/bin/python

# library with string manipulation routines

import sys
import subprocess
pymod=subprocess.check_output("echo $PYMOD",shell=True).strip().decode()
if pymod.startswith('local'):
    py_numpy=subprocess.check_output("echo $PY_NUMPY",shell=True).strip().decode()
    sys.path.insert(1,py_numpy)
import numpy as np


def fortran_str(numb):
    ''' convert a float to a Fortran-like string (with 'D' instead of 'e') '''
    return str(numb).replace('e','D');


def float_to_str(f):

    ''' convert a float into a string, replacing  '.' by 'p',
    and taking care of useless '0' '''
    s=str(f);
    if (s.endswith('.0')): s=s.rstrip("0").rstrip('.');

    return s.replace('.','p')


def invert_selection(namestot,namespart):
    ''' select all names in namestot (list of strings) that are not in namespart '''
    namesnew=[];
    for name in namestot:
        if not (name in namespart): namesnew.append(name);

    return namesnew;


def takeout_common(listname,flagcommon=False,nslash=2):

    ''' Takes out all the "words" that are common to all the
    names in listname:
    "words" are defined as being between two underscores "_" (or
    between the beginning and a underscore, or an underscore and the end).
    Gives in output a new list of names without those common words.

    Note: if there are "/"  characters, we take into account only the part
    of the name that is between the nslash-1 and nslash ones (counting from
    the end)
    for instance if a name is "blabla/toto/LHC_36b_csi4/data_prt.dat"
    and nslash=2, we take only "LHC_36b_csi4"

    if flagcommon=True, also output the common part of all names '''

    listname1=[];
    for iname,name in enumerate(listname):
        name1=name.split("/")
        if (len(name1)==1): listname1.append(name1[0]);
        elif (len(name1)>(nslash-1)): listname1.append(name1[-nslash]);
        else:
            print("Pb with nslash");sys.exit();

    # find the common words
    common=listname1[0].split("_")
    for name in listname1:
        tmp=name.split("_")
        for i,com in enumerate(common):

            flag=False;
            for item in tmp:
                flag=(flag or (com==item));

            if (not flag): common.pop(i);

    # take out the common words
    listnew=[];
    for name in listname1:
        for com in common:
            tmp=name.replace(com,"").replace("_"," ");
            name=tmp;

        listnew.append(name);


    if (not(flagcommon)): return listnew;
    else: return listnew,common;


def split_and_takeout_spaces(name):

    ''' split (according to spaces) and take out spaces from a string
    return a list with each "words" '''
    s=name.split(" ")
    slen=np.array([len(k) for k in s]);
    ind=np.where(slen>0);ind=ind[0];

    return [s[k] for k in ind];


def get_nice_string(listin):

    ''' to output a string from a list (for printing purposes) '''
    return " ".join( str(x) for x in listin)


def takeout_spaces(listname):

    ''' take out most of the spaces in each term of a list of strings
    (obtained from the 'takeout_common' routine) '''
    listnew=[];
    for name in listname:
        split_name=split_and_takeout_spaces(name);
        split_name=[nm+' ' for nm in split_name];
        newname=''.join(split_name);newname=newname[:-1];
        listnew.append(newname);

    return listnew;


def find_ind_names(namesref,names):
    ''' given 2 lists of strings namesref and names, find the indices ind
    such that names[ind]=namesref
    '=' means here that the longest of the two strings begins
    with the smallest.'''

    ind=-np.ones(len(namesref),dtype=int);
    for inamer,namer in enumerate(namesref):
        for iname,name in enumerate(names):

            if (len(namer)>len(name)):
                if namer.startswith(name):
                    if (ind[inamer]!=-1)or(any(ind-iname==0)): print("Pb in find_ind_names: twice same name!")
                    else: ind[inamer]=iname;
            else:
                if name.startswith(namer):
                    if (ind[inamer]!=-1)or(any(ind-iname==0)): print("Pb in find_ind_names: twice same name!")
                    else: ind[inamer]=iname;

    if any(ind==-1): print("Pb in find_ind_names: some non affected names!")

    return ind;
