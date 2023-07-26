#!/usr/bin/python

# library with parser routines (for command line options)

from optparse import OptionParser

def bunch_parse(option, opt_str, value, parser):

    ''' special bunch numbers parser (parse things like -b 1 -b 5 -b 7:9 -b 12:18:2) '''
    #print parser.rargs;
    value=getattr(parser.values,option.dest);
    #print value;
    if (value==None): value=[];
    arg=parser.rargs[0];
    if ":" not in arg: value.append(int(arg))
    else:
        s=arg.split(":")
        if (len(s)>=3):
            value.extend(list(range(int(s[0]),int(s[1])+1,int(s[2]))));
        else:
            value.extend(list(range(int(s[0]),int(s[1])+1)));
    del parser.rargs[:1]
    setattr(parser.values, option.dest, value)


def multistring_parse(option, opt_str, value, parser):

    ''' parse an option (variable number of string inputs) '''
    value=[];
    #print parser.rargs
    l=0;
    while (len(parser.rargs)>=(l+1))and(not(parser.rargs[l].startswith("-"))):
        value.append(parser.rargs[l]);
        l+=1;
    #print value;
    del parser.rargs[:l];
    setattr(parser.values, option.dest, value)
    #print parser.rargs,parser.values
    return;


def multifloatint_parse(option, opt_str, value, parser):

    ''' parse an option (variable number of float or integer inputs) '''
    value=[];
    #print parser.rargs
    l=0;
    while (len(parser.rargs)>=(l+1))and(not(parser.rargs[l].startswith("-"))):
        val=float(parser.rargs[l]);
        if val.is_integer(): val=int(parser.rargs[l]);
        value.append(val);
        l+=1;
    #print value;
    del parser.rargs[:l];
    setattr(parser.values, option.dest, value)
    #print parser.rargs,parser.values
    return;
