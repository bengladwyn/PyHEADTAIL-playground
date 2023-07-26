#!/usr/bin/python

# some basic physical constants

def electron_param():

    ''' gives basic electron parameters
    Output: e [C], mass [kg], speed of light [m/s], rest energy [J] '''

    e=1.602176487e-19; # elementary charge
    m0=9.10938e-31; # electron mass in kg
    c=299792458; # speed of light
    E0=m0*c**2 # rest energy

    return e,m0,c,E0


def proton_param():

    ''' gives basic proton parameters
    Output: e [C], mass [kg], speed of light [m/s], rest energy [J] '''

    e=1.602176487e-19; # elementary charge
    m0=1.6726216e-27; # proton mass in kg
    c=299792458; # speed of light
    E0=m0*c**2 # rest energy

    return e,m0,c,E0


def Pb54_param():

    ''' gives basic proton parameters
    Output: e [C], mass [kg], speed of light [m/s], rest energy [J] '''

    e=1.602176487e-19; # elementary charge
    m0= (207.947/208.)*1.6726216e-27; # Pb54 nucleon mass in kg
    c=299792458; # speed of light
    E0=m0*c**2 # rest energy

    return e,m0,c,E0


def Pb54_ion_param():

    ''' gives basic proton parameters
    Output: e [C], mass [kg], speed of light [m/s], rest energy [J] '''
    Z = 54
    e = Z*1.602176487e-19; # elementary charge
    m0 = 207.947*1.6726216e-27; # Pb54 nucleon mass in kg
    c = 299792458; # speed of light
    E0 = m0*c**2 # rest energy

    return e,m0,c,E0


def Pb82_param():

    ''' gives basic Pb82+ parameters
    Output: Z*e [C], mass [kg], speed of light [m/s], rest energy [J] '''
    Z = 82
    e = 1.602176487e-19 # elementary charge
    c = 299792458 # speed of light
    E0 = 193.68715*1e9*e # Pb (isotope 208) rest energy in J
    m0 = E0/(c**2) # mass in kg

    return Z*e,m0,c,E0
