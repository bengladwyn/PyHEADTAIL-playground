#!/usr/bin/python

# library for computation of impedance & wake models

import sys,subprocess,filecmp
# define user (for bsub command) (environment variable USER should be previously defined)
user=subprocess.check_output("echo $USER",shell=True).strip().decode()
if (len(user)>0): user_option=" -u "+user;
else: user_option='';
# define path for Yokoya factors file (environment variable YOK_PATH should be previously defined
#- see script_configure_bashrc.sh)
Yokoya_path=subprocess.check_output("echo $YOK_PATH",shell=True).strip().decode()
# define path for fiIW2D executables (environment variable IW2D_PATH should be previously defined
#- see ../../ImpedanceWake2D/script_configure_bashrc.sh)
IW2D_path=subprocess.check_output("echo $IW2D_PATH",shell=True).strip().decode()
# import local libraries if needed
pymod=subprocess.check_output("echo $PYMOD",shell=True).strip().decode()
if pymod.startswith('local'):
    py_numpy=subprocess.check_output("echo $PY_NUMPY",shell=True).strip().decode()
    sys.path.insert(1,py_numpy)
    py_matpl=subprocess.check_output("echo $PY_MATPL",shell=True).strip().decode()
    sys.path.insert(1,py_matpl)
    py_scipy=subprocess.check_output("echo $PY_SCIPY",shell=True).strip().decode()
    sys.path.insert(1,py_scipy)
from string import *
import numpy as np
import re,os
#import matplotlib
#from plot_lib import *
from io_lib import *
from tables_lib import *
from fourier_lib import *
from string_lib import invert_selection
from copy import deepcopy
from particle_param import *


class ImpedanceWakeFileNotFoundError(IOError):
    """
    Raised when impedance/wake files cannot
    be found. It's typically the case when the computation
    has been launched but is not finished yet.
    """
    pass


class ImpedanceWakeFileUnreadableError(IOError):
    """
    Raised upon impedance/wake file reading, 
    when the file is empty or unfinished.
    """
    pass


class layer(object):
    '''class for definition of layer properties (for resistive-wall code)'''

    def __init__(self,rhoDC=0.5e-5,tau=4.2e-12,epsb=1,mur=1,fmu=np.inf,thickness=np.inf):

        # default values for an infinite CFC layer (carbon-reinforced carbon). SI units here.
        self.rhoDC=rhoDC;
        self.tau=tau;
        self.epsb=epsb;
        self.mur=mur;
        self.fmu=fmu;
        self.thickness=thickness;


def CFC_layer(thickness=25e-3,Ageing_factor=1):
    '''define a layer of CFC (carbon fiber-reinforced carbon), type Tatsuno AC150
    see N. Mounet PhD thesis for references. It consider an ageing factor to make wort the resisitivity with time.'''

    return layer(rhoDC=5e-6*float(Ageing_factor),tau=4.2e-12,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def Ti6Al4V_layer(thickness=54e-3):
    '''define a layer of Ti6Al4V, used for TDIS'''

    return layer(rhoDC=1.7e-6,tau=0,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def graphite_layer(thickness=25e-3):
    '''define a layer of graphite, type SGL R4550
    see N. Mounet PhD thesis for references'''

    return layer(rhoDC=15e-6,tau=1.3e-12,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def carbon_layer(rhoDC=10e-6,thickness=25e-3):
    '''define a layer of carbon material, with user-defined resistivity
    relaxation time is changed accordingly (using parameters of CFC or graphite), proportionally to conductivity
    see N. Mounet PhD thesis for references.'''

    return layer(rhoDC=rhoDC,tau=4.2e-12*5e-6/rhoDC,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def carbon_3D_layer(rhoDC=30e-6,thickness=25e-3):
    '''define a layer of 3D carbon material, with user-defined resistivity
    relaxation time is changed accordingly (using parameters of CFC or graphite), proportionally to conductivity
    see N. Mounet PhD thesis for references.'''

    return layer(rhoDC=rhoDC,tau=0,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def aC_layer(thickness=0.5e-6):
    '''define a layer of amorphous carbon material,
    relaxation time is changed accordingly (using parameters of CFC or graphite), proportionally to conductivity
    see N. Mounet PhD thesis for references.'''
    rhoDC=1e-2; # pessimistic value from M. Taborelli

    return layer(rhoDC=rhoDC,tau=4.2e-12*5e-6/rhoDC,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def aC_HL_layer(thickness=0.5e-6):
    '''define a layer of amorphous carbon material, as the one used fo HL-LHC 
    coating (for e-cloud mitigation).
    Relaxation time is changed accordingly (using parameters of CFC or graphite), 
    proportionally to conductivity (see N. Mounet PhD thesis for references.)'''
    rhoDC = 1/400.
    epsb = 5.4
    # both values above from N. Biancacci - see talk by B. Salvant, Na Wang and 
    # Carlo Zannini, "Update on new triplet beam screen impedance", on 14/12/2015

    return layer(rhoDC=rhoDC,tau=4.2e-12*5e-6/rhoDC,epsb=epsb,mur=1,fmu=np.inf,thickness=thickness);


def Cu300K_layer(thickness=2e-3):
    '''define a layer of pure copper material at room temperature
    see N. Mounet PhD thesis for references.'''

    return layer(rhoDC=1.7e-8,tau=0.027e-12,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def Cu300K_in_TDI_layer(thickness=1e-6):
    '''define a layer of copper material, at room temperature, on top of the TDI jaw (most
    probably contaminated by other species, so more resistive - see some mesurements by S. Calatroni & W. Vollenberg)
    Also used for CU coating on CFC in TCDQ - see measurements by W. Vollenberg
    in https://edms.cern.ch/document/1337884/1
    '''

    return layer(rhoDC=2.6e-8,tau=0.027e-12/(2.6/1.7),epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def CuDCMS_layer(thickness=1e-6):
    '''define a layer of copper coating material on top of graphite (DCMS
    - Direct Current Magnetron Sputtering, done by CERN & used e.g. in
     high-rad mat. expts. HRTM35), at room temperature
     (see F. Carra, C. Accettura, N. Biancacci)
    '''

    return layer(rhoDC=1./14e6,tau=0.027e-12/(58.8/14),epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def CuHIPIMS_layer(thickness=1e-6):
    '''define a layer of copper coating material on top of graphite
    (HIPIMS - High-Power Impulse Magnetron Sputtering, done by DTI),
    at room temperature (see F. Carra, C. Accettura, N. Biancacci)
    '''

    return layer(rhoDC=1./45e6,tau=0.027e-12/(58.8/45),epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def CuCr1Zr_layer(thickness=1e-6):
    '''define a layer of Coper Zirconium material for the third block of the TDI. Conductivity is 80%IACS. Inigo Lamas Garcia 02/11/2015)
    '''

    return layer(rhoDC=2.15e-8,tau=0,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def W_layer(thickness=25e-3):
    '''define a layer of pure tungsten material
    see N. Mounet PhD thesis for references.'''

    return layer(rhoDC=5.4e-8,tau=0.005e-12,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def Inermet180_layer(thickness=25e-3):
    '''define a layer of tungsten alloy as it is in many LHC collimators
    see Federico Carra <Federico.Carra@cern.ch> and 
    Jorge Guardia Valenzuela <jorge.guardia.valenzuela@cern.ch> (email on 27/08/2018 to N. Mounet)
    Relaxation time is the one of pure tungsten.'''

    return layer(rhoDC=1./10.5e6,tau=0.005e-12,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def CuCD_layer(thickness=25e-3):
    '''define a layer of copper-diamond alloy as it is in some LHC collimators
    see Federico Carra <Federico.Carra@cern.ch> and 
    Jorge Guardia Valenzuela <jorge.guardia.valenzuela@cern.ch> (email on 18/12/2019 to N. Mounet)
    Relaxation time is neglected.'''

    return layer(rhoDC=1./9.3e6,tau=0.,epsb=1,mur=1,fmu=np.inf,thickness=thickness)


def TZM_layer(thickness=25e-3):
    '''define a layer of TZM (Mo-based) alloy as it is in some LHC collimators
    see F.-X. Nuiry and F. Carra (collimator materials table, WP5.2 meeting on 13/04/2022)
    Relaxation time is neglected.'''

    return layer(rhoDC=1./14.5e6,tau=0.,epsb=1,mur=1,fmu=np.inf,thickness=thickness)


def Alu_layer(rhoDC=2.65e-8,thickness=1e-2):
    '''define a layer of pure aluminum material.'''

    return layer(rhoDC=2.65e-8,tau=0,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def Ti_layer(thickness=3e-6):
    '''define a layer of pure titanium material
    see N. Mounet PhD thesis for references.'''

    return layer(rhoDC=4.3e-7,tau=0,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def Ti_in_TDI_layer(thickness=3e-6):
    '''define a layer of titanium material as it is in the LHC TDI (coating on hBN)
    see N. Mounet Evian 2011 paper for references.'''

    return layer(rhoDC=2.5e-6,tau=0,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def Al_layer(thickness=54e-3):
    '''define a layer of pure aluminum material
    see N. Mounet PhD thesis for references.'''

    return layer(rhoDC=2.7e-8,tau=0.008e-12,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def TiB2_layer(thickness=5e-6):
    ''' from F.Carra email 27-04-2015 11.1 MS/m (Mo is 20MS/m)'''
    return layer(rhoDC=9e-8,tau=0,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def TiN_layer(thickness=5e-6):
    ''' from F.Carra email 27-04-2015 2.5 MS/m (Mo is 20MS/m)'''
    return layer(rhoDC=4e-7,tau=0,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def Mo_layer(thickness=50e-6):
    '''define a layer of pure molybdenum material
    see A. Bertarelli...'''

    return layer(rhoDC=5.35e-8,tau=0.01e-12,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def Mo3_layer(thickness=50e-6):
    '''define a layer of pure molybdenum material
    with a factor 3 on resistivity from a possible effect of irradiation.'''

    return layer(rhoDC=1./6e6,tau=0.01e-12,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def Mo_on_MoC_layer(thickness=50e-6):
    '''define a layer of pure molybdenum material, as measured on a layer of Mo-graphite.
    See Jorge Guardia Valenzuela, Impedance meeting 24/08/2018.
    Relaxation time is the one of pure Mo.'''

    return layer(rhoDC=6.5e-8,tau=0.01e-12,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def Mo_on_graphite_layer(thickness=50e-6):
    '''define a layer of pure molybdenum material, as measured on a layer of graphite.
    See Carlotta Accettura.
    Relaxation time is the one of pure Mo.'''

    return layer(rhoDC=1/8e6,tau=0.01e-12,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def MoC_layer(thickness=25e-3):
    '''define a layer of molybdenum-graphite material
    mail from Luca Gentilini on 2-12-2014:
    Ici les proprietes :
    NB: (L)=longitudinal, (T)=transversal
    Densite: 2.65 g/cm3
    Conductivite electrique: 0.81 MS/m (L), 0.4 MS/m (T)
    Coefficient dilatation thermique: 3.0e-6 K-1 (L), 20e-6 K-1 (T)
    Module de Young: 62.7 MPa (L), 6.9 MPa (T)
    Tensile strength: 73 MPa (L), 17.5 MPa (T)
    We assume only L conductivity for the moment.

    From: Alessandro Bertarelli
    Sent: 05 January 2015 15:43
    To: Benoit Salvant
    Subject: Measurement of Electrical conductivity
    Conductivity of MoGr is ~1MS/m


    '''

    return layer(rhoDC=1e-6,tau=0.,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def hBN_layer(thickness=54e-3):
    '''define a layer of pure hBN material
    see N. Mounet PhD thesis for references.'''

    return layer(rhoDC=4e12,tau=0,epsb=4,mur=1,fmu=np.inf,thickness=thickness);


def ss304L_layer(thickness=np.inf):
    '''define a layer of stainless steel material (304L) at room temperature
    see N. Mounet PhD thesis for references.'''

    return layer(rhoDC=7.2e-7,tau=0,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def Torlon_layer(thickness=1.5e-6):
    '''define a layer of Torlon. Info from Benoit on 02 May 2017'''

    return layer(rhoDC=1e13,tau=0,epsb=5,mur=1,fmu=np.inf,thickness=thickness);


def AlMg4_layer(thickness=1.5e-6):
    '''define a layer of AlMg4. Info from Benoit on 02 May 2017'''

    return layer(rhoDC=5.9e-8,tau=0,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def NEG_layer(thickness=1.5e-6):
    '''define a layer of NEG material at room temperature e.g. see
    http://accelconf.web.cern.ch/accelconf/IPAC2014/papers/wepme050.pdf
    (rhoDC=1e-6, er=1)
    Update following B.Salvant email on 02 May 2017
    '''

    return layer(rhoDC=2.5e-5,tau=0,epsb=10,mur=1,fmu=np.inf,thickness=thickness);


def NEG_in_TDI_layer(rhoDC=6.e-6,thickness=1.e-6):
    '''define a layer of NEG material at room temperature coated on top of a porous hBN layer
    (lot of outgassing -> degrade electric properties, but not known by how much - cf. Sergio Calatroni)'''

    return layer(rhoDC=rhoDC,tau=0,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def ssP506_20K_layer(thickness=np.inf):
    '''define a layer of stainless steel material (P506) at 20K
    see N. Mounet PhD thesis for references.'''

    return layer(rhoDC=6e-7,tau=0,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def ss316LN_layer(thickness=np.inf):
    '''define a layer of stainless steel material (316LN) as for the PS chamber (see also Phys. Rev. Accel. Beams 19, 041001 ).'''

    return layer(rhoDC=75e-8,tau=0,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def Inconel_X750_layer(thickness=np.inf):
    '''define a layer of Inconel X750 alloy as in the PS vacuum chamber (see also Phys. Rev. Accel. Beams 19, 041001 ).'''

    return layer(rhoDC=120e-8,tau=0,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def Cu20K_layer(thickness=50e-6,RRR=70,B=7000/(2803.9*0.299792458)):
    '''define a layer of pure copper material at 20K
    with magnetoresistance taken into account (RRR=residual resistivity ratio,
    B=magnetic field). Default B given here is for the LHC (2803.9 m radius of curvature) at 7TeV.
    see N. Mounet PhD thesis for references.
    WARNING: RRR not properly taken into account - it should also be used to
    compute the B=0 value of resistivity (here 2.4e-10)'''

    rhoDC=2.4e-10*(1.0048+0.0038*RRR*B) # from C. Rathjen, EDMS 329882

    # tauAC formula from Ascroft-Mermin (Z=1 for copper), Drude model (used also for other
    # materials defined above) with parameters from CRC - Handbook of Chem. and Phys.
    me=9.10938e-31;e=1.60218e-19; # electron mass and elementary charge
    rhom=9.06 # volumic mass in g/cm3
    A=63.546; # atomic mass in g/mol
    Z=1;# number of valence electrons
    n=6.022e23*Z*rhom*1e6/A;
    tauAC=me/(n*rhoDC*e**2); # relaxation time (s)

    return layer(rhoDC=rhoDC,tau=tauAC,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def Cu50K_layer(thickness=50e-6,RRR=70,B=7000/(2803.9*0.299792458)):
    '''define a layer of pure copper material at 50K
    with magnetoresistance taken into account (RRR=residual resistivity ratio,
    B=magnetic field). Default B given here is for the LHC (2803.9 m radius of curvature) at 7TeV.
    We use "Kohler law" (F. Caspers et al, LHC project report 307).
    WARNING: RRR not taken into account!
    WARNING: Kohler law seems to underestimate magnetoresistance of real beam screen.'''

    rhoDC_B0=7.5e-10; # resistivity for B=0 at 50 K
    rhoDC_B0_273K=1.5e-8; # resistivity for B=0 at 273 K

    rhoDC=rhoDC_B0*(1+ 10.**(-2.69)*(B*rhoDC_B0_273K/rhoDC_B0)**1.055);

    # tauAC formula from Ascroft-Mermin (Z=1 for copper), Drude model (used also for other
    # materials defined above) with parameters from CRC - Handbook of Chem. and Phys.
    me=9.10938e-31;e=1.60218e-19; # electron mass and elementary charge
    rhom=9.06 # volumic mass in g/cm3
    A=63.546; # atomic mass in g/mol
    Z=1;# number of valence electrons
    n=6.022e23*Z*rhom*1e6/A;
    tauAC=me/(n*rhoDC*e**2); # relaxation time (s)

    return layer(rhoDC=rhoDC,tau=tauAC,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def Cu100K_layer(thickness=50e-6,RRR=70,B=7000/(2803.9*0.299792458)):
    '''define a layer of pure copper material at 100K
    with magnetoresistance taken into account (RRR=residual resistivity ratio,
    B=magnetic field). Default B given here is for the LHC (2803.9 m radius of curvature) at 7TeV.
    We use "Kohler law" (F. Caspers et al, LHC project report 307).
    WARNING: RRR not taken into account!
    WARNING: Kohler law seems to underestimate magnetoresistance of real beam screen.'''

    rhoDC_B0=0.33e-8; # resistivity for B=0 at 100 K
    rhoDC_B0_273K=1.5e-8; # resistivity for B=0 at 273 K

    rhoDC=rhoDC_B0*(1+ 10.**(-2.69)*(B*rhoDC_B0_273K/rhoDC_B0)**1.055);

    # tauAC formula from Ascroft-Mermin (Z=1 for copper), Drude model (used also for other
    # materials defined above) with parameters from CRC - Handbook of Chem. and Phys.
    me=9.10938e-31;e=1.60218e-19; # electron mass and elementary charge
    rhom=9.06 # volumic mass in g/cm3
    A=63.546; # atomic mass in g/mol
    Z=1;# number of valence electrons
    n=6.022e23*Z*rhom*1e6/A;
    tauAC=me/(n*rhoDC*e**2); # relaxation time (s)

    return layer(rhoDC=rhoDC,tau=tauAC,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def Cu_layer(thickness=50e-6,T=20,RRR=70,B=7000/(2803.9*0.299792458)):
    '''
    Define a layer of pure copper material at any temperature, any B field
    and any RRR.
    
    This function is meant to generalize and replace the various, half-wrong,
    copper layer functions (in particular Cu20K_layer, Cu50K_layer and Cu100K_layer).
    
    We use a magnetoresistance law fitted from the UPPER curve of the plot in 
    NIST, "Properties of copper and copper alloys at cryogenic temperatures",
    by Simon, Crexler and Reed, 1992 (p. 8-27, Fig. 8-14).
    The upper curve was chosen, as it is in relative agreement with C. Rathjen
    measurements and actual LHC beam screens (CERN EDMS document Nr. 329882).
    
    The resistivity vs. temperature and RRR is found from Hust & Lankford, 
    "Thermal conductivity of aluminum, copper, iron, and tungsten fortemperatures
    from 1K to the melting point", National Bureau of Standards, 1984, 
    from Eqs. (1.2.3) to (1.2.6) p. 8, with parameters from p. 22, EXCEPT P4 
    for which we took the opposite value (to get the right behaviour at high 
    temperature - most probably there is a typo in the reference).
    
    The law vs. temperature was found in good agreement with the NIST reference
    above (p. 8-5, Fig. 8-2).
    
    :param thickness: material thickness in m
    :param T: temperature in K
    :param RRR: residual resistivity ratio, i.e. rho(273K)/rho(0K) at B=0
    :param B: magnetic field in T
    
    :return: a layer object
    '''

    # magnetoresistance law (Kohler-like) drho/rho = f(B*rho_273K(B=0)/rho_T(B=0))
    def Kohler(BS):
        p = [0.029497104404715, 0.905633738689341, -2.361415783729567]
        if BS==0.:
            return 0.
        else:
            return 10.**np.polyval(p,np.log10(BS))
    
    # rho at B=0 vs. temperature and RRR
    def rho_vs_T(T,RRR):
        rhoc = 0. # residual (deviation from the law - unknown, we neglect it here)
        rho273K = 15.5*1e-9 # resistivity at 273K, in Ohm.m
        P1 = 1.171e-17
        P2 = 4.49
        P3 = 3.841e10
        P4 = -1.14
        P5 = 50.
        P6 = 6.428
        P7 = 0.4531
        rho0 = rho273K/RRR
        rhoi = P1*T**P2/(1.+P1*P3*T**(P2+P4)*np.exp(-(P5/T)**P6)) + rhoc
        rhoi0 = P7*rhoi*rho0/(rhoi+rho0)
        return (rho0+rhoi+rhoi0)        
    
    rhoDC_B0 = rho_vs_T(T,RRR) # resistivity for B=0
    Sratio = rho_vs_T(273.,RRR)/rhoDC_B0
    rhoDC = float('{:.3g}'.format(rhoDC_B0*(1.+Kohler(B*Sratio)))) # we round it to 3 significant digits

    # tauAC formula from Ascroft-Mermin (Z=1 for copper), Drude model (used also for other
    # materials defined above) with parameters from CRC - Handbook of Chem. and Phys.
    me = 9.10938e-31 # electron mass
    e = 1.60218e-19 # electron elementary charge
    rhom = 9.0 # Cu volumic mass in g/cm3 (it is actually 9.02 at 4K, 8.93 at 273K,
    # see https://www.copper.org/resources/properties/cryogenic/ )
    A = 63.546 # Cu atomic mass in g/mol
    Z = 1 # number of valence electrons
    n = 6.022e23*Z*rhom*1e6/A
    tauAC = float('{:.3g}'.format(me/(n*rhoDC*e**2))) # relaxation time (s) (3 significant digits)

    return layer(rhoDC=rhoDC,tau=tauAC,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def vacuum_layer(thickness=np.inf):
    '''define a vacuum layer'''

    return layer(rhoDC=np.inf,tau=0,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def Z_layer(rhoDC=1.7e-8,thickness=0.001):
    '''define a generic conducting layer'''

    return layer(rhoDC=rhoDC,tau=0,epsb=1,mur=1,fmu=np.inf,thickness=thickness);


def find_nmaterial(param_filename):
    ''' Find the number of materials in a parameter file.
    This is the number of columns for which the header begins by "[mM]aterial" '''

    nmat=0;flag=True;
    while flag:
        if nmat==0: res=read_ncol_file_identify_header(param_filename,'[mM]aterial',dispflag=False);
        else: res=read_ncol_file_identify_header(param_filename,'[mM]aterial'+str(nmat+1),dispflag=False);
        flag=(len(res)!=0);
        nmat+=flag

    if (nmat==0): print("Pb in find_nmaterial: no material column found !");sys.exit();

    return nmat;


def read_materials(param_filename,ind):
    ''' Read materials columns in a parameter file.
    These are columns for which the header begins by "[mM]aterial"
    ind is an array to re-order the materials extracted. '''

    # find the number of materials in the file
    nmat=find_nmaterial(param_filename)
    # read materials and thicknesses
    material=[];thick=[];
    for n in range(1,nmat+1):
        if nmat==1: strn='';
        else: strn=str(n);
        material.append(read_ncol_file_identify_header(param_filename,'[mM]aterial'+strn));
        material[-1]=[material[-1][i] for i in ind];
        if n<nmat:
            thick.append(read_ncol_file_identify_header(param_filename,'[tT]hickness'+strn));
            thick[-1]=[thick[-1][i] for i in ind]

    return material,thick,nmat;


def construct_layers(layers,thickness=[],Bfield=0,Age=[],rhoDC=[]):
    ''' Construct list of layer ojects from a list of either layer object
    OR material names (if material names then the corresponding
    thicknesses should be given). List elements can also be 'None' or None
    (then no layer).
    Bfield should be given in case of Cu20K (magnetoresistant).
    Age indicates an ageing factor for the primary and secondaries in CFC.
    Rho indicated if specific conductivities have to be changed. '''

    if not Age: Age=np.ones(len(thickness));
    layers_iw=[];
    for ilay,lay in enumerate(layers):
        if isinstance(lay,layer): layers_iw.append(lay);
        elif lay.startswith('Cu20K'):
            # cold copper case: take into account magnetoresistance
            eval('layers_iw.append('+lay+'_layer(thickness=thickness[ilay],B=Bfield))');
        elif lay.startswith('CFC'):
            eval('layers_iw.append('+lay+'_layer(thickness=thickness[ilay],Ageing_factor=Age[ilay]))')
        elif lay.startswith('Z'):
            eval('layers_iw.append('+lay+'_layer(thickness=thickness[ilay],rhoDC=rhoDC[ilay]))');
        elif (lay.startswith('None'))or(lay==None):
            # do nothing
            pass;
        else:
            # assume lay is a string with the name of the material
            eval('layers_iw.append('+lay+'_layer(thickness=thickness[ilay]))');

    return layers_iw;


def construct_layers_generalized(layers,thickness=[],**params):
    ''' Construct list of layer ojects from a list of either layer object
    OR material names (if material names then the corresponding
    thicknesses should be given). List elements can also be 'None' or None
    (then no layer).
    Additional keyword parameters can be added in params (ignored
    for any given layer that raises an exception when we pass them
    to the corresponding function). '''

    layers_iw=[]
    
    for ilay,lay in enumerate(layers):
        if isinstance(lay,layer):
            layers_iw.append(lay)
        elif (lay.startswith('None')) or (lay is None):
            # do nothing
            pass
        else:
            try:
                eval('layers_iw.append('+lay+'_layer(thickness=thickness[ilay],**params))')
            except TypeError:
                try:
                    # try to assume that each keyword param is actually a list
                    # to be associated to the list of different layers
                    the_params = {key: value[ilay] for key,value in list(params.items())}
                    eval('layers_iw.append('+lay+'_layer(thickness=thickness[ilay],**the_params))')
                except (TypeError,IndexError):
                    # if all the previous has failed, we just ignore the additional parameters
                    eval('layers_iw.append('+lay+'_layer(thickness=thickness[ilay]))')
    
    print("Layers: {}".format([(l.rhoDC,l.thickness) for l in layers_iw]))
    
    return layers_iw


class freq_param(object):
    '''class for definition of frequency sampling (all in Hz)'''

    def __init__(self,fmin=10,fmax=1e13,ftypescan=2,fsamplin=1e8,nflog=20,
        fminrefine=1e11,fmaxrefine=5e12,nrefine=5000,fadded=[1e-2,0.1,1,1e15]):

        self.fmin=fmin;
        self.fmax=fmax;
        self.ftypescan=ftypescan;
        self.fsamplin=fsamplin;
        self.nflog=nflog;
        self.fminrefine=fminrefine;
        self.fmaxrefine=fmaxrefine;
        self.nrefine=nrefine;
        self.fadded=fadded;


class z_param(object):
    '''class for definition of z (distance) sampling (all in m)'''

    def __init__(self,zmin=0.01,zmax=1e6,ztypescan=2,nzlog=100,
        zminrefine=0.5e-5,zmaxrefine=0.01,zsamplin=0.5e-5,zadded=[]):

        self.zmin=zmin;
        self.zmax=zmax;
        self.ztypescan=ztypescan;
        self.nzlog=nzlog;
        self.zminrefine=zminrefine;
        self.zmaxrefine=zmaxrefine;
        self.zsamplin=zsamplin;
        self.zadded=zadded;


class impedance_wake_input(object):
    '''class for impedance/wake analytical computations with ImpedanceWake2D'''

    def __init__(self,machine='LHC',gamma=7460.52,length=1,b=[2e-3],
        layers=[layer()],fpar=freq_param(),zpar=z_param(),
        longfactor=100.,waketol=1.e12,freqlinbisect=1.e11,
        geometry='round',Yokoya=[1,1,1,0,0],comment=''):

        '''all in SI units here.

        geometry can be 'round' or 'flat'

        Yokoya factors are used only for round geometry, to transform it into another
        geometry (flat, elliptical or rectangular). Order: long., xdip, ydip, xquad, yquad.

        b is the radius (for round) or half-gap (for flat).
        b can have one or two values: if two values and flat geometry, then top-bottom
        symmetry is 'no', and the first value is for the first upper layers, the second one
        for the first lower layer. Both values can still be equal.
        If one of the values in b is zero, the corresponding part (upper or lower) is empty
        (i.e. single plate case, like e.g. the LHC TCDQ).

        layers contain the layer properties (see above class 'layer'), first the
        upper layers then the lower ones (in the flat case). The last upper layer
        must be characterized   by an infinite thickness.

        longfactor: longitudinal factor, waketol: tolerance for wake calculation,
        freqlinbisect: frequency above which the mesh bisecting is linear.'''

        self.machine=machine;
        self.gamma=gamma;
        self.length=length;
        self.b=b;
        self.layers=layers;
        self.fpar=fpar;
        self.zpar=zpar;
        self.longfactor=longfactor;
        self.waketol=waketol;
        self.freqlinbisect=freqlinbisect;
        self.geometry=geometry;
        self.Yokoya=Yokoya;
        self.comment=comment;


def write_line(file1,sentence,value):
    '''write in file 'file1' (already opened) a line with the string 'sentence'
    followed by a tab and the number 'value'.
    if value=inf, 'Infinity' is printed.'''

    if (value==np.inf): val='Infinity'
    elif (np.isscalar(value)): val=str(value);
    else:
        val=str(np.array(value)).rstrip(']').lstrip('[');
        while (val.count('  ')>0): val=val.replace('  ',' ');
        if val.startswith(' '): val=val.lstrip(' ');

    file1.write(sentence+"\t"+val+"\n");


def write_layer(file1,n,layer):
    '''write the layer properties on the file (already opened) 'file1'
    layer has number 'n' (can be negative)'''

    write_line(file1,"Layer "+str(n)+" DC resistivity (Ohm.m):",layer.rhoDC);
    write_line(file1,"Layer "+str(n)+" relaxation time for resistivity (ps):",layer.tau*1e12);
    write_line(file1,"Layer "+str(n)+" real part of dielectric constant:",layer.epsb);
    write_line(file1,"Layer "+str(n)+" magnetic susceptibility:",layer.mur-1.);
    write_line(file1,"Layer "+str(n)+" relaxation frequency of permeability (MHz):",layer.fmu/1e6);
    write_line(file1,"Layer "+str(n)+" thickness in mm:",float(layer.thickness)*1e3);


def write_freq_param(file1,fpar):
    '''write the frequency sampling parameters on file 'file1' (already opened)'''

    write_line(file1,"start frequency exponent (10^) in Hz:",np.log10(fpar.fmin));
    write_line(file1,"stop frequency exponent (10^) in Hz:",np.log10(fpar.fmax));
    write_line(file1,"linear (1) or logarithmic (0) or both (2) frequency scan:",fpar.ftypescan);
    write_line(file1,"sampling frequency exponent (10^) in Hz (for linear):",np.log10(fpar.fsamplin));
    write_line(file1,"Number of points per decade (for log):",fpar.nflog);
    write_line(file1,"when both, fmin of the refinement (in THz):",fpar.fminrefine/1e12);
    write_line(file1,"when both, fmax of the refinement (in THz):",fpar.fmaxrefine/1e12);
    write_line(file1,"when both, number of points in the refinement:",fpar.nrefine);
    write_line(file1,"added frequencies (Hz):",fpar.fadded);


def write_z_param(file1,zpar):
    '''write the z sampling parameters on file 'file1' (already opened)'''

    write_line(file1,"linear (1) or logarithmic (0) or both (2) scan in z for the wake:",zpar.ztypescan);
    write_line(file1,"sampling distance in m for the linear sampling:",zpar.zsamplin);
    write_line(file1,"zmin in m of the linear sampling:",zpar.zminrefine);
    write_line(file1,"zmax in m of the linear sampling:",zpar.zmaxrefine);
    write_line(file1,"Number of points per decade for the logarithmic sampling:",zpar.nzlog);
    write_line(file1,"exponent (10^) of zmin (in m) of the logarithmic sampling:",np.log10(zpar.zmin));
    write_line(file1,"exponent (10^) of zmax (in m) of the logarithmic sampling:",np.log10(zpar.zmax));
    write_line(file1,"added z [m]:",zpar.zadded);


def write_impedance_wake_input(filename,iw_input):

    '''write an input file for ImpedanceWake2D (impedance or wake calc.)
    note that some lines might be useless (e.g. wake parameters for an impedance compuation)
    but this does not matter (those lines won't be read by ImpedanceWake2D).
    "iw_input" is an object of the class "impedance_wake_input".
    writes into file of name "filename" .'''

    file1=open(filename,'w');

    topbot='yes';flagtopbot=True; # default case is with top-bottom symmetry
    nup=len(iw_input.layers);nlow=0;

    write_line(file1,"Machine:",iw_input.machine);
    write_line(file1,"Relativistic Gamma:",iw_input.gamma);
    write_line(file1,"Impedance Length in m:",iw_input.length);

    if (iw_input.geometry.startswith('round')):
        # case of round geometry
        write_line(file1,"Number of layers:",len(iw_input.layers));
        write_line(file1,"Layer 1 inner radius in mm:",iw_input.b[0]*1e3);

    elif (iw_input.geometry.startswith('flat')):
        # case of flat geometry

        if (len(iw_input.b)>=2):
            # without top-bottom symmetry
            topbot='no';flagtopbot=False;
            # number of upper and lower layers
            if (iw_input.b[0]==0):
                nup=0;
                nlow=len(iw_input.layers);
            elif (iw_input.b[1]==0):
                nup=len(iw_input.layers);
                nlow=0;
            else:
                # find number of upper layers (last upper layer has infinite thickness)
                nup=0;
                while (getattr(iw_input.layers[nup],'thickness')!=np.inf): nup+=1;
                nup+=1;
                if (nup==len(iw_input.layers)):
                    print("Pb (flat case): both half-gaps are non zero and yet upper or lower part is empty!");sys.exit();
                # the rest are lower layers
                nlow=len(iw_input.layers)-nup;

        write_line(file1,"Number of upper layers in the chamber wall:",nup);
        if (nup>0): write_line(file1,"Layer 1 inner half gap in mm:",iw_input.b[0]*1e3);

    # write upper (or round) layers
    for n in range(nup): write_layer(file1,n+1,iw_input.layers[n]);

    if iw_input.geometry.startswith('flat'):
        # for flat geometry
        if not(flagtopbot):
            # case without top-bottom symmetry
            write_line(file1,"Number of lower layers in the chamber wall:",nlow);
            if (nlow>0): write_line(file1,"Layer -1 inner half gap in mm:",iw_input.b[1]*1e3);
            # write lower layers (if any)
            for n in range(nup,nup+nlow): write_layer(file1,nup-n-1,iw_input.layers[n]);

        write_line(file1,"Top bottom symmetry (yes or no):",topbot);

    # write frequency parameters
    write_freq_param(file1,iw_input.fpar)

    # write z (distance) parameters
    write_z_param(file1,iw_input.zpar)

    if (iw_input.geometry.startswith('round')):
        # for round geometry
        write_line(file1,"Yokoya factors long, xdip, ydip, xquad, yquad:",iw_input.Yokoya);

    write_line(file1,"factor weighting the longitudinal impedance error:",iw_input.longfactor);
    write_line(file1,"tolerance (in wake units) to achieve:",iw_input.waketol);
    write_line(file1,"frequency above which the mesh bisecting is linear [Hz]:",iw_input.freqlinbisect);
    write_line(file1,"Comments for the output files names:",iw_input.comment);

    file1.close();


def kicker_imp(a,b,d,L,material,gamma,fpar, er_input = '1', ur_input='1', sig_input='1e6'):

    from scipy.constants import c, e, epsilon_0
    import cmath
    imp_mod=[];

    freq=freqscan_from_fpar(fpar)
    Z0=120*np.pi;
    q=e;
    u0=Z0/c;
    e0=epsilon_0
    beta=np.sqrt(1.-1./gamma**2)
    v=beta*c;

    if 'ferrite' in material:
        rho=1.e6;
        ui=460.;
        kf=1./(20e6);
        media_ur='1+(ui/(1+(kf*kf*f*f)))-1j*((ui*kf*f)/(1+(kf*kf*f*f)))'
        media_er='12-1j*(1/(om*rho*e0))';
    elif 'graphite' in material:
        media_ur='1;';
        media_er='1-1j*1e5/(om*e0)';
    else:
        media_ur='%s'%ur_input;
        media_er='%s-1j*%s/(om*e0)'%(er_input,sig_input);

    print(media_ur, media_er)
    # longitudinal and quadrupolar impedance
    Nmax=150

    N=Z0*q/2/beta/a;
    An=np.zeros(Nmax);
    Zlong=[]; Zxquad=[]; Zyquad=[]

    for jj, f in enumerate(freq):

        n=np.arange(Nmax)
        g1n=(2*n+1)*np.pi/2/a;
        k1n=(2*n+1)*np.pi/2/a;

        om=2*np.pi*f;
        k0=om/c;
        k=om/beta/c;
        kr=om/(beta*c*gamma);

        ur=eval(media_ur);
        er=eval(media_er);

        k2n=np.array(list(map(cmath.sqrt,(-kr**2-k1n**2))))
        g2n=np.array(list(map(cmath.sqrt,(k0**2*ur*er-k**2-g1n**2))))
        x1=k2n*b
        x2=g2n*(b-d)
        Eyn=N*np.exp(1j*x1)

        An = Eyn*( k2n*np.cos(x1)*(1j*k2n*gamma**2*(k**2-k0**2*er*ur)-g2n*k**2*er*1./np.tan(x2)) + np.sin(x1)*\
        (k**2*(kr**2*er*ur+k1n**2*(1+gamma**2*(-1+er*ur))) +1j*g2n*k2n*kr**2*gamma**2*ur*np.tan(x2)) )/\
        (k*k2n*gamma**2*(-g2n*k2n*er*np.cos(x1)**2*1./np.tan(x2)+(k**2+(-k0**2+kr**2)*er*ur+k1n**2*(1+er*ur))*\
        np.cos(x1)*np.sin(x1)-g2n*k2n*ur*np.sin(x1)**2*np.tan(x2)));


        Zlong.append(-L/q*np.sum(An))
        Zxquad.append(L/k/q*np.sum(An*k1n**2))
        Zyquad.append(L/k/q*np.sum(An*k2n**2))

        # Dipolar Impedance Y

    YAn=np.zeros(Nmax);
    Zydip=[]

    for jj, f in enumerate(freq):

        n=np.arange(Nmax)

        k1n=(2*n+1)*np.pi/2/a;

        om=2*np.pi*f;
        k0=om/c;
        k=om/beta/c;
        kr=om/(beta*c*gamma);

        eval(media_ur);
        eval(media_er);

        k2n=np.array(list(map(cmath.sqrt,(-kr**2-k1n**2))));
        k3n=np.array(list(map(cmath.sqrt,(k0**2*ur*er-k**2-k1n**2))));
        x1=k2n*b;
        x2=k3n*(b-d);


        M1n=1j*k2n**2*k0**2*(1-er*ur)+k2n*kr**2*(1j*k2n-k3n*er*1./np.tan(x2))
        M2n=-k1n**2*k0**2*(1-er*ur)+k2n*kr**2*ur*(-k2n*er+1j*k3n*np.tan(x2));
        M3n=k2n*k3n*(er*np.sin(x1)**2*1./np.tan(x2)+ur*np.cos(x1)**2*np.tan(x2));
        M4n=(k0**2*(-1+er*ur)+k2n**2*(1+er*ur))*np.cos(x1)*np.sin(x1);

        YAn=(1j*Z0*k2n*L/2/a/beta/k)  *  k2n*np.exp(1j*k2n*b)*(-M1n*np.sin(x1)+M2n*np.cos(x1))/(k*k2n*(M3n-M4n));

        Zydip.append(np.sum(YAn))

        # Dipolar Impedance X

    XAn=np.zeros(Nmax);
    Zxdip=[]


    for jj, f in enumerate(freq):

        n=np.arange(Nmax)
        k1n=(n)*np.pi/a;

        om=2*np.pi*f;
        k0=om/c;
        k=om/beta/c;
        kr=om/(beta*c*gamma);

        eval(media_ur);
        eval(media_er);

        k2n=np.array(list(map(cmath.sqrt,(-kr**2-k1n**2))));
        k3n=np.array(list(map(cmath.sqrt,(k0**2*ur*er-k**2-k1n**2))));

        x1=k2n*b;
        x2=k3n*(b-d);


        M1n=1j*k2n**2*k0**2*(1-er*ur)+k2n*kr**2*(1j*k2n-k3n*er*1./np.tan(x2))
        M2n=-k1n**2*k0**2*(1-er*ur)+k2n*kr**2*ur*(-k2n*er+1j*k3n*np.tan(x2));
        M3n=-k2n*k3n*(er*np.cos(x1)**2*1./np.tan(x2)+ur*np.sin(x1)**2*np.tan(x2));
        M4n=(k0**2*(-1+er*ur)+k2n**2*(1+er*ur))*np.cos(x1)*np.sin(x1);


        XAn=(1j*Z0*L/2/a/beta/k)  *k1n*np.exp(1j*k2n*b)*(M1n*np.cos(x1)+M2n*np.sin(x1))/(k*(M3n-M4n));
        #XAn=-L*(-1i*Z0*k2n)./(2*a*beta)  .*k1n.*  (-1i*k1n.^2*(1-er*ur))./(k2n.*(M3n-M4n));
        Zxdip.append(np.sum(XAn))

    imp_mod.append(impedance_wake(a=0,b=0,c=0,d=0, plane ='z', var=freq, func=np.column_stack((np.real(Zlong),np.imag(Zlong))) ))
    imp_mod.append(impedance_wake(a=1,b=0,c=0,d=0, plane ='x', var=freq, func=np.column_stack((np.real(Zxdip),np.imag(Zxdip))) ))
    imp_mod.append(impedance_wake(a=0,b=1,c=0,d=0, plane ='y', var=freq, func=np.column_stack((np.real(Zydip),np.imag(Zydip))) ))
    imp_mod.append(impedance_wake(a=0,b=0,c=1,d=0, plane ='x', var=freq, func=np.column_stack((np.real(Zxquad),np.imag(Zxquad))) ))
    imp_mod.append(impedance_wake(a=0,b=0,c=0,d=1, plane ='y', var=freq, func=np.column_stack((np.real(Zyquad),np.imag(Zyquad))) ))

    return imp_mod


def transverse_imp_taper_round_Yokoya(a,b,tantheta):

    '''computes transverse R/Q (broad-band) in Ohm/m of one single round linear taper using
    Yokoya's formula (CERN SL/90-88)
    parameters: small radius a, large one b, slope of taper is tantheta=tan(theta)

    valid under the conditions (sigma_z = bunch length):
         b tan(theta)/sigma_z << 1, and either a/sigma_z >>1 or (b-a)/a<<1
    This is equivalent (with k=omega/c ~ 1/sigma_z at the main bunch spectrum frequency) to
         kb tan(theta) << 1, and either ka >> 1 or (b-a)/a<<1'''

    I=tantheta*(1./a-1./b);
    Z0=4e-7*np.pi*299792458; # free space impedance
    # shunt impedance /Q
    RoverQ=I*Z0/(2.*np.pi);

    return RoverQ;


def longitudinal_imp_taper_round_Yokoya(a,b,tantheta,fcutoff):

    '''computes transverse R/Q (broad-band) in Ohm of one single round linear taper using
    Yokoya's formula (CERN SL/90-88)
    parameters: small radius a, large one b, slope of taper is tantheta=tan(theta),
    beam pipe cutoff frequency fcutoff and relativistic velocity factor beta.

    valid under the conditions (sigma_z = bunch length):
        b tan(theta)/sigma_z << 1, and either a/sigma_z >>1 or (b-a)/a<<1
    This is equivalent (with k=omega/c ~ 1/sigma_z at the main bunch spectrum frequency) to
        kb tan(theta) << 1, and either ka >> 1 or (b-a)/a<<1'''

    I=tantheta*(b-a);
    mu0=4e-7*np.pi; # free space inductance
    # shunt impedance /Q
    RoverQ=I*mu0*fcutoff/2.;

    return RoverQ;


def broadband_imp_taper_flat_Stupakov(a,b,tantheta,w,fcutoff,approx=True,
        listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad']):

    '''computes R/Q (broad-band) in Ohm(/m if not long.) of one single rectangular
    linear taper using Stupakov's formulas (Phys. Rev. STAB 10, 094401 - 2007),
    multiplied by Z0*c/(4*pi) to convert to SI units.
    Taper is in the vertical plane.
    parameters: small vertical half-gap a, large one b, slope of taper is tantheta=tan(theta),
    half-width=w (constant over taper), fcutoff the cutoff frequency (used only
    for the longitudinal component), approx=True if one assumes small b/w,
    listcomp is a list with the names of the components for which ones computes the R/Q
    (ex: Zlong, Zydip, Zxquad, etc.), beta the relativistic velocity factor

    WARNING: we use here half apertures (half-gap and half-width) whereas
    Stupakov's paper is expressed with full apertures. This does not make
    any difference except for an additional factor 4 here for longitudinal
    impedance.

    returns a list of R/Q for each component

    valid under the conditions of low frequency and length of taper much larger than
    its transverse dimensions'''

    from scipy import special as sp

    mu0=4e-7*np.pi;
    Z0=mu0*299792458; # free space impedance
    RoverQ=[];

    for comp in listcomp:

        if comp.endswith('long'):
            cst=4.*mu0*fcutoff/2.; # factor 4 due to use of half-gaps here
            Gint=0; # G0=F used in integral
            gpow=0; # power of 1/g in integral
            I=7.*sp.zeta(3,1)/(2.*np.pi**2)*tantheta*(b-a); #  approx. integral (sp.zeta(3.,1.) is Riemann zeta function at x=3)
            #I=7.*1.202057/(2.*np.pi**2)*tantheta*(b-a);

        elif comp.endswith('ydip'):
            cst=Z0*w*np.pi/4.;
            Gint=1; # G1 used in integral
            gpow=3; # power of 1/g in integral
            I=tantheta*(1./(a**2)-1./(b**2))/(2.*np.pi); # approx. integral

        elif comp.endswith('xquad'):
            cst=-Z0*np.pi/4.;
            Gint=2; # G2 used in integral
            gpow=2; # power of 1/g in integral
            I=tantheta*(1./a-1./b)/(np.pi**2); # approx. integral

        elif comp.endswith('xdip'):
            cst=Z0*np.pi/4.;
            Gint=3; # G3 used in integral
            gpow=2; # power of 1/g in integral
            I=tantheta*(1./a-1./b)/(np.pi**2); # approx. integral

        elif comp.endswith('yquad'):
            cst=Z0*np.pi/4.;
            Gint=2; # G2 used in integral
            gpow=2; # power of 1/g in integral
            I=tantheta*(1./a-1./b)/(np.pi**2); # approx. integral

        if not(approx):
            # computes numerically the integral instead of using its approximation
            from scipy import integrate as inte
            I,err=inte.quadrature(integrand_Stupakov,a,b,args=(w,Gint,gpow),tol=1.e-3,maxiter=200);
            I *= tantheta # put back g' factor that was dropped

        # shunt impedance /Q
        RoverQ.append(cst*I);

    return np.array(RoverQ);


def integrand_Stupakov(g,w,Gint,gpow):
    '''computes integrand for Stupakov integral for a rectangular linear taper at a given
    vertical half-gap(s) g (can be an array)
    Note that (g')^2 has been taken out of the integral.

    w is the width of the taper.
    Gint is the indice of the G function to use (0 is for G0=F in Stupakov's paper)
    gpow is the power to which 1/g is taken
    See Phys. Rev. STAB 10, 094401 (2007)'''

    return G_Stupakov(g/w,Gint)/(g**gpow);


def G_Stupakov(x,Gint):
    '''computes F (when Gint=0) or G_[Gint] function from Stupakov's formulas
    for a rectangular linear taper at a given x=g/w ratio
    x can be an array
    See Phys. Rev. STAB 10, 094401 (2007)'''

    x=np.array(create_list(x));
    G=np.zeros(len(x));

    for ix,xelem in enumerate(x):

        oldG=1.;
        eps=1e-5; # relative precision of the summation
        incr=10; # increment for sum computation
        Garray=np.zeros(incr);
        m=np.arange(incr);
        mlimit=1e6;

        while (abs(G[ix]-oldG)>eps*abs(G[ix]))and(m[-1]<mlimit):
            Garray=G_element(m,xelem,Gint);
            oldG=G[ix];
            G[ix]+=np.sum(Garray);
            m+=incr;

        #print Gint,m[-1]
        if (m[-1]>=mlimit): print("Pb in G_Stupakov: maximum number of elements reached !",m[-1],xelem,", err=",abs((G[ix]-oldG)/G[ix]))

    return G;


def G_element(m,x,Gint):
    '''computes each element of series defining F (when Gint=0) or G_[Gint] function
    from Stupakov's formulas for a rectangular linear taper at a given x=g/w ratio
    m is an array of integers where we compute those monomials
    See Phys. Rev. STAB 10, 094401 (2007)'''

    if Gint==0:
        val = (2*m+1)*np.pi*x/2
        res2 = ((1 - np.exp(-2*val))/(1 + np.exp(-2*val)) * (2*(np.exp(-val))/(1 + np.exp(-2*val)))**2    )/(2*m+1)

        return res2;

    elif Gint==1:
        val = (2*m+1)*np.pi*x/2
        res2 = x**3*(2*m+1) * (4*np.exp(-2*val)) * (1+np.exp(-2*val)) /  ( (1-np.exp(-2*val))**3)
        return res2;

    elif Gint==2:
        val = (2*m+1)*np.pi*x/2
        res2 = x**2*(2*m+1)*(((1 - np.exp(-2*val))/(1 + np.exp(-2*val)) * (2*(np.exp(-val))/(1 + np.exp(-2*val)))**2    ))
        return res2;

    elif Gint==3:
        val = m*np.pi*x
        res2 = x**2*(2*m)*((1 - np.exp(-2*val))/(1 + np.exp(-2*val)) * (2*(np.exp(-val))/(1 + np.exp(-2*val)))**2    )
        return res2

    else:
        print("Pb in G_element for Stupakov's formulas: Gint not 0, 1, 2 or 3 !")


def transverse_broadband_holes_in_round_pipe_Kurennoy(Lh,Wh,b,eta):

    ''' compute transverse (broad-band) R/Q of holes in a round pipe,
    per unit length in longitudinal (so in Ohm/m^2).
    From Kurennoy's formula in Part. Acc. 1995, vol. 50, p 167-175
    Note: it is different from Gluckstern's formula (CERN SL/92-06 or
    Phys. ReV. A 46 (2), p. 1110-1115, 1992) by a factor 2.

    Thanks to A. Mostacci, La Sapienza, Rome

    Lh: hole length, Wh: hole width,
    b: pipe radius, eta: proportion of total surface covered by holes.'''

    Z0=4e-7*np.pi*299792458; # free space impedance
    A=Wh*(Lh-Wh)+np.pi*Wh**2/4; # hole area

    # polarizabilities for rounded rectangular holes (A. Mostacci)
    # NOTE (from A. Mostacci): thickness SHOULD NOT be taken into account for
    # transverse impedance; elec. and mag. moments of holes on the INSIDE of the pipe
    alphae,alpham=polarizabilities_hole_Mostacci(Lh,Wh,0,Cm=1.,Ce=1.);

    # shunt impedance /Q /L
    RoverQ=2.*Z0*eta*(alphae+alpham)/(np.pi*A*b**3);

    return RoverQ;


def longitudinal_broadband_holes_in_round_pipe_Kurennoy(Lh,Wh,b,eta,fr=5e9):

    '''compute longitudinal (broad-band) R/Q of pumping holes in a round pipe,
    per unit length in longitudinal (so in Ohm/m).
    From Kurennoy's formula in Part. Acc. 1995, vol. 50, p 167-175

    Thanks to A. Mostacci, La Sapienza, Rome

    Lh: hole length, Wh: hole width,
    b: pipe radius, eta: proportion of total surface covered by holes,
    fr: cutoff frequency of broad-band resonator.'''

    omegar=2.*np.pi*fr;
    mu0=4e-7*np.pi; # free space permeability
    A=Wh*(Lh-Wh)+np.pi*Wh**2/4; # hole area

    # polarizabilities for rounded rectangular holes (A. Mostacci)
    # NOTE (from A. Mostacci): thickness SHOULD NOT be taken into account for
    # broad-band longitudinal impedance; elec. and mag. moments of holes on the INSIDE of the pipe
    alphae,alpham=polarizabilities_hole_Mostacci(Lh,Wh,0,Cm=1.,Ce=1.);

    # shunt impedance /Q /L
    RoverQ=omegar*mu0*eta*(alphae+alpham)/(2.*np.pi*A*b);

    return RoverQ;


def longitudinal_imp_holes_in_round_pipe_Mostacci(freq,Lh,Wh,T,b,d,eta,rhob,rhod,length,fcutoff=5e9,
        nb_holes_per_crosssection=8,Cm=1.,Ce=1.):

    '''computes the total (real + imag.) longitudinal impedance for holes
    in a round beam pipe.
    From A. Mostacci PhD thesis: take approximate formula for real part (eq. 1.35)
    and first term of eq. 1.16 for imaginary part (same as Kurennoy's inductive impedance,
    see Part. Acc. 1995, vol. 50, p 167-175).

    NOTE: real part of impedance not linear with length (due to interference between holes -> goes
    as N^2 with N number of holes, up to the attenuation length - then it is linear)

    freq: frequencies at which we compute impedance, Lh: hole length, Wh: hole width,
    T: pipe wall thickness, b: pipe inner radius,
    d: radius of outer shell (behind the pipe of radius b and thickness T),
    eta: proportion of total surface covered by holes, rhob: conductivity of pipe wall,
    rhod: conductivity of outer shell, length: length of the device,
    nb_holes_per_crosssection: number of holes per cross-section of the device (default=value for LHC beam screen),
    fcutoff: cutoff frequency (above-> zero impedance) (default = LHC one, for r=18.4mm aperture).
    Cm and Ce are constants found from simulations, here taken as 1 (pessimistic).'''

    omega=2.*np.pi*freq;

    clight=299792458; # speed of light
    mu0=4e-7*np.pi; # free space permeability
    Z0=mu0*clight; # free space impedance
    Area=Wh*(Lh-Wh)+np.pi*Wh**2/4; # hole area
    b2=b+T; # pipe outer radius

    a=(np.sqrt(rhob)/b2+np.sqrt(rhod)/d)*np.sqrt(mu0/2.)/(2.*Z0*np.log(d/b2)); # attenuation constant over sqrt(omega)
    alpha=a*np.sqrt(omega); # attenuation constant

    N=eta*2*np.pi*b2*length/Area; # number of holes
    D=nb_holes_per_crosssection*Area/(2.*np.pi*b2*eta); # longitudinal spacing between holes

    # polarizabilities for rounded rectangular holes (Kurennoy, also in A. Mostacci) on the INSIDE
    # of the pipe (used for constant inductive impedance only)
    alphae0,alpham0=polarizabilities_hole_Mostacci(Lh,Wh,0,Cm=1.,Ce=1.);

    # polarizabilities for rounded rectangular holes (A. Mostacci) on the OUTSIDE
    # of the pipe (used for the real part of the impedance); take into account finite
    # pipe thickness T.
    alphae,alpham=polarizabilities_hole_Mostacci(Lh,Wh,T,Cm=Cm,Ce=Ce);

    Zl=np.zeros((len(freq),2));
    ind=np.where(freq<=fcutoff);
    # real part of longitudinal impedance, from A. Mostacci PhD, eq. 1.35 (approx. formula)
    # (due to propagating TEM coaxial mode outside the pipe)
    Zl[ind,0]=Z0*omega[ind]**2*(alpham+alphae)**2*(N/2.+N/(alpha[ind]*D)+(np.exp(-alpha[ind]*D*N)-1)/(alpha[ind]*D)**2)/(16*np.pi**3*b**2*b2**2*np.log(d/b2)*clight**2);

    # add also broad-band model, from Kurennoy (~constant imag(Z/n) )
    #Zl[ind,1]=mu0*omega[ind]*(alpham0+alphae0)*N/(4.*np.pi**2*b**2);
    # computes shunt impedance
    R=length*longitudinal_broadband_holes_in_round_pipe_Kurennoy(Lh,Wh,b,eta,fr=fcutoff)
    Zl += resonator_long_impedance(R,fcutoff,1.,freq);

    return Zl;


def transverse_broadband_single_hole_in_round_pipe_Kurennoy(Lh,Wh,b):

    '''compute transverse (broad-band) R/Q of a single hole in a round pipe (in Ohm/m).
    From Kurennoy's formula in Part. Acc. 1995, vol. 50, p 167-175
    We assume the hole is right in front of the beam azimuthally (cosine factor in Kurennoy's formula =1)
    Note: it is different from Gluckstern's formula (CERN SL/92-06 or
    Phys. ReV. A 46 (2), p. 1110-1115, 1992) by a factor 2.

    Thanks to A. Mostacci, La Sapienza, Rome

    Lh: hole length, Wh: hole width, b: pipe radius.'''

    Z0=4e-7*np.pi*299792458; # free space impedance

    # polarizabilities for rounded rectangular hole (A. Mostacci)
    # NOTE (from A. Mostacci): thickness SHOULD NOT be taken into account for
    # transverse impedance; elec. and mag. moments of hole on the INSIDE of the pipe
    alphae,alpham=polarizabilities_hole_Mostacci(Lh,Wh,0,Cm=1.,Ce=1.);

    # shunt impedance /Q
    RoverQ=Z0*(alphae+alpham)/(np.pi**2*b**4);

    return RoverQ;


def longitudinal_broadband_single_hole_in_round_pipe_Kurennoy(Lh,Wh,b,fr=5e9):

    '''compute longitudinal (broad-band) R/Q of a single hole in a round pipe (in Ohm),
    From Kurennoy's formula in Part. Acc. 1995, vol. 50, p 167-175

    Thanks to A. Mostacci, La Sapienza, Rome

    Lh: hole length, Wh: hole width, b: pipe radius,
    fr: cutoff frequency of broad-band resonator.'''

    omegar=2.*np.pi*fr;
    mu0=4e-7*np.pi; # free space permeability

    # polarizabilities for rounded rectangular hole (A. Mostacci)
    # NOTE (from A. Mostacci): thickness SHOULD NOT be taken into account for
    # broad-band longitudinal impedance; elec. and mag. moments of hole on the INSIDE of the pipe
    alphae,alpham=polarizabilities_hole_Mostacci(Lh,Wh,0,Cm=1.,Ce=1.);

    # shunt impedance /Q
    RoverQ=omegar*mu0*(alphae+alpham)/(4.*np.pi**2*b**2);

    return RoverQ;


def polarizabilities_hole_Mostacci(Lh,Wh,T,Cm=1.,Ce=1.):

    '''compute polarizabilities (electric and magnetic) of small (w.r.t. wavelength) rounded
    rectangular holes (i.e. rectangle with semi-circles at the ends of the largest sides)
    in a thick pipe.
    From A. Mostacci, PhD thesis, LHC project note 195, LHC project report 199
    and Mathematica notebook: wwwslap.cern.ch/collective/mostacci/slots/note/slots.nb

    Lh: hole length, Wh: hole width, T: wall thickness,
    Cm and Ce are constants found from simulations, here taken as 1 (pessimistic).'''

    # length of hole of rectangular shape with same area as rounded rectangular hole
    Lr=Lh-(1-np.pi/4)*Wh;

    x=Wh/Lh;
    V=Lh*Wh**2;

    alphae=-Ce*V*np.pi/16.*(1.-0.765*x+0.1894*x**2)*np.exp(-np.pi*T*np.sqrt(1./Lr**2+1./Wh**2));
    alpham=Cm*V*np.pi/16.*(1.-0.0857*x-0.0654*x**2)*np.exp(-np.pi*T/Wh);

    return alphae,alpham;


def transverse_broadband_half_ellipsoidal_protrusion_in_round_pipe_Kurennoy(Le,We,He,b):

    '''compute transverse (broad-band) R/Q of a single small (w.r.t. wavelength) half-ellipsoid
    obstacle of semi-axes Le (along longitudinal direction), We (azimuthal)
    and He (radial - it's simply the height of the obstacle), protruding in a round beam pipe of radius b.
    Assumes Le,We,He << b.
    From S. Kurennoy, Phys. Rev. E 55 (3), pp. 3529-3532 (1997). In his notations Le=a, We=c, He=b, b=R.'''

    Z0=4e-7*np.pi*299792458; # free space impedance

    # polarizabilities for half-ellipsoid
    alphae,alpham=polarizabilities_half_ellipsoidal_protrusion_Kurennoy(Le,We,He);

    # shunt impedance /Q
    RoverQ=Z0*(alphae+alpham)/(np.pi**2*b**4);

    return RoverQ;


def longitudinal_broadband_half_ellipsoidal_protrusion_in_round_pipe_Kurennoy(Le,We,He,b,fr=5e9):

    '''compute longitudinal (broad-band) R/Q of a single small (w.r.t. wavelength) half-ellipsoid
    obstacle of semi-axes Le (along longitudinal direction), We (azimuthal)
    and He (radial - it's simply the height of the obstacle), protruding in a round beam pipe of radius b.
    Assumes Le,We,He << b.
    From S. Kurennoy, Phys. Rev. E 55 (3), pp. 3529-3532 (1997). In his notations Le=a, We=c, He=b, b=R.'''

    # fr: cutoff frequency of broad-band resonator.

    omegar=2.*np.pi*fr;
    mu0=4e-7*np.pi; # free space permeability

    # polarizabilities for half-ellipsoid
    alphae,alpham=polarizabilities_half_ellipsoidal_protrusion_Kurennoy(Le,We,He);

    # shunt impedance /Q
    RoverQ=omegar*mu0*(alphae+alpham)/(4.*np.pi**2*b**2);

    return RoverQ;


def polarizabilities_half_ellipsoidal_protrusion_Kurennoy(Le,We,He):

    '''compute polarizabilities (electric and magnetic) of a single small (w.r.t. wavelength) half-ellipsoid
    (obstacle) protruding in a round beam pipe, of semi-axes Le (along longitudinal direction), We (azimuthal)
    and He (radial - it's simply the height of the obstacle).
    Assumes Le,We,He << b (radius of beam pipe).
    From S. Kurennoy, Phys. Rev. E 55 (3), pp. 3529-3532 (1997). In his notations Le=a, We=c and He=b.'''

    # integrals in Kurennoy's paper
    Ib=integral_Kurennoy(He,Le,We)*Le*We*He/2.;
    Ic=integral_Kurennoy(We,Le,He)*Le*We*He/2.;

    alphae=2.*np.pi*Le*We*He/(3.*Ib);
    alpham=2.*np.pi*Le*We*He/(3.*(Ic-1.));

    return alphae,alpham;


def integral_Kurennoy(b,a,c):

    '''computes integral from 0 to + infinity of
    1 / ( (s+b^2)^(3/2)*(s+a^2)^(1/2)*(s+c^2)^(1/2) )
    Assumes a>0, b>0 and c>0.'''

    from scipy import special as sp
    from scipy import sqrt,arcsin # generalized functions for complex numbers
    from elliptic import ellipeinc_gen

    # particular cases first (from S. Kurennoy, Phys. Rev. E 55 (3), pp. 3529-3532 (1997) ), that I re-demonstrated
    if (a==b):
        return sp.hyp2f1(1,0.5,2.5,1.-a**2/c**2)*2./(3.*a*b*c);

    if (b==c):
        return sp.hyp2f1(1,0.5,2.5,1.-c**2/a**2)*2./(3.*a*b*c);

    # General case when a different from b, and b different from c
    # Formula directly from Mathematica (not re-demonstrated) + in some cases of
    # the parameters some formulas from Abramowitz-Stegun (pp. 592-593)

    # compute first EllipticE(Arcsin(sqrt(a^2-b^2)/a),(a^2-c^2)/(a^2-b^2))
    # (in all possible sign cases)
    E=ellipeinc_gen(arcsin(sqrt(a**2-b**2)/a),(a**2-c**2)/(a**2-b**2));

    result=2.*( (b**2-a**2)*c + a*b*sqrt(a**2-b**2)*E ) / ( a*b*(a**2-b**2)*(b**2-c**2));
    # check with numerical integration
    result_num=integral_Kurennoy_num(b,a,c);
    if (np.abs((result-result_num)/result)>1e-3):
        print(" Warning in integral_Kurennoy: numerical integration of Kurennoy's integral does not give the same result",result,result_num)

    return np.real(result);


def integrand_Kurennoy(s,b,a,c):

    '''computes 1 / ( (s+b^2)^(3/2)*(s+a^2)^(1/2)*(s+c^2)^(1/2) )
    Assumes a>0, b>0 and c>0. s can be an array.'''

    return 1. / ( (s+b**2)**(3./2.)*(s+a**2)**(1./2.)*(s+c**2)**(1./2.) );


def integrand_Kurennoy2(u,b,a,c):

    '''computes 1 / ( (1/u^2+a^2-b^2)^(1/2)*(1/u^2+c^2-b^2)^(1/2) )
    Assumes a>0, b>0 and c>0. u can be an array.'''

    return 1. / ( (1./u**2+a**2-b**2)**(1./2.)*(1./u**2+c**2-b**2)**(1./2.) );


def integral_Kurennoy_num(b,a,c):

    '''computes integral from 0 to + infinity of
    1 / ( (s+b^2)^(3/2)*(s+a^2)^(1/2)*(s+c^2)^(1/2) )
    Assumes a>0, b>0 and c>0.
    This one does it numerically.'''

    from scipy import integrate as inte
    #result,err=inte.quad(integrand_Kurennoy,0,np.inf,args=(b,a,c),limit=1000);
    # previous one (FORTRAN QUADPACK) does not work (compared to Mathematica)

    # next one uses change of variable to get non-infinite domain and use
    # Gaussian quadrature
    result,err=inte.quadrature(integrand_Kurennoy2,0,1./b,args=(b,a,c),maxiter=100);
    result *=2.;

    return result;


def longitudinal_imp_4striplinesBPM_Ng(l,angle,freq,Zc=50):

    '''compute longitudinal impedance from Ng formula as quoted in
    Handbook of Accelerator Physics and Engr., sec 3.2 (Impedances and wakes functions)

    l = length of strips (=electrodes), angle = azimuthal angle over which each of them is seen
    Zc = characterisitc impedance of transmission lines from strips (usually
    50 Ohm) and freq the frequencies over which to compute impedance
    NB: only valid for 2 pairs of striplines in symmetric position.
    '''

    k=2.*np.pi*freq/299792458;

    Zl=np.zeros((len(freq),2));

    Zl[:,0]=(np.sin(k*l))**2;
    Zl[:,1]=np.sin(k*l)*np.cos(k*l);
    Zl *= Zc*(2*angle/(np.pi))**2;

    return Zl;


def transverse_imp_4striplinesBPM_Ng(l,angle,b,freq,Zc=50):

    '''compute transverse impedance from Ng formula as quoted in
    Handbook of Accelerator Physics and Engr., sec 3.2 (Impedances and wakes functions)

    l = length of strips, angle = azimuthal angle over which each of them is seen,
    b = pipe radius (or half distance between strips),
    Zc = characterisitc impedance of transmission lines from strips (usually
    50 Ohm) and freq the frequencies over which to compute impedance
    NB: only valid for 2 pair of striplines in symmetric position.
    '''

    k=2.*np.pi*freq/299792458;

    Zt=np.zeros((len(freq),2));

    Zt[:,0] =  8.*Zc / np.pi**2 / b**2 * 1/k * np.sin(angle)**2 * np.sin(k*l)**2
    Zt[:,1] =  8.*Zc / np.pi**2 / b**2 * 1/k * np.sin(angle)**2 * np.cos(k*l) * np.sin(k*l)

    return Zt;


def transverse_wake_4striplinesBPM_Ng(l,angle,b,z,Zc=50):

    '''compute transverse wake function from Ng formula as quoted in
    Handbook of Accelerator Physics and Engr., sec 3.2 (Impedances and wakes functions)

    l = length of strips, angle = azimuthal angle over which each of them is seen,
    b = pipe radius (or half distance between strips),
    Zc = characterisitc impedance of transmission lines from strips (usually
    50 Ohm) and freq the frequencies over which to compute impedance
    NB: only valid for 2 pairs of striplines in symmetric position.

    TO BE CHECKED
    '''

    c=299792458;

    Wt=np.zeros((len(z),2));
    ind=np.where((z-2*l)*z<0)[0];

    Wt[ind,0] = 8*Zc*c * (np.sin(angle)/(b*np.pi))**2;

    return Wt;


def longitudinal_imp_striplineBPM_Ng(l,angle,freq,Zc=50):

    '''compute longitudinal impedance from Ng formula as quoted in
    Handbook of Accelerator Physics and Engr., sec 3.2 (Impedances and wakes functions)

    l = length of strips (=electrodes), angle = azimuthal angle over which each of them is seen
    Zc = characterisitc impedance of transmission lines from strips (usually
    50 Ohm) and freq the frequencies over which to compute impedance
    NB: only valid for a pair of striplines in symmetric position.
    '''

    k=2.*np.pi*freq/299792458;

    Zl=np.zeros((len(freq),2));

    Zl[:,0]=2*(np.sin(k*l))**2;
    Zl[:,1]=np.sin(2*k*l);
    Zl *= Zc*(angle/(2*np.pi))**2;

    return Zl;


def transverse_imp_striplineBPM_Ng(l,angle,b,freq,Zc=50):

    '''compute transverse impedance from Ng formula as quoted in
    Handbook of Accelerator Physics and Engr., sec 3.2 (Impedances and wakes functions)

    l = length of strips, angle = azimuthal angle over which each of them is seen,
    b = pipe radius (or half distance between strips),
    Zc = characterisitc impedance of transmission lines from strips (usually
    50 Ohm) and freq the frequencies over which to compute impedance
    NB: only valid for a pair of striplines in symmetric position.
    '''

    k=2.*np.pi*freq/299792458;

    Zl=longitudinal_imp_striplineBPM_Ng(l,angle,freq,Zc=50);

    Zt=np.zeros((len(freq),2));

    for i in range(2): Zt[:,i] = Zl[:,i] * (4*np.sin(angle/2.)/(b*angle))**2 / k;

    return Zt;


def transverse_wake_striplineBPM_Ng(l,angle,b,z,Zc=50):

    '''compute transverse wake function from Ng formula as quoted in
    Handbook of Accelerator Physics and Engr., sec 3.2 (Impedances and wakes functions)

    l = length of strips, angle = azimuthal angle over which each of them is seen,
    b = pipe radius (or half distance between strips),
    Zc = characterisitc impedance of transmission lines from strips (usually
    50 Ohm) and freq the frequencies over which to compute impedance
    NB: only valid for a pair of striplines in symmetric position.
    '''

    c=299792458;

    Wt=np.zeros((len(z),2));
    ind=np.where((z-2*l)*z<0)[0];

    Wt[ind,0] = 8*Zc*c * (np.sin(angle/2.)/(b*np.pi))**2;

    return Wt;


def resonator_impedance(R,fr,Q,freq,save=None):

    '''resonator model impedance (transverse)
    R=shunt impedance (Ohm/m), fr=resonant frequency, Q=quality factor
    computes impedance at the frequencies freq (2 columns array)
    if save different from None and contains a string, write imp. to this file name (3 columns)'''

    Z=np.zeros((len(freq),2));
    Zt=(R*fr/freq)/(1.-1j*Q*(fr/freq-freq/fr));
    Z[:,0]=np.real(Zt);Z[:,1]=np.imag(Zt);

    if (save!=None):
        data=np.hstack((freq.reshape((-1,1)),Z));
        write_ncol_file(save,data,header="Frequency[Hz]\tRe(Z)[Ohm/m]\tIm(Z)[Ohm/m]");

    return Z;


def resonator_wake(R,fr,Q,tau,save=None):

    '''resonator model wake function (transverse)
    R=shunt impedance (Ohm/m), fr=resonant frequency, Q=quality factor
    computes wake at the times tau (in ns) (must be different from 0) (tau>0 behind the source)
    NOTE: Q must be different from 0.5 !
    if save different from None and contains a string, write wake to this file name (2 columns)'''

    from scipy import sqrt,sin

    W=np.zeros((len(tau),1));tau=tau.reshape((-1,1))*1e-9;
    omegar=2*np.pi*fr;omegarbar=omegar*sqrt(1.-1./((2.*Q)**2));alpha=omegar/(2.*Q);
    ind=np.where(tau>0);
    W[ind]=(R*omegar**2/(Q*omegarbar))*np.exp(-alpha*tau[ind])*sin(omegarbar*tau[ind]);

    if (save!=None):
        data=np.hstack((tau*1e9,1.e-15*W));
        write_ncol_file(save,data,header="Tau[ns]\tWake[V/pC/mm]");

    return W;


def resonator_wake_potential(R,fr,Q,sigmaz,z,beta=1):

    '''resonator model wake potential (transverse) for a Gaussian bunch of line density
    lambda_z=exp(-z^2/(2*sigz^2)) / (sqrt(2*pi)*sigz)
    Formula from Chao's Handbook (p. 237, sec 3.2), re-derived by hand
    R=shunt impedance (Ohm/m), fr=resonant frequency, Q=quality factor, sigmaz=rms bunch length
    beta is the relativistic velocity factor
    computes wake at the distances z (in m) (must be different from 0) (z>0 behind the source)
    NOTE: Q must be different from 0.5 !
    R, fr and Q can be scalar or lists of the same length (then add all resonators)
    '''

    from scipy import special as sp
    from scipy import sqrt,exp

    c=299792458;

    W=np.zeros(len(z));

    Rlist=create_list(R,n=1);
    frlist=create_list(fr,n=1);
    Qlist=create_list(Q,n=1);

    for i in range(len(Rlist)):

        omegar=2*np.pi*frlist[i];omegarbar=omegar*sqrt(1.-1./((2.*Qlist[i])**2));
        alpha=omegar/(2.*Qlist[i]);
        kr=omegarbar/(beta*c);alphar=alpha/(beta*c);
        cst=Rlist[i]*omegar**2/(2.*Qlist[i]*omegarbar)*exp((alphar**2-kr**2)*sigmaz**2/2.);

        # intermediate z-dependent quantities
        erf=sp.erf((z-alphar*sigmaz**2+1j*kr*sigmaz**2)/(np.sqrt(2.)*sigmaz));
        expo1=exp(-alphar*z);
        expo2=exp(1j*(kr*z-alphar*kr*sigmaz**2));

        W += cst*expo1*np.imag(expo2*(1.+erf));

    return W;


def resonator_long_impedance(R,fr,Q,freq,save=None):

    '''resonator model impedance (longitudinal)
    R=shunt impedance (Ohm), fr=resonant frequency, Q=quality factor
    computes impedance at the frequencies freq (2 columns array)
    if save different from None and contains a string, write imp. to this file name (3 columns)'''

    Z=np.zeros((len(freq),2));
    Zt=R/(1.-1j*Q*(fr/freq-freq/fr));
    Z[:,0]=np.real(Zt);Z[:,1]=np.imag(Zt);

    if (save!=None):
        data=np.hstack((freq.reshape((-1,1)),Z));
        write_ncol_file(save,data,header="Frequency[Hz]\tRe(Z)[Ohm]\tIm(Z)[Ohm]");

    return Z;


def resonator_long_wake(R,fr,Q,tau,save=None):

    '''resonator model wake function (longitudinal)
    R=shunt impedance (Ohm), fr=resonant frequency, Q=quality factor
    computes wake at the times tau (in ns) (must be different from 0) (tau>0 behind the source)
    NOTE: Q must be different from 0.5 !
    if save different from None and contains a string, write wake to this file name (2 columns)'''

    from scipy import sqrt,cos,sin

    W=np.zeros((len(tau),1));tau=tau.reshape((-1,1))*1e-9;
    omegar=2*np.pi*fr;omegarbar=omegar*sqrt(1.-1./((2.*Q)**2));alpha=omegar/(2.*Q);
    ind=np.where(tau>0);
    W[ind]=(R*omegar/Q)*np.exp(-alpha*tau[ind])*(cos(omegarbar*tau[ind])-(alpha/omegarbar)*sin(omegarbar*tau[ind]));

    if (save!=None):
        data=np.hstack((tau*1e9,1.e-15*W));
        write_ncol_file(save,data,header="Tau[ns]\tWake[kV/pC]");

    return W;


def resonator_long_wake_potential(R,fr,Q,sigmaz,z,beta=1):

    '''resonator model wake potential (longitudinal) for a Gaussian bunch of line density
    lambda_z=exp(-z^2/(2*sigz^2)) / (sqrt(2*pi)*sigz)
    Formula from Chao's Handbook (p. 237, sec 3.2), partly re-derived by hand
    R=shunt impedance (Ohm), fr=resonant frequency, Q=quality factor, sigmaz=rms bunch length
    beta is the relativistic velocity factor
    computes wake at the distances z (in m) (must be different from 0) (z>0 behind the source)
    NOTE: Q must be different from 0.5 !
    R, fr and Q can be scalar or lists of the same length (then add all resonators)
    '''

    from scipy import special as sp

    c=299792458;

    W=np.zeros(len(z));

    Rlist=create_list(R,n=1);
    frlist=create_list(fr,n=1);
    Qlist=create_list(Q,n=1);

    for i in range(len(Rlist)):

        omegar=2*np.pi*frlist[i];omegarbar=omegar*np.sqrt(1.-1./((2.*Qlist[i])**2));
        alpha=omegar/(2.*Qlist[i]);
        kr=omegarbar/(beta*c);alphar=alpha/(beta*c);
        cstsin=-Rlist[i]*omegar**2/(4.*Qlist[i]**2*omegarbar);
        cstcos=Rlist[i]*omegar/(2.*Qlist[i]);
        cst=np.exp((alphar**2-kr**2)*sigmaz**2/2.);

        # intermediate z-dependent quantities
        erf=sp.erf((z-alphar*sigmaz**2+1j*kr*sigmaz**2)/(np.sqrt(2.)*sigmaz));
        expo1=np.exp(-alphar*z);
        expo2=np.exp(1j*(kr*z-alphar*kr*sigmaz**2));

        W += cst*expo1*(cstsin*np.imag(expo2*(1.+erf))+cstcos*np.real(expo2*(1.+erf)));

    return W;


def add_resonator_wake(R,fr,Q,wakefile,plane='x',save=None):

    '''add a transverse wake function from a resonator model to a wake contained in the file
    named "wakefile" (without header, first column is tau in ns, transverse wake components
    are in V/pC/mm)
    Only dipolar component of plane 'plane' is affected (second or thrid column)
    R=shunt impedance (Ohm/m), fr=resonant frequency, Q=quality factor
    if 'save' different from None and is a string, write resulting wake
    to a file with the name given by 'save'.'''

    # read initial wake file
    wake_old=read_ncol_file(wakefile);

    # compute resonator wake function
    tau=wake_old[:,0];
    wake_res=1.e-15*resonator_wake(R,fr,Q,tau);

    # add resonator wake to initial wake
    wake_new=wake_old;
    if plane.startswith('x'): iplane=1;
    else: iplane=2;
    #print "plane=",plane,", iplane=",iplane;
    wake_new[:,iplane]=wake_new[:,iplane]+wake_res.reshape((1,-1));

    if (save!=None):
        # save in file
        write_ncol_file(save,wake_new);

    return wake_new;


class impedance_wake(object):
    '''class for impedance or wake fucntion model (impedances components vs. frequency,
    or wake function vs. distance z)'''

    def __init__(self,a=1,b=0,c=0,d=0,plane='x',var=10**np.arange(1,13.1,0.1),func=np.zeros(121)):

        '''all in SI units here.

        a,b,c,d are the powers the source and test coordinates in the force (or kick) due to impedance/wake
        ex:
         - a=b=c=d=0, plane ='z' : longitudinal constant impedance/wake
         - a=b=c=d=0, plane ='x' : horizontal constant impedance/wake
         - a=b=c=d=0, plane ='y' : vertical constant impedance/wake
         - a=1, b=c=d=0, plane ='x' : horizontal dipolar impedance/wake
         - b=1, a=c=d=0, plane ='y' : vertical dipolar impedance/wake
         - c=1, a=b=d=0, plane ='x' : horizontal quadrupolar impedance/wake
         - d=1, a=b=c=0, plane ='y' : vertical quadrupolar impedance/wake
         - b=1, a=c=d=0, plane ='x' : horizontal dipolar coupled-term impedance/wake
         - d=1, a=b=c=0, plane ='x' : horizontal quadrupolar coupled-term impedance/wake
         - a=1, b=c=d=0, plane ='y' : vertical dipolar coupled-term impedance/wake
         - c=1, a=b=d=0, plane ='y' : vertical quadrupolar coupled-term impedance/wake
        etc...'''

        self.a=a; # power of x_s (source horizontal position)
        self.b=b; # power of y_s (source vertical position)
        self.c=c; # power of x_t (test horizontal position)
        self.d=d; # power of y_t (test vertical position)
        self.plane=plane; # plane of the impedance/wake ('x' for horizontal impedance/wake, 'y' for vertical, 'z' for longitudinal)
        self.var=deepcopy(var); # table with varying variable (frequencies [Hz] for impedance, z [m] for wake function)
        self.func=deepcopy(func); # 2-columns table with the corresponding impedance or wake function - real and
        # imag. parts (Ohm/m^(a+b+c+d) for impedances, V/C/m^(a+b+c+d) for wakes).
        # Usually imaginary part of wake is zero, but could be non zero (for e.g. feedback wakes).


def freqscan_from_fpar(fpar):
    '''gives frequency scan associated with a freq_param object'''

    if (np.mod(fpar.ftypescan,2)==0):
        delta=1./float(fpar.nflog);
        freq=10**(np.arange(np.log10(fpar.fmin),np.log10(fpar.fmax)+delta,delta));
        if fpar.ftypescan==2:
            delta=(fpar.fmaxrefine-fpar.fminrefine)/fpar.nrefine;
            freq2=np.arange(fpar.fminrefine,fpar.fmaxrefine+delta,delta);
            freq=np.concatenate((freq,freq2));
    else: # fpar.ftypescan==1
        freq=np.arange(fpar.fmin,fpar.fmax+fpar.fsamplin,fpar.fsamplin);

    freq=np.concatenate((freq,fpar.fadded));

    return sort_and_delete_duplicates(freq);


def zscan_from_zpar(zpar):
    '''gives z (distance) scan associated with a z_param object'''

    if (np.mod(zpar.ztypescan,2)==0):
        delta=1./float(zpar.nzlog);
        z=10**(np.arange(np.log10(zpar.zmin),np.log10(zpar.zmax)+delta,delta));
        if zpar.ztypescan==2:
            z2=np.arange(zpar.zminrefine,zpar.zmaxrefine+zpar.zsamplin,zpar.zsamplin);
            z=np.concatenate((z,z2));
    else: # zpar.ztypescan==1
        z=np.arange(zpar.zminrefine,zpar.zmaxrefine+zpar.zsamplin,zpar.zsamplin);

    z=np.concatenate((z,zpar.zadded));

    return sort_and_delete_duplicates(z);


def test_impedance_wake_comp(iw,a,b,c,d,plane):
    '''test if the impedance or wake in iw has these a,b,c,d and plane'''
    flag=((iw.a==a)and(iw.b==b))and((iw.c==c)and(iw.d==d));

    return (flag)and(iw.plane==plane);


def add_impedance_wake(iw_model,iw_added,weightx,weighty,print_warning=True):

    '''add the model "iw_added" to an impedance or wake model in "iw_model".
    both iw_added and iw_model are impedance or wake models, i.e. a list of objects from the class
    impedance_wake. The two lists are not necessarily of same length, and the frequency/z scans
    are not necessarily all the same.
    weightx and weighty are resp. betax/avbetax and betay/avbetay (i.e. without the square root)
    at the place of the added impedance "iw_added".

    It does this "in-place" (on iw_model itself)'''

    # loop over all impedance objects in iw_added
    for iw in iw_added:

        # compute the impedance weight
        powx=iw.a/2.+iw.c/2.;
        powy=iw.b/2.+iw.d/2.;
        if iw.plane.startswith('x'): powx +=1./2.;
        elif iw.plane.startswith('y'): powy +=1./2.;

        weight=(weightx**powx) * (weighty**powy);

        # determine if in iw_model_new there is already the same kind of component as in iw.
        flag=False;
        for kiw_mod,iw_mod in enumerate(iw_model):
            if test_impedance_wake_comp(iw_mod,iw.a,iw.b,iw.c,iw.d,iw.plane):
                flag=True;ksamecomp=kiw_mod;

        if flag:
            # same kind of component as in imp already present in iw_model_new
            iw_mod=iw_model[ksamecomp];
            # concatenate the 2 var scans
            vartot=sort_and_delete_duplicates(np.concatenate((iw.var,iw_mod.var)));
            if print_warning and (vartot[-1]>iw.var[-1])and(np.abs((iw.func[-1,0]+1j*iw.func[-1,1])/np.average(iw.func[:,0]+1j*iw.func[:,1]))>1e-2):
                print(" Warning in add_impedance_wake: added model extrapolated above",iw.var[-1],"by zeros; func[-1]=",iw.func[-1,0],iw.func[-1,1]);
            if print_warning and (vartot[-1]>iw_mod.var[-1])and(np.abs((iw_mod.func[-1,0]+1j*iw_mod.func[-1,1])/np.average(iw_mod.func[:,0]+1j*iw_mod.func[:,1]))>1e-2):
                print(" Warning in add_impedance_wake: initial model extrapolated above",iw_mod.var[-1],"by zeros; func[-1]=",iw_mod.func[-1,0],iw_mod.func[-1,1]);
            # sum impedances
            functot=np.zeros((len(vartot),2),dtype=float);
            for icol in range(2):
                functot[:,icol]=np.interp(vartot,iw_mod.var,iw_mod.func[:,icol],right=0.)+weight*np.interp(vartot,iw.var,iw.func[:,icol],right=0.);
            # include it in the new model
            iw_mod.var=vartot;
            iw_mod.func=functot;

        else:
            # same kind of component as in iw is not already present:
            # append to iw_model_new a new impedance/wake with this component (with the correct weight)
            iw2=deepcopy(iw);iw2.func=weight*iw.func;
            iw_model.append(iw2);

    return iw_model;


def multiply_impedance_wake(iw_model,factor):

    '''multiply all components of the model "iw_model" by the factor "factor"
    It does this "in-place" (on iw_model itself)'''

    # loop over all impedance objects in iw_model
    for iw in iw_model:

        iw.func *= factor;

    return iw_model;


def multiply_impedance_wake_file(iw_model,factor_filename):

    '''multiply components of the model "iw_model" by the factors
    (freq. or distance dependent) in the file "factor_filename"
    headers of the file should provide info on the component to multiply
    (see format in function 'identify_component') and real / imag. part
    (at the end of headers: re or Re / im or Im)
    NOTE: frequency or distance scan should be LARGER in the factor file than
    in the imp. or wake model

    It does it "in-place" (on iw_model itself)'''

    # read headers
    fid=open(factor_filename);
    headers=fid.readline();heads=headers.split();
    if not ((heads[0].startswith('freq'))or(heads[0].startswith('Freq')))or((heads[0].startswith('dist'))or(heads[0].startswith('Dist'))):
        print("Pb in multiply_impedance_wake_file: frequencies or distances not found in 1st column!");
        sys.exit();
    fid.close();

    # read file
    factor=read_ncol_file(factor_filename,ignored_rows=1);

    # apply factors for each component
    for icomp,comp in enumerate(heads[1:]):
        # identify component type
        a,b,c,d,plane,wakeflag=identify_component(comp);
        # identify if real or imag part (default=real)
        if (comp.endswith('re'))or(comp.endswith('Re')): ir=0;
        elif (comp.endswith('im'))or(comp.endswith('Im')): ir=1;
        else: i=0;

        # loop over all impedance objects in iw_model and apply factor (interpolated)
        for iw in iw_model:

            if test_impedance_wake_comp(iw,a,b,c,d,plane):
                iw.func[:,ir] *= np.interp(iw.var,factor[:,0],factor[:,icomp+1]);


    return iw_model;


def rotate_imp_wake(iw_mod,theta):

    ''' Rotate in the transverse plane an impedance or wake model by the angle theta
    Initial model is in iw_mod (list of impedance/wake objects from class impedance_wake),
    given on the initial basis (x',y').
    Output is an impedance/wake model (list of impedance_wake objects)
    of the same element in the axes (x,y) rotated by theta such that
    theta = angle between e_x and e_x' (counted positively anticlockwise
    from x to x')
    Note: for LHC collimators, theta=pi/2-alpha, with alpha the skew
    angle of LHC collimators.'''
    from scipy.misc import comb

    iw_mod_rotated=[];
    costheta=np.cos(theta);sintheta=np.sin(theta);

    for iw in iw_mod:
        # loop over all components of model iw_mod
        if (iw.plane.startswith('x')): coefx=costheta;coefy=sintheta;
        elif (iw.plane.startswith('y')): coefx=-sintheta;coefy=costheta;

        for i in range(iw.a+1):
            for j in range(iw.b+1):
                for k in range(iw.c+1):
                    for l in range(iw.d+1):
                        # product of binomial coefficients
                        binprod=int(comb(iw.a,i,exact=1)*comb(iw.b,j,exact=1)*comb(iw.c,k,exact=1)*comb(iw.d,l,exact=1));
                        # multiply by power of cos and sin
                        coef=(-1)**(j+l)*binprod*costheta**(i+iw.b-j+k+iw.d-l)*sintheta**(iw.a-i+j+iw.c-k+l);
                        # add impedance/wake term to model
                        if (iw.plane.startswith('z')):
                            # case when initial imp/wake is longitudinal
                            add_impedance_wake(iw_mod_rotated,[impedance_wake(a=i+j,b=iw.a-i+iw.b-j,
                                    c=k+l,d=iw.c-k+iw.d-l,plane='z',var=iw.var,func=coef*iw.func)],1,1);
                        else:
                            # case when initial imp/wake is transverse
                            add_impedance_wake(iw_mod_rotated,[impedance_wake(a=i+j,b=iw.a-i+iw.b-j,
                                    c=k+l,d=iw.c-k+iw.d-l,plane='x',var=iw.var,func=coefx*coef*iw.func)],1,1);
                            add_impedance_wake(iw_mod_rotated,[impedance_wake(a=i+j,b=iw.a-i+iw.b-j,
                                    c=k+l,d=iw.c-k+iw.d-l,plane='y',var=iw.var,func=coefy*coef*iw.func)],1,1);

    return iw_mod_rotated;


def Yokoya_elliptic(small_semiaxis,large_semiaxis,Yokoyafilename=Yokoya_path+'/Yokoya_elliptic_from_Elias_USPAS.dat'):

    '''compute Yokoya factor for an elliptic shape (semiaxis given in input),
    assuming the small axis is in the vertical direction.'''

    # read Yokoya factors interpolation file
    # BEWARE: columns are ratio, dipy, dipx, quady, quadx
    yokoya_file=read_ncol_file(Yokoyafilename,ignored_rows=1);
    ratio_col=yokoya_file[:,0];
    # compute semi-axes ratio (first column of this file)
    ratio=(large_semiaxis-small_semiaxis)/(large_semiaxis+small_semiaxis);

    # interpolate Yokoya file at the correct ratio
    yoklong=1.;
    yokydip=np.interp(ratio,ratio_col,yokoya_file[:,1]);
    yokxdip=np.interp(ratio,ratio_col,yokoya_file[:,2]);
    yokyquad=np.interp(ratio,ratio_col,yokoya_file[:,3]);
    yokxquad=np.interp(ratio,ratio_col,yokoya_file[:,4]);

    return [yoklong,yokxdip,yokydip,yokxquad,yokyquad];


def apply_Yokoya_elliptic(iw_model,b,w,orientation='V'):

    '''apply Yokoya factors to a impedance or wake model
    small semi-axis is b, large one w, and orientation 'H' (horizontal)
    or 'V' (vertical) defined as the axis along which the semi-axis 'b' is.
    NOTE: initial model is assumed to be "round", i.e. only longitudinal
    and dipolar terms are taken into account (and xdip assumed to be equal to
    ydip). All other terms are ignored (in particular those of higher order).'''

    if (b==w): return iw_model;
    else:
        iw_model_new=[];
        yok=Yokoya_elliptic(b,w);
        if (orientation=='H')or(orientation=='h'):
            # invert planes in Yokoya factors
            yokold2=yok[2];yokold4=yok[4];
            yok[2]=yok[1];
            yok[1]=yokold2;
            yok[4]=yok[3];
            yok[3]=yokold4;

        # longitudinal
        for iw in iw_model:
            if test_impedance_wake_comp(iw,0,0,0,0,'z'):
                iw_model_new.append(impedance_wake(a=0,b=0,c=0,d=0,plane='z',var=iw.var,func=yok[0]*iw.func))
        for iw in iw_model:
            if test_impedance_wake_comp(iw,1,0,0,0,'x'):
                # x-dipolar
                iw_model_new.append(impedance_wake(a=1,b=0,c=0,d=0,plane='x',var=iw.var,func=yok[1]*iw.func))
                # y-dipolar
                iw_model_new.append(impedance_wake(a=0,b=1,c=0,d=0,plane='y',var=iw.var,func=yok[2]*iw.func))
                # x-quadrupolar
                iw_model_new.append(impedance_wake(a=0,b=0,c=1,d=0,plane='x',var=iw.var,func=yok[3]*iw.func))
                # y-quadrupolar
                iw_model_new.append(impedance_wake(a=0,b=0,c=0,d=1,plane='y',var=iw.var,func=yok[4]*iw.func))


    return iw_model_new;


def readZ(filename):
    '''read impedance file (3 columns: frequency, real and imag. parts - there can be a header which is ignored)'''
    freq=[];Z=[];
    fh=open(filename,"r")
    for l in fh.readlines():
        if (not l.startswith('Freq'))and(not l.startswith('freq')):
            ll=l.strip().split();
            if (len(freq)==0)or(float(ll[0])!=freq[-1]): # avoid duplicates
                freq.append(float(ll[0]));
                Z.append(float(ll[1]));
                Z.append(float(ll[2]));
    fh.close()
    Z=np.array(Z);
    return np.array(freq),Z.reshape((-1,2)); # two columns format for Z


def imp_model_from_IW2D(iw_input,wake_calc=False,path=IW2D_path,
        flagrm=True,lxplusbatch=None,queue='1nh',dire='', job_flavour=None):

    '''create an impedance model (without beta functions weight) from an ImpedanceWake2D (IW2D) computation.
    input parameters are in iw_input (object of the class impedance_wake_input).
    if wakecalc=True, do a wake+imp calculation  (otherwise computes only imp.)
    if flagrm=True, remove output files automatically
    path contains the path to the ImpedaceWake2D executables.
    lxplusbatch: if None, no use of lxplus batch system
                   if 'launch' -> launch calculation on lxplus on queue 'queue'
                   if 'retrieve' -> retrieve outputs
    dire contains the directory name where to put the outputs (default='./'=directory of IW2D)
    NOTE: input should be clearly identified by its 'comment' '''

    def getstatusoutput(cmd):
        '''
        This emulates the obsolete function "commands.getstatusoutput"
        '''
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        out, _ = process.communicate()
        return (process.returncode, out.strip().decode())

    def create_submission_file_HTCondor(submission_filename='batch.job',
                    executable='',initialdir='',input_filename='',job_flavour='tomorrow'):
        """
        Create a submission file for HTCondor
        """
        here = os.getcwd()
        with open(submission_filename,'w') as filejob:
            os.system("cp "+os.path.join(here,input_filename)+" "+dire)
            filejob.write("executable = "+executable+"\n")
            filejob.write("input = "+input_filename+"\n")
            filejob.write("ID = $(Cluster).$(Process)\n")
            filejob.write("output = $(ID).out\n")
            filejob.write("error = $(ID).err\n")
            filejob.write("log = $(Cluster).log\n")
            filejob.write("universe = vanilla\n")
            filejob.write("initialdir = "+initialdir+"\n")
            # filejob.write("should_transfer_file = YES\n")
            filejob.write("when_to_transfer_output = ON_EXIT\n")
            # filejob.write("transfer_input_file = "+filename+"\n")
            filejob.write("+JobFlavour =\""+job_flavour+"\"\n")
            filejob.write("queue\n")

        return

    imp_model = []
    wake_model = []

    # list of different components
    listname = ['long','xdip','ydip','xquad','yquad','ycst']
    listplane = ['z','x','y','x','y','y']
    lista = [0,1,0,0,0,0]
    listb = [0,0,1,0,0,0]
    listc = [0,0,0,1,0,0]
    listd = [0,0,0,0,1,0]

    # find number of upper layers (last upper layer has infinite thickness)
    nup = 0
    while (getattr(iw_input.layers[nup],'thickness')!=np.inf):
        nup += 1
    nup += 1

    if (len(iw_input.b)>=2):
        # include cst y component (no top-bottom symmetry)
        ncomp = 6
    else:
        # top-bottom symmetry -> Zycst=Wycst=0
        ncomp = 5

    # save path and go to executable path
    oldpath = os.getcwd()
    # commands.getoutput("pwd")
    os.chdir(path)

    # write input file
    filename = 'IW2D_input'+iw_input.comment+'.txt'
    write_impedance_wake_input(filename, iw_input)

    if wake_calc:
        wakestr = 'wake_'
    else:
        wakestr=''

    os.system("mkdir -p "+dire)

    if lxplusbatch and job_flavour is None:
        # Equivalent on HTCondor to LSF queues
        if queue in ['8nm','1nh','8nh','1nd']:
            job_flavour = 'tomorrow'
        elif queue in ['2nd', '1nw','2nw']:
            job_flavour = 'nextweek'
        else:
            print('Error in LSF queue, can not find HTCondor equivalent, use nextweek as default')
            job_flavour = 'nextweek'


    if (lxplusbatch==None):
        # launch computation
        os.system("./"+wakestr+iw_input.geometry+"chamber.x < "+filename+" > outtmp")
        os.system(" mv *"+iw_input.comment+"*.dat "+dire)

    elif lxplusbatch.startswith('launch'):

        # Launch calculation on HTCondor batch system
        # First create the submission file
        create_submission_file_HTCondor(submission_filename='batch'+iw_input.comment+'.job',
                                        executable=wakestr+iw_input.geometry+"chamber.x",
                                        initialdir=dire,input_filename=filename,job_flavour=job_flavour)

        # Submit the job to HTCondor
        # os.system("chmod 744 batch"+iw_input.comment+".job")
        # os.system("bsub"+user_option+" -e error1.out -q "+queue+" batch"+iw_input.comment+".job")
        os.system("condor_submit batch"+iw_input.comment+".job")

    if (lxplusbatch==None) or (lxplusbatch.startswith('retrieve')):
        # extract output suffix
        deltab = [0.,0.005,-0.005,0.01,-0.01]
        status = 1
        ib = 0

        # scan closed-by radii to the one looked for
        while ((status!=0) and (ib<len(deltab))):
            fb = iw_input.b[0]*1e3 + deltab[ib]
            strb="%.2lf" % fb
            #print "ls "+dire+"InputData*"+str(nup)+'layer*'+strb+"mm"+iw_input.comment+".dat"
            status,name = getstatusoutput("ls "+dire+"InputData*"+str(nup)+'layer*'+strb+"mm"+iw_input.comment+".dat")
            if status == 0:
                # additional check that input files are really identical
                status = 0 if filecmp.cmp(name,filename) else 1
            #print strb,status,name
            #status,name=getstatusoutput("ls "+dire+"InputData*"+str(nup)+'layer*'+iw_input.comment+".dat")
            ib += 1
        
        if (status!=0):
            print(" No impedance file -> Probably computation failed... ")
            if lxplusbatch is not None:
                print(" We try to relaunch")
                create_submission_file_HTCondor(submission_filename='batch'+iw_input.comment+'.job',
                                                     executable=wakestr+iw_input.geometry+"chamber.x",
                                                     initialdir=dire,input_filename=filename,job_flavour=job_flavour)
                os.system("condor_submit batch"+iw_input.comment+".job")
            # go back to previous path and raise an error
            os.chdir(oldpath)
            raise ImpedanceWakeFileNotFoundError("File",path+'/'+dire+"InputData*"+str(nup)+'layer*'+iw_input.comment+".dat not found !")

        # print iw_input.b[0]*1e3, fb, strb, status, name
        suffix1 = name.replace(dire+'InputData', '').split('\n')[0] # take the first one when there are several...

        # change output file name
        if (lxplusbatch==None):
            os.system("mv outtmp "+dire+"out"+suffix1)

        # use most precise calculation in case of wake
        if wake_calc:
            suffix = suffix1.replace('.dat','_precise.dat')
        else:
            suffix = suffix1

        # read impedance outputs and put them in model(s)
        # loop over the different components (long., x dip., etc.)
        for icomp in range(ncomp):
            namecomp = listname[icomp]
            freq, Z = readZ(dire+"Z"+namecomp+suffix)

            if len(Z)==0:
                print("    Impedance file empty -> Probably computation failed... ")
                if lxplusbatch is not None:
                    print("    We try to relaunch")
                    create_submission_file_HTCondor(submission_filename='batch'+iw_input.comment+'.job',
                                                   executable=wakestr+iw_input.geometry+"chamber.x",
                                                   initialdir=dire,input_filename=filename,job_flavour=job_flavour)
                    os.system("condor_submit batch"+iw_input.comment+".job")
                # go back to previous path and raise an error
                os.chdir(oldpath)
                raise ImpedanceWakeFileUnreadableError('File Z'+namecomp+suffix+' has no data!')

            imp_model.append(impedance_wake(a=lista[icomp],
                                            b=listb[icomp],
                                            c=listc[icomp],
                                            d=listd[icomp],
                                            plane=listplane[icomp],
                                            var=freq,func=Z))

        if wake_calc:
        # read wake outputs and put them in model(s)
            # loop over the different components (long., x dip., etc.)
            for icomp in range(ncomp):
                namecomp = listname[icomp]
                s = read_ncol_file(dire+"W"+namecomp+suffix, ignored_rows=1)
                W = np.zeros((len(s[:,0]),2))
                W[:,0] = s[:,1]
                wake_model.append(impedance_wake(a=lista[icomp],
                                                 b=listb[icomp],
                                                 c=listc[icomp],
                                                 d=listd[icomp],
                                                 plane=listplane[icomp],
                                                 var=s[:,0],func=W))

        # print "rm *"+suffix1.replace(".dat","*.dat")
        if flagrm:
            # os.system("rm -f "+dire+"*"+suffix1.replace(".dat","*.dat"))
            os.system("rm -rf LSFJOB_* error1.out batch"+iw_input.comment+".job "+filename)


    #elif (lxplusbatch.startswith('relaunch')):
        ## extract output suffix
        #deltab = [0.,0.005,-0.005,0.01,-0.01]
        #status = 1
        #ib = 0

        ## scan closed-by radii to the one looked for
        #while ((status!=0)and(ib<len(deltab))):
            #fb = iw_input.b[0]*1e3 + deltab[ib]
            #strb = "%.2lf" % fb
            ## status,name = commands.getstatusoutput("ls "+dire+"InputData*"+str(nup)+'layer*'+strb+"mm"+iw_input.comment+".dat")
            #status,name = commands.getstatusoutput("ls "+dire+"Zxquad*"+str(nup)+'layer*'+iw_input.comment+".dat")
            #ib += 1

        #if (status!=0):
            #print "\n Error: File",path+'/'+dire+"InputData*"+str(nup)+'layer*'+iw_input.comment+".dat not found !"
            #print " Probably computation failed; try with a longer queue (current queue:",queue+")"
            #os.system("bsub"+user_option+" -e error1.out -q "+queue+" batch"+iw_input.comment+".job")


    # go back to previous path
    os.chdir(oldpath)

    return imp_model, wake_model


def imp_model_elliptic(iw_input,w,orientation='V',wake_calc=False,
        flagrm=True,lxplusbatch=None,queue='1nh',dire='', job_flavour=None):

    '''small wrapper to compute wall impedance of an elliptic element
    of smaller semi-axis iw_input.b, larger one w, and orientation 'H' (horizontal)
    or 'V' (vertical) defined as the axis along which the semi-axis 'b' is.
    see other parameters in function imp_model_from_IW2D
    use round chamber theory (ImpedanceWake2D) + Yokoya factors
    dire contains the directory name where to put the outputs (default='./'=directory of IW2D)'''

    iw_input.geometry='round';
    if (iw_input.b[0]==w): iw_input.Yokoya=[1,1,1,0,0];
    else:
        yok=Yokoya_elliptic(iw_input.b[0],w);
        if (orientation=='H')or(orientation=='h'):
            # invert planes in Yokoya factors
            yokold2=yok[2];yokold4=yok[4];
            yok[2]=yok[1];
            yok[1]=yokold2;
            yok[4]=yok[3];
            yok[3]=yokold4;
        iw_input.Yokoya=yok;

    #print iw_input.b,w,iw_input.Yokoya;
    imp_model,wake_model=imp_model_from_IW2D(iw_input,wake_calc=wake_calc,flagrm=flagrm,lxplusbatch=lxplusbatch,queue=queue,dire=dire,job_flavour=job_flavour)

    return imp_model,wake_model;


def imp_model_transverse_resonator(R,fr,Q,beta=1,wake_calc=False,
        fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param(),plane=''):

    '''define an imp/wake model for a transverse resonator (SI units)
    R, fr, Q: shunt impedance, resonance frequency and quality factor
    beta: relativistic velocity factor (for z to tau conversion
    if wake_calc==True, calculate wake also
    assume axisymmetry (purely dipolar)
    plane can be 'x' (only to dipolar x component), 'y' (only to dipolar y component),
    or '' (both x and y dipolar components)'''

    imp_mod=[];wake_mod=[];

    freq=freqscan_from_fpar(fpar);
    Zt=resonator_impedance(R,fr,Q,freq);
    if (plane=='')or(plane=='x'): imp_mod.append(impedance_wake(a=1,b=0,c=0,d=0,plane='x',var=freq,func=Zt));
    if (plane=='')or(plane=='y'): imp_mod.append(impedance_wake(a=0,b=1,c=0,d=0,plane='y',var=freq,func=Zt));

    if wake_calc:
        z=zscan_from_zpar(zpar);
        tau=z*1e9/(beta*299792458); # conversion in ns
        Wt=np.zeros((len(z),2));
        tmp=resonator_wake(R,fr,Q,tau);
        Wt[:,0]=tmp[:,0];
        if (plane=='')or(plane=='x'): wake_mod.append(impedance_wake(a=1,b=0,c=0,d=0,plane='x',var=z,func=Wt));
        if (plane=='')or(plane=='y'): wake_mod.append(impedance_wake(a=0,b=1,c=0,d=0,plane='y',var=z,func=Wt));

    return imp_mod,wake_mod;


def imp_model_longitudinal_resonator(R,fr,Q,beta=1,wake_calc=False,
        fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param()):

    '''define an imp/wake model for a longtudinal resonator (SI units)
    R, fr, Q: shunt impedance, resonance frequency and quality factor
    beta: relativistic velocity factor (for z to tau conversion
    if wake_calc==True, calculate wake also'''

    imp_mod=[];wake_mod=[];

    freq=freqscan_from_fpar(fpar);
    Zl=resonator_long_impedance(R,fr,Q,freq);
    imp_mod.append(impedance_wake(a=0,b=0,c=0,d=0,plane='z',var=freq,func=Zl));

    if wake_calc:
        z=zscan_from_zpar(zpar);
        tau=z*1e9/(beta*299792458); # conversion in ns
        Wl=np.zeros((len(z),2));
        tmp=resonator_long_wake(R,fr,Q,tau);
        Wl[:,0]=tmp[:,0]
        wake_mod.append(impedance_wake(a=0,b=0,c=0,d=0,plane='z',var=z,func=Wl));

    return imp_mod,wake_mod;


def imp_model_resonator(Rlist,frlist,Qlist,beta=1,wake_calc=False,
        fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param(),listcomp=['Zxdip','Zydip']):

    '''function doing the same as the two previous ones but in a more general way.

    listcomp contains the component list in the format defined by 'identify_component';
    the model will contain a resonator with the same parameters for each component of the list.
    The first letter (W or Z) of each component does not matter: wake computation
    is determined by the wake_calc flag.
    It includes the longitudinal case (if component begins by Zl or Wl).

    Rlist, frlist and Qlist can be scalar or lists (of the size of listcomp); if they are lists then
    each component will have a different R, fr and Q.
    all units are SI.'''

    imp_mod=[];wake_mod=[];

    freq=freqscan_from_fpar(fpar);
    z=zscan_from_zpar(zpar);
    tau=z*1e9/(beta*299792458); # conversion in ns

    # test if Rlist, frlist and Qlist are list or scalar; if they are scalar replace them
    # by lists (if length that of listcomp) with repeated elements
    n=len(listcomp);
    Rlist=create_list(Rlist,n);
    frlist=create_list(frlist,n);
    Qlist=create_list(Qlist,n);

    if ((len(Rlist)!=n)or(len(frlist)!=n))or(len(Qlist)!=n):
        print("imp_model_resonator: problem in size of lists Rlist, frlist or Qlist",len(Rlist),len(frlist),len(Qlist),n)

    for icomp,comp in enumerate(listcomp):
        a,b,c,d,plane,wakeflag=identify_component(comp);
        R=Rlist[icomp];fr=frlist[icomp];Q=Qlist[icomp];
        if Q!=0:
            if (comp[1]=='l'): strlong='_long';
            else: strlong='';

            Z=eval('resonator'+strlong+'_impedance(R,fr,Q,freq)');
            imp_mod.append(impedance_wake(a=a,b=b,c=c,d=d,plane=plane,var=freq,func=Z));

            if wake_calc:
                W=np.zeros((len(z),2));
                tmp=eval('resonator'+strlong+'_wake(R,fr,Q,tau)')
                W[:,0]=tmp[:,0]
                wake_mod.append(impedance_wake(a=a,b=b,c=c,d=d,plane=plane,var=z,func=W));

    return imp_mod,wake_mod;


def imp_model_holes_Mostacci(Lh,Wh,T,b,d,eta,rhob,rhod,length,
        fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param(),fcutoff=5e9,
        nb_holes_per_crosssection=8,Cm=1.,Ce=1.,beta=1.):

    '''define an impedance & wake model for rounded rectangular holes of length Lh and width Wh
    in a round beam pipe of length 'length', inner radius b, thickness T and resistivity rhob,
    itself inside a larger round pipe of radius d and resistivity rhod,
    with the eta the fraction of the total surface covered by holes,
    nb_holes_per_crosssection the number of "lines of holes", assuming
    the holes are uniformly distributed (pessimistic), and with a cutoff
    frequency for the impedance equal to fcutoff (default = LHC value for
    an aperture of 18.4mm).
    The frequency scan is defined by fpar, the distacne scan by zpar.
    Cm and Ce are constants found from simulations, here taken as 1 (pessimistic).
    From A. Mostacci PhD thesis + Kurennoy results (Part. Acc. 1995, vol. 50, p 167-175)
    for broadband (long. and transverse).
    beta is the relativistic velocity factor
    NOTE: real long. impedance due to power dissipated by the wave propagating
    behind the holes is not taken into account in the long. wake (but is in the impedance)'''

    imp_mod=[];wake_mod=[];
    freq=freqscan_from_fpar(fpar);
    z=zscan_from_zpar(zpar);
    tau=z*1e9/(beta*299792458.); # conversion in ns

    # longitudinal model from Mostacci (including broadband from Kurennoy for imag. part)
    Zl=longitudinal_imp_holes_in_round_pipe_Mostacci(freq,Lh,Wh,T,b,d,eta,rhob,rhod,length,
        fcutoff=fcutoff,nb_holes_per_crosssection=nb_holes_per_crosssection,Cm=Cm,Ce=Ce);

    imp_mod.append(impedance_wake(a=0,b=0,c=0,d=0,plane='z',var=freq,func=Zl));

    # longitudinal wake
    Rl=longitudinal_broadband_holes_in_round_pipe_Kurennoy(Lh,Wh,b,eta,fr=fcutoff)*length;
    Wl=np.zeros((len(z),2));
    tmp=resonator_long_wake(Rl,fcutoff,1.,tau)
    Wl[:,0]=tmp[:,0]

    wake_mod.append(impedance_wake(a=0,b=0,c=0,d=0,plane='z',var=z,func=Wl));

    # transverse broadband model from Kurennoy
    Rt=transverse_broadband_holes_in_round_pipe_Kurennoy(Lh,Wh,b,eta)*length;
    Zt=resonator_impedance(Rt,fcutoff,1.,freq);
    Wt=np.zeros((len(z),2));
    tmp=resonator_wake(Rt,fcutoff,1.,tau);
    Wt[:,0]=tmp[:,0]
    print(Rt)

    # transverse wake
    imp_mod.append(impedance_wake(a=1,b=0,c=0,d=0,plane='x',var=freq,func=Zt));
    imp_mod.append(impedance_wake(a=0,b=1,c=0,d=0,plane='y',var=freq,func=Zt));

    wake_mod.append(impedance_wake(a=1,b=0,c=0,d=0,plane='x',var=z,func=Wt));
    wake_mod.append(impedance_wake(a=0,b=1,c=0,d=0,plane='y',var=z,func=Wt));

    return imp_mod,wake_mod;


def imp_model_from_HOMfile(filename,beta=1,unit_freq='Hz',spread=False, spread_rng=3e6, seed=False,damp_fm=False,fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param()):

    '''Create an impedance or wake model from a file containing HOMs (high
    order modes) data, i.e. resonator models data. Impedance / wake is in
    the end of sum of all those modes.
    In the file each line is one mode, and each column should have a header indicating:
    - Rl: longitudinal shunt impedance
    - Rxd(ip): horizontal dipolar impedance shunt impedance
    - same with Ryd(ip), Rxq(uad), Ryq(uad)
    - same with Ql, Qxd(ip), etc. for each quality factor
    - same with fl, fxd(ip), etc. for each resonance frequency.
    Input:
    filename: the HOM filename to model.
    beta: indicates the relativistic velocity factor.
    unit_freq: units in the filename HOM frequency
    spread: if True applies a uniform spread within (spread_rng)
    spread_rng: defines max frequency spread of the HOM (uniform distribution +/- spread_rng)
    seed: seed for random generator (True or False)
    damp_fm: if True, damp the first mode (fundamental mode) putting R=R/Q and Q=1. It applies only in transverse.
    fpar: frequency range (it is then expanded with HOMs)
    zpar: position range.
    '''

    colname=['Rl','Rxd','Ryd','Rxq','Ryq'];
    Rlist=[];frlist=[];Qlist=[];complist=[];

    # read file
    for icol,col in enumerate(colname):


        R=read_ncol_file_identify_header(filename,col,dispflag=False);
        fr=read_ncol_file_identify_header(filename,col.replace('R','f'),dispflag=False);
        Q=read_ncol_file_identify_header(filename,col.replace('R','Q'),dispflag=False);
        if spread:
            if seed:
                np.random.seed(seed)
            deltaf=np.random.uniform(-spread_rng,spread_rng,len(fr))
            fr+=deltaf

        if (damp_fm) and (icol!=0) and (len(R)>0):
            print('Damp fundamental mode')
            R[0]=R[0]/Q[0]
            Q[0]=1

        if (len(R)*len(fr)*len(Q)>0):
            complist.append(col.replace('R','Z'));
            Rlist.append(R);frlist.append(fr);Qlist.append(Q);

    if (len(complist)==0): print("Pb in imp_model_from_HOMfile: no valid data in HOM file!");sys.exit();

    # construct model
    imp_mod=[];wake_mod=[];

    # sum all modes
    for i in range(len(Rlist[0])):
            # compute resonator model for each mode
        Rlistcomp=[Rlist[icomp][i] for icomp in range(len(complist))];
        Qlistcomp=[Qlist[icomp][i] for icomp in range(len(complist))];
        if unit_freq=='MHz':
            fact=1e6;
        elif unit_freq=='GHz':
            fact=1e9;
        else:
            fact=1


        frlistcomp=[frlist[icomp][i]*fact for icomp in range(len(complist))];

        # define frequency axis


        freq=pow(10,np.arange(1,15,1/10.0));
        for kk in range(len(frlistcomp)):
            fres=frlistcomp[kk];
            Q=Qlistcomp[kk];
            if Q!=0:
                D=fres/2/Q;
                N=50;
                freq_reso=np.array([fres]);
                for nn in 10**(np.arange(0,11,2.)/2):
                    frmin=np.min(fres);frmax=np.max(fres);
                    freq_reso=np.concatenate([freq_reso,np.arange(frmin-nn*D,frmin,(nn*D)/N)]);
                    freq_reso=np.concatenate([freq_reso,np.arange(frmax,frmax+nn*D,(nn*D)/N)]);
                freq_reso=np.sort(np.unique(freq_reso));
                freq=np.sort(np.concatenate([freq,freq_reso]));
                freq=np.unique(freq[freq>0])


        fpar=freq_param(ftypescan=1,fmin=1e6,fmax=1e6,fsamplin=1e6,fadded=freq);

        imp,wake=imp_model_resonator(Rlistcomp,frlistcomp,Qlistcomp,beta=beta,wake_calc=True,fpar=fpar,zpar=zpar,listcomp=complist);

        add_impedance_wake(imp_mod,imp,1,1);
        add_impedance_wake(wake_mod,wake,1,1);


    return imp_mod,wake_mod;


def imp_model_striplineBPM_Ng(l,angle,b,Zc=50.,beta=1,plane='x',wake_calc=False,
        fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param()):

    '''compute impedance and wake model from Ng formula as quoted in
    Handbook of Accelerator Physics and Engr., sec 3.2 (Impedances and wakes functions)
    (transverse and longitudinal)

    Note that the longitudinal wake is zero because it's expressed as delta function
    in Ng's formulas (to be checked).

    - l = length of strips, angle = azimuthal angle over which each of them is seen,
    - b = pipe radius (or half distance between strips),
    - Zc = characterisitc impedance of transmission lines from strips (usually
    50 Ohm)
    - fpar / zpar:frequencies / distances over which to compute impedance /wake
    NB: valid only for one pair of striplines  '''

    imp_mod=[];wake_mod=[];

    freq=freqscan_from_fpar(fpar);
    z=zscan_from_zpar(zpar);

    Zl=longitudinal_imp_striplineBPM_Ng(l,angle,freq,Zc=50);
    imp_mod.append(impedance_wake(a=0,b=0,c=0,d=0,plane='z',var=freq,func=Zl));

    Zt=transverse_imp_striplineBPM_Ng(l,angle,b,freq,Zc=50);
    if plane=='x':
        imp_mod.append(impedance_wake(a=1,b=0,c=0,d=0,plane='x',var=freq,func=Zt));
    else:
        imp_mod.append(impedance_wake(a=0,b=1,c=0,d=0,plane='y',var=freq,func=Zt));

    if wake_calc:
        Wt=transverse_wake_striplineBPM_Ng(l,angle,b,z,Zc=50);
        if plane=='x':
            wake_mod.append(impedance_wake(a=1,b=0,c=0,d=0,plane='x',var=z,func=Wt));
        else:
            wake_mod.append(impedance_wake(a=0,b=1,c=0,d=0,plane='y',var=z,func=Wt));

    return imp_mod,wake_mod;


def imp_model_4striplinesBPM_Ng(l,angle,b,Zc=50.,beta=1,wake_calc=False,
        fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param()):

    '''compute impedance and wake model from Ng formula as quoted in
    Handbook of Accelerator Physics and Engr., sec 3.2 (Impedances and wakes functions)
    (transverse and longitudinal)

    Note that the longitudinal wake is zero because it's expressed as delta function
    in Ng's formulas (to be checked).

    - l = length of strips, angle = azimuthal angle over which each of them is seen,
    - b = pipe radius (or half distance between strips),
    - Zc = characterisitc impedance of transmission lines from strips (usually
    50 Ohm)
    - fpar / zpar:frequencies / distances over which to compute impedance /wake
    NB: valid for two pairs of striplines in x/y directions '''

    imp_mod=[];wake_mod=[];

    freq=freqscan_from_fpar(fpar);
    z=zscan_from_zpar(zpar);

    Zl=longitudinal_imp_4striplinesBPM_Ng(l,angle,freq,Zc=50);
    imp_mod.append(impedance_wake(a=0,b=0,c=0,d=0,plane='z',var=freq,func=Zl));

    Zt=transverse_imp_4striplinesBPM_Ng(l,angle,b,freq,Zc=50);

    imp_mod.append(impedance_wake(a=1,b=0,c=0,d=0,plane='x',var=freq,func=Zt));
    imp_mod.append(impedance_wake(a=0,b=1,c=0,d=0,plane='y',var=freq,func=Zt));

    if wake_calc:
        Wt=transverse_wake_striplinesBPM_Ng(l,angle,b,z,Zc=50);
        wake_mod.append(impedance_wake(a=1,b=0,c=0,d=0,plane='x',var=z,func=Wt));
        wake_mod.append(impedance_wake(a=0,b=1,c=0,d=0,plane='y',var=z,func=Wt));

    return imp_mod,wake_mod;


def imp_model_from_file(filename,compname,ignored_rows=0,sign=1):
    '''Create an impedance or wake model from a file containing an impedance
    or wake function (typically from a simulation).
    File has 'ignored_rows' header lines.
    There should be one column with frequency [Hz] or distance [m] (first column)
    then one or two columns with real and imaginary part of impedance or wake,
    for the component given in 'compname' (component identification should respect
    the standard given in function "identify_component").
    Type (wake or impedance) is detected from first letter (W or Z) of 'compname'
    (but actually it plays no role here).
    sign is the sign convention: 1 if keep the same sign as in the file,
    -1 otherwise. (use -1 for GdFidl wakes for instance).'''

    a,b,c,d,plane,wakeflag=identify_component(compname);
    imp_mod=[];

    s=read_ncol_file(filename,ignored_rows=ignored_rows);

    func=np.zeros((len(s[:,0]),2));
    if (len(s[0,:])>=3):
        func=sign*s[:,1:3];
    elif (len(s[0,:])==2):
        func[:,0]=sign*s[:,1];
    else:
        print("imp_model_from_file: not enough columns in file "+filename);
        sys.exit();

    imp_mod.append(impedance_wake(a=a,b=b,c=c,d=d,plane=plane,var=s[:,0],func=func));

    return imp_mod;


def imp_model_from_files(filenamelist,scan,value,compname,ignored_rows=0,sign=1):
    '''Create an impedance or wake model from several files, each containing
    an impedance or wake function (typically from simulations).
    Each file has 'ignored_rows' header lines, and should have one column
    with frequency [Hz] or distance [m] (first column)
    then one or two columns with real and imaginary part of impedance or wake,
    for the component given in 'compname' (component identification should respect
    the standard given in function "identify_component").
    The different files represent a scan given by 'scan' (e.g. a half-gap scan)
    and we interpolate between them at the value 'value' (extrapolation is
    possible: we use the closest value in scan).
    Type (wake or impedance) is detected from first letter (W or Z) of 'compname'
    (but actually it plays no role here).
    sign is the sign convention: 1 if keep the same sign as in the files,
    -1 otherwise. (use -1 for GdFidl wakes for instance).'''

    from scipy import interpolate as itp

    a,b,c,d,plane,wakeflag=identify_component(compname);

    # build the total scan interpolation table
    for ifile,filename in enumerate(filenamelist):

        # read file
        s=read_ncol_file(filename,ignored_rows=ignored_rows);

        # choose first file to define the distance or frequency scan
        if ifile==0: var=s[:,0];inttable=np.zeros((len(filenamelist),len(var),2));

        func=np.zeros((len(s[:,0]),2));
        if (len(s[0,:])>=3):
            func=sign*s[:,1:3];
        elif (len(s[0,:])==2):
            func[:,0]=sign*s[:,1];
        else:
            print("imp_model_from_file: not enough columns in file "+filename);
            sys.exit();

        # construct final interpolation table
        inttable[ifile,:,0]=np.interp(var,s[:,0],func[:,0]);
        inttable[ifile,:,1]=np.interp(var,s[:,0],func[:,1]);

    # interpolated function
    fint=itp.interp1d(scan,inttable,axis=0);

    # first treat extrapolation cases
    if (value<scan[0]): funcfinal=inttable[0,:,:].squeeze();
    elif (value>scan[-1]): funcfinal=inttable[-1,:,:].squeeze();
    else: # general case
        funcfinal=fint(value);

    imp_mod=[];
    imp_mod.append(impedance_wake(a=a,b=b,c=c,d=d,plane=plane,var=var,func=funcfinal));

    return imp_mod;


def power_spectrum(f,sigz,powerspectrum='gaussian'):

    '''compute the power-spectrum at a frequency f (which can be an array)
    sigz is the RMS bunch length in m (total length ~ 4*sigz)
    powerspectrum=type of line distribution:
           'gaussian' or by default (if not present): Gaussian
              line distribution (in number of particles):
               lambda(omega)=Nb*exp(-omega^2*sigz^2/(2*c^2))
               or in z domain (in m^-1):
               lambda(z)=Nb/(sqrt(2*pi)*sigz)*exp(-z^2/(2*sigz^2))
           'parabolic': parabolic line distribution
           2-column table of numbers: any power spectrum defined by
           amplitude^2 (second column) vs frequency (1st column) (note: THIS
           IS NOT THE ANGULAR FREQUENCY)'''

    c=299792458; # speed of light

    if isinstance(powerspectrum,np.ndarray)and(len(powerspectrum[0,:])==2):
        # customized bunch spectrum
        spec=np.interp(f,powerspectrum[:,0],powerspectrum[:,1]);

    elif powerspectrum.startswith('parabolic'):
        # parabolic line bunch spectrum
        u=4.*np.pi*f*sigz/c;
        spec=9.*(np.sin(u)-u*np.cos(u))**2/u**6;

    elif powerspectrum.startswith('gaussian'):
        # Gaussian bunch spectrum
        spec=np.exp(-(2.*np.pi*f)**2*sigz**2/(c**2));

    else:
        print("power_spectrum: type of powerspectrum",powerspectrum,"not recognized!");
        spec=0;

    return spec;


def power_loss(imp_mod,sigz,gamma,Nb,M,circum,powerspectrum='gaussian',particle='proton'):

    '''compute the power loss from the longitudinal impedance present in the impedance
    model "imp_mod".
     - sigz is the RMS bunch length in m (total length ~ 4*sigz),
     - gamma: relativistic mass factor,
     - Nb is the number of particles per bunch
     - M is the number of bunches (assumed to be equidistant)
     - circum is the circumference of the machine
     - powerspectrum: type of line distribution:
           'gaussian' or by default (if not present): Gaussian
               line distribution (in number of particles):
               lambda(omega)=Nb*exp(-omega^2*sigz^2/(2*c^2))
               or in z domain (in m^-1):
               lambda(z)=Nb/(sqrt(2*pi)*sigz)*exp(-z^2/(2*sigz^2))
           'parabolic': parabolic line distribution
           2-column table of numbers: any power spectrum defined by
           amplitude^2 (second column) vs frequency (1st column) (note: THIS
           IS NOT THE ANGULAR FREQUENCY)

    see Giovanni Rumolo's USPAS 2009 course : Collective Effects in the
    Longitudinal Plane'''

    e,m0,c,E0=eval(particle+'_param()');
    beta=np.sqrt(1.-1./(gamma**2))
    f0=beta*c/circum; # revolution frequency
    omega0=2*np.pi*f0; # revolution angular frequency

    # find the longitudinal component in imp_mod
    k_long=-1;
    for kiw,iw in enumerate(imp_mod):
        if test_impedance_wake_comp(iw,0,0,0,0,'z'): k_long=kiw;

    if (k_long==-1): return 0; # power loss is zero if no longitudinal impedance
    else:
        freq=imp_mod[k_long].var;
        Zlong=imp_mod[k_long].func[:,0]; # real part of longitudinal impedance

        # determines at which p the line density will begin to be small
        pmax=np.ceil(beta*c/(sigz*omega0));
        pend=np.ceil(pmax/10.);
        s=1; # partial sum
        Ploss=0.;
        p=np.arange(1,pend+1);

        while (p[-1]<2.*pmax)or(np.abs(s/Ploss)>1.e-7):
            # real impedance at p*omega0
            if (p[-1]*M*f0>freq[-1]):
                print('power_loss: warning: interpolating above the upper limit of Zlong in frequency')
            Z=np.interp(p*M*f0,freq,Zlong);
            # weight by bunch power spectrum and sum
            spec=power_spectrum(p*M*f0,sigz,powerspectrum=powerspectrum);
            s=np.sum(spec*Z);

            # add to the sum
            Ploss += s;
            # increment
            p += pend;

        Ploss*=-2.*M**2*Nb**2*e**2*f0**2; # the factor 2 is for the symmetric sum (from p=-inf to p=-1)

        return Ploss;


def transverse_kick_factor(imp_mod,sigz,powerspectrum='gaussian',compname='Zxdip'):

    '''compute transverse kick factor for a certain RMS bunch length sigz and a certain powerspectrum:
    gaussian, parabolic or user defined by an array with 2 columns: power spectrum vs frequency
    (see function power_spectrum above)
    see e.g. "Effect of transverse impedance on luminosity measurement", talk by E. Metral at the ICE meeting,
    slide 3 (03/07/2013)

    It computes the kick factor only for the component 'compname' (identify_component format)'''

    # find the correct component in imp_mod
    a,b,c,d,plane,wakeflag=identify_component(compname);
    k_comp=-1;
    for kiw,iw in enumerate(imp_mod):
        if test_impedance_wake_comp(iw,a,b,c,d,plane): k_comp=kiw;

    if (k_comp==-1): return 0; # power loss is zero if no longitudinal impedance
    else:
        freq=imp_mod[k_comp].var;
        Zimag=imp_mod[k_comp].func[:,1]; # imag. part of impedance

        omega=2.*np.pi*freq;

        # bunch power spectrum
        spec=power_spectrum(freq,sigz,powerspectrum=powerspectrum);

        # integrate over omega (trapz method)
        kickfact=-np.trapz(spec*Zimag,omega)/np.pi;

        return kickfact;


def hmm(m,omega,taub,modetype='sinusoidal'):

    '''compute hmm power spectrum of Sacherer formula, for azimuthal mode number m,
    at angular frequency 'omega' (rad/s) (can be an arrray), for total bunch length
    'taub' (s), and for a kind of mode specified by 'modetype'
    (which can be 'Hermite' - leptons -  or 'sinusoidal' - protons).'''

    if (modetype.startswith('sinus'))or(modetype.startswith('Sinus')):
        # best for protons
        hmm=( (taub*(np.abs(m)+1.))**2/(2.*np.pi**4) )* ( 1.+(-1)**m * np.cos(omega*taub) ) / ( ( (omega*taub/np.pi)**2 - (np.abs(m)+1.)**2)**2);

    elif (modetype=='Hermite')or(modetype=='hermite'):
        # best for leptons
        hmm=(omega*taub/4)**(2*m) * np.exp(-(omega*taub/4.)**2);

    else:
        print("Pb in hmm: : kind of mode not recognized!");sys.exit();

    return hmm;


def hmmsum(m,omega0,M,offk,taub,omegaksi,eps=1e-5,omegas=0.,kmax=20,modetype='sinusoidal',Z=None,omegaZ=None,flagtrapz=False):

    '''compute sum of hmm functions (defined above), weighted or not by the impedance Z
    (table of complex impedances in Ohm given from negative to positive angular frequencies
    omegaZ in rad/s] If these are None, then only sum the hmm.
    Use the trapz integration method if flagtrapz==True.

     - m: azimuthal mode number,
     - omega0: angular revolution frequency in rad/s,
     - M: number of bunches,
     - offk: offset in k (typically nx+[Q] where nx is the coupled-bunch mode and [Q]
     the fractional part of the tune),
     - taub: total bunch length in s,
     - omegaksi: chromatic angular frequency,
     - eps: relative precision of the sum,
     - omegas: synchrotron frequency,
     - kmax: step in k between sums,
     - modetype: kind of mode for the hmm power spectrum ('sinusoidal', 'Hermite').

    In the end the sum runs over k with hmm taken at the angular frequencies
    (offk+k*M)*omega0+m*omegas-omegaksi
    but the impedance is taken at (offk+k*M)*omega0+m*omegas'''

    Momega0=M*omega0;

    if (np.any(omegaZ)==None):
        omega=np.arange(-100.01/taub,100.01/taub,0.01/taub);
        #pylab.plot(omega,hmm(m,omega,taub,modetype=modetype));pylab.show();
        #pylab.loglog(omega,hmm(m,omega,taub,modetype=modetype));pylab.show();
    else:
        omega=omegaZ;

    # sum initialization
    omegak=offk*omega0+m*omegas;
    #omegak=Qfrac*omega0+m*omegas;
    hmm_k=hmm(m,omegak-omegaksi,taub,modetype=modetype) ;

    if flagtrapz:
        # initialization of correcting term sum_i (with an integral instead of discrete sum)
        ind_i=np.where(np.sign(omega-omegak-Momega0)*np.sign(omega-1e15)==-1);
        ind_mi=np.where(np.sign(omega-omegak+Momega0)*np.sign(omega+1e15)==-1);
        omega_i=omega[ind_i];
        omega_mi=omega[ind_mi];
        hmm_i=hmm(m,omega_i-omegaksi,taub,modetype=modetype);
        hmm_mi=hmm(m,omega_mi-omegaksi,taub,modetype=modetype);
        if (Z is not None):
            Z_i=Z[ind_i];Z_mi=Z[ind_mi];
            sum_i=(np.trapz(Z_i*hmm_i,omega_i)+np.trapz(Z_mi*hmm_mi,omega_mi))/(Momega0);
        else:
            sum_i=(np.trapz(hmm_i,omega_i)+np.trapz(hmm_mi,omega_mi))/(Momega0);
    else:
        sum_i=0.;

    if (np.any(Z)!=None):
        Zpk=np.interp(omegak,omega,np.real(Z))+1j*np.interp(omegak,omega,np.imag(Z));
        sum1=Zpk*hmm_k+sum_i;
    else:
        sum1=hmm_k+sum_i;

    k=np.arange(1,kmax+1);oldsum1=10.*sum1;

    while ( (np.abs(np.real(sum1-oldsum1)))>eps*np.abs(np.real(sum1)) )or( (np.abs(np.imag(sum1-oldsum1)))>eps*np.abs(np.imag(sum1)) ):
        oldsum1=sum1;
        # omega_k^x and omega_-k^x in Elias's slides:
        omegak=(offk+k*M)*omega0+m*omegas;
        omegamk=(offk-k*M)*omega0+m*omegas;
        # power spectrum function h(m,m) for k and -k:
        hmm_k=hmm(m,omegak-omegaksi,taub,modetype=modetype) ;
        hmm_mk=hmm(m,omegamk-omegaksi,taub,modetype=modetype) ;

        if flagtrapz:
            # subtract correction (rest of the sum considered as integral -> should suppress redundant terms)
            ind_i=np.where(np.sign(omega-omegak[0])*np.sign(omega-omegak[-1]-Momega0)==-1);
            ind_mi=np.where(np.sign(omega-omegamk[0])*np.sign(omega-omegamk[-1]+Momega0)==-1);
            omega_i=omega[ind_i];
            omega_mi=omega[ind_mi];
            hmm_i=hmm(m,omega_i-omegaksi,taub,modetype=modetype);
            hmm_mi=hmm(m,omega_mi-omegaksi,taub,modetype=modetype);
            if (Z is not None):
                Z_i=Z[ind_i];Z_mi=Z[ind_mi];
                sum_i=(np.trapz(Z_i*hmm_i,omega_i)+np.trapz(Z_mi*hmm_mi,omega_mi))/(Momega0);
            else:
                sum_i=(np.trapz(hmm_i,omega_i)+np.trapz(hmm_mi,omega_mi))/(Momega0);
        else:
            sum_i=0.;

        if (np.any(Z)!=None):
            # impedances at omegak and omegamk
            Zpk=np.interp(omegak,omegaZ,np.real(Z))+1j*np.interp(omegak,omegaZ,np.imag(Z));
            Zpmk=np.interp(omegamk,omegaZ,np.real(Z))+1j*np.interp(omegamk,omegaZ,np.imag(Z));
            # sum
            sum1=sum1+np.sum(Zpk*hmm_k)+np.sum(Zpmk*hmm_mk)-sum_i;
        else:
            # sum
            sum1=sum1+np.sum(hmm_k)+np.sum(hmm_mk)-sum_i;
        k=k+kmax;

    #print k[-1],kmax,omegak[-1],omegaksi,m*omegas;

    return sum1;


def Zl_eff(imp_mod,machine,l,normalized_over_n = 0, revert_to_Z = 1, modeltype = 'Gaussian'):
    '''Computes the longitudinal effective Zl/n for a given impedance model and a given azimuthal mode 'l'.
    Zleff=Zl_eff(imp_mod,machine,l):
    - imp_mod is the impedance model, of which the 'Zlong' component will be used
    - machine is the class defined for the machine: the code relies in the definition of machine.taub (4*rms bunch length in seconds), the sigma_z will be recalculated.
    - l is the azimuthal mode of Hermite (or Gaussian) modes.

    Returns Zleff complex number.
    '''

    import scipy.constants
    clight=scipy.constants.c

    Qs=machine.Qs
    taub=machine.taub
    f0=machine.f0
    sigma_z=machine.taub*clight/4
    compname='Zlong'
    nperc=1e-16 # tolerance on spectrum




    # find the correct component in imp_mod
    a,b,c,d,plane,wakeflag=identify_component(compname);
    k_comp=-1;
    for kiw,iw in enumerate(imp_mod):
        if test_impedance_wake_comp(iw,a,b,c,d,plane): k_comp=kiw;

    if (k_comp==-1): print("Pb in sacherer: component",compname,"not found!");sys.exit();

    freq=imp_mod[k_comp].var;
    Zreal=imp_mod[k_comp].func[:,0]; # real part of impedance
    Zimag=imp_mod[k_comp].func[:,1]; # imag. part of impedance

    omegap=2.*np.pi*freq;omegam=-omegap[::-1];

    if normalized_over_n == 1:
        Zpcomp=Zreal+1j*Zimag;Zmcomp=-Zreal[::-1]+1j*Zimag[::-1]; # compute complex impedance and 'symmetrize' it for negative frequencies
    elif (normalized_over_n == 0) and (revert_to_Z == 1):
        Zreal = Zreal*freq/f0
        Zimag = Zimag * freq/f0
        Zpcomp=Zreal+1j*Zimag;Zmcomp=Zreal[::-1]-1j*Zimag[::-1]; # compute complex impedance and 'symmetrize' it for negative frequencies
    elif (normalized_over_n == 0) and (revert_to_Z == 0):
        Zreal = Zreal/(freq/f0)
        Zimag = Zimag /( freq/f0)

        Zpcomp=Zreal+1j*Zimag;Zmcomp=-Zreal[::-1]+1j*Zimag[::-1];

    omega=np.concatenate((omegam,omegap));
    Zcomp=np.concatenate((Zmcomp,Zpcomp));


    omega_0=2*np.pi*f0;
    omega_s=Qs*omega_0;



    if modeltype=='Gaussian':

        h=(omega*sigma_z/clight)**(2*l)*np.exp(-omega**2*sigma_z**2/clight**2);

        omega_part=omega[omega>0];
        h_part=h[omega>0];
        omega_max=omega_part[h_part.argmax()];
        h_part2=h_part[omega_part>omega_max];
        omega_part2=omega_part[omega_part>omega_max];
        h_part2,ind=np.unique(h_part2,return_index=True); # delete trailed zeros;
        omega_part2=omega_part2[ind];
        omega_extr=np.interp(nperc*h_part.max(),h_part2,omega_part2);

        Pmax=np.floor(omega_extr/omega_0);
        print('max number of line index Pmax=%d'%Pmax) # approx value for max intergation;
        p=np.arange(-Pmax,Pmax)
        p=p[p!=0]

        omega_p=p*omega_0+l*omega_s;

        omega_p=np.unique(omega_p);
        h_p=(omega_p*sigma_z/clight)**(2*l)*np.exp(-omega_p**2*sigma_z**2/clight**2);
        Z_p=np.interp(omega_p,omega,Zcomp.real)+1j*np.interp(omega_p,omega,Zcomp.imag);
        N=Z_p*h_p
        D=h_p;

    elif modeltype=='Parabolic':

        print(machine.taub)
        tmax = machine.taub/2
        h = (3*(np.sin(omega*tmax)-omega*tmax*np.cos(omega*tmax))/(omega*tmax)**3)**2

        omega_part=omega[omega>0];
        h_part=h[omega>0];
        omega_max=omega_part[h_part.argmax()];
        h_part2=h_part[omega_part>omega_max];
        omega_part2=omega_part[omega_part>omega_max];
        h_part2,ind=np.unique(h_part2,return_index=True); # delete trailed zeros;
        omega_part2=omega_part2[ind];
        omega_extr=np.interp(nperc*h_part.max(),h_part2,omega_part2);

        Pmax=np.floor(omega_extr/omega_0);
        print('max number of line index Pmax=%d'%Pmax) # approx value for max intergation;
        p=np.arange(-Pmax,Pmax)
        p=p[p!=0]

        omega_p=p*omega_0+l*omega_s;

        omega_p=np.unique(omega_p);
        h_p=(3*(np.sin(omega_p*tmax)-omega_p*tmax*np.cos(omega_p*tmax))/(omega_p*tmax)**3)**2
        Z_p=np.interp(omega_p,omega,Zcomp.real)+1j*np.interp(omega_p,omega,Zcomp.imag);
        N=Z_p*h_p
        D=h_p;

    Zleff=np.sum(N)/np.sum(D)
    kl = np.sum(N)* machine.f0

    return Zleff, kl, omega_p, Z_p, h_p


def sacherer(imp_mod,Qpscan,nxscan,Nbscan,omegasscan,M,omega0,Q,gamma,eta,taub,mmax,
        particle='proton',modetype='sinusoidal',compname='Zxdip',flagtrapz=None):

    '''computes frequency shift and effective impedance from Sacherer formula, in transverse, in the case of low
    intensity perturbations (no mode coupling), for modes of kind 'modetype'.
    It gives in output:
     - tuneshift_most: tune shifts for the most unstable multibunch mode and synchrotron modes
    sorted by ascending imaginary parts (most unstable synchrotron mode first).
    Array of dimensions len(Qpscan)*len(Nbscan)*len(omegasscan)*(2*mmax+1)
     - tuneshiftnx: tune shifts for all multibunch modes and synchrotron modes m.
    Array of dimensions len(Qpscan)*len(nxscan)*len(Nbscan)*len(omegasscan)*(2*mmax+1)
     - tuneshiftm0: tune shifts for the most unstable multibunch mode and synchrotron mode m=0.
    Array of dimensions len(Qpscan)*len(Nbscan)*len(omegasscan)
     - Zeff: effective impedance for different multibunch modes and synchrotron modes m.
    Array of dimensions len(Qpscan)*len(nxscan)*len(omegasscan)*(2*mmax+1)

    Input parameters are similar to DELPHI's ones:
     - imp_mod: impedance model (list of impedance-wake objects),
     - Qpscan: scan in Q' (DeltaQ*p/Deltap),
     - nxscan: scan in multibunch modes (from 0 to M-1),
     - Nbscan: scan in number of particles per bunch,
     - omegasscan: scan in synchrotron angular frequency (Qs*omega0),
     - M: number of bunches,
     - omega0: angular revolution frequency
     - Q: transverse betatron tune (integer part + fractional part),
     - gamma: relativistic mass factor,
     - eta: slip factor (Elias's convention, i.e. oppostie to Joel Le Duff),
     - taub: total bunch length in seconds,
     - mmax: azimuthal modes considered are from -mmax to mmax,
     - particle: 'proton' or 'electron',
     - modetype: 'sinusoidal' or 'Hermite': kind of modes in effective impedance,'
     - compname: component to extract from impedance model

     see Elias Metral's USPAS 2009 course : Bunched beams transverse coherent
     instabilities.
     
     NOTE: this is NOT the original Sacherer formula, which assumes an impedance normalized by beta
     (see E. Metral, USPAS 2009 lectures, or C. Zannini, 
     https://indico.cern.ch/event/766028/contributions/3179810/attachments/1737652/2811046/Z_definition.pptx)
     Here this formula is instead divided by beta (compared to Sacherer initial one),
     so is valid with our usual definition of impedance (not beta-normalized).
     This was corrected on April 15th, 2019. NM
     '''


    e,m0,c,E0=eval(particle+'_param()');

    # some parameters
    Z0=4.e-7*np.pi*c; # free space impedance: here mu0 c (SI unit - Ohm) or 4 pi/c (c.g.s)
    beta=np.sqrt(1.-1./(gamma**2)); # relativistic velocity factor
    f0=omega0/(2.*np.pi); # revolution angular frequency
    Ibscan=e*Nbscan*f0; # single-bunch intensity
    Qfrac=Q-np.floor(Q); # fractional part of the tune
    Lb=taub*beta*c; # full bunch length (in meters)

    # find the correct component in imp_mod
    a,b,c,d,plane,wakeflag=identify_component(compname);
    k_comp=-1;
    for kiw,iw in enumerate(imp_mod):
        if test_impedance_wake_comp(iw,a,b,c,d,plane): k_comp=kiw;

    if (k_comp==-1): print("Pb in sacherer: component",compname,"not found!");sys.exit();

    freq=imp_mod[k_comp].var;
    Zreal=imp_mod[k_comp].func[:,0]; # real part of impedance
    Zimag=imp_mod[k_comp].func[:,1]; # imag. part of impedance

    omegap=2.*np.pi*freq;omegam=-omegap[::-1];
    # compute complex impedance and 'symmetrize' it for negative frequencies
    Zpcomp=Zreal+1j*Zimag;Zmcomp=-Zpcomp[::-1].conjugate();
    omega=np.concatenate((omegam,omegap));Zcomp=np.concatenate((Zmcomp,Zpcomp));

    # first guess of the maximum k, and step for k (compute sums on the array
    # [0 kmax-1]+[Integer]*kmax)

    eps=1.e-5; # relative precision of the summations
    tuneshiftnx=np.zeros((len(Qpscan),len(nxscan),len(Nbscan),len(omegasscan),2*mmax+1),dtype=complex);
    tuneshift_most=np.zeros((len(Qpscan),len(Nbscan),len(omegasscan),2*mmax+1),dtype=complex);
    tuneshiftm0=np.zeros((len(Qpscan),len(Nbscan),len(omegasscan)),dtype=complex);
    Zeff=np.zeros((len(Qpscan),len(nxscan),len(omegasscan),2*mmax+1),dtype=complex);

    for iQp,Qp in enumerate(Qpscan):

        omegaksi=Qp*omega0/eta;
        if flagtrapz is None:
            flagtrapz=(np.ceil(100.*(4.*np.pi/taub+abs(omegaksi))/omega0/M)>1e9);
        #print np.ceil(100.*(4.*np.pi/taub+abs(omegaksi))/omega0/M),"flagtrapz=",flagtrapz

        for inx,nx in enumerate(nxscan): # coupled-bunch modes

            for iomegas,omegas in enumerate(omegasscan):

                for im,m in enumerate(range(-mmax,mmax+1)):

                    # consider each sychrotron mode individually
                    # sum power spectrum functions and computes effective impedance

                    # sum power functions
                    # BE CAREFUL: maybe for this "normalization sum" the sum should run
                    # on all single-bunch harmonics instead of only coupled-bunch
                    # harmonics (and then the frequency shift should be multiplied by
                    # M). This has to be checked.
                    sum1=hmmsum(m,omega0,M,nx+Qfrac,taub,omegaksi,eps=eps,omegas=omegas,kmax=20,
                        modetype=modetype,flagtrapz=flagtrapz);

                    # effective impedance
                    sum2=hmmsum(m,omega0,M,nx+Qfrac,taub,omegaksi,eps=eps,omegas=omegas,kmax=20,
                        modetype=modetype,Z=Zcomp,omegaZ=omega,flagtrapz=flagtrapz);

                    Zeff[iQp,inx,iomegas,im]=sum2/sum1;
                    freqshift =( 1j*e*Ibscan/(2.*(np.abs(m)+1.)*m0*gamma*Q*omega0*Lb) ) * sum2/sum1; #15/04/2019 NM: beta suppressed (was for a "beta-normalized" definition of impedance)
                    tuneshiftnx[iQp,inx,:,iomegas,im]=freqshift/omega0+m*omegas/omega0;

        # find the most unstable coupled-bunch mode
        for iomegas,omegas in enumerate(omegasscan):

            for im,m in enumerate(range(-mmax,mmax+1)):

                inx=np.argmin(np.imag(tuneshiftnx[iQp,:,-1,iomegas,im])); # check one intensity (last one) is enough
                #print "Sacherer: Qp=",Qp,", M=",M,", omegas=",omegas,", m=",m,", Most unstable coupled-bunch mode: ",nxscan[inx];
                tuneshift_most[iQp,:,iomegas,im]=tuneshiftnx[iQp,inx,:,iomegas,im];
                if (m==0): tuneshiftm0[iQp,:,iomegas]=tuneshiftnx[iQp,inx,:,iomegas,im];

            # sort tuneshift_most (most unstable modes first) (only one instensity - last one - is enough)
            ind=np.argmin(np.imag(tuneshift_most[iQp,-1,iomegas,:]));
            for iNb,Nb in enumerate(Nbscan): tuneshift_most[iQp,iNb,iomegas,:]=tuneshift_most[iQp,iNb,iomegas,ind];

    return tuneshift_most,tuneshiftnx,tuneshiftm0,Zeff;


def write_wake_HEADTAIL(wake_mod,filename,beta=1,ncomp=6):

    '''write a wake file, in the HEADTAIL format, from the wake model 'wake_mod'
    reminder: units are ns and kV/pC(/m)
    beta: relativistic beta factor (for unit conversion - m to ns)
    ncomp= number of components to take (correspond to flag i_waketb in HEADTAIL input file)
    NOTE: HEADTAIL assumes xy=yx for coupling terms; we also make this assumption here'''

    c=299792458; # speed of light

    listcomp=['Wlong','Wxdip','Wydip','Wxquad','Wyquad','Wxydip','Wxyquad','Wxcst','Wycst']
    lista=[0,1,0,0,0,0,0,0,0];listb=[0,0,1,0,0,1,0,0,0];
    listc=[0,0,0,1,0,0,0,0,0];listd=[0,0,0,0,1,0,1,0,0];
    listplane=['z','x','y','x','y','x','x','x','y'];

    # components to select (depends on ncomp=i_waketb; if its odd Wlong is in the list)
    ind=np.arange(1-np.mod(ncomp,2),1-np.mod(ncomp,2)+ncomp);

    k_comp=[];
    for icomp in ind:
        flagfound=False;
        for kw,w in enumerate(wake_mod):
            if test_impedance_wake_comp(w,lista[icomp],listb[icomp],listc[icomp],listd[icomp],listplane[icomp]): k_comp.append(kw);flagfound=True;
        if (not flagfound): print("Pb in write_wake_HEADTAIL: component "+listcomp[icomp]+" not found !");sys.exit();

    # collect all distances z into one array
    z=[];
    for k in k_comp: z=sort_and_delete_duplicates(np.concatenate((z,wake_mod[k].var)));

    # unit conversion and build the final table
    data=z.reshape((-1,1))*1e9/(beta*c);
    for k in k_comp:
        s=1e-15*np.interp(z,wake_mod[k].var,wake_mod[k].func[:,0]);
        data=np.hstack((data,s.reshape((-1,1))));

    # write in file
    write_ncol_file(filename,data);

    return;


def write_imp_wake_mod(imp_mod,name,listcomp=['Zxdip','Zydip'],dire='',header_flag=True):

    '''write files with components of an impedance or wake model given in imp_mod.
    filenames will be: [component name][name].dat
    listcomp= components to be written (each in a separate file, with 3 columns:
    frequency - or distance, real and imaginary parts)
    dire is the directory where to put the files (do not forget the "/")'''

    units=["","/m","/m^2","/m^3","/m^4"];

    for icomp,comp in enumerate(listcomp):

        a,b,c,d,plane,wakeflag=identify_component(comp);
        unit=units[a+b+c+d];

        if header_flag:
            if wakeflag: header="Distance[m] Re_"+comp+"[V/C"+unit+"] Im_"+comp+"[V/C"+unit+"]";
            else: header="Frequency[Hz] Re_"+comp+"[Ohm"+unit+"] Im_"+comp+"[Ohm"+unit+"]";
        else:
            header=None;

        for iw in imp_mod:
            flag=False;
            if test_impedance_wake_comp(iw,a,b,c,d,plane): Z=deepcopy(iw.func);freq=deepcopy(iw.var);flag=True;
            if flag:
                data=np.hstack((freq.reshape((-1,1)),Z.reshape((-1,2))));
                write_ncol_file(dire+comp+name+'.dat',data,header=header);

    return;


def identify_component(compname):
    '''identify a component from its name, and convert it to the impedance/wake class
    format in terms  of a, b, c, d, plane and wakeflag (wakeflag=True if this is a
    wake component).
    first letter of compname should be 'Z' (for impedance) or 'W' (for wake potential)
    for e.g. impedance, compname can be: Zl(ong), Zxd(ip), Zyd(ip), Zxq(uad), Zyq(uad), Zxyd(ip), Zyxd(ip),
    Zxyq(uad), Zyxq(uad), Zxc(st), Zyc(st), or a name directly indicating a, b, c, d
    and the plane, in the form e.g. "Z1000x" (Zxdip in this case)
    (for wakes, it is the same replacing 'Z' by 'W')'''

    if (not compname.startswith('Z'))and(not compname.startswith('W')):
        print("Pb in identify_component: component name does not begin by Z or W");
        sys.exit();

    else:
        first=compname[0]; # Z or W
        if first=='W': wakeflag=True;
        else: wakeflag=False;

        if compname.startswith(first+'l'): a=0;b=0;c=0;d=0;plane='z';
        elif compname.startswith(first+'xd'): a=1;b=0;c=0;d=0;plane='x';
        elif compname.startswith(first+'yd'): a=0;b=1;c=0;d=0;plane='y';
        elif compname.startswith(first+'xq'): a=0;b=0;c=1;d=0;plane='x';
        elif compname.startswith(first+'yq'): a=0;b=0;c=0;d=1;plane='y';
        elif compname.startswith(first+'xyd'): a=0;b=1;c=0;d=0;plane='x';
        elif compname.startswith(first+'yxd'): a=1;b=0;c=0;d=0;plane='y';
        elif compname.startswith(first+'xyq'): a=0;b=0;c=0;d=1;plane='x';
        elif compname.startswith(first+'yxq'): a=0;b=0;c=1;d=0;plane='y';
        elif compname.startswith(first+'xc'): a=0;b=0;c=0;d=0;plane='x';
        elif compname.startswith(first+'yc'): a=0;b=0;c=0;d=0;plane='y';

        else:
            try:
                a=int(compname[1]);b=int(compname[2]);
                c=int(compname[3]);d=int(compname[4]);
                plane=compname[5];
            except IndexError:
                print("Pb in identify_component: bad length for component name");
                sys.exit();
            except ValueError:
                print("Pb in identify_component: bad format for component name");
                sys.exit();

    return a,b,c,d,plane,wakeflag;


def plot_compare_imp_model(imp_mod_list,leglist,listcomp=['Zlong','Zxdip','Zydip'],
        figimp=None,aximp=None,figratio=None,axratio=None,saveimp='',saveratio='',
        xlim=[1e3,1e11],ylim=[1e5,1e9],ylim_ratio=None,yliml=[1,1e6],bounds=[20e6,2e9],beta=1.,
        legpos=0,plotpercent=False,legpercentpos=(0.8,0.6),maxpercent=110.,markimp=None):

    '''plot on the same plot various impedance models provided in imp_mod_list
    (list of impedance models, each of them being itself a list of components).
    comparison is done for a set of components given in listcomp (see format in
    function 'identify_component' just above).
    The impedance ratio is also plotted (w.r.t. to the first model in the list).
    - if figimp /aximp (resp. figratio / axratio) are provided and not None,
    imp. plots (resp. imp. ratio plots) are made on the corresponding axes
    (if they are lists of the length of listcomp, then each component is plotted
    on a separate plot).
    otherwise figures and axes are created here and all components are plotted
    separately.
    - leglist is the list of legends (same length as imp_mod_list)
    - if saveimp (resp. saveratio) is a string of non-zero length, save plot
    on this filename.
    - xlim, ylim, ylim_ratio and yliml are the axes limits (only for impedance plot for ylim,
    only for ratio plot for ylim_ratio - if None use [0,maximum ratio] - and only
    for longitudinal impedance for yliml)
    - bounds are the frequencies between which we compute the maximum ratio;
    they are converted to distances in the case of a wake (beta is used in
    that case, for the conversion)
    - If plotpercent is True, we plot (on a single plot per component) the percentages
    w.r.t to the first model given, in a "filled" manner, instead of curves (i.e. the
    plot is filled betweeen each individual curves),
    - legpercentpos: position (bbox_to_anchor) of legend for percentage plot.
    - maxpercent: maximum of y-axis for the percentage plot (in percent).
    - markimp: if not None, list of markers  (like 'x') for each curve of impedance plot
    (not used for ratio plot). Default is no marker.

    It also works for wakes, with component names beginning by "W" instead of "Z".
    Then xlim and ylim have to be changed accordingly (e.g. resp. [1e-5,1e6] and [1e8,1e19] for LHC total wakes)'''

    from matplotlib.patches import Rectangle
    from matplotlib.font_manager import FontProperties
    from plot_lib import init_figure,end_figure,build_colors,plot,fillplot_percent

    clight=299792458;

    # if figure and axes are not None, check if they are lists, if not replace them
    # by a list with the same length as listcomp (repeating the same element)
    for name in ['imp','ratio']:
        if (eval('fig'+name+'!=None'))or(eval('ax'+name+'!=None')):
            exec('fig'+name+'=create_list_for_figure(fig'+name+',n=len(listcomp))');
            exec('ax'+name+'=create_list_for_figure(ax'+name+',n=len(listcomp))');
        else:
            for b in ['fig','ax']: exec(b+name+'=[]');
            for i in range(len(listcomp)):
                for ir in range(plotpercent+1): # do twice more figures (one for real, one for imag.) if "filled percentage plot"
                    if (len(imp_mod_list)<=6): fig,ax=init_figure();
                    else:
                        fig,ax=init_figure(axes=[0.08,0.12,0.55,0.8],figsize=(20,12));
                        legpos=(1.05,-0.05);legpercentpos=(1.65,0.8);
                    for b in ['fig','ax']:eval(b+name+'.append('+b+')');

    pat=['-','--'];units=["","/m","/m$^2$","/m$^3$","/m$^4$"];
    #col=['b','r','g','m','k','c','y']; # colors
    col=build_colors(len(imp_mod_list),randomize=True);
    if markimp==None: markimp=create_list('',n=len(imp_mod_list));
    maxratio=np.zeros((len(listcomp),len(imp_mod_list)-1,2));
    partcomp=['real','imag'];

    for icomp,comp in enumerate(listcomp):
        # identify the kind of component in terms of a,b,c,d and plane
        a,b,c,d,plane,wakeflag=identify_component(comp);
        unit=units[a+b+c+d];

        if wakeflag:
            parts=[''];str1="wakes";str2="Wake";
            xlab="Distance behind the source [m]";ylab="Wake [V/C";
            boundmin=beta*clight/bounds[1];boundmax=beta*clight/bounds[0];
        else:
            parts=['Real part, ','Imag part, '];str1="impedances";str2="Imp.";
            xlab="Frequency [Hz]";ylab="Z [$\Omega$";
            boundmin=bounds[0];boundmax=bounds[1];

        if plotpercent:
            p=[];leglistless=[];
            for ir,r in enumerate(parts): p.append([]);leglistless.append([]);

        kiw0_comp=-1;
        for imod,imp_mod in enumerate(imp_mod_list):

            kiw_comp=-1;
            for kiw,iw in enumerate(imp_mod):
                if test_impedance_wake_comp(iw,a,b,c,d,plane): kiw_comp=kiw;

            if (kiw_comp>=0):

                if imod==0:
                    kiw0_comp=kiw_comp;
                    if plotpercent: # initialize sums for filled percent plot
                        xsum=imp_mod[kiw0_comp].var;ysum=np.zeros(len(xsum),dtype=complex);

                # impedance plot
                for ir,r in enumerate(parts):
                    plot(imp_mod[kiw_comp].var,np.abs(imp_mod[kiw_comp].func[:,ir]),r+leglist[imod],pat[ir]+markimp[imod],ylab+unit+"]",aximp[icomp],3,xlab=xlab,colr=col[imod]);

                # impedance ratio plot
                if (imod>0)and(kiw0_comp>=0):
                    for ir,r in enumerate(parts):
                        imp0=np.interp(imp_mod[kiw_comp].var,imp_mod_list[0][kiw0_comp].var,imp_mod_list[0][kiw0_comp].func[:,ir]);
                        ratio=imp_mod[kiw_comp].func[:,ir]/imp0;

                        if plotpercent:
                            sumratio=np.interp(imp_mod[kiw_comp].var,xsum,getattr(ysum,partcomp[ir]))/imp0;
                            fillplot_percent(imp_mod[kiw_comp].var,sumratio+ratio,sumratio,xlab,leglist[imod],col[imod],axratio[2*icomp+ir])
                            ysum += (1j)**ir *np.interp(xsum,imp_mod[kiw_comp].var,imp_mod[kiw_comp].func[:,ir])
                            # with fill_between we have to create the legend box ourselves
                            p[ir].append(Rectangle((0, 0), 1, 1,axes=axratio[2*icomp+ir],color=col[imod]));
                            leglistless[ir].append(leglist[imod]);

                        else: plot(imp_mod[kiw_comp].var,np.abs(ratio),r+leglist[imod],pat[ir],str2+" ratio w.r.t "+leglist[0],axratio[icomp],1,xlab=xlab,colr=col[imod]);

                        # find max ratio
                        ind=np.where((imp_mod[kiw_comp].var<boundmax)*(imp_mod[kiw_comp].var>boundmin))
                        maxratio[icomp,imod-1,ir]=np.max(ratio[ind])
                        print((r+comp+": max. ratio from {} kHz to {} GHz, between {} and {} {}: {}".format(
                            round(bounds[0]/1e3,3),round(bounds[1]/1e9,3),leglist[imod],
                            leglist[0],str1,round(maxratio[icomp,imod-1,ir],4))))

        # finish plots
        aximp[icomp].set_xlim(xlim);
        if (comp.find('l')!=-1): aximp[icomp].set_ylim(yliml);
        else: aximp[icomp].set_ylim(ylim);
        end_figure(figimp[icomp],aximp[icomp],save=(len(saveimp)>0)*(saveimp+'_'+comp),legpos=legpos);

        if imod>0:

            if not(plotpercent):
                axratio[icomp].set_xlim(xlim);
                if ylim_ratio==None: axratio[icomp].set_ylim([0,np.ceil(np.max(maxratio[icomp,:,:])*5)/5.]);
                else: axratio[icomp].set_ylim(ylim_ratio);
                end_figure(figratio[icomp],axratio[icomp],save=(len(saveratio)>0)*(saveratio+'_'+comp),legpos=legpos);

            else:
                fontP = FontProperties();
                if len(imp_mod_list)<=6: fontP.set_size('medium');
                else: fontP.set_size('small');
                for ir,r in enumerate(parts):
                    axratio[2*icomp+ir].legend(p[ir][::-1],leglistless[ir][::-1],bbox_to_anchor=legpercentpos,prop=fontP);
                    axratio[2*icomp+ir].set_xscale('log');
                    axratio[2*icomp+ir].set_xlim(xlim);
                    axratio[2*icomp+ir].set_ylim([0,maxpercent]);
                    if len(r)>0: part='_'+r;part=part.replace(', ','').replace(' ','_').replace('_part','')
                    else: part='';
                    end_figure(figratio[2*icomp+ir],axratio[2*icomp+ir],save=(len(saveratio)>0)*(saveratio+'_percent_'+comp+part),legpos=legpercentpos);

    return maxratio;


if __name__ == "__main__":

    write_impedance_wake_input('test.txt',impedance_wake_input())

    layer1=layer(thickness=25e-3);
    layer2=layer(rhoDC=7.2e-7,tau=0,epsb=1,mur=1,fmu=np.inf,thickness=np.inf);
    imp=impedance_wake_input(layers=[layer1,layer2],Yokoya=[1,np.pi**2/24,np.pi**2/12,-np.pi**2/24,np.pi**2/24]);
    write_impedance_wake_input('test.txt',imp)

    imp2=imp;imp2.geometry='flat';
    write_impedance_wake_input('test.txt',imp2)

    layer1=layer(rhoDC=15e-6,tau=1.3e-12,epsb=1,mur=1,fmu=np.inf,thickness=25e-3);
    imp3=impedance_wake_input(gamma=3730.26,length=3,b=[6.5e-3,0],layers=[layer1,layer2],geometry='flat');
    write_impedance_wake_input('test.txt',imp3)

    subprocess.Popen("pwd", shell=True)
    a=subprocess.check_output("pwd", shell=True).strip().decode()
    print(a)
