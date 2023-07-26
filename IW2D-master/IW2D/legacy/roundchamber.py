from typing import Dict

import numpy as np
from .utils import generate_frequency_range, make_1d_array, read_input, read_input_layer
from ..utils import py_complex
from ..constants import *
import shutil

from cppyy.gbl import amp, C, alphaTM


def roundchamber(input_dict: Dict[str, str], input_file: str):

    ######## Get input from input file ###########

    # Get strings for output filename. Note: Does not raise warning or error when missing
    machine         =   input_dict.get("Machine", "")
    commentoutput   =   input_dict.get("Comments for the output files names", "")

    # Get general parameters
    gamma           =   read_input(input_dict, fMP, "Relativistic Gamma", default="479.6")
    beta            =   amp.sqrt(oneMP-oneMP/amp.sqr(gamma))
    L               =   read_input(input_dict, float, "Impedance Length in m", default="1")
    N               =   read_input(input_dict, int, "Number of layers", default="1")
    yokoya          =   read_input(input_dict, str, "Yokoya factors long, xdip, ydip, xquad, yquad", default = "1 1 1 0 0")
    yokoya          =   [float(yokoya) for yokoya in yokoya.split()]

    # Get layer properties
    rho             =   read_input_layer(input_dict, N, fMP, "DC resistivity (Ohm.m)", default="1e5")
    tau             =   np.array(read_input_layer(input_dict, N, fMP, "relaxation time for resistivity (ps)", default="0"))
    tau             *=  1.e-12 # ps -> s
    epsb            =   read_input_layer(input_dict, N, fMP, "real part of dielectric constant",default="1")
    chi             =   read_input_layer(input_dict, N, fMP, "magnetic susceptibility", default="0")
    fmu             =   read_input_layer(input_dict, N, fMP, "relaxation frequency of permeability (MHz)", default="Infinity")
    fmu             *=  1.e6 # MHz -> Hz
    b0              =   read_input(input_dict, fMP, "Layer 1 inner radius in mm", default="2")
    b0              *=  1e-3 # mm -> m
    thick           =   read_input_layer(input_dict, N, fMP, "thickness in mm", default = "Infinity")
    thick           *=  1e-3 # mm -> m
    b               =   np.cumsum(np.concatenate([[b0], thick]))

    # Get frequency range and generate array
    freqs = generate_frequency_range(input_dict)

    # Set up output arrays
    Zlong  = np.zeros_like(freqs, dtype=complex)
    Zxdip  = np.zeros_like(freqs, dtype=complex)
    Zydip  = np.zeros_like(freqs, dtype=complex)
    Zxquad = np.zeros_like(freqs, dtype=complex)
    Zyquad = np.zeros_like(freqs, dtype=complex)

    ######## Loop over frequencies, calculate impedance #########
    for i, freq in enumerate(freqs):
        omega = twopiMP*freq # ampf
        k = omega/(beta*C) # ampf


        eps1 = np.zeros(N+1, dtype=cMP) # campf
        mu1  = np.zeros(N+1, dtype=cMP) # campf

        # Innermost layer always vacuum
        eps1[0] = cMP(1) # = 1 + 0j
        mu1[0] = cMP(1)

        # Calculate layer eps1 and mu1
        for j in range(N):
            if rho[j].isFiniteNumber():
                eps1[j+1] = epsb[j] + oneMP/(jimagMP*eps0*rho[j]*omega*(oneMP + jimagMP*omega*tau[j]))
            else:
                eps1[j+1] = epsb[j]
            mu1[j+1] = oneMP + chi[j]/(oneMP + jimagMP*omega/(twopiMP*fmu[j]))
        
        # Numpy vectorized version of calculation:
        # eps1 = epsb + one_over(rho*jimag*eps0*omega*(tau*jimag*omega + one))
        # mu1  = chi*one_over(one_over(fmu*amp.twopi[precision]())*rho*jimag*eps0*omega + one) + one
        # This works, but was discarded because numpy didn't interface too well with cppyy and amp, so it turned out less readable
        
        # Make ap.template_1d_arrays from the numpy arrays
        eps1_arr = make_1d_array((1,N+1), cMP, eps1)
        mu1_arr = make_1d_array((1,N+1), cMP, mu1)
        b_arr = make_1d_array((1,N+1), fMP, b)

        # Calculate alpha coefficients (in c++)
        alphaTM0=alphaTM(N+1,0,eps1_arr,mu1_arr,b_arr,k,beta) # complex
        alphaTM1=alphaTM(N+1,1,eps1_arr,mu1_arr,b_arr,k,beta) # complex   

        # Compute all impedances from the alphas
        cst = jimagMP*L*(k*Z0/(beta*amp.sqr(gamma)*twopiMP)) # campf
        Zlong[i] = yokoya[0]*alphaTM0*py_complex(cst) # complex
        cst = cst*(k/(amp.sqr(gamma)*twoMP)) # campf
        Zxdip[i] = yokoya[1]*py_complex(cst)*alphaTM1 # complex
        Zydip[i] = yokoya[2]*py_complex(cst)*alphaTM1 # complex
        if yokoya[0] == yokoya[1] == yokoya[2] == 1:
            Zxquad[i] = py_complex(cst)*alphaTM0 # complex
            Zyquad[i] = Zxquad[i] # complex
        else:
            Zxquad[i] = py_complex(cst)*alphaTM1*yokoya[3] # complex
            Zyquad[i] = py_complex(cst)*alphaTM1*yokoya[4] # complex


    ########## Save the results to .dat files ##########
    filename_base = f"W{machine}_{N}layers{1e3*b0.toDouble():.2f}mm{commentoutput}"
    shutil.copy(input_file, f"InputData{filename_base}.dat")

    for name, data, unit in zip(
        ["Zlong","Zxdip", "Zydip", "Zxquad", "Zyquad"],
        [Zlong, Zxdip, Zydip, Zxquad, Zyquad],
        ["Ohm","Ohm/m","Ohm/m","Ohm/m","Ohm/m"]
    ):
        np.savetxt(
            f"{name}{filename_base}.dat",
            np.transpose([freqs, np.real(data), np.imag(data)]),
            fmt="%13.8e",
            header=f"Frequency [Hz]\tRe({name}) [{unit}]\tIm({name}) [{unit}]",
            comments=""
        )
