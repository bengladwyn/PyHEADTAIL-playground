from array import array
from typing import Dict, Tuple

from .utils import make_1d_array, read_input, read_input_layer, generate_frequency_range
from ..constants import *
import numpy as np
from itertools import product
import shutil


from cppyy.gbl import amp, gsl_integration_workspace_alloc, gsl_integration_workspace_free, alphamn, C, get_value_from_interp, memorycontainer, coefmn

# constants
INTEGRATION_WORKSPACE_LIMIT = 1000
MAXMEM = 50000

# Pairs of inputs that are not accepted together
MUTUALLY_EXCLUSIVE_INPUTS = [
    ["frequency dependent conductivity", "DC resistivity (Ohm.m)"],
    ["frequency dependent conductivity", "relaxation time for resistivity (ps)"],
    ["frequency dependent conductivity", "dielectric tan(delta)"],
    ["frequency dependent conductivity", "frequency dependent relative complex permittivity"],
    ["frequency dependent relative complex permittivity", "DC resistivity (Ohm.m)"],
    ["frequency dependent relative complex permittivity", "relaxation time for resistivity (ps)"],
    ["frequency dependent relative complex permittivity", "dielectric tan(delta)"],
    ["frequency dependent relative complex permittivity", "real part of dielectric constant"],
    ["dielectric tan(delta)", "DC resistivity (Ohm.m)"],
    ["dielectric tan(delta)", "relaxation time for resistivity (ps)"],
    ["frequency dependent relative complex permeability", "magnetic susceptibility"],
    ["frequency dependent relative complex permeability", "relaxation frequency of permeability (MHz)"]
]



def flatchamber(input_dict: Dict[str, str], input_file: str):

    ######## Get input from input file ###########

    # Get strings for output filename. Note: Does not raise warning or error when missing
    machine         =   input_dict.get("Machine", "")
    commentoutput   =   input_dict.get("Comments for the output files names", "")

    # Get general parameters
    gamma           =   read_input(input_dict, fMP, "Relativistic Gamma", default="479.6")
    beta            =   amp.sqrt(oneMP-oneMP/amp.sqr(gamma))
    order           =   read_input(input_dict, int, "Maximum order of non-linear terms", default="0")
    L               =   read_input(input_dict, float, "Impedance Length in m", default="1")
    topbotsym_str   =   read_input(input_dict, str, "Top bottom symmetry (yes or no)", default="yes")
    flag_topbotsym  =   (topbotsym_str.lower() in ["yes", "y", "1"])
    N               =   read_input(input_dict, int, "Number of upper layers in the chamber wall", default="2")
    layer_numbers   =   list(range(1,N+1))

    if flag_topbotsym:
        M = N
    else:
        M               =   read_input(input_dict, int, "Number of lower layers in the chamber wall", default="2")
        layer_numbers   =   layer_numbers + list(range(-1, -(M+1), -1))
        
    # Check that no mutually exclusive inputs were given together for the same layer
    for p in layer_numbers:
        layer = f"Layer {p} "
        for entry1, entry2 in MUTUALLY_EXCLUSIVE_INPUTS:
            if ((layer + entry1) in input_dict.keys() and (layer+entry2) in input_dict.keys()):
                raise ValueError(f"Cannot use both '{entry1}' and '{entry2}'!")
    

    # Set up dictionaries to store layer paramteres
    eps_vs_freq_mode = {}
    eps_of_f = {}
    sigma_of_f = {}
    rho = {}
    tau = {}
    epsb = {}
    tandelta = {}
    mu_vs_freq_mode = {}
    chi = {}
    fmu = {}
    mu_of_f = {}

    for p in layer_numbers:

        layer = f"Layer {p} "

        # Find which input mode was used for epsilon
        # For tan(delta) or resistivity, the dictionary entry is just a number,
        # while for frequency dependent conductivity/permittivity, the dictionary entry is set to a numpy array
        # The information of which is used is stored in eps_vs_freq_mode
        if layer + "dielectric tan(delta)" in input_dict.keys():
            eps_vs_freq_mode[p] = 3
            tandelta[p] = read_input(input_dict,fMP, layer+"dielectric tan(delta)")
            epsb[p] = read_input(input_dict, fMP, layer+"real part of dielectric constant", default="1")

        elif layer + "frequency dependent conductivity" in input_dict.keys():
            eps_vs_freq_mode[p] = 2
            sigma_of_f[p] = np.loadtxt(input_dict[layer + "frequency dependent conductivity"], skiprows=1)
            epsb[p] = read_input(input_dict, fMP, layer+"real part of dielectric constant", default="1")
            
        elif layer + "frequency dependent relative complex permittivity" in input_dict.keys():
            eps_vs_freq_mode[p] = 1
            eps_of_f[p] = np.loadtxt(input_dict[layer + "frequency dependent relative complex permittivity"], skiprows=1)

        else:
            eps_vs_freq_mode[p] = 0
            rho[p] = read_input(input_dict,fMP, layer+"DC resistivity (Ohm.m)",default="1e5")
            tau[p] = read_input(input_dict, fMP, layer+"relaxation time for resistivity (ps)", default="0")
            tau[p] *= 1e-12 # ps -> s
            epsb[p] = read_input(input_dict, fMP, layer+"real part of dielectric constant", default="1")

        # Find which input mode was used for mu
        # The information of which mode is used is stored in mu_vs_freq_mode
        if layer + "frequency dependent relative complex permeability" in input_dict.keys():
            mu_vs_freq_mode[p] = 1
            mu_of_f[p] = np.loadtxt(input_dict[layer + "frequency dependent relative complex permeability"], skiprows=1)
        else:
            mu_vs_freq_mode[p] = 0
            chi[p] = read_input(input_dict, fMP, layer+"magnetic susceptibility", default="0")
            fmu[p] = read_input(input_dict, fMP, layer+"relaxation frequency of permeability (MHz)", default="Infinity")
            fmu[p] *= 1e6 # MHz -> Hz
    

    # Make dict to more easily pass all arguments as kwargs to a function
    layer_properties = {
    "eps_vs_freq_mode" : eps_vs_freq_mode,
    "eps_of_f" : eps_of_f,
    "sigma_of_f" : sigma_of_f,
    "rho" : rho,
    "tau" : tau,
    "epsb" : epsb,
    "tandelta" : tandelta,
    "mu_vs_freq_mode" : mu_vs_freq_mode,
    "chi" : chi,
    "fmu" : fmu,
    "mu_of_f" : mu_of_f
    }


    b0              =   read_input(input_dict, fMP, "Layer 1 inner half gap in mm", default="2")
    b0              *=  1e-3 # mm -> m
    thick           =   read_input_layer(input_dict, N, fMP, "thickness in mm", default="Infinity")
    thick           *=  1e-3 # mm -> m
    b               =   np.cumsum(np.concatenate([[b0], thick]))

    if not flag_topbotsym:
        b0m         =   - read_input(input_dict, fMP, "Layer -1 inner half gap in mm", default="2")
        b0m         *=  1e-3 # mm -> m
        thickm      =   - read_input_layer(input_dict, -M, fMP, "thickness in mm", default="Infinity")
        thickm      *=  1e-3 # mm -> m
        bm          =   np.cumsum(np.concatenate([[b0m], thickm]))


    # Get frequency range and generate array
    freqs = generate_frequency_range(input_dict)

    filename_base = f"W{machine}_{N}layersup_{M}layersdown{1e3*b0.toDouble():.2f}mm{commentoutput}"

    # Make copy of input file
    shutil.copy(input_file, f"InputData{filename_base}.dat")

    # Set up integration workspace
    w = gsl_integration_workspace_alloc(INTEGRATION_WORKSPACE_LIMIT)

    # Set up lookup tables for matrix multiplication results
    memory = memorycontainer(MAXMEM)

    if order <= 0: #This should really never be negative

        # Set up output arrays
        Zlong = np.zeros_like(freqs, dtype=complex)
        Zycst = np.zeros_like(freqs, dtype=complex)
        Zxdip = np.zeros_like(freqs, dtype=complex)
        Zydip = np.zeros_like(freqs, dtype=complex)
        Zxquad = np.zeros_like(freqs, dtype=complex)
        Zyquad = np.zeros_like(freqs, dtype=complex)


        ######## Loop over frequencies, calculate impedance #########
        for i, freq in enumerate(freqs):
            memory.mem = 0 # Reset memory for each frequency

            omega = twopiMP*freq # ampf
            k = omega/(beta*C) # ampf
            kovergamma = k/gamma # ampf

            # Calculate eps1 and mu1 from the layer properties
            eps1, mu1 = calculate_eps1_and_mu1(N, freq, omega, **layer_properties)

            if not flag_topbotsym:
                eps1m, mu1m = calculate_eps1_and_mu1(-M, freq, omega, **layer_properties)


            # Make ap::template_1d_arrays from the numpy arrays
            b_arr   = make_1d_array((1,N+1), fMP, b)
            eps1_arr= make_1d_array((1,N+1), cMP, eps1)
            mu1_arr = make_1d_array((1,N+1), cMP, mu1)

            if flag_topbotsym:
                bm_arr      = make_1d_array((1,N+1), fMP, -b)
                eps1m_arr   = make_1d_array((1,N+1), cMP, eps1)
                mu1m_arr    = make_1d_array((1,N+1), cMP, mu1)
            else:
                bm_arr = make_1d_array((1, M+1),fMP,bm)
                eps1m_arr = make_1d_array((1, M+1), cMP, eps1m)
                mu1m_arr = make_1d_array((1,M+1), cMP, mu1m)
            

            # Calculate the alpha coefficients (in c++)
            alpha00 = alphamn(flag_topbotsym,M,N,b_arr,bm_arr,beta,eps1_arr,eps1m_arr,mu1_arr,mu1m_arr,omega,k,kovergamma,0,0,INTEGRATION_WORKSPACE_LIMIT,w,memory) # complex
            alpha01 = alphamn(flag_topbotsym,M,N,b_arr,bm_arr,beta,eps1_arr,eps1m_arr,mu1_arr,mu1m_arr,omega,k,kovergamma,0,1,INTEGRATION_WORKSPACE_LIMIT,w,memory) # complex
            alpha02 = alphamn(flag_topbotsym,M,N,b_arr,bm_arr,beta,eps1_arr,eps1m_arr,mu1_arr,mu1m_arr,omega,k,kovergamma,0,2,INTEGRATION_WORKSPACE_LIMIT,w,memory) # complex
            alpha11 = alphamn(flag_topbotsym,M,N,b_arr,bm_arr,beta,eps1_arr,eps1m_arr,mu1_arr,mu1m_arr,omega,k,kovergamma,1,1,INTEGRATION_WORKSPACE_LIMIT,w,memory) # complex

            # Compute all impedances from the alphas
            cst = 1j*L*(k*Z0/(beta*amp.sqr(gamma)*twopiMP)).toDouble() # complex
            Zlong[i] = cst*alpha00 # complex
            Zycst[i] = cst*alpha01/(gamma.toDouble()) # complex
            cst = cst*(k/(twoMP*amp.sqr(gamma))).toDouble()  # complex
            Zxdip[i] = cst*(alpha02-alpha00) # complex
            Zydip[i] = 2*cst*alpha11 # complex
            Zxquad[i] = -Zxdip[i] # complex
            Zyquad[i] =cst*(alpha02+alpha00) # complex

            ####### Save output to files ########
            for name, data, unit in zip(
                ["Zlong","Zycst","Zxdip", "Zydip", "Zxquad", "Zyquad"],
                [Zlong,   Zycst,  Zxdip,   Zydip,   Zxquad,   Zyquad],
                ["Ohm",  "Ohm",  "Ohm/m", "Ohm/m", "Ohm/m",  "Ohm/m"]
            ):
                np.savetxt(
                    f"{name}{filename_base}.dat",
                    np.transpose([freqs, np.real(data), np.imag(data)]),
                    fmt="%13.8e",
                    header=f"Frequency [Hz]\tRe({name}) [{unit}]\tIm({name}) [{unit}]",
                    comments=""
                )

    else:
        # if order > 0

        # initialize impedance array. Shape: (#planes, #power_combination, #frequencies)
        # Note: this is much larger than needed, since many power/plane combinations are invalid
        Z = np.zeros((3, int((order+1)**4), len(freqs)), dtype=complex)

        # initialize list of combinations of plane and powers to actually save
        saving_info = []

        # loop over all frequencies, calculate impedances
        for i, freq in enumerate(freqs):

            # reset memory
            memory.mem=0

            alphas: Dict[Tuple[int, int], complex] = {}

            # calculate frequency dependent variables
            omega=twopiMP*freq # ampf
            k=omega/(beta*C) # ampf
            kovergamma=k/gamma # ampf

             # Calculate eps1 and mu1 from the layer properties
            eps1, mu1 = calculate_eps1_and_mu1(N, freq, omega, **layer_properties)

            if not flag_topbotsym:
                eps1m, mu1m = calculate_eps1_and_mu1(-M, freq, omega, **layer_properties)

            # make ap::template_1d_arrays from numpy arrays
            b_arr = make_1d_array((1,N+1), fMP, b)
            eps1_arr = make_1d_array((1,N+1), cMP, eps1)
            mu1_arr = make_1d_array((1,N+1), cMP, mu1)

            if flag_topbotsym:
                bm_arr = make_1d_array((1,N+1), fMP, -b)
                eps1m_arr = make_1d_array((1,N+1), cMP, eps1)
                mu1m_arr = make_1d_array((1,N+1), cMP, mu1)
            else:
                bm_arr = make_1d_array((1, M+1),fMP,bm)
                eps1m_arr = make_1d_array((1, M+1), cMP, eps1m)
                mu1m_arr = make_1d_array((1,M+1), cMP, mu1m)


            # Loop over all planes
            for plane, plane_char in enumerate(["l","x","y"]):

                if (plane==0):
                    cst=4j*L*(k*Z0/(beta*amp.sqr(gamma)*twopiMP)).toDouble() # complex
                else:
                    cst=4j*L*(Z0/(beta*amp.sqr(gamma)*twopiMP)).toDouble() # complex

                # Loop over powers of the coordinates, filter to only valid combinations
                for j, (power_a, power_b, power_c, power_d) in enumerate(product(range(order+1), repeat=4)):
                    sumpower = power_a + power_b + power_c + power_d
                    if (
                        ((power_a+power_c+(plane==1))%2 == 0) and
                        ((sumpower)<=order) and
                        (((sumpower+(plane>0))%2)*flag_topbotsym==0)
                    ):

                        aplusc2=(power_a+power_c+(plane==1))//2

                        # Loop over r and q
                        for r in range(power_b//2 + 1):
                            for q in range(aplusc2+(power_d+(plane==2))//2 + 1):

                                m = power_b-2*r
                                n = power_a+power_c+power_d+(plane>0)-2*q

                                if plane == 0:
                                    fact = 1
                                elif plane == 1:
                                    fact = power_c+1
                                else:
                                    fact = power_d+1
                                
                                coef = fact*coefmn(power_a,power_b,power_c+(plane==1),power_d+(plane==2),m,n,r,q)
                        
                                # Calculate alpha coefficienct (in c++)
                                if alphas.get((m,n), None) is None:
                                    alphas[(m,n)] = alphamn(flag_topbotsym,M,N,b_arr,bm_arr,beta,eps1_arr,eps1m_arr,mu1_arr,mu1m_arr,omega,k,kovergamma,m,n,INTEGRATION_WORKSPACE_LIMIT,w,memory)


                                # Calculate impedance contribution
                                Z[plane, j, i] = Z[plane, j, i] + alphas[(m,n)]*cst*coef*(kovergamma.toDouble()/2.)**(sumpower+(plane>0))


                        # prepare info for saving. Only done in first iteration
                        if i==0:
                            # Set units and file names
                            if sumpower == 0:
                                unit = "Ohm"
                                if plane==0:
                                    name = f"Zlong"
                                else:
                                    name = f"Zycst" # necessarily, plane==2 (it cannot be hor. because of the condition above (a+c+(plane==1))%2 == 0)
                            elif sumpower == 1:
                                unit = "Ohm/m"
                                if ( power_a==1 and plane==1 ):
                                    name = f"Zxdip"
                                elif ( power_b==1 and plane==2 ):
                                    name = f"Zydip"
                                elif ( power_c==1 and plane==1 ):
                                    name = f"Zxquad"
                                elif ( power_d==1 and plane==2 ):
                                    name = f"Zyquad"
                                else:
                                    unit = f"Ohm/m^{sumpower}"
                                    name = f"Z{plane_char}{power_a}{power_b}{power_c}{power_d}"
                            else:
                                unit = f"Ohm/m^{sumpower}"
                                name = f"Z{plane_char}{power_a}{power_b}{power_c}{power_d}"
                            
                            # Store the indices of the Z array that should be saved to files, as well as the name and unit
                            saving_info.append((plane,j,unit,name))

                            # in top/bottom symmetric case, still fill (with zeros) the Zycst term (for back compatibility)
                            if ( flag_topbotsym and sumpower==0 and plane==2 ):
                                np.savetxt(
                                    f"Zycst{filename_base}.dat",
                                    np.transpose([freqs, np.zeros_like(freqs), np.zeros_like(freqs)]),
                                    fmt="%13.8e",
                                    header=f"Frequency [Hz]\tRe(Zycst) [Ohm]\tIm(Zycst) [Ohm]",
                                    comments=""
                                )

        # save to .dat files   
        for plane, j, unit, name in saving_info:
            np.savetxt(
                f"{name}{filename_base}.dat",
                np.transpose([freqs, np.real(Z[plane,j,:]), np.imag(Z[plane,j,:])]),
                fmt="%13.8e",
                header=f"Frequency [Hz]\tRe({name if name.isalpha() else 'Z'}) [{unit}]\tIm({name if name.isalpha() else 'Z'}) [{unit}]",
                comments=""
            )                
    
    # free memory used for integration
    gsl_integration_workspace_free(w)


def calculate_eps1_and_mu1(
    N,
    freq,
    omega, 
    eps_vs_freq_mode,
    mu_vs_freq_mode,
    epsb, rho, tau, eps_of_f, sigma_of_f, tandelta,
    chi, fmu, mu_of_f) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate eps1 and mu1 for all layers, using the method specified in eps_vs_freq_mode and mu_vs_freq_mode

    Give a negative N to get the negative layers
    
    Returns two numpy arrays with the results
    """

    # Set up arrays
    eps1 = np.zeros(np.abs(N)+1, dtype=object)
    mu1 = np.zeros(np.abs(N)+1, dtype=object)

    # Inner layer is always vacuum
    eps1[0] = cMP(1)
    mu1[0] = cMP(1)

    for p in range(1,np.abs(N)+1):
        p_sign = p*np.sign(N)
        # Calculate eps1
        if eps_vs_freq_mode[p_sign] == 0:
            # DC Resistivity and relaxation time
            if rho[p_sign].isFiniteNumber():
                eps1[p]=epsb[p_sign]+oneMP/(jimagMP*eps0*rho[p_sign]*omega*(oneMP+jimagMP*omega*tau[p_sign]))
            else:
                eps1[p]=epsb[p_sign]

        elif eps_vs_freq_mode[p_sign] == 1:
            # Relative permittivity as function of frequency
            eps_frequencies, eps_real, eps_imag = eps_of_f[p_sign].T
            eps_frequencies = array("d",eps_frequencies)
            eps_values = eps_real + 1j*eps_imag
            eps1[p]=get_value_from_interp(freq, eps_frequencies, eps_values, len(eps_frequencies))
            
        elif eps_vs_freq_mode[p_sign] == 2:
            # Conductivity as function of frequency
            sigma_frequencies, sigma_real, sigma_imag = sigma_of_f[p_sign].T
            sigma_frequencies = array("d", sigma_frequencies)
            sigma_values = sigma_real + 1j*sigma_imag
            sigma = get_value_from_interp(freq, sigma_frequencies, sigma_values, len(sigma_frequencies))
            eps1[p] = epsb[p_sign] + sigma/(jimagMP*eps0*omega)
        
        else: # eps_vs_freq_mode[p] == 3:
            # tan(delta)
            eps1[p]=epsb[p_sign]*(oneMP-jimagMP*tandelta[p_sign])
        
        # Calculate mu1
        if mu_vs_freq_mode[p_sign] == 0:
            # Susceptibility and relaxation frequency
            mu1[p]=oneMP+chi[p_sign]/(oneMP+jimagMP*omega/(twopiMP*fmu[p_sign]))

        else:
            # Relative permeability as function of frequency
            mu_frequencies, mu_real, mu_imag = mu_of_f[p_sign].T
            mu_frequencies = array("d", mu_frequencies)
            mu_values = mu_real + 1j*mu_imag
            mu1[p]=get_value_from_interp(freq, mu_frequencies, mu_values, len(mu_frequencies))
        
    return eps1, mu1
