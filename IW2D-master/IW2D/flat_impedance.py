from .interface import FlatIW2DInput
from .constants import fMP, twopiMP, Z0, twoMP, C, oneMP
from .utils import make_layer_arrays
from typing import Any,  Dict, Tuple
import numpy as np
import pandas as pd
from numpy.typing import ArrayLike

from cppyy.gbl import amp, gsl_integration_workspace_alloc, memorycontainer, alphamn, gsl_integration_workspace_free


def _iw2d_flat_impedance_single_frequency(
    input_obj: FlatIW2DInput,
    frequency: float,
    integration_workspace_limit: int = 1000,
    maxmem: int = 50000
) -> np.ndarray:
    """### Description
    Calculate the impedance at the given frequency in a flat chamber setup defined in input_obj.

    ### Parameters
    input_obj : FlatIW2DInput
        An input object containing all relevant information about the setup
    frequency : float
        The frequency in Hz to evaluate the impedance at

    ### Returns
    np.ndarray (shape=(6,), dtype=complex)
        An array of the impedances Zlong, Zycst, Zxdip, Zydip, Zxquad, Zyquad
    """
    
    N = len(input_obj.top_layers)
    M = len(input_obj.bottom_layers) if input_obj.bottom_layers is not None else N
    gammaMP = fMP(input_obj.relativistic_gamma) # ampf
    beta = amp.sqrt(oneMP-oneMP/amp.sqr(gammaMP)) # ampf
    omega = twopiMP*fMP(frequency) # ampf
    k = omega/(beta*C) # ampf
    kovergamma = k/gammaMP # ampf

    w = gsl_integration_workspace_alloc(integration_workspace_limit)

    memory = memorycontainer(maxmem)


    # Calculate layer parameters for the top layers
    eps1, mu1, b = make_layer_arrays(input_obj.top_layers, input_obj.top_half_gap, frequency)

    # Calculate layer parameters for the bottom layers. Note that b0 must be negative
    if input_obj.top_bottom_symmetry:
        eps1m, mu1m, bm = make_layer_arrays(input_obj.top_layers, input_obj.top_half_gap, frequency, b_negative=True)
    else:
        eps1m, mu1m, bm = make_layer_arrays(input_obj.bottom_layers, input_obj.bottom_half_gap, frequency, b_negative=True) # type: ignore

    # Compute alphamn needed for first order impedances (in C++)
    alpha00 = alphamn(input_obj.top_bottom_symmetry,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,0,0,integration_workspace_limit,w,memory) # complex
    alpha01 = alphamn(input_obj.top_bottom_symmetry,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,0,1,integration_workspace_limit,w,memory) # complex
    alpha02 = alphamn(input_obj.top_bottom_symmetry,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,0,2,integration_workspace_limit,w,memory) # complex
    alpha11 = alphamn(input_obj.top_bottom_symmetry,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,1,1,integration_workspace_limit,w,memory) # complex

    # Compute all impedances from the alphas
    cst = 1j*input_obj.length*(k*Z0/(beta*amp.sqr(gammaMP)*twopiMP)).toDouble() # complex
    Zlong = cst*alpha00 # complex
    Zycst = cst*alpha01/(gammaMP.toDouble()) # complex
    cst = cst*(k/(twoMP*amp.sqr(gammaMP))).toDouble()  # complex
    Zxdip = cst*(alpha02-alpha00) # complex
    Zydip = 2*cst*alpha11 # complex
    Zxquad = -Zxdip # complex
    Zyquad =cst*(alpha02+alpha00) # complex

    # Free memory
    gsl_integration_workspace_free(w)

    return np.array([Zlong, Zycst, Zxdip, Zydip, Zxquad, Zyquad], dtype=complex)


def iw2d_flat_impedance(
    input_obj: FlatIW2DInput,
    frequencies: ArrayLike,
    integration_workspace_limit: int = 1000,
    maxmem: int = 50000
) -> Tuple[pd.DataFrame, Dict[str, Dict[str, Any]]]:
    """### Description
    Calculate impedances for the frequencies given as input for a flat chamber and return in a pandas.DataFrame.
    The impedances calculated are longitudinal, constant y, dipolar x and y, and quadrupolar x and y.

    Also returns a metadata dictionary with units, planes and exponents.

    ### Parameters
    input_obj : FlatIW2DInput
        A FlatIW2DInput object
    frequencies : ArrayLike
        A single frequency or an iterable of frequencies in Hz

    ### Returns
    Tuple[pandas.DataFrame, Dict[str, Dict[str, Any]]]
        First element is a DataFrame with frequencies (float) as index and impedances (complex) as values
        Second element is a dictionary with the DataFrame column names as keys,
        and as values dictionaries containing the metadata for that column (units, plane, exponents)
    """

    # make index series of frequencies. Works for both iterable and non-iterable frequencies
    frequency_series = pd.Series(data=frequencies, dtype=float, name="Frequency [Hz]")
    
    # initialize array of impedances
    impedances = np.zeros((len(frequency_series), 6), dtype=complex)
    
    # loop over frequencies and calculate
    for i, frequency in enumerate(frequency_series):
        impedances[i,:] = _iw2d_flat_impedance_single_frequency(input_obj, frequency, integration_workspace_limit, maxmem)
    
    # create dataframe
    impedance_dataframe = pd.DataFrame(
        data=impedances,
        index=frequency_series,
        columns=["Zlong", "Zycst", "Zxdip", "Zydip", "Zxquad", "Zyquad"]
    )

    # create metadata dictionary
    impedance_metadata = {
        "Zlong" : {"Units":"Ohm",   "Plane":"z", "Exponents":(0,0,0,0)},
        "Zycst" : {"Units":"Ohm",   "Plane":"y", "Exponents":(0,0,0,0)},
        "Zxdip" : {"Units":"Ohm/m", "Plane":"x", "Exponents":(1,0,0,0)},
        "Zydip" : {"Units":"Ohm/m", "Plane":"y", "Exponents":(0,1,0,0)},
        "Zxquad": {"Units":"Ohm/m", "Plane":"x", "Exponents":(0,0,1,0)},
        "Zyquad": {"Units":"Ohm/m", "Plane":"y", "Exponents":(0,0,0,1)},
    }

    return impedance_dataframe, impedance_metadata
