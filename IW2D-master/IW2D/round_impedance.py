from .interface import RoundIW2DInput
from .constants import fMP, C, twopiMP, oneMP, Z0, jimagMP, twoMP
import numpy as np
import pandas as pd
from .utils import make_layer_arrays, py_complex
from typing import Tuple, Any, Dict
from numpy.typing import ArrayLike


from cppyy.gbl import alphaTM, amp

def _iw2d_round_impedance_single_frequency(input_obj: RoundIW2DInput, frequency: float) -> np.ndarray:
    """Calculate the impedance at the given frequency in a round chamber setup defined in input_obj.

    :param input_obj: An input object containing all relevant information about the system
    :type input_obj: RoundIW2DInput
    :param frequency: The frequency in Hz to evaluate the impedance at
    :type frequency: float
    :return: An array of the impedances Zlong, Zxdip, Zydip, Zxquad, Zyquad
    :rtype: np.ndarray
    """

    N = len(input_obj.layers)
    gammaMP = fMP(input_obj.relativistic_gamma) # ampf
    beta = amp.sqrt(oneMP-oneMP/amp.sqr(gammaMP)) # ampf
    omega = twopiMP*frequency # ampf
    k = omega/(beta*C) # ampf

    eps1, mu1, b = make_layer_arrays(input_obj.layers, input_obj.inner_layer_radius, frequency)


    # Calculate alpha coefficients (in c++)
    alphaTM0=alphaTM(N+1,0,eps1,mu1,b,k,beta) # complex
    alphaTM1=alphaTM(N+1,1,eps1,mu1,b,k,beta) # complex   

    # Compute all impedances from the alphas
    cst = jimagMP*input_obj.length*(k*Z0/(beta*amp.sqr(gammaMP)*twopiMP)) # campf
    Zlong = input_obj.yokoya_Zlong*alphaTM0*py_complex(cst) # complex
    cst = cst*(k/(amp.sqr(gammaMP)*twoMP)) # campf
    Zxdip = input_obj.yokoya_Zxdip*py_complex(cst)*alphaTM1 # complex
    Zydip = input_obj.yokoya_Zydip*py_complex(cst)*alphaTM1 # complex
    if input_obj.yokoya_Zlong == input_obj.yokoya_Zxdip == input_obj.yokoya_Zydip == 1:
        # case of axisymmetric geometry -> very small, residual quadrupolar impedance
        Zxquad = py_complex(cst)*alphaTM0 # complex
        Zyquad = Zxquad # complex
    else:
        # neglecting the small quad. impedance above whenever the Yokoya factors correspond
        # to a non-axisymmetric stucture
        Zxquad = py_complex(cst)*alphaTM1*input_obj.yokoya_Zxquad # complex
        Zyquad = py_complex(cst)*alphaTM1*input_obj.yokoya_Zyquad # complex
    
    return np.array([Zlong, Zxdip, Zydip, Zxquad, Zyquad])

def iw2d_round_impedance(input_obj: RoundIW2DInput, frequencies: ArrayLike) -> Tuple[pd.DataFrame, Dict[str, Dict[str, Any]]]:
    """Calculate impedances for the frequencies given as input for a round chamber and return in a pandas.DataFrame.
    The impedances calculated are (in return order) longitudinal, dipolar x and y, and quadrupolar x and y.

    Also returns a metadata dictionary with units, planes and exponents.

    :param input_obj: An input object containing all relevant information about the system
    :type input_obj: RoundIW2DInput
    :param frequencies: A single frequency or an iterable of frequencies, in Hz, for which to calculate the impedance
    :type frequencies: ArrayLike
    :return: First element is a DataFrame with frequencies (float) as index and impedances (complex) as values
    Second element is a dictionary with the DataFrame column names as keys,
    and as values dictionaries containing the metadata for that column (units, plane, exponents)
    :rtype: Tuple[pd.DataFrame, Dict[str, Dict[str, Any]]]
    """

    # make index series of frequencies. Works for both iterable and non-iterable frequencies
    frequency_series = pd.Series(data=frequencies, dtype=float, name="Frequency [Hz]")
    
    # initialize array of impedances
    impedances = np.zeros((len(frequency_series), 5), dtype=complex)
    
    # loop over frequencies and calculate
    for i, frequency in enumerate(frequency_series):
        impedances[i,:] = _iw2d_round_impedance_single_frequency(input_obj, frequency)
    
    # create dataframe
    impedance_dataframe = pd.DataFrame(
        data=impedances,
        index=frequency_series,
        columns=["Zlong", "Zxdip", "Zydip", "Zxquad", "Zyquad"]
    )

    # create metadata dictionary
    impedance_metadata = {
        "Zlong" : {"Units":"Ohm",   "Plane":"z", "Exponents":(0,0,0,0)},
        "Zxdip" : {"Units":"Ohm/m", "Plane":"x", "Exponents":(1,0,0,0)},
        "Zydip" : {"Units":"Ohm/m", "Plane":"y", "Exponents":(0,1,0,0)},
        "Zxquad": {"Units":"Ohm/m", "Plane":"x", "Exponents":(0,0,1,0)},
        "Zyquad": {"Units":"Ohm/m", "Plane":"y", "Exponents":(0,0,0,1)},
    }

    return impedance_dataframe, impedance_metadata
