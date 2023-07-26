from typing import Any, Dict, Tuple, Callable, Union
import numpy as np
from warnings import warn
from math import ceil

from cppyy.gbl import ap


def read_input(input_dict: Dict[str, str], output_type: Callable[[str],Any], key: str, default: Union[str, None]=None) -> Any:
    """
    Finds the input from the line in the input file starting with 'key', and converts it to the type defined by 'out_type'.

    If the input file does not contain the specified key, it raises a warning and returns a default value, if one is specified.
    If no default is specified, an error is raised.
    """
    value = input_dict.get(key, None)

    if value is None and default is None:
        err_msg = f"No argument supplied for entry: '{key}' and no default exists."
        raise ValueError(err_msg)
    
    elif value is None:
        warn_msg = f"No argument supplied for entry: '{key}'\nUsing default: '{default}'"
        warn(warn_msg)
        return output_type(default)  # type: ignore
    
    else:
        return output_type(value)


def read_input_layer(input_dict: Dict[str, str], N_layers: int, output_type: Callable[[str], Any], key: str, default: Union[Any, None]=None) -> np.ndarray:
    """
    Reads all input file entries starting with "Layer x <key>", and returns them as a numpy array.
    If N_layers is positive, it will return the layers with numbers in [1,N_layers] (both ends inclusive).
    If N_layers is negative, the layers are instead [-N_layers, -1].

    If an expected line is not found it raises a warning and returns a default value, if one is specified.
    If no default is specified, an error is raised.
    """
    entries = []
    for p in range(1, abs(N_layers)+1):
        full_key = f"Layer {np.sign(N_layers)*p} " + key
        entries.append(read_input(input_dict, output_type, full_key, default=default))
    return np.array(entries)

def generate_frequency_range(input_dict: Dict[str, str]) -> np.ndarray:
    """
    Generate a frequency range using one of three modes (as decided by 'typescan') parameter:

    0 - logarithmic frequency range from 10^fminlog to 10^fmaxlog, with nflog freq's per decade\n
    1 - linear frequency range from 10^fminlog to 10^fmaxlog, with distance fsamplin between points\n
    2 - logarithmic frequency range from 10^fminlog to 10^fmaxlog, with nflog freq's per decade, and a linear refinement with nflin frequencies between fminlin and fmaxlin

    returns sorted numpy.array with the frequencies
    """
    typescan        =   read_input(input_dict, int, "linear (1) or logarithmic (0) or both (2) frequency scan", default="0")
    if typescan not in [0,1,2]:
        msg = f"Invalid input at entry: 'linear (1) or logarithmic (0) or both (2) frequency scan'\n\tExpected 0, 1, or 2, got {typescan}"
        raise ValueError(msg)

    fminlog         =   read_input(input_dict, float, "start frequency exponent (10^) in Hz", default="2")
    fmaxlog         =   read_input(input_dict, float, "stop frequency exponent (10^) in Hz", default="13")

    if typescan in [0,2]: # log or both
        nflog       =   read_input(input_dict, int, "Number of points per decade (for log)", default="10")
        nflogtot    =   ceil(nflog*(fmaxlog-fminlog)) + 1
    
    if typescan == 1: # lin
        fsamplin    =   read_input(input_dict, float, "sampling frequency exponent (10^) in Hz (for linear)", default="8")

    if typescan == 2: # both
        fminlin =   read_input(input_dict, float, "when both, fmin of the refinement (in THz)", default="1")
        fminlin *=  1.e12
        fmaxlin =   read_input(input_dict, float, "when both, fmax of the refinement (in THz)", default="2")
        fmaxlin *=  1.e12
        nflin   =   read_input(input_dict, int, "when both, number of points in the refinement", default="100") + 1

    added_freqs =   [float(freq) for freq in input_dict.get("added frequencies [Hz]","").split()]

    if typescan == 0:
        freqs = np.logspace(fminlog, fmaxlog, nflogtot)
    elif typescan == 1:
        # freqs = np.arange(10**fminlog, 10**fmaxlog + 10**(fsamplin-1), 10**fsamplin)
        # Manual "arange" because rounding is not precisely the same
        freq_list = []
        freq = 0
        n = 0
        while freq < 10**fmaxlog:
            freq = 10**fminlog + n*10**fsamplin
            freq_list.append(freq)
            n += 1
        freqs = np.array(freq_list)

    elif typescan == 2:
        freqs = np.concatenate((
            np.logspace(fminlog, fmaxlog, nflogtot),
            np.linspace(fminlin, fmaxlin, nflin)
        ))
    
    # Note: unique also sorts
    freqs = np.unique(np.concatenate((freqs, np.array(added_freqs))))
    return freqs


def make_1d_array(bounds: Tuple[int, int], output_type: Callable[[Any], Any], input_array: np.ndarray):
    """
    Converts a numpy array to an ap::template_1d_array with the specified type and bounds.
    """

    assert ( len(input_array) == bounds[1] - bounds[0] + 1)

    output_arr = ap.template_1d_array[output_type]()
    lower, upper = bounds
    output_arr.setbounds(lower, upper)
    for i in range(lower, upper+1):
        output_arr[i] = input_array[i-lower]
    return output_arr
