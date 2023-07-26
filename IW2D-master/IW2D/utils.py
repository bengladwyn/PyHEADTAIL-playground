from pathlib import Path
import numpy as np
from typing import Union, Any, Callable, Sequence, Optional
from .interface import IW2DLayer

from .constants import cMP, fMP
from cppyy.gbl import ap

def linear_interpolation_from_file(data_file: Union[Path, str], **loadtxt_kwargs) -> Callable[[float], complex]:
    """### Description
    Interpolate a function z(x): real -> complex from a file
    
    Read a file with three columns: input variable, Real(ouput variable) and Imag(output variable),
    and make a function that linearly interpolates these values.

    Uses numpy.loadtxt to load the file. Any options can be passed to this through **kwargs.

    ### Parameters
    data_file : Union[Path, str]
        The path to the file to make the interpolation from

    ### Returns
    Callable[[float], complex]
        A function that takes an input variable and returns the interpolated output that co
    """
    frequencies, values_real, values_imag = np.loadtxt(data_file, **loadtxt_kwargs).T
    def interpolation(freq: float) -> complex:
        return np.interp(freq, frequencies, values_real + 1j*values_imag)
    return interpolation

def make_cMP(z: Any):
    """### Description
    Create a multiprecision complex number (amp::campf[PRECISION])

    If z is NOT a python complex, the imaginary component of the output will be 0.

    ### Parameters
    z : Any
        Either a python complex or anything parsable by the amp::campf constructor

    ### Returns
    cMP
        z as a complex multiprecision float
    """
    if isinstance(z, complex):
        zMP = cMP()
        zMP.x = fMP(z.real)
        zMP.y = fMP(z.imag)
        return zMP
    else:
        return cMP(z)

def make_layer_arrays(layers: Sequence[IW2DLayer], b0: Any, frequency: float, b_negative: bool = False):
    """### Description
    Makes the ap::template_1d_arrays for eps1, mu1 and b (thickness) used by the c++ calculation scripts,
    for a list of layers and a given frequency.

    b is retrieved from the layers' thickness attributes, while eps1 and mu1 are generated with the layers' corresponding
    callable attributes using frequency as input.

    ### Parameters
    layers : Sequence[IW2DLayer]
        A sequence of IW2DLayers
    b0 : float, str (anything parsable by fMP constructor)
        The half-gap in m of the innermost vacuum region.
    frequency : float
        A frequency in Hz
    b_negative : bool, optional
        Set to true if the entries in b array should have a negative sign


    ### Returns
    ap.template_1d_array, ap.template_1d_array, ap.template_1d_array
        The arrays of eps1 (cMP), mu1 (cMP) and b (fMP) for all layers given in the input.
    """

    # get sign of b
    b_sign = -1 if b_negative else 1

    # Initialize arrays
    eps1 = ap.template_1d_array[cMP]()
    mu1  = ap.template_1d_array[cMP]()
    b    = ap.template_1d_array[fMP]()

    # Set bounds of arrays
    eps1.setbounds(1, len(layers)+1)
    mu1.setbounds(1, len(layers)+1)
    b.setbounds(1, len(layers)+1)

    # Set innermost vacuum layer
    eps1[1] = cMP(1)
    mu1[1]  = cMP(1)
    b[1]    = fMP(b_sign*b0)

    # loop over layers, fill arrays
    for i, layer in enumerate(layers):
        p = i+2 # +1 to compensate for innermost layer being vacuum, +1 to compensate for 1-indexation

        eps1[p] = make_cMP(layer.eps1(frequency))
        mu1[p]  = make_cMP(layer.mu1(frequency))
        b[p]    = b[p-1] + fMP(layer.thickness)*b_sign
    
    return eps1, mu1, b


def py_complex(z):
    """Convert an amp::campf complex number to the corresponding python complex."""
    return z.x.toDouble() + 1j*z.y.toDouble()
