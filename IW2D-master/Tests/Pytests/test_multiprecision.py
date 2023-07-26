from typing import List, Tuple, Callable
from _pytest.mark.structures import ParameterSet
import pytest
import numpy as np
from scipy.special import iv, kv, ive, kve
import cppyy
import IW2D


# Bring into namespace.
# Note: Only available here because multiprecision.cc is included in conftest
from cppyy.gbl import amp, cexpMP, csqrtMP, clogMP


######### Define constants ##########

# Number of mantissa bits in MP numbers
PRECISION = cppyy.gbl.precision

# Relative and absolute error allowed in float equality checks
# Note: double equality assertions check double_1 == double_2 +-max(rel_tol, abs_tol), not with min() as one might expect 
RELATIVE_TOLERANCE_DOUBLES = 2**-52
ABSOLUTE_TOLERANCE_DOUBLES = 0

RELATIVE_TOLERANCE_MP = 2**-(PRECISION-2)

###### Define utility functions #####

def complex_to_campf(z: complex) -> amp.campf[PRECISION]:
    """Make a amp::campf from a python complex number"""
    z_MP = amp.campf[PRECISION]()
    z_MP.x = amp.ampf[PRECISION](z.real)
    z_MP.y = amp.ampf[PRECISION](z.imag)
    return z_MP


def relative_error(z1: amp.campf[PRECISION], z2: amp.campf[PRECISION]) -> Tuple[float, float]:
    """Gives as a tuple the relative error of the real and imag components of two campf numbers z1 and z2.
    The relative error is given by |x1 - x2|/max(|x1|,|x2|), where x1, x2 are the real or imaginary components of z1, z2."""
    
    if z1.x.toDouble() == 0 and z2.x.toDouble() == 0:
        relative_error_real = 0.
    else:
        relative_error_real = (amp.abs(z1.x - z2.x)/amp.maximum(amp.abs(z1.x), amp.abs(z2.x))).toDouble()
    
    if z1.y.toDouble() == 0 and z2.y.toDouble() == 0:
        relative_error_imag = 0.
    else:
        relative_error_imag = (amp.abs(z1.y - z2.y)/amp.maximum(amp.abs(z1.y), amp.abs(z2.y))).toDouble()
    
    return (relative_error_real, relative_error_imag)
    

def probe_points_complex_plane(magnitudes: List[float] = [1.], only_first_quadrant: bool = False) -> List[ParameterSet]:
    """Generates 8 complex numbers for each specified magnitude pointing in the directions n*pi/4, n integer, 0 <= n < 8.
    If only_first_quadrant is set to True, n will be limited to {0,1,2}"""

    directions = 8
    if only_first_quadrant:
        directions = 3
    
    points = []
    for magnitude in magnitudes:
        for n in range(directions):
            x = magnitude*np.cos(n*np.pi/4)
            y = magnitude*np.sin(n*np.pi/4)
            z = x + 1j*y
            points.append(pytest.param(z, id=f"{magnitude:.2e} * exp({n}*pi/4*j)"))
    return points

########### Define tests ############

@pytest.mark.parametrize("input_z",probe_points_complex_plane(magnitudes = [1.e-100, 1., 20., 1.e100]))
def test_clogMP_against_numpy(input_z: complex) -> None:
    """Test that the clogMP function in multiprecision.cc gives the same result as numpy."""
    numpy_answer = np.log(input_z)
    mp_answer = clogMP(complex_to_campf(input_z))
    assert mp_answer.x.toDouble() == pytest.approx(np.real(numpy_answer), rel=RELATIVE_TOLERANCE_DOUBLES, abs=ABSOLUTE_TOLERANCE_DOUBLES)
    assert mp_answer.y.toDouble() == pytest.approx(np.imag(numpy_answer), rel=RELATIVE_TOLERANCE_DOUBLES, abs=ABSOLUTE_TOLERANCE_DOUBLES)

@pytest.mark.parametrize("input_z",probe_points_complex_plane(magnitudes=[0.01, 1., 100]))
def test_cexpMP_against_numpy(input_z: complex) -> None:
    """Test that the cexpMP function in multiprecision.cc gives the same result as numpy."""
    numpy_answer = np.exp(input_z)
    mp_answer = cexpMP(complex_to_campf(input_z))
    assert mp_answer.x.toDouble() == pytest.approx(np.real(numpy_answer), rel=RELATIVE_TOLERANCE_DOUBLES, abs=ABSOLUTE_TOLERANCE_DOUBLES)
    assert mp_answer.y.toDouble() == pytest.approx(np.imag(numpy_answer), rel=RELATIVE_TOLERANCE_DOUBLES, abs=ABSOLUTE_TOLERANCE_DOUBLES)

@pytest.mark.parametrize("input_z", probe_points_complex_plane(magnitudes=[1e-100, 1., 1e100]))
def test_csqrtMP_against_numpy(input_z: complex) -> None:
    """Test that the csqrtMP function in multiprecision.cc gives the same result as numpy."""
    numpy_answer = np.sqrt(input_z)
    mp_answer = csqrtMP(complex_to_campf(input_z))
    assert mp_answer.x.toDouble() == pytest.approx(np.real(numpy_answer), rel=RELATIVE_TOLERANCE_DOUBLES, abs=ABSOLUTE_TOLERANCE_DOUBLES)
    assert mp_answer.y.toDouble() == pytest.approx(np.imag(numpy_answer), rel=RELATIVE_TOLERANCE_DOUBLES, abs=ABSOLUTE_TOLERANCE_DOUBLES)

@pytest.mark.parametrize("input_z",[
    ("5", "0"),
    ("-1", "3"),
    ("14", "-1.58"),
    ("65", "0.8")
])
def test_log_of_exp(input_z: Tuple[str, str]) -> None:
    """Test that taking log(exp(z)), where z is a complex number, is the identity operation on z.
    input_z is a tuple on the format ("z_real", "z_imag")"""

    z = amp.campf[PRECISION]()
    z.x, z.y = input_z

    output_1 = z
    output_2 = clogMP(cexpMP(z)) # Note: requires -pi < imag(z) <= pi for equality

    relative_error_real, relative_error_imag = relative_error(output_1, output_2)
    assert relative_error_real < RELATIVE_TOLERANCE_MP
    assert relative_error_imag < RELATIVE_TOLERANCE_MP

@pytest.mark.parametrize("input_z",[
    ("5", "0"),
    ("-1", "3"),
    ("14", "-1.58"),
    ("65", "0.8")
])
def test_log_of_sqrt(input_z: Tuple[str, str]) -> None:
    """Test that taking log(sqrt(z)) gives 0.5*log(z) for a complex number z.
    input_z is a tuple on the format ("z_real", "z_imag")"""

    z = amp.campf[PRECISION]()
    z.x, z.y = input_z

    output_1 = clogMP(csqrtMP(z))
    output_2 = clogMP(z) / amp.ampf[PRECISION](2)

    relative_error_real, relative_error_imag = relative_error(output_1, output_2)
    assert relative_error_real < RELATIVE_TOLERANCE_MP
    assert relative_error_imag < RELATIVE_TOLERANCE_MP

@pytest.mark.parametrize("func,expected_output",[(cexpMP,1.),(csqrtMP,0.),(clogMP,"-@Inf@")])
def test_zero_input(func: Callable, expected_output) -> None:
    """Test that inputting 0+0j to the functions cexpMP, csqrtMP and clogMP gives the expected results."""

    # Declare complex zero = 0 + 0j
    complex_zero_MP = amp.campf[PRECISION]()

    output = func(complex_zero_MP)
    expected_output_MP = amp.campf[PRECISION](expected_output)

    assert output == expected_output_MP

@pytest.mark.parametrize("input_1,input_2",[
    pytest.param(("-0", "0"), ( "0", "0"), id="(+-0+0j)"),
    pytest.param(("-0","-0"), ( "0", "0"), id="(+-0-0j)"),
    pytest.param(( "0","-0"), ( "0", "0"), id="(+0+-0j)"),
    pytest.param(( "1", "0"), ( "1","-0"), id="(+1+-0j)"),
    pytest.param(("-1", "0"), ("-1","-0"), id="(-1+-0j)"),
    pytest.param(( "0", "1"), ("-0", "1"), id="(+-0+1j)"),
    pytest.param(( "0","-1"), ("-0","-1"), id="(+-0-1j)")
    ])
def test_cexpMP_negative_zeros(input_1: Tuple[str, str], input_2: Tuple[str, str]) -> None:
    """Test if cexpMP treats +0 and -0 the same in various contexts.
    input_1 and input_2 are tuples, each containing 2 strings used to initialize the
    real and imaginary part of a complex number, respectively."""

    input_1_mp = amp.campf[PRECISION]()
    input_1_mp.x, input_1_mp.y = input_1
    output_1_mp = cexpMP(input_1_mp)

    input_2_mp = amp.campf[PRECISION]()
    input_2_mp.x, input_2_mp.y = input_2
    output_2_mp = cexpMP(input_2_mp)

    assert output_1_mp == output_2_mp

@pytest.mark.parametrize("input_1,input_2",[
    pytest.param(("-0", "0"), ( "0", "0"), id="(+-0+0j)"),
    pytest.param(("-0","-0"), ( "0", "0"), id="(+-0-0j)"),
    pytest.param(( "0","-0"), ( "0", "0"), id="(+0+-0j)"),
    pytest.param(( "1", "0"), ( "1","-0"), id="(+1+-0j)"),
    pytest.param(("-1", "0"), ("-1","-0"), id="(-1+-0j)", marks=pytest.mark.xfail(strict=True, reason="Branch cut")),
    pytest.param(( "0", "1"), ("-0", "1"), id="(+-0+1j)"),
    pytest.param(( "0","-1"), ("-0","-1"), id="(+-0-1j)")
    ])
def test_csqrtMP_negative_zeros(input_1: Tuple[str, str], input_2: Tuple[str, str]) -> None:
    """Test if csqrtMP treats +0 and -0 the same in various contexts.
    input_1 and input_2 are tuples, each containing 2 strings used to initialize the
    real and imaginary part of a complex number, respectively."""

    input_1_mp = amp.campf[PRECISION]()
    input_1_mp.x, input_1_mp.y = input_1
    output_1_mp = csqrtMP(input_1_mp)

    input_2_mp = amp.campf[PRECISION]()
    input_2_mp.x, input_2_mp.y = input_2
    output_2_mp = csqrtMP(input_2_mp)

    assert output_1_mp == output_2_mp

# @pytest.mark.xfail # Marked as expected fail until the issue has been addressed.
@pytest.mark.parametrize("input_1,input_2",[
    pytest.param(( "1", "0"), ( "1","-0"), id="(+1+-0j)"),
    pytest.param(("-1", "0"), ("-1","-0"), id="(-1+-0j)", marks=pytest.mark.xfail(strict=True, reason="Branch cut")), # Test passes only if assertion fails
    pytest.param(( "0", "1"), ("-0", "1"), id="(+-0+1j)"),
    pytest.param(( "0","-1"), ("-0","-1"), id="(+-0-1j)")
    ]) # Probably don't need all of these
def test_clogMP_negative_zeros(input_1: Tuple[str, str], input_2: Tuple[str, str]) -> None:
    """Test if clogMP treats +0 and -0 the same in various contexts.
    input_1 and input_2 are tuples, each containing 2 strings used to initialize the
    real and imaginary part of a complex number, respectively."""

    input_1_mp = amp.campf[PRECISION]()
    input_1_mp.x, input_1_mp.y = input_1
    output_1_mp = clogMP(input_1_mp)

    input_2_mp = amp.campf[PRECISION]()
    input_2_mp.x, input_2_mp.y = input_2
    output_2_mp = clogMP(input_2_mp)

    assert output_1_mp == output_2_mp

@pytest.mark.parametrize("input_z", probe_points_complex_plane(magnitudes=[1.,10.,50.,70.], only_first_quadrant=False))
@pytest.mark.parametrize("besseltype", ["besi","besk"])
@pytest.mark.parametrize("input_order", [0,1,2,3,4])
def test_bessel_against_scipy(besseltype: str, input_z: complex, input_order: int) -> None:
    """Test that the bessel function implementation in multiprecision.cc gives the same result as the scipy implementation.
    The equality check only fails if the compared numbers differ by more than BOTH relative_tolerance and absolute_tolerance.
    
    arguments:
        besseltype (str)            :   which bessel function output to check, one of: "besi","besk"
        input_z (complex)           :   The argument z to input into the bessel function
        input_order (int)           :   The order of the bessel function"""
    
    # Set tolerances
    # Both checks must fail for test to fail
    relative_tolerance=1e-14
    absolute_tolerance=1e-14

    # Instantiate variables
    besi = amp.campf[PRECISION]()
    besk = amp.campf[PRECISION]()

    # Run bessel function. Outputs are inserted into besi, besk, besinorm and besknorm (which are passed by reference)
    cppyy.gbl.bessel(input_order, complex_to_campf(input_z), besi, besk)

    # Select the bessel function chosen for this test
    if besseltype == "besi":
        ans = besi
        expected = iv(input_order, input_z)

    elif besseltype == "besk":
        ans = besk
        expected = kv(input_order, input_z)

    else:
        raise ValueError(f"Got bessel type argument: {besseltype}, expected from ['besi','besk']")
    
    assert ans.x.toDouble() == pytest.approx(np.real(expected), rel=relative_tolerance, abs=absolute_tolerance)
    assert ans.y.toDouble() == pytest.approx(np.imag(expected), rel=relative_tolerance, abs=absolute_tolerance)
