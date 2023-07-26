from IW2D.legacy import utils
import pytest
import numpy as np


# Bring into namespace
# Note: precision is only available since multiprecision.cc is included in conftest
from cppyy.gbl import amp, precision

def test_make_1d_array():
    lower, upper = 3,7 # bounds, inclusive
    size = upper - lower + 1 # '+ 1' because inclusive
    content_type = amp.ampf[precision]

    # Generate numpy array
    numpy_array = np.arange(size)

    # Generate ap::template_1d_array with bounds (3,7) and same contents, with type 'output_type'
    ap_1d_array = utils.make_1d_array(
        bounds = (lower,upper),
        output_type = content_type,
        input_array = numpy_array
    )

    # Loop over elements
    for i in range(lower,upper+1):

        # Check that all elements are the same
        assert ap_1d_array[i].toDouble() == numpy_array[i-lower]

        # Check that all elements have the correct type
        assert type(ap_1d_array[i]) == content_type


def test_read_input():
    """Test that read_input function works as intended."""
    fMP = amp.ampf[precision]

    # Make a dict imitating the format of input dicts generated in IW2D.legacy.__main__
    input_dict = {"InputA":"123", "InputB":"456", "InputC":"789"}

    # Read A from input_dict, make str
    output_A = utils.read_input(input_dict, str, "InputA")
    expected_output_A = "123"

    # Read B from input_dict, make int
    output_B = utils.read_input(input_dict, int, "InputB")
    expected_output_B = 456

    # Read C from input dict, make multiprecision float
    output_C = utils.read_input(input_dict, fMP, "InputC")
    expected_output_C = fMP(789)

    # Check that all outputs have the correct type and value
    assert isinstance(output_A, str)
    assert output_A == expected_output_A

    assert isinstance(output_B, int)
    assert output_B == expected_output_B

    assert isinstance(output_C, fMP)
    assert output_C == expected_output_C



@pytest.mark.filterwarnings("ignore:No argument supplied")
def test_read_input_default():
    """Test supplying a default value in read_input."""
    # Make empty input dict
    input_dict = {}

    # Search for an entry which is not present in the input dict, should set it to supplied default
    output = utils.read_input(input_dict, float, "NOT PRESENT", default="3.14")
    expected_output = 3.14

    # Check that output has correct type and value
    assert isinstance(output, float)
    assert output == expected_output

def test_read_input_invalid():
    """Test that trying to fetch a non-present input without a default raises an error."""
    # Make empty input dict
    input_dict = {}

    with pytest.raises(ValueError):
        # Make sure this raises an error (No input found and no default)
        output = utils.read_input(input_dict, int, "NOT PRESENT")

@pytest.mark.filterwarnings("ignore:No argument supplied")
@pytest.mark.parametrize("sign", ["","-"])
def test_read_layer(sign:str):
    """Test that read_layer function works as intended"""
    # Make input dict. Note that no input is supplied for 'Layer 3'
    input_dict = {
        f"Layer {sign}1 input" : "3",
        f"Layer {sign}2 input" : "4",
        f"Layer {sign}4 input" : "6"
    }

    N_layers = int(f"{sign}4") # 4 or -4

    # Read layers, give non-present layers (here: Layer 3) the value int(5)
    output = utils.read_input_layer(input_dict, N_layers, int, "input", default="5")
    expected_output = np.array([3,4,5,6])

    # Check that all entries match
    assert (output == expected_output).all()
