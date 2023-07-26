from typing import List, Union
import IW2D

from cppyy.gbl import amp, precision

PRECISION = precision

def make_underline(left_campf: amp.campf[PRECISION], right_campf: amp.campf[PRECISION]) -> str:
    """Make a string to underline and thus highlight where the (decimal) digits of the complex number right_campf does not match those of left_campf."""

    underline_symbol = "^"

    # Generate python strings of the real and imag components of left and right
    left_real = str(left_campf.x.toDec())
    left_imag = str(left_campf.y.toDec())
    right_real = str(right_campf.x.toDec())
    right_imag = str(right_campf.y.toDec())

    for str_repr in [left_real, left_imag, right_real, right_imag]:
        if "@" in str_repr:
            # If NaN or infinity: skip making underline
            return ""

    # Split into mantissa and exponent
    left_real_mantissa, left_real_exponent = left_real.split("E")
    left_imag_mantissa, left_imag_exponent = left_imag.split("E")
    right_real_mantissa, right_real_exponent = right_real.split("E")
    right_imag_mantissa, right_imag_exponent = right_imag.split("E")

    # Check for mantissa mismatch
    mantissa_mismatch_digits = []

    for left, right in zip(
        [left_real_mantissa, left_imag_mantissa],
        [right_real_mantissa, right_imag_mantissa]):
        # First check real, then imag

        # Set default 0 mismatch
        mismatch_digits = 0

        for i in range(len(right)):
            # Loop over characters in string

            if left[i] != right[i]:
                # If mismatch, set that number of mismatching digits and exit loop.
                mismatch_digits = len(right)-i
                break
        
        mantissa_mismatch_digits.append(mismatch_digits)


    
    # Check for exponent mismatch
    exponent_mismatch_real = (left_real_exponent != right_real_exponent)*len(right_real_exponent)
    exponent_mismatch_imag = (left_imag_exponent != right_imag_exponent)*len(right_imag_exponent)


    # Make underline string
    underline = (
        (underline_symbol*mantissa_mismatch_digits[0]).rjust(len(right_real_mantissa)," ") + # Underline for real mantissa
        (underline_symbol*exponent_mismatch_real).rjust(len(right_real_exponent)+1, " ") +   # Underline for real exponent. +1 to compensate for E symbol
        "     " +                                                                            # 5 spaces to compensate for "  +  " in the numbers themselves
        (underline_symbol*mantissa_mismatch_digits[1]).rjust(len(right_imag_mantissa)," ") + # Underline for imag component
        (underline_symbol*exponent_mismatch_imag).rjust(len(right_imag_exponent)+1, " ")     # Underline for imag exponent
    ) 

    
    return underline





def pytest_assertrepr_compare(op, left, right) -> Union[List[str], None]:
    """Implement custom error messages to failed comparison assertion. Pytest calls this function automatically."""

    # Equality comparison of two amp::campfs
    if isinstance(left, amp.campf[PRECISION]) and isinstance(right, amp.campf[PRECISION]) and op == "==":
        return ["equality of two amp::campf instances",
            f"{left.x.toDec()}  +  {left.y.toDec()} j !=",
            f"{right.x.toDec()}  +  {right.y.toDec()} j",
            make_underline(left, right),
            "Comparison done on the numbers in binary representation. Conversion to decimal only done in error message."
        ]
    
    return None
