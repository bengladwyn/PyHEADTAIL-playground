import pytest
from typing import Dict, Any, Optional, Iterable, Tuple, Callable
from numpy.typing import ArrayLike
from IW2D.interface import IW2DLayer, FlatIW2DInput, RoundIW2DInput, Eps1FromCallableConductivity, Eps1FromResistivity, Eps1FromTandelta, Mu1FromSusceptibility
from IW2D.utils import make_layer_arrays, linear_interpolation_from_file, py_complex
from IW2D.flat_impedance import _iw2d_flat_impedance_single_frequency, iw2d_flat_impedance
from IW2D.round_impedance import _iw2d_round_impedance_single_frequency, iw2d_round_impedance
import numpy as np
import pandas as pd
from pathlib import Path

FREQUENCY_DEPENDENT_DATA_PATH = Path(__file__).parents[2] / "examples" / "input_data_files"

# TODO: consider moving all kwarg dicts to fixture
# Dictionaries used to initialize IW2DLayer and IW2DInput objects:

EPS1_RESISTIVITY = Eps1FromResistivity(
    dc_resistivity=2.5e-8,
    resistivity_relaxation_time=0.027e-12,
    re_dielectric_constant=1
)

EPS1_TANDELTA = Eps1FromTandelta(
    re_dielectric_constant=1,
    dielectric_tan_delta = 0.12
)

EPS1_CONDUCTIVITY = Eps1FromCallableConductivity(
    re_dielectric_constant=1,
    frequency_dependent_conductivity=linear_interpolation_from_file(FREQUENCY_DEPENDENT_DATA_PATH / "conductivity_vs_freq_hBN_example.dat", skiprows=1)
)

EPS1_FILEINTERP = linear_interpolation_from_file(FREQUENCY_DEPENDENT_DATA_PATH / "relative_permittivity_vs_freq_hBN_example.dat", skiprows=1)

MU1_SUSCEPTIBILITY_ZERO = Mu1FromSusceptibility(
    magnetic_susceptibility=0,
    permeability_relaxation_frequency=np.inf
)

MU1_SUSCEPTIBILITY_NONZERO = Mu1FromSusceptibility(
    magnetic_susceptibility=0.5, # increased from default 0 to get more meaningful test of mu1
    permeability_relaxation_frequency=1e-5 # decreased from inf
)

MU1_FILEINTERP = linear_interpolation_from_file(FREQUENCY_DEPENDENT_DATA_PATH / "ferrite_8C11_muprime_musecond.dat", skiprows=1)

LAYER_KWARGS_STANDARD = dict(
    thickness=np.inf,
    eps1=EPS1_RESISTIVITY,
    mu1=MU1_SUSCEPTIBILITY_NONZERO
)

# A layer using tan(delta) for eps1 calculation
LAYER_TANDELTA = IW2DLayer(
    thickness=np.inf,
    eps1=EPS1_TANDELTA,
    mu1=MU1_SUSCEPTIBILITY_ZERO
)

# A layer with eps1 interpolated from a file
LAYER_FREQEPS1 = IW2DLayer(
    thickness=np.inf, 
    eps1=EPS1_FILEINTERP,
    mu1=MU1_SUSCEPTIBILITY_ZERO
)

# A layer with eps1 calculated from sigma1, which is interpolated from a file
LAYER_FREQSIGMA = IW2DLayer(
    thickness=np.inf, 
    eps1=EPS1_CONDUCTIVITY,
    mu1=MU1_SUSCEPTIBILITY_ZERO
)

# A layer with mu1 interpolated from a file
LAYER_FREQMU1 = IW2DLayer(
    thickness=np.inf,
    eps1 = EPS1_RESISTIVITY,
    mu1=MU1_FILEINTERP
)

# Standard input parameters for a flat chamber setup. Mimics examples/input_files/FlatChamberInputFile.txt
FLAT_INPUT_KWARGS_STANDARD = dict(
    machine="LHC",
    length=9.0,
    relativistic_gamma=6927.62871617,
    top_bottom_symmetry=False,
    top_half_gap=4e-3,
    bottom_half_gap=2e-3,
    top_layers= tuple([
        IW2DLayer(**LAYER_KWARGS_STANDARD)
    ]),
    bottom_layers=[],
    comment="",
    calculate_wake=False
)

# Standard input parameters for a round chamber setup. Mimics examples/input_files/RoundChamberInputFile.txt
ROUND_INPUT_KWARGS_STANDARD = dict(
    machine="LHC",
    length=1,
    relativistic_gamma=479.605064966,
    calculate_wake=False,
    layers=tuple([
        IW2DLayer(
            thickness=np.inf,
            eps1=Eps1FromResistivity(
                dc_resistivity=5e-6,
                resistivity_relaxation_time=4.2e-12,
                re_dielectric_constant=1
            ),
            mu1=Mu1FromSusceptibility(
                magnetic_susceptibility=0.,
                permeability_relaxation_frequency=np.inf
            )
        )
    ]),
    inner_layer_radius=0.002,
)


def modify_dict(initial_dict: Dict[Any, Any], remove_keys: Optional[Iterable[Any]] = None, add_dict: Optional[Dict[Any, Any]] = None) -> Dict[Any, Any]:
    """Modify a dictionary by removing entries with keys specified by remove_keys, 
    and then add entries from add_dict.

    Mainly intended to modify the initialization kwarg dictionaries.

    Copies the input dictionary, so initial_dict is left unchanged.

    ### Parameters

    initial_dict : Dict[Any, Any]
        The dict to be modified.
    remove_keys : Iterable[Any], optional
        The keys of all entries to be removed, by default None
    add_dict : Dict[Any, Any], optional
        A dictionary with new entries to the dictionary, by default None

    ### Returns

    Dict[Any, Any]
        A copy of initial_dict where remove_keys have been removed and add_dict has been added.
    """
    new_dict = initial_dict.copy()

    if remove_keys is not None:
        for key in remove_keys:
            new_dict.pop(key)
    
    if add_dict is not None:
        new_dict.update(add_dict)
    
    return new_dict

@pytest.mark.parametrize(["eps1_method", "expected_eps1"], 
    [
        [EPS1_RESISTIVITY, -121975.17927628598-7190041433808.935j],
        [EPS1_TANDELTA, 1-0.12j],
        [EPS1_FILEINTERP , 2.796302001293914-0.029310057387198146j],
        [EPS1_CONDUCTIVITY, 1-0.00022903219348931287j]
    ]
)
def test_eps1(eps1_method: Callable[[float], complex], expected_eps1: complex):
    """Test that eps1 method gives the correct value"""

    output_eps1 = eps1_method(1e5) # 1e5 Hz, arbitrarily chosen

    assert output_eps1 == pytest.approx(expected_eps1)


@pytest.mark.parametrize(["mu1_method", "expected_mu1"], 
    [
        [MU1_SUSCEPTIBILITY_NONZERO, 1-5.000000000000001e-11j],
        [MU1_FILEINTERP, 1000-159.361j],
    ]
)
def test_mu1(mu1_method: Callable[[float], complex], expected_mu1: complex):
    """Test that mu1 method gives the correct value"""

    output_mu1 = mu1_method(1e5) # 1e5 Hz, arbitrarily chosen

    assert output_mu1 == pytest.approx(expected_mu1)


@pytest.mark.parametrize("b_negative", (False, True))
def test_make_layer_array(b_negative: bool):
    """Test if make_layer_arrays makes arrays with the correct size and values"""

    # make several layers
    input_layers = tuple([
        IW2DLayer(**modify_dict(LAYER_KWARGS_STANDARD, remove_keys=["thickness"], add_dict={"thickness":0.003})),
        LAYER_TANDELTA
    ])

    # make arrays
    output_eps1, output_mu1, output_b = make_layer_arrays(input_layers, b0=0.006, frequency=1e5, b_negative=b_negative)

    expected_eps1  = [1+0j, -121975.17927628598-7190041433808.935j, 1-0.12j]
    expected_mu1   = [1+0j, 1-5.000000000000001e-11j, 1+0j]
    if not b_negative:
        expected_b = [0.006, 0.009, np.inf]
    else:
        expected_b = [-0.006, -0.009, -np.inf]

    for p in range(len(input_layers)+1):
        assert py_complex(output_eps1[p+1]) == pytest.approx(expected_eps1[p])
        assert py_complex(output_mu1[p+1]) == pytest.approx(expected_mu1[p])
        assert output_b[p+1].toDouble() == pytest.approx(expected_b[p])


@pytest.mark.parametrize(["component", "expected_impedance"], [
    (0, 2.66666405e+00  + 2.66942000e+00j), # Zlong
    (1, 2.82696134e+01  + 2.82973329e+01j), # Zycst
    (2, 7.06442485e+03  + 7.07415538e+03j), # Zxdip
    (3, 7.06442485e+03  + 7.07415539e+03j), # Zydip
    (4, -7.06442485e+03 - 7.07415538e+03j), # Zxquad
    (5, 7.06442485e+03  + 7.07415539e+03j)  # Zyquad
])
def test_single_frequency_flat_impedance(component: int, expected_impedance: complex):
    """Test that iw2d_single_frequency_flat_impedance gives correct impedances (for a randomly selected frequency)"""

    input_frequency = 5.62341325e+08 # [Hz]
    input_obj = FlatIW2DInput(**FLAT_INPUT_KWARGS_STANDARD) # type: ignore

    output_impedances = _iw2d_flat_impedance_single_frequency(input_obj, input_frequency)
    output_impedance = output_impedances[component]

    assert output_impedance == pytest.approx(expected_impedance)


@pytest.mark.parametrize(["component", "expected_impedance"], [
    (0, 3.75369320e+01 + 5.97673564e+01j), # Zlong
    (1, 5.00925310e+04 + 7.92633035e+04j), # Zxdip
    (2, 5.00925310e+04 + 7.92633035e+04j), # Zydip
    (3, 3.04103463e-02 + 4.84202067e-02j), # Zxquad
    (4, 3.04103463e-02 + 4.84202067e-02j)  # Zyquad
])
def test_single_frequency_round_impedance(component: int, expected_impedance: complex):
    """Test that iw2d_single_frequency_flat_impedance gives correct impedances (for a randomly selected frequency)"""

    input_frequency = 1.77827941e+10 # [Hz]
    input_obj = RoundIW2DInput(**ROUND_INPUT_KWARGS_STANDARD) # type: ignore

    output_impedances = _iw2d_round_impedance_single_frequency(input_obj, input_frequency)
    output_impedance = output_impedances[component]

    assert output_impedance == pytest.approx(expected_impedance)


@pytest.mark.parametrize(("input_obj", "impedance_function"), (
    (FlatIW2DInput(**FLAT_INPUT_KWARGS_STANDARD), iw2d_flat_impedance), # type: ignore
    (RoundIW2DInput(**ROUND_INPUT_KWARGS_STANDARD), iw2d_round_impedance) # type: ignore
))
@pytest.mark.parametrize(("frequencies", "expected_length"), (
    (1e2,                       1),
    ([1e2,1e3,1e4],             3),
    (np.geomspace(1e2,1e4,10),  10),
    (pd.Series(data=[1e2,1e3]), 2)
))
def test_impedance_function_frequency_parsing(input_obj, impedance_function, frequencies: ArrayLike, expected_length: int):
    """iw2d_calculate_flat_impedance and iw2d_calculate_round_impedance work with both floats and iterables of floats as its 'frequencies' input"""
    
    output_dataframe, _ = impedance_function(input_obj, frequencies)
    
    assert len(output_dataframe) == expected_length


@pytest.mark.parametrize(
    ("input_obj", "impedance_function"),
    (
        (FlatIW2DInput(**FLAT_INPUT_KWARGS_STANDARD), iw2d_flat_impedance), # type: ignore
        (RoundIW2DInput(**ROUND_INPUT_KWARGS_STANDARD), iw2d_round_impedance) # type: ignore
    )
)
def test_complete_metadata(input_obj, impedance_function):
    """Each column in the impedance DataFrame should have a corresponding complete entry in the metadata dictionary"""
    frequencies = 1e5

    output_dataframe, output_metadata = impedance_function(input_obj, frequencies)

    assert set(output_dataframe.columns) == set(output_metadata.keys())

    for column_metadata in output_metadata.values():
        assert set(column_metadata.keys()) == {"Units","Plane","Exponents"}


@pytest.mark.parametrize(("layer", "expected_hash"),(
    (IW2DLayer(**LAYER_KWARGS_STANDARD),  "1828194f86d6dc3750a6de3f7fad4dd397c92d5a8049178abd6c2d74cd48ef62"),
    (LAYER_TANDELTA,  "c6850403a573aebec3e7a173d6063f961dd09257955a514a237f507c698fbdb0"),
    (LAYER_FREQEPS1,  "abc914e509d1a86fc452a94210cf397a9e8ca8a91b3c61ab5905bdd6feba8156"),
    (LAYER_FREQSIGMA, "eddccb1f33ddf63cc845c9ae93bf24a984697a5f3be9fb8a7d601fdca44291d7"),
    (LAYER_FREQMU1,   "8ffc864d84407a0b29ed63d8f4d3f38a18d7c3fbb8e151e791d90178bb3c62b2"),
))
def test_flatiw2dinput_hashes_reproducible(layer, expected_hash):
    """FlatIW2DInput should give reproducible hashes. Most important to check for layers with Callable's"""

    input_obj = FlatIW2DInput(**modify_dict(
        FLAT_INPUT_KWARGS_STANDARD,
        remove_keys=["top_layers"],
        add_dict={"top_layers":(layer,)}
    ))

    output_hash = input_obj.output_identification_hash()

    assert output_hash == expected_hash

@pytest.mark.parametrize(("layer", "expected_hash"),(
    (IW2DLayer(**LAYER_KWARGS_STANDARD),  "784f1cca52f2cec9fe8e7b38056bac9fb17f23b148a3a40714c3e3aea19672fa"),
    (LAYER_TANDELTA,  "c2072dc265412a3ace6e0d153b46e6826c30e95ed8aa663a8a9cd484b4b1ae58"),
    (LAYER_FREQEPS1,  "345fce112f2ea82f86cdecd800baaab8ef9bff382b55ed1b06a6f6069fce8d82"),
    (LAYER_FREQSIGMA, "218827f7fd229cefdc912df46155f15669055a52bf19680d9296df1c38d9a73e"),
    (LAYER_FREQMU1,   "ca39d08389a7a8a9b768a2fe398ee45e89f5e63ace99bc60d7e038f7ef3de836"),
))
def test_roundiw2dinput_hashes_reproducible(layer, expected_hash):
    """RoundIW2DInput should give reproducible hashes. Most important to check for layers with Callable's"""

    input_obj = RoundIW2DInput(**modify_dict(
        ROUND_INPUT_KWARGS_STANDARD,
        remove_keys=["layers"],
        add_dict={"layers":(layer,)}
    ))

    output_hash = input_obj.output_identification_hash()

    assert output_hash == expected_hash
