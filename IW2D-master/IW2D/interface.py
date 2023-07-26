from dataclasses import dataclass, field
from typing import Optional, Tuple, Callable, Iterable
import numpy as np
from scipy.constants import epsilon_0
from hashlib import sha256

@dataclass(frozen=True, eq=True)
class Eps1FromResistivity:
    """A function object for calculating the relative complex permittivity from the resistivity.

    :param dc_resistivity: The DC resistivity rho (Ohm*m). Can be infinity, e.g. for vacuum
    :type dc_resistivity: float
    :param resistivity_relaxation_time: The relaxation time tau (s) for the AC resistivity. 0 for vacuum.
    :type resistivity_relaxation_time: float
    :param re_dielectric_constant: The real part of the relative permittivity, epsb (unitless). 1 for vacuum or metals.
    :type re_dielectric_constant: float
    """

    dc_resistivity: float # [Ohm*m]
    resistivity_relaxation_time: float # [s]
    re_dielectric_constant: float # [1]

    def __post_init__(self):
        """Runs automatically as part of __init__, checks whether inputs are in valid range"""
        assert np.isfinite(self.re_dielectric_constant), "re_dielectric_constant must be finite number"
        assert (self.dc_resistivity > 0), "dc_resistivity must be positive (non-zero) number"


    def __call__(self, frequency: float) -> complex:
        """Calculate the relative complex permittivity.

        :param frequency: The frequency f in Hz
        :type frequency: float
        :return: The relative complex permittivity eps1 given by: eps1 = epsb+1/(1j*epsilon_0*rho*2*pi*f*(1+1j*2*pi*f*tau))
        :rtype: complex
        """

        rho = self.dc_resistivity
        tau = self.resistivity_relaxation_time
        epsb = self.re_dielectric_constant
        omega = 2*np.pi*frequency

        # this check isn't necessary for calculation with python float, but is necessary if we go back to multiprecision
        if np.isfinite(rho) and np.isfinite(tau):
            return epsb+1/(1j*epsilon_0*rho*omega*(1+1j*omega*tau))
        else:
            return epsb


@dataclass(frozen=True, eq=True)
class Eps1FromCallableConductivity:
    """A function object for calculating the relative complex permittivity from the frequency dependent conductivity.

    :param re_dielectric_constant: The real part of the relative permittivity, epsb (unitless). 1 for vacuum or metals.
    :type re_dielectric_constant: float
    :param frequency_dependent_conductivity: A function that takes one argument: the frequency in Hz, and returns the conductivity sigma in S/m
    :type frequency_dependent_conductivity: callable
    """

    re_dielectric_constant: float # [1]
    frequency_dependent_conductivity: Callable[[float], complex] # [Hz] -> [S/m]

    def __post_init__(self):
        """Runs automatically as part of __init__, checks whether inputs are in valid range"""
        assert np.isfinite(self.re_dielectric_constant), "re_dielectric_constant must be finite number"

    def __call__(self, frequency: float) -> complex:
        """Calculate the relative complex permittivity.

        :param frequency: The frequency f in Hz
        :type frequency: float
        :return: The relative complex permittivity eps1 given by: eps1 = epsb + sigma/(1j*epsilon_0*2*pi*f)
        :rtype: complex
        """

        epsb = self.re_dielectric_constant
        sigma = self.frequency_dependent_conductivity(frequency)
        omega = 2*np.pi*frequency
        return epsb + sigma/(1j*epsilon_0*omega)


@dataclass(frozen=True, eq=True)
class Eps1FromTandelta:
    """A function object for calculating the relative complex permittivity from the materials dielectric tan(delta)

    :param re_dielectric_constant: The real part of the relative permittivity, epsb (unitless)
    :type re_dielectric_constant: float
    :param dielectric_tan_delta: The dielectric tan(delta) (unitless)
    :type dielectric_tan_delta: float
    """

    re_dielectric_constant: float # [1]
    dielectric_tan_delta: float # [1]

    def __post_init__(self):
        """Runs automatically as part of __init__, checks whether inputs are in valid range"""
        assert np.isfinite(self.re_dielectric_constant), "re_dielectric_constant must be finite number"
        assert np.isfinite(self.dielectric_tan_delta), "dielectric_tan_delta must be finite number"

    def __call__(self, frequency: float) -> complex:
        """Calculate the relative complex permittivity.

        :param frequency: The frequency f in Hz (not used)
        :type frequency: float
        :return: The relative complex permittivity eps1 given by: eps1 = epsb*(1-1j*tan(delta))
        :rtype: complex
        """

        epsb = self.re_dielectric_constant
        tandelta = self.dielectric_tan_delta
        return epsb*(1-1j*tandelta)


@dataclass(frozen=True, eq=True)
class Mu1FromSusceptibility:
    """Function object for calculating the relative complex permeability from the material's magnetic susceptibility.

    :param magnetic_susceptibility: The magnetic susceptibility chi (unitless), equal to mu_r-1 if the permeability is real. 0 for vacuum.
    :type magnetic_susceptibility: float
    :param permeability_relaxation_frequency: The cut-off frequency f_mu (Hz) for the frequency-dependent complex permeability. Usually Infinity (vacuum, metals, and others).

    :type permeability_relaxation_frequency: float
    """

    magnetic_susceptibility: float # [1]
    permeability_relaxation_frequency: float # [Hz]

    def __post_init__(self):
        """Runs automatically as part of __init__, checks whether inputs are in valid range"""
        assert np.isfinite(self.magnetic_susceptibility), "magnetic_susceptibility must be finite number"
        assert self.permeability_relaxation_frequency > 0, "permeability_relaxation_frequency must be greater than 0"


    def __call__(self, frequency: float) -> complex:
        """Calculate the relative complex permeability.

        :param frequency: The frequency f in Hz
        :type frequency: float
        :return: The relative complex permeability mu1 given by: mu1 = 1+chi/(1+1j*f/f_mu)
        :rtype: complex
        """

        chi = self.magnetic_susceptibility
        f_mu = self.permeability_relaxation_frequency
        return 1+chi/(1+1j*frequency/f_mu)


@dataclass(frozen=True, eq=True)
class IW2DLayer:
    """Defines a a chamber wall layer, with all material and geometric properties used by IW2D calculations.

    A note on performance: The supplied functions eps1 and mu1 get called once for each layer for each frequency.
    As long as one call of these functions takes about 1 ms or less, it should not affect performance considerably.

    :param thickness: The layer thickness in m.
    :type thickness: float
    :param eps1: A callable that takes a frequency in Hz as its only argument, and returns the relative complex permittivity of the material for that frequency.
    :type eps1: callable
    :param mu1: A callable that takes a frequency in Hz as its only argument, and returns the relative complex permeability of the material for that frequency.
    :type mu1: callable
    :param sample_frequencies_for_output_repr: The frequencies in Hz at which the outputs of eps1 and mu1 will be evaluated when generating an output repr string.
    Defaults to tuple([10**i for i in range(2,17)])
    :type sample_frequencies_for_output_repr: Iterable[float]
    """

    thickness: float # [m]
    eps1: Callable[[float], complex] # [Hz] -> [1]
    mu1: Callable[[float], complex] # [Hz] -> [1]
    sample_frequencies_for_output_repr: Iterable[float] = tuple([10**i for i in range(2,17)]) # [Hz]

    def output_repr(self) -> str:
        """Generate a string with the thickness and the outputs of eps1, mu1 from the layer,
        where eps1 and mu1 are evaluated at every frequency in self.sample_frequencies_for_output_repr.

        This function is mainly meant for generating hashes from IW2DInput objects in a consistent manner,
        that works reproducibly for any function that might be given as eps1 and mu1

        :return: A reproducible string generated from the outputs of the layers's thickness, eps1, mu1
        :rtype: str
        """

        esp1_output_str = "".join([str(self.eps1(frequency)) for frequency in self.sample_frequencies_for_output_repr])
        mu1_output_str  = "".join([str(self.mu1(frequency))  for frequency in self.sample_frequencies_for_output_repr])

        return esp1_output_str + mu1_output_str + str(self.thickness)


# Define several dataclasses for IW2D input elements. We must split mandatory
# and optional arguments into private dataclasses to respect the resolution
# order. The public classes RoundIW2DInput and FlatIW2D input inherit from
# from the private classes.
# https://stackoverflow.com/questions/51575931/class-inheritance-in-python-3-7-dataclasses

@dataclass(frozen=True, eq=True)
class _IW2DInputBase:
    length: float
    relativistic_gamma: float
    calculate_wake: bool
    # f_params: Sampling


@dataclass(frozen=True, eq=True)
class _IW2DInputOptional:
    # z_params: Optional[Sampling] = None
    machine: str = ""
    long_factor: Optional[float] = None
    wake_tol: Optional[float] = None
    freq_lin_bisect: Optional[float] = None
    comment: str = ""


@dataclass(frozen=True, eq=True)
class _RoundIW2DInputBase(_IW2DInputBase):
    layers: Tuple[IW2DLayer, ...] = field(repr=False) # remove repr for hashing
    inner_layer_radius: float


@dataclass(frozen=True, eq=True)
class _RoundIW2DInputOptional(_IW2DInputOptional):
    yokoya_Zlong: float = 1.
    yokoya_Zxdip: float = 1.
    yokoya_Zydip: float = 1.
    yokoya_Zxquad: float = 0.
    yokoya_Zyquad: float = 0.


@dataclass(frozen=True, eq=True)
class RoundIW2DInput(_RoundIW2DInputOptional, _RoundIW2DInputBase):
    
    def output_identification_hash(self) -> str:
        """Generate a sha256 hash of the object, intended to be used to uniquely identify the output it generates.

        The hash uses the repr of the input object itself (note that top_layers and bottom_layers are excluded from the repr)
        and a string generated from the output eps1 and mu1 methods and the thickness of the layers.

        :return: A sha256 encoded hash from the objects repr and the outputs of the layers.
        :rtype: str
        """

        repr_with_layers = repr(self) + "".join([layer.output_repr() for layer in self.layers])
        return sha256(repr_with_layers.encode()).hexdigest()


@dataclass(frozen=True, eq=True)
class _FlatIW2DInputBase(_IW2DInputBase):
    top_bottom_symmetry: bool
    top_layers: Tuple[IW2DLayer, ...] = field(repr=False) # remove repr for hashing
    top_half_gap: float


@dataclass(frozen=True, eq=True)
class _FlatIW2DInputOptional(_IW2DInputOptional):
    bottom_layers: Optional[Tuple[IW2DLayer, ...]] = field(default=None, repr=False) # remove repr for hashing
    bottom_half_gap: Optional[float] = None


@dataclass(frozen=True, eq=True)
class FlatIW2DInput(_FlatIW2DInputOptional, _FlatIW2DInputBase):
    
    def output_identification_hash(self) -> str:
        """Generate a sha256 hash of the object, intended to be used to uniquely identify the output it generates.

        The hash uses the repr of the input object itself (note that top_layers and bottom_layers are excluded from the repr)
        and a string generated from the output eps1 and mu1 methods and the thickness of the layers.

        :return: A sha256 encoded hash from the objects repr and the outputs of the layers.
        :rtype: str
        """

        if self.bottom_layers is None:
            layers = self.top_layers
        else:
            layers = tuple(self.top_layers) + tuple(self.bottom_layers) # concatenate the tuples

        repr_with_layers = repr(self) + "".join([layer.output_repr() for layer in layers])
        return sha256(repr_with_layers.encode()).hexdigest()


def verify_wake_possible():
    """Should make a function to verify that when wake==True, all required optional parameters are passed"""
    # TODO: implement when wake is implemented
    pass
