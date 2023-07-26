#  ImpedanceWake2D

Code to compute the longitudinal and transverse beam coupling
impedances and wake functions of a portion of an infinitely long multilayer structure, either cylindrical
or flat and infinitely large. The main purpose is to evaluate the so-called wall (or resistive-wall) impedance.

Authors: Nicolas Mounet, Nicolò Biancacci, David Amorim, Eskil Vik, CERN, BE dpt, Geneva

Last modifications: 24/02/2023

Comments / remarks are welcome and can be sent to nicolas.mounet@cern.ch or nicolas.mounet@m4x.org.

This program is free software; you can redistribute it and/or
modify it under the terms version 3 of the GNU General Public License
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


## 1. INTRODUCTION

ImpedanceWake2D is a package containing four different codes that compute the longitudinal and transverse
beam coupling impedances and wake functions in a multilayer axisymmetric or flat structure that is two dimensional
(i.e. of infinite length, that is without endings in the longitudinal direction). The impedance is given
for a portion of the structure of given length. The number of layers
can be in principle anything, and each of them can be made of any linear material. One can use either 
general-purpose frequency-dependent expressions for the complex permittivity and permeability, depending
on 5 parameters (see below), or fully general frequency-dependent quantities defined through interpolation tables.
The last layer has an infinite thickness (hence there is no perfectly conducting layer in the end).
The code relies on the analytic computation of the electromagnetic fields created by a point-charge beam travelling
at any speed (not necessarily ultrarelativistic) in the whole structure.

For the impedances of a cylindrical structure, the original ideas are from B. Zotter [1-3]. They were previously 
implemented in code called LAWAT [2], later [4,5] converted to a Mathematica [6] notebook. A full 
review of the formalism is in Ref. [7]. It was also implemented by B. Salvant (CERN, BE dpt) in a Mathematica
notebook called ReWall.

For the impedances of a flat structure, the formalism was developped during N. Mounet's PhD thesis [8], and is fully 
described in Ref. [9].

The theory on both the cylindrical and flat cases is shortly described in Ref. [10], and extensively in Ref. [8].

To compute the wake functions (time domain counterparts of the impedances), an algorithm based on a Filon's kind 
of method is used, instead of an FFT. The main advantages are the accuracy and the fact that an uneven frequency 
sampling can be used, which is useful in the case of impedances because they exhibit different features for 
different frequency ranges, i.e. their detailed behaviour over a large number of decades is needed 
to compute accurately the wakes. The codes automatically refine the frequency sampling until a 
prescribed accuracy on the wake function is reached. The full method is also described in Ref. [8].
 

IMPORTANT NOTE: the impedance computed here is the so-called "wall impedance", which is sligthly different from the
"resistive-wall impedance": the wall impedance contains the indirect space charge term (perfect conducor impedance)
whereas the resistive-wall one doesn't. This is because the indirect space-charge term is crucial for the 
low-frequency behaviour (we cannot easily separate its effect from the resistive part). Details on this concept 
can be found in Ref. [11]. The same applies for the wake functions: they also include the indirect space charge. 

NOTE ON PARTICULAR USES: 
- the codes work for any relativistic mass factor gamma (or equivalently for any relativistic 
velocity factor beta). If one wants to use the code in the ultrarelativistic limit, the only possibility is to put a 
very large gamma in the input file. Note anyway that above a certain frequency, the computation will be probably different
from the ultrarelativistic case (this frequency limit certainly increases with gamma).
- there is no perfectly conducting layer (PEC) on the outside. If one wants to put one, one possibility is to put a very
low resistivity on the last layer (e.g. 10^-15 Ohm.m), but some numerical issues could occur, in which case it is advised
to increase the precision (see below).


## 2. DESCRIPTION OF THE PACKAGE

The repository is structured as follows:

- `IW2D/` : Python and c++ source code
    - `External_libs/` : Multiprecision external libraries
        - `gmp-6.0.0a.tar.gz` and `mpfr-3.1.2.tar.gz` : compressed source code for the libraries GMP and MPFR
        - `compile_ext_libs.sh` bash script for extracting and compiling the libraries
        - Note: Use is optional, as system or conda installation also works, see Installation section below
    - `legacy/` : python package to calculate impedances and wake field with the exact same input and output files as the C++ executables
        - `flatchamber.py`, `roundchamber.py`, `wake_flatchamber.py`, `wake_roundchamber.py` : Python scripts that compute resp. the impedances of a cylindrical chamber, the impedances and wake functions of a cylindrical chamber, the impedances of a flat chamber, and the impedances and wake functions of a flat chamber
        - `utils.py` : Python methods used in the scripts above
        - `__main__.py` : Python script called from command line to read inputs and execute the scripts above
    - `cpp/` : (Mainly) C++ source code and Makefiles to compile them
        - `roundchamber.x`, `wake_roundchamber.x`, `flatchamber.x` and `wake_flatchamber.x` : executables that compute resp. the impedances of a cylindrical chamber, the impedances and wake functions of a cylindrical chamber, the impedances of a flat chamber, and the impedances and wake functions of a flat chamber.
        - `plot_Imp_flat.py` and `plot_Wake_flat.py` : python routines to plot the impedances and wake functions in output (they work for both round and flat cases).
        - `plot_Imp_all_terms_in_one_plot.py` and `plot_Imp_vs_displacement.py` : python plotting routines specifically for the non-linear case (now implemented only for a flat chamber)
- `PYTHON_codes_and_scripts/` : Old python scripts, currently not supported
- `examples/` : example files
    - `input_files/` : example input files for executables/legacy scripts. Note that these are used as inputs to some tests.
    - `outputs/` : folders containing the expected outputs for the input files in `examples/input_files`, to use for testing/verification.
    - `input_data_files/` : examples of data files for frequency dependent material properties
- `Tests/` : Unit tests and full package tests
    - `Pytests/` : `pytest` unit tests of python and C++ functions (using `cppyy` for the latter)
    - `run_and_check_IW2D.sh` : Script that runs the C++ executables on the input files in `IW2D_examples/` and check that they produce the expected output.
    - `verify_backwards_python.py` : Script for comparing results of C++ executables and python scripts


## 3. INSTALLATION
 
The following instructions should be valid for a Linux operating system. They were checked under CERN CentOS 7, CERN Alma9 and under Ubuntu 22.04.

If experiencing other problems, refer to Section 7 - Troubleshooting below.

For other operating systems and/or compilers, the user is welcome to share his/her installation experience by sending an email to nicolas.mounet@cern.ch .

### 3.1 Installation of prerequisite libraries

The code requires the GSL (http://www.gnu.org/software/gsl/), MPFR (https://www.mpfr.org/), Arb (https://arblib.org/index.html) and cppyy (https://cppyy.readthedocs.io/en/latest/) libraries to run.

The simplest installation procedure for these libraries is to use `conda`, which can be installed from one of the following:

- Miniconda (https://docs.conda.io/en/latest/miniconda.html) - minimal installer for `conda`
- Anaconda (https://www.anaconda.com/products/distribution) - installer for `conda` with several other features

With conda installed, run the following commands to create and activate a new environment (here called "my_env"):

    conda create -n my_env python=3.8
    conda activate my_env

The option `python=3.8` sets the python version used by the environment to 3.8. It can be changed to another python version (>=3.7), or removed to use the system installation of python.

Then, install the GSL, MPFR, Arb and cppyy libraries from the conda-forge channel into the new environment:

    conda install -c conda-forge gsl mpfr arb cppyy

---

To install without using conda, one can use a system package manager such as apt-get to install GSL, MPFR and Arb to the system. Under Ubuntu:

    sudo apt-get install libgsl-dev libmpfr-dev libflint-dev libflint-arb-dev


while cppyy can be installed using pip:
    
    pip install cppyy

---

MPFR can also be compiled locally, using a custom install script:

    cd IW2D/External_libs  # either where you cloned the repository or in the pip installation folder
    ./compile_ext_libs.sh

### 3.1.1 A note about Arb

The Arb library may have different name on other distributions. While it on Ubuntu and Debian is called `libflint-arb`, on most other systems it is called `libarb`. More information is available at: https://arblib.org/setup.html

**PLEASE NOTE:** If arb is insalled under this name, the environment variable `IW2D_FLINT_ARB` must be set (to anything) for the Python code to work. Otherwise, the program attempts to load `libarb` by default. On linux systems, this can be done by running the following in the terminal:

    echo 'export IW2D_FLINT_ARB=1' >> ~/.bashrc && source ~/.bashrc

Alternatively, to define the environment variable only for the current shell session, run `export IW2D_FLINT_ARB=1`.

The dependency Flint (libflint-dev) should be installed with the same installer as the one used to install Arb.


### 3.2 Installation of python package

To install the python package, simply use

    pip install git+ssh://git@gitlab.cern.ch:7999/IRIS/IW2D.git

or, if you have cloned the repository locally,

    pip install -e path/to/repository/root_dir

where the `-e` option is optional and makes any changes you make to the source code immediately affect the installed package, without having to reinstall.

To also install testing requirements, run:
    
    pip install -e path/to/repository/root_dir[test]

### 3.3 Compilation of C++ code

The python package is designed as a complete replacement of the executables generated directly from the C++ files (`.x` files). If one still wishes or needs the executables, they can be generated by doing the following:

- If MPFR was compiled locally (see above):

        cd IW2D/cpp
        cp ./Makefile_local_GMP_MPFR Makefile
        make

- If MPFR is installed to system or conda environment:

        cd IW2D/cpp
        cp ./Makefile_system_GMP_MPFR Makefile
        make

    If you moved the `mpfr-3.1.2` and/or the `gmp-6.0.0` directory to another location (i.e. not `IW2D/External_libs/`), make sure to update the PATH_MPFR and/or PATH_GMP variables accordingly.

The Arb library might be installed under the name `libflint-arb.so`, in which case `-larb` should be switched for `-lflint-arb` in the Makefiles.

To remove all `.o` and `.x` files, run (from the `IW2D/cpp` directory):

    make clean

This procedure requires that the `g++` compiler is installed on the system. To use another compiler, simply change the `CC` variable in the `Makefile`.

### 3.4 Testing the library

The following script tests that the C++ executables (`.x` files) produce the expected output.

    cd Tests
    ./run_and_check_IW2D.sh
       
## 4. USING THE CODE



###   4.1 Calculation from input files

To calculate impedances or wake functions from input files (described below) run one of the following:

- To use the `IW2D` python package, run e.g.:

        python3 -m IW2D.legacy flatchamber ./FlatChamberInputFile.txt

    where `flatchamber` can be replaced by `roundchamber`, `wake_flatchamber` or `wake_roundchamber`, and `./FlatChamberInputFile.txt` is the input file.

- To use the C++ executables, run e.g.:

        ./flatchamber.x < ./FlatChamberInputFile.txt
    
    where `flatchamber.x` can be replaced with any of the other executables, and `./FlatChamberInputFile.txt` is again the input file.


NOTE: the wake computations codes take much more time than the impedance ones, in general.
       

###   4.2 How to use the Python libraries (currently not supported)
       
The library "libIW2D.so" is used in the Python programs located in 
../PYTHON_codes_and_scripts/General_Python_tools/fourier_lib.py
In order for those to work, you need to add a line to your .bashrc file to add the IW2D path
to your LD_LIBRARY_PATH.
Also, for the Python Impedance library located in
../PYTHON_codes_and_scripts/Impedance_lib_Python/Impedance.py
you need to affect the environment variable IW2D_PATH with the path to the current directory
(containing all executables needed by Impedance.py)

All this can be simply done by running once the script given in this directory: type
       
    ./script_configure_bashrc.sh
       
and then logout and log in again, or do

    source ~/.bashrc
            
- General_Python_tools: to be installed first,
    * if you want to use the Fourier library (fourier_lib.py), you need first to have
    installed ImpedanceWake2D (with its library) (see above)
    DO NOT FORGET TO PUT IN YOUR .bashrc the corresponding LD_LIBRARY_PATH
    (use ../ImpedanceWake2D/script_configure_bashrc.sh)
    
    * you need to add the path to this directory to your PYTHONPATH. The best
    is do to it in your .bashrc file. For this you can type in a terminal

        cd PYTHON_codes_and_scripts/General_Python_tools
        ./script_configure_bashrc.sh
    
- Impedance_lib_Python: Impedance python library:
    * you need to have installed ImpedanceWake2D first, as well as General_Python_tools (see above)
    * to be able to use the function 'imp_model_from_IW2D' with lxplusbatch='launch' or 'retrieve'
    (i.e. to launch parallel jobs on a cluster - typically lxplus at CERN), you need
    to have the HTContor queuing system installed (this is the case on lxplus at CERN),
    with typical commands "condor_q", "condor_submit", etc.
    Otherwise you can still use this routine but sequentially, with lxplusbatch=None.
    
    * you need to add the path to this directory to your PYTHONPATH, and define
    the environement variable YOK_PATH for the path to the Yokoya factors file.
    The simplest is to modify your .bashrc file, typing in a terminal
    
        cd PYTHON_codes_and_scripts/Impedance_lib_Python
        ./script_configure_bashrc.sh


## 5. THE INPUT FILES

IMPORTANT NOTE: there must ALWAYS be a tab between the parameter description and its value, and more generally the exact 
sentence of the parameter description should be kept identical.

If there are more layers described in the input file than indicated at the 
beginning of it, the additional ones are ignored.

Examples of input files are both found below, and in the `examples/input_files/` directory of the repository.

###   5.1 Input file for the round chamber impedance computation
       
- Input file example:

        Machine:    LHC
        Relativistic Gamma: 479.6
        Impedance Length in m:  1
        Number of layers:   2
        Layer 1 inner radius in mm: 4
        Layer 1 DC resistivity (Ohm.m): 5.4e-8
        Layer 1 relaxation time for resistivity (ps):   0.005
        Layer 1 real part of dielectric constant:   1
        Layer 1 magnetic susceptibility:    0
        Layer 1 relaxation frequency of permeability (MHz): Infinity
        Layer 1 thickness in mm:    25
        Layer 2 DC resistivity (Ohm.m): 7.2e-7
        Layer 2 relaxation time for resistivity (ps):   0
        Layer 2 real part of dielectric constant:   1
        Layer 2 magnetic susceptibility:    0
        Layer 2 relaxation frequency of permeability (MHz): Infinity
        Layer 2 thickness in mm:    Infinity
        start frequency exponent (10^) in Hz:   0
        stop  frequency exponent (10^) in Hz:   15
        linear (1) or logarithmic (0) or both (2) frequency scan:   0
        sampling frequency exponent (10^) in Hz (for linear):   8
        Number of points per decade (for log):  20
        when both, fmin of the refinement (in THz): 1
        when both, fmax of the refinement (in THz): 10
        when both, number of points in the refinement:  500
        added frequencies [Hz]: 1e-5 1e16
        Yokoya factors long, xdip, ydip, xquad, yquad:  1 0.411233516712057 0.822467033424113 -0.411233516712057 0.411233516712057
        Comments for the output files names:    _some_element


- Description of each parameter:

    - Machine: name of the machine in which the element is located. It's used
    only in the output files names.
    - Relativistic Gamma: gamma=1/sqrt(1-beta^2) for the point-charge creating the
    electromagnetic fields and impedances.
    - Impedance Length in m: length (along the beam orbit) of the element
    considered.
    - Number of layers: number of layers in the multilayer cylindrical structure
    considered. For example 1 layer means that around the vacuum region only
    one layer of material is considered, with infinite thickness. Typically,
    a pipe wall made of a coated material and surrounded by vacuum will have
    three layers (coating + main wall material + vacuum). Any number of layers can be treated.
    - Layer 1 inner radius in mm: inner radius of the structure.
    - Layer n DC resistivity (Ohm.m):  DC resistivity rho of the layer considered
    (can be Infinity. e.g. for vacuum).
    - Layer n relaxation time for resistivity (ps): relaxation time tau in
    picoseconds (BEWARE OF THE UNITS), for the AC conductivity (or
    resistivity). 0 for vacuum.
    - Layer n real part of dielectric constant: eps_b (1 for vacuum or metals).
    - Layer n magnetic susceptibility: it is chi_m=mu_r-1 if the permeability is
    real. 0 for vacuum.
    - Layer n relaxation frequency of permeability (MHz): f_mu, cut-off frequency for the
    frequency-dependent complex permeability. Usually Inifinity (vacuum,
    metals, and others).
    - Layer n thickness in mm: thickness of the layer n. Can be Infinity (actually it
    has to be Infinity for the last layer).

    The 6 previous items are repeated for all the layers, changing the number
    n from 1 to "Number of layers" (specified at the beginning). There can be 
    more layers than needed (i.e. than "Number of layers") in the input file, 
    it is not a problem.
    In summary, the total complex permittivity in each layer is:
    
        eps_c = eps_0 * eps_b + 1/(j*rho*2*pi*f*(1+j*2*pi*f*tau))
    
    and the permeability is:
    
        mu = mu_0 * (1 + chi_m/(1+j*f/f_mu) )
        
    where f is the frequency, eps_0 the vacuum permittivity and mu_0 its
    permeability.

    - start frequency exponent (10^) in Hz: fmin: the scan (see flag below) in 
    frequency for the impedance calculation begins at 10^fmin Hz
    - stop  frequency exponent (10^) in Hz: fmax: the scan (see flag below) in 
    frequency for the impedance calculation ends at 10^fmax Hz
    - linear (1) or logarithmic (0) or both (2) frequency scan: flag for the 
    frequency sampling to be used: if it's 1, a linear (evenly spaced) scan is performed,
    if it's 0 it's a logarithmically spaced scan, if it's 2, we use
    both (a logarithmic scan plus a linear refinement, usually to sample
    accurately the high-frequency peak).
    - sampling frequency exponent (10^) in Hz (for linear): when option 1 is
    selected above, the scan goes from 10^fmin to 10^fmax with a spacing of
    10^"sampling frequency exponent".
    - Number of points per decade (for log): when option 0 or 2 is selected
    above, the scan goes from 10^fmin to 10^fmax with this number of points
    per decades.
    - when both, fmin of the refinement (in THz):  fminref
    - when both, fmax of the refinement (in THz):  fmaxref
    - when both, number of points in the refinement:   nref
    When the option above is 2, the last three parameters define an additional
    linear scan of nref points from fminref to fmaxref (BEWARE of UNITS:
    THZ). It is added to the log scan previously defined, sorting all the 
    frequencies in ascending order (also, if some frequencies
    appear twice, the redundant ones are erased).
    - added frequencies [Hz]:  several individual frequencies (separated by 
    spaces in the input file) to be added to all the previously defined scan.
    It can be empty.
    - Yokoya factors long, xdip, ydip, xquad, yquad: 5 values (separated by
    spaces) giving the factors to apply to get the final impedances and wakes.
    Basically, the code computes the longitudinal and transverse dipolar impedance
    of a cylindrical geometry, and apply those factors to get the long., dip. and quad. impedance
    for another geometry. First value is for the longitudinal term, second one for the dipolar x
    transverse term, third one for the dipolar y transverse term,
    fourth one for the quadrupolar x transverse term, and fifth one for the
    quadrupolar y transverse terms. NOTE: in the cylindrical case (Yokoya factors="1 1 1 0 0"), we take into
    account the non-ultrarelativistic axisymmetric quadrupolar impedance (new
    quadrupolar impedance, see [7]), such that the quadrupolar impedance is not zero.
    - Comments for the output files names: additional string
    that will be added at the end of the output file name. It can be empty.


###   5.2 Input file for the round chamber wake computation
       
- Input file example:

        Machine:    LHC
        Relativistic Gamma: 7460.52
        Impedance Length in m:  1
        Number of layers:   2
        Layer 1 inner radius in mm: 18.375
        Layer 1 DC resistivity (Ohm.m): 7.7e-10
        Layer 1 relaxation time for resistivity (ps):   0.5
        Layer 1 real part of dielectric constant:   1
        Layer 1 magnetic susceptibility:    0
        Layer 1 relaxation frequency of permeability (MHz): Infinity
        Layer 1 thickness in mm:    0.05
        Layer 2 DC resistivity (Ohm.m): 6.e-7
        Layer 2 relaxation time for resistivity (ps):   0
        Layer 2 real part of dielectric constant:   1
        Layer 2 magnetic susceptibility:    0
        Layer 2 relaxation frequency of permeability (MHz): Infinity
        Layer 2 thickness in mm:    Infinity
        linear (1) or logarithmic (0) or both (2) scan in z for the wake:   2
        sampling distance in m for the linear sampling: 0.5e-5
        zmin in m of the linear sampling:   0.5e-5
        zmax in m of the linear sampling:   0.01
        Number of points per decade for the logarithmic sampling:   100
        exponent (10^) of zmin (in m) of the logarithmic sampling:  -2
        exponent (10^) of zmax (in m) of the logarithmic sampling:  6
        added z [m]:    
        Yokoya factors long, xdip, ydip, xquad, yquad:  1 0.651 0.898 -0.238 0.241
        factor weighting the longitudinal impedance error:  1.
        tolerance (in wake units) to achieve:   1.e11
        frequency above which the mesh bisecting is linear [Hz]:    1.e11
        Comments for the output files names:


- Description of each parameter:

    - The first parameters, up to the last layer thickness, are the same as for the round chamber
    impedance (see 4.1),
    - linear (1) or logarithmic (0) or both (2) scan in z for the wake: same
    kind of flag as above: we can perform either a purely linear scan (1), a
    purely logarithmic scan (0), or both (2). Note that the scan is over z, the
    distance (in m) BEHIND the travelling point-charge. Negative values can
    be entered in principle (meaning that we look for the wake in front of the beam).
    - sampling distance in m for the linear sampling:     distance between two
    points for the linear sampling.
    - zmin in m of the linear sampling: first point of the linear sampling.
    - zmax in m of the linear sampling: last point of the linear sampling.
    - Number of points per decade for the logarithmic sampling: self
    explanatory.
    - exponent (10^) of zmin (in m) of the logarithmic sampling: the log
    scan begins at 10^"this value".
    - exponent (10^) of zmax (in m) of the logarithmic sampling: the log
    scan ends at 10^"this value".
    - added z [m]: as above for frequencies: individual z values to be
    entered in the scan.
    As for the frequency scan, the final scan is sorted and duplicates are
    erased.
    - Yokoya factors long, xdip, ydip, xquad, yquad: see above in 4.1
    - factor weighting the longitudinal impedance error: to reach convergence
    on the wake functions, an error between the analytical impedances and its
    interpolation is integrated over the frequency range. This parameter is a factor
    weighting the error on the longitudinal impedance in the total error computation. It's useful to put a high value
    (100 or even 1000) if one wants to get an accurate longitudinal wake, because otherwise,
    since the transverse impedances are numerically usually much higher than the longitudinal one,
    the error on the longitudinal impedance is rather negligible in the integration, and therefore
    in the convergence process.
    - tolerance (in wake units) to achieve: this is the absolute error one wants to reach
    for the wake calculation. The code refines the frequency mesh used to compute the impedance
    until the error is less than this. Note that this is approximate only, the final wake functions
    error may still be above this parameter. Nevertheless, the lower it is, the more accurate will be the wake,
    and the longer it will take the code to reach convergence. A rule of thumb is to choose
    a few times below the minimum value the wake is likely to reach on the chosen z mesh. 
    - frequency above which the mesh bisecting is linear [Hz]: above this frequency, the refinement
    procedure bisects linearly the frequency intervals, whereas below it bisects them
    logarithmically. Changing this parameter does not matter too much. For ceramic structures, or more
    generally structures having lots of resonant peaks, one can try to put it lower (e.g. 1e9) because
    the resonant peaks are usually with equal frequencies between them.
    - Comments for the output files names: see above in 4.1.


###   5.3 Input file for the flat chamber impedance computation
       
- Input file example:

        Machine:    LHC
        Relativistic Gamma: 3730.26
        Impedance Length in m:  1.
        Maximum order of non-linear terms:  0
        Number of upper layers in the chamber wall: 2
        Layer 1 inner half gap in mm:   1.5
        Layer 1 DC resistivity (Ohm.m): 5.e-6
        Layer 1 relaxation time for resistivity (ps):   4.2
        Layer 1 real part of dielectric constant:   1
        Layer 1 magnetic susceptibility:    0
        Layer 1 relaxation frequency of permeability (MHz): Infinity
        Layer 1 thickness in mm:    25
        Layer 2 DC resistivity (Ohm.m): 7.2e-7
        Layer 2 relaxation time for resistivity (ps):   0
        Layer 2 real part of dielectric constant:   1
        Layer 2 magnetic susceptibility:    0
        Layer 2 relaxation frequency of permeability (MHz): Infinity
        Layer 2 thickness in mm:    Infinity
        Top bottom symmetry (yes or no):    yes
        start frequency exponent (10^) in Hz:   2
        stop frequency exponent (10^) in Hz:    9
        linear (1) or logarithmic (0) or both (2) frequency scan:   0
        sampling frequency exponent (10^) in Hz (for linear):   8
        Number of points per decade (for log):  20
        when both, fmin of the refinement (in THz): 5.5
        when both, fmax of the refinement (in THz): 7
        when both, number of points in the refinement:  100
        added frequencies (Hz):
        Comments for the output files names:


- Description of each parameter:

    See the round case above (4.1), except for:
    - Number of upper layers in the chamber wall: number of layers in the upper part of the chamber.
    - Layer 1 inner half gap in mm: position on the vertical axis of the chamber first upper layer, 
    with respect to the beam orbit.
    - Top bottom symmetry (yes or no): if "yes", the lower part of the chamber is exactly the same as
    the upper part, defined as in the round case (see 4.1) with the lines "Layer n ...". If "no", the
    chamber is asymetrical and one needs to define the lower part of the chamber, with the following
    lines:
    - Number of lower layers in the chamber wall: number of layers in the lower part of the chamber.
    - Layer -1 inner half gap in mm: position (in absolute value, i.e. >=0) on the vertical axis of the 
    chamber first lower layer, with respect to the beam orbit.
    - Layer -n DC resistivity (Ohm.m): DC resistivity rho of the n^th layer in the lower part 
    of the chamber (can be Infinity. e.g. for vacuum).
    - Layer -n relaxation time for resistivity (ps): relaxation time tau in
    picoseconds (BEWARE OF THE UNITS), for the AC conductivity (or
    resistivity). 0 for vacuum.
    - Layer -n real part of dielectric constant: eps_b (1 for vacuum or metals).
    - Layer -n magnetic susceptibility: it is chi_m=mu_r-1 if the permeability is
    real. 0 for vacuum.
    - Layer -n relaxation frequency of permeability (MHz): f_mu, cut-off frequency for the
    frequency-dependent complex permeability. Usually Inifinity (vacuum,
    metals, and others).
    - Layer -n thickness in mm: thickness (positive) of the n^th layer in the lower part 
    of the chamber. Can be Infinity (actually it has to be Infinity for the last layer).

Layer -1 to -n are working the same way as layers 1 to n (note that layers 1 and -1 are the closest
to the beam, and layers n and -n the farthest).
Note that the number of layers in the upper or lower part can be zero.

- There are also other options available for the layer propertes (flat impedance only):
For any layer:

        Layer 1 dielectric tan(delta):  0.1
        Layer 1 frequency dependent conductivity:       conductivity_vs_freq_hBN_example.dat
        Layer 1 frequency dependent relative complex permittivity:      relative_permittivity_vs_freq_hBN_example.dat
        Layer 1 frequency dependent relative complex permeability:      ferrite_8C11_muprime_musecond.dat

These are mutually exclusive: one cannot use a tan(delta) with a resistivity or relaxation time
(it will raise an error), nor a frequency dependent conductivity - resp. permittivity - 
with any other permittivity related quantity (resistivity, relaxation time, dielectric constant, 
tan(delta), freq. dependent permittivity - resp. freq. dependent conductivity).
The same goes for the permeability, but permeability/permittivity properties
are independent from each other so not mutually exclusive.

Frequency-dependent quantities are defined in files looking like (see also the example files
quoted in the lines above):

        Frequency [Hz]  eps_prime   eps_second
        0.      3.12385         -0.203868386315
        106.082 3.12385         -0.203868386315
        106.289 3.12356849182   -0.203850014556
        142.51  3.07431         -0.183599648141
        148.657 3.07207572348   -0.180577225441
        203.092 3.05229         -0.163971130987
        207.913 3.05058721854   -0.162512712659

(here for the relative complex permittivity). The first line is always
ignored (headers). For the conductivity, a complex sigma should be given,
but one can set to zero the imaginary part if needed (still, there should always be 3 columns).
Units are S.I (relatice complex permittivity/permeability are dimensionless,
conductivity is in S/m), and beware of the sign convension - eps'' and mu'' 
should be negative for positive freq. (see N. Mounet's PhD thesis and references therehin).

If "Maximum order of non-linear terms" is strictly higher than 0, it will also
give additional files of the same format, of name Z\<plane>\<a>\<b><c>\<d>\*.dat
with the general impedance term in front of x1^a y1^b x^c y^d , for all
orders such that a+b+c+d <= Max. order

NOTE: frequency-dependent quantities are linearly extrapolated ouside the frequency
range defined in the file, using the slope from the edges. So to avoid any
unphysical value, better extend the frequency range artificially with well
chosen values.



###   5.4 Input file for the flat chamber wake computation
       
- Input file example:

        Machine:    LHC
        Relativistic Gamma: 3730.26
        Impedance Length in m:  3.
        Number of upper layers in the chamber wall: 2
        Layer 1 inner half gap in mm:   6.5
        Layer 1 DC resistivity (Ohm.m): 1.5e-5
        Layer 1 relaxation time for resistivity (ps):   1.3
        Layer 1 real part of dielectric constant:   1
        Layer 1 magnetic susceptibility:    0
        Layer 1 relaxation frequency of permeability (MHz): Infinity
        Layer 1 thickness in mm:    25
        Layer 2 DC resistivity (Ohm.m): 7.2e-7
        Layer 2 relaxation time for resistivity (ps):   0
        Layer 2 real part of dielectric constant:   1
        Layer 2 magnetic susceptibility:    0
        Layer 2 relaxation frequency of permeability (MHz): Infinity
        Layer 2 thickness in mm:    Infinity
        Number of lower layers in the chamber wall: 0
        Top bottom symmetry (yes or no):    no
        linear (1) or logarithmic (0) or both (2) scan in z for the wake:   2
        sampling distance in m for the linear sampling: 0.5e-5
        zmin in m of the linear sampling:   0.5e-5
        zmax in m of the linear sampling:   0.01
        Number of points per decade for the logarithmic sampling:   100
        exponent (10^) of zmin (in m) of the logarithmic sampling:  -2
        exponent (10^) of zmax (in m) of the logarithmic sampling:  6
        added z [m]:    
        factor weighting the longitudinal impedance error:  1.
        tolerance (in wake units) to achieve:   1.e12
        frequency above which the mesh bisecting is linear [Hz]:    1.e10
        Comments for the output files names:    


- Description of each parameter:

    See the three sections above (in particular, section 4.3 for the parameters related to the 
    flat chamber properties, and section 4.2 for the z sampling parameters and the parameters related
    to the wake computation).


## 6. THE OUTPUT FILES

The output files (ascii files) are put in the directory from which the code is launched.
In each case, one of the output file (name begins with "InputData") is a repetition of the input file.
The first line of the impedance or wake output files contains columns headers with the description and unit
of each column. All units are SI units. The output files have all the .dat extension.
 
NOTE: the impedances (resp. wakes) are sorted with respect to ascending frequencies (resp. distances).
it is possible that some duplicate points remain in the output files (this can cause problems later when using
the files, for instance in an interpolation routine).

For the wake computation, the code prints on the standard output the iterations to get the prescribed accuracy on
the wake functions. This can be quite long.


###   6.1 Output files for the round chamber impedance computation
       
There are 5 output files, with three columns each: frequency, impedance real part,
impedance imaginary part.

Files are: longitudinal (begins with Zlong), x dipolar (Zxdip), y dipolar (Zydip), x quadrupolar
(Zxquad) and y quadrupolar (Zyquad) impedances.
       

###   6.2 Output files for the round chamber wake computation
       
There are 10 impedances output files, with three columns each: frequency, impedance real part,
impedance imaginary part. 5 files are with the (supposedly) converged frequency mesh,
and 5 files (with "_precise.dat" in the end) are with a finer mesh (twice more frequencies).

There are also 10 wake functions output files, with two columns each: distance behind
the source, and wake function. 5 files are computed with the (supposedly) converged frequency mesh,
and 5 files (with "_precise.dat" in the end) with a finer mesh (twice more frequencies).

For the impedances, the files are: longitudinal (Zlong), x dipolar (Zxdip), y dipolar (Zydip), 
x quadrupolar (Zxquad) and y quadrupolar (Zyquad) impedances.

For the wakes, the files are: longitudinal (Wlong), x dipolar (Wxdip), y dipolar (Wydip), 
x quadrupolar (Wxquad) and y quadrupolar (Wyquad) wakes.
       
       
###   6.3 Output files for the flat chamber impedance computation
       
There are 6 output files, with three columns each: frequency, impedance real part,
impedance imaginary part.

Files are: longitudinal (Zlong), x dipolar (Zxdip), y dipolar (Zydip), x quadrupolar
(Zxquad), y quadrupolar (Zyquad), and y contant (Zycst) impedances (the last one is non zero for 
an asymmetrical structure).
       

###   6.4 Output files for the flat chamber wake computation

There are 12 impedances output files, with three columns each: frequency, impedance real part,
impedance imaginary part. 6 files are with the (supposedly) converged frequency mesh,
and 6 files (with "_precise.dat" in the end) are with a finer mesh (twice more frequencies).

There are also 12 wake functions output files, with two columns each: distance behind
the source, and wake function. 6 files are computed with the (supposedly) converged frequency mesh,
and 6 files (with "_precise.dat" in the end) with a finer mesh (twice more frequencies).

For the impedances, the files are: longitudinal (Zlong), x dipolar (Zxdip), y dipolar (Zydip), 
x quadrupolar (Zxquad), y quadrupolar (Zyquad) impedances, and y contant impedances (Zycst).

For the wakes, the files are: longitudinal (Wlong), x dipolar (Wxdip), y dipolar (Wydip), 
x quadrupolar (Wxquad) and y quadrupolar (Wyquad), and y contant (Wycst) wakes.
       
       
###   6.5 How to plot the output files

You can use the python routines provided to plot the impedances and wakes (round or flat).
python2.6 with numpy and pylab need to be installed. The routines can be launched in the
command line.

Examples:

   ./plot_Imp_flat.py -f ZlongWLHC_2layersup_0layersdown6.50mm_precise.dat -g "first case" 
   -f ZlongWLHC_2layers18.38mm_precise.dat -g "second case" -l
       
This will plot in log-log scale all the impedances ending with "WLHC_2layersup_0layersdown6.50mm_precise.dat"
(long., x dip, y dip, x quad, y quad) , and all
ending with "WLHC_2layers18.38mm_precise.dat", putting on the same plot the longitudinal
impedances together, the x dipolar impedances together, etc., with the legends "first case" and 
"second case" on the graphs.
       
   ./plot_Wake_flat.py -f WlongWLHC_2layersup_0layersdown6.50mm_precise.dat -g "first case" 
   -f WlongWLHC_2layers18.38mm_precise.dat -g "second case" -l
       
This will plot in log-log scale all the wakes ending with "WLHC_2layersup_0layersdown6.50mm_precise.dat"
(long., x dip, y dip, x quad, y quad, and y cst if -c option is present) , and all
ending with "WLHC_2layers18.38mm_precise.dat", putting on the same plot the longitudinal
wakes together, the x dipolar wakes together, etc., with the legends "first case" and "second case" on
the graphs.
       

Short description of the options:

       -l (optional): plot in log-log scale (instead of semilogx for impedances and semilogy for wakes)
       -m (optional): units on the y axis are divided by meter (i.e. impedances & wakes are given per meter length)
       -c (optional): also plots the constant term (for flat asymmetrical chamber)
       -f: longitudinal file to be plotted (needs to begin by Zlong for plot_Imp_flat.py, or Wlong for
       plot_Wake_flat.py). There can be as many files as desired. The transverse impedances and wakes are loaded
       automatically (replacing e.g. Zlong by Zxdip, Zydip, etc.)
       -g (optional): legend to be associated with each file above (otherwise the name of the file is put 
       in the legends).
       -o (option): if used, all plots are directly put into files named e.g. Zlong[this_option].png, 
       Zlong[this_option].eps, Zxdip[this_option].png, etc., instead of printing them on the screen.
       This option is particularly useful if you use a local matplotlib without GUI (so without interactive plots)
       because then the default behaviour (printing on the screen) will fail.
       
   When computing the wakes, it is usually a good idea to compare the converged wake functions with the
   "precise" ones (ending with "_precise.dat"), to check the accuracy of the convergence: there should not be
   too much difference between the two wakes. Typically one would need to launch:
       
   ./plot_Wake_flat.py -f WlongWLHC_2layersup_0layersdown6.50mm_precise.dat -g "more precise wake" 
   -f WlongWLHC_2layersup_0layersdown6.50mm.dat -g "converged wake" -l
       
Note: when curves are superimposed so undistinguishable, it is easy to change one of the plot pattern
to take out the solid line and put instead for instance some crosses ('x') or circles ('o').
To do so, go in the code of the function plot_Imp_flat.py or plot_Wake_flat.py and change the line

color=['b','r','g','m','c','y','k'];

into (for example)

color=['b','xr','g','m','c','y','k'];


There is an example of running these routines in
command_plots_impedance.txt
       

## 7. TROUBLESHOOTING

### Installation
If there seems to be a problem loading GSL or MPFR (e.g. `cppyy.load_library(<package>)` raising an error), check for the existence of "`libgmp`", "`libmpfr`", "`libgsl`" and "`libgslcblas`" (with extensions `.a`, `.la` and/or `.so`) in `/usr/lib/` , `/usr/lib64/`, `/usr/local/lib/` or `<your_conda_installation_path>/envs/<your_env>/lib`
as well as the existence of "`gmp.h`", "`mpfr.h`" and the "`gsl`" directory in `/usr/include/`, `/usr/local/include/` or `<your_conda_path>/envs/<your_env>/include`.

---

The Arb library has a different name when downloaded with a Debian / Ubuntu package manager (`flint-arb`). Some cases where this causes problems may not have been handled correctly. If there seems to be a problem related to loading the Arb library, it is therefore likely caused by trying to load `arb` instead of `flint-arb` in `IW2D/loading.py` or `IW2D/cpp/Makefile`.

It may also happen that Arb is unable to find its dependency MPFR or GMP, in particular if Arb is installed using e.g. apt while MPFR is installed using conda. It is therefore advised to install all prerequisites using the same installer.

More troubleshooting advice for Arb is available at https://arblib.org/setup.html.

---

Both the python package and the Makefiles adds both your conda environment's and your system's library and header paths to the search path before compilation. If you have unpacked the local installation of GMP and MPFR, these are also automatically added to the library and include paths of the python package, whereas they are only added when compiling C++ executables when using the separate Makefile `Makefile_local_GMP_MPFR`.

This could potentially cause problems if there is a mismatch between two of these intallations. When debugging issues with these packages, it is therefore recommended to make sure you only have one copy of the libraries installed between your system, your conda environment, and locally in the repository.

---

If `cppyy` fails to install when running `pip install`, try installing it separately first. Installing with `conda` from `conda-forge` has been found to work. Do this by activating the conda environment you want to run `IW2D` from, and running:

    conda install -c conda-forge cppyy

Other installation procedures are available at: https://cppyy.readthedocs.io/en/latest/installation.html

### Using the code

The codes use floating point numbers with higher precision than doubles or even quadruple precision numbers.
This is required by the field matching problem resolution where very small and very high numbers appear.
This is also why we had to use the ALGLIB library (linear algebra with high precision floating point numbers).

In the round chamber computations, problems can arise due to the lack of accuracy when computing the 4*4
matrices, in particular due to a lack of precision of the modified Bessel functions.
In the flat case, some problems might arise when performing the numerical integration (with GSL).

In both cases, this is most often solved by increasing the precision of all the calculations. To do so, in
BOTH  ImpedanceWake2D/IW2D.h and ImpedanceWake2D/globals.h
   
   #define  Precision
   
by increasing the number after "Precision" (this is the number of bits for the mantissa of all numbers).

(Typically, 160 is fine. Putting more or equal to 200 for round, makes it
very very slow - it literally takes forever).
  
Note that the longitudinal wake function is often quite noisy (because it is quite small), as well as the 
quadrupolar wake functions for round chambers (when putting the Yokoya factors to "1 1 1 0 0"). Try to decrease
the tolerance on the wake convergence, in case this is an issue (see input files description).
 
  
## 8. LATEST MODIFICATIONS

- Restructured repository into python package
- Made python scripts to mimic behaviour of C++ executables
- Updated installation procedure
- New python interface
   
## 9. FUTURE DEVELOPMENTS (TO-DO LIST)
   
- Round chamber impedances and wakes: better implementation of the modified Bessel functions computation: now we use simply the
Taylor series for small arguments, and the asymptotic expansion for large arguments. A Chebyshev approximation
would be better in the intermediate range.
 
- Include the possibility to have the "perfect conductor" and "perfect magnet" kinds of boundary conditions
for the last layer.


## 10. REFERENCES

[1] B. Zotter. Longitudinal instabilities of charged particle beams inside cylindrical walls of finite thickness.
Part. Accel., 1:311-326, 1970.

[2] E. Keil and B. Zotter. The impedance of layered vacuum chambers. In 2nd European Particle Accelerator
Conference, Stockholm, Sweden, pages 963-965, 1998.

[3] B. Zotter. New results on the impedance of resistive metal walls of finite thickness, 2005. Report No.
CERN-AB-2005-043.

[4] E. Métral, B. Zotter, and B. Salvant. Resistive-wall impedance of an infinitely long multi-layer cylindrical
beam pipe. In 22nd Particle Accelerator Conference, Albuquerque, 2007.

[5] B. Salvant. Impedance model of the CERN SPS and aspects of LHC single-bunch stability, 
EPFL PhD Thesis 4585 (2010).

[6] Mathematica®, 2010. ©Wolfram Research, Inc., Champaign, Illinois , USA.

[7] N. Mounet and E. Métral. Electromagnetic Fields Created by a Beam in an
Infinitely Long and Axisymmetric Multilayer Resistive Pipe, Report No. CERN-BE-2009-039 (2009).

[8] N. Mounet. The LHC Transverse Coupled-Bunch Instability, PhD thesis 5305 (EPFL, 2012),
http://infoscience.epfl.ch/record/174672/files/EPFL_TH5305.pdf

[9] N. Mounet and E. Métral. Electromagnetic fields and beam coupling impedances in a multilayer flat chamber,
CERN-ATS-Note-2010-056 TECH. (2010).

[10] N.Mounet and E.Métral. Impedances of two dimensional multilayer cylindrical and flat chambers
in the non-ultrarelativistic case. In HB2010, 46th ICFA Advanced Beam Dynamics Workshop on
High-Intensity and High-Brightness Hadron Beams, Morschach, Switzerland, pages 353-357, 2010.

[11] F. Roncarolo, F. Caspers, T. Kroyer, E. Métral, N. Mounet, B. Salvant, and B. Zotter. Comparison between
laboratory measurements, simulations and analytical predictions of the transverse wall impedance at low
frequencies. Phys. Rev. ST Accel. Beams, 12(084401), 2009.
