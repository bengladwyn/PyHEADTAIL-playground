import subprocess
from pathlib import Path
import cppyy
import os

def install_external_libraries():
    """Unpack and compile the GMP and MPFR libraries from the tarballs located in 'IW2D/External_libs/'.
    This function is equivalent to executing the compilation script in the same folder directly.
    Refer to the README.md for complete installation procedure."""
    package_path = Path(__file__).resolve().parent
    subprocess.run(package_path / "External_libs" / "compile_ext_libs.sh")

def load_libraries():
    """Loads the libraries gsl, gslcblas, m, gmp and mpfr and adds their headers to the include path.
    Adds alglib to the include path and includes ap.cpp and amp.cpp.
    Adds The repository itself to the include path."""

    package_path = Path(__file__).parent.resolve()

    # If GMP is installed locally, add to path
    if (package_path / "External_libs/gmp-6.0.0").exists():
        cppyy.add_library_path(str(package_path / "External_libs/gmp-6.0.0/.libs"))
        cppyy.add_include_path(str(package_path / "External_libs/gmp-6.0.0")) # Header location
    
    # If MPFR is installed locally, add to path
    if (package_path / "External_libs/mpfr-3.1.2/src").exists():
        cppyy.add_library_path(str(package_path / "External_libs/mpfr-3.1.2/src/.libs"))
        cppyy.add_include_path(str(package_path / "External_libs/mpfr-3.1.2/src")) # Header location


    #### Add include path for cpp scripts, including ALGLIB #####
    cppyy.add_include_path(str(package_path / "cpp"))

    # Load libraries
    cppyy.load_library("gslcblas")
    cppyy.load_library("gsl")
    # cppyy.load_library("stdc++")
    cppyy.load_library("gmp")
    cppyy.load_library("mpfr")
    # cppyy.load_library("flint")
    # cppyy.load_library("m")
    
    # The Arb library is called "flint-arb" when installed with a Debian/Ubuntu package manager.
    # As a simple solution to this, we force the user to set the environment variable IW2D_FLINT_ARB
    # if this is the case, otherwise the more standard "arb" will be loaded as a default.
    # More information about arb setup: https://arblib.org/setup.html
    if os.environ.get("IW2D_FLINT_ARB", None) is not None:
        cppyy.load_library("flint-arb")
    else:
        cppyy.load_library("arb")


    # Include definitions from ap.cpp and amp.cpp (which is not a shared library)
    cppyy.include("ap.cpp")
    cppyy.include("amp.cpp")

def load_iw2d():
    cppyy.include("multiprecision.cc")
    cppyy.include("interp.cc")
    cppyy.include("input.cc")
    cppyy.include("flatchamber.cc")
    cppyy.include("flat.cc")
    cppyy.include("round.cc")

# Helper function to avoid errors if files are attempted loaded multiple times.
# Removed since everything is now included only in __init__.py, and therefore only once
# def safe_include(*source_files: Union[Path, str]):
#     for source_file in source_files:
#         try:
#             cppyy.include(source_file)
#         except ImportError as err:
#             print(f"Include {source_file} failed")
