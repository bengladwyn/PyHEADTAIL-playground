from setuptools import setup, find_packages
from pathlib import Path

# Read README.md and use as long description in package metadata
root_dir = Path(__file__).parent.absolute()
with (root_dir / "README.md").open("rt") as f:
    long_description = f.read().strip()


# TODO: Add version requirements?
requirements = {
    "core": [
        "numpy",
        "cppyy<3.0.0",
        "pandas",
        "scipy"
    ],
    "dev": ["pre-commit"],
    "test":["pytest"]
}

setup(
    name='IW2D',
    version='0.1.0',
    description='Codes for calculating impedances and wake functions in 2D for round and flat chambers',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://gitlab.cern.ch/IRIS/IW2D/',
    author='Nicolas Mounet, Nicolò Biancacci, David Amorim, Eskil Vik',
    packages=find_packages() # finds all the packages (i.e. all __init__.py files) in the repository
            + ["IW2D.cpp", "IW2D.External_libs"], # include packages that do not contain any python files
    python_requires=">=3.7", # >=3.7 might work, but only 3.8 has been used for development
    install_requires=requirements["core"],
    extras_require= requirements,
    include_package_data=True
)
