#!/bin/bash

###
# Bash script for ImpedanceWake2D installation on Linux systems. It works on LXPLUS at CERN (CC7 and SLC6) or under Ubuntu (versions 10.04, 12.04, 16.04, 20.04)
# This has to be executed after cloning the repository.
# It installs GMP and MPFR and compiles the code.
###

set -e

#####################
# Installation set-up
#####################

# Folder in which the installation will be performed
SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IW2D_PATH="$(cd "$SCRIPT_PATH/../" && pwd)"


###########################
# Extract and compile GMP
###########################
cd $IW2D_PATH/External_libs
tar -zxvf gmp-6.0.0a.tar.gz
cd $IW2D_PATH/External_libs/gmp-6.0.0
./configure
make

############################
# Extract and compile MPFR
############################
cd $IW2D_PATH/External_libs
tar -zxvf mpfr-3.1.2.tar.gz
cd $IW2D_PATH/External_libs/mpfr-3.1.2
./configure --with-gmp-lib=$IW2D_PATH/External_libs/gmp-6.0.0 --with-gmp-include=$IW2D_PATH/External_libs/gmp-6.0.0 --disable-thread-safe
make
