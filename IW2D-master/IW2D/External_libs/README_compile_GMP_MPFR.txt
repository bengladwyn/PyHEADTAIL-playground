# To compile GMP:
# copy somewhere (i.e. NOT in this repository) the file gmp-6.0.0a.tar.gz
# then
tar -xzvf gmp-6.0.0a.tar.gz
cd gmp-6.0.0/
./configure
make


# To compile MPFR:
# copy somewhere (i.e. NOT in this repository) the file mpfr-3.1.2.tar.gz
# then
tar -xzvf mpfr-3.1.2.tar.gz
cd mpfr-3.1.2/

# configure using the local GMP (put the absolute path to the GMP directory you just installed, after "--with-gmp")
./configure --with-gmp-lib=/home/inewton/GMP/gmp-6.0.0 --with-gmp-include=/home/inewton/GMP/gmp-6.0.0 --disable-thread-safe

# then change all options -O3 into -O2 in files "libtool", "Makefile" and "config.status"
# and finally do
make


Note: you might need to install something called 'm4' during the process. Under Ubuntu you can install it via:
sudo apt-get install m4
