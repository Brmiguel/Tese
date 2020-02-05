# Installations 

## GFortran

- $sudo apt-get install gfortran

## Library fftw3
-Go to http://www.fftw.org/download.html, download file fftw-3.3.8.tar.gz .
-Extract the folder and open it in terminal.
- $ export CC="gcc -m64" 
- $ export MPICC="mpicc -m64"
- $ export F77="gfortran -m64"
- $ ./configure --enable-mpi --enable-threads --enable-openmp --enable-avx 
- make 
- make install

## Library cfitsio

-Go to https://heasarc.gsfc.nasa.gov/fitsio/ , download file cfitsio-3.47.tar.gz .
-Extract the folder and open it in terminal.
- $ ./configure --prefix=/usr/local/ --enable-sse2 --enable-ssse3 
- $ make 
- $ install

## Install g++ Compiler

sudo apt install g++

## Install PGI Compiler

-Go to https://www.pgroup.com/support/download_community.php?file=pgi-community-linux-x64 and download PGI Community Edition Version 

- $ sudo .\install 
- $ export PGI=/opt/pgi;
- $ export PATH=/opt/pgi/linux86-64/19.4/bin:$PATH;
- $ export MANPATH=$MANPATH:/opt/pgi/linux86-64/19.4/man;


## Compiling code without parallelization

pgfortran diffraction_nompi.f90 -o diffraction_nompi -I/usr/local/include -lfftw3_threads -lfftw3 -lm  -lcfitsio -lnsl -L/lib/x86_64-linux-gnu


## Compiling code without parallelization

pgfortran diffraction_nompi.f90 -o diffraction_nompi -I/usr/local/include -lfftw3_threads -lfftw3 -lm  -lcfitsio -lnsl -L/lib/x86_64-linux-gnu -fast -ta=multicore
