WIPPS post-processes output from the EPOCH (University of Warwick) PIC code to explore Whistler-mode wave-particle interactions. It can create wave spectra, and calculate theory based growth rates and particle diffusion coefficients.

Requires:
* Make and a C++ 11 compatible compiler
* boost (special functions, filesystem etc)
* FFTW (Fourier transform functionality only)
* EPOCH's SDF file reader (included, SDF file inputs only)

This code includes a copy of the SDF file library from EPOCH, licensed under the 3-clause BSD license. The SDF development team do not endorse this code. 
A copy of valgrind.h (http://valgrind.org/) is included under BSD license for testing purposes. 
This code is under development: features and calculations may be incomplete, incorrect or limited in scope. 

