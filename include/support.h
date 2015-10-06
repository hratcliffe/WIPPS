//
//  support.h
//  
//
//  Created by Heather Ratcliffe on 24/09/2015.
//
//

#ifndef _support_h
#define _support_h

#define _USE_FLOAT
//temporary. This should go to makefile

#ifdef _USE_FLOAT

#define ADD_FFTW(x) fftwf_ ## x
#define cplx_type ADD_FFTW(complex)
#define my_type float
#define my_sdf_type SDF_DATATYPE_REAL4

#else

#define ADD_FFTW(x) fftw_ ## x
#define cplx_type ADD_FFTW(complex)
#define my_type double
#define my_sdf_type SDF_DATATYPE_REAL8

#endif
/** These set up our types so we can easily recompile to work with doubles or floats. First adds correct FFTW library prefix, adjust to float or normal. Next defines suitable complex type. Third s working data type, and the last is set to suitable SDF data type matching my_type. Lets be sane, and assume we want the f libraries to work with float data, and the double to work with doubles. So we don't have extraneous copying and false precision.
*/

const int MAX_SIZE = 10000;
const my_type io_verify = 3.0/32.0;
//An exactly binary representable my_type to verify we're reading what we're writing...

const float pi = 3.14159;
const float v0 = 2.997924e8; 

//q0 = 1.602176565d-19 ; C
//m0 = 9.10938291d-31  ; kg
//v0 = 2.99792458d8    ; m/s^2
//kb = 1.3806488d-23   ; J/K
//mu0 = 4.0d-7 * !dpi  ; N/A^2
//epsilon0 = 8.8541878176203899d-12 ; F/m
//h_planck = 6.62606957d-34 ; J s

struct deck_constants{
  float v_t;
  float omega_pe;
  float omega_ce;
  int ppc;
  

};

struct mpi_info_struc{

  int rank;
  int n_procs;
};

const int DEFAULT_N_ANG = 100;

void my_print(std::string text, int rank, int rank_to_write=0);

std::string mk_str(int i);/**<Converts int to string*/
std::string mk_str(bool b);/**<Converts bool to string*/
//std::string mk_str(size_t i){ return mk_str((int) i);} /**<Converts size_t to string*/


#endif
