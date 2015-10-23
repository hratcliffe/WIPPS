//
//  support.h
//  
//
//  Created by Heather Ratcliffe on 24/09/2015.
//
//

#ifndef _support_h
#define _support_h


#include <fstream>
#include <iostream>
#include <stdio.h>

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

#define calc_type double


const int MAX_SIZE = 10000;
const my_type io_verify = 3.0/32.0;
//An exactly binary representable my_type to verify we're reading what we're writing...

const calc_type pi = 3.14159;
const calc_type v0 = 2.997924e8; //m/s^2
const calc_type me = 9.10938291e-31; //kg
const calc_type mp = me*1836.15267; //kg
const calc_type q0 = 1.602176565e-19; // C
const calc_type eps0 =8.85418782e-12; //F/m

//kb = 1.3806488d-23   ; J/K
//mu0 = 4.0d-7 * !dpi  ; N/A^2
//epsilon0 = 8.8541878176203899d-12 ; F/m
//h_planck = 6.62606957d-34 ; J s

struct deck_constants{
  float v_t;
  float omega_pe;
  float omega_ce;
  float omega_ci;
  int ppc;
  

};

struct mpi_info_struc{

  int rank;
  int n_procs;
};

struct mu{
  calc_type mu;
  calc_type mug;
  calc_type dmudr;
  calc_type dmudlat;
  calc_type dmudtheta;
  calc_type dmudom;
  calc_type alpha;
  int err;

};

struct mu_dmudom{
  calc_type mu;
  calc_type dmudom;
  calc_type phi;
  int err;

};


const int DEFAULT_N_ANG = 100;

void my_print(std::string text, int rank, int rank_to_write=0);

std::string mk_str(int i);/**<Converts int to string*/
std::string mk_str(bool b);/**<Converts bool to string*/
//std::string mk_str(size_t i){ return mk_str((int) i);} /**<Converts size_t to string*/

calc_type integrator(calc_type * start, int len, calc_type * increment);


const int WAVE_WHISTLER = 1;
const int WAVE_PLASMA = 2;
const int WAVE_O = 3;

const int FUNCTION_DELTA = 1;
const int FUNCTION_GAUSS = 2;
const int FUNCTION_ISO = 3;


#endif
