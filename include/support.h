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
#include <vector>

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

#define tiny_calc_type 1e-12


const int MAX_SIZE = 10000;/**< Maximum array size allowed (per processor if MPI in use) */
const int MAX_FILENAME_DIGITS = 7;/**< Maximum number of digits in filename dump number string*/
const my_type io_verify = 3.0/32.0;/**< An exactly binary representable my_type to verify we're reading what we're writing.*/
const calc_type pi = 3.14159265359;/**< Pi */
const calc_type v0 = 2.997924e8; /**< Speed of light in m/s^2 */
const calc_type me = 9.10938291e-31; /**< Electron mass in kg */
const calc_type mp = me*1836.15267; /**< Proton mass in kg */
const calc_type q0 = 1.602176565e-19; /**< Electron charge in C*/
const calc_type eps0 =8.85418782e-12; /**< Epsilon_0 permittivity of free space in F/m*/

//kb = 1.3806488d-23   ; J/K
//mu0 = 4.0d-7 * !dpi  ; N/A^2
//epsilon0 = 8.8541878176203899d-12 ; F/m
//h_planck = 6.62606957d-34 ; J s

/** \brief Constants read from deck
*
*Holds run parameters extracted from deck file such as temperature, reference frequencies etc
*/
struct deck_constants{
  float v_t; /**< Electron thermal velocity */
  float omega_pe; /**< Plasma frequency (reference) */
  float omega_ce; /**< Electron cyclotron frequency (reference)*/
  float omega_ci; /**< Ion cyclotron frequency (reference)*/
  int ppc;/**< Particles per cell used */
  
};

/** \brief MPI information
*
*Holds info on MPI: processor ranks etc
*/
struct mpi_info_struc{

  int rank;
  int n_procs;
};

/** \brief Full refractive index
*
*Contains refractive index mu and all derivatives, plus error flag
*/
struct mu{
  calc_type mu; /**< Refractive index */
  calc_type mug; /**< */
  calc_type dmudr; /**< d mu /dr (radial distance) */
  calc_type dmudlat; /**< d mu/ d latitude */
  calc_type dmudtheta; /**< d mu / d theta (wave normal angle) */
  calc_type dmudom; /**< d mu / d omega (wave frequency)*/
  calc_type alpha; /**< */
  int err; /**< 0 if mu found successfully, 1 else*/

};

/** \brief Reduced refractive index
*
*Contains refractive index mu, the two derivative needed for diffusion coefficient calculation, the phi function of e.g. Albert 2005 and error flag
*/
struct mu_dmudom{
  calc_type mu; /**< Refractive index */
  calc_type dmudom; /**< d mu / d omega (wave frequency)*/
  calc_type dmudtheta; /**< d mu / d theta (wave normal angle) */
  calc_type phi;/**< Phi from Albert 2005*/
  int err; /**< 0 if mu found successfully, 1 else*/

};

const int DEFAULT_N_ANG = 100;/**< Default number of wave normal angles to consider*/

void my_print(std::string text, int rank, int rank_to_write=0);
void my_print(std::fstream * handle, std::string text, int rank, int rank_to_write=0);

std::string mk_str(int i);/**<Converts int to string*/
std::string mk_str(bool b);/**<Converts bool to string*/
//std::string mk_str(size_t i){ return mk_str((int) i);} /**<Converts size_t to string*/
std::string mk_str(double i);
std::string mk_str(float i);
std::string mk_str(long double i);

//int whereb(my_type * ax_ptr, int len, my_type target, int &cut,int sign=1.0);
int where(my_type * ax_ptr, int len, my_type target);


template<typename T> T interpolate(T* axis, T* vals, T target, int pts);
template<typename T> T integrator(T * start, int len, T * increment);
template<typename T> void inplace_boxcar_smooth(T * start, int len, int width, bool periodic = 0);
//void inplace_boxcar_smooth(calc_type * start, int len, int width, bool periodic = 0);
calc_type square_integrator(calc_type * start, int len, calc_type * increment);

std::vector<calc_type> cubic_solve(calc_type a, calc_type b, calc_type c);

const int WAVE_WHISTLER = 1; /**< Code to id wave as whistler mode */
const int WAVE_PLASMA = 2; /**< Code to id wave as plasma/Langmuir mode */
const int WAVE_O = 3; /**< Code to id wave as ordinary EM mode */

const int FUNCTION_NULL = 0; /**< Code to id spectral angular distribution as absent */
const int FUNCTION_DELTA = 1; /**< Code to id spectral angular distribution as delta function (with integral 1) */
const int FUNCTION_GAUSS = 2; /**< Code to id spectral angular distribution as gaussian (with integral 1) */
const int FUNCTION_ISO = 3; /**< Code to id spectral angular distribution as isotropic (with integral 1) */


#endif
